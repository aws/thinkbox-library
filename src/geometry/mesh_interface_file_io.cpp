// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/mesh_interface_file_io.hpp>

#include <boost/bind.hpp>
#include <boost/cstdint.hpp>
#ifndef FRANTIC_DISABLE_THREADS
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#endif

#include <tbb/atomic.h>

#include <frantic/files/compression_stream.hpp>
#include <frantic/files/files.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/rply.hpp>
#include <frantic/geometry/xmesh_writer.hpp>
#include <frantic/locale/locale.hpp>
#include <frantic/threads/synchronizedqueue.hpp>

namespace {
/**
 * An interface for managing temp files.
 * If close is called then it will finalize the file.
 * If the file manager falls out of scope and close was not called
 *    then it will delete the temp file
 */
class temp_file_manager {
    frantic::tstring m_targetFile;
    frantic::tstring m_tempFile;
    bool m_isClosed;

  public:
    temp_file_manager( const frantic::tstring& targetFile, bool useTempFile = true )
        : m_targetFile( targetFile )
        , m_isClosed( false ) {
        if( useTempFile ) {
            boost::filesystem::path tempPath =
                boost::filesystem::unique_path( targetFile + _T(".%%%%-%%%%-%%%%-%%%%.tmp") );
            m_tempFile = frantic::files::to_tstring( tempPath );
        } else {
            m_tempFile = m_targetFile;
        }
    }
    const frantic::tstring& get_path() const { return m_tempFile; }
    void close() {
        boost::filesystem::path tempFile( m_tempFile ), targetFile( m_targetFile );
        if( tempFile != targetFile ) {
            boost::system::error_code errcode;

            // Delete the target file (if it exists) so that we can rename our new file to the final name it needs.
            boost::filesystem::remove( targetFile, errcode );
            if( errcode )
                throw std::runtime_error( "temp_file_manager failed to remove existing file " +
                                          frantic::strings::to_string( m_targetFile ) + " : " +
                                          frantic::strings::to_string( errcode.message() ) );

            // TODO check if tempFile exists
            boost::filesystem::rename( tempFile, targetFile, errcode );
            if( errcode )
                throw std::runtime_error( "temp_file_manager fail to rename " +
                                          frantic::strings::to_string( m_tempFile ) + " to " +
                                          frantic::strings::to_string( m_targetFile ) + " : " +
                                          frantic::strings::to_string( errcode.message() ) );
        }
        m_isClosed = true;
    }
    ~temp_file_manager() {
        if( !m_isClosed ) {
            boost::filesystem::path tempFile( m_tempFile );
            boost::system::error_code errcode;
            // try to delete the temp file (if it exists) because the file manager fell out of scope without being
            // closed
            boost::filesystem::remove( tempFile, errcode );
            if( errcode )
                FF_LOG( warning ) << "temp_file_manager failed to remove " << m_tempFile << " : "
                                  << frantic::files::to_tstring( errcode.message() );
        }
    }
};
} // anonymous namespace

namespace frantic {
namespace geometry {

namespace {
#ifndef FRANTIC_DISABLE_THREADS
// A simple thread pool for saving xmesh channels in parallel.
class thread_pool {
    boost::asio::io_service m_service;
    boost::shared_ptr<boost::asio::io_service::work> m_work;
    boost::thread_group m_threads;
    boost::thread* m_thread;

  public:
    thread_pool( std::size_t threadCount )
        : m_service( static_cast<int>( std::max<std::size_t>( threadCount, 1 ) ) )
        , m_work( new boost::asio::io_service::work( m_service ) ) {
        threadCount = std::max<std::size_t>( threadCount, 1 );
        for( std::size_t i = 0; i < threadCount; ++i ) {
            m_thread = m_threads.create_thread( boost::bind( &boost::asio::io_service::run, &m_service ) );
        }
    }

    ~thread_pool() {
        m_work.reset();
        m_threads.join_all();
    }

    template <typename F>
    void schedule( F task ) {
        m_service.post( task );
    }
};
#else
class thread_pool {
  public:
    thread_pool( std::size_t /*threadCount*/ ) {}

    template <typename F>
    void schedule( F task ) {
        task();
    }
};
#endif

typedef frantic::threads::SynchronizedQueue<std::string> error_queue_t;

// thread safe way to update progress info
class progress_info_wrapper {
    frantic::logging::progress_logger& m_progressLogger;

    progress_info_wrapper( const progress_info_wrapper& );            // not implemented
    progress_info_wrapper& operator=( const progress_info_wrapper& ); // not implemented

  public:
    tbb::atomic<bool> abortRequested;
    tbb::atomic<std::size_t> progressElementsPassed;
    tbb::atomic<std::size_t> progressTotalCount;

    progress_info_wrapper( frantic::logging::progress_logger& progress )
        : m_progressLogger( progress ) {
        abortRequested = false;
        progressElementsPassed = 0;
        progressTotalCount = 0;
    }

    // increases the number of elements passed and checks if an abort has been requested
    void increment_progress( std::size_t amountPassed ) {
        progressElementsPassed += amountPassed;
#ifdef FRANTIC_DISABLE_THREADS
        // if threads are disabled it is safe to update the progress directly
        m_progressLogger.update_progress( progressElementsPassed, progressTotalCount );
#endif
        if( abortRequested )
            throw frantic::logging::progress_cancel_exception( "Operation cancelled in another thread" );
    }
};

// returns a sum of the number of verts, elements and faces in a mesh
std::size_t get_progress_total_count( const mesh_interface_ptr& mesh ) {
    const mesh_interface::mesh_channel_map& vertexChannels = mesh->get_vertex_channels();
    const mesh_interface::mesh_channel_map& faceChannels = mesh->get_face_channels();

    std::size_t progressTotalCount = 0;

    progressTotalCount += mesh->get_num_verts();
    progressTotalCount += mesh->get_num_faces();

    mesh_interface::mesh_channel_map::const_iterator it = vertexChannels.begin();
    mesh_interface::mesh_channel_map::const_iterator itEnd = vertexChannels.end();
    for( ; it != itEnd; ++it ) {
        progressTotalCount += it->second->get_num_elements();
        progressTotalCount += it->second->get_num_faces();
    }
    it = faceChannels.begin();
    itEnd = faceChannels.end();
    for( ; it != itEnd; ++it ) {
        progressTotalCount += it->second->get_num_elements();
    }
    return progressTotalCount;
}

void get_channel_element_data( const mesh_channel* channel, std::vector<char>& outData ) {
    outData.resize( channel->get_element_size() * channel->get_num_elements() );
    for( std::size_t i = 0; i < channel->get_num_elements(); ++i ) {
        channel->get_value( i, &outData[i * channel->get_element_size()] );
    }
}
bool has_custom_faces( const mesh_channel* channel, const mesh_interface_ptr& mesh ) {
    if( mesh->get_num_faces() != channel->get_num_faces() )
        return true;
    for( std::size_t i = 0; i < channel->get_num_faces(); ++i ) {
        if( mesh->get_num_face_verts( i ) != channel->get_num_face_verts( i ) )
            return true;
        for( std::size_t k = 0; k < channel->get_num_face_verts( i ); ++k ) {
            if( mesh->get_face_vert_index( i, k ) != channel->get_fv_index( i, k ) )
                return true;
        }
    }
    return false;
}

void process_geom_verts( xmesh_writer& fileWriter, const mesh_interface_ptr& mesh, progress_info_wrapper& progress,
                         error_queue_t& errorQueue ) {
    bool done = false;
    std::string errMsg;
    try {
        // check for cancellation before starting any work
        progress.increment_progress( 0 );
        // write geom verts
        std::vector<float> dataBuffer( 3 * mesh->get_num_verts() );
        for( std::size_t i = 0; i < mesh->get_num_verts(); ++i ) {
            mesh->get_vert( i, *reinterpret_cast<float( * )[3]>( &dataBuffer[3 * i] ) );
        }
        progress.increment_progress( 0 ); // check for cancellation
        fileWriter.write_vertex_channel( _T("verts"),
                                         dataBuffer.size() > 0 ? reinterpret_cast<char( * )>( &dataBuffer[0] ) : 0,
                                         frantic::channels::data_type_float32, 3, mesh->get_num_verts() );
        progress.increment_progress( mesh->get_num_verts() );

        done = true;
    } catch( std::exception& e ) {
        errMsg = e.what();
    } catch( ... ) {
        errMsg = "process_geom_verts_and_faces Error: An unknown exception occurred.";
    }
    if( !done ) {
        errorQueue.enter( errMsg );
    }
}
void process_geom_faces( xmesh_writer& fileWriter, const mesh_interface_ptr& mesh, progress_info_wrapper& progress,
                         error_queue_t& errorQueue ) {
    bool done = false;
    std::string errMsg;
    try {
        // check for cancellation before starting any work
        progress.increment_progress( 0 );
        // write geom faces
        if( mesh->get_num_faces() > 0 ) {
            std::size_t faceBufferSize = 0;
            for( std::size_t i = 0; i < mesh->get_num_faces(); ++i ) {
                faceBufferSize += mesh->get_num_face_verts( i );
            }
            std::vector<int> faceBuffer( faceBufferSize );

            std::size_t pos = 0;
            for( std::size_t currFace = 0; currFace < mesh->get_num_faces(); ++currFace ) {
                for( std::size_t currFaceVert = 0; currFaceVert < mesh->get_num_face_verts( currFace );
                     ++currFaceVert ) {
                    faceBuffer[pos] = (int)mesh->get_face_vert_index( currFace, currFaceVert );
                    ++pos;
                }
                // mark the last vertex in each face by making it negitive ( and -1 to handle the 0 case )
                faceBuffer[pos - 1] = -faceBuffer[pos - 1] - 1;
            }
            progress.increment_progress( 0 ); // check for cancellation
            fileWriter.write_vertex_channel_faces( _T("verts"), faceBuffer.size() > 0 ? &faceBuffer[0] : 0,
                                                   faceBuffer.size() );

            progress.increment_progress( mesh->get_num_faces() );
        }
        done = true;
    } catch( std::exception& e ) {
        errMsg = e.what();
    } catch( ... ) {
        errMsg = "process_geom_verts_and_faces Error: An unknown exception occurred.";
    }
    if( !done ) {
        errorQueue.enter( errMsg );
    }
}
void process_vert_channel( xmesh_writer& fileWriter, const mesh_channel* channel, progress_info_wrapper& progress,
                           error_queue_t& errorQueue ) {
    bool done = false;
    std::string errMsg;
    try {
        // check for cancellation before starting any work
        progress.increment_progress( 0 );
        std::vector<char> dataBuffer;
        get_channel_element_data( channel, dataBuffer );
        progress.increment_progress( 0 ); // check for cancellation
        fileWriter.write_vertex_channel( channel->get_name(), dataBuffer.size() > 0 ? &dataBuffer[0] : 0,
                                         channel->get_data_type(), channel->get_data_arity(),
                                         channel->get_num_elements() );

        progress.increment_progress( channel->get_num_elements() );
        done = true;
    } catch( std::exception& e ) {
        errMsg = e.what();
    } catch( ... ) {
        errMsg = "process_vert_channel Error: An unknown exception occurred.";
    }
    if( !done ) {
        errorQueue.enter( errMsg );
    }
}
void process_vert_channel_faces( xmesh_writer& fileWriter, const mesh_channel* channel, const mesh_interface_ptr& mesh,
                                 progress_info_wrapper& progress, error_queue_t& errorQueue ) {
    bool done = false;
    std::string errMsg;
    try {
        // check for cancellation before starting any work
        progress.increment_progress( 0 );
        if( has_custom_faces( channel, mesh ) ) {
            // There are custom faces. We need to translate them into the appropriate encoding.
            std::size_t faceBufferSize = 0;
            for( std::size_t i = 0; i < channel->get_num_faces(); ++i ) {
                faceBufferSize += channel->get_num_face_verts( i );
            }
            std::vector<int> faceBuffer( faceBufferSize );

            std::size_t pos = 0;
            for( std::size_t currFace = 0; currFace < channel->get_num_faces(); ++currFace ) {
                for( std::size_t currFaceVert = 0; currFaceVert < channel->get_num_face_verts( currFace );
                     ++currFaceVert ) {
                    faceBuffer[pos] = (int)channel->get_fv_index( currFace, currFaceVert );
                    ++pos;
                }
                // mark the last vertex in each face by making it negitive ( and -1 to handle the 0 case )
                faceBuffer[pos - 1] = -faceBuffer[pos - 1] - 1;
            }
            progress.increment_progress( 0 ); // check for cancellation
            fileWriter.write_vertex_channel_faces( channel->get_name(), faceBuffer.size() > 0 ? &faceBuffer[0] : 0,
                                                   faceBuffer.size() );
        }
        progress.increment_progress( channel->get_num_faces() );

        done = true;
    } catch( std::exception& e ) {
        errMsg = e.what();
    } catch( ... ) {
        errMsg = "process_geom_verts_and_faces Error: An unknown exception occurred.";
    }
    if( !done ) {
        errorQueue.enter( errMsg );
    }
}
void process_face_channel( xmesh_writer& fileWriter, const mesh_channel* channel, progress_info_wrapper& progress,
                           error_queue_t& errorQueue ) {
    bool done = false;
    std::string errMsg;
    try {
        // check for cancellation before starting any work
        progress.increment_progress( 0 );
        std::vector<char> dataBuffer;
        get_channel_element_data( channel, dataBuffer );
        progress.increment_progress( 0 ); // check for cancellation
        fileWriter.write_face_channel( channel->get_name(), dataBuffer.size() > 0 ? &dataBuffer[0] : 0,
                                       channel->get_data_type(), channel->get_data_arity(),
                                       channel->get_num_elements() );

        progress.increment_progress( channel->get_num_elements() );
        done = true;
    } catch( std::exception& e ) {
        errMsg = e.what();
    } catch( ... ) {
        errMsg = "process_face_channel Error: An unknown exception occurred.";
    }
    if( !done ) {
        errorQueue.enter( errMsg );
    }
}

void write_xmesh_mesh_file_worker( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                                   xmesh_metadata metadata, progress_info_wrapper& progress, std::size_t threadCount,
                                   error_queue_t& errorQueue ) {
    bool done = false;
    std::string errMsg;
    try {
        metadata.set_boundbox( compute_boundbox( mesh ) );

        xmesh_writer fileWriter( path );
        fileWriter.set_metadata( metadata );

        { // scope for threadPool
            thread_pool threadPool( std::max<std::size_t>( 1, threadCount ) );

            const mesh_interface::mesh_channel_map& vertexChannels = mesh->get_vertex_channels();
            const mesh_interface::mesh_channel_map& faceChannels = mesh->get_face_channels();

            threadPool.schedule( boost::bind( process_geom_verts, boost::ref( fileWriter ), mesh,
                                              boost::ref( progress ), boost::ref( errorQueue ) ) );
            threadPool.schedule( boost::bind( process_geom_faces, boost::ref( fileWriter ), mesh,
                                              boost::ref( progress ), boost::ref( errorQueue ) ) );

            // write vert channels
            mesh_interface::mesh_channel_map::const_iterator it = vertexChannels.begin();
            mesh_interface::mesh_channel_map::const_iterator itEnd = vertexChannels.end();
            for( ; it != itEnd; ++it ) {
                threadPool.schedule( boost::bind( process_vert_channel, boost::ref( fileWriter ), it->second,
                                                  boost::ref( progress ), boost::ref( errorQueue ) ) );
                threadPool.schedule( boost::bind( process_vert_channel_faces, boost::ref( fileWriter ), it->second,
                                                  mesh, boost::ref( progress ), boost::ref( errorQueue ) ) );
            }

            // write face channels
            it = faceChannels.begin();
            itEnd = faceChannels.end();
            for( ; it != itEnd; ++it ) {
                threadPool.schedule( boost::bind( process_face_channel, boost::ref( fileWriter ), it->second,
                                                  boost::ref( progress ), boost::ref( errorQueue ) ) );
            }
        } // scope for threadPool

        // All threads are joined now,
        // because the thread pool is destroyed (out of scope)

        if( errorQueue.size() == 0 ) {
            // only close the file if there was no error
            fileWriter.close();
        }

        done = true;
    } catch( std::exception& e ) {
        errMsg = e.what();
    } catch( ... ) {
        errMsg = "write_xmesh_mesh_file_worker Error: An unknown exception occurred.";
    }
    if( !done ) {
        errorQueue.enter( errMsg );
    }
}

} // anonymous namespace

void write_xmesh_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                            frantic::logging::progress_logger& progress, std::size_t threadCount ) {
    xmesh_metadata metadata;
    write_xmesh_mesh_file( path, mesh, metadata, progress, threadCount );
}

void write_xmesh_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                            const xmesh_metadata& metadata, frantic::logging::progress_logger& progress,
                            std::size_t threadCount ) {
    if( !mesh ) {
        throw std::runtime_error( "write_xmesh_mesh_file Error: mesh is NULL" );
    }

    progress_info_wrapper threadProgress( progress ); // thread safe
    threadProgress.progressTotalCount = get_progress_total_count( mesh );
    threadProgress.progressElementsPassed = 0;

    error_queue_t errorQueue;

#ifndef FRANTIC_DISABLE_THREADS
    boost::thread worker( boost::bind( write_xmesh_mesh_file_worker, path, mesh, metadata, boost::ref( threadProgress ),
                                       threadCount, boost::ref( errorQueue ) ) );
    while( worker.joinable() ) {
        worker.try_join_for( boost::chrono::milliseconds( 100 ) );
        try {
            progress.update_progress( threadProgress.progressElementsPassed, threadProgress.progressTotalCount );
        } catch( frantic::logging::progress_cancel_exception& ) {
            // tell threads to abort because user canceled
            threadProgress.abortRequested = true;
        }
    }
#else
    write_xmesh_mesh_file_worker( path, mesh, metadata, threadProgress, threadCount, boost::ref( errorQueue ) );
#endif

    // this will throw progress_cancel_exception if user canceled
    progress.update_progress( 100.f );

    // Check for thread errors
    if( errorQueue.size() ) {
        std::string errMsg;
        bool success = errorQueue.leave( errMsg );
        if( success ) {
            throw std::runtime_error( errMsg );
        } else {
            throw std::runtime_error( "An unknown error occurred while attempting to retrieve worker error message." );
        }
    }
}

namespace {
bool has_all_zero_z( const mesh_channel_cvt<frantic::graphics::vector3f>& acc ) {
    for( std::size_t i = 0, ie = acc.get_num_elements(); i != ie; ++i ) {
        frantic::graphics::vector3f v = acc.get_value( i );
        if( v.z != 0 ) {
            return false;
        }
    }
    return true;
}
} // anonymous namespace

void write_obj_mesh_file( const frantic::tstring& destFile, const mesh_interface_ptr& mesh,
                          frantic::logging::progress_logger& progress ) {
    std::size_t progressSkip = 256; // number of lines to write between progress updates
    using frantic::locale::fprintf_c;

    if( !mesh )
        throw std::runtime_error( "write_obj_mesh_file: mesh is NULL. Failed to write mesh to \"" +
                                  frantic::strings::to_string( destFile ) + "\"." );

    temp_file_manager fileManager( destFile );
    frantic::files::file_ptr fout( frantic::files::tfopen( fileManager.get_path().c_str(), _T("w") ) );

    if( !fout )
        throw std::runtime_error( "write_obj_mesh_file: Could not open file \"" +
                                  frantic::strings::to_string( destFile ) + "\" for writing." );

    fprintf_c( fout, "# OBJ polygon mesh file generated from Thinkbox Software mesh_interface object.\n\n" );

    const mesh_interface::mesh_channel_map& vertexChannels = mesh->get_vertex_channels();
    bool hasTverts = mesh->has_vertex_channel( _T("TextureCoord") );
    bool hasNormals = mesh->has_vertex_channel( _T("Normal") );

    int progressSectionsCount = 2;
    if( hasTverts ) {
        ++progressSectionsCount;
    }
    if( hasNormals ) {
        ++progressSectionsCount;
    }
    int progressSectionsPassed = 0;

    // write vertices
    fprintf_c( fout, "# %d verts\n", int( mesh->get_num_verts() ) );
    { // scope for the vertices progress_logger_subinterval_tracker
        float progressSubintervalStart = ( (float)progressSectionsPassed / (float)progressSectionsCount ) * 100.0f;
        float progressSubintervalEnd =
            ( (float)( progressSectionsPassed + 1 ) / (float)progressSectionsCount ) * 100.0f;
        frantic::logging::progress_logger_subinterval_tracker progressSubinterval( progress, progressSubintervalStart,
                                                                                   progressSubintervalEnd );

        for( std::size_t i = 0, ie = mesh->get_num_verts(); i != ie; ++i ) {
            frantic::graphics::vector3f vert;
            mesh->get_vert( i, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
            fprintf_c( fout, "v %g %g %g\n", vert.x, vert.y, vert.z );

            if( i % progressSkip == 0 )
                progress.update_progress( i, ie );
        }
    }
    ++progressSectionsPassed;

    // write texture vertices
    mesh_channel_cvt<frantic::graphics::vector3f> tverts;
    if( hasTverts ) {
        tverts = mesh_channel_cvt<frantic::graphics::vector3f>( vertexChannels, _T("TextureCoord") );
        fprintf_c( fout, "\n# %d texture verts\n", int( tverts.get_num_elements() ) );

        float progressSubintervalStart = ( (float)progressSectionsPassed / (float)progressSectionsCount ) * 100.0f;
        float progressSubintervalEnd =
            ( (float)( progressSectionsPassed + 1 ) / (float)progressSectionsCount ) * 100.0f;
        frantic::logging::progress_logger_subinterval_tracker progressSubinterval( progress, progressSubintervalStart,
                                                                                   progressSubintervalEnd );

        if( has_all_zero_z( tverts ) ) {
            for( std::size_t i = 0, ie = tverts.get_num_elements(); i != ie; ++i ) {
                frantic::graphics::vector3f tvert = tverts.get_value( i );
                fprintf_c( fout, "vt %g %g\n", tvert.x, tvert.y );

                if( i % progressSkip == 0 )
                    progress.update_progress( i, ie );
            }
        } else {
            for( std::size_t i = 0, ie = tverts.get_num_elements(); i != ie; ++i ) {
                frantic::graphics::vector3f tvert = tverts.get_value( i );
                fprintf_c( fout, "vt %g %g %g\n", tvert.x, tvert.y, tvert.z );

                if( i % progressSkip == 0 )
                    progress.update_progress( i, ie );
            }
        }
        ++progressSectionsPassed;
    }

    // write normals
    mesh_channel_cvt<frantic::graphics::vector3f> normals;
    if( hasNormals ) {
        normals = mesh_channel_cvt<frantic::graphics::vector3f>( vertexChannels, _T("Normal") );
        fprintf_c( fout, "\n# %d normals\n", int( normals.get_num_elements() ) );

        float progressSubintervalStart = ( (float)progressSectionsPassed / (float)progressSectionsCount ) * 100.0f;
        float progressSubintervalEnd =
            ( (float)( progressSectionsPassed + 1 ) / (float)progressSectionsCount ) * 100.0f;
        frantic::logging::progress_logger_subinterval_tracker progressSubinterval( progress, progressSubintervalStart,
                                                                                   progressSubintervalEnd );

        for( std::size_t i = 0, ie = normals.get_num_elements(); i != ie; ++i ) {
            frantic::graphics::vector3f normal = normals.get_value( i );
            fprintf_c( fout, "vn %g %g %g\n", normal.x, normal.y, normal.z );
            if( i % progressSkip == 0 )
                progress.update_progress( i, ie );
        }
        ++progressSectionsPassed;
    }

    // write faces
    fprintf_c( fout, "\ng geometry\n\n# %d faces\n", int( mesh->get_num_faces() ) );
    { // scope for the faces progress_logger_subinterval_tracker
        float progressSubintervalStart = ( (float)progressSectionsPassed / (float)progressSectionsCount ) * 100.0f;
        float progressSubintervalEnd =
            ( (float)( progressSectionsPassed + 1 ) / (float)progressSectionsCount ) * 100.0f;
        frantic::logging::progress_logger_subinterval_tracker progressSubinterval( progress, progressSubintervalStart,
                                                                                   progressSubintervalEnd );

        for( std::size_t i = 0, ie = mesh->get_num_faces(); i != ie; ++i ) {
            // note: compensating for 1-based indexing in OBJ format by adding 1
            std::size_t geomFace = mesh->get_num_face_verts( i );

            // Initially, we just exported faces when a channel had custom faces, but the 3ds max .obj loader
            // just ignores the texture coordinates in that case.  Thus, now we export faces always.
            if( tverts.is_valid() ) {
                if( normals.is_valid() ) {
                    fprintf_c( fout, "f" );
                    for( std::size_t igeom = 0; igeom != geomFace; ++igeom ) {
                        fprintf_c( fout, " %d/%d/%d", int( mesh->get_face_vert_index( i, igeom ) + 1 ),
                                   int( tverts.get_fv_index( i, igeom ) + 1 ),
                                   int( normals.get_fv_index( i, igeom ) + 1 ) );
                    }
                    fprintf_c( fout, "\n" );
                } else {
                    fprintf_c( fout, "f" );
                    for( std::size_t igeom = 0; igeom != geomFace; ++igeom ) {
                        fprintf( fout, " %u/%u", unsigned( mesh->get_face_vert_index( i, igeom ) + 1 ),
                                 unsigned( tverts.get_fv_index( i, igeom ) + 1 ) );
                    }
                    fprintf_c( fout, "\n" );
                }
            } else {
                if( normals.is_valid() ) {
                    fprintf_c( fout, "f" );
                    for( std::size_t igeom = 0; igeom != geomFace; ++igeom ) {
                        fprintf_c( fout, " %d//%d", int( mesh->get_face_vert_index( i, igeom ) + 1 ),
                                   int( normals.get_fv_index( i, igeom ) + 1 ) );
                    }
                    fprintf_c( fout, "\n" );
                } else {
                    fprintf_c( fout, "f" );
                    for( std::size_t igeom = 0; igeom != geomFace; ++igeom ) {
                        fprintf_c( fout, " %d", int( mesh->get_face_vert_index( i, igeom ) + 1 ) );
                    }
                    fprintf_c( fout, "\n" );
                }
            }
            if( i % progressSkip == 0 )
                progress.update_progress( i, ie );
        }
    }
    fprintf_c( fout, "\ng\n" );

    fout.close();
    fileManager.close();
}

void write_ascii_stl_mesh_file( const frantic::tstring& filename, const mesh_interface_ptr& mesh,
                                frantic::logging::progress_logger& progress ) {
    std::size_t progressSkip = 256; // number of lines to write between progress updates
    if( !mesh )
        throw std::runtime_error( "write_ascii_stl_mesh_file: mesh is NULL. Failed to write mesh to \"" +
                                  frantic::strings::to_string( filename ) + "\"." );

    temp_file_manager fileManager( filename );
    std::ofstream fout( fileManager.get_path().c_str() );

    if( !fout )
        throw std::runtime_error( "write_ascii_stl_mesh_file: Could not open file \"" +
                                  frantic::strings::to_string( filename ) + "\" for writing." );

    const mesh_interface::mesh_channel_map& faceChannels = mesh->get_face_channels();

    fout.precision( 6 );

    fout << "solid \n";

    mesh_channel_cvt<frantic::graphics::vector3f> normalAcc;
    if( mesh->has_face_channel( _T("Normal") ) ) {
        normalAcc = mesh_channel_cvt<frantic::graphics::vector3f>( faceChannels, _T("Normal") );
    }

    for( std::size_t i = 0, ie = mesh->get_num_faces(); i != ie; ++i ) {
        // get 3 vertices based on the face
        std::size_t geomFace = mesh->get_num_face_verts( i );
        int cornerCount = 0;
        frantic::graphics::vector3f faceVertices[3];
        for( std::size_t igeom = 0; igeom != geomFace; ++igeom ) {
            if( cornerCount >= 3 )
                throw std::runtime_error(
                    "write_ascii_stl_mesh_file: mesh has an invalid face(more than 3 vertices) for STL "
                    "file format. Failed to write mesh to \"" +
                    frantic::strings::to_string( filename ) + "\"." );
            std::size_t temp = mesh->get_face_vert_index( i, igeom );
            mesh->get_vert( temp, *reinterpret_cast<float( * )[3]>( &faceVertices[cornerCount][0] ) );
            ++cornerCount;
        }
        if( cornerCount < 3 )
            throw std::runtime_error(
                "write_ascii_stl_mesh_file: mesh has an invalid face(less than 3 vertices) for STL "
                "file format. Failed to write mesh to \"" +
                frantic::strings::to_string( filename ) + "\"." );

        frantic::graphics::vector3f normal;

        if( normalAcc.is_valid() ) {
            normal = normalAcc.get_value( i ); // in for loop for faces already
        } else {
            normal = frantic::graphics::triangle_normal( faceVertices[0], faceVertices[1], faceVertices[2] );
        }

        fout << "\tfacet normal " << std::scientific << normal.x << " " << normal.y << " " << normal.z << "\n";
        fout << "\t\touter loop\n";
        for( int corner = 0; corner < 3; ++corner ) {
            fout << "\t\t\tvertex ";
            fout << std::scientific << faceVertices[corner].x << " " << faceVertices[corner].y << " "
                 << faceVertices[corner].z << "\n";
        }
        fout << "\t\tendloop\n";
        fout << "\tendfacet\n";

        if( i % progressSkip == 0 )
            progress.update_progress( i, ie );
    }

    fout << "endsolid";

    fout.close();
    fileManager.close();
}

namespace {
class binary_stl_writer {
  public:
    binary_stl_writer( const frantic::tstring& filename )
        : m_file( frantic::files::tfopen( filename.c_str(), _T("wb") ) )
        , m_filename( filename ) {
        if( !m_file )
            throw std::runtime_error( "binary_stl_writer: Could not open file \"" +
                                      frantic::strings::to_string( filename ) + "\" for writing." );
    }

    void write( char* data, std::size_t size ) {
        if( !size ) {
            return;
        }

        assert( data );
        assert( m_file );

        const std::size_t writeCount = fwrite( data, size, 1, m_file );
        if( writeCount != 1 ) {
            throw std::runtime_error( "binary_stl_writer: Error writing to file \"" +
                                      frantic::strings::to_string( m_filename ) + "\"." );
        }
    }

    void close() { m_file.close(); }

  private:
    binary_stl_writer(); // not implemented

    frantic::files::file_ptr m_file;
    frantic::tstring m_filename;
};
} // anonymous namespace

void write_binary_stl_mesh_file( const frantic::tstring& filename, const mesh_interface_ptr& mesh, bool isSolidView,
                                 frantic::logging::progress_logger& progress ) {
    std::size_t progressSkip = 256; // number of lines to write between progress updates
    if( !mesh )
        throw std::runtime_error( "write_binary_stl_mesh_file: mesh is NULL. Failed to write mesh to \"" +
                                  frantic::strings::to_string( filename ) + "\"." );

    temp_file_manager fileManager( filename );
    binary_stl_writer writer( fileManager.get_path() );

    const mesh_interface::mesh_channel_map& vertexChannels = mesh->get_vertex_channels();
    const mesh_interface::mesh_channel_map& faceChannels = mesh->get_face_channels();

    std::vector<char> header( 80 );

    mesh_channel_cvt<frantic::graphics::vector3f> colorAcc;
    if( mesh->has_vertex_channel( _T("Color") ) ) {
        colorAcc = mesh_channel_cvt<frantic::graphics::vector3f>( vertexChannels, _T("Color") );

        // magics header string
        if( !isSolidView ) {
            std::string colorString =
                "COLOR=\0\0\0\0,MATERIAL=\0\0\0\0\0\0\0\0\0\0\0\0"; // whole object doesnt have color so just zero it
            memcpy( &header[0], &colorString[0], colorString.size() );
        }
    }
    writer.write( &header[0], 80 );

    mesh_channel_cvt<frantic::graphics::vector3f> normalAcc;
    if( mesh->has_face_channel( _T("Normal") ) ) {
        normalAcc = mesh_channel_cvt<frantic::graphics::vector3f>( faceChannels, _T("Normal") );
    }

    const std::size_t vertexCount = mesh->get_num_verts();
    const std::size_t faceCount = mesh->get_num_faces();

    const boost::int32_t facets = static_cast<boost::int32_t>( faceCount );

    writer.write( (char*)&facets, 4 );

    std::vector<frantic::graphics::vector3f> faceVertexColors;

    for( std::size_t i = 0, ie = faceCount; i != ie; ++i ) {
        // get 3 vertices based on the face
        std::size_t geomFace = mesh->get_num_face_verts( i );
        size_t cornerCount = 0;
        frantic::graphics::vector3f faceVertices[3];
        faceVertexColors.clear();
        for( std::size_t igeom = 0; igeom != geomFace; ++igeom ) {
            if( cornerCount >= 3 )
                throw std::runtime_error(
                    "write_binary_stl_mesh_file: mesh has an invalid face(more than 3 vertices) for STL "
                    "file format. Failed to write mesh to \"" +
                    frantic::strings::to_string( filename ) + "\"." );
            std::size_t temp = mesh->get_face_vert_index( i, igeom );
            mesh->get_vert( temp, *reinterpret_cast<float( * )[3]>( &faceVertices[cornerCount][0] ) );
            if( colorAcc.is_valid() ) {
                if( colorAcc.get_num_elements() == vertexCount ) {
                    std::size_t temp2 = colorAcc.get_fv_index( i, igeom );
                    faceVertexColors.push_back( colorAcc.get_value( temp2 ) );
                }
            }
            ++cornerCount;
        }
        if( cornerCount < 3 )
            throw std::runtime_error(
                "write_binary_stl_mesh_file: mesh has an invalid face(less than 3 vertices) for STL "
                "file format. Failed to write mesh to \"" +
                frantic::strings::to_string( filename ) + "\"." );

        frantic::graphics::vector3f normal;

        if( normalAcc.is_valid() ) {
            normal = normalAcc.get_value( i ); // in for loop for faces already
        } else {
            normal = frantic::graphics::triangle_normal( faceVertices[0], faceVertices[1], faceVertices[2] );
        }

        writer.write( (char*)&normal.x, 4 );
        writer.write( (char*)&normal.y, 4 );
        writer.write( (char*)&normal.z, 4 );

        // write the 3 vertices
        for( int corner = 0; corner < 3; ++corner ) {
            writer.write( (char*)&faceVertices[corner].x, 4 );
            writer.write( (char*)&faceVertices[corner].y, 4 );
            writer.write( (char*)&faceVertices[corner].z, 4 );
        }

        boost::uint16_t attribute = 0;
        // check if want to save colors
        // if yes check whether solidview or magics
        // if no color channel set not valid or per-object color
        if( colorAcc.is_valid() && faceVertexColors.size() == cornerCount ) {
            // average out vertices colors
            frantic::graphics::vector3f faceColor;
            for( size_t corner = 0; corner < cornerCount; ++corner ) {
                faceColor += faceVertexColors[corner];
            }
            faceColor /= static_cast<float>( cornerCount );
            boost::uint8_t red = (boost::uint8_t)frantic::math::round( faceColor.x * 31 );
            boost::uint8_t green = (boost::uint8_t)frantic::math::round( faceColor.y * 31 );
            boost::uint8_t blue = (boost::uint8_t)frantic::math::round( faceColor.z * 31 );
            // STL is assumed little endian so it is best to assume lsb bit numbering
            // Solidview
            if( isSolidView ) {
                attribute = 1;
                attribute <<= 5; // left shift the bits of attribute
                attribute |= red;
                attribute <<= 5;    // left shift the bits of attribute
                attribute |= green; // or attribute with green
                attribute <<= 5;
                attribute |= blue;
            } else {
                // Magics
                attribute = blue;
                attribute <<= 5;    // left shift the bits of attribute
                attribute |= green; // or attribute with green
                attribute <<= 5;
                attribute |= red;
            }
        } else {
            attribute = 0; // no color channel just set it as 0
        }

        writer.write( (char*)&attribute, 2 );

        if( i % progressSkip == 0 )
            progress.update_progress( i, ie );
    }
    writer.close();
    fileManager.close();
}

namespace {
// RAII wrapper for ply
// \note don't use the regular ply_create and ply_close when using this wrapper
class ply_file {
    p_ply m_ply;
    bool m_isClosed;

  public:
    ply_file( const frantic::tchar* name, e_ply_storage_mode storage_mode, p_ply_error_cb error_cb, long idata,
              void* pdata )
        : m_isClosed( false ) {
        m_ply = ply_create( name, storage_mode, error_cb, idata, pdata );
    }

    ~ply_file() {
        if( !m_isClosed )
            close();
    }

    operator p_ply() { return m_ply; }

    operator const p_ply() const { return m_ply; }

    int close() {
        if( m_isClosed )
            return 0;
        m_isClosed = true;
        return ply_close( m_ply );
    }
};
} // namespace

void write_ply_mesh_file( const frantic::tstring& filename, const mesh_interface_ptr& mesh,
                          frantic::logging::progress_logger& progress ) {
    if( !mesh )
        throw std::runtime_error( "write_ply_mesh_file: mesh is NULL. Failed to write mesh to \"" +
                                  frantic::strings::to_string( filename ) + "\"." );

    temp_file_manager fileManager( filename );
    // Write the file in little endian binary
    ply_file ply( fileManager.get_path().c_str(), PLY_LITTLE_ENDIAN, NULL, 0, NULL );

    if( !ply )
        throw std::runtime_error( "write_ply_mesh_file: Could not open file \"" +
                                  frantic::strings::to_string( filename ) + "\" for writing." );

    mesh_channel_cvt<frantic::graphics::vector3f> colors;
    mesh_channel_cvt<frantic::graphics::vector3f> textures;
    mesh_channel_cvt<frantic::graphics::vector3f> normals;

    const mesh_interface::mesh_channel_map& vertexChannels = mesh->get_vertex_channels();

    typedef long index_t;

    if( mesh->get_num_verts() > std::numeric_limits<index_t>::max() ) {
        throw std::runtime_error( "write_ply_mesh_file: While writing \"" + frantic::strings::to_string( filename ) +
                                  "\": "
                                  "too many vertices in mesh (" +
                                  boost::lexical_cast<std::string>( mesh->get_num_verts() ) + ")" );
    }

    if( mesh->get_num_faces() > std::numeric_limits<index_t>::max() ) {
        throw std::runtime_error( "write_ply_mesh_file: While writing \"" + frantic::strings::to_string( filename ) +
                                  "\": "
                                  "too many faces in mesh (" +
                                  boost::lexical_cast<std::string>( mesh->get_num_faces() ) + ")" );
    }

    size_t vertexCount = mesh->get_num_verts();
    size_t faceCount = boost::numeric_cast<index_t>( mesh->get_num_faces() );

    // Add comment
    if( !ply_add_comment( ply, "PLY polygon mesh file generated from Thinkbox Software mesh_interface object" ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write comment to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );

    // Add vertex element
    if( !ply_add_element( ply, "vertex", static_cast<index_t>( vertexCount ) ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write element vertex to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );

    // Add the vertex properties
    if( !ply_add_scalar_property( ply, "x", PLY_FLOAT ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write property x to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );
    if( !ply_add_scalar_property( ply, "y", PLY_FLOAT ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write property y to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );
    if( !ply_add_scalar_property( ply, "z", PLY_FLOAT ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write property z to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );

    // Add the vertex color properties if 'Color' channel exists
    if( mesh->has_vertex_channel( _T("Color") ) ) {
        colors = mesh_channel_cvt<frantic::graphics::vector3f>( vertexChannels, _T("Color") );
        if( colors.get_num_elements() == vertexCount ) {
            if( !ply_add_scalar_property( ply, "red", PLY_UCHAR ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property red to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
            if( !ply_add_scalar_property( ply, "green", PLY_UCHAR ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property green to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
            if( !ply_add_scalar_property( ply, "blue", PLY_UCHAR ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property blue to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
        }
    }

    // Add the vertex texture properties if 'TextureCoord' channel exists
    if( mesh->has_vertex_channel( _T("TextureCoord") ) ) {
        textures = mesh_channel_cvt<frantic::graphics::vector3f>( vertexChannels, _T("TextureCoord") );
        if( textures.get_num_elements() == vertexCount ) {
            if( !ply_add_scalar_property( ply, "u", PLY_FLOAT ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property u to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
            if( !ply_add_scalar_property( ply, "v", PLY_FLOAT ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property v to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
            if( !ply_add_scalar_property( ply, "t", PLY_FLOAT ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property t to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
        }
    }

    // Add the vertex normal properties if 'Normal' channel exists
    if( mesh->has_vertex_channel( _T("Normal") ) ) {
        normals = mesh_channel_cvt<frantic::graphics::vector3f>( vertexChannels, _T("Normal") );
        if( normals.get_num_elements() == vertexCount ) {
            if( !ply_add_scalar_property( ply, "nx", PLY_FLOAT ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property nx to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
            if( !ply_add_scalar_property( ply, "ny", PLY_FLOAT ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property ny to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
            if( !ply_add_scalar_property( ply, "nz", PLY_FLOAT ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write property nz to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
        }
    }

    // Add face element
    if( !ply_add_element( ply, "face", static_cast<index_t>( faceCount ) ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write element face to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );
    // Add property list containing vertex indices
    if( !ply_add_list_property( ply, "vertex_indices", PLY_UCHAR, PLY_INT ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write property list vertex_indices to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );
    // Write the header information to the file
    if( !ply_write_header( ply ) )
        throw std::runtime_error( "write_ply_mesh_file: Unable to write header to file \"" +
                                  frantic::strings::to_string( filename ) + "\"" );

    // Write vertices to the file
    frantic::logging::progress_logger_subinterval_tracker progressSubinterval( progress, 0.0f, 50.0f );
    for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
        frantic::graphics::vector3f vert;
        mesh->get_vert( vertexIndex, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
        if( !ply_write( ply, vert.x ) )
            throw std::runtime_error( "write_ply_mesh_file: Unable to write data for x to file \"" +
                                      frantic::strings::to_string( filename ) + "\"" );
        if( !ply_write( ply, vert.y ) )
            throw std::runtime_error( "write_ply_mesh_file: Unable to write data for y to file \"" +
                                      frantic::strings::to_string( filename ) + "\"" );
        if( !ply_write( ply, vert.z ) )
            throw std::runtime_error( "write_ply_mesh_file: Unable to write data for z to file \"" +
                                      frantic::strings::to_string( filename ) + "\"" );

        if( colors.is_valid() ) {
            if( colors.get_num_elements() == vertexCount ) {
                frantic::graphics::vector3f color = colors.get_value( vertexIndex );
                for( int i = 0; i < 3; ++i ) {
                    color[i] = 255 * frantic::math::clamp<float>( color[i], 0, 1 );
                }
                if( !ply_write( ply, color.x ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for red to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
                if( !ply_write( ply, color.y ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for green to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
                if( !ply_write( ply, color.z ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for blue to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
            }
        }

        if( textures.is_valid() ) {
            if( textures.get_num_elements() == vertexCount ) {
                frantic::graphics::vector3f texture = textures.get_value( vertexIndex );
                if( !ply_write( ply, texture.x ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for u to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
                if( !ply_write( ply, texture.y ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for v to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
                if( !ply_write( ply, texture.z ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for t to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
            }
        }

        if( normals.is_valid() ) {
            if( normals.get_num_elements() == vertexCount ) {
                frantic::graphics::vector3f normal = normals.get_value( vertexIndex );
                if( !ply_write( ply, normal.x ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for nx to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
                if( !ply_write( ply, normal.y ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for ny to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
                if( !ply_write( ply, normal.z ) )
                    throw std::runtime_error( "write_ply_mesh_file: Unable to write data for nz to file \"" +
                                              frantic::strings::to_string( filename ) + "\"" );
            }
        }

        if( ( vertexIndex & 0xffu ) == 0 )
            progress.update_progress( vertexIndex, vertexCount );
    }

    // Write faces to the file
    progressSubinterval.reset( 50.0f, 100.0f );
    std::vector<index_t> face;
    for( size_t i = 0; i < faceCount; ++i ) {
        std::size_t cornerCount = mesh->get_num_face_verts( i );
        face.clear();
        for( std::size_t corner = 0; corner != cornerCount; ++corner ) {
            face.push_back( (index_t)mesh->get_face_vert_index( i, corner ) );
        }
        // Write the length of the list into the file
        if( !ply_write( ply, (double)face.size() ) )
            throw std::runtime_error( "write_ply_mesh_file: Unable to write vertex_indices length to file \"" +
                                      frantic::strings::to_string( filename ) + "\"" );
        // Write the face into the file
        for( std::size_t corner = 0; corner < cornerCount; ++corner ) {
            if( !ply_write( ply, face[corner] ) )
                throw std::runtime_error( "write_ply_mesh_file: Unable to write vertex_indices values to file \"" +
                                          frantic::strings::to_string( filename ) + "\"" );
        }
        if( ( i & 0xff ) == 0 )
            progress.update_progress( i, faceCount );
    }
    if( !ply.close() )
        throw std::runtime_error( "write_ply_mesh_file: Could not close file \"" +
                                  frantic::strings::to_string( filename ) + "\" after writing." );

    fileManager.close();
}

} // namespace geometry
} // namespace frantic
