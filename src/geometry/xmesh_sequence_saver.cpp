// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/files/compression_stream.hpp>
#include <frantic/files/files.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/xmesh_reader.hpp>
#include <frantic/geometry/xmesh_sequence_saver.hpp>
#include <frantic/geometry/xmesh_writer.hpp>
#include <frantic/misc/string_functions.hpp>
#include <frantic/strings/to_string_classic.hpp>
#include <frantic/strings/utf8.hpp>
#include <frantic/threads/synchronizedqueue.hpp>
#include <frantic/tinyxml/frantic_tinyxml_graphics.hpp>
#include <frantic/tinyxml/frantic_tinyxml_utility.hpp>

#include <tinyxml2.h>

#include <boost/bind.hpp>
#include <boost/cstdint.hpp>
#ifndef FRANTIC_DISABLE_THREADS
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#endif

#include <fstream>

using namespace std;
using namespace frantic;
using namespace frantic::geometry;
using namespace frantic::files;
using namespace frantic::channels;

namespace frantic {
namespace geometry {

namespace {
// TODO: It's not obvious how char string channel names will be encoded.
// I think it should be UTF-8, but we have code that can produce different
// encodings.  For now, I assume that it's UTF-8 if it's valid UTF-8.

// TODO: These functions are shared with xmesh_writer.  Move them somewhere
// common.

inline std::wstring wstring_from_channel_name( const std::string& s ) {
    if( frantic::strings::is_valid_utf8( s ) ) {
        return frantic::strings::wstring_from_utf8( s );
    } else {
        return frantic::strings::wstring_from_utf8( frantic::strings::to_utf8( s ) );
    }
}

inline const std::wstring& wstring_from_channel_name( const std::wstring& s ) { return s; }

inline std::string utf8_from_channel_name( const std::string& s ) {
    if( frantic::strings::is_valid_utf8( s ) ) {
        return s;
    } else {
        return frantic::strings::to_utf8( s );
    }
}

inline std::string utf8_from_channel_name( const std::wstring& s ) { return frantic::strings::to_utf8( s ); }

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

// Wrapper for invoking write_xmesh_array_file in a worker thread.
// Catches exceptions and pushes them onto a queue.
inline void write_xmesh_array_file( const boost::filesystem::path::string_type& filename, const void* data,
                                    const std::string& dataType, std::size_t dataSize, std::size_t numData,
                                    int zlibCompressionLevel, error_queue_t& errorQueue,
                                    boost::filesystem::path::string_type& filenameOut ) {
    bool done = false;
    std::string errMsg;
    try {
        filenameOut.clear();
        frantic::geometry::write_xmesh_array_file( filename, data, dataType, dataSize, numData, zlibCompressionLevel );
        filenameOut.assign( boost::filesystem::path( filename ).filename().native() );
        done = true;
    } catch( std::exception& e ) {
        errMsg = e.what();
    } catch( ... ) {
        errMsg = "write_xmesh_array_file Error: An unknown exception occurred.";
    }
    if( !done ) {
        errorQueue.enter( errMsg );
    }
}

size_t element_size( data_type_t dataType, size_t arity ) {
    return frantic::channels::sizeof_channel_data_type( dataType ) * arity;
}

size_t element_size( const xmesh_sequence_saver::xmesh_cache_channel& cacheChannel ) {
    return element_size( cacheChannel.m_dataType, cacheChannel.m_arity );
}

size_t num_elements( const xmesh_sequence_saver::xmesh_cache_channel& cacheChannel ) {
    return cacheChannel.m_data.size() / element_size( cacheChannel.m_dataType, cacheChannel.m_arity );
}

frantic::graphics::boundbox3f get_bound_box( const xmesh_sequence_saver::xmesh_cache_channel& cacheChannel ) {
    const size_t elementSize = element_size( cacheChannel );
    const size_t numElements = num_elements( cacheChannel );

    frantic::graphics::boundbox3f result;

    for( size_t i = 0; i < numElements; ++i ) {
        const vector3f& position = *reinterpret_cast<const vector3f*>( cacheChannel.m_data.ptr_at( i * elementSize ) );
        result += position;
    }

    return result;
}

bool matches_cache_channel( const xmesh_sequence_saver::xmesh_cache_channel& cacheChannel, data_type_t dataType,
                            size_t arity, size_t elementCount, const void* data ) {
    if( cacheChannel.m_dataType != dataType || cacheChannel.m_arity != arity ||
        cacheChannel.m_data.size() / element_size( dataType, arity ) != elementCount )
        return false;

    if( memcmp( cacheChannel.m_data.ptr_at( 0 ), data, cacheChannel.m_data.size() ) != 0 )
        return false;

    return true;
}

bool matches_vertices( const xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel, const trimesh3& mesh ) {
    const void* data = mesh.vertex_count() > 0 ? reinterpret_cast<const void*>( &mesh.get_vertex( 0 ) ) : 0;
    return matches_cache_channel( cacheChannel, data_type_float32, 3, mesh.vertex_count(), data );
}

bool matches_vertex_channel( const xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                             const const_trimesh3_vertex_channel_general_accessor& trimeshChannel ) {
    return matches_cache_channel( cacheChannel, trimeshChannel.data_type(), trimeshChannel.arity(),
                                  trimeshChannel.size(), reinterpret_cast<const void*>( trimeshChannel.data( 0 ) ) );
}

bool matches_vertex_channel( const xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                             const polymesh3_const_vertex_accessor<void>& polymeshChannel ) {
    return matches_cache_channel( cacheChannel, polymeshChannel.get_type(), polymeshChannel.get_arity(),
                                  polymeshChannel.vertex_count(),
                                  reinterpret_cast<const void*>( polymeshChannel.get_vertex( 0 ) ) );
}

bool matches_vertex_channel_custom_faces( const xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                                          size_t indiceCount, const boost::int32_t* customFacesData ) {
    if( indiceCount == cacheChannel.m_customFaces.size() ) {
        if( cacheChannel.m_customFaces.size() > 0 )
            return ( memcmp( &cacheChannel.m_customFaces[0], customFacesData,
                             sizeof( boost::int32_t ) * cacheChannel.m_customFaces.size() ) == 0 );
        else
            return true;
    } else
        return false;
}

bool matches_indices( const xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel, const trimesh3& mesh ) {
    const boost::int32_t* faceIndices = mesh.face_count() > 0 ? &( mesh.faces_ref()[0][0] ) : 0;
    return matches_vertex_channel_custom_faces( cacheChannel, mesh.face_count() * 3, faceIndices );
}

bool matches_vertex_channel_custom_faces( const xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                                          const const_trimesh3_vertex_channel_general_accessor& trimeshChannel ) {
    if( trimeshChannel.has_custom_faces() ) {
        const boost::int32_t* faceIndices = trimeshChannel.face_count() > 0 ? &trimeshChannel.face( 0 )[0] : 0;
        return matches_vertex_channel_custom_faces( cacheChannel, trimeshChannel.face_count() * 3, faceIndices );
    } else {
        return cacheChannel.m_customFaces.size() == 0;
    }
}

bool matches_vertex_channel_custom_faces( const xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                                          const polymesh3_const_vertex_accessor<void>& polymeshChannel ) {
    if( polymeshChannel.has_custom_faces() && polymeshChannel.face_count() > 0 ) {
        size_t indicesCount =
            polymeshChannel.get_face( polymeshChannel.face_count() - 1 ).second - polymeshChannel.get_face( 0 ).first;
        return matches_vertex_channel_custom_faces( cacheChannel, indicesCount, polymeshChannel.get_face( 0 ).first );
    } else
        return cacheChannel.m_customFaces.size() == 0;
}

bool matches_face_channel( const xmesh_sequence_saver::xmesh_cache_face_channel& cacheChannel,
                           const const_trimesh3_face_channel_general_accessor& trimeshChannel ) {
    return matches_cache_channel( cacheChannel, trimeshChannel.data_type(), trimeshChannel.arity(),
                                  trimeshChannel.size(), reinterpret_cast<const void*>( trimeshChannel.data( 0 ) ) );
}

bool matches_face_channel( const xmesh_sequence_saver::xmesh_cache_face_channel& cacheChannel,
                           const polymesh3_const_face_accessor<void>& polymeshChannel ) {
    return matches_cache_channel( cacheChannel, polymeshChannel.get_type(), polymeshChannel.get_arity(),
                                  polymeshChannel.face_count(),
                                  reinterpret_cast<const void*>( polymeshChannel.get_face( 0 ) ) );
}

void copy_cache_channel( xmesh_sequence_saver::xmesh_cache_channel& cacheChannel, data_type_t dataType, size_t arity,
                         size_t elementCount, const void* data ) {
    cacheChannel.m_dataType = dataType;
    cacheChannel.m_arity = arity;
    cacheChannel.m_data.resize( elementCount * element_size( dataType, arity ) );
    memcpy( cacheChannel.m_data.ptr_at( 0 ), data, cacheChannel.m_data.size() );
}

void copy_vertices( xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel, const trimesh3& mesh ) {
    const void* srcData = mesh.vertex_count() > 0 ? reinterpret_cast<const void*>( &mesh.vertices_ref()[0] ) : 0;
    copy_cache_channel( cacheChannel, data_type_float32, 3, mesh.vertex_count(), srcData );
}

void copy_vertex_channel( xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                          const const_trimesh3_vertex_channel_general_accessor& trimeshChannel ) {
    copy_cache_channel( cacheChannel, trimeshChannel.data_type(), trimeshChannel.arity(), trimeshChannel.size(),
                        reinterpret_cast<const void*>( trimeshChannel.data( 0 ) ) );
}

void copy_vertex_channel( xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                          const polymesh3_const_vertex_accessor<void>& polymeshChannel ) {
    copy_cache_channel( cacheChannel, polymeshChannel.get_type(), polymeshChannel.get_arity(),
                        polymeshChannel.vertex_count(),
                        reinterpret_cast<const void*>( polymeshChannel.get_vertex( 0 ) ) );
}

void copy_vertex_channel_custom_faces( xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                                       size_t indiceCount, const boost::int32_t* customFaces ) {
    cacheChannel.m_customFaces.resize( indiceCount );
    if( cacheChannel.m_customFaces.size() > 0 ) {
        memcpy( &cacheChannel.m_customFaces[0], customFaces, indiceCount * sizeof( boost::int32_t ) );
    }
}

void copy_indices( xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel, const trimesh3& mesh ) {
    const boost::int32_t* faceIndices = mesh.face_count() > 0 ? &mesh.faces_ref()[0][0] : 0;
    copy_vertex_channel_custom_faces( cacheChannel, mesh.face_count() * 3, faceIndices );
}

void copy_vertex_channel_custom_faces( xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                                       const const_trimesh3_vertex_channel_general_accessor& trimeshChannel ) {
    if( trimeshChannel.has_custom_faces() ) {
        const boost::int32_t* faceIndices = trimeshChannel.face_count() > 0 ? &trimeshChannel.face( 0 )[0] : 0;
        copy_vertex_channel_custom_faces( cacheChannel, trimeshChannel.face_count() * 3, faceIndices );
    } else {
        cacheChannel.m_customFaces.resize( 0 );
    }
}

void copy_vertex_channel_custom_faces( xmesh_sequence_saver::xmesh_cache_vertex_channel& cacheChannel,
                                       const polymesh3_const_vertex_accessor<void>& polymeshChannel ) {
    if( polymeshChannel.has_custom_faces() ) {
        size_t indicesCount =
            polymeshChannel.get_face( polymeshChannel.face_count() - 1 ).second - polymeshChannel.get_face( 0 ).first;
        copy_vertex_channel_custom_faces( cacheChannel, indicesCount, polymeshChannel.get_face( 0 ).first );
    } else
        cacheChannel.m_customFaces.resize( 0 );
}

void copy_face_channel( xmesh_sequence_saver::xmesh_cache_face_channel& cacheChannel,
                        const const_trimesh3_face_channel_general_accessor& trimeshChannel ) {
    return copy_cache_channel( cacheChannel, trimeshChannel.data_type(), trimeshChannel.arity(), trimeshChannel.size(),
                               trimeshChannel.data( 0 ) );
}

void copy_face_channel( xmesh_sequence_saver::xmesh_cache_face_channel& cacheChannel,
                        const polymesh3_const_face_accessor<void>& polymeshChannel ) {
    return copy_cache_channel( cacheChannel, polymeshChannel.get_type(), polymeshChannel.get_arity(),
                               polymeshChannel.face_count(), polymeshChannel.get_face( 0 ) );
}

void write_xmesh_data_file( thread_pool& threadPool, const xmesh_sequence_saver::xmesh_cache_channel& cacheChannel,
                            const boost::filesystem::path& filePath, int compressionLevel, error_queue_t& errorQueue,
                            boost::filesystem::path::string_type& outSavedPath ) {
    const void* data = 0;

    if( num_elements( cacheChannel ) > 0 )
        data = static_cast<const void*>( cacheChannel.m_data.ptr_at( 0 ) );

    // HACK: do not use frantic::strings::to_string(
    threadPool.schedule( boost::bind( write_xmesh_array_file, filePath.native(), data,
                                      frantic::strings::to_string( frantic::channels::channel_data_type_str(
                                          cacheChannel.m_arity, cacheChannel.m_dataType ) ),
                                      element_size( cacheChannel ), num_elements( cacheChannel ), compressionLevel,
                                      boost::ref( errorQueue ), boost::ref( outSavedPath ) ) );
}

void write_xmesh_vertex_custom_faces_data_file( thread_pool& threadPool,
                                                xmesh_sequence_saver::xmesh_cache_vertex_channel& vertexChannel,
                                                const boost::filesystem::path& filePath, int compressionLevel,
                                                const std::vector<int>& faceEndOffsets, error_queue_t& errorQueue,
                                                boost::filesystem::path::string_type& outSavedPath ) {

    const void* data = 0;

    vertexChannel.m_serializationCustomFaces.resize( vertexChannel.m_customFaces.size() );

    if( vertexChannel.m_serializationCustomFaces.size() > 0 ) {
        if( faceEndOffsets.size() > 0 ) {
            vertexChannel.m_serializationCustomFaces = vertexChannel.m_customFaces;
            for( size_t i = 0; i < faceEndOffsets.size(); ++i ) {
                size_t actualFaceEnd = faceEndOffsets[i] - 1;
                vertexChannel.m_serializationCustomFaces[actualFaceEnd] =
                    -vertexChannel.m_serializationCustomFaces[actualFaceEnd] - 1;
            }
        } else {
            for( size_t i = 0; i < vertexChannel.m_serializationCustomFaces.size(); ++i ) {
                if( i % 3 == 2 )
                    vertexChannel.m_serializationCustomFaces[i] = -vertexChannel.m_customFaces[i] - 1;
                else
                    vertexChannel.m_serializationCustomFaces[i] = vertexChannel.m_customFaces[i];
            }
        }

        data = &( vertexChannel.m_serializationCustomFaces[0] );
    }

    // HACK: do not use frantic::strings::to_string(
    threadPool.schedule( boost::bind( write_xmesh_array_file, filePath.native(), data, "int32[1]",
                                      sizeof( boost::int32_t ), vertexChannel.m_serializationCustomFaces.size(),
                                      compressionLevel, boost::ref( errorQueue ), boost::ref( outSavedPath ) ) );
}

} // anonymous namespace

// Loads an xmesh file as saved by the xmesh_sequence_saver.  Existence checks are built into all
// the calls to the tinyxml wrapper, which will throw an exception if the requested tag isnt
// present in the xml file.
void load_xmesh_file( const frantic::tstring& srcFile, frantic::geometry::trimesh3& outMesh ) {
    using namespace frantic::channels;
    // using namespace frantic::graphics;

    // NOTE: This was causing all kinds of trouble because often people were passing in pre-intialized meshes
    //       that had the same channels but without the same counts, which isn't caught by many trimesh3
    //       routines. Make sure this mesh gets cleared at the start.
    outMesh.clear();

#if defined( _WIN32 ) || defined( _WIN64 )
    boost::filesystem::path srcPath( frantic::strings::to_wstring( srcFile ) );
#else
    boost::filesystem::path srcPath( srcFile );
#endif
    xmesh_reader fileReader( srcPath );

    const xmesh_vertex_channel& vertChannel = fileReader.get_vertex_channel( _T("verts") );
    std::size_t numVerts = vertChannel.get_vertex_count();
    std::size_t numFaceElements = vertChannel.get_face_count();
    std::size_t numPolygons = numFaceElements / 3;

    // If there are no verts, then there can be no faces. If there are no verts there can be no vertex channels
    // without custom faces, but there are no faces, so there cannot be any vertex channels with custom faces.
    // There can be no face channels either because there are no faces. So this mesh must be empty.
    if( numVerts == 0 )
        return;

    outMesh.vertices_ref().resize( numVerts );
    outMesh.faces_ref().resize( numPolygons );

    float* pVertData = &outMesh.vertices_ref()[0].x;
    int* pFaceData = numPolygons > 0 ? (int*)&outMesh.faces_ref()[0].x : 0;

    fileReader.load_vertex_channel( _T("verts"), (char*)pVertData, data_type_float32, 3, numVerts );
    if( numFaceElements > 0 )
        fileReader.load_vertex_channel_faces( _T("verts"), (char*)pFaceData, numFaceElements );

    // Every third index is negative due to the encoding used in the xmesh file. We now go through and undo this.
    // Also, we must check that the face indices are actually valid relative to the vertex array, or there will
    // probably be an access to invalid memory.
    for( std::size_t i = 0, iEnd = numFaceElements; i < iEnd; i += 3 ) {
        if( pFaceData[i] < 0 || pFaceData[i + 1] < 0 || pFaceData[i + 2] >= 0 )
            throw std::runtime_error( "load_xmesh_file() The file \"" + frantic::strings::to_string( srcFile ) +
                                      "\" was not a triangle mesh" );
        pFaceData[i + 2] = -pFaceData[i + 2] - 1;
        if( pFaceData[i] >= (int)numVerts || pFaceData[i + 1] >= (int)numVerts || pFaceData[i + 2] >= (int)numVerts )
            throw std::runtime_error( "load_xmesh_file() The file \"" + frantic::strings::to_string( srcFile ) +
                                      "\" had invalid geometry face indices" );
    }

    std::vector<frantic::tstring> channelNames;
    fileReader.get_vertex_channel_names( channelNames );

    for( std::size_t i = 0, iEnd = channelNames.size(); i < iEnd; ++i ) {
        const xmesh_vertex_channel& ch = fileReader.get_vertex_channel( channelNames[i] );
        xmesh_reader::channel_type_t chType = ch.get_vertex_type();
        std::size_t chVertexCount = ch.get_vertex_count();
        std::size_t chFaceCount = ch.get_face_count();

        // Check that this channel has the correct correspondence with the geometry channel. If there aren't any
        // custom faces, then this channel needs to have the same number of vertices as the geometry channel. If
        // there are custom faces, then the number of faces needs to match the geometry channel.
        if( chFaceCount > 0 ) {
            if( chFaceCount != numFaceElements )
                throw std::runtime_error( "load_xmesh_file() The channel \"" +
                                          frantic::strings::to_string( channelNames[i] ) + "\" from file \"" +
                                          frantic::strings::to_string( srcFile ) +
                                          "\" had a different number of custom faces from the geometry channel" );
        } else {
            if( chVertexCount != numVerts )
                throw std::runtime_error(
                    "load_xmesh_file() The channel \"" + frantic::strings::to_string( channelNames[i] ) +
                    "\" from file \"" + frantic::strings::to_string( srcFile ) +
                    "\" had a different number of vertices from the geometry channel, and no custom faces" );
        }

        outMesh.add_vertex_channel_raw( channelNames[i], chType.second, chType.first, chVertexCount,
                                        ( chFaceCount > 0 ) );

        frantic::geometry::trimesh3_vertex_channel_general_accessor chAcc =
            outMesh.get_vertex_channel_general_accessor( channelNames[i] );

        fileReader.load_vertex_channel( channelNames[i], chAcc.data( 0 ), chAcc.data_type(), chAcc.arity(),
                                        chAcc.size() );
        if( chFaceCount > 0 ) {
            int* pFaceData = (int*)&chAcc.face( 0 ).x;
            fileReader.load_vertex_channel_faces( channelNames[i], (char*)pFaceData, chFaceCount );

            for( std::size_t i = 0, iEnd = chFaceCount; i < iEnd; i += 3 ) {
                if( pFaceData[i] < 0 || pFaceData[i + 1] < 0 || pFaceData[i + 2] >= 0 )
                    throw std::runtime_error( "load_xmesh_file() The channel \"" +
                                              frantic::strings::to_string( channelNames[i] ) + "\" from file \"" +
                                              frantic::strings::to_string( srcFile ) +
                                              "\" did not have a triangle mapping" );
                pFaceData[i + 2] = -pFaceData[i + 2] - 1;
                if( pFaceData[i] >= (int)chVertexCount || pFaceData[i + 1] >= (int)chVertexCount ||
                    pFaceData[i + 2] >= (int)chVertexCount )
                    throw std::runtime_error( "load_xmesh_file() The channel \"" +
                                              frantic::strings::to_string( channelNames[i] ) + "\" from file \"" +
                                              frantic::strings::to_string( srcFile ) +
                                              "\" had invalid custom face indices" );
            }
        }
    }

    channelNames.clear();
    fileReader.get_face_channel_names( channelNames );

    for( std::size_t i = 0, iEnd = channelNames.size(); i < iEnd; ++i ) {
        const xmesh_face_channel& ch = fileReader.get_face_channel( channelNames[i] );
        xmesh_reader::channel_type_t chType = ch.get_face_type();
        std::size_t chFaceCount = ch.get_face_count();

        if( chFaceCount != numPolygons )
            throw std::runtime_error( "load_xmesh_file() The channel \"" +
                                      frantic::strings::to_string( channelNames[i] ) + "\" from file \"" +
                                      frantic::strings::to_string( srcFile ) +
                                      "\" did not have the same number of polygons as the geometry channel" );

        outMesh.add_face_channel_raw( channelNames[i], chType.second, chType.first );

        frantic::geometry::trimesh3_face_channel_general_accessor chAcc =
            outMesh.get_face_channel_general_accessor( channelNames[i] );

        fileReader.load_face_channel( channelNames[i], chAcc.data( 0 ), chAcc.data_type(), chAcc.arity(),
                                      chAcc.size() );
    }
}

// This loads the mesh data into the memory buffers, comparing it
// with the existing data.  If any mismatches are found, the corresponding
// filename variables are set to the empty string.
void xmesh_sequence_saver::retrieve_and_compare_mesh( const trimesh3& mesh ) {
    // Triangle meshes don't use m_cachedFaceOffsets
    m_cachedFaceOffsets.clear();

    /////
    // Check the mesh vertices
    /////

    // An initial check for vertex equality through size of the vertex arrays
    if( !matches_vertices( m_cachedGeometryChannel, mesh ) ) {
        m_savedVertsFile.clear();
        copy_vertices( m_cachedGeometryChannel, mesh );
    }

    /////
    // Check the mesh faces
    /////

    // Another initial size check for equality
    if( !matches_indices( m_cachedGeometryChannel, mesh ) ) {
        m_savedFacesFile.clear();
        copy_indices( m_cachedGeometryChannel, mesh );
    }

    /////
    // Check all the vertex channels
    /////

    std::vector<frantic::tstring> newVertChannelNames;
    mesh.get_vertex_channel_names( newVertChannelNames );

    // only do the extra channel stuff if there are actually extra channels in the new
    // mesh, otherwise just wipe out the old channel data
    if( !newVertChannelNames.empty() ) {

        for( unsigned channel = 0; channel < newVertChannelNames.size(); ++channel ) {
            retrieve_and_compare_vertex_channel( mesh, newVertChannelNames[channel] );
        }

        std::vector<frantic::tstring> currentChannels;
        // if we have extra channels in the saved mesh that aren't in the new one, get rid of them
        for( vertex_channels_cache_t::iterator it = m_cachedVertexChannels.begin(); it != m_cachedVertexChannels.end();
             ++it ) {
            currentChannels.push_back( it->first );
        }

        for( size_t i = 0; i < currentChannels.size(); ++i ) {
            if( !mesh.has_vertex_channel( currentChannels[i] ) ) {
                m_cachedVertexChannels.erase( currentChannels[i] );
                m_vertChannelVertDataChanged.erase( currentChannels[i] );
                m_vertChannelFaceDataChanged.erase( currentChannels[i] );
            }
        }
    } else {
        m_cachedVertexChannels.clear();
        m_vertChannelVertDataChanged.clear();
        m_vertChannelFaceDataChanged.clear();
    }

    /////
    // Check all the face channels
    /////

    std::vector<frantic::tstring> newFaceChannelNames;
    mesh.get_face_channel_names( newFaceChannelNames );

    // only do the extra channel stuff if there are actually extra channels in the new
    // mesh, otherwise just wipe out the old channel data
    if( !newFaceChannelNames.empty() ) {

        for( unsigned channel = 0; channel < newFaceChannelNames.size(); ++channel ) {
            retrieve_and_compare_face_channel( mesh, newFaceChannelNames[channel] );
        }

        std::vector<frantic::tstring> currentChannels;
        // if we have extra channels in the saved mesh that aren't in the new one, get rid of them
        for( face_channels_cache_t::iterator it = m_cachedFaceChannels.begin(); it != m_cachedFaceChannels.end();
             ++it ) {
            currentChannels.push_back( it->first );
        }

        for( size_t i = 0; i < currentChannels.size(); ++i ) {
            if( !mesh.has_face_channel( currentChannels[i] ) ) {
                m_cachedFaceChannels.erase( currentChannels[i] );
                m_faceChannelDataChanged.erase( currentChannels[i] );
            }
        }
    } else {
        m_cachedFaceChannels.clear();
        m_faceChannelDataChanged.clear();
    }
}

// Same as retrieve_and_compare_mesh, but does it for one map channel of a mesh against the saved mesh
void xmesh_sequence_saver::retrieve_and_compare_vertex_channel( const trimesh3& mesh,
                                                                const frantic::tstring& channelName ) {

    // Grab the vertex channel data
    const_trimesh3_vertex_channel_general_accessor vertChannel =
        mesh.get_vertex_channel_general_accessor( channelName );
    vertex_channels_cache_t::iterator found = m_cachedVertexChannels.find( channelName );

    // first check if the saved mesh even has that channel
    if( found == m_cachedVertexChannels.end() ) {
        // channel will have to be added
        m_vertChannelVertDataChanged[channelName].clear();
        m_vertChannelFaceDataChanged[channelName].clear();
        xmesh_cache_vertex_channel& cachedVertexChannel = m_cachedVertexChannels[channelName];
        copy_vertex_channel( cachedVertexChannel, vertChannel );
        copy_vertex_channel_custom_faces( cachedVertexChannel, vertChannel );
    } else {
        // the channel exists, so grab the corresponding saved vertex channel data too
        xmesh_cache_vertex_channel& cachedVertexChannel = found->second;

        // Check for matching vertex data within the channel, re-write if no match
        if( !matches_vertex_channel( cachedVertexChannel, vertChannel ) ) {
            m_vertChannelVertDataChanged[channelName].clear();
            copy_vertex_channel( cachedVertexChannel, vertChannel );
        }

        // Check for changes in the custom face channel, and set it accordingly as well
        if( !matches_vertex_channel_custom_faces( cachedVertexChannel, vertChannel ) ) {
            m_vertChannelFaceDataChanged[channelName].clear();
            copy_vertex_channel_custom_faces( cachedVertexChannel, vertChannel );
        }
    }
}

// Same as retrieve_and_compare_mesh, but does it for one map channel of a mesh against the saved mesh
void xmesh_sequence_saver::retrieve_and_compare_face_channel( const trimesh3& mesh,
                                                              const frantic::tstring& channelName ) {

    // Grab the vertex channel data
    const_trimesh3_face_channel_general_accessor faceChannel = mesh.get_face_channel_general_accessor( channelName );

    face_channels_cache_t::iterator found = m_cachedFaceChannels.find( channelName );

    // first check if the saved mesh even has that channel
    if( found == m_cachedFaceChannels.end() ) {
        // channel will have to be added
        m_faceChannelDataChanged[channelName].clear();
        xmesh_cache_face_channel& cachedFaceChannel = m_cachedFaceChannels[channelName];
        copy_face_channel( cachedFaceChannel, faceChannel );
    } else {
        xmesh_cache_face_channel& cachedFaceChannel = found->second;

        if( !matches_face_channel( cachedFaceChannel, faceChannel ) ) {
            m_faceChannelDataChanged[channelName].clear();
            copy_face_channel( cachedFaceChannel, faceChannel );
        }
    }
}

void xmesh_sequence_saver::retrieve_and_compare_mesh( const frantic::geometry::polymesh3_ptr polymesh ) {

    bool invalidateFaces = false;
    bool isIncomingTriangleMesh = polymesh->is_triangle_mesh();

    if( !isIncomingTriangleMesh ) {
        if( m_isCachedTriangleMesh ) {
            invalidateFaces = true;
            polymesh->get_face_end_offsets( m_cachedFaceOffsets );
        } else {
            std::vector<int> newFaceOffsets;
            polymesh->get_face_end_offsets( newFaceOffsets );

            if( newFaceOffsets != m_cachedFaceOffsets ) {
                m_cachedFaceOffsets.swap( newFaceOffsets );
                invalidateFaces = true;
            }
        }
    } else {
        if( !m_isCachedTriangleMesh ) {
            m_cachedFaceOffsets.clear();
            invalidateFaces = true;
        }
    }

    m_isCachedTriangleMesh = isIncomingTriangleMesh;

    /////
    // Check the mesh vertices
    /////

    polymesh3_const_vertex_accessor<void> geometryAccessor = polymesh->get_const_vertex_accessor( _T( "verts" ) );

    if( !matches_vertex_channel( m_cachedGeometryChannel, geometryAccessor ) ) {
        m_savedVertsFile.clear();
        copy_vertex_channel( m_cachedGeometryChannel, geometryAccessor );
    }

    /////
    // Check the mesh faces
    /////
    if( invalidateFaces || !matches_vertex_channel_custom_faces( m_cachedGeometryChannel, geometryAccessor ) ) {
        m_savedFacesFile.clear();
        copy_vertex_channel_custom_faces( m_cachedGeometryChannel, geometryAccessor );
    }

    /////
    // Check all the vertex channels
    /////

    std::vector<frantic::tstring> newVertChannelNames;
    polymesh->get_vertex_channel_names( newVertChannelNames );

    // only do the extra channel stuff if there are actually extra channels in the new
    // mesh, otherwise just wipe out the old channel data
    if( !newVertChannelNames.empty() ) {
        for( unsigned channel = 0; channel < newVertChannelNames.size(); ++channel ) {
            retrieve_and_compare_vertex_channel( polymesh, newVertChannelNames[channel], invalidateFaces );
        }

        // if we have extra channels in the saved mesh that aren't in the new one, get rid of them
        std::vector<frantic::tstring> currentChannels;
        for( vertex_channels_cache_t::iterator it = m_cachedVertexChannels.begin(); it != m_cachedVertexChannels.end();
             ++it ) {
            currentChannels.push_back( it->first );
        }

        for( size_t i = 0; i < currentChannels.size(); ++i ) {
            if( !polymesh->has_vertex_channel( currentChannels[i] ) ) {
                m_cachedVertexChannels.erase( currentChannels[i] );
                m_vertChannelVertDataChanged.erase( currentChannels[i] );
                m_vertChannelFaceDataChanged.erase( currentChannels[i] );
            }
        }
    } else {
        m_cachedVertexChannels.clear();
        m_vertChannelVertDataChanged.clear();
        m_vertChannelFaceDataChanged.clear();
    }

    /////
    // Check all the face channels
    /////

    std::vector<frantic::tstring> newFaceChannelNames;
    polymesh->get_face_channel_names( newFaceChannelNames );

    // only do the extra channel stuff if there are actually extra channels in the new
    // mesh, otherwise just wipe out the old channel data
    if( !newFaceChannelNames.empty() ) {
        for( unsigned channel = 0; channel < newFaceChannelNames.size(); ++channel ) {
            retrieve_and_compare_face_channel( polymesh, newFaceChannelNames[channel], invalidateFaces );
        }

        // if we have extra channels in the saved mesh that aren't in the new one, get rid of them
        std::vector<frantic::tstring> currentChannels;
        for( face_channels_cache_t::iterator it = m_cachedFaceChannels.begin(); it != m_cachedFaceChannels.end();
             ++it ) {
            currentChannels.push_back( it->first );
        }

        for( size_t i = 0; i < currentChannels.size(); ++i ) {
            if( !polymesh->has_face_channel( currentChannels[i] ) ) {
                m_cachedFaceChannels.erase( currentChannels[i] );
                m_faceChannelDataChanged.erase( currentChannels[i] );
            }
        }
    } else {
        m_cachedFaceChannels.clear();
        m_faceChannelDataChanged.clear();
    }
}

void xmesh_sequence_saver::retrieve_and_compare_vertex_channel( const frantic::geometry::polymesh3_ptr polymesh,
                                                                const frantic::tstring& channelName,
                                                                bool invalidateFaces ) {
    // Grab the vertex channel data
    polymesh3_const_vertex_accessor<void> vertChannel = polymesh->get_const_vertex_accessor( channelName );
    vertex_channels_cache_t::iterator found = m_cachedVertexChannels.find( channelName );

    // first check if the saved mesh even has that channel
    if( invalidateFaces || found == m_cachedVertexChannels.end() ) {
        // channel will have to be added
        m_vertChannelVertDataChanged[channelName].clear();
        m_vertChannelFaceDataChanged[channelName].clear();
        xmesh_cache_vertex_channel& cachedVertexChannel = m_cachedVertexChannels[channelName];
        copy_vertex_channel( cachedVertexChannel, vertChannel );
        copy_vertex_channel_custom_faces( cachedVertexChannel, vertChannel );
    } else {
        // the channel exists, so grab the corresponding saved vertex channel data too
        xmesh_cache_vertex_channel& cachedVertexChannel = found->second;

        // Check for matching vertex data within the channel, re-write if no match
        if( !matches_vertex_channel( cachedVertexChannel, vertChannel ) ) {
            m_vertChannelVertDataChanged[channelName].clear();
            copy_vertex_channel( cachedVertexChannel, vertChannel );
        }

        // Check for changes in the custom face channel, and set it accordingly as well
        if( !matches_vertex_channel_custom_faces( cachedVertexChannel, vertChannel ) ) {
            m_vertChannelFaceDataChanged[channelName].clear();
            copy_vertex_channel_custom_faces( cachedVertexChannel, vertChannel );
        }
    }
}

void xmesh_sequence_saver::retrieve_and_compare_face_channel( const frantic::geometry::polymesh3_ptr polymesh,
                                                              const frantic::tstring& channelName,
                                                              bool invalidateFaces ) {
    // Grab the vertex channel data
    polymesh3_const_face_accessor<void> faceChannel = polymesh->get_const_face_accessor( channelName );

    face_channels_cache_t::iterator found = m_cachedFaceChannels.find( channelName );

    // first check if the saved mesh even has that channel
    if( invalidateFaces || found == m_cachedFaceChannels.end() ) {
        // channel will have to be added
        m_faceChannelDataChanged[channelName].clear();
        xmesh_cache_face_channel& cachedFaceChannel = m_cachedFaceChannels[channelName];
        copy_face_channel( cachedFaceChannel, faceChannel );
    } else {
        xmesh_cache_face_channel& cachedFaceChannel = found->second;

        if( !matches_face_channel( cachedFaceChannel, faceChannel ) ) {
            m_faceChannelDataChanged[channelName].clear();
            copy_face_channel( cachedFaceChannel, faceChannel );
        }
    }
}

// For all the places where there's a flag set (ie, an empty filename), this saves the data
// to a filename derived from the provided one.
void xmesh_sequence_saver::write_new_channels( const boost::filesystem::path& filename ) {
    error_queue_t errorQueue;
    { // scope for threadPool
        thread_pool threadPool( std::max<std::size_t>( 1, m_threadCount ) );

        if( m_savedVertsFile.empty() ) {
            boost::filesystem::path vertsFile( files::replace_extension(
                files::add_before_sequence_number( filename.wstring(), L"_verts" ), L".xmdat" ) );
            write_xmesh_data_file( threadPool, m_cachedGeometryChannel, vertsFile, m_compressionLevel,
                                   boost::ref( errorQueue ), boost::ref( m_savedVertsFile ) );
        }

        if( m_savedFacesFile.empty() ) {
            boost::filesystem::path facesFile( files::replace_extension(
                files::add_before_sequence_number( filename.wstring(), L"_faces" ), L".xmdat" ) );
            write_xmesh_vertex_custom_faces_data_file( threadPool, m_cachedGeometryChannel, facesFile,
                                                       m_compressionLevel, m_cachedFaceOffsets,
                                                       boost::ref( errorQueue ), boost::ref( m_savedFacesFile ) );
        }

        //  Loop through the channel flags and write out any data that changed

        // vertex channel vert data
        std::map<frantic::tstring, boost::filesystem::path::string_type>::iterator i;
        for( i = m_vertChannelVertDataChanged.begin(); i != m_vertChannelVertDataChanged.end(); ++i ) {
            if( i->second.empty() ) {
                vertex_channels_cache_t::iterator channelIt = m_cachedVertexChannels.find( i->first );
                xmesh_cache_vertex_channel& vertexChannel = channelIt->second;
                boost::filesystem::path vertsFile( files::replace_extension(
                    files::add_before_sequence_number(
                        filename.wstring(), L"_channel_" + wstring_from_channel_name( i->first ) + L"_verts" ),
                    L".xmdat" ) );
                write_xmesh_data_file( threadPool, vertexChannel, vertsFile, m_compressionLevel,
                                       boost::ref( errorQueue ), boost::ref( i->second ) );
            }
        }
        // vertex channel custom face data
        for( i = m_vertChannelFaceDataChanged.begin(); i != m_vertChannelFaceDataChanged.end(); ++i ) {
            if( i->second.empty() ) {
                vertex_channels_cache_t::iterator channelIt = m_cachedVertexChannels.find( i->first );
                xmesh_cache_vertex_channel& vertexChannel = channelIt->second;

                if( vertexChannel.m_customFaces.size() > 0 ) {
                    boost::filesystem::path facesFile( files::replace_extension(
                        files::add_before_sequence_number(
                            filename.wstring(), L"_channel_" + wstring_from_channel_name( i->first ) + L"_faces" ),
                        L".xmdat" ) );
                    write_xmesh_vertex_custom_faces_data_file( threadPool, vertexChannel, facesFile, m_compressionLevel,
                                                               m_cachedFaceOffsets, boost::ref( errorQueue ),
                                                               boost::ref( i->second ) );
                }
            }
        }

        // face channel data
        for( i = m_faceChannelDataChanged.begin(); i != m_faceChannelDataChanged.end(); ++i ) {
            if( i->second.empty() ) {
                face_channels_cache_t::iterator channelIt = m_cachedFaceChannels.find( i->first );
                xmesh_cache_face_channel& facesChannel = channelIt->second;
                boost::filesystem::path facesFile( files::replace_extension(
                    files::add_before_sequence_number(
                        filename.wstring(), L"_channel_" + wstring_from_channel_name( i->first ) + L"_facedata" ),
                    L".xmdat" ) );
                write_xmesh_data_file( threadPool, facesChannel, facesFile, m_compressionLevel,
                                       boost::ref( errorQueue ), boost::ref( i->second ) );
            }
        }
    }

    // All threads are joined now,
    // because the thread pool is destroyed (out of scope)

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

    // Make sure an output file was assigned for every channel.
    // This check should be redundant: any errors that prevent the assignment
    // should be reported in the errorQueue.
    if( m_savedVertsFile.empty() ) {
        throw std::runtime_error( "Internal Error: missing output data file for vertices." );
    }
    if( m_savedFacesFile.empty() ) {
        throw std::runtime_error( "Internal Error: missing output data file for faces." );
    }

    std::map<frantic::tstring, boost::filesystem::path::string_type>::iterator i;
    for( i = m_vertChannelVertDataChanged.begin(); i != m_vertChannelVertDataChanged.end(); ++i ) {
        if( i->second.empty() ) {
            throw std::runtime_error( "Internal Error: missing output data file for vertex channel \"" +
                                      frantic::strings::to_string( i->first ) + "\"." );
        }
    }
    for( i = m_vertChannelFaceDataChanged.begin(); i != m_vertChannelFaceDataChanged.end(); ++i ) {
        if( i->second.empty() ) {
            vertex_channels_cache_t::iterator channelIt = m_cachedVertexChannels.find( i->first );
            xmesh_cache_vertex_channel& vertexChannel = channelIt->second;

            if( vertexChannel.m_customFaces.size() > 0 ) {
                throw std::runtime_error( "Internal Error: missing output face file for vertex channel \"" +
                                          frantic::strings::to_string( i->first ) + "\"." );
            }
        }
    }
    for( i = m_faceChannelDataChanged.begin(); i != m_faceChannelDataChanged.end(); ++i ) {
        if( i->second.empty() ) {
            throw std::runtime_error( "Internal Error: missing output data file for face channel \"" +
                                      frantic::strings::to_string( i->first ) + "\"." );
        }
    }
}

// Writes out all the tags and corresponding information needed to locate and read the
// data files.
void xmesh_sequence_saver::write_xmesh_xml_file( const boost::filesystem::path& filename,
                                                 const frantic::geometry::xmesh_metadata* metadata ) {
    tinyxml2::XMLDocument doc( filename.c_str() );
    doc.LinkEndChild( doc.NewDeclaration() );

    tinyxml2::XMLElement* xmesh = doc.NewElement( "xmesh" );
    doc.LinkEndChild( xmesh );

    tinyxml2::XMLElement* pVersion = xmesh->LinkEndChild( doc.NewElement( "version" ) )->ToElement();
    pVersion->LinkEndChild( doc.NewText( "1" ) );

    tinyxml2::XMLElement* vertCount = doc.NewElement( "vertCount" );
    xmesh->LinkEndChild( vertCount );
    vertCount->LinkEndChild(
        doc.NewText( frantic::strings::to_string_classic<std::string>( num_elements( m_cachedGeometryChannel ) ).c_str() ) );
    tinyxml2::XMLElement* verts = doc.NewElement( "verts" );
    xmesh->LinkEndChild( verts );
    verts->LinkEndChild(
        doc.NewText( frantic::strings::to_utf8( boost::filesystem::path( m_savedVertsFile ).wstring() ).c_str() ) );

    tinyxml2::XMLElement* faceCount = doc.NewElement( "faceCount" );
    xmesh->LinkEndChild( faceCount );
    // todo: switch based on whether its a trimesh or not
    faceCount->LinkEndChild( doc.NewText(
        frantic::strings::to_string_classic<std::string>( m_cachedGeometryChannel.m_customFaces.size() /*/ 3*/ ).c_str() ) );
    tinyxml2::XMLElement* faces = doc.NewElement( "faces" );
    xmesh->LinkEndChild( faces );
    faces->LinkEndChild(
        doc.NewText( frantic::strings::to_utf8( boost::filesystem::path( m_savedFacesFile ).wstring() ).c_str() ) );

    frantic::graphics::boundbox3f boundbox( get_bound_box( m_cachedGeometryChannel ) );
    tinyxml2::XMLElement* boundboxElement = doc.NewElement( "boundbox" );
    xmesh->LinkEndChild( boundboxElement );
    tinyxml2::XMLElement* boundboxMin = doc.NewElement( "minimum" );
    boundboxElement->LinkEndChild( boundboxMin );
    boundboxMin->LinkEndChild( doc.NewText( boundbox.minimum().str().c_str() ) );
    tinyxml2::XMLElement* boundboxMax = doc.NewElement( "maximum" );
    boundboxElement->LinkEndChild( boundboxMax );
    boundboxMax->LinkEndChild( doc.NewText( boundbox.maximum().str().c_str() ) );

    std::map<frantic::tstring, boost::filesystem::path::string_type>::const_iterator i;
    for( i = m_vertChannelVertDataChanged.begin(); i != m_vertChannelVertDataChanged.end(); ++i ) {
        tinyxml2::XMLElement* map = doc.NewElement( "vertexChannel" );
        xmesh->LinkEndChild( map );

        tinyxml2::XMLElement* dataType = doc.NewElement( "type" );
        map->LinkEndChild( dataType );
        vertex_channels_cache_t::iterator channelIt = m_cachedVertexChannels.find( i->first );
        xmesh_cache_vertex_channel& vertexChannel = channelIt->second;

        // trimesh3_vertex_channel_general_accessor vertChannel =
        // m_savedMeshData.get_vertex_channel_general_accessor(i->first);
        dataType->LinkEndChild( doc.NewText(
            frantic::strings::to_string( channel_data_type_str( vertexChannel.m_arity, vertexChannel.m_dataType ) ).c_str() ) );

        const std::string channelNameUTF8( utf8_from_channel_name( i->first ) );
        if( !frantic::strings::is_valid_utf8( channelNameUTF8 ) ) {
            throw std::runtime_error( "xmesh_sequence_saver::write_xmesh_xml_file Internal Error: While writing \"" +
                                      filename.string() +
                                      "\": Vertex channel name was not converted to valid UTF-8: \"" +
                                      frantic::strings::to_string( i->first ) + "\"." );
        }
        tinyxml2::XMLElement* channelName = doc.NewElement( "name" );
        map->LinkEndChild( channelName );
        channelName->LinkEndChild( doc.NewText( channelNameUTF8.c_str() ) );

        tinyxml2::XMLElement* vertCount = doc.NewElement( "vertCount" );
        map->LinkEndChild( vertCount );
        vertCount->LinkEndChild(
            doc.NewText( frantic::strings::to_string_classic<std::string>( num_elements( vertexChannel ) ).c_str() ) );

        const boost::filesystem::path vertDataPath( boost::filesystem::path( i->second ) );
        const std::wstring vertData( vertDataPath.wstring() );
        const std::string vertDataUTF8( frantic::strings::to_utf8( vertData ) );
        if( !frantic::strings::is_valid_utf8( vertDataUTF8 ) ) {
            throw std::runtime_error(
                "xmesh_sequence_saver::write_xmesh_xml_file Internal Error: While writing \"" + filename.string() +
                "\": Vertex data filename was not converted to valid UTF-8: \"" + vertDataPath.string() + "\"." );
        }
        tinyxml2::XMLElement* mapVerts = doc.NewElement( "vertData" );
        map->LinkEndChild( mapVerts );
        mapVerts->LinkEndChild( doc.NewText( vertDataUTF8.c_str() ) );

        /* Check the custom faces flag before adding the tag for faces */
        if( vertexChannel.m_customFaces.size() > 0 ) {
            tinyxml2::XMLElement* faceCount = doc.NewElement( "faceCount" );
            map->LinkEndChild( faceCount );

            // todo: switch based on whether its a trimesh or not
            faceCount->LinkEndChild( doc.NewText(
                frantic::strings::to_string_classic<std::string>( vertexChannel.m_customFaces.size() /*/ 3*/ ).c_str() ) );

            const boost::filesystem::path faceDataPath( m_vertChannelFaceDataChanged[i->first] );
            const std::wstring faceData( faceDataPath.wstring() );
            const std::string faceDataUTF8( frantic::strings::to_utf8( faceData ) );
            if( !frantic::strings::is_valid_utf8( faceDataUTF8 ) ) {
                throw std::runtime_error(
                    "xmesh_sequence_saver::write_xmesh_xml_file Internal Error: While writing \"" + filename.string() +
                    "\": Vertex face filename was not converted to valid UTF-8: \"" + faceDataPath.string() + "\"." );
            }
            tinyxml2::XMLElement* mapFaces = doc.NewElement( "faceData" );
            map->LinkEndChild( mapFaces );
            mapFaces->LinkEndChild( doc.NewText( faceDataUTF8.c_str() ) );
        }
    }

    for( i = m_faceChannelDataChanged.begin(); i != m_faceChannelDataChanged.end(); ++i ) {
        tinyxml2::XMLElement* map = doc.NewElement( "faceChannel" );
        xmesh->LinkEndChild( map );

        tinyxml2::XMLElement* dataType = doc.NewElement( "type" );
        map->LinkEndChild( dataType );
        face_channels_cache_t::iterator channelIt = m_cachedFaceChannels.find( i->first );
        xmesh_cache_face_channel& facesChannel = channelIt->second;

        // trimesh3_face_channel_general_accessor faceChannel = m_savedMeshData.get_face_channel_general_accessor;
        dataType->LinkEndChild( doc.NewText(
            frantic::strings::to_string( channel_data_type_str( facesChannel.m_arity, facesChannel.m_dataType ) ).c_str() ) );

        const std::string channelNameUTF8( utf8_from_channel_name( i->first ) );
        if( !frantic::strings::is_valid_utf8( channelNameUTF8 ) ) {
            throw std::runtime_error( "xmesh_sequence_saver::write_xmesh_xml_file Internal Error: While writing \"" +
                                      filename.string() + "\": Face channel name was not converted to valid UTF-8: \"" +
                                      frantic::strings::to_string( i->first ) + "\"." );
        }
        tinyxml2::XMLElement* channelName = doc.NewElement( "name" );
        map->LinkEndChild( channelName );
        channelName->LinkEndChild( doc.NewText( channelNameUTF8.c_str() ) );

        tinyxml2::XMLElement* faceCount = doc.NewElement( "faceCount" );
        map->LinkEndChild( faceCount );
        faceCount->LinkEndChild(
            doc.NewText( frantic::strings::to_string_classic<std::string>( num_elements( facesChannel ) ).c_str() ) );

        const boost::filesystem::path faceDataPath( i->second );
        const std::wstring faceData( faceDataPath.wstring() );
        const std::string faceDataUTF8( frantic::strings::to_utf8( faceData ) );
        if( !frantic::strings::is_valid_utf8( faceDataUTF8 ) ) {
            throw std::runtime_error(
                "xmesh_sequence_saver::write_xmesh_xml_file Internal Error: While writing \"" + filename.string() +
                "\": Face data filename was not converted to valid UTF-8: \"" + faceDataPath.string() + "\"." );
        }
        tinyxml2::XMLElement* mapFaces = doc.NewElement( "faceData" );
        map->LinkEndChild( mapFaces );
        mapFaces->LinkEndChild( doc.NewText( faceDataUTF8.c_str() ) );
    }

    if( metadata ) {
        write_xmesh_metadata( doc, *metadata );
    }

    frantic::files::file_ptr fout( frantic::files::fopen( filename, "w" ) );
    tinyxml2::XMLError result = doc.SaveFile( fout );
    if( result != tinyxml2::XMLError::XML_SUCCESS )
        throw std::runtime_error( "xmesh_sequence_saver::write_xmesh_xml_file: Failed to save to file \"" +
                                  filename.string() + "\"." );
}

xmesh_sequence_saver::xmesh_sequence_saver()
    : m_isCachedTriangleMesh( false )
    , m_compressionLevel( Z_DEFAULT_COMPRESSION )
    , m_threadCount( 1 ) {}

xmesh_sequence_saver::~xmesh_sequence_saver() {}

/**
 * Clears the currently saved mesh data.
 * <p>
 * This function clears the cached trimesh object and resets all the
 * stored file "pointers".
 *
 */
void xmesh_sequence_saver::clear() {
    m_cachedGeometryChannel.m_data.clear();
    m_cachedVertexChannels.clear();
    m_cachedFaceChannels.clear();
    m_savedVertsFile.clear();
    m_savedFacesFile.clear();
    m_vertChannelVertDataChanged.clear();
    m_vertChannelFaceDataChanged.clear();
    m_faceChannelDataChanged.clear();
}

/**
 * Writes a given trimesh3 to the given filename
 * using the xmesh format.
 * <p>
 * This function checks the given mesh against the last mesh written with
 * the xmesh_sequence_saver instantiation and saves out any new data.
 * @param  mesh		the trimesh3 to be saved out
 * @param  metadata the (optional) metadata to include in the xmesh file
 * @param  filename	the location/name of the destination xmesh file
 */
void xmesh_sequence_saver::write_xmesh( const trimesh3& mesh, const xmesh_metadata* metadata,
                                        const frantic::tstring& filename ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    boost::filesystem::path filepath( frantic::strings::to_wstring( filename ) );
#else
    boost::filesystem::path filepath( filename );
#endif

    // Load the mesh data into this class's internal structures.  For any channel
    // which wasn't changed by the update, the filename is left as is, but for ones which
    // did change the filename is reset to the empty string.
    retrieve_and_compare_mesh( mesh );

    // Save all the channel data that wasn't kept the same from previous saves.  This function
    // uses the provided filename as its template for creating the files, and inserts just the filename
    // into all the cached filename variables.
    write_new_channels( filepath );

    // Save the xml file with the references to all this data
    write_xmesh_xml_file( filepath, metadata );
}

void xmesh_sequence_saver::write_xmesh( const frantic::geometry::polymesh3_ptr polymesh, const xmesh_metadata* metadata,
                                        const frantic::tstring& filename ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    boost::filesystem::path filepath( frantic::strings::to_wstring( filename ) );
#else
    boost::filesystem::path filepath( filename );
#endif

    // Load the mesh data into this class's internal structures.  For any channel
    // which wasn't changed by the update, the filename is left as is, but for ones which
    // did change the filename is reset to the empty string.
    retrieve_and_compare_mesh( polymesh );

    // Save all the channel data that wasn't kept the same from previous saves.  This function
    // uses the provided filename as its template for creating the files, and inserts just the filename
    // into all the cached filename variables.
    write_new_channels( filepath );

    // Save the xml file with the references to all this data
    write_xmesh_xml_file( filepath, metadata );
}

/**
 * Writes a given trimesh3 to the given filename
 * using the xmesh format.
 * <p>
 * This function checks the given mesh against the last mesh written with
 * the xmesh_sequence_saver instantiation and saves out any new data.
 * @param  mesh		the trimesh3 to be saved out
 * @param  filename	the location/name of the destination xmesh file
 */
void xmesh_sequence_saver::write_xmesh( const trimesh3& mesh, const frantic::tstring& filename ) {
    write_xmesh( mesh, (xmesh_metadata*)NULL, filename );
}

/**
 * Writes a given trimesh3 to the given filename
 * using the xmesh format.
 * <p>
 * This function checks the given mesh against the last mesh written with
 * the xmesh_sequence_saver instantiation and saves out any new data.
 * @param  mesh		the trimesh3 to be saved out
 * @param  metadata the metadata to include in the xmesh file
 * @param  filename	the location/name of the destination xmesh file
 */
void xmesh_sequence_saver::write_xmesh( const trimesh3& mesh, const xmesh_metadata& metadata,
                                        const frantic::tstring& filename ) {
    write_xmesh( mesh, &metadata, filename );
}

void xmesh_sequence_saver::write_xmesh( const frantic::geometry::polymesh3_ptr polymesh,
                                        const frantic::tstring& filename ) {
    write_xmesh( polymesh, (xmesh_metadata*)NULL, filename );
}

void xmesh_sequence_saver::write_xmesh( const frantic::geometry::polymesh3_ptr polymesh,
                                        const frantic::geometry::xmesh_metadata& metadata,
                                        const frantic::tstring& filename ) {
    write_xmesh( polymesh, &metadata, filename );
}

void xmesh_sequence_saver::set_compression_level( int compressionLevel ) { m_compressionLevel = compressionLevel; }

void xmesh_sequence_saver::set_thread_count( std::size_t threadCount ) { m_threadCount = threadCount; }

} // namespace geometry
} // namespace frantic
