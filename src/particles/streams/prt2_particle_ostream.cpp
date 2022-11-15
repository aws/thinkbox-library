// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/prt2_particle_ostream.hpp>

#pragma warning( push )
#pragma warning( disable : 4100 4512 )
#include <tbb/pipeline.h>
#pragma warning( pop )

#include <tbb/tbb_thread.h>

using namespace std;
using namespace frantic;
using namespace frantic::prtfile;
using namespace frantic::particles::streams;
using frantic::prtfile::prt2_writer;

namespace {

bool try_pop( tbb::concurrent_queue<vector<char>>& q, std::vector<char>& buffer ) {
#if TBB_INTERFACE_VERSION_MAJOR >= 4
    return q.try_pop( buffer );
#else
    return q.pop_if_present( buffer );
#endif
}

} // anonymous namespace

namespace detail {
/**
 * A particle chunk generator which reads from a concurrent queue and creates particle chunks from it. This is
 * needed to translate the push interface of the ostream into the pull interface of the prt2_writer.
 */
class queued_chunk_generator : public tbb::filter, boost::noncopyable {
    typedef prt2_writer::particle_chunk chunk_type;

    std::size_t m_structureSize;
    tbb::concurrent_queue<vector<char>>& m_particleChunkQueue;
    volatile bool& m_closeRequested;
    bool m_sentTerminationChunk;

  public:
    queued_chunk_generator( std::size_t structureSize, tbb::concurrent_queue<vector<char>>& particleChunkQueue,
                            volatile bool& closeRequested )
        : tbb::filter( true )
        , m_structureSize( structureSize )
        , m_particleChunkQueue( particleChunkQueue )
        , m_closeRequested( closeRequested )
        , m_sentTerminationChunk( false ) {}

    /**
     * Pop the particle chunk queue and return the result. Waits for a chunk in the queue if none is available.
     */
    void* operator()( void* ) {
        std::vector<char> chunkBuffer;

        // If a chunk can be retrieved, pass it along.
        if( try_pop( m_particleChunkQueue, chunkBuffer ) ) {
            chunk_type* chunk = new chunk_type();
            chunk->particleCount = chunkBuffer.size() / m_structureSize;
            chunk->uncompressed = chunkBuffer;
            return chunk;
        }

        // If no close has been requested, send an ignore chunk. Otherwise, return NULL to terminate the generator.
        if( !m_closeRequested ) {
            return &prt2_writer::IGNORE_CHUNK;
        }

        // If we have not yet sent a termination chunk, send one before terminating.
        if( !m_sentTerminationChunk ) {
            m_sentTerminationChunk = true;
            return &prt2_writer::TERMINATION_CHUNK;
        }
        return NULL;
    }
};

/**
 * A task that is spawned so that the pipelined prt2 writing can happen in parallel with the ostream.
 */
class queued_write_task : public tbb::task {
    boost::shared_ptr<prt2_particle_ostream::modal_pipeline> m_chunkPipeline;
    boost::shared_ptr<prt2_writer::chunk_writer> m_chunkWriter;
    prt2_writer& m_prt2;

  public:
    queued_write_task( boost::shared_ptr<prt2_particle_ostream::modal_pipeline> chunkPipeline,
                       boost::shared_ptr<prt2_writer::chunk_writer> chunkWriter, frantic::prtfile::prt2_writer& prt2 )
        : m_chunkPipeline( chunkPipeline )
        , m_chunkWriter( chunkWriter )
        , m_prt2( prt2 ) {}

    tbb::task* execute() {
        m_prt2.write_particle_chunks( m_chunkPipeline, m_chunkWriter );
        return NULL;
    }
};
} // namespace detail

/**
 * prt2_particle_ostream::modal_pipeline
 */
void prt2_particle_ostream::modal_pipeline::add_filter( tbb::filter& filter ) {
    m_filters.push_back( &filter );
    m_autoPipeline.add_filter( filter );
}

void prt2_particle_ostream::modal_pipeline::clear() {
    m_data = this;
    m_filters.clear();
    m_autoPipeline.clear();
}

void prt2_particle_ostream::modal_pipeline::run( size_t maxLiveTokens ) {
    // If we have finished or are already running, just return.
    if( m_data == NULL || is_running() ) {
        return;
    }

    if( m_filters.empty() ) {
        throw std::runtime_error( "modal_pipeline::run() should not be called if empty!" );
    }

    // Atomically switch the state to automatic running from stopped. If we fail to switch, just return.
    State oldState = static_cast<State>( m_state.compare_and_swap( AutomaticRunning, Stopped ) );
    if( oldState != Stopped ) {
        return;
    }

    m_autoPipeline.run( maxLiveTokens );
    m_data = NULL;
    m_state = Complete;
}

bool prt2_particle_ostream::modal_pipeline::step() {
    if( m_state == AutomaticRunning ) {
        throw std::runtime_error( "modal_pipeline::step() should not be called when automatic running!" );
    }

    // Stepping after completion just does nothing.
    if( !m_data ) {
        return false;
    }

    // Stepping before starting will block.
    while( m_state == Stopped ) {
        tbb::this_tbb_thread::yield();
    }

    for( std::vector<tbb::filter*>::iterator i = m_filters.begin(), end = m_filters.end(); i != end; ++i ) {
        tbb::filter& filter = **i;
        m_data = filter( m_data );

        if( m_data == NULL ) {
            m_state = Complete;
            return false;
        }
    }
    return true;
}

void prt2_particle_ostream::modal_pipeline::manual_run_non_blocking() {
    // If we have finished or are already running, just return.
    if( m_data == NULL || is_running() ) {
        return;
    }

    if( m_filters.empty() ) {
        throw std::runtime_error( "modal_pipeline::manual_run_non_blocking() should not be called if empty!" );
    }

    m_state.compare_and_swap( ManualRunning, Stopped );
}

/**
 * prt2_particle_ostream
 */
prt2_particle_ostream::prt2_particle_ostream(
    const frantic::tstring& file, const frantic::channels::channel_map& particleChannelMap,
    const frantic::channels::channel_map& particleChannelMapForFile,
    frantic::prtfile::prt2_compression_t compressionScheme, bool useTempFile, const boost::filesystem::path& tempDir,
    const frantic::channels::property_map* globalMetadata,
    const std::map<frantic::tstring, frantic::channels::property_map>* channelMetadata,
    intptr_t desiredChunkSizeInBytes )
    : m_closeRequested( false )
    , m_chunkPipeline( new modal_pipeline() ) {
    // Open the output file, this writes the header and the initial 'Chan' chunk,
    // and initializes the file-layout channel_map inside of m_prt2
    m_prt2.open( file, particleChannelMapForFile, useTempFile, tempDir );

    // Write all the metadata
    if( globalMetadata != NULL ) {
        m_prt2.write_general_metadata_filechunks( *globalMetadata );
    }
    if( channelMetadata != NULL ) {
        for( std::map<frantic::tstring, frantic::channels::property_map>::const_iterator i = channelMetadata->begin(),
                                                                                         iend = channelMetadata->end();
             i != iend; ++i ) {
            m_prt2.write_channel_metadata_filechunks( i->first, i->second );
        }
    }

    set_channel_map( particleChannelMap );

    // The particle chunk buffer uses the channel map from the prt2_writer
    m_particleChunkBuffer.reserve( desiredChunkSizeInBytes );
    m_desiredChunkSizeInBytes = desiredChunkSizeInBytes;

    // Initialize the pipeline with its filters immediately, as we can't rely on the asynchronous task to do it if we
    // want to be able to switch to manual mode on the pipeline.
    boost::shared_ptr<::detail::queued_chunk_generator> chunkGenerator( new ::detail::queued_chunk_generator(
        m_prt2.get_channel_map().structure_size(), m_particleChunkQueue, m_closeRequested ) );
    m_chunkPipelineFilters.push_back( chunkGenerator );

    boost::shared_ptr<frantic::prtfile::prt2_writer::chunk_writer> chunkWriter;

    m_prt2.get_filters( 0, m_nullProgress, _T(""), false, compressionScheme, m_chunkPipelineFilters, chunkWriter );
    m_chunkPipelineFilters.push_back( chunkWriter );

    for( std::vector<boost::shared_ptr<tbb::filter>>::iterator i = m_chunkPipelineFilters.begin(),
                                                               end = m_chunkPipelineFilters.end();
         i != end; ++i ) {
        m_chunkPipeline->add_filter( **i );
    }

    // Initialize the dummy task which will hold the queued_write_task so it can run and we can wait for it.
    m_dummyTask = new( tbb::task::allocate_root() ) tbb::empty_task();
    m_dummyTask->set_ref_count( 2 );

    ::detail::queued_write_task& queuedWriteTask =
        *new( m_dummyTask->allocate_child() )::detail::queued_write_task( m_chunkPipeline, chunkWriter, m_prt2 );
    m_dummyTask->spawn( queuedWriteTask );

    // Sleep for a short period of time, and then make the modal_pipeline manual if it hasn't run.
    // 0.01 is a bit of a magic number. 0.001 is short enough that this consistently goes to manual mode, but 0.01 seems
    // to work well enough.
    tbb::this_tbb_thread::yield();
    tbb::this_tbb_thread::sleep( tbb::tick_count::interval_t( 0.01 ) );

    m_chunkPipeline->manual_run_non_blocking();
}

prt2_particle_ostream::~prt2_particle_ostream() { close(); }

void prt2_particle_ostream::close() {
    if( !m_closeRequested ) {
        // Finish the particle chunks
        if( !m_particleChunkBuffer.empty() ) {
            m_particleChunkQueue.push( m_particleChunkBuffer );
            m_particleChunkBuffer.clear();
        }
        m_closeRequested = true;

        // If in manual mode, step until we can't step any more.
        if( m_chunkPipeline->is_manually_running() ) {
            while( m_chunkPipeline->step() ) {
            }
        }

        // Wait until our child tasks are done.
        m_dummyTask->wait_for_all();
        m_dummyTask->destroy( *m_dummyTask );
        m_dummyTask = NULL;

        // Write the bounding box if we were accumulating it
        if( m_posAccessor.is_valid() ) {
            if( m_particleChannelMap[_T("Position")].data_type() == frantic::channels::data_type_float64 ) {
                m_prt2.write_channel_metadata_filechunk<frantic::graphics::boundbox3fd>( _T("Position"), _T("Extents"),
                                                                                         m_boundbox );
            } else {
                m_prt2.write_channel_metadata_filechunk<frantic::graphics::boundbox3f>(
                    _T("Position"), _T("Extents"), frantic::graphics::boundbox3f( m_boundbox ) );
            }
        }

        m_prt2.close();
    }
}

const boost::filesystem::path& prt2_particle_ostream::get_target_file() const { return m_prt2.get_target_file(); }

void prt2_particle_ostream::set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
    m_particleChannelMap = particleChannelMap;
    // Initialize the adaptor for converting the particle format to the one in the file
    m_pcmAdaptor.set( m_prt2.get_channel_map(), m_particleChannelMap );

    m_posAccessor.reset();
    if( m_particleChannelMap.has_channel( _T("Position") ) )
        m_posAccessor = m_particleChannelMap.get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
}

void prt2_particle_ostream::put_particle( const char* rawParticleData ) {
    if( m_closeRequested )
        throw std::runtime_error( "prt2_particle_ostream.put_particle: Tried to write to particle file \"" +
                                  frantic::strings::to_string( m_prt2.get_stream_name() ) + "\" after it was closed." );

    if( m_posAccessor.is_valid() )
        m_boundbox += m_posAccessor( rawParticleData );

    // Append the particle to the particle chunk buffer
    std::size_t nextI = m_particleChunkBuffer.size();
    m_particleChunkBuffer.resize( m_particleChunkBuffer.size() + m_pcmAdaptor.dest_size() );
    m_pcmAdaptor.copy_structure( &m_particleChunkBuffer[nextI], rawParticleData );

    // Flush the particle chunk buffer if it has reached the threshold
    boost::int64_t currentChunkSize =
        static_cast<boost::int64_t>( m_particleChunkBuffer.size() + m_prt2.get_channel_map().structure_size() );
    if( currentChunkSize > m_desiredChunkSizeInBytes ) {
        m_particleChunkQueue.push( m_particleChunkBuffer );
        m_particleChunkBuffer.clear();

        // If we are in manual mode, step the pipeline for every chunk we add.
        if( m_chunkPipeline->is_manually_running() ) {
            if( !m_chunkPipeline->step() ) {
                throw std::runtime_error(
                    "prt2_particle_ostream::put_particle() - step() should never return false here." );
            }
        }
    }
}
