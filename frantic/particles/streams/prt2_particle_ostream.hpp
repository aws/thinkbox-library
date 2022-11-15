// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/particle_ostream.hpp>
#include <frantic/prtfile/prt2_writer.hpp>

#include <tbb/task.h>
#include <tbb/task_scheduler_init.h>

#include <tbb/concurrent_queue.h>

namespace frantic {
namespace particles {
namespace streams {

class prt2_particle_ostream : public frantic::particles::streams::particle_ostream {
  public:
    /**
     * A pipeline that can be manually stepped through.
     */
// Allow us to compile with the deprecated tbb::pipeline for now...
#pragma warning( push )
#pragma warning( disable : 4996 )
    class modal_pipeline {
        enum State {
            Stopped = 0,
            AutomaticRunning = 1,
            ManualRunning = 2,
            Complete = 3,
        };

        tbb::atomic<boost::uint8_t> m_state;
        tbb::pipeline m_autoPipeline;        // The tbb pipeline that is used when this is in automatic mode.
        std::vector<tbb::filter*> m_filters; // The filters in the pipeline.
        void* m_data;                        // The data passed down the pipeline. NULL if the pipeline is done.

      public:
        modal_pipeline()
            : m_filters()
            , m_data( this ) {
            m_state = Stopped;
        }

        virtual ~modal_pipeline() {}

        bool is_manually_running() const { return m_state == ManualRunning; }
        bool is_running() const { return m_state == AutomaticRunning || m_state == ManualRunning; }

        void add_filter( tbb::filter& filter );
        void clear();

        /**
         * If the pipeline is not stopped, this does nothing. Otherwise, this runs the pipeline using TBB and blocks
         * until TBB is complete.
         */
        void run( size_t maxLiveTokens );

        /**
         * If the pipeline is stopped, this switches it to manually running but does not block until completion.
         */
        void manual_run_non_blocking();

        /**
         * If the pipeline is manual, this performs a single iteration of its filters. Otherwise, this throws an
         * exception.
         */
        bool step();
    };
#pragma warning( pop )

  private:
    frantic::prtfile::prt2_writer m_prt2;
    tbb::concurrent_queue<std::vector<char>> m_particleChunkQueue;

    // Holds the filters as the pipeline has a weak pointer.
    std::vector<boost::shared_ptr<tbb::filter>> m_chunkPipelineFilters;

    // The pipeline that does the heavy lifting of file writing.
    boost::shared_ptr<modal_pipeline> m_chunkPipeline;

    // Signalling variable to communicate with our child tasks.
    volatile bool m_closeRequested;

    // A dummy task which spawns our child task (the queued_write_task) and allows us to wait for its completion.
    tbb::task* m_dummyTask;

    // A null progress logger we can count on to be around as long as this instance is.
    frantic::logging::null_progress_logger m_nullProgress;

    frantic::channels::channel_map m_particleChannelMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;

    std::vector<char> m_particleChunkBuffer;
    boost::int64_t m_desiredChunkSizeInBytes;

    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> m_posAccessor;
    frantic::graphics::boundbox3fd m_boundbox;

    // Private copy constructor and assignment operator to disable copying
    prt2_particle_ostream( const prt2_particle_ostream& );            // not implemented
    prt2_particle_ostream& operator=( const prt2_particle_ostream& ); // not implemented

  public:
    prt2_particle_ostream(
        const frantic::tstring& file, const frantic::channels::channel_map& particleChannelMap,
        const frantic::channels::channel_map& particleChannelMapForFile,
        frantic::prtfile::prt2_compression_t compressionScheme = frantic::prtfile::prt2_compression_default,
        bool useTempFile = true, const boost::filesystem::path& tempDir = boost::filesystem::path(),
        const frantic::channels::property_map* generalMetadata = NULL,
        const std::map<frantic::tstring, frantic::channels::property_map>* channelMetadata = NULL,
        intptr_t desiredChunkSizeInBytes = 1000000 );

    virtual ~prt2_particle_ostream();

    void close();

    /** Get the file path where we are writing. */
    const boost::filesystem::path& get_target_file() const;

    /** Get the name of this ostream
     *  (this is not guaranteed to be the same as the output file.  Use get_target_file() instead for that) */
    const frantic::tstring& get_stream_name() const { return m_prt2.get_stream_name(); }

    /** This is the particle layout for particles provided to the stream. */
    const frantic::channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    /** Sets the particle layout for particles provided to the stream. */
    void set_channel_map( const frantic::channels::channel_map& particleChannelMap );

    /** This is how big one particle is */
    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    void put_particle( const char* rawParticleData );
};

} // namespace streams
} // namespace particles
} // namespace frantic
