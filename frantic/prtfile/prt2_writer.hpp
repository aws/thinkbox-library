// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <memory>
#include <set>

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/files/files.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/prtfile/prt2_common.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_exception.h>
#pragma warning( push )
#pragma warning( disable : 4100 4512 )
#include <tbb/pipeline.h>
#pragma warning( pop )

// Forward declaration
namespace tbb {
class filter;
}

namespace frantic {
namespace prtfile {

/**
 * An object to write to PRT2 files using a tbb::pipeline for multithreading.
 *
 * To use this class:
 * 1. Call open on the file or stream you want to write to.
 * 2. Call write_particle_chunks for each 'Part' or 'PrtO' file chunk you want to write, and call
 *    write_prt2_buffered_filechunk, write_metadata_filechunk, or write_channel_metadata_filechunk for each metadata
 *    file chunk to write. Metadata file chunks and particle file chunks can be written in any order.
 * 3. Call close.
 *
 */
class prt2_writer {
  public:
    struct particle_chunk_index_info {
        boost::int64_t PRTChunkSize, PRTChunkParticleCount;
    };

    /**
     * Writes a 'Part' or 'PrtO' file chunk to an output stream. This also handles deleting the chunks that came from
     * the generator.

     * This class is public as classes that use the prt2_writer and manually initialize their pipeline need to be aware
     * of it.
     */
    class chunk_writer : public tbb::filter, boost::noncopyable {
        prt2_writer& m_writer;

        std::ostream& m_outputStream;
        const frantic::tstring m_fileStreamName;
        const frantic::tstring m_particleStreamName;
        bool m_usePositionOffset;
        prt2_compression_t m_compressionScheme;

        // The positions within the file where the PRTs values started.
        boost::int64_t m_prtsChunkSizeSeek, m_prtsParticleCountSeek;

        boost::uint64_t m_particleCountTotal;
        logging::progress_logger& m_progress;
        bool m_isCancelled;
        bool m_hasWrittenHeader;

        // The number of particles and number of particle chunks written out so far, respectively.
        boost::uint64_t m_particleCount, m_particleChunkCount;

        // The particle chunk indexes written to this file chunk.
        std::vector<prt2_writer::particle_chunk_index_info> m_particleChunkIndex;

        // Write the 'Part' or 'PrtO' filechunk header.
        void begin_particle_chunks();

        // Fill in any missing information in the filechunk header and write the 'PInd' particle chunk index.
        void end_particle_chunks();

      public:
        chunk_writer( prt2_writer& writer, std::ostream& outputStream, const frantic::tstring& fileStreamName,
                      const frantic::tstring& particleStreamName, bool usePositionOffset,
                      prt2_compression_t compressionScheme, boost::uint64_t totalParticleCount,
                      logging::progress_logger& progress )
            : tbb::filter( true )
            , m_writer( writer )
            , m_outputStream( outputStream )
            , m_fileStreamName( fileStreamName )
            , m_particleStreamName( particleStreamName )
            , m_usePositionOffset( usePositionOffset )
            , m_compressionScheme( compressionScheme )
            , m_prtsChunkSizeSeek( -1 )
            , m_prtsParticleCountSeek( -1 )
            , m_particleCountTotal( totalParticleCount )
            , m_progress( progress )
            , m_isCancelled( false )
            , m_hasWrittenHeader( false )
            , m_particleCount( 0 )
            , m_particleChunkCount( 0 )
            , m_particleChunkIndex() {}

        void* operator()( void* item );

        bool is_cancelled() const { return m_isCancelled; }
    };

    /**
     * A particle chunk. A chunk generator must generate these so that they can be fed into the pipeline.
     */
    struct particle_chunk {
        particle_chunk()
            : particleCount( 0 )
            , uncompressed()
            , compressed()
            , positionOffset()
            , usePositionOffset( false ) {}

        // These fields are initialized by the chunk generator.
        boost::uint64_t particleCount;     // The number of particles in the chunk.
        std::vector<char> uncompressed;    // The uncompressed bytes of the particles in the chunk.
        graphics::vector3f positionOffset; // The position offset of the particles in the chunk.
        bool usePositionOffset;            // Whether or not to use the position index.

        // These fields are initialized by the compressor.
        std::vector<char> compressed; // The compressed bytes.
    };

    /**
     * A special particle chunk- the ignore chunk. If a chunk generator returns the address of this chunk, the pipeline
     * will ignore the chunk and the generator will continue to run.
     */
    static particle_chunk IGNORE_CHUNK;

    /**
     * A special particle chunk- the termination chunk. A chunk generator must return the address of this chunk once it
     * is done generating chunks so that the chunk_writer can write the end of the file chunk.
     */
    static particle_chunk TERMINATION_CHUNK;

    prt2_writer();
    ~prt2_writer();

    const channels::channel_map& get_channel_map() const { return m_fileParticleChannelMap; }

    /**
     * Gets the name of this stream.  Usually this is exactly the same as the path of the file being written to.
     * However, it can differ if it is being written to a temporary file or folder first (In this case, it will not be a
     * valid filename at all)
     */
    const frantic::tstring& get_stream_name() const { return m_streamname; }

    /**
     * Gets the name of the file that is ultimately being written to
     */
    const boost::filesystem::path& get_target_file() const { return m_targetFile; }

    // Writes a filechunk that was buffered in memory.
    void write_prt2_buffered_filechunk( const frantic::graphics::raw_byte_buffer& buf, boost::uint32_t chunkName );

    // Writes a 'Meta' filechunk for a general metadata property.
    template <class T>
    void write_general_metadata_filechunk( const frantic::tstring& propertyName, const T& value ) {
        write_metadata_filechunk( propertyName, value );
    }

    // Writes a 'Meta' filechunk for the channel metadata.
    template <class T>
    void write_channel_metadata_filechunk( const frantic::tstring& channelName, const frantic::tstring& propertyName,
                                           const T& value ) {
        write_metadata_filechunk( channelName + _T(".") + propertyName, value );
    }
    // Writes a 'Meta' filechunk for each property in the property_map.
    void write_general_metadata_filechunks( const frantic::channels::property_map& pm );

    // Writes a 'Meta' filechunk for each property in the property_map.
    void write_channel_metadata_filechunks( const frantic::tstring& channelName,
                                            const frantic::channels::property_map& pm );

    /**
     * Write the header and the 'Chan' filechunk.
     *
     * \note The provided particleChannelMap is NOT preserved as the layout of particles that
     *       must be provided when writing particle chunks. The prt2_writer creates a new channel map
     *       according to the PRT spec, and the caller must call prt2_writer::get_channel_map() to
     *       get the actual particle layout for writing. This is handled by the prt2_ostream, for
     *       example.
     */
    void open( const frantic::tstring& filename, const channels::channel_map& particleChannelMap,
               bool useTempFile = true, const boost::filesystem::path& tempDir = boost::filesystem::path() );
    void open( std::ostream* os, const frantic::tstring& filename, const channels::channel_map& particleChannelMap );

    void close();

    /**
     * Get the filters this PRT2 writer would add to its pipeline to produce the correct 'Part' or 'PrtO' file chunk.
     *
     * \param totalParticleCount The total number of particles to write. Only used for progress logging.
     * \param progress Where to log progress to.
     * \param compressionScheme The compression method for the particle chunks in this file chunk.
     * \param outParticleChunkIndex The particle chunk index that outChunkWriter will write to.
     * \param outFilters The intermidiary filters in the pipeline.
     * \param outChunkWriter The filter which writes to the file. This is provided explicitly as it is guaranteed to
     *                       exist and holds information about cancellation.
     */
    void get_filters( boost::uint64_t totalParticleCount, frantic::logging::progress_logger& progress,
                      const frantic::tstring& particleStreamName, bool usePositionOffset,
                      prt2_compression_t compressionScheme, std::vector<boost::shared_ptr<tbb::filter>>& outFilters,
                      boost::shared_ptr<chunk_writer>& outChunkWriter );

    /**
     * Run the chunk pipeline to write particle chunks to a file within a 'Part' or 'PrtO' file chunk.
     *
     * \param chunkPipeline The pipeline which generates chunks, modifies them, and writes them to the file.
     * \param chunkWriter The filter which is writing chunks to the file. Needed to check for cancellation.
     * \param particleChunkIndex The spatial index of this file chunk, which must be filled by the time the pipeline
     *                           finishes running and is written at the end of the file chunk.
     * \param particleStreamName The name of this file chunk.
     * \param usePositionOffset If true, this is a 'PrtO' file chunk. Otherwise, this is a 'Part' file chunk.
     * \param compressionScheme The compression method for the particle chunks in this file chunk.
     */
    template <typename Pipeline>
    void write_particle_chunks( boost::shared_ptr<Pipeline> chunkPipeline,
                                boost::shared_ptr<chunk_writer> chunkWriter ) {
        // Run the pipeline.
        try {
            const std::size_t tokenCount = tbb::task_scheduler_init::default_num_threads();
            chunkPipeline->run( tokenCount );
        } catch( tbb::tbb_exception& e ) {
            if( chunkWriter->is_cancelled() ) {
                throw frantic::logging::progress_cancel_exception( e.what() );
            } else {
                // Translate into runtime_error, because Sequoia is currently swallowing tbb_exception.
                throw std::runtime_error( e.what() );
            }
        }
    }

    /**
     * Write the particle chunks from the chunk generator to the file as a 'Part' or 'PrtO' file chunk.
     *
     * \param chunkGenerator The filter which creates the particle chunks for this file chunk. The generator passes
     *                       ownership of the chunk onto this prt2_writer.
     * \param totalParticleCount The number of particles to be written. Only used for progress logging.
     * \param progress Where to log progress to.
     * \param particleStreamName The name of this file chunk.
     * \param usePositionOffset If true, this is a 'PrtO' file chunk. Otherwise, this is a 'Part' file chunk.
     * \param compressionScheme The compression method for the particle chunks in this file chunk.
     */
    void write_particle_chunks( boost::shared_ptr<tbb::filter> chunkGenerator, boost::uint64_t totalParticleCount,
                                frantic::logging::progress_logger& progress, const frantic::tstring& particleStreamName,
                                bool usePositionOffset, prt2_compression_t compressionScheme );

    /**
     * Write the particle chunks from the chunk generator to the file as a 'Part' or 'PrtO' file chunk.
     *
     * Same parameter meanings as above, but sets particleStreamName to "".
     */
    void write_particle_chunks( boost::shared_ptr<tbb::filter> chunkGenerator, boost::uint64_t totalParticleCount,
                                frantic::logging::progress_logger& progress, bool usePositionOffset,
                                prt2_compression_t compressionScheme ) {
        write_particle_chunks( chunkGenerator, totalParticleCount, progress, _T(""), usePositionOffset,
                               compressionScheme );
    }

  private:
    void write_header();

    void initialize_temp_file_state( bool useTempFile, const frantic::tstring& filename,
                                     const boost::filesystem::path& tempDir );

    template <class T>
    void write_metadata_filechunk( const frantic::tstring& metaName, const T& value ) {
        frantic::graphics::raw_byte_buffer buf;
        serialize::write_metadata_filechunk( buf, metaName, value );
        write_prt2_buffered_filechunk( buf, PRT2_FILECHUNK_NAME_Meta );
    }

    // The channel map of the particles in the file.
    channels::channel_map m_fileParticleChannelMap;

    // Stream object of the file being saved.
    std::unique_ptr<std::ostream> m_file; // In C++11, unique_ptr is preferable
    frantic::tstring m_streamname;

    // The set of particle streams that have already been written. Used for catching duplicates.
    std::set<frantic::tstring> m_writtenParticleStreamNames;

    // It's generally good to write the PRT to a temp file, then move it into place.
    bool m_useTempFile;
    boost::filesystem::path m_tempDir;
    boost::filesystem::path m_targetFile;
    boost::filesystem::path m_tempLocalFile;
    boost::filesystem::path m_tempRemoteFile;
};

namespace detail {
void write_prt2_buffered_filechunk( std::ostream& out, const frantic::graphics::raw_byte_buffer& buf,
                                    boost::uint32_t chunkName, const frantic::tstring& streamName );
} // namespace detail

} // namespace prtfile
} // namespace frantic
