// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <memory>

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/particles/particle_array.hpp>

namespace frantic {
namespace prtfile {

/**
 * This class reads in the FileChunk structure of a PRT2 file, and provides
 * functionality to read and interpret the particle index, metadata, and
 * SPRT index.
 */
class prt2_reader {
  public:
    struct particle_chunk_index_info {
        boost::int64_t chunkSeek, particleIndex;
    };

    struct particle_chunk_named {
        /** The number of particles, or -1 for uninitialized/invalid PRT file */
        boost::int64_t particleCount, particleChunkCount;
        /** The compression scheme for the particles */
        boost::uint32_t compressionScheme;

        /** The location in the file of the 'Part' or 'PrtO' FileChunk */
        boost::int64_t chunkSeek, chunkSize;
        /** If this is true, it is a 'PrtO' FileChunk instead of 'Part' */
        bool isPrtO;
        /**
         * When isPrtO is true, this provides access to the Position channel of the file,
         * so the simpler particle reading functions can fix them up and hide the details
         * of the Position+Offset particle chunks.
         */
        channels::channel_accessor<graphics::vector3f> posAccessor;

        /**
         * This contains one more entry than the number of chunks, so that
         * the number of particles in chunk i is equal to
         * particleChunkIndex[i+1].particleIndex - particleChunkIndex[i].particleIndex
         */
        std::vector<particle_chunk_index_info> particleChunkIndex;

        void clear() {
            particleCount = -1;
            particleChunkCount = -1;
            compressionScheme = 0;
            chunkSeek = -1;
            chunkSize = -1;
            isPrtO = false;
            particleChunkIndex.clear();
        }
    };

    struct unrecognized_chunk_index_info {
        boost::uint32_t chunkType;
        boost::int64_t chunkSeek, chunkSize;
    };

  private:
    particle_chunk_named m_defaultPRTChunk;
    std::map<frantic::tstring, particle_chunk_named> m_namedPRTChunks;

    /** Holds the layout of one particle */
    channels::channel_map m_particleChannelMap;

    /** Combined global and per-channel metadata */
    frantic::particles::particle_file_metadata m_metadata;

    /** Array of unrecognized chunks */
    std::vector<unrecognized_chunk_index_info> m_unrecognizedFileChunks;

    /** Stream object of the file being loaded */
    std::unique_ptr<std::istream> m_file; // In C++11, unique_ptr is preferable
    frantic::tstring m_streamname;
    bool m_isFileStream;

    /** After m_file is initialized, call this to read the file structure */
    void read_file_structure();

    /**
     * Read one chunk of particles into a raw buffer of the required bytes_size
     *
     * \param pcn  The internal structure describing the named stream of particles to read.
     * \param i    The index of the particle chunk within `pcn`.
     * \param rawBuffer  The byte-buffer into which the particles should be read. The
     *                   caller must ensure it has enough space.
     * \param outPositionOffset  If NULL, the position offset from any 'PrtO' filechunk is
     *                           added into the particles by this routine. If non-NULL, this
     *                           is populated with the position offset for 'PrtO' filechunks
     *                           or the zero vector for 'Part' filechunks.
     */
    void read_particle_chunk( const particle_chunk_named& pcn, size_t i, char* rawBuffer,
                              graphics::vector3f* outPositionOffset = NULL );
    /** Read one chunk of particles into a raw buffer, resizing it as necessary */
    void read_particle_chunk( const particle_chunk_named& pcn, size_t i, std::vector<char>& rawBuffer,
                              graphics::vector3f* outPositionOffset = NULL );
    /**
     * Read one chunk of particles into a particle_array. Changes the output's channel map
     * to match the file's being read in.
     */
    void read_particle_chunk( const particle_chunk_named& pcn, size_t i, frantic::particles::particle_array& parray,
                              graphics::vector3f* outPositionOffset = NULL );
    /** Reads all the particles from the filechunk */
    void read_particles( const particle_chunk_named& pcn, frantic::particles::particle_array& parray );

    /** Returns the range of particle indexes for the provided particle chunk index */
    static std::pair<boost::int64_t, boost::int64_t> get_particle_chunk_extents( const particle_chunk_named& pcn,
                                                                                 size_t i ) {
        const particle_chunk_index_info* pcii = &pcn.particleChunkIndex[i];
        return std::make_pair( pcii[0].particleIndex, pcii[1].particleIndex );
    }
    /** Returns the number of particles for the provided particle chunk index */
    static boost::int64_t get_particle_chunk_particle_count( const particle_chunk_named& pcn, size_t i ) {
        const particle_chunk_index_info* pcii = &pcn.particleChunkIndex[i];
        return pcii[1].particleIndex - pcii[0].particleIndex;
    }
    /**
     * Returns the number of particles for the provided range of particle chunk indices
     * irange.first <= i < irange.second.
     */
    static boost::int64_t get_particle_chunk_particle_count( const particle_chunk_named& pcn,
                                                             const std::pair<boost::int64_t, boost::int64_t>& irange ) {
        return pcn.particleChunkIndex[std::size_t( irange.second )].particleIndex -
               pcn.particleChunkIndex[std::size_t( irange.first )].particleIndex;
    }
    /** How many bytes the specified particle chunk requires */
    boost::int64_t get_particle_chunk_uncompressed_bytes_size( const particle_chunk_named& pcn, size_t i ) const {
        return get_particle_chunk_particle_count( pcn, i ) * m_particleChannelMap.structure_size();
    }

    /** Retrieves the given particle_chunk_named.  Does not do bounds checking.  If name is blank, gets the default. */
    const particle_chunk_named& get_particle_chunk_named( const frantic::tstring& name ) const {
        if( name.empty() )
            return m_defaultPRTChunk;
        return m_namedPRTChunks.find( name )->second;
    }

    /**
     * re-open a closed file stream
     */
    void reopen_stream();

    /**
     * Check and open a closed file stream
     */
    void auto_open_stream();

  public:
    /** Initialize in empty state */
    prt2_reader() { m_defaultPRTChunk.clear(); }

    /** Initialize by opening a file */
    prt2_reader( const frantic::tstring& filename ) {
        m_defaultPRTChunk.clear();
        open( filename );
    }

    /** Reads in the file structure, without reading in any particles */
    void open( const frantic::tstring& filename );

    /**
     * Reads in the file structure, without reading in any particles
     * Ownership of the 'is' istream passes into the prt2_reader object.
     */
    void open( std::istream* is, const frantic::tstring& filename );

    void close();

    /**
     * Closes the stream without discarding loaded information about the file
     */
    void close_stream();

    const frantic::tstring& get_stream_name() const { return m_streamname; }

    /** Returns a reference to the vector of unrecognized filechunks */
    const std::vector<unrecognized_chunk_index_info>& get_unrecognized_filechunks() const {
        return m_unrecognizedFileChunks;
    }

    /** Removes all chunks matching 'chunkType' from the unrecognized chunk vector, appending them to 'outChunkInfo' */
    void remove_unrecognized_filechunks( boost::uint32_t chunkType,
                                         std::vector<unrecognized_chunk_index_info>& outChunkInfo );

    /** Reads the data for the  */
    void read_unrecognized_filechunk( const unrecognized_chunk_index_info& unrecognizedChunk, char* outBuffer );

    /** The total number of particles in the file in the given chunk */
    boost::int64_t get_particle_count( const frantic::tstring& name = _T("") ) const {
        return get_particle_chunk_named( name ).particleCount;
    }
    /** The number of named channels per particle */
    size_t get_channel_count() const { return m_particleChannelMap.channel_count(); }
    /** Whether the Position channel is stored in Position+Offset mode */
    bool is_position_offset() const;
    /** The channel map which holds the channel names, types, and layout */
    const channels::channel_map& get_channel_map() const { return m_particleChannelMap; }
    /** The combined general and per-channel metadata */
    const particles::particle_file_metadata& get_metadata() const;
    /** The property map holding the general metadata */
    const channels::property_map& get_general_metadata() const;
    /** Gets the property map for one channel by name (this is redundant, it's also in m_metadata) */
    const channels::property_map& get_channel_metadata( const frantic::tstring& channelName ) const;

    /** Gets the name of all the particle chunks and places them in the given vector */
    void get_particle_chunk_names( std::vector<frantic::tstring>& outNames ) const {
        for( std::map<frantic::tstring, particle_chunk_named>::const_iterator iter = m_namedPRTChunks.begin();
             iter != m_namedPRTChunks.end(); ++iter ) {
            outNames.push_back( iter->first );
        }
    }
    /** Returns true if the given chunk name is contained in this prt2 file */
    bool contains_particle_chunk_name( const frantic::tstring& name ) const {
        return m_namedPRTChunks.find( name ) != m_namedPRTChunks.end();
    }

    /** The total number of particle chunks in the default 'Part' file chunk */
    boost::int64_t get_particle_chunk_count( const frantic::tstring& name = _T("") ) const {
        return get_particle_chunk_named( name ).particleChunkCount;
    }
    /** Returns the range of particle indexes for the provided particle chunk index in the given chunk */
    std::pair<boost::int64_t, boost::int64_t>
    get_particle_chunk_extents( size_t i, const frantic::tstring& name = _T("") ) const {
        return get_particle_chunk_extents( get_particle_chunk_named( name ), i );
    }
    /** Returns the number of particles for the provided particle chunk index in the given chunk */
    boost::int64_t get_particle_chunk_particle_count( size_t i, const frantic::tstring& name = _T("") ) const {
        return get_particle_chunk_particle_count( get_particle_chunk_named( name ), i );
    }
    /**
     * Returns the number of particles for the provided particle chunk range of indices
     * irange.first <= i < irange.second.
     */
    boost::int64_t get_particle_chunk_particle_count( const std::pair<boost::int64_t, boost::int64_t>& irange,
                                                      const frantic::tstring& name = _T("") ) const {
        return get_particle_chunk_particle_count( get_particle_chunk_named( name ), irange );
    }

    /** How many bytes the specified particle chunk requires for the given chunk */
    boost::int64_t get_particle_chunk_uncompressed_bytes_size( size_t i, const frantic::tstring& name = _T("") ) const {
        return get_particle_chunk_uncompressed_bytes_size( get_particle_chunk_named( name ), i );
    }

    /** Read one chunk of particles into a raw buffer of the required bytes_size from the given chunk */
    void read_particle_chunk( size_t i, char* rawBuffer, graphics::vector3f* outPositionOffset = NULL,
                              const frantic::tstring& name = _T("") ) {
        read_particle_chunk( get_particle_chunk_named( name ), i, rawBuffer, outPositionOffset );
    }
    /** Read one chunk of particles into a raw buffer from the given chunk, resizing it as necessary */
    void read_particle_chunk( size_t i, std::vector<char>& rawBuffer, graphics::vector3f* outPositionOffset = NULL,
                              const frantic::tstring& name = _T("") ) {
        read_particle_chunk( get_particle_chunk_named( name ), i, rawBuffer, outPositionOffset );
    }
    /**
     * Read one chunk of particles into a particle_array from the given chunk.
     * Changes the output's channel map to match the file's being read in.
     */
    void read_particle_chunk( size_t i, frantic::particles::particle_array& parray,
                              graphics::vector3f* outPositionOffset = NULL, const frantic::tstring& name = _T("") ) {
        read_particle_chunk( get_particle_chunk_named( name ), i, parray, outPositionOffset );
    }

    /** Reads all the particles from the default 'Part' chunk in the given chunk */
    void read_particles( frantic::particles::particle_array& parray, const frantic::tstring& name = _T("") ) {
        read_particles( get_particle_chunk_named( name ), parray );
    }
};

} // namespace prtfile
} // namespace frantic
