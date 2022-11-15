// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/scoped_array.hpp>

#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/prtfile/prt2_reader.hpp>

namespace frantic {
namespace particles {
namespace streams {

class prt2_particle_istream : public frantic::particles::streams::particle_istream {
  public:
    struct chunk_entry {
        frantic::tstring m_chunkName;
        std::vector<std::pair<boost::int64_t, boost::int64_t>> m_chunkIndices;

        explicit chunk_entry( const frantic::tstring& chunkName = _T("") )
            : m_chunkName( chunkName ) {}

        static chunk_entry all( const frantic::tstring& chunkName = _T("") ) {
            chunk_entry entry( chunkName );
            entry.m_chunkIndices.push_back(
                std::make_pair( boost::int64_t( 0 ), std::numeric_limits<boost::int64_t>::max() ) );
            return entry;
        }

        void add_index( boost::int64_t index ) {
            if( !m_chunkIndices.empty() ) {
                if( m_chunkIndices.back().second + 1 == index ) {
                    m_chunkIndices.back().second = index;
                    return;
                }
            }
            m_chunkIndices.push_back( std::make_pair( index, index ) );
        }

        void add_indices( boost::int64_t index, boost::int64_t count ) {
            if( count <= 0 )
                return;
            if( !m_chunkIndices.empty() ) {
                if( m_chunkIndices.back().second + 1 == index ) {
                    m_chunkIndices.back().second = index + count - 1;
                    return;
                }
            }
            m_chunkIndices.push_back( std::make_pair( index, index + count - 1 ) );
        }
    };
    typedef std::vector<chunk_entry> chunk_list_t;

  private:
    boost::shared_ptr<frantic::prtfile::prt2_reader> m_prt2;
    frantic::tstring m_streamName;

    channels::channel_map m_particleChannelMap;
    channels::channel_map_adaptor m_pcmAdaptor;

    boost::scoped_array<char> m_defaultParticle;
    boost::int64_t m_currentParticleIndex;

    /** Buffer that holds one chunk read from the PRT2 file */
    std::vector<char> m_chunkBuffer;

    /** List of chunks to read */
    chunk_list_t m_chunksList;
    boost::uint64_t m_currentChunkIndex;
    boost::uint64_t m_currentChunkEntryIndex;

    /** Which chunk is next for m_chunkBuffer */
    boost::int64_t m_nextChunkIndex;
    /** How large m_chunkBuffer is*/
    boost::int64_t m_nextParticleIndexInChunk, m_chunkParticleCount;

    /** Precalculated particle count */
    boost::int64_t m_particleCount;

    /** In float64 output with Position+Offset input, we apply the offset in float64 */
    bool m_isFloat64PositionOffset;
    frantic::channels::channel_accessor<frantic::graphics::vector3fd> m_pos64Accessor;
    frantic::graphics::vector3f m_pos64Offset;

    /** If necessary, reads a new chunk. Returns true if particles are still available */
    bool chunk_entry_advance();
    /** If necessary, reads a new chunk. Returns true if particles are still available */
    bool chunk_advance();

    // Private copy constructor and assignment operator to disable copying
    prt2_particle_istream( const prt2_particle_istream& );            // not implemented
    prt2_particle_istream& operator=( const prt2_particle_istream& ); // not implemented

    void init();

    /**
     * Given the channel map from the file, possibly increase/decrease the Position precision depending on the
     * requested position type.
     */
    void set_channel_map_with_position_type( const frantic::channels::channel_map& pcm,
                                             frantic::channels::data_type_t positionTypeHint );

  public:
    prt2_particle_istream( const frantic::tstring& file,
                           frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid,
                           const frantic::tstring& chunkName = _T("") );
    prt2_particle_istream( const frantic::tstring& file, const frantic::channels::channel_map& particleChannelMap,
                           const frantic::tstring& chunkName = _T("") );

    prt2_particle_istream( const frantic::tstring& file, const chunk_entry& chunk,
                           frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid );
    prt2_particle_istream( const frantic::tstring& file, const frantic::channels::channel_map& particleChannelMap,
                           const chunk_entry& chunk );

    prt2_particle_istream( const frantic::tstring& file, const chunk_list_t& chunks,
                           frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid );
    prt2_particle_istream( const frantic::tstring& file, const frantic::channels::channel_map& particleChannelMap,
                           const chunk_list_t& chunks );

    // Constructors for reading from a single shared reader (to reduce the number of file handles required).
    // Assumes file is already opened for reading

    prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                           frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid,
                           const frantic::tstring& chunkName = _T("") );
    prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                           const frantic::channels::channel_map& particleChannelMap,
                           const frantic::tstring& chunkName = _T("") );

    prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file, const chunk_entry& chunk,
                           frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid );
    prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                           const frantic::channels::channel_map& particleChannelMap, const chunk_entry& chunk );
    prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file, const chunk_list_t& chunks,
                           frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid );

    virtual ~prt2_particle_istream();

    void close();

    //////////////////////
    // Information functions
    //////////////////////

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_streamName; }

    boost::int64_t particle_count() const { return m_particleCount; }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return m_particleCount - m_currentParticleIndex - 1; }

    boost::int64_t particle_progress_count() const { return m_particleCount; }

    boost::int64_t particle_progress_index() const { return m_currentParticleIndex; }

    const channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    // Access to the channel_map in the file
    const channels::channel_map& get_native_channel_map() const { return m_prt2->get_channel_map(); }

    const frantic::channels::property_map& get_general_metadata() const { return m_prt2->get_general_metadata(); }

    const frantic::channels::property_map* get_channel_metadata( const frantic::tstring& channelName ) const {
        if( m_prt2->get_channel_map().has_channel( channelName ) ) {
            return &m_prt2->get_channel_metadata( channelName );
        } else {
            return NULL;
        }
    }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    void set_channel_map( const channels::channel_map& particleChannelMap );

    void set_default_particle( char* rawParticleBuffer );

    /**
     * This reads a particle into a buffer matching the channel_map.  It returns true if a particle was read, false
     * otherwise.
     * IMPORTANT: Make sure the buffer you pass in is at least as big as particle_size() bytes.
     */
    bool get_particle( char* rawParticleBuffer );

    bool get_particles( char* particleBuffer, std::size_t& numParticles );
};

} // namespace streams
} // namespace particles
} // namespace frantic
