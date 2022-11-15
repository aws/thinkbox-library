// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/move/core.hpp>

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <frantic/graphics/raw_byte_buffer.hpp>

#include <frantic/volumetrics/voxel_coord_system.hpp>
#include <frantic/volumetrics/voxel_edge_stepper.hpp>

#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/particles/streams/particle_ostream.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/particle_cursor.hpp>

#include <frantic/logging/logging_level.hpp>
#include <frantic/logging/progress_logger.hpp>

#include <frantic/particles/particle_array_iterator.hpp>

namespace frantic {
namespace particles {

using frantic::graphics::raw_byte_buffer;

class particle_array {
    // Mark this class as copyable and movable.
    BOOST_COPYABLE_AND_MOVABLE( particle_array )

    frantic::channels::channel_map m_channelMap;

    raw_byte_buffer m_particleBuffer; // contains the particle raw data

    friend class particle_array_iterator;

  public:
    typedef particle_array_iterator iterator;
    typedef const_particle_array_iterator const_iterator;

    ///////////
    // Constructors
    ///////////

    // Construct an empty particle array, with a channel map consisting of just position.
    particle_array() {
        m_channelMap.define_channel( _T("Position"), 3, frantic::channels::data_type_float32 );
        m_channelMap.end_channel_definition();
        m_particleBuffer.reserve( 1 );
    }

    // Constructs an empty particle array, with the provided channel map
    explicit particle_array( const frantic::channels::channel_map& map ) {
        // if( map.channel_count() == 0 ) {
        //	throw std::runtime_error( "Cannot create a particle array with an empty channel map.");
        // }

        // DH 1/17/2013 I'm allowing an empty channel_map, though this is guaranteed to cause problems if you try to add
        // to it, or resize it.

        m_channelMap = map;
        if( m_channelMap.channel_count() > 0 )
            m_particleBuffer.reserve( 1 );
    }

    particle_array( const particle_array& other ) {
        m_channelMap = other.m_channelMap;
        m_particleBuffer = other.m_particleBuffer;
    }

    particle_array( BOOST_RV_REF( particle_array ) other ) { swap( other ); }

    particle_array& operator=( BOOST_COPY_ASSIGN_REF( particle_array ) other ) {
        particle_array copy( other );
        swap( copy );
        return *this;
    }

    particle_array& operator=( BOOST_RV_REF( particle_array ) other ) {
        swap( other );
        return *this;
    }

    const frantic::channels::channel_map& get_channel_map() const { return m_channelMap; }

    void reset( const frantic::channels::channel_map& map ) {
        m_channelMap = map;
        m_particleBuffer.clear();
    }

    bool empty() const { return m_particleBuffer.empty(); }

    void resize( const std::size_t numElements ) {
        m_particleBuffer.resize( numElements * m_channelMap.structure_size() );
    }

    void reserve( const std::size_t numElements ) {
        m_particleBuffer.reserve( numElements * m_channelMap.structure_size() );
    }

    void trim() { m_particleBuffer.trim(); }

    char* release() { return m_particleBuffer.release(); }

    void clear() { m_particleBuffer.clear(); }

    void swap( particle_array& rhs ) {
        m_particleBuffer.swap( rhs.m_particleBuffer );
        m_channelMap.swap( rhs.m_channelMap );
    }

    iterator begin() { return particle_array_iterator( m_particleBuffer.begin(), m_channelMap.structure_size() ); }

    iterator end() { return particle_array_iterator( m_particleBuffer.end(), m_channelMap.structure_size() ); }

    const_iterator begin() const {
        return const_particle_array_iterator( m_particleBuffer.begin(), m_channelMap.structure_size() );
    }

    const_iterator end() const {
        return const_particle_array_iterator( m_particleBuffer.end(), m_channelMap.structure_size() );
    }

    char* at( std::size_t i ) { return m_particleBuffer.ptr_at( i * m_channelMap.structure_size() ); }

    const char* at( std::size_t i ) const { return m_particleBuffer.ptr_at( i * m_channelMap.structure_size() ); }

    char* operator[]( std::size_t i ) { return m_particleBuffer.ptr_at( i * m_channelMap.structure_size() ); }

    const char* operator[]( std::size_t i ) const {
        return m_particleBuffer.ptr_at( i * m_channelMap.structure_size() );
    }

    bool equals( const particle_array& rhs ) const {
        if( m_channelMap != rhs.m_channelMap )
            return false;

        return m_particleBuffer == rhs.m_particleBuffer;
    }

    /**
     * Adds a particle to the end of the array.
     *
     * @param  rawParticleData  This is the memory for the particle to add to the end of the array.
     */
    void push_back( const char* rawParticleData ) {
        m_particleBuffer.add_element( rawParticleData, m_channelMap.structure_size() );
    }

    /**
     * Adds a particle to the end of the array, using the provided channel map adaptor to
     * convert to the particle_array's format.
     *
     * The caller must ensure that the destination of the pcmAdaptor is the same as the particle_array's channel_map.
     */
    void push_back( const char* rawParticleData, const frantic::channels::channel_map_adaptor& pcmAdaptor ) {
        pcmAdaptor.copy_structure( m_particleBuffer.add_element( m_channelMap.structure_size() ), rawParticleData );
    }

    void pop_back() { m_particleBuffer.resize( m_particleBuffer.size() - m_channelMap.structure_size() ); }

    /**
     * Returns the number of the particles contained within the array.
     */
    std::size_t particle_count() const { return m_particleBuffer.size() / m_channelMap.structure_size(); }

    /**
     * Returns the number of the particles contained within the array.
     */
    std::size_t size() const { return m_particleBuffer.size() / m_channelMap.structure_size(); }

    /**
     * Returns the size of the particle buffer in bytes.
     */
    std::size_t size_in_memory() const { return m_particleBuffer.size(); }

    /**
     * Returns the how large the array may grow to before reallocating the underlying memory
     */
    std::size_t capacity() const { return m_particleBuffer.capacity() / m_channelMap.structure_size(); }

    /**
     * Inserts a particle_istream's entire stream into a particle array
     */
    void insert_particles( boost::shared_ptr<particles::streams::particle_istream> pin ) {
        frantic::logging::null_progress_logger nullProgress;
        insert_particles( pin, NULL, nullProgress );
    }

    /**
     * Inserts a particle_istream's entire stream into a particle array
     * Overload with progress logging.
     */
    void insert_particles( boost::shared_ptr<particles::streams::particle_istream> pin,
                           frantic::logging::progress_logger& progress ) {
        insert_particles( pin, NULL, progress );
    }

    /**
     * Inserts a particle_istream's entire stream into a particle array
     * Overload with progress logging and default particle.
     * Useful to set default values for channels that do not exist within the stream.
     */
    void insert_particles( boost::shared_ptr<particles::streams::particle_istream> pin, char* defaultParticle,
                           frantic::logging::progress_logger& progress ) {
        static const std::size_t CHUNK_SIZE = 50000; // Arbitrarily process at most this many particles at a time.
        boost::int64_t streamCount = pin->particle_count_left();

        // Explicity set the stream's channel map to match this array's channel map.
        // NOTE: I feel justified changing the incoming stream object because we intend
        // to exhaust it during this call; hence we don't need its channel map to remain unchanged.
        pin->set_channel_map( m_channelMap );
        if( defaultParticle )
            pin->set_default_particle( defaultParticle );

        // If the stream knows how many particles it has left, reserve the memory for it all at once
        // Then process the particles in CHUNK_SIZE chunks.
        if( streamCount > 0 ) {
            std::size_t offset = size();
            std::size_t expectedCount = CHUNK_SIZE;

            resize( offset + (std::size_t)streamCount );
            while( pin->get_particles( at( offset ), expectedCount ) ) {
                offset += expectedCount;
                expectedCount = CHUNK_SIZE;

                progress.update_progress( offset, size() );
            }

            if( ( offset + expectedCount ) != size() ) {
                std::stringstream ss;
                ss << "particle_array.insert_particles() - Did not recieve the expected number: " << size()
                   << " of particles.\n";
                ss << "Instead it received: " << ( offset + expectedCount ) << " particles.\n";
                ss << "The stream still has: " << pin->particle_count_left() << " particles." << std::endl;
                throw std::runtime_error( ss.str() );
            }

            progress.update_progress( offset, size() );
        }
        // Otherwise, a negative number indicates the stream doesn't know how many particles it has, so copy them in
        // CHUNK_SIZE chunks.
        else if( streamCount < 0 ) {
            std::size_t offset = size();
            std::size_t expectedCount = CHUNK_SIZE;

            boost::int64_t progressVal;
            boost::int64_t progressTotal = pin->particle_progress_count();

            // We need to use add_element() because it has memory doubling behaviour.
            char* p = m_particleBuffer.add_element( m_channelMap.structure_size() * CHUNK_SIZE, 1 );
            while( pin->get_particles( p, expectedCount ) ) {
                offset += expectedCount;
                expectedCount = CHUNK_SIZE;

                // Truncate unused particles from the end, then grow the vector for the next batch
                resize( offset );
                p = m_particleBuffer.add_element( m_channelMap.structure_size() * CHUNK_SIZE );

                // Make sure that progressTotal is greater than progressVal. This is usually only false
                // when loading .CSV files.
                progressVal = pin->particle_progress_index();
                if( progressVal < progressTotal )
                    progress.update_progress( progressVal, progressTotal );
            }

            resize( offset + expectedCount );
        }
    }

    template <class ForwardParticleIterator>
    void insert_particles( const frantic::channels::channel_map& cm, ForwardParticleIterator iter,
                           ForwardParticleIterator iterEnd ) {
        // make a converter to go to from the incoming channel map to the array's channel map
        frantic::channels::channel_map_adaptor cma( m_channelMap, cm );
        std::size_t destinationSize = m_channelMap.structure_size();

        if( cma.is_identity() ) {
            for( ; iter != iterEnd; ++iter )
                memcpy( m_particleBuffer.add_element( destinationSize ), *iter, destinationSize );
        } else {
            char* particle;
            for( ; iter != iterEnd; ++iter ) {
                particle = m_particleBuffer.add_element( destinationSize );

                memset( particle, 0, destinationSize );
                cma.copy_structure( particle, *iter );
            }
        }
    }

    void write_particles( boost::shared_ptr<particles::streams::particle_ostream> pout ) const {
        pout->set_channel_map( m_channelMap );
        for( const_iterator i = begin(), endIter = end(); i != endIter; ++i )
            pout->put_particle( *i );
    }

    /**
     * Copies all of the particles from a given particle array.
     *
     * @param  sourceParticleArray This is the particle array from which we are copying.
     */
    void copy_particles_from( const frantic::particles::particle_array& sourceParticleArray ) {
        frantic::channels::channel_map_adaptor cma( get_channel_map(), sourceParticleArray.get_channel_map() );

        std::vector<char> defaultData( cma.dest_size(), 0 );
        m_particleBuffer.reserve( ( particle_count() + sourceParticleArray.particle_count() ) *
                                  get_channel_map().structure_size() );
        for( std::size_t i = 0, ie = sourceParticleArray.particle_count(); i != ie; ++i ) {
            cma.copy_structure( m_particleBuffer.add_element( cma.dest_size() ), sourceParticleArray.at( i ),
                                &defaultData[0] );
        }
    }

    /**
     * Checks whether this particle_array has the specified named channel.
     *
     * @param channelName The channel name to check for.
     * @return true if the particle_array has the specified named channel,
     *		or false otherwise.
     */
    bool has_channel( const frantic::tstring& channelName ) { return m_channelMap.has_channel( channelName ); }

    /**
     * Sets the channel map to the input channel map, and updates all of the particle data while maintaining strong
     * exception safety.
     *
     * @param  channelMap  This is the new channel map for the particle array.
     */
    void set_channel_map( const frantic::channels::channel_map& channelMap ) {
        if( channelMap.channel_count() == 0 ) {
            throw std::runtime_error( "Cannot create a particle array with an empty channel map." );
        }
        particle_array pArray( channelMap );
        if( size() > 0 ) {
            pArray.copy_particles_from( *this );
        }
        swap( pArray );
    }

    /**
     * Removes a channel from the particle_array by creating a new particle_array and copying over all of the data. This
     * is slow, but has strong exception safety.")
     *
     * @param  channelName This is the channel to be removed.
     */
    void remove_channel( const frantic::tstring& channelName ) {
        const channel_map& cmIn = get_channel_map();
        channel_map cmOut = get_channel_map();
        cmOut.delete_channel( channelName );

        if( cmOut.channel_count() == 0 ) {
            throw std::runtime_error(
                "Cannot remove the last channel (" + frantic::strings::to_string( channelName ) +
                ") from the particle array because empty channel maps are not allowed. If you would like to get rid of "
                "this "
                "channel, try using set_channel_map() to replace it with a (non-empty) channel_map." );
        }

        channels::channel_map_adaptor cma( cmOut, cmIn );
        particle_array tempArray( cmOut );
        tempArray.reserve( particle_count() );
        tempArray.copy_particles_from( *this );

        swap( tempArray );
    }

    /**
     * Appends `count` particles from the `source` buffer without doing any checking of the particle channel map.
     *
     * \pre The caller should have confirmed that the channel map of the data in `source` matches that of the
     *      `particle_array`.
     *
     * \param source  The buffer with particles to append.
     * \param count   The number of particles to append.
     */
    void unchecked_insert_particles( const char* source, size_t count ) {
        size_t size = m_channelMap.structure_size() * count;
        memcpy( m_particleBuffer.add_element( size ), source, size );
    }
};
} // namespace particles
} // namespace frantic
