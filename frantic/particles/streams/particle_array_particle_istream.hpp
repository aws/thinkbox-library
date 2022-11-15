// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once


#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <frantic/particles/particle_array.hpp>

namespace frantic {
namespace particles {
namespace streams {

using frantic::channels::channel_map;
using frantic::channels::channel_map_adaptor;

class particle_array_particle_istream : public particle_istream {
    frantic::tstring m_name;
    channel_map m_particleChannelMap;
    channel_map_adaptor m_pcmAdaptor;
    boost::int64_t m_currentParticleIndex;
    std::vector<char> m_defaultParticleBuffer;

    const frantic::particles::particle_array& m_particles;

    // Private copy constructor to disable copying
    particle_array_particle_istream( const particle_array_particle_istream& rhs )
        : m_particles( rhs.m_particles ) {}

    // Private assignment operator to disable assignment
    particle_array_particle_istream& operator=( const particle_array_particle_istream& ); // not implemented

  public:
    // This loads a particle file, and feeds back the results in the particle's native particle storage layout.
    particle_array_particle_istream( const frantic::particles::particle_array& particles,
                                     const frantic::tstring& name = _T("<particle array>") )
        : m_name( name )
        , m_particles( particles ) {
        // Start the particle index one before the first valid index
        m_currentParticleIndex = -1;
        // Initialize the particle channel map that the user of this class sees to match the native one provided by the
        // .csv
        set_channel_map( m_particles.get_channel_map() );
    }

    particle_array_particle_istream( const frantic::particles::particle_array& particles,
                                     const channel_map& particleChannelMap,
                                     const frantic::tstring& name = _T("<particle array>") )
        : m_name( name )
        , m_particles( particles ) {
        m_currentParticleIndex = -1;
        set_channel_map( particleChannelMap );
    }

    virtual ~particle_array_particle_istream() {}

    void close() {}

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_name; }

    boost::int64_t particle_count() const { return m_particles.size(); }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return particle_count() - particle_index() - 1; }

    boost::int64_t particle_progress_count() const { return m_particles.size(); }

    boost::int64_t particle_progress_index() const { return m_currentParticleIndex; }

    const channel_map& get_channel_map() const { return m_particleChannelMap; }

    // Access to the channel_map in the file
    const channel_map& get_native_channel_map() const { return m_particles.get_channel_map(); }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const channel_map& particleChannelMap ) {
        std::vector<char> newDefaultParticle( particleChannelMap.structure_size() );
        if( m_defaultParticleBuffer.size() > 0 ) {
            frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
            defaultAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticleBuffer[0] );
        } else
            memset( &newDefaultParticle[0], 0, particleChannelMap.structure_size() );
        m_defaultParticleBuffer.swap( newDefaultParticle );

        // Set the map and the adaptor
        m_particleChannelMap = particleChannelMap;
        m_pcmAdaptor.set( m_particleChannelMap, m_particles.get_channel_map() );
    }

    void set_default_particle( char* buffer ) {
        memcpy( &m_defaultParticleBuffer[0], buffer, m_particleChannelMap.structure_size() );
    }

    bool get_particle( char* rawParticleBuffer ) {
        // Detect the end of the array
        if( particle_index() + 1 >= particle_count() )
            return false;

        ++m_currentParticleIndex;

        if( m_pcmAdaptor.is_identity() ) {
            // If we can, dump the data straight to the output
            m_particles.get_channel_map().copy_structure( rawParticleBuffer,
                                                          m_particles.at( (std::size_t)m_currentParticleIndex ) );
        } else {
            // Otherwise use the adaptor to convert the particles to the appropriate layout, filling missing channels
            // with the defaults
            m_pcmAdaptor.copy_structure( rawParticleBuffer, m_particles[(std::size_t)m_currentParticleIndex],
                                         &m_defaultParticleBuffer[0] );
        }

        return true;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        for( std::size_t i = 0; i < numParticles; ++i, buffer += m_particleChannelMap.structure_size() ) {
            if( !get_particle( buffer ) ) {
                numParticles = i;
                return false;
            }
        }

        return true;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
