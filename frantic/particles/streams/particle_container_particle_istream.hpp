// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once


#include <boost/scoped_array.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <boost/smart_ptr.hpp>

namespace frantic {
namespace particles {
namespace streams {

template <class ForwardIterator>
class particle_container_particle_istream : public particle_istream {
    frantic::tstring m_name;
    frantic::channels::channel_map m_sourceMap;
    frantic::channels::channel_map m_outputMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    boost::scoped_array<char> m_defaultParticle;

    boost::int64_t m_particleIndex;
    boost::int64_t m_particleCount;
    ForwardIterator m_curIter, m_endIter;

    // Private copy constructor to disable copying
    particle_container_particle_istream( const particle_container_particle_istream& rhs );

    // Private assignment operator to disable assignment
    particle_container_particle_istream& operator=( const particle_container_particle_istream& );

  public:
    /**
     * @param itBegin	The iterator to start the stream at.
     * @param itEnd		The iterator which marks the end of the stream.
     * @param particleMap	The channel map for the stream to use.
     * @param name	The name of the particle stream.
     */
    particle_container_particle_istream( ForwardIterator itBegin, ForwardIterator itEnd,
                                         const frantic::channels::channel_map& particleMap,
                                         const frantic::tstring& name = _T("<particle container>") )
        : m_name( name )
        , m_curIter( itBegin )
        , m_endIter( itEnd )
        , m_sourceMap( particleMap )
        , m_outputMap( particleMap ) {
        // Start the particle index one before the first valid index
        m_particleIndex = -1;
        // This might not be valid for forward iterators, so it might have to be a parameter.
        m_particleCount = ( itEnd - itBegin );

        m_pcmAdaptor.set( m_outputMap, m_outputMap );
    }

    /**
     * Constructor that allows specifying the particle count for containers which cannot have their size inferred
     * from their iterators, such as a particle_grid_tree.
     *
     * @param itBegin	The iterator to start the stream at.
     * @param itEnd		The iterator which marks the end of the stream.
     * @param particleCount	The number of particles in the container.
     * @param particleMap	The channel map for the stream to use.
     * @param name	The name of the particle stream.
     */
    particle_container_particle_istream( ForwardIterator itBegin, ForwardIterator itEnd, boost::int64_t particleCount,
                                         const frantic::channels::channel_map& particleMap,
                                         const frantic::tstring& name = _T("<particle container>") )
        : m_name( name )
        , m_curIter( itBegin )
        , m_endIter( itEnd )
        , m_particleCount( particleCount )
        , m_sourceMap( particleMap )
        , m_outputMap( particleMap ) {
        // Start the particle index one before the first valid index
        m_particleIndex = -1;

        m_pcmAdaptor.set( m_outputMap, m_outputMap );
    }

    particle_container_particle_istream( ForwardIterator itBegin, ForwardIterator itEnd,
                                         const frantic::channels::channel_map& particleMap,
                                         const frantic::channels::channel_map& outputMap,
                                         const std::string& name = "<particle container>" )
        : m_name( name )
        , m_curIter( itBegin )
        , m_endIter( itEnd )
        , m_sourceMap( particleMap )
        , m_outputMap( outputMap ) {
        // Start the particle index one before the first valid index
        m_particleIndex = -1;
        // This might not be valid for forward iterators, so it might have to be a parameter.
        m_particleCount = ( itEnd - itBegin );

        m_pcmAdaptor.set( m_outputMap, m_sourceMap );
    }

    virtual ~particle_container_particle_istream() {}

    void close() {}

    std::size_t particle_size() const { return m_outputMap.structure_size(); }

    frantic::tstring name() const { return m_name; }

    boost::int64_t particle_count() const { return m_particleCount; }

    boost::int64_t particle_index() const { return m_particleIndex; }

    boost::int64_t particle_count_left() const { return m_particleCount - m_particleIndex - 1; }

    boost::int64_t particle_progress_count() const { return m_particleCount; }

    boost::int64_t particle_progress_index() const { return m_particleIndex; }

    const frantic::channels::channel_map& get_channel_map() const { return m_outputMap; }

    // Access to the channel_map in the file
    const frantic::channels::channel_map& get_native_channel_map() const { return m_sourceMap; }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        if( m_defaultParticle ) {
            boost::scoped_array<char> newDefault( new char[particleChannelMap.structure_size()] );
            memset( newDefault.get(), 0, particleChannelMap.structure_size() );

            frantic::channels::channel_map_adaptor tempAdaptor( m_outputMap, particleChannelMap );
            tempAdaptor.copy_structure( newDefault.get(), m_defaultParticle.get() );

            m_defaultParticle.swap( newDefault );
        }

        // Set the map and the adaptor
        m_outputMap = particleChannelMap;
        m_pcmAdaptor.set( m_outputMap, m_sourceMap );
    }

    void set_default_particle( char* buffer ) {
        if( !m_defaultParticle )
            m_defaultParticle.reset( new char[m_outputMap.structure_size()] );
        memcpy( m_defaultParticle.get(), buffer, m_outputMap.structure_size() );
    }

    bool get_particle( char* rawParticleBuffer ) {
        if( m_curIter == m_endIter )
            return false;
        if( m_pcmAdaptor.is_identity() ) {
            m_outputMap.copy_structure( rawParticleBuffer, *m_curIter );
        } else {
            if( !m_defaultParticle.get() ) {
                memset( rawParticleBuffer, 0, m_outputMap.structure_size() );
                m_pcmAdaptor.copy_structure( rawParticleBuffer, *m_curIter );
            } else {
                m_pcmAdaptor.copy_structure( rawParticleBuffer, *m_curIter, m_defaultParticle.get() );
            }
        }

        ++m_curIter;
        ++m_particleIndex;
        return true;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        for( std::size_t i = 0; i < numParticles; ++i, buffer += m_outputMap.structure_size() ) {
            if( !get_particle( buffer ) ) {
                numParticles = i;
                return false;
            }
        }

        return true;
    }
};

template <class ForwardIterator>
inline boost::shared_ptr<particle_istream>
make_particle_container_particle_istream( ForwardIterator itBegin, ForwardIterator itEnd,
                                          const frantic::channels::channel_map& sourceMap,
                                          const frantic::tstring& name = _T("<particle container>") ) {
    return boost::shared_ptr<particle_istream>(
        new particle_container_particle_istream<ForwardIterator>( itBegin, itEnd, sourceMap, name ) );
}

} // namespace streams
} // namespace particles
} // namespace frantic
