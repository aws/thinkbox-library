// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_kdtree.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

class instancing_particle_istream : public delegated_particle_istream {
    struct instance_particle {
        vector3f position;
        float scale;
        int id;

        instance_particle()
            : id( -1 ) {}
        operator const vector3f&() const { return position; }
        operator vector3f&() { return position; }
        float& operator[]( int i ) { return position[i]; }
        float operator[]( int i ) const { return position[i]; }
    };

    boost::scoped_array<char> m_delegateParticle;
    std::vector<instance_particle> m_instanceParticles;
    std::size_t m_instanceIndex;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccessor;

    void load_instancing_particles( boost::shared_ptr<particle_istream> pin, float scaleMinimum,
                                    float scaleMultiplier ) {
        frantic::channels::channel_map inMap;
        inMap.define_channel( _T("Position"), 3, frantic::channels::data_type_float32, 0 );
        inMap.define_channel( _T("Scale"), 1, frantic::channels::data_type_float32, 12 );
        inMap.define_channel( _T("ID"), 1, frantic::channels::data_type_int32, 16 );
        inMap.end_channel_definition( 1, true, false );
        pin->set_channel_map( inMap );

        const int CHUNK_SIZE = 50000;

        if( pin->particle_count() > 0 )
            m_instanceParticles.reserve( (std::size_t)pin->particle_count() );
        std::size_t counter = 0, expectedCount = 0;

        do {
            counter += expectedCount;
            expectedCount = CHUNK_SIZE;

            m_instanceParticles.resize( counter + CHUNK_SIZE );
        } while( pin->get_particles( (char*)&m_instanceParticles[counter], expectedCount ) );
        m_instanceParticles.resize( counter + expectedCount );

        if( scaleMultiplier > 0 ) {
            frantic::particles::particle_kdtree<instance_particle> kdTree;
            kdTree.set_particles_with_swap( m_instanceParticles );
            kdTree.balance_kdtree();

            struct not_equal_id {
                int id;
                not_equal_id( int _id )
                    : id( _id ) {}
                bool operator()( const instance_particle& p ) const { return p.id != id; }
            };

            // For each particle, set its id to prevent queries from returning itself. Then set
            // the particle's scale to a linear function of the distance to the nearest particle.
            for( std::size_t i = 0; i < kdTree.size(); ++i ) {
                const instance_particle& queryParticle = kdTree[i];
                const_cast<instance_particle&>( queryParticle ).id = (int)i;

                const instance_particle& closeParticle =
                    kdTree.locate_closest_particle_if( queryParticle.position, not_equal_id( queryParticle.id ) );
                const_cast<instance_particle&>( queryParticle ).scale =
                    scaleMinimum +
                    scaleMultiplier * vector3f::distance( queryParticle.position, closeParticle.position );
            }

            kdTree.set_particles_with_swap( m_instanceParticles );
        } else {
            for( std::size_t i = 0; i < m_instanceParticles.size(); ++i )
                m_instanceParticles[i].scale = scaleMinimum;
        }

        m_instanceIndex = m_instanceParticles.size();
    }

    void init_channel_map( const frantic::channels::channel_map& pcm ) {
        m_posAccessor = pcm.get_accessor<frantic::graphics::vector3f>( _T("Position") );
        m_delegateParticle.reset( new char[pcm.structure_size()] );
    }

  public:
    instancing_particle_istream( boost::shared_ptr<particle_istream> delegateStream,
                                 boost::shared_ptr<particle_istream> instanceStream, float scaleMinimum = 1.f,
                                 float scaleMultiplier = 0.f )
        : delegated_particle_istream( delegateStream )
        , m_instanceIndex( static_cast<std::size_t>( -1 ) ) {
        init_channel_map( delegateStream->get_channel_map() );
        load_instancing_particles( instanceStream, scaleMinimum, scaleMultiplier );
    }

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        init_channel_map( pcm );
        m_delegate->set_channel_map( pcm );
    }

    boost::int64_t particle_index() const {
        boost::int64_t delegateIndex = m_delegate->particle_index();
        if( delegateIndex < 0 )
            return delegateIndex;
        return delegateIndex * m_instanceParticles.size() + m_instanceIndex;
    }

    boost::int64_t particle_count() const { return m_delegate->particle_count() * m_instanceParticles.size(); }

    boost::int64_t particle_count_left() const {
        return ( m_delegate->particle_count_left() + 1 ) * m_instanceParticles.size() - m_instanceIndex;
    }

    bool get_particle( char* outBuffer ) {
        if( m_instanceIndex >= m_instanceParticles.size() ) {
            if( !m_delegate->get_particle( m_delegateParticle.get() ) )
                return false;
            m_instanceIndex = 0;
        }

        memcpy( outBuffer, m_delegateParticle.get(), m_delegate->particle_size() );
        m_posAccessor.get( outBuffer ) = m_instanceParticles[m_instanceIndex].scale * m_posAccessor.get( outBuffer ) +
                                         m_instanceParticles[m_instanceIndex].position;

        ++m_instanceIndex;
        return true;
    }

    bool get_particles( char* outBuffer, std::size_t& numParticles ) {
        std::size_t particleSize = m_delegate->particle_size();
        for( std::size_t i = 0; i < numParticles; ++i, outBuffer += particleSize ) {
            if( !get_particle( outBuffer ) ) {
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
