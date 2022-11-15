// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/implicitsurface/detail/xyzr_packet_array.hpp>

using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::simd;

using frantic::volumetrics::implicitsurface::detail::xyzr_packet_array;

xyzr_packet_array::xyzr_packet_array() { resize( 0 ); }

void xyzr_packet_array::clear() { resize( 0 ); }

void xyzr_packet_array::load( const char* const* particlesBegin, const char* const* particlesEnd,
                              const channel_accessor<vector3f>& positionAcc,
                              const channel_accessor<float>& radiusAcc ) {
    const std::size_t particleCount = particlesEnd - particlesBegin;

    resize( particleCount );

    if( particleCount == 0 ) {
        return;
    }

#ifdef FRANTIC_HAS_SSE2

    std::ptrdiff_t offset = reinterpret_cast<const char*>( &radiusAcc( particlesBegin[0] ) ) -
                            reinterpret_cast<const char*>( &positionAcc( particlesBegin[0] ).x );
    if( offset == 12 ) {
        load_adjacent_internal( particlesBegin, positionAcc );
    } else {
        std::size_t i = 0;
        std::size_t outIndex = 0;

        for( ; i < get_filled_particle_count(); i += float_v::static_size ) {
            float_v x = float_v::load3( &positionAcc( particlesBegin[i] ).x );
            float r0 = radiusAcc( particlesBegin[i] );
            float_v y = float_v::load3( &positionAcc( particlesBegin[i + 1] ).x );
            float r1 = radiusAcc( particlesBegin[i + 1] );
            float_v z = float_v::load3( &positionAcc( particlesBegin[i + 2] ).x );
            float r2 = radiusAcc( particlesBegin[i + 2] );
            float_v w = float_v::load3( &positionAcc( particlesBegin[i + 3] ).x );
            float r3 = radiusAcc( particlesBegin[i + 3] );
            _MM_TRANSPOSE4_PS( x.native(), y.native(), z.native(), w.native() );
            float_v r( r0, r1, r2, r3 );
            m_particles[outIndex].position.set( x, y, z );
            m_particles[outIndex].radius = r;

            ++outIndex;
        }

        if( m_remainderParticleCount > 0 ) {
            float_v x, y, z, dummy( 0 );
            float r0, r1, r2;
            x = float_v::load3( &positionAcc( particlesBegin[i] ).x );
            r0 = radiusAcc( particlesBegin[i] );
            if( m_remainderParticleCount > 1 ) {
                y = float_v::load3( &positionAcc( particlesBegin[i + 1] ).x );
                r1 = radiusAcc( particlesBegin[i + 1] );
            } else {
                y = 0;
                r1 = 0;
            }
            if( m_remainderParticleCount > 2 ) {
                z = float_v::load3( &positionAcc( particlesBegin[i + 2] ).x );
                r2 = radiusAcc( particlesBegin[i + 2] );
            } else {
                z = 0;
                r2 = 0;
            }
            _MM_TRANSPOSE4_PS( x.native(), y.native(), z.native(), dummy.native() );
            float_v r( r0, r1, r2, 0 );

            m_particles[outIndex].position.set( x, y, z );
            m_particles[outIndex].radius = r;
        }
    }

#else

    BOOST_STATIC_ASSERT( float_v::static_size == 1 );

    for( std::size_t i = 0; i < get_filled_particle_count(); ++i ) {
        const frantic::graphics::vector3f position = positionAcc( particlesBegin[i] );
        m_particles[i].position.set( position.x, position.y, position.z );
        m_particles[i].radius = radiusAcc( particlesBegin[i] );
    }

#endif
}

namespace {

float_v create_remainder_mask( std::size_t remainderCount ) {
    switch( remainderCount ) {
    case 0:
        return float_v( 0 );
#ifdef FRANTIC_HAS_SSE2
    case 1:
        return float_v::reinterpret( int_v( -1, 0, 0, 0 ) );
    case 2:
        return float_v::reinterpret( int_v( -1, -1, 0, 0 ) );
    case 3:
        return float_v::reinterpret( int_v( -1, -1, -1, 0 ) );
#endif
    default:
        throw std::runtime_error( "create_remainder_mask Error: unexpected remainder count" );
    }
}

} // anonymous namespace

void xyzr_packet_array::resize( std::size_t particleCount ) {
    m_particleCount = particleCount;
    const std::size_t packetCount = ( particleCount + float_v::static_size - 1 ) / float_v::static_size;
    m_filledPacketCount = particleCount / float_v::static_size;
    m_filledParticleCount = m_filledPacketCount * float_v::static_size;
    m_remainderParticleCount = particleCount - ( m_filledPacketCount * float_v::static_size );
    m_remainderParticleMask = create_remainder_mask( m_remainderParticleCount );
    m_particles.resize( packetCount );
}

#ifdef FRANTIC_HAS_SSE2

void xyzr_packet_array::load_adjacent_internal( const char* const* particlesBegin,
                                                const channel_accessor<vector3f>& positionAcc ) {
    std::size_t outIndex = 0;
    std::size_t i = 0;

    for( ; i < get_filled_particle_count(); i += float_v::static_size ) {
        float_v x = float_v::load( &positionAcc( particlesBegin[i] ).x );
        float_v y = float_v::load( &positionAcc( particlesBegin[i + 1] ).x );
        float_v z = float_v::load( &positionAcc( particlesBegin[i + 2] ).x );
        float_v r = float_v::load( &positionAcc( particlesBegin[i + 3] ).x );
        _MM_TRANSPOSE4_PS( x.native(), y.native(), z.native(), r.native() );
        m_particles[outIndex].position.set( x, y, z );
        m_particles[outIndex].radius = r;

        ++outIndex;
    }

    if( m_remainderParticleCount > 0 ) {
        float_v x, y, z, r( 0 );
        x = float_v::load( &positionAcc( particlesBegin[i] ).x );
        if( m_remainderParticleCount > 1 ) {
            y = float_v::load( &positionAcc( particlesBegin[i + 1] ).x );
        } else {
            y = 0;
        }
        if( m_remainderParticleCount > 2 ) {
            z = float_v::load( &positionAcc( particlesBegin[i + 2] ).x );
        } else {
            z = 0;
        }
        _MM_TRANSPOSE4_PS( x.native(), y.native(), z.native(), r.native() );
        m_particles[outIndex].position.set( x, y, z );
        m_particles[outIndex].radius = r;
    }
}

#endif
