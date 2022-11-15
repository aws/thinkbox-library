// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/particle_array.hpp>
#include <frantic/volumetrics/implicitsurface/detail/xyzr_packet_array.hpp>

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif

using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::particles;
using namespace frantic::simd;

using frantic::volumetrics::implicitsurface::detail::xyzr_packet_array;

class XYZRPacketArrayLoad : public ::testing::TestWithParam<bool> {
  protected:
    void create_particles_to_load( frantic::particles::particle_array& outParticles, std::vector<char*>& outPParticles,
                                   std::size_t particleCount ) {
        const bool adjacentChannels = GetParam();

        channel_map channelMap;
        channelMap.define_channel( _T("Position"), 3, data_type_float32, 0 );
        if( adjacentChannels ) {
            channelMap.define_channel( _T("Radius"), 1, data_type_float32, 12 );
            channelMap.define_channel( _T("Color"), 3, data_type_float32, 16 );
        } else {
            channelMap.define_channel( _T("Color"), 3, data_type_float32, 12 );
            channelMap.define_channel( _T("Radius"), 1, data_type_float32, 24 );
        }
        channelMap.end_channel_definition( 4, true, false );

        outParticles.reset( channelMap );

        std::vector<char> buffer( channelMap.structure_size() );

        channel_accessor<vector3f> positionAcc = channelMap.get_accessor<vector3f>( _T("Position") );
        channel_accessor<float> radiusAcc = channelMap.get_accessor<float>( _T("Radius") );

        for( std::size_t i = 0; i < particleCount; ++i ) {
            positionAcc( buffer ) = vector3f( static_cast<float>( i + 1 ), 2, 3 );
            radiusAcc( buffer ) = 2 * static_cast<float>( i + 1 );
            outParticles.push_back( &buffer[0] );
        }

        outPParticles.clear();
        outPParticles.reserve( outParticles.size() );
        for( std::size_t i = 0; i < outParticles.size(); ++i ) {
            outPParticles.push_back( outParticles[i] );
        }
    }
};

TEST_P( XYZRPacketArrayLoad, Load ) {
    particle_array particles;
    std::vector<char*> pParticles;

    create_particles_to_load( particles, pParticles, 8 );

    const channel_map& channelMap = particles.get_channel_map();
    channel_accessor<vector3f> positionAcc = channelMap.get_accessor<vector3f>( _T("Position") );
    channel_accessor<float> radiusAcc = channelMap.get_accessor<float>( _T("Radius") );

    {
        xyzr_packet_array result;
        result.load( &pParticles[0], &pParticles[0], positionAcc, radiusAcc );

        EXPECT_EQ( 0, result.get_particle_count() );
        EXPECT_EQ( 0, result.get_filled_packet_count() );
        EXPECT_EQ( 0, result.get_remainder_particle_count() );
        EXPECT_FALSE( result.has_remainder() );
        EXPECT_EQ( 0, result.get_particles().size() );
    }

#ifdef FRANTIC_HAS_SSE2

    {
        xyzr_packet_array result;
        result.load( &pParticles[0], &pParticles[1], positionAcc, radiusAcc );

        EXPECT_EQ( 1, result.get_particle_count() );
        EXPECT_EQ( 0, result.get_filled_packet_count() );
        EXPECT_EQ( 1, result.get_remainder_particle_count() );
        EXPECT_TRUE( result.has_remainder() );
        EXPECT_EQ( 1, result.get_particles().size() );
        EXPECT_EQ( int_v( -1, 0, 0, 0 ), int_v::reinterpret( result.get_remainder_mask() ) );

        EXPECT_EQ( 1, result.get_particles()[0].position.x[0] );
        EXPECT_EQ( 2, result.get_particles()[0].radius[0] );
    }

    {
        xyzr_packet_array result;
        result.load( &pParticles[0], &pParticles[2], positionAcc, radiusAcc );

        EXPECT_EQ( 2, result.get_particle_count() );
        EXPECT_EQ( 0, result.get_filled_packet_count() );
        EXPECT_EQ( 2, result.get_remainder_particle_count() );
        EXPECT_TRUE( result.has_remainder() );
        EXPECT_EQ( 1, result.get_particles().size() );
        EXPECT_EQ( int_v( -1, -1, 0, 0 ), int_v::reinterpret( result.get_remainder_mask() ) );

        EXPECT_EQ( 1, result.get_particles()[0].position.x[0] );
        EXPECT_EQ( 2, result.get_particles()[0].radius[0] );
        EXPECT_EQ( 2, result.get_particles()[0].position.x[1] );
        EXPECT_EQ( 4, result.get_particles()[0].radius[1] );
    }

    {
        xyzr_packet_array result;
        result.load( &pParticles[0], &pParticles[3], positionAcc, radiusAcc );

        EXPECT_EQ( 3, result.get_particle_count() );
        EXPECT_EQ( 0, result.get_filled_packet_count() );
        EXPECT_EQ( 3, result.get_remainder_particle_count() );
        EXPECT_TRUE( result.has_remainder() );
        EXPECT_EQ( 1, result.get_particles().size() );
        EXPECT_EQ( int_v( -1, -1, -1, 0 ), int_v::reinterpret( result.get_remainder_mask() ) );

        EXPECT_EQ( 1, result.get_particles()[0].position.x[0] );
        EXPECT_EQ( 2, result.get_particles()[0].radius[0] );
        EXPECT_EQ( 2, result.get_particles()[0].position.x[1] );
        EXPECT_EQ( 4, result.get_particles()[0].radius[1] );
        EXPECT_EQ( 3, result.get_particles()[0].position.x[2] );
        EXPECT_EQ( 6, result.get_particles()[0].radius[2] );
    }

    {
        xyzr_packet_array result;
        result.load( &pParticles[0], &pParticles[4], positionAcc, radiusAcc );

        EXPECT_EQ( 4, result.get_particle_count() );
        EXPECT_EQ( 1, result.get_filled_packet_count() );
        EXPECT_EQ( 0, result.get_remainder_particle_count() );
        EXPECT_FALSE( result.has_remainder() );
        EXPECT_EQ( 1, result.get_particles().size() );

        EXPECT_EQ( float_v( 1, 2, 3, 4 ), result.get_particles()[0].position.x );
        EXPECT_EQ( float_v( 2 ), result.get_particles()[0].position.y );
        EXPECT_EQ( float_v( 3 ), result.get_particles()[0].position.z );
        EXPECT_EQ( float_v( 2, 4, 6, 8 ), result.get_particles()[0].radius );
    }

    {
        xyzr_packet_array result;
        result.load( &pParticles[0], &pParticles[5], positionAcc, radiusAcc );

        EXPECT_EQ( 5, result.get_particle_count() );
        EXPECT_EQ( 1, result.get_filled_packet_count() );
        EXPECT_EQ( 1, result.get_remainder_particle_count() );
        EXPECT_TRUE( result.has_remainder() );
        EXPECT_EQ( 2, result.get_particles().size() );
        EXPECT_EQ( int_v( -1, 0, 0, 0 ), int_v::reinterpret( result.get_remainder_mask() ) );

        EXPECT_EQ( float_v( 1, 2, 3, 4 ), result.get_particles()[0].position.x );
        EXPECT_EQ( float_v( 2 ), result.get_particles()[0].position.y[0] );
        EXPECT_EQ( float_v( 3 ), result.get_particles()[0].position.z[0] );
        EXPECT_EQ( float_v( 2, 4, 6, 8 ), result.get_particles()[0].radius );

        EXPECT_EQ( 5, result.get_particles()[1].position.x[0] );
        EXPECT_EQ( 10, result.get_particles()[1].radius[0] );
    }

#else

    {
        xyzr_packet_array result;
        result.load( &pParticles[0], &pParticles[1], positionAcc, radiusAcc );

        EXPECT_EQ( 1, result.get_particle_count() );
        EXPECT_EQ( 1, result.get_filled_packet_count() );
        EXPECT_EQ( 0, result.get_remainder_particle_count() );
        EXPECT_FALSE( result.has_remainder() );
        EXPECT_EQ( 1, result.get_particles().size() );
        EXPECT_EQ( int_v( 0 ), int_v::reinterpret( result.get_remainder_mask() ) );

        EXPECT_EQ( 1, result.get_particles()[0].position.x );
        EXPECT_EQ( 2, result.get_particles()[0].radius );
    }

#endif
}

INSTANTIATE_TEST_CASE_P( XYZRPacketArrayLoad, XYZRPacketArrayLoad, ::testing::Values( true, false ) );
