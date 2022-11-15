// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( E57_AVAILABLE )

#include <algorithm>
#include <iostream>

#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/e57_particle_istream.hpp>
#include <frantic/particles/streams/e57_particle_ostream.hpp>

#include <boost/algorithm/clamp.hpp>
#include <boost/scoped_array.hpp>

using namespace frantic;
using namespace frantic::particles::streams;
using boost::algorithm::clamp;
using channels::channel_map;
using channels::data_type_float16;
using channels::data_type_float32;
using channels::data_type_float64;
using particles::particle_array;

TEST( E57ParticleOstream, Write ) {
    const boost::int64_t expectedParticleCount = 4096;
    const tstring fileName = _T("TestOutputs/Write.e57");
    channel_map particleChannelMap;
    particleChannelMap.define_channel( _T("Position"), 3, data_type_float32 );
    particleChannelMap.define_channel( _T("Intensity"), 1, data_type_float32 );
    particleChannelMap.end_channel_definition();

    e57_particle_ostream os( fileName, particleChannelMap, expectedParticleCount, NULL, NULL );
    float rawParticleData[expectedParticleCount * 4];

    for( int i = 0; i < expectedParticleCount; ++i ) {
        rawParticleData[i * 4 + 0] = 0.1f * i;                                      // x
        rawParticleData[i * 4 + 1] = 0.1f * i;                                      // y
        rawParticleData[i * 4 + 2] = 0.1f * i;                                      // z
        rawParticleData[i * 4 + 3] = std::max( std::min( 0.02f * i, 1.0f ), 0.0f ); // intensity
        os.put_particle( (char*)( rawParticleData + i * 4 ) );
    }

    os.close();

    e57_particle_istream is( fileName );
    is.set_channel_map( particleChannelMap );
    ASSERT_EQ( expectedParticleCount, is.particle_count() );

    float writtenData[expectedParticleCount * 4];
    std::fill_n( writtenData, expectedParticleCount * 4, 0.0f );
    size_t numParticles = expectedParticleCount;

    is.get_particles( (char*)writtenData, numParticles );

    ASSERT_EQ( expectedParticleCount, numParticles );
    for( int i = 0; i < expectedParticleCount; ++i ) {
        EXPECT_EQ( rawParticleData[i * 4], writtenData[i * 4] );         // x
        EXPECT_EQ( rawParticleData[i * 4 + 1], writtenData[i * 4 + 1] ); // y
        EXPECT_EQ( rawParticleData[i * 4 + 2], writtenData[i * 4 + 2] ); // z
        EXPECT_EQ( rawParticleData[i * 4 + 3], writtenData[i * 4 + 3] ); // intensity
    }
}

TEST( E57ParticleOstream, IntenseWrite ) {
    const boost::int64_t expectedParticleCount = 1000000;
    const tstring fileName = _T("TestOutputs/IntenseWrite.e57");
    channel_map particleChannelMap;
    const size_t totalArity = 7;
    particleChannelMap.define_channel( _T("Position"), 3, data_type_float64 );
    particleChannelMap.define_channel( _T("Intensity"), 1, data_type_float64 );
    particleChannelMap.define_channel( _T("Color"), 3, data_type_float64 );
    particleChannelMap.end_channel_definition();

    e57_particle_ostream os( fileName, particleChannelMap, expectedParticleCount, NULL, NULL );
    boost::scoped_array<double> rawParticleData( new double[expectedParticleCount * totalArity] );

    for( int i = 0; i < expectedParticleCount; ++i ) {
        rawParticleData[i * totalArity + 0] = 0.000001 * i + 1.0;                              // x
        rawParticleData[i * totalArity + 1] = 0.000001 * i;                                    // y
        rawParticleData[i * totalArity + 2] = 10.0 - 0.000001 * i;                             // z
        rawParticleData[i * totalArity + 3] = clamp( 0.02 * i, 0.0, 1.0 );                     // intensity
        rawParticleData[i * totalArity + 4] = 0.1;                                             // r
        rawParticleData[i * totalArity + 5] = clamp( 0.001 * i, 0.0, 1.0 );                    // g
        rawParticleData[i * totalArity + 6] = clamp( fmod( 1.0 - 0.001 * i, 1.0 ), 0.0, 1.0 ); // b
        os.put_particle( (char*)( &rawParticleData[i * totalArity] ) );
    }

    os.close();

    e57_particle_istream is( fileName, channels::data_type_float64 );
    is.set_channel_map( particleChannelMap );
    ASSERT_EQ( expectedParticleCount, is.particle_count() );

    boost::scoped_array<double> writtenData( new double[expectedParticleCount * totalArity] );
    size_t numParticles = expectedParticleCount;

    is.get_particles( (char*)writtenData.get(), numParticles );

    ASSERT_EQ( expectedParticleCount, numParticles );
    for( int i = 0; i < expectedParticleCount; ++i ) {
        EXPECT_EQ( rawParticleData[i * totalArity], writtenData[i * totalArity] );         // x
        EXPECT_EQ( rawParticleData[i * totalArity + 1], writtenData[i * totalArity + 1] ); // y
        EXPECT_EQ( rawParticleData[i * totalArity + 2], writtenData[i * totalArity + 2] ); // z
        // istream is imprecise for intensity/r/g/b
        EXPECT_FLOAT_EQ( (float)rawParticleData[i * totalArity + 3],
                         (float)writtenData[i * totalArity + 3] ); // intensity
        EXPECT_FLOAT_EQ( (float)rawParticleData[i * totalArity + 4], (float)writtenData[i * totalArity + 4] ); // r
        EXPECT_FLOAT_EQ( (float)rawParticleData[i * totalArity + 5], (float)writtenData[i * totalArity + 5] ); // g
        EXPECT_FLOAT_EQ( (float)rawParticleData[i * totalArity + 6], (float)writtenData[i * totalArity + 6] ); // b
    }
}

#endif
