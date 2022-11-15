// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <vector>

#include <boost/lexical_cast.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/channel_map_lerp.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/particles/particle_classes.hpp>
#include <frantic/particles/particle_kdtree.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>

#include <frantic/diagnostics/profiling_section.hpp>

#include <frantic/locale/locale.hpp>

using namespace std;
using namespace frantic::graphics;
using namespace frantic::particles;
using namespace frantic::channels;

TEST( ChannelMap, Creation ) {

    channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<vector3f>( _T("Velocity") );
    channelMap.define_channel<vector3f>( _T("Acceleration") );
    channelMap.define_channel<half>( _T("U") );
    channelMap.define_channel<half>( _T("V") );
    channelMap.define_channel<boost::uint32_t>( _T("ID") );

    EXPECT_THROW( channelMap.define_channel<vector3f>( _T("Position") ), std::runtime_error )
        << "Defining the Position channel twice";

    EXPECT_THROW( channel_accessor<vector3f> temp = channelMap.get_accessor<vector3f>( _T("Position") ),
                  std::runtime_error )
        << " Trying to get an accessor before calling end_channel_definition ";

    EXPECT_THROW( channel_cvt_accessor<vector3f> temp = channelMap.get_cvt_accessor<vector3f>( _T("Position") ),
                  std::runtime_error )
        << "Trying to get a cvt accessor before calling end_channel_definition ";

    channelMap.end_channel_definition();

    EXPECT_THROW( channelMap.define_channel<vector3f>( _T("Exceptional") );, std::runtime_error )
        << "Trying to access a non-existent channel";

    // Check that the offsets we got accord with our expectations.  The channels have been sorted in descending order of
    // type size, using a stable sort.
    EXPECT_EQ( channelMap.channel_offset( _T("Position") ), 0 );
    EXPECT_EQ( channelMap.channel_offset( _T("Velocity") ), 12 );
    EXPECT_EQ( channelMap.channel_offset( _T("Acceleration") ), 24 );
    EXPECT_EQ( channelMap.channel_offset( _T("ID") ), 36 );
    EXPECT_EQ( channelMap.channel_offset( _T("U") ), 40 );
    EXPECT_EQ( channelMap.channel_offset( _T("V") ), 42 );

    // Make sure that the total byte size is right
    EXPECT_EQ( channelMap.structure_size(), 44 );

    // Make sure the channel count is right
    EXPECT_EQ( channelMap.channel_count(), 6 );

    // When accessed by index, the channels should be in the order of their defined position
    EXPECT_TRUE( channelMap[0].name() == _T("Position") );
    EXPECT_TRUE( channelMap[1].name() == _T("Velocity") );
    EXPECT_TRUE( channelMap[2].name() == _T("Acceleration") );
    EXPECT_TRUE( channelMap[3].name() == _T("ID") );
    EXPECT_TRUE( channelMap[4].name() == _T("U") );
    EXPECT_TRUE( channelMap[5].name() == _T("V") );
}

TEST( ChannelMap, Accessor ) {

    struct {
        double density;
        vector3f position;
        half color[3];
    } testParticle;

    channel_map pcm;
    pcm.define_channel<double>( _T("Density") );
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel( _T("Color"), 3, data_type_float16 );
    pcm.end_channel_definition( 1, true );

    testParticle.position = vector3f( 1, 2, 3 );
    testParticle.color[0] = 4;
    testParticle.color[1] = 5;
    testParticle.color[2] = 6;
    testParticle.density = 7;

    char* testParticlePtr = reinterpret_cast<char*>( &testParticle );

    // Test that the raw particle_accessor class works
    channel_accessor<vector3f> pav = pcm.get_accessor<vector3f>( _T("Position") );
    EXPECT_TRUE( pav.get( testParticlePtr ) == testParticle.position );

    EXPECT_THROW( pav = pcm.get_accessor<vector3f>( _T("Color") ), std::exception )
        << "Creating an accessor with a non-matching data type ";

    channel_accessor<double> pad = pcm.get_accessor<double>( _T("Density") );
    EXPECT_EQ( pad.get( testParticlePtr ), testParticle.density );

    // Test that conversion accessor works
    channel_cvt_accessor<vector3f> pcav = pcm.get_cvt_accessor<vector3f>( _T("Position") );
    EXPECT_TRUE( pcav.get( testParticlePtr ) == testParticle.position );
    pcav = pcm.get_cvt_accessor<vector3f>( _T("Color") );
    EXPECT_TRUE( pcav.get( testParticlePtr ) == vector3f( 4, 5, 6 ) );

    channel_cvt_accessor<double> pcad = pcm.get_cvt_accessor<double>( _T("Density") );
    EXPECT_EQ( pcad.get( testParticlePtr ), testParticle.density );

    channel_cvt_accessor<float> pcaf = pcm.get_cvt_accessor<float>( _T("Density") );
    EXPECT_EQ( pcaf.get( testParticlePtr ), 7.f );

    EXPECT_THROW( pcaf = pcm.get_cvt_accessor<float>( _T("Color") ), std::exception )
        << "Creating a cvt accessor with a non-matching arity ";
}

TEST( ChannelMap, Adaptor ) {

    channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<vector3f>( _T("Velocity") );
    channelMap.define_channel<vector3f>( _T("Acceleration") );
    channelMap.end_channel_definition();

    // Test that the identity transform produces a single memcpy block
    channel_map_adaptor pcmaIdentity( channelMap, channelMap );
    EXPECT_EQ( pcmaIdentity.block_operation_count(), 1 );
    EXPECT_TRUE( pcmaIdentity.is_identity() );

    // Test that swapping two channels works
    channel_map channelMapCopyTest;
    channelMapCopyTest.define_channel( _T("Acceleration"), 3, data_type_float64 );
    channelMapCopyTest.define_channel( _T("Velocity"), 3, data_type_float32 );
    channelMapCopyTest.define_channel( _T("Position"), 3, data_type_float64 );
    channelMapCopyTest.end_channel_definition();

    channel_map_adaptor pcmaCopyForward( channelMapCopyTest, channelMap );
    channel_map_adaptor pcmaCopyBackward( channelMap, channelMapCopyTest );
    vector<char> inputParticle( channelMap.structure_size() ), outputParticle( channelMap.structure_size() ),
        copyTestParticle( channelMapCopyTest.structure_size() );
    channel_accessor<vector3f> pos = channelMap.get_accessor<vector3f>( _T("Position") );
    channel_accessor<vector3f> vel = channelMap.get_accessor<vector3f>( _T("Velocity") );
    channel_accessor<vector3f> acc = channelMap.get_accessor<vector3f>( _T("Acceleration") );
    pos.get( &inputParticle[0] ) = vector3f( 3, 7, 1.457f );
    vel.get( &inputParticle[0] ) = vector3f( 12, 5, 1.357f );
    acc.get( &inputParticle[0] ) = vector3f( 2.99f, 1e38f, 18.451f );
    pcmaCopyForward.copy_structure( &copyTestParticle[0], &inputParticle[0] );
    pcmaCopyBackward.copy_structure( &outputParticle[0], &copyTestParticle[0] );
    EXPECT_EQ( memcmp( &inputParticle[0], &outputParticle[0], channelMap.structure_size() ), 0 );

    // Test that converting the data types works
    channel_map channelMapConvert1;
    channelMapConvert1.define_channel( _T("Position"), 3, data_type_float64 );
    channelMapConvert1.define_channel( _T("Velocity"), 3, data_type_float16 );
    channelMapConvert1.define_channel( _T("Acceleration"), 3, data_type_float32 );
    channelMapConvert1.end_channel_definition();

    channel_map_adaptor pcmaConvert1( channelMapConvert1, channelMap );
    EXPECT_EQ( pcmaConvert1.block_operation_count(), 3 );

    // Test that it fails when there's an invalid data_type conversion

    channel_map channelMapConvert2;
    channelMapConvert2.define_channel( _T("Position"), 3, data_type_float16 );
    channelMapConvert2.define_channel( _T("Velocity"), 3, data_type_float16 );
    channelMapConvert2.define_channel( _T("Acceleration"), 3, data_type_float16 );
    channelMapConvert2.end_channel_definition();

    channel_map_adaptor pcmaConvert2( channelMapConvert2, channelMap );
    EXPECT_EQ( pcmaConvert2.block_operation_count(), 1 );

    // Test that it fails when there's an invalid data_type conversion
    channel_map channelMapConvert3;
    channelMapConvert3.define_channel( _T("Position"), 3, data_type_int8 );
    channelMapConvert3.define_channel( _T("Velocity"), 3, data_type_float16 );
    channelMapConvert3.define_channel( _T("Acceleration"), 3, data_type_float16 );
    channelMapConvert3.end_channel_definition();

    EXPECT_THROW( channel_map_adaptor pcmaConvert3( channelMapConvert3, channelMap ), std::exception )
        << "Creating a channel_map_adaptor with an arity change";

    // Test that it fails when there's an arity change
    channel_map channelMapConvert4;
    channelMapConvert4.define_channel( _T("Position"), 3, data_type_float32 );
    channelMapConvert4.define_channel( _T("Velocity"), 3, data_type_float16 );
    channelMapConvert4.define_channel( _T("Acceleration"), 4, data_type_float16 );
    channelMapConvert4.end_channel_definition();

    EXPECT_THROW( channel_map_adaptor pcmaConvert4( channelMapConvert4, channelMap ), std::exception );

    // Test that it throws when you try to use an incomplete channel map
    channel_map channelMapIncomplete;
    channelMapIncomplete.define_channel( _T("Position"), 3, data_type_float32 );

    EXPECT_THROW( channel_map_adaptor pcmaIncomplete( channelMapIncomplete, channelMap ), std::exception )
        << "Creating a channel_map_adaptor from an incomplete ";

    // Create an adaptor which uses all valid type conversions

    channel_map pcmAllTypesSource;
    pcmAllTypesSource.define_channel( _T("F16_F32"), 1, data_type_float16 );
    pcmAllTypesSource.define_channel( _T("F16_F64"), 1, data_type_float16 );
    pcmAllTypesSource.define_channel( _T("F32_F16"), 1, data_type_float32 );
    pcmAllTypesSource.define_channel( _T("F32_F64"), 1, data_type_float32 );
    pcmAllTypesSource.define_channel( _T("F64_F16"), 1, data_type_float64 );
    pcmAllTypesSource.define_channel( _T("F64_F32"), 1, data_type_float64 );
    pcmAllTypesSource.define_channel( _T("I8_I16"), 1, data_type_int8 );
    pcmAllTypesSource.define_channel( _T("I8_I32"), 1, data_type_int8 );
    pcmAllTypesSource.define_channel( _T("I8_I64"), 1, data_type_int8 );
    pcmAllTypesSource.define_channel( _T("I16_I8"), 1, data_type_int16 );
    pcmAllTypesSource.define_channel( _T("I16_I32"), 1, data_type_int16 );
    pcmAllTypesSource.define_channel( _T("I16_I64"), 1, data_type_int16 );
    pcmAllTypesSource.define_channel( _T("I32_I8"), 1, data_type_int32 );
    pcmAllTypesSource.define_channel( _T("I32_I16"), 1, data_type_int32 );
    pcmAllTypesSource.define_channel( _T("I32_I64"), 1, data_type_int32 );
    pcmAllTypesSource.define_channel( _T("I64_I8"), 1, data_type_int64 );
    pcmAllTypesSource.define_channel( _T("I64_I16"), 1, data_type_int64 );
    pcmAllTypesSource.define_channel( _T("I64_I32"), 1, data_type_int64 );

    pcmAllTypesSource.define_channel( _T("U8_U16"), 1, data_type_uint8 );
    pcmAllTypesSource.define_channel( _T("U8_U32"), 1, data_type_uint8 );
    pcmAllTypesSource.define_channel( _T("U8_U64"), 1, data_type_uint8 );
    pcmAllTypesSource.define_channel( _T("U16_U8"), 1, data_type_uint16 );
    pcmAllTypesSource.define_channel( _T("U16_U32"), 1, data_type_uint16 );
    pcmAllTypesSource.define_channel( _T("U16_U64"), 1, data_type_uint16 );
    pcmAllTypesSource.define_channel( _T("U32_U8"), 1, data_type_uint32 );
    pcmAllTypesSource.define_channel( _T("U32_U16"), 1, data_type_uint32 );
    pcmAllTypesSource.define_channel( _T("U32_U64"), 1, data_type_uint32 );
    pcmAllTypesSource.define_channel( _T("U64_U8"), 1, data_type_uint64 );
    pcmAllTypesSource.define_channel( _T("U64_U16"), 1, data_type_uint64 );
    pcmAllTypesSource.define_channel( _T("U64_U32"), 1, data_type_uint64 );

    pcmAllTypesSource.define_channel( _T("U8_I16"), 1, data_type_uint8 );
    pcmAllTypesSource.define_channel( _T("U8_I32"), 1, data_type_uint8 );
    pcmAllTypesSource.define_channel( _T("U8_I64"), 1, data_type_uint8 );
    pcmAllTypesSource.define_channel( _T("U16_I32"), 1, data_type_uint16 );
    pcmAllTypesSource.define_channel( _T("U16_I64"), 1, data_type_uint16 );
    pcmAllTypesSource.define_channel( _T("U32_I64"), 1, data_type_uint32 );

    pcmAllTypesSource.end_channel_definition();

    channel_map pcmAllTypesDest;
    pcmAllTypesDest.define_channel( _T("F16_F32"), 1, data_type_float32 );
    pcmAllTypesDest.define_channel( _T("F16_F64"), 1, data_type_float64 );
    pcmAllTypesDest.define_channel( _T("F32_F16"), 1, data_type_float16 );
    pcmAllTypesDest.define_channel( _T("F32_F64"), 1, data_type_float64 );
    pcmAllTypesDest.define_channel( _T("F64_F16"), 1, data_type_float16 );
    pcmAllTypesDest.define_channel( _T("F64_F32"), 1, data_type_float32 );

    pcmAllTypesDest.define_channel( _T("I8_I16"), 1, data_type_int16 );
    pcmAllTypesDest.define_channel( _T("I8_I32"), 1, data_type_int32 );
    pcmAllTypesDest.define_channel( _T("I8_I64"), 1, data_type_int64 );
    pcmAllTypesDest.define_channel( _T("I16_I8"), 1, data_type_int8 );
    pcmAllTypesDest.define_channel( _T("I16_I32"), 1, data_type_int32 );
    pcmAllTypesDest.define_channel( _T("I16_I64"), 1, data_type_int64 );
    pcmAllTypesDest.define_channel( _T("I32_I8"), 1, data_type_int8 );
    pcmAllTypesDest.define_channel( _T("I32_I16"), 1, data_type_int16 );
    pcmAllTypesDest.define_channel( _T("I32_I64"), 1, data_type_int64 );
    pcmAllTypesDest.define_channel( _T("I64_I8"), 1, data_type_int8 );
    pcmAllTypesDest.define_channel( _T("I64_I16"), 1, data_type_int16 );
    pcmAllTypesDest.define_channel( _T("I64_I32"), 1, data_type_int32 );

    pcmAllTypesDest.define_channel( _T("U8_U16"), 1, data_type_uint16 );
    pcmAllTypesDest.define_channel( _T("U8_U32"), 1, data_type_uint32 );
    pcmAllTypesDest.define_channel( _T("U8_U64"), 1, data_type_uint64 );
    pcmAllTypesDest.define_channel( _T("U16_U8"), 1, data_type_uint8 );
    pcmAllTypesDest.define_channel( _T("U16_U32"), 1, data_type_uint32 );
    pcmAllTypesDest.define_channel( _T("U16_U64"), 1, data_type_uint64 );
    pcmAllTypesDest.define_channel( _T("U32_U8"), 1, data_type_uint8 );
    pcmAllTypesDest.define_channel( _T("U32_U16"), 1, data_type_uint16 );
    pcmAllTypesDest.define_channel( _T("U32_U64"), 1, data_type_uint64 );
    pcmAllTypesDest.define_channel( _T("U64_U8"), 1, data_type_uint8 );
    pcmAllTypesDest.define_channel( _T("U64_U16"), 1, data_type_uint16 );
    pcmAllTypesDest.define_channel( _T("U64_U32"), 1, data_type_uint32 );

    pcmAllTypesDest.define_channel( _T("U8_I16"), 1, data_type_int16 );
    pcmAllTypesDest.define_channel( _T("U8_I32"), 1, data_type_int32 );
    pcmAllTypesDest.define_channel( _T("U8_I64"), 1, data_type_int64 );
    pcmAllTypesDest.define_channel( _T("U16_I32"), 1, data_type_int32 );
    pcmAllTypesDest.define_channel( _T("U16_I64"), 1, data_type_int64 );
    pcmAllTypesDest.define_channel( _T("U32_I64"), 1, data_type_int64 );

    pcmAllTypesDest.end_channel_definition();

    // Verify that both channel_map instances got the same number of channels
    EXPECT_EQ( pcmAllTypesSource.channel_count(), pcmAllTypesDest.channel_count() );

    EXPECT_NO_THROW( {
        channel_map_adaptor pcmaAllTypesConvert( pcmAllTypesDest, pcmAllTypesSource );
        // The number of block operations should match the number of channels
        EXPECT_EQ( pcmaAllTypesConvert.block_operation_count(), pcmAllTypesSource.channel_count() );
    } );
}

TEST( ChannelMap, Lerp ) {
    channel_map channelMap;
    channelMap.define_channel<half>( _T("Half") );
    channelMap.define_channel<float>( _T("Float") );
    channelMap.define_channel<double>( _T("Double") );
    channelMap.define_channel<int>( _T("Int") );
    channelMap.end_channel_definition();

    channel_map_lerp pcmLerp( channelMap );

    channel_accessor<half> accHalf = channelMap.get_accessor<half>( _T("Half") );
    channel_accessor<float> accFloat = channelMap.get_accessor<float>( _T("Float") );
    channel_accessor<double> accDouble = channelMap.get_accessor<double>( _T("Double") );
    channel_accessor<int> accInt = channelMap.get_accessor<int>( _T("Int") );

    std::vector<char> particleA( channelMap.structure_size() ), particleB( channelMap.structure_size() ),
        particleC( channelMap.structure_size() );

    accHalf.get( &particleA[0] ) = 1;
    accHalf.get( &particleB[0] ) = 5;
    accFloat.get( &particleA[0] ) = -1;
    accFloat.get( &particleB[0] ) = 3;
    accDouble.get( &particleA[0] ) = 2;
    accDouble.get( &particleB[0] ) = -2;
    accInt.get( &particleA[0] ) = 1;
    accInt.get( &particleB[0] ) = 5;

    // Test the lerp function where the inputs are different from the output

    // Test the extreme cases
    pcmLerp.lerp( &particleC[0], &particleA[0], &particleB[0], 0 );
    EXPECT_TRUE( memcmp( &particleC[0], &particleA[0], channelMap.structure_size() ) == 0 );
    pcmLerp.lerp( &particleC[0], &particleA[0], &particleB[0], 1 );
    EXPECT_TRUE( memcmp( &particleC[0], &particleB[0], channelMap.structure_size() ) == 0 );

    pcmLerp.lerp( &particleC[0], &particleA[0], &particleB[0], 0.5f );
    EXPECT_FLOAT_EQ( accHalf.get( &particleC[0] ), 3.f );
    EXPECT_FLOAT_EQ( accFloat.get( &particleC[0] ), 1.f );
    EXPECT_DOUBLE_EQ( accDouble.get( &particleC[0] ), 0.f );
    // The integer parameter could go either way on this, maybe we should define what 0.5 does?

    pcmLerp.lerp( &particleC[0], &particleA[0], &particleB[0], 0.25f );
    EXPECT_FLOAT_EQ( accHalf.get( &particleC[0] ), 2.f );
    EXPECT_FLOAT_EQ( accFloat.get( &particleC[0] ), 0.f );
    EXPECT_DOUBLE_EQ( accDouble.get( &particleC[0] ), 1.f );
    EXPECT_EQ( accInt.get( &particleC[0] ), 1 );

    pcmLerp.lerp( &particleC[0], &particleA[0], &particleB[0], 0.75f );
    EXPECT_FLOAT_EQ( accHalf.get( &particleC[0] ), 4.f );
    EXPECT_FLOAT_EQ( accFloat.get( &particleC[0] ), 2.f );
    EXPECT_DOUBLE_EQ( accDouble.get( &particleC[0] ), -1.f );
    EXPECT_EQ( accInt.get( &particleC[0] ), 5 );

    // Test the lerp function where the first input is also the output

    // Test the extreme cases
    particleC = particleA;
    pcmLerp.lerp( &particleC[0], &particleB[0], 0 );
    EXPECT_TRUE( memcmp( &particleC[0], &particleA[0], channelMap.structure_size() ) == 0 );
    particleC = particleA;
    pcmLerp.lerp( &particleC[0], &particleB[0], 1 );
    EXPECT_TRUE( memcmp( &particleC[0], &particleB[0], channelMap.structure_size() ) == 0 );

    particleC = particleA;
    pcmLerp.lerp( &particleC[0], &particleB[0], 0.5f );
    EXPECT_FLOAT_EQ( accHalf.get( &particleC[0] ), 3.f );
    EXPECT_FLOAT_EQ( accFloat.get( &particleC[0] ), 1.f );
    EXPECT_DOUBLE_EQ( accDouble.get( &particleC[0] ), 0. );
    // The integer parameter could go either way on this, maybe we should define what 0.5 does?

    particleC = particleA;
    pcmLerp.lerp( &particleC[0], &particleB[0], 0.25f );
    EXPECT_FLOAT_EQ( accHalf.get( &particleC[0] ), 2.f );
    EXPECT_FLOAT_EQ( accFloat.get( &particleC[0] ), 0.f );
    EXPECT_DOUBLE_EQ( accDouble.get( &particleC[0] ), 1.0 );
    EXPECT_EQ( accInt.get( &particleC[0] ), 1 );

    particleC = particleA;
    pcmLerp.lerp( &particleC[0], &particleB[0], 0.75f );
    EXPECT_FLOAT_EQ( accHalf.get( &particleC[0] ), 4.f );
    EXPECT_FLOAT_EQ( accFloat.get( &particleC[0] ), 2.f );
    EXPECT_DOUBLE_EQ( accDouble.get( &particleC[0] ), -1.0 );
    EXPECT_EQ( accInt.get( &particleC[0] ), 5 );
}

TEST( ChannelMap, String ) {

    // A layout with some numbers and a string
    channel_map channelMapA;
    channelMapA.define_channel<half>( _T("Half") );
    channelMapA.define_channel<float>( _T("Float") );
    channelMapA.define_channel<double>( _T("Double") );
    channelMapA.define_channel<int>( _T("Int") );
    channelMapA.define_channel<frantic::tstring>( _T("String") );
    channelMapA.end_channel_definition();

    // A different layout with the same string channel
    channel_map channelMapC;
    channelMapC.define_channel<frantic::tstring>( _T("String") );
    channelMapC.end_channel_definition();

    // A layout with a different string channel
    channel_map channelMapD;
    channelMapD.define_channel<frantic::tstring>( _T("DifferentString") );
    channelMapD.end_channel_definition();

    std::vector<char> particleA( channelMapA.structure_size() ), particleB( channelMapA.structure_size() ),
        particleC( channelMapC.structure_size() ), particleD( channelMapD.structure_size() );

    channelMapA.construct_structure( particleA );
    channelMapA.construct_structure( particleB );
    channelMapC.construct_structure( particleC );
    channelMapD.construct_structure( particleD );

    channel_accessor<frantic::tstring> strA = channelMapA.get_accessor<frantic::tstring>( _T("String") );
    channel_cvt_accessor<frantic::tstring> strCvtA = channelMapA.get_cvt_accessor<frantic::tstring>( _T("String") );
    channel_accessor<frantic::tstring> strC = channelMapC.get_accessor<frantic::tstring>( _T("String") );
    channel_accessor<frantic::tstring> strD = channelMapD.get_accessor<frantic::tstring>( _T("DifferentString") );

    // Check that it defaults to an empty string, and various ways to do string comparisons compile ok
    EXPECT_TRUE( strA( particleA ) == _T("") );
    EXPECT_TRUE( strCvtA( particleA ) == _T("") );
    EXPECT_TRUE( _T("") == strA( particleA ) );
    EXPECT_TRUE( _T("") == strCvtA( particleA ) );
    EXPECT_TRUE( strA( particleA ) == frantic::tstring() );
    EXPECT_TRUE( strCvtA( particleA ) == frantic::tstring() );
    EXPECT_TRUE( frantic::tstring() == strA( particleA ) );
    EXPECT_TRUE( frantic::tstring() == strCvtA( particleA ) );

    // Check that assignment works
    strA( particleA ) = _T("testing");
    EXPECT_TRUE( strA( particleA ) == _T("testing") );
    // Check that != operator works too
    EXPECT_TRUE( strA( particleA ) != _T("") );

    channelMapA.copy_structure( particleB, particleA );
    EXPECT_TRUE( strA( particleB ) == _T("testing") );
    strA( particleA ) = _T("more"); // This is to make sure there wasn't just a memcpy
    EXPECT_TRUE( strA( particleB ) == _T("testing") );

    // Assign the string a bunch of times, to check for a memory look. (i.e. run the test and look at task manager)
    // for( int i = 0; i < 10000000; ++i ) {
    //	strA(particleA) = "teshjmdasfa;lskjdflasjflkjsadlfkjaslkdfjslakjdflaskjdflaskjdflsakjdflaskjdflkt" +
    // boost::lexical_cast<std::string>(i);
    //}

    // Check that the channel_map_adaptor will work to copy a channel
    strA( particleA ) = _T("again");
    channel_map_adaptor cma( channelMapC, channelMapA );
    cma.copy_structure( particleC, particleA );
    strA( particleA ) = _T("editA"); // This is to make sure there wasn't just a memcpy
    EXPECT_TRUE( strC( particleC ) == _T("again") );

    // Make sure that it doesn't touch the DifferentString channel in particle D
    cma.set( channelMapD, channelMapA );
    strD( particleD ) = _T("d test");
    cma.copy_structure( particleD, particleA );
    strA( particleA ) = _T("a");
    EXPECT_TRUE( strD( particleD ) == _T("d test") );

    // Now check that copying a string channel from a default particle works ok
    cma.set( channelMapA, channelMapD );
    // cma.dump(cerr);
    cerr << endl;
    strA( particleB ) = _T("default");
    cma.copy_structure( particleA, particleD, particleB );
    strA( particleB ) = _T("no longer default");
    EXPECT_TRUE( strA( particleA ) == _T("default") );

    channelMapA.destruct_structure( particleA );
    channelMapA.destruct_structure( particleB );
    channelMapC.destruct_structure( particleC );
    channelMapD.destruct_structure( particleD );
}

TEST( ChannelMap, PropertyMap ) {

    channel_map channelMapA;
    channelMapA.define_channel<half>( _T("Half") );
    channelMapA.define_channel<float>( _T("Float") );
    channelMapA.define_channel<double>( _T("Double") );
    channelMapA.define_channel<int>( _T("Int") );
    channelMapA.define_channel<frantic::tstring>( _T("String") );
    channelMapA.end_channel_definition();

    channel_map channelMapB;
    channelMapB.define_channel<float>( _T("Float") );
    channelMapB.define_channel<int>( _T("Int") );
    channelMapB.end_channel_definition();

    channel_map channelMapC;
    channelMapC.define_channel<float>( _T("MoreFloat") );
    channelMapC.define_channel<int>( _T("MoreInt") );
    channelMapC.define_channel<frantic::tstring>( _T("Str") );
    channelMapC.end_channel_definition();

    property_map pm( channelMapA );

    // Check that setting the values works ok
    pm.get<float>( _T("Float") ) = 3;
    pm.get<double>( _T("Double") ) = 7.5f;
    pm.get<half>( _T("Half") ) = 6;
    pm.get<int>( _T("Int") ) = 55;
    pm.get<frantic::tstring>( _T("String") ) = _T("testing");
    EXPECT_TRUE( pm.get<float>( _T("Float") ) == 3 );
    EXPECT_TRUE( pm.get<double>( _T("Double") ) == 7.5f );
    EXPECT_TRUE( pm.get<half>( _T("Half") ) == 6 );
    EXPECT_TRUE( pm.get<int>( _T("Int") ) == 55 );
    EXPECT_TRUE( pm.get<frantic::tstring>( _T("String") ) == _T("testing") );
    // Check has_property with a non-existent property
    EXPECT_TRUE( !pm.has_property( _T("strinG") ) );

    // Verify that assignment works
    property_map pm2;
    pm2 = pm;
    EXPECT_TRUE( pm2.get<float>( _T("Float") ) == 3 );
    EXPECT_TRUE( pm2.get<double>( _T("Double") ) == 7.5f );
    EXPECT_TRUE( pm2.get<half>( _T("Half") ) == 6 );
    EXPECT_TRUE( pm2.get<int>( _T("Int") ) == 55 );
    EXPECT_TRUE( pm2.get<frantic::tstring>( _T("String") ) == _T("testing") );

    // Verify that data gets retained ok when setting the channel_map
    pm2.set_channel_map( channelMapB );
    EXPECT_TRUE( pm2.get<float>( _T("Float") ) == 3 );
    EXPECT_TRUE( !pm2.has_property( _T("Double") ) );
    EXPECT_TRUE( !pm2.has_property( _T("Half") ) );
    EXPECT_TRUE( pm2.get<int>( _T("Int") ) == 55 );
    EXPECT_TRUE( !pm2.has_property( _T("String") ) );

    // Check that swapping switches things around
    pm.swap( pm2 );
    EXPECT_TRUE( !pm.has_property( _T("Double") ) );
    EXPECT_TRUE( pm2.get<double>( _T("Double") ) == 7.5f );

    // Check that clearing works
    pm2.clear();
    EXPECT_TRUE( pm2.empty() );
    EXPECT_TRUE( !pm2.has_property( _T("Float") ) );

    // Check that merging two property maps works
    pm.set_channel_map( channelMapB );
    pm.get<float>( _T("Float") ) = 3.75f;
    pm.get<int>( _T("Int") ) = 77;
    pm2.set_channel_map( channelMapC );
    pm2.get<float>( _T("MoreFloat") ) = -1.5f;
    pm2.get<int>( _T("MoreInt") ) = -8;
    pm2.get<frantic::tstring>( _T("Str") ) = _T("string here");
    pm.merge_property_map( pm2 ); // MERGE
    pm2.get<frantic::tstring>( _T("Str") ) =
        _T("string no longer here"); // reset the value of the property in the original pm2
    EXPECT_TRUE( pm.get<float>( _T("Float") ) == 3.75f );
    EXPECT_TRUE( pm.get<int>( _T("Int") ) == 77 );
    EXPECT_TRUE( pm.get<float>( _T("MoreFloat") ) == -1.5f );
    EXPECT_TRUE( pm.get<int>( _T("MoreInt") ) == -8 );
    EXPECT_TRUE( pm.get<frantic::tstring>( _T("Str") ) == _T("string here") );
}

TEST( ChannelMap, parse_channel_value_from_string_with_german_locale ) {

#ifdef _WIN32
#if defined( _MSC_VER ) && _MSC_VER >= 1700
    const char* german = "de-DE";
#else
    const char* german = "German";
#endif
#else
    const char* german = "de_DE.UTF-8";
#endif
    frantic::locale::set_locale_in_scope setLocale( german );

    const std::string s = "12.3";

    float f;
    parse_channel_value_from_string( channel_data_type_traits<float>::data_type(), s, reinterpret_cast<char*>( &f ) );
    EXPECT_EQ( f, 12.3f );

    double d;
    parse_channel_value_from_string( channel_data_type_traits<double>::data_type(), s, reinterpret_cast<char*>( &d ) );
    EXPECT_EQ( d, 12.3 );
}

TEST( ChannelMap, FrostFriendly ) {

    channel_map channelMap1;
    channelMap1.define_channel<vector3f>( _T("Position") );
    channelMap1.define_channel<vector3f>( _T("Velocity") );
    channelMap1.define_channel<half>( _T("U") );
    channelMap1.define_channel<boost::uint32_t>( _T("ID") );
    channelMap1.define_channel<float>( _T("Radius") );
    channelMap1.define_channel<boost::uint64_t>( _T("uint64_2") );
    channelMap1.define_channel<double>( _T("DoublePrecision") );
    channelMap1.define_channel<boost::uint8_t>( _T("SmallInt") );
    channelMap1.define_channel<vector3f>( _T("Acceleration") );
    channelMap1.define_channel<boost::uint64_t>( _T("uint64_1") );
    channelMap1.define_channel<int>( _T("integer") );
    channelMap1.define_channel<half>( _T("V") );
    channelMap1.define_channel<float>( _T("float") );

    channel_map channelMap = frantic::volumetrics::implicitsurface::create_optimized_channel_map( channelMap1 );

    // Check that the offsets we got accord with our expectations.  The channels have been sorted in descending order of
    // type size, using a stable sort.
    EXPECT_EQ( channelMap.channel_offset( _T("Position") ), 0 );
    EXPECT_EQ( channelMap.channel_offset( _T("Radius") ), 12 );
    EXPECT_EQ( channelMap.channel_offset( _T("DoublePrecision") ), 16 );
    EXPECT_EQ( channelMap.channel_offset( _T("uint64_2") ), 24 );
    EXPECT_EQ( channelMap.channel_offset( _T("uint64_1") ), 32 );
    EXPECT_EQ( channelMap.channel_offset( _T("Velocity") ), 40 );
    EXPECT_EQ( channelMap.channel_offset( _T("Acceleration") ), 52 );
    EXPECT_EQ( channelMap.channel_offset( _T("float") ), 64 );
    EXPECT_EQ( channelMap.channel_offset( _T("ID") ), 68 );
    EXPECT_EQ( channelMap.channel_offset( _T("integer") ), 72 );
    EXPECT_EQ( channelMap.channel_offset( _T("U") ), 76 );
    EXPECT_EQ( channelMap.channel_offset( _T("V") ), 78 );
    EXPECT_EQ( channelMap.channel_offset( _T("SmallInt") ), 80 );

    // Make sure that the total byte size is right
    EXPECT_EQ( channelMap.structure_size(), 84 );

    // Make sure the channel count is right
    EXPECT_EQ( channelMap.channel_count(), 13 );

    // When accessed by index, the channels should be in the order of their defined position
    EXPECT_TRUE( channelMap[0].name() == _T("Position") );
    EXPECT_TRUE( channelMap[1].name() == _T("Radius") );
    EXPECT_TRUE( channelMap[2].name() == _T("DoublePrecision") );
    EXPECT_TRUE( channelMap[3].name() == _T("uint64_2") );
    EXPECT_TRUE( channelMap[4].name() == _T("uint64_1") );
    EXPECT_TRUE( channelMap[5].name() == _T("Velocity") );
    EXPECT_TRUE( channelMap[6].name() == _T("Acceleration") );
    EXPECT_TRUE( channelMap[7].name() == _T("float") );
    EXPECT_TRUE( channelMap[8].name() == _T("ID") );
    EXPECT_TRUE( channelMap[9].name() == _T("integer") );
    EXPECT_TRUE( channelMap[10].name() == _T("U") );
    EXPECT_TRUE( channelMap[11].name() == _T("V") );
    EXPECT_TRUE( channelMap[12].name() == _T("SmallInt") );
}
