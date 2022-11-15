// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/*
 *	Test e57 particle stream.
 */

// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( E57_AVAILABLE )

#include <frantic/channels/channel_accessor.hpp>
#include <frantic/channels/channel_map.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/e57_particle_istream.hpp>
#include <frantic/strings/tstring.hpp>

#include "UnitTests/gtest-helper.h"
#include "utilities/e57_generator.hpp"

using namespace std;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::particles;
using namespace frantic::particles::streams;
using frantic::particles::particle_array;
using frantic::particles::particle_file_stream_factory_object;

// Set generate_files to true to generate fresh files when tests are run. For use when e57_generator.cpp is updated, ie.
// to create new files or update existing ones.
// NOTE: The current implementation of the prt2_ostream class automatically adds interpretation metadata when generating
//       the prt file. However the .e57 being generated does not have this metadata, so if the files are generated, the
//       tests below will all fail until the prt2_ostream is changed.
bool generateFiles = false;

// Test particle stream where e57 file has Position, Intensity, and Color.
TEST( E57Stream, NoNormals ) {
    if( generateFiles )
        generate_files( "TestInputs/noNormals", true, true, false, false, 1 );
    EXPECT_PARTICLE_FILE_EQ( _T("TestInputs/noNormals.prt"), _T("TestInputs/noNormals.e57") );
}

// Files with Position, Intensity, Color, and Normals.
TEST( E57Stream, Normals ) {
    if( generateFiles )
        generate_files( "TestInputs/wnormals", true, true, true, false, 1 );
    EXPECT_PARTICLE_FILE_EQ( _T("TestInputs/wnormals.prt"), _T("TestInputs/wnormals.e57") );
}

// Test e57 istream when some position values have double precision.
// In order for position channel to have double precision, must set factory position_type_hint.
TEST( E57Stream, DoubleFloatPos ) {
    if( generateFiles )
        generate_files( "TestInputs/doubleFloat", true, true, true, true, 1 );
    EXPECT_PARTICLE_FILE_EQ_POS_HINT( _T("TestInputs/doubleFloat.prt"), _T("TestInputs/doubleFloat.e57"),
                                      data_type_float64 );
}

// Files with 2 scans
TEST( E57Stream, twoScans ) {
    if( generateFiles )
        generate_files( "TestInputs/twoScans", true, true, true, false, 2 );
    EXPECT_PARTICLE_FILE_EQ( _T("TestInputs/twoScans.prt"), _T("TestInputs/twoScans.e57") );
}

// Files where the two scans have different transforms
TEST( E57Stream, transforms ) {
    if( generateFiles )
        generate_files( "TestInputs/transforms", true, true, true, false, 2, false, true );
    EXPECT_PARTICLE_FILE_EQ( _T("TestInputs/transforms.prt"), _T("TestInputs/transforms.e57") );
}

// Scaled int node is used for .e57 positions
TEST( E57Stream, ScaledInteger ) {
    if( generateFiles )
        generate_files( "TestInputs/scaledInteger", true, true, true, false, 1, true );
    EXPECT_PARTICLE_FILE_EQ( _T("TestInputs/scaledInteger.prt"), _T("TestInputs/scaledInteger.e57") );
}

// Color is int from 0 to 255
TEST( E57Stream, IntColor ) {
    if( generateFiles )
        generate_files( "TestInputs/intColor", true, true, true, false, 1, false, false, true );
    EXPECT_PARTICLE_FILE_EQ( _T("TestInputs/intColor.prt"), _T("TestInputs/intColor.e57") );
}

// Load a single scan, untransformed
TEST( E57Stream, LoadSingleScanUntransformed ) {
    const std::size_t particleCount = 10;
    const frantic::tstring scanGuid = _T("04FDBC32-C52A-491D-8248-D8DE54469E12");
    // reference data
    const double cartesianX[particleCount] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    const double cartesianY[particleCount] = { 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1 };
    const double cartesianZ[particleCount] = { 1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2, 9.1, 10.2 };

    const frantic::tstring fileName = _T("TestInputs/loadSingleScan");
    generate_files( frantic::strings::to_string( fileName ), true, true, true, false, 4, false, true, true );
    channel_map map;
    map.define_channel<vector3fd>( _T("Position") );
    map.define_channel<boost::uint32_t>( _T("ScannerIndex") );
    map.end_channel_definition();

    std::size_t numParticles = particleCount;
    boost::scoped_array<char> particleBuffer( new char[map.structure_size() * numParticles] );
    e57_particle_istream is( fileName + _T(".e57"), scanGuid, false );
    is.set_channel_map( map );

    ASSERT_TRUE( is.get_particles( &particleBuffer[0], numParticles ) );
    ASSERT_TRUE( numParticles == particleCount );

    channel_const_cvt_accessor<vector3f> posAccessor = map.get_const_cvt_accessor<vector3f>( _T("Position") );
    channel_const_cvt_accessor<boost::uint32_t> scanIndexAccessor =
        map.get_const_cvt_accessor<boost::uint32_t>( _T("ScannerIndex") );
    for( int i = 0; i < numParticles; ++i ) {
        const vector3f position = posAccessor.get( &particleBuffer[i * map.structure_size()] );
        EXPECT_FLOAT_EQ( (float)cartesianX[i], position.x );
        EXPECT_FLOAT_EQ( (float)cartesianY[i], position.y );
        EXPECT_FLOAT_EQ( (float)cartesianZ[i], position.z );
        EXPECT_EQ( 0, scanIndexAccessor.get( &particleBuffer[i * map.structure_size()] ) );
    }

    frantic::particles::particle_file_metadata metadata = is.get_metadata();
    frantic::channels::property_map pm = metadata.get_general_metadata();
    vector<transform4fd> scannerTransforms;
    frantic::particles::prt::get_scanner_transforms( pm, scannerTransforms );
    EXPECT_EQ( 0, scannerTransforms.size() );
}

TEST( E57Stream, LoadSingleScanTransformed ) {
    const std::size_t particleCount = 10;
    const boost::uint32_t scanIndex = 3;
    const frantic::tstring scanGuid = _T("E47BBBF3-C4A5-43E4-BC59-963545A6AD8A");
    // reference data
    const double cartesianX[particleCount] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    const double cartesianY[particleCount] = { 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1 };
    const double cartesianZ[particleCount] = { 1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2, 9.1, 10.2 };
    const vector3fd translationTransform = vector3fd( 1.0, 0.2, 3.0 ) * ( scanIndex - 1 );

    const frantic::tstring fileName = _T("TestInputs/loadSingleScanT");
    generate_files( frantic::strings::to_string( fileName ), true, true, true, false, 4, false, true, true );
    channel_map map;
    map.define_channel<vector3fd>( _T("Position") );
    map.define_channel<boost::uint32_t>( _T("ScannerIndex") );
    map.end_channel_definition();

    std::size_t numParticles = particleCount;
    boost::scoped_array<char> particleBuffer( new char[map.structure_size() * numParticles] );
    e57_particle_istream is( fileName + _T(".e57"), scanGuid, true );
    is.set_channel_map( map );
    ASSERT_TRUE( is.get_particles( &particleBuffer[0], numParticles ) );
    ASSERT_TRUE( numParticles == particleCount );

    channel_const_cvt_accessor<vector3f> posAccessor = map.get_const_cvt_accessor<vector3f>( _T("Position") );
    channel_const_cvt_accessor<boost::uint32_t> scanIndexAccessor =
        map.get_const_cvt_accessor<boost::uint32_t>( _T("ScannerIndex") );
    for( int i = 0; i < numParticles; ++i ) {
        const vector3f position = posAccessor.get( &particleBuffer[i * map.structure_size()] );
        EXPECT_FLOAT_EQ( float( cartesianX[i] + translationTransform.x ), position.x );
        EXPECT_FLOAT_EQ( float( cartesianY[i] + translationTransform.y ), position.y );
        EXPECT_FLOAT_EQ( float( cartesianZ[i] + translationTransform.z ), position.z );
        EXPECT_EQ( 1u, scanIndexAccessor.get( &particleBuffer[i * map.structure_size()] ) );
    }

    frantic::particles::particle_file_metadata metadata = is.get_metadata();
    frantic::channels::property_map pm = metadata.get_general_metadata();
    vector<transform4fd> scannerTransforms;
    frantic::particles::prt::get_scanner_transforms( pm, scannerTransforms );
    EXPECT_EQ( 1, scannerTransforms.size() );
}

#endif
