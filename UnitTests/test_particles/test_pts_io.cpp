// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/channels/channel_map.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/pts_particle_istream.hpp>

#include "gtest-helper.h"

using namespace std;
using namespace boost;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using frantic::files::scoped_file_cleanup;
using frantic::particles::particle_array;
using frantic::particles::particle_file_stream_factory_object;
namespace fs = boost::filesystem;

TEST( PTS, ColorScaling ) {
    // Confirm that the RGB colors are scaled from 0-255 range down to 0.0-1.0 range

    particle_file_stream_factory_object factory;
    channel_map cm;
    cm.define_channel<vector3f>( _T("Position") );
    cm.define_channel<vector3f>( _T("Color") );
    cm.end_channel_definition();
    particle_array prt( cm );

    prt.insert_particles( factory.create_istream( _T("TestInputs/simple_xyzrgb_pts.pts") ) );
    ASSERT_EQ( 5u, prt.size() );
    const vector3f* data = reinterpret_cast<const vector3f*>( prt[0] );
    // Position
    EXPECT_EQ( vector3f( 0, 0, 0 ), data[0] );
    EXPECT_EQ( vector3f( 1, 0, 0 ), data[2] );
    EXPECT_EQ( vector3f( 0, 1, 0 ), data[4] );
    EXPECT_EQ( vector3f( 0, 0, 1 ), data[6] );
    EXPECT_EQ( vector3f( 2, 3, 4 ), data[8] );
    // Color
    EXPECT_VECTOR3F_NEAR( vector3f( 0, 0, 0 ) / 255.f, data[1], 0.0001f );
    EXPECT_VECTOR3F_NEAR( vector3f( 128, 255, 255 ) / 255.f, data[3], 0.0001f );
    EXPECT_VECTOR3F_NEAR( vector3f( 255, 128, 0 ) / 255.f, data[5], 0.0001f );
    EXPECT_VECTOR3F_NEAR( vector3f( 0, 255, 128 ) / 255.f, data[7], 0.0001f );
    EXPECT_VECTOR3F_NEAR( vector3f( 63, 128, 255 ) / 255.f, data[9], 0.0001f );
}

TEST( PTS, Float64Position ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.pts" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "1\n", f );
    std::fputs( "1.234567890123456 3.14159265358979 1.1e-200\n", f );
    f.close();

    // Test this through the factory object, set to request float64 Position data
    particle_file_stream_factory_object pfactory;
    frantic::particles::particle_istream_ptr pin;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> positionAcc;
    vector<char> buffer;
    vector3fd v;

    // Read the file, and verify that the full precision was retained
    pfactory.set_position_type_hint( data_type_float64 );
    pin = pfactory.create_istream( frantic::files::to_tstring( path ) );
    positionAcc = pin->get_channel_map().get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
    buffer.resize( pin->get_channel_map().structure_size() );
    ASSERT_TRUE( pin->get_particle( buffer ) );
    v = positionAcc.get( buffer );
    EXPECT_DOUBLE_EQ( 1.234567890123456, v.x );
    EXPECT_DOUBLE_EQ( 3.14159265358979, v.y );
    EXPECT_DOUBLE_EQ( 1.1e-200, v.z );

    // Read it again, but with float32 precision, and confirm that it loses the precision
    pfactory.set_position_type_hint( data_type_float32 );
    pin = pfactory.create_istream( frantic::files::to_tstring( path ) );
    positionAcc = pin->get_channel_map().get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
    buffer.resize( pin->get_channel_map().structure_size() );
    ASSERT_TRUE( pin->get_particle( buffer ) );
    v = positionAcc.get( buffer );
    EXPECT_NE( 1.234567890123456, v.x );
    EXPECT_NE( 3.14159265358979, v.y );
    EXPECT_NE( 1.1e-200, v.z );
    // they should be equal with a fuzzy float32 comparison though
    EXPECT_FLOAT_EQ( 1.234567890123456f, (float)v.x );
    EXPECT_FLOAT_EQ( 3.14159265358979f, (float)v.y );
    EXPECT_EQ( 0.f, v.z );
}

TEST( PTS, InconsistentFieldCount ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.pts" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "7\n", f );
    std::fputs( "1 2 3\n", f );
    std::fputs( "2 2 3 0.5\n", f );
    std::fputs( "3 2 3 0.5 255 0 255\n", f );
    std::fputs( "3 2 3 0.5 0.0 1.0 0.0\n", f );
    f.close();

    // Test this through the factory object, set to request float64 Position data
    particle_file_stream_factory_object pfactory;
    frantic::particles::particle_istream_ptr pin;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> positionAcc;
    pin = pfactory.create_istream( frantic::files::to_tstring( path ) );

    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );
    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Intensity") ) );
    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Color") ) );
    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Normal") ) );

    // TODO: Could test that the values come out right too.
}
