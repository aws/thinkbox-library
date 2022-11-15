// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/algorithm/string/join.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/ptx_particle_istream.hpp>

#include "gtest-helper.h"

using namespace std;
using namespace boost;
using namespace frantic;
using frantic::channels::channel_map;
using frantic::files::scoped_file_cleanup;
using frantic::graphics::vector3f;
using frantic::graphics::vector3fd;
using frantic::particles::particle_array;
using frantic::particles::particle_file_stream_factory_object;
namespace fs = boost::filesystem;

TEST( PTX, OneScan ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.ptx" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "1\n", f );
    std::fputs( "1\n", f );
    std::fputs( "0.0 0.0 0.0\n", f );
    std::fputs( "1.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0\n", f );
    std::fputs( "1.0 0.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 0.0 1.0\n", f );
    std::fputs( "1.2 3.4 5.6\n", f );
    f.close();

    frantic::particles::particle_istream_ptr pin(
        new frantic::particles::streams::ptx_particle_istream( frantic::files::to_tstring( path ) ) );

    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAcc(
        pin->get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") ) );
    std::vector<char> buffer( pin->get_channel_map().structure_size() );

    ASSERT_TRUE( pin->get_particle( buffer ) );
    EXPECT_EQ( positionAcc.get( buffer ), frantic::graphics::vector3f( 1.2f, 3.4f, 5.6f ) );
}

TEST( PTX, TransformedScan ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.ptx" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "1\n", f );
    std::fputs( "1\n", f );
    std::fputs( "0.0 0.0 0.0\n", f );
    std::fputs( "1.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0\n", f );
    std::fputs( "0.0 -1.0 0.0 0.0\n", f );
    std::fputs( "1.0 0.0 0.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0 1.0\n", f );
    std::fputs( "1.0 2.0 3.0\n", f );
    f.close();

    frantic::particles::particle_istream_ptr pin(
        new frantic::particles::streams::ptx_particle_istream( frantic::files::to_tstring( path ) ) );

    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAcc(
        pin->get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") ) );
    std::vector<char> buffer( pin->get_channel_map().structure_size() );

    ASSERT_TRUE( pin->get_particle( buffer ) );
    EXPECT_EQ( positionAcc.get( buffer ), frantic::graphics::vector3f( 2.f, -1.f, 4.f ) );
}

TEST( PTX, TwoScans ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.ptx" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "2\n", f );
    std::fputs( "1\n", f );
    std::fputs( "0.0 0.0 0.0\n", f );
    std::fputs( "1.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0\n", f );
    std::fputs( "1.0 0.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 0.0 1.0\n", f );
    std::fputs( "0.0 1.0 2.0\n", f );
    std::fputs( "3.0 4.0 5.0\n", f );
    std::fputs( "1\n", f );
    std::fputs( "1\n", f );
    std::fputs( "0.0 0.0 0.0\n", f );
    std::fputs( "1.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0\n", f );
    std::fputs( "1.0 0.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 0.0 1.0\n", f );
    std::fputs( "1.2 3.4 5.6\n", f );
    f.close();

    frantic::particles::particle_istream_ptr pin(
        new frantic::particles::streams::ptx_particle_istream( frantic::files::to_tstring( path ) ) );

    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAcc(
        pin->get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") ) );
    std::vector<char> buffer( pin->get_channel_map().structure_size() );

    ASSERT_TRUE( pin->get_particle( buffer ) );
    EXPECT_EQ( positionAcc.get( buffer ), frantic::graphics::vector3f( 0.f, 1.f, 2.f ) );

    ASSERT_TRUE( pin->get_particle( buffer ) );
    EXPECT_EQ( positionAcc.get( buffer ), frantic::graphics::vector3f( 3.f, 4.f, 5.f ) );

    ASSERT_TRUE( pin->get_particle( buffer ) );
    EXPECT_EQ( positionAcc.get( buffer ), frantic::graphics::vector3f( 1.2f, 3.4f, 5.6f ) );
}

TEST( PTX, Float64Position ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.ptx" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "1\n", f );
    std::fputs( "1\n", f );
    std::fputs( "0.0 0.0 0.0\n", f );
    std::fputs( "1.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0\n", f );
    std::fputs( "1.0 0.0 0.0 0.0\n", f );
    std::fputs( "0.0 1.0 0.0 0.0\n", f );
    std::fputs( "0.0 0.0 1.0 0.0\n", f );
    std::fputs( "0.0 0.0 0.0 1.0\n", f );
    std::fputs( "1.234567890123456 3.14159265358979 1.1e-200\n", f );
    f.close();

    // Test this through the factory object, set to request float64 Position data
    particles::particle_file_stream_factory_object pfactory;
    frantic::particles::particle_istream_ptr pin;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> positionAcc;
    vector<char> buffer;
    vector3fd v;

    // Read the file, and verify that the full precision was retained
    pfactory.set_position_type_hint( channels::data_type_float64 );
    pin = pfactory.create_istream( frantic::files::to_tstring( path ) );
    positionAcc = pin->get_channel_map().get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
    buffer.resize( pin->get_channel_map().structure_size() );
    ASSERT_TRUE( pin->get_particle( buffer ) );
    v = positionAcc.get( buffer );
    EXPECT_DOUBLE_EQ( 1.234567890123456, v.x );
    EXPECT_DOUBLE_EQ( 3.14159265358979, v.y );
    EXPECT_DOUBLE_EQ( 1.1e-200, v.z );

    // Read it again, but with float32 precision, and confirm that it loses the precision
    pfactory.set_position_type_hint( channels::data_type_float32 );
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

TEST( PTX, ColorScaling ) {
    // Confirm that the RGB colors are scaled from 0-255 range down to 0.0-1.0 range

    particle_file_stream_factory_object factory;
    channel_map cm;
    cm.define_channel<vector3f>( _T("Position") );
    cm.define_channel<vector3f>( _T("Color") );
    cm.end_channel_definition();
    particle_array prt( cm );

    prt.insert_particles( factory.create_istream( _T("TestInputs/simple_xyzrgb_ptx.ptx") ) );
    ASSERT_EQ( 4u, prt.size() );
    const vector3f* data = reinterpret_cast<const vector3f*>( prt[0] );
    // Position
    EXPECT_EQ( vector3f( 0, 0, 0 ), data[0] );
    EXPECT_EQ( vector3f( 1, 0, 0 ), data[2] );
    EXPECT_EQ( vector3f( 0, 1, 0 ), data[4] );
    EXPECT_EQ( vector3f( 0, 0, 1 ), data[6] );
    // Color
    EXPECT_VECTOR3F_NEAR( vector3f( 255, 0, 128 ) / 255.f, data[1], 0.001f );
    EXPECT_VECTOR3F_NEAR( vector3f( 64, 192, 255 ) / 255.f, data[3], 0.001f );
    EXPECT_VECTOR3F_NEAR( vector3f( 50, 90, 100 ) / 255.f, data[5], 0.001f );
    EXPECT_VECTOR3F_NEAR( vector3f( 12, 22, 250 ) / 255.f, data[7], 0.001f );
}
