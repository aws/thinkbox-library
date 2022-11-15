// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/algorithm/string/join.hpp>

#include <frantic/files/paths.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>

#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/ptg_particle_istream.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::particles::streams;
using namespace frantic::graphics;
using frantic::files::scoped_file_cleanup;
using frantic::files::to_tstring;
using frantic::strings::to_tstring;
namespace fs = boost::filesystem;

// A bit of helper code to create a one-particle PTG file
template <typename T>
static void write_to_file( T value, std::ofstream& fout ) {
    fout.write( reinterpret_cast<char*>( &value ), sizeof( value ) );
}
static void puts_to_ptg_file( std::string str, std::ofstream& fout ) {
    size_t temp = str.size() + 1;
    write_to_file<boost::int32_t>( (boost::int32_t)temp, fout );
    fout.write( str.c_str(), temp );
}
template <class VecType>
static void create_test_ptg_file( frantic::tstring fileName, const transform4fd& xform, const VecType& pos ) {
    std::ofstream fout( fileName.c_str(), std::ios::out | std::ios::binary );
    fout.write( "PTG", 4 );
    write_to_file<boost::int32_t>( 2458887111, fout );
    puts_to_ptg_file( "%%header_begin", fout );
    puts_to_ptg_file( "%%version", fout );
    write_to_file<boost::int32_t>( 1, fout );
    puts_to_ptg_file( "%%cols", fout );
    write_to_file<boost::int32_t>( 1, fout );
    puts_to_ptg_file( "%%rows", fout );
    write_to_file<boost::int32_t>( 1, fout );
    puts_to_ptg_file( "%%transform", fout );
    write_to_file<double>( xform[0], fout );
    write_to_file<double>( xform[1], fout );
    write_to_file<double>( xform[2], fout );
    write_to_file<double>( xform[3], fout );
    write_to_file<double>( xform[4], fout );
    write_to_file<double>( xform[5], fout );
    write_to_file<double>( xform[6], fout );
    write_to_file<double>( xform[7], fout );
    write_to_file<double>( xform[8], fout );
    write_to_file<double>( xform[9], fout );
    write_to_file<double>( xform[10], fout );
    write_to_file<double>( xform[11], fout );
    write_to_file<double>( xform[12], fout );
    write_to_file<double>( xform[13], fout );
    write_to_file<double>( xform[14], fout );
    write_to_file<double>( xform[15], fout );
    puts_to_ptg_file( "%%properties", fout );
    write_to_file<boost::int32_t>( ( sizeof( VecType ) == 12u ? ptg_particle_istream::PTG_POSITION_AS_FLOAT
                                                              : ptg_particle_istream::PTG_POSITION_AS_DOUBLE ) |
                                       ptg_particle_istream::PTG_INTENSITY | ptg_particle_istream::PTG_COLOR,
                                   fout );
    puts_to_ptg_file( "%%header_end", fout );
    int64_t colOffset0 = int64_t( fout.tellp() ) + 8;
    write_to_file<boost::int64_t>( colOffset0, fout );
    fout.put( 1 );
    write_to_file<typename VecType::float_type>( pos.x, fout );
    write_to_file<typename VecType::float_type>( pos.y, fout );
    write_to_file<typename VecType::float_type>( pos.z, fout );
    write_to_file<float>( 1.0f, fout );
    fout.put( (char)1 );
    fout.put( (char)255 );
    fout.put( (char)254 );
}

// creates a simple ptg file (that has a transfrom), reads it into ptg_particle_istream, get's the particles, asserts
// that the particles position is correct (i.e. the transform was applied correctly)
TEST( PTG, TransformedParticle ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.ptg" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    // transform matrix to rotate z 90 degrees
    create_test_ptg_file( to_tstring( path ), transform4fd( 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 2, 3, 1 ),
                          vector3f( 1, 2, 3 ) );

    frantic::particles::streams::ptg_particle_istream testStream( to_tstring( path ) );
    size_t numParticles = (size_t)testStream.particle_count_guess();
    std::vector<char> particleBuffer( numParticles * testStream.particle_size() );
    testStream.get_particles( &particleBuffer[0], numParticles );
    particleBuffer.resize( numParticles * testStream.particle_size() );
    // check particle
    const frantic::channels::channel_map& cm( testStream.get_native_channel_map() );
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector3f> positionAcc(
        cm.get_const_cvt_accessor<frantic::graphics::vector3f>( _T("Position") ) );
    const frantic::graphics::vector3f position = positionAcc( &particleBuffer[0] );
    EXPECT_EQ( position.x, 3 );
    EXPECT_EQ( position.y, 1 );
    EXPECT_EQ( position.z, 6 );
}

TEST( PTG, Float64Position ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.ptg" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    // Create a test file with identity transform and some double-precision values
    create_test_ptg_file( frantic::files::to_tstring( path ), transform4fd(),
                          vector3fd( 1.234567890123456, 3.14159265358979, 1.1e-200 ) );

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
