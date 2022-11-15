// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <fstream>

#include <boost/make_shared.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/las_particle_istream.hpp>

#include "gtest-helper.h"

using frantic::channels::channel_map;
using frantic::graphics::vector3f;
using frantic::particles::particle_array;
using frantic::particles::particle_file_stream_factory_object;

namespace {

template <class T>
void write( std::ofstream& out, const T& value ) {
    out.write( reinterpret_cast<const char*>( &value ), sizeof( value ) );
}

} // anonymous namespace

// Test an LAS file with extra per-point data
TEST( LASParticleIStream, ExtraBytes ) {
    using namespace boost::filesystem;
    using namespace frantic::channels;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::particles::streams;

    path tempDir = temp_directory_path() / unique_path();
    create_directory( tempDir );
    path filename = tempDir / "extraBytes.las";

    frantic::files::scoped_file_cleanup cleanup;
    cleanup.add( tempDir );

    { // scope for writing file
        const boost::uint16_t pointDataRecordLength = 40;
        const std::size_t particleCount = 2;

        std::ofstream out( frantic::files::to_tstring( filename ).c_str(), std::ios::binary );

        // Write LAS 1.2 header, using Point Record Format 3,
        // and with "extra bytes" in the point data record
        // (40 bytes vs 34 required by the Point Record Format).

        // signature
        out.write( "LASF", 4 );
        // file source ID
        write( out, boost::uint16_t( 0 ) );
        // global encoding
        write( out, boost::uint16_t( 0 ) );
        // project ID
        const char projectId[16] = {};
        out.write( projectId, sizeof( projectId ) );
        // version major
        write( out, boost::uint8_t( 1 ) );
        // version minor
        write( out, boost::uint8_t( 2 ) );
        // system identifier
        const char systemIdentifier[32] = {};
        out.write( systemIdentifier, sizeof( systemIdentifier ) );
        // generating software
        const char generatingSoftware[32] = {};
        out.write( generatingSoftware, sizeof( generatingSoftware ) );
        // file creation day of year
        write( out, boost::uint16_t( 0 ) );
        // file creation year
        write( out, boost::uint16_t( 0 ) );
        // header size
        write( out, boost::uint16_t( 227 ) );
        // offset to point data
        write( out, boost::uint32_t( 227 ) );
        // number of variable length records
        write( out, boost::uint32_t( 0 ) );
        // point data record format
        write( out, boost::uint8_t( 3 ) );
        // point data record length
        write( out, boost::uint16_t( pointDataRecordLength ) );
        // legacy number of point records
        write( out, boost::uint32_t( particleCount ) );
        // legacy number of points by return
        for( int i = 0; i < 5; ++i ) {
            write( out, boost::uint32_t( i == 0 ? particleCount : 0 ) );
        }
        // X scale factor
        write( out, double( 1 ) );
        // Y scale factor
        write( out, double( 1 ) );
        // Z scale factor
        write( out, double( 1 ) );
        // X offset
        write( out, double( 0 ) );
        // Y offset
        write( out, double( 0 ) );
        // Z offset
        write( out, double( 0 ) );
        // max X
        write( out, double( 1 ) );
        // min X
        write( out, double( 0 ) );
        // max Y
        write( out, double( 0 ) );
        // min Y
        write( out, double( 0 ) );
        // max Z
        write( out, double( 0 ) );
        // min Z
        write( out, double( 0 ) );

        std::vector<char> buffer( pointDataRecordLength );

        channel_map channelMap;
        channelMap.define_channel( _T("Position"), 3, data_type_int32, 0 );
        channelMap.define_channel( _T("Intensity"), 1, data_type_uint16, 12 );
        channelMap.define_channel( _T("Color"), 3, data_type_uint16, 28 );
        channelMap.end_channel_definition( 1, true, false );

        ASSERT_TRUE( pointDataRecordLength >= channelMap.structure_size() );

        channel_accessor<vector3> positionAcc = channelMap.get_accessor<vector3>( _T("Position") );
        channel_accessor<boost::uint16_t> intensityAcc = channelMap.get_accessor<boost::uint16_t>( _T("Intensity") );
        channel_general_accessor colorAcc = channelMap.get_general_accessor( _T("Color") );

        intensityAcc( &buffer[0] ) = 65535;
        {
            boost::uint16_t color[3] = { 65535, 0, 0 };
            colorAcc.set_channel_from_primitive( &buffer[0], reinterpret_cast<const char*>( color ) );
        }

        out.write( &buffer[0], buffer.size() );

        positionAcc( &buffer[0] ) = vector3( 1, 0, 0 );
        intensityAcc( &buffer[0] ) = 0;
        {
            boost::uint16_t color[3] = { 0, 65535, 0 };
            colorAcc.set_channel_from_primitive( &buffer[0], reinterpret_cast<const char*>( color ) );
        }

        out.write( &buffer[0], buffer.size() );
    }

    particle_istream_ptr pin = boost::make_shared<las_particle_istream>( frantic::files::to_tstring( filename ) );
    ASSERT_TRUE( pin != 0 );

    const channel_map& channelMap = pin->get_channel_map();

    channel_cvt_accessor<vector3f> positionAcc = channelMap.get_cvt_accessor<vector3f>( _T("Position") );
    channel_cvt_accessor<float> intensityAcc = channelMap.get_cvt_accessor<float>( _T("Intensity") );
    channel_cvt_accessor<vector3f> colorAcc = channelMap.get_cvt_accessor<vector3f>( _T("Color") );

    std::vector<char> buffer( channelMap.structure_size() );

    bool gotParticle = false;

    gotParticle = pin->get_particle( buffer );
    ASSERT_TRUE( gotParticle );

    EXPECT_EQ( vector3f( 0, 0, 0 ), positionAcc( buffer ) );
    EXPECT_EQ( 1, intensityAcc( buffer ) );
    EXPECT_EQ( vector3f( 1, 0, 0 ), colorAcc( buffer ) );

    gotParticle = pin->get_particle( buffer );
    ASSERT_TRUE( gotParticle );

    EXPECT_EQ( vector3f( 1, 0, 0 ), positionAcc( buffer ) );
    EXPECT_EQ( 0, intensityAcc( buffer ) );
    EXPECT_EQ( vector3f( 0, 1, 0 ), colorAcc( buffer ) );
}

TEST( LasIStream, ScaleOffset ) {
    // The scale_offset_las.las file includes a scale factor and offset in it

    particle_file_stream_factory_object factory;
    channel_map cm;
    cm.define_channel<vector3f>( _T("Position") );
    cm.end_channel_definition();
    particle_array prt( cm );

    prt.insert_particles( factory.create_istream( _T("TestInputs/scale_offset_las.las") ) );
    ASSERT_EQ( 5u, prt.size() );
    const vector3f* data = reinterpret_cast<const vector3f*>( prt[0] );
    EXPECT_EQ( vector3f( 0, 0, 0 ), data[0] );
    EXPECT_EQ( vector3f( 1, 0, 0 ), data[1] );
    EXPECT_EQ( vector3f( 0, 1, 0 ), data[2] );
    EXPECT_EQ( vector3f( 0, 0, 1 ), data[3] );
    EXPECT_EQ( vector3f( 2, 3, 4 ), data[4] );
}

/**
 * Some software incorrectly writes LAS color values in the range [0, 255]
 * instead of the proper [0, 65535].  As a workaround, we attempt to detect
 * such files and scale them into the full range.
 */
TEST( LasIStream, ScaleSmallColorRange ) {
    using namespace frantic::channels;
    using namespace frantic::particles;

    channel_map cm;
    cm.define_channel<vector3f>( _T("Position") );
    cm.define_channel<vector3f>( _T("Color") );
    cm.end_channel_definition();

    particle_array pa( cm );

    pa.insert_particles( particle_file_istream_factory( _T("TestInputs/rgb_magnitude_127.las") ) );

    channel_accessor<vector3f> colorAcc = cm.get_accessor<vector3f>( _T("Color") );

    ASSERT_EQ( 4, pa.size() );
    EXPECT_VECTOR3F_NEAR( vector3f( 0, 0, 0 ), colorAcc( pa[0] ), 0.01f );
    EXPECT_VECTOR3F_NEAR( vector3f( 0.5, 0, 0 ), colorAcc( pa[1] ), 0.01f );
    EXPECT_VECTOR3F_NEAR( vector3f( 0, 0.5, 0 ), colorAcc( pa[2] ), 0.01f );
    EXPECT_VECTOR3F_NEAR( vector3f( 0, 0, 0.5 ), colorAcc( pa[3] ), 0.01f );
}
