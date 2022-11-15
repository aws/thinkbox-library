// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <fstream>

#include "tbb/task_scheduler_init.h"

#include <frantic/files/files.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/prt2_particle_istream.hpp>
#include <frantic/particles/streams/prt2_particle_ostream.hpp>
#include <frantic/prtfile/prt2_reader.hpp>
#include <frantic/prtfile/prt2_writer.hpp>

using namespace std;
using namespace boost;
using namespace frantic::prtfile;
using namespace frantic::particles::streams;
using namespace frantic::files;
using frantic::channels::channel_cvt_accessor;
using frantic::channels::channel_map;
using frantic::channels::property_map;
using frantic::graphics::boundbox3f;
using frantic::graphics::color3h;
using frantic::graphics::vector3f;
using frantic::particles::particle_array;

TEST( PRT2Stream, EmptyFile ) {
    tbb::task_scheduler_init taskScheduler;

    frantic::channels::channel_map cm;
    cm.define_channel<frantic::graphics::vector3f>( _T("Position") );
    cm.end_channel_definition();

    frantic::tstring testFile = _T("Test.PRT2Stream.EmptyFile.prt2");
    if( file_exists( testFile ) ) {
        delete_file( testFile );
    }

    frantic::particles::particle_ostream_ptr out(
        new frantic::particles::streams::prt2_particle_ostream( testFile, cm, cm ) );
}

void test_PRT2Stream_RoundTrip( prt2_compression_t compressionScheme ) {
    // Write some particles to a prt2 output stream
    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<int32_t>( _T("ID") );
    pcm.end_channel_definition();

    channel_map gmcm;
    gmcm.define_channel<vector3f>( _T("SomeVector") );
    gmcm.define_channel<double>( _T("LengthUnitInMicrometers") );
    gmcm.end_channel_definition();

    channel_map idcm;
    idcm.define_channel<int32_t>( _T("FavoriteID") );
    idcm.end_channel_definition();

    property_map globalMetadata( gmcm );
    map<frantic::tstring, property_map> channelMetadata;
    globalMetadata.set_cvt<vector3f>( _T("SomeVector"), vector3f( 2, 4, 6 ) );
    globalMetadata.set_cvt<double>( _T("LengthUnitInMicrometers"), 1000.0 );
    channelMetadata.insert( pair<frantic::tstring, property_map>( _T("ID"), property_map( idcm ) ) );
    channelMetadata[_T("ID")].set_cvt<int32_t>( _T("FavoriteID"), 42 );

    frantic::tstring testFile = _T("Test.PRT2Stream.RoundTrip.prt2");
    if( file_exists( testFile ) ) {
        delete_file( testFile );
    }
    prt2_particle_ostream pout( testFile, pcm, pcm, compressionScheme, true, boost::filesystem::path(), &globalMetadata,
                                &channelMetadata, 3 * ( 12 + 4 ) );

    channel_cvt_accessor<vector3f> posAcc = pcm.get_cvt_accessor<vector3f>( _T("Position") );
    channel_cvt_accessor<int32_t> idAcc = pcm.get_cvt_accessor<int32_t>( _T("ID") );
    vector<char> buf( pcm.structure_size() );
    for( int i = 0; i < 13; ++i ) {
        posAcc.set( &buf[0], vector3f( (float)i, -(float)i, 1.5f ) );
        idAcc.set( &buf[0], i + 40 );
        pout.put_particle( &buf[0] );
    }

    pout.close();

    // Verify the particle file using the raw prt2_reader object
    prt2_reader prt2;
    prt2.open( testFile );
    EXPECT_EQ( 13, prt2.get_particle_count() );
    EXPECT_EQ( 5, prt2.get_particle_chunk_count() );

    // Verify the channel map
    EXPECT_EQ( 2u, prt2.get_channel_count() );
    EXPECT_TRUE( prt2.get_channel_map().has_channel( _T("Position") ) );
    EXPECT_TRUE( prt2.get_channel_map().has_channel( _T("ID") ) );

    // Verify the global metadata
    EXPECT_EQ( 2u, prt2.get_general_metadata().get_channel_map().channel_count() );
    EXPECT_EQ( vector3f( 2, 4, 6 ), prt2.get_general_metadata().get<vector3f>( _T("SomeVector") ) );
    EXPECT_EQ( 1000.f, prt2.get_general_metadata().get<double>( _T("LengthUnitInMicrometers") ) );

    // Verify the channel metadata
    EXPECT_EQ( 1u, prt2.get_channel_metadata( _T("Position") ).get_channel_map().channel_count() );
    boundbox3f bb = prt2.get_channel_metadata( _T("Position") ).get<boundbox3f>( _T("Extents") );
    EXPECT_EQ( bb.minimum(), vector3f( 0.f, -12.f, 1.5f ) );
    EXPECT_EQ( bb.maximum(), vector3f( 12.f, 0.f, 1.5f ) );
    EXPECT_EQ( 1u, prt2.get_channel_metadata( _T("ID") ).get_channel_map().channel_count() );
    EXPECT_EQ( 42, prt2.get_channel_metadata( _T("ID") ).get<int32_t>( _T("FavoriteID") ) );

    // Verify the combined metadata
    { // scope for metadata
        // The metadata pointers aren't required to be equal, but they are
        // in the current implementation.  I'm taking advantage of this so
        // I can compare the pointers instead of the individual data items.
        const frantic::particles::particle_file_metadata& metadata = prt2.get_metadata();
        EXPECT_EQ( &prt2.get_general_metadata(), &metadata.get_general_metadata() );
        EXPECT_EQ( &prt2.get_channel_metadata( _T("Position") ), metadata.get_channel_metadata( _T("Position") ) );
        EXPECT_EQ( &prt2.get_channel_metadata( _T("ID") ), metadata.get_channel_metadata( _T("ID") ) );
    }

    // Verify the particles
    particle_array parr;
    prt2.read_particles( parr );
    ASSERT_EQ( 13u, parr.size() );
    posAcc = parr.get_channel_map().get_cvt_accessor<vector3f>( _T("Position") );
    idAcc = parr.get_channel_map().get_cvt_accessor<int32_t>( _T("ID") );
    for( int i = 0; i < 13; ++i ) {
        EXPECT_EQ( vector3f( (float)i, -(float)i, 1.5f ), posAcc( parr.at( i ) ) );
        EXPECT_EQ( i + 40, idAcc( parr.at( i ) ) );
    }

    prt2.close();

    // Verify the particles again using a prt2 input stream
    boost::shared_ptr<prt2_particle_istream> pin( new prt2_particle_istream( testFile ) );
    EXPECT_EQ( 13, pin->particle_count() );
    EXPECT_EQ( 13, pin->particle_count_left() );

    // Verify the channel map
    EXPECT_EQ( 2u, pin->get_channel_map().channel_count() );
    EXPECT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );
    EXPECT_TRUE( pin->get_channel_map().has_channel( _T("ID") ) );

    // Verify the global metadata
    EXPECT_EQ( 2u, pin->get_general_metadata().get_channel_map().channel_count() );
    EXPECT_EQ( vector3f( 2, 4, 6 ), pin->get_general_metadata().get<vector3f>( _T("SomeVector") ) );
    EXPECT_EQ( 1000.f, pin->get_general_metadata().get<double>( _T("LengthUnitInMicrometers") ) );

    // Verify the channel metadata
    EXPECT_EQ( 1u, pin->get_channel_metadata( _T("Position") )->get_channel_map().channel_count() );
    bb = pin->get_channel_metadata( _T("Position") )->get<boundbox3f>( _T("Extents") );
    EXPECT_EQ( bb.minimum(), vector3f( 0.f, -12.f, 1.5f ) );
    EXPECT_EQ( bb.maximum(), vector3f( 12.f, 0.f, 1.5f ) );
    EXPECT_EQ( 1u, pin->get_channel_metadata( _T("ID") )->get_channel_map().channel_count() );
    EXPECT_EQ( 42, pin->get_channel_metadata( _T("ID") )->get<int32_t>( _T("FavoriteID") ) );

    // Verify the particles
    parr.clear();
    parr.insert_particles( pin );
    ASSERT_EQ( 13u, parr.size() );
    posAcc = parr.get_channel_map().get_cvt_accessor<vector3f>( _T("Position") );
    idAcc = parr.get_channel_map().get_cvt_accessor<int32_t>( _T("ID") );
    for( int i = 0; i < 13; ++i ) {
        EXPECT_EQ( vector3f( (float)i, -(float)i, 1.5f ), posAcc( parr.at( i ) ) );
        EXPECT_EQ( i + 40, idAcc( parr.at( i ) ) );
    }

    pin->close();
}

