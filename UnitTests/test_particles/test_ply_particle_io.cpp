// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/filesystem.hpp>
#include <boost/scope_exit.hpp>

#include <frantic/files/files.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_file_io.hpp>
#include <frantic/particles/streams/ply_particle_istream.hpp>
#include <frantic/particles/streams/ply_particle_ostream.hpp>

using frantic::channels::channel_accessor;
using frantic::graphics::vector3f;
using frantic::particles::particle_istream_ptr;
using frantic::particles::particle_ostream_ptr;
using frantic::particles::streams::ply_particle_istream;
using frantic::particles::streams::ply_particle_ostream;

namespace {

class temp_file {
  public:
    temp_file()
        : m_path( boost::filesystem::unique_path(
              ( boost::filesystem::temp_directory_path() / L"%%%%-%%%%-%%%%-%%%%.ply" ).wstring() ) ) {}

    ~temp_file() {
        boost::system::error_code ec;
        boost::filesystem::remove( m_path, ec );
    }

    frantic::tstring get_path() const { return frantic::files::to_tstring( m_path ); }

  private:
    boost::filesystem::path m_path;
};

} // anonymous namespace

TEST( PLYParticle, WriteAndReadNoParticles ) {
    temp_file file;

    { // scope for output stream
        frantic::channels::channel_map channelMap;
        channelMap.define_channel<vector3f>( _T( "Position" ) );
        channelMap.end_channel_definition();

        particle_ostream_ptr pout( new ply_particle_ostream( file.get_path(), channelMap, channelMap ) );
        pout->close();
    }

    particle_istream_ptr pin( new ply_particle_istream( file.get_path() ) );

    ASSERT_TRUE( pin != 0 );
    EXPECT_EQ( pin->particle_count(), 0 );
    EXPECT_TRUE( pin->get_channel_map().has_channel( _T( "Position" ) ) );

    std::vector<char> buffer( pin->get_channel_map().structure_size() );

    bool gotParticle = pin->get_particle( buffer );
    EXPECT_FALSE( gotParticle );
}

TEST( PLYParticle, WriteAndReadTwoPositions ) {
    temp_file file;

    { // scope for output stream
        frantic::channels::channel_map channelMap;
        channelMap.define_channel<vector3f>( _T( "Position" ) );
        channelMap.end_channel_definition();

        std::vector<char> buffer( channelMap.structure_size() );
        channel_accessor<vector3f> acc( channelMap.get_accessor<vector3f>( _T( "Position" ) ) );

        particle_ostream_ptr pout( new ply_particle_ostream( file.get_path(), channelMap, channelMap ) );

        acc( buffer ).set( 0, 1, 2 );
        pout->put_particle( buffer );

        acc( buffer ).set( 3, 4, 5 );
        pout->put_particle( buffer );

        pout->close();
    }

    particle_istream_ptr pin( new ply_particle_istream( file.get_path() ) );

    ASSERT_TRUE( pin != 0 );
    EXPECT_EQ( pin->particle_count(), 2 );
    EXPECT_TRUE( pin->get_channel_map().has_channel( _T( "Position" ) ) );

    std::vector<char> buffer( pin->get_channel_map().structure_size() );
    channel_accessor<vector3f> acc( pin->get_channel_map().get_accessor<vector3f>( _T( "Position" ) ) );

    bool gotParticle = pin->get_particle( buffer );
    EXPECT_TRUE( gotParticle );
    EXPECT_EQ( acc( buffer ), vector3f( 0, 1, 2 ) );

    gotParticle = pin->get_particle( buffer );
    EXPECT_TRUE( gotParticle );
    EXPECT_EQ( acc( buffer ), vector3f( 3, 4, 5 ) );

    gotParticle = pin->get_particle( buffer );
    EXPECT_FALSE( gotParticle );
}

TEST( PLYParticle, WriteAndReadAllChannels ) {
    temp_file file;

    { // scope for output stream
        frantic::channels::channel_map channelMap;
        channelMap.define_channel<vector3f>( _T( "Position" ) );
        channelMap.define_channel<vector3f>( _T( "Color" ) );
        channelMap.define_channel<vector3f>( _T( "TextureCoord" ) );
        channelMap.define_channel<vector3f>( _T( "Normal" ) );
        channelMap.end_channel_definition();

        std::vector<char> buffer( channelMap.structure_size() );
        channel_accessor<vector3f> positionAcc( channelMap.get_accessor<vector3f>( _T( "Position" ) ) );
        channel_accessor<vector3f> colorAcc( channelMap.get_accessor<vector3f>( _T( "Color" ) ) );
        channel_accessor<vector3f> textureAcc( channelMap.get_accessor<vector3f>( _T( "TextureCoord" ) ) );
        channel_accessor<vector3f> normalAcc( channelMap.get_accessor<vector3f>( _T( "Normal" ) ) );

        particle_ostream_ptr pout( new ply_particle_ostream( file.get_path(), channelMap, channelMap ) );

        positionAcc( buffer ).set( 0, 1, 2 );
        colorAcc( buffer ).set( 0, 1, 0 );
        textureAcc( buffer ).set( 6, 7, 8 );
        normalAcc( buffer ).set( 9, 10, 11 );

        pout->put_particle( buffer );

        pout->close();
    }

    particle_istream_ptr pin( new ply_particle_istream( file.get_path() ) );

    ASSERT_TRUE( pin != 0 );
    EXPECT_EQ( pin->particle_count(), 1 );

    std::vector<char> buffer( pin->get_channel_map().structure_size() );
    channel_accessor<vector3f> positionAcc( pin->get_channel_map().get_accessor<vector3f>( _T( "Position" ) ) );
    channel_accessor<vector3f> colorAcc( pin->get_channel_map().get_accessor<vector3f>( _T( "Color" ) ) );
    channel_accessor<vector3f> textureAcc( pin->get_channel_map().get_accessor<vector3f>( _T( "TextureCoord" ) ) );
    channel_accessor<vector3f> normalAcc( pin->get_channel_map().get_accessor<vector3f>( _T( "Normal" ) ) );

    bool gotParticle = pin->get_particle( buffer );
    EXPECT_TRUE( gotParticle );

    EXPECT_EQ( positionAcc( buffer ), vector3f( 0, 1, 2 ) );
    EXPECT_EQ( colorAcc( buffer ), vector3f( 0, 1, 0 ) );
    EXPECT_EQ( textureAcc( buffer ), vector3f( 6, 7, 8 ) );
    EXPECT_EQ( normalAcc( buffer ), vector3f( 9, 10, 11 ) );

    gotParticle = pin->get_particle( buffer );
    EXPECT_FALSE( gotParticle );
}

TEST( PLYParticle, WritePositionAndNormal ) {
    using frantic::geometry::polymesh3_vertex_accessor;

    temp_file file;

    const frantic::tstring pathString( file.get_path() );

    { // scope for ply_particle_ostream
        frantic::channels::channel_map channelMap;
        channelMap.define_channel<vector3f>( _T( "Position" ) );
        channelMap.define_channel<vector3f>( _T( "Normal" ) );
        channelMap.end_channel_definition();

        particle_ostream_ptr pout( new ply_particle_ostream( pathString, channelMap, channelMap ) );

        channel_accessor<vector3f> positionAcc( channelMap.get_accessor<vector3f>( _T( "Position" ) ) );
        channel_accessor<vector3f> normalAcc( channelMap.get_accessor<vector3f>( _T( "Normal" ) ) );

        std::vector<char> buffer( channelMap.structure_size() );

        positionAcc( buffer ).set( 0, 1, 2 );
        normalAcc( buffer ).set( 0, 0, 1 );

        pout->put_particle( &buffer[0] );

        pout->close();
    }

    frantic::geometry::polymesh3_ptr mesh = frantic::geometry::load_ply_polymesh_file( pathString );

    ASSERT_TRUE( mesh.get() != 0 );

    ASSERT_EQ( mesh->vertex_count(), 1 );
    EXPECT_EQ( mesh->face_count(), 0 );

    EXPECT_EQ( mesh->get_vertex( 0 ), vector3f( 0, 1, 2 ) );

    ASSERT_TRUE( mesh->has_vertex_channel( _T( "Normal" ) ) );

    polymesh3_vertex_accessor<vector3f> normalAcc = mesh->get_vertex_accessor<vector3f>( _T( "Normal" ) );
    ASSERT_TRUE( normalAcc.is_valid() );
    ASSERT_EQ( normalAcc.vertex_count(), 1 );
    EXPECT_EQ( normalAcc.get_vertex( 0 ), vector3f( 0, 0, 1 ) );
}
