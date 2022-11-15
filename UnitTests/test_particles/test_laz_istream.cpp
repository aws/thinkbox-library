// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#ifdef LASZIP_LIB_AVAILABLE

#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/laz_particle_istream.hpp>

using namespace frantic;
using frantic::channels::channel_map;
using graphics::vector3f;
using graphics::vector3fd;
using particles::particle_array;
using particles::streams::laz_particle_istream;
using particles::streams::particle_istream_ptr;

TEST( LazIStream, ParticleCount ) {
    particles::particle_file_stream_factory_object streamFactory;
    particle_istream_ptr lazStream = streamFactory.create_istream( _T("TestInputs/ot_35120B4116B_1.laz") );

    std::size_t counter = 0;

    std::vector<char> buffer;
    buffer.resize( lazStream->particle_size() );

    EXPECT_EQ( lazStream->particle_count(), lazStream->particle_count_left() );

    while( lazStream->get_particle( buffer ) ) {
        EXPECT_EQ( counter, lazStream->particle_index() );
        EXPECT_EQ( counter, lazStream->particle_progress_index() );
        ++counter;
    }

    EXPECT_EQ( 0, lazStream->particle_count_left() );
    EXPECT_EQ( counter, lazStream->particle_count() );
    std::cout << "Read in " << counter << " particles." << std::endl;
}

TEST( LazIStream, 32bit ) {
    // Read in the particles with double precision.
    particle_array laz64bitParticles;
    {
        particles::particle_file_stream_factory_object stream64bitFactory;
        stream64bitFactory.set_position_type_hint( channels::data_type_float64 );

        particle_istream_ptr laz64bitStream =
            stream64bitFactory.create_istream( _T("TestInputs/ot_35120B4116B_1.laz") );
        laz64bitParticles.reset( laz64bitStream->get_channel_map() );

        std::vector<char> buffer;
        buffer.resize( laz64bitStream->particle_size() );

        while( laz64bitStream->get_particle( buffer ) ) {
            laz64bitParticles.push_back( &buffer[0] );
        }
        laz64bitStream->close();
    }

    // Read in the particles with single precision.
    particle_array laz32bitParticles;
    {
        particles::particle_file_stream_factory_object stream32bitFactory;
        stream32bitFactory.set_position_type_hint( channels::data_type_float32 );

        particle_istream_ptr laz32bitStream =
            stream32bitFactory.create_istream( _T("TestInputs/ot_35120B4116B_1.laz") );
        laz32bitParticles.reset( laz32bitStream->get_channel_map() );

        std::vector<char> buffer;
        buffer.resize( laz32bitStream->particle_size() );

        while( laz32bitStream->get_particle( buffer ) ) {
            laz32bitParticles.push_back( &buffer[0] );
        }
        laz32bitStream->close();
    }

    // Verify the particles converted properly.
    ASSERT_EQ( laz64bitParticles.size(), laz32bitParticles.size() );
    channels::channel_accessor<vector3fd> pos64Accessor =
        laz64bitParticles.get_channel_map().get_accessor<vector3fd>( _T("Position") );
    channels::channel_accessor<vector3f> pos32Accessor =
        laz32bitParticles.get_channel_map().get_accessor<vector3f>( _T("Position") );

    for( std::size_t i = 0, iEnd = laz64bitParticles.size(); i < iEnd; ++i ) {
        vector3fd pos64 = pos64Accessor( laz64bitParticles[i] );
        vector3f pos32 = pos32Accessor( laz32bitParticles[i] );

        ASSERT_EQ( static_cast<float>( pos64.x ), pos32.x );
        ASSERT_EQ( static_cast<float>( pos64.y ), pos32.y );
        ASSERT_EQ( static_cast<float>( pos64.z ), pos32.z );
    }
}

TEST( LazIStream, ScaleOffset ) {
    // The scale_offset_laz.laz file includes a scale factor and offset in it

    particles::particle_file_stream_factory_object factory;
    channel_map cm;
    cm.define_channel<vector3f>( _T("Position") );
    cm.end_channel_definition();
    particle_array prt( cm );

    prt.insert_particles( factory.create_istream( _T("TestInputs/scale_offset_laz.laz") ) );
    ASSERT_EQ( 5u, prt.size() );
    const vector3f* data = reinterpret_cast<const vector3f*>( prt[0] );
    EXPECT_EQ( vector3f( 0, 0, 0 ), data[0] );
    EXPECT_EQ( vector3f( 1, 0, 0 ), data[1] );
    EXPECT_EQ( vector3f( 0, 1, 0 ), data[2] );
    EXPECT_EQ( vector3f( 0, 0, 1 ), data[3] );
    EXPECT_EQ( vector3f( 2, 3, 4 ), data[4] );
}

#endif // LASZIP_LIB_AVAILABLE
