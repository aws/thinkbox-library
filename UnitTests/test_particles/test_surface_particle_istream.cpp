// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 *	Creates a surface_particle_istream and tests the points given by the
 *	stream to ensure they fall on the mesh and have the expected normals.
 *
 *	Tests the number of particles given by the stream to ensure that it
 *	matches the number specified.
 */

// clang-format off
#include <stdafx.h>
// clang-format on

#include <frantic/graphics/graphics_utils.hpp>
#include <frantic/particles/streams/surface_particle_istream.hpp>
#include <utilities/mesh_generators.hpp>

using namespace frantic::particles::streams;
using namespace frantic::geometry;

TEST( SurfaceParticleIStream, ParticleCount ) {

    // set up the mesh
    boost::shared_ptr<trimesh3> mesh( new trimesh3() );
    mesh->add_vertex( 0, 0, 0 );
    mesh->add_vertex( 1, 1, 0 );
    mesh->add_vertex( 1, 0, 0 );
    mesh->add_face( 0, 1, 2 );

    int numParticles = 100;
    boost::uint32_t randomSeed = 132234223;

    frantic::channels::channel_map channelMap;
    channelMap.define_channel( _T( "Position" ), 3, frantic::channels::data_type_float32 );
    channelMap.define_channel( _T( "Normal" ), 3, frantic::channels::data_type_float32 );
    channelMap.end_channel_definition();

    particle_istream_ptr particleStream =
        create_surface_particle_istream_using_count( *mesh, numParticles, randomSeed );
    particleStream->set_channel_map( channelMap );
    frantic::channels::channel_map map = particleStream->get_channel_map();
    frantic::channels::channel_accessor<vector3f> pAccessor = map.get_accessor<vector3f>( _T( "Position" ) );
    frantic::channels::channel_accessor<vector3f> nAccessor = map.get_accessor<vector3f>( _T( "Normal" ) );
    std::vector<char> buffer( map.structure_size() );

    int count = 0;
    vector3f surfaceNormal( 0, 0, -1 );
    // run the test: grab points, test if they are on the triangle, test if they have the expected normal.
    while( particleStream->get_particle( &buffer[0] ) ) {
        vector3f position = pAccessor.get( buffer );
        vector3f normal = nAccessor.get( buffer );

        EXPECT_EQ( true, frantic::graphics::contained_in_triangle( position, mesh->get_vertex( 0 ),
                                                                   mesh->get_vertex( 1 ), mesh->get_vertex( 2 ) ) );
        EXPECT_EQ( normal, surfaceNormal );
        count++;
    }

    EXPECT_EQ( numParticles, count );
}

TEST( SurfaceParticleIStream, ParticleSpacing ) {

    // set up the mesh
    boost::shared_ptr<trimesh3> mesh( new trimesh3() );
    mesh->add_vertex( 0, 0, 0 );
    mesh->add_vertex( 1, 1, 0 );
    mesh->add_vertex( 1, 0, 0 );
    mesh->add_face( 0, 1, 2 );

    float averageSpacing = 0.0025f;

    frantic::channels::channel_map channelMap;
    channelMap.define_channel( _T( "Position" ), 3, frantic::channels::data_type_float32 );
    channelMap.define_channel( _T( "Normal" ), 3, frantic::channels::data_type_float32 );
    channelMap.end_channel_definition();

    particle_istream_ptr particleStream = create_surface_particle_istream_using_spacing( *mesh, averageSpacing );
    particleStream->set_channel_map( channelMap );
    frantic::channels::channel_map map = particleStream->get_channel_map();
    frantic::channels::channel_accessor<vector3f> pAccessor = map.get_accessor<vector3f>( _T( "Position" ) );
    frantic::channels::channel_accessor<vector3f> nAccessor = map.get_accessor<vector3f>( _T( "Normal" ) );
    std::vector<char> buffer( map.structure_size() );

    vector3f surfaceNormal( 0, 0, -1 );
    int count = 0;
    // run the test: grab points, test if they are on the triangle, test if they have the expected normal.
    while( particleStream->get_particle( &buffer[0] ) ) {
        vector3f position = pAccessor.get( buffer );
        vector3f normal = nAccessor.get( buffer );

        EXPECT_EQ( true, frantic::graphics::contained_in_triangle( position, mesh->get_vertex( 0 ),
                                                                   mesh->get_vertex( 1 ), mesh->get_vertex( 2 ) ) );
        EXPECT_EQ( normal, surfaceNormal );
        count++;
    }

    // Get an estimate of the particle count.
    float surfaceArea =
        frantic::graphics::triangle_area( mesh->get_vertex( 0 ), mesh->get_vertex( 1 ), mesh->get_vertex( 2 ) );
    float particleEstimate = ( surfaceArea / static_cast<float>( pow( averageSpacing, 2 ) ) );

    EXPECT_NEAR( particleEstimate, count, 1 );
}
