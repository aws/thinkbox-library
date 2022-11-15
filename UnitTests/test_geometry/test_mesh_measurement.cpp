// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"

#include <frantic/geometry/mesh_measurement.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/surface_particle_istream.hpp>

namespace {

// assuming edge length 2.0
const float expectedSurfaceAreaPerTetrahedronFace = std::sqrt( 3.0f );

} // namespace

TEST( MeshMeasurement, SingleComponentSurfaceArea ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    EXPECT_FLOAT_EQ( 1.0f, get_mesh_surface_area( mesh ) );

    make_regular_tetrahedron( mesh );

    const float expectedSurfaceArea = expectedSurfaceAreaPerTetrahedronFace * 4;

    EXPECT_FLOAT_EQ( expectedSurfaceArea, get_mesh_surface_area( mesh ) );

    make_cube_mesh( mesh );

    EXPECT_FLOAT_EQ( 2.0f * 2.0f * 6.0f, get_mesh_surface_area( mesh ) );
}

TEST( MeshMeasurement, SingleComponentSurfaceAreaNonTriangulated ) {
    using namespace frantic::geometry;

    polymesh3_ptr polymesh;
    polymesh = make_cube_polymesh();
    std::unique_ptr<polymesh3_interface> polyInterface = polymesh3_interface::create_instance( polymesh );

    EXPECT_FLOAT_EQ( 2.0f * 2.0f * 6.0f, get_mesh_surface_area( polyInterface.get() ) );
}

TEST( MeshMeasurement, ApproxSphereSurfaceArea ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    make_sphere_mesh( 250, mesh );

    // radius is 1.0
    const float expectedSurfaceArea = float( 4.0 * M_PI ); // *1.0*1.0

    EXPECT_NEAR( expectedSurfaceArea, get_mesh_surface_area( mesh ), 1e-2f );
}

TEST( MeshMeasurement, FaceSubsetSurfaceArea ) {

    using namespace frantic::geometry;

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    EXPECT_FLOAT_EQ( 1.0f, get_mesh_surface_area( mesh ) );

    make_regular_tetrahedron( mesh );

    size_t component[2] = {
        0,
        3,
    };

    EXPECT_EQ( expectedSurfaceAreaPerTetrahedronFace * 2,
               get_mesh_component_surface_area( mesh, component, component + 2 ) );

    make_cube_mesh( mesh );

    size_t cubeComponents[4] = {
        0,
        1,
        4,
        5,
    };

    EXPECT_FLOAT_EQ( 2.0f * 2.0f * 2.0f, get_mesh_component_surface_area( mesh, cubeComponents, cubeComponents + 4 ) );
}

TEST( MeshMeasurement, SingleComponentVolume ) {
    using namespace frantic::geometry;

    trimesh3 mesh;

    make_regular_tetrahedron( mesh );

    // volume of a regular tetrahedron
    const float expectedVolume = float( ( 2.0 * 2.0 * 2.0 ) / ( 6.0 * sqrt( 2.0 ) ) );

    EXPECT_FLOAT_EQ( expectedVolume, tetrahedron_volume( mesh.get_vertex( 0 ), mesh.get_vertex( 1 ),
                                                         mesh.get_vertex( 2 ), mesh.get_vertex( 3 ) ) );
    EXPECT_FLOAT_EQ( expectedVolume, get_mesh_volume( mesh ) );

    make_cube_mesh( mesh );

    EXPECT_FLOAT_EQ( 2.0f * 2.0f * 2.0f, get_mesh_volume( mesh ) );
}

TEST( MeshMeasurement, SingleComponentVolumeNonTriangulated ) {
    using namespace frantic::geometry;

    polymesh3_ptr polymesh;
    polymesh = make_cube_polymesh();
    std::unique_ptr<polymesh3_interface> polyInterface = polymesh3_interface::create_instance( polymesh );

    EXPECT_FLOAT_EQ( 2.0f * 2.0f * 2.0f, get_mesh_volume( polyInterface.get() ) );
}

TEST( MeshMeasurement, ApproxSphereVolume ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    make_sphere_mesh( 250, mesh );

    // radius is 1.0
    const float expectedVolume = float( ( 4.0 / 3.0 ) * M_PI ); // *1.0*1.0*1.0

    EXPECT_NEAR( expectedVolume, get_mesh_volume( mesh ), 1e-2f );
}

TEST( MeshMeasurement, Hausdorff ) {

    using namespace frantic::geometry;
    using namespace frantic::channels;
    using namespace frantic::particles;
    using namespace frantic::particles::streams;

    boost::shared_ptr<trimesh3> meshA( new trimesh3() );
    boost::shared_ptr<trimesh3> meshB( new trimesh3() );
    boost::shared_ptr<particle_array> particlesA( new particle_array() );
    boost::shared_ptr<particle_array> particlesB( new particle_array() );

    bool faces = false;

    const bool vertices = true;
    const size_t numSamples = 150;

    // Test on empty meshes.
    double distance = hausdorff_distance_two_sided( *meshA, *meshB, vertices, faces, numSamples );
    EXPECT_EQ( 0.0, distance );

    // Test same mesh, distance should be zero when only using vertices.
    make_sphere_mesh( 5, *meshA );
    make_sphere_mesh( 5, *meshB );

    distance = hausdorff_distance_two_sided( *meshA, *meshB, vertices, faces, numSamples );
    EXPECT_EQ( 0.0, distance );

    // Test known difference, should be exactly one when not using the faces of the mesh.
    make_quad_trimesh( *meshA );
    make_quad_trimesh( *meshB );

    vector3f translate( 0, 0, 1 );
    meshA->translate( translate );

    distance = hausdorff_distance_two_sided( *meshA, *meshB, vertices, faces, numSamples );
    EXPECT_EQ( 1.0, distance );

    // When using faces, distance should be very close to the expected.
    faces = true;
    double error_buffer_faces = 0.0000001;
    distance = hausdorff_distance_two_sided( *meshB, *meshA, vertices, faces, numSamples );
    EXPECT_NEAR( 1.0, distance, error_buffer_faces );

    // Test that the result is reasonably close when using opposite ordering of meshes.
    make_sphere_mesh( 10, *meshA );
    make_sphere_mesh( 5, *meshB );

    double distance1 = hausdorff_distance_two_sided( *meshA, *meshB, vertices, faces, numSamples );
    double distance2 =
        hausdorff_distance_two_sided( *meshB, *meshA, vertices, faces, numSamples, particlesA, particlesB );
    double error_buffer = 0.005;
    EXPECT_NEAR( distance1, distance2, error_buffer );

    // Test that the particle streams that are setup in the hausdorff function are valid and properly filled.
    channel_map map = particlesA->get_channel_map();
    channel_accessor<vector3f> positions = map.get_accessor<vector3f>( _T( "Position" ) );
    channel_accessor<double> distances = map.get_accessor<double>( _T( "Distance" ) );
    size_t particleCount = particlesA->particle_count();

    EXPECT_TRUE( particleCount <= numSamples * 2 );

    for( int i = 0; i < particleCount; ++i ) {
        char* buffer = particlesA->at( i );
        EXPECT_NO_THROW( positions.get( buffer ) );
        EXPECT_NO_THROW( distances.get( buffer ) );
    }
}
