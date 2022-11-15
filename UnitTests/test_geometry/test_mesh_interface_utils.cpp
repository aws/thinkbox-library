// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>

#include <frantic/particles/particle_array.hpp>

#include "utilities/mesh_generators.hpp"

TEST( MeshInterfaceUtils, ComputeBoundboxEmpty ) {
    frantic::geometry::trimesh3 mesh;

    frantic::geometry::mesh_interface::ptr_type i(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    frantic::geometry::boundbox3f bounds = frantic::geometry::compute_boundbox( i );

    ASSERT_TRUE( bounds.is_empty() );
}

TEST( MeshInterfaceUtils, ComputeBoundboxUnitCube ) {
    frantic::geometry::trimesh3 mesh;
    mesh.set_to_box( frantic::graphics::boundbox3f( 0, 1, 0, 1, 0, 1 ) );

    frantic::geometry::mesh_interface::ptr_type i(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    frantic::geometry::boundbox3f bounds = frantic::geometry::compute_boundbox( i );

    ASSERT_EQ( bounds, frantic::graphics::boundbox3f( 0, 1, 0, 1, 0, 1 ) );
}

TEST( MeshInterfaceUtils, GetFaceVertexCount ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 0 );
    mesh.add_vertex( 0, 1, 0 );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 2, 3 );

    frantic::geometry::mesh_interface::ptr_type m(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    EXPECT_EQ( frantic::geometry::get_face_vertex_count( m.get() ), 6 );
}

TEST( MeshInterfaceUtils, IsClosedManifoldOneTriangle ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 0 );

    mesh.add_face( 0, 1, 2 );

    frantic::geometry::mesh_interface::ptr_type i(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    EXPECT_FALSE( frantic::geometry::is_closed_manifold( i ) );
}

TEST( MeshInterfaceUtils, IsClosedManifoldUnitCube ) {
    frantic::geometry::trimesh3 mesh;
    mesh.set_to_box( frantic::graphics::boundbox3f( 0, 1, 0, 1, 0, 1 ) );

    frantic::geometry::mesh_interface::ptr_type i(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    EXPECT_TRUE( frantic::geometry::is_closed_manifold( i ) );
}

TEST( MeshInterfaceUtils, IsClosedManifoldOneSingularEdge ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 0, 0, 1 );
    mesh.add_vertex( -1, 0, 0 );
    mesh.add_vertex( 0, 1, 0 );
    mesh.add_vertex( 1, 0, 0 );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 1, 3 );
    mesh.add_face( 0, 1, 4 );

    frantic::geometry::mesh_interface::ptr_type i(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    EXPECT_FALSE( frantic::geometry::is_closed_manifold( i ) );
}

TEST( MeshInterfaceUtils, HasDegenerateFaceNormalTriangle ) {
    using namespace frantic::geometry;

    polymesh3_ptr mesh = make_triangle_polymesh();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );

    EXPECT_FALSE( has_degenerate_faces( meshInterface ) );
}

TEST( MeshInterfaceUtils, HasDegenerateFaceSliverTriangle ) {
    using namespace frantic::geometry;

    polymesh3_ptr mesh = make_sliver_polymesh();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );

    EXPECT_TRUE( has_degenerate_faces( meshInterface ) );
}

TEST( MeshInterfaceUtils, FixDegenerateFaceSliverTriangle ) {
    using namespace frantic::geometry;

    polymesh3_ptr mesh = make_sliver_polymesh();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    polymesh3_ptr fixed = fix_degenerate_faces( meshInterface );

    EXPECT_EQ( 0, fixed->face_count() );
}

TEST( MeshInterfaceUtils, FixDegenerateFaceDuplicateCorner ) {
    using namespace frantic::geometry;

    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, 1, 0 );
    const int points[] = { 0, 1, 1, 2 };
    builder.add_polygon( points, 4 );
    polymesh3_ptr mesh = builder.finalize();

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    polymesh3_ptr fixed = fix_degenerate_faces( meshInterface );

    ASSERT_EQ( 1, fixed->face_count() );
    polymesh3_const_vertex_accessor<void> geomAcc = fixed->get_const_vertex_accessor( _T("verts") );

    polymesh3_const_vertex_accessor<void>::const_face_range face = geomAcc.get_face( 0 );
    ASSERT_EQ( 3, face.second - face.first );
    EXPECT_EQ( 0, face.first[0] );
    EXPECT_EQ( 1, face.first[1] );
    EXPECT_EQ( 2, face.first[2] );
}

TEST( MeshInterfaceUtils, FixDegenerateFacePinchedPolygon ) {
    using namespace frantic::geometry;

    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, 0, 2 );
    builder.add_vertex( 0, 1, 2 );
    builder.add_vertex( 0, 1, 0 );
    const int points[] = { 0, 1, 2, 3, 1, 4 };
    builder.add_polygon( points, 6 );
    polymesh3_ptr mesh = builder.finalize();

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    polymesh3_ptr fixed = fix_degenerate_faces( meshInterface );

    ASSERT_EQ( 2, fixed->face_count() );
    polymesh3_const_vertex_accessor<void> geomAcc = fixed->get_const_vertex_accessor( _T("verts") );

    polymesh3_const_vertex_accessor<void>::const_face_range face = geomAcc.get_face( 0 );
    ASSERT_EQ( 3, face.second - face.first );
    EXPECT_EQ( 1, face.first[0] );
    EXPECT_EQ( 2, face.first[1] );
    EXPECT_EQ( 3, face.first[2] );

    face = geomAcc.get_face( 1 );
    ASSERT_EQ( 3, face.second - face.first );
    EXPECT_EQ( 0, face.first[0] );
    EXPECT_EQ( 1, face.first[1] );
    EXPECT_EQ( 4, face.first[2] );
}

TEST( MeshInterfaceUtils, FixDegenerateFacePinchedSliver ) {
    using namespace frantic::geometry;

    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, 0, 2 );
    builder.add_vertex( 0, 1, 0 );
    const int points[] = { 0, 1, 2, 1, 3 };
    builder.add_polygon( points, 5 );
    polymesh3_ptr mesh = builder.finalize();

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    polymesh3_ptr fixed = fix_degenerate_faces( meshInterface );

    ASSERT_EQ( 1, fixed->face_count() );
    polymesh3_const_vertex_accessor<void> geomAcc = fixed->get_const_vertex_accessor( _T("verts") );

    polymesh3_const_vertex_accessor<void>::const_face_range face = geomAcc.get_face( 0 );
    ASSERT_EQ( 3, face.second - face.first );
    EXPECT_EQ( 0, face.first[0] );
    EXPECT_EQ( 1, face.first[1] );
    EXPECT_EQ( 3, face.first[2] );
}

TEST( MeshInterfaceUtils, FixDegenerateFaceVertexChannel ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, 0, 2 );
    builder.add_vertex( 0, 1, 2 );
    builder.add_vertex( 0, 1, 0 );
    const int points[] = { 0, 1, 2, 3, 1, 4 };
    builder.add_polygon( points, 6 );
    polymesh3_ptr mesh = builder.finalize();

    float data[] = { 5, 6, 7, 8 };
    frantic::graphics::raw_byte_buffer buffer( data, sizeof( data ) );

    std::vector<int> polygon;
    polygon.push_back( 2 );
    polygon.push_back( 1 );
    polygon.push_back( 0 );
    polygon.push_back( 1 );
    polygon.push_back( 3 );
    polygon.push_back( 2 );

    mesh->add_vertex_channel( _T("Channel"), data_type_float32, 1, buffer, &polygon );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    polymesh3_ptr fixed = fix_degenerate_faces( meshInterface );

    polymesh3_const_vertex_accessor<void> chanAcc = fixed->get_const_vertex_accessor( _T("Channel") );
    ASSERT_EQ( 2, chanAcc.face_count() );

    polymesh3_const_vertex_accessor<void>::const_face_range face = chanAcc.get_face( 0 );
    ASSERT_EQ( 3, face.second - face.first );
    EXPECT_EQ( 1, face.first[0] );
    EXPECT_EQ( 0, face.first[1] );
    EXPECT_EQ( 1, face.first[2] );

    face = chanAcc.get_face( 1 );
    ASSERT_EQ( 3, face.second - face.first );
    EXPECT_EQ( 2, face.first[0] );
    EXPECT_EQ( 3, face.first[1] );
    EXPECT_EQ( 2, face.first[2] );
}

TEST( MeshInterfaceUtils, FixDegenerateFaceFaceChannel ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, 0, 2 );
    builder.add_vertex( 0, 1, 2 );
    builder.add_vertex( 0, 1, 0 );
    const int points[] = { 0, 1, 2, 3, 1, 4 };
    builder.add_polygon( points, 6 );
    polymesh3_ptr mesh = builder.finalize();

    float data[] = { 5 };
    frantic::graphics::raw_byte_buffer buffer( data, sizeof( data ) );
    mesh->add_face_channel( _T("Channel"), data_type_float32, 1, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    polymesh3_ptr fixed = fix_degenerate_faces( meshInterface );

    polymesh3_face_accessor<void> chanAcc = fixed->get_face_accessor( _T("Channel") );
    ASSERT_EQ( 2, chanAcc.face_count() );

    const float* face = reinterpret_cast<float*>( chanAcc.get_face( 0 ) );
    EXPECT_EQ( 5, *face );

    face = reinterpret_cast<float*>( chanAcc.get_face( 1 ) );
    EXPECT_EQ( 5, *face );
}

TEST( MeshInterfaceUtils, FixDegenerateFaceMultiplePolygons ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder builder;
    builder.add_vertex( 0, -2, 1 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, -1, 0 );
    builder.add_vertex( 0, -2, 0 );
    builder.add_vertex( 0, 2, 1 );
    builder.add_vertex( 0, 2, 0 );
    builder.add_vertex( 0, 1, 0 );

    const int first[4] = { 0, 1, 2, 3 };
    builder.add_polygon( first, 4 );
    const int second[4] = { 2, 1, 1, 6 };
    builder.add_polygon( second, 4 );
    const int third[4] = { 1, 4, 5, 6 };
    builder.add_polygon( third, 4 );

    polymesh3_ptr mesh = builder.finalize();

    float data[] = { 5, 6, 7 };
    frantic::graphics::raw_byte_buffer buffer( data, sizeof( data ) );
    mesh->add_face_channel( _T("Channel"), data_type_float32, 1, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    polymesh3_ptr fixed = fix_degenerate_faces( meshInterface );

    ASSERT_EQ( 3, fixed->face_count() );
    polymesh3_const_vertex_accessor<void> geomAcc = fixed->get_const_vertex_accessor( _T("verts") );

    polymesh3_const_vertex_accessor<void>::const_face_range face = geomAcc.get_face( 0 );
    ASSERT_EQ( 4, face.second - face.first );
    EXPECT_EQ( 0, face.first[0] );
    EXPECT_EQ( 1, face.first[1] );
    EXPECT_EQ( 2, face.first[2] );
    EXPECT_EQ( 3, face.first[3] );

    face = geomAcc.get_face( 1 );
    ASSERT_EQ( 3, face.second - face.first );
    EXPECT_EQ( 2, face.first[0] );
    EXPECT_EQ( 1, face.first[1] );
    EXPECT_EQ( 6, face.first[2] );

    face = geomAcc.get_face( 2 );
    ASSERT_EQ( 4, face.second - face.first );
    EXPECT_EQ( 1, face.first[0] );
    EXPECT_EQ( 4, face.first[1] );
    EXPECT_EQ( 5, face.first[2] );
    EXPECT_EQ( 6, face.first[3] );

    polymesh3_face_accessor<void> chanAcc = fixed->get_face_accessor( _T("Channel") );
    ASSERT_EQ( 3, chanAcc.face_count() );

    const float* chanFace = reinterpret_cast<float*>( chanAcc.get_face( 0 ) );
    EXPECT_EQ( 5, *chanFace );

    chanFace = reinterpret_cast<float*>( chanAcc.get_face( 1 ) );
    EXPECT_EQ( 6, *chanFace );

    chanFace = reinterpret_cast<float*>( chanAcc.get_face( 2 ) );
    EXPECT_EQ( 7, *chanFace );
}

TEST( MeshInterfaceUtils, ComputeVertexNormals ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 0 );

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 0, 1, 1 );
    mesh.add_vertex( 0, 0, 1 );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 3, 4, 5 );

    frantic::geometry::mesh_interface::ptr_type i(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    std::vector<frantic::graphics::vector3f> vertexNormals;

    frantic::geometry::compute_vertex_normals( i.get(), vertexNormals );

    ASSERT_EQ( vertexNormals.size(), 6 );

    EXPECT_EQ( vertexNormals[0], frantic::graphics::vector3f( 0, 0, 1 ) );
    EXPECT_EQ( vertexNormals[1], frantic::graphics::vector3f( 0, 0, 1 ) );
    EXPECT_EQ( vertexNormals[2], frantic::graphics::vector3f( 0, 0, 1 ) );

    EXPECT_EQ( vertexNormals[3], frantic::graphics::vector3f( 1, 0, 0 ) );
    EXPECT_EQ( vertexNormals[4], frantic::graphics::vector3f( 1, 0, 0 ) );
    EXPECT_EQ( vertexNormals[5], frantic::graphics::vector3f( 1, 0, 0 ) );
}

TEST( MeshInterfaceUtils, CreateVertexNormalChannelMeshInterface ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 0 );

    mesh.add_face( 0, 1, 2 );

    // add a Color channel, which should be propagated to the output mesh
    mesh.add_vertex_channel<frantic::graphics::vector3f>( _T("Color") );
    {
        frantic::geometry::trimesh3_vertex_channel_accessor<frantic::graphics::vector3f> acc(
            mesh.get_vertex_channel_accessor<frantic::graphics::vector3f>( _T("Color") ) );
        acc[0].set( 1, 0, 0 );
        acc[1].set( 0, 1, 0 );
        acc[2].set( 0, 0, 1 );
    }

    // add a FaceSelection channel, which should be propagated to the output mesh
    mesh.add_face_channel<boost::int32_t>( _T("FaceSelection") );
    {
        frantic::geometry::trimesh3_face_channel_accessor<boost::int32_t> acc(
            mesh.get_face_channel_accessor<boost::int32_t>( _T("FaceSelection") ) );
        acc[0] = 1;
    }

    // we'll generate normals in this meshCopy later, for the sake of comparison
    frantic::geometry::trimesh3 meshCopy( mesh );

    frantic::geometry::mesh_interface::ptr_type mBeforeNormal(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    frantic::geometry::mesh_interface::ptr_type m(
        frantic::geometry::create_vertex_normal_channel_mesh_interface( mBeforeNormal.get() ) );

    ASSERT_TRUE( bool( m ) );
    ASSERT_TRUE( m->is_valid() );

    ASSERT_TRUE( m->get_vertex_channels().has_channel( _T("Color") ) );
    ASSERT_TRUE( m->get_vertex_channels().has_channel( _T("Normal") ) );
    ASSERT_TRUE( m->get_face_channels().has_channel( _T("FaceSelection") ) );

    meshCopy.build_vertex_normals();

    frantic::geometry::mesh_interface::ptr_type mCopy(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( meshCopy ) ).release() );

    ASSERT_TRUE( frantic::geometry::is_equal( m, mCopy ) );
}

TEST( MeshInterfaceUtils, CreateVertexNormalChannelMeshInterfaceReplaceVertexChannel ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 0 );

    mesh.add_face( 0, 1, 2 );

    mesh.add_vertex_channel<int>( _T("Normal") );
    {
        frantic::geometry::trimesh3_vertex_channel_accessor<int> acc(
            mesh.get_vertex_channel_accessor<int>( _T("Normal") ) );
        for( std::size_t i = 0, ie = acc.size(); i < ie; ++i ) {
            acc[i] = -1;
        }
    }

    frantic::geometry::mesh_interface::ptr_type mBeforeNormal(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );
    frantic::geometry::mesh_interface::ptr_type m(
        frantic::geometry::create_vertex_normal_channel_mesh_interface( mBeforeNormal.get() ) );

    ASSERT_TRUE( bool( m ) );
    ASSERT_TRUE( m->is_valid() );

    const frantic::geometry::mesh_channel* normalChannel = m->get_vertex_channels().get_channel( _T("Normal" ) );

    ASSERT_TRUE( normalChannel != 0 );

    ASSERT_EQ( normalChannel->get_num_faces(), 1 );
    ASSERT_EQ( normalChannel->get_num_face_verts( 0 ), 3 );

    ASSERT_TRUE( frantic::geometry::mesh_channel::is_stored_at_vertex( normalChannel->get_channel_type() ) );

    ASSERT_EQ( normalChannel->get_data_type(), frantic::channels::data_type_float32 );
    ASSERT_EQ( normalChannel->get_data_arity(), 3 );

    frantic::geometry::mesh_channel_cvt<frantic::graphics::vector3f> acc( normalChannel );

    EXPECT_EQ( acc.get_value( acc.get_fv_index( 0, 0 ) ), frantic::graphics::vector3f( 0, 0, 1 ) );
    EXPECT_EQ( acc.get_value( acc.get_fv_index( 0, 1 ) ), frantic::graphics::vector3f( 0, 0, 1 ) );
    EXPECT_EQ( acc.get_value( acc.get_fv_index( 0, 2 ) ), frantic::graphics::vector3f( 0, 0, 1 ) );
}

TEST( MeshInterfaceUtils, CreateParticleArrayMeshInterface ) {
    using frantic::graphics::vector3f;

    std::unique_ptr<frantic::geometry::mesh_interface> mesh;

    { // scope for particles
        frantic::channels::channel_map channelMap;
        channelMap.define_channel<vector3f>( _T("Position") );
        channelMap.define_channel<vector3f>( _T("Color") );
        channelMap.end_channel_definition();

        frantic::particles::particle_array particles( channelMap );

        frantic::channels::channel_accessor<vector3f> positionAcc(
            particles.get_channel_map().get_accessor<vector3f>( _T("Position") ) );
        frantic::channels::channel_accessor<vector3f> colorAcc(
            particles.get_channel_map().get_accessor<vector3f>( _T("Color") ) );

        std::vector<char> buffer( particles.get_channel_map().structure_size() );

        positionAcc( buffer ).set( 0, 1, 2 );
        colorAcc( buffer ).set( 3, 4, 5 );

        particles.push_back( &buffer[0] );

        mesh = frantic::geometry::create_particle_array_mesh_interface( boost::move( particles ) );
    }

    ASSERT_TRUE( mesh.get() != 0 );
    ASSERT_TRUE( mesh->is_valid() );

    ASSERT_EQ( mesh->get_num_verts(), 1 );
    EXPECT_EQ( mesh->get_num_faces(), 0 );

    EXPECT_EQ( mesh->get_vert( 0 ), vector3f( 0, 1, 2 ) );

    EXPECT_FALSE( mesh->has_vertex_channel( _T("Position") ) );
    EXPECT_TRUE( mesh->has_vertex_channel( _T("Color") ) );

    const frantic::geometry::mesh_interface::mesh_channel_map& channels = mesh->get_vertex_channels();

    const frantic::geometry::mesh_channel* colorChannel = channels.get_channel( _T("Color") );

    ASSERT_TRUE( colorChannel != 0 );
    ASSERT_EQ( colorChannel->get_data_type(), frantic::channels::data_type_float32 );
    ASSERT_EQ( colorChannel->get_data_arity(), 3 );
    ASSERT_EQ( colorChannel->get_num_elements(), 1 );
    EXPECT_EQ( colorChannel->get_name(), _T("Color") );

    vector3f color;
    colorChannel->get_value( 0, &color[0] );
    EXPECT_EQ( color, vector3f( 3, 4, 5 ) );
}

TEST( MeshInterfaceUtils, FindFaceCorner ) {
    using namespace frantic::geometry;

    trimesh3 mesh;

    make_regular_tetrahedron( mesh );

    mesh_interface_ptr iMesh( trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    for( size_t faceId = 0; faceId < iMesh->get_num_faces(); ++faceId ) {
        for( size_t faceVertex = 0; faceVertex < 3; ++faceVertex ) {
            EXPECT_EQ( faceVertex,
                       find_face_corner( iMesh, faceId, iMesh->get_face_vert_index( faceId, faceVertex ) ) );
        }

        EXPECT_ANY_THROW( find_face_corner( iMesh, faceId, faceId ) );
    }
}

TEST( MeshInterfaceUtils, IsTriangleMesh ) {
    using namespace frantic::geometry;

    trimesh3 trimesh;
    make_quad_trimesh( trimesh );
    mesh_interface::ptr_type tri(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( trimesh ) ).release() );

    polymesh3_ptr polymesh = make_quad_polymesh();
    mesh_interface::ptr_type poly( frantic::geometry::polymesh3_interface::create_instance( polymesh ) );

    EXPECT_TRUE( is_triangle_mesh( tri.get() ) );
    EXPECT_FALSE( is_triangle_mesh( poly.get() ) );
}

TEST( MeshInterfaceUtils, AssertValidIndices ) {
    using namespace frantic::geometry;

    trimesh3 original;
    original.add_vertex( 0, 0, 0 );
    original.add_vertex( 1, 0, 0 );
    original.add_vertex( 1, 1, 0 );
    original.add_face( 0, 1, 2 );

    original.add_vertex_channel<float>( _T("Channel"), 1, true );

    EXPECT_NO_THROW( assert_valid_indices( trimesh3_interface::create_instance( &original ).get() ) );

    trimesh3 badFaceIndex( original );
    badFaceIndex.get_face( 0 )[0] = 3;

    EXPECT_ANY_THROW( assert_valid_indices( trimesh3_interface::create_instance( &badFaceIndex ).get() ) );

    trimesh3 badChannelIndex( original );
    trimesh3_vertex_channel_accessor<float> acc = badChannelIndex.get_vertex_channel_accessor<float>( _T("Channel") );
    acc.face( 0 )[2] = 1;

    EXPECT_ANY_THROW( assert_valid_indices( trimesh3_interface::create_instance( &badChannelIndex ).get() ) );
}
