// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>

#include <boost/range/irange.hpp>

namespace {

void check_adjacent_tetrahedron_faces( const frantic::graphics::vector3& lhs, const frantic::graphics::vector3& rhs,
                                       const size_t endpoints[2] ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    vector3 leftFace = lhs.to_sorted();
    vector3 rightFace = rhs.to_sorted();

    size_t sortedEndpoints[2] = { endpoints[0], endpoints[1] };
    std::sort( &sortedEndpoints[0], &sortedEndpoints[2] );

    size_t leftDiff[3];
    size_t* leftDiffEnd =
        std::set_difference( &leftFace[0], &leftFace[3], &sortedEndpoints[0], &sortedEndpoints[2], leftDiff );
    EXPECT_EQ( 1, leftDiffEnd - leftDiff );

    size_t rightDiff[3];
    size_t* rightDiffEnd =
        std::set_difference( &rightFace[0], &rightFace[3], &sortedEndpoints[0], &sortedEndpoints[2], rightDiff );
    EXPECT_EQ( 1, rightDiffEnd - rightDiff );

    EXPECT_NE( leftDiff[0], rightDiff[0] );
}

frantic::graphics::vector3 get_triangle_face( const frantic::geometry::mesh_interface_ptr meshPtr, size_t faceId ) {
    return frantic::graphics::vector3( (int)meshPtr->get_face_vert_index( (int)faceId, 0 ),
                                       (int)meshPtr->get_face_vert_index( (int)faceId, 1 ),
                                       (int)meshPtr->get_face_vert_index( (int)faceId, 2 ) );
}

void check_tetrahedron_adjacency( const frantic::geometry::mesh_interface_ptr meshPtr ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    meshPtr->init_adjacency();

    for( size_t vertexId = 0; vertexId < 4; ++vertexId ) {

        vertex_iterator vIt;
        EXPECT_TRUE( meshPtr->init_vertex_iterator( vIt, vertexId ) );

        std::vector<size_t> incidentVerts;

        do {
            incidentVerts.push_back( meshPtr->get_edge_endpoint( vIt ) );

            vector3 leftFace = get_triangle_face( meshPtr, meshPtr->get_edge_left_face( vIt ) ).to_sorted();
            vector3 rightFace = get_triangle_face( meshPtr, meshPtr->get_edge_right_face( vIt ) ).to_sorted();
            size_t endPoints[2] = { vertexId, incidentVerts.back() };

            check_adjacent_tetrahedron_faces( leftFace, rightFace, endPoints );

        } while( meshPtr->advance_vertex_iterator( vIt ) );

        EXPECT_EQ( 3, incidentVerts.size() );

        std::sort( incidentVerts.begin(), incidentVerts.end() );

        size_t vertexDiff[4];
        size_t* vertexDiffEnd = std::set_difference( boost::irange( 0, 4 ).begin(), boost::irange( 0, 4 ).end(),
                                                     incidentVerts.begin(), incidentVerts.end(), vertexDiff );

        EXPECT_EQ( 1, vertexDiffEnd - vertexDiff );
        EXPECT_EQ( vertexId, vertexDiff[0] );
    }

    for( size_t faceId = 0; faceId < 4; ++faceId ) {

        face_iterator fIt;

        meshPtr->init_face_iterator( fIt, faceId );

        std::vector<size_t> incidentFaces;

        vector3 thisFace = get_triangle_face( meshPtr, faceId ).to_sorted();

        do {
            incidentFaces.push_back( meshPtr->get_face_neighbor( fIt ) );

            vector3 otherFace = get_triangle_face( meshPtr, incidentFaces.back() ).to_sorted();
            size_t endPoints[2] = { meshPtr->get_face_next_vertex( fIt ), meshPtr->get_face_prev_vertex( fIt ) };

            check_adjacent_tetrahedron_faces( thisFace, otherFace, endPoints );

        } while( meshPtr->advance_face_iterator( fIt ) );

        EXPECT_EQ( 3, incidentFaces.size() );

        std::sort( incidentFaces.begin(), incidentFaces.end() );

        // TODO: Should be asserting something about this...
        // size_t faceDiff[4];
        // size_t* faceDiffEnd = std::set_difference( boost::irange( 0, 4 ).begin(), boost::irange( 0, 4 ).end(),
        //                                           incidentFaces.begin(), incidentFaces.end(), faceDiff );
    }
}

} // namespace

TEST( MeshInterface, TestTrimesh3InterfaceAdjacency ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    trimesh3 baseMesh;
    make_regular_tetrahedron( baseMesh );
    mesh_interface_ptr baseMeshInterface( trimesh3_interface::create_instance( &baseMesh ).release() );
    check_tetrahedron_adjacency( baseMeshInterface );
}

TEST( MeshInterface, TestPolymesh3InterfaceAdjacency ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( make_regular_tetrahedron() ).release() );
    check_tetrahedron_adjacency( meshInterface );
}

// Test read/write of elements when converting to an arity lower than the source data
TEST( MeshInterface, MeshChannelArityAdjustmentDown ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;
    using namespace frantic::channels;

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( make_regular_tetrahedron() ).release() );

    const frantic::tstring simpleChannelName( _T( "Simple" ) );

    mesh_channel* simpleChannel =
        meshInterface->add_vertex_channel( simpleChannelName, frantic::channels::data_type_float32, 3 );

    ASSERT_NO_THROW( mesh_channel_cvt<frantic::graphics2d::vector2f> testNoThrow( simpleChannel ) );

    mesh_channel_cvt<frantic::graphics::vector3f> testChannel3( simpleChannel );
    for( size_t i = 0; i < testChannel3.get_num_elements(); ++i ) {
        testChannel3.set_value( i, frantic::graphics::vector3f( -1.0f, -1.0f, -1.0f ) );
    }

    mesh_channel_cvt<frantic::graphics2d::vector2f> testChannel2( simpleChannel );

    frantic::graphics2d::vector2f outValue2;
    frantic::graphics::vector3f outValue3;

    outValue2 = testChannel2.get_value( 0 );
    EXPECT_FLOAT_EQ( -1.0f, outValue2.x );
    EXPECT_FLOAT_EQ( -1.0f, outValue2.y );

    // test setting/getting in the reduced arity
    testChannel2.set_value( 0, frantic::graphics2d::vector2f( 1.0f, 2.0f ) );
    outValue2 = testChannel2.get_value( 0 );
    EXPECT_FLOAT_EQ( 1.0f, outValue2.x );
    EXPECT_FLOAT_EQ( 2.0f, outValue2.y );

    // test that the 3rd value was truncated in the larger arity
    outValue3 = testChannel3.get_value( 0 );
    EXPECT_FLOAT_EQ( 1.0f, outValue3.x );
    EXPECT_FLOAT_EQ( 2.0f, outValue3.y );
    EXPECT_FLOAT_EQ( 0.0f, outValue3.z );
}
