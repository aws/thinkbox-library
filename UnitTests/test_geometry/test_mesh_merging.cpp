// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/geometry/mesh_merging.hpp>

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <utilities/mesh_generators.hpp>

TEST( MeshMerging, WeldVerticesExact ) {
    using namespace frantic::geometry;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_vertices( meshInterfaces, 0 );
    polymesh3_ptr goal = make_welded_quad();

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, WeldVerticesTolerance ) {
    using namespace frantic::geometry;

    { // Past tolerance
        std::vector<polymesh3_ptr> meshes = make_separated_quad( 0.005f );

        const std::size_t meshCount = meshes.size();
        std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
        for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
            mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
            meshInterfaces[meshIndex].swap( meshInterface );
        }

        polymesh3_ptr result = weld_vertices( meshInterfaces, 0 );
        polymesh3_ptr goal = make_welded_quad();

        ASSERT_FALSE( is_consistent_topology( result, goal ) );
    }

    { // Within tolerance
        std::vector<polymesh3_ptr> meshes = make_separated_quad( 0.005f );

        const std::size_t meshCount = meshes.size();
        std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
        for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
            mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
            meshInterfaces[meshIndex].swap( meshInterface );
        }

        // Note: Different vertex welding selections could be correct depending on the implementation
        polymesh3_ptr result = weld_vertices( meshInterfaces, 0.01f );
        polymesh3_ptr goal = make_welded_quad();

        ASSERT_TRUE( is_consistent_topology( result, goal ) );
    }
}

TEST( MeshMerging, WeldChannelsBasic ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    { // First mesh channels
        float faceData[] = { 1 };
        meshes[0]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4 };
        meshes[0]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 1, 2, 3 };
        int customFaces[] = { 0, 1, 1, 2 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[0]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    { // Second mesh channels
        float faceData[] = { 2 };
        meshes[1]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 5, 6, 7, 8 };
        meshes[1]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 4, 5, 6 };
        int customFaces[] = { 2, 1, 1, 0 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[1]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_vertices( meshInterfaces, 0 );
    polymesh3_ptr goal = make_welded_quad();

    { // Goal mesh channels
        float faceData[] = { 1, 2 };
        goal->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4, 7, 8 };
        goal->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 1, 2, 3, 4, 5, 6 };
        int customFaces[] = { 0, 1, 1, 2, 5, 4, 4, 3 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 8 );
        goal->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    ASSERT_TRUE( is_equal( result, goal ) );
}

TEST( MeshMerging, WeldChannelsMissing ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    { // First mesh channels
        float faceData[] = { 1 };
        meshes[0]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4 };
        meshes[0]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 1, 2, 3 };
        int customFaces[] = { 0, 1, 1, 2 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[0]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_vertices( meshInterfaces, 0 );
    polymesh3_ptr goal = make_welded_quad();

    { // Goal mesh channels
        float faceData[] = { 1, 0 };
        goal->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4, 0, 0 };
        goal->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 0, 1, 2, 3 };
        int customFaces[] = { 1, 2, 2, 3, 0, 0, 0, 0 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 8 );
        goal->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    ASSERT_TRUE( is_equal( result, goal ) );
}

TEST( MeshMerging, WeldChannelsTypeConversion ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    { // First mesh channels
        double faceData[] = { 1 };
        meshes[0]->add_face_channel( _T("FaceChannel"), faceData );

        double vertexData[] = { 1, 2, 3, 4 };
        meshes[0]->add_vertex_channel( _T("VertexChannel"), vertexData );

        double customData[] = { 1, 2, 3 };
        int customFaces[] = { 0, 1, 1, 2 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[0]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    { // Second mesh channels
        float faceData[] = { 2 };
        meshes[1]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 5, 6, 7, 8 };
        meshes[1]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 4, 5, 6 };
        int customFaces[] = { 2, 1, 1, 0 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[1]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_vertices( meshInterfaces, 0 );
    polymesh3_ptr goal = make_welded_quad();

    { // Goal mesh channels
        double faceData[] = { 1, 2 };
        goal->add_face_channel( _T("FaceChannel"), faceData );

        double vertexData[] = { 1, 2, 3, 4, 7, 8 };
        goal->add_vertex_channel( _T("VertexChannel"), vertexData );

        double customData[] = { 1, 2, 3, 4, 5, 6 };
        int customFaces[] = { 0, 1, 1, 2, 5, 4, 4, 3 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 8 );
        goal->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    ASSERT_TRUE( is_equal( result, goal ) );
}

TEST( MeshMerging, WeldBoundariesExact ) {
    using namespace frantic::geometry;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0 );
    polymesh3_ptr goal = make_welded_quad();

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, WeldBoundariesTolerance ) {
    using namespace frantic::geometry;

    { // Past tolerance
        std::vector<polymesh3_ptr> meshes = make_separated_quad( 0.005f );

        const std::size_t meshCount = meshes.size();
        std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
        for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
            mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
            meshInterfaces[meshIndex].swap( meshInterface );
        }

        polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0 );
        polymesh3_ptr goal = make_welded_quad();

        ASSERT_FALSE( is_consistent_topology( result, goal ) );
    }

    { // Within tolerance
        std::vector<polymesh3_ptr> meshes = make_separated_quad( 0.005f );

        const std::size_t meshCount = meshes.size();
        std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
        for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
            mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
            meshInterfaces[meshIndex].swap( meshInterface );
        }

        // Note: Different vertex welding selections could be correct depending on the implementation
        polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0.01f );
        polymesh3_ptr goal = make_welded_quad();

        ASSERT_TRUE( is_consistent_topology( result, goal ) );
    }
}

/**
 * One of the criteria for two boundary edges to be welded is that each edge's head is within tolerance of the other's
 * tail. This means that, if they are using 'duplicate' vertices, the two edges must have opposite winding order. This
 * is to prevent welding edges from faces which point in opposite directions.
 */
TEST( MeshMerging, WeldBoundariesWindingOrder ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes( 2 );
    polymesh3_ptr goal;

    {
        polymesh3_builder firstBuilder;

        firstBuilder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( -1.0, 1.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );

        int face[] = { 0, 1, 2, 3 };
        firstBuilder.add_polygon( face, 4 );

        meshes[0] = firstBuilder.finalize();
    }

    {
        polymesh3_builder secondBuilder;

        secondBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );

        int face[] = { 0, 1, 2, 3 };
        secondBuilder.add_polygon( face, 4 );

        meshes[1] = secondBuilder.finalize();
    }

    {
        polymesh3_builder goalBuilder;

        goalBuilder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( -1.0, 1.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );

        int firstFace[] = { 0, 1, 2, 3 };
        goalBuilder.add_polygon( firstFace, 4 );

        int secondFace[] = { 4, 5, 6, 7 };
        goalBuilder.add_polygon( secondFace, 4 );

        goal = goalBuilder.finalize();
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0 );

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, WeldBoundariesIntraMesh ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes( 1 );
    polymesh3_ptr goal;

    {
        polymesh3_builder builder;

        builder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
        builder.add_vertex( vector3f( -1.0, 1.0, 0.0 ) );
        builder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        builder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        builder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        builder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        builder.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
        builder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

        int firstFace[] = { 0, 1, 2, 3 };
        builder.add_polygon( firstFace, 4 );

        int secondFace[] = { 4, 5, 6, 7 };
        builder.add_polygon( secondFace, 4 );

        goal = builder.finalize();
        meshes[0] = goal;
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0 );

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, WeldBoundariesNonManifold ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes( 2 );
    polymesh3_ptr goal;

    {
        polymesh3_builder firstBuilder;

        firstBuilder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, 0.0, 1.0 ) );

        int firstFace[] = { 0, 1, 2 };
        firstBuilder.add_polygon( firstFace, 3 );

        int secondFace[] = { 2, 3, 1 };
        firstBuilder.add_polygon( secondFace, 3 );

        meshes[0] = firstBuilder.finalize();
    }

    {
        polymesh3_builder secondBuilder;

        secondBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );

        int face[] = { 0, 1, 2 };
        secondBuilder.add_polygon( face, 3 );

        meshes[1] = secondBuilder.finalize();
    }

    {
        polymesh3_builder goalBuilder;

        goalBuilder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 0.0, 1.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );

        int firstFace[] = { 0, 1, 2 };
        goalBuilder.add_polygon( firstFace, 3 );

        int secondFace[] = { 2, 3, 1 };
        goalBuilder.add_polygon( secondFace, 3 );

        int thirdFace[] = { 4, 5, 6 };
        goalBuilder.add_polygon( thirdFace, 3 );

        goal = goalBuilder.finalize();
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0 );

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, WeldBoundariesAdjacentEdges ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes( 2 );
    polymesh3_ptr goal;

    {
        polymesh3_builder firstBuilder;

        firstBuilder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        firstBuilder.add_vertex( vector3f( 0.0, -1.0, 0.0 ) );

        int firstFace[] = { 0, 1, 2 };
        firstBuilder.add_polygon( firstFace, 3 );

        int secondFace[] = { 0, 2, 3 };
        firstBuilder.add_polygon( secondFace, 3 );

        meshes[0] = firstBuilder.finalize();
    }

    {
        polymesh3_builder secondBuilder;

        secondBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );
        secondBuilder.add_vertex( vector3f( 0.0, -1.0, 0.0 ) );

        int firstFace[] = { 0, 1, 2 };
        secondBuilder.add_polygon( firstFace, 3 );

        int secondFace[] = { 0, 2, 3 };
        secondBuilder.add_polygon( secondFace, 3 );

        meshes[1] = secondBuilder.finalize();
    }

    {
        polymesh3_builder goalBuilder;

        goalBuilder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 0.0, -1.0, 0.0 ) );
        goalBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

        int firstFace[] = { 0, 1, 2 };
        goalBuilder.add_polygon( firstFace, 3 );

        int secondFace[] = { 0, 2, 3 };
        goalBuilder.add_polygon( secondFace, 3 );

        int thirdFace[] = { 2, 1, 4 };
        goalBuilder.add_polygon( thirdFace, 3 );

        int fourthFace[] = { 2, 4, 3 };
        goalBuilder.add_polygon( fourthFace, 3 );

        goal = goalBuilder.finalize();
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0 );

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, WeldBoundariesEmpty ) {
    using namespace frantic::geometry;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    polymesh3_builder emptyBuilder;
    meshes.push_back( emptyBuilder.finalize() );

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    meshInterfaces.push_back( mesh_interface_ptr() );

    polymesh3_ptr result = weld_boundary_edges( meshInterfaces, 0 );
    polymesh3_ptr goal = make_welded_quad();

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, CombineBasic ) {
    using namespace frantic::geometry;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = combine( meshInterfaces );
    polymesh3_ptr goal = make_combined_quad();

    ASSERT_TRUE( is_consistent_topology( result, goal ) );
}

TEST( MeshMerging, CombineChannelsBasic ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    { // First mesh channels
        float faceData[] = { 1 };
        meshes[0]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4 };
        meshes[0]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 1, 2, 3 };
        int customFaces[] = { 0, 1, 1, 2 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[0]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    { // Second mesh channels
        float faceData[] = { 2 };
        meshes[1]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 5, 6, 7, 8 };
        meshes[1]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 4, 5, 6 };
        int customFaces[] = { 2, 1, 1, 0 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[1]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = combine( meshInterfaces );
    polymesh3_ptr goal = make_combined_quad();

    { // Goal mesh channels
        float faceData[] = { 1, 2 };
        goal->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
        goal->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 1, 2, 3, 4, 5, 6 };
        int customFaces[] = { 0, 1, 1, 2, 5, 4, 4, 3 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 8 );
        goal->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    ASSERT_TRUE( is_equal( result, goal ) );
}

TEST( MeshMerging, CombineChannelsMissing ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    { // First mesh channels
        float faceData[] = { 1 };
        meshes[0]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4 };
        meshes[0]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 1, 2, 3 };
        int customFaces[] = { 0, 1, 1, 2 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[0]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = combine( meshInterfaces );
    polymesh3_ptr goal = make_combined_quad();

    { // Goal mesh channels
        float faceData[] = { 1, 0 };
        goal->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 1, 2, 3, 4, 0, 0, 0, 0 };
        goal->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 0, 1, 2, 3 };
        int customFaces[] = { 1, 2, 2, 3, 0, 0, 0, 0 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 8 );
        goal->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    ASSERT_TRUE( is_equal( result, goal ) );
}

TEST( MeshMerging, CombineChannelsTypeConversion ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<polymesh3_ptr> meshes = make_separated_quad( 0 );

    { // First mesh channels
        double faceData[] = { 1 };
        meshes[0]->add_face_channel( _T("FaceChannel"), faceData );

        double vertexData[] = { 1, 2, 3, 4 };
        meshes[0]->add_vertex_channel( _T("VertexChannel"), vertexData );

        double customData[] = { 1, 2, 3 };
        int customFaces[] = { 0, 1, 1, 2 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[0]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    { // Second mesh channels
        float faceData[] = { 2 };
        meshes[1]->add_face_channel( _T("FaceChannel"), faceData );

        float vertexData[] = { 5, 6, 7, 8 };
        meshes[1]->add_vertex_channel( _T("VertexChannel"), vertexData );

        float customData[] = { 4, 5, 6 };
        int customFaces[] = { 2, 1, 1, 0 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 4 );
        meshes[1]->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( meshes[meshIndex] ).release() );
        meshInterfaces[meshIndex].swap( meshInterface );
    }

    polymesh3_ptr result = combine( meshInterfaces );
    polymesh3_ptr goal = make_combined_quad();

    { // Goal mesh channels
        double faceData[] = { 1, 2 };
        goal->add_face_channel( _T("FaceChannel"), faceData );

        double vertexData[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
        goal->add_vertex_channel( _T("VertexChannel"), vertexData );

        double customData[] = { 1, 2, 3, 4, 5, 6 };
        int customFaces[] = { 0, 1, 1, 2, 5, 4, 4, 3 };
        std::vector<int> customFacesBuffer( &customFaces[0], &customFaces[0] + 8 );
        goal->add_vertex_channel( _T("CustomChannel"), customData, &customFacesBuffer );
    }

    ASSERT_TRUE( is_equal( result, goal ) );
}
