// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/geometry/trimesh3.hpp>

#include <frantic/geometry/trimesh3_interface.hpp>
#include <frantic/geometry/xmesh_standard_mesh_interface.hpp>

#include <frantic/graphics/vector4f.hpp>

using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::graphics2d;

TEST( XMeshStandardMeshInterface, Test ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 0 );

    mesh.add_face( 0, 1, 2 );

    // Standard has different type (float16 -> float32)
    mesh.add_vertex_channel_raw( _T("Color"), 3, frantic::channels::data_type_float16 );
    // Standard has greater arity (2 -> 3)
    mesh.add_vertex_channel_raw( _T("TextureCoord"), 2, frantic::channels::data_type_float32 );
    // Standard has lesser arity (4 -> 3)
    mesh.add_vertex_channel_raw( _T("Mapping2"), 4, frantic::channels::data_type_float32 );
    // Non-standard channel name (should be removed)
    mesh.add_vertex_channel_raw( _T("NonStandardVertex"), 3, frantic::channels::data_type_float32 );

    // Standard has different type (uint8 -> uint16)
    mesh.add_face_channel_raw( _T("MaterialID"), 1, frantic::channels::data_type_uint8 );
    // Non-standard channel name (should be removed)
    mesh.add_face_channel_raw( _T("NonStandardFace"), 1, frantic::channels::data_type_uint8 );

    trimesh3_vertex_channel_cvt_accessor<vector3f> colorAcc =
        mesh.get_vertex_channel_cvt_accessor<vector3f>( _T("Color") );
    colorAcc.set( 0, vector3f( 0, 1, 2 ) );

    trimesh3_vertex_channel_cvt_accessor<vector2f> textureCoordAcc =
        mesh.get_vertex_channel_cvt_accessor<vector2f>( _T("TextureCoord") );
    textureCoordAcc.set( 0, vector2f( 3, 4 ) );

    trimesh3_vertex_channel_cvt_accessor<vector4f> mapping2Acc =
        mesh.get_vertex_channel_cvt_accessor<vector4f>( _T("Mapping2") );
    mapping2Acc.set( 0, vector4f( 5, 6, 7, 8 ) );

    trimesh3_face_channel_accessor<boost::uint8_t> materialIdAcc =
        mesh.get_face_channel_accessor<boost::uint8_t>( _T("MaterialID") );
    materialIdAcc[0] = 9;

    std::unique_ptr<frantic::geometry::trimesh3_interface> originalInterface =
        frantic::geometry::trimesh3_interface::create_instance( &mesh );

    ASSERT_TRUE( originalInterface.get() != 0 );

    std::unique_ptr<mesh_interface> standardInterface = create_xmesh_standard_mesh_interface( originalInterface.get() );

    ASSERT_TRUE( standardInterface.get() != 0 );

    mesh_interface::mesh_channel_map& vertexChannels = standardInterface->get_vertex_channels();
    mesh_interface::mesh_channel_map& faceChannels = standardInterface->get_face_channels();

    EXPECT_TRUE( vertexChannels.has_channel( _T("Color") ) );
    const mesh_channel* colorChannel = vertexChannels.get_channel( _T("Color") );
    ASSERT_TRUE( colorChannel != 0 );
    EXPECT_EQ( colorChannel->get_data_arity(), 3 );
    EXPECT_EQ( colorChannel->get_data_type(), frantic::channels::data_type_float32 );
    vector3f color;
    colorChannel->get_value( 0, &color[0] );
    EXPECT_EQ( vector3f( 0, 1, 2 ), color );

    EXPECT_TRUE( vertexChannels.has_channel( _T("TextureCoord") ) );
    const mesh_channel* textureCoordChannel = vertexChannels.get_channel( _T("TextureCoord") );
    ASSERT_TRUE( textureCoordChannel != 0 );
    EXPECT_EQ( textureCoordChannel->get_data_arity(), 3 );
    EXPECT_EQ( textureCoordChannel->get_data_type(), frantic::channels::data_type_float32 );
    vector3f textureCoord( 1 );
    textureCoordChannel->get_value( 0, &textureCoord[0] );
    EXPECT_EQ( vector3f( 3, 4, 0 ), textureCoord );

    EXPECT_TRUE( vertexChannels.has_channel( _T("Mapping2") ) );
    const mesh_channel* mapping2Channel = vertexChannels.get_channel( _T("Mapping2") );
    ASSERT_TRUE( mapping2Channel != 0 );
    EXPECT_EQ( mapping2Channel->get_data_arity(), 3 );
    EXPECT_EQ( mapping2Channel->get_data_type(), frantic::channels::data_type_float32 );
    // we expect a vector3f, but I'm using a vector4f to verify that the fourth
    // component isn't modified
    vector4f mapping2( 1, 1, 1, 1 );
    mapping2Channel->get_value( 0, &mapping2[0] );
    EXPECT_EQ( vector4f( 5, 6, 7, 1 ), mapping2 );

    EXPECT_FALSE( vertexChannels.has_channel( _T("NonStandardVertex") ) );

    EXPECT_TRUE( faceChannels.has_channel( _T("MaterialID") ) );
    const mesh_channel* materialIdChannel = faceChannels.get_channel( _T("MaterialID") );
    ASSERT_TRUE( materialIdChannel != 0 );
    EXPECT_EQ( materialIdChannel->get_data_arity(), 1 );
    EXPECT_EQ( materialIdChannel->get_data_type(), frantic::channels::data_type_uint16 );
    boost::uint16_t materialId;
    materialIdChannel->get_value( 0, &materialId );
    EXPECT_EQ( 9, materialId );

    EXPECT_FALSE( faceChannels.has_channel( _T("NonStandardFace") ) );
}
