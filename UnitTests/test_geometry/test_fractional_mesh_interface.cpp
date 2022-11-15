// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"

#include <frantic/geometry/fractional_mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>

using namespace frantic::geometry;
namespace {

polymesh3_ptr create_test_mesh() {
    // Create the polymesh3 to test
    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );
    builder.add_vertex( 0, 1, 0 );
    builder.add_vertex( 2, 0, 0 );
    builder.add_vertex( 2, 1, 0 );
    int faceIndices[] = { 0, 1, 2, 3 };
    builder.add_polygon( faceIndices, 4 );
    int faceIndices2[] = { 1, 4, 5, 2 };
    builder.add_polygon( faceIndices2, 4 );

    polymesh3_ptr polymesh = builder.finalize();
    EXPECT_TRUE( bool( polymesh ) );

    polymesh->add_empty_vertex_channel( _T( "simple" ), frantic::channels::data_type_int32, 1 );
    polymesh3_vertex_accessor<boost::int32_t> simpleAcc(
        polymesh->get_vertex_accessor<boost::int32_t>( _T( "simple" ) ) );
    simpleAcc.get_vertex( 0 ) = 4;
    simpleAcc.get_vertex( 1 ) = 5;
    simpleAcc.get_vertex( 2 ) = 6;
    simpleAcc.get_vertex( 3 ) = 7;
    simpleAcc.get_vertex( 4 ) = 8;
    simpleAcc.get_vertex( 5 ) = 9;

    polymesh->add_empty_vertex_channel( _T( "custom" ), frantic::channels::data_type_float32, 3, 2 );
    polymesh3_vertex_accessor<frantic::graphics::vector3f> customAcc(
        polymesh->get_vertex_accessor<frantic::graphics::vector3f>( _T( "custom" ) ) );
    customAcc.get_vertex( 0 ).set( 2, 3, 4 );
    customAcc.get_vertex( 1 ).set( 5, 6, 7 );
    customAcc.get_face( 0 ).first[0] = 0;
    customAcc.get_face( 0 ).first[1] = 0;
    customAcc.get_face( 0 ).first[2] = 1;
    customAcc.get_face( 0 ).first[3] = 1;
    customAcc.get_face( 1 ).first[0] = 0;
    customAcc.get_face( 1 ).first[1] = 0;
    customAcc.get_face( 1 ).first[2] = 1;
    customAcc.get_face( 1 ).first[3] = 1;

    polymesh->add_empty_face_channel( _T( "face" ), frantic::channels::data_type_uint8, 1 );
    polymesh3_face_accessor<boost::uint8_t> faceAcc( polymesh->get_face_accessor<boost::uint8_t>( _T( "face" ) ) );
    faceAcc.get_face( 0 ) = 12;
    faceAcc.get_face( 1 ) = 13;

    return polymesh;
}
} // anonymous namespace

TEST( FractionalInterface, WholeFaces ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    const mesh_interface_ptr poly( polymesh3_interface::create_const_instance( m ).release() );
    const mesh_interface_ptr mesh =
        boost::shared_ptr<fractional_face_interface>( new fractional_face_interface( poly.get(), 1 ) );
    ASSERT_TRUE( mesh.get() != NULL );
    // is_valid()
    EXPECT_TRUE( mesh->is_valid() );

    EXPECT_TRUE( is_equal( mesh, poly ) );
}

TEST( FractionalInterface, WholeFacesLimited ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_face_interface>( new fractional_face_interface( poly.get(), 1, 1 ) ) );
    ASSERT_TRUE( mesh.get() != NULL );
    // is_valid()
    EXPECT_TRUE( mesh->is_valid() );

    // get_num_verts()
    EXPECT_EQ( mesh->get_num_verts(), 4 );

    // get_vert()
    frantic::graphics::vector3f vert;
    mesh->get_vert( 0, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 0 ) );
    mesh->get_vert( 1, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 0, 0 ) );
    mesh->get_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 1, 0 ) );
    mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 0, 1, 0 ) );
    EXPECT_THROW( mesh->get_vert( 4, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );
    EXPECT_THROW( mesh->get_vert( 5, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );

    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), mesh->get_vert( 3 ) );

    // get_num_faces()
    EXPECT_EQ( mesh->get_num_faces(), 1 );

    // get_num_face_verts()
    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 4 );
    EXPECT_THROW( mesh->get_num_face_verts( 1 ), out_of_range );

    // get_face_vert_index()
    EXPECT_EQ( mesh->get_face_vert_index( 0, 0 ), 0 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 1 ), 1 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 2 ), 2 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 3 ), 3 );
    EXPECT_THROW( mesh->get_face_vert_index( 1, 0 ), out_of_range );

    // get_face_vert_indices()
    std::size_t faceVertIndices[4];
    mesh->get_face_vert_indices( 0, faceVertIndices );
    EXPECT_EQ( faceVertIndices[0], 0 );
    EXPECT_EQ( faceVertIndices[1], 1 );
    EXPECT_EQ( faceVertIndices[2], 2 );
    EXPECT_EQ( faceVertIndices[3], 3 );
    EXPECT_THROW( mesh->get_face_vert_indices( 1, faceVertIndices ), out_of_range );

    // get_face_verts()
    float faceVerts[4][3];
    mesh->get_face_verts( 0, faceVerts );
    EXPECT_EQ( faceVerts[0][0], 0 );
    EXPECT_EQ( faceVerts[0][1], 0 );
    EXPECT_EQ( faceVerts[0][2], 0 );
    EXPECT_EQ( faceVerts[1][0], 1 );
    EXPECT_EQ( faceVerts[1][1], 0 );
    EXPECT_EQ( faceVerts[1][2], 0 );
    EXPECT_EQ( faceVerts[2][0], 1 );
    EXPECT_EQ( faceVerts[2][1], 1 );
    EXPECT_EQ( faceVerts[2][2], 0 );
    EXPECT_EQ( faceVerts[3][0], 0 );
    EXPECT_EQ( faceVerts[3][1], 1 );
    EXPECT_EQ( faceVerts[3][2], 0 );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
        std::set<frantic::tstring> vertexChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( vertexChannelMap.begin() ),
             ie( vertexChannelMap.end() );
             i != ie; ++i ) {
            vertexChannelNames.insert( i->first );
        }
        EXPECT_EQ( vertexChannelNames.size(), 2 );
        EXPECT_EQ( vertexChannelNames.count( _T( "simple" ) ), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "custom" ) ), 1 );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "simple" ) ) );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "custom" ) ) );

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_EQ( simpleChannel->get_name(), _T( "simple" ) );
            EXPECT_EQ( simpleChannel->get_channel_type(), mesh_channel::vertex );
            EXPECT_EQ( simpleChannel->get_num_elements(), 4 );
            EXPECT_EQ( simpleChannel->get_num_faces(), 1 );
            EXPECT_EQ( simpleChannel->get_data_type(), frantic::channels::data_type_int32 );
            EXPECT_EQ( simpleChannel->get_data_arity(), 1 );
            EXPECT_EQ( simpleChannel->get_element_size(),
                       frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_int32 ) );

            boost::int32_t i;
            simpleChannel->get_value( 0, &i );
            EXPECT_EQ( i, 4 );
            simpleChannel->get_value( 1, &i );
            EXPECT_EQ( i, 5 );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 6 );
            simpleChannel->get_value( 3, &i );
            EXPECT_EQ( i, 7 );

            EXPECT_EQ( simpleChannel->get_num_face_verts( 0 ), 4 );

            EXPECT_EQ( simpleChannel->get_fv_index( 0, 0 ), 0 );
            EXPECT_EQ( simpleChannel->get_fv_index( 0, 1 ), 1 );
            EXPECT_EQ( simpleChannel->get_fv_index( 0, 2 ), 2 );
            EXPECT_EQ( simpleChannel->get_fv_index( 0, 3 ), 3 );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel != NULL );
            EXPECT_EQ( customChannel->get_name(), _T( "custom" ) );
            EXPECT_EQ( customChannel->get_channel_type(), mesh_channel::face_vertex );
            EXPECT_EQ( customChannel->get_num_elements(), 2 );
            EXPECT_EQ( customChannel->get_num_faces(), 1 );
            EXPECT_EQ( customChannel->get_data_type(), frantic::channels::data_type_float32 );
            EXPECT_EQ( customChannel->get_data_arity(), 3 );
            EXPECT_EQ( customChannel->get_element_size(), sizeof( frantic::graphics::vector3f ) );

            frantic::graphics::vector3f v;
            customChannel->get_value( 0, &v[0] );
            EXPECT_EQ( v, frantic::graphics::vector3f( 2, 3, 4 ) );
            customChannel->get_value( 1, &v[0] );
            EXPECT_EQ( v, frantic::graphics::vector3f( 5, 6, 7 ) );

            EXPECT_EQ( customChannel->get_num_face_verts( 0 ), 4 );

            EXPECT_EQ( customChannel->get_fv_index( 0, 0 ), 0 );
            EXPECT_EQ( customChannel->get_fv_index( 0, 1 ), 0 );
            EXPECT_EQ( customChannel->get_fv_index( 0, 2 ), 1 );
            EXPECT_EQ( customChannel->get_fv_index( 0, 3 ), 1 );
        }
    }

    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
        std::set<frantic::tstring> faceChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( faceChannelMap.begin() ), ie( faceChannelMap.end() );
             i != ie; ++i ) {
            faceChannelNames.insert( i->first );
        }
        EXPECT_EQ( faceChannelNames.size(), 1 );
        EXPECT_EQ( faceChannelNames.count( _T( "face" ) ), 1 );
        EXPECT_TRUE( faceChannelMap.has_channel( _T( "face" ) ) );

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel != NULL );
        EXPECT_EQ( faceChannel->get_name(), _T( "face" ) );
        EXPECT_EQ( faceChannel->get_channel_type(), mesh_channel::face );
        EXPECT_EQ( faceChannel->get_num_elements(), 1 );
        EXPECT_EQ( faceChannel->get_num_faces(), 1 );
        EXPECT_EQ( faceChannel->get_data_type(), frantic::channels::data_type_uint8 );
        EXPECT_EQ( faceChannel->get_data_arity(), 1 );
        EXPECT_EQ( faceChannel->get_element_size(),
                   frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_uint8 ) );

        boost::uint8_t i;
        faceChannel->get_value( 0, &i );
        EXPECT_EQ( i, 12 );
    }
}

TEST( FractionalInterface, FacesWholeReadOnly ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_face_interface>( new fractional_face_interface( poly.get(), 1 ) ) );

    ASSERT_TRUE( mesh.get() != NULL );

    // is_read_only()
    EXPECT_TRUE( mesh->is_read_only() );

    // set_vert()
    frantic::graphics::vector3f vert;
    vert.set( -1 );
    EXPECT_THROW( mesh->set_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), runtime_error );
    // should keep its original value
    vert.set( 0 );
    mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), vert );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_FALSE( simpleChannel->is_writeable() );

            boost::int32_t i;

            // should keep its original value after set_value()
            const boost::int32_t eight = 8;
            simpleChannel->set_value( 2, &eight );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 6 );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel != NULL );
            EXPECT_FALSE( customChannel->is_writeable() );

            frantic::graphics::vector3f v;

            // should keep its original value after set_value()
            customChannel->set_value( 0, &( frantic::graphics::vector3f( 0, 1, 2 )[0] ) );
            customChannel->get_value( 0, &v[0] );
            EXPECT_EQ( v, frantic::graphics::vector3f( 2, 3, 4 ) );

            EXPECT_THROW( customChannel->set_fv_index( 0, 3, 0 ), runtime_error );
            EXPECT_EQ( customChannel->get_fv_index( 0, 3 ), 1 );
        }
    }

    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel != NULL );
        EXPECT_FALSE( faceChannel->is_writeable() );

        boost::uint8_t i;

        // should keep its original value after set_value()
        const boost::uint8_t eight = 8;
        faceChannel->set_value( 0, &eight );
        faceChannel->get_value( 0, &i );
        EXPECT_EQ( i, 12 );
    }
}

TEST( FractionalInterface, FacesFractional ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_face_interface>( new fractional_face_interface( poly.get(), .5 ) ) );
    ASSERT_TRUE( mesh.get() != NULL );

    // is_valid()
    EXPECT_TRUE( mesh->is_valid() );

    // get_num_verts()
    EXPECT_EQ( mesh->get_num_verts(), 4 );

    // get_vert()
    frantic::graphics::vector3f vert;
    mesh->get_vert( 0, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 0, 0 ) );
    mesh->get_vert( 1, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 1, 0 ) );
    mesh->get_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 2, 0, 0 ) );
    mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 2, 1, 0 ) );
    EXPECT_THROW( mesh->get_vert( 4, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );

    EXPECT_EQ( frantic::graphics::vector3f( 2, 1, 0 ), mesh->get_vert( 3 ) );

    // get_num_faces()
    EXPECT_EQ( mesh->get_num_faces(), 1 );

    // get_num_face_verts()
    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 4 );
    EXPECT_THROW( mesh->get_num_face_verts( 1 ), out_of_range );

    // get_face_vert_index()
    EXPECT_EQ( mesh->get_face_vert_index( 0, 0 ), 0 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 1 ), 2 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 2 ), 3 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 3 ), 1 );

    // get_face_vert_indices()
    std::size_t faceVertIndices[4];
    mesh->get_face_vert_indices( 0, faceVertIndices );
    EXPECT_EQ( faceVertIndices[0], 0 );
    EXPECT_EQ( faceVertIndices[1], 2 );
    EXPECT_EQ( faceVertIndices[2], 3 );
    EXPECT_EQ( faceVertIndices[3], 1 );

    // get_face_verts()
    float faceVerts[4][3];
    mesh->get_face_verts( 0, faceVerts );
    EXPECT_EQ( faceVerts[0][0], 1 );
    EXPECT_EQ( faceVerts[0][1], 0 );
    EXPECT_EQ( faceVerts[0][2], 0 );
    EXPECT_EQ( faceVerts[1][0], 2 );
    EXPECT_EQ( faceVerts[1][1], 0 );
    EXPECT_EQ( faceVerts[1][2], 0 );
    EXPECT_EQ( faceVerts[2][0], 2 );
    EXPECT_EQ( faceVerts[2][1], 1 );
    EXPECT_EQ( faceVerts[2][2], 0 );
    EXPECT_EQ( faceVerts[3][0], 1 );
    EXPECT_EQ( faceVerts[3][1], 1 );
    EXPECT_EQ( faceVerts[3][2], 0 );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
        std::set<frantic::tstring> vertexChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( vertexChannelMap.begin() ),
             ie( vertexChannelMap.end() );
             i != ie; ++i ) {
            vertexChannelNames.insert( i->first );
        }
        EXPECT_EQ( vertexChannelNames.size(), 2 );
        EXPECT_EQ( vertexChannelNames.count( _T( "simple" ) ), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "custom" ) ), 1 );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "simple" ) ) );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "custom" ) ) );

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_EQ( simpleChannel->get_name(), _T( "simple" ) );
            EXPECT_EQ( simpleChannel->get_channel_type(), mesh_channel::vertex );
            EXPECT_EQ( simpleChannel->get_num_elements(), 4 );
            EXPECT_EQ( simpleChannel->get_num_faces(), 1 );
            EXPECT_EQ( simpleChannel->get_data_type(), frantic::channels::data_type_int32 );
            EXPECT_EQ( simpleChannel->get_data_arity(), 1 );
            EXPECT_EQ( simpleChannel->get_element_size(),
                       frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_int32 ) );

            boost::int32_t i;
            simpleChannel->get_value( 0, &i );
            EXPECT_EQ( i, 5 );
            simpleChannel->get_value( 1, &i );
            EXPECT_EQ( i, 6 );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 8 );
            simpleChannel->get_value( 3, &i );
            EXPECT_EQ( i, 9 );

            EXPECT_EQ( simpleChannel->get_num_face_verts( 0 ), 4 );

            EXPECT_EQ( simpleChannel->get_fv_index( 0, 0 ), 0 );
            EXPECT_EQ( simpleChannel->get_fv_index( 0, 1 ), 2 );
            EXPECT_EQ( simpleChannel->get_fv_index( 0, 2 ), 3 );
            EXPECT_EQ( simpleChannel->get_fv_index( 0, 3 ), 1 );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel != NULL );
            EXPECT_EQ( customChannel->get_name(), _T( "custom" ) );
            EXPECT_EQ( customChannel->get_channel_type(), mesh_channel::face_vertex );
            EXPECT_EQ( customChannel->get_num_elements(), 2 );
            EXPECT_EQ( customChannel->get_num_faces(), 1 );
            EXPECT_EQ( customChannel->get_data_type(), frantic::channels::data_type_float32 );
            EXPECT_EQ( customChannel->get_data_arity(), 3 );
            EXPECT_EQ( customChannel->get_element_size(), sizeof( frantic::graphics::vector3f ) );

            frantic::graphics::vector3f v;
            customChannel->get_value( 0, &v[0] );
            EXPECT_EQ( v, frantic::graphics::vector3f( 2, 3, 4 ) );
            customChannel->get_value( 1, &v[0] );
            EXPECT_EQ( v, frantic::graphics::vector3f( 5, 6, 7 ) );

            EXPECT_EQ( customChannel->get_num_face_verts( 0 ), 4 );

            EXPECT_EQ( customChannel->get_fv_index( 0, 0 ), 0 );
            EXPECT_EQ( customChannel->get_fv_index( 0, 1 ), 0 );
            EXPECT_EQ( customChannel->get_fv_index( 0, 2 ), 1 );
            EXPECT_EQ( customChannel->get_fv_index( 0, 3 ), 1 );
        }
    }

    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
        std::set<frantic::tstring> faceChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( faceChannelMap.begin() ), ie( faceChannelMap.end() );
             i != ie; ++i ) {
            faceChannelNames.insert( i->first );
        }
        EXPECT_EQ( faceChannelNames.size(), 1 );
        EXPECT_EQ( faceChannelNames.count( _T( "face" ) ), 1 );
        EXPECT_TRUE( faceChannelMap.has_channel( _T( "face" ) ) );

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel != NULL );
        EXPECT_EQ( faceChannel->get_name(), _T( "face" ) );
        EXPECT_EQ( faceChannel->get_channel_type(), mesh_channel::face );
        EXPECT_EQ( faceChannel->get_num_elements(), 1 );
        EXPECT_EQ( faceChannel->get_num_faces(), 1 );
        EXPECT_EQ( faceChannel->get_data_type(), frantic::channels::data_type_uint8 );
        EXPECT_EQ( faceChannel->get_data_arity(), 1 );
        EXPECT_EQ( faceChannel->get_element_size(),
                   frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_uint8 ) );

        boost::uint8_t i;
        faceChannel->get_value( 0, &i );
        EXPECT_EQ( i, 13 );
    }
}

TEST( FractionalInterface, FacesAdjacency ) {

    const size_t numMeshes = 4;
    const double fraction = .5;

    std::unique_ptr<mesh_interface> meshes[numMeshes];
    // Simple mesh, just a single triangle
    meshes[0] = polymesh3_interface::create_const_instance( make_triangle_polymesh() );

    // A cube mesh with 6 quads
    meshes[1] = polymesh3_interface::create_const_instance( make_cube_polymesh() );

    // A cube mesh with 4 quads and two distinct holes
    std::set<size_t> holes;
    holes.insert( 0 );
    holes.insert( 1 );
    meshes[2] = polymesh3_interface::create_const_instance( make_cube_polymesh( holes ) );

    // A cube mesh with 4 quads and one hole
    holes.clear();
    holes.insert( 0 );
    holes.insert( 2 );
    meshes[3] = polymesh3_interface::create_const_instance( make_cube_polymesh( holes ) );

    for( size_t meshId = 0; meshId < numMeshes; ++meshId ) {
        mesh_interface_ptr mesh = mesh_interface_ptr( std::unique_ptr<fractional_face_interface>(
            new fractional_face_interface( meshes[meshId].get(), fraction ) ) );

        // Build up a map of all face adjacencies for testing
        typedef std::map<std::pair<size_t, size_t>, size_t> face_map_t;
        face_map_t leftFaces;

        for( size_t i = 0; i < mesh->get_num_faces(); ++i ) {
            const size_t faceSize = mesh->get_num_face_verts( i );
            for( size_t j = 0; j < faceSize; ++j ) {
                leftFaces[std::make_pair( mesh->get_face_vert_index( i, j ),
                                          mesh->get_face_vert_index( i, ( j + 1 ) % faceSize ) )] = i;
            }
        }

        mesh->init_adjacency();

        for( size_t i = 0; i < mesh->get_num_verts(); ++i ) {

            vertex_iterator vIt;
            ASSERT_TRUE( mesh->init_vertex_iterator( vIt, i ) );
            std::set<size_t> adjacents;

            do {
                size_t endpoint = mesh->get_edge_endpoint( vIt );

                face_map_t::iterator leftFace = leftFaces.find( std::make_pair( i, endpoint ) );
                // reverse the edge direction to get the other face
                face_map_t::iterator rightFace = leftFaces.find( std::make_pair( endpoint, i ) );

                EXPECT_TRUE( leftFace != leftFaces.end() || rightFace != leftFaces.end() )
                    << "Error, edge does not exist.";

                EXPECT_TRUE( adjacents.find( endpoint ) == adjacents.end() )
                    << "Adjacent vertices should only appear once.";

                if( leftFace != leftFaces.end() ) {
                    EXPECT_EQ( leftFace->second, mesh->get_edge_left_face( vIt ) )
                        << "Incorrect face to the left of this edge.";
                } else {
                    EXPECT_EQ( mesh_interface::HOLE_INDEX, mesh->get_edge_left_face( vIt ) )
                        << "Not labeled as a hole face.";
                }

                if( rightFace != leftFaces.end() ) {
                    EXPECT_EQ( rightFace->second, mesh->get_edge_right_face( vIt ) )
                        << "Incorrect face to the right of this edge.";
                } else {
                    EXPECT_EQ( mesh_interface::HOLE_INDEX, mesh->get_edge_right_face( vIt ) )
                        << "Not labeled as a hole face.";
                }

                EXPECT_EQ( leftFace == leftFaces.end() || rightFace == leftFaces.end(), mesh->is_edge_boundary( vIt ) )
                    << "Did not identify edge as boundary correctly.";

                adjacents.insert( endpoint );
            } while( mesh->advance_vertex_iterator( vIt ) );
        }

        for( size_t i = 0; i < mesh->get_num_faces(); ++i ) {
            face_iterator fIt;

            const size_t faceSize = mesh->get_num_face_verts( i );

            mesh->init_face_iterator( fIt, i );

            size_t currentVertex = 0;

            do {
                size_t curr = mesh->get_face_vert_index( i, currentVertex );
                size_t next = mesh->get_face_vert_index( i, ( currentVertex + 1 ) % faceSize );

                face_map_t::iterator oppositeFace = leftFaces.find( std::make_pair( next, curr ) );

                if( oppositeFace != leftFaces.end() ) {
                    EXPECT_EQ( oppositeFace->second, mesh->get_face_neighbor( fIt ) ) << "Incorrect neighbour face.";
                } else {
                    EXPECT_EQ( mesh_interface::HOLE_INDEX, mesh->get_face_neighbor( fIt ) ) << "Expected hole face.";
                }
                ++currentVertex;
            } while( mesh->advance_face_iterator( fIt ) );

            EXPECT_EQ( currentVertex, faceSize ) << "Incorrect number of face iterations.";
        }
    }
}

TEST( FractionalInterface, FacesRegressionTests ) {
    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );
    builder.add_vertex( 0, 1, 0 );
    builder.add_vertex( 2, 0, 0 );

    int faceIndices[] = { 0, 1, 2, 3 };
    builder.add_polygon( faceIndices, 4 );
    int faceIndices2[] = { 1, 4, 2 };
    builder.add_polygon( faceIndices2, 3 );

    polymesh3_ptr polymesh = builder.finalize();
    EXPECT_TRUE( bool( polymesh ) );

    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( polymesh );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_face_interface>( new fractional_face_interface( poly.get(), .5 ) ) );

    ASSERT_TRUE( mesh.get() != NULL );

    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 3 );

    // get_face_vert_indices()
    std::size_t faceVertIndices[3];
    mesh->get_face_vert_indices( 0, faceVertIndices );
    EXPECT_EQ( faceVertIndices[0], 0 );
    EXPECT_EQ( faceVertIndices[1], 2 );
    EXPECT_EQ( faceVertIndices[2], 1 );
}

TEST( FractionalInterface, WholeVertex ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_vertex_interface>( new fractional_vertex_interface( poly.get(), 1 ) ) );
    ASSERT_TRUE( mesh.get() != NULL );

    // is_valid()
    EXPECT_TRUE( mesh->is_valid() );

    // get_num_verts()
    EXPECT_EQ( mesh->get_num_verts(), 6 );

    // get_vert()
    frantic::graphics::vector3f vert;
    mesh->get_vert( 0, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 0 ) );
    mesh->get_vert( 1, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 0, 0 ) );
    mesh->get_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 1, 0 ) );
    mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 0, 1, 0 ) );
    mesh->get_vert( 4, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 2, 0, 0 ) );
    mesh->get_vert( 5, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 2, 1, 0 ) );

    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), mesh->get_vert( 3 ) );

    // get_num_faces()
    EXPECT_EQ( mesh->get_num_faces(), 0 );

    // get_num_face_verts()
    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 0 );
    EXPECT_EQ( mesh->get_num_face_verts( 1 ), 0 );

    // get_face_vert_index()
    EXPECT_THROW( mesh->get_face_vert_index( 0, 0 ), runtime_error );

    // get_face_vert_indices()
    std::size_t faceVertIndices[4];
    EXPECT_THROW( mesh->get_face_vert_indices( 0, faceVertIndices ), runtime_error );

    // get_face_verts()
    float faceVerts[4][3];
    EXPECT_THROW( mesh->get_face_verts( 0, faceVerts ), runtime_error );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
        std::set<frantic::tstring> vertexChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( vertexChannelMap.begin() ),
             ie( vertexChannelMap.end() );
             i != ie; ++i ) {
            vertexChannelNames.insert( i->first );
        }
        EXPECT_EQ( vertexChannelNames.size(), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "simple" ) ), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "custom" ) ), 0 );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "simple" ) ) );
        EXPECT_FALSE( vertexChannelMap.has_channel( _T( "custom" ) ) );

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_FALSE( simpleChannel == NULL );
            EXPECT_EQ( simpleChannel->get_name(), _T( "simple" ) );
            EXPECT_EQ( simpleChannel->get_channel_type(), mesh_channel::vertex );
            EXPECT_EQ( simpleChannel->get_num_elements(), 6 );
            EXPECT_EQ( simpleChannel->get_num_faces(), 0 );
            EXPECT_EQ( simpleChannel->get_data_type(), frantic::channels::data_type_int32 );
            EXPECT_EQ( simpleChannel->get_data_arity(), 1 );
            EXPECT_EQ( simpleChannel->get_element_size(),
                       frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_int32 ) );

            boost::int32_t i;
            simpleChannel->get_value( 0, &i );
            EXPECT_EQ( i, 4 );
            simpleChannel->get_value( 1, &i );
            EXPECT_EQ( i, 5 );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 6 );
            simpleChannel->get_value( 3, &i );
            EXPECT_EQ( i, 7 );

            EXPECT_EQ( simpleChannel->get_num_face_verts( 0 ), 0 );

            EXPECT_THROW( simpleChannel->get_fv_index( 0, 0 ), runtime_error );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel == NULL );
        }
    }

    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
        std::set<frantic::tstring> faceChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( faceChannelMap.begin() ), ie( faceChannelMap.end() );
             i != ie; ++i ) {
            faceChannelNames.insert( i->first );
        }
        EXPECT_EQ( faceChannelNames.size(), 0 );
        EXPECT_EQ( faceChannelNames.count( _T( "face" ) ), 0 );
        EXPECT_FALSE( faceChannelMap.has_channel( _T( "face" ) ) );

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel == NULL );
    }
}

TEST( FractionalInterface, WholeVertexLimited ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_vertex_interface>( new fractional_vertex_interface( poly.get(), 1, 3 ) ) );

    ASSERT_TRUE( mesh.get() != NULL );

    // is_valid()
    EXPECT_TRUE( mesh->is_valid() );

    // get_num_verts()
    EXPECT_EQ( mesh->get_num_verts(), 3 );

    // get_vert()
    frantic::graphics::vector3f vert;
    mesh->get_vert( 0, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 0 ) );
    mesh->get_vert( 1, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 0, 0 ) );
    mesh->get_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 1, 0 ) );
    EXPECT_THROW( mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );

    EXPECT_EQ( frantic::graphics::vector3f( 1, 1, 0 ), mesh->get_vert( 2 ) );

    // get_num_faces()
    EXPECT_EQ( mesh->get_num_faces(), 0 );

    // get_num_face_verts()
    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 0 );
    EXPECT_EQ( mesh->get_num_face_verts( 1 ), 0 );

    // get_face_vert_index()
    EXPECT_THROW( mesh->get_face_vert_index( 0, 0 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 0, 1 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 1, 0 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 1, 1 ), runtime_error );

    // get_face_vert_indices()
    std::size_t faceVertIndices[4];
    EXPECT_THROW( mesh->get_face_vert_indices( 0, faceVertIndices ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_indices( 1, faceVertIndices ), runtime_error );

    // get_face_verts()
    float faceVerts[4][3];
    EXPECT_THROW( mesh->get_face_verts( 0, faceVerts ), runtime_error );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
        std::set<frantic::tstring> vertexChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( vertexChannelMap.begin() ),
             ie( vertexChannelMap.end() );
             i != ie; ++i ) {
            vertexChannelNames.insert( i->first );
        }
        EXPECT_EQ( vertexChannelNames.size(), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "simple" ) ), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "custom" ) ), 0 );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "simple" ) ) );
        EXPECT_FALSE( vertexChannelMap.has_channel( _T( "custom" ) ) );

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_FALSE( simpleChannel == NULL );
            EXPECT_EQ( simpleChannel->get_name(), _T( "simple" ) );
            EXPECT_EQ( simpleChannel->get_channel_type(), mesh_channel::vertex );
            EXPECT_EQ( simpleChannel->get_num_elements(), 3 );
            EXPECT_EQ( simpleChannel->get_num_faces(), 0 );
            EXPECT_EQ( simpleChannel->get_data_type(), frantic::channels::data_type_int32 );
            EXPECT_EQ( simpleChannel->get_data_arity(), 1 );
            EXPECT_EQ( simpleChannel->get_element_size(),
                       frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_int32 ) );

            boost::int32_t i;
            simpleChannel->get_value( 0, &i );
            EXPECT_EQ( i, 4 );
            simpleChannel->get_value( 1, &i );
            EXPECT_EQ( i, 5 );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 6 );
            EXPECT_THROW( simpleChannel->get_value( 3, &i ), out_of_range );

            EXPECT_EQ( simpleChannel->get_num_face_verts( 0 ), 0 );

            EXPECT_THROW( simpleChannel->get_fv_index( 0, 0 ), runtime_error );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel == NULL );
        }
    }

    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
        std::set<frantic::tstring> faceChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( faceChannelMap.begin() ), ie( faceChannelMap.end() );
             i != ie; ++i ) {
            faceChannelNames.insert( i->first );
        }
        EXPECT_EQ( faceChannelNames.size(), 0 );
        EXPECT_EQ( faceChannelNames.count( _T( "face" ) ), 0 );
        EXPECT_FALSE( faceChannelMap.has_channel( _T( "face" ) ) );

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel == NULL );
    }
}
TEST( FractionalInterface, VertexFractional ) {
    using namespace std;
    using namespace boost;
    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_vertex_interface>( new fractional_vertex_interface( poly.get(), .5 ) ) );
    ASSERT_TRUE( mesh.get() != NULL );

    // is_valid()
    EXPECT_TRUE( mesh->is_valid() );

    // get_num_verts()
    EXPECT_EQ( mesh->get_num_verts(), 3 );

    // get_vert()
    frantic::graphics::vector3f vert;
    mesh->get_vert( 0, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 0, 0 ) );
    mesh->get_vert( 1, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 0, 1, 0 ) );
    mesh->get_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 2, 1, 0 ) );
    EXPECT_THROW( mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );
    EXPECT_THROW( mesh->get_vert( 4, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );
    EXPECT_THROW( mesh->get_vert( 5, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );

    EXPECT_EQ( frantic::graphics::vector3f( 2, 1, 0 ), mesh->get_vert( 2 ) );

    // get_num_faces()
    EXPECT_EQ( mesh->get_num_faces(), 0 );

    // get_num_face_verts()
    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 0 );
    EXPECT_EQ( mesh->get_num_face_verts( 1 ), 0 );

    // get_face_vert_index()
    EXPECT_THROW( mesh->get_face_vert_index( 0, 0 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 0, 1 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 1, 0 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 1, 1 ), runtime_error );

    // get_face_vert_indices()
    std::size_t faceVertIndices[4];
    EXPECT_THROW( mesh->get_face_vert_indices( 0, faceVertIndices ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_indices( 1, faceVertIndices ), runtime_error );

    // get_face_verts()
    float faceVerts[4][3];
    EXPECT_THROW( mesh->get_face_verts( 0, faceVerts ), runtime_error );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
        std::set<frantic::tstring> vertexChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( vertexChannelMap.begin() ),
             ie( vertexChannelMap.end() );
             i != ie; ++i ) {
            vertexChannelNames.insert( i->first );
        }
        EXPECT_EQ( vertexChannelNames.size(), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "simple" ) ), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "custom" ) ), 0 );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "simple" ) ) );
        EXPECT_FALSE( vertexChannelMap.has_channel( _T( "custom" ) ) );

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_EQ( simpleChannel->get_name(), _T( "simple" ) );
            EXPECT_EQ( simpleChannel->get_channel_type(), mesh_channel::vertex );
            EXPECT_EQ( simpleChannel->get_num_elements(), 3 );
            EXPECT_EQ( simpleChannel->get_num_faces(), 0 );
            EXPECT_EQ( simpleChannel->get_data_type(), frantic::channels::data_type_int32 );
            EXPECT_EQ( simpleChannel->get_data_arity(), 1 );
            EXPECT_EQ( simpleChannel->get_element_size(),
                       frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_int32 ) );

            boost::int32_t i;
            simpleChannel->get_value( 0, &i );
            EXPECT_EQ( i, 5 );
            simpleChannel->get_value( 1, &i );
            EXPECT_EQ( i, 7 );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 9 );
            EXPECT_THROW( simpleChannel->get_value( 3, &i ), out_of_range );

            EXPECT_EQ( simpleChannel->get_num_face_verts( 0 ), 0 );

            EXPECT_THROW( simpleChannel->get_fv_index( 0, 0 ), runtime_error );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel == NULL );
        }
    }

    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
        std::set<frantic::tstring> faceChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( faceChannelMap.begin() ), ie( faceChannelMap.end() );
             i != ie; ++i ) {
            faceChannelNames.insert( i->first );
        }
        EXPECT_EQ( faceChannelNames.size(), 0 );
        EXPECT_EQ( faceChannelNames.count( _T( "face" ) ), 0 );
        EXPECT_FALSE( faceChannelMap.has_channel( _T( "face" ) ) );

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel == NULL );
    }
}

TEST( FractionalInterface, VertexFractionalLimited ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_vertex_interface>( new fractional_vertex_interface( poly.get(), .5, 2 ) ) );
    ASSERT_TRUE( mesh.get() != NULL );

    // is_valid()
    EXPECT_TRUE( mesh->is_valid() );

    // get_num_verts()
    EXPECT_EQ( mesh->get_num_verts(), 2 );

    // get_vert()
    frantic::graphics::vector3f vert;
    mesh->get_vert( 0, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 1, 0, 0 ) );
    mesh->get_vert( 1, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( vert, frantic::graphics::vector3f( 0, 1, 0 ) );
    EXPECT_THROW( mesh->get_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );
    EXPECT_THROW( mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );
    EXPECT_THROW( mesh->get_vert( 4, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );
    EXPECT_THROW( mesh->get_vert( 5, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), out_of_range );

    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), mesh->get_vert( 1 ) );

    // get_num_faces()
    EXPECT_EQ( mesh->get_num_faces(), 0 );

    // get_num_face_verts()
    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 0 );
    EXPECT_EQ( mesh->get_num_face_verts( 1 ), 0 );

    // get_face_vert_index()
    EXPECT_THROW( mesh->get_face_vert_index( 0, 0 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 0, 1 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 1, 0 ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_index( 1, 1 ), runtime_error );

    // get_face_vert_indices()
    std::size_t faceVertIndices[4];
    EXPECT_THROW( mesh->get_face_vert_indices( 0, faceVertIndices ), runtime_error );
    EXPECT_THROW( mesh->get_face_vert_indices( 1, faceVertIndices ), runtime_error );

    // get_face_verts()
    float faceVerts[4][3];
    EXPECT_THROW( mesh->get_face_verts( 0, faceVerts ), runtime_error );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
        std::set<frantic::tstring> vertexChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( vertexChannelMap.begin() ),
             ie( vertexChannelMap.end() );
             i != ie; ++i ) {
            vertexChannelNames.insert( i->first );
        }
        EXPECT_EQ( vertexChannelNames.size(), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "simple" ) ), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T( "custom" ) ), 0 );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T( "simple" ) ) );
        EXPECT_FALSE( vertexChannelMap.has_channel( _T( "custom" ) ) );

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_EQ( simpleChannel->get_name(), _T( "simple" ) );
            EXPECT_EQ( simpleChannel->get_channel_type(), mesh_channel::vertex );
            EXPECT_EQ( simpleChannel->get_num_elements(), 2 );
            EXPECT_EQ( simpleChannel->get_num_faces(), 0 );
            EXPECT_EQ( simpleChannel->get_data_type(), frantic::channels::data_type_int32 );
            EXPECT_EQ( simpleChannel->get_data_arity(), 1 );
            EXPECT_EQ( simpleChannel->get_element_size(),
                       frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_int32 ) );

            boost::int32_t i;
            simpleChannel->get_value( 0, &i );
            EXPECT_EQ( i, 5 );
            simpleChannel->get_value( 1, &i );
            EXPECT_EQ( i, 7 );
            EXPECT_THROW( simpleChannel->get_value( 2, &i ), out_of_range );
            EXPECT_THROW( simpleChannel->get_value( 3, &i ), out_of_range );

            EXPECT_EQ( simpleChannel->get_num_face_verts( 0 ), 0 );

            EXPECT_THROW( simpleChannel->get_fv_index( 0, 0 ), runtime_error );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel == NULL );
        }
    }

    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
        std::set<frantic::tstring> faceChannelNames;
        for( mesh_interface::mesh_channel_map::const_iterator i( faceChannelMap.begin() ), ie( faceChannelMap.end() );
             i != ie; ++i ) {
            faceChannelNames.insert( i->first );
        }
        EXPECT_EQ( faceChannelNames.size(), 0 );
        EXPECT_EQ( faceChannelNames.count( _T( "face" ) ), 0 );
        EXPECT_FALSE( faceChannelMap.has_channel( _T( "face" ) ) );

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel == NULL );
    }
}

TEST( FractionalInterface, VertexReadOnly ) {
    using namespace std;
    using namespace boost;

    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_const_instance( m );
    mesh_interface_ptr mesh = mesh_interface_ptr(
        boost::shared_ptr<fractional_vertex_interface>( new fractional_vertex_interface( poly.get(), .5 ) ) );

    ASSERT_TRUE( mesh.get() != NULL );

    // is_read_only()
    EXPECT_TRUE( mesh->is_read_only() );

    // set_vert()
    frantic::graphics::vector3f vert;
    vert.set( -1 );
    EXPECT_THROW( mesh->set_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) ), runtime_error );
    // should keep its original value
    vert.set( 0 );
    mesh->get_vert( 2, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( frantic::graphics::vector3f( 2, 1, 0 ), vert );

    { // scope for vertexChannelMap
        mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();

        { // scope for simpleChannel
            const mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T( "simple" ) );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_FALSE( simpleChannel->is_writeable() );

            boost::int32_t i;

            // should keep its original value after set_value()
            const boost::int32_t eight = 8;
            simpleChannel->set_value( 2, &eight );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 9 );
        }

        { // scope for customChannel
            const mesh_channel* customChannel = vertexChannelMap.get_channel( _T( "custom" ) );
            ASSERT_TRUE( customChannel == NULL );
        }
    }
    { // scope for faceChannelMap
        mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();

        const mesh_channel* faceChannel = faceChannelMap.get_channel( _T( "face" ) );
        ASSERT_TRUE( faceChannel == NULL );
    }
}

TEST( FractionalInterface, VertexAdjacency ) {

    const size_t numMeshes = 4;

    std::unique_ptr<mesh_interface> meshes[numMeshes];
    const double fraction = 1;

    // Simple mesh, just a single triangle
    meshes[0] = polymesh3_interface::create_const_instance( make_triangle_polymesh() );

    // A cube mesh with 6 quads
    meshes[1] = polymesh3_interface::create_const_instance( make_cube_polymesh() );

    // A cube mesh with 4 quads and two distinct holes
    std::set<size_t> holes;
    holes.insert( 0 );
    holes.insert( 1 );
    meshes[2] = polymesh3_interface::create_const_instance( make_cube_polymesh( holes ) );

    // A cube mesh with 4 quads and one hole
    holes.clear();
    holes.insert( 0 );
    holes.insert( 2 );
    meshes[3] = polymesh3_interface::create_const_instance( make_cube_polymesh( holes ) );

    for( size_t meshId = 0; meshId < numMeshes; ++meshId ) {
        mesh_interface_ptr mesh = mesh_interface_ptr( std::unique_ptr<fractional_vertex_interface>(
            new fractional_vertex_interface( meshes[meshId].get(), fraction ) ) );

        mesh->init_adjacency();

        for( size_t i = 0; i < mesh->get_num_verts(); ++i ) {
            vertex_iterator vIt;
            EXPECT_FALSE( mesh->init_vertex_iterator( vIt, i ) );
        }
    }
}

TEST( FractionalInterface, AddingChannels ) {
    polymesh3_ptr m = create_test_mesh();
    std::unique_ptr<polymesh3_interface> poly = polymesh3_interface::create_instance( m );
    boost::shared_ptr<fractional_face_interface> face_mesh =
        boost::shared_ptr<fractional_face_interface>( new fractional_face_interface( poly.get(), 1 ) );
    boost::shared_ptr<fractional_vertex_interface> vert_mesh =
        boost::shared_ptr<fractional_vertex_interface>( new fractional_vertex_interface( poly.get(), 1 ) );

    ASSERT_TRUE( face_mesh.get() != NULL );
    ASSERT_TRUE( vert_mesh.get() != NULL );

    EXPECT_FALSE( face_mesh->request_channel( _T( "newFace" ), false, false, false ) );
    EXPECT_THROW( face_mesh->request_channel( _T( "newFace" ), false, false ), std::runtime_error );

    poly->add_face_channel( _T( "newFace" ), frantic::channels::data_type_int32, 1 );

    EXPECT_FALSE( face_mesh->get_face_channels().has_channel( _T( "newFace" ) ) );
    EXPECT_FALSE( vert_mesh->get_face_channels().has_channel( _T( "newFace" ) ) );

    EXPECT_TRUE( face_mesh->request_channel( _T( "newFace" ), false, false, false ) );
    EXPECT_FALSE( face_mesh->get_face_channels().has_channel( _T( "newFace" ) ) );

    EXPECT_FALSE( vert_mesh->request_channel( _T( "newFace" ), false, false, false ) );
    EXPECT_THROW( vert_mesh->request_channel( _T( "newFace" ), false, false ), std::runtime_error );

    EXPECT_FALSE( vert_mesh->request_channel( _T( "newVert" ), true, false, false ) );
    EXPECT_THROW( vert_mesh->request_channel( _T( "newVert" ), true, false ), std::runtime_error );
    EXPECT_FALSE( face_mesh->request_channel( _T( "newVert" ), true, false, false ) );
    EXPECT_THROW( face_mesh->request_channel( _T( "newVert" ), true, false ), std::runtime_error );

    poly->add_vertex_channel( _T( "newVert" ), frantic::channels::data_type_int32, 1 );

    EXPECT_TRUE( face_mesh->request_channel( _T( "newVert" ), true, false, false ) );
    EXPECT_FALSE( face_mesh->get_vertex_channels().has_channel( _T( "newVert" ) ) );

    EXPECT_TRUE( vert_mesh->request_channel( _T( "newVert" ), true, false, false ) );
    EXPECT_FALSE( vert_mesh->get_vertex_channels().has_channel( _T( "newVert" ) ) );

    face_mesh->reset_mesh();
    EXPECT_TRUE( face_mesh->get_vertex_channels().has_channel( _T( "newVert" ) ) );
    EXPECT_TRUE( face_mesh->get_face_channels().has_channel( _T( "newFace" ) ) );

    vert_mesh->reset_mesh();
    EXPECT_TRUE( vert_mesh->get_vertex_channels().has_channel( _T( "newVert" ) ) );
    EXPECT_FALSE( vert_mesh->get_face_channels().has_channel( _T( "newFace" ) ) );
}
