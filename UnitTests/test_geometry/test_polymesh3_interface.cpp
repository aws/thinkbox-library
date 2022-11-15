// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"

#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif

namespace {

frantic::geometry::polymesh3_ptr create_test_mesh() {
    // Create the polymesh3 to test
    frantic::geometry::polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );
    builder.add_vertex( 0, 1, 0 );
    int faceIndices[] = { 0, 1, 2, 3 };
    builder.add_polygon( faceIndices, 4 );

    frantic::geometry::polymesh3_ptr polymesh = builder.finalize();
    EXPECT_TRUE( bool( polymesh ) );

    polymesh->add_empty_vertex_channel( _T( "simple" ), frantic::channels::data_type_int32, 1 );
    frantic::geometry::polymesh3_vertex_accessor<boost::int32_t> simpleAcc(
        polymesh->get_vertex_accessor<boost::int32_t>( _T( "simple" ) ) );
    simpleAcc.get_vertex( 0 ) = 4;
    simpleAcc.get_vertex( 1 ) = 5;
    simpleAcc.get_vertex( 2 ) = 6;
    simpleAcc.get_vertex( 3 ) = 7;

    polymesh->add_empty_vertex_channel( _T( "custom" ), frantic::channels::data_type_float32, 3, 2 );
    frantic::geometry::polymesh3_vertex_accessor<frantic::graphics::vector3f> customAcc(
        polymesh->get_vertex_accessor<frantic::graphics::vector3f>( _T( "custom" ) ) );
    customAcc.get_vertex( 0 ).set( 2, 3, 4 );
    customAcc.get_vertex( 1 ).set( 5, 6, 7 );
    customAcc.get_face( 0 ).first[0] = 0;
    customAcc.get_face( 0 ).first[1] = 0;
    customAcc.get_face( 0 ).first[2] = 1;
    customAcc.get_face( 0 ).first[3] = 1;

    polymesh->add_empty_face_channel( _T( "face" ), frantic::channels::data_type_uint8, 1 );
    frantic::geometry::polymesh3_face_accessor<boost::uint8_t> faceAcc(
        polymesh->get_face_accessor<boost::uint8_t>( _T( "face" ) ) );
    faceAcc.get_face( 0 ) = 12;

    return polymesh;
}

frantic::geometry::mesh_interface_ptr create_test_mesh_interface( const std::string& type ) {
    using namespace frantic::geometry;

    polymesh3_ptr mesh = create_test_mesh();

    if( type == "ReadOnly" ) {
        return mesh_interface_ptr( polymesh3_interface::create_const_instance( mesh ) );
    } else if( type == "Writable" ) {
        return mesh_interface_ptr( polymesh3_interface::create_instance( mesh ) );
    } else {
        throw std::runtime_error( "create_test_mesh_interface Error: unknown type: \"" + type + "\"" );
    }
}

} // anonymous namespace

// Test behavior that is common to both ReadOnly and Writable interfaces
class Polymesh3InterfaceCommon : public ::testing::TestWithParam<std::string> {};

TEST_P( Polymesh3InterfaceCommon, Common ) {
    using namespace std;
    using namespace boost;

    frantic::geometry::mesh_interface_ptr mesh = create_test_mesh_interface( GetParam() );
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

    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), mesh->get_vert( 3 ) );

    // get_num_faces()
    EXPECT_EQ( mesh->get_num_faces(), 1 );

    // get_num_face_verts()
    EXPECT_EQ( mesh->get_num_face_verts( 0 ), 4 );

    // get_face_vert_index()
    EXPECT_EQ( mesh->get_face_vert_index( 0, 0 ), 0 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 1 ), 1 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 2 ), 2 );
    EXPECT_EQ( mesh->get_face_vert_index( 0, 3 ), 3 );

    // get_face_vert_indices()
    std::size_t faceVertIndices[4];
    mesh->get_face_vert_indices( 0, faceVertIndices );
    EXPECT_EQ( faceVertIndices[0], 0 );
    EXPECT_EQ( faceVertIndices[1], 1 );
    EXPECT_EQ( faceVertIndices[2], 2 );
    EXPECT_EQ( faceVertIndices[3], 3 );

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
        frantic::geometry::mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
        std::set<frantic::tstring> vertexChannelNames;
        for( frantic::geometry::mesh_interface::mesh_channel_map::const_iterator i( vertexChannelMap.begin() ),
             ie( vertexChannelMap.end() );
             i != ie; ++i ) {
            vertexChannelNames.insert( i->first );
        }
        EXPECT_EQ( vertexChannelNames.size(), 2 );
        EXPECT_EQ( vertexChannelNames.count( _T("simple") ), 1 );
        EXPECT_EQ( vertexChannelNames.count( _T("custom") ), 1 );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T("simple") ) );
        EXPECT_TRUE( vertexChannelMap.has_channel( _T("custom") ) );

        { // scope for simpleChannel
            const frantic::geometry::mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T("simple") );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_EQ( simpleChannel->get_name(), _T("simple") );
            EXPECT_EQ( simpleChannel->get_channel_type(), frantic::geometry::mesh_channel::vertex );
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
            const frantic::geometry::mesh_channel* customChannel = vertexChannelMap.get_channel( _T("custom") );
            ASSERT_TRUE( customChannel != NULL );
            EXPECT_EQ( customChannel->get_name(), _T("custom") );
            EXPECT_EQ( customChannel->get_channel_type(), frantic::geometry::mesh_channel::face_vertex );
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
        frantic::geometry::mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
        std::set<frantic::tstring> faceChannelNames;
        for( frantic::geometry::mesh_interface::mesh_channel_map::const_iterator i( faceChannelMap.begin() ),
             ie( faceChannelMap.end() );
             i != ie; ++i ) {
            faceChannelNames.insert( i->first );
        }
        EXPECT_EQ( faceChannelNames.size(), 1 );
        EXPECT_EQ( faceChannelNames.count( _T("face") ), 1 );
        EXPECT_TRUE( faceChannelMap.has_channel( _T("face") ) );

        const frantic::geometry::mesh_channel* faceChannel = faceChannelMap.get_channel( _T("face") );
        ASSERT_TRUE( faceChannel != NULL );
        EXPECT_EQ( faceChannel->get_name(), _T("face") );
        EXPECT_EQ( faceChannel->get_channel_type(), frantic::geometry::mesh_channel::face );
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

INSTANTIATE_TEST_CASE_P( Polymesh3InterfaceCommon, Polymesh3InterfaceCommon,
                         ::testing::Values( "ReadOnly", "Writable" ) );

TEST( Polymesh3Interface, ReadOnly ) {
    using namespace std;
    using namespace boost;

    frantic::geometry::mesh_interface_ptr mesh = create_test_mesh_interface( "ReadOnly" );
    ASSERT_TRUE( mesh.get() != NULL );

    // is_read_only()
    EXPECT_TRUE( mesh->is_read_only() );

    // set_vert()
    frantic::graphics::vector3f vert;
    vert.set( -1 );
    EXPECT_ANY_THROW( mesh->set_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) ) );
    // should keep its original value
    vert.set( 0 );
    mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), vert );

    { // scope for vertexChannelMap
        frantic::geometry::mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();

        { // scope for simpleChannel
            const frantic::geometry::mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T("simple") );
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
            const frantic::geometry::mesh_channel* customChannel = vertexChannelMap.get_channel( _T("custom") );
            ASSERT_TRUE( customChannel != NULL );
            EXPECT_FALSE( customChannel->is_writeable() );

            frantic::graphics::vector3f v;

            // should keep its original value after set_value()
            customChannel->set_value( 0, &( frantic::graphics::vector3f( 0, 1, 2 )[0] ) );
            customChannel->get_value( 0, &v[0] );
            EXPECT_EQ( v, frantic::graphics::vector3f( 2, 3, 4 ) );

            // should keep its original value after set_fv_index()
            EXPECT_ANY_THROW( customChannel->set_fv_index( 0, 3, 0 ) );
            EXPECT_EQ( customChannel->get_fv_index( 0, 3 ), 1 );
        }
    }

    { // scope for faceChannelMap
        frantic::geometry::mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();

        const frantic::geometry::mesh_channel* faceChannel = faceChannelMap.get_channel( _T("face") );
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

TEST( Polymesh3Interface, Writable ) {
    using namespace std;
    using namespace boost;

    frantic::geometry::mesh_interface_ptr mesh = create_test_mesh_interface( "Writable" );
    ASSERT_TRUE( mesh.get() != NULL );

    // is_read_only()
    EXPECT_FALSE( mesh->is_read_only() );

    // set_vert()
    frantic::graphics::vector3f vert;
    vert.set( -1 );
    mesh->set_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    vert.set( 0 );
    mesh->get_vert( 3, *reinterpret_cast<float( * )[3]>( &vert[0] ) );
    EXPECT_EQ( frantic::graphics::vector3f( -1 ), vert );

    mesh->set_vert( 3, frantic::graphics::vector3f( 0, 1, 0 ) );
    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), mesh->get_vert( 3 ) );

    { // scope for vertexChannelMap
        frantic::geometry::mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();

        { // scope for simpleChannel
            const frantic::geometry::mesh_channel* simpleChannel = vertexChannelMap.get_channel( _T("simple") );
            ASSERT_TRUE( simpleChannel != NULL );
            EXPECT_TRUE( simpleChannel->is_writeable() );

            boost::int32_t i;

            const boost::int32_t eight = 8;
            simpleChannel->set_value( 2, &eight );
            simpleChannel->get_value( 2, &i );
            EXPECT_EQ( i, 8 );
        }

        { // scope for customChannel
            const frantic::geometry::mesh_channel* customChannel = vertexChannelMap.get_channel( _T("custom") );
            ASSERT_TRUE( customChannel != NULL );
            EXPECT_TRUE( customChannel->is_writeable() );

            frantic::graphics::vector3f v;

            customChannel->set_value( 0, &( frantic::graphics::vector3f( 0, 1, 2 )[0] ) );
            customChannel->get_value( 0, &v[0] );
            EXPECT_EQ( v, frantic::graphics::vector3f( 0, 1, 2 ) );

            customChannel->set_fv_index( 0, 3, 0 );
            EXPECT_EQ( customChannel->get_fv_index( 0, 3 ), 0 );
        }
    }

    { // scope for faceChannelMap
        frantic::geometry::mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();

        const frantic::geometry::mesh_channel* faceChannel = faceChannelMap.get_channel( _T("face") );
        ASSERT_TRUE( faceChannel != NULL );
        EXPECT_TRUE( faceChannel->is_writeable() );

        boost::uint8_t i;

        const boost::uint8_t eight = 8;
        faceChannel->set_value( 0, &eight );
        faceChannel->get_value( 0, &i );
        EXPECT_EQ( i, 8 );
    }
}

TEST( Polymesh3Interface, Adjacency ) {
    using namespace frantic::geometry;

    const size_t numMeshes = 4;

    std::unique_ptr<polymesh3_interface> meshes[numMeshes];

    // Simple mesh, just a single triangle
    meshes[0] = polymesh3_interface::create_instance( make_triangle_polymesh() );

    // A cube mesh with 6 quads
    meshes[1] = polymesh3_interface::create_instance( make_cube_polymesh() );

    // A cube mesh with 4 quads and two distinct holes
    std::set<size_t> holes;
    holes.insert( 0 );
    holes.insert( 1 );
    meshes[2] = polymesh3_interface::create_instance( make_cube_polymesh( holes ) );

    // A cube mesh with 4 quads and one hole
    holes.clear();
    holes.insert( 0 );
    holes.insert( 2 );
    meshes[3] = polymesh3_interface::create_instance( make_cube_polymesh( holes ) );

    for( size_t meshId = 0; meshId < numMeshes; ++meshId ) {
        std::unique_ptr<polymesh3_interface> mesh = std::move( meshes[meshId] );

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
            mesh->init_vertex_iterator( vIt, i );

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
