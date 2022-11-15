// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/assign/list_of.hpp>

#include "utilities/mesh_generators.hpp"
#include <frantic/files/filename_sequence.hpp>
#include <frantic/files/files.hpp>
#include <frantic/geometry/mesh_file_io_factory.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>

using namespace std;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::channels;

static const frantic::tstring& test_temp_dir() {
    static frantic::tstring _static_test_temp_dir = _T( "./XMeshSequenceSaverTestScratchDirectory/" );
    return _static_test_temp_dir;
}

static void clean_temp_directory() {
    std::vector<frantic::tstring> fileListing;
    frantic::files::get_filenames_in_directory( test_temp_dir(), fileListing );

    for( size_t i = 0; i < fileListing.size(); ++i ) {
        frantic::files::delete_file( test_temp_dir() + fileListing[i] );
    }
}

static bool is_equal_geometry( frantic::geometry::polymesh3_ptr a, frantic::geometry::polymesh3_ptr b ) {
    const std::size_t vertexCount = a->vertex_count();
    if( vertexCount != b->vertex_count() ) {
        return false;
    }

    for( std::size_t i = 0; i < vertexCount; ++i ) {
        if( a->get_vertex( i ) != b->get_vertex( i ) ) {
            return false;
        }
    }

    return true;
}

class GeometryTestFixture : public ::testing::Test {
  public:
    virtual void SetUp() {
        if( frantic::files::directory_exists( test_temp_dir() ) ) {
            clean_temp_directory();
        } else {
            frantic::files::make_directory( test_temp_dir() );
            cout << "directory made" << endl;
        }
    }

    virtual void TearDown() { clean_temp_directory(); }
};

TEST_F( GeometryTestFixture, MeshFileIoFactory ) {

    const frantic::tstring objFilename( test_temp_dir() + _T("test.obj") );
    const frantic::tstring xmeshFilename( test_temp_dir() + _T("test.xmesh") );

    frantic::geometry::polymesh3_ptr yUpMesh;
    {
        frantic::geometry::polymesh3_builder builder;
        builder.add_vertex( 0, 0, 0 );
        builder.add_vertex( 1, 0, 0 );
        builder.add_vertex( 1, 0, -1 );
        std::vector<int> indices = boost::assign::list_of( 0 )( 1 )( 2 );
        builder.add_polygon( indices );
        yUpMesh = builder.finalize();
    }

    frantic::geometry::polymesh3_ptr zUpMesh;
    {
        frantic::geometry::polymesh3_builder builder;
        builder.add_vertex( 0, 0, 0 );
        builder.add_vertex( 1, 0, 0 );
        builder.add_vertex( 1, 1, 0 );
        std::vector<int> indices = boost::assign::list_of( 0 )( 1 )( 2 );
        builder.add_polygon( indices );
        zUpMesh = builder.finalize();
    }

    frantic::geometry::mesh_file_io_factory factory;
    factory.set_coordinate_system( frantic::graphics::coordinate_system::right_handed_yup );
    factory.set_default_file_coordinate_system( frantic::graphics::coordinate_system::right_handed_zup );

    factory.write( objFilename, yUpMesh );
    factory.write( xmeshFilename, yUpMesh );

    // we should read the same yUpMesh back in again
    frantic::geometry::polymesh3_ptr mesh = factory.read_polymesh3( objFilename );
    ASSERT_TRUE( is_equal_geometry( mesh, yUpMesh ) );

    mesh = factory.read_polymesh3( xmeshFilename );
    ASSERT_TRUE( is_equal_geometry( mesh, yUpMesh ) );

    // change coordinate system to right-handed z-up
    // now we should get the zUpMesh
    factory.set_coordinate_system( frantic::graphics::coordinate_system::right_handed_zup );

    mesh = factory.read_polymesh3( objFilename );
    ASSERT_TRUE( is_equal_geometry( mesh, zUpMesh ) );

    mesh = factory.read_polymesh3( xmeshFilename );
    ASSERT_TRUE( is_equal_geometry( mesh, zUpMesh ) );
}

TEST_F( GeometryTestFixture, MeshInterfaceFileIoXmeshBoundbox ) {
    const frantic::tstring xmeshFilename( test_temp_dir() + _T("test.xmesh") );

    ASSERT_TRUE( !frantic::files::file_exists( xmeshFilename ) );

    frantic::geometry::trimesh3 mesh;
    mesh.set_to_box( frantic::graphics::boundbox3f( 1, 2, 3, 4, 5, 6 ) );

    frantic::geometry::mesh_interface::ptr_type meshInterface(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );
    frantic::logging::null_progress_logger progress;
    frantic::geometry::write_xmesh_mesh_file( xmeshFilename, meshInterface, progress );

    frantic::geometry::xmesh_metadata metadata;
    frantic::geometry::read_xmesh_metadata( xmeshFilename, metadata );

    ASSERT_TRUE( metadata.has_boundbox() );
    ASSERT_EQ( metadata.get_boundbox(), frantic::graphics::boundbox3f( 1, 2, 3, 4, 5, 6 ) );
}

TEST( GeometryTest, MeshInterfaceAddAndEraseCustomIndexVertexChannel ) {

    const frantic::tstring channelName = _T("MyChannelName");

    const size_t numInterfaces = 2;
    mesh_interface* interfaces[numInterfaces];

    std::unique_ptr<polymesh3_interface> iPolymesh = polymesh3_interface::create_instance( make_cube_polymesh() );
    trimesh3 triangleMesh;
    make_regular_tetrahedron( triangleMesh );
    std::unique_ptr<trimesh3_interface> iTrimesh = trimesh3_interface::create_instance( boost::move( triangleMesh ) );

    interfaces[0] = iPolymesh.get();
    interfaces[1] = iTrimesh.get();

    for( size_t iteration = 0; iteration < numInterfaces; ++iteration ) {
        mesh_interface* mesh = interfaces[iteration];

        mesh->add_vertex_channel_custom_faces( channelName, data_type_float32, 3, 15 );

        ASSERT_TRUE( mesh->has_vertex_channel( channelName ) ) << "Did not create vertex channel.";

        const mesh_channel* customChannel = mesh->get_vertex_channels().get_channel( channelName );

        ASSERT_TRUE( customChannel ) << "Returned a null channel.";

        size_t currentIndex = 0;

        for( size_t i = 0; i < mesh->get_num_faces(); ++i ) {
            for( size_t j = 0; j < mesh->get_num_face_verts( i ); ++j ) {
                customChannel->set_fv_index( i, j, currentIndex % 15 );
                ++currentIndex;
            }
        }

        currentIndex = 0;

        for( size_t i = 0; i < mesh->get_num_faces(); ++i ) {
            for( size_t j = 0; j < mesh->get_num_face_verts( i ); ++j ) {
                ASSERT_EQ( currentIndex % 15, customChannel->get_fv_index( i, j ) );
                ++currentIndex;
            }
        }

        mesh->erase_vertex_channel( channelName );

        ASSERT_TRUE( !mesh->has_vertex_channel( channelName ) ) << "Did not destroy vertex channel.";
        ASSERT_EQ( static_cast<mesh_channel*>( NULL ), mesh->get_vertex_channels().get_channel( channelName ) )
            << "Returned an erased channel.";

        // attempting to double remove should throw an error
        EXPECT_ANY_THROW( mesh->erase_vertex_channel( channelName ) )
            << "Attempting to erase a missing channel did not throw an error";
    }
}

TEST( GeometryTest, MeshInterfaceAddAndEraseFaceChannel ) {
    const frantic::tstring channelName = _T("MyChannelName");

    const size_t numInterfaces = 2;
    mesh_interface* interfaces[numInterfaces];

    std::unique_ptr<polymesh3_interface> iPolymesh = polymesh3_interface::create_instance( make_cube_polymesh() );
    trimesh3 triangleMesh;
    make_regular_tetrahedron( triangleMesh );
    std::unique_ptr<trimesh3_interface> iTrimesh = trimesh3_interface::create_instance( boost::move( triangleMesh ) );

    interfaces[0] = iPolymesh.get();
    interfaces[1] = iTrimesh.get();

    for( size_t iteration = 0; iteration < numInterfaces; ++iteration ) {
        mesh_interface* mesh = interfaces[iteration];

        mesh->add_face_channel( channelName, data_type_float32, 3 );

        ASSERT_TRUE( mesh->has_face_channel( channelName ) ) << "Did not create vertex channel.";
        const mesh_channel* customChannel = mesh->get_face_channels().get_channel( channelName );

        ASSERT_TRUE( customChannel ) << "Returned a null channel.";

        mesh->erase_face_channel( channelName );

        ASSERT_TRUE( !mesh->has_face_channel( channelName ) ) << "Did not destroy face channel.";
        ASSERT_EQ( static_cast<mesh_channel*>( NULL ), mesh->get_face_channels().get_channel( channelName ) )
            << "Returned an erased channel.";

        // attempting to double remove should throw an error
        EXPECT_ANY_THROW( mesh->erase_face_channel( channelName ) )
            << "Attempting to erase a missing channel did not throw an error";
    }
}
