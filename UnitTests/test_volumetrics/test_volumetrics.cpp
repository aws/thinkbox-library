// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/assign/std/vector.hpp>

#include <boost/assign/list_of.hpp>

#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/strings/tstring.hpp>

#include <frantic/graphics/boundbox3.hpp>

#include <frantic/geometry/trimesh3_degeneracy_removal.hpp>
#include <frantic/geometry/trimesh3_file_io.hpp>
#include <frantic/misc/functor.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_rle_level_set.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>
#include <frantic/volumetrics/levelset/geometry_to_levelset.hpp>
#include <frantic/volumetrics/levelset/rle_level_set_file_io.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

#include <frantic/files/files.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>

using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace frantic;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::fluids;
using namespace frantic::files;
using namespace frantic::strings;
using namespace frantic::volumetrics::levelset;
using namespace frantic::volumetrics::implicitsurface;
using namespace frantic::volumetrics;
using namespace frantic::particles;
using namespace frantic::channels;

static void remove_out_of_bounds_faces( trimesh3& mcMesh ) {
    std::vector<char> facesArr;
    std::vector<bool> vertMask;
    const std::vector<vector3f>& vertices = mcMesh.vertices_ref();
    boundbox3f bounds( vector3f( 0.5f ), vector3f( 1.5f ) );
    for( size_t i = 0; i < mcMesh.vertex_count(); ++i ) {
        vertMask.push_back( bounds.contains( vertices[i] ) );
    }

    const std::vector<vector3>& faces = mcMesh.faces_ref();

    for( size_t i = 0; i < mcMesh.face_count(); ++i ) {
        if( !vertMask[faces[i][0]] || !vertMask[faces[i][1]] || !vertMask[faces[i][2]] )
            facesArr.push_back( char( false ) );
        else
            facesArr.push_back( char( true ) );
    }
    remove_faces( mcMesh, frantic::make_bool_array_functor( facesArr ) );
    remove_dead_vertices( mcMesh );
}

TEST( Volumetrics, VoxelCoordSystem ) {

    voxel_coord_system vcs( vector3f( 7, -4.2f, 12.f ), 0.76f );

    transform4f identity = vcs.voxel_to_world_transform() * vcs.world_to_voxel_transform();

    ASSERT_TRUE( identity.is_identity( 0.0001f ) );
}

TEST( Volumetrics, CachedTrilerp ) {

    const voxel_coord_system vcs;

    std::vector<vector3> coordinateArray;
    coordinateArray += vector3( 1, 1, 1 ), vector3( 2, 1, 1 ), vector3( 1, 2, 1 ), vector3( 2, 2, 1 );
    coordinateArray += vector3( 1, 1, 2 ), vector3( 2, 1, 2 ), vector3( 1, 2, 2 ), vector3( 2, 2, 2 );

    rle_index_spec ris;
    ris.build_from_voxel_array( coordinateArray );

    std::vector<float> distanceData;
    distanceData += 0.f, 1.f, 2.f, 4.f, 8.f, 16.f, 32.f, 64.f;

    rle_level_set ls( vcs, ris, distanceData, 7.f, 7.f );

    ls.add_channel<float>( _T("float") );
    rle_channel_accessor<float> floatAcc = ls.get_channel_accessor<float>( _T("float") );
    rle_channel_general_accessor generalFloatAcc = ls.get_channel_general_accessor( _T("float") );
    for( size_t i = 0; i < 8; ++i ) {
        floatAcc[i] = distanceData[i];
    }

    ls.add_channel<vector3f>( _T("vector3f") );
    rle_channel_accessor<vector3f> vecAcc = ls.get_channel_accessor<vector3f>( _T("vector3f") );
    rle_channel_general_accessor generalVecAcc = ls.get_channel_general_accessor( _T("vector3f") );
    for( size_t i = 0; i < 8; ++i ) {
        vecAcc[i] = vector3f( distanceData[i], 0.5f * distanceData[i], 2.f * distanceData[i] );
    }

    cached_trilerp trilerp( ris );

    bool success = false;
    float floatOut = -1.f;
    vector3f vecOut = vector3f( -1.f );

    // float lookups
    success = trilerp.get( &floatAcc[0], vector3f( 1.5f ), floatOut );
    EXPECT_TRUE( success );
    EXPECT_EQ( floatOut, 0.f );

    success = trilerp.get( &floatAcc[0], vector3f( 3.f ), floatOut ); // not fully defined
    EXPECT_TRUE( success );
    EXPECT_EQ( floatOut, 64.f );

    success = trilerp.get( floatAcc, vector3f( 2.f, 1.5f, 1.5f ), floatOut );
    EXPECT_TRUE( success );
    EXPECT_FLOAT_EQ( floatOut, 0.5f );

    success = trilerp.get( const_rle_channel_accessor<float>( floatAcc ), vector3f( 3.f, 2.f, 2.f ),
                           floatOut ); // not fully defined
    EXPECT_TRUE( success );
    EXPECT_FLOAT_EQ( floatOut, 21.25f );

    success = trilerp.get( generalFloatAcc, vector3f( 1.5f, 2.f, 1.5f ), reinterpret_cast<char*>( &floatOut ) );
    EXPECT_TRUE( success );
    EXPECT_FLOAT_EQ( floatOut, 1.f );

    success = trilerp.get( generalFloatAcc, vector3f( 5.f ), reinterpret_cast<char*>( &floatOut ) );
    ASSERT_FALSE( success );

    // level set lookups
    floatOut = trilerp.get( ls, vector3f( 2.5f, 2.5f, 1.5f ) );
    EXPECT_FLOAT_EQ( floatOut, 4.f );

    floatOut = trilerp.get( ls, vector3f( 2.f, 2.f, 3.f ) ); // partially outside
    EXPECT_FLOAT_EQ( floatOut, 30.f );

    floatOut = trilerp.get( ls, vector3f( -1.f ) ); // completely outside
    EXPECT_FLOAT_EQ( floatOut, 8.f );

    // vector3f lookups
    success = trilerp.get( &vecAcc[0], vector3f( 1.5f, 1.5f, 2.5f ), vecOut );
    EXPECT_TRUE( success );
    EXPECT_EQ( vecOut, vector3f( 8.f, 4.f, 16.f ) );

    success = trilerp.get( vecAcc, vector3f( 2.f, 1.5f, 2.5f ), vecOut );
    EXPECT_TRUE( success );
    EXPECT_EQ( vecOut, vector3f( 12.f, 6.f, 24.f ) );

    success = trilerp.get( generalVecAcc, vector3f( 1.5f, 2.f, 2.5f ), reinterpret_cast<char*>( &vecOut ) );
    EXPECT_TRUE( success );
    EXPECT_EQ( vecOut, vector3f( 20.f, 10.f, 40.f ) );

    success = trilerp.get( generalVecAcc, vector3f( 3.f, 1.f, 2.f ),
                           reinterpret_cast<char*>( &vecOut ) ); // not fully defined
    EXPECT_TRUE( success );
    EXPECT_EQ( vecOut, vector3f( 8.5f, 4.25f, 17.f ) );
}

TEST( Volumetrics, RleIndexFill ) {

    voxel_coord_system vcs( vector3f(), 0.5f );
    boundbox3 voxelBounds( vector3( 0 ), vector3( 4, 4, 6 ) );
    boundbox3f baseBox( vector3f( -1 ), vector3f( 4.f ) );
    trimesh3 box;
    box.set_to_box( baseBox );

    rle_level_set rls( vcs );

    convert_geometry_to_levelset( box, -1.f, 1.f, voxelBounds, rls );

    EXPECT_TRUE( rls.get_rle_index_spec().check_consistency( cout ) ) << "Consistency check failed.";

    rle_index_spec velSpec;
    FF_LOG( debug ) << "Original # of defined voxels=" << rls.size() << endl;
    velSpec.build_by_filling( rls.get_rle_index_spec() );
    FF_LOG( debug ) << "Filled # of defined voxels=" << velSpec.data_size() << endl;

    EXPECT_TRUE( velSpec.check_consistency( cout ) ) << "Consistency check failed";
}

TEST( Volumetrics, RleVoxelFieldIO ) {
    using namespace frantic::files;

    boost::filesystem::path tempPath = boost::filesystem::temp_directory_path();
    frantic::tstring filename = to_tstring( tempPath / "unit_voxel_field_test.rls" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( filename );

    boundbox3 bounds( vector3(), vector3( 100 ) );
    voxel_coord_system vcs( vector3f(), 1.f );

    rle_index_spec ris;
    ris.build_from_random( bounds, 4.f );
    rle_voxel_field voxelField( vcs, ris );

    voxelField.add_channel<vector3f>( _T("Velocity") );
    voxelField.add_channel<int>( _T("Fun") );

    FF_LOG( debug ) << "VoxelCount:" << voxelField.size() << endl;

    rle_channel_accessor<int> fun = voxelField.get_channel_accessor<int>( _T("Fun") );
    rle_channel_accessor<vector3f> vel = voxelField.get_channel_accessor<vector3f>( _T("Velocity") );

    for( rle_defined_iterator i = ris.begin(), ie = ris.end(); i != ie; ++i ) {
        fun[i.get_data_index()] = rand();
        vel[i.get_data_index()] = vector3f::from_random();
    }

    levelset::write_rle_voxel_field_file( filename, voxelField );

    rle_voxel_field newVoxelField( vcs );
    levelset::read_rle_voxel_field_file( filename, newVoxelField );

    EXPECT_TRUE( newVoxelField.size() == voxelField.size() ) << "Voxel Field Sizes don't match";

    rle_channel_accessor<int> newFun = newVoxelField.get_channel_accessor<int>( _T("Fun") );
    rle_channel_accessor<vector3f> newVel = newVoxelField.get_channel_accessor<vector3f>( _T("Velocity") );

    for( rle_defined_iterator i = ris.begin(), ie = ris.end(); i != ie; ++i ) {
        int index = i.get_data_index();
        EXPECT_TRUE( fun[index] == newFun[index] ) << "Channel Fun doesn't match";

        EXPECT_TRUE( vel[index] == newVel[index] ) << "Channel Velocity doesn't match";
    }
}

class VolumetricsMarchingCubes : public ::testing::Test {
  protected:
    void SetUp() {
        singleCube.push_back( vector3( 0, 0, 0 ) ); // 0
        singleCube.push_back( vector3( 1, 0, 0 ) ); // 1
        singleCube.push_back( vector3( 0, 1, 0 ) ); // 2
        singleCube.push_back( vector3( 1, 1, 0 ) ); // 3
        singleCube.push_back( vector3( 0, 0, 1 ) ); // 4
        singleCube.push_back( vector3( 1, 0, 1 ) ); // 5
        singleCube.push_back( vector3( 0, 1, 1 ) ); // 6
        singleCube.push_back( vector3( 1, 1, 1 ) ); // 7
        ris.build_from_voxel_array( singleCube );
        data.resize( 8 );
        rls = rle_level_set( vcs, ris, data, 2.f, 2.f );
    }

    voxel_coord_system vcs;
    rle_index_spec ris;
    vector<vector3> singleCube;
    trimesh3 mcMesh;
    vector<float> data;
    rle_level_set rls;
};

// Go through all the 14 unique cube cases to confirm it's generating something reasonable

TEST_F( VolumetricsMarchingCubes, CubeCase1 ) {

    // Pick the cube case where everything is outside except [0,0,0] (cube case 1)
    rls[0] = -1;
    rls[1] = 1;
    rls[2] = 1;
    rls[3] = 1;
    rls[4] = 1;
    rls[5] = 1;
    rls[6] = 1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 3 );
    EXPECT_EQ( mcMesh.face_count(), 1 );
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( 1, 1, 1 ) ) > 0 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase1Inverse ) {

    // Now try an inverse cube of that with a rotation, where everything is inside except [1,0,0] (cube case 1)
    rls[0] = -1;
    rls[1] = 1;
    rls[2] = -1;
    rls[3] = -1;
    rls[4] = -1;
    rls[5] = -1;
    rls[6] = -1;
    rls[7] = -1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 3 );
    ASSERT_EQ( mcMesh.face_count(), 1 );
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( 1, -1, -1 ) ) > 0 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase2 ) {

    // Now try a cube case where just two adjacent samples are inside, [1,0,0] and [1,1,0] (cube case 2)
    rls[0] = 1;
    rls[1] = -1;
    rls[2] = 1;
    rls[3] = -1;
    rls[4] = 1;
    rls[5] = 1;
    rls[6] = 1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 4 );
    EXPECT_EQ( mcMesh.face_count(), 2 );
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( -1, 0, 1 ) ) > 0 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase3 ) {

    // Now try a cube case where just two samples opposite on one face are inside, [0,1,1] and [1,1,0] (cube case 3)
    rls[0] = 1;
    rls[1] = 1;
    rls[2] = 1;
    rls[3] = -1;
    rls[4] = 1;
    rls[5] = 1;
    rls[6] = -1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 6 );
    EXPECT_TRUE( mcMesh.face_count() == 2 || mcMesh.face_count() == 4 ); // this case is ambiguous
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( 0, -1, 0 ) ) > 0 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase4 ) {

    // Now try a cube case where just two samples opposite on the cube are inside, [0,0,1] and [1,1,0] (cube case 4)
    rls[0] = 1;
    rls[1] = 1;
    rls[2] = 1;
    rls[3] = -1;
    rls[4] = -1;
    rls[5] = 1;
    rls[6] = 1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct
    EXPECT_EQ( mcMesh.vertex_count(), 6 );
    EXPECT_TRUE( mcMesh.face_count() ==
                 2 ); // this case could be ambiguous, but we choose not to connect the opposite corners
}

TEST_F( VolumetricsMarchingCubes, CubeCase5 ) {

    // Now try a cube case where three out of four corners on one cube face are inside, [0,0,1], [1,0,1] and [0,1,1]  (
    // cube case 5 )
    rls[0] = 1;
    rls[1] = 1;
    rls[2] = 1;
    rls[3] = 1;
    rls[4] = -1;
    rls[5] = -1;
    rls[6] = -1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 5 );
    EXPECT_EQ( mcMesh.face_count(), 3 );
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( 0.5f, 0.5f, -1 ) ) > 0 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase6 ) {

    // Now try cube case 6, [0,0,1], [0,1,1] and [1,1,0] (cube case 6)
    rls[0] = 1;
    rls[1] = 1;
    rls[2] = 1;
    rls[3] = -1;
    rls[4] = -1;
    rls[5] = 1;
    rls[6] = -1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 7 );
    EXPECT_TRUE( mcMesh.face_count() == 3 || mcMesh.face_count() == 5 ); // this case is ambiguous
}

TEST_F( VolumetricsMarchingCubes, CubeCase7 ) {

    // Now try cube case 7, [1,0,0], [0,1,0] and [1,1,1] (cube case 7)
    rls[0] = 1;
    rls[1] = -1;
    rls[2] = -1;
    rls[3] = 1;
    rls[4] = 1;
    rls[5] = 1;
    rls[6] = 1;
    rls[7] = -1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct
    EXPECT_EQ( mcMesh.vertex_count(), 9 );
    EXPECT_TRUE( mcMesh.face_count() == 3 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase8 ) {

    // Now try a cube case where all the corners on one face are inside, [0,0,0], [0,0,1], [0,1,0] and [0,1,1] (cube
    // case 8)
    rls[0] = -1;
    rls[1] = 1;
    rls[2] = -1;
    rls[3] = 1;
    rls[4] = -1;
    rls[5] = 1;
    rls[6] = -1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 4 );
    EXPECT_EQ( mcMesh.face_count(), 2 );
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( 1, 0, 0 ) ) > 0 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase9 ) {

    // Now try cube case 9, [0,1,0], [1,0,0], [1,1,0] and [1,1,1] (cube case 9)
    rls[0] = 1;
    rls[1] = -1;
    rls[2] = -1;
    rls[3] = -1;
    rls[4] = 1;
    rls[5] = 1;
    rls[6] = 1;
    rls[7] = -1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 6 );
    EXPECT_EQ( mcMesh.face_count(), 4 );
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( -1, -1, 1 ) ) > 0 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase10 ) {

    // Now try cube case 10, [1,0,0], [1,0,1], [0,1,0] and [0,1,1] (cube case 10)
    rls[0] = 1;
    rls[1] = -1;
    rls[2] = -1;
    rls[3] = 1;
    rls[4] = 1;
    rls[5] = -1;
    rls[6] = -1;
    rls[7] = 1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct
    EXPECT_EQ( mcMesh.vertex_count(), 8 );
    EXPECT_TRUE( mcMesh.face_count() == 4 ); // This case is ambiguous, but both have 4 faces.
}

TEST_F( VolumetricsMarchingCubes, CubeCase11 ) {

    // Now try cube case 11, [0,0,0], [1,0,1], [1,1,0] and [1,1,1] (cube case 11)
    rls[0] = -1;
    rls[1] = 1;
    rls[2] = 1;
    rls[3] = -1;
    rls[4] = 1;
    rls[5] = -1;
    rls[6] = 1;
    rls[7] = -1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct
    EXPECT_EQ( mcMesh.vertex_count(), 8 );
    EXPECT_TRUE( mcMesh.face_count() == 4 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase12 ) {

    // Now try cube case 12, [1,0,0], [0,1,0], [0,0,1] and [1,1,1] (cube case 12)
    rls[0] = 1;
    rls[1] = -1;
    rls[2] = -1;
    rls[3] = 1;
    rls[4] = -1;
    rls[5] = 1;
    rls[6] = 1;
    rls[7] = -1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct
    EXPECT_EQ( mcMesh.vertex_count(), 12 );
    EXPECT_TRUE( mcMesh.face_count() == 4 );
}

TEST_F( VolumetricsMarchingCubes, CubeCase13 ) {

    // Now try cube case 13, [1,0,0], [1,1,0], [0,1,1] and [1,1,1] (cube case 13)
    rls[0] = 1;
    rls[1] = -1;
    rls[2] = 1;
    rls[3] = -1;
    rls[4] = 1;
    rls[5] = 1;
    rls[6] = -1;
    rls[7] = -1;
    convert_levelset_to_trimesh3( rls, mcMesh );
    remove_out_of_bounds_faces( mcMesh );
    // Check that the face and vert counts are correct, and that the face normals are pointing in the expected
    // direction.
    EXPECT_EQ( mcMesh.vertex_count(), 6 );
    EXPECT_TRUE( mcMesh.face_count() == 4 );
    for( int i = 0; i < (int)mcMesh.face_count(); ++i )
        EXPECT_TRUE( vector3f::dot( mcMesh.compute_face_normal( i ), vector3f( -1, -1, 0 ) ) > 0 );
}

TEST( Volumetrics, BuildLevelSetFromDirectLinearRLSISPolicy ) {
    unsigned int seed = 12345;
    srand( seed );

    stringstream ss;

    // set test bounds for the random level sets
    int bound = 10;
    boundbox3 testBounds( -bound, bound, -bound, bound, -bound, bound );
    for( int c = 0; c < 10000; c++ ) {

        // create a random test level set
        rle_index_spec ris;
        std::vector<float> rlsData;
        ris.build_from_random( testBounds, 4, 2 );

        // generate random data for the index spec
        for( rle_defined_iterator i = ris.begin(); i != ris.end(); ++i )
            rlsData.push_back( 10.f * ( (float)rand() / RAND_MAX - 0.5f ) );

        ASSERT_TRUE( ris.check_consistency( ss ) )
            << ss.str() +
                   "testBuildLevelSetFromDirectLinearRLSISPolicy() - Random RLE index spec failed consistency check";

        float interfaceWidth = 2.0; // not really important for the purposes of this testing
        rle_level_set inRLS( voxel_coord_system(), ris, rlsData, interfaceWidth, interfaceWidth );
        inRLS.add_channel( _T("TestData"), 3, data_type_float32 );
        // generate some random channel data
        rle_channel_accessor<vector3f> channelAccessor = inRLS.get_channel_accessor<vector3f>( _T("TestData") );
        for( size_t i = 0; i < inRLS.size(); ++i )
            channelAccessor[i] = vector3f( 10.f * ( (float)rand() / RAND_MAX - 0.5f ) );

        direct_linear_rle_level_set_is_policy lisp( inRLS, 3 ); // -3 is specific to the ISP internals

        rle_level_set outRLS;
        build_rle_level_set( lisp, outRLS );

        // check for level set equality (vcs, index spec, distance data)
        ASSERT_TRUE( ( inRLS == outRLS ) ) << "input and output rle level sets are not equal";

        // check that the named channel data is equal
        ASSERT_TRUE( outRLS.has_channel( _T("TestData") ) ) << "output level set failed to copy the named channel";

        rle_channel_accessor<vector3f> inTestData = inRLS.get_channel_accessor<vector3f>( _T("TestData") );
        rle_channel_accessor<vector3f> outTestData = outRLS.get_channel_accessor<vector3f>( _T("TestData") );
        for( size_t i = 0; i < inRLS.size(); ++i ) {
            EXPECT_TRUE( inTestData[i].is_equal( outTestData[i] ) )
                << "The copied data for the named channel is not equal to the original data at location " +
                       boost::lexical_cast<std::string>( i );
        }
    }
}
