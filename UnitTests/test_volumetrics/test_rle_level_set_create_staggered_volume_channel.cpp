// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/assign/std/vector.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>
#include <frantic/volumetrics/levelset/geometry_to_levelset.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

using namespace std;
using namespace boost;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::geometry;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::fluids;

TEST( RleLevelSetCreateStaggeredVolumeChannel, HalfOpen ) {
    // First build a level set to work with
    voxel_coord_system vcs( vector3f( 0.127f, 0.12f, 0.032f ), 0.478f );
    trimesh3 plane;
    boundbox3f worldBounds( -2, 2, -2, 2, -2, 2 );

    // First do the YZ plane, and a slightly offset X coordinate
    // float planeX = 1.67f;
    float planeX = -0.39f;
    plane.set_to_rectangular_grid( 2, 2, boundrect2f( -10, 10, -10, 10 ) );
    plane.transform( transform4f::from_y_rotation( float( M_PI ) / 2.f ) );
    plane.translate( vector3f( planeX, 0, 0 ) );
    rle_level_set rls( vcs );
    convert_geometry_to_levelset( plane, -2, 2, worldBounds, rls );

    // Now create the voxel field for the staggered velocities
    rle_index_spec ris;
    ris.build_by_filling( rls.get_rle_index_spec() );
    rle_voxel_field rvf( vcs, ris );

    // Create the staggered volume channel
    rls.create_volume_fraction_channels( rvf, _T("CntVol"), _T("StgVol") );

    rle_channel_accessor<float> cntVolAcc = rvf.get_channel_accessor<float>( _T("CntVol") );
    rle_channel_accessor<vector3f> stgVolAcc = rvf.get_channel_accessor<vector3f>( _T("StgVol") );

    // Go through and validate the volumes
    for( rle_defined_iterator i = rvf.get_rle_index_spec().begin(), ie = rvf.get_rle_index_spec().end(); i != ie;
         ++i ) {
        // FF_LOG(debug) << "Validating voxel coordinate " << i.get_coord() << endl;
        float cntVol = cntVolAcc[i.get_data_index()];
        vector3f stgVol = stgVolAcc[i.get_data_index()];
        float xDistanceCent = vcs.get_world_x_coord( i.get_coord().x + 0.5f ) - planeX,
              xDistanceStg = vcs.get_world_x_coord( (float)i.get_coord().x ) - planeX;
        float cntVolComputed;
        vector3f stgVolComputed;
        if( xDistanceStg < -vcs.voxel_length() / 2 )
            stgVolComputed.x = 1;
        else if( xDistanceStg < vcs.voxel_length() / 2 )
            stgVolComputed.x = 0.5f - xDistanceStg / vcs.voxel_length();
        else
            stgVolComputed.x = 0;

        if( xDistanceCent < -vcs.voxel_length() / 2 )
            cntVolComputed = stgVolComputed.y = stgVolComputed.z = 1;
        else if( xDistanceCent < vcs.voxel_length() / 2 )
            cntVolComputed = stgVolComputed.y = stgVolComputed.z = 0.5f - xDistanceCent / vcs.voxel_length();
        else
            cntVolComputed = stgVolComputed.y = stgVolComputed.z = 0;

        // Check that the staggered volumes from the create function match the ones calculated here
        //				FF_LOG(debug) << "vox coord: " << i.get_coord() << endl;
        EXPECT_NEAR( cntVol, cntVolComputed, 0.01f );
        EXPECT_NEAR( stgVol.x, stgVolComputed.x, 0.01f );
        EXPECT_NEAR( stgVol.y, stgVolComputed.y, 0.01f );
        EXPECT_NEAR( stgVol.z, stgVolComputed.z, 0.01f );
    }
}

TEST( RleLevelSetCreateStaggeredVolumeChannel, RotatedBox ) {
    ///////
    // This test checks the volumes of a rotated box
    ///////

    // First build a level set to work with
    voxel_coord_system vcs( vector3f( 0.127f, 0.12f, 0.032f ), 0.478f );
    boundbox3f baseBox( vector3f( -1.2f, -1.33f, -1.71f ), vector3f( 4.f, 3.37f, 2.9f ) );
    trimesh3 box;
    box.set_to_box( baseBox );
    box.transform( transform4f::from_euler_xyz_rotation( vector3f( 1.2f, 2.7f, 0.5f ) ) );
    rle_level_set rls( vcs );
    convert_geometry_to_levelset( box, -3, 3, rls );
    rls.dilate_defined_voxels( 2, _T("") );
    rls.reinitialize_signed_distance( _T("Populated") );
    rls.trim_to_populated( _T("Populated") );

    // Now create the voxel field for the staggered velocities
    rle_index_spec ris;
    ris.build_by_filling( rls.get_rle_index_spec() );
    rle_voxel_field rvf( vcs, ris );

    // Create the staggered volume channel
    rls.create_volume_fraction_channels( rvf, _T("CntVol"), _T("StgVol") );

    rle_channel_accessor<float> cntVolAcc = rvf.get_channel_accessor<float>( _T("CntVol") );
    rle_channel_accessor<vector3f> stgVolAcc = rvf.get_channel_accessor<vector3f>( _T("StgVol") );

    // Accumulate all the staggered volumes.  On completion, they should be equal.
    float totalCntVol = 0;
    vector3f totalStgVol( 0 );

    // Go through and validate the volumes
    for( rle_defined_iterator i = rvf.get_rle_index_spec().begin(), ie = rvf.get_rle_index_spec().end(); i != ie;
         ++i ) {
        vector3f voxCoordFlt( (float)i.get_coord().x, (float)i.get_coord().y, (float)i.get_coord().z );
        float cntVol = cntVolAcc[i.get_data_index()];
        vector3f stgVol = stgVolAcc[i.get_data_index()];

        float cntVolCalculated = 0;
        vector3f stgVolCalculated( 0 );
        // Calculate the centered volume estimate
        vector3f cubeCenter = vcs.get_world_coord( voxCoordFlt + vector3f( 0.5f, 0.5f, 0.5f ) );
        for( int c = 0; c < 2; ++c ) {
            for( int b = 0; b < 2; ++b ) {
                for( int a = 0; a < 2; ++a ) {
                    float sd = rls.trilerp_signed_distance( cubeCenter + 0.5f * vcs.voxel_length() *
                                                                             vector3f( a - 0.5f, b - 0.5f, c - 0.5f ) );
                    cntVolCalculated += .5f - .5f * max( -1.f, min( 1.f, 4 * sd / vcs.voxel_length() ) );
                }
            }
        }

        // Calculate the staggered volume using the linear interpolation function using the same formula used in the
        // create routine,
        // to validate that it's getting the linear interpolations right.
        cubeCenter = vcs.get_world_coord( voxCoordFlt + vector3f( 0, 0.5f, 0.5f ) );
        for( int c = 0; c < 2; ++c ) {
            for( int b = 0; b < 2; ++b ) {
                for( int a = 0; a < 2; ++a ) {
                    float sd = rls.trilerp_signed_distance( cubeCenter + 0.5f * vcs.voxel_length() *
                                                                             vector3f( a - 0.5f, b - 0.5f, c - 0.5f ) );
                    stgVolCalculated.x += .5f - .5f * max( -1.f, min( 1.f, 4 * sd / vcs.voxel_length() ) );
                }
            }
        }

        cubeCenter = vcs.get_world_coord( voxCoordFlt + vector3f( 0.5f, 0, 0.5f ) );
        for( int c = 0; c < 2; ++c ) {
            for( int b = 0; b < 2; ++b ) {
                for( int a = 0; a < 2; ++a ) {
                    float sd = rls.trilerp_signed_distance( cubeCenter + 0.5f * vcs.voxel_length() *
                                                                             vector3f( a - 0.5f, b - 0.5f, c - 0.5f ) );
                    stgVolCalculated.y += .5f - .5f * max( -1.f, min( 1.f, 4 * sd / vcs.voxel_length() ) );
                }
            }
        }

        cubeCenter = vcs.get_world_coord( voxCoordFlt + vector3f( 0.5f, 0.5f, 0 ) );
        for( int c = 0; c < 2; ++c ) {
            for( int b = 0; b < 2; ++b ) {
                for( int a = 0; a < 2; ++a ) {
                    float sd = rls.trilerp_signed_distance( cubeCenter + 0.5f * vcs.voxel_length() *
                                                                             vector3f( a - 0.5f, b - 0.5f, c - 0.5f ) );
                    stgVolCalculated.z += .5f - .5f * max( -1.f, min( 1.f, 4 * sd / vcs.voxel_length() ) );
                }
            }
        }

        cntVolCalculated *= 0.125f;
        stgVolCalculated *= 0.125f;

        // Check that the staggered volumes from the create function match the ones calculated here
        EXPECT_NEAR( cntVol, cntVolCalculated, 0.01f );
        EXPECT_NEAR( stgVol.x, stgVolCalculated.x, 0.01f );
        EXPECT_NEAR( stgVol.y, stgVolCalculated.y, 0.01f );
        EXPECT_NEAR( stgVol.z, stgVolCalculated.z, 0.01f );

        totalCntVol += cntVol;
        totalStgVol += stgVol;
    }

    // Verify that all the staggered volumes summed to roughly the same thing
    EXPECT_NEAR( totalStgVol.x, totalCntVol, max( fabs( totalStgVol.x ), fabs( totalCntVol ) ) * 0.0001f );
    EXPECT_NEAR( totalStgVol.x, totalStgVol.y, max( fabs( totalStgVol.x ), fabs( totalStgVol.y ) ) * 0.0001f );
    EXPECT_NEAR( totalStgVol.z, totalStgVol.y, max( fabs( totalStgVol.z ), fabs( totalStgVol.y ) ) * 0.0001f );

    // Verify that the staggered volume is roughly equal volume of the box that was input
    float voxelVolume = vcs.voxel_length() * vcs.voxel_length() * vcs.voxel_length();
    EXPECT_NEAR( voxelVolume * totalStgVol.x, baseBox.volume(),
                 max( fabs( voxelVolume * totalStgVol.x ), baseBox.volume() ) * 0.06f );

    float calculatedVolume = rls.compute_volume();
    EXPECT_NEAR( voxelVolume * totalStgVol.x, calculatedVolume,
                 max( fabs( voxelVolume * totalStgVol.x ), calculatedVolume ) * 0.06f );

    FF_LOG( debug ) << "Total centered volume (in cubic voxels for X, Y, Z): " << totalCntVol << endl;
    FF_LOG( debug ) << "Total staggered volume: " << totalStgVol << endl;
    FF_LOG( debug ) << "Original box volume: " << baseBox.volume() / voxelVolume << endl;
    FF_LOG( debug ) << "compute_volume volume: " << calculatedVolume / voxelVolume << endl;
}
