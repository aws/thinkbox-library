// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/assign/std/vector.hpp>

#include <frantic/files/files.hpp>
#include <frantic/strings/tstring.hpp>
#include <frantic/volumetrics/levelset/level_set_fast_marching_reinitialization.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

#define RLS_PRINT_2D( rls )                                                                                            \
    do {                                                                                                               \
        stringstream ss;                                                                                               \
        rls.print2d( ss );                                                                                             \
        FF_LOG( debug ) << frantic::strings::to_tstring( ss.str() );                                                   \
    } while( 0 )

using namespace std;
using namespace boost::assign;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::volumetrics::levelset::detail;
using namespace frantic::files;
using namespace frantic::channels;
using namespace frantic::graphics;

class ReinitializeSignedDistance : public ::testing::Test {
  protected:
    void SetUp() {

        vcs = voxel_coord_system( vector3f( 0 ), 1 );
        rls = rle_level_set( vcs );
        rls_reinit = rle_level_set( vcs );

        frantic::volumetrics::levelset::detail::build_test1( rls );

        frantic::logging::set_logging_level( 5 );
        RLS_PRINT_2D( rls );

        // Count is sizeof(correct_distanceData)/sizeof(correct_distanceData[0])
        originalData += 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f,
            -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f, 1.0f;
        // Count is sizeof(correct_distanceData)/sizeof(correct_distanceData[0])
        bandNotUpdatedData += 1.70711f, 1.0f, 1.0f, 1.0f, 1.0f, 1.70711f, 1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 1.0f, 1.0f,
            -1.0f, -1.70711f, -1.70711f, -1.0f, 1.0f, 1.0f, -1.0f, -1.70711f, -1.70711f, -1.0f, 1.0f, 1.0f, -1.0f,
            -1.0f, -1.0f, -1.0f, 1.0f, 1.70711f, 1.0f, 1.0f, 1.0f, 1.0f, 1.70711f;
        // Count is sizeof(correct_distanceData)/sizeof(correct_distanceData[0])
        bandUpdatedData += 1.20711f, 0.5f, 0.5f, 0.5f, 0.5f, 1.20711f, 0.5f, -0.353553f, -0.5f, -0.5f, -0.353553f, 0.5f,
            0.5f, -0.5f, -1.20711f, -1.20711f, -0.5f, 0.5f, 0.5f, -0.5f, -1.20711f, -1.20711f, -0.5f, 0.5f, 0.5f,
            -0.353553f, -0.5f, -0.5f, -0.353553f, 0.5f, 1.20711f, 0.5f, 0.5f, 0.5f, 0.5f, 1.20711f;
        // Count is sizeof(correct_distanceData)/sizeof(correct_distanceData[0])
        zeroStoppingWidthBandUpdatedData += 1.0f, 0.5f, 0.5f, 0.5f, 0.5f, 1.0f, 0.5f, -0.353553f, -0.5f, -0.5f,
            -0.353553f, 0.5f, 0.5f, -0.5f, -1.0f, -1.0f, -0.5f, 0.5f, 0.5f, -0.5f, -1.0f, -1.0f, -0.5f, 0.5f, 0.5f,
            -0.353553f, -0.5f, -0.5f, -0.353553f, 0.5f, 1.0f, 0.5f, 0.5f, 0.5f, 0.5f, 1.0f;
    }

    rle_level_set rls, rls_reinit;
    voxel_coord_system vcs;
    vector<float> originalData;
    vector<float> bandNotUpdatedData;
    vector<float> bandUpdatedData;
    vector<float> zeroStoppingWidthBandUpdatedData;
};

TEST_F( ReinitializeSignedDistance, ZeroStoppingDistance ) {

    // Note: this should only update the initial band distances

    rls_reinit = rls;
    rls_reinit.reinitialize_signed_distance( _T(""), 0, 0 );

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], zeroStoppingWidthBandUpdatedData[i], 1e-5 );
    }

    rls_reinit = rls;
    rls_reinit.reinitialize_signed_distance( _T("Populated"), 0, 0 );

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], zeroStoppingWidthBandUpdatedData[i], 1e-5 );
    }
}

TEST_F( ReinitializeSignedDistance, InfinityStoppingDistance ) {

    // Note: this should update all of the distances

    rls_reinit = rls;
    rls_reinit.reinitialize_signed_distance( _T(""), -std::numeric_limits<float>::infinity(),
                                             std::numeric_limits<float>::infinity() );

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], bandUpdatedData[i], 1e-5 );
    }

    rls_reinit = rls;
    rls_reinit.tag_interface_voxels( _T("Populated") );
    rls_reinit.reinitialize_signed_distance( _T("Populated"), -std::numeric_limits<float>::infinity(),
                                             std::numeric_limits<float>::infinity() );

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], bandUpdatedData[i], 1e-5 );
    }
}

TEST_F( ReinitializeSignedDistance, DefaultInputs ) {

    rls_reinit = rls;
    rls_reinit.tag_interface_voxels( _T("Populated") );
    rls_reinit.reinitialize_signed_distance();

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], bandUpdatedData[i], 1e-5 );
    }

    rls_reinit = rls;
    rls_reinit.reinitialize_signed_distance();

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], bandUpdatedData[i], 1e-5 );
    }
}

TEST_F( ReinitializeSignedDistance, DistanceLesserThanOne ) {

    // Test reinitialize_signed_distance() to ensure that the signed distance value at each voxel is no more
    // than one voxel length different than each of its adjacent neighbours

    rls_reinit = rls;

    const ris_adjacency& adj = rls_reinit.get_rle_index_spec().get_cached_adjacency();

    rls_reinit.tag_interface_voxels( _T("Populated") );
    rls_reinit.reinitialize_signed_distance();

    float voxelLength = rls_reinit.get_voxel_coord_system().voxel_length();
    stringstream ss;
    rls_reinit.dump( ss );
    FF_LOG( debug ) << frantic::strings::to_tstring( ss.str() );
    // for each defined voxel
    for( int iterIndex = 0, iterIndexEnd = (int)adj.data_size(); iterIndex != iterIndexEnd; ++iterIndex ) {

        // find the neighbour indices
        int xposDataIndex = adj[iterIndex].x_pos;
        int xnegDataIndex = adj[iterIndex].x_neg;
        int yposDataIndex = adj[iterIndex].y_pos;
        int ynegDataIndex = adj[iterIndex].y_neg;
        int zposDataIndex = adj[iterIndex].z_pos;
        int znegDataIndex = adj[iterIndex].z_neg;

        // ensure that none of the distance values of the neighbours are more than one voxel length (plus a small
        // epsilon for rounding error) from the current distance value
        if( xposDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[xposDataIndex], voxelLength + 1e-6 );
        }
        if( xnegDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[xnegDataIndex], voxelLength + 1e-6 );
        }
        if( yposDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[yposDataIndex], voxelLength + 1e-6 );
        }
        if( ynegDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[ynegDataIndex], voxelLength + 1e-6 );
        }
        if( zposDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[zposDataIndex], voxelLength + 1e-6 );
        }
        if( znegDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[znegDataIndex], voxelLength + 1e-6 );
        }
    }
}

TEST_F( ReinitializeSignedDistance, FromPopulatedZeroStoppingDistance ) {

    // Test reinitialize_signed_distance_from_populated() with a stopping distance of zero
    // Note: this should only update the initial band distances

    rls_reinit = rls;
    rls_reinit.tag_interface_voxels( _T("Populated") );
    rls_reinit.reinitialize_signed_distance_from_populated( _T("Populated"), 0, 0 );

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], originalData[i], 1e-5 );
    }
}

TEST_F( ReinitializeSignedDistance, FromPopulatedInfinityStoppingDistance ) {

    // Test reinitialize_signed_distance() with a stopping distance of infinity
    // Note: this should update all of the non-band distances

    rls_reinit = rls;
    rls_reinit.tag_interface_voxels( _T("Populated") );
    rls_reinit.reinitialize_signed_distance_from_populated( _T("Populated"), -std::numeric_limits<float>::infinity(),
                                                            std::numeric_limits<float>::infinity() );

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], bandNotUpdatedData[i], 1e-5 );
    }
}

TEST_F( ReinitializeSignedDistance, FromPopulatedDefaultInputs ) {

    rls_reinit = rls;
    rls_reinit.tag_interface_voxels( _T("Populated") );
    rls_reinit.reinitialize_signed_distance_from_populated( _T("Populated") );

    for( unsigned int i = 0; i < rls.size(); ++i ) {
        // compare with the expected output (allow for small rounding errors in the hardcoded data values)
        EXPECT_NEAR( rls_reinit[i], bandNotUpdatedData[i], 1e-5 );
    }
}

TEST_F( ReinitializeSignedDistance, FromPopulatedDistanceLesserThanOne ) {
    ////////////////////////////////////////////////////////
    // Test reinitialize_signed_distance_from_populated() to ensure that the signed distance value at each voxel is no
    // more
    // than one voxel length different than each of its adjacent neighbours
    ////////////////////////////////////////////////////////

    frantic::volumetrics::levelset::detail::build_test1( rls_reinit );

    // The level set this test operates on has been altered to a more random level set.
    // A coin flip is used to determine whether the level set will be biased to inside/outside.
    // Then random distance values are generated towards the bias.  Then a couple
    // of values of the opposite sign are generated.  It's done this way to ensure
    // that we don't general mostly (or only) interface voxels, so that some significant
    // marching takes place.
    unsigned seed = (unsigned)time( NULL );
    srand( seed );
    FF_LOG( debug ) << "seed:" << seed << std::endl;

    float bias;
    if( ( float( rand() ) / RAND_MAX ) - 0.5f < 0.f )
        bias = -1.f;
    else
        bias = 1.f;

    for( int i = 0; i < 6; ++i )
        for( int j = 0; j < 6; ++j )
            rls_reinit[i * 6 + j] = bias * ( float( rand() ) / RAND_MAX ) * 0.5f;

    for( int i = 0; i < 4; ++i ) {
        int r = (int)( ( float( rand() ) / RAND_MAX ) * rls_reinit.size() );
        rls_reinit[r] = -bias * ( float( rand() ) / RAND_MAX ) * 0.5f;
    }

    const ris_adjacency& adj2 = rls_reinit.get_rle_index_spec().get_cached_adjacency();

    rls_reinit.tag_interface_voxels( _T("Populated") );
    rls_reinit.reinitialize_signed_distance_from_populated( _T("Populated") );

    float voxelLength = rls_reinit.get_voxel_coord_system().voxel_length();

    // for each defined voxel
    for( int iterIndex = 0, iterIndexEnd = (int)adj2.data_size(); iterIndex != iterIndexEnd; ++iterIndex ) {

        // find the neighbour indices
        int xposDataIndex = adj2[iterIndex].x_pos;
        int xnegDataIndex = adj2[iterIndex].x_neg;
        int yposDataIndex = adj2[iterIndex].y_pos;
        int ynegDataIndex = adj2[iterIndex].y_neg;
        int zposDataIndex = adj2[iterIndex].z_pos;
        int znegDataIndex = adj2[iterIndex].z_neg;

        // ensure that none of the distance values of the neighbours are more than one voxel length (plus a small
        // epsilon for rounding error) from the current distance value
        if( xposDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[xposDataIndex], voxelLength + 1e-6 );
        }
        if( xnegDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[xnegDataIndex], voxelLength + 1e-6 );
        }
        if( yposDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[yposDataIndex], voxelLength + 1e-6 );
        }
        if( ynegDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[ynegDataIndex], voxelLength + 1e-6 );
        }
        if( zposDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[zposDataIndex], voxelLength + 1e-6 );
        }
        if( znegDataIndex > 0 ) {
            EXPECT_NEAR( rls_reinit[iterIndex], rls_reinit[znegDataIndex], voxelLength + 1e-6 );
        }
    }
}

TEST( ReinitializationRle, IgnoreFlag ) {
    // input populated :
    //  2  0  1  1  1  0  2
    // input signed distance :
    // 64  0 -1  0  1  0 32

    // expected output signed distance :
    // 64 -2 -1  0  1  2 32

    // Note: the populated flag's post condition is ill-defined right now --
    // the documentation is not consistent with the actual behaviour.
    {
        const int voxelCount = 7;

        frantic::volumetrics::voxel_coord_system vcs( frantic::graphics::vector3f( 0 ), 1.f );

        std::vector<frantic::graphics::vector3> voxelArray;
        voxelArray.reserve( voxelCount );
        for( boost::int32_t i = 0; i < voxelCount; ++i ) {
            voxelArray.push_back( frantic::graphics::vector3( i, 0, 0 ) );
        }

        frantic::volumetrics::levelset::rle_index_spec ris;
        ris.build_from_voxel_array( voxelArray );

        const float initialPhi[voxelCount] = { 64.f, 0, -1.f, 0, 1.f, 0, 32.f };
        std::vector<float> phi( initialPhi, initialPhi + voxelCount );
        EXPECT_EQ( phi.size(), voxelCount );

        frantic::volumetrics::levelset::rle_level_set ls( vcs, ris, phi, 5.f, 5.f );

        ls.add_channel<boost::uint8_t>( _T("Populated") );
        ls.zero_channel( _T("Populated") );
        const boost::uint8_t initialPop[voxelCount] = { 2, 0, 1, 1, 1, 0, 2 };
        frantic::volumetrics::levelset::rle_channel_accessor<boost::uint8_t> pop =
            ls.get_channel_accessor<boost::uint8_t>( _T("Populated") );
        memcpy( &pop[0], &initialPop[0], voxelCount * sizeof( boost::uint8_t ) );

        ls.reinitialize_signed_distance_from_populated( _T("Populated") );

        EXPECT_FLOAT_EQ( ls[0], 64.f );
        EXPECT_FLOAT_EQ( ls[1], -2.f );
        EXPECT_FLOAT_EQ( ls[2], -1.f );
        EXPECT_FLOAT_EQ( ls[3], 0 );
        EXPECT_FLOAT_EQ( ls[4], 1.f );
        EXPECT_FLOAT_EQ( ls[5], 2.f );
        EXPECT_FLOAT_EQ( ls[6], 32.f );
    }
}
TEST( ReinitializationRle, IgnoreFlagOnInterface ) {
    // Ignored voxels must not be used to change the initial band.
    //
    // input populated :
    //  2  2
    //  1  1
    // input signed distance :
    // 10 10
    // -1  1
    // expected output signed distance :
    //  10   10
    //  -0.5  0.5
    {
        const int voxelCount = 4;

        const float voxelLength = 1.f;

        frantic::volumetrics::voxel_coord_system vcs( frantic::graphics::vector3f( 0 ), voxelLength );

        std::vector<frantic::graphics::vector3> voxelArray;
        voxelArray.reserve( voxelCount );
        for( int j = 0; j < 2; ++j ) {
            for( int i = 0; i < 2; ++i ) {
                voxelArray.push_back( frantic::graphics::vector3( i, j, 0 ) );
            }
        }

        frantic::volumetrics::levelset::rle_index_spec ris;
        ris.build_from_voxel_array( voxelArray );

        const float initialPhi[voxelCount] = { 10.f, 10.f, -voxelLength, voxelLength };
        std::vector<float> phi( initialPhi, initialPhi + voxelCount );
        EXPECT_EQ( phi.size(), voxelCount );

        frantic::volumetrics::levelset::rle_level_set ls( vcs, ris, phi, 5.f, 5.f );

        ls.add_channel<boost::uint8_t>( _T("Populated") );
        ls.zero_channel( _T("Populated") );
        const boost::uint8_t initialPop[voxelCount] = { 2, 2, 1, 1 };
        frantic::volumetrics::levelset::rle_channel_accessor<boost::uint8_t> pop =
            ls.get_channel_accessor<boost::uint8_t>( _T("Populated") );
        memcpy( &pop[0], &initialPop[0], voxelCount * sizeof( boost::uint8_t ) );

        frantic::logging::null_progress_logger nullProgressLogger;

        fast_marching_reinitialization( nullProgressLogger, ls.get_rle_index_spec(), &pop[0], &ls[0], voxelLength, true,
                                        5.f * voxelLength, 5.f * voxelLength );

        EXPECT_FLOAT_EQ( ls[0], 10.f );
        EXPECT_FLOAT_EQ( ls[1], 10.f );
        EXPECT_FLOAT_EQ( ls[2], -0.5f * voxelLength );
        EXPECT_FLOAT_EQ( ls[3], 0.5f * voxelLength );
    }
}
