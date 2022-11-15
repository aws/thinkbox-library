// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/assign/std/vector.hpp>
#include <frantic/volumetrics/levelset/rle_defined_and_adj_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/levelset/rle_pairwise_run_iterator.hpp>

using namespace std;
using namespace boost;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::fluids;

class RleIndexSpec : public ::testing::Test {
  protected:
    virtual void SetUp() {
        coordinateArray.push_back( vector3( -1, 3, 5 ) );
        coordinateArray.push_back( vector3( 3, 3, 5 ) );
        coordinateArray.push_back( vector3( 2, 3, 8 ) );
        coordinateArray.push_back( vector3( -1, 2, 8 ) );
        coordinateArray.push_back( vector3( 3, 2, 8 ) );
        coordinateArray.push_back( vector3( 1, 4, 6 ) );
        coordinateArray.push_back( vector3( 2, 3, 5 ) );
        coordinateArray.push_back( vector3( 2, 3, 6 ) );
        coordinateArray.push_back( vector3( 2, 3, 7 ) );
        ris.build_from_voxel_array( coordinateArray );
    }
    rle_index_spec ris;
    std::vector<vector3> coordinateArray;
};

TEST_F( RleIndexSpec, Creation ) {

    rle_index_spec ris;

    // Test creation based on an array of voxel coordinates
    std::vector<vector3> coordinateArray;

    // Check the consistency, and check that the XYZtoDataIndex function finds all the vectors
    ASSERT_TRUE( ris.check_consistency( cout ) );
    ASSERT_EQ( coordinateArray.size(), ris.data_size() );
    for( unsigned i = 0; i < coordinateArray.size(); ++i ) {
        int dataIndex = ris.XYZtoDataIndex( coordinateArray[i] );
        EXPECT_EQ( dataIndex, i );
    }
}

TEST_F( RleIndexSpec, Get2x2x2DataIndexBlockFunction ) {

    rle_index_spec ris;

    boundbox3 bounds = ris.outer_bounds();
    boost::int32_t dataIndicesF2x2x2DIB[8];
    boost::int32_t dataIndicesFDIB[8];
    for( int z = bounds.minimum().z - 2; z < bounds.maximum().z + 2; ++z ) {
        for( int y = bounds.minimum().y - 2; y < bounds.maximum().y + 2; ++y ) {
            for( int x = bounds.minimum().x - 2; x < bounds.maximum().x + 2; ++x ) {
                // Get the block of data indices
                ris.fill_2x2x2_data_index_box( vector3( x, y, z ), dataIndicesF2x2x2DIB );
                // Get it again with the general function
                ris.fill_data_index_box( boundbox3( x, x + 1, y, y + 1, z, z + 1 ), dataIndicesFDIB );
                // Compare the values retrieved with the block function versus using the XYZtoDataIndex function
                for( int dz = 0; dz < 2; ++dz ) {
                    for( int dy = 0; dy < 2; ++dy ) {
                        for( int dx = 0; dx < 2; ++dx ) {
                            int rawDataIndex = ris.XYZtoDataIndex( vector3( x + dx, y + dy, z + dz ) );
                            int blockDataIndexF2x2x2DIB = dataIndicesF2x2x2DIB[dx + 2 * dy + 4 * dz];
                            int blockDataIndexFDIB = dataIndicesFDIB[dx + 2 * dy + 4 * dz];
                            EXPECT_EQ( rawDataIndex, blockDataIndexF2x2x2DIB )
                                << "outer bounds: " << bounds << "\n"
                                << "outer position: (" << x << "," << y << "," << z << ")\n"
                                << "inner position: (" << dx << "," << dy << "," << dz << "), index "
                                << dx + 2 * dy + 4 * dz << "\n"
                                << "raw data index: " << rawDataIndex << "\n"
                                << "2x2x2 block data index: " << blockDataIndexF2x2x2DIB << "\n"
                                << "general block data index: " << blockDataIndexFDIB << "\n";

                            EXPECT_EQ( rawDataIndex, blockDataIndexFDIB )
                                << "outer bounds: " << bounds << "\n"
                                << "outer position: (" << x << "," << y << "," << z << ")\n"
                                << "inner position: (" << dx << "," << dy << "," << dz << "), index "
                                << dx + 2 * dy + 4 * dz << "\n"
                                << "raw data index: " << rawDataIndex << "\n"
                                << "2x2x2 block data index: " << blockDataIndexF2x2x2DIB << "\n"
                                << "general block data index: " << blockDataIndexFDIB << "\n";
                        }
                    }
                }
            }
        }
    }
}

TEST_F( RleIndexSpec, GetDataIndexBlock ) {

    rle_index_spec ris;

    boundbox3 randomTestBounds = ris.outer_bounds();
    randomTestBounds.expand( 4 );
    for( int randomTest = 0; randomTest < 1000; ++randomTest ) {
        // Create a random box within the randomTestBounds
        boundbox3 testBox;
        /* The following lines were giving me an error that I cant get a radnom vector from an empty box. Commenting
           them out made the test pass. It might be because the LHS operand is not of the same type as the RHS operand
         */
        // testBox += randomTestBounds.random_vector();
        // testBox += randomTestBounds.random_vector();

        vector<boost::int32_t> dataIndices( testBox.get_volume() );

        // Get the indexes
        ris.fill_data_index_box( testBox, dataIndices.size() ? &dataIndices[0] : NULL );

        for( int z = testBox.minimum().z; z < testBox.maximum().z; ++z ) {
            for( int y = testBox.minimum().y; y < testBox.maximum().y; ++y ) {
                for( int x = testBox.minimum().x; x < testBox.maximum().x; ++x ) {
                    // Compare the values retrieved with the block function versus using the XYZtoDataIndex function
                    int rawDataIndex = ris.XYZtoDataIndex( vector3( x, y, z ) );
                    int blockDataIndex =
                        dataIndices[testBox.size().get_index( vector3( x, y, z ) - testBox.minimum() )];
                    EXPECT_EQ( rawDataIndex, blockDataIndex )
                        << "outer bounds: " << randomTestBounds << "\n"
                        << "data index box: " << testBox << "\n"
                        << "outer position: (" << x << "," << y << "," << z << ")\n"
                        << "raw data index: " << rawDataIndex << "\n"
                        << "general block data index: " << blockDataIndex << "\n";
                }
            }
        }
    }
}

TEST_F( RleIndexSpec, AxisPermutationfunction ) {

    rle_index_spec permuteRisA, permuteRisB;

    // Test that the identity (0,1,2) works
    permuteRisA.build_from_axis_permutation( vector3( 0, 1, 2 ), ris );
    ASSERT_TRUE( permuteRisA == ris );

    // Test that the transposition (1,0,2) works
    permuteRisA.build_from_axis_permutation( vector3( 1, 0, 2 ), ris );
    ASSERT_TRUE( permuteRisA != ris );
    permuteRisB.build_from_axis_permutation( vector3( 1, 0, 2 ), permuteRisA );
    ASSERT_TRUE( permuteRisB == ris );

    // Test that the transposition (0,2,1) works
    permuteRisA.build_from_axis_permutation( vector3( 0, 2, 1 ), ris );
    ASSERT_TRUE( permuteRisA != ris );
    permuteRisB.build_from_axis_permutation( vector3( 0, 2, 1 ), permuteRisA );
    ASSERT_TRUE( permuteRisB == ris );

    // Test that the transposition (2,1,0) works
    permuteRisA.build_from_axis_permutation( vector3( 2, 1, 0 ), ris );
    ASSERT_TRUE( permuteRisA != ris );
    permuteRisB.build_from_axis_permutation( vector3( 2, 1, 0 ), permuteRisA );
    ASSERT_TRUE( permuteRisB == ris );

    // Test that the 3-cycle (1,2,0) works
    permuteRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), ris );
    ASSERT_TRUE( permuteRisA != ris );
    permuteRisB.build_from_axis_permutation( vector3( 1, 2, 0 ), permuteRisA );
    ASSERT_TRUE( permuteRisB != ris );
    ASSERT_TRUE( permuteRisB != permuteRisA );
    permuteRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), permuteRisB );
    ASSERT_TRUE( permuteRisA == ris );

    // Test that the 3-cycle (2,0,1) works
    permuteRisA.build_from_axis_permutation( vector3( 2, 0, 1 ), ris );
    ASSERT_TRUE( permuteRisA != ris );
    permuteRisB.build_from_axis_permutation( vector3( 2, 0, 1 ), permuteRisA );
    ASSERT_TRUE( permuteRisB != ris );
    ASSERT_TRUE( permuteRisB != permuteRisA );
    permuteRisA.build_from_axis_permutation( vector3( 2, 0, 1 ), permuteRisB );
    ASSERT_TRUE( permuteRisA == ris );
}

TEST_F( RleIndexSpec, AdjacencyStructure ) {
    // Make sure the is_neighbor_index_direction_positive function works correctly
    EXPECT_TRUE( !is_neighbor_index_direction_positive( rae_index_x_neg ) );
    EXPECT_TRUE( is_neighbor_index_direction_positive( rae_index_x_pos ) );
    EXPECT_TRUE( !is_neighbor_index_direction_positive( rae_index_y_neg ) );
    EXPECT_TRUE( is_neighbor_index_direction_positive( rae_index_y_pos ) );
    EXPECT_TRUE( !is_neighbor_index_direction_positive( rae_index_z_neg ) );
    EXPECT_TRUE( is_neighbor_index_direction_positive( rae_index_z_pos ) );

    // Make sure the neighbor_index_axis function works correctly
    EXPECT_EQ( neighbor_index_axis( rae_index_x_neg ), 0 );
    EXPECT_EQ( neighbor_index_axis( rae_index_x_pos ), 0 );
    EXPECT_EQ( neighbor_index_axis( rae_index_y_neg ), 1 );
    EXPECT_EQ( neighbor_index_axis( rae_index_y_pos ), 1 );
    EXPECT_EQ( neighbor_index_axis( rae_index_z_neg ), 2 );
    EXPECT_EQ( neighbor_index_axis( rae_index_z_pos ), 2 );

    ris_adjacency adj;
    adj.compute( ris );
    for( rle_defined_iterator i = ris.begin(), ie = ris.end(); i != ie; ++i ) {
        // cout << "Testing data index " << i.get_data_index() << endl;
        EXPECT_EQ( adj[i.get_data_index()].x_neg, ris.XYZtoDataIndex( i.get_coord() - vector3( 1, 0, 0 ) ) );
        EXPECT_EQ( adj[i.get_data_index()].x_pos, ris.XYZtoDataIndex( i.get_coord() + vector3( 1, 0, 0 ) ) );
        EXPECT_EQ( adj[i.get_data_index()].y_neg, ris.XYZtoDataIndex( i.get_coord() - vector3( 0, 1, 0 ) ) );
        EXPECT_EQ( adj[i.get_data_index()].y_pos, ris.XYZtoDataIndex( i.get_coord() + vector3( 0, 1, 0 ) ) );
        EXPECT_EQ( adj[i.get_data_index()].z_neg, ris.XYZtoDataIndex( i.get_coord() - vector3( 0, 0, 1 ) ) );
        EXPECT_EQ( adj[i.get_data_index()].z_pos, ris.XYZtoDataIndex( i.get_coord() + vector3( 0, 0, 1 ) ) );
    }
}

TEST( RleIndexSpecTest, FillDataIndexBox ) {

    unsigned seed = (unsigned)time( NULL );
    srand( seed );
    FF_LOG( debug ) << "seed:" << seed << std::endl;

    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 20, boxTestCount = 20;
    for( int test = 0; test < testCount; ++test ) {
        rle_index_spec ris;
        ris.build_from_random( testBounds, 4 );

        // For each rle_index_spec, test it with a number of random boxes
        for( int boxTest = 0; boxTest < boxTestCount; ++boxTest ) {
            // Make a randomized outer_bounds
            boundbox3 dataIndexBounds( testBounds.random_vector() );
            dataIndexBounds += testBounds.random_vector();

            vector<boost::int32_t> indices( dataIndexBounds.get_volume() ),
                indicesBoundaryDuplicated( dataIndexBounds.get_volume() );
            ris.fill_data_index_box( dataIndexBounds, &indices[0] );
            ris.fill_data_index_box_boundary_duplicated( dataIndexBounds, &indicesBoundaryDuplicated[0] );

            vector3 coord;
            int index = 0;
            for( coord.z = dataIndexBounds.zminimum(); coord.z <= dataIndexBounds.zmaximum(); ++coord.z ) {
                for( coord.y = dataIndexBounds.yminimum(); coord.y <= dataIndexBounds.ymaximum(); ++coord.y ) {
                    for( coord.x = dataIndexBounds.xminimum(); coord.x <= dataIndexBounds.xmaximum();
                         ++coord.x, ++index ) {
                        vector3 clampedCoord = ris.outer_bounds().clamp( coord );
                        boost::int32_t dataIndex = ris.XYZtoDataIndex( coord ),
                                       dataIndexBoundaryDuplicated = ris.XYZtoDataIndex( clampedCoord );
                        EXPECT_EQ( dataIndex, indices[index] );
                        if( dataIndexBoundaryDuplicated != indicesBoundaryDuplicated[index] )
                            FF_LOG( debug ) << "coord: " << coord << ", outer bound: " << ris.outer_bounds()
                                            << ", exterior code: " << ris.get_exterior_region_code() << endl;
                        EXPECT_EQ( dataIndexBoundaryDuplicated, indicesBoundaryDuplicated[index] );
                    }
                }
            }
        }
    }
}

TEST( RleIndexSpecTest, DefinedIterator ) {
    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 1000;
    for( int test = 0; test < testCount; ++test ) {
        rle_index_spec ris;
        ris.build_from_random( testBounds, 4 );

        // Validate the defined iterator
        int index = 0;
        for( rle_defined_iterator i = ris.begin(), ie = ris.end(); i != ie; ++i, ++index ) {
            // cout << "coord: " << i.get_coord() << ", data index: " << i.get_data_index() << endl;
            EXPECT_EQ( i.get_data_index(), index );
            EXPECT_EQ( i.get_data_index(), ris.XYZtoDataIndex( i.get_coord() ) );
        }
        // Validate that we visited every defined voxel
        EXPECT_EQ( ris.data_size(), index );
    }
}

TEST( RleIndexSpecTest, DefinedAndAdjIterator ) {
    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 1000;
    for( int test = 0; test < testCount; ++test ) {
        rle_index_spec ris;
        ris.build_from_random( testBounds, 4 );

        // Validate the defined and adjacent iterator
        int index = 0;
        for( rle_defined_and_adj_iterator i( ris ), ie( ris, true ); i != ie; ++i, ++index ) {
            EXPECT_EQ( i.get_center_data_index(), index );
            EXPECT_EQ( i.get_center_data_index(), ris.XYZtoDataIndex( i.get_coord() ) );
            EXPECT_EQ( i.get_adjacent_data_index( rae_index_x_neg ),
                       ris.XYZtoDataIndex( i.get_coord() + vector3( -1, 0, 0 ) ) );
            EXPECT_EQ( i.get_adjacent_data_index( rae_index_x_pos ),
                       ris.XYZtoDataIndex( i.get_coord() + vector3( 1, 0, 0 ) ) );
            EXPECT_EQ( i.get_adjacent_data_index( rae_index_y_neg ),
                       ris.XYZtoDataIndex( i.get_coord() + vector3( 0, -1, 0 ) ) );
            EXPECT_EQ( i.get_adjacent_data_index( rae_index_y_pos ),
                       ris.XYZtoDataIndex( i.get_coord() + vector3( 0, 1, 0 ) ) );
            EXPECT_EQ( i.get_adjacent_data_index( rae_index_z_neg ),
                       ris.XYZtoDataIndex( i.get_coord() + vector3( 0, 0, -1 ) ) );
            EXPECT_EQ( i.get_adjacent_data_index( rae_index_z_pos ),
                       ris.XYZtoDataIndex( i.get_coord() + vector3( 0, 0, 1 ) ) );
        }
        // Validate that we visited every defined voxel
        EXPECT_EQ( ris.data_size(), index );
    }
}

TEST( RleIndexSpecTest, RunIterator ) {

    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 1000;
    for( int test = 0; test < testCount; ++test ) {
        rle_index_spec ris;
        ris.build_from_random( testBounds, 4 );

        int definedVoxelCount = 0;
        // Validate the run iterator
        boundbox3 outerBounds = ris.outer_bounds();
        for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
            for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
                rle_run_iterator i( ris, ris.y_to_b( y ), ris.z_to_c( z ) ), ie;
                int previousRunEndX = i.get_xmin();
                // Iterate through all the runs in this scanline
                for( ; i != ie; ++i ) {
                    if( i.get_data_index() >= 0 )
                        definedVoxelCount += i.get_xsize();
                    pair<int, int> runIndexRange =
                        ris.x_interval_to_run_index_interval( i.get_xmin(), i.get_xmax(), y, z );
                    EXPECT_EQ( runIndexRange.first, runIndexRange.second );
                    EXPECT_EQ( i.get_data_index(), ris.run_index_to_data_index( runIndexRange.first ) );
                    // one past the last run should always equal the start of this run
                    EXPECT_EQ( previousRunEndX, i.get_xmin() );
                    previousRunEndX = i.get_xmin() + i.get_xsize();
                }
                // As a special case, once the run iterator is done, it should contain the X value one past the scanline
                // end
                EXPECT_EQ( previousRunEndX, i.get_xmin() );
            }
        }
        // Validate that we visited every defined voxel
        EXPECT_EQ( ris.data_size(), definedVoxelCount );
    }
}

TEST( RleIndexSpecTest, FillDataIndexMap ) {

    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 1000;
    for( int test = 0; test < testCount; ++test ) {
        rle_index_spec ris1, ris2;
        ris1.build_from_random( testBounds, 4 );
        ris2.build_from_random( testBounds, 4 );

        // Create the data index map
        vector<boost::int32_t> diMap( ris1.data_size() );
        if( ris1.data_size() > 0 )
            ris1.fill_data_index_map( ris2, &diMap[0] );

        // Validate the data index map
        for( rle_defined_iterator i = ris1.begin(), ie = ris1.end(); i != ie; ++i ) {
            EXPECT_EQ( diMap[i.get_data_index()], ris2.XYZtoDataIndex( i.get_coord() ) );
        }
    }
}

TEST( RleIndexSpecTest, BuildWithTrimBounds ) {

    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 1000;
    for( int test = 0; test < testCount; ++test ) {
        // cout << "Test " << test << " of " << testCount << endl;
        rle_index_spec ris;
        ris.build_from_random( testBounds, 4 );

        // Create a random trimBounds bounding box
        boundbox3 trimBounds( testBounds.random_vector() );
        trimBounds += testBounds.random_vector();

        rle_index_spec risTrim;
        risTrim.build_with_trim_bounds( ris, trimBounds );

        ASSERT_TRUE( risTrim.check_consistency( std::cout ) );

        // Ensure that the outer bounds are inside the trim bounds
        ASSERT_TRUE( risTrim.outer_bounds().is_empty() || trimBounds.contains( risTrim.outer_bounds() ) );

        ASSERT_EQ( ris.get_exterior_region_code(), risTrim.get_exterior_region_code() );

        // Ensure that all the voxels inside the trim bounds matched up correctly
        vector3 coord;
        for( coord.z = trimBounds.zminimum(); coord.z <= trimBounds.zmaximum(); ++coord.z ) {
            for( coord.y = trimBounds.yminimum(); coord.y <= trimBounds.ymaximum(); ++coord.y ) {
                for( coord.x = trimBounds.xminimum(); coord.x <= trimBounds.xmaximum(); ++coord.x ) {
                    int risDataIndex = ris.XYZtoDataIndex( coord );
                    int risTrimDataIndex = risTrim.XYZtoDataIndex( coord );
                    // Either both should be positive, or they should both be negative and equal
                    if( risDataIndex < 0 || risTrimDataIndex < 0 )
                        EXPECT_EQ( risDataIndex, risTrimDataIndex );
                }
            }
        }
    }
}

TEST( RleIndexSpecTest, PairwiseRunIterator ) {

    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 1000;
    for( int test = 0; test < testCount; ++test ) {
        rle_index_spec risA, risB;
        risA.build_from_random( testBounds, 8 );
        risB.build_from_random( testBounds, 8 );

        // These should be all set to 1 during the iteration
        vector<char> touchedA( risA.data_size(), 0 ), touchedB( risB.data_size(), 0 );

        int subIntervalCount = 0;

        for( int z = testBounds.zminimum(); z <= testBounds.zmaximum(); ++z ) {
            for( int y = testBounds.yminimum(); y <= testBounds.ymaximum(); ++y ) {
                rle_pairwise_run_iterator i( risA, risA.y_to_b( y ), risA.z_to_c( z ), risB, risB.y_to_b( y ),
                                             risB.z_to_c( z ) ),
                    ie;
                int x = i.get_xmin();
                for( ; i != ie; ++i ) {
                    ++subIntervalCount;
                    // Ensure that all the intervals we process are contiguous
                    EXPECT_EQ( x, i.get_xmin() );
                    x += i.get_xsize();
                    // Check that the run for risA matches all the data indices correctly
                    if( i.get_first_data_index() >= 0 ) {
                        for( int d = 0; d < i.get_xsize(); ++d ) {
                            vector3 vox( i.get_xmin() + d, y, z );
                            if( i.get_first_data_index() + d != risA.XYZtoDataIndex( vox ) )
                                cout << "vox: " << vox << "\n";
                            EXPECT_EQ( i.get_first_data_index() + d, risA.XYZtoDataIndex( vox ) );
                            touchedA[i.get_first_data_index() + d] = 1;
                        }
                    } else {
                        for( int d = 0; d < i.get_xsize(); ++d ) {
                            EXPECT_EQ( i.get_first_data_index(),
                                       risA.XYZtoDataIndex( vector3( i.get_xmin() + d, y, z ) ) );
                        }
                    }
                    // Check that the run for risB matches all the data indices correctly
                    if( i.get_second_data_index() >= 0 ) {
                        for( int d = 0; d < i.get_xsize(); ++d ) {
                            EXPECT_EQ( i.get_second_data_index() + d,
                                       risB.XYZtoDataIndex( vector3( i.get_xmin() + d, y, z ) ) );
                            touchedB[i.get_second_data_index() + d] = 1;
                        }
                    } else {
                        for( int d = 0; d < i.get_xsize(); ++d ) {
                            EXPECT_EQ( i.get_second_data_index(),
                                       risB.XYZtoDataIndex( vector3( i.get_xmin() + d, y, z ) ) );
                        }
                    }
                }
                // After iteration, get_xmin() should equal the ending X.
                EXPECT_EQ( x, i.get_xmin() );
            }
        }

        // Verify that this iteration caught all of the defined voxels in both A and B
        for( size_t i = 0; i < touchedA.size(); ++i )
            EXPECT_EQ( touchedA[i], 1 );
        for( size_t i = 0; i < touchedB.size(); ++i )
            EXPECT_EQ( touchedB[i], 1 );
    }
}

TEST( RleIndexSpecTest, CombineForBlend ) {

    // testBounds is the bounds within which this whole test occurs
    boundbox3 testBounds( -20, 20, -20, 20, -20, 20 );
    int testCount = 100;
    for( int test = 0; test < testCount; ++test ) {
        rle_index_spec risA, risB;
        risA.build_from_random( testBounds, 8, 2 );
        risB.build_from_random( testBounds, 8, 2 );

        // The combine_for_blend function requires that the exterior region codes match, so let's force that
        risB.set_exterior_region_code( risA.get_exterior_region_code() );

        // Do the combine operation
        rle_index_spec ris;
        ris.combine_for_blend( risA, risB );

        // Ensure that all the voxels inside the test bounds are combined appropriately in the combined ris
        vector3 coord;
        for( coord.z = testBounds.zminimum(); coord.z <= testBounds.zmaximum(); ++coord.z ) {
            for( coord.y = testBounds.yminimum(); coord.y <= testBounds.ymaximum(); ++coord.y ) {
                for( coord.x = testBounds.xminimum(); coord.x <= testBounds.xmaximum(); ++coord.x ) {
                    int risDataIndex = ris.XYZtoDataIndex( coord );
                    int risDataIndexA = risA.XYZtoDataIndex( coord );
                    int risDataIndexB = risB.XYZtoDataIndex( coord );
                    if( risDataIndexA >= 0 || risDataIndexB >= 0 || risDataIndexA != risDataIndexB ) {
                        // the combined result should be defined
                        EXPECT_LE( 0, risDataIndex );
                    } else {
                        // the combined result should be undefined and match risDataIndexA
                        EXPECT_EQ( risDataIndex, risDataIndexA );
                        // the combined result should be undefined and match risDataIndexB
                        EXPECT_EQ( risDataIndex, risDataIndexB );
                    }
                }
            }
        }
    }
}
