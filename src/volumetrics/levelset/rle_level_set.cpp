// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/level_set_fast_marching_reinitialization.hpp>
#include <frantic/volumetrics/levelset/level_set_marching_extrapolation.hpp>
#include <frantic/volumetrics/levelset/rle_defined_and_adj_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_defined_box_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/levelset/rle_pairwise_run_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_run_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

#include <frantic/volumetrics/rle_weno_interpolation.hpp>

#pragma warning( push )
#pragma warning( disable : 4512 4100 )
#include <tbb/task_scheduler_init.h>
#pragma warning( pop )

using namespace std;
using namespace frantic;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::geometry;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::fluids; // for rle_voxel_field

namespace frantic {
namespace volumetrics {
namespace levelset {

namespace detail {

/*
 * This function creates a specific 10x10 rle level set.
 * Note: This is used inside the unit tests, specifically in testReinitializationRLE, so any changes to
 * the output of this function should be updated there
 *
 * @param  levelSet  The level set where the 10x10 test level set will be stored.
 */
void build_test1( rle_level_set& levelSet ) {
    rle_index_spec& ris = levelSet.m_rleIndex;

    // set up the bounding box
    vector3 orig( 0, 0, 0 );
    size3 coordSize( 10, 10, 1 );
    ris.m_abcCoordOrigin = orig;
    ris.m_abcCoordSize = coordSize;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    // set the exterior region
    ris.m_exteriorRegionCode = -1;

    // write the run data
    int bcIndex;

    // Iterate over all the runs
    // for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c) {
    for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {

        bcIndex = b; //+ c * ris.m_abcCoordSize.ysize();
        ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();
        int a;
        if( b >= 2 && b <= ris.m_abcCoordSize.ysize() - 3 ) {
            for( a = 2; a <= ris.m_abcCoordSize.xsize() - 3; ++a ) {
                // insert run
                if( a == 2 ) {
                    ris.m_runData.push_back( run_data( a, (int)levelSet.m_distanceData.size() ) );
                    levelSet.m_distanceData.push_back( 1.f );
                } else if( a == ris.m_abcCoordSize.xsize() - 3 ) {
                    levelSet.m_distanceData.push_back( 1.f );
                } else if( b == 2 || b == ris.m_abcCoordSize.ysize() - 3 ) {
                    levelSet.m_distanceData.push_back( 1.f );
                } else {
                    levelSet.m_distanceData.push_back( -1.f );
                }
            }
            ris.m_runData.push_back( run_data( a, ris.m_exteriorRegionCode ) );
        } else {
            // make a zero-sized exterior run
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
        }
        //}
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = levelSet.m_distanceData.size();
}

// Build an RLE level set that we can use for automatic testing in the test suite
void build_testSuite( rle_level_set& levelSet, vector3 orig, const size3 coordSize, const int exteriorRegionCode ) {
    rle_index_spec& ris = levelSet.m_rleIndex;

    // set up the bounding box
    ris.m_abcCoordOrigin = orig;
    ris.m_abcCoordSize = coordSize;

    // set the exterior region
    ris.m_exteriorRegionCode = exteriorRegionCode;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    // Iterate over all the runs
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {

            int bcIndex = b + c * ris.m_abcCoordSize.ysize();

            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            // empty run
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = levelSet.m_distanceData.size();
}

/*
 * This function creates a specific 10x10 rle level set.
 * Note: This is used inside the unit tests, so any changes to the output of this function should be updated there
 *
 * @param  levelSet  The level set where the 10x10 test level set will be stored.
 */
void build_test1_complement( rle_level_set& levelSet ) {
    rle_index_spec& ris = levelSet.m_rleIndex;

    // set up the bounding box
    vector3 orig( 0, 0, 0 );
    size3 coordSize( 10, 10, 1 );
    ris.m_abcCoordOrigin = orig;
    ris.m_abcCoordSize = coordSize;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    // set the exterior region
    ris.m_exteriorRegionCode = -2;

    // write the run data
    int bcIndex;

    // Iterate over all the runs
    for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {

        bcIndex = b;
        ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();
        int a;
        if( b >= 2 && b <= ris.m_abcCoordSize.ysize() - 3 ) {
            for( a = 2; a <= ris.m_abcCoordSize.xsize() - 3; ++a ) {
                // insert run
                if( a == 2 ) {
                    ris.m_runData.push_back( run_data( a, (int)levelSet.m_distanceData.size() ) );
                    levelSet.m_distanceData.push_back( -1.0f );
                } else if( a == ris.m_abcCoordSize.xsize() - 3 ) {
                    levelSet.m_distanceData.push_back( -1.0f );
                } else if( b == 2 || b == ris.m_abcCoordSize.ysize() - 3 ) {
                    levelSet.m_distanceData.push_back( -1.0f );
                } else {
                    levelSet.m_distanceData.push_back( 1.0f );
                }
            }
            ris.m_runData.push_back( run_data( a, ris.m_exteriorRegionCode ) );
        } else {
            // make a zero-sized exterior run
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = levelSet.m_distanceData.size();
}

/*
 * This function creates a specific 10x10 rle level set.
 * Note: This is used inside the unit tests, so any changes to the output of this function should be updated there
 *
 * @param  levelSet  The level set where the 10x10 test level set will be stored.
 */
void build_test2( rle_level_set& levelSet ) {
    rle_index_spec& ris = levelSet.m_rleIndex;

    // set up the bounding box
    vector3 orig( 0, 0, 0 );
    size3 coordSize( 10, 10, 1 );
    ris.m_abcCoordOrigin = orig;
    ris.m_abcCoordSize = coordSize;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    // set the exterior region
    ris.m_exteriorRegionCode = -1;

    // write the run data
    int bcIndex;

    // Iterate over all the runs
    for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {

        bcIndex = b;
        ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();
        int a;
        if( b >= 4 && b <= ris.m_abcCoordSize.ysize() - 4 ) {
            for( a = 4; a <= ris.m_abcCoordSize.xsize() - 4; ++a ) {
                // insert run
                if( a == 4 ) {
                    ris.m_runData.push_back( run_data( a, (int)levelSet.m_distanceData.size() ) );
                    levelSet.m_distanceData.push_back( 1.0f );
                } else if( a == ris.m_abcCoordSize.xsize() - 4 ) {
                    levelSet.m_distanceData.push_back( 1.0f );
                } else if( b == 4 || b == ris.m_abcCoordSize.ysize() - 4 ) {
                    levelSet.m_distanceData.push_back( 1.0f );
                } else {
                    levelSet.m_distanceData.push_back( -1.0f );
                }
            }
            ris.m_runData.push_back( run_data( a, ris.m_exteriorRegionCode ) );
        } else {
            // make a zero-sized exterior run
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
            ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = levelSet.m_distanceData.size();
}

/*
 * This function creates a specific 10x10 rle level set.
 * Note: This is used inside the unit tests, so any changes to the output of this function should be updated there
 *
 * @param  levelSet  The level set where the 10x10 test level set will be stored.
 */
void build_test12_union( rle_level_set& levelSet ) { build_test1( levelSet ); }

/*
 * This function creates a specific 10x10 rle level set.
 * Note: This is used inside the unit tests, so any changes to the output of this function should be updated there
 *
 * @param  levelSet  The level set where the 10x10 test level set will be stored.
 */
void build_test12_intersect( rle_level_set& levelSet ) { build_test2( levelSet ); }
void build_plane_test( rle_level_set& levelSet ) {
    rle_index_spec& ris = levelSet.m_rleIndex;

    // set up the bounding box
    vector3 orig( 0, 0, 0 );
    size3 coordSize( 10, 10, 10 );
    ris.m_abcCoordOrigin = orig;
    ris.m_abcCoordSize = coordSize;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    // set the exterior region
    ris.m_exteriorRegionCode = -1;

    // write the run data
    int bcIndex;

    // Iterate over all the runs
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {

            bcIndex = b + c * ris.m_abcCoordSize.ysize();
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();
            // int a;

            // insert inside run
            ris.m_runData.push_back( run_data( 0, -2 ) );

            // insert outside run
            ris.m_runData.push_back( run_data( 5, -1 ) );

            // end the runs
            ris.m_runData.push_back( run_data( 10, -1 ) );
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = levelSet.m_distanceData.size();
}

} // namespace detail

// this prints out 2d rle level sets
void rle_level_set::print2d( std::ostream& out ) const {
    const rle_index_spec& ris = m_rleIndex;

    // check that zsize = 1 (we don't want to print 3d)
    if( !( ris.m_abcCoordSize.zsize() == 1 ) ) {
        out << "\tsize is " << m_rleIndex.m_abcCoordSize.xsize() << "x" << m_rleIndex.m_abcCoordSize.ysize() << "x"
            << m_rleIndex.m_abcCoordSize.zsize() << endl;
        throw std::runtime_error( "This function is for 2D level sets only." );
    }

    // output some standard RLE stuff
    out << "\nPrinting 2d level set:\n";

    out << "\tsize is " << m_rleIndex.m_abcCoordSize.xsize() << "x" << m_rleIndex.m_abcCoordSize.ysize() << "x"
        << m_rleIndex.m_abcCoordSize.zsize() << endl;

    out << "\torigin is "
        << m_rleIndex.m_abcCoordOrigin; //.xsize() << "x" << m_rleIndex.m_abcCoordOrigin.ysize() << "x"
                                        //<< m_rleIndex.m_abcCoordOrigin.zsize() << endl;

    out << "\n\tdistanceData size is " << m_distanceData.size() << endl;

    out << "\trunIndexData size is " << m_rleIndex.m_runData.size() << endl;

    out << "exterior is ";

    if( m_rleIndex.m_exteriorRegionCode == -1 ) {
        out << "outside\n";
    } else {
        out << "inside\n";
    }

    // for each run in the y dimension
    for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {

        out << b << ": ";

        int bcIndex = b;
        int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;
        int nextXToFill = ris.m_abcCoordOrigin.x, lastXToFill = nextXToFill + ris.m_abcCoordSize.xsize() - 1;

        int xStart, xPastTheEnd = ris.m_runData[runRangeStart].x;

        // go through the run
        for( int run = runRangeStart; run != runRangeEnd; ++run ) {
            xStart = xPastTheEnd;

            xPastTheEnd = ris.m_runData[run + 1].x;

            // If there are values before xStart, fill them up with "outside" values
            while( nextXToFill < xStart && nextXToFill <= lastXToFill ) {

                // print + for outside, - for inside
                if( ris.m_exteriorRegionCode == -1 ) {
                    out << "+ ";
                } else {
                    out << "- ";
                }
                ++nextXToFill;
            }

            int dataIndex = ris.m_runData[run].dataIndex;

            // print the run data
            if( dataIndex >= 0 ) {
                // If the run is defined, copy the distance values
                while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                    if( m_distanceData[dataIndex] < 0 ) {
                        out << m_distanceData[dataIndex++] << " ";
                    } else {
                        out << " " << m_distanceData[dataIndex++] << " ";
                    }
                    ++nextXToFill;
                }
            } else {
                // If the run is undefined, copy +/-levelSetOutsideDistance
                if( dataIndex == -1 ) {
                    while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                        out << "+ ";
                        ++nextXToFill;
                    }
                } else {
                    while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                        out << "- ";
                        ++nextXToFill;
                    }
                }
            }
        }
        // If there are any values to fill after the last run, fill them in with the exterior value
        while( nextXToFill <= lastXToFill ) {

            // print + for outside, - for inside
            if( ris.m_exteriorRegionCode == -1 ) {
                out << "+ ";
            } else {
                out << "- ";
            }
            ++nextXToFill;
        }
        // end of the run, print on a new line
        out << endl;
    }
}

void rle_level_set::dump( std::ostream& out ) const {
    out << "DUMPING RLE LEVEL SET\n";
    out << "Defined voxel count: " << m_rleIndex.data_size() << "\n";
    out << "Outer bounds: " << m_rleIndex.outer_bounds() << "\n";
    // m_rleIndex.dump(out);

    out << "Voxel Coordinate System: " << m_voxelCoordSystem << "\n";
    out << "Interface Voxel Width Inside: " << m_interfaceVoxelWidthInside << "\n";
    out << "Interface Voxel Width Outside: " << m_interfaceVoxelWidthOutside << "\n";
    out << "Outside Distance: " << m_outsideDistance << "\n";
    out << "Inside Distance: " << m_insideDistance << "\n";
    out << "\n";
    out << "Number of named channels: " << m_namedChannels.size() << "\n";
    int index = 0;
    for( map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.begin(); i != m_namedChannels.end();
         ++i ) {
        out << "Channel " << index++ << ", \"" << frantic::strings::to_string( i->first ) << "\"\n";
        out << "Data type: " << frantic::strings::to_string( i->second.type_str() ) << "\n";
        out << "Channel size: " << i->second.size() << "\n";
        // out << "Data: ";
        // channels::channel_data_type_print( out, ",", i->second.arity() * i->second.size(), i->second.data_type(),
        // i->second.data() ); out << "\n";
    }
    out << "FINISHED DUMPING RLE LEVEL SET" << endl;
}

void rle_level_set::getYZ( const rle_level_set& rle, int y, int z ) {
    rle_index_spec& ris = m_rleIndex;

    // set up the bounding box
    vector3 orig( rle.m_rleIndex.m_abcCoordOrigin.x, 0, 0 );
    size3 coordSize( rle.m_rleIndex.m_abcCoordSize.xsize(), 1, 1 );
    ris.m_abcCoordOrigin = orig;
    ris.m_abcCoordSize = coordSize;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
    ris.m_runData.clear();
    m_distanceData.clear();

    // set the exterior region
    ris.m_exteriorRegionCode = rle.m_rleIndex.m_exteriorRegionCode;

    // write the run data
    int bcIndex = ( y - rle.m_rleIndex.m_abcCoordOrigin.y ) +
                  ( z - rle.m_rleIndex.m_abcCoordOrigin.z ) * rle.m_rleIndex.m_abcCoordSize.ysize();
    ris.m_bcToRunIndex[0] = (int)ris.m_runData.size();

    int runRangeStart = rle.m_rleIndex.m_bcToRunIndex[bcIndex],
        runRangeEnd = rle.m_rleIndex.m_bcToRunIndex[bcIndex + 1] - 1;
    int nextXToFill = rle.m_rleIndex.m_abcCoordOrigin.x,
        lastXToFill = nextXToFill + rle.m_rleIndex.m_abcCoordSize.xsize() - 1;

    int xStart, xPastTheEnd = rle.m_rleIndex.m_runData[runRangeStart].x;

    // go through the run
    for( int run = runRangeStart; run != runRangeEnd + 1; ++run ) {

        if( run == runRangeEnd ) {
            ris.m_runData.push_back( run_data( rle.m_rleIndex.m_runData[run].x, ris.m_exteriorRegionCode ) );
            break;
        }

        xStart = xPastTheEnd;

        xPastTheEnd = rle.m_rleIndex.m_runData[run + 1].x;

        int dataIndex = rle.m_rleIndex.m_runData[run].dataIndex;
        if( dataIndex >= 0 ) {
            // If the run is defined, copy the distance values
            ris.m_runData.push_back( run_data( xStart, (int)m_distanceData.size() ) );
            while( xStart < xPastTheEnd && nextXToFill <= lastXToFill ) {
                m_distanceData.push_back( rle.m_distanceData[dataIndex++] );
                ++nextXToFill;
                ++xStart;
            }
        } else {
            // If the run is undefined, copy +/-levelSetOutsideDistance
            ris.m_runData.push_back( run_data( rle.m_rleIndex.m_runData[run].x, dataIndex ) );
            ++nextXToFill;
            //++dataIndex;
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = m_distanceData.size();
}

void rle_level_set::apply_axis_permutation( const vector3& axisPermutation ) {
    // Skip any processing if the permutation is the identity
    if( axisPermutation.x != 0 || axisPermutation.y != 1 || axisPermutation.z != 2 ) {
        // NOTE: this is not a vector, because vectors could copy the data and corrupt it
        deque<raw_byte_buffer> tempProcessingChannelData;
        vector<rle_index_spec_channel_copying_data> channelsToRemap;
        // Go through all the channels, and set up the channels to remap during the axis permutation
        for( map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.begin(); i != m_namedChannels.end();
             ++i ) {
            // Allocate a new buffer for this channel, and swap it so the destination is within the named channel data
            // structure
            tempProcessingChannelData.push_back( raw_byte_buffer() );
            tempProcessingChannelData.back().resize( i->second.m_data.size() );
            tempProcessingChannelData.back().swap( i->second.m_data );

            // Set up the entry for remapping this channel
            channelsToRemap.push_back( rle_index_spec_channel_copying_data() );
            channelsToRemap.back().primitiveSize = i->second.primitive_size();
            channelsToRemap.back().inputData = tempProcessingChannelData.back().begin();
            channelsToRemap.back().outputData = i->second.m_data.begin();

            //			cout << "channelToRemap[" << index++ << "] is the " << i->first << " channel " << endl;
            //			cout << "It has size " << i->second.size() << endl;
            //			cout << "Input has range " << (void*)tempProcessingChannelData.back().begin() << " to "
            //<< (void*)tempProcessingChannelData.back().end() << endl;
        }

        // Create the permutation for the distance data
        vector<float> tempDistanceData( m_distanceData.size() );
        if( m_distanceData.size() > 0 ) {
            tempDistanceData.swap( m_distanceData );
            channelsToRemap.push_back( rle_index_spec_channel_copying_data() );
            channelsToRemap.back().primitiveSize = sizeof( float );
            channelsToRemap.back().inputData = reinterpret_cast<char*>( &tempDistanceData[0] );
            channelsToRemap.back().outputData = reinterpret_cast<char*>( &m_distanceData[0] );
        }

        //		cout << "distance array has size " << tempDistanceData.size() << endl;

        // Now apply the axis permutation, providing the data to remap
        if( channelsToRemap.empty() )
            m_rleIndex.apply_axis_permutation( axisPermutation );
        else
            m_rleIndex.apply_axis_permutation( axisPermutation, &channelsToRemap[0],
                                               &channelsToRemap[0] + channelsToRemap.size() );
    }

    // run a consistency check (debugging code, leave it in for now)
    stringstream sout;
    if( !m_rleIndex.check_consistency( sout ) ) {
        throw std::runtime_error(
            "rle_level_set.apply_axis_permutation: RLE Index Spec consistency check failed (permutation=" +
            axisPermutation.str() + "):\n" + sout.str() );
    }
}

void rle_level_set::dilate_defined_voxels( int dilationVoxels, const frantic::tstring& populatedChannelName ) {
    //	boundbox3 originalOuterBounds = m_rleIndex.outer_bounds();
    rle_index_spec ris;
    // Create the dilated rle_index_spec
    ris.build_from_dilation( m_rleIndex, dilationVoxels );
    // Switch this level set to use the new rle_index_spec
    switch_rle_index_spec_with_swap( ris, populatedChannelName );
    // Copy channel data from the nearest voxel within the old outer bounds into the newly defined voxels
    //	extend_channels_from_original_outer_bounds(originalOuterBounds);
}

void rle_level_set::bounded_dilate_defined_voxels( int dilationVoxels, const frantic::graphics::boundbox3& voxelBounds,
                                                   const frantic::tstring& populatedChannelName ) {
    //	boundbox3 originalOuterBounds = m_rleIndex.outer_bounds();
    rle_index_spec ris0, ris1;
    // Create the dilated and trimmed rle_index_spec
    ris0.build_from_dilation( m_rleIndex, dilationVoxels );
    ris1.build_with_trim_bounds( ris0, voxelBounds );
    // Switch this level set to use the new rle_index_spec
    switch_rle_index_spec_with_swap( ris1, populatedChannelName );
    // Copy channel data from the nearest voxel within the old outer bounds into the newly defined voxels
    //	extend_channels_from_original_outer_bounds(originalOuterBounds);
}

void rle_level_set::dilate_defined_voxels( const boundbox3& dilationBox,
                                           const frantic::tstring& populatedChannelName ) {
    //	boundbox3 originalOuterBounds = m_rleIndex.outer_bounds();
    rle_index_spec ris;
    // Create the dilated rle_index_spec
    ris.build_from_dilation( m_rleIndex, dilationBox );
    // Switch this level set to use the new rle_index_spec
    switch_rle_index_spec_with_swap( ris, populatedChannelName );
    // Copy channel data from the nearest voxel within the old outer bounds into the newly defined voxels
    //	extend_channels_from_original_outer_bounds(originalOuterBounds);
}

void rle_level_set::extend_channels_from_original_outer_bounds( const boundbox3& originalOuterBounds ) {
    // If the new outer bounds aren't bigger than the old ones, there's nothing to process
    if( originalOuterBounds.contains( m_rleIndex.outer_bounds() ) )
        return;

    // TODO: This function should copy the defined voxel values from the edge of the original outer bounds outwards
    //       to the newly defined voxels.  Three stages:
    //       1) Go through all the runs within the original bounds YZ plane, and if a run crosses the edge of the
    //          original outer bounds, copy the value on the edge outwards for all channels.  It's unclear what to do
    //          if there's a gap and then more defined voxels outside the bounding box extents.  Probably want to keep
    //          copying.
    //       2) For all the runs within the original bounds Z extents but outside the Y extents, use a pairwise
    //          run iterator to copy all the channel data from the nearest scanline within the bound box along the Y
    //          direction.
    //       3) For all the runs that are outside the original bounds Z extents, use a pairwise
    //          run iterator to copy all the channel data from the nearest scanline within the bound box along the Z
    //          direction.
    //
    //       This should suffice to extend the channel data into the newly dilated voxels that extended beyond the
    //       bounding box, or that were revealed based on switching the rle index spec.  The idea here is that we can
    //       use a strategy of extending a voxel based on the nearest voxel in the previous bounding box rather than
    //       just defaulting to the exterior region code.
    //
    //       Another useful function to write is one that grows the outer bounds of a level set, extending defined
    //       voxels or undefined voxels at the edge, then using this function to fill in those newly defined voxels.
    //       This would allow for the LS Mesher to have an "intersect with bounds" mode, where you extend the bounding
    //       box with this nearest voxel copying strategy, then create a level set of the voxel bounds, do a level set
    //       intersection, and finally mesh the result.  This will produce a much higher quality result than the current
    //       result when you uncheck "Clip Bounds."

    throw runtime_error( "rle_level_set.extend_channels_from_original_outer_bounds() - TODO: implement me!" );
}

void rle_level_set::trim_to_populated( const frantic::tstring& populatedChannelName ) {
    if( m_rleIndex.data_size() == 0 )
        return;

    // Make sure the populated channel exists
    if( !has_channel( populatedChannelName ) )
        throw runtime_error(
            "rle_level_set::trim_to_populated() - There must be a \"" +
            frantic::strings::to_string( populatedChannelName ) +
            "\" channel defined in the level set specifying which voxels are populated and which should be trimmed." );

    // Make sure its a uint8 channel
    rle_channel_general_accessor populatedAccessor = get_channel_general_accessor( populatedChannelName );
    if( populatedAccessor.arity() != 1 || populatedAccessor.data_type() != channels::data_type_uint8 )
        throw runtime_error( "rle_level_set::trim_to_populated() - The  \"" +
                             frantic::strings::to_string( populatedChannelName ) +
                             "\" channel (for the 'Populated' data) must be a uint8 channel, it is set to " +
                             populatedAccessor.type_str() + "." );

    // Get a pointer to the channel data
    unsigned char* populatedChannel = reinterpret_cast<unsigned char*>( populatedAccessor.data( 0 ) );

    rle_index_spec ris;
    // Use the short-circuiting to do nothing if the populatedChannel indicates nothing is unpopulated.
    if( ris.build_from_populated( m_rleIndex, populatedChannel, &m_distanceData[0],
                                  /*shortCircuitFullyPopulated=*/true ) ) {
        // No point in keeping the populated channel around, it will generally just be all ones now
        erase_channel( populatedChannelName );
        switch_rle_index_spec_with_swap( ris, populatedChannelName );
    } else {
        // No point in keeping the populated channel around, it will generally just be all ones now
        erase_channel( populatedChannelName );
    }
}

void rle_level_set::trim_to_bounds( const boundbox3& trimVoxelBounds ) {
    // Only do the trimming if the outer bounds are bigger than the trim bounds
    if( !m_rleIndex.outer_bounds().is_empty() && !trimVoxelBounds.contains( m_rleIndex.outer_bounds() ) ) {
        // To do the trimming, create a temp trimmed rle index spec, then switch the level set using
        // a swap on the rle index spec for efficiency.
        rle_index_spec risTemp;
        risTemp.build_with_trim_bounds( m_rleIndex, trimVoxelBounds );
        switch_rle_index_spec_with_swap( risTemp );
    }
}

void rle_level_set::compute_upwind_gradient( const frantic::tstring& channelName,
                                             const frantic::tstring& gradientChannelToCreate,
                                             const frantic::tstring& populatedChannelToCreate, int upwindDir ) {
    if( upwindDir != -1 && upwindDir != +1 )
        throw runtime_error( "rle_level_set::compute_upwind_gradient() - The upwind direction must be -1 or +1,"
                             "the value provided is " +
                             boost::lexical_cast<string>( upwindDir ) + "." );
    // Convert it to a float.
    float distanceSignMultiplier = (float)upwindDir;

    const_rle_channel_general_accessor chanAcc = get_channel_general_accessor( channelName );

    if( !channels::is_channel_data_type_float( chanAcc.data_type() ) )
        throw runtime_error( "rle_level_set::compute_upwind_gradient() - The channel \"" +
                             frantic::strings::to_string( channelName ) + "\" for the gradient has data type " +
                             chanAcc.type_str() + ", which is not a floating point type." );

    // Create the output gradient channel
    add_channel( gradientChannelToCreate, chanAcc.arity() * 3, chanAcc.data_type() );
    rle_channel_general_accessor gradChanAcc = get_channel_general_accessor( gradientChannelToCreate );

    rle_channel_accessor<unsigned char> populatedAcc;
    if( !populatedChannelToCreate.empty() ) {
        add_channel( populatedChannelToCreate, 1, data_type_uint8 );
        populatedAcc = get_channel_accessor<unsigned char>( populatedChannelToCreate );
    }

    // If there are no defined voxels, then exit after creating the output field, there is no data to compute.
    if( m_rleIndex.data_size() == 0 )
        return;

    channel_weighted_sum_combine_function_t combineFn = chanAcc.get_weighted_sum_combine_function();
    size_t arity = chanAcc.arity();
    size_t primitiveSize = chanAcc.primitive_size();

    float weights[2] = { -1.f / m_voxelCoordSystem.voxel_length(), +1.f / m_voxelCoordSystem.voxel_length() };
    const char* data[2] = { 0, 0 };

    const float* signedDistance = &m_distanceData[0];

    // This mask is used to accumulate which axes actually have any partial derivatives.  If a field is two dimensional,
    // one of its bits will not be set, then we can use it to set voxels to populated which have partial derivatives for
    // the other two directions.
    unsigned char fullPopulatedMask = 0;

    for( rle_defined_and_adj_iterator i( m_rleIndex ), ie( m_rleIndex, true ); i != ie; ++i ) {
        const char* centerData = chanAcc.data( i.get_center_data_index() );
        char* gradientData = gradChanAcc.data( i.get_center_data_index() );
        float centerDistance = distanceSignMultiplier * signedDistance[i.get_center_data_index()];
        unsigned char populated = 7;
        // X partial derivative
        if( i.get_x_neg_data_index() >= 0 &&
            distanceSignMultiplier * signedDistance[i.get_x_neg_data_index()] >= centerDistance ) {
            if( i.get_x_pos_data_index() >= 0 &&
                distanceSignMultiplier * signedDistance[i.get_x_pos_data_index()] >=
                    distanceSignMultiplier * signedDistance[i.get_x_neg_data_index()] ) {
                // upwind is towards positive X
                data[0] = centerData;
                data[1] = chanAcc.data( i.get_x_pos_data_index() );
                combineFn( weights, data, 2, arity, gradientData );
            } else {
                // upwind is towards negative X
                data[0] = chanAcc.data( i.get_x_neg_data_index() );
                data[1] = centerData;
                combineFn( weights, data, 2, arity, gradientData );
            }
        } else if( i.get_x_pos_data_index() >= 0 &&
                   distanceSignMultiplier * signedDistance[i.get_x_pos_data_index()] >= centerDistance ) {
            // upwind is towards positive X
            data[0] = centerData;
            data[1] = chanAcc.data( i.get_x_pos_data_index() );
            combineFn( weights, data, 2, arity, gradientData + 0 );
        } else {
            // there is no upwind direction in X
            memset( gradientData, 0, primitiveSize );
            populated -= 1;
        }
        // Move the gradient data pointer from the X to the Y component
        gradientData += primitiveSize;
        // Y partial derivative
        if( i.get_y_neg_data_index() >= 0 &&
            distanceSignMultiplier * signedDistance[i.get_y_neg_data_index()] >= centerDistance ) {
            if( i.get_y_pos_data_index() >= 0 &&
                distanceSignMultiplier * signedDistance[i.get_y_pos_data_index()] >=
                    distanceSignMultiplier * signedDistance[i.get_y_neg_data_index()] ) {
                // upwind is towards positive Y
                data[0] = centerData;
                data[1] = chanAcc.data( i.get_y_pos_data_index() );
                combineFn( weights, data, 2, arity, gradientData );
            } else {
                // upwind is towards negative Y
                data[0] = chanAcc.data( i.get_y_neg_data_index() );
                data[1] = centerData;
                combineFn( weights, data, 2, arity, gradientData );
            }
        } else if( i.get_y_pos_data_index() >= 0 &&
                   distanceSignMultiplier * signedDistance[i.get_y_pos_data_index()] >= centerDistance ) {
            // upwind is towards positive Y
            data[0] = centerData;
            data[1] = chanAcc.data( i.get_y_pos_data_index() );
            combineFn( weights, data, 2, arity, gradientData );
        } else {
            // there is no upwind direction in Y
            memset( gradientData, 0, primitiveSize );
            populated -= 2;
        }
        // Move the gradient data pointer from the Y to the Z component
        gradientData += primitiveSize;
        // Z partial derivative
        if( i.get_z_neg_data_index() >= 0 &&
            distanceSignMultiplier * signedDistance[i.get_z_neg_data_index()] >= centerDistance ) {
            if( i.get_z_pos_data_index() >= 0 &&
                distanceSignMultiplier * signedDistance[i.get_z_pos_data_index()] >=
                    distanceSignMultiplier * signedDistance[i.get_z_neg_data_index()] ) {
                // upwind is towards positive Z
                data[0] = centerData;
                data[1] = chanAcc.data( i.get_z_pos_data_index() );
                combineFn( weights, data, 2, arity, gradientData );
            } else {
                // upwind is towards negative Z
                data[0] = chanAcc.data( i.get_z_neg_data_index() );
                data[1] = centerData;
                combineFn( weights, data, 2, arity, gradientData );
            }
        } else if( i.get_z_pos_data_index() >= 0 &&
                   distanceSignMultiplier * signedDistance[i.get_z_pos_data_index()] >= centerDistance ) {
            // upwind is towards positive Z
            data[0] = centerData;
            data[1] = chanAcc.data( i.get_z_pos_data_index() );
            combineFn( weights, data, 2, arity, gradientData );
        } else {
            // there is no upwind direction in Z
            memset( gradientData, 0, primitiveSize );
            populated -= 4;
        }

        // Set populated to to the mask indicating which directions had a valid gradient
        if( populatedAcc.valid() ) {
            populatedAcc[i.get_center_data_index()] = populated;
            fullPopulatedMask |= populated;
        }
    }

    // Now go through the populated channel, and mark as populated any voxel that had the maximal number of partial
    // derivatives we saw. This way, a 1D or 2D field will still get the appropriate cells marked as populated.
    if( populatedAcc.valid() ) {
        for( size_t i = 0, ie = populatedAcc.size(); i != ie; ++i ) {
            if( populatedAcc[i] == fullPopulatedMask )
                populatedAcc[i] = 1;
            else
                populatedAcc[i] = 0;
        }
    }
}

void rle_level_set::fill_plane( const boundrect2& voxelXYExtents, int voxelZ, float* outVoxelCornerValues ) const {
    fill_box( boundbox3( vector3( voxelXYExtents.minimum().x, voxelXYExtents.minimum().y, voxelZ ),
                         vector3( voxelXYExtents.maximum().x, voxelXYExtents.maximum().y, voxelZ ) ),
              outVoxelCornerValues );
}

void rle_level_set::fill_plane( const boundrect2& voxelXYExtents, int voxelZ,
                                std::vector<float>& outVoxelCornerValues ) const {
    fill_box( boundbox3( vector3( voxelXYExtents.minimum().x, voxelXYExtents.minimum().y, voxelZ ),
                         vector3( voxelXYExtents.maximum().x, voxelXYExtents.maximum().y, voxelZ ) ),
              outVoxelCornerValues );
}

/**
 *	This function fills a z plane sparsely with data from the level set, and a rle plane
 *	structure to sparsely access it.
 *
 *	@param	voxelXYExtents	The xy extents of the plane to be populated.
 *	@param	voxelZ			The z value of the plane.
 *	@param	channelNames	A vector of strings of channel names to be fetched.
 *	@param	channelData		A vector of char* to arrays of channel data to be populated.
 *	@param	outRLP			The rle plane structure to access the output data sparsely.
 */
void rle_level_set::fill_sparse_plane_channel_data( const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                                                    std::vector<frantic::tstring>& channelNames,
                                                    std::vector<char*>& channelData,
                                                    frantic::volumetrics::rle_plane& outRLP ) const {
    const rle_index_spec& ris = m_rleIndex;

    // build an array of rle_channel accessors and the output channels they match to
    std::vector<std::pair<const_rle_channel_general_accessor, char*>> channelPairs;
    int signedDistanceChannel = -1;
    for( size_t i = 0; i < channelNames.size(); ++i ) {
        if( channelNames[i] == _T("SignedDistance") )
            signedDistanceChannel = (int)i;
        else
            channelPairs.push_back( std::pair<const_rle_channel_general_accessor, char*>(
                get_channel_general_accessor( channelNames[i] ), channelData[i] ) );
    }

    // rest of the channels
    outRLP.reset( voxelXYExtents );

    // This is the run vector that we build for each extent
    vector<pair<int, int>> runs;
    vector<int> runCodes;

    // Get the B coordinate range and the C coordinate
    int bMin = voxelXYExtents.minimum().y - ris.m_abcCoordOrigin.y;
    int bMax = voxelXYExtents.maximum().y - ris.m_abcCoordOrigin.y;
    int c = voxelZ - ris.m_abcCoordOrigin.z; // the c coord is just the provided z coord
    int minX = voxelXYExtents.minimum().x, maxX = voxelXYExtents.maximum().x;

    // If the C coordinate or the B coordinate interval is outside the range of the bounding box,
    // or the x min/max are outside of the range of the run data.  There will be no defined run data,
    // and we can just fill the plane with the exterior region code and return.
    if( bMax < 0 || bMin >= ris.m_abcCoordSize.ysize() || c < 0 || c >= ris.m_abcCoordSize.zsize() ||
        maxX < ris.m_abcCoordOrigin.x || minX >= ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() ) {
        for( int i = 0; i < voxelXYExtents.ysize(); i++ )
            outRLP.append_runs_by_extent( runs, i, ris.m_exteriorRegionCode );
        return;
    }

    // These are the number of extents to fill in before copying values within the defined B range of the bounding box
    int bMinExteriorCount = 0, bMaxExteriorCount = 0;

    // If the array desired extends in front of the defined B range, we will have to fill the first part with exterior
    // values.
    if( bMin < 0 ) {
        bMinExteriorCount = -bMin;
        bMin = 0;
    }

    // Fill in any leading empty scanlines
    int extent = 0;
    for( extent = 0; extent < bMinExteriorCount; ++extent )
        outRLP.append_runs_by_extent( runs, extent, ris.m_exteriorRegionCode );

    // If the array desired extends beyond the defined B range, we will have to fill the last part with
    // exterior values
    if( bMax >= ris.m_abcCoordSize.ysize() ) {
        bMaxExteriorCount = bMax - ris.m_abcCoordSize.ysize() + 1;
        bMax = m_rleIndex.m_abcCoordSize.ysize() - 1;
    }

    // Now iterate through all the B coordinates and copy the scanlines
    for( int b = bMin; b <= bMax; ++b ) {

        // Get the run range
        int bcIndex = b + c * ris.m_abcCoordSize.ysize();
        int nextXToFill = minX, lastXToFill = maxX;
        int runRangeStart = ris.BCIndexXtoRunIndex( bcIndex, nextXToFill );
        int runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

        // This is current x-extent (y coordinate) of the plane we are traversing
        extent = ris.m_abcCoordOrigin.y + b - voxelXYExtents.minimum().y;
        int extentStartIndex = extent * voxelXYExtents.xsize() - voxelXYExtents.minimum().x;
        // Handle the empty scanline case.
        if( ris.m_runData[runRangeStart].x == ris.m_runData[runRangeEnd].x ||
            ris.m_runData[runRangeStart].x > lastXToFill || nextXToFill >= ris.m_runData[runRangeEnd].x ) {
            outRLP.append_runs_by_extent( runs, extent, ris.m_exteriorRegionCode );
            continue;
        }

        int xStart, xPastTheEnd = ris.m_runData[runRangeStart].x;
        int dataIndex = 0;

        for( int runIndex = runRangeStart; runIndex < runRangeEnd; ++runIndex ) {

            // Get the x range of this run
            xStart = xPastTheEnd;
            xPastTheEnd = ris.m_runData[runIndex + 1].x;

            dataIndex = ris.m_runData[runIndex].dataIndex;

            if( nextXToFill < xStart )
                nextXToFill = xStart;

            if( dataIndex >= 0 ) {
                // If the run is defined, copy the channel data values and add the run index info
                dataIndex += nextXToFill - ris.m_runData[runIndex].x;
                std::pair<int, int> run( nextXToFill + extentStartIndex,
                                         min( xPastTheEnd - 1, lastXToFill ) + extentStartIndex );
                int runLength = run.second - run.first + 1;

                // signed distance data
                if( signedDistanceChannel >= 0 )
                    memcpy( channelData[signedDistanceChannel] + sizeof( float ) * run.first,
                            &m_distanceData[dataIndex], runLength * sizeof( float ) );

                // rest of the channel data
                for( size_t i = 0; i < channelPairs.size(); ++i )
                    memcpy( channelPairs[i].second + channelPairs[i].first.primitive_size() * run.first,
                            channelPairs[i].first.data( dataIndex ),
                            runLength * channelPairs[i].first.primitive_size() );
                nextXToFill += runLength;

                runs.push_back( run );
                runCodes.push_back( 0 );
            } else {
                runCodes.push_back( dataIndex );
                runs.push_back( std::pair<int, int>( nextXToFill + extentStartIndex,
                                                     min( lastXToFill, xPastTheEnd - 1 ) + extentStartIndex ) );
            }

            // Break out if our current run goes past the maxX of the plane we are populating
            if( xPastTheEnd > lastXToFill )
                break;
        }

        // If we had only a single undefined run, insert an empty extent with the proper code
        if( runs.size() == 0 && dataIndex < 0 ) {
            runs.clear();
            runCodes.clear();
            outRLP.append_runs_by_extent( runs, extent, dataIndex );
        }

        // Otherwise add the runs for this extent to the indexing structure
        else {
            outRLP.append_runs_by_extent( runs, runCodes, extent );
            runs.clear();
            runCodes.clear();
        }
    }

    // Fill in any following empty scanlines
    for( int i = 0; i < bMaxExteriorCount; i++ )
        outRLP.append_runs_by_extent( runs, extent + i + 1, ris.m_exteriorRegionCode );
}

// This function has only been tested as used through fill_plane so far.  Need to validate that it works in more general
// settings.
void rle_level_set::fill_box( const frantic::graphics::boundbox3& voxelExtents, float* outVoxelCornerValues ) const {
    if( !voxelExtents.is_empty() && outVoxelCornerValues == NULL ) {
        throw std::runtime_error( "rle_level_set.fill_box Error: the output array is NULL." );
    }

    ////////////////
    // INITIALIZE
    ////////////////

    const rle_index_spec& ris = m_rleIndex;

    // The distance we use for "outside" values is 1 more than the number of voxels the interface radius is
    float levelSetOutsideDistance = m_outsideDistance;
    // The rle_index_spec specifies a region code for all the voxels exterior to the surface.  This determines the level
    // set value to copy for that case.
    float levelSetExteriorDistance =
        ( ris.m_exteriorRegionCode == -1 ) ? levelSetOutsideDistance : -levelSetOutsideDistance;

    // Get the B coordinate range and the C coordinate
    int bMin = voxelExtents.minimum().y - ris.m_abcCoordOrigin.y,
        bMax = voxelExtents.maximum().y - ris.m_abcCoordOrigin.y;
    int cMin = voxelExtents.minimum().z - ris.m_abcCoordOrigin.z,
        cMax = voxelExtents.maximum().z - ris.m_abcCoordOrigin.z;

    // If the C coordinate interval or the B coordinate interval is outside the range of the bounding box, set the whole
    // array to "outside"
    if( bMax < 0 || bMin >= ris.m_abcCoordSize.ysize() || cMax < 0 || cMin >= ris.m_abcCoordSize.zsize() ) {
        for( int i = 0, ie = voxelExtents.get_volume(); i != ie; ++i ) {
            outVoxelCornerValues[i] = levelSetExteriorDistance;
        }
        return;
    }

    // These are the number of voxels to fill in before copying values withing the defined B range of the bounding box
    int bMinExteriorCount = 0, bMaxExteriorCount = 0;

    // If the array desired extends in front of the defined B range, we will have to fill the first part with exterior
    // values.
    if( bMin < 0 ) {
        bMinExteriorCount = -bMin * voxelExtents.xsize();
        bMin = 0;
    }

    // If the array desired extends beyond the defined B range, we will have to fill the last part with exterior values
    if( bMax >= ris.m_abcCoordSize.ysize() ) {
        bMaxExteriorCount = ( bMax - ris.m_abcCoordSize.ysize() + 1 ) * voxelExtents.xsize();
        bMax = m_rleIndex.m_abcCoordSize.ysize() - 1;
    }

    ////////////////
    // COPY LEVEL SET VALUES
    ////////////////

    // This is the value we use to step through the outVoxelCornerValues
    unsigned outputIndex = 0;

    // If the array desired extends in front of the defined C range, fill the first part with exterior values.
    if( cMin < 0 ) {
        int blankValueCount = -cMin * voxelExtents.xsize() * voxelExtents.ysize();
        for( int i = 0; i < blankValueCount; ++i )
            outVoxelCornerValues[outputIndex++] = levelSetExteriorDistance;
        cMin = 0;
    }

    // Trim the C values we iterate over to within the range
    if( cMax >= ris.m_abcCoordSize.zsize() ) {
        cMax = m_rleIndex.m_abcCoordSize.zsize() - 1;
    }

    // Now iterate through all the C coordinates and copy the planes
    for( int c = cMin; c <= cMax; ++c ) {
        // Fill in the starting B out-of-bounds scanlines with exterior values
        for( int i = 0; i < bMinExteriorCount; ++i )
            outVoxelCornerValues[outputIndex++] = levelSetExteriorDistance;

        // Now iterate through all the B coordinates and copy the scanlines
        for( int b = bMin; b <= bMax; ++b ) {
            // Get the run range
            int bcIndex = b + c * ris.m_abcCoordSize.ysize();
            int nextXToFill = voxelExtents.minimum().x, lastXToFill = voxelExtents.maximum().x;
            int runRangeStart = ris.BCIndexXtoRunIndex( bcIndex, nextXToFill ),
                runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

            int xStart, xPastTheEnd = ris.m_runData[runRangeStart].x;
            for( int run = runRangeStart; run != runRangeEnd; ++run ) {
                // Get the x range of this run
                xStart = xPastTheEnd;
                xPastTheEnd = ris.m_runData[run + 1].x;

                // If there are values before xStart, fill them up with "outside" values
                while( nextXToFill < xStart && nextXToFill <= lastXToFill ) {
                    outVoxelCornerValues[outputIndex++] = levelSetExteriorDistance;
                    ++nextXToFill;
                }

                int dataIndex = ris.m_runData[run].dataIndex;
                if( dataIndex >= 0 ) {
                    dataIndex += nextXToFill - ris.m_runData[run].x;
                    // If the run is defined, copy the distance values
                    while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                        outVoxelCornerValues[outputIndex++] = m_distanceData[dataIndex++];
                        ++nextXToFill;
                    }
                } else {
                    // If the run is undefined, copy +/-levelSetOutsideDistance
                    while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                        // Region code -1 indicates "outside", anything else indicates "inside".
                        outVoxelCornerValues[outputIndex++] =
                            dataIndex == -1 ? levelSetOutsideDistance : -levelSetOutsideDistance;
                        ++nextXToFill;
                    }
                }
            }
            // If there are any values to fill after the last run, fill them in with the exterior value
            while( nextXToFill <= lastXToFill ) {
                outVoxelCornerValues[outputIndex++] = levelSetExteriorDistance;
                ++nextXToFill;
            }
        }
        // Fill in the ending B out-of-bounds scanlines with exterior values
        for( int i = 0; i < bMaxExteriorCount; ++i )
            outVoxelCornerValues[outputIndex++] = levelSetExteriorDistance;
    }

    // If the array extended beyond in the C direction, we need to fill the rest with exterior values
    const unsigned int outputIndexEnd = voxelExtents.get_volume();
    // while( outputIndex < outVoxelCornerValues.size() )
    while( outputIndex < outputIndexEnd )
        outVoxelCornerValues[outputIndex++] = levelSetExteriorDistance;
}

// This function has only been tested as used through fill_plane so far.  Need to validate that it works in more general
// settings.
void rle_level_set::fill_box( const frantic::graphics::boundbox3& voxelExtents,
                              std::vector<float>& outVoxelCornerValues ) const {
    ////////////////
    // INITIALIZE
    ////////////////

    // Make sure the output array is the right size
    outVoxelCornerValues.resize( voxelExtents.get_volume() );

    if( outVoxelCornerValues.size() > 0 ) {
        fill_box( voxelExtents, &outVoxelCornerValues[0] );
    }
}

// A functor which orders vector3 values in the rle_index_spec sort order.
// NOTE: This is the same sorted order as in rle_index_spec.cpp rle_index_spec_sort_order class.
class rle_index_spec_sort_order_pair {
  public:
    bool operator()( const pair<vector3, float>& lhs, const pair<vector3, float>& rhs ) const {
        if( lhs.first.z != rhs.first.z )
            return lhs.first.z < rhs.first.z;
        else if( lhs.first.y != rhs.first.y )
            return lhs.first.y < rhs.first.y;
        else
            return lhs.first.x < rhs.first.x;
    }
};

// This sets the rle level set to the given data
void rle_level_set::set( const voxel_coord_system& vcs, const rle_index_spec& ris,
                         const std::vector<float>& distanceData, float interfaceVoxelWidthInside,
                         float interfaceVoxelWidthOutside ) {
    if( ris.data_size() != distanceData.size() )
        throw runtime_error( "rle_level_set.set: The size of the data provided, " +
                             boost::lexical_cast<std::string>( distanceData.size() ) +
                             ", doesn't match the data size of the RLE Index Spec provided, " +
                             boost::lexical_cast<std::string>( ris.data_size() ) + "." );

    m_rleIndex = ris;
    m_voxelCoordSystem = vcs;
    m_interfaceVoxelWidthInside = interfaceVoxelWidthInside;
    m_interfaceVoxelWidthOutside = interfaceVoxelWidthOutside;
    m_outsideDistance = ( m_interfaceVoxelWidthOutside + 1 ) * m_voxelCoordSystem.voxel_length();
    m_insideDistance = -( m_interfaceVoxelWidthInside + 1 ) * m_voxelCoordSystem.voxel_length();
    m_distanceData = distanceData;
    m_namedChannels.clear();
}

// This sets the rle level set to the given data, using swap functions to set the data efficiently.  Note that
// the variables you pass in will have arbitrary values in them after this function is called.
void rle_level_set::set_with_swap( voxel_coord_system& vcs, rle_index_spec& ris, std::vector<float>& distanceData,
                                   float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside ) {
    if( ris.data_size() != distanceData.size() )
        throw runtime_error( "rle_level_set.set_with_swap: The size of the data provided, " +
                             boost::lexical_cast<std::string>( distanceData.size() ) +
                             ", doesn't match the data size of the RLE Index Spec provided, " +
                             boost::lexical_cast<std::string>( ris.data_size() ) + "." );

    m_rleIndex.swap( ris );
    m_voxelCoordSystem.swap( vcs );
    m_interfaceVoxelWidthInside = interfaceVoxelWidthInside;
    m_interfaceVoxelWidthOutside = interfaceVoxelWidthOutside;
    m_outsideDistance = ( m_interfaceVoxelWidthOutside + 1 ) * m_voxelCoordSystem.voxel_length();
    m_insideDistance = -( m_interfaceVoxelWidthInside + 1 ) * m_voxelCoordSystem.voxel_length();
    m_distanceData.swap( distanceData );
    m_namedChannels.clear();
}

// Convenience method for creating a level set, based on a map
void rle_level_set::set( const std::map<vector3, float>& mapLevelSet, float interfaceVoxelWidthInside,
                         float interfaceVoxelWidthOutside ) {
    // Copy the pairs into a vector, and sort them in the rle_index_spec sort order
    vector<pair<vector3, float>> valuesCopy;
    std::copy( mapLevelSet.begin(), mapLevelSet.end(), std::back_inserter( valuesCopy ) );
    sort( valuesCopy.begin(), valuesCopy.end(), rle_index_spec_sort_order_pair() );

    // Copy the voxel coordinates into their own array, and build the index spec from them
    vector<vector3> voxels( valuesCopy.size() );
    for( unsigned i = 0; i < valuesCopy.size(); ++i )
        voxels[i] = valuesCopy[i].first;
    m_rleIndex.build_from_voxel_array( voxels );

    // Copy the data coordinates into the distances array
    m_distanceData.resize( valuesCopy.size() );
    for( unsigned i = 0; i < valuesCopy.size(); ++i )
        m_distanceData[i] = valuesCopy[i].second;

    m_interfaceVoxelWidthInside = interfaceVoxelWidthInside;
    m_interfaceVoxelWidthOutside = interfaceVoxelWidthOutside;
    m_outsideDistance = ( m_interfaceVoxelWidthOutside + 1 ) * m_voxelCoordSystem.voxel_length();
    m_insideDistance = -( m_interfaceVoxelWidthInside + 1 ) * m_voxelCoordSystem.voxel_length();

    // Clear all the named channels
    m_namedChannels.clear();
}

// Set the current run to empty
void rle_level_set::set_to_empty( boost::int32_t exteriorRegionCode ) {
    m_rleIndex.clear();

    // set the exterior region
    m_rleIndex.m_exteriorRegionCode = exteriorRegionCode;

    // Clear all the channels to zero-size
    m_distanceData.clear();
    for( map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.begin(); i != m_namedChannels.end(); ++i ) {
        i->second.m_data.clear();
    }

    // Indicate that we're really far away from the surface (possibly use 1e38 instead?)
    m_outsideDistance = 10000 * m_voxelCoordSystem.voxel_length();
    m_insideDistance = -10000 * m_voxelCoordSystem.voxel_length();
}

float rle_level_set::compute_volume() const {
    if( m_rleIndex.get_exterior_region_code() != -1 )
        throw runtime_error(
            "rle_level_set::compute_volume() - The exterior region code is \"inside\", so the volume is "
            "effectively infinite." );

    if( size() == 0 ||
        m_rleIndex.m_runData.size() < 3 ) { // if there aren't enough runs to have at least one defined run
        return 0.0f;
    }

    float voxelLength = m_voxelCoordSystem.voxel_length();
    float volumeOfVoxel = voxelLength * voxelLength * voxelLength;

    float totalVolume = 0;
    float halfVoxelLength = voxelLength * 0.5f;

    boundbox3 outerBounds = m_rleIndex.outer_bounds();
    for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
        for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
            // Iterate through all the runs in this scanline
            for( rle_run_iterator i( m_rleIndex, m_rleIndex.y_to_b( y ), m_rleIndex.z_to_c( z ) ), ie; i != ie; ++i ) {
                if( i.get_data_index() >= 0 ) {
                    // Loop through all the distance values in this defined run
                    for( const float *dist = &m_distanceData[i.get_data_index()], *distEnd = dist + i.get_xsize();
                         dist != distEnd; ++dist ) {
                        float phi = *dist;
                        if( phi < -halfVoxelLength ) { // the voxel is completely inside
                            totalVolume += volumeOfVoxel;
                        } else if( phi < halfVoxelLength ) { // the voxel is partially inside
                            totalVolume += ( halfVoxelLength - phi ) * voxelLength * voxelLength;
                        } // otherwise there isn't any fluid in the voxel
                    }
                } else if( i.get_data_index() < -1 ) {
                    // An "inside" run, all voxels are inside
                    totalVolume += volumeOfVoxel * i.get_xsize();
                }
            }
        }
    }

    return totalVolume;
}

/**
 * This function resamples a given rle level set into the current voxel coordinate system. It is buggy, See Comment
 *Notes!
 **/
void rle_level_set::resample( rle_level_set& rleSource,
                              frantic::graphics::transform4f transformNoScaleSourceToCurrent ) {
    // ensure that the rotation matrix does not scale (allow for round-off errors though)
    float xScale = sqrtf( transformNoScaleSourceToCurrent.get( 0, 0 ) * transformNoScaleSourceToCurrent.get( 0, 0 ) +
                          transformNoScaleSourceToCurrent.get( 1, 0 ) * transformNoScaleSourceToCurrent.get( 1, 0 ) +
                          transformNoScaleSourceToCurrent.get( 2, 0 ) * transformNoScaleSourceToCurrent.get( 2, 0 ) );
    if( xScale < 1 - 1e-4 || xScale > 1 + 1e-4 ) {
        throw std::runtime_error( "rle_level_set.resample: the input transform matrix contains a scaling factor of " +
                                  boost::lexical_cast<string>( xScale ) +
                                  ". Please input only a rotation transform matrix (scale can be defined in the "
                                  "voxel_length of the current voxel coordinate system)." );
    }
    float yScale = sqrtf( transformNoScaleSourceToCurrent.get( 0, 1 ) * transformNoScaleSourceToCurrent.get( 0, 1 ) +
                          transformNoScaleSourceToCurrent.get( 1, 1 ) * transformNoScaleSourceToCurrent.get( 1, 1 ) +
                          transformNoScaleSourceToCurrent.get( 2, 1 ) * transformNoScaleSourceToCurrent.get( 2, 1 ) );
    if( yScale < 1 - 1e-4 || yScale > 1 + 1e-4 ) {
        throw std::runtime_error( "rle_level_set.resample: the input transform matrix contains a scaling factor of " +
                                  boost::lexical_cast<string>( yScale ) +
                                  ". Please input only a rotation transform matrix (scale can be defined in the "
                                  "voxel_length of the current voxel coordinate system)." );
    }
    float zScale = sqrtf( transformNoScaleSourceToCurrent.get( 0, 2 ) * transformNoScaleSourceToCurrent.get( 0, 2 ) +
                          transformNoScaleSourceToCurrent.get( 1, 2 ) * transformNoScaleSourceToCurrent.get( 1, 2 ) +
                          transformNoScaleSourceToCurrent.get( 2, 2 ) * transformNoScaleSourceToCurrent.get( 2, 2 ) );
    if( zScale < 1 - 1e-4 || zScale > 1 + 1e-4 ) {
        throw std::runtime_error( "rle_level_set.resample: the input transform matrix contains a scaling factor of " +
                                  boost::lexical_cast<string>( zScale ) +
                                  ". Please input only a rotation transform matrix (scale can be defined in the "
                                  "voxel_length of the current voxel coordinate system)." );
    }

    rle_index_spec& ris = m_rleIndex;

    //	float voxelLength = m_voxelCoordSystem.voxel_length();
    vector3f worldOrigin = m_voxelCoordSystem.world_origin();

    // The relative transforms don't properly take into account that the voxel centers are the sample locations,
    // but the voxel coordinates are at the corner of the voxel. This results in the output levelset being offset
    // by half a voxel length in every direction.
    // We tried to fix this by creating a voxel center to world transform and the inverse of that. But they don't seem
    // to work
    transform4f xfrmCurentVoxelToWorld = m_voxelCoordSystem.voxel_center_to_world_transform();

    transform4f xfrmSourceWorldToVoxel = rleSource.get_voxel_coord_system().world_voxel_center_to_voxel_transform();
    // cout << "transformSourceWorldToVoxel:\n" << xfrmSourceWorldToVoxel << endl;
    // cout << "transformCurrentVoxelToWorld:\n" << xfrmCurentVoxelToWorld << endl;

    transform4f transformCurrentToSource =
        xfrmSourceWorldToVoxel * transformNoScaleSourceToCurrent * xfrmCurentVoxelToWorld;

    // cout << "transformSourceWorldToVoxel^-1:\n" << xfrmSourceWorldToVoxel.to_inverse() << endl;

    // cout << "transformCurrentToSource:\n" << transformCurrentToSource << endl;

    transform4f transformSourceToCurrent = transformCurrentToSource.to_inverse();

    // convert the sourceBounds to floats
    boundbox3f sourceBounds;
    sourceBounds.set( rleSource.get_rle_index_spec().outer_bounds().minimum(),
                      rleSource.get_rle_index_spec().outer_bounds().maximum() + vector3( 1 ) );

    // transform the bounds
    boundbox3f transformedBounds;
    transformedBounds = transformSourceToCurrent * sourceBounds;

    // cout << "sourceBounds: " << sourceBounds << endl;
    // cout << "transformedBounds: " << transformedBounds << endl;

    ris.m_abcCoordOrigin = vector3( (int)transformedBounds.minimum().x, (int)transformedBounds.minimum().y,
                                    (int)transformedBounds.minimum().z );
    ris.m_abcCoordSize = size3( (int)transformedBounds.size().xsize(), (int)transformedBounds.size().ysize(),
                                (int)transformedBounds.size().zsize() );

    // Reset the recorded data size to 0, so when we create the named channels they start with a size of zero
    ris.m_dataSize = 0;

    // Get the vertex channel names in the mesh
    vector<frantic::tstring> channelNames;
    rleSource.get_channel_names( channelNames );
    rle_channel_general_accessor x;
    vector<rle_channel_general_accessor> inputChannels;
    vector<rle_channel_general_accessor> outputChannels;

    // Build the corresponding arrays of input and output channels.
    for( unsigned i = 0; i < channelNames.size(); ++i ) {
        // A small hard-coded list of channels we don't want to copy into the level set.  For example, we
        // don't want a channel of normals...
        if( channelNames[i] != _T("Normal") ) {
            inputChannels.push_back( rleSource.get_channel_general_accessor( channelNames[i] ) );
            add_channel( channelNames[i], inputChannels.back().arity(), inputChannels.back().data_type() );
            outputChannels.push_back( get_channel_general_accessor( channelNames[i] ) );
        }
    }

    // If the exterior region code was forced, use the one specified
    // TODO: can we default the exterior to outside? Do Flood and Flood:Spray do that?
    if( rleSource.get_rle_index_spec().m_exteriorRegionCode < 0 ) {
        ris.m_exteriorRegionCode = rleSource.get_rle_index_spec().m_exteriorRegionCode;
    }

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );

    ris.m_runData.clear();
    m_distanceData.clear();

    // find the largest primitive_size in the input channels and create a vector of them of size 8

    size_t largestPrimitiveSize = 1;
    for( size_t i = 0; i < inputChannels.size(); ++i ) {
        if( inputChannels[i].primitive_size() > largestPrimitiveSize ) {
            largestPrimitiveSize = inputChannels[i].primitive_size();
        }
    }

    // Create a temp buffer to hold one intermediate named channel value.
    // and collect a set of pointers into each voxel block
    vector<char> data( 8 * largestPrimitiveSize );
    char* dataPTR[8];

    for( int i = 0; i < 8; ++i ) {
        dataPTR[i] = &data[i * largestPrimitiveSize];
    }

    // rleSource.dump( cout );

    vector3 xyz;
    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            xyz.y = b + ris.m_abcCoordOrigin.y;
            xyz.z = c + ris.m_abcCoordOrigin.z;

            int bcIndex = b + c * ris.m_abcCoordSize.ysize();
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            // The default state is an undefined region of the exterior region code
            bool creatingUndefinedRun = true;
            int creatingRegionCode = ris.m_exteriorRegionCode;
            // Scan along X all the way.
            int xEnd = ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize();
            float levelSetDistance = 0;
            for( xyz.x = ris.m_abcCoordOrigin.x; xyz.x < xEnd; ++xyz.x ) {

                // get the voxel coordinate xyz position in the rleSource vcs
                vector3f xyzSource = transformCurrentToSource * xyz;

                // store the weights and indices for later, also we can use the indices to check for undefined voxels
                float trilerpWeights[8] = { 0, 0, 0, 0, 0, 0, 0, 0 }; // TODO: This should only be of size 3
                boost::int32_t trilerpDataIndices[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
                get_trilerp_indices( rleSource.get_rle_index_spec(), xyzSource, trilerpWeights, trilerpDataIndices );

                bool voxelIsDefined = true;

                // if the dataIndices are negative (and all the same), we are undefined
                if( trilerpDataIndices[0] < 0 && trilerpDataIndices[0] == trilerpDataIndices[1] &&
                    trilerpDataIndices[0] == trilerpDataIndices[2] && trilerpDataIndices[0] == trilerpDataIndices[3] &&
                    trilerpDataIndices[0] == trilerpDataIndices[4] && trilerpDataIndices[0] == trilerpDataIndices[5] &&
                    trilerpDataIndices[0] == trilerpDataIndices[6] && trilerpDataIndices[0] == trilerpDataIndices[7] ) {
                    voxelIsDefined = false;
                }

                if( voxelIsDefined ) {

                    levelSetDistance = 0;
                    for( int i = 0; i < 8; ++i ) {
                        if( trilerpDataIndices[i] > 0 ) {
                            // TODO: trilerpWeights will only have first 3 instantiated.
                            //		 Is this supposed to be trilerpMultipliers?
                            levelSetDistance += trilerpWeights[i] * rleSource[trilerpDataIndices[i]];
                        }
                    }

                    // This is a defined voxel
                    // Flag the start of a new defined run if we switched over from an undefined run
                    if( creatingUndefinedRun ) {
                        ris.m_runData.push_back( run_data( xyz.x, (int)m_distanceData.size() ) );
                        creatingUndefinedRun = false;
                    }
                    // Add the distance value
                    m_distanceData.push_back( levelSetDistance );

                    // Add the values for all the additional named channels by doing an appropriate trilinear
                    // interpolation (using the weights and indices from above)
                    // TODO: This will not do a trilinear interpolation because get_weights only returns 3 deltas.
                    for( unsigned i = 0; i < outputChannels.size(); ++i ) {

                        size_t primitiveSize = inputChannels[i].primitive_size();
                        for( int j = 0; j < 8; ++j ) {
                            if( trilerpDataIndices[j] < 0 ) {
                                memset( &data[j * primitiveSize], 0, primitiveSize );
                            } else {
                                memcpy( &data[j * primitiveSize], inputChannels[i].data( trilerpDataIndices[j] ),
                                        primitiveSize );
                            }
                        }

                        inputChannels[i].get_trilinear_interpolated( trilerpWeights, trilerpDataIndices,
                                                                     outputChannels[i].add_element() );
                    }
                } else {
                    int regionCode = trilerpDataIndices[0]; // -1 means "outside", -2 means "inside"
                    // This is an undefined voxel
                    if( !creatingUndefinedRun ) {
                        ris.m_runData.push_back( run_data( xyz.x, regionCode ) );
                        creatingUndefinedRun = true;
                        creatingRegionCode = regionCode;
                    }
                }
            }

            // If any runs were flagged, finish off the last run.  Otherwise add a zero-sized run
            if( (int)ris.m_runData.size() > ris.m_bcToRunIndex[bcIndex] ) {
                // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to optimize
                // the storage and access.
                if( !( creatingUndefinedRun && creatingRegionCode == ris.m_exteriorRegionCode ) )
                    ris.m_runData.push_back( run_data( xEnd, -1 ) );
            } else {
                ris.m_runData.push_back( run_data( 0, -1 ) );
                ris.m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }
    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = m_distanceData.size();

    // dump(cout);
}

void rle_level_set::switch_rle_index_spec_with_swap( rle_index_spec& ris,
                                                     const frantic::tstring& populatedChannelToCreate ) {
    //	ofstream fout("c:\\temp\\rls.txt", ios::app | ios::out);

    //	fout << "outside distance: " << m_outsideDistance << endl;
    //	fout << "inside distance: " << m_insideDistance << endl;

    // Initialize the level set we're computing from this
    rle_level_set result( m_voxelCoordSystem );
    result.m_insideDistance = m_insideDistance;
    result.m_outsideDistance = m_outsideDistance;
    result.m_interfaceVoxelWidthInside = m_interfaceVoxelWidthInside;
    result.m_interfaceVoxelWidthOutside = m_interfaceVoxelWidthOutside;

    // Swap in the provided rle_index_spec
    result.m_rleIndex.swap( ris );

    // Create all the data channels we need in the result
    result.m_distanceData.resize( result.m_rleIndex.data_size() );
    for( std::map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.begin();
         i != m_namedChannels.end(); ++i ) {
        result.add_channel( i->second.name(), i->second.arity(), i->second.data_type() );
        // Initialize each new channel to all 0 bytes
        std::map<frantic::tstring, rle_channel>::iterator iResult = result.m_namedChannels.find( i->second.name() );
        memset( iResult->second.m_data.begin(), 0, iResult->second.m_data.size() );
    }

    // Create the populated channel.  By initializing all the values in the input to 1, and all the values
    // in the output to 0, we will automatically end up with a result that has 1 wherever there was data defined
    // and 0 whereever a new defined voxel was created.
    if( !populatedChannelToCreate.empty() ) {
        // Add the populated channel to the source level set, setting all its values to 1
        add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( populatedChannelToCreate );
        memset( i->second.m_data.begin(), 1, i->second.m_data.size() );

        // Add the populated channel to the result level set, setting all its values to 0
        result.add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        i = result.m_namedChannels.find( populatedChannelToCreate );
        memset( i->second.m_data.begin(), 0, i->second.m_data.size() );
    }

    // Build the accessors to all the named channels
    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;
    for( std::map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.begin();
         i != m_namedChannels.end(); ++i ) {
        inputAccessors.push_back( get_channel_general_accessor( i->first ) );
        outputAccessors.push_back( result.get_channel_general_accessor( i->first ) );
    }
    // Fill the channel copying data array
    vector<rle_index_spec_channel_copying_data> channelCopyData( inputAccessors.size() + 1 );
    for( size_t i = 0, ie = inputAccessors.size(); i != ie; ++i ) {
        rle_index_spec_channel_copying_data& cd = channelCopyData[i];
        cd.inputData = inputAccessors[i].data( 0 );
        cd.outputData = outputAccessors[i].data( 0 );
        cd.primitiveSize = inputAccessors[i].primitive_size();
    }
    // Also add the m_distanceData channel
    channelCopyData.back().inputData = m_distanceData.empty() ? 0 : reinterpret_cast<const char*>( &m_distanceData[0] );
    channelCopyData.back().outputData =
        result.m_distanceData.empty() ? 0 : reinterpret_cast<char*>( &result.m_distanceData[0] );
    channelCopyData.back().primitiveSize = sizeof( float );
    // Use the rle_index_spec method to copy all the defined data.
    rle_index_spec::copy_data_channels( result.m_rleIndex, m_rleIndex, &channelCopyData[0],
                                        &channelCopyData[0] + channelCopyData.size() );

    //	fout << "exterior region code: " << m_rleIndex.m_exteriorRegionCode << endl;
    //	fout << "exterior distance: " << exteriorDistanceToSet << endl;

    //	fout << "Input abcCoordOrigin: " << m_rleIndex.m_abcCoordOrigin << endl;
    //	fout << "Input abcCoordSize: " << m_rleIndex.m_abcCoordSize << endl;
    //	fout << "B range: [" << bBegin << ", " << bEnd << ")" << endl;
    //	fout << "C range: [" << cBegin << ", " << cEnd << ")" << endl;

    // Now, all the locations that are defined in the output, but weren't defined in the output, don't
    // have a good signed distance value.  This goes through and sets those values.
    for( int c = 0; c < result.m_rleIndex.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < result.m_rleIndex.m_abcCoordSize.ysize(); ++b ) {
            // Get the bcIndex for the input rle index spec
            int bInput = b + result.m_rleIndex.m_abcCoordOrigin.y - m_rleIndex.m_abcCoordOrigin.y;
            int cInput = c + result.m_rleIndex.m_abcCoordOrigin.z - m_rleIndex.m_abcCoordOrigin.z;

            rle_pairwise_run_iterator i( result.m_rleIndex, b, c, m_rleIndex, bInput, cInput ), ie;

            for( ; i != ie; ++i ) {
                // If this sub-interval is defined in result, and not in self, set the level set values.
                if( i.get_first_data_index() >= 0 && i.get_second_data_index() < 0 ) {
                    int dataIndex = i.get_first_data_index();
                    int dataIndexEnd = dataIndex + i.get_xsize();
                    if( i.get_second_data_index() == -1 ) {
                        for( ; dataIndex != dataIndexEnd; ++dataIndex )
                            result.m_distanceData[dataIndex] = m_outsideDistance;
                    } else {
                        for( ; dataIndex != dataIndexEnd; ++dataIndex )
                            result.m_distanceData[dataIndex] = m_insideDistance;
                    }
                }
            }
        }
    }

    // Swap the answer back into this
    result.swap( *this );
}

void rle_level_set::create_data_index_map_channel( const frantic::tstring& dataIndexMapChannelToCreate,
                                                   const rle_index_spec& risMappingTarget ) {
    if( !channels::is_valid_channel_name( dataIndexMapChannelToCreate ) )
        throw std::runtime_error( "rle_level_set.create_data_index_map_channel() - Tried to create channel \"" +
                                  frantic::strings::to_string( dataIndexMapChannelToCreate ) +
                                  "\".  The specified channel name is not valid for a channel name, channel names must "
                                  "start with a letter or "
                                  "underscore, and must contain only letters, numbers, and underscores." );

    add_channel( dataIndexMapChannelToCreate, 1, channels::data_type_int32 );

    rle_channel_accessor<boost::int32_t> dataIndexMapChannelAccessor =
        get_channel_accessor<boost::int32_t>( dataIndexMapChannelToCreate );
    if( dataIndexMapChannelAccessor.size() > 0 )
        m_rleIndex.fill_data_index_map( risMappingTarget, &dataIndexMapChannelAccessor[0] );
}

namespace detail {

/**
 * Computes the three staggered volumes for the faces of the given voxel coordinates.  It returns true if it actually
 * computed the volumes, and false otherwise, in which case a simpler estimate should be used.
 *
 */
bool compute_volume_fractions( boost::int32_t x, boost::int32_t y, boost::int32_t z, float inverseVoxelLength,
                               const rle_index_spec& ris, const float* distanceChannel, float* outCenteredVolume,
                               vector3f* outStaggeredVolumes ) {
    // Initialize the result to 0
    if( outCenteredVolume )
        *outCenteredVolume = 0;
    if( outStaggeredVolumes )
        outStaggeredVolumes->set( 0 );
    // In this case, we estimate the amount of volume inside each of the three face-centered cubes
    // using the formula from Bridson.  First get the 3x3x3 array of data indices for this interpolation.
    boost::int32_t adjacentIndices[3 * 3 * 3];
    boundbox3 adjacentBox( x - 1, x + 1, y - 1, y + 1, z - 1, z + 1 );
    ris.fill_data_index_box_boundary_duplicated( adjacentBox, adjacentIndices );
    // The first diagram is the middle slice of voxel centers where the distances array corresponds to.
    // The * and . are the X and Y faces for which we are estimating their surrounding cube volumes.
    // The second diagram shows the X cube filled with *s.
    // Note that we don't need index 17 to do the estimation, but we calculate it any way to keep the
    // code simpler.
    //
    // 15 | 16  | 17    +  |  +  |  +
    // ___|_____|___    ___|_____|___
    //    |     |        *****   |
    // 12 * 13  | 14    +*****+  |  +
    // ___|__.__|___    _*****___|___
    //    |     |          |     |
    // 9  | 10  | 11    +  |  +  |  +
    // Now convert the indexes into distances for the interpolation.
    // This 3x3x3 grid is positioned at voxel centers.
    float distances[3 * 3 * 3];
    for( int j = 0; j < 3 * 3 * 3; ++j ) {
        int index = adjacentIndices[j];
        if( index >= 0 )
            distances[j] = distanceChannel[index];
        else
            return false;
    }
    //		cout << "Distances (3 XY planes):" << endl;
    //		for( int b = 0; b < 3; ++b ) {
    //			for( int c = 0; c < 3; ++c ) {
    //				for( int a = 0; a < 3; ++a ) {
    //					cout << distances[a + (b + c * 3) * 3] << " ";
    //				}
    //				cout << "   ";
    //			}
    //			cout << endl;
    //		}
    // Now use the distances to compute sub-voxel volume estimates.  In the
    // diagram below, A, B, ..., I are a 2D projection of the sub-voxel centers
    // where we linearly interpolate the distance function to compute the estimates
    //
    // _____|_________|_
    //      |         |
    //   G  |  H   I  |
    // 12   |    13   |
    //   D  |  E   F  |
    // _____|_________|_
    //      |         |
    //   A  |  B   C  |
    // 9    |    10   |
    //
    float aAlpha, bAlpha, cAlpha;
    int aMin, bMin, cMin;
    for( int c = 0; c < 3; ++c ) {
        // If c is 0 or 1, we interpolate from distance samples Z=0 to 1.  If c is 2, it's from Z=1 to 2.
        if( c < 2 ) {
            cAlpha = c * 0.5f + 0.25f;
            cMin = 0;
        } else {
            cAlpha = 0.25f;
            cMin = 1;
        }
        for( int b = 0; b < 3; ++b ) {
            // If b is 0 or 1, we interpolate from distance samples Y=0 to 1.  If c is 2, it's from Y=1 to 2.
            if( b < 2 ) {
                bAlpha = b * 0.5f + 0.25f;
                bMin = 0;
            } else {
                bAlpha = 0.25f;
                bMin = 1;
            }
            for( int a = 0; a < 3; ++a ) {
                bool centered = outCenteredVolume && a > 0 && b > 0 && c > 0,
                     staggeredX = outStaggeredVolumes && a < 2 && b > 0 && c > 0,
                     staggeredY = outStaggeredVolumes && a > 0 && b < 2 && c > 0,
                     staggeredZ = outStaggeredVolumes && a > 0 && b > 0 && c < 2;

                if( centered || staggeredX || staggeredY || staggeredZ ) {
                    if( a < 2 ) {
                        aAlpha = a * 0.5f + 0.25f;
                        aMin = 0;
                    } else {
                        aAlpha = 0.25f;
                        aMin = 1;
                    }

                    int minIndex = aMin + ( bMin + cMin * 3 ) * 3;

                    // Compute the signed distance value at the center of this sub-voxel
                    float dist = ( 1 - cAlpha ) * ( ( 1 - bAlpha ) * ( ( 1 - aAlpha ) * distances[minIndex] +
                                                                       aAlpha * distances[minIndex + 1] ) +
                                                    bAlpha * ( ( 1 - aAlpha ) * distances[minIndex + 3] +
                                                               aAlpha * distances[minIndex + 3 + 1] ) ) +
                                 cAlpha * ( ( 1 - bAlpha ) * ( ( 1 - aAlpha ) * distances[minIndex + 9] +
                                                               aAlpha * distances[minIndex + 9 + 1] ) +
                                            bAlpha * ( ( 1 - aAlpha ) * distances[minIndex + 9 + 3] +
                                                       aAlpha * distances[minIndex + 9 + 3 + 1] ) );

                    //						cout << "X?: " << staggeredX << " Y?: " << staggeredY << " Z?: " << staggeredZ <<
                    //endl; 						cout << "aAlpha: " << aAlpha << ", bAlpha: " << bAlpha << ", cAlpha: " << cAlpha << endl;
                    //						cout << "interpolated distance for sub-volume " << a << "," << b << "," << c << ": "
                    //<< dist << endl;

                    // volume fraction = 1/2 - 1/2 * clamp[-1,1](2 * phi / (0.5 * deltaX))
                    float subVolume = 0.125f * ( 0.5f - 2 * dist * inverseVoxelLength );
                    if( subVolume < 0 )
                        subVolume = 0;
                    else if( subVolume > 0.125f )
                        subVolume = 0.125f;

                    //						cout << "estimated sub-volume: " << subVolume << endl;

                    // Add the subvolume estimate to the relevant volumes
                    if( centered )
                        *outCenteredVolume += subVolume;
                    if( staggeredX )
                        outStaggeredVolumes->x += subVolume;
                    if( staggeredY )
                        outStaggeredVolumes->y += subVolume;
                    if( staggeredZ )
                        outStaggeredVolumes->z += subVolume;
                }
            }
        }
    }
    return true;
}
} // namespace detail

namespace detail {
template <class rle_channel_provider>
struct rle_channel_provider_traits;

template <>
struct rle_channel_provider_traits<rle_level_set> {
    static std::string class_name() { return "rle_level_set"; }
};

template <>
struct rle_channel_provider_traits<frantic::fluids::rle_voxel_field> {
    static std::string class_name() { return "rle_voxel_field"; }
};

template <class rle_channel_provider>
void create_volume_fraction_channels( const rle_level_set& inputLS, rle_channel_provider& outputField,
                                      const frantic::tstring& outputCenteredVolumeChannelName,
                                      const frantic::tstring& outputStaggeredVolumeChannelName ) {
    // Make sure the voxel coordinate systems are the same
    if( !inputLS.get_voxel_coord_system().equals( outputField.get_voxel_coord_system() ) )
        throw runtime_error( "rle_level_set::create_staggered_volume_channel() - The input "
                             "rle_level_set and output " +
                             rle_channel_provider_traits<rle_channel_provider>::class_name() +
                             " provided have differing voxel coordinate systems (" +
                             inputLS.get_voxel_coord_system().str() + " versus " +
                             outputField.get_voxel_coord_system().str() + ")" );

    rle_channel_accessor<float> centeredVolumeChannelAcc;
    rle_channel_accessor<vector3f> staggeredVolumeChannelAcc;

    bool createCentered = false, createStaggered = false;
    // Create the output channels, and get accessors to them
    if( !outputCenteredVolumeChannelName.empty() ) {
        outputField.add_channel( outputCenteredVolumeChannelName, 1, channels::data_type_float32 );
        centeredVolumeChannelAcc = outputField.template get_channel_accessor<float>( outputCenteredVolumeChannelName );
        memset( &centeredVolumeChannelAcc[0], 0xa0, 4 * centeredVolumeChannelAcc.size() );
        createCentered = true;
    }
    if( !outputStaggeredVolumeChannelName.empty() ) {
        outputField.add_channel( outputStaggeredVolumeChannelName, 3, channels::data_type_float32 );
        staggeredVolumeChannelAcc =
            outputField.template get_channel_accessor<vector3f>( outputStaggeredVolumeChannelName );
        createStaggered = true;
    }

    // If neither field is being created, don't do anything
    if( !createCentered && !createStaggered )
        return;

    const rle_index_spec& risDest = outputField.get_rle_index_spec();
    const rle_index_spec& risLS = inputLS.get_rle_index_spec();

    const float voxelLength = inputLS.get_voxel_coord_system().voxel_length();

    float inverseVoxelLength = 1.f / voxelLength;
    // If the voxel center is more than this distance away from the interface, the staggered volumes will be 1.0 or 0.0.
    float distanceThreshold = 1.8f * voxelLength;

    boundbox3 outerBounds = risDest.outer_bounds();
    for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
        for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
            rle_pairwise_run_iterator i( risDest, risDest.y_to_b( y ), risDest.z_to_c( z ), risLS, risLS.y_to_b( y ),
                                         risLS.z_to_c( z ) ),
                ie;
            // Go through all the sub-intervals shared between the destination voxel field and the level set
            for( ; i != ie; ++i ) {
                // Process everything where the destination sub-interval is defined
                if( i.get_first_data_index() >= 0 ) {
                    // If the level set is undefined here, all the volume estimates are set to either 1 or 0
                    if( i.get_second_data_index() < 0 ) {
                        float volumeFraction = ( i.get_second_data_index() < -1 ) ? 1.f : 0.f;
                        if( createCentered ) {
                            for( float *fc = &centeredVolumeChannelAcc[i.get_first_data_index()],
                                       *fcEnd = &centeredVolumeChannelAcc[i.get_first_data_index() + i.get_xsize()];
                                 fc != fcEnd; ++fc )
                                *fc = volumeFraction;
                        }
                        if( createStaggered ) {
                            for( vector3f *vc = &staggeredVolumeChannelAcc[i.get_first_data_index()],
                                          *vcEnd = &staggeredVolumeChannelAcc[i.get_first_data_index() + i.get_xsize()];
                                 vc != vcEnd; ++vc )
                                vc->set( volumeFraction );
                        }
                    }
                    // Otherwise, the level set is defined here, so we have to look at the signed distance values
                    else {
                        float* fc = createCentered ? &centeredVolumeChannelAcc[i.get_first_data_index()] : 0;
                        vector3f* vc = createStaggered ? &staggeredVolumeChannelAcc[i.get_first_data_index()] : 0;
                        int x = i.get_xmin();
                        for( int index = i.get_second_data_index(),
                                 indexEnd = i.get_second_data_index() + i.get_xsize();
                             index != indexEnd; ++index, ++x ) {
                            float signedDistance = inputLS[index];
                            //							cout << "Processing voxel coordinate " << x << "," << y << "," << z << " -
                            //world voxel center " << m_voxelCoordSystem.get_world_voxel_center(vector3(x,y,z)) << endl;
                            // cout << "voxel length: "
                            //<< m_voxelCoordSystem.voxel_length() << endl;
                            if( signedDistance < -distanceThreshold || signedDistance > distanceThreshold ||
                                !detail::compute_volume_fractions(
                                    x, y, z, inverseVoxelLength, inputLS.get_rle_index_spec(), &inputLS[0], fc, vc ) ) {
                                // Otherwise, set them to full or empty
                                if( signedDistance < 0 ) {
                                    if( fc )
                                        *fc = 1;
                                    if( vc )
                                        vc->set( 1 );
                                } else {
                                    if( fc )
                                        *fc = 0;
                                    if( vc )
                                        vc->set( 0 );
                                }
                            }
                            if( fc )
                                ++fc;
                            if( vc )
                                ++vc;
                        }
                    }
                }
            }
        }
    }
}
} // namespace detail

void rle_level_set::create_volume_fraction_channels( rle_level_set& outputLS,
                                                     const frantic::tstring& outputCenteredVolumeChannelName,
                                                     const frantic::tstring& outputStaggeredVolumeChannelName ) const {
    detail::create_volume_fraction_channels( *this, outputLS, outputCenteredVolumeChannelName,
                                             outputStaggeredVolumeChannelName );
}

void rle_level_set::create_volume_fraction_channels( fluids::rle_voxel_field& outputField,
                                                     const frantic::tstring& outputCenteredVolumeChannelName,
                                                     const frantic::tstring& outputStaggeredVolumeChannelName ) const {
    detail::create_volume_fraction_channels( *this, outputField, outputCenteredVolumeChannelName,
                                             outputStaggeredVolumeChannelName );
}

void rle_level_set::duplicate_channel( const frantic::tstring& destChannelName,
                                       const frantic::tstring& sourceChannelName ) {
    std::map<frantic::tstring, rle_channel>::iterator sourceChannel = m_namedChannels.find( sourceChannelName );
    if( sourceChannel == m_namedChannels.end() )
        throw std::runtime_error( "rle_level_set.duplicate_channel() - Tried to duplicate channel \"" +
                                  frantic::strings::to_string( sourceChannelName ) +
                                  "\", but no channel of that name exists." );

    if( !channels::is_valid_channel_name( destChannelName ) )
        throw std::runtime_error( "rle_level_set.duplicate_channel() - Tried to create channel \"" +
                                  frantic::strings::to_string( destChannelName ) +
                                  "\".  The specified channel name is not valid for a channel name, channel names must "
                                  "start with a letter or "
                                  "underscore, and must contain only letters, numbers, and underscores." );

    // Create the destination channel with the appropriate data type
    std::map<frantic::tstring, rle_channel>::iterator destChannel = m_namedChannels.find( destChannelName );
    if( destChannel == m_namedChannels.end() ) {
        m_namedChannels.insert(
            std::make_pair( destChannelName, rle_channel( destChannelName, sourceChannel->second.arity(),
                                                          sourceChannel->second.data_type() ) ) );
        destChannel = m_namedChannels.find( destChannelName );
        // Make sure the vertex array count matches that of the rle_index_spec's data size
        destChannel->second.m_data.resize( m_rleIndex.data_size() * destChannel->second.primitive_size() );
    } else if( destChannel->second.data_type() != sourceChannel->second.data_type() ||
               destChannel->second.arity() != sourceChannel->second.arity() ) {
        destChannel->second.set( destChannelName, sourceChannel->second.arity(), sourceChannel->second.data_type() );
        // Make sure the vertex array count matches that of the rle_index_spec's data size
        destChannel->second.m_data.resize( m_rleIndex.data_size() * destChannel->second.primitive_size() );
    }

    // Copy the memory to the duplicated channel
    memcpy( destChannel->second.m_data.begin(), sourceChannel->second.m_data.begin(),
            m_rleIndex.data_size() * destChannel->second.primitive_size() );
}

void rle_level_set::duplicate_signed_distance_channel( const frantic::tstring& destChannelName ) {
    if( !channels::is_valid_channel_name( destChannelName ) )
        throw std::runtime_error( "rle_level_set.duplicate_signed_distance_channel() - Tried to create channel \"" +
                                  frantic::strings::to_string( destChannelName ) +
                                  "\".  The specified channel name is not valid for a channel name, channel names must "
                                  "start with a letter or "
                                  "underscore, and must contain only letters, numbers, and underscores." );

    // arity and data type of the signed distance
    const std::size_t arity = 1;
    const data_type_t dataType = frantic::channels::data_type_float32;

    // Create the destination channel with the appropriate data type
    std::map<frantic::tstring, rle_channel>::iterator destChannel = m_namedChannels.find( destChannelName );
    if( destChannel == m_namedChannels.end() ) {
        m_namedChannels.insert( std::make_pair( destChannelName, rle_channel( destChannelName, arity, dataType ) ) );
        destChannel = m_namedChannels.find( destChannelName );
        if( destChannel == m_namedChannels.end() )
            throw std::runtime_error( "rle_level_set.duplicate_signed_distance_channel() - Created the channel \"" +
                                      frantic::strings::to_string( destChannelName ) +
                                      "\" but it could not be found afterwards." );
        // Make sure the vertex array count matches that of the rle_index_spec's data size
        destChannel->second.m_data.resize( m_rleIndex.data_size() * destChannel->second.primitive_size() );
    } else if( destChannel->second.data_type() != dataType || destChannel->second.arity() != arity ) {
        destChannel->second.set( destChannelName, arity, dataType );
        // Make sure the vertex array count matches that of the rle_index_spec's data size
        destChannel->second.m_data.resize( m_rleIndex.data_size() * destChannel->second.primitive_size() );
    }

    if( m_distanceData.size() != m_rleIndex.data_size() ) {
        throw std::runtime_error( "rle_level_set.duplicate_signed_distance_channel() - the size of the source and "
                                  "destination channels do not match." );
    }

    if( m_distanceData.size() > 0 ) {
        if( sizeof( m_distanceData[0] ) != destChannel->second.primitive_size() ) {
            throw std::runtime_error( "rle_level_set.duplicate_signed_distance_channel() - the size of the source and "
                                      "destination channel elements do not match." );
        }

        memcpy( destChannel->second.m_data.begin(), &( m_distanceData[0] ),
                m_distanceData.size() * destChannel->second.primitive_size() );
    }
}

void rle_level_set::add_channel( const frantic::tstring& channelName, std::size_t arity, data_type_t dataType ) {
    if( !channels::is_valid_channel_name( channelName ) )
        throw std::runtime_error(
            "rle_level_set.add_channel() - Tried to add channel \"" + frantic::strings::to_string( channelName ) +
            "\", with data type " + frantic::strings::to_string( channel_data_type_str( arity, dataType ) ) +
            ".  "
            "The specified channel name is not valid for a channel name, channel names must start "
            "with a letter or underscore, "
            "and must contain only letters, numbers, and underscores." );

    std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( channelName );
    if( i == m_namedChannels.end() ) {
        m_namedChannels.insert( std::make_pair( channelName, rle_channel( channelName, arity, dataType ) ) );
        i = m_namedChannels.find( channelName );
        // Make sure the voxel data array count matches that of the rle_index_spec's data size
        i->second.m_data.resize( m_rleIndex.data_size() * i->second.primitive_size() );
    } else {
        // If the existing channel doesn't match arity and data type, then we tweak it so it does
        if( i->second.arity() != arity || i->second.data_type() != dataType ) {
            i->second.set( channelName, arity, dataType );
            // Make sure the voxel data array count matches that of the rle_index_spec's data size
            i->second.m_data.resize( m_rleIndex.data_size() * i->second.primitive_size() );
        }
    }
}

void rle_level_set::zero_channel( const frantic::tstring& channelName ) {
    std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( channelName );
    if( i == m_namedChannels.end() ) {
        throw runtime_error( "rle_level_set.zero_channel() - Tried to set channel \"" +
                             frantic::strings::to_string( channelName ) +
                             "\" "
                             "to all zeros, but no such channel exists in the level set." );
    } else {
        // memset the buffer for this channel to all zeros
        size_t size = i->second.m_data.size();
        if( size > 0 )
            memset( i->second.m_data.begin(), 0, size );
    }
}

float rle_level_set::get_channel_max_norm( const frantic::tstring& channelName ) const {
    if( !has_channel( channelName ) )
        throw runtime_error( "rle_level_set.get_channel_max_norm() - Tried to get the maximum L2 norm of channel \"" +
                             frantic::strings::to_string( channelName ) +
                             "\", "
                             "but no such channel exists." );

    const_rle_channel_general_accessor chanAcc = get_channel_general_accessor( channelName );
    if( chanAcc.data_type() != channels::data_type_float32 )
        throw runtime_error( "rle_level_set.get_channel_max_norm() - Tried to get the maximum L2 norm of channel \"" +
                             frantic::strings::to_string( channelName ) +
                             "\", "
                             "but only float32[] data types are supported currently.  Its data type is " +
                             chanAcc.type_str() + "." );

    float maximumNorm = 0;
    for( size_t i = 0, ie = chanAcc.size(); i != ie; ++i ) {
        // Compute the sum of the squares of the components (L2 norm squared)
        float norm = 0;
        const float* data = reinterpret_cast<const float*>( chanAcc.data( i ) );
        for( size_t j = 0, je = chanAcc.arity(); j != je; ++j ) {
            norm += ( *data ) * ( *data );
            ++data;
        }
        if( norm > maximumNorm )
            maximumNorm = norm;
    }

    // Take the square root to return the norm
    return sqrtf( maximumNorm );
}

void rle_level_set::copy_channels( const rle_level_set& inputRLS, const std::vector<frantic::tstring>& channelsToCopy,
                                   const frantic::tstring& populatedChannelToCreate ) {
    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;

    // The channel copying data array
    vector<rle_index_spec_channel_copying_data> channelCopyData;

    // First make sure that all the necesary channels are defined
    for( size_t i = 0, ie = channelsToCopy.size(); i != ie; ++i ) {
        if( !inputRLS.has_channel( channelsToCopy[i] ) )
            throw runtime_error(
                "rle_level_set::copy_channels() - The input rle_level_set doesn't have a channel named \"" +
                frantic::strings::to_string( channelsToCopy[i] ) + "\", which was requested as input" );
        inputAccessors.push_back( inputRLS.get_channel_general_accessor( channelsToCopy[i] ) );

        // If the output already had this data channel, then try to reuse it.
        if( has_channel( channelsToCopy[i] ) ) {
            rle_channel_general_accessor outputAccessor = get_channel_general_accessor( channelsToCopy[i] );
            if( outputAccessor.data_type() == inputAccessors.back().data_type() &&
                outputAccessor.arity() == inputAccessors.back().arity() ) {
                outputAccessors.push_back( outputAccessor );
            } else {
                erase_channel( channelsToCopy[i] );
                add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
                outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
            }
        } else {
            add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
            outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
        }

        // Initialize the output channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = inputAccessors.back().data( 0 );
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = inputAccessors.back().primitive_size();
    }

    // If the caller requested that the 'Populated' channel be made, do so
    vector<char> populatedChannelInput( populatedChannelToCreate.empty() ? 0
                                                                         : inputRLS.get_rle_index_spec().data_size() );
    if( !populatedChannelToCreate.empty() ) {
        if( !populatedChannelInput.empty() ) {
            // Initialize the input 'Populated' data to all ones
            memset( &populatedChannelInput[0], 1, populatedChannelInput.size() );
        }
        add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        outputAccessors.push_back( get_channel_general_accessor( populatedChannelToCreate ) );
        // Initialize the output 'Populated' channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = populatedChannelInput.empty() ? 0 : &populatedChannelInput[0];
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = 1;
    }

    // Use the rle_index_spec method to copy all the defined data.
    if( !channelCopyData.empty() )
        rle_index_spec::copy_data_channels( m_rleIndex, inputRLS.get_rle_index_spec(), &channelCopyData[0],
                                            &channelCopyData[0] + channelCopyData.size() );
}

void rle_level_set::copy_channels( const rle_voxel_field& inputRVF, const std::vector<frantic::tstring>& channelsToCopy,
                                   const frantic::tstring& populatedChannelToCreate ) {
    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;

    // The channel copying data array
    vector<rle_index_spec_channel_copying_data> channelCopyData;

    // First make sure that all the necesary channels are defined
    for( size_t i = 0, ie = channelsToCopy.size(); i != ie; ++i ) {
        if( !inputRVF.has_channel( channelsToCopy[i] ) )
            throw runtime_error(
                "rle_level_set::copy_channels() - The input rle_voxel_field doesn't have a channel named \"" +
                frantic::strings::to_string( channelsToCopy[i] ) + "\", which was requested as input" );
        inputAccessors.push_back( inputRVF.get_channel_general_accessor( channelsToCopy[i] ) );

        // If the output already had this data channel, then try to reuse it.
        if( has_channel( channelsToCopy[i] ) ) {
            rle_channel_general_accessor outputAccessor = get_channel_general_accessor( channelsToCopy[i] );
            if( outputAccessor.data_type() == inputAccessors.back().data_type() &&
                outputAccessor.arity() == inputAccessors.back().arity() ) {
                outputAccessors.push_back( outputAccessor );
            } else {
                erase_channel( channelsToCopy[i] );
                add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
                outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
            }
        } else {
            add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
            outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
        }

        // Initialize the output channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = inputAccessors.back().data( 0 );
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = inputAccessors.back().primitive_size();
    }

    // If the caller requested that the 'Populated' channel be made, do so
    vector<char> populatedChannelInput( populatedChannelToCreate.empty() ? 0
                                                                         : inputRVF.get_rle_index_spec().data_size() );
    if( !populatedChannelToCreate.empty() ) {
        if( !populatedChannelInput.empty() ) {
            // Initialize the input 'Populated' data to all ones
            memset( &populatedChannelInput[0], 1, populatedChannelInput.size() );
        }
        add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        outputAccessors.push_back( get_channel_general_accessor( populatedChannelToCreate ) );
        // Initialize the output 'Populated' channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = populatedChannelInput.empty() ? 0 : &populatedChannelInput[0];
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = 1;
    }

    // Use the rle_index_spec method to copy all the defined data.
    if( !channelCopyData.empty() )
        rle_index_spec::copy_data_channels( m_rleIndex, inputRVF.get_rle_index_spec(), &channelCopyData[0],
                                            &channelCopyData[0] + channelCopyData.size() );
}

class rle_scanline_level_set_advect_rk3_linear {
    const boundbox3& m_bounds;
    const size3& m_boxSize;
    const rle_index_spec &m_velRIS, m_srcRIS;
    const voxel_coord_system &m_srcVcs, m_velVcs; //, m_resultVcs;
    vector3 m_centerOffset;
    int m_cellBoxIndex;

    const_rle_channel_accessor<vector3f>& m_velocityAcc;
    const_rle_channel_accessor<boost::uint8_t>& m_populatedChannelAcc;
    const std::vector<float>& m_srcDistanceData;
    std::vector<float>& m_resultDistanceData;

    float m_dt;

    rle_scanline_level_set_advect_rk3_linear&
    operator=( const rle_scanline_level_set_advect_rk3_linear& ); // not implemented

  public:
    rle_scanline_level_set_advect_rk3_linear( const boundbox3& bounds, const size3& boxSize,
                                              const rle_index_spec& velRIS, const rle_index_spec& srcRIS,
                                              const voxel_coord_system& srcVcs, const voxel_coord_system& velVcs,
                                              const_rle_channel_accessor<vector3f>& velocityAcc,
                                              const_rle_channel_accessor<boost::uint8_t>& populatedChannelAcc,
                                              const std::vector<float>& srcDistanceData,
                                              std::vector<float>& resultDistanceData, float dt )
        : m_bounds( bounds )
        , m_boxSize( boxSize )
        , m_velRIS( velRIS )
        , m_srcRIS( srcRIS )
        , m_srcVcs( srcVcs )
        , m_velVcs( velVcs )
        , m_velocityAcc( velocityAcc )
        , m_populatedChannelAcc( populatedChannelAcc )
        , m_srcDistanceData( srcDistanceData )
        , m_resultDistanceData( resultDistanceData )
        , m_dt( dt ) {
        m_centerOffset = vector3( boxSize.xsize() / 2, boxSize.ysize() / 2, boxSize.zsize() / 2 );
        m_cellBoxIndex = m_centerOffset.x + m_centerOffset.y * m_boxSize.xsize() +
                         m_centerOffset.z * m_boxSize.xsize() * m_boxSize.ysize();
    }

    static int get_required_dilation( float maxVoxelMotion ) {
        return 3 + static_cast<int>( max<float>( 0, ceilf( maxVoxelMotion - 0.49f ) ) );
        // return 4;
    }

    static size3 get_index_box_size( float maxVoxelMotion ) {
        const int maxWindowSize = 7;
        const int voxelMotionWindowSize = 3 + 2 * static_cast<int>( max<float>( 0, ceilf( maxVoxelMotion - 0.49f ) ) );
        return size3( std::min<int>( maxWindowSize, voxelMotionWindowSize ) );
    }

    void operator()( const tbb::blocked_range2d<size_t>& r ) const {
        const std::vector<boost::int32_t>& bcToRunIndex = m_srcRIS.get_bc_to_run_index_vector();
        const std::vector<run_data>& runIndexData = m_srcRIS.get_run_index_data_vector();
        const vector3& boundsMin = m_bounds.minimum();

        int ysize = m_bounds.ysize(); //, zsize = m_bounds.zsize();

        for( size_t c = r.rows().begin(); c != r.rows().end(); ++c ) {
            int z = static_cast<int>( c ) + boundsMin.z;
            size_t cOffset = c * ysize;
            for( size_t b = r.cols().begin(); b != r.cols().end(); ++b ) {
                int y = static_cast<int>( b ) + boundsMin.y;
                size_t bcIndex = b + cOffset;

                int runRangeStart = bcToRunIndex[bcIndex];
                int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

                // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the
                // runs with the iterator
                if( runRangeStart != runRangeEnd ||
                    runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                    // construct a block iterator for this scanline

                    // for velocity advection we need one voxel in the negative direction, and two in hte positive
                    // to make sure we have all the data indices that we will need to process the lookup. This means we
                    // will end up with more overall indices, but should be faster than carefully choosing and finding
                    // only the ones we need int xmin=boundsMin.x-centerOffset.x, xmax = boundsMin.x+centerOffset.x;

                    // construct our iterator box
                    boundbox3 box(
                        vector3( boundsMin.x - m_centerOffset.x, y - m_centerOffset.y, z - m_centerOffset.z ),
                        m_boxSize );

                    for( levelset::rle_defined_box_iterator boxIter( m_srcRIS, box );
                         !boxIter.is_xplane_finished( m_centerOffset.x ); ++boxIter ) {

                        const boost::int32_t* const dataIndices = boxIter.get_indices();

                        /*		for (int i=0; i!= boxSize.volume(); ++i ){
                              logging::error << dataIndices[i] << endl;
                            }*/
                        boost::int32_t index = dataIndices[m_cellBoxIndex];

                        // if the center cell is not actually defined then we can skip
                        if( index < 0 ) {
                            // logging::error << "undef: "  << currentMin + centerOffset << std::endl;
                            continue;
                        }

                        if( m_populatedChannelAcc[index] ) {
                            const boundbox3& currentBox = boxIter.current_box();
                            const vector3& currentMin = currentBox.minimum();

                            const vector3 cellVoxelCoord = currentMin + m_centerOffset;
                            const vector3f lookup = m_srcVcs.get_world_voxel_center( cellVoxelCoord );
                            vector3f finalPosition;
                            fluids::trace_runge_kutta_3( m_velVcs, m_velRIS, m_velocityAcc, lookup, -m_dt,
                                                         finalPosition );

                            const vector3f finalPositionVoxelCoordFloat = m_srcVcs.get_voxel_coord( finalPosition );
                            const vector3 finalPositionVoxelCoord =
                                vector3::from_floor( finalPositionVoxelCoordFloat - vector3f( 0.5f ) );
                            const boundbox3 definedBox( currentBox.minimum(), currentBox.maximum() - vector3( 1 ) );

                            if( definedBox.contains( finalPositionVoxelCoord ) ) {
                                // TODO don't silently ignore the returned error code !
                                trilerp_float_from_indices( &m_srcDistanceData[0], finalPositionVoxelCoordFloat,
                                                            dataIndices, currentBox, &m_resultDistanceData[index] );
                            } else {
                                // TODO don't silently ignore the returned error code !
                                trilerp_float( m_srcRIS, &m_srcDistanceData[0], finalPositionVoxelCoordFloat,
                                               &m_resultDistanceData[index] );
                            }
                        }
                    }
                }
            }
        }
    }
};

class rle_scanline_level_set_advect_rk3_weno3 {
    const boundbox3& m_bounds;
    const size3& m_boxSize;
    const rle_index_spec &m_velRIS, m_srcRIS;
    const voxel_coord_system &m_srcVcs, m_velVcs; //, m_resultVcs;
    vector3 m_centerOffset;
    int m_cellBoxIndex;

    const_rle_channel_accessor<vector3f>& m_velocityAcc;
    const_rle_channel_accessor<boost::uint8_t>& m_populatedChannelAcc;
    const std::vector<float>& m_srcDistanceData;
    std::vector<float>& m_resultDistanceData;

    float m_dt;

    rle_scanline_level_set_advect_rk3_weno3&
    operator=( const rle_scanline_level_set_advect_rk3_weno3& ); // not implemented

  public:
    rle_scanline_level_set_advect_rk3_weno3( const boundbox3& bounds, const size3& boxSize,
                                             const rle_index_spec& velRIS, const rle_index_spec& srcRIS,
                                             const voxel_coord_system& srcVcs, const voxel_coord_system& velVcs,
                                             const_rle_channel_accessor<vector3f>& velocityAcc,
                                             const_rle_channel_accessor<boost::uint8_t>& populatedChannelAcc,
                                             const std::vector<float>& srcDistanceData,
                                             std::vector<float>& resultDistanceData, float dt )
        : m_bounds( bounds )
        , m_boxSize( boxSize )
        , m_velRIS( velRIS )
        , m_srcRIS( srcRIS )
        , m_srcVcs( srcVcs )
        , m_velVcs( velVcs )
        , m_velocityAcc( velocityAcc )
        , m_populatedChannelAcc( populatedChannelAcc )
        , m_srcDistanceData( srcDistanceData )
        , m_resultDistanceData( resultDistanceData )
        , m_dt( dt )

    {
        m_centerOffset = vector3( boxSize.xsize() / 2, boxSize.ysize() / 2, boxSize.zsize() / 2 );
        m_cellBoxIndex = m_centerOffset.x + m_centerOffset.y * m_boxSize.xsize() +
                         m_centerOffset.z * m_boxSize.xsize() * m_boxSize.ysize();
    }

    static int get_required_dilation( float maxVoxelMotion ) {
        return 3 + static_cast<int>( max<float>( 0, ceilf( maxVoxelMotion - 0.49f ) ) );
    }

    static size3 get_index_box_size( float maxVoxelMotion ) {
        const int maxWindowSize = 7;
        const int voxelMotionWindowSize = 5 + 2 * static_cast<int>( max<float>( 0, ceilf( maxVoxelMotion - 0.49f ) ) );
        return size3( std::min<int>( maxWindowSize, voxelMotionWindowSize ) );
    }

    void operator()( const tbb::blocked_range2d<size_t>& r ) const {

        const std::vector<boost::int32_t>& bcToRunIndex = m_srcRIS.get_bc_to_run_index_vector();
        const std::vector<run_data>& runIndexData = m_srcRIS.get_run_index_data_vector();
        const vector3& boundsMin = m_bounds.minimum();

        int ysize = m_bounds.ysize(); //, zsize = m_bounds.zsize();

        for( size_t c = r.rows().begin(); c != r.rows().end(); ++c ) {
            int z = static_cast<int>( c ) + boundsMin.z;
            size_t cOffset = c * ysize;
            for( size_t b = r.cols().begin(); b != r.cols().end(); ++b ) {
                int y = static_cast<int>( b ) + boundsMin.y;
                size_t bcIndex = b + cOffset;

                int runRangeStart = bcToRunIndex[bcIndex];
                int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

                // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the
                // runs with the iterator
                if( runRangeStart != runRangeEnd ||
                    runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                    // construct a block iterator for this scanline

                    // for velocity advection we need one voxel in the negative direction, and two in hte positive
                    // to make sure we have all the data indices that we will need to process the lookup. This means we
                    // will end up with more overall indices, but should be faster than carefully choosing and finding
                    // only the ones we need int xmin=boundsMin.x-centerOffset.x, xmax = boundsMin.x+centerOffset.x;

                    // construct our iterator box
                    boundbox3 box(
                        vector3( boundsMin.x - m_centerOffset.x, y - m_centerOffset.y, z - m_centerOffset.z ),
                        m_boxSize );

                    for( levelset::rle_defined_box_iterator boxIter( m_srcRIS, box );
                         !boxIter.is_xplane_finished( m_centerOffset.x ); ++boxIter ) {

                        const boost::int32_t* const dataIndices = boxIter.get_indices();

                        /*		for (int i=0; i!= boxSize.volume(); ++i ){
                              logging::error << dataIndices[i] << endl;
                            }*/
                        boost::int32_t index = dataIndices[m_cellBoxIndex];

                        // if the center cell is not actually defined then we can skip
                        if( index < 0 ) {
                            // logging::error << "undef: "  << currentMin + centerOffset << std::endl;
                            continue;
                        }

                        if( m_populatedChannelAcc[index] ) {
                            const boundbox3& currentBox = boxIter.current_box();
                            const vector3& currentMin = currentBox.minimum();

                            // logging::error << "Current Index: " << index << endl;
                            // logging::error << "Current Box : " << currentBox << endl;

                            vector3 cellVoxelCoord = currentMin + m_centerOffset;

                            //		if( i.get_coord().z == bounds.minimum().z)
                            //			logging::set_logging_level(5);

                            const vector3f lookup = m_srcVcs.get_world_voxel_center( cellVoxelCoord );
                            vector3f finalPosition;
                            fluids::trace_runge_kutta_3( m_velVcs, m_velRIS, m_velocityAcc, lookup, -m_dt,
                                                         finalPosition );

                            const vector3f finalPositionVoxelCoordFloat = m_srcVcs.get_voxel_coord( finalPosition );
                            // This definedBox is used to test whether all of
                            // the indices required by weno3_signed_distance_lookup
                            // are inside the advection window.
                            // TODO: simplify this test
                            const vector3 finalPositionVoxelCoord =
                                vector3::from_floor( finalPositionVoxelCoordFloat - vector3f( 0.5f ) );
                            const boundbox3 definedBox( currentBox.minimum() + vector3( 1 ),
                                                        currentBox.maximum() - vector3( 2 ) );

                            if( definedBox.contains( finalPositionVoxelCoord ) ) {
                                m_resultDistanceData[index] =
                                    weno3_signed_distance_lookup( dataIndices, m_boxSize, currentMin, m_srcDistanceData,
                                                                  finalPositionVoxelCoordFloat );
                            } else { // outside advection window; need rle lookup
                                boost::int32_t finalDataIndices[64];
                                const boundbox3 finalVoxelExtents( finalPositionVoxelCoord - vector3( 1 ),
                                                                   finalPositionVoxelCoord + vector3( 2 ) );
                                m_srcRIS.fill_data_index_box( finalVoxelExtents, finalDataIndices );
                                m_resultDistanceData[index] = weno3_signed_distance_lookup(
                                    finalDataIndices, size3( 4 ), finalVoxelExtents.minimum(), m_srcDistanceData,
                                    finalPositionVoxelCoordFloat );
                            }
                        }
                    }
                }
            }
        }
    }
};

// Helper functions which are called by semi_lagrangian_advect_staggered()
namespace detail {

/**
 *  Call f(index) for each index that is on the interior surface of the wall.
 */
template <class UnaryFunction>
void for_each_index_on_box_interior_surface( const rle_index_spec& ris, const boundbox3& wall, UnaryFunction f ) {

    const boundbox3& outerBounds = ris.outer_bounds();
    for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
        if( z < wall.minimum().z || z > wall.maximum().z )
            continue;
        for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
            if( y < wall.minimum().y || y > wall.maximum().y )
                continue;

            rle_run_iterator runIter( ris, ris.y_to_b( y ), ris.z_to_c( z ) ), runIterEnd;
            for( ; runIter != runIterEnd; ++runIter ) {
                const boost::int32_t dataIndex = runIter.get_data_index();
                if( dataIndex >= 0 ) {
                    const int xMin = runIter.get_xmin();
                    const int xMax = runIter.get_xmax();
                    if( y == wall.minimum().y || y == wall.maximum().y || z == wall.minimum().z ||
                        z == wall.maximum().z ) {
                        const int xBegin = std::max<int>( xMin, wall.minimum().x );
                        const int xEnd = 1 + std::min<int>( xMax, wall.maximum().x );
                        const int iBegin = xBegin - xMin;
                        const int iEnd = xEnd - xMin;
                        for( int i = iBegin; i < iEnd; ++i ) {
                            f( dataIndex + i );
                        }
                    } else {
                        if( xMin <= wall.minimum().x && wall.minimum().x <= xMax ) {
                            const int i = wall.minimum().x - xMin;
                            f( dataIndex + i );
                        }
                        if( wall.minimum().x != wall.maximum().x && xMin <= wall.maximum().x &&
                            wall.maximum().x <= xMax ) {
                            const int i = wall.maximum().x - xMin;
                            f( dataIndex + i );
                        }
                    }
                }
            }
        }
    }
}

/**
 *  Call f(index) for each index that is outside of the wall.
 */
template <class UnaryFunction>
void for_each_index_outside_box( const rle_index_spec& ris, const boundbox3& wall, UnaryFunction f ) {

    const boundbox3& outerBounds = ris.outer_bounds();
    for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
        for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
            rle_run_iterator runIter( ris, ris.y_to_b( y ), ris.z_to_c( z ) ), runIterEnd;
            for( ; runIter != runIterEnd; ++runIter ) {
                const boost::int32_t dataIndex = runIter.get_data_index();
                if( dataIndex >= 0 ) {
                    const int xMin = runIter.get_xmin();
                    const int xMax = runIter.get_xmax();
                    if( y < wall.minimum().y || y > wall.maximum().y || z < wall.minimum().z || z > wall.maximum().z ) {
                        for( int i = 0; i < runIter.get_xsize(); ++i ) {
                            f( dataIndex + i );
                        }
                    } else {
                        if( xMin < wall.minimum().x ) {
                            int xEnd;
                            if( xMax < wall.minimum().x )
                                xEnd = xMax + 1;
                            else
                                xEnd = wall.minimum().x;
                            const int iEnd = xEnd - xMin;
                            for( int i = 0; i < iEnd; ++i ) {
                                f( dataIndex + i );
                            }
                        }
                        if( xMax > wall.maximum().x ) {
                            int xBegin;
                            if( xMin <= wall.maximum().x )
                                xBegin = wall.maximum().x + 1;
                            else
                                xBegin = xMin;
                            const int iBegin = xBegin - xMin;
                            const int iEnd = runIter.get_xsize();
                            for( int i = iBegin; i < iEnd; ++i ) {
                                f( dataIndex + i );
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 *  Set the populated channel to ignore for the specified index.
 */
struct set_populated_to_ignore {
  private:
    set_populated_to_ignore& operator=( set_populated_to_ignore& ); // not implemented
  public:
    rle_channel_accessor<boost::uint8_t>& m_populated;
    set_populated_to_ignore( rle_channel_accessor<boost::uint8_t>& populated )
        : m_populated( populated ) {}
    void operator()( boost::int32_t index ) { m_populated[index] = 2; }
};

/**
 *  Start extrapolation from populated indices into their unpopulated
 * neighbours.
 *
 *  If the index is populated(1), then for each unpopulated(0) neighbour:
 * set the neighbour's populated flag to 3, and append the neighbour's
 * index to the index vector.
 */
struct add_unpopulated_neighbours_to_near {
  private:
    const ris_adjacency& m_adj;
    rle_channel_accessor<boost::uint8_t>& m_populated;
    std::vector<boost::int32_t>& m_indexVector;
    add_unpopulated_neighbours_to_near& operator=( add_unpopulated_neighbours_to_near& ); // not implemented
  public:
    add_unpopulated_neighbours_to_near( const ris_adjacency& adj, rle_channel_accessor<boost::uint8_t>& populated,
                                        std::vector<boost::int32_t>& indexVector )
        : m_adj( adj )
        , m_populated( populated )
        , m_indexVector( indexVector ) {}
    void operator()( boost::int32_t centerIndex ) {
        if( m_populated[centerIndex] == 1 ) {
            for( int face = 0; face < 6; ++face ) {
                const boost::int32_t adjIndex = m_adj[centerIndex][face];
                if( adjIndex >= 0 && m_populated[adjIndex] == 0 ) {
                    m_populated[adjIndex] = 3;
                    m_indexVector.push_back( adjIndex );
                }
            }
        }
    }
};

/**
 *  Continue extrapolation for one step.  Update the index, and add its
 * unpopulated neighbours.
 *
 *  Note that this function does not update the populated flag for the
 * index.
 *
 *  Update the signed distance of the specified index to the average
 * signed distance of its populated neighbours.  Also, for each
 * unpopulated(0) neighbour, add the unpopulated neighbour's index
 * to the nextWavefront vector, and set its populated flag to 3.
 */
void update_center_and_add_unpopulated_neighbours_to_near( const ris_adjacency& adj, float* phi,
                                                           rle_channel_accessor<boost::uint8_t>& populatedAcc,
                                                           boost::int32_t centerIndex,
                                                           std::vector<boost::int32_t>& nextWavefront ) {
    float phiAccum = 0;
    float weight = 0;
    for( int face = 0; face < 6; ++face ) {
        const boost::int32_t adjIndex = adj[centerIndex][face];
        if( adjIndex >= 0 ) {
            if( populatedAcc[adjIndex] == 1 ) {
                phiAccum += phi[adjIndex];
                weight += 1.f;
            } else if( populatedAcc[adjIndex] == 0 ) {
                nextWavefront.push_back( adjIndex );
                populatedAcc[adjIndex] = 3;
            }
        }
    }
    phi[centerIndex] = phiAccum / weight;
}

/**
 *  Extrapolate the signed distance in the level set away from populated
 * voxels on the inside surface of the bounding box into neighbouring
 * unpopulated voxels.
 *
 *  This is intended for use when the voxels inside the bounding box are
 * all populated, and the voxels outside the bounding box are all
 * unpopulated, but this is not enforced.
 *
 * @note this will allocate ris_adjacency data for the level set.
 *
 * @param[in,out] ls the level set which holds the signed distance to
 *		extrapolate.
 * @param box extrapolation will start from populated voxels that are
 *		on the inside surface of this boundbox.
 * @param populatedChannelName the name of a uint8 populated channel
 *		in the level set.  Voxels with a populated value 1 are populated
 *		and may be used as an extrapolation source, while voxels with
 *		a populated value of 0 may be filled by the extrapolation.
 *
 */
void extrapolate_signed_distance_away_from_populated_box_interior_surface(
    rle_level_set& ls, const frantic::graphics::boundbox3& box, const frantic::tstring& populatedChannelName ) {
    typedef boost::uint8_t populated_t;

    if( !ls.has_channel( populatedChannelName ) ) {
        throw std::runtime_error(
            "extrapolate_signed_distance_away_from_populated_box_interior_surface Error: the level "
            "set is missing the specified \'" +
            frantic::strings::to_string( populatedChannelName ) + "\' populated channel." );
    } else {
        rle_channel_general_accessor tempAcc = ls.get_channel_general_accessor( populatedChannelName );
        typedef frantic::channels::channel_data_type_traits<populated_t> required;
        if( tempAcc.arity() != required::arity() || tempAcc.data_type() != required::data_type() ) {
            throw std::runtime_error( "extrapolate_signed_distance_away_from_populated_box_interior_surface Error: the "
                                      "specified populated channel \'" +
                                      frantic::strings::to_string( populatedChannelName ) + "\' must be of type " +
                                      frantic::strings::to_string( required::type_str() ) +
                                      ", but instead it is of type " +
                                      frantic::strings::to_string( tempAcc.type_str() ) + "." );
        }
    }

    if( ls.size() == 0 )
        return;

    const rle_index_spec& ris = ls.get_rle_index_spec();
    const ris_adjacency& adj = ris.get_cached_adjacency();
    rle_channel_accessor<populated_t> populatedAcc =
        ls.get_channel_accessor<populated_t>( _T("ReinitializationPopulated") );

    std::vector<boost::int32_t> currentWavefront;
    std::vector<boost::int32_t> nextWavefront;

    detail::for_each_index_on_box_interior_surface(
        ris, box, detail::add_unpopulated_neighbours_to_near( adj, populatedAcc, currentWavefront ) );
    while( !currentWavefront.empty() ) {
        for( std::vector<boost::int32_t>::iterator i( currentWavefront.begin() ); i != currentWavefront.end(); ++i ) {
            detail::update_center_and_add_unpopulated_neighbours_to_near( adj, &ls[0], populatedAcc, *i,
                                                                          nextWavefront );
        }
        // update the populated flag of the current wavefront
        // this is done after updating the signed distance of the current
        // wavefront so we don't update using inputs from the same
        // wavefront
        for( std::vector<boost::int32_t>::iterator i( currentWavefront.begin() ); i != currentWavefront.end(); ++i ) {
            populatedAcc[*i] = 1;
        }
        std::swap( currentWavefront, nextWavefront );
        nextWavefront.clear();
    }
}

float get_maximum_voxel_motion_from_staggered_velocity( const rle_voxel_field& field,
                                                        const frantic::tstring& staggeredVelocityChannelName,
                                                        const float dt ) {
    const_rle_channel_accessor<vector3f> velAcc = field.get_channel_accessor<vector3f>( staggeredVelocityChannelName );

    float maxVelocityMagnitude = 0;
    for( size_t i = 0; i < velAcc.size(); ++i ) {
        for( int axis = 0; axis < 3; ++axis ) {
            maxVelocityMagnitude = std::max<float>( maxVelocityMagnitude, fabsf( velAcc[i][axis] ) );
        }
    }

    return maxVelocityMagnitude * fabsf( dt ) / field.get_voxel_coord_system().voxel_length();
}

}; // namespace detail

void rle_level_set::semi_lagrangian_advect_staggered( const boundbox3& voxelBounds,
                                                      const frantic::fluids::rle_voxel_field& velocityField,
                                                      const frantic::tstring& staggeredVelocityChannelName, float dt ) {
    tbb::task_scheduler_init taskSchedulerInit;

    typedef rle_scanline_level_set_advect_rk3_weno3 advector_t;

    if( size() == 0 )
        return;

    const float maxVoxelMotion =
        detail::get_maximum_voxel_motion_from_staggered_velocity( velocityField, staggeredVelocityChannelName, dt );
    const int dilationWidth = advector_t::get_required_dilation( maxVoxelMotion );

    // first lets extrapolate the distance values once and then use that to complete the weno lookups.
    rle_level_set tempRLS( m_voxelCoordSystem, m_rleIndex, m_distanceData, m_insideDistance, m_outsideDistance );
    tempRLS.dilate_defined_voxels( dilationWidth, _T("OriginalPopulatedChannel") );
    tempRLS.duplicate_channel( _T("ReinitializationPopulated"), _T("OriginalPopulatedChannel") );

    const rle_index_spec& tempRIS = tempRLS.get_rle_index_spec();
    // Reinitialize without reading or modifying voxels that are outside of the
    // voxelBounds.
    { // scope for tempReinitializationPopulatedAcc
        rle_channel_accessor<boost::uint8_t> tempReinitializationPopulatedAcc =
            tempRLS.get_channel_accessor<boost::uint8_t>( _T("ReinitializationPopulated") );
        detail::for_each_index_outside_box( tempRIS, voxelBounds,
                                            detail::set_populated_to_ignore( tempReinitializationPopulatedAcc ) );
    }
    tempRLS.reinitialize_signed_distance_from_populated( _T("ReinitializationPopulated") );

    // Assume that, after reinitialization, the populated channel is
    // 0 or 1 everywhere: ignore voxels (2) have been set to unpopulated (0).
    // This is the current behaviour of
    // reinitialize_signed_distance_from_populated, but its specification
    // was recently changed.

    // Now extrapolate away from the populated voxels inside the voxelBounds
    // into unpopulated voxels outside of the voxelBounds.
    // This may be unnecessary or even ill-advised -- we may be better off
    // filling the outside voxels using a constant +halfVoxelLength or similar.
    detail::extrapolate_signed_distance_away_from_populated_box_interior_surface( tempRLS, voxelBounds,
                                                                                  _T("ReinitializationPopulated") );

    const_rle_channel_accessor<boost::uint8_t> populatedChannelAccessor =
        tempRLS.get_channel_accessor<boost::uint8_t>( _T("OriginalPopulatedChannel") );

    std::vector<float> newDistanceData( tempRLS.size() );

    const voxel_coord_system& velVCS = velocityField.get_voxel_coord_system();
    const rle_index_spec& ris = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velAcc =
        velocityField.get_channel_accessor<vector3f>( staggeredVelocityChannelName );

    const size3 boxSize( advector_t::get_index_box_size( maxVoxelMotion ) );

    // create a data index channel for copying the advected values back into the current level set
    create_data_index_map_channel( _T("AdvectedLevelSetLookup"), tempRIS );
    const_rle_channel_accessor<boost::int32_t> indexAccessor =
        get_channel_accessor<boost::int32_t>( _T("AdvectedLevelSetLookup") );

    boundbox3 bounds = tempRIS.outer_bounds();
    boundbox3 originalBounds = get_rle_index_spec().outer_bounds();

    advector_t advector( bounds, boxSize, ris, tempRIS, m_voxelCoordSystem, velVCS, velAcc, populatedChannelAccessor,
                         tempRLS.get_distance_data(), newDistanceData, dt );

    // tbb::blocked_range2d<size_t> range( 0, bounds.zsize(), 0, bounds.ysize() );
    tbb::blocked_range2d<size_t> range(
        originalBounds.zminimum() - bounds.zminimum(), originalBounds.zmaximum() - bounds.zminimum() + 1,
        originalBounds.yminimum() - bounds.yminimum(), originalBounds.ymaximum() - bounds.yminimum() + 1 );

    // logging::error << "threading block range: " << bounds.ysize() << ", " << bounds.zsize() << endl;
    tbb::parallel_for( range, advector, tbb::auto_partitioner() );

    // could probably paralellize this update loop as well
    for( std::size_t i = 0; i != m_distanceData.size(); ++i ) {
        m_distanceData[i] = newDistanceData[indexAccessor[i]];
    }

    // clean up the lookup channel
    erase_channel( _T("AdvectedLevelSetLookup") );
    ////m_distanceData.swap(newDistanceData);

    // std::unique_ptr<rle_level_set> debugRLS( new rle_level_set( tempRLS ) );
    // return debugRLS.release();
}

void rle_level_set::semi_lagrangian_advect_staggered( const frantic::fluids::rle_voxel_field& velocityField,
                                                      const frantic::tstring& staggeredVelocityChannelName, float dt ) {
    typedef rle_scanline_level_set_advect_rk3_weno3 advector_t;

    const float maxVoxelMotion =
        detail::get_maximum_voxel_motion_from_staggered_velocity( velocityField, staggeredVelocityChannelName, dt );
    if( maxVoxelMotion > 12.f ) {
        throw std::runtime_error(
            "rle_level_set::semi_lagrangian_advect_staggered Error: the maximum velocity in the velocity field (" +
            boost::lexical_cast<std::string>( maxVoxelMotion ) + " voxel lengths in channel \'" +
            frantic::strings::to_string( staggeredVelocityChannelName ) +
            "\') exceeds the maximum limit (12 voxel lengths)." );
    }
    const int indexBoxLength = advector_t::get_index_box_size( maxVoxelMotion ).xsize();

    // first lets extrapolate the distance values once and then use that to complete the weno lookups.
    rle_level_set tempRLS( m_voxelCoordSystem, m_rleIndex, m_distanceData, m_insideDistance, m_outsideDistance );
    tempRLS.dilate_defined_voxels( indexBoxLength / 2 + 1, _T("OriginalPopulatedChannel") );
    tempRLS.duplicate_channel( _T("ReinitializationPopulated"), _T("OriginalPopulatedChannel") );

    tempRLS.reinitialize_signed_distance_from_populated( _T("ReinitializationPopulated") );

    // tempRLS.extrapolate_channels( "ReinitializationPopulated" );

    const_rle_channel_accessor<boost::uint8_t> populatedChannelAccessor =
        tempRLS.get_channel_accessor<boost::uint8_t>( _T("OriginalPopulatedChannel") );

    std::vector<float> newDistanceData( tempRLS.size() );

    const voxel_coord_system& velVCS = velocityField.get_voxel_coord_system();
    const rle_index_spec& ris = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velAcc =
        velocityField.get_channel_accessor<vector3f>( staggeredVelocityChannelName );

    size3 boxSize( indexBoxLength );

    const rle_index_spec& tempRIS = tempRLS.get_rle_index_spec();

    // create a data index channel for copying the advected values back into the current level set
    create_data_index_map_channel( _T("AdvectedLevelSetLookup"), tempRIS );

    boundbox3 bounds = tempRIS.outer_bounds();

    advector_t advector( bounds, boxSize, ris, tempRIS, m_voxelCoordSystem, velVCS, velAcc, populatedChannelAccessor,
                         tempRLS.get_distance_data(), newDistanceData, dt );

    tbb::blocked_range2d<size_t> range( 0, bounds.zsize(), 0, bounds.ysize() );

    // logging::error << "threading block range: " << bounds.ysize() << ", " << bounds.zsize() << endl;
    tbb::parallel_for( range, advector, tbb::auto_partitioner() );

    const_rle_channel_accessor<boost::int32_t> indexAccessor =
        get_channel_accessor<boost::int32_t>( _T("AdvectedLevelSetLookup") );

    // could probably paralellize this update loop as well
    for( std::size_t i = 0; i != m_distanceData.size(); ++i ) {
        m_distanceData[i] = newDistanceData[indexAccessor[i]];
    }

    // clean up the lookup channel
    erase_channel( _T("AdvectedLevelSetLookup") );
    ////m_distanceData.swap(newDistanceData);
}

void rle_level_set::serial_semi_lagrangian_advect_staggered( const frantic::fluids::rle_voxel_field& velocityField,
                                                             const frantic::tstring& staggeredVelocityChannelName,
                                                             float dt ) {

    // first lets extrapolate the distance values once and then use that to complete the weno lookups.
    rle_level_set tempRLS( m_voxelCoordSystem, m_rleIndex, m_distanceData, m_insideDistance, m_outsideDistance );
    tempRLS.dilate_defined_voxels( 3, _T("OriginalPopulatedChannel") );
    tempRLS.duplicate_channel( _T("ReinitializationPopulated"), _T("OriginalPopulatedChannel") );

    tempRLS.reinitialize_signed_distance_from_populated( _T("ReinitializationPopulated") );

    // tempRLS.extrapolate_channels( "ReinitializationPopulated" );

    const_rle_channel_accessor<boost::uint8_t> populatedChannelAccessor =
        tempRLS.get_channel_accessor<boost::uint8_t>( _T("OriginalPopulatedChannel") );

    std::vector<float> newDistanceData( tempRLS.size() );

    const voxel_coord_system& velVCS = velocityField.get_voxel_coord_system();
    const rle_index_spec& ris = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velAcc =
        velocityField.get_channel_accessor<vector3f>( staggeredVelocityChannelName );
    // precompute the index to the center cell (2,2,2) in terms of the local 5x5x5 box
    // this is written out for clarity, but could be collapsed to save a few cyclces...
    size3 boxSize( 7, 7, 7 );
    vector3 centerOffset( 3 );

    const rle_index_spec& tempRIS = tempRLS.get_rle_index_spec();

    // create a data index channel for copying the advected values back into the current level set
    create_data_index_map_channel( _T("AdvectedLevelSetLookup"), tempRIS );

    boundbox3 bounds = tempRIS.outer_bounds();

    const vector3& boundsMin = bounds.minimum();
    const vector3& boundsMax = bounds.maximum();

    int ysize = bounds.ysize();

    FF_LOG( error ) << "Original Bounds: " << m_rleIndex.outer_bounds() << endl;
    FF_LOG( error ) << "Dilated Bounds: " << bounds << endl;

    const vector<boost::int32_t>& bcToRunIndex = tempRIS.get_bc_to_run_index_vector();
    const std::vector<run_data>& runIndexData = tempRIS.get_run_index_data_vector();

    vector3f finalPosition;

    int cellBoxIndex =
        centerOffset.x + centerOffset.y * boxSize.xsize() + centerOffset.z * boxSize.xsize() * boxSize.ysize();

    rle_defined_iterator i = tempRIS.begin(), ie = tempRIS.end();

    // int count = 0;
    // for( ; i!=ie; ++i ) {
    //	if( populatedChannelAccessor[i.get_data_index()] ) {
    //		++count;
    //	}
    // }

    // logging::error << "Populated Count: " << count << "\tTotal: " << tempRLS.size() << "\tOrignal Size: " <<
    // m_distanceData.size() << endl;

    for( int z = boundsMin.z; z <= boundsMax.z; ++z ) {
        int cIndex = tempRIS.z_to_c( z );
        for( int y = boundsMin.y; y <= boundsMax.y; ++y ) {
            // compute the bcIndex of the scanline
            int bcIndex = tempRIS.y_to_b( y ) + cIndex * ysize;

            // logging::error << "Computing Scanline: [" << y << ", " << z << "]: bcIndex:" << bcIndex << endl;
            //  get the run range of the scanline
            int runRangeStart = bcToRunIndex[bcIndex];
            int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

            // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the runs
            // with the iterator
            if( runRangeStart != runRangeEnd || runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                // construct a block iterator for this scanline

                // for velocity advection we need one voxel in the negative direction, and two in hte positive
                // to make sure we have all the data indices that we will need to process the lookup. This means we will
                // end up with more overall indices, but should be faster than carefully choosing and finding only the
                // ones we need
                boundbox3 box( vector3( boundsMin.x - centerOffset.x, y - centerOffset.y, z - centerOffset.z ),
                               boxSize );

                for( levelset::rle_defined_box_iterator boxIter( tempRIS, box );
                     !boxIter.is_xplane_finished( centerOffset.x ); ++boxIter ) {
                    const boost::int32_t* const dataIndices = boxIter.get_indices();

                    /*		for (int i=0; i!= boxSize.volume(); ++i ){
                          logging::error << dataIndices[i] << endl;
                        }*/
                    boost::int32_t index = dataIndices[cellBoxIndex];

                    // if the center cell is not actually defined then we can skip
                    if( index < 0 ) {
                        // logging::error << "undef: "  << currentMin + centerOffset << std::endl;
                        continue;
                    }

                    if( populatedChannelAccessor[index] ) {
                        const boundbox3& currentBox = boxIter.current_box();
                        const vector3& currentMin = currentBox.minimum();

                        vector3 cellVoxelCoord = currentMin + centerOffset;

                        vector3f finalPosition;

                        //		if( i.get_coord().z == bounds.minimum().z)
                        //			logging::set_logging_level(5);

                        vector3f lookup = m_voxelCoordSystem.get_world_voxel_center( cellVoxelCoord );

                        // logging::error  << "\nInitial: (coord) " << i.get_coord() << " (phi) "  <<
                        // m_distanceData[i.get_data_index()] << " step:" << -dt << endl; logging::error  << "Lookup
                        // Voxel Coord: "
                        // << cellVoxelCoord << endl; logging::error  << "Lookup World Coord: " << lookup << endl;
                        fluids::trace_runge_kutta_4( velVCS, ris, velAcc, lookup, -dt, finalPosition );
                        // logging::error  << "Final Voxel Coord: " << m_voxelCoordSystem.get_voxel_coord(finalPosition)
                        // << endl; logging::debug << "lookup index: " << index << endl; newDistanceData[index] =
                        // trilerp_signed_distance( finalPosition );

                        newDistanceData[index] =
                            weno3_signed_distance_lookup( dataIndices, boxSize, currentMin, tempRLS.get_distance_data(),
                                                          m_voxelCoordSystem.get_voxel_coord( finalPosition ) );

                        // logging::debug << "trilerp: " << newDistanceData[index] << endl;
                        // FF_LOG(debug) << "Final: (coord) " << m_voxelCoordSystem.get_voxel_coord( finalPosition );
                        // FF_LOG(debug)<< " (phi) " <<  newDistanceData[i.get_data_index()] << endl;

                        // if( i.get_coord().z == bounds.minimum().z)
                        //	logging::set_logging_level(3);
                    }
                }
            }
        }
    }

    const_rle_channel_accessor<boost::int32_t> indexAccessor =
        get_channel_accessor<boost::int32_t>( _T("AdvectedLevelSetLookup") );

    i = m_rleIndex.begin(), ie = m_rleIndex.end();
    for( ; i != ie; ++i ) {
        boost::int32_t index = i.get_data_index();
        m_distanceData[index] = newDistanceData[indexAccessor[index]];
    }

    // clean up the lookup channel
    erase_channel( _T("AdvectedLevelSetLookup") );
    // m_distanceData.swap(newDistanceData);
}

void rle_level_set::semi_lagrangian_advect( const rle_level_set& velocityField,
                                            const frantic::tstring& velocityChannelName, float dt ) {

    rle_defined_iterator i = m_rleIndex.begin();
    rle_defined_iterator iEnd = m_rleIndex.end();

    std::vector<float> newDistanceData( m_distanceData.size() );

    boundbox3 bounds = m_rleIndex.outer_bounds();

    const voxel_coord_system& velVCS = velocityField.get_voxel_coord_system();
    const rle_index_spec& ris = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velAcc = velocityField.get_channel_accessor<vector3f>( velocityChannelName );

    for( ; i != iEnd; ++i ) {
        vector3f finalPosition;

        vector3f lookup = m_voxelCoordSystem.get_world_voxel_center( i.get_coord() );

        fluids::trace_runge_kutta_2_unstaggered( velVCS, ris, velAcc, lookup, -dt, finalPosition );
        newDistanceData[i.get_data_index()] = trilerp_signed_distance( finalPosition );

        if( i.get_coord().z == bounds.minimum().z )
            logging::set_logging_level( 3 );
    }

    m_distanceData.swap( newDistanceData );
}

void rle_level_set::semi_lagrangian_advect_channels_staggered(
    const rle_level_set& inputRLS, const fluids::rle_voxel_field& inputStaggeredVelocityField, float timeStep,
    const std::vector<frantic::tstring>& channelsToAdvect, const frantic::tstring& targetVoxelsChannelName,
    const frantic::tstring& populatedChannelName ) {
    bool useTargetVoxelsChannel = !targetVoxelsChannelName.empty();
    bool createPopulatedChannel = !populatedChannelName.empty();

    if( useTargetVoxelsChannel && createPopulatedChannel && targetVoxelsChannelName == populatedChannelName )
        throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels_staggered() - The specified target voxels "
                             "channel and output populated channel were requested to be the same name, \"" +
                             frantic::strings::to_string( populatedChannelName ) + "\", these must be different." );

    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;
    const_rle_channel_accessor<boost::uint8_t> inputTargetVoxelsAccessor;
    rle_channel_accessor<boost::uint8_t> outputPopulatedAccessor;

    // This is a memory buffer holding memory set to all 0 values.
    vector<char> zero( 32 );

    // First make sure that all the necesary channels are defined
    for( size_t i = 0, ie = channelsToAdvect.size(); i != ie; ++i ) {
        // Disallow advection of the 'Velocity' channel
        if( channelsToAdvect[i] == _T("Velocity") )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels_staggered() - This function isn't able to advect the "
                "requested \"Velocity\" channel, it needs this channel to remain constant during the advection "
                "process." );
        if( createPopulatedChannel && channelsToAdvect[i] == populatedChannelName )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels_staggered() - This function isn't able to "
                "advect the requested \"" +
                frantic::strings::to_string( populatedChannelName ) +
                "\" channel, because a new such channel is being created for 'Populated' data." );
        // Make sure that the input level set has all the requested channels.
        if( !inputRLS.has_channel( channelsToAdvect[i] ) )
            throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels_staggered() - The input rle_level_set "
                                 "doesn't have a channel named \"" +
                                 frantic::strings::to_string( channelsToAdvect[i] ) +
                                 "\", which was requested as input" );
        inputAccessors.push_back( inputRLS.get_channel_general_accessor( channelsToAdvect[i] ) );

        // If the output already had this data channel, this will reuse it.
        add_channel( channelsToAdvect[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
        outputAccessors.push_back( get_channel_general_accessor( channelsToAdvect[i] ) );

        // If the zero buffer isn't big enough for this channel, make it bigger.
        size_t primitiveSize = inputAccessors.back().primitive_size();
        if( primitiveSize > zero.size() )
            zero.resize( primitiveSize );

        // No need to initialize the channel to all 0's, because all the values get set in the main loop below
    }

    if( useTargetVoxelsChannel ) {
        if( !has_channel( targetVoxelsChannelName ) )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels_staggered() - A target voxels flag channel named \"" +
                frantic::strings::to_string( targetVoxelsChannelName ) +
                "\" was indicated, but the target level set doesn't have such a channel." );
        inputTargetVoxelsAccessor = get_channel_accessor<boost::uint8_t>( targetVoxelsChannelName );
    }

    if( createPopulatedChannel ) {
        // add_channel will create the channel if it doesn't exist,
        // or make sure it has the right arity and type if it doesn't
        add_channel( populatedChannelName, 1, data_type_uint8 );
        outputPopulatedAccessor = get_channel_accessor<boost::uint8_t>( populatedChannelName );
        // No need to initialize the channel to all 0's, because all the values get set in the main loop below
    }

    const voxel_coord_system& vcs = m_voxelCoordSystem;

    float deltas[3];
    float weights[8];
    boost::int32_t dataIndices[8];
    const char* dataPointers[8];
    char* zeroPointer = &zero[0];
    bool isPopulated;

    if( !inputStaggeredVelocityField.has_channel( _T("StaggeredVelocity") ) )
        throw runtime_error(
            "rle_level_set::semi_lagrangian_advect_channels_staggered() - The input velocity field must "
            "have a \"StaggeredVelocity\" channel in order to do a staggered semi-lagrangian advection." );

    const_rle_channel_accessor<vector3f> velocityAccessor =
        inputStaggeredVelocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    const rle_index_spec& velRIS = inputStaggeredVelocityField.get_rle_index_spec();
    const voxel_coord_system& velVCS = inputStaggeredVelocityField.get_voxel_coord_system();

    // Iterate through all the defined voxels, and do the semi-lagrangian backtrace for each one.
    for( rle_defined_iterator i = m_rleIndex.begin(), ie = m_rleIndex.end(); i != ie; ++i ) {
        int dataIndex = i.get_data_index();
        if( !useTargetVoxelsChannel || inputTargetVoxelsAccessor[dataIndex] != 0 ) {
            // vector3f velocity = velocityAccessor[dataIndex];
            vector3f voxelCenter = vcs.get_world_voxel_center( i.get_coord() );
            // Backtrace to sourceCoord.  Note that this simple linear back trace could be made selectably more
            // sophisticated. vector3f sourceCoord = voxelCenter - timeStep * velocity;
            vector3f sourceCoord;
            fluids::trace_runge_kutta_2( velVCS, velRIS, velocityAccessor, voxelCenter, -timeStep, sourceCoord );

            get_trilerp_indices( inputRLS.get_rle_index_spec(),
                                 inputRLS.get_voxel_coord_system().get_voxel_coord( sourceCoord ), deltas,
                                 dataIndices );

            // convert the deltas into our 8 interpolation weights
            get_trilerp_weights( deltas, weights );

            // If requested, store the 'Populated' channel value
            if( createPopulatedChannel ) {
                isPopulated = true;
                for( int j = 0; j != 8; ++j ) {
                    if( dataIndices[j] < 0 ) {
                        isPopulated = false;
                        break;
                    }
                }
                outputPopulatedAccessor[dataIndex] = isPopulated ? 1 : 0;
            }

            // Go through all the channels, and copy the linear interpolated value
            for( size_t j = 0, je = inputAccessors.size(); j != je; ++j ) {
                const_rle_channel_general_accessor& inputAccessor = inputAccessors[j];
                // Get all the input data pointers for this channel
                for( int k = 0; k != 8; ++k ) {
                    int inputDataIndex = dataIndices[k];
                    if( inputDataIndex < 0 ) {
                        dataPointers[k] = zeroPointer;
                    } else {
                        dataPointers[k] = inputAccessor.data( inputDataIndex );
                    }
                }
                // Do the linear interpolation for this channel
                inputAccessor.get_weighted_sum_combine_function()( weights, dataPointers, 8, inputAccessor.arity(),
                                                                   outputAccessors[j].data( dataIndex ) );
            }
        } else {
            // If this wasn't a target voxel, but we're creating a populated channel, then indicate that
            // we didn't populate this channel with data.
            if( createPopulatedChannel )
                outputPopulatedAccessor[dataIndex] = 0;
        }
    }
}

void rle_level_set::semi_lagrangian_advect_channels( const rle_level_set& inputRLS, float timeStep,
                                                     const std::vector<frantic::tstring>& channelsToAdvect,
                                                     const frantic::tstring& targetVoxelsChannelName,
                                                     const frantic::tstring& populatedChannelName ) {
    bool useTargetVoxelsChannel = !targetVoxelsChannelName.empty();
    bool createPopulatedChannel = !populatedChannelName.empty();

    if( useTargetVoxelsChannel && createPopulatedChannel && targetVoxelsChannelName == populatedChannelName )
        throw runtime_error(
            "rle_level_set::semi_lagrangian_advect_channels() - The specified target voxels channel and "
            "output populated channel were requested to be the same name, \"" +
            frantic::strings::to_string( populatedChannelName ) + "\", these must be different." );

    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;
    const_rle_channel_accessor<boost::uint8_t> inputTargetVoxelsAccessor;
    rle_channel_accessor<boost::uint8_t> outputPopulatedAccessor;
    rle_channel_accessor<vector3f> velocityAccessor;

    if( !has_channel( _T("Velocity") ) )
        throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels() - The output level set must have a "
                             "\"Velocity\" channel in order to do a semi-lagrangian advection." );

    velocityAccessor = get_channel_accessor<vector3f>( _T("Velocity") );

    // This is a memory buffer holding memory set to all 0 values.
    vector<char> zero( 32 );

    // First make sure that all the necesary channels are defined
    for( size_t i = 0, ie = channelsToAdvect.size(); i != ie; ++i ) {
        // Disallow advection of the 'Velocity' channel
        if( channelsToAdvect[i] == _T("Velocity") )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels() - This function isn't able to advect the requested "
                "\"Velocity\" channel, it needs this channel to remain constant during the advection process." );
        if( createPopulatedChannel && channelsToAdvect[i] == populatedChannelName )
            throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels() - This function isn't able to "
                                 "advect the requested \"" +
                                 frantic::strings::to_string( populatedChannelName ) +
                                 "\" channel, because a new such channel is being created for 'Populated' data." );
        // Make sure that the input level set has all the requested channels.
        if( !inputRLS.has_channel( channelsToAdvect[i] ) )
            throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels() - The input rle_level_set doesn't "
                                 "have a channel named \"" +
                                 frantic::strings::to_string( channelsToAdvect[i] ) +
                                 "\", which was requested as input" );
        inputAccessors.push_back( inputRLS.get_channel_general_accessor( channelsToAdvect[i] ) );

        // If the output already had this data channel, this will reuse it.
        add_channel( channelsToAdvect[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
        outputAccessors.push_back( get_channel_general_accessor( channelsToAdvect[i] ) );

        // If the zero buffer isn't big enough for this channel, make it bigger.
        size_t primitiveSize = inputAccessors.back().primitive_size();
        if( primitiveSize > zero.size() )
            zero.resize( primitiveSize );

        // No need to initialize the channel to all 0's, because all the values get set in the main loop below
    }

    if( useTargetVoxelsChannel ) {
        if( !has_channel( targetVoxelsChannelName ) )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels() - A target voxels flag channel named \"" +
                frantic::strings::to_string( targetVoxelsChannelName ) +
                "\" was indicated, but the target level set doesn't have such a channel." );
        inputTargetVoxelsAccessor = get_channel_accessor<boost::uint8_t>( targetVoxelsChannelName );
    }

    if( createPopulatedChannel ) {
        // add_channel will create the channel if it doesn't exist,
        // or make sure it has the right arity and type if it doesn't
        add_channel( populatedChannelName, 1, data_type_uint8 );
        outputPopulatedAccessor = get_channel_accessor<boost::uint8_t>( populatedChannelName );
        // No need to initialize the channel to all 0's, because all the values get set in the main loop below
    }

    const voxel_coord_system& vcs = m_voxelCoordSystem;

    float weights[8], deltas[3];
    boost::int32_t dataIndices[8];
    const char* dataPointers[8];
    char* zeroPointer = &zero[0];
    bool isPopulated;

    // Iterate through all the defined voxels, and do the semi-lagrangian backtrace for each one.
    for( rle_defined_iterator i = m_rleIndex.begin(), ie = m_rleIndex.end(); i != ie; ++i ) {
        int dataIndex = i.get_data_index();
        if( !useTargetVoxelsChannel || inputTargetVoxelsAccessor[dataIndex] != 0 ) {
            vector3f velocity = velocityAccessor[dataIndex];
            vector3f voxelCenter = vcs.get_world_voxel_center( i.get_coord() );
            // Backtrace to sourceCoord.  Note that this simple linear back trace could be made selectably more
            // sophisticated.
            vector3f sourceCoord = voxelCenter - timeStep * velocity;

            get_trilerp_indices( inputRLS.get_rle_index_spec(),
                                 inputRLS.get_voxel_coord_system().get_voxel_coord( sourceCoord ), deltas,
                                 dataIndices );

            // convert the deltas into our 8 interpolation weights
            get_trilerp_weights( deltas, weights );

            // If requested, store the 'Populated' channel value
            if( createPopulatedChannel ) {
                isPopulated = true;
                for( int j = 0; j != 8; ++j ) {
                    if( dataIndices[j] < 0 ) {
                        isPopulated = false;
                        break;
                    }
                }
                outputPopulatedAccessor[dataIndex] = isPopulated ? 1 : 0;
            }

            // Go through all the channels, and copy the linear interpolated value
            for( size_t j = 0, je = inputAccessors.size(); j != je; ++j ) {
                const_rle_channel_general_accessor& inputAccessor = inputAccessors[j];
                // Get all the input data pointers for this channel
                for( int k = 0; k != 8; ++k ) {
                    int inputDataIndex = dataIndices[k];
                    if( inputDataIndex < 0 ) {
                        dataPointers[k] = zeroPointer;
                    } else {
                        dataPointers[k] = inputAccessor.data( inputDataIndex );
                    }
                }
                // Do the linear interpolation for this channel
                inputAccessor.get_weighted_sum_combine_function()( weights, dataPointers, 8, inputAccessor.arity(),
                                                                   outputAccessors[j].data( dataIndex ) );
            }

        } else {
            // If this wasn't a target voxel, but we're creating a populated channel, then indicate that
            // we didn't populate this channel with data.
            if( createPopulatedChannel )
                outputPopulatedAccessor[dataIndex] = 0;
        }
    }
}

void rle_level_set::semi_lagrangian_advect_channels_with_undefined_default(
    const rle_level_set& inputRLS, const rle_level_set& defaultRLS, float timeStep,
    const std::vector<frantic::tstring>& channelsToAdvect, const frantic::tstring& targetVoxelsChannelName,
    const frantic::tstring& populatedChannelName ) {
    bool useTargetVoxelsChannel = !targetVoxelsChannelName.empty();
    bool createPopulatedChannel = !populatedChannelName.empty();

    if( useTargetVoxelsChannel && createPopulatedChannel && targetVoxelsChannelName == populatedChannelName )
        throw runtime_error(
            "rle_level_set::semi_lagrangian_advect_channels() - The specified target voxels channel and "
            "output populated channel were requested to be the same name, \"" +
            frantic::strings::to_string( populatedChannelName ) + "\", these must be different." );

    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<const_rle_channel_general_accessor> defaultAccessors;
    vector<rle_channel_general_accessor> outputAccessors;
    const_rle_channel_accessor<boost::uint8_t> inputTargetVoxelsAccessor;
    rle_channel_accessor<boost::uint8_t> outputPopulatedAccessor;
    rle_channel_accessor<vector3f> velocityAccessor;

    if( !has_channel( _T("Velocity") ) )
        throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels() - The output level set must have a "
                             "\"Velocity\" channel in order to do a semi-lagrangian advection." );

    velocityAccessor = get_channel_accessor<vector3f>( _T("Velocity") );

    // This is a memory buffer holding memory set to all 0 values.
    vector<char> zero( 32 );

    // First make sure that all the necesary channels are defined
    for( size_t i = 0, ie = channelsToAdvect.size(); i != ie; ++i ) {
        // Disallow advection of the 'Velocity' channel
        if( channelsToAdvect[i] == _T("Velocity") )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels() - This function isn't able to advect the requested "
                "\"Velocity\" channel, it needs this channel to remain constant during the advection process." );
        if( createPopulatedChannel && channelsToAdvect[i] == populatedChannelName )
            throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels() - This function isn't able to "
                                 "advect the requested \"" +
                                 frantic::strings::to_string( populatedChannelName ) +
                                 "\" channel, because a new such channel is being created for 'Populated' data." );
        // Make sure that the input level set has all the requested channels.
        if( !inputRLS.has_channel( channelsToAdvect[i] ) )
            throw runtime_error( "rle_level_set::semi_lagrangian_advect_channels() - The input rle_level_set doesn't "
                                 "have a channel named \"" +
                                 frantic::strings::to_string( channelsToAdvect[i] ) +
                                 "\", which was requested as input" );
        inputAccessors.push_back( inputRLS.get_channel_general_accessor( channelsToAdvect[i] ) );

        if( !defaultRLS.has_channel( channelsToAdvect[i] ) )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels() - The default rle_level_set doesn't have "
                "a channel named \"" +
                frantic::strings::to_string( channelsToAdvect[i] ) + "\", which was requested as input" );
        defaultAccessors.push_back( defaultRLS.get_channel_general_accessor( channelsToAdvect[i] ) );

        // If the output already had this data channel, then try to reuse it.
        if( has_channel( channelsToAdvect[i] ) ) {
            rle_channel_general_accessor outputAccessor = get_channel_general_accessor( channelsToAdvect[i] );
            if( outputAccessor.data_type() == inputAccessors.back().data_type() &&
                outputAccessor.arity() == inputAccessors.back().arity() ) {
                outputAccessors.push_back( outputAccessor );
            } else {
                erase_channel( channelsToAdvect[i] );
                add_channel( channelsToAdvect[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
                outputAccessors.push_back( get_channel_general_accessor( channelsToAdvect[i] ) );
            }
        } else {
            add_channel( channelsToAdvect[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
            outputAccessors.push_back( get_channel_general_accessor( channelsToAdvect[i] ) );
        }

        // If the zero buffer isn't big enough for this channel, make it bigger.
        size_t primitiveSize = inputAccessors.back().primitive_size();
        if( primitiveSize > zero.size() )
            zero.resize( primitiveSize );

        // No need to initialize the channel to all 0's, because all the values get set in the main loop below
    }

    if( useTargetVoxelsChannel ) {
        if( !has_channel( targetVoxelsChannelName ) )
            throw runtime_error(
                "rle_level_set::semi_lagrangian_advect_channels() - A target voxels flag channel named \"" +
                frantic::strings::to_string( targetVoxelsChannelName ) +
                "\" was indicated, but the target level set doesn't have such a channel." );
        inputTargetVoxelsAccessor = get_channel_accessor<boost::uint8_t>( targetVoxelsChannelName );
    }

    if( createPopulatedChannel ) {
        // Create the 'Populated' channel, or reuse it if it already exists
        add_channel( populatedChannelName, 1, channels::data_type_uint8 );
        outputPopulatedAccessor = get_channel_accessor<boost::uint8_t>( populatedChannelName );
        // No need to initialize the channel to all 0's, because all the values get set in the main loop below
    }

    const voxel_coord_system& vcs = m_voxelCoordSystem;

    float weights[8], deltas[3];
    boost::int32_t dataIndices[8];
    const char* dataPointers[8];
    char* zeroPointer = &zero[0];
    bool isPopulated;

    // Iterate through all the defined voxels, and do the semi-lagrangian backtrace for each one.
    for( rle_defined_iterator i = m_rleIndex.begin(), ie = m_rleIndex.end(); i != ie; ++i ) {
        int dataIndex = i.get_data_index();
        if( !useTargetVoxelsChannel || inputTargetVoxelsAccessor[dataIndex] != 0 ) {
            vector3f velocity = velocityAccessor[dataIndex];
            vector3f voxelCenter = vcs.get_world_voxel_center( i.get_coord() );
            vector3 voxelCoord = i.get_coord();
            // Backtrace to sourceCoord.  Note that this simple linear back trace could be made selectably more
            // sophisticated.
            vector3f sourceCoord = voxelCenter - timeStep * velocity;

            get_trilerp_indices( inputRLS.get_rle_index_spec(),
                                 inputRLS.get_voxel_coord_system().get_voxel_coord( sourceCoord ), deltas,
                                 dataIndices );

            // convert the deltas into our 8 interpolation weights
            get_trilerp_weights( deltas, weights );

            // If requested, store the 'Populated' channel value
            if( createPopulatedChannel ) {
                isPopulated = true;
                for( int j = 0; j != 8; ++j ) {
                    if( dataIndices[j] < 0 ) {
                        isPopulated = false;
                        break;
                    }
                }
                outputPopulatedAccessor[dataIndex] = isPopulated ? 1 : 0;
            }

            // Go through all the channels, and copy the linear interpolated value
            for( size_t j = 0, je = inputAccessors.size(); j != je; ++j ) {
                const_rle_channel_general_accessor& inputAccessor = inputAccessors[j];
                bool zeroFound = false;
                // Get all the input data pointers for this channel
                for( int k = 0; k != 8; ++k ) {
                    int inputDataIndex = dataIndices[k];
                    if( inputDataIndex < 0 ) {
                        dataPointers[k] = zeroPointer;
                        zeroFound = true;
                    } else {
                        dataPointers[k] = inputAccessor.data( inputDataIndex );
                    }
                }

                if( zeroFound ) {
                    const_rle_channel_general_accessor& defaultAccessor = defaultAccessors[j];
                    int defaultDataIndex = defaultRLS.XYZtoDataIndex( voxelCoord );

                    if( defaultDataIndex > 0 ) {

                        memcpy( outputAccessors[j].data( dataIndex ), defaultAccessor.data( defaultDataIndex ),
                                defaultAccessor.primitive_size() );

                    } else {

                        // Do the linear interpolation for this channel
                        inputAccessor.get_weighted_sum_combine_function()(
                            weights, dataPointers, 8, inputAccessor.arity(), outputAccessors[j].data( dataIndex ) );
                    }
                } else {

                    // Do the linear interpolation for this channel
                    inputAccessor.get_weighted_sum_combine_function()( weights, dataPointers, 8, inputAccessor.arity(),
                                                                       outputAccessors[j].data( dataIndex ) );
                }
            }
        } else {
            // If this wasn't a target voxel, but we're creating a populated channel, then indicate that
            // we didn't populate this channel with data.
            if( createPopulatedChannel )
                outputPopulatedAccessor[dataIndex] = 0;
        }
    }
}

static int baseCounter = 0;

void rle_level_set::extrapolate_channels( const std::vector<frantic::tstring>& channelsToExtrapolate,
                                          const std::vector<frantic::tstring>& correspondingGradientChannelsToMatch,
                                          const frantic::tstring& populatedChannelName ) {
    if( channelsToExtrapolate.size() != correspondingGradientChannelsToMatch.size() )
        throw runtime_error( "rle_level_set::extrapolate_channels() - The provided channelsToExtrapolate and"
                             " correspondingGradientChannelsToMatch arrays had different"
                             " sizes, they must be the same size." );

    // If there are no defined voxels, then there's definitely not going to be any extrapolation to do.
    if( m_rleIndex.data_size() <= 1 )
        return;

    // Confirm that there is a 'Populated' channel.
    if( !has_channel( populatedChannelName ) )
        throw runtime_error(
            "rle_level_set::extrapolate_channels() - There must be a \"" +
            frantic::strings::to_string( populatedChannelName ) +
            "\" channel defined in the level set specifying which voxels are extrapolation sources and "
            "which are targets." );

    rle_channel_general_accessor populatedAccessor = get_channel_general_accessor( populatedChannelName );
    // Make sure its a uint8 channel
    if( populatedAccessor.arity() != 1 || populatedAccessor.data_type() != channels::data_type_uint8 )
        throw runtime_error( "rle_level_set::extrapolate_channels() - The  \"" +
                             frantic::strings::to_string( populatedChannelName ) +
                             "\" channel (for the 'Populated' data) must be a uint8 channel, it is set to " +
                             populatedAccessor.type_str() + "." );

    vector<marching_extrapolation_channel> extrapChannels;
    // First make sure that all the necesary channels are defined, and build up the extrapChannels vector
    for( size_t i = 0, ie = channelsToExtrapolate.size(); i != ie; ++i ) {
        if( channelsToExtrapolate[i] == populatedChannelName )
            throw runtime_error( "rle_level_set::extrapolate_channels() - This function isn't able to extrapolate"
                                 " the requested \"" +
                                 frantic::strings::to_string( populatedChannelName ) +
                                 "\" channel, this channel is used to define how extrapolation occurs." );
        if( !has_channel( channelsToExtrapolate[i] ) )
            throw runtime_error(
                "rle_level_set::extrapolate_channels() - The rle_level_set doesn't have a channel named \"" +
                frantic::strings::to_string( channelsToExtrapolate[i] ) +
                "\", for which extrapolation was requested." );

        extrapChannels.push_back( marching_extrapolation_channel() );
        marching_extrapolation_channel& mec = extrapChannels.back();
        rle_channel_general_accessor accessor = get_channel_general_accessor( channelsToExtrapolate[i] );
        mec.arity = accessor.arity();
        mec.data = accessor.data( 0 );
        mec.primitiveSize = accessor.primitive_size();
        mec.weightedSumFn = accessor.get_weighted_sum_combine_function();
        if( correspondingGradientChannelsToMatch[i].empty() ) {
            mec.gradientToMatch = 0;
        } else {
            // Validate that the channel exists
            if( !has_channel( correspondingGradientChannelsToMatch[i] ) )
                throw runtime_error( "rle_level_set::extrapolate_channels() - The rle_level_set doesn't have a"
                                     " channel named \"" +
                                     frantic::strings::to_string( correspondingGradientChannelsToMatch[i] ) +
                                     "\", requested as a gradient to match for channel \"" +
                                     frantic::strings::to_string( channelsToExtrapolate[i] ) + "\"." );

            // Only float types can be extrapolated while matching a gradient
            if( !channels::is_channel_data_type_float( accessor.data_type() ) )
                throw runtime_error( "rle_level_set::extrapolate_channels() - The rle_level_set"
                                     " channel named \"" +
                                     frantic::strings::to_string( channelsToExtrapolate[i] ) +
                                     "\"specified to be extrapolated matching gradient channel \"" +
                                     frantic::strings::to_string( correspondingGradientChannelsToMatch[i] ) +
                                     "\","
                                     " does not have a floating point data type.  Its data type is " +
                                     accessor.type_str() + "." );

            rle_channel_general_accessor gradAcc =
                get_channel_general_accessor( correspondingGradientChannelsToMatch[i] );

            // Validate that the requested gradient channel is compatible with the channel being extrapolated
            if( gradAcc.arity() != 3 * mec.arity )
                throw runtime_error( "rle_level_set::extrapolate_channels() - The rle_level_set"
                                     " channel named \"" +
                                     frantic::strings::to_string( correspondingGradientChannelsToMatch[i] ) +
                                     "\", has a different arity than 3 times the arity of channel \"" +
                                     frantic::strings::to_string( channelsToExtrapolate[i] ) +
                                     "\","
                                     " for which it was requested as a gradient to match."
                                     " (" +
                                     boost::lexical_cast<string>( gradAcc.arity() ) + " instead of " +
                                     boost::lexical_cast<string>( 3 * mec.arity ) + ")" );
            if( gradAcc.data_type() != accessor.data_type() )
                throw runtime_error(
                    "rle_level_set::extrapolate_channels() - The rle_level_set"
                    " channel named \"" +
                    frantic::strings::to_string( correspondingGradientChannelsToMatch[i] ) +
                    "\" has a different data type than the channel it was requested as a gradient to match, \"" +
                    frantic::strings::to_string( channelsToExtrapolate[i] ) +
                    "\"."
                    " (" +
                    frantic::strings::to_string( channels::channel_data_type_str( gradAcc.data_type() ) ) +
                    " instead of " +
                    frantic::strings::to_string( channels::channel_data_type_str( accessor.data_type() ) ) + ")" );

            // Get a pointer to the gradient data
            mec.gradientToMatch = gradAcc.data( 0 );
        }
    }

    const_rle_channel_general_accessor temp = get_channel_general_accessor( populatedChannelName );

    extrapolation_debug_info debugInfo;

    debugInfo.set( m_voxelCoordSystem, m_rleIndex, temp, 0.5f );

    if( logging::is_logging_debug() ) {
        debugInfo.counter = baseCounter;
        // std::cout << " Extrapolation Debug Info -- Counter:\t" << debugInfo.counter << "\tBase Counter: " <<
        // baseCounter
        // <<  std::endl;
    }

    marching_extrapolation_iterative( m_rleIndex.get_cached_adjacency(), m_voxelCoordSystem.voxel_length(),
                                      reinterpret_cast<unsigned char*>( populatedAccessor.data( 0 ) ),
                                      &m_distanceData[0], extrapChannels, debugInfo );

    if( logging::is_logging_debug() ) {
        baseCounter = debugInfo.counter;
        // std::cout << " Extrapolation Debug Info -- Counter:\t" << debugInfo.counter << "\tBase Counter: " <<
        // baseCounter
        // <<  std::endl;
    }

    // debugInfo.reset();
}

void rle_level_set::extrapolate_channel( const frantic::tstring& channelToExtrapolate,
                                         const frantic::tstring& gradientChannelToMatch,
                                         const frantic::tstring& populatedChannelName ) {
    std::vector<frantic::tstring> channelsToExtrapolate( 1, channelToExtrapolate );
    std::vector<frantic::tstring> gradientChannelsToMatch( 1, gradientChannelToMatch );
    extrapolate_channels( channelsToExtrapolate, gradientChannelsToMatch, populatedChannelName );
}

void rle_level_set::extrapolate_channels( const std::vector<frantic::tstring>& channelsToExtrapolate,
                                          const frantic::tstring& populatedChannelName ) {
    // If there are no defined voxels, then there's definitely not going to be any extrapolation to do.
    if( m_rleIndex.data_size() <= 1 )
        return;

    // Confirm that there is a 'Populated' channel.
    if( !has_channel( populatedChannelName ) )
        throw runtime_error(
            "rle_level_set::extrapolate_channels() - There must be a \"" +
            frantic::strings::to_string( populatedChannelName ) +
            "\" channel defined in the level set specifying which voxels are extrapolation sources and "
            "which are targets." );
    rle_channel_general_accessor populatedAccessor = get_channel_general_accessor( populatedChannelName );
    // Make sure its a uint8 channel
    if( populatedAccessor.arity() != 1 || populatedAccessor.data_type() != channels::data_type_uint8 )
        throw runtime_error( "rle_level_set::extrapolate_channels() - The  \"" +
                             frantic::strings::to_string( populatedChannelName ) +
                             "\" channel (for the 'Populated' data) must be a uint8 channel, it is set to " +
                             populatedAccessor.type_str() + "." );

    vector<marching_extrapolation_channel> extrapChannels;
    // First make sure that all the necesary channels are defined, and build up the extrapChannels vector
    for( size_t i = 0, ie = channelsToExtrapolate.size(); i != ie; ++i ) {
        if( channelsToExtrapolate[i] == populatedChannelName )
            throw runtime_error(
                "rle_level_set::extrapolate_channels() - This function isn't able to extrapolate the requested \"" +
                frantic::strings::to_string( populatedChannelName ) +
                "\" channel, this channel is used to define how extrapolation occurs." );
        if( !has_channel( channelsToExtrapolate[i] ) )
            throw runtime_error(
                "rle_level_set::extrapolate_channels() - The rle_level_set doesn't have a channel named \"" +
                frantic::strings::to_string( channelsToExtrapolate[i] ) +
                "\", for which extrapolation was requested." );
        extrapChannels.push_back( marching_extrapolation_channel() );
        marching_extrapolation_channel& mec = extrapChannels.back();
        rle_channel_general_accessor accessor = get_channel_general_accessor( channelsToExtrapolate[i] );
        mec.arity = accessor.arity();
        mec.data = accessor.data( 0 );
        mec.gradientToMatch = 0;
        mec.primitiveSize = accessor.primitive_size();
        mec.weightedSumFn = accessor.get_weighted_sum_combine_function();
    }

    const_rle_channel_general_accessor temp = get_channel_general_accessor( populatedChannelName );

    extrapolation_debug_info debugInfo;
    debugInfo.set( m_voxelCoordSystem, m_rleIndex, temp, 0.5f );

    // add_channel<std::size_t>( "PopulatedOrder" );
    // zero_channel( "PopulatedOrder" );
    // rle_channel_accessor<std::size_t> populatedOrderAccessor = get_channel_accessor<std::size_t>( "PopulatedOrder" );
    marching_extrapolation_iterative( m_rleIndex.get_cached_adjacency(), m_voxelCoordSystem.voxel_length(),
                                      reinterpret_cast<unsigned char*>( populatedAccessor.data( 0 ) ),
                                      &m_distanceData[0], extrapChannels, debugInfo );

    debugInfo.reset();
}

void rle_level_set::extrapolate_channel( const frantic::tstring& channelToExtrapolate,
                                         const frantic::tstring& populatedChannelName ) {
    std::vector<frantic::tstring> channelsToExtrapolate( 1, channelToExtrapolate );
    extrapolate_channels( channelsToExtrapolate, populatedChannelName );
}

void rle_level_set::extrapolate_channels( const frantic::tstring& populatedChannelName ) {
    // If there are no defined voxels, then there's definitely not going to be any extrapolation to do.
    if( m_rleIndex.data_size() <= 1 )
        return;

    // Confirm that there is a 'Populated' channel.
    if( !has_channel( populatedChannelName ) )
        throw runtime_error(
            "rle_level_set::extrapolate_channels() - There must be a \"" +
            frantic::strings::to_string( populatedChannelName ) +
            "\" channel defined in the level set specifying which voxels are extrapolation sources and "
            "which are targets." );
    rle_channel_general_accessor populatedAccessor = get_channel_general_accessor( populatedChannelName );
    // Make sure its a uint8 channel
    if( populatedAccessor.arity() != 1 || populatedAccessor.data_type() != channels::data_type_uint8 )
        throw runtime_error( "rle_level_set::extrapolate_channels() - The  \"" +
                             frantic::strings::to_string( populatedChannelName ) +
                             "\" channel (for the 'Populated' data) must be a uint8 channel, it is set to " +
                             frantic::strings::to_string( channels::channel_data_type_str(
                                 populatedAccessor.arity(), populatedAccessor.data_type() ) ) +
                             "." );

    vector<marching_extrapolation_channel> extrapChannels;
    // First make sure that all the necesary channels are defined, and build up the extrapChannels vector
    for( std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.begin(), ie = m_namedChannels.end();
         i != ie; ++i ) {
        if( i->first != populatedChannelName ) {
            extrapChannels.push_back( marching_extrapolation_channel() );
            marching_extrapolation_channel& mec = extrapChannels.back();
            mec.arity = i->second.arity();
            mec.data = i->second.data();
            mec.gradientToMatch = 0;
            mec.primitiveSize = i->second.primitive_size();
            mec.weightedSumFn = channel_weighted_sum_combine_function( i->second.data_type() );
        }
    }

    const_rle_channel_general_accessor temp = get_channel_general_accessor( populatedChannelName );

    extrapolation_debug_info debugInfo;
    debugInfo.set( m_voxelCoordSystem, m_rleIndex, temp, 0.5f );

    marching_extrapolation_iterative( m_rleIndex.get_cached_adjacency(), m_voxelCoordSystem.voxel_length(),
                                      reinterpret_cast<unsigned char*>( populatedAccessor.data( 0 ) ),
                                      &m_distanceData[0], extrapChannels, debugInfo );

    debugInfo.reset();
}

void rle_level_set::extrapolate_staggered_channel( rle_voxel_field& rvf,
                                                   const frantic::tstring& rvfStaggeredChannelName, int extrapDirection,
                                                   const frantic::tstring& rlsStaggeredPopulatedChannelName,
                                                   const frantic::tstring& rlsDataIndexMapChannelName ) {
    // Both fields must have defined voxels for this to make any sense
    if( m_rleIndex.data_size() <= 1 || rvf.get_rle_index_spec().data_size() <= 1 )
        return;

    // Make sure the voxel coordinate systems are the same
    if( !m_voxelCoordSystem.equals( rvf.get_voxel_coord_system() ) )
        throw runtime_error(
            "rle_level_set::extrapolate_staggered_channel() - The rle_level_set and rle_voxel_field provided "
            "have differing voxel coordinate systems (" +
            m_voxelCoordSystem.str() + " versus " + rvf.get_voxel_coord_system().str() + ")" );

    if( !rvf.has_channel( rvfStaggeredChannelName ) )
        throw runtime_error( "rle_level_set::extrapolate_staggered_channel() - The rle_voxel_field provided "
                             "did not have the staggered field specified, \"" +
                             frantic::strings::to_string( rvfStaggeredChannelName ) + "\"." );

    if( !has_channel( rlsStaggeredPopulatedChannelName ) )
        throw runtime_error( "rle_level_set::extrapolate_staggered_channel() - The rle_level_set "
                             "did not have the staggered populated channel specified, \"" +
                             frantic::strings::to_string( rlsStaggeredPopulatedChannelName ) + "\"." );

    rle_channel_general_accessor staggeredChannelAccessor;
    staggeredChannelAccessor = rvf.get_channel_general_accessor( rvfStaggeredChannelName );
    if( staggeredChannelAccessor.arity() != 3 || staggeredChannelAccessor.data_type() != data_type_float32 )
        throw runtime_error( "rle_level_set::extrapolate_staggered_channel() - The channel specified for the staggered "
                             "extrapolation, \"" +
                             frantic::strings::to_string( rvfStaggeredChannelName ) + "\", has the incorrect type " +
                             staggeredChannelAccessor.type_str() + ", its type should instead be float32[3]." );
    float* staggeredExtrapChannel = reinterpret_cast<float*>( staggeredChannelAccessor.data( 0 ) );

    rle_channel_general_accessor staggeredPopulatedChannelAccessor;
    staggeredPopulatedChannelAccessor = get_channel_general_accessor( rlsStaggeredPopulatedChannelName );
    if( staggeredPopulatedChannelAccessor.arity() != 1 ||
        staggeredPopulatedChannelAccessor.data_type() != data_type_uint8 )
        throw runtime_error( "rle_level_set::extrapolate_staggered_channel() - The channel specified for the staggered "
                             "populated channel \"" +
                             frantic::strings::to_string( rlsStaggeredPopulatedChannelName ) +
                             "\", has the incorrect type " + staggeredPopulatedChannelAccessor.type_str() +
                             ", its type should instead be uint8." );
    boost::uint8_t* staggeredPopulatedChannel =
        reinterpret_cast<boost::uint8_t*>( staggeredPopulatedChannelAccessor.data( 0 ) );

    boost::int32_t* dataIndexMapChannel;
    rle_channel_general_accessor dataIndexMapChannelAccessor;
    // Get the dataIndexMapChannel, either by creating it now if the specified channel is ""
    // or grabbing the channel otherwise.
    vector<boost::int32_t> dataIndexMapChannelVector;
    if( rlsDataIndexMapChannelName.empty() ) {
        dataIndexMapChannelVector.resize( m_rleIndex.data_size() );
        dataIndexMapChannel = &dataIndexMapChannelVector[0];
        m_rleIndex.fill_data_index_map( rvf.get_rle_index_spec(), dataIndexMapChannel );
    } else {
        if( !has_channel( rlsDataIndexMapChannelName ) )
            throw runtime_error(
                "rle_level_set::extrapolate_staggered_channel() - The channel specified for the data index map \"" +
                frantic::strings::to_string( rlsDataIndexMapChannelName ) +
                "\", did not exist in this rle level set." );
        dataIndexMapChannelAccessor = get_channel_general_accessor( rlsDataIndexMapChannelName );
        if( dataIndexMapChannelAccessor.arity() != 1 || dataIndexMapChannelAccessor.data_type() != data_type_int32 )
            throw runtime_error(
                "rle_level_set::extrapolate_staggered_channel() - The channel specified for the data index map \"" +
                frantic::strings::to_string( rlsDataIndexMapChannelName ) + "\", has the incorrect type " +
                dataIndexMapChannelAccessor.type_str() + ", its type should instead be int32." );
        dataIndexMapChannel = reinterpret_cast<boost::int32_t*>( dataIndexMapChannelAccessor.data( 0 ) );
    }

    // logging::error << "**** (pre LS stagg extrapolate) staggered velocity channel size=" <<
    // staggeredChannelAccessor.size() << " ****" << endl;

    // FF_LOG(error)  << "*** (pre LS stagg extrapolate)  *** field size = " <<  rvf.size() << " ***\n";
    // vector<string> names;
    // rvf.get_channel_names(names);
    // for( size_t i =0 ;i<names.size(); ++i) {
    //	const_rle_channel_general_accessor acc = rvf.get_channel_general_accessor(names[i]);
    //	FF_LOG(error)  << "*** (pre LS stagg extrapolate)  *** "<< names[i] <<  " size = " <<  acc.size() << " ***\n";
    // }
    // FF_LOG(error) << endl;

    // if( !rvf.check_consistency( cout ) )
    //	throw runtime_error("staggered field extrapolation() - (pre) Voxel Field Failed consistency check");

    //	const_rle_channel_general_accessor temp = get_channel_general_accessor(rlsStaggeredPopulatedChannelName);

    // Finally, call the staggered field marching extrapolation function with the collected arguments
    staggered_field_marching_extrapolation_iterative( m_rleIndex.get_cached_adjacency(), &m_distanceData[0],
                                                      staggeredPopulatedChannel, extrapDirection, dataIndexMapChannel,
                                                      staggeredExtrapChannel );

    // if( !rvf.check_consistency( cout ) )
    //	throw runtime_error("staggered field extrapolation() - (post) Voxel Field Failed consistency check");

    // FF_LOG(error)  << "*** (post LS stagg extrapolate)  *** field size = " <<  rvf.size() << " ***\n";
    ////vector<string> names;
    // names.clear();
    // rvf.get_channel_names(names);
    // for( size_t i =0 ;i<names.size(); ++i) {
    //	const_rle_channel_general_accessor acc = rvf.get_channel_general_accessor(names[i]);
    //	FF_LOG(error)  << "*** (post LS stagg extrapolate)  *** "<< names[i] <<  " size = " <<  acc.size() << " ***\n";
    // }
    // FF_LOG(error) << endl;
    erase_channel( _T("StaggeredFlagChannel") );
}

/**
 * Low level function which tags the interface voxels, given raw pointers to the signed distance data and
 * populated channel output data.
 *
 * The 'Populated' channel stores a value of 1 for points that are next to the interface, and 0 for all other points.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 */
void tag_interface_voxels( const ris_adjacency& adj, const float* signedDistanceChannel,
                           unsigned char* populatedChannel ) {
    // Initialize all the populated values to 0
    memset( populatedChannel, 0, adj.data_size() );

    // Iterate and set populated channel to 1 if there is an adjacent neighbour on the other side of the interface
    for( size_t i = 0, ie = adj.data_size(); i != ie; ++i ) {
        const ris_adj_entry& rae = adj[i];
        if( signedDistanceChannel[i] >= 0 ) {
            if( rae.x_pos >= 0 && signedDistanceChannel[rae.x_pos] <= 0 ) {
                populatedChannel[rae.x_pos] = 1;
                populatedChannel[i] = 1;
            }
            if( rae.y_pos >= 0 && signedDistanceChannel[rae.y_pos] <= 0 ) {
                populatedChannel[rae.y_pos] = 1;
                populatedChannel[i] = 1;
            }
            if( rae.z_pos >= 0 && signedDistanceChannel[rae.z_pos] <= 0 ) {
                populatedChannel[rae.z_pos] = 1;
                populatedChannel[i] = 1;
            }
        } else {
            if( rae.x_pos >= 0 && signedDistanceChannel[rae.x_pos] >= 0 ) {
                populatedChannel[rae.x_pos] = 1;
                populatedChannel[i] = 1;
            }
            if( rae.y_pos >= 0 && signedDistanceChannel[rae.y_pos] >= 0 ) {
                populatedChannel[rae.y_pos] = 1;
                populatedChannel[i] = 1;
            }
            if( rae.z_pos >= 0 && signedDistanceChannel[rae.z_pos] >= 0 ) {
                populatedChannel[rae.z_pos] = 1;
                populatedChannel[i] = 1;
            }
        }
    }
}

void rle_level_set::tag_interface_voxels( const frantic::tstring& populatedChannelName ) {
    // Create the 'Populated' channel, or reuse it if it already exists
    add_channel( populatedChannelName, 1, channels::data_type_uint8 );

    // If there are no defined voxels, then skip out early so we can assume there are some later in the function.
    if( m_rleIndex.data_size() == 0 )
        return;

    rle_channel_general_accessor populatedAccessor;

    populatedAccessor = get_channel_general_accessor( populatedChannelName );

    // Get a raw pointer to the populated channel
    unsigned char* populatedChannel = reinterpret_cast<unsigned char*>( populatedAccessor.data( 0 ) );
    const float* signedDistanceChannel = &m_distanceData[0];

    ::tag_interface_voxels( m_rleIndex.get_cached_adjacency(), signedDistanceChannel, populatedChannel );
}

void rle_level_set::tag_distance_range( const frantic::tstring& populatedChannelName, float minVoxelDistance,
                                        float maxVoxelDistance ) {
    // Create the 'Populated' channel, or reuse it if it already exists
    add_channel( populatedChannelName, 1, channels::data_type_uint8 );

    // If there are no defined voxels, then skip out early so we can assume there are some later in the function.
    if( m_rleIndex.data_size() == 0 )
        return;

    rle_channel_general_accessor populatedAccessor;

    populatedAccessor = get_channel_general_accessor( populatedChannelName );

    // Convert the voxel distances into world distances for comparison with the signed distance channel
    float minDistance = minVoxelDistance * m_voxelCoordSystem.voxel_length();
    float maxDistance = maxVoxelDistance * m_voxelCoordSystem.voxel_length();

    // Get a raw pointer to the populated channel, so we can access it faster
    unsigned char* populatedChannel = reinterpret_cast<unsigned char*>( populatedAccessor.data( 0 ) );
    const float* signedDistanceChannel = &m_distanceData[0];

    for( size_t i = 0, ie = m_rleIndex.data_size(); i != ie; ++i ) {
        float distance = signedDistanceChannel[i];
        if( minDistance <= distance && distance <= maxDistance )
            populatedChannel[i] = 1;
        else
            populatedChannel[i] = 0;
    }
}

void rle_level_set::tag_distance_range_staggered( const frantic::tstring& staggeredPopulatedChannelName,
                                                  float minVoxelDistance, float maxVoxelDistance ) {

    bool minIsInfinite = false;
    if( math::is_infinite( minVoxelDistance ) ) {
        minIsInfinite = ( minVoxelDistance < 0 );
    }

    bool maxIsInfinite = false;
    if( math::is_infinite( maxVoxelDistance ) ) {
        maxIsInfinite = ( maxVoxelDistance > 0 );
    }

    // std::cout << "minDistance: " << minVoxelDistance << " inf: " << minIsInfinite << " maxDistance: " <<
    // maxVoxelDistance << " inf: " << maxIsInfinite << std::endl;

    // Create the 'StaggeredPopulated' channel, or reuse it if it already exists
    add_channel( staggeredPopulatedChannelName, 1, channels::data_type_uint8 );

    // If there are no defined voxels, then skip out early so we can assume there are some later in the function.
    if( m_rleIndex.data_size() == 0 )
        return;

    rle_channel_general_accessor staggeredPopulatedAccessor;
    staggeredPopulatedAccessor = get_channel_general_accessor( staggeredPopulatedChannelName );

    // Convert the voxel distances into world distances for comparison with the signed distance channel
    float minDistance = minVoxelDistance * m_voxelCoordSystem.voxel_length();
    float maxDistance = maxVoxelDistance * m_voxelCoordSystem.voxel_length();

    // Get a raw pointer to the populated channel, so we can access it faster
    unsigned char* staggeredPopulatedChannel = reinterpret_cast<unsigned char*>( staggeredPopulatedAccessor.data( 0 ) );
    const float* signedDistanceChannel = &m_distanceData[0];

    const ris_adjacency& adj = m_rleIndex.get_cached_adjacency();

    for( size_t i = 0, ie = m_rleIndex.data_size(); i != ie; ++i ) {
        float iDistance = signedDistanceChannel[i];
        //		const ris_adj_entry& rae = adj[i];
        boost::uint8_t staggeredPopulated = 0;
        // X face
        int iNeighbor = adj[i].x_neg;
        if( iNeighbor >= 0 ) {
            float distance = 0.5f * ( iDistance + signedDistanceChannel[iNeighbor] );
            if( minDistance <= distance && distance <= maxDistance )
                staggeredPopulated |= 0x01;
        } else if( iNeighbor < -1 ) { // inside region code
            if( minIsInfinite ) {
                staggeredPopulated |= 0x01;
            }
        } else { // outside region code
            if( maxIsInfinite ) {
                staggeredPopulated |= 0x01;
            }
        }

        // Y face
        iNeighbor = adj[i].y_neg;
        if( iNeighbor >= 0 ) {
            float distance = 0.5f * ( iDistance + signedDistanceChannel[iNeighbor] );
            if( minDistance <= distance && distance <= maxDistance )
                staggeredPopulated |= 0x02;
        } else if( iNeighbor < -1 ) { // inside region code
            if( minIsInfinite ) {
                staggeredPopulated |= 0x02;
            }
        } else { // outside region code
            if( maxIsInfinite ) {
                staggeredPopulated |= 0x02;
            }
        }

        // Z face
        iNeighbor = adj[i].z_neg;
        if( iNeighbor >= 0 ) {
            float distance = 0.5f * ( iDistance + signedDistanceChannel[iNeighbor] );
            if( minDistance <= distance && distance <= maxDistance )
                staggeredPopulated |= 0x04;
        } else if( iNeighbor < -1 ) { // inside region code
            if( minIsInfinite ) {
                staggeredPopulated |= 0x04;
            }
        } else { // outside region code
            if( maxIsInfinite ) {
                staggeredPopulated |= 0x04;
            }
        }

        // now check the undefined regions in the positive direction ( this is because we only look at the neg adjacent
        // voxels above)
        //( not the right flags... ! )
        if( minIsInfinite ) {
            if( adj[i].x_pos < -1 )
                staggeredPopulated |= 0x01;
            if( adj[i].y_pos < -1 )
                staggeredPopulated |= 0x02;
            if( adj[i].z_pos < -1 )
                staggeredPopulated |= 0x04;
        }

        if( maxIsInfinite ) {
            if( adj[i].x_pos == -1 )
                staggeredPopulated |= 0x01;
            if( adj[i].y_pos == -1 )
                staggeredPopulated |= 0x02;
            if( adj[i].z_pos == -1 )
                staggeredPopulated |= 0x04;
        }

        // Set the channel value.
        staggeredPopulatedChannel[i] = staggeredPopulated;
    }
}

void rle_level_set::reinitialize_signed_distance( frantic::logging::progress_logger& progressLogger,
                                                  const frantic::tstring& populatedChannelToCreate,
                                                  const float insideVoxelDistance, const float outsideVoxelDistance ) {
    if( insideVoxelDistance > 0 )
        throw runtime_error( "rle_level_set::reinitialize_signed_distance() - The insideVoxelDistance value must be "
                             "negative or zero, but instead its value is " +
                             boost::lexical_cast<string>( insideVoxelDistance ) + "." );
    if( outsideVoxelDistance < 0 )
        throw runtime_error( "rle_level_set::reinitialize_signed_distance() - The outsideVoxelDistance value must be "
                             "positive or zero, but instead its value is " +
                             boost::lexical_cast<string>( outsideVoxelDistance ) + "." );

    if( !populatedChannelToCreate.empty() ) {
        // find the initial band, and put it into the 'Populated' channel
        tag_interface_voxels( populatedChannelToCreate );

        // If there are no defined voxels, then there's definitely not going to be any reinitialization to do.
        // NOTE: We're still doing this test after creating the populated channel, so that that channel exists
        //       when the function returns.
        if( m_rleIndex.data_size() == 0 )
            return;

        // tag_interface_voxels created a uint8_t channel
        rle_channel_accessor<boost::uint8_t> populatedAccessor =
            get_channel_accessor<boost::uint8_t>( populatedChannelToCreate );

        // reinitialize from that channel
        fast_marching_reinitialization( progressLogger, m_rleIndex, &populatedAccessor[0], &m_distanceData[0],
                                        m_voxelCoordSystem.voxel_length(), true,
                                        -insideVoxelDistance * m_voxelCoordSystem.voxel_length(),
                                        outsideVoxelDistance * m_voxelCoordSystem.voxel_length() );

    } else {
        // If there are no defined voxels, then there's definitely not going to be any reinitialization to do.
        if( m_rleIndex.data_size() == 0 )
            return;

        const ris_adjacency& adj = m_rleIndex.get_cached_adjacency();

        // find the initial band, and put it into a vector
        std::vector<unsigned char> populatedChannelVector( adj.data_size() );
        ::tag_interface_voxels( adj, &m_distanceData[0], &populatedChannelVector[0] );

        // reintialize with that vector as the populated channel
        fast_marching_reinitialization( progressLogger, m_rleIndex, &populatedChannelVector[0], &m_distanceData[0],
                                        m_voxelCoordSystem.voxel_length(), true,
                                        -insideVoxelDistance * m_voxelCoordSystem.voxel_length(),
                                        outsideVoxelDistance * m_voxelCoordSystem.voxel_length() );
    }
}

void rle_level_set::reinitialize_signed_distance_from_populated( frantic::logging::progress_logger& progressLogger,
                                                                 const frantic::tstring& populatedChannelName,
                                                                 const float insideVoxelDistance,
                                                                 const float outsideVoxelDistance ) {
    // If there are no defined voxels, then there's definitely not going to be any extrapolation to do.
    if( m_rleIndex.data_size() == 0 )
        return;

    if( insideVoxelDistance > 0 )
        throw runtime_error(
            "rle_level_set::reinitialize_signed_distance_from_populated() - The insideVoxelDistance value "
            "must be negative or zero, but instead its value is " +
            boost::lexical_cast<string>( insideVoxelDistance ) + "." );
    if( outsideVoxelDistance < 0 )
        throw runtime_error( "rle_level_set::reinitialize_signed_distance_from_populated() - The outsideVoxelDistance "
                             "value must be positive or zero, but instead its value is " +
                             boost::lexical_cast<string>( outsideVoxelDistance ) + "." );

    rle_channel_general_accessor populatedAccessor = get_channel_general_accessor( populatedChannelName );
    fast_marching_reinitialization(
        progressLogger, m_rleIndex, reinterpret_cast<unsigned char*>( populatedAccessor.data( 0 ) ), &m_distanceData[0],
        m_voxelCoordSystem.voxel_length(), false, -insideVoxelDistance * m_voxelCoordSystem.voxel_length(),
        outsideVoxelDistance * m_voxelCoordSystem.voxel_length() );
}

bool rle_level_set::find_defined_times_along_ray( const ray3f& ray, double tMin, double tMax,
                                                  unsigned int numBisectionSteps, double& tMinOut, double& tMaxOut ) {
    double tMinTemp = tMin;
    double tMaxTemp = tMax;
    int dataIndexMin = XYZtoDataIndex(
        vector3( int( ray.at( tMin ).x + 0.5 ), int( ray.at( tMin ).y + 0.5 ), int( ray.at( tMin ).z + 0.5 ) ) );
    int dataIndexMax = XYZtoDataIndex(
        vector3( int( ray.at( tMax ).x + 0.5 ), int( ray.at( tMax ).y + 0.5 ), int( ray.at( tMax ).z + 0.5 ) ) );

    // if the two points are in the same type of undefined region, assume that there is not a surface between them.
    if( dataIndexMin < 0 && dataIndexMin == dataIndexMax ) {
        return false;
    }

    // if tMin is in an undefined region, do a binary search to find a point on the same side of the surface that is
    // defined
    if( dataIndexMin < 0 ) {

        bool isInside = dataIndexMin == -2;
        bool isFound = false;

        double deltaT = ( tMax - tMin ) / 2;
        tMinTemp = tMin + deltaT;
        int dataIndexTemp;

        for( unsigned i = 0; i < numBisectionSteps && !isFound; ++i ) {
            dataIndexTemp =
                XYZtoDataIndex( vector3( int( ray.at( tMinTemp ).x + 0.5 ), int( ray.at( tMinTemp ).y + 0.5 ),
                                         int( ray.at( tMinTemp ).z + 0.5 ) ) );
            deltaT = deltaT / 2;

            int caseNum = isInside ? 1 : 0;       // is the original point inside?
            caseNum += dataIndexTemp > 0 ? 2 : 0; // are we in a defined region?
            if( dataIndexTemp > 0 ) {
                // pointTemp is defined
                caseNum += m_distanceData[dataIndexTemp] < 0 ? 4 : 0; // is pointTemp inside?
            } else {
                // pointTemp is undefined
                caseNum += dataIndexTemp == -2 ? 4 : 0; // is pointTemp inside?
            }

            switch( caseNum ) {
            case 0: // temp is outside undefined, started at an outside point
                tMinTemp += deltaT;
                break;
            case 1: // temp is outside undefined, started at an inside point
                tMinTemp -= deltaT;
                break;
            case 2: // temp is outside defined, started at an outside point
                isFound = true;
                break;
            case 3: // temp is outside defined, started at an inside point
                tMinTemp -= deltaT;
                break;
            case 4: // temp is inside undefined, started at an outside point
                tMinTemp -= deltaT;
                break;
            case 5: // temp is inside undefined, started at an inside point
                tMinTemp += deltaT;
                break;
            case 6: // temp is inside defined, started at an outside point
                tMinTemp -= deltaT;
                break;
            case 7: // temp is inside defined, started at an inside point
                isFound = true;
                break;
            }
        }
        if( !isFound ) {
            return false;
        }
    }
    tMinOut = tMinTemp;

    // if tMax is in an undefined region, do a binary search to find a point on the same side of the surface that is
    // defined
    if( dataIndexMax < 0 ) {

        bool isInside = dataIndexMax == -2;
        bool isFound = false;

        double deltaT = ( tMax - tMin ) / 2;
        tMaxTemp = tMin + deltaT;
        int dataIndexTemp;

        for( unsigned i = 0; i < numBisectionSteps && !isFound; ++i ) {
            dataIndexTemp =
                XYZtoDataIndex( vector3( int( ray.at( tMaxTemp ).x + 0.5 ), int( ray.at( tMaxTemp ).y + 0.5 ),
                                         int( ray.at( tMaxTemp ).z + 0.5 ) ) );
            deltaT = deltaT / 2;

            int caseNum = isInside ? 1 : 0;       // is the original point inside?
            caseNum += dataIndexTemp > 0 ? 2 : 0; // are we in a defined region?
            if( dataIndexTemp > 0 ) {
                // pointTemp is defined
                caseNum += m_distanceData[dataIndexTemp] < 0 ? 4 : 0; // is pointTemp inside?
            } else {
                // pointTemp is undefined
                caseNum += dataIndexTemp == -2 ? 4 : 0; // is pointTemp inside?
            }

            switch( caseNum ) {
            case 0: // temp is outside undefined, started at an outside point
                tMaxTemp -= deltaT;
                break;
            case 1: // temp is outside undefined, started at an inside point
                tMaxTemp += deltaT;
                break;
            case 2: // temp is outside defined, started at an outside point
                isFound = true;
                break;
            case 3: // temp is outside defined, started at an inside point
                tMaxTemp += deltaT;
                break;
            case 4: // temp is inside undefined, started at an outside point
                tMaxTemp += deltaT;
                break;
            case 5: // temp is inside undefined, started at an inside point
                tMaxTemp -= deltaT;
                break;
            case 6: // temp is inside defined, started at an outside point
                tMaxTemp += deltaT;
                break;
            case 7: // temp is inside defined, started at an inside point
                isFound = true;
                break;
            }
        }
        if( !isFound ) {
            return false;
        }
    }
    tMaxOut = tMaxTemp;

    return true;
}

bool rle_level_set::intersect_ray( const ray3f& ray, double tMin, double tMax, double& tOut, vector3f& /*outNormal*/ ) {
    // todo: replace this function with a kd-tree raytrace function
    double t0, t1;
    unsigned int numSteps = 100;

    // find t0 and t1, such that they exist in defined regions of the level set on opposite sides of the surface
    if( find_defined_times_along_ray( ray, tMin, tMax, numSteps, t0, t1 ) ) {

        float phi0 = trilerp_signed_distance( ray.at( t0 ) );
        float phi1 = trilerp_signed_distance( ray.at( t1 ) );

        // if both points are on the same side of the surface
        if( ( phi0 > 0 && phi1 > 0 ) || ( phi0 < 0 && phi1 < 0 ) ) {
            return false;
        }

        vector3f outPosition;

        // find the surface location using Neubauer's Method (repeated linear interpolation)
        // apparently 2-3 is a good number of times to repeat the linear interpolation solve,
        // according to "Fast and Acurate Ray-Voxel Intersection Techniques for Iso-Surface Ray Tracing" by Marmitt et
        // al. They also claim that this method is accurate provided there is only one surface crossing along the ray,
        // which is about all we can guarantee without a kd-tree raytrace function (see the "todo" above).
        for( unsigned i = 0; i < 3; ++i ) {
            tOut = t0 + ( t1 - t0 ) * ( -phi0 / ( phi1 - phi0 ) );

            outPosition = ray.at( tOut );
            float outPhi = trilerp_signed_distance( outPosition );

            if( ( outPhi > 0 && phi0 > 0 ) || ( outPhi < 0 && phi0 < 0 ) ) {
                t0 = tOut;
                phi0 = outPhi;
            } else {
                t1 = tOut;
                phi1 = outPhi;
            }
        }

        tOut = t0 + ( t1 - t0 ) * ( -phi0 / ( phi1 - phi0 ) );
        outPosition = ray.at( tOut );

        return true;
    } else {
        return false;
    }
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
