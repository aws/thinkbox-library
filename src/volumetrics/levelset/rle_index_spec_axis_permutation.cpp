// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::graphics;

namespace frantic {
namespace volumetrics {
namespace levelset {

namespace detail {

/**
 * This struct represents one entry in a linked list of run starts for the axis permutation algorithm which
 * switches the compression axis.
 */
struct bfap_run_start_ll {
    // dataIndex     Is either the undefined region code, or the input data index of the first voxel in this run.
    // targetX       Is the X coordinate of this run start in the target voxel coordinates.
    // nextRunStart  Is a link to the next run start in this target scanline.
    int dataIndex, targetX, nextRunStart;
};

struct bfap_run_starts_for_plane {
    /**
     * This is a 1D array, with one value for each target scanline, which contains the head of
     * the linked list of run starts for that target scanline.  All its values start as -1,
     * and if a value is non-negative, it is an index into the runStartLinkedLists array for this plane.
     * The dimension of this array is the same as the output compression axis, because it is for planes
     * which contains both the input and the ouput compression axes.
     */
    vector<int32_t> linkedListHeadArray;
    /**
     * This is a pointer to an array of integers, which should at the end be the number of runs
     * in each scanline.  It is used to do a prefix sum so that we know the layout of the bcToRunIndex
     * array before we start filling in the final RLE Index Spec.
     */
    int32_t* runStartCounts;
    /**
     * This is a pointer to an array of integers, which should at the end be the number of defined voxels
     * in each scanline.  It is used to do a prefix sum so that we know the layout of the data array
     * before we start filling in the final RLE Index Spec.
     */
    int32_t* definedVoxelsCounts;
    /**
     * This array contains all the data for the linked lists of the run starts for target scanlines.
     * We store it in a flat array like this in order to minimize heap allocations, and thus improve
     * performance.
     */
    vector<bfap_run_start_ll> runStartLinkedLists;

    void init( size_t size, int32_t* runStartCountArray_, int32_t* definedVoxelsCountArray_ ) {
        linkedListHeadArray.resize( size, -1 );
        runStartCounts = runStartCountArray_;
        definedVoxelsCounts = definedVoxelsCountArray_;
    }
};

} // namespace detail

// This applies a permutation to the voxel coordinates of the input rle_index_spec, producing a copy.
// It can be used to swap the X and Z coordinates, for instance, to recompress the structure on a different axis.
void rle_index_spec::build_from_axis_permutation(
    // void rle_index_spec::build_from_axis_permutation_unthreaded(
    const vector3& axisPermutation, const rle_index_spec& ris,
    const rle_index_spec_channel_copying_data* channelsToRemapBegin,
    const rle_index_spec_channel_copying_data* channelsToRemapEnd ) {
    if( !( ( axisPermutation.x == 0 || axisPermutation.y == 0 || axisPermutation.z == 0 ) &&
           ( axisPermutation.x == 1 || axisPermutation.y == 1 || axisPermutation.z == 1 ) &&
           ( axisPermutation.x == 2 || axisPermutation.y == 2 || axisPermutation.z == 2 ) ) )
        throw std::runtime_error( "rle_index_spec.build_from_axis_permutation: The permutation provided, " +
                                  axisPermutation.str() + ", is not valid as a permutation." );

    // If the bounding box of the rle_index_spec is empty, then just clear ourselves because there's nothing to do.
    if( ris.m_abcCoordSize.volume() == 0 ) {
        clear();
        m_exteriorRegionCode = ris.m_exteriorRegionCode;
        return;
    }

    // TODO: Should we turn off checks like this once we're confident everything works great?
    /*
    std::stringstream msg;
    if( !ris.check_consistency(msg) ) {
      //ofstream fout("c:\\temp\\blah.txt");
      //ris.dump(fout);
      throw std::runtime_error( "rle_index_spec.build_from_axis_permutation: The rle_index_spec provided is not
    valid.\n"
    + msg.str() );
    }
    */

    /*
    static int seq = 0;
    ofstream fout(files::replace_sequence_number("c:\\temp\\blah.txt",seq).c_str());
    ++seq;
    */

    // Because we're modifying this rle_index_spec, we need to reset the cached ris_adjacency
    free_cached_adjacency();

    m_exteriorRegionCode = ris.m_exteriorRegionCode;

    if( axisPermutation.x == 0 ) {
        // If the compression axis doesn't change, it's simpler, and the permutation is either the identity or (y z)
        if( axisPermutation.y == 1 ) {
            ////// If it's the identity permutation, it's a straight copy //////
            *this = ris;

            // Do a straight copy of the channels to remap
            for( const rle_index_spec_channel_copying_data* i = channelsToRemapBegin; i != channelsToRemapEnd; ++i ) {
                memcpy( i->outputData, i->inputData, i->primitiveSize * data_size() );
            }

        } else {
            ////// This permutation swaps the y and z axes. //////
            m_abcCoordOrigin = vector3( ris.m_abcCoordOrigin.x, ris.m_abcCoordOrigin.z, ris.m_abcCoordOrigin.y );
            m_abcCoordSize =
                size3( ris.m_abcCoordSize.xsize(), ris.m_abcCoordSize.zsize(), ris.m_abcCoordSize.ysize() );
            m_dataSize = ris.m_dataSize;

            // Reserve the correct amount of space
            m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
            m_runData.clear();
            m_runData.reserve( ris.m_runData.size() );

            int currentDataIndex = 0;
            // Iterate through the voxel coordinates and build up the arrays.
            for( int c = 0; c < m_abcCoordSize.zsize(); ++c ) {
                for( int b = 0; b < m_abcCoordSize.ysize(); ++b ) {
                    int bcIndex = b + c * m_abcCoordSize.ysize();
                    m_bcToRunIndex[bcIndex] = (int)m_runData.size();

                    int inputbcIndex = c + b * ris.m_abcCoordSize.ysize();
                    int inputRunRangeStart = ris.m_bcToRunIndex[inputbcIndex],
                        inputRunRangeEnd = ris.m_bcToRunIndex[inputbcIndex + 1] - 1;
                    if( ris.m_runData[inputRunRangeStart].x == ris.m_runData[inputRunRangeStart + 1].x ) {
                        // Create a zero-sized run, necessary for the coordinate queries to work properly
                        m_runData.push_back( run_data( 0, -1 ) );
                        m_runData.push_back( run_data( 0, -1 ) );
                    } else {
                        // Copy the scanline, adjusting the data indices for the new order
                        for( int inputRun = inputRunRangeStart; inputRun != inputRunRangeEnd; ++inputRun ) {
                            if( ris.m_runData[inputRun].dataIndex < 0 ) {
                                // For undefined runs, make an exact copy of the run index data
                                m_runData.push_back( ris.m_runData[inputRun] );
                            } else {
                                // For defined runs, we have to remap the data indices, and copy chunks of data in the
                                // channels being remapped
                                m_runData.push_back( run_data( ris.m_runData[inputRun].x, currentDataIndex ) );

                                int currentRunLength = ris.m_runData[inputRun + 1].x - ris.m_runData[inputRun].x;
                                int sourceDataIndex = ris.m_runData[inputRun].dataIndex;

                                // Do a straight copy of the data for this defined scanline
                                for( const rle_index_spec_channel_copying_data* i = channelsToRemapBegin;
                                     i != channelsToRemapEnd; ++i ) {
                                    memcpy( i->outputData + i->primitiveSize * currentDataIndex,
                                            i->inputData + i->primitiveSize * sourceDataIndex,
                                            i->primitiveSize * currentRunLength );
                                }

                                currentDataIndex += currentRunLength;
                            }
                        }
                        // Copy the data value at the end of the scanline, all that's used is its X value
                        m_runData.push_back( ris.m_runData[inputRunRangeEnd] );
                    }
                }
            }
            // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range
            // of the last scanline
            m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();
        }
    } else {
        // Compute the inverse permutation for the value lookups
        vector3 inverseAxisPermutation;
        inverseAxisPermutation[axisPermutation.x] = 0;
        inverseAxisPermutation[axisPermutation.y] = 1;
        inverseAxisPermutation[axisPermutation.z] = 2;

        //		cout << "axis permutation: " << axisPermutation << endl;
        //		cout << "inv axis permutation: " << inverseAxisPermutation << endl;

        /////////
        // Initialize the target rle_index_spec with some basic values
        /////////

        // The compression axis changes, so we have to recompress along another axis
        m_abcCoordOrigin = vector3( ris.m_abcCoordOrigin[axisPermutation.x], ris.m_abcCoordOrigin[axisPermutation.y],
                                    ris.m_abcCoordOrigin[axisPermutation.z] );
        m_abcCoordSize = size3( ris.m_abcCoordSize[axisPermutation.x], ris.m_abcCoordSize[axisPermutation.y],
                                ris.m_abcCoordSize[axisPermutation.z] );
        m_dataSize = ris.m_dataSize;

        // Reserve the correct amount of space, and clear the arrays
        m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
        m_runData.clear();

        /////////
        // The following data structures are used to compute the layout of the data array (m_bcToRunIndex is
        // computed similarly) so that each scanline can be processed independently.
        /////////

        // This is temporary, it will directly use m_bcToRunIndex in the final implementation
        vector<int32_t> bcToRunIndexTemp( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
        // This is a mapping from bc index to the first data index in its scanline
        vector<int32_t> bcToDataIndex( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );

        /////////
        // The following data structures are used to compute the run starts and data index adjacency needed
        // to efficiently traverse the runs in the target compression axis direction.
        /////////

        // This is the count of run starts for each target scanline, which will be run through the prefix-sum algorithm
        // to fill m_bcToRunIndex in the middle phase of the algorithm.  Initialize the counts with the value 1, for
        // the terminator run start entry each scanline must have.
        vector<int32_t> runStartCounts( m_abcCoordSize.ysize() * m_abcCoordSize.zsize(), 1 );

        // This is the count of defined voxels for each target scanline, which will be run through the prefix sum
        // algorithm to determine the per-scanline layout of the defined voxels.
        vector<int32_t> definedVoxelCounts( m_abcCoordSize.ysize() * m_abcCoordSize.zsize(), 0 );

        // This is an array which maps an input data index to the adjacent data index in the positive direction
        // along the target compression axis.  When processing a defined run, this lets us efficiently visit all
        // the data along that run as we build up the output.
        vector<int32_t> nextInputDataIndex( ris.data_size() );

        // This is the number of XY or XZ planes there are
        int compressionAxisPlanesCount;
        if( axisPermutation.x == 1 )
            compressionAxisPlanesCount = ris.m_abcCoordSize.zsize(); // XY planes
        else
            compressionAxisPlanesCount = ris.m_abcCoordSize.ysize(); // XZ planes

        // This holds one plane of run start linked list and run start counts data for each Z coordinate.
        vector<detail::bfap_run_starts_for_plane> runStartsForPlane( compressionAxisPlanesCount );
        for( size_t i = 0, ie = runStartsForPlane.size(); i != ie; ++i )
            runStartsForPlane[i].init( ris.m_abcCoordSize.xsize(), &runStartCounts[i * ris.m_abcCoordSize.xsize()],
                                       &definedVoxelCounts[i * ris.m_abcCoordSize.xsize()] );

        if( axisPermutation.x == 1 ) {
            // In this case, the permutation is either (0 1) which corresponds to a vector3(1,0,2), or
            // it's (0 1 2) which corresponds to a vector3(1,2,0).  In either case, the input compression
            // axis is X and the output compression axis is Y, so Z is the free axis along which the first
            // step can be multi-threaded.

            /////////
            // INITIAL PHASE: Compute these data structures by iterating over pairs of scanlines in the input
            // compression axis direction.
            /////////

            int exteriorRegionCode = ris.m_exteriorRegionCode;

            // Deal with the boundary case for the nextInputDataIndex array
            for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {

                // Get the run starts data for this plane
                detail::bfap_run_starts_for_plane& rsfp = runStartsForPlane[c];

                // the trailing Y's
                int bcIndex = ( c + 1 ) * ris.m_abcCoordSize.ysize() - 1;
                int runRangeStart = ris.m_bcToRunIndex[bcIndex];
                int runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x, a = rd->x - ris.m_abcCoordOrigin.x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex, ++a ) {
                            nextInputDataIndex[dataIndex] = -1;
                            // Count this defined voxel for the scanline
                            ++rsfp.definedVoxelsCounts[a];
                        }
                    }
                }
            }
            // Deal with the interior, going from maximum B to minimum B, so that the runStartLinkedLists end up
            // in the right order for efficient forward traversal later
            int bsize = ris.m_abcCoordSize.ysize();
            for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {

                // Get the run starts data for this plane
                detail::bfap_run_starts_for_plane& rsfp = runStartsForPlane[c];

                // At each value of b, we process the two scanlines at [*,b-1,c] and [*,b,c] together
                for( int b = bsize - 1; b != 0; --b ) {
                    int bcIndexA = b - 1 + c * bsize;
                    int bcIndexB = bcIndexA + 1;
                    int runRangeStartA = ris.m_bcToRunIndex[bcIndexA],
                        runRangeEndA = ris.m_bcToRunIndex[bcIndexA + 1] - 1;
                    int runRangeStartB = runRangeEndA + 1, runRangeEndB = ris.m_bcToRunIndex[bcIndexB + 1] - 1;

                    int dataIndexA = exteriorRegionCode, dataIndexB = exteriorRegionCode;

                    const run_data* rdA = &ris.m_runData[runRangeStartA];
                    const run_data* rdB = &ris.m_runData[runRangeStartB];
                    int runA = runRangeStartA, runB = runRangeStartB;
                    int xA = rdA->x, xB = rdB->x;
                    // Deal with empty scanlines
                    if( xA == ( rdA + 1 )->x ) {
                        runA = runRangeEndA;
                        xA = numeric_limits<int>::max();
                    }
                    if( xB == ( rdB + 1 )->x ) {
                        runB = runRangeEndB;
                        xB = numeric_limits<int>::max();
                    }
                    do {
                        int x;
                        // Increment to the next sub-run
                        if( xA < xB ) {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                        } else if( xA > xB ) {
                            x = xB;
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        } else {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        }

                        // Figure out the size of this sub-run
                        int size = min( xA, xB ) - x;

                        //
                        if( dataIndexA >= 0 ) {
                            if( dataIndexB >= 0 ) {
                                // In this case we have two defined runs
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    nextInputDataIndex[dataIndexA] = dataIndexB;
                                    // Count this defined voxel for the scanline
                                    ++rsfp.definedVoxelsCounts[a];
                                    ++dataIndexA;
                                    ++dataIndexB;
                                }
                            } else {
                                // in this case just A is defined
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    nextInputDataIndex[dataIndexA] = -1;
                                    // Count this defined voxel for the scanline
                                    ++rsfp.definedVoxelsCounts[a];
                                    ++dataIndexA;
                                }
                                // there is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.y + b;
                                a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    //								cout << i << "/" << size << ", coord " << a << "/" <<
                                    //ris.m_abcCoordSize.xsize() <<
                                    //", " << c << ", in size " << rsfp.linkedListHeadArray.size() << endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                                    rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                                    rsfp.runStartLinkedLists.push_back( rsl );
                                    // If this isn't the very last run in the scanline, or if it is, and the region code
                                    // doesn't match the exterior region code, then we increment the run start count for
                                    // this scanline
                                    if( rsl.nextRunStart != -1 || rsl.dataIndex != exteriorRegionCode ) {
                                        ++rsfp.runStartCounts[a];
                                    }
                                }
                            }
                        } else {
                            if( dataIndexB >= 0 ) {
                                // In this case just B is defined, so there is the start
                                // of a defined run for each voxel in this sub-run along scanline B.
                                detail::bfap_run_start_ll rsl;
                                rsl.targetX = ris.m_abcCoordOrigin.y + b;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    //									cout << "coord " << a << ", " << c << ", in size " <<
                                    //linkedListHeadArray.size()
                                    //<< endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.dataIndex = dataIndexB;
                                    ++dataIndexB;
                                    rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                                    rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                                    rsfp.runStartLinkedLists.push_back( rsl );
                                    // Since this run is defined, we always increment the run start count for this
                                    // scanline
                                    ++rsfp.runStartCounts[a];
                                }
                            } else if( dataIndexA != dataIndexB ) {
                                // In this case both are undefined, but the region code changed.
                                // There is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.y + b;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    //									cout << "coord " << a << ", " << c << ", in size " <<
                                    //linkedListHeadArray.size()
                                    //<< endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                                    rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                                    rsfp.runStartLinkedLists.push_back( rsl );
                                    // If this isn't the very last run in the scanline, or if it is, and the region code
                                    // doesn't match the exterior region code, then we increment the run start count for
                                    // this scanline
                                    if( rsl.nextRunStart != -1 || rsl.dataIndex != exteriorRegionCode ) {
                                        ++rsfp.runStartCounts[a];
                                    }
                                }
                            }
                        }
                    } while( xB != numeric_limits<int>::max() || xA != numeric_limits<int>::max() );
                }
            }
            // Deal with the boundary case for the linkedListHeadArray array
            for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {

                // Get the run starts data for this plane
                detail::bfap_run_starts_for_plane& rsfp = runStartsForPlane[c];

                // the leading Y's
                int bcIndex = 0 + c * ris.m_abcCoordSize.ysize();
                int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                detail::bfap_run_start_ll rsl;
                rsl.targetX = ris.m_abcCoordOrigin.y + 0;
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                        int a = rd->x - ris.m_abcCoordOrigin.x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex, ++a ) {
                            // Hook this run start into the front of the scanline linked list
                            rsl.dataIndex = dataIndex;
                            rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                            rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                            rsfp.runStartLinkedLists.push_back( rsl );
                            // Since this run is defined, we always increment the run start count for this scanline
                            ++rsfp.runStartCounts[a];
                        }
                    } else if( dataIndex != exteriorRegionCode ) {
                        int a = rd->x - ris.m_abcCoordOrigin.x, aEnd = a + ( rd + 1 )->x - rd->x;
                        rsl.dataIndex = dataIndex;
                        for( ; a != aEnd; ++a ) {
                            // Hook this run start into the front of the scanline linked list
                            rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                            rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                            rsfp.runStartLinkedLists.push_back( rsl );
                            // Since this run is not of the exterior region code, we always increment the run start
                            // count for this scanline
                            ++rsfp.runStartCounts[a];
                        }
                    }
                }
            }

            /////////
            // PREFIX SUM PHASE: Compute the per-scanline layout of the run index data and the defined voxels.
            /////////

            // First fix up the runStartCounts array so that any empty scanlines have count 2 instead of 1
            for( vector<int32_t>::iterator i = runStartCounts.begin(), ie = runStartCounts.end(); i != ie; ++i ) {
                if( *i == 1 )
                    *i = 2;
            }

            // Compute the prefix sum to fill bcToRunIndex, and the prefix sum to get the per-scanline data index layout
            bcToRunIndexTemp[0] = 0;
            bcToDataIndex[0] = 0;
            int32_t accumRun = 0, accumData = 0;
            // Depending on whether Z maps to X, or Z maps to Z, the indexing order of runStartCounts is different
            if( axisPermutation.z == 0 ) {
                size_t i = 0;
                for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
                    for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b, ++i ) {
                        int j = c + csize * b;
                        accumRun += runStartCounts[j];
                        bcToRunIndexTemp[i + 1] = accumRun;
                        accumData += definedVoxelCounts[j];
                        bcToDataIndex[i + 1] = accumData;
                    }
                }
            } else {
                for( size_t i = 0, ie = runStartCounts.size(); i != ie; ++i ) {
                    accumRun += runStartCounts[i];
                    bcToRunIndexTemp[i + 1] = accumRun;
                    accumData += definedVoxelCounts[i];
                    bcToDataIndex[i + 1] = accumData;
                }
            }

        } else { // axisPermutation.x == 2
            // In this case, the permuation is either (0 2) which corresponds to a vector3(2,1,0), or
            // it's (0 2 1) which corresponds to a vector3(2,0,1).  In either case, the input compression
            // axis is X and the output compression axis is Z, so Y is the free axis along which the first
            // step can be multi-threaded.

            /////////
            // INITIAL PHASE: Compute these data structures by iterating over pairs of scanlines in the input
            // compression axis direction.
            /////////

            int exteriorRegionCode = ris.m_exteriorRegionCode;

            // Deal with the boundary case for the nextInputDataIndex array
            for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {

                // Get the run starts data for this plane
                detail::bfap_run_starts_for_plane& rsfp = runStartsForPlane[b];

                // the trailing Z's
                int bcIndex = b + ( ris.m_abcCoordSize.zsize() - 1 ) * ris.m_abcCoordSize.ysize();
                int runRangeStart = ris.m_bcToRunIndex[bcIndex];
                int runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x, a = rd->x - ris.m_abcCoordOrigin.x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex, ++a ) {
                            nextInputDataIndex[dataIndex] = -1;
                            // Count this defined voxel for the scanline
                            ++rsfp.definedVoxelsCounts[a];
                        }
                    }
                }
            }
            // Deal with the interior, going from maximum C to minimum C, so that the runStartLinkedLists end up
            // in the right order for efficient forward traversal later
            int csize = ris.m_abcCoordSize.zsize();
            for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {

                // Get the run starts data for this plane
                detail::bfap_run_starts_for_plane& rsfp = runStartsForPlane[b];

                // At each value of c, we process the two scanlines at [*,b,c-1] and [*,b,c] together
                for( int c = csize - 1; c != 0; --c ) {
                    int bcIndexA = b + ( c - 1 ) * bsize;
                    int bcIndexB = bcIndexA + bsize;
                    int runRangeStartA = ris.m_bcToRunIndex[bcIndexA],
                        runRangeEndA = ris.m_bcToRunIndex[bcIndexA + 1] - 1;
                    int runRangeStartB = ris.m_bcToRunIndex[bcIndexB],
                        runRangeEndB = ris.m_bcToRunIndex[bcIndexB + 1] - 1;

                    int dataIndexA = exteriorRegionCode, dataIndexB = exteriorRegionCode;

                    const run_data* rdA = &ris.m_runData[runRangeStartA];
                    const run_data* rdB = &ris.m_runData[runRangeStartB];
                    int runA = runRangeStartA, runB = runRangeStartB;
                    int xA = rdA->x, xB = rdB->x;
                    // Deal with empty scanlines
                    if( xA == ( rdA + 1 )->x ) {
                        runA = runRangeEndA;
                        xA = numeric_limits<int>::max();
                    }
                    if( xB == ( rdB + 1 )->x ) {
                        runB = runRangeEndB;
                        xB = numeric_limits<int>::max();
                    }
                    do {
                        int x;
                        // Increment to the next sub-run
                        if( xA < xB ) {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                        } else if( xA > xB ) {
                            x = xB;
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        } else {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        }

                        // Figure out the size of this sub-run
                        int size = min( xA, xB ) - x;

                        //
                        if( dataIndexA >= 0 ) {
                            if( dataIndexB >= 0 ) {
                                // In this case we have two defined runs
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    nextInputDataIndex[dataIndexA] = dataIndexB;
                                    // Count this defined voxel for the scanline
                                    ++rsfp.definedVoxelsCounts[a];
                                    ++dataIndexA;
                                    ++dataIndexB;
                                }
                            } else {
                                // in this case just A is defined
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    nextInputDataIndex[dataIndexA] = -1;
                                    // Count this defined voxel for the scanline
                                    ++rsfp.definedVoxelsCounts[a];
                                    ++dataIndexA;
                                }
                                // there is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.z + c;
                                a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    //									cout << "coord " << a << ", " << b << ", in size " <<
                                    //linkedListHeadArray.size()
                                    //<< endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                                    rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                                    rsfp.runStartLinkedLists.push_back( rsl );
                                    // If this isn't the very last run in the scanline, or if it is, and the region code
                                    // doesn't match the exterior region code, then we increment the run start count for
                                    // this scanline
                                    if( rsl.nextRunStart != -1 || rsl.dataIndex != exteriorRegionCode ) {
                                        ++rsfp.runStartCounts[a];
                                    }
                                }
                            }
                        } else {
                            if( dataIndexB >= 0 ) {
                                // In this case just B is defined, so there is the start
                                // of a defined run for each voxel in this sub-run along scanline B.
                                detail::bfap_run_start_ll rsl;
                                rsl.targetX = ris.m_abcCoordOrigin.z + c;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    //									cout << "coord " << a << ", " << b << ", in size " <<
                                    //linkedListHeadArray.size()
                                    //<< endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.dataIndex = dataIndexB;
                                    ++dataIndexB;
                                    rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                                    rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                                    rsfp.runStartLinkedLists.push_back( rsl );
                                    // Since this run is defined, we always increment the run start count for this
                                    // scanline
                                    ++rsfp.runStartCounts[a];
                                }
                            } else if( dataIndexA != dataIndexB ) {
                                // In this case both are undefined, but the region code changed.
                                // There is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.z + c;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    //									cout << "coord " << a << ", " << b << ", in size " <<
                                    //linkedListHeadArray.size()
                                    //<< endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                                    rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                                    rsfp.runStartLinkedLists.push_back( rsl );
                                    // If this isn't the very last run in the scanline, or if it is, and the region code
                                    // doesn't match the exterior region code, then we increment the run start count for
                                    // this scanline
                                    if( rsl.nextRunStart != -1 || rsl.dataIndex != exteriorRegionCode ) {
                                        ++rsfp.runStartCounts[a];
                                    }
                                }
                            }
                        }
                    } while( xB != numeric_limits<int>::max() || xA != numeric_limits<int>::max() );
                }
            }
            // Deal with the boundary case for the linkedListHeadArray array
            for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {

                // Get the run starts data for this plane
                detail::bfap_run_starts_for_plane& rsfp = runStartsForPlane[b];

                // the leading Z's
                int bcIndex = b;
                int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                detail::bfap_run_start_ll rsl;
                rsl.targetX = ris.m_abcCoordOrigin.z + 0;
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                        int a = rd->x - ris.m_abcCoordOrigin.x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex, ++a ) {
                            // Hook this run start into the front of the scanline linked list
                            rsl.dataIndex = dataIndex;
                            rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                            rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                            rsfp.runStartLinkedLists.push_back( rsl );
                            // Since this run is defined, we always increment the run start count for this scanline
                            ++rsfp.runStartCounts[a];
                        }
                    } else if( dataIndex != exteriorRegionCode ) {
                        int a = rd->x - ris.m_abcCoordOrigin.x, aEnd = a + ( rd + 1 )->x - rd->x;
                        rsl.dataIndex = dataIndex;
                        for( ; a != aEnd; ++a ) {
                            // Hook this run start into the front of the scanline linked list
                            rsl.nextRunStart = rsfp.linkedListHeadArray[a];
                            rsfp.linkedListHeadArray[a] = (int)rsfp.runStartLinkedLists.size();
                            rsfp.runStartLinkedLists.push_back( rsl );
                            // Since this run is not of the exterior region code, we always increment the run start
                            // count for this scanline
                            ++rsfp.runStartCounts[a];
                        }
                    }
                }
            }

            /////////
            // PREFIX SUM PHASE: Compute the per-scanline layout of the run index data and the defined voxels.
            /////////

            // First fix up the runStartCounts array so that any empty scanlines have count 2 instead of 1
            for( vector<int32_t>::iterator i = runStartCounts.begin(), ie = runStartCounts.end(); i != ie; ++i ) {
                if( *i == 1 )
                    *i = 2;
            }

            // Compute the prefix sum to fill bcToRunIndex, and the prefix sum to get the per-scanline data index layout
            bcToRunIndexTemp[0] = 0;
            bcToDataIndex[0] = 0;
            int32_t accumRun = 0, accumData = 0;
            // Depending on whether Z maps to X, or Z maps to Z, the indexing order of runStartCounts is different
            if( axisPermutation.z == 0 ) {
                size_t i = 0;
                for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
                    for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b, ++i ) {
                        int j = c + csize * b;
                        accumRun += runStartCounts[j];
                        bcToRunIndexTemp[i + 1] = accumRun;
                        accumData += definedVoxelCounts[j];
                        bcToDataIndex[i + 1] = accumData;
                    }
                }
            } else {
                for( size_t i = 0, ie = runStartCounts.size(); i != ie; ++i ) {
                    accumRun += runStartCounts[i];
                    bcToRunIndexTemp[i + 1] = accumRun;
                    accumData += definedVoxelCounts[i];
                    bcToDataIndex[i + 1] = accumData;
                }
            }
        }

        /////////
        // FINAL PHASE: Now use the data structures created to build up the runs and copy the data for the target.
        /////////

        //		fout << "Perm: " << axisPermutation << "\n";
        //		fout << "Found " << runStartLinkedLists.size() << " run starts\n";
        //		fout << "Outer bounds: " << outer_bounds() << "\n";

        int nextOutputDataIndex = 0;
        vector3 xyz;
        // Iterate through the voxel coordinates and build up the arrays.
        for( int c = 0; c < m_abcCoordSize.zsize(); ++c ) {
            for( int b = 0; b < m_abcCoordSize.ysize(); ++b ) {
                xyz.y = b + m_abcCoordOrigin.y;
                xyz.z = c + m_abcCoordOrigin.z;

                //				fout << "BC coordinate (" << b << ", " << c << ")\n";

                int bcIndex = b + c * m_abcCoordSize.ysize();
                m_bcToRunIndex[bcIndex] = (int)m_runData.size();

                // The linkedListHeadArray was indexed based on the input, not the target, so we have to check how to
                // get its index here.
                int runStartLLIndex;

                detail::bfap_run_starts_for_plane* rsfp = 0;

                // Test whether to swich the role of b and c in the linkedListHeadArray
                if( axisPermutation.z == 0 ) {
                    rsfp = &runStartsForPlane[b];
                    runStartLLIndex = rsfp->linkedListHeadArray[c];
                } else {
                    rsfp = &runStartsForPlane[c];
                    runStartLLIndex = rsfp->linkedListHeadArray[b];
                }

                // The default state is an undefined region of the exterior region code
                int creatingRegionCode = m_exteriorRegionCode;
                while( runStartLLIndex != -1 ) {
                    const detail::bfap_run_start_ll& rsll = rsfp->runStartLinkedLists[runStartLLIndex];
                    int inputDataIndex = rsll.dataIndex;

                    //					fout << "rsll.targetX: " << rsll.targetX << "\n";
                    //					fout << "rsll.dataIndex: " << rsll.dataIndex << "\n";
                    //					fout << "rsll.nextRunStart: " << rsll.nextRunStart << "\n";

                    if( inputDataIndex != creatingRegionCode ) {
                        creatingRegionCode = inputDataIndex;
                        if( inputDataIndex < 0 ) {
                            m_runData.push_back( run_data( rsll.targetX, inputDataIndex ) );
                        } else {
                            m_runData.push_back( run_data( rsll.targetX, nextOutputDataIndex ) );

                            do {
                                // Copy this value for each of the channels being remapped
                                for( const rle_index_spec_channel_copying_data* i = channelsToRemapBegin;
                                     i != channelsToRemapEnd; ++i ) {
                                    //							cout << "channel " << i << ", inputDataIndex " << inputDataIndex <<
                                    //", nextOutputDataIndex " << nextOutputDataIndex << endl;
                                    // cout << "channelsToRemap[i].outputData + channelsToRemap[i].primitiveSize *
                                    // nextOutputDataIndex " << (void*)(channelsToRemap[i].outputData +
                                    // channelsToRemap[i].primitiveSize * nextOutputDataIndex) << endl; 							cout
                                    // << "channelsToRemap[i].inputData + channelsToRemap[i].primitiveSize *
                                    //inputDataIndex " << (void*)(channelsToRemap[i].inputData +
                                    //channelsToRemap[i].primitiveSize * inputDataIndex) << endl; 							cout <<
                                    //"channelsToRemap[i].primitiveSize " << channelsToRemap[i].primitiveSize << endl;
                                    // cout << "isbadreadptr: " << IsBadReadPtr(channelsToRemap[i].outputData +
                                    // channelsToRemap[i].primitiveSize * nextOutputDataIndex,
                                    // channelsToRemap[i].primitiveSize) << endl; 							cout << "isbadreadptr: " <<
                                    //IsBadReadPtr(channelsToRemap[i].inputData + channelsToRemap[i].primitiveSize *
                                    // inputDataIndex, channelsToRemap[i].primitiveSize) << endl;
                                    memcpy( i->outputData + i->primitiveSize * nextOutputDataIndex,
                                            i->inputData + i->primitiveSize * inputDataIndex, i->primitiveSize );
                                }
                                ++nextOutputDataIndex;
                                inputDataIndex = nextInputDataIndex[inputDataIndex];
                            } while( inputDataIndex >= 0 );
                        }
                    }

                    runStartLLIndex = rsll.nextRunStart;
                }
                // If any runs were flagged, finish off the last run.  Otherwise add a zero-sized run
                if( (int)m_runData.size() > m_bcToRunIndex[bcIndex] ) {
                    // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to
                    // optimize the storage and access.
                    if( creatingRegionCode != m_exteriorRegionCode )
                        m_runData.push_back( run_data( m_abcCoordOrigin.x + m_abcCoordSize.xsize(), -1 ) );
                } else {
                    m_runData.push_back( run_data( 0, -1 ) );
                    m_runData.push_back( run_data( 0, -1 ) );
                }
            }
        }
        // Finish off the m_bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of
        // the last scanline
        m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();

        // TEMP DEBUG CODE: The bcToRunIndexTemp should match m_bcToRunIndexTemp
        for( size_t i = 0, ie = bcToRunIndexTemp.size(); i != ie; ++i ) {
            if( bcToRunIndexTemp[i] != m_bcToRunIndex[i] ) {
                std::ofstream fout( "c:\\debug.txt" );
                ris.dump( fout );
                for( size_t j = 0, je = bcToRunIndexTemp.size(); j != je; ++j ) {
                    fout << m_bcToRunIndex[j] << "\t" << bcToRunIndexTemp[j] << "\n";
                }
                throw runtime_error( "axis permutation (Z->X): bc to run index " + lexical_cast<string>( i ) +
                                     " doesn't match." );
            }
        }
        // TEMP DEBUG CODE: The bcToDataIndex should match up with the info in m_runData
        for( size_t i = 0, ie = bcToDataIndex.size() - 1; i != ie; ++i ) {
            // Find the first defined run
            int runRangeStart = m_bcToRunIndex[i], runRangeEnd = m_bcToRunIndex[i + 1] - 1;
            for( int run = runRangeStart; run != runRangeEnd; ++run ) {
                if( m_runData[run].dataIndex >= 0 ) {
                    // In the first defined run, check that bcToDataIndex[i] matches
                    if( bcToDataIndex[i] != m_runData[run].dataIndex ) {
                        std::ofstream fout( "c:\\debug.txt" );
                        this->dump( fout );
                        for( size_t j = 0, je = bcToDataIndex.size() - 1; j != je; ++j ) {
                            int runRangeStart = m_bcToRunIndex[j], runRangeEnd = m_bcToRunIndex[j + 1] - 1;
                            int run = runRangeStart;
                            while( run < runRangeEnd - 1 && m_runData[m_bcToRunIndex[run]].dataIndex < 0 )
                                ++run;
                            // This dump isn't quite the first defined data index, but we can make it more precise if
                            // necessary
                            fout << m_runData[run].dataIndex << "\t" << bcToDataIndex[j] << "\n";
                        }
                        fout << endl << "run index range: " << runRangeStart << " " << runRangeEnd << endl;
                        fout << "run: " << run << endl;
                        fout << endl << "m_runData[run].dataIndex: " << m_runData[run].dataIndex << endl;
                        fout << "bcToDataIndex[i]: " << bcToDataIndex[i] << endl;
                        throw runtime_error( "axis permutation (Z->X): bc to data index " + lexical_cast<string>( i ) +
                                             " doesn't match." );
                    }
                    break;
                }
            }
        }
    }

    /*
    std::stringstream msg2;
    if( !check_consistency(msg2) ) {
      fout << "\n\nERROR!!!: " << msg2.str() << "\n\n";
      ris.dump(fout);
      fout << "\n\n\n\n\n\n";
      dump(fout);
      throw std::runtime_error( "rle_index_spec.build_from_axis_permutation: The rle_index_spec produced is not
    valid.\n"
    + msg2.str() );
    }
    */
}

/////////////////////////////////////////////////////////////////////////
// The following is to be deleted:
/////////////////////////////////////////////////////////////////////////

// This applies a permutation to the voxel coordinates of the input rle_index_spec, producing a copy.
// It can be used to swap the X and Z coordinates, for instance, to recompress the structure on a different axis.
void rle_index_spec::build_from_axis_permutation_unthreaded(
    // void rle_index_spec::build_from_axis_permutation(
    const vector3& axisPermutation, const rle_index_spec& ris,
    const rle_index_spec_channel_copying_data* channelsToRemapBegin,
    const rle_index_spec_channel_copying_data* channelsToRemapEnd ) {
    if( !( ( axisPermutation.x == 0 || axisPermutation.y == 0 || axisPermutation.z == 0 ) &&
           ( axisPermutation.x == 1 || axisPermutation.y == 1 || axisPermutation.z == 1 ) &&
           ( axisPermutation.x == 2 || axisPermutation.y == 2 || axisPermutation.z == 2 ) ) )
        throw std::runtime_error( "rle_index_spec.build_from_axis_permutation: The permutation provided, " +
                                  axisPermutation.str() + ", is not valid as a permutation." );

    // TODO: Should we turn off checks like this once we're confident everything works great?
    /*
    std::stringstream msg;
    if( !ris.check_consistency(msg) ) {
      //ofstream fout("c:\\temp\\blah.txt");
      //ris.dump(fout);
      throw std::runtime_error( "rle_index_spec.build_from_axis_permutation: The rle_index_spec provided is not
    valid.\n"
    + msg.str() );
    }
    */

    /*
    static int seq = 0;
    ofstream fout(files::replace_sequence_number("c:\\temp\\blah.txt",seq).c_str());
    ++seq;
    */

    // Because we're modifying this rle_index_spec, we need to reset the cached ris_adjacency
    free_cached_adjacency();

    m_exteriorRegionCode = ris.m_exteriorRegionCode;

    if( axisPermutation.x == 0 ) {
        // If the compression axis doesn't change, it's simpler, and the permutation is either the identity or (y z)
        if( axisPermutation.y == 1 ) {
            ////// If it's the identity permutation, it's a straight copy //////
            *this = ris;

            // Do a straight copy of the channels to remap
            for( const rle_index_spec_channel_copying_data* i = channelsToRemapBegin; i != channelsToRemapEnd; ++i ) {
                memcpy( i->outputData, i->inputData, i->primitiveSize * data_size() );
            }

        } else {
            ////// This permutation swaps the y and z axes. //////
            m_abcCoordOrigin = vector3( ris.m_abcCoordOrigin.x, ris.m_abcCoordOrigin.z, ris.m_abcCoordOrigin.y );
            m_abcCoordSize =
                size3( ris.m_abcCoordSize.xsize(), ris.m_abcCoordSize.zsize(), ris.m_abcCoordSize.ysize() );
            m_dataSize = ris.m_dataSize;

            // Reserve the correct amount of space
            m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
            m_runData.clear();
            m_runData.reserve( ris.m_runData.size() );

            int currentDataIndex = 0;
            // Iterate through the voxel coordinates and build up the arrays.
            for( int c = 0; c < m_abcCoordSize.zsize(); ++c ) {
                for( int b = 0; b < m_abcCoordSize.ysize(); ++b ) {
                    int bcIndex = b + c * m_abcCoordSize.ysize();
                    m_bcToRunIndex[bcIndex] = (int)m_runData.size();

                    int inputbcIndex = c + b * ris.m_abcCoordSize.ysize();
                    int inputRunRangeStart = ris.m_bcToRunIndex[inputbcIndex],
                        inputRunRangeEnd = ris.m_bcToRunIndex[inputbcIndex + 1] - 1;
                    if( ris.m_runData[inputRunRangeStart].x == ris.m_runData[inputRunRangeStart + 1].x ) {
                        // Create a zero-sized run, necessary for the coordinate queries to work properly
                        m_runData.push_back( run_data( 0, -1 ) );
                        m_runData.push_back( run_data( 0, -1 ) );
                    } else {
                        // Copy the scanline, adjusting the data indices for the new order
                        for( int inputRun = inputRunRangeStart; inputRun != inputRunRangeEnd; ++inputRun ) {
                            if( ris.m_runData[inputRun].dataIndex < 0 ) {
                                // For undefined runs, make an exact copy of the run index data
                                m_runData.push_back( ris.m_runData[inputRun] );
                            } else {
                                // For defined runs, we have to remap the data indices, and copy chunks of data in the
                                // channels being remapped
                                m_runData.push_back( run_data( ris.m_runData[inputRun].x, currentDataIndex ) );

                                int currentRunLength = ris.m_runData[inputRun + 1].x - ris.m_runData[inputRun].x;
                                int sourceDataIndex = ris.m_runData[inputRun].dataIndex;

                                // Do a straight copy of the data for this defined scanline
                                for( const rle_index_spec_channel_copying_data* i = channelsToRemapBegin;
                                     i != channelsToRemapEnd; ++i ) {
                                    memcpy( i->outputData + i->primitiveSize * currentDataIndex,
                                            i->inputData + i->primitiveSize * sourceDataIndex,
                                            i->primitiveSize * currentRunLength );
                                }

                                currentDataIndex += currentRunLength;
                            }
                        }
                        // Copy the data value at the end of the scanline, all that's used is its X value
                        m_runData.push_back( ris.m_runData[inputRunRangeEnd] );
                    }
                }
            }
            // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range
            // of the last scanline
            m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();
        }
    } else {
        // Compute the inverse permutation for the value lookups
        vector3 inverseAxisPermutation;
        inverseAxisPermutation[axisPermutation.x] = 0;
        inverseAxisPermutation[axisPermutation.y] = 1;
        inverseAxisPermutation[axisPermutation.z] = 2;

        //		cout << "axis permutation: " << axisPermutation << endl;
        //		cout << "inv axis permutation: " << inverseAxisPermutation << endl;

        /////////
        // Initialize the target rle_index_spec with some basic values
        /////////

        // The compression axis changes, so we have to recompress along another axis
        m_abcCoordOrigin = vector3( ris.m_abcCoordOrigin[axisPermutation.x], ris.m_abcCoordOrigin[axisPermutation.y],
                                    ris.m_abcCoordOrigin[axisPermutation.z] );
        m_abcCoordSize = size3( ris.m_abcCoordSize[axisPermutation.x], ris.m_abcCoordSize[axisPermutation.y],
                                ris.m_abcCoordSize[axisPermutation.z] );
        m_dataSize = ris.m_dataSize;

        // Reserve the correct amount of space, and clear the arrays
        m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
        m_runData.clear();

        /////////
        // The following data structures are used to compute the run starts and data index adjacency needed
        // to efficiently traverse the runs in the target compression axis direction.
        /////////

        // This is a 2D array, with one value for each target scanline, which contains the head of
        // the linked list of run starts for that target scanline.  All its values start as -1,
        // and if a value is non-negative, it is an index into runStartLinkedLists.
        // In the input ABC coordinates, the coordinates into this array will be either [a,b] or [a,c].
        // This array is based on that ordering, with an array of arrays along the A direction.  This may
        // differ from the canonical 2D arrays in the target index spec.
        vector<int> linkedListHeadArray( m_abcCoordSize.ysize() * m_abcCoordSize.zsize(), -1 );
        //		cout << "input abc coord size: " << ris.m_abcCoordSize << endl;
        //		cout << "target abc coord size: " << m_abcCoordSize << endl;
        // This array contains all the data for the linked lists of the run starts for target scanlines.
        // We store it in a flat array like this in order to minimize heap allocations, and thus improve
        // performance.
        vector<detail::bfap_run_start_ll> runStartLinkedLists;
        // Reserve a bit more storage than there are runs, to try and avoid extra memory allocations
        runStartLinkedLists.reserve( ris.m_runData.size() * 12 / 8 );

        // This is an array which maps an input data index to the adjacent data index in the positive direction
        // along the target compression axis.  When processing a defined run, this lets us efficiently visit all
        // the data along that run as we build up the output.
        vector<int> nextInputDataIndex( ris.data_size() );

        /////////
        // Compute these data structures by iterating over pairs of scanlines in the input compression axis direction.
        /////////

        int exteriorRegionCode = ris.m_exteriorRegionCode;
        // If the target compression axis is equivalent to the input Y axis
        if( axisPermutation.x == 1 ) {
            // Deal with the boundary case for the nextInputDataIndex array
            for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {
                // the trailing Y's
                int bcIndex = ( c + 1 ) * ris.m_abcCoordSize.ysize() - 1;
                int runRangeStart = ris.m_bcToRunIndex[bcIndex];
                int runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex ) {
                            nextInputDataIndex[dataIndex] = -1;
                        }
                    }
                }
            }
            // Deal with the interior, going from maximum B to minimum B, so that the runStartLinkedLists end up
            // in the right order for efficient forward traversal later
            int bsize = ris.m_abcCoordSize.ysize();
            for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {
                // At each value of b, we process the two scanlines at [*,b-1,c] and [*,b,c] together
                for( int b = bsize - 1; b != 0; --b ) {
                    int bcIndexA = b - 1 + c * bsize;
                    int bcIndexB = bcIndexA + 1;
                    int runRangeStartA = ris.m_bcToRunIndex[bcIndexA],
                        runRangeEndA = ris.m_bcToRunIndex[bcIndexA + 1] - 1;
                    int runRangeStartB = runRangeEndA + 1, runRangeEndB = ris.m_bcToRunIndex[bcIndexB + 1] - 1;

                    int dataIndexA = exteriorRegionCode, dataIndexB = exteriorRegionCode;

                    const run_data* rdA = &ris.m_runData[runRangeStartA];
                    const run_data* rdB = &ris.m_runData[runRangeStartB];
                    int runA = runRangeStartA, runB = runRangeStartB;
                    int xA = rdA->x, xB = rdB->x;
                    // Deal with empty scanlines
                    if( xA == ( rdA + 1 )->x ) {
                        runA = runRangeEndA;
                        xA = numeric_limits<int>::max();
                    }
                    if( xB == ( rdB + 1 )->x ) {
                        runB = runRangeEndB;
                        xB = numeric_limits<int>::max();
                    }
                    do {
                        int x;
                        // Increment to the next sub-run
                        if( xA < xB ) {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                        } else if( xA > xB ) {
                            x = xB;
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        } else {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        }

                        // Figure out the size of this sub-run
                        int size = min( xA, xB ) - x;

                        //
                        if( dataIndexA >= 0 ) {
                            if( dataIndexB >= 0 ) {
                                // In this case we have two defined runs
                                for( int i = 0; i != size; ++i ) {
                                    nextInputDataIndex[dataIndexA] = dataIndexB;
                                    ++dataIndexA;
                                    ++dataIndexB;
                                }
                            } else {
                                // in this case just A is defined
                                for( int i = 0; i != size; ++i ) {
                                    nextInputDataIndex[dataIndexA] = -1;
                                    ++dataIndexA;
                                }
                                // there is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.y + b;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    // Get the index in the linked list head array
                                    int index = a + c * ris.m_abcCoordSize.xsize();
                                    //									cout << "index " << index << ", coord " << a << ", " << c << ",
                                    //in size " << linkedListHeadArray.size() << endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = linkedListHeadArray[index];
                                    linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                                    runStartLinkedLists.push_back( rsl );
                                }
                            }
                        } else {
                            if( dataIndexB >= 0 ) {
                                // In this case just B is defined, so there is the start
                                // of a defined run for each voxel in this sub-run along scanline B.
                                detail::bfap_run_start_ll rsl;
                                rsl.targetX = ris.m_abcCoordOrigin.y + b;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    // Get the index in the linked list head array
                                    int index = a + c * ris.m_abcCoordSize.xsize();
                                    //									cout << "index " << index << ", coord " << a << ", " << c << ",
                                    //in size " << linkedListHeadArray.size() << endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.dataIndex = dataIndexB;
                                    ++dataIndexB;
                                    rsl.nextRunStart = linkedListHeadArray[index];
                                    linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                                    runStartLinkedLists.push_back( rsl );
                                }
                            } else if( dataIndexA != dataIndexB ) {
                                // In this case both are undefined, but the region code changed.
                                // There is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.y + b;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    // Get the index in the linked list head array
                                    int index = a + c * ris.m_abcCoordSize.xsize();
                                    //									cout << "index " << index << ", coord " << a << ", " << c << ",
                                    //in size " << linkedListHeadArray.size() << endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = linkedListHeadArray[index];
                                    linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                                    runStartLinkedLists.push_back( rsl );
                                }
                            }
                        }
                    } while( xB != numeric_limits<int>::max() || xA != numeric_limits<int>::max() );
                }
            }
            // Deal with the boundary case for the linkedListHeadArray array
            for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {
                // the leading Y's
                int bcIndex = 0 + c * ris.m_abcCoordSize.ysize();
                int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                detail::bfap_run_start_ll rsl;
                rsl.targetX = ris.m_abcCoordOrigin.y + 0;
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                        int a = rd->x - ris.m_abcCoordOrigin.x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex, ++a ) {
                            // Get the index in the linked list head array
                            int index = a + c * ris.m_abcCoordSize.xsize();
                            // Hook this run start into the front of the scanline linked list
                            rsl.dataIndex = dataIndex;
                            rsl.nextRunStart = linkedListHeadArray[index];
                            linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                            runStartLinkedLists.push_back( rsl );
                        }
                    } else if( dataIndex != exteriorRegionCode ) {
                        int a = rd->x - ris.m_abcCoordOrigin.x, aEnd = a + ( rd + 1 )->x - rd->x;
                        rsl.dataIndex = dataIndex;
                        for( ; a != aEnd; ++a ) {
                            // Get the index in the linked list head array
                            int index = a + c * ris.m_abcCoordSize.xsize();
                            // Hook this run start into the front of the scanline linked list
                            rsl.nextRunStart = linkedListHeadArray[index];
                            linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                            runStartLinkedLists.push_back( rsl );
                        }
                    }
                }
            }
        }
        // Otherwise the target compression axis is equivalent to the input Z axis
        else {
            // Deal with the boundary case for the nextInputDataIndex array
            for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {
                // the trailing Z's
                int bcIndex = b + ( ris.m_abcCoordSize.zsize() - 1 ) * ris.m_abcCoordSize.ysize();
                int runRangeStart = ris.m_bcToRunIndex[bcIndex];
                int runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex ) {
                            nextInputDataIndex[dataIndex] = -1;
                        }
                    }
                }
            }
            // Deal with the interior, going from maximum C to minimum C, so that the runStartLinkedLists end up
            // in the right order for efficient forward traversal later
            int csize = ris.m_abcCoordSize.zsize();
            // At each value of c, we process the two scanlines at [*,b,c-1] and [*,b,c] together
            for( int c = csize - 1; c != 0; --c ) {
                for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {
                    int bcIndexA = b + ( c - 1 ) * bsize;
                    int bcIndexB = bcIndexA + bsize;
                    int runRangeStartA = ris.m_bcToRunIndex[bcIndexA],
                        runRangeEndA = ris.m_bcToRunIndex[bcIndexA + 1] - 1;
                    int runRangeStartB = ris.m_bcToRunIndex[bcIndexB],
                        runRangeEndB = ris.m_bcToRunIndex[bcIndexB + 1] - 1;

                    int dataIndexA = exteriorRegionCode, dataIndexB = exteriorRegionCode;

                    const run_data* rdA = &ris.m_runData[runRangeStartA];
                    const run_data* rdB = &ris.m_runData[runRangeStartB];
                    int runA = runRangeStartA, runB = runRangeStartB;
                    int xA = rdA->x, xB = rdB->x;
                    // Deal with empty scanlines
                    if( xA == ( rdA + 1 )->x ) {
                        runA = runRangeEndA;
                        xA = numeric_limits<int>::max();
                    }
                    if( xB == ( rdB + 1 )->x ) {
                        runB = runRangeEndB;
                        xB = numeric_limits<int>::max();
                    }
                    do {
                        int x;
                        // Increment to the next sub-run
                        if( xA < xB ) {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                        } else if( xA > xB ) {
                            x = xB;
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        } else {
                            x = xA;
                            if( runA != runRangeEndA ) {
                                dataIndexA = rdA->dataIndex;
                                ++rdA;
                                ++runA;
                                xA = rdA->x;
                            } else {
                                dataIndexA = exteriorRegionCode;
                                xA = numeric_limits<int>::max();
                            }
                            if( runB != runRangeEndB ) {
                                dataIndexB = rdB->dataIndex;
                                ++rdB;
                                ++runB;
                                xB = rdB->x;
                            } else {
                                dataIndexB = exteriorRegionCode;
                                xB = numeric_limits<int>::max();
                            }
                        }

                        // Figure out the size of this sub-run
                        int size = min( xA, xB ) - x;

                        //
                        if( dataIndexA >= 0 ) {
                            if( dataIndexB >= 0 ) {
                                // In this case we have two defined runs
                                for( int i = 0; i != size; ++i ) {
                                    nextInputDataIndex[dataIndexA] = dataIndexB;
                                    ++dataIndexA;
                                    ++dataIndexB;
                                }
                            } else {
                                // in this case just A is defined
                                for( int i = 0; i != size; ++i ) {
                                    nextInputDataIndex[dataIndexA] = -1;
                                    ++dataIndexA;
                                }
                                // there is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.z + c;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    // Get the index in the linked list head array
                                    int index = a + b * ris.m_abcCoordSize.xsize();
                                    //									cout << "index " << index << ", coord " << a << ", " << b << ",
                                    //in size " << linkedListHeadArray.size() << endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = linkedListHeadArray[index];
                                    linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                                    runStartLinkedLists.push_back( rsl );
                                }
                            }
                        } else {
                            if( dataIndexB >= 0 ) {
                                // In this case just B is defined, so there is the start
                                // of a defined run for each voxel in this sub-run along scanline B.
                                detail::bfap_run_start_ll rsl;
                                rsl.targetX = ris.m_abcCoordOrigin.z + c;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    // Get the index in the linked list head array
                                    int index = a + b * ris.m_abcCoordSize.xsize();
                                    //									cout << "index " << index << ", coord " << a << ", " << b << ",
                                    //in size " << linkedListHeadArray.size() << endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.dataIndex = dataIndexB;
                                    ++dataIndexB;
                                    rsl.nextRunStart = linkedListHeadArray[index];
                                    linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                                    runStartLinkedLists.push_back( rsl );
                                }
                            } else if( dataIndexA != dataIndexB ) {
                                // In this case both are undefined, but the region code changed.
                                // There is the start of an undefined run for each voxel in this sub-run along scanline
                                // B
                                detail::bfap_run_start_ll rsl;
                                rsl.dataIndex = dataIndexB;
                                rsl.targetX = ris.m_abcCoordOrigin.z + c;
                                int a = x - ris.m_abcCoordOrigin.x;
                                for( int i = 0; i != size; ++i, ++a ) {
                                    // Get the index in the linked list head array
                                    int index = a + b * ris.m_abcCoordSize.xsize();
                                    //									cout << "index " << index << ", coord " << a << ", " << b << ",
                                    //in size " << linkedListHeadArray.size() << endl;
                                    // Hook this run start into the front of the scanline linked list
                                    rsl.nextRunStart = linkedListHeadArray[index];
                                    linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                                    runStartLinkedLists.push_back( rsl );
                                }
                            }
                        }
                    } while( xB != numeric_limits<int>::max() || xA != numeric_limits<int>::max() );
                }
            }
            // Deal with the boundary case for the linkedListHeadArray array
            for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {
                // the leading Z's
                int bcIndex = b;
                int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

                const run_data* rd = &ris.m_runData[runRangeStart];
                detail::bfap_run_start_ll rsl;
                rsl.targetX = ris.m_abcCoordOrigin.z + 0;
                for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                    int dataIndex = rd->dataIndex;
                    if( dataIndex >= 0 ) {
                        int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                        int a = rd->x - ris.m_abcCoordOrigin.x;
                        for( ; dataIndex != dataIndexEnd; ++dataIndex, ++a ) {
                            // Get the index in the linked list head array
                            int index = a + b * ris.m_abcCoordSize.xsize();
                            // Hook this run start into the front of the scanline linked list
                            rsl.dataIndex = dataIndex;
                            rsl.nextRunStart = linkedListHeadArray[index];
                            linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                            runStartLinkedLists.push_back( rsl );
                        }
                    } else if( dataIndex != exteriorRegionCode ) {
                        int a = rd->x - ris.m_abcCoordOrigin.x, aEnd = a + ( rd + 1 )->x - rd->x;
                        rsl.dataIndex = dataIndex;
                        for( ; a != aEnd; ++a ) {
                            // Get the index in the linked list head array
                            int index = a + b * ris.m_abcCoordSize.xsize();
                            // Hook this run start into the front of the scanline linked list
                            rsl.nextRunStart = linkedListHeadArray[index];
                            linkedListHeadArray[index] = (int)runStartLinkedLists.size();
                            runStartLinkedLists.push_back( rsl );
                        }
                    }
                }
            }
        }

        /////////
        // Now use the data structures created to build up the runs and copy the data for the target.
        /////////

        //		fout << "Perm: " << axisPermutation << "\n";
        //		fout << "Found " << runStartLinkedLists.size() << " run starts\n";
        //		fout << "Outer bounds: " << outer_bounds() << "\n";

        int nextOutputDataIndex = 0;
        vector3 xyz;
        // Iterate through the voxel coordinates and build up the arrays.
        for( int c = 0; c < m_abcCoordSize.zsize(); ++c ) {
            for( int b = 0; b < m_abcCoordSize.ysize(); ++b ) {
                xyz.y = b + m_abcCoordOrigin.y;
                xyz.z = c + m_abcCoordOrigin.z;

                //				fout << "BC coordinate (" << b << ", " << c << ")\n";

                int bcIndex = b + c * m_abcCoordSize.ysize();
                m_bcToRunIndex[bcIndex] = (int)m_runData.size();

                // The linkedListHeadArray was indexed based on the input, not the target, so we have to check how to
                // get its index here.
                int runStartLLIndex;

                // Test whether to swich the role of b and c in the linkedListHeadArray
                if( axisPermutation.z == 0 )
                    runStartLLIndex = linkedListHeadArray[c + b * m_abcCoordSize.zsize()];
                else
                    runStartLLIndex = linkedListHeadArray[bcIndex];

                // The default state is an undefined region of the exterior region code
                int creatingRegionCode = m_exteriorRegionCode;
                while( runStartLLIndex != -1 ) {
                    const detail::bfap_run_start_ll& rsll = runStartLinkedLists[runStartLLIndex];
                    int inputDataIndex = rsll.dataIndex;

                    //					fout << "rsll.targetX: " << rsll.targetX << "\n";
                    //					fout << "rsll.dataIndex: " << rsll.dataIndex << "\n";
                    //					fout << "rsll.nextRunStart: " << rsll.nextRunStart << "\n";

                    if( inputDataIndex != creatingRegionCode ) {
                        creatingRegionCode = inputDataIndex;
                        if( inputDataIndex < 0 ) {
                            m_runData.push_back( run_data( rsll.targetX, inputDataIndex ) );
                        } else {
                            m_runData.push_back( run_data( rsll.targetX, nextOutputDataIndex ) );

                            do {
                                // Copy this value for each of the channels being remapped
                                for( const rle_index_spec_channel_copying_data* i = channelsToRemapBegin;
                                     i != channelsToRemapEnd; ++i ) {
                                    //							cout << "channel " << i << ", inputDataIndex " << inputDataIndex <<
                                    //", nextOutputDataIndex " << nextOutputDataIndex << endl;
                                    // cout << "channelsToRemap[i].outputData + channelsToRemap[i].primitiveSize *
                                    // nextOutputDataIndex " << (void*)(channelsToRemap[i].outputData +
                                    // channelsToRemap[i].primitiveSize * nextOutputDataIndex) << endl; 							cout
                                    // << "channelsToRemap[i].inputData + channelsToRemap[i].primitiveSize *
                                    //inputDataIndex " << (void*)(channelsToRemap[i].inputData +
                                    //channelsToRemap[i].primitiveSize * inputDataIndex) << endl; 							cout <<
                                    //"channelsToRemap[i].primitiveSize " << channelsToRemap[i].primitiveSize << endl;
                                    // cout << "isbadreadptr: " << IsBadReadPtr(channelsToRemap[i].outputData +
                                    // channelsToRemap[i].primitiveSize * nextOutputDataIndex,
                                    // channelsToRemap[i].primitiveSize) << endl; 							cout << "isbadreadptr: " <<
                                    //IsBadReadPtr(channelsToRemap[i].inputData + channelsToRemap[i].primitiveSize *
                                    // inputDataIndex, channelsToRemap[i].primitiveSize) << endl;
                                    memcpy( i->outputData + i->primitiveSize * nextOutputDataIndex,
                                            i->inputData + i->primitiveSize * inputDataIndex, i->primitiveSize );
                                }
                                ++nextOutputDataIndex;
                                inputDataIndex = nextInputDataIndex[inputDataIndex];
                            } while( inputDataIndex >= 0 );
                        }
                    }

                    runStartLLIndex = rsll.nextRunStart;
                }
                // If any runs were flagged, finish off the last run.  Otherwise add a zero-sized run
                if( (int)m_runData.size() > m_bcToRunIndex[bcIndex] ) {
                    // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to
                    // optimize the storage and access.
                    if( creatingRegionCode != m_exteriorRegionCode )
                        m_runData.push_back( run_data( m_abcCoordOrigin.x + m_abcCoordSize.xsize(), -1 ) );
                } else {
                    m_runData.push_back( run_data( 0, -1 ) );
                    m_runData.push_back( run_data( 0, -1 ) );
                }
            }
        }
        // Finish off the m_bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of
        // the last scanline
        m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();
    }

    /*
    std::stringstream msg2;
    if( !check_consistency(msg2) ) {
      fout << "\n\nERROR!!!: " << msg2.str() << "\n\n";
      ris.dump(fout);
      fout << "\n\n\n\n\n\n";
      dump(fout);
      throw std::runtime_error( "rle_index_spec.build_from_axis_permutation: The rle_index_spec produced is not
    valid.\n"
    + msg2.str() );
    }
    */
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
