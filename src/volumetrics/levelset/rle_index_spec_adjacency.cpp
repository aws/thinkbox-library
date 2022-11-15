// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_adjacency.hpp>

using namespace std;

namespace frantic {
namespace volumetrics {
namespace levelset {

void ris_adjacency::compute( const rle_index_spec& ris ) {
    if( m_adjacencyList ) {
        delete[] m_adjacencyList;
        m_adjacencyList = 0;
    }
    m_dataSize = ris.data_size();

    // If there is no data, then don't allocate anything
    if( m_dataSize == 0 )
        return;

    m_adjacencyList = new ris_adj_entry[m_dataSize];
    int exteriorRegionCode = ris.m_exteriorRegionCode;

    //////////
    // First compute all the x_pos and x_neg values
    //////////
    for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

            int previousRunRegionCode = exteriorRegionCode;
            const run_data* rd = &ris.m_runData[runRangeStart];
            for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
                int dataIndex = rd->dataIndex;
                if( dataIndex >= 0 ) {
                    int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                    if( dataIndex != dataIndexEnd ) {
                        m_adjacencyList[dataIndex].x_neg = previousRunRegionCode;
                        while( dataIndex != dataIndexEnd - 1 ) {
                            m_adjacencyList[dataIndex].x_pos = dataIndex + 1;
                            m_adjacencyList[dataIndex + 1].x_neg = dataIndex;
                            ++dataIndex;
                        }
                        if( run == runRangeEnd - 1 )
                            m_adjacencyList[dataIndex].x_pos = exteriorRegionCode;
                        else
                            m_adjacencyList[dataIndex].x_pos = ( rd + 1 )->dataIndex;
                    }
                } else {
                    previousRunRegionCode = dataIndex;
                }
            }
        }
    }

    //////////
    // Now compute all the y_pos and y_neg values
    //////////

    // Deal with the boundary cases
    for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {
        // the leading Y's
        int bcIndex = 0 + c * ris.m_abcCoordSize.ysize();
        int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

        const run_data* rd = &ris.m_runData[runRangeStart];
        for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
            int dataIndex = rd->dataIndex;
            if( dataIndex >= 0 ) {
                int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                for( ; dataIndex != dataIndexEnd; ++dataIndex ) {
                    m_adjacencyList[dataIndex].y_neg = exteriorRegionCode;
                }
            }
        }

        // the trailing Y's
        bcIndex = ( c + 1 ) * ris.m_abcCoordSize.ysize() - 1;
        runRangeStart = ris.m_bcToRunIndex[bcIndex];
        runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

        rd = &ris.m_runData[runRangeStart];
        for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
            int dataIndex = rd->dataIndex;
            if( dataIndex >= 0 ) {
                int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                for( ; dataIndex != dataIndexEnd; ++dataIndex ) {
                    m_adjacencyList[dataIndex].y_pos = exteriorRegionCode;
                }
            }
        }
    }
    // Deal with the interior
    for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize - 1; ++b ) {
            int bcIndexA = b + c * bsize;
            int bcIndexB = bcIndexA + 1;
            int runRangeStartA = ris.m_bcToRunIndex[bcIndexA], runRangeEndA = ris.m_bcToRunIndex[bcIndexA + 1] - 1;
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
                            m_adjacencyList[dataIndexA].y_pos = dataIndexB;
                            m_adjacencyList[dataIndexB].y_neg = dataIndexA;
                            ++dataIndexA;
                            ++dataIndexB;
                        }
                    } else {
                        // in this case just A is defined
                        for( int i = 0; i != size; ++i ) {
                            m_adjacencyList[dataIndexA].y_pos = dataIndexB;
                            ++dataIndexA;
                        }
                    }
                } else {
                    if( dataIndexB >= 0 ) {
                        // in this case just B is defined
                        for( int i = 0; i != size; ++i ) {
                            m_adjacencyList[dataIndexB].y_neg = dataIndexA;
                            ++dataIndexB;
                        }
                    }
                }
            } while( xB != numeric_limits<int>::max() || xA != numeric_limits<int>::max() );
        }
    }

    //////////
    // Now compute all the z_pos and z_neg values
    //////////

    // Deal with the boundary cases
    for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {
        // the leading Z's
        int bcIndex = b;
        int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

        const run_data* rd = &ris.m_runData[runRangeStart];
        for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
            int dataIndex = rd->dataIndex;
            if( dataIndex >= 0 ) {
                int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                for( ; dataIndex != dataIndexEnd; ++dataIndex ) {
                    m_adjacencyList[dataIndex].z_neg = exteriorRegionCode;
                }
            }
        }

        // the trailing Z's
        bcIndex = b + ( ris.m_abcCoordSize.zsize() - 1 ) * ris.m_abcCoordSize.ysize();
        runRangeStart = ris.m_bcToRunIndex[bcIndex];
        runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;

        rd = &ris.m_runData[runRangeStart];
        for( int run = runRangeStart; run != runRangeEnd; ++run, ++rd ) {
            int dataIndex = rd->dataIndex;
            if( dataIndex >= 0 ) {
                int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                for( ; dataIndex != dataIndexEnd; ++dataIndex ) {
                    m_adjacencyList[dataIndex].z_pos = exteriorRegionCode;
                }
            }
        }
    }
    // Deal with the interior
    for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize - 1; ++c ) {
        for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndexA = b + c * bsize;
            int bcIndexB = bcIndexA + bsize;
            int runRangeStartA = ris.m_bcToRunIndex[bcIndexA], runRangeEndA = ris.m_bcToRunIndex[bcIndexA + 1] - 1;
            int runRangeStartB = ris.m_bcToRunIndex[bcIndexB], runRangeEndB = ris.m_bcToRunIndex[bcIndexB + 1] - 1;

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
                            m_adjacencyList[dataIndexA].z_pos = dataIndexB;
                            m_adjacencyList[dataIndexB].z_neg = dataIndexA;
                            ++dataIndexA;
                            ++dataIndexB;
                        }
                    } else {
                        // in this case just A is defined
                        for( int i = 0; i != size; ++i ) {
                            m_adjacencyList[dataIndexA].z_pos = dataIndexB;
                            ++dataIndexA;
                        }
                    }
                } else {
                    if( dataIndexB >= 0 ) {
                        // in this case just B is defined
                        for( int i = 0; i != size; ++i ) {
                            m_adjacencyList[dataIndexB].z_neg = dataIndexA;
                            ++dataIndexB;
                        }
                    }
                }
            } while( xB != numeric_limits<int>::max() || xA != numeric_limits<int>::max() );
        }
    }
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
