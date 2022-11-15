// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

namespace detail {
/**
 * This class contains the data needed for one rle index spec cursor advancing along X.
 */
struct rle_index_spec_block_iterator_x_data {
    // dataIndex - This is the current data index of this particular iterator.
    // nextX - This is the next X which this cursor can skip to.  Within a defined run, this increments
    //         by one, within an undefined run, this will be the start of the next run.
    // nextRunXStart - Only used when processing a defined run.
    int dataIndex, nextX, nextRunXStart;
    int runIndex, runIndexEnd;
};
} // namespace detail

/**
 * This class manages a 2D grid of cursors in the YZ plane, stepping along the X axis
 * in the positive direction.  It provides convenient ways to jump over large swathes of
 * undefined voxels, so as to be able to take advantage of the rle compression.
 */
class rle_index_spec_block_iterator_x {
    const rle_index_spec& m_ris;
    std::vector<detail::rle_index_spec_block_iterator_x_data> m_cursorBlock;
    int m_yMin, m_yMax, m_zMin, m_zMax;
    // The current X coordinate the block iterator is at, as well as the next one it can skip to
    int m_currentX, m_nextX;

    // Disable the assignment operator
    rle_index_spec_block_iterator_x& operator=( const rle_index_spec_block_iterator_x& );

  public:
    rle_index_spec_block_iterator_x( const rle_index_spec& ris, int xStart, int yMin, int yMax, int zMin, int zMax )
        : m_ris( ris )
        , m_cursorBlock( ( yMax - yMin + 1 ) * ( zMax - zMin + 1 ) )
        , m_yMin( yMin )
        , m_yMax( yMax )
        , m_zMin( zMin )
        , m_zMax( zMax ) {
        m_currentX = xStart;
        m_nextX = ( std::numeric_limits<int>::max )();

        int index = 0;
        for( int z = zMin; z <= zMax; ++z ) {
            for( int y = yMin; y <= yMax; ++y, ++index ) {

                detail::rle_index_spec_block_iterator_x_data& data = m_cursorBlock[index];
                int b = y - ris.m_abcCoordOrigin.y, c = z - ris.m_abcCoordOrigin.z;

                if( (unsigned)b < (unsigned)ris.m_abcCoordSize.ysize() &&
                    (unsigned)c < (unsigned)ris.m_abcCoordSize.zsize() ) {
                    int bcIndex = b + c * ris.m_abcCoordSize.ysize();

                    data.runIndexEnd = ris.m_bcToRunIndex[bcIndex + 1] - 2;
                    // Try to find the run which contains xStart.  If no such run
                    // exists, then the run X range will not contain xStart, so we
                    // have to check that later.
                    data.runIndex = ris.BCIndexXtoRunIndex( bcIndex, xStart );

                    // If the start and end run indexes are equal, we have one of 3 cases.
                    // 1)	Empty scanline.
                    // 2)	Full undefined scanline.
                    // 3)	Full defined scanline.
                    // We don't care about case 1, so toss that one out.
                    if( data.runIndex == data.runIndexEnd ) {
                        const run_data* rd = &ris.m_runData[data.runIndex];
                        if( rd->x == ( rd + 1 )->x ) {
                            data.dataIndex = ris.m_exteriorRegionCode;
                            data.nextRunXStart = ( std::numeric_limits<int>::max )();
                            data.nextX = ( std::numeric_limits<int>::max )();
                            data.runIndexEnd = 0;
                            data.runIndex = 1;
                            continue;
                        }
                    }

                    int runXStart = ris.m_runData[data.runIndex].x;
                    data.nextRunXStart = ris.m_runData[data.runIndex + 1].x;

                    if( xStart < runXStart ) {

                        // In this case, the xStart is before the first run in the scanline.
                        data.nextX = runXStart;
                        data.dataIndex = ris.m_exteriorRegionCode;
                        // Point the run towards one before the first run, so that when it is incremented next we get to
                        // the first run
                        --data.runIndex;
                        data.nextRunXStart = runXStart;
                    } else if( xStart >= data.nextRunXStart ) {
                        // In this case, the xStart is after the last run in the scanline.
                        data.nextX = ( std::numeric_limits<int>::max )();
                        data.dataIndex = ris.m_exteriorRegionCode;
                    } else {
                        data.dataIndex = ris.m_runData[data.runIndex].dataIndex;
                        if( data.dataIndex >= 0 ) {
                            data.nextX = xStart + 1;
                            data.dataIndex += xStart - runXStart;
                        } else {
                            data.nextX = data.nextRunXStart;
                        }
                    }

                    // Compute the minimum nextX
                    if( data.nextX < m_nextX )
                        m_nextX = data.nextX;
                } else {
                    // For invalid BC coordinates, it always returns the exterior region code.
                    data.dataIndex = ris.m_exteriorRegionCode;
                    data.nextRunXStart = ( std::numeric_limits<int>::max )();
                    data.nextX = ( std::numeric_limits<int>::max )();
                    data.runIndexEnd = 0;
                    data.runIndex = 1;
                }
            }
        }
    }

    /**
     * This function returns the X value this block iterator is currently sitting at.
     */
    int current_x() const { return m_currentX; }

    /**
     * This function returns the X value this block iterator can skip to in one step.  All data from
     * current_x() to next_x()-1 is the same.  The purpose of this structure is to be able to gain performance
     * from the rle compression even though we're tracking a whole block of scanlines at once.
     */
    int next_x() const { return m_nextX; }

    /**
     * This function increments current_x() by one.
     *
     * @return True if the cursor advanced to a new set of data indexes, false if not.
     */
    bool increment_x() {
        if( m_currentX + 1 < m_nextX ) {
            ++m_currentX;
            return false;
        } else {
            jump_to_next_x();
            return true;
        }
    }

    /**
     * This function advances to the requested x.
     *
     * @return True if the cursor advanced to a new set of data indexes, false if not.
     */
    bool jump_to_x( int x ) {
        if( x >= m_nextX ) {
            do {
                jump_to_next_x();
            } while( x >= m_nextX );
            m_currentX = x;
            return true;
        } else {
            m_currentX = x;
            return false;
        }
    }

    /**
     * This function advances current_x() to the next_x() value.  This is what you typically want to use.
     */
    void jump_to_next_x() {
        if( m_currentX < ( std::numeric_limits<int>::max )() ) {
            m_currentX = m_nextX;
            m_nextX = ( std::numeric_limits<int>::max )();
            for( std::vector<detail::rle_index_spec_block_iterator_x_data>::iterator i = m_cursorBlock.begin(),
                                                                                     iterEnd = m_cursorBlock.end();
                 i != iterEnd; ++i ) {
                // Increment this particular iterator if necessary
                if( i->nextX <= m_currentX ) {
                    // Increment the iterator to the next run if indicated.
                    // NOTE: Could also test that (i->dataIndex < 0 || i->nextRunXStart <= m_currentX), but based on
                    //       the invariants of this iterator, that shouldn't be necessary.
                    if( i->nextRunXStart <= m_currentX ) {
                        ++i->runIndex;
                        if( i->runIndex <= i->runIndexEnd ) {
                            const run_data* rd = &m_ris.m_runData[i->runIndex];
                            i->dataIndex = rd->dataIndex;
                            i->nextRunXStart = ( rd + 1 )->x;
                            if( i->dataIndex >= 0 )
                                i->nextX = m_currentX + 1;
                            else
                                i->nextX = i->nextRunXStart;
                        } else {
                            // When it goes passed the end, the data index is the exterior region code
                            i->nextX = ( std::numeric_limits<int>::max )();
                            i->dataIndex = m_ris.m_exteriorRegionCode;
                        }
                    } else {
                        // Just increment the data index and nextX values.  Note that i->dataIndex >= 0 in this case,
                        // even though that isn't explicitly tested for above.
                        ++i->dataIndex;
                        ++i->nextX;
                    }
                }
                if( i->nextX < m_nextX )
                    m_nextX = i->nextX;
            }
        }
    }

    /**
     * This function fills the parameter outIndexes with the indexes currently in the block iterator.
     *
     * @param  outIndexes  This is where the indexes go.  It must point to an array of ints large enough to hold
     *                     all the indexes.
     */
    void get_indexes( int* outIndexes ) const {
        const detail::rle_index_spec_block_iterator_x_data* cursorBlockData = &m_cursorBlock[0];
        int zMax = m_zMax, yMax = m_yMax;
        for( int z = m_zMin; z <= zMax; ++z ) {
            for( int y = m_yMin; y <= yMax; ++y, ++outIndexes, ++cursorBlockData ) {
                *outIndexes = cursorBlockData->dataIndex;
            }
        }
    }

    /**
     * This function dumps the current iterator data to an outstream.
     *
     * @param  o  The outstream to dump to.
     */
    void dump( std::ostream& o ) {

        o << std::endl << "---Dumping RLE Index Spec Block Iterator X" << std::endl;
        o << "Current x: " << m_currentX << "  Next x: " << m_nextX << std::endl;
        o << "Iterators in block: " << std::endl;
        int index = 0;
        for( int z = m_zMin; z <= m_zMax; ++z ) {
            for( int y = m_yMin; y <= m_yMax; ++y, ++index ) {
                detail::rle_index_spec_block_iterator_x_data& i = m_cursorBlock[index];
                o << "YX coord (" << y << "," << z << ") runIndex: " << i.runIndex << "  nextX: " << i.nextX
                  << std::endl;
            }
        }
        o << "---Finished Dumping RLE Index Spec Block Iterator X" << std::endl;
    }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
