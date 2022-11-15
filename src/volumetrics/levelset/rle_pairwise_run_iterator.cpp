// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_pairwise_run_iterator.hpp>

using namespace std;
using namespace boost;

namespace frantic {
namespace volumetrics {
namespace levelset {

rle_pairwise_run_iterator::rle_pairwise_run_iterator( const rle_index_spec& risFirst, int bFirst, int cFirst,
                                                      const rle_index_spec& risSecond, int bSecond, int cSecond )
    : m_exteriorRegionCode( risFirst.m_exteriorRegionCode, risSecond.m_exteriorRegionCode ) {
    //	cout << "BC coord first: " << bFirst << ", " << cFirst << "\n";
    //	cout << "BC coord second: " << bSecond << ", " << cSecond << "\n";

    // Figure out m_rdFirst and m_rdFirstEnd
    if( (unsigned)bFirst < (unsigned)risFirst.m_abcCoordSize.ysize() &&
        (unsigned)cFirst < (unsigned)risFirst.m_abcCoordSize.zsize() ) {
        int bcIndex = bFirst + cFirst * risFirst.m_abcCoordSize.ysize();
        //		cout << "BC Index first: " << bcIndex << "\n";
        int32_t runIndexStart = risFirst.m_bcToRunIndex[bcIndex];
        const run_data* rd = &risFirst.m_runData[runIndexStart];
        if( rd->x != ( rd + 1 )->x ) {
            m_rdFirst = rd;
            m_rdFirstEnd = rd + ( risFirst.m_bcToRunIndex[bcIndex + 1] - 1 - runIndexStart );
        } else {
            // It's an empty scanline
            m_rdFirst = 0;
            m_rdFirstEnd = 0;
        }
    } else {
        // The BC coordinates are out of range
        m_rdFirst = 0;
        m_rdFirstEnd = 0;
        //		cout << "No first\n";
    }

    // Figure out m_rdSecond and m_rdSecondEnd
    if( (unsigned)bSecond < (unsigned)risSecond.m_abcCoordSize.ysize() &&
        (unsigned)cSecond < (unsigned)risSecond.m_abcCoordSize.zsize() ) {
        int bcIndex = bSecond + cSecond * risSecond.m_abcCoordSize.ysize();
        //		cout << "BC Index second: " << bcIndex << "\n";
        int runIndexStart = risSecond.m_bcToRunIndex[bcIndex];
        const run_data* rd = &risSecond.m_runData[runIndexStart];
        if( rd->x != ( rd + 1 )->x ) {
            m_rdSecond = rd;
            m_rdSecondEnd = rd + ( risSecond.m_bcToRunIndex[bcIndex + 1] - 1 - runIndexStart );
        } else {
            // It's an empty scanline
            m_rdSecond = 0;
            m_rdSecondEnd = 0;
        }
    } else {
        // The BC coordinates are out of range
        m_rdSecond = 0;
        m_rdSecondEnd = 0;
        //		cout << "No second\n";
    }

    // Now set up the iteration parameters
    if( m_rdFirst != 0 ) {
        if( m_rdSecond != 0 ) {
            // Both are non-empty scanlines
            int32_t xFirst = m_rdFirst->x, xSecond = m_rdSecond->x;
            if( xFirst < xSecond ) {
                // The first sub-interval is defined in risFirst, and exterior region code in risSecond
                m_xStart = xFirst;
                m_dataIndex.first = m_rdFirst->dataIndex;
                m_dataIndex.second = m_exteriorRegionCode.second;
                ++m_rdFirst;
                // This sub-interval stops at the end of this risFirst run or the start of the next risSecond run
                m_xSize = std::min( xSecond, m_rdFirst->x ) - xFirst;
            } else if( xSecond < xFirst ) {
                // The first sub-interval is defined in risSecond, and exterior region code in risFirst
                m_xStart = xSecond;
                m_dataIndex.first = m_exteriorRegionCode.first;
                m_dataIndex.second = m_rdSecond->dataIndex;
                ++m_rdSecond;
                // This sub-interval stops at the end of this risSecond run or the start of the next risFirst run
                m_xSize = std::min( xFirst, m_rdSecond->x ) - xSecond;
            } else {
                // The first sub-interval is defined in both risFirst and risSecond
                m_xStart = xFirst;
                m_dataIndex.first = m_rdFirst->dataIndex;
                m_dataIndex.second = m_rdSecond->dataIndex;
                ++m_rdFirst;
                ++m_rdSecond;
                // This sub-interval stops at the end of whichever run is shorter
                m_xSize = std::min( m_rdFirst->x, m_rdSecond->x ) - xFirst;
            }
        } else {
            // Just m_rdFirst is a non-empty scanline
            m_dataIndex.first = m_rdFirst->dataIndex;
            m_dataIndex.second = m_exteriorRegionCode.second;
            m_xStart = m_rdFirst->x;
            ++m_rdFirst;
            m_xSize = m_rdFirst->x - m_xStart;
        }
    } else {
        if( m_rdSecond != 0 ) {
            // Just m_rdSecond is a non-empty scanline
            m_dataIndex.first = m_exteriorRegionCode.first;
            m_dataIndex.second = m_rdSecond->dataIndex;
            m_xStart = m_rdSecond->x;
            ++m_rdSecond;
            m_xSize = m_rdSecond->x - m_xStart;
        } else {
            // Both are empty scanlines
            m_dataIndex.first = m_exteriorRegionCode.first;
            m_dataIndex.second = m_exteriorRegionCode.second;
            m_xStart = 0;
            m_xSize = 0;
        }
    }
}

// Past-the end iterator
rle_pairwise_run_iterator::rle_pairwise_run_iterator() {
    m_exteriorRegionCode.first = 0;
    m_exteriorRegionCode.second = 0;
    m_rdFirst = 0;
    m_rdFirstEnd = 0;
    m_rdSecond = 0;
    m_rdSecondEnd = 0;
    m_dataIndex.first = -1;
    m_dataIndex.second = -1;
    m_xStart = 0;
    m_xSize = 0;
}

rle_pairwise_run_iterator& rle_pairwise_run_iterator::operator++() {
    //	cout << "rdFirst: " << m_rdFirst << ", rdFirstEnd: " << m_rdFirstEnd << "\n";
    //	cout << "rdSecond: " << m_rdSecond << ", rdSecondEnd: " << m_rdSecondEnd << "\n";
    //	cout << "xstart: " << m_xStart << ", xsize: " << m_xSize << "\n";

    // Increment to the next sub-interval
    if( m_rdFirst != 0 ) {
        if( m_rdSecond != 0 ) {
            // Both still have runs to process
            int32_t xFirst = m_rdFirst->x, xSecond = m_rdSecond->x;
            if( xFirst < xSecond ) {
                // The first sub-interval is defined in risFirst, and a continuation of the previous in risSecond
                m_xStart = xFirst;
                // If the continued run within risSecond is defined, need to add to its data index
                if( m_dataIndex.second >= 0 )
                    m_dataIndex.second += m_xSize;
                if( m_rdFirst != m_rdFirstEnd ) {
                    m_dataIndex.first = m_rdFirst->dataIndex;
                    ++m_rdFirst;
                    // This sub-interval stops at the end of this risFirst run or the start of the next risSecond run
                    m_xSize = std::min( xSecond, m_rdFirst->x ) - xFirst;
                } else {
                    m_dataIndex.first = m_exteriorRegionCode.first;
                    m_rdFirst = 0;
                    m_rdFirstEnd = 0;
                    m_xSize = xSecond - xFirst;
                }
            } else if( xSecond < xFirst ) {
                // The first sub-interval is defined in risSecond, and a continuation of the previous in risFirst
                m_xStart = xSecond;
                // If the continued run within risFirst is defined, need to add to its data index
                if( m_dataIndex.first >= 0 )
                    m_dataIndex.first += m_xSize;
                if( m_rdSecond != m_rdSecondEnd ) {
                    m_dataIndex.second = m_rdSecond->dataIndex;
                    ++m_rdSecond;
                    // This sub-interval stops at the end of this risSecond run or the start of the next risFirst run
                    m_xSize = std::min( xFirst, m_rdSecond->x ) - xSecond;
                } else {
                    m_dataIndex.second = m_exteriorRegionCode.second;
                    m_rdSecond = 0;
                    m_rdSecondEnd = 0;
                    m_xSize = xFirst - xSecond;
                }
            } else {
                // The next sub-interval is defined in both risFirst and risSecond
                if( m_rdFirst != m_rdFirstEnd ) {
                    m_dataIndex.first = m_rdFirst->dataIndex;
                    ++m_rdFirst;
                    if( m_rdSecond != m_rdSecondEnd ) {
                        m_dataIndex.second = m_rdSecond->dataIndex;
                        ++m_rdSecond;
                        m_xStart = xFirst;
                        // This sub-interval stops at the end of whichever run is shorter
                        m_xSize = std::min( m_rdFirst->x, m_rdSecond->x ) - xFirst;
                    } else {
                        m_dataIndex.second = m_exteriorRegionCode.second;
                        m_rdSecond = 0;
                        m_rdSecondEnd = 0;
                        // With rdSecond out of the picture, just rdFirst determines the next subinterval
                        m_xStart = xFirst;
                        m_xSize = m_rdFirst->x - xFirst;
                    }
                } else {
                    m_dataIndex.first = m_exteriorRegionCode.first;
                    m_rdFirst = 0;
                    m_rdFirstEnd = 0;
                    if( m_rdSecond != m_rdSecondEnd ) {
                        m_dataIndex.second = m_rdSecond->dataIndex;
                        ++m_rdSecond;
                        // With rdFirst out of the picture, just rdSecond determines the next subinterval
                        m_xStart = xSecond;
                        m_xSize = m_rdSecond->x - xSecond;
                    } else {
                        // Both scanlines have finished, so the iterator is now past-the-end
                        m_dataIndex.second = m_exteriorRegionCode.second;
                        m_rdSecond = 0;
                        m_rdSecondEnd = 0;
                        // When we go past-the-end, get_xmin() becomes equal to the end X coordinate of the scanline we
                        // processed.
                        m_xStart += m_xSize;
                        m_xSize = 0;
                    }
                }
            }
        } else if( m_rdFirst != m_rdFirstEnd ) {
            // Just m_rdFirst is a non-empty scanline
            m_dataIndex.first = m_rdFirst->dataIndex;
            m_xStart = m_rdFirst->x;
            ++m_rdFirst;
            m_xSize = m_rdFirst->x - m_xStart;
        } else {
            // Both have gone past the end
            // When we go past-the-end, get_xmin() becomes equal to the end X coordinate of the scanline we processed.
            m_xStart += m_xSize;
            m_xSize = 0;
            // Set the run_data pointers to 0 so that this iterator compares equal with the past-the-end iterator
            m_rdFirst = 0;
            m_rdFirstEnd = 0;
            m_rdSecond = 0;
            m_rdSecondEnd = 0;
        }
    } else {
        if( m_rdSecond != m_rdSecondEnd ) {
            // Just m_rdSecond is a non-empty scanline
            m_dataIndex.second = m_rdSecond->dataIndex;
            m_xStart = m_rdSecond->x;
            ++m_rdSecond;
            m_xSize = m_rdSecond->x - m_xStart;
        } else {
            // Both have gone past the end
            // When we go past-the-end, get_xmin() becomes equal to the end X coordinate of the scanline we processed.
            m_xStart += m_xSize;
            m_xSize = 0;
            // Set the run_data pointers to 0 so that this iterator compares equal with the past-the-end iterator
            m_rdFirst = 0;
            m_rdFirstEnd = 0;
            m_rdSecond = 0;
            m_rdSecondEnd = 0;
        }
    }

    //	cout << "after xstart: " << m_xStart << ", xsize: " << m_xSize << "\n";

    return *this;
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
