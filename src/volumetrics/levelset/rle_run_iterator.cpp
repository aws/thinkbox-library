// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_run_iterator.hpp>

using namespace std;
using namespace boost;

namespace frantic {
namespace volumetrics {
namespace levelset {

rle_run_iterator::rle_run_iterator( const rle_index_spec& ris, int b, int c ) {
    m_exteriorRegionCode = ris.get_exterior_region_code();

    if( (unsigned)b < (unsigned)ris.m_abcCoordSize.ysize() && (unsigned)c < (unsigned)ris.m_abcCoordSize.zsize() ) {
        int bcIndex = b + c * ris.m_abcCoordSize.ysize();
        int32_t runIndexStart = ris.m_bcToRunIndex[bcIndex];
        const run_data* rd = &ris.m_runData[runIndexStart];
        m_xStart = rd->x;
        m_xSize = ( rd + 1 )->x - m_xStart;
        if( m_xSize != 0 ) {
            m_dataIndex = rd->dataIndex;
            m_rd = rd;
            m_rdEnd = rd + ( ris.m_bcToRunIndex[bcIndex + 1] - 1 - runIndexStart );
            // In this case, there is a non-empty run, so we're good to go
            return;
        }
    }

    // If we didn't initialize the data and return above, initialize to past-the-end.
    m_xStart = 0;
    m_xSize = 0;
    m_dataIndex = m_exteriorRegionCode;
    m_rd = 0;
    m_rdEnd = 0;
}

rle_run_iterator::rle_run_iterator() {
    // the exterior region code doesnt really matter for the end iterator
    m_exteriorRegionCode = -1;

    // Initialize to past-the-end.
    m_xStart = 0;
    m_xSize = 0;
    m_dataIndex = m_exteriorRegionCode;
    m_rd = 0;
    m_rdEnd = 0;
}

rle_run_iterator& rle_run_iterator::operator++() {
    if( m_rd != 0 ) { // if we are past the end, don't move anything
        if( ++m_rd != m_rdEnd ) {
            m_xStart = m_rd->x;
            m_xSize = ( m_rd + 1 )->x - m_xStart;
            m_dataIndex = m_rd->dataIndex;
        } else {
            // get_xmin() gives the one-past-the-end X value after iteration is complete
            m_xStart = m_rd->x;
            m_xSize = 0;
            m_dataIndex = m_exteriorRegionCode;
            m_rd = 0;
            m_rdEnd = 0;
        }
    }
    return *this;
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
