// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>

using namespace std;
using namespace boost;

namespace frantic {
namespace volumetrics {
namespace levelset {

rle_defined_iterator::rle_defined_iterator( const rle_index_spec& ris, bool pastTheEnd )
    : m_ris( &ris ) {
    if( !pastTheEnd ) {
        // This loop searches for the first defined voxel
        int32_t csize = m_ris->m_abcCoordSize.zsize(), bsize = m_ris->m_abcCoordSize.ysize();
        for( m_c = 0; m_c != csize; ++m_c ) {
            m_coord.z = m_ris->m_abcCoordOrigin.z + m_c;
            for( m_b = 0; m_b != bsize; ++m_b ) {
                m_coord.y = m_ris->m_abcCoordOrigin.y + m_b;
                int bcIndex = m_b + m_c * bsize;
                const int32_t* bcToRI = &m_ris->m_bcToRunIndex[bcIndex];
                const run_data* rd = &m_ris->m_runData[*bcToRI];
                const run_data* rdEnd = rd + *( bcToRI + 1 ) - *bcToRI - 1;
                // If the init function finds a defined voxel in this scanline, we're done
                if( m_rsip.init( *m_ris, rd, rdEnd, m_coord.x, m_dataIndex ) )
                    return;
            }
        }
    }
    // If no defined voxels were found, or pastTheEnd was specified, set up the end iterator
    m_coord = m_ris->m_abcCoordOrigin + m_ris->m_abcCoordSize;
    m_b = m_ris->m_abcCoordSize.ysize();
    m_c = m_ris->m_abcCoordSize.zsize();
    m_rsip.clear();
    m_dataIndex = m_ris->get_exterior_region_code();
}

// Preincrement
rle_defined_iterator& rle_defined_iterator::operator++() {
    // Call postincrement
    operator++( 0 );
    return *this;
}

// Postincrement
void rle_defined_iterator::operator++( int ) {
    // if we are passed the end, don't move anything
    if( !m_rsip.finished() ) {
        // Increment x, and if there are
        if( !m_rsip.increment_to_next_defined_x( *m_ris, m_coord.x, m_dataIndex ) ) {
            int32_t csize = m_ris->m_abcCoordSize.zsize(), bsize = m_ris->m_abcCoordSize.ysize();
            // No defined run was found, So we need to find the next BC coordinate
            ++m_b;
            if( m_b == bsize ) {
                m_b = 0;
                ++m_c;
            }
            for( ; m_c != csize; ++m_c ) {
                m_coord.z = m_ris->m_abcCoordOrigin.z + m_c;
                for( ; m_b != bsize; ++m_b ) {
                    m_coord.y = m_ris->m_abcCoordOrigin.y + m_b;

                    int bcIndex = m_b + m_c * bsize;
                    const int32_t* bcToRI = &m_ris->m_bcToRunIndex[bcIndex];
                    const run_data* rd = &m_ris->m_runData[*bcToRI];
                    const run_data* rdEnd = rd + *( bcToRI + 1 ) - *bcToRI - 1;
                    // If the init function finds a defined voxel in this scanline, we're done
                    if( m_rsip.init( *m_ris, rd, rdEnd, m_coord.x, m_dataIndex ) )
                        return;
                }
                m_b = 0;
            }
            // All things are complete, so this iterator is now past the end
            m_rsip.clear();
            m_coord = m_ris->m_abcCoordOrigin + m_ris->m_abcCoordSize;
            m_dataIndex = m_ris->get_exterior_region_code();
            m_b = m_ris->m_abcCoordSize.ysize();
            m_c = m_ris->m_abcCoordSize.zsize();
        }
    }
}

int rle_defined_iterator::get_bc_index() const { return m_b + m_c * m_ris->m_abcCoordSize.ysize(); }

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
