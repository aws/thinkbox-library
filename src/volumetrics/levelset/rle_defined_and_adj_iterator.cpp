// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_defined_and_adj_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_adjacency.hpp>

using namespace std;
using namespace boost;

namespace frantic {
namespace volumetrics {
namespace levelset {

rle_defined_and_adj_iterator::rle_defined_and_adj_iterator( const rle_index_spec& ris, bool pastTheEnd )
    : m_ris( &ris ) {
    if( !pastTheEnd ) {
        // This loop searches for the first defined voxel
        int32_t csize = ris.m_abcCoordSize.zsize(), bsize = ris.m_abcCoordSize.ysize();
        for( m_c = 0; m_c != csize; ++m_c ) {
            m_coord.z = m_ris->m_abcCoordOrigin.z + m_c;
            for( m_b = 0; m_b != bsize; ++m_b ) {
                m_coord.y = m_ris->m_abcCoordOrigin.y + m_b;
                int bcIndex = m_b + m_c * bsize;
                const int32_t* bcToRI = &m_ris->m_bcToRunIndex[bcIndex];
                const run_data* rd = &m_ris->m_runData[*bcToRI];
                const run_data* rdEnd = rd + *( bcToRI + 1 ) - *bcToRI - 1;
                // If the init function finds a defined voxel in this scanline, we're done
                if( m_rsipCenter.init( *m_ris, rd, rdEnd, m_coord.x, m_dataIndexNeighbors[rae_index_x_neg],
                                       m_dataIndexCenter, m_dataIndexNeighbors[rae_index_x_pos] ) ) {
                    m_rsipYNeg.init( *m_ris, m_b - 1, m_c, m_coord.x, m_dataIndexNeighbors[rae_index_y_neg] );
                    m_rsipYPos.init( *m_ris, m_b + 1, m_c, m_coord.x, m_dataIndexNeighbors[rae_index_y_pos] );
                    m_rsipZNeg.init( *m_ris, m_b, m_c - 1, m_coord.x, m_dataIndexNeighbors[rae_index_z_neg] );
                    m_rsipZPos.init( *m_ris, m_b, m_c + 1, m_coord.x, m_dataIndexNeighbors[rae_index_z_pos] );
                    return;
                }
            }
        }
    }
    // If no defined voxels were found, or pastTheEnd was specified, set up the end iterator
    m_coord = m_ris->m_abcCoordOrigin + m_ris->m_abcCoordSize;
    m_b = m_ris->m_abcCoordSize.ysize();
    m_c = m_ris->m_abcCoordSize.zsize();
    m_rsipCenter.clear();
    m_dataIndexCenter = m_ris->get_exterior_region_code();
}

// Preincrement
rle_defined_and_adj_iterator& rle_defined_and_adj_iterator::operator++() {
    // Call postincrement
    operator++( 0 );
    return *this;
}

// Postincrement
void rle_defined_and_adj_iterator::operator++( int ) {
    // if we are passed the end, don't move anything
    if( !m_rsipCenter.finished() ) {
        // Increment x, and if there are
        int32_t previousX = m_coord.x;
        if( m_rsipCenter.increment_to_next_defined_x( *m_ris, m_coord.x, m_dataIndexNeighbors[rae_index_x_neg],
                                                      m_dataIndexCenter, m_dataIndexNeighbors[rae_index_x_pos] ) ) {
            m_rsipYNeg.increment_to_x( previousX, m_coord.x, *m_ris, m_dataIndexNeighbors[rae_index_y_neg] );
            m_rsipYPos.increment_to_x( previousX, m_coord.x, *m_ris, m_dataIndexNeighbors[rae_index_y_pos] );
            m_rsipZNeg.increment_to_x( previousX, m_coord.x, *m_ris, m_dataIndexNeighbors[rae_index_z_neg] );
            m_rsipZPos.increment_to_x( previousX, m_coord.x, *m_ris, m_dataIndexNeighbors[rae_index_z_pos] );
        } else {
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
                    if( m_rsipCenter.init( *m_ris, rd, rdEnd, m_coord.x, m_dataIndexNeighbors[rae_index_x_neg],
                                           m_dataIndexCenter, m_dataIndexNeighbors[rae_index_x_pos] ) ) {
                        m_rsipYNeg.init( *m_ris, m_b - 1, m_c, m_coord.x, m_dataIndexNeighbors[rae_index_y_neg] );
                        m_rsipYPos.init( *m_ris, m_b + 1, m_c, m_coord.x, m_dataIndexNeighbors[rae_index_y_pos] );
                        m_rsipZNeg.init( *m_ris, m_b, m_c - 1, m_coord.x, m_dataIndexNeighbors[rae_index_z_neg] );
                        m_rsipZPos.init( *m_ris, m_b, m_c + 1, m_coord.x, m_dataIndexNeighbors[rae_index_z_pos] );
                        return;
                    }
                }
                m_b = 0;
            }
            // All things are complete, so this iterator is now past the end
            m_coord = m_ris->m_abcCoordOrigin + m_ris->m_abcCoordSize;
            m_b = m_ris->m_abcCoordSize.ysize();
            m_c = m_ris->m_abcCoordSize.zsize();
            m_rsipCenter.clear();
            m_dataIndexCenter = m_ris->get_exterior_region_code();
        }
    }
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
