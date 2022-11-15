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

///////////////////////////
// rle_scanline_defined_iteration_primitive implementation
///////////////////////////

bool rle_scanline_defined_iteration_primitive::init( const rle_index_spec& /*ris*/, const run_data* rdBegin,
                                                     const run_data* rdEnd, boost::int32_t& outXCurrent,
                                                     boost::int32_t& outDataIndex ) {
    for( ; rdBegin != rdEnd; ++rdBegin ) {
        if( rdBegin->dataIndex >= 0 ) {
            // If we find a defined run, initialize the primitive, and return true
            m_rd = rdBegin;
            m_rdEnd = rdEnd;
            m_xEnd = ( rdBegin + 1 )->x;
            outXCurrent = rdBegin->x;
            outDataIndex = rdBegin->dataIndex;
            return true;
        }
    }
    // If no defined run was found, clear to an invalid state, and return false
    clear();
    return false;
}

bool rle_scanline_defined_iteration_primitive::init( const rle_index_spec& ris, const run_data* rdBegin,
                                                     const run_data* rdEnd, boost::int32_t& outXCurrent,
                                                     boost::int32_t& outXNegDataIndex, boost::int32_t& outDataIndex,
                                                     boost::int32_t& outXPosDataIndex ) {
    int32_t dataIndexPrevious = ris.get_exterior_region_code();
    for( ; rdBegin != rdEnd; ++rdBegin ) {
        if( rdBegin->dataIndex >= 0 ) {
            // If we find a defined run, initialize the primitive, and return true
            m_rd = rdBegin;
            m_rdEnd = rdEnd;
            m_xEnd = ( rdBegin + 1 )->x;
            outXCurrent = rdBegin->x;
            outXNegDataIndex = dataIndexPrevious;
            outDataIndex = rdBegin->dataIndex;
            if( outXCurrent + 1 < m_xEnd ) {
                // If the run has more than one voxel in it, its just the next data index
                outXPosDataIndex = outDataIndex + 1;
            } else {
                // otherwise look at the next run, and if it's not past-the-end, its the starting data
                // index of that run, otherwise it's the exterior region code
                ++rdBegin;
                if( rdBegin != rdEnd )
                    outXPosDataIndex = rdBegin->dataIndex;
                else
                    outXPosDataIndex = ris.get_exterior_region_code();
            }
            return true;
        } else {
            // Save the undefined region code for the data index in the negative X direction
            dataIndexPrevious = rdBegin->dataIndex;
        }
    }
    // If no defined run was found, clear to an invalid state, and return false
    clear();
    return false;
}

void rle_scanline_defined_iteration_primitive::clear() {
    m_rd = 0;
    m_rdEnd = 0;
    m_xEnd = numeric_limits<boost::int32_t>::min();
}

bool rle_scanline_defined_iteration_primitive::jump_to_next_defined_x( const rle_index_spec& /*ris*/,
                                                                       boost::int32_t& xTarget,
                                                                       boost::int32_t& xCurrent,
                                                                       boost::int32_t& dataIndex ) {
    // cout << "xCurrent: " << xCurrent << ", m_dataIndex: " << m_dataIndex << endl;
    //++xCurrent;

    if( xTarget < m_xEnd ) {
        // we can just move directly there
        // dataIndex is always assumed to be nonnegative for this iteration primitive
        dataIndex += xTarget - xCurrent;
        return true;
    } else if( m_rd != m_rdEnd ) {
        ++m_rd;
        for( ; m_rd != m_rdEnd; ++m_rd ) {
            if( m_rd->dataIndex >= 0 ) {
                // If we find a defined run, set the necessary variables return true
                dataIndex = m_rd->dataIndex;
                m_xEnd = ( m_rd + 1 )->x;
                xCurrent = m_rd->x;
                // if we overshoot, then we stop and report where we actually are
                if( xCurrent > xTarget ) {
                    return false;
                } else if( xCurrent == xTarget ) {
                    // we have found the target
                    return true;
                }
                // otherwise we need to jump to it
                else if( xTarget < m_xEnd ) {
                    // we can just move directly there
                    // dataIndex is always assumed to be nonnegative for this iteration primitive
                    dataIndex += xTarget - xCurrent;
                    return true;
                }
                // else keep searching
            }
        }
    }
    // If we found no defined X, set rd and rdEnd to NULL to signal that this->finished() should return true.
    m_rd = 0;
    m_rdEnd = 0;
    return false;
}

bool rle_scanline_defined_iteration_primitive::increment_to_next_defined_x( const rle_index_spec& /*ris*/,
                                                                            boost::int32_t& xCurrent,
                                                                            boost::int32_t& dataIndex ) {
    // cout << "xCurrent: " << xCurrent << ", m_dataIndex: " << m_dataIndex << endl;
    ++xCurrent;
    if( xCurrent < m_xEnd ) {
        // dataIndex is always assumed to be nonnegative for this iteration primitive
        ++dataIndex;
        return true;
    } else if( m_rd != m_rdEnd ) {
        ++m_rd;
        for( ; m_rd != m_rdEnd; ++m_rd ) {
            if( m_rd->dataIndex >= 0 ) {
                // If we find a defined run, set the necessary variables return true
                dataIndex = m_rd->dataIndex;
                m_xEnd = ( m_rd + 1 )->x;
                xCurrent = m_rd->x;
                return true;
            }
        }
    }
    // If we found no defined X, set rd and rdEnd to NULL to signal that this->finished() should return true.
    m_rd = 0;
    m_rdEnd = 0;
    return false;
}

bool rle_scanline_defined_iteration_primitive::increment_to_next_defined_x( const rle_index_spec& ris,
                                                                            boost::int32_t& xCurrent,
                                                                            boost::int32_t& outXNegDataIndex,
                                                                            boost::int32_t& dataIndex,
                                                                            boost::int32_t& outXPosDataIndex ) {
    // cout << "xCurrent: " << xCurrent << ", m_dataIndex: " << m_dataIndex << endl;
    ++xCurrent;
    if( xCurrent < m_xEnd ) {
        // dataIndex is always assumed to be nonnegative for this iteration primitive
        outXNegDataIndex = dataIndex;
        ++dataIndex;
        if( xCurrent + 1 < m_xEnd ) {
            // If the next voxel is still in the same run, its just the next data index
            outXPosDataIndex = dataIndex + 1;
        } else {
            // otherwise look at the next run, and if it's not past-the-end, its the starting data
            // index of that run, otherwise it's the exterior region code
            const run_data* rd = m_rd + 1;
            if( rd != m_rdEnd )
                outXPosDataIndex = rd->dataIndex;
            else
                outXPosDataIndex = ris.get_exterior_region_code();
        }
        return true;
    } else if( m_rd != m_rdEnd ) {
        int32_t dataIndexPrevious = ris.get_exterior_region_code();
        ++m_rd;
        for( ; m_rd != m_rdEnd; ++m_rd ) {
            if( m_rd->dataIndex >= 0 ) {
                // If we find a defined run, set the necessary variables return true
                m_xEnd = ( m_rd + 1 )->x;
                xCurrent = m_rd->x;
                outXNegDataIndex = dataIndexPrevious;
                dataIndex = m_rd->dataIndex;
                if( xCurrent + 1 < m_xEnd ) {
                    // If the run has more than one voxel in it, its just the next data index
                    outXPosDataIndex = dataIndex + 1;
                } else {
                    // otherwise look at the next run, and if it's not past-the-end, its the starting data
                    // index of that run, otherwise it's the exterior region code
                    const run_data* rd = m_rd + 1;
                    if( rd != m_rdEnd )
                        outXPosDataIndex = rd->dataIndex;
                    else
                        outXPosDataIndex = ris.get_exterior_region_code();
                }
                return true;
            } else {
                // Save the undefined region code for the data index in the negative X direction
                dataIndexPrevious = m_rd->dataIndex;
            }
        }
    }
    // If we found no defined X, set rd and rdEnd to NULL to signal that this->finished() should return true.
    m_rd = 0;
    m_rdEnd = 0;
    return false;
}

///////////////////////////
// rle_scanline_iteration_primitive implementation
///////////////////////////

void rle_scanline_iteration_primitive::init( const rle_index_spec& ris, boost::int32_t b, boost::int32_t c,
                                             boost::int32_t xCurrent, boost::int32_t& outDataIndex ) {
    int32_t bsize = ris.m_abcCoordSize.ysize(), csize = ris.m_abcCoordSize.zsize();
    if( (uint32_t)b < (uint32_t)bsize && (uint32_t)c < (uint32_t)csize ) {
        int32_t bcIndex = b + c * bsize;
        const int32_t* bcToRI = &ris.m_bcToRunIndex[bcIndex];
        const run_data* rd = &ris.m_runData[*bcToRI];
        const run_data* rdEnd = rd + *( bcToRI + 1 ) - *bcToRI - 1;
        if( xCurrent < rd->x ) {
            // If this is before the first run, set up xEnd to start at the first run in the range
            outDataIndex = ris.get_exterior_region_code();
            m_rd = rd - 1;
            m_rdEnd = rdEnd;
            m_xEnd = rd->x;
        } else {
            while( rd != rdEnd && xCurrent >= ( rd + 1 )->x )
                ++rd;
            if( rd != rdEnd ) {
                outDataIndex = rd->dataIndex;
                if( outDataIndex >= 0 ) {
                    outDataIndex += xCurrent - rd->x;
                }
                m_rd = rd;
                m_rdEnd = rdEnd;
                m_xEnd = ( rd + 1 )->x;
            } else {
                outDataIndex = ris.get_exterior_region_code();
                m_rd = 0;
                m_rdEnd = 0;
            }
        }
    } else {
        outDataIndex = ris.get_exterior_region_code();
        m_rd = 0;
        m_rdEnd = 0;
    }
}

void rle_scanline_iteration_primitive::clear() {
    m_rd = 0;
    m_rdEnd = 0;
}

void rle_scanline_iteration_primitive::increment_to_x( boost::int32_t xPrevious, boost::int32_t xNew,
                                                       const rle_index_spec& ris, boost::int32_t& dataIndex ) {
    // Once m_rd is zero, the data index never changes from the exterior region code
    if( m_rd != 0 ) {
        if( xNew < m_xEnd ) {
            // If we're within the same run, then increment dataIndex if it's a defined run
            if( dataIndex >= 0 )
                dataIndex += xNew - xPrevious;
        } else {
            // Need to find a new run
            while( m_rd != m_rdEnd && xNew >= ( m_rd + 1 )->x )
                ++m_rd;
            if( m_rd != m_rdEnd ) {
                dataIndex = m_rd->dataIndex;
                if( dataIndex >= 0 ) {
                    dataIndex += xNew - m_rd->x;
                }
                m_xEnd = ( m_rd + 1 )->x;
            } else {
                // We've passed the end of all the runs, so from now on everything is the exterior region code
                dataIndex = ris.get_exterior_region_code();
                m_rd = 0;
                m_rdEnd = 0;
            }
        }
    }
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
