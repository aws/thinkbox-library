// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * This iterator class allows you to iterate over two separate rle_index_spec objects
 * simultaneously, going over all the common sub-intervals of the runs from both.  This
 * is intended for use by operations such as CSG union and intersection, and linear
 * interpolation.
 *
 * Here's an example usage:
 *
 * <pre>
 *  // Make i the begin iterator for the scanline being processed.
 *  rle_pairwise_run_iterator i(risA, risA.y_to_b(y), risA.z_to_c(z), risB, risB.y_to_b(y), risB.z_to_c(z));
 *  // Make ie the past-the-end iterator using the default constructor.
 *  rle_pairwise_run_iterator ie;
 *  for( ; i != ie; ++i ) {
 *    // i.get_xmin() is the starting X coordinate
 *    // i.get_xsize() is the length of the sub-interval
 *    // i.get_first_data_index() is the data index or region code from risA
 *    // i.get_second_data_index() is the data index or region code from risB
 *  }
 *  // At loop end, i.get_xmin() is the X coordinate at the very end of the processed scanline (assuming a non-empty
 * scanline)
 *  // (NOTE that ie.get_xmin() is *NOT* that same X coordinate, note how it was created without knowledge of the
 * scanline)
 * </pre>
 */
class rle_pairwise_run_iterator {
    std::pair<boost::int32_t, boost::int32_t> m_exteriorRegionCode;
    const run_data *m_rdFirst, *m_rdFirstEnd;
    const run_data *m_rdSecond, *m_rdSecondEnd;
    std::pair<boost::int32_t, boost::int32_t> m_dataIndex;

    int m_xStart, m_xSize;

  public:
    /**
     * Constructs an iterator which will iterate over the specified scanlines in risFirst and risSecond in lock-step.
     * If either BC coordinate is out of range, it is treated as having the exterior region code, so the code creating
     * the iterator doesn't have to special-case that situation.
     *
     * If you want to compare two scanlines from the same rle_index_spec, you can use this iterator by passing the
     * same rle_index_spec to risFirst and risSecond.
     *
     * @param  risFirst  The first rle_index_spec for the pairwise iteration.
     * @param  bFirst    The b coordinate for the scanline of the first rle_index_spec.
     * @param  cFirst    The c coordinate for the scanline of the first rle_index_spec.
     * @param  risSecond  The first rle_index_spec for the pairwise iteration.
     * @param  bSecond    The b coordinate for the scanline of the first rle_index_spec.
     * @param  cSecond    The c coordinate for the scanline of the first rle_index_spec.
     */
    rle_pairwise_run_iterator( const rle_index_spec& risFirst, int bFirst, int cFirst, const rle_index_spec& risSecond,
                               int bSecond, int cSecond );

    /**
     * Constructs the 'past-the-end' iterator.
     */
    rle_pairwise_run_iterator();

    /**
     * Pre-increments the iterator.
     *
     * @todo Also do post-increment.
     */
    rle_pairwise_run_iterator& operator++();

    /**
     * For equality comparison, we don't check everything.  But, the subset of what we
     * check is enough to determine equality.  Note that we're checking the things most likely
     * to be different first, to be able to return false as quickly as possible in the common
     * cases.
     */
    bool operator==( const rle_pairwise_run_iterator& iter ) const {
        return m_rdFirst == iter.m_rdFirst && m_rdSecond == iter.m_rdSecond;
    }

    bool operator!=( const rle_pairwise_run_iterator& iter ) const { return !operator==( iter ); }

    /**
     * This gets the data index or the region code of the current sub-interval for risFirst.
     */
    boost::int32_t get_first_data_index() const { return m_dataIndex.first; }

    /**
     * This gets the data index or the region code of the current sub-interval for risSecond.
     */
    boost::int32_t get_second_data_index() const { return m_dataIndex.second; }

    /**
     * This returns the minimum X coordinate of the sub-interval.
     */
    boost::int32_t get_xmin() const { return m_xStart; }

    /**
     * This returns the maximum X coordinate of the sub-interval.
     */
    boost::int32_t get_xmax() const { return m_xStart + m_xSize - 1; }

    /**
     * This returns the size of the sub-interval in the X dimension.
     */
    int get_xsize() const { return m_xSize; }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
