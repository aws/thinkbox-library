// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frantic {
namespace volumetrics {
namespace levelset {

// Forward declaration of the rle_index_spec classes.
class rle_index_spec;
struct run_data;

/**
 * This rle_run_iterator is designed to provide convenient iteration over all the runs for
 * one scanline in an rle_index_spec.
 *
 * Here's an example usage for one run:
 *
 * <pre>
 *  // Make i the begin iterator for the scanline being processed.
 *  rle_run_iterator i(ris, ris.y_to_b(y), ris.z_to_c(z));
 *  // Make ie the past-the-end iterator using the default constructor.
 *  rle_run_iterator ie;
 *  for( ; i != ie; ++i ) {
 *    // i.get_xmin() is the starting X coordinate
 *    // i.get_xsize() is the length of the run
 *    // i.get_data_index() is the data index or region code for the run
 *  }
 *  // At loop end, i.get_xmin() is the X coordinate at the very end of the processed scanline (assuming a non-empty
 * scanline)
 *  // (NOTE that ie.get_xmin() is *NOT* that same X coordinate, note how it was created without knowledge of the
 * scanline)
 * </pre>
 *
 * Here's how you would typically use it to visit all runs in an rle_index_spec:
 *
 * <pre>
 *  const rle_index_spec& ris = ...;
 *  boundbox3 outerBounds = ris.outer_bounds();
 *  for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
 *    for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
 *      rle_run_iterator i(ris, ris.y_to_b(y), ris.z_to_c(z)), ie;
 *      for( ; i != ie; ++i ) {
 *        // i.get_xmin() is the starting X coordinate
 *        // i.get_xsize() is the length of the run
 *        // i.get_data_index() is the data index or region code for the run
 *      }
 *      // Here, i.get_xmin() is the X coordinate at the very end of the processed scanline (assuming a non-empty
 * scanline)
 *      // (NOTE that ie.get_xmin() is *NOT* that same X coordinate, note how it was created without knowledge of the
 * scanline)
 *    }
 *  }
 * </pre>
 */
class rle_run_iterator {
    const run_data *m_rd, *m_rdEnd;
    int m_exteriorRegionCode;
    boost::int32_t m_dataIndex;

    boost::int32_t m_xStart, m_xSize;

  public:
    /**
     * Constructs an iterator which will iterate over the runs in the scanline of the given BC coordinates.
     *
     * @param  ris  The rle_index_spec for run iteration.
     * @param  b    The b coordinate for the scanline.
     * @param  c    The c coordinate for the scanline.
     */
    rle_run_iterator( const rle_index_spec& ris, int b, int c );

    /**
     * Constructs the 'past-the-end' iterator.
     */
    rle_run_iterator();

    /**
     * Pre-increments the iterator.
     *
     * @todo Do both pre- and post-increment correctly
     */
    rle_run_iterator& operator++();

    /**
     * For equality comparison, we don't check everything.  But, the subset of what we
     * check is enough to determine equality.  Note that we're checking the things most likely
     * to be different first, to be able to return false as quickly as possible in the common
     * cases.
     */
    bool operator==( const rle_run_iterator& iter ) const { return m_rd == iter.m_rd; }

    bool operator!=( const rle_run_iterator& iter ) const { return !operator==( iter ); }

    /**
     * This gets the data index of the run, or the region code of the run if it
     * is an undefined run.
     */
    boost::int32_t get_data_index() const { return m_dataIndex; }

    /**
     * This returns the minimum X coordinate of the run.
     */
    boost::int32_t get_xmin() const { return m_xStart; }

    /**
     * This returns the maximum X coordinate of the run.
     */
    boost::int32_t get_xmax() const { return m_xStart + m_xSize - 1; }

    /**
     * This returns the size of the run in the X dimension.
     */
    boost::int32_t get_xsize() const { return m_xSize; }

    /**
     * This gets the data index of the next run, or the region code of the run if it
     * is an undefined run.  It's more costly than retrieving the current index due
     * to the extra indirection.
     */
    boost::int32_t get_next_run_data_index() const {

        if( m_rd == m_rdEnd )
            return m_exteriorRegionCode;

        const run_data* next_rd = m_rd + 1;
        if( next_rd != m_rdEnd )
            return next_rd->dataIndex;
        else
            return m_exteriorRegionCode;
    }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
