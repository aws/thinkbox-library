// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/vector3.hpp>
#include <frantic/volumetrics/levelset/rle_scanline_iteration_primitives.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

// Forward declaration of the rle_index_spec classes.
class rle_index_spec;
struct run_data;

/**
 * This rle_index_spec iterator is designed to provide convenient iteration over all the defined
 * voxels in the rle_index_spec.
 */
class rle_defined_iterator {
    rle_scanline_defined_iteration_primitive m_rsip;
    const rle_index_spec* m_ris;
    frantic::graphics::vector3 m_coord;
    boost::int32_t m_dataIndex, m_b, m_c;

  public:
    /**
     * Constructor which creates an iterator starting at the beginning.
     *
     * @param  ris  The rle_index_spec for iteration.
     * @param  pastTheEnd If false, creates an iterator at the start, otherwise creates an iterator past the end.
     */
    rle_defined_iterator( const rle_index_spec& ris, bool pastTheEnd = false );

    /**
     * Pre-increments the iterator.
     */
    rle_defined_iterator& operator++();

    /**
     * Post-increments the iterator.  Void return value as allowed by input iterator concept.
     */
    void operator++( int );

    /**
     * For equality comparison, we don't check everything.  But, the subset of what we
     * check is enough to determine equality.  Note that we're checking the things most likely
     * to be different first, to be able to return false as quickly as possible in the common
     * cases.
     */
    bool operator==( const rle_defined_iterator& iter ) const { return m_dataIndex == iter.m_dataIndex; }

    bool operator!=( const rle_defined_iterator& iter ) const { return !operator==( iter ); }

    /**
     * Returns the voxel coordinate of our current iterator.
     * Remember, this coordinate defines the bottom, left, front of the voxel (not center).
     */
    const frantic::graphics::vector3& get_coord() const { return m_coord; }

    int get_bc_index() const;

    int get_data_index() const { return m_dataIndex; }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
