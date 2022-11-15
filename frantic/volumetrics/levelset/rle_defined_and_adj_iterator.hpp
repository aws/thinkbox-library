// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/vector3.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_adjacency.hpp>
#include <frantic/volumetrics/levelset/rle_scanline_iteration_primitives.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

// Forward declaration of the rle_index_spec classes.
class rle_index_spec;
struct run_data;

/**
 * This rle_index_spec iterator is designed to provide convenient iteration over all the defined
 * voxels in the rle_index_spec, while simultaneously providing the data index/region code for each of the 6
 * immediate neighbors.
 */
class rle_defined_and_adj_iterator {
    rle_scanline_defined_iteration_primitive m_rsipCenter;
    rle_scanline_iteration_primitive m_rsipYNeg, m_rsipYPos, m_rsipZNeg, m_rsipZPos;
    boost::int32_t m_dataIndexCenter, m_dataIndexNeighbors[6];
    const rle_index_spec* m_ris;
    frantic::graphics::vector3 m_coord;
    boost::int32_t m_b, m_c;

  public:
    /**
     * Constructor which creates an iterator starting at the beginning.
     *
     * @param  ris  The rle_index_spec for iteration.
     * @param  pastTheEnd If false, creates an iterator at the start, otherwise creates an iterator past the end.
     */
    rle_defined_and_adj_iterator( const rle_index_spec& ris, bool pastTheEnd = false );

    /**
     * Pre-increments the iterator.
     */
    rle_defined_and_adj_iterator& operator++();

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
    bool operator==( const rle_defined_and_adj_iterator& iter ) const {
        return m_dataIndexCenter == iter.m_dataIndexCenter;
    }

    bool operator!=( const rle_defined_and_adj_iterator& iter ) const { return !operator==( iter ); }

    /**
     * Gets the coordinate of the current voxel in iteration.
     */
    const frantic::graphics::vector3& get_coord() const { return m_coord; }

    /**
     * Gets the data index of the current voxel in iteration.
     */
    boost::int32_t get_center_data_index() const { return m_dataIndexCenter; }

    /**
     * Gets the data index of one of the 6 neighbor voxels.
     *
     * @param  neighborIndex  The index of the neighbor.  For now, use the rae_index enumeration in
     * rle_index_spec_adjacency.hpp. (rae_index_x_neg, rae_index_x_pos, ...)
     */
    boost::int32_t get_adjacent_data_index( size_t neighborIndex ) const { return m_dataIndexNeighbors[neighborIndex]; }

    boost::int32_t get_x_neg_data_index() const { return m_dataIndexNeighbors[rae_index_x_neg]; }

    boost::int32_t get_x_pos_data_index() const { return m_dataIndexNeighbors[rae_index_x_pos]; }

    boost::int32_t get_y_neg_data_index() const { return m_dataIndexNeighbors[rae_index_y_neg]; }

    boost::int32_t get_y_pos_data_index() const { return m_dataIndexNeighbors[rae_index_y_pos]; }

    boost::int32_t get_z_neg_data_index() const { return m_dataIndexNeighbors[rae_index_z_neg]; }

    boost::int32_t get_z_pos_data_index() const { return m_dataIndexNeighbors[rae_index_z_pos]; }

    boost::int32_t get_b_coord() const { return m_b; }

    boost::int32_t get_c_coord() const { return m_c; }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
