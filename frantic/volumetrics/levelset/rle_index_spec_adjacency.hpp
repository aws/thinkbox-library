// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

// Forward declaration of the rle_index_spec class
class rle_index_spec;

/**
 * @todo  Recommended to refactor into the ordering x_neg, x_pos, y_neg, y_pos, z_neg, z_pos.  This will make it
 *        more consistent with some other things, and seems more intuitive to me.  -Mark Wiebe
 */
enum rae_index { rae_index_x_pos, rae_index_x_neg, rae_index_y_pos, rae_index_y_neg, rae_index_z_pos, rae_index_z_neg };

/**
 * This returns true if the neighbor index is rae_index_x_pos, rae_index_y_pos or rae_index_z_pos (assuming its within
 * bounds).
 *
 * @todo  If the ordering is refactored as recommended, this test will flip to !=.
 */
inline bool is_neighbor_index_direction_positive( int neighborIndex ) { return ( neighborIndex & 0x1 ) == 0; }

/**
 * This returns the axis (0, 1 or 2 for X, Y Z respectively) of the given neighbor index.
 */
inline int neighbor_index_axis( int neighborIndex ) { return neighborIndex >> 1; }

/**
 *  This returns the neighbor rae_index for the specified axis and sign.
 *
 * @param axis the axis of the desired neighbor index. Acceptable values
 *		are 0, 1 or 2 for X, Y Z respectively.
 * @param direction +1 to get the index in the positive direction
 *		along the axis (rae_index_x_pos, rae_index_y_pos, or rae_index_z_pos),
 *		or -1 to get the index in the negative direction instead
 *		(rae_index_x_neg, rae_index_y_neg, or rae_index_z_neg).
 * @return the frantic::volumetrics::levelset::rae_index corresponding
 *		to the specified axis and direction.
 */
inline int get_rae_neighbor_index( const int axis, const int direction ) {
    switch( axis ) {
    case 0:
        return direction > 0 ? rae_index_x_pos : rae_index_x_neg;
    case 1:
        return direction > 0 ? rae_index_y_pos : rae_index_y_neg;
    case 2:
        return direction > 0 ? rae_index_z_pos : rae_index_z_neg;
    default:
        throw( "get_rae_neighbor_index Error: axis must be 0, 1, or 2, but instead it was " +
               boost::lexical_cast<std::string>( axis ) + "." );
    }
}

/**
 *  This returns the rae_index neighbor along the specified axis
 * in the positive direction.
 *
 * @param axis the axis along which to find the rae_index neighbor
 *		in the positive direction.  Acceptable values are 0, 1, or 2
 *		for X, Y, and Z respectively.
 * @return the frantic::volumetrics::levelset::rae_index corresponding
 *		to the specified axis in the positive direction.
 */
inline int get_positive_rae_neighbor_index( const int axis ) {
    const int positiveIndices[] = { rae_index_x_pos, rae_index_y_pos, rae_index_z_pos };
    return positiveIndices[axis];
}

/**
 *  This returns the rae_index neighbor along the specified axis
 * in the negative direction.
 *
 * @param axis the axis along which to find the rae_index neighbor
 *		in the negative direction.  Acceptable values are 0, 1, or 2
 *		for X, Y, and Z respectively.
 * @return the frantic::volumetrics::levelset::rae_index corresponding
 *		to the specified axis in the negative direction.
 */
inline int get_negative_rae_neighbor_index( const int axis ) {
    const int negativeIndices[] = { rae_index_x_neg, rae_index_y_neg, rae_index_z_neg };
    return negativeIndices[axis];
}

/**
 *  Return the rae_index neighbour along the same axis, but in
 * the opposite direction.  For example, the opposite of
 * rae_index_x_pos is rae_index_x_neg.
 *
 * @param face the rae index to get the opposite of.
 * @return the rae index that is along the same axis as the
 *		specified neighborIndex, but in the opposite direction.
 */
inline int get_opposite_rae_neighbor_index( const int face ) {
    switch( face ) {
    case rae_index_x_pos:
        return rae_index_x_neg;
    case rae_index_x_neg:
        return rae_index_x_pos;
    case rae_index_y_pos:
        return rae_index_y_neg;
    case rae_index_y_neg:
        return rae_index_y_pos;
    case rae_index_z_pos:
        return rae_index_z_neg;
    case rae_index_z_neg:
        return rae_index_z_pos;
    default:
        throw std::runtime_error( "get_opposite_rae_neighbor_index Error: invalid face number" );
    }
}

/**
 * Each element specifies the data index of the voxel one unit in
 * the positive or negative direction, along the specified axis.
 *
 * @todo  Recommended to refactor into the ordering x_neg, x_pos, y_neg, y_pos, z_neg, z_pos.  This will make it
 *        more consistent with some other things, and seems more intuitive to me.  -Mark Wiebe
 */
struct ris_adj_entry {
    boost::int32_t x_pos, x_neg;
    boost::int32_t y_pos, y_neg;
    boost::int32_t z_pos, z_neg;

    /**
     * This indexing operator lets you do loops for i = 0 to 5.  The order
     * is specified as X positive/negative, Y positive/negative, Z positive/negative.
     */
    const boost::int32_t& operator[]( std::size_t i ) const { return ( &x_pos )[i]; }
};

/**
 * This class stores the adjacency graph of the defined voxels in an rle_index_spec (or possibly other
 * voxel structure).  For each defined voxel, as represented by a data index, it can give you the
 * data index of all 6 immediately adjacent voxels, providing you with the region code if the
 * adjacent voxel is not defined.
 *
 * @todo  I would suggest renaming this to voxel_adjacency_graph, because it is, indeed, providing
 *        the adjacency graph of the defined voxels.  There is nothing rle-specific about it, other
 *        than the creation code, and an ris_adjacency/voxel_adjacency_graph could just as easily be created for a dense
 *        voxel structure as well.  -Mark Wiebe
 */
class ris_adjacency {
    // This is the array of entries
    ris_adj_entry* m_adjacencyList;

    // For some minimal validation by users of this class to make sure it's the right size.
    std::size_t m_dataSize;

    // Don't implement these.  Copy construction and assignment are disabled.
    ris_adjacency( const ris_adjacency& );
    ris_adjacency& operator=( const ris_adjacency& );

  public:
    ris_adjacency() {
        m_adjacencyList = 0;
        m_dataSize = 0;
    }

    ~ris_adjacency() {
        if( m_adjacencyList ) {
            delete[] m_adjacencyList;
            m_adjacencyList = 0;
        }
    }

    const ris_adj_entry& operator[]( std::size_t dataIndex ) const { return m_adjacencyList[dataIndex]; }

    /**
     * This function computes the full adjacency information from the given rle_index_spec.
     */
    void compute( const rle_index_spec& ris );

    std::size_t data_size() const { return m_dataSize; }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
