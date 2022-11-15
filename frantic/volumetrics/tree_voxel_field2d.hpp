// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file tree_voxel_field2d.hpp
 * This file provides the class tree_voxel_field2d which is a 2-dimensional voxel field
 * compressed using a heirarchical grid.
 */

#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/volumetrics/tree_voxel_field2d_block_iterator.hpp>
#include <vector>

namespace frantic {
namespace volumetrics {

namespace detail {
// These templates run a meta-program that calculates the largest power of N that is less than 2^32-1
// To get the value, do:
//    static const unsigned long result = max_power<N>::value;
//
template <unsigned long N, unsigned long V = 1, bool C = ( V * N < 0xFFFFFFFFul / N )>
struct max_32bit_power {
    static const unsigned long value = max_32bit_power<N, V * N, ( V * N < 0xFFFFFFFFul / N )>::value;
};

template <unsigned long N, unsigned long V>
struct max_32bit_power<N, V, false> {
    static const unsigned long value = V;
};

/**
 * This struct is used as an individual node in a tree_voxel_field2d.
 */
struct tree_voxel_field2d_node {
    unsigned x, y;
    unsigned length;
    unsigned dataIndex;
};
} // namespace detail

template <unsigned Subdivs = 5>
class tree_voxel_field2d {
    typedef detail::tree_voxel_field2d_node node_type;

    // NOTE: The maximum width is the closest power of Subdivs that is less than 2^32-1
    static const unsigned MAXWIDTH = detail::max_32bit_power<Subdivs>::value;

    // TODO: Could a std::deque be better?
    std::vector<node_type> m_nodeData;

    // std::size_t m_dataItemSize;
    frantic::channels::channel_map m_leafDataMap;
    frantic::graphics::raw_byte_buffer m_leafData;

    inline static bool is_leaf( const node_type& node ) { return ( node.length == Subdivs ); }
    inline static bool is_internal( const node_type& node ) { return ( node.length > Subdivs ); }
    inline static bool has_child( const node_type& node ) { return node.dataIndex != 0xFFFFFFFF; }

    inline static bool contains( const node_type& node, unsigned x, unsigned y ) {
        return ( node.x <= x && node.y <= y && node.x + node.length > x && node.y + node.length > y );
    }

  public:
    friend class tree_voxel_field2d_block_iterator<Subdivs>;
    typedef tree_voxel_field2d_block_iterator<Subdivs> iterator;

  private:
    /**
     * This function will generate the root node for an empty tree, such that the tree
     * contains the voxel [x,y]. It will generate 1 root, leaf node that is not pointing
     * at any data yet.
     * @param x The x-coord for the voxel to contain.
     * @param y The y-coord for the voxel to contain.
     */
    void construct_root_at( unsigned x, unsigned y );

    /**
     * This function will iteratively grow the tree upwards until it contains the
     * specified voxel [x,y].
     * @param x The x-coord for the voxel to contain.
     * @param y The y-coord for the voxel to contain.
     */
    void ensure_root_contains( unsigned x, unsigned y );

    /**
     * This function will generate the next level of nodes below this internal node and
     * update its dataIndex.
     * @note This method takes an index because it grows the node array, which may
     *        invalidate a ptr or ref if they were used instead.
     * @param nodeIndex The index of the internal node for which to populate its children
     */
    void populate_internal_children( unsigned nodeIndex );

    /**
     * Will allocate Subdivs * Subdivs voxels and set the leaf node's dataIndex field to point
     * at the first of these newly allocated voxels.
     *
     * @param nodeIndex The index of the leaf node to allocate voxel data for.
     */
    void populate_leaf_children( unsigned nodeIndex );

    /**
     * This function will return the leaf node containing the position specified. If
     * the leaf and ancestor nodes do not already exist, NULL is returned.
     *
     * @result A pointer to a leaf node, or NULL if no leaf containing that voxel coordinate exists.
     *
     * @param x The x position of the voxel
     * @param y The y position of the voxel
     * @return A ptr to the leaf node containing the voxel [x,y], or NULL if the is no such leaf.
     */
    const node_type* get_leaf_containing( unsigned x, unsigned y ) const;

    /**
     * This function will return the deepest node containing the specified position. If the deepest
     * node is an internal node, or a leaf node without data then outIsValidLeaf is left as false.
     * If the deepest node is a leaf with allocated data then outIsValidLeaf is set to true.
     *
     * @param x The x position of the voxel
     * @param y The y position of the voxel
     * @param outIsValidLeaf Must be set to false on before calling this. Is set to true if the result
     *                        node is a leaf with a valid data ptr.
     * @return A ptr to the deepest node containing the voxel [x,y], or NULL if the root node does not
     *         contain the point.
     */
    const node_type* get_deepest_containing( unsigned x, unsigned y, bool& outIsValidLeaf ) const;

    /**
     * This function will return the leaf node containing the position specified. If
     * the leaf and ancestor nodes do not already exist, they are created and outIsNew is set to true.
     *
     * @note outIsNew will not be modified if the leaf already exists, so it to to false before calling.
     *
     * @param x The x position of the voxel
     * @param y The y position of the voxel
     * @param outIsNew Is set to true if the leaf was created on this call.
     * @return A reference to the leaf node containing the voxel [x,y]
     */
    node_type& get_or_make_leaf_containing( unsigned x, unsigned y, bool& outIsNew );

    /**
     * This function will return the data index of the specified voxel, or -1 if it does not exist.
     *
     * @param x The x position of the voxel
     * @param y The y position of the voxel
     * @return The data index of the voxel, or -1 if it does not exist.
     */
    int get_data_index( unsigned x, unsigned y ) const;

    /**
     * This function will return the data index of the specified voxel. If it does not currently
     * exist, the tree will be modified such that the specified voxel is defined, and the appropriate
     * index is returned.
     *
     * @param x The x position of the voxel
     * @param y The y position of the voxel
     * @return The data index of the voxel
     */
    int get_or_make_data_index( unsigned x, unsigned y );

    /**
     * This function will fill an array of data indices to voxels in the box with corners [x,y]
     * and [x+filterWidth-1,y+filterWidth-1]. Any undefined voxels will have a -1 index in the
     * output array.
     *
     * @param x The x position of the box
     * @param y The y position of the box
     * @param filterWidth The width of the box of voxels to retrieve
     * @param outIndices The array to fill with indices of the voxels. Must be of
     *                    size (filterWidth * filterWidth)
     */
    void get_data_indices( unsigned x, unsigned y, unsigned filterWidth, int outIndices[] ) const;

    /**
     * This function will fill an array of data indices to voxels in the box with corners [x,y]
     * and [x+filterWidth-1,y+filterWidth-1]. Any leaf nodes which are not defined but overlap this
     * box will be created so that all indices are to valid voxels.
     *
     * @param x The x position of the box
     * @param y The y position of the box
     * @param filterWidth The width of the box of voxels to retrieve
     * @param outIndices The array to fill with indices of the voxels. Must be of
     *                    size (filterWidth * filterWidth)
     * @param outCreatedBounds Is set to the bounding box of all new leaf nodes created by this call
     */
    void get_or_make_data_indices( unsigned x, unsigned y, unsigned filterWidth, int outIndices[],
                                   frantic::graphics2d::boundrect2& outCreatedBounds );

    /**
     * This function will collect the 4 indices required for bilinear interpolation of the
     * 2x2 box with corner [x,y]. It is provided as an optimization over get_data_indices(), which
     * is a more general version of this function.
     *
     * @param x The x position of the 2x2 box
     * @param y The y position of the 2x2 box
     * @param outIndices The 4 element array to be filled with the bi-lerp indices.
     */
    void get_bilerp_indices( unsigned x, unsigned y, int outIndices[] ) const;

  public:
    tree_voxel_field2d() {}

    tree_voxel_field2d( const frantic::channels::channel_map& voxelMap )
        : m_leafDataMap( voxelMap ) {}

    void clear() {
        m_nodeData.clear();
        m_leafData.clear();
    }

    void reset( const frantic::channels::channel_map& voxelMap ) {
        clear();
        m_leafDataMap = voxelMap;
    }

    bool is_empty() const { return m_leafData.empty(); }

    iterator begin() const { return iterator( *this ); }

    iterator end() const { return iterator(); }

    char* operator[]( int i ) { return m_leafData.ptr_at( (std::size_t)i * m_leafDataMap.structure_size() ); }

    const char* operator[]( int i ) const {
        return m_leafData.ptr_at( (std::size_t)i * m_leafDataMap.structure_size() );
    }

    const frantic::channels::channel_map& get_channel_map() const { return m_leafDataMap; }

    /**
     * @overload
     * @see tree_voxel_field2d::get_data_index( unsigned, unsigned )
     */
    inline int get_data_index( int x, int y ) const {
        unsigned offsetX = static_cast<unsigned>( x + (int)( MAXWIDTH / 2 ) );
        unsigned offsetY = static_cast<unsigned>( y + (int)( MAXWIDTH / 2 ) );
        return get_data_index( offsetX, offsetY );
    }

    /**
     * @overload
     * @see tree_voxel_field2d::get_or_make_data_index( unsigned, unsigned )
     */
    inline int get_or_make_data_index( int x, int y ) {
        unsigned offsetX = static_cast<unsigned>( x + (int)( MAXWIDTH / 2 ) );
        unsigned offsetY = static_cast<unsigned>( y + (int)( MAXWIDTH / 2 ) );
        return get_or_make_data_index( offsetX, offsetY );
    }

    /**
     * @overload
     * @see tree_voxel_field2d::get_data_indices( unsigned, unsigned, unsigned, int[] )
     */
    inline void get_data_indices( int x, int y, unsigned filterWidth, int outIndices[] ) const {
        unsigned offsetX = static_cast<unsigned>( x + (int)( MAXWIDTH / 2 ) );
        unsigned offsetY = static_cast<unsigned>( y + (int)( MAXWIDTH / 2 ) );
        get_data_indices( offsetX, offsetY, filterWidth, outIndices );
    }

    /**
     * @overload
     * @see tree_voxel_field2d::get_or_make_data_indices( unsigned, unsigned, unsigned, int[],
     * frantic::graphics2d::boundrect2& )
     */
    inline void get_or_make_data_indices( int x, int y, unsigned filterWidth, int outIndices[],
                                          frantic::graphics2d::boundrect2& outCreatedBounds ) {
        unsigned offsetX = static_cast<unsigned>( x + (int)( MAXWIDTH / 2 ) );
        unsigned offsetY = static_cast<unsigned>( y + (int)( MAXWIDTH / 2 ) );
        get_or_make_data_indices( offsetX, offsetY, filterWidth, outIndices, outCreatedBounds );

        outCreatedBounds.minimum().x -= ( MAXWIDTH / 2 );
        outCreatedBounds.minimum().y -= ( MAXWIDTH / 2 );
        outCreatedBounds.maximum().x -= ( MAXWIDTH / 2 );
        outCreatedBounds.maximum().y -= ( MAXWIDTH / 2 );
    }

    /**
     * @overload
     */
    inline void get_or_make_data_indices( int x, int y, unsigned filterWidth, int outIndices[] ) {
        frantic::graphics2d::boundrect2 tempBounds;
        unsigned offsetX = static_cast<unsigned>( x + (int)( MAXWIDTH / 2 ) );
        unsigned offsetY = static_cast<unsigned>( y + (int)( MAXWIDTH / 2 ) );
        get_or_make_data_indices( offsetX, offsetY, filterWidth, outIndices, tempBounds );
    }

    /**
     * @overload
     */
    inline void get_bilerp_indices( int x, int y, int outIndices[] ) const {
        unsigned offsetX = static_cast<unsigned>( x + (int)( MAXWIDTH / 2 ) );
        unsigned offsetY = static_cast<unsigned>( y + (int)( MAXWIDTH / 2 ) );
        get_bilerp_indices( offsetX, offsetY, outIndices );
    }

    template <class ParticleContainer>
    void get_voxels_as_particles( ParticleContainer& outParticles ) const {
        using frantic::graphics::vector3f;

        frantic::channels::channel_map particleMap;
        particleMap.define_channel<vector3f>( _T("Position") );
        particleMap.union_channel_map( m_leafDataMap );
        particleMap.end_channel_definition();

        frantic::channels::channel_map_adaptor adaptor( particleMap, m_leafDataMap );
        frantic::channels::channel_accessor<vector3f> posAccessor =
            particleMap.get_accessor<vector3f>( _T("Position") );

        outParticles.reset( particleMap );

        // Allocate temp space for the particle.
        char* particle = (char*)alloca( particleMap.structure_size() );

        for( std::vector<node_type>::const_iterator it = m_nodeData.begin(), itEnd = m_nodeData.end(); it != itEnd;
             ++it ) {
            if( !is_leaf( *it ) || !has_child( *it ) )
                continue;

            int signedX = (int)it->x - (int)( MAXWIDTH / 2 );
            int signedY = (int)it->y - (int)( MAXWIDTH / 2 );

            for( int y = 0, c = 0; y < Subdivs; ++y ) {
                for( int x = 0; x < Subdivs; ++x, ++c ) {
                    adaptor.copy_structure( particle, ( *this )[it->dataIndex + c] );
                    posAccessor.get( particle ).set( (float)( signedX + x ), (float)( signedY + y ), 0 );

                    outParticles.push_back( particle );
                }
            }
        }
    }
};

} // namespace volumetrics
} // namespace frantic
