// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/volumetrics/tree_voxel_field2d.hpp>

using frantic::graphics2d::boundrect2;
using frantic::graphics2d::vector2;

namespace frantic {
namespace volumetrics {

template <unsigned Subdivs>
void tree_voxel_field2d<Subdivs>::construct_root_at( unsigned x, unsigned y ) {
    unsigned rootX = ( x / Subdivs ) * Subdivs;
    unsigned rootY = ( y / Subdivs ) * Subdivs;

    m_nodeData.resize( 1 );
    m_nodeData[0].x = rootX;
    m_nodeData[0].y = rootY;
    m_nodeData[0].length = Subdivs;
    m_nodeData[0].dataIndex = 0xFFFFFFFF;
}

template <unsigned Subdivs>
void tree_voxel_field2d<Subdivs>::ensure_root_contains( unsigned x, unsigned y ) {
    while( !contains( m_nodeData[0], x, y ) ) {
        unsigned newIndex = static_cast<unsigned>( m_nodeData.size() );
        m_nodeData.resize( Subdivs * Subdivs + newIndex );

        unsigned currentLength = m_nodeData[0].length;
        unsigned newRootLength = currentLength * Subdivs;
        unsigned rootX = ( m_nodeData[0].x / newRootLength ) * newRootLength;
        unsigned rootY = ( m_nodeData[0].y / newRootLength ) * newRootLength;

        if( newRootLength < currentLength )
            throw std::runtime_error(
                "tree_array2D.ensure_root_contains() - Cannot grow to contain the specified point." );

        unsigned index = newIndex;
        for( unsigned j = 0; j < Subdivs; ++j ) {
            for( unsigned i = 0; i < Subdivs; ++i, ++index ) {
                m_nodeData[index].x = rootX + i * currentLength;
                m_nodeData[index].y = rootY + j * currentLength;
                m_nodeData[index].length = currentLength;
                m_nodeData[index].dataIndex = 0xFFFFFFFF;
            }
        }

        unsigned offsetX = ( m_nodeData[0].x - rootX ) / currentLength;
        unsigned offsetY = ( m_nodeData[0].y - rootY ) / currentLength;
        m_nodeData[newIndex + offsetX + ( offsetY * Subdivs )].dataIndex = m_nodeData[0].dataIndex;
        m_nodeData[0].x = rootX;
        m_nodeData[0].y = rootY;
        m_nodeData[0].length = newRootLength;
        m_nodeData[0].dataIndex = newIndex;
    }
}

template <unsigned Subdivs>
void tree_voxel_field2d<Subdivs>::populate_internal_children( unsigned nodeIndex ) {
    unsigned childLength = m_nodeData[nodeIndex].length / Subdivs;

    unsigned newIndex = static_cast<unsigned>( m_nodeData.size() );
    m_nodeData.resize( newIndex + ( Subdivs * Subdivs ) );

    unsigned index = newIndex;
    for( unsigned j = 0; j < Subdivs; ++j ) {
        for( unsigned i = 0; i < Subdivs; ++i, ++index ) {
            m_nodeData[index].x = m_nodeData[nodeIndex].x + i * childLength;
            m_nodeData[index].y = m_nodeData[nodeIndex].y + j * childLength;
            m_nodeData[index].length = childLength;
            m_nodeData[index].dataIndex = 0xFFFFFFFF;
        }
    }

    m_nodeData[nodeIndex].dataIndex = newIndex;
}

template <unsigned Subdivs>
void tree_voxel_field2d<Subdivs>::populate_leaf_children( unsigned nodeIndex ) {
    std::size_t newIndex = m_leafData.size() / m_leafDataMap.structure_size();
    std::size_t allocSize = ( Subdivs * Subdivs ) * m_leafDataMap.structure_size();

    char* pNewData = m_leafData.add_element( allocSize, ( Subdivs * Subdivs ) );
    memset( pNewData, 0, allocSize );

    m_nodeData[nodeIndex].dataIndex = (unsigned)newIndex;
}

template <unsigned Subdivs>
const typename tree_voxel_field2d<Subdivs>::node_type*
tree_voxel_field2d<Subdivs>::get_leaf_containing( unsigned x, unsigned y ) const {
    if( !contains( m_nodeData[0], x, y ) )
        return NULL;

    unsigned curID = 0;
    node_type curNode = m_nodeData[0];
    while( is_internal( curNode ) ) {
        if( !has_child( curNode ) )
            return NULL;
        unsigned childLength = curNode.length / Subdivs;
        unsigned offsetX = ( x - curNode.x ) / childLength;
        unsigned offsetY = ( y - curNode.y ) / childLength;
        curID = curNode.dataIndex + offsetX + ( offsetY * Subdivs );
        curNode = m_nodeData[curID];
    }

    if( !has_child( curNode ) )
        return NULL;

    return &m_nodeData[curID];
}

template <unsigned Subdivs>
const typename tree_voxel_field2d<Subdivs>::node_type*
tree_voxel_field2d<Subdivs>::get_deepest_containing( unsigned x, unsigned y, bool& outIsValidLeaf ) const {
    outIsValidLeaf = false;

    if( !contains( m_nodeData[0], x, y ) )
        return NULL;

    unsigned curID = 0;
    node_type curNode = m_nodeData[0];
    while( is_internal( curNode ) ) {
        if( !has_child( curNode ) )
            return &m_nodeData[curID];
        unsigned childLength = curNode.length / Subdivs;
        unsigned offsetX = ( x - curNode.x ) / childLength;
        unsigned offsetY = ( y - curNode.y ) / childLength;
        curID = curNode.dataIndex + offsetX + ( offsetY * Subdivs );
        curNode = m_nodeData[curID];
    }

    if( has_child( curNode ) )
        outIsValidLeaf = true;

    return &m_nodeData[curID];
}

template <unsigned Subdivs>
typename tree_voxel_field2d<Subdivs>::node_type&
tree_voxel_field2d<Subdivs>::get_or_make_leaf_containing( unsigned x, unsigned y, bool& outIsNew ) {
    if( m_nodeData.size() == 0 )
        construct_root_at( x, y );
    else
        ensure_root_contains( x, y );

    unsigned curNode = 0;
    while( is_internal( m_nodeData[curNode] ) ) {
        if( !has_child( m_nodeData[curNode] ) )
            populate_internal_children( curNode );

        // Set length to be the child node's length.
        unsigned childLength = m_nodeData[curNode].length / Subdivs;
        unsigned offsetX = ( x - m_nodeData[curNode].x ) / childLength;
        unsigned offsetY = ( y - m_nodeData[curNode].y ) / childLength;
        curNode = m_nodeData[curNode].dataIndex + offsetX + ( offsetY * Subdivs );
    }

    // We have arrived at the relevant leaf node.
    if( !has_child( m_nodeData[curNode] ) ) {
        populate_leaf_children( curNode );
        outIsNew = true;
    }
    return m_nodeData[curNode];
}

template <unsigned Subdivs>
void tree_voxel_field2d<Subdivs>::get_or_make_data_indices( unsigned x, unsigned y, unsigned filterWidth,
                                                            int outIndices[], boundrect2& outCreatedBounds ) {
    // Store the filter's base coordinates, so that x and y can be modified to be the current
    // inspected coordinates as we move through the leaves.
    unsigned filterX = x;
    unsigned filterY = y;

    while( y < filterY + filterWidth ) {
        bool isNewLeaf = false;
        node_type* pCurrentNode = &get_or_make_leaf_containing( x, y, isNewLeaf );

        if( isNewLeaf ) {
            outCreatedBounds += vector2( pCurrentNode->x, pCurrentNode->y );
            outCreatedBounds +=
                vector2( pCurrentNode->x + pCurrentNode->length, pCurrentNode->y + pCurrentNode->length );
        }

        // Calculate the offsets in this leaf node, and in the output array.
        unsigned outXOff = ( x - filterX );
        unsigned outYOff = ( y - filterY );
        unsigned leafXOff = ( x - pCurrentNode->x );
        unsigned leafYOff = ( y - pCurrentNode->y );

        // Calculate the index in the output array, and into the leaf's data indices.
        unsigned outIndex = outXOff + filterWidth * outYOff;
        unsigned dataIndex = pCurrentNode->dataIndex + leafXOff + Subdivs * leafYOff;

        // Calculate the portion of this leaf overlapped by the filter.
        unsigned xWidth = std::min( Subdivs - leafXOff, filterWidth - outXOff );
        unsigned yWidth = std::min( Subdivs - leafYOff, filterWidth - outYOff );

        // Run through the part of the filter that fits into this leaf node, filling in the
        // appropriate output elements.
        for( unsigned r = 0; r < yWidth; ++r ) {
            for( unsigned c = 0; c < xWidth; ++c )
                outIndices[outIndex + c] = dataIndex + c;
            dataIndex += Subdivs;
            outIndex += filterWidth;
        }

        // Try the next leaf in the X direction. If it is beyond our filter width, then
        // try the next leaf in the Y direction. If it is beyond out filter's extent, we are done.
        x += xWidth;
        if( x >= filterX + filterWidth ) {
            x = filterX;
            y += yWidth;
        }
    }
}

template <unsigned Subdivs>
void tree_voxel_field2d<Subdivs>::get_data_indices( unsigned x, unsigned y, unsigned filterWidth,
                                                    int outIndices[] ) const {
    if( m_nodeData.size() == 0 ) {
        for( int i = 0, iEnd = filterWidth * filterWidth; i < iEnd; ++i )
            outIndices[i] = -1;
        return;
    }

    unsigned maxYStep = m_nodeData[0].length;

    // Store the filter's base coordinates, so that x and y can be modified to be the current
    // inspected coordinates as we move through the leaves.
    unsigned filterX = x;
    unsigned filterY = y;

    while( y < filterY + filterWidth ) {
        unsigned xWidth, yWidth;

        bool isValidLeaf = false;
        const node_type* pCurrentNode = get_deepest_containing( x, y, isValidLeaf );

        if( isValidLeaf ) {
            // Calculate the offsets in this leaf node, and in the output array.
            unsigned outXOff = ( x - filterX );
            unsigned outYOff = ( y - filterY );
            unsigned leafXOff = ( x - pCurrentNode->x );
            unsigned leafYOff = ( y - pCurrentNode->y );

            // Calculate the index in the output array, and into the leaf's data indices.
            unsigned outIndex = outXOff + filterWidth * outYOff;
            unsigned dataIndex = pCurrentNode->dataIndex + leafXOff + Subdivs * leafYOff;

            // Calculate the portion of this leaf overlapped by the filter.
            xWidth = std::min( Subdivs - leafXOff, filterWidth - outXOff );
            yWidth = std::min( Subdivs - leafYOff, filterWidth - outYOff );

            // Run through the part of the filter that fits into this leaf node, filling in the
            // appropriate output elements.
            for( unsigned r = 0; r < yWidth; ++r ) {
                for( unsigned c = 0; c < xWidth; ++c )
                    outIndices[outIndex + c] = dataIndex + c;
                dataIndex += Subdivs;
                outIndex += filterWidth;
            }
        } else if( pCurrentNode ) {
            // Calculate the offsets in this non-leaf node, and in the output array.
            unsigned outXOff = ( x - filterX );
            unsigned outYOff = ( y - filterY );
            unsigned leafXOff = ( x - pCurrentNode->x );
            unsigned leafYOff = ( y - pCurrentNode->y );

            // Calculate the index in the output array.
            unsigned outIndex = outXOff + filterWidth * outYOff;

            // Calculate the portion of this node overlapped by the filter.
            xWidth = std::min( pCurrentNode->length - leafXOff, filterWidth - outXOff );
            yWidth = std::min( pCurrentNode->length - leafYOff, filterWidth - outYOff );

            // Run through the part of the filter that fits into this leaf node, filling in the
            // appropriate output elements.
            for( unsigned r = 0; r < yWidth; ++r ) {
                for( unsigned c = 0; c < xWidth; ++c )
                    outIndices[outIndex + c] = -1;
                outIndex += filterWidth;
            }
        } else {
            // There is no node that contains point since pCurrentNode == NULL.
            // Calculate which part of the filter is not overlapped by the root node.
            unsigned outXOff = ( x - filterX );
            unsigned outYOff = ( y - filterY );

            // Calculate the index in the output array.
            unsigned outIndex = outXOff + filterWidth * outYOff;

            // Calculate the portion of this non-node overlapped by the filter.
            xWidth = filterWidth - outXOff;
            if( x < m_nodeData[0].x )
                xWidth = std::min( m_nodeData[0].x - x, xWidth );

            yWidth = filterWidth - outYOff;
            if( y < m_nodeData[0].y )
                yWidth = std::min( m_nodeData[0].y - y, yWidth );

            // Run through the part of the filter that fits into this leaf node, filling in the
            // appropriate output elements.
            for( unsigned r = 0; r < yWidth; ++r ) {
                for( unsigned c = 0; c < xWidth; ++c )
                    outIndices[outIndex + c] = -1;
                outIndex += filterWidth;
            }
        }

        if( yWidth < maxYStep )
            maxYStep = yWidth;

        // Try the next leaf in the X direction. If it is beyond our filter's extent, then
        // try the next leaf in the Y direction. If it is also beyond out filter's extent, we are done.
        x += xWidth;
        if( x >= filterX + filterWidth ) {
            x = filterX;
            // y += yWidth;
            /* If a filter overlaps regions like this it gets tricky to traverse, so I restrict y-movement to the
             smallest amount.
             *---*---*---*---*
             | 0 |           |
             *---*           *
             | 2 |     1     |
             *---*           *
             |   |           |
             *---*---*---*---*
            */
            y += maxYStep;
            maxYStep = m_nodeData[0].length;
        }
    }
}

template <unsigned Subdivs>
int tree_voxel_field2d<Subdivs>::get_data_index( unsigned x, unsigned y ) const {
    if( m_nodeData.size() == 0 || !contains( m_nodeData[0], x, y ) )
        return -1;

    unsigned curNode = 0;
    while( is_internal( m_nodeData[curNode] ) ) {
        if( !has_child( m_nodeData[curNode] ) )
            return -1;

        // Set length to be the child node's length.
        unsigned childLength = m_nodeData[curNode].length / Subdivs;
        unsigned offsetX = ( x - m_nodeData[curNode].x ) / childLength;
        unsigned offsetY = ( y - m_nodeData[curNode].y ) / childLength;
        curNode = m_nodeData[curNode].dataIndex + offsetX + ( offsetY * Subdivs );
    }

    // We have arrived at the relevant leaf node.
    if( !has_child( m_nodeData[curNode] ) )
        return -1;

    unsigned offsetX = ( x - m_nodeData[curNode].x );
    unsigned offsetY = ( y - m_nodeData[curNode].y );
    return m_nodeData[curNode].dataIndex + offsetX + ( offsetY * Subdivs );
}

template <unsigned Subdivs>
int tree_voxel_field2d<Subdivs>::get_or_make_data_index( unsigned x, unsigned y ) {
    if( m_nodeData.size() == 0 )
        construct_root_at( x, y );
    else
        ensure_root_contains( x, y );

    unsigned curNode = 0;
    while( is_internal( m_nodeData[curNode] ) ) {
        if( !has_child( m_nodeData[curNode] ) )
            populate_internal_children( curNode );

        // Set length to be the child node's length.
        unsigned childLength = m_nodeData[curNode].length / Subdivs;
        unsigned offsetX = ( x - m_nodeData[curNode].x ) / childLength;
        unsigned offsetY = ( y - m_nodeData[curNode].y ) / childLength;
        curNode = m_nodeData[curNode].dataIndex + offsetX + ( offsetY * Subdivs );
    }

    // We have arrived at the relevant leaf node.
    if( !has_child( m_nodeData[curNode] ) )
        populate_leaf_children( curNode );

    unsigned offsetX = ( x - m_nodeData[curNode].x );
    unsigned offsetY = ( y - m_nodeData[curNode].y );
    return m_nodeData[curNode].dataIndex + offsetX + ( offsetY * Subdivs );
}

template <unsigned Subdivs>
void tree_voxel_field2d<Subdivs>::get_bilerp_indices( unsigned x, unsigned y, int outIndices[] ) const {
    if( m_nodeData.size() == 0 ) {
        outIndices[0] = outIndices[1] = outIndices[2] = outIndices[3] = -1;
        return;
    }

    unsigned offsetX = x % Subdivs;
    unsigned offsetY = y % Subdivs;
    if( offsetX < ( Subdivs - 1 ) ) {
        if( offsetY < ( Subdivs - 1 ) ) {
            // One node stores all values
            if( const node_type* pLeaf = get_leaf_containing( x, y ) ) {
                unsigned index = pLeaf->dataIndex + offsetX + ( offsetY * Subdivs );
                outIndices[0] = index;
                outIndices[1] = index + 1;
                outIndices[2] = index + Subdivs;
                outIndices[3] = index + Subdivs + 1;
            } else {
                outIndices[0] = outIndices[1] = outIndices[2] = outIndices[3] = -1;
            }
        } else {
            // Two nodes to access
            if( const node_type* pLeaf = get_leaf_containing( x, y ) ) {
                unsigned index = pLeaf->dataIndex + offsetX + ( offsetY * Subdivs );
                outIndices[0] = index;
                outIndices[1] = index + 1;
            } else {
                outIndices[0] = outIndices[1] = -1;
            }

            if( const node_type* pLeaf = get_leaf_containing( x, y + 1 ) ) {
                unsigned index = pLeaf->dataIndex + offsetX; // We can assume that y+1 % Subdivs == 0
                outIndices[2] = index;
                outIndices[3] = index + 1;
            } else {
                outIndices[2] = outIndices[3] = -1;
            }
        }
    } else {
        if( offsetY < ( Subdivs - 1 ) ) {
            // Two nodes to access
            if( const node_type* pLeaf = get_leaf_containing( x, y ) ) {
                unsigned index = pLeaf->dataIndex + offsetX + ( offsetY * Subdivs );
                outIndices[0] = index;
                outIndices[2] = index + Subdivs;
            } else {
                outIndices[0] = outIndices[2] = -1;
            }

            if( const node_type* pLeaf = get_leaf_containing( x + 1, y ) ) {
                unsigned index = pLeaf->dataIndex + ( offsetY * Subdivs ); // We can assume that x+1 % Subdivs == 0
                outIndices[1] = index;
                outIndices[3] = index + Subdivs;
            } else {
                outIndices[1] = outIndices[3] = -1;
            }
        } else {
            // Four nodes to access
            // TODO: Try getting the node above the leaf level and trying there.
            outIndices[0] = get_data_index( x, y );
            outIndices[1] = get_data_index( x + 1, y );
            outIndices[2] = get_data_index( x, y + 1 );
            outIndices[3] = get_data_index( x + 1, y + 1 );
        }
    }
}

template class tree_voxel_field2d<3>;
template class tree_voxel_field2d<4>;
template class tree_voxel_field2d<5>;

} // namespace volumetrics
} // namespace frantic
