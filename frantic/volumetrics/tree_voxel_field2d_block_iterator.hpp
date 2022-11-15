// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frantic {
namespace volumetrics {

template <unsigned>
class tree_voxel_field2d;

// TODO: Generalize, to non-2x2 blocks
template <unsigned Subdivs = 4>
class tree_voxel_field2d_block_iterator {
    typedef tree_voxel_field2d<Subdivs> tree_array;

    enum { BLOCK_SIZE = Subdivs * Subdivs };

    const tree_array* m_pTreeArray;
    unsigned m_voxelX, m_voxelY;
    unsigned m_leafDataIndex;

  private:
    // Assume m_voxelX and m_voxelY point to the voxel space coords of desired leaf,
    // change them if that position doesn't exist.
    bool reposition_recursive( unsigned rootIndex = 0 ) {
        typename tree_array::node_type currentRoot = m_pTreeArray->m_nodeData[rootIndex];
        if( !m_pTreeArray->is_leaf( currentRoot ) ) {
            unsigned childLength = currentRoot.length / Subdivs;
            unsigned localOffsetX = ( m_voxelX - currentRoot.x ) / childLength;
            unsigned localOffsetY = ( m_voxelY - currentRoot.y ) / childLength;
            unsigned localIndex = localOffsetX + localOffsetY * Subdivs;

            typename tree_array::node_type nextNode = m_pTreeArray->m_nodeData[currentRoot.dataIndex + localIndex];
            if( m_pTreeArray->has_child( nextNode ) ) {
                if( reposition_recursive( currentRoot.dataIndex + localIndex ) )
                    return true;
            }

            while( ++localIndex < BLOCK_SIZE ) {
                nextNode = m_pTreeArray->m_nodeData[currentRoot.dataIndex + localIndex];
                m_voxelX = nextNode.x;
                m_voxelY = nextNode.y;

                if( m_pTreeArray->has_child( nextNode ) ) {
                    if( reposition_recursive( currentRoot.dataIndex + localIndex ) )
                        return true;
                }
            }

            return false;
        } else {
            m_leafDataIndex = currentRoot.dataIndex;
            return true;
        }
    }

    void init() {
        if( m_pTreeArray->m_nodeData.size() > 0 ) {
            typename tree_array::node_type rootNode = m_pTreeArray->m_nodeData[0];
            m_voxelX = rootNode.x;
            m_voxelY = rootNode.y;
            if( !reposition_recursive() )
                m_voxelX = m_voxelY = UINT_MAX;
        } else {
            m_voxelX = m_voxelY = UINT_MAX;
        }
    }

    void next() {
        unsigned length = Subdivs;
        unsigned childLength = 1;
        do {
            m_voxelX += childLength;
            if( ( m_voxelX % length ) == 0 ) { // wrapped in X
                m_voxelX -= length;
                m_voxelY += childLength;
                if( ( m_voxelY % length ) == 0 ) { // wrapped in Y, go up a level
                    m_voxelY -= length;
                    childLength = length;
                    length *= Subdivs;
                } else
                    break;
            } else
                break;
        } while( 1 );

        // If we changed leaf nodes, we need to reposition m_leafDataIndex, skipping empty space.
        if( childLength > 1 ) {
            if( !reposition_recursive() )
                m_voxelX = m_voxelY = UINT_MAX;
        }
        return;
    }

  public:
    void get_data_indices( int outIndices[4] ) const {
        unsigned offsetX = ( m_voxelX % Subdivs );
        unsigned offsetY = ( m_voxelY % Subdivs );

        if( offsetX + 1 < Subdivs ) {
            if( offsetY + 1 < Subdivs ) {
                unsigned index = m_leafDataIndex + offsetX + Subdivs * offsetY;
                outIndices[0] = index;
                outIndices[1] = index + 1;
                outIndices[2] = index + Subdivs;
                outIndices[3] = index + Subdivs + 1;
            } else {
                // Different Y Node
                unsigned index = m_leafDataIndex + offsetX + Subdivs * offsetY;
                outIndices[0] = index;
                outIndices[1] = index + 1;
                if( const typename tree_array::node_type* pLeaf =
                        m_pTreeArray->get_leaf_containing( m_voxelX, m_voxelY + 1 ) ) {
                    unsigned index = pLeaf->dataIndex + offsetX; // We can assume that y+1 % Subdivs == 0
                    outIndices[2] = index;
                    outIndices[3] = index + 1;
                } else {
                    outIndices[2] = -1;
                    outIndices[3] = -1;
                }
            }
        } else if( offsetY + 1 < Subdivs ) {
            // Different X Node
            unsigned index = m_leafDataIndex + offsetX + Subdivs * offsetY;
            outIndices[0] = index;
            outIndices[2] = index + Subdivs;
            if( const typename tree_array::node_type* pLeaf =
                    m_pTreeArray->get_leaf_containing( m_voxelX + 1, m_voxelY ) ) {
                unsigned index = pLeaf->dataIndex + Subdivs * offsetY; // We can assume that x+1 % Subdivs == 0
                outIndices[1] = index;
                outIndices[3] = index + Subdivs;
            } else {
                outIndices[1] = -1;
                outIndices[3] = -1;
            }
        } else {
            // All different leaf nodes
            outIndices[0] = m_leafDataIndex + offsetX + Subdivs * offsetY;
            outIndices[1] = m_pTreeArray->get_data_index( m_voxelX + 1, m_voxelY );
            outIndices[2] = m_pTreeArray->get_data_index( m_voxelX, m_voxelY + 1 );
            outIndices[3] = m_pTreeArray->get_data_index( m_voxelX + 1, m_voxelY + 1 );
        }
    }

    int get_x() const { return (int)( m_voxelX - ( tree_array::MAXWIDTH / 2 ) ); }

    int get_y() const { return (int)( m_voxelY - ( tree_array::MAXWIDTH / 2 ) ); }

    tree_voxel_field2d_block_iterator& operator++() {
        next();
        return *this;
    }

    bool operator==( const tree_voxel_field2d_block_iterator& rhs ) {
        return ( m_voxelX == rhs.m_voxelX ) && ( m_voxelY == rhs.m_voxelY );
    }

    bool operator!=( const tree_voxel_field2d_block_iterator& rhs ) {
        return ( m_voxelX != rhs.m_voxelX ) || ( m_voxelY != rhs.m_voxelY );
    }

  public:
    tree_voxel_field2d_block_iterator()
        : m_pTreeArray( NULL )
        , m_voxelX( UINT_MAX )
        , m_voxelY( UINT_MAX ) {}

    tree_voxel_field2d_block_iterator( const tree_array& treeArray )
        : m_pTreeArray( &treeArray ) {
        init();
    }
};

} // namespace volumetrics
} // namespace frantic
