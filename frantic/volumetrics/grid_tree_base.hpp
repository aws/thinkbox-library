// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/integer/static_log2.hpp>

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/misc/exception_stream.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

namespace frantic {
namespace volumetrics {

/**
 * grid_tree_base_node
 * This is a struct for the nodes in the grid_tree. It represents a length^3 (of its length member)
 * node starting at location nodeCoord. It is templated on the leaf data class
 * and a side length (to avoid extra heap indirection for child lookups, etc.)
 */
template <class LEAFDATA, int SL>
struct grid_tree_base_node {
    // Node position and length.
    frantic::graphics::vector3 coord;
    boost::uint32_t length;

    // Pointers to child nodes.
    grid_tree_base_node<LEAFDATA, SL>* children[SL * SL * SL];

    // Leaf data pointer. It is NULL if this node is not a leaf node.
    // It is up to the specialized class of grid_tree_base to create this data.
    LEAFDATA* leafData;

    grid_tree_base_node() {
        for( int i = 0; i < SL * SL * SL; ++i )
            children[i] = NULL;
        leafData = NULL;
    }

    ~grid_tree_base_node() {
        for( size_t i = 0; i < SL * SL * SL; ++i )
            if( children[i] )
                delete children[i];
        if( leafData )
            delete leafData;
    }
};

/**
 * grid_tree_base
 * This is a base class for any 3d grid tree object (for example the particle_grid_tree).
 * It is templated on the leaf data class and side length. The functionality this class
 * provides is traversing to, and creating new leaf nodes in a tree representing a 3d area.
 */
template <class LEAFDATA, int SL>
class grid_tree_base {
  private:
    // tree quadrant roots.
    // the grid tree is separated into 8 quadrants due to the way we divide the tree up hierachically.
    // see the wiki page for more information.
    grid_tree_base_node<LEAFDATA, SL>* m_rootNodeQuadrants[8];

  protected:
    /**
     * Default constructor.
     */
    grid_tree_base();

    /**
     * Destructor. Frees all the children and leaf data.
     */
    ~grid_tree_base();

    /**
     * Frees all the children and leaf data.
     */
    void clear();

    /**
     * Swaps the tree with a different tree.
     * @param  rhs  The tree to swap with.
     */
    void swap( grid_tree_base<LEAFDATA, SL>& rhs );

    /**
     * Tree traversal function.
     * @return  The leaf node for a given voxel coordinate, or NULL if it does not exist.
     */
    grid_tree_base_node<LEAFDATA, SL>* navigate_to_leaf( const frantic::graphics::vector3& leafCoord ) const;

    /**
     * Tree traversal function. It will create the node and add it to the tree if it is not already there.
     * @return  The leaf node for a given voxel coordinate.
     */
    grid_tree_base_node<LEAFDATA, SL>* navigate_to_leaf_with_create( const frantic::graphics::vector3& leafCoord );

    /**
     * Tree traversal function. It will return the deepest node in the tree for which the given voxel
     * is a child. The return value can be a leaf node (length == 1), an internal node, or NULL if the
     * query point is outside the tree bounds.
     * @return The deepest node containing the given coordinate, or NULL if not possible.
     */
    grid_tree_base_node<LEAFDATA, SL>* navigate_to_deepest_node( const frantic::graphics::vector3& leafCoord ) const;

  private:
    /// Internal function called by navigate_to_leaf and navigate_to_leaf_with_create. This function does the actual
    /// work of tree traversal.
    grid_tree_base_node<LEAFDATA, SL>* navigate_to_leaf_internal( const frantic::graphics::vector3& leafCoord,
                                                                  bool createIfNeeded );
};

//
//
// grid_tree_base implementation
//
//

template <class LEAFDATA, int SL>
grid_tree_base<LEAFDATA, SL>::grid_tree_base() {
    for( int i = 0; i < 8; ++i )
        m_rootNodeQuadrants[i] = NULL;
}

template <class LEAFDATA, int SL>
grid_tree_base<LEAFDATA, SL>::~grid_tree_base() {
    clear();
}

template <class LEAFDATA, int SL>
void grid_tree_base<LEAFDATA, SL>::clear() {
    // this will delete all the nodes in the tree
    // this function is called on destructor
    for( int i = 0; i < 8; ++i ) {
        if( m_rootNodeQuadrants[i] ) {
            delete m_rootNodeQuadrants[i];
            m_rootNodeQuadrants[i] = NULL;
        }
    }
}

template <class LEAFDATA, int SL>
void grid_tree_base<LEAFDATA, SL>::swap( grid_tree_base<LEAFDATA, SL>& rhs ) {
    for( int i = 0; i < 8; ++i )
        std::swap( m_rootNodeQuadrants[i], rhs.m_rootNodeQuadrants[i] );
}

template <class LEAFDATA, int SL>
grid_tree_base_node<LEAFDATA, SL>*
grid_tree_base<LEAFDATA, SL>::navigate_to_leaf( const frantic::graphics::vector3& leafCoord ) const {
    // call the internal tranversal function.
    // we use a const_cast because the one function does both the "with node create" and "without" versions based on the
    // flag. return const_cast< grid_tree_base<LEAFDATA,SL>* >( this )->navigate_to_leaf_internal( leafCoord, false );
    grid_tree_base_node<LEAFDATA, SL>* pNode =
        const_cast<grid_tree_base<LEAFDATA, SL>*>( this )->navigate_to_leaf_internal( leafCoord, false );
    if( !pNode || pNode->length != 1 )
        return NULL;
    return pNode;
}

template <class LEAFDATA, int SL>
grid_tree_base_node<LEAFDATA, SL>*
grid_tree_base<LEAFDATA, SL>::navigate_to_deepest_node( const frantic::graphics::vector3& leafCoord ) const {
    return const_cast<grid_tree_base<LEAFDATA, SL>*>( this )->navigate_to_leaf_internal( leafCoord, false );
}

template <class LEAFDATA, int SL>
grid_tree_base_node<LEAFDATA, SL>*
grid_tree_base<LEAFDATA, SL>::navigate_to_leaf_with_create( const frantic::graphics::vector3& leafCoord ) {
    // call the internal tranversal function.
    return navigate_to_leaf_internal( leafCoord, true );
}

namespace detail {

template <int SL, class Enable = void>
class grid_tree_base_traversal {
  public:
    grid_tree_base_traversal( boost::uint32_t parentLength )
        : m_parentLength( parentLength ) {}

    boost::uint32_t get_parent_length() const { return m_parentLength; }

    int get_child_node_index( int leafCoord ) const {
        return (int)( ( leafCoord % m_parentLength ) * SL / m_parentLength );
    }

    void step() { m_parentLength /= SL; }

    int get_parent_coord_offset( int childNodeIndex ) const { return childNodeIndex * m_parentLength; }

  private:
    boost::uint32_t m_parentLength;
};

// Special case for when SL is a power of two, so we can use bitwise
// operations instead of / and % when working with parentLength
template <int SL>
class grid_tree_base_traversal<SL, typename boost::enable_if_c<( SL & ( SL - 1 ) ) == 0>::type> {
  public:
    grid_tree_base_traversal( boost::uint32_t parentLength )
        : m_parentLength( parentLength ) {

        if( !frantic::math::is_power_of_two( parentLength ) ) {
            throw frantic::exception_stream() << "grid_tree_base_traversal Internal Error: "
                                              << "parent length (" << parentLength << ") is not a power of two.";
        }

        m_parentLengthBitShift = frantic::math::log2_uint32( m_parentLength );
    }

    boost::uint32_t get_parent_length() const { return m_parentLength; }

    int get_child_node_index( int leafCoord ) const {
        // Note that ( leafCoord & ( parentLength - 1 ) ) == ( leafCoord % parentLength )
        // given that parentLength is a power of two.
        return (int)( ( leafCoord & ( m_parentLength - 1 ) ) * SL >> m_parentLengthBitShift );
    }

    void step() {
        m_parentLengthBitShift -= boost::static_log2<SL>::value;
        m_parentLength /= SL;
    }

    int get_parent_coord_offset( int childNodeIndex ) const { return childNodeIndex << m_parentLengthBitShift; }

  private:
    boost::uint32_t m_parentLength;

    // How many bits do we need to shift by to multiply or divide by the
    // parentLength?
    boost::uint32_t m_parentLengthBitShift;
};

} // namespace detail

template <class LEAFDATA, int SL>
grid_tree_base_node<LEAFDATA, SL>*
grid_tree_base<LEAFDATA, SL>::navigate_to_leaf_internal( const frantic::graphics::vector3& inputLeafCoord,
                                                         bool createIfNeeded ) {
    // determine which quadrant this leaf will go in
    bool xQuad = inputLeafCoord.x < 0;
    bool yQuad = inputLeafCoord.y < 0;
    bool zQuad = inputLeafCoord.z < 0;

    // the index into our m_rootNodeQuadrants array that this leaf belongs to
    int rootNodeIndex = ( ( xQuad ) ? 1 : 0 ) + ( ( yQuad ) ? 2 : 0 ) + ( ( zQuad ) ? 4 : 0 );

    // get the root node for this quadrant
    grid_tree_base_node<LEAFDATA, SL>* rootNode = m_rootNodeQuadrants[rootNodeIndex];

    // handle empty tree
    if( !rootNode ) {
        // create new root node if tree is empty
        if( createIfNeeded ) {
            rootNode = new grid_tree_base_node<LEAFDATA, SL>;
            rootNode->coord = inputLeafCoord;
            rootNode->length = 1;
            m_rootNodeQuadrants[rootNodeIndex] = rootNode;
        }
        return rootNode;
    }

    // get the leaf coordinate in the positive quadrant (for tree navigation purposes)
    boost::uint32_t leafX = ( boost::uint32_t )( ( xQuad ) ? ~inputLeafCoord.x : inputLeafCoord.x );
    boost::uint32_t leafY = ( boost::uint32_t )( ( yQuad ) ? ~inputLeafCoord.y : inputLeafCoord.y );
    boost::uint32_t leafZ = ( boost::uint32_t )( ( zQuad ) ? ~inputLeafCoord.z : inputLeafCoord.z );

    // get the parent coordiante in the positive quadrant (for tree navigation purposes)
    boost::uint32_t parentX = ( boost::uint32_t )( ( xQuad ) ? ~rootNode->coord.x : rootNode->coord.x );
    boost::uint32_t parentY = ( boost::uint32_t )( ( yQuad ) ? ~rootNode->coord.y : rootNode->coord.y );
    boost::uint32_t parentZ = ( boost::uint32_t )( ( zQuad ) ? ~rootNode->coord.z : rootNode->coord.z );
    boost::uint32_t parentLength = rootNode->length;

    // reparent the root node if need be
    if( createIfNeeded ) {

        // this loop creates root nodes until the root node encapsulates our new leaf node
        while( leafX < parentX || leafX > parentX + parentLength - 1 || leafY < parentY ||
               leafY > parentY + parentLength - 1 || leafZ < parentZ || leafZ > parentZ + parentLength - 1 )

        {
            // reparent the tree
            boost::uint32_t newParentLength = parentLength * SL;
            boost::uint32_t newParentX = ( parentX / newParentLength ) * newParentLength;
            boost::uint32_t newParentY = ( parentY / newParentLength ) * newParentLength;
            boost::uint32_t newParentZ = ( parentZ / newParentLength ) * newParentLength;

            // link new root node to our old root
            int childNodeIndexX = (int)( ( parentX - newParentX ) / parentLength );
            int childNodeIndexY = (int)( ( parentY - newParentY ) / parentLength );
            int childNodeIndexZ = (int)( ( parentZ - newParentZ ) / parentLength );

            // create the new parent and put its coordinates in its own correct quadrant
            grid_tree_base_node<LEAFDATA, SL>* parent = new grid_tree_base_node<LEAFDATA, SL>;
            parent->coord.x = ( boost::int32_t )( ( xQuad ) ? ~newParentX : newParentX );
            parent->coord.y = ( boost::int32_t )( ( yQuad ) ? ~newParentY : newParentY );
            parent->coord.z = ( boost::int32_t )( ( zQuad ) ? ~newParentZ : newParentZ );
            parent->length = newParentLength;
            parent->children[childNodeIndexX + ( childNodeIndexY + childNodeIndexZ * SL ) * SL] = rootNode;

            // update loop varibles for our next iteration
            parentX = newParentX;
            parentY = newParentY;
            parentZ = newParentZ;
            parentLength = newParentLength;

            // reroot the tree
            rootNode = parent;
            m_rootNodeQuadrants[rootNodeIndex] = rootNode;
        }

    } else { //! createIfNeeded

        // check if this is in the bounds of our tree if we are not to reparent
        if( leafX < parentX || leafX > parentX + parentLength - 1 || leafY < parentY ||
            leafY > parentY + parentLength - 1 || leafZ < parentZ || leafZ > parentZ + parentLength - 1 ) {
            // this leaf is out of the bounds of our tree
            return NULL;
        }
    }

    // traverse the tree to our node.
    grid_tree_base_node<LEAFDATA, SL>* currNode = rootNode;

    detail::grid_tree_base_traversal<SL> traversal( parentLength );

    while( traversal.get_parent_length() > 1 ) {

        // get this node's appropriate child node
        int childNodeIndexX = traversal.get_child_node_index( leafX );
        int childNodeIndexY = traversal.get_child_node_index( leafY );
        int childNodeIndexZ = traversal.get_child_node_index( leafZ );

        int index = childNodeIndexX + ( childNodeIndexY + childNodeIndexZ * SL ) * SL;
        grid_tree_base_node<LEAFDATA, SL>* newNode = currNode->children[index];

        // update parentCoord and parentLength
        traversal.step();
        parentX += traversal.get_parent_coord_offset( childNodeIndexX );
        parentY += traversal.get_parent_coord_offset( childNodeIndexY );
        parentZ += traversal.get_parent_coord_offset( childNodeIndexZ );

        // if no child is found
        if( !newNode ) {

            if( !createIfNeeded )
                return currNode;
            // return NULL;

            // create a new node and link it up to our old parent
            // and put its coordinates in its own correct quadrant
            newNode = new grid_tree_base_node<LEAFDATA, SL>;
            newNode->coord.x = ( boost::int32_t )( ( xQuad ) ? ~parentX : parentX );
            newNode->coord.y = ( boost::int32_t )( ( yQuad ) ? ~parentY : parentY );
            newNode->coord.z = ( boost::int32_t )( ( zQuad ) ? ~parentZ : parentZ );
            newNode->length = traversal.get_parent_length();
            currNode->children[index] = newNode;
        }

        currNode = newNode;
    }

    return currNode;
}

} // namespace volumetrics
} // namespace frantic
