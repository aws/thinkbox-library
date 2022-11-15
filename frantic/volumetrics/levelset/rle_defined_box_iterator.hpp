// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <limits>
#include <vector>

#include <frantic/volumetrics/levelset/rle_general_block_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * This class manages a 3D grid of cursors in the rle index spec, stepping along the X axis
 * in the positive direction.  It provides convenient ways to jump over large swathes of
 * undefined voxels, so as to be able to take advantage of the rle compression. And do access
 * neighboring information efficiently
 */
class rle_defined_box_iterator {
    const rle_index_spec& m_ris;

    // the leading slice of the box
    rle_general_block_iterator m_blockIterator;

    // When the m_box's minimum x reaches this value, we are done.
    boost::int32_t m_finishAtX;

    // maintaining the xsize because its used so often
    boost::uint32_t m_xSize;
    boost::uint32_t m_blockSize;

    // data indices of the box, stored in the order of the box, so m_dataIndices[0] is (xMin,YMin,ZMin) and
    // (xMin,yMax,zMin) would be stored at m_dataIndices[yMax] (Y + z*ysize + x*ysize*zsize ) with x,y,z relative the
    // current box. This is to allow block copying when retreiving the indices from the block primitives

    // I think the above comment is obsolete, and now :
    // ( x + y * xSize + z * xSize * ySize )
    std::vector<boost::int32_t> m_dataIndices;

    // Pre-allocate temporary arrays
    // for the full box:
    std::vector<boost::int32_t> m_tempBoxIndices;
    // for the leading slice of the box:
    std::vector<boost::int32_t> m_tempSliceIndices;

    // the actual current box
    frantic::graphics::boundbox3 m_box;

    // Disable the assignment operator
    rle_defined_box_iterator& operator=( const rle_defined_box_iterator& );

  public:
    /**
     * Constructs a box iterator for an rle index spec using the provided box for dimensions as well as starting
     * location
     */
    rle_defined_box_iterator( const rle_index_spec& ris, const frantic::graphics::boundbox3& startBox )
        : m_ris( ris )
        , m_blockIterator( ris, startBox.minimum().x, startBox.minimum().y, startBox.maximum().y, startBox.minimum().z,
                           startBox.maximum().z )
        , m_finishAtX( std::numeric_limits<boost::int32_t>::max() - startBox.xsize() )
        , m_box( startBox ) {
        if( m_box.is_empty() )
            std::runtime_error( "rle_defined_box_iterator() - The provided box to iterate with is empty." );

        frantic::graphics::vector3& boxMin = m_box.minimum();

        const int xMin = boxMin.x;

        m_xSize = m_box.xsize();
        m_blockSize = m_box.ysize() * m_box.zsize();

        m_dataIndices.resize( m_box.get_volume() );

        m_tempBoxIndices.resize( m_box.get_volume() );
        m_tempSliceIndices.resize( m_blockSize );

        // for each plane get the indices for this x
        for( boost::uint32_t i = 0; i < m_xSize; ++i ) {
            if( i == 0 ) {
                m_blockIterator.jump_to_x( xMin );
                m_blockIterator.get_indexes( &m_tempSliceIndices[0] );
            } else {
                m_blockIterator.increment_x_and_get_indexes( &m_tempSliceIndices[0] );
            }

            for( boost::uint32_t j = 0; j < m_blockSize; ++j ) {
                const boost::int32_t yzOffset = static_cast<boost::int32_t>( m_xSize * j );
                m_dataIndices[i + yzOffset] = m_tempSliceIndices[j];
            }

            if( m_blockIterator.is_finished() ) {
                m_finishAtX = std::min<boost::int32_t>( xMin + i, m_finishAtX );
            }
        }
    }

    /**
     * Returns the current box of the iterator
     */
    const frantic::graphics::boundbox3& current_box() const { return m_box; }

    /**
     * Returns the current minimum coordinate of the iterator
     */
    const frantic::graphics::vector3& current_min() const { return m_box.minimum(); }

    /**
     * Returns the current maximum coordinate of the iterator
     */
    const frantic::graphics::vector3& current_max() const { return m_box.maximum(); }

    /**
     * This function advances the box in the x direction
     */
    rle_defined_box_iterator& operator++() {

        // if we aren't already done, jump the 'front' block to the next plane and reset the front
        if( m_box.minimum().x < m_finishAtX ) {
            // move the box
            ++m_box.minimum().x;
            ++m_box.maximum().x;
            m_blockIterator.increment_x_and_get_indexes( &m_tempSliceIndices[0] );

            if( m_blockIterator.is_finished() ) {
                m_finishAtX = std::min<boost::int32_t>( m_box.maximum().x, m_finishAtX );
            }

            // shuffle forward the current indices

            // Shift forward by 1 in x dimension ( move the data back by 1 ).
            // This places incorrect data in the last x position, which
            // should be overwritten in the loop below.
            if( ( m_dataIndices.size() ) > 1 )
                memcpy( &m_tempBoxIndices[0], &m_dataIndices[1],
                        sizeof( boost::int32_t ) * ( m_dataIndices.size() - 1 ) );

            // silliness to avoid secure scl checks
            // todo: remove this when you're no longer using msvc8/9 for flood
            boost::int32_t* tempBoxIndices = &m_tempBoxIndices[0];
            boost::int32_t* tempSliceIndices = &m_tempSliceIndices[0];

            // fill in the last (y,z) plane
            const int lastx = m_xSize - 1;
            int boxIndex = lastx;
            for( boost::uint32_t i = 0; i < m_blockSize; ++i ) {
                tempBoxIndices[boxIndex] = tempSliceIndices[i];
                boxIndex += m_xSize;
            }

            m_tempBoxIndices.swap( m_dataIndices );
        }

        return *this;
    }

    bool is_xplane_finished( unsigned int xoffset ) {
        if( xoffset > m_xSize ) {
            throw std::runtime_error(
                "is_xplane_finished() - x offset provided is outside of the iterator range in x: " +
                boost::lexical_cast<std::string>( m_xSize ) );
        }

        return ( m_box.minimum().x + static_cast<boost::int32_t>( xoffset ) ) >= m_finishAtX;
    }

    /**
     * Check if this block iterator passed the end.
     *
     * Return true if it passed the end. Return false otherwise.
     *
     */
    bool is_finished() { return m_box.minimum().x >= m_finishAtX; }

    /**
     *  Return a pointer to the current box indices.
     *
     * @note The index array is of length current_box().get_volume().  This
     *		is the same as the volume of the box which was passed to the
     *		constructor.
     *
     * @note The array is indexed using [x + y * xSize + z * xSize * ySize],
     *		where (x, y, z) is the offset from the current_min() corner of
     *		the box.
     *
     * @return a pointer to the current box indices.
     */
    const boost::int32_t* const get_indices() const { return &m_dataIndices[0]; }

    /**
     * This function dumps the current iterator data to an outstream.
     *
     * @param  o  The outstream to dump to.
     */
    void dump( std::ostream& o ) const {
        o << std::endl << "---Dumping RLE Index Spec Box Iterator" << std::endl;
        o << "Current box: " << m_box << std::endl;
        o << "Iteration will end when the box\'s minimum x reaches: " << m_finishAtX << std::endl;
        // o << "Iterators in block: " << m_blockIterators.size() << std::endl;
        // unsigned int index = 0;

        /*for ( index = 0; index < m_cursorBlocks.size(); ++index ) {
          m_cursorBlocks[index]->dump( o );
        }*/

        // o << "Current Front:\t" << m_currentStart << std::endl;
        o << "Data Indices: \n";

        // int xsize = m_box.xsize();
        int ysize = m_box.ysize();
        int zsize = m_box.zsize();
        // reshuffle the indices into an expected permutation for outside use
        for( int x = 0; x < (int)m_xSize; ++x ) {
            for( int z = 0; z < (int)zsize; ++z ) {
                for( int y = 0; y < (int)ysize; ++y ) {
                    // for ( int x = m_currentStart, i=0 ; i < m_xSize ; x = (x+1)% m_xSize, ++i){
                    //  I think the indexing order has changed to this ?
                    o << frantic::graphics::vector3( x, y, z ) + m_box.minimum() << ": "
                      << m_dataIndices[x + y * m_xSize + z * m_xSize * ysize] << std::endl;
                    // o << vector3(x,y,z) + m_box.minimum() << ": " << m_dataIndices[y + z*ysize + x*ysize*zsize] <<
                    // std::endl; o << vector3(i,y,z) << ": " << m_dataIndices[y + z*ysize + x*ysize*zsize] <<
                    // std::endl;
                }
            }
        }
    }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
