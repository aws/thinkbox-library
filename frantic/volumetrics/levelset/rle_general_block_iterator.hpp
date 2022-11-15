// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_run_iterator.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * This class manages a 2D grid of cursors in the YZ plane, stepping along the X axis
 * in the positive direction.  It provides convenient ways to jump over large swathes of
 * undefined voxels, so as to be able to take advantage of the rle compression.
 */
class rle_general_block_iterator {
    const rle_index_spec& m_ris;

    // min and max values of the block on y-axis and z-axis
    int m_yMin, m_yMax, m_zMin, m_zMax;

    // the scanline iterators, data indices and current x positions of each iterator for the block,
    // one for each yz coord
    std::vector<rle_run_iterator> m_iterBlock;
    std::vector<boost::int32_t> m_dataIndices, m_iterOffsetX;
    rle_run_iterator m_end;

    // the current X and next X for the whole block
    int m_currentX;

    // is the iterator past all the run data
    bool m_finished;

    // Disable the assignment operator
    rle_general_block_iterator& operator=( const rle_general_block_iterator& );

    int end_x() const { return ( std::numeric_limits<int>::max )(); }

    /**
     * This function will find the next closest run's starting x value.  It considers both defined and
     * undefined runs.
     */
    int find_next_run() const {

        int minX = ( std::numeric_limits<int>::max )();
        int minIter = -1;
        for( int i = 0; i < (int)m_iterBlock.size(); ++i ) {
            int x = m_iterBlock[i].get_xmax();
            if( x < minX ) {
                minX = x;
                minIter = i;
            }
        }

        if( minIter < 0 )
            return end_x();
        else
            return minX + 1;
    }

    /**
     * This function will find the next x value in the block.  This is either the next defined piece of data, or
     * the next run start, whichever is defined first.
     */
    int find_next_x() const {

        // check if any of the iterators are in a defined run, because if they are, the next x
        // is just m_currentX+1
        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {
            if( m_iterBlock[i].get_data_index() >= 0 && m_iterOffsetX[i] >= 0 &&
                m_iterBlock[i].get_xsize() >= m_iterOffsetX[i] )
                return m_currentX + 1;
        }

        // if the above isn't the case, find the next run start
        int minX = std::numeric_limits<int>::max();
        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {
            if( m_iterBlock[i] != m_end ) {
                if( m_iterOffsetX[i] < 0 && m_iterBlock[i].get_xmin() < minX )
                    minX = m_iterBlock[i].get_xmin();
                else if( m_iterBlock[i].get_xmax() + 1 < minX )
                    minX = m_iterBlock[i].get_xmax() + 1;
            }
        }

        // if you couldn't find any, return the end x value
        if( minX == std::numeric_limits<int>::max() )
            return end_x();
        else
            return minX;
    }

    /**
     * This function will find the next defined x value in the block.
     */
    int find_next_defined_x() const {

        // check if any of the iterators are in a defined run, because if they are, the next defined x
        // is just m_currentX+1
        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {
            if( m_iterBlock[i].get_data_index() >= 0 && m_iterOffsetX[i] >= 0 &&
                m_iterBlock[i].get_xsize() >= m_iterOffsetX[i] )
                return m_currentX + 1;
        }

        // otherwise, make a copy of the iterators and advance them all until they hit a defined x
        std::vector<rle_run_iterator> tempIterBlock( m_iterBlock );
        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {
            while( tempIterBlock[i] != m_end && tempIterBlock[i].get_data_index() < 0 )
                ++tempIterBlock[i];
        }

        // take the smallest one
        int minX = std::numeric_limits<int>::max();
        for( size_t i = 0; i < tempIterBlock.size(); ++i ) {
            int iterX = tempIterBlock[i].get_xmin();
            if( tempIterBlock[i] != m_end && iterX < minX )
                minX = iterX;
        }

        // if you couldn't find any, return the end x value
        if( minX == std::numeric_limits<int>::max() )
            return end_x();
        else
            return minX;
    }

  public:
    /**
     *
     * This function will check if there is indexes for defined data.
     *
     * Deprecated: use has_defined_index() instead.
     *
     */
    bool is_currentX_has_defined_indexes() const {
        std::vector<boost::int32_t> indexes;
        get_indexes( indexes );
        for( unsigned int i = 0; i < indexes.size(); i++ ) {
            if( indexes[i] >= 0 )
                return true;
        }
        return false;
    }

    rle_general_block_iterator( const rle_index_spec& ris, int xStart, int yMin, int yMax, int zMin, int zMax )
        : m_ris( ris )
        , m_yMin( yMin )
        , m_yMax( yMax )
        , m_zMin( zMin )
        , m_zMax( zMax ) {

        if( yMin > yMax || zMin > zMax )
            std::runtime_error(
                "rle_general_block_iterator::init - One or both of the block max coords is less than the "
                "block min coords.  (y_min,z_min): (" +
                boost::lexical_cast<std::string>( yMin ) + "," + boost::lexical_cast<std::string>( zMin ) +
                ")   (y_max, z_max): (" + boost::lexical_cast<std::string>( yMax ) + "," +
                boost::lexical_cast<std::string>( zMax ) + ")" );

        // allocate the required space
        size_t blockSize = ( yMax - yMin + 1 ) * ( zMax - zMin + 1 );
        m_iterBlock.resize( blockSize );
        m_iterOffsetX.resize( blockSize );

        // init all the iterators in the block
        m_finished = true;
        size_t i = 0;
        for( int z = zMin; z <= zMax; ++z ) {
            int c = m_ris.z_to_c( z );
            for( int y = yMin; y <= yMax; ++y ) {

                m_iterBlock[i] = rle_run_iterator( m_ris, m_ris.y_to_b( y ), c );

                // advance the iterator to the start point
                while( m_iterBlock[i] != m_end && m_iterBlock[i].get_xmax() < xStart )
                    ++m_iterBlock[i];

                if( m_iterBlock[i] != m_end )
                    m_finished = false;

                m_iterOffsetX[i] = xStart - m_iterBlock[i].get_xmin();
                ++i;
            } // for y
        }     // for z

        m_currentX = xStart;
    }

    /**
     * This function returns the X value this block iterator is currently sitting at.
     */
    int current_x() const { return m_currentX; }

    /**
     * This function returns the X value this block iterator will step to next.
     */
    int next_x() const { return find_next_x(); }

    /**
     * This function increments the block in the x direction by one.
     *
     * @return	True if the cursor advanced to a new set of data indexes, false if not.
     */
    bool increment_x() {

        if( m_finished )
            return false;

        bool newData = false;
        m_finished = true;
        ++m_currentX;
        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {

            if( m_iterBlock[i] == m_end )
                continue;

            m_finished = false;

            // if we're at the end of the run, bump into the next run and reset the offset
            if( m_iterBlock[i].get_xmax() < m_currentX ) {
                ++m_iterBlock[i];
                m_iterOffsetX[i] = 0;
                newData = true;
            }
            // if we're in a run, just bump up the x offset
            else {
                ++m_iterOffsetX[i];
                // if it's a defined run, or we just entered an undefined run from an exterior
                // undefined region, set the new data flag
                if( m_iterBlock[i].get_data_index() >= 0 || m_iterOffsetX[i] == 0 )
                    newData = true;
            }
        }

        return newData;
    }

    /**
     * This function increments the block in the x direction by one, and gets
     * the block iterator's indices after the increment.
     *
     * @note this is a combination of jump_to_x and get_indexes.
     *
     * @param outIndexes This is where the indexes go.  It must point to vector of ints to hold all the indexes.
     * @return	True if the cursor advanced to a new set of data indexes, false if not.
     */
    bool increment_x_and_get_indexes( boost::int32_t* outIndexes ) {

        if( m_finished ) {
            for( size_t i = 0; i < m_iterBlock.size(); ++i ) {
                outIndexes[i] = m_ris.m_exteriorRegionCode;
            }
            return false;
        }

        // silliness to avoid secure scl checks
        // TODO: get rid of this when you're no longer using msvc8/9 for flood
        rle_run_iterator* iterBlock = m_iterBlock.size() > 0 ? &m_iterBlock[0] : 0;
        boost::int32_t* iterOffsetX = m_iterOffsetX.size() > 0 ? &m_iterOffsetX[0] : 0;

        bool newData = false;
        m_finished = true;
        ++m_currentX;
        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {
            // jump along runs until you find the one containing the x
            if( iterBlock[i] != m_end && iterBlock[i].get_xmax() < m_currentX ) {
                ++iterBlock[i];
                newData = true;
            }

            if( iterBlock[i] == m_end ) {
                outIndexes[i] = m_ris.m_exteriorRegionCode;
                continue;
            }

            m_finished = false;

            // if we are in a defined run, then there is new data.  if we moved into an
            // undefined run from an exterior region, then there's new data.
            const int newOffset = m_currentX - iterBlock[i].get_xmin();
            if( iterBlock[i].get_data_index() >= 0 || iterOffsetX[i] * newOffset <= 0 )
                newData = true;

            // update the offset
            iterOffsetX[i] = newOffset;

            // if the current x is before the start of the first run,
            // then assign the exterior region code
            if( newOffset < 0 )
                outIndexes[i] = m_ris.m_exteriorRegionCode;
            else {
                boost::int32_t newDataIndex = iterBlock[i].get_data_index();
                // if we're in an undefined run, return the undefined code
                // and if we're in a defined run, offset the data index appropriately
                if( newDataIndex >= 0 )
                    newDataIndex += newOffset;
                outIndexes[i] = newDataIndex;
            }
        }
        return newData;
    }

    /**
     * This function advances the block iterator to the requested x.
     *
     * @return True if the cursor advanced to a new set of data indexes, false if not.
     */
    bool jump_to_x( int x ) {

        if( x < m_currentX )
            throw std::runtime_error(
                "rle_general_block_iterator::jump_to_x() - Cannot jump to the given x location, " +
                boost::lexical_cast<std::string>( x ) + ".  The block is already at x location " +
                boost::lexical_cast<std::string>( m_currentX ) + " and jumping backwards is not supported." );

        if( x == m_currentX || m_finished )
            return false;

        // avoid secure scl checks
        // TODO: get rid of this when you're no longer using msvc8/9
        rle_run_iterator* iterBlock = m_iterBlock.size() > 0 ? &m_iterBlock[0] : 0;
        boost::int32_t* iterOffsetX = m_iterOffsetX.size() > 0 ? &m_iterOffsetX[0] : 0;

        bool newData = false;
        m_finished = true;
        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {
            // jump along runs until you find the one containing the x
            while( iterBlock[i] != m_end && iterBlock[i].get_xmax() < x ) {
                ++iterBlock[i];
                newData = true;
            }

            if( iterBlock[i] == m_end )
                continue;

            m_finished = false;

            // if we are in a defined run, then there is new data.  if we moved into an
            // undefined run from an exterior region, then there's new data.
            int newOffset = x - iterBlock[i].get_xmin();
            if( iterBlock[i].get_data_index() >= 0 || iterOffsetX[i] * newOffset <= 0 )
                newData = true;

            // update the offset
            iterOffsetX[i] = newOffset;
        }
        m_currentX = x;
        return newData;
    }

    /**
     * This function advances the block to the next defined x value, or the next run start,
     * whichever is first.
     */
    void jump_to_next_x() { jump_to_x( find_next_x() ); }

    /**
     * This function advances the block to the next defined x value, or the next run start,
     * whichever is first.
     */
    rle_general_block_iterator& operator++() {
        jump_to_x( find_next_x() );
        return *this;
    }

    /**
     * Check if this block iterator passed the end of the defined region of the index spec.
     *
     * Return true if it passed the end. Return false otherwise.
     *
     */
    bool is_finished() const { return m_finished; }

    /**
     * This function fills the parameter outIndexes with the data indexes currently in the block iterator.
     *
     * @param  outIndexes  This is where the indexes go.  It must point to vector of ints to hold all the indexes.
     */
    void get_indexes( boost::int32_t* outIndexes ) const {

        rle_run_iterator end;

        // avoid secure scl checks
        // TODO: get rid of this when you're no longer using msvc8/9
        const rle_run_iterator* iterBlock = m_iterBlock.size() > 0 ? &m_iterBlock[0] : 0;
        const boost::int32_t* iterOffsetX = m_iterOffsetX.size() > 0 ? &m_iterOffsetX[0] : 0;

        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {

            // if this iterator is past the end, or the current x is before the start of the first run,
            // then assign the exterior region code
            if( iterBlock[i] == end || m_currentX < iterBlock[i].get_xmin() || iterOffsetX[i] < 0 )
                outIndexes[i] = m_ris.m_exteriorRegionCode;

            else {
                outIndexes[i] = iterBlock[i].get_data_index();
                // if we're in an undefined run, return the undefined code
                // and if we're in a defined run, offset the data index appropriately
                if( outIndexes[i] >= 0 )
                    outIndexes[i] += iterOffsetX[i];
            }
        }
    }

    /**
     * This function fills the parameter outIndexes with the data indexes currently in the block iterator.
     *
     * @param  outIndexes  This is where the indexes go.  It must point to vector of ints to hold all the indexes.
     */
    void get_indexes( std::vector<boost::int32_t>& outIndexes ) const {

        if( outIndexes.size() != m_iterBlock.size() )
            outIndexes.resize( m_iterBlock.size() );

        if( outIndexes.size() != 0 )
            get_indexes( &outIndexes[0] );
    }

    /**
     *
     * This function will check if there there is at least one index for defined data at the current iterator position.
     *
     */
    bool has_defined_index() const {
        // const bool expectedResult = is_currentX_has_defined_indexes();

        rle_run_iterator end;

        for( size_t i = 0; i < m_iterBlock.size(); ++i ) {

            // if this iterator is past the end, or the current x is before the start of the first run,
            // then assign the exterior region code
            if( m_iterBlock[i] != end && m_currentX >= m_iterBlock[i].get_xmin() && m_iterOffsetX[i] >= 0 ) {
                if( m_iterBlock[i].get_data_index() >= 0 ) {
                    return true;
                }
            }
        }

        return false;
    }

    /**
     * This function dumps the current iterator data to an outstream.
     *
     * @param  o  The outstream to dump to.
     */
    void dump( std::ostream& o ) const {
        o << std::endl << "---Dumping RLE Index Spec Block Iterator X---" << std::endl;
        o << "Current x: " << m_currentX << std::endl;
        o << "Next x: " << find_next_x() << std::endl;
        o << "Iterators in block: " << std::endl;

        std::vector<boost::int32_t> dataIndexes( m_iterBlock.size() );
        get_indexes( &dataIndexes[0] );
        unsigned int i = 0;
        rle_run_iterator end;
        for( int z = m_zMin; z <= m_zMax; ++z ) {
            for( int y = m_yMin; y <= m_yMax; ++y ) {

                o << "\nYZ coord: (" << y << "," << z << ")   ";

                if( m_iterBlock[i] == end )
                    o << "xpos: END   ";
                else
                    o << "xpos: " << m_iterBlock[i].get_xmin() + m_iterOffsetX[i] << "   ";

                o << "runStart: " << m_iterBlock[i].get_xmin() << "   xOffset: " << m_iterOffsetX[i]
                  << "   dataIndex: " << dataIndexes[i];

                ++i;
            }
        }
        o << "\n---Finished Dumping RLE Index Spec Block Iterator X---" << std::endl;
    }
};
} // namespace levelset
} // namespace volumetrics
} // namespace frantic
