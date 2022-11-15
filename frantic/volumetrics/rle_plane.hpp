// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <cstdlib>
#include <frantic/graphics2d/boundrect2.hpp>

namespace frantic {
namespace volumetrics {

class rle_plane {

    std::vector<std::pair<int, int>> m_runs;     // pairs of run indices
    std::vector<int> m_runCodes;                 // codes for each run
    std::vector<int> m_xStarts;                  // index into m_runs for the first run in each x extent
    frantic::graphics2d::boundrect2 m_xyExtents; // the plane bounds

  public:
    /**
     *	Constructor takes the xyextents for the plane to be run length encoded.
     *
     *	@param	xyExtents	A boundrect2 of the min/max points on the plane.
     */
    rle_plane( const frantic::graphics2d::boundrect2& xyExtents )
        : m_xyExtents( xyExtents ) {}

    /**
     *	Dumps the contents of the rle plane to an outstream.
     *
     *	@param	fout	The ostream to dump to.
     */
    void dump( std::ostream& fout ) const {
        fout << std::endl << "DUMP of rle_plane" << std::endl;
        fout << "xyExtents: " << m_xyExtents.str() << std::endl;
        fout << "xyExtents.xsize: " << m_xyExtents.xsize() << "  xyExtents.ysize: " << m_xyExtents.ysize() << std::endl;
        fout << "number of runs: " << m_runs.size() << std::endl;
        fout << "number of defined extents: " << m_xStarts.size() << std::endl;
        fout << "Defined runs, codes:";
        for( int i = 0; i < (int)m_xStarts.size(); ++i ) {
            fout << std::endl << "runs for y extent " << i << ":" << std::endl;

            int start = get_y_extent_first_run( i );
            int end = get_y_extent_last_run( i );

            if( start < 0 ) {
                fout << " code: " << start;
                continue;
            }

            std::cout << " first run index: " << start << "  last run index: " << end << " code: " << m_xStarts[i]
                      << std::endl
                      << " ";
            for( int j = start; j <= end; ++j )
                fout << "([" << m_runs[j].first << " " << m_runs[j].second << "], " << m_runCodes[j] << "), ";
        }
        fout << std::endl << "FINISHED DUMP" << std::endl;
    }

    /**
     *	Swaps the contents of the rle plane with the given rle plane.
     *	@param	rlp	The rle plane to swap contents with.
     */
    void swap( rle_plane& rlp ) {
        m_runs.swap( rlp.m_runs );
        m_runCodes.swap( rlp.m_runCodes );
        m_xStarts.swap( rlp.m_xStarts );
        frantic::graphics2d::boundrect2 tempExtents( rlp.m_xyExtents );
        rlp.m_xyExtents = m_xyExtents;
        m_xyExtents = tempExtents;
    }

    /**
     *	Checks if the rle plane has any run data in it.
     *
     *	@return		True if there is any run data, false if there isn't.
     */
    bool is_empty() const { return ( m_runs.empty() ); }

    /**
     *	Gets the "size" of the rle plane.
     *
     *	@return		The number of run index pairs in the rle plane.
     */
    size_t size() const { return m_runs.size(); }

    /**
     *	Clears the run information from the rle plane.  Leaves the xyExtents intact.
     */
    void clear() {
        m_runs.clear();
        m_runCodes.clear();
        m_xStarts.clear();
    }

    /**
     *	Resets the plane, clearing the run information and resetting the xyextents.
     *	@param	xyExtents	The new xyExtents for the plane.
     */
    void reset( const frantic::graphics2d::boundrect2& xyExtents ) {
        m_runs.clear();
        m_runCodes.clear();
        m_xStarts.clear();
        m_xyExtents = xyExtents;
    }

    /**
     *	Resets the plane, clearing the run information.
     */
    void reset() {
        m_runs.clear();
        m_runCodes.clear();
        m_xStarts.clear();
    }
    /**
     *	Gets the xyExtents of the rle plane.
     *
     *	@return		A boundrect2 holding the max/min extents of the rle plane.
     */
    frantic::graphics2d::boundrect2 get_xy_extents() const { return m_xyExtents; }

    /**
     *	Gets the i-th run in the rle plane.
     *
     *	@param	i	The index to retrieve a run from.
     *	@return		A pair of run indices.
     */
    std::pair<int, int> operator[]( const int i ) const { return m_runs[i]; }

    /**
     *	Gets the i-th run in the rle plane.
     *
     *	@param	i	The index to retrieve a run from.
     *	@return		A pair of run indices.
     */
    std::pair<int, int> get_run( const int i ) const { return m_runs[i]; }

    /**
     *	Gets the code associated with the i-th run in the rle plane.
     *
     *	@param	i	The index to retrieve a run code from.
     *	@return		An integer code associated with the run.
     */
    int get_run_code( const int i ) const { return m_runCodes[i]; }

    /**
     *	Gets the index of the first run of the given x extent in the rle plane.  If there are no defined
     *  runs for the requested extent, returns -1.
     *
     *	@param	index	The x extent that you want the first run of.
     *	@return			An the index in the rle plane that contains the first run of the given x extent.
     */
    int get_y_extent_first_run( const int index ) const {
        if( index >= (int)m_xyExtents.ysize() )
            throw std::runtime_error( "rle_plane::get_y_extent_first_run() - Extent index " +
                                      boost::lexical_cast<std::string>( index ) +
                                      " is larger than the largest y extent of the plane (" +
                                      boost::lexical_cast<std::string>( m_xyExtents.ysize() ) + ")" );
        if( index >= (int)m_xStarts.size() )
            throw std::runtime_error( "rle_plane::get_y_extent_first_run() - Extent index " +
                                      boost::lexical_cast<std::string>( index ) +
                                      " is larger than the number of y extents currently in the plane (" +
                                      boost::lexical_cast<std::string>( m_xStarts.size() ) + ")" );

        return m_xStarts[index];
    }

    /**
     *	Gets the index of the last run of the given x extent in the rle plane.  If there are no defined
     *  runs for the requested extent, returns -1.
     *
     *	@param	index	The x extent that you want the last run of.
     *	@return			An the index in the rle plane that contains the last run of the given x extent.
     */
    int get_y_extent_last_run( int index ) const {
        if( index >= (int)m_xyExtents.ysize() )
            throw std::runtime_error( "rle_plane::get_y_extent_first_run() - Extent index " +
                                      boost::lexical_cast<std::string>( index ) +
                                      " is larger than the number of extents in the plane (" +
                                      boost::lexical_cast<std::string>( m_xyExtents.ysize() ) + ")" );

        // extent has no defined data
        if( index >= (int)m_xStarts.size() || m_xStarts[index] < 0 )
            return m_xStarts[index];

        do
            index++;
        while( index < (int)m_xStarts.size() && m_xStarts[index] < 0 );

        if( index == (int)m_xStarts.size() )
            return (int)m_runs.size() - 1;

        return m_xStarts[index] - 1;
    }

    /**
     *	Appends all the runs for an extent of the plane.  A relatively quick operation since it doesn't
     *	check for consistency of run information when added.  A check_run_consistency function exists to
     *	do a single pass to check all run information.  It can be used after a series of calls to this
     *	function.  If you try to add runs for an extent which is already defined in the structure,
     *  an exception will be thrown.  Run codes for each run will default to 0 when runs are appended
     *	by this function.  If particular run codes are desired, use the other overload of this function.
     *
     *	@param	runs	A vector of run pairs to add.
     *	@param	extent	The extent for which to add runs.
     *	@param	code	The negative code for the extent if it is undefined, defaults to -1
     */
    void append_runs_by_extent( const std::vector<std::pair<int, int>>& runs, const int extent, const int code = -1 ) {

        if( extent >= m_xyExtents.ysize() )
            throw std::runtime_error( "rle_plane::append_runs_by_extent() - Cannot append runs for extent " +
                                      boost::lexical_cast<std::string>( extent ) +
                                      ".  This exceeds the maximum y extent (" +
                                      boost::lexical_cast<std::string>( m_xyExtents.ysize() - 1 ) + ")." );
        if( m_xyExtents.get_area() == 0 )
            throw std::runtime_error(
                "rle_plane::append_runs_by_extent() - Cannot append runs, the rle plane bounds are empty." );
        if( extent < 0 )
            throw std::runtime_error( "rle_plane::append_runs_by_extent() - Cannot append runs for extent " +
                                      boost::lexical_cast<std::string>( extent ) +
                                      ".  The RLE plane indexes extents from 0." );
        if( (int)m_xStarts.size() > extent )
            throw std::runtime_error( "rle_plane::append_runs_by_extent() - Cannot append runs for extent " +
                                      boost::lexical_cast<std::string>( extent ) +
                                      ", the RLE plane already has extents defined up to " +
                                      boost::lexical_cast<std::string>( m_xStarts.size() - 1 ) );

        if( runs.empty() ) {
            if( code >= 0 )
                throw std::runtime_error(
                    "rle_plane::append_runs_by_extent() - Given extent code (" +
                    boost::lexical_cast<std::string>( code ) +
                    ") must be negative, so that it doesn't interfere with the non-negative indexing." );
            while( (int)m_xStarts.size() <= extent )
                m_xStarts.push_back( code );
            return;
        } else {
            while( (int)m_xStarts.size() < extent )
                m_xStarts.push_back( code );
            m_xStarts.push_back( int( m_runs.size() ) );
            m_runs.resize( m_runs.size() + runs.size() );
            memcpy( &m_runs[m_runs.size() - runs.size()], &runs[0], sizeof( std::pair<int, int> ) * runs.size() );
            m_runCodes.resize( m_runCodes.size() + runs.size() );
            memset( &m_runCodes[m_runCodes.size() - runs.size()], 0, sizeof( int ) * runs.size() );
        }
    }

    /**
     *	Appends all the runs for an extent of the plane.  A relatively quick operation since it doesn't
     *	check for consistency of run information when added.  A check_run_consistency function exists to
     *	do a single pass to check all run information.  It can be used after a series of calls to this
     *	function.  If you try to add runs for an extent which is already defined in the structure,
     *  an exception will be thrown.
     *
     *	@param	runs		A vector of run pairs to add.
     *	@param	runCodes	A vector of codes corresponding to the above runs to be added.
     *	@param	extent		The extent for which to add runs.
     *	@param	code		The negative code for the extent if it is undefined, defaults to -1
     */
    void append_runs_by_extent( const std::vector<std::pair<int, int>>& runs, const std::vector<int>& runCodes,
                                const int extent, const int code = -1 ) {

        if( extent >= m_xyExtents.ysize() )
            throw std::runtime_error( "rle_plane::append_runs_by_extent() - Cannot append runs for extent " +
                                      boost::lexical_cast<std::string>( extent ) +
                                      ".  This exceeds the maximum y extent (" +
                                      boost::lexical_cast<std::string>( m_xyExtents.ysize() - 1 ) + ")." );
        if( m_xyExtents.get_area() == 0 )
            throw std::runtime_error(
                "rle_plane::append_runs_by_extent() - Cannot append runs, the rle plane bounds are empty." );
        if( runs.size() != runCodes.size() )
            throw std::runtime_error(
                "rle_plane::append_runs_by_extent() - Cannot append runs, the length of the runs vector (" +
                boost::lexical_cast<std::string>( runs.size() ) +
                ") differs from the length of the run codes vector (" +
                boost::lexical_cast<std::string>( runCodes.size() ) + ")." );
        if( extent < 0 )
            throw std::runtime_error( "rle_plane::append_runs_by_extent() - Cannot append runs for extent " +
                                      boost::lexical_cast<std::string>( extent ) +
                                      ".  The RLE plane indexes extents from 0." );
        if( (int)m_xStarts.size() > extent )
            throw std::runtime_error( "rle_plane::append_runs_by_extent() - Cannot append runs for extent " +
                                      boost::lexical_cast<std::string>( extent ) +
                                      ", the RLE plane already has extents defined up to " +
                                      boost::lexical_cast<std::string>( m_xStarts.size() - 1 ) );

        if( runs.empty() ) {
            if( code >= 0 )
                throw std::runtime_error(
                    "rle_plane::append_runs_by_extent() - Given extent code (" +
                    boost::lexical_cast<std::string>( code ) +
                    ") must be negative, so that it doesn't interfere with the non-negative indexing." );
            while( (int)m_xStarts.size() <= extent )
                m_xStarts.push_back( code );
            return;
        } else {

            while( (int)m_xStarts.size() < extent )
                m_xStarts.push_back( code );

            m_xStarts.push_back( int( m_runs.size() ) );

            m_runs.resize( m_runs.size() + runs.size() );
            memcpy( &m_runs[m_runs.size() - runs.size()], &runs[0], sizeof( std::pair<int, int> ) * runs.size() );

            m_runCodes.resize( m_runCodes.size() + runCodes.size() );
            memcpy( &m_runCodes[m_runCodes.size() - runCodes.size()], &runCodes[0], sizeof( int ) * runCodes.size() );
        }
    }

    /**
     *	Merges the two defined runs of the given RLE planes into this one.  This assumes that both
     *  planes are consistent, otherwise the merge of the two planes is not likely to result in a
     *  consistent rle plane.  For disagreements in undefined information, the merge function
     *	always takes the lesser of the codes.  The function only works on planes that contain
     *	only defined runs.
     *
     *	@param	rlp		The rle plane to be merged with this one.
     *
    void merge_defined_runs( frantic::volumetrics::rle_plane& firstRLP, frantic::volumetrics::rle_plane& secondRLP){

      // We can't merge planes with differing dimensions
      if ( firstRLP.get_xy_extents().maximum().x != secondRLP.get_xy_extents().maximum().x ||
         firstRLP.get_xy_extents().maximum().y != secondRLP.get_xy_extents().maximum().y ||
         firstRLP.get_xy_extents().minimum().x != secondRLP.get_xy_extents().minimum().x ||
         firstRLP.get_xy_extents().minimum().y != secondRLP.get_xy_extents().minimum().y )
        throw std::runtime_error("rle_plane::merge_defined_runs() - Cannot merge the provided rle planes, because they
    differ in xy extent size.  First plane has extents: " + firstRLP.get_xy_extents().str() + " and second has extents:
    "
    + secondRLP.get_xy_extents().str());

      // If one plane is empty, just make a copy of the other one
      if ( firstRLP.is_empty() ) {
        m_runs = secondRLP.m_runs;
        m_runCodes = secondRLP.m_runCodes;
        m_xStarts = secondRLP.m_xStarts;
        m_xyExtents = secondRLP.m_xyExtents;
        return;
      }
      if ( secondRLP.is_empty() ) {
        m_runs = firstRLP.m_runs;
        m_runCodes = firstRLP.m_runCodes;
        m_xStarts = firstRLP.m_xStarts;
        m_xyExtents = firstRLP.m_xyExtents;
        return;
      }

      // Set up the merged plane's data members
      m_xyExtents = firstRLP.get_xy_extents();
      m_runs.clear();
      m_runCodes.clear();
      m_xStarts.reserve(m_xyExtents.ysize());

      // loop through the runs extent by extent and create the combined rle plane
      int firstRun = 0, secondRun = 0;
      for ( int currentExtent = 0; currentExtent < firstRLP.get_xy_extents().ysize(); ++currentExtent ) {

        //std::cout << "extent: " << currentExtent << std::endl;

        firstRun = firstRLP.get_y_extent_first_run(currentExtent);
        secondRun = secondRLP.get_y_extent_first_run(currentExtent);

        if ( firstRun < 0 && secondRun < 0 ) {
          //std::cout << "empty" << std::endl;
          m_xStarts.push_back(std::min(firstRun,secondRun));
          continue;
        }

        m_xStarts.push_back((int)m_runs.size());

        //std::cout << "beginning run traversal" << std::endl;
        while ( firstRun >= 0 && firstRun <= firstRLP.get_y_extent_last_run(currentExtent) &&
              secondRun >= 0 && secondRun <= secondRLP.get_y_extent_last_run(currentExtent) ) {

          // push on a new run with the smaller start index and set the plane traversal flag
          bool inFirstPlane;
          if ( firstRLP[firstRun].first < secondRLP[secondRun].first ) {
            //std::cout << "starting run in first plane at " << firstRLP[firstRun].first << std::endl;
            inFirstPlane = true;
            m_runs.push_back( std::pair<int,int>(firstRLP[firstRun].first, 0) );
          }
          else {
            //std::cout << "starting run in second plane at " << secondRLP[secondRun].first << std::endl;
            inFirstPlane = false;
            m_runs.push_back( std::pair<int,int>(secondRLP[secondRun].first, 0) );
          }

          // Traverse runs in both planes until you hit a space where they don't overlap.
          bool endRun = false;
          while ( !endRun ) {
            //std::cout << "current run in first plane: " << firstRLP[firstRun].first << "," <<
    firstRLP[firstRun].second
    << std::endl;
            //std::cout << "current run in second plane: " << secondRLP[secondRun].first << "," <<
    secondRLP[secondRun].second  << std::endl; if ( inFirstPlane && firstRLP[firstRun].second <=
    secondRLP[secondRun].second && firstRLP[firstRun].second >= secondRLP[secondRun].first) {
              //std::cout << "moving to run in second plane" << std::endl;
              inFirstPlane = false;
              while ( ++firstRun < firstRLP.size() && firstRLP[firstRun].second <= secondRLP[secondRun].second );
              if ( firstRun == firstRLP.size() )
                endRun = true;
            }
            else if ( !inFirstPlane && secondRLP[secondRun].second <= firstRLP[firstRun].second &&
    secondRLP[secondRun].second >= firstRLP[firstRun].first) {
              //std::cout << "moving to run in first plane" << std::endl;
              inFirstPlane = true;
              while ( ++secondRun < secondRLP.size() && secondRLP[secondRun].second <= firstRLP[firstRun].second );
              if ( secondRun == secondRLP.size() )
                endRun = true;
            }
            else {
              endRun = true;
              if ( inFirstPlane )
                while ( ++secondRun < secondRLP.size() && secondRLP[secondRun].second <= firstRLP[firstRun].second );
              else
                while ( ++firstRun < firstRLP.size() && firstRLP[firstRun].second <= secondRLP[secondRun].second );
            }
          }

          // Wherever I am right now, I have to stop the run
          if ( inFirstPlane ) {
            //std::cout << "ending run in first plane at " << firstRLP[firstRun].second << std::endl;
            m_runs[m_runs.size()-1].second = firstRLP[firstRun++].second;
          }
          else {
            //std::cout << "ending run in second plane at " << secondRLP[secondRun].second << std::endl;
            m_runs[m_runs.size()-1].second = secondRLP[secondRun++].second;
          }
        }

        // we've run out of runs in at least one extent.  copy the remaining runs over
        // for this extent
        while ( firstRun >= 0 && firstRun <= firstRLP.get_y_extent_last_run(currentExtent) )
          m_runs.push_back(firstRLP[firstRun++]);
        while ( secondRun >= 0 && secondRun <= secondRLP.get_y_extent_last_run(currentExtent) )
          m_runs.push_back(secondRLP[secondRun++]);

      }

      // resize the run code array to match
      m_runCodes.resize(m_runs.size(),0);

    }
*/
    /**
     *	Checks the consistency of the run data currently in the rle plane.  If any run information is out
     *	of order, the function returns false.  Otherwise true.  Prints any errors to the given ostream.
     *
     *	@param		o	An ostream to print errors to.
     *	@return			True if the run information in the plane is all consistent, false if not.
     */
    bool check_consistency( std::ostream& o ) const {
        std::pair<int, int> lastRun( -1, -1 );
        int currentExtent = 0;

        for( int i = 0; i < (int)m_runs.size(); i++ ) {
            while( m_runs[i].first >= currentExtent * m_xyExtents.xsize() )
                currentExtent++;
            if( m_runs[i].first > m_runs[i].second ) { // run length is -ve
                o << "Negative run length at run " << i << "(" << m_runs[i].first << "," << m_runs[i].second << ")"
                  << std::endl;
                return false;
            }
            if( m_runs[i].first <= lastRun.second ) { // out of order/overlapping run
                o << "Start of run " << i << " (" << m_runs[i].first << "," << m_runs[i].second
                  << ") is before end of previous run (" << lastRun.first << "," << lastRun.second << ")" << std::endl;
                return false;
            }
            if( m_runs[i].second >= currentExtent * m_xyExtents.xsize() ) { // run overruns extent
                o << "Run " << i << " (" << m_runs[i].first << "," << m_runs[i].second << ") overruns extent "
                  << currentExtent << " which begins at " << currentExtent * m_xyExtents.xsize() << " and ends at "
                  << ( currentExtent + 1 ) * m_xyExtents.xsize() - 1 << std::endl;
                return false;
            }
            lastRun = m_runs[i];
        }
        return true;
    }

    /**
     *	Checks the consistency of the run data currently in the rle plane.  If any run information is out
     *	of order, the function returns false.  Otherwise true.
     *
     *	@return		True if the run information in the plane is all consistent, false if not.
     */
    bool check_consistency() const {
        std::pair<int, int> lastRun( -1, -1 );
        int currentExtent = 0;

        for( int i = 0; i < (int)m_runs.size(); i++ ) {
            while( m_runs[i].first >= currentExtent * m_xyExtents.xsize() )
                currentExtent++;
            if( m_runs[i].first > m_runs[i].second ) { // run length is -ve
                return false;
            }
            if( m_runs[i].first <= lastRun.second ) { // out of order/overlapping run
                return false;
            }
            if( m_runs[i].second >= currentExtent * m_xyExtents.xsize() ) { // run overruns extent
                return false;
            }
            lastRun = m_runs[i];
        }
        return true;
    }
};

class rle_plane_defined_iterator {
    const rle_plane* m_rlp;
    int m_extentIndex;
    int m_runIndex;
    int m_dataIndex;

  public:
    rle_plane_defined_iterator( const rle_plane* rlp, bool pastTheEnd = false ) {

        if( !pastTheEnd && !rlp->is_empty() ) {

            m_rlp = rlp;
            m_extentIndex = 0;
            m_runIndex = -1;

            while( m_runIndex < 0 && m_extentIndex < rlp->get_xy_extents().ysize() ) {
                m_runIndex = rlp->get_y_extent_first_run( m_extentIndex++ );
            }

            if( m_runIndex >= 0 ) {
                while( m_runIndex < (int)( m_rlp->size() ) && m_rlp->get_run_code( m_runIndex ) < 0 )
                    m_runIndex++;

                if( m_runIndex == (int)( m_rlp->size() ) ) {
                    m_rlp = rlp;
                    m_extentIndex = m_rlp->get_xy_extents().ysize();
                    m_runIndex = (int)m_rlp->size();
                    m_dataIndex = boost::numeric::bounds<int>::highest();
                    return;
                }

                m_dataIndex = ( *m_rlp )[m_runIndex].first;
                --m_extentIndex;
            } else
                m_runIndex = (int)m_rlp->size();
        } else {
            m_rlp = rlp;
            m_extentIndex = m_rlp->get_xy_extents().ysize();
            m_runIndex = (int)m_rlp->size();
            m_dataIndex = boost::numeric::bounds<int>::highest();
        }
    };

    rle_plane_defined_iterator& operator++() {

        if( m_dataIndex == boost::numeric::bounds<int>::highest() )
            return *this;

        if( m_dataIndex >= ( *m_rlp )[m_runIndex].second ) {
            while( ++m_runIndex < (int)m_rlp->size() && m_rlp->get_run_code( m_runIndex ) < 0 )
                ;
            while( m_runIndex >= ( m_extentIndex + 1 ) * m_rlp->get_xy_extents().xsize() )
                m_extentIndex++;

            if( m_runIndex == (int)m_rlp->size() )
                m_dataIndex = boost::numeric::bounds<int>::highest();
            else
                m_dataIndex = ( *m_rlp )[m_runIndex].first;
        } else
            ++m_dataIndex;

        return *this;
    };

    /**
     *	Jumps forward to the given x, if it is defined.  If not, it jumps to the next defined x after it.
     *	If the given x is less than the x the iterator is currently at, it doesn't move.
     */
    void jump_to( int x ) {

        // std::cout << "x: " << x << std::endl;
        // std::cout << "m_runIndex: " << m_runIndex << std::endl;
        // std::cout << "m_dataIndex: " << m_dataIndex << std::endl;
        // std::cout << "m_extentIndex: " << m_extentIndex << std::endl;

        // if we are already past the index or end, don't go anywhere
        if( m_dataIndex == boost::numeric::bounds<int>::highest() || x <= m_dataIndex )
            return;

        // jump by extent first
        // std::cout << "jumping by extent" << std::endl;
        while( m_extentIndex < (int)m_rlp->get_xy_extents().ysize() &&
               x >= ( m_extentIndex + 1 ) * (int)m_rlp->get_xy_extents().xsize() )
            m_extentIndex++;
        if( m_extentIndex >= (int)m_rlp->get_xy_extents().ysize() ) {
            m_dataIndex = boost::numeric::bounds<int>::highest();
            return;
        }

        // jump by runs next
        // std::cout << "jumping by runs" << std::endl;
        while( m_runIndex < (int)m_rlp->size() &&
               ( ( *m_rlp ).get_run_code( m_runIndex ) < 0 || // ignore undefined runs
                 x > ( *m_rlp )[m_runIndex].second ) ) {
            m_runIndex++;
        }
        if( m_runIndex >= (int)m_rlp->size() ) {
            m_dataIndex = boost::numeric::bounds<int>::highest();
            return;
        }

        // If the run we jumped to starts after the requested index, set the iterator to the
        // start of the run.
        // std::cout << "jumping by dataindex" << std::endl;
        if( ( *m_rlp )[m_runIndex].first > x )
            m_dataIndex = ( *m_rlp )[m_runIndex].first;
        // Otherwise, the request index is in the run, so just set the iterator to it
        else
            m_dataIndex = x;
    }

    rle_plane_defined_iterator& operator=( const rle_plane_defined_iterator& iter ) {

        m_rlp = iter.m_rlp;
        m_runIndex = iter.m_runIndex;
        return *this;
    }

    bool operator==( const rle_plane_defined_iterator& iter ) const {
        return m_dataIndex == iter.m_dataIndex && m_runIndex == iter.m_runIndex && m_rlp == iter.m_rlp;
    }

    bool operator!=( const rle_plane_defined_iterator& iter ) const { return !operator==( iter ); }

    int get_extent_index() const { return m_extentIndex; }

    int get_run_index() const { return m_runIndex; }

    int get_data_index() const { return m_dataIndex; }
};

} // namespace volumetrics
}; // namespace frantic
