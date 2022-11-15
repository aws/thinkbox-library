// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#include <cfloat>
#include <limits>

#include <frantic/files/filename_sequence.hpp>

using namespace std;
using namespace boost;
using namespace frantic::files;
using namespace frantic;

namespace frantic {
namespace files {

// constructors
/**
 * frame_set constructor.  This is the default constructor.
 */
frame_set::frame_set() {}

/**
 * frame_set constructor.  This constructor initializes the frame_set to the frame numbers
 * contained in the supplied std::vector<double>.
 *
 * @param	frames	The frame numbers which the frame_set will initially contain.
 */
frame_set::frame_set( const std::vector<double>& frames ) { m_frames.insert( frames.begin(), frames.end() ); }

/**
 * frame_set constructor.  This constructor initializes the frame_set to the frame numbers
 * contained in the supplied std::vector<int>.
 *
 * @param	wholeframes	The frame numbers which the frame_set will initially contain.
 */
frame_set::frame_set( const std::vector<int>& wholeframes ) {
    m_frames.insert( wholeframes.begin(), wholeframes.end() );
}

/**
 * frame_set constructor.  This constructor initializes the frame_set to include all integer
 * frame numbers between startFrame and endFrame, inclusive.
 *
 * @param	startFrame	The starting frame number of the included range.
 * @param	endFrame	The ending frame number of the included range.
 */
frame_set::frame_set( int startFrame, int endFrame ) {
    int start = std::min( startFrame, endFrame ), end = std::max( startFrame, endFrame );
    for( int i = start; i <= end; ++i )
        m_frames.insert( i );
}

frame_set::frame_set( const frame_set& rhs ) { m_frames = rhs.m_frames; }

// set operations
/**
 * Adds the specified integer frame number to the frame_set.  No error occurs if the frame number is
 * already in the frame_set.
 *
 * @param	wholeframe	The frame number to be included in the frame_set.
 */
int frame_set::add_frame( int wholeframe ) {
    m_frames.insert( double( wholeframe ) );
    return wholeframe;
}

/**
 * Adds the specified frame number to the frame_set.  No error occurs if the frame number is
 * already in the frame_set.
 *
 * @param	frame	The frame number to be included in the frame_set.
 */
double frame_set::add_frame( double frame ) {
    m_frames.insert( frame );
    return frame;
}

/**
 * Adds all the integer frame numbers in the supplied std::vector<int> to the frame_set.
 * No error occurs if frames in the vector already exist in the set.
 *
 * @param	wholeframes	The frame numbers to be added to the frame_set.
 */
void frame_set::add_frames( const std::vector<int>& wholeframes ) {
    m_frames.insert( wholeframes.begin(), wholeframes.end() );
}

/**
 * Adds all frame numbers in the supplied std::vector<double> to the frame_set.  No error
 * occurs if frames in the vector already exist in the set.
 *
 * @param	frames	The frame numbers to be added to the frame_set.
 */
void frame_set::add_frames( const std::vector<double>& frames ) { m_frames.insert( frames.begin(), frames.end() ); }

/**
 * Removes the specified integer frame number from the frame_set.  No error occurs if the number
 * is not in the set.
 *
 * @param wholeframe	The frame number to be removed from the frame_set.
 */
int frame_set::remove_frame( int wholeframe ) {
    m_frames.erase( double( wholeframe ) );
    return wholeframe;
}

/**
 * Removes the specified frame number from the frame_set.  No error occurs if the number
 * is not in the set.
 *
 * @param	frame	The frame number to be removed from the frame_set.
 */
double frame_set::remove_frame( double frame ) {
    m_frames.erase( frame );
    return frame;
}

/**
 * Removes the integer frame numbers in the supplied vector from the frame_set.  No error occurs
 * if numbers in the vector do not exist in the set.
 *
 * @param	wholeframes	The integer frame numbers to be removed from the frame_set.
 */
void frame_set::remove_frames( const std::vector<int>& wholeframes ) {
    for( vector<int>::size_type i = 0; i < wholeframes.size(); ++i )
        m_frames.erase( double( wholeframes[i] ) );
}

/**
 * Removes the frame numbers in the supplied vector from the frame_set.  No error occurs
 * if numbers in the vector do not exist in the set.
 *
 * @param	frames	The frame numbers to be removed from the frame_set.
 */
void frame_set::remove_frames( const std::vector<double>& frames ) {
    for( vector<double>::size_type i = 0; i < frames.size(); ++i )
        m_frames.erase( frames[i] );
}

/**
 * Clears all frame numbers from the frame_set.
 */
void frame_set::clear() { m_frames.clear(); }

// set inspection
/**
 * Indicates whether the frame_set is empty.
 *
 * @return	Returns true if the frame_set is empty.
 */
bool frame_set::empty() const { return m_frames.empty(); }

/**
 * Indicates whether the frame_set contains any (non-integer) subframe numbers.
 *
 * @return	Returns true if the frame_set contains (non-integer) subframe numbers.
 */
bool frame_set::has_subframes() const {
    for( set<double>::const_iterator itr = m_frames.begin(); itr != m_frames.end(); ++itr )
        if( floor( *itr ) != *itr )
            return true;

    return false;
}

/**
 * Indicates whether the frame_set contains any integer frame numbers.
 *
 * @return Returns true if the frame_set contains any integer frame numbers.
 */
bool frame_set::has_wholeframes() const {
    for( set<double>::const_iterator itr = m_frames.begin(); itr != m_frames.end(); ++itr )
        if( floor( *itr ) == *itr )
            return true;

    return false;
}

/**
 * Indicates whether the specified frame number exists in the frame_set.
 *
 * @param	frame	A frame number to check for existence in the frame_set.
 * @return	Returns true if the frame number exists in the frame_set.
 */
bool frame_set::frame_exists( double frame ) const {
    // const set<double>::iterator itr = m_frames.find( frame );
    return ( m_frames.find( frame ) != m_frames.end() );
}

/**
 * Indicates the number of frame numbers contained in the frame_set.
 *
 * @return	The number of frame numbers contained in the frame_set.
 */
frantic::files::frame_set::size_type frame_set::size() const { return m_frames.size(); }

/**
 * Indicates the range of all frame numbers in the frame_set.
 *
 * @return	A std::pair that contains the lowest and highest frame numbers in the frame_set.
 */
std::pair<double, double> frame_set::allframe_range() const {
    if( m_frames.size() == 0 )
        throw std::runtime_error( "frame_set::allframe_range():  The frame set is empty." );

    return std::pair<double, double>( *( m_frames.begin() ), *( m_frames.rbegin() ) );
}

/**
 * Indicates the range of subrame (non-integer) framenumbers in the frame_set.
 *
 * @return A std::pair that contains the lowers and highest subframe (non-integer) numbers in the frame_set.
 */
std::pair<double, double> frame_set::subframe_range() const {
    if( m_frames.size() == 0 )
        throw std::runtime_error( "frame_set::subframe_range():  The frame set is empty." );

    double front = -DBL_MAX;

    for( set<double>::const_iterator itr = m_frames.begin(); itr != m_frames.end(); ++itr )
        if( floor( *itr ) != *itr ) {
            front = *itr;
            break;
        }

    if( front == -DBL_MAX )
        throw std::runtime_error( "frame_set::subframe_range():  The frame set does not contain any subframes." );

    double back = front;

    for( set<double>::const_reverse_iterator ritr = m_frames.rbegin(); ritr != m_frames.rend(); ++ritr )
        if( floor( *ritr ) != *ritr ) {
            back = *ritr;
            break;
        }

    return std::pair<double, double>( front, back );
}

/**
 * Indicates the range of subrame (non-integer) framenumbers in the frame_set.
 *
 * @return A std::pair that contains the lowers and highest subframe (non-integer) numbers in the frame_set.
 */
std::pair<int, int> frame_set::wholeframe_range() const {
    if( m_frames.size() == 0 )
        throw std::runtime_error( "frame_set::wholeframe_range():  The frame set is empty." );

    int front = -INT_MAX;

    for( set<double>::const_iterator itr = m_frames.begin(); itr != m_frames.end(); ++itr )
        if( floor( *itr ) == *itr ) {
            front = int( *itr );
            break;
        }

    if( front == -INT_MAX )
        throw std::runtime_error( "frame_set::wholeframe_range():  The frame set does not contain any whole frames." );

    int back = front;

    for( set<double>::const_reverse_iterator ritr = m_frames.rbegin(); ritr != m_frames.rend(); ++ritr )
        if( floor( *ritr ) == *ritr ) {
            back = int( *ritr );
            break;
        }

    return std::pair<int, int>( front, back );
}

/**
 * Finds the missing integer frame numbers between the lowest and highest integer frame numbers in the frame_set.
 *
 * @param	missingWholeFrames	This out parameter will hold the missing integer frame numbers in the frame_set.
 */
void frame_set::missing_wholeframes( std::vector<int>& missingWholeFrames ) const {
    std::pair<int, int> range;

    try {
        range = wholeframe_range();
    } catch( const std::exception& ) {
        throw std::runtime_error(
            "frame_set::missing_wholeframes():  The frame set is empty or does not contain any whole frames." );
    }

    missingWholeFrames.clear();
    if( range.first == range.second )
        return;

    set<double>::const_iterator itr = m_frames.find( range.first ), itr_end = m_frames.find( range.second );
    ++itr_end; // move stop point to one-past range.second
    int checkingFor = int( ( *itr ) + 1 );

    while( itr != itr_end ) {
        if( *itr > checkingFor ) {
            missingWholeFrames.push_back( checkingFor );
            ++checkingFor;
        } else {
            if( *itr == checkingFor )
                ++checkingFor;

            ++itr;
        }
    }
}

/**
 * Finds the missing integer frame numbers in the user-specified, inclusive range.
 *
 * @param	lowFrameNumber		The low frame number in the range to be checked.
 * @param	highFrameNumber		The high frame number in the range to be checked.
 * @param	missingWholeFrames	This out parameter will hold the missing integer frame numbers in the frame_set.
 */
void frame_set::missing_wholeframes( int lowFrameNumber, int highFrameNumber,
                                     std::vector<int>& missingWholeFrames ) const {
    missingWholeFrames.clear();
    int loFrame = std::min( lowFrameNumber, highFrameNumber ), hiFrame = std::max( lowFrameNumber, highFrameNumber );

    if( m_frames.size() == 0 || hiFrame < *( m_frames.begin() ) ||
        loFrame > *( m_frames.rbegin() ) ) { // all frames in specified range are missing
        missingWholeFrames.reserve( hiFrame - loFrame + 1 );
        for( int i = loFrame; i <= hiFrame; ++i )
            missingWholeFrames.push_back( i );

        return;
    }

    set<double>::const_iterator itr = m_frames.begin(), itr_end = m_frames.end();

    int checkingFor = loFrame;
    // add missing frames before set
    for( ; checkingFor < *itr; ++checkingFor )
        missingWholeFrames.push_back( checkingFor );

    // add missing frames within set
    while( itr != itr_end && *itr <= hiFrame ) {
        if( *itr > checkingFor ) {
            missingWholeFrames.push_back( checkingFor );
            ++checkingFor;
        } else {
            if( *itr == checkingFor )
                ++checkingFor;

            ++itr;
        }
    }

    // add missing frames after set
    for( ; checkingFor <= hiFrame; ++checkingFor )
        missingWholeFrames.push_back( checkingFor );

    return;
}

/**
 * Indicates whether the whole frame numbers are continguous between the lowest and highest integer frame numbers
 * in the frame_set.  By convention, this function returns false if the frame_set is empty.
 *
 * @return	Returns true if the whole frame numbers are continguous between the lowest and highest integer
 * frame numbers in the frame_set.
 */
bool frame_set::contiguous() const {
    if( m_frames.empty() )
        return false;

    vector<int> missingFrames;
    try {
        missing_wholeframes( missingFrames );
    } catch( const std::exception& ) {
        return false;
    }

    return ( missingFrames.size() == 0 ); // the set is continguous if there are no missing whole frames.
}

/**
 * Indicates whether the whole frame numbers are continguous in the user-specified, inclusive range.  By convention,
 * this function returns false if the frame_set is empty.
 *
 * @param	lowFrameNumber	The low frame number in the range to be checked.
 * @param	highFrameNumber	The high frame number in the range to be checked.
 * @return	Returns true the whole frame numbers are continguous in the specified range.
 */
bool frame_set::contiguous( int lowFrameNumber, int highFrameNumber ) const {
    if( m_frames.empty() )
        return false;

    vector<int> missingFrames;
    try {
        missing_wholeframes( lowFrameNumber, highFrameNumber, missingFrames );
    } catch( const std::exception& ) {
        return false;
    }

    return ( missingFrames.size() == 0 ); // the set is continguous over the range if there are no missing whole frames.
}

// retrieval
/**
 * Retrieves all the frame numbers in the frame_set.
 *
 * @param	allframeNumbers	This out parameter will hold the frame numbers.
 */
void frame_set::allframe_numbers( std::vector<double>& allframeNumbers ) const {
    allframeNumbers.clear();
    allframeNumbers.insert( allframeNumbers.begin(), m_frames.begin(), m_frames.end() );
}

/**
 * Retrieves all subframe(non-integer) frame numbers in the frame_set.
 *
 * @param	subframeNumbers	This out parameter will hold the subframe (non-integer) frame numbers.
 */
void frame_set::subframe_numbers( std::vector<double>& subframeNumbers ) const {
    subframeNumbers.clear();
    for( set<double>::const_iterator itr = m_frames.begin(); itr != m_frames.end(); ++itr )
        if( floor( *itr ) != *itr )
            subframeNumbers.push_back( *itr );
}

/**
 * Retrieves all integer frame numbers in the frame_set.
 *
 * @param	wholeframeNumbers	This out parameter will hold the integer frame numbers.
 */
void frame_set::wholeframe_numbers( std::vector<int>& wholeframeNumbers ) const {
    wholeframeNumbers.clear();
    for( set<double>::const_iterator itr = m_frames.begin(); itr != m_frames.end(); ++itr )
        if( floor( *itr ) == *itr )
            wholeframeNumbers.push_back( int( *itr ) );
}

/**
 * Retrieves all integer frame numbers in the frame_set.
 *
 * @param	wholeframeNumbers	This out parameter will hold the integer frame numbers.
 */
void frame_set::wholeframe_numbers( std::vector<double>& wholeframeNumbers ) const {
    wholeframeNumbers.clear();
    for( set<double>::const_iterator itr = m_frames.begin(); itr != m_frames.end(); ++itr )
        if( floor( *itr ) == *itr )
            wholeframeNumbers.push_back( int( *itr ) );
}

/**
 * Retrieves the nearest bracketing whole framenumbers to frame in the frame set.
 * Also returns the fractional alpha value indicating how far from the start of the interval
 * that frame lies.
 * Returns false if unable to find bracketing frames and doesn't touch any of the return data.
 *
 * @param	frame		the frame number to find the bracket for
 * @param	outInterval	the return interval of closest frames surrounding the given frame
 * @param	outAlpha	the alpha value, between 0 and 1, indicating the location of the given frame in the
 * interval
 */
bool frame_set::get_nearest_wholeframe_interval( const double frame, std::pair<double, double>& outInterval,
                                                 float& outAlpha ) const {
    double startFrame = floor( frame );
    double endFrame = ceil( frame );

    std::set<double>::const_iterator end( m_frames.end() );

    if( m_frames.empty() || m_frames.find( startFrame ) == end || m_frames.find( endFrame ) == end )
        return false;

    if( startFrame == frame ) {
        outInterval.first = frame;
        outInterval.second = frame;

        outAlpha = 0.0f;
    } else {
        outInterval.first = startFrame;
        outInterval.second = endFrame;

        outAlpha = float( ( frame - outInterval.first ) / ( outInterval.second - outInterval.first ) );
    }

    return true;
}

/**
 * Retrieves the nearest bracketing framenumbers to frame in the frame set, including subframes.
 * Also returns the fractional alpha value indicating how far from the start of the interval
 * that frame lies.
 * Returns false if unable to find bracketing frames and doesn't touch any of the return data.
 *
 * @param	frame		the frame number to find the bracket for
 * @param	outInterval	the return interval of closest frames surrounding the given frame
 * @param	outAlpha	the alpha value, between 0 and 1, indicating the location of the given frame in the
 * interval
 */
bool frame_set::get_nearest_subframe_interval( const double frame, std::pair<double, double>& outInterval,
                                               float& outAlpha ) const {
    set<double>::const_iterator intervalEnd = m_frames.lower_bound( frame );

    if( m_frames.empty() || ( intervalEnd == m_frames.begin() && frame < *intervalEnd ) ||
        intervalEnd == m_frames.end() )
        return false;

    if( *intervalEnd > frame ) {
        outInterval.second = *intervalEnd;
        outInterval.first = *( --intervalEnd );
        outAlpha = float( ( frame - outInterval.first ) / ( outInterval.second - outInterval.first ) );
        return true;
    } else {
        outInterval.first = frame;
        outInterval.second = frame;
        outAlpha = 0.0f;
        return true;
    }
}

/**
 * Retrieves the nearest whole framenumber to frame in the frame set.
 *
 * Return false if unable to find the frame.
 *
 * @param	givenFrame		the frame number to find
 * @param	outFrame		the return frame found at the rounded given frame
 */
bool frame_set::get_nearest_wholeframe( const double givenFrame, double& outFrame ) const {
    set<double>::const_iterator itr = m_frames.find( floor( givenFrame + 0.5f ) );

    if( itr == m_frames.end() )
        return false;
    else {
        outFrame = *itr;
        return true;
    }
}

/**
 * Retrieves the nearest subframe in the set to the given frame.
 *
 * Return false if unable to find bracketing frames.
 *
 * @param	givenFrame		the frame number to find
 * @param	outFrame		the return number of closest frames surrounding the given frame
 */
bool frame_set::get_nearest_subframe( const double givenFrame, double& outFrame ) const {
    if( m_frames.size() == 0 )
        return false;

    // lower_bound returns an iterator to the first element in a set with a key that
    // is equal to or greater than a specified key.
    set<double>::const_iterator itr = m_frames.lower_bound( givenFrame );

    // the end of m_frames is a dummy node
    if( itr == m_frames.end() ) {
        outFrame = *( --itr );
        return true;
    }
    double lowerBoundNumber = *itr;

    if( itr == m_frames.begin() )
        outFrame = *itr;
    else {
        --itr;
        // checking either lowerboundnumber or the number right before is nearest to givenFrame
        if( std::abs( *(itr)-givenFrame ) < std::abs( lowerBoundNumber - givenFrame ) )
            outFrame = *itr;
        else
            outFrame = lowerBoundNumber;
    }
    return true;
}

// iterator support
/**
 * @return	Returns a frame_set::iterator that points to the first frame number in the frame_set.
 */
frantic::files::frame_set::iterator frame_set::begin() { return m_frames.begin(); }

/**
 * @return	Returns a frame_set::const_iterator that points to the first frame number in the frame_set.
 */
frantic::files::frame_set::const_iterator frame_set::begin() const { return m_frames.begin(); }

/**
 * @return	Returns a frame_set::iterator that points to one past the last frame number in the frame_set.
 */
frantic::files::frame_set::iterator frame_set::end() { return m_frames.end(); }

/**
 * @return	Returns a frame_set::const_iterator that points to one past the last frame number in the frame_set.
 */
frantic::files::frame_set::const_iterator frame_set::end() const { return m_frames.end(); }

/**
 * @return	Returns a frame_set::reverse_iterator that points to the first frame number in the reverse frame_set.
 */
frantic::files::frame_set::reverse_iterator frame_set::rbegin() { return m_frames.rbegin(); }

/**
 * @return	Returns a frame_set::const_reverse_iterator that points to the first frame number in the reverse
 * frame_set.
 */
frantic::files::frame_set::const_reverse_iterator frame_set::rbegin() const { return m_frames.rbegin(); }

/**
 * @return	Returns a frame_set::reverse_iterator that points to one past the last frame number in the reverse
 * frame_set.
 */
frantic::files::frame_set::reverse_iterator frame_set::rend() { return m_frames.rend(); }

/**
 * @return	Returns a frame_set::const_reverse_iterator that points to one past the last frame number in the reverse
 * frame_set.
 */
frantic::files::frame_set::const_reverse_iterator frame_set::rend() const { return m_frames.rend(); }
} // namespace files
} // namespace frantic
