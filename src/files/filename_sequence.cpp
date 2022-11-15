// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

//#include <ShlObj.h>
#include <string>

#include <frantic/files/filename_sequence.hpp>
#include <frantic/files/files.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/misc/string_functions.hpp>

using namespace std;
using namespace boost;
using namespace frantic::files;
using namespace frantic;

namespace frantic {
namespace files {

// filename_sequence implementation
// constructors

/**
 * filename_sequence constructor.  This is the default constructor.
 *
 */
filename_sequence::filename_sequence()
    : m_pattern()
    , m_set() {
    ;
}

/**
 * filename_sequence constructor.  For this construtor, the supplied path should contain
 * an example file (i.e., the path cannot be just a folder).  The file portion can be a pattern,
 * containing #### for the sequence numbers.
 *
 * @param	path	An example path on which the pattern will be based.
 */
filename_sequence::filename_sequence( const frantic::tstring& path )
    : m_pattern( path )
    , m_set() {
    ;
}

/**
 * filename_sequence constructor.   For this construtor, the supplied path should contain
 * an example file (i.e., the path cannot be just a folder).  The file portion can be a pattern,
 * containing #### for the sequence numbers.
 *
 * @param	path				An example path on which the pattern will be based.
 * @param	negativeConvention	The negative convention to be used for the pattern.
 */
filename_sequence::filename_sequence( const frantic::tstring& path, int negativeConvention )
    : m_pattern( path, negativeConvention )
    , m_set() {
    ;
}

/**
 * filename_sequence constructor.  This constructor accepts the directory and filename pattern portions separately.
 * The direcotry can be an empty string.  The file portion can be a pattern, containing #### for the sequence numbers.
 *
 * @param	dir				The directory portion of the path.
 * @param	filenamePattern	The filename portion of the path.
 */
filename_sequence::filename_sequence( const frantic::tstring& dir, const frantic::tstring& filenamePattern )
    : m_pattern( frantic::tstring( ensure_trailing_pathseparator( backward_slashes( dir ) ) + filenamePattern ) )
    , m_set() {
    ;
}

/**
 * filename_sequence constructor.  This constructor accepts the directory and filename pattern portions separately.
 * The direcotry can be an empty string.  The file portion can be a pattern, containing #### for the sequence numbers.
 * Additionally, the negative convention can be directly specified.
 *
 * @param	dir				The directory portion of the path.
 * @param	filenamePattern	The filename portion of the path.
 * @param	negativeConvention	The convention to use for creating strings from negative numbers.
 */
filename_sequence::filename_sequence( const frantic::tstring& dir, const frantic::tstring& filenamePattern,
                                      int negativeConvention )
    : m_pattern( frantic::tstring( ensure_trailing_pathseparator( backward_slashes( dir ) ) + filenamePattern ),
                 negativeConvention )
    , m_set() {
    ;
}

/**
 * filename_sequence constructor.  This constructor builds up the pattern from the individual components.  The negative
 * convention is set to the default FS_REPLACE_LEADING_ZERO convention.
 *
 * @param	dir					The directory for the pattern.  This can be an empty string.
 * @param	prefix				The filename prefix.
 * @param	numDigits			The minimum number of digits to use for the sequence numbers.
 * @param	extension			The filename extension.
 */
filename_sequence::filename_sequence( const frantic::tstring& dir, const frantic::tstring& prefix, int numDigits,
                                      const frantic::tstring& extension )
    : m_pattern( dir, prefix, numDigits, extension, FS_REPLACE_LEADING_ZERO )
    , m_set() {
    ;
}

/**
 * filename_sequence constructor.  This constructor builds up the pattern from the individual components.
 *
 * @param	dir					The directory for the pattern.  This can be an empty string.
 * @param	prefix				The filename prefix.
 * @param	numDigits			The minimum number of digits to use for the sequence numbers.
 * @param	extension			The filename extension.
 * @param	negativeConvention	The negative convention to be used for the pattern.
 */
filename_sequence::filename_sequence( const frantic::tstring& dir, const frantic::tstring& prefix, int numDigits,
                                      const frantic::tstring& extension, int negativeConvention )
    : m_pattern( dir, prefix, numDigits, extension, negativeConvention )
    , m_set() {
    ;
}

// filename_pattern shortcuts

/**
 * Generates a sequenced filename based on the supplied wholeframeNumber.
 *
 * @param	wholeframeNumber	The framenumber for which the sequenced filename will be generated.
 * @return	A sequenced filename based on the supplied wholeframeNumber.
 */
frantic::tstring filename_sequence::operator[]( int wholeframeNumber ) const { return m_pattern[wholeframeNumber]; }

/**
 * Generates a sequenced filename based on the supplied subframeNumber.
 *
 * @param	subframeNumber	The framenumber for which the sequenced filename will be generated.
 * @return	A sequenced filename based on the supplied wholeframeNumber.
 */
frantic::tstring filename_sequence::operator[]( double subframeNumber ) const { return m_pattern[subframeNumber]; }

// frame_set shortcuts

/**
 * Retrieves all the frame numbers in the frame_set.
 *
 * @param	allframeNumbers	This out parameter will hold the frame numbers.
 */
void filename_sequence::allframe_numbers( std::vector<double>& allframeNumbers ) const {
    m_set.allframe_numbers( allframeNumbers );
}

/**
 * Retrieves all subframe (non-integer) frame numbers in the frame_set.
 *
 * @param	subframeNumbers	This out parameter will hold the subframe (non-integer) frame numbers.
 */
void filename_sequence::subframe_numbers( std::vector<double>& subframeNumbers ) const {
    m_set.subframe_numbers( subframeNumbers );
}

/**
 * Retrieves all integer frame numbers in the frame_set.
 *
 * @param	wholeFrameNumbers	This out parameter will hold the integer frame numbers.
 */
void filename_sequence::wholeframe_numbers( std::vector<int>& wholeFrameNumbers ) const {
    m_set.wholeframe_numbers( wholeFrameNumbers );
}

/**
 * Retrieves all integer frame numbers in the frame_set.
 *
 * @param	wholeFrameNumbers	This out parameter will hold the integer frame numbers.
 */
void filename_sequence::wholeframe_numbers( std::vector<double>& wholeFrameNumbers ) const {
    m_set.wholeframe_numbers( wholeFrameNumbers );
}

/**
 * Retrieves the nearest bracketing whole framenumbers to frame in the member frame set.
 * Also returns the fractional alpha value indicating how far from the start of the interval
 * that frame lies.
 * Returns false if unable to find bracketing frames and doesn't touch any of the return data.
 * On success, if alpha is negative, that indicates that no interpolation was required, and that
 * the interval is a single value.
 *
 * @param		frame			the frame number to find the bracket for
 * @param[out]	outInterval		the return interval of closest frames surrounding the given frame
 * @param[out]	outAlpha		the alpha value, between 0 and 1, indicating the location of the given frame in the
 * interval
 */
bool filename_sequence::get_nearest_wholeframe_interval( const double frame, std::pair<double, double>& outInterval,
                                                         float& outAlpha ) const {
    return m_set.get_nearest_wholeframe_interval( frame, outInterval, outAlpha );
}

/**
 * Retrieves the nearest bracketing framenumbers to frame in the member frame set, including subframes.
 * Also returns the fractional alpha value indicating how far from the start of the interval
 * that frame lies.
 * Returns false if unable to find bracketing frames and doesn't touch any of the return data.
 * On success, if alpha is negative, that indicates that no interpolation was required, and that
 * the interval is a single value.
 *
 * @param	frame			the frame number to find the bracket for
 * @param	outInterval		the return interval of closest frames surrounding the given frame
 * @param	outAlpha		the alpha value, between 0 and 1, indicating the location of the given frame in the
 * interval
 */
bool filename_sequence::get_nearest_subframe_interval( const double frame, std::pair<double, double>& outInterval,
                                                       float& outAlpha ) const {
    return m_set.get_nearest_subframe_interval( frame, outInterval, outAlpha );
}

/**
 * Retrieves the nearest whole framenumber to frame in the frame set.
 *
 * Return false if unable to find the frame.
 *
 * @param	givenFrame		the frame number to find
 * @param	outFrame		the return frame found at the rounded given frame
 */
bool filename_sequence::get_nearest_wholeframe( const double givenFrame, double& outFrame ) const {
    return m_set.get_nearest_wholeframe( givenFrame, outFrame );
}

/**
 * Retrieves the nearest subframe in the set to the given frame.
 *
 * Return false if unable to find bracketing frames.
 *
 * @param	givenFrame		the frame number to find
 * @param	outFrame		the return number of closest frames surrounding the given frame
 */
bool filename_sequence::get_nearest_subframe( const double givenFrame, double& outFrame ) const {
    return m_set.get_nearest_subframe( givenFrame, outFrame );
}

// utility functions
/**
 * Retrieves all sequenced filenames for all frames in the frame_set.
 *
 * @param	filenames	This out parameter will hold the sequenced filenames.
 */
void filename_sequence::allframe_filenames( std::vector<frantic::tstring>& filenames ) const {
    filenames.clear();
    for( frame_set::const_iterator itr = m_set.begin(); itr != m_set.end(); ++itr )
        filenames.push_back( m_pattern[*itr] );
}

/**
 * Retrieves the sequenced filenames for all subframe (non-integer) frame numbers in the frame_set.
 *
 * @param	filenames	This out parameter will hold the subframe (non-integer) sequenced filenames.
 */
void filename_sequence::subframe_filenames( std::vector<frantic::tstring>& filenames ) const {
    filenames.clear();
    vector<double> subframeNumbers;
    m_set.subframe_numbers( subframeNumbers );
    for( vector<double>::const_iterator itr = subframeNumbers.begin(); itr != subframeNumbers.end(); ++itr )
        filenames.push_back( m_pattern[*itr] );
}

/**
 * Retrieves the sequenced filenames for all integer frame numbers in the frame_set.
 *
 * @param	filenames	This out parameter will hold the integer sequenced filenames.
 */
void filename_sequence::wholeframe_filenames( std::vector<frantic::tstring>& filenames ) const {
    filenames.clear();
    vector<int> wholeframeNumbers;
    m_set.wholeframe_numbers( wholeframeNumbers );
    for( vector<int>::const_iterator itr = wholeframeNumbers.begin(); itr != wholeframeNumbers.end(); ++itr )
        filenames.push_back( m_pattern[*itr] );
}

// file system related

/**
 * Inidcates whether the filename_pattern's directory exists.
 *
 * @return	Returns true if the directory exists.
 */
bool filename_sequence::directory_exists() const {
    if( m_pattern.get_directory( false ).length() == 0 )
        return false;

    return files::directory_exists( m_pattern.get_directory( false ) );
}

/**
 * Creates the filename_pattern's directory if it does not alreay exist.
 */
void filename_sequence::create_directory() const {
    if( m_pattern.get_directory( false ).length() == 0 )
        throw std::runtime_error(
            "filename_sequence::create_directory():  Cannot create directory from a filename_pattern "
            "with empty string for the directory." );

    if( !directory_exists() )
        files::make_directory( m_pattern.get_directory( false ) );
}

/**
 * This function searches the filename_pattern's directory for sequenced filenames that match the filename_pattern.
 * If a file's name matches the filename pattern, its frame number is added to the frame_set.
 */
void filename_sequence::sync_frame_set() {
    if( m_pattern.get_directory( false ).length() == 0 )
        throw std::runtime_error( "filename_sequence::sync_frame_set():  Cannot sync to non-existent directory." );

    // trim off the trailing backslash and check for existence.
    if( !files::directory_exists( m_pattern.get_directory( false ) ) )
        throw std::runtime_error( "filename_sequence::sync_frame_set():  Cannot sync to non-existent directory: " +
                                  frantic::strings::to_string( m_pattern.get_directory( false ) ) );

    // get the files in the directory
    vector<frantic::tstring> fileListing;
    files::get_filenames_in_directory( m_pattern.get_directory( false ), fileListing );

    m_set.clear();
    bool isPatternMatch, isValidFrameNumber;
    double frameNumber;
    for( vector<frantic::tstring>::size_type i = 0; i < fileListing.size(); ++i ) {
        isPatternMatch = m_pattern.matches_pattern(
            frantic::tstring( m_pattern.get_directory( true ) + fileListing[i] ), isValidFrameNumber, frameNumber );
        if( isPatternMatch && isValidFrameNumber )
            m_set.add_frame( frameNumber );
    }
}

/**
 * This function deletes the sequenced file from disk indicated by the frame number
 * and also removes the frame number from the frame_set.
 */
void filename_sequence::delete_frame( double framenumber ) {
    files::delete_file( m_pattern[framenumber] );
    m_set.remove_frame( framenumber );
}

// member objects
/**
 * @return	Returns the address of filename_pattern member object of this filename_sequence.
 */
filename_pattern& filename_sequence::get_filename_pattern() { return m_pattern; }
const filename_pattern& filename_sequence::get_filename_pattern() const { return m_pattern; }

/**
 * @return	Returns the address of the frame_set member object of this filename_sequence.
 */
frame_set& filename_sequence::get_frame_set() { return m_set; }
const frame_set& filename_sequence::get_frame_set() const { return m_set; }

} // namespace files
} // namespace frantic
