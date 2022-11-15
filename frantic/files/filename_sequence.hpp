// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <set>
#include <string>
#include <vector>

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace files {

const static int FS_REPLACE_LEADING_ZERO = 0; // myFile####.ext becomes myFile-###.ext
const static int FS_PREFIX = 1;               // myFile####.ext becomes myFile-####.ext
const static int FS_NON_NEGATIVE = 2;         // specifies that there are no negative frame numbers

const static int FS_NEG_CONV_MIN = 0;       // the lowest negative convention id
const static int FS_NEG_CONV_MAX = 2;       // the highest negative convention id
const static int FS_FRAMENUMLENGTH_MIN = 1; // the minimum number of digits in the whole portion of the frame number
// const static int FS_FRAMENUMLENGTH_MAX = 8;   // the maximum number of digits in the whole portion of the frame
// number
const static int FS_DECIMALLENGTH_MAX =
    7; // the maximum number of digits in the decimal portion of the frame number (not counting the '.')

const static frantic::tstring FS_DEFAULT_DIRECTORY( _T("") );
const static frantic::tstring FS_DEFAULT_PREFIX( _T("") );
const static int FS_DEFAULT_NUMDIGITS = 4;
const static frantic::tstring FS_DEFAULT_EXTENSION( _T(".") );

const static frantic::tstring
    FS_INVALID_IDENTIFIERS( _T("/\\\t?\":><|\n") ); // Do not include '.' character in this list.

/**
 * A filename_pattern object manages the naming pattern for sequenced filenames. The filename_pattern class
 * is intended primarily to support the filename_sequence class. However, in cases where the set of frames
 * does not need to be managed, a filename_pattern can be used on its own. Once constructed, the filename_pattern
 * can be used to easily construct sequenced filenames and to identify files that match the pattern.
 *
 * @author	James Coulter
 */
class filename_pattern {
  public:
    // construction
    filename_pattern();
    filename_pattern( const frantic::tstring& path );
    filename_pattern( const frantic::tstring& path, int negativeConvention );
    filename_pattern( const frantic::tstring& dir, const frantic::tstring& prefix, int numDigits,
                      const frantic::tstring extension, int negativeConvention );

    frantic::tstring get_pattern() const;
    bool matches_pattern( const frantic::tstring& path ) const;
    bool matches_pattern( const frantic::tstring& path, bool& isValidFrameNumber, double& frameNumber ) const;

    void set( const frantic::tstring& path );

    frantic::tstring get_directory( bool includeTrailingSeparator ) const;
    void set_directory( const frantic::tstring& directory );

    frantic::tstring get_prefix() const;
    void set_prefix( const frantic::tstring& prefix );

    int get_negative_convention() const;
    void set_negative_convention( int negativeConvention );

    int get_num_digits() const;
    void set_num_digits( int numDigits );

    frantic::tstring get_extension() const;
    void set_extension( const frantic::tstring& extension );

    frantic::tstring operator[]( int wholeframeNumber ) const;
    frantic::tstring operator[]( double subframeNumber ) const;

    /**
     * Given a path that is acceptable to a filename_pattern object, add an additional string to the filename prefix
     * before the sequence number.
     * @param path The path to modify
     * @param addition The string to insert into the path
     * @param negativeConvention The convention to use for creating strings from negative numbers.
     * @result The modified path with the additional string added before the sequence number.
     */
    static frantic::tstring add_before_sequence_number( const frantic::tstring& path, const frantic::tstring& addition,
                                                        int negativeConvention = FS_REPLACE_LEADING_ZERO );

  private:
    int m_negConvention, m_numDigits;
    frantic::tstring m_directory;
    frantic::tstring m_prefix;
    frantic::tstring m_extension;

    static void separate_directory_and_filename( const frantic::tstring& path, frantic::tstring& outDirectory,
                                                 frantic::tstring& outFilename );
    static void decompose_path( const frantic::tstring& path, int negativeConvention, frantic::tstring& outDirectory,
                                frantic::tstring& outPrefix, frantic::tstring& outSequence, int& outNumDigits,
                                frantic::tstring& outExtension );

    frantic::tstring build_whole_frame_seq_string( int wholeframeNumber, bool makeNegative ) const;
    frantic::tstring get_invalid_identifiers( const frantic::tstring& s ) const;
    bool test_patternmatch( const frantic::tstring& path, frantic::tstring& outSequencePart ) const;
};

/**
 * A frame_set object manages a set of frame numbers. The frame_set class is intended primarily
 * to support the filename_sequence class. It is essentially a fancy wrapper for std::set<double>
 * aimed at frame number management. Management of frames is broken down into three categories for
 * most operations. These are:
 * <ul>
 * <li>allframes: This indicates all frames in the set.
 * <li>subframes: This indicates only the non-integer frames in the set.
 * <li>wholeframes: This indicates only the integer frames in the set.
 * </ul>
 * <p>
 *
 * @author	James Coulter
 */
class frame_set {
  public:
    typedef std::set<double>::size_type size_type;
    typedef std::set<double>::iterator iterator;
    typedef std::set<double>::const_iterator const_iterator;
    typedef std::set<double>::reverse_iterator reverse_iterator;
    typedef std::set<double>::const_reverse_iterator const_reverse_iterator;

    // construction
    frame_set();
    frame_set( const std::vector<double>& frames );
    frame_set( const std::vector<int>& wholeframes );
    frame_set( int startFrame, int endFrame );
    frame_set( const frame_set& rhs );

    // set operations
    int add_frame( int wholeframe );
    double add_frame( double frame );
    void add_frames( const std::vector<int>& wholeframes );
    void add_frames( const std::vector<double>& frames );
    int remove_frame( int wholeframe );
    double remove_frame( double frame );
    void remove_frames( const std::vector<int>& wholeframes );
    void remove_frames( const std::vector<double>& frames );
    void clear();

    // set inspection
    bool empty() const;
    bool has_subframes() const;
    bool has_wholeframes() const;
    bool frame_exists( double frame ) const;
    size_type size() const; // number of frames in the set

    std::pair<double, double> allframe_range() const;
    std::pair<double, double> subframe_range() const;
    std::pair<int, int> wholeframe_range() const;

    void missing_wholeframes( std::vector<int>& missingWholeFrames ) const;
    void missing_wholeframes( int lowFrameNumber, int highFrameNumber, std::vector<int>& missingWholeFrames ) const;
    bool contiguous() const;
    bool contiguous( int lowFrameNumber, int highFrameNumber ) const;

    // retrieval
    void allframe_numbers( std::vector<double>& allframeNumbers ) const;
    void subframe_numbers( std::vector<double>& subframeNumbers ) const;
    void wholeframe_numbers( std::vector<double>& wholeframeNumbers ) const;
    void wholeframe_numbers( std::vector<int>& wholeframeNumbers ) const;

    bool get_nearest_wholeframe_interval( const double frame, std::pair<double, double>& outInterval,
                                          float& outAlpha ) const;
    bool get_nearest_subframe_interval( const double frame, std::pair<double, double>& outInterval,
                                        float& outAlpha ) const;

    bool get_nearest_wholeframe( const double givenFrame, double& outFrame ) const;
    bool get_nearest_subframe( const double givenFrame, double& outFrame ) const;

    // iterator support
    iterator begin(); // points to first element
    const_iterator begin() const;
    iterator end(); // points to one-past-last element
    const_iterator end() const;

    reverse_iterator rbegin(); // points to first element of reverse sequence
    const_reverse_iterator rbegin() const;
    reverse_iterator rend(); // point to one-past-last element of reverse sequence
    const_reverse_iterator rend() const;

  private:
    std::set<double> m_frames;
};

/**
 * The filename_sequence class is intended to be useful for both regular file sequences and also
 * file sequences with subframe data, as is often found with simulation data. The filename_sequence
 * class has two member objects which provide most of the functionality. The filename_pattern object
 * is used to recognize and generate sequenced filenames, and the frame_set object is used to
 * manage sets of frame numbers. There are additional functions in filename_sequence which tie
 * the members together to round out the functionality.
 * <p>
 * With the exception of a set of convenience functions, a filename_sequence object does not bind to
 * the file system, nor does it assume how the sequence is being used. This places the
 * responsibility / flexibility of disk file management in the hands of the consumer.
 * <p>
 *
 * @author	James Coulter
 */
class filename_sequence {
  public:
    // constructors
    filename_sequence();
    filename_sequence( const frantic::tstring& path );
    filename_sequence( const frantic::tstring& path, int negativeConvention );
    filename_sequence( const frantic::tstring& dir, const frantic::tstring& filenamePattern );
    filename_sequence( const frantic::tstring& dir, const frantic::tstring& filenamePattern, int negativeConvention );
    filename_sequence( const frantic::tstring& dir, const frantic::tstring& prefix, int numDigits,
                       const frantic::tstring& extension );
    filename_sequence( const frantic::tstring& dir, const frantic::tstring& prefix, int numDigits,
                       const frantic::tstring& extension, int negativeConvention );

    // copy construction follows standard semantics so copy constructors are not explicitly stated.

    // shortcut functions to member objects
    // filename_pattern shortcuts
    frantic::tstring operator[]( int wholeframeNumber ) const;
    frantic::tstring operator[]( double subframeNumber ) const;

    // frame_set shortcuts
    void allframe_numbers( std::vector<double>& allframeNumbers ) const;
    void subframe_numbers( std::vector<double>& subframeNumbers ) const;
    void wholeframe_numbers( std::vector<int>& wholeFrameNumbers ) const;
    void wholeframe_numbers( std::vector<double>& wholeFrameNumbers ) const;

    // utility functions
    void allframe_filenames( std::vector<frantic::tstring>& filenames ) const;
    void subframe_filenames( std::vector<frantic::tstring>& filenames ) const;
    void wholeframe_filenames( std::vector<frantic::tstring>& filenames ) const;
    bool get_nearest_wholeframe_interval( const double frame, std::pair<double, double>& outInterval,
                                          float& outAlpha ) const;
    bool get_nearest_subframe_interval( const double frame, std::pair<double, double>& outInterval,
                                        float& outAlpha ) const;
    bool get_nearest_wholeframe( const double givenFrame, double& outFrame ) const;
    bool get_nearest_subframe( const double givenFrame, double& outFrame ) const;

    // file system related
    bool directory_exists() const;
    void create_directory() const;
    void sync_frame_set();
    void delete_frame( double framenumber );

    // member objects (both const and non-const)
    filename_pattern& get_filename_pattern();
    const filename_pattern& get_filename_pattern() const;
    frame_set& get_frame_set();
    const frame_set& get_frame_set() const;

  private:
    frantic::files::filename_pattern m_pattern;
    frantic::files::frame_set m_set;
};

} // namespace files
} // namespace frantic
