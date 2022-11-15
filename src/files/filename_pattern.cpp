// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

//#include <ShlObj.h>
#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <string>

#include <frantic/files/filename_sequence.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/misc/string_functions.hpp>

//#pragma warning(push)
//#pragma warning( disable : 4103 4244 )
//#include <boost/filesystem/operations.hpp>
//
//#if defined(_MSC_VER) && _MSC_VER > 1200 // Only use regex here for Visual Studio 7 or higher
//#define BOOST_REGEX_STATIC_LINK
//#include <boost/regex.hpp>
//#endif

//#pragma warning(pop)

using namespace std;
using namespace boost;
using namespace frantic;

namespace frantic {
namespace files {

// filename_pattern implementation
// constructors

/**
 * filename_pattern constructor.  This is the default constructor.
 */
filename_pattern::filename_pattern() {
    this->m_negConvention = FS_REPLACE_LEADING_ZERO;

    this->m_directory = FS_DEFAULT_DIRECTORY;
    this->m_prefix = FS_DEFAULT_PREFIX;
    this->m_numDigits = FS_DEFAULT_NUMDIGITS;
    this->m_extension = FS_DEFAULT_EXTENSION;
}

/**
 * filename_pattern constructor.  For this construtor, the supplied path should contain
 * an example file (i.e., the path cannot be just a folder).  The file portion can be a pattern,
 * containing #### for the sequence numbers.
 *
 * @param	path	An example path on which the pattern will be based.
 */
filename_pattern::filename_pattern( const frantic::tstring& path ) { set( path ); }

/**
 * filename_pattern constructor.  For this construtor, the supplied path should contain
 * an example file (i.e., the path cannot be just a folder).  The file portion can be a pattern,
 * containing #### for the sequence numbers.
 *
 * @param	path	An example path on which the pattern will be based.
 * @param	negativeConvention	The negative convention to be used for the pattern.
 */
filename_pattern::filename_pattern( const frantic::tstring& path, int negativeConvention ) {
    set_negative_convention( negativeConvention );

    frantic::tstring dir, prefix, sequence, extension;
    int numDigits;
    decompose_path( path, negativeConvention, dir, prefix, sequence, numDigits, extension );

    numDigits = ( numDigits < FS_FRAMENUMLENGTH_MIN ) ? FS_DEFAULT_NUMDIGITS : numDigits;

    set_directory( dir );
    set_prefix( prefix );
    set_num_digits( numDigits );
    set_extension( extension );
}

/**
 * filename_pattern constructor.  This constructor builds up the pattern from the individual components.
 *
 * @param	dir			The directory for the pattern.  This can be an empty string.
 * @param	prefix		The filename prefix.
 * @param	numDigits	The minimum number of digits to use for the sequence numbers.
 * @param	extension	The filename extension.
 * @param	negativeConvention	The negative convention to be used for the pattern.
 */
filename_pattern::filename_pattern( const frantic::tstring& dir, const frantic::tstring& prefix, int numDigits,
                                    const frantic::tstring extension, int negativeConvention ) {
    set_directory( dir );
    set_prefix( prefix );
    set_num_digits( numDigits );
    set_extension( extension );
    set_negative_convention( negativeConvention );
}

// ==== PUBLIC FUNCTIONS ====

/**
 * Retrieves the pattern string for this filename_pattern object.
 *
 * @return The pattern string, placing number signs #### where the sequence digits will be.
 */
frantic::tstring filename_pattern::get_pattern() const {
    return ( m_directory + m_prefix + frantic::tstring( m_numDigits, '#' ) + m_extension );
}

/**
 * Indicates whether the supplied filename matches the pattern of this filename_pattern object.
 * The supplied filename can be a typical sequenced filename, or a pattern containing #### signs
 * for the digits.
 *
 * @param	path	The fully specified filename to be tested.
 * @return	Returns true if the filename matches the pattern.
 */
bool filename_pattern::matches_pattern( const frantic::tstring& path ) const {

    frantic::tstring sequencePart;
    bool matchesPattern = test_patternmatch( path, sequencePart );

    return matchesPattern;
}

/**
 * Indicates whether the supplied filename matches the pattern of this filename_pattern object.
 * The supplied filename can be a typical sequenced filename, or a pattern containing #### signs
 * for the digits.  The additional parameters are out parameters which will indicate whether the
 * filename contains and valid frame number, and the actual frame number if so.
 *
 * @param	path	The fully specified filename to be tested.
 * @param	isValidFrameNumber	An out parameter that will be true if the filename constians a valid frame
 * number.
 * @param	frameNumber			An out parameter that will contian the frame number, if isValidFrameNumber is
 * true.
 * @return	Returns true if the filename matches the pattern.
 */
bool filename_pattern::matches_pattern( const frantic::tstring& path, bool& isValidFrameNumber,
                                        double& frameNumber ) const {
    // test for pattern match and get the sequencePart.
    frantic::tstring sequencePart;
    if( !test_patternmatch( path, sequencePart ) ) {
        isValidFrameNumber = false;
        frameNumber = 0;
        return false;
    }

    // check for '#' characters.
    if( sequencePart.find( '#' ) != string::npos ) {
        // pattern is valid, but framenumber is not because it has # in it.
        isValidFrameNumber = false;
        frameNumber = 0;
        return true;
    }

    // replace the comma that's used for decimal separation with a dot.
    string::size_type commaPos = sequencePart.find( ',' );
    if( commaPos != string::npos )
        sequencePart[commaPos] = '.';

    // do the conversion
    std::basic_stringstream<frantic::tchar> ss( sequencePart );
    ss.imbue( std::locale::classic() ); // "C" locale with '.' decimal separator
    if( ss >> frameNumber ) {
        isValidFrameNumber = true;
    } else {
        // TODO: what?  We used to call strtod(), which would return 0 on failure
        frameNumber = 0;
        isValidFrameNumber = false;
    }
    return true;
}

/**
 * Sets the filename_pattern object based on the example input file name path.
 *
 * @param	path	An example path on which the pattern will be based.
 */
void filename_pattern::set( const frantic::tstring& path ) {
    // Assign a default negative convention
    m_negConvention = FS_REPLACE_LEADING_ZERO;

    frantic::tstring dir, prefix, sequence, extension;
    int numDigits;
    decompose_path( path, m_negConvention, dir, prefix, sequence, numDigits, extension );

    numDigits = ( numDigits < FS_FRAMENUMLENGTH_MIN ) ? FS_DEFAULT_NUMDIGITS : numDigits;

    set_directory( dir );
    set_prefix( prefix );
    set_num_digits( numDigits );
    set_extension( extension );
}

/**
 * Retrieves the directory for this filename_pattern object.
 *
 * @param	includeTrailingSeparator	Set to true to include a trailing separator.
 * @return	The directory for this filename_pattern object.
 */
frantic::tstring filename_pattern::get_directory( bool includeTrailingSeparator ) const {
    if( includeTrailingSeparator || ( m_directory.length() <= 1 ) ) // preserve if the directory is just "\".
        return m_directory;
    else
        return m_directory.substr( 0, m_directory.length() - 1 ); // trim trailing separator
}

/**
 * Sets the directory for this filename_pattern object.  The directory can be an empty string.
 *
 * @param directory The directory to be set.
 */
void filename_pattern::set_directory( const frantic::tstring& directory ) {
    if( directory.length() == 0 ) {
        m_directory = FS_DEFAULT_DIRECTORY;
        return;
    }

    frantic::tstring checkDir =
        ensure_trailing_pathseparator( strict_path_name( backward_slashes( directory ) ) ); // the order is important!

    // FILESYSTEM_DEPENDENCY
    frantic::tstring invalidChars = get_invalid_identifiers( checkDir );
    if( invalidChars.find( '\\' ) != string::npos )
        invalidChars.replace( invalidChars.find( '\\' ), 1, _T("") ); // remove '\\' from returned list
    if( invalidChars.find( '/' ) != string::npos )
        invalidChars.replace( invalidChars.find( '/' ), 1, _T("") ); // remove '/' from returned list
    if( invalidChars.find( ':' ) != string::npos )
        invalidChars.replace( invalidChars.find( ':' ), 1, _T("") ); // remove ':' from returned list
    if( invalidChars.find( ' ' ) != string::npos )
        invalidChars.replace( invalidChars.find( ' ' ), 1, _T("") ); // remove ' ' from returned list

    if( invalidChars.length() > 0 ) {
        string message( "filename_pattern::set_directory():  directory contains invalid characters " );
        message += string( "[" ) + frantic::strings::to_string( invalidChars ) + string( "]." );
        throw std::runtime_error( message.c_str() );
    }

    m_directory = checkDir;
}

/**
 * Retrieves the prefix for this filename_pattern object.
 *
 * @return	The prefix for this filename_pattern object.
 */
frantic::tstring filename_pattern::get_prefix() const { return m_prefix; }

/**
 * Sets the prefix for this filename_pattern object.
 */
void filename_pattern::set_prefix( const frantic::tstring& prefix ) {
    if( prefix.length() == 0 ) {
        m_prefix = FS_DEFAULT_PREFIX;
        return;
    }

    frantic::tstring invalidChars = get_invalid_identifiers( prefix );
    if( invalidChars.length() > 0 ) {
        string message( "filename_pattern::set_prefix(): " + frantic::strings::to_string( prefix ) +
                        " contains invalid characters " );
        message += string( "[" ) + frantic::strings::to_string( invalidChars ) + string( "]." );
        throw std::runtime_error( message.c_str() );
    }

    static frantic::tstring checkNumbers( _T("0123456789#") );
    if( prefix.length() > 0 && checkNumbers.rfind( prefix[prefix.length() - 1] ) != string::npos )
        throw std::runtime_error( std::string( "filename_pattern::set_prefix(): " ) +
                                  frantic::strings::to_string( prefix ) +
                                  std::string( " cannot contain trailing numbers or '#' character." ) );

    m_prefix = prefix;
}

/**
 * Retrieves the negative convention for this filename_pattern object.
 *
 * @return	The negative convention for this filename_pattern object.
 */
int filename_pattern::get_negative_convention() const { return m_negConvention; }

/**
 * Sets the negative convention for this filename_pattern object.
 *
 * @param	negativeConvention	The negative convention to be set.
 */
void filename_pattern::set_negative_convention( int negativeConvention ) {
    // Check the range on m_negConvention
    if( negativeConvention < FS_NEG_CONV_MIN || negativeConvention > FS_NEG_CONV_MAX )
        throw std::runtime_error(
            "filename_pattern::set_negative_convention():  The supplied negativeConvention is out of range." );

    m_negConvention = negativeConvention;
}

/**
 * Retrieves the number of digits being used for the sequence portion of this filename_pattern object.
 *
 * @return	The number of digits being used for the sequence portion of this filename_pattern object.
 */
int filename_pattern::get_num_digits() const { return m_numDigits; }

/**
 * Sets the number of digits to be used for the sequence portion of this filename_pattern object.
 *
 * @param	numDigits	the number of digits to be used for the sequence portion of this filename_pattern
 * object.
 */
void filename_pattern::set_num_digits( int numDigits ) {
    if( numDigits < FS_FRAMENUMLENGTH_MIN )
        throw std::runtime_error(
            "filename_pattern::set_num_digits():  numDigits cannot be less than FS_FRAMENUMLENGTH_MIN." );

    m_numDigits = numDigits;
}

/**
 * Retrieves the filename extension for this filename_pattern object.  The dot is included as part of the extension.
 *
 * @return	The filename extension for this filename_pattern object.
 */
frantic::tstring filename_pattern::get_extension() const { return m_extension; }

/**
 * Sets the filename extension for this filename_pattern object.  The dot will be prefixed if it is missing from the
 * supplied string.
 *
 * @param	extension	the filename extension for this fileanme_pattern object.
 */
void filename_pattern::set_extension( const frantic::tstring& extension ) {
    if( extension.length() == 0 ) {
        m_extension = FS_DEFAULT_EXTENSION;
        return;
    }

    frantic::tstring invalidChars = get_invalid_identifiers( extension );
    if( invalidChars.length() > 0 ) {
        string message( "filename_pattern::set_extension():  extension contains invalid characters " );
        message += string( "[" ) + frantic::strings::to_string( invalidChars ) + string( "]." );
        throw std::runtime_error( message.c_str() );
    }

    string::size_type dotPos = extension.rfind( '.' );
    if( dotPos == string::npos ) {
        m_extension = _T(".") + extension; // prepend the missing dot.
        return;
    } else {
        if( dotPos != 0 ) // is there a dot in the middle somewhere?
            throw std::runtime_error(
                "filename_pattern::set_extension(): extension cannot contain dot in character position other than 0." );
    }

    m_extension = extension;
}

/**
 * Generates a sequenced filename based on the supplied wholeframeNumber.
 *
 * @param	wholeframeNumber	The framenumber for which the sequenced filename will be generated.
 * @return	A sequenced filename based on the supplied wholeframeNumber.
 */
frantic::tstring filename_pattern::operator[]( int wholeframeNumber ) const {
    bool isNeg = ( wholeframeNumber < 0 );

    if( isNeg && m_negConvention == FS_NON_NEGATIVE )
        throw std::runtime_error(
            "filename_pattern::operator[] (int& frameNumber) cannot form negative filename when the "
            "negative convention is FS_NON_NEGATIVE." );

    frantic::tstring seqString( build_whole_frame_seq_string( wholeframeNumber, false ) );

    return ( m_directory + m_prefix + seqString + m_extension );
}

/**
 * Generates a sequenced filename based on the supplied subframeNumber.
 *
 * @param	subframeNumber	The framenumber for which the sequenced filename will be generated.
 * @return	A sequenced filename based on the supplied subframeNumber.
 */
frantic::tstring filename_pattern::operator[]( double subframeNumber ) const {
    bool isNeg = ( subframeNumber < 0 );

    if( isNeg && m_negConvention == FS_NON_NEGATIVE )
        throw std::runtime_error(
            "filename_pattern::operator[] (double& frameNumber) cannot form negative filename when "
            "the negative convention is FS_NON_NEGATIVE." );

    frantic::tstring seqString( build_whole_frame_seq_string( int( subframeNumber ), isNeg ) );
    // if (seqString.length() + 1 > FS_FRAMENUMLENGTH_MAX)
    //	throw std::runtime_error("filename_pattern::operator[] (double& frameNumber) integer portion of supplied double
    // exceeds maximum allowed digits.");

    // check for whole frame
    if( subframeNumber - int( subframeNumber ) == 0 )
        return ( m_directory + m_prefix + seqString + m_extension );

    // build a string of the decimal part
    const double subframeDecimal = fabs( subframeNumber ) - floor( fabs( subframeNumber ) );

    std::basic_stringstream<frantic::tchar> ss;
    ss.imbue( std::locale::classic() ); // Make sure we're using the "C" locale, with a '.' decimal separator.
    ss.precision( FS_DECIMALLENGTH_MAX );
    ss.setf( ios::fixed ); // fixed point representation
    ss << subframeDecimal;
    const frantic::tstring stringRep = ss.str();

    frantic::tstring decString; // subframe number after the decimal point
    if( stringRep == _T("0") ) {
        // Check for rounding to a whole frame number.
        // This empty string is detected and an error is thrown later on.
        decString = _T("");
    } else if( starts_with( stringRep, _T("0.") ) && stringRep.length() > 2 ) {
        // This should be the usual decimal case, for example 0.5
        decString = stringRep.substr( 2, min( frantic::tstring::size_type( FS_DECIMALLENGTH_MAX ),
                                              stringRep.size() ) ); // truncate to max decimal length
    } else if( starts_with( stringRep, _T("1") ) ) {
        // The decimal part was rounded up to 1, bump it back down to 0.999...
        decString = frantic::tstring( frantic::tstring::size_type( FS_DECIMALLENGTH_MAX ), '9' );
    } else {
        throw std::runtime_error( "filename_pattern::operator[] Error: unexpected result \'" +
                                  frantic::strings::to_string( stringRep ) +
                                  "\' while creating decimal string for frame number " +
                                  boost::lexical_cast<std::string>( subframeNumber ) + "." );
    }

    // TODO:  How should we handle rounding when truncation is necessary?  For example, rounding up could cause the
    // number to roll to the next
    //  whole frame number, and that would be bad.
    //  - for now such numbers are rounded back down to 0.999...

    // trim off trailing zeroes
    frantic::tstring::size_type nonZeroPos = decString.find_last_not_of( _T("0") );

    if( nonZeroPos == string::npos )
        decString = _T("");
    else
        decString = decString.substr( 0, nonZeroPos + 1 );

    // check for rounding to whole frame number.
    if( decString.length() == 0 )
        throw std::runtime_error(
            "filename_pattern::operator[] (double& frameNumber) decimal portion of supplied double '" +
            boost::lexical_cast<std::string>( subframeNumber ) + " is too small to represent." );

    seqString = seqString + _T(",") + decString;

    return ( m_directory + m_prefix + seqString + m_extension );
}

frantic::tstring filename_pattern::add_before_sequence_number( const frantic::tstring& path,
                                                               const frantic::tstring& addition,
                                                               int negativeConvention ) {
    int numSequenceDigits;
    frantic::tstring dirPart, filePart, sequencePart, extPart;

    decompose_path( path, negativeConvention, dirPart, filePart, sequencePart, numSequenceDigits, extPart );

    return dirPart + filePart + addition + sequencePart + _T(".") + extPart;
}

// == Private Methods ==

/**
 * This private function decomposes the supplied path into the directory and filename components. The directory
 * portion retains the trailing separator.
 *
 * @param	path	The path to be decomposed.
 * @param	outDirectory	An output parameter to hold the directory portion of the path.
 * @param	outFilename		An output parameter to hold the filename portion of the path.
 */
void filename_pattern::separate_directory_and_filename( const frantic::tstring& path, frantic::tstring& outDirectory,
                                                        frantic::tstring& outFilename ) {
    frantic::tstring strictPath = strict_path_name( path );
#if defined( _WIN32 )
    frantic::tstring::size_type posSlash = strictPath.rfind( '\\' );
#else
    frantic::tstring::size_type posSlash = strictPath.rfind( '/' );
#endif

    // determine where to split between the pathPart and the filePart.
    if( posSlash == string::npos ) { // all file and no path
        outDirectory = _T("");
        outFilename = strictPath;
    } else if( posSlash == strictPath.length() - 1 ) {              // all path and no file
        outDirectory = strictPath.substr( 0, strictPath.length() ); // keep the trailing slash
        outFilename = _T("");
    } else {
        // slash in there somewhere
        outDirectory = strictPath.substr( 0, posSlash + 1 ); // up to and inlcuding the last backslash
        outFilename = strictPath.substr( posSlash + 1,
                                         strictPath.length() - posSlash ); // from after the last backslash to the end;
    }
}

/**
 * This private function fully decomposes the supplied path into the directory, prefix, sequence & numDigits, and
 * extension components.  Note that the supplied negativeConvention affect how the path is decomposed.
 *
 * @param	path			The path to be decomposed.
 * @param	negativeConvention	The negative convention on which to base the decomposition.
 * @param	outDirectory	An out parameter that will hold the directory portion of the pattern.
 * @param	outPrefix		An out paramter that will hold the prefix portion of the pattern.
 * @param	outSequence		An out parameter that will hold the sequence portion of the pattern.
 * @param	outNumDigits	An out parameter that will hold the number of sequence digits in the pattern.
 * @param	outExtension	An out parameter that will hold the filename extension portion of the pattern.
 */
void filename_pattern::decompose_path( const frantic::tstring& path, int negativeConvention,
                                       frantic::tstring& outDirectory, frantic::tstring& outPrefix,
                                       frantic::tstring& outSequence, int& outNumDigits,
                                       frantic::tstring& outExtension ) {
    // Clean the path
    frantic::tstring strictPath = strict_path_name( path );

    // Check for empty path
    if( strictPath.length() == 0 ) {
        outDirectory = _T("");
        outPrefix = _T("");
        outSequence = _T("");
        outNumDigits = 0;
        outExtension = _T("");
        return;
    }

    frantic::tstring filename, prefix, sequence, ext;
    separate_directory_and_filename( strictPath, outDirectory, filename );

    enum DecomposeState { UNKNOWN = 0, EXTENSION, WHOLE_OR_SUB_FRAME, WHOLE_FRAME, FOUND_LETTERS, DONE };

    int prefixLength = 0;
    int letterIndex = (int)filename.length() - 1;
    int extensionStart = (int)filename.length();
    int extensionLength = 0;
    DecomposeState currState = EXTENSION;
    int subFrameDigits = 0;
    int wholeFrameDigits = 0;
    bool isNegative = false;

    // begin crawl
    currState = EXTENSION;

    while( letterIndex >= 0 && currState != DONE ) {
        if( currState == EXTENSION ) {
            if( filename[letterIndex] == '.' ) {
                extensionStart = letterIndex + 1;
                extensionLength = (int)filename.length() - extensionStart;
                currState = WHOLE_OR_SUB_FRAME;
            }
        } else if( ( filename[letterIndex] >= '0' && filename[letterIndex] <= '9' ) ||
                   ( filename[letterIndex] == '#' ) ) {
            switch( currState ) {
            case WHOLE_OR_SUB_FRAME:
            case WHOLE_FRAME:
                wholeFrameDigits++;
                break;
            default:
                currState = DONE;
            }
        } else {
            switch( filename[letterIndex] ) {
            case ',':
                switch( currState ) {
                case WHOLE_OR_SUB_FRAME:
                    subFrameDigits = wholeFrameDigits;
                    wholeFrameDigits = 0;
                    currState = WHOLE_FRAME;
                    break;

                default:
                    currState = DONE;
                }
                break;

            case '-':
                isNegative = negativeConvention != FS_NON_NEGATIVE;
                currState = DONE;
                break;

            default:
                currState = DONE;
            }
        }

        --letterIndex;
    }

    // end crawl

    prefixLength = letterIndex + 1;

    // If the current state is done then the prefix length has been decremented
    // beyond the end of the prefix
    if( currState == DONE ) {
        if( isNegative ) {
            if( negativeConvention == FS_REPLACE_LEADING_ZERO )
                wholeFrameDigits++;
        } else
            prefixLength++;
    } else {
        currState = DONE;
    }

    // check for a weird case like myfile,0000.txt
    if( ( wholeFrameDigits == 0 ) && ( subFrameDigits > 0 ) ) {
        prefixLength++;
        wholeFrameDigits = subFrameDigits;
        subFrameDigits = 0;
    }

    // make final assignments
    if( prefixLength <= 0 )
        prefix = _T("");
    else
        prefix = filename.substr( 0, prefixLength );

    if( wholeFrameDigits <= 0 )
        sequence = _T("");
    else if( subFrameDigits > 0 )
        sequence = filename.substr( prefixLength, wholeFrameDigits + subFrameDigits + 1 );
    else
        sequence = filename.substr( prefixLength, wholeFrameDigits );

    if( extensionLength <= 0 )
        ext = _T("");
    else
        ext = filename.substr( extensionStart, filename.length() - extensionStart );

    // std::stringstream strstm;
    // strstm << prefix << " + " << sequence << " + " << ext;
    // throw std::runtime_error( strstm.str() );

    // Fill the out params.
    outPrefix = prefix;
    outSequence = sequence;
    outNumDigits = wholeFrameDigits;
    outExtension = ext;
}

/**
 * This private function builds the wholeframe portion of the sequence string for the supplied wholeframeNumber.
 * The makeNegative parameter instructs the fucntion to make a zero frame negative.  This is used, for example,
 * by operator[] (double) when it is given a frame like -0.25, where the whole frame part is zero, but must be negative.
 *
 * @param	wholeframeNumber	The integer frame number for which the wholeframe portion of the sequence string will
 * be built.
 * @param	makeNegative		Indicates that the returned string should be negative.
 * @return	The wholeframe portion of the sequence string for the supplied wholeframeNumber.
 */
frantic::tstring filename_pattern::build_whole_frame_seq_string( int wholeframeNumber, bool makeNegative ) const {
    bool isNeg = ( wholeframeNumber < 0 );

    frantic::tstring fileNumStr( frantic::strings::to_tstring( 
        std::to_string( std::abs( wholeframeNumber ) )
    ) );

    // pad to full length
    if( int( fileNumStr.length() ) < m_numDigits )
        fileNumStr = frantic::tstring( m_numDigits - fileNumStr.length(), '0' ) + fileNumStr;

    // handle negatives
    if( isNeg || makeNegative ) {
        switch( m_negConvention ) {
        case FS_REPLACE_LEADING_ZERO:
            fileNumStr = ( fileNumStr[0] == '0' ) ? ( _T("-") + fileNumStr.substr( 1, fileNumStr.length() - 1 ) )
                                                  : ( _T("-") + fileNumStr );
            break;

        case FS_PREFIX:
            fileNumStr = _T("-") + fileNumStr;
            break;

        default:
            throw std::runtime_error( "filename_pattern::operator[] unhandled negative convention." );
        }
    }

    return fileNumStr;
}

/**
 * This private function returns a string of invalid identifiers contained in the supplied string.
 *
 * @param	s	The string to be checked for invalid indentifiers.
 * @return	A string of invalid identifiers contained in the supplied string.
 */
frantic::tstring filename_pattern::get_invalid_identifiers( const frantic::tstring& s ) const {
    if( s.length() == 0 )
        return frantic::tstring();

    frantic::tstring invalids;

    for( frantic::tstring::size_type i = 0; i < FS_INVALID_IDENTIFIERS.length(); ++i )
        if( s.find( FS_INVALID_IDENTIFIERS[i] ) != string::npos )
            invalids += FS_INVALID_IDENTIFIERS[i];

    return invalids;
}

/**
 * This private function tests the supplied path for a match to the pattern of this filename_pattern object.
 * If the path is a match, the sequence portion of the string is returned via the outSequencePart out parameter.
 *
 * @param	path			The path to be tested for a match to the pattern of this filename_pattern
 * object.
 * @param	outSequencePart	An output paramter to hold the sequence portion if the path is a match.
 * @return	Returns true if the supplied patch is a match to the pattern of this filename_pattern object.
 */
bool filename_pattern::test_patternmatch( const frantic::tstring& path, frantic::tstring& outSequencePart ) const {

    // FILESYSTEM_DEPENDENCY - some filesystems may be case sensitive.
    // Use defines to change the comparible versions accordingly.
    frantic::tstring comparable_Path = strings::to_lower( path );
    frantic::tstring comparable_mDir = strings::to_lower( m_directory );
    frantic::tstring comparable_mPrefix = strings::to_lower( m_prefix );
    frantic::tstring comparable_mExt = strings::to_lower( m_extension );

    frantic::tstring dirPart, filePart;
    separate_directory_and_filename( comparable_Path, dirPart, filePart );

    // check the directory
    if( dirPart != comparable_mDir ) {
        outSequencePart = _T("");
        return false;
    }

    // check the extension and prefix parts
    frantic::tstring leftOfExt, ext;
    frantic::tstring::size_type posDot = filePart.rfind( '.' );
    if( posDot == string::npos ) {
        ext = _T("");
        leftOfExt = filePart;
    } else {
        ext = filePart.substr( posDot, filePart.length() - posDot );
        leftOfExt = filePart.substr( 0, posDot );
    }

    // check the extension, length of the prefix, and content of the prefix
    if( ext != comparable_mExt || leftOfExt.length() <= comparable_mPrefix.length() ||
        leftOfExt.substr( 0, comparable_mPrefix.length() ) != comparable_mPrefix ) {
        outSequencePart = _T("");
        return false;
    }

    // Isolate the sequencePart.
    outSequencePart = leftOfExt.substr( comparable_mPrefix.length(), leftOfExt.length() - comparable_mPrefix.length() );

    if( outSequencePart.length() == 0 )
        return false;

    // check that the sequence part contains only valid characters. (This part might be better handled wtih regex.)
    frantic::tstring validChars( _T("-0123456789,#") );
    int dashCount = 0, commaCount = 0, commaPos = -1;
    for( string::size_type i = 0; i < outSequencePart.length(); ++i ) {
        if( outSequencePart[i] == '-' )
            ++dashCount;

        if( outSequencePart[i] == ',' ) {
            ++commaCount;
            commaPos = (int)i;
        }

        if( validChars.find( outSequencePart[i] ) == string::npos || commaCount > 1 || dashCount > 1 ) {
            return false;
        }
    }

    // If there is a negative symbol, it can only be on the far left.
    if( dashCount == 1 && outSequencePart[0] != '-' )
        return false;

    frantic::tstring wholeframePart;

    // subframe handling
    if( commaCount == 1 ) {
        // trailing zero or comma not permitted due to round-trip invariant rule.
        if( outSequencePart[outSequencePart.length() - 1] == '0' ||
            outSequencePart[outSequencePart.length() - 1] == ',' )
            return false;

        wholeframePart = outSequencePart.substr( 0, commaPos );
    } else
        wholeframePart = outSequencePart;

    // check for sufficient total digits for round-trip invariant rule.  Applies to all negative convention cases.
    int wholeframePartLen = (int)wholeframePart.length();
    if( wholeframePartLen < m_numDigits )
        return false;

    // negative convention handling
    if( dashCount == 1 ) {
        switch( m_negConvention ) {
        case FS_REPLACE_LEADING_ZERO:
            if( wholeframePartLen == m_numDigits || wholeframePartLen == m_numDigits + 1 )
                return true;

            return ( wholeframePart[1] != '0' ); // wholeframe has extra digits, so first digit after negative symbol
                                                 // cannot be zero due to round-trip invariant rule.

        case FS_PREFIX:
            if( wholeframePartLen < m_numDigits + 1 )
                return false;

            if( wholeframePartLen == m_numDigits + 1 )
                return true;

            return ( wholeframePart[1] != '0' ); // wholeframe has extra digits, so first digit after negative symbol
                                                 // cannot be zero due to round-trip invariant rule.

        case FS_NON_NEGATIVE:
            return false;

        default:
            throw std::runtime_error( "test_patternmatch():  unhandled negative convention." );
        }
    }

    // non-negative handling
    if( wholeframePartLen == m_numDigits )
        return true;

    // if the digit count is greater than m_numDigits, the leading digit cannot be zero due to round-trip invariant
    // rule.
    return ( wholeframePart[0] != '0' );
}

} // namespace files
} // namespace frantic
