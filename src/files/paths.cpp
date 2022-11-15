// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( _WIN32 ) || defined( _WIN64 )
#include <ShlObj.h>
#endif

#include <algorithm>

#include <frantic/files/paths.hpp>
#include <frantic/misc/string_functions.hpp>
#include <frantic/strings/tstring.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>

#if( defined( _MSC_VER ) && _MSC_VER > 1200 ) ||                                                                       \
    defined( linux ) // Only use regex here for Visual Studio 7 or higher and Linux
#define BOOST_REGEX_STATIC_LINK
#include <boost/regex.hpp>
#endif

using namespace std;
using namespace boost;
using namespace frantic;
namespace fs = boost::filesystem;

namespace {
class ff_set_default_name_check {
  public:
    ff_set_default_name_check() {
        // Set the path checking for windows paths
        // Deprecated.
        // boost::filesystem::path::default_name_check( boost::filesystem::windows_name );
    }
};

ff_set_default_name_check ffsdnc;
} // namespace

namespace frantic {
namespace files {

std::string forward_slashes( std::string Path ) {
    for( string::size_type i = 0; i < Path.size(); ++i )
        if( Path[i] == '\\' )
            Path[i] = '/';

    return Path;
}

std::string backward_slashes( std::string Path ) {
    for( string::size_type i = 0; i < Path.size(); ++i )
        if( Path[i] == '/' )
            Path[i] = '\\';

    return Path;
}

std::wstring forward_slashes( std::wstring Path ) {
    for( string::size_type i = 0; i < Path.size(); ++i )
        if( Path[i] == '\\' )
            Path[i] = '/';

    return Path;
}

std::wstring backward_slashes( std::wstring Path ) {
    for( string::size_type i = 0; i < Path.size(); ++i )
        if( Path[i] == '/' )
            Path[i] = '\\';

    return Path;
}

namespace detail {
template <class CharType>
std::basic_string<CharType> ensure_trailing_pathseparator( std::basic_string<CharType> Path ) {
    // Empty strings should stay as empty strings
    if( Path.size() == 0 )
        return std::basic_string<CharType>();

#if defined( _WIN32 )
    if( Path[Path.size() - 1] != '\\' && Path[Path.size() - 1] != '\\' )
        Path += '\\';
#else
    if( Path[Path.size() - 1] != '/' && Path[Path.size() - 1] != '/' )
        Path += '/';
#endif

    return Path;
}
} // namespace detail

std::string ensure_trailing_pathseparator( std::string Path ) { return detail::ensure_trailing_pathseparator( Path ); }

std::wstring ensure_trailing_pathseparator( std::wstring Path ) {
    return detail::ensure_trailing_pathseparator( Path );
}

namespace detail {
template <class CharType>
std::basic_string<CharType> normalized_directory_name( std::basic_string<CharType> name ) {
    const std::basic_string<CharType> slash( 1, '/' );
    const std::basic_string<CharType> backslash( 1, '\\' );

#if defined( _WIN32 )
    const std::basic_string<CharType> doubleBackslash( 2, '\\' );

    // Switch to DOS style separators
    name = strings::string_replace( name, slash, backslash );

    // Protect the UNC path at the start
    bool UNC = false;
    if( name.size() > 1 && name[0] == '\\' && name[1] == '\\' ) {
        UNC = true;
        name[1] = 'X';
    }

    // Remove doubled up separators
    name = strings::string_replace( name, doubleBackslash, backslash );

    // Restore the UNC path at the start
    if( UNC )
        name[1] = '\\';
#else
    const std::basic_string<CharType> doubleSlash( 2, '/' );

    // Switch to Linux style separators
    name = strings::string_replace( name, backslash, slash );

    // Remove doubled up separators
    name = strings::string_replace( name, doubleSlash, slash );
#endif

    return name;
}
} // namespace detail

// Converts paths like "C:\/windows", "c:/amaretto/\shaders", etc to proper paths
std::string normalized_directory_name( const std::string& name ) { return detail::normalized_directory_name( name ); }

std::wstring normalized_directory_name( const std::wstring& name ) { return detail::normalized_directory_name( name ); }

std::string clean_directory_name( std::string name, std::string replacement ) {
    name = strings::string_replace( name, "/", replacement );
    name = strings::string_replace( name, "\\", replacement );
    name = strings::string_replace( name, " ", replacement );
    name = strings::string_replace( name, ".", replacement );
    name = strings::string_replace( name, "\t", replacement );
    name = strings::string_replace( name, "?", replacement );
    name = strings::string_replace( name, "\"", replacement );
    name = strings::string_replace( name, ":", replacement );
    name = strings::string_replace( name, ">", replacement );
    name = strings::string_replace( name, "<", replacement );
    name = strings::string_replace( name, "|", replacement );
    name = strings::string_replace( name, "\n", replacement );

    return name;
}

// cleans up a path name, removing ", trimming whitespaces, converting all slashes to backslashes and normalizing the
// path.
// Good for commandline arguments that have to play nice with boost::filesystem::path. Not very pretty or efficient.
std::string strict_path_name( std::string const& path ) {
    return normalized_directory_name( frantic::strings::trim( frantic::strings::string_replace( path, "\"", "" ) ) );
}

std::wstring strict_path_name( std::wstring const& path ) {
    return normalized_directory_name( frantic::strings::trim( frantic::strings::string_replace( path, L"\"", L"" ) ) );
}

std::string replace_directory( const string& Path, string ReplacementDir = "" ) {
    return ( boost::filesystem::path( files::normalized_directory_name( ReplacementDir ) ) /
             filename_from_path( Path ) )
        .string();
}

std::string replace_filename( const string& Path, const string& ReplacementName = "" ) {
    return ( boost::filesystem::path( files::normalized_directory_name( Path ) ).branch_path() / ReplacementName )
        .string();
}

std::wstring replace_filename( const std::wstring& Path, const std::wstring& ReplacementName = L"" ) {
    return ( boost::filesystem::path( files::normalized_directory_name( Path ) ).branch_path() / ReplacementName )
        .wstring();
}

boost::filesystem::path replace_extension_p( const boost::filesystem::path& Path, const string& ReplacementExt ) {
    return Path.branch_path() / ( basename_from_path_p( Path ) + ReplacementExt );
}

std::string replace_extension( const string& Path, const string& ReplacementExt ) {
    return replace_extension_p( boost::filesystem::path( files::normalized_directory_name( Path ) ), ReplacementExt )
        .string();
}

// forward declaration
std::wstring basename_from_path_p_w( const boost::filesystem::path& thePath );

boost::filesystem::path replace_extension_p( const boost::filesystem::path& Path, const std::wstring& ReplacementExt ) {
    return Path.branch_path() / ( basename_from_path_p_w( Path ) + ReplacementExt );
}

std::wstring replace_extension( const std::wstring& Path, const std::wstring& ReplacementExt ) {
    return replace_extension_p( boost::filesystem::path( files::normalized_directory_name( Path ) ), ReplacementExt )
        .wstring();
}

std::string directory_from_path( const string& Path ) {
    return boost::filesystem::path( files::normalized_directory_name( Path ) ).branch_path().string();
}

std::wstring directory_from_path( const wstring& Path ) {
    return boost::filesystem::path( files::normalized_directory_name( Path ) ).branch_path().wstring();
}

std::string filename_from_path( const string& Path ) {
    return boost::filesystem::path( files::normalized_directory_name( Path ) ).leaf().string();
}

std::wstring filename_from_path( const std::wstring& Path ) {
    return boost::filesystem::path( files::normalized_directory_name( Path ) ).leaf().wstring();
}

std::string basename_from_path_p( const boost::filesystem::path& thePath ) {
    // Grab the filename
    string filename = thePath.leaf().string();

    // Split it into the basename and the extension
    size_t extIndex = filename.rfind( '.' );
    if( extIndex != string::npos ) {
        return filename.substr( 0, extIndex );
    } else {
        return filename;
    }
}

std::wstring basename_from_path_p_w( const boost::filesystem::path& thePath ) {
    // Grab the filename
    std::wstring filename = thePath.leaf().wstring();

    // Split it into the basename and the extension
    size_t extIndex = filename.rfind( '.' );
    if( extIndex != string::npos ) {
        return filename.substr( 0, extIndex );
    } else {
        return filename;
    }
}

std::string basename_from_path( const string& Path ) {
    return basename_from_path_p( boost::filesystem::path( files::normalized_directory_name( Path ) ) );
}

std::wstring basename_from_path( const std::wstring& Path ) {
    return basename_from_path_p_w( boost::filesystem::path( files::normalized_directory_name( Path ) ) );
}

std::string extension_from_path( const string& Path ) {
    // Extract the extension
    size_t extIndex = Path.rfind( '.' );
    if( extIndex != string::npos ) {
        return Path.substr( extIndex );
    } else {
        return "";
    }
}

std::wstring extension_from_path( const std::wstring& path ) {
    return boost::filesystem::path( files::normalized_directory_name( path ) ).extension().wstring();
}

bool extension_iequals( const std::string& Path, const std::string& Ext ) {
    std::string pathExtLower = extension_from_path( Path );
    boost::algorithm::to_lower( pathExtLower );
    std::string extLower = Ext;
    boost::algorithm::to_lower( extLower );
    return pathExtLower == extLower;
}

bool extension_iequals( const std::wstring& Path, const std::wstring& Ext ) {
    std::wstring pathExtLower = extension_from_path( Path );
    boost::algorithm::to_lower( pathExtLower );
    std::wstring extLower = Ext;
    boost::algorithm::to_lower( extLower );
    return pathExtLower == extLower;
}

namespace {
// TODO: Roughly translated from latest boost tree, should just use it within boost when we switch to 1.60
fs::path lexically_normal( const fs::path& p ) {
    if( p.empty() )
        return p;

    fs::path temp;
    fs::path::iterator start( p.begin() );
    fs::path::iterator last( p.end() );
    fs::path::iterator stop( last-- );
    for( fs::path::iterator itr( start ); itr != stop; ++itr ) {
        // ignore "." except at start and last
        if( itr->native().size() == 1 && ( itr->native() )[0] == '.' && itr != start && itr != last )
            continue;

        // ignore a name and following ".."
        if( !temp.empty() && itr->native().size() == 2 && ( itr->native() )[0] == '.' &&
            ( itr->native() )[1] == '.' ) // dot dot
        {
            fs::path::string_type lf( temp.filename().native() );
            if( lf.size() > 0 && ( lf.size() != 1 || ( lf[0] != '.' && lf[0] != fs::path::preferred_separator ) ) &&
                ( lf.size() != 2 || ( lf[0] != '.' && lf[1] != '.'
#ifdef BOOST_WINDOWS_API
                                      && lf[1] != ':'
#endif
                                      ) ) ) {
                temp.remove_filename();
                //// if not root directory, must also remove "/" if any
                // if (temp.native().size() > 0
                //  && temp.native()[temp.native().size()-1]
                //    == separator)
                //{
                //  string_type::size_type rds(
                //    root_directory_start(temp.native(), temp.native().size()));
                //  if (rds == string_type::npos
                //    || rds != temp.native().size()-1)
                //  {
                //    temp.m_pathname.erase(temp.native().size()-1);
                //  }
                //}

                fs::path::iterator next( itr );
                if( temp.empty() && ++next != stop && next == last && *last == fs::path( _T( "." ) ) ) {
                    temp /= fs::path( _T( "." ) );
                }
                continue;
            }
        }

        temp /= *itr;
    };

    if( temp.empty() )
        temp /= fs::path( _T( "." ) );
    return temp;
}

// Case insensitive compare on Windows, case sensitive on Unix. This is in general a mess, and this
// rule is not "correct" by any means. Maybe it should be case insensitive on Unix. Maybe it should
// always use OS case normalization, then case sensitive compare. I'm making it like this right now,
// and if we need to we can change it...
bool compare_path_components( const fs::path& lhs, const fs::path& rhs ) {
#ifdef _WIN32
    return boost::algorithm::iequals( lhs.wstring(), rhs.wstring() );
#else
    return lhs == rhs;
#endif
}

// C++14 provide a mismatch algorithm with four iterator arguments(), but earlier
// standard libraries didn't, so provide this needed functionality.
inline std::pair<fs::path::iterator, fs::path::iterator> mismatch( fs::path::iterator it1, fs::path::iterator it1end,
                                                                   fs::path::iterator it2, fs::path::iterator it2end ) {
    for( ; it1 != it1end && it2 != it2end && compare_path_components( *it1, *it2 ); ) {
        ++it1;
        ++it2;
    }
    return std::make_pair( it1, it2 );
}

fs::path lexically_relative( const fs::path& p, const fs::path& base ) {
    std::pair<fs::path::iterator, fs::path::iterator> mm = mismatch( p.begin(), p.end(), base.begin(), base.end() );
    if( mm.first == p.begin() && mm.second == base.begin() )
        return fs::path();
    if( mm.first == p.end() && mm.second == base.end() )
        return fs::path( _T( "." ) );
    fs::path tmp;
    for( ; mm.second != base.end(); ++mm.second )
        tmp /= fs::path( _T( ".." ) );
    for( ; mm.first != p.end(); ++mm.first )
        tmp /= *mm.first;
    return tmp;
}
} // anonymous namespace

bool make_relative_path( const std::string& rootPath, std::string& inoutPath ) {
    fs::path rp = lexically_normal( rootPath ), p = lexically_normal( inoutPath );
    fs::path rel = lexically_relative( p, rp );
    if( !rel.empty() ) {
        inoutPath = rel.string();
        return true;
    } else {
        return false;
    }
}
bool make_relative_path( const std::wstring& rootPath, std::wstring& inoutPath ) {
    fs::path rp = lexically_normal( rootPath ), p = lexically_normal( inoutPath );
    fs::path rel = lexically_relative( p, rp );
    if( !rel.empty() ) {
        inoutPath = rel.wstring();
        return true;
    } else {
        return false;
    }
}

bool has_sequence_number( const std::string& path ) {
    string::size_type dotPos = path.rfind( '.' );
    if( dotPos != string::npos ) {
        // Find the last directory separator, to make sure the '.' character is in the actual filename
        string::size_type slashPos = path.rfind( '\\' );
        string::size_type forwardslashPos = path.rfind( '/' );
        if( forwardslashPos != string::npos ) {
            if( slashPos == string::npos )
                slashPos = forwardslashPos;
            else
                slashPos = ( std::max )( slashPos, forwardslashPos );
        }

        if( slashPos == string::npos || dotPos > slashPos ) {
            if( dotPos > 0 )
                return isdigit( path[dotPos - 1] ) != 0;
            else
                return false;
        }
    }

    if( path.empty() )
        return false;
    else
        return isdigit( path[path.size() - 1] ) != 0;
}

boost::filesystem::path replace_sequence_number( const boost::filesystem::path& thePath, long newNumber,
                                                 int numDigits ) {
    // Grab the filename
    frantic::tstring filename = to_tstring( thePath.leaf() );

    // Split it into the basename and the extension
    frantic::tstring basename, extension;
    size_t extIndex = filename.rfind( '.' );
    if( extIndex != string::npos ) {
        basename = filename.substr( 0, extIndex );
        extension = filename.substr( extIndex );
    } else {
        basename = filename;
    }

    // Now remove the existing digits from basename if required

    if( numDigits < 0 ) {
        int digitsIndex = (int)basename.size() - 1;
        numDigits = 0;
        while( digitsIndex >= 0 && isdigit( basename[digitsIndex] ) ) {
            --digitsIndex;
            numDigits++;
        }

        // Allow for negative numbers
        if( digitsIndex >= 0 && numDigits > 0 && basename[digitsIndex] == '-' ) {
            --digitsIndex;
            numDigits++;
        }

        if( numDigits > 0 ) {
            basename = basename.substr( 0, digitsIndex + 1 );
        } else {
            // Default to 4 digits
            numDigits = 4;
        }
    }

    // Generate the number and copy it into the destination string
    if( numDigits > 0 ) {
        basename += frantic::strings::zero_pad( newNumber, numDigits );
    }

    return thePath.branch_path() / ( basename + extension );
}

std::string replace_sequence_number( const string& thePath, long newNumber, int numDigits ) {
    return replace_sequence_number( boost::filesystem::path( files::normalized_directory_name( thePath ) ), newNumber,
                                    numDigits )
        .string();
}

std::wstring replace_sequence_number( const std::wstring& thePath, long newNumber, int numDigits ) {
    return replace_sequence_number( boost::filesystem::path( files::normalized_directory_name( thePath ) ), newNumber,
                                    numDigits )
        .wstring();
}

std::string increment_sequence_number( const std::string& thePath ) {
    return increment_sequence_number( boost::filesystem::path( files::normalized_directory_name( thePath ) ) ).string();
}

boost::filesystem::path increment_sequence_number( const boost::filesystem::path& thePath ) {
    return replace_sequence_number( thePath, extract_sequence_number( thePath ) + 1 );
}

int extract_sequence_number( const std::string& thePath ) {
    return extract_sequence_number( boost::filesystem::path( files::normalized_directory_name( thePath ) ) );
}

int extract_sequence_number( const std::wstring& thePath ) {
    return extract_sequence_number( boost::filesystem::path( files::normalized_directory_name( thePath ) ) );
}

namespace detail {

inline int isdigit( int c ) { return std::isdigit( c ); }
#if defined( _WIN32 ) || defined( _WIN64 ) ||                                                                          \
    defined( linux ) // apparently wint_t is the same as int on OSX, which results in a redefinition error
inline int isdigit( wint_t c ) { return std::iswdigit( c ); }
#endif
inline int atoi( const char* s ) { return ::atoi( s ); }

inline int atoi( const wchar_t* s ) {
    wchar_t* endPtr = 0;
    return static_cast<int>( std::wcstol( s, &endPtr, 10 ) );
}

template <class CharType>
void split_sequence_path( const std::basic_string<CharType>& thePath, std::basic_string<CharType>& outPrefix,
                          int& outDigitCount, int& outSequenceNumber, std::basic_string<CharType>& outPostfix ) {
    // Extract the extension
    typename std::basic_string<CharType>::size_type extIndex = thePath.rfind( '.' );
    if( extIndex != string::npos ) {
        outPostfix = thePath.substr( extIndex );
    } else {
        outPostfix.clear();
        extIndex = thePath.size();
    }

    // Now isolate the digits string
    int digitsIndex = (int)extIndex - 1;
    outDigitCount = 0;
    while( digitsIndex >= 0 && isdigit( thePath[digitsIndex] ) ) {
        --digitsIndex;
        outDigitCount++;
    }

    if( digitsIndex >= 0 && outDigitCount > 0 && thePath[digitsIndex] == '-' ) {
        --digitsIndex;
        outDigitCount++;
    }

    outSequenceNumber = atoi( thePath.substr( digitsIndex + 1, outDigitCount ).c_str() );
    outPrefix = thePath.substr( 0, digitsIndex + 1 );
}
} // namespace detail

// This splits an input path into three components, a prefix, sequence digit block (represented by the count of digits),
// and a postfix.
// A path can then be reassembled by outPrefix + frantic::strings::zero_pad(outFrameNumber, outDigitCount) + outPostfix
void split_sequence_path( const std::string& thePath, std::string& outPrefix, int& outDigitCount,
                          int& outSequenceNumber, std::string& outPostfix ) {
    detail::split_sequence_path( thePath, outPrefix, outDigitCount, outSequenceNumber, outPostfix );
}

void split_sequence_path( const boost::filesystem::path& thePath, boost::filesystem::path& outPrefix,
                          int& outDigitCount, int& outSequenceNumber, std::string& outPostfix ) {
    std::string prefix;
    split_sequence_path( thePath.string(), prefix, outDigitCount, outSequenceNumber, outPostfix );
    outPrefix = boost::filesystem::path( prefix );
}

void split_sequence_path( const std::wstring& thePath, std::wstring& outPrefix, int& outDigitCount,
                          int& outSequenceNumber, std::wstring& outPostfix ) {
    detail::split_sequence_path( thePath, outPrefix, outDigitCount, outSequenceNumber, outPostfix );
}

void split_sequence_path( const boost::filesystem::path& thePath, boost::filesystem::path& outPrefix,
                          int& outDigitCount, int& outSequenceNumber, std::wstring& outPostfix ) {
    std::wstring prefix;
    split_sequence_path( thePath.wstring(), prefix, outDigitCount, outSequenceNumber, outPostfix );
    outPrefix = boost::filesystem::path( prefix );
}

int extract_sequence_number( const boost::filesystem::path& thePath ) {
    // Grab the filename
    frantic::tstring filename = frantic::files::to_tstring( thePath.leaf() );

    // Remove the extension
    size_t extIndex = filename.rfind( '.' );
    if( extIndex != string::npos ) {
        filename = filename.substr( 0, extIndex );
    }

    // Now isolate the digits string
    int digitsIndex = (int)filename.size() - 1;
    int numDigits = 0;
    while( digitsIndex >= 0 && detail::isdigit( filename[digitsIndex] ) ) {
        --digitsIndex;
        numDigits++;
    }

    if( digitsIndex >= 0 && numDigits > 0 && filename[digitsIndex] == '-' ) {
        --digitsIndex;
        numDigits++;
    }

    if( numDigits > 0 ) {
        filename = filename.substr( digitsIndex + 1 );
    } else {
        throw std::runtime_error(
            "extract_sequence_number: Unsuccessfully tried to extract a sequence number from the filename \"" +
            thePath.string() + "\"" );
    }

    return detail::atoi( filename.c_str() );
}

frantic::tstring to_unix_style_search_path( const frantic::tstring& searchPath ) {
    std::vector<frantic::tstring> splitPaths;
    frantic::tstring txt( _T("") );
    strings::split_search_path( searchPath, splitPaths );
    for( unsigned int i = 0; i < splitPaths.size(); ++i ) {
        if( i < splitPaths.size() - 1 )
            txt += splitPaths[i] + _T(":");
        else
            txt += splitPaths[i];
    }

    return txt;
}

std::vector<frantic::tstring> get_program_files_folders() {
    std::vector<frantic::tstring> result;
#if defined( _WIN32 )
    TCHAR buffer[MAX_PATH];
#if defined( CSIDL_PROGRAM_FILES )
    if( SUCCEEDED( SHGetFolderPath( NULL, CSIDL_PROGRAM_FILES, NULL, 0, buffer ) ) )
        result.push_back( files::ensure_trailing_pathseparator( buffer ) );
#endif
#if defined( CSIDL_PROGRAM_FILESX86 )
    if( SUCCEEDED( SHGetFolderPath( NULL, CSIDL_PROGRAM_FILESX86, NULL, 0, buffer ) ) )
        result.push_back( files::ensure_trailing_pathseparator( buffer ) );
#endif
#endif
    return result;
}

// the following two functions are not valid on Linux
#if defined( _WIN32 ) || defined( _WIN64 )
frantic::tstring get_application_data_folder( bool allUsers ) {
    int flags = CSIDL_FLAG_CREATE;
    if( allUsers )
        flags |= CSIDL_COMMON_APPDATA;
    else
        flags |= CSIDL_APPDATA;

    TCHAR buffer[MAX_PATH];
    if( SUCCEEDED( SHGetFolderPath( NULL, flags, NULL, 0, buffer ) ) )
        return files::ensure_trailing_pathseparator( buffer );

    return _T("");
}

frantic::tstring get_local_application_data_folder() {
    TCHAR buffer[MAX_PATH];
    if( SUCCEEDED( SHGetFolderPath( NULL, CSIDL_FLAG_CREATE | CSIDL_LOCAL_APPDATA, NULL, 0, buffer ) ) )
        return files::ensure_trailing_pathseparator( buffer );

    return _T("");
}
#endif

#if( defined( _MSC_VER ) && _MSC_VER > 1200 ) ||                                                                       \
    defined( linux ) // Regex is only defined in visual studio 7 or higher and Linux
// Parses the filename to find the tile information.  The form is like
// ...._TILE_#x#_#x#... where the first pair of numbers are the 1-based tile coordinate, and
// the second pair are the # of tiles along x and y.
void get_tile_from_filename( const std::string& filename, int& outX, int& outY, int& outWidth, int& outHeight ) {
    static boost::regex findTile( "_tile_(\\d+)x(\\d+)_(\\d+)x(\\d+)_", boost::regex::perl | boost::regex::icase );

    boost::smatch what;
    if( boost::regex_search( filename, what, findTile ) ) {
        outX = boost::lexical_cast<int>( what[1].str() );
        outY = boost::lexical_cast<int>( what[2].str() );
        outWidth = boost::lexical_cast<int>( what[3].str() );
        outHeight = boost::lexical_cast<int>( what[4].str() );
    } else {
        throw std::runtime_error( "get_tile_from_filename: Filename \"" + filename +
                                  "\" did not contain a _tile_ tag." );
    }
}

std::string replace_tile_in_filename( const std::string& filename, int x, int y ) {
    static boost::regex replaceTile( "_tile_\\d+x\\d+_", boost::regex::perl | boost::regex::icase );
    return boost::regex_replace( filename, replaceTile,
                                 "_tile_" + boost::lexical_cast<std::string>( x ) + "x" +
                                     boost::lexical_cast<std::string>( y ) + "_" );
}

std::string remove_tile_tag_from_filename( const std::string& filename ) {
    static boost::regex eraseTile( "_tile_\\d+x\\d+_\\d+x\\d+_", boost::regex::perl | boost::regex::icase );
    return boost::regex_replace( filename, eraseTile, "" );
}

void get_part_from_filename( const frantic::tstring& filename, int& outIndex, int& outCount ) {
    static boost::basic_regex<frantic::tchar> findPart( _T(".*_part(\\d+)of(\\d+)_"),
                                                        boost::regex::perl | boost::regex::icase );

    boost::match_results<frantic::tstring::const_iterator> what;
    if( boost::regex_search( filename, what, findPart ) ) {
        outIndex = boost::lexical_cast<int>( what[1].str() );
        outCount = boost::lexical_cast<int>( what[2].str() );
    } else {
        throw std::runtime_error( "get_part_from_filename: Filename \"" + frantic::strings::to_string( filename ) +
                                  "\" did not contain a valid _part tag." );
    }
}

frantic::tstring replace_part_in_filename( const frantic::tstring& filename, int index ) {
    static boost::basic_regex<frantic::tchar> findPart( _T("^(.*_part)(\\d+)of(\\d+)(_.*)$"),
                                                        boost::regex::perl | boost::regex::icase );

    boost::match_results<frantic::tstring::const_iterator> what;
    if( boost::regex_search( filename, what, findPart ) ) {
        if( what[2].str().size() == what[3].str().size() ) {
            frantic::tstring part = frantic::strings::zero_pad( index, (unsigned)what[3].str().size() );
            return what[1].str() + part + _T("of") + what[3].str() + what[4].str();
        } else {
            frantic::tstring part = boost::lexical_cast<frantic::tstring>( index );
            return what[1].str() + part + _T("of") + what[3].str() + what[4].str();
        }
    } else {
        throw std::runtime_error( "replace_part_in_filename: Filename \"" + frantic::strings::to_string( filename ) +
                                  "\" did not contain a valid _part tag." );
    }
}
#else
void get_tile_from_filename( const std::string& filename, int& outX, int& outY, int& outWidth, int& outHeight ) {
    int tileIndex = filename.find( "_tile_" );
    if( tileIndex >= 0 ) {
        int outXIndex = filename.find( "_", tileIndex + 1 );
        int outYIndex = filename.find( "x", outXIndex + 1 );
        int outWidthIndex = filename.find( "_", outYIndex + 1 );
        int outHeightIndex = filename.find( "x", outWidthIndex + 1 );
        int lastIndex = filename.find( "_", outHeightIndex + 1 );

        outX = boost::lexical_cast<int>( filename.substr( outXIndex + 1, outYIndex - outXIndex - 1 ) );
        outY = boost::lexical_cast<int>( filename.substr( outYIndex + 1, outWidthIndex - outYIndex - 1 ) );
        outWidth = boost::lexical_cast<int>( filename.substr( outWidthIndex + 1, outHeightIndex - outWidthIndex - 1 ) );
        outHeight = boost::lexical_cast<int>( filename.substr( outHeightIndex + 1, lastIndex - outHeightIndex - 1 ) );
    } else
        throw std::runtime_error( "get_tile_from_filename: Filename \"" + filename +
                                  "\" did not contain a _tile_ tag." );
}

std::string replace_tile_in_filename( const std::string& filename, int x, int y ) {
    int tileIndex = filename.find( "_tile_" );
    if( tileIndex >= 0 ) {
        int outXIndex = filename.find( "_", tileIndex + 1 );
        int outYIndex = filename.find( "x", outXIndex + 1 );
        int outWidthIndex = filename.find( "_", outYIndex + 1 );
        std::string tempFilename = filename;
        return tempFilename.replace( tileIndex, outWidthIndex - tileIndex + 1,
                                     "_tile_" + boost::lexical_cast<std::string>( x ) + "x" +
                                         boost::lexical_cast<std::string>( y ) + "_" );
    }
    return filename;
}

std::string remove_tile_tag_from_filename( const std::string& filename ) {
    int tileIndex = filename.find( "_tile_" );
    if( tileIndex >= 0 ) {
        int outXIndex = filename.find( "_", tileIndex + 1 );
        int outYIndex = filename.find( "x", outXIndex + 1 );
        int outWidthIndex = filename.find( "_", outYIndex + 1 );
        std::string tempFilename = filename;
        return tempFilename.replace( tileIndex, outWidthIndex - tileIndex + 1, "" );
    }
    return filename;
}
#endif

frantic::tstring get_universal_name( const frantic::tstring& path ) {
#if defined( _WIN32 ) || defined( _WIN64 )

    if( path.length() < 2 || path[1] != ':' )
        return path;

    DWORD dwRetVal;
    WCHAR Buffer[1024];
    DWORD dwBufferLength = 1024;
    UNIVERSAL_NAME_INFO* unameinfo;
    // REMOTE_NAME_INFO *remotenameinfo;

    unameinfo = (UNIVERSAL_NAME_INFO*)&Buffer;
    dwRetVal = WNetGetUniversalName( path.c_str(), UNIVERSAL_NAME_INFO_LEVEL, (LPVOID)unameinfo, &dwBufferLength );

    if( dwRetVal == NO_ERROR ) {
        return frantic::tstring( unameinfo->lpUniversalName );
    } else {
        std::string errStr0( "get_universal_name:  Failed for InfoLevel=UNIVERSAL_NAME_INFO_LEVEL with error: " );
        std::string errStr1 = errStr0 + boost::lexical_cast<std::string>( dwRetVal );
        throw std::runtime_error( errStr1.c_str() );
    }
#else
    throw std::runtime_error( "files::get_universal_name():  This method is not yet implemented!" );
#endif
}

frantic::tstring to_tstring( const boost::filesystem::path& path ) {
#ifdef FRANTIC_USE_WCHAR
    return path.wstring();
#else
    return path.string();
#endif
}
} // namespace files
} // namespace frantic
