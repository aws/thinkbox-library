// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/strings/tstring.hpp>
#include <string>

#include <boost/filesystem/path.hpp>

// The following is needed by  function.
#if defined( _WIN32 ) || defined( _WIN64 )
#pragma comment( lib, "shell32.lib" ) // needed for SHGetFolderPath() in get_local_application_data_folder()
#pragma comment( lib, "Mpr.lib" )     // needed for WNetGetUniversalName() in get_universal_name()
#endif

namespace frantic {
namespace files {

std::string forward_slashes( std::string Path );
std::string backward_slashes( std::string Path );
std::wstring forward_slashes( std::wstring Path );
std::wstring backward_slashes( std::wstring Path );
// Ensures '/' for unix, '\\' for windows
std::string ensure_trailing_pathseparator( std::string Path );
std::wstring ensure_trailing_pathseparator( std::wstring Path );

// Fixes paths so they have the proper path separators for DOS paths
std::string normalized_directory_name( const std::string& name );
std::wstring normalized_directory_name( const std::wstring& name );

// removes slashes, spaces, dots, etc from a string. Offending characters will be
// replaced with the replacement string.
std::string clean_directory_name( std::string name, std::string replacement = std::string() );

std::string replace_directory( const std::string& Path, std::string ReplacementDir );
std::string replace_filename( const std::string& Path, const std::string& ReplacementName );
std::wstring replace_filename( const std::wstring& Path, const std::wstring& ReplacementName );
// ReplacementExt should be like ".ext"
std::string replace_extension( const std::string& Path, const std::string& ReplacementExt );
boost::filesystem::path replace_extension_p( const boost::filesystem::path& Path, const std::string& ReplacementExt );
std::wstring replace_extension( const std::wstring& Path, const std::wstring& ReplacementExt );
boost::filesystem::path replace_extension_p( const boost::filesystem::path& Path, const std::wstring& ReplacementExt );

// If path is "c:/winnt/system.dll"

// this returns "c:/winnt/"
std::string directory_from_path( const std::string& Path );
std::wstring directory_from_path( const std::wstring& Path );

// "system.dll"
std::string filename_from_path( const std::string& Path );
std::wstring filename_from_path( const std::wstring& Path );

// "system"
std::string basename_from_path( const std::string& Path );
std::wstring basename_from_path( const std::wstring& Path );
std::string basename_from_path_p( const boost::filesystem::path& thePath );

// ".dll"
std::string extension_from_path( const std::string& Path );
std::wstring extension_from_path( const std::wstring& Path );

// Returns true if the extension of `Path` is equal to `Ext` with
// a case-insensitive comparison.
bool extension_iequals( const std::string& Path, const std::string& Ext );
bool extension_iequals( const std::wstring& Path, const std::wstring& Ext );

/**
 * If `inoutPath` and `rootPath` share the same root, creates a relative path
 * such that `rootPath / relPath` gives back the original `inoutPath`,
 * and return true. Otherwise leave `inoutPath` unchanged, and return false.
 *
 * Example:
 *    string rootPath = "C:\\MyDir\\Subdir", path = "C:\\MyDir\\Other\\File.txt";
 *    assert( make_relative_path(rootPath, path) );
 *    assert( path == "..\\Other\\File.txt" );
 */
bool make_relative_path( const std::string& rootPath, std::string& inoutPath );
bool make_relative_path( const std::wstring& rootPath, std::wstring& inoutPath );

// Replaces the # in a filename like "output0003.tif"
// If numDigits is matchExistingSequence, the output of ReplaceSequenceNumber( inputstring, 192 ) is
//		output.tif		->	output0192.tif  (Defaults to 4 digits when there are none in the given
// path) 		output3.tif		->	output2.tif 		output2732.tif	->	output0192.tif
// If numDigits is >= 0, then it will append that number of digits before the extension. numDigits == 4:
//		output.tif		->  output0192.tif
//      ouput3.tif		->  output30192.tif
enum { matchExistingSequence = -1 };
std::string replace_sequence_number( const std::string& path, long newNumber,
                                     int numDigits = files::matchExistingSequence );
inline std::string replace_sequence_number( const char* path, long newNumber,
                                            int numDigits = files::matchExistingSequence ) {
    return replace_sequence_number( std::string( path ), newNumber, numDigits );
}
std::wstring replace_sequence_number( const std::wstring& path, long newNumber,
                                      int numDigits = files::matchExistingSequence );
boost::filesystem::path replace_sequence_number( const boost::filesystem::path& path, long newNumber,
                                                 int numDigits = files::matchExistingSequence );

int extract_sequence_number( const boost::filesystem::path& path );
int extract_sequence_number( const std::string& path );
int extract_sequence_number( const std::wstring& path );

/**
 * This function checks whether the given filename has a sequence number.
 *
 * @param  path  The filename (possibly including the full path) to test for a sequence number.
 * @return  True if the filename has a sequence number, false otherwise.
 */
bool has_sequence_number( const std::string& path );

std::string increment_sequence_number( const std::string& path );
boost::filesystem::path increment_sequence_number( const boost::filesystem::path& path );

void split_sequence_path( const std::string& path, std::string& outPrefix, int& outDigitCount, int& outSequenceNumber,
                          std::string& outPostfix );
void split_sequence_path( const boost::filesystem::path& path, boost::filesystem::path& outPrefix, int& outDigitCount,
                          int& outSequenceNumber, std::string& outPostfix );

void split_sequence_path( const std::wstring& path, std::wstring& outPrefix, int& outDigitCount, int& outSequenceNumber,
                          std::wstring& outPostfix );
void split_sequence_path( const boost::filesystem::path& path, boost::filesystem::path& outPrefix, int& outDigitCount,
                          int& outSequenceNumber, std::wstring& outPostfix );

std::string to_unix_style_search_path( const std::string& searchPath );

// Returns the 32-bit and 64-bit program files folders, as they're available
std::vector<frantic::tstring> get_program_files_folders();

frantic::tstring get_application_data_folder( bool allUsers );
frantic::tstring get_local_application_data_folder();

// cleans up a path name, removing ", trimming whitespaces, converting all slashes to backslashes and normalizing the
// path.
// Good for commandline arguments that have to play nice with boost::filesystem::path. Not very pretty or efficient.
std::string strict_path_name( std::string const& path );
std::wstring strict_path_name( std::wstring const& path );

namespace detail {
// Adds a string before the sequence number in the given filename
// For example add_before_sequence_number( "test0000.tif, "_diffuse" ) produces "test_diffuse0000.tif"
template <class CharType>
inline std::basic_string<CharType> add_before_sequence_number( const std::basic_string<CharType>& path,
                                                               const std::basic_string<CharType>& addition ) {
    typename std::basic_string<CharType>::size_type prefixPosition = path.rfind( '.' );

    if( prefixPosition == std::basic_string<CharType>::npos && !isdigit( path[path.size() - 1] ) )
        return path + addition;

    if( prefixPosition == std::basic_string<CharType>::npos )
        prefixPosition = path.size() - 1;

    while( prefixPosition > 0 && isdigit( path[prefixPosition - 1] ) )
        --prefixPosition;

    // If we found the subframe divider (',') move past it
    if( prefixPosition > 1 && path[prefixPosition - 1] == ',' && isdigit( path[prefixPosition - 2] ) ) {

        --prefixPosition;

        while( prefixPosition > 0 && isdigit( path[prefixPosition - 1] ) )
            --prefixPosition;
    }

    return path.substr( 0, prefixPosition ) + addition + path.substr( prefixPosition );
}
} // namespace detail

// Adds a string before the sequence number in the given filename
// For example add_before_sequence_number( L"test0000.tif, L"_diffuse" ) produces L"test_diffuse0000.tif"
inline std::wstring add_before_sequence_number( const std::wstring& path, const std::wstring& addition ) {
    return detail::add_before_sequence_number<wchar_t>( path, addition );
}

#ifndef FRANTIC_USE_WCHAR
// Adds a string before the sequence number in the given filename
// For example add_before_sequence_number( "test0000.tif, "_diffuse" ) produces "test_diffuse0000.tif"
inline std::string add_before_sequence_number( const std::string& path, const std::string& addition ) {
    return detail::add_before_sequence_number<char>( path, addition );
}
#endif

// Parses the filename to find the tile information.  The form is like
// ...._TILE_#x#_#x#... where the first pair of numbers are the 1-based tile coordinate, and
// the second pair are the # of tiles along x and y.
void get_tile_from_filename( const std::string& filename, int& outX, int& outY, int& outWidth, int& outHeight );
std::string replace_tile_in_filename( const std::string& filename, int x, int y );
std::string remove_tile_tag_from_filename( const std::string& filename );

// Parses the filename to find the partition information.  The form is like
// ..._part#of#_... where the first number is the index and the second is the count.
void get_part_from_filename( const frantic::tstring& filename, int& outIndex, int& outCount );
frantic::tstring replace_part_in_filename( const frantic::tstring& filename, int index );

// Returns the UNC equivalent of the supplied path
frantic::tstring get_universal_name( const frantic::tstring& path );

frantic::tstring to_tstring( const boost::filesystem::path& path );
} // namespace files
} // namespace frantic
