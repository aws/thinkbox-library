// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#if defined( _WIN32 ) || defined( _WIN64 )
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <frantic/win32/utility.hpp>
#else
#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#endif

#include <fstream>
#include <iostream>
#include <string>

#include <boost/filesystem/operations.hpp>

#include <frantic/files/paths.hpp>
#include <frantic/misc/string_functions.hpp>
#include <frantic/os/environment.hpp>
#include <frantic/strings/tstring.hpp>

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif

namespace frantic {
namespace files {

/**
 * This provides a wrapper for a FILE* object. The constructor opens the file and the destructor closes it.
 * This should used to help ensure exception safety when working with FILE*s.
 */
class file_ptr {
    FILE* m_file;

    // Disable copy construction/assignment
    file_ptr( const file_ptr& );            // do not implement
    file_ptr& operator=( const file_ptr& ); // do not implement
  public:
    /**
     * Default constructor.
     */
    file_ptr()
        : m_file( 0 ) {}

    /**
     * Constructor accepts a FILE* to initialize the class.
     *
     * @param file - the FILE struct opened with std::fopen
     */
    file_ptr( FILE* file )
        : m_file( file ) {}

    /**
     * The file object is closed on destruction
     */
    ~file_ptr() { close(); }

    /**
     * The reset function closes any existing open file, and replaces its FILE* with the given variable.
     *
     * @param file - the FILE struct opened with std::fopen
     */
    void reset( FILE* file ) {
        close();
        m_file = file;
    }

    int close() {
        int result = 0;
        if( m_file ) {
            result = fclose( m_file );
            m_file = 0;
        }
        return result;
    }

    /**
     * Provides the file pointer
     */
    FILE* get() { return m_file; }

    const FILE* get() const { return m_file; }

    operator FILE*() { return m_file; }

    operator const FILE*() const { return m_file; }

    FILE* operator->() { return m_file; }
};

/**
 * file_exists() - Tests to see if the specified file exists
 *
 * @param filename - The full path of the file to be detected.
 */
bool file_exists( const std::string& filename );

bool file_exists( const std::wstring& filename );

/*
 * delete_file() deletes the specified filename.
 *
 * @param filename The full path of the file to be deleted.
 * @return        void
 */
void delete_file( const std::string& filename );

/*
 * delete_file() deletes the specified filename.
 *
 * @param filename The full path of the file to be deleted.
 * @return        void
 */
void delete_file( const std::wstring& filename );

#if defined( _WIN32 ) || defined( _WIN64 )
/**
 * Rename the specified file.
 *
 * @param from The name of the file to rename.
 * @param to The new name for the file.
 */
int rename_file( const std::wstring& from, const std::wstring& to );
#endif

/**
 * Rename the specified file.
 *
 * @param from The name of the file to rename.
 * @param to The new name for the file.
 */
int rename_file( const std::string& from, const std::string& to );

#if defined( _WIN32 ) || defined( _WIN64 )
bool set_hidden_attribute( const std::string& filename, bool hide );

std::string file_modification_time_iso( const char* filename );

bool are_modified_times_equal( std::string filename1, std::string filename2 );

std::string file_modification_time_iso( const std::string& filename );

bool directory_exists( const std::string& dirname );

bool directory_exists( const std::wstring& dirname );
#else
bool directory_exists( const std::string& dirname );
#endif

// This function reates a directory by recursively creating the chain of directories necessary
void make_directory( const std::string& directory );
#if defined( _WIN32 ) || defined( _WIN64 )
void make_directory( const std::wstring& directory );
#endif

int fseek64( FILE* stream, boost::int64_t offset, int origin );

boost::int64_t ftell64( FILE* stream );

boost::int64_t file_size( const std::string& filename );

#if defined( _WIN32 ) || defined( _WIN64 )
bool files_different( WIN32_FIND_DATA& wfd1, WIN32_FIND_DATA& wfd2 );

bool files_different( WIN32_FIND_DATA& wfd1, WIN32_FILE_ATTRIBUTE_DATA& wfad2 );

bool files_different( WIN32_FILE_ATTRIBUTE_DATA& wfad1, WIN32_FILE_ATTRIBUTE_DATA& wfad2 );

bool files_different( WIN32_FIND_DATA& wfd1, const frantic::tstring& fileName );

// Returns true if the last modified time and the file sizes match, false otherwise
bool files_different( const frantic::tstring& first, const frantic::tstring& second );
#else
bool files_different( const std::string& /*first*/, const std::string& /*second*/ );
#endif

bool files_different( const boost::filesystem::path& first, const boost::filesystem::path& second );

std::string first_available_numbered_file( const std::string& fileName, int digits );

std::wstring first_available_numbered_file( const std::wstring& fileName, int digits );

std::string first_available_numbered_file( const std::string& fileName );

void rename_to_backup( const std::string& file, const std::string& backupTag );

#if defined( _WIN32 ) || defined( _WIN64 )
void rename_to_backup( const std::wstring& file, const std::wstring& backupTag );
#endif

// Returns the first file in a semicolon separated list of files that exists.  Returns
// an empty string if nothing matched
frantic::tstring search_file_list( const frantic::tstring& semicolon_file_list );

#ifdef _WIN32
frantic::tstring search_file_list_for_64bit( const frantic::tstring& semicolon_file_list );

frantic::tstring search_file_list_for_32bit( const frantic::tstring& semicolon_file_list );
#endif

// Returns the first directory in a semicolon separated list of files that exists.  Returns
// an empty string if nothing matched
frantic::tstring search_directory_list_for_directory( const frantic::tstring& semicolon_directory_list );

frantic::tstring find_file_in_search_path( const frantic::tstring& fileName, const frantic::tstring& searchPath );

/**
 * search_directory_list() - Returns the first file in a semicolon separated list of files that exists.
 *
 * @param semicolon_directory_list - the semicolon seperated list of directories
 * @param file - the filename to search
 * @returns the first matching file (with path) found or an empty string if nothing matched.
 */
frantic::tstring search_directory_list( const frantic::tstring& semicolon_directory_list,
                                        const frantic::tstring& file );

frantic::tstring search_path( const frantic::tstring& file );

/**
 * add_to_path() adds a directory to the environment variable PATH
 *
 * @param newDirectory - the full path of the directory to add
 */
void add_to_path( const frantic::tstring& newDirectory );

/*
 * get_filenames_in_directory() retrieves a listing of unhidden files in the
 * specified directory.  The listing does not include sub directories.
 * An exception is thrown if the directory does not exist.
 *
 * @param directory The directory fromm which to get the filenames.
 * @param fileListing An out parameter to hold the filenames.
 * @return        void
 */
void get_filenames_in_directory( const frantic::tstring& directory, std::vector<frantic::tstring>& fileListing );

class file_open_error : public std::runtime_error {
  public:
    file_open_error( const std::string& message )
        : std::runtime_error( message ) {}
};

/**
 * Create a temporary file for `mainOutputFile` in `directory`.
 * @param directory The directory to create the file in.
 * @param prefix The prefix to use for the temporary file name.
 */
inline boost::filesystem::path create_temp_path( const boost::filesystem::path& directory,
                                                 const boost::filesystem::path& prefix ) {
    frantic::tstring out = frantic::files::to_tstring( prefix.filename() );
    const boost::filesystem::path tempPath(
        boost::filesystem::unique_path( directory / ( out + _T(".%%%%-%%%%-%%%%-%%%%.tmp") ) ) );
    return tempPath;
}

#if defined( FRANTIC_USE_WCHAR ) && ( defined( _WIN32 ) || defined( _WIN64 ) )
inline FILE* tfopen( const wchar_t* filename, const wchar_t* mode ) { return ::_wfopen( filename, mode ); }
#else
inline FILE* tfopen( const char* filename, const char* mode ) { return ::fopen( filename, mode ); }
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
inline FILE* fopen( const boost::filesystem::path& filename, const char* mode ) {
    return ::_wfopen( filename.c_str(), frantic::strings::to_wstring( mode ).c_str() );
}
#else
inline FILE* fopen( const boost::filesystem::path& filename, const char* mode ) {
    return ::fopen( filename.c_str(), mode );
}
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
inline int remove_file( const wchar_t* filename ) { return ::_wremove( filename ); }
#endif

inline int remove_file( const char* filename ) { return std::remove( filename ); }
} // namespace files
} // namespace frantic

#if defined( _MSC_VER )
#pragma warning( pop )
#endif
