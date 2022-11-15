// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// This file was created to avoid putting the DWORD typedef in the global namespace for those that include files.hpp. It
// was conflicting with other headers.
// The rest of the code remains in files.hpp.

// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/files/files.hpp>

#ifndef _WIN32
typedef unsigned char BYTE;
typedef unsigned int UINT;
typedef unsigned long ULONG;
typedef ULONG DWORD;
#endif

using namespace frantic;
namespace fs = boost::filesystem;

namespace {

#if defined( _WIN32 ) || defined( _WIN64 )
inline bool ffGetFileAttributesEx( const std::string& filename, GET_FILEEX_INFO_LEVELS fInfoLevelId,
                                   LPVOID lpFileInformation ) {
    return GetFileAttributesExA( filename.c_str(), fInfoLevelId, lpFileInformation ) != 0;
}

inline bool ffGetFileAttributesEx( const std::wstring& filename, GET_FILEEX_INFO_LEVELS fInfoLevelId,
                                   LPVOID lpFileInformation ) {
    return GetFileAttributesExW( filename.c_str(), fInfoLevelId, lpFileInformation ) != 0;
}
#endif

template <typename CharType>
bool ff_file_exists( const std::basic_string<CharType>& filename ) {
#ifdef _WIN32
    WIN32_FILE_ATTRIBUTE_DATA fileAD;
    memset( &fileAD, 0, sizeof( fileAD ) );

    bool result = true;

    if( !ffGetFileAttributesEx( filename, GetFileExInfoStandard, &fileAD ) ) {
        // std::cout << "File \"" << fileName << "\" does not exist because " << win32::GetLastErrorMessage() <<
        // std::endl;
        result = false;
        if( GetLastError() != ERROR_FILE_NOT_FOUND && GetLastError() != ERROR_PATH_NOT_FOUND &&
            GetLastError() != ERROR_BAD_NETPATH && GetLastError() != ERROR_BAD_NET_NAME &&
            GetLastError() != ERROR_INVALID_NAME ) {
            throw std::runtime_error( "file_exists: Could not get the file attributes of \"" +
                                      frantic::strings::to_string( filename ) +
                                      "\": " + win32::GetLastErrorMessageA() );
        }
    } else if( fileAD.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY ) {
        result = false;
    }

    return result;
#else
    // The following did not work consistently on Win32, Boost 1.32.0...
    // TODO It might be reasonable to try this with the newer boost. Though for now, the above is suitable.
    fs::path thePath = fs::path( filename );
    return fs::exists( thePath ) && !fs::is_directory( thePath );
#endif
}
} // anonymous namespace

bool files::file_exists( const std::string& filename ) { return ff_file_exists( filename ); }

bool files::file_exists( const std::wstring& filename ) { return ff_file_exists( filename ); }

void files::delete_file( const std::string& filename ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    DeleteFileA( filename.c_str() );
#else
    fs::remove( fs::path( filename ) );
#endif
}

void files::delete_file( const std::wstring& filename ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    DeleteFileW( filename.c_str() );
#else
    fs::remove( fs::path( filename ) );
#endif
}

#if defined( _WIN32 ) || defined( _WIN64 )
int files::rename_file( const std::wstring& from, const std::wstring& to ) {
    return ::_wrename( from.c_str(), to.c_str() );
}
#endif

int files::rename_file( const std::string& from, const std::string& to ) {
    return ::rename( from.c_str(), to.c_str() );
}

#if defined( _WIN32 ) || defined( _WIN64 )
bool files::set_hidden_attribute( const std::string& filename, bool hide ) {
    if( file_exists( filename ) ) {
        DWORD fileAttributes = GetFileAttributesA( filename.c_str() );
        if( hide ) {
            fileAttributes = fileAttributes | FILE_ATTRIBUTE_HIDDEN;
        } else {
            fileAttributes = fileAttributes & ( ~FILE_ATTRIBUTE_HIDDEN );
        }
        return SetFileAttributesA( filename.c_str(), fileAttributes ) != 0;
    } else
        return false;
}

std::string files::file_modification_time_iso( const char* filename ) {
    fs::path filePath( filename );
    if( fs::exists( filePath ) ) {
        time_t tUnstructured = fs::last_write_time( filePath );

        tm* t = localtime( &tUnstructured );
        char buffer[64];
        // tm_year is the year minus 1900
        // tm_mon is the month where 0 = Jan.
        sprintf( buffer, "%d-%02d-%02d %02d:%02d:%02d%s", t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour,
                 t->tm_min, t->tm_sec, frantic::win32::getTimeZoneOffset().c_str() );
        return buffer;
    } else {
        return "";
    }
}

bool files::are_modified_times_equal( std::string filename1, std::string filename2 ) {
    bool result = false;
    fs::path file1Path( filename1 );
    fs::path file2Path( filename2 );
    if( fs::exists( file1Path ) && fs::exists( file2Path ) ) {
        result = ( fs::last_write_time( file1Path ) == fs::last_write_time( file2Path ) );
    }
    return result;
}

std::string files::file_modification_time_iso( const std::string& filename ) {
    fs::path filePath( filename );
    if( fs::exists( filePath ) ) {
        time_t tUnstructured = fs::last_write_time( filePath );

        tm* t = localtime( &tUnstructured );
        char buffer[64];
        // tm_year is the year minus 1900
        // tm_mon is the month where 0 = Jan.
        sprintf( buffer, "%d-%02d-%02d %02d:%02d:%02d%s", t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour,
                 t->tm_min, t->tm_sec, frantic::win32::getTimeZoneOffset().c_str() );
        return buffer;
    } else {
        return "";
    }
}

bool files::directory_exists( const std::string& dirname ) {
    fs::path p( dirname );
    return fs::is_directory( p );
}

bool files::directory_exists( const std::wstring& dirname ) {
    fs::path p( dirname );
    return fs::is_directory( p );
}
#else
bool files::directory_exists( const std::string& dirname ) {
    DIR* dir = opendir( dirname.c_str() );
    if( dir == 0 )
        return false;
    else {
        closedir( dir );
        return true;
    }
}
#endif

namespace {
#if defined( _WIN32 ) || defined( _WIN64 )
inline void ffCreateDirectory( const std::string& name, LPSECURITY_ATTRIBUTES lpSecurityAttributes ) {
    CreateDirectoryA( name.c_str(), lpSecurityAttributes );
}
inline void ffCreateDirectory( const std::wstring& name, LPSECURITY_ATTRIBUTES lpSecurityAttributes ) {
    CreateDirectoryW( name.c_str(), lpSecurityAttributes );
}
#endif

template <class CharType>
void ff_make_directory( std::basic_string<CharType> directory ) {
    if( directory.size() == 0 )
        return;

    if( directory[directory.size() - 1] == '\\' || directory[directory.size() - 1] == '/' )
        directory.resize( directory.size() - 1 );

    if( files::directory_from_path( directory ).size() >= directory.size() )
        return;

    // if( files::directory_exists( directory ) )
    //	return;

    files::make_directory( files::directory_from_path( directory ) );

#if defined( _WIN32 ) || defined( _WIN64 )
    ffCreateDirectory( directory, NULL );
#else
    mkdir( directory.c_str(), 0777 ); // not 100% sure on the mode_t value to use
#endif
}

} // anonymous namespace

// This function reates a directory by recursively creating the chain of directories necessary
void files::make_directory( const std::string& directory ) { ff_make_directory( directory ); }
#if defined( _WIN32 ) || defined( _WIN64 )
void files::make_directory( const std::wstring& directory ) { ff_make_directory( directory ); }
#endif

int files::fseek64( FILE* stream, boost::int64_t offset, int origin ) {
#ifdef _WIN32
    return _fseeki64( stream, offset, origin );
#elif defined( __APPLE__ )
    return fseeko( stream, offset, origin );
#else
    return fseeko64( stream, offset, origin );
#endif
}

boost::int64_t files::ftell64( FILE* stream ) {
#ifdef _WIN32
    return _ftelli64( stream );
#elif defined( __APPLE__ )
    return ftello( stream );
#else
    return ftello64( stream );
#endif
}

boost::int64_t files::file_size( const std::string& filename ) {
    if( file_exists( filename ) ) {
        std::ifstream ftest( filename.c_str(), std::ios::binary );
        if( !ftest )
            return -1;
        else {
            ftest.seekg( 0, std::ios::end );
            return ftest.tellg();
        }
    } else {
        return -1;
    }
}

#if defined( _WIN32 ) || defined( _WIN64 )
bool files::files_different( WIN32_FIND_DATA& wfd1, WIN32_FIND_DATA& wfd2 ) {
    if( wfd1.nFileSizeLow != wfd2.nFileSizeLow || wfd1.nFileSizeHigh != wfd2.nFileSizeHigh )
        return true;

    if( wfd1.ftLastWriteTime.dwLowDateTime != wfd2.ftLastWriteTime.dwLowDateTime ||
        wfd1.ftLastWriteTime.dwHighDateTime != wfd2.ftLastWriteTime.dwHighDateTime )
        return true;

    return false;
}

bool files::files_different( WIN32_FIND_DATA& wfd1, WIN32_FILE_ATTRIBUTE_DATA& wfad2 ) {
    bool result = false;

    if( wfd1.nFileSizeLow != wfad2.nFileSizeLow || wfd1.nFileSizeHigh != wfad2.nFileSizeHigh )
        result = true;

    if( wfd1.ftLastWriteTime.dwLowDateTime != wfad2.ftLastWriteTime.dwLowDateTime ||
        wfd1.ftLastWriteTime.dwHighDateTime != wfad2.ftLastWriteTime.dwHighDateTime )
        result = true;

    //*
    if( result ) {
        std::cout << "files_different: File Size A: " << wfd1.nFileSizeLow << ", " << wfd1.nFileSizeHigh << "\n";
        std::cout << "files_different: File Size B: " << wfad2.nFileSizeLow << ", " << wfad2.nFileSizeHigh << "\n";
        std::cout << "files_different: Last Write A: " << wfd1.ftLastWriteTime.dwLowDateTime << ", "
                  << wfd1.ftLastWriteTime.dwHighDateTime << "\n";
        std::cout << "files_different: Last Write B: " << wfad2.ftLastWriteTime.dwLowDateTime << ", "
                  << wfad2.ftLastWriteTime.dwHighDateTime << "\n";
    }
    //*/

    return result;
}

bool files::files_different( WIN32_FILE_ATTRIBUTE_DATA& wfad1, WIN32_FILE_ATTRIBUTE_DATA& wfad2 ) {
    bool result = false;

    if( wfad1.nFileSizeLow != wfad2.nFileSizeLow || wfad1.nFileSizeHigh != wfad2.nFileSizeHigh )
        result = true;

    if( wfad1.ftLastWriteTime.dwLowDateTime != wfad2.ftLastWriteTime.dwLowDateTime ||
        wfad1.ftLastWriteTime.dwHighDateTime != wfad2.ftLastWriteTime.dwHighDateTime )
        result = true;

    //*
    if( result ) {
        std::cout << "files_different: File Size A: " << wfad1.nFileSizeLow << ", " << wfad1.nFileSizeHigh << "\n";
        std::cout << "files_different: File Size B: " << wfad2.nFileSizeLow << ", " << wfad2.nFileSizeHigh << "\n";
        std::cout << "files_different: Last Write A: " << wfad1.ftLastWriteTime.dwLowDateTime << ", "
                  << wfad1.ftLastWriteTime.dwHighDateTime << "\n";
        std::cout << "files_different: Last Write B: " << wfad2.ftLastWriteTime.dwLowDateTime << ", "
                  << wfad2.ftLastWriteTime.dwHighDateTime << "\n";
    }
    //*/

    return result;
}

bool files::files_different( WIN32_FIND_DATA& wfd1, const frantic::tstring& fileName ) {
    WIN32_FILE_ATTRIBUTE_DATA fileAD;

    if( !GetFileAttributesEx( fileName.c_str(), GetFileExInfoStandard, &fileAD ) )
        throw std::runtime_error( "files_different: Could not get file attributes of file \"" +
                                  frantic::strings::to_string( fileName ) + "\": " + win32::GetLastErrorMessageA() );

    bool result = files_different( wfd1, fileAD );

    //*
    if( result ) {
        std::cout << "files_different: Files " << frantic::strings::to_string( wfd1.cFileName ) << " and "
                  << frantic::strings::to_string( fileName ) << " differed" << std::endl;
    }
    //*/

    return result;
}

// Returns true if the last modified time and the file sizes match, false otherwise
bool files::files_different( const frantic::tstring& first, const frantic::tstring& second ) {
    WIN32_FILE_ATTRIBUTE_DATA fileAD1, fileAD2;

    if( !GetFileAttributesEx( first.c_str(), GetFileExInfoStandard, &fileAD1 ) )
        throw std::runtime_error( "files_different: Could not get file attributes of file \"" +
                                  frantic::strings::to_string( first ) + "\": " + win32::GetLastErrorMessageA() );

    if( !GetFileAttributesEx( second.c_str(), GetFileExInfoStandard, &fileAD2 ) )
        throw std::runtime_error( "files_different: Could not get file attributes of file \"" +
                                  frantic::strings::to_string( second ) + "\": " + win32::GetLastErrorMessageA() );

    bool result = files_different( fileAD1, fileAD2 );

    //*
    if( result ) {
        std::cout << "files_different: Files " << frantic::strings::to_string( first ) << " and "
                  << frantic::strings::to_string( second ) << " differed" << std::endl;
    }
    //*/

    return result;
}
#else
bool files::files_different( const std::string& /*first*/, const std::string& /*second*/ ) {
    throw std::runtime_error( "files_different not implemented in unix yet." );
}
#endif

bool files::files_different( const fs::path& first, const fs::path& second ) {
    return files_different( frantic::files::to_tstring( first ), frantic::files::to_tstring( second ) );
}

namespace {
template <class CharType>
inline std::basic_string<CharType> ff_first_available_numbered_file( const std::basic_string<CharType>& fileName,
                                                                     int digits ) {
    int currDigits = digits;
    int index = 0;
    std::basic_string<CharType> testFile;
    do {
        // testFile = files::ReplaceSequenceNumber( fileName, index, digits );
        testFile = files::replace_sequence_number( fileName, index, currDigits );
        ++index;

        // increment digit count if index has more digits than digit count
        char temp[1000]; // TODO: replace this with an appropriate MAXPATH value
        sprintf( temp, "%d", index );
        if( (int)strlen( temp ) > currDigits ) {
            ++currDigits;
        }
    } while( files::file_exists( testFile ) );
    return testFile;
}
} // anonymous namespace

std::string files::first_available_numbered_file( const std::string& fileName, int digits ) {
    return ff_first_available_numbered_file( fileName, digits );
}

std::wstring files::first_available_numbered_file( const std::wstring& fileName, int digits ) {
    return ff_first_available_numbered_file( fileName, digits );
}

std::string files::first_available_numbered_file( const std::string& fileName ) {
    return first_available_numbered_file( fileName, 0 );
}

void files::rename_to_backup( const std::string& file, const std::string& backupTag ) {
    std::string taggedFile = files::directory_from_path( file ) + "\\" + files::basename_from_path( file ) + "_" +
                             backupTag + files::extension_from_path( file );
#if defined( _WIN32 ) || defined( _WIN64 )
    MoveFileA( file.c_str(), first_available_numbered_file( taggedFile, 0 ).c_str() );
#else
    rename( file.c_str(), first_available_numbered_file( taggedFile, 0 ).c_str() );
#endif
}

#if defined( _WIN32 ) || defined( _WIN64 )
void files::rename_to_backup( const std::wstring& file, const std::wstring& backupTag ) {
    std::wstring taggedFile = files::directory_from_path( file ) + L"\\" + files::basename_from_path( file ) + L"_" +
                              backupTag + files::extension_from_path( file );
    MoveFileW( file.c_str(), first_available_numbered_file( taggedFile, 0 ).c_str() );
}
#endif

frantic::tstring files::search_file_list( const frantic::tstring& semicolon_file_list ) {
    std::vector<frantic::tstring> separated;
    frantic::strings::split( semicolon_file_list, separated, _T( ';' ) );
    for( unsigned i = 0; i < separated.size(); ++i ) {
        if( files::file_exists( separated[i] ) )
            return separated[i];
    }
    return _T("");
}

#ifdef _WIN32
frantic::tstring files::search_file_list_for_64bit( const frantic::tstring& semicolon_file_list ) {
    std::vector<frantic::tstring> separated;
    frantic::strings::split( semicolon_file_list, separated, _T( ';' ) );
    for( unsigned i = 0; i < separated.size(); ++i ) {
        if( files::file_exists( separated[i] ) ) {
            if( win32::IsFile64BitDLLorEXE( separated[i] ) )
                return separated[i];
        }
    }
    return _T("");
}

frantic::tstring files::search_file_list_for_32bit( const frantic::tstring& semicolon_file_list ) {
    std::vector<frantic::tstring> separated;
    frantic::strings::split( semicolon_file_list, separated, _T( ';' ) );
    for( unsigned i = 0; i < separated.size(); ++i ) {
        if( files::file_exists( separated[i] ) ) {
            if( !win32::IsFile64BitDLLorEXE( separated[i] ) )
                return separated[i];
        }
    }
    return _T("");
}
#endif

frantic::tstring files::search_directory_list_for_directory( const frantic::tstring& semicolon_directory_list ) {
    std::vector<frantic::tstring> separated;
    frantic::strings::split( semicolon_directory_list, separated, _T( ';' ) );
    for( unsigned i = 0; i < separated.size(); ++i ) {
        frantic::tstring directory = separated[i];

        // remove trailing forward slash on directory
        if( directory[directory.size() - 1] == _T( '/' ) || directory[directory.size() - 1] == _T( '\\' ) )
            directory = directory.substr( 0, directory.size() - 1 );

        if( files::directory_exists( directory ) )
            return directory;
    }
    return _T("");
}

frantic::tstring files::find_file_in_search_path( const frantic::tstring& fileName,
                                                  const frantic::tstring& searchPath ) {
    std::vector<frantic::tstring> searchList;
    strings::split_search_path( searchPath, searchList );
    for( unsigned i = 0; i < searchList.size(); ++i ) {
        fs::path possibleMatch = fs::path( files::normalized_directory_name( searchList[i] ) ) / fileName;
        if( fs::exists( possibleMatch ) ) {
            return frantic::files::to_tstring( possibleMatch );
        }
    }
    return _T("");
}

frantic::tstring files::search_directory_list( const frantic::tstring& semicolon_directory_list,
                                               const frantic::tstring& file ) {
    std::vector<frantic::tstring> separated;
    frantic::strings::split( semicolon_directory_list, separated, _T( ';' ) );
    for( unsigned i = 0; i < separated.size(); ++i ) {
        frantic::tstring candidateFile = files::ensure_trailing_pathseparator( separated[i] ) + file;
        if( files::file_exists( candidateFile ) )
            return candidateFile;
    }
    return _T("");
}

frantic::tstring files::search_path( const frantic::tstring& file ) {
    frantic::tstring path = frantic::os::get_environment_variable( _T("PATH") );
    if( path.empty() )
        return _T("");
    return search_directory_list( path, file );
}

void files::add_to_path( const frantic::tstring& newDirectory ) {
    frantic::tstring path = frantic::os::get_environment_variable( _T("PATH") );

    if( path != _T("") )
#if defined( _WIN32 )
        path += _T(";");
#else
        path += _T(":");
#endif

    path += newDirectory;
    frantic::os::set_environment_variable( _T("PATH"), path );
}

void files::get_filenames_in_directory( const frantic::tstring& directory,
                                        std::vector<frantic::tstring>& fileListing ) {
#if defined( _WIN32 ) || defined( _WIN64 )

    // Win32 implementation
    WIN32_FIND_DATA FindFileData;
    HANDLE hFind = INVALID_HANDLE_VALUE;
    DWORD gfid_dwError;

    frantic::tstring strictDir = ensure_trailing_pathseparator( strict_path_name( backward_slashes( directory ) ) ) +
                                 frantic::tstring( _T("*") );

    hFind = FindFirstFile( strictDir.c_str(), &FindFileData );
    if( hFind == INVALID_HANDLE_VALUE ) {
        gfid_dwError = GetLastError();
        std::string errorMessage = std::string( "files::get_filenames_in_directory():  Error code :  " ) +
                                   frantic::win32::GetErrorMessageA( gfid_dwError );
        throw std::runtime_error( errorMessage.c_str() );
    }

    fileListing.clear();

    if( ( FindFileData.dwFileAttributes != FILE_ATTRIBUTE_DIRECTORY ) &&
        ( FindFileData.dwFileAttributes != FILE_ATTRIBUTE_HIDDEN ) )
        fileListing.push_back( frantic::tstring( FindFileData.cFileName ) );

    while( FindNextFile( hFind, &FindFileData ) != 0 )
        if( ( FindFileData.dwFileAttributes != FILE_ATTRIBUTE_DIRECTORY ) &&
            ( FindFileData.dwFileAttributes != FILE_ATTRIBUTE_HIDDEN ) )
            fileListing.push_back( frantic::tstring( FindFileData.cFileName ) );

    gfid_dwError = GetLastError();
    FindClose( hFind );

    if( gfid_dwError != ERROR_NO_MORE_FILES ) {
        std::string errorMessage =
            std::string( "files::get_filenames_in_directory():  Unexpected error in FindNextFile loop, code:  " ) +
            frantic::win32::GetErrorMessageA( gfid_dwError );
        throw std::runtime_error( errorMessage.c_str() );
    }
#else
    fileListing.clear();

    if( !fs::exists( directory ) ) {
        throw std::runtime_error( "files::get_filenames_in_directory(): Directory does not exist: " + directory );
    }

    for( fs::directory_iterator i( directory ), ie; i != ie; ++i ) {
        bool accept = false;

        if( fs::exists( *i ) ) {
            if( fs::is_regular_file( *i ) ) {
                accept = true;
            } else if( fs::is_symlink( *i ) ) {
                // TODO: switch to canonical when we update boost
                // fs::path p( fs::canonical( *i ) );
                fs::path p( fs::read_symlink( *i ) );
                if( fs::exists( p ) && fs::is_regular_file( p ) ) {
                    accept = true;
                }
            }
        }

        if( accept ) {
            fileListing.push_back( i->path().filename().c_str() );
        }
    }
#endif
}