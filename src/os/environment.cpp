// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/os/environment.hpp>

#if defined( _WIN32 )
#include <frantic/win32/utility.hpp>
#include <windows.h>

// This requires Userenv.lib to be linked to compile
// Used by current_os_environment_variable_reader
#include <UserEnv.h>
#pragma comment( lib, "Userenv.lib" )

// Only used in Windows build right now
#include <frantic/logging/logging_level.hpp>
#endif

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif

namespace frantic {
namespace os {

#if defined( _WIN32 )
namespace {

/**
 * Wrapper class for reading from the set of environment variables from the OS/User rather than the process.
 */
class current_os_environment_variable_reader {
    // Raw list of variables
    LPVOID m_envVariables;

  public:
    current_os_environment_variable_reader()
        : m_envVariables( NULL ) {
        // Note that according to https://msdn.microsoft.com/en-us/library/windows/desktop/ms683179%28v=vs.85%29.aspx
        // the handle returned from GetCurrentProcess does not need to be closed
        HANDLE userToken = NULL;
        const BOOL useUserToken = OpenProcessToken( GetCurrentProcess(), TOKEN_QUERY, &userToken );
        if( !useUserToken ) {
            FF_LOG( warning ) << _T( "current_os_environment_variable_reader: " )
                                 _T( "Unable to get user token to read environment variables from: " )
                              << win32::GetLastErrorMessage() << std::endl;
            userToken = NULL;
        }

        const BOOL useEnvVariables = CreateEnvironmentBlock( &m_envVariables, userToken, FALSE );
        if( !useEnvVariables ) {
            FF_LOG( warning ) << _T( "current_os_environment_variable_reader: " )
                                 _T( "Unable to create environment variable block for reading: " )
                              << win32::GetLastErrorMessage() << std::endl;
            m_envVariables = NULL;
        }

        if( useUserToken ) {
            const BOOL ok = CloseHandle( userToken );
            if( !ok ) {
                FF_LOG( warning ) << _T( "current_os_environment_variable_reader: " )
                                     _T( "Unable to close user token: " )
                                  << win32::GetLastErrorMessage() << std::endl;
            }
        }
    }

    ~current_os_environment_variable_reader() {
        if( m_envVariables ) {
            const BOOL ok = DestroyEnvironmentBlock( m_envVariables );
            if( !ok ) {
                FF_LOG( warning ) << _T( "current_os_environment_variable_reader: " )
                                     _T( "Unable to destroy environment variable block: " )
                                  << win32::GetLastErrorMessage() << std::endl;
            }
        }
    }

    tstring read_value( const tstring& key ) const {
        // According to https://msdn.microsoft.com/en-us/library/windows/desktop/bb762270%28v=vs.85%29.aspx,
        // the environemnt variable list is of type unicode which are wchars
        WCHAR* envVariableList = static_cast<WCHAR*>( m_envVariables );
        WCHAR* envVariableListIter = envVariableList;
        const std::wstring wKey = strings::to_wstring( key );

        // According to https://msdn.microsoft.com/en-us/library/windows/desktop/ms682429%28v=vs.85%29.aspx,
        // the environment list is assumed to be of the format
        //
        // <name1>=<value1>\0
        // <name2>=<value2>\0
        // ...
        // \0
        //
        // Each entry is delimited by a null character with the last entry marked by two null characters
        while( envVariableListIter[0] ) {
            const std::wstring entry = envVariableListIter;

            // Extract the key
            const std::size_t indexOfEqual = entry.find( L'=' );
            if( entry.substr( 0, indexOfEqual ) == wKey ) {
                // Extract the value
                return strings::to_tstring( entry.substr( indexOfEqual + 1 ) );
            }

            // Move to next entry
            envVariableListIter += entry.size() + 1;
        }

        return _T( "" );
    }

    bool is_valid() const { return m_envVariables != NULL; }
};

} // anonymous namespace
#endif

void set_environment_variable( const frantic::tstring& key, const frantic::tstring& value ) {
#if defined( _WIN32 )
    int result = _tputenv( ( key + _T("=") + value ).c_str() );
    if( result < 0 ) {
        std::stringstream ss;
        ss << "Failed to set environment variable \"" + frantic::strings::to_string( key ) + "\" because:\n\n";
        ss << "Error number " << errno << "\n";
        ss << strerror( errno );
        throw std::runtime_error( ss.str() );
    }
#else
    int nResult = setenv( key.c_str(), value.c_str(), 1 );
    if( 0 != nResult ) {
        std::stringstream ss;
        ss << "Failed to set environment variable \"" + key + "\" because:\n\n";
        ss << "Error number " << errno << "\n";
        ss << strerror( errno );
        throw std::runtime_error( ss.str() );
    }
#endif
}

frantic::tstring get_environment_variable( const frantic::tstring& key ) {
#if defined( _WIN32 )
    DWORD buffSize = GetEnvironmentVariable( key.c_str(), 0, 0 );

    if( buffSize == 0 )
        return _T("");

    std::vector<frantic::strings::tchar> buffer( buffSize );
    GetEnvironmentVariable( key.c_str(), &buffer[0], buffSize );
    return frantic::strings::to_tstring( &buffer[0] );
#else
    frantic::tchar* cvalue = getenv( key.c_str() );
    frantic::tstring value;
    if( cvalue != 0 )
        value = cvalue;
    return value;
#endif
}

frantic::tstring get_environment_variable_from_os_scope( const frantic::tstring& key, bool* outIsFromOS ) {
#if defined( _WIN32 )
    current_os_environment_variable_reader envReader;
    if( envReader.is_valid() ) {
        if( outIsFromOS ) {
            *outIsFromOS = true;
        }
        return envReader.read_value( key );
    }
#endif
    // TODO: implement for Linux and Mac
    if( outIsFromOS ) {
        *outIsFromOS = false;
    }
    return get_environment_variable( key );
}

} // namespace os
} // namespace frantic

#if defined( _MSC_VER )
#pragma warning( pop )
#endif
