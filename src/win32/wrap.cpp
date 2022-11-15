// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( _WIN32 )

#include <frantic/win32/wrap.hpp>

namespace frantic {
namespace win32 {

void GetRegSubKeyNames( HKEY hKey, std::vector<frantic::tstring>& out ) {
    out.clear();

    LONG result;

    // set by RegQueryInfoKey()
    DWORD maxSubKeyLen = 0; // in unicode not including NULL terminator
    DWORD subKeyCount = 0;

    result = RegQueryInfoKey( hKey, 0, 0, 0, &subKeyCount, &maxSubKeyLen, 0, 0, 0, 0, 0, 0 );
    if( result != ERROR_SUCCESS ) {
        // failure
        // return value is a system error code
        // throw std::runtime_error( "GetRegSubKeyNames: error getting key info" );
        return;
    }

    std::vector<TCHAR> subKeyName( maxSubKeyLen + 2 );

    for( DWORD i = 0; i < subKeyCount; ++i ) {
        DWORD subKeyNameSize = maxSubKeyLen + 1;
        memset( &subKeyName[0], 0, subKeyName.size() * sizeof( TCHAR ) );
        result = RegEnumKeyEx( hKey, i, &subKeyName[0], &subKeyNameSize, NULL, NULL, NULL, NULL );

        if( result == ERROR_SUCCESS ) {
            out.push_back( frantic::tstring( subKeyName.begin(), subKeyName.end() ) );
        } else if( result == ERROR_MORE_DATA ) {
            break;
        } else if( result == ERROR_NO_MORE_ITEMS ) {
            break;
        } else {
            // failure
            // return value is a system error code
            // throw std::runtime_error( "GetRegSubKeyNames: error getting key name" );
            break;
        }
    }
}

frantic::tstring GetINISetting( const frantic::tstring& section, const frantic::tstring& key,
                                const frantic::tstring& defaultValue, const frantic::tstring& iniFile ) {
    const DWORD bufferSize = 1024;
    TCHAR answer[bufferSize];
    DWORD valueCount = ::GetPrivateProfileString( section.c_str(), key.c_str(), defaultValue.c_str(), answer,
                                                  bufferSize, iniFile.c_str() );

    return frantic::tstring( answer, valueCount );
}

bool SetINISetting( const frantic::tstring& section, const frantic::tstring& key, const frantic::tstring& value,
                    const frantic::tstring& iniFile ) {
    return ::WritePrivateProfileString( section.c_str(), key.c_str(), value.c_str(), iniFile.c_str() ) != 0;
}

bool RemoveINISetting( const frantic::tstring& section, const frantic::tstring& key, const frantic::tstring& iniFile ) {
    return ::WritePrivateProfileString( section.c_str(), key.c_str(), NULL, iniFile.c_str() ) != 0;
}

} // namespace win32
} // namespace frantic

#endif // defined( _WIN32 )
