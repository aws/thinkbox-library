// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <commdlg.h>
#include <windows.h>

#include <stdexcept>
#include <string>
#include <vector>

#pragma warning( push )
#pragma warning( disable : 4512 )
#include <frantic/diagnostics/assert_macros.hpp>

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace win32 {

// This helper function will build a dropdown menu and add it to an existing menubar.
// The intended usage is to pass it the window's menu bar (ie. Topmost menu), and it will
// fill out a new drop down menu at end, with the specified buttons.
// For example it will add the File dropwdown with Save, Open, and Close options;
inline HMENU AppendDropDownMenu( HMENU menuBar, const frantic::tstring& name,
                                 const std::vector<std::pair<frantic::tstring, UINT>>& buttons ) {
    HMENU theMenu = CreatePopupMenu();

    MENUITEMINFO mii;
    memset( &mii, 0, sizeof( MENUITEMINFO ) );
    mii.cbSize = sizeof( MENUITEMINFO );
    mii.fMask = MIIM_STRING | MIIM_FTYPE | MIIM_SUBMENU;
    mii.fType = MFT_STRING;
    mii.dwTypeData = (TCHAR*)name.c_str();
    mii.hSubMenu = theMenu;

    InsertMenuItem( menuBar, GetMenuItemCount( menuBar ), TRUE, &mii );

    for( std::size_t i = 0; i < buttons.size(); ++i ) {
        memset( &mii, 0, sizeof( MENUITEMINFO ) );
        mii.cbSize = sizeof( MENUITEMINFO );

        if( buttons[i].first == _T("-") && buttons[i].second == 0 ) {
            mii.fMask = MIIM_FTYPE;
            mii.fType = MFT_SEPARATOR;
        } else {
            mii.fMask = MIIM_STRING | MIIM_FTYPE | MIIM_ID;
            mii.fType = MFT_STRING;
            mii.dwTypeData = (TCHAR*)buttons[i].first.c_str();
            mii.wID = buttons[i].second;
        }

        InsertMenuItem( theMenu, (UINT)i, TRUE, &mii );
    }

    return theMenu;
}

class AutoRegKey {
  private:
    HKEY m_key;

  public:
    static HKEY create_key( const HKEY parent, const frantic::tstring& keyName, bool bForce64Bit = false ) {
        HKEY key;
        DWORD result = RegCreateKeyEx( parent, keyName.c_str(), 0, NULL, REG_OPTION_NON_VOLATILE,
                                       KEY_ALL_ACCESS | ( bForce64Bit ? KEY_WOW64_64KEY : 0 ), NULL, &key, NULL );

        FRANTIC_ASSERT_THROW( ERROR_SUCCESS == result,
                              "Failed to create the key \"" + frantic::strings::to_string( keyName ) + "\"" );

        return key;
    }

    static HKEY open_key( const HKEY parent, const frantic::tstring& keyName, bool bForce64Bit = false ) {
        HKEY key;
        DWORD result =
            RegOpenKeyEx( parent, keyName.c_str(), 0, KEY_ALL_ACCESS | ( bForce64Bit ? KEY_WOW64_64KEY : 0 ), &key );

        FRANTIC_ASSERT_THROW( ERROR_SUCCESS == result,
                              "Failed to open the key \"" + frantic::strings::to_string( keyName ) + "\"" );

        return key;
    }

    AutoRegKey()
        : m_key( (HKEY)-1 ) {}
    explicit AutoRegKey( const HKEY k )
        : m_key( k ) {}
    explicit AutoRegKey( const HKEY parent, const frantic::tstring& keyName, bool bForce64Bit = false )
        : m_key( create_key( parent, keyName, bForce64Bit ) ) {}

    ~AutoRegKey() { close(); }

    void reset( HKEY key ) {
        if( m_key != key )
            close();
        m_key = key;
    }

    void set_value_string( const frantic::tstring& valueName, const frantic::tstring& data ) {
        DWORD result = RegSetValueEx( m_key, valueName.c_str(), 0, REG_SZ, (const BYTE*)data.c_str(),
                                      (DWORD)( ( data.size() + 1 ) * sizeof( TCHAR ) ) );
        FRANTIC_ASSERT_THROW( ERROR_SUCCESS == result, "Could not set the value: name:\"" +
                                                           frantic::strings::to_string( valueName ) + "\", data\"" +
                                                           frantic::strings::to_string( data ) + "\"" );
    }

    std::string get_value_string( const std::string& valueName ) {
        // Call it once to get the size needed.
        DWORD buffSize = 0;
        DWORD result = RegQueryValueExA( m_key, valueName.c_str(), 0, NULL, NULL, &buffSize );
        if( ERROR_SUCCESS == result && buffSize > 0 ) {
            // Set the buffer size and retrieve the data
            std::string buffer( buffSize, '\0' );
            result = RegQueryValueExA( m_key, valueName.c_str(), 0, NULL, (BYTE*)&buffer[0], &buffSize );

            if( ERROR_SUCCESS == result )
                return buffer.substr( 0, buffSize - 1 ); // Strip the trailing NULL
        }

        // throw std::runtime_error("AutoRegKey::get_value_string(): failed to get value");
        return "";
    }

    std::wstring get_value_string( const std::wstring& valueName ) {
        // Call it once to get the size needed.
        DWORD buffSize = 0;
        DWORD result = RegQueryValueExW( m_key, valueName.c_str(), 0, NULL, NULL, &buffSize );
        if( ERROR_SUCCESS == result && buffSize > 0 ) {
            // Set the buffer size and retrieve the data
            // Round up to nearest multiple of wchar.  This may be unnecessary.
            std::vector<TCHAR> buffer( ( buffSize + sizeof( wchar_t ) - 1 ) / sizeof( wchar_t ), 0 );
            buffSize = static_cast<DWORD>( buffer.size() * sizeof( wchar_t ) );
            result = RegQueryValueExW( m_key, valueName.c_str(), 0, NULL, (BYTE*)&buffer[0], &buffSize );

            if( ERROR_SUCCESS == result ) {
                // strip the trailing NULL
                if( buffer.size() > 0 ) {
                    if( buffer[buffer.size() - 1] == 0 ) {
                        buffer.resize( buffer.size() - 1 );
                    }
                    return std::wstring( buffer.begin(), buffer.end() );
                } else {
                    return L"";
                }
            }
        }

        // throw std::runtime_error("AutoRegKey::get_value_string(): failed to get value");
        return L"";
    }

    operator HKEY() const { return m_key; }
    void close() {
        if( m_key != (HKEY)-1 ) {
            RegCloseKey( m_key );
            m_key = (HKEY)-1;
        }
    }
};

class AutoHandle {
  private:
    HANDLE m_handle;

  public:
    explicit AutoHandle( HANDLE h ) { m_handle = h; }

    ~AutoHandle() { CloseHandle( m_handle ); }

    bool operator==( const HANDLE& handle ) const { return m_handle == handle; }

    operator HANDLE() { return m_handle; }
    operator const HANDLE() const { return m_handle; }
};

inline frantic::tstring GetModuleFileName( HMODULE hModule = NULL ) {
    TCHAR lpDLLFileName[MAX_PATH];
    if( ::GetModuleFileName( hModule, lpDLLFileName, MAX_PATH ) == 0 )
        return _T("");
    else
        return lpDLLFileName;
}

inline frantic::tstring GetWindowsDirectory() {
    TCHAR lpWindowsDirectory[MAX_PATH] = { 0 };
    if( ::GetWindowsDirectory( lpWindowsDirectory, MAX_PATH ) == 0 )
        return _T("");
    else
        return lpWindowsDirectory;
}

frantic::tstring GetINISetting( const frantic::tstring& section, const frantic::tstring& key,
                                const frantic::tstring& defaultValue, const frantic::tstring& iniFile );
bool SetINISetting( const frantic::tstring& section, const frantic::tstring& key, const frantic::tstring& value,
                    const frantic::tstring& iniFile );
bool RemoveINISetting( const frantic::tstring& section, const frantic::tstring& key, const frantic::tstring& iniFile );

inline frantic::tstring ffGetRegKeyString( HKEY parentKey, const frantic::tstring& subKey,
                                           const frantic::tstring& valueName, bool b64BitRegistry = false ) {
    DWORD lRet;
    HKEY tempKey;
    frantic::tstring answer;

    lRet = RegOpenKeyEx( parentKey, subKey.c_str(), 0,
                         KEY_QUERY_VALUE | ( b64BitRegistry ? KEY_WOW64_64KEY : KEY_WOW64_32KEY ), &tempKey );
    if( ERROR_SUCCESS == lRet ) {
        DWORD buffSize = 0;
        lRet = RegQueryValueEx( tempKey, valueName.c_str(), 0, NULL, NULL, &buffSize );
        if( lRet == ERROR_SUCCESS ) {
            std::vector<TCHAR> buffer( buffSize + 1, 0 );
            lRet = RegQueryValueEx( tempKey, valueName.c_str(), 0, NULL, (BYTE*)&buffer[0], &buffSize );

            if( ERROR_SUCCESS == lRet )
                answer = &buffer[0];
        }

        RegCloseKey( tempKey );
    }

    return answer;
}

inline void ffSetRegKeyString( HKEY hkey, const frantic::tstring& subKey, const frantic::tstring& valueName,
                               const frantic::tstring& value, bool b64BitRegistry = false ) {
    long lRet;
    HKEY hKey;
    lRet = RegOpenKeyEx( hkey, subKey.c_str(), 0,
                         KEY_SET_VALUE | ( b64BitRegistry ? KEY_WOW64_64KEY : KEY_WOW64_32KEY ), &hKey );
    if( lRet != ERROR_SUCCESS )
        throw std::runtime_error( "Could not open registry key " + frantic::strings::to_string( subKey ) +
                                  ", value name " + frantic::strings::to_string( valueName ) + ", to " +
                                  frantic::strings::to_string( value ) );

    lRet =
        RegSetValueEx( hKey, valueName.c_str(), NULL, REG_SZ, (const BYTE*)value.c_str(), (DWORD)( value.size() + 1 ) );

    if( lRet != ERROR_SUCCESS )
        throw std::runtime_error( "Could not set registry key " + frantic::strings::to_string( subKey ) +
                                  ", value name " + frantic::strings::to_string( valueName ) + ", to " +
                                  frantic::strings::to_string( value ) );
    RegCloseKey( hKey );
}

/**
 *  Get the subkey names of a specified key.
 *
 * @todo how to report failure?
 *
 * @param hKey the parent key whose subkeys should be found.
 * @param out[out] the name of immediate descendant subkeys of hKey.  On error,
 *		this may be empty or missing some subkeys.
 */
void GetRegSubKeyNames( HKEY hKey, std::vector<frantic::tstring>& out );

inline std::string GetTempPathA() {
    char lpTempPath[MAX_PATH];
    if( ::GetTempPathA( sizeof( lpTempPath ), lpTempPath ) == 0 )
        return "";
    else
        return lpTempPath;
}

// conditional because GetTempPath is a macro that resolves to GetTempPathA
// which is already defined when we're not using unicode
#if defined( UNICODE ) || defined( _UNICODE )
inline frantic::tstring GetTempPath() {
    // Get the length of the path
    DWORD len = ::GetTempPath( 0, _T("") );
    if( len == 0 )
        return _T("");

    // Temporary buffer, the `len` returned includes the NUL terminator
    std::vector<TCHAR> tempPath( len );
    len = ::GetTempPath( (DWORD)tempPath.size(), &tempPath[0] );
    if( len == 0 || len > tempPath.size() )
        return _T("");

    return frantic::tstring( &tempPath[0], &tempPath[0] + len );
}
#endif

inline std::string ffGetTempFilenameA() {
    char lpTempFilename[MAX_PATH];
    if( ::GetTempFileNameA( frantic::win32::GetTempPathA().c_str(), "tmp", 0, lpTempFilename ) == 0 )
        return "";
    else
        return lpTempFilename;
}

inline std::basic_string<TCHAR> ffGetTempFilename() {
    TCHAR lpTempFilename[MAX_PATH];
    if( ::GetTempFileName( GetTempPath().c_str(), _T("tmp"), 0, lpTempFilename ) == 0 )
        return _T("");
    else
        return lpTempFilename;
}

inline int ffMessageBoxA( HWND hWnd, const std::string& lpText, const std::string& lpCaption, ::UINT uType = MB_OK ) {
    return ::MessageBoxA( hWnd, lpText.c_str(), lpCaption.c_str(), uType );
}

inline int ffMessageBox( HWND hWnd, const frantic::tstring& lpText, const frantic::tstring& lpCaption,
                         ::UINT uType = MB_OK ) {
    return ::MessageBox( hWnd, lpText.c_str(), lpCaption.c_str(), uType );
}

inline frantic::tstring
ffLoadDialog( HWND hwnd, const frantic::tstring& title, const frantic::tstring& initialFile,
              const std::vector<std::pair<frantic::tstring, frantic::tstring>>& filterStrings ) {
    OPENFILENAME ofn;
    memset( &ofn, 0, sizeof( OPENFILENAMEA ) );

    // create a buffer for windows to put the filename in
    TCHAR buffer[MAX_PATH] = { 0 };
    if( initialFile.size() < MAX_PATH )
        memcpy( buffer, initialFile.c_str(), initialFile.size() );

    // create a windows-style filter string
    frantic::tstring filter = _T("");
    for( std::size_t i = 0; i < filterStrings.size(); ++i )
        filter += filterStrings[i].first + TCHAR( 0 ) + filterStrings[i].second + TCHAR( 0 );

    ofn.lStructSize = sizeof( OPENFILENAME );
    ofn.hwndOwner = hwnd;
    ofn.lpstrFile = buffer;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrFilter = filter.c_str();
    ofn.nMaxCustFilter = (DWORD)filter.size() + 1;
    ofn.nFilterIndex = 1;
    ofn.lpstrTitle = title.c_str();
    ofn.Flags = OFN_FILEMUSTEXIST;

    if( GetOpenFileName( &ofn ) )
        return buffer;
    return _T("");
}

inline frantic::tstring
ffSaveDialog( HWND hwnd, const frantic::tstring& title, const frantic::tstring& initialFile,
              const frantic::tstring defaultExtension,
              const std::vector<std::pair<frantic::tstring, frantic::tstring>>& filterStrings ) {
    OPENFILENAME ofn;
    memset( &ofn, 0, sizeof( OPENFILENAME ) );

    // create a buffer for windows to put the filename in
    TCHAR buffer[MAX_PATH] = { 0 };
    if( initialFile.size() < MAX_PATH )
        memcpy( buffer, initialFile.c_str(), initialFile.size() );

    // create a windows-style filter string
    frantic::tstring filter = _T("");
    for( std::size_t i = 0; i < filterStrings.size(); ++i )
        filter += filterStrings[i].first + TCHAR( 0 ) + filterStrings[i].second + TCHAR( 0 );

    ofn.lStructSize = sizeof( OPENFILENAME );
    ofn.hwndOwner = hwnd;
    ofn.lpstrFile = buffer;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrFilter = filter.c_str();
    ofn.nMaxCustFilter = (DWORD)filter.size() + 1;
    ofn.nFilterIndex = 1;
    ofn.lpstrTitle = title.c_str();
    ofn.lpstrDefExt = defaultExtension.c_str();

    if( GetSaveFileName( &ofn ) )
        return buffer;
    return _T("");
}

inline HMODULE LoadLibrary( const frantic::tstring& lpLibFileName ) { return ::LoadLibrary( lpLibFileName.c_str() ); }

#pragma warning( pop )

} // namespace win32
} // namespace frantic
