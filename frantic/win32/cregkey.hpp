// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once
/////////////////////////////////////////////////////////////////////////////
// CRegKey

// This is the microsoft code removed from the ATL components so it will run
// in not ATL apps.
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <assert.h>
#include <windows.h>

class CRegKey {
  public:
    CRegKey();
    ~CRegKey();

    // Attributes
  public:
    operator HKEY() const;
    HKEY m_hKey;

    // Operations
  public:
    LONG SetValue( DWORD dwValue, LPCTSTR lpszValueName );
    LONG QueryValue( DWORD& dwValue, LPCTSTR lpszValueName );
    LONG QueryValue( LPTSTR szValue, LPCTSTR lpszValueName, DWORD* pdwCount );
    LONG SetValue( LPCTSTR lpszValue, LPCTSTR lpszValueName = NULL );

    LONG SetKeyValue( LPCTSTR lpszKeyName, LPCTSTR lpszValue, LPCTSTR lpszValueName = NULL );
    static LONG WINAPI SetValue( HKEY hKeyParent, LPCTSTR lpszKeyName, LPCTSTR lpszValue,
                                 LPCTSTR lpszValueName = NULL );

    LONG Create( HKEY hKeyParent, LPCTSTR lpszKeyName, LPTSTR lpszClass = REG_NONE,
                 DWORD dwOptions = REG_OPTION_NON_VOLATILE, REGSAM samDesired = KEY_ALL_ACCESS,
                 LPSECURITY_ATTRIBUTES lpSecAttr = NULL, LPDWORD lpdwDisposition = NULL );
    LONG Open( HKEY hKeyParent, LPCTSTR lpszKeyName, REGSAM samDesired = KEY_ALL_ACCESS );
    LONG Close();
    HKEY Detach();
    void Attach( HKEY hKey );
    LONG DeleteSubKey( LPCTSTR lpszSubKey );
    LONG RecurseDeleteKey( LPCTSTR lpszKey );
    LONG DeleteValue( LPCTSTR lpszValue );
};

inline CRegKey::CRegKey() { m_hKey = NULL; }

inline CRegKey::~CRegKey() { Close(); }

inline CRegKey::operator HKEY() const { return m_hKey; }

inline HKEY CRegKey::Detach() {
    HKEY hKey = m_hKey;
    m_hKey = NULL;
    return hKey;
}

inline void CRegKey::Attach( HKEY hKey ) {
    assert( m_hKey == NULL );
    m_hKey = hKey;
}

inline LONG CRegKey::DeleteSubKey( LPCTSTR lpszSubKey ) {
    assert( m_hKey != NULL );
    return RegDeleteKey( m_hKey, lpszSubKey );
}

inline LONG CRegKey::DeleteValue( LPCTSTR lpszValue ) {
    assert( m_hKey != NULL );
    return RegDeleteValue( m_hKey, (LPTSTR)lpszValue );
}

inline LONG CRegKey::Close() {
    LONG lRes = ERROR_SUCCESS;
    if( m_hKey != NULL ) {
        lRes = RegCloseKey( m_hKey );
        m_hKey = NULL;
    }
    return lRes;
}

inline LONG CRegKey::Create( HKEY hKeyParent, LPCTSTR lpszKeyName, LPTSTR lpszClass, DWORD dwOptions, REGSAM samDesired,
                             LPSECURITY_ATTRIBUTES lpSecAttr, LPDWORD lpdwDisposition ) {
    assert( hKeyParent != NULL );
    DWORD dw;
    HKEY hKey = NULL;
    LONG lRes = RegCreateKeyEx( hKeyParent, lpszKeyName, 0, lpszClass, dwOptions, samDesired, lpSecAttr, &hKey, &dw );
    if( lpdwDisposition != NULL )
        *lpdwDisposition = dw;
    if( lRes == ERROR_SUCCESS ) {
        lRes = Close();
        m_hKey = hKey;
    }
    return lRes;
}

inline LONG CRegKey::Open( HKEY hKeyParent, LPCTSTR lpszKeyName, REGSAM samDesired ) {
    assert( hKeyParent != NULL );
    HKEY hKey = NULL;
    LONG lRes = RegOpenKeyEx( hKeyParent, lpszKeyName, 0, samDesired, &hKey );
    if( lRes == ERROR_SUCCESS ) {
        lRes = Close();
        assert( lRes == ERROR_SUCCESS );
        m_hKey = hKey;
    }
    return lRes;
}

inline LONG CRegKey::QueryValue( DWORD& dwValue, LPCTSTR lpszValueName ) {
    DWORD dwType = NULL;
    DWORD dwCount = sizeof( DWORD );
    LONG lRes = RegQueryValueEx( m_hKey, (LPTSTR)lpszValueName, NULL, &dwType, (LPBYTE)&dwValue, &dwCount );
    assert( ( lRes != ERROR_SUCCESS ) || ( dwType == REG_DWORD ) );
    assert( ( lRes != ERROR_SUCCESS ) || ( dwCount == sizeof( DWORD ) ) );
    return lRes;
}

inline LONG CRegKey::QueryValue( LPTSTR szValue, LPCTSTR lpszValueName, DWORD* pdwCount ) {
    assert( pdwCount != NULL );
    DWORD dwType = NULL;
    LONG lRes = RegQueryValueEx( m_hKey, (LPTSTR)lpszValueName, NULL, &dwType, (LPBYTE)szValue, pdwCount );
    assert( ( lRes != ERROR_SUCCESS ) || ( dwType == REG_SZ ) || ( dwType == REG_MULTI_SZ ) ||
            ( dwType == REG_EXPAND_SZ ) );
    return lRes;
}

inline LONG WINAPI CRegKey::SetValue( HKEY hKeyParent, LPCTSTR lpszKeyName, LPCTSTR lpszValue, LPCTSTR lpszValueName ) {
    assert( lpszValue != NULL );
    CRegKey key;
    LONG lRes = key.Create( hKeyParent, lpszKeyName );
    if( lRes == ERROR_SUCCESS )
        lRes = key.SetValue( lpszValue, lpszValueName );
    return lRes;
}

inline LONG CRegKey::SetKeyValue( LPCTSTR lpszKeyName, LPCTSTR lpszValue, LPCTSTR lpszValueName ) {
    assert( lpszValue != NULL );
    CRegKey key;
    LONG lRes = key.Create( m_hKey, lpszKeyName );
    if( lRes == ERROR_SUCCESS )
        lRes = key.SetValue( lpszValue, lpszValueName );
    return lRes;
}

inline LONG CRegKey::SetValue( DWORD dwValue, LPCTSTR lpszValueName ) {
    assert( m_hKey != NULL );
    return RegSetValueEx( m_hKey, lpszValueName, NULL, REG_DWORD, (BYTE* const)&dwValue, sizeof( DWORD ) );
}

inline HRESULT CRegKey::SetValue( LPCTSTR lpszValue, LPCTSTR lpszValueName ) {
    assert( lpszValue != NULL );
    assert( m_hKey != NULL );
    return RegSetValueEx( m_hKey, lpszValueName, NULL, REG_SZ, (BYTE* const)lpszValue,
                          ( lstrlen( lpszValue ) + 1 ) * sizeof( TCHAR ) );
}

// RecurseDeleteKey is necessary because on NT RegDeleteKey doesn't work if the
// specified key has subkeys
inline LONG CRegKey::RecurseDeleteKey( LPCTSTR lpszKey ) {
    CRegKey key;
    LONG lRes = key.Open( m_hKey, lpszKey );
    if( lRes != ERROR_SUCCESS )
        return lRes;
    FILETIME time;
    TCHAR szBuffer[256];
    DWORD dwSize = 256;
    while( RegEnumKeyEx( key.m_hKey, 0, szBuffer, &dwSize, NULL, NULL, NULL, &time ) == ERROR_SUCCESS ) {
        lRes = key.RecurseDeleteKey( szBuffer );
        if( lRes != ERROR_SUCCESS )
            return lRes;
        dwSize = 256;
    }
    key.Close();
    return DeleteSubKey( lpszKey );
}
