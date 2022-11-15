// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once
// Header file to provide small utility functions for win32


// Moved stuff which needs winsock2 into simple_socket.hpp
//#include "Winsock2.h"
//#include "Winsock.h"

// Only include this in visual studio 2003 and higher
#if _MSC_VER >= 1300
#include <psapi.h>
#endif

#include <fstream>
#include <stdexcept>
#include <vector>

#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>

#pragma warning( push )
#pragma warning( disable : 4511 4512 )
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#pragma warning( pop )

#include <boost/algorithm/string/predicate.hpp>

#include <frantic/misc/string_functions.hpp>
#include <frantic/strings/tstring.hpp>
#include <frantic/win32/cregkey.hpp>
#include <frantic/win32/wrap.hpp>

namespace frantic {
namespace win32 {

// Disable some VC2005 warnings.
#pragma warning( push )
#pragma warning( disable : 4996 )

inline frantic::tstring GetCurrentUserName() {
    DWORD size = 2048;
    TCHAR buffer[2048] = { 0 };
    if( GetUserName( (LPTSTR)buffer, &size ) ) {
        return buffer;
    }
    return _T("");
}

/**
 * Class used to obtain the systems UUID. Equivalent of doing "wmic path win32_computersystemproduct get uuid"
 */
class get_win32_computer_system_uuid {
  public:
    /**
     * Gets the systems UUID and throws an exception if unable to retrieve the UUID.
     * Virtual such that it can be overwritten while integration testing
     * @return systems UUID in a wide string
     */
    virtual frantic::tstring get() const;
};

// Returns the FormatMessage() result from the error code provided.
inline frantic::tstring GetErrorMessage( DWORD errorCode ) {
    LPVOID lpMsgBuf = 0;
    FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL,
                   errorCode, MAKELANGID( LANG_NEUTRAL, SUBLANG_DEFAULT ), // Default language
                   (LPTSTR)&lpMsgBuf, 0, NULL );
    if( lpMsgBuf == 0 )
        return _T("");
    else {
        frantic::tstring result = (const TCHAR*)lpMsgBuf;
        while( result.size() > 0 &&
               ( result[result.size() - 1] == TCHAR( '\n' ) || result[result.size() - 1] == TCHAR( '\r' ) ) )
            result.resize( result.size() - 1 );
        return result;
    }
}

inline std::string GetErrorMessageA( DWORD errorCode ) {
    return frantic::strings::to_string( GetErrorMessage( errorCode ) );
}

// Returns the FormatMessage() result from the GetLastError() code
inline frantic::tstring GetLastErrorMessage() { return GetErrorMessage( GetLastError() ); }

inline std::string GetLastErrorMessageA() { return frantic::strings::to_string( GetLastErrorMessage() ); }

inline bool RegKeyExists( HKEY hKey, const frantic::tstring& keyToSearch, REGSAM accessRights ) {
    CRegKey regKey;
    bool keyExists = true;
    LONG lres = regKey.Open( hKey, keyToSearch.c_str(), accessRights );
    if( lres != ERROR_SUCCESS ) {
        keyExists = false;
    }
    return keyExists;
}

// Waits for a process handle to terminate, returns true if it terminated, false if it timed out
// timeOut of 0 makes it wait forever
inline bool WaitForProcessToTerminate( HANDLE hProcess, DWORD timeOut ) {
    DWORD endTime = GetTickCount() + timeOut;
    long timeLeft = timeOut;
    DWORD returnCode = STILL_ACTIVE;

    GetExitCodeProcess( hProcess, &returnCode );
    timeLeft = long( endTime - GetTickCount() );

    while( returnCode == STILL_ACTIVE && ( timeOut == 0 || timeLeft > 0 ) ) {
        // Block for a while
        WaitForSingleObject( hProcess, timeLeft );
        // Check if it really exited
        GetExitCodeProcess( hProcess, &returnCode );
        timeLeft = long( endTime - GetTickCount() );
    }

    return returnCode != STILL_ACTIVE;
}

inline DWORD GetExitCodeProcess( HANDLE hProcess ) {
    DWORD exitCode = 0;

    if( !::GetExitCodeProcess( hProcess, &exitCode ) )
        throw std::runtime_error( "GetExitCodeProcess: " + GetLastErrorMessageA() );

    return exitCode;
}

// Only do this in visual studio 2003 and higher
#if _MSC_VER >= 1300

#pragma comment( lib, "psapi.lib" )
// Returns the process id of a given executable, or 0
inline DWORD FindProcess( const frantic::tstring& executable ) {
#define MAX_PROCESSES 300
#define MAX_MODULES 300
    HMODULE phModules[MAX_MODULES];
    HMODULE* phModule;
    HMODULE hModule;
    HANDLE hProcess;
    DWORD lpidProcesses[MAX_PROCESSES];
    DWORD* pdw;
    DWORD dwRet;
    DWORD cbNeeded, cbModulesNeeded;
    BOOL bRet, bFoundDFRNode;
    TCHAR szFilename[MAX_PATH];
    int numProcesses, numModules;

    //	std::ofstream fout("c:/proc.txt");

    // see if it's already running first
    bFoundDFRNode = FALSE;
    bRet = EnumProcesses( lpidProcesses, sizeof( lpidProcesses ), &cbNeeded );
    if( !bRet )
        return 0;

    // for every process found, query it for module info
    numProcesses = cbNeeded / sizeof( DWORD );
    for( int iProcessNum = 0; iProcessNum < numProcesses; iProcessNum++ ) {
        pdw = lpidProcesses;
        pdw += iProcessNum;
        hProcess = OpenProcess( PROCESS_ALL_ACCESS, FALSE, *pdw );
        if( hProcess ) {
            bRet = EnumProcessModules( hProcess, phModules, sizeof( phModules ), &cbModulesNeeded );
            if( bRet ) {
                numModules = cbModulesNeeded / sizeof( HMODULE );
                for( int iModuleNum = 0; iModuleNum < numModules; iModuleNum++ ) {
                    phModule = phModules;
                    phModule += iModuleNum;
                    hModule = (HMODULE)*phModule;

                    dwRet = GetModuleBaseName( hProcess, hModule, szFilename, sizeof( szFilename ) / sizeof( TCHAR ) );
                    if( dwRet ) {
                        //						fout << szFilename << std::endl;
                        if( ::frantic::strings::to_lower( szFilename ) == frantic::strings::to_lower( executable ) ) {
                            CloseHandle( hProcess );
                            //							fout << "Returning pid " << *pdw <<
                            // std::endl;
                            return *pdw;
                        }
                    }
                }
            }
            CloseHandle( hProcess );
        }
    }

    return 0;
}

#endif // _MSC_VER >= 1300

namespace detail {
typedef struct _PROCESS_BASIC_INFORMATION {
    PVOID Reserved1;
    /*PPEB*/ void* PebBaseAddress;
    PVOID Reserved2[2];
    unsigned long UniqueProcessId;
    PVOID Reserved3;
} PROCESS_BASIC_INFORMATION;

typedef long NTSTATUS;

typedef enum _PROCESSINFOCLASS {
    ProcessBasicInformation,
    ProcessQuotaLimits,
    ProcessIoCounters,
    ProcessVmCounters,
    ProcessTimes,
    ProcessBasePriority,
    ProcessRaisePriority,
    ProcessDebugPort,
    ProcessExceptionPort,
    ProcessAccessToken,
    ProcessLdtInformation,
    ProcessLdtSize,
    ProcessDefaultHardErrorMode,
    ProcessIoPortHandlers, // Note: this is kernel mode only
    ProcessPooledUsageAndLimits,
    ProcessWorkingSetWatch,
    ProcessUserModeIOPL,
    ProcessEnableAlignmentFaultFixup,
    ProcessPriorityClass,
    ProcessWx86Information,
    ProcessHandleCount,
    ProcessAffinityMask,
    ProcessPriorityBoost,
    MaxProcessInfoClass
} PROCESSINFOCLASS;
} // namespace detail

// TODO: move this somewhere else
inline std::string getTimeZoneOffset() {
    _tzset();
    int hours = -( _timezone / 3600 );
    int mins = _timezone + ( hours * 3600 );

    char buf[16];
    sprintf( buf, "%+03d%02d", hours, mins );
    return buf;
}

inline std::string GetLocalIsoTime() {
    SYSTEMTIME t;
    GetLocalTime( &t );
    char result[64];
    sprintf( result, "%d-%02d-%02d %02d:%02d:%02d%s", t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond,
             getTimeZoneOffset().c_str() );
    return result;
    //	return (boost::format("%d-%02d-%02d %02d:%02d:%02d") % t.wYear % t.wMonth % t.wDay % t.wHour % t.wMinute %
    // t.wSecond).str();
}

inline std::string GetLocalCleanedIsoTime() {
    std::string result = GetLocalIsoTime();
    for( unsigned i = 0; i < result.size(); ++i ) {
        if( result[i] == ' ' || result[i] == ':' )
            result[i] = '_';
    }
    return result;
}

// TODO: move this somewhere else
inline std::string GetLocalIsoDate() {
    SYSTEMTIME t;
    GetLocalTime( &t );
    char result[64];
    sprintf( result, "%d-%02d-%02d", t.wYear, t.wMonth, t.wDay );
    return result;
    //	return (boost::format("%d-%02d-%02d %02d:%02d:%02d") % t.wYear % t.wMonth % t.wDay % t.wHour % t.wMinute %
    // t.wSecond).str();
}

namespace GetChildWindows_detail {
BOOL CALLBACK EnumChildProc( HWND hwnd, LPARAM lParam );
}

inline void GetChildWindows( HWND parentWindow, std::vector<HWND>& children ) {
    if( !EnumChildWindows( parentWindow, GetChildWindows_detail::EnumChildProc, (LPARAM)&children ) ) {
        // Sometimes it says it succeeded, but returns false.
        if( GetLastError() != ERROR_SUCCESS )
            throw std::runtime_error( "GetChildWindows: Error enumerating child windows: " + GetLastErrorMessageA() );
    }
}

inline frantic::tstring ffGetWindowText( HWND hWnd ) {
    int textLength = GetWindowTextLength( hWnd );
    if( textLength <= 0 )
        return _T("");
    else {
        boost::shared_array<TCHAR> text( new TCHAR[textLength + 1] );
        if( GetWindowText( hWnd, text.get(), textLength + 1 ) > 0 )
            return text.get();
        else
            return _T("");
    }
}

inline frantic::tstring ffGetWindowClassName( HWND hWnd ) {
    int textLength = 1024;
    if( textLength == 0 )
        return _T("");
    else {
        boost::shared_array<TCHAR> text( new TCHAR[textLength + 1] );
        if( GetClassName( hWnd, text.get(), textLength + 1 ) > 0 )
            return text.get();
        else
            return _T("");
    }
}

inline int ffGetListBoxItemCount( HWND hWnd ) { return (int)SendMessage( hWnd, LB_GETCOUNT, 0, 0 ); }

inline frantic::tstring ffGetListBoxItemText( HWND hWnd, int itemIndex ) {
    int textLength = (int)SendMessage( hWnd, LB_GETTEXTLEN, itemIndex, 0 );
    if( textLength <= 0 )
        return _T("");
    else {
        boost::shared_array<TCHAR> text( new TCHAR[textLength + 1] );
        if( SendMessage( hWnd, LB_GETTEXT, itemIndex, (LPARAM)text.get() ) > 0 )
            return text.get();
        else
            return _T("");
    }
}

inline frantic::tstring ffGetPopupSummary( HWND hWnd ) {
    std::vector<HWND> childHwnds;
    try {
        win32::GetChildWindows( hWnd, childHwnds );
    } catch( const std::exception& ) {
        // TODO: what's with the eating the exception here???
    }
    frantic::tstring mostLikelyMessage = _T("");
    for( size_t j = 0; j < childHwnds.size(); ++j ) {
        frantic::tstring childText = win32::ffGetWindowText( childHwnds[j] );
        // Only consider non-empty visible child windows
        if( ( GetWindowLong( childHwnds[j], GWL_STYLE ) & WS_VISIBLE ) ) {
            // Consider the longest message to be the likely error text
            if( childText.size() > mostLikelyMessage.size() )
                mostLikelyMessage = childText;
        }
    }
    return _T("Title \"") + win32::ffGetWindowText( hWnd ) + _T("\", Message \"") + mostLikelyMessage + _T("\"");
}

inline int ffListBox_GetSelCount( HWND hwnd ) {
    int selCount = (int)SendMessage( hwnd, LB_GETSELCOUNT, 0, 0 );
    if( selCount == LB_ERR ) {
        return 0; // what to do ?
    } else {
        return selCount;
    }
}

inline void ffListBox_GetSelItems( HWND hwnd, std::vector<int>& outIndices ) {
    outIndices.clear();

    int selCount = (int)SendMessage( hwnd, LB_GETSELCOUNT, 0, 0 );
    if( selCount != LB_ERR && selCount > 0 ) {
        std::vector<int> temp( selCount );
        selCount = (int)SendMessage( hwnd, LB_GETSELITEMS, selCount, LPARAM( &temp[0] ) );
        if( selCount != LB_ERR ) {
            temp.resize( selCount );
            temp.swap( outIndices );
        }
    }
}

inline void ffListBox_SetSelItems( HWND hwnd, const std::vector<int>& indices ) {
    SendMessage( hwnd, LB_SETSEL, FALSE, -1 );

    const LRESULT listBoxItemCount = SendMessage( hwnd, LB_GETCOUNT, 0, 0 );
    if( listBoxItemCount != LB_ERR ) {
        for( std::vector<int>::const_iterator i = indices.begin(); i != indices.end(); ++i ) {
            const int index = *i;
            if( index >= 0 && (LRESULT)index < listBoxItemCount ) {
                SendMessage( hwnd, LB_SETSEL, TRUE, (LPARAM)index );
            }
        }
    }
}

inline void ffListBox_Clear( HWND hwnd ) { SendMessage( hwnd, LB_RESETCONTENT, 0, 0 ); }

inline void ffListBox_AddStrings( HWND hwnd, const std::vector<frantic::tstring>& strings ) {
    for( std::size_t i = 0; i < strings.size(); ++i ) {
        SendMessage( hwnd, LB_INSERTSTRING, WPARAM( -1 ), LPARAM( strings[i].c_str() ) );
    }
}

inline void ffListBox_SetStrings( HWND hwnd, const std::vector<frantic::tstring>& strings ) {
    ffListBox_Clear( hwnd );

    std::size_t prealloc = 0;
    for( std::size_t i = 0; i < strings.size(); ++i ) {
        prealloc += 1 + strings[i].size();
    }

    SendMessage( hwnd, LB_INITSTORAGE, strings.size() * sizeof( TCHAR ), prealloc );

    ffListBox_AddStrings( hwnd, strings );
}

inline frantic::tstring ffListBox_GetItemString( HWND hwnd, int index ) {
    LRESULT len = SendMessage( hwnd, LB_GETTEXTLEN, index, 0 );
    if( len == LB_ERR ) {
        return _T("");
    }
    std::vector<TCHAR> buffer( 1 + len, 0 );
    len = SendMessage( hwnd, LB_GETTEXT, index, LPARAM( &buffer[0] ) );
    if( len == LB_ERR ) {
        return _T("");
    }
    buffer.resize( 1 + len );
    return frantic::tstring( &buffer[0] );
}

inline void ffComboBox_Clear( HWND hwnd ) { SendMessage( hwnd, CB_RESETCONTENT, 0, 0 ); }

inline int ffComboBox_AddString( HWND hwnd, const frantic::tstring& s ) {
    return (int)SendMessage( hwnd, CB_ADDSTRING, 0L, (LPARAM)s.c_str() );
}

inline void ffComboBox_SetStrings( HWND hwnd, const std::vector<frantic::tstring>& strings ) {
    ffComboBox_Clear( hwnd );

    std::size_t prealloc = 0;

    for( std::size_t i = 0; i < strings.size(); ++i ) {
        prealloc += 1 + strings[i].size();
    }

    SendMessage( hwnd, CB_INITSTORAGE, strings.size() * sizeof( TCHAR ), prealloc );

    for( std::size_t i = 0; i < strings.size(); ++i ) {
        SendMessage( hwnd, CB_ADDSTRING, 0, LPARAM( strings[i].c_str() ) );
    }
}

inline bool ffComboBox_GetItemString( HWND hwnd, int i, frantic::tstring& out ) {
    LRESULT expectLen = SendMessage( hwnd, CB_GETLBTEXTLEN, i, 0 );
    if( expectLen != CB_ERR && expectLen >= 0 ) {
        std::vector<TCHAR> buffer( expectLen + 2 );
        LRESULT gotLen = SendMessage( hwnd, CB_GETLBTEXT, i, LPARAM( &buffer[0] ) );
        if( gotLen != CB_ERR ) {
            if( gotLen > expectLen ) {
                throw std::runtime_error( "ComboBox_GetItemString Error: string length exceeds expected length" );
            }
            buffer.resize( gotLen );
            out = frantic::tstring( buffer.begin(), buffer.end() );
            return true;
        }
    }
    return false;
}

inline frantic::tstring ffComboBox_GetItemString( HWND hwnd, int i ) {
    frantic::tstring out;
    bool success = ffComboBox_GetItemString( hwnd, i, out );
    if( success ) {
        return out;
    } else {
        return _T("");
    }
}

// case insensitive search
// returns -1 if the item is not found
inline int ffComboBox_IFindItem( HWND hwnd, const frantic::tstring& s ) {
    frantic::tstring itemString;
    const int itemCount = static_cast<int>( SendMessage( hwnd, CB_GETCOUNT, 0, 0 ) );
    if( itemCount == CB_ERR ) {
        return -1;
    }
    for( int i = 0; i < itemCount; ++i ) {
        bool gotString = ffComboBox_GetItemString( hwnd, i, itemString );
        if( gotString ) {
            if( boost::algorithm::iequals( s, itemString ) ) {
                return i;
            }
        }
    }
    return -1;
}

// Makes sure windows messages keep passing around
inline void PumpMessages() {
    MSG msg;
    while( ::PeekMessage( &msg, NULL, NULL, NULL, PM_REMOVE ) ) {
        ::TranslateMessage( &msg );
        ::DispatchMessage( &msg );
    }
}

#pragma comment( lib, "Version.lib" )

// Returns the version of a dll or exe as a 64 bit integer, where the major, minor, etc. version are encoded in
// 16bit chunks with the most significant 4 bytes being the major version number.
// ex.
// Version 2.5.0.0 = 0x0002000500000000
inline __int64 GetVersion( const frantic::tstring& path ) {
    DWORD garbage = 0;
    DWORD infoSize = GetFileVersionInfoSize( path.c_str(), &garbage );
    if( infoSize == 0 )
        throw std::runtime_error( "GetVersion: GetFileVersionInfoSize(" + frantic::strings::to_string( path ) +
                                  ") failed because:\n" + frantic::win32::GetLastErrorMessageA() );

    boost::scoped_array<char> pBuffer( new char[infoSize] );
    if( !GetFileVersionInfo( path.c_str(), NULL, infoSize, pBuffer.get() ) )
        throw std::runtime_error( "GetVersion: GetFileVersionInfo(" + frantic::strings::to_string( path ) +
                                  ") failed because:\n" + frantic::win32::GetLastErrorMessageA() );

    VS_FIXEDFILEINFO* ffi;
    UINT ffiSize;

    if( !VerQueryValue( pBuffer.get(), _T("\\"), (LPVOID*)&ffi, &ffiSize ) )
        throw std::runtime_error( "GetVersion: VerQueryValue() failed" );

    return ( (__int64)ffi->dwProductVersionMS << 32 ) | (__int64)ffi->dwProductVersionLS;
}

inline frantic::tstring GetExecutableVersion( const frantic::tstring& fileName ) {
    // Get version info for this app and put it into dialog box member variables
    frantic::tstring versionString = _T("unknown");

    // find size of ver info block
    DWORD dwVerHnd = 0;

    // Create temporary buffer ( so that arg1 is not const char* )
    TCHAR* tBuf;
    tBuf = (TCHAR*)malloc( (size_t)( fileName.length() + 1 ) * sizeof( TCHAR ) );
    _tcscpy( tBuf, fileName.c_str() );

    DWORD dwVerInfoSize = ::GetFileVersionInfoSize( tBuf, &dwVerHnd );
    if( ( dwVerInfoSize > 0 ) && ( dwVerHnd == 0 ) ) {
        // allocate memory to hold it
        HGLOBAL hMem = ::GlobalAlloc( GMEM_MOVEABLE, dwVerInfoSize );
        LPVOID pVerInfo = ::GlobalLock( hMem );

        // get the block
        if( ::GetFileVersionInfo( tBuf, dwVerHnd, dwVerInfoSize, pVerInfo ) ) {

            struct SLangCP {
                WORD wLang;
                WORD wCP;
            } * pTrans;

            ::UINT nVersionLen = 0;
            // Read the list of languages and code pages.  We will use only the	first one, for now.
            ::VerQueryValue( pVerInfo, _T("\\VarFileInfo\\Translation"), (LPVOID*)&pTrans, &nVersionLen );

            TCHAR szLookup[64];

            // Get the file version and if valid, write it to the dialog member	variable
            wsprintf( szLookup, _T("\\StringFileInfo\\%04x%04x\\ProductVersion"), pTrans[0].wLang, pTrans[0].wCP );
            LPTSTR pszVersion = NULL;
            BOOL bRetCode = ::VerQueryValue( pVerInfo, szLookup, (LPVOID*)&pszVersion, &nVersionLen );
            if( bRetCode && ( nVersionLen > 0 ) && ( pszVersion != NULL ) ) {
                versionString = pszVersion;
            }
        }
        ::GlobalUnlock( hMem );
        ::GlobalFree( hMem );
        delete( tBuf );
    }
    return versionString;
}

// Uses types from <winnt.h>
//
// The purpose of this routine is to parse a Windows "Portable Executable" file (a .dll or .exe) and decide if
// it is 32-bit or 64-bit.
inline bool IsFile64BitDLLorEXE( const frantic::tstring& filename ) {
    namespace fs = boost::filesystem;

    fs::path p( filename );
    if( !fs::exists( p ) )
        throw std::runtime_error( "IsFile64BitDLLorEXE: File \"" + frantic::strings::to_string( filename ) +
                                  "\" does not exist." );

    fs::ifstream fin;
    fin.open( p, std::ios::in | std::ios::binary );

    // Verify this is an exe or dll by checking for the magic signature of the first chunk, which is
    // the legacy DOS header.
    std::vector<char> buffer( sizeof( IMAGE_DOS_HEADER ), '\0' );
    fin.read( &buffer[0], sizeof( IMAGE_DOS_HEADER ) );

    IMAGE_DOS_HEADER* pDosHeader = reinterpret_cast<IMAGE_DOS_HEADER*>( &buffer[0] );
    if( pDosHeader->e_magic != IMAGE_DOS_SIGNATURE )
        throw std::runtime_error( "IsFile64BitDLLorEXE: \"" + frantic::strings::to_string( filename ) +
                                  "\" doesn't have a valid DOS signature value, so is not a PE or not a dll/exe" );

    // The DOS header has a relative pointer to the NT header.
    fin.seekg( pDosHeader->e_lfanew, std::ios::beg );

    buffer.resize( sizeof( IMAGE_NT_HEADERS ), '\0' );
    fin.read( &buffer[0], sizeof( IMAGE_NT_HEADERS ) );

    if( !fin )
        throw std::runtime_error( "IsFile64BitDLLorEXE: \"" + frantic::strings::to_string( filename ) +
                                  "\" doesn't have a valid NT signature value, so is not a PE" );

    fin.close();

    // Check if the NT header's magic number matches the 64bit header indicator.
    IMAGE_NT_HEADERS* pNTHeader = reinterpret_cast<IMAGE_NT_HEADERS*>( &buffer[0] );
    return pNTHeader->OptionalHeader.Magic == IMAGE_NT_OPTIONAL_HDR64_MAGIC;
}

#pragma warning( pop )

} // namespace win32
} // namespace frantic
