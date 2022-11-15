// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( _WIN32 )

#include <frantic/diagnostics/timeout_tracker.hpp>
#include <frantic/files/files.hpp>
#include <frantic/win32/child_process.hpp>
#include <frantic/win32/utility.hpp>
#include <frantic/win32/wrap.hpp>

using namespace std;
using namespace frantic;

namespace frantic {
namespace process {

child_process::child_process( bool hideWindow, bool terminateOnExit, bool controlStdOut ) {
    m_hProcess = 0;
    m_hReadHandle = 0;
    m_hWriteHandle = 0;
    m_hJobObject = 0;

    reset();

    m_terminateOnExit = terminateOnExit;
    m_controlStdOut = controlStdOut;
    m_hideWindow = hideWindow;
    m_binaryStdOut = false;

    m_aggressiveBuffer = false;

    // Default 4K buffer
    m_stdoutBuffer.resize( 4096 );
}

child_process::~child_process() { reset(); }

void child_process::reset() {
    if( m_hProcess != 0 ) {
        if( m_terminateOnExit )
            this->terminate( 1234 );

        CloseHandle( m_hProcess );
        CloseHandle( m_hReadHandle );
        CloseHandle( m_hWriteHandle );
        if( m_hJobObject != 0 )
            CloseHandle( m_hJobObject );
    }
    m_stdoutBufferStart = 0;
    m_stdoutBufferCount = 0;
    m_stdoutBufferSkipEmptyLine = false;

    m_hProcess = 0;
    m_hReadHandle = 0;
    m_hWriteHandle = 0;
    m_hJobObject = 0;

    m_pid = 0;

    m_useJobObject = false;
    m_createNewConsole = false;
}

void child_process::set_aggressiveBuffer( bool value ) { m_aggressiveBuffer = value; }

void child_process::set_bufferSize( int size ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::set_bufferSize: called after process was already launched." );

    m_stdoutBuffer.resize( size );
}

void child_process::set_terminateOnExit( bool value ) { m_terminateOnExit = value; }

void child_process::set_controlStdOut( bool value ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::set_controlStdOut: called after process was already launched." );

    m_controlStdOut = value;
}

void child_process::set_controlStdIn( bool value ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::set_controlStdIn: called after process was already launched." );

    m_controlStdIn = value;
}

void child_process::set_binaryStdOut( bool value ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::set_controlStdOut: called after process was already launched." );

    m_binaryStdOut = value;
}

void child_process::set_hideWindow( bool value ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::set_hideWindow: called after process was already launched." );

    m_hideWindow = value;
}

void child_process::set_useJobObject( bool value ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::set_useJobObject: called after process was already launched." );

    m_useJobObject = value;
}

void child_process::set_createNewConsole( bool value ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::set_createNewConsole: called after process was already launched." );

    m_createNewConsole = value;
}

void child_process::launch_as( const std::string& executable, const std::string& arguments,
                               const std::string& startupdir, const std::string& username, const std::string& domain,
                               const std::string& password ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::launch_as: called after a process was already launched." );

    if( m_controlStdOut ) {
        throw runtime_error( "child_process::launch_as: cannot control stdout when calling this function" );
    }

    if( m_controlStdIn ) {
        throw runtime_error( "child_process::launch_as: cannot control stdin when calling this function" );
    }

    if( m_useJobObject ) {
        m_hJobObject = CreateJobObject( NULL, NULL );
        if( m_hJobObject == NULL )
            throw runtime_error( "child_process::launch_as: Error creating JobObject: " +
                                 win32::GetLastErrorMessageA() );
    }

    wchar_t wUsername[256];
    wchar_t wDomain[256];
    wchar_t wPassword[256];
    wchar_t wCmdline[256];
    wchar_t wStartupdir[256];

    STARTUPINFOW StartupInfo;
    PROCESS_INFORMATION ProcessInfo;

    //	Initialize the structures to zero
    memset( &StartupInfo, 0, sizeof( StartupInfo ) );
    memset( &ProcessInfo, 0, sizeof( ProcessInfo ) );

    //	Initialize the startupinfo
    StartupInfo.cb = sizeof( StartupInfo );
    StartupInfo.dwFlags = STARTF_USESHOWWINDOW;
    StartupInfo.wShowWindow = m_hideWindow ? SW_HIDE : SW_SHOW;

    //	Set up command line and startup directory
    string cmdline;
    if( executable == "" )
        cmdline = arguments;
    else
        cmdline = "\"" + executable + "\" " + arguments;

    string launchStartupDir = startupdir;
    if( launchStartupDir == "" && executable != "" )
        launchStartupDir = files::directory_from_path( executable );

    //	Convert c strings to wide char arrays required by CreateProcessWithLogonW function
    MultiByteToWideChar( CP_ACP, 0, username.c_str(), -1, wUsername, sizeof( wUsername ) / sizeof( wUsername[0] ) );
    MultiByteToWideChar( CP_ACP, 0, domain.c_str(), -1, wDomain, sizeof( wDomain ) / sizeof( wDomain[0] ) );
    MultiByteToWideChar( CP_ACP, 0, password.c_str(), -1, wPassword, sizeof( wPassword ) / sizeof( wPassword[0] ) );
    MultiByteToWideChar( CP_ACP, 0, cmdline.c_str(), -1, wCmdline, sizeof( wCmdline ) / sizeof( wCmdline[0] ) );
    MultiByteToWideChar( CP_ACP, 0, launchStartupDir.c_str(), -1, wStartupdir,
                         sizeof( wStartupdir ) / sizeof( wStartupdir[0] ) );

    //	Set required flags
    DWORD dwLogonFlags = LOGON_WITH_PROFILE;
    DWORD dwCreationFlags = CREATE_SUSPENDED + ( m_createNewConsole ? CREATE_NEW_CONSOLE : 0 );

    //	Create the child process
    if( !CreateProcessWithLogonW( wUsername, wDomain, wPassword, dwLogonFlags, NULL, wCmdline, dwCreationFlags, NULL,
                                  wStartupdir, &StartupInfo, &ProcessInfo ) ) {
        if( m_hJobObject != 0 )
            CloseHandle( m_hJobObject );
        m_hJobObject = 0;
        throw runtime_error( "child_process::launch_as: Error starting \"" +
                             ( executable == "" ? arguments : executable ) + "\" for user " + username + ": " +
                             win32::GetLastErrorMessageA() );
    }

    m_hProcess = ProcessInfo.hProcess;
    m_pid = ProcessInfo.dwProcessId;

    //	Put the process in the job object
    if( m_useJobObject ) {
        if( !AssignProcessToJobObject( m_hJobObject, m_hProcess ) ) {
            CloseHandle( m_hJobObject );
            CloseHandle( ProcessInfo.hThread );
            m_hJobObject = 0;
            throw runtime_error( "child_process::launch_as: Failed to assign process \"" +
                                 ( executable == "" ? arguments : executable ) +
                                 "\" to job object: " + win32::GetLastErrorMessageA() );
        }
    }

    //	Make the child process inherit the priority class
    setpriorityclass( GetPriorityClass( GetCurrentProcess() ) );

    ResumeThread( ProcessInfo.hThread );
    CloseHandle( ProcessInfo.hThread );
}

void child_process::launch_with_environment( const frantic::tstring& executable, const frantic::tstring& arguments,
                                             const frantic::tstring& startupdir,
                                             const std::map<frantic::tstring, frantic::tstring>* environment ) {
    if( m_hProcess != 0 )
        throw runtime_error( "child_process::launch: called again after a process was already launched." );

    if( m_useJobObject ) {

        m_hJobObject = CreateJobObject( NULL, NULL );
        if( m_hJobObject == NULL )
            throw runtime_error( "child_process::launch: Error creating JobObject: " + win32::GetLastErrorMessageA() );
    }

    HANDLE PipeReadHandle = 0;
    HANDLE PipeWriteHandle = 0;
    HANDLE ChildReadHandle = 0;
    HANDLE ChildWriteHandle = 0;
    PROCESS_INFORMATION ProcessInfo;
    SECURITY_ATTRIBUTES SecurityAttributes;
    STARTUPINFO StartupInfo;

    // Initialize the structures to zero
    memset( &StartupInfo, 0, sizeof( StartupInfo ) );
    memset( &ProcessInfo, 0, sizeof( ProcessInfo ) );
    memset( &SecurityAttributes, 0, sizeof( SecurityAttributes ) );

    //	Make a pipe to connect to the renderer's stdout
    SecurityAttributes.nLength = sizeof( SECURITY_ATTRIBUTES );
    SecurityAttributes.bInheritHandle = TRUE;
    SecurityAttributes.lpSecurityDescriptor = 0;

    if( m_controlStdOut ) {
        // Create a pipe with the default number of reserved bytes
        if( !CreatePipe( &PipeReadHandle, &PipeWriteHandle, &SecurityAttributes, 0 ) ) {
            std::string errorMessage = win32::GetLastErrorMessageA();
            if( m_hJobObject != 0 )
                CloseHandle( m_hJobObject );
            m_hJobObject = 0;
            throw runtime_error( "Error creating pipe for process i/o: " + errorMessage );
        }

        if( PipeWriteHandle == 0 || PipeWriteHandle == INVALID_HANDLE_VALUE ) {
            std::string errorMessage = win32::GetLastErrorMessageA();
            if( m_hJobObject != 0 )
                CloseHandle( m_hJobObject );
            m_hJobObject = 0;
            if( PipeReadHandle != 0 )
                CloseHandle( PipeReadHandle );
            if( PipeWriteHandle != 0 )
                CloseHandle( PipeWriteHandle );
            throw runtime_error( "CreatePipe returned success, but produced invalid handles: " + errorMessage );
        }
    }

    if( m_controlStdIn ) {
        // Create a pipe to allow writing to the stdin of the child process
        if( !CreatePipe( &ChildReadHandle, &ChildWriteHandle, &SecurityAttributes, 0 ) ) {
            std::string errorMessage = win32::GetLastErrorMessageA();
            if( m_hJobObject != 0 )
                CloseHandle( m_hJobObject );
            m_hJobObject = 0;
            throw runtime_error( "Error creating pipe for process i/o: " + errorMessage );
        }
        if( ChildWriteHandle == 0 || ChildWriteHandle == INVALID_HANDLE_VALUE ) {
            std::string errorMessage = win32::GetLastErrorMessageA();
            if( ChildWriteHandle != 0 )
                CloseHandle( ChildWriteHandle );
            if( ChildReadHandle != 0 )
                CloseHandle( ChildReadHandle );
            throw runtime_error( "CreatePipe returned success, but produced invalid handles: " + errorMessage );
        }

        SetHandleInformation( ChildWriteHandle, HANDLE_FLAG_INHERIT, 0 );
    }

    //	Initialize the startupinfo
    StartupInfo.cb = sizeof( STARTUPINFO );
    StartupInfo.dwFlags = STARTF_USESHOWWINDOW | ( m_controlStdOut ? STARTF_USESTDHANDLES : 0 );
    StartupInfo.wShowWindow = m_hideWindow ? SW_HIDE : SW_SHOW;
    if( m_controlStdOut ) {
        StartupInfo.hStdOutput = PipeWriteHandle;
        StartupInfo.hStdError = PipeWriteHandle;
    }
    if( m_controlStdIn ) {
        StartupInfo.hStdInput = ChildReadHandle;
    }

    //----------------------------------------------------------------------------
    //	Create the child process.
    //----------------------------------------------------------------------------
    frantic::tstring cmdline;
    if( executable == _T("") )
        cmdline = arguments;
    else
        cmdline = _T("\"") + executable + _T("\" ") + arguments;

    frantic::tstring launchStartupDir = startupdir;

    if( launchStartupDir == _T("") && executable != _T("") )
        launchStartupDir = files::directory_from_path( executable );

    // Prepare the child process environment variable block - if needed
    LPTSTR lpszEnvironmentBlock = NULL;
    frantic::tstring envBuf = _T("");
    if( environment != NULL ) {
        for( std::map<frantic::tstring, frantic::tstring>::const_iterator it = environment->begin();
             it != environment->end(); it++ )
            envBuf = it->first + _T("=") + it->second + _T("\0");
        lpszEnvironmentBlock = const_cast<LPTSTR>( envBuf.c_str() );
    }

    if( !CreateProcess( NULL, const_cast<LPTSTR>( cmdline.c_str() ),
                        NULL,                                                               // lpProcessAttributes
                        NULL,                                                               // lpThreadAttributes
                        m_controlStdOut ? TRUE : FALSE,                                     // bInheritHandles
                        CREATE_SUSPENDED + ( m_createNewConsole ? CREATE_NEW_CONSOLE : 0 ), // dwCreationFlags
                        lpszEnvironmentBlock, // use parent's environment block
                        launchStartupDir == _T("") ? NULL : launchStartupDir.c_str(), &StartupInfo, &ProcessInfo ) ) {
        if( PipeReadHandle != 0 )
            CloseHandle( PipeReadHandle );
        if( PipeWriteHandle != 0 )
            CloseHandle( PipeWriteHandle );
        if( ChildWriteHandle != 0 )
            CloseHandle( ChildWriteHandle );
        if( ChildReadHandle != 0 )
            CloseHandle( ChildReadHandle );
        if( m_hJobObject != 0 )
            CloseHandle( m_hJobObject );
        m_hJobObject = 0;
        throw runtime_error(
            "child_process::launch: Error starting \"" +
            ( frantic::strings::to_string( executable ) == "" ? frantic::strings::to_string( arguments )
                                                              : frantic::strings::to_string( executable ) ) +
            "\" in \"" + frantic::strings::to_string( launchStartupDir ) + "\" : " + win32::GetLastErrorMessageA() );
    }

    if( PipeWriteHandle != 0 )
        CloseHandle( PipeWriteHandle );
    if( ChildReadHandle != 0 )
        CloseHandle( ChildReadHandle );

    m_hProcess = ProcessInfo.hProcess;
    m_pid = ProcessInfo.dwProcessId;
    m_hReadHandle = PipeReadHandle;
    m_hWriteHandle = ChildWriteHandle;

    // Put the process in the job object
    if( m_useJobObject ) {
        if( !AssignProcessToJobObject( m_hJobObject, m_hProcess ) ) {
            TerminateProcess( m_hProcess, 0 );
            CloseHandle( m_hProcess );
            m_hProcess = 0;
            CloseHandle( m_hJobObject );
            CloseHandle( ProcessInfo.hThread );
            m_hJobObject = 0;
            throw runtime_error( "child_process::launch: Failed to assign process \"" +
                                 ( frantic::strings::to_string( executable ) == ""
                                       ? frantic::strings::to_string( arguments )
                                       : frantic::strings::to_string( executable ) ) +
                                 "\" to job object: " + win32::GetLastErrorMessageA() );
        }
    }

    // Make the child process inherit the priority class
    setpriorityclass( GetPriorityClass( GetCurrentProcess() ) );

    ResumeThread( ProcessInfo.hThread );
    CloseHandle( ProcessInfo.hThread );
}

void child_process::launch( const frantic::tstring& executable, const frantic::tstring& arguments,
                            const frantic::tstring& startupdir ) {
    launch_with_environment( executable, arguments, startupdir, NULL );
}

void child_process::terminate( int exitCode ) {
    if( m_useJobObject )
        TerminateJobObject( m_hJobObject, exitCode );
    else
        TerminateProcess( m_hProcess, exitCode );
}

void child_process::close_windows() { post_wm_message( WM_CLOSE, 0, 0 ); }

void child_process::post_wm_message( UINT Msg, WPARAM wParam, LPARAM lParam ) {
    vector<HWND> windows;
    get_windows( windows );
    for( size_t i = 0; i < windows.size(); ++i ) {
        PostMessage( windows[i], Msg, wParam, lParam );
    }
}

namespace detail {
struct child_processWindowAdder {
    set<DWORD> pids;
    std::vector<HWND>* windows;
};
} // namespace detail

BOOL CALLBACK child_processGetWindowsProc( HWND hwnd, LPARAM lParam ) {
    detail::child_processWindowAdder* data = (detail::child_processWindowAdder*)lParam;

    DWORD windowProcessId;
    if( GetWindowThreadProcessId( hwnd, &windowProcessId ) ) {

        if( data->pids.find( windowProcessId ) != data->pids.end() ) {
            data->windows->push_back( hwnd );
            return TRUE;
        }
    }
    return TRUE;
}

void child_process::get_windows( std::vector<HWND>& windows ) const {
    if( m_pid != 0 ) {
        detail::child_processWindowAdder data;
        get_processes( std::inserter( data.pids, data.pids.end() ) );
        data.windows = &windows;
        EnumWindows( child_processGetWindowsProc, reinterpret_cast<LPARAM>( &data ) );
    }
}

void child_process::get_popups( std::vector<HWND>& popups ) const {
    popups.clear();
    std::vector<HWND> hwnds;
    get_windows( hwnds );

    for( size_t i = 0; i < hwnds.size(); ++i ) {
        frantic::tstring text = win32::ffGetWindowText( hwnds[i] ),
                         windowClass = win32::ffGetWindowClassName( hwnds[i] );

        long style = GetWindowLong( hwnds[i], GWL_STYLE );

        // The dialog window class is "#32770", so consider any window with this class,
        // a non-empty title and WS_POPUP or WS_OVERLAPPED to be a blocking popup window.
        if( text == _T("") || windowClass != _T("#32770") || !( style & WS_POPUP || style & WS_OVERLAPPED ) ||
            ( style & WS_CHILD ) )
            continue;

        vector<HWND> childHwnds;
        try {
            win32::GetChildWindows( hwnds[i], childHwnds );
        } catch( const std::exception& ) {
            // LogMessage( "Window \"" + text + "\": " + e.what(), ll_warning );
            continue;
        }

        // Apply a simple heuristic to the child windows to determine whether it's a popup dialog
        if( childHwnds.size() > 20 ) // At most 20 controls
            continue;
        int buttonCount = 0, messageCount = 0, otherCount = 0;
        for( size_t j = 0; j < childHwnds.size(); ++j ) {
            if( ( GetWindowLong( childHwnds[j], GWL_STYLE ) & WS_VISIBLE ) ) {
                frantic::tstring childClass = win32::ffGetWindowClassName( childHwnds[j] );
                if( childClass == _T("Button") )
                    buttonCount++;
                else if( childClass == _T("Static") || childClass == _T("RICHEDIT") )
                    messageCount++;
                else
                    otherCount++;
            }
        }
        // Need a button and a message to consider it a popup
        if( buttonCount == 0 || messageCount == 0 )
            continue;

        // This is considered a popup, so add it to the list
        popups.push_back( hwnds[i] );
    }
}

// Waits for the process to exit.  Returns true if it did exit, false if it timed out
bool child_process::waitforexit( int timeout ) { return win32::WaitForProcessToTerminate( m_hProcess, timeout ); }

int child_process::getexitcode() { return win32::GetExitCodeProcess( m_hProcess ); }

bool child_process::running() { return m_hProcess != 0 && win32::GetExitCodeProcess( m_hProcess ) == STILL_ACTIVE; }

bool child_process::any_running() {
    bool running = false;
    vector<DWORD> pids;
    get_processes( std::back_inserter( pids ) );
    for( unsigned i = 0; i < pids.size(); ++i ) {
        HANDLE hProcess = OpenProcess( PROCESS_QUERY_INFORMATION, false, pids[i] );
        if( hProcess != 0 && win32::GetExitCodeProcess( hProcess ) == STILL_ACTIVE )
            running = true;

        if( hProcess != 0 )
            CloseHandle( hProcess );

        if( running )
            break;
    }

    return running;
}

void child_process::get_process_pids( std::vector<DWORD>& pids ) { get_processes( std::back_inserter( pids ) ); }

bool child_process::setpriorityclass( DWORD dwPriorityClass ) {
    return SetPriorityClass( m_hProcess, dwPriorityClass ) != 0;
}

int child_process::getpriorityclass() { return GetPriorityClass( m_hProcess ); }

bool child_process::getprocessaffinitymask( DWORD_PTR& dwProcessAffinityMask ) {
    DWORD_PTR dwSystemAffinityMask;
    return GetProcessAffinityMask( m_hProcess, &dwProcessAffinityMask, &dwSystemAffinityMask ) != 0;
}

bool child_process::setprocessaffinitymask( DWORD_PTR dwProcessAffinityMask ) {
    return SetProcessAffinityMask( m_hProcess, dwProcessAffinityMask ) != 0;
}

bool child_process::is_stdout_data() {
    if( m_hReadHandle == 0 )
        return false;

    if( m_stdoutBufferCount == m_stdoutBufferStart )
        GetDataForBuffer( 0 );

    if( m_stdoutBufferCount == m_stdoutBufferStart )
        return false;

    if( m_stdoutBufferSkipEmptyLine ) {
        if( m_stdoutBufferCount - m_stdoutBufferStart < 2 )
            return false;
        // If there's an end of line in the buffer, we will return a line
        return ( memchr( &m_stdoutBuffer[m_stdoutBufferStart + 1], '\r',
                         m_stdoutBufferCount - m_stdoutBufferStart - 1 ) != NULL ||
                 memchr( &m_stdoutBuffer[m_stdoutBufferStart + 1], '\n',
                         m_stdoutBufferCount - m_stdoutBufferStart - 1 ) != NULL );
    } else {
        // If there's an end of line in the buffer, we will return a line
        return (
            memchr( &m_stdoutBuffer[m_stdoutBufferStart], '\r', m_stdoutBufferCount - m_stdoutBufferStart ) != NULL ||
            memchr( &m_stdoutBuffer[m_stdoutBufferStart], '\n', m_stdoutBufferCount - m_stdoutBufferStart ) != NULL );
    }
}

bool child_process::getstdoutbytes( int timeout, char* bytes, int bytesCount ) {
    if( m_hReadHandle == 0 )
        throw runtime_error(
            "child_process::getstdoutbytes: Should not call this function when not controlling stdout" );

    if( !m_binaryStdOut )
        throw runtime_error( "child_process::getstdoutbytes: Binary stdout is disabled, use getstdoutline instead" );

    if( ExtractDataFromBuffer( bytes, bytesCount ) )
        return true;

    GetDataForBuffer( timeout );

    return ExtractDataFromBuffer( bytes, bytesCount );
}

bool child_process::getstdoutline( int timeout, std::string& outString ) {
    if( m_hReadHandle == 0 )
        throw runtime_error(
            "child_process::getstdoutline: Should not call this function when not controlling stdout" );

    if( m_binaryStdOut )
        throw runtime_error( "child_process::getstdoutline: Binary stdout is enabled, use getstdoutbytes instead" );

    if( ExtractLineFromBuffer( outString ) )
        return true;

    GetDataForBuffer( timeout );

    return ExtractLineFromBuffer( outString );
}

void child_process::GetDataForBuffer( int timeout ) {
    diagnostics::timeout_tracker theTimeout( timeout );

    // Reset the buffer start to 0
    if( m_stdoutBufferStart != 0 ) {
        memcpy( &m_stdoutBuffer[0], &m_stdoutBuffer[m_stdoutBufferStart], m_stdoutBufferCount - m_stdoutBufferStart );
        m_stdoutBufferCount -= m_stdoutBufferStart;
        m_stdoutBufferStart = 0;
    }

    DWORD bytesAvail = 0, bytesRead = 0, bytesLeftThisMessage = 0;

    while( bytesAvail == 0 && !theTimeout.timed_out() ) {
        PeekNamedPipe( m_hReadHandle, NULL, 0, &bytesRead, &bytesAvail, &bytesLeftThisMessage );
        theTimeout.sleep( 50 );
    }

    if( bytesAvail == 0 )
        return;

    if( m_aggressiveBuffer ) {
        while( m_stdoutBufferCount < (int)m_stdoutBuffer.size() * 4 / 5 ) {
            if( !ReadFile( m_hReadHandle, &m_stdoutBuffer[m_stdoutBufferCount],
                           (DWORD)m_stdoutBuffer.size() - m_stdoutBufferCount, &bytesRead, NULL ) )
                return;
            m_stdoutBufferCount += bytesRead;
        }
    } else {
        if( !ReadFile( m_hReadHandle, &m_stdoutBuffer[m_stdoutBufferCount],
                       (DWORD)m_stdoutBuffer.size() - m_stdoutBufferCount, &bytesRead, NULL ) )
            return;
        m_stdoutBufferCount += bytesRead;
    }
}

// writes a string to the standard in handle of the process
DWORD child_process::WriteToStdIn( std::string line ) {
    if( !m_controlStdIn )
        throw runtime_error( "child_process::WriteToStdIn called without enabling control std in on this process." );
    DWORD numBytesWritten = 0;
    if( m_hWriteHandle != 0 ) {
        line.append( "\n" ); // Appends a newline so the stdin is read by a console application.
        if( !WriteFile( m_hWriteHandle, line.c_str(), (DWORD)line.length(), &numBytesWritten, NULL ) ) {
            throw runtime_error( "child_process::WriteToStdIn: There was an error writing to the child process. " );
        }
    }
    return numBytesWritten;
}

bool child_process::ExtractDataFromBuffer( char* bytes, int bytesCount ) {
    if( m_stdoutBufferCount - m_stdoutBufferStart >= bytesCount ) {
        memcpy( bytes, &m_stdoutBuffer[m_stdoutBufferStart], bytesCount );
        m_stdoutBufferStart += bytesCount;
        return true;
    } else
        return false;
}

// Extracts a string from the stdout buffer
bool child_process::ExtractLineFromBuffer( std::string& outString ) {
    // If there was a half-EOL at the end of the string, which is empty, try to remove that empty line.
    // If it can not remove the empty line, do not get a string
    if( m_stdoutBufferSkipEmptyLine ) {
        m_stdoutBufferSkipEmptyLine = false;
        ExtractLineFromBuffer( outString );
        if( m_stdoutBufferSkipEmptyLine )
            return false;
    }

    // Find the first carriage return
    void* LF = memchr( &m_stdoutBuffer[m_stdoutBufferStart], '\n', m_stdoutBufferCount - m_stdoutBufferStart );
    void* CR = memchr( &m_stdoutBuffer[m_stdoutBufferStart], '\r', m_stdoutBufferCount - m_stdoutBufferStart );

    if( LF == 0 && CR == 0 )
        return false;

    // Get the distance to the closest end of line character
    int size = 0;
    if( LF == 0 )
        size = int( static_cast<char*>( CR ) - &m_stdoutBuffer[m_stdoutBufferStart] );
    else if( CR == 0 )
        size = int( static_cast<char*>( LF ) - &m_stdoutBuffer[m_stdoutBufferStart] );
    else
        size = min( int( static_cast<char*>( CR ) - &m_stdoutBuffer[m_stdoutBufferStart] ),
                    int( static_cast<char*>( LF ) - &m_stdoutBuffer[m_stdoutBufferStart] ) );

    // Half a linefeed at the end of the buffer is bad
    if( size + 1 == m_stdoutBufferCount - m_stdoutBufferStart )
        m_stdoutBufferSkipEmptyLine = true;

    // Get the string
    outString = string( &m_stdoutBuffer[m_stdoutBufferStart], size );

    // Move the buffer to exclude that line
    if( !m_stdoutBufferSkipEmptyLine ) {
        size++;
        if( size < m_stdoutBufferCount ) { // Count "\r\n" and "\n\r" as end of lines too
            if( ( m_stdoutBuffer[m_stdoutBufferStart + size] == '\r' &&
                  m_stdoutBuffer[m_stdoutBufferStart + size - 1] == '\n' ) ||
                ( m_stdoutBuffer[m_stdoutBufferStart + size] == '\n' &&
                  m_stdoutBuffer[m_stdoutBufferStart + size - 1] == '\r' ) )
                size++;
        }
    }

    m_stdoutBufferStart += size;

    return true;
}

} // namespace process
} // namespace frantic

#endif // defined( _WIN32 )
