// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#if !defined( _WIN32 )

#include <fstream>
#include <stdexcept>
#include <vector>

#include <fcntl.h>
#include <sys/time.h>

#include <frantic/diagnostics/timeout_tracker.hpp>
#include <frantic/misc/string_functions.hpp>
#include <frantic/unix/child_process.hpp>

using std::string;

using namespace std;
using namespace frantic;

// ----------------------------------------------------------------------------------
// Begin Child Process Class

namespace frantic {
namespace process {

child_process::child_process() {
    m_terminateOnExit = false;
    m_controlStdOut = true;
    m_hideWindow = false;
    m_binaryStdOut = false;

    m_aggressiveBuffer = false;

    // Default 4K buffer
    m_stdoutBuffer.resize( 4096 );

    m_dataPipePtoC[0] = -1;
    m_dataPipePtoC[1] = -1;
    m_dataPipeCtoP[0] = -1;
    m_dataPipeCtoP[1] = -1;
    m_childStatus = -1;
    reset();
}

child_process::~child_process() { reset(); }

void child_process::reset() {
    m_stdoutBufferStart = 0;
    m_stdoutBufferCount = 0;
    m_stdoutBufferSkipEmptyLine = false;

    if( m_dataPipePtoC[0] != -1 ) {
        close( m_dataPipePtoC[0] );
        m_dataPipePtoC[0] = -1;
    }
    if( m_dataPipePtoC[1] != -1 ) {
        close( m_dataPipePtoC[1] );
        m_dataPipePtoC[1] = -1;
    }
    if( m_dataPipeCtoP[0] != -1 ) {
        close( m_dataPipeCtoP[0] );
        m_dataPipeCtoP[0] = -1;
    }
    if( m_dataPipeCtoP[1] != -1 ) {
        close( m_dataPipeCtoP[1] );
        m_dataPipeCtoP[1] = -1;
    }

    m_pid = -1;

    m_useJobObject = false;
    m_createNewConsole = false;
}

int child_process::Get_pID() { return m_pid; }

void child_process::set_aggressiveBuffer( bool value ) { m_aggressiveBuffer = value; }

void child_process::set_bufferSize( int size ) {
    if( m_pid > 0 )
        perror( "child_process::set_bufferSize: called after process was already launched." );

    m_stdoutBuffer.resize( size );
    // cout << "ERROR: Resize doens't work right now." << endl;
}

void child_process::set_terminateOnExit( bool value ) {
    if( m_pid > 0 )
        perror( "child_process::set_terminateOnExit: called after process was already launched." );

    m_terminateOnExit = value;
}

void child_process::set_controlStdOut( bool value ) {
    if( m_pid > 0 )
        perror( "child_process::set_controlStdOut: called after process was already launched." );

    m_controlStdOut = value;
}

void child_process::set_binaryStdOut( bool value ) {
    if( m_pid > 0 )
        perror( "child_process::set_controlStdOut: called after process was already launched." );

    m_binaryStdOut = value;
}

void child_process::set_hideWindow( bool value ) {
    if( m_pid > 0 )
        perror( "child_process::set_hideWindow: called after process was already launched." );

    m_hideWindow = value;
}

void child_process::set_useJobObject( bool value ) {
    if( m_pid > 0 )
        perror( "child_process::set_useJobObject: called after process was already launched." );

    m_useJobObject = value;
}

void child_process::set_createNewConsole( bool value ) {
    if( m_pid > 0 )
        perror( "child_process::set_createNewConsole: called after process was already launched." );

    m_createNewConsole = value;
}

void child_process::launch_as( const std::string& executable, const std::string& arguments,
                               const std::string& startupdir, const std::string& username, const std::string& domain,
                               const std::string& password ) {
    perror( "child_process::launch_as: not supported in Unix envirmonment." );

    /*if( m_pid != 0 )
      perror( "child_process::launch_as: called after a process was already launched." );

    if( m_controlStdOut ) {
      perror( "child_process::launch_as: cannot control stdout when calling this function" );
    }

    if( m_useJobObject ) {
      m_hJobObject = CreateJobObject( NULL, NULL );
      if( m_hJobObject == NULL )
        perror( "child_process::launch_as: Error creating JobObject: " + win32::GetLastErrorMessage() );
    }

    wchar_t					wUsername[256];
    wchar_t					wDomain[256];
    wchar_t					wPassword[256];
    wchar_t					wCmdline[256];
    wchar_t					wStartupdir[256];

    STARTUPINFOW			StartupInfo;
    PROCESS_INFORMATION		ProcessInfo;

    //	Initialize the structures to zero
    memset( &StartupInfo, 0, sizeof(StartupInfo) );
    memset( &ProcessInfo, 0, sizeof(ProcessInfo) );

    //	Initialize the startupinfo
    StartupInfo.cb           = sizeof(StartupInfo);
    StartupInfo.dwFlags      = STARTF_USESHOWWINDOW;
    StartupInfo.wShowWindow  = m_hideWindow ? SW_HIDE : SW_SHOW;

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
    MultiByteToWideChar( CP_ACP, 0, username.c_str(), -1, wUsername, sizeof(wUsername) / sizeof(wUsername[0]) );
    MultiByteToWideChar( CP_ACP, 0, domain.c_str(), -1, wDomain, sizeof(wDomain) / sizeof(wDomain[0]) );
    MultiByteToWideChar( CP_ACP, 0, password.c_str(), -1, wPassword, sizeof(wPassword) / sizeof(wPassword[0]) );
    MultiByteToWideChar( CP_ACP, 0, cmdline.c_str(), -1, wCmdline, sizeof(wCmdline) / sizeof(wCmdline[0]) );
    MultiByteToWideChar( CP_ACP, 0, launchStartupDir.c_str(), -1, wStartupdir, sizeof(wStartupdir) /
    sizeof(wStartupdir[0]) );

    //	Set required flags
    DWORD dwLogonFlags = LOGON_WITH_PROFILE;
    DWORD dwCreationFlags = CREATE_SUSPENDED + (m_createNewConsole ? CREATE_NEW_CONSOLE : 0);

    //	Create the child process
    if( !CreateProcessWithLogonW( wUsername, wDomain, wPassword, dwLogonFlags, NULL, wCmdline,
      dwCreationFlags, NULL, wStartupdir, &StartupInfo, &ProcessInfo ) ) {
      if( m_hJobObject != 0 )
        CloseHandle( m_hJobObject );
      m_hJobObject = 0;
      throw runtime_error( "child_process::launch_as: Error starting \"" + (executable == "" ? arguments : executable) +
    "\" for user " + username + ": " + win32::GetLastErrorMessage() );
    }

    m_hProcess = ProcessInfo.hProcess;
    m_pid = ProcessInfo.dwProcessId;

    //	Put the process in the job object
    if( m_useJobObject ) {
      if( !AssignProcessToJobObject( m_hJobObject, m_hProcess ) ) {
        CloseHandle( m_hJobObject );
        CloseHandle( ProcessInfo.hThread );
        m_hJobObject = 0;
        throw runtime_error( "child_process::launch_as: Failed to assign process \"" + (executable == "" ? arguments :
    executable) + "\" to job object: " + win32::GetLastErrorMessage() );
      }
    }

    //	Make the child process inherit the priority class
    setpriorityclass( GetPriorityClass( GetCurrentProcess() ) );

    ResumeThread( ProcessInfo.hThread );
    CloseHandle( ProcessInfo.hThread );*/
}

int child_process::launch( const std::string& executable, const std::string& arguments,
                           const std::string& startupdir ) // = "./" )
{
    cout << "Entered launch." << endl;
    cout << "exe: " << executable << endl;
    cout << "arg: " << arguments << endl;

    //  We should not have a processs id yet.
    if( m_pid > 0 )
        perror( "child_process::launch: called again after a process was already launched." );

    /*
//~ if( m_useJobObject ) {
//~ m_hJobObject = CreateJobObject( NULL, NULL );
//~ if( m_hJobObject == NULL )
//~ throw runtime_error( "child_process::launch: Error creating JobObject: " + win32::GetLastErrorMessage() );
//~ } */

    // Initialize the structures to zero
    //~ memset( &StartupInfo, 0, sizeof(StartupInfo) );
    //~ memset( &ProcessInfo, 0, sizeof(ProcessInfo) );
    //~ memset( &SecurityAttributes, 0, sizeof(SecurityAttributes) );

    // Setup up the pipes for data transfer between parent and child processes
    int rc_PtoC;
    int rc_CtoP;

    if( m_controlStdOut ) {
        // cout << "Creating Pipes.\n" << endl;
        // Create the pipe.
        rc_PtoC = pipe( m_dataPipePtoC );
        if( rc_PtoC == -1 ) {
            perror( "Error creating Parent to Child pipe for process i/o: " );
            exit( 1 );
        }

        rc_CtoP = pipe( m_dataPipeCtoP );
        if( rc_CtoP == -1 ) {
            perror( "Error creating Parent to Child pipe for process i/o: " );
            exit( 1 );
        }
    }

    //~ if( launchStartupDir == "" && executable != "" )
    //~ launchStartupDir = files::directory_from_path( executable );

    // Fork the process into parent and child and get the process id
    m_pid = fork();

    if( m_pid == -1 ) { //  Was there an error?
        // cout << "fork error.\n" << endl;
        perror( "fork error" );
    } else if( m_pid == 0 ) { //  We're in the child process

        // ofstream fout( "debug.txt" );

        std::vector<std::string> argVector;
        argVector.push_back( executable );
        frantic::strings::split( arguments, argVector );

        std::vector<char*> argArray;
        for( size_t i = 0; i < argVector.size(); i++ ) {
            argArray.push_back( const_cast<char*>( argVector[i].c_str() ) );
        }

        argArray.push_back( 0 );

        if( m_controlStdOut ) {
            //  close stdin
            close( m_dataPipePtoC[1] );
            m_dataPipePtoC[1] = -1;
            dup2( m_dataPipePtoC[0], STDIN_FILENO );

            //  close stdout
            close( m_dataPipeCtoP[0] );
            m_dataPipeCtoP[0] = -1;
            // dup2( data_pipe_CtoP[ 1 ], STDERR_FILENO );
            dup2( m_dataPipeCtoP[1], STDOUT_FILENO );
        }

        execvp( executable.c_str(), &argArray[0] );

        close( m_dataPipePtoC[0] );
        m_dataPipePtoC[0] = -1;
        close( m_dataPipeCtoP[1] );
        m_dataPipeCtoP[1] = -1;

        perror( "Child Process Error" );
        raise( SIGKILL );
        exit( 1 );

    } else { //  We're in the parent

        if( m_controlStdOut ) {
            // Close the unnecessary pipe ends
            close( m_dataPipeCtoP[1] );
            m_dataPipeCtoP[1] = -1;
            close( m_dataPipePtoC[0] );
            m_dataPipePtoC[0] = -1;
        }

        frantic::diagnostics::sleep_milliseconds( 50 );

        int status = -1;
        int pid = waitpid( m_pid, &status, WNOHANG );

        m_childStatus = status;

        if( pid == m_pid && WIFSIGNALED( status ) ) {
            throw std::runtime_error(
                "Failed to launch process \"" + executable +
                "\".  Check that this file exists and that it has proper execute permissions for the current user." );
        } else if( pid < 0 ) {
            throw std::runtime_error( "Failed to get status from waitpid when checking the status of \"" + executable +
                                      "\" execution.  Got return code: " + boost::lexical_cast<string>( pid ) );
        }
    }

    return m_pid;
}

void child_process::terminate( int exitCode ) {
    cout << "Terminating Child from Parent.\n" << endl;
    kill( m_pid, SIGKILL );
    // m_pid = 0;
    // sleep(1);
}

/*
void child_process::close_windows()
{
  post_wm_message( WM_CLOSE, 0, 0 );
}

void child_process::post_wm_message( UINT Msg, WPARAM wParam, LPARAM lParam )
{
  vector<HWND> windows;
  get_windows(windows);
  for( size_t i = 0; i < windows.size(); ++i ) {
    PostMessage( windows[i], Msg, wParam, lParam );
  }
}


namespace detail {
  struct child_processWindowAdder {
    set<DWORD> pids;
    std::vector<HWND>* windows;
  };
}

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

void child_process::get_windows( std::vector<HWND>& windows ) const
{
  if( m_pid != 0 ) {
    detail::child_processWindowAdder data;
    get_processes( std::inserter( data.pids, data.pids.end() ) );
    data.windows = &windows;
    EnumWindows( child_processGetWindowsProc, (LPARAM)&data );
  }
}*/

/*void child_process::get_popups( std::vector<HWND>& popups ) const
{
  popups.clear();
  std::vector<HWND> hwnds;
  get_windows(hwnds);

  for( size_t i = 0; i < hwnds.size(); ++i ) {
    string	text = win32::ffGetWindowText( hwnds[i] ),
        windowClass = win32::ffGetWindowClassName( hwnds[i] );

    long style = GetWindowLong( hwnds[i], GWL_STYLE );

    // The dialog window class is "#32770", so consider any window with this class,
    // a non-empty title and WS_POPUP or WS_OVERLAPPED to be a blocking popup window.
    if( text == "" || windowClass != "#32770" || !(style & WS_POPUP || style & WS_OVERLAPPED) || (style & WS_CHILD) )
      continue;

    vector<HWND> childHwnds;
    try {
      win32::GetChildWindows( hwnds[i], childHwnds );
    } catch( std::runtime_error& ) {
      //LogMessage( "Window \"" + text + "\": " + e.what(), ll_warning );
      continue;
    }

    // Apply a simple heuristic to the child windows to determine whether it's a popup dialog
    if( childHwnds.size() > 20 ) // At most 20 controls
      continue;
    int buttonCount = 0, messageCount = 0, otherCount = 0;
    for( size_t j = 0; j < childHwnds.size(); ++j ) {
      if( (GetWindowLong( childHwnds[j], GWL_STYLE ) & WS_VISIBLE) ) {
        string childClass = win32::ffGetWindowClassName( childHwnds[j] );
        if( childClass == "Button" )
          buttonCount++;
        else if( childClass == "Static" || childClass == "RICHEDIT" )
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
}*/

// Waits for the process to exit.
bool child_process::waitforexit( int timeout ) {
    // Should use waitforpid
    //  return true;
    //	return win32::WaitForProcessToTerminate( m_hProcess, timeout );
    if( !running() ) {
        if( m_childStatus == 0 ) {
            return true;
        }
        int status = -1;
        // This should reap the zombie process
        int pid = wait( &status );
        m_childStatus = status;
        if( pid == -1 ) {
            perror( "waitforexit:: The process was completed. wait error:" );
            return false;
        }
        return true;
    }
    if( timeout <= 0 ) {
        int status = -1;
        while( true ) {
            int pid = waitpid( m_pid, &status, 0 );
            if( pid == -1 ) {
                perror( "waitforexit::waitpid error" );
                return false;
            } else if( pid == m_pid ) {
                break;
            } else {
                continue;
            }
        }
        m_childStatus = status;
        return true;
    }
    using frantic::diagnostics::get_milliseconds;
    int64_t stop_at = get_milliseconds() + timeout;
    float delay = 0.01;
    bool notTimedOut = true;
    while( true ) {
        int status = -1;
        int pid = waitpid( m_pid, &status, WNOHANG );
        m_childStatus = status;
        if( pid == -1 ) {
            perror( "waitforexit::timeout::waitpid error" );
            return false;
        } else if( pid == m_pid ) {
            m_childStatus = status;
            break;
        }
        if( static_cast<int64_t>( get_milliseconds() ) > stop_at ) {
            notTimedOut = false;
            break;
        }
        sleep( delay );
    }
    return notTimedOut;
}

int child_process::getexitcode() { return WEXITSTATUS( m_childStatus ); }

bool child_process::running() {
    // cout << kill( m_pid, 0 ) << endl;
    return kill( m_pid, 0 ) != -1;
}

/*bool child_process::any_running()
{
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
}*/

bool child_process::setpriorityclass( PROCESS_PRIORITY priorityClass ) {
    //  Fork and execute a command line in Child Process
    // signal( SIGCHLD, child_signal );    //  Set up termination signal so we catch termination properly

    char priorityChar[10];
    sprintf( priorityChar, "%i", priorityClass );
    // cout << priorityClass << endl;

    char pIDchar[10];
    sprintf( pIDchar, "%i", m_pid );
    // cout << m_pid <<endl;

    // printf( executable );

    pid_t tempID;

    tempID = fork();

    if( tempID == -1 ) {
        perror( "fork error" );
    } else if( tempID == 0 ) { //  We're in the child process
        // cout << "Changing Child Process Priority." << endl;

        //  execute our target
        // cout << "Using Linux (renice).\n" << endl;
        int executeCode = execlp( "renice", "-n", (char*)priorityChar, (char*)pIDchar, (char*)0 );

        cout << executeCode << endl;
        cout << "Shouldn't get here - Child Process Priority." << endl;
        exit( 0 );
    } else {
        // printf ("At Parent after fork.\n");
    }

    return true; // +dwPriorityClass 987 -u daemon root -p 32;
}

bool child_process::fork_to( char* executeString ) {
    //  Fork and execute a command line in Child Process
    // signal( SIGCHLD, child_signal );    //  Set up termination signal so we catch termination properly

    pid_t tempID;

    tempID = fork();

    if( tempID == -1 ) {
        perror( "fork error" );
    } else if( tempID == 0 ) { //  We're in the child process
        cout << "Executing String in Child Process." << endl;
        //  execute our target
        int executeCode = execlp( executeString, executeString, (char*)0 );

        cout << executeCode;
        cout << "Shouldn't get here - Child Process Priority." << endl;
        exit( 0 );
    } else {
        // printf ("At Parent after fork.\n");
    }

    return true; // +dwPriorityClass 987 -u daemon root -p 32;
}

bool child_process::is_stdout_data() {
    if( !m_controlStdOut )
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

    if( !m_controlStdOut ) {
        perror( "child_process::getstdoutbytes: Should not call this function when not controlling stdout" );
        return false;
    }

    if( !m_binaryStdOut ) {
        perror( "child_process::getstdoutbytes: Binary stdout is disabled, use getstdoutline instead" );
        return false;
    }

    if( ExtractDataFromBuffer( bytes, bytesCount ) )
        return true;

    GetDataForBuffer( timeout );

    return ExtractDataFromBuffer( bytes, bytesCount );
}

bool child_process::getstdoutline( int timeout, std::string& outString ) {
    if( !m_controlStdOut ) {
        perror( "child_process::getstdoutbytes: Should not call this function when not controlling stdout" );
        return false;
    }

    if( m_binaryStdOut ) {
        perror( "child_process::getstdoutline: Binary stdout is enabled, use getstdoutbytes instead" );
        return false;
    }

    if( ExtractLineFromBuffer( outString ) )
        return true;

    GetDataForBuffer( timeout );

    return ExtractLineFromBuffer( outString );
}

void child_process::GetDataForBuffer( int timeout ) {
    // diagnostics::timeout_tracker theTimeout( timeout );

    // Reset the buffer start to 0
    if( m_stdoutBufferStart != 0 ) {
        memcpy( &m_stdoutBuffer[0], &m_stdoutBuffer[m_stdoutBufferStart], m_stdoutBufferCount - m_stdoutBufferStart );
        m_stdoutBufferCount -= m_stdoutBufferStart;
        m_stdoutBufferStart = 0;
    }
    if( timeout < 0 ) {
        fd_set rfds;
        FD_ZERO( &rfds );
        FD_SET( m_dataPipeCtoP[0], &rfds );

        int n = 0;
        if( FD_ISSET( m_dataPipeCtoP[0], &rfds ) ) {
            if( select( m_dataPipeCtoP[0] + 1, &rfds, NULL, NULL, NULL ) <= 0 )
                return;
            n = read( m_dataPipeCtoP[0], &m_stdoutBuffer[m_stdoutBufferCount],
                      (unsigned long)( m_stdoutBuffer.size() ) - m_stdoutBufferCount );
            if( n > 0 ) {
                m_stdoutBufferCount += n;
            }
        }
    } else {
        struct timeval tv;
        tv.tv_sec = 0;
        tv.tv_usec = timeout;

        fd_set rfds;
        FD_ZERO( &rfds );
        FD_SET( m_dataPipeCtoP[0], &rfds );

        int n = 0;
        if( FD_ISSET( m_dataPipeCtoP[0], &rfds ) ) {
            if( select( m_dataPipeCtoP[0] + 1, &rfds, NULL, NULL, &tv ) <= 0 )
                return;
            n = read( m_dataPipeCtoP[0], &m_stdoutBuffer[m_stdoutBufferCount],
                      (unsigned long)m_stdoutBuffer.size() - m_stdoutBufferCount );
            // std::cerr << "read " << n << " bytes into buffer. Starts with " << frantic::strings::to_hex_string(
            // &m_stdoutBuffer[m_stdoutBufferCount], std::min(n,12) ) << std::endl;
            if( n > 0 ) {
                m_stdoutBufferCount += n;
            }
        }
    }
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

#endif // !defined( _WIN32 )
