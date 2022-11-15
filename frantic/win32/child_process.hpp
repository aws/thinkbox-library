// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <map>
#include <string>
#include <vector>
#include <windows.h>

#include <frantic/strings/tstring.hpp>

#define DESKTOP_ALL                                                                                                    \
    ( DESKTOP_READOBJECTS | DESKTOP_CREATEWINDOW | DESKTOP_CREATEMENU | DESKTOP_HOOKCONTROL | DESKTOP_JOURNALRECORD |  \
      DESKTOP_JOURNALPLAYBACK | DESKTOP_ENUMERATE | DESKTOP_WRITEOBJECTS | DESKTOP_SWITCHDESKTOP |                     \
      STANDARD_RIGHTS_REQUIRED )

#define WINSTA_ALL                                                                                                     \
    ( WINSTA_ENUMDESKTOPS | WINSTA_READATTRIBUTES | WINSTA_ACCESSCLIPBOARD | WINSTA_CREATEDESKTOP |                    \
      WINSTA_WRITEATTRIBUTES | WINSTA_ACCESSGLOBALATOMS | WINSTA_EXITWINDOWS | WINSTA_ENUMERATE | WINSTA_READSCREEN |  \
      STANDARD_RIGHTS_REQUIRED )

#define GENERIC_ACCESS ( GENERIC_READ | GENERIC_WRITE | GENERIC_EXECUTE | GENERIC_ALL )

namespace frantic {
namespace process {

class child_process {
    std::vector<char> m_stdoutBuffer; // The buffer for getting the next line from stdout
    int m_stdoutBufferStart;          // Where in the buffer we currently are
    int m_stdoutBufferCount;          // How many bytes of the above buffer are valid
    bool m_stdoutBufferSkipEmptyLine; // A flag to indicate the string on the newline was already retrieved (
    bool m_createNewConsole;          // A flag indicating to pass CREATE_NEW_CONSOLE to the process creation

    HANDLE m_hProcess, m_hReadHandle, m_hWriteHandle, m_hJobObject;
    DWORD m_pid;
    bool m_terminateOnExit, m_controlStdOut, m_controlStdIn, m_hideWindow, m_useJobObject, m_binaryStdOut;

    // If enabled, assumes stdout will never block, and aggressively fills the buffer as much as possible
    // each time.  Used for reading Spore cmdline output.
    bool m_aggressiveBuffer;

    bool ExtractDataFromBuffer( char* bytes, int bytesCount );
    bool ExtractLineFromBuffer( std::string& outString );
    void GetDataForBuffer( int timeout );

    // Don't allow copies
    child_process( const child_process& rhs );

  public:
    child_process( bool hideWindow = false, bool terminateOnExit = false, bool controlStdOut = false );
    ~child_process();

    void reset();

    // Can only be called before launch() occured
    void set_terminateOnExit( bool value );
    void set_controlStdOut( bool value );
    void set_controlStdIn( bool value );
    void set_binaryStdOut( bool value );
    void set_hideWindow( bool value );
    void set_useJobObject( bool value );
    void set_createNewConsole( bool value );
    DWORD WriteToStdIn( std::string line );

    // Aggressive buffering - tries to fill the buffer as much as possible.  This hasn't been
    // checked to guarantee it won't block, because it is for use with the spore, which continually
    // generates output.
    void set_aggressiveBuffer( bool value );
    // Buffer size - A large buffer + aggressive buffering increases utilisation of the stdout connection
    // significantly.
    void set_bufferSize( int size );

    void launch_as( const std::string& executable, const std::string& arguments = "",
                    const std::string& startupdir = "", const std::string& username = "",
                    const std::string& domain = "", const std::string& password = "" );
    void launch( const frantic::tstring& executable, const frantic::tstring& arguments = _T(""),
                 const frantic::tstring& startupdir = _T("") );
    void launch_with_environment( const frantic::tstring& executable, const frantic::tstring& arguments,
                                  const frantic::tstring& startupdir,
                                  const std::map<frantic::tstring, frantic::tstring>* environment = NULL );

    // terminates the process
    void terminate( int exitCode );

    // sends a WM_CLOSE message to all the windows of this process
    void close_windows();

    void post_wm_message( UINT Msg, WPARAM wParam, LPARAM lParam );

    // Gets all the top level windows of the process (or if a job object is used, all the processes
    // in the job object)
    void get_windows( std::vector<HWND>& windows ) const;

    // Gets all the top level windows which look like popup dialogs
    void get_popups( std::vector<HWND>& windows ) const;

    // Waits for the process to exit.  Returns true if it did exit, false if it timed out
    bool waitforexit( int timeout = 0 );

    int getexitcode();

    // Whether the process is active
    bool running();

    // Whether any process in the job object is active
    bool any_running();

    void get_process_pids( std::vector<DWORD>& pids );

    bool setpriorityclass( DWORD dwPriorityClass );
    int getpriorityclass();

    bool getprocessaffinitymask( DWORD_PTR& dwProcessAffinityMask );
    bool setprocessaffinitymask( DWORD_PTR dwProcessAffinityMask );

    /**
     * Gets a line from the child process's stdout, into outString
     *
     * @param timeout The time (in milliseconds) it should wait for. Negative number means it will wait for infinity.
     * @param[out] outString string to get the line into.
     *
     * @return true if it got a line, false if it timed out.
     */
    bool getstdoutline( int timeout, std::string& outString );
    // Gets a block of binary data from the child process's stdout
    bool getstdoutbytes( int timeout, char* bytes, int bytesCount );
    // Returns true if getstdoutline has more data to pass along
    bool is_stdout_data();

    // Gets all the processes in the job object (or just the one pid if no job object is uesd)
    template <class InsertIter>
    void get_processes( InsertIter pidIns ) const {
        if( m_pid != 0 ) {
#ifndef _WIN64
            if( m_useJobObject ) {
                int bufferSize = 1024;
                char* buffer = new char[bufferSize];
                JOBOBJECT_BASIC_PROCESS_ID_LIST* plist = (JOBOBJECT_BASIC_PROCESS_ID_LIST*)buffer;

                if( !QueryInformationJobObject( m_hJobObject, JobObjectBasicProcessIdList, (void*)plist, bufferSize,
                                                NULL ) ) {
                    delete[] buffer;
                    throw runtime_error( "child_process::get_processes: Error in QueryInformationJobObject: " +
                                         win32::GetLastErrorMessageA() );
                }
                // Increase the buffer until all the pids fit in the list
                while( plist->NumberOfAssignedProcesses != plist->NumberOfProcessIdsInList ) {
                    delete[] buffer;
                    bufferSize *= 2;
                    buffer = new char[bufferSize];
                    plist = (JOBOBJECT_BASIC_PROCESS_ID_LIST*)buffer;
                    if( !QueryInformationJobObject( m_hJobObject, JobObjectBasicProcessIdList, (void*)plist, bufferSize,
                                                    NULL ) ) {
                        delete[] buffer;
                        throw runtime_error( "child_process::get_processes: Error in QueryInformationJobObject: " +
                                             win32::GetLastErrorMessageA() );
                    }
                }
                for( unsigned int i = 0; i < plist->NumberOfProcessIdsInList; ++i )
                    *pidIns++ = plist->ProcessIdList[i];

                delete[] buffer;
            } else {
#endif
                *pidIns++ = m_pid;
#ifndef _WIN64
            }
#endif
        }
    }
};

} // namespace process
} // namespace frantic
