// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once


#include <string>
#include <vector>

#include "./process_priority.hpp"
#include <iostream>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

namespace frantic {
namespace process {

class child_process {
    // char m_stdoutBuffer[65536];             // = new char[4096];
    std::vector<char> m_stdoutBuffer;
    int m_stdoutBufferStart;          // Where in the buffer we currently are
    int m_stdoutBufferCount;          // How many bytes of the above buffer are valid
    bool m_stdoutBufferSkipEmptyLine; // A flag to indicate the string on the newline was already retrieved (
    bool m_createNewConsole;          // A flag indicating to pass CREATE_NEW_CONSOLE to the process creation

    //	HANDLE m_hProcess, m_hReadHandle, m_hJobObject;
    //	DWORD m_pid;

    /* storage place for the pid of the child process, and its exit status. */
    pid_t m_pid;
    int m_childStatus;

    int m_dataPipePtoC[2]; /* an array to store the file descriptors of the parent to child pipe. */
    int m_dataPipeCtoP[2]; /* an array to store the file descriptors of the child to parent pipe. */

    bool m_terminateOnExit, m_controlStdOut, m_hideWindow, m_useJobObject, m_binaryStdOut;

    // If enabled, assumes stdout will never block, and aggressively fills the buffer as much as possible
    // each time.  Used for reading Spore cmdline output.
    bool m_aggressiveBuffer;

    bool ExtractDataFromBuffer( char* bytes, int bytesCount );
    bool ExtractLineFromBuffer( std::string& outString );
    void GetDataForBuffer( int timeout );

    bool fork_to( char* executeString );

    // Don't allow copies
    child_process( const child_process& rhs );

  public:
    child_process();
    ~child_process();

    void reset();

    int Get_pID();

    char* s3;

    // Can only be called before launch() occured
    void set_terminateOnExit( bool value );
    void set_controlStdOut( bool value );
    void set_binaryStdOut( bool value );
    void set_hideWindow( bool value );
    void set_useJobObject( bool value );
    void set_createNewConsole( bool value );

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
    // int launch( std::string executable,  char *const arguments[], const char* startupdir );
    // int launch( std::string executable,  std::string arguments[], int numArgs, const char* startupdir );
    int launch( const std::string& executable, const std::string& arguments, const std::string& startupdir = "./" );
    // terminates the process
    void terminate( int exitCode );

    // sends a WM_CLOSE message to all the windows of this process
    //	void close_windows();

    //	void post_wm_message( UINT Msg, WPARAM wParam, LPARAM lParam );

    // Gets all the top level windows of the process (or if a job object is used, all the processes
    // in the job object)
    //	void get_windows( std::vector<HWND>& windows ) const;

    // Gets all the top level windows which look like popup dialogs
    //	void get_popups( std::vector<HWND>& windows ) const;

    /**
     *  Waits for the process to exit.  Returns true if it did exit, false if it timed out
     *
     * @param timeout Maximum waiting time in milliseconds
     */
    bool waitforexit( int timeout = 0 );

    int getexitcode();

    // Whether the process is active
    bool running();

    // Whether any process in the job object is active
    bool any_running();

    bool setpriorityclass( PROCESS_PRIORITY priorityClass );

    /**
     * Gets a line from the child process's stdout, into outString.  Returns true if it got a line, false if it timed
     * out.
     *
     * @param timeout Maximum waiting time in milliseconds
     * @param[out] outString The string in which the output will be stored
     */
    bool getstdoutline( int timeout, std::string& outString );

    /**
     * Gets a block of binary data from the child process's stdout
     *
     * @param timeout Maximum waiting time in milliseconds
     * @param[out] bytes Char array in which the stdout will be stored
     * @param bytesCount Number of bytes to be stored
     */
    bool getstdoutbytes( int timeout, char* bytes, int bytesCount );
    // Returns true if getstdoutline has more data to pass along
    bool is_stdout_data();
};

} // namespace process
} // namespace frantic
