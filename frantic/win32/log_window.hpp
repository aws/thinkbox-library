// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <string>

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace win32 {

class log_window {
    HINSTANCE m_hInst;
    HWND m_hWnd;
    HWND m_hEdit;
    frantic::tstring m_name;
    int m_horizontalExtent; // scrollable width of the window, measured in pixels

    // from the msdn manual:
    // Windows NT/2000/XP: No window classes registered by a DLL registers are unregistered when the .dll is unloaded.
    // so, it seems we have to register the window class when the first instance is created and
    // reference count to unregister when the last instance is destroyed.
    static int m_ref_count;

  private:
    friend LRESULT CALLBACK WndProcedure( HWND, UINT, WPARAM, LPARAM );
    void fire_menu_button( UINT uID, UINT uMsg );

    // Gets text from the log window, dumping it into the provided streams
    void get_selected_text( std::basic_ostream<TCHAR>& out );
    void get_all_text( std::basic_ostream<TCHAR>& out );

    // Gets text from the log window, returning it as a string
    frantic::tstring get_selected_text();
    frantic::tstring get_all_text();

    // Window creation may be delayed from the construction of this class.  The reason for the lag is that you
    // can't switch a window's owner, and the log window needs to be owned by the main 3ds Max window.
    void create_window( HWND hWndOwner );

    // Set the scrollable width of the window
    void set_horizontal_extent( int widthInPixels );
    // Set the scrollable width of the window to the maximum of its current width or the specified widthInPixels
    void union_horizontal_extent( int widthInPixels );

  public:
    /**
     * Creates a new log_window object.  Optionally register the window class and create the window.
     * @note Previous versions of this did not take an HINSTANCE and instead used GetModuleHandle(NULL) internally. This
     *caused the window class to be registered with the process handle, not the DLL handle which means that WNDCLASS was
     *registered only once even when used in multiple DLLs. This was bad because global data (for example
     *frantic::logging::get_logging_level()) was only being used from one DLL (ie. the first to register) which was not
     *the intended design. By passing the HINSTANCE of the DLL, the class will be separately registered by each DLL.
     *
     * @param windowName The name to display in the titlebar of the log window
     * @param hDLLInst The HINSTANCE of the DLL or exe that owns this window. Should be the value passed to DLLMain(),
     *not GetModuleHandle(NULL) which is the exe's handle. If deferCreation is true, then this is ignored.  This can be
     *0 only if deferCreation is true.
     * @param owner The owner window of this log window. Will make sure the window is always on top. This is not used if
     *deferCreation is false.
     * @param deferCreation If true, the window will not be automatically created. You will need to call
     *create_window(HWND) at a later time.
     */
    log_window( const frantic::tstring& windowName, HINSTANCE hDLLInst = 0, HWND owner = 0,
                bool deferCreation = false );
    virtual ~log_window();

    /**
     *  Register the window class and create the window.
     *
     * @param hDLLInst The HINSTANCE of the DLL or exe that owns this window. Should be the value passed to DLLMain(),
     * not GetModuleHandle(NULL) which is the exe's handle.
     * @param hWndOwner The owner window of this log window. Will make sure the window is always on top.
     *
     * @note Previous versions of this did not take an HINSTANCE and instead used GetModuleHandle(NULL) internally. This
     * caused the window class to be registered with the process handle, not the DLL handle which means that WNDCLASS
     * was registered only once even when used in multiple DLLs. This was bad because global data (for example
     * frantic::logging::get_logging_level()) was only being used from one DLL (ie. the first to register) which was not
     * the intended design. By passing the HINSTANCE of the DLL, the class will be separately registered by each DLL.
     */
    void init( HINSTANCE hDLLInst, HWND hWndOwner );

    void log( const std::string& message );
    void log( const std::wstring& message );
    void show( bool s = true );
    bool is_visible();
    HWND handle() { return m_hWnd; }
};

} // namespace win32
} // namespace frantic
