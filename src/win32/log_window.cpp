// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( _WIN32 )

#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/string_functions.hpp>
#include <frantic/win32/log_window.hpp>
#include <frantic/win32/wrap.hpp>

// include this since WIN32_LEAN_AND_MEAN prevents it otherwise
#include <Commdlg.h>

#define IDC_MAIN_EDIT 101
#define FRANTIC_LOG_WINDOW_CLASS_NAME _T("Frantic_Log_Window_Class")

using namespace std;

namespace {

BOOL ffGetTextExtentPoint32( HDC dc, const std::string& s, SIZE* pTextSize ) {
    return GetTextExtentPoint32A( dc, s.c_str(), static_cast<int>( s.size() ), pTextSize );
}

BOOL ffGetTextExtentPoint32( HDC dc, const std::wstring& s, SIZE* pTextSize ) {
    return GetTextExtentPoint32W( dc, s.c_str(), static_cast<int>( s.size() ), pTextSize );
}

template <class CharType>
int get_string_width_in_pixels( HWND hEdit, const std::basic_string<CharType>& s ) {
    HDC hDC = GetDC( hEdit );
    HFONT hFont = (HFONT)SendMessage( hEdit, WM_GETFONT, 0, 0 );
    HGDIOBJ hOldFont = SelectObject( hDC, hFont );

    SIZE textSize;
    BOOL success = ffGetTextExtentPoint32( hDC, s, &textSize );

    SelectObject( hDC, hOldFont );
    ReleaseDC( hEdit, hDC );

    if( success ) {
        return textSize.cx;
    } else {
        return 0;
    }
}

} // anonymous namespace

namespace frantic {
namespace win32 {
// The ids for menu buttons so we can know which one fired.
enum Events {
    ID_FILE = 4000,
    ID_FILE_SAVEAS,
    ID_FILE_CLOSE,
    ID_EDIT = 4010,
    ID_EDIT_COPY,
    ID_EDIT_CLEAR,
    ID_EDIT_LOGGING = 4020,
    ID_EDIT_LOGGING_NONE,
    ID_EDIT_LOGGING_ERROR,
    ID_EDIT_LOGGING_WARNING,
    ID_EDIT_LOGGING_PROGRESS,
    ID_EDIT_LOGGING_STATS,
    ID_EDIT_LOGGING_DEBUG
};

// Static member that handles registering and unregistering the window class.
int log_window::m_ref_count = 0;

LRESULT CALLBACK WndProcedure( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam ) {
    HWND hEdit;
    switch( uMsg ) {
    case WM_SIZE:
        RECT rcClient;

        GetClientRect( hWnd, &rcClient );

        hEdit = GetDlgItem( hWnd, IDC_MAIN_EDIT );
        SetWindowPos( hEdit, NULL, 0, 0, rcClient.right, rcClient.bottom, SWP_NOZORDER );
        break;
    case WM_COMMAND:
        if( HIWORD( wParam ) == 0 && hWnd ) { // Check if its a menu-button notification
            log_window* pLogWnd = reinterpret_cast<log_window*>( GetWindowLongPtr( hWnd, GWLP_USERDATA ) );
            pLogWnd->fire_menu_button( LOWORD( wParam ), uMsg );
        } else {
            return DefWindowProc( hWnd, uMsg, wParam, lParam );
        }
        break;
    case WM_ENTERMENULOOP: {
        const int currentLogID = ID_EDIT_LOGGING_NONE + frantic::logging::get_logging_level();
        const HMENU hMenu = GetMenu( hWnd );
        CheckMenuRadioItem( hMenu, ID_EDIT_LOGGING_NONE, ID_EDIT_LOGGING_DEBUG, currentLogID, MF_BYCOMMAND );
    }
        return DefWindowProc( hWnd, uMsg, wParam, lParam );
    case WM_CLOSE:
        ShowWindow( hWnd, SW_HIDE );
        break;
    default:
        return DefWindowProc( hWnd, uMsg, wParam, lParam );
    }
    return TRUE;
}

void log_window::fire_menu_button( UINT uID, UINT uMsg ) {
    try {
        switch( uID ) {
        case ID_FILE_SAVEAS: {

            std::vector<std::pair<frantic::tstring, frantic::tstring>> filter( 2 );
            filter[0] = std::make_pair( _T("Log Files (*.log)"), _T("*.log") );
            filter[1] = std::make_pair( _T("All Files (*.*)"), _T("*.*") );
            frantic::tstring filename =
                frantic::win32::ffSaveDialog( m_hWnd, _T("Save As"), _T(""), _T("log"), filter );

            if( filename != _T("") ) {
                std::basic_ofstream<TCHAR> fout;
                fout.open( filename.c_str(), std::ios::out );
                get_all_text( fout );
                fout << _T("\n");
                fout.close();
            }

            break;
        }
        case ID_FILE_CLOSE:
            ShowWindow( m_hWnd, SW_HIDE );
            break;
        case ID_EDIT_COPY: {
            frantic::tstring theText = get_selected_text(); // Get the text
            if( OpenClipboard( m_hWnd ) ) {
                EmptyClipboard();

                if( theText.size() == 0 )
                    break;

                HGLOBAL hText = GlobalAlloc( GMEM_MOVEABLE, theText.size() * sizeof( frantic::tstring::value_type ) );

                TCHAR* clipboardData = (TCHAR*)GlobalLock( hText );
                memcpy( clipboardData, theText.c_str(), theText.size() * sizeof( frantic::tstring::value_type ) );
                clipboardData[theText.size() - 2] = TCHAR( '\0' );
                GlobalUnlock( hText );

#ifdef _UNICODE
                const UINT clipboardFormat = CF_UNICODETEXT;
#else
                const UINT clipboardFormat = CF_TEXT;
#endif
                SetClipboardData( clipboardFormat, hText );
                CloseClipboard();
            }
            break;
        }
        case ID_EDIT_CLEAR:
            set_horizontal_extent( 0 );
            SendMessage( m_hEdit, LB_RESETCONTENT, 0, 0 );
            break;
        case ID_EDIT_LOGGING_NONE:
        case ID_EDIT_LOGGING_ERROR:
        case ID_EDIT_LOGGING_WARNING:
        case ID_EDIT_LOGGING_PROGRESS:
        case ID_EDIT_LOGGING_STATS:
        case ID_EDIT_LOGGING_DEBUG:
            frantic::logging::set_logging_level( uID - ID_EDIT_LOGGING_NONE ); // NOTE: Requires ordered enums
            CheckMenuRadioItem( GetMenu( m_hWnd ), ID_EDIT_LOGGING_NONE, ID_EDIT_LOGGING_DEBUG, uID, MF_BYCOMMAND );
            break;
        default:
            std::stringstream ss;
            ss << "Got a message: " << uMsg << " from an unknown button: " << uID << "\n";

            MessageBoxA( m_hWnd, ss.str().c_str(), "Log Window Error!", MB_OK );
            break;
        }
    } catch( const std::exception& e ) {
        MessageBoxA( m_hWnd, e.what(), "Log Window Error!", MB_OK );
    }
}

void log_window::get_all_text( std::basic_ostream<TCHAR>& out ) {
    if( m_hWnd ) {
        int nLines = (int)SendMessage( m_hEdit, LB_GETCOUNT, NULL, NULL );

        std::vector<TCHAR> buffer;
        for( int i = 0; i < nLines; ++i ) {
            int buffSize = (int)SendMessage( m_hEdit, LB_GETTEXTLEN, i, NULL );
            buffer.resize( buffSize + 1 );

            SendMessage( m_hEdit, LB_GETTEXT, i, (LPARAM)&buffer[0] );
            out << &buffer[0] << _T("\n");
        }
    }
}

frantic::tstring log_window::get_all_text() {
    std::basic_stringstream<TCHAR> ss;
    get_all_text( ss );
    return ss.str();
}

void log_window::get_selected_text( std::basic_ostream<TCHAR>& out ) {
    if( m_hWnd ) {
        int numSelected = (int)SendMessage( m_hEdit, LB_GETSELCOUNT, NULL, NULL );
        if( numSelected <= 0 )
            return;

        std::vector<int> selectedIndices( numSelected, 0 );
        SendMessage( m_hEdit, LB_GETSELITEMS, numSelected, (LPARAM)&selectedIndices[0] );

        for( int i = 0; i < numSelected; ++i ) {
            int buffSize = (int)SendMessage( m_hEdit, LB_GETTEXTLEN, selectedIndices[i], NULL );

            std::vector<TCHAR> buff( buffSize + 1, TCHAR( '\0' ) ); // buffSize doesn't include the NULL
            SendMessage( m_hEdit, LB_GETTEXT, selectedIndices[i], (LPARAM)&buff[0] );

            out << ( &buff[0] ) << _T("\r\n");
        }
    }
}

frantic::tstring log_window::get_selected_text() {
    std::basic_stringstream<TCHAR> ss;
    get_selected_text( ss );
    return ss.str();
}

log_window::log_window( const frantic::tstring& windowName, HINSTANCE hDLLInst, HWND owner, bool deferCreation )
    : m_name( windowName )
    , m_hInst( hDLLInst )
    , m_hWnd( 0 )
    , m_hEdit( 0 )
    , m_horizontalExtent( 0 ) {
    if( !deferCreation )
        init( hDLLInst, owner );
}

log_window::~log_window() {
    if( m_hWnd ) {
        DestroyWindow( m_hWnd );
        m_hWnd = 0;
        m_hEdit = 0;
    }

    // Decrement the reference count and unregister the class if the count has gone to zero.
    if( m_ref_count > 0 ) {
        m_ref_count--;
        if( m_ref_count == 0 ) {
            UnregisterClass( FRANTIC_LOG_WINDOW_CLASS_NAME, m_hInst );
            m_hInst = 0;
        }
    }
}

void log_window::init( HINSTANCE hDLLInst, HWND hWndOwner ) {
    if( !hDLLInst ) {
        throw std::runtime_error( "Frantic Log Window Error: DLL module handle is NULL" );
    }

    if( !m_ref_count ) {
        m_hInst = hDLLInst;

        WNDCLASSEX WndClsEx;

        // Populate the WNDCLASSEX structure
        WndClsEx.cbSize = sizeof( WNDCLASSEX );
        WndClsEx.style = CS_HREDRAW | CS_VREDRAW;
        WndClsEx.lpfnWndProc = WndProcedure; // Need to define a global function with the WNDPROC interface
        WndClsEx.cbClsExtra = 0;
        WndClsEx.cbWndExtra = sizeof( LONG_PTR ); // Need to store a per-hWnd pointer to the log_window
        WndClsEx.hIcon = LoadIcon( NULL, IDI_QUESTION );
        WndClsEx.hCursor = LoadCursor( NULL, IDC_ARROW );
        WndClsEx.hbrBackground = (HBRUSH)GetStockObject( WHITE_BRUSH );
        WndClsEx.lpszMenuName = NULL;
        WndClsEx.lpszClassName = FRANTIC_LOG_WINDOW_CLASS_NAME;
        WndClsEx.hInstance = m_hInst;
        WndClsEx.hIconSm = LoadIcon( NULL, IDI_QUESTION );

        // Register the class
        RegisterClassEx( &WndClsEx );
    }
    m_ref_count++;

    if( !m_hWnd )
        create_window( hWndOwner );
}

void log_window::show( bool s ) {
    if( m_hWnd ) {
        // Show the window
        ShowWindow( m_hWnd, s ? SW_SHOWNORMAL : SW_HIDE );
        UpdateWindow( m_hWnd );
    }
}

void log_window::create_window( HWND hWndOwner ) {
    if( !m_hWnd ) {
        // Build the log window's menu, remember to delete this when the destructor is called
        HMENU theMenuBar = CreateMenu();
        HMENU editMenu, logLevelMenu;

        std::vector<std::pair<frantic::tstring, UINT>> fileButtons;
        fileButtons.push_back( std::make_pair( _T("Save &As..."), ID_FILE_SAVEAS ) );
        fileButtons.push_back( std::make_pair( _T("&Close"), ID_FILE_CLOSE ) );
        win32::AppendDropDownMenu( theMenuBar, _T("&File"), fileButtons );

        std::vector<std::pair<frantic::tstring, UINT>> editButtons;
        editButtons.push_back( std::make_pair( _T("&Copy\tCtrl+c"), ID_EDIT_COPY ) );
        editButtons.push_back( std::make_pair( _T("-"), 0 ) );
        editButtons.push_back( std::make_pair( _T("Clea&r"), ID_EDIT_CLEAR ) );
        editMenu = win32::AppendDropDownMenu( theMenuBar, _T("&Edit"), editButtons );

        std::vector<std::pair<frantic::tstring, UINT>> loggingButtons;
        loggingButtons.push_back( std::make_pair( _T("Logging &None"), ID_EDIT_LOGGING_NONE ) );
        loggingButtons.push_back( std::make_pair( _T("Logging &Errors"), ID_EDIT_LOGGING_ERROR ) );
        loggingButtons.push_back( std::make_pair( _T("Logging &Warnings"), ID_EDIT_LOGGING_WARNING ) );
        loggingButtons.push_back( std::make_pair( _T("Logging &Progress"), ID_EDIT_LOGGING_PROGRESS ) );
        loggingButtons.push_back( std::make_pair( _T("Logging &Statistics"), ID_EDIT_LOGGING_STATS ) );
        loggingButtons.push_back( std::make_pair( _T("Logging &Debug"), ID_EDIT_LOGGING_DEBUG ) );
        logLevelMenu = win32::AppendDropDownMenu( editMenu, _T("&Logging Level"), loggingButtons );

        int currentLogID = frantic::logging::get_logging_level();
        CheckMenuRadioItem( logLevelMenu, 0, GetMenuItemCount( logLevelMenu ), currentLogID, MF_BYPOSITION );

        // After the static initialization, we can create a new window.
        m_hWnd = CreateWindowEx( 0, FRANTIC_LOG_WINDOW_CLASS_NAME, m_name.c_str(), WS_OVERLAPPEDWINDOW, CW_USEDEFAULT,
                                 CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, hWndOwner, theMenuBar, m_hInst, NULL );

        // Give the window a pointer back to this log_window object.
        SetWindowLongPtr( m_hWnd, GWLP_USERDATA, (LONG_PTR)this );

        // Now add the text edit
        HFONT hfDefault;
        m_hEdit = CreateWindowEx( WS_EX_CLIENTEDGE, _T("LISTBOX"), _T(""),
                                  WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_HSCROLL | LBS_EXTENDEDSEL | LBS_USETABSTOPS,
                                  0, 0, 100, 100, m_hWnd, (HMENU)IDC_MAIN_EDIT, m_hInst, NULL );
        if( m_hEdit == NULL )
            throw std::runtime_error(
                "Frantic Log Window was unable to create message box" ); // probably could improve this message?

        hfDefault = (HFONT)GetStockObject( ANSI_FIXED_FONT );
        SendMessage( m_hEdit, WM_SETFONT, (WPARAM)hfDefault, MAKELPARAM( FALSE, 0 ) );
    }
}

void log_window::set_horizontal_extent( int widthInPixels ) {
    SendMessage( m_hEdit, LB_SETHORIZONTALEXTENT, widthInPixels, 0 );
    m_horizontalExtent = widthInPixels;
}

void log_window::union_horizontal_extent( int widthInPixels ) {
    if( widthInPixels > m_horizontalExtent ) {
        set_horizontal_extent( widthInPixels );
    }
}

void log_window::log( const std::string& msg ) {
    if( m_hWnd ) {
        UINT newItem = (UINT)SendMessageA( m_hEdit, LB_ADDSTRING, 0, (LPARAM)msg.c_str() );
        SendMessageA( m_hEdit, LB_SETTOPINDEX, newItem, NULL );
        union_horizontal_extent( get_string_width_in_pixels( m_hEdit, msg ) );
    }
}

void log_window::log( const std::wstring& msg ) {
    if( m_hWnd ) {
        UINT newItem = (UINT)SendMessageW( m_hEdit, LB_ADDSTRING, 0, (LPARAM)msg.c_str() );
        SendMessageW( m_hEdit, LB_SETTOPINDEX, newItem, NULL );
        union_horizontal_extent( get_string_width_in_pixels( m_hEdit, msg ) );
    }
}

bool log_window::is_visible() {
    if( m_hWnd ) {
        return IsWindowVisible( m_hWnd ) != 0;
    } else {
        return false;
    }
}

} // namespace win32
} // namespace frantic

#endif // defined( _WIN32 )
