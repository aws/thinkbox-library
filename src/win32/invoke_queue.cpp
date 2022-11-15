// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( _WIN32 )

#include <frantic/win32/invoke_queue.hpp>

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>

#define WM_INVOKE ( WM_USER + 53 )
#define WM_INVOKE_FANCY ( WM_USER + 54 )

namespace frantic {
namespace win32 {

LRESULT invoke_queue::WndProc( HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam ) {
    switch( uMsg ) {
    case WM_INVOKE:
        reinterpret_cast<void ( * )( void* )>( lParam )( reinterpret_cast<void*>( wParam ) );
        break;
    case WM_INVOKE_FANCY: {
        std::unique_ptr<const boost::function<void( void )>> pFn(
            reinterpret_cast<const boost::function<void( void )>*>( lParam ) );

        ( *pFn )();
    } break;
    default:
        return DefWindowProc( hwnd, uMsg, wParam, lParam );
    }

    return 0;
}

invoke_queue::invoke_queue()
    : m_hQueueWnd( NULL )
    , m_hModuleInst( NULL ) {}

invoke_queue::~invoke_queue() {
    if( m_hQueueWnd != NULL )
        DestroyWindow( m_hQueueWnd );

    m_hQueueWnd = NULL;

    // TODO: We end up not having the Window Class unregistered... Is that a problem?
}

void invoke_queue::initialize( HINSTANCE hModuleInstance ) {
    assert( !m_hQueueWnd ); // Don't call initialize more than once.

    m_hModuleInst = hModuleInstance;

    WNDCLASSEX wndClass;

    // Create the window class if it isn't created already.
    if( !GetClassInfoEx( m_hModuleInst, TEXT( "ThinkboxInvokeQueueClass" ), &wndClass ) ) {
        ZeroMemory( &wndClass, sizeof( WNDCLASSEX ) );

        wndClass.cbSize = sizeof( WNDCLASSEX );
        wndClass.style = 0;
        wndClass.lpfnWndProc = &invoke_queue::WndProc;
        wndClass.cbClsExtra = 0;
        wndClass.cbWndExtra = 0;
        wndClass.hInstance = m_hModuleInst;
        wndClass.hIcon = NULL;
        wndClass.hCursor = NULL;
        wndClass.hbrBackground = NULL;
        wndClass.lpszMenuName = NULL;
        wndClass.lpszClassName = TEXT( "ThinkboxInvokeQueueClass" );
        wndClass.hIconSm = NULL;

        RegisterClassEx( &wndClass );
    }

    m_hQueueWnd = CreateWindow( TEXT( "ThinkboxInvokeQueueClass" ), TEXT( "ThinkboxInvokeQueue" ), 0, 0, 0, 0, 0,
                                HWND_MESSAGE, NULL, m_hModuleInst, NULL );
}

bool invoke_queue::is_owning_thread() const {
    return GetWindowThreadProcessId( m_hQueueWnd, NULL ) == GetCurrentThreadId();
}

void invoke_queue::invoke( void ( *pCallback )( void* ), void* pData ) {
    if( m_hQueueWnd == NULL )
        throw std::logic_error( "invoke_queue not initialized" );

    PostMessage( m_hQueueWnd, WM_INVOKE, reinterpret_cast<WPARAM>( pData ), reinterpret_cast<LPARAM>( pCallback ) );
}

void invoke_queue::invoke_impl( const boost::function<void( void )>& fn ) {
    if( m_hQueueWnd == NULL )
        throw std::logic_error( "invoke_queue not initialized" );

    std::unique_ptr<const boost::function<void( void )>> pFn( new boost::function<void( void )>( fn ) );

    PostMessage( m_hQueueWnd, WM_INVOKE_FANCY, NULL, reinterpret_cast<LPARAM>( pFn.get() ) );

    pFn.release();
}

} // namespace win32
} // namespace frantic

#endif // defined( _WIN32 )
