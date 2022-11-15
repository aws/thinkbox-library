// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/files/pipelined_ifstream.hpp>

#ifndef _WIN32
#define FRANTIC_DISABLE_THREADS
#endif

namespace frantic {
namespace files {

namespace {
enum { ERR_CLOSED, ERR_IO, ERR_UNEXPECTED = 0xFFFF };
}

#ifndef FRANTIC_DISABLE_THREADS
namespace {
const char* err_msg[] = {
    "The thread was requested to close",
    "The thread had an I/O error",
};
}

static const char* get_err_msg( unsigned long e ) {
    if( e >= sizeof( err_msg ) / sizeof( const char* ) )
        return err_msg[e];
    return "The thread had an unexpected error";
}
#endif

pipelined_ifstream::pipelined_ifstream()
    : m_buffer( NULL ) {
#ifdef _WIN32
    m_threadHandle = NULL;
#endif
}

pipelined_ifstream::pipelined_ifstream( FILE* f, std::size_t bufferSize )
    : m_buffer( NULL ) {
#ifdef _WIN32
    m_threadHandle = NULL;
#endif
    open( f, bufferSize );
}

pipelined_ifstream::~pipelined_ifstream() { close(); }

void pipelined_ifstream::close() {
#ifndef FRANTIC_DISABLE_THREADS
    if( m_threadHandle ) {
        SetEvent( m_threadCloseEvent );
        WaitForSingleObject( m_threadHandle, INFINITE );

        CloseHandle( m_threadHandle );
        CloseHandle( m_threadCloseEvent );
        CloseHandle( m_bufferAvailableEvent );
        CloseHandle( m_outputAvailableEvent );

        m_threadHandle = NULL;
    }
#endif

    if( m_buffer ) {
        delete[] m_buffer;
        m_buffer = NULL;
    }
}

#ifndef FRANTIC_DISABLE_THREADS

#pragma warning( push )
#pragma warning( disable : 4127 )
unsigned long pipelined_ifstream::thread_proc( void* pData ) {
    pipelined_ifstream& pi = *reinterpret_cast<pipelined_ifstream*>( pData );
    HANDLE hArray[2];
    hArray[0] = pi.m_bufferAvailableEvent;
    hArray[1] = pi.m_threadCloseEvent;

    while( 1 ) {
        DWORD sObj = WaitForMultipleObjects( 2, hArray, FALSE, INFINITE );
        if( sObj == WAIT_OBJECT_0 ) {
            pi.m_writeBufferSize = std::fread( (void*)pi.m_writeBuffer, 1, pi.m_bufferSize, pi.m_fin );
            SetEvent( pi.m_outputAvailableEvent );
            if( ferror( pi.m_fin ) )
                return ERR_IO;
        } else if( sObj == WAIT_OBJECT_0 + 1 )
            return ERR_CLOSED;
        else
            return ERR_UNEXPECTED;
    }
}
#pragma warning( pop )

void pipelined_ifstream::open( FILE* f, std::size_t bufferSize ) {
    close(); // Close if already open

    m_fin = f;
    m_buffer = new char[2 * bufferSize];
    m_bufferSize = bufferSize;

    m_threadCloseEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
    m_bufferAvailableEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
    m_outputAvailableEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
    if( !m_threadCloseEvent || !m_bufferAvailableEvent || !m_outputAvailableEvent )
        throw std::runtime_error( "Failed to create an I/O thread event" );

    m_threadHandle = CreateThread( NULL, 0, &pipelined_ifstream::thread_proc, this, 0, NULL );
    if( !m_threadHandle )
        throw std::runtime_error( "Failed to create I/O thread" );

    m_writeBuffer = m_buffer;
    m_writeBufferSize = 0;

    SetEvent( m_bufferAvailableEvent ); // Get the I/O thread started
}

std::pair<char*, std::size_t> pipelined_ifstream::read() {
    char* result;
    std::size_t resultSize;

    HANDLE handleArray[2];
    handleArray[0] = m_outputAvailableEvent;
    handleArray[1] = m_threadHandle;

    DWORD r = WaitForMultipleObjects( 2, handleArray, FALSE, INFINITE );
    if( r == WAIT_OBJECT_0 ) {
        result = (char*)m_writeBuffer;
        resultSize = m_writeBufferSize;

        // Set the thread's next buffers
        m_writeBuffer = ( m_writeBuffer == m_buffer ) ? m_buffer + m_bufferSize : m_buffer;
        m_writeBufferSize = 0;
        SetEvent( m_bufferAvailableEvent );
    } else if( r == WAIT_OBJECT_0 + 1 ) {
        // The thread closed while we were expecting output ...
        DWORD threadRetCode;
        if( !GetExitCodeThread( m_threadHandle, &threadRetCode ) )
            throw std::runtime_error( "pipelined_ifstream.read: Failed to check the exit code for the i/o thread" );
        throw std::runtime_error( "pipelined_ifstream.read: " + std::string( get_err_msg( threadRetCode ) ) );
    } else
        throw std::runtime_error( "pipelined_ifstream.read: Unexpected result from WaitForMultipleObjects" );

    return std::make_pair( result, resultSize );
}

#else
// FRANTIC_DISABLE_THREADS, provides an alternate no-threaded version
#pragma message( "WARNING: FranticLibrary threads have been disabled for this build" )

// TODO: Check if the file is already open
void pipelined_ifstream::open( FILE* f, std::size_t bufferSize ) {
    close();

    m_fin = f;
    m_buffer = new char[2 * bufferSize];
    m_bufferSize = bufferSize;

    m_writeBuffer = m_buffer;
    m_writeBufferSize = 0;
}

std::pair<char*, std::size_t> pipelined_ifstream::read() {
    char* result;
    std::size_t resultSize;

    result = (char*)m_writeBuffer;
    resultSize = std::fread( result, 1, m_bufferSize, m_fin );

    // Swap the buffers next buffers
    m_writeBuffer = ( m_writeBuffer == m_buffer ) ? m_buffer + m_bufferSize : m_buffer;
    m_writeBufferSize = 0;

    return std::make_pair( result, resultSize );
}
#endif

} // namespace files
} // namespace frantic
