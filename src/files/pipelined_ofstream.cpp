// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/files/pipelined_ofstream.hpp>

#ifndef _WIN32
#define FRANTIC_DISABLE_THREADS
#endif

namespace frantic {
namespace files {

pipelined_ofstream::pipelined_ofstream()
    : m_fout( NULL )
    , m_buffer( NULL )
    , m_bufferSize( 0 )
    , m_ioError( 0 )
    , m_ioSysError( 0 ) {
#ifdef _WIN32
    m_threadHandle = 0;
    m_threadCloseEvent = 0;
    m_outputAvailableEvent = 0;
    m_bufferAvailableEvent = 0;
#endif
}

pipelined_ofstream::pipelined_ofstream( FILE* f, std::size_t bufferSize ) { open( f, bufferSize ); }

pipelined_ofstream::~pipelined_ofstream() { close(); }

#ifdef _WIN32

#pragma warning( push )
#pragma warning( disable : 4127 )
unsigned long pipelined_ofstream::thread_proc( void* pData ) {
    pipelined_ofstream& pi = *reinterpret_cast<pipelined_ofstream*>( pData );
    HANDLE hArray[2];
    hArray[0] = pi.m_outputAvailableEvent;
    hArray[1] = pi.m_threadCloseEvent;

    std::size_t numWritten;

    while( 1 ) {
        DWORD sObj = WaitForMultipleObjects( 2, hArray, FALSE, INFINITE );
        if( sObj == WAIT_OBJECT_0 ) {
            numWritten = fwrite( (void*)pi.m_writeBuffer, 1, pi.m_writeBufferSize, pi.m_fout );
            if( numWritten != pi.m_writeBufferSize || ferror( pi.m_fout ) ) {
                pi.m_ioError = errno;
                pi.m_ioSysError = _doserrno;
                return 1;
            } else
                SetEvent( pi.m_bufferAvailableEvent );
        } else if( sObj == WAIT_OBJECT_0 + 1 ) {
            return 0;
        } else {
            return 2; // This was an unexpected result so return 2 as an error
        }
    }
}
#pragma warning( pop )

#endif

// TODO: Check if the file is already open
void pipelined_ofstream::open( FILE* f, std::size_t bufferSize ) {
    m_fout = f;
    m_buffer = new char[2 * bufferSize];
    m_bufferSize = bufferSize;

    m_ioError = 0;
    m_ioSysError = 0;

    m_writeBuffer = m_buffer;
    m_writeBufferSize = 0;

#ifndef FRANTIC_DISABLE_THREADS
    m_threadCloseEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
    m_bufferAvailableEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
    m_outputAvailableEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
    if( !m_threadCloseEvent || !m_bufferAvailableEvent || !m_outputAvailableEvent )
        throw std::runtime_error( "Failed to create an I/O thread event" );

    m_threadHandle = CreateThread( NULL, 0, &pipelined_ofstream::thread_proc, this, 0, NULL );
    if( !m_threadHandle )
        throw std::runtime_error( "Failed to create I/O thread" );

    SetEvent( m_bufferAvailableEvent );
#endif
}

void pipelined_ofstream::close() {
#ifndef FRANTIC_DISABLE_THREADS
    if( m_threadHandle ) {
        SetEvent( m_threadCloseEvent );
        WaitForSingleObject( m_threadHandle, INFINITE );

        CloseHandle( m_threadHandle );
        m_threadHandle = NULL;
    }

    if( m_threadCloseEvent ) {
        CloseHandle( m_threadCloseEvent );
        m_threadCloseEvent = NULL;
    }

    if( m_bufferAvailableEvent ) {
        CloseHandle( m_bufferAvailableEvent );
        m_bufferAvailableEvent = NULL;
    }

    if( m_outputAvailableEvent ) {
        CloseHandle( m_outputAvailableEvent );
        m_outputAvailableEvent = NULL;
    }
#endif

    if( m_buffer ) {
        delete[] m_buffer;
        m_buffer = NULL;
    }
}

#ifndef FRANTIC_DISABLE_THREADS

void pipelined_ofstream::write( std::size_t numBytes ) {
    if( numBytes > m_bufferSize )
        throw std::out_of_range( "pipelined_ofstream::write() - The write buffer size supplied was too large" );

    HANDLE hArray[2];
    hArray[0] = m_threadHandle;
    hArray[1] = m_bufferAvailableEvent;

    DWORD sObj = WaitForMultipleObjects( 2, hArray, FALSE, INFINITE );
    if( sObj == WAIT_OBJECT_0 ) {
        CloseHandle( m_threadHandle );
        m_threadHandle = NULL; // The thread has exited and we don't want close() to try and touch it.

        // The thread has exited due to some error. Time to toss an exception.
        std::stringstream ss;
        ss << "pipelined_ofstream::write() The I/O thread was closed due to an error\n";
        ss << "\tError number: " << m_ioError << "\n";
        ss << "\tOS Error number: " << m_ioSysError << "\n";
        ss << "\tError message: " << strerror( m_ioError ) << std::endl;

        throw std::runtime_error( ss.str() );
    } else if( sObj == WAIT_OBJECT_0 + 1 ) {
        // There is an available buffer to swap with the one the user has marked as being filled.
        m_writeBuffer = get_write_buffer();
        m_writeBufferSize = numBytes;
        SetEvent( m_outputAvailableEvent );
    } else
        throw std::runtime_error( "pipelined_ofstream::write() WaitForMultipleObjects returned unexpected object." );
}

#else

void pipelined_ofstream::write( std::size_t numBytes ) { fwrite( (const void*)m_writeBuffer, 1, numBytes, m_fout ); }

#endif

#ifndef FRANTIC_DISABLE_THREADS

char* pipelined_ofstream::get_write_buffer() {
    return ( m_writeBuffer == m_buffer ) ? m_buffer + m_bufferSize : m_buffer;
}

#else

char* pipelined_ofstream::get_write_buffer() { return m_buffer; }

#endif

} // namespace files
} // namespace frantic
