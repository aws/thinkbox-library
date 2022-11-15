// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdio>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#include <frantic/diagnostics/profiling_section.hpp>

namespace frantic {
namespace files {

class pipelined_ofstream {
  private:
    FILE* m_fout;
    char* m_buffer;
    std::size_t m_bufferSize;

    volatile char* m_writeBuffer;
    volatile std::size_t m_writeBufferSize;

    volatile int m_ioError;
    volatile int m_ioSysError;

#ifdef _WIN32
    HANDLE m_threadHandle;
    HANDLE m_threadCloseEvent;
    HANDLE m_outputAvailableEvent; // Indicates that the output buffer has been filled
    HANDLE m_bufferAvailableEvent; // Indicates that the output has completed

    static unsigned long __stdcall thread_proc( void* pData );
#endif

  public:
    pipelined_ofstream();
    pipelined_ofstream( FILE* f, std::size_t bufferSize );
    ~pipelined_ofstream();

    void open( FILE* f, std::size_t bufferSize );
    void close();

    // Writes the user buffer (ie. the one returned by get_write_buffer) to file
    void write( std::size_t numBytes );

    // Returns the current buffer that is available to the user.
    char* get_write_buffer();

    // Returns the size of the buffer available to the user.
    std::size_t get_write_buffer_size() const { return m_bufferSize; }
};

} // namespace files
} // namespace frantic
