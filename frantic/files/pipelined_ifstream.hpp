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

class pipelined_ifstream {
  private:
    FILE* m_fin;
    char* m_buffer;
    std::size_t m_bufferSize;

    volatile char* m_writeBuffer;
    volatile std::size_t m_writeBufferSize;

#ifdef _WIN32
    HANDLE m_threadHandle;
    HANDLE m_threadCloseEvent;
    HANDLE m_bufferAvailableEvent;
    HANDLE m_outputAvailableEvent;

    static unsigned long __stdcall thread_proc( void* pData );
#endif

  public:
    pipelined_ifstream();
    pipelined_ifstream( FILE* f, std::size_t bufferSize );
    ~pipelined_ifstream();

    void open( FILE* f, std::size_t bufferSize );
    void close();

    std::pair<char*, std::size_t> read();
};

} // namespace files
} // namespace frantic
