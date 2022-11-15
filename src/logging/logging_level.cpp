// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/logging/logging_level.hpp>

// 0 is no logging, 1 is errors only, 2 is warnings, 3 is with progress, 4 is with stats, 5 is debugging
static int g_loggingLevel = 4;

#ifdef FRANTIC_USE_WCHAR
#define FRANTIC_TCLOG ( std::wclog )
#define FRANTIC_TCERR ( std::wcerr )
#define FRANTIC_TCOUT ( std::wcout )
#else
#define FRANTIC_TCLOG ( std::clog )
#define FRANTIC_TCERR ( std::cerr )
#define FRANTIC_TCOUT ( std::cout )
#endif

namespace frantic {
namespace logging {

std::basic_ostream<frantic::tchar> debug( FRANTIC_TCLOG.rdbuf() );
std::basic_ostream<frantic::tchar> stats( FRANTIC_TCLOG.rdbuf() );
std::basic_ostream<frantic::tchar> progress( FRANTIC_TCLOG.rdbuf() );
std::basic_ostream<frantic::tchar> warning( FRANTIC_TCLOG.rdbuf() );
std::basic_ostream<frantic::tchar> error( FRANTIC_TCERR.rdbuf() );

std::string logging_level_as_string( int level ) {
    switch( level ) {
    case level::none:
        return "0 - No Logging";
    case level::error:
        return "1 - Errors";
    case level::warning:
        return "2 - Warnings";
    case level::progress:
        return "3 - Progress";
    case level::stats:
        return "4 - Stats";
    default:
        return "5 - Debug";
    }
}
std::basic_ostream<frantic::tchar>& get_logging_stream( int streamLevel ) {
    switch( streamLevel ) {
    case level::none:
        return FRANTIC_TCOUT;
    case level::error:
        return error;
    case level::warning:
        return warning;
    case level::progress:
        return progress;
    case level::stats:
        return stats;
    case level::debug:
        return debug;
    default:
        return FRANTIC_TCOUT;
    }
}

void redirect_all_streams( std::basic_streambuf<frantic::tchar>* newBuff ) {
    debug.rdbuf( newBuff );
    stats.rdbuf( newBuff );
    progress.rdbuf( newBuff );
    warning.rdbuf( newBuff );
    error.rdbuf( newBuff );
    debug.rdbuf( newBuff );
    stats.rdbuf( newBuff );
    progress.rdbuf( newBuff );
    warning.rdbuf( newBuff );
    error.rdbuf( newBuff );
}

void reset_default_streams() {
    debug.rdbuf( FRANTIC_TCLOG.rdbuf() );
    stats.rdbuf( FRANTIC_TCLOG.rdbuf() );
    progress.rdbuf( FRANTIC_TCLOG.rdbuf() );
    warning.rdbuf( FRANTIC_TCLOG.rdbuf() );
    error.rdbuf( FRANTIC_TCERR.rdbuf() );
}

// 0 is no logging, 1 is errors only, 2 is warnings, 3 is with progress, 4 is with stats, 5 is debugging
void set_logging_level( int level ) { g_loggingLevel = level; }

std::string get_logging_level_string() { return logging_level_as_string( g_loggingLevel ); }

int get_logging_level() { return g_loggingLevel; }

set_logging_level_in_scope::set_logging_level_in_scope( int level )
    : m_oldLevel( get_logging_level() ) {
    set_logging_level( level );
}

set_logging_level_in_scope::~set_logging_level_in_scope() { set_logging_level( m_oldLevel ); }

bool is_logging_errors() { return g_loggingLevel > 0; }

bool is_logging_warnings() { return g_loggingLevel > 1; }

bool is_logging_progress() { return g_loggingLevel > 2; }

bool is_logging_stats() { return g_loggingLevel > 3; }

bool is_logging_debug() { return g_loggingLevel > 4; }

} // namespace logging
} // namespace frantic
