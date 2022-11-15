// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/strings/tstring.hpp>
#include <sstream>

// This is implemented as a macro because it will prevent functions from being called if they piped into a disabled
// stream. This was desired to prevent expensive operations from executing, only to have the results discarded. Usage:
//		FF_LOG(debug) << "This is a test: " << expensive_function() << std::endl;
#define FF_LOG( stream )                                                                                               \
    if( frantic::logging::get_logging_level() >= frantic::logging::level::stream && frantic::logging::stream.rdbuf() ) \
    frantic::logging::stream

namespace frantic {
namespace logging {

// 0 is no logging, 1 is errors only, 2 is warnings, 3 is with progress, 4 is with stats, 5 is debugging
enum logging_level { LOG_NONE = 0, LOG_ERRORS, LOG_WARNINGS, LOG_PROGRESS, LOG_STATS, LOG_DEBUG, LOG_CUSTOM = -1 };

// By default all the log streams write directly to stderr or stdlog. Use std::ostream::rdbuf( newStreamBuf ) to
// redirect them to other locations.
extern std::basic_ostream<frantic::tchar> debug;
extern std::basic_ostream<frantic::tchar> stats;
extern std::basic_ostream<frantic::tchar> progress;
extern std::basic_ostream<frantic::tchar> warning;
extern std::basic_ostream<frantic::tchar> error;

// This alternate enum exists to allow the macro FFLog() to operate correctly. Perhaps it should replace the other one?
namespace level {
enum { none = 0, error, warning, progress, stats, debug, custom = -1 };
}

inline void redirect_stream( std::basic_ostream<frantic::tchar>& stream,
                             std::basic_streambuf<frantic::tchar>* newBuff ) {
    stream.rdbuf( newBuff );
}

// Helper function to automagically redirect all of the streams to a specific stream buffer.
void redirect_all_streams( std::basic_streambuf<frantic::tchar>* newBuff );
void reset_default_streams();

// TODO: make this use the logging level as defined above.
void set_logging_level( int level );
int get_logging_level();
std::string get_logging_level_string();

class set_logging_level_in_scope {
  public:
    set_logging_level_in_scope( int level );
    ~set_logging_level_in_scope();

  private:
    int m_oldLevel;
};

// returns a string for the logging level
std::string logging_level_as_string( int level );

// a little function to easily provide a stream for a logging level without any need for knowledge
// of how its set up
std::basic_ostream<frantic::tchar>& get_logging_stream( int streamLevel );

bool is_logging_errors();
bool is_logging_warnings();
bool is_logging_progress();
bool is_logging_stats();
bool is_logging_debug();

// This functor is used in ffstreambuf to send strings one line at a time
// to another functor that handles them.
template <class Func>
class for_each_line {
    Func m_func;

  public:
    explicit for_each_line( const Func& f )
        : m_func( f ) {}
    void operator()( frantic::tchar* p ) {
        frantic::tchar* pLast = p;
        for( ; *p != '\0'; ++p ) {
            if( *p == '\n' ) {
                *p = '\0';
                m_func( pLast );
                pLast = p + 1;
            }
        }

        if( pLast != p )
            m_func( pLast );
    }
};

// This is a specialization of std::stringbuf that will send strings one line at a time
// to a functor when it is flushed. This is intended for use w/ std::ostreams that write
// to non-standard locations like the maxscript listener or a frantic::win32::log_window.
// Use std::ostream::rdbuf( newStreamBuf ) to set a stream's buffer.
// TODO: Put this in a better header file.
template <class Func>
class ffstreambuf : public std::basic_stringbuf<frantic::tchar> {
    for_each_line<Func> m_func;

  public:
    explicit ffstreambuf( const Func& f )
        : m_func( f ) {}

    int sync() {
        sputc( '\0' );
        m_func( pbase() );
        setp( pbase(), epptr() );
        return std::basic_stringbuf<frantic::tchar>::sync();
    }
};

// This helper function removes the need to explicitly set the type of the ffstreambuf's functor object.
template <class Func>
inline ffstreambuf<Func>* new_ffstreambuf( const Func& f ) {
    return new ffstreambuf<Func>( f );
}

} // namespace logging
} // namespace frantic
