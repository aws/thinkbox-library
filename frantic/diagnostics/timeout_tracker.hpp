// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#if defined( _WIN32 ) || defined( _WIN64 )
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <boost/chrono.hpp>
#include <windows.h>
#endif

#include <string>
#if defined( __GNUC__ ) && __GNUC__ < 3
#include <limits.h>
#else
#include <limits>
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
#include <time.h>
#elif defined( __GNUC__ )
#include <sys/time.h>
#endif

#include <algorithm>

namespace frantic {
namespace diagnostics {

#if defined( _WIN32 ) || defined( _WIN64 )
// Returns a counter, calibrated in milliseconds
inline double get_milliseconds_precise() {
    LARGE_INTEGER liNow, liPrecision;
    QueryPerformanceCounter( &liNow );
    QueryPerformanceFrequency( &liPrecision );

    if( liPrecision.QuadPart == 0 ) {
        return (double)liNow.QuadPart;
    } else {
        // The below always evaluates to 1000 because a number divided by itself is ussually == 1!
        // return (double)liNow.QuadPart / liNow.QuadPart * 1000;
        // I think you mean this:
        return ( double( liNow.QuadPart ) / double( liPrecision.QuadPart ) ) * 1000.0;
    }
}
#endif

// Returns a counter, calibrated in milliseconds
inline unsigned long get_milliseconds() {
#if defined( _WIN32 ) || defined( _WIN64 )
    return (unsigned long)( boost::chrono::duration_cast<boost::chrono::milliseconds>(
                                boost::chrono::system_clock::now().time_since_epoch() )
                                .count() );
#elif defined( __APPLE__ )
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return tv.tv_sec * 1000 + ( tv.tv_usec / 1000 );
#else
    struct timespec ts;
#if __GNUC__ < 3
    clock_gettime( CLOCK_REALTIME, &ts );
#else
    clock_gettime( CLOCK_MONOTONIC, &ts );
#endif
    return ts.tv_sec * 1000 + ( ts.tv_nsec / 1000000L );
#endif
}

// Sleeps for a while
inline void sleep_milliseconds( long msec ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    Sleep( msec );
#else
    struct timeval tv;

    tv.tv_sec = 0;
    tv.tv_usec = msec * 1000;
    select( 0, NULL, NULL, NULL, &tv );
#endif
}

// Given a timeout value, it keeps track of when the timeout should expire
// If timeout is negative, the timeout is considered infinite
// Timeout is always in milliseconds.
class timeout_tracker {
    bool m_isTimeOut;
    long m_whenToTimeOut;

  public:
    timeout_tracker() {
        m_isTimeOut = false;
        m_whenToTimeOut = get_milliseconds();
    }
    timeout_tracker( int timeout ) {
        m_isTimeOut = timeout >= 0;
        m_whenToTimeOut = get_milliseconds() + timeout;
    }
    timeout_tracker( const timeout_tracker& rhs )
        : m_isTimeOut( rhs.m_isTimeOut )
        , m_whenToTimeOut( rhs.m_whenToTimeOut ) {}

    // Restarts the timeout period
    void restart_timeout( int timeout ) {
        m_isTimeOut = timeout >= 0;
        m_whenToTimeOut = get_milliseconds() + timeout;
    }

    int time_left() const {
        if( m_isTimeOut )
            return ( std::max )( 0, int( m_whenToTimeOut - get_milliseconds() ) );
        else
            return ( std::numeric_limits<int>::max )();
    }

    bool timed_out() const {
        if( m_isTimeOut )
            return ( time_left() == 0 );
        else
            return false;
    }

    // if true, then it will never time out
    bool is_infinite() const { return !m_isTimeOut; }

    void sleep( int time ) const {
        int actualSleep = time;
        if( m_isTimeOut ) {
            int timeLeft = time_left();
            if( timeLeft < actualSleep )
                actualSleep = timeLeft;
        }

        sleep_milliseconds( actualSleep );
    }

    void wait() const { sleep_milliseconds( time_left() ); }
};

} // namespace diagnostics
} // namespace frantic
