// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <sstream>
#include <string>

#include <frantic/diagnostics/timeout_tracker.hpp>
#include <frantic/misc/string_functions.hpp>
#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace diagnostics {

#ifdef FRANTIC_DISABLE_PROFILING

class profiling_section {
    profiling_section() {}
    profiling_section( const frantic::tstring& ) {}

    frantic::tstring name() const { return _T("(profiling disabled)"); }

    unsigned long totalTime() const { return 0; }

    void reset() {}

    void enter() {}
    void exit() {}

    frantic::tstring str() const { return _T("(profiling disabled)"); }

    template <class CharType>
    friend std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& o, const profiling_section& ps );
};

template <class CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& o, const profiling_section& ) {
    o << _T("(profiling disabled)");
    return o;
}

class scoped_profile {
  public:
    scoped_profile( profiling_section& ) {}

    void exit() {}
};

#else

/**
 * profiling_section is a class that provides timing support and data collection for profiling code. It keeps
 * track of the time within enter() and exit() function calls - works with nesting, but not multithreaded
 * use scoped_profile with this class. The class also supports outputing to standard C++ strings.
 */
class profiling_section {
    frantic::tstring m_name;
    unsigned long m_startTime, m_totalTime, m_lastTiming;
    long m_calledCount;
    bool m_nestedCalls;

    // Should probably use a counted locking primitive
    int m_insideCount;

  public:
    profiling_section() {
        m_name = _T("(unnamed)");
        m_startTime = m_totalTime = 0;
        m_insideCount = 0;
        m_calledCount = 0;
        m_lastTiming = 0;
        m_nestedCalls = false;
    }
    profiling_section( const frantic::tstring& name ) {
        m_name = name;
        m_startTime = m_totalTime = 0;
        m_insideCount = 0;
        m_calledCount = 0;
        m_lastTiming = 0;
        m_nestedCalls = false;
    }

    const frantic::tstring& name() const { return m_name; }
    void set_name( const frantic::tstring& newName ) { m_name = newName; }

    void add_manual_timing( unsigned long timeAmountMS ) {
        m_calledCount++;
        m_lastTiming = timeAmountMS;
        m_totalTime += timeAmountMS;
    }

    unsigned long last_timing() const { return m_lastTiming; }

    frantic::tstring last_timing_seconds() const { return frantic::strings::ms_to_string( m_lastTiming ); }

    unsigned long total_time() const {
        if( m_insideCount > 0 )
            return m_totalTime + ( get_milliseconds() - m_startTime );
        else
            return m_totalTime;
    }

    void reset() {
        m_startTime = m_totalTime = 0;
        m_insideCount = 0;
        m_calledCount = 0;
        m_lastTiming = 0;
        m_nestedCalls = false;
    }

    void enter() {
        if( m_insideCount == 0 )
            m_startTime = get_milliseconds();
        else
            m_nestedCalls = true;

        ++m_insideCount;
        ++m_calledCount;
    }

    void exit() {
        if( m_insideCount > 0 ) {
            --m_insideCount;

            if( m_insideCount == 0 ) {
                m_lastTiming = ( get_milliseconds() - m_startTime );
                m_totalTime += m_lastTiming;
            }
        }
    }

    frantic::tstring str() const;

    template <class CharType>
    friend std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& o, const profiling_section& ps );
};

template <class CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& o, const profiling_section& ps ) {
    using frantic::strings::ms_to_basic_string;
    o << "Section \"" << ps.name().c_str() << "\": \n";
    if( !ps.m_nestedCalls ) {
        o << "\tTotal " << ms_to_basic_string<CharType>( ps.total_time() ) << "\t Called " << ps.m_calledCount
          << " times";
        if( ps.m_calledCount > 0 )
            o << "\t Avg "
              << ms_to_basic_string<CharType>( (unsigned long)( (float)ps.total_time() / ps.m_calledCount ) );
    } else
        o << "\tTot " << ms_to_basic_string<CharType>( ps.total_time() ) << "\t Called " << ps.m_calledCount
          << " times, with nesting";
    return o;
}

inline frantic::tstring profiling_section::str() const {
    std::basic_stringstream<frantic::tchar> s;

    s << *this;

    return s.str();
}

// Provides a class to encapsulate using scope to enter and exit a profiling section
class scoped_profile {
    profiling_section& m_ps;
    bool m_didExit;

    // Don't allow copying
    scoped_profile& operator=( const scoped_profile& ) { return *this; }

  public:
    scoped_profile( profiling_section& ps )
        : m_ps( ps )
        , m_didExit( false ) {
        m_ps.enter();
    }
    ~scoped_profile() {
        if( !m_didExit )
            m_ps.exit();
    }

    void exit() {
        if( !m_didExit ) {
            m_ps.exit();
            m_didExit = true;
        }
    }
};

#endif // FRANTIC_DISABLE_PROFILING

} // namespace diagnostics
} // namespace frantic
