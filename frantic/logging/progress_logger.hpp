// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/logging/logging_level.hpp>
#include <frantic/strings/tstring.hpp>
#include <iostream>
#include <vector>

namespace frantic {
namespace logging {

class progress_cancel_exception : public std::runtime_error {
  public:
    progress_cancel_exception( const std::string& operation )
        : std::runtime_error( "progress_cancel_exception: Canceled operation \"" + operation + "\"" ) {}
};

// This is the virtual base class for progress logging.  Programs can accept progress loggers
// which either print out progress to the console or update a progress bar through this mechanism.
class progress_logger {
  protected:
    std::vector<std::pair<float, float>> m_progressStack;
    // Percentage range for this progress display
    float m_progressStart, m_progressEnd;

    float get_adjusted_progress( float progressPercent ) {
        return ( ( 100.f - progressPercent ) * m_progressStart + progressPercent * m_progressEnd ) / 100.f;
    }

  public:
    progress_logger( float progressStart = 0, float progressEnd = 100 ) {
        m_progressStart = progressStart;
        m_progressEnd = progressEnd;
    }
    virtual ~progress_logger() {}

    void push_progress( float progressStart, float progressEnd ) {
        // Save the current progress range
        m_progressStack.push_back( std::make_pair( m_progressStart, m_progressEnd ) );
        float oldProgressStart = m_progressStart, oldProgressEnd = m_progressEnd;
        // Compute the new progress range
        m_progressStart = ( 1 - 0.01f * progressStart ) * oldProgressStart + ( 0.01f * progressStart ) * oldProgressEnd;
        m_progressEnd = ( 1 - 0.01f * progressEnd ) * oldProgressStart + ( 0.01f * progressEnd ) * oldProgressEnd;
    }

    void pop_progress() {
        if( m_progressStack.size() > 0 ) {
            m_progressStart = m_progressStack.back().first;
            m_progressEnd = m_progressStack.back().second;
            m_progressStack.pop_back();
        } else if( is_logging_warnings() ) {
            std::cerr << "WARNING: pop_progress was called when the progress stack was empty." << std::endl;
        }
    }

    void reset( float progressStart, float progressEnd ) {
        m_progressStack.clear();
        m_progressStart = progressStart;
        m_progressEnd = progressEnd;
    }

    virtual void set_title( const frantic::tstring& title ) = 0;

    virtual void update_progress( long long completed, long long maximum ) = 0;
    virtual void update_progress( float percent ) = 0;

    /**
     * Check if the operation should be aborted. It will throw a progress_cancel_exception if the operation should be
     * aborted. This is intended for use if there is no reasonable way to update the progress but you want to allow for
     * cancellation.
     */
    virtual void check_for_abort() {
        // Default implementation never aborts. Subclasses should override this to do something meaningful.
    }
};

// An RAII helper class used to ensure correct usage of the stack in the progress tracker
class progress_logger_subinterval_tracker {
    progress_logger& m_pl;

    // Disable copy constructor and asignment operator by declaring them private
    progress_logger_subinterval_tracker( const progress_logger_subinterval_tracker& src )
        : m_pl( src.m_pl ) {}
    progress_logger_subinterval_tracker& operator=( const progress_logger_subinterval_tracker& ) { return *this; }

  public:
    progress_logger_subinterval_tracker( progress_logger& pl, float progressStart, float progressEnd )
        : m_pl( pl ) {
        m_pl.push_progress( progressStart, progressEnd );
    }
    ~progress_logger_subinterval_tracker() { m_pl.pop_progress(); }
    void reset( float progressStart, float progressEnd ) {
        m_pl.pop_progress();
        m_pl.push_progress( progressStart, progressEnd );
    }
};

// A null logger for situations where logging is just not required.
class null_progress_logger : public progress_logger {
  public:
    null_progress_logger() {}

    void set_title( const frantic::tstring& /*title*/ ) {}

    void update_progress( long long /*completed*/, long long /*maximum*/ ) {}
    void update_progress( float /*percent*/ ) {}
};

} // namespace logging
} // namespace frantic
