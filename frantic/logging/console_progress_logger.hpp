// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/diagnostics/timeout_tracker.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace logging {

// Prints out progress information to the console.
class console_progress_logger : public progress_logger {
    frantic::diagnostics::timeout_tracker m_progressTimeout;
    std::string m_title;

    void print_progress( long long completed, long long maximum, float progress ) {
        std::cout << "\r" << m_title << "Progress: ";
        if( completed == 0 && maximum == 0 ) {
            std::cout << get_adjusted_progress( progress ) << "%                     ";
        } else {
            std::cout << completed << "/" << maximum << " (" << get_adjusted_progress( progress )
                      << "%)                     ";
        }
        std::cout.flush();
    }

  public:
    console_progress_logger( float progressStart = 0, float progressEnd = 100 )
        : progress_logger( progressStart, progressEnd )
        , m_title( "" ) {
        // print_progress( 0, 0, 0 );
        //  Print progress 5 times a second
        m_progressTimeout.restart_timeout( 200 );
    }
    virtual ~console_progress_logger() {
        if( is_logging_progress() )
            std::cout << std::endl;
    }

    // This won't become active until the next progress update
    void set_title( const frantic::tstring& title ) {
        if( title != _T("") )
            m_title = frantic::strings::to_string( title ) + " - ";
        else
            m_title = "";
    }

    void update_progress( long long completed, long long maximum ) {
        if( is_logging_progress() && ( completed == maximum || m_progressTimeout.timed_out() ) ) {
            if( maximum == 0 )
                print_progress( completed, maximum, 0 );
            else
                print_progress( completed, maximum, (float)completed * 100.f / maximum );
            // Print progress 5 times a second
            m_progressTimeout.restart_timeout( 200 );
        }
    }

    void update_progress( float progressPercent ) {
        if( is_logging_progress() && ( progressPercent >= 100 || m_progressTimeout.timed_out() ) ) {
            print_progress( 0, 0, progressPercent );
            // Print progress 5 times a second
            m_progressTimeout.restart_timeout( 200 );
        }
    }
};

} // namespace logging
} // namespace frantic
