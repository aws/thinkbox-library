// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/logging/progress_logger.hpp>

/**
 * This macro provides an efficient way to log global progress. By making a call such as:
 * FF_GLOBAL_PROGRESS().update_progress( 123, 999 );
 * This will only call update_progress(long,long) if global progress is actually active. Compared to:
 * get_global_progress_logger().update_progress( 123, 999 );
 * which has to call a virtual function even if logging is disabled (albeit to a null_progress_logger which just
 * returns.)
 */
#define FF_GLOBAL_PROGRESS()                                                                                           \
    if( frantic::logging::is_global_progress_active() )                                                                \
    frantic::logging::detail::get_global_progress_impl()

namespace frantic {
namespace logging {

namespace detail {
// This global object is exposed for returning from get_global_progress_logger() when progress is disabled.
extern frantic::logging::null_progress_logger g_nullProgress;

// This global bool is exposed for use by is_global_progress_active()
extern bool g_isGlobalProgressActive;

/**
 * Returns a progress_logger object for use by clients. The result of this is undefined if g_isGlobalProgressActive is
 * false
 * @return A reference to a thread-safe global progress logger.
 */
frantic::logging::progress_logger& get_global_progress_impl();

void set_global_progress_impl( frantic::logging::progress_logger& progressLogger, bool isThreadSafe );
void reset_global_progress_impl();
} // namespace detail

inline bool is_global_progress_active() { return detail::g_isGlobalProgressActive; }

inline frantic::logging::progress_logger& get_global_progress_logger() {
    if( is_global_progress_active() )
        return detail::get_global_progress_impl();
    return detail::g_nullProgress;
}

inline void check_global_abort() { FF_GLOBAL_PROGRESS().check_for_abort(); }

class scoped_global_progress_section {
  public:
    scoped_global_progress_section( frantic::logging::progress_logger& progressLogger, bool isThreadSafe = false ) {
        detail::set_global_progress_impl( progressLogger, isThreadSafe );
    }

    ~scoped_global_progress_section() { detail::reset_global_progress_impl(); }
};

} // namespace logging
} // namespace frantic
