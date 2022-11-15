// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/thread.hpp>
#include <frantic/logging/global_progress_logger.hpp>

namespace frantic {
namespace logging {
namespace detail {

frantic::logging::null_progress_logger g_nullProgress;
bool g_isGlobalProgressActive = false;

namespace {
frantic::logging::progress_logger* g_currentProgress = NULL;

class serializing_progress_logger : public frantic::logging::progress_logger {
    boost::thread::id m_mainThreadID;
    frantic::logging::progress_logger* m_delegateLogger;

    // TODO: Use a tbb::atomic<bool> for this, since that is more correct. However we can't tolerate being forced to
    // link to tbb, so ...
    volatile bool m_abortRequested;

  public:
    serializing_progress_logger() {}

    virtual ~serializing_progress_logger() {}

    void reset( frantic::logging::progress_logger& delegateLogger ) {
        m_delegateLogger = &delegateLogger;
        m_abortRequested = false;
        m_mainThreadID = boost::this_thread::get_id();
    }

    virtual void update_progress( float percentage ) {
        if( m_mainThreadID == boost::this_thread::get_id() ) {
            try {
                m_delegateLogger->update_progress( percentage );
            } catch( const frantic::logging::progress_cancel_exception& ) {
                m_abortRequested = true;
                throw;
            }
        } else if( m_abortRequested )
            throw frantic::logging::progress_cancel_exception( "Operation cancelled in another thread" );
    }

    virtual void update_progress( long long completed, long long maximum ) {
        if( m_mainThreadID == boost::this_thread::get_id() ) {
            try {
                m_delegateLogger->update_progress( completed, maximum );
            } catch( const frantic::logging::progress_cancel_exception& ) {
                m_abortRequested = true;
                throw;
            }
        } else if( m_abortRequested )
            throw frantic::logging::progress_cancel_exception( "Operation cancelled in another thread" );
    }

    virtual void set_title( const frantic::tstring& title ) {
        if( m_mainThreadID == boost::this_thread::get_id() )
            m_delegateLogger->set_title( title );
    }

    virtual void check_for_abort() {
        if( m_mainThreadID == boost::this_thread::get_id() ) {
            try {
                m_delegateLogger->check_for_abort();
            } catch( const frantic::logging::progress_cancel_exception& ) {
                m_abortRequested = true;
                throw;
            }
        } else if( m_abortRequested )
            throw frantic::logging::progress_cancel_exception( "Operation cancelled in another thread" );
    }
} g_serializingLogger;
} // namespace

frantic::logging::progress_logger& get_global_progress_impl() { return *g_currentProgress; }

void set_global_progress_impl( frantic::logging::progress_logger& progressLogger, bool isThreadSafe ) {
    g_isGlobalProgressActive = true;

    if( isThreadSafe ) {
        g_currentProgress = &progressLogger;
    } else {
        g_serializingLogger.reset( progressLogger );
        g_currentProgress = &g_serializingLogger;
    }
}

void reset_global_progress_impl() {
    g_isGlobalProgressActive = false;
    g_currentProgress = &g_nullProgress;
}

} // namespace detail
} // namespace logging
} // namespace frantic
