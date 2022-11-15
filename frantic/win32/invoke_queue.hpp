// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <WinDef.h>

#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/function.hpp>

namespace frantic {
namespace win32 {

/**
 * Invokes callable objects (ex. function pointers, functors, std::function<void(void)>, etc.) on the thread which owns
 * the queue. This is typically used to invoke function calls on the main thread, or a specific UI thread.
 */
class invoke_queue {
  public:
    /**
     * Creates a new, unbound queue. The "owning" thread must call initialize() before other threads call invoke.
     */
    invoke_queue();

    ~invoke_queue();

    /**
     * Initializes ownership of the queue. When other threads call invoke(), the passed callable object will be invoked
     * in the context of the thread that owns this queue. \param hModuleInstance The handle for the module hosting this
     * invoke_queue instance. \note The calling thread must implement a Windows message loop (ie. calling
     * GetMessage/PeekMessage, TranslateMessage & DispatchMessage) for this object to function.
     */
    void initialize( HINSTANCE hModuleInstance );

    /**
     * \return True if the calling thread is the owner (ie. the thread that called initialize) of the queue.
     */
    bool is_owning_thread() const;

    /**
     * Invokes the function pointer in the context of the thread that initialized this invoke_queue.
     * \param pCallback A function pointer to call. This is a pointer to a non-member function with signature 'void
     * some_func(void*)' \param pData Optional pointer to some data interpreted by by the called function.
     */
    void invoke( void ( *pCallback )( void* ), void* pData = NULL );

    /**
     * Invokes the callable object in the context of the thread that initialized this invoke_queue.
     * \tparam Callable A functor implementing 'void operator()() const'
     * \param fn The functor to invoke.
     */
    template <class Callable>
    void invoke( const Callable& fn );

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    // This is less useful than it seemed since Boost.Function doesn't support move for function objects. Gotta switch
    // to std::function.
    template <class Callable>
    void invoke( Callable&& fn );
#endif

    /**
     * Invokes the callable object in the context of the thread that initialized this invoke_queue.
     * \tparam Callable A functor implementing 'void operator()( const ParamType& p1 ) const'.
     * \param fn The functor to invoke.
     * \param p1 The parameter to pass when invoking 'fn'. A copy will be stored, so use boost::ref(p1) if you want a
     * reference stored and you can guarantee the lifetime of the referred object.
     */
    template <class Callable, class ParamType>
    void invoke( const Callable& fn, const ParamType& p1 );

    /**
     * Invokes the callable object in the context of the thread that initialized this invoke_queue.
     * \tparam Callable A functor implementing 'void operator()( const ParamType& p1, const ParamType& p2 ) const'.
     * \param fn The functor to invoke.
     * \param p1 The first parameter to pass when invoking 'fn'. A copy will be stored, so use boost::ref(p1) if you
     * want a reference stored and you can guarantee the lifetime of the referred object. \param p1 The second parameter
     * to pass when invoking 'fn'. A copy will be stored, so use boost::ref(p2) if you want a reference stored and you
     * can guarantee the lifetime of the referred object.
     */
    template <class Callable, class ParamType1, class ParamType2>
    void invoke( const Callable& fn, const ParamType1& p1, const ParamType2& p2 );

  private:
    static LRESULT CALLBACK WndProc( HWND, UINT, WPARAM, LPARAM );

    void invoke_impl( const boost::function<void( void )>& fn );

  private:
    HINSTANCE m_hModuleInst;
    HWND m_hQueueWnd;
};

template <class Callable>
inline void invoke_queue::invoke( const Callable& fn ) {
    this->invoke_impl( fn );
}

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template <class Callable>
inline void invoke_queue::invoke( Callable&& fn ) {
    this->invoke_impl( std::forward<Callable>( fn ) );
}
#endif

template <class Callable, class ParamType>
inline void invoke_queue::invoke( const Callable& fn, const ParamType& p1 ) {
    this->invoke_impl( boost::bind( fn, p1 ) );
}

template <class Callable, class ParamType1, class ParamType2>
inline void invoke_queue::invoke( const Callable& fn, const ParamType1& p1, const ParamType2& p2 ) {
    this->invoke_impl( boost::bind( fn, p1, p2 ) );
}

} // namespace win32
} // namespace frantic
