// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file
 * @author Stephen Kiazyk
 *
 * A collection of some functors I've found useful at some point or another
 */
#pragma once

#if defined( _MSC_VER )
#pragma warning( push )
// Disable warning C4800: '<type>' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning( disable : 4800 )
#endif

namespace frantic {

template <class A, class R>
class array_functor {
  private:
    array_functor& operator=( const array_functor& other );
    const A& m_array;

  public:
    array_functor( const A& a )
        : m_array( a ) {}
    array_functor( const array_functor& other )
        : m_array( other.m_array ) {}
    R operator()( size_t i ) const { return static_cast<R>( m_array[i] ); }
};

template <class A>
array_functor<A, bool> make_bool_array_functor( const A& a ) {
    return array_functor<A, bool>( a );
}

} // namespace frantic

#if defined( _MSC_VER )
#pragma warning( pop )
#endif
