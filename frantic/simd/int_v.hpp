// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

#include <frantic/simd/detect.hpp>

#ifdef FRANTIC_HAS_SSE2
#include <emmintrin.h>
#endif

namespace frantic {
namespace simd {
namespace FRANTIC_SIMD_NAMESPACE {

class float_v;

class int_v {
  public:
    typedef boost::int32_t value_type;
    typedef std::size_t size_type;

#ifdef FRANTIC_HAS_SSE2
    typedef __m128i native_type;
    enum { static_size = 4 };
#else
    typedef boost::int32_t native_type;
    enum { static_size = 1 };
#endif

    int_v();

    int_v( value_type a );

#ifdef FRANTIC_HAS_SSE2
    int_v( value_type a, value_type b, value_type c, value_type d );

    int_v( const native_type& a );
#endif // #ifdef FRANTIC_HAS_SSE2

    explicit int_v( const float_v& a );

    value_type operator[]( std::size_t i ) const;

    const native_type& native() const;

    int_v& operator&=( const int_v& a );
    int_v& and_not( const int_v& a );

    static int_v reinterpret( const float_v& a );

    static int_v select( const int_v& a, const int_v& b, const int_v& mask );

  private:
    native_type m_native;
};

int_v operator&( const int_v& a, const int_v& b );

int_v operator+( const int_v& a, const int_v& b );
int_v operator-( const int_v& a, const int_v& b );

int_v operator<( const int_v& a, const int_v& b );
int_v operator>( const int_v& a, const int_v& b );
int_v operator<=( const int_v& a, const int_v& b );
int_v operator>=( const int_v& a, const int_v& b );

bool operator==( const int_v& a, const int_v& b );

template <typename CharType>
std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const int_v& v );

} // namespace FRANTIC_SIMD_NAMESPACE

using namespace FRANTIC_SIMD_NAMESPACE;

} // namespace simd
} // namespace frantic

namespace std {

frantic::simd::int_v max( const frantic::simd::int_v& a, const frantic::simd::int_v& b );
frantic::simd::int_v min( const frantic::simd::int_v& a, const frantic::simd::int_v& b );

} // namespace std

#include <frantic/simd/int_v_impl.hpp>
