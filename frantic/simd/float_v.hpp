// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

#include <frantic/simd/detect.hpp>

#ifdef FRANTIC_HAS_SSE2
#include <emmintrin.h>
#endif

// Forward declarations

class half;

namespace frantic {
namespace graphics {

template <class T>
class vector3t;

typedef vector3t<float> vector3f;

} // namespace graphics
} // namespace frantic

namespace frantic {
namespace simd {
namespace FRANTIC_SIMD_NAMESPACE {

class int_v;

class float_v {
  public:
    typedef float value_type;
    typedef std::size_t size_type;

#ifdef FRANTIC_HAS_SSE2
    typedef __m128 native_type;
    enum { static_size = 4 };
#else
    typedef float native_type;
    enum { static_size = 1 };
#endif

    float_v();

    float_v( value_type a );

#ifdef FRANTIC_HAS_SSE2
    float_v( value_type a, value_type b, value_type c, value_type d );

    float_v( const native_type& a );

    explicit float_v( const frantic::graphics::vector3f& a );
#endif // #ifdef FRANTIC_HAS_SSE2

    explicit float_v( const int_v& a );

    native_type& native();
    const native_type& native() const;

    value_type operator[]( size_type i ) const;

    float_v& operator&=( const float_v& a );

    float_v& operator+=( const float_v& a );
    float_v& operator*=( const float_v& a );

    value_type sum() const;
#ifdef FRANTIC_HAS_SSE2
    value_type sum3() const;
#endif // #ifdef FRANTIC_HAS_SSE2

    static float_v reinterpret( const int_v& a );

    static float_v select( const float_v& a, const float_v& b, const float_v& mask );

    static float_v load( const value_type* a );
    static float_v load1( const value_type* a );
#ifdef FRANTIC_HAS_SSE2
    static float_v load2( const value_type* a );
    static float_v load3( const value_type* a );
    static float_v load4( const value_type* a );
#endif // #ifdef FRANTIC_HAS_SSE2

    static float_v load1( const half* a );
#ifdef FRANTIC_HAS_SSE2
    static float_v load2( const half* a );
    static float_v load3( const half* a );
    static float_v load4( const half* a );
#endif // #ifdef FRANTIC_HAS_SSE2

    void store( value_type* a ) const;
    void store1( value_type* a ) const;
#ifdef FRANTIC_HAS_SSE2
    void store2( value_type* a ) const;
    void store3( value_type* a ) const;
    void store4( value_type* a ) const;
#endif // #ifdef FRANTIC_HAS_SSE2

    void store1( half* a ) const;
#ifdef FRANTIC_HAS_SSE2
    void store2( half* a ) const;
    void store3( half* a ) const;
    void store4( half* a ) const;
#endif // #ifdef FRANTIC_HAS_SSE2

  private:
    native_type m_native;
};

float_v operator&( const float_v& a, const float_v& b );

float_v operator+( const float_v& a, const float_v& b );
float_v operator-( const float_v& a, const float_v& b );
float_v operator*( const float_v& a, const float_v& b );
float_v operator/( const float_v& a, const float_v& b );

float_v operator<( const float_v& a, const float_v& b );
float_v operator>( const float_v& a, const float_v& b );
float_v operator<=( const float_v& a, const float_v& b );
float_v operator>=( const float_v& a, const float_v& b );

bool operator==( const float_v& a, const float_v& b );

template <typename CharType>
std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const float_v& v );

} // namespace FRANTIC_SIMD_NAMESPACE

using namespace FRANTIC_SIMD_NAMESPACE;

} // namespace simd
} // namespace frantic

namespace std {

frantic::simd::float_v floor( const frantic::simd::float_v& a );
frantic::simd::float_v max( const frantic::simd::float_v& a, const frantic::simd::float_v& b );
frantic::simd::float_v min( const frantic::simd::float_v& a, const frantic::simd::float_v& b );
frantic::simd::float_v sqrt( const frantic::simd::float_v& a );

} // namespace std

#include <frantic/simd/float_v_impl.hpp>
