// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/simd/float_v.hpp>

namespace frantic {
namespace simd {
namespace FRANTIC_SIMD_NAMESPACE {

inline int_v::int_v() {}

inline int_v::int_v( int_v::value_type a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_set1_epi32( a );
#else
    m_native = a;
#endif
}

#ifdef FRANTIC_HAS_SSE2

inline int_v::int_v( value_type a, value_type b, value_type c, value_type d ) {
    m_native = _mm_set_epi32( d, c, b, a );
}

inline int_v::int_v( const native_type& a ) { m_native = a; }

#endif // #ifdef FRANTIC_HAS_SSE2

inline int_v::int_v( const float_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_cvtps_epi32( a.native() );
#else
    m_native = static_cast<native_type>( a.native() );
#endif
}

inline int_v::value_type int_v::operator[]( std::size_t i ) const {
    assert( i < int_v::static_size );

#ifdef FRANTIC_HAS_SSE2
#if defined( _WIN32 ) && !defined( __clang__ )
    return m_native.m128i_i32[i];
#else
    union {
        __m128i vector;
        value_type scalar[4];
    } u;
    u.vector = m_native;
    return u.scalar[i];
#endif
#else
    (void)i; // suppress "unreferenced formal parameter"
    return m_native;
#endif
}

inline const int_v::native_type& int_v::native() const { return m_native; }

inline int_v& int_v::operator&=( const int_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_and_si128( m_native, a.m_native );
#else
    m_native &= a.native();
#endif
    return *this;
}

inline int_v& int_v::and_not( const int_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_andnot_si128( a.native(), m_native );
#else
    m_native &= ~a.native();
#endif
    return *this;
}

inline int_v int_v::reinterpret( const float_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_castps_si128( a.native() );
#else
    union {
        float f;
        boost::int32_t i;
    } u;
    u.f = a.native();
    return u.i;
#endif
}

inline int_v int_v::select( const int_v& a, const int_v& b, const int_v& mask ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_or_si128( _mm_andnot_si128( mask.m_native, a.m_native ), _mm_and_si128( mask.m_native, b.m_native ) );
#else
    return ( ~mask.native() & a.native() ) | ( mask.native() & b.native() );
#endif
}

inline int_v operator&( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return int_v( _mm_and_si128( a.native(), b.native() ) );
#else
    return a.native() & b.native();
#endif
}

inline int_v operator+( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return int_v( _mm_add_epi32( a.native(), b.native() ) );
#else
    return a.native() + b.native();
#endif
}

inline int_v operator-( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return int_v( _mm_sub_epi32( a.native(), b.native() ) );
#else
    return a.native() - b.native();
#endif
}

inline int_v operator<( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return int_v( _mm_cmplt_epi32( a.native(), b.native() ) );
#else
    return a.native() < b.native() ? -1 : 0;
#endif
}

inline int_v operator>( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return int_v( _mm_cmpgt_epi32( a.native(), b.native() ) );
#else
    return a.native() > b.native() ? -1 : 0;
#endif
}

inline int_v operator<=( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_or_si128( _mm_cmpeq_epi32( a.native(), b.native() ), _mm_cmplt_epi32( a.native(), b.native() ) );
#else
    return a.native() <= b.native() ? -1 : 0;
#endif
}

inline int_v operator>=( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_or_si128( _mm_cmpeq_epi32( a.native(), b.native() ), _mm_cmpgt_epi32( a.native(), b.native() ) );
#else
    return a.native() >= b.native() ? -1 : 0;
#endif
}

inline bool operator==( const int_v& a, const int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_movemask_epi8( _mm_cmpeq_epi32( a.native(), b.native() ) ) == 0xFFFF;
#else
    return a.native() == b.native();
#endif
}

template <typename CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const int_v& v ) {
    out << "[";
    if( int_v::static_size > 0 ) {
        out << v[0];
    }
    for( int_v::size_type i = 1; i < int_v::static_size; ++i ) {
        out << ", " << v[i];
    }
    out << "]";
    return out;
}

} // namespace FRANTIC_SIMD_NAMESPACE
} // namespace simd
} // namespace frantic

namespace std {

inline frantic::simd::int_v max( const frantic::simd::int_v& a, const frantic::simd::int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    __m128i mask = _mm_cmpgt_epi32( a.native(), b.native() );
    return _mm_or_si128( _mm_and_si128( mask, a.native() ), _mm_andnot_si128( mask, b.native() ) );
    // If you have SSE4.1:
    // return int_v( _mm_max_epi32( a.m_native, b.m_native ) );
#else
    return std::max( a.native(), b.native() );
#endif
}

inline frantic::simd::int_v min( const frantic::simd::int_v& a, const frantic::simd::int_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    __m128i mask = _mm_cmplt_epi32( a.native(), b.native() );
    return _mm_or_si128( _mm_and_si128( mask, a.native() ), _mm_andnot_si128( mask, b.native() ) );
    // If you have SSE4.1:
    // return int_v( _mm_min_epi32( a.m_native, b.m_native ) );
#else
    return std::min( a.native(), b.native() );
#endif
}

} // namespace std
