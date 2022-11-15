// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <locale>
#include <sstream>

#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#include <boost/predef/compiler.h>

#include <half.h>

#include <frantic/graphics/vector3f.hpp>
#include <frantic/simd/int_v.hpp>

#ifdef FRANTIC_HAS_SSE2
#include <frantic/simd/sse/floor.hpp>
#endif

namespace frantic {
namespace simd {
namespace FRANTIC_SIMD_NAMESPACE {

#ifndef FRANTIC_HAS_SSE2

// Some functions to support using a bitwise mask
// TODO: change to use a mask with a single bool per float instead?

namespace detail {

inline float float_from_bits( boost::uint32_t i ) {
    union {
        float f;
        boost::uint32_t i;
    } u;
    u.i = i;
    return u.f;
}

inline boost::uint32_t bits_from_float( float f ) {
    union {
        float f;
        boost::uint32_t i;
    } u;
    u.f = f;
    return u.i;
}

inline float float_all_bits_set() { return float_from_bits( std::numeric_limits<boost::uint32_t>::max() ); }

} // namespace detail

#endif // ifndef FRANTIC_HAS_SSE2

inline float_v::float_v() {}

inline float_v::float_v( value_type a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_set1_ps( a );
#else
    m_native = a;
#endif
}

#ifdef FRANTIC_HAS_SSE2

inline float_v::float_v( value_type a, value_type b, value_type c, value_type d ) {
    m_native = _mm_set_ps( d, c, b, a );
}

inline float_v::float_v( const native_type& a ) { m_native = a; }

inline float_v::float_v( const frantic::graphics::vector3f& a ) { *this = float_v::load3( &a.x ); }

#endif // #ifdef FRANTIC_HAS_SSE2

inline float_v::float_v( const int_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_cvtepi32_ps( a.native() );
#else
    m_native = static_cast<native_type>( a.native() );
#endif
}

inline float_v::native_type& float_v::native() { return m_native; }

inline const float_v::native_type& float_v::native() const { return m_native; }

inline float_v::value_type float_v::operator[]( size_type i ) const {
    assert( i < float_v::static_size );

#ifdef FRANTIC_HAS_SSE2
#if defined( _WIN32 ) && !defined( __clang__ )
    return m_native.m128_f32[i];
#elif BOOST_COMP_GNUC && BOOST_COMP_GNUC < BOOST_VERSION_NUMBER( 5, 0, 0 )
    union {
        __m128 vector;
        value_type scalar[4];
    } u;
    u.vector = m_native;
    return u.scalar[i];
#else
    return m_native[i];
#endif
#else
    (void)i; // suppress "unreferenced formal parameter"
    return m_native;
#endif
}

inline float_v& float_v::operator&=( const float_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_and_ps( m_native, a.m_native );
#else
    *this = *this & a;
#endif
    return *this;
}

inline float_v& float_v::operator+=( const float_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_add_ps( m_native, a.m_native );
#else
    m_native += a.m_native;
#endif
    return *this;
}

inline float_v& float_v::operator*=( const float_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    m_native = _mm_mul_ps( m_native, a.m_native );
#else
    m_native *= a.m_native;
#endif
    return *this;
}

inline float_v::value_type float_v::sum() const {
#ifdef FRANTIC_HAS_SSE2
    __m128 zw = _mm_movehl_ps( m_native, m_native );
    __m128 sum = _mm_add_ps( m_native, zw );
    return _mm_cvtss_f32( _mm_add_ps( sum, _mm_shuffle_ps( sum, sum, 1 ) ) );
#else
    return m_native;
#endif
}

#ifdef FRANTIC_HAS_SSE2

inline float_v::value_type float_v::sum3() const {
    __m128 z = _mm_movehl_ps( m_native, m_native );
    __m128 y = _mm_shuffle_ps( m_native, m_native, 1 );
    return _mm_cvtss_f32( _mm_add_ps( _mm_add_ps( m_native, z ), y ) );
}

#endif // #ifdef FRANTIC_HAS_SSE2

inline float_v float_v::reinterpret( const int_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_castsi128_ps( a.native() );
#else
    return detail::float_from_bits( static_cast<boost::uint32_t>( a.native() ) );
#endif
}

inline float_v float_v::select( const float_v& a, const float_v& b, const float_v& mask ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_or_ps( _mm_andnot_ps( mask.m_native, a.m_native ), _mm_and_ps( mask.m_native, b.m_native ) );
#else
    boost::uint32_t aBits = detail::bits_from_float( a.native() );
    boost::uint32_t bBits = detail::bits_from_float( b.native() );
    boost::uint32_t maskBits = detail::bits_from_float( mask.native() );
    boost::uint32_t result = bBits & maskBits | aBits & ~maskBits;
    return detail::float_from_bits( result );
#endif
}

namespace detail {

template <std::size_t N>
struct float_v_memory;

template <>
struct float_v_memory<1> {
    static float_v load( const float_v::value_type* a ) {
        assert( a );

#ifdef FRANTIC_HAS_SSE2
        return _mm_load_ss( a );
#else
        return a[0];
#endif
    }

    static void store( float_v::value_type* dest, const float_v& src ) {
        assert( dest );

#ifdef FRANTIC_HAS_SSE2
        _mm_store_ss( dest, src.native() );
#else
        dest[0] = src.native();
#endif
    }
};

#ifdef FRANTIC_HAS_SSE2

template <>
struct float_v_memory<2> {
    static float_v load( const float_v::value_type* a ) {
        assert( a );

        return _mm_loadl_pi( _mm_setzero_ps(), reinterpret_cast<const __m64*>( a ) );
    }

    static void store( float_v::value_type* dest, const float_v& src ) {
        assert( dest );

        _mm_storel_pi( reinterpret_cast<__m64*>( dest ), src.native() );
    }
};

template <>
struct float_v_memory<3> {
    static float_v load( const float_v::value_type* a ) {
        assert( a );

        __m128 x = _mm_load_ss( a );
        __m128 y = _mm_load_ss( a + 1 );
        __m128 z = _mm_load_ss( a + 2 );
        __m128 xy = _mm_movelh_ps( x, y );
        return _mm_shuffle_ps( xy, z, _MM_SHUFFLE( 2, 0, 2, 0 ) );
    }

    static void store( float_v::value_type* dest, const float_v& src ) {
        assert( dest );

        dest[0] = src[0];
        dest[1] = src[1];
        dest[2] = src[2];
    }
};

template <>
struct float_v_memory<4> {
    static float_v load( const float_v::value_type* a ) {
        assert( a );

        return _mm_loadu_ps( a );
    }

    static void store( float_v::value_type* dest, const float_v& src ) {
        assert( dest );

        _mm_storeu_ps( dest, src.native() );
    }
};

#endif

} // namespace detail

inline float_v float_v::load( const value_type* a ) { return detail::float_v_memory<float_v::static_size>::load( a ); }

inline float_v float_v::load1( const value_type* a ) { return detail::float_v_memory<1>::load( a ); }

#ifdef FRANTIC_HAS_SSE2

inline float_v float_v::load2( const value_type* a ) { return detail::float_v_memory<2>::load( a ); }

inline float_v float_v::load3( const value_type* a ) { return detail::float_v_memory<3>::load( a ); }

inline float_v float_v::load4( const value_type* a ) { return detail::float_v_memory<4>::load( a ); }

#endif // #ifdef FRANTIC_HAS_SSE2

inline float_v float_v::load1( const half* a ) {
    assert( a );

#ifdef FRANTIC_HAS_SSE2
    return _mm_set_ss( a[0] );
#else
    return static_cast<float>( a[0] );
#endif
}

#ifdef FRANTIC_HAS_SSE2

inline float_v float_v::load2( const half* a ) {
    assert( a );
    __m128 x = _mm_set_ss( a[0] );
    __m128 y = _mm_set_ss( a[1] );
    return _mm_unpacklo_ps( x, y );
}

inline float_v float_v::load3( const half* a ) {
    assert( a );
    __m128 x = _mm_set_ss( a[0] );
    __m128 y = _mm_set_ss( a[1] );
    __m128 z = _mm_set_ss( a[2] );
    return _mm_unpacklo_ps( _mm_unpacklo_ps( x, z ), y );
}

inline float_v float_v::load4( const half* a ) {
    assert( a );
    return _mm_setr_ps( a[0], a[1], a[2], a[3] );
}

#endif // #ifdef FRANTIC_HAS_SSE2

inline void float_v::store( value_type* a ) const { detail::float_v_memory<float_v::static_size>::store( a, *this ); }

inline void float_v::store1( value_type* a ) const { detail::float_v_memory<1>::store( a, *this ); }

#ifdef FRANTIC_HAS_SSE2

inline void float_v::store2( value_type* a ) const { detail::float_v_memory<2>::store( a, *this ); }

inline void float_v::store3( value_type* a ) const { detail::float_v_memory<3>::store( a, *this ); }

inline void float_v::store4( value_type* a ) const { detail::float_v_memory<4>::store( a, *this ); }

#endif // #ifdef FRANTIC_HAS_SSE2

inline void float_v::store1( half* a ) const {
    assert( a );
    a[0] = operator[]( 0 );
}

#ifdef FRANTIC_HAS_SSE2

inline void float_v::store2( half* a ) const {
    assert( a );
    a[0] = operator[]( 0 );
    a[1] = operator[]( 1 );
}

inline void float_v::store3( half* a ) const {
    assert( a );
    a[0] = operator[]( 0 );
    a[1] = operator[]( 1 );
    a[2] = operator[]( 2 );
}

inline void float_v::store4( half* a ) const {
    assert( a );
    a[0] = operator[]( 0 );
    a[1] = operator[]( 1 );
    a[2] = operator[]( 2 );
    a[3] = operator[]( 3 );
}

#endif // #ifdef FRANTIC_HAS_SSE2

inline float_v operator&( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_and_ps( a.native(), b.native() ) );
#else
    boost::uint32_t aBits = detail::bits_from_float( a.native() );
    boost::uint32_t bBits = detail::bits_from_float( b.native() );
    boost::uint32_t result = aBits & bBits;
    return detail::float_from_bits( result );
#endif
}

inline float_v operator+( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_add_ps( a.native(), b.native() ) );
#else
    return a.native() + b.native();
#endif
}

inline float_v operator-( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_sub_ps( a.native(), b.native() ) );
#else
    return a.native() - b.native();
#endif
}

inline float_v operator*( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_mul_ps( a.native(), b.native() ) );
#else
    return a.native() * b.native();
#endif
}

inline float_v operator/( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_div_ps( a.native(), b.native() ) );
#else
    return a.native() / b.native();
#endif
}

inline float_v operator<( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_cmplt_ps( a.native(), b.native() ) );
#else
    return a.native() < b.native() ? detail::float_all_bits_set() : 0;
#endif
}

inline float_v operator>( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_cmpgt_ps( a.native(), b.native() ) );
#else
    return a.native() > b.native() ? detail::float_all_bits_set() : 0;
#endif
}

inline float_v operator<=( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_cmple_ps( a.native(), b.native() ) );
#else
    return a.native() <= b.native() ? detail::float_all_bits_set() : 0;
#endif
}

inline float_v operator>=( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return float_v( _mm_cmpge_ps( a.native(), b.native() ) );
#else
    return a.native() >= b.native() ? detail::float_all_bits_set() : 0;
#endif
}

inline bool operator==( const float_v& a, const float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return _mm_movemask_ps( _mm_cmpneq_ps( a.native(), b.native() ) ) == 0;
#else
    return a.native() == b.native();
#endif
}

template <typename CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const float_v& v ) {
    std::basic_stringstream<CharType> ss;
    // nonfinite_num_put for a standard representation of non-finite numbers,
    // for example, "inf" instead of the "1.#INF" in Visual Studio
    ss.imbue( std::locale( std::locale::classic(), new boost::math::nonfinite_num_put<CharType> ) );
    ss << "[";
    if( float_v::static_size > 0 ) {
        ss << v[0];
    }
    for( float_v::size_type i = 1; i < float_v::static_size; ++i ) {
        ss << ", " << v[i];
    }
    ss << "]";
    out << ss.str();
    return out;
}

} // namespace FRANTIC_SIMD_NAMESPACE
} // namespace simd
} // namespace frantic

namespace std {

inline frantic::simd::float_v floor( const frantic::simd::float_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    return frantic::simd::float_v( frantic::simd::sse::floor( a.native() ) );
#else
    return std::floor( a.native() );
#endif
}

inline frantic::simd::float_v max( const frantic::simd::float_v& a, const frantic::simd::float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return frantic::simd::float_v( _mm_max_ps( a.native(), b.native() ) );
#else
    return std::max( a.native(), b.native() );
#endif
}

inline frantic::simd::float_v min( const frantic::simd::float_v& a, const frantic::simd::float_v& b ) {
#ifdef FRANTIC_HAS_SSE2
    return frantic::simd::float_v( _mm_min_ps( a.native(), b.native() ) );
#else
    return std::min( a.native(), b.native() );
#endif
}

inline frantic::simd::float_v sqrt( const frantic::simd::float_v& a ) {
#ifdef FRANTIC_HAS_SSE2
    return frantic::simd::float_v( _mm_sqrt_ps( a.native() ) );
#else
    return std::sqrt( a.native() );
#endif
}

} // namespace std
