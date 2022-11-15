// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/simd/float_v.hpp>
#include <frantic/simd/int_v.hpp>

#include <frantic/graphics/vector3f.hpp>

#ifdef FRANTIC_HAS_SSE2

using frantic::simd::float_v;

TEST( SIMD, FloatVConstructDefault ) {
    float_v v;
    EXPECT_EQ( 4, v.static_size );
}

TEST( SIMD, FloatVConstruct1 ) {
    float_v v( 1 );
    EXPECT_EQ( 4, v.static_size );
    for( std::size_t i = 0; i < v.static_size; ++i ) {
        EXPECT_EQ( 1, v[i] );
    }
}

TEST( SIMD, FloatVConstruct4 ) {
    float_v v( 1, 2, 3, 4 );

    EXPECT_EQ( 1, v[0] );
    EXPECT_EQ( 2, v[1] );
    EXPECT_EQ( 3, v[2] );
    EXPECT_EQ( 4, v[3] );
}

TEST( SIMD, FloatVConstructNative ) {
    __m128 native = _mm_setr_ps( 1, 2, 3, 4 );

    float_v v( native );

    EXPECT_EQ( 1, v[0] );
    EXPECT_EQ( 2, v[1] );
    EXPECT_EQ( 3, v[2] );
    EXPECT_EQ( 4, v[3] );
}

TEST( SIMD, FloatVConstructVector3f ) {
    frantic::graphics::vector3f v( 1, 2, 3 );
    float_v r( v );
    EXPECT_EQ( float_v( 1, 2, 3, 0 ), r );
}

TEST( SIMD, FloatVConstructIntV ) {
    frantic::simd::int_v vi( 1 );
    float_v vf( vi );
    EXPECT_EQ( 4, vf.static_size );
    for( std::size_t i = 0; i < vf.static_size; ++i ) {
        EXPECT_EQ( 1, vf[i] );
    }
}

TEST( SIMD, FloatVNative ) {
    float_v a;
    a.native() = _mm_setr_ps( 1, 2, 3, 4 );
    EXPECT_EQ( float_v( 1, 2, 3, 4 ), a );
}

TEST( SIMD, FloatVNativeConst ) {
    const float_v a( 1, 2, 3, 4 );
    float_v r( a.native() );
    EXPECT_EQ( float_v( 1, 2, 3, 4 ), r );
}

TEST( SIMD, FloatVBitwiseAndAssignment ) {
    float_v mask = float_v( 0, 0, 1, 1 ) <= float_v( 0 );
    float_v r( 1, 2, 3, 4 );
    r &= mask;
    EXPECT_EQ( float_v( 1, 2, 0, 0 ), r );
}

TEST( SIMD, FloatVPlusAssignment ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 1 );
    a += b;
    EXPECT_EQ( float_v( 1, 2, 3, 4 ), a );
}

TEST( SIMD, FloatVTimesAssignment ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 2 );
    a *= b;
    EXPECT_EQ( float_v( 0, 2, 4, 6 ), a );
}

TEST( SIMD, FloatVSum ) {
    float_v a( 1, 2, 4, 8 );

    EXPECT_EQ( 15, a.sum() );
}

TEST( SIMD, FloatVSum3 ) {
    float_v a( 1, 2, 4, 8 );

    EXPECT_EQ( 7, a.sum3() );
}

TEST( SIMD, FloatVReinterpretIntV ) {
    frantic::simd::int_v a( 0x3f800000 );
    float_v r = float_v::reinterpret( a );
    EXPECT_EQ( float_v( 1 ), r );
}

TEST( SIMD, FloatVSelect ) {
    float_v mask( float_v( 0, 1, 3, 4 ) > float_v( 2 ) );

    float_v a( 1, 2, 3, 4 );
    float_v b( 4, 5, 6, 7 );

    float_v r = float_v::select( a, b, mask );

    EXPECT_EQ( float_v( 1, 2, 6, 7 ), r );
}

TEST( SIMD, FloatVLoad ) {
    float v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load( v );

    EXPECT_EQ( float_v( 1, 2, 3, 4 ), r );
}

TEST( SIMD, FloatVLoad1 ) {
    float v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load1( v );

    EXPECT_EQ( float_v( 1, 0, 0, 0 ), r );
}

TEST( SIMD, FloatVLoad2 ) {
    float v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load2( v );

    EXPECT_EQ( float_v( 1, 2, 0, 0 ), r );
}

TEST( SIMD, FloatVLoad3 ) {
    frantic::graphics::vector3f v( 1, 2, 3 );

    float_v r = float_v::load3( &v.x );

    EXPECT_EQ( float_v( 1, 2, 3, 0 ), r );
}

TEST( SIMD, FloatVLoad4 ) {
    float v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load4( v );

    EXPECT_EQ( float_v( 1, 2, 3, 4 ), r );
}

TEST( SIMD, FloatVLoad1Half ) {
    half v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load1( v );

    EXPECT_EQ( float_v( 1, 0, 0, 0 ), r );
}

TEST( SIMD, FloatVLoad2Half ) {
    half v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load2( v );

    EXPECT_EQ( float_v( 1, 2, 0, 0 ), r );
}

TEST( SIMD, FloatVLoad3Half ) {
    half v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load3( v );

    EXPECT_EQ( float_v( 1, 2, 3, 0 ), r );
}

TEST( SIMD, FloatVLoad4Half ) {
    half v[] = { 1, 2, 3, 4 };

    float_v r = float_v::load4( v );

    EXPECT_EQ( float_v( 1, 2, 3, 4 ), r );
}

TEST( SIMD, FloatVStore ) {
    float r[] = { -1, -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( 2, r[1] );
    EXPECT_EQ( 3, r[2] );
    EXPECT_EQ( 4, r[3] );
    EXPECT_EQ( -1, r[4] );
}

TEST( SIMD, FloatVStore1 ) {
    float r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store1( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( -1, r[1] );
    EXPECT_EQ( -1, r[2] );
    EXPECT_EQ( -1, r[3] );
}

TEST( SIMD, FloatVStore2 ) {
    float r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store2( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( 2, r[1] );
    EXPECT_EQ( -1, r[2] );
    EXPECT_EQ( -1, r[3] );
}

TEST( SIMD, FloatVStore3 ) {
    float r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store3( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( 2, r[1] );
    EXPECT_EQ( 3, r[2] );
    EXPECT_EQ( -1, r[3] );
}

TEST( SIMD, FloatVStore4 ) {
    float r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store4( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( 2, r[1] );
    EXPECT_EQ( 3, r[2] );
    EXPECT_EQ( 4, r[3] );
}

TEST( SIMD, FloatVStore1Half ) {
    half r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store1( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( -1, r[1] );
    EXPECT_EQ( -1, r[2] );
    EXPECT_EQ( -1, r[3] );
}

TEST( SIMD, FloatVStore2Half ) {
    half r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store2( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( 2, r[1] );
    EXPECT_EQ( -1, r[2] );
    EXPECT_EQ( -1, r[3] );
}

TEST( SIMD, FloatVStore3Half ) {
    half r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store3( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( 2, r[1] );
    EXPECT_EQ( 3, r[2] );
    EXPECT_EQ( -1, r[3] );
}

TEST( SIMD, FloatVStore4Half ) {
    half r[] = { -1, -1, -1, -1 };

    float_v v( 1, 2, 3, 4 );
    v.store4( r );

    EXPECT_EQ( 1, r[0] );
    EXPECT_EQ( 2, r[1] );
    EXPECT_EQ( 3, r[2] );
    EXPECT_EQ( 4, r[3] );
}

TEST( SIMD, FloatVBitwiseAnd ) {
    float_v mask = float_v( 0, 0, 1, 1 ) <= float_v( 0 );
    float_v a( 1, 2, 3, 4 );
    float_v r = a & mask;
    EXPECT_EQ( float_v( 1, 2, 0, 0 ), r );
}

TEST( SIMD, FloatVPlus ) {
    float_v a( 1 );
    float_v b( 2 );
    float_v r = a + b;
    EXPECT_EQ( 4, r.static_size );
    for( std::size_t i = 0; i < r.static_size; ++i ) {
        EXPECT_EQ( 3, r[i] );
    }
}

TEST( SIMD, FloatVMinus ) {
    float_v a( 1 );
    float_v b( 2 );
    float_v r = a - b;
    EXPECT_EQ( 4, r.static_size );
    for( std::size_t i = 0; i < r.static_size; ++i ) {
        EXPECT_EQ( -1, r[i] );
    }
}

TEST( SIMD, FloatVTimes ) {
    float_v a( 2 );
    float_v b( 3 );
    float_v r = a * b;
    EXPECT_EQ( 4, r.static_size );
    for( std::size_t i = 0; i < r.static_size; ++i ) {
        EXPECT_EQ( 6, r[i] );
    }
}

TEST( SIMD, FloatDivide ) {
    float_v a( 10 );
    float_v b( 2 );
    float_v r = a / b;
    EXPECT_EQ( 4, r.static_size );
    for( std::size_t i = 0; i < r.static_size; ++i ) {
        EXPECT_EQ( 5, r[i] );
    }
}

TEST( SIMD, FloatVCompareLT ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 2 );

    float_v r = a < b;

    EXPECT_NE( 0, r[0] );
    EXPECT_NE( 0, r[1] );
    EXPECT_EQ( 0, r[2] );
    EXPECT_EQ( 0, r[3] );
}

TEST( SIMD, FloatVCompareGT ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 2 );

    float_v r = a > b;

    EXPECT_EQ( 0, r[0] );
    EXPECT_EQ( 0, r[1] );
    EXPECT_EQ( 0, r[2] );
    EXPECT_NE( 0, r[3] );
}

TEST( SIMD, FloatVCompareLE ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 2 );

    float_v r = a <= b;

    EXPECT_NE( 0, r[0] );
    EXPECT_NE( 0, r[1] );
    EXPECT_NE( 0, r[2] );
    EXPECT_EQ( 0, r[3] );
}

TEST( SIMD, FloatVCompareGE ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 2 );

    float_v r = a >= b;

    EXPECT_EQ( 0, r[0] );
    EXPECT_EQ( 0, r[1] );
    EXPECT_NE( 0, r[2] );
    EXPECT_NE( 0, r[3] );
}

TEST( SIMD, FloatVEquality ) {
    EXPECT_TRUE( float_v( 0 ) == float_v( 0 ) );
    EXPECT_TRUE( float_v( 1 ) == float_v( 1 ) );
    EXPECT_FALSE( float_v( 0, 0, 0, 1 ) == float_v( 0 ) );
    EXPECT_FALSE( float_v( 0, 0, 1, 0 ) == float_v( 0 ) );
    EXPECT_FALSE( float_v( 0, 1, 0, 0 ) == float_v( 0 ) );
    EXPECT_FALSE( float_v( 1, 0, 0, 0 ) == float_v( 0 ) );
}

TEST( SIMD, FloatVOstreamInsert ) {
    std::stringstream ss;
    ss << float_v( 0, 0.5, 1, 2 );
    EXPECT_EQ( "[0, 0.5, 1, 2]", ss.str() );
}

TEST( SIMD, FloatVFloor ) {
    float_v a( -1.1f, 0, 1, 1.1f );

    float_v r = std::floor( a );

    EXPECT_EQ( float_v( -2, 0, 1, 1 ), r );
}

TEST( SIMD, FloatVMax ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 2, 2, 2, 2 );

    float_v r = std::max( a, b );

    EXPECT_EQ( float_v( 2, 2, 2, 3 ), r );
}

TEST( SIMD, FloatVMin ) {
    float_v a( 0, 1, 2, 3 );
    float_v b( 2, 2, 2, 2 );

    float_v r = std::min( a, b );

    EXPECT_EQ( float_v( 0, 1, 2, 2 ), r );
}

TEST( SIMD, FloatVSqrt ) {
    float_v a( 16 );

    float_v r = std::sqrt( a );

    EXPECT_EQ( float_v( 4 ), r );
}

#endif
