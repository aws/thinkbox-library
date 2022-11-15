// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#define FRANTIC_DISABLE_SIMD

#include <frantic/simd/float_v.hpp>

using frantic::simd::float_v;

TEST( SIMDScalar, FloatVConstructDefault ) {
    float_v v;
    EXPECT_EQ( 1, v.static_size );
}

TEST( SIMDScalar, FloatVConstruct1 ) {
    float_v v( 1 );
    EXPECT_EQ( 1, v[0] );
}

TEST( SIMDScalar, FloatVConstructIntV ) {
    frantic::simd::int_v vi( 1 );
    float_v vf( vi );
    EXPECT_EQ( 1, vf[0] );
}

TEST( SIMDScalar, FloatVNative ) {
    float_v v( 2 );
    EXPECT_EQ( 2, v.native() );
}

TEST( SIMDScalar, FloatVSubscript ) {
    float_v v( 1 );
    EXPECT_EQ( 1, v[0] );
}

TEST( SIMDScalar, FloatVBitwiseAndAssignment ) {
    {
        float_v mask = float_v( 1 ) > float_v( 0 );
        float_v a( 2 );
        a &= mask;
        EXPECT_EQ( 2, a[0] );
    }
    {
        float_v mask = float_v( 0 ) > float_v( 1 );
        float_v a( 2 );
        a &= mask;
        EXPECT_EQ( 0, a[0] );
    }
}

TEST( SIMDScalar, FloatVPlusAssignment ) {
    float_v a( 1 );
    float_v b( 2 );
    a += b;
    EXPECT_EQ( 3, a[0] );
}

TEST( SIMDScalar, FloatVTimesAssignment ) {
    float_v a( 2 );
    float_v b( 3 );
    a *= b;
    EXPECT_EQ( 6, a[0] );
}

TEST( SIMDScalar, FloatVSum ) {
    float_v a( 2 );
    EXPECT_EQ( 2, a.sum() );
}

TEST( SIMDScalar, FloatVReinterpretIntV ) {
    frantic::simd::int_v a( 0x3f800000 );
    float_v r = float_v::reinterpret( a );
    EXPECT_EQ( float_v( 1 ), r );
}

TEST( SIMDScalar, FloatVSelectA ) {
    float_v mask = float_v( 0 ) > float_v( 1 );
    float_v a( 1 );
    float_v b( 2 );
    float_v r = float_v::select( a, b, mask );

    EXPECT_EQ( 1, r[0] );
}

TEST( SIMDScalar, FloatVSelectB ) {
    float_v mask = float_v( 1 ) > float_v( 0 );
    float_v a( 1 );
    float_v b( 2 );
    float_v r = float_v::select( a, b, mask );

    EXPECT_EQ( 2, r[0] );
}

TEST( SIMDScalar, FloatVLoad ) {
    float a[] = { 4 };
    float_v r = float_v::load( a );
    EXPECT_EQ( 4, r[0] );
}

TEST( SIMDScalar, FloatVLoad1 ) {
    float a[] = { 4 };
    float_v r = float_v::load1( a );
    EXPECT_EQ( 4, r[0] );
}

TEST( SIMDScalar, FloatVLoad1Half ) {
    half v[] = { 4 };

    float_v r = float_v::load1( v );

    EXPECT_EQ( 4, r[0] );
}

TEST( SIMDScalar, FloatVStore ) {
    float_v a = float_v( 2 );
    float r[] = { -1, -1 };
    a.store( r );

    EXPECT_EQ( 2, r[0] );
    EXPECT_EQ( -1, r[1] );
}

TEST( SIMDScalar, FloatVStore1 ) {
    float_v a = float_v( 2 );
    float r[] = { -1, -1 };
    a.store1( r );

    EXPECT_EQ( 2, r[0] );
    EXPECT_EQ( -1, r[1] );
}

TEST( SIMDScalar, FloatVStore1Half ) {
    float_v a = float_v( 2 );
    half r[] = { -1, -1 };
    a.store1( r );

    EXPECT_EQ( 2, r[0] );
    EXPECT_EQ( -1, r[1] );
}

TEST( SIMDScalar, FloatVBitwiseAnd ) {
    {
        float_v mask = float_v( 1 ) > float_v( 0 );
        float_v a( 2 );
        float_v r = a & mask;
        EXPECT_EQ( 2, r[0] );
    }
    {
        float_v mask = float_v( 1 ) > float_v( 0 );
        float_v a( 2 );
        float_v r = mask & a;
        EXPECT_EQ( 2, r[0] );
    }
    {
        float_v mask = float_v( 0 ) > float_v( 1 );
        float_v a( 2 );
        float_v r = a & mask;
        EXPECT_EQ( 0, r[0] );
    }
}

TEST( SIMDScalar, FloatVPlus ) {
    float_v a( 1 );
    float_v b( 2 );
    float_v r = a + b;
    EXPECT_EQ( 3, r[0] );
}

TEST( SIMDScalar, FloatVMinus ) {
    float_v a( 1 );
    float_v b( 2 );
    float_v r = a - b;
    EXPECT_EQ( -1, r[0] );
}

TEST( SIMDScalar, FloatVTimes ) {
    float_v a( 2 );
    float_v b( 4 );
    float_v r = a * b;
    EXPECT_EQ( 8, r[0] );
}

TEST( SIMDScalar, FloatVDivide ) {
    float_v a( 8 );
    float_v b( 2 );
    float_v r = a / b;
    EXPECT_EQ( 4, r[0] );
}

TEST( SIMDScalar, FloatVCompareLT ) {
    float_v a( 1 );
    float_v b( 2 );

    EXPECT_NE( 0, ( a < b ).native() );
    EXPECT_EQ( 0, ( b < a ).native() );
    EXPECT_EQ( 0, ( a < a ).native() );
}

TEST( SIMDScalar, FloatVCompareGT ) {
    float_v a( 1 );
    float_v b( 2 );

    EXPECT_EQ( 0, ( a > b ).native() );
    EXPECT_NE( 0, ( b > a ).native() );
    EXPECT_EQ( 0, ( a > a ).native() );
}

TEST( SIMDScalar, FloatVCompareLE ) {
    float_v a( 1 );
    float_v b( 2 );

    EXPECT_NE( 0, ( a <= b ).native() );
    EXPECT_EQ( 0, ( b <= a ).native() );
    EXPECT_NE( 0, ( a <= a ).native() );
}

TEST( SIMDScalar, FloatVCompareGE ) {
    float_v a( 1 );
    float_v b( 2 );

    EXPECT_EQ( 0, ( a >= b ).native() );
    EXPECT_NE( 0, ( b >= a ).native() );
    EXPECT_NE( 0, ( a >= a ).native() );
}

TEST( SIMDScalar, FloatVEquality ) {
    float_v a( 1 );
    float_v b( 2 );

    EXPECT_TRUE( a == a );
    EXPECT_FALSE( a == b );
}

TEST( SIMDScalar, FloatVOstreamInsert ) {
    std::stringstream ss;
    ss << float_v( 2 );
    EXPECT_EQ( "[2]", ss.str() );
}

TEST( SIMDScalar, FloatVMax ) {
    float_v a( 1 );
    float_v b( 2 );
    float_v r = std::max( a, b );
    EXPECT_EQ( 2, r[0] );
}

TEST( SIMDScalar, FloatVMin ) {
    float_v a( 1 );
    float_v b( 2 );
    float_v r = std::min( a, b );
    EXPECT_EQ( 1, r[0] );
}

TEST( SIMDScalar, FloatSqrt ) {
    float_v a( 16 );
    float_v r = std::sqrt( a );
    EXPECT_EQ( 4, r[0] );
}
