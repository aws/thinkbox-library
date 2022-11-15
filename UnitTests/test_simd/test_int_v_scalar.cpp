// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#define FRANTIC_DISABLE_SIMD

#include <frantic/simd/int_v.hpp>

using frantic::simd::float_v;
using frantic::simd::int_v;

TEST( SIMDScalar, IntVConstructDefault ) {
    int_v v;
    EXPECT_EQ( 1, v.static_size );
}

TEST( SIMDScalar, IntVConstruct1 ) {
    int_v v( 2 );
    EXPECT_EQ( 2, v[0] );
}

TEST( SIMDScalar, IntVConstructFloatV ) {
    float_v a( 4 );
    int_v r( a );
    EXPECT_EQ( 4, r[0] );
}

TEST( SIMDScalar, IntVSubscript ) {
    int_v r( 4 );
    EXPECT_EQ( 4, r[0] );
}

TEST( SIMDScalar, IntVNative ) {
    int_v r( 2 );
    EXPECT_EQ( 2, r.native() );
}

TEST( SIMDScalar, IntVBitwiseAndAssignment ) {
    int_v a( 0xFFFF0000 );
    int_v b( 0xAAAAAAAA );
    a &= b;
    EXPECT_EQ( 0xAAAA0000, a );
}

TEST( SIMDScalar, IntVBitwiseAndNotAssignment ) {
    int_v a( 0xFFFF0000 );
    int_v b( 0xAAAAAAAA );
    a.and_not( b );
    EXPECT_EQ( 0x55550000, a );
}

TEST( SIMDScalar, IntVReinterpretFloat ) {
    int_v a( 0xABCDEF12 );
    float_v b = float_v::reinterpret( a );
    int_v r = int_v::reinterpret( b );
    EXPECT_EQ( 0xABCDEF12, r.native() );
}

TEST( SIMDScalar, IntVSelectA ) {
    int_v mask = int_v( 0 ) > int_v( 1 );
    int_v a( 0xAAAAAAAA );
    int_v b( 0xBBBBBBBB );
    int_v r = int_v::select( a, b, mask );

    EXPECT_EQ( 0xAAAAAAAA, r[0] );
}

TEST( SIMDScalar, IntVSelectB ) {
    int_v mask = int_v( 1 ) > int_v( 0 );
    int_v a( 0xAAAAAAAA );
    int_v b( 0xBBBBBBBB );
    int_v r = int_v::select( a, b, mask );

    EXPECT_EQ( 0xBBBBBBBB, r[0] );
}

TEST( SIMDScalar, IntVBitwiseAnd ) {
    int_v a( 0x0000FFFF );
    int_v b( 0x55555555 );
    int_v r = a & b;
    EXPECT_EQ( 0x00005555, r[0] );
}

TEST( SIMDScalar, IntVPlus ) {
    int_v a( 1 );
    int_v b( 2 );
    int_v r = a + b;
    EXPECT_EQ( 3, r[0] );
}

TEST( SIMDScalar, IntVMinus ) {
    int_v a( 1 );
    int_v b( 2 );
    int_v r = a - b;
    EXPECT_EQ( -1, r[0] );
}

TEST( SIMDScalar, IntVCompareLT ) {
    int_v a( 1 );
    int_v b( 2 );

    EXPECT_NE( 0, ( a < b ).native() );
    EXPECT_EQ( 0, ( b < a ).native() );
    EXPECT_EQ( 0, ( a < a ).native() );
}

TEST( SIMDScalar, IntVCompareGT ) {
    int_v a( 1 );
    int_v b( 2 );

    EXPECT_EQ( 0, ( a > b ).native() );
    EXPECT_NE( 0, ( b > a ).native() );
    EXPECT_EQ( 0, ( a > a ).native() );
}

TEST( SIMDScalar, IntVCompareLE ) {
    int_v a( 1 );
    int_v b( 2 );

    EXPECT_NE( 0, ( a <= b ).native() );
    EXPECT_EQ( 0, ( b <= a ).native() );
    EXPECT_NE( 0, ( a <= a ).native() );
}

TEST( SIMDScalar, IntVCompareGE ) {
    int_v a( 1 );
    int_v b( 2 );

    EXPECT_EQ( 0, ( a >= b ).native() );
    EXPECT_NE( 0, ( b >= a ).native() );
    EXPECT_NE( 0, ( a >= a ).native() );
}

TEST( SIMDScalar, IntVEquality ) {
    int_v a( 1 );
    int_v a2( 1 );
    int_v b( 2 );

    EXPECT_TRUE( a == a2 );
    EXPECT_FALSE( a == b );
}

TEST( SIMDScalar, IntVOstreamInsert ) {
    std::stringstream ss;
    ss << int_v( 3 );
    EXPECT_EQ( "[3]", ss.str() );
}

TEST( SIMDScalar, IntVMax ) {
    int_v a( 1 );
    int_v b( 2 );
    int_v r = std::max( a, b );
    EXPECT_EQ( 2, r[0] );
}

TEST( SIMDScalar, IntVMin ) {
    int_v a( 1 );
    int_v b( 2 );
    int_v r = std::min( a, b );
    EXPECT_EQ( 1, r[0] );
}
