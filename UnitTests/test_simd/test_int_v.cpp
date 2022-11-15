// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/simd/int_v.hpp>

#ifdef FRANTIC_HAS_SSE2

using frantic::simd::float_v;
using frantic::simd::int_v;

TEST( SIMD, IntVStaticSize ) {
    int_v v;
    EXPECT_EQ( 4, v.static_size );
}

TEST( SIMD, IntVConstruct1 ) {
    int_v v( 1 );
    EXPECT_EQ( 4, v.static_size );
    for( std::size_t i = 0; i < v.static_size; ++i ) {
        EXPECT_EQ( 1, v[i] );
    }
}

TEST( SIMD, IntVConstruct4 ) {
    int_v v( 1, 2, 3, 4 );

    EXPECT_EQ( 4, v.static_size );

    EXPECT_EQ( 1, v[0] );
    EXPECT_EQ( 2, v[1] );
    EXPECT_EQ( 3, v[2] );
    EXPECT_EQ( 4, v[3] );
}

TEST( SIMD, IntVConstructNative ) {
    int_v v( _mm_setr_epi32( 1, 2, 3, 4 ) );

    EXPECT_EQ( 4, v.static_size );

    EXPECT_EQ( 1, v[0] );
    EXPECT_EQ( 2, v[1] );
    EXPECT_EQ( 3, v[2] );
    EXPECT_EQ( 4, v[3] );
}

TEST( SIMD, IntVConstructFloatV ) {
    frantic::simd::float_v vf( 0, 1, 2, 3 );
    int_v vi( vf );
    EXPECT_EQ( int_v( 0, 1, 2, 3 ), vi );
}

TEST( SIMD, IntVBitwiseAndAssignment ) {
    int_v a( 0, 0, -1, -1 );
    int_v b( 0, -1, 0, -1 );

    a &= b;

    EXPECT_EQ( int_v( 0, 0, 0, -1 ), a );
}

TEST( SIMD, IntVBitwiseAndNotAssignment ) {
    int_v a( 0, 0, -1, -1 );
    int_v b( 0, -1, 0, -1 );

    a.and_not( b );

    EXPECT_EQ( int_v( 0, 0, -1, 0 ), a );
}

TEST( SIMD, IntVReinterpretFloatV ) {
    float_v mask( float_v( 1 ) > float_v( 0 ) );
    int_v r = int_v::reinterpret( mask );
    EXPECT_EQ( int_v( -1 ), r );
}

TEST( SIMD, IntVSelect ) {
    int_v mask( int_v( 0, 1, 3, 4 ) > int_v( 2 ) );

    int_v a( 1, 2, 3, 4 );
    int_v b( 4, 5, 6, 7 );

    int_v r = int_v::select( a, b, mask );

    EXPECT_EQ( int_v( 1, 2, 6, 7 ), r );
}

TEST( SIMD, IntVBitwiseAnd ) {
    int_v a( 0, 0, -1, -1 );
    int_v b( 0, -1, 0, -1 );

    int_v r = a & b;

    EXPECT_EQ( int_v( 0, 0, 0, -1 ), r );
}

TEST( SIMD, IntVPlus ) {
    int_v a( 1 );
    int_v b( 2 );

    int_v r = a + b;

    EXPECT_EQ( 3, r[0] );
}

TEST( SIMD, IntVMinus ) {
    int_v a( 1 );
    int_v b( 2 );

    int_v r = a - b;

    EXPECT_EQ( -1, r[0] );
}

TEST( SIMD, IntVCompareLT ) {
    int_v a( 0, 1, 2, 3 );
    int_v b( 2 );

    int_v r = a < b;

    EXPECT_EQ( -1, r[0] );
    EXPECT_EQ( -1, r[1] );
    EXPECT_EQ( 0, r[2] );
    EXPECT_EQ( 0, r[3] );
}

TEST( SIMD, IntVCompareGT ) {
    int_v a( 0, 1, 2, 3 );
    int_v b( 2 );

    int_v r = a > b;

    EXPECT_EQ( 0, r[0] );
    EXPECT_EQ( 0, r[1] );
    EXPECT_EQ( 0, r[2] );
    EXPECT_EQ( -1, r[3] );
}

TEST( SIMD, IntVCompareLE ) {
    int_v a( 0, 1, 2, 3 );
    int_v b( 2 );

    int_v r = a <= b;

    EXPECT_EQ( int_v( -1, -1, -1, 0 ), r );
}

TEST( SIMD, IntVCompareGE ) {
    int_v a( 0, 1, 2, 3 );
    int_v b( 2 );

    int_v r = a >= b;

    EXPECT_EQ( int_v( 0, 0, -1, -1 ), r );
}

TEST( SIMD, IntVEquality ) {
    EXPECT_TRUE( int_v( 0 ) == int_v( 0 ) );
    EXPECT_TRUE( int_v( 1 ) == int_v( 1 ) );
    EXPECT_FALSE( int_v( 0, 0, 0, 1 ) == int_v( 0 ) );
    EXPECT_FALSE( int_v( 0, 0, 1, 0 ) == int_v( 0 ) );
    EXPECT_FALSE( int_v( 0, 1, 0, 0 ) == int_v( 0 ) );
    EXPECT_FALSE( int_v( 1, 0, 0, 0 ) == int_v( 0 ) );
    EXPECT_FALSE( int_v( 0, 0, 0, 128 ) == int_v( 0 ) );
    EXPECT_FALSE( int_v( 0, 0, 0, 32768 ) == int_v( 0 ) );
}

TEST( SIMD, IntVOstreamInsert ) {
    std::stringstream ss;
    ss << int_v( 0, 1, 2, 3 );
    EXPECT_EQ( "[0, 1, 2, 3]", ss.str() );
}

TEST( SIMD, IntVMax ) {
    int_v a( 0, 1, 2, 3 );
    int_v b( 2, 2, 2, 2 );
    EXPECT_EQ( int_v( 2, 2, 2, 3 ), std::max( a, b ) );
}

TEST( SIMD, IntVMin ) {
    int_v a( 0, 1, 2, 3 );
    int_v b( 2, 2, 2, 2 );
    EXPECT_EQ( int_v( 0, 1, 2, 2 ), std::min( a, b ) );
}

#endif
