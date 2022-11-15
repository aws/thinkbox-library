// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics/ray3f.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::math;

TEST( Ray3f, Basics ) {
    vector3f v0( 0, 0, 0 ), v1( 0, 1, 0 ), v2( 1, 0, 0 );
    ray3f ray( vector3f( 0.45f, 0.45f, 1.f ), vector3f( 0, 0, -1.f ) );
    double t;
    vector3f baryCentric;
    EXPECT_TRUE( ray.intersect_with_triangle( v0, v1, v2, t, baryCentric ) );
    EXPECT_NEAR( t, 1, 0.0001 );
    EXPECT_NEAR( baryCentric.x, 0.45f, 0.001f );
    EXPECT_NEAR( baryCentric.y, 0.45f, 0.001f );
    ray.set( vector3f( 0.55f, 0.55f, 1.f ), vector3f( 0, 0, -1.f ) );
    EXPECT_FALSE( ray.intersect_with_triangle( v0, v1, v2, t, baryCentric ) );
}

TEST( Ray3f, NearestPointParametersParallel ) {
    vector3f p0( 10.f, 11.f, 10.f ), v0( 3.f, 1.1f, 7.f );
    vector3f p1( 2.2f, 1.f, 5.f ), v1( 3.f, 1.1f, 7.f );
    ray3f r0( p0, v0 ), r1( p1, v1 );
    std::pair<float, float> soln = r0.nearest_point_parameters( r1 );
    EXPECT_TRUE( soln.first != soln.first );   // NaN
    EXPECT_TRUE( soln.second != soln.second ); // NaN
}

TEST( Ray3f, NearestPointParametersAntiparallel ) {
    vector3f p0( 46.f, -4.135f, 1.f ), v0( 10.7f, 3.f, -2.6f );
    vector3f p1( 2.f, -100.f, 6.6f ), v1( -10.7f, -3.f, 2.6f );
    ray3f r0( p0, v0 ), r1( p1, v1 );
    std::pair<float, float> soln = r0.nearest_point_parameters( r1 );
    EXPECT_TRUE( soln.first != soln.first );   // NaN
    EXPECT_TRUE( soln.second != soln.second ); // NaN
}

TEST( Ray3f, NearestPointParametersNegative ) {
    // Looks like this:
    //
    //        ^
    //        |
    //        |
    //        .
    //
    // <--------------.
    vector3f p0( 4.f, 0.f, -2.5f ), v0( -1.f, 0.f, 0.f );
    vector3f p1( 0.f, 1.f, 153.17f ), v1( 0.f, 1.f, 0.f );
    ray3f r0( p0, v0 ), r1( p1, v1 );
    std::pair<float, float> soln = r0.nearest_point_parameters( r1 );
    EXPECT_EQ( 4.f, soln.first );
    EXPECT_EQ( -1.f, soln.second );
}

TEST( Ray3f, NearestPointParametersBothPositiveNonintersecting ) {
    vector3f p0( -1.f, 1.f, 1.f ), v0( 1.f, 1.f, 0.f );
    vector3f p1( 2.f, 1.f, -1.f ), v1( -1.f, 1.f, 0.f );
    ray3f r0( p0, v0 ), r1( p1, v1 );
    std::pair<float, float> soln = r0.nearest_point_parameters( r1 );
    EXPECT_EQ( 1.5f, soln.first );
    EXPECT_EQ( 1.5f, soln.second );
}

TEST( Ray3f, NearestPointParametersIntersecting ) {
    vector3f p0( -1.f, 1.f, 4.3f ), v0( 1.f, 1.f, 0.f );
    vector3f p1( 3.f, 1.f, 4.3f ), v1( -1.f, 1.f, 0.f );
    ray3f r0( p0, v0 ), r1( p1, v1 );
    std::pair<float, float> soln = r0.nearest_point_parameters( r1 );
    EXPECT_EQ( 2.f, soln.first );
    EXPECT_EQ( 2.f, soln.second );
}
