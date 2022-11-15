// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics2d/vector2f.hpp>

TEST( Vector2f, TriangleMethods ) {
    using namespace std;
    using namespace frantic::graphics2d;
    vector2f a, b, c;
    a = vector2f( 0, 0 );
    b = vector2f( 1, 0 );
    c = vector2f( 0, 1 );
    EXPECT_LT( fabsf( 0.5f - vector2f::triangle_area( a, b, c ) ), 0.0001f );
    EXPECT_LT( 0, vector2f::triangle_curvature( a, b, c ) );
    EXPECT_LT( vector2f::triangle_curvature( c, b, a ), 0 );
}
