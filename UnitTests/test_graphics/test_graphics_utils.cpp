// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/math/constants/constants.hpp>

#include <frantic/graphics/graphics_utils.hpp>

using namespace frantic::graphics;

TEST( GraphicsUtils, GetTriangleAngle ) {
    const vector3f a( 0, 1, 0 );
    const vector3f b( 0, 0, 0 );
    const vector3f c( 1, 0, 0 );

    const float angle = get_triangle_angle( a, b, c );

    EXPECT_EQ( boost::math::constants::pi<float>() / 2, angle );
}
