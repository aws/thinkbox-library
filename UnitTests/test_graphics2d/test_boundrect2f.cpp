// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics2d/boundrect2f.hpp>

TEST( Boundrect2f, AreaIntersectingLineSegment ) {
    using namespace std;
    using namespace frantic::graphics2d;

    boundrect2f br( vector2f( 0, 0 ), vector2f( 1, 1 ) );

    EXPECT_TRUE( br.is_area_intersecting_line_segment( vector2f( -100, 1 ), vector2f( 100, 0.99999f ) ) );
    EXPECT_TRUE( br.is_area_intersecting_line_segment( vector2f( -10, 0.5f ), vector2f( 10, 0.5f ) ) );
    EXPECT_TRUE( br.is_area_intersecting_line_segment( vector2f( 10, 0.5f ), vector2f( -10, 0.5f ) ) );
    EXPECT_TRUE( br.is_area_intersecting_line_segment( vector2f( 0, -10 ), vector2f( 1, 10 ) ) );
    EXPECT_TRUE( br.is_area_intersecting_line_segment( vector2f( 0.4999f, 10 ), vector2f( 0.5f, -10 ) ) );

    EXPECT_TRUE( br.is_area_intersecting_line_segment( vector2f( 0, 1.9f ), vector2f( 2, 0 ) ) );
    EXPECT_TRUE( !br.is_area_intersecting_line_segment( vector2f( 0, 2.1f ), vector2f( 2, 0 ) ) );
    EXPECT_TRUE( !br.is_area_intersecting_line_segment( vector2f( -100, .5f ), vector2f( 0.5f, 1.01f ) ) );
    EXPECT_TRUE( br.is_area_intersecting_line_segment( vector2f( -100, .5f ), vector2f( 0.5f, 1.001f ) ) );
}
