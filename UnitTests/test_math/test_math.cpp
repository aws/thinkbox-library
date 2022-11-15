// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/random.hpp>

#include <frantic/math/utils.hpp>

using namespace std;

TEST( Math, Utils ) {
    double zero = 0.0;
    double inf = numeric_limits<double>::infinity();
    double nan = numeric_limits<double>::quiet_NaN();
    // is_finite
    EXPECT_TRUE( frantic::math::is_finite( zero ) );
    EXPECT_FALSE( frantic::math::is_finite( inf ) );
    EXPECT_FALSE( frantic::math::is_finite( -inf ) );
    EXPECT_FALSE( frantic::math::is_infinite( nan ) );
    // is_infinite
    EXPECT_FALSE( frantic::math::is_infinite( zero ) );
    EXPECT_TRUE( frantic::math::is_infinite( inf ) );
    EXPECT_TRUE( frantic::math::is_infinite( -inf ) );
    EXPECT_FALSE( frantic::math::is_infinite( nan ) );
    // is_nan
    EXPECT_FALSE( frantic::math::is_nan( zero ) );
    EXPECT_FALSE( frantic::math::is_nan( inf ) );
    EXPECT_FALSE( frantic::math::is_nan( -inf ) );
    EXPECT_TRUE( frantic::math::is_nan( nan ) );
}
