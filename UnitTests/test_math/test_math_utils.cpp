// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/math/utils.hpp>

using namespace frantic::math;

TEST( MathUtils, Log2Uint32 ) {
    const boost::uint32_t lsb = 1;
    EXPECT_EQ( 0, log2_uint32( lsb ) );

    const boost::uint32_t msb = 0xFFFFFFFF;
    EXPECT_EQ( 31, log2_uint32( msb ) );
}

TEST( MathUtils, IsPowerOfTwo ) {
    EXPECT_FALSE( is_power_of_two( 0 ) );
    EXPECT_TRUE( is_power_of_two( 1 ) );
    EXPECT_TRUE( is_power_of_two( 2 ) );
    EXPECT_FALSE( is_power_of_two( 3 ) );
    EXPECT_TRUE( is_power_of_two( 4 ) );
}
