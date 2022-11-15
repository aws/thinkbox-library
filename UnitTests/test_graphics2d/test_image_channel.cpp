// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics2d/image_channel.hpp>

using namespace frantic::graphics2d;

struct test_function {
    inline float operator()( const float& left, const float& right ) { return left * 10.f + right; }
};

TEST( ImageChannel, ApplyFunction ) {
    image_channel<float> leftChannel;
    leftChannel.resize( size2( 2, 2 ) );
    leftChannel.set_pixel( 0, 0, 0.f );
    leftChannel.set_pixel( 0, 1, 1.f );
    leftChannel.set_pixel( 1, 0, 2.f );
    leftChannel.set_pixel( 1, 1, 3.f );

    image_channel<float> rightChannel;
    rightChannel.resize( size2( 2, 2 ) );
    rightChannel.set_pixel( 0, 0, 4.f );
    rightChannel.set_pixel( 0, 1, 5.f );
    rightChannel.set_pixel( 1, 0, 6.f );
    rightChannel.set_pixel( 1, 1, 7.f );

    leftChannel.apply_function( rightChannel, test_function() );

    EXPECT_EQ( 4.f, leftChannel.get_pixel( 0, 0 ) );
    EXPECT_EQ( 15.f, leftChannel.get_pixel( 0, 1 ) );
    EXPECT_EQ( 26.f, leftChannel.get_pixel( 1, 0 ) );
    EXPECT_EQ( 37.f, leftChannel.get_pixel( 1, 1 ) );
}
