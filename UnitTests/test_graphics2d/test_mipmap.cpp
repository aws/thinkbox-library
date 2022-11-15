// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics2d/mipmap.hpp>

frantic::graphics2d::image_channel<frantic::graphics::color3f> generate_base() {
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;

    image_channel<color3f> base;
    base.set_size( size2( 4, 4 ) );

    base.set_pixel( 0, 0, color3f::blue() );
    base.set_pixel( 0, 1, color3f::blue() );
    base.set_pixel( 1, 0, color3f::blue() );
    base.set_pixel( 1, 1, color3f::blue() );

    base.set_pixel( 2, 0, color3f::red() );
    base.set_pixel( 3, 0, color3f::green() );
    base.set_pixel( 2, 1, color3f::blue() );
    base.set_pixel( 3, 1, color3f::white() );

    base.set_pixel( 0, 2, color3f::black() );
    base.set_pixel( 0, 3, color3f::black() );
    base.set_pixel( 1, 2, color3f::white() );
    base.set_pixel( 1, 3, color3f::white() );

    base.set_pixel( 2, 2, color3f::red() );
    base.set_pixel( 2, 3, color3f::red() );
    base.set_pixel( 3, 2, color3f::red() );
    base.set_pixel( 3, 3, color3f::blue() );

    return base;
}

TEST( Mipmap, generate_mipmap ) {
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;

    const image_channel<color3f> base = generate_base();
    const std::vector<image_channel<color3f>> mipmap = generate_mipmap<color3f, four_pixel_mean<color3f>>( base );

    ASSERT_EQ( 3, mipmap.size() );

    // Test 4x4
    for( int x = 0; x < mipmap.back().width(); ++x ) {
        for( int y = 0; y < mipmap.back().height(); ++y ) {
            ASSERT_EQ( mipmap.back().at( x, y ), base.at( x, y ) );
        }
    }

    // Test 2x2
    const static color3f upperLeft = color3f::blue();
    const static color3f upperRight( 0.5f, 0.5f, 0.5f );
    const static color3f& lowerLeft = upperRight;
    const static color3f lowerRight( 0.75f, 0.0f, 0.25f );

    const image_channel<color3f>& twoByTwo = mipmap[1];

    EXPECT_EQ( upperLeft, twoByTwo.at( 0, 0 ) );
    EXPECT_EQ( upperRight, twoByTwo.at( 1, 0 ) );
    EXPECT_EQ( lowerLeft, twoByTwo.at( 0, 1 ) );
    EXPECT_EQ( lowerRight, twoByTwo.at( 1, 1 ) );

    // Test 1x1
    const static color3f finalExpected( 0.4375f, 0.25f, 0.5625f );
    EXPECT_EQ( finalExpected, mipmap.front().at( 0, 0 ) );

    image_channel<int> intImage;
    intImage.set_size( size2( 1, 1 ) );
    const std::vector<image_channel<int>> intMipmap = generate_mipmap<int, four_pixel_mean<int>>( intImage );

    EXPECT_EQ( intMipmap.size(), 1 );
}

TEST( Mipmap, nearest_mipmap_resolution ) {
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;

    typedef image_channel<color3f> image_t;

    const image_channel<color3f> base = generate_base();
    const std::vector<image_channel<color3f>> mipmap = generate_mipmap<color3f, four_pixel_mean<color3f>>( base );

    const static color3f zeroUpperLeft( 0.4375f, 0.25f, 0.5625f );
    const static color3f oneUpperLeft = color3f::blue();
    const static color3f& twoUpperLeft = oneUpperLeft;

    const color3f& zero = nearest_mipmap_resolution<color3f>( mipmap, 0 ).at( 0, 0 );
    const color3f& one = nearest_mipmap_resolution<color3f>( mipmap, 1 ).at( 0, 0 );
    const color3f& two = nearest_mipmap_resolution<color3f>( mipmap, 2 ).at( 0, 0 );
    const color3f& three = nearest_mipmap_resolution<color3f>( mipmap, 3 ).at( 0, 0 );
    const color3f& four = nearest_mipmap_resolution<color3f>( mipmap, 4 ).at( 0, 0 );
    const color3f& oneTwentyEight = nearest_mipmap_resolution<color3f>( mipmap, 128 ).at( 0, 0 );

    EXPECT_EQ( zero, zeroUpperLeft );
    EXPECT_EQ( one, zeroUpperLeft );
    EXPECT_EQ( two, oneUpperLeft );
    EXPECT_EQ( three, twoUpperLeft );
    EXPECT_EQ( four, twoUpperLeft );
    EXPECT_EQ( oneTwentyEight, twoUpperLeft );

    EXPECT_THROW( nearest_mipmap_resolution<color3f>( std::vector<image_t>(), 128 ), std::runtime_error );
    EXPECT_THROW( nearest_mipmap_resolution<color3f>( std::vector<image_t>(), 0 ), std::runtime_error );

    image_t bigImage;
    bigImage.set_size( size2( 256, 256 ) );
    std::vector<image_t> malformed;
    malformed.push_back( base );
    malformed.push_back( base );
    malformed.push_back( bigImage );

    EXPECT_THROW( nearest_mipmap_resolution<color3f>( malformed, 128 ), std::runtime_error );
}
