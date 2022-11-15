// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics2d/framebuffer.hpp>

#include <vector>

TEST( Framebuffer, Channels ) {
    using namespace std;
    using namespace frantic::graphics2d;
    using namespace frantic::graphics;

    vector<string> names1;
    vector<string> names2;
    vector<string> names3;

    names1.push_back( "aaa" );
    names2.push_back( "bbb" );
    names1.push_back( "ccc" );
    names2.push_back( "ddd" );
    names1.push_back( "eee" );
    names2.push_back( "eee" );
    names2.push_back( "fff" );
    names2.push_back( "ggg" );
    names1.push_back( "hhh" );

    names3.push_back( "aaa" );
    names3.push_back( "bbb" );
    names3.push_back( "ccc" );
    names3.push_back( "ddd" );
    names3.push_back( "eee" );
    names3.push_back( "fff" );
    names3.push_back( "ggg" );
    names3.push_back( "hhh" );

    framebuffer<color3f> fb1;

    fb1.set_channels( names1 );

    {
        vector<string> channels( fb1.channel_names() );
        sort( channels.begin(), channels.end() );
        for( int i = 0; i < fb1.channel_count(); i++ )
            EXPECT_TRUE( channels[i] == names1[i] );
    }

    fb1.ensure_channels( names2 );
    {
        vector<string> channels( fb1.channel_names() );
        sort( channels.begin(), channels.end() );
        for( int i = 0; i < fb1.channel_count(); i++ )
            EXPECT_TRUE( channels[i] == names3[i] );
    }
}

TEST( Framebuffer, FloodFill ) {
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;

    color3f clearColor( 0.0f, 0.0f, 0.0f );

    framebuffer<color3f> baseBuffer( size2( 640, 480 ) );
    baseBuffer.fill( clearColor );

    std::vector<vector2> polyline;

    polyline.push_back( vector2( 10, 10 ) );
    polyline.push_back( vector2( 30, 10 ) );
    polyline.push_back( vector2( 30, 40 ) );
    polyline.push_back( vector2( 80, 40 ) );
    polyline.push_back( vector2( 80, 20 ) );
    polyline.push_back( vector2( 150, 20 ) );
    polyline.push_back( vector2( 230, 100 ) );
    polyline.push_back( vector2( 180, 200 ) );
    polyline.push_back( vector2( 380, 250 ) );
    polyline.push_back( vector2( 5, 260 ) );

    color3f boundaryColor( 1.0f, 0.0f, 0.0f );

    for( size_t i = 0; i < polyline.size(); ++i ) {
        baseBuffer.draw_line( polyline[i], polyline[( i + 1 ) % polyline.size()], boundaryColor );
    }

    std::vector<vector2> pointsInside;

    pointsInside.push_back( vector2( 20, 20 ) );
    pointsInside.push_back( vector2( 90, 30 ) );
    pointsInside.push_back( vector2( 220, 100 ) );
    pointsInside.push_back( vector2( 170, 200 ) );
    pointsInside.push_back( vector2( 370, 249 ) );
    pointsInside.push_back( vector2( 9, 251 ) );

    color3f fillColor( 0.0f, 1.0f, 0.0 );

    for( size_t i = 0; i < pointsInside.size(); ++i ) {
        framebuffer<color3f> copyBuffer( baseBuffer );

        copyBuffer.flood_fill( pointsInside[i], fillColor, boundaryColor );

        for( size_t j = 0; j < pointsInside.size(); ++j ) {
            color3f found = copyBuffer.get_pixel( pointsInside[j] );

            EXPECT_EQ( fillColor.r, found.r );
            EXPECT_EQ( fillColor.g, found.g );
            EXPECT_EQ( fillColor.b, found.b );
        }

        for( size_t y = 0; y < copyBuffer.height(); ++y ) {

            color3f outsideLeft = copyBuffer.get_pixel( vector2( 0, vector2::value_type( y ) ) );

            EXPECT_EQ( clearColor.r, outsideLeft.r );
            EXPECT_EQ( clearColor.g, outsideLeft.g );
            EXPECT_EQ( clearColor.b, outsideLeft.b );

            color3f outsideRight = copyBuffer.get_pixel(
                vector2( vector2::value_type( copyBuffer.width() - 1 ), vector2::value_type( y ) ) );

            EXPECT_EQ( clearColor.r, outsideRight.r );
            EXPECT_EQ( clearColor.g, outsideRight.g );
            EXPECT_EQ( clearColor.b, outsideRight.b );
        }

        for( size_t x = 0; x < copyBuffer.width(); ++x ) {
            color3f outsideTop = copyBuffer.get_pixel( vector2( vector2::value_type( x ), 0 ) );

            EXPECT_EQ( clearColor.r, outsideTop.r );
            EXPECT_EQ( clearColor.g, outsideTop.g );
            EXPECT_EQ( clearColor.b, outsideTop.b );

            color3f outsideBottom = copyBuffer.get_pixel(
                vector2( vector2::value_type( x ), vector2::value_type( copyBuffer.height() - 1 ) ) );

            EXPECT_EQ( clearColor.r, outsideBottom.r );
            EXPECT_EQ( clearColor.g, outsideBottom.g );
            EXPECT_EQ( clearColor.b, outsideBottom.b );
        }
    }
}
