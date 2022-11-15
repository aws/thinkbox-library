// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics/color3f.hpp>

using namespace frantic::graphics;

// Taken from the table here: https://en.wikipedia.org/wiki/HSL_and_HSV#Examples
// clang-format off
const static float colors[] = {
	//R     G       B       H       C       V       L       SV      SL
	1.000f, 1.000f, 1.000f, 0.0f,	0.000f, 1.000f, 1.000f, 0.000f, 0.000f,
	0.500f, 0.500f, 0.500f, 0.0f,	0.000f, 0.500f, 0.500f, 0.000f, 0.000f,
	0.000f, 0.000f, 0.000f, 0.0f,	0.000f, 0.000f, 0.000f, 0.000f, 0.000f,
	1.000f, 0.000f, 0.000f, 0.0f,	1.000f, 1.000f, 0.500f, 1.000f, 1.000f,
	0.750f, 0.750f, 0.000f, 60.0f,	0.750f, 0.750f, 0.375f, 1.000f, 1.000f,
	0.000f, 0.500f, 0.000f, 120.0f, 0.500f, 0.500f, 0.250f, 1.000f, 1.000f,
	0.500f, 1.000f, 1.000f, 180.0f, 0.500f, 1.000f, 0.750f, 0.500f, 1.000f,
	0.500f, 0.500f, 1.000f, 240.0f, 0.500f, 1.000f, 0.750f, 0.500f, 1.000f,
	0.750f, 0.250f, 0.750f, 300.0f, 0.500f, 0.750f, 0.500f, 0.667f, 0.500f,
	0.628f, 0.643f, 0.142f, 61.8f,	0.501f, 0.643f, 0.393f, 0.779f, 0.638f,
	0.255f, 0.104f, 0.918f, 251.1f, 0.814f, 0.918f, 0.511f, 0.887f, 0.832f,
	0.116f, 0.675f, 0.255f, 134.9f, 0.559f, 0.675f, 0.396f, 0.828f, 0.707f,
	0.941f, 0.785f, 0.053f, 49.5f,	0.888f, 0.941f, 0.497f, 0.944f, 0.893f,
	0.704f, 0.187f, 0.897f, 283.7f, 0.710f, 0.897f, 0.542f, 0.792f, 0.775f,
	0.931f, 0.463f, 0.316f, 14.3f,	0.615f, 0.931f, 0.624f, 0.661f, 0.817f,
	0.998f, 0.974f, 0.532f, 56.9f,	0.466f, 0.998f, 0.765f, 0.467f, 0.991f,
	0.099f, 0.795f, 0.591f, 162.4f, 0.696f, 0.795f, 0.447f, 0.875f, 0.779f,
	0.211f, 0.149f, 0.597f, 248.3f, 0.448f, 0.597f, 0.373f, 0.750f, 0.601f,
	0.495f, 0.493f, 0.721f, 240.5f, 0.228f, 0.721f, 0.607f, 0.316f, 0.290f
};
// clang-format on

const static std::size_t numColors = 19;

std::size_t row_offset( std::size_t i ) { return i * 9; }

color3f get_color( std::size_t i ) {
    const std::size_t rowOffset = row_offset( i );
    return color3f( colors[rowOffset], colors[rowOffset + 1], colors[rowOffset + 2] );
}

float get_hue( std::size_t i ) { return colors[row_offset( i ) + 3]; }

float get_chroma( std::size_t i ) { return colors[row_offset( i ) + 4]; }

float get_value( std::size_t i ) { return colors[row_offset( i ) + 5]; }

float get_lightness( std::size_t i ) { return colors[row_offset( i ) + 6]; }

float get_hsv_saturation( std::size_t i ) { return colors[row_offset( i ) + 7]; }

float get_hsl_saturation( std::size_t i ) { return colors[row_offset( i ) + 8]; }

float round1( float in ) { return std::floor( in * 10.0f + 0.5f ) / 10.0f; }

float round3( float in ) { return std::floor( in * 1000.0f + 0.5f ) / 1000.0f; }

color3f round3( const color3f& in ) {
    color3f result;
    result.r = round3( in.r );
    result.g = round3( in.g );
    result.b = round3( in.b );
    return result;
}

bool within_rounding_error( const color3f& expected, const color3f& actual ) {
    return std::abs( actual.r - expected.r ) <= 0.001f && std::abs( actual.g - expected.g ) <= 0.001f &&
           std::abs( actual.b - expected.b ) <= 0.001f;
}

TEST( color3f, HSLHSV ) {
    for( std::size_t i = 0; i < 19; ++i ) {
        const color3f expected = get_color( i );

        // Test component calculators

        EXPECT_EQ( round1( expected.hue() ), get_hue( i ) );
        EXPECT_EQ( round3( expected.chroma() ), get_chroma( i ) );
        EXPECT_EQ( round3( expected.value() ), get_value( i ) );
        EXPECT_EQ( round3( expected.lightness() ), get_lightness( i ) );
        EXPECT_EQ( round3( expected.saturation() ), get_hsv_saturation( i ) );
        EXPECT_EQ( round3( expected.hsl_saturation() ), get_hsl_saturation( i ) );

        // Test generators

        const color3f fromHSV = round3( color3f::from_hsv( get_hue( i ), get_hsv_saturation( i ), get_value( i ) ) );
        const color3f fromHSL =
            round3( color3f::from_hsl( get_hue( i ), get_hsl_saturation( i ), get_lightness( i ) ) );

        EXPECT_TRUE( within_rounding_error( expected, fromHSV ) );
        EXPECT_TRUE( within_rounding_error( expected, fromHSL ) );
    }

    // Test Input Fixing

    EXPECT_EQ( color3f::from_hsv( 0.0f, 0.0f, 0.0f ), color3f::from_hsv( -1.0f, -1.0f, -1.0f ) );
    EXPECT_EQ( color3f::from_hsv( 0.0f, 1.0f, 1.0f ), color3f::from_hsv( 360.0f, 2.0f, 2.0f ) );

    EXPECT_EQ( color3f::from_hsl( 0.0f, 0.0f, 0.0f ), color3f::from_hsl( -1.0f, -1.0f, -1.0f ) );
    EXPECT_EQ( color3f::from_hsl( 0.0f, 1.0f, 1.0f ), color3f::from_hsl( 360.0f, 2.0f, 2.0f ) );
}

TEST( color3f, HSVClamping ) {
    EXPECT_EQ( color3f( 1.0f, 0.0f, 0.0f ), color3f::from_hsv( 0.0f, 1.0f, 2.0f, color3f::clamped ) );
    EXPECT_EQ( color3f( 2.0f, 0.0f, 0.0f ), color3f::from_hsv( 0.0f, 1.0f, 2.0f, color3f::unclamped ) );
}

TEST( color3f, to_inverse ) {
    const color3f original( 0.25f, 0.5f, 0.75f );
    EXPECT_EQ( color3f( 0.75f, 0.5f, 0.25f ), original.to_inverse() );
}
