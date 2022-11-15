// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <algorithm>
#include <cmath>
#include <stdexcept>

// The FRANTIC_DISABLE_IOSTREAM macro disables << functionality and the .str() function.
#ifndef FRANTIC_DISABLE_IOSTREAM
#include <iostream>
#endif

#include <boost/algorithm/clamp.hpp>
#include <boost/cstdint.hpp>
#include <frantic/misc/string_functions.hpp>

namespace frantic {
namespace graphics {

// Pack vector3f with an alignment of 2, so they're guaranteed to be adjacent
#pragma pack( push, 4 )

class color3f {
  public:
    float r, g, b;

    enum { Channels = 3 };

    //////////////
    // Constructors
    //////////////

    color3f( float R, float G, float B ) {
        r = R;
        g = G;
        b = B;
    }

    explicit color3f( float I ) {
        r = I;
        g = I;
        b = I;
    }

    explicit color3f( int I ) {
        r = (float)I;
        g = (float)I;
        b = (float)I;
    }

    explicit color3f( const float* col ) {
        r = col[0];
        g = col[1];
        b = col[2];
    }

    color3f() {
        r = 0;
        g = 0;
        b = 0;
    }

#ifdef _COLOR_H
    color3f( const Color& c ) {
        r = c.r;
        g = c.g;
        b = c.b;
    }
#endif
    static color3f white() { return color3f( 1, 1, 1 ); }

    static color3f red() { return color3f( 1, 0, 0 ); }

    static color3f green() { return color3f( 0, 1, 0 ); }

    static color3f blue() { return color3f( 0, 0, 1 ); }

    static color3f black() { return color3f( 0 ); }

    static color3f from_RGBA( boost::uint32_t c ) {
        return color3f( ( c & 0xff ) / 255.f, ( ( c & 0xff00 ) >> 8 ) / 255.f, ( ( c & 0xff0000 ) >> 16 ) / 255.f );
    }

    boost::uint32_t to_RGBA() const {
        int cr = ( std::max )( 0, ( std::min )( 255, static_cast<int>( floorf( r * 255.f + 0.5f ) ) ) ),
            cg = ( std::max )( 0, ( std::min )( 255, static_cast<int>( floorf( g * 255.f + 0.5f ) ) ) ),
            cb = ( std::max )( 0, ( std::min )( 255, static_cast<int>( floorf( b * 255.f + 0.5f ) ) ) );
        return cr | ( cg << 8 ) | ( cb << 16 );
    }

    float& operator[]( int i ) { return ( &r )[i]; }

    float operator[]( int i ) const { return ( &r )[i]; }

    //////////////
    // Queries
    //////////////

    float max_abs_component() const { return ( std::max )( fabsf( r ), ( std::max )( fabsf( g ), fabsf( b ) ) ); }

    float component_sum() const { return r + g + b; }

    static bool equals_relative_error( const color3f& a, const color3f& b, float maxError ) {
        float absValue = ( std::max )( a.max_abs_component(), b.max_abs_component() );
        maxError *= absValue;
        if( maxError < 0.0001f ) // Limit colors max abs error to 10^-4
            maxError = 0.0001f;

        float difference =
            sqrtf( ( a.r - b.r ) * ( a.r - b.r ) + ( a.g - b.g ) * ( a.g - b.g ) + ( a.b - b.b ) * ( a.b - b.b ) );
        return difference < maxError;
    }

    static color3f abs( const color3f& c ) { return color3f( fabsf( c.r ), fabsf( c.g ), fabsf( c.b ) ); }

    //////////////
    // Operators
    //////////////

    void set( float r_, float g_, float b_ ) {
        r = r_;
        g = g_;
        b = b_;
    }

    // There's no alpha channel, so just do an additive blend
    void blend_over( const color3f& color ) { operator+=( color ); }

    // There's no alpha channel, so just do an additive blend
    void blend_under( const color3f& color ) { operator+=( color ); }

    void apply_gamma( float gamma, float one = 1 ) {
        float exponent = 1 / gamma;
        r = one * pow( r / one, exponent );
        g = one * pow( g / one, exponent );
        b = one * pow( b / one, exponent );
    }

    void clamp( float minValue = 0.0f, float maxValue = 1.0f ) {
        r = ( std::max )( minValue, ( std::min )( r, maxValue ) );
        g = ( std::max )( minValue, ( std::min )( g, maxValue ) );
        b = ( std::max )( minValue, ( std::min )( b, maxValue ) );
    }

    color3f to_clamped( float minValue = 0.0f, float maxValue = 1.0f ) const {
        return color3f( ( std::max )( minValue, ( std::min )( r, maxValue ) ),
                        ( std::max )( minValue, ( std::min )( g, maxValue ) ),
                        ( std::max )( minValue, ( std::min )( b, maxValue ) ) );
    }

    float hue() const {
        float minValue = ( std::min )( r, ( std::min )( g, b ) );
        float maxValue = ( std::max )( r, ( std::max )( g, b ) );
        float delta = maxValue - minValue;
        if( delta == 0.0f || maxValue == 0.0f )
            return 0.0f;

        float h = 0.0f;

        if( r == maxValue )
            h = ( g - b ) / delta;
        else if( g == maxValue )
            h = 2.0f + ( b - r ) / delta;
        else if( b == maxValue )
            h = 4.0f + ( r - g ) / delta;

        h *= 60.0f;
        if( h < 0.0f )
            h += 360.0f;

        return h;
    }

    // HSV saturation
    float saturation() const {
        float minValue = ( std::min )( r, ( std::min )( g, b ) );
        float maxValue = ( std::max )( r, ( std::max )( g, b ) );
        float delta = maxValue - minValue;
        return ( maxValue == 0.0f ) ? 0.0f : delta / maxValue;
    }

    // some people call this "value"
    float max_component() const { return ( std::max )( r, ( std::max )( g, b ) ); }

    inline float min_component() const { return ( std::min )( r, ( std::min )( g, b ) ); }

    inline float chroma() const {
        const float max = max_component();
        const float min = min_component();
        return max - min;
    }

    inline float value() const { return max_component(); }

    inline float lightness() const { return 0.5f * ( min_component() + max_component() ); }

    inline float hsl_saturation() const {
        const float c = chroma();
        if( c == 0.0f ) {
            return 0.0f;
        }

        const float denominator = 1 - std::abs( 2 * lightness() - 1 );
        return c / denominator;
    }

    enum hsv_value_clamping { unclamped, clamped };

    static color3f from_hsv( float hue, float saturation, float value, hsv_value_clamping clampValue = clamped ) {
        const float h = hue < 0.0f || hue >= 360.0f ? std::fmod( hue, 360.0f ) : hue;
        const float s = boost::algorithm::clamp( saturation, 0.0f, 1.0f );
        const float v = clampValue == clamped ? boost::algorithm::clamp( value, 0.0f, 1.0f ) : value;

        const float chroma = v * s;
        const float huePrime = h / 60.0f;
        const float x = chroma * ( 1.0f - std::abs( std::fmod( huePrime, 2.0f ) - 1.0f ) );

        color3f result;

        if( huePrime < 1.0f ) {
            result = color3f( chroma, x, 0.0f );
        } else if( huePrime < 2.0f ) {
            result = color3f( x, chroma, 0.0f );
        } else if( huePrime < 3.0f ) {
            result = color3f( 0.0f, chroma, x );
        } else if( huePrime < 4.0f ) {
            result = color3f( 0.0f, x, chroma );
        } else if( huePrime < 5.0f ) {
            result = color3f( x, 0.0f, chroma );
        } else if( huePrime < 6.0f ) {
            result = color3f( chroma, 0.0f, x );
        } else {
            result = color3f( 0.0f, 0.0f, 0.0f );
        }

        const float min = v - chroma;

        result.r += min;
        result.g += min;
        result.b += min;

        return result;
    }

    static color3f from_hsl( float hue, float saturation, float lightness ) {
        const float h = hue < 0.0f || hue >= 360.0f ? std::fmod( hue, 360.0f ) : hue;
        const float s = boost::algorithm::clamp( saturation, 0.0f, 1.0f );
        const float l = boost::algorithm::clamp( lightness, 0.0f, 1.0f );

        const float chroma = ( 1.0f - std::abs( 2.0f * l - 1.0f ) ) * s;
        const float huePrime = h / 60.0f;
        const float x = chroma * ( 1.0f - std::abs( std::fmod( huePrime, 2.0f ) - 1.0f ) );

        color3f result;

        if( huePrime < 1.0f ) {
            result = color3f( chroma, x, 0.0f );
        } else if( huePrime < 2.0f ) {
            result = color3f( x, chroma, 0.0f );
        } else if( huePrime < 3.0f ) {
            result = color3f( 0.0f, chroma, x );
        } else if( huePrime < 4.0f ) {
            result = color3f( 0.0f, x, chroma );
        } else if( huePrime < 5.0f ) {
            result = color3f( x, 0.0f, chroma );
        } else if( huePrime < 6.0f ) {
            result = color3f( chroma, 0.0f, x );
        } else {
            result = color3f( 0.0f, 0.0f, 0.0f );
        }

        const float min = l - 0.5f * chroma;

        result.r += min;
        result.g += min;
        result.b += min;

        return result;
    }

    color3f fabs() const { return color3f( ::fabs( r ), ::fabs( g ), ::fabs( b ) ); }

    color3f operator-() const { return color3f( -r, -g, -b ); }

    color3f& operator+=( const color3f& a ) {
        r += a.r;
        g += a.g;
        b += a.b;
        return *this;
    }

    color3f& operator-=( const color3f& a ) {
        r -= a.r;
        g -= a.g;
        b -= a.b;
        return *this;
    }

    color3f& operator*=( const color3f& a ) {
        r *= a.r;
        g *= a.g;
        b *= a.b;
        return *this;
    }

    color3f& operator*=( float x ) {
        r *= x;
        g *= x;
        b *= x;
        return *this;
    }

    color3f& operator/=( float x ) {
        r /= x;
        g /= x;
        b /= x;
        return *this;
    }

    bool operator==( const color3f& c ) const { return r == c.r && g == c.g && b == c.b; }

    bool operator!=( const color3f& c ) const { return r != c.r || g != c.g || b != c.b; }

#ifdef _COLOR_H
    operator Color() const { return Color( r, g, b ); }
#endif

    std::string str() const;

    static color3f parse( const frantic::tstring& input ) {
        std::vector<frantic::tstring> values;
        frantic::strings::split( frantic::strings::string_replace( input, _T("color "), _T("") ), values,
                                 _T("()[], ") );
        if( values.size() == 1 ) {
            float X;
            try {
                X = boost::lexical_cast<float>( values[0] );
            } catch( boost::bad_lexical_cast& ) {
                throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                          frantic::strings::to_string( input ) + "\"" );
            }
            return color3f( X );
        } else if( values.size() == 3 ) {
            float X, Y, Z;
            try {
                X = boost::lexical_cast<float>( values[0] );
                Y = boost::lexical_cast<float>( values[1] );
                Z = boost::lexical_cast<float>( values[2] );
            } catch( boost::bad_lexical_cast& ) {
                throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                          frantic::strings::to_string( input ) + "\"" );
            }

            return color3f( X, Y, Z );
        } else {
            throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                      frantic::strings::to_string( input ) + "\"" );
        }
    }

    color3f to_inverse() const { return color3f( 1.f - r, 1.f - g, 1.f - b ); }
};

#pragma pack( pop )

inline color3f operator*( const color3f& a, float k ) { return color3f( a.r * k, a.g * k, a.b * k ); }

inline color3f operator*( const color3f& a, const color3f& b ) { return color3f( a.r * b.r, a.g * b.g, a.b * b.b ); }

inline color3f operator*( float k, const color3f& a ) { return color3f( k * a.r, k * a.g, k * a.b ); }

inline color3f operator/( const color3f& a, float k ) { return color3f( a.r / k, a.g / k, a.b / k ); }

inline color3f operator/( const color3f& a, const color3f& b ) { return color3f( a.r / b.r, a.g / b.g, a.b / b.b ); }

inline color3f operator+( const color3f& a, const color3f& b ) { return color3f( a.r + b.r, a.g + b.g, a.b + b.b ); }

inline color3f operator-( const color3f& a, const color3f& b ) { return color3f( a.r - b.r, a.g - b.g, a.b - b.b ); }

#ifndef FRANTIC_DISABLE_IOSTREAM
template <class ElemType>
inline std::basic_ostream<ElemType>& operator<<( std::basic_ostream<ElemType>& out, const color3f& a ) {
    out << "(color " << a.r << ", " << a.g << ", " << a.b << " )";
    return out;
}

inline std::string color3f::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

#endif

} // namespace graphics
} // namespace frantic
