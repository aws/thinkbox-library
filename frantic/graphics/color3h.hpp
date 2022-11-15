// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/color3f.hpp>

#if( defined( _MSC_VER ) && _MSC_VER > 1200 ) || __GNUC__

#include <half.h>

namespace frantic {
namespace graphics {

// Pack vector3h with an alignment of 2, so they're guaranteed to be adjacent
#pragma pack( push, 2 )

class color3h {
  public:
    half r, g, b;

    //////////////
    // Constructors
    //////////////

    color3h( half R, half G, half B ) {
        r = R;
        g = G;
        b = B;
    }

    color3h( float R, float G, float B ) {
        r = R;
        g = G;
        b = B;
    }

    explicit color3h( half I ) {
        r = I;
        g = I;
        b = I;
    }

    explicit color3h( float I ) {
        r = I;
        g = I;
        b = I;
    }

    explicit color3h( int I ) {
        r = (float)I;
        g = (float)I;
        b = (float)I;
    }

    explicit color3h( const half* col ) {
        r = col[0];
        g = col[1];
        b = col[2];
    }

    explicit color3h( const float* col ) {
        r = col[0];
        g = col[1];
        b = col[2];
    }

    color3h() {
        r = 0;
        g = 0;
        b = 0;
    }

    color3h( const frantic::graphics::color3f& c ) {
        r = c.r;
        g = c.g;
        b = c.b;
    }

#ifdef _COLOR_H
    color3h( const Color& c ) {
        r = c.r;
        g = c.g;
        b = c.b;
    }
#endif

    static color3h white() { return color3h( 1, 1, 1 ); }

    static color3h red() { return color3h( 1, 0, 0 ); }

    static color3h green() { return color3h( 0, 1, 0 ); }

    static color3h blue() { return color3h( 0, 0, 1 ); }

    static color3h black() { return color3h( 0 ); }

    half& operator[]( int i ) { return ( &r )[i]; }

    const half operator[]( int i ) const { return ( &r )[i]; }

    //////////////
    // Queries
    //////////////

    float max_abs_component() const { return ( std::max )( fabsf( r ), ( std::max )( fabsf( g ), fabsf( b ) ) ); }

    float component_sum() const { return (float)r + (float)g + (float)b; }

    //////////////
    // Operators
    //////////////

    // There's no alpha channel, so just do an additive blend
    void blend_over( const color3f& color ) { operator+=( color ); }

    // There's no alpha channel, so just do an additive blend
    void blend_under( const color3f& color ) { operator+=( color ); }

    void clamp( float minValue = 0.0f, float maxValue = 1.0f ) {
        r = ( std::max )( minValue, ( std::min )( (float)r, maxValue ) );
        g = ( std::max )( minValue, ( std::min )( (float)g, maxValue ) );
        b = ( std::max )( minValue, ( std::min )( (float)b, maxValue ) );
    }

    color3f to_clamped( float minValue = 0.0f, float maxValue = 1.0f ) const {
        return color3f( ( std::max )( minValue, ( std::min )( (float)r, maxValue ) ),
                        ( std::max )( minValue, ( std::min )( (float)g, maxValue ) ),
                        ( std::max )( minValue, ( std::min )( (float)b, maxValue ) ) );
    }

    // some people call this "value"
    float max_component() const { return ( std::max )( r, ( std::max )( g, b ) ); }

    color3h operator-() const { return color3h( -r, -g, -b ); }

    color3h& operator+=( const color3f& a ) {
        r += a.r;
        g += a.g;
        b += a.b;
        return *this;
    }

    color3h& operator-=( const color3f& a ) {
        r -= a.r;
        g -= a.g;
        b -= a.b;
        return *this;
    }

    color3h& operator*=( const color3f& a ) {
        r *= a.r;
        g *= a.g;
        b *= a.b;
        return *this;
    }

    color3h& operator*=( float x ) {
        r *= x;
        g *= x;
        b *= x;
        return *this;
    }

    color3h& operator/=( float x ) {
        r /= x;
        g /= x;
        b /= x;
        return *this;
    }

    bool operator==( const color3h& c ) { return r == c.r && g == c.g && b == c.b; }

    bool operator!=( const color3h& c ) { return r != c.r || g != c.g || b != c.b; }

    operator color3f() const { return color3f( r, g, b ); }

#ifdef _COLOR_H
    operator Color() const { return Color( r, g, b ); }
#endif
};

#pragma pack( pop )

// Convert results to color3f so that intermediate calculations are done with floats as much as possible
inline color3f operator*( const color3h& a, float k ) {
    return color3f( (float)a.r * k, (float)a.g * k, (float)a.b * k );
}

inline color3f operator*( float k, const color3h& a ) {
    return color3f( k * (float)a.r, k * (float)a.g, k * (float)a.b );
}

inline color3f operator*( const color3h& a, const color3h& b ) {
    return color3f( (float)a.r * (float)b.r, (float)a.g * (float)b.g, (float)a.b * (float)b.b );
}

inline color3f operator*( const color3h& a, const color3f& b ) {
    return color3f( (float)a.r * b.r, (float)a.g * b.g, (float)a.b * b.b );
}

inline color3f operator*( const color3f& a, const color3h& b ) {
    return color3f( a.r * (float)b.r, a.g * (float)b.g, a.b * (float)b.b );
}

inline color3f operator/( const color3h& a, float k ) {
    return color3h( (float)a.r / k, (float)a.g / k, (float)a.b / k );
}

inline color3f operator+( const color3h& a, const color3h& b ) {
    return color3f( (float)a.r + (float)b.r, (float)a.g + (float)b.g, (float)a.b + (float)b.b );
}

inline color3f operator+( const color3h& a, const color3f& b ) {
    return color3f( (float)a.r + b.r, (float)a.g + b.g, (float)a.b + b.b );
}

inline color3f operator+( const color3f& a, const color3h& b ) {
    return color3f( a.r + (float)b.r, a.g + (float)b.g, a.b + (float)b.b );
}

inline color3f operator-( const color3h& a, const color3h& b ) {
    return color3f( (float)a.r - (float)b.r, (float)a.g - (float)b.g, (float)a.b - (float)b.b );
}

inline color3f operator-( const color3h& a, const color3f& b ) {
    return color3f( (float)a.r - b.r, (float)a.g - b.g, (float)a.b - b.b );
}

inline color3f operator-( const color3f& a, const color3h& b ) {
    return color3f( a.r - (float)b.r, a.g - (float)b.g, a.b - (float)b.b );
}

#ifndef FRANTIC_DISABLE_IOSTREAM
inline std::ostream& operator<<( std::ostream& out, const color3h& a ) {
    out << "(color " << (float)a.r << ", " << (float)a.g << ", " << (float)a.b << " )";
    return out;
}
#endif

} // namespace graphics
} // namespace frantic

#endif // defined(_MSC_VER) && _MSC_VER > 1200
