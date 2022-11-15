// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>

#include <frantic/graphics/color3f.hpp>

namespace frantic {
namespace graphics {

// This represents and rgb opacity.  Adding alpha3f's together is the same as occluding the alpha
// values in front of each other.  Multiplying alpha3f's by floating point values scales the transparency
class alpha3f {
  public:
    float ar, ag, ab;

    enum { Channels = 3 };

    //////////////
    // Constructors
    //////////////

    alpha3f( float aR, float aG, float aB ) {
        ar = aR;
        ag = aG;
        ab = aB;
    }

    explicit alpha3f( float I ) {
        ar = I;
        ag = I;
        ab = I;
    }

    explicit alpha3f( int I ) {
        ar = (float)I;
        ag = (float)I;
        ab = (float)I;
    }

    explicit alpha3f( const float* col ) {
        ar = col[0];
        ag = col[1];
        ab = col[2];
    }

    alpha3f() {
        ar = 0;
        ag = 0;
        ab = 0;
    }

    explicit alpha3f( const color3f& c ) {
        ar = c.r;
        ag = c.g;
        ab = c.b;
    }

#ifdef _COLOR_H
    explicit alpha3f( const Color& c ) {
        ar = c.r;
        ag = c.g;
        ab = c.b;
    }
#endif

    //////////////
    // Queries
    //////////////

    float max_abs_component() const { return ( std::max )( fabs( ar ), ( std::max )( fabs( ag ), fabs( ab ) ) ); }

    static bool equals_relative_error( const alpha3f& a, const alpha3f& b, float maxError ) {
        float absValue = ( std::max )( a.max_abs_component(), b.max_abs_component() );
        maxError *= absValue;
        if( maxError < 0.0001f ) // Limit colors max abs error to 10^-4
            maxError = 0.0001f;

        float difference = sqrtf( ( a.ar - b.ar ) * ( a.ar - b.ar ) + ( a.ag - b.ag ) * ( a.ag - b.ag ) +
                                  ( a.ab - b.ab ) * ( a.ab - b.ab ) );
        return difference < maxError;
    }

    float component_sum() const { return ar + ag + ab; }

    static alpha3f abs( const alpha3f& c ) { return alpha3f( fabsf( c.ar ), fabsf( c.ag ), fabsf( c.ab ) ); }

    //////////////
    // Operators
    //////////////

    float to_float() const {
        // Not sure what the best formula to use is
        return 0.333333f * ( ar + ag + ab );
    }

    color3f occlude( const color3f& c ) const {
        return color3f( ( 1 - ar ) * c.r, ( 1 - ag ) * c.g, ( 1 - ab ) * c.b );
    }

    alpha3f occlude( const alpha3f& a ) const {
        return alpha3f( ( 1 - ar ) * a.ar, ( 1 - ag ) * a.ag, ( 1 - ab ) * a.ab );
    }

    color3f premultiply( const color3f& c ) const { return color3f( ar * c.r, ag * c.g, ab * c.b ); }

    void blend_over( const alpha3f& a ) {
        ar = 1 - ( 1 - ar ) * ( 1 - a.ar );
        ag = 1 - ( 1 - ag ) * ( 1 - a.ag );
        ab = 1 - ( 1 - ab ) * ( 1 - a.ab );
    }

    void blend_under( const alpha3f& a ) {
        ar = 1 - ( 1 - ar ) * ( 1 - a.ar );
        ag = 1 - ( 1 - ag ) * ( 1 - a.ag );
        ab = 1 - ( 1 - ab ) * ( 1 - a.ab );
    }

    alpha3f& operator+=( const alpha3f& a ) {
        ar += a.ar;
        ag += a.ag;
        ab += a.ab;
        return *this;
    }

    alpha3f& operator*=( float x ) {
        ar *= x;
        ag *= x;
        ab *= x;
        return *this;
    }

    alpha3f& operator/=( float x ) {
        ar /= x;
        ag /= x;
        ab /= x;
        return *this;
    }

    alpha3f& operator*=( const alpha3f& rhs ) {
        ar *= rhs.ar;
        ag *= rhs.ag;
        ab *= rhs.ab;
        return *this;
    }

    std::string str() const;

    bool operator!=( const alpha3f& other ) const { return ar != other.ar || ag != other.ag || ab != other.ab; }

    bool operator==( const alpha3f& other ) const { return ar == other.ar && ag == other.ag && ab == other.ab; }

    alpha3f to_inverse() const { return alpha3f( 1.f - ar, 1.f - ag, 1.f - ab ); }
};

inline alpha3f operator+( const alpha3f& a, const alpha3f& b ) {
    return alpha3f( a.ar + b.ar, a.ag + b.ag, a.ab + b.ab );
}

inline alpha3f operator-( const alpha3f& a, const alpha3f& b ) {
    return alpha3f( a.ar - b.ar, a.ag - b.ag, a.ab - b.ab );
}

inline alpha3f operator*( const alpha3f& a, const alpha3f& b ) {
    return alpha3f( a.ar * b.ar, a.ag * b.ag, a.ab * b.ab );
}

inline alpha3f operator*( const alpha3f& a, float k ) { return alpha3f( a.ar * k, a.ag * k, a.ab * k ); }

inline alpha3f operator*( float k, const alpha3f& a ) { return alpha3f( k * a.ar, k * a.ag, k * a.ab ); }

inline alpha3f operator/( const alpha3f& a, float k ) { return alpha3f( a.ar / k, a.ag / k, a.ab / k ); }

inline std::ostream& operator<<( std::ostream& out, const alpha3f& a ) {
    out << "(alpha " << a.ar << ", " << a.ag << ", " << a.ab << " )";
    return out;
}

inline std::string alpha3f::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

} // namespace graphics
} // namespace frantic
