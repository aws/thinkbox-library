// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>

#include <frantic/graphics/color3f.hpp>

namespace frantic {
namespace graphics {

// This represents and rgb opacity.  Adding alpha1f's together is the same as occluding the alpha
// values in front of each other.  Multiplying alpha1f's by floating point values scales the transparency
class alpha1f {
  public:
    float a;

    enum { Channels = 1 };

    //////////////
    // Constructors
    //////////////

    alpha1f( float alpha ) { a = alpha; }

    explicit alpha1f( int I ) { a = (float)I; }

    alpha1f() { a = 0; }

    //////////////
    // Queries
    //////////////

    float component_sum() const { return a; }

    static alpha1f abs( const alpha1f& c ) { return alpha1f( fabsf( c.to_float() ) ); }

    //////////////
    // Operators
    //////////////

    float to_float() const { return a; }

    color3f occlude( const color3f& c ) const { return color3f( ( 1 - a ) * c.r, ( 1 - a ) * c.g, ( 1 - a ) * c.b ); }

    color3f premultiply( const color3f& c ) const { return color3f( a * c.r, a * c.g, a * c.b ); }

    void blend_over( const alpha1f& alpha ) { a = 1 - ( 1 - a ) * ( 1 - alpha.a ); }

    void blend_under( const alpha1f& alpha ) { a = 1 - ( 1 - a ) * ( 1 - alpha.a ); }

    alpha1f& operator+=( const alpha1f& alpha ) {
        a += alpha.a;
        return *this;
    }

    alpha1f& operator*=( float x ) {
        a *= x;
        return *this;
    }

    alpha1f& operator/=( float x ) {
        a /= x;
        return *this;
    }

    bool operator!=( const alpha1f& other ) { return a != other.a; }

    bool operator==( const alpha1f& other ) { return a == other.a; }
};

inline alpha1f operator+( const alpha1f& alpha, const alpha1f& beta ) { return alpha1f( alpha.a + beta.a ); }

inline alpha1f operator*( const alpha1f& alpha, float k ) { return alpha1f( alpha.a * k ); }

inline alpha1f operator*( float k, const alpha1f& alpha ) { return alpha1f( k * alpha.a ); }

inline alpha1f operator/( const alpha1f& alpha, float k ) { return alpha1f( alpha.a / k ); }

inline std::ostream& operator<<( std::ostream& out, const alpha1f& alpha ) {
    out << "(alpha " << alpha.a << " )";
    return out;
}

} // namespace graphics
} // namespace frantic
