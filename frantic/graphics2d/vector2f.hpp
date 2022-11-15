// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

#include <frantic/graphics2d/vector2t.hpp>

namespace frantic {
namespace graphics2d {

class vector2f : public vector2t<float, vector2f> {
  public:
    typedef value_type float_type;

    vector2f()
        : vector2t<float, vector2f>( 0.0f, 0.0f ) {}
    vector2f( float X, float Y )
        : vector2t<float, vector2f>( X, Y ) {}
    vector2f( const std::pair<float, float>& p )
        : vector2t<float, vector2f>( p.first, p.second ) {}
    explicit vector2f( float X )
        : vector2t<float, vector2f>( X, X ) {}
    explicit vector2f( float* vec )
        : vector2t<float, vector2f>( vec[0], vec[1] ) {}
    vector2f( const vector2f& v )
        : vector2t<float, vector2f>( v.x, v.y ) {}

    // TODO: We should use the boost random number generator for high quality and fast random numbers
    static vector2f from_random() { return vector2f( (float)rand() / RAND_MAX, (float)rand() / RAND_MAX ); }

    template <class RandomNumberGenerator>
    static vector2f from_unit_disk_random( RandomNumberGenerator& rng ) {
        // There's probably a better approach which uses fewer rng() calls
        vector2f result;
        result.x = 2 * rng() - 1;
        result.y = 2 * rng() - 1;

        while( result.get_magnitude_squared() > 1 ) {
            result.x = 2 * rng() - 1;
            result.y = 2 * rng() - 1;
        }

        return result;
    }

    // TODO: We should use the boost random number generator for high quality and fast random numbers
    static vector2f from_random_gaussian();

    template <class RandomNumberGenerator>
    static vector2f from_random_gaussian( RandomNumberGenerator& rng ) {
        // NOTE: The polar method of the Box-Muller transformation produced really bad
        //       artifacts.  That is why this code does not use it.

        // Do the box-muller transformation
        // NOTE: This method produces significantly better results than the polar form of the transformation
        float x = rng(), y = rng();
        double coefficient = sqrt( -2 * log( x ) );
        return vector2f( float( coefficient * cos( 2 * M_PI * y ) ), float( coefficient * sin( 2 * M_PI * y ) ) );
    }

    // TODO: We should use the boost random number generator for high quality and fast random numbers
    static vector2f from_unit_random() {
        vector2f result = vector2f::from_random_gaussian();
        // Normalize to unit distance
        result.normalize();
        return result;
    }

    static float triangle_area( const vector2f& a, const vector2f& b, const vector2f& c ) {
        float ad = vector2f::distance( a, b ), bd = vector2f::distance( b, c ), cd = vector2f::distance( c, a );
        float s = ( ad + bd + cd ) / 2;
        return sqrtf( s * ( s - ad ) * ( s - bd ) * ( s - cd ) );
    }

    static float triangle_curvature( const vector2f& a, const vector2f& b, const vector2f& c ) {
        double ad = sqrt( (double)vector2f::distance_squared( a, b ) ),
               bd = sqrt( (double)vector2f::distance_squared( b, c ) ),
               cd = sqrt( (double)vector2f::distance_squared( c, a ) );
        double s = ( ad + bd + cd ) / 2;
        double area = sqrt( s * ( s - ad ) * ( s - bd ) * ( s - cd ) );
        double absCurvature = 4 * area / ( ad * bd * cd );
        double signTest = ( (double)a.x - b.x ) * ( (double)c.y - b.y ) - ( (double)c.x - b.x ) * ( (double)a.y - b.y );
        if( signTest > 0 )
            return float( -absCurvature );
        else
            return float( absCurvature );
    }

    // Returns the z magnitude of the cross product, which is |a| * |b| * sin(theta)

    static vector2f xaxis() { return vector2f( 1, 0 ); }

    static vector2f yaxis() { return vector2f( 0, 1 ); }

// exclude this for visual studio 6 :P
#if !defined( _MSC_VER ) || _MSC_VER > 1200
    float max_abs_component() const { return ( std::max )( fabs( x ), fabs( y ) ); }
#endif

    void normalize() {
        float magnitudeSquared = get_magnitude_squared();
        if( magnitudeSquared < 0.00000000001f || ( magnitudeSquared > 0.9999999f && magnitudeSquared < 1.0000001f ) )
            return;
        float length = sqrtf( magnitudeSquared );
        x /= length;
        y /= length;
    }

    vector2f to_normalized() const {
        vector2f result = *this;
        result.normalize();
        return result;
    }

    void from_frantic_to_fusion() {
        x -= .5f;
        y -= .5f;
    }

    void from_frantic_to_fusion4( float heightOfImage ) {
        x -= .5f;
        y = heightOfImage - y - .5f;
    }

    void from_fusion_to_frantic() {
        x += .5f;
        y += .5f;
    }

    void from_fusion4_to_frantic( float heightOfImage ) {
        x += .5f;
        y = heightOfImage - y - .5f;
    }

    static vector2f normalize( vector2f v ) {
        v.normalize();
        return v;
    }

    static vector2f abs( vector2f v ) { return vector2f( fabs( v.x ), fabs( v.y ) ); }
};

inline vector2f vector2f::from_random_gaussian() {
    // Use the non-polar form of the box-muller transformation
    float x = (float)rand() / RAND_MAX, y = (float)rand() / RAND_MAX;
    double coefficient = sqrt( -2 * log( x ) );
    return vector2f( float( coefficient * cos( 2 * M_PI * y ) ), float( coefficient * sin( 2 * M_PI * y ) ) );
}

// The compiler wasn't finding these particular operators in the vector2t template class.
inline vector2f operator*( float val, const vector2f& b ) { return vector2f( val * b.x, val * b.y ); }

inline vector2f operator*( const vector2f& a, float val ) { return vector2f( a.x * val, a.y * val ); }

inline vector2f operator/( const vector2f& a, float val ) { return vector2f( a.x / val, a.y / val ); }

// LessThan operator used for specialized algorithms regarding sorted points.
inline bool operator<( const vector2f& a, const vector2f& b ) {
    return ( ( a.x < b.x ) || ( ( a.x == b.x ) && ( a.y < b.y ) ) );
}

} // namespace graphics2d
} // namespace frantic
