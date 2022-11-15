// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <iostream>
#include <sstream>

#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/size2.hpp>
#include <frantic/graphics2d/size2f.hpp>
#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

namespace frantic {
namespace graphics {

// This should never be done, importing a namespace into another is a bad idea.
// using namespace graphics2d;

template <class FloatType>
class spherical_coords_3t {
  public:
    typedef FloatType float_type;
    typedef vector3t<float_type> vector3f_type;

    // It was found that double accuracy is necessary to resolve the difference between the coordinates
    // 320,240 and adjacent pixels in a 640x480 perspective projection.
    // Theta is in [0,360), Phi is in [0,180]
    float_type theta, phi;

    /// Constructors

    spherical_coords_3t() {
        theta = 0;
        phi = 0;
    }

    spherical_coords_3t( float_type thetaValue, float_type phiValue ) {
        theta = thetaValue;
        phi = phiValue;
    }

    // Constructs the spherical coordinates, assuming we're
    // looking towards positive z, y is up, and x is rightwards.
    // This is looking towards the north pole of the sphere, where phi is 0.
    explicit spherical_coords_3t( const vector3f_type& v ) {
        theta = atan2( v.y, v.x );
        phi = acos( v.z / sqrt( v.x * v.x + v.y * v.y + v.z * v.z ) );
    }

    static spherical_coords_3t from_longlat( float_type longitude, float_type latitude ) {
        float_type theta = fmod( longitude, 1 ) * 2 * float_type( M_PI );
        float_type phi = fmod( latitude, 1 ) * float_type( M_PI );
        return spherical_coords_3t( theta, phi );
    }

    static spherical_coords_3t from_random() { return spherical_coords_3t( vector3f_type::from_random_gaussian() ); }

    vector3f to_vector3f() const {
        return vector3f( float( sin( phi ) * cos( theta ) ), float( sin( phi ) * sin( theta ) ), float( cos( phi ) ) );
    }

    vector3f_type to_vector3t() const {
        return vector3f_type( float_type( sin( phi ) * cos( theta ) ), float_type( sin( phi ) * sin( theta ) ),
                              float_type( cos( phi ) ) );
    }

    frantic::graphics2d::vector2f to_longlat() const {
        frantic::graphics2d::vector2f result( (float)fmod( theta / float_type( ( 2 * M_PI ) ), 1 ),
                                              (float)fmod( phi / float_type( M_PI ), 1 ) );
        // fmod returns a result that's the same sign as the original, not quite what we want.
        if( result.x < 0 )
            result.x += 1;
        if( result.y < 0 )
            result.y += 1;
        return result;
    }

    frantic::graphics2d::vector2f to_fisheye( float fov ) const {
        // TODO: Why do we convert to vector3f?  It should just be a direct manipulation of the
        //       theta and phi values.
        vector3f_type v = to_vector3t();

        float_type phi = acos( -v.z );
        float_type theta = atan2( -v.y, v.x );

        // Fov used to be in degrees
        // return vector2f((float)(cos(theta) * phi * 2 / math::degrees_to_radians(fov)), (float)(sin(theta) * phi * 2 /
        // (fov * M_PI / 180)));

        return frantic::graphics2d::vector2f( (float)( cos( theta ) * phi * 2 / fov ),
                                              (float)( sin( theta ) * phi * 2 / fov ) );
    }

    static float_type angle_between( spherical_coords_3t a, spherical_coords_3t b ) {
        if( a == b )
            return 0.f;

        vector3f_type va = a.to_vector3f(), vb = b.to_vector3f();
        return float( acos( vector3f_type::dot( va, vb ) ) );
    }

    std::string str() const;

    bool operator==( spherical_coords_3t const& b ) const { return theta == b.theta && phi == b.phi; }
};

template <class FloatType>
inline std::ostream& operator<<( std::ostream& out, const spherical_coords_3t<FloatType>& v ) {
    out << "(spherical theta=" << v.theta << " phi=" << v.phi << " )";
    return out;
}

template <class FloatType>
inline std::string spherical_coords_3t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

typedef spherical_coords_3t<float> spherical_coords_3f;
typedef spherical_coords_3t<float> spherical_coords_3fd;

// Retain backwards compatibility with existing users of this class
typedef spherical_coords_3t<double> spherical_coords;

} // namespace graphics
} // namespace frantic
