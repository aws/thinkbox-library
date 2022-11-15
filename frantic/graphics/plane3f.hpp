// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <iostream>

#include <assert.h>

#include <frantic/graphics/ray3f.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace graphics {

template <class FloatType>
class plane3t {
    vector3t<FloatType> m_normal;
    FloatType m_constant;

  public:
    typedef FloatType float_type;

    plane3t( float_type X, float_type Y, float_type Z, float_type K )
        : m_normal( X, Y, Z )
        , m_constant( K ) {}

    plane3t( const vector3t<float_type>& N, float_type K )
        : m_normal( N )
        , m_constant( K ) {}

    plane3t() {
        m_normal = vector3t<float_type>( 0, 0, 1 );
        m_constant = 0;
    }

    // Constructs a plane with normal N, passing through P
    static plane3t from_normal_and_point( vector3t<float_type> N, const vector3t<float_type>& P ) {
        // Normalize N, make it a default Z axis if it's a 0 vector
        float_type length = N.get_magnitude();
        if( length == 0 )
            N = vector3t<float_type>::from_zaxis();
        else
            N /= length;

        return plane3t( N, vector3t<float_type>::dot( N, P ) );
    }

    static plane3t from_triangle( const vector3t<float_type>& a, const vector3t<float_type>& b,
                                  const vector3t<float_type>& c ) {
        return plane3t::from_normal_and_point( vector3t<float_type>::cross( c - b, a - b ), b );
    }

    static plane3t from_triangle( const vector3t<float_type>* verts ) {
        return plane3t::from_triangle( verts[0], verts[1], verts[2] );
    }

    /**
     *  Given a triangle with vertices onPlaneA, onPlaneB, and offPlaneC,
     * create a plane which touches onPlaneA and onPlaneB, and with a normal
     * that is orthogonal to the triangle's normal.
     *
     * @param onPlaneA a point on the plane.
     * @param onPlaneB a second point on the plane.
     * @param offPlaneC a point which is not on the plane, used to determine
     *		the face's normal.
     * @return a plane which touches the points onPlaneA and onPlaneB, and
     *		which has a normal orthogonal to the normal of the triangle
     *		(onPlaneA, onPlaneB, offPlaneC).
     */
    static plane3t from_triangle_edge( const vector3t<float_type>& onPlaneA, const vector3t<float_type>& onPlaneB,
                                       const vector3t<float_type>& offPlaneC ) {
        const vector3t<float_type> edgeAB = onPlaneA - onPlaneB;
        const vector3t<float_type> edgeCB = offPlaneC - onPlaneB;

        vector3t<float_type> planeNormal;

        const float_type magnitudeSquaredAB = edgeAB.get_magnitude_squared();
        if( magnitudeSquaredAB > 0 ) {
            const vector3t<float_type> projectedCBonAB =
                edgeAB * ( vector3t<float_type>::dot( edgeAB, edgeCB ) / edgeAB.get_magnitude_squared() );
            planeNormal = edgeCB - projectedCBonAB;
        } else {
            planeNormal = edgeCB;
        }

        return plane3t::from_normal_and_point( planeNormal, onPlaneA );
    }

    static plane3t xaxis() { return plane3t( vector3t<float_type>::from_xaxis(), 0 ); }

    static plane3t yaxis() { return plane3t( vector3t<float_type>::from_yaxis(), 0 ); }

    static plane3t zaxis() { return plane3t( vector3t<float_type>::from_zaxis(), 0 ); }

    vector3t<float_type> normal() const { return m_normal; }

    float_type constant() const { return m_constant; }

    // This returns the perpendicular distance from the point P to the plane.
    // The regular convention is that positive values are "outside" and negative values are "inside", same as with
    // our level sets.
    float_type get_signed_distance_to_plane( const vector3t<float_type>& P ) const {
        return vector3t<float_type>::dot( P, m_normal ) - m_constant;
    }

    double get_signed_distance_to_plane_double( const vector3t<float_type>& P ) const {
        return vector3t<float_type>::dot_double( P, m_normal ) - static_cast<double>( m_constant );
    }

    int get_sign( const vector3t<float_type>& P ) const {
        using namespace std;

        float_type PdotNormal = vector3t<float_type>::dot( P, m_normal );
        float_type d = PdotNormal - m_constant;

        // The reason for using this as the relative error bounds is that the subtraction between these two quantities
        // to compute d is going to be the most common cause of precision loss.  Thus, the greatest absolute value of
        // the two values being subtracted will roughly be the "magnitude" for error bounds purposes. It also turned out
        // that there were cases where P and m_normal were much larger values than PdotNormal and m_constant, so also
        // needed to be included in this error term.
        float_type relativeErrorBounds =
            0.000001f * ( std::max )( ( std::max )( P.max_abs_component(), abs( PdotNormal ) ),
                                      ( std::max )( m_normal.max_abs_component(), abs( m_constant ) ) );

        if( d < -relativeErrorBounds )
            return -1;
        else if( d > relativeErrorBounds )
            return 1;
        else
            return 0;
    }

    static int get_sign_exact( const vector3t<float_type>& a, const vector3t<float_type>& b,
                               const vector3t<float_type>& c, const vector3t<float_type>& d );

    vector3t<float_type> project_onto_plane( const vector3t<float_type>& P ) const {
        return P + ( m_constant - vector3t<float_type>::dot( P, m_normal ) ) * m_normal;
    }

    // returns true if the line segment from a to b intersects the plane
    bool is_intersecting( const vector3t<float_type>& a, const vector3t<float_type>& b ) const {
        int s0 = get_sign( a );
        int s1 = get_sign( b );

        return s0 == 0 || s1 == 0 || ( s0 == -1 && s1 == 1 ) || ( s0 == 1 && s1 == -1 );
    }

    bool is_collinear( const vector3t<float_type>& a, const vector3t<float_type>& b ) const {
        int s0 = get_sign( a );
        int s1 = get_sign( b );
        return s0 == 0 && s1 == 0;
    }

    bool is_coplanar( const vector3t<float_type>& p ) const { return get_sign( p ) == 0; }

    // This gets the intersection of the given line segment and the plane.  Fills outIntersectionPoint and returns true
    // when the point exists and is unique, otherwise returns false.
    bool get_intersection( const vector3t<float_type>& start, const vector3t<float_type>& end,
                           vector3t<float_type>& outIntersectionPoint ) const {
        float_type aDist = get_signed_distance_to_plane( start ), bDist = get_signed_distance_to_plane( end );
        if( aDist > 0 ) {
            if( bDist <= 0 ) {
                float_type normalizationFactor = 1.f / ( aDist - bDist );
                outIntersectionPoint = ( aDist * normalizationFactor ) * end - ( bDist * normalizationFactor ) * start;
                return true;
            }
        } else {
            if( bDist > 0 ) {
                float_type normalizationFactor = 1.f / ( bDist - aDist );
                outIntersectionPoint = ( bDist * normalizationFactor ) * start - ( aDist * normalizationFactor ) * end;
                return true;
            }
        }
        return false;
    }

    double get_distance_to_intersection( const ray3t<FloatType>& ray ) const {
        using namespace std;

        double normalDotRayDirection = vector3t<float_type>::dot_double( m_normal, ray.direction() );
        // TODO: Probably shouldn't be using an absolute epsilon value!
        if( abs( normalDotRayDirection ) < 0.0001f ) {
            if( is_coplanar( ray.origin() ) )
                return 0;
            else
                return -1;
        }
        double t =
            ( (double)m_constant - vector3t<float_type>::dot_double( m_normal, ray.origin() ) ) / normalDotRayDirection;
        return t >= 0 ? t : -1;
    }

    double get_distance_to_intersection( const ray3t<float_type>& ray, bool& outIsCoplanar ) {
        using namespace std;

        double normalDotRayDirection = vector3t<float_type>::dot_double( m_normal, ray.direction() );
        // TODO: Probably shouldn't be using an absolute epsilon value!
        outIsCoplanar = false;
        if( abs( normalDotRayDirection ) < 0.0001f ) {
            if( is_coplanar( ray.origin() ) ) {
                outIsCoplanar = true;
                return 0;
            } else
                return -1;
        }
        double t =
            ( (double)m_constant - vector3t<float_type>::dot_double( m_normal, ray.origin() ) ) / normalDotRayDirection;
        return t >= 0 ? t : -1;
    }

    // same as get_distance_to_intersection(ray), only limits the search to the interval [t0, t1).
    double get_distance_to_intersection( const ray3t<float_type>& ray, float_type t0, float_type t1 ) const {
        double distance = get_distance_to_intersection( ray );
        if( distance >= t0 && distance < t1 )
            return distance;
        else
            return -1.0;
    }

    // Clip the given ray segment to only the parts that are on the -ve side of the plane.
    // Returns false if the ray is totally culled, otherwise true with the new segment parameters set in t0/t1.
    bool clip_ray_segment_to_negative_region( const ray3t<float_type>& ray, float_type& t0, float_type& t1 ) const {
        // See which side of the plane the ray starts on. +ve is above, -ve is below
        float_type originDistance = get_signed_distance_to_plane( ray.at( t0 ) );

        // Ray intersects plane surface. see if we're above plane or below
        float_type hitDistance = (float_type)get_distance_to_intersection( ray, t0, t1 );
        if( hitDistance > 0.0f ) {
            if( originDistance >= 0.0f ) // We're above, so t0 gets pushed to the plane surface
                t0 = hitDistance;
            else // We're below, so t1 gets clipped to the plane surface
                t1 = hitDistance;
        } else {
            // We don't intersect the water plane, so if we're above the plane the entire
            // segment gets culled
            if( originDistance >= 0.0f )
                return false;
        }
        return true;
    }

    void get_basis( vector3t<float_type>& u, vector3t<float_type>& v ) const {
        vector3t<float_type> P = project_onto_plane( vector3t<float_type>( 0 ) );

        u = project_onto_plane( P + vector3t<float_type>::from_xaxis() ) - P;
        u = vector3t<float_type>::max_magnitude( u, project_onto_plane( P + vector3t<float_type>::from_yaxis() ) - P );
        u = vector3t<float_type>::max_magnitude( u, project_onto_plane( P + vector3t<float_type>::from_zaxis() ) - P );

        v = vector3t<float_type>::cross( u, m_normal );

        u.normalize();
        v.normalize();
    }

    /**
     * Returns an affine transformation that will take points from this plane into the xy plane.
     * Note: this is not a projection, this will preserve distances between points
     *
     * @param originOnPlane  The projection of this point onto the plane is mapped to the origin in the resulting
     *                       transform.
     * @return a 1-1 transformation such that given any point on plane, it maps to a unique point with z = 0, preserving
     *         euclidean distance
     */
    transform4t<float_type>
    get_planar_transform( const vector3t<float_type>& originOnPlane = vector3t<float_type>( 0 ) ) {
        vector3t<float_type> u;
        vector3t<float_type> v;

        get_basis( u, v );

        // form a rotation matrix to restore a z-up orientation
        transform4t<float_type> xform( u, v, m_normal );

        // add in a translation down to the z=0 plane
        return xform.to_transpose() * transform4t<float_type>::from_translation( -project_onto_plane( originOnPlane ) );
    }

    //////////////////
    // Operators
    //////////////////

    plane3t operator-() { return plane3t( -m_normal, -m_constant ); }

    void flip() {
        m_normal = -m_normal;
        m_constant = -m_constant;
    }

    std::string str() const;
};

/**
 * Computes the 6 supporting planes of a bounding box. The generated planes
 * will have normals pointing away from the bounding box.
 * @remark The reason for this function's peculiar inclusion is due to
 *   file inclusion dependencies.
 * @param[in] box the bounding box to use
 * @param[out] outPlanes where to put the resulting planes
 * @remark the planes will be placed into the ouput according to the
 *  indices specified by the `cube_face::default_cube_face` enum
 */
template <typename FloatType>
void get_bounding_box_planes( const boundbox3t<FloatType>& box, plane3t<FloatType>* outPlanes ) {
    outPlanes[int( cube_face::CF_X_NEG )] =
        plane3t<FloatType>( FloatType( -1.0 ), FloatType( 0.0 ), FloatType( 0.0 ), -box.minimum().x );
    outPlanes[int( cube_face::CF_X_POS )] =
        plane3t<FloatType>( FloatType( 1.0 ), FloatType( 0.0 ), FloatType( 0.0 ), box.maximum().x );
    outPlanes[int( cube_face::CF_Y_NEG )] =
        plane3t<FloatType>( FloatType( 0.0 ), FloatType( -1.0 ), FloatType( 0.0 ), -box.minimum().y );
    outPlanes[int( cube_face::CF_Y_POS )] =
        plane3t<FloatType>( FloatType( 0.0 ), FloatType( 1.0 ), FloatType( 0.0 ), box.maximum().y );
    outPlanes[int( cube_face::CF_Z_NEG )] =
        plane3t<FloatType>( FloatType( 0.0 ), FloatType( 0.0 ), FloatType( -1.0 ), -box.minimum().z );
    outPlanes[int( cube_face::CF_Z_POS )] =
        plane3t<FloatType>( FloatType( 0.0 ), FloatType( 0.0 ), FloatType( 1.0 ), box.maximum().z );
}

} // namespace graphics
} // namespace frantic

namespace frantic {
namespace graphics {

template <class FloatType>
inline std::ostream& operator<<( std::ostream& out, const plane3t<FloatType>& v ) {
    vector3t<FloatType> N = v.normal();
    out << "(plane3f " << N.x << ", " << N.y << ", " << N.z << ", " << v.constant() << " )";
    return out;
}

template <class FloatType>
inline std::string plane3t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

typedef plane3t<float> plane3f;
typedef plane3t<double> plane3fd;
} // namespace graphics
} // namespace frantic
