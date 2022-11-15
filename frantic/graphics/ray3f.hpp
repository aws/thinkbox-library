// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#if defined( __GNUC__ ) && __GNUC__ < 3
#include <limits.h>
#else
#include <limits>
#endif
#include <iostream>

#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace graphics {

/**
 * This class represents a ray.  The direction does not need to be normalized, however note that distances along the ray
 * are in multiples of the length of the m_direction vector.  The reason for not normalizing it is so that we transform
 * a ray in any way, including assymmetric skews and scales, then compute distances to geometry intersections within the
 * transformed space.
 */
template <typename FloatType>
class ray3t {
  public:
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;
    typedef size3t<FloatType> size3f_type;

  private:
    vector3f_type m_origin, m_direction;

  public:
    ///////////////
    // Constructors
    ///////////////

    ray3t()
        : m_origin( 0 )
        , m_direction( vector3f_type::from_xaxis() ) {}

    ray3t( const vector3f_type& origin, const vector3f_type& direction )
        : m_origin( origin )
        , m_direction( direction ) {}

    static ray3t from_random( boundbox3t<float_type>& box ) {
        return ray3t( box.random_vector(), vector3f_type::from_unit_random() );
    }

    // TODO: We should use the boost random number generator for high quality and fast random numbers
    static ray3t from_random_towards_box( boundbox3t<float_type>& box ) {
        if( box.is_empty() )
            throw std::runtime_error( "ray3t::from_random_towards_box: Cannot get a random ray from an empty box." );
        // Get a vector outside the box
        vector3f_type start = 10 * ( box.random_vector() - box.center() ) + box.center();
        while( box.contains( start ) ) {
            start = 10 * ( box.random_vector() - box.center() ) + box.center();
        }
        // Get a vector inside the box
        vector3f_type end = box.random_vector();
        return ray3t( start, ( end - start ).to_normalized() );
    }

    void set( const vector3f_type& origin, const vector3f_type& direction ) {
        m_origin = origin;
        m_direction = direction;
    }

    void set_origin( const vector3f_type& origin ) { m_origin = origin; }

    void set_direction( const vector3f_type& direction ) { m_direction = direction; }

    vector3f_type origin() const { return m_origin; }

    vector3f_type direction() const { return m_direction; }

    static ray3t from_line_segment( const vector3f_type& start, const vector3f_type& end ) {
        return ray3t( start, end - start );
    }

    static size3t<FloatType> from_cube_side_length( float sideLength ) {
        return size3t<FloatType>( sideLength, sideLength, sideLength );
    }

    vector3f_type at( float d ) const {
        return vector3f_type( m_origin.x + d * m_direction.x, m_origin.y + d * m_direction.y,
                              m_origin.z + d * m_direction.z );
    }

    vector3f_type at( double d ) const {
        return vector3f_type( static_cast<float>( m_origin.x + d * m_direction.x ),
                              static_cast<float>( m_origin.y + d * m_direction.y ),
                              static_cast<float>( m_origin.z + d * m_direction.z ) );
    }

    vector3f_type operator[]( float d ) const { return at( d ); }

    vector3f_type operator[]( double d ) const { return at( d ); }

    bool contains_point( vector3f_type p, float_type tolerance = 0.00001f ) const {
        // TODO: Make tolerance depend on float_type
        vector3f_type delta = p - m_origin;

        float_type dTolerance = ( std::max )( delta.get_magnitude_squared(), p.get_magnitude_squared() );
        ;

        delta -= vector3f_type::project( delta, m_direction );
        float_type dSquared = delta.get_magnitude_squared();
        if( dTolerance < 0.0000001f ) // Handle the case where p is approximately 0
            dTolerance = 0.0000001f;
        else
            dTolerance *= tolerance * tolerance;
        // Check whether the relative error is within tolerance
        return dSquared < dTolerance;
    }

    /**
     * Return the parametric location of the nearest point to `p`
     * along the current ray.
     * @remark this will also include the negative parameter space
     *
     * @param p the point to project onto the ray
     * @return the t-value along the ray of the nearest point
     */
    float_type nearest_point_parameter( const vector3f_type& p ) const {
        vector3f_type offset = p - m_origin;
        return vector3f_type::dot( m_direction, offset ) / m_direction.get_magnitude_squared();
    }

    /**
     * Get the nearest point to `p` along the current ray.
     * @remark this will also include the negative direction of the ray
     *
     * @param p the point ot project onto the ray
     * @return the nearest point in terms of perpendicular distance
     */
    vector3f_type nearest_point( const vector3f_type& p ) const { return at( nearest_point_parameter( p ) ); }

    /**
     * Get the time at which the two given rays are closest to each other. Solutions may be negative.
     *
     * @return The parametric equation parameters which result in a minimal distance between the points if there is a
     *         unique solution, or (NaN, NaN) if there are infinite solutions (rays are parallel or antiparallel).
     */
    std::pair<float_type, float_type>
    nearest_point_parameters( const frantic::graphics::ray3t<FloatType>& other ) const {
        // This involves solving the following system of two equations:
        // For o_1 = m_origin, o_2 = other.origin(), d_1 = m_direction, d_2 = other.direction()
        // [ (o_1 - o_2) dot d_1 ] = [  (d_1 dot d_1) t1 - (d_1 dot d_2) t2 ]
        // [ (o_2 - o_1) dot d_2 ] = [ -(d_2 dot d_1) t1 + (d_2 dot d_2) t2 ]

        const vector3f_type deltaStart = m_origin - other.origin();
        const float_type const1 = vector3f_type::dot( -deltaStart, m_direction );
        const float_type const2 = vector3f_type::dot( deltaStart, other.direction() );
        const float_type ray1ray1 = vector3f_type::dot( m_direction, m_direction );
        const float_type ray1ray2 = vector3f_type::dot( m_direction, other.direction() );
        const float_type ray2ray2 = vector3f_type::dot( other.direction(), other.direction() );

        const float_type a11 = ray1ray1;
        const float_type a12 = -ray1ray2;
        const float_type a21 = -ray1ray2;
        const float_type a22 = ray2ray2;

        const float_type determinate = ( a11 * a22 - a12 * a21 );

        const float_type t1 = ( const1 * a22 - const2 * a12 );
        const float_type t2 = ( -const1 * a21 + const2 * a11 );

        return std::pair<float_type, float_type>( t1 / determinate, t2 / determinate );
    }

    float_type distance_from_point( vector3f_type p ) const {
        vector3f_type delta = p - m_origin;
        delta -= vector3f_type::dot( delta, m_direction ) * m_direction;
        return delta.get_magnitude();
    }

    float_type distance_squared_from_point( vector3f_type p ) const {
        vector3f_type delta = p - m_origin;
        delta -= vector3f_type::dot( delta, m_direction ) * m_direction;
        return delta.get_magnitude_squared();
    }

    std::string str() const;

    // Moves the current ray to the far edge of the box.  Returns false if no movement occurred
    bool move_to_far_box_intersection( const boundbox3t<FloatType>& box ) {
        double nearIntersection = 0, farIntersection = 0;
        if( get_intersection( box, nearIntersection, farIntersection ) ) {
            if( farIntersection > 0 ) {
                m_origin = at( farIntersection );
                return true;
            } else
                return false;
        } else
            return false;
    }

    // Moves the current ray to the near edge of the box.  Returns false if no movement occurred
    bool move_to_near_box_intersection( const boundbox3t<FloatType>& box ) {
        double nearIntersection = 0, farIntersection = 0;
        if( get_intersection( box, nearIntersection, farIntersection ) ) {
            if( nearIntersection > 0 ) {
                m_origin = at( nearIntersection );
                return true;
            } else
                return false;
        } else
            return false;
    }

    // Moves the current ray to the next edge of the box.  Returns false if no movement occurred
    bool move_to_next_box_intersection( const boundbox3f& box ) {
        if( box.contains( m_origin ) )
            return move_to_far_box_intersection( box );
        else
            return move_to_near_box_intersection( box );
    }

    bool get_intersection( const vector3f_type& boxMin, const vector3f_type& boxMax, double& outEntryLocation,
                           double& outExitLocation ) const {
        vector3f_type entryNormal, exitNormal;
        return intersect_with_box( boxMin, boxMax, outEntryLocation, outExitLocation, &entryNormal, &exitNormal );
    }

    bool get_intersection( const boundbox3t<FloatType>& box, double& outEntryLocation, double& outExitLocation ) const {
        vector3f_type entryNormal, exitNormal;
        return intersect_with_box( box.minimum(), box.maximum(), outEntryLocation, outExitLocation, &entryNormal,
                                   &exitNormal );
    }

    bool intersect_with_box( const boundbox3t<FloatType>& box, double& outEntryLocation, double& outExitLocation,
                             vector3f_type* outEntryNormal, vector3f_type* outExitNormal ) const {
        return intersect_with_box( box.minimum(), box.maximum(), outEntryLocation, outExitLocation, outEntryNormal,
                                   outExitNormal );
    }

    // Intersects the ray with the specified triangle, returning the intersection point
    bool intersect_with_triangle( const vector3f_type& vert0, const vector3f_type& vert1, const vector3f_type& vert2,
                                  double& outIntersection, vector3f_type& outBaryCoord ) const {
        // From "Fast, Minimum Storage Ray/Triangle Intersection" by Tomas Akenine-Moller and Ben Trumbore
        // This method uses Cramer's rule to solve for the ray intersection.
        // Suppose vert0 is at the origin, to simplify the discussion.  Set E1 = vert1, E2 = vert2, B = ray origin, and
        // D = ray direction Then, we want to solve the following 3x3 linear system: (E1 E2 D)*(u) = B
        //           (v)
        //           (t)
        // In this system, u and v are the barycentric coordinates, and t is the distance along the ray direction D.
        vector3f_type edge1 = vert1 - vert0, edge2 = vert2 - vert0;
        vector3f_type pVec = vector3f_type::cross( m_direction, edge2 );
        double det = vector3f_type::dot_double( edge1, pVec );

        // TODO: what value to use for this epsilon, that's not hardcoded like this?
        if( det > -0.000001 && det < 0.000001 )
            return false;

        double invDet = 1 / det;
        vector3f_type tVec = m_origin - vert0;

        outBaryCoord.x = (float_type)( vector3f_type::dot_double( tVec, pVec ) * invDet );
        if( outBaryCoord.x < 0 || outBaryCoord.x > 1 )
            return false;

        vector3f_type qVec = vector3f_type::cross( tVec, edge1 );

        outBaryCoord.y = (float_type)( vector3f_type::dot_double( m_direction, qVec ) * invDet );
        if( outBaryCoord.y < 0 || outBaryCoord.x + outBaryCoord.y > 1 )
            return false;

        outIntersection = vector3f_type::dot_double( edge2, qVec ) * invDet;

        outBaryCoord.z = 1.f - outBaryCoord.x - outBaryCoord.y;

        return true;
    }

    // Intersects the ray with the specified triangle, returning the intersection point
    bool intersect_with_triangle( const vector3f_type& vert0, const vector3f_type& vert1, const vector3f_type& vert2,
                                  vector3f_type& outIntersection, vector3f_type& outBaryCoord ) const {
        double t = 0;
        if( intersect_with_triangle( vert0, vert1, vert2, t, outBaryCoord ) ) {
            outIntersection = at( t );
            return true;
        } else {
            return false;
        }
    }

    // Intersects the ray with the box.  Uses the range [0, infinity) for t.
    // ray.At(t0) is the first intersection found, with the given normal
    // ray.At(t1) is the second intersection
    bool intersect_with_box( const vector3f_type& boxMin, const vector3f_type& boxMax, double& outEntryLocation,
                             double& outExitLocation, vector3f_type* outEntryNormal,
                             vector3f_type* outExitNormal ) const;

    bool is_intersecting_box_surface( const boundbox3t<FloatType>& box, float_type start, float_type end ) const {
        bool startInside = box.contains_as_open_set( at( start ) ), endInside = box.contains_as_open_set( at( end ) );
        if( startInside && endInside )
            return false;
        if( startInside || endInside )
            return true;
        return is_intersecting_box_volume( box, start, end );
    }

    bool is_intersecting_box_volume( const boundbox3t<FloatType>& box, float_type start, float_type end ) const;

    /**
     * This function will clamp tMin and tMax, such that both this->at(tMin) and this->at(tMax)
     * are within the given box.
     * @param box The box to clamp to
     * @param tMin The input/output minimum ray value. After this function this->at(tMin) will be in 'box'
     * @param tMax The input/output maximum ray value. After this function this->at(tMax) will be in 'box'
     * @return False if the ray does not intersect the box, or if it does so outside of the range [tMin,tMax]
     */
    bool clamp_to_box( const boundbox3t<FloatType>& box, double& tMin, double& tMax ) const {
        if( m_direction.x > 0 ) {
            double xTMin = ( box.minimum().x - m_origin.x ) / m_direction.x;
            if( xTMin > tMin ) {
                if( xTMin > tMax )
                    return false;
                tMin = xTMin;
            }

            double xTMax = ( box.maximum().x - m_origin.x ) / m_direction.x;
            if( xTMax < tMax ) {
                if( xTMax < tMin )
                    return false;
                tMax = xTMax;
            }
        } else if( m_direction.x < 0 ) {
            double xTMin = ( box.maximum().x - m_origin.x ) / m_direction.x;
            if( xTMin > tMin ) {
                if( xTMin > tMax )
                    return false;
                tMin = xTMin;
            }

            double xTMax = ( box.minimum().x - m_origin.x ) / m_direction.x;
            if( xTMax < tMax ) {
                if( xTMax < tMin )
                    return false;
                tMax = xTMax;
            }
        } else if( m_origin.x < box.minimum().x || m_origin.x > box.maximum().x )
            return false;

        if( m_direction.y > 0 ) {
            double yTMin = ( box.minimum().y - m_origin.y ) / m_direction.y;
            if( yTMin > tMin ) {
                if( yTMin > tMax )
                    return false;
                tMin = yTMin;
            }

            double yTMax = ( box.maximum().y - m_origin.y ) / m_direction.y;
            if( yTMax < tMax ) {
                if( yTMax < tMin )
                    return false;
                tMax = yTMax;
            }
        } else if( m_direction.y < 0 ) {
            double yTMin = ( box.maximum().y - m_origin.y ) / m_direction.y;
            if( yTMin > tMin ) {
                if( yTMin > tMax )
                    return false;
                tMin = yTMin;
            }

            double yTMax = ( box.minimum().y - m_origin.y ) / m_direction.y;
            if( yTMax < tMax ) {
                if( yTMax < tMin )
                    return false;
                tMax = yTMax;
            }
        } else if( m_origin.y < box.minimum().y || m_origin.y > box.maximum().y )
            return false;

        if( m_direction.z > 0 ) {
            double zTMin = ( box.minimum().z - m_origin.z ) / m_direction.z;
            if( zTMin > tMin ) {
                if( zTMin > tMax )
                    return false;
                tMin = zTMin;
            }

            double zTMax = ( box.maximum().z - m_origin.z ) / m_direction.z;
            if( zTMax < tMax ) {
                if( zTMax < tMin )
                    return false;
                tMax = zTMax;
            }
        } else if( m_direction.z < 0 ) {
            double zTMin = ( box.maximum().z - m_origin.z ) / m_direction.z;
            if( zTMin > tMin ) {
                if( zTMin > tMax )
                    return false;
                tMin = zTMin;
            }

            double zTMax = ( box.minimum().z - m_origin.z ) / m_direction.z;
            if( zTMax < tMax ) {
                if( zTMax < tMin )
                    return false;
                tMax = zTMax;
            }
        } else if( m_origin.z < box.minimum().z || m_origin.z > box.maximum().z )
            return false;

        return true;
    }

    /////////////////////////
    // Operators
    /////////////////////////

    // Jitters the ray by a very very tiny amount.  Used as a stop-gap solution for when
    // the voxel edge stepper screws up (which is rare, but even 1 in 100 million rays is bad).
    // The cases where it screws up are detected, however, and throw an exception.  Thus we
    // can jitter and try again.
    void micro_jitter() {
        float originMagnitude = ( std::max )( m_origin.get_magnitude(), 1.f );
        m_origin += 0.00001f * originMagnitude * vector3f_type::from_random_gaussian();
        m_direction += 0.00001f * vector3f_type::from_random_gaussian();
    }

    ray3t& operator+=( const vector3f_type& a ) {
        m_origin += a;
        return *this;
    }

    ray3t& operator-=( const vector3f_type& a ) {
        m_origin -= a;
        return *this;
    }

    bool operator==( const ray3t& rhs ) const {
        return ( m_origin == rhs.m_origin ) && ( m_direction == rhs.m_direction );
    }
};

/////////////////////////
// Operators
/////////////////////////

template <class FloatType>
std::ostream& operator<<( std::ostream& out, const ray3t<FloatType>& r ) {
    out << "(ray " << r.origin() << " " << r.direction() << " )";
    return out;
}

template <class FloatType>
std::string ray3t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

typedef ray3t<float> ray3f;
typedef ray3t<double> ray3fd;

namespace detail {
// A mechanism to get the more complicated/slow to compile code into the .cpp file
bool intersect_with_box( const ray3f& ray, const vector3f& boxMin, const vector3f& boxMax, double& outEntryLocation,
                         double& outExitLocation, vector3f* outEntryNormal, vector3f* outExitNormal );
bool intersect_with_box( const ray3fd& ray, const vector3fd& boxMin, const vector3fd& boxMax, double& outEntryLocation,
                         double& outExitLocation, vector3fd* outEntryNormal, vector3fd* outExitNormal );
bool is_intersecting_box_volume( const ray3f& ray, const boundbox3f& box, float start, float end );
bool is_intersecting_box_volume( const ray3fd& ray, const boundbox3fd& box, double start, double end );
} // namespace detail

template <class FloatType>
inline bool ray3t<FloatType>::intersect_with_box( const vector3f_type& boxMin, const vector3f_type& boxMax,
                                                  double& outEntryLocation, double& outExitLocation,
                                                  vector3f_type* outEntryNormal, vector3f_type* outExitNormal ) const {
    return detail::intersect_with_box( *this, boxMin, boxMax, outEntryLocation, outExitLocation, outEntryNormal,
                                       outExitNormal );
}

template <class FloatType>
bool ray3t<FloatType>::is_intersecting_box_volume( const boundbox3t<FloatType>& box, float_type start,
                                                   float_type end ) const {
    return detail::is_intersecting_box_volume( *this, box, start, end );
}
} // namespace graphics
} // namespace frantic
