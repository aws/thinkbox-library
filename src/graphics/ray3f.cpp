// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// vector3f.cpp
//
// Ray functions.

// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics/ray3f.hpp>

using namespace frantic::graphics;

// Intersects the ray with the box.  Uses the range [0, infinity) for t.
// ray.At(t0) is the first intersection found, with the given normal
// ray.At(t1) is the second intersection
template <class FloatType>
bool intersect_with_box_templ( const ray3t<FloatType>& ray, const vector3t<FloatType>& boxMin,
                               const vector3t<FloatType>& boxMax, double& outEntryLocation, double& outExitLocation,
                               vector3t<FloatType>* outEntryNormal, vector3t<FloatType>* outExitNormal ) {
    typedef vector3t<FloatType> vector3f_type;

    // TODO: Make another version of this that doesn't return the normal directions, to speed up those cases.
    double tNear = -( std::numeric_limits<double>::max )(), tFar = ( std::numeric_limits<double>::max )();
    bool swappedBounds = false;

    if( ray.direction().x == 0 ) {
        if( ( boxMin.x > ray.origin().x ) || ( ray.origin().x > boxMax.x ) ) {
            outEntryLocation = outExitLocation = 0;
            return false;
        }
    } else {
        double invDirX = 1.0 / ray.direction().x;
        double tMin = ( (double)boxMin.x - ray.origin().x ) * invDirX,
               tMax = ( (double)boxMax.x - ray.origin().x ) * invDirX;
        swappedBounds = false;
        if( tMin > tMax ) {
            std::swap( tMin, tMax );
            swappedBounds = true;
        }
        if( tMin > tNear ) { // largest tNear
            tNear = tMin;
            if( outEntryNormal )
                *outEntryNormal = swappedBounds ? vector3f_type::from_xaxis() : -vector3f_type::from_xaxis();
        }
        if( tMax < tFar ) { // closest tFar
            tFar = tMax;
            if( outExitNormal )
                *outExitNormal = ( !swappedBounds ) ? vector3f_type::from_xaxis() : -vector3f_type::from_xaxis();
        }
        if( tFar < 0 ) {
            outEntryLocation = outExitLocation = 0;
            return false;
        }
    }

    if( ray.direction().y == 0 ) {
        if( ( boxMin.y > ray.origin().y ) || ( ray.origin().y > boxMax.y ) ) {
            outEntryLocation = outExitLocation = 0;
            return false;
        }
    } else {
        double invDirY = 1.0 / ray.direction().y;
        double tMin = ( (double)boxMin.y - ray.origin().y ) * invDirY,
               tMax = ( (double)boxMax.y - ray.origin().y ) * invDirY;
        swappedBounds = false;
        if( tMin > tMax ) {
            std::swap( tMin, tMax );
            swappedBounds = true;
        }
        if( tMin > tNear ) { // largest tNear
            tNear = tMin;
            if( outEntryNormal )
                *outEntryNormal = swappedBounds ? vector3f_type::from_yaxis() : -vector3f_type::from_yaxis();
        }
        if( tMax < tFar ) { // closest tFar
            tFar = tMax;
            if( outExitNormal )
                *outExitNormal = ( !swappedBounds ) ? vector3f_type::from_yaxis() : -vector3f_type::from_yaxis();
        }
        if( tFar < 0 ) {
            outEntryLocation = outExitLocation = 0;
            return false;
        }
    }

    if( ray.direction().z == 0 ) {
        if( ( boxMin.z > ray.origin().z ) || ( ray.origin().z > boxMax.z ) ) {
            outEntryLocation = outExitLocation = 0;
            return false;
        }
    } else {
        double invDirZ = 1.0 / ray.direction().z;
        double tMin = ( (double)boxMin.z - ray.origin().z ) * invDirZ,
               tMax = ( (double)boxMax.z - ray.origin().z ) * invDirZ;
        swappedBounds = false;
        if( tMin > tMax ) {
            std::swap( tMin, tMax );
            swappedBounds = true;
        }
        if( tMin > tNear ) { // largest tNear
            tNear = tMin;
            if( outEntryNormal )
                *outEntryNormal = swappedBounds ? vector3f::from_zaxis() : -vector3f::from_zaxis();
        }
        if( tMax < tFar ) { // closest tFar
            tFar = tMax;
            if( outExitNormal )
                *outExitNormal = ( !swappedBounds ) ? vector3f::from_zaxis() : -vector3f::from_zaxis();
        }
        if( tFar < 0 ) {
            outEntryLocation = outExitLocation = 0;
            return false;
        }
    }

    outEntryLocation = tNear;
    outExitLocation = tFar;
    // The comparison between tNear and tFar is <= because we're considering the bound box as a closed set.  This case
    // happens when a ray is just grazing the edge.
    return ( tNear <= tFar ) && ( tFar >= 0 );
}

bool detail::intersect_with_box( const ray3f& ray, const vector3f& boxMin, const vector3f& boxMax,
                                 double& outEntryLocation, double& outExitLocation, vector3f* outEntryNormal,
                                 vector3f* outExitNormal ) {
    return intersect_with_box_templ<float>( ray, boxMin, boxMax, outEntryLocation, outExitLocation, outEntryNormal,
                                            outExitNormal );
}
bool detail::intersect_with_box( const ray3fd& ray, const vector3fd& boxMin, const vector3fd& boxMax,
                                 double& outEntryLocation, double& outExitLocation, vector3fd* outEntryNormal,
                                 vector3fd* outExitNormal ) {
    return intersect_with_box_templ<double>( ray, boxMin, boxMax, outEntryLocation, outExitLocation, outEntryNormal,
                                             outExitNormal );
}

template <class FloatType>
bool is_intersecting_box_volume_templ( const ray3t<FloatType>& ray, const boundbox3t<FloatType>& box, FloatType start,
                                       FloatType end ) {
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;

    const vector3f_type startPoint( ray.at( start ) );
    const vector3f_type endPoint( ray.at( end ) );
    if( box.is_empty() ) {
        return false;
    } else if( box.contains( startPoint ) || box.contains( endPoint ) ) {
        return true;
    } else if( ( startPoint.x < box.minimum().x && endPoint.x < box.minimum().x ) ||
               ( startPoint.y < box.minimum().y && endPoint.y < box.minimum().y ) ||
               ( startPoint.z < box.minimum().z && endPoint.z < box.minimum().z ) ||
               ( startPoint.x > box.maximum().x && endPoint.x > box.maximum().x ) ||
               ( startPoint.y > box.maximum().y && endPoint.y > box.maximum().y ) ||
               ( startPoint.z > box.maximum().z && endPoint.z > box.maximum().z ) ) {
        // A number of special cases of guaranteed non-intersection, which were very numerically
        // unstable in the test below.
        return false;
    } else {
        float_type minimumOfFarDistances = ( std::numeric_limits<float_type>::max )();
        float_type maximumOfNearDistances = ( std::numeric_limits<float_type>::min )();

        if( ray.direction().x != 0 ) {
            float_type invDirX = 1.f / ray.direction().x;
            if( invDirX >= 0 ) {
                maximumOfNearDistances = ( box.minimum().x - ray.origin().x ) * invDirX;
                minimumOfFarDistances = ( box.maximum().x - ray.origin().x ) * invDirX;
            } else {
                maximumOfNearDistances = ( box.maximum().x - ray.origin().x ) * invDirX;
                minimumOfFarDistances = ( box.minimum().x - ray.origin().x ) * invDirX;
            }
        }

        if( ray.direction().y != 0 ) {
            float_type invDirY = 1.f / ray.direction().y;
            if( invDirY >= 0 ) {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.minimum().y - ray.origin().y ) * invDirY );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.maximum().y - ray.origin().y ) * invDirY );
            } else {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.maximum().y - ray.origin().y ) * invDirY );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.minimum().y - ray.origin().y ) * invDirY );
            }
        }

        if( ray.direction().z != 0 ) {
            float_type invDirZ = 1.f / ray.direction().z;
            if( invDirZ >= 0 ) {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.minimum().z - ray.origin().z ) * invDirZ );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.maximum().z - ray.origin().z ) * invDirZ );
            } else {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.maximum().z - ray.origin().z ) * invDirZ );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.minimum().z - ray.origin().z ) * invDirZ );
            }
        }

        if( maximumOfNearDistances >= end || minimumOfFarDistances <= start ) {
            return false;
        }

        return ( maximumOfNearDistances <= minimumOfFarDistances );
    }
}

bool detail::is_intersecting_box_volume( const ray3f& ray, const boundbox3f& box, float start, float end ) {
    return is_intersecting_box_volume_templ<float>( ray, box, start, end );
}
bool detail::is_intersecting_box_volume( const ray3fd& ray, const boundbox3fd& box, double start, double end ) {
    return is_intersecting_box_volume_templ<double>( ray, box, start, end );
}
