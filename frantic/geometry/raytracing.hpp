// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/ray3f.hpp>

namespace frantic {
namespace geometry {

using frantic::graphics::ray3f;
using frantic::graphics::vector3f;

// Computes the greedy SAH cost of a given split.
inline float greedy_sah_cost( float boundsSurfaceArea, float leftBoundsSurfaceArea, float rightBoundsSurfaceArea,
                              int leftPrimitiveCount, int perpPrimitiveCount, int rightPrimitiveCount,
                              int& outBestPerpSide ) {
    // C(T) = K_T + P[V1|V]C(V_l) + P[V2|V]C(V_r)
    //      = K_T + SA(leftBounds)/SA(bounds) * leftPrimitiveCount + SA(rightBounds)/SA(bounds) * rightPrimitiveCount
    float probabilityLeft = leftBoundsSurfaceArea / boundsSurfaceArea,
          probabilityRight = rightBoundsSurfaceArea / boundsSurfaceArea;

    // Setting K_T to 0 for now...
    float costLeftPerp =
        ( leftPrimitiveCount + perpPrimitiveCount ) * probabilityLeft + rightPrimitiveCount * probabilityRight;
    float costRightPerp =
        leftPrimitiveCount * probabilityLeft + ( rightPrimitiveCount + perpPrimitiveCount ) * probabilityRight;

    // Empty space heuristic - empty nodes get bonus points
    //		if( leftPrimitiveCount + perpPrimitiveCount == 0 || rightPrimitiveCount == 0 )
    //			costLeftPerp *= 0.85f;
    //		if( rightPrimitiveCount == 0 || rightPrimitiveCount + perpPrimitiveCount == 0 )
    //			costLeftPerp *= 0.85f;

    // Choose the best of the right or left split
    if( costLeftPerp < costRightPerp ) {
        outBestPerpSide = -1;
        return costLeftPerp;
    } else {
        outBestPerpSide = +1;
        return costRightPerp;
    }
}

// This class stores the result of a trimesh3 ray intersection.
struct raytrace_intersection {
    ray3f ray;
    double distance;
    vector3f position;
    vector3f geometricNormal;
    int faceIndex;
    vector3f barycentricCoords;
    int primitiveIndex;

    /// Currently set only by the particle_collision_detector, not by the regular ray intersections
    vector3f motionDuringTimeStep;

    /**
     * Comparison based on distance to the intersection, so that we can sort intersections by distance.
     */
    bool operator<( const raytrace_intersection& rhs ) const { return distance < rhs.distance; }
};

// This class stores the result of a trimesh3 nearest point search.
struct nearest_point_search_result {
    float distance;
    vector3f position;
    vector3f geometricNormal;
    int faceIndex;
    vector3f barycentricCoords;
    int primitiveIndex; // added this because I need to access the primitive index in my searches

    bool operator<( const nearest_point_search_result& rhs ) const { return distance < rhs.distance; }
};

} // namespace geometry
} // namespace frantic
