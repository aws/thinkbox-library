// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace geometry {

class trimesh3; // forward declaration

/**
 * This class holds the data for one intersection
 *
 */
struct scan_conversion_intersection {
    // The z coordinate of this intersection
    float z;
    // The face index and barycentric coordinates of this intersection
    int faceIndex;
    graphics::vector3f barycentricCoord;
    // Whether the normal of this face is pointing towards positive Z
    bool normalFacingZPositive;

    /**
     * First sort by ascending Z, and if the Z values are equal, consider
     * the ones with normal facing towards negative Z as smaller than the ones with normal
     * facing towards positive Z.
     *
     */
    bool operator<( const scan_conversion_intersection& rhs ) const {
        if( z == rhs.z )
            return !normalFacingZPositive && ( normalFacingZPositive != rhs.normalFacingZPositive );
        else
            return z < rhs.z;
    }
};

/**
 * Scan converts an input triangle mesh into a series of surface depths.  All the triangles are converted into samples
 * on the grid defined by the input parameters, and added to the appropriate elements of the outIntersectionDepths
 * array. <p> Once the scan conversion is complete, each subarray is sorted in the ascending Z direction.
 *
 * @param  xform      The transform matrix which should place the input mesh in the correct position for scan
 * conversion.
 * @param  mesh       The input mesh to be scan converted.
 * @param  dimensions The XY dimensions of the output intersection depths.  The scan conversion operates in the space
 *                    from [0,0,-infinity] to [dimensions.xsize,dimensions.ysize,+infinity].  The actual coordinates
 * sampled are at the middle of the squares, so for instance index 0 in the output array refers to the line
 *                    [0.5,0.5,-infinity] to [0,5,0,5,+infinity].
 * @param  outIntersectionDepths  This output parameter is where all the scan converted samples are placed.  When the
 * function is called, it should already be initialized to the right size.  Each element of this array consists of the
 *                                Z coordinate of the surface, as well as true if the normal of the surface is pointing
 * upwards.
 */
void trimesh3_scan_convert( const graphics::transform4f& xform, const trimesh3& mesh, graphics2d::size2 dimensions,
                            std::vector<std::vector<scan_conversion_intersection>>& outIntersectionDepths );

} // namespace geometry
} // namespace frantic
