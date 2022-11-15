// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/vector2f.hpp>

namespace frantic {
namespace graphics2d {

namespace detail {
/**
 * Sort the triangle points according to y coordinate
 * @param points an array (of size 3) containing the vertex coordinates in a triangle
 * @param perm an array (of size 3) containing the vertex indices of triangle. At the end of the function
 *		the indices will be sorted so that the perm[0] indicates the point
 *		with the lowest y value and perm[2] indicates the point with the
 *		highest y value
 */
inline void sort_triangle_points( frantic::graphics2d::vector2f points[3], size_t perm[3] ) {
    if( points[perm[0]].y > points[perm[1]].y ) {
        std::swap( perm[0], perm[1] );
    }
    if( points[perm[0]].y > points[perm[2]].y ) {
        std::swap( perm[0], perm[2] );
    }
    if( points[perm[1]].y > points[perm[2]].y ) {
        std::swap( perm[1], perm[2] );
    }
}

/**
 * Perform conservative rasterization on the left edge
 * @param pixelY scanline being considered
 * @param v0 triangle point with the most extreme y coordinate
 * @param v1 right most triangle point
 * @param v2 left most triangle point
 * @param leftEdge vector from v0 to v2
 * @param leftGradient gradient of leftEdge
 * @param scans overall starting, middle and ending scanline for overall triangle
 * @return the left most pixel on scanline that need to be colored
 */
inline int left_conservative( const int pixelY, const frantic::graphics2d::vector2f v0,
                              const frantic::graphics2d::vector2f v1, const frantic::graphics2d::vector2f v2,
                              const vector2f leftEdge, const float leftGradient, const int scans[3] ) {
    using namespace frantic::graphics2d;

    vector2f leftPixelVal;
    // leftPixelVal is the cartesian coordinates of
    // the left most lower right corner or
    // upper right corner (depending on the gradient of left
    // edge) such that leftPixelVal lies to the right of the
    // left edge
    leftPixelVal[1] = pixelY + ( leftGradient > 0 ? 0.f : 1.f );
    leftPixelVal[0] = floorf( ( ( leftPixelVal[1] - v2.y ) * ( -leftEdge.x ) ) / ( -leftEdge.y ) + v2.x + 1.0f );

    // Account for horizontal line
    if( leftEdge.y == 0 ) {
        leftPixelVal[0] = floorf( v2.x );
    }

    // At the starting and ending scan lines, just use the relevant point as
    // leftPixelVal[0]. This is important if 2 edges in the triangle
    // are almost horizontal
    if( pixelY == scans[0] ) {
        if( leftGradient > 0 ) {
            leftPixelVal[0] = v0.y < v2.y ? ceilf( v0.x ) : ceilf( v2.x );
        }
    }
    if( pixelY == scans[2] - 1 ) {
        if( leftGradient < 0 ) {
            leftPixelVal[0] = v0.y < v2.y ? ceilf( v2.x ) : ceilf( v0.x );
        }
    }

    // If the triangle is split such that the two halves
    // of the triangle have different left edges then grab the
    // middle point for the middle scanline if it is higher than
    // the current leftPixelVal[0]. This is important for almost horizontal
    // line
    if( pixelY == scans[1] - 1 ) {
        if( ( v0.y <= v2.y && v2.y <= v1.y ) || ( v0.y >= v2.y && v2.y >= v1.y ) ) {
            if( v0.x > v2.x ) {
                leftPixelVal[0] = std::max( ceilf( v2.x ), leftPixelVal[0] );
            } else {
                leftPixelVal[0] = std::max( ceilf( v0.x ), leftPixelVal[0] );
            }
        }
    }

    int leftX = static_cast<int>( leftPixelVal[0] ) - 1;
    return leftX;
}

/**
 * Perform interior rasterization for the left edge
 * @param pixelY scanline being considered
 * @param v2 left most triangle point
 * @param leftEdge the left edge of the triangle
 * @param conservativeVertex boolean to indicate if the v2 is part of an edge conservatively rasterized
 * @return the left most pixel on scanline that need to be colored
 */
inline int left_interior( const int pixelY, const vector2f v2, const vector2f leftEdge,
                          const bool conservativeVertex ) {

    using namespace frantic::graphics2d;

    vector2f leftPixelVal;

    // leftPixelVal is the pixel whose center lies on or
    // to the right of the leftEdge
    leftPixelVal[1] = static_cast<float>( pixelY ) + 0.5f;
    leftPixelVal[0] = ( ( leftPixelVal[1] - v2.y ) * ( -leftEdge.x ) ) / ( -leftEdge.y ) + v2.x;
    leftPixelVal[0] = ceilf( leftPixelVal[0] - 0.5f );

    // Take care of horizontal lines
    if( leftEdge.y == 0 ) {
        leftPixelVal[0] = conservativeVertex ? floorf( v2.x ) : ceilf( v2.x - 0.5f );
    }
    int leftX = static_cast<int>( leftPixelVal[0] );
    return leftX;
}

/**
 * Perform conservative rasterization on the right edge
 * @param pixelY scanline being considered
 * @param v0 triangle point with the most extreme y coordinate
 * @param v1 right most triangle point
 * @param v2 left most triangle point
 * @param rightEdge vector from v0 to v1
 * @param rightGradient gradient of rightEdge
 * @param scans overall starting, middle and ending scanline for overall triangle
 * @return the right most pixel on scanline that need to be colored + 1
 */
inline int right_conservative( const int pixelY, const frantic::graphics2d::vector2f v0,
                               const frantic::graphics2d::vector2f v1, const frantic::graphics2d::vector2f v2,
                               const vector2f rightEdge, const float rightGradient, const int scans[3] ) {

    using namespace frantic::graphics2d;

    vector2f rightPixelVal;

    // rightPixelVal is the cartesian coordinates of
    // the right most lower left corner or
    // lower right corner (depending on the gradient of right
    // edge) such that rightPixelVal lies to the left of the
    // right edge
    rightPixelVal[1] = pixelY + ( rightGradient < 0 ? 0.f : 1.f );
    rightPixelVal[0] = ceilf( ( ( rightPixelVal[1] - v1.y ) * ( -rightEdge.x ) ) / ( -rightEdge.y ) + v1.x - 1.0f );

    // Account for horizontal lines
    if( rightEdge.y == 0 ) {
        rightPixelVal[0] = floorf( v1.x );
    }

    // At the starting and ending scan lines, just use the relevant point as
    // rightPixelVal[0]. This is important if 2 edges in the triangle
    // are almost horizontal
    if( pixelY == scans[0] ) {
        if( rightGradient < 0 ) {
            rightPixelVal[0] = v0.y < v1.y ? floorf( v0.x ) : floorf( v1.x );
        }
    }
    if( pixelY == scans[2] - 1 ) {
        if( rightGradient > 0 ) {
            rightPixelVal[0] = v0.y < v1.y ? floorf( v1.x ) : floorf( v0.x );
        }
    }

    // If the triangle is split such that the two halves
    // of the triangle have different right edges then grab the
    // middle point for the middle scanline if it is higher than
    // the current rightPixelVal[0]. This is important for almost horizontal
    // line
    if( pixelY == scans[1] - 1 ) {
        if( ( v0.y <= v2.y && v2.y >= v1.y ) || ( v0.y >= v2.y && v2.y <= v1.y ) ) {
            if( v0.x < v1.x ) {
                rightPixelVal[0] = std::min( rightPixelVal[0], floorf( v1.x ) );
            } else {
                rightPixelVal[0] = std::min( rightPixelVal[0], floorf( v0.x ) );
            }
        }
    }
    int rightX = static_cast<int>( rightPixelVal[0] ) + 1;
    return rightX;
}

/**
 * Perform interior rasterization for the right edge
 * @param pixelY scanline being considered
 * @param v1 right most triangle point
 * @param rightEdge the right edge of the triangle
 * @param conservativeVertex boolean to indicate if the v1 is part of an edge conservatively rasterized
 * @return the right most pixel on scanline that need to be colored + 1
 */
inline int right_interior( const int pixelY, const vector2f v1, const vector2f rightEdge,
                           const bool conservativeVertex ) {

    using namespace frantic::graphics2d;

    vector2f rightPixelVal;

    // rightPixelVal is the pixel whose center lies on or
    // to the left of the rightEdge
    rightPixelVal[1] = static_cast<float>( pixelY ) + 0.5f;
    rightPixelVal[0] =
        frantic::math::round( ( ( rightPixelVal[1] - v1.y ) * ( -rightEdge.x ) ) / ( -rightEdge.y ) + v1.x );
    if( rightEdge.y == 0 ) {
        rightPixelVal[0] = conservativeVertex ? ceilf( v1.x ) : floorf( v1.x + 0.5f );
    }
    int rightX = static_cast<int>( rightPixelVal[0] );
    return rightX;
}

/**
 * Collects information about rasterization of a triangle into a buffer
 * (specifically for scanlines starting from startScan to endScan (not inclusive))
 * This is done with the help of edge functions as shown in
 * https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/rasterization-stage
 * @param verts coordinates of triangle vertices
 * @param perm indices of triangle vertices
 *		perm[0] the highest point in the y direction
 *		perm[1] the right side leg of the triangle half
 *		perm[2] the left side leg of the triangle half
 * @param startScan starting scanline for triangle half
 * @param endScan ending scanline (not inclusive) for triangle half
 * @param scans starting, middle and ending scanline for overall triangle
 * @param imageWidth the width of the buffer
 * @param conservativeEdges an array of booleans indicating if an edge needs
  to be conservatively rasterized or not
    conservativeEdges[0] indicates edge between vertex 0 and 1
    conservativeEdges[1] indicates edge between vertex 1 and 2
    conservativeEdges[2] indicates edge between vertex 2 and 0
 * @param outRuns vector of structs holding the relevant rasterization information
 */
template <typename RunType>
inline void rasterize_triangle_half( const frantic::graphics2d::vector2f verts[3], size_t perm[3], const int startScan,
                                     const int endScan, const int scans[3], const int imageWidth,
                                     const bool conservativeEdges[3], std::vector<RunType>& outRuns ) {

    using namespace frantic::graphics2d;
    vector2f v0 = verts[perm[0]], v1 = verts[perm[1]], v2 = verts[perm[2]];

    // Figure out the left edge
    vector2f leftEdge = v2 - v0;
    float leftGradient = leftEdge.y / leftEdge.x;
    size_t minLeftVertInd = std::min( perm[0], perm[2] );
    size_t maxLeftVertInd = std::max( perm[0], perm[2] );

    // Figure out the right edge
    vector2f rightEdge = v1 - v0;
    float rightGradient = rightEdge.y / rightEdge.x;
    size_t minRightVertInd = std::min( perm[0], perm[1] );
    size_t maxRightVertInd = std::max( perm[0], perm[1] );

    // Find out if the left edge or the right edges are conservative edges
    bool leftConservative = false, rightConservative = false;
    if( ( conservativeEdges[minLeftVertInd] && maxLeftVertInd == ( minLeftVertInd + 1 ) % 3 ) ||
        ( conservativeEdges[maxLeftVertInd] && minLeftVertInd == ( maxLeftVertInd + 1 ) % 3 ) ) {
        leftConservative = true;
    }
    if( ( conservativeEdges[minRightVertInd] && maxRightVertInd == ( minRightVertInd + 1 ) % 3 ) ||
        ( conservativeEdges[maxRightVertInd] && minRightVertInd == ( maxRightVertInd + 1 ) % 3 ) ) {
        rightConservative = true;
    }

    // Store information about vertices that are part of
    // edges that undergoes conservative rasterization
    bool conservativeVertices[3] = { false, false, false };
    for( int b = 0; b < 3; ++b ) {
        if( conservativeEdges[b] ) {
            conservativeVertices[b] = true;
            conservativeVertices[( b + 1 ) % 3] = true;
        }
    }

    for( int pixelY = startScan; pixelY < endScan; ++pixelY ) {
        int leftX, rightX;

        vector2f leftPixelVal, rightPixelVal;

        if( leftConservative ) {
            leftX = left_conservative( pixelY, v0, v1, v2, leftEdge, leftGradient, scans );

        } else {
            leftX = left_interior( pixelY, v2, leftEdge, conservativeVertices[perm[2]] );
        }
        if( rightConservative ) {
            rightX = right_conservative( pixelY, v0, v1, v2, rightEdge, rightGradient, scans );

        } else {
            rightX = right_interior( pixelY, v1, rightEdge, conservativeVertices[perm[1]] );
        }

        if( leftX < 0 ) {
            leftX = 0;
        }
        if( rightX > imageWidth ) {
            rightX = imageWidth;
        }

        vector2f offset( ( leftX + 0.5f - v0.x ), ( static_cast<float>( pixelY ) + 0.5f - v0.y ) );
        outRuns.push_back( RunType( leftX, rightX, offset, leftEdge, rightEdge, perm ) );
    }
}
} // namespace detail

/**
 * Data structure to hold rasterization information - specifically,
 * the range (from left to right) of columns in a specific row in a buffer touched by a triangle
 */
struct simple_run {

    // Holds the range of column indices in a specific row touched
    // by a triangle in the buffer
    std::pair<int, int> columnRange;

    /**
     * @param leftX left most pixel in the specific row touched by triangle in the buffer
     * @param rightX right most pixel in the specific row touched by triangle in the buffer + 1
     * @param  pixelOffset(not needed for this run)
     * @param leftEdge (not needed for this run)
     * @param rightEdge (not needed for this run)
     * @param perm (not needed for this run)
     */
    simple_run( const int leftX, const int rightX, const vector2f /*pixelOffset*/, const vector2f /*leftEdge*/,
                const vector2f /*rightEdge*/, const size_t* /*perm*/ )
        : columnRange( leftX, rightX ) {}
};

/**
 * Data structure to hold rasterization information - specifically,
 * the range (from left to right) of columns in a specific row in a buffer touched by a triangle
 * and the relevant barycentric coordinate information. The barycentric coordinate information
 * will contain the barycentric coordinate at the left column and the the difference between
 * barycentric coordinates per pixel in the specific row.
 */
struct barycentric_run {
    /* Holds the range of column indices */
    std::pair<int, int> columnRange;

    /* The starting barycentric coordinate. */
    frantic::graphics::vector3f firstBarycentricCoord;

    /*The difference between the barycentric coordinates per pixel in the specific row */
    frantic::graphics::vector3f deltaBarycentricCoord;

    /**
     * @param leftX left most pixel in the specific row touched by triangle in the buffer
     * @param rightX right most pixel in the specific row touched by triangle in the buffer + 1
     * @param pixelOffset edge from the most extreme point to the middle of the pixel being considered
     * @param leftEdge the left edge of the triangle
     * @param rightEdge the right edge of the triangle
     * @param perm indices of triangle vertices
     *		the structure for permutations is
     *		perm[0] the highest point in the y direction
     *		perm[1] the right side leg of the triangle half
     *		perm[2] the left side leg of the triangle half
     */
    barycentric_run( const int leftX, const int rightX, const vector2f pixelOffset, const vector2f leftEdge,
                     const vector2f rightEdge, const size_t perm[3] ) {
        using namespace frantic::graphics2d;
        using namespace frantic::graphics;

        columnRange = std::make_pair( leftX, rightX );

        float leftRightCross = rightEdge.x * leftEdge.y - leftEdge.x * rightEdge.y;
        vector3f leftBaryCoord;
        leftBaryCoord[perm[1]] = ( pixelOffset[0] * leftEdge.y - leftEdge.x * pixelOffset[1] ) / leftRightCross;
        leftBaryCoord[perm[2]] = ( rightEdge.x * pixelOffset[1] - pixelOffset[0] * rightEdge.y ) / leftRightCross;
        leftBaryCoord[perm[0]] = 1 - leftBaryCoord[perm[1]] - leftBaryCoord[perm[2]];

        vector3f baryCoordDiff;
        baryCoordDiff[perm[1]] = leftEdge.y / leftRightCross;
        baryCoordDiff[perm[2]] = -rightEdge.y / leftRightCross;
        baryCoordDiff[perm[0]] = -baryCoordDiff[perm[1]] - baryCoordDiff[perm[2]];

        firstBarycentricCoord = leftBaryCoord;
        deltaBarycentricCoord = baryCoordDiff;
    }
};

/**
 * Rasterize a triangle on a buffer and collect the resulting rasterization information
 *
 * The edges of the triangle can be rasterized using either conservative rasterization
 * or interior rasterization.
 * In conservative rasterization, if an edge intersects a pixel then that pixel is considered to
 * be the most extreme pixel for a specific scanline.
 * In interior rasterization, the most exteme pixel whose center lies to the right (for the left edge)
 * and to the left (for the right edge) is considered to be the most extreme pixel for a
 * specific scanline.
 * For example, if the left edge of the triangle is conservatively rasterized then any pixel
 * touched by the left edge will be the left most pixel for every scan line. But if the right
 * edge of the triangle undergoes interior rasterization, then only the rightmost pixel whose center
 * is to the left of the edge will be considered to be the right most pixel for every scan line.
 *
 * @param p0 a point on the triangle
 * @param p1 a point on the triangle
 * @param p2 a point on the triangle
 * @param imageWidth width of the buffer
 * @param imageHeight height of the buffer
 * @param conservativeEdges
 *		an array of booleans indicating if an edge of the triangle
 *		needs to be conservatively rasterized or not
 *		conservativeEdges[0] indicates edge between vertex 0 and 1
 *		conservativeEdges[1] indicates edge between vertex 1 and 2
 *		conservativeEdges[2] indicates edge between vertex 2 and 0
 * @param outFirstRow the lowest row in the buffer touched by the triangle
 * @param outRuns vector of structs holding the relevant rasterization information - the runs indicating pixels filled
 *		and possibly the barycentric coordinates
 */
template <typename RunType>
void rasterize_triangle( const frantic::graphics2d::vector2f& p0, const frantic::graphics2d::vector2f& p1,
                         const frantic::graphics2d::vector2f& p2, const int imageWidth, const int imageHeight,
                         const bool conservativeEdges[3], int& outFirstRow, std::vector<RunType>& outRuns ) {

    using namespace frantic::graphics2d;
    vector2f points[3] = { p0, p1, p2 };
    size_t perm[3] = { 0, 1, 2 };
    bool conservativeVerts[3] = { false, false, false };

    for( int b = 0; b < 3; ++b ) {
        if( conservativeEdges[b] ) {
            conservativeVerts[b] = true;
            conservativeVerts[( b + 1 ) % 3] = true;
        }
    }

    detail::sort_triangle_points( points, perm );

    // Round the y values of the points appropriately
    int roundedY[3];

    roundedY[perm[0]] = static_cast<int>( conservativeVerts[perm[0]] ? floorf( points[perm[0]].y )
                                                                     : ceilf( points[perm[0]].y - 0.5f ) );

    roundedY[perm[1]] = static_cast<int>( conservativeVerts[perm[1]] ? ceilf( points[perm[1]].y )
                                                                     : frantic::math::round( points[perm[1]].y ) );

    roundedY[perm[2]] = static_cast<int>( conservativeVerts[perm[2]] ? ceilf( points[perm[2]].y )
                                                                     : frantic::math::round( points[perm[2]].y ) );

    // get the breaking points of the triangle's upper/lower pieces
    int yStart = frantic::math::clamp( roundedY[perm[0]], 0, imageHeight );
    int yMiddle = frantic::math::clamp( roundedY[perm[1]], 0, imageHeight );
    int yEnd = frantic::math::clamp( roundedY[perm[2]], 0, imageHeight );

    const vector2f e01 = points[perm[1]] - points[perm[0]];
    const vector2f e02 = points[perm[2]] - points[perm[0]];

    float cross = vector2f::cross( e01, e02 );
    // the structure for permutations is
    // perm[0] the extreme point in the y direction
    // perm[1] the right side leg of the triangle half
    // perm[2] the left side leg of the triangle half

    // degenerate, don't render
    if( cross == 0.0f ) {
        return;
    } else if( cross < 0.0f ) {
        // triangle is CW, swap left and right ends
        std::swap( perm[1], perm[2] );
    }

    outFirstRow = yStart;
    int scans[3] = { yStart, yMiddle, yEnd };
    if( yStart < yMiddle ) {
        detail::rasterize_triangle_half<RunType>( points, perm, yStart, yMiddle, scans, imageWidth, conservativeEdges,
                                                  outRuns );
    }
    if( cross < 0.0f ) {
        std::swap( perm[0], perm[1] );
    } else {
        std::swap( perm[0], perm[2] );
    }

    if( yMiddle < yEnd ) {
        detail::rasterize_triangle_half<RunType>( points, perm, yMiddle, yEnd, scans, imageWidth, conservativeEdges,
                                                  outRuns );
    }
}

} // namespace graphics2d
} // namespace frantic
