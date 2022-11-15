// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/ray3f.hpp>
#include <frantic/graphics/vector3f.hpp>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace frantic {
namespace graphics {

inline float triangle_area( const vector3f& a, const vector3f& b, const vector3f& c ) {
    float ad = vector3f::distance( a, b ), bd = vector3f::distance( b, c ), cd = vector3f::distance( c, a );
    float s = ( ad + bd + cd ) / 2;
    float areaSquared = s * ( s - ad ) * ( s - bd ) * ( s - cd );
    if( areaSquared > 0 )
        return sqrtf( areaSquared );
    return 0;
}

inline vector3f triangle_normal( const vector3f& a, const vector3f& b, const vector3f& c ) {
    vector3f result = vector3f::cross( c - b, a - b );
    float length = result.get_magnitude();
    if( length != 0 )
        result /= length;
    return result;
}

inline vector3f get_triangle_angles( const vector3f& a, const vector3f& b, const vector3f& c ) {
    vector3f b_a = b - a, c_b = c - b, a_c = a - c;

    // Normalize the vectors
    b_a.normalize();
    c_b.normalize();
    a_c.normalize();

    // Get the cosine of the angle
    float cosAlphaA = -vector3f::dot( b_a, a_c ), cosAlphaB = -vector3f::dot( b_a, c_b ),
          cosAlphaC = -vector3f::dot( c_b, a_c );

    // Just in case
    cosAlphaA = frantic::math::clamp( cosAlphaA, -1.f, 1.f );
    cosAlphaB = frantic::math::clamp( cosAlphaB, -1.f, 1.f );
    cosAlphaC = frantic::math::clamp( cosAlphaC, -1.f, 1.f );

    // Return the arc cosine
    return vector3f( acos( cosAlphaA ), acos( cosAlphaB ), acos( cosAlphaC ) );
}

/**
 * Compute the angle between two edges in a triangle.
 *
 * @return The angle of the corner formed by the edges from current to
 *         previous and from current to next.  Measured in radians.
 */
inline float get_triangle_angle( const frantic::graphics::vector3f& previous,
                                 const frantic::graphics::vector3f& current, const frantic::graphics::vector3f& next ) {
    frantic::graphics::vector3f next_current = next - current;
    frantic::graphics::vector3f previous_current = previous - current;

    next_current.normalize();
    previous_current.normalize();

    float cosAlpha = frantic::graphics::vector3f::dot( previous_current, next_current );

    cosAlpha = frantic::math::clamp( cosAlpha, -1.f, 1.f );

    return acos( cosAlpha );
}

inline bool intersect_line_segment_with_triangle( const vector3f& lineStart, const vector3f& lineEnd,
                                                  const vector3f& vert0, const vector3f& vert1, const vector3f& vert2,
                                                  vector3f& outIntersection ) {
    ray3f ray( lineStart, lineEnd - lineStart );
    vector3f baryCoord;
    double t;
    if( ray.intersect_with_triangle( vert0, vert1, vert2, t, baryCoord ) ) {
        if( t >= 0 && t <= 1 ) {
            outIntersection = ray.at( t );
            return true;
        }
    }
    return false;
}

template <typename FloatType>
inline vector3t<FloatType> linear_interpolate( const vector3t<FloatType>& a, const vector3t<FloatType>& b,
                                               FloatType t ) {
    FloatType one_t = 1 - t;
    return vector3t<FloatType>( one_t * a.x + t * b.x, one_t * a.y + t * b.y, one_t * a.z + t * b.z );
}

// Given the triangle a,b,c, this computes the two major axes along which
// to compute barycentric coordinates, and the inverse determinant of the 2x2
// matrix used to do that computation.
inline void compute_barycentric_helpers( const vector3f& a, const vector3f& b, const vector3f& c,
                                         char& outBarycentric0Axis, char& outBarycentric1Axis,
                                         float& outBarycentricInverseDeterminant ) {
    vector3f normal = triangle_normal( a, b, c );
    vector3f edge0 = b - a, edge1 = c - a;

    // Discard the axis along which the normal has greatest absolute value.
    // Store the remaining two axes in i0 and i1.
    int i0 = 1, i1 = 2;
    float maxAbsComponent = fabs( normal.x );
    if( fabs( normal.y ) > maxAbsComponent ) {
        maxAbsComponent = fabs( normal.y );
        i0 = 0;
        i1 = 2;
    }
    if( fabs( normal.z ) > maxAbsComponent ) {
        i0 = 0;
        i1 = 1;
    }

    // The resulting 2x2 linear equation is what we need to solve, so to help
    // that be done later, we compute the inverse determinant of that matrix
    outBarycentricInverseDeterminant = 1.f / ( edge0[i0] * edge1[i1] - edge1[i0] * edge0[i1] );
    outBarycentric0Axis = (char)i0;
    outBarycentric1Axis = (char)i1;
}

inline vector3f compute_barycentric_coordinates_with_helpers( const vector3f& position, const vector3f& a,
                                                              const vector3f& b, const vector3f& c,
                                                              int barycentric0Axis, int barycentric1Axis,
                                                              float barycentricInverseDeterminant ) {
    vector3f edge0 = b - a, edge1 = c - a;
    vector3f relativePosition = position - a;

    // Solve the resulting 2x2 linear equation
    float barycentricB =
        barycentricInverseDeterminant * ( edge1[barycentric1Axis] * relativePosition[barycentric0Axis] -
                                          edge1[barycentric0Axis] * relativePosition[barycentric1Axis] );
    float barycentricC =
        barycentricInverseDeterminant * ( edge0[barycentric0Axis] * relativePosition[barycentric1Axis] -
                                          edge0[barycentric1Axis] * relativePosition[barycentric0Axis] );

    return vector3f( 1.f - barycentricB - barycentricC, barycentricB, barycentricC );
}

// Slower than using the above cached version
inline vector3f compute_barycentric_coordinates( const vector3f& position, const vector3f& a, const vector3f& b,
                                                 const vector3f& c ) {
    char barycentric0Axis = 0, barycentric1Axis = 0;
    float barycentricInverseDeterminant = 0;
    compute_barycentric_helpers( a, b, c, barycentric0Axis, barycentric1Axis, barycentricInverseDeterminant );
    return compute_barycentric_coordinates_with_helpers( position, a, b, c, barycentric0Axis, barycentric1Axis,
                                                         barycentricInverseDeterminant );
}

// Returns vector from point to line perpendicular to the line
inline vector3f get_distance_vector_to_line( const vector3f& v, const vector3f& lineStart, const vector3f& lineEnd ) {
    vector3f lineVector = lineEnd - lineStart;
    float t = vector3f::dot( lineVector, ( v - lineStart ) ) / lineVector.get_magnitude_squared();
    return ( v - ( lineStart + t * lineVector ) );
}

// Computes the square of the distance from this point to a given line segment.
inline float get_distance_squared_to_line_segment( const vector3f& v, const vector3f& lineStart,
                                                   const vector3f& lineEnd ) {
    vector3f lineVector = lineEnd - lineStart;
    float t = vector3f::dot( lineVector, ( v - lineStart ) ) / lineVector.get_magnitude_squared();
    if( t <= 0 ) {
        return vector3f::distance_squared( v, lineStart );
    } else if( t >= 1 ) {
        return vector3f::distance_squared( v, lineEnd );
    } else {
        return vector3f::distance_squared( v, lineStart + t * lineVector );
    }
}

inline float get_distance_to_line_segment( const vector3f& v, const vector3f& lineStart, const vector3f& lineEnd ) {
    return sqrtf( get_distance_squared_to_line_segment( v, lineStart, lineEnd ) );
}

inline bool contained_in_triangle( const vector3f& v, const vector3f& a, const vector3f& b, const vector3f& c,
                                   float relativeTolerance = 0.0001f ) {
    // First make sure that the point is on the same plane as the triangle
    vector3f triangleNormal = triangle_normal( a, b, c );
    float distanceFromPlane = vector3f::dot( triangleNormal, v ) - vector3f::dot( triangleNormal, b );
    float absoluteMagnitude = a.max_abs_component();
    absoluteMagnitude = ( std::max )( absoluteMagnitude, b.max_abs_component() );
    absoluteMagnitude = ( std::max )( absoluteMagnitude, c.max_abs_component() );
    absoluteMagnitude = ( std::max )( absoluteMagnitude, v.max_abs_component() );
    if( distanceFromPlane > relativeTolerance * absoluteMagnitude )
        return false;
    // Then ensure that the point is actually inside the triangle
    vector3f baryCoords = compute_barycentric_coordinates( v, a, b, c );
    return baryCoords.x >= 0 && baryCoords.y >= 0 && baryCoords.z >= 0;
}

inline float distance_squared( const vector3f& a, const vector3f& b ) { return ( a - b ).get_magnitude_squared(); }

/**
 * This function takes vector3f sample inputs, and outputs a trilinearly interpolated value.
 * <p>
 * This function takes as input a size 3 array and a size 8 array.  The delta weights array consists of the distance
 * from the target interpolation point from the minimum of the sample cube, and the samples array are the corresponding
 * vector field samples.
 *
 * @param  trilerpDeltas  This array must be of length 3, and should consist of delta weighting factors for a trilinear
 * interpolation.
 * @param  vectorFieldSamples  This array must be of length 8, and consists of the values from the vector field.
 */
inline vector3f trilerp_value( const float* trilerpDeltas, const vector3f* vectorFieldSamples ) {
    float dx = trilerpDeltas[0], dy = trilerpDeltas[1], dz = trilerpDeltas[2];
    float dxAlt = 1 - dx, dyAlt = 1 - dy, dzAlt = 1 - dz;
    return dzAlt * ( dyAlt * ( dxAlt * vectorFieldSamples[0] + dx * vectorFieldSamples[1] ) +
                     dy * ( dxAlt * vectorFieldSamples[2] + dx * vectorFieldSamples[3] ) ) +
           dz * ( dyAlt * ( dxAlt * vectorFieldSamples[4] + dx * vectorFieldSamples[5] ) +
                  dy * ( dxAlt * vectorFieldSamples[6] + dx * vectorFieldSamples[7] ) );
}

/**
 * This function takes flot sample inputs, and outputs the gradient of the continous function created by trilinearly
 * interpolating the samples. <p> This function takes as input a size 3 array and a size 8 array.  The delta weights
 * array consists of the distance from the target interpolation point from the minimum of the sample cube, and the
 * samples array are the corresponding scalar field samples. The voxel length is required to take first derivatives.
 *
 * @param  trilerpDeltas  This array must be of length 3, and should consist of delta weighting factors for a trilinear
 * interpolation.
 * @param  scalarFieldSamples  This array must be of length 8, and consists of the values from the scalar field.
 * @param  voxelLength  This is the side length of the cube representing the 8 samples.
 */
inline vector3f trilerp_gradient( const float* trilerpDeltas, const float* scalarFieldSamples, float voxelLength ) {
    float dx = trilerpDeltas[0], dy = trilerpDeltas[1], dz = trilerpDeltas[2];
    float dxAlt = 1 - dx, dyAlt = 1 - dy, dzAlt = 1 - dz;

    vector3f result;

    // d F(x,y,z) / dx
    result.x = dzAlt * ( dyAlt * ( scalarFieldSamples[1] - scalarFieldSamples[0] ) +
                         dy * ( scalarFieldSamples[3] - scalarFieldSamples[2] ) ) +
               dz * ( dyAlt * ( scalarFieldSamples[5] - scalarFieldSamples[4] ) +
                      dy * ( scalarFieldSamples[7] - scalarFieldSamples[6] ) );

    // d F(x,y,z) / dy
    result.y = dzAlt * ( dxAlt * ( scalarFieldSamples[2] - scalarFieldSamples[0] ) +
                         dy * ( scalarFieldSamples[3] - scalarFieldSamples[1] ) ) +
               dz * ( dxAlt * ( scalarFieldSamples[6] - scalarFieldSamples[4] ) +
                      dy * ( scalarFieldSamples[7] - scalarFieldSamples[5] ) );

    // d F(x,y,z) / dz
    result.z = dyAlt * ( dxAlt * ( scalarFieldSamples[4] - scalarFieldSamples[0] ) +
                         dy * ( scalarFieldSamples[5] - scalarFieldSamples[1] ) ) +
               dy * ( dxAlt * ( scalarFieldSamples[6] - scalarFieldSamples[2] ) +
                      dy * ( scalarFieldSamples[7] - scalarFieldSamples[3] ) );

    result /= voxelLength;
    return result;
}

/**
 * This function takes vector3f sample inputs, and outputs a trilinearly interpolated curl value.

 * This function takes as input a size 3 array and a size 8 array.  The delta weights array consists of the distance
 * from the target interpolation point from the minimum of the sample cube, and the samples array are the corresponding
 * vector field samples. The voxel length is required to take first derivatives.
 *
 * @param  trilerpDeltas  This array must be of length 3, and should consist of delta weighting factors for a trilinear
 * interpolation.
 * @param  vectorFieldSamples  This array must be of length 8, and consists of the values from the vector field.
 * @param  voxelLength  This is the side length of the cube representing the 8 samples.
 */
inline vector3f trilerp_curl( const float* trilerpDeltas, const vector3f* vectorFieldSamples, float voxelLength ) {
    float dx = trilerpDeltas[0], dy = trilerpDeltas[1], dz = trilerpDeltas[2];
    float dxAlt = 1 - dx, dyAlt = 1 - dy, dzAlt = 1 - dz;
    float inverseVoxelLength = 1 / voxelLength;

    vector3f result;
    float partialA, partialB;

    //////////
    // Do X, (dv_z/dy - dv_y/dz)
    //////////
    // voxelLength * dv_z/dy
    partialA = dzAlt * ( ( dxAlt * vectorFieldSamples[2].z + dx * vectorFieldSamples[3].z ) -
                         ( dxAlt * vectorFieldSamples[0].z + dx * vectorFieldSamples[1].z ) ) +
               dz * ( ( dxAlt * vectorFieldSamples[6].z + dx * vectorFieldSamples[7].z ) -
                      ( dxAlt * vectorFieldSamples[4].z + dx * vectorFieldSamples[5].z ) );
    // voxelLength * dv_y/dz
    partialB = ( dyAlt * ( dxAlt * vectorFieldSamples[4].y + dx * vectorFieldSamples[5].y ) +
                 dy * ( dxAlt * vectorFieldSamples[6].y + dx * vectorFieldSamples[7].y ) ) -
               ( dyAlt * ( dxAlt * vectorFieldSamples[0].y + dx * vectorFieldSamples[1].y ) +
                 dy * ( dxAlt * vectorFieldSamples[2].y + dx * vectorFieldSamples[3].y ) );
    result.x = inverseVoxelLength * ( partialA - partialB );

    //////////
    // Do Y, (dv_x/dz - dv_z/dx)
    //////////
    // voxelLength * dv_x/dz
    partialA = ( dyAlt * ( dxAlt * vectorFieldSamples[4].x + dx * vectorFieldSamples[5].x ) +
                 dy * ( dxAlt * vectorFieldSamples[6].x + dx * vectorFieldSamples[7].x ) ) -
               ( dyAlt * ( dxAlt * vectorFieldSamples[0].x + dx * vectorFieldSamples[1].x ) +
                 dy * ( dxAlt * vectorFieldSamples[2].x + dx * vectorFieldSamples[3].x ) );
    // voxelLength * dv_z/dx
    partialB = dzAlt * ( dyAlt * ( vectorFieldSamples[1].z - vectorFieldSamples[0].z ) +
                         dy * ( vectorFieldSamples[3].z - vectorFieldSamples[2].z ) ) +
               dz * ( dyAlt * ( vectorFieldSamples[5].z - vectorFieldSamples[4].z ) +
                      dy * ( vectorFieldSamples[7].z - vectorFieldSamples[6].z ) );
    result.y = inverseVoxelLength * ( partialA - partialB );

    //////////
    // Do Z, (dv_y/dx - dv_x/dy)
    //////////
    // voxelLength * dv_y/dx
    partialA = dzAlt * ( dyAlt * ( vectorFieldSamples[1].y - vectorFieldSamples[0].y ) +
                         dy * ( vectorFieldSamples[3].y - vectorFieldSamples[2].y ) ) +
               dz * ( dyAlt * ( vectorFieldSamples[5].y - vectorFieldSamples[4].y ) +
                      dy * ( vectorFieldSamples[7].y - vectorFieldSamples[6].y ) );
    // voxelLength * dv_x/dy
    partialB = dzAlt * ( ( dxAlt * vectorFieldSamples[2].x + dx * vectorFieldSamples[3].x ) -
                         ( dxAlt * vectorFieldSamples[0].x + dx * vectorFieldSamples[1].x ) ) +
               dz * ( ( dxAlt * vectorFieldSamples[6].x + dx * vectorFieldSamples[7].x ) -
                      ( dxAlt * vectorFieldSamples[4].x + dx * vectorFieldSamples[5].x ) );
    result.z = inverseVoxelLength * ( partialA - partialB );

    return result;
}

template <class InputIterator>
typename std::iterator_traits<InputIterator>::value_type get_centroid( InputIterator begin, InputIterator end ) {
    typedef typename std::iterator_traits<InputIterator>::value_type vector_type;
    typedef typename vector_type::float_type float_type;

    size_t count = 0;
    vector_type acculm;

    for( InputIterator it = begin; it != end; ++it ) {
        acculm += *it;
        ++count;
    }

    return acculm / float_type( count );
}

} // namespace graphics
} // namespace frantic
