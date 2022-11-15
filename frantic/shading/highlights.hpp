// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once


#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/math/utils.hpp>

namespace frantic {
namespace shading {

// This function returns a unit vector tangent to "normal". It does extra work if
// the normal is more vertical than horizontal in order to avoid the degenerate case of
// a normal of (0,0,1).
// NOTE: This assumes normal is a unit vector.
inline frantic::graphics::vector3f compute_tangent( const frantic::graphics::vector3f& normal ) {
    using frantic::graphics::vector3f;

    if( std::abs( normal.z ) <= std::abs( normal.x ) ) {
        return vector3f::normalize( vector3f::cross( normal, vector3f( 0, 0, 1 ) ) );
    } else {
        return vector3f::cross( vector3f::normalize( vector3f::cross( normal, vector3f( 1, 0, 0 ) ) ), normal );
    }
}

// This function returns unit vectors tangent and binormal to "normal". It does extra work if
// the normal is more vertical than horizontal in order to avoid the degenerate case of
// a normal of (0,0,1).
// NOTE: This assumes normal is a unit vector.
inline void compute_tangent_binormal( const frantic::graphics::vector3f& normal, frantic::graphics::vector3f& tangent,
                                      frantic::graphics::vector3f& binormal ) {
    using frantic::graphics::vector3f;

    if( std::abs( normal.z ) <= std::abs( normal.x ) ) {
        tangent = vector3f::normalize( vector3f::cross( normal, vector3f( 0, 0, 1 ) ) );
        binormal = vector3f::cross( normal, tangent );
    } else {
        binormal = vector3f::normalize( vector3f::cross( normal, vector3f( 1, 0, 0 ) ) );
        tangent = vector3f::cross( binormal, normal );
    }
}

// This function returns a vector with magnitude equal to incidentDir, with direction reflected off the a
// surface "normal".
// NOTE: This assumes normal is a unit vector.
//  @param incidentDir Points towards the surface from the viewing location.
//  @param normal Points away from the surface.
inline frantic::graphics::vector3f compute_reflection( const frantic::graphics::vector3f& incidentDir,
                                                       const frantic::graphics::vector3f& normal ) {
    using frantic::graphics::vector3f;
    return incidentDir - 2 * vector3f::dot( incidentDir, normal ) * normal;
}

// This an isotropic specular highlight function created by Phong and Blinn.
// See: http://en.wikipedia.org/wiki/Specular_highlight#Phong_distribution
// NOTE: Only valid if dot(normalize(normal), normalize(half)) > 0
inline float phong_isotropic_highlight( const frantic::graphics::vector3f& normal, // normalize(normal)
                                        const frantic::graphics::vector3f& half,   // normalize(toLight + toEye)
                                        float glossiness ) {
    using frantic::graphics::vector3f;

    float NdotH = vector3f::dot( normal, half );
    return ( NdotH > 0 ) ? std::pow( NdotH, glossiness ) : 0.f;
}

// This is an anisotropic specular highlight function created by Ward.
// See: http://en.wikipedia.org/wiki/Specular_highlight#Ward_anisotropic_distribution
// NOTE: Only valid if NDotL > 0 && NDotV > 0.
inline float ward_anisotropic_highlight( const frantic::graphics::vector3f& normal,  // normalize(normal)
                                         const frantic::graphics::vector3f& tangent, // normalize(tangent)
                                         const frantic::graphics::vector3f& half,    // normalize(toLight + toEye)
                                         float NdotL, // dot(normalize(Normal), normalize(toLight))
                                         float NdotV, // dot(normalize(Normal), normalize(toEye))
                                         float tRoughness, float bRoughness ) {
    using frantic::graphics::vector3f;

    vector3f binormal = vector3f::normalize( vector3f::cross( normal, tangent ) );

    float TdotH = vector3f::dot( tangent, half );
    float BdotH = vector3f::dot( binormal, half );
    float NdotH = vector3f::dot( normal, half );
    float rho =
        std::exp( -2 * ( math::square( TdotH / tRoughness ) + math::square( BdotH / bRoughness ) ) / ( 1.f + NdotH ) );

    return rho * NdotL / ( 4 * tRoughness * bRoughness * std::sqrt( NdotL * NdotV ) );
}

} // namespace shading
} // namespace frantic
