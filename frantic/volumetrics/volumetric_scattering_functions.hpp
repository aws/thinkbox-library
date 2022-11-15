// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace volumetrics {

using frantic::graphics::color3f;
using frantic::graphics::vector3f;

inline float isotropic_scattering() { return 1 / ( 4.f * (float)M_PI ); }

// Henyey-Greenstein phase function.
// The incoming and outgoing vectors should be normalized.
inline float phase_hg( const vector3f& w, const vector3f& wp, float g ) {
    float costheta = vector3f::dot( w, wp );
    return 1.0f / ( 4.0f * (float)M_PI ) * ( 1.0f - g * g ) / powf( 1.0f + g * g - 2.0f * g * costheta, 1.5f );
}

inline color3f phase_hg( const vector3f& w, const vector3f& wp, const color3f& g ) {
    return color3f( phase_hg( w, wp, g.r ), phase_hg( w, wp, g.g ), phase_hg( w, wp, g.b ) );
}

// Schlick phase function (approximates Henyey-Greenstein)
// The incoming and outgoing vectors should be normalized.
inline float phase_schlick( const vector3f& w, const vector3f& wp, float g ) {
    float k = 1.55f * g - .55f * g * g * g;
    float kcostheta = k * vector3f::dot( w, wp );
    return 1.0f / ( 4.0f * (float)M_PI ) * ( 1.0f - k * k ) / ( ( 1.0f - kcostheta ) * ( 1.0f - kcostheta ) );
}

inline color3f phase_schlick( const vector3f& w, const vector3f& wp, const color3f& g ) {
    return color3f( phase_schlick( w, wp, g.r ), phase_schlick( w, wp, g.g ), phase_schlick( w, wp, g.b ) );
}

} // namespace volumetrics
} // namespace frantic
