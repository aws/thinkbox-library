// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/motion_blurred_transform.hpp>
#include <frantic/rendering/lights/directedlight_interface.hpp>

namespace frantic {
namespace rendering {
namespace lights {

class directlight : public directedlight_interface {
  public:
    // Constructor from XML
    directlight(
        // lightinterface
        const frantic::tstring& name, const frantic::graphics::motion_blurred_transform<float>& xform,
        const frantic::graphics::color3f& flux, int decayExponent, bool bShadowsEnabled, float shadowDensity,
        int shadowMapWidth,

        bool bUseNearAtten, bool bUseFarAtten, float nearAttenuationStart, float nearAttenuationEnd,
        float farAttenuationStart, float farAttenuationEnd,

        // directedlight_interface
        LIGHT_SHAPE lightShape, float lightAspect, float innerRectRadius, float outerRectRadius,
        float decayOffset = 1.f )
        : directedlight_interface( name, flux, decayExponent, bShadowsEnabled, shadowDensity,

                                   bUseNearAtten, bUseFarAtten, nearAttenuationStart, nearAttenuationEnd,
                                   farAttenuationStart, farAttenuationEnd,

                                   frantic::graphics::camera<float>(
                                       frantic::graphics::projection_mode::orthographic, xform, 2 * outerRectRadius,
                                       frantic::graphics2d::size2( shadowMapWidth ), 0.001f, 1e+10, lightAspect ),
                                   lightShape, lightAspect, innerRectRadius, outerRectRadius, decayOffset ) {}

    virtual ~directlight() {}

    virtual void write_xml( std::basic_ostream<frantic::tchar>& out, const frantic::tstring& prefix ) const {
        out << prefix << _T("<type>directlight</type>\n");
        directedlight_interface::write_xml( out, prefix );
    }

    // Calculates the irradiance from the light at a given point in space
    // bIsValid is true if the light shines on the given shadingPosition
    virtual frantic::graphics::color3f irradiance( const frantic::graphics::vector3f& shadingPosition,
                                                   float motionSegmentTime, bool& bIsValid ) const {
        using frantic::graphics::vector3f;
        vector3f posInCamera = get_camera().world_transform_inverse( motionSegmentTime ) * shadingPosition;
        float parallelDist = -posInCamera.z;
        float perpDist;

        if( m_lightShape == LIGHT_SHAPE_ROUND )
            perpDist = sqrtf( posInCamera.x * posInCamera.x + posInCamera.y * posInCamera.y );
        else {
            // LIGHT_SHAPE_SQUARE
            perpDist = std::max( fabsf( posInCamera.x ), fabsf( posInCamera.y ) * m_lightAspect );
        }

        // Return black for points outside the beam or behind(positive z) the light position.
        if( perpDist > m_outerRadius || parallelDist < 0 ) {
            bIsValid = false;
            return frantic::graphics::color3f();
        }

        float multiplier = 1.f;
        if( perpDist > m_innerRadius )
            multiplier -= frantic::math::smoothstep( perpDist, m_innerRadius, m_outerRadius );

        if( use_near_attenuation() ) {
            if( parallelDist < m_nearAttenuationStart ) {
                bIsValid = false;
                return frantic::graphics::color3f();
            } else if( parallelDist < m_nearAttenuationEnd )
                multiplier *= frantic::math::smoothstep( parallelDist, m_nearAttenuationStart, m_nearAttenuationEnd );
        }

        if( use_far_attenuation() ) {
            if( parallelDist > m_farAttenuationEnd ) {
                bIsValid = false;
                return frantic::graphics::color3f();
            } else if( parallelDist > m_farAttenuationStart )
                multiplier *=
                    1.f - frantic::math::smoothstep( parallelDist, m_farAttenuationStart, m_farAttenuationEnd );
        }

        const frantic::graphics::color3f weightedFlux = m_flux_over_4pi * multiplier;

        const float divisor = parallelDist < m_decayOffset - 1.f ? m_decayOffset : parallelDist;
        switch( m_decayExponent ) {
        case 3:
            return weightedFlux / ( divisor * divisor * divisor );
        case 2:
            return weightedFlux / ( divisor * divisor );
        case 1:
            return weightedFlux / divisor;
        default:
            return weightedFlux;
        }
    }

    // Calcualates the irradiance on a surface at a given point with a given normal
    // This is irradiance(position) * cos(angle between normal and direction to light)
    virtual frantic::graphics::color3f radiance( const frantic::graphics::vector3f& shadingPosition,
                                                 float motionSegmentTime, const frantic::graphics::vector3f& normal,
                                                 bool& bIsValid ) const {
        float multiplier = frantic::graphics::vector3f::dot(
            normal, get_camera().world_transform( motionSegmentTime ).get_column( 3 ) /*-viewDirection*/ );
        return ( multiplier > 0 ) ? ( multiplier * irradiance( shadingPosition, motionSegmentTime, bIsValid ) )
                                  : frantic::graphics::color3f();
    }
};

} // namespace lights
} // namespace rendering
} // namespace frantic
