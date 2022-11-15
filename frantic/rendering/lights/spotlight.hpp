// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/motion_blurred_transform.hpp>
#include <frantic/rendering/lights/directedlight_interface.hpp>

namespace frantic {
namespace rendering {
namespace lights {

class spotlight : public directedlight_interface {
  public:
    spotlight(
        // lightinterface
        const frantic::tstring& name, const frantic::graphics::motion_blurred_transform<float>& xform,
        const frantic::graphics::color3f& flux, int decayExponent, bool bShadowsEnabled, float shadowDensity,
        int shadowMapWidth,

        bool bUseNearAtten, bool bUseFarAtten, float nearAttenuationStart, float nearAttenuationEnd,
        float farAttenuationStart, float farAttenuationEnd,

        // directedlight_interface
        LIGHT_SHAPE lightShape, float lightAspect,
        float innerConeHalfAngle, // IMPORTANT: do not pass in the whole angle, half of it only.
        float outerConeHalfAngle, // IMPORTANT: do not pass in the whole angle, half of it only.
        float decayOffset = 1.f )
        : directedlight_interface( name, flux, decayExponent, bShadowsEnabled, shadowDensity,

                                   bUseNearAtten, bUseFarAtten, nearAttenuationStart, nearAttenuationEnd,
                                   farAttenuationStart, farAttenuationEnd,

                                   frantic::graphics::camera<float>(
                                       frantic::graphics::projection_mode::perspective, xform,
                                       2.0f * ::atanf( ::tanf( outerConeHalfAngle ) * sqrt( lightAspect ) ),
                                       frantic::graphics2d::size2( shadowMapWidth ), 0.001f, 1e+10, lightAspect ),
                                   lightShape, lightAspect, innerConeHalfAngle, outerConeHalfAngle, decayOffset ) {}

    virtual ~spotlight() {}

    virtual void write_xml( std::basic_ostream<frantic::tchar>& out, const frantic::tstring& prefix ) const {
        out << prefix << _T("<type>spotlight</type>\n");
        directedlight_interface::write_xml( out, prefix );
    }

    // Calculates the irradiance from the light at a given point in space
    // bIsValid is true if the light shines on the given shadingPosition
    virtual frantic::graphics::color3f irradiance( const frantic::graphics::vector3f& shadingPosition,
                                                   float motionSegmentTime, bool& bIsValid ) const {
        using frantic::graphics::color3f;
        using frantic::graphics::vector3f;

        vector3f camSpacePos = get_camera().world_transform_inverse( motionSegmentTime ) * shadingPosition;
        float rayDist = camSpacePos.get_magnitude();
        float angle;

        if( m_lightShape == LIGHT_SHAPE_ROUND ) {
            angle = acosf( -camSpacePos.z / rayDist );
        } else {
            // LIGHT_SHAPE_SQUARE
            // This code normalizes the angles based on the light aspect. Effectively, this
            // divides out the aspect so that the x, and y angles can be compared to m_outerRadius
            // directly.
            angle = std::max( atanf( fabsf( camSpacePos.x / camSpacePos.z / m_sqrtLightAspect ) ),
                              atanf( fabsf( camSpacePos.y / camSpacePos.z / m_sqrtLightAspect * m_lightAspect ) ) );
        }

        if( angle > m_outerRadius || camSpacePos.z > 0 ) {
            bIsValid = false;
            return frantic::graphics::color3f();
        }

        float falloffWeight = 1.f;
        if( angle > m_innerRadius )
            falloffWeight -= frantic::math::smoothstep( angle, m_innerRadius, m_outerRadius );

        if( use_near_attenuation() ) {
            if( rayDist < m_nearAttenuationStart ) {
                bIsValid = false;
                return color3f();
            } else if( rayDist < m_nearAttenuationEnd )
                falloffWeight *= frantic::math::smoothstep( rayDist, m_nearAttenuationStart, m_nearAttenuationEnd );
        }

        if( use_far_attenuation() ) {
            if( rayDist > m_farAttenuationEnd ) {
                bIsValid = false;
                return frantic::graphics::color3f();
            } else if( rayDist > m_farAttenuationStart )
                falloffWeight *= 1.f - frantic::math::smoothstep( rayDist, m_farAttenuationStart, m_farAttenuationEnd );
        }

        const frantic::graphics::color3f weightedFlux = m_flux_over_4pi * falloffWeight;

        const float divisor = rayDist < m_decayOffset - 1.f ? m_decayOffset : rayDist;
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
        using frantic::graphics::color3f;
        using frantic::graphics::vector3f;

        vector3f vectorToLight = directedlight_interface::position( motionSegmentTime ) - shadingPosition;

        float cosTheta = vector3f::dot( normal, vectorToLight ) / vectorToLight.get_magnitude();
        return ( cosTheta > 0 ) ? cosTheta * irradiance( shadingPosition, motionSegmentTime, bIsValid ) : color3f();
    }
};

} // namespace lights
} // namespace rendering
} // namespace frantic
