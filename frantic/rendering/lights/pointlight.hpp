// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/motion_blurred_transform.hpp>
#include <frantic/rendering/lights/lightinterface.hpp>

namespace frantic {
namespace rendering {
namespace lights {

class pointlight : public lightinterface {
  protected:
    frantic::graphics::motion_blurred_transform<float> m_transform;
    frantic::graphics2d::size2 m_shadowMapSize;
    // TODO: I don't like the name "Shadow Density" for this, it should be something more like "shadow scaling factor"
    // or some such
    float m_shadowDensity;

  public:
    pointlight( const frantic::tstring& name, const frantic::graphics::motion_blurred_transform<float>& xform,
                const frantic::graphics::color3f& flux, int decayExponent = 2, bool bShadowsEnabled = false,
                float shadowDensity = 1.f, int shadowMapWidth = 512, bool bUseNearAtten = false,
                bool bUseFarAtten = false, float nearAttenuationStart = 0.f, float nearAttenuationEnd = 0.f,
                float farAttenuationStart = std::numeric_limits<float>::max(),
                float farAttenuationEnd = std::numeric_limits<float>::max(), float decayOffset = 1.f

                )
        : lightinterface( name, flux, decayExponent, bShadowsEnabled, bUseNearAtten, bUseFarAtten, nearAttenuationStart,
                          nearAttenuationEnd, farAttenuationStart, farAttenuationEnd, decayOffset )
        , m_transform( xform )
        , m_shadowMapSize( frantic::graphics2d::size2( shadowMapWidth, shadowMapWidth ) )
        , m_shadowDensity( shadowDensity ) {}

    virtual ~pointlight() {}

    virtual bool is_directional_light() const { return false; }
    virtual bool is_area_light() const { return false; }
    virtual void write_xml( std::basic_ostream<frantic::tchar>& out, const frantic::tstring& prefix ) const {
        out << prefix << _T("<type>pointlight</type>\n");
        lightinterface::write_xml( out, prefix );
    }

    virtual frantic::graphics2d::size2 shadow_map_size() const { return m_shadowMapSize; }
    virtual float shadow_density() const { return m_shadowDensity; }
    virtual const frantic::graphics::transform4f& transform_matrix() const { return m_transform.get_transform(); }
    virtual const frantic::graphics::vector3f& position() const { return m_transform.get_transform().get_column( 3 ); }

    virtual frantic::graphics::transform4f transform_matrix( float motionSegmentTime ) const {
        return m_transform.get_transform( motionSegmentTime );
    }
    virtual frantic::graphics::vector3f position( float motionSegmentTime ) const {
        return m_transform.get_transform( motionSegmentTime ).get_column( 3 );
    }
    virtual frantic::graphics::vector3f position( const frantic::graphics::vector3f& /*relativeTo*/ ) const {
        return m_transform.get_transform().get_column( 3 );
    }
    virtual frantic::graphics::vector3f position( float motionSegmentTime,
                                                  const frantic::graphics::vector3f& /*relativeTo*/ ) const {
        return m_transform.get_transform( motionSegmentTime ).get_column( 3 );
    }

    // Calculates the irradiance from the light at a given point in space
    // bIsValid remains true if the light shines on the given shadingPosition.  For the pointlight, it's always true.
    virtual frantic::graphics::color3f irradiance( const frantic::graphics::vector3f& shadingPosition,
                                                   float motionSegmentTime, bool& /*bIsValid*/ ) const {
        float falloffWeight = 1.f;
        float distSquared =
            frantic::graphics::vector3f::distance_squared( position( motionSegmentTime ), shadingPosition );
        float dist = std::sqrt( distSquared );
        if( use_near_attenuation() ) {
            if( dist < m_nearAttenuationStart ) {
                return frantic::graphics::color3f();
            } else if( dist < m_nearAttenuationEnd )
                falloffWeight *= frantic::math::smoothstep( dist, m_nearAttenuationStart, m_nearAttenuationEnd );
        }

        if( use_far_attenuation() ) {
            if( dist > m_farAttenuationEnd ) {
                return frantic::graphics::color3f();
            } else if( dist > m_farAttenuationStart )
                falloffWeight *= 1.f - frantic::math::smoothstep( dist, m_farAttenuationStart, m_farAttenuationEnd );
        }

        const frantic::graphics::color3f weightedFlux = m_flux_over_4pi * falloffWeight;

        const float divisor = dist < m_decayOffset - 1.f ? m_decayOffset : dist;
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
    // This function doesn't touch bIsValid, because the pointlight always has output everywhere.
    virtual frantic::graphics::color3f radiance( const frantic::graphics::vector3f& shadingPosition,
                                                 float motionSegmentTime, const frantic::graphics::vector3f& normal,
                                                 bool& /*bIsValid*/ ) const {
        frantic::graphics::vector3f vectorToLight = position( motionSegmentTime ) - shadingPosition;
        float distanceSquared = vectorToLight.get_magnitude_squared();
        float distance = std::sqrt( distanceSquared );
        float cosTheta = frantic::graphics::vector3f::dot( normal, vectorToLight ) / distance;
        if( cosTheta > 0 ) {
            if( m_decayExponent == 3 ) {
                return m_flux_over_4pi * cosTheta / ( distanceSquared * distance );
            } else if( m_decayExponent == 2 ) {
                return m_flux_over_4pi * cosTheta / distanceSquared;
            } else if( m_decayExponent == 1 ) {
                return m_flux_over_4pi * cosTheta / distance;
            } else {
                return m_flux_over_4pi * cosTheta;
            }
        } else
            return frantic::graphics::color3f();
    }
};

} // namespace lights
} // namespace rendering
} // namespace frantic
