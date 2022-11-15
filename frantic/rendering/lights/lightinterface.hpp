// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/size2.hpp>
#include <frantic/math/radian.hpp>
#include <frantic/rendering/deep_attenuation_loaders.hpp>
#include <frantic/rendering/deep_attenuation_savers.hpp>

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace rendering {
namespace lights {

class lightinterface {
  protected:
    frantic::tstring m_name;
    bool m_bShadowsEnabled;
    bool m_bUseNearAtten;
    bool m_bUseFarAtten;

    frantic::graphics::color3f m_flux_over_4pi;
    int m_decayExponent;

    // values for the near/far light attenuation that can be used with any light
    float m_nearAttenuationStart;
    float m_nearAttenuationEnd;
    float m_farAttenuationStart;
    float m_farAttenuationEnd;

    // interface to save light attenuation to a file.
    boost::shared_ptr<frantic::rendering::atten_saver> m_attenSaver;

    // interface to include pre-rendered light attenuation from a file.
    boost::shared_ptr<frantic::rendering::atten_loader> m_attenLoader;

    float m_decayOffset;

  public:
    lightinterface( const frantic::tstring& name, const frantic::graphics::color3f& flux, int decayExponent,
                    bool bShadowsEnabled, bool bUseNearAtten = false, bool bUseFarAtten = false,
                    float nearAttenuationStart = 0.f, float nearAttenuationEnd = 0.f,
                    float farAttenuationStart = std::numeric_limits<float>::max(),
                    float farAttenuationEnd = std::numeric_limits<float>::max(), float decayOffset = 1.f )
        : m_name( name )
        , m_flux_over_4pi( flux / ( 2.f * (float)M_2PI ) )
        , m_decayExponent( decayExponent )
        , m_bShadowsEnabled( bShadowsEnabled )
        , m_bUseNearAtten( bUseNearAtten )
        , m_bUseFarAtten( bUseFarAtten )
        , m_nearAttenuationStart( nearAttenuationStart )
        , m_nearAttenuationEnd( nearAttenuationEnd )
        , m_farAttenuationStart( farAttenuationStart )
        , m_farAttenuationEnd( farAttenuationEnd )
        , m_decayOffset( decayOffset ) {}

    virtual ~lightinterface() {}

    virtual const frantic::tstring& name() const { return m_name; }
    virtual bool is_shadows_enabled() const { return m_bShadowsEnabled; }
    virtual bool use_near_attenuation() const { return m_bUseNearAtten; }
    virtual bool use_far_attenuation() const { return m_bUseFarAtten; }

    // Output attenuation saving stuff
    virtual void set_attenuation_saver( boost::shared_ptr<frantic::rendering::atten_saver> attenSaver ) {
        m_attenSaver = attenSaver;
    }
    virtual boost::shared_ptr<frantic::rendering::atten_saver> get_attenuation_saver() { return m_attenSaver; }

    // Input attenuation loading stuff
    virtual void set_attenuation_loader( boost::shared_ptr<frantic::rendering::atten_loader> attenLoader ) {
        m_attenLoader = attenLoader;
    }
    virtual boost::shared_ptr<frantic::rendering::atten_loader> get_attenuation_loader() { return m_attenLoader; }

    virtual void write_xml( std::basic_ostream<frantic::tchar>& out, const frantic::tstring& prefix ) const {
        out << prefix << _T("<name>") << m_name << _T("</name>\n");
        out << prefix << _T("<transform>") << transform_matrix() << _T("</transform>\n");
        out << prefix << _T("<flux>") << ( m_flux_over_4pi * 2.f * (float)M_2PI ) << _T("</flux>\n");
        out << prefix << _T("<decayexponent>") << m_decayExponent << _T("</decayexponent>\n");
        out << prefix << _T("<castshadows>") << ( m_bShadowsEnabled ? _T("true") : _T("false") )
            << _T("</castshadows>\n");
        out << prefix << _T("<shadowmapwidth>") << shadow_map_size().xsize << _T("</shadowmapwidth>\n");
        // TODO: I don't like the name "Shadow Density" for this, it should be something more like "shadow scaling
        // factor" or some such
        out << prefix << _T("<shadowdensity>") << shadow_density() << _T("</shdowdensity>\n");
    }

    ///---------------------///
    /// Interface functions ///
    ///---------------------///
    // Some basic properties that external processors can use to determine how this light will function
    virtual bool is_directional_light() const = 0;
    virtual bool is_area_light() const = 0;

    virtual frantic::graphics2d::size2 shadow_map_size() const = 0;
    // TODO: I don't like the name "Shadow Density" for this, it should be something more like "shadow scaling factor"
    // or some such
    virtual float shadow_density() const = 0;

    virtual const frantic::graphics::transform4f& transform_matrix() const = 0;
    virtual const frantic::graphics::vector3f& position() const = 0;

    virtual frantic::graphics::transform4f transform_matrix( float motionSegmentTime ) const = 0;
    virtual frantic::graphics::vector3f position( float motionSegmentTime ) const = 0;
    virtual frantic::graphics::vector3f position( const frantic::graphics::vector3f& relativeTo ) const = 0;
    virtual frantic::graphics::vector3f position( float motionSegmentTime,
                                                  const frantic::graphics::vector3f& relativeTo ) const = 0;

    /**
     * Calculates the irradiance from the light at a given point in world space
     * @note bIsValid is true if the light shines on the given shadingPosition
     */
    virtual frantic::graphics::color3f irradiance( const frantic::graphics::vector3f& shadingPosition,
                                                   float motionSegmentTime, bool& bIsValid ) const = 0;

    /**
     * @overload
     */
    frantic::graphics::color3f irradiance( const frantic::graphics::vector3f& shadingPosition, bool& bIsValid ) const {
        return this->irradiance( shadingPosition, 0.5f, bIsValid );
    }

    /**
     * Calcualates the irradiance on a surface at a given point with a given normal
     * This is irradiance(position) * cos(angle between normal and direction to light)
     */
    virtual frantic::graphics::color3f radiance( const frantic::graphics::vector3f& shadingPosition,
                                                 float motionSegmentTime, const frantic::graphics::vector3f& normal,
                                                 bool& bIsValid ) const = 0;

    /**
     * @overload
     */
    frantic::graphics::color3f radiance( const frantic::graphics::vector3f& shadingPosition,
                                         const frantic::graphics::vector3f& normal, bool& bIsValid ) const {
        return this->radiance( shadingPosition, 0.5f, normal, bIsValid );
    }
};

} // namespace lights
} // namespace rendering
} // namespace frantic
