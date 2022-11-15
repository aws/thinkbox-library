// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/camera.hpp>
#include <frantic/rendering/lights/lightinterface.hpp>

namespace frantic {
namespace rendering {
namespace lights {

enum LIGHT_SHAPE { LIGHT_SHAPE_SQUARE, LIGHT_SHAPE_ROUND };

class directedlight_interface : public lightinterface {
  protected:
    frantic::graphics::camera<float> m_camera;

    // TODO: I don't like the name "Shadow Density" for this, it should be something more like "shadow scaling factor"
    // or some such
    float m_shadowDensity;

    LIGHT_SHAPE m_lightShape;
    float m_lightAspect;
    float m_sqrtLightAspect;
    float m_innerRadius;
    float m_outerRadius;

  public:
    directedlight_interface(
        // lightinterface
        const frantic::tstring& name, const frantic::graphics::color3f& flux, int decayExponent, bool bShadowsEnabled,
        float shadowDensity,

        bool bUseNearAtten, bool bUseFarAtten, float nearAttenuationStart, float nearAttenuationEnd,
        float farAttenuationStart, float farAttenuationEnd,

        // directedlight_interface
        const frantic::graphics::camera<float>& camera, LIGHT_SHAPE lightShape, float lightAspect, float innerRadius,
        float outerRadius, float decayOffset = 1.f )
        : lightinterface( name, flux, decayExponent, bShadowsEnabled, bUseNearAtten, bUseFarAtten, nearAttenuationStart,
                          nearAttenuationEnd, farAttenuationStart, farAttenuationEnd, decayOffset )
        , m_camera( camera )
        , m_lightShape( lightShape )
        , m_lightAspect( lightAspect )
        , m_sqrtLightAspect( sqrt( lightAspect ) )
        , m_innerRadius( innerRadius )
        , m_outerRadius( outerRadius )
        , m_shadowDensity( shadowDensity )

    {}

    virtual ~directedlight_interface() {}

    /**
     * From lightinterface
     */

    // Some basic properties that external processors can use to determine how this light will function
    virtual bool is_directional_light() const { return true; }
    virtual bool is_area_light() const { return false; }
    virtual frantic::graphics2d::size2 shadow_map_size() const { return m_camera.output_size(); }
    virtual float shadow_density() const { return m_shadowDensity; }

    virtual const frantic::graphics::vector3f& position() const { return m_camera.camera_position(); }
    virtual const frantic::graphics::transform4f& transform_matrix() const { return m_camera.world_transform(); }

    virtual frantic::graphics::vector3f position( float motionSegmentTime ) const {
        return m_camera.camera_position( motionSegmentTime );
    }
    virtual frantic::graphics::transform4f transform_matrix( float motionSegmentTime ) const {
        return m_camera.world_transform( motionSegmentTime );
    }
    virtual frantic::graphics::vector3f position( const frantic::graphics::vector3f& relativeTo ) const {
        return m_camera.camera_position( relativeTo );
    }
    virtual frantic::graphics::vector3f position( float motionSegmentTime,
                                                  const frantic::graphics::vector3f& relativeTo ) const {
        return m_camera.camera_position( motionSegmentTime, relativeTo );
    }
    virtual void write_xml( std::basic_ostream<frantic::tchar>& out, const frantic::tstring& prefix ) const {
        lightinterface::write_xml( out, prefix );
        out << prefix << _T("<lightshape>") << m_lightShape << _T("</lightshape>\n");
        out << prefix << _T("<lightaspect>") << m_lightAspect << _T("</lightaspect>\n");
        out << prefix << _T("<innerradius>") << m_innerRadius << _T("</innerradius>\n");
        out << prefix << _T("<outerradius>") << m_outerRadius << _T("</outerradius>\n");
    }

    /**
     * directedlight_interface specifc functions
     */
    virtual frantic::graphics::vector3f direction() const { return m_camera.view_direction(); }
    virtual frantic::graphics::vector3f direction( float motionSegmentTime ) const {
        return m_camera.view_direction( motionSegmentTime );
    }
    virtual const frantic::graphics::camera<float>& get_camera() const { return m_camera; }
};

} // namespace lights
} // namespace rendering
} // namespace frantic
