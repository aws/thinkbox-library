// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

#include <boost/shared_ptr.hpp>

#include <frantic/rendering/lights/directlight.hpp>
#include <frantic/rendering/lights/lightinterface.hpp>
#include <frantic/rendering/lights/pointlight.hpp>
#include <frantic/rendering/lights/spotlight.hpp>
#include <frantic/rendering/lights_accessor.hpp>

#ifdef MAX_VERSION
#include <frantic/max3d/parameter_extraction.hpp>
#define MR_OMNI_LIGHT_CLASS_ID Class_ID( 112233, 554433 )
#define MR_SPOT_LIGHT_CLASS_ID Class_ID( 112233, 554434 )
#endif

namespace frantic {
namespace rendering {
namespace lights {

class light_list : public frantic::rendering::lights_accessor {
    std::vector<boost::shared_ptr<lightinterface>> m_lights;

  public:
    void add_light( const boost::shared_ptr<lightinterface>& light ) { m_lights.push_back( light ); }

#ifdef MAX_VERSION
    // All 3ds Max standard lights are subclassed from GenLight, and in order to get at the decay
    // parameters we need to be able to know if we have a GenLight. The SDK does not offer the ability to
    // do so, so this will have to suffice.
    bool IsSubClassOf_GenLight( Object* obj ) {
        return obj->IsSubClassOf( Class_ID( OMNI_LIGHT_CLASS_ID, 0x0 ) ) ||
               obj->IsSubClassOf( Class_ID( SPOT_LIGHT_CLASS_ID, 0x0 ) ) ||
               obj->IsSubClassOf( Class_ID( FSPOT_LIGHT_CLASS_ID, 0x0 ) ) ||
               obj->IsSubClassOf( Class_ID( DIR_LIGHT_CLASS_ID, 0x0 ) ) ||
               obj->IsSubClassOf( Class_ID( TDIR_LIGHT_CLASS_ID, 0x0 ) ) ||
               obj->IsSubClassOf( MR_OMNI_LIGHT_CLASS_ID ) || obj->IsSubClassOf( MR_SPOT_LIGHT_CLASS_ID );
    }

    void add_light( INode* lightNode, TimeValue t, float mblurInterval = 0.5f, float mblurBias = 0.f ) {
        using namespace frantic::graphics;

        Object* obj = lightNode->GetObjectRef();
        if( !obj || obj->SuperClassID() != LIGHT_CLASS_ID )
            return;

        LightObject* lightObj = static_cast<LightObject*>( obj );
        LightState lightState;

        Interval ivalid = FOREVER;
        if( lightObj->EvalLightState( t, ivalid, &lightState ) != REF_SUCCEED )
            return;
        if( !lightState.on || lightState.intens <= 0.f )
            return;

        int decayExponent = 2; // Default based off of real world.
        float unitDecayDistance = 40.f;

        if( IsSubClassOf_GenLight( lightObj ) ) {
            decayExponent = ( (GenLight*)obj )->GetDecayType();
            unitDecayDistance = ( (GenLight*)obj )->GetDecayRadius( t );
        }

        float shadowDensity = max3d::mxs::expression( _T("try(obj.shadowMultiplier)catch(1.0)") )
                                  .bind( _T("obj"), lightObj )
                                  .at_time( t )
                                  .evaluate<float>();
        int shadowMapSize = max3d::mxs::expression( _T("try(obj.mapsize)catch(512)") )
                                .bind( _T("obj"), lightObj )
                                .at_time( t )
                                .evaluate<int>();

        float multiplier = lightState.intens * 4 * (float)M_PI;
        if( decayExponent == 1 )
            multiplier *= unitDecayDistance;
        else if( decayExponent == 2 )
            multiplier *= unitDecayDistance * unitDecayDistance;

        frantic::graphics::color3f color = frantic::graphics::color3f( lightState.color ) * multiplier;
        motion_blurred_transform<float> lightTransform( lightNode, t, mblurInterval, mblurBias, NULL );

        frantic::tstring name = lightNode->GetName();
        bool shadowsEnabled = ( lightState.shadow != FALSE );

        // for near/far threshold light attenuation
        bool useNearAtten = ( lightState.useNearAtten != FALSE );
        bool useFarAtten = ( lightState.useAtten != FALSE );
        float nearAttenuationStart = lightState.nearAttenStart;
        float nearAttenuationEnd = lightState.nearAttenEnd;
        float farAttenuationStart = lightState.attenStart;
        float farAttenuationEnd = lightState.attenEnd;
        if( farAttenuationEnd < farAttenuationStart )
            farAttenuationEnd = farAttenuationStart;
        if( nearAttenuationEnd < nearAttenuationStart )
            nearAttenuationEnd = nearAttenuationStart;

        switch( lightState.type ) {
        case OMNI_LGT: {
            boost::shared_ptr<lightinterface> newLight( new pointlight(
                name, lightTransform, color, decayExponent, shadowsEnabled, shadowDensity, shadowMapSize, useNearAtten,
                useFarAtten, nearAttenuationStart, nearAttenuationEnd, farAttenuationStart, farAttenuationEnd ) );
            add_light( newLight );

            FF_LOG( debug ) << "Added omnilight: " << name << std::endl;
            break;
        }
        case SPOT_LGT: {
            lights::LIGHT_SHAPE shape = (lights::LIGHT_SHAPE)lightState.shape;
            float aspect = ( shape == LIGHT_SHAPE_SQUARE ) ? lightState.aspect : 1.f;
            float innerRadius = math::degrees_to_radians(
                lightState.hotsize / 2 ); // divide by two since the constructor takes the half angle
            float outerRadius = math::degrees_to_radians(
                lightState.fallsize / 2 ); // divide by two since the constructor takes the half angle

            if( outerRadius < innerRadius )
                outerRadius = innerRadius;

            boost::shared_ptr<lightinterface> newLight(
                new spotlight( name, lightTransform, color, decayExponent, shadowsEnabled, shadowDensity, shadowMapSize,
                               useNearAtten, useFarAtten, nearAttenuationStart, nearAttenuationEnd, farAttenuationStart,
                               farAttenuationEnd, shape, aspect, innerRadius, outerRadius ) );
            add_light( newLight );

            FF_LOG( debug ) << "Added spotlight: " << name << std::endl;
            break;
        }
        case DIRECT_LGT: {
            lights::LIGHT_SHAPE shape = (lights::LIGHT_SHAPE)lightState.shape;
            float aspect = ( shape == LIGHT_SHAPE_SQUARE ) ? lightState.aspect : 1.f;
            float innerRadius = lightState.hotsize;
            float outerRadius = lightState.fallsize;

            if( outerRadius < innerRadius )
                outerRadius = innerRadius;

            boost::shared_ptr<lightinterface> newLight(
                new directlight( name, lightTransform, color, decayExponent, shadowsEnabled, shadowDensity,
                                 shadowMapSize, useNearAtten, useFarAtten, nearAttenuationStart, nearAttenuationEnd,
                                 farAttenuationStart, farAttenuationEnd, shape, aspect, innerRadius, outerRadius ) );
            add_light( newLight );

            FF_LOG( debug ) << "Added directlight: " << name << std::endl;
            break;
        }
        default:
            // mprintf("Light: %s has an unsupported type.\n", name.c_str());
            break;
        }
    }
#endif // 3ds max headers

    int size() const { return (int)m_lights.size(); }

    lightinterface& light( int index ) { return *m_lights[index].get(); }

    const lightinterface& light( int index ) const { return *m_lights[index].get(); }

    boost::shared_ptr<lightinterface> light_ptr( int index ) { return m_lights[index]; }

    void illuminate( int lightIndex, const frantic::graphics::vector3f& position,
                     const frantic::graphics::vector3f& normal, frantic::graphics::color3f& outCl ) const {
        bool bResult = true;
        outCl = m_lights[lightIndex]->radiance( position, normal, bResult );
    }

    void illuminate( int lightIndex, const frantic::graphics::vector3f& position,
                     frantic::graphics::color3f& outCl ) const {
        bool bResult = true;
        outCl = m_lights[lightIndex]->irradiance( position, bResult );
    }
};
} // namespace lights
} // namespace rendering
} // namespace frantic
