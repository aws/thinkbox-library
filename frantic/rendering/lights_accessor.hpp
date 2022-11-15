// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/rendering/lights/lightinterface.hpp>
#include <frantic/rendering/shadow_generator.hpp>

// The lights_accessor class is intended to serve as a wrapper around light list access.
// The idea is that by creating a virtual interface, we will be able to use the same code
// for rendering within 3ds max, Brazil, mental ray, Gelato or any other render in which
// we want to put our code.

namespace frantic {
namespace rendering {
// using frantic::graphics::color3f;
// using frantic::graphics::vector3f;

class lights_accessor {
  public:
    virtual ~lights_accessor() {}

    virtual int size() const = 0;

    virtual lights::lightinterface& light( int index ) = 0;
    virtual const lights::lightinterface& light( int index ) const = 0;
    virtual void illuminate( int lightIndex, const frantic::graphics::vector3f& position,
                             const frantic::graphics::vector3f& normal, frantic::graphics::color3f& outCl ) const = 0;

    /*
    void initialize_shadows( const shadow_combiner& globalShadows, boost::shared_ptr<const
    frantic::geometry::raytraced_geometry_collection>& geometry ) { std::cerr << "Initializing shadows for " << size()
    << " lights" << std::endl; for( int i = 0; i < size(); ++i ) {
        //light(i).initialize_shadows( globalShadows, geometry );	is not implemented for any lights
      }
    }
    */
};

/* //Doesn't compile .... I have no idea what shadows() was supposed to do!
template< class lights_accessor_class, class volume_scattering_function >
frantic::graphics::color3f volumetric_inscattered_light( lights_accessor_class& lights, const
frantic::graphics::vector3f& samplePoint, const frantic::graphics::vector3f& rayDirection, volume_scattering_function
vsf ) { using namespace frantic::graphics;

  color3f result;

  int lightCount = lights.size();
  for( int lightIndex = 0; lightIndex < lightCount; ++lightIndex ) {
    vector3f lightPosition = lights.light_position( lightIndex );

    float visibility = shadows().visibility( samplePoint, lightPosition );

    if( visibility > 0 ) {

      vector3f sampleToLightDirection = lightPosition - samplePoint;
      sampleToLightDirection.normalize();

      color3f lightColor;
      if( lights.illuminate( lightIndex, samplePoint, sampleToLightDirection, HALFPI, lightColor ) ) {
        // The illuminate call incorporated all shadow calculation, so no need to do that here.

        // todo: this prolly does NOT belong here
        color3f phase = vsf( rayDirection, sampleToLightDir );
        // todo: prolly wrong.
        phase *= 4.0f * PI;
        result += visibility * phase * lightColor;
      }
    }
  }

  return result;
}
*/

} // namespace rendering
} // namespace frantic
