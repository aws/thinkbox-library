// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics/camera.hpp>
#include <frantic/graphics/color_with_alpha.hpp>
#include <frantic/graphics2d/framebuffer.hpp>
#include <frantic/rendering/framebuffer_cubeface.hpp>

namespace frantic {
namespace rendering {

// Virtual base class for providing environment lookups
template <class ColorType>
class environment_map_provider {
  public:
    typedef ColorType color_type;
    typedef frantic::graphics::alpha3f alpha_type;
    typedef frantic::graphics::color_with_alpha<color_type, alpha_type> pixel_type;

    virtual ~environment_map_provider() {}

    // Simple environment lookup with no filter width.
    virtual ColorType lookup_environment( const frantic::graphics::vector3f& direction ) const = 0;

    virtual void render_background( const frantic::graphics::camera<float>& cam,
                                    frantic::graphics2d::framebuffer<pixel_type>& outImage, float mblurTime = 0.5f );

    virtual void render_cube_background( const frantic::graphics::transform4f& toWorld,
                                         frantic::rendering::framebuffer_cubeface<pixel_type>& outImage );
};

template <class ColorType>
void environment_map_provider<ColorType>::render_background( const frantic::graphics::camera<float>& cam,
                                                             frantic::graphics2d::framebuffer<pixel_type>& outImage,
                                                             float mblurTime ) {
    frantic::graphics::transform4f toWorld = cam.world_transform( mblurTime );

    // TODO: This doesn't do any sort of antialiasing. It probably should. Maybe. Ya know. Also this should definitely
    // be done in parallel.
    for( int y = 0, yEnd = outImage.height(); y < yEnd; ++y ) {
        for( int x = 0, xEnd = outImage.width(); x < xEnd; ++x ) {
            bool isValid = true;
            frantic::graphics2d::vector2f p( (float)x + 0.5f, (float)y + 0.5f );
            frantic::graphics::vector3f dir =
                toWorld.transform_no_translation( cam.to_cameraspace_direction( p, isValid ) );

            ColorType c = this->lookup_environment( dir );

            outImage.blend_over( x, y, pixel_type( c, alpha_type( 0.f ) ) );
        }
    }
}

template <class ColorType>
void environment_map_provider<ColorType>::render_cube_background(
    const frantic::graphics::transform4f& toWorld, frantic::rendering::framebuffer_cubeface<pixel_type>& outImage ) {
    frantic::graphics2d::size2 cubeSize( outImage.size() );

    // TODO: This doesn't do any sort of antialiasing. It probably should. Maybe. Ya know. Also this should definitely
    // be done in parallel.
    for( int z = 0; z < 6; ++z ) {
        frantic::graphics::cube_face::default_cube_face cubeFace =
            static_cast<frantic::graphics::cube_face::default_cube_face>( z );

        for( int y = 0, yEnd = cubeSize.xsize; y < yEnd; ++y ) {
            for( int x = 0, xEnd = cubeSize.ysize; x < xEnd; ++x ) {
                frantic::graphics2d::vector2f p( (float)x + 0.5f, (float)y + 0.5f );

                // Convert to [0,1] coordinates.
                p.x /= static_cast<float>( outImage.size() );
                p.y /= static_cast<float>( outImage.size() );

                frantic::graphics::vector3f dir =
                    toWorld.transform_no_translation( frantic::graphics::from_cube_face_coordinate( p, cubeFace ) );

                ColorType c = this->lookup_environment( dir );

                outImage.blend_over( z, x, y, pixel_type( c, alpha_type( 0.f ) ) );
            }
        }
    }
}

} // namespace rendering
} // namespace frantic
