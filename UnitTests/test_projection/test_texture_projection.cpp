// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#ifdef OIIO_LIB_AVAILABLE

#include <gtest/gtest.h>

#include <frantic/projection/texture_projection.hpp>

#include <frantic/graphics/camera.hpp>
#include <frantic/graphics/color_rgba_f.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <frantic/graphics/spherical_coords.hpp>

#include <frantic/graphics/cubeface.hpp>

#include <OpenImageIO/imagebuf.h>
#include <OpenImageIO/typedesc.h>

#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>

#include <cmath>

TEST( TextureProjection, SingleCamera ) {
    using namespace frantic::projection;
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;
    using namespace OpenImageIO;

    const size_t imageHeight = 512;
    const size_t imageWidth = 512;

    projection_mode::projection_mode testModes[] = {
        projection_mode::perspective,
        projection_mode::spherical,
        projection_mode::panoramic,
    };

    for( size_t projType = 0; projType < 3; ++projType ) {

        boost::shared_ptr<ImageBuf> projectionImage(
            new ImageBuf( ImageSpec( imageWidth, imageHeight, 4, TypeDesc::FLOAT ) ) );

        for( int y = 0; y < imageHeight; ++y ) {
            for( int x = 0; x < imageWidth; ++x ) {
                const color_rgba_f pixelColor( float( x ), float( y ), 0.0f, 1.0f );
                projectionImage->setpixel( x, y, (float*)&pixelColor );
            }
        }

        frantic::graphics::camera<double> projectionCamera;
        projectionCamera.set_projection_mode( testModes[projType] );
        projectionCamera.set_transform( frantic::graphics::transform4fd::from_angle_axis(
                                            1.4, frantic::graphics::vector3fd( 3.0, 4.0, 5.0 ).to_normalized() ) *
                                        transform4f::from_translation( -65.0, 12.0, -4.0 ) );
        projectionCamera.set_output_size( frantic::graphics2d::size2( imageWidth, imageHeight ) );

        texture_projection proj;
        proj.set_cubic( false );
        proj.set_projection_camera( projectionCamera );

        proj.set_projection_texture( projectionImage );

        for( int y = 0; y < imageHeight; ++y ) {
            for( int x = 0; x < imageHeight; ++x ) {
                const color_rgba_f expectedColor( float( x ), float( y ), 0.0f, 1.0f );
                bool outValid;
                const vector3fd outDirection = projectionCamera.to_worldspace_direction(
                    vector2f( float( x ) + 0.5f, float( y ) + 0.5f ), outValid );
                const color_rgba_f actualColor =
                    proj.get_projection_color( projectionCamera.camera_position() + outDirection );

                EXPECT_NEAR( expectedColor.get_r(), actualColor.get_r(), 0.1f );
                EXPECT_NEAR( expectedColor.get_g(), actualColor.get_g(), 0.1f );
                EXPECT_NEAR( expectedColor.get_b(), actualColor.get_b(), 0.1f );
                EXPECT_NEAR( 1.0, actualColor.get_a(), 0.0001f );
            }
        }
    }
}

TEST( TextureProjection, Cubic ) {
    using namespace frantic::projection;
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;
    using namespace OpenImageIO;

    boost::shared_ptr<ImageBuf> projectionImages[6];

    const size_t imageHeight = 1024;
    const size_t imageWidth = 1024;

    frantic::graphics::camera<double> projectionCamera;
    projectionCamera.set_transform(
        transform4fd::from_angle_axis( 0.9, frantic::graphics::vector3fd( 5.0, 4.0, 3.0 ).to_normalized() ) *
        transform4f::from_translation( 9.0, -15.0, -3.0 ) );
    projectionCamera.set_output_size( size2( imageWidth, imageHeight ) );

    texture_projection proj;
    proj.set_cubic( true );
    proj.set_projection_camera( projectionCamera );

    for( size_t i = 0; i < 6; ++i ) {
        const cube_face::default_cube_face currentFace = cube_face::default_cube_face( i );

        projectionImages[i].reset( new ImageBuf( ImageSpec( imageWidth, imageHeight, 4, TypeDesc::FLOAT ) ) );

        for( int y = 0; y < imageHeight; ++y ) {
            for( int x = 0; x < imageWidth; ++x ) {
                const vector2f faceCoord( float( x ) / float( imageWidth ), float( y ) / float( imageHeight ) );
                const vector3f position = from_cube_face_coordinate( faceCoord, currentFace ).to_normalized();
                const color_rgba_f pixelColor( position[0], position[1], position[2], 1.0f );
                projectionImages[i]->setpixel( x, y, (float*)&pixelColor );
            }
        }

        proj.set_projection_texture( projectionImages[i], currentFace );
    }

    for( size_t theta = 0; theta < 1000; ++theta ) {
        for( size_t phi = 0; phi < 1000; ++phi ) {
            const spherical_coords sphCoord( double( theta ) / ( 2.0 * M_PI ), double( phi ) / M_PI );
            const vector3fd cameraDirection = sphCoord.to_vector3t();

            const vector3fd worldPosition = projectionCamera.world_transform() * cameraDirection;

            const color_rgba_f actualColor = proj.get_projection_color( worldPosition );

            EXPECT_NEAR( cameraDirection[0], actualColor.get_r(), 0.1f );
            EXPECT_NEAR( cameraDirection[1], actualColor.get_g(), 0.1f );
            EXPECT_NEAR( cameraDirection[2], actualColor.get_b(), 0.1f );
            EXPECT_NEAR( 1.0, actualColor.get_a(), 0.000001f );
        }
    }
}

#endif // OIIO_LIB_AVAILABLE
