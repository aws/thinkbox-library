// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest-helper.h"

#include <frantic/graphics/camera.hpp>
#include <frantic/graphics/vector3f.hpp>

using namespace frantic::graphics;
using namespace frantic::graphics2d;

// from_perspective_projection is actually a property of vector3f but based on the comments there,
// it belongs better in camera instead
TEST( CameraFd, PerspectiveDirection ) {
    // Test center ray direction
    const float fov = float( M_PI / 2 );
    const float tanHalfFov = float( std::tan( fov / 2 ) );
    const size2f size( 128, 256 );
    const vector2f checkPoint( size.get_xsize() / 2, size.get_ysize() / 2 );
    const double aspect = 1;

    vector3fd resultDirection = vector3fd::from_perspective_projection( checkPoint, size, tanHalfFov, aspect );
    vector3fd expectedDirection( 0, 0, -1 );
    EXPECT_VECTOR3FD_EQ( expectedDirection, resultDirection );
}

namespace {

template <typename FloatType>
void test_points_match_up( const camera<FloatType>& camera, const vector3t<FloatType>& testPointWorld,
                           const vector2f& testPointScreen ) {
    const int width = camera.get_output_size().get_xsize();
    const int height = camera.get_output_size().get_ysize();

    bool validOutput = true;
    const ray3t<FloatType> ray = camera.get_cameraspace_ray( testPointScreen, validOutput );

    EXPECT_TRUE( validOutput );
    EXPECT_TRUE( ray.contains_point( testPointWorld ) );

    const transform4t<FloatType> trans = camera.get_projection_transform();
    vector3t<FloatType> testPointWorldOnScreen = trans * testPointWorld;
    testPointWorldOnScreen.x = ( testPointWorldOnScreen.x + 1 ) / 2 * width;
    testPointWorldOnScreen.y = ( testPointWorldOnScreen.y + 1 ) / 2 * height;

    EXPECT_FLOAT_EQ( (float)testPointScreen.x, (float)testPointWorldOnScreen.x );
    EXPECT_FLOAT_EQ( (float)testPointScreen.y, (float)testPointWorldOnScreen.y );
}

template <typename FloatType>
void test_points_match_up( const camera<FloatType>& camera, const vector2f& testPointScreen, FloatType rayTime ) {
    bool validOutput = true;
    const ray3t<FloatType> ray = camera.get_cameraspace_ray( testPointScreen, validOutput );

    test_points_match_up( camera, ray.at( rayTime ), testPointScreen );
}

} // anonymous namespace

TEST( CameraFd, ProjectionTransform ) {
    // Test if camera ray direction and camera to screen transforms are consistent with each other

    camera<double> camera;
    camera.set_near( 1 );
    camera.set_far( 16 );
    camera.set_horizontal_fov( M_PI / 2 );
    camera.set_orthographic_width( 2 );

    // Test is only valid for orthographic and perspective projections
    std::vector<projection_mode::projection_mode> modes;
    modes.push_back( projection_mode::perspective );
    modes.push_back( projection_mode::orthographic );

    // Test some various image sizes
    std::vector<size2> imageSizes;
    imageSizes.push_back( size2( 640, 480 ) );
    imageSizes.push_back( size2( 480, 640 ) );
    imageSizes.push_back( size2( 500, 500 ) );

    // Test some various pixel aspects
    std::vector<double> pixelAspects;
    pixelAspects.push_back( 1 );
    pixelAspects.push_back( 2 );
    pixelAspects.push_back( 0.5 );

    for( std::vector<projection_mode::projection_mode>::const_iterator projectionIter = modes.begin();
         projectionIter != modes.end(); ++projectionIter ) {
        camera.set_projection_mode( *projectionIter );

        for( std::vector<size2>::const_iterator sizeIter = imageSizes.begin(); sizeIter != imageSizes.end();
             ++sizeIter ) {
            const size2& size = *sizeIter;
            camera.set_output_size( size );

            for( std::vector<double>::const_iterator aspectIter = pixelAspects.begin();
                 aspectIter != pixelAspects.end(); ++aspectIter ) {
                camera.set_pixel_aspect( *aspectIter );

                test_points_match_up( camera, vector3fd( 0, 0, 8 ), vector2f( size.xsize / 2.f, size.ysize / 2.f ) );
                test_points_match_up( camera, vector2f( 32, 32 ), 16.0 );
                test_points_match_up( camera, vector2f( 582, 389 ), 16.0 );
                test_points_match_up( camera, vector2f( 43, 132 ), 16.0 );
                test_points_match_up( camera, vector2f( 476, 159 ), 16.0 );
            }
        }
    }
}

namespace {
void checkFOVValid( double expectedFOV, const camera<double>& camera ) {
    EXPECT_DOUBLE_EQ( expectedFOV, camera.fov() );
    const double pixelAspect = camera.pixel_aspect();

    const size2& size = camera.output_size();
    if( size.xsize > size.ysize * pixelAspect ) {
        // Wide image
        EXPECT_DOUBLE_EQ( expectedFOV, camera.vertical_fov() );
        EXPECT_LT( expectedFOV, camera.horizontal_fov() );

    } else if( size.xsize < size.ysize * pixelAspect ) {
        // Tall image
        EXPECT_DOUBLE_EQ( expectedFOV, camera.horizontal_fov() );
        EXPECT_LT( expectedFOV, camera.vertical_fov() );

    } else {
        // Square image
        EXPECT_DOUBLE_EQ( expectedFOV, camera.horizontal_fov() );
        EXPECT_DOUBLE_EQ( expectedFOV, camera.vertical_fov() );
    }
}

void checkFOVValid( double expectedFOV, bool isHorizontal, const camera<double>& camera ) {
    const size2& size = camera.output_size();
    const double pixelAspect = camera.pixel_aspect();

    if( size.xsize > size.ysize * pixelAspect ) {
        // Wide image
        if( isHorizontal ) {
            EXPECT_DOUBLE_EQ( expectedFOV, camera.horizontal_fov() );
            EXPECT_GT( expectedFOV, camera.vertical_fov() );
            EXPECT_GT( expectedFOV, camera.fov() );
        } else {
            EXPECT_DOUBLE_EQ( expectedFOV, camera.fov() );
            EXPECT_DOUBLE_EQ( expectedFOV, camera.vertical_fov() );
            EXPECT_LT( expectedFOV, camera.horizontal_fov() );
        }

    } else if( size.xsize < size.ysize * pixelAspect ) {
        // Tall image
        if( isHorizontal ) {
            EXPECT_DOUBLE_EQ( expectedFOV, camera.fov() );
            EXPECT_DOUBLE_EQ( expectedFOV, camera.horizontal_fov() );
            EXPECT_LT( expectedFOV, camera.vertical_fov() );
        } else {
            EXPECT_DOUBLE_EQ( expectedFOV, camera.vertical_fov() );
            EXPECT_GT( expectedFOV, camera.horizontal_fov() );
            EXPECT_GT( expectedFOV, camera.fov() );
        }

    } else {
        // Square image
        EXPECT_DOUBLE_EQ( expectedFOV, camera.fov() );
        EXPECT_DOUBLE_EQ( expectedFOV, camera.horizontal_fov() );
        EXPECT_DOUBLE_EQ( expectedFOV, camera.vertical_fov() );
    }
}

void checkOrthoValid( double expectedOrtho, const camera<double>& camera ) {
    EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_size() );

    const size2& size = camera.output_size();
    const double pixelAspect = camera.pixel_aspect();
    if( size.xsize > size.ysize * pixelAspect ) {
        // Wide image
        EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_height() );
        EXPECT_LT( expectedOrtho, camera.orthographic_width() );

    } else if( size.xsize < size.ysize * pixelAspect ) {
        // Tall image
        EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_width() );
        EXPECT_LT( expectedOrtho, camera.orthographic_height() );

    } else {
        // Square image
        EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_width() );
        EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_height() );
    }
}

void checkOrthoValid( double expectedOrtho, bool isHorizontal, const camera<double>& camera ) {
    const size2& size = camera.output_size();
    const double pixelAspect = camera.pixel_aspect();
    if( size.xsize > size.ysize * pixelAspect ) {
        // Wide image
        if( isHorizontal ) {
            EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_width() );
            EXPECT_GT( expectedOrtho, camera.orthographic_height() );
            EXPECT_GT( expectedOrtho, camera.orthographic_size() );
        } else {
            EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_size() );
            EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_height() );
            EXPECT_LT( expectedOrtho, camera.orthographic_width() );
        }

    } else if( size.xsize < size.ysize * pixelAspect ) {
        // Tall image
        if( isHorizontal ) {
            EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_size() );
            EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_width() );
            EXPECT_LT( expectedOrtho, camera.orthographic_height() );
        } else {
            EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_height() );
            EXPECT_GT( expectedOrtho, camera.orthographic_width() );
            EXPECT_GT( expectedOrtho, camera.orthographic_size() );
        }

    } else {
        // Square image
        EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_size() );
        EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_width() );
        EXPECT_DOUBLE_EQ( expectedOrtho, camera.orthographic_height() );
    }
}
} // anonymous namespace

TEST( CameraFd, FieldOfView ) {
    std::vector<double> anglesAndSizes;
    anglesAndSizes.push_back( 1 );
    anglesAndSizes.push_back( M_PI / 2 );
    anglesAndSizes.push_back( 2 * M_PI / 3 );
    anglesAndSizes.push_back( M_PI / 4 );

    std::vector<size2> sizes;
    sizes.push_back( size2( 64, 64 ) );
    sizes.push_back( size2( 640, 480 ) );
    sizes.push_back( size2( 800, 600 ) );
    sizes.push_back( size2( 600, 800 ) );
    sizes.push_back( size2( 1024, 768 ) );
    sizes.push_back( size2( 300, 300 ) );
    sizes.push_back( size2( 480, 640 ) );
    sizes.push_back( size2( 768, 1024 ) );

    std::vector<double> aspects;
    aspects.push_back( 1 );
    aspects.push_back( 4 );
    aspects.push_back( 0.25 );

    // Make sure changing the image size after setting the field of view / ortho size does not change it
    for( std::vector<double>::const_iterator valueIter = anglesAndSizes.begin(); valueIter != anglesAndSizes.end();
         ++valueIter ) {
        camera<double> camera;
        camera.set_pixel_aspect( 1 );

        const double fovOrSize = *valueIter;
        camera.set_fov( fovOrSize );
        camera.set_orthographic_size( fovOrSize );

        for( std::vector<size2>::const_iterator sizeIter = sizes.begin(); sizeIter != sizes.end(); ++sizeIter ) {
            const size2& size = *sizeIter;

            camera.set_output_size( size );

            checkFOVValid( fovOrSize, camera );
            checkOrthoValid( fovOrSize, camera );
        }
    }

    // Make sure changing the pixel aspect after setting the field of view / ortho size does not change it
    for( std::vector<double>::const_iterator valueIter = anglesAndSizes.begin(); valueIter != anglesAndSizes.end();
         ++valueIter ) {
        camera<double> camera;
        camera.set_output_size( size2( 256, 256 ) );

        const double fovOrSize = *valueIter;
        camera.set_fov( fovOrSize );
        camera.set_orthographic_size( fovOrSize );

        for( std::vector<double>::const_iterator aspectIter = aspects.begin(); aspectIter != aspects.end();
             ++aspectIter ) {
            camera.set_pixel_aspect( *aspectIter );

            checkFOVValid( fovOrSize, camera );
            checkOrthoValid( fovOrSize, camera );
        }
    }

    // Test if setting / getting the field of view and ortho size of the minimum image dimension works as expected
    for( std::vector<size2>::const_iterator sizeIter = sizes.begin(); sizeIter != sizes.end(); ++sizeIter ) {
        const size2& size = *sizeIter;
        for( std::vector<double>::const_iterator aspectIter = aspects.begin(); aspectIter != aspects.end();
             ++aspectIter ) {
            camera<double> camera;
            camera.set_pixel_aspect( *aspectIter );
            camera.set_output_size( size );

            for( std::vector<double>::const_iterator valueIter = anglesAndSizes.begin();
                 valueIter != anglesAndSizes.end(); ++valueIter ) {
                const double fovOrSize = *valueIter;

                // Make sure setting the field of view and ortho sizes updates the horizontal and vertical values
                // properly
                camera.set_fov( fovOrSize );
                checkFOVValid( fovOrSize, camera );

                camera.set_orthographic_size( fovOrSize );
                checkOrthoValid( fovOrSize, camera );

                // Make sure setting the horizontal or vertical fovs/ortho sizes sets the actual fovs/ortho sizes
                // properly
                camera.set_horizontal_fov( fovOrSize );
                checkFOVValid( fovOrSize, true, camera );

                camera.set_vertical_fov( fovOrSize );
                checkFOVValid( fovOrSize, false, camera );

                camera.set_orthographic_width( fovOrSize );
                checkOrthoValid( fovOrSize, true, camera );

                camera.set_orthographic_height( fovOrSize );
                checkOrthoValid( fovOrSize, false, camera );
            }
        }
    }
}

TEST( CameraFd, SetPositionOrientation ) {
    camera<double> camera;
    const vector3fd translation( -30., 1.3, 4.4446 );
    camera.set_position( translation );
    EXPECT_VECTOR3FD_EQ( translation, camera.world_transform().translation() );

    const quat4fd rotation( 0.185, -0.003, -0.932, 0.311 );
    camera.set_orientation( rotation );
    EXPECT_QUAT4FD_NEAR( rotation, -camera.camera_orientation(), 0.0003 );
    EXPECT_VECTOR3FD_EQ( translation, camera.world_transform().translation() );

    const vector3fd translation2( -2., 153.436, -5.925 );
    camera.set_position( translation2 );
    EXPECT_QUAT4FD_NEAR( rotation, -camera.camera_orientation(), 0.0003 );
    EXPECT_VECTOR3FD_EQ( translation2, camera.world_transform().translation() );
}

TEST( CameraFd, RotateCamera ) {
    const vector3fd worldUp( 0., 1., 0. );
    {
        camera<double> camera;
        camera.rotate_camera( 2 * M_PI, 2 * M_PI, worldUp );
        EXPECT_QUAT4FD_NEAR( quat4fd(), camera.camera_orientation(), 0.0003 );
    }
    {
        camera<double> camera;
        camera.rotate_camera( 0., M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -1. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( -1., 0., 0. ), camera.view_direction(), 1e-15 );
    }
    {
        camera<double> camera;
        camera.rotate_camera( 0., -M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_right(), 1e-15 );
    }
    {
        camera<double> camera;
        camera.rotate_camera( M_PI / 2, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
    }
    {
        camera<double> camera;
        camera.rotate_camera( -M_PI / 2, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -1. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
    }
    {
        camera<double> camera;
        camera.rotate_camera( M_PI, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
        // Test upside down rotation reversal
        camera.rotate_camera( 0., M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_right(), 1e-15 );
    }
    {
        camera<double> camera;
        camera.rotate_camera( 0., M_PI, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( -1., 0., 0. ), camera.view_right(), 1e-15 );
    }
    {
        camera<double> camera;
        camera.rotate_camera( M_PI / 2, M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -1. ), camera.view_right(), 1e-15 );
    }
    {
        camera<double> camera;
        const vector3fd position( 20., 40., 75. );
        camera.set_position( position );

        camera.rotate_camera( M_PI, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
        // Test upside down rotation reversal
        camera.rotate_camera( 0., M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_right(), 1e-15 );

        EXPECT_VECTOR3FD_NEAR( position, camera.camera_position(), 1e-15 );
    }
}

TEST( CameraFd, OrbitCameraAroundCentre ) {
    // Almost same as RotateCamera
    const vector3fd worldUp( 0., 1., 0. );
    const vector3fd pivot( 7.777777, 100., -12.3 );
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, 2 * M_PI, 2 * M_PI, worldUp );
        EXPECT_QUAT4FD_NEAR( quat4fd(), camera.camera_orientation(), 0.0003 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, 0., M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -1. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( -1., 0., 0. ), camera.view_direction(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, 0., -M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, M_PI / 2, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, -M_PI / 2, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -1. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, M_PI, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
        // Test upside down rotation reversal
        camera.orbit_camera( pivot, 0., M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, 0., M_PI, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( -1., 0., 0. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
    {
        camera<double> camera;
        camera.set_position( pivot );
        camera.orbit_camera( pivot, M_PI / 2, M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -1. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_EQ( pivot, camera.camera_position() );
    }
}

TEST( CameraFd, OrbitCameraAroundExternal ) {
    const vector3fd worldUp( 0., 1., 0. );
    const vector3fd pivot;
    const vector3fd startingPos( 0., 0., 5. );
    {
        camera<double> camera;
        camera.set_position( startingPos );
        camera.orbit_camera( pivot, 2 * M_PI, 2 * M_PI, worldUp );
        EXPECT_QUAT4FD_NEAR( quat4fd(), camera.camera_orientation(), 0.0003 );
        EXPECT_VECTOR3FD_NEAR( startingPos, camera.camera_position(), 1e-14 );
    }
    {
        camera<double> camera;
        camera.set_position( startingPos );
        camera.orbit_camera( pivot, 0., M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -1. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 5., 0., 0. ), camera.camera_position(), 1e-14 );
    }
    {
        camera<double> camera;
        camera.set_position( startingPos );
        camera.orbit_camera( pivot, M_PI / 2, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -5., 0. ), camera.camera_position(), 1e-14 );
    }
    {
        camera<double> camera;
        camera.set_position( startingPos );
        camera.orbit_camera( pivot, M_PI, 0., worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 1., 0., 0. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -5. ), camera.camera_position(), 1e-14 );
        // Test upside down rotation reversal
        camera.orbit_camera( pivot, 0., M_PI / 2, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., -1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., 1. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 5., 0., 0. ), camera.camera_position(), 1e-14 );
    }
    {
        camera<double> camera;
        camera.set_position( startingPos );
        camera.orbit_camera( pivot, 0., M_PI, worldUp );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 1., 0. ), camera.view_up(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( -1., 0., 0. ), camera.view_right(), 1e-15 );
        EXPECT_VECTOR3FD_NEAR( vector3fd( 0., 0., -5. ), camera.camera_position(), 1e-14 );
    }
}

namespace {
bool isVisible( const camera<double>& camera, const vector3fd& point ) {
    vector3fd v = camera.get_projection_transform() * camera.world_transform_inverse() * point;
    return ( -1. <= v.x ) && ( v.x <= 1. ) && ( -1. <= v.y ) && ( v.y <= 1. ) && ( v.z <= 1. );
}

bool isCompletelyVisible( const camera<double>& camera, const boundbox3fd& box ) {
    for( int i = 0; i < 8; ++i ) {
        if( !isVisible( camera, box.get_corner( i ) ) )
            return false;
    }
    return true;
}
} // namespace

TEST( CameraFd, GetZoomDistancePoint ) {
    {
        camera<double> camera;
        vector3fd point( 0., 0.1, -2. );
        EXPECT_TRUE( isVisible( camera, point ) );
        double distance = camera.get_suggested_depth( point );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isVisible( camera, point ) );
    }
    {
        camera<double> camera;
        vector3fd point( 0., 5., -2. );
        EXPECT_FALSE( isVisible( camera, point ) );
        double distance = camera.get_suggested_depth( point );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isVisible( camera, point ) );
    }
    {
        camera<double> camera;
        vector3fd point( 0., -6., -2.2 );
        EXPECT_FALSE( isVisible( camera, point ) );
        double distance = camera.get_suggested_depth( point );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isVisible( camera, point ) );
    }
    {
        camera<double> camera;
        vector3fd point( 0., -1., 2.2 );
        EXPECT_FALSE( isVisible( camera, point ) );
        double distance = camera.get_suggested_depth( point );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isVisible( camera, point ) );
    }
    {
        camera<double> camera;
        vector3fd point( 0., 0., 3. );
        EXPECT_FALSE( isVisible( camera, point ) );
        double distance = camera.get_suggested_depth( point );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isVisible( camera, point ) );
    }
    {
        camera<double> camera;
        camera.set_fov_degrees( 2.0 );
        vector3fd point( -1., 3., 12. );
        EXPECT_FALSE( isVisible( camera, point ) );
        double distance = camera.get_suggested_depth( point );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isVisible( camera, point ) );
    }
    {
        camera<double> camera;
        camera.set_fov_degrees( 120.0 );
        vector3fd point( -1., 3., 12. );
        EXPECT_FALSE( isVisible( camera, point ) );
        double distance = camera.get_suggested_depth( point );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isVisible( camera, point ) );
    }
    {
        camera<double> camera( projection_mode::orthographic, motion_blurred_transform<double>(), 2.0 );
        vector3fd point( 500., 0., -2. );
        try {
            camera.get_suggested_depth( point );
            GTEST_FAIL();
        } catch( const std::runtime_error& ) {
        }
    }
}

TEST( CameraFd, GetZoomDistanceBBox ) {
    {
        camera<double> camera;
        boundbox3fd bbox( -0.1, 0.1, -0.1, 0.1, -3.0, -2.0 );
        EXPECT_TRUE( isCompletelyVisible( camera, bbox ) );
        double distance = camera.get_suggested_depth( bbox );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isCompletelyVisible( camera, bbox ) );
    }
    {
        camera<double> camera;
        boundbox3fd bbox( -10., 20., -70., 0.1, -3.0, -2.0 );
        EXPECT_FALSE( isCompletelyVisible( camera, bbox ) );
        double distance = camera.get_suggested_depth( bbox );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isCompletelyVisible( camera, bbox ) );
    }
    {
        camera<double> camera;
        boundbox3fd bbox( -0.1, 0.1, -0.1, 0.1, 7.0, 7.5 );
        EXPECT_FALSE( isCompletelyVisible( camera, bbox ) );
        double distance = camera.get_suggested_depth( bbox );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isCompletelyVisible( camera, bbox ) );
    }
    {
        camera<double> camera;
        camera.set_fov_degrees( 2.0 );
        boundbox3fd bbox( -0.1, 0.1, -0.1, 0.1, 7.0, 7.5 );
        EXPECT_FALSE( isCompletelyVisible( camera, bbox ) );
        double distance = camera.get_suggested_depth( bbox );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isCompletelyVisible( camera, bbox ) );
    }
    {
        camera<double> camera;
        camera.set_fov_degrees( 144. );
        boundbox3fd bbox( -0.1, 0.1, -0.1, 0.1, 7.0, 7.5 );
        EXPECT_FALSE( isCompletelyVisible( camera, bbox ) );
        double distance = camera.get_suggested_depth( bbox );
        camera.set_position( camera.camera_position() - distance * camera.view_direction() );
        EXPECT_TRUE( isCompletelyVisible( camera, bbox ) );
    }
}

TEST( CameraFd, CreateCubicCameras ) {
    std::vector<camera<float>> cameraList;

    {
        frantic::graphics::camera<float> camera;
        camera.set_far( 400.0f );
        camera.set_near( 10.0f );
        camera.set_position( vector3f( 5.0f, -15.0f, 33.33f ) );

        cameraList.push_back( camera );
    }

    {
        frantic::graphics::camera<float> camera;
        camera.set_fov( 13.5f );
        camera.set_fstop( 20.0f );
        camera.set_pixel_aspect( 3.0f );
        camera.set_position( vector3f( 1000.0f, 2000.0f, -111.0f ) );

        cameraList.push_back( camera );
    }

    {
        frantic::graphics::camera<float> camera;
        cameraList.push_back( camera );
    }

    {
        frantic::graphics::camera<float> camera;
        camera.set_far( 2.0f );
        camera.set_near( 0.01f );
        camera.set_fstop( 1.39f );
        camera.set_subpixel_offset( vector2f( 19.5f, -3.2f ) );
        camera.set_fov( 100.0f );
        camera.set_focal_distance( 50.0f );
        camera.set_focal_length( 9.992f );
        camera.set_pixel_aspect( 15.0f );

        cameraList.push_back( camera );
    }

    for( std::vector<frantic::graphics::camera<float>>::const_iterator currCamera = cameraList.begin();
         currCamera != cameraList.end(); ++currCamera ) {
        const std::vector<frantic::graphics::camera<float>> cubicCameras = currCamera->get_cubic_cameras( 500 );

        EXPECT_EQ( 6, cubicCameras.size() );

        const vector3f dirXPos( 1.0f, 0.0f, 0.0f );
        const vector3f dirXNeg( -1.0f, 0.0f, 0.0f );
        const vector3f dirYPos( 0.0f, 1.0f, 0.0f );
        const vector3f dirYNeg( 0.0f, -1.0f, 0.0f );
        const vector3f dirZNeg( 0.0f, 0.0f, -1.0f );
        const vector3f dirZPos( 0.0f, 0.0f, 1.0f );

        EXPECT_VECTOR3F_EQ( dirZPos, cubicCameras[0].view_direction() );
        EXPECT_VECTOR3F_EQ( dirZNeg, cubicCameras[1].view_direction() );
        EXPECT_VECTOR3F_EQ( dirYNeg, cubicCameras[2].view_direction() );
        EXPECT_VECTOR3F_EQ( dirYPos, cubicCameras[3].view_direction() );
        EXPECT_VECTOR3F_EQ( dirXNeg, cubicCameras[4].view_direction() );
        EXPECT_VECTOR3F_EQ( dirXPos, cubicCameras[5].view_direction() );

        EXPECT_FLOAT_EQ( 90.0f, cubicCameras[0].fov_degrees() );
        EXPECT_FLOAT_EQ( 90.0f, cubicCameras[1].fov_degrees() );
        EXPECT_FLOAT_EQ( 90.0f, cubicCameras[2].fov_degrees() );
        EXPECT_FLOAT_EQ( 90.0f, cubicCameras[3].fov_degrees() );
        EXPECT_FLOAT_EQ( 90.0f, cubicCameras[4].fov_degrees() );
        EXPECT_FLOAT_EQ( 90.0f, cubicCameras[5].fov_degrees() );

        for( std::vector<frantic::graphics::camera<float>>::const_iterator it = cubicCameras.begin();
             it != cubicCameras.end(); ++it ) {
            EXPECT_EQ( size2( 500, 500 ), it->get_output_size() );
            EXPECT_FLOAT_EQ( 1.0f, it->pixel_aspect() );
            EXPECT_FLOAT_EQ( it->near_distance(), currCamera->near_distance() );
            EXPECT_FLOAT_EQ( it->far_distance(), currCamera->far_distance() );
            EXPECT_FLOAT_EQ( it->fstop(), currCamera->fstop() );
            EXPECT_FLOAT_EQ( it->focal_length(), currCamera->focal_length() );
            EXPECT_FLOAT_EQ( it->focal_distance(), currCamera->focal_distance() );
            EXPECT_EQ( it->get_subpixel_offset(), currCamera->get_subpixel_offset() );
        }

        EXPECT_VECTOR3F_EQ( currCamera->camera_position(), cubicCameras[0].camera_position() );
        EXPECT_VECTOR3F_EQ( currCamera->camera_position(), cubicCameras[1].camera_position() );
        EXPECT_VECTOR3F_EQ( currCamera->camera_position(), cubicCameras[2].camera_position() );
        EXPECT_VECTOR3F_EQ( currCamera->camera_position(), cubicCameras[3].camera_position() );
        EXPECT_VECTOR3F_EQ( currCamera->camera_position(), cubicCameras[4].camera_position() );
        EXPECT_VECTOR3F_EQ( currCamera->camera_position(), cubicCameras[5].camera_position() );
    }
}
