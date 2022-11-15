// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/camera.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/graphics/size3.hpp>
#include <frantic/graphics/spherical_coords.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/volumetrics/voxel_edge_stepper.hpp>

using namespace std;
using namespace boost;
using namespace frantic::graphics;
using namespace frantic::volumetrics;
using namespace frantic::graphics2d;
using namespace frantic::math;

TEST( GraphicsTest, Vector3Cubeface ) {

    ASSERT_EQ( vector3( 1, 0, 0 ).neighbor_vector_to_cube_face(), 0 );
    ASSERT_EQ( vector3( -1, 0, 0 ).neighbor_vector_to_cube_face(), 1 );
    ASSERT_EQ( vector3( 0, 1, 0 ).neighbor_vector_to_cube_face(), 2 );
    ASSERT_EQ( vector3( 0, -1, 0 ).neighbor_vector_to_cube_face(), 3 );
    ASSERT_EQ( vector3( 0, 0, 1 ).neighbor_vector_to_cube_face(), 4 );
    ASSERT_EQ( vector3( 0, 0, -1 ).neighbor_vector_to_cube_face(), 5 );

    ASSERT_EQ( vector3( 0, 0, 0 ).neighbor_vector_to_cube_face_nothrow(), -1 );
    ASSERT_EQ( vector3( 0, 1, 1 ).neighbor_vector_to_cube_face_nothrow(), -1 );
    ASSERT_EQ( vector3( 0, -2, 0 ).neighbor_vector_to_cube_face_nothrow(), -1 );
    ASSERT_EQ( vector3( -1, 1, 0 ).neighbor_vector_to_cube_face_nothrow(), -1 );
}

TEST( GraphicsTest, Vector3Index ) {
    size3 size( 10, 11, 7 );
    int index = 0;
    for( int z = 0; z < size.zsize(); ++z ) {
        for( int y = 0; y < size.ysize(); ++y ) {
            for( int x = 0; x < size.xsize(); ++x, ++index ) {
                vector3 v( x, y, z );
                ASSERT_EQ( size.get_index( v ), index );
                ASSERT_TRUE( size.get_coord( index ) == v );
            }
        }
    }
}

TEST( GraphicsTest, BoundBox3Basics ) {

    ////// Test construction //////
    boundbox3 boxFromDefaultConstructor;
    ASSERT_TRUE( boxFromDefaultConstructor.is_empty() );

    boundbox3 boxFromOneVector( vector3( 3, 2, 7 ) );
    ASSERT_EQ( boxFromOneVector.minimum(), vector3( 3, 2, 7 ) );
    ASSERT_EQ( boxFromOneVector.maximum(), vector3( 3, 2, 7 ) );

    boundbox3 boxFromTwoVectors( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    ASSERT_EQ( boxFromTwoVectors.minimum(), vector3( 1, 7, 2 ) );
    ASSERT_EQ( boxFromTwoVectors.maximum(), vector3( 2, 9, 5 ) );

    boundbox3 boxFromVectorAndSize( vector3( -3, 2, 2 ), size3( 7, 1, 9 ) );
    ASSERT_EQ( boxFromVectorAndSize.minimum(), vector3( -3, 2, 2 ) );
    ASSERT_EQ( boxFromVectorAndSize.maximum(), vector3( 3, 2, 10 ) );

    int intArray[] = { 3, 7, 9, 12, 9, 18 };
    boundbox3 boxFromArray( intArray );
    ASSERT_EQ( boxFromArray.minimum(), vector3( 3, 9, 9 ) );
    ASSERT_EQ( boxFromArray.maximum(), vector3( 7, 12, 18 ) );

    ////// Test set functions //////
    boundbox3 box;
    box.set( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    ASSERT_EQ( box.minimum(), vector3( 1, 7, 2 ) );
    ASSERT_EQ( box.maximum(), vector3( 2, 9, 5 ) );

    box.set_to_empty();
    ASSERT_TRUE( box.is_empty() );

    box.set_to_point( vector3( 3, 2, 7 ) );
    ASSERT_EQ( box.minimum(), vector3( 3, 2, 7 ) );
    ASSERT_EQ( box.maximum(), vector3( 3, 2, 7 ) );

    ASSERT_TRUE( boundbox3::from_empty().is_empty() );

    ////// Test various operations //////
    box.set( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    ASSERT_EQ( box.xsize(), 2 );
    ASSERT_EQ( box.ysize(), 3 );
    ASSERT_EQ( box.zsize(), 4 );
    ASSERT_EQ( box.size(), size3( 2, 3, 4 ) );
    ASSERT_TRUE( !box.is_empty() );
    ASSERT_TRUE( !box.is_cube() );
    ASSERT_EQ( box.get_max_dimension(), 4 );
    ASSERT_EQ( box.get_max_dimension_axis(), 2 );
    ASSERT_EQ( box.get_min_dimension(), 2 );
    ASSERT_LT( vector3f::distance( box.center(), vector3f( 2, 8.5, 4 ) ), 0.0000001f );
    ASSERT_EQ( box.get_volume(), 24 );

    box.set( vector3( 1, 9, 4 ), vector3( 4, 12, 7 ) );
    ASSERT_TRUE( box.is_cube() );

    // Validate the 8 corners.  They should be in 'binary' order.
    box.set( vector3( 0 ), vector3( 1 ) );
    for( int i = 0; i < 8; ++i ) {
        vector3 corner = box.get_corner( i );
        ASSERT_EQ( i, corner.x * 4 + corner.y * 2 + corner.z );
    }

    ////// Test contains //////
    box.set( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    ASSERT_TRUE( box.contains( vector3( 2, 9, 5 ) ) );
    ASSERT_TRUE( box.contains( vector3( 2, 7, 5 ) ) );
    ASSERT_TRUE( box.contains( vector3( 1, 8, 3 ) ) );
    ASSERT_TRUE( !box.contains( vector3( 0, 8, 3 ) ) );

    // For floating point containment, the box is the half-open set defined by the union of all its voxels.
    ASSERT_TRUE( box.contains( vector3f( 1, 7, 2 ) ) );
    ASSERT_TRUE( box.contains( vector3f( 2.9f, 9.9f, 5.9f ) ) );
    ASSERT_TRUE( !box.contains( vector3f( 3, 10, 11 ) ) );

    boundbox3 otherBox( vector3( 2, 7, 3 ), vector3( 2, 8, 5 ) );
    ASSERT_TRUE( box.contains( otherBox ) );
    ASSERT_TRUE( !otherBox.contains( box ) );
    otherBox.set( vector3( 2, 7, 3 ), vector3( 2, 10, 5 ) );
    ASSERT_TRUE( !box.contains( otherBox ) );
    ASSERT_TRUE( !otherBox.contains( box ) );

    ////// Test clamp //////
    box.set( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    ASSERT_EQ( box.clamp( vector3( -4, 8, 9 ) ), vector3( 1, 8, 5 ) );

    ////// Test expand //////
    box.set( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    box.expand( 3 );
    ASSERT_EQ( box, boundbox3( vector3( -2, 4, -1 ), vector3( 5, 12, 8 ) ) );
    box.expand( vector3( 1, 2, 3 ) );
    ASSERT_EQ( box, boundbox3( vector3( -3, 2, -4 ), vector3( 6, 14, 11 ) ) );

    ////// Test parsing //////
    box.set( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    ASSERT_EQ( boundbox3::parse( "boundbox3([1,7,2],[2,9,5])" ), box );
    ASSERT_EQ( boundbox3::parse( box.str() ), box );

    ////// Test union and intersection //////
    box.set( vector3( 1, 7, 2 ), vector3( 2, 9, 5 ) );
    box += vector3( 5, 8, 1 );
    ASSERT_EQ( box, boundbox3( vector3( 1, 7, 1 ), vector3( 5, 9, 5 ) ) );
    box += boundbox3( vector3( 2, 3, 1 ), vector3( 4, 12, 6 ) );
    ASSERT_EQ( box, boundbox3( vector3( 1, 3, 1 ), vector3( 5, 12, 6 ) ) );
    box.intersect_with( boundbox3( vector3( 2, 6, 4 ), vector3( 4, 19, 8 ) ) );
    ASSERT_EQ( box, boundbox3( vector3( 2, 6, 4 ), vector3( 4, 12, 6 ) ) );
}

TEST( GraphicsTest, VoxelEdgeStepper ) {

    boundbox3f startingBox = boundbox3f::from_center_and_width( vector3f(), 1000 );

    vector<pair<vector3f, vector3f>> stepperTests;
    /////// generate test cases
    // Some preset tests
    stepperTests.push_back(
        make_pair( vector3f( -265.77f, -494.14f, 334.864f ), vector3f( -0.892006f, 8.552e-005f, 0.452024f ) ) );
    stepperTests.push_back(
        make_pair( vector3f( 189.566f, 167.928f, 81.7438f ), vector3f( -3.58724e-005f, 0.748325f, -0.663333f ) ) );
    stepperTests.push_back( make_pair( vector3f( 100 ), vector3f( -0.545058f, 0.838398f, -1.44654e-005f ) ) );
    stepperTests.push_back( make_pair( vector3f(), vector3f( 1, 0.5f, 0.25f ).to_normalized() ) );
    stepperTests.push_back( make_pair( vector3f(), vector3f( .8f, .123f, 0.01423f ).to_normalized() ) );
    stepperTests.push_back( make_pair( vector3f(), vector3f( 1, 1, 1 ).to_normalized() ) );
    // A bunch of random tests
    for( int i = 0; i < 100; ++i )
        stepperTests.push_back( make_pair( startingBox.random_vector(), vector3f::from_unit_random() ) );

    for( unsigned stepperTestIndex = 0; stepperTestIndex < stepperTests.size(); ++stepperTestIndex ) {
        ray3f ray( stepperTests[stepperTestIndex].first, stepperTests[stepperTestIndex].second );

        voxel_edge_stepper ves( ray );

        for( int i = 0; i < 10000; ++i ) {
            ves.step();
            // NOTE: both of these checks are -very- lenient, because the stepper can produce
            // very inaccurate results when it stabilizes the cases of the lines and corners of the grid.
            if( !ray.contains_point( ves.isect_position(), 0.001f ) ) {
                string message = "voxel_edge_stepper_regressions: The stepper along ray " + ray.str() +
                                 " stepped away from the ray, to point " + ves.isect_position().str();
                cout << message << endl;
                cout << "point: " << ves.isect_position() << endl;
                cout << "ray: " << ray << endl;
                cout << "point-ray distance: " << ray.distance_from_point( ves.isect_position() ) << endl;
                cout << "step: " << i << endl;
                ves.debug_print( cout );
                FAIL() << message;
            }
            // TODO: Could likely better characterize this as accumulated error per step.
            if( abs( vector3f::distance( ray.origin(), ves.isect_position() ) - ves.isect_distance() ) /
                    ves.isect_position().max_abs_component() >
                0.005f ) {
                string message =
                    "voxel_edge_stepper_regressions: The stepper along ray " + ray.str() +
                    " reported an invalid distance. reported: " + lexical_cast<string>( ves.isect_distance() ) +
                    ", actual " + lexical_cast<string>( vector3f::distance( ray.origin(), ves.isect_position() ) );
                cout << message << endl;
                cout << "step: " << i << endl;
                ves.debug_print( cout );
                FAIL() << message;
            }
            //			out << "P: " << ves.isect_position() << ", axis: " << ves.isect_axis() << ", dist: " <<
            // ves.isect_distance() << ", next voxel: " << ves.isect_voxel_being_entered() << endl;
        }
    }
}

TEST( GraphicsTest, SphericalCoords ) {

    frantic::graphics2d::size2 s( 640, 480 );
    frantic::graphics2d::vector2 p;
    bool isValid;
    for( p.y = 0; p.y < s.ysize; ++p.y ) {
        for( p.x = 0; p.x < s.xsize; ++p.x ) {
            frantic::graphics::spherical_coords sc = frantic::graphics::spherical_coords(
                vector3f::from_perspective_projection( p, s, tan( 60 * float( M_PI ) / 180 / 2 ) ) );
            isValid = true;
            frantic::graphics2d::vector2 pEquiv =
                sc.to_vector3f().to_perspective_projection( s, tan( 60 * float( M_PI ) / 180 / 2 ), isValid );
            ASSERT_TRUE( isValid );
            if( p != pEquiv )
                cout << "p: " << p << ", sc: " << sc << ", pEquiv: " << pEquiv << endl;
            ASSERT_EQ( p, pEquiv );
        }
    }
    frantic::graphics2d::vector2 a, b;
    frantic::graphics2d::size2f sf( s );
    isValid = true;
    a = frantic::graphics2d::vector2( 0, 240 );
    b = frantic::graphics2d::vector2::round(
        frantic::graphics::spherical_coords( float( M_PI ), 150 * float( M_PI ) / 180 )
            .to_vector3f()
            .to_perspective_projection( sf, (float)tan( 60 * M_PI / 180 / 2 ), isValid ) );
    ASSERT_TRUE( isValid );
    if( a != b )
        cout << "a: " << a << ", b: " << b << endl;
    ASSERT_EQ( a, b );
    isValid = true;
    a = frantic::graphics2d::vector2( 640, 240 );
    b = frantic::graphics2d::vector2::round(
        frantic::graphics::spherical_coords( 0, 150 * float( M_PI ) / 180 )
            .to_vector3f()
            .to_perspective_projection( sf, (float)tan( 60 * M_PI / 180 / 2 ), isValid ) );
    ASSERT_TRUE( isValid );
    if( a != b )
        cout << "a: " << a << ", b: " << b << endl;
    ASSERT_EQ( a, b );
    isValid = true;
    a = frantic::graphics2d::vector2( 320, -80 );
    b = frantic::graphics2d::vector2::round(
        frantic::graphics::spherical_coords( 3 * float( M_PI ) / 2, 150 * float( M_PI ) / 180 )
            .to_vector3f()
            .to_perspective_projection( sf, (float)tan( 60 * M_PI / 180 / 2 ), isValid ) );
    ASSERT_TRUE( isValid );
    if( a != b )
        cout << "a: " << a << ", b: " << b << endl;
    ASSERT_EQ( a, b );
    isValid = true;
    a = frantic::graphics2d::vector2( 320, 560 );
    b = frantic::graphics2d::vector2::round(
        frantic::graphics::spherical_coords( float( M_PI ) / 2, 150 * float( M_PI ) / 180 )
            .to_vector3f()
            .to_perspective_projection( sf, (float)tan( 60 * M_PI / 180 / 2 ), isValid ) );
    ASSERT_TRUE( isValid );
    if( a != b )
        cout << "a: " << a << ", b: " << b << endl;
    ASSERT_EQ( a, b );

    // Test latitude-longitude usage
    frantic::graphics2d::vector2f c = frantic::graphics::spherical_coords::from_longlat( 0.0, 0.0 ).to_longlat();
    ASSERT_TRUE( fabs( c.x - 0 ) < 0.0001f );
    ASSERT_TRUE( fabs( c.y - 0 ) < 0.0001f );

    c = frantic::graphics::spherical_coords::from_longlat( 0.5, 0.5 ).to_longlat();
    ASSERT_TRUE( fabs( c.x - 0.5f ) < 0.0001f );
    ASSERT_TRUE( fabs( c.y - 0.5f ) < 0.0001f );
}

TEST( GraphicsTest, RawByteBuffer ) {

    raw_byte_buffer rbb;
    // On creation, the buffer should consist of a single NULL pointer, with no data or capacity
    ASSERT_TRUE( rbb.empty() );
    ASSERT_EQ( rbb.size(), 0 );
    ASSERT_EQ( rbb.capacity(), 0 );
    ASSERT_EQ( rbb.begin(), (char*)0 );
    ASSERT_EQ( rbb.end(), (char*)0 );
    ASSERT_LT( 2000000000, rbb.max_size() );

    // After a reserve, it should still be empty, but now should be pointing inside an allocated buffer
    rbb.reserve( 100 );
    ASSERT_TRUE( rbb.empty() );
    ASSERT_EQ( rbb.size(), 0 );
    ASSERT_EQ( rbb.capacity(), 100 );
    ASSERT_NE( rbb.begin(), (char*)0 );
    ASSERT_NE( rbb.end(), (char*)0 );

    // Adding the first element (without triggering a reallocation) should point to the begin pointer
    ASSERT_EQ( rbb.add_element( 10 ), rbb.begin() );
    ASSERT_NE( rbb.add_element( 10 ), rbb.begin() );

    // A resize should change the size of the buffer
    rbb.resize( 0 );
    ASSERT_TRUE( rbb.empty() );
    rbb.resize( 50 );
    ASSERT_EQ( rbb.size(), 50 );
    rbb.resize( 500 );
    ASSERT_EQ( rbb.size(), 500 );
    ASSERT_LE( rbb.size(), rbb.capacity() );

    // Clearing the buffer should free the memory and return it back to a single NULL pointer
    rbb.clear();
    ASSERT_TRUE( rbb.empty() );
    ASSERT_EQ( rbb.size(), 0 );
    ASSERT_EQ( rbb.capacity(), 0 );
    ASSERT_EQ( rbb.begin(), (char*)0 );
    ASSERT_EQ( rbb.end(), (char*)0 );

    // Adding an element to a buffer with a NULL pointer should allocate some extra room to start
    rbb.add_element( 10 );
    ASSERT_LT( rbb.size(), rbb.capacity() );

    // Fill the values so we're not comparing uninitialized memory in the following test
    memcpy( rbb.begin(), "abcdefghij", 10 );

    raw_byte_buffer rbb2;
    // They shouldn't match, because they're different sizes
    ASSERT_TRUE( rbb != rbb2 );
    ASSERT_TRUE( !( rbb == rbb2 ) );
    // The should match after an assignment though
    rbb2 = rbb;
    ASSERT_TRUE( rbb == rbb2 );
    ASSERT_TRUE( !( rbb != rbb2 ) );
    // If we force a reallocation using reserve, they should still match
    rbb.reserve( 1000 );
    ASSERT_TRUE( rbb == rbb2 );
    ASSERT_TRUE( !( rbb != rbb2 ) );
    // They shouldn't match if the values differ
    ( *rbb2.begin() )++;
    ASSERT_TRUE( rbb != rbb2 );
    ASSERT_TRUE( !( rbb == rbb2 ) );
    // A buffer should match after copy construction
    raw_byte_buffer rbb3( rbb );
    ASSERT_TRUE( rbb == rbb3 );

    // A buffer with a NULL pointer and a zero-sized  buffer with an allocation should be equal
    rbb.clear();
    rbb2.clear();
    rbb2.reserve( 50 );
    ASSERT_TRUE( rbb == rbb2 );
    ASSERT_TRUE( rbb2 == rbb );
    ASSERT_TRUE( !( rbb != rbb2 ) );
    ASSERT_TRUE( !( rbb2 != rbb ) );

    // Adding an element with provided data should work
    char testData[] = { 1, 2, 3 };
    rbb.add_element( testData, sizeof( testData ) );
    ASSERT_EQ( sizeof( testData ), rbb.size() );
    ASSERT_EQ( memcmp( testData, rbb.begin(), sizeof( testData ) ), 0 );
    // Adding an element by getting the pointer, then copying in the data should work
    memcpy( rbb.add_element( sizeof( testData ) ), testData, sizeof( testData ) );
    ASSERT_EQ( 2 * sizeof( testData ), rbb.size() );
    ASSERT_EQ( memcmp( testData, rbb.begin() + sizeof( testData ), sizeof( testData ) ), 0 );

    { // scope for new raw_byte_buffer
        // test resize_with_exponential_growth

        raw_byte_buffer rbb;
        rbb.resize_with_exponential_growth( 0 );
        ASSERT_EQ( rbb.size(), 0 );

        rbb.resize_with_exponential_growth( 1 );
        ASSERT_EQ( rbb.size(), 1 );
        ASSERT_TRUE( rbb.capacity() >= 1 );

        rbb.resize_with_exponential_growth( 100 );
        ASSERT_EQ( rbb.size(), 100 );
        // At the time of implementation, this would pass with
        // capacity() >= 200 .  I'm comparing against a smaller value in
        // case you decrease the exponent.
        ASSERT_TRUE( rbb.capacity() >= 119 );

        // This should cause an exponential increase in capacity
        rbb.resize_with_exponential_growth( 101 );
        ASSERT_EQ( rbb.size(), 101 );
        ASSERT_TRUE( rbb.capacity() >= 119 );

        // This should keep the same size
        rbb.resize_with_exponential_growth( 101 );
        ASSERT_EQ( rbb.size(), 101 );
        ASSERT_TRUE( rbb.capacity() >= 101 );

        // This should reduce the size to 100
        rbb.resize_with_exponential_growth( 100 );
        ASSERT_EQ( rbb.size(), 100 );
        ASSERT_TRUE( rbb.capacity() >= 100 );

        // This should reduce the size to 0
        rbb.resize_with_exponential_growth( 0 );
        ASSERT_EQ( rbb.size(), 0 );
    }

    { // Create from zero bytes of data
        const int data = 0;

        raw_byte_buffer buffer( &data, 0 );

        EXPECT_EQ( 0, buffer.size() );
    }

    { // Create from one byte of data
        const char data = 'a';

        raw_byte_buffer buffer( &data, sizeof( data ) );

        ASSERT_EQ( 1, buffer.size() );
        EXPECT_EQ( data, buffer.begin()[0] );
    }
}

TEST( GraphicsTest, PanoramicCameraAxes ) {

    size2f size( 1024, 512 );
    float hfov = float( 2.0 * M_PI );

    vector3f axes[6] = { vector3f( 0, 0, -1 ), vector3f( 0, 0, 1 ),  vector3f( 0, -1, 0 ),
                         vector3f( 0, 1, 0 ),  vector3f( -1, 0, 0 ), vector3f( 1, 0, 0 ) };

    vector2f results2d[6];
    vector2f results2dNorm[6];
    vector3f results3d[6];

    for( size_t i = 0; i < 6; ++i ) {
        results2d[i] = to_panoramic_yup( axes[i], size, hfov );
        results2dNorm[i] = vector2f( results2d[i].x / size.xsize, results2d[i].y / size.ysize );
        results3d[i] = from_panoramic_yup( results2d[i], size, hfov );
    }

    ASSERT_NEAR( results2dNorm[0].x, 0.5f, 0.00001f );
    ASSERT_NEAR( results2dNorm[0].y, 0.5f, 0.00001f );
    ASSERT_NEAR( results2dNorm[1].x, 1.0f, 0.00001f );
    ASSERT_NEAR( results2dNorm[1].y, 0.5f, 0.00001f );
    ASSERT_NEAR( results2dNorm[2].x, 0.5f, 0.00001f );
    ASSERT_NEAR( results2dNorm[2].y, 0.0f, 0.00001f );
    ASSERT_NEAR( results2dNorm[3].x, 0.5f, 0.00001f );
    ASSERT_NEAR( results2dNorm[3].y, 1.0f, 0.00001f );
    ASSERT_NEAR( results2dNorm[4].x, 0.25f, 0.00001f );
    ASSERT_NEAR( results2dNorm[4].y, 0.5f, 0.00001f );
    ASSERT_NEAR( results2dNorm[5].x, 0.75f, 0.00001f );
    ASSERT_NEAR( results2dNorm[5].y, 0.5f, 0.00001f );

    for( int i = 0; i < 6; ++i ) {
        for( int j = 0; j < 3; ++j ) {
            ASSERT_NEAR( results3d[i][j], axes[i][j], 0.00001f );
        }
    }
}

TEST( GraphicsTest, PanoramicCameraScan ) {

    size2f size( 1024, 512 );
    float hfov = float( 2.0 * M_PI );

    size_t yDropOff = 25;
    size_t xQuarter = size_t( size.xsize / 4 );

    float yMonotoneCheck = -1.00001f;

    // avoid testing too high or low, since the top/bottom rows will likely be nearly uniform x-coords
    for( size_t y = yDropOff; y < ( size.ysize - yDropOff ); ++y ) {

        float xMonotoneCheck = 0.00001f;

        for( size_t x = 0; x < xQuarter; ++x ) {
            vector3f conv = from_panoramic_yup( vector2f( float( x ), float( y ) ), size, hfov );

            ASSERT_LT( conv.x, xMonotoneCheck );
            if( xMonotoneCheck > conv.x ) {
                xMonotoneCheck = conv.x;
            } else {
                break;
            }

            ASSERT_LT( yMonotoneCheck, conv.y );
            ASSERT_LE( -0.00001f, conv.z );

            vector2f backOut = to_panoramic_yup( conv, size, hfov );

            if( x != 0 ) {
                ASSERT_EQ( x, frantic::math::round( backOut.x ) );
            } else {
                ASSERT_EQ( size.xsize, frantic::math::round( backOut.x ) );
            }

            ASSERT_EQ( y, frantic::math::round( backOut.y ) );
        }

        ASSERT_LT( -1.00001f, xMonotoneCheck );

        xMonotoneCheck = -1.00001f;
        for( size_t x = xQuarter; x < 3 * xQuarter; ++x ) {
            vector3f conv = from_panoramic_yup( vector2f( float( x ), float( y ) ), size, hfov );

            ASSERT_LT( xMonotoneCheck, conv.x );
            if( xMonotoneCheck < conv.x ) {
                xMonotoneCheck = conv.x;
            } else {
                break;
            }

            ASSERT_LT( yMonotoneCheck, conv.y );
            ASSERT_LE( conv.z, 0.00001f );

            vector2f backOut = to_panoramic_yup( conv, size, hfov );

            if( x != 0 ) {
                ASSERT_EQ( x, frantic::math::round( backOut.x ) );
            } else {
                ASSERT_EQ( size.xsize, frantic::math::round( backOut.x ) );
            }

            ASSERT_EQ( y, frantic::math::round( backOut.y ) );
        }

        ASSERT_LT( xMonotoneCheck, 1.00001f );

        xMonotoneCheck = 1.00001f;

        for( size_t x = 3 * xQuarter; x < size.xsize; ++x ) {
            vector3f conv = from_panoramic_yup( vector2f( (float)x, (float)y ), size, hfov );

            ASSERT_LT( conv.x, xMonotoneCheck );
            if( xMonotoneCheck > conv.x ) {
                xMonotoneCheck = conv.x;
            } else {
                break;
            }

            ASSERT_LT( yMonotoneCheck, conv.y );
            ASSERT_LE( -0.00001f, conv.z );

            vector2f backOut = to_panoramic_yup( conv, size, hfov );

            if( x != 0 ) {
                ASSERT_EQ( x, frantic::math::round( backOut.x ) );
            } else {
                ASSERT_EQ( size.xsize, frantic::math::round( backOut.x ) );
            }

            ASSERT_EQ( y, frantic::math::round( backOut.y ) );
        }

        ASSERT_LT( xMonotoneCheck, 1.00001f );

        yMonotoneCheck = from_panoramic_yup( vector2f( xQuarter * 2.f, (float)y ), size, hfov ).y;
    }
    ASSERT_LT( yMonotoneCheck, 1.00001f );
}
