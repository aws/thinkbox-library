// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics/plane3f.hpp>
#include <frantic/graphics/spherical_coords.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <half.h>

#include "gtest-helper.h"

using namespace std;
using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::math;

TEST( Transform4f, Basics ) {
    transform4f xform;

    // Ensure that the default constructor produces the identity matrix
    EXPECT_TRUE( xform.is_identity() );
    EXPECT_EQ( xform.determinant(), 1 );
    EXPECT_EQ( xform.determinant_double(), 1 );
    EXPECT_EQ( xform.determinant_3x3(), 1 );

    // Test the zero matrix
    xform.set_to_zero();
    EXPECT_TRUE( xform.is_zero() );
    EXPECT_EQ( xform.determinant(), 0 );

    // Test the inverse function
    xform.set_to_identity();
    transform4f xformInverse = xform.to_inverse();
    EXPECT_TRUE( xformInverse.is_identity() );
    for( int decomposeTestNumber = 0; decomposeTestNumber < 500; ++decomposeTestNumber ) {
        float data[16];
        for( int i = 0; i < 16; ++i )
            data[i] = 2 * ( float( rand() ) / RAND_MAX ) - 1;
        xform.set_to_array( data );

        xformInverse = xform.to_inverse();

        // Have to choose the appropriate epsilon for floating point precision
        float comparisonEpsilon =
            0.000001f * ( std::max )( xformInverse.max_abs_component(), xform.max_abs_component() );

        transform4f product = xform * xformInverse;
        EXPECT_TRUE( product.is_identity( comparisonEpsilon ) )
            << "Produced incorrect inverse matrix:" << endl
            << " xform: " << xform << endl
            << " inverse: " << xformInverse << endl
            << " product: " << product << endl
            << " xform determinant: " << xform.determinant() << endl;
    }

    // TODO Test equivalence of vector3f::get_cube_face_coordinate to transform4f::from_cubeface
    transform4f transCubeFace[] = {
        transform4f::from_cubeface( cube_face::CF_X_POS ).to_inverse(), // cube_face::CF_RIGHT  is 0
        transform4f::from_cubeface( cube_face::CF_X_NEG ).to_inverse(), // cube_face::CF_LEFT   is 1
        transform4f::from_cubeface( cube_face::CF_Y_POS ).to_inverse(), // cube_face::CF_TOP    is 2
        transform4f::from_cubeface( cube_face::CF_Y_NEG ).to_inverse(), // cube_face::CF_BOTTOM is 3
        transform4f::from_cubeface( cube_face::CF_Z_POS ).to_inverse(), // cube_face::CF_REAR   is 4
        transform4f::from_cubeface( cube_face::CF_Z_NEG ).to_inverse(), // cube_face::CF_FRONT  is 5
    };
    transform4f transPerspective =
        transform4f::from_fov_perspective( float( M_PI ) / 2.0f, float( M_PI ) / 2.0f, 1.0f, 100.0f );

    for( int cubefaceTestNumber = 0; cubefaceTestNumber < 500; ++cubefaceTestNumber ) {
        float comparisonEpsilon = 0.000001f;

        // Generate random vectors with x, y, z between [-1, 1]
        vector3f v = vector3f::from_unit_random();

        // Calculate the cube face of current vector
        cube_face::default_cube_face cubeFace = get_cube_face( v );

        // Use transform4f::from_cubeface to rotate the vector so that it is facing forward
        vector3f transCube = ( transCubeFace[cubeFace] * v );
        // Calculate the projection normalized coordinate in 2d space
        bool isValid = true;
        frantic::graphics2d::vector2f transResult = transCube.to_perspective_projection(
            frantic::graphics2d::size2f( 1.0f, 1.0f ), tan( float( M_PI ) / 2.0f / 2 ), isValid );
        EXPECT_TRUE( isValid );

        // Use vector3f::get_cube_face_coordinate to calculate the normalized coordinate on its
        // corresponding cube face
        frantic::graphics2d::vector2f v3fResult = get_cube_face_coordinate( v, cubeFace );
        // Normalize from [-1, 1] to [0, 1] to compare with transform4f::from_cubeface result
        v3fResult.x = ( ( v3fResult.x + 1.0f ) / 2.0f );
        v3fResult.y = ( ( v3fResult.y + 1.0f ) / 2.0f );

        // DEBUGGING OUTPUT
        /*printf("\nV = ( %f, %f, %f )\n", v.x, v.y, v.z );
        printf("transResult = ( %f, %f )\n", transResult.x, transResult.y );
        printf("v3fResult = ( %f, %f )\n", v3fResult.x, v3fResult.y );*/

        // Calculate x and y difference between results
        float xdiff = fabsf( transResult.x - v3fResult.x );
        float ydiff = fabsf( transResult.y - v3fResult.y );

        // Ensure results are acceptable
        EXPECT_LT( xdiff, comparisonEpsilon );
        EXPECT_LT( ydiff, comparisonEpsilon );
    }

    transform4f ap = transform4f::from_axis_permutation( vector3( 0, 1, 2 ) );
    EXPECT_TRUE( ap.is_identity() );
    ap = transform4f::from_axis_permutation( vector3( 1, 0, 2 ) );
    EXPECT_EQ( vector3f( 1, 0, 2 ), ap * vector3f( 0, 1, 2 ) );
    ap = transform4f::from_axis_permutation( vector3( 1, 2, 0 ) );
    EXPECT_EQ( vector3f( 2, 0, 1 ), ap * vector3f( 0, 1, 2 ) );
    ap = transform4f::from_axis_permutation( vector3( 0, 2, 1 ) );
    EXPECT_EQ( vector3f( 0, 2, 1 ), ap * vector3f( 0, 1, 2 ) );
}

TEST( Transform4f, ConvertFromSmallerFloatType ) {
    transform4t<half> halfTransform;
    for( int i = 0; i < 16; ++i ) {
        halfTransform[i] = static_cast<float>( i );
    }
    transform4fd doubleTransform( halfTransform );
    // Can't use EXPECT_EQ( doubleTransform, halfTransform ), because that
    // comparison would also go through the transform4fd( transform4f<half> )
    // constructor.
    for( int i = 0; i < 16; ++i ) {
        EXPECT_EQ( doubleTransform[i], halfTransform[i] );
    }
}

TEST( Transform4f, LookAt ) {
    transform4f xform;

    xform = transform4f::from_look_dir( vector3f( 0 ), vector3f( 0, 1, 0 ), vector3f( 0, 0, 1 ) );
    EXPECT_EQ( xform, transform4f( 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1 ) );
    xform = transform4f::from_look_at( vector3f( 0 ), vector3f( 0, 3, 0 ), vector3f( 0, 0, 1 ) );
    EXPECT_EQ( xform, transform4f( 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1 ) );
}

TEST( Transform4f, FromNormal ) {
    transform4f xform;
    // Test the transform4f::from_normal function
    xform = transform4f::from_normal_z( vector3f::from_xaxis() );
    EXPECT_EQ( xform.determinant(), 1 );
    EXPECT_TRUE( xform.is_orthogonal() );
    xform = transform4f::from_normal_z( vector3f::from_yaxis() );
    EXPECT_EQ( xform.determinant(), 1 );
    EXPECT_TRUE( xform.is_orthogonal() );
    xform = transform4f::from_normal_z( vector3f::from_zaxis() );
    EXPECT_TRUE( xform.is_identity() );
    EXPECT_TRUE( xform.is_orthogonal() );
    for( int i = 0; i < 30; ++i ) {
        vector3f normal = vector3f::from_unit_random();
        xform = transform4f::from_normal_z( normal );
        float det = xform.determinant();
        if( fabs( det - 1 ) >= 0.0001f ) {
            cerr << endl << "transform4f::from_normal( " << normal << " ) == " << xform << endl;
            EXPECT_LT( fabsf( det - 1 ), 0.0001f );
        }
        if( !xform.is_orthogonal( 0.00005f ) ) {
            cerr << endl << "transform4f::from_normal( " << normal << " ) == " << xform << endl;
            EXPECT_TRUE( xform.is_orthogonal( 0.00005f ) );
        }
    }
}

TEST( Transform4f, FromRotations ) {
    transform4f xform;
    // Test transform4f::from_*_rotation
    xform = transform4f::from_x_rotation( math::degrees_to_radians( 90.f ) );
    EXPECT_TRUE( xform.is_orthogonal() );
    EXPECT_EQ( xform.determinant(), 1 );
    EXPECT_EQ( xform * vector3f::from_xaxis(), vector3f::from_xaxis() );
    EXPECT_TRUE( ( xform * vector3f::from_yaxis() ).is_equal( vector3f::from_zaxis(), 0.000001f ) );
    EXPECT_TRUE( ( xform * vector3f::from_zaxis() ).is_equal( -vector3f::from_yaxis(), 0.000001f ) );
    xform = transform4f::from_y_rotation( math::degrees_to_radians( 90.f ) );
    EXPECT_TRUE( xform.is_orthogonal() );
    EXPECT_EQ( xform.determinant(), 1 );
    EXPECT_TRUE( ( xform * vector3f::from_xaxis() ).is_equal( -vector3f::from_zaxis(), 0.000001f ) );
    EXPECT_EQ( xform * vector3f::from_yaxis(), vector3f::from_yaxis() );
    EXPECT_TRUE( ( xform * vector3f::from_zaxis() ).is_equal( vector3f::from_xaxis(), 0.000001f ) );
    xform = transform4f::from_z_rotation( math::degrees_to_radians( 90.f ) );
    EXPECT_TRUE( xform.is_orthogonal() );
    EXPECT_EQ( xform.determinant(), 1 );
    EXPECT_TRUE( ( xform * vector3f::from_xaxis() ).is_equal( -vector3f::from_yaxis(), 0.000001f ) );
    EXPECT_TRUE( ( xform * vector3f::from_yaxis() ).is_equal( vector3f::from_xaxis(), 0.000001f ) );
    EXPECT_EQ( xform * vector3f::from_zaxis(), vector3f::from_zaxis() );
}

TEST( Transform4f, Decompose ) {
    unsigned seed = (unsigned)time( NULL );
    srand( seed );
    std::cout << "Running test testTransform4f_decompose with seed " << seed << std::endl;
    transform4f xform;
    // Test transform4f::decompose
    for( int decomposeTestNumber = 0; decomposeTestNumber < 500; ++decomposeTestNumber ) {
        float data[16];
        for( int i = 0; i < 16; ++i )
            data[i] = 10 * ( float( rand() ) / RAND_MAX ) - 5;
        xform.set_to_array( data );

        transform4f persp, rotate, stretch;
        vector3f translate;
        xform.decompose( persp, translate, rotate, stretch );

        transform4f xformCompare = persp * transform4f::from_translation( translate ) * rotate * stretch;

        // Have to choose the appropriate epsilon for floating point precision
        float comparisonEpsilon = persp.max_abs_component();
        comparisonEpsilon = ( std::max )( comparisonEpsilon, translate.max_abs_component() );
        comparisonEpsilon = ( std::max )( comparisonEpsilon, rotate.max_abs_component() );
        comparisonEpsilon = ( std::max )( comparisonEpsilon, stretch.max_abs_component() );

        comparisonEpsilon *= 0.0001f;

        if( !xform.equals( xformCompare, comparisonEpsilon ) || !rotate.is_orthogonal( 0.000001f ) ) {
            cerr << endl << "xform: " << xform << endl;
            cerr << "xformCompare: " << xformCompare << endl;
            cerr << "persp: " << persp << endl;
            cerr << "translate: " << translate << endl;
            cerr << "rotate: " << rotate << endl;
            cerr << "stretch: " << stretch << endl;
            cerr << "comparisonEpsilon: " << comparisonEpsilon << endl;
        }

        EXPECT_TRUE( rotate.is_orthogonal( 0.000001f ) );

        EXPECT_TRUE( xform.equals( xformCompare, comparisonEpsilon ) );
    }
}

TEST( Transform4f, DecomposeSimplified ) {
    std::vector<vector3f> positions;
    positions.push_back( vector3f( 0, 0, 0 ) );
    positions.push_back( vector3f( 2, 0, 0 ) );
    positions.push_back( vector3f( 6, 3, -4 ) );

    std::vector<quat4f> rotations;
    rotations.push_back( quat4f() );
    rotations.push_back( quat4f( quat4fd::from_euler_angles( 0, 0, M_PI / 4 ) ) );
    rotations.push_back( quat4f( quat4fd::from_euler_angles( 0, -M_PI / 2, 0 ) ) );
    rotations.push_back( quat4f( quat4fd::from_euler_angles( M_PI / 2, -M_PI / 3, M_PI / 4 ) ) );

    std::vector<vector3f> scales;
    scales.push_back( vector3f( 1, 1, 1 ) );
    scales.push_back( vector3f( 1, 1, 4 ) );
    scales.push_back( vector3f( 1, 0.375f, 1 ) );
    scales.push_back( vector3f( 3, 8, -0.75f ) );

    for( std::vector<vector3f>::const_iterator iter1 = positions.begin(); iter1 != positions.end(); ++iter1 ) {
        const vector3f& expectedTranslation = *iter1;
        for( std::vector<quat4f>::const_iterator iter2 = rotations.begin(); iter2 != rotations.end(); ++iter2 ) {
            const quat4f& expectedRotation = *iter2;

            // No Scale
            {
                const transform4f baseTransform = transform4f::to_transform( expectedTranslation, expectedRotation );

                vector3f actualTranslation;
                quat4f actualRotation;
                baseTransform.decompose( actualTranslation, actualRotation );

                EXPECT_VECTOR3F_EQ( expectedTranslation, actualTranslation );
                EXPECT_QUAT4F_NEAR( expectedRotation, actualRotation, 0.000001f * actualRotation.max_abs_component() );
            }

            // With Scale
            for( std::vector<vector3f>::const_iterator iter3 = scales.begin(); iter3 != scales.end(); ++iter3 ) {
                const vector3f& expectedScale = *iter3;

                const transform4f baseTransform =
                    transform4f::to_transform( expectedTranslation, expectedRotation, expectedScale );

                vector3f actualTranslation;
                quat4f actualRotation;
                vector3f actualScale;
                baseTransform.decompose( actualTranslation, actualRotation, actualScale );

                EXPECT_VECTOR3F_EQ( expectedTranslation, actualTranslation );
                EXPECT_QUAT4F_NEAR( expectedRotation, actualRotation, 0.000001f * actualRotation.max_abs_component() );
                EXPECT_VECTOR3F_EQ( expectedScale, actualScale );
            }
        }
    }
}

TEST( Transform4f, PlanarTransform ) {
    std::srand( 8458973 );

    boundbox3f bounds( vector3f( -100.0f, -100.0f, -100.0f ), vector3f( 100.0f, 100.0f, 100.0f ) );

    std::vector<vector3f> points( 20 );

    // first, take 3 points to establish a plane
    for( size_t i = 0; i < 3; ++i ) {
        points[i] = bounds.random_vector();
    }

    plane3f plane( plane3f::from_triangle( points[0], points[1], points[2] ) );

    for( size_t i = 3; i < points.size(); ++i ) {
        points[i] = plane.project_onto_plane( bounds.random_vector() );
    }

    transform4f xform = plane.get_planar_transform();

    std::vector<vector3f> resultPoints( points.size() );

    transform4f iXform = xform.to_inverse();

    for( size_t i = 0; i < points.size(); ++i ) {
        EXPECT_NEAR( 0.0f, plane.get_signed_distance_to_plane( points[i] ), 0.0001f );

        resultPoints[i] = xform * points[i];

        // check that the transformation brings it down to the z = 0 plane
        EXPECT_NEAR( 0.0f, resultPoints[i].z, 1.0e-4f );

        // check that the transformation preserves distances
        for( size_t j = 0; j < i; ++j ) {
            double actualDistance = vector3f::distance_double( points[i], points[j] );
            double planarDistance = vector3f::distance_double( resultPoints[i], resultPoints[j] );
            EXPECT_NEAR( actualDistance, planarDistance, 1.0e-4 );
        }

        // check that the transformation inverts correctly
        vector3f backAgain = iXform * resultPoints[i];

        EXPECT_NEAR( points[i].x, backAgain.x, 1.0e-4f );
        EXPECT_NEAR( points[i].y, backAgain.y, 1.0e-4f );
        EXPECT_NEAR( points[i].z, backAgain.z, 1.0e-4f );
    }
}

TEST( Transform4f, ToAxisAngleConversion ) {
    std::srand( 6982334 );

    for( size_t i = 0; i < 1; ++i ) { // FIXME: to 100
        float angle = linearConvertRange( rand() / float( RAND_MAX ), 0.0f, 1.0f, 0.0f, float( M_PI * 2.0f ) );
        vector3f axis = spherical_coords::from_random().to_vector3f().to_normalized();

        transform4f rotated = transform4f::from_angle_axis( angle, axis );
        std::pair<float, vector3f> backOut = rotated.to_angle_axis();
        transform4f rotated2 = transform4f::from_angle_axis( backOut.first, backOut.second );

        vector3f solvedAxis = backOut.second;
        float solvedAngle = backOut.first;

        // The angle-axis representation is not unique, rather, the
        // method may have solved for the inverted axis.
        if( vector3f::dot( axis, backOut.second ) < 0.0f ) {
            solvedAxis = -solvedAxis;
            solvedAngle = float( M_PI * 2.0 ) - solvedAngle;
        }

        EXPECT_NEAR( angle, solvedAngle, 0.00001f );
        for( int a = 0; a < 3; ++a ) {
            EXPECT_NEAR( axis[a], solvedAxis[a], 0.00001f );
        }

        for( int e = 0; e < 16; ++e ) {
            EXPECT_NEAR( rotated[e], rotated2[e], 0.00001f );
        }
    }
}

TEST( Transform4f, EulerAngles ) {
    vector3fd angles( M_PI_4, M_PI_4, M_PI_4 );
    transform4fd fromEuler( transform4fd::from_euler_xyz_rotation( angles ) );
    vector3fd eulerAngles( fromEuler.get_euler_angles() );
    EXPECT_VECTOR3FD_EQ( angles, eulerAngles );
}

TEST( Transform4f, EulerAnglesPI_2 ) {
    vector3fd angles( M_PI_2, M_PI_4, 0 );
    transform4fd fromEuler( transform4fd::from_euler_xyz_rotation( angles ) );
    vector3fd eulerAngles( fromEuler.get_euler_angles() );
    EXPECT_VECTOR3FD_EQ( angles, eulerAngles );
}

TEST( Transform4f, EulerAnglesMinusPI_2 ) {
    vector3fd angles( -M_PI_2, M_PI_4, 0 );
    transform4fd fromEuler( transform4fd::from_euler_xyz_rotation( angles ) );
    vector3fd eulerAngles( fromEuler.get_euler_angles() );
    EXPECT_VECTOR3FD_EQ( angles, eulerAngles );
}
