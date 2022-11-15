// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest-helper.h"

#include <frantic/graphics/cubeface.hpp>
#include <frantic/graphics/vector3f.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::math;

TEST( Vector3f, MemoryLayout ) {
    // Make sure it has the memory layout we expect (and the channel code depends on)
    vector3f v;
    EXPECT_EQ( 12u, sizeof( v ) );
    EXPECT_EQ( 0, ( (char*)&v.x ) - ( (char*)&v ) );
    EXPECT_EQ( 4, ( (char*)&v.y ) - ( (char*)&v ) );
    EXPECT_EQ( 8, ( (char*)&v.z ) - ( (char*)&v ) );
    EXPECT_EQ( 0, ( (char*)&v[0] ) - ( (char*)&v ) );
    EXPECT_EQ( 4, ( (char*)&v[1] ) - ( (char*)&v ) );
    EXPECT_EQ( 8, ( (char*)&v[2] ) - ( (char*)&v ) );
}

TEST( Vector3fd, MemoryLayout ) {
    // Make sure it has the memory layout we expect (and the channel code depends on)
    vector3fd v;
    EXPECT_EQ( 24u, sizeof( v ) );
    EXPECT_EQ( 0, ( (char*)&v.x ) - ( (char*)&v ) );
    EXPECT_EQ( 8, ( (char*)&v.y ) - ( (char*)&v ) );
    EXPECT_EQ( 16, ( (char*)&v.z ) - ( (char*)&v ) );
    EXPECT_EQ( 0, ( (char*)&v[0] ) - ( (char*)&v ) );
    EXPECT_EQ( 8, ( (char*)&v[1] ) - ( (char*)&v ) );
    EXPECT_EQ( 16, ( (char*)&v[2] ) - ( (char*)&v ) );
}

TEST( Vector3f, Basics ) {
    vector3f x = vector3f::from_xaxis();

    EXPECT_EQ( x, vector3f( 1, 0, 0 ) );
    EXPECT_EQ( x.get_magnitude(), 1 );
    EXPECT_EQ( x.get_largest_axis(), 0 );
    EXPECT_EQ( vector3f::cross( x, vector3f::from_yaxis() ), vector3f::from_zaxis() );
    EXPECT_FALSE( x.is_nan() );
    EXPECT_EQ( vector3f::parse( "[1,2,3]" ), vector3f( 1, 2, 3 ) );
}

TEST( Vector3f, Assignment ) {
    vector3f x( 1, 2, 3 );

    x = vector3f( 4, 5, 6 );
    EXPECT_EQ( vector3f( 4, 5, 6 ), x );

    vector3fd y( 1, 2, 3 );

    // There was a regression in this operation on vector3f -> vector3t<> conversion
    y = vector3fd( 4, 5, 6 );
    EXPECT_EQ( vector3fd( 4, 5, 6 ), y );
}

TEST( Vector3fd, EqualityComparison ) {
    // Verify that operator== and operator!= behave as expected
    EXPECT_TRUE( frantic::graphics::vector3fd( 1, 2, 3 ) == frantic::graphics::vector3fd( 1, 2, 3 ) );
    EXPECT_FALSE( frantic::graphics::vector3fd( 1, 2, 3 ) == frantic::graphics::vector3fd( 4, 2, 3 ) );
    EXPECT_FALSE( frantic::graphics::vector3fd( 1, 2, 3 ) == frantic::graphics::vector3fd( 1, 4, 3 ) );
    EXPECT_FALSE( frantic::graphics::vector3fd( 1, 2, 3 ) == frantic::graphics::vector3fd( 1, 2, 4 ) );

    EXPECT_FALSE( frantic::graphics::vector3fd( 1, 2, 3 ) != frantic::graphics::vector3fd( 1, 2, 3 ) );
    EXPECT_TRUE( frantic::graphics::vector3fd( 1, 2, 3 ) != frantic::graphics::vector3fd( 4, 2, 3 ) );
    EXPECT_TRUE( frantic::graphics::vector3fd( 1, 2, 3 ) != frantic::graphics::vector3fd( 1, 4, 3 ) );
    EXPECT_TRUE( frantic::graphics::vector3fd( 1, 2, 3 ) != frantic::graphics::vector3fd( 1, 2, 4 ) );
}

TEST( Vector3f, NonFiniteString ) {
    const frantic::graphics::vector3f out( std::numeric_limits<float>::infinity(),
                                           std::numeric_limits<float>::quiet_NaN(), 0 );

    const std::string s = out.str();

    const frantic::graphics::vector3f in = frantic::graphics::vector3f::parse( s );

    EXPECT_TRUE( frantic::math::is_infinite( in.x ) );
    EXPECT_TRUE( frantic::math::is_nan( in.y ) );
    EXPECT_EQ( in.z, 0 );
}

TEST( Vector3f, CubeFace ) {
    // Make sure all six directions go to the right cube face indices
    EXPECT_EQ( get_cube_face( vector3f( 1, 0, 0 ) ), 0 );
    EXPECT_EQ( get_cube_face( vector3f( -1, 0, 0 ) ), 1 );
    EXPECT_EQ( get_cube_face( vector3f( 0, 1, 0 ) ), 2 );
    EXPECT_EQ( get_cube_face( vector3f( 0, -1, 0 ) ), 3 );
    EXPECT_EQ( get_cube_face( vector3f( 0, 0, 1 ) ), 4 );
    EXPECT_EQ( get_cube_face( vector3f( 0, 0, -1 ) ), 5 );

    // Make sure that the thresholded version catches the correct cube faces
    cube_face::default_cube_face faces[6];
    int faceCount = get_cube_faces( vector3f( 1, 1, 0 ), 0.9f, faces );
    EXPECT_EQ( faceCount, 2 );
    if( faceCount == 2 ) {
        EXPECT_TRUE( faces[0] == cube_face::CF_X_POS || faces[1] == cube_face::CF_X_POS );
        EXPECT_TRUE( faces[0] == cube_face::CF_Y_POS || faces[1] == cube_face::CF_Y_POS );
    } else {
        cerr << endl << "Got Faces: ";
        for( int i = 0; i < faceCount; ++i ) {
            cerr << faces[i] << " ";
        }
        cerr << endl;
    }
    faceCount = get_cube_faces( vector3f( 1, -1, -1 ), 0.9f, faces );
    EXPECT_EQ( faceCount, 3 );
    if( faceCount == 3 ) {
        EXPECT_TRUE( faces[0] == cube_face::CF_X_POS || faces[1] == cube_face::CF_X_POS ||
                     faces[2] == cube_face::CF_X_POS );
        EXPECT_TRUE( faces[0] == cube_face::CF_Y_NEG || faces[1] == cube_face::CF_Y_NEG ||
                     faces[2] == cube_face::CF_Y_NEG );
        EXPECT_TRUE( faces[0] == cube_face::CF_Z_NEG || faces[1] == cube_face::CF_Z_NEG ||
                     faces[2] == cube_face::CF_Z_NEG );
    } else {
        cerr << endl << "Got Faces: ";
        for( int i = 0; i < faceCount; ++i ) {
            cerr << faces[i] << " ";
        }
        cerr << endl;
    }
}

TEST( Vector3f, Normalize ) {
    vector3f v;
    bool success;

    // x axis
    v.set( 2, 0, 0 );
    success = v.normalize();
    EXPECT_TRUE( success );
    EXPECT_EQ( vector3f( 1, 0, 0 ), v );

    // negative x axis
    v.set( -2, 0, 0 );
    success = v.normalize();
    EXPECT_TRUE( success );
    EXPECT_EQ( vector3f( -1, 0, 0 ), v );

    // y axis
    v.set( 0, 2, 0 );
    success = v.normalize();
    EXPECT_TRUE( success );
    EXPECT_EQ( vector3f( 0, 1, 0 ), v );

    // z axis
    v.set( 0, 0, 2 );
    success = v.normalize();
    EXPECT_TRUE( success );
    EXPECT_EQ( vector3f( 0, 0, 1 ), v );

    // multiple axes
    v.set( 1, 1, 1 );
    success = v.normalize();
    EXPECT_TRUE( success );
    EXPECT_FLOAT_EQ( 1 / sqrt( 3.f ), v.x );
    EXPECT_FLOAT_EQ( 1 / sqrt( 3.f ), v.y );
    EXPECT_FLOAT_EQ( 1 / sqrt( 3.f ), v.z );

    // "small" value that didn't work previously
    v.set( 0.000001f, 0, 0 );
    success = v.normalize();
    EXPECT_TRUE( success );
    EXPECT_EQ( vector3f( 1, 0, 0 ), v );

    v.set( 0.000001f, 0.000001f, 0 );
    success = v.normalize();
    EXPECT_TRUE( success );
    EXPECT_FLOAT_EQ( 1 / sqrt( 2.f ), v.x );
    EXPECT_FLOAT_EQ( 1 / sqrt( 2.f ), v.y );
    EXPECT_FLOAT_EQ( 0, v.z );

    // The following values currently don't succeed, so we don't care about
    // the exact output value.  However, due to how we defined the function,
    // we do expect the output to be a unit vector.

    v.set( 0, 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    // "large" value that previously resulted in a zero vector
    v.set( 1e20f, 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( std::numeric_limits<float>::min(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( -std::numeric_limits<float>::min(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( std::numeric_limits<float>::max(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( -std::numeric_limits<float>::max(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( std::numeric_limits<float>::infinity(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( -std::numeric_limits<float>::infinity(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( std::numeric_limits<float>::denorm_min(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( std::numeric_limits<float>::quiet_NaN(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );

    v.set( std::numeric_limits<float>::signaling_NaN(), 0, 0 );
    success = v.normalize();
    EXPECT_FALSE( success );
    EXPECT_FLOAT_EQ( 1, v.get_magnitude() );
}

TEST( Vector3fd, Normalize ) {
    vector3fd v( 0.000001, 0, 0 );
    v.normalize();
    EXPECT_EQ( vector3fd( 1, 0, 0 ), v );
}

TEST( Vector3f, Dot ) {
    EXPECT_FLOAT_EQ( 1, vector3f::dot( vector3f( 1, 0, 0 ), vector3f( 1, 0, 0 ) ) );
    EXPECT_FLOAT_EQ( 1, vector3f::dot( vector3f( 1, 2, 0 ), vector3f( 1, 0, 4 ) ) );
    EXPECT_FLOAT_EQ( 0, vector3f::dot( vector3f( 1, 1, 0 ), vector3f( 0, 0, 1 ) ) );
}

TEST( Vector3f, Cross ) {
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 0 ), vector3f::cross( vector3f( 1, 0, 0 ), vector3f( 1, 0, 0 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 1, 0, 0 ), vector3f::cross( vector3f( 0, 1, 0 ), vector3f( 0, 0, 1 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( -1, 0, 0 ), vector3f::cross( vector3f( 0, 0, 1 ), vector3f( 0, 1, 0 ) ) );
}

TEST( Vector3f, Angle ) {
    EXPECT_FLOAT_EQ( 0, vector3f::angle( vector3f( 1, 0, 0 ), vector3f( 1, 0, 0 ) ) );
    EXPECT_FLOAT_EQ( float( M_PI ), vector3f::angle( vector3f( 0, 1, 0 ), vector3f( 0, -1, 0 ) ) );
    EXPECT_FLOAT_EQ( float( M_PI / 2 ), vector3f::angle( vector3f( 0, -1, 0 ), vector3f( 0, 0, 1 ) ) );

    EXPECT_FLOAT_EQ( float( M_PI / 2 ),
                     vector3f::angle( vector3f( 0, 0, 0 ), vector3f( 0, 0, 7 ), vector3f( 0, -5, 7 ) ) );
    EXPECT_FLOAT_EQ( float( M_PI / 4 ),
                     vector3f::angle( vector3f( 0, 0, 0 ), vector3f( 0, 1, 0 ), vector3f( -1, 0, 0 ) ) );
}

TEST( Vector3f, AngleTriangleSum ) {
    // Triangle Test: make sure angles add up to 180 degrees
    std::vector<vector3f> p1s;
    p1s.push_back( vector3f( 0, 3, 5 ) );
    p1s.push_back( vector3f( 0, 0, -4 ) );
    p1s.push_back( vector3f( 10, 20, 30 ) );

    std::vector<vector3f> p2s;
    p2s.push_back( vector3f( 4, -6, 3 ) );
    p2s.push_back( vector3f( -9, 17, 1 ) );
    p2s.push_back( vector3f( -3, -2, -5 ) );

    std::vector<vector3f> p3s;
    p3s.push_back( vector3f( -9, 1, -2 ) );
    p3s.push_back( vector3f( -8, -5, 2 ) );
    p3s.push_back( vector3f( 7, 0, 0 ) );

    for( std::vector<vector3f>::const_iterator iter1 = p1s.begin(); iter1 != p1s.end(); ++iter1 ) {
        const vector3f& p1 = *iter1;
        for( std::vector<vector3f>::const_iterator iter2 = p2s.begin(); iter2 != p2s.end(); ++iter2 ) {
            const vector3f& p2 = *iter2;
            for( std::vector<vector3f>::const_iterator iter3 = p3s.begin(); iter3 != p3s.end(); ++iter3 ) {
                const vector3f& p3 = *iter3;

                const float a2 = vector3f::angle( p1, p2, p3 );
                const float a3 = vector3f::angle( p2, p3, p1 );
                const float a1 = vector3f::angle( p3, p1, p2 );

                EXPECT_FLOAT_EQ( float( M_PI ), a1 + a2 + a3 );
            }
        }
    }
}

TEST( Vector3f, SelfProjection ) {
    // Project vector onto itself
    std::vector<vector3f> toTest;
    toTest.push_back( vector3f( 0, 0, 1 ) );
    toTest.push_back( vector3f( 0, 3, 1 ) );
    toTest.push_back( vector3f( -4, 2, 1 ) );

    for( std::vector<vector3f>::const_iterator iter = toTest.begin(); iter != toTest.end(); ++iter ) {
        const vector3f& expected = *iter;
        EXPECT_VECTOR3F_EQ( expected, vector3f::project( expected, expected ) );
        EXPECT_VECTOR3F_EQ( expected, vector3f::project( expected, expected * 13 ) );
        EXPECT_VECTOR3F_EQ( expected * 13, vector3f::project( expected * 13, expected ) );
    }
}

TEST( Vector3f, Projection ) {
    // Project parallel / anti parallel vectors
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 1 ), vector3f::project( vector3f( 0, 0, 1 ), vector3f( 0, 0, -1 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 2, 0 ), vector3f::project( vector3f( 0, 2, 0 ), vector3f( 0, 3, 0 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 0, -2, 0 ), vector3f::project( vector3f( 0, -2, 0 ), vector3f( 0, 1, 0 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 0, -5, 0 ), vector3f::project( vector3f( 0, -5, 0 ), vector3f( 0, 6, 0 ) ) );

    // Project orthogonal vectors
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 0 ), vector3f::project( vector3f( 0, 0, 1 ), vector3f( 0, 1, 0 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 0 ), vector3f::project( vector3f( 0, 0, 3 ), vector3f( 0, -4, 0 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 0 ), vector3f::project( vector3f( -4, 9, 0 ), vector3f( 0, 0, -3 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 0 ), vector3f::project( vector3f( -4, 9, 0 ), vector3f( 0, 0, -3 ) ) );

    // Arbitrary vectors
    EXPECT_VECTOR3F_EQ( vector3f( 0, 5, 0 ), vector3f::project( vector3f( 5, 5, 0 ), vector3f( 0, 2, 0 ) ) );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 5, 0 ), vector3f::project( vector3f( -5, 5, 5 ), vector3f( 0, 2, 0 ) ) );
}

TEST( Vector3f, OrthogonalProjection ) {
    // Orthogonality / Direction test
    // Make sure projection result goes in the appropriate direction and that it is possible to decompose a vector into
    // 2 orthogonal components
    std::vector<vector3f> toTestSource;
    toTestSource.push_back( vector3f( 0, 0, 1 ) );
    toTestSource.push_back( vector3f( 2, 5, -9 ) );
    toTestSource.push_back( vector3f( -4, -7, 6 ) );

    std::vector<vector3f> toTestDestination;
    toTestSource.push_back( vector3f( 0, 1, 0 ) );
    toTestSource.push_back( vector3f( -9, 4, 3 ) );
    toTestSource.push_back( vector3f( 2, 0, -5 ) );

    for( std::vector<vector3f>::const_iterator iter1 = toTestSource.begin(); iter1 != toTestSource.end(); ++iter1 ) {
        const vector3f& testSource = *iter1;
        for( std::vector<vector3f>::const_iterator iter2 = toTestDestination.begin(); iter2 != toTestDestination.end();
             ++iter2 ) {
            const vector3f& testDestination = *iter2;

            const vector3f testProjection = vector3f::project( testSource, testDestination );
            const vector3f testPerpendicular = testSource - testProjection;

            EXPECT_FLOAT_EQ( 0, vector3f::dot( testProjection, testPerpendicular ) );
            EXPECT_FLOAT_EQ( float( M_PI / 2 ), vector3f::angle( testProjection, testPerpendicular ) );

            if( testProjection.get_magnitude_squared() != 0 ) {
                EXPECT_VECTOR3F_EQ( testDestination.to_normalized(), testProjection.to_normalized() );
            }
        }
    }
}
