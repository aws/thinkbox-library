// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "gtest-helper.h"

#include <frantic/graphics/quat4f.hpp>

#include <boost/random.hpp>
#include <vector>

using namespace std;
using namespace frantic::graphics;

TEST( Quat4f, MetaFunctionHelpers ) {
    // Confirm the metafunction is producing the bigger type
    EXPECT_EQ( 4u, sizeof( frantic::graphics::detail::biggest_type<float, float>::type ) );
    EXPECT_EQ( 8u, sizeof( frantic::graphics::detail::biggest_type<float, double>::type ) );
    EXPECT_EQ( 8u, sizeof( frantic::graphics::detail::biggest_type<double, float>::type ) );
    EXPECT_EQ( 8u, sizeof( frantic::graphics::detail::biggest_type<double, double>::type ) );
}

TEST( Quat4f, BasicRotations ) {
    vector3f v( 1, 0, 0 );
    quat4f q;
    transform4f t;

    // Identity
    q.set_to_identity();
    q.as_transform4f( t );
    EXPECT_EQ( transform4f(), t );
    EXPECT_EQ( v, rotate_point( q, v ) );
    // Other simple rotations
    q = quat4f( 0, 1, 0, 0 );
    q.as_transform4f( t );
    EXPECT_EQ( t * v, rotate_point( q, v ) );
    q = quat4f( 0, 0, 1, 0 );
    q.as_transform4f( t );
    EXPECT_EQ( t * v, rotate_point( q, v ) );
    q = quat4f( 0, 0, 0, 1 );
    q.as_transform4f( t );
    EXPECT_EQ( t * v, rotate_point( q, v ) );
}

TEST( Quat4f, Transform4fConversions ) {
    // Identity
    const transform4f identity = transform4f::identity();

    EXPECT_QUAT4F_EQ( quat4f::from_identity(), quat4f::from_transform4f( identity ) );

    // Random quaternions, test round trip conversion
    {
        const float w = 1;
        const float x = 2;
        const float y = 3;
        const float z = 4;
        const float dist = sqrt( w * w + x * x + y * y + z * z );

        quat4f testQuat( w / dist, x / dist, y / dist, z / dist );
        transform4f testMatrix;
        testQuat.as_transform4f( testMatrix );
        quat4f resultQuat = quat4f::from_transform4f( testMatrix );

        EXPECT_QUAT4F_EQ( testQuat, resultQuat );
    }

    {
        const float w = 3;
        const float x = -2;
        const float y = 1;
        const float z = -4;
        const float dist = sqrt( w * w + x * x + y * y + z * z );

        quat4f testQuat( w / dist, x / dist, y / dist, z / dist );
        transform4f testMatrix;
        testQuat.as_transform4f( testMatrix );
        quat4f resultQuat = quat4f::from_transform4f( testMatrix );

        EXPECT_QUAT4F_EQ( testQuat, resultQuat );
    }
}

TEST( Quat4f, FromAngleAxisIdentities ) {
    // To Identity
    const quat4f identityQuat = quat4f::from_identity();

    EXPECT_QUAT4F_EQ( identityQuat, quat4f::from_angle_axis( 0, vector3f( 0, 0, 1 ) ) );
    EXPECT_QUAT4F_EQ( identityQuat, quat4f::from_angle_axis( 0, vector3f( 0, 1, 0 ) ) );
    EXPECT_QUAT4F_EQ( identityQuat, quat4f::from_angle_axis( 0, vector3f( -1 / sqrt( 2.f ), 1 / sqrt( 2.f ), 0 ) ) );

    // From Identity
    {
        float readAngle;
        vector3f readAxis;
        quat4f quat4;
        quat4.as_angle_axis( readAngle, readAxis );
        EXPECT_FLOAT_EQ( 0, readAngle );

        // Axis does not matter but code assumes z-axis
        EXPECT_VECTOR3F_EQ( vector3( 0, 0, 1 ), readAxis );
    }
}

TEST( Quat4f, FromAngleAxisRoundTrip ) {
    // Angle-Axis to Quaternion to Angle-Axis round trip conversion
    {
        const float expectedAngle = 2.56f;
        const vector3f expectedAxis( 1, 0, 0 );

        const quat4f quat = quat4f::from_angle_axis( expectedAngle, expectedAxis );

        float actualAngle;
        vector3f actualAxis;
        quat.as_angle_axis( actualAngle, actualAxis );

        EXPECT_FLOAT_EQ( expectedAngle, actualAngle );
        EXPECT_VECTOR3F_EQ( expectedAxis, actualAxis );
    }

    {
        const float expectedAngle = 0.675f;
        const vector3f expectedAxis( -1 / sqrt( 3.f ), -1 / sqrt( 3.f ), 1 / sqrt( 3.f ) );

        const quat4f quat = quat4f::from_angle_axis( expectedAngle, expectedAxis );

        float actualAngle;
        vector3f actualAxis;
        quat.as_angle_axis( actualAngle, actualAxis );

        EXPECT_FLOAT_EQ( expectedAngle, actualAngle );
        EXPECT_VECTOR3F_EQ( expectedAxis, actualAxis );
    }

    {
        const float expectedAngle = float( M_PI * 3 / 2 );
        const vector3f expectedAxis( 1, float( -1 / sqrt( 2.f ) ), float( -1 / sqrt( 2.f ) ) );

        const quat4f quat = quat4f::from_angle_axis( expectedAngle, expectedAxis );

        float actualAngle;
        vector3f actualAxis;
        quat.as_angle_axis( actualAngle, actualAxis );

        EXPECT_FLOAT_EQ( expectedAngle, actualAngle );
        EXPECT_VECTOR3F_EQ( expectedAxis, actualAxis );
    }

    // Corner Case
    {
        const float expectedAngle = float( M_PI );
        vector3f expectedAxis( 0, 0, 1 );

        quat4f quat = quat4f::from_angle_axis( expectedAngle, expectedAxis );

        float actualAngle;
        vector3f actualAxis;
        quat.as_angle_axis( actualAngle, actualAxis );

        // Angle sign and z axis direction doesn't matter
        actualAngle = std::abs( actualAngle );
        actualAxis.z = std::abs( actualAxis.z );

        EXPECT_FLOAT_EQ( expectedAngle, actualAngle );
        EXPECT_VECTOR3F_EQ( expectedAxis, actualAxis );
    }
}
