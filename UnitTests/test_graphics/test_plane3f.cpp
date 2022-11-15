// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics/plane3f.hpp>
#include <frantic/graphics/ray3f.hpp>
#include <frantic/graphics/transform4f.hpp>

#include <boost/random.hpp>
#include <vector>

using namespace std;
using namespace frantic::graphics;

TEST( Plane3f, PlaneDistanceQuery ) {
    plane3f p = plane3f::from_normal_and_point( vector3f( 0, 0, 1 ), vector3f() );

    // Make sure that we can hit the plane properly from a number of different angles
    ray3f ray( vector3f( 1, 0, 1 ), vector3f( 0, 0, -1 ).to_normalized() );
    float distance, desiredDistance;
    distance = (float)p.get_distance_to_intersection( ray );
    desiredDistance = 1;
    EXPECT_TRUE( p.is_coplanar( ray.at( distance ) ) );
    if( fabs( distance - desiredDistance ) > 0.00001f ) {
        cerr << "distance (" << desiredDistance << " desired): " << distance << endl;
        EXPECT_TRUE( fabs( distance - desiredDistance ) < 0.00001f );
    }
    ray.set( vector3f( 1, 0, 1 ), vector3f( 1, 0, -1 ).to_normalized() );
    distance = (float)p.get_distance_to_intersection( ray );
    desiredDistance = (float)sqrt( 2.f );
    EXPECT_TRUE( p.is_coplanar( ray.at( distance ) ) );
    if( fabs( distance - desiredDistance ) > 0.00001f ) {
        cerr << "distance (" << desiredDistance << " desired): " << distance << endl;
        EXPECT_TRUE( fabs( distance - desiredDistance ) < 0.00001f );
    }

    ray.set( vector3f( 1, 1, 1 ), vector3f( 1, 1, -1 ).to_normalized() );
    distance = (float)p.get_distance_to_intersection( ray );
    desiredDistance = (float)sqrt( 3.f );
    EXPECT_TRUE( p.is_coplanar( ray.at( distance ) ) );
    if( fabs( distance - desiredDistance ) > 0.00001f ) {
        cerr << "distance (" << desiredDistance << " desired): " << distance << endl;
        EXPECT_TRUE( fabs( distance - desiredDistance ) < 0.00001f );
    }

    p = plane3f( -3.86641e-010f, 0.999391f, -0.0348996f, 8.73886f );
    desiredDistance = 31.6226f;
    ray.set( vector3f( 0.031356f, 2.97882f, -48.448f ),
             vector3f( -0.0644832f, 0.0941238f, -0.99347f ).to_normalized() );
    distance = (float)p.get_distance_to_intersection( ray );
    EXPECT_TRUE( p.is_coplanar( ray.at( distance ) ) );
    if( fabs( distance - desiredDistance ) > 0.0001f ) {
        cerr << "distance (" << desiredDistance << " desired): " << distance << endl;
        EXPECT_TRUE( fabs( distance - desiredDistance ) < 0.0001f );
    }

    ray.set( vector3f( 4.3797f, 7.57111f, -52.4579f ),
             vector3f( -0.167559f, -0.00126226f, -0.985861f ).to_normalized() );
    distance = (float)p.get_distance_to_intersection( ray );
    desiredDistance = -1;
    if( fabs( distance - desiredDistance ) > 0.00001f ) {
        cerr << "distance (" << desiredDistance << " desired): " << distance << endl;
        EXPECT_TRUE( fabs( distance - desiredDistance ) < 0.00001f );
    }

    ray.set( vector3f( 5.03853f, 8.71001f, -60.3489f ),
             vector3f( -0.219909f, -0.0323859f, -0.974983f ).to_normalized() );
    distance = (float)p.get_distance_to_intersection( ray );
    desiredDistance = -1;
    // cerr << "ray.at(distance): " << ray.at(distance) << ", p.get_signed_distance_to_plane(ray.at(distance)): " <<
    // p.get_signed_distance_to_plane(ray.at(distance)) << endl;
    if( fabs( distance - desiredDistance ) > 0.0001f ) {
        cerr << "distance (" << desiredDistance << " desired): " << distance << endl;
        EXPECT_TRUE( fabs( distance - desiredDistance ) < 0.0001f );
    }

    ray.set( vector3f( 5.36794f, 9.27946f, -64.2944f ),
             vector3f( -0.254549f, -0.0531907f, -0.965596f ).to_normalized() );
    distance = (float)p.get_distance_to_intersection( ray );
    desiredDistance = 142.8f;
    EXPECT_TRUE( p.is_coplanar( ray.at( distance ) ) );
    if( fabs( distance - desiredDistance ) > 0.001f ) {
        cerr << "distance (" << desiredDistance << " desired): " << distance << endl;
        EXPECT_TRUE( fabs( distance - desiredDistance ) < 0.001f );
    }

    // ray: (ray (vec3f  ) (vec3f  ) ), plane: (plane3f -3.86641e-010, 0.999391, -0.0348996, 8.73886 ), attenDistance:
    // -1, waterPlaneDistance: 36.9008 ray: (ray (vec3f  ) (vec3f  ) ), plane: (plane3f -3.86641e-010, 0.999391,
    // -0.0348996, 8.73886 ), attenDistance: 142.801, waterPlaneDistance: 33.1734
}

TEST( Plane3f, PlaneLineSegmentIntersection ) {
    plane3f p = plane3f::from_normal_and_point( vector3f( 0, 0, 1 ), vector3f() );

    // Test getting the intersection with line segments
    vector3f isect;
    EXPECT_TRUE( p.get_intersection( vector3f( 0, 0, -1 ), vector3f( 0, 0, 1 ), isect ) );
    EXPECT_EQ( isect, vector3f() );

    EXPECT_TRUE( p.get_intersection( vector3f( 5, 0, -1 ), vector3f( 7, 0, 1 ), isect ) );
    EXPECT_EQ( isect, vector3f( 6, 0, 0 ) );
}

TEST( Plane3f, PlanarTransform ) {
    boost::mt19937 generator( 8458973 );
    boost::uniform_real<float> uniformDistro( 0.0f, 1.0f );
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float>> uniformRandom( generator, uniformDistro );

    std::srand( 8458973 );

    boundbox3f bounds( vector3f( -100.0f, -100.0f, -100.0f ), vector3f( 100.0f, 100.0f, 100.0f ) );

    std::vector<vector3f> points( 20 );

    // first, take 3 points to establish a plane
    for( size_t i = 0; i < 3; ++i ) {
        points[i] = bounds.random_vector( uniformRandom );
    }

    plane3f plane( plane3f::from_triangle( points[0], points[1], points[2] ) );

    for( size_t i = 3; i < points.size(); ++i ) {
        points[i] = plane.project_onto_plane( bounds.random_vector( uniformRandom ) );
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
