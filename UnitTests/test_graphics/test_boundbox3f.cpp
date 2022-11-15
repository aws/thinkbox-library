// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <gtest/gtest.h>

#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/plane3f.hpp>
#include <frantic/graphics/ray3f.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::math;

TEST( BoundBox3f, Basics ) {
    boundbox3f box( vector3f( 0, 0, 0 ), vector3f( 1, 1, 1 ) );

    boundbox3f voxelBox = box.get_voxel_bounds( vector3( 0, 1, 2 ), size3( 8, 8, 8 ) );

    EXPECT_EQ( voxelBox.minimum(), vector3f( 0, 1 / 8.f, 2 / 8.f ) );
    EXPECT_EQ( voxelBox.maximum(), vector3f( 1 / 8.f, 2 / 8.f, 3 / 8.f ) );

    EXPECT_TRUE( voxelBox.is_cube() );

    // Test <count> random rays
    int count = 1000;
    for( int i = 0; i < count; ++i ) {
        ray3f ray = ray3f::from_random_towards_box( box );

        double startDist = 0, endDist = 0;
        vector3f startNormal, endNormal;
        ray.intersect_with_box( box, startDist, endDist, &startNormal, &endNormal );
        vector3f start = ray.at( startDist ), end = ray.at( endDist );

        // Both points should have landed on the box's surface
        EXPECT_LT( fabs( box.distance_function( start ) ), 0.001f );
        EXPECT_LT( fabs( box.distance_function( end ) ), 0.001f );
    }

    // Validate the 8 corners.  They should be in 'binary' order, with x changing every index increment.
    box.set( vector3f( 0 ), vector3f( 1 ) );
    for( int i = 0; i < 8; ++i ) {
        vector3f corner = box.get_corner( i );
        EXPECT_EQ( i, int( corner.x + corner.y * 2 + corner.z * 4 ) );
    }

    // test that the support function works correctly
    EXPECT_EQ( box.support( vector3f( -1.47f, -1.f, -0.2f ) ), box.get_corner( 0 ) );
    EXPECT_EQ( box.support( vector3f( 1, -0.76f, -0.2f ) ), box.get_corner( 1 ) );
    EXPECT_EQ( box.support( vector3f( -0.001f, 0.0021f, -0.2f ) ), box.get_corner( 2 ) );
    EXPECT_EQ( box.support( vector3f( 0.1f, 0.1f, -0.2f ) ), box.get_corner( 3 ) );
    EXPECT_EQ( box.support( vector3f( -2, -0.54f, 0.2f ) ), box.get_corner( 4 ) );
    EXPECT_EQ( box.support( vector3f( 1, -2434.4f, 0.2f ) ), box.get_corner( 5 ) );
    EXPECT_EQ( box.support( vector3f( -7, 0.1f, 0.2f ) ), box.get_corner( 6 ) );
    EXPECT_EQ( box.support( vector3f( 1, 0.1f, 0.2f ) ), box.get_corner( 7 ) );

    EXPECT_TRUE( boundbox3f( 0, 1, 1, 2, 0, 1 ).is_intersecting( boundbox3f( 0.9f, 2, 1, 1.5f, -1, 0 ) ) );
    EXPECT_FALSE( boundbox3f( 0, 1, 1, 2, 0, 1 ).is_intersecting( boundbox3f( 1.01f, 2, 1, 1.5f, -1, 0 ) ) );
    EXPECT_FALSE( boundbox3f( 0, 1, 1, 2, 0, 1 ).is_intersecting( boundbox3f( 0.9f, 2, 1, 1.5f, -1, -0.01f ) ) );
    EXPECT_FALSE( boundbox3f( 0, 1, 1, 2, 0, 1 ).is_intersecting( boundbox3f( 1.1f, 2, 2.01f, 2.5f, -1, -0.1f ) ) );
}

TEST( BoundBox3f, Modulus ) {
    boundbox3f theBox( vector3f( 0.f ), vector3f( 1.f ) );

    // test the basic mod vector function
    for( size_t attempts = 0; attempts < 1000; ++attempts ) {
        vector3f testPoint = theBox.mod_point( vector3f( (float)rand(), (float)rand(), (float)rand() ) );
        EXPECT_TRUE( theBox.contains( testPoint ) );
    }

    // test the periodic bounds helper function
    vector<vector3f> thePoints;

    // x
    theBox.compute_mod_boundary_points( vector3f( 0.5f, 0.f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.ysize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.5f, 1.f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.ysize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.5f, 0.f, 1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.ysize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.5f, 1.f, 1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.ysize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    // y
    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.f, 0.5f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 1.f, 0.5f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.f, 0.5f, 1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 1.f, 0.5f, 1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    // z
    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.f, 0.f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 1.f, 0.0f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.f, 1.f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 1.f, 1.f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ), sqrt( theBox.xsize() + theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    // test the 8 corners
    for( int i = 0; i < 8; ++i ) {
        vector3f testCorner = theBox.get_corner( i );
        thePoints.clear();

        theBox.compute_mod_boundary_points( testCorner, 0.01f, thePoints );

        for( int j = 0; j < 8; ++j ) {
            vector3f realCorner = theBox.get_corner( j );
            bool found = false;
            for( int p = 0; p < 8; ++p ) {
                if( thePoints[p] == realCorner ) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE( found );
        }

    } //*/

    // TEST rectangular bounds
    theBox = boundbox3f( vector3f( -1.f ), vector3f( 2.f, 5.f, 0.f ) );

    // test the basic mod vector function
    for( size_t attempts = 0; attempts < 1000; ++attempts ) {
        vector3f testPoint = theBox.mod_point( vector3f( (float)rand(), (float)rand(), (float)rand() ) );
        EXPECT_TRUE( theBox.contains( testPoint ) );
    }

    // test the periodic bounds helper function

    // x
    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.5f, -1.f, -1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    // cout << "\n P: " << thePoints[1] << " P: " << thePoints[3];
    // cout <<  "distance: " << vector3f::distance(thePoints[1], thePoints[3]) << " ysize: " << theBox.ysize() << "
    // zsize: " << theBox.zsize() << endl;

    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.ysize() * theBox.ysize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.5f, 5.f, -1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.ysize() * theBox.ysize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.5f, -1.f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.ysize() * theBox.ysize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 0.5f, 5.f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.ysize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.ysize() * theBox.ysize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    // y
    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( -1.f, 0.5f, -1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 2.f, 0.5f, -1.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( -1.f, 0.5f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 2.f, 0.5f, 0.f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.zsize() * theBox.zsize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.zsize(), 0.0001f );

    // z
    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( -1.f, -1.f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.ysize() * theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( -1.f, 5.0f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.ysize() * theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 2.f, -1.f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.ysize() * theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    thePoints.clear();
    theBox.compute_mod_boundary_points( vector3f( 2.f, 5.f, 0.5f ), 0.01f, thePoints );
    EXPECT_TRUE( thePoints.size() == 4 );
    EXPECT_NEAR( vector3f::distance( thePoints[0], thePoints[3] ), theBox.xsize(), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[1], thePoints[3] ),
                 sqrt( theBox.xsize() * theBox.xsize() + theBox.ysize() * theBox.ysize() ), 0.0001f );
    EXPECT_NEAR( vector3f::distance( thePoints[2], thePoints[3] ), theBox.ysize(), 0.0001f );

    //*
    // test the 8 corners
    for( int i = 0; i < 8; ++i ) {
        vector3f testCorner = theBox.get_corner( i );
        thePoints.clear();

        theBox.compute_mod_boundary_points( testCorner, 0.01f, thePoints );

        for( int j = 0; j < 8; ++j ) {
            vector3f realCorner = theBox.get_corner( j );
            bool found = false;
            for( int p = 0; p < 8; ++p ) {
                if( thePoints[p] == realCorner ) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE( found );
        }

    } //*/
}

namespace {
struct vector3f_triple {
    frantic::graphics::vector3f a, b, c;
    vector3f_triple( const frantic::graphics::vector3f& ai, const frantic::graphics::vector3f& bi,
                     const frantic::graphics::vector3f& ci ) {
        a = ai;
        b = bi;
        c = ci;
    }
};
} // anonymous namespace

// Tests the boundbox3f::is_intersecting_triangle function
TEST( BoundBox3f, IsIntersectingTriangle ) {
    // An arbitrary box for our tests.
    boundbox3f box( vector3f( 10, 1, 2 ), vector3f( 15, 5, 4 ) );

    // Some of the tests below were created by making a box of the above dimensions in 3ds Max, then moving around a
    // triangle into intersecting
    // and non-intersecting positions.  The following maxscript was used to print out the triangle's vertices:
    //   for i = 1 to (getnumverts $) do (v = (getVert $ i);format "vector3f(%f,%f,%f), " v.x v.y v.z ); format "\n"

    std::vector<vector3f_triple> trianglesIntersect, trianglesNotIntersect;

    // Some triangles which are exactly touching the boundary of the box
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 5, 0, 2 ), vector3f( 20, 0, 2 ), vector3f( 12, 1.1f, 2 ) ) ); // XY plane, Z minimum
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 15, 5, 4 ), vector3f( 20, 2, 4 ), vector3f( 10, 0, 4 ) ) ); // XY plane, Z maximum
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 11, 1, 2 ), vector3f( 14, 1, 3 ), vector3f( 14, 1, 2 ) ) ); // XZ plane, Y minimum
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 9, 5, 3 ), vector3f( 11, 5, 1 ), vector3f( 9, 5, 0 ) ) ); // XZ plane, Y maximum
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 10, 2, 3 ), vector3f( 10, 5, 4 ), vector3f( 10, 3, 1 ) ) ); // YZ plane, X minimum
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 15, 3, 5 ), vector3f( 15, 1, 0 ), vector3f( 15, 4, 1 ) ) ); // YZ plane, X maximum
    // Some triangles which are just outside the boundary of the box
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 5, 0, 1.999f ), vector3f( 20, 0, 1.999f ),
                                                      vector3f( 12, 1.1f, 1.999f ) ) ); // XY plane, Z minimum
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 15, 5, 4.001f ), vector3f( 20, 2, 4.001f ),
                                                      vector3f( 10, 0, 4.001f ) ) ); // XY plane, Z maximum
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 11, 0.999f, 2 ), vector3f( 14, 0.999f, 3 ),
                                                      vector3f( 14, 0.999f, 2 ) ) ); // XZ plane, Y minimum
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 9, 5.001f, 3 ), vector3f( 11, 5.001f, 1 ),
                                                      vector3f( 9, 5.001f, 0 ) ) ); // XZ plane, Y maximum
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 9.999f, 2, 3 ), vector3f( 9.999f, 5, 4 ),
                                                      vector3f( 9.999f, 3, 1 ) ) ); // YZ plane, X minimum
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 15.001f, 3, 5 ), vector3f( 15.001f, 1, 0 ),
                                                      vector3f( 15.001f, 4, 1 ) ) ); // YZ plane, X maximum
    // Triangles which are grazing the four Z edges
    trianglesIntersect.push_back( vector3f_triple( vector3f( 16, 2, 2 ), vector3f( 16, 2, 4 ), vector3f( 14, 5, 3 ) ) );
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 13, 12, 2 ), vector3f( 13, 12, 4 ), vector3f( 8, -2, 3 ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 1, 8, 2 ), vector3f( 1, 8, 4 ), vector3f( 18, -4, 3 ) ) );
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 19, 6, 2 ), vector3f( 19, 6, 4 ), vector3f( 12, -2, 3 ) ) );
    // Triangles which are just missing the four Z edges
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 16, 3, 2 ), vector3f( 16, 3, 4 ), vector3f( 14, 8, 3 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 12, 13, 2 ), vector3f( 12, 13, 4 ), vector3f( 8, -1, 3 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 1, 8, 2 ), vector3f( 1, 8, 4 ), vector3f( 18, -6, 3 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 19, 4, 2 ), vector3f( 19, 4, 4 ), vector3f( 12, -2, 3 ) ) );
    // Triangles which are grazing the four Y edges
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 8, 1, -4 ), vector3f( 8, 5, -4 ), vector3f( 12, 3, 12 ) ) );
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 19, 1, -3 ), vector3f( 19, 5, -3 ), vector3f( 12, 3, 9 ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 12, 1, 1 ), vector3f( 12, 5, 1 ), vector3f( 9, 3, 3 ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 7, 1, -5 ), vector3f( 7, 5, -5 ), vector3f( 20, 3, 7 ) ) );
    // Triangles which are just missing the four Y edges
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 8, 1, -4 ), vector3f( 8, 5, -4 ), vector3f( 12, 3, 13 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 19, 1, -3 ), vector3f( 19, 5, -3 ), vector3f( 12, 3, 10 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 12, 1, 1 ), vector3f( 12, 5, 1 ), vector3f( 7, 3, 3 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 7, 1, -5 ), vector3f( 7, 5, -5 ), vector3f( 20, 3, 6 ) ) );
    // Triangles which are grazing the four X edges
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 12, -6, -2 ), vector3f( 12, 10, 11 ), vector3f( 12, 0, 12 ) ) );
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 12, 0, -7 ), vector3f( 12, 6, 4 ), vector3f( 12, 10, -2 ) ) );
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 12, 0, 5 ), vector3f( 12, 4, -6 ), vector3f( 12, -3, 0 ) ) );
    trianglesIntersect.push_back(
        vector3f_triple( vector3f( 12, 1, 11 ), vector3f( 12, 9, -3 ), vector3f( 12, 10, 10 ) ) );
    // Triangles which are just missing the four X edges
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 12, -6, -1 ), vector3f( 12, 10, 11 ), vector3f( 12, 0, 12 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 12, 1, -7 ), vector3f( 12, 6, 4 ), vector3f( 12, 10, -2 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 12, 0, 5 ), vector3f( 12, 4, -8 ), vector3f( 12, -3, 0 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( 12, 1, 12 ), vector3f( 12, 9, -3 ), vector3f( 12, 10, 10 ) ) );
    // Triangles which are grazing the corners
    trianglesIntersect.push_back( vector3f_triple( vector3f( 8.4675f, 1.47506f, 5.91476f ),
                                                   vector3f( 12.7314f, -6.28364f, -1.58767f ),
                                                   vector3f( 18.0826f, 4.15162f, 5.83094f ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 8.4675f, 1.47506f, 3.93348f ),
                                                   vector3f( 7.10301f, -6.28364f, -1.58767f ),
                                                   vector3f( 14.9383f, 4.15162f, 7.17738f ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 6.9999f, -1.3363f, 6.50172f ),
                                                   vector3f( 9.26569f, 4.73f, -1.58767f ),
                                                   vector3f( 14.9383f, -1.29998f, 1.63947f ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 10.4688f, -6.34769f, 6.50172f ),
                                                   vector3f( 14.6013f, 2.73863f, -1.58767f ),
                                                   vector3f( 18.509f, 4.53375f, 2.43517f ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 22.4509f, 2.14957f, 5.95894f ),
                                                   vector3f( 13.4377f, 3.11324f, -1.58767f ),
                                                   vector3f( 11.1875f, 8.41777f, 2.43517f ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 17.0432f, 9.55174f, 5.95894f ),
                                                   vector3f( 13.4377f, 3.11324f, -1.58767f ),
                                                   vector3f( 2.4229f, 3.42399f, 2.43517f ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 17.0432f, 9.55174f, 5.95894f ),
                                                   vector3f( 11.2697f, 7.78706f, -2.05579f ),
                                                   vector3f( 2.4229f, -1.35744f, 6.10678f ) ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( 18.6213f, 3.679f, 5.95894f ),
                                                   vector3f( 14.6339f, 5.99776f, -2.05579f ),
                                                   vector3f( 11.3685f, 5.62751f, 6.10678f ) ) );
    // Triangles which are just missing the corners
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 8.4675f, 1.47506f, 5.91476f ),
                                                      vector3f( 12.7314f, -6.28364f, -1.58767f ),
                                                      vector3f( 18.0826f, 4.15162f, 6.3791f ) ) );
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 8.4675f, 1.47506f, 4.08239f ),
                                                      vector3f( 7.10301f, -6.28364f, -1.58767f ),
                                                      vector3f( 14.9383f, 4.15162f, 7.17738f ) ) );
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 6.9999f, -1.3363f, 6.50172f ),
                                                      vector3f( 9.26569f, 4.73f, -1.58767f ),
                                                      vector3f( 14.9383f, -1.75229f, 1.63947f ) ) );
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 10.4688f, -6.34769f, 5.95894f ),
                                                      vector3f( 14.6013f, 2.73863f, -1.58767f ),
                                                      vector3f( 18.509f, 4.53375f, 2.43517f ) ) );
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 22.4509f, 2.14957f, 5.95894f ),
                                                      vector3f( 13.4377f, 3.11324f, -1.58767f ),
                                                      vector3f( 11.4385f, 8.6271f, 2.43517f ) ) );
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 17.0432f, 9.55174f, 5.95894f ),
                                                      vector3f( 13.4377f, 3.11324f, -2.05579f ),
                                                      vector3f( 2.4229f, 3.42399f, 2.43517f ) ) );
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 16.7589f, 9.55174f, 5.95894f ),
                                                      vector3f( 11.2697f, 7.78706f, -2.05579f ),
                                                      vector3f( 2.4229f, -1.35744f, 6.10678f ) ) );
    trianglesNotIntersect.push_back( vector3f_triple( vector3f( 18.6213f, 3.679f, 5.95894f ),
                                                      vector3f( 14.6339f, 5.99776f, -2.05579f ),
                                                      vector3f( 11.3685f, 5.62751f, 7.43034f ) ) );

    // Test all the triangles that are intersecting
    for( unsigned i = 0; i < trianglesIntersect.size(); ++i ) {
        bool isIntersecting =
            box.is_intersecting_triangle( trianglesIntersect[i].a, trianglesIntersect[i].b, trianglesIntersect[i].c );
        if( !isIntersecting ) {
            std::cout << "is_intersecting_triangle " << trianglesIntersect[i].a << ", " << trianglesIntersect[i].b
                      << ", " << trianglesIntersect[i].c << std::endl;
            EXPECT_TRUE( isIntersecting );
        }
        boundbox3f intersection;
        isIntersecting = box.intersect_with_triangle( trianglesIntersect[i].a, trianglesIntersect[i].b,
                                                      trianglesIntersect[i].c, intersection );
        if( !isIntersecting ) {
            std::cout << "intersect_with_triangle " << trianglesIntersect[i].a << ", " << trianglesIntersect[i].b
                      << ", " << trianglesIntersect[i].c << std::endl;
            EXPECT_TRUE( isIntersecting );
        }
    }
    // Test all the triangles that are not intersecting
    for( unsigned i = 0; i < trianglesNotIntersect.size(); ++i ) {
        bool isIntersecting = box.is_intersecting_triangle( trianglesNotIntersect[i].a, trianglesNotIntersect[i].b,
                                                            trianglesNotIntersect[i].c );
        if( isIntersecting ) {
            std::cout << "is_intersecting_triangle " << trianglesNotIntersect[i].a << ", " << trianglesNotIntersect[i].b
                      << ", " << trianglesNotIntersect[i].c << std::endl;
            EXPECT_FALSE( isIntersecting );
        }
        boundbox3f intersection;
        isIntersecting = box.intersect_with_triangle( trianglesNotIntersect[i].a, trianglesNotIntersect[i].b,
                                                      trianglesNotIntersect[i].c, intersection );
        if( isIntersecting ) {
            std::cout << "intersect_with_triangle " << trianglesNotIntersect[i].a << ", " << trianglesNotIntersect[i].b
                      << ", " << trianglesNotIntersect[i].c << std::endl;
            EXPECT_FALSE( isIntersecting );
        }
    }

    // Test with a box-plane
    trianglesIntersect.clear();
    trianglesNotIntersect.clear();
    box.set( vector3f( -2, 5, 8 ), vector3f( 4, 12, 8 ) );
    trianglesIntersect.push_back( vector3f_triple( vector3f( -2, 5, 8 ), vector3f( 4, 12, 8 ), vector3f( 4, 5, 8 ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( -2, 5, 8.00001f ), vector3f( 4, 12, 8.00001f ), vector3f( 4, 5, 8.00001f ) ) );
    trianglesNotIntersect.push_back(
        vector3f_triple( vector3f( -2, 5, 7.99999f ), vector3f( 4, 12, 7.99999f ), vector3f( 4, 5, 7.99999f ) ) );

    // Test all the triangles that are intersecting
    for( unsigned i = 0; i < trianglesIntersect.size(); ++i ) {
        bool isIntersecting =
            box.is_intersecting_triangle( trianglesIntersect[i].a, trianglesIntersect[i].b, trianglesIntersect[i].c );
        if( !isIntersecting ) {
            std::cout << "is_intersecting_triangle " << trianglesIntersect[i].a << ", " << trianglesIntersect[i].b
                      << ", " << trianglesIntersect[i].c << std::endl;
            EXPECT_TRUE( isIntersecting );
        }
        boundbox3f intersection;
        isIntersecting = box.intersect_with_triangle( trianglesIntersect[i].a, trianglesIntersect[i].b,
                                                      trianglesIntersect[i].c, intersection );
        if( !isIntersecting ) {
            std::cout << "intersect_with_triangle " << trianglesIntersect[i].a << ", " << trianglesIntersect[i].b
                      << ", " << trianglesIntersect[i].c << std::endl;
            EXPECT_TRUE( isIntersecting );
        }
    }
    // Test all the triangles that are not intersecting
    for( unsigned i = 0; i < trianglesNotIntersect.size(); ++i ) {
        bool isIntersecting = box.is_intersecting_triangle( trianglesNotIntersect[i].a, trianglesNotIntersect[i].b,
                                                            trianglesNotIntersect[i].c );
        if( isIntersecting ) {
            std::cout << "is_intersecting_triangle " << trianglesNotIntersect[i].a << ", " << trianglesNotIntersect[i].b
                      << ", " << trianglesNotIntersect[i].c << std::endl;
            EXPECT_FALSE( isIntersecting );
        }
        boundbox3f intersection;
        isIntersecting = box.intersect_with_triangle( trianglesNotIntersect[i].a, trianglesNotIntersect[i].b,
                                                      trianglesNotIntersect[i].c, intersection );
        if( isIntersecting ) {
            std::cout << "intersect_with_triangle " << trianglesNotIntersect[i].a << ", " << trianglesNotIntersect[i].b
                      << ", " << trianglesNotIntersect[i].c << std::endl;
            EXPECT_FALSE( isIntersecting );
        }
    }
}

TEST( BoundBox3f, IntersectWithTriangle ) {
    // An arbitrary box for our tests.
    boundbox3f box( vector3f( 10, 1, 2 ), vector3f( 15, 5, 4 ) );

    boundbox3f intersectedBox;

    // Intersect various triangles, and verify that the resulting intersection boxes are correct
    EXPECT_TRUE( box.intersect_with_triangle( vector3f( 14, 2, 4 ), vector3f( 18, 0, 4 ), vector3f( 16, -1, 4 ),
                                              intersectedBox ) );
    EXPECT_TRUE( intersectedBox.minimum() == vector3f( 14, 1, 4 ) );
    EXPECT_TRUE( intersectedBox.maximum() == vector3f( 15, 2, 4 ) );
    EXPECT_TRUE( box.intersect_with_triangle( vector3f( 14, 2, 3 ), vector3f( 18, 0, 3 ), vector3f( 16, -1, 3 ),
                                              intersectedBox ) );
    EXPECT_TRUE( intersectedBox.minimum() == vector3f( 14, 1, 3 ) );
    EXPECT_TRUE( intersectedBox.maximum() == vector3f( 15, 2, 3 ) );
    EXPECT_TRUE( box.intersect_with_triangle( vector3f( 14, 2, 2 ), vector3f( 18, 0, 2 ), vector3f( 16, -1, 2 ),
                                              intersectedBox ) );
    EXPECT_TRUE( intersectedBox.minimum() == vector3f( 14, 1, 2 ) );
    EXPECT_TRUE( intersectedBox.maximum() == vector3f( 15, 2, 2 ) );

    EXPECT_TRUE( box.intersect_with_triangle( vector3f( 12, 3, 3 ), vector3f( 12, 5, 5 ), vector3f( 12, 1, 5 ),
                                              intersectedBox ) );
    EXPECT_TRUE( intersectedBox.minimum() == vector3f( 12, 2, 3 ) );
    EXPECT_TRUE( intersectedBox.maximum() == vector3f( 12, 4, 4 ) );

    EXPECT_TRUE( box.intersect_with_triangle( vector3f( 12, 2, 2.5f ), vector3f( 11, 3, 3 ), vector3f( 13, 4, 3 ),
                                              intersectedBox ) );
    EXPECT_TRUE( intersectedBox.minimum() == vector3f( 11, 2, 2.5f ) );
    EXPECT_TRUE( intersectedBox.maximum() == vector3f( 13, 4, 3 ) );
}

TEST( BoundBox3f, GetBoundBoxPlanes ) {
    using namespace frantic::graphics;

    vector3f coords[2] = {
        vector3f( 3, 4, 5 ),
        vector3f( 9, 16, 25 ),
    };
    boundbox3f testBox( coords[0], coords[1] );

    plane3f planes[6];
    get_bounding_box_planes( testBox, planes );

    for( size_t z = 0; z < 2; ++z ) {
        for( size_t y = 0; y < 2; ++y ) {
            for( size_t x = 0; x < 2; ++x ) {
                vector3f position( coords[x].x, coords[y].y, coords[z].z );

                if( x == 0 ) {
                    EXPECT_NEAR( 0.0f, planes[cube_face::CF_X_NEG].get_signed_distance_to_plane( position ), 0.0001f );
                    EXPECT_LT( 0.0001f,
                               std::abs( planes[cube_face::CF_X_POS].get_signed_distance_to_plane( position ) ) );
                } else {
                    EXPECT_NEAR( 0.0f, planes[cube_face::CF_X_POS].get_signed_distance_to_plane( position ), 0.0001f );
                    EXPECT_LT( 0.0001f,
                               std::abs( planes[cube_face::CF_X_NEG].get_signed_distance_to_plane( position ) ) );
                }

                if( y == 0 ) {
                    EXPECT_NEAR( 0.0f, planes[cube_face::CF_Y_NEG].get_signed_distance_to_plane( position ), 0.0001f );
                    EXPECT_LT( 0.0001f,
                               std::abs( planes[cube_face::CF_Y_POS].get_signed_distance_to_plane( position ) ) );
                } else {
                    EXPECT_NEAR( 0.0f, planes[cube_face::CF_Y_POS].get_signed_distance_to_plane( position ), 0.0001f );
                    EXPECT_LT( 0.0001f,
                               std::abs( planes[cube_face::CF_Y_NEG].get_signed_distance_to_plane( position ) ) );
                }

                if( z == 0 ) {
                    EXPECT_NEAR( 0.0f, planes[cube_face::CF_Z_NEG].get_signed_distance_to_plane( position ), 0.0001f );
                    EXPECT_LT( 0.0001f,
                               std::abs( planes[cube_face::CF_Z_POS].get_signed_distance_to_plane( position ) ) );
                } else {
                    EXPECT_NEAR( 0.0f, planes[cube_face::CF_Z_POS].get_signed_distance_to_plane( position ), 0.0001f );
                    EXPECT_LT( 0.0001f,
                               std::abs( planes[cube_face::CF_Z_NEG].get_signed_distance_to_plane( position ) ) );
                }
            }
        }
    }
}

TEST( BoundBox3f, TestInfiniteBox ) {
    using namespace frantic::graphics;

    boundbox3f infinityBox = boundbox3f::infinite();

    EXPECT_TRUE( infinityBox.is_infinite() );
    EXPECT_TRUE( infinityBox.contains( vector3f( std::numeric_limits<float>::max() ) ) );
    EXPECT_TRUE( infinityBox.contains( vector3f( std::numeric_limits<float>::min() ) ) );
}
