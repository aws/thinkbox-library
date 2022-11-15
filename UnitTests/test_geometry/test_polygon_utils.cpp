// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/geometry/polygon_utils.hpp>

#include <frantic/misc/utility.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/random.hpp>

#include <cmath>
#include <map>
#include <utility>
#include <vector>

TEST( PolygonUtils, PolygonArea ) {
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;

    // clang-format off
	vector3f simplePolygon[4] = {
		vector3f( -1.0f, -1.0f, -1.0f ),
		vector3f(  1.0f, -1.0f,  1.0f ),
		vector3f(  1.0f,  1.0f,  1.0f ),
		vector3f( -1.0f,  1.0f, -1.0f ),
	};
    // clang-format on

    const double area = polygon3_area( simplePolygon, simplePolygon + 4 );
    const double expectedArea = vector3f::distance( simplePolygon[0], simplePolygon[1] ) *
                                vector3f::distance( simplePolygon[1], simplePolygon[2] );

    EXPECT_NEAR( expectedArea, area, 0.0001f );
}

TEST( PolygonUtils, SimplePointInPolygonInt ) {
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;

    std::vector<vector2> polygon;
    polygon.push_back( vector2( 0, 0 ) );
    polygon.push_back( vector2( 19, 30 ) );
    polygon.push_back( vector2( 40, 6 ) );
    polygon.push_back( vector2( 31, -31 ) );

    // some easy interior points
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 20, 0 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 10, 10 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 30, -10 ) ) );

    // some easy exterior points
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 45, 3 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 20, -40 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -10, -5 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 10, 25 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 30, 24 ) ) );

    // some edge interior points
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 1, 0 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 19, 28 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 30, -28 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 39, 6 ) ) );

    // some edge exterior points
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -1, 0 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 41, 6 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 20, 30 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 18, 30 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 32, -31 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 30, -31 ) ) );
}

TEST( PolygonUtils, ComplexPointInPolygonInt ) {
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;

    std::vector<vector2> polygon;
    polygon.push_back( vector2( 0, 0 ) );
    polygon.push_back( vector2( -10, 20 ) );
    polygon.push_back( vector2( 10, 10 ) );
    polygon.push_back( vector2( 25, 20 ) );
    polygon.push_back( vector2( 14, -2 ) );
    polygon.push_back( vector2( 20, -15 ) );
    polygon.push_back( vector2( 12, -8 ) );
    polygon.push_back( vector2( -12, -21 ) );

    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 7, 0 ) ) );

    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 1, 0 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -1, 0 ) ) );

    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -9, 20 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -11, 20 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -10, 21 ) ) );

    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 10, 11 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 10, 9 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 9, 10 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 11, 10 ) ) );

    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 26, 20 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 24, 20 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 25, 21 ) ) );

    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 14, -3 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 15, -2 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 13, -2 ) ) );

    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 20, -16 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 21, -15 ) ) );

    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( 12, -9 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 12, -7 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 13, -8 ) ) );
    EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( 11, -8 ) ) );

    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -12, -22 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -12, -20 ) ) );
    EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( -13, -21 ) ) );
}

TEST( PolygonUtils, EdgePointInPolygonInt ) {
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;

    std::vector<vector2> polygon;
    polygon.push_back( vector2( 0, 0 ) );
    polygon.push_back( vector2( 0, 10 ) );
    polygon.push_back( vector2( 10, 10 ) );
    polygon.push_back( vector2( 10, 0 ) );

    for( int i = polygon[0].y; i < polygon[1].y; ++i ) {
        EXPECT_EQ( boundary_relation::boundary, point_in_polygon( polygon, vector2( polygon[0].x, i ) ) );
        EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( polygon[0].x - 1, i ) ) );
        if( i > polygon[0].y ) {
            EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( polygon[0].x + 1, i ) ) );
        }
    }

    for( int i = polygon[1].x; i < polygon[2].x; ++i ) {
        EXPECT_EQ( boundary_relation::boundary, point_in_polygon( polygon, vector2( i, polygon[1].y ) ) );
        EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( i, polygon[1].y + 1 ) ) );
        if( i > polygon[1].x ) {
            EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( i, polygon[1].y - 1 ) ) );
        }
    }

    for( int i = polygon[2].y; i > polygon[3].y; --i ) {
        EXPECT_EQ( boundary_relation::boundary, point_in_polygon( polygon, vector2( polygon[2].x, i ) ) );
        EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( polygon[2].x + 1, i ) ) );
        if( i < polygon[2].y ) {
            EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( polygon[2].x - 1, i ) ) );
        }
    }

    for( int i = polygon[3].x; i > polygon[0].x; --i ) {
        EXPECT_EQ( boundary_relation::boundary, point_in_polygon( polygon, vector2( i, polygon[3].y ) ) );
        EXPECT_EQ( boundary_relation::outside, point_in_polygon( polygon, vector2( i, polygon[3].y - 1 ) ) );
        if( i < polygon[3].x ) {
            EXPECT_EQ( boundary_relation::inside, point_in_polygon( polygon, vector2( i, polygon[3].y + 1 ) ) );
        }
    }
}

TEST( PolygonUtils, MinimumEnclosingRectangle ) {
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;

    std::srand( 20495346 );

    std::vector<vector2f> points;

    const size_t rectY = 10;
    const size_t rectX = 20;

    for( size_t y = 0; y < rectY; ++y ) {
        for( size_t x = 0; x < rectX; ++x ) {
            points.push_back( vector2f( float( x ), float( y ) ) );
        }
    }

    float theta = float( ( (float)rand() / RAND_MAX ) * ( 2.0f * M_PI ) );
    float xOff = float( ( (float)rand() / RAND_MAX ) * 100.0f );
    float yOff = float( ( (float)rand() / RAND_MAX ) * 100.0f );

    for( size_t i = 0; i < points.size(); ++i ) {
        points[i] = vector2f( cosf( theta ) * points[i].x + xOff, sinf( theta ) * points[i].y + yOff );
    }

    std::vector<vector2f> outRect;

    minimum_enclosing_rectange( points, outRect );

    EXPECT_EQ( 4, outRect.size() );

    bool found[4] = { false, false, false, false };

    size_t indices[4] = { 0, rectX - 1, ( rectY - 1 ) * rectX, ( rectX * rectY ) - 1 };

    for( size_t i = 0; i < outRect.size(); ++i ) {
        for( size_t j = 0; j < 4; ++j ) {
            if( vector2f::distance( outRect[i], points[indices[j]] ) < 1e-4f ) {
                EXPECT_TRUE( !found[j] )
                    << "More than one point on the enclosing rectangle matches a single expected result.";
                found[j] = true;
            }
        }
    }

    for( size_t j = 0; j < 4; ++j ) {
        EXPECT_TRUE( found[j] ) << "Not all expected points were found in the enclosing rectangle.";
    }
}

void check_triangluation( size_t n, frantic::graphics::vector3* triangles ) {
    using namespace frantic;

    typedef std::pair<size_t, size_t> edge_index_t;

    std::map<edge_index_t, size_t> edgeCounts;

    for( size_t i = 0; i < n; ++i ) {
        edgeCounts[make_sorted_pair( i, ( i + 1 ) % n )] = 1;
    }

    for( size_t i = 0; i < n - 2; ++i ) {
        for( size_t j = 0; j < 3; ++j ) {
            edge_index_t e = make_sorted_pair( triangles[i][j], triangles[i][( j + 1 ) % 3] );
            EXPECT_NE( e.first, e.second );

            std::map<edge_index_t, size_t>::iterator it = edgeCounts.find( e );

            if( it == edgeCounts.end() ) {
                edgeCounts[e] = 1;
            } else {
                EXPECT_EQ( 1, it->second ) << "Found more than two triangles sharing an edge.";
                ++it->second;
            }
        }
    }

    for( std::map<edge_index_t, size_t>::iterator it = edgeCounts.begin(); it != edgeCounts.end(); ++it ) {
        EXPECT_EQ( it->second, 2 )
            << "An edge of the triangulation either did not have exactly two internal edges, or did "
               "not have exactly one boundary edge.";
    }
}

TEST( PolygonUtils, Triangulate ) {
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;
    using namespace frantic::math;

    const size_t numPolygons = 3;

    std::vector<vector2f> polygons[numPolygons];

    polygons[0].push_back( vector2f( 0.0f, 0.0f ) );
    polygons[0].push_back( vector2f( 1.0f, 0.0f ) );
    polygons[0].push_back( vector2f( 1.0f, 1.0f ) );
    polygons[0].push_back( vector2f( 0.0f, 1.0f ) );

    polygons[1].push_back( vector2f( -1.0f, 1.0f ) );
    polygons[1].push_back( vector2f( 0.0f, 1.0f ) );
    polygons[1].push_back( vector2f( 1.0f, 1.0f ) );
    polygons[1].push_back( vector2f( 0.0f, 0.5f ) );

    polygons[2].push_back( vector2f( 0.0f, 0.0f ) );
    polygons[2].push_back( vector2f( 1.0f, 0.75f ) );
    polygons[2].push_back( vector2f( 1.0f, 0.0f ) );
    polygons[2].push_back( vector2f( 2.0f, 1.0f ) );
    polygons[2].push_back( vector2f( 1.25f, 0.25f ) );
    polygons[2].push_back( vector2f( 1.25f, 1.0f ) );
    polygons[2].push_back( vector2f( 0.25f, 0.25f ) );
    polygons[2].push_back( vector2f( 0.0f, 1.0f ) );

    for( size_t i = 0; i < numPolygons; ++i ) {
        std::vector<vector3> outTriangles;
        triangulate_polygon( &polygons[i][0], polygons[i].size(), std::back_inserter( outTriangles ) );
        check_triangluation( polygons[i].size(), &outTriangles[0] );
    }

    const size_t numPlanes = 5;
    plane3f planes[numPlanes] = {
        plane3f( 1.0f, 0.0f, 0.0f, 4.0f ),
        plane3f( 0.0f, 1.0f, 0.0f, -3.0f ),
        plane3f( 0.0f, 0.0f, 1.0f, 22.0f ),
        plane3f( vector3f( 3.0f, 4.0f, 5.0f ).to_normalized(), 19.0f ),
        plane3f( vector3f( 0.0f, 2.0f, 1.0f ).to_normalized(), -53.0f ),
    };

    for( size_t i = 0; i < numPlanes; ++i ) {
        transform4f transformTo3dPlane = planes[i].get_planar_transform().to_inverse();

        for( size_t j = 0; j < numPolygons; ++j ) {
            std::vector<vector3f> polygon3;

            for( size_t k = 0; k < polygons[j].size(); ++k ) {
                polygon3.push_back( vector3f( polygons[j][k].x, polygons[j][k].y, 0.0f ) );
            }

            std::vector<vector3> outTriangles;
            triangulate_polygon( &polygon3[0], polygon3.size(), std::back_inserter( outTriangles ) );
            check_triangluation( polygon3.size(), &outTriangles[0] );
        }
    }
}

/**
 * These are some example polygons that are non-coplanar, and thus triangulation fails with them
 * This is just to check that our detection of such cases is working
 */
TEST( PolygonUtils, OffPlaneTriangulation ) {
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;
    using namespace frantic::math;

    // clang-format off
	std::vector<vector3f> originalPoly1 = boost::assign::list_of
		( vector3f( -59.0461769f, 16.6899071f, -10.2295694f ) )
		( vector3f( -59.0450516f, 16.6974659f, -10.2161627f ) )
		( vector3f( -59.0402374f, 16.7018642f, -10.2243195f ) )
		( vector3f( -59.0358849f, 16.6921577f, -10.2237978f ) );
    // clang-format on

    std::vector<vector3> outTriangles1;
    const char* outError1 = NULL;
    EXPECT_FALSE( triangulate_polygon_nothrow( &originalPoly1[0], originalPoly1.size(),
                                               std::back_inserter( outTriangles1 ), outError1 ) );

    // clang-format off
	std::vector<vector3f> originalPoly2 = boost::assign::list_of
		( vector3f( -59.4944038f, 16.7060242f, -10.2015619f ) )
		( vector3f( -59.5037155f, 16.7089767f, -10.1945639f ) )
		( vector3f( -59.4987755f, 16.7183495f, -10.1956663f ) )
		( vector3f( -59.4933128f, 16.7135830f, -10.1881523f ) );
    // clang-format on

    std::vector<vector3> outTriangles2;
    const char* outError2 = NULL;
    EXPECT_FALSE( triangulate_polygon_nothrow( &originalPoly1[0], originalPoly1.size(),
                                               std::back_inserter( outTriangles2 ), outError2 ) );
}
