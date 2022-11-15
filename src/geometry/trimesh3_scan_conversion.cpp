// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_scan_conversion.hpp>
#include <frantic/math/utils.hpp>

using namespace std;
using namespace frantic::math;

namespace frantic {
namespace geometry {

void trimesh3_scan_convert( const graphics::transform4f& xform, const trimesh3& mesh, graphics2d::size2 dimensions,
                            std::vector<std::vector<scan_conversion_intersection>>& outIntersectionDepths ) {
    if( (int)outIntersectionDepths.size() != dimensions.get_area() )
        throw runtime_error(
            "trimesh3_scan_convert: The output intersection depth array provided had the incorrect size." );

    // First transform all the vertices into the destination space
    std::vector<vector3f> transformedVertices;
    transformedVertices.reserve( mesh.vertex_count() );
    for( unsigned i = 0; i < mesh.vertex_count(); ++i ) {
        transformedVertices.push_back( xform * mesh.get_vertex( i ) );
    }

    scan_conversion_intersection sci;
    vector<pair<int, int>> scanlines( dimensions.ysize );
    int yStart, yEnd;
    vector3f a, b, c;
    for( sci.faceIndex = 0; sci.faceIndex < (int)mesh.face_count(); ++sci.faceIndex ) {
        vector3 face = mesh.get_face( sci.faceIndex );
        a = transformedVertices[face.x];
        b = transformedVertices[face.y];
        c = transformedVertices[face.z];

        // If there's a weird vertex, skip the drawing of this face.
        if( !a.is_finite() || !b.is_finite() || !c.is_finite() )
            continue;

        sci.normalFacingZPositive = ( ( c.x - b.x ) * ( a.y - b.y ) - ( c.y - b.y ) * ( a.x - b.x ) ) >= 0;

        //		cout << "Face " << sci.faceIndex << endl;
        //		cout << "normal facing Z pos: " << sci.normalFacingZPositive << endl;
        //		cout << "normal " << vector3f::from_triangle_normal(a,b,c) << endl;
        //		cout << a << ", " << b << ", " << c << endl;

        //////////////
        // Sort the vertices in ascending Y
        //////////////

        int vertPerm[3] = { 0, 1, 2 };

        if( b.y < a.y && b.y <= c.y ) {
            swap( a, b );
            swap( vertPerm[0], vertPerm[1] );
        } else if( c.y < a.y && c.y <= b.y ) {
            swap( a, c );
            swap( vertPerm[0], vertPerm[2] );
        }
        if( c.y < b.y ) {
            swap( b, c );
            swap( vertPerm[1], vertPerm[2] );
        }

        // If the whole triangle is outside of the bounds, skip to the next one
        if( c.y < 0.5f || a.y > dimensions.ysize - 0.5f )
            continue;
        if( ( a.x < 0.5f && b.x < 0.5f && c.x < 0.5f ) ||
            ( a.x > dimensions.xsize - 0.5f && b.x > dimensions.xsize - 0.5f && c.x > dimensions.xsize - 0.5f ) )
            continue;

        // Get the y range of values
        yStart = (int)ceilf( max( 0.f, a.y - 0.5f ) );
        yEnd = (int)floorf( min( dimensions.ysize - 1.f, c.y - 0.5f ) );

        //		cout << " y range: " << yStart << ", " << yEnd << endl;

        //////////////
        // Convert the triangle into a set of X-intervals in the scanlines array
        //////////////

        // Determine whether the middle vertex is to the right or to the left of the line
        // connecting the top and bottom vertices.
        if( ( ( c.x - b.x ) * ( a.y - b.y ) - ( c.y - b.y ) * ( a.x - b.x ) ) < 0 ) {
            // Put the two segments in the first variable
            int yMiddle = (int)floorf( b.y );
            if( yMiddle < yStart )
                yMiddle = yStart;
            if( yMiddle > yEnd )
                yMiddle = yEnd;

            //			cout << " y middle: " << yMiddle << endl;

            float lineInverseSlope;
            if( b.y > a.y ) {
                lineInverseSlope = ( b.x - a.x ) / ( b.y - a.y );
                for( int y = yStart; y <= yMiddle; ++y ) {
                    float x = a.x + ( y + 0.5f - a.y ) * lineInverseSlope - 0.5f;
                    scanlines[y].first = (int)ceilf( x );
                    //					cout << "f: (" << scanlines[y].first << ", " << y << ")" <<
                    //endl;
                }
            } else {
                // Have to initialize scanlines[yMiddle].first to a small value, so the max operation below picks the
                // right one
                scanlines[yMiddle].first = 0;
            }
            if( c.y > b.y ) {
                lineInverseSlope = ( c.x - b.x ) / ( c.y - b.y );
                // Take the maximum of the value from the two lines for the middle y coordinate
                int x = (int)ceilf( b.x + ( yMiddle + 0.5f - b.y ) * lineInverseSlope - 0.5f );
                if( x > scanlines[yMiddle].first )
                    scanlines[yMiddle].first = x;

                //				cout << "f: (" << scanlines[yMiddle].first << ", " << yMiddle << ")" <<
                //endl;

                for( int y = yMiddle + 1; y <= yEnd; ++y ) {
                    x = (int)ceilf( b.x + ( y + 0.5f - b.y ) * lineInverseSlope - 0.5f );
                    scanlines[y].first = x;
                    //					cout << "f: (" << scanlines[y].first << ", " << y << ")" <<
                    //endl;
                }
            }

            // Put the third segment in the second variable
            if( c.y > a.y ) {
                lineInverseSlope = ( c.x - a.x ) / ( c.y - a.y );
                for( int y = yStart; y <= yEnd; ++y ) {
                    float x = a.x + ( y + 0.5f - a.y ) * lineInverseSlope - 0.5f;
                    scanlines[y].second = (int)floorf( x );
                    //					cout << "s: (" << scanlines[y].second << ", " << y << ")" <<
                    //endl;
                }
            }

        } else {
            // Put the two segments in the second variable
            int yMiddle = (int)floorf( b.y );
            if( yMiddle < yStart )
                yMiddle = yStart;
            if( yMiddle > yEnd )
                yMiddle = yEnd;

            //			cout << " y middle: " << yMiddle << endl;

            float lineInverseSlope;
            if( b.y > a.y ) {
                lineInverseSlope = ( b.x - a.x ) / ( b.y - a.y );
                for( int y = yStart; y <= yMiddle; ++y ) {
                    float x = a.x + ( y + 0.5f - a.y ) * lineInverseSlope - 0.5f;
                    scanlines[y].second = (int)floorf( x );
                    //					cout << "s: (" << scanlines[y].second << ", " << y << ")" <<
                    //endl;
                }
            } else {
                // Have to initialize scanlines[yMiddle].second to a large value, so the min operation below picks the
                // right one
                scanlines[yMiddle].second = dimensions.xsize;
                //				cout << "s skip" << endl;
            }
            if( c.y > b.y ) {
                lineInverseSlope = ( c.x - b.x ) / ( c.y - b.y );
                // Take the minimum of the value from the two lines for the middle y coordinate
                int x = (int)floorf( b.x + ( yMiddle + 0.5f - b.y ) * lineInverseSlope - 0.5f );
                if( x < scanlines[yMiddle].second )
                    scanlines[yMiddle].second = x;

                //				cout << "s: (" << scanlines[yMiddle].second << ", " << yMiddle << ")" <<
                //endl;

                for( int y = yMiddle + 1; y <= yEnd; ++y ) {
                    x = (int)floorf( b.x + ( y + 0.5f - b.y ) * lineInverseSlope - 0.5f );
                    scanlines[y].second = x;
                    //					cout << "s: (" << scanlines[y].second << ", " << y << ")" <<
                    //endl;
                }
            }

            // Put the third segment in the first variable
            if( c.y > a.y ) {
                lineInverseSlope = ( c.x - a.x ) / ( c.y - a.y );
                for( int y = yStart; y <= yEnd; ++y ) {
                    float x = a.x + ( y + 0.5f - a.y ) * lineInverseSlope - 0.5f;
                    scanlines[y].first = (int)ceilf( x );
                    //					cout << "f: (" << scanlines[y].first << ", " << y << ")" <<
                    //endl;
                }
            }
        }

        //////////////
        // Convert each X-interval into samples
        //////////////

        // Get the linear equation which maps (x,y) -> z for this plane
        // Note that we're baking the +0.5 for converting from integer indexes to
        // float sample position into this linear equation.
        vector3f N = vector3f::cross( b - a, c - a );
        N /= N.z;
        float planeA = -N.x, planeB = -N.y, planeC = vector3f::dot( N, a ) - 0.5f * ( N.x + N.y );

        // Compute the barycentric coordinate matrix
        float baryInvMatrix[4], baryMatrix[4];
        baryInvMatrix[0] = a.x - c.x;
        baryInvMatrix[1] = a.y - c.y;
        baryInvMatrix[2] = b.x - c.x;
        baryInvMatrix[3] = b.y - c.y;
        float baryInvMatrixDet = baryInvMatrix[0] * baryInvMatrix[3] - baryInvMatrix[2] * baryInvMatrix[1];
        if( baryInvMatrixDet != 0 ) {
            // Invert the matrix, to get the mapping from (x,y) to barycentric coordinates
            baryMatrix[0] = baryInvMatrix[3] / baryInvMatrixDet;
            baryMatrix[1] = -baryInvMatrix[1] / baryInvMatrixDet;
            baryMatrix[2] = -baryInvMatrix[2] / baryInvMatrixDet;
            baryMatrix[3] = baryInvMatrix[0] / baryInvMatrixDet;

            for( int y = yStart; y <= yEnd; ++y ) {
                int xStart = scanlines[y].first, xEnd = scanlines[y].second;
                if( xStart < 0 )
                    xStart = 0;
                if( xEnd >= dimensions.xsize )
                    xEnd = dimensions.xsize - 1;
                for( int x = xStart; x <= xEnd; ++x ) {
                    // Compute the Z of the intersection
                    sci.z = planeA * x + planeB * y + planeC;
                    // Somehow, some non-finite values were slipping into here.  Cull out any intersections that gave a
                    // non-finite value.
                    if( is_finite( sci.z ) ) {
                        // Compute the barycentric coordinates
                        float fx = x - c.x + 0.5f, fy = y - c.y + 0.5f;
                        float baryX = baryMatrix[0] * fx + baryMatrix[2] * fy,
                              baryY = baryMatrix[1] * fx + baryMatrix[3] * fy;
                        sci.barycentricCoord[vertPerm[0]] = baryX;
                        sci.barycentricCoord[vertPerm[1]] = baryY;
                        sci.barycentricCoord[vertPerm[2]] = 1 - baryX - baryY;
                        // Save the intersection
                        outIntersectionDepths[dimensions.get_index( x, y )].push_back( sci );
                    }
                }
            }
        }
    }

    // Finally, sort all the intersections into ascending order
    // TODO: This part is easy to multithread
    for( vector<vector<scan_conversion_intersection>>::iterator i = outIntersectionDepths.begin(),
                                                                iterEnd = outIntersectionDepths.end();
         i != iterEnd; ++i )
        sort( i->begin(), i->end() );

    /*
    cout << "\n";
    for( int y = 0; y < dimensions.ysize; ++y ) {
      for( int x = 0; x < dimensions.xsize; ++x ) {
        vector<scan_conversion_intersection>& sci = outIntersectionDepths[dimensions.get_index(x,y)];
        //cout << samples[dim.get_index(x,y)].size();
        if( sci.empty() )
          cout << "XX";
        else
          //cout << sci[0].barycentricCoord.get_largest_axis() << " ";
          //cout << (int)sci.back().z << "";
          //cout << sci.size() << "";
          cout << sci.front().normalFacingZPositive << sci.back().normalFacingZPositive;
      }
      cout << "\n";
    }
    cout << endl;
    */
}

} // namespace geometry
} // namespace frantic
