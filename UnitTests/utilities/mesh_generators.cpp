// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "utilities/mesh_generators.hpp"

#include <frantic/geometry/mesh_conversion.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>

#include <frantic/graphics/graphics_utils.hpp>
#include <frantic/graphics/spherical_coords.hpp>
#include <frantic/graphics/vector3.hpp>

#include <frantic/misc/indexer.hpp>

#include <boost/random.hpp>

#include <cmath>

static void add_trimesh_quad( frantic::geometry::trimesh3& mesh, size_t x0, size_t x1, size_t x2, size_t x3 ) {
    mesh.add_face( (int)x0, (int)x1, (int)x2 );
    mesh.add_face( (int)x0, (int)x2, (int)x3 );
}

void make_single_triangle_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0, 0.0, 0.0 );
    mesh.add_vertex( 0.0, 1.0, 0.0 );
    mesh.add_vertex( 1.0, 1.0, 0.0 );

    mesh.add_face( 0, 1, 2 );
}

void make_quad_trimesh( frantic::geometry::trimesh3& outMesh ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    outMesh.clear();

    outMesh.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    outMesh.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
    outMesh.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
    outMesh.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

    outMesh.add_face( vector3( 0, 1, 2 ) );
    outMesh.add_face( vector3( 0, 2, 3 ) );
}

void make_regular_tetrahedron( frantic::geometry::trimesh3& mesh ) {
    const float sqrt2inv = 1.0f / std::sqrt( 2.0f );

    mesh.clear();

    mesh.add_vertex( 1.0f, 0.0f, -sqrt2inv );
    mesh.add_vertex( -1.0f, 0.0f, -sqrt2inv );
    mesh.add_vertex( 0.0f, 1.0f, sqrt2inv );
    mesh.add_vertex( 0.0f, -1.0f, sqrt2inv );

    // This face-creation order ensure that face 'i' is the only face
    // without vertex 'i'. This is a very convenient way to associate
    // faces and vertices in a simplex.
    mesh.add_face( 1, 3, 2 );
    mesh.add_face( 0, 2, 3 );
    mesh.add_face( 0, 3, 1 );
    mesh.add_face( 0, 1, 2 );
}

frantic::geometry::polymesh3_ptr make_regular_tetrahedron() {
    frantic::geometry::trimesh3 inMesh;
    make_regular_tetrahedron( inMesh );
    return frantic::geometry::trimesh3_to_polymesh3( inMesh );
}

void make_triangle_beveled_vertex_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );          // 0
    mesh.add_vertex( 2.0f, 0.0f, 0.0f );          // 1
    mesh.add_vertex( 1.0f, sqrtf( 3.0f ), 0.0f ); // 2

    frantic::graphics::vector3f centroid( 1.0f, 1.0f / sqrtf( 3.0f ), 0.0f );

    mesh.add_vertex( centroid ); // 3
    mesh.add_vertex( centroid ); // 4
    mesh.add_vertex( centroid ); // 5

    mesh.add_face( 0, 1, 3 );
    mesh.add_face( 1, 4, 3 );
    mesh.add_face( 1, 2, 4 );
    mesh.add_face( 2, 5, 4 );
    mesh.add_face( 0, 3, 5 );
    mesh.add_face( 0, 5, 2 );
    mesh.add_face( 3, 4, 5 );
}

void make_square_beveled_vertex_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0f, 0.0f, 0.0f ); // 0
    mesh.add_vertex( 0.0f, 1.0f, 0.0f ); // 1
    mesh.add_vertex( 1.0f, 1.0f, 0.0f ); // 2
    mesh.add_vertex( 1.0f, 0.0f, 0.0f ); // 3

    mesh.add_vertex( 0.5f, 0.5f, 0.0f ); // 4
    mesh.add_vertex( 0.5f, 0.5f, 0.0f ); // 5
    mesh.add_vertex( 0.5f, 0.5f, 0.0f ); // 6

    mesh.add_face( 0, 1, 4 );
    mesh.add_face( 1, 5, 4 );
    mesh.add_face( 1, 2, 5 );
    mesh.add_face( 2, 6, 5 );
    mesh.add_face( 2, 3, 6 );
    mesh.add_face( 3, 0, 6 );
    mesh.add_face( 0, 4, 6 );
    mesh.add_face( 4, 5, 6 );
}

void make_tetrahedron_beveled_edge_mesh( frantic::geometry::trimesh3& mesh ) {

    float sqrt2inv = 1.0f / std::sqrt( 2.0f );
    mesh.clear();
    mesh.add_vertex( 1.0f, 0.0f, -sqrt2inv );  // 0
    mesh.add_vertex( 1.0f, 0.0f, -sqrt2inv );  // 1
    mesh.add_vertex( -1.0f, 0.0f, -sqrt2inv ); // 2
    mesh.add_vertex( -1.0f, 0.0f, -sqrt2inv ); // 3
    mesh.add_vertex( 0.0f, 1.0f, sqrt2inv );   // 4
    mesh.add_vertex( 0.0f, 1.0f, sqrt2inv );   // 5
    mesh.add_vertex( 0.0f, -1.0f, sqrt2inv );  // 6
    mesh.add_vertex( 0.0f, -1.0f, sqrt2inv );  // 7

    mesh.add_face( 1, 3, 4 );
    mesh.add_face( 0, 3, 1 );
    mesh.add_face( 0, 2, 3 );
    mesh.add_face( 0, 6, 2 );
    mesh.add_face( 2, 6, 3 );
    mesh.add_face( 6, 7, 3 );
    mesh.add_face( 0, 5, 6 );
    mesh.add_face( 6, 5, 7 );
    mesh.add_face( 5, 4, 7 );
    mesh.add_face( 3, 7, 4 );
    mesh.add_face( 0, 1, 4 );
    mesh.add_face( 5, 0, 4 );
}

void make_chained_combined_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );  // 0
    mesh.add_vertex( 1.0f, -1.0f, 0.0f ); // 1
    mesh.add_vertex( 2.0f, 0.0f, 0.0f );  // 2
    mesh.add_vertex( 1.0f, 1.0f, 0.0f );  // 3

    mesh.add_vertex( 0.5f, 0.0f, 0.0f ); // 4
    mesh.add_vertex( 0.5f, 0.0f, 0.0f ); // 5
    mesh.add_vertex( 0.5f, 0.0f, 0.0f ); // 6
    mesh.add_vertex( 0.5f, 0.0f, 0.0f ); // 7

    mesh.add_face( 0, 1, 4 );
    mesh.add_face( 1, 5, 4 );
    mesh.add_face( 1, 6, 5 );
    mesh.add_face( 1, 2, 6 );
    mesh.add_face( 2, 3, 6 );
    mesh.add_face( 3, 7, 6 );
    mesh.add_face( 3, 4, 7 );
    mesh.add_face( 0, 4, 3 );
}

void make_wheel_mesh( frantic::geometry::trimesh3& mesh, size_t wheelCount ) {
    mesh.clear();

    if( wheelCount < 3 ) {
        throw std::runtime_error( "make_wheel_mesh: Error, must have at least 3 exterior vertices." );
    }

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );

    for( int i = 0; i < wheelCount; ++i ) {
        double rotation = ( double( i ) / double( wheelCount ) ) * ( 2.0 * M_PI );
        mesh.add_vertex( float( std::cos( rotation ) ), float( std::sin( rotation ) ), 0.0f );
    }

    for( int i = 0; i < wheelCount; ++i ) {
        mesh.add_face( 0, i + 1, ( ( i + 1 ) % wheelCount ) + 1 );
    }
}

void make_octofan_mesh( frantic::geometry::trimesh3& mesh ) { make_wheel_mesh( mesh, 8 ); }

void make_doubly_covered_triangle( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    mesh.add_vertex( 1.0f, 0.0f, 0.0f );
    mesh.add_vertex( 0.0f, 1.0f, 0.0f );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 2, 1 );
}

void make_saddle_vertex_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( -0.0f, -1.0f, 0.0f );
    mesh.add_vertex( 1.0f, -1.0f, -0.5f );
    mesh.add_vertex( 1.0f, 1.0f, 0.0f );
    mesh.add_vertex( -1.0f, 1.0f, -0.5f );

    mesh.add_vertex( 0.0f, 0.0f, 1.0f );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 2, 3 );

    mesh.add_face( 0, 1, 4 );
    mesh.add_face( 1, 2, 4 );
    mesh.add_face( 2, 3, 4 );
    mesh.add_face( 3, 0, 4 );
}

void make_tetrahedral_cap( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( -2.0f, -1.0f, 0.0f );
    mesh.add_vertex( 2.0f, -1.0f, 0.0f );
    mesh.add_vertex( 0.0f, 2.0f, 0.0f );

    mesh.add_face( 0, 1, 2 );

    add_spikey_face( mesh, 0, 0.5f, 1.0f );
}

void make_sphere_mesh( size_t tesselationFactorIn, frantic::geometry::trimesh3& mesh ) {
    using namespace frantic::geometry;

    mesh.clear();

    const trimesh3::index_t northPoleId = 0;
    const trimesh3::index_t southPoleId = 1;

    mesh.add_vertex( 0.0f, 0.0f, 1.0f );
    mesh.add_vertex( 0.0f, 0.0f, -1.0f );

    const trimesh3::index_t tesselationFactor = trimesh3::index_t( tesselationFactorIn + 3 );

    float phiIncrem = float( M_PI / ( tesselationFactor - 1 ) );
    float thetaIncrem = float( ( 2 * M_PI ) / tesselationFactor );

    for( int phiCounter = 1; phiCounter < ( tesselationFactor - 1 ); ++phiCounter ) {
        for( int thetaCounter = 0; thetaCounter < tesselationFactor; ++thetaCounter ) {
            frantic::graphics::spherical_coords coord( thetaIncrem * thetaCounter, phiIncrem * phiCounter );
            mesh.add_vertex( coord.to_vector3f() );
        }
    }

    trimesh3::index_t firstRoundOffset = 2;
    trimesh3::index_t lastRoundOffset =
        firstRoundOffset + trimesh3::index_t( tesselationFactor * ( tesselationFactor - 3 ) );

    for( int thetaCounter = 0; thetaCounter < tesselationFactor; ++thetaCounter ) {
        mesh.add_face( northPoleId, firstRoundOffset + thetaCounter,
                       firstRoundOffset + ( ( thetaCounter + 1 ) % tesselationFactor ) );
    }

    for( int thetaCounter = 0; thetaCounter < tesselationFactor; ++thetaCounter ) {
        mesh.add_face( lastRoundOffset + thetaCounter, southPoleId,
                       lastRoundOffset + ( ( thetaCounter + 1 ) % tesselationFactor ) );
    }

    trimesh3::index_t currentRoundOffset = firstRoundOffset;
    trimesh3::index_t nextRoundOffset = firstRoundOffset + tesselationFactor;

    for( int phiCounter = 1; phiCounter < ( tesselationFactor - 2 ); ++phiCounter ) {
        for( int thetaCounter = 0; thetaCounter < tesselationFactor; ++thetaCounter ) {
            int thetaCounterNext = ( thetaCounter + 1 ) % tesselationFactor;
            mesh.add_face( currentRoundOffset + thetaCounter, nextRoundOffset + thetaCounter,
                           nextRoundOffset + thetaCounterNext );
            mesh.add_face( currentRoundOffset + thetaCounter, nextRoundOffset + thetaCounterNext,
                           currentRoundOffset + thetaCounterNext );
        }
        currentRoundOffset = nextRoundOffset;
        nextRoundOffset += tesselationFactor;
    }
}

void make_hourglass_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    // apex
    mesh.add_vertex( 0.0f, 0.0f, 0.0f );

    float sqrt3Over2 = std::sqrt( 3.0f ) / 2.0f;

    // 'top' pyramid vertices
    mesh.add_vertex( -1.0f, -sqrt3Over2, 1.0f ); // 1
    mesh.add_vertex( 1.0f, -sqrt3Over2, 1.0f );  // 2
    mesh.add_vertex( 0.0f, sqrt3Over2, 1.0f );   // 3

    // 'bottom' pyramid vertices
    mesh.add_vertex( -1.0f, -sqrt3Over2, -1.0f ); // 4
    mesh.add_vertex( 1.0f, -sqrt3Over2, -1.0f );  // 5
    mesh.add_vertex( 0.0f, sqrt3Over2, -1.0f );   // 6

    // 'top' pyramid
    mesh.add_face( 0, 2, 1 );
    mesh.add_face( 0, 3, 2 );
    mesh.add_face( 0, 3, 2 );
    mesh.add_face( 1, 2, 3 );

    // 'bottom' pyramid
    mesh.add_face( 0, 4, 5 );
    mesh.add_face( 0, 5, 6 );
    mesh.add_face( 0, 6, 4 );
    mesh.add_face( 4, 6, 5 );
}

void make_cube_mesh( frantic::geometry::trimesh3& mesh ) {

    mesh.clear();

    mesh.add_vertex( -1.0f, -1.0f, -1.0f ); // 0
    mesh.add_vertex( 1.0f, -1.0f, -1.0f );  // 1
    mesh.add_vertex( -1.0f, 1.0f, -1.0f );  // 2
    mesh.add_vertex( 1.0f, 1.0f, -1.0f );   // 3
    mesh.add_vertex( -1.0f, -1.0f, 1.0f );  // 4
    mesh.add_vertex( 1.0f, -1.0f, 1.0f );   // 5
    mesh.add_vertex( -1.0f, 1.0f, 1.0f );   // 6
    mesh.add_vertex( 1.0f, 1.0f, 1.0f );    // 7

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 1, 3, 2 );
    mesh.add_face( 1, 5, 3 );
    mesh.add_face( 5, 7, 3 );
    mesh.add_face( 5, 4, 7 );
    mesh.add_face( 4, 6, 7 );
    mesh.add_face( 4, 0, 6 );
    mesh.add_face( 0, 2, 6 );
    mesh.add_face( 2, 3, 6 );
    mesh.add_face( 3, 7, 6 );
    mesh.add_face( 4, 5, 0 );
    mesh.add_face( 5, 1, 0 );
}

void make_sharkfin_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    mesh.add_vertex( 1.0f, 0.0f, 0.0f );
    mesh.add_vertex( 1.0f, 1.0f, 0.0f );
    mesh.add_vertex( 0.0f, 1.0f, 0.0f );
    mesh.add_vertex( 0.5f, 0.5f, 0.5f );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 2, 3 );
    mesh.add_face( 0, 2, 4 );
}

void make_bubble_straw_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0f, 0.0f, 0.0f ); // 0
    mesh.add_vertex( 1.0f, 1.0f, 0.0f ); // 1
    mesh.add_vertex( 1.0f, 1.0f, 0.0f ); // 2
    mesh.add_vertex( 0.0f, 1.0f, 0.0f ); // 3

    mesh.add_vertex( -1.0f, 0.0f, -1.0f ); // 4
    mesh.add_vertex( -1.0f, 0.0f, 1.0f );  // 5

    mesh.add_vertex( 2.0f, 0.0f, -1.0f ); // 6
    mesh.add_vertex( 2.0f, 0.0f, 1.0f );  // 7

    // the 'straw'
    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 2, 3 );
    mesh.add_face( 0, 3, 1 );
    mesh.add_face( 1, 3, 2 );

    // the left 'bubble'
    mesh.add_face( 0, 3, 5 );
    mesh.add_face( 3, 4, 5 );
    mesh.add_face( 0, 5, 4 );
    mesh.add_face( 0, 4, 3 );

    // the right 'bubble'
    mesh.add_face( 1, 7, 2 );
    mesh.add_face( 1, 6, 7 );
    mesh.add_face( 2, 7, 6 );
    mesh.add_face( 1, 2, 6 );
}

void make_genus_mesh( size_t genus, frantic::geometry::trimesh3& mesh ) {
    using namespace frantic::graphics;

    const size_t xSize = genus * 2 + 2;
    const size_t ySize = 4;

    mesh.clear();

    frantic::indexer<size_t, 3> index = frantic::make_indexer( xSize, ySize, size_t( 2 ) );

    mesh.set_vertex_count( index.address_space() );

    for( size_t layer = 0; layer < 2; ++layer ) {
        for( size_t y = 0; y < ySize; ++y ) {
            for( size_t x = 0; x < xSize; ++x ) {
                mesh.get_vertex( index.address( x, y, layer ) ) = vector3f( float( x ), float( y ), float( layer ) );
            }
        }
    }

    // add the grid 'tiles'
    for( size_t x = 0; x < xSize - 1; ++x ) {
        for( size_t y = 0; y < ySize - 1; ++y ) {
            if( x % 2 == 1 && y % 2 == 1 ) {
                add_trimesh_quad( mesh, index.address( x, y, 0 ), index.address( x, y + 1, 0 ),
                                  index.address( x, y + 1, 1 ), index.address( x, y, 1 ) );
                add_trimesh_quad( mesh, index.address( x + 1, y + 1, 0 ), index.address( x + 1, y, 0 ),
                                  index.address( x + 1, y, 1 ), index.address( x + 1, y + 1, 1 ) );
                add_trimesh_quad( mesh, index.address( x + 1, y, 0 ), index.address( x, y, 0 ),
                                  index.address( x, y, 1 ), index.address( x + 1, y, 1 ) );
                add_trimesh_quad( mesh, index.address( x, y + 1, 0 ), index.address( x + 1, y + 1, 0 ),
                                  index.address( x + 1, y + 1, 1 ), index.address( x, y + 1, 1 ) );
            } else {
                add_trimesh_quad( mesh, index.address( x, y, 0 ), index.address( x, y + 1, 0 ),
                                  index.address( x + 1, y + 1, 0 ), index.address( x + 1, y, 0 ) );
                add_trimesh_quad( mesh, index.address( x, y, 1 ), index.address( x + 1, y, 1 ),
                                  index.address( x + 1, y + 1, 1 ), index.address( x, y + 1, 1 ) );
            }
        }
    }

    // add the exteriors
    for( size_t x = 0; x < xSize - 1; ++x ) {
        add_trimesh_quad( mesh, index.address( x, 0, 0 ), index.address( x + 1, 0, 0 ), index.address( x + 1, 0, 1 ),
                          index.address( x, 0, 1 ) );
        add_trimesh_quad( mesh, index.address( x + 1, ySize - 1, 0 ), index.address( x, ySize - 1, 0 ),
                          index.address( x, ySize - 1, 1 ), index.address( x + 1, ySize - 1, 1 ) );
    }

    for( size_t y = 0; y < ySize - 1; ++y ) {
        add_trimesh_quad( mesh, index.address( 0, y + 1, 0 ), index.address( 0, y, 0 ), index.address( 0, y, 1 ),
                          index.address( 0, y + 1, 1 ) );
        add_trimesh_quad( mesh, index.address( xSize - 1, y, 0 ), index.address( xSize - 1, y + 1, 0 ),
                          index.address( xSize - 1, y + 1, 1 ), index.address( xSize - 1, y, 1 ) );
    }
}

void add_spikey_face( frantic::geometry::trimesh3& mesh, size_t faceId, float offsetFactor, float faceHeight ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    vector3 innerFace = add_hole_face( mesh, faceId, offsetFactor );

    const float oneThird = 1.0f / 3.0f;

    vector3f innerFaceVerts[3] = {
        mesh.get_vertex( innerFace[0] ),
        mesh.get_vertex( innerFace[1] ),
        mesh.get_vertex( innerFace[2] ),
    };

    vector3f triangleCenter =
        oneThird * innerFaceVerts[0] + oneThird * innerFaceVerts[1] + oneThird * innerFaceVerts[2];

    plane3f plane = plane3f::from_triangle( innerFaceVerts[0], innerFaceVerts[1], innerFaceVerts[2] );

    vector3::value_type peakVertex = vector3::value_type( mesh.vertex_count() );

    mesh.add_vertex( triangleCenter + plane.normal() * faceHeight );

    for( size_t fVert = 0; fVert < 3; ++fVert ) {
        vector3 peakFace = innerFace;
        peakFace[fVert] = peakVertex;
        mesh.add_face( peakFace );
    }
}

frantic::graphics::vector3 add_hole_face( frantic::geometry::trimesh3& mesh, size_t faceId, float offsetFactor ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    const float oneThird = 1.0f / 3.0f;
    const float baseOffset = frantic::math::clamp( offsetFactor, 0.0f, 1.0f ) * oneThird;
    const float offsetComplement = ( 1.0f - baseOffset ) / 2.0f;

    const vector3 face = mesh.get_face( faceId );

    vector3f faceVerts[3] = {
        mesh.get_vertex( face[0] ),
        mesh.get_vertex( face[1] ),
        mesh.get_vertex( face[2] ),
    };

    vector3 innerFace;

    for( size_t fVert = 0; fVert < 3; ++fVert ) {
        float weights[3] = { offsetComplement, offsetComplement, offsetComplement };
        weights[fVert] = baseOffset;
        vector3f location = weights[0] * faceVerts[0] + weights[1] * faceVerts[1] + weights[2] * faceVerts[2];
        innerFace[fVert] = vector3::value_type( mesh.vertex_count() );
        mesh.add_vertex( location );
    }

    bool first = true;

    for( size_t fVert = 0; fVert < 3; ++fVert ) {
        vector3 edgeFace = face;
        edgeFace[fVert] = innerFace[fVert];

        // replace the existing face first, then add the additional faces
        if( first ) {
            mesh.get_face( faceId ) = edgeFace;
            first = false;
        } else {
            mesh.add_face( edgeFace );
        }

        vector3 pointFace = innerFace;
        pointFace[fVert] = face[fVert];
        std::swap( pointFace[0], pointFace[1] );

        mesh.add_face( pointFace );
    }

    return innerFace;
}

frantic::geometry::trimesh3::index_t add_face_vertex( frantic::geometry::trimesh3& mesh, size_t faceId,
                                                      const frantic::graphics::vector3f& baryCoords ) {
    using namespace frantic::graphics;

    vector3 face = mesh.get_face( faceId );

    vector3f faceLocation = mesh.get_vertex( face[0] ) * baryCoords[0] + mesh.get_vertex( face[1] ) * baryCoords[1] +
                            mesh.get_vertex( face[2] ) * baryCoords[2];

    int newVertexId = static_cast<int>( mesh.vertex_count() );

    mesh.add_vertex( faceLocation );

    mesh.get_face( faceId )[2] = newVertexId;

    mesh.add_face( face[1], face[2], newVertexId );
    mesh.add_face( face[2], face[0], newVertexId );

    return newVertexId;
}

void make_spikey_tetrahedron( frantic::geometry::trimesh3& mesh, float offsetFactor, float faceHeight ) {
    make_regular_tetrahedron( mesh );

    for( size_t i = 0; i < 4; ++i ) {
        add_spikey_face( mesh, i, offsetFactor, faceHeight );
    }
}

void make_holed_tetrahedron( frantic::geometry::trimesh3& mesh, float offsetFactor ) {
    make_regular_tetrahedron( mesh );

    for( size_t i = 0; i < 4; ++i ) {
        add_hole_face( mesh, i, offsetFactor );
    }
}

void make_t_vertex_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.add_vertex( -1.0f, 0.0f, 0.0f );
    mesh.add_vertex( 0.0f, -1.0f, 0.0f );
    mesh.add_vertex( 1.0f, 0.0f, 0.0f );
    mesh.add_vertex( 0.0f, 0.0f, 0.0f );

    mesh.add_face( 0, 1, 3 );
    mesh.add_face( 1, 2, 3 );
    mesh.add_face( 0, 3, 2 );
}

void make_plane_mesh( frantic::geometry::trimesh3& mesh, int xSize, int ySize, float xDelta, float yDelta ) {
    frantic::indexer<int, 2> vertexIndexer( frantic::make_indexer( ySize + 1, xSize + 1 ) );

    for( int y = 0; y < ySize + 1; ++y ) {
        for( int x = 0; x < xSize + 1; ++x ) {
            mesh.add_vertex( xDelta * x, yDelta * y, 0.0f );
        }
    }

    for( int y = 0; y < ySize; ++y ) {
        for( int x = 0; x < xSize; ++x ) {
            mesh.add_face( vertexIndexer( y, x ), vertexIndexer( y + 1, x ), vertexIndexer( y, x + 1 ) );
            mesh.add_face( vertexIndexer( y, x + 1 ), vertexIndexer( y + 1, x ), vertexIndexer( y + 1, x + 1 ) );
        }
    }
}

void make_terrain_mesh( frantic::geometry::trimesh3& mesh, int xSize, int ySize, float xDelta, float yDelta,
                        float zDelta, boost::uint32_t zSeed ) {
    make_plane_mesh( mesh, xSize, ySize, xDelta, yDelta );

    boost::random::mt19937 randomEngine( zSeed );
    boost::random::uniform_real_distribution<float> distro( -zDelta, zDelta );

    for( size_t i = 0; i < mesh.vertex_count(); ++i ) {
        mesh.get_vertex( i ).z = distro( randomEngine );
    }
}

void make_elongated_tetrahedron_mesh( frantic::geometry::trimesh3& mesh ) {
    mesh.clear();

    mesh.add_vertex( 0.0692248f, -0.136369f, -0.154853f );
    mesh.add_vertex( -0.0825064f, 0.107311f, 0.17459f );
    mesh.add_vertex( 0.00620087f, -0.0687743f, -0.0374749f );
    mesh.add_vertex( 0.0291289f, -0.0527215f, -0.0177045f );

    mesh.add_face( 1, 3, 0 );
    mesh.add_face( 3, 2, 0 );
    mesh.add_face( 1, 0, 2 );
    mesh.add_face( 3, 1, 2 );
}

frantic::geometry::polymesh3_ptr make_sliver_polymesh() {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder polymeshBuilder;

    polymeshBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

    int polygon[] = { 0, 1, 1 };

    polymeshBuilder.add_polygon( polygon, 3 );

    return polymeshBuilder.finalize();
}

frantic::geometry::polymesh3_ptr make_triangle_polymesh() {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder polymeshBuilder;

    polymeshBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );

    int polygon[] = { 0, 1, 2 };

    polymeshBuilder.add_polygon( polygon, 3 );

    return polymeshBuilder.finalize();
}

frantic::geometry::polymesh3_ptr make_quad_polymesh() {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder polymeshBuilder;

    polymeshBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

    int polygon[] = { 0, 3, 2, 1 };

    polymeshBuilder.add_polygon( polygon, 4 );

    return polymeshBuilder.finalize();
}

frantic::geometry::polymesh3_ptr make_split_quad_polymesh() {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder polymeshBuilder;

    polymeshBuilder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
    polymeshBuilder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

    int polygon[] = { 0, 1, 2 };
    polymeshBuilder.add_polygon( polygon, 3 );

    polygon[1] = 2;
    polygon[2] = 3;
    polymeshBuilder.add_polygon( polygon, 3 );

    return polymeshBuilder.finalize();
}

frantic::geometry::polymesh3_ptr make_cube_polymesh( const std::set<size_t>& removedFaces ) {
    using namespace frantic::graphics;
    using namespace frantic::geometry;

    polymesh3_builder builder;

    builder.add_vertex( -1.0f, -1.0f, -1.0f );
    builder.add_vertex( 1.0f, -1.0f, -1.0f );
    builder.add_vertex( -1.0f, 1.0f, -1.0f );
    builder.add_vertex( 1.0f, 1.0f, -1.0f );
    builder.add_vertex( -1.0f, -1.0f, 1.0f );
    builder.add_vertex( 1.0f, -1.0f, 1.0f );
    builder.add_vertex( -1.0f, 1.0f, 1.0f );
    builder.add_vertex( 1.0f, 1.0f, 1.0f );

    int ids[4];

#define __ADD_TEST_FACE__( f, i0, i1, i2, i3 )                                                                         \
    do {                                                                                                               \
        if( removedFaces.count( f ) == 0 ) {                                                                           \
            ids[0] = i0;                                                                                               \
            ids[1] = i1;                                                                                               \
            ids[2] = i2;                                                                                               \
            ids[3] = i3;                                                                                               \
            builder.add_polygon( ids, 4 );                                                                             \
        }                                                                                                              \
    } while( 0 )
    __ADD_TEST_FACE__( 0, 1, 5, 7, 3 ); // X+ face
    __ADD_TEST_FACE__( 1, 2, 6, 4, 0 ); // X- face
    __ADD_TEST_FACE__( 2, 3, 7, 6, 2 ); // Y+ face
    __ADD_TEST_FACE__( 3, 0, 4, 5, 1 ); // Y-face
    __ADD_TEST_FACE__( 4, 4, 6, 7, 5 ); // Z+ face
    __ADD_TEST_FACE__( 5, 1, 3, 2, 0 ); // Z- face
#undef __ADD_TEST_FACE__

    return builder.finalize();
}

frantic::geometry::polymesh3_ptr make_sphere_polymesh( size_t tesselationFactor ) {
    using namespace frantic::geometry;
    trimesh3 mesh;
    make_sphere_mesh( tesselationFactor, mesh );
    return frantic::geometry::trimesh3_to_polymesh3( mesh );
}

std::vector<frantic::geometry::polymesh3_ptr> make_separated_quad( float separation ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<frantic::geometry::polymesh3_ptr> meshes( 2 );

    polymesh3_builder firstBuilder;

    firstBuilder.add_vertex( vector3f( -1.0f, 0.0f, 0.0f ) );
    firstBuilder.add_vertex( vector3f( -1.0f, 1.0f, 0.0f ) );
    firstBuilder.add_vertex( vector3f( 0.0f, 1.0f, 0.0f ) );
    firstBuilder.add_vertex( vector3f( 0.0f, 0.0f, 0.0f ) );

    int face[] = { 0, 1, 2, 3 };
    firstBuilder.add_polygon( face, 4 );
    meshes[0] = firstBuilder.finalize();

    polymesh3_builder secondBuilder;

    secondBuilder.add_vertex( vector3f( 0.0f + separation, 0.0f, 0.0f ) );
    secondBuilder.add_vertex( vector3f( 0.0f + separation, 1.0f, 0.0f ) );
    secondBuilder.add_vertex( vector3f( 1.0f + separation, 1.0f, 0.0f ) );
    secondBuilder.add_vertex( vector3f( 1.0f + separation, 0.0f, 0.0f ) );

    secondBuilder.add_polygon( face, 4 );
    meshes[1] = secondBuilder.finalize();

    return meshes;
}

frantic::geometry::polymesh3_ptr make_welded_quad() {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder builder;

    builder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
    builder.add_vertex( vector3f( -1.0, 1.0, 0.0 ) );
    builder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
    builder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    builder.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
    builder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

    int firstFace[] = { 0, 1, 2, 3 };
    builder.add_polygon( firstFace, 4 );

    int secondFace[] = { 3, 2, 4, 5 };
    builder.add_polygon( secondFace, 4 );

    return builder.finalize();
}

frantic::geometry::polymesh3_ptr make_combined_quad() {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder builder;

    builder.add_vertex( vector3f( -1.0, 0.0, 0.0 ) );
    builder.add_vertex( vector3f( -1.0, 1.0, 0.0 ) );
    builder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
    builder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    builder.add_vertex( vector3f( 0.0, 0.0, 0.0 ) );
    builder.add_vertex( vector3f( 0.0, 1.0, 0.0 ) );
    builder.add_vertex( vector3f( 1.0, 1.0, 0.0 ) );
    builder.add_vertex( vector3f( 1.0, 0.0, 0.0 ) );

    int firstFace[] = { 0, 1, 2, 3 };
    builder.add_polygon( firstFace, 4 );
    int secondFace[] = { 4, 5, 6, 7 };
    builder.add_polygon( secondFace, 4 );

    return builder.finalize();
}