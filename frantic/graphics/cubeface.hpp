// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/algorithm/clamp.hpp>

#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/size2.hpp>
#include <frantic/graphics2d/size2f.hpp>
#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace frantic {
namespace graphics {

// TODO: mw - This hemisphere face thing seems like a really bad idea to me, should use the cube face stuff instead to
//            implement it
namespace hemisphere_face {
enum default_hemisphere_face { HF_X_POS = 0, HF_X_NEG = 1, HF_Y_POS = 2, HF_Y_NEG = 3, HF_Z_NEG = 4, HF_NO_FACE = 5 };
}

namespace cube_face {
enum default_cube_face { CF_X_POS = 0, CF_X_NEG = 1, CF_Y_POS = 2, CF_Y_NEG = 3, CF_Z_POS = 4, CF_Z_NEG = 5 };
struct cube_face_mapping {
    int xTileCount;       // Number of tiles on x axis
    int yTileCount;       // Number of tiles on y axis
    bool rotatePositiveZ; // Wether or not to rotate positive z tile (needed for vertical cross)
    int xTile[6];         // X tile number for each cube face (added to optimize from_cube_face_map)
    int yTile[6];         // Y tile number for each cube face (added to optimize from_cube_face_map)
    float xOffset[6];     // Normalized [0,1] x offset for each cube face (see default_cube_face)
    float yOffset[6];     // Normalized [0,1] y offset for each cube face (see default_cube_face)
};

const cube_face_mapping CROSS_VERTICAL = {

    // MAPPING
    //
    //         ______
    //        |      |
    //        |  Y+  |
    //  ______|______|______
    // |      |      |      |
    // |  X-  |  Z-  |  X+  |
    // |______|______|______|
    //        |      |
    //        |  Y-  |
    //        |______|
    //        |      |
    //        |  +Z  | <---- Rotated 180 degrees
    //        |______|

    3,    // Number of x-tiles
    4,    // Number of y-tiles
    true, // Rotate Positive Z tile
          // ---------- Tile numbers -----------
          // Face					X+	X-	Y+	Y-	Z+	Z-
    /* X tile number*/ { 2, 0, 1, 1, 1, 1 },
    /* Y tile number*/ { 2, 2, 3, 1, 0, 2 },
    // --------- X tile offsets ----------
    {
        2.f / 3.f, // X-Positive
        0.f / 3.f, // X-Negative
        1.f / 3.f, // Y-Positive
        1.f / 3.f, // Y-Negative
        1.f / 3.f, // Z-Positive
        1.f / 3.f  // Z-Negative
    },
    // --------- Y tile offsets ----------
    {
        2.f / 4.f, // X-Positive
        2.f / 4.f, // X-Negative
        3.f / 4.f, // Y-Positive
        1.f / 4.f, // Y-Negative
        0.f / 4.f, // Z-Positive
        2.f / 4.f  // Z-Negative
    } };

const cube_face_mapping CROSS_HORIZONTAL = {

    // MAPPING
    //
    //         ______
    //        |      |
    //        |  Y+  |
    //  ______|______|_____________
    // |      |      |      |      |
    // |  X-  |  Z-  |  X+  |  Z+  |
    // |______|______|______|______|
    //        |      |
    //        |  Y-  |
    //        |______|

    4,     // Number of x-tiles
    3,     // Number of y-tiles
    false, // Do not rotate Positive Z tile
           // ---------- Tile numbers -----------
           // Face					X+	X-	Y+	Y-	Z+	Z-
    /* X tile number*/ { 2, 0, 1, 1, 3, 1 },
    /* Y tile number*/ { 1, 1, 2, 0, 1, 1 },
    // --------- X tile offsets ----------
    {
        2.f / 4.f, // X-Positive
        0.f / 4.f, // X-Negative
        1.f / 4.f, // Y-Positive
        1.f / 4.f, // Y-Negative
        3.f / 4.f, // Z-Positive
        1.f / 4.f  // Z-Negative
    },
    // --------- Y tile offsets ----------
    {
        1.f / 3.f, // X-Positive
        1.f / 3.f, // X-Negative
        2.f / 3.f, // Y-Positive
        0.f / 3.f, // Y-Negative
        1.f / 3.f, // Z-Positive
        1.f / 3.f  // Z-Negative
    } };

const cube_face_mapping RECTANGLE_HORIZONTAL = {

    // MAPPING
    //  ____________________
    // |      |      |      |
    // |  X+  |  Y+  |  Z-  |
    // |______|______|______|
    // |      |      |      |
    // |  X-  |  Y-  |  Z+  |  <--- Z+ is the front in Gelato
    // |______|______|______|

    3,     // Number of x-tiles
    2,     // Number of y-tiles
    false, // Do not rotate Positive Z tile

    // The Z-Positive and Z-Negative tiles are
    // flipped because Gelato uses the left handed
    // coordinate system

    // ---------- Tile numbers -----------
    // Face					X+	X-	Y+	Y-	Z+	Z-
    /* X tile number*/ { 0, 0, 1, 1, 2, 2 },
    /* Y tile number*/ { 1, 0, 1, 0, 0, 1 },
    // --------- X tile offsets ----------
    {
        0.f / 3.f, // X-Positive
        0.f / 3.f, // X-Negative
        1.f / 3.f, // Y-Positive
        1.f / 3.f, // Y-Negative
        2.f / 3.f, // Z-Positive
        2.f / 3.f  // Z-Negative
    },
    // --------- Y tile offsets ----------
    {
        1.f / 2.f, // X-Positive
        0.f / 2.f, // X-Negative
        1.f / 2.f, // Y-Positive
        0.f / 2.f, // Y-Negative
        0.f / 2.f, // Z-Positive
        1.f / 2.f  // Z-Negative
    } };

const cube_face_mapping STRIP_VERTICAL = {

    // MAPPING
    //  ______
    // |      |
    // |  X+  |
    // |______|
    // |      |
    // |  X-  |
    // |______|
    // |      |
    // |  Y+  |
    // |______|
    // |      |
    // |  Y-  |
    // |______|
    // |      |
    // |  Z-  | <----
    // |______|      \______ Swapped because Gelato uses
    // |      |      /       left hand coordinate system
    // |  Z+  | <----
    // |______|

    1,     // Number of x-tiles
    6,     // Number of y-tiles
    false, // Do not rotate Positive Z tile
           // ---------- Tile numbers -----------
           // Face					X+	X-	Y+	Y-	Z+	Z-
    /* X tile number*/ { 0, 0, 0, 0, 0, 0 },
    /* Y tile number*/ { 5, 4, 3, 2, 0, 1 },
    // --------- X tile offsets ----------
    {
        0.f / 1.f, // X-Positive
        0.f / 1.f, // X-Negative
        0.f / 1.f, // Y-Positive
        0.f / 1.f, // Y-Negative
        0.f / 1.f, // Z-Positive
        0.f / 1.f  // Z-Negative
    },
    // --------- Y tile offsets ----------
    {
        5.f / 6.f, // X-Positive
        4.f / 6.f, // X-Negative
        3.f / 6.f, // Y-Positive
        2.f / 6.f, // Y-Negative
        0.f / 6.f, // Z-Positive
        1.f / 6.f  // Z-Negative
    } };

} // namespace cube_face

// This returns which cube face the vector is pointing at
inline cube_face::default_cube_face get_cube_face( const frantic::graphics::vector3f& v ) {
    float ax = fabsf( v.x ), ay = fabsf( v.y ), az = fabsf( v.z );
    if( ay > ax ) {     // faces 2, 3, 4, 5
        if( ay > az ) { // faces 2, 3
            if( v.y > 0 )
                return cube_face::CF_Y_POS;
            else
                return cube_face::CF_Y_NEG;
        } else { // faces 4, 5
            if( v.z > 0 )
                return cube_face::CF_Z_POS;
            else
                return cube_face::CF_Z_NEG;
        }
    } else {            // faces 0, 1, 4, 5
        if( ax > az ) { // faces 0, 1
            if( v.x > 0 )
                return cube_face::CF_X_POS;
            else
                return cube_face::CF_X_NEG;
        } else { // faces 4, 5
            if( v.z > 0 )
                return cube_face::CF_Z_POS;
            else
                return cube_face::CF_Z_NEG;
        }
    }
}

// This returns the 2D coordinate that the vector is pointing at withing the given cube face
// The returned (x,y) coordinates are in the [-1,1] interval for each coordinate.
// The top left is at (-1, -1)
inline frantic::graphics2d::vector2f get_cube_face_coordinate( const vector3f& v,
                                                               cube_face::default_cube_face cubeFace ) {
    switch( cubeFace ) {
    case cube_face::CF_X_POS: // (1, 0, 0)
        return frantic::graphics2d::vector2f( v.z / v.x, v.y / v.x );
    case cube_face::CF_X_NEG: // (-1, 0, 0)
        return frantic::graphics2d::vector2f( v.z / v.x, -v.y / v.x );
    case cube_face::CF_Y_POS: // (0, 1, 0)
        return frantic::graphics2d::vector2f( v.x / v.y, v.z / v.y );
    case cube_face::CF_Y_NEG: // (0, -1, 0)
        return frantic::graphics2d::vector2f( -v.x / v.y, v.z / v.y );
    case cube_face::CF_Z_POS: // (0, 0, 1)
        return frantic::graphics2d::vector2f( -v.x / v.z, v.y / v.z );
    case cube_face::CF_Z_NEG: // (0, 0, -1)
        return frantic::graphics2d::vector2f( -v.x / v.z, -v.y / v.z );
    default:
        throw std::runtime_error( "vector3f::get_cube_face_coordinate: Provided invalid cube face value " +
                                  boost::lexical_cast<std::string>( cubeFace ) );
    }
}

// This returns the 2D coordinate that the vector is pointing at withing the given cube face, and the z-depth in that
// cubeface direction.
// The returned (x,y) coordinates are in the [-1,1] interval for each coordinate. The returned z-depth is in
// cameraspace.
// The top left is at (-1, -1)
inline frantic::graphics2d::vector2f
get_cube_face_coordinate_and_zdepth( const vector3f& v, cube_face::default_cube_face cubeFace, float& outZDepth ) {
    switch( cubeFace ) {
    case cube_face::CF_X_POS: // (1, 0, 0)
        outZDepth = v.x;
        return frantic::graphics2d::vector2f( v.z / v.x, v.y / v.x );
    case cube_face::CF_X_NEG: // (-1, 0, 0)
        outZDepth = -v.x;
        return frantic::graphics2d::vector2f( v.z / v.x, -v.y / v.x );
    case cube_face::CF_Y_POS: // (0, 1, 0)
        outZDepth = v.y;
        return frantic::graphics2d::vector2f( v.x / v.y, v.z / v.y );
    case cube_face::CF_Y_NEG: // (0, -1, 0)
        outZDepth = -v.y;
        return frantic::graphics2d::vector2f( -v.x / v.y, v.z / v.y );
    case cube_face::CF_Z_POS: // (0, 0, 1)
        outZDepth = v.z;
        return frantic::graphics2d::vector2f( -v.x / v.z, v.y / v.z );
    case cube_face::CF_Z_NEG: // (0, 0, -1)
        outZDepth = -v.z;
        return frantic::graphics2d::vector2f( -v.x / v.z, -v.y / v.z );
    default:
        throw std::runtime_error( "vector3f::get_cube_face_coordinate: Provided invalid cube face value " +
                                  boost::lexical_cast<std::string>( cubeFace ) );
    }
}

//							 Y		  *
//  _________                |   /
// |  *      |      \        | /
// |    X+   |  ----->  -----+----- X
// |_________|      /       /|
//                        /  |
//						Z
//
// TODO: Should take in a size and absolute coordinate instead of normalized
// v is a normalized coordinate between 0 and 1 where (0,0) is bottom left
// cubeFace is the cube face that the point is on
// Returns the 3d view vector that the input corresponds to
inline frantic::graphics::vector3f from_cube_face_coordinate( const frantic::graphics2d::vector2f& v,
                                                              const cube_face::default_cube_face cubeFace ) {
    float y = v.y * 2 - 1;
    float x = v.x * 2 - 1;

    switch( cubeFace ) {
    case cube_face::CF_X_POS:
        return vector3f( 1, y, x );
    case cube_face::CF_X_NEG:
        return vector3f( -1, y, -x );
    case cube_face::CF_Y_POS:
        return vector3f( x, 1, y );
    case cube_face::CF_Y_NEG:
        return vector3f( x, -1, -y );
    case cube_face::CF_Z_POS:
        return vector3f( -x, y, 1 );
    case cube_face::CF_Z_NEG:
        return vector3f( x, y, -1 );
    default:
        throw std::runtime_error( "vector3f::from_cube_face_coordinate: Invalid cubeface" );
    }
}

inline frantic::graphics::vector3f from_cube_face_map( const frantic::graphics2d::vector2f& v,
                                                       const frantic::graphics2d::size2f& size,
                                                       const cube_face::cube_face_mapping& cfMap, bool& withinMap ) {
    float xPixelsPerTile = size.xsize / cfMap.xTileCount;
    float yPixelsPerTile = size.ysize / cfMap.yTileCount;
    // Get tile of v
    int tileX = (int)( v.x / xPixelsPerTile );
    int tileY = (int)( v.y / yPixelsPerTile );

    // Calculate the normalized [0,1] position in the tile
    frantic::graphics2d::vector2f posInTile( ( v.x - tileX * xPixelsPerTile ) / xPixelsPerTile,
                                             ( v.y - tileY * yPixelsPerTile ) / yPixelsPerTile );

    withinMap = true;

    for( int i = 0; i < 6; i++ ) {
        cube_face::default_cube_face cubeFace = (cube_face::default_cube_face)i;
        if( tileX == cfMap.xTile[cubeFace] && tileY == cfMap.yTile[cubeFace] ) {
            // Rotate the Z tile if specified in the cube face map
            if( cfMap.rotatePositiveZ && cubeFace == cube_face::CF_Z_POS ) {
                posInTile.x = 1 - posInTile.x;
                posInTile.y = 1 - posInTile.y;
            }
            return from_cube_face_coordinate( posInTile, cubeFace );
        }
    }

    withinMap = false;
    return vector3f( 0 );
}

inline frantic::graphics::vector3f from_cube_face_map( const frantic::graphics2d::vector2& v,
                                                       const frantic::graphics2d::size2& size,
                                                       const cube_face::cube_face_mapping& cfMap, bool& withinMap ) {
    return from_cube_face_map( frantic::graphics2d::vector2f( v.x + 0.5f, v.y + 0.5f ),
                               frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ), cfMap,
                               withinMap );
}

inline frantic::graphics::vector3f from_axis( int axis ) {
    switch( axis ) {
    case 0:
        return vector3f( 1.f, 0.f, 0.f );
    case 1:
        return vector3f( 0.f, 1.f, 0.f );
    case 2:
        return vector3f( 0.f, 0.f, 1.f );
    default:
        return vector3f();
    }
}

inline frantic::graphics2d::vector2f to_cube_face_map( const frantic::graphics::vector3f& v,
                                                       const frantic::graphics2d::size2f& size,
                                                       const cube_face::cube_face_mapping& cfMap, float edgeCut ) {
    // Get the cube face and coordinate
    cube_face::default_cube_face cubeFace = get_cube_face( v );
    frantic::graphics2d::vector2f result = get_cube_face_coordinate( v, cubeFace );

    // Rotate the Z-Positive tile if necessary
    if( cfMap.rotatePositiveZ && cubeFace == cube_face::CF_Z_POS ) {
        result.x = -result.x;
        result.y = -result.y;
    }

    // Calculates the scaling factor for the edge cut
    if( edgeCut > 0 ) {
        float xPixels = size.xsize / cfMap.xTileCount;
        float yPixels = size.ysize / cfMap.yTileCount;
        float xRatio = ( xPixels - 2.f * edgeCut ) / xPixels;
        float yRatio = ( yPixels - 2.f * edgeCut ) / yPixels;

        result.x *= xRatio;
        result.y *= yRatio;
    }

    // Normalize coordinates to [0,1]
    result.x = ( result.x + 1.f ) / 2.f;
    result.y = ( result.y + 1.f ) / 2.f;

    // Scale the result to the cube face map
    result.x /= cfMap.xTileCount;
    result.y /= cfMap.yTileCount;

    // Translate the coordinate to its position in the cube face map
    result.x += cfMap.xOffset[cubeFace];
    result.y += cfMap.yOffset[cubeFace];

    // Scale result to image size
    result.x *= size.xsize - 1.f;
    result.y *= size.ysize;
    result.y = boost::algorithm::clamp( result.y, 0.f, size.ysize - 1.f );

    return result;
}

inline frantic::graphics2d::vector2 to_cube_face_map( const vector3f& vec, const frantic::graphics2d::size2& size,
                                                      const cube_face::cube_face_mapping& cfMap, float edgeCut ) {
    frantic::graphics2d::vector2f v = to_cube_face_map(
        vec, frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ), cfMap, edgeCut );
    return frantic::graphics2d::vector2( (int)floorf( v.x ), (int)floorf( v.y ) );
}

// This returns which hemisphere face the vector is pointing at
inline hemisphere_face::default_hemisphere_face get_hemisphere_face( const vector3f& v ) {
    float ax = fabsf( v.x ), ay = fabsf( v.y ), az = fabsf( v.z );
    if( ay > ax ) {       // faces HF_Y_POS, HF_Y_NEG, HF_Z_NEG, HF_NO_FACE
        if( ay > az ) {   // faces HF_Y_POS, HF_Y_NEG, HF_NO_FACE
            if( v.z > 0 ) // if pointed more than half way up the side
                return hemisphere_face::HF_NO_FACE;
            else if( v.y > 0 )
                return hemisphere_face::HF_Y_POS;
            else
                return hemisphere_face::HF_Y_NEG;
        } else { // faces HF_Z_NEG, HF_NO_FACE
            if( v.z > 0 )
                return hemisphere_face::HF_NO_FACE;
            else
                return hemisphere_face::HF_Z_NEG;
        }
    } else {              // faces HF_X_POS, HF_X_NEG, HF_Z_NEG, HF_NO_FACE
        if( ax > az ) {   // faces HF_X_POS, HF_X_NEG, HF_NO_FACE
            if( v.z > 0 ) // if pointed more than half way up the side
                return hemisphere_face::HF_NO_FACE;
            else if( v.x > 0 )
                return hemisphere_face::HF_X_POS;
            else
                return hemisphere_face::HF_X_NEG;
        } else { // faces HF_Z_NEG, HF_NO_FACE
            if( v.z > 0 )
                return hemisphere_face::HF_NO_FACE;
            else
                return hemisphere_face::HF_Z_NEG;
        }
    }
}

// This gets all the hemisphere faces you need to touch if you care about values
// within the given ratioTolerance. If your draw_point affects a radius 4 region,
// and the hemisphere full face width is 512. This means that your ratioTolerance
// should be (512-2*4)/512. outHemisphereFaces should be an array of at least 5 ints.
// The return value is the number of hemisphere faces set in the output array.
inline int get_hemisphere_faces( const vector3f& v, float ratioTolerance,
                                 hemisphere_face::default_hemisphere_face* outHemisphereFaces ) {
    int hemisphereFaceCount = 0;
    float ax = fabsf( v.x ), ay = fabsf( v.y ), az = fabsf( v.z );

    if( ay >= ratioTolerance * ax ) {     // faces HF_Y_POS, HF_Y_NEG, HF_Z_NEG, HF_NO_FACE
        if( ay >= ratioTolerance * az ) { // faces HF_Y_POS, HF_Y_NEG, HF_NO_FACE
            if( v.z > 0 )                 // if pointed more than half way up the side
            {                             // no face to add
            } else if( v.y > 0 )
                outHemisphereFaces[hemisphereFaceCount++] = hemisphere_face::HF_Y_POS;
            else
                outHemisphereFaces[hemisphereFaceCount++] = hemisphere_face::HF_Y_NEG;
        }
        if( az >= ratioTolerance * ay ) { // faces HF_Z_NEG, HF_NO_FACE
            if( v.z > 0 ) {               // no face to add
            } else
                outHemisphereFaces[hemisphereFaceCount++] = hemisphere_face::HF_Z_NEG;
        }
    }

    if( ax >= ratioTolerance * ay ) {     // faces HF_X_POS, HF_X_NEG, HF_Z_NEG, HF_NO_FACE
        if( ax >= ratioTolerance * az ) { // faces HF_X_POS, HF_X_NEG, HF_NO_FACE
            if( v.z > 0 )                 // if pointed more than half way up the side
            {                             // no face to add
            } else if( v.x > 0 ) {
                outHemisphereFaces[hemisphereFaceCount++] = hemisphere_face::HF_X_POS;
            } else
                outHemisphereFaces[hemisphereFaceCount++] = hemisphere_face::HF_X_NEG;
        }
        if( az >= ratioTolerance * ax ) { // faces HF_Z_NEG, HF_NO_FACE
            if( v.z > 0 ) {               // no face to add
            } else if( std::find( outHemisphereFaces, outHemisphereFaces + hemisphereFaceCount,
                                  hemisphere_face::HF_Z_NEG ) == outHemisphereFaces + hemisphereFaceCount )
                outHemisphereFaces[hemisphereFaceCount++] = hemisphere_face::HF_Z_NEG;
        }
    }

    return hemisphereFaceCount;
}

// This gets all the cube faces you need to touch if you care about values
// within the given ratioTolerance.  Let's say your draw_point affects a radius 4 region,
// and the cube face width is 512.  This means that ratioTolerance should be (512-2*4)/512.
// outCubeFaces should be an array of at least 6 ints.
// The return value is the number of cube facess set in the output array
inline int get_cube_faces( const vector3f& v, float ratioTolerance, cube_face::default_cube_face* outCubeFaces ) {
    int cubeFaceCount = 0;
    float ax = fabsf( v.x ), ay = fabsf( v.y ), az = fabsf( v.z );
    if( ay >= ratioTolerance * ax ) {     // faces 2, 3, 4, 5
        if( ay >= ratioTolerance * az ) { // faces 2, 3
            if( v.y > 0 )
                outCubeFaces[cubeFaceCount++] = cube_face::CF_Y_POS;
            else
                outCubeFaces[cubeFaceCount++] = cube_face::CF_Y_NEG;
        }
        if( az >= ratioTolerance * ay ) { // faces 4, 5
            if( v.z > 0 )
                outCubeFaces[cubeFaceCount++] = cube_face::CF_Z_POS;
            else
                outCubeFaces[cubeFaceCount++] = cube_face::CF_Z_NEG;
        }
    }
    if( ax >= ratioTolerance * ay ) {     // faces 0, 1, 4, 5
        if( ax >= ratioTolerance * az ) { // faces 0, 1
            if( v.x > 0 )
                outCubeFaces[cubeFaceCount++] = cube_face::CF_X_POS;
            else
                outCubeFaces[cubeFaceCount++] = cube_face::CF_X_NEG;
        }
        if( az >= ratioTolerance * ax ) { // faces 4, 5
            // In this case, the face may already have been found by the previous test
            if( v.z > 0 ) {
                if( std::find( outCubeFaces, outCubeFaces + cubeFaceCount, cube_face::CF_Z_POS ) ==
                    outCubeFaces + cubeFaceCount )
                    outCubeFaces[cubeFaceCount++] = cube_face::CF_Z_POS;
            } else {
                if( std::find( outCubeFaces, outCubeFaces + cubeFaceCount, cube_face::CF_Z_NEG ) ==
                    outCubeFaces + cubeFaceCount )
                    outCubeFaces[cubeFaceCount++] = cube_face::CF_Z_NEG;
            }
        }
    }
    return cubeFaceCount;
}

// This returns the 2D coordinate that the vector is pointing at within the given cube face.
// The returned (x,y) coordinates are in [-1,1].
// The top left is at coord=(-1,-1).
// If the vector does not actually point at the given hemisphere face, the output is not usable.
// Ensure that you get the correct hemisphereFace by calling get_hemisphere_face().
// TODO: The required call to get_hemisphere_face could be encapsulated in this call?
inline frantic::graphics2d::vector2f
get_hemisphere_face_coordinate( const vector3f& v, hemisphere_face::default_hemisphere_face hemisphereFace ) {
    switch( hemisphereFace ) {
    case hemisphere_face::HF_X_POS:
        return frantic::graphics2d::vector2f( v.z / v.x, v.y / v.x );
    case hemisphere_face::HF_X_NEG:
        return frantic::graphics2d::vector2f( v.z / v.x, -v.y / v.x );
    case hemisphere_face::HF_Y_POS:
        return frantic::graphics2d::vector2f( v.x / v.y, v.z / v.y );
    case hemisphere_face::HF_Y_NEG:
        return frantic::graphics2d::vector2f( -v.x / v.y, v.z / v.y );
    case hemisphere_face::HF_Z_NEG:
        return frantic::graphics2d::vector2f( -v.x / v.z, -v.y / v.z );
    default:
        throw std::runtime_error( "vector2f::get_hemisphere_face_coordinate: Provided invalid hemisphere face value " +
                                  boost::lexical_cast<std::string>( hemisphereFace ) );
    }
}
} // namespace graphics
} // namespace frantic
