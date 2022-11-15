// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frantic {
namespace graphics2d {
namespace image {

// -------------------------------------------------------------------
//
// get_pixel functions operate in two formats:
//
//		integer: relating to specific x, y coordinates.
//		(0,height-1)
//		  ^
//		  |
//		  |
//		  |
//		  |
//		  |
//		  |
//		  |
//		  |
//		  |
//		(0,0)---------------------->(width-1,0)
//
//
//		normalized (float): relating to sub-pixel coordinates
//		              (0,1)
//						^
//						|
//						|
//						|
//						|
//		(-1,0)<-------(0,0)------->(1,0)
//						|
//						|
//						|
//						|
//						v
//					 (-1,0)
//
//
// Take care to know which one of these formats you are working in
//
// -------------------------------------------------------------------

// An optimized version of the special case b-spline cubic
struct BSplineCubicFilter {
    float radius() { return 2.0; }
    float operator()( float x ) {
        float xx = x * x;
        if( x < 0 )
            x = -x;
        if( x < 1 ) {
            return ( 1.f / 6.f ) * ( 3 * x * xx - 6 * xx + 4 );
        } else if( x < 2.0 ) {
            return ( 1.f / 6.f ) * ( -x * xx + 6 * xx - 12 * x + 8 );
        }
        return 0;
    }
};

/**
 * This is an SSE version of the BSplineCubicFilter. It is slightly faster,
 * but probably not enough to warrant its use.
 */
#ifdef _FVEC_H_INCLUDED

#pragma warning( push )
#pragma warning( disable : 4239 )
struct BSplineCubicFilter4 {
    float radius() { return 2.0; }
    void operator()( F32vec4& x ) {
        x = simd_max( x, x * F32vec4( -1.0f ) );

        F32vec4 xx = x * x;
        F32vec4 xxx = x * xx;

        F32vec4 lt1 = cmplt( x, F32vec4( 1.0f ) );
        F32vec4 lt2 = cmplt( x, F32vec4( 2.0f ) );

        lt2 ^= lt1;

        F32vec4 r1( 4.f );
        r1 += F32vec4( 3.0f ) * xxx;
        r1 -= F32vec4( 6.f ) * xx;

        F32vec4 r2( 8.f );
        r2 -= xxx;
        r2 += F32vec4( 6.f ) * xx;
        r2 -= F32vec4( 12.f ) * x;

        x = F32vec4( 1 / 6.f ) * ( ( r1 & lt1 ) | ( r2 & lt2 ) );
    }
};
#pragma warning( pop )

#endif

template <typename ChannelDataType>
void set_pixel( vector2 pixel, const ChannelDataType& value, std::vector<ChannelDataType>& channel,
                const size2& size ) {
    if( size.is_coord_valid( pixel ) )
        channel[size.get_index( pixel )] = value;
}

template <typename ChannelDataType>
void set_pixel( int x, int y, const ChannelDataType& value, std::vector<ChannelDataType>& channel, const size2& size ) {
    set_pixel( vector2( x, y ), value, channel, size );
}

template <typename ChannelDataType>
const ChannelDataType& get_pixel( const vector2& pixel, const std::vector<ChannelDataType>& channel,
                                  const size2& size ) {
    if( size.is_coord_valid( pixel ) ) {
        return channel[size.get_index( pixel )];
    } else
        throw std::runtime_error( "image.get_pixel: Requested pixel index, " + pixel.str() +
                                  ", out of bounds of image size, " + size.str() );
}

template <typename ChannelDataType>
const ChannelDataType& get_pixel( int x, int y, const std::vector<ChannelDataType>& channel, const size2& size ) {
    return get_pixel( vector2( x, y ), channel, size );
}

template <typename ChannelDataType>
const ChannelDataType& get_pixel_fast( const vector2& pixel, const std::vector<ChannelDataType>& channel,
                                       const size2& size ) {
    return channel[size.get_index( pixel )];
}

template <typename ChannelDataType>
const ChannelDataType& get_pixel_fast( int x, int y, const std::vector<ChannelDataType>& channel, const size2& size ) {
    return get_pixel_fast( vector2( x, y ), channel, size );
}

// Gets a pixel value from normalized coordinates which are in [-1, 1]
// This wraps along the x direction, and mirrors along the y direction
template <typename ChannelDataType>
ChannelDataType get_pixel_bilinear_wrap_mirror( vector2f normalizedCoords, const std::vector<ChannelDataType>& channel,
                                                const size2 size ) {
    normalizedCoords.x = ( normalizedCoords.x + 1 ) * 0.5f;
    normalizedCoords.y = ( normalizedCoords.y + 1 ) * 0.5f;

    // First get these coords into the [0,1] range
    normalizedCoords.x = fmodf( normalizedCoords.x, 1 );
    normalizedCoords.y = fmodf( normalizedCoords.y, 2 );
    if( normalizedCoords.x < 0 )
        normalizedCoords.x += 1;
    if( normalizedCoords.y < 0 )
        normalizedCoords.y += 2;

    // Mirror the y
    if( normalizedCoords.y > 1 )
        normalizedCoords.y = 2 - normalizedCoords.y;

    normalizedCoords.x *= size.xsize;
    normalizedCoords.y *= size.ysize;

    // The + 2 - 2 is so that it rounds down
    int xMin = (int)( normalizedCoords.x + 2 ) - 2, yMin = (int)( normalizedCoords.y + 2 ) - 2;
    int xSize = 2, ySize = 2;
    float xAlpha[2], yAlpha[2];

    // Compute the blending values
    xAlpha[1] = normalizedCoords.x - xMin;
    xAlpha[0] = 1 - xAlpha[1];
    yAlpha[1] = normalizedCoords.y - yMin;
    yAlpha[0] = 1 - yAlpha[1];

    ChannelDataType result;
    if( xMin < 0 || xMin >= size.xsize - 1 || yMin < 0 || yMin >= size.ysize - 1 ) {
        // Slower lookup when wrapping must occur
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                // Get the wrapped coordinate within the bitmap
                result +=
                    xAlpha[x] * yAlpha[y] * get_pixel_wrap_duplicate( vector2( xMin + x, yMin + y ), channel, size );
            }
        }
    } else {
        // Faster indexed lookup
        int index = size.get_index( xMin, yMin );
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                result += xAlpha[x] * yAlpha[y] * channel[index + x];
            }
            index += size.xsize;
        }
    }
    return result;
}

// Gets a pixel value from normalized coordinates which are in [-1, 1]
// This method wraps along both the x and y directions
template <typename ChannelDataType>
ChannelDataType get_pixel_bilinear_wrap_wrap( vector2f normalizedCoords, const std::vector<ChannelDataType>& channel,
                                              const size2& size ) {
    normalizedCoords.x = ( normalizedCoords.x + 1 ) * 0.5f;
    normalizedCoords.y = ( normalizedCoords.y + 1 ) * 0.5f;

    // First get these coords into the [0,1] range
    normalizedCoords.x = fmodf( normalizedCoords.x, 1 );
    normalizedCoords.y = fmodf( normalizedCoords.y, 1 );
    if( normalizedCoords.x < 0 )
        normalizedCoords.x += 1;
    if( normalizedCoords.y < 0 )
        normalizedCoords.y += 1;

    normalizedCoords.x *= size.xsize;
    normalizedCoords.y *= size.ysize;

    // The + 2 - 2 is so that it rounds down
    int xMin = (int)( normalizedCoords.x + 2 ) - 2, yMin = (int)( normalizedCoords.y + 2 ) - 2;
    int xSize = 2, ySize = 2;
    float xAlpha[2], yAlpha[2];

    // Compute the blending values
    xAlpha[1] = normalizedCoords.x - xMin;
    xAlpha[0] = 1 - xAlpha[1];
    yAlpha[1] = normalizedCoords.y - yMin;
    yAlpha[0] = 1 - yAlpha[1];

    ChannelDataType result;
    if( xMin < 0 || xMin >= size.xsize - 1 || yMin < 0 || yMin >= size.ysize - 1 ) {
        // Slower lookup when wrapping must occur
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                // Get the wrapped coordinate within the bitmap
                result += xAlpha[x] * yAlpha[y] * get_pixel_wrap_wrap( vector2( xMin + x, yMin + y ), channel, size );
            }
        }
    } else {
        // Faster indexed lookup
        int index = size.get_index( xMin, yMin );
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                result += xAlpha[x] * yAlpha[y] * channel[index + x];
            }
            index += size.xsize;
        }
    }
    return result;
}

// Gets a pixel value from normalized coordinates which are in [-1, 1]
template <typename ChannelDataType>
const ChannelDataType& get_pixel_nearest_neighbor( float x, float y, const std::vector<ChannelDataType>& channel,
                                                   const size2& size ) {
    x = ( x + 1.0f ) * 0.5f * size.xsize;
    y = ( y + 1.0f ) * 0.5f * size.ysize;
    // This is the correct conversion to tile the pixels all the way from -1 to 1
    int ix = (int)( x + 0.5f ), iy = (int)( y + 0.5f );
    if( ix < 0 )
        ix = 0;
    if( ix >= size.xsize )
        ix = size.xsize - 1;
    if( iy < 0 )
        iy = 0;
    if( iy >= size.ysize )
        iy = size.ysize - 1;
    return channel[size.get_index( ix, iy )];
}

// Gets a pixel value from normalized coordinates which are in [-1, 1]
template <typename ChannelDataType>
const ChannelDataType& get_pixel_nearest_neighbor( vector2f normalizedCoords,
                                                   const std::vector<ChannelDataType>& channel, const size2& size ) {
    return get_pixel_nearest_neighbor( normalizedCoords.x, normalizedCoords.y, channel, size );
}

template <typename ChannelDataType>
const ChannelDataType& get_pixel_nearest_neighbor_int( int x, int y, const std::vector<ChannelDataType>& channel,
                                                       const size2& size ) {
    x = std::max( 0, x );
    y = std::max( 0, y );
    x = std::min( x, size.xsize - 1 );
    y = std::min( y, size.ysize - 1 );
    return channel[size.get_index( x, y )];
}

// Gets a pixel, wrapping along x and duplicating the edge along y.
template <typename ChannelDataType>
const ChannelDataType& get_pixel_wrap_duplicate( vector2 coords, const std::vector<ChannelDataType>& channel,
                                                 const size2& size ) {
    coords.x += size.xsize;
    coords.x = coords.x % size.xsize;

    if( coords.y < 0 )
        coords.y = 0;
    else if( coords.y >= size.ysize )
        coords.y = size.ysize - 1;
    return channel[size.get_index( coords )];
}

// Gets a pixel, wrapping along both x and y.
template <typename ChannelDataType>
const ChannelDataType& get_pixel_wrap_wrap( vector2 coords, const std::vector<ChannelDataType>& channel,
                                            const size2& size ) {
    // Takes care of negatives
    coords.x += size.xsize;
    coords.y += size.ysize;

    coords.x = coords.x % size.xsize;
    coords.y = coords.y % size.ysize;

    return channel[size.get_index( coords )];
}

template <typename ChannelDataType>
ChannelDataType get_pixel_bilinear_duplicate_duplicate( vector2f normalizedCoords,
                                                        const std::vector<ChannelDataType>& channel,
                                                        const size2& size ) {
    if( size.xsize == 0 || size.ysize == 0 )
        return ChannelDataType();

    normalizedCoords.x = ( normalizedCoords.x + 1 ) * 0.5f * size.xsize;
    normalizedCoords.y = ( normalizedCoords.y + 1 ) * 0.5f * size.ysize;
    normalizedCoords.x -= 0.5f;
    normalizedCoords.y -= 0.5f;

    // The + 2 - 2 is so that it rounds down
    int xMin = (int)( normalizedCoords.x + 2 ) - 2, yMin = (int)( normalizedCoords.y + 2 ) - 2;
    float xAlpha[2], yAlpha[2];

    // Compute the blending values
    xAlpha[1] = normalizedCoords.x - xMin;
    xAlpha[0] = 1 - xAlpha[1];
    yAlpha[1] = normalizedCoords.y - yMin;
    yAlpha[0] = 1 - yAlpha[1];

    if( xMin < 0 ) {
        if( yMin < 0 )
            return channel[0]; // Bottom left corner
        if( yMin + 1 >= size.ysize )
            return channel[( size.ysize - 1 ) * size.xsize]; // Top left corner

        int offset = yMin * size.xsize;
        return yAlpha[0] * channel[offset] + yAlpha[1] * channel[offset + size.xsize];
    } else if( xMin + 1 >= size.xsize ) {
        if( yMin < 0 )
            return channel[size.xsize - 1]; // Bottom right corner
        if( yMin + 1 >= size.ysize )
            return channel[size.xsize * size.ysize - 1]; // Top right corner

        int offset = ( yMin + 1 ) * size.xsize - 1;
        return yAlpha[0] * channel[offset] + yAlpha[1] * channel[offset + size.xsize];
    } else if( yMin < 0 ) {
        int offset = xMin;
        return xAlpha[0] * channel[offset] + xAlpha[1] * channel[offset + 1];
    } else if( yMin + 1 >= size.ysize ) {
        int offset = ( size.ysize - 1 ) * size.xsize + xMin;
        return xAlpha[0] * channel[offset] + xAlpha[1] * channel[offset + 1];
    }

    int index = size.get_index( xMin, yMin );
    return xAlpha[0] * yAlpha[0] * channel[index] + xAlpha[1] * yAlpha[0] * channel[index + 1] +
           xAlpha[0] * yAlpha[1] * channel[index + size.xsize] +
           xAlpha[1] * yAlpha[1] * channel[index + size.xsize + 1];
}

// Gets a pixel value from normalized coordinates which are in [-1, 1]
// This method also interpolates with black along the edges
template <typename ChannelDataType>
ChannelDataType get_pixel_bilinear( vector2f normalizedCoords, const std::vector<ChannelDataType>& channel,
                                    const size2& size, float* outWeightAccumulator = 0 ) {
    normalizedCoords.x = ( normalizedCoords.x + 1 ) * 0.5f * size.xsize;
    normalizedCoords.y = ( normalizedCoords.y + 1 ) * 0.5f * size.ysize;
    normalizedCoords.x -= 0.5f;
    normalizedCoords.y -= 0.5f;

    // The + 2 - 2 is so that it rounds down
    int xMin = (int)( normalizedCoords.x + 2 ) - 2, yMin = (int)( normalizedCoords.y + 2 ) - 2;
    int xSize = 2, ySize = 2;
    float xAlpha[2], yAlpha[2];

    // Compute the blending values
    xAlpha[1] = normalizedCoords.x - xMin;
    xAlpha[0] = 1 - xAlpha[1];
    yAlpha[1] = normalizedCoords.y - yMin;
    yAlpha[0] = 1 - yAlpha[1];

    if( xMin < 0 ) {
        if( xMin < -1 )
            return ChannelDataType();
        xMin = 0;
        --xSize;
        xAlpha[0] = xAlpha[1];
    }
    if( xMin >= size.xsize - 1 ) {
        if( xMin > size.xsize - 1 )
            return ChannelDataType();
        --xSize;
    }
    if( yMin < 0 ) {
        if( yMin < -1 )
            return ChannelDataType();
        yMin = 0;
        --ySize;
        yAlpha[0] = yAlpha[1];
    }
    if( yMin >= size.ysize - 1 ) {
        if( yMin > size.ysize - 1 )
            return ChannelDataType();
        --ySize;
    }

    ChannelDataType result = ChannelDataType();
    int index = size.get_index( xMin, yMin );
    float weight;
    // Do this test outside of the inner loops for a little bit of efficiency
    if( outWeightAccumulator == 0 ) {
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                weight = xAlpha[x] * yAlpha[y];
                result += weight * channel[index + x];
            }
            index += size.xsize;
        }
    } else {
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                weight = xAlpha[x] * yAlpha[y];
                result += weight * channel[index + x];
                *outWeightAccumulator += weight;
            }
            index += size.xsize;
        }
    }
    return result;
}

/**
 * Bicubic bSpline filter over a 4x4 grid.
 * This version duplicates the edges.
 */
template <typename ChannelDataType>
ChannelDataType get_pixel_bicubic_duplicate_duplicate( vector2f normalizedCoords,
                                                       const std::vector<ChannelDataType>& channel,
                                                       const size2& size ) {
    BSplineCubicFilter bSplineCubic;

    normalizedCoords.x = ( normalizedCoords.x + 1 ) * 0.5f * size.xsize;
    normalizedCoords.y = ( normalizedCoords.y + 1 ) * 0.5f * size.ysize;
    normalizedCoords.x -= 0.5f;
    normalizedCoords.y -= 0.5f;

    // The + 2 - 3 is so that it rounds down, and sticks it in bottom left corner of the filter
    int xMin = (int)( normalizedCoords.x + 2 ) - 3, yMin = (int)( normalizedCoords.y + 2 ) - 3;
    int xOffsets[4], yOffsets[4];
    float xWeights[4], yWeights[4];

    xOffsets[0] = math::clamp( xMin, 0, size.xsize - 1 );
    xOffsets[1] = math::clamp( xMin + 1, 0, size.xsize - 1 );
    xOffsets[2] = math::clamp( xMin + 2, 0, size.xsize - 1 );
    xOffsets[3] = math::clamp( xMin + 3, 0, size.xsize - 1 );

    yOffsets[0] = math::clamp( yMin, 0, size.ysize - 1 ) * size.xsize;
    yOffsets[1] = math::clamp( yMin + 1, 0, size.ysize - 1 ) * size.xsize;
    yOffsets[2] = math::clamp( yMin + 2, 0, size.ysize - 1 ) * size.xsize;
    yOffsets[3] = math::clamp( yMin + 3, 0, size.ysize - 1 ) * size.xsize;

    xWeights[0] = bSplineCubic( xMin - normalizedCoords.x );
    xWeights[1] = bSplineCubic( xMin + 1 - normalizedCoords.x );
    xWeights[2] = bSplineCubic( xMin + 2 - normalizedCoords.x );
    xWeights[3] = bSplineCubic( xMin + 3 - normalizedCoords.x );

    yWeights[0] = bSplineCubic( yMin - normalizedCoords.y );
    yWeights[1] = bSplineCubic( yMin + 1 - normalizedCoords.y );
    yWeights[2] = bSplineCubic( yMin + 2 - normalizedCoords.y );
    yWeights[3] = bSplineCubic( yMin + 3 - normalizedCoords.y );

    /*for(int y = 0; y < 4; ++y){
      for(int x = 0; x < 4; ++x)
        result += (xWeights[x] * yWeights[y]) * channel[xOffsets[x] + yOffsets[y]];
    }*/
    return ( xWeights[0] * yWeights[0] ) * channel[xOffsets[0] + yOffsets[0]] +
           ( xWeights[1] * yWeights[0] ) * channel[xOffsets[1] + yOffsets[0]] +
           ( xWeights[2] * yWeights[0] ) * channel[xOffsets[2] + yOffsets[0]] +
           ( xWeights[3] * yWeights[0] ) * channel[xOffsets[3] + yOffsets[0]] +

           ( xWeights[0] * yWeights[1] ) * channel[xOffsets[0] + yOffsets[1]] +
           ( xWeights[1] * yWeights[1] ) * channel[xOffsets[1] + yOffsets[1]] +
           ( xWeights[2] * yWeights[1] ) * channel[xOffsets[2] + yOffsets[1]] +
           ( xWeights[3] * yWeights[1] ) * channel[xOffsets[3] + yOffsets[1]] +

           ( xWeights[0] * yWeights[2] ) * channel[xOffsets[0] + yOffsets[2]] +
           ( xWeights[1] * yWeights[2] ) * channel[xOffsets[1] + yOffsets[2]] +
           ( xWeights[2] * yWeights[2] ) * channel[xOffsets[2] + yOffsets[2]] +
           ( xWeights[3] * yWeights[2] ) * channel[xOffsets[3] + yOffsets[2]] +

           ( xWeights[0] * yWeights[3] ) * channel[xOffsets[0] + yOffsets[3]] +
           ( xWeights[1] * yWeights[3] ) * channel[xOffsets[1] + yOffsets[3]] +
           ( xWeights[2] * yWeights[3] ) * channel[xOffsets[2] + yOffsets[3]] +
           ( xWeights[3] * yWeights[3] ) * channel[xOffsets[3] + yOffsets[3]];
}

// Gets a pixel value from normalized coordinates which are in [-1, 1]
// This method also interpolates with black along the edges
template <typename ChannelDataType>
ChannelDataType get_pixel_bicubic( vector2f normalizedCoords, const std::vector<ChannelDataType>& channel,
                                   const size2& size, float* outWeightAccumulator = 0 ) {
    BSplineCubicFilter bSplineCubic;

    normalizedCoords.x = ( normalizedCoords.x + 1 ) * 0.5f * size.xsize;
    normalizedCoords.y = ( normalizedCoords.y + 1 ) * 0.5f * size.ysize;
    normalizedCoords.x -= 0.5f;
    normalizedCoords.y -= 0.5f;

    /*
    A query point at x uses the surrounding
    sixteen sample points (@) for the bicubic. The bicubic filter
    gives a zero weight to everything outside +-2.

      |			|			|			|			|
    ----@-----------@-----------@-----------@----------------
      |			|			|			|			|
      |			|			|			|			|
      |			|			|			|			|
    ----@-----------@-----------@-----------@----------------
      |			|			|			|			|
      |			|	x		|			|			|
      |			|			|			|			|
    ----@-----------@-----------@-----------@----------------
      |			|			|			|			|
      |			|			|			|			|
      |			|			|			|			|
ymin----@-----------@-----------@-----------@----------------
      |			|			|			|			|
      |			|			|			|			|
      |			|			|			|			|
    ---------------------------------------------------------
      |			|			|			|			|
      xmin

    */

    // The + 2 - 3 is so that it rounds down
    int xMin = (int)( normalizedCoords.x + 2 ) - 3, yMin = (int)( normalizedCoords.y + 2 ) - 3;
    int xSize = 4, ySize = 4;

    if( xMin < 0 ) {
        if( xMin < -3 )
            return ChannelDataType();
        xSize += xMin;
        xMin = 0;
    }
    if( xMin >= size.xsize - 3 ) {
        if( xMin > size.xsize - 1 )
            return ChannelDataType();
        xSize = ( size.xsize - xMin );
    }
    if( yMin < 0 ) {
        if( yMin < -3 )
            return ChannelDataType();
        ySize += yMin;
        yMin = 0;
    }
    if( yMin >= size.ysize - 3 ) {
        if( yMin > size.ysize - 1 )
            return ChannelDataType();
        ySize = ( size.ysize - yMin );
    }

    float xWeight[4];
    float yWeight[4];
    for( int y = 0; y < ySize; ++y ) {
        float offset = ( yMin + y ) - normalizedCoords.y;
        yWeight[y] = bSplineCubic( offset );
    }
    for( int x = 0; x < xSize; ++x ) {
        float offset = ( xMin + x ) - normalizedCoords.x;
        xWeight[x] = bSplineCubic( offset );
    }

    ChannelDataType result = ChannelDataType();
    int index = size.get_index( xMin, yMin );
    float weight;
    // Do this test outside of the inner loops for a little bit of efficiency
    if( outWeightAccumulator == 0 ) {
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                weight = yWeight[y] * xWeight[x];
                result += weight * channel[index + x];
            }
            index += size.xsize;
        }
    } else {
        for( int y = 0; y < ySize; ++y ) {
            for( int x = 0; x < xSize; ++x ) {
                weight = yWeight[y] * xWeight[x];
                result += weight * channel[index + x];
                *outWeightAccumulator += weight;
            }
            index += size.xsize;
        }
    }

    return result;
}

template <typename ChannelDataType>
void get_pixel_gradient( vector2 coords, const std::vector<ChannelDataType>& channel, const size2& size,
                         float& outGradient, float& outGradientDirection ) {
    get_pixel_gradient( coords.x, coords.y, channel, size, outGradient, outGradientDirection );
}

template <typename ChannelDataType, typename AccumType>
void get_pixel_gradient( int x, int y, const std::vector<ChannelDataType>& channel, const size2& size,
                         AccumType& outGradient, float& outGradientDirection ) {
    // Apply a sobel filter to the pixel in question.

    // SobelX
    // [ -1  0  1 ]
    // [ -2  0  2 ]
    // [ -1  0  1 ]

    // SobelY
    // [ -1 -2 -1 ]
    // [  0  0  0 ]
    // [  1  2  1 ]

    AccumType sobelX = 0;
    AccumType sobelY = 0;
    ChannelDataType p;

    // Top Left
    p = get_pixel_nearest_neighbor_int( x - 1, y - 1, channel, size );
    sobelX += -p;
    sobelY += -p;

    // Top Middle
    p = get_pixel_nearest_neighbor_int( x, y - 1, channel, size );
    sobelY += -( p + p );

    // Top Right
    p = get_pixel_nearest_neighbor_int( x + 1, y - 1, channel, size );
    sobelX += p;
    sobelY += -p;

    // Middle Left
    p = get_pixel_nearest_neighbor_int( x - 1, y, channel, size );
    sobelX += -( p + p );

    // Middle Right
    p = get_pixel_nearest_neighbor_int( x + 1, y, channel, size );
    sobelX += p + p;

    // Bottom Left
    p = get_pixel_nearest_neighbor_int( x - 1, y + 1, channel, size );
    sobelX += -p;
    sobelY += p;

    // Bottom Middle
    p = get_pixel_nearest_neighbor_int( x, y + 1, channel, size );
    sobelY += p + p;

    // Bottom Right
    p = get_pixel_nearest_neighbor_int( x + 1, y + 1, channel, size );
    sobelX += p;
    sobelY += p;

    outGradient = fabs( sobelX ) + fabs( sobelY );
    outGradientDirection = (float)atan2( sobelX, sobelY );
}

template <typename ChannelDataType>
ChannelDataType get_x_partial_derivative( const vector2& coord, const std::vector<ChannelDataType>& channel,
                                          const size2& size ) {
    if( coord.x > 0 ) {
        if( coord.x < size.xsize - 1 ) {
            // Use the centered gradient to reduce bias
            ChannelDataType cPos = get_pixel( vector2( coord.x + 1, coord.y ), channel, size );
            ChannelDataType cNeg = get_pixel( vector2( coord.x - 1, coord.y ), channel, size );
            return 0.5f * ( cPos - cNeg );
        } else {
            // gradient at x==width-1
            return get_pixel( vector2( size.xsize - 1, coord.y ), channel, size ) -
                   get_pixel( vector2( size.xsize - 2, coord.y ), channel, size );
        }
    } else {
        // gradient at x==0
        return get_pixel( vector2( 1, coord.y ), channel, size ) - get_pixel( vector2( 0, coord.y ), channel, size );
    }
}

template <typename ChannelDataType>
ChannelDataType get_x_partial_derivative_squared( const vector2& coord, const std::vector<ChannelDataType>& channel,
                                                  const size2& size ) {
    if( coord.x > 0 ) {
        if( coord.x < size.xsize - 1 ) {
            // Use the product of the forward and backwards derivatives to reduce bias
            ChannelDataType cMid = get_pixel( coord, channel, size );
            ChannelDataType cPos = get_pixel( vector2( coord.x + 1, coord.y ), channel, size );
            ChannelDataType cNeg = get_pixel( vector2( coord.x - 1, coord.y ), channel, size );
            return ( cPos - cMid ) * ( cMid - cNeg );
        } else {
            // gradient at x==width-1
            ChannelDataType c = get_pixel( vector2( size.xsize - 1, coord.y ), channel, size ) -
                                get_pixel( vector2( size.xsize - 2, coord.y ), channel, size );
            return c * c;
        }
    } else {
        // gradient at x==0
        ChannelDataType c =
            get_pixel( vector2( 1, coord.y ), channel, size ) - get_pixel( vector2( 0, coord.y ), channel, size );
        return c * c;
    }
}

template <typename ChannelDataType>
ChannelDataType get_y_partial_derivative( const vector2& coord, const std::vector<ChannelDataType>& channel,
                                          const size2& size ) {
    if( coord.y > 0 ) {
        if( coord.y < size.ysize - 1 ) {
            // Use the centered gradient to reduce bias
            ChannelDataType cPos = get_pixel( vector2( coord.x, coord.y + 1 ), channel, size );
            ChannelDataType cNeg = get_pixel( vector2( coord.x, coord.y - 1 ), channel, size );
            return 0.5f * ( cPos - cNeg );
        } else {
            // gradient at y==height-1
            return get_pixel( vector2( coord.x, size.ysize - 1 ), channel, size ) -
                   get_pixel( vector2( coord.x, size.ysize - 2 ), channel, size );
        }
    } else {
        // gradient at y==0
        return get_pixel( vector2( coord.x, 1 ), channel, size ) - get_pixel( vector2( coord.x, 0 ), channel, size );
    }
}

template <typename ChannelDataType>
ChannelDataType get_y_partial_derivative_squared( const vector2& coord, const std::vector<ChannelDataType>& channel,
                                                  const size2& size ) {
    if( coord.y > 0 ) {
        if( coord.y < size.ysize - 1 ) {
            // Use the product of the forward and backwards derivatives to reduce bias
            ChannelDataType cMid = get_pixel( coord, channel, size );
            ChannelDataType cPos = get_pixel( vector2( coord.x, coord.y + 1 ), channel, size );
            ChannelDataType cNeg = get_pixel( vector2( coord.x, coord.y - 1 ), channel, size );
            return ( cPos - cMid ) * ( cMid - cNeg );
        } else {
            // gradient at y==height-1
            ChannelDataType c = get_pixel( vector2( coord.x, size.ysize - 1 ), channel, size ) -
                                get_pixel( vector2( coord.x, size.ysize - 2 ), channel, size );
            return c * c;
        }
    } else {
        // gradient at y==0
        ChannelDataType c =
            get_pixel( vector2( coord.x, 1 ), channel, size ) - get_pixel( vector2( coord.x, 0 ), channel, size );
        return c * c;
    }
}
} // namespace image
} // namespace graphics2d
} // namespace frantic
