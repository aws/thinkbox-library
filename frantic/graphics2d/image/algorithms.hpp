// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stack>

namespace frantic {
namespace graphics2d {
namespace image {

// This algorithm doe
template <typename ColorType>
void flood_fill( vector2 start, const ColorType& fill, const ColorType& boundary, bool /*overwriteBoundary*/,
                 std::vector<ColorType>& data, const size2 size ) {
    using std::stack;

    if( start.x < 0 || start.y < 0 || start.x >= size.xsize || start.y >= size.ysize )
        return;

    stack<vector2> scanStart;
    scanStart.push( start );

    while( !scanStart.empty() ) {
        vector2 curScan = scanStart.top();
        scanStart.pop();

        // Fill the scanline
        int x = curScan.x;

        ColorType c = data[curScan.y * size.xsize + curScan.x];

        bool originalFill = ( c != boundary );

        // Left
        while( x >= 0 && c != boundary ) {
            data[curScan.y * size.xsize + x] = fill;

            // Whenever we hit a spot where the upper pixel is "free" but the
            // pixel to the left of it is not, push the location
            if( curScan.y > 0 ) {
                ColorType cNear = data[( curScan.y - 1 ) * size.xsize + x];
                if( cNear != boundary && cNear != fill ) {
                    if( x == 0 )
                        scanStart.push( vector2( 0, curScan.y - 1 ) );
                    else if( data[( curScan.y - 1 ) * size.xsize + x - 1] == boundary )
                        scanStart.push( vector2( x, curScan.y - 1 ) );
                }
            }

            if( curScan.y < size.ysize - 1 ) {
                ColorType cNear = data[( curScan.y + 1 ) * size.xsize + x];
                if( cNear != boundary && cNear != fill ) {
                    if( x == 0 )
                        scanStart.push( vector2( 0, curScan.y + 1 ) );
                    else if( data[( curScan.y + 1 ) * size.xsize + x - 1] == boundary )
                        scanStart.push( vector2( x, curScan.y + 1 ) );
                }
            }

            x--;
            if( x >= 0 )
                c = data[curScan.y * size.xsize + x];
        }

        if( x > 0 && x < size.xsize - 1 && c == boundary ) {
            if( curScan.y > 0 ) {
                ColorType cNear = data[( curScan.y - 1 ) * size.xsize + ( x + 1 )];
                if( cNear != boundary && cNear != fill ) {
                    scanStart.push( vector2( x + 1, curScan.y - 1 ) );
                }
            }

            if( curScan.y < size.ysize - 1 ) {
                ColorType cNear = data[( curScan.y + 1 ) * size.xsize + ( x + 1 )];
                if( cNear != boundary && cNear != fill ) {
                    scanStart.push( vector2( x + 1, curScan.y + 1 ) );
                }
            }
        }

        // Right side
        x = curScan.x;
        c = data[curScan.y * size.xsize + curScan.x];
        while( ( x == curScan.x && originalFill ) || ( x < size.xsize && c != boundary ) ) {
            data[curScan.y * size.xsize + x] = fill;
            // Whenever we hit a spot where the upper pixel is "free" but the
            // pixel to the left of it is not, push the location
            if( curScan.y > 0 ) {
                ColorType cNear = data[( curScan.y - 1 ) * size.xsize + x];
                if( cNear != boundary && cNear != fill ) {
                    if( x == 0 )
                        scanStart.push( vector2( 0, curScan.y - 1 ) );
                    else if( data[( curScan.y - 1 ) * size.xsize + x - 1] == boundary )
                        scanStart.push( vector2( x, curScan.y - 1 ) );
                }
            }

            if( curScan.y < size.ysize - 1 ) {
                ColorType cNear = data[( curScan.y + 1 ) * size.xsize + x];
                if( cNear != boundary && cNear != fill ) {
                    if( x == 0 )
                        scanStart.push( vector2( 0, curScan.y + 1 ) );
                    else if( data[( curScan.y + 1 ) * size.xsize + x - 1] == boundary )
                        scanStart.push( vector2( x, curScan.y + 1 ) );
                }
            }

            x++;
            if( x < size.xsize )
                c = data[curScan.y * size.xsize + x];
        }
    }
}

template <class Image>
void resize( Image& src, const size2& size ) {
    Image dest( size );
    dest.set_channels( src.channel_names() );

    int channelCount = src.channel_count();
    for( int i = 0; i < channelCount; ++i ) {
        resize_up( src.channel_data( i ), src.size(), dest.channel_data( src.channel_name( i ) ), dest.size() );
    }

    src.swap( dest );
}

template <typename ColorType>
void resize_up( std::vector<ColorType>& dest, const size2 destSize, const std::vector<ColorType>& src,
                const size2 srcSize ) {
    int yAdd = srcSize.ysize + srcSize.ysize;
    int xAdd = srcSize.xsize + srcSize.xsize;
    int xAccum = destSize.xsize;
    int yAccum = destSize.ysize;
    int xAccumMax = destSize.xsize + destSize.xsize;
    int yAccumMax = destSize.ysize + destSize.ysize;

    int srcX = 0;
    int srcY = 0;
    int destIndex = 0;

    ColorType oldTop = src[0];
    ColorType oldBottom = src[0];
    ColorType topLeft = src[0];
    ColorType topRight = src[0];
    ColorType bottomLeft = src[0];
    ColorType bottomRight = src[0];

    ColorType topLeftWAlpha = src[0];
    ColorType topRightWAlpha = src[0];
    ColorType bottomLeftWAlpha = src[0];
    ColorType bottomRightWAlpha = src[0];

    float alphaRight = 0;
    float alphaLeft = 0;
    float alphaBottom = 0;
    float alphaTop = 0;

    // For each Y in the image
    for( int y = 0; y < destSize.ysize; ++y, yAccum += yAdd ) {
        if( yAccum > yAccumMax ) {
            yAccum -= yAccumMax;
            oldTop = oldBottom;
            srcY++;
            oldBottom = src[std::min( srcY, srcSize.ysize - 1 ) * srcSize.xsize];
        }

        topRight = topLeft = oldTop;
        bottomRight = bottomLeft = oldBottom;

        alphaBottom = (float)yAccum / yAccumMax;
        alphaTop = 1 - alphaBottom;

        topLeftWAlpha = topLeft * alphaTop;
        topRightWAlpha = topRight * alphaTop;
        bottomLeftWAlpha = bottomLeft * alphaBottom;
        bottomRightWAlpha = bottomRight * alphaBottom;

        srcX = 0;
        xAccum = destSize.xsize;
        for( int x = 0; x < destSize.xsize; ++x, xAccum += xAdd ) {
            if( xAccum > xAccumMax ) {
                xAccum -= xAccumMax;
                ++srcX;
                topLeft = topRight;
                topRight = src[std::min( srcX, srcSize.xsize - 1 ) + std::max( srcY - 1, 0 ) * srcSize.xsize];
                bottomLeft = bottomRight;
                bottomRight =
                    src[std::min( srcX, srcSize.xsize - 1 ) + std::min( srcY, srcSize.ysize - 1 ) * srcSize.xsize];

                topLeftWAlpha = topLeft * alphaTop;
                topRightWAlpha = topRight * alphaTop;
                bottomLeftWAlpha = bottomLeft * alphaBottom;
                bottomRightWAlpha = bottomRight * alphaBottom;
            }

            alphaRight = (float)xAccum / xAccumMax;
            alphaLeft = 1 - alphaRight;

            dest[destIndex++] = alphaLeft * topLeftWAlpha + alphaRight * topRightWAlpha + alphaLeft * bottomLeftWAlpha +
                                alphaRight * bottomRightWAlpha;
        }
    }
}
} // namespace image
} // namespace graphics2d
} // namespace frantic
