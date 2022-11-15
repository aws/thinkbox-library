// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics2d/size2f.hpp>
#include <vector>

#include <frantic/math/filters.hpp>

namespace frantic {
namespace graphics2d {
namespace image {

// A collection of image filter utilities that operate on a vector
class filter {
  public:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    // Apply filters to vectors

    // Applies a filter to the image in the x direction
    template <class ElementType>
    static void apply_filter_x( const std::vector<float>& filter, std::vector<ElementType>& outMap,
                                const size2& mapSize ) {
        // Get the middle index, which lands on the pixel being processed
        int middle = int( filter.size() ) / 2;
        // Allocate a temporary scanline buffer
        std::vector<ElementType> lineBuffer( mapSize.xsize );
        for( int y = 0; y < mapSize.ysize; ++y ) {
            int lineOffset = y * mapSize.xsize;
            // copy the scanline into the scanline buffer
            for( int i = 0; i < mapSize.xsize; ++i )
                lineBuffer[i] = outMap[lineOffset + i];
            for( int x = 0; x < mapSize.xsize; ++x ) {
                ElementType mapValue = ElementType();
                for( int f = 0; f < int( filter.size() ); ++f ) {
                    int xf = x + f - middle;
                    if( xf < 0 )
                        xf = 0;
                    if( xf >= mapSize.xsize )
                        xf = mapSize.xsize - 1;
                    mapValue += filter[f] * lineBuffer[xf];
                }
                outMap[lineOffset + x] = mapValue;
            }
        }
    }

    // Applies a filter to the image in the y direction
    template <class ElementType>
    static void apply_filter_y( const std::vector<float>& filter, std::vector<ElementType>& outMap,
                                const size2& mapSize ) {
        // Get the middle index, which lands on the pixel being processed
        int middle = int( filter.size() ) / 2;
        // Allocate a temporary vertical line buffer
        std::vector<ElementType> lineBuffer( mapSize.ysize );
        for( int x = 0; x < mapSize.xsize; ++x ) {
            // copy the vertical line into the line buffer
            for( int i = 0; i < mapSize.ysize; ++i )
                lineBuffer[i] = outMap[x + i * mapSize.xsize];
            for( int y = 0; y < mapSize.ysize; ++y ) {
                ElementType mapValue = ElementType(); // Color always defaults to zero
                for( int f = 0; f < int( filter.size() ); ++f ) {
                    int yf = y + f - middle;
                    if( yf < 0 )
                        yf = 0;
                    if( yf >= mapSize.ysize )
                        yf = mapSize.ysize - 1;
                    mapValue += filter[f] * lineBuffer[yf];
                }
                outMap[x + y * mapSize.xsize] = mapValue;
            }
        }
    }

    // Applies a separable filter to the image, by applying it once along x and
    // once along y.
    template <class ElementType>
    static void apply_separable_filter( const std::vector<float>& filter, std::vector<ElementType>& outMap,
                                        const size2& mapSize ) {
        apply_filter_x( filter, outMap, mapSize );
        apply_filter_y( filter, outMap, mapSize );
    }

    // Applies a square filter to the image
    template <class ElementType>
    static void apply_filter( const std::vector<float>& filter, const unsigned int filterSize,
                              std::vector<ElementType>& outMap, const size2& mapSize ) {

        if( filterSize % 2 != 1 )
            throw std::runtime_error( "image_filter_utils.apply_filter: Filter size must be odd." );

        int fSize = filterSize / 2;

        // Get a copy of the image
        std::vector<ElementType> mapBuffer;
        mapBuffer = outMap;

        for( int y = 0; y < mapSize.ysize; y++ ) {
            for( int x = 0; x < mapSize.xsize; x++ ) {
                ElementType mapValue = ElementType(); // Color always defaults to zero

                for( int fy = std::max( 0, y - fSize ); fy <= std::min( mapSize.ysize - 1, y + fSize ); fy++ ) {
                    for( int fx = std::max( 0, x - fSize ); fx <= std::min( mapSize.xsize - 1, x + fSize ); fx++ ) {
                        // fx, fy is in bounds
                        mapValue += filter[( fy - y + fSize ) * filterSize + ( fx - x + fSize )] *
                                    mapBuffer[fy * mapSize.xsize + fx];
                    }
                }
                outMap[x + y * mapSize.xsize] = mapValue;
            }
        }
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    // Common functions with filters

    // Applies a gaussian blur with standard deviation sigma to the image.
    // TODO: Applying multiple small gaussian blurs makes a bigger gaussian blur.  We
    //       can optimize this further by taking advantage of that fact.  Just be careful
    //       that you use a bigger scanlineBuffer in the x and y directions to prevent
    //       the truncation of intermediate results.  On second thought, this does not
    //       appear to provide any speedup whatsoever.
    template <class ElementType>
    static void apply_gaussian( float sigma, std::vector<ElementType>& outMap, const size2& mapSize ) {
        if( sigma > 0 ) {
            std::vector<float> filter;
            frantic::math::filters::gaussian_filter( sigma, filter );
            apply_separable_filter( filter, outMap, mapSize );
        }
    }

    template <class ElementType>
    static void apply_laplacian( int size, std::vector<ElementType>& outMap, const size2& mapSize ) {
        std::vector<float> filter;
        frantic::math::filters::laplacian_filter( size, filter );
        apply_filter( filter, size, outMap, mapSize );
    }

    template <class ElementType>
    static void apply_sobel( std::vector<ElementType>& outMap, const size2& mapSize ) {
        std::vector<float> filterX;
        std::vector<float> filterY;
    }

    static void apply_non_maximal_suppression( const std::vector<float>& inMap, const std::vector<float>& inDirection,
                                               std::vector<float>& outMap, const size2& mapSize,
                                               const float& angleThreshold = M_PI ) {
        // Create three filters for the four different angles we recognize
        std::vector<float> f0( 9, 1 );
        std::vector<float> f45( 9, 1 );
        std::vector<float> f90( 9, 1 );
        std::vector<float> f135( 9, 1 );

        // Horizontal
        f0[1] = 0;
        f0[7] = 0;
        // Upward sloping diagonal
        f45[0] = 0;
        f45[8] = 0;
        // Vertical
        f90[3] = 0;
        f90[5] = 0;
        // Downward sloping diagonal
        f135[2] = 0;
        f135[6] = 0;

        for( int y = 0; y < mapSize.ysize; y++ ) {
            for( int x = 0; x < mapSize.xsize; x++ ) {

                float val = inMap[y * mapSize.xsize + x];
                float direction = inDirection[y * mapSize.xsize + x];
                std::vector<float>* filter = 0;
                if( direction < 0.3926990f )
                    filter = &f0;
                else if( direction < 1.1780972f )
                    filter = &f45;
                else if( direction < 1.9634954f )
                    filter = &f90;
                else if( direction < 2.7488935f )
                    filter = &f135;
                else
                    filter = &f0;

                int index = 0;
                for( int fy = -1; fy <= 1; fy++ ) {
                    for( int fx = -1; fx <= 1; fx++ ) {
                        if( fx + x > 0 && fx + x < mapSize.xsize && fy + y > 0 && fy + y < mapSize.ysize ) {
                            float curDir = inDirection[( fy + y ) * mapSize.xsize + ( fx + x )];
                            if( ( *filter )[index] != 1 &&
                                inMap[( fy + y ) * mapSize.xsize + ( fx + x )] > val
                                // Threshold the angle (if the angles are "close")
                                && ( curDir > direction - angleThreshold && curDir < direction + angleThreshold ) )
                                val = 0;
                        }
                        index++;
                    }
                }
                outMap[y * mapSize.xsize + x] = val;
            }
        }
    }
};
} // namespace image
} // namespace graphics2d
} // namespace frantic
