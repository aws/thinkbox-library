// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cassert>

// File: utils.hpp
// Contains a set of useful mathematical values and tools.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>

#include <frantic/math/utils.hpp>

namespace frantic {
namespace math {
namespace filters {

inline void gaussian_filter( float sigma, std::vector<float>& outFilter ) {
    if( sigma > 0 ) {
        // compute the kernel out to about 3 standard deviations
        int middle = (int)ceilf( 3 * sigma );
        // Make sure the filter width is at least 3.
        if( middle == 0 )
            middle = 1;
        int newVectorSize = 2 * middle + 1;
        outFilter.resize( newVectorSize );
        // fill the filter with a gaussian kernel

        // Instead of using the theoretical normalization constant, we compute the sum of all the
        // values so that the result is actually normalized
        //			float normalizationConstant = 0.3989422804014326779f / sigma; // (1/(sqrt(2*M_PI)*sigma)

        float sumOfAllValues = 1.f;
        outFilter[middle] = 1.f;

        for( int i = 0; i < middle; ++i ) {
            float value = exp( -frantic::math::square( middle - i ) / ( 2 * sigma * sigma ) );
            outFilter[i] = value;
            outFilter[outFilter.size() - 1 - i] = value;
            sumOfAllValues += 2 * value;
        }

        // normalize all values to have a sum of 1.0
        for( int normalize = 0; normalize < newVectorSize; ++normalize )
            outFilter[normalize] /= sumOfAllValues;
    } else {
        outFilter.resize( 1 );
        outFilter[0] = 1;
    }
}

inline void sobel_perpendicular_filter( std::vector<float>& outFilter ) {
    outFilter.resize( 3 );
    outFilter[0] = -1;
    outFilter[1] = 0;
    outFilter[2] = 1;
}

inline void sobel_parallel_filter( std::vector<float>& outFilter ) {
    outFilter.resize( 3 );
    outFilter[0] = 1;
    outFilter[1] = 2;
    outFilter[2] = 1;
}

// This is a 2D laplacian filter, stored in a square buffer of size by size.
inline void laplacian_filter( int size, std::vector<float>& outFilter ) {
    if( size % 2 != 1 )
        throw std::runtime_error( "image_filter_utils::laplacian_filter: Filter size must be odd." );

    outFilter.resize( size * size, -1.0f );
    outFilter[outFilter.size() / 2] = (float)outFilter.size() - 1;
}

/////////////////////////////////////////
// Functions to apply the filters to arrays
/////////////////////////////////////////

template <class ElementType>
inline void apply_filter( const std::vector<float>& filter, std::vector<ElementType>& outMap ) {
    // Get the middle index, which lands on the pixel being processed (Assuming an odd filter size)
    int middle = int( filter.size() ) / 2, mapSize = (int)outMap.size();
    // Create a buffer for output
    std::vector<ElementType> tempBuffer( outMap.size() );
    for( int x = 0; x < mapSize; ++x ) {
        ElementType mapValue = ElementType();
        for( int f = 0; f < int( filter.size() ); ++f ) {
            int xf = x + f - middle;
            if( xf < 0 )
                xf = 0;
            if( xf >= mapSize )
                xf = mapSize - 1;
            mapValue += filter[f] * outMap[xf];
        }
        tempBuffer[x] = mapValue;
    }

    // Put the filtered version back into the parameter vector
    outMap.swap( tempBuffer );
}

template <class ElementType>
inline void apply_gaussian_filter( float sigma, std::vector<ElementType>& outMap ) {
    if( sigma > 0 ) {
        std::vector<float> filter;
        gaussian_filter( sigma, filter );
        apply_filter( filter, outMap );
    }
}

} // namespace filters
} // namespace math
} // namespace frantic
