// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/math/differentiable_noise_generator.hpp>
#include <frantic/math/utils.hpp>

namespace frantic {
namespace math {

/**
 * This class provides a Perlin noise function.
 *
 * This class can be used to generate random or pseudorandom noise in
 * 1,2,3, or 4 dimensions.
 *
 * It works by picking a gradient at each of points surronding the input point.
 * Therefore the number of gradients used will be 2^n where n is the number of dimensions.
 * Using the values at each of these points we can perform 2^n - 1 lerps to obtain the noise value.
 */
class perlin_noise_generator : public differentiable_noise_generator {
  private:
    int m_permCount;
    std::vector<int> m_permutations;

  public:
    perlin_noise_generator( int octaves = 5, float persistence = 0.5f )
        : differentiable_noise_generator( octaves, persistence ) {
        generate_permutations( 256 );
    }

    /**
     * A 4 dimensional version of the noise function
     * @param x - x value
     * @param y - y value
     * @param z - z value
     * @param t - t value
     */
    float noise4d( float x, float y, float z, float t );

    /**
     * A 4 dimensional noise function which also calculates
     * partial derivatives of each input
     * @param x - x value
     * @param y - y value
     * @param z - z value
     * @param t - t value
     * @param dx - pointer to store x derivative
     * @param dy - pointer to store y derivative
     * @param dz - pointer to store z derivative
     * @param dt - pointer to store t derivative
     */
    float dnoise4d( float x, float y, float z, float t, float* dx, float* dy, float* dz, float* dt );

    /**
     * Regenerate the set of pseudorandom permutations used to create the noise
     * The count will be rounded up to the nearest power of 2
     * Different seeds can be used to generate different permutations
     *
     * @param count - the number of permuations to generate
     * @param seed - a seed for rand
     */
    void generate_permutations( unsigned int count, int seed = 2 );

  private:
    float lerp( float m, float a, float b );

    float grad( int hash );
    float grad( int hash, float x, float y, float z, float t );

    float dfade( float x );
    float fade( float x );

    float dot( int* grad, float x, float y );
    float dot( int* grad, float x, float y, float z );
    float dot( int* grad, float x, float y, float z, float t );

    float get_noise( float x );
    float get_noise( float x, float y );
    float get_noise( float x, float y, float z );

    float get_dnoise( float x, float* dx );
    float get_dnoise( float x, float y, float* dx, float* dy );
    float get_dnoise( float x, float y, float z, float* dx, float* dy, float* dz );

    float pnoise( float x, float* dx = 0 );
    float pnoise( float x, float y, float* dx = 0, float* dy = 0 );
    float pnoise( float x, float y, float z, float* dx = 0, float* dy = 0, float* dz = 0 );
    float pnoise( float x, float y, float z, float t, float* dx = 0, float* dy = 0, float* dz = 0, float* dt = 0 );
};
} // namespace math
} // namespace frantic
