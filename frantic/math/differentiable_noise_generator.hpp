// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/math/master_noise_generator.hpp>

namespace frantic {
namespace math {
/**
 * An abstract class for noise functions that are differentiable.
 * It calculates and returns the noise value just as a regular noise generator
 * but additionally it will store the partial derivatives for each
 * input at an supplied initialized pointer.
 *
 * When there are multiple octaves the derivatives are summed the same way
 * noise values are in the master noise generator.
 *
 */

class differentiable_noise_generator : public master_noise_generator {
  public:
    differentiable_noise_generator( int octaves, float persistence )
        : master_noise_generator( octaves, persistence ) {}

    /**
     * The dnoise function for 1 dimension
     * @param x - x input
     * @param dx - pointer to store x derivative
     */
    float dnoise( float x, float* dx ) {
        float total = 0.0f;
        float dx_total = 0.0f;
        float freq = powf( m_lacunarity, static_cast<float>( m_startingOctave ) );

        for( int i = m_startingOctave; i < ( m_startingOctave + m_octaves ); i++ ) {
            float amplitude = powf( m_persistence, static_cast<float>( i - m_startingOctave ) );

            total += get_dnoise( x * freq, dx ) * amplitude;

            if( dx != 0 )
                dx_total += *dx * amplitude;

            freq *= m_lacunarity;
        }

        if( dx != 0 )
            *dx = dx_total;

        return total;
    }

    /**
     * The dnoise function for 2 dimensions
     * @param x - x input
     * @param y - y input
     * @param dx - pointer to store x derivative
     * @param dy - pointer to store y derivative
     */
    float dnoise( float x, float y, float* dx, float* dy ) {
        float total = 0.0f;
        float dx_total = 0.0f;
        float dy_total = 0.0f;
        float freq = powf( m_lacunarity, static_cast<float>( m_startingOctave ) );

        for( int i = m_startingOctave; i < ( m_startingOctave + m_octaves ); i++ ) {
            float amplitude = powf( m_persistence, static_cast<float>( i - m_startingOctave ) );

            total += get_dnoise( x * freq, y * freq, dx, dy ) * amplitude;

            if( dx != 0 && dy != 0 ) {
                dx_total += *dx * amplitude;
                dy_total += *dy * amplitude;
            }

            freq *= m_lacunarity;
        }

        if( dx != 0 && dy != 0 ) {
            *dx = dx_total;
            *dy = dy_total;
        }

        return total;
    }

    /**
     * The dnoise function for 3 dimensions
     * @param x - x input
     * @param y - y input
     * @param z - z input
     * @param dx - pointer to store x derivative
     * @param dy - pointer to store y derivative
     * @param dz - pointer to store z derivative
     */
    float dnoise( float x, float y, float z, float* dx, float* dy, float* dz ) {
        float total = 0.0f;
        float dx_total = 0.0f;
        float dy_total = 0.0f;
        float dz_total = 0.0f;
        float freq = powf( m_lacunarity, static_cast<float>( m_startingOctave ) );

        for( int i = m_startingOctave; i < ( m_startingOctave + m_octaves ); i++ ) {
            float amplitude = powf( m_persistence, static_cast<float>( i - m_startingOctave ) );

            total += get_dnoise( x * freq, y * freq, z * freq, dx, dy, dz ) * amplitude;

            if( dx != 0 && dy != 0 && dz != 0 ) {
                dx_total += *dx * amplitude;
                dy_total += *dy * amplitude;
                dz_total += *dz * amplitude;
            }

            freq *= m_lacunarity;
        }

        if( dx != 0 && dy != 0 && dz != 0 ) {
            *dx = dx_total;
            *dy = dy_total;
            *dz = dz_total;
        }

        return total;
    }

  protected:
    virtual float get_dnoise( float x, float* dx ) = 0;
    virtual float get_dnoise( float x, float y, float* dx, float* dy ) = 0;
    virtual float get_dnoise( float x, float y, float z, float* dx, float* dy, float* dz ) = 0;
};
} // namespace math
} // namespace frantic
