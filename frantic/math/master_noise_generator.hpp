// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>

#define FASTFLOOR( x ) ( ( ( x ) > 0 ) ? ( (int)x ) : ( (int)x - 1 ) )

namespace frantic {
namespace math {

/**
 * An abstract base class for noise generators.
 * Maintains information about the number of octaves and persistence.
 * Handles combining noise values of multiple octaves.
 *
 * The final noise value is calculated as
 *
 * for octave=startingOctave to startingOctave + numOctaves
 *     noiseSum += noise(x * lacunarity ^ octave) * persistence ^ (octave - startingOctave)
 *
 * where noise is a function supplied by a child of this class.
 **/
class master_noise_generator {
  protected:
    int m_octaves;
    float m_persistence;
    float m_lacunarity;
    int m_startingOctave;

  public:
    master_noise_generator( int octaves = 5, float persistence = 1.0f )
        : m_octaves( octaves )
        , m_persistence( persistence )
        , m_lacunarity( 2.f )
        , m_startingOctave( 0 ) {}

    virtual ~master_noise_generator() {}

    int get_octaves() { return m_octaves; }
    float get_persistence() { return m_persistence; }
    int get_starting_octave() { return m_startingOctave; }
    float get_lacunarity() { return m_lacunarity; }

    /**
     * Sets the number of octaves to obtain noise values from
     * @param octaves - the number of octaves
     */
    void set_octaves( int octaves ) { m_octaves = std::max<int>( 1, octaves ); }

    /**
     * Sets the persistence which controls the influence given to each additional octave.
     * Usually but not limited to between 0 and 1.
     * @param persistence - the persistence
     */
    void set_persistence( float persistence ) { m_persistence = persistence; }

    /**
     * Sets the octave to starting octave.  This is normally 0.
     * @param octave - the starting octave
     */
    void set_starting_octave( int octave ) { m_startingOctave = std::max<int>( 0, octave ); }

    /**
     *  Set the lacunarity.
     *
     * @param lacunarity the multiplicative increase between successive
     *		frequencies.  Use 2.0 for octaves.
     */
    void set_lacunarity( float lacunarity ) { m_lacunarity = lacunarity; }

    /**
     * Generates a noise value for a 1 dimensional input
     * @param x - x value
     */
    float noise( float x ) {
        float total = 0.0f;
        float freq = powf( m_lacunarity, static_cast<float>( m_startingOctave ) );

        for( int i = m_startingOctave; i < ( m_startingOctave + m_octaves ); i++ ) {
            float amplitude = powf( m_persistence, static_cast<float>( i - m_startingOctave ) );

            total += get_noise( x * freq ) * amplitude;

            freq *= m_lacunarity;
        }

        return total;
    }

    /**
     * Generates a noise value for a 2 dimensional input
     * @param x - x value
     * @param y - y value
     */
    float noise( float x, float y ) {
        float total = 0.0f;
        float freq = powf( m_lacunarity, static_cast<float>( m_startingOctave ) );

        for( int i = m_startingOctave; i < ( m_startingOctave + m_octaves ); i++ ) {
            float amplitude = powf( m_persistence, static_cast<float>( i - m_startingOctave ) );

            total += get_noise( x * freq, y * freq ) * amplitude;

            freq *= m_lacunarity;
        }

        return total;
    }

    /**
     * Generates a noise value for a 3 dimensional input
     * @param x - x value
     * @param y - y value
     * @param z - z value
     */
    float noise( float x, float y, float z ) {
        float total = 0.0f;
        float freq = powf( m_lacunarity, static_cast<float>( m_startingOctave ) );

        for( int i = m_startingOctave; i < ( m_startingOctave + m_octaves ); i++ ) {
            float amplitude = powf( m_persistence, static_cast<float>( i - m_startingOctave ) );

            total += get_noise( x * freq, y * freq, z * freq ) * amplitude;

            freq *= m_lacunarity;
        }

        return total;
    }

    float turbulence( float x, float y, float z ) {
        float result = 0;
        float amplitude = 1.f;

        if( m_startingOctave > 0 ) {
            float frequency = 1.f;
            for( int i = 0; i < m_startingOctave; ++i ) {
                amplitude *= m_persistence;
                frequency *= m_lacunarity;
            }
            x *= frequency;
            y *= frequency;
            z *= frequency;
        }
        for( int i = 0; i < m_octaves; ++i ) {
            result += amplitude * std::abs( get_noise( x, y, z ) );
            amplitude *= m_persistence;
            x *= m_lacunarity;
            y *= m_lacunarity;
            z *= m_lacunarity;
        }
        return result;
    }

  protected:
    virtual float get_noise( float x ) = 0;
    virtual float get_noise( float x, float y ) = 0;
    virtual float get_noise( float x, float y, float z ) = 0;
};
} // namespace math
} // namespace frantic
