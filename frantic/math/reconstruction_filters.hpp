// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frantic {
namespace math {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double sinc( double x ) {
    x *= M_PI;
    if( x != 0 )
        return ( sin( x ) / x );
    return ( 1.0 );
}

struct gaussian_filter {
    double m_sigma;
    double m_negHalfInvSigmaSq;
    double m_normalizer;
    gaussian_filter( double sigma )
        : m_sigma( sigma )
        , m_negHalfInvSigmaSq( -0.5 / ( sigma * sigma ) )
        , m_normalizer( 1.0 / sqrt( 2.0 * M_PI ) * m_sigma ) {}
    double radius() { return 4.0 * m_sigma; }
    double operator()( double t ) const { return m_normalizer * exp( t * t * m_negHalfInvSigmaSq ); }
    std::string name() { return "gaussian"; }
};

struct vanilla_filter {
    double radius() { return 1.0; }
    double operator()( double t ) const {
        /* f(t) = 2|t|^3 - 3|t|^2 + 1, -1 <= t <= 1 */
        if( t < 0.0 )
            t = -t;
        if( t < 1.0 )
            return ( ( 2.0 * t - 3.0 ) * t * t + 1.0 );
        return ( 0.0 );
    }
    std::string name() { return "vanilla"; }
};

struct box_filter {
    double radius() { return 0.5; }
    double operator()( double t ) const {
        if( ( t > -0.5 ) && ( t <= 0.5 ) )
            return ( 1.0 );
        return ( 0.0 );
    }
    std::string name() { return "box"; }
};

struct triangle_filter {
    double radius() { return 1.0; }
    double operator()( double t ) const {
        if( t < 0.0 )
            t = -t;
        if( t < 1.0 )
            return ( 1.0 - t );
        return ( 0.0 );
    }
    std::string name() { return "triangle"; }
};

struct bell_filter {
    double radius() { return 1.5; }
    double operator()( double t ) const {
        if( t < 0 )
            t = -t;

        if( t < .5 )
            return ( .75 - ( t * t ) );

        if( t < 1.5 ) {
            t = ( t - 1.5 );
            return ( .5 * ( t * t ) );
        }
        return ( 0.0 );
    }
    std::string name() { return "bell"; }
};

struct b_spline_filter {
    double radius() { return 2.0; }
    double operator()( double t ) const {
        double tt;

        if( t < 0 )
            t = -t;

        if( t < 1 ) {
            tt = t * t;
            return ( ( .5 * tt * t ) - tt + ( 2.0 / 3.0 ) );
        } else if( t < 2 ) {
            t = 2 - t;
            return ( ( 1.0 / 6.0 ) * ( t * t * t ) );
        }
        return ( 0.0 );
    }
    std::string name() { return "b_spline"; }
};

struct lanczos_3_filter {
    double radius() { return 3.0; }
    double operator()( double t ) const {
        if( t < 0 )
            t = -t;

        if( t < 3.0 )
            return ( sinc( t ) * sinc( t / 3.0 ) );

        return ( 0.0 );
    }
    std::string name() { return "lanczos 3"; }
};

struct mitchell_filter {
    double radius() { return 2.0; }
    double operator()( double t ) const {
        double tt;
        double B = 1.0 / 3.0;
        double C = B;

        tt = t * t;

        if( t < 0 )
            t = -t;

        if( t < 1.0 ) {
            t = ( ( ( 12.0 - 9.0 * B - 6.0 * C ) * ( t * tt ) ) + ( ( -18.0 + 12.0 * B + 6.0 * C ) * tt ) +
                  ( 6.0 - 2 * B ) );
            return t / 6.0;
        } else if( t < 2.0 ) {
            t = ( ( ( -1.0 * B - 6.0 * C ) * ( t * tt ) ) + ( ( 6.0 * B + 30.0 * C ) * tt ) +
                  ( ( -12.0 * B - 48.0 * C ) * t ) + ( 8.0 * B + 24 * C ) );
            return t / 6.0;
        }
        return ( 0.0 );
    }
    std::string name() { return "mitchell"; }
};

} // namespace math
} // namespace frantic
