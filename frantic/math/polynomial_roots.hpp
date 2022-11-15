// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stdexcept>

#ifndef M_PI
#define M_PI 3.141592653589
#endif

namespace frantic {
namespace math {
namespace polynomial_roots {

// Computes the roots of the quadratic formula a*x^2 + b*x + c = 0
// Returns the number of roots
template <class T>
inline int get_quadratic_roots( T a, T b, T c, T& outRoot1, T& outRoot2 ) {
    if( a == 0 )
        throw std::runtime_error( "frantic::math::polynomial_roots::get_quadratic_roots: Invalid parameter: a = 0." );

    T discriminant = b * b - 4 * a * c;
    if( discriminant < 0 ) {
        outRoot1 = outRoot2 = 0;
        return 0; // No roots;
    }

    if( discriminant == 0 ) {
        // For stability reasons, choose the largest denominator
        if( fabs( 2 * a ) > fabs( b ) ) {
            outRoot1 = outRoot2 = -b / ( 2 * a );
        } else {
            outRoot1 = outRoot2 = ( 2 * c ) / ( -b );
        }
        return 1;
    } else {
        T sqrtDiscriminant = sqrt( discriminant );

        // The reason we discriminate between positive and negative values of B is
        // because we want B to be the same sign as sqrtDiscriminant. This is so that
        // we do not have rounding errors when B is close to sqrtDiscriminant.

        if( b >= 0 ) {
            // Get the stable roots for positive b
            outRoot1 = ( 2 * c ) / ( -b - sqrtDiscriminant );
            outRoot2 = ( -b - sqrtDiscriminant ) / ( 2 * a );
        } else {
            // Get the stable roots for negative b
            outRoot1 = ( -b + sqrtDiscriminant ) / ( 2 * a );
            outRoot2 = ( 2 * c ) / ( -b + sqrtDiscriminant );
        }
    }

    return 2;
}

// Computes a root of the quadratic formula a*x^2 + b*x + c = 0, corresponding to (-b + sqrt(b^2-4ac))/2a.
// Returns true if it found the root.
template <class T>
inline bool get_quadratic_larger_root( T a, T b, T c, T& outRoot ) {
    if( a == 0 )
        throw std::runtime_error( "frantic::math::polynomial_roots::get_quadratic_roots: Invalid parameter: a = 0." );

    T discriminant = b * b - 4 * a * c;
    if( discriminant < 0 ) {
        outRoot = 0;
        return false; // No roots;
    }

    if( discriminant == 0 ) {
        // For stability reasons, choose the largest denominator
        if( fabs( 2 * a ) > fabs( b ) ) {
            outRoot = -b / ( 2 * a );
        } else {
            outRoot = ( 2 * c ) / ( -b );
        }
        return true;
    } else {
        T sqrtDiscriminant = sqrt( discriminant );

        // The reason we discriminate between positive and negative values of B is
        // because we want B to be the same sign as sqrtDiscriminant. This is so that
        // we do not have rounding errors when B is close to sqrtDiscriminant.

        if( b >= 0 ) {
            // Get the stable roots for positive b
            outRoot = ( 2 * c ) / ( -b - sqrtDiscriminant );
        } else {
            // Get the stable roots for negative b
            outRoot = ( -b + sqrtDiscriminant ) / ( 2 * a );
        }
    }

    return true;
}

// Cubic function x^3 + a*x^2 + b*x + c = 0
// Returns the number of roots.
inline int get_cubic_roots( double a, double b, double c, double& outRoot1, double& outRoot2, double& outRoot3 ) {
    double Q = ( a * a - 3 * b ) / 9, R = ( 2 * a * a * a - 9 * a * b + 27 * c ) / 54;

    double Qcubed = Q * Q * Q, Rsquared = R * R;
    ;

    if( Rsquared < Qcubed ) { // The equation has 3 real roots.
        double theta = acos( R / sqrt( Qcubed ) );
        double sqrtQ = sqrt( Q );
        outRoot1 = -2 * sqrtQ * cos( theta / 3 ) - a / 3;
        outRoot2 = -2 * sqrtQ * cos( ( theta + 2 * M_PI ) / 3 ) - a / 3;
        outRoot3 = -2 * sqrtQ * cos( ( theta - 2 * M_PI ) / 3 ) - a / 3;
        return 3;
    } else {
        double A = ( R > 0 ? -1 : ( R < 0 ? 1 : 0 ) ) * pow( fabs( R ) + sqrt( Rsquared - Qcubed ), 1 / 3.0 );
        double B = A == 0 ? 0 : Q / A;
        // std::cout << "Q: " << Q << ", R: " << R << ", A: " << A << ", B: " << B << std::endl;
        outRoot1 = A + B - a / 3;
        outRoot2 = 0;
        outRoot3 = 0;
        return 1;
    }
}
} // namespace polynomial_roots
} // namespace math
} // namespace frantic
