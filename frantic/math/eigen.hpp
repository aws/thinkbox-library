// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics2d/vector2f.hpp>
#include <frantic/math/polynomial_roots.hpp>

namespace frantic {
namespace math {
namespace linearalgebra {

using frantic::graphics2d::vector2f;
using namespace frantic::math::polynomial_roots;

// Gets the eigenvalues of the 2x2 matrix
// [ a b ]
// [ c d ]
// and returns the number of eigenvalues
inline int get_eigenvalues( float a, float b, float c, float d, float& outEigen1, float& outEigen2 ) {
    // det(A - lambda I) = 0

    // det | a - lambda         b     |
    //     |     c         d - lambda |

    // (a-lambda)(d-lambda) - bc = 0

    return get_quadratic_roots( 1.0f, -( a + d ), ( a * d - b * c ), outEigen1, outEigen2 );
}

// helper functions for get_eigenvalues_symmetric_3x3
namespace {
int get_i_pp( int axis ) {
    int lut[] = { 3, 5, 0 };
    return lut[axis];
}
int get_i_qq( int axis ) {
    int lut[] = { 5, 0, 3 };
    return lut[axis];
}
int get_i_pq( int axis ) {
    int lut[] = { 4, 2, 1 };
    return lut[axis];
}
int get_i_rp( int axis ) {
    int lut[] = { 1, 4, 2 };
    return lut[axis];
}
int get_i_rq( int axis ) {
    int lut[] = { 2, 1, 4 };
    return lut[axis];
}
int get_i_p( int axis ) {
    int lut[] = { 1, 2, 0 };
    return lut[axis];
}
int get_i_q( int axis ) {
    int lut[] = { 2, 0, 1 };
    return lut[axis];
}

template <class T>
int get_max_axis( T a[6] ) {
    const T x = fabs( a[4] );
    const T y = fabs( a[2] );
    const T z = fabs( a[1] );

    if( x > y ) {
        if( x > z ) {
            return 0;
        } else {
            return 2;
        }
    } else {
        if( y > z ) {
            return 1;
        } else {
            return 2;
        }
    }
}

template <class T>
T get_t( T a[6], int axis ) {
    const T a_qq = a[get_i_qq( axis )];
    const T a_pp = a[get_i_pp( axis )];
    const T a_pq = a[get_i_pq( axis )];

    const T diff = a_qq - a_pp;

    const T abs_pq = fabs( a_pq );
    const T abs_diff = fabs( diff );

    if( abs_diff + abs_pq == abs_diff ) {
        return a_pq / diff;
    } else {
        const T theta = diff / ( 2 * a_pq );
        const T t = 1 / ( fabs( theta ) + sqrt( theta * theta + 1 ) );
        return ( theta < 0 ) ? -t : t;
    }
}
} // namespace

/**
 *  Get eigenvalues of a symmetric 3x3 matrix.
 *
 *	This function uses the Jacobi method.  See:
 *
 *		Press et al., "Numerical Recipes in C++".  2nd ed.
 *
 * @param m the coefficients of a symmetric 3x3 matrix, stored in the order:
 *			[ 0 1 2 ]
 *			[   3 4 ]
 *			[     5 ]
 * @param[out] outEigen1 an eigenvalue of the matrix.
 * @param[out] outEigen2 an eigenvalue of the matrix.
 * @param[out] outEigen3 an eigenvalue of the matrix.
 */
template <class T>
void get_eigenvalues_symmetric_3x3( const T m[6], T& outEigen1, T& outEigen2, T& outEigen3 ) {
    const int iterMax = 50;
    const T thresh = 0;
    T a[6];
    memcpy( a, m, 6 * sizeof( T ) );

    for( int iter = 0; iter < iterMax; ++iter ) {
        const int axis = get_max_axis( a );

        const T a_pq = a[get_i_pq( axis )];

        if( fabs( a_pq ) <= thresh ) {
            break;
        }

        const T t = get_t( a, axis );
        const T c = 1 / sqrt( t * t + 1 );
        const T s = t * c;
        const T tau = s / ( 1 + c );
        const T a_rp = a[get_i_rp( axis )];
        const T a_rq = a[get_i_rq( axis )];

        a[get_i_pq( axis )] = 0;
        a[get_i_pp( axis )] -= t * a_pq;
        a[get_i_qq( axis )] += t * a_pq;
        a[get_i_rp( axis )] -= s * ( a_rq + tau * a_rp );
        a[get_i_rq( axis )] += s * ( a_rp - tau * a_rq );
    }

    outEigen1 = a[0];
    outEigen2 = a[3];
    outEigen3 = a[5];
}

/**
 *  Get eigenvectors and eigenvalues of a symmetric 3x3 matrix.
 *
 *	This function uses the Jacobi method.  See:
 *
 *		Press et al., "Numerical Recipes in C++".  2nd ed.
 *
 * @param m the coefficients of a symmetric 3x3 matrix, stored in the order:
 *			[ 0 1 2 ]
 *			[   3 4 ]
 *			[     5 ]
 * @param[out] outVector1 an eigenvector of the matrix.
 * @param[out] outVector2 an eigenvector of the matrix.
 * @param[out] outVector3 an eigenvector of the matrix.
 * @param[out] outValues eigenvalues of the matrix.
 */
template <class T>
void eigendecompose_symmetric_3x3( const T m[6], T outVector1[3], T outVector2[3], T outVector3[3], T outValues[3] ) {
    const int iterMax = 50;
    const T thresh = 0;
    T a[6];
    memcpy( a, m, 6 * sizeof( T ) );

    T eigenvec[3][3];
    memset( eigenvec, 0, 9 * sizeof( T ) );
    for( int i = 0; i < 3; ++i ) {
        eigenvec[i][i] = T( 1 );
    }

    for( int iter = 0; iter < iterMax; ++iter ) {
        const int axis = get_max_axis( a );

        const T a_pq = a[get_i_pq( axis )];

        if( fabs( a_pq ) <= thresh ) {
            break;
        }

        const T t = get_t( a, axis );
        const T c = 1 / sqrt( t * t + 1 );
        const T s = t * c;
        const T tau = s / ( 1 + c );
        const T a_rp = a[get_i_rp( axis )];
        const T a_rq = a[get_i_rq( axis )];

        a[get_i_pq( axis )] = 0;
        a[get_i_pp( axis )] -= t * a_pq;
        a[get_i_qq( axis )] += t * a_pq;
        a[get_i_rp( axis )] -= s * ( a_rq + tau * a_rp );
        a[get_i_rq( axis )] += s * ( a_rp - tau * a_rq );

        const int p = get_i_p( axis );
        const int q = get_i_q( axis );
        for( int i = 0; i < 3; ++i ) {
            const T e_ip = eigenvec[i][p];
            const T e_iq = eigenvec[i][q];
            eigenvec[i][p] -= s * ( e_iq + tau * e_ip );
            eigenvec[i][q] += s * ( e_ip - tau * e_iq );
        }
    }

    for( int i = 0; i < 3; ++i ) {
        outVector1[i] = eigenvec[i][0];
        outVector2[i] = eigenvec[i][1];
        outVector3[i] = eigenvec[i][2];
    }

    outValues[0] = a[0];
    outValues[1] = a[3];
    outValues[2] = a[5];
}

} // namespace linearalgebra
} // namespace math
} // namespace frantic
