// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file iterative_root_solve.hpp
 * @author Darcy Harrison
 *
 * This file contains routines for iteratively solving functions that are not known to be polynomials.
 * For root solving polynomials there are better functions, try polynomial_roots.hpp
 */
#pragma once

#include <limits>

namespace frantic {
namespace math {

/**
 * This function will solve for the root of a function of one variable. It uses a bisection method
 * that is guaranteed to converge. The user is responsible for supplying a lower and upper bound
 * on a region that contains a root. Essentially it must be true that: fn(lb) <= 0 <= fn(ub).
 *
 * @note This implementation will linearly solve for the result once the bisection range is smaller
 *        than the supplied tolerance. This means that functions can use a largish tolerance value
 *        if the expect to behave linearly in small regions.
 *
 * @tparam Fn This class must have operator()(float) that is the numerical function to find the roots for.
 * @param fn The function to solve the roots for.
 * @param lb The lower bound on the region containing the root.
 * @param ub The upper bound on the region containing the root.
 * @param tol The amount of absolute error to expect in the result.
 * @return The root of the function, within the given tolerance.
 */
template <class Fn>
inline float find_root_bisection( const Fn& fn, float lb, float ub, float tol = 1e-3f ) {
    float diff, mid, val;

    for( ;; ) {
        diff = ub - lb;
        if( diff < tol ) {
            // Within this small interval, assume the function is linear. Consider the case where the
            // true solution is on one endpoint of the bisection range. If the function were assumed
            // constant in this range, it would never correctly choose one of the endpoints.
            float lbVal = fn( lb );
            float ubVal = fn( ub );
            float valDiff = ( ubVal - lbVal );
            // If both endpoints have a similar value, choose the midpoint in order to avoid
            // the divide by almost zero. This should always be positive unless badness happened.
            if( valDiff < 1e-5f || lbVal > 0 || ubVal < 0 )
                return lb + 0.5f * diff;
            return lb - lbVal * ( diff / valDiff );
        }
        mid = lb + 0.5f * diff;
        val = fn( mid );
        if( val < 0 )
            lb = mid;
        else
            ub = mid;
    }
}

/**
 * This function will solve for the root of a function of one variable. It uses a newton-raphson iteration, but trys
 * to avoid stepping outside of the known bracket. It takes a bisection step in order to prevent this. The user is
 * responsible for supplying a lower and upper bound on a region that contains a root. Essentially it must
 * be true that: fn(lb) <= 0 <= fn(ub).
 *
 * @tparam Fn  This class must have operator()(float,float&) that is the numerical function to find the roots for. The
 *             derivative of the function is passed as an output reference.
 * @param fn The function to solve the roots for.
 * @param lb The lower bound on the region containing the root.
 * @param ub The upper bound on the region containing the root.
 * @param guess An initial estimate of the root. The better this is, the faster the function will converge.
 * @param tol The amount of absolute error to expect in the result.
 * @return The root of the function, within the given tolerance.
 */
template <class Fn>
inline float find_root_newton_raphson( const Fn& fn, float lb, float ub,
                                       float guess = ( std::numeric_limits<float>::max )(), float tol = 1e-3f ) {
    float x;
    float val, derivative;

    x = guess != ( std::numeric_limits<float>::max )() ? guess : 0.5f * ( lb + ub );
    do {
        val = fn( x, derivative );

        // Is the current value effectively zero?
        if( fabsf( val ) < 1e-5f )
            return x;
        if( val < 0 )
            lb = x;
        else // val > 0
            ub = x;

        // If the newton-raphson step will take us outside the known bounds, do a bisection step.
        // This should also handle the case of a zero derivative, since nextX will become infinite
        // and this line will produce false.
        float nextX = x - val / derivative;
        if( nextX > lb && nextX < ub ) {
            x = nextX;
        } else {
            x = lb + 0.5f * ( ub - lb );
        }
        // Has the bracketing range for the root converged?
    } while( ub - lb > tol );

    // TODO: Another option is to do a linear solve in [lb,ub] for the root. This is done in find_root_bisection().
    return x;
}
} // namespace math
} // namespace frantic
