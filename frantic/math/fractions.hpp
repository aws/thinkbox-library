// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

namespace frantic {
namespace math {

/**
 * This function computes a rational number that is guaranteed to be a better approximation for the given
 * floating point argument than any other rational number with a smaller denominator.
 * It uses continued fractions to iteratively generate the rational number and stops once the numerator or
 * denominator exceeds the supplied value (minNumDenom) or 100 continued fraction terms, whichever occurs first.
 * The result is a rational number, returned as a numerator/denominator pair.
 */
inline std::pair<boost::int64_t, boost::int64_t>
get_rational_representation( double d, boost::int64_t minNumDenom = ( 1 << 24 ) ) {
    boost::int64_t pPrev = 0, pCur = 1, qPrev = 1, qCur = 0;
    double val = std::abs( d );
    for( int i = 0; i < 100; ++i ) {
        boost::int64_t a = static_cast<boost::int64_t>( val );

        // Make sure to handle overflows gracefully. For small input values, a can become very large so (a * qCur)
        // overflows.
        // TODO: 'a' can only be 0 on the first iteration. If this loop is unrolled, this check can be eliminated.
        // HACK: This does not prevent an integer overflow from the addition in the recurrence relation!!!!
        if( a != 0 ) {
            boost::int64_t maxVal = std::numeric_limits<boost::int64_t>::max() / a;
            if( qCur > maxVal || pCur > maxVal )
                break;
        }

        boost::int64_t pNext = a * pCur + pPrev;
        boost::int64_t qNext = a * qCur + qPrev;

        pPrev = pCur;
        pCur = pNext;
        qPrev = qCur;
        qCur = qNext;

        // Break if the desired denominator accuracy has been reached.
        // Also, break if the numerator is large to prevent overflow.
        if( pCur > minNumDenom || qCur > minNumDenom )
            break;

        double trunc = ( val - std::floor( val ) );
        if( trunc <= 0.0 ) // Should this use an epsilon?
            break;

        val = 1.0 / trunc;
    }

    if( d < 0 )
        return std::make_pair( -pCur, qCur );
    return std::make_pair( pCur, qCur );
}

} // namespace math
} // namespace frantic
