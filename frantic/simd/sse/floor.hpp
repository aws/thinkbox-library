#pragma once

#include <mmintrin.h>

namespace frantic {
namespace simd {
namespace sse {

/**
 *  Round down to the nearest integer.
 *
 *  Implementation based on:
 *
 *  Stephanie Rancourt, "Pre SSE 4.1 floor/ceil/round functions (+modulo bonus)"
 *  http://dss.stephanierct.com/DevBlog/?p=8
 *  Licensed under Boost Software License
 */
inline __m128 floor( const __m128& a ) {
    // Keep the original value if it's too large to round
    __m128 aAbs = _mm_and_ps( a, _mm_castsi128_ps( _mm_set1_epi32( 0x7fffffff ) ) );
    __m128 largeMask = _mm_set1_ps( 8388608.f );
    __m128 mask = _mm_cmple_ps( aAbs, largeMask );

    __m128 aRoundTrip = _mm_cvtepi32_ps( _mm_cvttps_epi32( a ) );
    __m128 roundTripIncreased = _mm_cmpgt_ps( aRoundTrip, a );
    __m128 correction = _mm_cvtepi32_ps( _mm_castps_si128( roundTripIncreased ) );
    __m128 aFloor = _mm_add_ps( aRoundTrip, correction );

    return _mm_or_ps( _mm_and_ps( mask, aFloor ), _mm_andnot_ps( mask, a ) );
}

} // namespace sse
} // namespace simd
} // namespace frantic
