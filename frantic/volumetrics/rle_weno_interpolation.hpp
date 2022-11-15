// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/fluids/velocity_advection.hpp>
#include <frantic/volumetrics/levelset/rle_defined_box_iterator.hpp>

namespace frantic {
namespace volumetrics {

inline float weno_indicator1( float f1, float f2, float f3 ) {
    return ( 26.f * f3 * f1 - 52.f * f2 * f1 - 76.f * f3 * f2 + 25.f * f3 * f3 + 64.f * f2 * f2 + 13.f * f1 * f1 ) /
           12.f;
}

inline float weno_indicator2( float f2, float f3, float f4 ) {
    return ( 26.f * f4 * f2 - 52.f * f4 * f3 - 76.f * f3 * f2 + 25.f * f2 * f2 + 64.f * f3 * f3 + 13.f * f4 * f4 ) /
           12.f;
}

inline float weno_polynomial1( float f1, float f2, float f3, float w ) {
    // FF_LOG(debug)  << "\tf1:" << f[0] << "\tf2:" << f[1] << "\tf3:" << f[2]
    return f2 + ( f3 - f1 ) * w * 0.5f + ( f3 - 2 * f2 + f1 ) * w * w * 0.5f;
}

inline float weno_polynomial2( float f2, float f3, float f4, float w ) {
    return f2 + ( -f4 + 4 * f3 - 3 * f2 ) * w * 0.5f + ( f4 - 2 * f3 + f2 ) * w * w * 0.5f;
}

/**
 * Computes a WENO3 interpolation of the 4 data points provided given that w = (x -x_i) and s1, s2 are precomputed
 * Smoothness Indicators
 *
 * @param f - array of 4 values f_i-1, f_i, f_i+1, fi+2
 * @param w - w = (x-x_i), where x is the interpolation point that lies between f_i and f_i+1
 * @param s1 - left hand smoothness indicator
 * @param s2 - right hand smoothness indicator
 * @returns the interpolated value
 */
inline float weno3_interpolant( const float* f, float w, float s1, float s2 ) {
    float p1 = weno_polynomial1( f[0], f[1], f[2], w );
    float p2 = weno_polynomial2( f[1], f[2], f[3], w );

    // FF_LOG(debug)  << "\tw:" << w << std::endl;
    // FF_LOG(debug)  << "\tf1:" << f[0] << "\tf2:" << f[1] << "\tf3:" << f[2] << "\tf4:" << f[3] << std::endl;
    // FF_LOG(debug)  << "\tp1:" << p1 << "\tp2:" << p2 << std::endl;
    // FF_LOG(debug)  << "\ts1:" << s1 << "\ts2:" << s2 << std::endl;

    float c1 = ( 2.f - w ) / 3.f;
    float c2 = ( 1.f + w ) / 3.f;
    // FF_LOG(debug) << "\tc1:" << c1 << "\tc2:" << c2 << std::endl;

    float a1 = static_cast<float>( c1 / ( ( 1e-6 + s1 ) * ( 1e-6 + s1 ) ) );
    float a2 = static_cast<float>( c2 / ( ( 1e-6 + s2 ) * ( 1e-6 + s2 ) ) );

    float asum = a1 + a2;
    float w1 = a1 / asum;
    float w2 = a2 / asum;
    // FF_LOG(debug)  << "\ta1:" << a1 << "\ta2:" << a2 << std::endl;
    // FF_LOG(debug)  << "\tw1:" << w1 << "\tw2:" << w2 << std::endl;
    // FF_LOG(debug)  << "\tr:" << w1*p1 + w2*p2 << std::endl;
    return w1 * p1 + w2 * p2;
}

/**
 * Computes a WENO3 interpolation of the 4 data points provided given that w = (x -x_i)
 * The Smoothness Indicators for the interpolationg are calculated directly from the
 * provided values.
 *
 * @param f - array of 4 values f_i-1, f_i, f_i+1, fi+2
 * @param w - w = (x-x_i), where x is the interpolation point that lies between f_i and f_i+1
 * @returns the interpolated value
 */
inline float weno3_interpolant( const float* f, float w ) {

    float s1 = weno_indicator1( f[0], f[1], f[2] );
    float p1 = weno_polynomial1( f[0], f[1], f[2], w );
    float s2 = weno_indicator2( f[1], f[2], f[3] );
    float p2 = weno_polynomial2( f[1], f[2], f[3], w );

    // FF_LOG(debug)  << "\tw:" << w << std::endl;
    // FF_LOG(debug)  << "\tf1:" << f[0] << "\tf2:" << f[1] << "\tf3:" << f[2] << "\tf4:" << f[3] << std::endl;
    // FF_LOG(debug)  << "\tp1:" << p1 << "\tp2:" << p2 << std::endl;
    // FF_LOG(debug)  << "\ts1:" << s1 << "\ts2:" << s2 << std::endl;

    float c1 = ( 2.f - w ) / 3.f;
    float c2 = ( 1.f + w ) / 3.f;
    // FF_LOG(debug) << "\tc1:" << c1 << "\tc2:" << c2 << std::endl;

    float a1 = static_cast<float>( c1 / ( ( 1e-6 + s1 ) * ( 1e-6 + s1 ) ) );
    float a2 = static_cast<float>( c2 / ( ( 1e-6 + s2 ) * ( 1e-6 + s2 ) ) );

    float asum = a1 + a2;
    float w1 = a1 / asum;
    float w2 = a2 / asum;
    // FF_LOG(debug)  << "\ta1:" << a1 << "\ta2:" << a2 << std::endl;
    // FF_LOG(debug)  << "\tw1:" << w1 << "\tw2:" << w2 << std::endl;
    // FF_LOG(debug)  << "\tr:" << w1*p1 + w2*p2 << std::endl;
    return w1 * p1 + w2 * p2;
}

float weno3_signed_distance_lookup( const boost::int32_t* const dataIndices, const frantic::graphics::size3& boxSize,
                                    const frantic::graphics::vector3& currentXYZMin,
                                    const std::vector<float>& distanceData, frantic::graphics::vector3f voxelLookup );

void staggered_weno3_debug_dump(
    const boost::int32_t* const dataIndices, const frantic::graphics::size3& boxSize,
    const frantic::graphics::vector3& currentXYZMin,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& indicatorAccessor1,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& indicatorAccessor2,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& velocityAccessor,
    frantic::graphics::vector3f voxelLookup );

frantic::graphics::vector3f staggered_weno3_lookup(
    const boost::int32_t* const dataIndices, const frantic::graphics::size3& boxSize,
    const frantic::graphics::vector3& currentXYZMin,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& indicatorXAccessor1,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& indicatorXAccessor2,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& velocityAccessor,
    frantic::graphics::vector3f voxelLookup );

void create_staggered_smoothness_indicator_x_channel( frantic::fluids::rle_voxel_field& field,
                                                      const frantic::tstring& StaggeredChannelName,
                                                      const frantic::tstring& indicatorChannelOneName,
                                                      const frantic::tstring& indicatorChannelTwoName );

} // namespace volumetrics
} // namespace frantic
