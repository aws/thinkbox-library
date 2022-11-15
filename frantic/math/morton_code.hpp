// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/vector2f.hpp>

namespace frantic {
namespace math {

using frantic::graphics::vector3f;
using frantic::graphics2d::vector2f;

namespace detail {

/**
 * Separates the bits of x evenly with two empty bits inbetween. Used for 3D version of Morton Code.
 * Only the first 21 bits of x are used, anything bigger is truncated.
 */
inline boost::uint64_t part_by_2( boost::uint64_t x ) {
    x &= UINT64_C( 0x00000000001fffff ); // x = -------- -------- -------- -------- -------- ---xxxxx xxxxxxxx xxxxxxxx
    x = ( x ^ ( x << 32 ) ) &
        UINT64_C( 0x00ff00000000ffff ); // x = -------- ---xxxxx -------- -------- -------- -------- xxxxxxxx xxxxxxxx
    x = ( x ^ ( x << 16 ) ) &
        UINT64_C( 0x00ff0000ff0000ff ); // x = -------- ---xxxxx -------- -------- xxxxxxxx -------- -------- xxxxxxxx
    x = ( x ^ ( x << 8 ) ) &
        UINT64_C( 0x100f00f00f00f00f ); // x = ---x---- ----xxxx -------- xxxx---- ----xxxx -------- xxxx---- ----xxxx
    x = ( x ^ ( x << 4 ) ) &
        UINT64_C(
            0x10c30c30c30c30c3 ); // x = ---x---- xx----xx ----xx-- --xx---- xx----xx	----xx-- --xx---- xx----xx
    x = ( x ^ ( x << 2 ) ) &
        UINT64_C( 0x1249249249249249 ); // x = ---x--x- -x--x--x --x--x-- x--x--x- -x--x--x --x--x-- x--x--x- -x--x--x
    return x;
}

/**
 * Separates the bits of x evenly with one empty bit inbetween. Used for 2D version of Morton Code.
 * Only the first 32 bits of x are used, anything bigger is truncated.
 */
inline boost::uint64_t part_by_1( boost::uint64_t x ) {
    x &= 0x00000000ffffffffLL; // x = -------- -------- -------- -------- xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
    x = ( x ^ ( x << 16 ) ) &
        UINT64_C( 0x0000ffff0000ffff ); // x = -------- -------- xxxxxxxx xxxxxxxx -------- -------- xxxxxxxx xxxxxxxx
    x = ( x ^ ( x << 8 ) ) &
        UINT64_C( 0x00ff00ff00ff00ff ); // x = -------- xxxxxxxx -------- xxxxxxxx -------- xxxxxxxx -------- xxxxxxxx
    x = ( x ^ ( x << 4 ) ) &
        UINT64_C( 0x0f0f0f0f0f0f0f0f ); // x = ----xxxx ----xxxx ----xxxx ----xxxx ----xxxx ----xxxx ----xxxx ----xxxx
    x = ( x ^ ( x << 2 ) ) &
        UINT64_C( 0x3333333333333333 ); // x = --xx--xx --xx--xx --xx--xx --xx--xx --xx--xx --xx--xx --xx--xx --xx--xx
    x = ( x ^ ( x << 1 ) ) &
        UINT64_C( 0x5555555555555555 ); // x = -x-x-x-x -x-x-x-x -x-x-x-x -x-x-x-x -x-x-x-x -x-x-x-x -x-x-x-x -x-x-x-x
    return x;
}

} // namespace detail

/**
 * Returns the morton code for position in relation to the boundingBoxCorner (origin). Morton code interleaves the bits
 * from the y and x values of position.
 */
inline boost::uint64_t morton_code( const vector2f& boundingBoxCorner, float voxelLength, const vector2f& position ) {
    vector2f v = ( position - boundingBoxCorner ) / voxelLength;
    return ( detail::part_by_1( static_cast<boost::uint64_t>( v.y ) ) << 1 ) +
           detail::part_by_1( static_cast<boost::uint64_t>( v.x ) );
}

/**
 * Returns the morton code for position in relation to the boundingBoxCorner (origin). Morton code interleaves the bits
 * from the z, y and x values of position.
 */
inline boost::uint64_t morton_code( const vector3f& boundingBoxCorner, float voxelLength, const vector3f& position ) {
    vector3f v = ( position - boundingBoxCorner ) / voxelLength;
    return ( detail::part_by_2( static_cast<boost::uint64_t>( v.z ) ) << 2 ) +
           ( detail::part_by_2( static_cast<boost::uint64_t>( v.y ) ) << 1 ) +
           detail::part_by_2( static_cast<boost::uint64_t>( v.x ) );
}

} // namespace math
} // namespace frantic
