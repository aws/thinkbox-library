// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>

namespace frantic {
namespace channels {

namespace detail {
struct cml_lerp_block {
    // The position and element count to lerp, in bytes.
    unsigned position, count;
    // The function which knows how to convert such a block.
    channel_weighted_sum_combine_function_t lerpCombine;

    cml_lerp_block( unsigned position_, unsigned count_, channel_weighted_sum_combine_function_t lerpCombine_ )
        : position( position_ )
        , count( count_ )
        , lerpCombine( lerpCombine_ ) {}
};
} // namespace detail

/**
 * The channel_map_lerp class can be used to do linear interpolations between two particles of the same type.
 *
 * For example:
 *
 *  channel_map cm;
 *  ...
 *  channel_map_lerp cml(cm);
 *  cml.lerp( destParticle, particleA, particleB, 0.25f );
 */
class channel_map_lerp {
    // These are all the runs of values that can be lerp'd
    std::vector<detail::cml_lerp_block> m_lerpBlocks;
    std::size_t m_structureSize;

  public:
    /**
     * Constructs an empty channel_map_lerp.
     */
    channel_map_lerp() { m_structureSize = 0; }

    /**
     * Constructs a channel_map_lerp which does linear interpolations between structures of the provided channel_map.
     */
    channel_map_lerp( const channel_map& channelMap ) { set( channelMap ); }

    /**
     * Clears the channel_map_lerp so that it operates on a zero-sized structure.
     */
    void clear() {
        m_lerpBlocks.clear();
        m_structureSize = 0;
    }

    /**
     * This is the size of the structure that this channel_map_lerp operates on.
     */
    std::size_t structre_size() const { return m_structureSize; }

    /**
     * This function resets the existing state of the channel_map_lerp instance, and configures it to
     * do linear interpolations between particles of the provided structure.
     */
    void set( const channel_map& channelMap );

    /**
     * This does a linear interpolation between the two source particles, into the destination particle.
     */
    void lerp( void* destData, const void* sourceDataA, const void* sourceDataB, float t );

    /**
     * This does a linear interpolation of the destination particle with the source particle.
     */
    void lerp( void* destDataA, const void* sourceDataB, float t );

    /**
     * This does a trilinear interpoation of the 8 source particles, into the dest particle.
     */
    void trilerp( void* dest, const void* sourceData[8], float ( &offsets )[3] );
};

} // namespace channels
} // namespace frantic
