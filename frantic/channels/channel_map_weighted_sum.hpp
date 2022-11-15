// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_propagation_policy.hpp>

namespace frantic {
namespace channels {

/**
 * Class for grouping channels that have the same data type before calculating the
 * weighted sum.
 */
class channel_map_weighted_sum {
  public:
    channel_map_weighted_sum();
    channel_map_weighted_sum( const channel_map& channelMap, const channel_propagation_policy& policy );

    /**
     * Get the weighted sum from the channel_map by calling the weighted sum functions on groups of channels.
     *
     * @param weights    An array of floating point weights.
     * @param outData    A pointer to a data element for output.
     * @param particles  An array of pointers, each element corresponding to a weight.
     */
    void channel_weighted_sum( const std::vector<float>& weights, const std::vector<char*>& particles,
                               std::vector<char>& outData ) const;

    /**
     * Overload for working with raw data points.
     * Get the weighted sum from the channel_map by calling the weighted sum functions on groups of channels.
     *
     * @param weights    An array of floating point weights.
     * @param outData    A pointer to a data element for output.
     * @param particles  An array of pointers, each element corresponding to a weight.
     * @param count      The number of elements in each array.
     */
    void channel_weighted_sum( const float* weights, const char* const* particles, const size_t count, char* outData );

    /**
     * Get the size of the structure that the channel_map_weighted_sum takes up.
     */
    size_t structure_size() const;

  private:
    /**
     * Contains the data required, from the channel_map and channel_propagation_policy that are passed into
     * the constructor, for performing the weighted sum calculation on groups of channels.
     */
    struct weighted_sum_data {
        weighted_sum_data() {}
        weighted_sum_data( offset_input_channel_weighted_sum_combine_function_t weightedSum, std::size_t offset,
                           std::size_t arity )
            : m_weightedSumFunc( weightedSum )
            , m_offset( offset )
            , m_arity( arity ) {}

        offset_input_channel_weighted_sum_combine_function_t m_weightedSumFunc;
        size_t m_offset;
        size_t m_arity;
    }; // weighted_sum_data

    /**
     * Determine which offset_input_channel_weighted_sum_combine_functions will be needed and
     * group together the channels with the same data types in lists.
     *
     * Grouping the channels together before getting the weighted sum is a performance optimization.
     *
     * The functions will have the same index as the corresponding group of channels will.
     *
     * @param channelMap  The input channel to use when creating the lists.
     */
    void set_functions_and_adjacent_channels( const channel_map& channelMap, const channel_propagation_policy& policy );

    /**
     * Determine if a given channel index on the channel_map is included according to the
     * channel_propagation_policy. Does not matter if policy is an include or exclude policy.
     *
     * @param index       The index into the channel_map to try to find in the policy.
     * @param channelMap  The channel_map to get the channel at the given index from.
     */
    bool is_channel_included( std::size_t index, const channel_map& channelMap,
                              const channel_propagation_policy& policy ) const;

    /**
     * Determine if a channel_map's layout reflects the memory locations of channels.
     * ie: Consecutive channels in the map are consecutive in memory.
     *
     * @param channelMap  The channel_map to check.
     */
    bool channels_are_contiguous( const channel_map& channelMap, size_t channelOne, size_t channelTwo ) const;

  protected:
    std::vector<weighted_sum_data> m_weightedSumData;
    std::size_t m_structureSize;
};

} // namespace channels
} // namespace frantic
