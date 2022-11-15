// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/channel_map_weighted_sum.hpp>

using namespace frantic::channels;

channel_map_weighted_sum::channel_map_weighted_sum()
    : m_structureSize( 0 ) {}

channel_map_weighted_sum::channel_map_weighted_sum( const channel_map& channelMap,
                                                    const channel_propagation_policy& policy ) {
    set_functions_and_adjacent_channels( channelMap, policy );
    m_structureSize = channelMap.structure_size();
}

bool channel_map_weighted_sum::is_channel_included( std::size_t index, const channel_map& channelMap,
                                                    const channel_propagation_policy& policy ) const {
    return policy.is_channel_included( channelMap[index].name() );
}

void channel_map_weighted_sum::set_functions_and_adjacent_channels( const channel_map& channelMap,
                                                                    const channel_propagation_policy& policy ) {
    for( size_t i = 0; i < channelMap.channel_count(); ) {
        if( is_channel_included( i, channelMap, policy ) ) {
            size_t offset = channelMap[i].offset();
            size_t arity = channelMap[i].arity();

            size_t currIndex = i + 1;
            while( currIndex < channelMap.channel_count() &&
                   channelMap[i].data_type() == channelMap[currIndex].data_type() &&
                   is_channel_included( currIndex, channelMap, policy ) &&
                   channels_are_contiguous( channelMap, currIndex - 1, currIndex ) ) {
                arity += channelMap[currIndex].arity();
                ++currIndex;
            }

            m_weightedSumData.push_back( weighted_sum_data(
                offset_input_channel_weighted_sum_combine_function( channelMap[i].data_type() ), offset, arity ) );
            i = currIndex;
        } else {
            ++i;
        }
    }
}

void channel_map_weighted_sum::channel_weighted_sum( const std::vector<float>& weights,
                                                     const std::vector<char*>& particles,
                                                     std::vector<char>& outData ) const {
    for( size_t i = 0; i < m_weightedSumData.size(); ++i ) {
        weighted_sum_data weightedSumData = m_weightedSumData[i];
        weightedSumData.m_weightedSumFunc( &weights[0], weightedSumData.m_offset, &particles[0], particles.size(),
                                           weightedSumData.m_arity, &outData[weightedSumData.m_offset] );
    }
}

void channel_map_weighted_sum::channel_weighted_sum( const float* weights, const char* const* particles,
                                                     const size_t count, char* outData ) {
    for( size_t i = 0; i < m_weightedSumData.size(); ++i ) {
        weighted_sum_data weightedSumData = m_weightedSumData[i];
        weightedSumData.m_weightedSumFunc( weights, weightedSumData.m_offset, particles, count, weightedSumData.m_arity,
                                           outData + weightedSumData.m_offset );
    }
}

bool channel_map_weighted_sum::channels_are_contiguous( const channel_map& channelMap, size_t channelOne,
                                                        size_t channelTwo ) const {
    if( channelMap[channelTwo].offset() ==
        ( channelMap[channelOne].offset() + channelMap[channelOne].primitive_size() ) ) {
        return true;
    } else {
        return false;
    }
}

size_t channel_map_weighted_sum::structure_size() const { return m_structureSize; }
