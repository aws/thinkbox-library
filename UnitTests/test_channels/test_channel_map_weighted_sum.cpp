// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <gtest/gtest.h>

#include <frantic/channels/channel_map_weighted_sum.hpp>
#include <frantic/channels/named_channel_data.hpp>

using namespace frantic::channels;
using namespace frantic::graphics;

// This class extends the channel_map_weighted_sum class and is merely for test purposes
namespace {
class channel_map_weighted_sum_accessor : public channel_map_weighted_sum {
  public:
    channel_map_weighted_sum_accessor( const channel_map& channelMap, const channel_propagation_policy& cpp )
        : channel_map_weighted_sum( channelMap, cpp ) {}

    std::size_t size_of_weighted_sum_data() { return m_weightedSumData.size(); }

    std::vector<std::size_t> arity_list() {
        std::vector<size_t> arity;
        for( size_t i = 0; i < m_weightedSumData.size(); ++i ) {
            arity.push_back( m_weightedSumData[i].m_arity );
        }
        return arity;
    }

    std::vector<std::size_t> offset_list() {
        std::vector<size_t> offset;
        for( size_t i = 0; i < m_weightedSumData.size(); ++i ) {
            offset.push_back( m_weightedSumData[i].m_offset );
        }
        return offset;
    }
}; // channel_map_weighted_sum_accessor

// Particle for running through channelMapWeightedSum.
struct test_particle {
    vector3f pos;
    float r;
    vector3f vel;
    vector3 color;
    vector3 color2;

    test_particle( vector3f pos, float r, vector3f vel, vector3 color, vector3 color2 )
        : pos( pos )
        , r( r )
        , vel( vel )
        , color( color )
        , color2( color2 ) {}
};

} // Anonymous namespace

// Test channels that are contiguous in the channel_map, but not contiguous in memory.
TEST( ChannelMapWeightedSum, NonContiguousMemory ) {
    channel_map map;
    map.define_channel( _T("Position"), 3, data_type_float32, 0 );
    map.define_channel( _T("Radius"), 1, data_type_float32, 12 );
    map.define_channel( _T("Velocity"), 3, data_type_float32, 28 );
    map.define_channel( _T("Color"), 3, data_type_int32, 40 );
    map.end_channel_definition( 4, true, false );

    channel_propagation_policy cpp( false );

    channel_map_weighted_sum_accessor channelMapWeightedSum( map, cpp );

    std::vector<size_t> arity = channelMapWeightedSum.arity_list();
    EXPECT_EQ( 4, arity[0] );
    EXPECT_EQ( 3, arity[1] );
    EXPECT_EQ( 3, arity[2] );

    std::vector<size_t> offset = channelMapWeightedSum.offset_list();
    EXPECT_EQ( 0, offset[0] );
    EXPECT_EQ( 28, offset[1] );
    EXPECT_EQ( 40, offset[2] );

    EXPECT_EQ( 3, channelMapWeightedSum.size_of_weighted_sum_data() );
}

// Test using a channel_map that has adjacent channels adjacent in memory.
TEST( ChannelMapWeightedSum, ContinguousMemory ) {
    std::vector<char*> particles;
    test_particle particle( vector3f( 1, 0, 0 ), 1, vector3f( 1, 1, 1 ), vector3( 1, 1, 1 ), vector3( 3, 2, 4 ) );
    particles.push_back( reinterpret_cast<char*>( &particle ) );

    channel_map map;
    map.define_channel( _T("Position"), 3, data_type_float32, 0 );
    map.define_channel( _T("Radius"), 1, data_type_float32, 12 );
    map.define_channel( _T("Velocity"), 3, data_type_float32, 16 );
    map.define_channel( _T("Color"), 3, data_type_int32, 28 );
    map.define_channel( _T("Color2"), 3, data_type_int32, 40 );
    map.end_channel_definition( 4, true, false );
    channel_propagation_policy cpp( false );

    channel_map_weighted_sum_accessor channelMapWeightedSum( map, cpp );

    std::vector<size_t> arity = channelMapWeightedSum.arity_list();
    EXPECT_EQ( 7, arity[0] );
    EXPECT_EQ( 6, arity[1] );

    std::vector<size_t> offset = channelMapWeightedSum.offset_list();
    EXPECT_EQ( 0, offset[0] );
    EXPECT_EQ( 28, offset[1] );

    EXPECT_EQ( 2, channelMapWeightedSum.size_of_weighted_sum_data() );

    std::vector<float> weights;
    weights.push_back( 1.0098f );
    std::vector<char> data( map.structure_size() );
    EXPECT_NO_THROW( channelMapWeightedSum.channel_weighted_sum( weights, particles, data ) );
}

// Test using a policy that does not include all of the channels in the map.
TEST( ChannelMapWeightedSum, ExcludingChannels ) {
    std::vector<char*> particles;
    test_particle particle( vector3f( 1, 0, 0 ), 1, vector3f( 1, 1, 1 ), vector3( 1, 1, 1 ), vector3( 5, 5, 5 ) );
    particles.push_back( reinterpret_cast<char*>( &particle ) );

    channel_map map;
    map.define_channel( _T("Position"), 3, data_type_float32, 0 );
    map.define_channel( _T("Radius"), 1, data_type_float32, 12 );
    map.define_channel( _T("Velocity"), 3, data_type_float32, 16 );
    map.define_channel( _T("Color"), 3, data_type_int32, 28 );
    map.define_channel( _T("Color2"), 3, data_type_int32, 40 );
    map.end_channel_definition( 4, true, false );

    channel_propagation_policy cpp( true );
    cpp.add_channel( _T("Position") );
    cpp.add_channel( _T("Velocity") );
    cpp.add_channel( _T("Color2") );

    channel_map_weighted_sum_accessor channelMapWeightedSum( map, cpp );

    std::vector<size_t> arity = channelMapWeightedSum.arity_list();
    EXPECT_EQ( 3, arity[0] );
    EXPECT_EQ( 3, arity[1] );
    EXPECT_EQ( 3, arity[2] );

    std::vector<size_t> offset = channelMapWeightedSum.offset_list();
    EXPECT_EQ( 0, offset[0] );
    EXPECT_EQ( 16, offset[1] );
    EXPECT_EQ( 40, offset[2] );

    EXPECT_EQ( 3, channelMapWeightedSum.size_of_weighted_sum_data() );

    std::vector<float> weights;
    weights.push_back( 1.0098f );
    std::vector<char> data( map.structure_size() );
    EXPECT_NO_THROW( channelMapWeightedSum.channel_weighted_sum( weights, particles, data ) );
}

TEST( ChannelMapWeightedSum, RawDataPointers ) {
    channel_map map;
    map.define_channel( _T("Position"), 3, data_type_float32, 0 );
    map.define_channel( _T("Radius"), 1, data_type_float32, 12 );
    map.define_channel( _T("Velocity"), 3, data_type_float32, 16 );
    map.define_channel( _T("Color"), 3, data_type_int32, 28 );
    map.define_channel( _T("Color2"), 3, data_type_int32, 40 );
    map.end_channel_definition( 4, true, false );
    channel_propagation_policy cpp( false );

    channel_map_weighted_sum_accessor channelMapWeightedSum( map, cpp );

    std::vector<size_t> arity = channelMapWeightedSum.arity_list();
    EXPECT_EQ( 7, arity[0] );
    EXPECT_EQ( 6, arity[1] );

    std::vector<size_t> offset = channelMapWeightedSum.offset_list();
    EXPECT_EQ( 0, offset[0] );
    EXPECT_EQ( 28, offset[1] );

    EXPECT_EQ( 2, channelMapWeightedSum.size_of_weighted_sum_data() );

    std::vector<float> weights;
    weights.push_back( 1.0098f );
    std::vector<char> data( map.structure_size() );
    std::vector<char*> particles;
    test_particle particle( vector3f( 1, 0, 0 ), 1, vector3f( 1, 1, 1 ), vector3( 1, 1, 1 ), vector3( 3, 2, 4 ) );
    particles.push_back( reinterpret_cast<char*>( &particle ) );
    EXPECT_NO_THROW( channelMapWeightedSum.channel_weighted_sum( &weights[0], &particles[0], 1, &data[0] ) );
}

TEST( ChannelMapWeightedSum, DataOutput ) {
    struct small_particle {
        vector3f pos;
        float radius;
        vector3 color;

        small_particle( vector3f pos, float r, vector3 color )
            : pos( pos )
            , radius( r )
            , color( color ) {}
    }; // small_particle

    // Check to ensure that all of the channels have weighted sum data.
    channel_map map;
    map.define_channel( _T("Position"), 3, data_type_float32, 0 );
    map.define_channel( _T("Radius"), 1, data_type_float32, 12 );
    map.define_channel( _T("Color"), 3, data_type_int32, 16 );
    map.end_channel_definition( 4, true, false );

    channel_propagation_policy cpp;
    channel_map_weighted_sum weightedSumObject( map, cpp );

    std::vector<char*> particles;
    small_particle particle1( vector3f( 200, 0, 50 ), 4, vector3( 4, 1, 7 ) );
    small_particle particle2( vector3f( 10, 50, 100 ), 2, vector3( 9, 15, 55 ) );
    particles.push_back( reinterpret_cast<char*>( &particle1 ) );
    particles.push_back( reinterpret_cast<char*>( &particle2 ) );

    std::vector<float> weights;
    weights.push_back( 0.5 );
    weights.push_back( 1 );

    // Fill in expectedData with what should be in the weighted sum
    std::vector<char> data( weightedSumObject.structure_size() );
    std::vector<char> expectedData;
    expectedData.resize( map.structure_size() );

    // Get the weighted sum of the floats.
    weighted_sum_float_combine_general<float>( &weights[0], &particles[0], 2, 4, &expectedData[0] );
    // Weighted sum of integer fields will fill with the highest weighted value.
    expectedData[16] = 9;
    expectedData[20] = 15;
    expectedData[24] = 55;

    weightedSumObject.channel_weighted_sum( weights, particles, data );

    for( size_t i = 0; i < data.size(); ++i ) {
        EXPECT_EQ( expectedData[i], data[i] );
    }

    // Check to ensure that channels that are skipped do not have weighted sum data.
    cpp.add_channel( _T( "Radius" ) );

    expectedData.clear();
    expectedData.resize( map.structure_size() );
    data.clear();
    data.resize( map.structure_size() );
    // refill expectedData.
    // Get the weighted sum of the floats.
    weighted_sum_float_combine_general<float>( &weights[0], &particles[0], 2, 3, &expectedData[0] );
    // Weighted sum of integer fields will fill with the highest weighted value.
    expectedData[16] = 9;
    expectedData[20] = 15;
    expectedData[24] = 55;

    weightedSumObject = channel_map_weighted_sum( map, cpp );
    weightedSumObject.channel_weighted_sum( weights, particles, data );
    for( size_t i = 0; i < data.size(); ++i ) {
        EXPECT_EQ( expectedData[i], data[i] );
    }
}
