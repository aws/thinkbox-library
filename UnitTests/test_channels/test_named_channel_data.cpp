// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/named_channel_data.hpp>

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif

using namespace frantic::channels;

template <typename T>
class OffsetInputWeightedSumCombine : public ::testing::Test {};

typedef ::testing::Types<boost::int8_t, boost::int16_t, boost::int32_t, boost::int64_t, boost::uint8_t, boost::uint16_t,
                         boost::uint32_t, boost::uint64_t, half, float, double>
    NumericChannelDataTypes;

TYPED_TEST_CASE( OffsetInputWeightedSumCombine, NumericChannelDataTypes );

TYPED_TEST( OffsetInputWeightedSumCombine, Offset1Arity1 ) {
    const std::size_t arity = 1;

    float weights[2] = { 0, 1 };

    TypeParam a[arity + 2] = { 1, 2, 3 };
    TypeParam b[arity + 2] = { 4, 5, 6 };
    const char* data[] = { reinterpret_cast<const char*>( a ), reinterpret_cast<const char*>( b ) };

    TypeParam out[arity + 1] = { 0, 0 };

    const data_type_t dataType = channel_data_type_traits<TypeParam>::data_type();

    offset_input_channel_weighted_sum_combine_function_t f =
        offset_input_channel_weighted_sum_combine_function( dataType );

    f( weights, sizeof( TypeParam ), data, 2, arity, reinterpret_cast<char*>( out ) );

    EXPECT_EQ( 5, out[0] );
    EXPECT_EQ( 0, out[1] );
}

TYPED_TEST( OffsetInputWeightedSumCombine, Offset1Arity2 ) {
    const std::size_t arity = 2;

    float weights[2] = { 1, 0 };

    TypeParam a[arity + 2] = { 1, 2, 3, 4 };
    TypeParam b[arity + 2] = { 5, 6, 7, 8 };
    const char* data[] = { reinterpret_cast<const char*>( a ), reinterpret_cast<const char*>( b ) };

    TypeParam out[arity + 1] = { 0, 0, 0 };

    const data_type_t dataType = channel_data_type_traits<TypeParam>::data_type();

    offset_input_channel_weighted_sum_combine_function_t f =
        offset_input_channel_weighted_sum_combine_function( dataType );

    f( weights, sizeof( TypeParam ), data, 2, arity, reinterpret_cast<char*>( out ) );

    EXPECT_EQ( 2, out[0] );
    EXPECT_EQ( 3, out[1] );
    EXPECT_EQ( 0, out[2] );
}
