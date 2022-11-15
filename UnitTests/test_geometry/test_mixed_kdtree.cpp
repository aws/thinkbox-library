// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/mixed_kdtree.hpp>
#include <frantic/graphics/vector3f.hpp>

#include "gtest/gtest.h"

#include <algorithm>
#include <array>

TEST(GeometryTest, MixedKDTreePointData) {
    // We're interested in testing mixed_kdtree_point_data::set_my_data_from
    // which is called by the overload of the constructor that takes raw data.
    // We're primarily interested in testing against buffer overflows caused by off-by-one errors
    // when the mixed_kdtree_point_data is in its static data mode (size <= 80).

    // Create an std::vector<char> of length 80 that starts with 'a', ends with 'z', and is
    // otherwise filled with 'b'.
    constexpr std::size_t size{ 80 };
    std::array<char, size> data; 
    data[0] = 'a';
    std::fill( data.begin() + 1, data.end(), 'b' );
    data[size - 1] =  'z';

    // Initialize a mixed_kdtree_point_data such that set_my_data_from is called
    frantic::graphics::vector3f origin;
    frantic::geometry::mixed_kdtree_point_data point_data( origin, 0., origin, nullptr, data.data(), data.size() );

    // Get the primitive data
    const char* primitive_data = point_data.get_primitive_data();
    const std::size_t primitive_data_size = point_data.get_primitive_data_size();

    // Check to make sure the data we get back is as expected.
    ASSERT_EQ( primitive_data[0],'a' );
    for( std::size_t i = 1ull; i < primitive_data_size - 1; ++i ) {
        ASSERT_EQ( primitive_data[i], 'b' );
    }
    ASSERT_EQ( primitive_data[size - 1], 'z' );
}

