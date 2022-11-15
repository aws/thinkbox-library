// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest-helper.h"
#include "gtest/gtest.h"

#include <frantic/particles/prt_metadata.hpp>

TEST( PropertyMap, DeleteProperty ) {
    frantic::channels::channel_map channelMap;
    channelMap.define_channel<frantic::graphics::vector3f>( _T("Mean") );
    channelMap.define_channel<frantic::graphics::boundbox3f>( _T("Extents") );
    channelMap.end_channel_definition();

    frantic::channels::property_map propertyMap;
    propertyMap.set_channel_map( channelMap );
    ASSERT_TRUE( propertyMap.has_property( _T("Mean") ) );
    ASSERT_TRUE( propertyMap.has_property( _T("Extents") ) );

    propertyMap.delete_property( _T("Extents") );
    ASSERT_FALSE( propertyMap.has_property( _T("Extents") ) );
    ASSERT_TRUE( propertyMap.has_property( _T("Mean") ) );

    propertyMap.delete_property( _T("Mean") );
    ASSERT_FALSE( propertyMap.has_property( _T("Mean") ) );
}

TEST( PropertyMap, AddProperty ) {
    frantic::channels::channel_map channelMap;
    channelMap.define_channel<frantic::graphics::vector3f>( _T("Mean") );
    channelMap.define_channel<frantic::graphics::boundbox3f>( _T("Extents") );
    channelMap.end_channel_definition();

    frantic::channels::property_map propertyMap;
    propertyMap.set_channel_map( channelMap );
    ASSERT_TRUE( propertyMap.has_property( _T("Mean") ) );
    ASSERT_TRUE( propertyMap.has_property( _T("Extents") ) );

    propertyMap.set_property<frantic::graphics::vector3f>( _T("StdDev"), frantic::graphics::vector3f( 1.0f ) );
    ASSERT_TRUE( propertyMap.has_property( _T("StdDev") ) );
    EXPECT_VECTOR3F_EQ( propertyMap.get_cvt<frantic::graphics::vector3f>( _T("StdDev") ),
                        frantic::graphics::vector3f( 1.0f ) );

    propertyMap.set_property<frantic::tstring>( _T("Interpretation"), _T("Point") );
    ASSERT_TRUE( propertyMap.has_property( _T("Interpretation") ) );
    EXPECT_EQ( propertyMap.get_cvt<frantic::tstring>( _T("Interpretation") ), _T("Point") );
}

TEST( PropertyMap, AddExistingProperty ) {
    frantic::channels::channel_map channelMap;
    channelMap.define_channel<frantic::graphics::vector3f>( _T("Mean") );
    channelMap.define_channel<frantic::graphics::boundbox3f>( _T("Extents") );
    channelMap.end_channel_definition();

    frantic::channels::property_map propertyMap;
    propertyMap.set_channel_map( channelMap );
    propertyMap.set_property<frantic::graphics::vector3f>( _T("StdDev"), frantic::graphics::vector3f( 1.0f ) );
    ASSERT_TRUE( propertyMap.has_property( _T("StdDev") ) );
    EXPECT_VECTOR3F_EQ( propertyMap.get_cvt<frantic::graphics::vector3f>( _T("StdDev") ),
                        frantic::graphics::vector3f( 1.0f ) );

    propertyMap.set_property<frantic::graphics::vector3f>( _T("StdDev"), frantic::graphics::vector3f( 2.0f ) );
    ASSERT_TRUE( propertyMap.has_property( _T("StdDev") ) );
    EXPECT_VECTOR3F_EQ( propertyMap.get_cvt<frantic::graphics::vector3f>( _T("StdDev") ),
                        frantic::graphics::vector3f( 2.0f ) );
}

TEST( PropertyMap, AddExistingPropertyThrows ) {
    frantic::channels::channel_map channelMap;
    channelMap.define_channel<frantic::graphics::vector3f>( _T("Mean") );
    channelMap.define_channel<frantic::graphics::boundbox3f>( _T("Extents") );
    channelMap.end_channel_definition();

    frantic::channels::property_map propertyMap;
    propertyMap.set_channel_map( channelMap );
    propertyMap.set_property<frantic::graphics::vector3f>( _T("StdDev"), frantic::graphics::vector3f( 1.0f ) );
    ASSERT_TRUE( propertyMap.has_property( _T("StdDev") ) );
    EXPECT_VECTOR3F_EQ( propertyMap.get_cvt<frantic::graphics::vector3f>( _T("StdDev") ),
                        frantic::graphics::vector3f( 1.0f ) );

    EXPECT_THROW( propertyMap.set_property<int>( _T("StdDev"), 2 ), std::runtime_error );
}