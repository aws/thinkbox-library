// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest-helper.h"
#include "gtest/gtest.h"

#include <frantic/particles/prt_metadata.hpp>

TEST( PRTMetadata, GetScannerTransformsNone ) {
    frantic::channels::property_map propertyMap;

    std::vector<frantic::graphics::transform4fd> scannerTransforms;

    frantic::particles::prt::get_scanner_transforms( propertyMap, scannerTransforms );

    EXPECT_TRUE( scannerTransforms.empty() );
}

TEST( PRTMetadata, GetScannerTransformsOneIdentity ) {
    frantic::channels::channel_map channelMap;
    channelMap.define_channel( _T("ScannerTransforms"), 16, frantic::channels::data_type_float32 );
    channelMap.end_channel_definition();

    frantic::channels::property_map propertyMap;
    propertyMap.set_channel_map( channelMap );
    char* charBuffer = propertyMap.get_channel_buffer( _T("ScannerTransforms") );
    float* buffer = reinterpret_cast<float*>( charBuffer );
    const float matrix[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    memcpy( buffer, matrix, sizeof( matrix ) );

    std::vector<frantic::graphics::transform4fd> scannerTransforms;

    frantic::particles::prt::get_scanner_transforms( propertyMap, scannerTransforms );

    ASSERT_EQ( scannerTransforms.size(), 1 );
    EXPECT_TRUE( scannerTransforms[0].is_identity() );
}

TEST( PRTMetadata, GetScannerTransformsTwoIdentity ) {
    frantic::channels::channel_map channelMap;
    channelMap.define_channel( _T("ScannerTransforms"), 32, frantic::channels::data_type_float64 );
    channelMap.end_channel_definition();

    frantic::channels::property_map propertyMap;
    propertyMap.set_channel_map( channelMap );
    char* charBuffer = propertyMap.get_channel_buffer( _T("ScannerTransforms") );
    double* buffer = reinterpret_cast<double*>( charBuffer );
    const double matrix[32] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
                                1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    memcpy( buffer, matrix, sizeof( matrix ) );

    std::vector<frantic::graphics::transform4fd> scannerTransforms;

    frantic::particles::prt::get_scanner_transforms( propertyMap, scannerTransforms );

    ASSERT_EQ( scannerTransforms.size(), 2 );
    EXPECT_TRUE( scannerTransforms[0].is_identity() );
    EXPECT_TRUE( scannerTransforms[1].is_identity() );
}

TEST( PRTMetadata, SetScannerTransformsNone ) {
    frantic::channels::property_map propertyMap;

    std::vector<frantic::graphics::transform4fd> inScannerTransforms;

    frantic::particles::prt::set_scanner_transforms( propertyMap, inScannerTransforms );

    std::vector<frantic::graphics::transform4fd> scannerTransforms;

    frantic::particles::prt::get_scanner_transforms( propertyMap, scannerTransforms );
    EXPECT_TRUE( scannerTransforms.empty() );
}

TEST( PRTMetadata, SetScannerTransformsClear ) {
    frantic::channels::channel_map channelMap;
    channelMap.define_channel( _T("ScannerTransforms"), 16, frantic::channels::data_type_float64 );
    channelMap.end_channel_definition();

    frantic::channels::property_map propertyMap;
    propertyMap.set_channel_map( channelMap );
    char* charBuffer = propertyMap.get_channel_buffer( _T("ScannerTransforms") );
    double* buffer = reinterpret_cast<double*>( charBuffer );
    const double matrix[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    memcpy( buffer, matrix, sizeof( matrix ) );

    std::vector<frantic::graphics::transform4fd> inScannerTransforms;
    frantic::particles::prt::set_scanner_transforms( propertyMap, inScannerTransforms );

    std::vector<frantic::graphics::transform4fd> scannerTransforms;

    frantic::particles::prt::get_scanner_transforms( propertyMap, scannerTransforms );
    EXPECT_TRUE( scannerTransforms.empty() );
}

TEST( PRTMetadata, SetScannerTransformsIdentity ) {
    frantic::channels::property_map propertyMap;

    std::vector<frantic::graphics::transform4fd> inScannerTransforms;
    inScannerTransforms.push_back( frantic::graphics::transform4fd::identity() );

    frantic::particles::prt::set_scanner_transforms( propertyMap, inScannerTransforms );

    std::vector<frantic::graphics::transform4fd> scannerTransforms;

    frantic::particles::prt::get_scanner_transforms( propertyMap, scannerTransforms );

    ASSERT_EQ( scannerTransforms.size(), 1 );
    EXPECT_TRUE( scannerTransforms[0].is_identity() );
}

TEST( PRTMetadata, SetScannerTransformsTwo ) {
    frantic::channels::property_map propertyMap;

    std::vector<frantic::graphics::transform4fd> inScannerTransforms;
    inScannerTransforms.push_back( frantic::graphics::transform4fd::identity() );
    inScannerTransforms.push_back( 2.0 * frantic::graphics::transform4fd::identity() );

    frantic::particles::prt::set_scanner_transforms( propertyMap, inScannerTransforms );

    std::vector<frantic::graphics::transform4fd> scannerTransforms;

    frantic::particles::prt::get_scanner_transforms( propertyMap, scannerTransforms );

    ASSERT_EQ( scannerTransforms.size(), 2 );
    EXPECT_TRUE( scannerTransforms[0].is_identity() );
    EXPECT_EQ( scannerTransforms[1], 2.0 * frantic::graphics::transform4fd::identity() );
}
