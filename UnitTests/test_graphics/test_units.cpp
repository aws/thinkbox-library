// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics/units.hpp>

using std::numeric_limits;
using namespace frantic::graphics;
using frantic::graphics::length_unit::format_autoscale;

TEST( Units, FormatAutoScale ) {
    double nan = numeric_limits<double>::quiet_NaN();
    double inf = numeric_limits<double>::infinity();
    // Infinity/NaN cases
    EXPECT_EQ( _T("nan"), format_autoscale( nan, length_unit::invalid, 2 ) );
    EXPECT_EQ( _T("inf"), format_autoscale( inf, length_unit::invalid, 2 ) );
    EXPECT_EQ( _T("-inf"), format_autoscale( -inf, length_unit::invalid, 2 ) );
    EXPECT_EQ( _T("nan"), format_autoscale( nan, length_unit::meters, 2 ) );
    EXPECT_EQ( _T("inf"), format_autoscale( inf, length_unit::meters, 2 ) );
    EXPECT_EQ( _T("-inf"), format_autoscale( -inf, length_unit::meters, 2 ) );
    // Without a unit specified
    EXPECT_EQ( _T("3.14"), format_autoscale( 3.141, length_unit::invalid, 3 ) );
    EXPECT_EQ( _T("3.141"), format_autoscale( 3.141, length_unit::invalid, 4 ) );
    // Reasonable values should stick to the specified unit
    EXPECT_EQ( _T("3.14 mm"), format_autoscale( 3.141, length_unit::millimeters, 3 ) );
    EXPECT_EQ( _T("3.14 cm"), format_autoscale( 3.141, length_unit::centimeters, 3 ) );
    EXPECT_EQ( _T("3.14 in"), format_autoscale( 3.141, length_unit::inches, 3 ) );
    EXPECT_EQ( _T("3.14 ft"), format_autoscale( 3.141, length_unit::feet, 3 ) );
    EXPECT_EQ( _T("3.14 yd"), format_autoscale( 3.141, length_unit::yards, 3 ) );
    // Large/small values should autoscale to different unit
    EXPECT_EQ( _T("3.14 m"), format_autoscale( 3141, length_unit::millimeters, 3 ) );
    EXPECT_EQ( _T("3.14 km"), format_autoscale( 3141000, length_unit::millimeters, 3 ) );
    EXPECT_EQ( _T("3.14 mm"), format_autoscale( 0.003141, length_unit::meters, 3 ) );
    EXPECT_EQ( _T("3.14 m"), format_autoscale( 0.003141, length_unit::kilometers, 3 ) );
    EXPECT_EQ( _T("3.14 mm"), format_autoscale( 0.000003141, length_unit::kilometers, 3 ) );
    EXPECT_EQ( _T("3.14 m"), format_autoscale( 3141, length_unit::millimeters, 3 ) );
    // note: The std::ostream formatting doesn't quite give enough control, might want to roll our own...
    EXPECT_EQ( _T("1 in"), format_autoscale( 1.0 / 12.0, length_unit::feet, 2 ) );
    EXPECT_EQ( _T("3 in"), format_autoscale( 3.0 / 12.0 / 3.0, length_unit::yards, 2 ) );
}