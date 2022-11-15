// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/scope_exit.hpp>

#include "gtest/gtest.h"

#include <frantic/locale/locale.hpp>
#include <frantic/misc/string_functions.hpp>
#include <frantic/strings/tstring.hpp>

TEST( ZeroPad, PositiveInteger ) {
    // The given integer should be converted to a string and padded with zeros until it contains
    // at least pad_to characters.
    const int pad_to = 4;
    const int number = 1;

    const frantic::tstring result = frantic::strings::zero_pad( number, pad_to );

    EXPECT_EQ( _T( "0001" ), result );
}

TEST( ZeroPad, BiggerThanPadding ) {
    // A number that is larger than the desired padding should not be truncated.
    const int pad_to = 4;
    const int number = 99999;

    const frantic::tstring result = frantic::strings::zero_pad( number, pad_to );

    EXPECT_EQ( _T( "99999" ), result );
}

TEST( ZeroPad, NegativeInteger ) {
    // When padding a negative number, we want the sign character to be included in the padding count.
    const int pad_to = 4;
    const int number = -1;

    const frantic::tstring result = frantic::strings::zero_pad( number, pad_to );

    EXPECT_EQ( _T( "-001" ), result );
}

TEST( ZeroPad, NegativeIntegerZeroSize ) {
    // Test the edge case where padding is zero.
    // This is mainly to ensure bugs aren't introduced due to swapping out
    // singed for unsigned types in the implementation and arithmetic operations then causing overflow
    const int pad_to = 0;
    const int number = -1;

    const frantic::tstring result = frantic::strings::zero_pad( number, pad_to );

    EXPECT_EQ( _T( "-1" ), result );
}
