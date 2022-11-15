// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest-helper.h"

#include <frantic/os/environment.hpp>

using namespace frantic::os;

TEST( EnvironmentVariables, SetEnvironmentVariable ) {
    frantic::tstring key = _T("MY_KEY");
    frantic::tstring value = _T("MY_VALUE");

    set_environment_variable( key, value );
    set_environment_variable( key, value );
}

TEST( EnvironmentVariables, GetEnvironmentVariable ) {
    frantic::tstring firstKey = _T("MY_FIRST_KEY");
    frantic::tstring secondKey = _T("MY_SECOND_KEY");
    frantic::tstring firstValue = _T("MY_FIRST_VALUE");
    frantic::tstring secondValue = _T("MY_SECOND_VALUE");

    set_environment_variable( firstKey, firstValue );
    EXPECT_EQ( firstValue, get_environment_variable( firstKey ) );

    set_environment_variable( firstKey, secondValue );
    EXPECT_EQ( secondValue, get_environment_variable( firstKey ) );

    EXPECT_EQ( _T(""), get_environment_variable( secondKey ) );
}
