// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/graphics/alpha3f.hpp>

using frantic::graphics::alpha3f;

TEST( alpha3f, to_inverse ) {
    const alpha3f original( 0.25f, 0.5f, 0.75f );
    EXPECT_EQ( alpha3f( 0.75f, 0.5f, 0.25f ), original.to_inverse() );
}

TEST( alpha3f, subtraction ) {
    const alpha3f a( 1.f );
    const alpha3f b( 0.5f );

    EXPECT_EQ( b, a - b );
}
