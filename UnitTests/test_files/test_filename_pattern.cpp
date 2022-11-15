
// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/files/filename_sequence.hpp>

#include "gtest/gtest.h"

TEST( FilesTest, FilenamePattern ) {
    // This is just a regression test to make sure the int to string conversion
    // in frantic::files::filename_pattern::build_whole_frame_seq_string
    // still works correctly after changing it to use a safer implementation.
    frantic::files::filename_pattern pattern( _T( "seq_####.exr" ) ); 
    ASSERT_EQ(pattern[0], _T( "seq_0000.exr" ) );
    ASSERT_EQ(pattern[1], _T( "seq_0001.exr" ) );
    ASSERT_EQ(pattern[10], _T( "seq_0010.exr" ) );
    ASSERT_EQ(pattern[10000], _T( "seq_10000.exr" ) );
}

