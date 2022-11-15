// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/tinyxml/frantic_tinyxml_utility.hpp>

TEST( TinyXMLUtility, InfiniteLoop ) {
    // Test to make sure calling remove_unsupported_xml_elements() doesn't cause an infinite loop.
    // This was know to happen previously when the first character was a '?' or '!' and there was no
    //'<' in the file.
    boost::filesystem::path filename( _T("TestInputs/infinite_loop.rsp") );
    // If the test fails, this call will cause an infinite loop.
    // The test will either timeout if run in the context of CI, or
    // otherwise need to be manually stopped.
    ASSERT_THROW( frantic::tinyxml::remove_unsupported_xml_elements( filename ), std::runtime_error );
}
