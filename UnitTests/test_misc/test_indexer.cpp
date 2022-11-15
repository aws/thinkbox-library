// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/misc/indexer.hpp>

TEST( Indexer, Simple2dIndexer ) {
    using namespace frantic;

    indexer<size_t, 2> idx = make_indexer<size_t>( 5, 5 );

    std::vector<bool> coverage( 5 * 5, false );

    size_t currentIndex = 0;

    for( size_t x0 = 0; x0 < 5; ++x0 ) {
        for( size_t x1 = 0; x1 < 5; ++x1 ) {
            size_t index = idx.address( x0, x1 );
            ASSERT_GT( coverage.size(), index );
            EXPECT_EQ( currentIndex, index );
            EXPECT_FALSE( coverage[index] );
            coverage[index] = true;
            ++currentIndex;
        }
    }
}

TEST( Indexer, Simple3dIndexer ) {
    using namespace frantic;

    indexer<size_t, 3> idx = make_indexer<size_t>( 5, 5, 5 );

    std::vector<bool> coverage( 5 * 5 * 5, false );

    size_t currentIndex = 0;

    for( size_t x0 = 0; x0 < 5; ++x0 ) {
        for( size_t x1 = 0; x1 < 5; ++x1 ) {
            for( size_t x2 = 0; x2 < 5; ++x2 ) {
                size_t index = idx.address( x0, x1, x2 );
                ASSERT_GT( coverage.size(), index );
                EXPECT_EQ( currentIndex, index );
                EXPECT_FALSE( coverage[index] );
                coverage[index] = true;
                ++currentIndex;
            }
        }
    }
}
