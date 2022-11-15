// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/misc/range_segmentation.hpp>

#include <boost/range/irange.hpp>

#include <vector>

TEST( MeshSegmentation, InitializeSingleChart ) {
    using namespace frantic;

    const size_t numFaces = 50;

    // all faces in single subset
    std::vector<size_t> segments( numFaces, 0 );

    range_segmentation seg( segments.begin(), segments.end(), 1 );

    EXPECT_EQ( 1, seg.get_num_subsets() );
    EXPECT_EQ( numFaces, seg.get_subset_size( 0 ) );
}

namespace {

void make_modulus_segmentation( size_t numFaces, size_t numCharts, frantic::range_segmentation& segmentation ) {
    std::vector<size_t> segments( numFaces, 0 );

    for( size_t i = 0; i < segments.size(); ++i ) {
        segments[i] = i % numCharts;
    }

    segmentation.assign( segments.begin(), segments.end(), numCharts );
}

} // namespace

TEST( MeshSegmentation, InitializeMultiChart ) {
    using namespace frantic;

    const size_t numFaces = 64;
    const size_t numCharts = 4;

    range_segmentation seg;
    make_modulus_segmentation( numFaces, numCharts, seg );

    EXPECT_EQ( 4, seg.get_num_subsets() );
    for( size_t subset = 0; subset < numCharts; ++subset ) {
        EXPECT_EQ( numFaces / numCharts, seg.get_subset_size( subset ) );
    }

    for( size_t i = 0; i < numFaces; ++i ) {
        EXPECT_EQ( i % numCharts, seg.get_face_subset( i ) );
    }
}

TEST( MeshSegmentation, MergeSubChart ) {
    using namespace frantic;

    const size_t numFaces = 64;
    const size_t numSubCharts = 4;
    const size_t numChartsPerSubChart = 4;

    range_segmentation inputSeg;
    make_modulus_segmentation( numFaces, numSubCharts, inputSeg );

    range_segmentation resultSeg( numFaces );

    for( size_t i = 0; i < numSubCharts; ++i ) {
        range_segmentation subSeg;
        make_modulus_segmentation( numFaces / numSubCharts, numChartsPerSubChart, subSeg );

        resultSeg.merge_from_subrange_segmentation( subSeg, inputSeg.subset_begin( i ), inputSeg.subset_end( i ) );
    }

    for( size_t i = 0; i < numFaces; ++i ) {
        EXPECT_EQ( ( ( i % numSubCharts ) * 4 ) + ( ( i / numSubCharts ) % numChartsPerSubChart ),
                   resultSeg.get_face_subset( i ) );
    }
}
