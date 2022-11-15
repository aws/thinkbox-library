// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"
#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/dcel_mesh_interface.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>

using namespace frantic::geometry;
namespace {

polymesh3_ptr make_single_vertex_polymesh() {
    // Create the polymesh3 to test
    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );

    polymesh3_ptr polymesh = builder.finalize();
    EXPECT_TRUE( bool( polymesh ) );

    return polymesh;
}
} // anonymous namespace

TEST( DCELInterface, VertAdjacency ) {
    using namespace boost;

    polymesh3_ptr triangle = make_triangle_polymesh();
    polymesh3_ptr vertex = make_single_vertex_polymesh();
    dcel adjacency;
    polymesh3_to_dcel( triangle, adjacency );
    mesh_interface_ptr connected =
        boost::shared_ptr<dcel_mesh_interface>( new dcel_mesh_interface( boost::move( adjacency ) ) );
    frantic::geometry::vertex_iterator vIt;
    EXPECT_TRUE( connected->init_vertex_iterator( vIt, 0 ) );

    polymesh3_to_dcel( vertex, adjacency );
    mesh_interface_ptr unconnected =
        boost::shared_ptr<dcel_mesh_interface>( new dcel_mesh_interface( boost::move( adjacency ) ) );

    EXPECT_FALSE( unconnected->init_vertex_iterator( vIt, 0 ) );
}
