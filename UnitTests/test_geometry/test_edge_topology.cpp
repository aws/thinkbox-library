// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/geometry/edge_topology.hpp>

#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <utilities/mesh_generators.hpp>

TEST( EdgeTopology, GetBoundaryEdgesManifold ) {
    using namespace frantic::geometry;

    polymesh3_ptr mesh = make_cube_polymesh();

    std::vector<topology::half_edge> boundaries;
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    get_boundary_edges( meshInterface.get(), boundaries );

    ASSERT_EQ( 0, boundaries.size() );
}

TEST( EdgeTopology, GetBoundaryEdgesManifoldWithBoundaries ) {
    using namespace frantic::geometry;

    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, 1, 0 );
    builder.add_vertex( 1, 0, 0 );

    const int first[3] = { 0, 2, 1 };
    builder.add_polygon( first, 3 );
    const int second[3] = { 1, 2, 3 };
    builder.add_polygon( second, 3 );
    const int third[3] = { 0, 3, 2 };
    builder.add_polygon( third, 3 );

    polymesh3_ptr mesh = builder.finalize();

    std::vector<topology::half_edge> boundaries;
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    get_boundary_edges( meshInterface.get(), boundaries );

    ASSERT_EQ( 3, boundaries.size() );
}

TEST( EdgeTopology, GetBoundaryEdgesNonManifold ) {
    using namespace frantic::geometry;

    polymesh3_builder builder;
    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 0, 0, 1 );
    builder.add_vertex( 0, 1, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 0, 1, 1 );

    const int first[3] = { 0, 2, 1 };
    builder.add_polygon( first, 3 );
    const int second[3] = { 1, 2, 3 };
    builder.add_polygon( second, 3 );
    const int third[3] = { 0, 3, 2 };
    builder.add_polygon( third, 3 );
    const int fourth[3] = { 1, 4, 2 };
    builder.add_polygon( fourth, 3 );

    polymesh3_ptr mesh = builder.finalize();

    std::vector<topology::half_edge> boundaries;
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    get_boundary_edges( meshInterface.get(), boundaries );

    ASSERT_EQ( 5, boundaries.size() );
}

TEST( EdgeTopology, GetBoundaryEdgesEmpty ) {
    using namespace frantic::geometry;

    polymesh3_builder builder;
    polymesh3_ptr mesh = builder.finalize();

    std::vector<topology::half_edge> boundaries;
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
    get_boundary_edges( meshInterface.get(), boundaries );

    ASSERT_EQ( 0, boundaries.size() );
}
