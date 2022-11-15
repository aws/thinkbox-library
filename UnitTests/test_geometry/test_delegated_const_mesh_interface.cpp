// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/delegated_const_mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>

TEST( DelegatedConstMeshInterface, OneTriangle ) {
    frantic::geometry::trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 1 );

    frantic::geometry::mesh_interface::ptr_type a(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    frantic::geometry::mesh_interface::ptr_type b( new frantic::geometry::delegated_const_mesh_interface( a.get() ) );

    ASSERT_TRUE( bool( b ) );
    ASSERT_TRUE( b->is_valid() );

    EXPECT_TRUE( frantic::geometry::is_equal( a, b ) );
}
