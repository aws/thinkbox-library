// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/mesh_conversion.hpp>

#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>

namespace frantic {
namespace geometry {

polymesh3_ptr trimesh3_to_polymesh3( const trimesh3& inMesh ) {
    // its simpler to just use the existing method
    std::unique_ptr<trimesh3_interface> reference( trimesh3_interface::create_instance( &inMesh ) );
    return create_polymesh3( reference.get() );
}

} // namespace geometry
} // namespace frantic
