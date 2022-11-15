// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/trimesh3.hpp>

namespace frantic {
namespace geometry {

/**
 * Builds a polymesh3 from a trimesh.
 *
 * @param inMesh a trimesh to use as the template
 * @return a newly constructed polymesh3
 */
polymesh3_ptr trimesh3_to_polymesh3( const trimesh3& inMesh );

} // namespace geometry
} // namespace frantic
