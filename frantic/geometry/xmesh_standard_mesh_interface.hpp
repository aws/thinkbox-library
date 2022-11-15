// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/mesh_interface.hpp>

namespace frantic {
namespace geometry {

/**
 *  Create a read-only mesh_interface that contains only standard XMesh
 * channels.
 *
 *  To handle non-standard mesh channels:
 *
 * - Channels with a non-standard data type will be converted to a
 *   standard data type (if such type conversion is not possible, an
 *   exception will be thrown; TODO: something else?).
 * - Channels with a non-standard arity will have their values
 *   zero-padded or truncated to the standard arity.
 * - Channels with unknown names will be removed.
 *
 */
std::unique_ptr<mesh_interface> create_xmesh_standard_mesh_interface( const mesh_interface* meshInterface );

} // namespace geometry
} // namespace frantic
