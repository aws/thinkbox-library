// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/mesh_interface.hpp>

namespace frantic {
namespace geometry {

/*
 * Transforms a mesh with a given transform4f
 * create_instance will return a transformed_mesh_interface of a mesh to access the transformed data
 *    The returned mesh is read only
 * Warning: changing the base mesh may cause undefined behaviour in the translated mesh as
 *    it does not make a deep copy of the base mesh
 */
class transformed_mesh_interface : public mesh_interface {
  public:
    static std::unique_ptr<transformed_mesh_interface> create_instance( frantic::geometry::mesh_interface_ptr mesh,
                                                                        const frantic::graphics::transform4f& xform );

    virtual ~transformed_mesh_interface() {}
};
} // namespace geometry
} // namespace frantic
