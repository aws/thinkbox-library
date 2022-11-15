// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polymesh3.hpp>

namespace frantic {
namespace geometry {

class polymesh3_interface : public mesh_interface {
  public:
    // TODO: take ownership of the mesh?
    static std::unique_ptr<polymesh3_interface> create_instance( frantic::geometry::polymesh3_ptr mesh );

    /**
     * @param mesh the mesh to create a read-only interface for.
     * @return a read-only interface for the mesh.
     */
    static std::unique_ptr<polymesh3_interface> create_const_instance( frantic::geometry::const_polymesh3_ptr mesh );

    virtual ~polymesh3_interface() {}
};

} // namespace geometry
} // namespace frantic
