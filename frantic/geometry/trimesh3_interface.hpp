// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/trimesh3.hpp>

#include <boost/move/move.hpp>

namespace frantic {
namespace geometry {

class trimesh3_interface : public mesh_interface {
  public:
    static std::unique_ptr<trimesh3_interface> create_instance( const trimesh3* mesh );
    static std::unique_ptr<trimesh3_interface> create_instance( BOOST_RV_REF( trimesh3 ) mesh );

    virtual ~trimesh3_interface() {}

    virtual const trimesh3& get_trimesh() const = 0;
};

} // namespace geometry
} // namespace frantic
