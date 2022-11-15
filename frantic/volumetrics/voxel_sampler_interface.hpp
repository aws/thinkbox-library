// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <boost/smart_ptr.hpp>

namespace frantic {
namespace volumetrics {

class voxel_sampler_interface {
  public:
    typedef boost::shared_ptr<voxel_sampler_interface> ptr_type;

    virtual void update_for_voxel( const frantic::graphics::vector3& voxel ) = 0;
    virtual bool get_next_position( frantic::graphics::vector3f& outPosition, float& outCompensationFactor ) = 0;
};

typedef voxel_sampler_interface::ptr_type voxel_sampler_interface_ptr;

} // namespace volumetrics
} // namespace frantic
