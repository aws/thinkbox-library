// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/camera.hpp>
#include <frantic/graphics/ray3f.hpp>
#include <frantic/graphics2d/framebuffer.hpp>
#include <frantic/rendering/framebuffer_cubeface.hpp>

namespace frantic {
namespace geometry {

class matte_geometry_interface {
  public:
    virtual void rebuild( float motionSegmentTime ) = 0;
    virtual bool is_visible( const frantic::graphics::ray3f& ray, double tMin, double tMax ) const = 0;
    virtual bool is_visible( const frantic::graphics::ray3f& ray, double tMin, double tMax, int& outLayer ) const = 0;
    virtual void generate_depth_map( const frantic::graphics::camera<float>& cam, float motionSegmentTime,
                                     frantic::graphics2d::framebuffer<float>& outBuffer, int numThreads = 0,
                                     bool forLight = false ) const = 0;
    virtual void generate_depth_map( const frantic::graphics::transform4f& camTM, float motionSegmentTime,
                                     frantic::rendering::framebuffer_cubeface<float>& outBuffer, int numThreads = 0,
                                     bool forLight = true ) const = 0;
};

} // namespace geometry
} // namespace frantic
