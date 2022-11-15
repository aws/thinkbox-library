// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/dcel.hpp>

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/parameterization/parameterization_chart.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/trimesh3.hpp>

#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace geometry {

void trimesh3_to_dcel( const trimesh3& mesh, dcel& out );

void dcel_to_trimesh3( const dcel& input, const std::vector<frantic::graphics::vector3f>& vertices, trimesh3& outMesh );

void polymesh3_to_dcel( const_polymesh3_ptr polymesh, dcel& out );

void mesh_interface_to_dcel( const mesh_interface_ptr meshInterface, dcel& out );

void mesh_interface_to_dcel( const mesh_interface* meshInterface, dcel& out );

void mesh_channel_to_dcel( const mesh_channel* meshChannel, dcel& out );

void parameterization_chart_to_dcel( const parameterization::parameterization_chart& chart, dcel& out );

} // namespace geometry
} // namespace frantic
