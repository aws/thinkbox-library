// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace geometry {
/**
 * Merge all duplicate vertices in a set of meshes, adjusting the face topology as necessary.
 *
 * Vertices are considered to be duplicates if the distance between them is not greater than the given tolerance.
 * Channels which are newly created due to the combination have their values defaulted to zero.
 *
 * @param meshes The meshes to weld.
 * @param tolerance Maximum distance between vertices which are considered duplicates.
 * @param progress Logger to track the progress of the weld.
 * @return The welded mesh.
 */
polymesh3_ptr weld_vertices( const std::vector<mesh_interface_ptr>& meshes, float tolerance );
polymesh3_ptr weld_vertices( const std::vector<const mesh_interface*>& meshes, float tolerance );
polymesh3_ptr weld_vertices( const std::vector<mesh_interface_ptr>& meshes, float tolerance,
                             frantic::logging::progress_logger& progress );
polymesh3_ptr weld_vertices( const std::vector<const mesh_interface*>& meshes, float tolerance,
                             frantic::logging::progress_logger& progress );

/**
 * Merge all matching boundary edges between the set of meshes, adjusting the face topology as necessary.
 *
 * Boundary edges are considered matching if their endpoints are within the given tolerance of one another, and they
 * have opposite winding orders. Each boundary edge can only be matched with one other edge. If the input meshes are all
 * manifold, then the resulting welded mesh will also be manifold. Channels which are newly created due to the
 * combination have their values defaulted to zero.
 *
 * @param meshes The meshes to weld.
 * @param tolerance Maximum distance between edge endpoints which are considered duplicates.
 * @param progress Logger to track the progress of the weld.
 * @return The welded mesh.
 */
polymesh3_ptr weld_boundary_edges( const std::vector<mesh_interface_ptr>& meshes, float tolerance );
polymesh3_ptr weld_boundary_edges( const std::vector<const mesh_interface*>& meshes, float tolerance );
polymesh3_ptr weld_boundary_edges( const std::vector<mesh_interface_ptr>& meshes, float tolerance,
                                   frantic::logging::progress_logger& progress );
polymesh3_ptr weld_boundary_edges( const std::vector<const mesh_interface*>& meshes, float tolerance,
                                   frantic::logging::progress_logger& progress );

/**
 * Combine all of the faces and vertices of the meshes into a single mesh.
 *
 * Channels which are newly created due to the combination have their values defaulted to zero.
 *
 * @param meshes The meshes to combine.
 * @param progress Logger to track the progress of the weld.
 * @return The combined mesh.
 */
polymesh3_ptr combine( const std::vector<mesh_interface_ptr>& meshes );
polymesh3_ptr combine( const std::vector<const mesh_interface*>& meshes );
polymesh3_ptr combine( const std::vector<mesh_interface_ptr>& meshes, frantic::logging::progress_logger& progress );
polymesh3_ptr combine( const std::vector<const mesh_interface*>& meshes, frantic::logging::progress_logger& progress );
} // namespace geometry
} // namespace frantic
