// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/xmesh_metadata.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace geometry {

/**
 * This function will write a mesh to disk. All channels in the mesh will be written.
 * @param path The filesystem location to write the xmesh files. The directory containing this path will get
 *              all the associated .xmdata files that store the channel information.
 * @param mesh The mesh to write to disk.
 * @param progress The progress logger to update while writing the file.
 */
void write_xmesh_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                            frantic::logging::progress_logger& progress, std::size_t threadCount = 1 );

/**
 * This function will write a mesh to disk. All channels in the mesh will be written.
 * @param path The filesystem location to write the xmesh files. The directory containing this path will get
 *              all the associated .xmdata files that store the channel information.
 * @param mesh The mesh to write to disk.
 * @param metadata The metadata to include in the xmesh file.
 * @param progress The progress logger to update while writing the file.
 */
void write_xmesh_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                            const xmesh_metadata& metadata, frantic::logging::progress_logger& progress,
                            std::size_t threadCount = 1 );

/**
 *  Write a mesh to the specified OBJ file.  This will write the
 * "TextureCoord" and "Normal" channels if they are present.
 *
 * @param path The filesystem location to write the OBJ file to.
 * @param mesh The mesh to write to disk.
 * @param progress The progress logger to update while writing the file.
 */
void write_obj_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                          frantic::logging::progress_logger& progress );

/**
 * Write a mesh to the specified ASCII STL file.  This will write the
 *    "Normal" channel.
 *
 * @param path The filesystem location to write the STL file to.
 * @param mesh The mesh to write to disk.
 * @param progress The progress logger to update while writing the file.
 */
void write_ascii_stl_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                                frantic::logging::progress_logger& progress );

/**
 * Write a mesh to the specified binary STL file.  This will write the
 *    "Normal" channel, and "Color" channel if present.
 *
 * @param path The filesystem location to write the STL file to.
 * @param mesh The mesh to write to disk.
 * @param isSolidView To tell if the color should be written in SolidView format or Magics.
 * @param progress The progress logger to update while writing the file.
 */
void write_binary_stl_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh, bool isSolidView,
                                 frantic::logging::progress_logger& progress );

/**
 * Write a mesh to the specified binary PLY file.  This will write the
 *    "TextureCoord" and "Color" channels if they are present.
 *
 * @param path The filesystem location to write the PLY file to.
 * @param mesh The mesh to write to disk.
 * @param progress The progress logger to update while writing the file.
 */
void write_ply_mesh_file( const frantic::tstring& path, const mesh_interface_ptr& mesh,
                          frantic::logging::progress_logger& progress );

} // namespace geometry
} // namespace frantic
