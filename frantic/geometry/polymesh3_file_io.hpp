// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/xmesh_metadata.hpp>
#include <frantic/geometry/xmesh_reader.hpp>
#include <frantic/graphics/units.hpp>

namespace frantic {
namespace geometry {

/**
 * This function will load a polymesh from the specified path.  It will detect the mesh file format based
 * on the path's extension.
 *
 * It will load all of the channels in the file, and it will be guaranteed to have a 'verts' channel,
 * probably with custom faces.
 *
 * @param path The filesystem location to load the mesh from
 * @return A new polymesh3 object with all the channels loaded.
 */
polymesh3_ptr load_polymesh_file( const frantic::tstring& path );

/**
 * This function will load a polymesh from the specified xmesh file. It will load all of the channels
 * in the file, and it will be guaranteed to have a 'verts' channel, probably with custom faces.
 * @param path The filesystem location to load the xmesh from
 * @return A new polymesh3 object with all the channels loaded.
 */
polymesh3_ptr load_xmesh_polymesh_file( const frantic::tstring& path );

/**
 * This function will load a polymesh from the specified xmesh_reader. It will load the channels
 * accepted by the channel_propagation_policy. It will be guaranteed to have a 'verts' channel.
 * @param reader The xmesh_reader to load the xmesh from
 * @param cpp The policy that determines which channels should be read.
 * @param loadFaces Whether to load the mesh channel faces.
 * @return A new polymesh3 object holding the loaded mesh data.
 */
polymesh3_ptr load_xmesh_polymesh_file( const frantic::geometry::xmesh_reader& reader,
                                        const frantic::channels::channel_propagation_policy& cpp,
                                        bool loadFaces = true );

/**
 * Load a Wavefront OBJ file from the specified path.  It will load the 'TextureCoord' and 'Normal'
 * channels if they are present.
 *
 * The resulting mesh is guaranteed to have a 'verts' channel, probably with custom faces.
 *
 * @param path The filesystem location to load the OBJ file from.
 * @return A new polymesh3 object with all the channels loaded.
 */
polymesh3_ptr load_obj_polymesh_file( const frantic::tstring& path );

/**
 * Load a STL file from the specified path.
 *
 * The resulting mesh is guaranteed to have a 'verts' and 'Normal' channel.
 *
 * @param path The filesystem location to load the STL file from.
 * @return A new polymesh3 object with all the channels loaded.
 */
polymesh3_ptr load_stl_polymesh_file( const frantic::tstring& path );

/**
 * Load a PLY file from the specified path.
 *
 * The resulting mesh is guaranteed to have a 'verts' channel.
 *
 * @param path The filesystem location to load the PLY file from.
 * @return A new polymesh3 object with all the channels loaded.
 */
polymesh3_ptr load_ply_polymesh_file( const frantic::tstring& path );

/**
 *  Write a polymesh3 to the specified path.  It will detect the file format
 * based on the path's extension.
 *
 * @param path The filesystem location to write the files to.
 * @param polymesh The polygonal mesh to write to disk.
 */
void write_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh );

/**
 *  Write a polymesh3 to the specified path.  It will detect the file format
 * based on the path's extension.
 *
 * @param path The filesystem location to write the files to.
 * @param polymesh The polygonal mesh to write to disk.
 * @param metadata The metadata to include in the file.
 */
void write_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh, const xmesh_metadata& metadata );

/**
 * This function will write a polymesh3 to disk. All channels in the polymesh will be written.
 * @param path The filesystem location to write the xmesh files. The directory containing this path will get
 *              all the associated .xmdata files that store the channel information.
 * @param polymesh The polygonal mesh to write to disk.
 */
void write_xmesh_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh );

/**
 * This function will write a polymesh3 to disk. All channels in the polymesh will be written.
 * @param path The filesystem location to write the xmesh files. The directory containing this path will get
 *              all the associated .xmdata files that store the channel information.
 * @param polymesh The polygonal mesh to write to disk.
 * @param metadata The metadata to include in the xmesh file.
 */
void write_xmesh_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh,
                                const xmesh_metadata& metadata );

/**
 *  Write a polymesh3 to the specified OBJ file.  This will write the
 * "TextureCoord" and "Normal" channels if they are present.
 *
 * @param path The filesystem location to write the OBJ file to.
 * @param polymesh The polygonal mesh to write to disk.
 */
void write_obj_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh );

/**
 *  Write a polymesh3 to the specified ASCII STL file.  This will write the
 * "Normal" channel.
 *
 *  @param path The filesystem location to write the STL file to.
 *  @param polymesh The polygonal mesh to write to disk.
 */
void write_ascii_stl_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh );

/**
 *  Write a polymesh3 to the specified binary STL file.  This will write the
 *  "Normal" channel, and "Color" channel if present.
 *
 *  @param path The filesystem location to write the STL file to.
 *  @param polymesh The polygonal mesh to write to disk.
 *  @param isSolidView To tell if the color should be written in SolidView format or Magics.
 */
void write_binary_stl_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh, bool isSolidView );

/**
 *  Write a polymesh3 to the specified binary PLY file.  This will write the
 * "TextureCoord" and "Color" channels if they are present.
 *
 *  @param path The filesystem location to write the PLY file to.
 *  @param polymesh The polygonal mesh to write to disk.
 */
void write_ply_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh );

} // namespace geometry
} // namespace frantic
