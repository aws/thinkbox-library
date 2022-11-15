// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/geometry/trimesh3.hpp>

namespace frantic {
namespace geometry {

// This loads an .obj file into a trimesh3
void load_obj_mesh_file( const frantic::tstring& objFile, trimesh3& mesh );
void write_obj_mesh_file( const frantic::tstring& destFile, const trimesh3& mesh );

// Detect extension/format
void load_mesh_file( const frantic::tstring& srcFile, trimesh3& mesh );
void write_mesh_file( const frantic::tstring& destFile, const trimesh3& mesh );

/**
 *  Attempt to load an approximate bounding box from the specified mesh file.
 * This will succeed only with .xmesh files which have a \<boundbox\> element.
 * Such bounding boxes are intended for use as a proxy to avoid loading the
 * entire mesh or its vertices.
 *
 * @param srcFile name of the file from which to attempt to load an approximate
 *		bounding box.
 * @param[out] outBoundBox if the function returns true, then this is the
 *		bounding box which was loaded from the specified filename.
 * @param showWarnings if true and the boundbox data could not be interpreted
 *		correctly, then a warning will be sent to FF_LOG(warning).
 * @return true if a bounding box metadata could be read from the file,
 *		otherwise false.
 */
bool try_load_boundbox_metadata( const frantic::tstring& srcFile, boundbox3f& outBoundBox,
                                 const bool showWarnings = false );

} // namespace geometry
} // namespace frantic
