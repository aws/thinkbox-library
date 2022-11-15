// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/misc/range_segmentation.hpp>

#include <vector>

namespace frantic {
namespace geometry {

/**
 * Assocates to each face a positive integer, such that faces in the same connected component will share the same value
 *
 * @param mesh The mesh to run the search on
 * @param outFaceLabels output array containing the component id of each face
 */
void get_face_connected_components( const frantic::geometry::mesh_interface* mesh, std::vector<size_t>& outFaceLabels );
void get_face_connected_components( const frantic::geometry::mesh_interface_ptr mesh,
                                    std::vector<size_t>& outFaceLabels );

/**
 * Given a set of face component labels (as computed by 'get_face_connected_components'), gather all connected
 *components' faces
 * together into a set of lists.
 *
 * @param faceLabels the list of connected component labels, one per face
 * @param outFaceLists storage for the list of component faces
 * @param outListOffsets contains the offsets for each component's set of faces, such that
 *	outFaceLists[outListOffsets[c]..outListOffsets[c+1]] will contain the faces for component 'c'
 * @remark the length of outListOffsets will be 1 plus the number of connected components
 */
void group_face_components( const std::vector<size_t>& faceLabels, std::vector<size_t>& outFaceLists,
                            std::vector<size_t>& outListOffsets );

/**
 * Modify a set of face labels such that all subsets are entirely face-connected
 *
 * @param edgeStructure a DCEL representing the desired topology
 * @param faceLabels array to be modified containing the segment id of each face
 */
void separate_segmentation_components( const dcel& edgeStructure, std::vector<size_t>& faceLabels );

/**
 * Modify a face segmentation such that all subsets are entirely face-connected
 *
 * @param edgeStructure a DCEL representing the desired topology
 * @param initialSegmentation input face segmentation
 * @param outSegmentation output face segmentation
 */
void separate_segmentation_components( const mesh_interface* mesh, const range_segmentation& initialSegmentation,
                                       range_segmentation& outSegmentation );
void separate_segmentation_components( const mesh_interface* mesh, range_segmentation& inoutSegmentation );

} // namespace geometry
} // namespace frantic
