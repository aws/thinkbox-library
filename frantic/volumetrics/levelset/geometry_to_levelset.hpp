// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * Converts the provided triangle mesh into a level set.  The output level set should be initialized to the desired
 * voxel space before this function is called.
 *
 * @note The result it produces is not an accurate distance function.  You will want to call a reinitialization function
 * if that's what you need.
 *
 * @param  geometry  The input triangle mesh which is converted into a level set.
 * @param  interfaceVoxelDistanceInside  The negative voxel distance from the interface towards the inside of the
 * object, completely to be made of defined voxels.  This value must be negative.
 * @param  interfaceVoxelDistanceOutside  The positive voxel distance from the interface towards the outside of the
 * object, completely to be made of defined voxels.  This value must be positive.
 * @param voxelBounds  The bounding box, in voxel space, to use as the outer bounds of the level set.  If this is empty,
 * this is derived from the input geometry.
 * @param  outLevelSet  The output level set.  When this function is called, this level set should already be
 * initialized with the desired voxel coordinate system.
 */
void convert_geometry_to_levelset( const geometry::trimesh3& geometry, float interfaceVoxelDistanceInside,
                                   float interfaceVoxelDistanceOutside, const frantic::graphics::boundbox3& voxelBounds,
                                   rle_level_set& outLevelSet );

/**
 * @overload
 *
 * This overload allows you to provide a world space bounding box instead of the voxel space bounding box.
 */
inline void convert_geometry_to_levelset( const geometry::trimesh3& geometry, float interfaceVoxelDistanceInside,
                                          float interfaceVoxelDistanceOutside,
                                          const frantic::graphics::boundbox3f& worldBounds,
                                          rle_level_set& outLevelSet ) {
    convert_geometry_to_levelset( geometry, interfaceVoxelDistanceInside, interfaceVoxelDistanceOutside,
                                  outLevelSet.get_voxel_coord_system().get_voxel_bounds( worldBounds ), outLevelSet );
}

/**
 * @overload
 *
 * This overload doesn't require a bounding box to be passed in.
 */
inline void convert_geometry_to_levelset( const geometry::trimesh3& geometry, float interfaceVoxelDistanceInside,
                                          float interfaceVoxelDistanceOutside, rle_level_set& outLevelSet ) {
    convert_geometry_to_levelset( geometry, interfaceVoxelDistanceInside, interfaceVoxelDistanceOutside,
                                  frantic::graphics::boundbox3::from_empty(), outLevelSet );
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
