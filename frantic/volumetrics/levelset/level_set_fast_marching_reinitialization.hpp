// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_adjacency.hpp>

#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * This function reintializes a level set.
 *
 * @note This function allocates and uses the cached ris_adjacency in the rle_index_spec.
 *
 * @param  progressLogger   the logger will be used to display progress
 * @param  ris  A filled in rle_index_spec structure which corresponds to the populatedChannel and
 *signedDistanceChannel.
 * @param[in,out]  populatedChannel  The raw memory for the 'Populated' channel.
 *		This should be set to 1 for voxels with correct signed distance that can
 *		be used as a reinitialization source; set to 0 for voxels with invalid
 *		distances that should be reinitialized; and set to 2 for voxels that
 *		should be ignored during reinitialization (their distance is not used
 *		as a source and is not changed during reinitialization).
 *		On return, this will be 1 for voxels that were successfully reinitialized
 *		or were initially populated, and 0 for voxels that were not reinitialized
 *		or that were set to ignore.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  voxelLength  The length of one voxel.
 * @param  updatePopulatedInterfaceDistances  Determines whether the intial band should updated the distances of points
 *that are populated and on the interface.
 * @param  insideStoppingDistance  The width at which to stop reinitializing inside the interface (defaults to infinity,
 *which means no stopping width)
 * @param  outsideStoppingDistance  The width at which to stop reinitializing outside the interface (defaults to
 *infinity, which means no stopping width)
 */
void fast_marching_reinitialization( frantic::logging::progress_logger& progressLogger, const rle_index_spec& ris,
                                     unsigned char* populatedChannel, float* signedDistanceChannel,
                                     const float voxelLength, const bool updatePopulatedInterfaceDistances,
                                     const float insideStoppingDistance = std::numeric_limits<float>::infinity(),
                                     const float outsideStoppingDistance = std::numeric_limits<float>::infinity() );

/**
 * Extrapolate an unsigned distance field away from the initially populated
 * voxels.
 *
 * This function is based on fast_marching_reinitialization and
 * uses nearly the same procedure, but here the extrapolation is from the
 * initially populated voxels, and there is no check for zero crossings or
 * any changes in the initially populated band.
 *
 * The unsigned distance channel is expected to hold non-negative distances
 * in its populated voxels.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 *		   This should be set to 1 for voxels with correct distance that can
 *		   be used as an extrapolation source; set to 0 for voxels with invalid
 *		   distances that should be set by extrapolation; and set to 2 for voxels that
 *		   should be ignored during extrapolation (their distance is not used
 *		   as a source and is not changed during extrapolation).
 *		   All other values are reserved.
 *		   On return, this will be 1 for voxels that were successfully extrapolated into
 *		   or were initially populated, and 0 for voxels that were not extrapolated into
 *		   or that were set to ignore.
 * @param  signedDistanceChannel  The raw memory for the unsigned distance channel.
 * @param  voxelLength  The length of one voxel.
 * @param  stoppingDistance  The distance at which to stop reinitializing the distance field (defaults to infinity,
 *which means no stopping distance)
 */
void extrapolate_unsigned_distance( const frantic::volumetrics::levelset::ris_adjacency& adj,
                                    unsigned char* populatedChannel, float* unsignedDistanceChannel,
                                    const float voxelLength,
                                    const float outsideStoppingDistance = std::numeric_limits<float>::infinity() );

/**
 * Solve the quadratic equation used to calculate the first-order fast
 * marching solution to the isotropic eikonal equation.
 *
 * @param  voxelLength  The length of one voxel.
 * @param  phi1  The distance from the given voxel to the closest point on the interface along the X axis.  Use
 * std::numeric_limits<float>::infinity() if there is no close point.
 * @param  phi2  The distance from the given voxel to the closest point on the interface along the Y axis.
 * @param  phi3  The distance from the given voxel to the closest point on the interface along the Z axis.
 */
float solve_fast_marching_first_order_quad( const float voxelLength, const float phi1, const float phi2,
                                            const float phi3 );

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
