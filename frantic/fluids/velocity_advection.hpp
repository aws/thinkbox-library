// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/logging/logging_level.hpp>

#include <frantic/fluids/rle_staggered_vel_accessor.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

namespace frantic {
namespace fluids {

void ready_profiling();
void dump_profiling( std::basic_ostream<frantic::tchar>& out );

/**
 * This function type describes a function does ODE numerical integration on a velocity field. This is used for tracing
 * the velocity in a velocity field to a final position and querying the velocity at that location
 *
 * @param vcs the voxel coordinate system to convert to and from world locations
 * @param ris the rle index spec for the velocity field
 * @param staggeredVelocityAcc the channel accessor for the Staggered Velocity channel of the voxel field
 * @param worldLocation the starting world location
 * @param dt the time step to take
 * @param outPosition the final position after the trace
 */
typedef void ( *velocity_trace_function_t )(
    const frantic::volumetrics::voxel_coord_system& vcs, const frantic::volumetrics::levelset::rle_index_spec& ris,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& staggeredVelocityAcc,
    const frantic::graphics::vector3f& worldLocation, float dt, frantic::graphics::vector3f& finalPosition );

/**
 * Performs a 2nd order Runge Kutta Method trace in the velocity and returns the velocity at the final location and the
 * final location.
 *
 * @param vcs the voxel coordinate system to convert to and from world locations
 * @param ris the rle index spec for the velocity field
 * @param staggeredVelocityAcc the channel accessor for the Staggered Velocity channel of the voxel field
 * @param worldLocation the starting world location
 * @param dt the time step to take
 * @param outPosition the final position after the trace
 */
void trace_runge_kutta_2(
    const frantic::volumetrics::voxel_coord_system& vcs, const frantic::volumetrics::levelset::rle_index_spec& ris,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& staggeredVelocityAcc,
    const frantic::graphics::vector3f& worldLocation, float dt, frantic::graphics::vector3f& finalPosition );

/**
 * Performs a 2nd order Runge Kutta Method trace in the velocity and returns the velocity at the final location and the
 * final location.
 *
 * @param vcs the voxel coordinate system to convert to and from world locations
 * @param ris the rle index spec for the velocity field
 * @param velocityAcc the channel accessor for the Velocity channel of the voxel field
 * @param worldLocation the starting world location
 * @param dt the time step to take
 * @param outPosition the final position after the trace
 */
void trace_runge_kutta_2_unstaggered(
    const frantic::volumetrics::voxel_coord_system& vcs, const frantic::volumetrics::levelset::rle_index_spec& ris,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& velocityAcc,
    const frantic::graphics::vector3f& worldLocation, float dt, frantic::graphics::vector3f& finalPosition );

/**
 * Performs a 3nd order Runge Kutta Method trace in the velocity and returns the velocity at the final location and the
 * final location.
 *
 * @param vcs the voxel coordinate system to convert to and from world locations
 * @param ris the rle index spec for the velocity field
 * @param staggeredVelocityAcc the channel accessor for the Staggered Velocity channel of the voxel field
 * @param worldLocation the starting world location
 * @param dt the time step to take
 * @param outPosition the final position after the trace
 */
void trace_runge_kutta_3(
    const frantic::volumetrics::voxel_coord_system& vcs, const frantic::volumetrics::levelset::rle_index_spec& ris,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& staggeredVelocityAcc,
    const frantic::graphics::vector3f& worldLocation, float dt, frantic::graphics::vector3f& finalPosition );

/**
 * Performs a 4th order Runge Kutta Method trace in the velocity and returns the velocity at the final location and the
 * final location.
 *
 * @param vcs the voxel coordinate system to convert to and from world locations
 * @param ris the rle index spec for the velocity field
 * @param staggeredVelocityAcc the channel accessor for the Staggered Velocity channel of the voxel field
 * @param worldLocation the starting world location
 * @param dt the time step to take
 * @param outPosition the final position after the trace
 */
void trace_runge_kutta_4(
    const frantic::volumetrics::voxel_coord_system& vcs, const frantic::volumetrics::levelset::rle_index_spec& ris,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& staggeredVelocityAcc,
    const frantic::graphics::vector3f& worldLocation, float dt, frantic::graphics::vector3f& finalPosition );

void compute_curl_channel( rle_voxel_field& source, const std::string srcChannel,
                           const std::string& curlChannelToCreate );
void compute_convective_derivative_channel( rle_voxel_field& source, const std::string srcChannel,
                                            const std::string& curlChannelToCreate );

/**
 * This function advects a velocity field using a WEN03 lookup and RK3 for the backtrace. It returns a new field with
 * the velocities advected based on the time step.
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement)
 * @param[out] result the resultling velocity field after the advection
 * @param velocityField the velocity field in which the source will be advected.
 * @param dt the time step taken
 */

void velocity_weno_advect( rle_voxel_field& source, rle_voxel_field& result, rle_voxel_field& velocityField, float dt );

/**
 *  This function advects a velocity field using a WEN03 lookup and RK3 for the backtrace. It returns a new field with
 * the velocities advected based on the time step. It does not perform any multithreading and used mostly for debugging.
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement)
 * @param[out] result the resultling velocity field after the advection
 * @param velocityField the velocity field in which the source will be advected.
 * @param dt the time step taken
 */
void serial_velocity_weno_advect( rle_voxel_field& source, rle_voxel_field& result, rle_voxel_field& velocityField,
                                  float dt );

/**
 * This function advects a velocity field, returning a new field with the velocities advected based on the time step. It
 * does not perform any mutlithreading and used mostly for debugging.
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement)
 * @param[out] result the resultling velocity field after the advection
 * @param velocityField the velocity field in which the source will be advected.
 * @param dt the time step taken
 */
void serial_velocity_advect( const rle_voxel_field& source, rle_voxel_field& result,
                             const rle_voxel_field& velocityField, float dt );

/**
 * This function advects a velocity field, returning a new field with the velocities advected based on the time step.
 *
 * This performs the same operation as the other velocity_advect function, but this version is faster because it does
 *not search for the maximum velocity in velocityField.
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement)
 * @param[out] result the resulting velocity field after the advection
 * @param velocityField the velocity field in which the source will be advected.
 * @param maxVoxelMotion the maximum motion of the velocityField during the
 *		time step ( = max( abs( StaggeredVelocity ) ) * dt ), measured in voxel
 *		lengths.  This is the CFL condition number for the velocityField with
 *		time step dt.  Using a value that is too small will cause a runtime
 *		error.  Using a value that is too large will slow the function down.
 * @param dt the time step taken
 */
void velocity_advect( const rle_voxel_field& source, rle_voxel_field& result, const rle_voxel_field& velocityField,
                      float maxVoxelMotion, float dt );

/**
 * This function advects a velocity field, returning a new field with the velocities advected based on the time step
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement)
 * @param[out] result the resultling velocity field after the advection
 * @param velocityField the velocity field in which the source will be advected.
 * @param dt the time step taken
 */
void velocity_advect( const rle_voxel_field& source, rle_voxel_field& result, const rle_voxel_field& velocityField,
                      float dt );

/**
 *  Advect a velocity field using MacCormack advection.
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement).
 * @param[out] result the resulting velocity field after advection.
 * @param velocityField the velocity field in which the source will be advected.
 * @param dt the time step taken.
 */
void velocity_maccormack_advect( const rle_voxel_field& source, rle_voxel_field& result,
                                 const rle_voxel_field& velocityField, float dt );

/**
 *  Advect a velocity field using BFECC advection.
 *
 *  This is intended to implement the BFECC algorithm described in:
 *  Andrew Selle et al.  "An Unconditionally Stable MacCormack Method".
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement).
 * @param[out] result the resulting velocity field after advection.
 * @param velocityField the velocity field in which the source will be advected.
 * @param dt the time step taken.
 */
void velocity_bfecc_advect( const rle_voxel_field& source, rle_voxel_field& result,
                            const rle_voxel_field& velocityField, float dt );

/**
 *  Advect a velocity field using BFECC advection.
 *
 *  I suspect this implementation was changed in a copy-paste mistake,
 * but people like its look, so I'm leaving it around.
 *
 * @param source the source velocity field (typically divergence-free, though that is not a requirement).
 * @param[out] result the resulting velocity field after advection.
 * @param velocityField the velocity field in which the source will be advected.
 * @param dt the time step taken.
 */
void velocity_modified_bfecc_advect( const rle_voxel_field& source, rle_voxel_field& result,
                                     const rle_voxel_field& velocityField, float dt );

} // namespace fluids
} // namespace frantic
