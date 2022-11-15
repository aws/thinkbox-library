// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/particle_array.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

/**
 * @param originalChannelMap a base channel_map to copy channels from.
 * @param volumeChannelName name to use for the volume channel.
 * @return a new channel_map that includes all of the channels in the
 *         originalChannelMap, plus the channels necessary for anisotropy
 *         calculations and anisotropic meshing.
 */
frantic::channels::channel_map
create_channel_map_with_anisotropy_channels( const frantic::channels::channel_map& originalChannelMap,
                                             const frantic::tstring& volumeChannelName = _T("__Volume") );

// these calls take a progress_logger parameter

void calculate_anisotropy( frantic::particles::particle_array& particles, float compactSupportScale,
                           float anisotropyWindowScale, float kr, std::size_t nEps,
                           frantic::logging::progress_logger& progressLogger );

void calculate_volume_with_anisotropic_kernel( frantic::particles::particle_array& particles,
                                               const frantic::tstring& volumeChannelName,
                                               frantic::logging::progress_logger& progressLogger );

void smooth_particle_positions( frantic::particles::particle_array& particles, float effectRadiusScale, float lambda,
                                frantic::logging::progress_logger& progressLogger );

// these calls take a shared_progress_logger_proxy parameter

void calculate_anisotropy( frantic::particles::particle_array& particles, float compactSupportScale,
                           float anisotropyWindowScale, float kr, std::size_t nEps,
                           frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );

void calculate_volume_with_anisotropic_kernel(
    frantic::particles::particle_array& particles, const frantic::tstring& volumeChannelName,
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );

void smooth_particle_positions( frantic::particles::particle_array& particles, float effectRadiusScale, float lambda,
                                frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );

struct smooth_particle_positions_params {
    frantic::particles::particle_array& particles;
    float effectRadiusScale;
    float lambda;
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger;

    smooth_particle_positions_params& operator=( const smooth_particle_positions_params& ); // not implemented

    smooth_particle_positions_params(
        frantic::particles::particle_array& particles, float effectRadiusScale, float lambda,
        frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );
};

void smooth_particle_positions( smooth_particle_positions_params& params );

struct calculate_volume_with_anisotropic_kernel_params {
    frantic::particles::particle_array& particles;
    frantic::tstring volumeChannelName;
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger;

    calculate_volume_with_anisotropic_kernel_params&
    operator=( const calculate_volume_with_anisotropic_kernel_params& ); // not implemented

    calculate_volume_with_anisotropic_kernel_params(
        frantic::particles::particle_array& particles, const frantic::tstring& volumeChannelName,
        frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );
};

void calculate_volume_with_anisotropic_kernel( calculate_volume_with_anisotropic_kernel_params& params );

struct calculate_anisotropy_params {
    frantic::particles::particle_array& particles;
    float effectRadiusScale;
    float anisotropyWindowRadiusScale;
    float kr;
    std::size_t ne;
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger;

    calculate_anisotropy_params& operator=( const calculate_anisotropy_params& );

    calculate_anisotropy_params( frantic::particles::particle_array& particles, float effectRadiusScale,
                                 float anisotropyWindowRadiusScale, float kr, std::size_t ne,
                                 frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );
};

void calculate_anisotropy( calculate_anisotropy_params& params );

} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic
