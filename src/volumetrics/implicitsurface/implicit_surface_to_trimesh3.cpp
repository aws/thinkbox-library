// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#pragma warning( push )
#pragma warning( disable : 4512 4100 )
#include <tbb/parallel_reduce.h>
#pragma warning( pop )

#include <boost/bind.hpp>

#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/geometry/trimesh3_degeneracy_removal.hpp>
#include <frantic/math/reconstruction_filters.hpp>
#include <frantic/particles/particle_grid_tree.hpp>
#include <frantic/particles/sph/sph_kernel_functions.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>
#include <frantic/volumetrics/implicitsurface/level_set_implicit_surface_policies.hpp>
#include <frantic/volumetrics/implicitsurface/marching_cubes_table.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>

#include <frantic/volumetrics/levelset/rle_run_iterator.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::geometry;
using namespace frantic::particles;
using namespace frantic::volumetrics::levelset;

void frantic::volumetrics::levelset::convert_levelset_to_ris_debug_trimesh3( const rle_level_set& levelSet,
                                                                             float cubeSize, trimesh3& outMesh ) {
    frantic::logging::null_progress_logger nullLogger;
    frantic::volumetrics::levelset::convert_levelset_to_ris_debug_trimesh3( levelSet, cubeSize, outMesh, nullLogger );
}

void frantic::volumetrics::levelset::convert_levelset_to_ris_debug_trimesh3(
    const rle_level_set& levelSet, float cubeSize, trimesh3& outMesh,
    frantic::logging::progress_logger& progressLogger ) {
    trimesh3 result, tempMesh;
    frantic::tstring colorChannelName = _T("Color");

    // Initialize the tempMesh to a box.  For performance reasons, this box is never reallocated, it just has its
    // vertices moved about.
    tempMesh.set_to_box();

    // Add a color channel, and get an accessor to it
    tempMesh.add_vertex_channel<color3f>( colorChannelName );
    trimesh3_vertex_channel_accessor<color3f> ca = tempMesh.get_vertex_channel_accessor<color3f>( colorChannelName );

    const voxel_coord_system& vcs = levelSet.get_voxel_coord_system();
    const rle_index_spec& ris = levelSet.get_rle_index_spec();

    // Do a red-purple-blue gradient across 5 voxels
    float insideDistForColor = -2.5f * vcs.voxel_length();
    float outsideDistForColor = 2.5f * vcs.voxel_length();

    boundbox3 outerBounds = ris.outer_bounds();
    int zMin = outerBounds.zminimum(), zMax = outerBounds.zmaximum(), yMin = outerBounds.yminimum(),
        yMax = outerBounds.ymaximum();
    int numIntervals = ( yMax - yMin + 1 ) * ( zMax - zMin + 1 );
    int currInterval = 0;

    for( int z = zMin; z <= zMax; ++z ) {
        for( int y = yMin; y <= yMax; ++y ) {
            for( rle_run_iterator i( ris, ris.y_to_b( y ), ris.z_to_c( z ) ), ie; i != ie; ++i ) {
                int dataIndex = i.get_data_index();
                if( dataIndex < 0 ) {
                    // Get a bounding box to represent this run
                    boundbox3f bounds( vcs.get_world_voxel_center( vector3( i.get_xmin(), y, z ) ),
                                       vcs.get_world_voxel_center( vector3( i.get_xmax(), y, z ) ) );
                    bounds.expand( 0.5f * cubeSize * vcs.voxel_length() );

                    tempMesh.set_existing_box_verts( bounds );
                    color3f col( 0 );
                    if( dataIndex == -1 ) {
                        col.r = 1;
                        col.g = 1;
                    } else {
                        col.g = 1;
                    }
                    for( size_t index = 0, indexEnd = ca.size(); index != indexEnd; ++index ) {
                        ca[index] = col;
                    }

                    result.combine( tempMesh );
                } else {
                    vector3 coord( i.get_xmin(), y, z );
                    int xEnd = coord.x + i.get_xsize();
                    for( ; coord.x != xEnd; ++coord.x, ++dataIndex ) {
                        boundbox3f bounds( vcs.get_world_voxel_center( coord ) );
                        bounds.expand( 0.5f * cubeSize * vcs.voxel_length() );

                        float d = levelSet[dataIndex];

                        tempMesh.set_existing_box_verts( bounds );
                        color3f col( 0 );
                        if( d < 0 ) {
                            col.r = ( insideDistForColor - d ) / insideDistForColor;
                            if( col.r < 0 )
                                col.r = 0;
                            col.b = 1;
                        } else {
                            col.r = 1;
                            col.b = ( outsideDistForColor - d ) / outsideDistForColor;
                            if( col.b < 0 )
                                col.b = 0;
                        }
                        for( size_t index = 0, indexEnd = ca.size(); index != indexEnd; ++index ) {
                            ca[index] = col;
                        }
                        result.combine( tempMesh );
                    }
                }
            }
            ++currInterval;
            progressLogger.update_progress( currInterval, numIntervals );
        }
    }

    result.swap( outMesh );
}

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

// without progress logger
void union_of_spheres_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                     const frantic::channels::channel_propagation_policy& cpp,
                                                     float maximumParticleRadius,
                                                     float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                     const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                     frantic::geometry::trimesh3& outMesh ) {
    particle_union_of_spheres_is_policy isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale,
                                             implicitThreshold, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh );
}

// with progress logger
void union_of_spheres_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                     const frantic::channels::channel_propagation_policy& cpp,
                                                     float maximumParticleRadius,
                                                     float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                     const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                     frantic::geometry::trimesh3& outMesh,
                                                     frantic::logging::progress_logger& progressLogger ) {
    particle_union_of_spheres_is_policy isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale,
                                             implicitThreshold, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh, progressLogger );
}

void union_of_spheres_convert_sparse_particles_to_trimesh3(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maximumParticleRadius, float particleRadiusToEffectRadiusScale, float implicitThreshold,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    frantic::logging::progress_logger& progressLogger ) {
    frantic::volumetrics::implicitsurface::shared_progress_logger_adapter uiAdapter( progressLogger );
    frantic::volumetrics::implicitsurface::union_of_spheres_sparse_params params(
        particles, cpp, maximumParticleRadius, particleRadiusToEffectRadiusScale, implicitThreshold, meshingVCS,
        vertexRefinement, outMesh, *uiAdapter.get_proxy() );
    boost::function<void( void )> f(
        boost::bind( frantic::volumetrics::implicitsurface::union_of_spheres_convert_sparse_particles_to_trimesh3,
                     boost::ref( params ) ) );
    uiAdapter.run( f );
}

void union_of_spheres_convert_sparse_particles_to_trimesh3(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maximumParticleRadius, float particleRadiusToEffectRadiusScale, float implicitThreshold,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    shared_progress_logger_proxy& progressLogger ) {
    particle_union_of_spheres_is_policy isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale,
                                             implicitThreshold, meshingVCS, vertexRefinement );
    convert_particle_implicit_surface_to_trimesh3( isp, cpp, outMesh, progressLogger );
}

union_of_spheres_sparse_params::union_of_spheres_sparse_params(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maximumParticleRadius, float particleRadiusToEffectRadiusScale, float implicitThreshold,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    shared_progress_logger_proxy& progressLogger )
    : particles( particles )
    , cpp( cpp )
    , maximumParticleRadius( maximumParticleRadius )
    , particleRadiusToEffectRadiusScale( particleRadiusToEffectRadiusScale )
    , implicitThreshold( implicitThreshold )
    , meshingVCS( meshingVCS )
    , vertexRefinement( vertexRefinement )
    , outMesh( outMesh )
    , progressLogger( progressLogger ) {}

void union_of_spheres_convert_sparse_particles_to_trimesh3( union_of_spheres_sparse_params& params ) {
    union_of_spheres_convert_sparse_particles_to_trimesh3(
        params.particles, params.cpp, params.maximumParticleRadius, params.particleRadiusToEffectRadiusScale,
        params.implicitThreshold, params.meshingVCS, params.vertexRefinement, params.outMesh, params.progressLogger );
}

// without progress logger
void metaball_convert_particles_to_trimesh3( particle_grid_tree& particles,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             float maximumParticleRadius, float particleRadiusToEffectRadiusScale,
                                             float implicitThreshold, const voxel_coord_system& meshingVCS,
                                             int vertexRefinement, trimesh3& outMesh ) {
    particle_metaball_is_policy isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale,
                                     implicitThreshold, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh );
}

// with progress logger
void metaball_convert_particles_to_trimesh3( particle_grid_tree& particles,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             float maximumParticleRadius, float particleRadiusToEffectRadiusScale,
                                             float implicitThreshold, const voxel_coord_system& meshingVCS,
                                             int vertexRefinement, trimesh3& outMesh,
                                             frantic::logging::progress_logger& progressLogger ) {
    particle_metaball_is_policy isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale,
                                     implicitThreshold, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh, progressLogger );
}

void metaball_convert_sparse_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                    const frantic::channels::channel_propagation_policy& cpp,
                                                    float maximumParticleRadius,
                                                    float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                    const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                    frantic::geometry::trimesh3& outMesh,
                                                    frantic::logging::progress_logger& progressLogger ) {
    frantic::volumetrics::implicitsurface::shared_progress_logger_adapter uiAdapter( progressLogger );
    frantic::volumetrics::implicitsurface::metaball_sparse_params params(
        particles, cpp, maximumParticleRadius, particleRadiusToEffectRadiusScale, implicitThreshold, meshingVCS,
        vertexRefinement, outMesh, *uiAdapter.get_proxy() );
    boost::function<void( void )> f( boost::bind(
        frantic::volumetrics::implicitsurface::metaball_convert_sparse_particles_to_trimesh3, boost::ref( params ) ) );
    uiAdapter.run( f );
}

void metaball_convert_sparse_particles_to_trimesh3( particle_grid_tree& particles,
                                                    const frantic::channels::channel_propagation_policy& cpp,
                                                    float maximumParticleRadius,
                                                    float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                    const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                    trimesh3& outMesh, shared_progress_logger_proxy& progressLogger ) {
    particle_metaball_is_policy isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale,
                                     implicitThreshold, meshingVCS, vertexRefinement );
    convert_particle_implicit_surface_to_trimesh3( isp, cpp, outMesh, progressLogger );
}

metaball_sparse_params::metaball_sparse_params( frantic::particles::particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float maximumParticleRadius, float particleRadiusToEffectRadiusScale,
                                                float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                                shared_progress_logger_proxy& progressLogger )
    : particles( particles )
    , cpp( cpp )
    , maximumParticleRadius( maximumParticleRadius )
    , particleRadiusToEffectRadiusScale( particleRadiusToEffectRadiusScale )
    , implicitThreshold( implicitThreshold )
    , meshingVCS( meshingVCS )
    , vertexRefinement( vertexRefinement )
    , outMesh( outMesh )
    , progressLogger( progressLogger ) {}

void metaball_convert_sparse_particles_to_trimesh3( metaball_sparse_params& params ) {
    metaball_convert_sparse_particles_to_trimesh3(
        params.particles, params.cpp, params.maximumParticleRadius, params.particleRadiusToEffectRadiusScale,
        params.implicitThreshold, params.meshingVCS, params.vertexRefinement, params.outMesh, params.progressLogger );
}

// without progress logger
void zhu_bridson_convert_particles_to_trimesh3( particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float maxParticleRadius, float effectRadius,
                                                float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
                                                const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                trimesh3& outMesh ) {
    particle_zhu_bridson_is_policy isp( particles, maxParticleRadius, effectRadius, lowDensityTrimmingDensity,
                                        lowDensityTrimmingStrength, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh );
}

// with progress logger
void zhu_bridson_convert_particles_to_trimesh3( particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float maxParticleRadius, float effectRadius,
                                                float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
                                                const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                trimesh3& outMesh, frantic::logging::progress_logger& progressLogger ) {
    particle_zhu_bridson_is_policy isp( particles, maxParticleRadius, effectRadius, lowDensityTrimmingDensity,
                                        lowDensityTrimmingStrength, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh, progressLogger );
}

void zhu_bridson_convert_sparse_particles_to_trimesh3(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maxParticleRadius, float effectRadiusScale, float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    frantic::logging::progress_logger& progressLogger ) {
    frantic::volumetrics::implicitsurface::shared_progress_logger_adapter uiAdapter( progressLogger );
    frantic::volumetrics::implicitsurface::zhu_bridson_sparse_params params(
        particles, cpp, maxParticleRadius, effectRadiusScale, lowDensityTrimmingDensity, lowDensityTrimmingStrength,
        meshingVCS, vertexRefinement, outMesh, *uiAdapter.get_proxy() );
    boost::function<void( void )> f(
        boost::bind( frantic::volumetrics::implicitsurface::zhu_bridson_convert_sparse_particles_to_trimesh3,
                     boost::ref( params ) ) );
    uiAdapter.run( f );
}

void zhu_bridson_convert_sparse_particles_to_trimesh3(
    particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
    float effectRadius, float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
    const voxel_coord_system& meshingVCS, int vertexRefinement, trimesh3& outMesh,
    shared_progress_logger_proxy& progressLogger ) {
    particle_zhu_bridson_is_policy isp( particles, maxParticleRadius, effectRadius, lowDensityTrimmingDensity,
                                        lowDensityTrimmingStrength, meshingVCS, vertexRefinement );
    convert_particle_implicit_surface_to_trimesh3( isp, cpp, outMesh, progressLogger );
}

zhu_bridson_sparse_params::zhu_bridson_sparse_params(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maximumParticleRadius, float particleRadiusToEffectRadiusScale, float lowDensityTrimmingDensity,
    float lowDensityTrimmingStrength, const voxel_coord_system& meshingVCS, int vertexRefinement,
    frantic::geometry::trimesh3& outMesh, shared_progress_logger_proxy& progressLogger )
    : particles( particles )
    , cpp( cpp )
    , maximumParticleRadius( maximumParticleRadius )
    , particleRadiusToEffectRadiusScale( particleRadiusToEffectRadiusScale )
    , lowDensityTrimmingDensity( lowDensityTrimmingDensity )
    , lowDensityTrimmingStrength( lowDensityTrimmingStrength )
    , meshingVCS( meshingVCS )
    , vertexRefinement( vertexRefinement )
    , outMesh( outMesh )
    , progressLogger( progressLogger ) {}

void zhu_bridson_convert_sparse_particles_to_trimesh3( zhu_bridson_sparse_params& params ) {
    zhu_bridson_convert_sparse_particles_to_trimesh3(
        params.particles, params.cpp, params.maximumParticleRadius, params.particleRadiusToEffectRadiusScale,
        params.lowDensityTrimmingDensity, params.lowDensityTrimmingStrength, params.meshingVCS, params.vertexRefinement,
        params.outMesh, params.progressLogger );
}

/*
void zhu_bridson_convert_particles_to_trimesh3_sparse(
  frantic::particles::particle_grid_tree& particles,
  frantic::channels::channel_propagation_policy cpp,
  float maxParticleRadius,
  float effectRadius,
  float lowDensityTrimmingDensity,
  float lowDensityTrimmingStrength,
  const voxel_coord_system& meshingVCS,
  int vertexRefinement,
  frantic::geometry::trimesh3& outMesh)
{
  particle_zhu_bridson_is_policy isp( particles, maxParticleRadius, effectRadius, lowDensityTrimmingDensity,
lowDensityTrimmingStrength, meshingVCS, vertexRefinement ); convert_sparse_implicit_surface_to_trimesh3( cpp, isp,
outMesh );
}
*/

void anisotropic_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                                frantic::logging::progress_logger& progressLogger ) {
    particle_anisotropic_is_policy isp( particles, implicitThreshold, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh, progressLogger );
}

void anisotropic_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                int vertexRefinement, frantic::geometry::trimesh3& outMesh ) {
    particle_anisotropic_is_policy isp( particles, implicitThreshold, meshingVCS, vertexRefinement );
    convert_implicit_surface_to_trimesh3( cpp, isp, outMesh );
}

void anisotropic_convert_sparse_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                       const frantic::channels::channel_propagation_policy& cpp,
                                                       float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                       int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                                       frantic::logging::progress_logger& progressLogger ) {
    frantic::volumetrics::implicitsurface::shared_progress_logger_adapter uiAdapter( progressLogger );
    frantic::volumetrics::implicitsurface::anisotropic_sparse_params params(
        particles, cpp, implicitThreshold, meshingVCS, vertexRefinement, outMesh, *uiAdapter.get_proxy() );
    boost::function<void( void )> f(
        boost::bind( frantic::volumetrics::implicitsurface::anisotropic_convert_sparse_particles_to_trimesh3,
                     boost::ref( params ) ) );
    uiAdapter.run( f );
}

anisotropic_sparse_params::anisotropic_sparse_params( frantic::particles::particle_grid_tree& particles,
                                                      const frantic::channels::channel_propagation_policy& cpp,
                                                      float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                      int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                                      shared_progress_logger_proxy& progressLogger )
    : particles( particles )
    , cpp( cpp )
    , implicitThreshold( implicitThreshold )
    , meshingVCS( meshingVCS )
    , vertexRefinement( vertexRefinement )
    , outMesh( outMesh )
    , progressLogger( progressLogger ) {}

void anisotropic_convert_sparse_particles_to_trimesh3( anisotropic_sparse_params& params ) {
    particle_anisotropic_is_policy isp( params.particles, params.implicitThreshold, params.meshingVCS,
                                        params.vertexRefinement );
    convert_particle_implicit_surface_to_trimesh3( isp, params.cpp, params.outMesh, params.progressLogger );
}

void convert_levelset_to_trimesh3( const rle_level_set& levelSet, frantic::geometry::trimesh3& outMesh,
                                   int clipBounds ) {
    frantic::channels::channel_propagation_policy cpp;
    convert_levelset_to_trimesh3( levelSet, cpp, outMesh, clipBounds );
}

void convert_levelset_to_trimesh3( const rle_level_set& levelSet, frantic::channels::channel_propagation_policy& cpp,
                                   frantic::geometry::trimesh3& outMesh, int clipBounds ) {
    frantic::logging::null_progress_logger nullLogger;
    convert_levelset_to_trimesh3( levelSet, cpp, outMesh, nullLogger, clipBounds );
}

void convert_levelset_to_trimesh3( const rle_level_set& levelSet, frantic::channels::channel_propagation_policy& cpp,
                                   frantic::geometry::trimesh3& outMesh,
                                   frantic::logging::progress_logger& /*progressLogger*/, int clipBounds ) {
    // Create the policy, and call the implicit surface conversion function
    direct_linear_rle_level_set_is_policy mcp( levelSet, clipBounds );
    frantic::volumetrics::implicitsurface::convert_dense_implicit_surface_to_trimesh3( cpp, mcp,
                                                                                       outMesh /*, progressLogger*/ );
}

/*
void convert_levelset_to_trimesh3_sparse( const rle_level_set& levelSet, frantic::geometry::trimesh3& outMesh, int
clipBounds )
{
  // Create the policy, and call the implicit surface conversion function
  direct_linear_rle_level_set_is_policy isp( levelSet, clipBounds );
  frantic::channels::channel_propagation_policy cpp;
  frantic::volumetrics::implicitsurface::convert_sparse_implicit_surface_to_trimesh3( cpp, isp, outMesh );
}
*/

/*
void convert_levelset_to_trimesh3_sparse( const rle_level_set& levelSet, frantic::channels::channel_propagation_policy&
cpp, frantic::geometry::trimesh3& outMesh, int clipBounds )
{
  // Create the policy, and call the implicit surface conversion function
  direct_linear_rle_level_set_is_policy isp( levelSet, clipBounds );
  frantic::volumetrics::implicitsurface::convert_sparse_implicit_surface_to_trimesh3( cpp, isp, outMesh );
}
*/

// TODO: Also specify which reconstruction filter to use
void convert_levelset_to_trimesh3( const rle_level_set& levelSet, const voxel_coord_system& meshingVCS,
                                   int vertexRefinement, frantic::geometry::trimesh3& outMesh, int clipBounds ) {
    frantic::channels::channel_propagation_policy cpp;
    convert_levelset_to_trimesh3( levelSet, cpp, meshingVCS, vertexRefinement, outMesh, clipBounds );
}

void convert_levelset_to_trimesh3( const rle_level_set& levelSet, frantic::channels::channel_propagation_policy& cpp,
                                   const voxel_coord_system& meshingVCS, int vertexRefinement,
                                   frantic::geometry::trimesh3& outMesh, int clipBounds ) {
    frantic::logging::null_progress_logger nullLogger;
    convert_levelset_to_trimesh3( levelSet, cpp, meshingVCS, vertexRefinement, outMesh, nullLogger, clipBounds );
}

void convert_levelset_to_trimesh3( const rle_level_set& levelSet, frantic::channels::channel_propagation_policy& cpp,
                                   const voxel_coord_system& meshingVCS, int vertexRefinement,
                                   frantic::geometry::trimesh3& outMesh,
                                   frantic::logging::progress_logger& progressLogger, int clipBounds ) {
    // Create the policy, and call the implicit surface conversion function
    reconstruction_filtered_rle_level_set_is_policy<math::b_spline_filter> mcp( levelSet, meshingVCS, vertexRefinement,
                                                                                math::b_spline_filter(), clipBounds );
    frantic::volumetrics::implicitsurface::convert_implicit_surface_to_trimesh3( cpp, mcp, outMesh, progressLogger );
}

void fix_marching_cubes_topology( frantic::geometry::trimesh3& mesh,
                                  frantic::logging::progress_logger& progressLogger ) {
    //  Removes faces in which all vertices are part of a non-manifold
    // edge.
    //
    //  This is intended to remove a specific non-manifold feature that
    // we have observed in the output mesh.  The feature is a double-sided
    // quad connecting two surfaces together -- imagine a flattened
    // drinking straw connecting two inflated balloons.
    frantic::geometry::remove_duplicate_halfedge_faces( mesh, progressLogger, true );
}

/*
void convert_levelset_to_trimesh3_sparse(
  const rle_level_set& levelSet,
  frantic::channels::channel_propagation_policy& cpp,
  const voxel_coord_system& meshingVCS,
  int vertexRefinement,
  frantic::geometry::trimesh3& outMesh,
  int clipBounds
 )
{
  // Create the policy, and call the implicit surface conversion function
  reconstruction_filtered_rle_level_set_is_policy<math::b_spline_filter> mcp( levelSet, meshingVCS, vertexRefinement,
math::b_spline_filter(), clipBounds );
  frantic::volumetrics::implicitsurface::convert_sparse_implicit_surface_to_trimesh3( cpp, mcp, outMesh );
}
*/

namespace detail {

// Given a slice of density samples, this function fills an array of cube case values for initial setup.  It fills the
// cube cases as if the voxel corner densities were duplicated.
//
// NOTE: The size of the output cube case values array will equal the size of the voxel corner density arrays, because
// the first values along both x and y are generated based on taking the first corner density twice.  This is done this
// way to deal with open edges properly.
void fill_initial_cube_case_values( size2 size, const std::vector<float>& voxelCornerDensities,
                                    std::vector<unsigned char>& outCubeCases ) {
    if( voxelCornerDensities.empty() )
        return;

    int i = 0;

    unsigned char cubeCase = 0;

    // Special case the (0,0) corner, which is for 8 copies of the corner density (both in the first row and the first
    // column)
    outCubeCases[i++] = ( voxelCornerDensities[0] < 0 ) ? (unsigned char)0xff : (unsigned char)0x00;

    // Special case the first row
    for( int x = 1; x < size.xsize; ++x ) {
        cubeCase = 0;
        if( voxelCornerDensities[x - 1] < 0 )
            cubeCase |= 0x55;
        if( voxelCornerDensities[x] < 0 )
            cubeCase |= 0xaa;
        outCubeCases[i++] = cubeCase;
    }

    for( int y = 1; y < size.ysize; ++y ) {
        int yoffset = y * size.xsize;

        // Special case the first column
        cubeCase = 0;
        if( voxelCornerDensities[yoffset - size.xsize] < 0 )
            cubeCase |= 0x33;
        if( voxelCornerDensities[yoffset] < 0 )
            cubeCase |= 0xcc;
        outCubeCases[i++] = cubeCase;

        // Now do the general case for the rest of this row
        for( int x = 1; x < size.xsize; ++x ) {
            cubeCase = 0;
            if( voxelCornerDensities[yoffset + x - size.xsize - 1] < 0 )
                cubeCase |= 0x11;
            if( voxelCornerDensities[yoffset + x - size.xsize] < 0 )
                cubeCase |= 0x22;
            if( voxelCornerDensities[yoffset + x - 1] < 0 )
                cubeCase |= 0x44;
            if( voxelCornerDensities[yoffset + x] < 0 )
                cubeCase |= 0x88;
            outCubeCases[i++] = cubeCase;
        }
    }

    //	cout << "fill_initial_cube_case_values( " << size <<  " );" << endl;
    //	cout << "  results: ";
    //	for( unsigned i = 0; i < outCubeCases.size(); ++i )
    //		cout << strings::to_hex_string(outCubeCases[i]) << " ";
    //	cout << endl;
}

void fill_initial_cube_case_values( frantic::graphics2d::size2 size, const std::vector<float>& voxelCornerDensities,
                                    std::vector<unsigned char>& outCubeCases,
                                    std::vector<boost::int32_t>& outNewVertexCount ) {
    if( voxelCornerDensities.empty() )
        return;

    if( (int)outNewVertexCount.size() != size.ysize ) {
        throw std::runtime_error(
            "fill_initial_cube_case_values Error: vertex count array size does not match ysize (" +
            boost::lexical_cast<std::string>( outNewVertexCount.size() ) +
            " != " + boost::lexical_cast<std::string>( size.ysize ) + ")" );
    }

    int i = 0;

    unsigned char cubeCase = 0;

    boost::int32_t newVertexCount = 0;

    // Special case the (0,0) corner, which is for 8 copies of the corner density (both in the first row and the first
    // column)
    cubeCase = ( voxelCornerDensities[0] < 0 ) ? (unsigned char)0xff : (unsigned char)0x00;
    outCubeCases[i++] = cubeCase;

    // Special case the first row
    for( int x = 1; x < size.xsize; ++x ) {
        cubeCase = 0;
        if( voxelCornerDensities[x - 1] < 0 )
            cubeCase |= 0x55;
        if( voxelCornerDensities[x] < 0 )
            cubeCase |= 0xaa;
        outCubeCases[i++] = cubeCase;
        if( cubeCase != 0 && cubeCase != 0xff ) {
            newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
        }
    }
    outNewVertexCount[0] = newVertexCount;

    for( int y = 1; y < size.ysize; ++y ) {
        int yoffset = y * size.xsize;
        newVertexCount = 0;

        // Special case the first column
        cubeCase = 0;
        if( voxelCornerDensities[yoffset - size.xsize] < 0 )
            cubeCase |= 0x33;
        if( voxelCornerDensities[yoffset] < 0 )
            cubeCase |= 0xcc;
        outCubeCases[i++] = cubeCase;
        if( cubeCase != 0 && cubeCase != 0xff ) {
            newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
        }

        // Now do the general case for the rest of this row
        for( int x = 1; x < size.xsize; ++x ) {
            cubeCase = 0;
            if( voxelCornerDensities[yoffset + x - size.xsize - 1] < 0 )
                cubeCase |= 0x11;
            if( voxelCornerDensities[yoffset + x - size.xsize] < 0 )
                cubeCase |= 0x22;
            if( voxelCornerDensities[yoffset + x - 1] < 0 )
                cubeCase |= 0x44;
            if( voxelCornerDensities[yoffset + x] < 0 )
                cubeCase |= 0x88;
            outCubeCases[i++] = cubeCase;
            if( cubeCase != 0 && cubeCase != 0xff ) {
                newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
            }
        }
        outNewVertexCount[y] = newVertexCount;
    }
    detail::cumulative_sum( outNewVertexCount );
}

/**
 *	This function just populates the vert indicies in the plane where there are adjacent defined verts,
 *	because a vert may need to be generated there later, and it will be necessary to check against an
 *	initial value to see if one was already generated.
 *
 *	@param	voxelCornerDensities	The voxel corner samples for the plane.
 *	@param	currentRLP				The rle plane for traversing the samples.
 *	@param	outVertexIndices		The vertex index array to be flagged in the appropriate places.
 */
void fill_sparse_initial_vertex_indices( const boost::shared_array<float>& /*voxelCornerDensities*/,
                                         const rle_plane& currentRLP, boost::shared_array<vector3>& outVertexIndices ) {
    if( currentRLP.is_empty() )
        return;

    const int xsize = currentRLP.get_xy_extents().xsize();
    const int ysize = currentRLP.get_xy_extents().ysize();

    rle_plane_defined_iterator prevRowIter( &currentRLP );

    // On first extent, only adjacent x defined corners after the starts of runs
    int run = currentRLP.get_y_extent_first_run( 0 );
    if( run >= 0 ) {
        for( ; run <= currentRLP.get_y_extent_last_run( 0 ); ++run ) {
            for( int x = currentRLP[run].first + 1; x <= currentRLP[run].second; ++x ) {
                outVertexIndices[x].x = -1;
            }
        }
    }

    // Go through each defined corner sample indicated in the rle plane
    for( int y = 1; y < ysize; ++y ) {
        int run = currentRLP.get_y_extent_first_run( y );

        // Skip any empty extents
        if( run < 0 )
            continue;

        // never adjacent x defined corners on first of run, but maybe y
        int x = currentRLP[run].first;
        prevRowIter.jump_to( x - xsize );
        if( prevRowIter.get_data_index() == x - xsize )
            outVertexIndices[x].y = -1;

        // any internal run location will have a defined pair of x corners, but have to check for y
        for( ; run <= currentRLP.get_y_extent_last_run( y ); ++run ) {
            for( x++; x <= currentRLP[run].second; ++x ) {
                outVertexIndices[x].x = -1;
                prevRowIter.jump_to( x - xsize );
                if( prevRowIter.get_data_index() == x - xsize )
                    outVertexIndices[x].y = -1;
            }
        }
    }
}

void fill_cube_case_values( frantic::graphics2d::size2 size, const std::vector<unsigned char>& previousCubeCases,
                            const std::vector<float>& voxelCornerDensities, std::vector<unsigned char>& outCubeCases ) {
    int i = 0;
    unsigned char cubeCase = 0;

    // Special case the (0,0) corner
    outCubeCases[i++] = ( ( previousCubeCases[0] & 0xf0 ) >> 4 ) |
                        ( ( voxelCornerDensities[0] < 0 ) ? (unsigned char)0xf0 : (unsigned char)0x00 );

    // Special case the first row
    for( int x = 1; x < size.xsize; ++x ) {
        cubeCase = ( previousCubeCases[i] & 0xf0 ) >> 4;
        if( voxelCornerDensities[x - 1] < 0 )
            cubeCase |= 0x50;
        if( voxelCornerDensities[x] < 0 )
            cubeCase |= 0xa0;
        outCubeCases[i++] = cubeCase;
    }

    for( int y = 1; y < size.ysize; ++y ) {
        int yoffset = y * size.xsize;

        // Special case the first column
        cubeCase = ( previousCubeCases[i] & 0xf0 ) >> 4;
        if( voxelCornerDensities[yoffset - size.xsize] < 0 )
            cubeCase |= 0x30;
        if( voxelCornerDensities[yoffset] < 0 )
            cubeCase |= 0xc0;
        outCubeCases[i++] = cubeCase;

        // Now do the general case for the rest of this row
        for( int x = 1; x < size.xsize; ++x ) {
            outCubeCases[i] = (unsigned char)marching_cubes_table::combine_cube_cases(
                outCubeCases[i - 1], outCubeCases[i - size.xsize], previousCubeCases[i],
                voxelCornerDensities[yoffset + x] < 0 );
            ++i;
        }
    }
    //	cout << "fill_cube_case_values( " << size <<  " );" << endl;
    //	cout << "  results: ";
    //	for( unsigned i = 0; i < outCubeCases.size(); ++i )
    //		cout << strings::to_hex_string(outCubeCases[i]) << " ";
    //	cout << endl;
}

void fill_cube_case_values( const frantic::graphics::size3& size, const std::vector<float>& voxelCornerDensities_,
                            std::vector<boost::uint8_t>& outCubeCases_ ) {
    using boost::uint8_t;
    using std::size_t;

    const float* voxelCornerDensities = voxelCornerDensities_.size() > 0 ? &voxelCornerDensities_[0] : 0;
    boost::uint8_t* outCubeCases = outCubeCases_.size() > 0 ? &outCubeCases_[0] : 0;

    const size_t xsize = size.xsize();
    const size_t ysize = size.ysize();
    const size_t zsize = size.zsize();
    const size_t xyArea = size.xsize() * size.ysize();

    size_t i = 0;
    uint8_t cubeCase = 0;

    // Special case the (0,0) corner, which is for 8 copies of the corner density (both in the first row and the first
    // column)
    outCubeCases[i++] = ( voxelCornerDensities[0] < 0 ) ? (uint8_t)0xff : (uint8_t)0x00;

    // Special case the first row
    for( size_t x = 1; x < xsize; ++x ) {
        cubeCase = 0;
        if( voxelCornerDensities[i - 1] < 0 )
            cubeCase |= 0x55;
        if( voxelCornerDensities[i] < 0 )
            cubeCase |= 0xaa;
        outCubeCases[i++] = cubeCase;
    }

    for( size_t y = 1; y < ysize; ++y ) {
        cubeCase = 0;
        if( voxelCornerDensities[i - xsize] < 0 )
            cubeCase |= 0x33;
        if( voxelCornerDensities[i] < 0 )
            cubeCase |= 0xcc;
        outCubeCases[i++] = cubeCase;

        for( size_t x = 1; x < xsize; ++x ) {
            cubeCase = 0;
            if( voxelCornerDensities[i - xsize - 1] < 0 )
                cubeCase |= 0x11;
            if( voxelCornerDensities[i - xsize] < 0 )
                cubeCase |= 0x22;
            if( voxelCornerDensities[i - 1] < 0 )
                cubeCase |= 0x44;
            if( voxelCornerDensities[i] < 0 )
                cubeCase |= 0x88;
            outCubeCases[i++] = cubeCase;
        }
    }

    for( size_t z = 1; z < zsize; ++z ) {
        // Special case the (0,0) corner
        outCubeCases[i] = ( ( outCubeCases[i - xyArea] & 0xf0 ) >> 4 ) |
                          ( ( voxelCornerDensities[i] < 0 ) ? (uint8_t)0xf0 : (uint8_t)0x00 );
        ++i;

        // Special case the first row
        for( size_t x = 1; x < xsize; ++x ) {
            cubeCase = ( outCubeCases[i - xyArea] & 0xf0 ) >> 4;
            if( voxelCornerDensities[i - 1] < 0 )
                cubeCase |= 0x50;
            if( voxelCornerDensities[i] < 0 )
                cubeCase |= 0xa0;
            outCubeCases[i++] = cubeCase;
        }

        for( size_t y = 1; y < ysize; ++y ) {
            // Special case the first column
            cubeCase = ( outCubeCases[i - xyArea] & 0xf0 ) >> 4;
            if( voxelCornerDensities[i - xsize] < 0 )
                cubeCase |= 0x30;
            if( voxelCornerDensities[i] < 0 )
                cubeCase |= 0xc0;
            outCubeCases[i++] = cubeCase;

            // Now do the general case for the rest of this row
            for( size_t x = 1; x < xsize; ++x ) {
                outCubeCases[i] = (unsigned char)marching_cubes_table::combine_cube_cases(
                    outCubeCases[i - 1], outCubeCases[i - xsize], outCubeCases[i - xyArea],
                    voxelCornerDensities[i] < 0 );
                // outCubeCases[i) = get_cube_case( outCubeCases[i-1-xsize-xyArea), outCubeCases[i-xsize-xyArea),
                // outCubeCases[i-1-xyArea), outCubeCases[i-xyArea), outCubeCases[i-1-xsize), outCubeCases[i-xsize),
                // outCubeCases[i-1), voxelCornerDensities[i) < 0 ); std::cout << int( outCubeCases[i) ) << " " << int(
                // get_cube_case( voxelCornerDensities[i-1-xsize-xyArea), voxelCornerDensities[i-xsize-xyArea),
                // voxelCornerDensities[i-1-xyArea), voxelCornerDensities[i-xyArea), voxelCornerDensities[i-1-xsize),
                // voxelCornerDensities[i-xsize), voxelCornerDensities[i-1), voxelCornerDensities[i) < 0 ) ) <<
                // std::endl; outCubeCases[i) = get_cube_case( voxelCornerDensities[i-1-xsize-xyArea),
                // voxelCornerDensities[i-xsize-xyArea), voxelCornerDensities[i-1-xyArea),
                // voxelCornerDensities[i-xyArea), voxelCornerDensities[i-1-xsize), voxelCornerDensities[i-xsize),
                // voxelCornerDensities[i-1), voxelCornerDensities[i) );
                ++i;
            }
        }
    }
}

struct FillCubeCaseValuesBody {
    const frantic::graphics2d::size2& m_size;
    const std::vector<boost::uint8_t>& m_previousCubeCases;
    const std::vector<float>& m_voxelCornerDensities;
    std::vector<boost::uint8_t>& m_outCubeCases;

#pragma warning( push )
#pragma warning( disable : 4822 ) // local class member function does not have a body
    FillCubeCaseValuesBody& operator=( const FillCubeCaseValuesBody& ); // not implemented
#pragma warning( pop )

    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        const std::size_t yStart = r.begin();
        const std::size_t yEnd = r.end();

        unsigned char cubeCase = 0;

        std::size_t y = yStart;

        // avoid secure scl checks
        // TODO: get rid of this when you're no longer using msvc8/9
        const float* voxelCornerDensities = m_voxelCornerDensities.size() > 0 ? &m_voxelCornerDensities[0] : 0;
        const boost::uint8_t* previousCubeCases = m_previousCubeCases.size() > 0 ? &m_previousCubeCases[0] : 0;
        boost::uint8_t* outCubeCases = m_outCubeCases.size() > 0 ? &m_outCubeCases[0] : 0;

        if( y == 0 ) {
            // Special case the (0,0) corner
            cubeCase = ( ( previousCubeCases[0] & 0xf0 ) >> 4 ) |
                       ( ( voxelCornerDensities[0] < 0 ) ? (unsigned char)0xf0 : (unsigned char)0x00 );
            outCubeCases[0] = cubeCase;

            // Special case the first row
            for( int x = 1; x < m_size.xsize && x < m_outCubeCases.size(); ++x ) {
                cubeCase = ( previousCubeCases[x] & 0xf0 ) >> 4;
                if( voxelCornerDensities[x - 1] < 0 )
                    cubeCase |= 0x50;
                if( voxelCornerDensities[x] < 0 )
                    cubeCase |= 0xa0;
                outCubeCases[x] = cubeCase;
            }

            ++y;
        }

        std::size_t i = y * m_size.xsize;
        for( ; y < yEnd; ++y ) {

            // Special case the first column
            cubeCase = ( previousCubeCases[i] & 0xf0 ) >> 4;
            if( voxelCornerDensities[i - m_size.xsize] < 0 )
                cubeCase |= 0x30;
            if( voxelCornerDensities[i] < 0 )
                cubeCase |= 0xc0;
            outCubeCases[i++] = cubeCase;

            // Now do the general case for the rest of this row
            if( y == yStart ) {
                for( int x = 1; x < m_size.xsize; ++x ) {
                    cubeCase = (unsigned char)marching_cubes_table::combine_cube_cases(
                        outCubeCases[i - 1], voxelCornerDensities[i - m_size.xsize] < 0, previousCubeCases[i],
                        voxelCornerDensities[i] < 0 );
                    outCubeCases[i] = cubeCase;
                    ++i;
                }
            } else {
                for( int x = 1; x < m_size.xsize; ++x ) {
                    cubeCase = (unsigned char)marching_cubes_table::combine_cube_cases(
                        outCubeCases[i - 1], outCubeCases[i - m_size.xsize], previousCubeCases[i],
                        voxelCornerDensities[i] < 0 );
                    outCubeCases[i] = cubeCase;
                    ++i;
                }
            }
        }
    }

    FillCubeCaseValuesBody( const frantic::graphics2d::size2& size, const std::vector<unsigned char>& previousCubeCases,
                            const std::vector<float>& voxelCornerDensities, std::vector<unsigned char>& outCubeCases )
        : m_size( size )
        , m_previousCubeCases( previousCubeCases )
        , m_voxelCornerDensities( voxelCornerDensities )
        , m_outCubeCases( outCubeCases ) {}
};

void fill_cube_case_values_mt( tbb::affinity_partitioner& partitioner, frantic::graphics2d::size2 size,
                               const std::vector<unsigned char>& previousCubeCases,
                               const std::vector<float>& voxelCornerDensities,
                               std::vector<unsigned char>& outCubeCases ) {
    FillCubeCaseValuesBody body( size, previousCubeCases, voxelCornerDensities, /*outAddedVertexCount,*/ outCubeCases );
    tbb::blocked_range<std::size_t> range( 0, size.ysize );
#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, partitioner );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif
}

struct FillCubeCaseValuesBody2 {
    const frantic::graphics2d::size2& m_size;
    const std::vector<boost::uint8_t>& m_previousCubeCases;
    const std::vector<float>& m_voxelCornerDensities;
    std::vector<boost::uint8_t>& m_outCubeCases;
    std::vector<boost::int32_t>& m_outNewVertexCount;

#pragma warning( push )
#pragma warning( disable : 4822 ) // local class member function does not have a body
    FillCubeCaseValuesBody2& operator=( const FillCubeCaseValuesBody2& ); // not implemented
#pragma warning( pop )

    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        const std::size_t yStart = r.begin();
        const std::size_t yEnd = r.end();

        unsigned char cubeCase = 0;

        boost::int32_t newVertexCount = 0;

        std::size_t y = yStart;

        // avoid secure scl checks
        // TODO: get rid of this when you're no longer using msvc8/9
        const float* voxelCornerDensities = m_voxelCornerDensities.size() > 0 ? &m_voxelCornerDensities[0] : 0;
        const boost::uint8_t* previousCubeCases = m_previousCubeCases.size() > 0 ? &m_previousCubeCases[0] : 0;
        boost::uint8_t* outCubeCases = m_outCubeCases.size() > 0 ? &m_outCubeCases[0] : 0;
        boost::int32_t* outNewVertexCount = m_outNewVertexCount.size() > 0 ? &m_outNewVertexCount[0] : 0;

        if( y == 0 ) {
            // Special case the (0,0) corner
            cubeCase = ( ( previousCubeCases[0] & 0xf0 ) >> 4 ) |
                       ( ( voxelCornerDensities[0] < 0 ) ? (unsigned char)0xf0 : (unsigned char)0x00 );
            outCubeCases[0] = cubeCase;
            if( cubeCase != 0 && cubeCase != 0xff ) {
                newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
            }

            // Special case the first row
            for( int x = 1; x < m_size.xsize && x < m_outCubeCases.size(); ++x ) {
                cubeCase = ( previousCubeCases[x] & 0xf0 ) >> 4;
                if( voxelCornerDensities[x - 1] < 0 )
                    cubeCase |= 0x50;
                if( voxelCornerDensities[x] < 0 )
                    cubeCase |= 0xa0;
                outCubeCases[x] = cubeCase;
                if( cubeCase != 0 && cubeCase != 0xff ) {
                    newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
                }
            }
            outNewVertexCount[y] = newVertexCount;

            ++y;
        }

        std::size_t i = y * m_size.xsize;
        for( ; y < yEnd; ++y ) {
            newVertexCount = 0;

            // Special case the first column
            cubeCase = ( previousCubeCases[i] & 0xf0 ) >> 4;
            if( voxelCornerDensities[i - m_size.xsize] < 0 )
                cubeCase |= 0x30;
            if( voxelCornerDensities[i] < 0 )
                cubeCase |= 0xc0;
            outCubeCases[i++] = cubeCase;
            if( cubeCase != 0 && cubeCase != 0xff ) {
                newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
            }

            // Now do the general case for the rest of this row
            if( y == yStart ) {
                for( int x = 1; x < m_size.xsize; ++x ) {
                    cubeCase = (unsigned char)marching_cubes_table::combine_cube_cases(
                        outCubeCases[i - 1], voxelCornerDensities[i - m_size.xsize] < 0, previousCubeCases[i],
                        voxelCornerDensities[i] < 0 );
                    outCubeCases[i] = cubeCase;
                    if( cubeCase != 0 && cubeCase != 0xff ) {
                        newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
                    }
                    ++i;
                }
            } else {
                for( int x = 1; x < m_size.xsize; ++x ) {
                    cubeCase = (unsigned char)marching_cubes_table::combine_cube_cases(
                        outCubeCases[i - 1], outCubeCases[i - m_size.xsize], previousCubeCases[i],
                        voxelCornerDensities[i] < 0 );
                    outCubeCases[i] = cubeCase;
                    if( cubeCase != 0 && cubeCase != 0xff ) {
                        newVertexCount += frantic::volumetrics::marching_cubes_table::get_added_vert_count( cubeCase );
                    }
                    ++i;
                }
            }
            outNewVertexCount[y] = newVertexCount;
        }
    }

    FillCubeCaseValuesBody2( const frantic::graphics2d::size2& size,
                             const std::vector<unsigned char>& previousCubeCases,
                             const std::vector<float>& voxelCornerDensities, std::vector<unsigned char>& outCubeCases,
                             std::vector<boost::int32_t>& outNewVertexCount )
        : m_size( size )
        , m_previousCubeCases( previousCubeCases )
        , m_voxelCornerDensities( voxelCornerDensities )
        , m_outCubeCases( outCubeCases )
        , m_outNewVertexCount( outNewVertexCount ) {}
};

void fill_cube_case_values_mt( tbb::affinity_partitioner& partitioner, frantic::graphics2d::size2 size,
                               const std::vector<unsigned char>& previousCubeCases,
                               const std::vector<float>& voxelCornerDensities, std::vector<unsigned char>& outCubeCases,
                               std::vector<boost::int32_t>& outNewVertexCount ) {
    if( (int)outNewVertexCount.size() != size.ysize ) {
        throw std::runtime_error( "fill_cube_case_values_mt Error: vertex count array size does not match ysize (" +
                                  boost::lexical_cast<std::string>( outNewVertexCount.size() ) +
                                  " != " + boost::lexical_cast<std::string>( size.ysize ) + ")" );
    }

    FillCubeCaseValuesBody2 body( size, previousCubeCases, voxelCornerDensities, outCubeCases, outNewVertexCount );
    tbb::blocked_range<std::size_t> range( 0, size.ysize );
#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, partitioner );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif
    detail::cumulative_sum( outNewVertexCount );
}

/**
 *	This function generates a new set of cube cases, using the previous cube cases
 *	and the new voxel corner density samples.  It does so sparsely using rle plane structures
 *	for each plane of voxelCornerDensities.  This means that, unfortunately, we can't blindly use
 *	previous data as the dense version did.  We have to check the previous corners that belong to each
 *  new cube case being examined for their defined/undefined status, as well as their inside/ouside value.
 *	This requires the maintenance of iterators for the previous plane (current and previous rows) and the
 *	current plane (previous row).
 *
 *	We treat each sparse corner sample as the upper back corner of a cube case (the 7th vert).  As a result,
 *	there are the same number of cube cases as there are corners in the plane, and the rle plane structure
 *	applies to both the voxel corner samples and the cube cases.  However, not all of the cube
 *	cases will be fully defined.  A new rlp structure is created to reflect fully defined cube cases which
 *  will eventually need to be processed by the later stages of the algorithm.
 *
 *	It also initializes the outVertexIndices array for the current plane of samples.  That is, it puts a -1
 *	flag value in for any pair of defined samples that it finds on the leading edges of the associated cube
 *	case.  This array is then used later on to determine if a vert has been generated for that edge when
 *	processing the sparse fully defined cube cases.  Any candidate location to be checked will have been
 *	initialized with a -1 value.
 *
 *	@param	prevVoxelCornerDensities	Corner samples from the previous plane.
 *	@param	prevRLP						rle plane for traversing the samples from the previous
 *plane
 *	@param	currentVoxelCornerDensities	Corner samples from the current plane.
 *	@param	currentRLP					rle plane for traversing the samples from the current
 *plane
 *	@param	outDefinedCubeCases			the fully defined cube cases for this plane
 *	@param	outDefinedCubeCasesRLP		the rlp structure for the defined cube cases
 *	@param	outVertexIndices			the vertex index flags/values for this plane
 */
void fill_sparse_cube_case_values( const float /*defaultOutsideDensity*/,
                                   const boost::shared_array<float>& prevVoxelCornerDensities, const rle_plane& prevRLP,
                                   const boost::shared_array<float>& voxelCornerDensities, const rle_plane& currentRLP,
                                   boost::shared_array<unsigned char>& outDefinedCubeCases,
                                   rle_plane& outDefinedCubeCasesRLP, boost::shared_array<vector3>& outVertexIndices ) {

    outDefinedCubeCasesRLP.reset( currentRLP.get_xy_extents() );

    if( prevRLP.is_empty() && currentRLP.is_empty() )
        return;

    const int xsize = currentRLP.get_xy_extents().xsize();
    const int ysize = currentRLP.get_xy_extents().ysize();

    // Create iterators for the previous plane and previous row
    rle_plane_defined_iterator prevPlanePrevRowIter( &prevRLP ), prevPlaneCurrentRowIter( &prevRLP ),
        prevRowIter( &currentRLP );

    unsigned char cubeCase = 0, definedCornerCount = 0; // cubeCaseFlags = 0;

    // There won't be any defined cube cases for the first extent, being that they are all on the edge.
    std::vector<std::pair<int, int>> runs;
    outDefinedCubeCasesRLP.append_runs_by_extent( runs, 0 );

    // However, we still have to check for defined corners so we can initialize the vertex indices array
    // On the first extent, only adjacent x defined corners are possible, and only after the start of a run
    int run = currentRLP.get_y_extent_first_run( 0 );
    if( run >= 0 ) {
        for( ; run <= currentRLP.get_y_extent_last_run( 0 ); ++run ) {
            for( int x = currentRLP[run].first + 1; x <= currentRLP[run].second; ++x ) {
                outVertexIndices[x].x = -1;
            }
        }
    }

    // Go through each defined corner sample indicated in the rle plane and define cube cases and
    // flags where appropriate.
    for( int y = 1; y < ysize; ++y ) {
        runs.clear();
        runs.reserve( xsize / 2 ); // upper bound on number of runs, to avoid on the fly resizing
        int run = currentRLP.get_y_extent_first_run( y );

        // Skip any empty extents
        if( run < 0 ) {
            outDefinedCubeCasesRLP.append_runs_by_extent( runs, y );
            continue;
        }

        for( ; run <= currentRLP.get_y_extent_last_run( y ); ++run ) {

            int x = currentRLP[run].first;

            // Check for defined vert pairs.  Won't be needing an x vert, since it's isolated
            // in that direction, might need a y,z
            prevRowIter.jump_to( x - xsize );
            if( prevRowIter.get_data_index() == x - xsize )
                outVertexIndices[x].y = -1;
            prevPlaneCurrentRowIter.jump_to( x );
            if( prevPlaneCurrentRowIter.get_data_index() == x )
                outVertexIndices[x].z = -1;

            // If the run is only one long, skip it, no defined cube case
            if( currentRLP[run].first == currentRLP[run].second ) {
                continue;
            }

            for( x++; x <= currentRLP[run].second; ++x ) {

                cubeCase = 0;
                definedCornerCount = 0;

                // we might need an x vert here at some point, so initialize the vert index
                outVertexIndices[x].x = -1;

                // Set up the cube corners that lie in the previous plane

                int prevCornerIndex = x - xsize - 1;
                prevPlanePrevRowIter.jump_to( prevCornerIndex );
                // Check if the previous column/row corner was defined
                if( prevPlanePrevRowIter.get_data_index() == prevCornerIndex ) {
                    if( prevVoxelCornerDensities[prevCornerIndex] < 0 )
                        cubeCase |= 0x01;
                    definedCornerCount++;
                }

                ++prevCornerIndex;
                prevPlanePrevRowIter.jump_to( prevCornerIndex );

                // Check if the previous row corner was defined
                if( prevPlanePrevRowIter.get_data_index() == prevCornerIndex ) {
                    if( prevVoxelCornerDensities[prevCornerIndex] < 0 )
                        cubeCase |= 0x02;
                    definedCornerCount++;
                }

                prevCornerIndex = x - 1;
                prevPlaneCurrentRowIter.jump_to( prevCornerIndex );
                // Check if the previous column corner was defined
                if( prevPlaneCurrentRowIter.get_data_index() == prevCornerIndex ) {
                    if( prevVoxelCornerDensities[prevCornerIndex] < 0 )
                        cubeCase |= 0x04;
                    definedCornerCount++;
                }

                // Check if the previous current corner was defined
                ++prevCornerIndex;
                prevPlaneCurrentRowIter.jump_to( prevCornerIndex );
                if( prevPlaneCurrentRowIter.get_data_index() == x ) {
                    if( prevVoxelCornerDensities[prevCornerIndex] < 0 )
                        cubeCase |= 0x08;
                    definedCornerCount++;

                    // might need a vert between these defined corners, so initialize the vert index
                    outVertexIndices[x].z = -1;
                }

                // Set up the current plane's cube corners

                // Move up the previous row iterator
                prevCornerIndex = x - xsize - 1;
                prevRowIter.jump_to( prevCornerIndex );
                // Check if the previous column/row corner was defined
                if( prevRowIter.get_data_index() == prevCornerIndex ) {
                    if( voxelCornerDensities[prevCornerIndex] < 0 )
                        cubeCase |= 0x10;
                    definedCornerCount++;
                }

                ++prevCornerIndex;
                prevRowIter.jump_to( prevCornerIndex );
                // Check if the previous row corner was defined
                if( prevRowIter.get_data_index() == prevCornerIndex ) {
                    if( voxelCornerDensities[prevCornerIndex] < 0 )
                        cubeCase |= 0x20;
                    definedCornerCount++;

                    // could need a y vert between these defined corners, initialize the vert index
                    outVertexIndices[x].y = -1;
                }

                prevCornerIndex = x - 1;
                // Check the previous x position
                if( voxelCornerDensities[prevCornerIndex] < 0 )
                    cubeCase |= 0x40;

                // Define the current x position
                if( voxelCornerDensities[x] < 0 )
                    cubeCase |= 0x80;

                // If the cube case is fully defined, add it to the fully defined cube case rle plane
                if( definedCornerCount == 6 ) {

                    // first run
                    if( runs.size() == 0 ) {
                        runs.push_back( std::pair<int, int>( x, x ) );
                    }

                    // If we moved more than one since the last defined position, need a new run
                    if( x - runs.back().second > 1 )
                        runs.push_back( std::pair<int, int>( x, x ) );

                    runs.back().second = x;

                    // add the defined cube case
                    outDefinedCubeCases[x] = cubeCase;
                }
            }
        }

        // append the runs of defined cube cases for this extent
        outDefinedCubeCasesRLP.append_runs_by_extent( runs, y );
    }
}

void generate_faces_for_plane( const marching_cubes_table& mct, size2 size,
                               const vector<unsigned char>& previousCubeCases,
                               const vector<unsigned char>& currentCubeCases, const vector<int>& previousVertexIndices,
                               const vector<int>& currentVertexIndices, trimesh3& outMesh ) {
    //	cout << "generate_faces_for_plane( " << size << ") " << endl;

    // Allocate a vector for the 12 cube vertices, and the faces
    vector<boost::int32_t> cubeVerts( 12 );

    // Skip the first row
    int index = size.xsize;
    for( int y = 1; y < size.ysize; ++y ) {
        // Skip the first column
        ++index;
        for( int x = 1; x < size.xsize; ++x ) {
            unsigned char cubeCase = currentCubeCases[index];
            if( cubeCase != 0x00 && cubeCase != 0xff ) {

                // Retrieve all the vertices from the edges of this cube
                mct.fill_verts( cubeCase, currentVertexIndices[index], currentCubeCases[index - 1],
                                currentCubeCases[index - size.xsize], previousCubeCases[index],
                                currentCubeCases[index - size.xsize - 1], previousCubeCases[index - 1],
                                previousCubeCases[index - size.xsize], currentVertexIndices[index - 1],
                                currentVertexIndices[index - size.xsize], previousVertexIndices[index],
                                currentVertexIndices[index - size.xsize - 1], previousVertexIndices[index - 1],
                                previousVertexIndices[index - size.xsize], cubeVerts );

                const vector<vector3>& cubeFaces = mct.get_cubecase_faces( cubeCase );

                for( unsigned i = 0; i < cubeFaces.size(); ++i ) {
                    vector3 face = cubeFaces[i];
                    outMesh.add_face( vector3( cubeVerts[face.x], cubeVerts[face.y], cubeVerts[face.z] ) );
                }
            }
            ++index;
        }
    }
}

/**
 *	Generates the faces for the defined cube cases in a sparse plane by traversing the currentCubeCase
 *	array sparsely using the currentCubeCasesRLP rle plane structure.
 *	The vertices in the mesh to be used for computed cube case faces are all included in the appropriate entries
 *	of the previousVertexIndices, currentVertexIndices.
 *
 *	@param	mct						The marching cubes table used to determine cube case
 *faces.
 *	@param	currentCubeCases		The defined cube cases to be processed.
 *	@param	currentCubeCasesRLP		The rle plane for traversing the defined cube cases.
 *	@param	previousVertexIndices	Indices from the previous plane of the verts in the mesh to be used for cube
 *case faces.
 *	@param	currentVertexIndices	Indices from the current plane of the verts in the mesh to be used for cube case
 *faces.
 *	@param	outMesh					The output mesh to add faces to, should already have the verts
 *required for the faces.
 */
void generate_faces_for_sparse_plane( const marching_cubes_table& mct,
                                      const boost::shared_array<unsigned char> currentCubeCases,
                                      const rle_plane& currentCubeCasesRLP,
                                      const boost::shared_array<vector3>& previousVertexIndices,
                                      const boost::shared_array<vector3>& currentVertexIndices, trimesh3& outMesh ) {

    // if the plane is empty, no work to be done
    if( currentCubeCasesRLP.is_empty() )
        return;

    int xsize = currentCubeCasesRLP.get_xy_extents().xsize();

    // Allocate a vector for the 12 cube vertices, and the faces
    vector<int> cubeVerts( 12, 0 );

    std::pair<int, int> run;
    for( int runIndex = 0; runIndex < (int)currentCubeCasesRLP.size(); ++runIndex ) {
        run = currentCubeCasesRLP[runIndex];
        for( int i = run.first; i <= run.second; ++i ) {
            unsigned char cubeCase = currentCubeCases[i];
            if( cubeCase != 0x00 && cubeCase != 0xff ) {

                if( previousVertexIndices[i - xsize].x >= 0 ) { // edge 01
                    cubeVerts[0] = previousVertexIndices[i - xsize].x;
                }
                if( previousVertexIndices[i - 1].y >= 0 ) { // edge 02
                    cubeVerts[1] = previousVertexIndices[i - 1].y;
                }
                if( currentVertexIndices[i - xsize - 1].z >= 0 ) { // edge 04
                    cubeVerts[2] = currentVertexIndices[i - xsize - 1].z;
                }
                if( previousVertexIndices[i].y >= 0 ) { // edge 13
                    cubeVerts[3] = previousVertexIndices[i].y;
                }
                if( currentVertexIndices[i - xsize].z >= 0 ) { // edge 15
                    cubeVerts[4] = currentVertexIndices[i - xsize].z;
                }
                if( previousVertexIndices[i].x >= 0 ) { // edge 23
                    cubeVerts[5] = previousVertexIndices[i].x;
                }
                if( currentVertexIndices[i - 1].z >= 0 ) { // edge 26
                    cubeVerts[6] = currentVertexIndices[i - 1].z;
                }
                if( currentVertexIndices[i].z >= 0 ) { // edge 37
                    cubeVerts[7] = currentVertexIndices[i].z;
                }
                if( currentVertexIndices[i - xsize].x >= 0 ) { // edge 45
                    cubeVerts[8] = currentVertexIndices[i - xsize].x;
                }
                if( currentVertexIndices[i - 1].y >= 0 ) { // edge 46
                    cubeVerts[9] = currentVertexIndices[i - 1].y;
                }
                if( currentVertexIndices[i].y >= 0 ) { // edge 57
                    cubeVerts[10] = currentVertexIndices[i].y;
                }
                if( currentVertexIndices[i].x >= 0 ) { // edge 67
                    cubeVerts[11] = currentVertexIndices[i].x;
                }

                // hand it to the table to get the faces
                const vector<vector3>& cubeFaces = mct.get_cubecase_faces( cubeCase );

                // add the faces to the mesh
                for( unsigned i = 0; i < cubeFaces.size(); ++i ) {
                    vector3 face = cubeFaces[i];
                    outMesh.add_face( vector3( cubeVerts[face.x], cubeVerts[face.y], cubeVerts[face.z] ) );
                }
            }
        }
    }
}

// this table doesn't really work
static const boost::int8_t blockColorFetchBoundaryCase[8][26] = {
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 12, -1, 13, 12, -1, 13, 12, -1, 13, 12, -1, 13, 12, 13, 12, -1, 13, 12, -1, 13, 12, -1, 13, 12, -1, 13 },
    { 10, 10, 10, -1, -1, -1, 15, 15, 15, 10, 10, 10, -1, -1, 15, 15, 15, 10, 10, 10, -1, -1, -1, 15, 15, 15 },
    { 9, 10, 11, 12, -1, 13, 14, 15, 16, 9, 10, 11, 12, 13, 14, 15, 16, 9, 10, 11, 12, -1, 13, 14, 15, 16 },
    { 4, 4, 4, 4, 4, 4, 4, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, 21, 21, 21, 21, 21, 21, 21, 21, 21 },
    { 3, 4, 5, 3, 4, 5, 3, 4, 5, 12, -1, 13, 12, 13, 12, -1, 13, 20, 21, 22, 20, 21, 22, 20, 21, 22 },
    { 1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10, -1, -1, 15, 15, 15, 18, 18, 18, 21, 21, 21, 24, 24, 24 },
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 } };

// bool must_fetch_block_boundary( int blockColor, int boundaryCase ) {
// return blockColorFetchBoundaryCase[blockColor][boundaryCase] != 0;
//}

static const boost::uint8_t blockIsBoundaryWriter[26] = { 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0,
                                                          1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1 };

bool is_boundary_writer( int /*blockColor*/, int boundaryCase ) { return blockIsBoundaryWriter[boundaryCase] != 0; }

int get_boundary_owner_direction_code( int blockColor, int boundaryCase ) {
    return blockColorFetchBoundaryCase[blockColor][boundaryCase];
}

static const vector3 boundaryCaseNeighborOffset[26] = { vector3( -1, -1, -1 ), // 0
                                                        vector3( 0, -1, -1 ),  // 1
                                                        vector3( 1, -1, -1 ),  // 2
                                                        vector3( -1, 0, -1 ),  // 3
                                                        vector3( 0, 0, -1 ),   // 4
                                                        vector3( 1, 0, -1 ),   // 5
                                                        vector3( -1, 1, -1 ),  // 6
                                                        vector3( 0, 1, -1 ),   // 7
                                                        vector3( 1, 1, -1 ),   // 8

                                                        vector3( -1, -1, 0 ), // 9
                                                        vector3( 0, -1, 0 ),  // 10
                                                        vector3( 1, -1, 0 ),  // 11
                                                        vector3( -1, 0, 0 ),  // 12
                                                        vector3( 1, 0, 0 ),   // 13
                                                        vector3( -1, 1, 0 ),  // 14
                                                        vector3( 0, 1, 0 ),   // 15
                                                        vector3( 1, 1, 0 ),   // 16

                                                        vector3( -1, -1, 1 ), // 17
                                                        vector3( 0, -1, 1 ),  // 18
                                                        vector3( 1, -1, 1 ),  // 19
                                                        vector3( -1, 0, 1 ),  // 20
                                                        vector3( 0, 0, 1 ),   // 21
                                                        vector3( 1, 0, 1 ),   // 22
                                                        vector3( -1, 1, 1 ),  // 23
                                                        vector3( 0, 1, 1 ),   // 24
                                                        vector3( 1, 1, 1 ) };

frantic::graphics::vector3 get_direction_neighbor_offset( int direction ) {
    return boundaryCaseNeighborOffset[direction];
}

void copy_color_vertices( const std::vector<exposed_voxel_vertices>& colorExposedVertices,
                          voxel_vertex_t& exposedVoxelVertices ) {
    BOOST_FOREACH( const exposed_voxel_vertices& x, colorExposedVertices ) {
        std::pair<voxel_vertex_t::iterator, bool> p =
            exposedVoxelVertices.insert( std::pair<vector3, vector3>( x.voxelCoord, x.vertices ) );
        if( !p.second ) {
            // if this key was already in the table..
            voxel_vertex_t::iterator i = p.first;

            vector3 vertices = i->second;
            bool hasUpdate = false;
            for( int axis = 0; axis < 3; ++axis ) {
                const boost::int32_t oldVertex = vertices[axis];
                boost::int32_t newVertex = x.vertices[axis];
                if( !is_vertex_number( oldVertex ) && is_vertex_number( newVertex ) ) {
                    vertices[axis] = newVertex;
                    hasUpdate = true;
                }
            }
            if( hasUpdate ) {
                i->second = vertices;
            }
        }
    }
}

} // namespace detail

namespace {

struct channels_descending_order_optimized {
    bool operator()( const channel& lhs, const channel& rhs ) {
        if( lhs.name() == _T( "Position" ) )
            return true;
        else if( rhs.name() == _T( "Position" ) )
            return false;
        else if( lhs.name() == _T( "Radius" ) )
            return true;
        else if( rhs.name() == _T( "Radius" ) )
            return false;

        size_t lhsSize = sizeof_channel_data_type( lhs.data_type() );
        size_t rhsSize = sizeof_channel_data_type( rhs.data_type() );

        if( lhsSize == rhsSize ) {
            if( lhs.data_type() == rhs.data_type() ) {
                return lhs.arity() > rhs.arity();
            } else {
                return lhs.data_type() > rhs.data_type();
            }
        } else {
            return lhsSize > rhsSize;
        }
    }
}; // channel_descending_order_optimized

} // Anonymous namespace

channel_map create_optimized_channel_map( const channel_map& channelMap, std::size_t alignmentPadding ) {

    channel_map optimizedChannels;
    optimizedChannels.end_channel_definition( alignmentPadding );
    std::vector<channel> channelsToAdd;

    // Get all the channels in the other channel map so that they can be ordered.
    for( size_t i = 0; i < channelMap.channel_count(); ++i ) {
        channelsToAdd.push_back( channelMap[i] );
    }

    std::stable_sort( channelsToAdd.begin(), channelsToAdd.end(), channels_descending_order_optimized() );

    // Add the channels to the optimized channel_map
    for( std::size_t i = 0; i < channelsToAdd.size(); ++i ) {
        optimizedChannels.append_channel( channelsToAdd[i].name(), channelsToAdd[i].arity(),
                                          channelsToAdd[i].data_type() );
    }

    return optimizedChannels;
}

} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic
