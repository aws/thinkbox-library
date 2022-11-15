// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#ifndef NOMINMAX
#define NOMINMAX
#endif

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif

#include <tbb/task_scheduler_init.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/make_shared.hpp>

#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/mesh_measurement.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>
#include <frantic/particles/streams/particle_array_particle_istream.hpp>
#include <frantic/volumetrics/implicitsurface/calculate_particle_anisotropic_params.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>

#include "gtest-helper.h"
#include "utilities/mesh_generators.hpp"

using namespace frantic::channels;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::logging;
using namespace frantic::math;
using namespace frantic::particles;
using namespace frantic::particles::streams;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::implicitsurface;

using boost::algorithm::ends_with;
using boost::algorithm::starts_with;

namespace {

std::vector<std::string> get_all_meshing_modes() {
    std::vector<std::string> result;
    result.push_back( "union_of_spheres_plane" );
    result.push_back( "union_of_spheres_block" );
    result.push_back( "metaball_plane" );
    result.push_back( "metaball_block" );
    result.push_back( "zhu_bridson_plane" );
    result.push_back( "zhu_bridson_block" );
    result.push_back( "anisotropic_plane" );
    result.push_back( "anisotropic_block" );
    return result;
}

struct meshing_parameters {
    meshing_parameters()
        : vertexRefinement( 0 )
        , createDensityChannel( false ) {}

    std::string meshingMode;
    int vertexRefinement;
    bool createDensityChannel;
};

template <class ImplicitSurfacePolicy>
void populate_density_channel( const ImplicitSurfacePolicy& isp, trimesh3& mesh ) {
    const frantic::tstring channelName( _T("Density") );

    mesh.add_vertex_channel<float>( channelName );

    trimesh3_vertex_channel_accessor<float> acc = mesh.get_vertex_channel_accessor<float>( channelName );

    for( std::size_t i = 0, ie = mesh.vertex_count(); i < ie; ++i ) {
        const vector3f position = mesh.get_vertex( i );
        const float density = isp.get_density( position );
        acc[i] = density;
    }
}

void create_red_particles_near_origin( particle_array& out, std::size_t particleCount ) {
    const float boxEdgeLength = 0.01f;

    channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<vector3f>( _T("Color") );
    channelMap.end_channel_definition();

    out.reset( channelMap );

    std::vector<char> buffer( channelMap.structure_size() );

    channel_accessor<vector3f> positionAcc = channelMap.get_accessor<vector3f>( _T("Position") );
    channel_accessor<vector3f> colorAcc = channelMap.get_accessor<vector3f>( _T("Color") );

    positionAcc( buffer ) = vector3f( 0 );
    colorAcc( buffer ) = vector3f( 1, 0, 0 );

    out.push_back( &buffer[0] );

    if( particleCount > 1 ) {
        // Additional particles, closely spaced so that vertex refinement
        // uses the SIMD code path.
        boundbox3f box( vector3f( -boxEdgeLength / 2 ), vector3f( boxEdgeLength / 2 ) );
        for( std::size_t i = 0, ie = particleCount - 1; i < ie; ++i ) {
            positionAcc( buffer ) = box.get_corner( static_cast<int>( i ) );
            out.push_back( &buffer[0] );
        }
    }
}

void create_unit_radius_particle_mesh( const meshing_parameters& params, const particle_array& particleArray,
                                       trimesh3& outMesh ) {
    outMesh.clear();

    const std::string& meshingMode = params.meshingMode;

    float radius = 1;
    if( starts_with( meshingMode, "zhu_bridson" ) ) {
        // Need to increase the radius to compensate for low density trimming.
        // I just eyeballed this in Frost MX.
        radius *= 1.45f;
    }

    channel_map channelMap = particleArray.get_channel_map();
    if( !channelMap.has_channel( _T("Radius") ) ) {
        channelMap.append_channel<float>( _T("Radius") );
    }

    const float voxelLength = 0.1f;
    const int vertexRefinement = params.vertexRefinement;

    particle_grid_tree particles( channelMap );

    particle_istream_ptr pin = boost::make_shared<particle_array_particle_istream>( particleArray );
    particles.insert_particles( pin );

    channel_accessor<float> radiusAcc = channelMap.get_accessor<float>( _T("Radius") );
    for( particle_grid_tree::iterator i = particles.begin(), ie = particles.end(); i != ie; ++i ) {
        radiusAcc( *i ) = radius;
    }

    channel_propagation_policy cpp;

    const voxel_coord_system meshingVCS( vector3f( 0 ), voxelLength );

    null_progress_logger progressLogger;

    if( starts_with( meshingMode, "union_of_spheres" ) ) {
        // want to get at least one voxel with meaningful density outside
        // of the particle
        const float particleRadiusToEffectRadiusScale = 1 + 2 * voxelLength / radius;
        const float implicitThreshold = 2 * voxelLength;

        if( ends_with( meshingMode, "_plane" ) ) {
            union_of_spheres_convert_particles_to_trimesh3( particles, cpp, radius, particleRadiusToEffectRadiusScale,
                                                            implicitThreshold, meshingVCS, vertexRefinement, outMesh );
        } else if( ends_with( meshingMode, "_block" ) ) {
            union_of_spheres_convert_sparse_particles_to_trimesh3(
                particles, cpp, radius, particleRadiusToEffectRadiusScale, implicitThreshold, meshingVCS,
                vertexRefinement, outMesh, progressLogger );
        } else {
            throw std::runtime_error( "create_unit_radius_particle_mesh Error: unknown meshing mode: " + meshingMode );
        }

        if( params.createDensityChannel ) {
            particle_union_of_spheres_is_policy policy( particles, radius, particleRadiusToEffectRadiusScale,
                                                        implicitThreshold, meshingVCS, vertexRefinement );
            populate_density_channel( policy, outMesh );
        }
    } else if( starts_with( meshingMode, "metaball" ) ) {
        const float effectRadiusScale = 1 + 2 * voxelLength / radius;
        // based on metaball_function(); intended to get unit radius
        const float implicitThreshold = 1.5f * square( 1 - 1 / effectRadiusScale );

        if( ends_with( meshingMode, "_plane" ) ) {
            metaball_convert_particles_to_trimesh3( particles, cpp, radius, effectRadiusScale, implicitThreshold,
                                                    meshingVCS, vertexRefinement, outMesh );
        } else if( ends_with( meshingMode, "_block" ) ) {
            metaball_convert_sparse_particles_to_trimesh3( particles, cpp, radius, effectRadiusScale, implicitThreshold,
                                                           meshingVCS, vertexRefinement, outMesh, progressLogger );
        } else {
            throw std::runtime_error( "create_unit_radius_particle_mesh Error: unknown meshing mode: " + meshingMode );
        }

        if( params.createDensityChannel ) {
            particle_metaball_is_policy policy( particles, radius, effectRadiusScale, implicitThreshold, meshingVCS,
                                                vertexRefinement );
            populate_density_channel( policy, outMesh );
        }
    } else if( starts_with( meshingMode, "zhu_bridson" ) ) {
        const float effectRadiusScale = 1.7f;
        const float lowDensityTrimmingDensity = 1;
        const float lowDensityTrimmingStrength = 15;

        if( ends_with( meshingMode, "_plane" ) ) {
            zhu_bridson_convert_particles_to_trimesh3( particles, cpp, radius, effectRadiusScale,
                                                       lowDensityTrimmingDensity, lowDensityTrimmingStrength,
                                                       meshingVCS, vertexRefinement, outMesh );
        } else if( ends_with( meshingMode, "_block" ) ) {
            zhu_bridson_convert_sparse_particles_to_trimesh3( particles, cpp, radius, effectRadiusScale,
                                                              lowDensityTrimmingDensity, lowDensityTrimmingStrength,
                                                              meshingVCS, vertexRefinement, outMesh, progressLogger );
        } else {
            throw std::runtime_error( "create_unit_radius_particle_mesh Error: unknown meshing mode: " + meshingMode );
        }

        if( params.createDensityChannel ) {
            particle_zhu_bridson_is_policy policy( particles, radius, effectRadiusScale, lowDensityTrimmingDensity,
                                                   lowDensityTrimmingStrength, meshingVCS, vertexRefinement );
            populate_density_channel( policy, outMesh );
        }
    } else if( starts_with( meshingMode, "anisotropic" ) ) {
        // I just eyeballed this in Frost MX
        const float compactSupportScale = 6.6f;
        // Default parameters in Frost
        const std::size_t minNeighborCount = 25;
        const float implicitThreshold = 0.5;
        const float maxAnisotropy = 4;
        const float anisotropyWindowScale = 2;

        // Need to rebuild the particles with anisotropy parameters
        const frantic::tstring volumeChannelName = _T("__Volume");
        channel_map anisoChannelMap = create_channel_map_with_anisotropy_channels( channelMap, volumeChannelName );

        particle_array particleArray( anisoChannelMap );
        particleArray.insert_particles( particles.get_channel_map(), particles.begin(), particles.end() );

        calculate_anisotropy( particleArray, compactSupportScale, compactSupportScale * anisotropyWindowScale,
                              maxAnisotropy, minNeighborCount, progressLogger );

        calculate_volume_with_anisotropic_kernel( particleArray, volumeChannelName, progressLogger );

        particles.reset( anisoChannelMap, 1 );
        particle_istream_ptr pin( new particle_array_particle_istream( particleArray ) );
        particles.insert_particles( pin );

        if( ends_with( meshingMode, "_plane" ) ) {
            anisotropic_convert_particles_to_trimesh3( particles, cpp, implicitThreshold, meshingVCS, vertexRefinement,
                                                       outMesh );
        } else if( ends_with( meshingMode, "_block" ) ) {
            anisotropic_convert_sparse_particles_to_trimesh3( particles, cpp, implicitThreshold, meshingVCS,
                                                              vertexRefinement, outMesh, progressLogger );
        } else {
            throw std::runtime_error( "create_unit_radius_particle_mesh Error: unknown meshing mode: " + meshingMode );
        }

        if( params.createDensityChannel ) {
            particle_anisotropic_is_policy policy( particles, implicitThreshold, meshingVCS, vertexRefinement );
            populate_density_channel( policy, outMesh );
        }
    } else {
        throw std::runtime_error( "create_unit_radius_particle_mesh Error: unknown meshing mode: " + meshingMode );
    }
}

} // anonymous namespace

class ConvertParticlesToTrimesh3 : public ::testing::TestWithParam<std::string> {};

TEST_P( ConvertParticlesToTrimesh3, SingleParticleUnitRadiusMesh ) {
    tbb::task_scheduler_init taskScheduler;

    trimesh3 mesh;

    meshing_parameters params;
    params.meshingMode = GetParam();
    params.vertexRefinement = 10;

    particle_array particles;
    create_red_particles_near_origin( particles, 1 );

    create_unit_radius_particle_mesh( params, particles, mesh );

    mesh_interface::ptr_type meshInterface( trimesh3_interface::create_instance( &mesh ).release() );

    EXPECT_TRUE( mesh.face_count() > 0 );

    EXPECT_TRUE( is_closed_manifold( meshInterface ) );

    ASSERT_TRUE( mesh.has_vertex_channel( _T("Color" ) ) );

    trimesh3_vertex_channel_accessor<vector3f> colorAcc = mesh.get_vertex_channel_accessor<vector3f>( _T("Color") );
    for( std::size_t i = 0; i < colorAcc.size(); ++i ) {
        EXPECT_VECTOR3F_EQ( vector3f( 1, 0, 0 ), colorAcc[i] );
    }

    // Expect all of the normal vectors to point outward
    mesh.build_vertex_normals();
    trimesh3_vertex_channel_accessor<vector3f> normalAcc = mesh.get_vertex_channel_accessor<vector3f>( _T("Normal") );
    for( std::size_t i = 0; i < normalAcc.size(); ++i ) {
        EXPECT_TRUE( vector3f::dot( normalAcc[i], mesh.get_vertex( i ) ) > 0 );
    }

    trimesh3 sphereMesh;
    make_sphere_mesh( 50, sphereMesh );

    EXPECT_LT( hausdorff_distance_two_sided( sphereMesh, mesh, true, true, 1000 ), 0.01 );
}

INSTANTIATE_TEST_CASE_P( ConvertParticlesToTrimesh3, ConvertParticlesToTrimesh3,
                         ::testing::ValuesIn( get_all_meshing_modes() ) );

class TestMeshingModeAndParticleCount : public ::testing::TestWithParam<std::tuple<std::string, std::size_t>> {};

namespace {

void get_sorted_density( const frantic::geometry::trimesh3& mesh, std::vector<float>& outDensity ) {
    outDensity.clear();

    const_trimesh3_vertex_channel_accessor<float> acc = mesh.get_vertex_channel_accessor<float>( _T("Density") );

    for( std::size_t i = 0, ie = mesh.vertex_count(); i < ie; ++i ) {
        outDensity.push_back( std::abs( acc[i] ) );
    }

    std::sort( outDensity.begin(), outDensity.end() );
}

} // anonymous namespace

// Test whether vertex refinement moves the mesh closer to zero
// density, that is, the true surface of the level set.
TEST_P( TestMeshingModeAndParticleCount, VertRefinement ) {
    tbb::task_scheduler_init taskScheduler;

    meshing_parameters params;
    params.meshingMode = std::get<0>( GetParam() );
    params.createDensityChannel = true;

    const std::size_t particleCount = std::get<1>( GetParam() );

    particle_array particles;
    create_red_particles_near_origin( particles, particleCount );

    trimesh3 mesh0;
    params.vertexRefinement = 0;
    create_unit_radius_particle_mesh( params, particles, mesh0 );

    trimesh3 mesh1;
    params.vertexRefinement = 1;
    create_unit_radius_particle_mesh( params, particles, mesh1 );

    // Currently the mesher can produce different vertex orders
    // between runs.  So, instead of comparing the densities directly,
    // I'm sorting them first, and then comparing the values after
    // sorting.
    std::vector<float> density0, density1;
    get_sorted_density( mesh0, density0 );
    get_sorted_density( mesh1, density1 );

    ASSERT_EQ( density0.size(), density1.size() );
    EXPECT_TRUE( density0.size() > 0 );
    EXPECT_TRUE( density0.back() > 0 );

    for( std::size_t i = 0; i < density0.size(); ++i ) {
        if( density0[0] != 0 ) {
            EXPECT_TRUE( density1[i] < density0[i] );
        }
    }
}

// Expect that the channel values from coincident particles are blended
// equally to form the output channel value.
TEST_P( TestMeshingModeAndParticleCount, SinglePositionChannelBlending ) {
    const std::string meshingMode = std::get<0>( GetParam() );
    const std::size_t particleCount = std::get<1>( GetParam() );

    channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel( _T("Data"), particleCount, data_type_float32 );
    channelMap.end_channel_definition();

    particle_array particles( channelMap );
    {
        std::vector<char> buffer( channelMap.structure_size() );
        channel_general_accessor dataAcc = channelMap.get_general_accessor( _T("Data") );

        // Set Data[particleIndex] to 1, so we can easily tell if the particle
        // contributes to the output channel value.
        for( std::size_t i = 0; i < particleCount; ++i ) {
            float* data = reinterpret_cast<float*>( dataAcc.get_channel_data_pointer( buffer ) );
            memset( data, 0, dataAcc.primitive_size() );
            data[i] = 1;
            particles.push_back( &buffer[0] );
        }
    }

    meshing_parameters params;
    params.meshingMode = meshingMode;
    params.createDensityChannel = true;

    trimesh3 mesh;
    params.vertexRefinement = 0;
    create_unit_radius_particle_mesh( params, particles, mesh );

    EXPECT_TRUE( mesh.vertex_count() > 0 );

    trimesh3_vertex_channel_general_accessor dataAcc = mesh.get_vertex_channel_general_accessor( _T("Data") );
    if( starts_with( meshingMode, "union_of_spheres" ) ) {
        // Union of spheres doesn't do any blending.  Instead, it copies the
        // channel values from one particle.  So, we should expect to find
        // a single 1 entry in the data, while all other entries are 0.
        for( std::size_t vertexIndex = 0; vertexIndex < mesh.vertex_count(); ++vertexIndex ) {
            std::size_t oneCount = 0;
            float* f = reinterpret_cast<float*>( dataAcc.data( vertexIndex ) );
            for( std::size_t i = 0; i < particleCount; ++i ) {
                if( f[i] == 1 ) {
                    ++oneCount;
                }
            }
            EXPECT_EQ( 1, oneCount );
        }
    } else {
        // The output value should be an equal blend of all input values.
        for( std::size_t vertexIndex = 0; vertexIndex < mesh.vertex_count(); ++vertexIndex ) {
            float* f = reinterpret_cast<float*>( dataAcc.data( vertexIndex ) );
            for( std::size_t i = 0; i < particleCount; ++i ) {
                EXPECT_NEAR( 1.f / particleCount, f[i], 0.001f );
            }
        }
    }
}

// Expect that nearby particles are weighted more heavily in blended channels.
TEST_P( TestMeshingModeAndParticleCount, TwoPositionChannelBlending ) {
    const std::string meshingMode = std::get<0>( GetParam() );
    const std::size_t particleCount = std::get<1>( GetParam() );

    channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<vector3f>( _T("Color") );
    channelMap.end_channel_definition();

    particle_array particles( channelMap );
    {
        std::vector<char> buffer( channelMap.structure_size() );
        channel_accessor<vector3f> positionAcc = channelMap.get_accessor<vector3f>( _T("Position") );
        channel_accessor<vector3f> colorAcc = channelMap.get_accessor<vector3f>( _T("Color") );

        // Red particles on the left
        positionAcc( &buffer[0] ).set( -0.1f, 0, 0 );
        colorAcc( &buffer[0] ).set( 1, 0, 0 );
        for( std::size_t i = 0; i < particleCount; ++i ) {
            particles.push_back( &buffer[0] );
        }

        // Green particles on the right
        positionAcc( &buffer[0] ).set( 0.1f, 0, 0 );
        colorAcc( &buffer[0] ).set( 0, 1, 0 );
        for( std::size_t i = 0; i < particleCount; ++i ) {
            particles.push_back( &buffer[0] );
        }
    }

    meshing_parameters params;
    params.meshingMode = meshingMode;
    params.createDensityChannel = true;

    trimesh3 mesh;
    params.vertexRefinement = 0;
    create_unit_radius_particle_mesh( params, particles, mesh );

    EXPECT_TRUE( mesh.vertex_count() > 0 );

    std::size_t leftVertexIndex = 0;
    float leftVertexPosition = std::numeric_limits<float>::max();
    std::size_t rightVertexIndex = 0;
    float rightVertexPosition = -std::numeric_limits<float>::max();

    for( std::size_t vertexIndex = 0; vertexIndex < mesh.vertex_count(); ++vertexIndex ) {
        const float x = mesh.get_vertex( vertexIndex ).x;
        if( x < leftVertexPosition ) {
            leftVertexIndex = vertexIndex;
            leftVertexPosition = x;
        }
        if( x > rightVertexPosition ) {
            rightVertexIndex = vertexIndex;
            rightVertexPosition = x;
        }
    }

    trimesh3_vertex_channel_accessor<vector3f> colorAcc = mesh.get_vertex_channel_accessor<vector3f>( _T("Color") );

    // Expect the left vertex to be more red
    EXPECT_TRUE( colorAcc[leftVertexIndex].x > colorAcc[rightVertexIndex].x );
    // And the right vertex to be more green
    EXPECT_TRUE( colorAcc[leftVertexIndex].y < colorAcc[rightVertexIndex].y );
}

INSTANTIATE_TEST_CASE_P( ConvertParticlesToTrimesh3MeshingModeAndParticleCount, TestMeshingModeAndParticleCount,
                         ::testing::Combine( ::testing::ValuesIn( get_all_meshing_modes() ),
                                             // Testing 1 particle to hit the scalar code path,
                                             // and 5 to hit the SIMD code path.
                                             ::testing::Values( 1, 5 ) ) );
