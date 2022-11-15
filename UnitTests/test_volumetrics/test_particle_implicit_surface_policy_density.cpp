// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/algorithm/string/predicate.hpp>

#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/particles/streams/particle_array_particle_istream.hpp>
#include <frantic/volumetrics/implicitsurface/calculate_particle_anisotropic_params.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/levelset/rle_level_set_file_io.hpp>

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif

using namespace frantic::channels;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::logging;
using namespace frantic::math;
using namespace frantic::particles;
using namespace frantic::particles::streams;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::implicitsurface;
using namespace frantic::volumetrics::levelset;

namespace {

void populate_particle_grid_tree( std::size_t particleCount, frantic::particles::particle_grid_tree& outParticles ) {
    const float radius = 1;

    channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<float>( _T("Radius") );
    channelMap.end_channel_definition();

    outParticles.reset( channelMap, 1 );

    std::vector<char> buffer( channelMap.structure_size() );

    channel_accessor<vector3f> positionAcc = channelMap.get_accessor<vector3f>( _T("Position") );
    channel_accessor<float> radiusAcc = channelMap.get_accessor<float>( _T("Radius") );

    if( particleCount >= 1 ) {
        positionAcc( buffer ) = vector3f( -0.9f );
        radiusAcc( buffer ) = radius;

        outParticles.insert( &buffer[0] );
    }

    if( particleCount >= 2 ) {
        positionAcc( buffer ) = vector3f( 0, 0.5f, 1.f );
        radiusAcc( buffer ) = 0.75f * radius;

        outParticles.insert( &buffer[0] );
    }

    if( particleCount > 2 ) {
        throw std::runtime_error( "populate_particle_grid_tree Error: unsupported particle count:" +
                                  boost::lexical_cast<std::string>( particleCount ) );
    }
}

template <class ImplicitSurfacePolicy>
void compute_density_impl( ImplicitSurfacePolicy& isp, const std::string& samplingMode, const boundbox3& voxelExtents,
                           std::vector<float>& outDensity ) {
    typename ImplicitSurfacePolicy::sparse_voxel_corner_density_workspace_t data;

    if( samplingMode == "fill_sparse_voxel_corner_densities_2d" ) {
        if( voxelExtents.zminimum() != voxelExtents.zmaximum() ) {
            throw std::runtime_error( "compute_density_impl Error: " + samplingMode +
                                      " requires zminimum == zmaximum" );
        }

        const boundrect2 xyExtents( vector2( voxelExtents.xminimum(), voxelExtents.yminimum() ),
                                    vector2( voxelExtents.xmaximum(), voxelExtents.ymaximum() ) );
        const int z = voxelExtents.zminimum();
        isp.fill_sparse_voxel_corner_densities( xyExtents, z, &outDensity[0], data );
    } else if( samplingMode == "fill_voxel_corner_densities_3d" ) {
        isp.fill_voxel_corner_densities( voxelExtents, &outDensity[0], data );
    } else {
        throw std::runtime_error( "compute_density_impl Error: unknown sampling mode: " + samplingMode );
    }
}

void compute_density( const std::string& meshingMode, const std::string& samplingMode, std::size_t particleCount,
                      const voxel_coord_system& meshingVCS, const boundbox3& voxelExtents,
                      std::vector<float>& outDensity, float& outDefaultOutsideDensity ) {
    const float voxelLength = meshingVCS.voxel_length();
    const float maximumParticleRadius = 1;
    const int vertexRefinement = 10;
    const float radius = 1;

    particle_grid_tree particles;
    populate_particle_grid_tree( particleCount, particles );

    outDensity.resize( voxelExtents.get_volume() );

    if( meshingMode == "union_of_spheres" ) {
        const float particleRadiusToEffectRadiusScale = 1.5;
        const float implicitThreshold = 0.3f;

        particle_union_of_spheres_is_policy policy( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale,
                                                    implicitThreshold, meshingVCS, vertexRefinement );

        compute_density_impl( policy, samplingMode, voxelExtents, outDensity );

        outDefaultOutsideDensity = policy.get_default_outside_distance();
    } else if( meshingMode == "metaball" ) {
        const float effectRadiusScale = 1 + 2 * voxelLength / radius;
        const float implicitThreshold = 1.5f * square( 1 - 1 / effectRadiusScale );

        particle_metaball_is_policy policy( particles, maximumParticleRadius, effectRadiusScale, implicitThreshold,
                                            meshingVCS, vertexRefinement );

        compute_density_impl( policy, samplingMode, voxelExtents, outDensity );

        outDefaultOutsideDensity = policy.get_default_outside_distance();
    } else if( meshingMode == "zhu_bridson" ) {
        const float effectRadiusScale = 1.7f;
        const float lowDensityTrimmingDensity = 1;
        const float lowDensityTrimmingStrength = 15;

        particle_zhu_bridson_is_policy policy( particles, maximumParticleRadius, effectRadiusScale,
                                               lowDensityTrimmingDensity, lowDensityTrimmingStrength, meshingVCS,
                                               vertexRefinement );

        compute_density_impl( policy, samplingMode, voxelExtents, outDensity );

        outDefaultOutsideDensity = policy.get_default_outside_distance();
    } else if( meshingMode == "anisotropic" ) {
        const float compactSupportScale = 4;
        // Default parameters in Frost
        const std::size_t minNeighborCount = 25;
        const float implicitThreshold = 0.5;
        const float maxAnisotropy = 4;
        const float anisotropyWindowScale = 2;

        // Need to rebuild the particles with anisotropy parameters
        const frantic::tstring volumeChannelName = _T("__Volume");
        channel_map anisoChannelMap =
            create_channel_map_with_anisotropy_channels( particles.get_channel_map(), volumeChannelName );

        particle_array particleArray( anisoChannelMap );
        particleArray.insert_particles( particles.get_channel_map(), particles.begin(), particles.end() );

        null_progress_logger progressLogger;

        calculate_anisotropy( particleArray, compactSupportScale, compactSupportScale * anisotropyWindowScale,
                              maxAnisotropy, minNeighborCount, progressLogger );

        calculate_volume_with_anisotropic_kernel( particleArray, volumeChannelName, progressLogger );

        particles.reset( anisoChannelMap, 1 );
        particle_istream_ptr pin( new particle_array_particle_istream( particleArray ) );
        particles.insert_particles( pin );

        particle_anisotropic_is_policy policy( particles, implicitThreshold, meshingVCS, vertexRefinement );

        compute_density_impl( policy, samplingMode, voxelExtents, outDensity );

        outDefaultOutsideDensity = policy.get_default_outside_distance();
    } else {
        throw std::runtime_error( "compute_density Error: unknown meshing mode: " + meshingMode );
    }
}

frantic::tstring get_reference_filename( const std::string& meshingMode, std::size_t particleCount ) {
    return _T("TestInputs/density_") + frantic::strings::to_tstring( meshingMode ) + _T("_") +
           boost::lexical_cast<frantic::tstring>( particleCount ) + _T(".rls");
}

void save_reference_density( const std::string& meshingMode, std::size_t particleCount,
                             const voxel_coord_system& meshingVCS, const boundbox3& voxelExtents,
                             const std::vector<float>& density, float defaultOutsideDensity ) {
    std::vector<vector3> voxelCoords;
    for( boost::int32_t z = voxelExtents.zminimum(); z <= voxelExtents.zmaximum(); ++z ) {
        for( boost::int32_t y = voxelExtents.yminimum(); y <= voxelExtents.ymaximum(); ++y ) {
            for( boost::int32_t x = voxelExtents.xminimum(); x <= voxelExtents.xmaximum(); ++x ) {
                voxelCoords.push_back( vector3( x, y, z ) );
            }
        }
    }

    rle_index_spec denseRis;
    denseRis.build_from_voxel_array( voxelCoords );

    std::vector<boost::uint8_t> populated( density.size() );
    std::vector<float> populatedDensity;

    for( std::size_t i = 0; i < density.size(); ++i ) {
        if( density[i] != defaultOutsideDensity ) {
            populated[i] = 1;
            populatedDensity.push_back( density[i] );
        }
    }

    rle_index_spec ris;
    ris.build_from_populated( denseRis, &populated[0], &density[0], false );

    const float interfaceWidthInside = std::numeric_limits<float>::max();
    // Want the level set's outside distance to equal the defaultOutsideDensity
    // so that we'll get the desired value if we
    const float interfaceWidthOutside = defaultOutsideDensity / meshingVCS.voxel_length() - 1;
    rle_level_set rls( meshingVCS, ris, populatedDensity, interfaceWidthInside, interfaceWidthOutside );

    ASSERT_EQ( defaultOutsideDensity, rls.get_outside_distance() );

    write_rls_rle_level_set_file( get_reference_filename( meshingMode, particleCount ), rls );
}

void load_reference_density( const std::string& meshingMode, std::size_t particleCount, const boundbox3& voxelExtents,
                             std::vector<float>& outDensity ) {
    rle_level_set rls;
    read_rle_level_set_file( get_reference_filename( meshingMode, particleCount ), rls );
    rls.fill_box( voxelExtents, outDensity );
}

boundbox3 get_voxel_extents( const std::string& samplingMode, boost::int32_t extentsWidth ) {
    using boost::algorithm::ends_with;

    const boundbox3 extents3d( vector3( -8 ), size3( extentsWidth ) );
    if( ends_with( samplingMode, "2d" ) ) {
        const boost::int32_t z = ( extents3d.zminimum() + extents3d.zmaximum() ) / 2;

        boundbox3 result( extents3d );
        result.minimum().z = z;
        result.maximum().z = z;

        return result;
    } else {
        return extents3d;
    }
}

} // anonymous namespace

class ParticleImplicitSurfacePolicyDensity
    : public ::testing::TestWithParam<std::tuple<std::string, std::string, std::size_t, boost::int32_t>> {};

INSTANTIATE_TEST_CASE_P(
    ParticleImplicitSurfacePolicyDensity, ParticleImplicitSurfacePolicyDensity,
    ::testing::Combine( ::testing::Values( "union_of_spheres", "metaball", "zhu_bridson", "anisotropic" ),
                        ::testing::Values( "fill_voxel_corner_densities_3d", "fill_sparse_voxel_corner_densities_2d" ),
                        ::testing::Range<std::size_t>( 0, 3 ), ::testing::Values<boost::int32_t>( 16, 15 ) ) );

TEST_P( ParticleImplicitSurfacePolicyDensity, CompareWithReference ) {
    const std::string meshingMode = std::get<0>( GetParam() );
    const std::string samplingMode = std::get<1>( GetParam() );
    const std::size_t particleCount = std::get<2>( GetParam() );
    const boost::int32_t extentsWidth = std::get<3>( GetParam() );

    const float voxelLength = 0.35f;
    const voxel_coord_system meshingVCS( vector3f( 0 ), voxelLength );

    const boundbox3 voxelExtents = get_voxel_extents( samplingMode, extentsWidth );

    std::vector<float> density;
    float defaultOutsideDensity;
    compute_density( meshingMode, samplingMode, particleCount, meshingVCS, voxelExtents, density,
                     defaultOutsideDensity );

    // Uncomment this if you want to save new reference data.
    /*
    if( samplingMode == "fill_voxel_corner_densities_3d" && extentsWidth == 16 ) {
      save_reference_density(
        meshingMode, particleCount, meshingVCS, voxelExtents, density, defaultOutsideDensity );
    }
    */

    bool foundDefault = true;
    bool foundNonDefault = false;
    for( std::size_t i = 0; i < density.size(); ++i ) {
        if( density[i] == defaultOutsideDensity ) {
            foundDefault = true;
        } else {
            foundNonDefault = true;
        }
    }

    if( particleCount == 0 ) {
        EXPECT_TRUE( foundDefault );
        EXPECT_FALSE( foundNonDefault );
    } else {
        // We don't necessarily need to have default outside density values in
        // the test, but currently I do expect to find them.
        EXPECT_TRUE( foundDefault );
        EXPECT_TRUE( foundNonDefault );
    }

    std::vector<float> referenceDensity;
    load_reference_density( meshingMode, particleCount, voxelExtents, referenceDensity );

    ASSERT_EQ( density.size(), referenceDensity.size() );
    EXPECT_TRUE( density.size() > 0 );

    std::size_t i = 0;
    for( boost::int32_t z = voxelExtents.zminimum(); z <= voxelExtents.zmaximum(); ++z ) {
        for( boost::int32_t y = voxelExtents.yminimum(); y <= voxelExtents.ymaximum(); ++y ) {
            for( boost::int32_t x = voxelExtents.xminimum(); x <= voxelExtents.xmaximum(); ++x ) {
                const float error = std::abs( density[i] - referenceDensity[i] );
                const float scale = std::abs( density[i] ) + std::abs( referenceDensity[i] );

                const float relativeError = error / scale;

                EXPECT_LE( relativeError, 1e-3f )
                    << "Differing density in voxel " << vector3( x, y, z ) << ".  "
                    << "Expected " << referenceDensity[i] << " but got " << density[i] << " instead.";

                ++i;
            }
        }
    }
    EXPECT_EQ( density.size(), i );
}
