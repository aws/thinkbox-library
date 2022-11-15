// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/volumetrics/implicitsurface/level_set_implicit_surface_policies.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/rle_plane.hpp>

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

namespace detail {
/**
 * This class is used as a workaround for error C2888: symbol cannot be defined within namespace 'levelset'
 */
class brls_ris_friend {
  public:
    /**
     *	This function takes an implicit surface policy and populates the level set with
     *  the data from the policy.  The channel propagation policy dictates which channels
     *  are propagated from the is policy to the level set.
     */
    template <class ImplicitSurfacePolicy>
    static void build_rle_level_set( const frantic::channels::channel_propagation_policy& cpp,
                                     ImplicitSurfacePolicy& isp,
                                     frantic::volumetrics::levelset::rle_level_set& outLevelSet,
                                     frantic::logging::progress_logger& logger ) {
        if( !cpp.is_channel_included( _T("SignedDistance") ) )
            throw std::runtime_error(
                "implicitsurface::build_rle_level_set() - The provided channel propagation policy does not include a "
                "'SignedDistance' channel.  This is a required channel to construct a level set." );

        logger.check_for_abort();

        // Get the bounds of the policy, we'll use these to build the rle index spec.
        // These are voxel bounds in the meshingVCS of the isp, which is what we want.
        frantic::graphics::boundbox3 bounds = isp.get_voxel_bounds();
        frantic::graphics2d::boundrect2 xyExtents( bounds.xminimum(), bounds.xmaximum(), bounds.yminimum(),
                                                   bounds.ymaximum() );
        // std::cout << bounds.str() << std::endl;
        // std::cout << bounds.size().str() << std::endl;

        // Make an rle index spec and initialize the required info
        frantic::volumetrics::levelset::rle_index_spec ris;

        // We're creating a level set that matches the bounds of the particle system, so the bounds of the
        // index spec should match the bounds of the level set
        ris.m_abcCoordOrigin = bounds.minimum(); // same as the voxel origin of vcs
        ris.m_abcCoordSize = bounds.size();      // same as the bounds of the particle system

        // Initialize the recorded data size to 0, so when we create the named channels they start with a size of zero
        ris.m_dataSize = 0;

        // Set the exterior region code
        ris.m_exteriorRegionCode = isp.get_exterior_region_code();

        // std::cout << "ris.m_abcCoordSize.ysize() " << ris.m_abcCoordSize.ysize() << std::endl;
        // std::cout << "ris.m_abcCoordSize.zsize() " << ris.m_abcCoordSize.zsize() << std::endl;

        // Reserve the correct amount of space in the bc to run index mapping, and clear the run index data array
        ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
        ris.m_runData.clear();

        // Need a voxel value array to populate, and an rle_plane to be populated with the run indexing data
        // for each plane from the isp
        boost::shared_array<float> distanceData( new float[xyExtents.get_area()] );
        frantic::volumetrics::rle_plane rlp( xyExtents );

        // The interface widths.
        float interfaceVoxelWidthInside, interfaceVoxelWidthOutside;
        isp.get_interface_widths( interfaceVoxelWidthInside, interfaceVoxelWidthOutside );

        // std::cout << "fetching channel info" << std::endl;
        //  the channel info
        std::vector<frantic::tstring> channelName;
        std::vector<frantic::channels::data_type_t> channelType;
        std::vector<size_t> channelArity;
        std::vector<size_t> channelPrimitiveSize;
        isp.get_channel_info( cpp, channelName, channelType, channelArity, channelPrimitiveSize );

        // setup the data arrays for the channel data using the channel info

        std::vector<boost::shared_array<char>> channelData;
        std::vector<char*> channelDataPointers;
        int signedDistanceChannel = -1;
        // std::cout << channelName.size() << " channels" << std::endl;
        for( size_t i = 0; i < channelName.size(); ++i ) {
            // std::cout << "initializing channel " << channelName[i] <<  std::endl;
            // std::cout << "type: " << channelType[i] << std::endl;
            // std::cout << "arity: " << channelArity[i] << std::endl;
            // std::cout << "primitive size: " << channelPrimitiveSize[i] << std::endl;
            //  allocate space for the channel data we'll retrieve
            channelData.push_back(
                boost::shared_array<char>( new char[xyExtents.get_area() * channelPrimitiveSize[i]] ) );
            channelDataPointers.push_back( channelData[i].get() );

            if( channelName[i] == _T("SignedDistance") ) {
                signedDistanceChannel = (int)i;
            }
        }
        // std::cout << "done channel initialization" << std::endl;

        // for( size_t i = 0; i < channelDataPointers.size(); ++i ) {
        //	std::cout << (int)channelDataPointers[i] << std::endl;
        //	std::cout << (int)channelData[i].get() << std::endl;
        // }

        // A data arrays for the level set
        std::vector<float> rlsDistanceData;
        std::vector<std::vector<char>> rlsChannelData( channelName.size() );

        // Iterate through the z voxel coordinates in the isp vcs and build up the level set index
        // spec.  I will do this once, and build the distance data.  Then I'll loop over the channels
        // to build each of the channel data arrays.
        for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {

            // std::cout << "c: " << c << std::endl;

            // get the sparse voxel corner densities for the plane
            isp.fill_sparse_channel_data( xyExtents, c + ris.m_abcCoordOrigin.z, channelName, channelDataPointers,
                                          rlp );

            logger.update_progress( c, ris.m_abcCoordSize.zsize() );

            // rlp.check_consistency(std::cout);
            // rlp.dump(std::cout);

            // populate the index spec using the rle plane generated
            for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
                // std::cout << "  b: " << b << std::endl;
                //  calculate the offset for the run indices
                int offset = b * ris.m_abcCoordSize.xsize() - ris.m_abcCoordOrigin.x;

                // for each b (y) extent, put the current runindex in the array
                int bcIndex = b + c * ris.m_abcCoordSize.ysize();

                // Indicate where in the m_runData array the runs for this BC coordinate begin
                ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

                // Assign the run data
                int rlpRunIndex = rlp.get_y_extent_first_run( b );

                // If the extent has no defined runs
                if( rlpRunIndex < 0 ) {
                    // std::cout << "    no defined runs" << std::endl;
                    //  It's one of two cases.  It's either an empty scanline, or a full undefined run.
                    //  If the data matches the exterior region code, empty run, otherwise, full undefined.
                    if( rlpRunIndex == ris.m_exteriorRegionCode ) {
                        // std::cout << "    inserting empty extent" << std::endl;
                        ris.m_runData.push_back(
                            frantic::volumetrics::levelset::run_data( 0, ris.m_exteriorRegionCode ) );
                        ris.m_runData.push_back(
                            frantic::volumetrics::levelset::run_data( 0, ris.m_exteriorRegionCode ) );
                    } else {
                        // std::cout << "    inserting full undefined run" << std::endl;
                        ris.m_runData.push_back(
                            frantic::volumetrics::levelset::run_data( ris.m_abcCoordOrigin.x, rlpRunIndex ) );
                        ris.m_runData.push_back( frantic::volumetrics::levelset::run_data(
                            ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize(), rlpRunIndex ) );
                    }

                } else {
                    // std::cout << "first run of extent: " << rlp[rlpRunIndex].first << "," << rlp[rlpRunIndex].second
                    // << std::endl;

                    // Index data for this extent
                    int lastRunOfExtent = rlp.get_y_extent_last_run( b );
                    // std::cout << "last run of extent: " << rlp[lastRunOfExtent].first << "," <<
                    // rlp[lastRunOfExtent].second << std::endl;

                    // loop through and add all the run data.
                    while( rlpRunIndex <= lastRunOfExtent ) {
                        // std::cout << "  currentrun: " << rlp[rlpRunIndex].first << "," << rlp[rlpRunIndex].second <<
                        // std::endl;

                        if( rlp.get_run_code( rlpRunIndex ) < 0 ) {
                            ris.m_runData.push_back( frantic::volumetrics::levelset::run_data(
                                rlp[rlpRunIndex].first - offset, rlp.get_run_code( rlpRunIndex ) ) );
                        } else {

                            // index data
                            int runStart = rlp[rlpRunIndex].first - offset;
                            // int runEnd = rlp[rlpRunIndex].second - offset;
                            int currentSize = (int)rlsDistanceData.size();

                            // std::cout << "    runStart: " << runStart << "  runEnd: " << runEnd << "  currentSize: "
                            // << currentSize
                            // << std::endl;

                            // add the defined run index data
                            ris.m_runData.push_back(
                                frantic::volumetrics::levelset::run_data( runStart, currentSize ) );
                            // std::cout << "    ris.m_runData.size(): " << ris.m_runData.size() << std::endl;

                            // copy in the actual distance data for the defined run
                            int runLength = rlp[rlpRunIndex].second - rlp[rlpRunIndex].first + 1;
                            rlsDistanceData.resize( currentSize + runLength );

                            // std::cout << "runLength: " << runLength << std::endl;
                            // memcpy(&rlsDistanceData[currentSize], &distanceData[rlp[rlpRunIndex].first],
                            // sizeof(float)*runLength);
                            memcpy( &rlsDistanceData[currentSize],
                                    &channelData[signedDistanceChannel][rlp[rlpRunIndex].first * sizeof( float )],
                                    sizeof( float ) * runLength );
                            // update the data size to reflect the allocation of the new defined run
                            ris.m_dataSize += runLength;
                        }
                        ++rlpRunIndex;
                        // for ( int i = 0; i < ris.m_runData.size() ;  ++i)
                        // std::cout << "x: " << ris.m_runData[i].x << "  dataIndex: " << ris.m_runData[i].dataIndex <<
                        // std::endl;
                    }

                    int terminalRunIndex = rlp[rlpRunIndex - 1].second + 1 - offset;

                    // Add in the terminal run at the end to "close off the scanline",
                    // doesnt matter what's in it.
                    ris.m_runData.push_back(
                        frantic::volumetrics::levelset::run_data( terminalRunIndex, ris.m_exteriorRegionCode ) );
                    // std::cout << "  terminating run  " << terminalRunIndex << "," << ris.m_exteriorRegionCode <<
                    // std::endl; std::cout << "  ris.m_runData.size(): " << ris.m_runData.size() << std::endl;
                }
            }

            // Use the rlp to go through the channel wise data and add it to the existing
            // channel data arrays.
            // std::cout << "ris.m_dataSize: " << ris.m_dataSize << std::endl;
            for( int channelNum = 0; channelNum < (int)( channelName.size() ); ++channelNum ) {
                // std::cout << c << " " << channelName[channelNum] << std::endl;

                // skip the channel if it is signed distance
                if( channelNum == signedDistanceChannel )
                    continue;

                size_t currentSize = rlsChannelData[channelNum].size();

                // if there was no new data, just skip
                if( currentSize != ris.m_dataSize * channelPrimitiveSize[channelNum] ) {
                    // std::cout << "currentSize: " << float(currentSize)/channelPrimitiveSize[channelNum] << std::endl;
                    rlsChannelData[channelNum].resize( ris.m_dataSize * channelPrimitiveSize[channelNum] );
                    // std::cout << "newSize: " <<
                    // float(rlsChannelData[channelNum].size())/channelPrimitiveSize[channelNum] << std::endl; std::cout
                    // << "processing channel: " << channelName[channelNum] << std::endl;
                    char* dest = &( rlsChannelData[channelNum][currentSize] );
                    char* source = channelData[channelNum].get();
                    // std::cout << "dest: " << (int)dest << "  source: " << (int)source << std::endl;
                    //  iterate through the rlp and copy data over
                    for( int runIndex = 0; runIndex < (int)rlp.size(); ++runIndex ) {
                        if( rlp.get_run_code( runIndex ) >= 0 ) {
                            size_t dataLength =
                                channelPrimitiveSize[channelNum] * ( rlp[runIndex].second - rlp[runIndex].first + 1 );
                            memcpy( dest, source + rlp[runIndex].first * channelPrimitiveSize[channelNum], dataLength );
                            // std::cout << ((vector3f*)(source +
                            // rlp[runIndex].first*channelPrimitiveSize[channelNum]))->str();
                            dest += dataLength;
                        }
                    }
                }
            }
        }

        // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of
        // the last scanline
        ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();

        // std::cout << "rundata: " << ris.m_runData[63].dataIndex << "  " << ris.m_runData[63].x << std::endl;
        // std::cout << "data: " << distanceData[ris.m_runData[63].dataIndex] << std::endl;

        // for ( int i = 0; i < ris.m_runData.size() ;  ++i)
        //	std::cout << "x: " << ris.m_runData[i].x << "  dataIndex: " << ris.m_runData[i].dataIndex << std::endl;
        // ris.dump(std::cout);
        // if ( !ris.check_consistency(std::cout) )
        //	throw std::runtime_error("build_rle_level_set() - Level set failed consistency check.");

        // Build a temp level set using the retrieved data
        frantic::volumetrics::levelset::rle_level_set result( isp.get_voxel_coord_system(), ris, rlsDistanceData,
                                                              interfaceVoxelWidthInside, interfaceVoxelWidthOutside );
        // std::cout << "level set size: " << result.size() << std::endl;
        //  Add in the channel data, unless there isn't any
        for( size_t i = 1; i < channelName.size(); ++i ) {
            // std::cout << channelName[i] << " channel size: " << rlsChannelData[i].size() << std::endl;
            result.add_channel( channelName[i], channelArity[i], channelType[i] );
            // if there is data to copy, copy it
            if( result.size() != 0 ) {
                frantic::volumetrics::levelset::rle_channel_general_accessor gca =
                    result.get_channel_general_accessor( channelName[i] );
                memcpy( gca.data( 0 ), &( rlsChannelData[i][0] ), channelPrimitiveSize[i] * ris.data_size() );
            }
        }
        // Swap it into the output level set
        outLevelSet.swap( result );

        // for( size_t i = 0; i < channelDataPointers.size(); ++i ){
        //	std::cout << (int)channelDataPointers[i] << std::endl;
        //	std::cout << (int)channelData[i].get() << std::endl;
        // }
    }
};
} // namespace detail

template <class ImplicitSurfacePolicy>
inline void build_rle_level_set( const frantic::channels::channel_propagation_policy& cpp, ImplicitSurfacePolicy& isp,
                                 frantic::volumetrics::levelset::rle_level_set& outLevelSet,
                                 frantic::logging::progress_logger& logger ) {
    detail::brls_ris_friend::build_rle_level_set( cpp, isp, outLevelSet, logger );
}

template <class ImplicitSurfacePolicy>
inline void build_rle_level_set( ImplicitSurfacePolicy& isp, frantic::volumetrics::levelset::rle_level_set& outLevelSet,
                                 frantic::logging::progress_logger& logger ) {
    detail::brls_ris_friend::build_rle_level_set( frantic::channels::channel_propagation_policy(), isp, outLevelSet,
                                                  logger );
}

template <class ImplicitSurfacePolicy>
inline void build_rle_level_set( ImplicitSurfacePolicy& isp,
                                 frantic::volumetrics::levelset::rle_level_set& outLevelSet ) {
    frantic::channels::channel_propagation_policy cpp;
    frantic::logging::null_progress_logger logger;
    build_rle_level_set( cpp, isp, outLevelSet, logger );
}

/**
 * Converts a set of particles into a level set using a union of spheres implicit
 * surface.  Ensure that the voxel length of the provided particle_grid_tree
 * is at least equal to the particle effect radius.  The conversion is done sparsely.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as the maximum particle's effect radius.
 * @param  maximumParticleRadius This is as big or larger than the radius of the largest particle in the system.
 * @param  particleRadiusToEffectRadiusScale	This is multiplied with each particle radius to determine the range of
 * it's effect.
 * @param  implicitThreshold     The isosurface value which defines the implicit surface.
 * @param  levelsetVCS           The voxel coordinate system for the destination level set.
 * @param  outRLS                The output parameter where the resulting mesh will go.
 * @param  logger                The progress logger for this function
 */
inline void union_of_spheres_convert_particles_to_level_set(
    frantic::particles::particle_grid_tree& particles, frantic::channels::channel_propagation_policy cpp,
    float maximumParticleRadius, float particleRadiusToEffectRadiusScale, float implicitThreshold,
    const voxel_coord_system& levelsetVCS, frantic::volumetrics::levelset::rle_level_set& outRLS,
    frantic::logging::progress_logger& logger ) {
    frantic::volumetrics::implicitsurface::particle_union_of_spheres_is_policy isp(
        particles, maximumParticleRadius, particleRadiusToEffectRadiusScale, implicitThreshold, levelsetVCS, 0 );

    frantic::volumetrics::implicitsurface::build_rle_level_set( cpp, isp, outRLS, logger );
};

inline void union_of_spheres_convert_particles_to_level_set(
    frantic::particles::particle_grid_tree& particles, float maximumParticleRadius,
    float particleRadiusToEffectRadiusScale, float implicitThreshold, const voxel_coord_system& levelsetVCS,
    frantic::volumetrics::levelset::rle_level_set& outRLS, frantic::logging::progress_logger& logger ) {
    frantic::volumetrics::implicitsurface::particle_union_of_spheres_is_policy isp(
        particles, maximumParticleRadius, particleRadiusToEffectRadiusScale, implicitThreshold, levelsetVCS, 0 );

    frantic::channels::channel_propagation_policy cpp;
    frantic::volumetrics::implicitsurface::build_rle_level_set( cpp, isp, outRLS, logger );
};

/**
 * Converts a set of particles into a level set using a metaballs implicit
 * surface.  Ensure that the voxel length of the provided particle_grid_tree
 * is at least equal to the particle effect radius.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as the maximum particle's effect radius.
 * @param  maximumParticleRadius Max radius of the particles in the above grid tree.
 * @param  particleRadiusToEffectRadiusScale	This is multiplied with each particle radius to determine the range of
 * it's effect.
 * @param  implicitThreshold     The isosurface value which defines the implicit surface.
 * @param  levelsetVCS           The voxel coordinate system for the destination level set.
 * @param  outRLS                The output parameter where the resulting level set will go.
 * @param  logger                Records the progress made by this function
 */
inline void metaball_convert_particles_to_level_set( frantic::particles::particle_grid_tree& particles,
                                                     frantic::channels::channel_propagation_policy& cpp,
                                                     float maximumParticleRadius,
                                                     float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                     const voxel_coord_system& levelsetVCS,
                                                     frantic::volumetrics::levelset::rle_level_set& outRLS,
                                                     frantic::logging::progress_logger& logger ) {
    if( particles.get_channel_map().has_channel( _T("SignedDistance") ) )
        throw std::runtime_error(
            "metaball_convert_particles_to_level_set() - The particle_grid_tree may not contain a "
            "SignedDistance channel " );
    frantic::volumetrics::implicitsurface::particle_metaball_is_policy isp(
        particles, maximumParticleRadius, particleRadiusToEffectRadiusScale, implicitThreshold, levelsetVCS, 0 );

    frantic::volumetrics::implicitsurface::build_rle_level_set( cpp, isp, outRLS, logger );
};

/**
 * Converts a set of particles into a level set using the implicit
 * surface described in the paper by Zhu and Bridson, "Animating Sand as a Fluid" from SIGGRAPH 2005.
 * Ensure that the voxel length of the provided particle_grid_tree is at least equal to the kernel compact support.
 * The radius of the particle must be stored in a float "Radius" channel.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as particleEffectRadius.
 * @param  cpp		A policy that determines which channels should be propagated from the particles
 *				to the output.
 * @param  maximumParticleRadius The max radius of particles in the system.
 * @param  effectRadius			 This is multiplied with each particle radius to determine the range of it's
 *effect.
 * @param  lowDensityTrimmingDensity    An attempt to trim away low density regions which would still be inside the
 *surface due to the formulation.  This parameter defines the density at which trimming begins.
 * @param  lowDensityTrimmingStrength   This is a multiplier which scales the magnitude of the low density trimming.
 * @param  levelsetVCS           The voxel coordinate system for the destination level set.
 * @param  outRLS                The output parameter where the resulting mesh will go.
 */
inline void zhu_bridson_convert_particles_to_level_set(
    frantic::particles::particle_grid_tree& particles, frantic::channels::channel_propagation_policy& cpp,
    float maximumParticleRadius, float effectRadius, float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
    const voxel_coord_system& levelsetVCS, frantic::volumetrics::levelset::rle_level_set& outRLS,
    frantic::logging::progress_logger& logger ) {
    frantic::volumetrics::implicitsurface::particle_zhu_bridson_is_policy isp(
        particles, maximumParticleRadius, effectRadius, lowDensityTrimmingDensity, lowDensityTrimmingStrength,
        levelsetVCS, 0 );

    frantic::volumetrics::implicitsurface::build_rle_level_set( cpp, isp, outRLS, logger );
};

inline void zhu_bridson_convert_particles_to_level_set(
    frantic::particles::particle_grid_tree& particles, float maximumParticleRadius, float effectRadius,
    float lowDensityTrimmingDensity, float lowDensityTrimmingStrength, const voxel_coord_system& levelsetVCS,
    frantic::volumetrics::levelset::rle_level_set& outRLS, frantic::logging::progress_logger& logger ) {
    frantic::volumetrics::implicitsurface::particle_zhu_bridson_is_policy isp(
        particles, maximumParticleRadius, effectRadius, lowDensityTrimmingDensity, lowDensityTrimmingStrength,
        levelsetVCS, 0 );
    frantic::channels::channel_propagation_policy cpp;
    frantic::volumetrics::implicitsurface::build_rle_level_set( cpp, isp, outRLS, logger );
};
} // namespace implicitsurface
} // namespace volumetrics
}; // namespace frantic
