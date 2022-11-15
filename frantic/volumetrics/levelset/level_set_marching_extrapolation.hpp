// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_adjacency.hpp>

#include <frantic/files/filename_sequence.hpp>
#include <frantic/geometry/trimesh3_file_io.hpp>
#include <frantic/volumetrics/voxel_field_visualization.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

struct extrapolation_debug_info {
    const frantic::volumetrics::voxel_coord_system* vcs;
    const frantic::volumetrics::levelset::rle_index_spec* ris;
    const_rle_channel_general_accessor* channelAccessor;
    float meshingRes;
    int counter;

    inline void set( const frantic::volumetrics::voxel_coord_system& v,
                     const frantic::volumetrics::levelset::rle_index_spec& r, const_rle_channel_general_accessor& c,
                     float res ) {
        vcs = &v;
        ris = &r;
        channelAccessor = &c;
        counter = 0;
        meshingRes = res;
    }

    inline void reset() {
        vcs = 0;
        ris = 0;
        channelAccessor = 0;
        counter = 0;
    }
    inline void create_debug_mesh() {
        if( vcs != 0 && ris != 0 ) {
            frantic::geometry::trimesh3 debugMesh;
            volumetrics::visualization::debug_mesh_populated_color_policy policy;
            volumetrics::visualization::convert_voxel_channel_to_debug_mesh( *vcs, *ris, meshingRes, *channelAccessor,
                                                                             policy, debugMesh );
            frantic::files::filename_pattern fp( _T("C:\\temp\\extrap\\extrapolation_debug_0000.xmesh") );
            write_mesh_file( fp[counter], debugMesh );
            ++counter;
        } else
            std::cout << "COULD NOT SAVE OUT DEBUG MESH BECAUSE THINGS ARE FUCKING NULL....." << std::endl;
    }
};

/**
 * This struct is used to pass data about each channel to extrapolate
 * to the marching_extrapolation algorithm.
 */
struct marching_extrapolation_channel {
    /// This is the raw channel data to extrapolate
    char* data;
    /// If non-zero, this is the gradient of the raw channel data to try and match.  It must have the same
    /// data type, and if non-zero, the arity of gradient channel must be 3 times the arity fo the raw data
    /// channel.  For data with arity != 1, the gradient is laid out as
    /// [all d(data)/dx components], [all d(data)/dy components], [all d(data)/dz components].
    const char* gradientToMatch;
    std::size_t primitiveSize, arity;
    channels::channel_weighted_sum_combine_function_t weightedSumFn;
};

/**
 * This function implements extrapolation of level set named channels.
 *
 * This is the algorithm "extrapolation_starting_points" given in pseudo-code there.
 *
 * Note that this algorithm does both constant extrapolation and extrapolation matching a gradient, based on
 * whether the gradientToMatch pointers in the marching_extrapolation_channel structures are set.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  voxelLength  The length of one voxel.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  extrapChannels  Data about the channels that need to be extrapolated.
 * @param  db  A class which records debugging information.
 */
void marching_extrapolation( const ris_adjacency& adj, float voxelLength, unsigned char* populatedChannel,
                             const float* signedDistanceChannel,
                             const std::vector<marching_extrapolation_channel>& extrapChannels,
                             extrapolation_debug_info& db );

/**
 * This function implements extrapolation of a single staggered channel, which must have type float32[3].
 * The staggered channel is expected to be in a different level set or voxel field, and the
 * int32* dataIndexMapChannel is a data index map (created by the create_data_index_map_channel
 * function, for instance) which provides the mapping from the indices for signedDistanceChannel to indices
 * for staggeredExtrapChannel.
 *
 * @note it looks like this function leaves the staggeredPopulatedChannel
 *		unmodified, unlike the centred channel extrapolation.
 *		Was this intentional ?  I'm not sure if this is a good idea, because
 *		the user may want to know if some voxels are missed during
 *		extrapolation.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information
 *              for the signedDistanceChannel (but *NOT* for staggeredExtrapChannel!).
 * @param  signedDistanceChannel  The raw memory for the float32 'SignedDistance' channel.
 * @param  staggeredPopulatedChannel  The raw memory for the uint8 'StaggeredPopulated' channel.  This field corresponds
 *                                    to the signedDistanceChannel.  Each voxel gets 3 bits, bit 0 indicating X
 *direction, bit 1 for Y, and bit 2 for Z.  For example, a value of 5 = 0b0101 indicates that X and Z are populated, but
 *Y is not.  Any populated voxel is used as a source for extrapolation, and any unpopulated voxel is used as a
 *destination.
 * @param  extrapDirection  This indicates which way extrapolation occurs.  +1 indicates from negative signed distance
 *to positive (inside to outside), and -1 indicates the opposite.
 * @param  dataIndexMapChannel  A mapping which goes from the data indices of signedDistanceChannel to
 *                              data indices of staggeredExtrapChannel.
 * @param  staggeredExtrapChannel  The raw memory for the float32[3] staggeredExtrapChannel.
 */
void staggered_field_marching_extrapolation( const ris_adjacency& adj, const float* signedDistanceChannel,
                                             const boost::uint8_t* staggeredPopulatedChannel, int extrapDirection,
                                             boost::int32_t* dataIndexMapChannel, float* staggeredExtrapChannel );

/**
 *  This function implements extrapolation of level set named channels.
 *
 *  This is a reimplementation of marching_extrapolation using a queue-based
 * iterative algorithm.  This implementation is no doubt slower but it will
 * hopefully be easier to get working correctly.
 *
 *  This function extrapolates channel values away from the phi = 0 isosurface.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  voxelLength  The length of one voxel.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  extrapChannels  The channels to extrapolate.
 * @param  db  A class which records debugging information.
 */
void marching_extrapolation_iterative( const ris_adjacency& adj, const float voxelLength,
                                       unsigned char* populatedChannel, const float* signedDistanceChannel,
                                       const std::vector<marching_extrapolation_channel>& extrapChannels,
                                       extrapolation_debug_info& db
                                       /*std::size_t* populatedOrderDebugChannel = 0*/ );

/**
 * This function implements extrapolation of a single staggered channel, which must have type float32[3].
 * The staggered channel is expected to be in a different level set or voxel field, and the
 * int32* dataIndexMapChannel is a data index map (created by the create_data_index_map_channel
 * function, for instance) which provides the mapping from the indices for signedDistanceChannel to indices
 * for staggeredExtrapChannel.
 *
 *  This is a reimplementation of staggered_field_marching_extrapolation
 * using a queue-based iterative algorithm.  This implementation is no doubt
 * slower but it will hopefully be easier to get working correctly.
 *
 *  This function does not change the values in the staggeredPopulatedChannel.
 *
 * @note it looks like this function leaves the staggeredPopulatedChannel
 *		unmodified, unlike the centred channel extrapolation.
 *		Was this intentional ?  I'm not sure if this is a good idea, because
 *		the user may want to know if some voxels are missed during
 *		extrapolation.
 * @todo consider changing the staggeredPopulatedChannel behaviour.  I'm
 *		pretty sure we want to change it to reflect the newly populated
 *		voxels, like we do in the centred channel extrapolation.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information
 *              for the signedDistanceChannel (but *NOT* for staggeredExtrapChannel!).
 * @param  signedDistanceChannel  The raw memory for the float32 'SignedDistance' channel.
 * @param  staggeredPopulatedChannel  The raw memory for the uint8 'StaggeredPopulated' channel.  This field corresponds
 *                                    to the signedDistanceChannel.  Bits [0..2] are used to
 *									  flag whether a face is populated, bit 0 indicating
 *X direction, bit 1 for Y, and bit 2 for Z.  For example, a value of 5 = 0b0101 indicates that X and Z are populated,
 *but Y is not.  Any populated voxel is used as a source for extrapolation, and any unpopulated voxel is used as a
 *destination.
 * @param  extrapDirection  This indicates which way extrapolation occurs.  +1 indicates from negative signed distance
 *to positive (inside to outside), and -1 indicates the opposite.
 * @param  dataIndexMapChannel  A mapping which goes from the data indices of signedDistanceChannel to
 *                              data indices of staggeredExtrapChannel.
 * @param  staggeredExtrapChannel  The raw memory for the float32[3] staggeredExtrapChannel.
 */
void staggered_field_marching_extrapolation_iterative( const ris_adjacency& adj, const float* signedDistanceChannel,
                                                       const boost::uint8_t* staggeredPopulatedChannel,
                                                       int extrapDirection, boost::int32_t* dataIndexMapChannel,
                                                       float* staggeredExtrapChannel );

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
