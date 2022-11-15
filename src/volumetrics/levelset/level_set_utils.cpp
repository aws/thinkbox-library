// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/volumetrics/levelset/level_set_utils.hpp>

using namespace frantic::graphics;

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * For a collection of levelsets, this function unions all of them together where there are any bounding box
 * intersections returns of a set of level sets that have mutually exclusive inside regions. This function allocates new
 * level sets where necessary to hold the union results.
 *
 * @param levelsets The collection of level sets to union
 * @returns unionResults a collection of level sets whose inside regions are exclusive (i.e have no intersections with
 * each other)
 */
void multiple_csg_union( std::vector<boost::shared_ptr<rle_level_set>>& levelSets,
                         std::vector<boost::shared_ptr<rle_level_set>>& unionResults ) {
    using namespace std;
    using namespace boost;

    typedef std::pair<size_t, boost::shared_ptr<rle_level_set>> levelset_pair_t;

    vector<boost::shared_ptr<rle_level_set>> seedLevelSets;
    priority_queue<levelset_pair_t, vector<levelset_pair_t>, rle_size_compare<boost::shared_ptr<rle_level_set>>>
        objsToUnion;

    voxel_coord_system vcs = levelSets[0]->get_voxel_coord_system();
    vector<boundbox3f> bboxes( levelSets.size() );
    for( unsigned obj = 0; obj < levelSets.size(); ++obj ) {
        bboxes[obj] = vcs.get_world_bounds( levelSets[obj]->get_rle_index_spec().outer_bounds() );
    }
    float voxelLength = vcs.voxel_length();
    float volumeOfVoxel = voxelLength * voxelLength * voxelLength;
    // perform any necessary unions of the levelsets
    for( unsigned obj = 0; obj < levelSets.size(); ++obj ) {
        // boundbox3f bounds(bboxes[obj]);
        boost::shared_ptr<rle_level_set> seedLS( new rle_level_set( vcs ) );

        // bounds.expand( 2*seederData.particleRadius );

        vector<int> toUnionIndices;

        // check the bound boxes to see if a union is necessary
        for( unsigned i = 0; i < levelSets.size(); ++i ) {
            if( i == obj )
                continue;

            bboxes[obj].intersect_with( bboxes[i] );

            // add based on the size of the intersected bounds
            if( bboxes[obj].volume() > volumeOfVoxel ) {
                objsToUnion.push( levelset_pair_t( levelSets[i]->size(), levelSets[i] ) );
            }
        }

        // there are some intersection objects so we have to perform at least one union
        if( objsToUnion.size() > 0 ) {
            objsToUnion.push( levelset_pair_t( levelSets[obj]->size(), levelSets[obj] ) );

            boost::shared_ptr<rle_level_set> theUnion;
            // perform any unions with the
            while( objsToUnion.size() > 1 ) {
                boost::shared_ptr<rle_level_set> first = objsToUnion.top().second;
                objsToUnion.pop();
                boost::shared_ptr<rle_level_set> second = objsToUnion.top().second;
                objsToUnion.pop();

                theUnion.reset( new rle_level_set( vcs ) );
                // fout << "Prepared Result Level Set " << endl;

                // psLevelSetUnion.enter();

                theUnion->csg_union( *first, *second );

                // psLevelSetUnion.exit();

                objsToUnion.push( levelset_pair_t( theUnion->size(), theUnion ) );

                // fout << "Done Union, queue size: " << objsToUnion.size() <<  endl;
            }

            // fout << "\tLevel Set pre-processing complete" << endl;
            // fout << "\tUnioned Region Level Set Size = " << objsToUnion.top().second->size() << endl;
            //  we know have the final level set and we can intersect it with the surface

            unionResults.push_back( objsToUnion.top().second );

        }
        // we can use the levelset directly since it doesn't intersect any objects
        else {
            unionResults.push_back( levelSets[obj] );
        }
    }
}

void self_advect_rle_level_set( rle_level_set& levelSet, float frameOffset ) {
    // Make a copy of the level set to advect into, so that we have all the required channels of appropriate size.
    rle_level_set rlsTemp = rle_level_set( levelSet );
    self_advect_rle_level_set( levelSet, rlsTemp, frameOffset );
}

void self_advect_rle_level_set( rle_level_set& levelSet, const rle_level_set& tempLevelSetCopy, float frameOffset ) {
    // advect based on velocity
    if( levelSet.has_channel( _T("Velocity") ) ) {

        // Fetch all the channel names to be advected
        std::vector<frantic::tstring> channelNames;
        levelSet.get_channel_names( channelNames );

        // Get rid of the velocity channel.  It can't be among the channels advected.
        for( size_t c = 0; c < channelNames.size(); ++c ) {
            if( channelNames[c] == _T("Velocity") ) {
                channelNames[c] = channelNames[channelNames.size() - 1];
                channelNames.pop_back();
                break;
            }
        }

        levelSet.semi_lagrangian_advect_channels( tempLevelSetCopy, frameOffset, channelNames );
        levelSet.semi_lagrangian_advect( tempLevelSetCopy, _T("Velocity"), frameOffset );
        levelSet.reinitialize_signed_distance( _T("ReinitializePopulated"),
                                               -levelSet.get_interface_voxel_width_inside(),
                                               levelSet.get_interface_voxel_width_outside() );
        levelSet.trim_to_populated( _T("ReinitializePopulated") );
        levelSet.erase_channel( _T("ReinitializePopulated") );
    }
}

void expand_rle_level_set_to_max_velocity( rle_level_set& levelSet, float timeStep ) {
    if( levelSet.has_channel( _T("Velocity") ) ) {

        // Fetch all the channel names to be advected
        std::vector<frantic::tstring> channelNames;
        levelSet.get_channel_names( channelNames );

        // resize to make sure we have the defined data for the advection
        float maxVel = levelSet.get_channel_max_norm( _T("Velocity") );

        float voxelLength = levelSet.get_voxel_coord_system().voxel_length();
        levelSet.dilate_defined_voxels( 1 + (int)( maxVel * fabsf( timeStep ) / voxelLength ), _T("DilatePopulated") );

        levelSet.duplicate_channel( _T("ExtrapolatePopulated"), _T("DilatePopulated") );
        levelSet.reinitialize_signed_distance_from_populated( _T("DilatePopulated") );
        levelSet.extrapolate_channels( channelNames, _T("ExtrapolatePopulated") );
        levelSet.erase_channel( _T("DilatePopulated") );
        levelSet.erase_channel( _T("ExtrapolatePopulated") );
    }
}

namespace detail {

bool has_channel( const frantic::fluids::rle_voxel_field& rvf, const frantic::tstring& channelName ) {
    return rvf.has_channel( channelName );
}

void add_channel( frantic::fluids::rle_voxel_field& rvf, const frantic::tstring& channelName, std::size_t arity,
                  data_type_t dataType ) {
    rvf.add_channel( channelName, arity, dataType );
}

std::pair<std::size_t, data_type_t> get_channel_type( const frantic::fluids::rle_voxel_field& rvf,
                                                      const frantic::tstring& channelName ) {
    if( !rvf.has_channel( channelName ) )
        throw std::runtime_error( "get_channel_type Error: the rle_voxel_field does not have the specified channel \'" +
                                  frantic::strings::to_string( channelName ) + "\'." );

    const_rle_channel_general_accessor acc = rvf.get_channel_general_accessor( channelName );

    return std::pair<std::size_t, data_type_t>( acc.arity(), acc.data_type() );
}

bool has_channel( const frantic::volumetrics::levelset::rle_level_set& rls, const frantic::tstring& channelName ) {
    if( channelName == _T("SignedDistance") )
        return true;
    else
        return rls.has_channel( channelName );
}

void add_channel( frantic::volumetrics::levelset::rle_level_set& rls, const frantic::tstring& channelName,
                  std::size_t arity, data_type_t dataType ) {
    if( channelName == _T("SignedDistance") ) {
        if( arity != frantic::channels::channel_data_type_traits<float>::arity() ||
            dataType != frantic::channels::channel_data_type_traits<float>::data_type() ) {
            throw std::runtime_error(
                "add_channel Error: requested adding a \'SignedDistance\' channel of type " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str( arity, dataType ) ) +
                " to a rle_level_set, but this channel must be of type " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<float>::type_str() ) + "." );
        }
    } else {
        rls.add_channel( channelName, arity, dataType );
    }
}

std::pair<std::size_t, data_type_t> get_channel_type( const frantic::volumetrics::levelset::rle_level_set& rls,
                                                      const frantic::tstring& channelName ) {
    if( channelName == _T("SignedDistance") ) {
        return std::pair<std::size_t, data_type_t>( 1, frantic::channels::data_type_float32 );
    } else {
        if( !rls.has_channel( channelName ) )
            throw std::runtime_error(
                "get_channel_type Error: the rle_level_set does not have the specified channel \'" +
                frantic::strings::to_string( channelName ) + "\'." );

        const_rle_channel_general_accessor acc = rls.get_channel_general_accessor( channelName );

        return std::pair<std::size_t, data_type_t>( acc.arity(), acc.data_type() );
    }
}

template <>
const float* get_signed_distance_data_pointer<float>( const frantic::volumetrics::levelset::rle_level_set& rls ) {
    if( rls.size() == 0 )
        return NULL;
    return &rls[0];
}

template <>
float* get_signed_distance_data_pointer<float>( frantic::volumetrics::levelset::rle_level_set& rls ) {
    if( rls.size() == 0 )
        return NULL;
    return &rls[0];
}

} // namespace detail

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
