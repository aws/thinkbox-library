// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/voxel_grid_tree.hpp>

using namespace frantic::volumetrics;

voxel_grid_tree_leaf_data::voxel_grid_tree_leaf_data( std::size_t dataSize ) {
    m_data.reset( new char[dataSize] );
    memset( m_data.get(), 0, dataSize );
}

void voxel_grid_tree_leaf_data::reset( std::size_t dataSize ) {
    m_data.reset( new char[dataSize] );
    memset( m_data.get(), 0, dataSize );
}

char* voxel_grid_tree_leaf_data::get_data() { return m_data.get(); }

const char* voxel_grid_tree_leaf_data::get_data() const { return m_data.get(); }

voxel_grid_tree::voxel_grid_tree( const frantic::channels::channel_map& channels ) { m_channelMap = channels; }

void voxel_grid_tree::reset( const frantic::channels::channel_map& channels ) {
    clear();
    m_channelMap = channels;
}

const frantic::channels::channel_map& voxel_grid_tree::get_channel_map() const { return m_channelMap; }

char* voxel_grid_tree::get_or_make_voxel_ptr( const voxel_coord& vc ) {
    const int LEAFSIZE = 5;

    voxel_coord quotient = frantic::graphics::size3( LEAFSIZE ).get_quotient_coord( vc );
    voxel_coord kernel =
        voxel_coord( vc.x - LEAFSIZE * quotient.x, vc.y - LEAFSIZE * quotient.y, vc.z - LEAFSIZE * quotient.z );

    node_type* node = navigate_to_leaf_with_create( quotient );
    if( node->leafData == NULL )
        node->leafData =
            new voxel_grid_tree_leaf_data( ( LEAFSIZE * LEAFSIZE * LEAFSIZE ) * m_channelMap.structure_size() );

    return node->leafData->get_data() +
           ( kernel.x + LEAFSIZE * ( kernel.y + LEAFSIZE * kernel.z ) ) * m_channelMap.structure_size();
}

const char* voxel_grid_tree::get_voxel_ptr( const voxel_coord& vc ) const {
    const int LEAFSIZE = 5;

    voxel_coord quotient = frantic::graphics::size3( LEAFSIZE ).get_quotient_coord( vc );
    voxel_coord kernel =
        voxel_coord( vc.x - LEAFSIZE * quotient.x, vc.y - LEAFSIZE * quotient.y, vc.z - LEAFSIZE * quotient.z );

    node_type* node = navigate_to_leaf( quotient );
    if( !node || node->leafData == NULL )
        return NULL;

    return node->leafData->get_data() +
           ( kernel.x + LEAFSIZE * ( kernel.y + LEAFSIZE * kernel.z ) ) * m_channelMap.structure_size();
}
