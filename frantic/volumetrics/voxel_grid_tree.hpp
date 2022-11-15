// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/shared_array.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/volumetrics/grid_tree_base.hpp>

namespace frantic {
namespace volumetrics {

const int VGT_SL = 3;

class voxel_grid_tree_leaf_data {
    boost::shared_array<char> m_data;

  public:
    voxel_grid_tree_leaf_data() {}

    voxel_grid_tree_leaf_data( std::size_t dataSize );

    void reset( std::size_t dataSize );

    char* get_data();

    const char* get_data() const;
};

class voxel_grid_tree : public frantic::volumetrics::grid_tree_base<voxel_grid_tree_leaf_data, VGT_SL> {
    frantic::channels::channel_map m_channelMap;

  private:
    typedef frantic::volumetrics::grid_tree_base_node<voxel_grid_tree_leaf_data, VGT_SL> node_type;

  public:
    typedef frantic::graphics::vector3 voxel_coord;

  public:
    voxel_grid_tree() {}

    voxel_grid_tree( const frantic::channels::channel_map& channels );

    void reset( const frantic::channels::channel_map& channels );

    const frantic::channels::channel_map& get_channel_map() const;

    char* get_or_make_voxel_ptr( const voxel_coord& vc );

    const char* get_voxel_ptr( const voxel_coord& vc ) const;
};

} // namespace volumetrics
} // namespace frantic