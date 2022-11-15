// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/unordered_map.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/graphics/size3.hpp>
#include <frantic/volumetrics/voxel_grid_tree.hpp>

namespace frantic {
namespace graphics {
inline std::size_t hash_value( const frantic::graphics::vector3& vec ) {
    const std::size_t p1 = 73856093;
    const std::size_t p2 = 19349663;
    const std::size_t p3 = 83492791;

    return ( vec.x * p1 ) ^ ( vec.y * p2 ) ^ ( vec.z * p3 );
}
} // namespace graphics
} // namespace frantic

namespace frantic {
namespace volumetrics {

class hash_grid {
    frantic::channels::channel_map m_channelMap;

    typedef frantic::graphics::vector3 voxel_coord;
    typedef boost::unordered_map<voxel_coord, voxel_grid_tree_leaf_data> hashtable_type;

    hashtable_type m_hashTable;

  public:
    hash_grid() {}

    hash_grid( const frantic::channels::channel_map& channels );

    void reset( const frantic::channels::channel_map& channels );

    const frantic::channels::channel_map& get_channel_map() const;

    char* get_or_make_voxel_ptr( const voxel_coord& vc );

    void get_or_make_voxel_ptrs( const voxel_coord& vc, int filterWidth, char* outPtrs[] );

    const char* get_voxel_ptr( const voxel_coord& vc ) const;

    void get_voxel_ptrs( const voxel_coord& vc, int filterWidth, const char* outPtrs[] ) const;
};
} // namespace volumetrics
} // namespace frantic
