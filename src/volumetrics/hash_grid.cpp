// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/hash_grid.hpp>

using namespace frantic::volumetrics;

hash_grid::hash_grid( const frantic::channels::channel_map& channels )
    : m_channelMap( channels ) {}

void hash_grid::reset( const frantic::channels::channel_map& channels ) {
    m_channelMap = channels;
    m_hashTable.clear();
}

const frantic::channels::channel_map& hash_grid::get_channel_map() const { return m_channelMap; }

char* hash_grid::get_or_make_voxel_ptr( const voxel_coord& vc ) {
    const int LEAFSIZE = 5;

    voxel_coord quotient = frantic::graphics::size3( LEAFSIZE ).get_quotient_coord( vc );
    voxel_coord kernel =
        voxel_coord( vc.x - LEAFSIZE * quotient.x, vc.y - LEAFSIZE * quotient.y, vc.z - LEAFSIZE * quotient.z );

    hashtable_type::mapped_type& voxels = m_hashTable[quotient];
    if( voxels.get_data() == NULL )
        voxels.reset( ( LEAFSIZE * LEAFSIZE * LEAFSIZE ) * m_channelMap.structure_size() );

    return voxels.get_data() +
           ( kernel.x + LEAFSIZE * ( kernel.y + LEAFSIZE * kernel.z ) ) * m_channelMap.structure_size();
}

void hash_grid::get_or_make_voxel_ptrs( const voxel_coord& vc, int filterWidth, char* outPtrs[] ) {
    const int LEAFSIZE = 5;

    voxel_coord curVC = vc;

    while( curVC.z < vc.z + filterWidth ) {
        voxel_coord quotient = frantic::graphics::size3( LEAFSIZE ).get_quotient_coord( curVC );
        voxel_coord kernel = voxel_coord( curVC.x - LEAFSIZE * quotient.x, curVC.y - LEAFSIZE * quotient.y,
                                          curVC.z - LEAFSIZE * quotient.z );

        hashtable_type::mapped_type& voxels = m_hashTable[quotient];
        if( voxels.get_data() == NULL )
            voxels.reset( ( LEAFSIZE * LEAFSIZE * LEAFSIZE ) * m_channelMap.structure_size() );

        // Calculate the offsets in this leaf node, and in the output array.
        frantic::graphics::vector3 off = curVC - vc;

        // Calculate the index in the output array, and into the leaf's data indices.
        unsigned outIndex = off.x + ( filterWidth * ( off.y + off.z * filterWidth ) );
        unsigned dataIndex = kernel.x + ( LEAFSIZE * ( kernel.y + kernel.z * LEAFSIZE ) );

        // Calculate the portion of this leaf overlapped by the filter.
        unsigned xWidth = std::min( LEAFSIZE - kernel.x, filterWidth - off.x );
        unsigned yWidth = std::min( LEAFSIZE - kernel.y, filterWidth - off.y );
        unsigned zWidth = std::min( LEAFSIZE - kernel.z, filterWidth - off.z );

        // Run through the part of the filter that fits into this leaf node, filling in the
        // appropriate output elements.
        for( unsigned s = 0; s < zWidth; ++s ) {
            for( unsigned r = 0; r < yWidth; ++r ) {
                for( unsigned c = 0; c < xWidth; ++c )
                    outPtrs[outIndex + filterWidth * r + c] =
                        voxels.get_data() + ( dataIndex + LEAFSIZE * r + c ) * m_channelMap.structure_size();
            }

            dataIndex += LEAFSIZE * LEAFSIZE;
            outIndex += filterWidth * filterWidth;
        }

        // Try the next leaf in the X direction. If it is beyond our filter width, then
        // try the next leaf in the Y direction. If it is beyond our filter's extent, then
        // try the next lead in the Z direction. If it is beyond our filter's extent, then we are done.
        curVC.x += xWidth;
        if( curVC.x >= vc.x + filterWidth ) {
            curVC.x = vc.x;
            curVC.y += yWidth;
            if( curVC.y >= vc.y + filterWidth ) {
                curVC.y = vc.y;
                curVC.z += zWidth;
            }
        }
    }
}

const char* hash_grid::get_voxel_ptr( const voxel_coord& vc ) const {
    const int LEAFSIZE = 5;

    voxel_coord quotient = frantic::graphics::size3( LEAFSIZE ).get_quotient_coord( vc );
    voxel_coord kernel =
        voxel_coord( vc.x - LEAFSIZE * quotient.x, vc.y - LEAFSIZE * quotient.y, vc.z - LEAFSIZE * quotient.z );

    hashtable_type::const_iterator it = m_hashTable.find( quotient );
    if( it == m_hashTable.end() )
        return NULL;

    return it->second.get_data() +
           ( kernel.x + LEAFSIZE * ( kernel.y + LEAFSIZE * kernel.z ) ) * m_channelMap.structure_size();
}

void hash_grid::get_voxel_ptrs( const voxel_coord& vc, int filterWidth, const char* outPtrs[] ) const {
    const int LEAFSIZE = 5;

    voxel_coord curVC = vc;

    while( curVC.z < vc.z + filterWidth ) {
        voxel_coord quotient = frantic::graphics::size3( LEAFSIZE ).get_quotient_coord( curVC );
        voxel_coord kernel = voxel_coord( curVC.x - LEAFSIZE * quotient.x, curVC.y - LEAFSIZE * quotient.y,
                                          curVC.z - LEAFSIZE * quotient.z );

        hashtable_type::const_iterator it = m_hashTable.find( quotient );
        const hashtable_type::mapped_type* pData = ( it != m_hashTable.end() ) ? &it->second : NULL;

        // Calculate the offsets in this leaf node, and in the output array.
        frantic::graphics::vector3 off = curVC - vc;

        // Calculate the index in the output array, and into the leaf's data indices.
        unsigned outIndex = off.x + ( filterWidth * ( off.y + off.z * filterWidth ) );
        unsigned dataIndex = kernel.x + ( LEAFSIZE * ( kernel.y + kernel.z * LEAFSIZE ) );

        // Calculate the portion of this leaf overlapped by the filter.
        unsigned xWidth = std::min( LEAFSIZE - kernel.x, filterWidth - off.x );
        unsigned yWidth = std::min( LEAFSIZE - kernel.y, filterWidth - off.y );
        unsigned zWidth = std::min( LEAFSIZE - kernel.z, filterWidth - off.z );

        // Run through the part of the filter that fits into this leaf node, filling in the
        // appropriate output elements.
        for( unsigned s = 0; s < zWidth; ++s ) {
            for( unsigned r = 0; r < yWidth; ++r ) {
                for( unsigned c = 0; c < xWidth; ++c )
                    outPtrs[outIndex + filterWidth * r + c] =
                        pData ? pData->get_data() + ( dataIndex + LEAFSIZE * r + c ) * m_channelMap.structure_size()
                              : NULL;
            }

            dataIndex += LEAFSIZE * LEAFSIZE;
            outIndex += filterWidth * filterWidth;
        }

        // Try the next leaf in the X direction. If it is beyond our filter width, then
        // try the next leaf in the Y direction. If it is beyond our filter's extent, then
        // try the next lead in the Z direction. If it is beyond our filter's extent, then we are done.
        curVC.x += xWidth;
        if( curVC.x >= vc.x + filterWidth ) {
            curVC.x = vc.x;
            curVC.y += yWidth;
            if( curVC.y >= vc.y + filterWidth ) {
                curVC.y = vc.y;
                curVC.z += zWidth;
            }
        }
    }
}