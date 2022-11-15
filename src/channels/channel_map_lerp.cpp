// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/channels/channel_map_lerp.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

using namespace std;
using namespace frantic::graphics;

namespace frantic {
namespace channels {

struct individual_element_lerp {
    // The source and destination positions, in bytes.
    unsigned position;
    data_type_t dataType;

    individual_element_lerp( unsigned position_, data_type_t dataType_ )
        : position( position_ )
        , dataType( dataType_ ) {}
};

// Prepares the class to do lerp operations for this particle channel map structure
void channel_map_lerp::set( const channel_map& channelMap ) {
    // Reset the existing values just in case
    clear();

    if( !channelMap.channel_definition_complete() )
        throw runtime_error( "channel_map_lerp.set: The channel_map definition process was not complete.  This must be "
                             "finished before a lerp accelerator can be created." );

    m_structureSize = channelMap.structure_size();

    std::vector<individual_element_lerp> elementLerps;

    // First convert the lerping into a element-by-element lerp representation
    for( unsigned channelNum = 0; channelNum < channelMap.channel_count(); ++channelNum ) {
        const channel& pc = channelMap[channelNum];

        unsigned dataSize = (unsigned)sizeof_channel_data_type( pc.data_type() );
        for( unsigned element = 0; element < pc.arity(); ++element ) {
            elementLerps.push_back(
                individual_element_lerp( unsigned( pc.offset() + element * dataSize ), pc.data_type() ) );
        }
    }

    // Then run-length encode the resulting lerps into a sequence of lerp blocks
    unsigned position = 0, lerpElementCount = 0;
    data_type_t dataType = data_type_invalid;
    std::size_t dataTypeSize = 0;
    for( unsigned i = 0; i < elementLerps.size(); ++i ) {
        // If we can combine this into an existing lerp block, do it
        if( lerpElementCount > 0 && dataType == elementLerps[i].dataType &&
            position + lerpElementCount * dataTypeSize == elementLerps[i].position ) {
            ++lerpElementCount;
        }
        // Otherwise close off any previous blocks, and start a new lerp block
        else {
            if( lerpElementCount > 0 ) {
                m_lerpBlocks.push_back( detail::cml_lerp_block( position, lerpElementCount,
                                                                channel_weighted_sum_combine_function( dataType ) ) );
                lerpElementCount = 0;
            }
            // Start a new lerp block
            position = elementLerps[i].position;
            lerpElementCount = 1;
            dataType = elementLerps[i].dataType;
            dataTypeSize = sizeof_channel_data_type( dataType );
        }
    }
    // Close off the last block
    if( lerpElementCount > 0 ) {
        m_lerpBlocks.push_back(
            detail::cml_lerp_block( position, lerpElementCount, channel_weighted_sum_combine_function( dataType ) ) );
    }
}

// This does a linear interpolation between the two source particles, into the destination particle.
void channel_map_lerp::lerp( void* destData, const void* sourceDataA, const void* sourceDataB, float t ) {
    if( t == 0 ) {
        if( destData != sourceDataA ) {
            memcpy( destData, sourceDataA, m_structureSize );
        }
    } else if( t == 1 ) {
        if( destData != sourceDataB ) {
            memcpy( destData, sourceDataB, m_structureSize );
        }
    } else {
        float weights[2];
        const char* data[2];
        weights[0] = 1 - t;
        weights[1] = t;
        for( vector<detail::cml_lerp_block>::const_iterator i = m_lerpBlocks.begin(); i != m_lerpBlocks.end(); ++i ) {
            data[0] = reinterpret_cast<const char*>( sourceDataA ) + i->position;
            data[1] = reinterpret_cast<const char*>( sourceDataB ) + i->position;
            i->lerpCombine( weights, data, 2, i->count, reinterpret_cast<char*>( destData ) + i->position );
        }
    }
}

// This does a linear interpolation of the destination particle with the source particle.
void channel_map_lerp::lerp( void* destDataA, const void* sourceDataB, float t ) {
    if( t == 1 ) {
        if( destDataA != sourceDataB ) {
            memcpy( destDataA, sourceDataB, m_structureSize );
        }
    } else if( t != 0 ) {
        float weights[2];
        const char* data[2];
        weights[0] = 1 - t;
        weights[1] = t;
        for( vector<detail::cml_lerp_block>::const_iterator i = m_lerpBlocks.begin(); i != m_lerpBlocks.end(); ++i ) {
            data[0] = reinterpret_cast<const char*>( destDataA ) + i->position;
            data[1] = reinterpret_cast<const char*>( sourceDataB ) + i->position;
            i->lerpCombine( weights, data, 2, i->count, reinterpret_cast<char*>( destDataA ) + i->position );
        }
    }
}

void channel_map_lerp::trilerp( void* _dest, const void* _src[8], float ( &offsets )[3] ) {
    float weights[8];
    const char* src[8];
    char* dest;

    frantic::volumetrics::levelset::get_trilerp_weights( offsets, weights );

    for( vector<detail::cml_lerp_block>::const_iterator i = m_lerpBlocks.begin(); i != m_lerpBlocks.end(); ++i ) {
        dest = reinterpret_cast<char*>( _dest ) + i->position;
        src[0] = reinterpret_cast<const char*>( _src[0] ) + i->position;
        src[1] = reinterpret_cast<const char*>( _src[1] ) + i->position;
        src[2] = reinterpret_cast<const char*>( _src[2] ) + i->position;
        src[3] = reinterpret_cast<const char*>( _src[3] ) + i->position;
        src[4] = reinterpret_cast<const char*>( _src[4] ) + i->position;
        src[5] = reinterpret_cast<const char*>( _src[5] ) + i->position;
        src[6] = reinterpret_cast<const char*>( _src[6] ) + i->position;
        src[7] = reinterpret_cast<const char*>( _src[7] ) + i->position;
        i->lerpCombine( weights, src, 8, i->count, dest );
    }
}

} // namespace channels
} // namespace frantic
