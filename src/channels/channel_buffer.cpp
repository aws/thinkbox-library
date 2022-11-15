// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "frantic/channels/channel_buffer.hpp"

using namespace frantic;
using namespace graphics;
using namespace channels;

/**
 * Creates a channel buffer without a channel map defined
 *
 */
channel_buffer::channel_buffer() {}

/**
 * Creates a channel buffer with a custom channel map
 *
 * @param	map			The user generated channel map
 */
channel_buffer::channel_buffer( const channel_map& map ) {
    m_channelMap = map;
    m_size = size3();
}

/**
 * Creates a channel buffer with a custom channel map and
 * a defined size
 *
 * @param	map			The user generated channel map
 * @param	size			The length of the channel
 */
channel_buffer::channel_buffer( const channel_map& map, int size ) {
    m_channelMap = map;
    resize( size );
}

/**
 * Creates a channel buffer with a custom channel map and
 * a defined size
 *
 * @param	map			The user generated channel map
 * @param	size			The width and height of the image
 */
channel_buffer::channel_buffer( const channel_map& map, size2 size ) {
    m_channelMap = map;
    resize( size );
}

/**
 * Creates a channel buffer with a custom channel map and
 * a defined size
 *
 * @param	map			The user generated channel map
 * @param	size			The width, height, and depth of the image
 */
channel_buffer::channel_buffer( const channel_map& map, size3 size ) {
    m_channelMap = map;
    resize( size );
}

void channel_buffer::clear() {
    m_channelBuffer.clear();
    m_size = size3( 0, 0, 0 );
}

void channel_buffer::set_channel_map( const channel_map& map ) {
    m_channelBuffer.clear();
    m_channelMap = map;
    resize( m_size );
}

/**
 * Resizes the buffer to the defined size.  This operation
 * will release the current buffer and claim a new buffer.
 */
void channel_buffer::resize( int size ) { resize( size3( size, 1, 1 ) ); }

/**
 * Resizes the buffer to the defined size.  This operation
 * will release the current buffer and claim a new buffer.
 */
void channel_buffer::resize( size2 size ) { resize( size3( size.xsize, size.ysize, 1 ) ); }

/**
 * Resizes the buffer to the defined size.  This operation
 * will release the current buffer and claim a new buffer.
 *
 * \throw runtime_error Exception thrown when a resize is requested
 *						 without a valid channel_map
 */
void channel_buffer::resize( size3 size ) {
    if( m_channelMap.channel_definition_complete() ) {
        m_size = size;
        if( size.volume() * m_channelMap.structure_size() > 0 )
            m_channelBuffer.resize( size.volume() * m_channelMap.structure_size() );
        else
            m_channelBuffer.clear();
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Cannot resize buffer without completed channel map";
        throw std::runtime_error( strstm.str() );
    }
}
