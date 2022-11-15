// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "frantic/graphics/color_rgba_h.hpp"

using namespace std;
using namespace frantic;
using namespace channels;
using namespace graphics;

channel_map color_rgba_h::get_channel_map() {
    assert( sizeof( color_rgba_h ) == 2 * 4 );

    channel_map channelMap;

    channelMap.define_channel( get_name(), 4, frantic::channels::data_type_float16 );
    channelMap.end_channel_definition();

    return channelMap;
}

void color_rgba_h::add_to_channel_map( channel_map* channelMap ) {
    assert( sizeof( color_rgba_h ) == 2 * 4 );

    if( !channelMap->has_channel( get_name() ) )
        channelMap->define_channel( get_name(), 4, frantic::channels::data_type_float16 );
}
