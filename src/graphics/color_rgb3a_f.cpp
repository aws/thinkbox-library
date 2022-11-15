// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "frantic/graphics/color_rgb3a_f.hpp"

using namespace std;
using namespace frantic;
using namespace channels;
using namespace graphics;

channel_map color_rgb3a_f::get_channel_map() {
    assert( sizeof( color_rgb3a_f ) == 4 * 6 );

    channel_map channelMap;

    channelMap.define_channel( get_name(), 6, frantic::channels::data_type_float32 );
    channelMap.end_channel_definition();

    return channelMap;
}

void color_rgb3a_f::add_to_channel_map( channel_map* channelMap ) {
    assert( sizeof( color_rgb3a_f ) == 4 * 6 );

    if( !channelMap->has_channel( get_name() ) )
        channelMap->define_channel( get_name(), 6, frantic::channels::data_type_float32 );
}
