// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "frantic/graphics/color_rgba_f.hpp"

using namespace std;
using namespace frantic;
using namespace channels;
using namespace graphics;

channel_map color_rgba_f::get_channel_map( frantic::tstring name ) {
    assert( sizeof( color_rgba_f ) == 4 * 4 );

    channel_map channelMap;

    channelMap.define_channel( name, 4, frantic::channels::data_type_float32 );
    channelMap.end_channel_definition();

    return channelMap;
}

void color_rgba_f::add_to_channel_map( channel_map* channelMap, frantic::tstring name ) {
    assert( sizeof( color_rgba_f ) == 4 * 4 );

    if( !channelMap->has_channel( name ) )
        channelMap->define_channel( name, 4, frantic::channels::data_type_float32 );
}
