// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "frantic/channels/named_channel_data.hpp"

namespace frantic {
namespace graphics {

using frantic::channels::data_type_t;
/**
 * This class handles the conversion from one pixel type to another, and is primarily used
 * by frantic::channels::generic_channel_buffer_iterator.  There are some caveats to using the
 * generic interface for pixel data that you must be aware of.
 *
 * frantic::graphics::color_rgb_f do not alter the alpha if the destination has an alpha channel.
 * Instead, the existing alpha value is premultiplied into the rgb values, and the resulting value
 * is stored.  The opposite applies an alpha division, and retains the existing alpha channel.
 *
 */
class color_converter {
  public:
    static void ( *get_conversion_function( const frantic::tstring& s, data_type_t sType, const frantic::tstring& d,
                                            data_type_t dType ) )( void*, void* );
};

} // namespace graphics
} // namespace frantic
