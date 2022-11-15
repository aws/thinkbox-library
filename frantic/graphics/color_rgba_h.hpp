// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "frantic/channels/channel_map.hpp"
#include "frantic/channels/named_channel_data.hpp"
#include "frantic/graphics/color_base.hpp"
#include <half.h>

namespace frantic {
namespace graphics {

// using frantic::channels::channel_map;
using frantic::channels::data_type_t;

/**
 * A pixel with a Red, Green, Blue, and Alpha channel store in an Imf half
 *
 * @author   Brian McKinnon
 * @since    May 2, 2007
 */
class color_rgba_h : public color_base<color_rgba_h, half, 3, 1> {
  public:
    color_rgba_h( half c = 0.0f, half a = 0.0f ) {
        set_r( c );
        set_g( c );
        set_b( c );
        set_a( a );
    }
    color_rgba_h( half r, half g, half b, half a = 0.0f ) {
        set_r( r );
        set_g( g );
        set_b( b );
        set_a( a );
    }

    half get_r() const { return this->get_color().get_data( 0 ); }
    half get_g() const { return this->get_color().get_data( 1 ); }
    half get_b() const { return this->get_color().get_data( 2 ); }
    half get_a() const { return this->get_alpha().get_data( 0 ); }

    void set_r( half r ) { this->get_color().set_data( 0, r ); }
    void set_g( half g ) { this->get_color().set_data( 1, g ); }
    void set_b( half b ) { this->get_color().set_data( 2, b ); }
    void set_a( half a ) { this->get_alpha().set_data( 0, a ); }

    static channel_map get_channel_map();
    static void add_to_channel_map( channel_map* channelMap );
    static data_type_t get_data_type() { return frantic::channels::data_type_float16; }
    static frantic::tstring get_name() { return _T("R_G_B_A"); }
};

inline std::ostream& operator<<( std::ostream& out, const color_rgba_h& p ) {
    out << "color_rgba_h( " << p.get_r() << ", " << p.get_g() << ", " << p.get_b() << ", " << p.get_a() << " )";
    return out;
}

} // namespace graphics
} // namespace frantic
