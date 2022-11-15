// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "frantic/channels/channel_map.hpp"
#include "frantic/channels/named_channel_data.hpp"
#include "frantic/graphics/color_base.hpp"

namespace frantic {
namespace graphics {

// using frantic::channels::channel_map;
using frantic::channels::data_type_t;

/**
 * A pixel with a Red, Green, Blue, and an Alpha channel for each colour
 *
 * @author   Brian McKinnon
 * @since    Apr 25, 2007
 */
class color_rgb3a_f : public color_base<color_rgb3a_f, float, 3, 3> {
  public:
    color_rgb3a_f( float c = 0.0f, float a = 1.0f ) {
        set_r( c );
        set_g( c );
        set_b( c );
        set_ar( a );
        set_ag( a );
        set_ab( a );
    }
    color_rgb3a_f( float r, float g, float b, float ar = 1.0f, float ag = 1.0f, float ab = 1.0f ) {
        set_r( r );
        set_g( g );
        set_b( b );
        set_ar( ar );
        set_ag( ag );
        set_ab( ab );
    }

    float get_r() const { return this->get_color().get_data( 0 ); }
    float get_g() const { return this->get_color().get_data( 1 ); }
    float get_b() const { return this->get_color().get_data( 2 ); }
    float get_ar() const { return this->get_alpha().get_data( 0 ); }
    float get_ag() const { return this->get_alpha().get_data( 1 ); }
    float get_ab() const { return this->get_alpha().get_data( 2 ); }

    void set_r( float r ) { this->get_color().set_data( 0, r ); }
    void set_g( float g ) { this->get_color().set_data( 1, g ); }
    void set_b( float b ) { this->get_color().set_data( 2, b ); }
    void set_ar( float ar ) { this->get_alpha().set_data( 0, ar ); }
    void set_ag( float ag ) { this->get_alpha().set_data( 1, ag ); }
    void set_ab( float ab ) { this->get_alpha().set_data( 2, ab ); }

    static channel_map get_channel_map();
    static void add_to_channel_map( channel_map* channelMap );
    static data_type_t get_data_type() { return frantic::channels::data_type_float32; }
    static frantic::tstring get_name() { return _T("R_G_B_AR_AG_AB"); }
};

inline std::ostream& operator<<( std::ostream& out, const color_rgb3a_f& p ) {
    out << "color_rgb3a_f( " << p.get_r() << ", " << p.get_g() << ", " << p.get_b() << ", " << p.get_ar() << ", "
        << p.get_ag() << ", " << p.get_ab() << " )";
    return out;
}

} // namespace graphics
} // namespace frantic
