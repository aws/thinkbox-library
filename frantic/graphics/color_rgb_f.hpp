// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "frantic/channels/channel_map.hpp"
#include "frantic/channels/named_channel_data.hpp"
#include "frantic/graphics/color_base.hpp"

namespace frantic {
namespace graphics {

/**
 * A pixel with a Red, Green, Blue channel
 *
 * @author   Brian McKinnon
 * @since    Apr 25, 2007
 */
class color_rgb_f : public color_base<color_rgb_f, float, 3, 0> {
  public:
    color_rgb_f( float c = 0.0f ) {
        set_r( c );
        set_g( c );
        set_b( c );
    }
    color_rgb_f( float r, float g, float b ) {
        set_r( r );
        set_g( g );
        set_b( b );
    }

    float get_r() const { return this->get_color().get_data( 0 ); }
    float get_g() const { return this->get_color().get_data( 1 ); }
    float get_b() const { return this->get_color().get_data( 2 ); }

    void set_r( float r ) { this->get_color().set_data( 0, r ); }
    void set_g( float g ) { this->get_color().set_data( 1, g ); }
    void set_b( float b ) { this->get_color().set_data( 2, b ); }

    static frantic::channels::channel_map get_channel_map();
    static void add_to_channel_map( frantic::channels::channel_map* channelMap );
    static frantic::channels::data_type_t get_data_type() { return frantic::channels::data_type_float32; }
    static frantic::tstring get_name() { return _T("R_G_B"); }
};

inline std::ostream& operator<<( std::ostream& out, const color_rgb_f& p ) {
    out << "color_rgb_f( " << p.get_r() << ", " << p.get_g() << ", " << p.get_b() << " )";
    return out;
}

inline color_rgb_f operator*( float s, const color_rgb_f& c ) {
    return color_rgb_f( s * c.get_r(), s * c.get_g(), s * c.get_b() );
}

inline color_rgb_f operator*( const color_rgb_f& c, float s ) {
    return color_rgb_f( s * c.get_r(), s * c.get_g(), s * c.get_b() );
}

} // namespace graphics
} // namespace frantic

namespace frantic {
namespace channels {

template <>
struct channel_data_type_traits<frantic::graphics::color_rgb_f> {
    typedef frantic::graphics::color_rgb_f value_type;
    inline static size_t arity() { return 3; }
    inline static data_type_t data_type() { return data_type_float32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
};

} // namespace channels
} // namespace frantic
