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
 * A pixel with a Red, Green, Blue, and Alpha channel
 *
 * @author   Brian McKinnon
 * @since    Apr 25, 2007
 */
class color_rgba_f : public color_base<color_rgba_f, float, 3, 1> {
  public:
    color_rgba_f( float c = 0.0f, float a = 0.0f ) {
        set_r( c );
        set_g( c );
        set_b( c );
        set_a( a );
    }
    color_rgba_f( float r, float g, float b, float a = 0.0f ) {
        set_r( r );
        set_g( g );
        set_b( b );
        set_a( a );
    }

    float get_r() const { return this->get_color().get_data( 0 ); }
    float get_g() const { return this->get_color().get_data( 1 ); }
    float get_b() const { return this->get_color().get_data( 2 ); }
    float get_a() const { return this->get_alpha().get_data( 0 ); }

    void set_r( float r ) { this->get_color().set_data( 0, r ); }
    void set_g( float g ) { this->get_color().set_data( 1, g ); }
    void set_b( float b ) { this->get_color().set_data( 2, b ); }
    void set_a( float a ) { this->get_alpha().set_data( 0, a ); }

    color_rgba_f operator*( float scale ) const {
        return color_rgba_f( get_r() * scale, get_g() * scale, get_b() * scale, get_a() * scale );
    }

    color_rgba_f operator+( const color_rgba_f& rhs ) const {
        return color_rgba_f( get_r() + rhs.get_r(), get_g() + rhs.get_g(), get_b() + rhs.get_b(),
                             get_a() + rhs.get_a() );
    }

    static channel_map get_channel_map( frantic::tstring name = _T("R_G_B_A") );
    static void add_to_channel_map( channel_map* channelMap, frantic::tstring name = _T("R_G_B_A") );
    static data_type_t get_data_type() { return frantic::channels::data_type_float32; }
    static frantic::tstring get_name() { return _T("R_G_B_A"); }
};

inline color_rgba_f operator*( float lhs, const color_rgba_f& rhs ) { return rhs * lhs; }

inline std::ostream& operator<<( std::ostream& out, const color_rgba_f& p ) {
    out << "color_rgba_f( " << p.get_r() << ", " << p.get_g() << ", " << p.get_b() << ", " << p.get_a() << " )";
    return out;
}

} // namespace graphics
} // namespace frantic

namespace frantic {
namespace channels {

template <>
struct channel_data_type_traits<frantic::graphics::color_rgba_f> {
    typedef frantic::graphics::color_rgba_f value_type;
    inline static size_t arity() { return 4; }
    inline static data_type_t data_type() { return data_type_float32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords,
                                                  const value_type& a, const value_type& b, const value_type& c ) {
        return barycentricCoords[0] * a + barycentricCoords[1] * b + barycentricCoords[2] * c;
    }
};

} // namespace channels
} // namespace frantic
