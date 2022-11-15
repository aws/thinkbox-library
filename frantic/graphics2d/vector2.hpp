// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#if defined( __GNUC__ ) && __GNUC__ < 3
#include <limits.h>
#else
#include <limits>
#endif
#include <stdexcept>
#include <string>

#include <frantic/graphics2d/vector2f.hpp>
#include <frantic/graphics2d/vector2t.hpp>

namespace frantic {
namespace graphics2d {

class vector2 : public vector2t<int, vector2> {
  public:
    vector2() {
        x = 0;
        y = 0;
    }

    vector2( int X, int Y ) {
        x = X;
        y = Y;
    }

    explicit vector2( int X ) {
        x = X;
        y = X;
    }

    static vector2 from_axis( int axis ) {
        switch( axis ) {
        case 0:
            return vector2( 1, 0 );
        case 1:
            return vector2( 0, 1 );
        default:
            return vector2();
        }
    }

    static vector2 floor( const vector2f& vec );
    static vector2 round( const vector2f& vec );

    static vector2 from_neighbor_index( int index ) {
        switch( index ) {
        case 0:
            return vector2( 1, 0 );
        case 1:
            return vector2( -1, 0 );
        case 2:
            return vector2( 0, 1 );
        case 3:
            return vector2( 0, -1 );
        default:
            throw std::runtime_error( "vector2::from_neighbor_index: Requested invalid neighbor index." );
        }
    }

    // Image pixel indices and coordinates are slightly different, these two functions convert
    // between them.
    vector2f to_image_coord() const;
    static vector2 from_image_coord( vector2f v );

    static vector2 minvalue() { return vector2( ( std::numeric_limits<int>::min )() ); }

    static vector2 maxvalue() { return vector2( ( std::numeric_limits<int>::max )() ); }

    // A valid neighbor had absolute value 1 and has two 0 coordinates
    bool is_valid_neighbor() const { return ( x * x + y * y ) == 1; }

    bool is_valid_extended_neighbor() const { return x * x <= 1 && y * y <= 1; }

    int to_neighbor_index() const {
        if( !is_valid_neighbor() )
            throw std::runtime_error( std::string( "vector2.to_neighbor_index: Tried to convert vector " ) + str() +
                                      std::string( " to a neighbor index" ) );
        if( x == 1 )
            return 0;
        if( x == -1 )
            return 1;
        if( y == 1 )
            return 2;
        if( y == -1 )
            return 3;
        throw std::runtime_error(
            std::string( "vector2.to_neighbor_index: Unexpected error, should never have come here.  Vector was " ) +
            str() );
    }
};

inline vector2f vector2::to_image_coord() const { return vector2f( x + 0.5f, y + 0.5f ); }

inline static vector2 from_image_coord( vector2f v ) { return vector2( (int)floorf( v.x ), (int)floorf( v.y ) ); }

inline vector2f operator*( float k, const vector2& a ) { return vector2f( k * a.x, k * a.y ); }

inline vector2f operator*( const vector2& a, float k ) { return vector2f( a.x * k, a.y * k ); }

inline vector2 vector2::floor( const vector2f& vec ) { return vector2( (int)floorf( vec.x ), (int)floorf( vec.y ) ); }

inline vector2 vector2::round( const vector2f& vec ) {
    return vector2( (int)floorf( vec.x + 0.5f ), (int)floorf( vec.y + 0.5f ) );
}

} // namespace graphics2d
} // namespace frantic
