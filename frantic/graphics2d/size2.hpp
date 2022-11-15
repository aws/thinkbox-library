// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics2d/size2t.hpp>
#include <frantic/graphics2d/vector2.hpp>

namespace frantic {
namespace graphics2d {

class size2 : public size2t<int, size2> {
  public:
    size2()
        : size2t<int, size2>( 0, 0 ) {}
    explicit size2( int w )
        : size2t<int, size2>( w, w ) {}
    size2( int w, int h )
        : size2t<int, size2>( w, h ) {}
    template <typename T, typename SizeType>
    explicit size2( const size2t<T, SizeType>& s )
        : size2t<int, size2>( (T)s.xsize, (T)s.ysize ) {}

    typedef vector2 vector_type;

    bool contains( const vector_type& coord ) const {
        // Do all 4 comparisons in 2 by using unsigned numbers.  If coord.x < 0, then its unsigned representation will
        // be bigger than xsize.  Note that this only works for integers.
        return (unsigned)coord.x < (unsigned)xsize && (unsigned)coord.y < (unsigned)ysize;
    }

    bool contains( int x, int y ) const {
        // Do all 4 comparisons in 2 by using unsigned numbers.  If x < 0, then its unsigned representation will
        // be bigger than xsize.  Note that this only works for integers.
        return (unsigned)x < (unsigned)xsize && (unsigned)y < (unsigned)ysize;
    }

    bool is_coord_valid( const vector_type& coord ) const { return contains( coord ); }

    bool is_index_valid( int ix, int iy ) const { return contains( ix, iy ); }

    // Gets the index corresponding to the vector
    int get_index( const vector_type& vec ) const {
        if( !contains( vec.x, vec.y ) )
            throw std::runtime_error( "size2.get_index: Tried to get the index of an out of bounds vector: " +
                                      vec.str() + ", " + this->str() );
        return vec.x + xsize * vec.y;
    }

    // Gets the index corresponding to the vector
    int get_index( int x, int y ) const {
        if( !contains( x, y ) )
            throw std::runtime_error( "size2.get_index: Tried to get the index of an out of bounds vector: " +
                                      vector2( x, y ).str() + ", " + this->str() );
        return x + xsize * y;
    }

    // Gets the vector corresponding to the index
    vector2 get_coord( int index ) const { return vector2( index % xsize, ( index / xsize ) % ysize ); }

    int area() const { return xsize * ysize; }
};

inline size2 operator*( const size2& lhs, int rhs ) { return size2( lhs.xsize * rhs, lhs.ysize * rhs ); }

inline size2 operator*( int lhs, const size2& rhs ) { return size2( rhs.xsize * lhs, rhs.ysize * lhs ); }

} // namespace graphics2d
} // namespace frantic
