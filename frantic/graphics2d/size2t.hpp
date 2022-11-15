// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics2d/vector2t.hpp>
#include <frantic/misc/string_functions.hpp>
#include <stdexcept>

namespace frantic {
namespace graphics2d {

template <typename T, typename SizeType>
class size2t {
  public:
    T xsize, ysize;

    typedef T value_type;

    /**************
     * Constructors
     **************/
    size2t() {
        xsize = 0;
        ysize = 0;
    }
    size2t( T w, T h ) {
        xsize = w;
        ysize = h;
    }

    //	template<typename U>
    //	size2t( size2t<U> & s ) {
    //		xsize = (T)s.xsize; ysize = (T)s.ysize;
    //	}

    /**************
     * Mutators
     **************/
    void set_xsize( T w ) { xsize = w; }
    void set_ysize( T h ) { ysize = h; }

    /**************
     * Accessors
     **************/
    T get_xsize() const { return xsize; }
    T get_ysize() const { return ysize; }
    T get_area() const { return xsize * ysize; }

    static SizeType from_cube_side_length( T sideLength ) { return SizeType( sideLength, sideLength ); }

    static SizeType from_square_side_length( T sideLength ) { return SizeType( sideLength, sideLength ); }

    static SizeType empty() { return SizeType( 0, 0 ); }

    T get_max_dimension() const { return ( std::max )( xsize, ysize ); }
    bool is_square() const { return xsize == ysize; }

    T get_square_length() const {
        if( !is_square() ) {
            throw std::runtime_error( "size2::get_square_length: not a cube" );
        }
        return xsize;
    }

    template <class VectorType>
    VectorType clamp( const vector2t<T, VectorType>& vec ) const {
        return VectorType( ( std::max )( 0, ( std::min )( xsize - 1, vec.x ) ),
                           ( std::max )( 0, ( std::min )( ysize - 1, vec.y ) ) );
    }

    // Assuming we have a grid of blocks the size of this size2, this function returns
    // the coordinate in that higher level grid
    template <class VectorType>
    VectorType get_quotient_coord( const vector2t<T, VectorType>& vec ) const {
        // Can't just use /, because then negative numbers go wrong
        VectorType result;
        if( vec.x >= 0 )
            result.x = vec.x / xsize;
        else
            result.x = -1 - ( ( -1 - vec.x ) / xsize );
        if( vec.y >= 0 )
            result.y = vec.y / ysize;
        else
            result.y = -1 - ( ( -1 - vec.y ) / ysize );

        return result;
    }

    // Assuming we have a grid of blocks the size of this size2, this function returns
    // the coordinate in the voxel that contains the vec
    template <class VectorType>
    VectorType get_kernel_coord( const vector2t<T, VectorType> vec ) const {
        VectorType quotientCoord = get_quotient_coord( vec );
        return VectorType( vec.x - xsize * quotientCoord.x, vec.y - ysize * quotientCoord.y );
    }

    static SizeType parse( const std::string& input ) {
        std::vector<std::string> values;
        frantic::strings::split( input, values, "()[], " );
        if( values.size() != 2 )
            throw std::runtime_error( "could not parse size2: " + input );
        T w, h;
        try {
            w = boost::lexical_cast<T>( values[0] );
            h = boost::lexical_cast<T>( values[1] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "could not parse size2: " + input );
        }

        return SizeType( w, h );
    }

    SizeType operator+( const size2t<T, SizeType>& a ) { return SizeType( xsize + a.xsize, ysize + a.ysize ); }

    SizeType operator-( const size2t<T, SizeType>& a ) { return SizeType( xsize - a.xsize, ysize - a.ysize ); }

    SizeType operator*( const size2t<T, SizeType>& a ) { return SizeType( xsize * a.xsize, ysize * a.ysize ); }

    SizeType& operator+=( const size2t<T, SizeType>& a ) {
        xsize += a.xsize;
        ysize += a.ysize;
        return *static_cast<SizeType*>( this );
    }

    SizeType& operator-=( const size2t<T, SizeType>& a ) {
        xsize -= a.xsize;
        ysize -= a.ysize;
        return *static_cast<SizeType*>( this );
    }

    SizeType& operator*=( const size2t<T, SizeType>& a ) {
        xsize *= a.xsize;
        ysize *= a.ysize;
        return *static_cast<SizeType*>( this );
    }

    SizeType& operator+=( T a ) {
        xsize += a;
        ysize += a;
        return *static_cast<SizeType*>( this );
    }

    SizeType& operator-=( T a ) {
        xsize -= a;
        ysize -= a;
        return *static_cast<SizeType*>( this );
    }

    SizeType& operator*=( T a ) {
        xsize *= a;
        ysize *= a;
        return *static_cast<SizeType*>( this );
    }

    SizeType& operator/=( T a ) {
        xsize /= a;
        ysize /= a;
        return *static_cast<SizeType*>( this );
    }

    std::string str() const;
};

//////////////
// Operators
//////////////
template <typename T, typename SizeType>
inline SizeType operator+( T x, const size2t<T, SizeType>& a ) {
    return SizeType( x + a.xsize, x + a.ysize );
}

template <typename T, typename SizeType, typename VectorType>
inline VectorType operator+( const vector2t<T, VectorType>& a, const size2t<T, SizeType>& b ) {
    return VectorType( a.x + b.xsize, a.y + b.ysize );
}

template <typename T, typename SizeType, typename VectorType>
inline VectorType operator+( const size2t<T, SizeType>& a, const vector2t<T, VectorType>& b ) {
    return VectorType( a.xsize + b.x, a.ysize + b.y );
}

template <typename T, typename SizeType>
inline SizeType operator-( T x, const size2t<T, SizeType>& a ) {
    return SizeType( x - a.xsize, x - a.ysize );
}

template <typename T, typename SizeType, typename VectorType>
inline VectorType operator-( const vector2t<T, VectorType>& a, const size2t<T, SizeType>& b ) {
    return VectorType( a.x - b.xsize, a.y - b.ysize );
}

template <typename T, typename SizeType, typename VectorType>
inline VectorType operator-( const size2t<T, SizeType>& a, const vector2t<T, VectorType>& b ) {
    return VectorType( a.xsize - b.x, a.ysize - b.y );
}

template <typename T, typename SizeType>
inline SizeType operator*( T x, const size2t<T, SizeType>& a ) {
    return SizeType( x * a.xsize, x * a.ysize );
}

template <typename T, typename SizeType, typename VectorType>
inline VectorType operator*( const vector2t<T, VectorType>& a, const size2t<T, SizeType>& b ) {
    return VectorType( a.x * b.xsize, a.y * b.ysize );
}

template <typename T, typename SizeType, typename VectorType>
inline VectorType operator*( const size2t<T, SizeType>& a, const vector2t<T, VectorType>& b ) {
    return VectorType( a.xsize * b.x, a.ysize * b.y );
}

template <typename T, typename SizeType>
inline bool operator==( const size2t<T, SizeType>& a, const size2t<T, SizeType>& b ) {
    return a.ysize == b.ysize && a.xsize == b.xsize;
}

template <typename T, typename SizeType>
inline bool operator!=( const size2t<T, SizeType>& a, const size2t<T, SizeType>& b ) {
    return a.ysize != b.ysize || a.xsize != b.xsize;
}

template <typename T, typename SizeType, typename CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const size2t<T, SizeType>& v ) {
    out << "(size2t " << v.xsize << ", " << v.ysize << " )";
    return out;
}

template <typename T, typename SizeType>
inline std::string size2t<T, SizeType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

} // namespace graphics2d
} // namespace frantic
