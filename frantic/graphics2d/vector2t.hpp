// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <limits>
#include <sstream>
#include <string>

namespace frantic {
namespace graphics2d {

template <typename T, typename VectorType>
class vector2t {
  public:
    T x, y;

    typedef T value_type;

    /**************
     * Constructors
     **************/
    vector2t()
        : x( 0 )
        , y( 0 ) {}
    vector2t( T X )
        : x( X )
        , y( X ) {}
    vector2t( T X, T Y )
        : x( X )
        , y( Y ) {}
    explicit vector2t( const std::pair<T, T>& p )
        : x( p.first )
        , y( p.second ) {}
    //		explicit template<typename U>
    //		vector2t( const vector2t<U> & vec ) : x( (T)vec.x ), y( (T)vec.y )	{}

    /**************
     * Mutators
     **************/
    void set( T X, T Y ) {
        x = X;
        y = Y;
    }
    void set_x( T X ) { x = X; }
    void set_y( T Y ) { y = Y; }

    /**************
     * Accessors
     **************/
    const T& get_x() const { return x; }
    T& get_x() { return x; }
    const T& get_y() const { return y; }
    T& get_y() { return y; }

    static VectorType minvalue() { return VectorType( ( std::numeric_limits<T>::min )() ); }
    static VectorType maxvalue() { return VectorType( ( std::numeric_limits<T>::max )() ); }

    static T cross( const VectorType& a, const VectorType& b ) { return a.x * b.y - b.x * a.y; }
    static T dot( const VectorType& a, const VectorType& b ) { return a.x * b.x + a.y * b.y; }
    static double dot_double( const VectorType& a, const VectorType& b ) {
        return double( a.x ) * double( b.x ) + double( a.y ) * double( b.y );
    }

    T get_magnitude_squared() const { return x * x + y * y; }
    float get_magnitude() const { return sqrtf( (float)get_magnitude_squared() ); }

    static T manhattan_distance( const VectorType& a, const VectorType& b ) {
        return (T)fabs( (float)( a.x - b.x ) ) + (T)fabs( (float)( a.y - b.y ) );
    }

    static T distance_squared( const VectorType& a, const VectorType& b ) {
        return ( a.x - b.x ) * ( a.x - b.x ) + ( a.y - b.y ) * ( a.y - b.y );
    }

    static float distance( const VectorType& a, const VectorType& b ) {
        return sqrtf( (float)vector2t::distance_squared( a, b ) );
    }

    static T signed_area( const VectorType& p1, const VectorType& p2, const VectorType& p3 );

    std::pair<T, T> to_pair() const { return std::make_pair( x, y ); }

    template <class VectorTypeU>
    VectorTypeU convert() {
        return VectorTypeU( (typename VectorTypeU::value_type)get_x(), (typename VectorTypeU::value_type)get_y() );
    }

    static inline bool find_intersection( const VectorType& p1, const VectorType& p2, const VectorType& p3,
                                          const VectorType& p4, VectorType* pi );

    VectorType left_normal() const { return VectorType( -y, x ); }

    VectorType right_normal() const { return VectorType( y, -x ); }

    /////////////////////
    // Operators
    /////////////////////

    T& operator[]( int i ) { return ( &x )[i]; }
    T& operator[]( size_t i ) { return ( &x )[i]; }
    const T& operator[]( int i ) const { return ( &x )[i]; }
    const T& operator[]( size_t i ) const { return ( &x )[i]; }

    VectorType& operator+=( T a ) {
        x += a;
        y += a;
        return *static_cast<VectorType*>( this );
    }

    VectorType& operator-=( T a ) {
        x -= a;
        y -= a;
        return *static_cast<VectorType*>( this );
    }

    VectorType& operator*=( T a ) {
        x *= a;
        y *= a;
        return *static_cast<VectorType*>( this );
    }

    VectorType& operator/=( T a ) {
        x /= a;
        y /= a;
        return *static_cast<VectorType*>( this );
    }

    VectorType& operator+=( const VectorType& a ) {
        x += a.x;
        y += a.y;
        return *static_cast<VectorType*>( this );
    }

    VectorType& operator-=( const VectorType& a ) {
        x -= a.x;
        y -= a.y;
        return *static_cast<VectorType*>( this );
    }

    VectorType& operator*=( const VectorType& a ) {
        x *= a.x;
        y *= a.y;
        return *static_cast<VectorType*>( this );
    }

    bool operator==( const VectorType& a ) const { return ( x == a.x ) && ( y == a.y ); }

    bool operator!=( const VectorType& a ) const { return ( x != a.x ) || ( y != a.y ); }

    std::string str() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }
};

template <typename T, typename VectorType>
inline std::ostream& operator<<( std::ostream& out, const vector2t<T, VectorType>& v ) {
    out << "(vec2t " << v.x << ", " << v.y << " )";
    return out;
}

template <typename T, typename VectorType>
inline VectorType operator+( const vector2t<T, VectorType>& a, const vector2t<T, VectorType>& b ) {
    return VectorType( a.x + b.x, a.y + b.y );
}

template <typename T, typename VectorType>
inline VectorType operator-( const vector2t<T, VectorType>& a, const vector2t<T, VectorType>& b ) {
    return VectorType( a.x - b.x, a.y - b.y );
}

template <typename T, typename VectorType>
inline VectorType operator*( const vector2t<T, VectorType>& a, const vector2t<T, VectorType>& b ) {
    return VectorType( a.x * b.x, a.y * b.y );
}

template <typename T, typename VectorType>
inline VectorType operator+( const vector2t<T, VectorType>& a, const std::pair<T, T>& p ) {
    return VectorType( a.x + p.first, a.y + p.second );
}

template <typename T, typename VectorType>
inline VectorType operator-( const vector2t<T, VectorType>& a, const std::pair<T, T>& p ) {
    return VectorType( a.x - p.first, a.y - p.second );
}

template <typename T, typename VectorType>
inline VectorType operator*( const vector2t<T, VectorType>& a, const std::pair<T, T>& p ) {
    return VectorType( a.x * p.first, a.y * p.second );
}

template <typename T, typename VectorType>
inline VectorType operator+( T val, const vector2t<T, VectorType>& b ) {
    return VectorType( val + b.x, val + b.y );
}

template <typename T, typename VectorType>
inline VectorType operator-( T val, const vector2t<T, VectorType>& b ) {
    return VectorType( val - b.x, val - b.y );
}

template <typename T, typename VectorType>
inline VectorType operator*( T val, const vector2t<T, VectorType>& b ) {
    return VectorType( val * b.x, val * b.y );
}

template <typename T, typename VectorType>
inline VectorType operator-( const vector2t<T, VectorType>& a ) {
    return VectorType( -( a.x ), -( a.y ) );
}

template <typename T, typename VectorType>
inline T distance_squared( const vector2t<T, VectorType>& a, const vector2t<T, VectorType>& b ) {
    return ( a - b ).get_magnitude_squared();
}

template <typename T, typename VectorType>
inline float distance( const vector2t<T, VectorType>& a, const vector2t<T, VectorType>& b ) {
    return ( a - b ).get_magnitude();
}

/*! Find the intersection point of the line segments defined by end point p1<->p2 and p3<->p4.
 *  Returns false if the line segments do not intersect.
 */
template <typename T, typename VectorType>
bool vector2t<T, VectorType>::find_intersection( const VectorType& p1, const VectorType& p2, const VectorType& p3,
                                                 const VectorType& p4, VectorType* pi ) {
    double denom;
    double ua;

    denom = ( ( p4.get_y() - p3.get_y() ) * ( p2.get_x() - p1.get_x() ) ) -
            ( ( p4.get_x() - p3.get_x() ) * ( p2.get_y() - p1.get_y() ) );

    if( denom == 0.0 )
        return false;

    ua = ( ( p4.get_x() - p3.get_x() ) * ( p1.get_y() - p3.get_y() ) ) -
         ( ( p4.get_y() - p3.get_y() ) * ( p1.get_x() - p3.get_x() ) );
    ua /= denom;

    pi->set_x( (T)( p1.get_x() + ( p2.get_x() - p1.get_x() ) * ua ) );
    pi->set_y( (T)( p1.get_y() + ( p2.get_y() - p1.get_y() ) * ua ) );

    return true;
}

/*! Compute the signed for the vector v1,v2,v3.  If v3 is left of the line passing through v1,v2 then the
 *	area will be positive
 */
template <typename T, typename VectorType>
T vector2t<T, VectorType>::signed_area( const VectorType& p1, const VectorType& p2, const VectorType& p3 ) {
    return ( p2.get_x() - p1.get_x() ) * ( p3.get_y() - p1.get_y() ) -
           ( p3.get_x() - p1.get_x() ) * ( p2.get_y() - p1.get_y() );
}

} // namespace graphics2d
} // namespace frantic
