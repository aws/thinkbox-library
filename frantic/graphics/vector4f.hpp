// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/transform4f.hpp>

#ifdef MAX_VERSION
#include <Quat.h>
#endif

namespace frantic {
namespace graphics {

template <typename FloatType>
class vector4t {
  public:
    FloatType x, y, z, w;
    typedef FloatType float_type;

    vector4t()
        : x( 0 )
        , y( 0 )
        , z( 0 )
        , w( 0 ) {}

    vector4t( float_type _x, float_type _y, float_type _z, float_type _w )
        : x( _x )
        , y( _y )
        , z( _z )
        , w( _w ) {}

    explicit vector4t( const float_type* ptr )
        : x( ptr[0] )
        , y( ptr[1] )
        , z( ptr[2] )
        , w( ptr[3] ) {}

#ifdef MAX_VERSION
    vector4t( const Quat& p ) {
        x = p.x;
        y = p.y;
        z = p.z;
        w = p.w;
    }

    vector4t& operator=( const Quat& p ) {
        x = p.x;
        y = p.y;
        z = p.z;
        w = p.w;
        return *this;
    }

    vector4t( const AngAxis& a ) {
        x = a.axis.x;
        y = a.axis.y;
        z = a.axis.z;
        w = a.angle;
    }

    vector4t& operator=( const AngAxis& a ) {
        x = a.axis.x;
        y = a.axis.y;
        z = a.axis.z;
        w = a.angle;
        return *this;
    }
#endif

    static vector4t from_random() {
        float_type m = 1.f / RAND_MAX;
        return vector4t( rand() * m - 0.5f, rand() * m - 0.5f, rand() * m - 0.5f, rand() * m - 0.5f );
    }

    float_type& operator[]( int i ) { return ( &x )[i]; }
    const float_type& operator[]( int i ) const { return ( &x )[i]; }

    float_type magnitude_squared() const { return x * x + y * y + z * z + w * w; }
    float_type magnitude() const {
        using namespace std;
        return sqrt( x * x + y * y + z * z + w * w );
    }

    bool is_zero() const { return x == 0 && y == 0 && z == 0 && w == 0; }

    vector4t operator-() const { return vector4t( -x, -y, -z, -w ); }

    vector4t& operator+=( const vector4t& v ) {
        x += v.x;
        y += v.y;
        z += v.z;
        w += v.w;
        return *this;
    }
    vector4t& operator-=( const vector4t& v ) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        w -= v.w;
        return *this;
    }
    vector4t& operator*=( float_type f ) {
        x *= f;
        y *= f;
        z *= f;
        w *= f;
        return *this;
    }
    vector4t& operator/=( float_type f ) {
        x /= f;
        y /= f;
        z /= f;
        w /= f;
        return *this;
    }
    bool operator==( const vector4t& v ) const { return x == v.x && y == v.y && z == v.z && w == v.w; }
    bool operator!=( const vector4t& v ) const { return x != v.x || y != v.y || z != v.z || w != v.w; }

    static vector4t cross( const vector4t& v1, const vector4t& v2, const vector4t& c );
    static float_type dot( const vector4t& v1, const vector4t& v2 );
    static vector4t component_mul( const vector4t& lhs, const vector4t& rhs );

    std::string str() const;

    // TODO: Convert anything using this to use quat4f/quat4fd, and remove from here!
    // Quaternion specific functions
    vector3t<float_type> quaternion_basis_vector( int i ) const;
    transform4t<float_type> quaternion_to_matrix() const;
    static vector4t<float_type> quaternion_from_rotation( const transform4t<float_type>& R );
    static vector4t<float_type> slerp( const vector4t<float_type>& qa, const vector4t<float_type>& qb, float_type t );
};

template <typename FloatType>
inline vector4t<FloatType> operator+( const vector4t<FloatType>& v1, const vector4t<FloatType>& v2 ) {
    return vector4t<FloatType>( v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w );
}
template <typename FloatType>
inline vector4t<FloatType> operator-( const vector4t<FloatType>& v1, const vector4t<FloatType>& v2 ) {
    return vector4t<FloatType>( v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w );
}
template <typename FloatType>
inline vector4t<FloatType> operator*( const vector4t<FloatType>& v, FloatType f ) {
    return vector4t<FloatType>( v ) *= f;
}
template <typename FloatType>
inline vector4t<FloatType> operator*( FloatType f, const vector4t<FloatType>& v ) {
    return v * f;
}
template <typename FloatType>
inline vector4t<FloatType> operator/( const vector4t<FloatType>& v, FloatType f ) {
    return vector4t<FloatType>( v ) /= f;
}

/**
 * Returns a column from the matrix equivalent to this quaternion.
 */
template <typename FloatType>
inline vector3t<FloatType> vector4t<FloatType>::quaternion_basis_vector( int i ) const {
    FloatType m = magnitude_squared();
    if( m <= 1e-6f ) // TODO: What should it actually do?
        return vector3t<FloatType>( 0 );

    switch( i ) {
    case 0:
        return vector3t<FloatType>( 1 - 2 * ( y * y + z * z ) / m, 2 * ( x * y + w * z ) / m,
                                    2 * ( x * z - w * y ) / m );
    case 1:
        return vector3t<FloatType>( 2 * ( x * y - w * z ) / m, 1 - 2 * ( x * x + z * z ) / m,
                                    2 * ( y * z + w * x ) / m );
    case 2:
        return vector3t<FloatType>( 2 * ( x * z + w * y ) / m, 2 * ( y * z - w * x ) / m,
                                    1 - 2 * ( x * x + y * y ) / m );
    default:
        throw std::range_error( "vector4f.quaternion_column: Index " + boost::lexical_cast<std::string>( i ) +
                                " out of range" );
    }
}

template <typename FloatType>
inline vector4t<FloatType> vector4t<FloatType>::cross( const vector4t<FloatType>& a, const vector4t<FloatType>& b,
                                                       const vector4t<FloatType>& c ) {
    // 4D Cross Product, yields a vector orthogonal to the other three vectors.
    //| i  j  k  h|
    //|x1 y1 z1 w1|				For example, the w component is:	|x1 y1 y1|
    //|x2 y2 z2 w2|													|x2
    //y2 z2|
    //|x3 y3 z3 w3|													|x3
    //y3 z3|
    FloatType temp1 = ( b.z * c.w - b.w * c.z );
    FloatType temp2 = ( a.z * c.w - a.w * c.z );
    FloatType temp3 = ( a.z * b.w - a.w * b.z );
    FloatType x = a.y * temp1 - b.y * temp2 + c.y * temp3;

    FloatType y = a.x * temp1 - b.x * temp2 + c.x * temp3;

    temp1 = ( b.x * c.y - b.y * c.x );
    temp2 = ( a.x * c.y - a.y * c.x );
    temp3 = ( a.x * b.y - a.y * b.x );
    FloatType z = a.w * temp1 - b.w * temp2 + c.w * temp3;

    FloatType w = a.z * temp1 - b.z * temp2 + c.z * temp3;

    return vector4t<FloatType>( -x, -y, -z, -w );
}

template <typename FloatType>
inline vector4t<FloatType> vector4t<FloatType>::component_mul( const vector4t<FloatType>& lhs,
                                                               const vector4t<FloatType>& rhs ) {
    return vector4t<FloatType>( lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z, lhs.w * rhs.w );
}

// See "Visualizing Quaternions", pg. 148 by Andrew J. Hanson for details.
template <typename FloatType>
inline vector4t<FloatType> vector4t<FloatType>::quaternion_from_rotation( const transform4t<FloatType>& R ) {
    using namespace std;

    FloatType X, Y, Z, W;
    FloatType trace = R[0] + R[5] + R[10];
    if( trace > 0 ) {
        W = 0.5f * sqrt( 1 + trace );

        X = 0.25f * ( R[6] - R[9] ) / W;
        Y = 0.25f * ( R[8] - R[2] ) / W;
        Z = 0.25f * ( R[1] - R[4] ) / W;
    } else {
        if( R[0] > R[5] && R[0] > R[10] ) {
            X = 0.5f * sqrt( 1 + R[0] - R[5] - R[10] );

            W = 0.25f * ( R[6] - R[9] ) / X;
            Y = 0.25f * ( R[1] + R[4] ) / X;
            Z = 0.25f * ( R[2] + R[8] ) / X;
        } else if( R[5] > R[10] ) {
            Y = 0.5f * sqrt( 1 + R[5] - R[0] - R[10] );

            W = 0.25f * ( R[8] - R[2] ) / Y;
            X = 0.25f * ( R[1] + R[4] ) / Y;
            Z = 0.25f * ( R[6] + R[9] ) / Y;
        } else {
            Z = 0.5f * sqrt( 1 + R[10] - R[0] - R[5] );

            W = 0.25f * ( R[1] - R[4] ) / Z;
            X = 0.25f * ( R[2] + R[8] ) / Z;
            Y = 0.25f * ( R[6] + R[9] ) / Z;
        }
    }

    return vector4t( X, Y, Z, W );
}

template <typename FloatType>
inline transform4t<FloatType> vector4t<FloatType>::quaternion_to_matrix() const {
    FloatType m = magnitude_squared();
    if( m <= 1e-6f ) // TODO: What should it actually do?
        return transform4t<FloatType>::zero();

    FloatType xx = x * x, yy = y * y, zz = z * z;
    FloatType xy = x * y, xz = x * z, yz = y * z;
    FloatType wx = w * x, wy = w * y, wz = w * z;

    return transform4t<FloatType>( 1 - 2 * ( yy + zz ) / m, // Column 1
                                   2 * ( xy + wz ) / m, 2 * ( xz - wy ) / m, 0,

                                   2 * ( xy - wz ) / m, // Column 2
                                   1 - 2 * ( xx + zz ) / m, 2 * ( yz + wx ) / m, 0,

                                   2 * ( xz + wy ) / m, // Colum 3
                                   2 * ( yz - wx ) / m, 1 - 2 * ( xx + yy ) / m, 0,

                                   0, 0, 0, 1 );
}

template <typename FloatType>
inline vector4t<FloatType> vector4t<FloatType>::slerp( const vector4t<FloatType>& qa, const vector4t<FloatType>& qb,
                                                       FloatType t ) {
    using namespace std;
    // quaternion to return
    vector4t<FloatType> qm;
    // Calculate angle between them.
    FloatType cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;
    // if qa=qb or qa=-qb then theta = 0 and we can return qa
    if( abs( cosHalfTheta ) >= 1 ) {
        qm.w = qa.w;
        qm.x = qa.x;
        qm.y = qa.y;
        qm.z = qa.z;
        return qm;
    }
    // Calculate temporary values.
    FloatType halfTheta = acos( cosHalfTheta );
    FloatType sinHalfTheta = sqrt( 1 - cosHalfTheta * cosHalfTheta );
    // if theta = 180 degrees then result is not fully defined
    // we could rotate around any axis normal to qa or qb
    if( abs( sinHalfTheta ) < 1e-5f ) { // fabs is floating point absolute
        qm.w = ( qa.w * 0.5f + qb.w * 0.5f );
        qm.x = ( qa.x * 0.5f + qb.x * 0.5f );
        qm.y = ( qa.y * 0.5f + qb.y * 0.5f );
        qm.z = ( qa.z * 0.5f + qb.z * 0.5f );
        return qm;
    }
    FloatType ratioA = sin( ( 1 - t ) * halfTheta ) / sinHalfTheta;
    FloatType ratioB = sin( t * halfTheta ) / sinHalfTheta;
    // calculate Quaternion.
    qm.w = ( qa.w * ratioA + qb.w * ratioB );
    qm.x = ( qa.x * ratioA + qb.x * ratioB );
    qm.y = ( qa.y * ratioA + qb.y * ratioB );
    qm.z = ( qa.z * ratioA + qb.z * ratioB );
    return qm;
}

template <typename FloatType>
inline static FloatType dot( const vector4t<FloatType>& v1, const vector4t<FloatType>& v2 ) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
}

template <class CharType, typename FloatType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& o, const vector4t<FloatType>& v ) {
    return o << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "]";
}

// Only matrix multiplication with a vector on the right is defined, because we are using
// the convention of column vectors.
template <typename FloatType>
inline vector4t<FloatType> operator*( const transform4t<FloatType>& m, const vector4t<FloatType>& v ) {
    return vector4t<FloatType>( ( v.x * m[0] + v.y * m[4] + v.z * m[8] + v.w * m[12] ),
                                ( v.x * m[1] + v.y * m[5] + v.z * m[9] + v.w * m[13] ),
                                ( v.x * m[2] + v.y * m[6] + v.z * m[10] + v.w * m[14] ),
                                ( v.x * m[3] + v.y * m[7] + v.z * m[11] + v.w * m[15] ) );
}

template <typename FloatType>
std::string vector4t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

typedef vector4t<float> vector4f;
typedef vector4t<double> vector4fd;
} // namespace graphics
} // namespace frantic
