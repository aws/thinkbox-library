// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>

namespace frantic {
namespace graphics {

namespace detail {
// Metafunction which returns the first or second, depending on Which
template <class S, class T, bool Which>
struct choose_type;
template <class S, class T>
struct choose_type<S, T, false> {
    typedef S type;
};
template <class S, class T>
struct choose_type<S, T, true> {
    typedef T type;
};

// Metafunction which returns the bigger of the two types
template <class S, class T>
struct biggest_type : public choose_type<S, T, sizeof( S ) < sizeof( T )> {};
} // namespace detail

template <class FloatType>
class quat4t {
  public:
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;
    typedef size3t<FloatType> size3f_type;

    // (w) is the real part, m_data[0-2] (x,y,z) is the imaginary part. W is last to match the Max layout
    float_type x, y, z, w;

  public:
    quat4t()
        : x( 0 )
        , y( 0 )
        , z( 0 )
        , w( 1 ) {}

    quat4t( float_type realPart, float_type i, float_type j, float_type k )
        : x( i )
        , y( j )
        , z( k )
        , w( realPart ) {}

    quat4t( float_type realPart, const vector3f_type& imaginaryPart )
        : x( imaginaryPart.x )
        , y( imaginaryPart.y )
        , z( imaginaryPart.z )
        , w( realPart ) {}

    /** Implicit conversion from smaller float type is allowed. */
    template <class OtherFloatType>
    quat4t( const quat4t<OtherFloatType>& t,
            typename boost::enable_if_c<( sizeof( OtherFloatType ) < sizeof( FloatType ) ), int*>::type = 0 )
        : w( t.w )
        , x( t.x )
        , y( t.y )
        , z( t.z ) {}

    /** Explicit cast required for conversion from larger float type. */
    template <class OtherFloatType>
    explicit quat4t( const quat4t<OtherFloatType>& t,
                     typename boost::enable_if_c<( sizeof( OtherFloatType ) > sizeof( FloatType ) ), int*>::type = 0 )
        : w( float_type( t.w ) )
        , x( float_type( t.x ) )
        , y( float_type( t.y ) )
        , z( float_type( t.z ) ) {}

    float_type real_part() const { return w; }

    vector3f_type vector_part() const { return vector3f_type( x, y, z ); }

    float_type magnitude_squared() const { return x * x + y * y + z * z + w * w; }

    float_type magnitude() const { return sqrt( magnitude_squared() ); }

    float_type max_abs_component() const { return ( std::max )( x, ( std::max )( y, ( std::max )( z, w ) ) ); }

    void normalize() {
        float_type mag = magnitude();
        if( mag != 0 ) {
            x /= mag;
            y /= mag;
            z /= mag;
            w /= mag;
        } else {
            set_to_identity();
        }
    }

    quat4t to_normalized() const {
        quat4t result = *this;
        result.normalize();
        return result;
    }

    quat4t get_inverse() const { return quat4t( *this ).invert(); }

    quat4t& invert() {
        const float_type m = magnitude_squared();
        x = -x / m;
        y = -y / m;
        z = -z / m;
        w /= m;
        return *this;
    }

    static quat4t conjugate( const quat4t& val ) { return quat4t( val.w, -val.x, -val.y, -val.z ); }

    static quat4t from_identity() { return quat4t().set_to_identity(); }

    quat4t& set_to_identity() {
        x = 0;
        y = 0;
        z = 0;
        w = 1.f;
        return *this;
    }

    static float_type dot( const quat4t& lhs, const quat4t& rhs ) {
        return lhs.w * rhs.w + lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    static double dot_double( const quat4t& lhs, const quat4t& rhs ) {
        return static_cast<double>( lhs.w ) * static_cast<double>( rhs.w ) +
               static_cast<double>( lhs.x ) * static_cast<double>( rhs.x ) +
               static_cast<double>( lhs.y ) * static_cast<double>( rhs.y ) +
               static_cast<double>( lhs.z ) * static_cast<double>( rhs.z );
    }

    /**
     * Interpolates between endpoint quat4ts using spherical linear interpolation. This will interpolate with constant
     * angular velocity, but is not commutative.
     *
     * @note lhs and rhs should be unit length or else something broken will happen.
     *
     * @param lhs The unit quat4t at the 0 side of the interpolation
     * @param rhs The unit quat4t at the 1 side of the interpolation
     * @param t The interpolation parameter [0,1]
     * @return A unit quat4t the is interpolated between lhs & rhs
     */
    static quat4t slerp( const quat4t& lhs, const quat4t& rhs, float_type t );

    /**
     * A normalized lerp of the quat4t. Corresponds to a valid interpolation between endpoints, but with non-constant
     * angular velocity. What this version does give you is A) Performance and B) Commutivity compared to slerp. Because
     * its commutative, we can use it for trilinear interpolation and get a consistent result.
     *
     * nlerp(a,b,t) = lerp(a,b,t) / ||lerp(a,b,t)||
     *
     * @note lhs and rhs should be unit length or else something broken will happen.
     *
     * @param lhs unit The quat4t at the 0 side of the interpolation
     * @param rhs unit The quat4t at the 1 side of the interpolation
     * @param t The interpolation parameter [0,1]
     * @return A unit quat4t the is interpolated between lhs & rhs
     */
    static quat4t nlerp( const quat4t& lhs, const quat4t& rhs, float_type t );

    /**
     * This version of nlerp does not check whether t is in [0,1]. Otherwise it is exactly like nlerp()
     */
    static quat4t nlerp_nocheck( const quat4t& lhs, const quat4t& rhs, float_type t );

    /**
     * Creates a quat4t by composing three rotations around the x-axis, then y-axis, then z-axis.
     */
    static quat4t from_euler_angles( float_type xAngleRads, float_type yAngleRads, float_type zAngleRads ) {
        // We compose a rotation around X, then Y, then Z. From the from_angle_axis( a, [nx,ny,nz] ) equation we see
        // this is straightforward for a single axis. quat( cos(a/2), nx*sin(a/2), ny*sin(a/2), nz*sin(a/2) )
        //
        // Multiplying two quat4ts will compose them:
        // result.w = lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z;
        // result.x = lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y;
        // result.y = lhs.w * rhs.y - lhs.x * rhs.z + lhs.y * rhs.w + lhs.z * rhs.x;
        // result.z = lhs.w * rhs.z + lhs.x * rhs.y - lhs.y * rhs.x + lhs.z * rhs.w;

        // Therefore:
        // float hX = 0.5f * x, hY = 0.5f * y, hZ = 0.5f * z;
        // float w = cos(hZ)*cos(hY);
        // float x = -sin(hZ)*sin(hY);
        // float y = cos(hZ)*sin(hY);
        // float z = sin(hZ)*cos(hY);
        // w = w*cos(hX) - x*sin(hX);
        // x = x*cos(hX) + w*sin(hX);
        // y = y*cos(hX) + z*sin(hX);
        // z = z*cos(hX) - y*sin(hX);

        float_type hX = 0.5f * xAngleRads, cX = cos( hX ), sX = sin( hX );
        float_type hY = 0.5f * yAngleRads, cY = cos( hY ), sY = sin( hY );
        float_type hZ = 0.5f * zAngleRads, cZ = cos( hZ ), sZ = sin( hZ );

        float_type qw = cZ * cY;
        float_type qx = sZ * sY; // I omitted the sign here and moved it into the final calculation.
        float_type qy = cZ * sY;
        float_type qz = sZ * cY;

        return quat4t( qw * cX + qx * sX, qw * sX - qx * cX, qy * cX + qz * sX, qz * cX - qy * sX );
    }

    /**
     * Creates a quat4t from an angle and a normalized direction vector.
     */
    static quat4t from_angle_axis( float_type angleRadians, const vector3f_type& normalizedAxis ) {
        float_type halfAngle = 0.5f * angleRadians;
        float_type cosHalfAngle = cos( halfAngle );
        float_type sinHalfAngle = sin( halfAngle );

        return quat4t( cosHalfAngle, sinHalfAngle * normalizedAxis.x, sinHalfAngle * normalizedAxis.y,
                       sinHalfAngle * normalizedAxis.z );
        // This quat is normalized (therefore a valid orientation) because:
        // magnitude = sqrt( cos(a)^2 + sin(a)^2*x*x + sin(a)^2*y*y + sin(a)^2*z*z )
        //           = sqrt( cos(a)^2 + sin(a)^2 * (x*x + y*y + z*z) )
        //           = sqrt( cos(a)^2 + sin(a)^2 * 1.0 )
        //           = sqrt( cos(a)^2 + sin(a)^2 ) = sqrt( 1.0 ) = 1.0
    }

    /**
     * Creates a quat4t from a rotation matrix.  Does not check if rotation matrix is valid
     */
    static quat4t from_transform4f( const transform4t<float_type>& transform ) {
        float_type m00 = transform[0];
        float_type m10 = transform[1];
        float_type m20 = transform[2];
        float_type m01 = transform[4];
        float_type m11 = transform[5];
        float_type m21 = transform[6];
        float_type m02 = transform[8];
        float_type m12 = transform[9];
        float_type m22 = transform[10];

        quat4t result;

        const float_type trace = m00 + m11 + m22;
        if( trace > 0 ) {
            const float_type s = sqrt( trace + 1 ) * 2;
            result = quat4t( s * 0.25f, ( m21 - m12 ) / s, ( m02 - m20 ) / s, ( m10 - m01 ) / s );
        } else if( m00 > m11 && m00 > m22 ) {
            const float_type s = sqrt( m00 - m11 - m22 + 1 ) * 2;
            result = quat4t( ( m21 - m12 ) / s, s * 0.25f, ( m01 + m10 ) / s, ( m02 + m20 ) / s );
        } else if( m11 > m22 ) {
            const float_type s = sqrt( m11 - m00 - m22 + 1 ) * 2;
            result = quat4t( ( m02 - m20 ) / s, ( m01 + m10 ) / s, s * 0.25f, ( m12 + m21 ) / s );
        } else {
            const float_type s = sqrt( m22 - m00 - m11 + 1 ) * 2;
            result = quat4t( ( m10 - m01 ) / s, ( m02 + m20 ) / s, ( m12 + m21 ) / s, s * 0.25f );
        }

        // TODO: Found this normalize was necessary, maybe a bug in the above?
        result.normalize();
        return result;
    }

    /**
     * Creates the quat4t representation of a 3D coordinate system. Assumes the input vectors are orthonormal.
     */
    static quat4t from_coord_sys( const vector3f_type& xAxis, const vector3f_type& yAxis, const vector3f_type& zAxis ) {
        float_type q0, q1, q2, q3;
        float_type denom;

        float_type trace = xAxis.x + yAxis.y + zAxis.z;
        if( trace > 0 ) {
            q0 = sqrt( 1 + trace );
            denom = q0 + q0;

            q0 = 0.5f * q0;
            q1 = ( yAxis.z - zAxis.y ) / denom;
            q2 = ( zAxis.x - xAxis.z ) / denom;
            q3 = ( xAxis.y - yAxis.x ) / denom;
        } else if( xAxis.x >= yAxis.y && xAxis.x >= zAxis.z ) {
            q1 = sqrt( 1 + xAxis.x - yAxis.y - zAxis.z );
            denom = q1 + q1;

            q0 = ( yAxis.z - zAxis.y ) / denom;
            q1 = 0.5f * q1;
            q2 = ( xAxis.y + yAxis.x ) / denom;
            q3 = ( zAxis.x + xAxis.z ) / denom;
        } else if( yAxis.y >= zAxis.z ) {
            q2 = sqrt( 1 - xAxis.x + yAxis.y - zAxis.z );
            denom = q2 + q2;

            q0 = ( zAxis.x - xAxis.z ) / denom;
            q1 = ( yAxis.x + xAxis.y ) / denom;
            q2 = 0.5f * q2;
            q3 = ( zAxis.y + yAxis.z ) / denom;
        } else {
            q3 = sqrt( 1 - xAxis.x - yAxis.y + zAxis.z );
            denom = q3 + q3;

            q0 = ( xAxis.y - yAxis.x ) / denom;
            q1 = ( xAxis.z + zAxis.x ) / denom;
            q2 = ( zAxis.y + yAxis.z ) / denom;
            q3 = 0.5f * q3;
        }

        return quat4t( q0, q1, q2, q3 );
    }

    template <typename OtherFloatType>
    void as_transform4f( transform4t<OtherFloatType>& outTransform ) const {
        // Use the bigger float type of the quat4f/transform4f as the precision to do the computation
        typedef typename detail::biggest_type<FloatType, OtherFloatType>::type BigFloatType;
        BigFloatType twoX = x + x;
        BigFloatType twoY = y + y;
        BigFloatType twoZ = z + z;

        BigFloatType x_twoX = x * twoX;
        BigFloatType x_twoY = x * twoY;
        BigFloatType x_twoZ = x * twoZ;
        BigFloatType y_twoY = y * twoY;
        BigFloatType y_twoZ = y * twoZ;
        BigFloatType z_twoZ = z * twoZ;

        BigFloatType w_twoX = w * twoX;
        BigFloatType w_twoY = w * twoY;
        BigFloatType w_twoZ = w * twoZ;

        outTransform[0] = 1.f - ( y_twoY + z_twoZ );
        outTransform[1] = x_twoY + w_twoZ;
        outTransform[2] = x_twoZ - w_twoY;
        outTransform[3] = 0.f;

        outTransform[4] = x_twoY - w_twoZ;
        outTransform[5] = 1.f - ( x_twoX + z_twoZ );
        outTransform[6] = y_twoZ + w_twoX;
        outTransform[7] = 0.f;

        outTransform[8] = x_twoZ + w_twoY;
        outTransform[9] = y_twoZ - w_twoX;
        outTransform[10] = 1.f - ( x_twoX + y_twoY );
        outTransform[11] = 0.f;

        outTransform[12] = 0.f;
        outTransform[13] = 0.f;
        outTransform[14] = 0.f;
        outTransform[15] = 1.f;

        // TODO: something seems wonky about either this method or from_angle_axis.
        //  Calculating transform from angle-axis to quaternion seems to require outTransform
        //  to be transposed.  Need to determine what's going on.
    }

    transform4t<FloatType> to_transform4f() const {
        transform4t<FloatType> result;
        as_transform4f( result );
        return result;
    }

    /**
     * Converts the quat4t into an angle and normalized axis representation.
     * Notes:
     *   * Does NOT check if the quaternion itself is normalized.
     *   * If the quaternion is an identity (abs of real component is 1), the axis is assumed to be in the z-axis with
     * the sign equal to the real component (In theory, any axis is valid for this case since there's no rotation
     * anyway)
     */
    template <typename OtherFloatType>
    void as_angle_axis( OtherFloatType& outAngle, vector3t<OtherFloatType>& outAxis,
                        const float_type& epsilon = 0 ) const {
        // Use the bigger float type of the quat4f/transform4f as the precision to do the computation
        typedef typename detail::biggest_type<FloatType, OtherFloatType>::type BigFloatType;

        using std::abs;
        using std::acos;
        using std::sqrt;

        if( abs( w ) <= epsilon ) {
            // angle at 180 degrees
            outAngle = static_cast<OtherFloatType>( M_PI );
            outAxis = vector3t<OtherFloatType>( x, y, z );

        } else if( abs( w - 1 ) <= epsilon ) {
            // angle at 0 degrees: any axis is fine
            outAngle = 0;
            outAxis = vector3t<OtherFloatType>( 0, 0, 1 );

        } else if( abs( w + 1 ) <= epsilon ) {
            // angle at 0 degrees: any axis is fine
            outAngle = 0;
            outAxis = vector3t<OtherFloatType>( 0, 0, -1 );

        } else {
            // All other cases
            const BigFloatType scale = sqrt( 1.0f - w * w );

            outAngle = acos( w ) * 2;
            outAxis = vector3t<OtherFloatType>( x / scale, y / scale, z / scale );
        }
    }

    quat4t& operator*=( float_type f ) {
        w *= f;
        x *= f;
        y *= f;
        z *= f;
        return ( *this );
    }

    quat4t& operator+=( const quat4t& rhs ) {
        w += rhs.w;
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return ( *this );
    }

    quat4t& operator-=( const quat4t& rhs ) {
        w -= rhs.w;
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return ( *this );
    }

    std::string str() const;
};

template <class FloatType>
inline quat4t<FloatType> operator*( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs ) {
    // return quat4t(
    //	lhs.real_part()*rhs.real_part() - vector3f::dot(lhs.vector_part(),rhs.vector_part()),
    //	lhs.real_part()*rhs.vector_part() +
    //	rhs.real_part()*lhs.vector_part() +
    //	vector3f::cross(lhs.vector_part(), rhs.vector_part()) );

    quat4t<FloatType> result;
    result.w = lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z;
    result.x = lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y;
    result.y = lhs.w * rhs.y - lhs.x * rhs.z + lhs.y * rhs.w + lhs.z * rhs.x;
    result.z = lhs.w * rhs.z + lhs.x * rhs.y - lhs.y * rhs.x + lhs.z * rhs.w;
    return result;
}

template <class FloatType>
inline const quat4t<FloatType> operator+( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs ) {
    return quat4t<FloatType>( lhs ) += rhs;
}

template <class FloatType>
inline const quat4t<FloatType> operator-( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs ) {
    return quat4t<FloatType>( lhs ) -= rhs;
}

template <class FloatType>
inline const quat4t<FloatType> operator*( FloatType f, const quat4t<FloatType>& q ) {
    return quat4t<FloatType>( q ) *= f;
}

template <class FloatType>
inline const quat4t<FloatType> operator*( const quat4t<FloatType>& q, FloatType f ) {
    return f * q;
}

template <class FloatType>
inline bool operator==( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs ) {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
}

template <class FloatType>
inline bool operator!=( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs ) {
    return lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w;
}

template <class FloatType>
inline quat4t<FloatType> operator-( const quat4t<FloatType>& q ) {
    return quat4t<FloatType>( -q.w, -q.x, -q.y, -q.z );
}

/**
 * Rotates 'rhs' by the rotation/orientation described by the unit quaterion 'lhs'.
 * @note lhs must be a unit quat4t or badness happens
 * @param lhs The unit quat4t describing the rotation
 * @param rhs The point to rotate
 * @return The point 'rhs' rotated by 'lhs'
 */
template <class FloatType>
inline vector3t<FloatType> rotate_point( const quat4t<FloatType>& lhs, const vector3t<FloatType>& rhs ) {
    // TODO Test this!
    return ( ( lhs * quat4t<FloatType>( 0.f, rhs.x, rhs.y, rhs.z ) ) * quat4t<FloatType>::conjugate( lhs ) )
        .vector_part();

    // TODO Expand this and collapse shared terms.
}

template <class FloatType>
inline quat4t<FloatType> quat4t<FloatType>::slerp( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs,
                                                   float_type t ) {
    using std::acos;
    using std::cos;
    using std::sin;

    if( t <= 0.f )
        return lhs;
    if( t >= 1.f )
        return rhs;

    // Compute the cosine of the angle between the two vectors.
    float_type dotVal = quat4t::dot( lhs, rhs );

    // If the inputs are too close for comfort, switch to nlerp
    // TODO WHY?
    const float_type DOT_THRESHOLD = float_type( 0.9995 );
    if( dotVal > DOT_THRESHOLD )
        return nlerp_nocheck( lhs, rhs, t ); // Use nocheck version since we've checked upon entering slerp().

    dotVal = frantic::math::clamp( dotVal, float_type( -1.0 ),
                                   float_type( 1.0 ) ); // Robustness: Stay within domain of acos()
    float_type theta_0 = acos( dotVal );                // theta_0 = angle between input vectors
    float_type theta = theta_0 * t;                     // theta = angle between lhs and result

    quat4t<FloatType> newRhs = rhs - lhs * static_cast<float_type>( dotVal );
    newRhs.normalize(); // { lhs, newRhs } is now an orthonormal basis

    return lhs * static_cast<float_type>( cos( theta ) ) + newRhs * static_cast<float_type>( sin( theta ) );
}

// TODO Move this to cpp file
template <class FloatType>
inline quat4t<FloatType> quat4t<FloatType>::nlerp_nocheck( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs,
                                                           float_type t ) {
    quat4t<FloatType> result( lhs.w + t * ( rhs.w - lhs.w ), lhs.x + t * ( rhs.x - lhs.x ),
                              lhs.y + t * ( rhs.y - lhs.y ), lhs.z + t * ( rhs.z - lhs.z ) );

    float_type m = result.magnitude();

    // TODO Handle the case where m is ~ 0. Maybe swtich to slerp in that case?
    result.w /= m;
    result.x /= m;
    result.y /= m;
    result.z /= m;

    return result;
}

template <class FloatType>
inline quat4t<FloatType> quat4t<FloatType>::nlerp( const quat4t<FloatType>& lhs, const quat4t<FloatType>& rhs,
                                                   float_type t ) {
    if( t <= 0 )
        return lhs;
    if( t >= 1.f )
        return rhs;

    return nlerp_nocheck( lhs, rhs, t );
}

template <class CharType, class FloatType>
std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const quat4t<FloatType>& q ) {
    out << "[ (" << q.x << ", " << q.y << ", " << q.z << "), " << q.w << ']';
    return out;
}

template <class FloatType>
std::string quat4t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

/**
 * Calculates the logarithm of the quaternion
 */
template <class FloatType>
quat4t<FloatType> log( const quat4t<FloatType>& quat ) {
    using std::acos;
    using std::log;

    const vector3t<FloatType> vectorPart( quat.vector_part() );
    const FloatType realPart( quat.real_part() );
    const FloatType magnitudePart( quat.magnitude() );
    const FloatType vectorPartMagnitude( vectorPart.get_magnitude() );

    return quat4t<FloatType>( log( magnitudePart ),
                              vectorPartMagnitude > 0
                                  ? vectorPart * acos( realPart / magnitudePart ) / vectorPartMagnitude
                                  : vector3t<FloatType>( 0 ) );
}

/**
 * Calculates the exponential of the quaternion
 */
template <class FloatType>
quat4t<FloatType> exp( const quat4t<FloatType>& quat ) {
    using std::cos;
    using std::exp;
    using std::sin;

    const vector3t<FloatType> vectorPart( quat.vector_part() );
    const FloatType expReal( exp( quat.real_part() ) );
    const FloatType magnitudePart( quat.magnitude() );
    const FloatType vectorPartMagnitude( vectorPart.get_magnitude() );

    return quat4t<FloatType>( expReal * cos( vectorPartMagnitude ),
                              vectorPartMagnitude > 0
                                  ? vectorPart * expReal * sin( vectorPartMagnitude ) / vectorPartMagnitude
                                  : vector3t<FloatType>( 0 ) );
}

typedef quat4t<float> quat4f;
typedef quat4t<double> quat4fd;
} // namespace graphics
} // namespace frantic
