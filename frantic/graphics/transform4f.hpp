// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stdexcept>

#include <frantic/graphics/cubeface.hpp>
#include <frantic/graphics/graphics_utils.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace graphics {
// Forward Declare quaternion
// (quat4f.hpp already #includes this file)
template <class FloatType>
class quat4t;
} // namespace graphics
} // namespace frantic

#include <boost/lexical_cast.hpp>

#ifdef MAX_VERSION
#if MAX_RELEASE >= 25000
  #include <geom/Matrix3.h>
#else
  #include <Matrix3.h>
#endif
#endif

namespace frantic {
namespace graphics {

enum frustum_style { frustum_style_max, frustum_style_opengl };

// Class: transform4t<FloatType>/transform4f/transform4fd
// This class represents a homogeneous 3D transform that is compatible with
// OpenGL's native format in memory organization.
//> |  0  4  8  12  |   |  RS R  R  T  |
//> |  1  5  9  13  |   |  R  RS R  T  |
//> |  2  6  10 14  | = |  R  R  RS T  |
//> |  3  7  11 15  |   |  0  0  0  1  |
template <class FloatType>
class transform4t {

  protected:
    /// <summary>
    /// Represents a homogeneous 3D transform that is compatible with
    /// OpenGL's native format in memory organization.
    /// </summary>

    // |  0  4  8  12  |   |  RS R  R  T  |
    // |  1  5  9  13  |   |  R  RS R  T  |
    // |  2  6  10 14  | = |  R  R  RS T  |
    // |  3  7  11 15  |   |  0  0  0  1  |
    FloatType m_elements[16];

    // ================================================================================

  public:
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;
    typedef quat4t<FloatType> quat4f_type;

    // ================================================================================
    // Constructors
    // ================================================================================

    /**
     * This is the default constructor which sets this
     * transformation matrix to the identity matrix.
     */
    transform4t() { set_to_identity(); }

    /** Implicit conversion from smaller float type is allowed. */
    template <class OtherFloatType>
    transform4t( const transform4t<OtherFloatType>& t,
                 typename boost::enable_if_c<( sizeof( OtherFloatType ) < sizeof( FloatType ) ), int*>::type = 0 ) {
        const OtherFloatType* el = &t[0];
        for( int i = 0; i < 16; i++ )
            m_elements[i] = el[i];
    }

    /** Explicit cast required for conversion from larger float type. */
    template <class OtherFloatType>
    explicit transform4t(
        const transform4t<OtherFloatType>& t,
        typename boost::enable_if_c<( sizeof( OtherFloatType ) > sizeof( FloatType ) ), int*>::type = 0 ) {
        const OtherFloatType* el = &t[0];
        for( int i = 0; i < 16; i++ )
            m_elements[i] = (FloatType)el[i];
    }

    // Constructor: transform4
    // A constructor that creates the matrix and initializes its values
    //  from the float array specified.
    transform4t( const FloatType* elements ) {
        for( int i = 0; i < 16; i++ )
            m_elements[i] = elements[i];
    }

    // Constructor: transform4
    // A constructor that creates the matrix with the specified values.
    transform4t( FloatType e11, FloatType e21, FloatType e31, FloatType e41, FloatType e12, FloatType e22,
                 FloatType e32, FloatType e42, FloatType e13, FloatType e23, FloatType e33, FloatType e43,
                 FloatType e14, FloatType e24, FloatType e34, FloatType e44 ) {
        m_elements[0] = e11;
        m_elements[1] = e21;
        m_elements[2] = e31;
        m_elements[3] = e41;
        m_elements[4] = e12;
        m_elements[5] = e22;
        m_elements[6] = e32;
        m_elements[7] = e42;
        m_elements[8] = e13;
        m_elements[9] = e23;
        m_elements[10] = e33;
        m_elements[11] = e43;
        m_elements[12] = e14;
        m_elements[13] = e24;
        m_elements[14] = e34;
        m_elements[15] = e44;
    }

    // Constructor: transform4t
    // A constructor that creates the matrix with the columns initialized to the given column vectors.
    // The last entries in each column are set as in the identity.
    transform4t( const vector3f_type& column1, const vector3f_type& column2, const vector3f_type& column3,
                 const vector3f_type& column4 ) {
        m_elements[0] = column1.x;
        m_elements[1] = column1.y;
        m_elements[2] = column1.z;
        m_elements[3] = FloatType( 0.0 );
        m_elements[4] = column2.x;
        m_elements[5] = column2.y;
        m_elements[6] = column2.z;
        m_elements[7] = FloatType( 0.0 );
        m_elements[8] = column3.x;
        m_elements[9] = column3.y;
        m_elements[10] = column3.z;
        m_elements[11] = FloatType( 0.0 );
        m_elements[12] = column4.x;
        m_elements[13] = column4.y;
        m_elements[14] = column4.z;
        m_elements[15] = FloatType( 1.0 );
    }

    transform4t( const vector3f_type& column1, const vector3f_type& column2, const vector3f_type& column3 ) {
        m_elements[0] = column1.x;
        m_elements[1] = column1.y;
        m_elements[2] = column1.z;
        m_elements[3] = FloatType( 0.0 );
        m_elements[4] = column2.x;
        m_elements[5] = column2.y;
        m_elements[6] = column2.z;
        m_elements[7] = FloatType( 0.0 );
        m_elements[8] = column3.x;
        m_elements[9] = column3.y;
        m_elements[10] = column3.z;
        m_elements[11] = FloatType( 0.0 );
        m_elements[12] = FloatType( 0.0 );
        m_elements[13] = FloatType( 0.0 );
        m_elements[14] = FloatType( 0.0 );
        m_elements[15] = FloatType( 1.0 );
    }

    // 3ds Max compatibility
#ifdef MAX_VERSION
    // NOTE: Matrix3 "row"s are the same as transform3f "column"s.  We are using
    //       column vectors, while 3ds Max is using row vectors.
    transform4t( const Matrix3& maxmat ) {
        const Point3& column1 = maxmat[0];
        const Point3& column2 = maxmat[1];
        const Point3& column3 = maxmat[2];
        const Point3& column4 = maxmat[3];
        m_elements[0] = column1.x;
        m_elements[1] = column1.y;
        m_elements[2] = column1.z;
        m_elements[3] = 0.0f;
        m_elements[4] = column2.x;
        m_elements[5] = column2.y;
        m_elements[6] = column2.z;
        m_elements[7] = 0.0f;
        m_elements[8] = column3.x;
        m_elements[9] = column3.y;
        m_elements[10] = column3.z;
        m_elements[11] = 0.0f;
        m_elements[12] = column4.x;
        m_elements[13] = column4.y;
        m_elements[14] = column4.z;
        m_elements[15] = 1.0f;
    }
#endif

    // Function: zero
    // Gets the zero matrix, the matrix with all values of 0.
    //
    // Returns:
    // A transform4t matrix with all zero values.
    static transform4t zero() { return transform4t( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ); }

    // Function: identity
    // Gets the identity matrix, the matrix with ones in the diagonal and zero everywhere else.
    //
    // Returns:
    // A transform4t matrix with ones in the diagonal and zero everywhere else.
    static transform4t identity() { return transform4t( 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 ); }

    static transform4t from_scale( FloatType x, FloatType y, FloatType z ) {
        return transform4t( x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1 );
    }

    /**
     * Create a set of normalized basis vectors for a coordinate system.
     *
     * @param  normal		The normal from which to create the rotated coord system
     * @param  outNormal	The outNormal of the coord system.
     * @param  outBinormal	The outBinormal of the coord system.
     * @param  outTangent	The outTangent of the coord system.
     */
    static void orthonormal_system_from_normal( const vector3f_type& normal, vector3f_type& outNormal,
                                                vector3f_type& outBinormal, vector3f_type& outTangent ) {
        outNormal = normal;
        outNormal.normalize();

        outTangent = vector3f( -outNormal.y, outNormal.x, FloatType( 0.0 ) );
        // Normalize the tangent vector.
        // If its magnitude is zero, then set it to the x unit vector instead.
        const FloatType tangentMagnitude = outTangent.get_magnitude();
        if( tangentMagnitude > 0 ) {
            outTangent /= tangentMagnitude;
        } else {
            outTangent = vector3f_type( FloatType( 1 ), FloatType( 0 ), FloatType( 0 ) );
        }

        outBinormal = vector3f_type::cross( outTangent, outNormal );
    }

    ///// Constructors /////

    /**
     * Create a matrix with its z axis aligned to the given normal vector.
     *
     * @param  normal  The normal from which to create the rotated coord system
     */
    static transform4t from_normal_z( vector3f_type normal ) {

        vector3f_type newNormal, newBinormal, newTangent;
        orthonormal_system_from_normal( normal, newNormal, newBinormal, newTangent );

        return transform4t( newTangent, -newBinormal, newNormal );
    }

    /**
     * Create a matrix with its x axis aligned to the given normal vector.
     *
     * @param  normal  The normal from which to create the rotated coord system
     */
    static transform4t from_normal_x( vector3f_type normal ) {

        vector3f_type newNormal, newBinormal, newTangent;
        orthonormal_system_from_normal( normal, newNormal, newBinormal, newTangent );

        return transform4t( newNormal, newBinormal, newTangent );
    }

    /**
     * Create a matrix with its y axis aligned to the given normal vector.
     *
     * @param  normal  The normal from which to create the rotated coord system
     */
    static transform4t from_normal_y( vector3f_type normal ) {

        vector3f_type newNormal, newBinormal, newTangent;
        orthonormal_system_from_normal( normal, newNormal, newBinormal, newTangent );

        return transform4t( newTangent, newNormal, newBinormal );
    }

    /**
     * Creates a world to camera space matrix.  Negative z looks in the specified direction.
     * Positive y is constrained to be "up" as much as possible.
     * This is similar to gluLookAt from the OpenGL documentation.
     *
     * @param  pos  This is the position, in world space, where the camera will should be placed.
     * @param  at   This is the target at which the camera will be looking.
     * @param  up   This is the direction for "up".  Typically, up means towards positive Z.
     */
    static transform4t from_look_at( const vector3f_type& pos, const vector3f_type& at,
                                     const vector3f_type& up = vector3f_type::from_zaxis() ) {
        return from_look_dir( pos, vector3f_type::normalize( at - pos ), up );
    }

    /**
     * Creates a world to camera space matrix.  Negative z looks in the specified direction.
     * The input parameter dir is assumed to be normalized. Positive y is constrained to be "up" as much as possible.
     * This is similar to gluLookAt from the OpenGL documentation.
     *
     * @param  pos  This is the position, in world space, where the camera will should be placed.
     * @param  dir  This is the direction the camera will be looking.  ASSUMED TO BE NORMALIZED.
     * @param  up   This is the direction for "up".  Typically, up means towards positive Z.
     */
    static transform4t from_look_dir( const vector3f_type& pos, const vector3f_type& dir,
                                      const vector3f_type& up = vector3f_type::from_zaxis() ) {
        vector3f_type s = vector3f_type::cross( dir, up );
        s.normalize();

        vector3f_type u = vector3f_type::cross( s, dir );
        u.normalize();

        transform4t result( s.x, u.x, -dir.x, 0, s.y, u.y, -dir.y, 0, s.z, u.z, -dir.z, 0,
                            -vector3f_type::dot( s, pos ), -vector3f_type::dot( u, pos ),
                            vector3f_type::dot( dir, pos ), 1 );

        return result;
    }

    static transform4t from_scale( vector3f_type pt ) { return transform4t::from_scale( pt.x, pt.y, pt.z ); }

    static transform4t from_scale( FloatType scale ) { return transform4t::from_scale( scale, scale, scale ); }

    static transform4t from_translation( FloatType x, FloatType y, FloatType z ) {
        transform4t m;
        m.m_elements[12] = x;
        m.m_elements[13] = y;
        m.m_elements[14] = z;
        return m;
    }

    static transform4t from_translation( const vector3f_type& pt ) {
        return transform4t::from_translation( pt.x, pt.y, pt.z );
    }

    static transform4t from_transpose( const transform4t& a ) {
        return transform4t( a.m_elements[0], a.m_elements[4], a.m_elements[8], a.m_elements[12], a.m_elements[1],
                            a.m_elements[5], a.m_elements[9], a.m_elements[13], a.m_elements[2], a.m_elements[6],
                            a.m_elements[10], a.m_elements[14], a.m_elements[3], a.m_elements[7], a.m_elements[11],
                            a.m_elements[15] );
    }

    // These matrices assume that we are in the right handed coordinate system.
    // They also assume that the forward direction is facing -z to be consistent with 3dsmax cameras
    // To use these matrices, multiply a camera view
    static transform4t from_cubeface( cube_face::default_cube_face cubeFace ) {
        switch( cubeFace ) {
        case cube_face::CF_X_POS: // Looking towards (1,0,0)
            return transform4t( vector3f_type( 0, 0, 1 ), vector3f_type( 0, 1, 0 ), vector3f_type( -1, 0, 0 ),
                                vector3f_type( 0 ) );
            //  0  0 -1
            //  0  1  0
            //  1  0  0
        case cube_face::CF_X_NEG: // Looking towards (-1,0,0)
            return transform4t( vector3f_type( 0, 0, -1 ), vector3f_type( 0, 1, 0 ), vector3f_type( 1, 0, 0 ),
                                vector3f_type( 0 ) );
            //  0  0  1
            //  0  1  0
            // -1  0  0
        case cube_face::CF_Y_POS: // Looking towards (0,1,0)
            return transform4t( vector3f_type( 1, 0, 0 ), vector3f_type( 0, 0, 1 ), vector3f_type( 0, -1, 0 ),
                                vector3f_type( 0 ) );
            //  1  0  0
            //  0  0 -1
            //  0  1  0
        case cube_face::CF_Y_NEG: // Looking towards (0,-1,0)
            return transform4t( vector3f_type( 1, 0, 0 ), vector3f_type( 0, 0, -1 ), vector3f_type( 0, 1, 0 ),
                                vector3f_type( 0 ) );
            //  1  0  0
            //  0  0  1
            //  0 -1  0
        case cube_face::CF_Z_POS: // Looking towards (0,0,1)
            return transform4t( vector3f_type( -1, 0, 0 ), vector3f_type( 0, 1, 0 ), vector3f_type( 0, 0, -1 ),
                                vector3f_type( 0 ) );
            // -1  0  0
            //  0  1  0
            //  0  0 -1
        case cube_face::CF_Z_NEG: // Looking towards (0,0,-1)
            return transform4t( vector3f_type( 1, 0, 0 ), vector3f_type( 0, 1, 0 ), vector3f_type( 0, 0, 1 ),
                                vector3f_type( 0 ) );
            //  1  0  0
            //  0  1  0
            //  0  0  1
        default:
            throw std::runtime_error( "transform4t::from_cubeface: Requested invalid cube face index." );
        }
    }

    static transform4t from_axis_permutation( frantic::graphics::vector3 axisPermutation ) {
        const FloatType ZERO( 0.0 );
        const FloatType ONE( 1.0 );

        if( !( ( axisPermutation.x == 0 || axisPermutation.y == 0 || axisPermutation.z == 0 ) &&
               ( axisPermutation.x == 1 || axisPermutation.y == 1 || axisPermutation.z == 1 ) &&
               ( axisPermutation.x == 2 || axisPermutation.y == 2 || axisPermutation.z == 2 ) ) )
            throw std::runtime_error( "transform4t::from_axis_permutation: The permutation provided, " +
                                      axisPermutation.str() + ", is not valid as a permutation." );

        return transform4t( axisPermutation.x == 0 ? ONE : ZERO, axisPermutation.x == 1 ? ONE : ZERO,
                            axisPermutation.x == 2 ? ONE : ZERO, ZERO, axisPermutation.y == 0 ? ONE : ZERO,
                            axisPermutation.y == 1 ? ONE : ZERO, axisPermutation.y == 2 ? ONE : ZERO, ZERO,
                            axisPermutation.z == 0 ? ONE : ZERO, axisPermutation.z == 1 ? ONE : ZERO,
                            axisPermutation.z == 2 ? ONE : ZERO, ZERO, ZERO, ZERO, ZERO, ONE );
    }

    ///// Operations /////

    static transform4t linear_interpolate( const transform4t& a, const transform4t& b, FloatType t ) {
        if( t <= FloatType( 0.0 ) )
            return a;
        else if( t >= FloatType( 1.0 ) )
            return b;

        transform4t result;

        FloatType one_t = FloatType( 1.0 ) - t;
        for( int i = 0; i < 16; i++ )
            result[i] = one_t * a[i] + t * b[i];

        return result;
        ;
    }

    static transform4t from_frustum( FloatType left, FloatType right, FloatType bottom, FloatType top,
                                     FloatType nearDist, FloatType farDist, frustum_style frustumStyle ) {

        transform4t m;

        FloatType xInvDelta;
        FloatType yInvDelta;
        FloatType zInvDelta;
        switch( frustumStyle ) {
        case frustum_style_opengl:
            xInvDelta = FloatType( 1.0 ) / ( left - right );
            yInvDelta = FloatType( 1.0 ) / ( top - bottom );
            zInvDelta = FloatType( 1.0 ) / ( farDist - nearDist );
            break;
        case frustum_style_max:
            xInvDelta = FloatType( 1.0 ) / ( right - left );
            yInvDelta = FloatType( 1.0 ) / ( top - bottom );
            zInvDelta = FloatType( 1.0 ) / ( farDist - nearDist );
            break;
        default:
            throw std::runtime_error( "requested frustumStyle not handled" );
        }

        m[0] = ( FloatType( 2 ) * nearDist ) * xInvDelta;
        m[4] = FloatType( 0 );
        m[8] = ( right + left ) * xInvDelta;
        m[12] = FloatType( 0 );

        m[1] = FloatType( 0 );
        m[5] = ( FloatType( 2 ) * nearDist ) * yInvDelta;
        m[9] = ( top + bottom ) * yInvDelta;
        m[13] = FloatType( 0 );

        m[2] = FloatType( 0 );
        m[6] = FloatType( 0 );
        m[10] = -( farDist + nearDist ) * zInvDelta;
        m[14] = -( FloatType( 2 ) * farDist * nearDist ) * zInvDelta;

        m[3] = FloatType( 0 );
        m[7] = FloatType( 0 );
        m[11] = FloatType( -1 );
        m[15] = FloatType( 0 );

        return m;
    }

    static transform4t from_perspective( FloatType yFovDegrees, FloatType aspect, FloatType nearDist, FloatType farDist,
                                         frustum_style frustumStyle = frustum_style_max ) {

        FloatType top =
            (FloatType)( tan( frantic::math::degrees_to_radians( yFovDegrees ) / FloatType( 2 ) ) * nearDist );
        FloatType bottom = -top;
        FloatType left = aspect * top;
        FloatType right = aspect * bottom;

        return from_frustum( left, right, bottom, top, nearDist, farDist, frustumStyle );
    }

    static transform4t from_fov_perspective( FloatType xFov, FloatType yFov, FloatType nearDist, FloatType farDist,
                                             frustum_style frustumStyle = frustum_style_max ) {

        FloatType top = (FloatType)( tan( yFov / 2 ) * nearDist );
        FloatType bottom = -top;
        FloatType left = -(FloatType)( tan( xFov / 2 ) * nearDist );
        FloatType right = -left;

        return from_frustum( left, right, bottom, top, nearDist, farDist, frustumStyle );
    }

    static transform4t from_x_rotation( FloatType thetaRadians ) {
        transform4t mat; // Defaults to the identity
        FloatType sinTheta = sin( thetaRadians ), cosTheta = cos( thetaRadians );
        mat.set( 1, 1, cosTheta );
        mat.set( 2, 1, -sinTheta );
        mat.set( 1, 2, sinTheta );
        mat.set( 2, 2, cosTheta );
        return mat;
    }

    static transform4t from_y_rotation( FloatType thetaRadians ) {
        transform4t mat; // Defaults to the identity
        FloatType sinTheta = sin( thetaRadians ), cosTheta = cos( thetaRadians );
        mat.set( 0, 0, cosTheta );
        mat.set( 2, 0, sinTheta );
        mat.set( 0, 2, -sinTheta );
        mat.set( 2, 2, cosTheta );
        return mat;
    }

    // mw: I think this one may be rotating the wrong way around Z
    // Brian: I think the negative on the sin values were backward
    // mw: switched it back, the unit tests caught this change. Need solid mathematical rationale if it is to be
    // changed.
    static transform4t from_z_rotation( FloatType thetaRadians ) {
        transform4t mat; // Defaults to the identity
        FloatType sinTheta = sin( thetaRadians ), cosTheta = cos( thetaRadians );
        mat.set( 0, 0, cosTheta );
        mat.set( 1, 0, sinTheta );
        mat.set( 0, 1, -sinTheta );
        mat.set( 1, 1, cosTheta );
        return mat;
    }

    static transform4t from_euler_xyz_rotation( const vector3f_type& angles ) {
        return from_z_rotation( angles.z ) * from_y_rotation( angles.y ) * from_x_rotation( angles.x );
    }

    static transform4t from_euler_xyz_rotation_degrees( const vector3f_type& angles ) {
        return from_z_rotation( frantic::math::degrees_to_radians( angles.z ) ) *
               from_y_rotation( frantic::math::degrees_to_radians( angles.y ) ) *
               from_x_rotation( frantic::math::degrees_to_radians( angles.x ) );
    }

    /**
     * Creates a transformation that rotates a vector by a given number of radians (angleRads)
     * around a given axis (normalAxis).
     *
     * @param  angleRads  The angle in radians to rotate through.
     * @param  normalAxis The NORMALIZED axis about which to rotate.
     */
    static transform4t from_angle_axis( FloatType angleRads, const vector3f_type& normalAxis ) {
        const FloatType ONE( 1.0 );
        const FloatType ZERO( 0.0 );
        FloatType s = sin( angleRads );
        FloatType c = cos( angleRads );
        FloatType x = normalAxis.x, y = normalAxis.y, z = normalAxis.z;
        return transform4t( ONE + ( ONE - c ) * ( x * x - ONE ), -z * s + ( ONE - c ) * x * y,
                            y * s + ( ONE - c ) * x * z, ZERO, z * s + ( ONE - c ) * x * y,
                            ONE + ( ONE - c ) * ( y * y - ONE ), -x * s + ( ONE - c ) * y * z, ZERO,
                            -y * s + ( ONE - c ) * x * z, x * s + ( ONE - c ) * y * z,
                            ONE + ( ONE - c ) * ( z * z - ONE ), ZERO, ZERO, ZERO, ZERO, ONE );
    }

    /// <summary>
    /// An OpenGL compatible linear indexer to the matrix elements
    /// |  0  4  8  12  |   |  RS R  R  T  |
    /// |  1  5  9  13  |   |  R  RS R  T  |
    /// |  2  6  10 14  | = |  R  R  RS T  |
    /// |  3  7  11 15  |   |  0  0  0  1  |
    /// </summary>
    FloatType& operator[]( int i ) {
        if( i < 0 || i > 15 )
            throw std::runtime_error( "transform4t: operator[] index must be in the range [0,15]" );

        return m_elements[i];
    }

    const FloatType& operator[]( int i ) const {
        if( i < 0 || i > 15 )
            throw std::runtime_error( "transform4t: operator[] index must be in the range [0,15]" );

        return m_elements[i];
    }

    FloatType& operator[]( size_t i ) {
        if( i > 15 )
            throw std::runtime_error( "transform4t: operator[] index must be in the range [0,15]" );

        return m_elements[i];
    }

    const FloatType& operator[]( size_t i ) const {
        if( i > 15 )
            throw std::runtime_error( "transform4t: operator[] index must be in the range [0,15]" );

        return m_elements[i];
    }

    // Function: get
    // Gets the value at row "row" and column "column"
    //
    // Parameters:
    // column	- the column the value is in.
    // row		- the row the value is in.
    //
    // Returns:
    // The float that is at (row, column).
    FloatType get( int column, int row ) const { return m_elements[row + column * 4]; }

    FloatType& get( int column, int row ) { return m_elements[row + column * 4]; }

    const vector3f_type& get_column( int column ) const {
        // column *= 4;
        // return vector3f_type( m_elements[column], m_elements[column+1], m_elements[column+2] );
        return *reinterpret_cast<const vector3f_type*>( &m_elements[column << 2] );
    }

    vector3f_type& get_column( int column ) {
        // column *= 4;
        // return vector3f_type( m_elements[column], m_elements[column+1], m_elements[column+2] );
        return *reinterpret_cast<vector3f_type*>( &m_elements[column << 2] );
    }

    vector3f_type get_row( int row ) const {
        return vector3f_type( m_elements[row], m_elements[row + 4], m_elements[row + 8] );
    }

    /**
     * Returns the euler angles from the transformation matrix.
     *
     * @note For all matrices except for the ones with theta values pi/2 and -pi/2, there are 2 correct solutions. To
     * break the ambiguity, we restrict the domain of theta values to be [-pi/2,pi/2]. Any theta value which does not
     * lie in this domain is not considered and thus the corresponding values of the other 2 angles are not considered
     * either.
     * @note The negative sign for the calculation of the z angle is due to the fact that the from_z_rotation function
     * has a tansposed matrix as compared to the one in the paper
     */
    vector3f_type get_euler_angles() const {
        if( m_elements[2] == -1 ) {
            FloatType theta = M_PI_2;
            FloatType phi = 0;
            FloatType psi = atan2( m_elements[4], m_elements[8] );
            return vector3f_type( psi, theta, phi );
        } else if( m_elements[2] == 1 ) {
            FloatType theta = -M_PI_2;
            FloatType phi = 0;
            FloatType psi = atan2( -m_elements[4], -m_elements[8] );
            return vector3f_type( psi, theta, phi );
        } else {
            FloatType theta_1 = -asin( m_elements[2] );
            FloatType theta_2 = M_PI - theta_1;
            FloatType psi_1 = atan2( m_elements[6] / cos( theta_1 ), m_elements[10] / cos( theta_1 ) );
            FloatType psi_2 = atan2( m_elements[6] / cos( theta_2 ), m_elements[10] / cos( theta_2 ) );
            FloatType phi_1 = atan2( -m_elements[1] / cos( theta_1 ), m_elements[0] / cos( theta_1 ) );
            FloatType phi_2 = atan2( -m_elements[1] / cos( theta_2 ), m_elements[0] / cos( theta_2 ) );

            if( theta_1 > M_PI_2 || theta_1 < -M_PI_2 ) {
                return vector3f_type( psi_2, theta_2, phi_2 );
            } else {
                return vector3f_type( psi_1, theta_1, phi_1 );
            }
        }
    }

    // Function: set
    // Sets the value at the specified row and column to the
    //  specified value.
    //
    // Parameters:
    // column	- the column of the value to set.
    // row		- the row of the value to set.
    // value	- the new value at (row, column).
    void set( int column, int row, FloatType value ) { m_elements[row + column * 4] = value; }

    // Function: set_to_zero
    // Zeros this transformation matrix, that is, all values
    //  of this matrix after this operation will be zero.
    void set_to_zero() {
        for( int i = 0; i < 16; ++i )
            m_elements[i] = FloatType( 0.0 );
    }

    // Function: set_to_identity
    // Sets this transformation matrix to the identity matrix,
    // that is, all values along the (top-left to bottom-right) diagonal
    // will be 1, all other values will be 0.
    void set_to_identity() {
        const FloatType ZERO( 0.0 );
        const FloatType ONE( 1.0 );
        m_elements[0] = ONE;
        m_elements[1] = ZERO;
        m_elements[2] = ZERO;
        m_elements[3] = ZERO;
        m_elements[4] = ZERO;
        m_elements[5] = ONE;
        m_elements[6] = ZERO;
        m_elements[7] = ZERO;
        m_elements[8] = ZERO;
        m_elements[9] = ZERO;
        m_elements[10] = ONE;
        m_elements[11] = ZERO;
        m_elements[12] = ZERO;
        m_elements[13] = ZERO;
        m_elements[14] = ZERO;
        m_elements[15] = ONE;
    }

    void set_to_scale( FloatType scale ) {
        const FloatType ZERO( 0.0 );
        m_elements[0] = scale;
        m_elements[1] = ZERO;
        m_elements[2] = ZERO;
        m_elements[3] = ZERO;
        m_elements[4] = ZERO;
        m_elements[5] = scale;
        m_elements[6] = ZERO;
        m_elements[7] = ZERO;
        m_elements[8] = ZERO;
        m_elements[9] = ZERO;
        m_elements[10] = scale;
        m_elements[11] = ZERO;
        m_elements[12] = ZERO;
        m_elements[13] = ZERO;
        m_elements[14] = ZERO;
        m_elements[15] = FloatType( 1.0 );
    }

    // Function: set_to_array
    // Sets the values of this array from the
    //  specified array.
    //
    // Parameters:
    // elements - a float array to set this tranformation matrix to.
    void set_to_array( FloatType* elements ) {
        for( int i = 0; i < 16; ++i )
            m_elements[i] = elements[i];
    }

    //---------------------------------------------------------------------------------------------------------

    bool equals( const transform4t& a ) const {
        for( int i = 0; i < 16; ++i ) {
            if( a.m_elements[i] != m_elements[i] )
                return false;
        }
        return true;
    }

    bool equals( const transform4t& a, FloatType tolerance ) const {
        FloatType errorSquared = 0;
        for( int i = 0; i < 16; ++i ) {
            errorSquared += frantic::math::square( m_elements[i] - a.m_elements[i] );
        }
        return errorSquared <= tolerance * tolerance;
    }

    bool equals_relative( const transform4t& a, FloatType relativeTolerance ) const {
        FloatType errorSquared = 0, absSumSquared = 0;
        for( int i = 0; i < 16; ++i ) {
            errorSquared += frantic::math::square( m_elements[i] - a.m_elements[i] );
            absSumSquared += frantic::math::square( m_elements[i] );
            absSumSquared += frantic::math::square( a.m_elements[i] );
        }
        // std::cout << "relative error " << sqrtf(errorSquared/absSumSquared) << std::endl;
        return errorSquared <= absSumSquared * relativeTolerance * relativeTolerance;
    }

    // Function: ==
    // Compares the specified transform matrix to this matrix and
    //  determines if they have exactly the
    //  same values.
    //
    // Parameters:
    // a - the transformation matrix to compare this matrix to.
    //
    // Returns:
    // True if the specified matrix is the same as this matrix.  False otherwise.
    bool operator==( const transform4t& a ) const { return equals( a ); }

    // Function: !=
    // Compares the specified transform matrix to this matrix
    //  and determines if they are different, that is, not equal.
    //
    // Parameters:
    // a - the transformation matrix to compare this matrix to.
    //
    // Returns:
    // true if the specified matrix is not the same as this matrix, false otherwise.
    bool operator!=( const transform4t& a ) const {
        for( int i = 0; i < 16; ++i )
            if( a.m_elements[i] != m_elements[i] )
                return true;
        return false;
    }

    bool is_zero() const {
        for( unsigned i = 0; i < 16; ++i ) {
            if( m_elements[i] != FloatType( 0.0 ) )
                return false;
        }
        return true;
    }

    bool is_zero( FloatType tolerance ) const {
        FloatType errorSquared = 0;
        for( int i = 0; i < 16; ++i ) {
            errorSquared += frantic::math::square( m_elements[i] );
        }
        return errorSquared <= tolerance * tolerance;
    }

    // Whether this matrix is the identity
    bool is_identity() const {
        const FloatType ZERO( 0.0 );
        const FloatType ONE( 1.0 );
        return m_elements[0] == ONE && m_elements[1] == ZERO && m_elements[2] == ZERO && m_elements[3] == ZERO &&
               m_elements[4] == ZERO && m_elements[5] == ONE && m_elements[6] == ZERO && m_elements[7] == ZERO &&
               m_elements[8] == ZERO && m_elements[9] == ZERO && m_elements[10] == ONE && m_elements[11] == ZERO &&
               m_elements[12] == ZERO && m_elements[13] == ZERO && m_elements[14] == ZERO && m_elements[15] == ONE;
    }

    // Whether this matrix is the identity, within a given tolerance
    bool is_identity( FloatType tolerance ) const {
        const FloatType ZERO( 0.0 );
        const FloatType ONE( 1.0 );
        FloatType errorSquared( 0.0 );

        errorSquared += frantic::math::square( m_elements[0] - ONE );
        errorSquared += frantic::math::square( m_elements[1] - ZERO );
        errorSquared += frantic::math::square( m_elements[2] - ZERO );
        errorSquared += frantic::math::square( m_elements[3] - ZERO );

        errorSquared += frantic::math::square( m_elements[4] - ZERO );
        errorSquared += frantic::math::square( m_elements[5] - ONE );
        errorSquared += frantic::math::square( m_elements[6] - ZERO );
        errorSquared += frantic::math::square( m_elements[7] - ZERO );

        errorSquared += frantic::math::square( m_elements[8] - ZERO );
        errorSquared += frantic::math::square( m_elements[9] - ZERO );
        errorSquared += frantic::math::square( m_elements[10] - ONE );
        errorSquared += frantic::math::square( m_elements[11] - ZERO );

        errorSquared += frantic::math::square( m_elements[12] - ZERO );
        errorSquared += frantic::math::square( m_elements[13] - ZERO );
        errorSquared += frantic::math::square( m_elements[14] - ZERO );
        errorSquared += frantic::math::square( m_elements[15] - ONE );

        //		std::cerr << "is_identity error term: " << sqrt(errorSquared) << std::endl;

        return errorSquared <= tolerance * tolerance;
    }

    // Whether this matrix is orthogonal (A * transpose(A) == I)
    bool is_orthogonal() const;
    bool is_orthogonal( FloatType tolerance ) const;

    // Whether this matrix is symmetric (A == transpose(A))
    bool is_symmetric() const;
    bool is_symmetric( FloatType tolerance ) const;

    // Whether this matrix does some kind of reflection of space
    bool is_orientation_inverting() const { return determinant() < FloatType( 0 ); }

    // Whether the matrix has any perspective transform effect, which
    // is indicated by the bottom row not being [0,0,0,1]
    bool is_perspective() const {
        return m_elements[3] != FloatType( 0 ) || m_elements[7] != FloatType( 0 ) || m_elements[11] != FloatType( 0 ) ||
               m_elements[15] != FloatType( 1 );
    }

    //---------------------------------------------------------------------------------------------------------

    // Decomposition of the matrix into perspective, translation, rotation and stretch.
    // Returns true if the decomposition was successful, false otherwise.
    //   The perspective matrix looks like this:
    //     [ 1 0 0 0   ]
    //     [ 0 1 0 0   ]
    //     [ 0 0 1 0   ]
    //     [ x y z w+1 ]
    //   The translation matrix (returned as a vector) looks like this:
    //     [ 1 0 0 x ]
    //     [ 0 1 0 y ]
    //     [ 0 0 1 z ]
    //     [ 0 0 0 1 ]
    //
    // To reconstruct the original matrix, use: outPerspective * transform4t::from_translation(outTranslation) *
    // outRotation * outStretch
    void decompose( transform4t& outPerspective, vector3f_type& outTranslation, transform4t& outRotation,
                    transform4t& outStretch ) const;

    // This does a polar decomposition of the matrix.  If the matrix is not a perspective
    // transform, and has no translation component, outQ and outS are a "rotation" and a "stretch".
    // This uses the method described by the paper linked to from the decompose method,
    // which is credited as [Higham 86] there.  It apparently has quadratic convergence
    // when Qi is nearly orthogonal
    void decompose_polar( transform4t& outQ, transform4t& outS, FloatType tolerance ) const {
        transform4t Qprev;
        outQ = *this;
        int maximumIterations = 200;
        do {
            Qprev = outQ;
            transform4t QinverseTranspose = Qprev.to_inverse();
            QinverseTranspose.transpose();
            outQ = Qprev + QinverseTranspose;
            outQ *= 0.5f;
        } while( !Qprev.equals( outQ, tolerance ) && --maximumIterations > 0 );
        outS = outQ.to_inverse() * ( *this );

        // Note that while polar decomposition gets us an orthogonal matrix, it is not necessarily the rotation matrix
        // This orthogonal matrix may contain a reflection component which should go with the scale matrix instead
        // (mathematically, polar decomposition gets us an orthogonal matrix but not necessarily a "special" orthogonal
        // matrix)
        // TODO: is there a better work around for this?
        if( outQ.is_orientation_inverting() ) {
            transform4t reflection;

            // Try reflecting over the "most negative" axis
            // It shouldn't matter which axis we reflect.
            // If we reflect the "wrong" one, the rotation matrix will compensate by including a 180 degree rotation.
            const float_type a00 = outS.get( 0, 0 );
            const float_type a11 = outS.get( 1, 1 );
            const float_type a22 = outS.get( 2, 2 );
            const float_type a33 = outS.get( 3, 3 );
            if( a00 < a11 && a00 < a22 && a00 < a33 ) {
                reflection.set( 0, 0, -1 );
            } else if( a11 < a22 && a11 < a33 ) {
                reflection.set( 1, 1, -1 );
            } else if( a22 < a33 ) {
                reflection.set( 2, 2, -1 );
            } else {
                reflection.set( 3, 3, -1 );
            }

            outS = reflection * outS;
            outQ = outQ * reflection;
        }
    }

    /**
     * Decomposes the matrix into a translation, a rotation, and x-y-z scale component
     * To reconstruct the original matrix, use: transform4t::from_translation(outTranslation) *
     * outRotation.to_transform4f() * transform4t::from_scale(outScale) If the transform contains a perspective
     * component, it is dropped.
     */
    void decompose( vector3f_type& outTranslation, quat4f_type& outRotation, vector3f_type& outScale ) const {
        transform4t perspective;
        transform4t rotation;
        transform4t stretch;
        decompose( perspective, outTranslation, rotation, stretch );
        outRotation = quat4f_type::from_transform4f( rotation );

        outScale.x = stretch.get( 0, 0 );
        outScale.y = stretch.get( 1, 1 );
        outScale.z = stretch.get( 2, 2 );
    }

    /**
     * Concatenates the given translation, rotation, and scale into a transform4f
     */
    static transform4t to_transform( const vector3f_type& translation, const quat4f_type& rotation,
                                     const vector3f_type& scale ) {
        return from_translation( translation ) * rotation.to_transform4f() * from_scale( scale );
    }

    /**
     * Decomposes the matrix into a translation and a rotation component
     * To reconstruct the original matrix, use: transform4t::from_translation(outTranslation) *
     * outRotation.to_transform4f() If the transform contains a perspective or scale component, they are dropped.
     */
    void decompose( vector3f_type& outTranslation, quat4f_type& outRotation ) const {
        vector3f_type scale;
        decompose( outTranslation, outRotation, scale );
    }

    /**
     * Concatenates the given translation and rotation into a transform4f
     */
    static transform4t to_transform( const vector3f_type& translation, const quat4f_type& rotation ) {
        return from_translation( translation ) * rotation.to_transform4f();
    }

    // Does an LU decomposition.  When it's done, (*this) == outL*outU,
    // L is lower triangular and U is upper triangular.
    // void decompose_lu( transform4t& outLU ) const {
    //	outU.set_to_identity();
    //	outL = *this;
    //}

    // One step of the polar decomposition iteration can be used to renormalize a rotation matrix.
    // This only works for matrices which are already close to a rotation matrix.
    void renormalize_as_rotation() {
        transform4t inverseTranspose = to_inverse();
        inverseTranspose.transpose();
        ( *this ) += inverseTranspose;
        ( *this ) *= 0.5f;
    }

    //---------------------------------------------------------------------------------------------------------

    void set_translation( vector3f_type translation ) {
        m_elements[12] = translation.x;
        m_elements[13] = translation.y;
        m_elements[14] = translation.z;
    }

    void translate( const vector3f_type& v ) {
        m_elements[12] += v.x;
        m_elements[13] += v.y;
        m_elements[14] += v.z;
    }

    vector3f_type translation() const { return vector3f_type( m_elements[12], m_elements[13], m_elements[14] ); }

    // Effect is result = transform4t::from_scale( scale,scale,scale ) * result;
    void scale( FloatType scale ) {
        m_elements[0] *= scale;
        m_elements[1] *= scale;
        m_elements[2] *= scale;
        m_elements[4] *= scale;
        m_elements[5] *= scale;
        m_elements[6] *= scale;
        m_elements[8] *= scale;
        m_elements[9] *= scale;
        m_elements[10] *= scale;
        m_elements[12] *= scale;
        m_elements[13] *= scale;
        m_elements[14] *= scale;
    }

    // Effect is result = transform4t::from_scale( scaleX,scaleY,scaleZ ) * result;
    void scale( FloatType scaleX, FloatType scaleY, FloatType scaleZ ) {
        m_elements[0] *= scaleX;
        m_elements[1] *= scaleY;
        m_elements[2] *= scaleZ;
        m_elements[4] *= scaleX;
        m_elements[5] *= scaleY;
        m_elements[6] *= scaleZ;
        m_elements[8] *= scaleX;
        m_elements[9] *= scaleY;
        m_elements[10] *= scaleZ;
        m_elements[12] *= scaleX;
        m_elements[13] *= scaleY;
        m_elements[14] *= scaleZ;
    }

    //---------------------------------------------------------------------------------------------------------

    FloatType max_abs_component() const {
        FloatType result = 0;
        for( unsigned i = 0; i < 16; ++i ) {
            FloatType absComponent = fabs( m_elements[i] );
            if( absComponent > result )
                result = absComponent;
        }
        return result;
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------

    /// <summary>
    /// transpose the matrix
    /// </summary>
    void transpose() {
        std::swap( m_elements[1], m_elements[4] );
        std::swap( m_elements[2], m_elements[8] );
        std::swap( m_elements[3], m_elements[12] );
        std::swap( m_elements[6], m_elements[9] );
        std::swap( m_elements[7], m_elements[13] );
        std::swap( m_elements[11], m_elements[14] );
    }

    transform4t to_transpose() const {
        transform4t result( *this );
        result.transpose();
        return result;
    }

    transform4t to_translation( const vector3f_type& v ) {
        transform4t result( *this );
        result.set_translation( v );
        return result;
    }

    transform4t to_scale_free() const {
        transform4t perspective, rotation, stretch;
        vector3f_type translation;
        decompose( perspective, translation, rotation, stretch );
        return perspective * transform4t::from_translation( translation ) * rotation;
    }

    //
    transform4t to_scale_free_not_skew() const {
        FloatType col0M =
            sqrt( m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1] + m_elements[2] * m_elements[2] );
        FloatType col1M =
            sqrt( m_elements[4] * m_elements[4] + m_elements[5] * m_elements[5] + m_elements[6] * m_elements[6] );
        FloatType col2M =
            sqrt( m_elements[8] * m_elements[8] + m_elements[9] * m_elements[9] + m_elements[10] * m_elements[10] );

        return transform4t( m_elements[0] / col0M, m_elements[1] / col0M, m_elements[2] / col0M, m_elements[3],
                            m_elements[4] / col1M, m_elements[5] / col1M, m_elements[6] / col1M, m_elements[7],
                            m_elements[8] / col2M, m_elements[9] / col2M, m_elements[10] / col2M, m_elements[11],
                            m_elements[12], m_elements[13], m_elements[14], m_elements[15] );
    }

    // Function: determinant
    // Calculates the determinant of this matrix.

    // e[0]  e[4]  e[8]  e[12]
    // e[1]  e[5]  e[9]  e[13]
    // e[2]  e[6]  e[10] e[14]
    // e[3]  e[7]  e[11] e[15]

    FloatType determinant() const {
        FloatType result( 0 );
        // Expand along the bottom row, because that will typically be [0,0,0,1]
        if( m_elements[3] != FloatType( 0 ) ) {
            result -= m_elements[3] * det3x3( m_elements[4], m_elements[5], m_elements[6], m_elements[8], m_elements[9],
                                              m_elements[10], m_elements[12], m_elements[13], m_elements[14] );
        }
        if( m_elements[7] != FloatType( 0 ) ) {
            result += m_elements[7] * det3x3( m_elements[0], m_elements[1], m_elements[2], m_elements[8], m_elements[9],
                                              m_elements[10], m_elements[12], m_elements[13], m_elements[14] );
        }
        if( m_elements[11] != FloatType( 0 ) ) {
            result -= m_elements[7] * det3x3( m_elements[0], m_elements[1], m_elements[2], m_elements[4], m_elements[5],
                                              m_elements[6], m_elements[12], m_elements[13], m_elements[14] );
        }
        if( m_elements[15] != FloatType( 0 ) ) {
            result +=
                m_elements[15] * det3x3( m_elements[0], m_elements[1], m_elements[2], m_elements[4], m_elements[5],
                                         m_elements[6], m_elements[8], m_elements[9], m_elements[10] );
        }
        return result;
    }

    //
    // Returns:
    // The determinant of this matrix as a double.
    double determinant_double() const {
        double result = 0.0;
        // Expand along the bottom row, because that will typically be [0,0,0,1]
        if( m_elements[3] != FloatType( 0 ) ) {
            result -= m_elements[3] * det3x3_double( m_elements[4], m_elements[5], m_elements[6], m_elements[8],
                                                     m_elements[9], m_elements[10], m_elements[12], m_elements[13],
                                                     m_elements[14] );
        }
        if( m_elements[7] != FloatType( 0 ) ) {
            result += m_elements[7] * det3x3_double( m_elements[0], m_elements[1], m_elements[2], m_elements[8],
                                                     m_elements[9], m_elements[10], m_elements[12], m_elements[13],
                                                     m_elements[14] );
        }
        if( m_elements[11] != FloatType( 0 ) ) {
            result -= m_elements[7] * det3x3_double( m_elements[0], m_elements[1], m_elements[2], m_elements[4],
                                                     m_elements[5], m_elements[6], m_elements[12], m_elements[13],
                                                     m_elements[14] );
        }
        if( m_elements[15] != FloatType( 0 ) ) {
            result += m_elements[15] * det3x3_double( m_elements[0], m_elements[1], m_elements[2], m_elements[4],
                                                      m_elements[5], m_elements[6], m_elements[8], m_elements[9],
                                                      m_elements[10] );
        }
        return result;
    }

    // Function: det3x3
    // Determines the determinant of a 3x3 matrix formed by 9 values. That is, the determinant of:
    // >		a11 a12 a13      e[0]  e[4]  e[8]
    // >		a21 a22 a23      e[1]  e[5]  e[9]
    // >		a31 a32 a33      e[2]  e[6]  e[10]
    //
    // Returns:
    // The determinant of the matrix formed by the specified values, as a float.

    static FloatType det3x3( FloatType a11, FloatType a21, FloatType a31, FloatType a12, FloatType a22, FloatType a32,
                             FloatType a13, FloatType a23, FloatType a33 ) {
        return a11 * ( a22 * a33 - a32 * a23 ) - a21 * ( a12 * a33 - a32 * a13 ) + a31 * ( a12 * a23 - a22 * a13 );
    }

    static double det3x3_double( double a11, double a21, double a31, double a12, double a22, double a32, double a13,
                                 double a23, double a33 ) {
        return a11 * ( a22 * a33 - a32 * a23 ) - a21 * ( a12 * a33 - a32 * a13 ) + a31 * ( a12 * a23 - a22 * a13 );
    }

    FloatType determinant_3x3() const {
        return m_elements[0] * ( m_elements[5] * m_elements[10] - m_elements[6] * m_elements[9] ) -
               m_elements[1] * ( m_elements[4] * m_elements[10] - m_elements[6] * m_elements[8] ) +
               m_elements[2] * ( m_elements[4] * m_elements[9] - m_elements[5] * m_elements[8] );
    }

    // Function: adjoint
    // Finds the adjoint of this transform matrix.
    //
    // Returns:
    // Another transformation matrix, which is the adjoint of this matrix.
    transform4t adjoint() const {
        transform4t result;
        FloatType a1, a2, a3, a4, b1, b2, b3, b4;
        FloatType c1, c2, c3, c4, d1, d2, d3, d4;

        /* assign to individual variable names to aid  */
        /* selecting correct values  */

        a1 = get( 0, 0 );
        b1 = get( 0, 1 );
        c1 = get( 0, 2 );
        d1 = get( 0, 3 );

        a2 = get( 1, 0 );
        b2 = get( 1, 1 );
        c2 = get( 1, 2 );
        d2 = get( 1, 3 );

        a3 = get( 2, 0 );
        b3 = get( 2, 1 );
        c3 = get( 2, 2 );
        d3 = get( 2, 3 );

        a4 = get( 3, 0 );
        b4 = get( 3, 1 );
        c4 = get( 3, 2 );
        d4 = get( 3, 3 );

        /* row column labeling reversed since we transpose rows & columns */

        result.set( 0, 0, det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4 ) );
        result.set( 1, 0, -det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4 ) );
        result.set( 2, 0, det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4 ) );
        result.set( 3, 0, -det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4 ) );

        result.set( 0, 1, -det3x3( b1, b3, b4, c1, c3, c4, d1, d3, d4 ) );
        result.set( 1, 1, det3x3( a1, a3, a4, c1, c3, c4, d1, d3, d4 ) );
        result.set( 2, 1, -det3x3( a1, a3, a4, b1, b3, b4, d1, d3, d4 ) );
        result.set( 3, 1, det3x3( a1, a3, a4, b1, b3, b4, c1, c3, c4 ) );

        result.set( 0, 2, det3x3( b1, b2, b4, c1, c2, c4, d1, d2, d4 ) );
        result.set( 1, 2, -det3x3( a1, a2, a4, c1, c2, c4, d1, d2, d4 ) );
        result.set( 2, 2, det3x3( a1, a2, a4, b1, b2, b4, d1, d2, d4 ) );
        result.set( 3, 2, -det3x3( a1, a2, a4, b1, b2, b4, c1, c2, c4 ) );

        result.set( 0, 3, -det3x3( b1, b2, b3, c1, c2, c3, d1, d2, d3 ) );
        result.set( 1, 3, det3x3( a1, a2, a3, c1, c2, c3, d1, d2, d3 ) );
        result.set( 2, 3, -det3x3( a1, a2, a3, b1, b2, b3, d1, d2, d3 ) );
        result.set( 3, 3, det3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3 ) );

        return result;
    }

    // ================================================================================
    // Function: *
    // Performs the matrix multiplication of [this Matrix] x [Matrix b]
    //
    // Parameters:
    // b - the matrix to multiply this matrix to.
    //
    // Results:
    // A transformation matrix which is the matrix multiplication of this and b.
    transform4t operator*( const transform4t& b ) const {
        transform4t mResult;

        // make things easier to read
        const FloatType* e0 = &this->m_elements[0];
        const FloatType* e1 = &b.m_elements[0];
        FloatType* result = &mResult.m_elements[0];
        if( !is_perspective() && !b.is_perspective() ) {
            // Special case the affine transform, because it is common.
            result[0] = e0[0] * e1[0] + e0[4] * e1[1] + e0[8] * e1[2];
            result[4] = e0[0] * e1[4] + e0[4] * e1[5] + e0[8] * e1[6];
            result[8] = e0[0] * e1[8] + e0[4] * e1[9] + e0[8] * e1[10];
            result[12] = e0[0] * e1[12] + e0[4] * e1[13] + e0[8] * e1[14] + e0[12];

            result[1] = e0[1] * e1[0] + e0[5] * e1[1] + e0[9] * e1[2];
            result[5] = e0[1] * e1[4] + e0[5] * e1[5] + e0[9] * e1[6];
            result[9] = e0[1] * e1[8] + e0[5] * e1[9] + e0[9] * e1[10];
            result[13] = e0[1] * e1[12] + e0[5] * e1[13] + e0[9] * e1[14] + e0[13];

            result[2] = e0[2] * e1[0] + e0[6] * e1[1] + e0[10] * e1[2];
            result[6] = e0[2] * e1[4] + e0[6] * e1[5] + e0[10] * e1[6];
            result[10] = e0[2] * e1[8] + e0[6] * e1[9] + e0[10] * e1[10];
            result[14] = e0[2] * e1[12] + e0[6] * e1[13] + e0[10] * e1[14] + e0[14];

            result[3] = 0;
            result[7] = 0;
            result[11] = 0;
            result[15] = 1;
        } else {
            // Do a full 4x4 matrix multiplication
            result[0] = e0[0] * e1[0] + e0[4] * e1[1] + e0[8] * e1[2] + e0[12] * e1[3];
            result[4] = e0[0] * e1[4] + e0[4] * e1[5] + e0[8] * e1[6] + e0[12] * e1[7];
            result[8] = e0[0] * e1[8] + e0[4] * e1[9] + e0[8] * e1[10] + e0[12] * e1[11];
            result[12] = e0[0] * e1[12] + e0[4] * e1[13] + e0[8] * e1[14] + e0[12] * e1[15];

            result[1] = e0[1] * e1[0] + e0[5] * e1[1] + e0[9] * e1[2] + e0[13] * e1[3];
            result[5] = e0[1] * e1[4] + e0[5] * e1[5] + e0[9] * e1[6] + e0[13] * e1[7];
            result[9] = e0[1] * e1[8] + e0[5] * e1[9] + e0[9] * e1[10] + e0[13] * e1[11];
            result[13] = e0[1] * e1[12] + e0[5] * e1[13] + e0[9] * e1[14] + e0[13] * e1[15];

            result[2] = e0[2] * e1[0] + e0[6] * e1[1] + e0[10] * e1[2] + e0[14] * e1[3];
            result[6] = e0[2] * e1[4] + e0[6] * e1[5] + e0[10] * e1[6] + e0[14] * e1[7];
            result[10] = e0[2] * e1[8] + e0[6] * e1[9] + e0[10] * e1[10] + e0[14] * e1[11];
            result[14] = e0[2] * e1[12] + e0[6] * e1[13] + e0[10] * e1[14] + e0[14] * e1[15];

            result[3] = e0[3] * e1[0] + e0[7] * e1[1] + e0[11] * e1[2] + e0[15] * e1[3];
            result[7] = e0[3] * e1[4] + e0[7] * e1[5] + e0[11] * e1[6] + e0[15] * e1[7];
            result[11] = e0[3] * e1[8] + e0[7] * e1[9] + e0[11] * e1[10] + e0[15] * e1[11];
            result[15] = e0[3] * e1[12] + e0[7] * e1[13] + e0[11] * e1[14] + e0[15] * e1[15];
        }

        return mResult;
    }

    transform4t& operator*=( FloatType x ) {
        for( int i = 0; i < 16; ++i )
            m_elements[i] *= x;
        return *this;
    }

    // Function: /
    // Reduces this matrix by a factor of b, that is,
    //  divides each value of this matrix by the specified float.
    // Note: This throws an exception if b is 0.0
    //
    // Parameters:
    // b - the factor to reduce the matrix by.
    //

    // Returns:
    // A matrix such that b * [Matrix result] = [this Matrix]
    transform4t operator/( FloatType b ) {
        if( b == FloatType( 0.0 ) ) {
            throw std::runtime_error( "can not divide matrix by zero" );
        }
        FloatType invB = FloatType( 1.0 ) / b;
        transform4t result = *this;
        for( int i = 0; i < 16; i++ ) {
            result.m_elements[i] *= invB;
        }
        return result;
    }

    transform4t operator+( const transform4t& b ) const {
        transform4t result;
        for( int i = 0; i < 16; ++i )
            result.m_elements[i] = m_elements[i] + b.m_elements[i];
        return result;
    }

    transform4t operator-( const transform4t& b ) const {
        transform4t result;
        for( int i = 0; i < 16; ++i )
            result.m_elements[i] = m_elements[i] - b.m_elements[i];
        return result;
    }

    transform4t& operator+=( const transform4t& b ) {
        for( int i = 0; i < 16; ++i )
            m_elements[i] += b.m_elements[i];
        return *this;
    }

    transform4t& operator-=( const transform4t& b ) {
        for( int i = 0; i < 16; ++i )
            m_elements[i] -= b.m_elements[i];
        return *this;
    }

#ifdef _MATRIX3_H
    // NOTE: This chops off the last row of data (or, the last element of each column).
    // NOTE: Matrix3 "row"s are the same as transform3f "column"s.  We are using
    //       column vectors, while 3ds Max is using row vectors.
    operator Matrix3() {
        // Point3's are always defined if Matrix3's are.
        Point3 row;
        Matrix3 mat;

        row.x = m_elements[0];
        row.y = m_elements[1];
        row.z = m_elements[2];
        mat.SetRow( 0, row );
        row.x = m_elements[4];
        row.y = m_elements[5];
        row.z = m_elements[6];
        mat.SetRow( 1, row );
        row.x = m_elements[8];
        row.y = m_elements[9];
        row.z = m_elements[10];
        mat.SetRow( 2, row );
        row.x = m_elements[12];
        row.y = m_elements[13];
        row.z = m_elements[14];
        mat.SetRow( 3, row );

        return mat;
    }
#endif

    // Elementary row operation: interchange two rows.
    void elementary_row_interchange( int row1, int row2 ) {
        for( int columnOffset = 0; columnOffset < 16; columnOffset += 4 )
            std::swap( m_elements[columnOffset + row1], m_elements[columnOffset + row2] );
    }

    // Elementary row operation: Multiply a row by a non-zero constant.
    void elementary_row_multiply( int row, FloatType multiplier ) {
        for( int columnOffset = 0; columnOffset < 16; columnOffset += 4 )
            m_elements[columnOffset + row] *= multiplier;
    }

    // Elementary row operation: Add a multiple of one row to another.
    void elementary_row_add_multiple( int rowSource, FloatType multiplier, int rowDest ) {
        for( int columnOffset = 0; columnOffset < 16; columnOffset += 4 )
            m_elements[columnOffset + rowDest] += multiplier * m_elements[columnOffset + rowSource];
    }

    // Will take pretend the transform is a 3x3 matrix of coefficients and will solve for conditions in vec.
    // Result is stored in vec, and the transform4t is no longer valid.
    // Much can be changed about this.
    bool solve_by_elimination_3x3( vector3f_type& vec ) {
        FloatType* nMat = &m_elements[0]; // Useful for traversing the elements with pointers instead of indices
        FloatType* rows[3];               // Used for row swapping

        // Augment this matrix with the passed values.
        nMat[12] = vec[0];
        nMat[13] = vec[1];
        nMat[14] = vec[2];

        // The virtual row array for simple, quick row swaps.
        // Points to first member of each row;
        rows[0] = nMat;
        rows[1] = nMat + 1;
        rows[2] = nMat + 2;

        int j = 0;
        for( int i = 0; i < 3; i++ ) {
            // Since this algorithm operates on rows and the consecutive data elements are column members,
            // we need to use this value to index in [row, column] order.
            int iColAddr = i << 2;

            // Find the largest possible pivot point
            FloatType max = abs( rows[i][iColAddr] );
            int imax = i;
            for( j = i + 1; j < 3; j++ ) {
                FloatType val = abs( rows[j][iColAddr] );
                if( val > max ) {
                    max = val;
                    imax = j;
                }
            }

            // Set the largest to the pivot point, or return false if the matrix is degenerate
            if( max == FloatType( 0.0 ) )
                return false;
            if( imax != i ) {
                FloatType* t = rows[i];
                rows[i] = rows[imax];
                rows[imax] = t;
            }

            // Inspect each pivot
            for( j = i + 1; j < 3; j++ ) {
                FloatType factor = rows[j][iColAddr] / rows[i][iColAddr];

                // Subtract to zero out the spots below the pivot
                // This traverses the row so we increment by 4, since adjacent row members are 4 elements apart
                for( int k = iColAddr + 4; k < 16; k += 4 )
                    rows[j][k] = rows[j][k] - factor * rows[i][k];
            }
        }

        // Back substitution
        vec.z = rows[2][12] / rows[2][8];
        vec.y = ( rows[1][12] - ( vec.z * rows[1][8] ) ) / rows[1][4];
        vec.x = ( rows[0][12] - ( vec.z * rows[0][8] ) - ( vec.y * rows[0][4] ) ) / rows[0][0];

        return true;
    }

    // Function: inverse
    // Calculates the inverse of this matrix.
    //
    // Returns:
    // The matrix inverse of this matrix.  That is, [Inverse Matrix]x[Matrix] = [Identity Matrix]
    transform4t to_inverse() const {
        using std::abs;

        // Defaults to the identity
        transform4t result;
        transform4t working = *this;
        // Do a row reduction, with the identity matrix tagging along for the ride.  Once
        // this matrix has been reduced to the identity, the result matrix has become the inverse.
        for( int column = 0; column < 4; ++column ) {
            // == Get the candidate value with the largest absolute value as the pivot ==
            // This value is chosen for numerical stability
            int pivotRow = column;
            int row = 0;
            FloatType absPivot = abs( working.get( column, pivotRow ) );
            for( row = pivotRow + 1; row < 4; ++row ) {
                FloatType testPivot = abs( working.get( column, row ) );
                if( testPivot > absPivot ) {
                    absPivot = testPivot;
                    pivotRow = row;
                }
            }
            if( absPivot == FloatType( 0.0 ) )
                throw std::runtime_error( "Cannot invert matrix " + str() + ", it is only of rank " +
                                          boost::lexical_cast<std::string>( column ) );

            // == Now use row operations to eliminate all the other values in this column ==
            // First, normalize the row to 1.
            FloatType pivotValue = working.get( column, pivotRow );
            result.elementary_row_multiply( pivotRow, FloatType( 1.0 ) / pivotValue );
            working.elementary_row_multiply( pivotRow, FloatType( 1.0 ) / pivotValue );
            // Then, subtract the correct multiple
            for( row = 0; row < 4; ++row ) {
                if( row != pivotRow ) {
                    FloatType multiplier = working.get( column, row );
                    if( multiplier != FloatType( 0 ) ) {
                        result.elementary_row_add_multiple( pivotRow, -multiplier, row );
                        working.elementary_row_add_multiple( pivotRow, -multiplier, row );
                    }
                }
            }
            // Finally, interchange the pivot row to its correct location if necessary.
            if( pivotRow != column ) {
                result.elementary_row_interchange( pivotRow, column );
                working.elementary_row_interchange( pivotRow, column );
            }
        }

        // Now the working matrix has been converted to the identity, and the result matrix
        // has been converted to this matrix's inverse.

        return result;
    }

    /**
     * Convert the upper 3x3 part of this matrix into an axis angle representation
     * (this is assumming that the transformation is a rotation matrix)
     *
     * \return a pair of the angle (in radians) and the axis which will produce this object's rotation component.
     */
    std::pair<FloatType, vector3f_type> to_angle_axis( FloatType epsilon = FloatType( 0.00001 ) ) const {
        using namespace frantic::math;

        const FloatType ZERO( 0.0 );
        const FloatType ONE( 1.0 );
        const FloatType TWO( 2.0 );
        const FloatType FOUR( 4.0 );

        FloatType angle;
        FloatType x;
        FloatType y;
        FloatType z;

        // It is important to check if we might induce a division by zero. This will happen when the rotation angle is
        // near 0 or pi. Unfortunately this special case makes up the bulk of the function.
        if( ( get_absolute( get( 0, 1 ) - get( 1, 0 ) ) < epsilon ) &&
            ( get_absolute( get( 0, 2 ) - get( 2, 0 ) ) < epsilon ) &&
            ( get_absolute( get( 2, 1 ) - get( 1, 2 ) ) < epsilon ) ) {

            // Check if this is the identity matrix (within some epsilon)
            if( is_identity( epsilon ) ) {
                return std::make_pair( ZERO, vector3f_type( ONE, ZERO, ZERO ) );
            }

            // otherwise we assume this arose from angle = pi, determine the axis
            angle = FloatType( M_PI );
            FloatType xx = ( get( 0, 0 ) + ONE ) / TWO;
            FloatType yy = ( get( 1, 1 ) + ONE ) / TWO;
            FloatType zz = ( get( 2, 2 ) + ONE ) / TWO;
            FloatType xy = ( get( 0, 1 ) + get( 1, 0 ) ) / FOUR;
            FloatType xz = ( get( 0, 2 ) + get( 2, 0 ) ) / FOUR;
            FloatType yz = ( get( 1, 2 ) + get( 2, 1 ) ) / FOUR;

            // Find the largest term in in the diagonal, and use it to determine the axis
            if( ( xx > yy ) && ( xx > zz ) ) { // m[0][0] is the largest diagonal term
                if( xx < epsilon ) {
                    x = ZERO;
                    y = FloatType( 0.7071067811865476 ); // 1.0 / sqrt(2.0);
                    z = FloatType( 0.7071067811865476 ); // 1.0 / sqrt(2.0);
                } else {
                    x = sqrt( xx );
                    y = xy / x;
                    z = xz / x;
                }
            } else if( yy > zz ) { // m[1][1] is the largest diagonal term
                if( yy < epsilon ) {
                    x = FloatType( 0.7071067811865476 ); // 1.0 / std::sqrt(2.0);
                    y = ZERO;
                    z = FloatType( 0.7071067811865476 ); // 1.0 / std::sqrt(2.0);
                } else {
                    y = sqrt( yy );
                    x = xy / y;
                    z = yz / y;
                }
            } else { // m[2][2] is the largest diagonal term so base result on this
                if( zz < epsilon ) {
                    x = FloatType( 0.7071067811865476 ); // 1.0 / sqrt(2.0);
                    y = FloatType( 0.7071067811865476 ); // 1.0 / sqrt(2.0);
                    z = ZERO;
                } else {
                    z = sqrt( zz );
                    x = xz / z;
                    y = yz / z;
                }
            }

            return std::make_pair( angle, vector3f_type( x, y, z ) );
        }

        // general case, solve for cos(angle), angle, and the axis
        FloatType s = sqrt( square( get( 2, 1 ) - get( 1, 2 ) ) + square( get( 0, 2 ) - get( 2, 0 ) ) +
                            square( get( 1, 0 ) - get( 0, 1 ) ) );

        // final measure to prevent division by zero, should not happen if the
        // original matrix was valid
        if( get_absolute( s ) < epsilon ) {
            s = FloatType( 1.0 );
        }

        angle = acos( ( get( 0, 0 ) + get( 1, 1 ) + get( 2, 2 ) - ONE ) / TWO );

        x = ( get( 2, 1 ) - get( 1, 2 ) ) / s;
        y = ( get( 0, 2 ) - get( 2, 0 ) ) / s;
        z = ( get( 1, 0 ) - get( 0, 1 ) ) / s;

        return std::make_pair( angle, vector3f_type( x, y, z ) );
    }

    vector3f_type projection_transform( vector3f_type v ) const {
        FloatType x = v.x, y = v.y, z = v.z;

        FloatType x2 = x * m_elements[0] + y * m_elements[4] + z * m_elements[8] + m_elements[12];
        FloatType y2 = x * m_elements[1] + y * m_elements[5] + z * m_elements[9] + m_elements[13];
        FloatType z2 = x * m_elements[2] + y * m_elements[6] + z * m_elements[10] + m_elements[14];
        FloatType w2 = x * m_elements[3] + y * m_elements[7] + z * m_elements[11] + m_elements[15];

        if( w2 != FloatType( 0 ) )
            return vector3f_type( x2, y2, z2 ) / w2;

        return vector3f_type();
    }

    static transform4t parse( const std::string& input ) {
        std::string s = input;

        s = frantic::strings::string_replace( s, "matrix3", "" );
        s = frantic::strings::string_replace( s, "matrix", "" );

        transform4t result;

        std::vector<std::string> tokens;
        frantic::strings::split( s, tokens, "()[],\n\r " );

        try {
            if( tokens.size() == 12 ) {
                int i = 0;
                for( int col = 0; col < 4; col++ ) {
                    for( int row = 0; row < 3; row++ ) {
                        result.set( col, row, boost::lexical_cast<FloatType>( tokens[i++] ) );
                    }
                }
            } else if( tokens.size() == 16 ) {
                int i = 0;
                for( int col = 0; col < 4; col++ ) {
                    for( int row = 0; row < 4; row++ ) {
                        result.set( col, row, boost::lexical_cast<FloatType>( tokens[i++] ) );
                    }
                }
            } else {
                throw std::runtime_error( "Could not parse matrix due to incorrect number of tokens:" + input );
            }
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "could not parse matrix due to invalid tokens: " + input );
        }

        return result;
    }

    vector3f_type transform_no_translation( const vector3f_type& v ) const {
        FloatType x = v.x, y = v.y, z = v.z;
        return vector3f_type( x * m_elements[0] + y * m_elements[4] + z * m_elements[8],
                              x * m_elements[1] + y * m_elements[5] + z * m_elements[9],
                              x * m_elements[2] + y * m_elements[6] + z * m_elements[10] );
    }

    vector3f_type transpose_transform_no_translation( const vector3f_type& v ) const {
        FloatType x = v.x, y = v.y, z = v.z;
        return vector3f_type( x * m_elements[0] + y * m_elements[1] + z * m_elements[2],
                              x * m_elements[4] + y * m_elements[5] + z * m_elements[6],
                              x * m_elements[8] + y * m_elements[9] + z * m_elements[10] );
    }

    std::string to_readable_string() {
        std::stringstream ss;
        ss << m_elements[0] << " " << m_elements[4] << " " << m_elements[8] << " " << m_elements[12] << std::endl;
        ss << m_elements[1] << " " << m_elements[5] << " " << m_elements[9] << " " << m_elements[13] << std::endl;
        ss << m_elements[2] << " " << m_elements[6] << " " << m_elements[10] << " " << m_elements[14] << std::endl;
        ss << m_elements[3] << " " << m_elements[7] << " " << m_elements[11] << " " << m_elements[15] << std::endl;
        return ss.str();
    }

    std::string str() const;
};

// todo: maybe make a friend so we can access the array data directly? batty, 2005/05/05
template <class CharType, typename FloatType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out,
                                                 const transform4t<FloatType>& xfrm ) {
    out << '[' << xfrm[0] << ", " << xfrm[1] << ", " << xfrm[2] << ", " << xfrm[3] << "],\n";
    out << '[' << xfrm[4] << ", " << xfrm[5] << ", " << xfrm[6] << ", " << xfrm[7] << "],\n";
    out << '[' << xfrm[8] << ", " << xfrm[9] << ", " << xfrm[10] << ", " << xfrm[11] << "],\n";
    out << '[' << xfrm[12] << ", " << xfrm[13] << ", " << xfrm[14] << ", " << xfrm[15] << ']';
    return out;
}

template <typename FloatType>
inline std::string transform4t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

template <typename FloatType>
inline transform4t<FloatType> operator*( FloatType b, const transform4t<FloatType>& a ) {
    transform4t<FloatType> result = a;
    for( int i = 0; i < 16; i++ )
        result[i] *= b;
    return result;
}

template <typename FloatType>
inline transform4t<FloatType> operator*( const transform4t<FloatType>& a, FloatType b ) {
    return b * a;
}

// Only matrix multiplication with a vector on the right is defined, because we are using
// the convention of column vectors.
template <typename FloatType>
inline typename transform4t<FloatType>::vector3f_type
operator*( const transform4t<FloatType>& m, const typename transform4t<FloatType>::vector3f_type& v ) {
    typedef typename transform4t<FloatType>::vector3f_type vector3f_type;
    FloatType x = v.x, y = v.y, z = v.z;
    if( m.is_perspective() ) {
        FloatType w = x * m[3] + y * m[7] + z * m[11] + m[15];
        // We need the case without dividing through by w==0 when we're multiplying
        // a position (which has w=1) by a derivative matrix (which has all zeros in
        // the bottom row), to get a vector (which has w=0, but we don't store it so that
        // information disappears).  For example, in the transformed_particle_istream class
        // when modifying the velocity channel.
        if( w != FloatType( 0 ) )
            return vector3f_type( ( x * m[0] + y * m[4] + z * m[8] + m[12] ) / w,
                                  ( x * m[1] + y * m[5] + z * m[9] + m[13] ) / w,
                                  ( x * m[2] + y * m[6] + z * m[10] + m[14] ) / w );
        else
            return vector3f_type( ( x * m[0] + y * m[4] + z * m[8] + m[12] ),
                                  ( x * m[1] + y * m[5] + z * m[9] + m[13] ),
                                  ( x * m[2] + y * m[6] + z * m[10] + m[14] ) );
    } else {
        return vector3f_type( x * m[0] + y * m[4] + z * m[8] + m[12], x * m[1] + y * m[5] + z * m[9] + m[13],
                              x * m[2] + y * m[6] + z * m[10] + m[14] );
    }
}

// Whether this matrix is orthogonal (A * transpose(A) == I)
template <typename FloatType>
inline bool transform4t<FloatType>::is_orthogonal() const {
    transform4t<FloatType> testMatrix = ( *this ) * transform4t<FloatType>::from_transpose( *this );
    //	std::cerr << "is_orthogonal test matrix: " << testMatrix << std::endl;
    return testMatrix.is_identity();
}

// Whether this matrix is orthogonal (A * transpose(A) == I) within a given tolerance
template <typename FloatType>
inline bool transform4t<FloatType>::is_orthogonal( FloatType tolerance ) const {
    transform4t<FloatType> testMatrix = ( *this ) * transform4t<FloatType>::from_transpose( *this );
    //	std::cerr << "is_orthogonal test matrix: " << testMatrix << std::endl;
    return testMatrix.is_identity( tolerance );
}

// Whether this matrix is symmetric (A == transpose(A))
template <typename FloatType>
inline bool transform4t<FloatType>::is_symmetric() const {
    return m_elements[1] == m_elements[4] && m_elements[2] == m_elements[8] && m_elements[3] == m_elements[12] &&
           m_elements[6] == m_elements[9] && m_elements[7] == m_elements[13] && m_elements[11] == m_elements[14];
}

// Whether this matrix is symmetric (A == transpose(A)) within a given tolerance
template <typename FloatType>
inline bool transform4t<FloatType>::is_symmetric( FloatType tolerance ) const {
    FloatType errorSquared = frantic::math::square( m_elements[1] - m_elements[4] ) +
                             frantic::math::square( m_elements[2] - m_elements[8] ) +
                             frantic::math::square( m_elements[3] - m_elements[12] ) +
                             frantic::math::square( m_elements[6] - m_elements[9] ) +
                             frantic::math::square( m_elements[7] - m_elements[13] ) +
                             frantic::math::square( m_elements[11] - m_elements[14] );
    return errorSquared <= tolerance * tolerance;
}

template <typename FloatType>
void transform4t<FloatType>::decompose( transform4t<FloatType>& outPerspective, vector3f_type& outTranslation,
                                        transform4t<FloatType>& outRotation,
                                        transform4t<FloatType>& outStretch ) const {
    transform4t<FloatType> partialXform = *this;

    // Extract the perspective transform component first
    outPerspective.set_to_identity();
    if( is_perspective() ) {
        // Remove the perspective portion from the partial transform and save it in
        // the perspective matrix.
        partialXform.m_elements[3] = 0;
        partialXform.m_elements[7] = 0;
        partialXform.m_elements[11] = 0;
        partialXform.m_elements[15] = 1;
        // This is a very inefficient way to solve for the perspective portion, but in
        // comparison to the polar decomposition iteration it is insignificant.
        transform4t<FloatType> perspExtract = partialXform.to_inverse();
        perspExtract.set_translation( vector3f( 0 ) );
        perspExtract.transpose();
        vector3f_type perspectiveCoords = perspExtract * vector3f_type( m_elements[3], m_elements[7], m_elements[11] );
        outPerspective.m_elements[3] = perspectiveCoords.x;
        outPerspective.m_elements[7] = perspectiveCoords.y;
        outPerspective.m_elements[11] = perspectiveCoords.z;
        outPerspective.m_elements[15] = m_elements[15] - m_elements[12] * perspectiveCoords.x -
                                        m_elements[13] * perspectiveCoords.y - m_elements[14] * perspectiveCoords.z;
    }

    // Extract the translation, and set it to zero in the partial transform.
    outTranslation.x = m_elements[12];
    outTranslation.y = m_elements[13];
    outTranslation.z = m_elements[14];
    partialXform.m_elements[12] = 0;
    partialXform.m_elements[13] = 0;
    partialXform.m_elements[14] = 0;

    partialXform.decompose_polar( outRotation, outStretch, 0.000001f );
}

template <typename FloatType>
inline ray3t<FloatType> operator*( const transform4t<FloatType>& m, const ray3t<FloatType>& r ) {
    return ray3t<FloatType>( m * r.origin(), m.transform_no_translation( r.direction() ) );
}

template <typename FloatType>
inline boundbox3t<FloatType> operator*( const transform4t<FloatType>& t, const boundbox3t<FloatType>& b ) {
    boundbox3t<FloatType> result( vector3t<FloatType>( t * b.get_corner( 0 ) ) );
    for( int i = 1; i < 8; ++i ) {
        result += vector3t<FloatType>( t * b.get_corner( i ) );
    }
    return result;
}

/**
 * For matrix multiplication, when the second argument's float type is bigger, we want the result
 * of the multiplcation to have that bigger size. This overload enables that.
 */
template <typename FloatType1, typename FloatType2>
inline typename boost::enable_if_c<( sizeof( FloatType1 ) < sizeof( FloatType2 ) ), transform4t<FloatType2>>::type
operator*( const transform4t<FloatType1>& t1, const transform4t<FloatType2>& t2 ) {
    return transform4t<FloatType2>( t1 ) * t2;
}

typedef transform4t<float> transform4f;
typedef transform4t<double> transform4fd;

} // namespace graphics
} // namespace frantic
