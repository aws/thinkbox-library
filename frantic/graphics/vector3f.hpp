// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <cstring>
#include <ostream>

#include <frantic/math/utils.hpp>

#include <frantic/graphics2d/size2.hpp>
#include <frantic/graphics2d/size2f.hpp>
#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

namespace frantic {
namespace graphics {

template <typename FloatType>
class vector3t {

  public:
    FloatType x, y, z;
    typedef FloatType float_type;

    vector3t( const float_type& X, const float_type& Y, const float_type& Z )
        : x( X )
        , y( Y )
        , z( Z ) {}

    explicit vector3t( const float_type& X )
        : x( X )
        , y( X )
        , z( X ) {}

    explicit vector3t( int X )
        : x( (float_type)X )
        , y( (float_type)X )
        , z( (float_type)X ) {}

    explicit vector3t( const float_type vec[3] )
        : x( vec[0] )
        , y( vec[1] )
        , z( vec[2] ) {}

    vector3t()
        : x( 0.f )
        , y( 0.f )
        , z( 0.f ) {
        // There was an issue with y being set to a 0 integer, it was coming out as a REALLY small float, but not 0.
        // I changed these to floats and it worked.
    }

    /** Implicit conversion from smaller float type is allowed. */
    template <class OtherFloatType>
    vector3t( const vector3t<OtherFloatType>& t,
              typename boost::enable_if_c<( sizeof( OtherFloatType ) < sizeof( FloatType ) ), int*>::type = 0 )
        : x( t.x )
        , y( t.y )
        , z( t.z ) {}

    /** Explicit cast required for conversion from larger float type. */
    template <class OtherFloatType>
    explicit vector3t( const vector3t<OtherFloatType>& t,
                       typename boost::enable_if_c<( sizeof( OtherFloatType ) > sizeof( FloatType ) ), int*>::type = 0 )
        : x( (float_type)t.x )
        , y( (float_type)t.y )
        , z( (float_type)t.z ) {}

    static vector3t from_floor( const vector3t& v ) { return vector3t( floor( v.x ), floor( v.y ), floor( v.z ) ); }

    // Returns a uniform random point from within the specified tetrahedron
    vector3t from_random_in_tetrahedron( const vector3t& v0, const vector3t& v1, const vector3t& v2,
                                         const vector3t& v3 ) {
        float_type r = (float_type)rand() / RAND_MAX;
        float_type s = (float_type)rand() / RAND_MAX;
        float_type t = (float_type)rand() / RAND_MAX;

        // fold the cube into a prism
        if( r + s > 1.f ) {
            r = 1.f - s;
            s = 1.f - t;
        }

        // fold the prism into a tetrahedron
        if( s + t > 1.f ) {
            float_type tmp = t;
            t = 1.f - r - s;
            s = 1.f - tmp;
        } else if( r + s + t > 1.f ) {
            float_type tmp = t;
            t = r + s + t - 1.f;
            r = 1 - s - tmp;
        }
        // This results in barycentric coordinates for creating the random point in the tet
        return ( 1 - r - s - t ) * v0 + r * v1 + s * v2 + t * v3;
    }

    // Generates a random vector in the unit cube [0,0,0] to [1,1,1]
    static vector3t from_random() {
        return vector3t( (float_type)rand() / RAND_MAX, (float_type)rand() / RAND_MAX, (float_type)rand() / RAND_MAX );
    }

    // Generates a random vector in the unit cube [0,0,0] to [1,1,1] (as long as that's what rng does)
    template <class RandomNumberGenerator>
    static vector3t from_random( RandomNumberGenerator& rng ) {
        return vector3t( rng(), rng(), rng() );
    }

    // Generates a random vector with each component having a gaussian distribution with variance 1
    static vector3t from_random_gaussian() {
        using std::cos;
        using std::log;
        using std::sin;
        using std::sqrt;

        bool goodResult = false;
        vector3t result;
        while( !goodResult ) {
            // Use the non-polar form of the box-muller transformation
            float_type x = (float_type)rand() / RAND_MAX, y = (float_type)rand() / RAND_MAX,
                       z = (float_type)rand() / RAND_MAX, w = (float_type)rand() / RAND_MAX;
            double coefficient = sqrt( -2 * log( x ) );
            result.set( float_type( coefficient * cos( 2 * M_PI * y ) ),
                        float_type( coefficient * sin( 2 * M_PI * y ) ),
                        float_type( sqrt( -2 * log( z ) ) * cos( 2 * M_PI * w ) ) );
            // If all the numbers aren't NaN, it's good
            goodResult = frantic::math::is_finite( result.x ) && frantic::math::is_finite( result.y ) &&
                         frantic::math::is_finite( result.z );
        }
        return result;
    }

    // Generates a random vector with each component having a gaussian distribution with variance 1
    // NOTE: rng should generate random values uniformly distributed in [0,1]
    template <class RandomNumberGenerator>
    static vector3t from_random_gaussian( RandomNumberGenerator& rng ) {
        using std::cos;
        using std::log;
        using std::sin;
        using std::sqrt;

        bool goodResult = false;
        vector3t result;
        while( !goodResult ) {
            // Use the non-polar form of the box-muller transformation
            float_type x = rng(), y = rng(), z = rng(), w = rng();
            double coefficient = sqrt( -2 * log( x ) );
            result.set( float_type( coefficient * cos( 2 * M_PI * y ) ),
                        float_type( coefficient * sin( 2 * M_PI * y ) ),
                        float_type( sqrt( -2 * log( z ) ) * cos( 2 * M_PI * w ) ) );
            // If all the numbers aren't NaN, it's good
            goodResult = result.is_finite();
        }
        return result;
    }

    // TODO: We should use the boost random number generator for high quality and fast random numbers
    static vector3t from_unit_random() {
        vector3t result = from_random_gaussian();
        // Normalize to unit distance
        result.normalize();
        return result;
    }

    template <typename RandomNumberGenerator>
    static vector3t from_unit_random( RandomNumberGenerator& rng ) {
        vector3t result = from_random_gaussian( rng );
        // Normalize to unit distance
        result.normalize();
        return result;
    }

    /**
     *  Return a random vector within a sphere of the given radius.
     * The specified random number generator must produce random numbers
     * uniformly distributed in [0,1].  With such a random number generator,
     * the output vectors will be uniformly distributed within the
     * sphere.
     *
     * @param rng a random number generator.  This must generate uniformly
     *		distributed numbers in the range [0,1].
     * @param radius the radius of a sphere in which to produce random
     *		vectors.
     * @return vector3f a vector within the specified sphere.
     */
    template <typename RandomNumberGenerator>
    static vector3t from_random_in_sphere( RandomNumberGenerator& rng, const float_type radius = 1.f ) {
        const vector3t randomDirection( from_random_gaussian( rng ) );
        float_type randomRadius( radius * pow( rng(), 1.f / 3 ) );
        float_type magnitudeSquared = randomDirection.get_magnitude_squared();

        if( magnitudeSquared > 0 ) {
            float_type f = randomRadius / sqrt( magnitudeSquared );
            return f * randomDirection;
        } else {
            return vector3t();
        }
    }

    /**
     * @overload
     */
    static vector3t from_random_in_sphere( const float_type radius = 1.f ) {
        const vector3t randomDirection( from_random_gaussian() );
        float_type randomRadius( radius * pow( rand() / (float_type)RAND_MAX, 1.f / 3 ) );
        float_type magnitudeSquared = randomDirection.get_magnitude_squared();

        if( magnitudeSquared > 0 ) {
            float_type f = randomRadius / sqrt( magnitudeSquared );
            return f * randomDirection;
        } else {
            return vector3t();
        }
    }

    /**
     * Returns a random point within a sphere of the specified radius. This method uses a rejection approach where
     * points are generated in the enclosing box and rejected if they fall outside the sphere. Given the ration of
     * volumes, we expect this to iterate ~ 2 times on average before selecting a point. The ratio of volumes is (6 /
     * pi) ~ 1.909. With a suitably performant generator, this should be faster than the direct approach in
     * from_random_in_sphere() which makes use of sqrt, cos, and others. \tparam Generator Must model the 'Generator'
     * concept from Boost.Random. See
     * http://www.boost.org/doc/libs/1_52_0/doc/html/boost_random/reference.html#boost_random.reference.generators
     * \param gen A reference to a Generator object for producing random bits.
     * \param radius The radius of the sphere to generate in.
     * \return A uninform random point in the sphere.
     */
    template <class Generator>
    static vector3t from_random_in_sphere_rejection( Generator& gen, const float_type radius = 1.f ) {
        boost::variate_generator<Generator&, boost::uniform_real<float_type>> rnd(
            gen, boost::uniform_real<float_type>( -1.f, 1.f ) );

        vector3t result;
        float_type magnitude;

        // Generate a random point in a [-1,-1,-1] to [1,1,1] box, discarding those outside the sphere. I expect this to
        // take on average 2 tries before a point is accepted.
        do {
            result.x = rnd();
            result.y = rnd();
            result.z = rnd();

            magnitude = result.get_magnitude_squared();
        } while( magnitude >= 1.f );

        return radius * result;
    }

    // TODO: Move all projection functions to camera class
    // ****************************************************************************
    // CONVERSION FROM 2D COORDINATE SYSTEMS
    //
    // All of the to_<coordinate system> functions share the following properties:
    // -The input coordinates are not normalized
    // -The origin (0,0) is the bottom left of the image
    //
    //   y
    //   ^
    //   |
    //   |
    //   |
    // (0,0)-----> x
    //
    // ****************************************************************************

    //
    // v is an coordinate (not normalized)
    // s is the size of the image
    // horizontalFovRadians is the angle from the left side of the image to the right side
    //
    // This perspective projection is projected onto x,y
    // The reason this is operating towards negative z is that max cameras operate
    // that way.  See also the version of this function in spherical_coords.
    static vector3t from_perspective_projection( const frantic::graphics2d::vector2f& v,
                                                 const frantic::graphics2d::size2f& size,
                                                 const float_type& tanHalfHorizontalFov,
                                                 const float_type& pixelAspect = 1.0f ) {
        const double tanHalfVerticalFov = tanHalfHorizontalFov * pixelAspect * size.ysize / size.xsize;

        const double xfactor = 2.f * tanHalfHorizontalFov;
        const double yfactor = 2.f * tanHalfVerticalFov;

        const float_type nx = float_type( xfactor * ( double( v.x ) / size.xsize - 0.5 ) );
        // const float_type ny = float_type(yfactor * (double(v.y) / size.xsize - 0.5 * size.ysize / size.xsize));
        const float_type ny = float_type( yfactor * ( double( v.y ) / size.ysize - 0.5 ) );

        // std::cout << "<dbg: " << v << "  nx: " << nx << "  ny: " << ny << std::endl;
        return vector3t( nx, ny, -1.f );
    }

    // v takes regular coordinates (not normalized)
    // TODO: pass around tan(half-fov) for efficiency
    static vector3t from_perspective_projection( const frantic::graphics2d::vector2& v,
                                                 const frantic::graphics2d::size2& size,
                                                 const float_type& tanHalfHorizontalFov,
                                                 const float_type& pixelAspect = 1.0f ) {
        return from_perspective_projection( v.to_image_coord(),
                                            frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ),
                                            tanHalfHorizontalFov, pixelAspect );
    }

    // TODO: move to camera class
    // Computes the 3dimensional view vector which corresponds to the latitude / longitude values
    // Longitude corresponds to the angle on the ZX plane starting from Z and going in the -X dir
    // Latitude corresponds to the angle from positive Y to negative Y
    //
    // v is the coordinate on the 2d longlat image where v.x is longitude and v.y is latitude
    // s is the size of the 2d longlat image
    static vector3t from_longlat_yup( const frantic::graphics2d::vector2f& v,
                                      const frantic::graphics2d::size2f& size ) {
        double theta = 2 * M_PI * ( v.x / size.xsize );
        double phi = M_PI * ( ( size.ysize - v.y ) / size.ysize );

        return vector3t( -float_type( sin( phi ) * sin( theta ) ), float_type( cos( phi ) ),
                         float_type( sin( phi ) * cos( theta ) ) );
    }

    // TODO: move to camera class
    static vector3t from_longlat_yup( frantic::graphics2d::vector2 v, frantic::graphics2d::size2 size ) {
        return from_longlat_yup( frantic::graphics2d::vector2f( v.x + 0.5f, v.y + 0.5f ),
                                 frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ) );
    }

    // TODO: move to camera class
    // v is the coordinate (not normalized) in a 2d fisheye image
    // (0, 0) is the bottom left corner
    // size is the size of the image
    // fov is in radians and corresponds to the X-axis
    // Z-up projection
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    static vector3t from_fisheye( const frantic::graphics2d::vector2f& v, const frantic::graphics2d::size2f& size,
                                  const float_type& horizontalFovRadians, bool& outIsValidInput,
                                  const float_type& pixelAspect = 1.0 ) {
        float_type x = v.x / size.xsize * 2 - 1;
        // float_type y = (v.y / size.ysize * 2 - 1) * (size.ysize / size.xsize);
        float_type y = ( v.y / size.ysize * 2 - 1 ) * ( pixelAspect * size.ysize / size.xsize );

        // Get angles
        double r = sqrt( x * x + y * y );
        double theta = atan2( y, x );
        double phi = r * ( horizontalFovRadians / 2 );

        // The input direction is valid only if the angle is less than 180 degrees, i.e. A
        // full 360 degree fisheye.
        outIsValidInput = outIsValidInput && ( phi <= M_PI );

        return vector3t( float_type( sin( phi ) * cos( theta ) ), float_type( sin( phi ) * sin( theta ) ),
                         -float_type( cos( phi ) ) );
    }

    // TODO: move to camera class
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    static vector3t from_fisheye( const frantic::graphics2d::vector2& v, const frantic::graphics2d::size2& size,
                                  const float_type& horizontalFovRadians, bool& outIsValidInput,
                                  const float_type& pixelAspect = 1.0 ) {
        return from_fisheye( frantic::graphics2d::vector2f( v.x + 0.5f, v.y + 0.5f ),
                             frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ),
                             horizontalFovRadians, outIsValidInput, pixelAspect );
    }

    // TODO: move to camera class
    static vector3t from_cylindrical( const frantic::graphics2d::vector2f& v, const frantic::graphics2d::size2f& size,
                                      const float_type& verticalFovRadians ) {
        double theta = 2 * M_PI * ( v.x / size.xsize );
        double phi = verticalFovRadians * ( ( size.ysize - v.y ) / size.ysize ) + ( M_PI - verticalFovRadians ) / 2;

        return vector3t( -float_type( sin( phi ) * sin( theta ) ), float_type( cos( phi ) ),
                         float_type( sin( phi ) * cos( theta ) ) );
    }

    // TODO: move to camera class
    static vector3t from_cylindrical( const frantic::graphics2d::vector2& v, const frantic::graphics2d::size2& size,
                                      const float_type& verticalFovRadians ) {
        return from_cylindrical( frantic::graphics2d::vector2f( v.x + 0.5f, v.y + 0.5f ),
                                 frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ),
                                 verticalFovRadians );
    }

    /**
     * This is useful for chaining operations, ie. vector3f v; v.copy().replace_element(3,1.f);
     * The copy-constructor does not return an l-value, so that's why I made this function.
     */
    inline vector3t copy() const { return vector3t( *this ); }

    inline vector3t& replace_element( int i, float_type f ) {
        ( &x )[i] = f;
        return *this;
    }

    // TODO: Move to camera class
    // ****************************************************************************
    // CONVERSION TO 2D COORDINATE SYSTEMS
    //
    // All of the to_<coordinate system> functions share the following properties:
    // -The returned coordinates are dependant on the size passed to the function
    // -The origin (0,0) is the bottom left of the image
    //
    //   y
    //   ^
    //   |
    //   |
    //   |
    // (0,0)-----> x
    //
    // ****************************************************************************

    // TODO: Move to camera class
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2f to_perspective_projection( frantic::graphics2d::size2f size,
                                                             float tanHalfHorizontalFov, bool& outIsValidInput,
                                                             float pixelAspect = 1.0f ) const {
        outIsValidInput = outIsValidInput && ( z < 0 );
        // May 7, 2009
        // I switched this because just slightly negative z values were projecting really far
        // away from just slightly positive z values. Neither solution is great, but it sort of
        // looks better this way.
        // double planeWidth = -z * tanHalfHorizontalFov;
        double planeWidth = std::abs( z ) * tanHalfHorizontalFov;
        double planeHeight = planeWidth * pixelAspect * size.ysize / size.xsize;
        double xVal = 0.5 * size.xsize * ( double( x ) / planeWidth + 1.0 );
        double yVal = 0.5 * size.ysize * ( double( y ) / planeHeight + 1.0 );
        return frantic::graphics2d::vector2f( float( xVal ), float( yVal ) );
    }

    // TODO: Move to camera class
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2 to_perspective_projection( frantic::graphics2d::size2 size, float tanHalfHorizontalFov,
                                                            bool& outIsValidInput, float pixelAspect = 1.0f ) const {
        frantic::graphics2d::vector2f v = to_perspective_projection(
            frantic::graphics2d::size2f( float_type( size.xsize ), float_type( size.ysize ) ), tanHalfHorizontalFov,
            outIsValidInput, pixelAspect );
        return frantic::graphics2d::vector2( (int)floorf( v.x ), (int)floorf( v.y ) );
    }

    // TODO: Move to camera class
    frantic::graphics2d::vector2f to_fisheye( const frantic::graphics2d::size2f& size, double horizontalFovRadians,
                                              float pixelAspect = 1.0f ) const {
        // This vector may not be normalized, so we have to divide by its norm to get the cosine.
        double phi = acos( -z / sqrt( x * x + y * y + z * z ) );
        double theta = atan2( y, x );

        double normX = cos( theta ) * phi * 2 / horizontalFovRadians;
        double normY = sin( theta ) * phi * 2 / ( horizontalFovRadians * pixelAspect * size.ysize ) * size.xsize;

        return frantic::graphics2d::vector2f( static_cast<float>( size.xsize * ( normX + 1 ) / 2.f ),
                                              static_cast<float>( size.ysize * ( normY + 1 ) / 2.f ) );
    }

    // TODO: Move to camera class
    frantic::graphics2d::vector2 to_fisheye( const frantic::graphics2d::size2& size, double horizontalFovRadians,
                                             float pixelAspect = 1.0f ) const {
        frantic::graphics2d::vector2f v =
            to_fisheye( frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ), horizontalFovRadians,
                        pixelAspect );
        return frantic::graphics2d::vector2( (int)floorf( v.x ), (int)floorf( v.y ) );
    }

    // TODO: Move to camera class
    frantic::graphics2d::vector2f to_longlat_yup( const frantic::graphics2d::size2f& size ) const {
        vector3t normalized = *this;
        normalized.normalize();

        double phi = M_PI - acos( normalized.y );
        double theta = atan2( -normalized.x, normalized.z );
        const double TWO_PI = 2 * M_PI;

        if( theta < 0 )
            theta += TWO_PI;

        return frantic::graphics2d::vector2f( float( theta / TWO_PI ) * ( size.xsize - 1 ) + .5f,
                                              float( phi / M_PI ) * ( size.ysize - 1 ) + .5f );
    }

    // TODO: Move to camera class
    frantic::graphics2d::vector2 to_longlat_yup( const frantic::graphics2d::size2& size ) const {
        frantic::graphics2d::vector2f v =
            to_longlat_yup( frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ) );
        return frantic::graphics2d::vector2( (int)floorf( v.x ), (int)floorf( v.y ) );
    }

    // TODO: Move to camera class
    frantic::graphics2d::vector2f to_cylindrical( const frantic::graphics2d::size2f& size,
                                                  double verticalFovRadians ) const {
        vector3t normalized = *this;
        normalized.normalize();

        double phi = ( M_PI - acos( normalized.y ) ) - ( M_PI - verticalFovRadians ) / 2;
        double theta = atan2( -normalized.x, normalized.z );
        const double TWO_PI = 2 * M_PI;

        if( theta < 0 )
            theta += TWO_PI;

        return frantic::graphics2d::vector2f( float( theta / TWO_PI ) * size.xsize,
                                              float( phi / verticalFovRadians ) * size.ysize );
    }

    // TODO: Move to camera class
    frantic::graphics2d::vector2 to_cylindrical( const frantic::graphics2d::size2& size,
                                                 double verticalFovRadians ) const {
        frantic::graphics2d::vector2f v = to_cylindrical(
            frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ), verticalFovRadians );
        return frantic::graphics2d::vector2( (int)floorf( v.x ), (int)floorf( v.y ) );
    }

    frantic::graphics2d::vector2f project_xy() const { return frantic::graphics2d::vector2f( float( x ), float( y ) ); }

    static vector3t from_axis( int axis ) {
        switch( axis ) {
        case 0:
            return vector3t( 1, 0, 0 );
        case 1:
            return vector3t( 0, 1, 0 );
        case 2:
            return vector3t( 0, 0, 1 );
        default:
            return vector3t();
        }
    }

    static vector3t from_xaxis() { return vector3t( 1, 0, 0 ); }

    static vector3t from_yaxis() { return vector3t( 0, 1, 0 ); }

    static vector3t from_zaxis() { return vector3t( 0, 0, 1 ); }

    static float_type distance_squared( const vector3t& a, const vector3t& b ) {
        return ( a.x - b.x ) * ( a.x - b.x ) + ( a.y - b.y ) * ( a.y - b.y ) + ( a.z - b.z ) * ( a.z - b.z );
    }

    static float_type distance( const vector3t& a, const vector3t& b ) {
        return sqrt( vector3t::distance_squared( a, b ) );
    }

    static double distance_squared_double( const vector3t& a, const vector3t& b ) {
        double x = double( a.x ) - b.x;
        double y = double( a.y ) - b.y;
        double z = double( a.z ) - b.z;
        return x * x + y * y + z * z;
    }

    static double distance_double( const vector3t& a, const vector3t& b ) {
        return sqrt( distance_squared_double( a, b ) );
    }

    float_type get_magnitude_squared() const { return x * x + y * y + z * z; }

    float_type get_magnitude() const {
        using std::sqrt;
        return sqrt( get_magnitude_squared() );
    }

    static vector3t max_magnitude( const vector3t& a, const vector3t& b ) {
        return a.get_magnitude_squared() >= b.get_magnitude_squared() ? a : b;
    }

    float_type max_abs_component() const {
        return ( std::max )( std::abs( x ), ( std::max )( std::abs( y ), std::abs( z ) ) );
    }

    int get_largest_axis() const {
        float absX = std::abs( x ), absY = std::abs( y ), absZ = std::abs( z );
        if( absX >= absY && absX >= absZ )
            return 0;
        else if( absY >= absZ )
            return 1;
        else
            return 2;
    }

    static float_type dot( const vector3t& a, const vector3t& b ) { return a.x * b.x + a.y * b.y + a.z * b.z; }

    static double dot_double( const vector3t& a, const vector3t& b ) {
        return (double)a.x * (double)b.x + (double)a.y * (double)b.y + (double)a.z * (double)b.z;
    }

    static vector3t cross( const vector3t& a, const vector3t& b ) {
        return vector3t( a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x );
    }

    /**
     * Gets the angle in radians between the two vectors
     */
    static float_type angle( const vector3t& a, const vector3t& b ) {
        using std::acos;
        using std::sqrt;

        return acos( dot( a, b ) / sqrt( a.get_magnitude_squared() * b.get_magnitude_squared() ) );
    }

    /**
     * Gets the angle in radians of two given line segments originating at "center" and going to "first" / "second"
     */
    static float_type angle( const vector3t& first, const vector3t& center, const vector3t& second ) {
        return angle( first - center, second - center );
    }

    /**
     * Gets the projection of the vector "source" onto "destination"
     */
    static vector3t project( const vector3t& source, const vector3t& destination ) {
        return destination * ( dot( source, destination ) / destination.get_magnitude_squared() );
    }

    // TODO: where does this go?
    static float_type tetrahedron_volume( const vector3t& pt1, const vector3t& pt2, const vector3t& pt3,
                                          const vector3t& pt4 ) {
        // return std::abs( vector3t::dot( vector3t::cross( (pt2 - pt1),(pt3 - pt1) ), (pt1 - pt4) ) / 6.f );
        //
        // Changed to this since the operator- is not defined yet.
        return ( std::abs( ( ( pt2.y - pt1.y ) * ( pt3.z - pt1.z ) - ( pt2.z - pt1.z ) * ( pt3.y - pt1.y ) ) *
                               ( pt1.x - pt4.x ) +
                           ( ( pt2.z - pt1.z ) * ( pt3.x - pt1.x ) - ( pt2.x - pt1.x ) * ( pt3.z - pt1.z ) ) *
                               ( pt1.y - pt4.y ) +
                           ( ( pt2.x - pt1.x ) * ( pt3.y - pt1.y ) - ( pt2.y - pt1.y ) * ( pt3.x - pt1.x ) ) *
                               ( pt1.z - pt4.z ) ) ) /
               6.f;
    }

    // TODO: remove one of these ( to_normalized )
    static vector3t normalize( const vector3t& v ) { return v.to_normalized(); }

    /**
     *  Normalize the vector, so it points in the same direction, but with
     * magnitude 1.
     *
     *  If the vector cannot be normalized, then it is set to (0, 0, 1), and
     * this function returns false.  This can happen if the vector is zero;
     * very large or very small, such that its get_magnitude() is not finite;
     * or if it is non-finite.
     *
     * @return true if the vector could be normalized, and false otherwise.
     */
    bool normalize() {
        float_type magnitude = get_magnitude();
        *this /= magnitude;

        // In addition to checking that the result is finite, also check that
        // the magnitude is finite.  Non-finite magnitude implies that we
        // could not normalize the vector correctly.
        if( is_finite() && frantic::math::is_finite( magnitude ) ) {
            return true;
        } else {
            set( 0, 0, 1 );
            return false;
        }
    }

    vector3t to_normalized() const {
        vector3t result = *this;
        result.normalize();
        return result;
    }

    vector3t to_floor() const { return vector3t( floor( x ), floor( y ), floor( z ) ); }

    static vector3t to_floor( const vector3t& v ) { return vector3t( floor( v.x ), floor( v.y ), floor( v.z ) ); }

    vector3t to_ceil() const { return vector3t( ceil( x ), ceil( y ), ceil( z ) ); }

    static vector3t to_ceil( const vector3t& v ) { return vector3t( ceil( v.x ), ceil( v.y ), ceil( v.z ) ); }

    bool is_finite() const {
        return frantic::math::is_finite( x ) && frantic::math::is_finite( y ) && frantic::math::is_finite( z );
    }

    bool is_nan() const {
        return frantic::math::is_nan( x ) || frantic::math::is_nan( y ) || frantic::math::is_nan( z );
    }

    bool is_inf() const {
        return frantic::math::is_infinite( x ) || frantic::math::is_infinite( y ) || frantic::math::is_infinite( z );
    }

    static vector3t component_multiply( const vector3t& a, const vector3t& b ) {
        return ( vector3t( a.x * b.x, a.y * b.y, a.z * b.z ) );
    }

    static vector3t component_divide( const vector3t& a, const vector3t& b ) {
        return vector3t( a.x / b.x, a.y / b.y, a.z / b.z );
    }

    static vector3t component_max( const vector3t& a, const vector3t& b ) {
        using std::max;
        return vector3t( max( a.x, b.x ), max( a.y, b.y ), max( a.z, b.z ) );
    }

    static vector3t component_min( const vector3t& a, const vector3t& b ) {
        using std::min;
        return vector3t( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ) );
    }

    //////////////////
    // Operators
    //////////////////

    void set( float_type X, float_type Y, float_type Z ) {
        x = X;
        y = Y;
        z = Z;
    }

    void set( float_type v ) {
        x = v;
        y = v;
        z = v;
    }

    vector3t& operator=( const vector3t& v ) {
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }

    float_type& operator[]( size_t i ) { return ( &x )[i]; }
    const float_type& operator[]( size_t i ) const { return ( &x )[i]; }

    float_type& at( const std::size_t i ) { return ( &x )[i]; }

    const float_type& at( const std::size_t i ) const { return ( &x )[i]; }

    bool is_equal( const vector3t& a ) const { return x == a.x && y == a.y && z == a.z; }

    bool is_equal( const vector3t& a, float tolerance ) const {
        float errorSquared = vector3t::distance_squared( *this, a );
        return errorSquared <=
               tolerance * tolerance * ( std::max )( get_magnitude_squared(), a.get_magnitude_squared() );
    }

    bool is_zero() const { return x == 0 && y == 0 && z == 0; }

    bool is_zero( float_type tolerance ) const { return get_magnitude_squared() <= tolerance * tolerance; }

    bool operator==( const vector3t& a ) const { return x == a.x && y == a.y && z == a.z; }

    bool operator!=( const vector3t& a ) const { return x != a.x || y != a.y || z != a.z; }

    vector3t operator+( const vector3t& b ) const { return vector3t<FloatType>( x + b.x, y + b.y, z + b.z ); }

    vector3t operator-( const vector3t& b ) const { return vector3t<FloatType>( x - b.x, y - b.y, z - b.z ); }

    vector3t operator*( const vector3t& b ) { return vector3t<FloatType>( x * b.x, y * b.y, z * b.z ); }

    vector3t operator/( const vector3t& b ) { return vector3t<FloatType>( x / b.x, y / b.y, z / b.z ); }

    vector3t operator-() const { return vector3t( -x, -y, -z ); }

    vector3t& operator+=( const vector3t& in ) {
        x += in.x;
        y += in.y;
        z += in.z;
        return *this;
    }

    vector3t& operator-=( const vector3t& in ) {
        x -= in.x;
        y -= in.y;
        z -= in.z;
        return *this;
    }

    vector3t& operator*=( const vector3t& in ) {
        x *= in.x;
        y *= in.y;
        z *= in.z;
        return *this;
    }

    vector3t& operator/=( const vector3t& in ) {
        x /= in.x;
        y /= in.y;
        z /= in.z;
        return *this;
    }

    vector3t& operator*=( const float_type& a ) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }

    vector3t& operator/=( const float_type& a ) {
        x /= a;
        y /= a;
        z /= a;
        return *this;
    }

    vector3t operator*( const float_type& k ) const { return vector3t( x * k, y * k, z * k ); }

    vector3t operator/( const float_type& k ) const { return vector3t( x / k, y / k, z / k ); }

    std::string str() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

    static vector3t parse( const std::string& input ) {
        std::vector<std::string> values;
        frantic::strings::split( input, values, "()[], " );
        if( values.size() != 3 )
            throw std::runtime_error( "could not parse vector3t: " + input );
        float_type X, Y, Z;
        try {
            X = boost::lexical_cast<float_type>( values[0] );
            Y = boost::lexical_cast<float_type>( values[1] );
            Z = boost::lexical_cast<float_type>( values[2] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "could not parse vector3t: " + input );
        }

        return vector3t( X, Y, Z );
    }

    friend frantic::graphics::vector3t<FloatType> operator*( const FloatType& k,
                                                             const frantic::graphics::vector3t<FloatType>& a ) {
        return frantic::graphics::vector3t<FloatType>( k * a.x, k * a.y, k * a.z );
    }
};

/**
 * For multiplication, when the second argument's float type is bigger, we want the result
 * of the multiplcation to have that bigger size. This overload enables that.
 */
template <typename FloatType1, typename FloatType2>
typename boost::enable_if_c<( sizeof( FloatType1 ) < sizeof( FloatType2 ) ), vector3t<FloatType2>>::type inline
operator*( const vector3t<FloatType1> t1, const FloatType2& t2 ) {
    return vector3t<FloatType2>( t1 ) * t2;
}

template <typename FloatType1, typename FloatType2>
typename boost::enable_if_c<( sizeof( FloatType2 ) < sizeof( FloatType1 ) ), vector3t<FloatType1>>::type inline
operator*( const FloatType1& t1, const vector3t<FloatType1>& t2 ) {
    return t1 * vector3t<FloatType1>( t2 );
}

template <typename FloatType>
inline FloatType distance( const vector3t<FloatType>& a, const vector3t<FloatType>& b ) {
    return ( a - b ).get_magnitude();
}

template <typename CharType, typename FloatType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const vector3t<FloatType>& v ) {
    std::basic_stringstream<CharType> ss;
    // nonfinite_num_put for a standard representation of non-finite numbers,
    // for example, "inf" instead of the "1.#INF" in Visual Studio
    ss.imbue( std::locale( std::locale::classic(), new boost::math::nonfinite_num_put<CharType> ) );
    ss << "[" << v.x << ", " << v.y << ", " << v.z << "]";
    out << ss.str();
    return out;
}

typedef vector3t<float> vector3f;
typedef vector3t<double> vector3fd;

} // namespace graphics
} // namespace frantic
