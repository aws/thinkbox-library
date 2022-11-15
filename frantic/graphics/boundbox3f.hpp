// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <iostream>

#include <boost/lexical_cast.hpp>

#include <frantic/graphics/size3.hpp>
#include <frantic/graphics/size3f.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace graphics {

template <typename FloatType>
class boundbox3t {
  public:
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;
    typedef size3t<FloatType> size3f_type;

  private:
    vector3f_type m_minimum, m_maximum;

  public:
    //////////////
    // Constructors
    //////////////

    boundbox3t() { set_to_empty(); }

    explicit boundbox3t( const vector3f_type& position )
        : m_minimum( position )
        , m_maximum( position ) {}

    boundbox3t( const vector3f_type& minimum, const vector3f_type& maximum )
        : m_minimum( minimum )
        , m_maximum( maximum ) {}

    boundbox3t( const vector3& minimum, const vector3& maximum )
        : m_minimum( (float_type)minimum.x, (float_type)minimum.y, (float_type)minimum.z )
        , m_maximum( (float_type)maximum.x, (float_type)maximum.y, (float_type)maximum.z ) {}

    boundbox3t( const vector3f_type& minimum, const size3f_type& size )
        : m_minimum( minimum )
        , m_maximum( minimum + size ) {}

    // This array is {xmin, xmax, ymin, ymax, zmin, zmax}
    explicit boundbox3t( const float_type* boundsArray )
        : m_minimum( boundsArray[0], boundsArray[2], boundsArray[4] )
        , m_maximum( boundsArray[1], boundsArray[3], boundsArray[5] ) {}

    boundbox3t( float_type xmin, float_type xmax, float_type ymin, float_type ymax, float_type zmin, float_type zmax )
        : m_minimum( xmin, ymin, zmin )
        , m_maximum( xmax, ymax, zmax ) {}

    /** Implicit conversion from smaller float type is allowed. */
    template <class OtherFloatType>
    boundbox3t( const boundbox3t<OtherFloatType>& rhs,
                typename boost::enable_if_c<( sizeof( OtherFloatType ) < sizeof( FloatType ) ), int*>::type = 0 )
        : m_minimum( rhs.minimum().x, rhs.minimum().y, rhs.minimum().z )
        , m_maximum( rhs.maximum().x, rhs.maximum().y, rhs.maximum().z ) {}

    /** Explicit cast required for conversion from larger float type. */
    template <class OtherFloatType>
    explicit boundbox3t(
        const boundbox3t<OtherFloatType>& rhs,
        typename boost::enable_if_c<( sizeof( OtherFloatType ) > sizeof( FloatType ) ), int*>::type = 0 )
        : m_minimum( (float_type)rhs.minimum().x, (float_type)rhs.minimum().y, (float_type)rhs.minimum().z )
        , m_maximum( (float_type)rhs.maximum().x, (float_type)rhs.maximum().y, (float_type)rhs.maximum().z ) {}

    void set( const vector3f_type& minimum, const vector3f_type& maximum ) {
        m_minimum = minimum;
        m_maximum = maximum;
    }

    void set_to_empty() {
        m_minimum = vector3f_type( ( std::numeric_limits<float_type>::max )() );
        m_maximum = vector3f_type( -( std::numeric_limits<float_type>::max )() );
    }

    void set_to_point( const vector3f_type& position ) {
        m_minimum = position;
        m_maximum = position;
    }

    void set_to_infinite() {
        m_minimum = vector3f_type( -std::numeric_limits<float_type>::infinity() );
        m_maximum = vector3f_type( std::numeric_limits<float_type>::infinity() );
    }

    static boundbox3t from_empty() {
        return boundbox3t( vector3f_type( ( std::numeric_limits<float_type>::max )() ),
                           vector3f_type( -( std::numeric_limits<float_type>::max )() ) );
    }

    static boundbox3t from_center_and_width( vector3f_type center, float_type width ) {
        vector3f_type widthVectorHalf( .5f * width );
        return boundbox3t( center - widthVectorHalf, center + widthVectorHalf );
    }

    static boundbox3t from_random( const boundbox3t& insideWhat ) {
        boundbox3t result = boundbox3t::from_empty();
        result += insideWhat.random_vector();
        result += insideWhat.random_vector();
        return result;
    }

    static boundbox3t infinite() {
        boundbox3t box;
        box.set_to_infinite();
        return box;
    }

    boundbox3t to_expanded_fractional( float_type expandAmount ) const {
        boundbox3t result( *this );
        result.expand_fractional( expandAmount );
        return result;
    }

    vector3f_type get_mapped_from_unit_cube( const vector3f_type& cubeCoord ) const {
        return vector3f_type( m_minimum.x * ( 1 - cubeCoord.x ) + m_maximum.x * cubeCoord.x,
                              m_minimum.y * ( 1 - cubeCoord.y ) + m_maximum.y * cubeCoord.y,
                              m_minimum.z * ( 1 - cubeCoord.z ) + m_maximum.z * cubeCoord.z );
    }

    //********
    // TODO Rename these functions to not use the word voxel, and hence not be confusing with other global voxel space
    // definitions.

    // Gets the voxel coordinate containing a given world position
    vector3f_type get_voxel_coord( vector3f_type coord, size3 voxelSize ) const {
        return vector3f_type( ( coord.x - m_minimum.x ) * voxelSize.xsize() / xsize(),
                              ( coord.y - m_minimum.y ) * voxelSize.ysize() / ysize(),
                              ( coord.z - m_minimum.z ) * voxelSize.zsize() / zsize() );
    }

    vector3f_type from_voxel_coord( vector3f_type voxelCoord, size3 voxelSize ) const {
        return vector3f_type( voxelCoord.x * xsize() / voxelSize.xsize() + m_minimum.x,
                              voxelCoord.y * ysize() / voxelSize.ysize() + m_minimum.y,
                              voxelCoord.z * zsize() / voxelSize.zsize() + m_minimum.z );
    }

    // Gets the voxel bounding box for a given voxel coordinate
    boundbox3t get_voxel_bounds( vector3 coord, size3 voxelSize ) const {
        // Calculate the min and max with the same formula for numerical stability
        return boundbox3t(
            get_mapped_from_unit_cube( vector3f_type( (float_type)coord.x / voxelSize.xsize(),
                                                      (float_type)coord.y / voxelSize.ysize(),
                                                      (float_type)coord.z / voxelSize.zsize() ) ),
            get_mapped_from_unit_cube( vector3f_type( (float_type)( coord.x + 1 ) / voxelSize.xsize(),
                                                      (float_type)( coord.y + 1 ) / voxelSize.ysize(),
                                                      (float_type)( coord.z + 1 ) / voxelSize.zsize() ) ) );
    }
    //********//

    //////////////
    // Member properties
    //////////////

    const vector3f_type& minimum() const { return m_minimum; }

    vector3f_type& minimum() { return m_minimum; }

    const vector3f_type& maximum() const { return m_maximum; }

    vector3f_type& maximum() { return m_maximum; }

    float_type xminimum() const { return m_minimum.x; }

    float_type yminimum() const { return m_minimum.y; }

    float_type zminimum() const { return m_minimum.z; }

    float_type xmaximum() const { return m_maximum.x; }

    float_type ymaximum() const { return m_maximum.y; }

    float_type zmaximum() const { return m_maximum.z; }

    float_type xsize() const { return m_maximum.x - m_minimum.x; }

    float_type ysize() const { return m_maximum.y - m_minimum.y; }

    float_type zsize() const { return m_maximum.z - m_minimum.z; }

    float_type size( int axis ) const { return m_maximum[axis] - m_minimum[axis]; }

    size3f_type size() const { return size3f_type::from_bounds( m_minimum, m_maximum ); }

    bool is_empty() const {
        return m_minimum.x > m_maximum.x || m_minimum.y > m_maximum.y || m_minimum.z > m_maximum.z;
    }

    bool is_infinite() const {
        for( size_t i = 0; i < 3; ++i ) {
            if( m_minimum[i] != -std::numeric_limits<float_type>::infinity() ||
                m_maximum[i] != std::numeric_limits<float_type>::infinity() ) {
                return false;
            }
        }

        return true;
    }

    bool
    is_cube( float_type tolerance = 0.000005f ) const { // TODO: Make tolerance depend on the precision of float_type
        float_type err = fabs( xsize() - ysize() ) + fabs( xsize() - zsize() ) + fabs( ysize() - zsize() );
        //		std::cerr << "Error metric: " << err << std::endl;
        //		std::cerr << "Comparison: " << 0.000001f * get_max_dimension() << std::endl;
        return err < tolerance * get_max_dimension();
    }

    vector3f_type center() const { return 0.5f * ( m_minimum + m_maximum ); }

    float_type get_max_dimension() const { return size().get_max_dimension(); }

    int get_max_dimension_axis() const {
        int axis = 2;

        if( ( m_maximum.x - m_minimum.x ) > ( m_maximum.z - m_minimum.y ) &&
            ( m_maximum.x - m_minimum.x ) > ( m_maximum.z - m_minimum.z ) )
            axis = 0;
        else if( ( m_maximum.y - m_minimum.y ) > ( m_maximum.z - m_minimum.z ) )
            axis = 1;

        return axis;
    }

    float_type get_min_dimension() const { return size().get_min_dimension(); }

    vector3f_type get_diagonal() const { return vector3f_type( xsize(), ysize(), zsize() ); }

    float_type get_diagonal_length() const { return get_diagonal().get_magnitude(); }

    float_type get_volume() const {
        if( is_empty() )
            return 0;
        else
            return xsize() * ysize() * zsize();
    }

    float_type get_surface_area() const {
        if( is_empty() )
            return 0;
        else {
            float_type w = xsize(), h = ysize(), d = zsize();
            return 2 * ( w * h + w * d + h * d );
        }
    }

    vector3f_type get_corner( int i ) const {
        switch( i ) {
        case 0:
            return m_minimum;
        case 1:
            return vector3f_type( m_maximum.x, m_minimum.y, m_minimum.z );
        case 2:
            return vector3f_type( m_minimum.x, m_maximum.y, m_minimum.z );
        case 3:
            return vector3f_type( m_maximum.x, m_maximum.y, m_minimum.z );
        case 4:
            return vector3f_type( m_minimum.x, m_minimum.y, m_maximum.z );
        case 5:
            return vector3f_type( m_maximum.x, m_minimum.y, m_maximum.z );
        case 6:
            return vector3f_type( m_minimum.x, m_maximum.y, m_maximum.z );
        case 7:
            return m_maximum;
        default:
            throw std::runtime_error(
                "boundbox3t.get_corner: Tried to get an invalid corner (must be between 0 and 7)" );
        }
    }

    vector3f_type get_corner_0() const { return m_minimum; }
    vector3f_type get_corner_1() const { return vector3f_type( m_maximum.x, m_minimum.y, m_minimum.z ); }
    vector3f_type get_corner_2() const { return vector3f_type( m_minimum.x, m_maximum.y, m_minimum.z ); }
    vector3f_type get_corner_3() const { return vector3f_type( m_maximum.x, m_maximum.y, m_minimum.z ); }
    vector3f_type get_corner_4() const { return vector3f_type( m_minimum.x, m_minimum.y, m_maximum.z ); }
    vector3f_type get_corner_5() const { return vector3f_type( m_maximum.x, m_minimum.y, m_maximum.z ); }
    vector3f_type get_corner_6() const { return vector3f_type( m_minimum.x, m_maximum.y, m_maximum.z ); }
    vector3f_type get_corner_7() const { return m_maximum; }

    // This tests whether the point is in the bounding box, considered as a closed set.
    template <class OtherVecType>
    bool contains( const OtherVecType& worldPos ) const {
        return m_minimum.x <= worldPos.x && worldPos.x <= m_maximum.x && m_minimum.y <= worldPos.y &&
               worldPos.y <= m_maximum.y && m_minimum.z <= worldPos.z && worldPos.z <= m_maximum.z;
    }

    // This tests whether the point is in the bounding box, considered as an open set.
    template <class OtherVecType>
    bool contains_as_open_set( const OtherVecType& worldPos ) const {
        return m_minimum.x < worldPos.x && worldPos.x < m_maximum.x && m_minimum.y < worldPos.y &&
               worldPos.y < m_maximum.y && m_minimum.z < worldPos.z && worldPos.z < m_maximum.z;
    }

    bool contains( const boundbox3t& worldBox ) const {
        return contains( worldBox.minimum() ) && contains( worldBox.maximum() );
    }

    vector3f_type random_vector() const {
        if( is_empty() )
            throw std::runtime_error( "boundbox3t.random_vector: Cannot get a random vector from an empty box." );
        return minimum() + vector3f_type::from_random() * size();
    }

    template <class RandomNumberGenerator>
    vector3f_type random_vector( RandomNumberGenerator& rng ) {
        if( is_empty() )
            throw std::runtime_error( "boundbox3t.random_vector: Cannot get a random vector from an empty box." );
        return minimum() + vector3f_type::from_random( rng ) * size();
        /*
        // TODO: Use this instead?
        vector3f_type unitRandom = vector3f_type::from_random(rng);
        return vector3f_type( (1-unitRandom.x)*m_minimum.x + unitRandom.x*m_maximum.x,
          (1-unitRandom.y)*m_minimum.y + unitRandom.y*m_maximum.y,
          (1-unitRandom.z)*m_minimum.z + unitRandom.z*m_maximum.z );
        */
    }

    // This computes the support of this bounding box in the direction d.  The
    // support is the (not necessarily unique) point in the bounding box which
    // is furthest along in the direction d.  That is,
    // (p in this) such that { dot(p,d) >= p' for all (p' in this) }.
    // Mathematically, the support would be the set of all such points, but we just
    // return one.
    vector3f_type support( const vector3f_type& d ) const {
        return vector3f_type( d.x > 0 ? m_maximum.x : m_minimum.x, d.y > 0 ? m_maximum.y : m_minimum.y,
                              d.z > 0 ? m_maximum.z : m_minimum.z );
    }

    vector3f_type clamp( vector3f_type v ) const {
        if( is_empty() ) {
            throw std::runtime_error( "boundbox3t: Cannot clamp vector " + boost::lexical_cast<std::string>( v ) +
                                      " because box is empty" );
        }
        if( v.x < m_minimum.x )
            v.x = m_minimum.x;
        if( v.y < m_minimum.y )
            v.y = m_minimum.y;
        if( v.z < m_minimum.z )
            v.z = m_minimum.z;
        if( v.x > m_maximum.x )
            v.x = m_maximum.x;
        if( v.y > m_maximum.y )
            v.y = m_maximum.y;
        if( v.z > m_maximum.z )
            v.z = m_maximum.z;
        return v;
    }

    // Clamps the vector without checking for emptiness. Note that when the boundbox3t is empty, the result
    // given won't be well-defined.
    vector3f_type clamp_nothrow( vector3f_type v ) const {
        if( v.x < m_minimum.x )
            v.x = m_minimum.x;
        if( v.y < m_minimum.y )
            v.y = m_minimum.y;
        if( v.z < m_minimum.z )
            v.z = m_minimum.z;
        if( v.x > m_maximum.x )
            v.x = m_maximum.x;
        if( v.y > m_maximum.y )
            v.y = m_maximum.y;
        if( v.z > m_maximum.z )
            v.z = m_maximum.z;
        return v;
    }

    // returns a distance function for this box (whose level set is the edge of the box)
    float_type distance_function( const vector3f_type& v ) const {
        if( is_empty() ) {
            throw std::runtime_error( "boundbox3f: Cannot get distance function for vector " +
                                      boost::lexical_cast<std::string>( v ) + " because box is empty" );
        }
        if( contains( v ) ) {
            float_type d = v.x - m_minimum.x;
            d = ( std::min )( d, v.y - m_minimum.y );
            d = ( std::min )( d, v.z - m_minimum.z );
            d = ( std::min )( d, m_maximum.x - v.x );
            d = ( std::min )( d, m_maximum.y - v.y );
            d = ( std::min )( d, m_maximum.z - v.z );
            return -d;
        } else {
            return ( clamp( v ) - v ).get_magnitude();
        }
    }

    vector3f_type get_nearest_surface_point( const vector3f_type& v ) const {
        if( is_empty() ) {
            throw std::runtime_error( "boundbox3t.get_nearest_surface_point: Cannot get the nearest point to vector " +
                                      boost::lexical_cast<std::string>( v ) + " because box is empty" );
        }
        if( contains( v ) ) {
            vector3f_type result;
            float_type d;
            if( v.x - m_minimum.x > m_maximum.x - v.x ) {
                result = vector3f_type( m_maximum.x, v.y, v.z );
                d = m_maximum.x - v.x;
            } else {
                result = vector3f_type( m_minimum.x, v.y, v.z );
                d = v.x - m_minimum.x;
            }
            if( v.y - m_minimum.y > m_maximum.y - v.y ) {
                if( m_maximum.y - v.y < d ) {
                    result = vector3f_type( v.x, m_maximum.y, v.z );
                    d = m_maximum.y - v.y;
                }
            } else {
                if( v.y - m_minimum.y < d ) {
                    result = vector3f_type( v.x, m_minimum.y, v.z );
                    d = v.y - m_minimum.y;
                }
            }
            if( v.z - m_minimum.z > m_maximum.z - v.z ) {
                if( m_maximum.z - v.z < d ) {
                    result = vector3f_type( v.x, v.y, m_maximum.z );
                    d = m_maximum.z - v.z;
                }
            } else {
                if( v.z - m_minimum.z < d ) {
                    result = vector3f_type( v.x, v.y, m_minimum.z );
                    d = v.z - m_minimum.z;
                }
            }

            return result;
        } else {
            return clamp( v );
        }
    }

    // This intersects a line segment with the bounding box, and returns both the intersection points and the number of
    // intersections (0, 1, or 2)
    int intersect_with_line_segment( const vector3f_type& start, const vector3f_type& end,
                                     vector3f_type& outFirstIntersection, vector3f_type& outSecondIntersection ) const;

    // Tests whether the box is intersecting another box.
    bool is_intersecting( const boundbox3t& box ) const {
        return !( m_minimum.x > box.m_maximum.x || m_maximum.x < box.m_minimum.x || m_minimum.y > box.m_maximum.y ||
                  m_maximum.y < box.m_minimum.y || m_minimum.z > box.m_maximum.z || m_maximum.z < box.m_minimum.z );
    }

    // Triangle Box intersection test as described by Tomas Akenine-Moller in his 2001 paper, "Fast 3D Triangle-Box
    // Overlap Testing."
    bool is_intersecting_triangle( const vector3f_type& vert0, const vector3f_type& vert1,
                                   const vector3f_type& vert2 ) const;

    // This function finds the intersection of the triangle and the bound box, returning the result in another bound
    // box. It is used by the kd-tree construction code to compute perfect splits.  An example of where a more naive
    // algorithm would erroneously include a triangle when it shouldn't is as follows:
    //
    // |-|   /
    // | |  /
    // |-| /
    //    /
    //   /
    //  /
    //
    // NOTE: This is much slower than the is_intersecting_triangle function, and you may want to use that function to
    //       first test for intersection before calling this one.
    bool intersect_with_triangle( const vector3f_type& vert0, const vector3f_type& vert1, const vector3f_type& vert2,
                                  boundbox3t& outIntersectedBounds ) const;

    float_type volume() const { return xsize() * ysize() * zsize(); }

    bool is_volume_intersecting_line_segment( const vector3f_type& start, const vector3f_type& end ) const;

    //////////////
    // Operators
    //////////////

    void expand_fractional( float_type fraction ) {
        float_type expandAmount = .5f * fraction * get_max_dimension();
        expand( expandAmount );
    }

    void expand( float_type worldDistance ) {
        if( !is_empty() ) {
            m_minimum.x -= worldDistance;
            m_minimum.y -= worldDistance;
            m_minimum.z -= worldDistance;
            m_maximum.x += worldDistance;
            m_maximum.y += worldDistance;
            m_maximum.z += worldDistance;
        }
    }

    void expand( const vector3f_type& worldDistance ) {
        if( !is_empty() ) {
            m_minimum -= worldDistance;
            m_maximum += worldDistance;
        }
    }

    void expand_to_cube() {
        vector3f_type c = center();
        float_type w = get_max_dimension() / 2;
        m_minimum = vector3f_type( c.x - w, c.y - w, c.z - w );
        m_maximum = vector3f_type( c.x + w, c.y + w, c.z + w );
    }

    boundbox3t& operator+=( const boundbox3t& box ) {
        if( is_empty() ) {
            m_minimum = box.minimum();
            m_maximum = box.maximum();
        } else if( !box.is_empty() ) {
            if( m_minimum.x > box.m_minimum.x )
                m_minimum.x = box.m_minimum.x;
            if( m_maximum.x < box.m_maximum.x )
                m_maximum.x = box.m_maximum.x;
            if( m_minimum.y > box.m_minimum.y )
                m_minimum.y = box.m_minimum.y;
            if( m_maximum.y < box.m_maximum.y )
                m_maximum.y = box.m_maximum.y;
            if( m_minimum.z > box.m_minimum.z )
                m_minimum.z = box.m_minimum.z;
            if( m_maximum.z < box.m_maximum.z )
                m_maximum.z = box.m_maximum.z;
        }
        return *this;
    }

    boundbox3t& operator+=( const vector3f_type& worldPos ) {
        if( is_empty() ) {
            m_minimum = worldPos;
            m_maximum = worldPos;
        } else {
            if( m_minimum.x > worldPos.x )
                m_minimum.x = worldPos.x;
            else if( m_maximum.x < worldPos.x )
                m_maximum.x = worldPos.x;
            if( m_minimum.y > worldPos.y )
                m_minimum.y = worldPos.y;
            else if( m_maximum.y < worldPos.y )
                m_maximum.y = worldPos.y;
            if( m_minimum.z > worldPos.z )
                m_minimum.z = worldPos.z;
            else if( m_maximum.z < worldPos.z )
                m_maximum.z = worldPos.z;
        }
        return *this;
    }

    void intersect_with( const boundbox3t& bounds ) {
        m_maximum.x = ( std::min )( m_maximum.x, bounds.m_maximum.x );
        m_maximum.y = ( std::min )( m_maximum.y, bounds.m_maximum.y );
        m_maximum.z = ( std::min )( m_maximum.z, bounds.m_maximum.z );
        m_minimum.x = ( std::max )( m_minimum.x, bounds.m_minimum.x );
        m_minimum.y = ( std::max )( m_minimum.y, bounds.m_minimum.y );
        m_minimum.z = ( std::max )( m_minimum.z, bounds.m_minimum.z );
    }

    // performs a modulus operation on the given point to place it within the bound box
    vector3f_type mod_point( const vector3f_type& point ) const {
        vector3f_type divisor = m_maximum - m_minimum;

        vector3f_type result;
        result.x = (float_type)fmod( point.x, divisor.x );
        result.y = (float_type)fmod( point.y, divisor.y );
        result.z = (float_type)fmod( point.z, divisor.z );

        return result;
    }

    std::pair<boundbox3t, boundbox3t> split( int axis, float_type kPlane ) const {
        if( m_minimum[axis] > kPlane || kPlane > m_maximum[axis] )
            throw std::runtime_error( "boundbox3t::split() - Split plane: " + boost::lexical_cast<std::string>( axis ) +
                                      "," + boost::lexical_cast<std::string>( kPlane ) +
                                      " did not intersect the boundbox: " + str() );

        return std::make_pair( boundbox3t( m_minimum, m_maximum.copy().replace_element( axis, kPlane ) ),
                               boundbox3t( m_minimum.copy().replace_element( axis, kPlane ), m_maximum ) );
    }

    // takes a point and performs a modulus operation to place the point within the boundbox
    // the return container will then have any points that are modulo permutations of the original point
    // that are within the given threshold of the faces of the boundbox
    // The modded original point is also returned and will be the last point in the container
    // This could be used to create cyclic boundary conditions
    void compute_mod_boundary_points( const vector3f_type& point, float_type boundaryThreshold,
                                      std::vector<vector3f_type>& outBoundaryPoints ) const;

    std::string str() const;

    static boundbox3t parse( const std::string& input ) {
        std::string s = input;

        s = frantic::strings::string_replace( s, "boundbox3", "" ); // magic boundbox3 into boundbox3f, enjoy
        s = frantic::strings::string_replace( s, "boundbox", "" );

        std::vector<std::string> tokens;
        frantic::strings::split( s, tokens, "()[],\n\r " );

        try {
            if( tokens.size() == 6 ) {
                float_type bounds[6];

                for( int i = 0; i < 6; ++i )
                    bounds[i] = boost::lexical_cast<float_type>( tokens[i] );

                return boundbox3t( bounds );
            } else {
                throw std::runtime_error( "Could not parse bounding box due to incorrect number of tokens:" + input );
            }
        } catch( const boost::bad_lexical_cast& ) {
            throw std::runtime_error( "could not parse bounding box due to invalid tokens: " + input );
        }
    }
};

template <class FloatType>
inline bool operator==( const boundbox3t<FloatType>& lhs, const boundbox3t<FloatType>& rhs ) {
    return lhs.minimum() == rhs.minimum() && lhs.maximum() == rhs.maximum();
}

template <class FloatType>
inline bool operator!=( const boundbox3t<FloatType>& lhs, const boundbox3t<FloatType>& rhs ) {
    return !( lhs == rhs );
}

template <class CharType, class FloatType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const boundbox3t<FloatType>& box ) {
    out << "(box3f " << box.minimum() << " " << box.maximum() << " )";
    return out;
}

template <class FloatType>
std::string boundbox3t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

typedef boundbox3t<float> boundbox3f;
typedef boundbox3t<double> boundbox3fd;

namespace detail {
// A mechanism to get the more complicated/slow to compile code into the .cpp file
bool is_intersecting_triangle( const boundbox3f& box, const vector3f& vert0, const vector3f& vert1,
                               const vector3f& vert2 );
bool is_intersecting_triangle( const boundbox3fd& box, const vector3fd& vert0, const vector3fd& vert1,
                               const vector3fd& vert2 );
bool intersect_with_triangle( const boundbox3f& box, const vector3f& vert0, const vector3f& vert1,
                              const vector3f& vert2, boundbox3f& outIntersectedBounds );
bool intersect_with_triangle( const boundbox3fd& box, const vector3fd& vert0, const vector3fd& vert1,
                              const vector3fd& vert2, boundbox3fd& outIntersectedBounds );
int intersect_with_line_segment( const boundbox3f& box, const vector3f& start, const vector3f& end,
                                 vector3f& outFirstIntersection, vector3f& outSecondIntersection );
int intersect_with_line_segment( const boundbox3fd& box, const vector3fd& start, const vector3fd& end,
                                 vector3fd& outFirstIntersection, vector3fd& outSecondIntersection );
bool is_volume_intersecting_line_segment( const boundbox3f& box, const vector3f& start, const vector3f& end );
bool is_volume_intersecting_line_segment( const boundbox3fd& box, const vector3fd& start, const vector3fd& end );
void compute_mod_boundary_points( const boundbox3f& box, const vector3f& point, float boundaryThreshold,
                                  std::vector<vector3f>& outBoundaryPoints );
void compute_mod_boundary_points( const boundbox3fd& box, const vector3fd& point, double boundaryThreshold,
                                  std::vector<vector3fd>& outBoundaryPoints );
} // namespace detail

template <class FloatType>
inline bool boundbox3t<FloatType>::is_intersecting_triangle( const vector3f_type& vert0, const vector3f_type& vert1,
                                                             const vector3f_type& vert2 ) const {
    return detail::is_intersecting_triangle( *this, vert0, vert1, vert2 );
}

template <class FloatType>
inline bool boundbox3t<FloatType>::intersect_with_triangle( const vector3f_type& vert0, const vector3f_type& vert1,
                                                            const vector3f_type& vert2,
                                                            boundbox3t& outIntersectedBounds ) const {
    return detail::intersect_with_triangle( *this, vert0, vert1, vert2, outIntersectedBounds );
}

template <class FloatType>
inline int boundbox3t<FloatType>::intersect_with_line_segment( const vector3f_type& start, const vector3f_type& end,
                                                               vector3f_type& outFirstIntersection,
                                                               vector3f_type& outSecondIntersection ) const {
    return detail::intersect_with_line_segment( *this, start, end, outFirstIntersection, outSecondIntersection );
}

template <class FloatType>
inline bool boundbox3t<FloatType>::is_volume_intersecting_line_segment( const vector3f_type& start,
                                                                        const vector3f_type& end ) const {
    return detail::is_volume_intersecting_line_segment( *this, start, end );
}

template <class FloatType>
inline void boundbox3t<FloatType>::compute_mod_boundary_points( const vector3f_type& point,
                                                                float_type boundaryThreshold,
                                                                std::vector<vector3f_type>& outBoundaryPoints ) const {
    return detail::compute_mod_boundary_points( *this, point, boundaryThreshold, outBoundaryPoints );
}

} // namespace graphics
} // namespace frantic
