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
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3.hpp>

namespace frantic {
namespace graphics {

/**
 * A 3-dimensional integer bounding box.
 */
class boundbox3 {
    vector3 m_minimum, m_maximum;

  public:
    //////////////
    // Constructors
    //////////////

    /**
     * The default boundbox3 constructor initializes to an empty boundbox3.
     */
    boundbox3() { set_to_empty(); }

    explicit boundbox3( const vector3& position ) {
        m_minimum = position;
        m_maximum = position;
    }

    boundbox3( const vector3& minimum, const vector3& maximum ) {
        m_minimum = minimum;
        m_maximum = maximum;
    }

    boundbox3( const vector3& minimum, const size3& size ) {
        m_minimum = minimum;
        m_maximum = minimum + size - vector3( 1 );
    }

    boundbox3( int xmin, int xmax, int ymin, int ymax, int zmin, int zmax )
        : m_minimum( xmin, ymin, zmin )
        , m_maximum( xmax, ymax, zmax ) {}

    /**
     * This constructor intializes to the data in the provided array.
     *
     * @param  boundsArray  The input bounds, as an array with order {xmin, xmax, ymin, ymax, zmin, zmax}.
     */
    explicit boundbox3( const int* boundsArray )
        : m_minimum( boundsArray[0], boundsArray[2], boundsArray[4] )
        , m_maximum( boundsArray[1], boundsArray[3], boundsArray[5] ) {}

    void set( const vector3& minimum, const vector3& maximum ) {
        m_minimum = minimum;
        m_maximum = maximum;
    }

    void set( int xmin, int xmax, int ymin, int ymax, int zmin, int zmax ) {
        m_minimum.x = xmin;
        m_minimum.y = ymin;
        m_minimum.z = zmin;
        m_maximum.x = xmax;
        m_maximum.y = ymax;
        m_maximum.z = zmax;
    }

    void set_to_empty() {
        m_minimum = vector3( 1, 1, 1 );
        m_maximum = vector3( 0, 0, 0 );
    }

    void set_to_point( const vector3& position ) {
        m_minimum = position;
        m_maximum = position;
    }

    static boundbox3 from_empty() { return boundbox3( vector3( 1, 1, 1 ), vector3( 0, 0, 0 ) ); }

    vector3 random_vector() const {
        if( is_empty() )
            throw std::runtime_error( "boundbox3.random_vector: Cannot get a random vector from an empty box." );
        return vector3( m_minimum.x + rand() % xsize(), m_minimum.y + rand() % ysize(),
                        m_minimum.z + rand() % zsize() );
    }

    //////////////
    // Member properties
    //////////////

    const vector3& minimum() const { return m_minimum; }

    vector3& minimum() { return m_minimum; }

    const vector3& maximum() const { return m_maximum; }

    vector3& maximum() { return m_maximum; }

    int xminimum() const { return m_minimum.x; }

    int yminimum() const { return m_minimum.y; }

    int zminimum() const { return m_minimum.z; }

    int xmaximum() const { return m_maximum.x; }

    int ymaximum() const { return m_maximum.y; }

    int zmaximum() const { return m_maximum.z; }

    int xsize() const { return m_maximum.x - m_minimum.x + 1; }

    int ysize() const { return m_maximum.y - m_minimum.y + 1; }

    int zsize() const { return m_maximum.z - m_minimum.z + 1; }

    size3 size() const { return size3( xsize(), ysize(), zsize() ); }

    bool is_empty() const {
        return ( m_minimum.x > m_maximum.x ) || ( m_minimum.y > m_maximum.y ) || ( m_minimum.z > m_maximum.z );
    }

    bool is_cube() const { return xsize() == ysize() && ysize() == zsize(); }

    // The center is defined to be the center of the set of voxels spanned by this box.  That means
    // for a bounding box where the minimum and maximum are both (0,0,0), the center will be at (0.5,0.5,0.5).
    vector3f center() const { return 0.5f * vector3f( m_minimum + m_maximum + vector3( 1 ) ); }

    // This returns the maximum size along any of the three axes
    int get_max_dimension() const { return size().get_max_dimension(); }

    // This returns which of the three axes is biggest
    int get_max_dimension_axis() const {
        int axis = 2;

        if( ( m_maximum.x - m_minimum.x ) > ( m_maximum.z - m_minimum.y ) &&
            ( m_maximum.x - m_minimum.x ) > ( m_maximum.z - m_minimum.z ) )
            axis = 0;
        else if( ( m_maximum.y - m_minimum.y ) > ( m_maximum.z - m_minimum.z ) )
            axis = 1;

        return axis;
    }

    int get_min_dimension() const { return size().get_min_dimension(); }

    int get_volume() const {
        if( is_empty() )
            return 0;
        else
            return xsize() * ysize() * zsize();
    }

    vector3 get_corner( int i ) const {
        switch( i ) {
        case 0:
            return vector3( m_minimum.x, m_minimum.y, m_minimum.z );
        case 1:
            return vector3( m_minimum.x, m_minimum.y, m_maximum.z );
        case 2:
            return vector3( m_minimum.x, m_maximum.y, m_minimum.z );
        case 3:
            return vector3( m_minimum.x, m_maximum.y, m_maximum.z );
        case 4:
            return vector3( m_maximum.x, m_minimum.y, m_minimum.z );
        case 5:
            return vector3( m_maximum.x, m_minimum.y, m_maximum.z );
        case 6:
            return vector3( m_maximum.x, m_maximum.y, m_minimum.z );
        case 7:
            return vector3( m_maximum.x, m_maximum.y, m_maximum.z );
        default:
            throw std::runtime_error(
                "boundbox3.get_corner: Tried to get an invalid corner (must be between 0 and 7)" );
        }
    }

    // The bounding box goes from integer coordinate m_minimum to m_maximum.
    // In a similar manner to how we treat images, the center of "voxel" (0,0,0) is (0.5,0.5,0.5), and its extents are
    // half-open from (0,0,0) to (1,1,1).
    bool contains( const vector3& coord ) const {
        return m_minimum.x <= coord.x && coord.x <= m_maximum.x && m_minimum.y <= coord.y && coord.y <= m_maximum.y &&
               m_minimum.z <= coord.z && coord.z <= m_maximum.z;
    }

    // For floating point coordinates, we consider "voxel" (m,n,r) as being the half-open cube from (m,n,r) to
    // (m+1,n+1,r+1). This way, if we consider the union of all the voxels, we exactly reconstruct the cartesion space.
    // For example, coordinate (0.5,0.5,0.0) belongs to voxel (0,0,0) and is not in voxel (0,0,-1).
    // NOTE: If we consider the mapping from integer coordinates to voxel centers (i.e. add vector3f(0.5)),
    // this contains operator is identical to the above vector3 version.
    bool contains( const vector3f& coord ) const {
        return m_minimum.x <= coord.x && coord.x < m_maximum.x + 1 && m_minimum.y <= coord.y &&
               coord.y < m_maximum.y + 1 && m_minimum.z <= coord.z && coord.z < m_maximum.z + 1;
    }

    bool contains( const boundbox3& worldBox ) const {
        return contains( worldBox.minimum() ) && contains( worldBox.maximum() );
    }

    // Tests whether the box is intersecting another box.
    bool is_intersecting( const boundbox3& box ) const {
        return !( m_minimum.x > box.m_maximum.x || m_maximum.x < box.m_minimum.x || m_minimum.y > box.m_maximum.y ||
                  m_maximum.y < box.m_minimum.y || m_minimum.z > box.m_maximum.z || m_maximum.z < box.m_minimum.z );
    }

    vector3 clamp( vector3 v ) const {
        if( is_empty() ) {
            throw std::runtime_error( "boundbox3: Cannot clamp vector " + boost::lexical_cast<std::string>( v ) +
                                      " because box is empty" );
        }
        if( v.x < m_minimum.x )
            v.x = m_minimum.x;
        if( v.y < m_minimum.y )
            v.y = m_minimum.y;
        if( v.z < m_minimum.z )
            v.z = m_minimum.z;
        if( v.x >= m_maximum.x )
            v.x = m_maximum.x;
        if( v.y >= m_maximum.y )
            v.y = m_maximum.y;
        if( v.z >= m_maximum.z )
            v.z = m_maximum.z;
        return v;
    }

    //////////////
    // Operators
    //////////////

    std::string str() const;

    static boundbox3 parse( const std::string& input ) {
        std::string s = input;

        s = frantic::strings::string_replace( s, "boundbox3", "" );
        s = frantic::strings::string_replace( s, "boundbox", "" );
        s = frantic::strings::string_replace( s, "box3", "" );
        s = frantic::strings::string_replace( s, "vec3", "" );

        std::vector<std::string> tokens;
        frantic::strings::split( s, tokens, "()[],\n\r " );

        try {
            if( tokens.size() == 6 ) {
                int bounds[6];

                for( int i = 0; i < 6; ++i )
                    bounds[i] = boost::lexical_cast<int>( tokens[i] );

                return boundbox3( vector3( &bounds[0] ), vector3( &bounds[3] ) );
            } else {
                throw std::runtime_error( "Could not parse bounding box due to incorrect number of tokens:" + input );
            }
        } catch( const boost::bad_lexical_cast& ) {
            throw std::runtime_error( "could not parse bounding box due to invalid tokens: " + input );
        }
    }

    void expand( int worldDistance ) { expand( vector3( worldDistance ) ); }

    void expand( const vector3& worldDistance ) {
        if( !is_empty() ) {
            m_minimum -= worldDistance;
            m_maximum += worldDistance;
        }
    }

    boundbox3& operator+=( const boundbox3& box ) {
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

    boundbox3& operator+=( const vector3& worldPos ) {
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

    /**
     * Expand the bounding box so that it includes this X interval.
     *
     * @param  xminimum  The minimum value of the input X interval.
     * @param  xmaximum  The maximum value of the input X interval.
     */
    void include_x_interval( int xminimum, int xmaximum ) {
        if( m_minimum.x > m_maximum.x ) {
            m_minimum.x = xminimum;
            m_maximum.x = xmaximum;
        } else if( xminimum <= xmaximum ) {
            if( xminimum < m_minimum.x )
                m_minimum.x = xminimum;
            if( xmaximum > m_maximum.x )
                m_maximum.x = xmaximum;
        }
    }

    /**
     * Expand the bounding box so that it includes this X value.
     *
     * @param  x  The X value to include.
     */
    void include_x( int x ) {
        if( m_minimum.x > m_maximum.x ) {
            m_minimum.x = x;
            m_maximum.x = x;
        } else if( x < m_minimum.x ) {
            m_minimum.x = x;
        } else if( x > m_maximum.x ) {
            m_maximum.x = x;
        }
    }

    /**
     * Expand the bounding box so that it includes this Y interval.
     *
     * @param  yminimum  The minimum value of the input Y interval.
     * @param  ymaximum  The maximum value of the input Y interval.
     */
    void include_y_interval( int yminimum, int ymaximum ) {
        if( m_minimum.y > m_maximum.y ) {
            m_minimum.y = yminimum;
            m_maximum.y = ymaximum;
        } else if( yminimum <= ymaximum ) {
            if( yminimum < m_minimum.y )
                m_minimum.y = yminimum;
            if( ymaximum > m_maximum.y )
                m_maximum.y = ymaximum;
        }
    }

    /**
     * Expand the bounding box so that it includes this Y value.
     *
     * @param  y  The Y value to include.
     */
    void include_y( int y ) {
        if( m_minimum.y > m_maximum.y ) {
            m_minimum.y = y;
            m_maximum.y = y;
        } else if( y < m_minimum.y ) {
            m_minimum.y = y;
        } else if( y > m_maximum.y ) {
            m_maximum.y = y;
        }
    }

    /**
     * Expand the bounding box so that it includes this Z interval.
     *
     * @param  zminimum  The minimum value of the input Z interval.
     * @param  zmaximum  The maximum value of the input Z interval.
     */
    void include_z_interval( int zminimum, int zmaximum ) {
        if( m_minimum.z > m_maximum.z ) {
            m_minimum.z = zminimum;
            m_maximum.z = zmaximum;
        } else if( zminimum <= zmaximum ) {
            if( zminimum < m_minimum.z )
                m_minimum.z = zminimum;
            if( zmaximum > m_maximum.z )
                m_maximum.z = zmaximum;
        }
    }

    /**
     * Expand the bounding box so that it includes this Z value.
     *
     * @param  z  The Z value to include.
     */
    void include_z( int z ) {
        if( m_minimum.z > m_maximum.z ) {
            m_minimum.z = z;
            m_maximum.z = z;
        } else if( z < m_minimum.z ) {
            m_minimum.z = z;
        } else if( z > m_maximum.z ) {
            m_maximum.z = z;
        }
    }

    void intersect_with( const boundbox3& bounds ) {
        if( bounds.m_maximum.x < m_maximum.x )
            m_maximum.x = bounds.m_maximum.x;
        if( bounds.m_maximum.y < m_maximum.y )
            m_maximum.y = bounds.m_maximum.y;
        if( bounds.m_maximum.z < m_maximum.z )
            m_maximum.z = bounds.m_maximum.z;
        if( bounds.m_minimum.x > m_minimum.x )
            m_minimum.x = bounds.m_minimum.x;
        if( bounds.m_minimum.y > m_minimum.y )
            m_minimum.y = bounds.m_minimum.y;
        if( bounds.m_minimum.z > m_minimum.z )
            m_minimum.z = bounds.m_minimum.z;
    }

    bool operator==( const boundbox3& box ) const { return m_minimum == box.m_minimum && m_maximum == box.m_maximum; }

    bool operator!=( const boundbox3& box ) const { return m_minimum != box.m_minimum || m_maximum != box.m_maximum; }
};

template <class CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const boundbox3& box ) {
    out << "(box3 [" << box.minimum().x << ", " << box.minimum().y << ", " << box.minimum().z << "] ["
        << box.maximum().x << ", " << box.maximum().y << ", " << box.maximum().z << "] )";
    return out;
}

inline std::string boundbox3::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

} // namespace graphics
} // namespace frantic
