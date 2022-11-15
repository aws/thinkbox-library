// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics2d/boundrect2t.hpp>
#include <frantic/graphics2d/size2f.hpp>
#include <frantic/graphics2d/vector2f.hpp>

namespace frantic {
namespace graphics2d {

class boundrect2f : public boundrect2t<vector2f, size2f, boundrect2f> {
  public:
    //////////////
    // Constructors
    //////////////

    boundrect2f()
        : boundrect2t<vector2f, size2f, boundrect2f>( vector2f( 0, 0 ), vector2f( 0, 0 ) ) {}
    boundrect2f( float xmin, float xmax, float ymin, float ymax )
        : boundrect2t<vector2f, size2f, boundrect2f>( vector2f( xmin, ymin ), vector2f( xmax, ymax ) ) {}
    explicit boundrect2f( const vector2f& position )
        : boundrect2t<vector2f, size2f, boundrect2f>( position ) {}
    boundrect2f( const vector2f& min, const vector2f& max )
        : boundrect2t<vector2f, size2f, boundrect2f>( min, max ) {}
    boundrect2f( const vector2f& min, const size2f& size )
        : boundrect2t<vector2f, size2f, boundrect2f>(
              min, vector2f( min.get_x() + size.get_xsize(), min.get_y() + size.get_ysize() ) ) {}

    value_type xsize() const { return m_maximum.x - m_minimum.x; }
    value_type ysize() const { return m_maximum.y - m_minimum.y; }

    void set_to_point( const vector2f& position ) {
        m_minimum = position;
        m_maximum = position;
    }
    void set_empty() {
        m_minimum = vector2f( 0, 0 );
        m_maximum = vector2f( 0, 0 );
    }

    bool is_area_intersecting_line_segment( const vector2f& start, const vector2f& end ) const {
        if( is_empty() ) {
            return false;
        } else if( contains( start ) || contains( end ) ) {
            return true;
        } else if( ( start.x < m_minimum.x && end.x < m_minimum.x ) ||
                   ( start.y < m_minimum.y && end.y < m_minimum.y ) ||
                   ( start.x > m_maximum.x && end.x > m_maximum.x ) ||
                   ( start.y > m_maximum.y && end.y > m_maximum.y ) ) {
            // A number of special cases of guaranteed non-intersection, which were very numerically
            // unstable in the test below.
            return false;
        } else if( start.x >= m_minimum.x && start.x <= m_maximum.x && end.x >= m_minimum.x && end.x <= m_maximum.x ) {
            // The test below was also missing true intersections like these.
            return ( ( start.y >= m_minimum.y && end.y <= m_maximum.y ) ||
                     ( end.y >= m_minimum.y && start.y <= m_maximum.y ) );
        } else if( start.y >= m_minimum.y && start.y <= m_maximum.y && end.y >= m_minimum.y && end.y <= m_maximum.y ) {
            // The test below was also missing true intersections like these.
            return ( ( start.x >= m_minimum.x && end.x <= m_maximum.x ) ||
                     ( end.x >= m_minimum.x && start.x <= m_maximum.x ) );
        } else {
            vector2f delta = end - start;

            // If the delta is numerically unstable, return that it doesn't intersect
            //			float absoluteMagnitude =  ( (std::max)( start.max_abs_component(), end.max_abs_component()
            //)
            //); 			if( delta.max_abs_component() < 0.001f * absoluteMagnitude ) 				return
            //false;

            float maximumOfNearDistances = 0, minimumOfFarDistances = 1;

            bool reasonableDeltaX = fabsf( delta.x ) > ( std::max )( fabsf( start.x ), fabsf( end.x ) ) * 0.00001f;
            bool reasonableDeltaY = fabsf( delta.y ) > ( std::max )( fabsf( start.y ), fabsf( end.y ) ) * 0.00001f;

            if( !reasonableDeltaX && !reasonableDeltaY )
                return false;

            // Use a relative bounds to determine whether it's worth testing the delta x
            // values.  When it's really vertical, the previous test will catch it.
            if( reasonableDeltaX ) {
                if( delta.x > 0 ) {
                    maximumOfNearDistances = ( m_minimum.x - start.x ) / delta.x;
                    minimumOfFarDistances = ( m_maximum.x - start.x ) / delta.x;
                } else {
                    maximumOfNearDistances = ( m_maximum.x - start.x ) / delta.x;
                    minimumOfFarDistances = ( m_minimum.x - start.x ) / delta.x;
                }
            }

            // Use a relative bounds to determine whether it's worth testing the delta y
            // values.  When it's really horizontal, the previous test will catch it.
            if( reasonableDeltaY ) {
                if( delta.y > 0 ) {
                    maximumOfNearDistances =
                        ( std::max )( ( m_minimum.y - start.y ) / delta.y, maximumOfNearDistances );
                    minimumOfFarDistances = ( std::min )( ( m_maximum.y - start.y ) / delta.y, minimumOfFarDistances );
                } else {
                    maximumOfNearDistances =
                        ( std::max )( ( m_maximum.y - start.y ) / delta.y, maximumOfNearDistances );
                    minimumOfFarDistances = ( std::min )( ( m_minimum.y - start.y ) / delta.y, minimumOfFarDistances );
                }
            }

            if( maximumOfNearDistances >= 1 || minimumOfFarDistances <= 0 ) {
                return false;
            }

            // if( _isnan( maximumOfNearDistances ) || _isnan( minimumOfFarDistances ) )
            //	throw std::runtime_error( "is_area_intersecting_line_segment: got a NaN" );

            return ( maximumOfNearDistances <= minimumOfFarDistances );
        }
    }

    vector2f get_mapped_from_unit_square( const vector2f& squareCoord ) const {
        return vector2f( m_minimum.x * ( 1 - squareCoord.x ) + m_maximum.x * squareCoord.x,
                         m_minimum.y * ( 1 - squareCoord.y ) + m_maximum.y * squareCoord.y );
    }

    bool is_square( float tolerance = 0.000005f ) const {
        float err = fabs( xsize() - ysize() );
        return err < tolerance * size().get_max_dimension();
    }

    //////////////
    // Operators
    //////////////

    void expand_fractional( value_type fraction ) {
        float expandAmount = .5f * fraction * size().get_max_dimension();
        expand( expandAmount );
    }
};

} // namespace graphics2d
} // namespace frantic
