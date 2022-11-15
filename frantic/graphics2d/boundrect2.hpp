// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics2d/boundrect2t.hpp>
#include <frantic/graphics2d/size2.hpp>
#include <frantic/graphics2d/vector2.hpp>

namespace frantic {
namespace graphics2d {

class boundrect2 : public boundrect2t<vector2, size2, boundrect2> {
  public:
    //////////////////
    // Constructors //
    //////////////////

    boundrect2()
        : boundrect2t<vector2, size2, boundrect2>( vector2( 1, 1 ), vector2( 0, 0 ) ) {}
    explicit boundrect2( const vector2& position )
        : boundrect2t<vector2, size2, boundrect2>( position ) {}
    boundrect2( int xmin, int xmax, int ymin, int ymax )
        : boundrect2t<vector2, size2, boundrect2>( xmin, xmax, ymin, ymax ) {}
    boundrect2( const vector2& min, const vector2& max )
        : boundrect2t<vector2, size2, boundrect2>( min, max ) {}
    boundrect2( const vector2& min, const size2& size )
        : boundrect2t<vector2, size2, boundrect2>(
              min, vector2( min.get_x() + size.get_xsize() - 1, min.get_y() + size.get_ysize() - 1 ) ) {}

    value_type xsize() const { return m_maximum.x - m_minimum.x + 1; }
    value_type ysize() const { return m_maximum.y - m_minimum.y + 1; }

    void set_to_point( const vector2& position ) {
        m_minimum = position;
        m_maximum = position;
    }
    void set_empty() {
        m_minimum = vector2( 1, 1 );
        m_maximum = vector2( 0, 0 );
    }

    bool is_square() const { return xsize() == ysize(); }
};

} // namespace graphics2d
} // namespace frantic
