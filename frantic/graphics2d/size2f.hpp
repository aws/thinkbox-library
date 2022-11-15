// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics2d/size2t.hpp>
#include <frantic/graphics2d/vector2f.hpp>

namespace frantic {
namespace graphics2d {

class size2f : public size2t<float, size2f> {
  public:
    size2f()
        : size2t<float, size2f>( 0.0f, 0.0f ) {}
    size2f( float w )
        : size2t<float, size2f>( w, w ) {}
    size2f( float w, float h )
        : size2t<float, size2f>( w, h ) {}
    template <typename T, typename SizeType>
    explicit size2f( const size2t<T, SizeType>& s )
        : size2t<float, size2f>( (float)s.xsize, (float)s.ysize ) {}

    typedef vector2f vector_type;

    bool contains( const vector_type& coord ) const {
        return 0 <= coord.x && coord.x < xsize && 0 <= coord.y && coord.y < ysize;
    }

    static size2f from_bounds( const vector2f& minimum, const vector2f& maximum ) {
        size2f result( maximum.x - minimum.x, maximum.y - minimum.y );
        if( result.xsize < 0 || result.ysize < 0 ) {
            result.xsize = 0;
            result.ysize = 0;
        }
        return result;
    }
};

} // namespace graphics2d
} // namespace frantic
