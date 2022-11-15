// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace graphics {

template <class FloatType>
class size3t {
  public:
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;

  private:
    float_type m_xsize, m_ysize, m_zsize;

  public:
    size3t() {
        m_xsize = 0;
        m_ysize = 0;
        m_zsize = 0;
    }

    size3t( float_type w, float_type h, float_type d ) {
        m_xsize = w;
        m_ysize = h;
        m_zsize = d;
    }

    static size3t from_cube_side_length( float_type sideLength ) {
        return size3t( sideLength, sideLength, sideLength );
    }

    static size3t empty() { return size3t( 0, 0, 0 ); }

    static size3t from_bounds( const vector3t<float_type>& minimum, const vector3t<float_type>& maximum ) {
        size3t result( maximum.x - minimum.x, maximum.y - minimum.y, maximum.z - minimum.z );
        if( result.m_xsize < 0 || result.m_ysize < 0 || result.m_zsize < 0 ) {
            result.m_xsize = 0;
            result.m_ysize = 0;
            result.m_zsize = 0;
        }
        return result;
    }

    float_type xsize() const { return m_xsize; }

    float_type ysize() const { return m_ysize; }

    float_type zsize() const { return m_zsize; }

    float_type get_max_dimension() const { return ( std::max )( m_xsize, ( std::max )( m_ysize, m_zsize ) ); }

    float_type get_min_dimension() const { return ( std::min )( m_xsize, ( std::min )( m_ysize, m_zsize ) ); }

    std::string str() const;

    //////////////
    // Operators
    //////////////

    friend size3t<FloatType> operator*( const typename size3t<FloatType>::vector3f_type& a,
                                        const size3t<FloatType>& b ) {
        return size3t<FloatType>( a.x * b.xsize(), a.y * b.ysize(), a.z * b.zsize() );
    }

    friend size3t<FloatType> operator*( const size3t<FloatType>& a,
                                        const typename size3t<FloatType>::vector3f_type& b ) {
        return b * a;
    }
};

//////////////
// Operators
//////////////

template <class FloatType>
inline size3t<FloatType> operator*( FloatType a, const size3t<FloatType>& b ) {
    return size3t<FloatType>( a * b.xsize(), a * b.ysize(), a * b.zsize() );
}

template <class FloatType>
inline size3t<FloatType> operator*( const size3t<FloatType>& a, FloatType b ) {
    return b * a;
}

template <class FloatType>
inline std::ostream& operator<<( std::ostream& out, const size3t<FloatType>& v ) {
    out << "(size3f " << v.xsize() << ", " << v.ysize() << ", " << v.zsize() << " )";
    return out;
}

template <class FloatType>
inline std::string size3t<FloatType>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

template <class FloatType>
inline typename size3t<FloatType>::vector3f_type operator+( const typename size3t<FloatType>::vector3f_type& a,
                                                            const size3t<FloatType>& b ) {
    return typename size3t<FloatType>::vector3f_type( a.x + b.xsize(), a.y + b.ysize(), a.z + b.zsize() );
}

typedef size3t<float> size3f;
typedef size3t<double> size3fd;

} // namespace graphics
} // namespace frantic
