// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stdexcept>

#include <frantic/graphics/vector3.hpp>

namespace frantic {
namespace graphics {

// Ensure proper alignment for the data members
// see operator[]
#pragma pack( push, 4 )
class size3 {
    // Formerly width, height, depth.
    int m_xsize, m_ysize, m_zsize;

  public:
    explicit size3( int w ) { m_xsize = m_ysize = m_zsize = w; }

    size3( int xsize, int ysize, int zsize ) {
        m_xsize = xsize;
        m_ysize = ysize;
        m_zsize = zsize;
    }

    size3() {
        m_xsize = 0;
        m_ysize = 0;
        m_zsize = 0;
    }

    void set( int size ) {
        m_xsize = size;
        m_ysize = size;
        m_zsize = size;
    }

    void set( int xsize, int ysize, int zsize ) {
        m_xsize = xsize;
        m_ysize = ysize;
        m_zsize = zsize;
    }

    // The box represented by the minimum and maximum values includes both minimum and maximum.
    // This is the same convention that is used in boundbox3.
    // NOTE: This is different from the floating point size, where the size intervals are treated as half-open
    // intervals.
    static size3 from_bounds( const vector3& minimum, const vector3& maximum ) {
        size3 result( maximum.x - minimum.x + 1, maximum.y - minimum.y + 1, maximum.z - minimum.z + 1 );
        if( result.xsize() < 0 || result.ysize() < 0 || result.zsize() < 0 ) {
            result.m_xsize = 0;
            result.m_ysize = 0;
            result.m_zsize = 0;
        }
        return result;
    }

    static size3 from_cube_side_length( int sideLength ) { return size3( sideLength, sideLength, sideLength ); }

    ///////////////////
    // Accessors
    ///////////////////

    int& xsize() { return m_xsize; }

    int xsize() const { return m_xsize; }

    int& ysize() { return m_ysize; }

    int ysize() const { return m_ysize; }

    int& zsize() { return m_zsize; }

    int zsize() const { return m_zsize; }

    bool is_cube() const { return m_xsize == m_ysize && m_ysize == m_zsize; }

    int get_cube_length() const {
        if( !is_cube() ) {
            throw std::runtime_error( "size3::get_cube_length: not a cube" );
        }
        return m_xsize;
    }

    vector3 clamp( const vector3& vec ) const {
        return vector3( (std::max<boost::int32_t>)( 0, (std::min<boost::int32_t>)( m_xsize - 1, vec.x ) ),
                        (std::max<boost::int32_t>)( 0, (std::min<boost::int32_t>)( m_ysize - 1, vec.y ) ),
                        (std::max<boost::int32_t>)( 0, (std::min<boost::int32_t>)( m_zsize - 1, vec.z ) ) );
    }

    bool contains( const vector3& coord ) const {
        // Do all 6 comparisons in 3 by using unsigned numbers.  If coord.x < 0, then its unsigned representation will
        // be bigger than m_xsize.  Note that this only works for integers.
        return ( (unsigned)coord.x < (unsigned)m_xsize ) && ( (unsigned)coord.y < (unsigned)m_ysize ) &&
               ( (unsigned)coord.z < (unsigned)m_zsize );
    }

    bool is_coord_valid( vector3 coord ) const { return contains( coord ); }

    /**
     * Gets the index corresponding to the vector.  If the vector is out of range, it throws an exception.
     *
     * @param   vec  The input vector coordinate.
     */
    int get_index( vector3 vec ) const {
        if( !contains( vec ) )
            throw std::runtime_error( "size3.get_index: Tried to get the index of an out of bounds vector: " +
                                      vec.str() + ", " + this->str() );
        return vec.x + ( vec.y + vec.z * m_ysize ) * m_xsize;
    }

    /**
     * Gets the index corresponding to the vector.  If the vector is out of range, it returns -1.
     *
     * @param   vec  The input vector coordinate.
     */
    int get_index_nothrow( vector3 vec ) const {
        if( !contains( vec ) )
            return -1;
        return vec.x + ( vec.y + vec.z * m_ysize ) * m_xsize;
    }

    /**
     * Gets the index corresponding to the vector.  If the vector is out of range, the behavior is unspecified, i.e. may
     * be valid or invalid. This is intended to be faster, because it doesn't do a bounds check at all.
     *
     * @param   vec  The input vector coordinate.
     */
    int get_index_nocheck( vector3 vec ) const { return vec.x + ( vec.y + vec.z * m_ysize ) * m_xsize; }

    // Gets the vector corresponding to the index
    vector3 get_coord( int index ) const {
        return vector3( index % m_xsize, ( index / m_xsize ) % m_ysize, ( index / ( m_xsize * m_ysize ) ) % m_zsize );
    }

    // Assuming we have a grid of blocks the size of this size3, this function returns
    // the coordinate in that higher level grid
    vector3 get_quotient_coord( vector3 vec ) const {
        // Can't just use /, because then negative numbers go wrong
        vector3 result;
        if( vec.x >= 0 )
            result.x = vec.x / m_xsize;
        else
            result.x = -1 - ( ( -1 - vec.x ) / m_xsize );
        if( vec.y >= 0 )
            result.y = vec.y / m_ysize;
        else
            result.y = -1 - ( ( -1 - vec.y ) / m_ysize );
        if( vec.z >= 0 )
            result.z = vec.z / m_zsize;
        else
            result.z = -1 - ( ( -1 - vec.z ) / m_zsize );

        return result;
    }

    // Assuming we have a grid of blocks the size of this size3, this function returns
    // the coordinate in the voxel that contains the vector3 vec
    vector3 get_kernel_coord( vector3 vec ) const {
        vector3 quotientCoord = get_quotient_coord( vec );
        return vector3( vec.x - m_xsize * quotientCoord.x, vec.y - m_ysize * quotientCoord.y,
                        vec.z - m_zsize * quotientCoord.z );
    }

    int volume() const { return m_xsize * m_ysize * m_zsize; }

    int get_max_dimension() const { return ( std::max )( m_xsize, ( std::max )( m_ysize, m_zsize ) ); }

    int get_min_dimension() const { return ( std::min )( m_xsize, ( std::min )( m_ysize, m_zsize ) ); }

    bool operator==( const size3& rhs ) const {
        return m_xsize == rhs.m_xsize && m_ysize == rhs.m_ysize && m_zsize == rhs.m_zsize;
    }

    bool operator!=( const size3& rhs ) const {
        return m_xsize != rhs.m_xsize || m_ysize != rhs.m_ysize || m_zsize != rhs.m_zsize;
    }

    // Uses member alignment ensured by #pragma pack
    int operator[]( int index ) const { return ( &m_xsize )[index]; }

    std::string str() const;
};
#pragma pack( pop )

inline vector3 operator+( const vector3& a, const size3& b ) {
    return vector3( a.x + b.xsize(), a.y + b.ysize(), a.z + b.zsize() );
}

inline vector3 operator*( const vector3& a, const size3& b ) {
    return vector3( a.x * b.xsize(), a.y * b.ysize(), a.z * b.zsize() );
}

inline size3 operator*( int a, const size3& b ) { return size3( a * b.xsize(), a * b.ysize(), a * b.zsize() ); }

inline vector3 operator*( const size3& a, const vector3& b ) {
    return vector3( a.xsize() * b.x, a.ysize() * b.y, a.zsize() * b.z );
}

template <class CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const size3& v ) {
    out << "(size3 " << v.xsize() << ", " << v.ysize() << ", " << v.zsize() << " )";
    return out;
}

inline std::string size3::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

} // namespace graphics
} // namespace frantic
