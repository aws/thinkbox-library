// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace graphics {

// Ensure proper alignment for the data members
// see operator[]
#pragma pack( push, 4 )
class vector3 {
  public:
    typedef boost::int32_t value_type;

  public:
    value_type x, y, z;

    vector3() {
        x = 0;
        y = 0;
        z = 0;
    }

    vector3( value_type X, value_type Y, value_type Z ) {
        x = X;
        y = Y;
        z = Z;
    }

    explicit vector3( const value_type* vec ) {
        x = vec[0];
        y = vec[1];
        z = vec[2];
    }

    explicit vector3( value_type X ) {
        x = X;
        y = X;
        z = X;
    }

    static vector3 xaxis() { return vector3( 1, 0, 0 ); }

    static vector3 yaxis() { return vector3( 0, 1, 0 ); }

    static vector3 zaxis() { return vector3( 0, 0, 1 ); }

    static vector3 from_axis( int axis ) {
        switch( axis ) {
        case 0:
            return vector3( 1, 0, 0 );
        case 1:
            return vector3( 0, 1, 0 );
        case 2:
            return vector3( 0, 0, 1 );
        default:
            return vector3();
        }
    }

    static vector3 from_floor( vector3f vec ) {
        return vector3( (value_type)floorf( vec.x ), (value_type)floorf( vec.y ), (value_type)floorf( vec.z ) );
    }

    static vector3 from_ceil( vector3f vec ) {
        return vector3( (value_type)ceilf( vec.x ), (value_type)ceilf( vec.y ), (value_type)ceilf( vec.z ) );
    }

    static vector3 from_cube_face( int index ) {
        switch( index ) {
        case 0:
            return vector3( 1, 0, 0 );
        case 1:
            return vector3( -1, 0, 0 );
        case 2:
            return vector3( 0, 1, 0 );
        case 3:
            return vector3( 0, -1, 0 );
        case 4:
            return vector3( 0, 0, 1 );
        case 5:
            return vector3( 0, 0, -1 );
        default:
            throw std::runtime_error( "vector3::from_cube_face: Requested invalid cube face index." );
        }
    }

    // builds vector from the minimum of component of the given vector3s
    static vector3 from_min( vector3 a, vector3 b ) {
        return vector3( ( std::min )( a.x, b.x ), ( std::min )( a.y, b.y ), ( std::min )( a.z, b.z ) );
    }

    // builds vector from the maximum of component of the given vector3s
    static vector3 from_max( vector3 a, vector3 b ) {
        return vector3( ( std::max )( a.x, b.x ), ( std::max )( a.y, b.y ), ( std::max )( a.z, b.z ) );
    }

    static vector3 minvalue() { return vector3( ( std::numeric_limits<value_type>::min )() ); }

    static vector3 maxvalue() { return vector3( ( std::numeric_limits<value_type>::max )() ); }

    // A valid neighbor had absolute value 1 and has two 0 coordinates
    bool is_valid_neighbor() const { return ( x * x + y * y + z * z ) == 1; }

    // An extended neighbor includes the all diagonal neighbors and the
    // current position itself.
    bool is_valid_extended_neighbor() const { return x * x <= 1 && y * y <= 1 && z * z <= 1; }

    // Converts this vector into a cube face index, assuming the vector is a valid neighbor vector.  If
    // it's not a valid neighbor vector, it returns -1.
    int neighbor_vector_to_cube_face_nothrow() const {
        if( x != 0 ) {
            if( y == 0 && z == 0 ) {
                if( x == 1 )
                    return 0;
                if( x == -1 )
                    return 1;
            }
        } else {
            if( z == 0 ) {
                if( y == 1 )
                    return 2;
                if( y == -1 )
                    return 3;
            } else if( y == 0 ) {
                if( z == 1 )
                    return 4;
                if( z == -1 )
                    return 5;
            }
        }
        // If none of the above tests succeeded, the coordinate is invalid
        return -1;
    }

    int neighbor_vector_to_cube_face() const {
        int result = neighbor_vector_to_cube_face_nothrow();
        // If none of the above tests succeeded, the coordinate is invalid
        if( result < 0 )
            throw std::runtime_error(
                "vector3.neighbor_vector_to_cube_face: Tried to convert invalid neighbor vector " + str() +
                " to a neighbor index" );
        return result;
    }

    /////////////////////
    // Operators
    /////////////////////

    void set( value_type V ) {
        x = V;
        y = V;
        z = V;
    }

    void set( value_type X, value_type Y, value_type Z ) {
        x = X;
        y = Y;
        z = Z;
    }

    vector3& operator+=( const vector3& a ) {
        x += a.x;
        y += a.y;
        z += a.z;
        return *this;
    }

    vector3& operator-=( const vector3& a ) {
        x -= a.x;
        y -= a.y;
        z -= a.z;
        return *this;
    }

    operator vector3f() { return vector3f( (float)x, (float)y, (float)z ); }

    operator const vector3f() const { return vector3f( (float)x, (float)y, (float)z ); }

    operator vector3fd() { return vector3fd( (double)x, (double)y, (double)z ); }

    operator const vector3fd() const { return vector3fd( (double)x, (double)y, (double)z ); }

    // Uses member alignment ensured by #pragma pack
    const value_type& operator[]( std::size_t index ) const { return ( &x )[index]; }

    value_type& operator[]( std::size_t index ) { return ( &x )[index]; }

    /**
     * Returns this vector3 as a sorted triple.
     * This comes up ad-hoc often enough that I think
     * we can stand to make it a method.
     *
     * @return A vector3 with the same components as this one sorted in non-decreasing order
     */
    vector3 to_sorted() const {
        if( x < y )
            if( y < z ) {
                return vector3( x, y, z );
            } else if( x < z ) {
                return vector3( x, z, y );
            } else {
                return vector3( z, x, y );
            }
        else if( x < z ) {
            return vector3( y, x, z );
        } else if( y < z ) {
            return vector3( y, z, x );
        } else {
            return vector3( z, y, x );
        }
    }

    std::string str() const;
};
#pragma pack( pop )

inline bool operator==( const vector3& a, const vector3& b ) { return a.x == b.x && a.y == b.y && a.z == b.z; }

inline bool operator<( const vector3& a, const vector3& b ) {
    if( a.x == b.x )
        if( a.y == b.y )
            return a.z < b.z;
        else
            return a.y < b.y;
    else
        return a.x < b.x;
}

inline bool operator!=( const vector3& a, const vector3& b ) { return a.x != b.x || a.y != b.y || a.z != b.z; }

inline vector3f operator*( float k, const vector3& a ) { return vector3f( k * a.x, k * a.y, k * a.z ); }

inline vector3f operator*( const vector3& a, float k ) { return vector3f( a.x * k, a.y * k, a.z * k ); }

inline vector3 operator*( vector3::value_type k, const vector3& a ) { return vector3( k * a.x, k * a.y, k * a.z ); }

inline vector3 operator*( const vector3& a, vector3::value_type k ) { return vector3( a.x * k, a.y * k, a.z * k ); }

inline vector3 operator/( const vector3& a, vector3::value_type k ) { return vector3( a.x / k, a.y / k, a.z / k ); }

inline vector3 operator+( const vector3& a, const vector3& b ) { return vector3( a.x + b.x, a.y + b.y, a.z + b.z ); }

inline vector3f operator+( const vector3f& a, const vector3& b ) { return vector3f( a.x + b.x, a.y + b.y, a.z + b.z ); }

inline vector3f operator+( const vector3& a, const vector3f& b ) { return vector3f( a.x + b.x, a.y + b.y, a.z + b.z ); }

inline vector3 operator-( const vector3& a, const vector3& b ) { return vector3( a.x - b.x, a.y - b.y, a.z - b.z ); }

inline vector3 operator-( const vector3& a ) { return vector3( -a.x, -a.y, -a.z ); }

template <class CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& out, const vector3& v ) {
    out << "( vec3 [" << v.x << ", " << v.y << ", " << v.z << "] )";
    return out;
}

inline std::string vector3::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

/**
 *  A hash function for a vector3, suitable for use with unordered_map.
 */
class vector3_hasher {
  public:
    /**
     *  Return a hash value for a voxel coordinate.
     *
     *  Hash from Teschner et al. "Optimized Spatial Hashing for Collision
     * Detection of Deformable Objects".
     *
     * @param voxelCoord the voxel coordinate to generate a hash value for.
     * @return a hash value for the voxel coordinate.
     */
    static size_t hash( const frantic::graphics::vector3& voxelCoord ) {
        const size_t p1 = 73856093;
        const size_t p2 = 19349663;
        const size_t p3 = 83492791;

        return ( voxelCoord.x * p1 ) ^ ( voxelCoord.y * p2 ) ^ ( voxelCoord.z * p3 );
    }

    size_t operator()( const vector3& a ) const { return hash( a ); }
};

// A hash function for a vector3 where the order of the x, y, and z components should be ignored.
// An example is if you are hashing triangles, where the x, y, and z components are indexes into
// the vertex array.  See tetramesh3.hpp for some usage of this class.
class unordered_vector3_hasher {
  public:
    size_t operator()( const vector3& a ) const {
        // Needs a hash that does not depend on vertex order, thus we sort the components.
        vector3 vec = a.to_sorted();
        return vector3_hasher::hash( vec );
    }

    bool operator()( const vector3& a, const vector3& b ) const {
        // Sort the components of both vectors, to make the hashing independent of the component order
        return a.to_sorted() == b.to_sorted();
    }
};

} // namespace graphics
} // namespace frantic
