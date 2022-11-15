// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <vector>

namespace frantic {
namespace volumetrics {

class voxel_coord_system {
    float m_voxelLength;                       // length of a voxel
    frantic::graphics::vector3f m_worldOrigin; // origin of voxel space in world
  public:
    voxel_coord_system() {
        // By default, just do an identity transform
        m_worldOrigin = frantic::graphics::vector3f( 0 );
        m_voxelLength = 1.0f;
    }

    voxel_coord_system( const frantic::graphics::vector3f& worldOrigin, float voxelLength ) {
        m_voxelLength = voxelLength;
        m_worldOrigin = worldOrigin;
    }

    frantic::graphics::vector3f world_origin() const { return m_worldOrigin; }

    float voxel_length() const { return m_voxelLength; }

    void set_voxel_length( float length ) { m_voxelLength = length; }

    void set_world_origin( const frantic::graphics::vector3f& origin ) { m_worldOrigin = origin; }

    void set( const frantic::graphics::vector3f& worldOrigin, float voxelLength ) {
        m_voxelLength = voxelLength;
        m_worldOrigin = worldOrigin;
    }

    void swap( voxel_coord_system& rhs ) {
        using std::swap;
        swap( m_voxelLength, rhs.m_voxelLength );
        swap( m_worldOrigin, rhs.m_worldOrigin );
    }

    // TODO: This probably needs to be adjusted, based on the precise meaning of
    //       a voxel.  For instance a voxel bounds whose min and max are the same would
    //       become a world bounds whose min and max are one voxel length apart.
    frantic::graphics::boundbox3f get_world_bounds( const frantic::graphics::boundbox3& voxelBounds ) const {
        return frantic::graphics::boundbox3f(
            get_world_coord( voxelBounds.minimum() ),
            get_world_coord( voxelBounds.maximum() + frantic::graphics::vector3( 1 ) ) );
    }

    /**
     * This returns the world-space bounding box containing all the voxel centers of the voxels
     * in the provided voxelBounds.
     */
    frantic::graphics::boundbox3f
    get_world_bounds_of_voxel_centers( const frantic::graphics::boundbox3& voxelBounds ) const {
        return frantic::graphics::boundbox3f(
            get_world_coord( voxelBounds.minimum() + frantic::graphics::vector3f( 0.5f ) ),
            get_world_coord( voxelBounds.maximum() + frantic::graphics::vector3f( 0.5f ) ) );
    }

    /**
     * This returns the smallest voxel bounding box which contains the given world bounds.
     */
    frantic::graphics::boundbox3 get_voxel_bounds( const frantic::graphics::boundbox3f& worldBounds ) const {
        if( worldBounds.is_empty() ) {
            return frantic::graphics::boundbox3::from_empty();
        } else {
            frantic::graphics::vector3f minimum = get_voxel_coord( worldBounds.minimum() ),
                                        maximum = get_voxel_coord( worldBounds.maximum() );
            return frantic::graphics::boundbox3(
                frantic::graphics::vector3( (int)floorf( minimum.x ), (int)floorf( minimum.y ),
                                            (int)floorf( minimum.z ) ),
                frantic::graphics::vector3( (int)ceilf( maximum.x ), (int)ceilf( maximum.y ),
                                            (int)ceilf( maximum.z ) ) );
        }
    }

    /**
     * Converts the provided world-space coordinate into voxel space.
     */
    frantic::graphics::vector3f get_voxel_coord( const frantic::graphics::vector3f& worldLocation ) const {
        return ( worldLocation - m_worldOrigin ) / m_voxelLength;
    }

    /**
     * Converts the provided world-space X coordinate into voxel space.
     */
    float get_voxel_x_coord( float worldX ) const { return ( worldX - m_worldOrigin.x ) / m_voxelLength; }

    /**
     * Converts the provided world-space Y coordinate into voxel space.
     */
    float get_voxel_y_coord( float worldY ) const { return ( worldY - m_worldOrigin.y ) / m_voxelLength; }

    /**
     * Converts the provided world-space Z coordinate into voxel space.
     */
    float get_voxel_z_coord( float worldZ ) const { return ( worldZ - m_worldOrigin.z ) / m_voxelLength; }

    /**
     * Converts the provided voxel-space coordinate into world space.
     */
    frantic::graphics::vector3f get_world_coord( const frantic::graphics::vector3f& voxelCoord ) const {
        return voxelCoord * m_voxelLength + m_worldOrigin;
    }

    /**
     * Converts the provided voxel-space coordinate into world space.
     */
    frantic::graphics::vector3f get_world_coord( const frantic::graphics::vector3& voxelCoord ) const {
        return voxelCoord * m_voxelLength + m_worldOrigin;
    }

    /**
     * Converts the provided voxel-space x coordinate into world space.
     */
    float get_world_x_coord( float voxelX ) const { return voxelX * m_voxelLength + m_worldOrigin.x; }

    /**
     * Converts the provided voxel-space y coordinate into world space.
     */
    float get_world_y_coord( float voxelY ) const { return voxelY * m_voxelLength + m_worldOrigin.y; }

    /**
     * Converts the provided voxel-space z coordinate into world space.
     */
    float get_world_z_coord( float voxelZ ) const { return voxelZ * m_voxelLength + m_worldOrigin.z; }

    frantic::graphics::vector3f get_world_voxel_center( const frantic::graphics::vector3& voxelIndex ) const {
        return ( voxelIndex + frantic::graphics::vector3f( 0.5f ) ) * m_voxelLength + m_worldOrigin;
    }

    float get_world_x_voxel_center( boost::int32_t x ) const { return ( x + 0.5f ) * m_voxelLength + m_worldOrigin.x; }

    float get_world_y_voxel_center( boost::int32_t y ) const { return ( y + 0.5f ) * m_voxelLength + m_worldOrigin.y; }

    float get_world_z_voxel_center( boost::int32_t z ) const { return ( z + 0.5f ) * m_voxelLength + m_worldOrigin.z; }

    /**
     * Returns the face center in World Coordinates. The face indices follow
     * (x,x+1,y,y+1,z,z+1)
     *
     *@param voxelCoord the coord of the voxel
     *@param faceIndex the index of the face
     */
    frantic::graphics::vector3f get_world_face_center( const frantic::graphics::vector3& voxelCoord,
                                                       int faceIndex ) const {

        frantic::graphics::vector3f v = get_world_voxel_center( voxelCoord );

        switch( faceIndex ) {
        // X
        case 0:
            v.x -= 0.5f * m_voxelLength;
            break;
        case 1:
            v.x += 0.5f * m_voxelLength;
            break;

        // Y
        case 2:
            v.y -= 0.5f * m_voxelLength;
            break;
        case 3:
            v.y += 0.5f * m_voxelLength;
            break;

        // 4 z neg
        case 4:
            v.z -= 0.5f * m_voxelLength;
            break;
        case 5:
            v.z += 0.5f * m_voxelLength;
            break;
        default:
            throw std::runtime_error( "get_world_face_center() -  Face Index ( " +
                                      boost::lexical_cast<std::string>( faceIndex ) +
                                      " ) is not valid. It must be in [0,5]" );
        }
        return v;
    }

    /**
     * Returns the face center in Voxel Coordinates. The face indices follow
     * (x,x+1,y,y+1,z,z+1)
     *
     *@param voxelCoord the coord of the voxel
     *@param faceIndex the index of the face
     */
    frantic::graphics::vector3f get_voxel_face_center( const frantic::graphics::vector3& voxelCoord,
                                                       int faceIndex ) const {
        frantic::graphics::vector3f v( voxelCoord );
        switch( faceIndex ) {
        // X
        case 0:
            v.y += 0.5f;
            v.z += 0.5f;
            break;
        case 1:
            v.y += 0.5f;
            v.z += 0.5f;
            v.x += 1.0f;
            break;

        // Y
        case 2:
            v.x += 0.5f;
            v.z += 0.5f;
            break;
        case 3:
            v.x += 0.5f;
            v.y += 1.0f;
            v.z += 0.5f;
            break;

        // 4 z neg
        case 4:
            v.x += 0.5f;
            v.y += 0.5f;
            break;
        case 5:
            v.x += 0.5f;
            v.y += 0.5f;
            v.z += 1.0f;
            break;
        default:
            throw std::runtime_error( "get_voxel_face_center() -  Face Index ( " +
                                      boost::lexical_cast<std::string>( faceIndex ) +
                                      " ) is not valid. It must be in [0,5]" );
        }
        return v;
    }

    std::string str() const {
        return std::string() + "WorldOrigin: " + m_worldOrigin.str() +
               ", length: " + boost::lexical_cast<std::string>( m_voxelLength );
    }

    /**
     * Returns a transform matrix which converts a voxel coordinate to the world coordinate of the voxel center (May be
     * buggy - to be reviewed)
     */
    graphics::transform4f voxel_center_to_world_transform() const {
        float half_length = m_voxelLength / 2.0f;

        return graphics::transform4f( m_voxelLength, 0, 0, 0, 0, m_voxelLength, 0, 0, 0, 0, m_voxelLength, 0,
                                      m_worldOrigin.x + half_length, m_worldOrigin.y + half_length,
                                      m_worldOrigin.z + half_length, 1 );
    }

    /**
     * Returns a transform matrix which converts world coordinate of a voxel center to voxel coordinates. (May be buggy
     * - to be reviewed)
     */
    graphics::transform4f world_voxel_center_to_voxel_transform() const {
        float inverseVoxelLength = 1 / m_voxelLength;
        return graphics::transform4f( inverseVoxelLength, 0, 0, 0, 0, inverseVoxelLength, 0, 0, 0, 0,
                                      inverseVoxelLength, 0, -m_worldOrigin.x * inverseVoxelLength - 0.5f,
                                      -m_worldOrigin.y * inverseVoxelLength - 0.5f,
                                      -m_worldOrigin.z * inverseVoxelLength - 0.5f, 1 );
    }

    /**
     * Returns a transform matrix which converts voxel coordinates to world coordinates.
     */
    graphics::transform4f voxel_to_world_transform() const {
        return graphics::transform4f( m_voxelLength, 0, 0, 0, 0, m_voxelLength, 0, 0, 0, 0, m_voxelLength, 0,
                                      m_worldOrigin.x, m_worldOrigin.y, m_worldOrigin.z, 1 );
    }

    /**
     * Returns a transform matrix which converts world coordinates to voxel coordinates.
     */
    graphics::transform4f world_to_voxel_transform() const {
        float inverseVoxelLength = 1 / m_voxelLength;
        return graphics::transform4f( inverseVoxelLength, 0, 0, 0, 0, inverseVoxelLength, 0, 0, 0, 0,
                                      inverseVoxelLength, 0, -m_worldOrigin.x * inverseVoxelLength,
                                      -m_worldOrigin.y * inverseVoxelLength, -m_worldOrigin.z * inverseVoxelLength, 1 );
    }

    /* TODO
    possible implementation interface which will be more natural template type interfaces:
    <type> voxel_to_world( <type> )
    <type> world_to_voxel( <type> )
    //*/

    bool operator==( const voxel_coord_system& rhs ) const {
        return m_voxelLength == rhs.m_voxelLength && m_worldOrigin == rhs.m_worldOrigin;
    }

    bool operator!=( const voxel_coord_system& rhs ) const {
        return m_voxelLength != rhs.m_voxelLength || m_worldOrigin != rhs.m_worldOrigin;
    }

    bool equals( const voxel_coord_system& rhs, float tolerance = 1e-5 ) const;
};

inline std::ostream& operator<<( std::ostream& out, const voxel_coord_system& vcs ) {
    out << "voxel coord system( " << vcs.world_origin() << ", " << vcs.voxel_length() << " )";
    return out;
}

inline bool voxel_coord_system::equals( const voxel_coord_system& rhs, float tolerance ) const {
    bool result =
        fabsf( m_voxelLength - rhs.m_voxelLength ) < ( std::max )( m_voxelLength, rhs.m_voxelLength ) * tolerance &&
        m_worldOrigin.is_equal( rhs.m_worldOrigin, tolerance );
    /*
      if( !result ) {
        std::ofstream fout("c:\\debugvcs.txt");
        fout << *this << std::endl;
        fout << rhs << std::endl;
        fout << tolerance << std::endl;
        fout << fabsf(m_voxelLength - rhs.m_voxelLength) << " < " << (std::max)(m_voxelLength, rhs.m_voxelLength) *
      tolerance << std::endl; fout << frantic::graphics::vector3f::distance_squared(m_worldOrigin, rhs.m_worldOrigin) <<
      " < " << tolerance * tolerance * (std::max)(m_worldOrigin.get_magnitude_squared(),
      rhs.m_worldOrigin.get_magnitude_squared()) << std::endl;
      }
    */
    return result;
}

} // namespace volumetrics
} // namespace frantic
