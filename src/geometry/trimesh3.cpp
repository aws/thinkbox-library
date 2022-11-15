// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/trimesh3.hpp>

#include <frantic/geometry/mesh_interface_utils.hpp>

#include <frantic/geometry/trimesh3_degeneracy_removal.hpp>

#include <frantic/geometry/polygon_utils.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>

#include <boost/foreach.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::channels;
using namespace frantic::graphics2d;

namespace frantic {
namespace geometry {

trimesh3::trimesh3( const trimesh3& rhs )
    : m_vertices( rhs.m_vertices )
    , m_faces( rhs.m_faces )
    , m_namedVertexChannels( rhs.m_namedVertexChannels )
    , m_namedFaceChannels( rhs.m_namedFaceChannels )
    , m_facePlanes( rhs.m_facePlanes )
    , m_faceBoundBoxes( rhs.m_faceBoundBoxes )
    , m_barycentric0Axis( rhs.m_barycentric0Axis )
    , m_barycentric1Axis( rhs.m_barycentric1Axis )
    , m_barycentricInverseDeterminant( rhs.m_barycentricInverseDeterminant ) {}

void trimesh3::swap( trimesh3& rhs ) {
    m_vertices.swap( rhs.m_vertices );
    m_faces.swap( rhs.m_faces );
    m_facePlanes.swap( rhs.m_facePlanes );
    m_faceBoundBoxes.swap( rhs.m_faceBoundBoxes );
    m_namedVertexChannels.swap( rhs.m_namedVertexChannels );
    m_namedFaceChannels.swap( rhs.m_namedFaceChannels );
    m_barycentric0Axis.swap( rhs.m_barycentric0Axis );
    m_barycentric1Axis.swap( rhs.m_barycentric1Axis );
    m_barycentricInverseDeterminant.swap( rhs.m_barycentricInverseDeterminant );
}

void trimesh3::dump( std::ostream& out ) const {
    out << "===================\n";
    out << "Dumping trimesh3\n";
    out << "\n";
    out << "Vertex count: " << m_vertices.size() << "\n";
    out << "Face count: " << m_faces.size() << "\n";
    out << "Additional vertex channel count: " << m_namedVertexChannels.size() << "\n";
    out << "Vertex channel names:\n";
    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.begin();
         i != m_namedVertexChannels.end(); ++i )
        out << "  " << frantic::strings::to_string( i->first ) << "\n";
    out << "Additional face channel count: " << m_namedFaceChannels.size() << "\n";
    out << "Face channel names:\n";
    for( std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.begin();
         i != m_namedFaceChannels.end(); ++i )
        out << "  " << frantic::strings::to_string( i->first ) << "\n";

    out << "\n";
    out << "Vertices:\n";
    for( unsigned i = 0; i < m_vertices.size(); ++i )
        out << m_vertices[i] << "\n";

    out << "\n";
    out << "Faces:\n";
    for( unsigned i = 0; i < m_faces.size(); ++i )
        out << m_faces[i] << "\n";

    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.begin();
         i != m_namedVertexChannels.end(); ++i ) {
        out << "\n\n";
        out << "Named Vertex Channel \"" << frantic::strings::to_string( i->first ) << "\":\n";
        out << " Has custom faces: " << ( i->second.has_custom_faces() ? "true" : "false" ) << "\n";
        out << " Vertex count: " << i->second.size() << "\n";
        if( i->second.has_custom_faces() )
            out << " Face count: " << i->second.face_count() << "\n";
        out << "Channel vertices:\n";
        const_trimesh3_vertex_channel_general_accessor channel = i->second.get_general_accessor( &m_faces );

        channel.dump( out );

        for( unsigned v = 0; v < channel.size(); ++v ) {
            channel.print( out, v );
            out << "\n";
        }
        if( i->second.has_custom_faces() ) {
            out << "\n";
            out << "Channel faces:\n";
            for( unsigned f = 0; f < channel.face_count(); ++f ) {
                out << channel.face( f ) << "\n";
            }
        }
    }

    for( std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.begin();
         i != m_namedFaceChannels.end(); ++i ) {
        out << "\n\n";
        out << "Named Face Channel \"" << frantic::strings::to_string( i->first ) << "\":\n";
        out << " Face count: " << i->second.size() << "\n";
        out << "Channel data:\n";
        const_trimesh3_face_channel_general_accessor channel = i->second.get_general_accessor();

        channel.dump( out );

        for( unsigned f = 0; f < channel.size(); ++f ) {
            channel.print( out, f );
            out << "\n";
        }
    }
    out << "===================\n";
}

void trimesh3::dump_channel_data( std::ostream& out ) const {
    out << "===================\n";
    out << "Dumping trimesh3 channel data\n";
    out << "\n";
    out << "Vertex count: " << m_vertices.size() << "\n";
    out << "Face count: " << m_faces.size() << "\n";
    out << "Additional vertex channel count: " << m_namedVertexChannels.size() << "\n";
    out << "Vertex channel names:\n";
    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.begin();
         i != m_namedVertexChannels.end(); ++i )
        out << "  " << frantic::strings::to_string( i->first ) << "\n";
    out << "Additional face channel count: " << m_namedFaceChannels.size() << "\n";
    out << "Face channel names:\n";
    for( std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.begin();
         i != m_namedFaceChannels.end(); ++i )
        out << "  " << frantic::strings::to_string( i->first ) << "\n";

    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.begin();
         i != m_namedVertexChannels.end(); ++i ) {
        out << "\n\n";
        out << "Named Vertex Channel \"" << frantic::strings::to_string( i->first ) << "\":\n";
        out << " Has custom faces: " << ( i->second.has_custom_faces() ? "true" : "false" ) << "\n";
        out << " Vertex count: " << i->second.size() << "\n";
        if( i->second.has_custom_faces() )
            out << " Face count: " << i->second.face_count() << "\n";
        const_trimesh3_vertex_channel_general_accessor channel = i->second.get_general_accessor( &m_faces );
        channel.dump( out );
    }

    for( std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.begin();
         i != m_namedFaceChannels.end(); ++i ) {
        out << "\n\n";
        out << "Named Face Channel \"" << frantic::strings::to_string( i->first ) << "\":\n";
        out << " Face count: " << i->second.size() << "\n";
        out << "Channel data:\n";
        const_trimesh3_face_channel_general_accessor channel = i->second.get_general_accessor();
        channel.dump( out );
    }
    out << "===================\n";
}

void trimesh3::set_existing_exploded_box_verts( const boundbox3f& box ) {
    if( m_vertices.size() != 24 )
        throw runtime_error(
            "trimesh3::set_existing_exploded_box_verts() - This function can only be called on a trimesh3 "
            "with 24 vertices, this mesh has " +
            lexical_cast<string>( m_vertices.size() ) + " vertices." );

    // xneg
    m_vertices[0] = box.get_corner( 0 );
    m_vertices[1] = box.get_corner( 2 );
    m_vertices[2] = box.get_corner( 4 );
    m_vertices[3] = box.get_corner( 6 );
    // xpos
    m_vertices[4] = box.get_corner( 1 );
    m_vertices[5] = box.get_corner( 3 );
    m_vertices[6] = box.get_corner( 5 );
    m_vertices[7] = box.get_corner( 7 );
    // yneg
    m_vertices[8] = box.get_corner( 0 );
    m_vertices[9] = box.get_corner( 1 );
    m_vertices[10] = box.get_corner( 4 );
    m_vertices[11] = box.get_corner( 5 );
    // ypos
    m_vertices[12] = box.get_corner( 2 );
    m_vertices[13] = box.get_corner( 3 );
    m_vertices[14] = box.get_corner( 6 );
    m_vertices[15] = box.get_corner( 7 );

    // zneg
    m_vertices[16] = box.get_corner( 0 );
    m_vertices[17] = box.get_corner( 1 );
    m_vertices[18] = box.get_corner( 2 );
    m_vertices[19] = box.get_corner( 3 );
    // zpos
    m_vertices[20] = box.get_corner( 4 );
    m_vertices[21] = box.get_corner( 5 );
    m_vertices[22] = box.get_corner( 6 );
    m_vertices[23] = box.get_corner( 7 );
}

void trimesh3::set_to_exploded_box( const boundbox3f& box ) {
    clear();

    m_vertices.resize( 24 );

    // xneg
    m_vertices[0] = box.get_corner( 0 );
    m_vertices[1] = box.get_corner( 2 );
    m_vertices[2] = box.get_corner( 4 );
    m_vertices[3] = box.get_corner( 6 );
    // xpos
    m_vertices[4] = box.get_corner( 1 );
    m_vertices[5] = box.get_corner( 3 );
    m_vertices[6] = box.get_corner( 5 );
    m_vertices[7] = box.get_corner( 7 );
    // yneg
    m_vertices[8] = box.get_corner( 0 );
    m_vertices[9] = box.get_corner( 1 );
    m_vertices[10] = box.get_corner( 4 );
    m_vertices[11] = box.get_corner( 5 );
    // ypos
    m_vertices[12] = box.get_corner( 2 );
    m_vertices[13] = box.get_corner( 3 );
    m_vertices[14] = box.get_corner( 6 );
    m_vertices[15] = box.get_corner( 7 );

    // zneg
    m_vertices[16] = box.get_corner( 0 );
    m_vertices[17] = box.get_corner( 1 );
    m_vertices[18] = box.get_corner( 2 );
    m_vertices[19] = box.get_corner( 3 );
    // ypos
    m_vertices[20] = box.get_corner( 4 );
    m_vertices[21] = box.get_corner( 5 );
    m_vertices[22] = box.get_corner( 6 );
    m_vertices[23] = box.get_corner( 7 );

    m_faces.resize( 12 );

    // xneg
    m_faces[0] = vector3( 1, 0, 2 );
    m_faces[1] = vector3( 2, 3, 1 );
    // xpos
    m_faces[2] = vector3( 4, 5, 7 );
    m_faces[3] = vector3( 7, 6, 4 );
    // yneg
    m_faces[4] = vector3( 8, 9, 11 );
    m_faces[5] = vector3( 11, 10, 8 );
    // ypos
    m_faces[6] = vector3( 13, 12, 14 );
    m_faces[7] = vector3( 14, 15, 13 );
    // zneg
    m_faces[8] = vector3( 18, 19, 17 );
    m_faces[9] = vector3( 17, 16, 18 );
    // zpos
    m_faces[10] = vector3( 20, 21, 23 );
    m_faces[11] = vector3( 23, 22, 20 );
}

void trimesh3::set_to_box( const boundbox3f& box ) {
    clear();

    m_vertices.resize( 8 );
    for( int i = 0; i < 8; ++i )
        m_vertices[i] = box.get_corner( i );

    m_faces.resize( 12 );
    m_faces[0] = vector3( 1, 0, 2 );
    m_faces[1] = vector3( 2, 3, 1 );
    m_faces[2] = vector3( 6, 4, 5 );
    m_faces[3] = vector3( 5, 7, 6 );
    m_faces[4] = vector3( 5, 4, 0 );
    m_faces[5] = vector3( 0, 1, 5 );
    m_faces[6] = vector3( 2, 6, 7 );
    m_faces[7] = vector3( 7, 3, 2 );
    m_faces[8] = vector3( 0, 4, 6 );
    m_faces[9] = vector3( 6, 2, 0 );
    m_faces[10] = vector3( 7, 5, 1 );
    m_faces[11] = vector3( 1, 3, 7 );
}

void trimesh3::set_existing_box_verts( const boundbox3f& box ) {
    if( m_vertices.size() != 8 )
        throw runtime_error(
            "trimesh3::set_existing_box_verts() - This function can only be called on a trimesh3 with 8 "
            "vertices, this mesh has " +
            lexical_cast<string>( m_vertices.size() ) + " vertices." );

    for( int i = 0; i < 8; ++i )
        m_vertices[i] = box.get_corner( i );
}

void trimesh3::set_to_tetrahedron( const vector3f& centre, float radius ) {
    clear();

    const float f = radius / sqrtf( 3 );

    m_vertices.resize( 4 );
    m_vertices[0] = centre + vector3f( f, f, f );
    m_vertices[1] = centre + vector3f( -f, -f, f );
    m_vertices[2] = centre + vector3f( -f, f, -f );
    m_vertices[3] = centre + vector3f( f, -f, -f );

    m_faces.resize( 4 );
    m_faces[0].set( 0, 2, 1 );
    m_faces[1].set( 0, 1, 3 );
    m_faces[2].set( 0, 3, 2 );
    m_faces[3].set( 1, 2, 3 );
}

void trimesh3::set_to_icosahedron( const vector3f& centre, float radius ) {
    // Construction based on the golden ratio.  Vertices are permutations of [0,+/-1,+/-phi].
    float phi = 0.5f * ( 1 + sqrtf( 5 ) );
    // Normalize the radius by the length of the vector [0,1,phi] to minimize multiplies and divides later.
    radius /= sqrt( 1 + phi * phi );

    clear();

    m_vertices.resize( 12 );
    m_vertices[0] = centre + vector3f( 0, -radius, -radius * phi );
    m_vertices[1] = centre + vector3f( 0, -radius, radius * phi );
    m_vertices[2] = centre + vector3f( 0, radius, -radius * phi );
    m_vertices[3] = centre + vector3f( 0, radius, radius * phi );

    m_vertices[4] = centre + vector3f( -radius, -radius * phi, 0 );
    m_vertices[5] = centre + vector3f( -radius, radius * phi, 0 );
    m_vertices[6] = centre + vector3f( radius, -radius * phi, 0 );
    m_vertices[7] = centre + vector3f( radius, radius * phi, 0 );

    m_vertices[8] = centre + vector3f( -radius * phi, 0, -radius );
    m_vertices[9] = centre + vector3f( radius * phi, 0, -radius );
    m_vertices[10] = centre + vector3f( -radius * phi, 0, radius );
    m_vertices[11] = centre + vector3f( radius * phi, 0, radius );

    m_faces.resize( 20 );
    m_faces[0].set( 6, 4, 0 );
    m_faces[1].set( 4, 6, 1 );
    m_faces[2].set( 7, 5, 3 );
    m_faces[3].set( 5, 7, 2 );
    m_faces[4].set( 8, 10, 5 );
    m_faces[5].set( 10, 8, 4 );
    m_faces[6].set( 9, 11, 6 );
    m_faces[7].set( 11, 9, 7 );
    m_faces[8].set( 11, 3, 1 );
    m_faces[9].set( 9, 2, 7 );
    m_faces[10].set( 0, 2, 9 );
    m_faces[11].set( 1, 6, 11 );
    m_faces[12].set( 1, 10, 4 );
    m_faces[13].set( 3, 10, 1 );
    m_faces[14].set( 10, 3, 5 );
    m_faces[15].set( 8, 5, 2 );
    m_faces[16].set( 8, 2, 0 );
    m_faces[17].set( 4, 8, 0 );
    m_faces[18].set( 11, 7, 3 );
    m_faces[19].set( 9, 6, 0 );
}

void trimesh3::set_to_rectangular_grid( int horizontalSegments, int verticalSegments,
                                        const ::frantic::graphics2d::boundrect2f& bounds ) {
    clear();
    // First create all the vertices
    m_vertices.reserve( ( horizontalSegments + 1 ) * ( verticalSegments + 1 ) );
    for( int vert_y = 0; vert_y <= verticalSegments; ++vert_y ) {
        for( int vert_x = 0; vert_x <= horizontalSegments; ++vert_x ) {
            vector2f coord = bounds.get_mapped_from_unit_square(
                vector2f( (float)vert_x / horizontalSegments, (float)vert_y / verticalSegments ) );
            m_vertices.push_back( vector3f( coord.x, coord.y, 0 ) );
        }
    }
    // Then create all the faces
    m_faces.reserve( 2 * horizontalSegments * verticalSegments );
    for( int face_y = 0; face_y < verticalSegments; ++face_y ) {
        for( int face_x = 0; face_x < horizontalSegments; ++face_x ) {
            int corners[4];
            // Clockwise around the quad
            corners[0] = face_y * ( horizontalSegments + 1 ) + face_x;
            corners[1] = corners[0] + 1;
            corners[2] = corners[1] + ( horizontalSegments + 1 );
            corners[3] = corners[2] - 1;
            m_faces.push_back( vector3( corners[0], corners[1], corners[3] ) );
            m_faces.push_back( vector3( corners[1], corners[2], corners[3] ) );
        }
    }
}

void trimesh3::set_to_hexagonal_grid( int horizontalSegments, int verticalSegments,
                                      const ::frantic::graphics2d::boundrect2f& bounds ) {
    using namespace ::frantic::graphics2d;

    clear();
    int vertexCount;
    if( verticalSegments % 2 == 0 ) {
        // With an even # (2*n) of vertical segments, it looks something like this:
        //   0   1   2   3
        //     4   5   6
        //   7   8   9   a
        // So the formula is n*m + (n+1)*(m+1)
        vertexCount = ( verticalSegments / 2 ) * horizontalSegments +
                      ( ( verticalSegments / 2 ) + 1 ) * ( horizontalSegments + 1 );
    } else {
        // With an odd # (2*n+1) of vertical segments, it looks something like this:
        //   0   1   2
        //     3   4
        //   5   6   7
        //   a 8   9 b
        // So the formula is (n+1)*m + (n+1)*(m+1) + 2
        vertexCount = ( ( verticalSegments / 2 ) + 1 ) * horizontalSegments +
                      ( ( verticalSegments / 2 ) + 1 ) * ( horizontalSegments + 1 ) + 2;
    }
    // First create all the vertices
    m_vertices.reserve( vertexCount );
    for( int vert_y = 0; vert_y <= verticalSegments; ++vert_y ) {
        if( vert_y % 2 == 0 ) { // On even lines there's 1 more vertex
            for( int vert_x = 0; vert_x <= horizontalSegments; ++vert_x ) {
                vector2f coord = bounds.get_mapped_from_unit_square(
                    vector2f( (float)vert_x / horizontalSegments, (float)vert_y / verticalSegments ) );
                m_vertices.push_back( vector3f( coord.x, coord.y, 0 ) );
            }
        } else {
            for( int vert_x = 0; vert_x < horizontalSegments; ++vert_x ) {
                vector2f coord = bounds.get_mapped_from_unit_square(
                    vector2f( (float)( vert_x + 0.5f ) / horizontalSegments, (float)vert_y / verticalSegments ) );
                m_vertices.push_back( vector3f( coord.x, coord.y, 0 ) );
            }
        }
    }
    if( verticalSegments % 2 != 0 ) {
        // There's two more vertices in the odd verticalSegments case
        vector2f coord = bounds.get_mapped_from_unit_square( vector2f( 0, 1 ) );
        m_vertices.push_back( vector3f( coord.x, coord.y, 0 ) );
        coord = bounds.get_mapped_from_unit_square( vector2f( 1, 1 ) );
        m_vertices.push_back( vector3f( coord.x, coord.y, 0 ) );
    }

    // Then create all the faces
    int faceCount = ( 2 * horizontalSegments - 3 ) * ( verticalSegments - 1 ) + 2 * ( ( verticalSegments + 1 ) / 2 );
    m_faces.reserve( faceCount );
    for( int face_y = 0; face_y < verticalSegments; ++face_y ) {
        // This is the bulk interior of the mesh
        bool evenSegment = ( face_y % 2 == 0 );
        int shortSegmentStart, longSegmentStart;
        longSegmentStart = ( ( face_y + 1 ) / 2 ) * ( 2 * horizontalSegments + 1 );
        if( evenSegment ) {
            shortSegmentStart = longSegmentStart + horizontalSegments + 1;
        } else {
            shortSegmentStart = longSegmentStart - horizontalSegments;
        }
        for( int face_x = 0; face_x < horizontalSegments; ++face_x ) {

            int corners[4];
            // Clockwise around the quad
            corners[0] = longSegmentStart + face_x;
            corners[1] = longSegmentStart + face_x + 1;
            corners[2] = shortSegmentStart + face_x + 1;
            corners[3] = shortSegmentStart + face_x;
            if( !evenSegment )
                std::swap( corners[1], corners[3] );
            m_faces.push_back( vector3( corners[0], corners[1], corners[3] ) );
            if( face_x < horizontalSegments - 1 )
                m_faces.push_back( vector3( corners[1], corners[2], corners[3] ) );
        }
    }
    for( int y = 0; y < verticalSegments - 1; y += 2 ) {
        int longSegmentStart = ( ( y + 1 ) / 2 ) * ( 2 * horizontalSegments + 1 );
        // This is the left and right edge faces to square it off
        m_faces.push_back( vector3( longSegmentStart, longSegmentStart + horizontalSegments + 1,
                                    longSegmentStart + 2 * horizontalSegments + 1 ) );
        m_faces.push_back( vector3( longSegmentStart + horizontalSegments,
                                    longSegmentStart + 3 * horizontalSegments + 1,
                                    longSegmentStart + 2 * horizontalSegments ) );
    }
    if( verticalSegments % 2 != 0 ) {
        int longSegmentStart = ( verticalSegments / 2 ) * ( 2 * horizontalSegments + 1 );
        // These are the bottom left and right corners
        m_faces.push_back( vector3( longSegmentStart, longSegmentStart + horizontalSegments + 1, vertexCount - 2 ) );
        m_faces.push_back( vector3( longSegmentStart + horizontalSegments, vertexCount - 1,
                                    longSegmentStart + 2 * horizontalSegments ) );
    }
}

void trimesh3::add_vertex_and_resize_channels( const vector3f& vertex ) {
    add_vertex( vertex );

    typedef std::map<frantic::tstring, trimesh3_vertex_channel> named_vertex_channels_t;

    // resize vertex channels that do not have custom faces
    for( named_vertex_channels_t::iterator i = m_namedVertexChannels.begin(); i != m_namedVertexChannels.end(); ++i ) {
        if( !i->second.m_hasCustomFaces ) {
            if( i->second.m_data.size() + 1 == vertex_count() ) {
                i->second.m_data.add_element( i->second.m_primitiveSize );
            } else {
                i->second.m_data.resize( m_vertices.size() * i->second.m_primitiveSize );
            }
        }
    }
}

void trimesh3::add_face_and_resize_channels( const vector3& face ) {
    add_face( face );

    typedef std::map<frantic::tstring, trimesh3_face_channel> named_face_channels_t;
    typedef std::map<frantic::tstring, trimesh3_vertex_channel> named_vertex_channels_t;

    // resize faces
    for( named_face_channels_t::iterator i = m_namedFaceChannels.begin(); i != m_namedFaceChannels.end(); ++i ) {
        if( i->second.m_data.size() + 1 == face_count() ) {
            i->second.m_data.add_element( i->second.m_primitiveSize );
        } else {
            i->second.m_data.resize( m_faces.size() * i->second.m_primitiveSize );
        }
    }

    // resize vertex channels that have custom faces
    for( named_vertex_channels_t::iterator i = m_namedVertexChannels.begin(); i != m_namedVertexChannels.end(); ++i ) {
        if( i->second.m_hasCustomFaces ) {
            if( i->second.m_faces.size() + 1 == face_count() ) {
                i->second.m_faces.push_back( vector3( 0 ) );
            } else {
                i->second.m_faces.resize( m_faces.size(), vector3( 0 ) );
            }
        }
    }
}

bool trimesh3::operator==( const trimesh3& rhs ) const {
    if( m_vertices.size() != rhs.m_vertices.size() )
        return false;
    if( m_faces.size() != rhs.m_faces.size() )
        return false;
    if( m_namedVertexChannels.size() != rhs.m_namedVertexChannels.size() )
        return false;
    if( m_namedFaceChannels.size() != rhs.m_namedFaceChannels.size() )
        return false;

    // vertices
    if( m_vertices.size() > 0 ) {
        if( memcmp( &m_vertices[0], &rhs.m_vertices[0], m_vertices.size() * sizeof( vector3f ) ) )
            return false;
    }

    // faces
    if( m_faces.size() > 0 ) {
        if( memcmp( &m_faces[0], &rhs.m_faces[0], m_faces.size() * sizeof( vector3 ) ) )
            return false;
    }

    // vertex channel properties
    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.begin(),
                                                                             ie = m_namedVertexChannels.end();
         i != ie; ++i ) {
        const frantic::tstring& vertexChannelName = i->first;

        std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator irhs =
            rhs.m_namedVertexChannels.find( vertexChannelName );
        if( irhs == rhs.m_namedVertexChannels.end() )
            return false;

        const trimesh3_vertex_channel& left = i->second;
        const trimesh3_vertex_channel& right = irhs->second;
        if( left.arity() != right.arity() )
            return false;
        if( left.data_type() != right.data_type() )
            return false;
        if( left.has_custom_faces() != right.has_custom_faces() )
            return false;
        if( left.face_count() != right.face_count() )
            return false;
        if( left.size() != right.size() )
            return false;
    }

    // face channel properties
    for( std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.begin(),
                                                                           ie = m_namedFaceChannels.end();
         i != ie; ++i ) {
        const frantic::tstring& faceChannelName = i->first;

        std::map<frantic::tstring, trimesh3_face_channel>::const_iterator irhs =
            rhs.m_namedFaceChannels.find( faceChannelName );
        if( irhs == rhs.m_namedFaceChannels.end() )
            return false;

        const trimesh3_face_channel& left = i->second;
        const trimesh3_face_channel& right = irhs->second;
        if( left.arity() != right.arity() )
            return false;
        if( left.data_type() != right.data_type() )
            return false;
        if( left.size() != right.size() )
            return false;
    }

    // vertex channel data
    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.begin(),
                                                                             ie = m_namedVertexChannels.end();
         i != ie; ++i ) {
        const frantic::tstring& vertexChannelName = i->first;
        const_trimesh3_vertex_channel_general_accessor left = get_vertex_channel_general_accessor( vertexChannelName );
        const_trimesh3_vertex_channel_general_accessor right =
            rhs.get_vertex_channel_general_accessor( vertexChannelName );
        if( left.size() ) {
            if( memcmp( left.data( 0 ), right.data( 0 ), left.size() * left.primitive_size() ) )
                return false;
        }
        if( left.has_custom_faces() ) {
            for( std::size_t faceNumber = 0; faceNumber < left.face_count(); ++faceNumber ) {
                if( left.face( faceNumber ) != right.face( faceNumber ) )
                    return false;
            }
        }
    }

    // face channel data
    for( std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.begin(),
                                                                           ie = m_namedFaceChannels.end();
         i != ie; ++i ) {
        const frantic::tstring& faceChannelName = i->first;
        const_trimesh3_face_channel_general_accessor left = get_face_channel_general_accessor( faceChannelName );
        const_trimesh3_face_channel_general_accessor right = rhs.get_face_channel_general_accessor( faceChannelName );
        if( left.size() ) {
            if( memcmp( left.data( 0 ), right.data( 0 ), left.size() * left.primitive_size() ) )
                return false;
        }
    }

    return true;
}

void trimesh3::translate( const vector3f& offset, const vector3f& offsetTimeDerivative ) {
    if( offset.x != 0 || offset.y != 0 || offset.z != 0 ) {
        m_facePlanes.clear();
        m_faceBoundBoxes.clear();

        for( size_t i = 0, count = m_vertices.size(); i < count; ++i )
            m_vertices[i] += offset;
    }

    if( ( offsetTimeDerivative.x != 0 || offsetTimeDerivative.y != 0 || offsetTimeDerivative.z != 0 ) &&
        has_vertex_channel( _T("Velocity") ) ) {
        trimesh3_vertex_channel_accessor<vector3f> velocityChannel =
            get_vertex_channel_accessor<vector3f>( _T("Velocity") );

        for( size_t i = 0, count = velocityChannel.size(); i < count; ++i )
            velocityChannel[i] += offsetTimeDerivative;
    }
}

// Scales the trimesh3 using the given factor and its time derivative
void trimesh3::scale( const vector3f& factor, const vector3f& factorTimeDerivative ) {
    if( factor.x != 1 || factor.y != 1 || factor.z != 1 ) {
        m_facePlanes.clear();
        m_faceBoundBoxes.clear();

        // Get the velocity channel
        trimesh3_vertex_channel_accessor<vector3f> velocityChannel;
        bool hasVelocityChannel = has_vertex_channel( _T("Velocity") );
        if( hasVelocityChannel ) {
            velocityChannel = get_vertex_channel_accessor<vector3f>( _T("Velocity") );

            if( velocityChannel.has_custom_faces() )
                throw runtime_error(
                    "trimesh3::scale() - The Velocity channel in the provided trimesh3 had custom faces.  The "
                    "velocities should correspond 1-1 to the vertices, which means they should share the "
                    "primary geometry faces." );

            if( velocityChannel.size() != vertex_count() )
                throw runtime_error(
                    "trimesh3::scale() - The Velocity channel had a different number of entries than the "
                    "vertex count.  The velocities should correspond 1-1 to the vertices." );
        }

        // Get the normal channel
        trimesh3_vertex_channel_accessor<vector3f> normalChannel;
        bool hasNormalChannel = has_vertex_channel( _T("Normal") );
        if( hasNormalChannel ) {
            normalChannel = get_vertex_channel_accessor<vector3f>( _T("Normal") );
        }

        // Apply the scaling to the positions and velocities
        if( hasVelocityChannel ) {
            for( size_t i = 0, count = m_vertices.size(); i < count; ++i ) {
                // Must scale the velocity before the vertex position, because it depends on the vertex position prior
                // to being scaled.
                velocityChannel[i] = vector3f::component_multiply( factor, velocityChannel[i] );
                velocityChannel[i] += vector3f::component_multiply( factorTimeDerivative, m_vertices[i] );
                m_vertices[i] = vector3f::component_multiply( m_vertices[i], factor );
            }
        } else {
            for( size_t i = 0, count = m_vertices.size(); i < count; ++i ) {
                m_vertices[i] = vector3f::component_multiply( m_vertices[i], factor );
            }
        }

        // Apply the scaling to the normals
        if( hasNormalChannel ) {
            vector3f inverseFactor( 1 / factor.x, 1 / factor.y, 1 / factor.z );
            for( size_t i = 0, count = normalChannel.size(); i < count; ++i ) {
                normalChannel[i] = vector3f::component_multiply( normalChannel[i], inverseFactor );
            }
        }

        if( factor.x * factor.y * factor.z < 0 )
            reverse_face_winding();
    } else if( ( factorTimeDerivative.x != 0 || factorTimeDerivative.y != 0 || factorTimeDerivative.z != 0 ) &&
               has_vertex_channel( _T("Velocity") ) ) {
        // If the scale is one, but the derivative of the factor is non-zero, the velocity channel will
        // still be affected
        trimesh3_vertex_channel_accessor<vector3f> velocityChannel =
            get_vertex_channel_accessor<vector3f>( _T("Velocity") );

        if( velocityChannel.has_custom_faces() )
            throw runtime_error( "trimesh3.scale: The Velocity channel in the provided trimesh3 had custom faces.  The "
                                 "velocities should correspond 1-1 to the vertices, which means they should share the "
                                 "primary geometry faces." );

        if( velocityChannel.size() != vertex_count() )
            throw runtime_error(
                "trimesh3.scale: The Velocity channel had a different number of entries than the vertex "
                "count.  The velocities should correspond 1-1 to the vertices." );

        for( size_t i = 0, count = velocityChannel.size(); i < count; ++i ) {
            velocityChannel[i] += vector3f::component_multiply( factorTimeDerivative, m_vertices[i] );
        }
    }
}

void trimesh3::transform( const transform4f& xform, const transform4f& xformTimeDerivative ) {
    if( !xform.is_identity() ) {

        m_facePlanes.clear();
        m_faceBoundBoxes.clear();

        // Get the velocity channel
        trimesh3_vertex_channel_general_accessor velocityChannel;
        bool hasVelocityChannel = has_vertex_channel( _T("Velocity") );
        if( hasVelocityChannel ) {
            velocityChannel = get_vertex_channel_general_accessor( _T("Velocity") );

            if( velocityChannel.has_custom_faces() )
                throw runtime_error(
                    "trimesh3::transform() - The Velocity channel in the provided trimesh3 had custom faces.  "
                    "The velocities should correspond 1-1 to the vertices, which means they should share the "
                    "primary geometry faces." );

            if( velocityChannel.size() != vertex_count() )
                throw runtime_error(
                    "trimesh3::transform() - The Velocity channel of the mesh has a different number of "
                    "entries than the vertex count.  The velocities should correspond 1-1 to the vertices." );

            if( velocityChannel.arity() != 3 )
                throw runtime_error( "trimesh3::transform() - The Velocity channel of the mesh has arity " +
                                     boost::lexical_cast<string>( velocityChannel.arity() ) +
                                     ", it should have arity 3." );
        }

        // Get the normal channel
        trimesh3_vertex_channel_general_accessor normalChannel;
        bool hasNormalChannel = has_vertex_channel( _T("Normal") );
        if( hasNormalChannel ) {
            normalChannel = get_vertex_channel_general_accessor( _T("Normal") );

            if( normalChannel.arity() != 3 )
                throw runtime_error( "trimesh3:transform() - The Normal channel of the mesh has arity " +
                                     boost::lexical_cast<string>( normalChannel.arity() ) +
                                     ", it should have arity 3." );
        }

        // Apply the transform to the positions and velocities
        if( hasVelocityChannel ) {
            // Get functions for converting to/from the velocity channel
            channel_type_convertor_function_t convertFromChannel =
                get_channel_type_convertor_function( velocityChannel.data_type(), data_type_float32, _T("Velocity") );
            channel_type_convertor_function_t convertToChannel =
                get_channel_type_convertor_function( data_type_float32, velocityChannel.data_type(), _T("Velocity") );
            vector3f velocity;
            for( size_t i = 0, count = m_vertices.size(); i < count; ++i ) {
                // Must transform the velocity before the vertex position, because it depends on the vertex position
                // prior to transformation
                // TODO: Should make a cvt vertex channel accessor for the trimesh3.
                convertFromChannel( reinterpret_cast<char*>( &velocity ), velocityChannel.data( i ), 3 );
                velocity = xform.transform_no_translation( velocity ) + xformTimeDerivative * m_vertices[i];
                convertToChannel( velocityChannel.data( i ), reinterpret_cast<char*>( &velocity ), 3 );
                m_vertices[i] = xform * m_vertices[i];
            }
        } else {
            for( size_t i = 0, count = m_vertices.size(); i < count; ++i ) {
                m_vertices[i] = xform * m_vertices[i];
            }
        }

        // Apply the transform to the normals
        if( hasNormalChannel ) {
            // Get functions for converting to/from the normal channel
            channel_type_convertor_function_t convertFromChannel =
                get_channel_type_convertor_function( normalChannel.data_type(), data_type_float32, _T("Normal") );
            channel_type_convertor_function_t convertToChannel =
                get_channel_type_convertor_function( data_type_float32, normalChannel.data_type(), _T("Normal") );
            transform4f xformInverse = xform.to_inverse();
            vector3f normal;
            for( size_t i = 0, count = normalChannel.size(); i < count; ++i ) {
                // TODO: Should make a cvt vertex channel accessor for the trimesh3.
                convertFromChannel( reinterpret_cast<char*>( &normal ), normalChannel.data( i ), 3 );
                normal = xformInverse.transpose_transform_no_translation( normal );
                // TODO: Should we renormalize?
                convertToChannel( normalChannel.data( i ), reinterpret_cast<char*>( &normal ), 3 );
            }
        }

        if( xform.is_orientation_inverting() )
            reverse_face_winding();
    } else if( !xformTimeDerivative.is_zero() && has_vertex_channel( _T("Velocity") ) ) {
        // If the matrix is the identity, but the derivative of the matrix is non-zero, the velocity channel will
        // still be affected
        trimesh3_vertex_channel_general_accessor velocityChannel =
            get_vertex_channel_general_accessor( _T("Velocity") );

        if( velocityChannel.has_custom_faces() )
            throw runtime_error(
                "trimesh3:transform() - The Velocity channel in the provided trimesh3 had custom faces.  "
                "The velocities should correspond 1-1 to the vertices, which means they should share the "
                "primary geometry faces." );

        if( velocityChannel.size() != vertex_count() )
            throw runtime_error(
                "trimesh3:transform() - The Velocity channel had a different number of entries than the "
                "vertex count.  The velocities should correspond 1-1 to the vertices." );

        if( velocityChannel.arity() != 3 )
            throw runtime_error( "trimesh3:transform() - The Velocity channel of the mesh has arity " +
                                 boost::lexical_cast<string>( velocityChannel.arity() ) + ", it should have arity 3." );

        // Get functions for converting to/from the velocity channel
        channel_type_convertor_function_t convertFromChannel =
            get_channel_type_convertor_function( velocityChannel.data_type(), data_type_float32, _T("Velocity") );
        channel_type_convertor_function_t convertToChannel =
            get_channel_type_convertor_function( data_type_float32, velocityChannel.data_type(), _T("Velocity") );

        vector3f velocity;
        for( size_t i = 0, count = velocityChannel.size(); i < count; ++i ) {
            // TODO: Should make a cvt vertex channel accessor for the trimesh3.
            convertFromChannel( reinterpret_cast<char*>( &velocity ), velocityChannel.data( i ), 3 );
            velocity = xform.transform_no_translation( velocity ) + xformTimeDerivative * m_vertices[i];
            convertToChannel( velocityChannel.data( i ), reinterpret_cast<char*>( &velocity ), 3 );
        }
    }
}

/**
 * Builds a vector of adjacenies indexed on the vertex indices. Each element contains the set of
 * of adjacent verts
 *
 * @returns vector of adjacent vertex sets
 */
void trimesh3::get_adjacency_list( std::vector<std::set<int>>& outAdjList ) const {
    outAdjList.clear();
    outAdjList.resize( m_vertices.size() );

    for( unsigned i = 0; i < m_faces.size(); ++i ) {
        vector3 face = m_faces[i];

        outAdjList[face.x].insert( face.y );
        outAdjList[face.x].insert( face.z );

        outAdjList[face.y].insert( face.x );
        outAdjList[face.y].insert( face.z );

        outAdjList[face.z].insert( face.x );
        outAdjList[face.z].insert( face.y );
    }
}

float trimesh3::get_minimum_edge_length() const {
    float minimumEdgeLengthSquared = ( std::numeric_limits<float>::max )();
    for( unsigned i = 0; i < m_faces.size(); ++i ) {
        vector3 face = m_faces[i];

        vector3f a = m_vertices[face.x];
        vector3f b = m_vertices[face.y];
        vector3f c = m_vertices[face.z];

        minimumEdgeLengthSquared = ( std::min )( minimumEdgeLengthSquared, vector3f::distance_squared( a, b ) );
        minimumEdgeLengthSquared = ( std::min )( minimumEdgeLengthSquared, vector3f::distance_squared( a, c ) );
        minimumEdgeLengthSquared = ( std::min )( minimumEdgeLengthSquared, vector3f::distance_squared( b, c ) );
    }

    return sqrtf( minimumEdgeLengthSquared );
}

// TODO: We could implement this one separately instead of calling the transformed case, for performance.
void trimesh3::combine( const trimesh3& mesh ) {
    std::vector<std::pair<frantic::tstring, vector_type>> transformChannels;
    transformChannels.push_back( std::pair<frantic::tstring, vector_type>( _T("Velocity"), VECTOR ) );
    transformChannels.push_back( std::pair<frantic::tstring, vector_type>( _T("Normal"), NORMAL ) );
    combine( transform4f::identity(), transform4f::zero(), mesh, transformChannels );
}

void trimesh3::combine( const transform4f& xform, const trimesh3& mesh ) {
    std::vector<std::pair<frantic::tstring, vector_type>> transformChannels;
    transformChannels.push_back( std::pair<frantic::tstring, vector_type>( _T("Velocity"), VECTOR ) );
    transformChannels.push_back( std::pair<frantic::tstring, vector_type>( _T("Normal"), NORMAL ) );
    combine( xform, transform4f::zero(), mesh, transformChannels );
}

void trimesh3::combine( const transform4f& xform, const trimesh3& mesh,
                        const std::vector<std::pair<frantic::tstring, vector_type>> transformChannels ) {
    combine( xform, transform4f::zero(), mesh, transformChannels );
}

void trimesh3::combine( const transform4f& xform, const transform4f& xformTimeDerivative, const trimesh3& mesh ) {
    std::vector<std::pair<frantic::tstring, vector_type>> transformChannels;
    transformChannels.push_back( std::pair<frantic::tstring, vector_type>( _T("Velocity"), VECTOR ) );
    transformChannels.push_back( std::pair<frantic::tstring, vector_type>( _T("Normal"), NORMAL ) );
    combine( xform, xformTimeDerivative, mesh, transformChannels );
}

// This merges all the faces and vertices of the given mesh into the current mesh, applying the provided transform to
// the mesh on the fly.
void trimesh3::combine( const transform4f& xform, const transform4f& xformTimeDerivative, const trimesh3& mesh,
                        const std::vector<std::pair<frantic::tstring, vector_type>> transformChannels ) {
    if( mesh.is_empty() )
        return;

    // First go through the input mesh and make sure every vertex channel there has a corresponding vertex channel in
    // this mesh.
    for( map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = mesh.m_namedVertexChannels.begin(),
                                                                        iterEnd = mesh.m_namedVertexChannels.end();
         i != iterEnd; ++i ) {
        if( !has_vertex_channel( i->first ) ) {
            add_vertex_channel_raw( i->first, i->second.arity(), i->second.data_type() );
            trimesh3_vertex_channel_general_accessor acc = get_vertex_channel_general_accessor( i->first );
            memset( acc.data( 0 ), 0, m_vertices.size() * i->second.m_primitiveSize );
        }
    }

    // First go through the input mesh and make sure every face channel there has a corresponding face channel in this
    // mesh.
    for( map<frantic::tstring, trimesh3_face_channel>::const_iterator i = mesh.m_namedFaceChannels.begin(),
                                                                      iterEnd = mesh.m_namedFaceChannels.end();
         i != iterEnd; ++i ) {
        if( !has_face_channel( i->first ) ) {
            add_face_channel_raw( i->first, i->second.arity(), i->second.data_type() );
            trimesh3_face_channel_general_accessor acc = get_face_channel_general_accessor( i->first );
            memset( acc.data( 0 ), 0, m_faces.size() * i->second.m_primitiveSize );
        }
    }

    // Then go through this mesh and copy the data of each vertex channel one by one
    for( map<frantic::tstring, trimesh3_vertex_channel>::iterator i = m_namedVertexChannels.begin(),
                                                                  iterEnd = m_namedVertexChannels.end();
         i != iterEnd; ++i ) {
        // If the input mesh doesn't have this channel, then just fill it out to all zeros
        if( !mesh.has_vertex_channel( i->first ) ) {
            int vertexOffset = (int)i->second.size();
            // Add a block of memory of the right size, and initialize it to all zero
            size_t dataSizeToAdd = i->second.primitive_size() * mesh.vertex_count();
            memset( i->second.m_data.add_element( dataSizeToAdd ), 0, dataSizeToAdd );
            // If the channel has custom faces, we need to add to them as well, otherwise leave them as the main faces
            // are dealt with later
            if( i->second.has_custom_faces() ) {
                // Add faces matching the faces of the input mesh
                for( size_t face = 0, count = mesh.face_count(); face < count; ++face ) {
                    i->second.m_faces.push_back( mesh.m_faces[face] + vector3( vertexOffset ) );
                }
            }
        }
        // If the input mesh does have this channel, have to make sure it's compatible, and then copy it.
        else {
            int type = -1;

            for( size_t j = 0; j < transformChannels.size(); j++ )
                if( i->first == transformChannels[j].first )
                    type = transformChannels[j].second;

            switch( type ) {
            case VECTOR: {
                const_trimesh3_vertex_channel_general_accessor inputVelocityChannel =
                    mesh.get_vertex_channel_general_accessor( i->first );
                trimesh3_vertex_channel_general_accessor outputVelocityChannel =
                    get_vertex_channel_general_accessor( i->first );

                if( inputVelocityChannel.has_custom_faces() || outputVelocityChannel.has_custom_faces() )
                    throw runtime_error(
                        "trimesh3::combine() - The " + frantic::strings::to_string( i->first ) +
                        " channel in the provided trimesh3 had custom faces.  The velocities should correspond "
                        "1-1 to the vertices, which means they should share the primary geometry faces." );

                if( inputVelocityChannel.size() != mesh.vertex_count() ||
                    outputVelocityChannel.size() != vertex_count() )
                    throw runtime_error(
                        "trimesh3::combine() - The " + frantic::strings::to_string( i->first ) +
                        " channel had a different number of entries than the vertex count.  The velocities "
                        "should correspond 1-1 to the vertices." );

                if( inputVelocityChannel.arity() != 3 || outputVelocityChannel.arity() != 3 )
                    throw runtime_error( "trimesh3::combine() - The " + frantic::strings::to_string( i->first ) +
                                         " channel of one of the inputs had an arity not equal to 3." );

                // Get functions for converting to/from the velocity channel
                channel_type_convertor_function_t convertFromChannel = get_channel_type_convertor_function(
                    inputVelocityChannel.data_type(), data_type_float32, i->first );
                channel_type_convertor_function_t convertToChannel = get_channel_type_convertor_function(
                    data_type_float32, outputVelocityChannel.data_type(), i->first );
                vector3f velocity;
                for( size_t j = 0, count = inputVelocityChannel.size(); j < count; ++j ) {
                    // TODO: Should make a cvt vertex channel accessor for the trimesh3.
                    convertFromChannel( reinterpret_cast<char*>( &velocity ), inputVelocityChannel.data( j ), 3 );
                    velocity = xform.transform_no_translation( velocity ) + xformTimeDerivative * mesh.get_vertex( j );
                    convertToChannel( outputVelocityChannel.add_vertex(), reinterpret_cast<char*>( &velocity ), 3 );
                }
            } break;
            case NORMAL: {
                const_trimesh3_vertex_channel_accessor<vector3f> inputNormalChannel =
                    mesh.get_vertex_channel_accessor<vector3f>( i->first );

                // Make sure that the channel we're copying to has custom faces if necessary.
                if( inputNormalChannel.has_custom_faces() && !i->second.has_custom_faces() )
                    set_vertex_channel_custom_faces( i->first, true );

                // Get the output accessor.  This has to be after we set the custom faces flag, because changing custom
                // faces invalidates any accessors pointing to that channel.
                trimesh3_vertex_channel_general_accessor outputNormalChannel =
                    get_vertex_channel_general_accessor( i->first );

                if( inputNormalChannel.arity() != 3 || outputNormalChannel.arity() != 3 )
                    throw runtime_error( "trimesh3::combine() - The " + frantic::strings::to_string( i->first ) +
                                         " channel of one of the inputs had an arity not equal to 3." );

                int vertexOffset = (int)outputNormalChannel.size();

                // Get functions for converting to/from the velocity channel
                channel_type_convertor_function_t convertFromChannel =
                    get_channel_type_convertor_function( inputNormalChannel.data_type(), data_type_float32, i->first );
                channel_type_convertor_function_t convertToChannel =
                    get_channel_type_convertor_function( data_type_float32, outputNormalChannel.data_type(), i->first );
                vector3f normal;
                transform4f xformInverse = xform.to_inverse();
                for( size_t j = 0, je = inputNormalChannel.size(); j != je; ++j ) {
                    // TODO: Should make a cvt vertex channel accessor for the trimesh3.
                    convertFromChannel( reinterpret_cast<char*>( &normal ), inputNormalChannel.data( j ), 3 );
                    normal = xformInverse.transpose_transform_no_translation( normal );
                    // TODO: Should we renormalize?
                    convertToChannel( outputNormalChannel.add_vertex(), reinterpret_cast<char*>( &normal ), 3 );
                }

                // If the channel has custom faces, we need to add to them as well, otherwise leave them as the main
                // faces are dealt with later
                if( i->second.has_custom_faces() ) {
                    // Add faces matching the faces of the input mesh
                    for( size_t face = 0, faceEnd = mesh.face_count(); face != faceEnd; ++face ) {
                        i->second.m_faces.push_back( inputNormalChannel.face( face ) + vector3( vertexOffset ) );
                    }
                }
            } break;
            default: {
                const_trimesh3_vertex_channel_general_accessor inputChannel =
                    mesh.get_vertex_channel_general_accessor( i->first );

                // If the arity of the input and output don't match, throw an exception
                if( inputChannel.arity() != i->second.arity() )
                    throw runtime_error(
                        "trimesh3::combine() - The \"" + frantic::strings::to_string( i->first ) +
                        "\" channel of the provided trimesh3 had a different arity than in the destination trimesh3." );

                // Get a convertor function for transferring the data
                channels::channel_type_convertor_function_t convertType = channels::get_channel_type_convertor_function(
                    inputChannel.data_type(), i->second.data_type(), i->first );

                // Make sure that the channel we're copying to has custom faces if necessary.
                if( inputChannel.has_custom_faces() && !i->second.has_custom_faces() )
                    set_vertex_channel_custom_faces( i->first, true );

                int vertexOffset = (int)i->second.size();
                size_t numData;
                if( inputChannel.has_custom_faces() )
                    numData = inputChannel.size();
                else
                    numData = mesh.vertex_count();

                // Add a block of memory of the right size, and copy the channel data
                size_t bytesToAdd = i->second.primitive_size() * numData, elementsToAdd = i->second.arity() * numData;
                convertType( i->second.m_data.add_element( bytesToAdd ), inputChannel.data( 0 ), elementsToAdd );

                // If the channel has custom faces, we need to add to them as well, otherwise leave them as the main
                // faces are dealt with later
                if( i->second.has_custom_faces() ) {
                    // Add faces matching the faces of the input mesh
                    for( size_t face = 0, faceEnd = mesh.face_count(); face != faceEnd; ++face ) {
                        i->second.m_faces.push_back( inputChannel.face( face ) + vector3( vertexOffset ) );
                    }
                }
            } break;
            }
        }
    }

    // Now go through this mesh and copy the data of each face channel one by one
    for( map<frantic::tstring, trimesh3_face_channel>::iterator i = m_namedFaceChannels.begin(),
                                                                iterEnd = m_namedFaceChannels.end();
         i != iterEnd; ++i ) {
        // If the input mesh doesn't have this channel, then just fill it out to all zeros
        if( !mesh.has_face_channel( i->first ) ) {
            // Add a block of memory of the right size, and initialize it to all zero
            size_t dataSizeToAdd = i->second.primitive_size() * mesh.face_count();
            memset( i->second.m_data.add_element( dataSizeToAdd ), 0, dataSizeToAdd );
        }
        // If the input mesh does have this channel, have to make sure it's compatible, and then copy it.
        else {
            const_trimesh3_face_channel_general_accessor inputChannel =
                mesh.get_face_channel_general_accessor( i->first );

            // If the arity of the input and output don't match, throw an exception
            if( inputChannel.arity() != i->second.arity() )
                throw runtime_error( "trimesh3::combine() - The \"" + frantic::strings::to_string( i->first ) +
                                     "\" face channel of the provided trimesh3 had a different arity than in the "
                                     "destination trimesh3." );

            // Get a convertor function for transferring the data
            channels::channel_type_convertor_function_t convertType = channels::get_channel_type_convertor_function(
                inputChannel.data_type(), i->second.data_type(), i->first );

            // Add a block of memory of the right size, and copy the channel data
            size_t bytesToAdd = i->second.primitive_size() * mesh.face_count(),
                   elementsToAdd = i->second.arity() * mesh.face_count();
            convertType( i->second.m_data.add_element( bytesToAdd ), inputChannel.data( 0 ), elementsToAdd );
        }
    }

    // NOTE: We call reserve() only when adding this mesh increases the vert or face count significantly.  This
    //       is so the normal resizing happens when adding small meshes to a big one, like building a particle
    //       mesh from lots of small particle meshes.

    // Clear the cached face data
    m_facePlanes.clear();
    m_faceBoundBoxes.clear();

    int vertexOffset = (int)m_vertices.size();
    int totalVertexCount = int( m_vertices.size() + mesh.m_vertices.size() );

    // First append all the vertices from the source mesh onto the destination mesh
    if( (std::size_t)totalVertexCount > m_vertices.size() * 14 / 10 )
        m_vertices.reserve( totalVertexCount );
    if( xform.is_identity() ) {
        for( unsigned i = 0; i < mesh.m_vertices.size(); ++i )
            m_vertices.push_back( mesh.m_vertices[i] );
    } else {
        for( unsigned i = 0; i < mesh.m_vertices.size(); ++i )
            m_vertices.push_back( xform * mesh.m_vertices[i] );
    }

    // Combine the faces, adjusting the vertex indices in the combined mesh
    std::size_t totalFaceCount = m_faces.size() + mesh.m_faces.size();
    if( (std::size_t)totalFaceCount > m_faces.size() * 14 / 10 )
        m_faces.reserve( totalFaceCount );
    for( unsigned face = 0; face < mesh.m_faces.size(); ++face ) {
        m_faces.push_back( mesh.m_faces[face] + vector3( vertexOffset ) );
    }
}

// Builds the vertex normals based on incident angle.  By default it puts it in a channel called "Normal"
// TODO: Eventually, when named face channels exist, we can use a "SmoothingGroup" channel to determine which edges are
// smoothed.
//       This would be done by generating faces based on the smoothing groups, then running the algorithm below as-is.
void trimesh3::build_vertex_normals( const frantic::tstring& vertexChannelName, bool normalizeNormals ) {
    // Add this channel if it's not already defined
    if( !has_vertex_channel( vertexChannelName ) )
        add_vertex_channel<vector3f>( vertexChannelName );

    trimesh3_vertex_channel_accessor<vector3f> normals = get_vertex_channel_accessor<vector3f>( vertexChannelName );

    // Make sure all the normals start as 0
    for( unsigned i = 0; i < vertex_count(); ++i )
        normals[i] = vector3f( 0.0f );

    // Add the contributions of all the faces to the normals array
    for( size_t faceIndex = 0; faceIndex < face_count(); ++faceIndex ) {
        // Get the triangle
        const vector3& face = get_face( faceIndex );
        const vector3f& a = get_vertex( face.x );
        const vector3f& b = get_vertex( face.y );
        const vector3f& c = get_vertex( face.z );

        // Get the normal, and use the triangle angles as the weights
        vector3f geoNormal = triangle_normal( a, b, c );
        vector3f weights = get_triangle_angles( a, b, c );
        normals[face.x] += weights.x * geoNormal;
        normals[face.y] += weights.y * geoNormal;
        normals[face.z] += weights.z * geoNormal;
    }

    // Normalize all the resulting vectors
    if( normalizeNormals ) {
        for( unsigned i = 0; i < normals.size(); ++i )
            normals[i].normalize();
    }
}

void trimesh3::build_face_normals( const frantic::tstring& faceChannelName ) {
    if( !has_face_channel( faceChannelName ) ) {
        add_face_channel<vector3f>( faceChannelName );
    }

    trimesh3_face_channel_accessor<vector3f> normals = get_face_channel_accessor<vector3f>( faceChannelName );

    for( size_t faceIndex = 0; faceIndex < face_count(); ++faceIndex ) {
        const vector3& face = get_face( faceIndex );
        const vector3f& a = get_vertex( face.x );
        const vector3f& b = get_vertex( face.y );
        const vector3f& c = get_vertex( face.z );

        normals[faceIndex] = triangle_normal( a, b, c );
    }
}

void trimesh3::linear_interpolation( const trimesh3& mesh1, const trimesh3& mesh2, const float& alpha ) {

    if( alpha < 0 || alpha > 1 )
        throw std::runtime_error(
            "trimesh3::linear_interpolation - Cannot interpolate, alpha value should bein the range [0,1]" );

    // Check that the topology is consistent
    if( mesh1.vertex_count() != mesh2.vertex_count() )
        throw std::runtime_error(
            "trimesh3::linear_interpolation - Cannot interpolate, inconsistent vertex count.  First mesh has " +
            lexical_cast<string>( mesh1.vertex_count() ) + " verts and second mesh has " +
            lexical_cast<string>( mesh2.vertex_count() ) + " verts." );
    if( mesh1.face_count() != mesh2.face_count() )
        throw std::runtime_error(
            "trimesh3::linear_interpolation - Cannot interpolate, inconsistent face count.  First mesh has " +
            lexical_cast<string>( mesh1.face_count() ) + " faces and second mesh has " +
            lexical_cast<string>( mesh2.face_count() ) + " faces." );
    if( mesh1.face_count() > 0 ) {
        if( std::memcmp( &( mesh1.faces_ref()[0] ), &( mesh2.faces_ref()[0] ),
                         mesh1.face_count() * sizeof( vector3 ) ) != 0 )
            throw std::runtime_error( "trimesh3::linear_interpolation - Cannot interpolate, inconsistent face data." );
    }

    trimesh3 result;
    result.set_vertex_count( mesh1.vertex_count() );
    result.set_face_count( mesh1.face_count() );

    // interpolate the verts
    for( unsigned i = 0; i < result.vertices_ref().size(); i++ )
        result.vertices_ref()[i] =
            mesh1.vertices_ref()[i] + alpha * ( mesh2.vertices_ref()[i] - mesh1.vertices_ref()[i] );

    // copy the faces
    if( mesh1.face_count() > 0 ) {
        memcpy( &( result.faces_ref()[0] ), &( mesh1.faces_ref()[0] ), mesh1.face_count() * sizeof( vector3 ) );
    }

    std::vector<frantic::tstring> channels;
    mesh1.get_vertex_channel_names( channels );

    // Interpolate channels.  Interpolation occurs only between channels which exist in both meshes.
    for( size_t i = 0; i < channels.size(); i++ ) {
        // if both meshes don't have the channel, then skip it
        if( !mesh2.has_vertex_channel( channels[i] ) )
            continue; // throw runtime_error( "todo write me" );

        const_trimesh3_vertex_channel_general_accessor channel1 =
            mesh1.get_vertex_channel_general_accessor( channels[i] );
        const_trimesh3_vertex_channel_general_accessor channel2 =
            mesh2.get_vertex_channel_general_accessor( channels[i] );

        if( channel1.arity() != channel2.arity() )
            throw std::runtime_error(
                "trimesh3::linear_interpolation - Cannot interpolate vertex channel " +
                frantic::strings::to_string( channels[i] ) + ".  The first mesh's channel has arity " +
                lexical_cast<string>( channel1.arity() ) + " and the second mesh's channel has arity " +
                lexical_cast<string>( channel2.arity() ) );
        if( channel1.data_type() != channel2.data_type() )
            throw std::runtime_error(
                "trimesh3::linear_interpolation - Cannot interpolate vertex channel " +
                frantic::strings::to_string( channels[i] ) + ".  The first mesh's channel is of type " +
                frantic::strings::to_string( channels::channel_data_type_str( channel1.data_type() ) ) +
                " and the second mesh's channel is of type " +
                frantic::strings::to_string( channels::channel_data_type_str( channel1.data_type() ) ) );

        if( channel1.has_custom_faces() != channel2.has_custom_faces() )
            throw runtime_error( "trimesh3::linear_interpolation - Cannot interpolate vertex channel " +
                                 frantic::strings::to_string( channels[i] ) +
                                 ".  The first mesh's channel has_custom_faces() = " +
                                 lexical_cast<string>( channel1.has_custom_faces() ) +
                                 " and the second mesh's chanel has_custom_faces() = " +
                                 lexical_cast<string>( channel2.has_custom_faces() ) );

        if( channel1.has_custom_faces() ) {
            if( channel1.size() != channel2.size() )
                throw std::runtime_error( "trimesh3::linear_interpolation - Cannot interpolate vertex channel " +
                                          frantic::strings::to_string( channels[i] ) +
                                          " with custom faces, inconsistent vertex count.  First mesh has " +
                                          lexical_cast<string>( channel2.size() ) + " verts and second mesh has " +
                                          lexical_cast<string>( channel2.size() ) + " verts." );
            if( channel1.face_count() > 0 ) {
                if( std::memcmp( &( channel1.face( 0 ) ), &( channel2.face( 0 ) ),
                                 channel1.face_count() * sizeof( vector3 ) ) != 0 )
                    throw std::runtime_error( "trimesh3::linear_interpolation - Cannot interpolate vertex channel " +
                                              frantic::strings::to_string( channels[i] ) +
                                              ".  Inconsistent face data." );
            }
        }

        // add the channel from the first mesh and get the accessor
        result.add_vertex_channel_raw( channels[i], channel1.arity(), channel1.data_type(), channel1.size(),
                                       channel1.has_custom_faces() );
        trimesh3_vertex_channel_general_accessor newChannel = result.get_vertex_channel_general_accessor( channels[i] );

        if( newChannel.has_custom_faces() ) {
            if( channel1.face_count() > 0 ) {
                memcpy( &( newChannel.face( 0 ) ), &( channel1.face( 0 ) ), channel1.face_count() * sizeof( vector3 ) );
            }
        }

        if( newChannel.size() > 0 ) {
            // use the weighted sum function of the accessor to interpolate between the two given channels
            channel_weighted_sum_combine_function_t ws = newChannel.m_weightedSumFunction;
            const char* data[2];
            data[0] = channel1.data( 0 );
            data[1] = channel2.data( 0 );
            float weights[2];
            weights[0] = 1.f - alpha;
            weights[1] = alpha;
            ws( weights, data, 2, newChannel.size() * newChannel.arity(), newChannel.data( 0 ) );
        }
    }

    std::vector<frantic::tstring> faceChannels;
    mesh1.get_face_channel_names( faceChannels );

    for( std::size_t i = 0; i < faceChannels.size(); ++i ) {
        const frantic::tstring& channelName = faceChannels[i];

        if( !mesh2.has_face_channel( channelName ) )
            continue;

        const_trimesh3_face_channel_general_accessor channel1 = mesh1.get_face_channel_general_accessor( channelName );
        const_trimesh3_face_channel_general_accessor channel2 = mesh2.get_face_channel_general_accessor( channelName );

        if( channel1.arity() != channel2.arity() )
            throw std::runtime_error(
                "trimesh3::linear_interpolation - Cannot interpolate face channel " +
                frantic::strings::to_string( channelName ) + ".  The first mesh's channel has arity " +
                lexical_cast<string>( channel1.arity() ) + " and the second mesh's channel has arity " +
                lexical_cast<string>( channel2.arity() ) );
        if( channel1.data_type() != channel2.data_type() )
            throw std::runtime_error(
                "trimesh3::linear_interpolation - Cannot interpolate face channel " +
                frantic::strings::to_string( channelName ) + ".  The first mesh's channel is of type " +
                frantic::strings::to_string( channels::channel_data_type_str( channel1.data_type() ) ) +
                " and the second mesh's channel is of type " +
                frantic::strings::to_string( channels::channel_data_type_str( channel1.data_type() ) ) );

        result.add_face_channel_raw( channelName, channel1.arity(), channel1.data_type() );
        trimesh3_face_channel_general_accessor newChannel = result.get_face_channel_general_accessor( channelName );

        if( newChannel.size() > 0 ) {
            channel_weighted_sum_combine_function_t ws = newChannel.m_weightedSumFunction;
            const char* data[2];
            data[0] = channel1.data( 0 );
            data[1] = channel2.data( 0 );
            float weights[2];
            weights[0] = 1.f - alpha;
            weights[1] = alpha;
            ws( weights, data, 2, newChannel.size() * newChannel.arity(), newChannel.data( 0 ) );
        }
    }

    // swap in the new mesh
    swap( result );
}

void trimesh3::velocity_offset( float timeOffset ) {

    // move the vertices based on the velocity
    if( timeOffset != 0 && has_vertex_channel( _T("Velocity") ) ) {
        const_trimesh3_vertex_channel_general_accessor velocity = get_vertex_channel_general_accessor( _T("Velocity") );
        if( velocity.arity() != 3 )
            throw std::runtime_error(
                "trimesh3::velocity_offset: The velocity channel from the input mesh had an arity "
                "different from 3. The arity is " +
                boost::lexical_cast<string>( velocity.arity() ) + "." );
        if( velocity.has_custom_faces() )
            throw std::runtime_error(
                "trimesh3::velocity_offset: The velocity channel of the input mesh has custom faces, "
                "which means it can't be applied for motion blur." );

        if( velocity.data_type() == data_type_float32 ) {
            for( size_t i = 0; i < velocity.size(); ++i ) {
                m_vertices[i] += timeOffset * ( *(vector3f*)( velocity.data( i ) ) );
            }
        } else {
            vector3f velocityVector;
            channel_type_convertor_function_t convertFromChannel =
                get_channel_type_convertor_function( velocity.data_type(), data_type_float32, _T("Velocity") );
            for( size_t i = 0; i < velocity.size(); ++i ) {
                convertFromChannel( reinterpret_cast<char*>( &velocityVector ), velocity.data( i ), 3 );

                m_vertices[i] += timeOffset * velocityVector;
            }
        }
    }
}

size_t trimesh3::count_named_vertex_channels_with_custom_faces() const {
    size_t count = 0;
    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator it = m_namedVertexChannels.begin();
         it != m_namedVertexChannels.end(); ++it ) {
        if( it->second.has_custom_faces() ) {
            ++count;
        }
    }
    return count;
}

void trimesh3::erase_all_named_channels() {
    m_namedVertexChannels.clear();
    m_namedFaceChannels.clear();
}

bool trimesh3::check_channel_sizes( std::ostream& out, frantic::logging::progress_logger* logger ) const {

    frantic::logging::null_progress_logger nullLogger;

    if( !logger ) {
        logger = &nullLogger;
    }

    typedef std::map<frantic::tstring, trimesh3_vertex_channel> named_vertex_channel_map_t;

    BOOST_FOREACH( const named_vertex_channel_map_t::value_type& channelIt, m_namedVertexChannels ) {
        const frantic::tstring& channelName = channelIt.first;
        const trimesh3_vertex_channel& channel = channelIt.second;

        if( channel.has_custom_faces() && channel.face_count() != face_count() ) {
            out << "Vertex channel \"" << frantic::strings::to_string( channelName )
                << "\" with custom faces has an incorrect face count (was " << channel.face_count() << ", expected "
                << face_count() << ")" << std::endl;
            return false;
        } else if( channel.size() != vertex_count() ) {
            out << "Vertex channel \"" << frantic::strings::to_string( channelName )
                << "\" has an incorrect vertex count (was " << channel.size() << ", expected " << vertex_count() << ")"
                << std::endl;
            return false;
        }
    }

    typedef std::map<frantic::tstring, trimesh3_face_channel> named_face_channel_map_t;

    BOOST_FOREACH( const named_face_channel_map_t::value_type& channelIt, m_namedFaceChannels ) {
        const frantic::tstring& channelName = channelIt.first;
        const trimesh3_face_channel& channel = channelIt.second;

        if( channel.size() != face_count() ) {
            out << "Face channel \"" << frantic::strings::to_string( channelName )
                << "\" has an incorrect face count (was " << channel.size() << ", expected " << face_count() << ")"
                << std::endl;
            return false;
        }
    }

    return true;
}

// The original implementation was quite slow due to using a map of all edges. This is my attempt to
// optimize the method by using a sorted array, and parallelizing some parts
bool trimesh3::check_duplicate_edges( std::ostream& out, frantic::logging::progress_logger* logger ) const {
    tbb::task_scheduler_init taskScheduleInit;

    frantic::logging::null_progress_logger nullLogger;
    if( !logger ) {
        logger = &nullLogger;
    }

    std::vector<edge_with_face> edgeSet( face_count() * 3 );

    const size_t updateInterval = ( edgeSet.size() / 100 ) + 1;

    logger->set_title( _T( "Checking duplicate edges" ) );

    logger->update_progress( 0.0f );

    const size_t customFaceChannelCount = count_named_vertex_channels_with_custom_faces() + 1;
    const float channelProgressStride = 100.0f / customFaceChannelCount;

    {
        frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, 0.0f, channelProgressStride );

        collect_edges_with_faces( *this, edgeSet );

        logger->update_progress( 30.0f );

        tbb::parallel_sort( edgeSet.begin(), edgeSet.end() );

        logger->update_progress( 60.0f );

        {
            frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, 60.f, 100.f );

            for( size_t i = 1; i < edgeSet.size(); ++i ) {
                if( edgeSet[i].m_edge == edgeSet[i - 1].m_edge ) {
                    out << "Duplicate edge (" << edgeSet[i].m_edge.first << "," << edgeSet[i].m_edge.second
                        << ") found bewteen faces " << edgeSet[i - 1].m_face << " and " << edgeSet[i].m_face << "."
                        << std::endl;
                    return false;
                }

                if( i % updateInterval == 0 ) {
                    logger->update_progress( float( i ) / 100.0f );
                }
            }
        }
    }

    size_t currentChannel = 1;

    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator ichan = m_namedVertexChannels.begin(),
                                                                             ichane = m_namedVertexChannels.end();
         ichan != ichane; ++ichan ) {

        const_trimesh3_vertex_channel_general_accessor channel = ichan->second.get_general_accessor( &m_faces );

        if( channel.has_custom_faces() ) {
            ++currentChannel;
            frantic::logging::progress_logger_subinterval_tracker subinterval(
                *logger, channelProgressStride * currentChannel, channelProgressStride * ( currentChannel + 1 ) );

            logger->update_progress( 0.0f );

            collect_edges_with_faces( *this, edgeSet );

            logger->update_progress( 30.0f );

            tbb::parallel_sort( edgeSet.begin(), edgeSet.end() );

            logger->update_progress( 60.0f );

            size_t duplicateCount = 0;

            {
                frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, 60.f, 100.f );

                for( size_t i = 1; i < edgeSet.size(); ++i ) {
                    if( edgeSet[i].m_edge == edgeSet[i - 1].m_edge ) {
                        ++duplicateCount;
                    } else {
                        duplicateCount = 1;
                    }

                    if( duplicateCount >= 2 ) {
                        out << "In channel \"" << frantic::strings::to_string( ichan->first ) << "\" : duplicate edge ("
                            << edgeSet[i].m_edge.first << "," << edgeSet[i].m_edge.second << ") found bewteen faces "
                            << edgeSet[i - 1].m_face << " and " << edgeSet[i].m_face << "." << std::endl;
                        return false;
                    }

                    if( i % updateInterval == 0 ) {
                        logger->update_progress( float( i ) / 100.0f );
                    }
                }
            }
        }
    }

    return true;
}

bool trimesh3::check_finite_vertices( std::ostream& out, frantic::logging::progress_logger* logger ) const {

    frantic::logging::null_progress_logger nullLogger;
    if( !logger ) {
        logger = &nullLogger;
    }

    logger->set_title( _T( "Checking for infinite vertices" ) );

    const size_t vertexUpdateInterval = ( vertex_count() / 20 ) + 1;

    // Check that all the vertices are finite
    for( size_t i = 0, ie = m_vertices.size(); i != ie; ++i ) {
        if( !m_vertices[i].is_finite() ) {
            out << "Vertex " << i << " has a non-finite coordinate value." << endl;
            return false;
        }

        if( i % vertexUpdateInterval == 0 ) {
            logger->update_progress( ( 100.0f * i ) / vertex_count() );
        }
    }

    return true;
}

bool trimesh3::check_index_ranges( std::ostream& out, frantic::logging::progress_logger* logger ) const {
    frantic::logging::null_progress_logger nullLogger;

    if( !logger ) {
        logger = &nullLogger;
    }

    logger->set_title( _T( "Checking mesh index ranges" ) );

    const size_t faceUpdateInterval = ( face_count() / 20 ) + 1;
    const size_t customFaceChannelCount = count_named_vertex_channels_with_custom_faces() + 1;
    const float channelProgressStride = 100.0f / customFaceChannelCount;

    {
        frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, 0.0f, channelProgressStride );

        // Check that all the indexes in the faces are within range
        for( size_t i = 0, ie = m_faces.size(); i != ie; ++i ) {
            vector3 face = m_faces[i];
            for( int j = 0; j < 3; ++j ) {
                if( (unsigned)face[j] >= m_vertices.size() ) {
                    out << "Face " << i << " has an invalid vertex index " << face[j] << endl;
                    return false;
                }
            }

            if( i % faceUpdateInterval == 0 ) {
                logger->update_progress( ( i * 100.0f ) / face_count() );
            }
        }
    }

    size_t currentChannel = 1;

    // Do the same checks for all the vertex channels
    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator ichan = m_namedVertexChannels.begin(),
                                                                             ichane = m_namedVertexChannels.end();
         ichan != ichane; ++ichan ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, channelProgressStride * currentChannel, channelProgressStride * ( currentChannel + 1 ) );
        ++currentChannel;

        const_trimesh3_vertex_channel_general_accessor channel = ichan->second.get_general_accessor( &m_faces );

        size_t size = channel.size();
        if( channel.has_custom_faces() ) {
            // Check that all the indexes in the faces are within range
            for( size_t i = 0, ie = channel.face_count(); i != ie; ++i ) {
                vector3 face = channel.face( i );
                for( int j = 0; j < 3; ++j ) {
                    if( (unsigned)face[j] >= size ) {
                        out << "In channel \"" << frantic::strings::to_string( ichan->first ) << "\", face " << i
                            << " has an invalid vertex index " << face[j] << endl;
                        return false;
                    }
                }

                if( i % faceUpdateInterval == 0 ) {
                    logger->update_progress( ( i * 100.0f ) / face_count() );
                }
            }
        } else if( size != m_vertices.size() ) {
            out << "Channel \"" << frantic::strings::to_string( ichan->first )
                << "\" does not have custom faces, but its data element count, " << size
                << ", doesn't match the mesh vertex count, " << m_vertices.size() << "." << endl;
            return false;
        }
    }

    return true;
}

bool trimesh3::check_zero_area_faces( std::ostream& out, frantic::logging::progress_logger* logger,
                                      float epsilon ) const {
    frantic::logging::null_progress_logger nullLogger;

    if( epsilon < 0.0f ) {
        throw std::runtime_error( "trimesh3::check_zero_area_faces: Error, epsilon must be a positive number." );
    }

    if( !logger ) {
        logger = &nullLogger;
    }

    logger->set_title( _T( "Checking zero area faces" ) );

    const size_t faceUpdateInterval = ( face_count() / 20 ) + 1;

    for( size_t faceId = 0; faceId < face_count(); ++faceId ) {
        const vector3& face = get_face( faceId );

        const vector3f& a = get_vertex( face[0] );
        const vector3f& b = get_vertex( face[1] );
        const vector3f& c = get_vertex( face[2] );

        float area = frantic::graphics::triangle_area( a, b, c );

        if( area < epsilon ) {
            out << "Face " << faceId << " had zero area." << std::endl;
            return false;
        }

        if( faceId % faceUpdateInterval == 0 ) {
            logger->update_progress( ( faceId * 100.0f ) / face_count() );
        }
    }

    return true;
}

bool trimesh3::check_duplicate_face_indices( std::ostream& out, frantic::logging::progress_logger* logger ) const {
    frantic::logging::null_progress_logger nullLogger;

    if( !logger ) {
        logger = &nullLogger;
    }

    logger->set_title( _T( "Checking duplicate face indices" ) );

    const size_t faceUpdateInterval = ( face_count() / 20 ) + 1;
    const size_t customFaceChannelCount = count_named_vertex_channels_with_custom_faces() + 1;
    const float channelProgressStride = 100.0f / customFaceChannelCount;

    {
        frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, 0.0f, channelProgressStride );

        for( size_t faceId = 0; faceId < face_count(); ++faceId ) {
            const vector3& face = get_face( faceId );

            for( size_t vertexIndex = 0; vertexIndex < 3; ++vertexIndex ) {
                size_t currV = face[vertexIndex];
                size_t nextV = face[( vertexIndex + 1 ) % 3];

                if( currV == nextV ) {
                    out << "Duplicate vertex index " << currV << " detected on face " << faceId << std::endl;
                    return false;
                }
            }

            if( faceId % faceUpdateInterval == 0 ) {
                logger->update_progress( ( faceId * 100.0f ) / face_count() );
            }
        }
    }

    size_t currentChannel = 1;

    // Do the same checks for all the vertex channels
    for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator ichan = m_namedVertexChannels.begin(),
                                                                             ichane = m_namedVertexChannels.end();
         ichan != ichane; ++ichan ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, channelProgressStride * currentChannel, channelProgressStride * ( currentChannel + 1 ) );
        ++currentChannel;

        const_trimesh3_vertex_channel_general_accessor channel = ichan->second.get_general_accessor( &m_faces );

        if( channel.has_custom_faces() ) {
            // Check that all the indexes in the faces are within range
            for( size_t faceId = 0; faceId < channel.face_count(); ++faceId ) {
                vector3 face = channel.face( faceId );
                for( size_t vertexIndex = 0; vertexIndex < 3; ++vertexIndex ) {
                    size_t currV = face[vertexIndex];
                    size_t nextV = face[( vertexIndex + 1 ) % 3];

                    if( currV == nextV ) {
                        out << "In channel \"" << frantic::strings::to_string( ichan->first )
                            << "\", duplicate vertex index " << currV << " detected on face " << faceId << std::endl;
                        return false;
                    }
                }

                if( faceId % faceUpdateInterval == 0 ) {
                    logger->update_progress( ( faceId * 100.0f ) / face_count() );
                }
            }
        }
    }

    return true;
}

// TODO: move these checks outside of trimesh3, and parallelize them
bool trimesh3::check_consistency( std::ostream& out, consistency_check_flag flags,
                                  frantic::logging::progress_logger* logger ) const {
    frantic::logging::null_progress_logger nullLogger;

    if( !logger ) {
        logger = &nullLogger;
    }

    const size_t totalCount = ( flags & bad_channel_sizes ? 1 : 0 ) + ( flags & out_of_range_indices ? 1 : 0 ) +
                              ( flags & infinite_vertices ? 1 : 0 ) + ( flags & duplicate_edges ? 1 : 0 ) +
                              ( flags & duplicate_face_indices ? 1 : 0 ) + ( flags & zero_area_faces ? 1 : 0 );
    const float progressIntervalStride = 100.0f / totalCount;

    size_t currentCheck = 0;

    if( flags & bad_channel_sizes ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, progressIntervalStride * currentCheck, progressIntervalStride * ( currentCheck + 1 ) );
        ++currentCheck;
        if( !check_channel_sizes( out, logger ) ) {
            return false;
        }
    }

    if( flags & out_of_range_indices ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, progressIntervalStride * currentCheck, progressIntervalStride * ( currentCheck + 1 ) );
        ++currentCheck;
        if( !check_index_ranges( out, logger ) ) {
            return false;
        }
    }

    if( flags & infinite_vertices ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, progressIntervalStride * currentCheck, progressIntervalStride * ( currentCheck + 1 ) );
        ++currentCheck;
        if( !check_finite_vertices( out, logger ) ) {
            return false;
        }
    }

    if( flags & duplicate_face_indices ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, progressIntervalStride * currentCheck, progressIntervalStride * ( currentCheck + 1 ) );
        ++currentCheck;
        if( !check_duplicate_face_indices( out, logger ) ) {
            return false;
        }
    }

    if( flags & duplicate_edges ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, progressIntervalStride * currentCheck, progressIntervalStride * ( currentCheck + 1 ) );
        ++currentCheck;
        if( !check_duplicate_edges( out, logger ) ) {
            return false;
        }
    }

    if( flags & zero_area_faces ) {
        frantic::logging::progress_logger_subinterval_tracker subinterval(
            *logger, progressIntervalStride * currentCheck, progressIntervalStride * ( currentCheck + 1 ) );
        ++currentCheck;
        if( !check_zero_area_faces( out, logger ) ) {
            return false;
        }
    }

    return true;
}

namespace {

void copy_vertices( frantic::geometry::trimesh3& outMesh, const frantic::geometry::mesh_interface* inMesh ) {
    assert( inMesh );

    const std::size_t vertexCount = inMesh->get_num_verts();

    outMesh.set_vertex_count( vertexCount );

    for( std::size_t i = 0; i < vertexCount; ++i ) {
        float buffer[3];
        inMesh->get_vert( i, buffer );
        outMesh.get_vertex( i ).set( buffer[0], buffer[1], buffer[2] );
    }
}

frantic::graphics::vector3 to_vector3( std::size_t indexBuffer[3] ) {
    frantic::graphics::vector3 result;
    for( int i = 0; i < 3; ++i ) {
        try {
            result[i] = boost::numeric_cast<boost::int32_t>( indexBuffer[i] );
        } catch( boost::bad_numeric_cast& e ) {
            throw std::runtime_error( "to_vector3 Error: " + std::string( e.what() ) );
        }
    }
    return result;
}

frantic::graphics::vector3 get_face( const frantic::geometry::mesh_interface* mesh, const std::size_t faceIndex ) {
    assert( mesh );

    std::size_t indexBuffer[3];
    const std::size_t cornerCount = mesh->get_num_face_verts( faceIndex );
    if( cornerCount != 3 ) {
        throw std::runtime_error( "copy_faces Error: Can only copy triangular faces, but attempted to copy face with " +
                                  boost::lexical_cast<std::string>( cornerCount ) + " corners." );
    }

    mesh->get_face_vert_indices( faceIndex, indexBuffer );

    return to_vector3( indexBuffer );
}

std::size_t get_triangulated_face_count( const frantic::geometry::mesh_interface* mesh ) {
    std::size_t total = 0;
    const std::size_t totalFaces = mesh->get_num_faces();
    for( std::size_t i = 0; i < totalFaces; ++i ) {
        const std::size_t cornerCount = mesh->get_num_face_verts( i );
        if( cornerCount > 2 ) {
            total += cornerCount - 2;
        }
    }
    return total;
}

void get_triangulated_faces( const frantic::geometry::mesh_interface* mesh, const std::size_t faceIndex,
                             std::vector<frantic::graphics::vector3>& outCornerIndices,
                             std::vector<frantic::graphics::vector3>& outFaces ) {
    assert( mesh );
    const std::size_t cornerCount = mesh->get_num_face_verts( faceIndex );

    if( cornerCount <= 2 ) {
        return;
    }

    // Get the vertex indices of the polygon
    std::vector<std::size_t> polygonIndices( cornerCount );
    mesh->get_face_vert_indices( faceIndex, &polygonIndices[0] );

    // Get the actual polygon vertex positions
    std::vector<vector3f> polygonVertices( cornerCount );
    for( std::size_t i = 0; i < cornerCount; ++i ) {
        polygonVertices[i] = mesh->get_vert( polygonIndices[i] );
    }

    // Triangulate the polygon
    outFaces.resize( cornerCount - 2 );
    outCornerIndices.reserve( cornerCount - 2 );

    const char* error = NULL;
    const bool ok = triangulate_polygon_robust( &polygonVertices[0], cornerCount, outCornerIndices, error );
    if( !ok ) {
        // "Infinite Loop Detected" error.  Try to ignore it and move on
        FF_LOG( warning ) << "Error triangulating polymesh: " << error << std::endl;
    }

    // Set the proper final faceindices
    for( std::size_t i = 0; i < cornerCount - 2; ++i ) {
        outFaces[i].x = static_cast<int>( polygonIndices[outCornerIndices[i].x] );
        outFaces[i].y = static_cast<int>( polygonIndices[outCornerIndices[i].y] );
        outFaces[i].z = static_cast<int>( polygonIndices[outCornerIndices[i].z] );
    }
}

void copy_faces( frantic::geometry::trimesh3& outMesh, const frantic::geometry::mesh_interface* inMesh ) {
    assert( inMesh );

    const std::size_t faceCount = inMesh->get_num_faces();
    outMesh.set_face_count( faceCount );

    for( std::size_t i = 0; i < faceCount; ++i ) {
        outMesh.get_face( i ) = get_face( inMesh, i );
    }
}

void copy_triangulated_faces( frantic::geometry::trimesh3& outMesh, const frantic::geometry::mesh_interface* inMesh,
                              std::vector<std::vector<vector3>>& outCorners ) {
    assert( inMesh );

    const std::size_t faceCount = inMesh->get_num_faces();
    outCorners.resize( faceCount );

    const std::size_t finalFaceCount = get_triangulated_face_count( inMesh );
    outMesh.set_face_count( finalFaceCount );

    std::size_t n = 0;
    std::vector<vector3> currentFace;
    for( std::size_t i = 0; i < faceCount; ++i ) {
        std::vector<vector3>& currentCorner = outCorners[i];
        get_triangulated_faces( inMesh, i, currentCorner, currentFace );
        for( std::size_t j = 0; j < currentFace.size(); ++j ) {
            outMesh.get_face( n ) = currentFace[j];
            ++n;
        }
    }
}

frantic::graphics::vector3 get_face( const frantic::geometry::mesh_channel* channel, const std::size_t faceIndex ) {
    assert( channel );

    const std::size_t cornerCount = channel->get_num_face_verts( faceIndex );
    if( cornerCount != 3 ) {
        throw std::runtime_error( "get_face Error: expected a triangle, but got a face with " +
                                  boost::lexical_cast<std::string>( cornerCount ) + " corners instead" );
    }

    std::size_t indexBuffer[3];
    for( std::size_t cornerIndex = 0; cornerIndex < 3; ++cornerIndex ) {
        indexBuffer[cornerIndex] = channel->get_fv_index( faceIndex, cornerIndex );
    }

    return to_vector3( indexBuffer );
}

void get_triangulated_faces( const frantic::geometry::mesh_channel* channel, const std::size_t faceIndex,
                             const std::vector<vector3>& readCorners, std::vector<vector3>& outFaces ) {
    assert( channel );

    const std::size_t cornerCount = channel->get_num_face_verts( faceIndex );
    std::vector<std::size_t> indexBuffer( cornerCount );
    for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
        indexBuffer[cornerIndex] = channel->get_fv_index( faceIndex, cornerIndex );
    }

    outFaces.resize( readCorners.size() );
    for( std::size_t i = 0; i < readCorners.size(); ++i ) {
        std::size_t maxCorner = std::max( readCorners[i].x, std::max( readCorners[i].y, readCorners[i].z ) );
        if( maxCorner >= cornerCount ) {
            throw std::runtime_error( "get_triangulated_faces requested a corner index of " +
                                      boost::lexical_cast<std::string>( maxCorner ) + " when there are " +
                                      boost::lexical_cast<std::string>( cornerCount ) + " corners" );
        }
        outFaces[i].x = static_cast<int>( indexBuffer[readCorners[i].x] );
        outFaces[i].y = static_cast<int>( indexBuffer[readCorners[i].y] );
        outFaces[i].z = static_cast<int>( indexBuffer[readCorners[i].z] );
    }
}

bool has_custom_faces( const frantic::geometry::mesh_channel* channel ) {
    if( channel && channel->get_channel_type() == frantic::geometry::mesh_channel::face_vertex ) {
        return true;
    } else {
        return false;
    }
}

template <class AccessorType>
void copy_channel_data( AccessorType& outAccessor, const frantic::geometry::mesh_channel* channel ) {
    assert( channel );

    const std::size_t elementSize = channel->get_element_size();
    const std::size_t elementCount = channel->get_num_elements();

    if( elementSize == 0 ) {
        throw std::runtime_error( "copy_channel_data Error: element size is 0" );
    }
    if( elementCount != outAccessor.size() ) {
        throw std::runtime_error( "copy_channel_data Error: mismatch between input and output channel size" );
    }
    if( elementSize != outAccessor.primitive_size() ) {
        throw std::runtime_error( "copy_channel_data Error: mismatch between input and output element size" );
    }

    boost::scoped_array<char> elementBuffer( new char[elementSize] );
    for( std::size_t element = 0; element < elementCount; ++element ) {
        channel->get_value( element, elementBuffer.get() );
        memcpy( outAccessor.data( element ), elementBuffer.get(), elementSize );
    }
}

void copy_channel_faces( frantic::geometry::trimesh3_vertex_channel_general_accessor& outAccessor,
                         const frantic::geometry::mesh_channel* channel ) {
    assert( channel );

    const std::size_t faceCount = channel->get_num_faces();

    if( faceCount != outAccessor.face_count() ) {
        throw std::runtime_error( "copy_channel_faces Error: mismatch between input and output face count" );
    }

    for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        outAccessor.face( faceIndex ) = get_face( channel, faceIndex );
    }
}

void copy_triangulated_channel_faces( frantic::geometry::trimesh3_vertex_channel_general_accessor& outAccessor,
                                      const frantic::geometry::mesh_channel* channel,
                                      const std::vector<std::vector<vector3>>& faceCornerMap ) {
    assert( channel );

    const std::size_t faceCount = channel->get_num_faces();

    std::size_t n = 0;
    std::vector<vector3> faces;
    for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        get_triangulated_faces( channel, faceIndex, faceCornerMap[faceIndex], faces );
        for( std::size_t outFaceIndex = 0; outFaceIndex < faces.size(); ++outFaceIndex ) {
            outAccessor.face( n ) = faces[outFaceIndex];
            ++n;
        }
    }
}

template <class AccessorType>
void copy_triangulated_channel_data( AccessorType& outAccessor, const frantic::geometry::mesh_channel* channel,
                                     const std::vector<std::vector<vector3>>& indexRemap ) {
    assert( channel );

    const std::size_t elementSize = channel->get_element_size();
    const std::size_t elementCount = channel->get_num_elements();

    if( elementSize == 0 ) {
        throw std::runtime_error( "copy_triangulated_channel_data Error: element size is 0" );
    }
    if( elementSize != outAccessor.primitive_size() ) {
        throw std::runtime_error(
            "copy_triangulated_channel_data Error: mismatch between input and output element size" );
    }

    boost::scoped_array<char> elementBuffer( new char[elementSize] );

    // Polymesh has been broken up into triangles, repeat the face data for each subtriangle
    std::size_t n = 0;
    for( std::size_t element = 0; element < elementCount; ++element ) {
        channel->get_value( element, elementBuffer.get() );
        for( std::size_t faceParts = 0; faceParts < indexRemap[element].size(); ++faceParts ) {
            if( n >= outAccessor.size() ) {
                throw std::runtime_error(
                    "copy_triangulated_channel_data Error: mismatch between input and output channel size" );
            }
            memcpy( outAccessor.data( n ), elementBuffer.get(), elementSize );
            ++n;
        }
    }
}

void copy_vertex_channels( frantic::geometry::trimesh3& outMesh, const frantic::geometry::mesh_interface* inMesh ) {
    assert( inMesh );

    using frantic::geometry::mesh_interface;

    const mesh_interface::mesh_channel_map& vertexChannels = inMesh->get_vertex_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = vertexChannels.begin(), ie = vertexChannels.end();
         i != ie; ++i ) {
        const frantic::geometry::mesh_channel* channel = i->second;
        if( !channel ) {
            throw std::runtime_error( "copy_vertex_channels Error: channel is NULL" );
        }

        const bool hasCustomFaces = has_custom_faces( channel );
        outMesh.add_vertex_channel_raw( channel->get_name(), channel->get_data_arity(), channel->get_data_type(),
                                        channel->get_num_elements(), hasCustomFaces );

        frantic::geometry::trimesh3_vertex_channel_general_accessor acc =
            outMesh.get_vertex_channel_general_accessor( channel->get_name() );
        copy_channel_data( acc, channel );

        if( hasCustomFaces ) {
            copy_channel_faces( acc, channel );
        }
    }
}

void copy_triangulated_vertex_channels( frantic::geometry::trimesh3& outMesh,
                                        const frantic::geometry::mesh_interface* inMesh,
                                        const std::vector<std::vector<vector3>>& faceCornerMap ) {
    assert( inMesh );

    using frantic::geometry::mesh_interface;

    const mesh_interface::mesh_channel_map& vertexChannels = inMesh->get_vertex_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = vertexChannels.begin(), ie = vertexChannels.end();
         i != ie; ++i ) {
        const frantic::geometry::mesh_channel* channel = i->second;
        if( !channel ) {
            throw std::runtime_error( "copy_vertex_channels Error: channel is NULL" );
        }

        const bool hasCustomFaces = has_custom_faces( channel );
        outMesh.add_vertex_channel_raw( channel->get_name(), channel->get_data_arity(), channel->get_data_type(),
                                        channel->get_num_elements(), hasCustomFaces );

        frantic::geometry::trimesh3_vertex_channel_general_accessor acc =
            outMesh.get_vertex_channel_general_accessor( channel->get_name() );
        copy_channel_data( acc, channel );

        if( hasCustomFaces ) {
            copy_triangulated_channel_faces( acc, channel, faceCornerMap );
        }
    }
}

void copy_face_channels( frantic::geometry::trimesh3& outMesh, const frantic::geometry::mesh_interface* inMesh ) {
    assert( inMesh );

    using frantic::geometry::mesh_interface;

    const mesh_interface::mesh_channel_map& faceChannels = inMesh->get_face_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = faceChannels.begin(), ie = faceChannels.end(); i != ie;
         ++i ) {
        const frantic::geometry::mesh_channel* channel = i->second;
        if( !channel ) {
            throw std::runtime_error( "copy_face_channels Error: channel is NULL" );
        }

        if( channel->get_num_elements() != outMesh.face_count() ) {
            throw std::runtime_error( "copy_face_channels Error: mismatch between number of faces in geometry (" +
                                      boost::lexical_cast<std::string>( outMesh.face_count() ) +
                                      ") and elements in face channel \"" +
                                      frantic::strings::to_string( channel->get_name() ) + "\" (" +
                                      boost::lexical_cast<std::string>( channel->get_num_elements() ) + ")" );
        }

        outMesh.add_face_channel_raw( channel->get_name(), channel->get_data_arity(), channel->get_data_type() );

        frantic::geometry::trimesh3_face_channel_general_accessor acc =
            outMesh.get_face_channel_general_accessor( channel->get_name() );
        copy_channel_data( acc, channel );
    }
}

void copy_triangulated_face_channels( frantic::geometry::trimesh3& outMesh,
                                      const frantic::geometry::mesh_interface* inMesh,
                                      const std::vector<std::vector<vector3>>& faceCornerMap ) {
    assert( inMesh );

    using frantic::geometry::mesh_interface;

    const mesh_interface::mesh_channel_map& faceChannels = inMesh->get_face_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = faceChannels.begin(), ie = faceChannels.end(); i != ie;
         ++i ) {
        const frantic::geometry::mesh_channel* channel = i->second;
        if( !channel ) {
            throw std::runtime_error( "copy_triangulated_face_channels Error: channel is NULL" );
        }

        outMesh.add_face_channel_raw( channel->get_name(), channel->get_data_arity(), channel->get_data_type() );

        frantic::geometry::trimesh3_face_channel_general_accessor acc =
            outMesh.get_face_channel_general_accessor( channel->get_name() );
        copy_triangulated_channel_data( acc, channel, faceCornerMap );
    }
}

} // anonymous namespace

void copy_to_trimesh3( const frantic::geometry::mesh_interface* meshInterface, frantic::geometry::trimesh3& outMesh,
                       bool autoTriangulate ) {
    if( !meshInterface ) {
        throw std::runtime_error( "copy_to_trimesh3 Error: meshInterface is NULL" );
    }

    frantic::geometry::trimesh3 result;

    copy_vertices( result, meshInterface );
    if( !autoTriangulate || is_triangle_mesh( meshInterface ) ) {
        copy_faces( result, meshInterface );

        copy_vertex_channels( result, meshInterface );
        copy_face_channels( result, meshInterface );
    } else {
        std::vector<std::vector<vector3>> corners;
        copy_triangulated_faces( result, meshInterface, corners );

        copy_triangulated_vertex_channels( result, meshInterface, corners );
        copy_triangulated_face_channels( result, meshInterface, corners );
    }

    outMesh.swap( result );
}
} // namespace geometry
} // namespace frantic
