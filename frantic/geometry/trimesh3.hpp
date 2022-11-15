// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include <boost/move/move.hpp>

#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/plane3f.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/boundrect2f.hpp>

#include <frantic/geometry/trimesh3_named_channels.hpp>

#include <frantic/logging/progress_logger.hpp>

#include <frantic/misc/utility.hpp>

namespace frantic {
namespace geometry {
using ::frantic::graphics::boundbox3f;
using ::frantic::graphics::plane3f;
using ::frantic::graphics::transform4f;
using ::frantic::graphics::vector3;
using ::frantic::graphics::vector3f;

// forward declaration
class mesh_interface;

class trimesh3 {
  public:
    typedef vector3::value_type index_t;

  private:
    // Mark this class as copyable and movable.
    BOOST_COPYABLE_AND_MOVABLE( trimesh3 )

    // The mesh starts with vertices and faces
    std::vector<vector3f> m_vertices;
    std::vector<vector3> m_faces;

    // The set of named vertex channels, supporting additional per-vertex data.
    std::map<frantic::tstring, trimesh3_vertex_channel> m_namedVertexChannels;

    // The set of named face channels, supporting additional per-face data.
    std::map<frantic::tstring, trimesh3_face_channel> m_namedFaceChannels;

    // Helpers to assist with the faces
    // TODO: Probably remove these members?
    std::vector<plane3f> m_facePlanes;
    std::vector<boundbox3f> m_faceBoundBoxes;

    // Helpers to assist with computing barycentric coordinates on a face.  These have the same
    // indices as the faces
    // TODO: Probably remove these members?
    std::vector<char> m_barycentric0Axis;
    std::vector<char> m_barycentric1Axis;
    std::vector<float> m_barycentricInverseDeterminant;

  public:
    enum vector_type { VECTOR, NORMAL };

    enum consistency_check_flag {
        no_checks = 0,
        infinite_vertices = 1 << 0,
        out_of_range_indices = 1 << 1,
        duplicate_edges = 1 << 2,
        duplicate_face_indices = 1 << 3,
        zero_area_faces = 1 << 4,
        bad_channel_sizes = 1 << 5,

        all_checks = 0xFFFFFFFF,
    };

    ///////////////
    // Constructors
    ///////////////

    trimesh3() {}

    trimesh3( const trimesh3& rhs );

    ~trimesh3() {}

    trimesh3( BOOST_RV_REF( trimesh3 ) rhs ) { swap( rhs ); }

    void swap( trimesh3& rhs );

    /**
     * This function sets the trimesh3 to a box.
     *
     * @param  box  The bounding box to use for creating the box.
     */
    void set_to_box( const boundbox3f& box = boundbox3f( 0, 1, 0, 1, 0, 1 ) );

    /**
     * This function sets the trimesh3 to a box, with disjoint faces. So there are 24 vertices and 12 faces.
     * The verts are added so that each set of 4 verts correspond to face ordering
     * xneg, xpos, yneg, ypos, zneg, zpos
     *
     * @param  box  The bounding box to use for creating the box.
     */
    void set_to_exploded_box( const boundbox3f& box = boundbox3f( 0, 1, 0, 1, 0, 1 ) );

    void set_existing_exploded_box_verts( const boundbox3f& box );

    /**
     * After a box has been created, this function moves the vertices to reflect a new bounding box.
     *
     * @param  box  The bounding box to use for setting the vert positions.
     */
    void set_existing_box_verts( const boundbox3f& box );

    void set_to_tetrahedron( const vector3f& centre = vector3f( 0 ), float radius = 1.f );

    void set_to_icosahedron( const vector3f& centre = vector3f( 0 ), float radius = 1 );

    void set_to_rectangular_grid( int horizontalSegments, int verticalSegments,
                                  const ::frantic::graphics2d::boundrect2f& bounds );

    void set_to_hexagonal_grid( int horizontalSegments, int verticalSegments,
                                const ::frantic::graphics2d::boundrect2f& bounds );

    ////////////////////////////////////////
    // Queries
    ////////////////////////////////////////

    bool is_empty() const { return m_vertices.size() == 0 || m_faces.size() == 0; }

    void dump( std::ostream& out ) const;
    void dump_channel_data( std::ostream& out ) const;

    ////////////////////////////////////////
    // Vertex methods
    ////////////////////////////////////////

    std::vector<vector3f>& vertices_ref() { return m_vertices; }

    const std::vector<vector3f>& vertices_ref() const { return m_vertices; }

    std::size_t vertex_count() const { return m_vertices.size(); }

    void set_vertex_count( std::size_t vertexCount ) { m_vertices.resize( vertexCount ); }

    const vector3f& get_vertex( std::size_t v ) const { return m_vertices[v]; }

    vector3f& get_vertex( std::size_t v ) { return m_vertices[v]; }

    void add_vertex( const vector3f& vertex ) { m_vertices.push_back( vertex ); }

    void add_vertex( float x, float y, float z ) { m_vertices.push_back( vector3f( x, y, z ) ); }

    /**
     *  Add a vertex to the mesh, and resize the named channels to match
     * the new vertex count.
     *
     * @param v the vertex coordinate to add to the mesh.
     */
    void add_vertex_and_resize_channels( const vector3f& v );

    void add_vertices( const std::vector<vector3f>& vertices ) {

        if( vertices.empty() )
            return;

        size_t currentNumVerts = m_vertices.size();
        m_vertices.resize( currentNumVerts + vertices.size() );

        memcpy( &m_vertices[currentNumVerts], &vertices[0], vertices.size() * sizeof( vector3f ) );
    }

    /**
     *	Adds n default verts to the mesh.
     */
    void add_vertices( size_t n ) { m_vertices.resize( m_vertices.size() + n ); }

    /**
     *	Reserves n vertices
     */
    void reserve_vertices( size_t n ) { m_vertices.reserve( n ); }
    ////////////////////////////////////////
    // Face methods
    ////////////////////////////////////////

    std::vector<vector3>& faces_ref() { return m_faces; }

    const std::vector<vector3>& faces_ref() const { return m_faces; }

    std::size_t face_count() const { return m_faces.size(); }

    void set_face_count( std::size_t faceCount ) { m_faces.resize( faceCount ); }

    vector3 get_face( std::size_t f ) const { return m_faces[f]; }

    vector3& get_face( std::size_t f ) { return m_faces[f]; }

    // NOTE: When code adds a face here, it is its responsibility to also add the corresponding face to any
    //		 named vertex channels with custom faces and any named face channels.
    //       Failure to do so will result in an exception when the accessor for a channel is requested.

    void add_face( const vector3& face ) { m_faces.push_back( face ); }

    void add_face( int v0, int v1, int v2 ) { m_faces.push_back( vector3( v0, v1, v2 ) ); }

    // Adds an arbitrary face to the mesh by triangulating it.
    // NOTE: It's important that this triangulates the face *exactly* like the add_face function in
    // trimesh3_channels.hpp
    void add_face( const std::vector<int>& vertexIndices ) {
        if( vertexIndices.size() > 2 ) {
            int iForward = 0, iBackward = (int)vertexIndices.size() - 1;
            while( iBackward - iForward > 1 ) {
                // Alternate incrementing iForward and decrementing iBackward to
                // avoid having all triangles emanate from the same vertex.
                if( ( iBackward - iForward ) % 2 == 0 ) {
                    m_faces.push_back(
                        vector3( vertexIndices[iForward], vertexIndices[iForward + 1], vertexIndices[iBackward] ) );
                    iForward++;
                } else {
                    m_faces.push_back(
                        vector3( vertexIndices[iBackward - 1], vertexIndices[iBackward], vertexIndices[iForward] ) );
                    iBackward--;
                }
            }
        }
    }

    /**
     *  Add n faces to the mesh.
     */
    void add_faces( size_t n ) { m_faces.resize( m_faces.size() + n ); }

    /**
     *	Reserves n faces
     */
    void reserve_faces( size_t n ) { m_faces.reserve( n ); }

    /**
     *  Add a face to the mesh, and resize the named channels to match
     * the new face count.
     *
     * @param face the vertex numbers of the face to add to the mesh.
     */
    void add_face_and_resize_channels( const vector3& face );

    ////////////////////////////////////////
    // Named vertex channel methods
    ////////////////////////////////////////

    bool has_vertex_channel( const frantic::tstring& name ) const {
        return m_namedVertexChannels.find( name ) != m_namedVertexChannels.end();
    }

    void erase_vertex_channel( const frantic::tstring& name ) { m_namedVertexChannels.erase( name ); }

    // This retrieves the names of all the vertex channels
    void get_vertex_channel_names( std::vector<frantic::tstring>& outNames ) const {
        outNames.reserve( outNames.size() + m_namedVertexChannels.size() );
        for( std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.begin(),
                                                                                 ie = m_namedVertexChannels.end();
             i != ie; ++i ) {
            outNames.push_back( i->first );
        }
    }

    template <class DataType>
    void add_vertex_channel( const frantic::tstring& name ) {
        add_vertex_channel_raw( name, frantic::channels::channel_data_type_traits<DataType>::arity(),
                                frantic::channels::channel_data_type_traits<DataType>::data_type() );
    }

    // TODO this should be changed to accept a bool for is custom data or not. The custom face count should be exactly
    // the same as the base face count
    //  if there are custom faces, or it will be empty and nothing else.
    template <class DataType>
    void add_vertex_channel( const frantic::tstring& name, std::size_t vertexCount, bool hasCustomFaces = false ) {
        add_vertex_channel_raw( name, frantic::channels::channel_data_type_traits<DataType>::arity(),
                                frantic::channels::channel_data_type_traits<DataType>::data_type(), vertexCount,
                                hasCustomFaces );
    }

    void add_vertex_channel_raw( const frantic::tstring& name, std::size_t arity, data_type_t dataType ) {
        std::map<frantic::tstring, trimesh3_vertex_channel>::iterator i = m_namedVertexChannels.find( name );
        if( i == m_namedVertexChannels.end() ) {
            m_namedVertexChannels.insert( std::make_pair( name, trimesh3_vertex_channel( name, arity, dataType ) ) );
            i = m_namedVertexChannels.find( name );
            // Make sure the vertex array count matches that of the mesh's vertices
            i->second.m_data.resize( m_vertices.size() * i->second.m_primitiveSize );
            i->second.m_hasCustomFaces = false;
        } else {
            // TODO: Should this rather erase the existing channel and replace it with a new one?
            if( i->second.arity() != arity || i->second.data_type() != dataType )
                throw std::runtime_error(
                    "trimesh3.add_vertex_channel: Tried to add vertex channel \"" +
                    frantic::strings::to_string( name ) + "\", with data type " +
                    frantic::strings::to_string( channel_data_type_str( arity, dataType ) ) +
                    ".  This could not be done, because the channel already exists with data type " +
                    frantic::strings::to_string( channel_data_type_str( i->second.arity(), i->second.data_type() ) ) +
                    "." );
        }
    }

    void add_vertex_channel_raw( const frantic::tstring& name, std::size_t arity, data_type_t dataType,
                                 std::size_t vertexCount, bool hasCustomFaces = false ) {
        std::map<frantic::tstring, trimesh3_vertex_channel>::iterator i = m_namedVertexChannels.find( name );
        if( i == m_namedVertexChannels.end() ) {
            m_namedVertexChannels.insert( std::make_pair( name, trimesh3_vertex_channel( name, arity, dataType ) ) );
            i = m_namedVertexChannels.find( name );
            // Reset the vertices and faces arrays to the requested size
            i->second.m_data.resize( vertexCount * i->second.m_primitiveSize );
            i->second.m_hasCustomFaces = hasCustomFaces;
            if( hasCustomFaces ) {
                i->second.m_faces.resize( m_faces.size() );
            } else {
                i->second.m_faces.clear();
            }
        } else {
            // TODO: Should this rather erase the existing channel and replace it with a new one?
            if( i->second.arity() != arity || i->second.data_type() != dataType )
                throw std::runtime_error(
                    "trimesh3.add_vertex_channel: Tried to add vertex channel \"" +
                    frantic::strings::to_string( name ) + "\", with data type " +
                    frantic::strings::to_string( channel_data_type_str( arity, dataType ) ) +
                    ".  This could not be done, because the channel already exists with data type " +
                    frantic::strings::to_string( channel_data_type_str( i->second.arity(), i->second.data_type() ) ) +
                    "." );

            // Reset the vertices and faces arrays to the requested size
            i->second.m_data.resize( vertexCount * i->second.m_primitiveSize );
            i->second.m_faces.clear(); // Clear to ensure the faces are all zeroed to start
            i->second.m_hasCustomFaces = hasCustomFaces;
            if( hasCustomFaces ) {
                i->second.m_faces.resize( m_faces.size() );
            }
        }
    }

    // NOTE: This invalidates any accessors pointing to this channel.
    void set_vertex_channel_custom_faces( const frantic::tstring& name, bool hasCustomFaces ) {
        std::map<frantic::tstring, trimesh3_vertex_channel>::iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            if( hasCustomFaces ) {
                if( !i->second.m_hasCustomFaces ) {
                    // Start off the faces matching the mesh's faces
                    i->second.m_faces = m_faces;
                    i->second.m_hasCustomFaces = true;
                }
            } else {
                if( i->second.m_hasCustomFaces ) {
                    // Remove the existing face array
                    i->second.m_faces.clear();
                    // Make sure the vertex array count matches that of the mesh's vertices
                    if( i->second.m_data.size() < m_vertices.size() )
                        i->second.m_data.resize( m_vertices.size() );
                    i->second.m_hasCustomFaces = false;
                }
            }
        } else {
            throw std::runtime_error(
                "trimesh3.set_vertex_channel_face_support: Tried to modify the face support of vertex channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    bool has_vertex_channel_custom_faces( const frantic::tstring& name ) const {
        std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            return i->second.has_custom_faces();
        } else {
            return false;
        }
    }

    size_t count_named_vertex_channels_with_custom_faces() const;

    trimesh3_vertex_channel_general_accessor get_vertex_channel_general_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, trimesh3_vertex_channel>::iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            return i->second.get_general_accessor( &m_faces );
        } else {
            throw std::runtime_error(
                "trimesh3.get_vertex_channel_general_accessor: Tried to retrieve an accessor for vertex channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    const_trimesh3_vertex_channel_general_accessor
    get_vertex_channel_general_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            return i->second.get_general_accessor( &m_faces );
        } else {
            throw std::runtime_error(
                "trimesh3.get_vertex_channel_general_accessor: Tried to retrieve an accessor for vertex channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    trimesh3_vertex_channel_accessor<DataType> get_vertex_channel_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, trimesh3_vertex_channel>::iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            return i->second.get_accessor<DataType>( &m_faces );
        } else {
            throw std::runtime_error(
                "trimesh3.get_vertex_channel_accessor: Tried to retrieve an accessor for vertex channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    const_trimesh3_vertex_channel_accessor<DataType> get_vertex_channel_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            return i->second.get_accessor<DataType>( &m_faces );
        } else {
            throw std::runtime_error(
                "trimesh3.get_vertex_channel_accessor: Tried to retrieve an accessor for vertex channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    trimesh3_vertex_channel_cvt_accessor<DataType> get_vertex_channel_cvt_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, trimesh3_vertex_channel>::iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            return i->second.get_cvt_accessor<DataType>( &m_faces );
        } else {
            throw std::runtime_error(
                "trimesh3.get_vertex_channel_cvt_accessor: Tried to retrieve an accessor for vertex channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    const_trimesh3_vertex_channel_cvt_accessor<DataType>
    get_vertex_channel_cvt_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, trimesh3_vertex_channel>::const_iterator i = m_namedVertexChannels.find( name );
        if( i != m_namedVertexChannels.end() ) {
            return i->second.get_cvt_accessor<DataType>( &m_faces );
        } else {
            throw std::runtime_error(
                "trimesh3.get_vertex_channel_cvt_accessor: Tried to retrieve an accessor for vertex channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    ////////////////////////////////////////
    // Named face channel methods
    ////////////////////////////////////////

    bool has_face_channel( const frantic::tstring& name ) const {
        return m_namedFaceChannels.find( name ) != m_namedFaceChannels.end();
    }

    void erase_face_channel( const frantic::tstring& name ) { m_namedFaceChannels.erase( name ); }

    // This retrieves the names of all the face channels
    void get_face_channel_names( std::vector<frantic::tstring>& outNames ) const {
        outNames.reserve( outNames.size() + m_namedFaceChannels.size() );
        for( std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.begin(),
                                                                               ie = m_namedFaceChannels.end();
             i != ie; ++i ) {
            outNames.push_back( i->first );
        }
    }

    template <class DataType>
    void add_face_channel( const frantic::tstring& name ) {
        add_face_channel_raw( name, frantic::channels::channel_data_type_traits<DataType>::arity(),
                              frantic::channels::channel_data_type_traits<DataType>::data_type() );
    }

    void add_face_channel_raw( const frantic::tstring& name, std::size_t arity, data_type_t dataType ) {
        std::map<frantic::tstring, trimesh3_face_channel>::iterator i = m_namedFaceChannels.find( name );
        if( i == m_namedFaceChannels.end() ) {
            m_namedFaceChannels.insert( std::make_pair( name, trimesh3_face_channel( name, arity, dataType ) ) );
            i = m_namedFaceChannels.find( name );
            // Make sure the data array count matches that of the mesh's faces
            i->second.m_data.resize( m_faces.size() * i->second.m_primitiveSize );
        } else {
            // TODO: Should this rather erase the existing channel and replace it with a new one?
            if( i->second.arity() != arity || i->second.data_type() != dataType )
                throw std::runtime_error(
                    "trimesh3.add_face_channel: Tried to add face channel \"" + frantic::strings::to_string( name ) +
                    "\", with data type " + frantic::strings::to_string( channel_data_type_str( arity, dataType ) ) +
                    ".  This could not be done, because the face channel already exists with data type " +
                    frantic::strings::to_string( channel_data_type_str( i->second.arity(), i->second.data_type() ) ) +
                    "." );

            // We potentially need to resize the channel in order to match the geometry channel. Otherwise this mesh
            // will become invalid because a per-face channel will not have one element per face.
            i->second.m_data.resize( m_faces.size() * i->second.m_primitiveSize );
        }
    }

    trimesh3_face_channel_general_accessor get_face_channel_general_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, trimesh3_face_channel>::iterator i = m_namedFaceChannels.find( name );
        if( i != m_namedFaceChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "trimesh3.get_face_channel_general_accessor: Tried to retrieve an accessor for face channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    const_trimesh3_face_channel_general_accessor
    get_face_channel_general_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.find( name );
        if( i != m_namedFaceChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "trimesh3.get_face_channel_general_accessor: Tried to retrieve an accessor for face channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    trimesh3_face_channel_accessor<DataType> get_face_channel_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, trimesh3_face_channel>::iterator i = m_namedFaceChannels.find( name );
        if( i != m_namedFaceChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "trimesh3.get_face_channel_accessor: Tried to retrieve an accessor for face channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    const_trimesh3_face_channel_accessor<DataType> get_face_channel_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.find( name );
        if( i != m_namedFaceChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "trimesh3.get_face_channel_accessor: Tried to retrieve an accessor for face channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    trimesh3_face_channel_cvt_accessor<DataType> get_face_channel_cvt_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, trimesh3_face_channel>::iterator i = m_namedFaceChannels.find( name );
        if( i != m_namedFaceChannels.end() ) {
            return i->second.get_cvt_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "trimesh3.get_face_channel_cvt_accessor: Tried to retrieve an accessor for face channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    template <class DataType>
    const_trimesh3_face_channel_cvt_accessor<DataType>
    get_face_channel_cvt_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, trimesh3_face_channel>::const_iterator i = m_namedFaceChannels.find( name );
        if( i != m_namedFaceChannels.end() ) {
            return i->second.get_cvt_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "trimesh3.get_face_channel_cvt_accessor: Tried to retrieve an accessor for face channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this trimesh3." );
        }
    }

    ////////////////////////////////////////
    // The face planes can be cached, for instance to speed up
    // ray intersections.
    ////////////////////////////////////////

    bool has_face_planes() const { return m_facePlanes.size() == m_faces.size(); }

    const plane3f& get_face_plane( int f ) const { return m_facePlanes[f]; }

    void build_face_planes( bool rebuild = false ) {
        if( !has_face_planes() || rebuild ) {
            m_facePlanes.resize( m_faces.size() );
            for( unsigned i = 0; i < m_facePlanes.size(); ++i )
                m_facePlanes[i] = calculate_face_plane( i );
        }
    }

    ////////////////////////////////////////
    // The face bounding boxes can be cached, for instance to speed up
    // ray intersections.
    ////////////////////////////////////////

    bool has_face_bound_boxes() const { return m_faceBoundBoxes.size() == m_faces.size(); }

    const boundbox3f& get_face_bound_box( int f ) const { return m_faceBoundBoxes[f]; }

    void build_face_bound_boxes( bool rebuild = false ) {
        if( !has_face_bound_boxes() || rebuild ) {
            m_faceBoundBoxes.resize( m_faces.size() );
            for( unsigned i = 0; i < m_faceBoundBoxes.size(); ++i )
                calculate_face_bound_box( i, m_faceBoundBoxes[i] );
        }
    }

    ////////////////////////////////////////

    boundbox3f compute_bound_box() const {
        if( m_vertices.size() > 0 ) {
            boundbox3f result( m_vertices[0] );
            for( unsigned i = 1; i < m_vertices.size(); ++i ) {
                result += m_vertices[i];
            }
            return result;
        } else {
            return boundbox3f::from_empty();
        }
    }

    boundbox3f compute_bound_box( const transform4f& xform ) const {
        if( m_vertices.size() > 0 ) {
            boundbox3f result( xform * m_vertices[0] );
            for( unsigned i = 1; i < m_vertices.size(); ++i ) {
                result += xform * m_vertices[i];
            }
            return result;
        } else {
            return boundbox3f::from_empty();
        }
    }

    vector3f compute_face_normal( std::size_t faceIndex ) const {
        vector3 face = get_face( faceIndex );
        return triangle_normal( get_vertex( face.x ), get_vertex( face.y ), get_vertex( face.z ) );
    }

    ////////////////////////////////////////
    // Operators
    ////////////////////////////////////////

    trimesh3& operator=( BOOST_COPY_ASSIGN_REF( trimesh3 ) rhs ) {
        m_vertices = rhs.m_vertices;
        m_faces = rhs.m_faces;
        m_facePlanes = rhs.m_facePlanes;
        m_faceBoundBoxes = rhs.m_faceBoundBoxes;
        m_namedVertexChannels = rhs.m_namedVertexChannels;
        m_namedFaceChannels = rhs.m_namedFaceChannels;
        m_barycentric0Axis = rhs.m_barycentric0Axis;
        m_barycentric1Axis = rhs.m_barycentric1Axis;
        m_barycentricInverseDeterminant = rhs.m_barycentricInverseDeterminant;

        return *this;
    }

    trimesh3& operator=( BOOST_RV_REF( trimesh3 ) rhs ) {
        if( this != &rhs ) {
            clear_and_deallocate();
            swap( rhs );
        }

        return *this;
    }

    /**
     *  Compare meshes for equality.
     *
     *  This comparison examines the mesh geometry and channel data.
     *
     * @param  rhs  the mesh to compare this mesh to.
     * @return  true if the meshes are the same, including geometry and
     *		channel data, and false otherwise.
     */
    bool operator==( const trimesh3& rhs ) const;

    /**
     * This clears all the data in the trimesh3 but doesn't deallocate.
     *
     */
    void clear() {
        m_vertices.clear();
        m_faces.clear();
        m_facePlanes.clear();
        m_faceBoundBoxes.clear();
        m_barycentric0Axis.clear();
        m_barycentric1Axis.clear();
        m_barycentricInverseDeterminant.clear();
        m_namedVertexChannels.clear();
        m_namedFaceChannels.clear();
    }

    /**
     * This clears all the data in the trimesh3.  It also deallocates the memory it had allocated.
     *
     */
    void clear_and_deallocate() {
        // Do clear_with_swap's so that the memory gets deallocated
        frantic::clear_with_swap( m_vertices );
        frantic::clear_with_swap( m_faces );
        frantic::clear_with_swap( m_facePlanes );
        frantic::clear_with_swap( m_faceBoundBoxes );
        frantic::clear_with_swap( m_barycentric0Axis );
        frantic::clear_with_swap( m_barycentric1Axis );
        frantic::clear_with_swap( m_barycentricInverseDeterminant );
        m_namedVertexChannels.clear();
        m_namedFaceChannels.clear();
    }

    /**
     * This clears all the named vertex channels.
     */
    void clear_vertex_channels() { m_namedVertexChannels.clear(); }

    /**
     * This clears all the named face channels.
     */
    void clear_face_channels() { m_namedFaceChannels.clear(); }

    /**
     * This translates the trimesh3 in place, using the provided translation offset and its time derivative.  It
     * affects the vertex positions as well as the Velocity channel.
     * <p>
     * A possible future extension would be to add a offsetTimeSecondDerivative parameter, which would
     * affect the Acceleration channel together with the given factors.  Currently the Acceleration channel
     * is left untouched.
     *
     * @param  offset                The translation applied to the mesh
     * @param  offsetTimeDerivative  The time derivative of the translation, which affects the Velocity channel.
     */
    void translate( const vector3f& offset, const vector3f& offsetTimeDerivative = vector3f( 0 ) );

    /**
     * This scales the trimesh3 in place, using the provided scale factor and its time derivative.  It
     * affects the vertex positions as well as the Velocity channel.
     * <p>
     * A possible future extension would be to add a factorTimeSecondDerivative parameter, which would
     * affect the Acceleration channel together with the given factors.  Currently the Acceleration channel
     * is left untouched.
     *
     * @param  factor                The scaling factor applied to the mesh
     * @param  factorTimeDerivative  The time derivative of the scaling factor, which affects the Velocity channel.
     */
    void scale( const vector3f& factor, const vector3f& factorTimeDerivative = vector3f( 0 ) );

    /**
     * This scales the trimesh3 in place, using the provided scale factor and its time derivative.  It
     * affects the vertex positions as well as the Velocity channel.
     * <p>
     * It calls the vector3f version of the scale function.
     *
     * @param  factor                The scaling factor applied to the mesh
     * @param  factorTimeDerivative  The time derivative of the scaling factor, which affects the Velocity channel.
     */
    void scale( float factor, float factorTimeDerivative = 0 ) {
        scale( vector3f( factor ), vector3f( factorTimeDerivative ) );
    }

    /**
     * This transforms the trimesh3 in place, using the provided transform matrix and time-derivative transform matrix.
     * It affects the vertex positions, as well as the Velocity and Normal channels if they exist.
     * <p>
     * A possible future extension would be to add a xformTimeSecondDerivative parameter, and use it in concert
     * with the other transforms to affect the Acceleration channel.  Currently this function doesn't modify the
     * Acceleration channel at all.
     *
     * @param  xform                The transform matrix to apply to the trimesh3.
     * @param  xformTimeDerivative  The time derivative of xform, which affects the Velocity channel.
     */
    void transform( const transform4f& xform, const transform4f& xformTimeDerivative = transform4f::zero() );

    void get_adjacency_list( std::vector<std::set<int>>& outAdjList ) const;

    float get_minimum_edge_length() const;

    /**
     * This merges all the faces and vertices of the given mesh into the current mesh.
     * <p>
     * Channels which are newly created due to the combination have their values defaulted to 0.
     *
     * @param  mesh  The input mesh to combine.
     */
    void combine( const trimesh3& mesh );

    /**
     * This merges all the faces and vertices of the given mesh into the current mesh.  It applies the provided
     * transform to the mesh on the fly.
     * <p>
     * Channels which are newly created due to the combination have their values defaulted to 0.
     *
     * @param  xform  The transform to apply to the input mesh.
     * @param  mesh  The input mesh to combine.
     */
    void combine( const transform4f& xform, const trimesh3& mesh );

    /**
     * This merges all the faces and vertices of the given mesh into the current mesh.  It applies the provided
     * transform to the mesh on the fly.
     * <p>
     * Channels which are newly created due to the combination have their values defaulted to 0.
     * <p>
     * Any channels that should not be transformed should be exluded from the transformChannels list.
     *
     * @param  xform  The transform to apply to the input mesh.
     * @param  mesh  The input mesh to combine.
     * @param  transformChannels  A list of channels that should be affected by the transform, and the type of transform
     * to use.
     */
    void combine( const transform4f& xform, const trimesh3& mesh,
                  std::vector<std::pair<frantic::tstring, vector_type>> transformChannels );

    /**
     * This merges all the faces and vertices of the given mesh into the current mesh.  It applies the provided
     * transform to the mesh on the fly.  The transform time derivative is used to affect the Velocity channel.
     * <p>
     * Channels which are newly created due to the combination have their values defaulted to 0.
     *
     * @param  xform  The transform to apply to the input mesh.
     * @param  xformTimeDerivative  The time derivative of the transform.  This affects the Velocity channel.
     * @param  mesh  The input mesh to combine.
     */
    void combine( const transform4f& xform, const transform4f& xformTimeDerivative, const trimesh3& mesh );

    /**
     * This merges all the faces and vertices of the given mesh into the current mesh.  It applies the provided
     * transform to the mesh on the fly.  The transform time derivative is used to affect the Velocity channel.
     * <p>
     * Channels which are newly created due to the combination have their values defaulted to 0.
     * <p>
     * Any channels that should not be transformed should be exluded from the transformChannels list.
     *
     * @param  xform  The transform to apply to the input mesh.
     * @param  xformTimeDerivative  The time derivative of the transform.  This affects the Velocity channel.
     * @param  mesh  The input mesh to combine.
     * @param  transformChannels  A list of channels that should be affected by the xform, and the type of transform to
     * use.
     */
    void combine( const transform4f& xform, const transform4f& xformTimeDerivative, const trimesh3& mesh,
                  std::vector<std::pair<frantic::tstring, vector_type>> transformChannels );

    void calculate_face_bound_box( std::size_t faceIndex, boundbox3f& outResult ) const {
        vector3 face = m_faces[faceIndex];
        outResult.set_to_point( m_vertices[face.x] );
        outResult += m_vertices[face.y];
        outResult += m_vertices[face.z];
        // Expand the bounds a bit to account for numerical error
        outResult.expand_fractional( 0.0001f );
    }

    /**
     *  Performs a linear interpolation between two trimeshes, including all
     *  vertex channels that are defined in both meshes.  If a channel is defined in one
     *	mesh, but not the other, the channel is not included in the final mesh.  Both
     *  meshes have to have the same topology (vertex count and face count must be identical).
     *
     * @param  mesh1	the first mesh operand
     * @param  mesh2	the second mesh operand
     * @param  alpha	the fractional distance between the two sample at which to interpolate the final sample
     */
    void linear_interpolation( const trimesh3& mesh1, const trimesh3& mesh2, const float& alpha );

    /**
     * Offsets a trimesh3, using a time offset in seconds to move the
     * vertices based on the Velocity channel. This is done by modifying the trimesh3 in place.
     *
     * @param  timeOffset  The time offset, used to move the vertices relative to their current
     *                     position, based on the Velocity channel of the input mesh.
     */
    void velocity_offset( float timeOffset );

    plane3f calculate_face_plane( std::size_t faceIndex ) const {
        vector3 face = m_faces[faceIndex];
        return plane3f::from_triangle( m_vertices[face.x], m_vertices[face.y], m_vertices[face.z] );
    }

    vector3f calculate_face_normal( std::size_t faceIndex ) const {
        vector3 face = m_faces[faceIndex];
        return triangle_normal( m_vertices[face.x], m_vertices[face.y], m_vertices[face.z] );
    }

    vector3f calculate_face_center( std::size_t faceIndex ) const {
        vector3 face = m_faces[faceIndex];
        return ( m_vertices[face.x] + m_vertices[face.y] + m_vertices[face.z] ) / 3.f;
    }

    bool has_barycentric_helpers() const {
        return m_barycentric0Axis.size() == m_faces.size() && m_barycentric1Axis.size() == m_faces.size() &&
               m_barycentricInverseDeterminant.size() == m_faces.size();
    }

    void build_barycentric_helpers( bool rebuild = false ) {
        if( !has_barycentric_helpers() || rebuild ) {
            m_barycentric0Axis.resize( m_faces.size() );
            m_barycentric1Axis.resize( m_faces.size() );
            m_barycentricInverseDeterminant.resize( m_faces.size() );

            for( unsigned i = 0; i < m_faces.size(); ++i ) {
                vector3 face = m_faces[i];

                compute_barycentric_helpers( m_vertices[face.x], m_vertices[face.y], m_vertices[face.z],
                                             m_barycentric0Axis[i], m_barycentric1Axis[i],
                                             m_barycentricInverseDeterminant[i] );
            }
        }
    }

    vector3f compute_barycentric_coordinates( int faceIndex, const vector3f& position ) const {
        vector3 face = m_faces[faceIndex];
        return compute_barycentric_coordinates_with_helpers(
            position, m_vertices[face.x], m_vertices[face.y], m_vertices[face.z], m_barycentric0Axis[faceIndex],
            m_barycentric1Axis[faceIndex], m_barycentricInverseDeterminant[faceIndex] );
    }

    vector3f get_barycentric_vertex( int faceIndex, const vector3f& barycentricCoordinates ) const {
        vector3 face = m_faces[faceIndex];
        return barycentricCoordinates.x * m_vertices[face.x] + barycentricCoordinates.y * m_vertices[face.y] +
               barycentricCoordinates.z * m_vertices[face.z];
    }

    //////////////////////////////////////////////////////
    // Miscellaneous operations

    void reverse_face_winding() {
        using std::swap;
        for( unsigned i = 0; i < m_faces.size(); ++i ) {
            swap( m_faces[i].x, m_faces[i].y );
        }
    }

    // Builds the vertex normals based on incident angle.
    void build_vertex_normals( const frantic::tstring& vertexChannelName = _T("Normal"), bool normalizeNormals = true );

    /**
     *  Build a channel which holds the face normals.
     *
     *  The resulting normals are usually unit vectors.  However, if a face has
     * coincident vertices, then its normal will have magnitude 0.
     *
     * @param faceChannelName the name of a float32[3] face channel to create
     *		and fill with the mesh's face normals.
     */
    void build_face_normals( const frantic::tstring& faceChannelName = _T("Normal") );

    /**
     * Clears out all named vertex channels, leaving only the raw geometry
     */
    void erase_all_named_channels();

    /// TODO: ideally these checks could be generalized into mesh_interface to
    /// make the accessible from all mesh types, something to consider?

    /**
     * Check that the sizes of the custom channels are consistent with the base mesh, i.e.
     * - All vertex channels without custom faces have the same number of vertices
     * - All vertex channels with custom faces have the same number of faces
     * - All face channels have the same number of faces
     *
     * @param out An output stream to print error messages to
     * @return false if a duplicate edge was detected
     */
    bool check_channel_sizes( std::ostream& out, frantic::logging::progress_logger* logger = NULL ) const;

    /**
     * Check the mesh for duplicate directed edges. Each pair of vertices should have at most one
     * directed edge between them.
     *
     * @param out An output stream to print error messages to
     * @return false if a duplicate edge was detected
     */
    bool check_duplicate_edges( std::ostream& out, frantic::logging::progress_logger* logger = NULL ) const;

    /**
     * Check that the mesh vertices all have non-infinite coordinates.
     *
     * @param out An output stream to print error messages to
     * @return false if an infinite vertex was detected
     */
    bool check_finite_vertices( std::ostream& out, frantic::logging::progress_logger* logger = NULL ) const;

    /**
     * Check that all mesh indices are within the range of available vertices.
     * Will also check all named channels for consistency as well.
     *
     * @param out An output stream to print error messages to
     * @return false if an out of range index was detected
     */
    bool check_index_ranges( std::ostream& out, frantic::logging::progress_logger* logger = NULL ) const;

    /**
     * Faces which report having zero area (typically degenerate slivers)
     * are problematic since they cause problems with surface normal generation.
     *
     * @param out An output stream to print error messages to
     * @return false if a zero area face was detected
     */
    bool check_zero_area_faces( std::ostream& out, frantic::logging::progress_logger* logger = NULL,
                                float epsilon = 1.0e-8f ) const;

    /**
     * Checks if any face has more than one of the same vertex index.
     * Applies to both the primary, and custom face channels.
     *
     * @param out An output stream to print error messages to
     * @return false if a bad face was detected
     */
    bool check_duplicate_face_indices( std::ostream& out, frantic::logging::progress_logger* logger = NULL ) const;

    /**
     * This function verifies that the mesh is ok.  Something that is not ok is a face with an index outside of its
     * valid index range.
     *
     * @param  out  The output stream where error messages are printed to.
     * @param flags The set of checks to perform
     * @return false if there was an error detected in the mesh
     */
    bool check_consistency( std::ostream& out,
                            consistency_check_flag flags = consistency_check_flag( infinite_vertices |
                                                                                   out_of_range_indices ),
                            frantic::logging::progress_logger* logger = NULL ) const;
};

/**
 *  Copy the data stored in the meshInterface into the outMesh.
 *
 *  an exception will be thrown if the meshInterface holds a non-triangle mesh and autoTriangulate is set to false.
 *
 * @param meshInterface the mesh data to copy.
 * @param[out] outMesh the mesh to copy the data into.
 * @param autoTriangulate true to triangulate the mesh if meshInterface is a non-triangle mesh
 */
void copy_to_trimesh3( const frantic::geometry::mesh_interface* meshInterface, frantic::geometry::trimesh3& outMesh,
                       bool autoTriangulate = true );

} // namespace geometry
} // namespace frantic
