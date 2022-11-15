// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/move/core.hpp>

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/misc/range_segmentation.hpp>

#include <boost/range/iterator_range.hpp>

#include <iterator>
#include <vector>

// Forward declaration
namespace frantic {
namespace particles {
class particle_array;
}
} // namespace frantic

namespace frantic {
namespace geometry {

/*
 * Compares all verts, faces and channel data
 * \return true if the two meshes are identical
 */
bool is_equal( const frantic::geometry::mesh_interface* mesh, const frantic::geometry::mesh_interface* otherMesh );
inline bool is_equal( const frantic::geometry::mesh_interface_ptr& mesh,
                      const frantic::geometry::mesh_interface_ptr& otherMesh ) {
    return is_equal( mesh.get(), otherMesh.get() );
}

/**
 *  Create a bounding box that contains all of the mesh's vertex positions.
 *
 * @param mesh the mesh to build a boundbox for.
 * @return a bounding box that contains all of the mesh's vertex positions.
 */
frantic::graphics::boundbox3f compute_boundbox( const frantic::geometry::mesh_interface* mesh );
inline frantic::graphics::boundbox3f compute_boundbox( const mesh_interface_ptr& mesh ) {
    return compute_boundbox( mesh.get() );
}

/**
 *  Return the total number of vertices summed over all faces in the mesh,
 *
 *  For example, a mesh with 2 triangles will return 3 vertices * 2 triangles = 6.
 *
 * @param mesh the mesh to examine.
 * @return the number of face-vertices in the mesh.
 */
std::size_t get_face_vertex_count( const frantic::geometry::mesh_interface* mesh );
inline std::size_t get_face_vertex_count( const mesh_interface_ptr& mesh ) {
    return get_face_vertex_count( mesh.get() );
}

/**
 *  Test if the mesh is a closed manifold.  This is determined by checking
 * if all edges are incident on exactly two faces.
 *
 * @param mesh the mesh to test.
 * @return true if the mesh is a closed manifold, and false otherwise.
 */
bool is_closed_manifold( const frantic::geometry::mesh_interface* mesh );
inline bool is_closed_manifold( const frantic::geometry::mesh_interface_ptr& mesh ) {
    return is_closed_manifold( mesh.get() );
}

/**
 * Test if the mesh has degenerate faces. Faces are considered degenerate if they contain multiple references to the
 * same vertex.
 *
 * @param mesh The mesh to test.
 * @return true if the mesh has degenerate faces, and false otherwise.
 */
bool has_degenerate_faces( const mesh_interface_ptr& mesh );

/**
 * Fix all degenerate faces in the mesh. Faces are considered degenerate if they contain multiple references to the same
 * vertex. Depending on the topolgy of the face with respect to the duplicate vertex, the face will either be
 * simplified, removed, or split into multiple faces. Additional channels with custom faces will be adjusted
 * appropriately to reflect the topology change.
 *
 * @param mesh The mesh to fix.
 * @return A deep copy of the mesh with degenerate faces fixed.
 */
polymesh3_ptr fix_degenerate_faces( const mesh_interface_ptr& mesh );

/**
 *  Compute the vertex normal for all vertices in the mesh.
 *
 *  The output vertex normals are normalized to have length 1.
 *
 * @param mesh the mesh for which to compute vertex normals.
 * @param[out] outVertexNormals the computed vertex normals.
 */
void compute_vertex_normals( const frantic::geometry::mesh_interface* mesh,
                             std::vector<frantic::graphics::vector3f>& outVertexNormals );

void compute_face_normals( const frantic::geometry::mesh_interface* mesh,
                           std::vector<frantic::graphics::vector3f>& outFaceNormals );

/**
 *  Return a read-only mesh_interface that includes all of the vertices,
 * faces, and channels of the input meshInterface, plus a "Normal" vertex
 * channel.
 *
 * @note the returned object holds an internal reference to the input
 *    meshInterface.  The caller must not delete or modify the input
 *    meshInterface or any of its channels during the lifetime of the
 *    returned object.
 *
 * @param meshInterface the mesh interface for which to create a "Normal"
 *    channel.
 * @return a new mesh_interface that includes all of the vertices, faces,
 *    and channels of the input meshInterface, plus a "Normal" vertex
 *    channel.
 */
std::unique_ptr<frantic::geometry::mesh_interface>
create_vertex_normal_channel_mesh_interface( const frantic::geometry::mesh_interface* meshInterface );

std::unique_ptr<frantic::geometry::mesh_interface>
create_face_normal_channel_mesh_interface( const frantic::geometry::mesh_interface* meshInterface );

/**
 *  Return a read-only mesh_interface that includes all of the vertices,
 * faces, and channels of the input meshInterface, plus an "FaceArea" face
 * channel, containing the surface area of each face
 *
 * @param mesh the mesh to override
 * @param dataType the type to store the area channel as, must be a valid floating-point type
 * @return a new mesh_interface with all of the existing channels, and the added "FaceArea" channel
 */
std::unique_ptr<mesh_interface>
create_face_area_channel_mesh_interface( const mesh_interface* mesh,
                                         frantic::channels::data_type_t dataType = frantic::channels::data_type_float32,
                                         const frantic::tstring& name = _T( "FaceArea" ) );

/**
 * Returns a 'shallow-copy' mesh interface that stores an independent set of vertex locations,
 * such that modifying their locations on the delegate will not modify the underlying mesh
 */
std::unique_ptr<frantic::geometry::mesh_interface>
create_modified_geometry_mesh_interface( const frantic::geometry::mesh_interface* meshInterface );

/**
 *  Create a vertex cloud mesh from a particle array.
 *
 *  The returned mesh_interface includes all of the particles as vertices, and
 * all of their channels as vertex channels.
 *
 * @note the returned mesh_interface is read-only.
 *
 * @return a vertex cloud mesh.
 * @param particles the particles to convert to a vertex cloud.
 */
std::unique_ptr<frantic::geometry::mesh_interface>
    create_particle_array_mesh_interface( BOOST_RV_REF( frantic::particles::particle_array ) particles );

std::unique_ptr<mesh_interface> make_submesh( const mesh_interface* principleMesh, const std::vector<size_t>& faces );

template <class InputIterator>
std::unique_ptr<mesh_interface> make_submesh( const mesh_interface* principleMesh, InputIterator facesBegin,
                                              InputIterator facesEnd ) {
    return make_submesh( principleMesh, std::vector<size_t>( facesBegin, facesEnd ) );
}

/**
 * Check to see if a given mesh is a triangle mesh.
 *
 * @param mesh the mesh to perform the check on.
 * @return true if the mesh is a triangle mesh, false if not.
 */
bool is_triangle_mesh( const mesh_interface* mesh );

/**
 *  Check if all of the indices in the mesh and its channels are in bounds.
 * Throw an exception if any index is out of bounds.
 *
 * @param mesh the mesh to check.
 * @param meshNameForErrorMessage the mesh's name, to help identify it in
 *                                error messages.
 */
void assert_valid_indices( const mesh_interface* mesh, const frantic::tstring& meshNameForErrorMessage = _T("mesh") );

/**
 * Provides a standards-compliant iterator through the vertices of a mesh_interface face
 *
 * TODO: make this a RandomAccessIterator
 */
class face_vertex_iterator {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef frantic::graphics::vector3f value_type;
    typedef ptrdiff_t difference_type;
    typedef frantic::graphics::vector3f* pointer;
    typedef frantic::graphics::vector3f& reference;

    enum tag_location {
        begin,
        end,
    };

  private:
    const mesh_interface* m_mesh;
    size_t m_face;
    size_t m_fVert;

  public:
    face_vertex_iterator( const mesh_interface* mesh, size_t face, tag_location location = begin )
        : m_mesh( mesh )
        , m_face( face ) {
        if( location == begin ) {
            m_fVert = 0;
        } else {
            m_fVert = m_mesh->get_num_face_verts( m_face );
        }
    }

    face_vertex_iterator( const mesh_interface* mesh, size_t face, size_t location )
        : m_mesh( mesh )
        , m_face( face )
        , m_fVert( location ) {}

    face_vertex_iterator( const face_vertex_iterator& other )
        : m_mesh( other.m_mesh )
        , m_face( other.m_face )
        , m_fVert( other.m_fVert ) {}

    face_vertex_iterator& operator=( const face_vertex_iterator& other ) {
        m_mesh = other.m_mesh;
        m_face = other.m_face;
        m_fVert = other.m_fVert;
        return *this;
    }

    bool related( const face_vertex_iterator& other ) const { return m_mesh == other.m_mesh && m_face == other.m_face; }

    bool operator==( const face_vertex_iterator& other ) const { return related( other ) && m_fVert == other.m_fVert; }

    bool operator!=( const face_vertex_iterator& other ) const { return !( *this == other ); }

    face_vertex_iterator& operator++() {
        ++m_fVert;
        return *this;
    }

    face_vertex_iterator operator++( int ) {
        face_vertex_iterator result = *this;
        ++m_fVert;
        return result;
    }

    ptrdiff_t operator-( const face_vertex_iterator& other ) const {
        if( related( other ) ) {
            return ptrdiff_t( m_fVert ) - ptrdiff_t( other.m_fVert );
        } else {
            throw std::runtime_error(
                "face_vertex_iterator::operator- : Error, cannot take the difference of unrelated iterators." );
        }
    }

    frantic::graphics::vector3f operator*() const { return m_mesh->get_face_vert( m_face, m_fVert ); }
};

typedef boost::iterator_range<face_vertex_iterator> face_vertex_range;

inline face_vertex_range make_face_vertex_range( const mesh_interface* mesh, size_t faceId ) {
    return face_vertex_range( face_vertex_iterator( mesh, faceId, face_vertex_iterator::begin ),
                              face_vertex_iterator( mesh, faceId, face_vertex_iterator::end ) );
}

/**
 * Provides a standards-compliant iterator through all vertices of a mesh_interface
 *
 * TODO: make this a RandomAccessIterator
 */
class mesh_vertex_iterator {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef frantic::graphics::vector3f value_type;
    typedef std::ptrdiff_t difference_type;
    typedef frantic::graphics::vector3f* pointer;
    typedef frantic::graphics::vector3f& reference;

    enum tag_location {
        begin,
        end,
    };

  private:
    const mesh_interface* m_mesh;
    size_t m_currentVert;

  public:
    mesh_vertex_iterator( const mesh_interface* mesh, tag_location location = begin )
        : m_mesh( mesh )
        , m_currentVert( location == begin ? 0 : mesh->get_num_verts() ) {}

    mesh_vertex_iterator( const mesh_vertex_iterator& other )
        : m_mesh( other.m_mesh )
        , m_currentVert( other.m_currentVert ) {}

    mesh_vertex_iterator& operator=( const mesh_vertex_iterator& other ) {
        m_mesh = other.m_mesh;
        m_currentVert = other.m_currentVert;
        return *this;
    }

    bool related( const mesh_vertex_iterator& other ) const { return m_mesh == other.m_mesh; }

    bool operator==( const mesh_vertex_iterator& other ) const {
        return related( other ) && m_currentVert == other.m_currentVert;
    }

    bool operator!=( const mesh_vertex_iterator& other ) const { return !( *this == other ); }

    mesh_vertex_iterator& operator++() {
        ++m_currentVert;
        return *this;
    }

    mesh_vertex_iterator operator++( int ) {
        mesh_vertex_iterator result = *this;
        ++m_currentVert;
        return result;
    }

    ptrdiff_t operator-( const mesh_vertex_iterator& other ) const {
        if( related( other ) ) {
            return ptrdiff_t( m_currentVert ) - ptrdiff_t( other.m_currentVert );
        } else {
            throw std::runtime_error(
                "mesh_vertex_iterator::operator- : Error, cannot take the difference of unrelated iterators." );
        }
    }

    frantic::graphics::vector3f operator*() const { return m_mesh->get_vert( m_currentVert ); }
};

class mesh_interface_geom_channel : public mesh_channel {
  private:
    mesh_interface* m_mesh;

  public:
    mesh_interface_geom_channel( mesh_interface* mesh );

    virtual void get_value( std::size_t index, void* outValue ) const;

    virtual void set_value( std::size_t index, const void* value ) const;

    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const;

    virtual std::size_t get_num_face_verts( size_t faceIndex ) const;
};

/**
 * Find the corresponding per-face vertex index on a face given a global vertex index.
 *
 * @param mesh the mesh on which to search
 * @param faceId the face of the mesh to search
 * @param vertexId the global index of the vertex
 * @return the face-index of the specified vertex on the given face
 * @throws std::runtime_error if that vertex was not on the given face
 */
inline size_t find_face_corner( mesh_interface_ptr mesh, size_t faceId, size_t vertexId ) {
    const size_t faceSize = mesh->get_num_face_verts( faceId );

    for( size_t i = 0; i < faceSize; ++i ) {
        if( mesh->get_face_vert_index( faceId, i ) == vertexId ) {
            return i;
        }
    }

    throw std::runtime_error( "Invalid vertex on face." );
}

/**
 * Return a segmentation such that each face-connected component is its own distinct subset
 *
 * @param iMesh the mesh
 * @param outSegmentation (output) stores the resulting segmentation
 */
void get_connected_components_segmentation( const mesh_interface* iMesh, range_segmentation& outSegmentation );

/**
 * Return the implicit segmentation of a given vertex channel
 *
 * @param mesh the mesh to operate on
 * @param customChannelName the name of the mesh channel to operate on
 * @param outSegmentation (output) stores the resulting segmentation found
 */
void get_channel_segmentation( const mesh_interface* mesh, const tstring& customChannelName,
                               range_segmentation& outSegmentation );

/**
 * Return the implicit segmentation of a given face-vertex channel
 *
 * @param mesh the mesh to operate on
 * @param channel the mesh channel to operate on
 * @param outSegmentation (output) stores the resulting segmentation found
 */
void get_channel_segmentation( const mesh_interface* mesh, const mesh_channel* channel,
                               range_segmentation& outSegmentation );

} // namespace geometry
} // namespace frantic
