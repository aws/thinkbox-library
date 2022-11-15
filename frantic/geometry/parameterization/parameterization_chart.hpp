// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/*
 * @file
 * @author Jeff Coukell
 */

#pragma once

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polygon_utils.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/vector2f.hpp>
#include <frantic/math/utils.hpp>
#include <frantic/misc/range_segmentation.hpp>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include <vector>

namespace frantic {
namespace geometry {
class dcel;
} // namespace geometry
} // namespace frantic

namespace frantic {
namespace geometry {
namespace parameterization {

/**
 * Stores the mapping from 3D vertices to 2D planar points for a piece of geometry or 'chart'. This is typically a
 * subset of a larger mesh whose combined charts are stored in a parameterization_atlas object.
 */
class parameterization_chart {
  public:
    typedef boost::shared_ptr<parameterization_chart> ptr_type;

    typedef std::vector<graphics::vector3f>::const_iterator vertex_iterator;
    typedef std::vector<graphics2d::vector2f>::const_iterator mapping_iterator;

    class face_vertex_iterator : public boost::iterator_facade<face_vertex_iterator, const graphics::vector3f,
                                                               std::random_access_iterator_tag> {
      private:
        const parameterization_chart* m_owner;
        size_t m_index;

      public:
        face_vertex_iterator( const parameterization_chart* chart, size_t index )
            : m_owner( chart )
            , m_index( index ) {}

      private:
        friend class boost::iterator_core_access;

        face_vertex_iterator::reference dereference() const {
            return m_owner->m_vertices[m_owner->m_faceVerts[m_index]];
        }

        bool equal( const face_vertex_iterator& other ) const {
            return m_owner == other.m_owner && m_index == other.m_index;
        }

        void increment() { ++m_index; }

        void decrement() { --m_index; }

        void advance( face_vertex_iterator::difference_type n ) { m_index += n; }

        face_vertex_iterator::difference_type distance_to( const face_vertex_iterator& other ) const {
            return other.m_index - m_index;
        }
    };

    class face_mapping_iterator : public boost::iterator_facade<face_mapping_iterator, const graphics2d::vector2f,
                                                                std::random_access_iterator_tag> {
      private:
        const parameterization_chart* m_owner;
        size_t m_index;

      public:
        face_mapping_iterator( const parameterization_chart* chart, size_t index )
            : m_owner( chart )
            , m_index( index ) {}

      private:
        friend class boost::iterator_core_access;

        face_mapping_iterator::reference dereference() const {
            return m_owner->m_mappings[m_owner->m_faceVerts[m_index]];
        }

        bool equal( const face_mapping_iterator& other ) const {
            return m_owner == other.m_owner && m_index == other.m_index;
        }

        void increment() { ++m_index; }

        void decrement() { --m_index; }

        void advance( face_mapping_iterator::difference_type n ) { m_index += n; }

        face_mapping_iterator::difference_type distance_to( const face_mapping_iterator& other ) const {
            return other.m_index - m_index;
        }
    };

  private:
    std::vector<graphics::vector3f> m_vertices;
    std::vector<graphics2d::vector2f> m_mappings;

    std::vector<size_t> m_faceVerts;
    // A list of positions in m_faceVerts marking the start of each face, this list will not be populated if the mesh is
    // triangulated
    std::vector<size_t> m_faceOffsets;
    // A mapping to face indices in the parent mesh for each face
    std::vector<size_t> m_originalFaces;

    bool hasDcel;
    dcel m_edgeStructure;

  public:
    parameterization_chart( const mesh_interface_ptr& mesh );

    /**
     * Construct a chart by including a subset of the faces from the input mesh. The indices of these faces should be
     * provided by the input iterator pair.
     */
    template <class InputIterator>
    parameterization_chart( const mesh_interface_ptr& mesh, InputIterator facesBegin, InputIterator facesEnd )
        : hasDcel( false ) {
        const size_t faceCount = facesEnd - facesBegin;

        m_vertices.reserve( faceCount ); // Conservative estimate
        m_originalFaces.reserve( faceCount );

        // Check if face offsets are necessary
        bool isTriangulated = true;
        for( size_t faceIndex = 0; faceIndex < faceCount && isTriangulated; ++faceIndex ) {
            isTriangulated = mesh->get_num_face_verts( faceIndex ) == 3;
        }

        // Populate face offsets and find number of face verts
        if( !isTriangulated ) {
            m_faceOffsets.reserve( faceCount + 1 );
            m_faceOffsets.push_back( 0 );

            for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                m_faceOffsets.push_back( mesh->get_num_face_verts( faceIndex ) + m_faceOffsets[faceIndex] );
            }

            m_faceVerts.reserve( m_faceOffsets.back() );
        } else {
            m_faceVerts.reserve( faceCount * 3 );
        }

        // Populate vertex data and face verts, keeping track of vertex indices which are already added
        boost::unordered_map<size_t, size_t> vertexRemap;
        for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
            m_originalFaces.push_back( facesBegin[faceIndex] );

            const size_t faceSize = mesh->get_num_face_verts( facesBegin[faceIndex] );
            for( size_t faceVert = 0; faceVert < faceSize; ++faceVert ) {
                const size_t vertexIndex = mesh->get_face_vert_index( facesBegin[faceIndex], faceVert );
                if( vertexRemap.find( vertexIndex ) == vertexRemap.end() ) {
                    vertexRemap[vertexIndex] = m_vertices.size();
                    m_vertices.push_back( mesh->get_vert( vertexIndex ) );
                }
                m_faceVerts.push_back( vertexRemap[vertexIndex] );
            }
        }

        m_mappings.resize( m_vertices.size(), graphics2d::vector2f( 0.0f, 0.0f ) );
    }

    /**
     * Construct a chart by including a subset of the faces from the input chart. The indices of these faces should be
     * provided by the input iterator pair.
     */
    template <class InputIterator>
    parameterization_chart( const parameterization_chart& chart, InputIterator facesBegin, InputIterator facesEnd )
        : hasDcel( false ) {
        const size_t faceCount = facesEnd - facesBegin;

        m_vertices.reserve( faceCount ); // Conservative estimate
        m_originalFaces.reserve( faceCount );

        // Populate face offsets and find number of face verts
        if( !chart.isTriangulated() ) {
            m_faceOffsets.reserve( faceCount + 1 );
            m_faceOffsets.push_back( 0 );

            for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                m_faceOffsets.push_back( chart.get_num_face_verts( faceIndex ) + m_faceOffsets[faceIndex] );
            }

            m_faceVerts.reserve( m_faceOffsets.back() );
        } else {
            m_faceVerts.reserve( faceCount * 3 );
        }

        // Populate vertex data and face verts, keeping track of vertex indices which are already added
        boost::unordered_map<size_t, size_t> vertexRemap;
        for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
            m_originalFaces.push_back( chart.get_original_face( facesBegin[faceIndex] ) );

            const size_t faceSize = chart.get_num_face_verts( facesBegin[faceIndex] );
            for( size_t faceVert = 0; faceVert < faceSize; ++faceVert ) {
                const size_t vertexIndex = chart.get_face_vert_index( facesBegin[faceIndex], faceVert );
                if( vertexRemap.find( vertexIndex ) == vertexRemap.end() ) {
                    vertexRemap[vertexIndex] = m_vertices.size();
                    m_vertices.push_back( chart.get_vertex( vertexIndex ) );
                    m_mappings.push_back( chart.get_mapping( vertexIndex ) );
                }
                m_faceVerts.push_back( vertexRemap[vertexIndex] );
            }
        }
    }

    /**
     * Divide this chart into several charts as specified by the given face segmentation and append the new charts onto
     * the output vector.
     */
    void split( const range_segmentation& segmentation, std::vector<ptr_type>& output ) const;

    /**
     * Convert the chart to a complete mesh object where the geometry is exactly copied and the 2D mapping is stored in
     * the given vertex channel.
     */
    mesh_interface_ptr to_mesh( const tstring& texCoordChannelName = _T("TextureCoord") ) const;

    /**
     * Output this chart as an OBJ mesh file with its mapping applied as the UV channel.
     */
    void output_as_obj( const tstring& filePath ) const;

    bool isTriangulated() const;

    size_t get_num_vertices() const;
    size_t get_num_faces() const;
    size_t get_num_face_verts( size_t face ) const;
    size_t get_face_vert_index( size_t face, size_t index ) const;
    size_t get_original_face( size_t face ) const;

    const graphics::vector3f& get_vertex( size_t vertex ) const;
    const graphics2d::vector2f& get_mapping( size_t mapping ) const;

    void set_mapping( size_t index, const graphics2d::vector2f& mapping );

    face_vertex_iterator face_begin( size_t face ) const;
    face_vertex_iterator face_end( size_t face ) const;

    face_mapping_iterator face_mapping_begin( size_t face ) const;
    face_mapping_iterator face_mapping_end( size_t face ) const;

    vertex_iterator vertices_begin() const;
    vertex_iterator vertices_end() const;

    mapping_iterator mappings_begin() const;
    mapping_iterator mappings_end() const;

    /**
     * If the chart already has a DCEL then return it, otherwise generate one and return it.
     */
    const dcel& get_dcel();

    float geometric_area() const;
    float parametric_area() const;

    float geometric_perimeter() const;
    float parametric_perimeter() const;

    graphics2d::vector2f parametric_centroid() const;
};

typedef parameterization_chart::ptr_type parameterization_chart_ptr;

} // namespace parameterization
} // namespace geometry
} // namespace frantic
