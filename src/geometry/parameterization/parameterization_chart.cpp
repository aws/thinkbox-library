// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/parameterization/parameterization_chart.hpp>

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/polygon_utils.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/misc/range_segmentation.hpp>

namespace frantic {
namespace geometry {
namespace parameterization {

using graphics2d::vector2f;

parameterization_chart::parameterization_chart( const mesh_interface_ptr& mesh )
    : hasDcel( false ) {
    const size_t vertexCount = mesh->get_num_verts();
    const size_t faceCount = mesh->get_num_faces();

    m_vertices.reserve( vertexCount );
    m_mappings.resize( vertexCount, vector2f( 0.0f, 0.0f ) );
    m_originalFaces.reserve( faceCount );

    // Populate vertex data
    for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
        m_vertices.push_back( mesh->get_vert( vertexIndex ) );
    }

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

    // Populate face verts
    for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        m_originalFaces.push_back( faceIndex );

        for( size_t faceVert = 0; faceVert < mesh->get_num_face_verts( faceIndex ); ++faceVert ) {
            m_faceVerts.push_back( mesh->get_face_vert_index( faceIndex, faceVert ) );
        }
    }
};

void parameterization_chart::split( const range_segmentation& segmentation, std::vector<ptr_type>& output ) const {
    const size_t chartCount = segmentation.get_num_subsets();
    output.reserve( output.size() + chartCount );
    if( chartCount == 1 ) {
        output.push_back( parameterization_chart_ptr( new parameterization_chart( *this ) ) );
    } else {
        for( size_t chart = 0; chart < chartCount; ++chart ) {
            output.push_back( parameterization_chart_ptr( new parameterization_chart(
                *this, segmentation.subset_begin( chart ), segmentation.subset_end( chart ) ) ) );
        }
    }
}

mesh_interface_ptr parameterization_chart::to_mesh( const tstring& texCoordChannelName ) const {
    using graphics::vector3f;

    polymesh3_builder builder;

    // Construct the mesh geometry
    const size_t vertexCount = m_vertices.size();
    for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
        builder.add_vertex( m_vertices[vertexIndex] );
    }

    std::vector<int> faceIndices;

    const size_t faceCount = get_num_faces();
    for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        const size_t faceSize = get_num_face_verts( faceIndex );

        if( faceSize > faceIndices.size() ) {
            faceIndices.resize( faceSize );
        }

        for( size_t faceVert = 0; faceVert < faceSize; ++faceVert ) {
            faceIndices[faceVert] = int( get_face_vert_index( faceIndex, faceVert ) );
        }

        builder.add_polygon( &faceIndices[0], faceSize );
    }

    polymesh3_ptr mesh = builder.finalize();

    // Create and add the mappings
    std::vector<vector3f> uvs;
    uvs.reserve( vertexCount );
    for( size_t mapping = 0; mapping < vertexCount; ++mapping ) {
        uvs.push_back( vector3f( m_mappings[mapping].x, m_mappings[mapping].y, 0.0f ) );
    }

    graphics::raw_byte_buffer buffer( &uvs[0], vertexCount * sizeof( vector3f ) );
    mesh->add_vertex_channel( texCoordChannelName, channels::data_type_float32, 3, buffer );

    return mesh_interface_ptr( polymesh3_interface::create_instance( mesh ).release() );
}

void parameterization_chart::output_as_obj( const tstring& filePath ) const {
    logging::null_progress_logger progress;
    write_obj_mesh_file( filePath, to_mesh( _T("TextureCoord") ), progress );
}

bool parameterization_chart::isTriangulated() const { return m_faceOffsets.size() == 0; }

size_t parameterization_chart::get_num_vertices() const { return m_vertices.size(); }

size_t parameterization_chart::get_num_faces() const { return m_originalFaces.size(); }

size_t parameterization_chart::get_num_face_verts( size_t face ) const {
    if( isTriangulated() ) {
        return 3;
    } else {
        return m_faceOffsets[face + 1] - m_faceOffsets[face];
    }
}

size_t parameterization_chart::get_face_vert_index( size_t face, size_t index ) const {
    if( isTriangulated() ) {
        return m_faceVerts[face * 3 + index];
    } else {
        return m_faceVerts[m_faceOffsets[face] + index];
    }
}

size_t parameterization_chart::get_original_face( size_t face ) const { return m_originalFaces[face]; }

const graphics::vector3f& parameterization_chart::get_vertex( size_t vertex ) const { return m_vertices[vertex]; }

const vector2f& parameterization_chart::get_mapping( size_t mapping ) const { return m_mappings[mapping]; }

void parameterization_chart::set_mapping( size_t index, const vector2f& mapping ) { m_mappings[index] = mapping; }

parameterization_chart::face_vertex_iterator parameterization_chart::face_begin( size_t face ) const {
    if( isTriangulated() ) {
        return face_vertex_iterator( this, face * 3 );
    } else {
        return face_vertex_iterator( this, m_faceOffsets[face] );
    }
}

parameterization_chart::face_vertex_iterator parameterization_chart::face_end( size_t face ) const {
    if( isTriangulated() ) {
        return face_vertex_iterator( this, ( face + 1 ) * 3 );
    } else {
        return face_vertex_iterator( this, m_faceOffsets[face + 1] );
    }
}

parameterization_chart::face_mapping_iterator parameterization_chart::face_mapping_begin( size_t face ) const {
    if( isTriangulated() ) {
        return face_mapping_iterator( this, face * 3 );
    } else {
        return face_mapping_iterator( this, m_faceOffsets[face] );
    }
}

parameterization_chart::face_mapping_iterator parameterization_chart::face_mapping_end( size_t face ) const {
    if( isTriangulated() ) {
        return face_mapping_iterator( this, ( face + 1 ) * 3 );
    } else {
        return face_mapping_iterator( this, m_faceOffsets[face + 1] );
    }
}

parameterization_chart::vertex_iterator parameterization_chart::vertices_begin() const { return m_vertices.begin(); }

parameterization_chart::vertex_iterator parameterization_chart::vertices_end() const { return m_vertices.end(); }

parameterization_chart::mapping_iterator parameterization_chart::mappings_begin() const { return m_mappings.begin(); }

parameterization_chart::mapping_iterator parameterization_chart::mappings_end() const { return m_mappings.end(); }

const dcel& parameterization_chart::get_dcel() {
    if( !hasDcel ) {
        parameterization_chart_to_dcel( *this, m_edgeStructure );
        hasDcel = true;
    }
    return m_edgeStructure;
}

float parameterization_chart::geometric_area() const {
    double area = 0.0f;

    const size_t faceCount = get_num_faces();
    for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        area += polygon3_area( face_begin( faceIndex ), face_end( faceIndex ) );
    }

    return static_cast<float>( area );
}

float parameterization_chart::parametric_area() const {
    double area = 0.0f;

    const size_t faceCount = get_num_faces();
    for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        area += polygon2_area( face_mapping_begin( faceIndex ), face_mapping_end( faceIndex ) );
    }

    return static_cast<float>( area );
}

float parameterization_chart::geometric_perimeter() const {
    if( !hasDcel ) {
        throw std::runtime_error( "parameterization_chart::geometric_perimeter Error: Chart has no initialized DCEL" );
    }

    float perimiter = 0.0f;

    const size_t boundaryCount = m_edgeStructure.boundary_count();
    for( size_t boundaryIndex = 0; boundaryIndex < boundaryCount; ++boundaryIndex ) {
        dcel::const_halfedge_handle edge = m_edgeStructure.get_boundary_halfedge( boundaryIndex );
        perimiter += graphics::vector3f::distance( m_vertices[edge.target_vertex()], m_vertices[edge.source_vertex()] );
    }

    return perimiter;
}

float parameterization_chart::parametric_perimeter() const {
    if( !hasDcel ) {
        throw std::runtime_error( "parameterization_chart::parametric_perimeter Error: Chart has no initialized DCEL" );
    }

    float perimiter = 0.0f;

    const size_t boundaryCount = m_edgeStructure.boundary_count();
    for( size_t boundaryIndex = 0; boundaryIndex < boundaryCount; ++boundaryIndex ) {
        dcel::const_halfedge_handle edge = m_edgeStructure.get_boundary_halfedge( boundaryIndex );
        perimiter += vector2f::distance( m_mappings[edge.target_vertex()], m_mappings[edge.source_vertex()] );
    }

    return perimiter;
}

vector2f parameterization_chart::parametric_centroid() const {
    vector2f positionSum( 0.0f );
    const size_t vertexCount = get_num_vertices();
    for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
        positionSum += m_mappings[vertexIndex];
    }
    return positionSum / float( vertexCount );
}

} // namespace parameterization
} // namespace geometry
} // namespace frantic
