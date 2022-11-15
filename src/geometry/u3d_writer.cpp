// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/u3d_writer.hpp>

using namespace frantic::geometry;
using namespace frantic::geometry::u3d;

u3d_writer::u3d_writer( const boost::filesystem::path& name )
    : m_filename( name ) {
    m_file.reset( new std::ofstream( m_filename.c_str(), std::ios::binary ) );
    if( !*m_file ) {
        throw std::runtime_error( "u3d_writer::u3d_writer Error: The file " + m_filename.string() +
                                  " failed to open " );
    }
    m_writeVertexColors = true;
    m_header.set_default_values();
    write_file_header();
}

void u3d_writer::write_file_header() { m_header.write( *m_file ); }

u3d_writer::~u3d_writer() { close(); }

void u3d_writer::close() { m_file.reset(); }

void u3d_writer::initialize_mesh_data( const mesh_interface* mesh ) {
    m_nodeModifierChain.set_default_values( 0 );
    m_nodeBlock.set_default_values();
    m_resourceModifierChain.set_default_values( 1 );
    m_meshDeclaration.set_default_values();
    m_mesh.set_default_values();
    m_shader.set_default_values( 0 );

    m_nodeModifierChain.name = _T("MeshNode");
    m_resourceModifierChain.name = _T("MeshResource");
    m_resourceModifierChain.padding = u3d::get_padding_guess( m_resourceModifierChain.name );
    m_nodeBlock.modelNodeName = m_nodeModifierChain.name;
    m_nodeBlock.modelResourceName = m_resourceModifierChain.name;
    m_meshDeclaration.name = m_resourceModifierChain.name;
    m_mesh.name = m_resourceModifierChain.name;

    boost::uint32_t numVertices = static_cast<boost::uint32_t>( mesh->get_num_verts() );
    boost::uint32_t numFaces = static_cast<boost::uint32_t>( mesh->get_num_faces() );

    m_meshDeclaration.minimumResolution = numVertices;
    m_meshDeclaration.finalMaximumResolution = numVertices;

    m_meshDeclaration.maxMeshDescription.faceCount = numFaces;
    m_meshDeclaration.maxMeshDescription.positionCount = numVertices;
    m_meshDeclaration.maxMeshDescription.shadingDescription[0] = m_shader;

    m_mesh.set_positions( mesh );
    m_mesh.set_faces( mesh );

    if( m_writeVertexColors && mesh->has_vertex_channel( _T("Color") ) ) {
        m_mesh.set_vertex_colors( mesh );
        m_meshDeclaration.maxMeshDescription.shadingDescription[0].shadingAttributes = 1;
        const frantic::geometry::mesh_channel* channel = mesh->get_vertex_channel( _T("Color") );
        m_meshDeclaration.maxMeshDescription.diffuseColorCount =
            static_cast<boost::uint32_t>( channel->get_num_elements() );
    }

    m_nodeModifierChain.dataSize +=
        4 + static_cast<boost::uint32_t>( m_nodeModifierChain.name.length() ) + m_nodeBlock.update_and_get_size() + 12;
    m_resourceModifierChain.dataSize += 2 + static_cast<boost::uint32_t>( m_resourceModifierChain.name.length() ) +
                                        m_resourceModifierChain.padding + m_meshDeclaration.update_and_get_size() + 12;
    m_mesh.update_and_get_size();
}

#ifdef OIIO_LIB_AVAILABLE
void u3d_writer::initialize_texture_data( const mesh_interface* mesh, const OIIO::ImageBuf* texture ) {
    m_textureDeclaration.set_default_values();
    m_texture.set_default_values();
    m_textureModifierChain.set_default_values( 2 );
    m_material.set_default_values();
    m_textureShader.set_default_values();
    m_shadingModifier.set_default_values();

    m_shadingModifier.name = _T("MeshNode");
    m_shadingModifier.shaderName[0] = _T("TextureShader");
    m_nodeModifierChain.modifierCount++;
    m_nodeModifierChain.dataSize += 12 + m_shadingModifier.update_and_get_size();

    m_material.name = _T("Material");
    m_material.dataSize += 2 + static_cast<boost::uint32_t>( m_material.name.length() );
    m_textureModifierChain.name = _T("Texture");
    m_textureDeclaration.name = _T("Texture");
    m_textureDeclaration.dataSize += 2 + static_cast<boost::uint32_t>( m_textureDeclaration.name.size() );
    // Include padding to dataSize (round up to a multiple of 4 bytes)
    m_textureDeclaration.dataSize = ( m_textureDeclaration.dataSize + 3 ) & ~0x03;
    m_texture.name = _T("Texture");
    m_texture.dataSize += 2 + static_cast<boost::uint32_t>( m_texture.name.size() );

    m_textureShader.name = _T("TextureShader");
    m_textureShader.materialName = _T("Material");
    m_textureShader.textureInformation.name = _T("Texture");
    m_textureShader.dataSize += 6 + static_cast<boost::uint32_t>( m_textureShader.name.length() ) +
                                static_cast<boost::uint32_t>( m_textureShader.materialName.length() ) +
                                static_cast<boost::uint32_t>( m_textureShader.textureInformation.name.length() );

    m_textureModifierChain.padding = u3d::get_padding_guess( m_textureModifierChain.name );
    m_textureModifierChain.dataSize += 4 + static_cast<boost::uint32_t>( m_textureModifierChain.name.length() ) +
                                       m_textureModifierChain.padding + m_textureDeclaration.dataSize + 12;

    m_meshDeclaration.maxMeshDescription.textureCoordCount =
        static_cast<boost::uint32_t>( mesh->get_vertex_channel( _T("TextureCoord") )->get_num_elements() );
    m_meshDeclaration.maxMeshDescription.shadingDescription[0].set_texture_layer_count( 1 );
    m_meshDeclaration.maxMeshDescription.shadingDescription[0].textureCoordDimensions[0] = 2;
    boost::uint32_t oldSize = m_meshDeclaration.dataSize;
    boost::uint32_t meshDeclarationSize = m_meshDeclaration.update_and_get_size();
    m_resourceModifierChain.dataSize += meshDeclarationSize - oldSize;
    m_mesh.set_textures( mesh );
    m_mesh.update_and_get_size();

    m_texture.set_image_data( texture );
    m_textureDeclaration.set_image_metadata( texture );
    m_textureDeclaration.continuationImageFormats[0].imageDataByteCount =
        m_texture.dataSize - ( 4 + 2 + static_cast<boost::uint32_t>( m_texture.name.size() ) );
}
#endif

void u3d_writer::write_mesh( const mesh_interface* mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "u3d_writer::write_mesh ERROR: mesh was Null" );
    }
    initialize_mesh_data( mesh );
    m_nodeModifierChain.write( *m_file );
    m_nodeBlock.write( *m_file );
    m_resourceModifierChain.write( *m_file );
    m_meshDeclaration.write( *m_file );
    m_mesh.write( *m_file );

    boost::uint32_t fileSize = static_cast<boost::uint32_t>( m_file->tellp() );
    m_file->seekp( 0, std::ios::beg );
    m_header.fileSize = fileSize - 1;
    m_header.declarationSize = fileSize - m_mesh.dataSize - 13;
    write_file_header();
}

#ifdef OIIO_LIB_AVAILABLE
void u3d_writer::write_mesh_with_texture( const mesh_interface* mesh, const OIIO::ImageBuf* texture ) {
    if( !mesh ) {
        throw std::runtime_error( "u3d_writer::write_mesh_with_texture ERROR: mesh was Null" );
    }
    if( !texture ) {
        throw std::runtime_error( "u3d_writer::write_mesh_with_texture : Texture was NULL." );
    }
    m_writeVertexColors = false;
    initialize_mesh_data( mesh );
    initialize_texture_data( mesh, texture );
    m_nodeModifierChain.write( *m_file );
    m_nodeBlock.write( *m_file );
    m_shadingModifier.write( *m_file );
    m_resourceModifierChain.write( *m_file );
    m_meshDeclaration.write( *m_file );
    m_textureModifierChain.write( *m_file );
    m_textureDeclaration.write( *m_file );
    m_textureShader.write( *m_file );
    m_material.write( *m_file );
    m_header.declarationSize = static_cast<boost::uint32_t>( m_file->tellp() );
    m_mesh.write( *m_file );
    m_texture.write( *m_file );

    boost::uint32_t fileSize = static_cast<boost::uint32_t>( m_file->tellp() );
    m_file->seekp( 0, std::ios::beg );
    m_header.fileSize = fileSize;
    write_file_header();
}
#endif
