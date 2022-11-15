// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <fstream>
#include <memory>

#include <boost/filesystem/path.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/u3d/u3d_common.hpp>
#include <frantic/geometry/u3d/u3d_shading.hpp>

#ifdef OIIO_LIB_AVAILABLE
#include <OpenImageIO/imagebuf.h>
#endif

namespace frantic {
namespace geometry {

/**
 * An object that writes a mesh_interface object along with a texture to a u3d file
 */
class u3d_writer {
    boost::filesystem::path m_filename;
    std::unique_ptr<std::ostream> m_file;

    u3d::u3d_file_header_block m_header;
    u3d::u3d_modifier_chain_block m_nodeModifierChain;
    u3d::u3d_model_node_block m_nodeBlock;
    u3d::u3d_modifier_chain_block m_resourceModifierChain;
    u3d::u3d_clod_mesh_declaration m_meshDeclaration;
    u3d::u3d_clod_base_mesh_continuation m_mesh;
    u3d::u3d_shading_description m_shader;
    u3d::u3d_shading_modifier_block m_shadingModifier;
    u3d::u3d_lit_texture_shader_block m_textureShader;
    u3d::u3d_material_resource_block m_material;
    u3d::u3d_modifier_chain_block m_textureModifierChain;
    u3d::u3d_texture_declaration_block m_textureDeclaration;
    u3d::u3d_texture_continuation_block m_texture;

    bool m_writeVertexColors;

    /**
     * Writes the file header block of the .U3D file
     */
    void write_file_header();

    /**
     * Initializes the data blocks related to the mesh that would be written to the .U3D file
     *
     * @param mesh Mesh that needs to be written to the file
     */
    void initialize_mesh_data( const mesh_interface* mesh );

#ifdef OIIO_LIB_AVAILABLE
    /**
     * Initializes the data blocks related to the mesh and the texture that would be written to the .U3D file
     *
     * @param mesh Mesh that needs to be written to the file
     * @param texture Texture that needs to be written to the file
     *
     * @note This function sets the flag writeVertexColors to false. This is beacuse writing both vertex colors and
     * texture coordinates together was resulting in a 3D data parsing error
     */
    void initialize_texture_data( const mesh_interface* mesh, const OIIO::ImageBuf* texture );
#endif

  public:
    /**
     * Creates a writer object that will write a .U3D file
     *
     * @param name Name of the output .U3D file
     */
    u3d_writer( const boost::filesystem::path& name );
    ~u3d_writer();

    /**
     * Writes a mesh object to a .U3D file
     *
     * @param mesh Mesh that needs to be written to the file
     */
    void write_mesh( const mesh_interface* mesh );

#ifdef OIIO_LIB_AVAILABLE
    /**
     * Writes a mesh object and a corresponding texture to a .U3D file
     *
     * @param mesh Mesh that needs to be written to the file
     * @param texture Texture that needs to be written to the file
     */
    void write_mesh_with_texture( const mesh_interface* mesh, const OIIO::ImageBuf* texture );
#endif

    /**
     * Closes the stream
     */
    void close();
};

} // namespace geometry
} // namespace frantic
