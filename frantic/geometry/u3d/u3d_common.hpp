// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/integer.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/color_rgba_f.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics/vector4f.hpp>
#include <frantic/strings/tstring.hpp>

#define U3D_FILE_HEADER_BLOCKTYPE 0x00443355
#define U3D_MODIFIER_CHAIN_BLOCKTYPE 0xFFFFFF14
#define U3D_CLOD_MESH_DECLARATION_BLOCKTYPE 0xFFFFFF31
#define U3D_CLOD_BASE_MESH_DATA_BLOCKTYPE 0xFFFFFF3B
#define U3D_MODEL_NODE_BLOCKTYPE 0xFFFFFF22

/**
 * This value corresponds to the MIBenum value for UTF-8.
 * For more details, please look at http://www.iana.org/assignments/character-sets
 */
#define U3D_UTF8_CHARACTER_ENCODING 106;

namespace frantic {
namespace geometry {
namespace u3d {

using namespace frantic::graphics;

/**
 * Writes a value of type T to a raw_byte_buffer
 *
 * @param buf The raw byte buffer to whcih the value must be written
 * @param value The value of type T to be written to the buffer
 */
template <typename T>
inline void write_value( raw_byte_buffer& buf, const T& value ) {
    memcpy( buf.add_element( sizeof( T ) ), reinterpret_cast<const char*>( &value ), sizeof( T ) );
}

/**
 * Writes a string to a raw_byte_buffer
 *
 * @param buf The raw byte buffer to which the string must be written
 * @param s   The string to be written
 */
void write_tstring( raw_byte_buffer& buf, const frantic::tstring& s );

/**
 * Writes a vector3f to a raw_byte_buffer
 *
 * @param buf The raw byte bufffer to which the vector3f must be written
 * @param v   The vector3f to be written to the file
 */
void write_vector3f( raw_byte_buffer& buf, const vector3f& v );

/**
 * Writes a vector4f to a raw_byte_buffer
 *
 * @param buf The raw byte bufffer to which the vector4f must be written
 * @param v   The vector4f to be written to the file
 */
void write_vector4f( raw_byte_buffer& buf, const vector4f& v );

/**
 * Writes a color3f to a raw_byte_buffer
 *
 * @param buf The raw byte bufffer to which the color3f must be written
 * @param c   The color3f to be written to the file
 */
void write_color3f( raw_byte_buffer& buf, const color3f& c );

/**
 * Writes a color_rgba_f to a raw_byte_buffer
 *
 * @param buf The raw byte bufffer to which the color3f must be written
 * @param c   The color_rgba_f to be written to the file
 */
void write_color_rgba_f( raw_byte_buffer& buf, const color_rgba_f& c );

/**
 * Writes a transform4f to a raw_byte_buffer
 *
 * @param buf The raw byte bufffer to which the transform4f must be written
 * @param t   The transform4f to be written to the file
 */
void write_transform4f( raw_byte_buffer& buf, const transform4f& t );

/**
 * Writes a raw byte buffer containing a u3d block to the u3d file
 *
 * @param stream The file to which the block must be written
 * @param buf The buffer containing the u3d block to be written
 * @param blockName The name of the u3d block which is being written
 */
void write_buffered_block( std::ostream& stream, const raw_byte_buffer& buf, const std::string& blockName );

/**
 * Returns the padding that needs to be added to the current position in the file for a certain alignment
 *
 * @param position Current position in the file
 * @param alignment The alginment that needs to be followed for the file ( Eg: 4-byte alignment )
 */
boost::uint32_t get_alignment_data_padding( const boost::uint32_t position, const boost::uint32_t alignment );

/**
 * Returns the padding that would be required to be added to maintain 32-bit alignment beacuse of a given string
 *
 * @param s String that may cause a need for alignment
 */
boost::uint32_t get_padding_guess( const frantic::tstring& s );

/**
 * Adds padding to the u3d file
 *
 * @param stream The file to which the padding must be added
 * @param paddingSize Amount of padding which needs to be added to the file
 */
void add_padding( std::ostream& stream, const boost::uint32_t paddingSize );

/**
 * Object which contains information about the basic structure of each block in the u3d file.
 *
 * @var blockType type of object associated with this block
 * @var dataSize size of the data section in bytes
 * @var metaDataSize size of the meta data section in bytes
 *
 * @note For this implementation of u3d we are assuming that no metadata is ever written. So this should always be 0.
 * If metadata is ever to be implemented, please read section 9.2.6 in the official file format
 */
struct u3d_block {
    boost::uint32_t blockType;
    boost::uint32_t dataSize;
    boost::uint32_t metaDataSize;
};

/**
 * Enum to hold the optional features sued by the u3d file.
 *
 * BASE_PROFILE no optional features used
 * EXTENSIBLE_PROFILE file may contain New Object Type blocks
 * NO_COMPRESSION_MODE file does not contain any compressed values
 * DEFINED_UNITS objects in the file are defined with units
 *
 * @note To combine these values, use the OR operator
 */
namespace U3D_PROFILE_IDENTIFIER {
enum { BASE_PROFILE, EXTENSIBLE_PROFILE = 2, NO_COMPRESSION_MODE = 4, DEFINED_UNITS = 8 };
}

/**
 * Object which holds the information for the file header block in a u3d file
 *
 * @var majorVersion version of the file format used to write this file. Currently 0.
 * @var minorVersion maybe used to indicate the version of the encodes used to write this file. Currently 0.
 * @var profileIdentifier a combination of values from the PROFILE_IDENTIFIER enum
 * @var delcarationSize number of bytes in the Declaration Block section of the file.( Includes the size of file header
 * block )
 * @var fileSize total number of bytes in the file
 * @var characterEncoding encoding used for the strings in the file. It is always set to U3D_DEFAULT_CHARACTER_ENCODING
 * @var unitsScalingFactor defines the units used in this file
 */
struct u3d_file_header_block : u3d_block {
    boost::int16_t majorVersion;
    boost::int16_t minorVersion;
    boost::uint32_t profileIdentifier;
    boost::uint32_t declarationSize;
    boost::uint64_t fileSize;
    boost::uint32_t characterEncoding;
    double unitsScalingFactor;

    /**
     * Sets the default values for file header block
     */
    void set_default_values() {
        blockType = U3D_FILE_HEADER_BLOCKTYPE;
        dataSize = 32;
        metaDataSize = 0;
        majorVersion = 0;
        minorVersion = 0;
        profileIdentifier = U3D_PROFILE_IDENTIFIER::NO_COMPRESSION_MODE | U3D_PROFILE_IDENTIFIER::DEFINED_UNITS;
        declarationSize = 32;
        fileSize = 32;
        characterEncoding = U3D_UTF8_CHARACTER_ENCODING;
        unitsScalingFactor = 1.0;
    }

    /**
     * Writes the file header block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_value<boost::int16_t>( buf, majorVersion );
        write_value<boost::int16_t>( buf, minorVersion );
        write_value<boost::uint32_t>( buf, profileIdentifier );
        write_value<boost::uint32_t>( buf, declarationSize );
        write_value<boost::uint64_t>( buf, fileSize );
        write_value<boost::uint32_t>( buf, characterEncoding );
        write_value<double>( buf, unitsScalingFactor );
        write_buffered_block( stream, buf, "File Header" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

/**
 * Object to hold the information about a model node. It is the first modifier in a node modifier chain.
 *
 * This must be a part of a node modifier chain
 *
 * @var modelNodeName Name of the model node
 * @var numParents Number of parents the node has. For our implementation it is always 0.
 * @var modelResourceName Name of the model resource chain used as input to the model node's modifier chain
 * @var modelVisibility Used to indicate whether the front facing or back facing surface should be drawn.
 *      0 : Not visible
 *      1 : Front Visible
 *      2 : Back Visible
 *      3 : Front and back visible
 */
struct u3d_model_node_block : u3d_block {
    frantic::tstring modelNodeName;
    boost::uint32_t numParents;
    frantic::tstring modelResourceName;
    boost::uint32_t modelVisibility;

    /**
     * Sets the default values for model node block
     */
    void set_default_values() {
        blockType = U3D_MODEL_NODE_BLOCKTYPE;
        dataSize = 8;
        metaDataSize = 0;
        numParents = 0;
        modelVisibility = 3;
    }

    /**
     * Updates and returns the dataSize stored by this block
     */
    boost::uint32_t update_and_get_size() {
        dataSize += 4 + static_cast<boost::uint32_t>( modelNodeName.length() ) +
                    static_cast<boost::uint32_t>( modelResourceName.length() );
        return dataSize;
    }

    /**
     * Writes the data of modifier node block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, modelNodeName );
        write_value<boost::uint32_t>( buf, numParents );
        write_tstring( buf, modelResourceName );
        write_value<boost::uint32_t>( buf, modelVisibility );
        write_buffered_block( stream, buf, "Root Node" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

/**
 * Object to represent modifier chain block in the u3d file
 *
 * @var name Name of the modifier chain and also the name of all the modifiers in the chain
 * @var type Indicates the type of modifier chain.
 *      0 : Node modifier chain ( also called instance modifier chain )
 *      1 : Model Resource modifier chain ( also called resource modifier chain )
 *      2 : Texture Resource modifier chain ( also called texture modifier chain )
 * @var attributes Indicates the presence of optional information about the modifier chain
 *      1 : Bounding sphere information present
 *      2 : Axis-aligned bounding box present
 * @var modifierCount Number of modifiers in the modifier chain
 */
struct u3d_modifier_chain_block : u3d_block {
    frantic::tstring name;
    boost::uint32_t type;
    boost::uint32_t attributes;
    boost::uint32_t modifierCount;
    boost::uint32_t padding;

    /**
     * Sets the default values for modifier chain block
     */
    void set_default_values( boost::uint32_t chainType ) {
        blockType = U3D_MODIFIER_CHAIN_BLOCKTYPE;
        dataSize = 12;
        metaDataSize = 0;
        attributes = 0;
        modifierCount = 1;
        padding = 0;
        type = chainType;
    }

    /**
     * Writes the data of modifier chain block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        std::string blockName;
        if( type == 1 ) {
            blockName = "Node Modifier Chain";
        } else if( type == 2 ) {
            blockName = "Resource Modifier Chain";
        } else {
            blockName = "Texture Modifier Chain";
        }
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, name );
        write_value<boost::uint32_t>( buf, type );
        write_value<boost::uint32_t>( buf, attributes );
        write_buffered_block( stream, buf, blockName );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
        stream.write( reinterpret_cast<const char*>( &modifierCount ), 4u );
    }
};

/**
 * Object to hold the shading description for 1 shading list used in the mesh. Indicates which per vertex attributes, in
 * addition to position and normal, are used by each shader list
 *
 * @var shadingAttributes Collection of flags to indicate the usage of per verte colors. Combine using the binary OR
 * operator 0 : Shader list uses neither diffuse colors nor specular colors 1 : Shader list uses per vertex diffuse
 * colors 2 : Shader list uses per vertex specular colors 3 : Shader list uses both diffuse and specular colors, per
 * vertex
 * @var textureLayerCount Number of texture layers used by this shader list
 * @var textureCoordDimensions Holds the number of dimensions in the texture coordinate vector for each texture layer
 * @var originalShadingID Original shading index for the shader list. Not needed in dcoding the file.
 */
struct u3d_shading_description {
  private:
    boost::uint32_t textureLayerCount;

  public:
    boost::uint32_t shadingAttributes;
    std::vector<boost::uint32_t> textureCoordDimensions;
    boost::uint32_t originalShadingID;

    /**
     * Returns the dataSize stored by this block
     */
    boost::uint32_t size() {
        boost::uint32_t dataSize = 12 + 4 * static_cast<boost::uint32_t>( textureCoordDimensions.size() );
        return dataSize;
    }

    void set_texture_layer_count( boost::uint32_t count ) {
        textureLayerCount = count;
        textureCoordDimensions.resize( textureLayerCount );
    }

    /**
     * Sets the default values for shading description
     *
     * @param shadingID id of the shader
     */
    void set_default_values( boost::uint32_t shadingID ) {
        textureLayerCount = 0;
        shadingAttributes = 0;
        originalShadingID = shadingID;
    }

    /**
     * Writes the data of shading description to a raw byte buffer
     *
     * @param buf The raw byte buffer to which the data must be written to.
     */
    void write( raw_byte_buffer& buf ) {
        write_value<boost::uint32_t>( buf, shadingAttributes );
        write_value<boost::uint32_t>( buf, textureLayerCount );
        for( boost::uint32_t i = 0; i < textureLayerCount; ++i ) {
            write_value<boost::uint32_t>( buf, textureCoordDimensions[i] );
        }
        write_value<boost::uint32_t>( buf, originalShadingID );
    }
};

/**
 * Object that defines the size of the mesh at full resolutiion. This can be used to allocate space for the mesh.
 *
 * @var meshAttributes Contains information that applies to the entire mesh. Only 1 flag defined currently.
 *      0 : Default. The faces in the mesh have a normal index at each corner
 *      1 : Exclude Normals. The faces in the mesh do not have a normal index at each corner.
 * @var faceCount Number of faces in the mesh
 * @var positionCount Number of positions in the position array
 * @var normalCount Number of normals in the normal array
 * @var diffuseColorCount Number of colors in the diffuse color array
 * @var specularColorCount Number of colors in the specular color array
 * @var textureCoordCount Number of texture coordinates in the texture coordinate array
 * @var shadingCount Number of shading descriptions used in the mesh
 * @var shadingDescription Holds the description for all the shader lists
 */
struct u3d_max_mesh_description {
  private:
    boost::uint32_t shadingCount;

  public:
    boost::uint32_t meshAttributes;
    boost::uint32_t faceCount;
    boost::uint32_t positionCount;
    boost::uint32_t normalCount;
    boost::uint32_t diffuseColorCount;
    boost::uint32_t specularColorCount;
    boost::uint32_t textureCoordCount;
    std::vector<u3d_shading_description> shadingDescription;

    /**
     * Returns the dataSize stored by this block
     */
    boost::uint32_t size() {
        boost::uint32_t dataSize = 32;
        for( size_t i = 0; i < shadingCount; ++i ) {
            dataSize += shadingDescription[i].size();
        }
        return dataSize;
    }

    /**
     * Sets default values for max mesh description
     */
    void set_default_values() {
        meshAttributes = 0x00000001;
        normalCount = 0;
        diffuseColorCount = 0;
        specularColorCount = 0;
        textureCoordCount = 0;
        shadingCount = 1;
        shadingDescription.resize( shadingCount );
        for( boost::uint32_t i = 0; i < shadingCount; ++i ) {
            shadingDescription[i].set_default_values( 0 );
        }
    }

    /**
     * Sets the number of shaders that will be used by this mesh
     */
    void set_shader_count( boost::uint32_t shaderCount ) {
        shadingCount = shaderCount;
        shadingDescription.resize( shadingCount );
        for( boost::uint32_t i = 0; i < shadingCount; ++i ) {
            shadingDescription[i].set_default_values( 0 );
        }
    }

    /**
     * Writes the data of max mesh description to a raw byte buffer
     *
     * @param buf The raw byte buffer to which the data must be written to.
     */
    void write( raw_byte_buffer& buf ) {
        write_value<boost::uint32_t>( buf, meshAttributes );
        write_value<boost::uint32_t>( buf, faceCount );
        write_value<boost::uint32_t>( buf, positionCount );
        write_value<boost::uint32_t>( buf, normalCount );
        write_value<boost::uint32_t>( buf, diffuseColorCount );
        write_value<boost::uint32_t>( buf, specularColorCount );
        write_value<boost::uint32_t>( buf, textureCoordCount );
        write_value<boost::uint32_t>( buf, shadingCount );
        for( boost::uint32_t i = 0; i < shadingCount; ++i ) {
            shadingDescription[i].write( buf );
        }
    }
};

/**
 * Object that contains the declaration information for a continuous level of detail mesh generator
 *
 * @var name Name of the CLOD mesh generator. Same as the resource modifier chain it is a part of
 * @var chainIndex position of the CLOD mesh generator in the model resource modifier chain. Always 0 for this blocktype
 * @var maxMeshDescription defines the size of the mesh at full resolution
 * @var minimumResolution Number of positions in the base mesh
 * @var finalMaximumResolution Number of positions in the Max Mesh Descriptions
 * @var positionQualityFactor Quality Factor associated with quantization of positions. Not needed for decoding
 * @var normalQualityFactor Quality Factor associated with quantization of normal vectors
 * @var textureCoordQualityFactor Quality Factor associated with  quantization of texture coordinates
 * @var positionInverseQuant Inverse Quantization factor used in the reconstruction of position vectors
 * @var normalInverseQuant Inverse Quantization factor used in the reconstruction of normal vectors
 * @var textureCoordInverseQuant Inverse Quantization factor used in the reconsturction of texture coordinates
 * @var diffuseColorInverseQuant Inverse Quantization factor used in the reconstruction of diffuse colors
 * @var specularColorInverseQuant Inverse Quantization factor used in the reconstruction of specular colors
 * @var normalCreaseParameter Read 9.6.1.1.5.3.1 of the documentation for details. We dont need this as we are not using
 normals with meshes.
 * @var normalUpdateParameter Read 9.6.1.1.5.3.2 of the documentation for details.Only for information. Not used during
 decoding. We dont need this as we are not using normals with meshes.
 * @var normalToleranceParameter Read 9.6.1.1.5.3.3 of the documentation for details. We dont need this as we are not
 using normals with meshes.
 * @var boneCount Number of bones associated with the mesh. This is always 0 in our case.
 *
 * @note If minimumResolution and finalMaximumResolutionare the same then the base mesh is the whole mesh and there is
 only
 *       1 level of detail.
 */
struct u3d_clod_mesh_declaration : u3d_block {
    frantic::tstring name;
    boost::uint32_t chainIndex;
    u3d_max_mesh_description maxMeshDescription;
    boost::uint32_t minimumResolution;
    boost::uint32_t finalMaximumResolution;
    boost::uint32_t positionQualityFactor;
    boost::uint32_t normalQualityFactor;
    boost::uint32_t textureCoordQualityFactor;
    float positionInverseQuant;
    float normalInverseQuant;
    float textureCoordInverseQuant;
    float diffuseColorInverseQuant;
    float specularColorInverseQuant;
    float normalCreaseParameter;
    float normalUpdateParameter;
    float normalToleranceParameter;
    boost::uint32_t boneCount;

    /**
     * Sets the default values for CLOD mesh declaration block
     */
    void set_default_values() {
        blockType = U3D_CLOD_MESH_DECLARATION_BLOCKTYPE;
        metaDataSize = 0;
        chainIndex = 0;
        maxMeshDescription.set_default_values();
        positionQualityFactor = 0;
        normalQualityFactor = 0;
        textureCoordQualityFactor = 0;
        positionInverseQuant = 1.0f;
        normalInverseQuant = 1.0f;
        textureCoordInverseQuant = 1.0f;
        diffuseColorInverseQuant = 1.0f;
        specularColorInverseQuant = 1.0f;
        normalCreaseParameter = 0.0f;
        normalUpdateParameter = 0.0f;
        normalToleranceParameter = 0.0f;
        boneCount = 0;
    }

    /**
     * Updates and returns the dataSize stored by this block
     */
    boost::uint32_t update_and_get_size() {
        dataSize = 60 + 2 + static_cast<boost::uint32_t>( name.length() ) +
                   static_cast<boost::uint32_t>( maxMeshDescription.size() );
        return dataSize;
    }

    /**
     * Writes the data of CLOD mesh declaration block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, name );
        write_value<boost::uint32_t>( buf, chainIndex );
        maxMeshDescription.write( buf );
        write_value<boost::uint32_t>( buf, minimumResolution );
        write_value<boost::uint32_t>( buf, finalMaximumResolution );
        write_value<boost::uint32_t>( buf, positionQualityFactor );
        write_value<boost::uint32_t>( buf, normalQualityFactor );
        write_value<boost::uint32_t>( buf, textureCoordQualityFactor );
        write_value<float>( buf, positionInverseQuant );
        write_value<float>( buf, normalInverseQuant );
        write_value<float>( buf, textureCoordInverseQuant );
        write_value<float>( buf, diffuseColorInverseQuant );
        write_value<float>( buf, specularColorInverseQuant );
        write_value<float>( buf, normalCreaseParameter );
        write_value<float>( buf, normalUpdateParameter );
        write_value<float>( buf, normalToleranceParameter );
        write_value<boost::uint32_t>( buf, boneCount );
        write_buffered_block( stream, buf, "Mesh Declaration" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

/**
 * Object that holds information about each face used in the mesh
 *
 * @var shadingID Index of the shader list descriptions used for this face
 * @var positionIndex Indices of all the vertices that desrcibe the face
 * @var normalIndex Indices of all the normals that descirbe the face. Absent if MaxMeshDescriptionn indicates Exclude
 * Normals.
 * @var diffuseColorIndex Indices of all the diffuseColors used by this face. Absent if shader list description doesn't
 * indicate it's presence.
 * @var specularColorIndex Indices of all specularColors used by this face. Absent if shader list description doesn't
 * indicate it's presence.
 * @var textureCoordIndex Indices of all textureCoords used by this face. In our case there is only 1 set of indices as
 * we have decided to have only 1 texture layer.
 */
struct u3d_base_mesh_face {
    boost::uint32_t shadingId;
    frantic::graphics::vector3 positionIndex;
    frantic::graphics::vector3 normalIndex;
    frantic::graphics::vector3 diffuseColorIndex;
    frantic::graphics::vector3 specularColorIndex;
    frantic::graphics::vector3 textureCoordIndex;
    bool hasTexture;
    bool hasColor;

    /**
     * Writes the data of a face of a mesh to a raw byte buffer
     *
     * @param buf The raw byte buffer to which the data must be written to.
     */
    void write( raw_byte_buffer& buf ) {
        write_value<boost::uint32_t>( buf, shadingId );
        if( hasColor ) {
            write_value<boost::uint32_t>( buf, positionIndex.x );
            write_value<boost::uint32_t>( buf, diffuseColorIndex.x );
            write_value<boost::uint32_t>( buf, positionIndex.y );
            write_value<boost::uint32_t>( buf, diffuseColorIndex.y );
            write_value<boost::uint32_t>( buf, positionIndex.z );
            write_value<boost::uint32_t>( buf, diffuseColorIndex.z );
        } else if( hasTexture ) {
            write_value<boost::uint32_t>( buf, positionIndex.x );
            write_value<boost::uint32_t>( buf, textureCoordIndex.x );
            write_value<boost::uint32_t>( buf, positionIndex.y );
            write_value<boost::uint32_t>( buf, textureCoordIndex.y );
            write_value<boost::uint32_t>( buf, positionIndex.z );
            write_value<boost::uint32_t>( buf, textureCoordIndex.z );
        } else {
            write_value<boost::uint32_t>( buf, positionIndex.x );
            write_value<boost::uint32_t>( buf, positionIndex.y );
            write_value<boost::uint32_t>( buf, positionIndex.z );
        }
    }
};

/**
 * Object which contains base mesh information for a continuouss level of detail mesh generator.
 * The base mesh is the minimum LOD mesh
 *
 * @var name Name of the CLOD mesh generator. Should be the same as the name of the resource modifier chain it is a part
 * of.
 * @var chainIndex Position of the CLOD mesh generator in the model resource modifier chain. It is 0 for this blocktype.
 * @var baseFaceCount Number of faces in the base mesh
 * @var basePositionCount Number of positions used by the base mesh in the position array.
 * @var baseNormalCount Number of normal used by the base mesh in the normal array
 * @var baseDiffuseColorCount Number of colors used by the base mesh in the diffuse color array.
 * @var baseSpecularColorCount Number of colors used by the base mesh in the specular color array.
 * @var baseTextureCoordCount Number of texture coodinates used by the base mesh in the texture coordinate array.
 * @var positions Array with all the positions( vertices )
 * @var normals Array with all the normals
 * @var diffuseColors Array with all the diffuse colors
 * @var specularColors Array with all the specular colors
 * @var textureCoords Array with all the 4D texture coordinates
 * @var faces A struct which holds all the information about the faces in the mesh
 */
struct u3d_clod_base_mesh_continuation : u3d_block {
  private:
    boost::uint32_t baseFaceCount;
    boost::uint32_t basePositionCount;
    boost::uint32_t baseNormalCount;
    boost::uint32_t baseDiffuseColorCount;
    boost::uint32_t baseSpecularColorCount;
    boost::uint32_t baseTextureCoordCount;

  public:
    frantic::tstring name;
    boost::uint32_t chainIndex;
    std::vector<frantic::graphics::vector3f> positions;
    std::vector<frantic::graphics::vector3f> normals;
    std::vector<frantic::graphics::color_rgba_f> diffuseColors;
    std::vector<frantic::graphics::color_rgba_f> specularColors;
    std::vector<frantic::graphics::vector4f> textureCoords;
    std::vector<u3d_base_mesh_face> faces;

    /**
     * Sets the default values for base mesh continuation block
     */
    void set_default_values() {
        blockType = U3D_CLOD_BASE_MESH_DATA_BLOCKTYPE;
        metaDataSize = 0;
        chainIndex = 0;
        baseNormalCount = 0;
        baseDiffuseColorCount = 0;
        baseSpecularColorCount = 0;
        baseTextureCoordCount = 0;
    }

    /**
     * Adds positions( vertices ) to the base mesh continuation block that will be written to the file
     *
     * @param mesh A mesh that contains all the positions
     */
    void set_positions( const mesh_interface* mesh ) {
        basePositionCount = static_cast<boost::uint32_t>( mesh->get_num_verts() );
        positions.resize( basePositionCount );
        for( size_t i = 0; i < basePositionCount; ++i ) {
            positions[i] = mesh->get_vert( i );
        }
    }

    /**
     * Adds faces to the base mesh continuation block that will be written to the file
     *
     * @param mesh A mesh that contains all the faces
     */
    void set_faces( const mesh_interface* mesh ) {
        baseFaceCount = static_cast<boost::uint32_t>( mesh->get_num_faces() );
        faces.resize( baseFaceCount );

        for( boost::uint32_t i = 0; i < baseFaceCount; ++i ) {
            u3d_base_mesh_face face;
            face.shadingId = 0;
            face.positionIndex.x = static_cast<boost::int32_t>( mesh->get_face_vert_index( i, 0 ) );
            face.positionIndex.y = static_cast<boost::int32_t>( mesh->get_face_vert_index( i, 1 ) );
            face.positionIndex.z = static_cast<boost::int32_t>( mesh->get_face_vert_index( i, 2 ) );
            faces[i] = face;
            faces[i].hasTexture = false;
            faces[i].hasColor = false;
        }
    }

    void set_vertex_colors( const mesh_interface* mesh ) {
        const frantic::geometry::mesh_channel* channel = mesh->get_vertex_channel( _T("Color") );
        baseDiffuseColorCount = static_cast<boost::uint32_t>( channel->get_num_elements() );
        mesh_channel_cvt<frantic::graphics::color3f> acc( channel );
        diffuseColors.resize( baseDiffuseColorCount );
        for( boost::uint32_t i = 0; i < baseDiffuseColorCount; ++i ) {
            color3f c( acc.get_value( i ) );
            diffuseColors[i].set_r( c.r );
            diffuseColors[i].set_g( c.g );
            diffuseColors[i].set_b( c.b );
            diffuseColors[i].set_a( 1.0f );
        }
        for( boost::uint32_t i = 0; i < baseFaceCount; ++i ) {
            faces[i].hasColor = true;
            faces[i].diffuseColorIndex.x = static_cast<boost::int32_t>( mesh->get_face_vert_index( i, 0 ) );
            faces[i].diffuseColorIndex.y = static_cast<boost::int32_t>( mesh->get_face_vert_index( i, 1 ) );
            faces[i].diffuseColorIndex.z = static_cast<boost::int32_t>( mesh->get_face_vert_index( i, 2 ) );
        }
    }

    /**
     * Sets the UV coordinates of the mesh and texture indices of the faces
     *
     * @param mesh Mesh that contains the UV coordinates and the texture indices for all the faces
     */
    void set_textures( const mesh_interface* mesh ) {
        const frantic::geometry::mesh_channel* channel = mesh->get_vertex_channel( _T("TextureCoord") );
        baseTextureCoordCount = static_cast<boost::uint32_t>( channel->get_num_elements() );
        mesh_channel_cvt<frantic::graphics2d::vector2f> acc( channel );
        textureCoords.resize( baseTextureCoordCount );
        for( size_t i = 0; i < baseTextureCoordCount; ++i ) {
            frantic::graphics2d::vector2f coord( acc.get_value( i ) );
            textureCoords[i].x = coord.x;
            textureCoords[i].y = coord.y;
            textureCoords[i].z = 0.0f;
            textureCoords[i].w = 0.0f;
        }

        for( size_t i = 0; i < baseFaceCount; ++i ) {
            faces[i].hasTexture = true;
            faces[i].textureCoordIndex.x = static_cast<boost::int32_t>( acc.get_fv_index( i, 0 ) );
            faces[i].textureCoordIndex.y = static_cast<boost::int32_t>( acc.get_fv_index( i, 1 ) );
            faces[i].textureCoordIndex.z = static_cast<boost::int32_t>( acc.get_fv_index( i, 2 ) );
        }
    }

    /**
     * Updates and returns the dataSize stored by this block
     */
    boost::uint32_t update_and_get_size() {
        dataSize = 28 + 12 * basePositionCount + 16 * baseFaceCount + 2 + static_cast<boost::uint32_t>( name.length() );
        if( baseTextureCoordCount > 0 ) {
            dataSize += 16 * baseTextureCoordCount + 12 * baseFaceCount;
        }
        if( baseDiffuseColorCount > 0 ) {
            dataSize += 16 * baseDiffuseColorCount + 12 * baseFaceCount;
        }
        return dataSize;
    }

    /**
     * Writes the data of CLOD base mesh continuation block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, name );
        write_value<boost::uint32_t>( buf, chainIndex );
        write_value<boost::uint32_t>( buf, baseFaceCount );
        write_value<boost::uint32_t>( buf, basePositionCount );
        write_value<boost::uint32_t>( buf, baseNormalCount );
        write_value<boost::uint32_t>( buf, baseDiffuseColorCount );
        write_value<boost::uint32_t>( buf, baseSpecularColorCount );
        write_value<boost::uint32_t>( buf, baseTextureCoordCount );

        for( boost::uint32_t i = 0; i < basePositionCount; ++i ) {
            write_vector3f( buf, positions[i] );
        }

        for( boost::uint32_t i = 0; i < baseDiffuseColorCount; ++i ) {
            write_color_rgba_f( buf, diffuseColors[i] );
        }

        for( boost::uint32_t i = 0; i < baseTextureCoordCount; ++i ) {
            write_vector4f( buf, textureCoords[i] );
        }

        for( boost::uint32_t i = 0; i < baseFaceCount; ++i ) {
            faces[i].write( buf );
        }

        write_buffered_block( stream, buf, "Mesh" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

} // namespace u3d
} // namespace geometry
} // namespace frantic
