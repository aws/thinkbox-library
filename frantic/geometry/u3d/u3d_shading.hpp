// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/filesystem/operations.hpp>
#include <frantic/geometry/u3d/u3d_common.hpp>

#ifdef OIIO_LIB_AVAILABLE
#include <OpenImageIO/imagebuf.h>
#endif

#define U3D_SHADING_MODIFIER_BLOCKTYPE 0XFFFFFF45
#define U3D_LIT_TEXTURE_SHADER_BLOCKTYPE 0xFFFFFF53
#define U3D_MATERIAL_RESOURCE_BLOCKTYPE 0xFFFFFF54
#define U3D_TEXTURE_RESOURCE_BLOCKTYPE 0xFFFFFF55
#define U3D_TEXTURE_CONTINUATION_BLOCKTYPE 0xFFFFFF5C

namespace frantic {
namespace geometry {
namespace u3d {

/**
 * Object that contains the information about the shading group that is used in the drawing of a renderable group
 *
 * @var name Name of the shading modifier. Same as the name of the modifier chain that contains this modifier.
 * @var chainIndex Position of this modifier in the modifier chain
 * @var shadingAttributes Collection of flags
 *      0x00000001 - Mesh: the shading group is applied to the renderable mesh group
 *      0x00000002 - Line: the shading group is applied to the renderable line group
 *      0x00000004 - Line: the shading group is applied to the renderable point group
 *      0x00000008 - Line: the shading group is applied to the gyph string
 * @var shaderListCount Number of shader lists in the shading group
 * @var shaderCount Vector that holds the number of shders in each shader list
 * @var shaderName Vector that holds the name of each shader in each shader list
 */
struct u3d_shading_modifier_block : u3d_block {
    frantic::tstring name;
    boost::uint32_t chainIndex;
    boost::uint32_t shadingAttributes;
    boost::uint32_t shaderListCount;
    std::vector<boost::uint32_t> shaderCount;
    std::vector<frantic::tstring> shaderName;

    /**
     * Sets the default values for this block
     */
    void set_default_values() {
        blockType = U3D_SHADING_MODIFIER_BLOCKTYPE;
        dataSize = 12;
        metaDataSize = 0;
        chainIndex = 1;
        shadingAttributes = 0x00000001;
        shaderListCount = 1;
        shaderCount.push_back( 1 );
        shaderName.resize( 1 );
    }

    /**
     * Updates and returns the dataSize stored by this block
     */
    boost::uint32_t update_and_get_size() {
        dataSize += 4 * static_cast<boost::uint32_t>( shaderCount.size() );
        dataSize += 2 + static_cast<boost::uint32_t>( name.size() );
        for( size_t i = 0; i < shaderName.size(); ++i ) {
            dataSize += 2 + static_cast<boost::uint32_t>( shaderName[i].size() );
        }
        // Include padding to dataSize (round up to a multiple of 4 bytes)
        dataSize = ( dataSize + 3 ) & ~0x03;
        return dataSize;
    }

    /**
     * Writes the data of shading modifier block to the file
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
        write_value<boost::uint32_t>( buf, shadingAttributes );
        write_value<boost::uint32_t>( buf, shaderListCount );
        write_value<boost::uint32_t>( buf, shaderCount[0] );
        write_tstring( buf, shaderName[0] );
        write_buffered_block( stream, buf, "Shading Modifier" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

/**
 * Holds the information about the texture used by a particular shader channel.
 *
 * @var name Name of the texture resource that is used for this texture layer
 * @var textureIntensity Scale factor applied to the color components of the texture
 * @var blendFunction Determines how the current texture layer is combined with the result from previous layers
 *       0 - Multiply: blended = current * previous
 *       1 - Add: blended = current + previous
 *       2 - Replace: blended = current
 *       3 - Blend: blended = current * currentAlpha + previous * ( 1 - currentAlpha )
 * @var blendSource Indicates whether the blending operation combines the current layer with the result from previous
 * layers
 *      0 - Alpha value of each pixel
 *      1 - Blending constant
 * @var blendConstant Used when combining the results of texture layers
 * @var textureMode Indicates the source of texture coordinates used to map the texture onto the model
 *      0x00: TM_NONE The shader does not generate texture coordinates
 *      0x01: TM_PLANAR The shader trasnforms the model by the inverse of the texture wrap transform and then performs a
 * planar x,y mapping of the texture onto the model
 *      0x02: TM_CYLINDRICAL The shader transofrms the model by the inverse of the texture wrap transform and then
 * performs a cylindrical mapping of the texture onto the model
 *            The Z-axis of the transformed model is the cylinder axis
 *      0x03: TM_SPHERICAL The shader transforms the model by the inverse of the texture wrap transform and then
 * performs a spherical mapping of the texture onto the model
 *            The Z-axis of the transformed model is the sphere's vertical axis
 *      0x04: TM_REFLECTION The shader performs a spherical reflection mapping. This is used to generate texture
 * coordinates for reflection mapping when using a specially designed
 *            spherical reflection texture.
 *      @note Only TM_NONE and TM_REFLECTION are supported by Acrobat
 * @var textureTransformMatrix Read section 9.8.3.10.7
 * @var textureWrapTransformMatrix Read section 9.8.3.10.8
 * @var textureRepeat Indicates whether or not the texture in the specified texture layer should be tiled beyond the
 * coordinate range.
 *      0x01 - Repeat in the direction of the first texture coordinate dimension
 *      0x02 - Repeat in the direction of the second texture coordinate dimension
 */
struct u3d_texture_information {
    frantic::tstring name;
    float textureIntensity;
    boost::uint8_t blendFunction;
    boost::uint8_t blendSource;
    float blendConstant;
    boost::uint8_t textureMode;
    frantic::graphics::transform4f textureTransformMatrix;
    frantic::graphics::transform4f textureWrapTransformMatrix;
    boost::uint8_t textureRepeat;

    /**
     * Sets the default values for this block
     */
    void set_default_values() {
        textureIntensity = 1.0f;
        blendFunction = 2;
        blendSource = 0;
        blendConstant = 0.5f;
        textureMode = 0x00;
        textureTransformMatrix.set_to_identity();
        textureWrapTransformMatrix.set_to_identity();
        textureRepeat = 0;
    }

    /**
     * Writes the data of texture information to a raw byte buffer
     *
     * @param buf The raw byte buffer to which the data must be written
     */
    void write( raw_byte_buffer& buf ) {
        write_tstring( buf, name );
        write_value<float>( buf, textureIntensity );
        write_value<boost::uint8_t>( buf, blendFunction );
        write_value<boost::uint8_t>( buf, blendSource );
        write_value<float>( buf, blendConstant );
        write_value<boost::uint8_t>( buf, textureMode );
        write_transform4f( buf, textureTransformMatrix );
        write_transform4f( buf, textureWrapTransformMatrix );
        write_value<boost::uint8_t>( buf, textureRepeat );
    }
};

/**
 * Contains the information needed to determine the appearance of a surface during rendering
 *
 * @var name Name that identifies this shader
 * @var textureAttributes Stores information about this shader. They are combined by a bitwise OR operation
 *      0x00000001: Lighting Enabled
 *      0x00000002: Alpha Test Enabled
 *      0x00000004: Use Vertex Color
 * @var alphaTestReference Value used in comparisons when alpha test is enabled
 * @var alphaTestFunction
 *      0x00000610: NEVER: The test never passes. No pixels are drawn
 *      0x00000611: LESS: The rendered alpha value must be less than the reference value
 *      0x00000612: GREATER: The rendered alpha value must be greater than the reference value
 *      0x00000613: EQUAL: The rendered alpha value must be equal to the reference value
 *      0x00000614: NOT_EQUAL: The rendered alpha value must not be equal to the reference value.
 *      0x00000615: LEQUAL: The rendered alpha value must be less than or equal to the reference value.
 *      0x00000616: GEQUAL: The rendered alpha value must be greater than or equal to the reference value.
 *      0x00000617: ALWAYS: The test always passes. No pixels are rejected.
 * @var colorBlendFunction Functiopn used to blend rendered pixels and the existing frame buffer
 *      0x00000604: FB_ADD: Add the RGB components into the framebuffer
 *      0x00000605: FB_MULTIPLY: Multiply the RGB components into the framebuffer
 *      0x00000606: FB_ALPHA_BLEND: Linear blend the RGB components into the framebuffer based on the rendered alpha
 *                                  value.
 *      0x00000607: FB_INV_ALPHA_BLEND: Linear blend the RGB components into framebuffer based on the inverse (1.0 - a)
 *                                      of the rendered alpha.
 * @var renderPassEnabledFlags Determines which passes this shader uses. Each bit (1<<n) in the flags determines if the
 *                             shader is used in pass n.
 * @var shaderChannels Bit field that determines which of the model's texture coordinate layers are used for this
 *                     shader. The least significant 8 bits are used
 *                     A layer is active if the corresponding bit is set.
 * @var alphaTextureChannels Bit field that determines which texture layers should use alpha component if an alpha
 *                           component exists. The least significant 8 bits are used
 *                           The Alpha Texture Channel bit shall not be set if the corresponding Shader Channel bit is
 *                           not set.
 * @var materialName Name of the material associated with this shader that determines how the shader appears when lit
 * @var textureInformation Identifies the texture used by a particular shader channel
 */
struct u3d_lit_texture_shader_block : u3d_block {
    frantic::tstring name;
    boost::uint32_t textureAttributes;
    float alphaTestReference;
    boost::uint32_t alphaTestFunction;
    boost::uint32_t colorBlendFunction;
    boost::uint32_t renderPassEnabledFlags;
    boost::uint32_t shaderChannels;
    boost::uint32_t alphaTextureChannels;
    frantic::tstring materialName;
    u3d_texture_information textureInformation;

    /**
     * Sets the default values for this block
     */
    void set_default_values() {
        blockType = U3D_LIT_TEXTURE_SHADER_BLOCKTYPE;
        metaDataSize = 0;
        dataSize = 168;
        textureAttributes = 0X00000001;
        alphaTestReference = 0.0f;
        alphaTestFunction = 0x00000617;
        colorBlendFunction = 0x00000606;
        renderPassEnabledFlags = 0X00000001;
        shaderChannels = 0x00000001;
        alphaTextureChannels = 0x00000000;
        textureInformation.set_default_values();
    }

    /**
     * Writes the data of texture shader block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, name );
        write_value<boost::uint32_t>( buf, textureAttributes );
        write_value<float>( buf, alphaTestReference );
        write_value<boost::uint32_t>( buf, alphaTestFunction );
        write_value<boost::uint32_t>( buf, colorBlendFunction );
        write_value<boost::uint32_t>( buf, renderPassEnabledFlags );
        write_value<boost::uint32_t>( buf, shaderChannels );
        write_value<boost::uint32_t>( buf, alphaTextureChannels );
        write_tstring( buf, materialName );
        textureInformation.write( buf );
        write_buffered_block( stream, buf, "Texture Shader" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

/**
 * Object that holds the information about how a material interacts with light in a scene.
 * @note A shader references a Material Resource to determine how surfaces will appear when rendered
 *
 * @var name Name of the material
 * @var materialAttributes Collection of flags that define which of the material attributes specified are enabled
 *      0x00000001 - Ambient
 *      0x00000002 - Diffuse
 *      0x00000004 - Specular
 *      0x00000008 - Emissive
 *      0x00000010 - Reflectivity
 *      0x00000020 - Opacity
 * @var ambientColor Defines the material's appearance in ambient light
 * @var diffuseColor Defines the material's appearance in diffuse light
 * @var specularColor Defines the material's appearance in specular light
 * @var emissiveColor Defines the light that material appears to give off
 * @var reflectivity Measures how shiny a material appears to be
 * @var opacity Measure of an object's transparency
 */
struct u3d_material_resource_block : u3d_block {
    frantic::tstring name;
    boost::uint32_t materialAttributes;
    frantic::graphics::color3f ambientColor;
    frantic::graphics::color3f diffuseColor;
    frantic::graphics::color3f specularColor;
    frantic::graphics::color3f emissiveColor;
    float reflectivity;
    float opacity;

    /**
     * Sets the default values for this block
     */
    void set_default_values() {
        blockType = U3D_MATERIAL_RESOURCE_BLOCKTYPE;
        metaDataSize = 0;
        dataSize = 60;
        diffuseColor = frantic::graphics::color3f::white();
        ambientColor = frantic::graphics::color3f::white();
        specularColor = frantic::graphics::color3f::white();
        emissiveColor = frantic::graphics::color3f::white();
        reflectivity = 0.0f;
        opacity = 1.0f;
        materialAttributes = 0x0000003F;
    }

    /**
     * Writes the data of material resource block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, name );
        write_value<boost::uint32_t>( buf, materialAttributes );
        write_color3f( buf, ambientColor );
        write_color3f( buf, diffuseColor );
        write_color3f( buf, specularColor );
        write_color3f( buf, emissiveColor );
        write_value<float>( buf, reflectivity );
        write_value<float>( buf, opacity );
        write_buffered_block( stream, buf, "Material Resource" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

/**
 * Object that contains information about the continuation image in the texture declaration
 *
 * @var compressionType Defines the scheme used to compress the Image Data in the texture continuation blocks
 *      0x01 - JPEG-24 (color, baseline profile)
 *      0x02 - PNG
 *      0x03 - JPEG-8 (greyscale, baseline profile)
 *      0x04 - TIFF
 *      @note TIFF is not supported by Acrobat
 * @var textureImageChannels Indicates which color channels of the texture image are composed using this image
 *      0x01: alpha channel
 *      0x02: blue channel
 *      0x04: green channel
 *      0x08: red channel
 *      0x10: luminance
 * @var continuationImageAttributes Contains additional information about the continuation image
 *      0x0000: default attributes
 *      0x0001: external continuation image file reference
 * @var imageDataByteCount sum of number of bytes of Image Data in all continuation blocks for this continuation image
 */
struct u3d_continuation_image_format {
    boost::uint8_t compressionType;
    boost::uint8_t textureImageChannels;
    boost::uint16_t continuationImageAttributes;
    boost::uint32_t imageDataByteCount;

    /**
     * Writes the data of the continuation image format to a raw byte buffer
     *
     * @param buf The raw byte buffer to which the data must be written to.
     */
    void write( raw_byte_buffer& buf ) {
        write_value<boost::uint8_t>( buf, compressionType );
        write_value<boost::uint8_t>( buf, textureImageChannels );
        write_value<boost::uint16_t>( buf, continuationImageAttributes );
        write_value<boost::uint32_t>( buf, imageDataByteCount );
    }
};

/**
 * Contains information for creating a texture image to be applied to geometry
 *
 * @var name Name used to identify the texture
 * @var textureHeight Height of texture in pixels
 * @var textureWidth Width of texture in pixels
 * @var textureImageType Identifies the color channels present in the texture image
 *      0x01 - alpha component
 *      0x0E - color RGB
 *      0x0F - color RGBA
 *      0x10 - luminance (greyscale)
 *      0x11 - luminance (greyscale and alpha)
 * @var continuationImageCount Number of images used to compose the texture image
 * @var continuationImageFormats Information about the continuation images used to make the texture declaration
 */
struct u3d_texture_declaration_block : u3d_block {
    frantic::tstring name;
    boost::uint32_t textureHeight;
    boost::uint32_t textureWidth;
    boost::uint8_t textureImageType;
    boost::uint32_t continuationImageCount;
    std::vector<u3d_continuation_image_format> continuationImageFormats;

    /**
     * Sets the default values for this block
     */
    void set_default_values() {
        blockType = U3D_TEXTURE_RESOURCE_BLOCKTYPE;
        dataSize = 21;
        metaDataSize = 0;
        continuationImageCount = 1;
        continuationImageFormats.resize( 1 );
        continuationImageFormats[0].continuationImageAttributes = 0;
    }

#ifdef OIIO_LIB_AVAILABLE
    /**
     * Sets the metadata about the image which will be embedded along with the mesh in the u3d file
     *
     * @param texture Image Buffer that contains the image
     */
    void set_image_metadata( const OIIO::ImageBuf* texture ) {
        OIIO::ImageSpec spec = texture->spec();
        textureHeight = spec.height;
        textureWidth = spec.width;
        continuationImageFormats[0].compressionType = 0x01;

        for( int i = 0; i < spec.nchannels; ++i ) {
            if( spec.channelnames[i] == "R" || spec.channelnames[i] == "G" || spec.channelnames[i] == "B" ) {
                textureImageType = 0x0E;
                continuationImageFormats[0].textureImageChannels = 0x02 | 0x04 | 0x08;
                break;
            } else {
                textureImageType = 0x10;
                continuationImageFormats[0].textureImageChannels = 0x10;
            }
        }

        if( spec.alpha_channel != -1 ) {
            textureImageType += 0x01;
            continuationImageFormats[0].textureImageChannels = continuationImageFormats[0].textureImageChannels | 0x01;
        }

        if( textureImageType == 0x10 || textureImageType == 0x11 ) {
            continuationImageFormats[0].compressionType += 0x02;
        }
    }
#endif

    /**
     * Writes the data of texture declaration block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, name );
        write_value<boost::uint32_t>( buf, textureHeight );
        write_value<boost::uint32_t>( buf, textureWidth );
        write_value<boost::uint8_t>( buf, textureImageType );
        write_value<boost::uint32_t>( buf, continuationImageCount );
        for( boost::uint32_t i = 0; i < continuationImageCount; ++i ) {
            continuationImageFormats[i].write( buf );
        }
        write_buffered_block( stream, buf, "Texture Declaration" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};

/**
 * Contains image data for a continuation image previously described in the texture declaration
 *
 * @var name Name of the texture resource with which this continuation block is associated
 * @var continuationImageIndex Index into the sequence of continuation image formats
 * @var imageData Format of the image data is indicated by Compression Type in the texture declaration
 */
struct u3d_texture_continuation_block : u3d_block {
    frantic::tstring name;
    boost::uint32_t continuationImageIndex;
    std::vector<char> imageData;

    /**
     * Sets the default values for this node
     */
    void set_default_values() {
        blockType = U3D_TEXTURE_CONTINUATION_BLOCKTYPE;
        dataSize = 4;
        metaDataSize = 0;
        continuationImageIndex = 0;
    }

#ifdef OIIO_LIB_AVAILABLE
    /**
     * Sets the imageData which will be embedded along with the mesh in the u3d file
     *
     * @param texture Image Buffer that contains the image
     */
    void set_image_data( const OIIO::ImageBuf* texture ) {
        namespace fs = boost::filesystem;
        fs::path tempFile = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.jpeg" );
        bool success = texture->write( tempFile.string(), "jpeg" );
        if( !success ) {
            throw std::runtime_error(
                "u3d_texture_continuation_block::set_image_data Error: Failed to write temporary jpeg file " +
                tempFile.string() + " with the texture" );
        }
        std::ifstream temp( tempFile.string().c_str(), std::ios::binary );
        if( temp.fail() ) {
            throw std::runtime_error( "u3d_texture_continuation_block::set_image_data Error: Temporary file " +
                                      tempFile.string() + " failed to open for reading" );
        }
        temp.seekg( 0, std::ios::end );
        boost::uint32_t size = (boost::uint32_t)temp.tellg();
        temp.seekg( 0, std::ios::beg );
        imageData.resize( size );
        temp.read( &imageData[0], size );
        temp.close();
        fs::remove( tempFile );

        dataSize += size;
    }
#endif

    /**
     * Writes the data of texture continuation block to the file
     *
     * @param stream The file to which the data must be written to.
     */
    void write( std::ostream& stream ) {
        raw_byte_buffer buf;
        write_value<boost::uint32_t>( buf, blockType );
        write_value<boost::uint32_t>( buf, dataSize );
        write_value<boost::uint32_t>( buf, metaDataSize );
        write_tstring( buf, name );
        write_value<boost::uint32_t>( buf, continuationImageIndex );
        memcpy( buf.add_element( imageData.size() ), &imageData[0], imageData.size() );
        write_buffered_block( stream, buf, "Texture" );
        add_padding( stream, get_alignment_data_padding( static_cast<boost::uint32_t>( stream.tellp() ), 4u ) );
    }
};
} // namespace u3d
} // namespace geometry
} // namespace frantic
