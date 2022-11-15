// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/scoped_array.hpp>

#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/framebuffer.hpp>
#include <frantic/strings/tstring.hpp>

// this file has the following class hierarchy:
//
// The class hierarchy is a little bit weird because originally it was just two classes. One for singleface, one for
// cubeface. However, we then added support for DTEX deep textures. At that point, the exr loading classes were split
// up, and started inheriting from a common EXR base. The DTEX classes have since been removed for licensing reasons, so
// it no longer exists in the code. The hierarchy remains as follows:
//
//     atten_loader
//         singleface_atten_loader
//	          singleface_atten_dtex_loader (has been removed due to licensing issues)
//		      singleface_atten_exr_loader*
//	      cubeface_atten_loader
//	          cubeface_atten_dtex_loader (has been removed due to licensing issues)
//		      cubeface_atten_exr_loader*
//
//  *exr loader classes also inheret from base_atten_exr_loader which provides common members/functions

namespace frantic {
namespace rendering {

class depthbuffer_singleface; // forward declare

/**
 * Abstract class has no interface, because of the generic interface between directional/cube lights in the lights code.
 * Normally this superclass is not useful. To use this top-level class, you must cast to either singleface_atten_loader
 * or cubeface_atten_loader because their interfaces are completely different.
 *
 */
class atten_loader {
  public:
    /**
     * Function to tells us what subclass it is.
     */
    virtual bool is_cubeface() const = 0;
};

/**
 * Abstract sampler classes for single face attenuation maps.
 * Create these objects by calling create_singleface_atten_loader.
 *
 */
class singleface_atten_loader : public atten_loader {
  protected:
    frantic::graphics2d::size2 m_mapDim;

  public:
    /**
     * Function to tells the base class (atten_loader) what subclass it is.
     */
    bool is_cubeface() const { return false; }

    /**
     * Gets the dimensions of the deep attenuation map
     */
    const frantic::graphics2d::size2& size() const { return m_mapDim; }

    /**
     * Gets the attenuation value at a given pixel for a provided z value.
     */
    virtual frantic::graphics::alpha3f get_sample( int x, int y, float z ) const = 0;
    virtual frantic::graphics::alpha3f get_sample( int x, int y, float z, float& outImgZ ) const = 0;

    /**
     * Gets the Z depth of the first non-zero alpha at the given pixel. This is useful for creating a z-depth pass.
     */
    virtual float get_zdepth( int x, int y ) const = 0;

    /**
     * Gets the attenuation value at a given pixel for a provided z value.
     * This is the bilinear interpolation of the "get_sample" function.
     * The attenuation value will be sampled from the neighbouring four pixels.
     */
    frantic::graphics::alpha3f get_sample_bilinear( float x, float y, float z ) const;

    /**
     * Gets the Z depth of the first non-zero alpha at the given pixel.
     * This function is useful for creating a z-depth pass.
     * This is the bilinear interpolation of the "get_zdepth" function.
     * The returned Z value will be sampled from the neighbouring four pixels, and the nearest Z value will be returned.
     */
    float get_zdepth_bilinear( float x, float y ) const;
};

/**
 * Abstract sampler classes for cube face attenuation maps.
 * Create these objects by calling create_cubeface_atten_loader.
 *
 */
class cubeface_atten_loader : public atten_loader {
  protected:
    int m_mapWidth;

  public:
    /**
     * Function to tells the base class (atten_loader) what subclass it is.
     */
    bool is_cubeface() const { return true; }

    /**
     * Gets the width of the deep attenuation cube map. Cube maps are always square.
     */
    int width() const { return m_mapWidth; }

    /**
     * Gets the attenuation value at a given pixel and cubeface index for a provided z value.
     */
    virtual frantic::graphics::alpha3f get_sample( int x, int y, int cubefaceIndex, float z ) const = 0;
    virtual frantic::graphics::alpha3f get_sample( int x, int y, int cubefaceIndex, float z, float& outImgZ ) const = 0;

    /**
     * Gets the Z depth of the first non-zero alpha at the given pixel and cubeface index. This is useful for creating a
     * z-depth pass.
     */
    virtual float get_zdepth( int x, int y, int cubefaceIndex ) const = 0;

    /**
     * Gets the attenuation value at a given position. The position is the sampling coordinate relative to the light.
     * This is the bilinear interpolation of the "get_sample" function.
     * The attenuation value will be sampled from the neighbouring four pixels.
     */
    frantic::graphics::alpha3f get_sample_bilinear( const frantic::graphics::vector3f& pos ) const;

    /**
     * Gets the Z depth of the first non-zero alpha at the given position. The position is the sampling coordinate
     * relative to the light. This function is useful for creating a z-depth pass. This is the bilinear interpolation of
     * the "get_zdepth" function. The returned Z value will be sampled from the neighbouring four pixels, and the
     * nearest Z value will be returned.
     */
    float get_zdepth_bilinear( const frantic::graphics::vector3f& pos ) const;
};

/**
 *
 * Common based class for exr loaders. Provides data members and function class for exr loaders.
 *
 */
class base_atten_exr_loader {
  protected:
    // Stores the Z depth of the beginning of the sampled data along each ray. Samples in front of this buffer are NOT
    // occluded at all.
    boost::shared_ptr<frantic::rendering::depthbuffer_singleface> m_depthMap;

    // Store each sample for every pixel. m_attenBuffers[i](x,y) is the attenuation at m_depthMap(x,y) + (i+1) *
    // m_sampleSpacing for pixel (x,y)
    boost::scoped_array<frantic::graphics2d::framebuffer<frantic::graphics::alpha3f>> m_attenBuffers;

    // Store the number of attenuation samples per-pixel. This is the size of the 'm_attenBuffers' array.
    int m_sampleCount;

    // Store the worldspace z-distance between samples. This is NOT ray distance, but camera Z distance.
    float m_sampleSpacing;

    // If false, samples are evenly spaced. If true, sample spacings double. ie. The n'th sample is at sum{i=0 to n-1}(
    // 2^i * m_sampleSpacing ) = m_sampleSpacing * (2^n - 1)
    bool m_exponentialSampleSpacing;

  protected:
    /// Helper function
    std::string get_exr_layer_name( int i );

    /**
     * Internal file loader function used by both singleface_atten_exr_loader and cubeface_atten_exr_loader
     */
    void load_exr_from_file( const frantic::tstring& depthMapPath, bool isCubeMap );

    /**
     * Internal data sampler function used by both singleface_atten_exr_loader and cubeface_atten_exr_loader
     */
    frantic::graphics::alpha3f get_sample_internal( int x, int y, int cubefaceIndex, float z, float& outImgZ ) const;

  public:
    virtual ~base_atten_exr_loader() {}
};

/**
 * Exr implementation for loading/sampling single face attenuation maps.
 * Instead of creating one of these directly, consider using create_singleface_atten_loader.
 *
 */
class singleface_atten_exr_loader : public singleface_atten_loader, public base_atten_exr_loader {
  public:
    /**
     * Constructor provides the filename of the single faced exr attenuation map to be loaded.
     */
    singleface_atten_exr_loader( const frantic::tstring& filename );

    /**
     * Gets the attenuation value at a given pixel for a provided z value.
     */
    virtual frantic::graphics::alpha3f get_sample( int x, int y, float z ) const;
    virtual frantic::graphics::alpha3f get_sample( int x, int y, float z, float& outImgZ ) const;

    /**
     * Gets the Z depth of the first non-zero alpha at the given pixel. This is useful for creating a z-depth pass.
     */
    virtual float get_zdepth( int x, int y ) const;
};

/**
 * Exr implementation for loading/sampling cube face attenuation maps.
 * Instead of creating one of these directly, consider using create_cubeface_atten_loader.
 *
 */
class cubeface_atten_exr_loader : public cubeface_atten_loader, public base_atten_exr_loader {
  public:
    /**
     * Constructor provides the filename of the cube faced exr attenuation map to be loaded.
     */
    cubeface_atten_exr_loader( const frantic::tstring& filename );

    /**
     * Gets the attenuation value at a given pixel and cubeface index for a provided z value.
     */
    virtual frantic::graphics::alpha3f get_sample( int x, int y, int cubefaceIndex, float z ) const;
    virtual frantic::graphics::alpha3f get_sample( int x, int y, int cubefaceIndex, float z, float& outImgZ ) const;

    /**
     * Gets the Z depth of the first non-zero alpha at the given pixel and cubeface index. This is useful for creating a
     * z-depth pass.
     */
    virtual float get_zdepth( int x, int y, int cubefaceIndex ) const;
};

//
//
// convenience factory function to create an exr singleface/cubeface loader
//
//

// Use this function to create single faced attenuation loaders. Could be extended to load different file types, but
// currently only loads our custom layered EXR files.
boost::shared_ptr<singleface_atten_loader> create_singleface_atten_loader( const frantic::tstring& filename );

// Use this function to create cube faced attenuation loaders. Could be extended to load different file types, but
// currently only loads our custom layered EXR files.
boost::shared_ptr<cubeface_atten_loader> create_cubeface_atten_loader( const frantic::tstring& filename );

} // namespace rendering
} // namespace frantic
