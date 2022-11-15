// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/scoped_array.hpp>

#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/framebuffer.hpp>
#include <string>

// this file has the following class hierarchy:
//
// The class hierarchy is a little bit weird because originally it was just two classes. One for singleface, one for
// cubeface. However, we then added support for DTEX deep textures. At that point, the exr loading classes were split
// up, and started inheriting from a common EXR base. The DTEX classes have since been removed for licensing reasons, so
// it no longer exists in the code. The hierarchy remains as follows:
//
//     atten_saver
//         singleface_atten_saver
//	          singleface_atten_dtex_saver (has been removed due to licensing issues)
//		      singleface_atten_exr_saver*
//	      cubeface_atten_saver
//	          cubeface_atten_dtex_saver (has been removed due to licensing issues)
//		      cubeface_atten_exr_saver*
//
//  *exr saver classes also inheret from base_atten_exr_saver which provides common members/functions

namespace frantic {
namespace rendering {

class depthbuffer_singleface; // forward declare

/**
 * Abstract class has no interface, because of the generic interface between directional/cube lights in the lights code.
 * Normally this superclass is not useful. To use this top-level class, you must cast to either singleface_atten_saver
 * or cubeface_atten_saver because their interfaces are completely different.
 *
 */
class atten_saver {
  public:
    /**
     * Function to tells us what subclass it is.
     */
    virtual bool is_cubeface() const = 0;

    /**
     * Function to finalize the map and write it to disk.
     * All input to this function (ie filename(s)) is provided on construction
     */
    virtual void write_file() = 0;
};

/**
 *
 * Abstract creator/saver for single face attenuation maps.
 *
 */
class singleface_atten_saver : public atten_saver {
  protected:
    frantic::graphics2d::size2 m_mapDim;

  public:
    /**
     * Function to tells the base class (atten_saver) what subclass it is.
     */
    bool is_cubeface() const { return false; }

    /**
     * Gets the dimensions of the deep attenuation map
     */
    const frantic::graphics2d::size2& size() const { return m_mapDim; }

    /**
     * Adds an attenuation value at a given pixel for a provided z value.
     */
    virtual void add_sample( int x, int y, float z, const frantic::graphics::alpha3f& atten ) = 0;

    /**
     * Adds an attenuation value at a given pixel for a provided z value.
     * This is the bilinear extrapolation of the "add_sample" function.
     * The attenuation value will be added to the neighbouring four pixels.
     */
    void add_sample_bilinear( float x, float y, float z, const frantic::graphics::alpha3f& atten );
};

/**
 *
 * Abstract creator/saver classes for cube face attenuation maps.
 *
 */
class cubeface_atten_saver : public atten_saver {
  protected:
    int m_mapWidth;

  public:
    /**
     * Function to tells the base class (atten_saver) what subclass it is.
     */
    bool is_cubeface() const { return true; }

    /**
     * Gets the width of the deep attenuation cube map. Note that the dimensions are always square.
     */
    int width() const { return m_mapWidth; }

    /**
     * Adds an attenuation value at a given pixel and cubeface for a provided z value.
     */
    virtual void add_sample( int x, int y, int cubefaceIndex, float z, const frantic::graphics::alpha3f& atten ) = 0;

    /**
     * Adds an attenuation value at a given position. The position is the sample's coordinate relative to the light.
     * This is the bilinear extrapolation of the "get_sample" function.
     * The attenuation value will be added to the neighbouring four pixels.
     */
    void add_sample_bilinear( const frantic::graphics::vector3f& pos, const frantic::graphics::alpha3f& atten );
};

/**
 *
 * Common based class for exr savers. Provides data members and function class for exr savers.
 *
 */
class base_atten_exr_saver {
  public:
    virtual ~base_atten_exr_saver() {}

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
    std::string get_exr_layer_name( int i ) const;

    /**
     * Internal initalize function used by the constructors of both singleface_atten_exr_saver and
     * cubeface_atten_exr_saver
     */
    void reset( frantic::graphics2d::size2 mapDim, int numSamples, float spacing, bool exponentialGrowth );

    /**
     * Internal file loader function used by both singleface_atten_exr_saver and cubeface_atten_exr_saver
     */
    void save_exr_to_file( const frantic::tstring& depthMapPath, bool isCubeMap );

    /**
     * Internal data sampler function used by both singleface_atten_exr_saver and cubeface_atten_exr_saver
     */
    void add_sample_internal( int x, int y, int cubefaceIndex, float z, const frantic::graphics::alpha3f& atten );
};

/**
 *
 * Exr implementation for creating/saving single face attenuation maps.
 *
 */
class singleface_atten_exr_saver : public singleface_atten_saver, public base_atten_exr_saver {
  private:
    frantic::tstring m_outputFilename;

  public:
    /**
     * Creates an Exr single-faced saver object. Final output filename is provided on construction so that write_file
     * can be abstract.
     */
    singleface_atten_exr_saver( const frantic::tstring& outputExrFilename, frantic::graphics2d::size2 mapDim,
                                int numSamples, float spacing, bool exponentialGrowth );

    /**
     * Adds an attenuation value at a given pixel for a provided z value.
     */
    virtual void add_sample( int x, int y, float z, const frantic::graphics::alpha3f& atten );

    /**
     * Function to finalize the map and write it to disk.
     * The input to this function (single exr filename) is provided on construction.
     */
    virtual void write_file();
};

/**
 *
 * Exr implementation for creating/saving cube face attenuation maps.
 *
 */
class cubeface_atten_exr_saver : public cubeface_atten_saver, public base_atten_exr_saver {
  private:
    frantic::tstring m_outputFilename;

  public:
    /**
     * Creates an Exr cube-faced saver object. Final output filename is provided on construction so that write_file can
     * be abstract.
     */
    cubeface_atten_exr_saver( const frantic::tstring& outputExrFilename, int mapWidth, int numSamples, float spacing,
                              bool exponentialGrowth );

    /**
     * Adds an attenuation value at a given pixel and cubeface for a provided z value.
     */
    virtual void add_sample( int x, int y, int cubefaceIndex, float z, const frantic::graphics::alpha3f& atten );

    /**
     * Function to finalize the map and write it to disk.
     * Output is one Exr file of dimensions [N,N*6] which stores all six sides of the N-lengthed cube map.
     * The input to this function (single exr filename) is provided on construction.
     */
    virtual void write_file();
};

} // namespace rendering
} // namespace frantic
