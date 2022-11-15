// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

// A port of some features from the "SphericalDistort" module of Awake
// Originally written by Mark Wiebe and Greg Koreman
// Ported to FranticLibrary by Evan Spearman

#include <vector>

#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics2d/image_channel.hpp>

namespace frantic {
namespace graphics2d {

enum projection_type { perspective, latlong, cubeface, last_projection_type };

struct image_data {
    virtual int num_channels() const = 0;
    virtual frantic::graphics2d::size2 size() const = 0;

    virtual std::vector<frantic::graphics::color6f>
    get_pixel_at( const frantic::graphics2d::vector2& pixelCoord ) const = 0;

    virtual void set_pixel_at( const frantic::graphics2d::vector2& pixelCoord,
                               const std::vector<frantic::graphics::color6f>& data ) = 0;
};

// OpenImageIO compatible data structure for storing input and output images.
struct openimageio_data : public image_data {
    std::vector<float> m_data;
    int m_numChannels;
    frantic::graphics2d::size2 m_imageSize;

  public:
    openimageio_data( const frantic::graphics2d::size2& imageSize, int numChannels )
        : m_numChannels( numChannels )
        , m_imageSize( imageSize )
        , m_data( numChannels * imageSize.area() ) {}

    int num_channels() const { return m_numChannels; }
    frantic::graphics2d::size2 size() const { return m_imageSize; }

    std::vector<frantic::graphics::color6f> get_pixel_at( const frantic::graphics2d::vector2& pixelCoord ) const;
    void set_pixel_at( const frantic::graphics2d::vector2& pixelCoord,
                       const std::vector<frantic::graphics::color6f>& data );

  private:
    inline int get_pixel_index( const frantic::graphics2d::vector2& pixelCoord ) const;
};

// contains the rendered image data for one side of the cube face
struct multiimage_container {
    frantic::graphics2d::size2 size;
    int imageCount;
    std::vector<boost::shared_ptr<frantic::graphics::color6f[]>> images;
};

class multiimage_data : public image_data {
    multiimage_container m_data[6];
    frantic::graphics2d::size2 m_size;
    int m_numElementsFilled;

  public:
    // creates multiple images from the given data
    multiimage_data( const std::vector<multiimage_container>& data );
    // initializes one instance of multiimage_struct with the given size and imageCount
    multiimage_data( const frantic::graphics2d::size2& size, int imageCount );

    int num_channels() const {
        return m_data[0].imageCount * ( sizeof( frantic::graphics::color6f ) / sizeof( float ) );
    }
    frantic::graphics2d::size2 size() const { return m_size; }

    // NOTE: assumes the six images are equally sized and  are arranged in a vertical strip with index 0 being the top
    std::vector<frantic::graphics::color6f> get_pixel_at( const frantic::graphics2d::vector2& pixelCoord ) const;
    void set_pixel_at( const frantic::graphics2d::vector2& pixelCoord,
                       const std::vector<frantic::graphics::color6f>& data );

    std::vector<boost::shared_ptr<frantic::graphics::color6f[]>> get_images_at( int imageIndex ) const;
};

class spherical_distortion {
    // The float  parameters
    float m_offX, m_offY; // Location of center point control (normalized coordinates)

    projection_type m_inputType;
    projection_type m_outputType;

    frantic::graphics::cube_face::cube_face_mapping
        m_cubeFaceTypeIn; // Cube face type of input (0 = Vertical cross, 1 = Horizontal cross)
    frantic::graphics::cube_face::cube_face_mapping
        m_cubeFaceTypeOut; // Cube face type of output (0 = Vertical cross, 1 = Horizontal cross)

    float m_fieldOfViewIn;
    float m_fieldOfViewOut;

    frantic::graphics::transform4f m_rotationMatrix; // Rotation matrix for viewing direction

  public:
    spherical_distortion(
        projection_type inputType, projection_type outputType,
        const frantic::graphics::cube_face::cube_face_mapping& cubeFaceTypeIn =
            frantic::graphics::cube_face::STRIP_VERTICAL,
        const frantic::graphics::cube_face::cube_face_mapping& cubeFaceTypeOut =
            frantic::graphics::cube_face::STRIP_VERTICAL,
        float fieldOfViewIn = 360.f, float fieldOfViewOut = 360.f,
        const frantic::graphics::transform4f& rotationMatrix = frantic::graphics::transform4f::identity(),
        float offX = 0.f, float offY = 0.f );
    void do_distortion( const image_data& input, image_data& output );

  private:
    static void copy_pixel( const image_data& input, const frantic::graphics2d::vector2f& inputPixel,
                            image_data& output, const frantic::graphics2d::vector2& outputPixel );
    static void make_black( image_data& output, const frantic::graphics2d::vector2& outputPixel );
    inline static frantic::graphics2d::vector2 floor( const frantic::graphics2d::vector2f& in );
};
} // namespace graphics2d
} // namespace frantic
