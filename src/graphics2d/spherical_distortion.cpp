// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics2d/spherical_distortion.hpp>
#include <frantic/logging/logging_level.hpp>

using namespace frantic::graphics;
using namespace frantic::graphics2d;

spherical_distortion::spherical_distortion( projection_type inputType, projection_type outputType,
                                            const frantic::graphics::cube_face::cube_face_mapping& cubeFaceTypeIn,
                                            const frantic::graphics::cube_face::cube_face_mapping& cubeFaceTypeOut,
                                            float fieldOfViewIn, float fieldOfViewOut,
                                            const frantic::graphics::transform4f& rotationMatrix, float offX,
                                            float offY )
    : m_offX( offX )
    , m_offY( offY )
    , m_inputType( inputType )
    , m_outputType( outputType )
    , m_cubeFaceTypeIn( cubeFaceTypeIn )
    , m_cubeFaceTypeOut( cubeFaceTypeOut )
    , m_fieldOfViewIn( fieldOfViewIn )
    , m_fieldOfViewOut( fieldOfViewOut )
    , m_rotationMatrix( rotationMatrix ) {}

void spherical_distortion::do_distortion( const image_data& input, image_data& output ) {
    if( input.num_channels() != output.num_channels() ) {
        throw std::runtime_error( "Input and Output must have the same number of channels." );
    }
    for( int x = 0; x < output.size().xsize; ++x ) {
        for( int y = 0; y < output.size().ysize; ++y ) {
            const vector2 outCoord = vector2( x, y );
            bool validInputDirection = true;
            const size2 outSize( output.size().xsize, output.size().ysize );
            frantic::graphics::vector3f viewDirection;
            switch( m_outputType ) {
            case perspective: // Perspective
                viewDirection = vector3f::from_perspective_projection(
                    outCoord, outSize, static_cast<float>( tan( m_fieldOfViewOut / 2.f ) ) );
                break;
            case latlong: // Lat Long
                viewDirection = vector3f::from_longlat_yup( outCoord, outSize );
                break;
            case cubeface: // Cube Face
            {
                viewDirection = from_cube_face_map( outCoord, outSize, m_cubeFaceTypeOut, validInputDirection );
            } break;
            default: {
                throw std::runtime_error( "Error: spherical_distortion: invalid output projection." );
            }
            }

            ///////////////////////////////////
            // Step 2: Apply the rotations
            ///////////////////////////////////
            viewDirection = m_rotationMatrix * viewDirection;

            ///////////////////////////////////
            // Step 3: Convert to source pixel coordinates
            ///////////////////////////////////

            size2f inSize( static_cast<float>( input.size().xsize ), static_cast<float>( input.size().ysize ) );

            // The 2D position on the IN image
            vector2f srcpos;

            switch( m_inputType ) {
            case perspective: // Perspective
            {
                viewDirection.normalize();
                srcpos = viewDirection.to_perspective_projection(
                    inSize, static_cast<float>( tan( m_fieldOfViewIn / 2.f ) ), validInputDirection );

                if( validInputDirection ) {
                    srcpos.x -= ( 0.5f - m_offX ) * static_cast<float>( input.size().xsize - 1 );
                    srcpos.y -= ( 0.5f - m_offY ) * static_cast<float>( input.size().ysize - 1 );

                    // Get colour
                    copy_pixel( input, srcpos, output, outCoord );
                } else {
                    make_black( output, outCoord );
                }
            } break;
            case latlong: // Lat Long
            {
                srcpos = viewDirection.to_longlat_yup( inSize );
                // Get colour
                copy_pixel( input, srcpos, output, outCoord );
            } break;
            case cubeface: // Cube Face
            {
                srcpos = to_cube_face_map( viewDirection, inSize, m_cubeFaceTypeIn, 0 );
                copy_pixel( input, srcpos, output, outCoord );
            } break;
            default: {
                throw std::runtime_error( "Error: spherical_distortion: invalid input projection." );
            }
            }
        }
    }
}

void spherical_distortion::copy_pixel( const image_data& input, const frantic::graphics2d::vector2f& inputPixel,
                                       image_data& output, const frantic::graphics2d::vector2& outputPixel ) {
    std::vector<color6f> pixelData = input.get_pixel_at( floor( inputPixel ) );
    output.set_pixel_at( outputPixel, pixelData );
}

void spherical_distortion::make_black( image_data& output, const frantic::graphics2d::vector2& outputPixel ) {
    std::vector<color6f> pixelData;
    static const color6f blackPixel;

    const int totalPixelCount = output.num_channels() / ( sizeof( color6f ) / sizeof( float ) );
    for( int i = 0; i < totalPixelCount; ++i ) {
        pixelData.push_back( blackPixel );
    }

    output.set_pixel_at( outputPixel, pixelData );
}

vector2 spherical_distortion::floor( const vector2f& in ) {
    return vector2( static_cast<int>( in.x ), static_cast<int>( in.y ) );
}

std::vector<color6f> openimageio_data::get_pixel_at( const frantic::graphics2d::vector2& pixelCoord ) const {
    std::vector<color6f> result;

    if( pixelCoord.x <= m_imageSize.xsize && pixelCoord.y <= m_imageSize.ysize ) {
        int pixelIndex = get_pixel_index( pixelCoord );
        const int totalIterations =
            static_cast<int>( ceil( static_cast<float>( m_numChannels ) / ( sizeof( color6f ) / sizeof( float ) ) ) );

        for( int iterationNum = 0; iterationNum < totalIterations; ++iterationNum ) {
            std::vector<float> colours( sizeof( color6f ) / sizeof( float ) );
            for( int i = 0; i < colours.size(); ++i ) {
                colours[i] = m_data[pixelIndex++];
            }

            color6f pixelData;
            pixelData.c.r = colours[0];
            pixelData.c.g = colours[1];
            pixelData.c.b = colours[2];
            pixelData.a.ar = colours[3];
            pixelData.a.ag = colours[4];
            pixelData.a.ab = colours[5];
            result.push_back( pixelData );
        }
    }

    return result;
}

void openimageio_data::set_pixel_at( const frantic::graphics2d::vector2& pixelCoord,
                                     const std::vector<color6f>& data ) {
    if( pixelCoord.x <= m_imageSize.xsize && pixelCoord.y <= m_imageSize.ysize ) {
        int pixelIndex = get_pixel_index( pixelCoord );
        static const int pixelChannels = sizeof( color6f ) / sizeof( float );
        const int channelsPerPixel = m_numChannels < pixelChannels ? m_numChannels : pixelChannels;

        for( std::vector<color6f>::const_iterator it = data.begin(); it != data.end(); ++it ) {
            std::vector<float> colours( pixelChannels );
            colours[0] = it->c.r;
            colours[1] = it->c.g;
            colours[2] = it->c.b;
            colours[3] = it->a.ar;
            colours[4] = it->a.ag;
            colours[5] = it->a.ab;

            for( int i = 0; i < channelsPerPixel; ++i ) {
                m_data[pixelIndex++] = colours[i];
            }
        }
    }
}

inline int openimageio_data::get_pixel_index( const frantic::graphics2d::vector2& pixelCoord ) const {
    return ( pixelCoord.x + pixelCoord.y * m_imageSize.xsize ) * m_numChannels;
}

multiimage_data::multiimage_data( const std::vector<multiimage_container>& data ) {
    const int numElements = sizeof( m_data ) / sizeof( multiimage_container );
    const int minimum = numElements < data.size() ? numElements : static_cast<int>( data.size() );
    for( int i = 0; i < minimum; ++i ) {
        m_data[i] = data[i];
    }

    m_size = frantic::graphics2d::size2( m_data[0].size.xsize, m_data[0].size.ysize * minimum );
    m_numElementsFilled = minimum;
}

multiimage_data::multiimage_data( const frantic::graphics2d::size2& size, int imageCount ) {
    m_size = size;

    m_data[0].size = size;
    m_data[0].imageCount = imageCount;

    for( int i = 0; i < imageCount; ++i ) {
        m_data[0].images.push_back( boost::shared_ptr<color6f[]>( new color6f[size.xsize * size.ysize] ) );
    }

    m_numElementsFilled = 1;
}

std::vector<color6f> multiimage_data::get_pixel_at( const frantic::graphics2d::vector2& pixelCoord ) const {
    std::vector<color6f> result;

    if( pixelCoord.x <= m_size.xsize && pixelCoord.y <= m_size.ysize ) {
        const int imageIndex = pixelCoord.y / m_data[0].size.ysize;
        const int y = pixelCoord.y % m_data[0].size.ysize;

        for( int i = 0; i < m_data[imageIndex].imageCount; ++i ) {
            color6f pixelData = m_data[imageIndex].images[i][pixelCoord.x + y * m_data[imageIndex].size.xsize];
            result.push_back( pixelData );
        }
    }

    return result;
}

void multiimage_data::set_pixel_at( const frantic::graphics2d::vector2& pixelCoord, const std::vector<color6f>& data ) {
    if( pixelCoord.x <= m_size.xsize && pixelCoord.y <= m_size.ysize && data.size() > 0 ) {
        const int imageIndex = pixelCoord.y / m_data[0].size.ysize;
        const int y = pixelCoord.y % m_data[0].size.ysize;

        for( int i = 0; i < m_data[imageIndex].imageCount; ++i ) {
            m_data[imageIndex].images[i][pixelCoord.x + y * m_data[imageIndex].size.xsize] = data[i];
        }
    }
}

std::vector<boost::shared_ptr<color6f[]>> multiimage_data::get_images_at( int imageIndex ) const {
    if( imageIndex < m_numElementsFilled ) {
        return m_data[imageIndex].images;
    } else {
        return std::vector<boost::shared_ptr<color6f[]>>();
    }
}
