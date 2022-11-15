// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/rendering/deep_attenuation_savers.hpp>

#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/framebuffer.hpp>
#include <frantic/rendering/depthbuffer_singleface.hpp>
#include <string>

#pragma warning( push, 3 )
#include <ImathBox.h>
#include <ImfChannelList.h>
#include <ImfEnvmap.h>
#include <ImfEnvmapAttribute.h>
#include <ImfFloatAttribute.h>
#include <ImfFrameBuffer.h>
#include <ImfHeader.h>
#include <ImfInputFile.h>
#include <ImfIntAttribute.h>
#include <ImfOutputFile.h>
#include <ImfStandardAttributes.h>
#include <half.h>
#pragma warning( pop )

using namespace frantic::graphics;

namespace frantic {
namespace rendering {

//
//
// abstract writer classes for single face attenuation maps
//
//

void singleface_atten_saver::add_sample_bilinear( float x, float y, float z, const frantic::graphics::alpha3f& atten ) {
    float fpixelX = floor( x - 0.5f );
    float fpixelY = floor( y - 0.5f );
    float alphaX = x - fpixelX - 0.5f;
    float alphaY = y - fpixelY - 0.5f;
    int pixelX = (int)fpixelX;
    int pixelY = (int)fpixelY;

    // TODO: This might be a performance issue, since we are skipping through multiple images for each pixel. Instead it
    // might be better
    //        to have the calculations for all 4 pixels done then skip through the image layers.
    if( (unsigned)pixelY < (unsigned)m_mapDim.ysize ) {
        if( (unsigned)pixelX < (unsigned)m_mapDim.xsize ) {
            float weight = alphaX * ( 1.f - alphaY );
            add_sample( pixelX, pixelY, z, atten * weight );
        }
        if( (unsigned)( pixelX + 1 ) < (unsigned)m_mapDim.xsize ) {
            float weight = alphaX * ( 1.f - alphaY );
            add_sample( pixelX + 1, pixelY, z, atten * weight );
        }
    }
    if( (unsigned)( pixelY + 1 ) < (unsigned)m_mapDim.ysize ) {
        if( (unsigned)pixelX < (unsigned)m_mapDim.xsize ) {
            float weight = ( 1.f - alphaX ) * alphaY;
            add_sample( pixelX, pixelY + 1, z, atten * weight );
        }
        if( (unsigned)( pixelX + 1 ) < (unsigned)m_mapDim.xsize ) {
            float weight = alphaX * alphaY;
            add_sample( pixelX + 1, pixelY + 1, z, atten * weight );
        }
    }
}

//
//
// abstract writer classes for cube face attenuation maps
//
//

void cubeface_atten_saver::add_sample_bilinear( const frantic::graphics::vector3f& pos,
                                                const frantic::graphics::alpha3f& atten ) {

    frantic::graphics::cube_face::default_cube_face cubeFaces[6];
    int numFaces = frantic::graphics::get_cube_faces( pos, (float)( m_mapWidth - 2 ) / (float)m_mapWidth, cubeFaces );

    frantic::graphics2d::vector2f coord = get_cube_face_coordinate( pos, cubeFaces[0] );

    // See framebuffer_cubeface for justification of this.
    float compensationFactor = ( coord.x * coord.x + coord.y * coord.y + 1 );
    compensationFactor *= sqrtf( compensationFactor );

    frantic::graphics::alpha3f compAtten = frantic::rendering::cubeface_compensation( atten, compensationFactor );
    const float dist = pos.get_magnitude();

    for( int i = 0; i < numFaces; ++i ) {
        coord = frantic::graphics::get_cube_face_coordinate( pos, cubeFaces[i] );
        coord += frantic::graphics2d::vector2f( 1, 1 );
        coord *= 0.5f * m_mapWidth;

        float floorX = floor( coord.x );
        float alphaX = coord.x - floorX;
        int x = (int)floorX;

        float floorY = floor( coord.y );
        float alphaY = coord.y - floorY;
        int y = (int)floorY;

        if( (unsigned)y < (unsigned)m_mapWidth ) {
            if( (unsigned)x < (unsigned)m_mapWidth ) {
                float weight = ( 1.f - alphaX ) * ( 1.f - alphaY );
                add_sample( x, y, cubeFaces[i], dist, compAtten * weight );
            }
            if( (unsigned)( x + 1 ) < (unsigned)m_mapWidth ) {
                float weight = alphaX * ( 1.f - alphaY );
                add_sample( x + 1, y, cubeFaces[i], dist, compAtten * weight );
            }
        }
        if( (unsigned)( y + 1 ) < (unsigned)m_mapWidth ) {
            if( (unsigned)x < (unsigned)m_mapWidth ) {
                float weight = ( 1.f - alphaX ) * alphaY;
                add_sample( x, y + 1, cubeFaces[i], dist, compAtten * weight );
            }
            if( (unsigned)( x + 1 ) < (unsigned)m_mapWidth ) {
                float weight = alphaX * alphaY;
                add_sample( x + 1, y + 1, cubeFaces[i], dist, compAtten * weight );
            }
        }
    }
}

//
//
// common data members and function class for exr savers
//
//

std::string base_atten_exr_saver::get_exr_layer_name( int i ) const {
    return "sample" + frantic::strings::to_string( frantic::strings::zero_pad( i, 4 ) );
}

void base_atten_exr_saver::reset( frantic::graphics2d::size2 mapDim, int numSamples, float spacing,
                                  bool exponentialGrowth ) {
    m_depthMap.reset( new frantic::rendering::depthbuffer_singleface( mapDim ) );
    m_depthMap->clear();

    m_sampleCount = numSamples;
    m_attenBuffers.reset( new frantic::graphics2d::framebuffer<frantic::graphics::alpha3f>[m_sampleCount] );
    for( int i = 0; i < m_sampleCount; ++i ) {
        m_attenBuffers[i].set_size( mapDim );
        m_attenBuffers[i].fill( frantic::graphics::alpha3f( 0 ) );
    }

    m_sampleSpacing = spacing;
    m_exponentialSampleSpacing = exponentialGrowth;
}

void base_atten_exr_saver::save_exr_to_file( const frantic::tstring& depthMapPath, bool isCubeMap ) {
    frantic::graphics2d::size2 imgSize = m_depthMap->size();

    Imf::Header fileHeader( imgSize.xsize, imgSize.ysize );
    fileHeader.insert( "attenSampleCount", Imf::IntAttribute( m_sampleCount ) );
    fileHeader.insert( "attenSampleSpacing", Imf::FloatAttribute( m_sampleSpacing ) );
    fileHeader.insert( "attenSampleSpacingType", Imf::IntAttribute( m_exponentialSampleSpacing ? 1 : 0 ) );

    if( isCubeMap ) {
        fileHeader.insert( "envmap", Imf::EnvmapAttribute( Imf::ENVMAP_CUBE ) );
    }

    fileHeader.channels().insert( "Z", Imf::Channel( Imf::FLOAT ) );

    for( int i = 0; i < m_sampleCount; ++i ) {
        std::string layerName = "sample" + frantic::strings::to_string( frantic::strings::zero_pad( i, 4 ) );
        fileHeader.channels().insert( ( layerName + ".AR" ).c_str(), Imf::Channel( Imf::FLOAT ) );
        fileHeader.channels().insert( ( layerName + ".AG" ).c_str(), Imf::Channel( Imf::FLOAT ) );
        fileHeader.channels().insert( ( layerName + ".AB" ).c_str(), Imf::Channel( Imf::FLOAT ) );
    }

    Imf::OutputFile fileOut( frantic::strings::to_string( depthMapPath ).c_str(), fileHeader );

    Imf::FrameBuffer fileFrame;
    fileFrame.insert( "Z", Imf::Slice( Imf::FLOAT, (char*)&m_depthMap->data()[imgSize.get_area() - imgSize.xsize],
                                       sizeof( float ), sizeof( float ) * -imgSize.xsize, 1, 1,
                                       std::numeric_limits<float>::max() ) );

    for( int i = 0; i < m_sampleCount; ++i ) {
        Imf::Slice layerSliceR( Imf::FLOAT, (char*)&m_attenBuffers[i].data()[imgSize.get_area() - imgSize.xsize].ar,
                                sizeof( frantic::graphics::alpha3f ),
                                sizeof( frantic::graphics::alpha3f ) * -imgSize.xsize );

        Imf::Slice layerSliceG( Imf::FLOAT, (char*)&m_attenBuffers[i].data()[imgSize.get_area() - imgSize.xsize].ag,
                                sizeof( frantic::graphics::alpha3f ),
                                sizeof( frantic::graphics::alpha3f ) * -imgSize.xsize );

        Imf::Slice layerSliceB( Imf::FLOAT, (char*)&m_attenBuffers[i].data()[imgSize.get_area() - imgSize.xsize].ab,
                                sizeof( frantic::graphics::alpha3f ),
                                sizeof( frantic::graphics::alpha3f ) * -imgSize.xsize );

        std::string layerName = get_exr_layer_name( i );
        fileFrame.insert( ( layerName + ".AR" ).c_str(), layerSliceR );
        fileFrame.insert( ( layerName + ".AG" ).c_str(), layerSliceG );
        fileFrame.insert( ( layerName + ".AB" ).c_str(), layerSliceB );
    }

    fileOut.setFrameBuffer( fileFrame );
    fileOut.writePixels( imgSize.ysize );
}

void base_atten_exr_saver::add_sample_internal( int x, int y, int cubefaceIndex, float z,
                                                const frantic::graphics::alpha3f& atten ) {

    int pixelX = x;
    int pixelY = y + m_depthMap->width() * cubefaceIndex;

    float imgZ = m_depthMap->get_depth_value( pixelX, pixelY );
    if( z < imgZ ) {
        // This should only happen if this pixel was previously un-initialized, otherwise we are writing out of order.
        if( imgZ != std::numeric_limits<float>::max() ) {
            // Clamp this sample to the first layer. This is kind of an error, since it means we were not writing in the
            // correct order!
            z = imgZ;
        } else {
            imgZ = z;
            m_depthMap->set_depth_value( pixelX, pixelY, z );
        }
    }

    int layer;
    if( m_exponentialSampleSpacing ) {
        layer = (int)floor( log( ( z - imgZ ) / m_sampleSpacing + 1 ) / M_LN2 );
    } else {
        layer = (int)floor( ( z - imgZ ) / m_sampleSpacing );
    }

    // Clamp to the final layer
    if( layer >= m_sampleCount )
        layer = m_sampleCount - 1;

    for( int i = layer; i < m_sampleCount; ++i )
        m_attenBuffers[i].blend_under( pixelX, pixelY, atten );
}

//
//
// exr implementation for creating/saving single face attenuation maps
//
//

singleface_atten_exr_saver::singleface_atten_exr_saver( const frantic::tstring& outputExrFilename,
                                                        frantic::graphics2d::size2 mapDim, int numSamples,
                                                        float spacing, bool exponentialGrowth ) {
    m_outputFilename = outputExrFilename;
    reset( mapDim, numSamples, spacing, exponentialGrowth );
    m_mapDim = mapDim;
}

void singleface_atten_exr_saver::add_sample( int x, int y, float z, const frantic::graphics::alpha3f& atten ) {
    add_sample_internal( x, y, 0, z, atten );
}

void singleface_atten_exr_saver::write_file() { save_exr_to_file( m_outputFilename, false ); }

//
//
// exr implementation for creating/saving cube face attenuation maps
//
//

cubeface_atten_exr_saver::cubeface_atten_exr_saver( const frantic::tstring& outputExrFilename, int mapWidth,
                                                    int numSamples, float spacing, bool exponentialGrowth ) {
    m_outputFilename = outputExrFilename;
    frantic::graphics2d::size2 mapDim( mapWidth, mapWidth * 6 );
    reset( mapDim, numSamples, spacing, exponentialGrowth );
    m_mapWidth = mapWidth;
}

void cubeface_atten_exr_saver::add_sample( int x, int y, int cubefaceIndex, float z,
                                           const frantic::graphics::alpha3f& atten ) {
    add_sample_internal( x, y, cubefaceIndex, z, atten );
}

void cubeface_atten_exr_saver::write_file() { save_exr_to_file( m_outputFilename, true ); }

} // namespace rendering
} // namespace frantic
