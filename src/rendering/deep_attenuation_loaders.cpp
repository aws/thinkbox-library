// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/rendering/deep_attenuation_loaders.hpp>

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
// abstract sampler classes for single face attenuation maps
//
//

frantic::graphics::alpha3f singleface_atten_loader::get_sample_bilinear( float x, float y, float z ) const {
    int fX = (int)floorf( x - 0.5f );
    int fY = (int)floorf( y - 0.5f );
    float alphaX = x - fX - 0.5f;
    float alphaY = y - fY - 0.5f;

    alpha3f a( 0.0f );
    if( (unsigned)fY < (unsigned)m_mapDim.ysize ) {
        if( (unsigned)fX < (unsigned)m_mapDim.xsize )
            a += get_sample( fX, fY, z ) * ( 1.f - alphaX ) * ( 1.f - alphaY );
        if( (unsigned)( fX + 1 ) < (unsigned)m_mapDim.xsize )
            a += get_sample( fX + 1, fY, z ) * alphaX * ( 1.f - alphaY );
    }
    if( (unsigned)( fY + 1 ) < (unsigned)m_mapDim.ysize ) {
        if( (unsigned)fX < (unsigned)m_mapDim.xsize )
            a += get_sample( fX, fY + 1, z ) * ( 1.f - alphaX ) * alphaY;
        if( (unsigned)( fX + 1 ) < (unsigned)m_mapDim.xsize )
            a += get_sample( fX + 1, fY + 1, z ) * alphaX * alphaY;
    }

    return a;
}

float singleface_atten_loader::get_zdepth_bilinear( float x, float y ) const {
    int fX = (int)floorf( x - 0.5f );
    int fY = (int)floorf( y - 0.5f );

    // gets the minimum z distance from all 4 neighbouring samples around this floating point sample.
    float minDepth = std::numeric_limits<float>::max();
    if( (unsigned)fY < (unsigned)m_mapDim.ysize ) {
        if( (unsigned)fX < (unsigned)m_mapDim.xsize )
            minDepth = std::min( minDepth, get_zdepth( fX, fY ) );
        if( (unsigned)( fX + 1 ) < (unsigned)m_mapDim.xsize )
            minDepth = std::min( minDepth, get_zdepth( fX + 1, fY ) );
    }
    if( (unsigned)( fY + 1 ) < (unsigned)m_mapDim.ysize ) {
        if( (unsigned)fX < (unsigned)m_mapDim.xsize )
            minDepth = std::min( minDepth, get_zdepth( fX, fY + 1 ) );
        if( (unsigned)( fX + 1 ) < (unsigned)m_mapDim.xsize )
            minDepth = std::min( minDepth, get_zdepth( fX + 1, fY + 1 ) );
    }
    return minDepth;
}

//
//
// abstract sampler classes for cube face attenuation maps
//
//

frantic::graphics::alpha3f cubeface_atten_loader::get_sample_bilinear( const frantic::graphics::vector3f& pos ) const {
    int imgSize = m_mapWidth;

    cube_face::default_cube_face cubeFaces[6];
    int numFaces = get_cube_faces( pos, (float)( imgSize - 2 ) / (float)imgSize, cubeFaces );

    frantic::graphics2d::vector2f coord;
    alpha3f result( 0.f );
    float totalWeight = 0.f;

    for( int i = 0; i < numFaces; ++i ) {
        const transform4f cubeXform = transform4f::from_cubeface( cubeFaces[i] );
        const float distCoeff = cubeFaces[i] >= 4 ? -1.f : 1.f;
        const float dist = distCoeff * ( cubeXform * pos ).z;

        coord = get_cube_face_coordinate( pos, cubeFaces[i] );
        coord += frantic::graphics2d::vector2f( 1, 1 );
        coord *= 0.5f * imgSize;

        float floorX = floor( coord.x );
        float alphaX = coord.x - floorX;
        int x = (int)floorX;

        float floorY = floor( coord.y );
        float alphaY = coord.y - floorY;
        int y = (int)floorY;

        if( (unsigned)y < (unsigned)imgSize ) {
            if( (unsigned)x < (unsigned)imgSize ) {
                float weight = ( 1.f - alphaX ) * ( 1.f - alphaY );
                result += get_sample( x, y, cubeFaces[i], dist ) * weight;
                totalWeight += weight;
            }
            if( (unsigned)( x + 1 ) < (unsigned)imgSize ) {
                float weight = alphaX * ( 1.f - alphaY );
                result += get_sample( x + 1, y, cubeFaces[i], dist ) * weight;
                totalWeight += weight;
            }
        }
        if( (unsigned)( y + 1 ) < (unsigned)imgSize ) {
            if( (unsigned)x < (unsigned)imgSize ) {
                float weight = ( 1.f - alphaX ) * alphaY;
                result += get_sample( x, y + 1, cubeFaces[i], dist ) * weight;
                totalWeight += weight;
            }
            if( (unsigned)( x + 1 ) < (unsigned)imgSize ) {
                float weight = alphaX * alphaY;
                result += get_sample( x + 1, y + 1, cubeFaces[i], dist ) * weight;
                totalWeight += weight;
            }
        }
    }

    return ( result / totalWeight );
}

float cubeface_atten_loader::get_zdepth_bilinear( const frantic::graphics::vector3f& pos ) const {
    int imgSize = m_mapWidth;

    cube_face::default_cube_face cubeFaces[6];
    int numFaces = get_cube_faces( pos, (float)( imgSize - 2 ) / (float)imgSize, cubeFaces );

    // gets the minimum z distance from all 4 neighbouring samples around this floating point sample, on ALL four faces.
    // do we need to make a non-filtered version of this function for anything?
    // using the same idea as get_sample_bilinear.
    frantic::graphics2d::vector2f coord;
    // float dist = pos.get_magnitude();

    float result = std::numeric_limits<float>::max();
    for( int i = 0; i < numFaces; ++i ) {
        coord = get_cube_face_coordinate( pos, cubeFaces[i] );
        coord += frantic::graphics2d::vector2f( 1, 1 );
        coord *= 0.5f * imgSize;

        int x = (int)floorf( coord.x );
        int y = (int)floorf( coord.y );

        if( (unsigned)y < (unsigned)imgSize ) {
            if( (unsigned)x < (unsigned)imgSize )
                result = std::min( result, get_zdepth( x, y, cubeFaces[i] ) );
            if( (unsigned)( x + 1 ) < (unsigned)imgSize )
                result = std::min( result, get_zdepth( x + 1, y, cubeFaces[i] ) );
        }
        if( (unsigned)( y + 1 ) < (unsigned)imgSize ) {
            if( (unsigned)x < (unsigned)imgSize )
                result = std::min( result, get_zdepth( x, y + 1, cubeFaces[i] ) );
            if( (unsigned)( x + 1 ) < (unsigned)imgSize )
                result = std::min( result, get_zdepth( x + 1, y + 1, cubeFaces[i] ) );
        }
    }

    return result;
}

//
//
// common data members and function class for exr loaders
//
//

std::string base_atten_exr_loader::get_exr_layer_name( int i ) {
    return "sample" + frantic::strings::to_string( frantic::strings::zero_pad( i, 4 ) );
}

void base_atten_exr_loader::load_exr_from_file( const frantic::tstring& depthMapPath, bool isCubeMap ) {

    m_depthMap.reset( new frantic::rendering::depthbuffer_singleface );

    Imf::InputFile inFile( frantic::strings::to_string( depthMapPath ).c_str() );

    m_sampleCount = inFile.header().typedAttribute<Imf::IntAttribute>( "attenSampleCount" ).value();
    m_sampleSpacing = inFile.header().typedAttribute<Imf::FloatAttribute>( "attenSampleSpacing" ).value();
    m_exponentialSampleSpacing = inFile.header().typedAttribute<Imf::IntAttribute>( "attenSampleSpacingType" ).value()
                                     ? true
                                     : false; // avoid MSVC bitching about conversion from int to bool

    Imath::Box2i dw = inFile.header().dataWindow();
    frantic::graphics2d::size2 imgSize( dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1 );

    if( isCubeMap ) {
        if( Imf::envmap( inFile.header() ) != Imf::ENVMAP_CUBE ) {
            std::stringstream ss;
            ss << "deep_attenuation_exr_base::load_from_file() The file \""
               << frantic::strings::to_string( depthMapPath ) << "\" was not an OpenEXR cubemap";
            throw std::runtime_error( ss.str() );
        }
        if( imgSize.ysize != 6 * imgSize.xsize ) {
            std::stringstream ss;
            ss << "deep_attenuation_exr_base::load_from_file() The file \""
               << frantic::strings::to_string( depthMapPath )
               << "\" had invalid dimensions for loading a cubemap: " << imgSize;
            throw std::runtime_error( ss.str() );
        }
    }

    m_depthMap->set_size( imgSize );
    m_depthMap->clear();

    m_attenBuffers.reset( new frantic::graphics2d::framebuffer<alpha3f>[m_sampleCount] );
    for( int i = 0; i < m_sampleCount; ++i ) {
        m_attenBuffers[i].set_size( imgSize );
        m_attenBuffers[i].fill( alpha3f( 0 ) );
    }

    Imf::FrameBuffer inBuffer;

    Imf::Slice sliceZ( Imf::FLOAT, (char*)&m_depthMap->data()[imgSize.get_area() - imgSize.xsize], sizeof( float ),
                       sizeof( float ) * -imgSize.xsize, 1, 1, std::numeric_limits<float>::max() );

    inBuffer.insert( "Z", sliceZ );

    for( int i = 0; i < m_sampleCount; ++i ) {
        Imf::Slice layerSliceR( Imf::FLOAT, (char*)&m_attenBuffers[i].data()[imgSize.get_area() - imgSize.xsize].ar,
                                sizeof( alpha3f ), sizeof( alpha3f ) * -imgSize.xsize );

        Imf::Slice layerSliceG( Imf::FLOAT, (char*)&m_attenBuffers[i].data()[imgSize.get_area() - imgSize.xsize].ag,
                                sizeof( alpha3f ), sizeof( alpha3f ) * -imgSize.xsize );

        Imf::Slice layerSliceB( Imf::FLOAT, (char*)&m_attenBuffers[i].data()[imgSize.get_area() - imgSize.xsize].ab,
                                sizeof( alpha3f ), sizeof( alpha3f ) * -imgSize.xsize );

        std::string layerName = get_exr_layer_name( i );
        inBuffer.insert( ( layerName + ".AR" ).c_str(), layerSliceR );
        inBuffer.insert( ( layerName + ".AG" ).c_str(), layerSliceG );
        inBuffer.insert( ( layerName + ".AB" ).c_str(), layerSliceB );
    }

    inFile.setFrameBuffer( inBuffer );
    inFile.readPixels( dw.min.y, dw.max.y );
}

frantic::graphics::alpha3f base_atten_exr_loader::get_sample_internal( int x, int y, int cubefaceIndex, float z,
                                                                       float& outImgZ ) const {
    int pixelX = x;
    int pixelY = y + m_depthMap->width() * cubefaceIndex;
    alpha3f outA( 0.0f );
    float imgZ = m_depthMap->get_depth_value_fast( pixelX, pixelY );
    outImgZ = imgZ;
    if( z > imgZ ) {
        if( m_exponentialSampleSpacing ) {
            // A depth layer starts at z = m_sampleSpacing * (2^n - 1), so a sample depth z falls into layer floor(
            // log(z / m_sampleSpacing + 1) / log(2) ).
            int layer = static_cast<int>( log( ( z - imgZ ) / m_sampleSpacing + 1 ) / M_LN2 );
            if( layer < 0 ) {
                layer = 0;
            }
            if( layer >= m_sampleCount )
                outA += m_attenBuffers[m_sampleCount - 1].get_pixel_fast( pixelX, pixelY );
            else {
                // The linear interpolation alpha value for depth 'z' in layer 'n' is ( z - m_sampleSpacing * (2^n - 1)
                // ) / ( m_sampleSpacing * 2^n )
                float layerWidth = m_sampleSpacing * float( 1 << layer ); //(1 << layer) == 2 ^ layer
                float alpha = ( z - imgZ - layerWidth + m_sampleSpacing ) / layerWidth;
                if( layer == 0 )
                    outA += alpha * m_attenBuffers[0].get_pixel_fast( pixelX, pixelY );
                else {
                    outA += ( 1.f - alpha ) * m_attenBuffers[layer - 1].get_pixel_fast( pixelX, pixelY );
                    outA += alpha * m_attenBuffers[layer].get_pixel_fast( pixelX, pixelY );
                }
            }
        } else {
            float layerDepth = ( z - imgZ ) / m_sampleSpacing;
            int layer = (int)layerDepth;

            if( layer < 0 ) {
                layer = 0;
            }

            if( layer >= m_sampleCount )
                outA += m_attenBuffers[m_sampleCount - 1].get_pixel_fast( pixelX, pixelY );
            else if( layer == 0 )
                outA += layerDepth * m_attenBuffers[0].get_pixel_fast( pixelX, pixelY );
            else {
                float alpha = ( layerDepth - layer );
                outA += ( 1.f - alpha ) * m_attenBuffers[layer - 1].get_pixel_fast( pixelX, pixelY );
                outA += alpha * m_attenBuffers[layer].get_pixel_fast( pixelX, pixelY );
            }
        }
    }
    return outA;
}

//
//
// exr implementation for loading/sampling single face attenuation maps
//
//

singleface_atten_exr_loader::singleface_atten_exr_loader( const frantic::tstring& filename ) {
    load_exr_from_file( filename, false );
    m_mapDim = m_depthMap->size();
}

frantic::graphics::alpha3f singleface_atten_exr_loader::get_sample( int x, int y, float z ) const {
    float garbage = 0;
    return get_sample_internal( x, y, 0, z, garbage );
}

frantic::graphics::alpha3f singleface_atten_exr_loader::get_sample( int x, int y, float z, float& outImgZ ) const {
    return get_sample_internal( x, y, 0, z, outImgZ );
}

float singleface_atten_exr_loader::get_zdepth( int x, int y ) const { return m_depthMap->get_depth_value_fast( x, y ); }

//
//
// exr implementation for loading/sampling cube face attenuation maps
//
//

cubeface_atten_exr_loader::cubeface_atten_exr_loader( const frantic::tstring& filename ) {
    load_exr_from_file( filename, true );
    m_mapWidth = m_depthMap->width();
}

frantic::graphics::alpha3f cubeface_atten_exr_loader::get_sample( int x, int y, int cubefaceIndex, float z ) const {
    float garbage = 0.f;
    return get_sample_internal( x, y, cubefaceIndex, z, garbage );
}

frantic::graphics::alpha3f cubeface_atten_exr_loader::get_sample( int x, int y, int cubefaceIndex, float z,
                                                                  float& outImgZ ) const {
    return get_sample_internal( x, y, cubefaceIndex, z, outImgZ );
}

float cubeface_atten_exr_loader::get_zdepth( int x, int y, int cubefaceIndex ) const {
    return m_depthMap->get_depth_value_fast( x, y + m_mapWidth * cubefaceIndex );
}

//
//
// convenience factory function to create a exr singleface/cubeface loader
//
//

boost::shared_ptr<singleface_atten_loader> create_singleface_atten_loader( const frantic::tstring& filename ) {
    boost::shared_ptr<singleface_atten_loader> attenLoader;
    frantic::tstring ext = frantic::strings::to_lower( frantic::files::extension_from_path( filename ) );
    // Only layered exrs currently supported, but could be extended to load different file types. For example, it used
    // to load DTEX files.
    //  .cxr files are corona's multichannel exr file format
    if( ext == _T(".exr") || ext == _T(".cxr") ) {
        attenLoader.reset( new frantic::rendering::singleface_atten_exr_loader( filename ) );
    } else {
        throw std::runtime_error( "Attenuation saving error: " + frantic::strings::to_string( ext ) +
                                  " is not a known extension type. Please use a .exr file." );
    }
    return attenLoader;
}

boost::shared_ptr<cubeface_atten_loader> create_cubeface_atten_loader( const frantic::tstring& filename ) {
    boost::shared_ptr<cubeface_atten_loader> attenLoader;
    frantic::tstring ext = frantic::strings::to_lower( frantic::files::extension_from_path( filename ) );
    // Only layered exrs currently supported, but could be extended to load different file types. For example, it used
    // to load a vector of 6 DTEX files.
    //  .cxr files are corona's multichannel exr file format
    if( ext == _T(".exr") || ext == _T(".cxr") )
        attenLoader.reset( new cubeface_atten_exr_loader( filename ) );

    if( !attenLoader )
        throw std::runtime_error( "Attenuation saving error: Please specify a cube-faced .exr file." );
    return attenLoader;
}

} // namespace rendering
} // namespace frantic
