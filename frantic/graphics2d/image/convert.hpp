// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <vector>

#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/color_with_alpha.hpp>
#include <frantic/graphics2d/image/image_file_io.hpp>

namespace frantic {
namespace graphics2d {
namespace image {

using frantic::graphics::color4f;

// If the fusion header file Image.h is included, add this functionality
#ifdef _IMAGE_H_ // Digital Fusion image header
template <class Buffer>
void from_FusionImage( Buffer& buffer, Image* image, bool addAllChannels ) {
    buffer.set_size( size2( image->Width, image->Height ) );

    std::vector<std::string> channels;
    if( addAllChannels ) {
        channels.push_back( "BgR" );
        channels.push_back( "BgG" );
        channels.push_back( "BgB" );
        channels.push_back( "BgA" );
        channels.push_back( "Z" );
        channels.push_back( "U" );
        channels.push_back( "V" );
        channels.push_back( "Coverage" );
        channels.push_back( "ObjectID" );
        channels.push_back( "MaterialID" );
        channels.push_back( "NormalX" );
        channels.push_back( "NormalY" );
        channels.push_back( "NormalZ" );
        channels.push_back( "VectorX" );
        channels.push_back( "VectorY" );

        buffer.ensure_channels( channels );
        vector<float>& BgR = buffer.channel( "BgR" ).data();
        vector<float>& BgG = buffer.channel( "BgG" ).data();
        vector<float>& BgB = buffer.channel( "BgB" ).data();
        vector<float>& BgA = buffer.channel( "BgA" ).data();

        vector<float>& Z = buffer.channel( "Z" ).data();
        vector<float>& U = buffer.channel( "U" ).data();
        vector<float>& V = buffer.channel( "V" ).data();

        vector<float>& Coverage = buffer.channel( "Coverage" ).data();
        vector<float>& ObjectID = buffer.channel( "ObjectID" ).data();
        vector<float>& MaterialID = buffer.channel( "MaterialID" ).data();

        vector<float>& NormalX = buffer.channel( "NormalX" ).data();
        vector<float>& NormalY = buffer.channel( "NormalY" ).data();
        vector<float>& NormalZ = buffer.channel( "NormalZ" ).data();

        vector<float>& VectorX = buffer.channel( "VectorX" ).data();
        vector<float>& VectorY = buffer.channel( "VectorY" ).data();

        int i = 0;
        PixPtr ptr( image );
        for( int y = 0; y < buffer.height(); y++ ) {
            ptr.GotoXY( 0, y );
            for( int x = 0; x < buffer.width(); x++ ) {
                FltPixel p;
                p >>= ptr;
                color4f c( p.R, p.G, p.B, p.A );
                buffer.data()[i] = c;
                BgR[i] = p.BgR;
                BgG[i] = p.BgG;
                BgB[i] = p.BgB;
                BgA[i] = p.BgA;

                Z[i] = p.Z;
                U[i] = p.U;
                V[i] = p.V;

                Coverage[i] = p.Coverage;
                ObjectID[i] = (float)p.ObjectID;
                MaterialID[i] = (float)p.MaterialID;

                NormalX[i] = p.NX;
                NormalY[i] = p.NY;
                NormalZ[i] = p.NZ;

                VectorX[i] = p.VectX;
                VectorY[i] = p.VectY;
                i++;
            }
        }
    } else {

        int i = 0;
        PixPtr ptr( image );
        for( int y = 0; y < buffer.height(); y++ ) {
            ptr.GotoXY( 0, y );
            for( int x = 0; x < buffer.width(); x++ ) {
                FltPixel p;
                p >>= ptr;
                color4f c( p.R, p.G, p.B, p.A );
                buffer.data()[i] = c;
                i++;
            }
        }
    }
}

template <class Buffer> // AKA Image, just avoiding name collision
void to_FusionImage( Image* img, const Buffer& buffer ) {
    bool zExists = false;
    bool zCoverageExists = false;

    if( buffer.xsize() != img->Width || buffer.ysize() != img->Height )
        throw std::runtime_error( "image::to_FusionImage: Tried to add images with different dimensions." );

    image_channel<float>& z = image_channel<float>();
    if( buffer.channel_exists( "Z" ) ) {
        zExists = true;
        z = buffer.channel( "Z" );
    }

    image_channel<float>& zCoverage = image_channel<float>();
    if( buffer.channel_exists( "Coverage" ) ) {
        zCoverageExists = true;
        zCoverage = buffer.channel( "Coverage" );
    }

    PixPtr ptr( img );
    for( int y = 0; y < buffer.size().ysize; y++ ) {
        ptr.GotoXY( 0, y );
        for( int x = 0; x < buffer.size().xsize; x++ ) {
            FltPixel p;
            color4f c = buffer.get_pixel( x, y );
            p.R = c.c.r;
            p.G = c.c.g;
            p.B = c.c.b;
            p.A = c.a.a;
            if( zExists )
                p.Z = z.get_pixel( x, y );
            if( zCoverageExists )
                p.Coverage = zCoverage.get_pixel( x, y );
            ptr >>= p;
        }
    }
}

#endif // Digital Fusion image header

#ifdef IMGCLASSID // 3ds Max bitmap header

// TODO: define this for other colour types
inline void from_3dsMaxBitmap( Bitmap* bm, std::vector<color4f>& data, size2& size ) {
    if( size.xsize != bm->Width() || size.ysize != bm->Height() )
        throw std::runtime_error(
            "image::from_3dsMaxBitmap: Tried to retrieve a 3ds Max bitmap with mismatched sizes." );
    // Allocate a scanline for processing
    std::vector<BMM_Color_fl> scanline( size.xsize );
    int dataOffset = 0;
    for( int y = 0; y < size.ysize; ++y ) {
        // Get one scanline of pixels
        if( !bm->GetPixels( 0, size.ysize - y - 1, size.xsize, &scanline[0] ) )
            throw std::runtime_error( "image_utils::from_3dsMaxBitmap: Failed to retrieve scanline " +
                                      boost::lexical_cast<std::string>( y ) + " from a 3ds Max Bitmap." );
        for( int x = 0; x < size.xsize; ++x ) {
            data[dataOffset].c.r = scanline[x].r;
            data[dataOffset].c.g = scanline[x].g;
            data[dataOffset].c.b = scanline[x].b;
            data[dataOffset].a = scanline[x].a;
            ++dataOffset;
        }
    }
}

template <class ColorType, class AlphaType>
void to_3dsMaxBitmap( Bitmap* bm, const std::vector<frantic::graphics::color_with_alpha<ColorType, AlphaType>>& data,
                      size2 size ) {
    if( size.xsize != bm->Width() || size.ysize != bm->Height() )
        throw std::runtime_error( "image::to_3dsMaxBitmap: Tried to set a 3ds Max bitmap with mismatched sizes." );

    BitmapStorage* bs = bm->Storage();
    if( bs->Type() == BMM_REALPIX_32 )
        throw std::runtime_error( "image::to_3dsMaxBitmap: Cannot write to RealPixel format." );

    std::vector<BMM_Color_fl> scanline( size.xsize );
    int dataOffset = 0;
    for( int y = 0; y < size.ysize; ++y ) {
        for( int x = 0; x < size.xsize; ++x ) {
            scanline[x].r = data[dataOffset].c.r;
            scanline[x].g = data[dataOffset].c.g;
            scanline[x].b = data[dataOffset].c.b;
            scanline[x].a = data[dataOffset].a.to_float();
            ++dataOffset;
        }
        // Set one scanline of pixels
        if( !bm->PutPixels( 0, size.ysize - y - 1, size.xsize, &scanline[0] ) )
            throw std::runtime_error( "image::from_3dsMaxBitmap: Failed to write scanline " +
                                      boost::lexical_cast<std::string>( y ) + " to a 3ds Max Bitmap." );
    }
}

inline void to_3dsMaxBitmap( Bitmap* bm, const std::vector<float>& data, size2 size ) {
    if( size.xsize != bm->Width() || size.ysize != bm->Height() )
        throw std::runtime_error( "image::to_3dsMaxBitmap: Tried to set a 3ds Max bitmap with mismatched sizes." );

    BitmapStorage* bs = bm->Storage();
    if( bs->Type() == BMM_REALPIX_32 )
        throw std::runtime_error( "image::to_3dsMaxBitmap: Cannot write to RealPixel format." );

    std::vector<BMM_Color_fl> scanline( size.xsize );
    int dataOffset = 0;
    for( int y = 0; y < size.ysize; ++y ) {
        for( int x = 0; x < size.xsize; ++x ) {
            scanline[x].r = data[dataOffset];
            scanline[x].g = data[dataOffset];
            scanline[x].b = data[dataOffset];
            scanline[x].a = 1.f;
            ++dataOffset;
        }
        // Set one scanline of pixels
        if( !bm->PutPixels( 0, size.ysize - y - 1, size.xsize, &scanline[0] ) )
            throw std::runtime_error( "image::from_3dsMaxBitmap: Failed to write scanline " +
                                      boost::lexical_cast<std::string>( y ) + " to a 3ds Max Bitmap." );
    }
}

#endif // 3ds Max bitmap header

} // namespace image

} // namespace graphics2d
} // namespace frantic
