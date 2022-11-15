// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

// TODO: Need to delete this file, and fully use the new image I/O stuff.

#include <ImfArray.h>
#include <ImfBoxAttribute.h>
#include <ImfChannelList.h>
#include <ImfIO.h>
#include <ImfInputFile.h>
#include <ImfLineOrder.h>
#include <ImfOutputFile.h>
#include <ImfRgbaFile.h>
#include <half.h>

#include <iostream>
#include <vector>

#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/color3h.hpp>
#include <frantic/graphics/color_with_alpha.hpp>
#include <frantic/graphics2d/size2.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::graphics;
using namespace Imf;
using namespace Imath;

namespace frantic {
namespace graphics2d {
namespace image {

void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::alpha3f>& data,
                       const frantic::graphics2d::size2& size ) {
    using namespace Imf;
    using namespace Imath;

    RgbaOutputFile file( path.c_str(), size.xsize, size.ysize, WRITE_RGB );

    std::vector<Rgba> pixels( size.get_area() );

    // I tried using line order of DECREASING_Y but that didn't change the output for some
    // reason, so I'm just going to reverse the data myself. -batty 2005/05/06
    // Same for reading

    for( int row = 0; row < size.ysize; row++ ) {
        for( int col = 0; col < size.xsize; col++ ) {

            int inindex = col + row * size.xsize;

            alpha3f source = data[inindex];

            Rgba dest;
            dest.r = source.ar;
            dest.g = source.ag;
            dest.b = source.ab;
            dest.a = 1;

            int outindex = col + ( size.ysize - row - 1 ) * size.xsize;
            pixels[outindex] = dest;
        }
    }

    file.setFrameBuffer( &pixels[0], 1, size.xsize );
    file.writePixels( size.ysize );
}

void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::color4f>& data,
                       const frantic::graphics2d::size2& size ) {
    using namespace Imf;
    using namespace Imath;

    RgbaOutputFile file( path.c_str(), size.xsize, size.ysize, WRITE_RGBA );

    std::vector<Rgba> pixels( size.get_area() );

    // I tried using line order of DECREASING_Y but that didn't change the output for some
    // reason, so I'm just going to reverse the data myself. -batty 2005/05/06
    // Same for reading

    for( int row = 0; row < size.ysize; row++ ) {
        for( int col = 0; col < size.xsize; col++ ) {

            int inindex = col + row * size.xsize;
            color4f source = data[inindex].to_color4f();

            Rgba dest;
            dest.r = source.color().r;
            dest.g = source.color().g;
            dest.b = source.color().b;
            dest.a = source.alpha().a;

            int outindex = col + ( size.ysize - row - 1 ) * size.xsize;
            pixels[outindex] = dest;
        }
    }

    file.setFrameBuffer( &pixels[0], 1, size.xsize );
    file.writePixels( size.ysize );
}

void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::color6f>& data,
                       const frantic::graphics2d::size2& size ) {
    RgbaOutputFile file( path.c_str(), size.xsize, size.ysize, WRITE_RGBA );

    std::vector<Rgba> pixels( size.get_area() );

    // I tried using line order of DECREASING_Y but that didn't change the output for some
    // reason, so I'm just going to reverse the data myself. -batty 2005/05/06
    // Same for reading

    for( int row = 0; row < size.ysize; row++ ) {
        for( int col = 0; col < size.xsize; col++ ) {

            int inindex = col + row * size.xsize;
            frantic::graphics::color4f source = data[inindex].to_color4f();

            Rgba dest;
            dest.r = source.color().r;
            dest.g = source.color().g;
            dest.b = source.color().b;
            dest.a = source.alpha().a;

            int outindex = col + ( size.ysize - row - 1 ) * size.xsize;
            pixels[outindex] = dest;
        }
    }

    file.setFrameBuffer( &pixels[0], 1, size.xsize );
    file.writePixels( size.ysize );
}

void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::color3f>& data,
                       const frantic::graphics2d::size2& size ) {
    RgbaOutputFile file( path.c_str(), size.xsize, size.ysize, WRITE_RGB );

    std::vector<Rgba> pixels( size.get_area() );

    // I tried using line order of DECREASING_Y but that didn't change the output for some
    // reason, so I'm just going to reverse the data myself. -batty 2005/05/06
    // Same for reading

    for( int row = 0; row < size.ysize; row++ ) {
        for( int col = 0; col < size.xsize; col++ ) {

            int inindex = col + row * size.xsize;

            frantic::graphics::color3f source = data[inindex];

            Rgba dest;
            dest.r = source.r;
            dest.g = source.g;
            dest.b = source.b;
            dest.a = 1;

            int outindex = col + ( size.ysize - row - 1 ) * size.xsize;
            pixels[outindex] = dest;
        }
    }

    file.setFrameBuffer( &pixels[0], 1, size.xsize );
    file.writePixels( size.ysize );
}

void write_to_OpenEXR( const std::string& path, const std::vector<float>& data,
                       const frantic::graphics2d::size2& size ) {
    Header header( size.xsize, size.ysize );
    header.channels().insert( "R", Channel( Imf::FLOAT ) );
    header.channels().insert( "G", Channel( Imf::FLOAT ) );
    header.channels().insert( "B", Channel( Imf::FLOAT ) );

    OutputFile file( path.c_str(), header );

    FrameBuffer frameBuffer;

    frameBuffer.insert( "R", Slice( Imf::FLOAT, (char*)( &data[( size.ysize - 1 ) * size.xsize] ),
                                    sizeof( data[0] ) * 1,           // xStride
                                    sizeof( data[0] ) * -size.xsize, // yStride
                                    1, 1, 0.0 ) );

    frameBuffer.insert( "G", Slice( Imf::FLOAT, (char*)( &data[( size.ysize - 1 ) * size.xsize] ),
                                    sizeof( data[0] ) * 1,           // xStride
                                    sizeof( data[0] ) * -size.xsize, // yStride
                                    1, 1, 0.0 ) );

    frameBuffer.insert( "B", Slice( Imf::FLOAT, (char*)( &data[( size.ysize - 1 ) * size.xsize] ),
                                    sizeof( data[0] ) * 1,           // xStride
                                    sizeof( data[0] ) * -size.xsize, // yStride
                                    1, 1, 0.0 ) );

    file.setFrameBuffer( frameBuffer );
    file.writePixels( size.ysize );
}

void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::alpha3f>& data,
                        frantic::graphics2d::size2& size ) {
    InputFile file( path.c_str() );
    Box2i dw = file.header().dataWindow();

    size.xsize = dw.max.x - dw.min.x + 1;
    size.ysize = dw.max.y - dw.min.y + 1;
    data.resize( size.get_area() );

    FrameBuffer frameBuffer;
    frameBuffer.insert( "R",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].ar - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( alpha3f ) * 1,           // xStride
                               sizeof( alpha3f ) * -size.xsize, // yStride
                               1, 1,                            // x/y sampling
                               0.0 ) );                         // fillValue

    frameBuffer.insert( "G",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].ag - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( alpha3f ) * 1,           // xStride
                               sizeof( alpha3f ) * -size.xsize, // yStride
                               1, 1,                            // x/y sampling
                               0.0 ) );                         // fillValue

    frameBuffer.insert( "B",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].ab - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( alpha3f ) * 1,           // xStride
                               sizeof( alpha3f ) * -size.xsize, // yStride
                               1, 1,                            // x/y sampling
                               0.0 ) );                         // fillValue

    file.setFrameBuffer( frameBuffer );
    file.readPixels( dw.min.y, dw.max.y );
}

void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::color3f>& data,
                        frantic::graphics2d::size2& size ) {
    InputFile file( path.c_str() );
    Box2i dw = file.header().dataWindow();

    size.xsize = dw.max.x - dw.min.x + 1;
    size.ysize = dw.max.y - dw.min.y + 1;
    data.resize( size.get_area() );

    FrameBuffer frameBuffer;
    frameBuffer.insert( "R",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].r - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( color3f ) * 1,           // xStride
                               sizeof( color3f ) * -size.xsize, // yStride
                               1, 1,                            // x/y sampling
                               0.0 ) );                         // fillValue

    frameBuffer.insert( "G",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].g - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( color3f ) * 1,           // xStride
                               sizeof( color3f ) * -size.xsize, // yStride
                               1, 1,                            // x/y sampling
                               0.0 ) );                         // fillValue

    frameBuffer.insert( "B",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].b - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( color3f ) * 1,           // xStride
                               sizeof( color3f ) * -size.xsize, // yStride
                               1, 1,                            // x/y sampling
                               0.0 ) );                         // fillValue

    file.setFrameBuffer( frameBuffer );
    file.readPixels( dw.min.y, dw.max.y );
}

void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::color3h>& outData,
                        frantic::graphics2d::size2& outSize ) {
    FrameBuffer frameBuffer;

    InputFile file( path.c_str() );

    Box2i dw = file.header().dataWindow();

    outSize.xsize = dw.max.x - dw.min.x + 1;
    outSize.ysize = dw.max.y - dw.min.y + 1;
    outData.resize( outSize.get_area() );

    frameBuffer.insert( "R",              // name
                        Slice( Imf::HALF, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &outData[outSize.ysize * outSize.xsize - outSize.xsize].r - // base
                                        dw.min.x - dw.min.y * outSize.xsize ),
                               sizeof( color3h ),                  // xStride
                               sizeof( color3h ) * -outSize.xsize, // yStride
                               1, 1,                               // x/y sampling
                               0.0 ) );                            // fillValue

    frameBuffer.insert( "G",              // name
                        Slice( Imf::HALF, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &outData[outSize.ysize * outSize.xsize - outSize.xsize].g - // base
                                        dw.min.x - dw.min.y * outSize.xsize ),
                               sizeof( color3h ),                  // xStride
                               sizeof( color3h ) * -outSize.xsize, // yStride
                               1, 1,                               // x/y sampling
                               0.0 ) );                            // fillValue

    frameBuffer.insert( "B",              // name
                        Slice( Imf::HALF, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &outData[outSize.ysize * outSize.xsize - outSize.xsize].b - // base
                                        dw.min.x - dw.min.y * outSize.xsize ),
                               sizeof( color3h ),                  // xStride
                               sizeof( color3h ) * -outSize.xsize, // yStride
                               1, 1,                               // x/y sampling
                               0.0 ) );                            // fillValue

    file.setFrameBuffer( frameBuffer );
    file.readPixels( dw.min.y, dw.max.y );
}

void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::color4f>& data,
                        frantic::graphics2d::size2& size ) {
    InputFile file( path.c_str() );
    Box2i dw = file.header().dataWindow();

    size.xsize = dw.max.x - dw.min.x + 1;
    size.ysize = dw.max.y - dw.min.y + 1;
    data.resize( size.get_area() );

    FrameBuffer frameBuffer;
    frameBuffer.insert( "R",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].c.r - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( frantic::graphics::color4f ) * 1,           // xStride
                               sizeof( frantic::graphics::color4f ) * -size.xsize, // yStride
                               1, 1,                                               // x/y sampling
                               0.0 ) );                                            // fillValue

    frameBuffer.insert( "G",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].c.g - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( frantic::graphics::color4f ) * 1,           // xStride
                               sizeof( frantic::graphics::color4f ) * -size.xsize, // yStride
                               1, 1,                                               // x/y sampling
                               0.0 ) );                                            // fillValue

    frameBuffer.insert( "B",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].c.b - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( frantic::graphics::color4f ) * 1,           // xStride
                               sizeof( frantic::graphics::color4f ) * -size.xsize, // yStride
                               1, 1,                                               // x/y sampling
                               0.0 ) );                                            // fillValue

    frameBuffer.insert( "A",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize].a.a - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( frantic::graphics::color4f ) * 1,           // xStride
                               sizeof( frantic::graphics::color4f ) * -size.xsize, // yStride
                               1, 1,                                               // x/y sampling
                               1.0 ) );                                            // fillValue

    file.setFrameBuffer( frameBuffer );
    file.readPixels( dw.min.y, dw.max.y );
}

void read_from_OpenEXR( const std::string& path, std::vector<float>& data, frantic::graphics2d::size2& size ) {
    InputFile file( path.c_str() );
    Box2i dw = file.header().dataWindow();

    size.xsize = dw.max.x - dw.min.x + 1;
    size.ysize = dw.max.y - dw.min.y + 1;
    data.resize( size.get_area() );

    FrameBuffer frameBuffer;
    frameBuffer.insert( "R",               // name
                        Slice( Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                               (char*)( &data[size.ysize * size.xsize - size.xsize] - // base
                                        dw.min.x - dw.min.y * size.xsize ),
                               sizeof( float ) * 1,           // xStride
                               sizeof( float ) * -size.xsize, // yStride
                               1, 1,                          // x/y sampling
                               0.0 ) );                       // fillValue
                                                              /*
                                                                frameBuffer.insert ("G", // name
                                                                  Slice (Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                                                                  (char *) (&data[size.ysize*size.xsize-size.xsize].g - // base
                                                                  dw.min.x - dw.min.y * size.xsize),
                                                                  sizeof (color3f) * 1, // xStride
                                                                  sizeof (color3f) * -size.xsize,// yStride
                                                                  1, 1, // x/y sampling
                                                                  0.0)); // fillValue
                                                          
                                                                frameBuffer.insert ("B", // name
                                                                  Slice (Imf::FLOAT, // type (must be fully qualified to disambiguate from windows.h)
                                                                  (char *) (&data[size.ysize*size.xsize-size.xsize].b - // base
                                                                  dw.min.x - dw.min.y * size.xsize),
                                                                  sizeof (color3f) * 1, // xStride
                                                                  sizeof (color3f) * -size.xsize,// yStride
                                                                  1, 1, // x/y sampling
                                                                  0.0)); // fillValue
                                                              */
    file.setFrameBuffer( frameBuffer );
    file.readPixels( dw.min.y, dw.max.y );
}

frantic::graphics2d::size2 get_exr_size( const std::string& filename ) {
    InputFile file( filename.c_str() );
    Box2i dw = file.header().dataWindow();

    return frantic::graphics2d::size2( dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1 );
}

frantic::graphics2d::size2 get_exr_display_size( const std::string& filename ) {
    InputFile file( filename.c_str() );
    Box2i dw = file.header().displayWindow();

    return frantic::graphics2d::size2( dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1 );
}

} // namespace image
} // namespace graphics2d
} // namespace frantic
