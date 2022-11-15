// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <algorithm>
#include <cmath>
#include <stack>
#include <vector>

#include <boost/shared_ptr.hpp>

// No OpenEXR for MSVC6
#if( defined( _MSC_VER ) && _MSC_VER > 1200 ) || __GNUC__
#include <ImfIO.h>
#endif

#if( defined( _WIN32 ) || defined( _WIN64 ) ) && !defined( _INC_WINDOWS )
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif

#if( defined( _WIN32 ) || defined( _WIN64 ) )
#include <frantic/win32/utility.hpp>
#include <frantic/win32/wrap.hpp>
#endif

#include <frantic/geometry/polygon_utils.hpp>
#include <frantic/graphics/color_with_alpha.hpp>
#include <frantic/graphics2d/file_io/dpx_file_io.hpp>
#include <frantic/graphics2d/image/convert.hpp>
#include <frantic/graphics2d/image/image_file_io.hpp>
#include <frantic/graphics2d/image_channel.hpp>
#include <frantic/graphics2d/rasterization.hpp>
#include <frantic/math/utils.hpp>

// TODO: This file is really long - we should refactor it somehow into multiple
//       files in a way which makes sense.

namespace frantic {
namespace graphics2d {

// using frantic::graphics::vector3f;

namespace draw_point_filter {
enum draw_point_filter_enum {
    nearest_neighbor,
    bilinear, // averages the pixel energy over the nearest 4 pixels
    bicubic   // averages the pixel energy using a 3x3 filter
};

inline draw_point_filter_enum from_string( const std::string& filter ) {
    if( filter == "Nearest Neighbor" ) {
        return nearest_neighbor;
    } else if( filter == "Bilinear" ) {
        return bilinear;
    } else if( filter == "Bicubic" ) {
        return bicubic;
    } else {
        throw std::runtime_error( "draw_point_filter::from_string() - Invalid draw point filter, \"" + filter + "\"." );
    }
}
} // namespace draw_point_filter

namespace pixel_lookup_filter {
enum pixel_lookup_filter_enum { nearest_neighbor_filter, bilinear_filter, bicubic_filter };
}

// Some extra functions to help with pixel blending
template <class ColorType>
inline void color_blend_under( ColorType& outDest, const ColorType& source ) {
    outDest.blend_under( source );
}
inline void color_blend_under( float& outDest, float source ) { outDest += source; }
template <class ColorType>
inline void color_blend_over( ColorType& outDest, const ColorType& source ) {
    outDest.blend_over( source );
}
inline void color_blend_over( float& outDest, float source ) { outDest += source; }

template <typename ColorType>
class framebuffer : public image_channel<ColorType> {
  public:
    // A member-function typedef for the get pixel function
    typedef void ( framebuffer<ColorType>::*set_pixel_function_type )( vector2f pixelLocation, const ColorType& color );
    typedef ColorType ( framebuffer<ColorType>::*get_pixel_function_type )( const vector2f& pixelLocation ) const;

  private:
    set_pixel_function_type m_drawPointFunction;
    get_pixel_function_type m_getPixelFunction;

  protected:
    std::vector<boost::shared_ptr<image_channel<float>>> m_channels;

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
  public:
    framebuffer() {
        m_drawPointFunction = &framebuffer<ColorType>::draw_point_bilinear;
        m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_bilinear;
    }
    framebuffer( int squareSize )
        : image_channel<ColorType>( size2( squareSize, squareSize ) ) {
        m_drawPointFunction = &framebuffer<ColorType>::draw_point_bilinear;
        m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_bilinear;
    }
    framebuffer( size2 size )
        : image_channel<ColorType>( size ) {
        m_drawPointFunction = &framebuffer<ColorType>::draw_point_bilinear;
        m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_bilinear;
    }
    framebuffer( std::vector<ColorType> data, size2 size_ )
        : image_channel<ColorType>( data, size_ ) {
        m_drawPointFunction = &framebuffer<ColorType>::draw_point_bilinear;
        m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_bilinear;
    }
    framebuffer( framebuffer<ColorType> data, size2 size_ )
        : image_channel<ColorType>( data.m_data, size_ ) {
        m_drawPointFunction = &framebuffer<ColorType>::draw_point_bilinear;
        m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_bilinear;
    }
    ~framebuffer() {}

    void clear() {
        image_channel<ColorType>::clear();
        for( unsigned i = 0; i < m_channels.size(); i++ )
            m_channels[i]->clear_helper();

        m_drawPointFunction = &framebuffer<ColorType>::draw_point_bilinear;
    }

    void set_size( size2 size ) {
        image_channel<ColorType>::set_size( size );
        for( unsigned i = 0; i < m_channels.size(); i++ )
            m_channels[i]->set_size_helper( size );
    }

    void swap( framebuffer<ColorType>& rhs ) {
        image_channel<ColorType>::swap( rhs );
        std::swap( m_channels, rhs.m_channels );
    }

    void set_draw_point_filter( draw_point_filter::draw_point_filter_enum mode ) {
        switch( mode ) {
        case draw_point_filter::nearest_neighbor:
            m_drawPointFunction = &framebuffer<ColorType>::draw_point_nearest_neighbor;
            break;
        case draw_point_filter::bilinear:
            m_drawPointFunction = &framebuffer<ColorType>::draw_point_bilinear;
            break;
        case draw_point_filter::bicubic:
            m_drawPointFunction = &framebuffer<ColorType>::draw_point_bicubic;
            break;
        default:
            throw std::runtime_error( "framebuffer.set_draw_point_filter: Invalid filter specified." );
        }
    }

    draw_point_filter::draw_point_filter_enum get_draw_point_filter() const {
        if( m_drawPointFunction == &framebuffer<ColorType>::draw_point_nearest_neighbor )
            return draw_point_filter::nearest_neighbor;
        else if( m_drawPointFunction == &framebuffer<ColorType>::draw_point_bilinear )
            return draw_point_filter::bilinear;
        else if( m_drawPointFunction == &framebuffer<ColorType>::draw_point_bicubic )
            return draw_point_filter::bicubic;
        else
            throw std::runtime_error(
                "framebuffer.get_draw_point_filter: The draw point filter was in an invalid state, "
                "this is likely a bug in the software somewhere." );
    }

    void set_pixel_lookup_filter( pixel_lookup_filter::pixel_lookup_filter_enum lookupFilter ) {
        switch( lookupFilter ) {
        case pixel_lookup_filter::nearest_neighbor_filter:
            m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_nearest_neighbor;
            break;
        case pixel_lookup_filter::bilinear_filter:
            m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_bilinear;
            break;
        case pixel_lookup_filter::bicubic_filter:
            m_getPixelFunction = &framebuffer<ColorType>::get_pixel_filter_bicubic;
            break;
        default:
            throw std::runtime_error( "framebuffer.set_pixel_lookup_filter: Invalid pixel lookup filter specified." );
        }
    }

    pixel_lookup_filter::pixel_lookup_filter_enum get_pixel_lookup_filter() const {
        if( m_getPixelFunction == &framebuffer<ColorType>::get_pixel_nearest_neighbor )
            return pixel_lookup_filter::nearest_neighbor_filter;
        else if( m_getPixelFunction == &framebuffer<ColorType>::get_pixel_bilinear )
            return pixel_lookup_filter::bilinear_filter;
        else if( m_getPixelFunction == &framebuffer<ColorType>::get_pixel_bicubic )
            return pixel_lookup_filter::bicubic_filter;
        else
            throw std::runtime_error(
                "framebuffer.get_pixel_lookup_filter: The pixel lookup filter was in an invalid state, "
                "this is likely a bug in the software somewhere." );
    }

    template <class OtherColorType>
    void copy_from( const framebuffer<OtherColorType>& source ) {
        image_channel<ColorType>::copy_from( source );

        // Delete all our channels
        clear_channels();
        for( int chan = 0; chan < source.channel_count(); ++chan ) {
            // Add channel using name and data from source
            add_channel( source.channel( chan ).name(), source.channel( chan ).data() );
        }
    }

    // TODO: This function should do some error checking about the sizes, image allocation, etc.
    void resize( const size2& newSize ) {
        framebuffer<ColorType> resized( newSize );
        resized.ensure_channels( channel_names() );

        image::resize_up<ColorType>( resized.data(), newSize, this->data(), this->size() );
        for( int chan = 0; chan < channel_count(); ++chan )
            image::resize_up<float>( resized.channel( channel( chan ).name() ).data(), newSize, channel( chan ).data(),
                                     this->size() );

        swap( resized );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Channel functions

    void clear_channels() {
        m_channels.clear();
        std::vector<boost::shared_ptr<image_channel<float>>> chan;
        m_channels.swap( chan );
    }

    void set_channels( const std::vector<std::string> names ) {
        // Remove all other channels
        clear_channels();

        for( unsigned int i = 0; i < names.size(); ++i )
            add_channel( names[i] );
    }

    void ensure_channels( const std::vector<std::string> namesList ) {
        // These are copies because we sort them
        std::vector<std::string> channelNames( channel_names() );
        std::vector<std::string> names( namesList );

        std::sort( names.begin(), names.end() );
        std::sort( channelNames.begin(), channelNames.end() );

        int channelNamesIndex = 0, namesIndex = 0;

        while( namesIndex < (int)names.size() ) {
            // If there are any channels left to look at, compare. Otherwise just add it
            int cmp = channelNamesIndex < (int)channelNames.size()
                          ? names[namesIndex].compare( channelNames[channelNamesIndex] )
                          : -1;
            if( cmp == 0 ) {
                namesIndex++;
                channelNamesIndex++;
            } else if( cmp < 0 ) {
                add_channel( names[namesIndex] );
                namesIndex++;
            } else {
                channelNamesIndex++;
            }
        }
    }

    std::vector<std::string> channel_names() {
        std::vector<std::string> names( channel_count() );
        int channelCount = channel_count();
        for( int i = 0; i < channelCount; ++i )
            names[i] = channel( i ).name();

        return names;
    }

    int add_channel( const std::string& name ) {
        int id = channel_id_nothrow( name );
        if( id >= 0 )
            return id;

        m_channels.push_back(
            boost::shared_ptr<image_channel<float>>( new image_channel<float>( (image_channel_base*)this, name ) ) );
        return (int)m_channels.size() - 1;
    }

    int add_channel( const std::string& name, const std::vector<float>& data ) {
        if( this->size().get_area() != (int)data.size() )
            throw std::runtime_error( "framebuffer.add_channel: Tried to add a channel of area " +
                                      boost::lexical_cast<std::string>( data.size() ) + " to an image of area " +
                                      boost::lexical_cast<std::string>( this->size().get_area() ) );

        int id = channel_id_nothrow( name );
        if( id >= 0 ) {
            image_channel<float>& chan = channel( id );
            chan.data().insert( chan.data().begin(), data.begin(), data.end() );
            return id;
        }

        m_channels.push_back( boost::shared_ptr<image_channel<float>>(
            new image_channel<float>( (image_channel_base*)this, name, data ) ) );
        return (int)m_channels.size() - 1;
    }

    bool channel_exists( const std::string& name ) const {
        for( int i = 0; i < channel_count(); ++i )
            if( channel( i ).name().compare( name ) == 0 )
                return true;
        return false;
    }

    int channel_id( const std::string& name ) const {
        for( int i = 0; i < channel_count(); ++i )
            if( channel( i ).name().compare( name ) == 0 )
                return i;

        throw std::runtime_error( "framebuffer.channel_id: Channel of name \"" + name + "\" does not exist." );
        // return -1;
    }

    int channel_id_nothrow( const std::string& name ) const {
        for( int i = 0; i < channel_count(); ++i )
            if( channel( i ).name().compare( name ) == 0 )
                return i;
        return -1;
    }

    int channel_count() const { return (int)m_channels.size(); }

    const image_channel<float>& channel( int id ) const { return *( m_channels[id] ); }

    image_channel<float>& channel( int id ) { return *( m_channels[id] ); }

    const image_channel<float>& channel( const std::string& name ) const { return channel( channel_id( name ) ); }

    image_channel<float>& channel( const std::string& name ) { return channel( channel_id( name ) ); }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Pixel Get and Set for channels

    // Three wrappers for the image_channel get_pixel functions. These wrappers all have the same signature
    // so they can be stored in the get_pixel_function function pointer.

    // point normalized to [-1,1]
    ColorType get_pixel_filter_nearest_neighbor( const vector2f& normalizedCoord ) const {
        return this->get_pixel_nearest_neighbor( normalizedCoord );
    }

    // point normalized to [-1,1]
    ColorType get_pixel_filter_bilinear( const vector2f& normalizedCoord ) const {
        return image::get_pixel_bilinear_duplicate_duplicate( normalizedCoord, this->m_data, this->m_size );
    }

    // point normalized to [-1,1]
    ColorType get_pixel_filter_bicubic( const vector2f& normalizedCoord ) const {
        return image::get_pixel_bicubic_duplicate_duplicate( normalizedCoord, this->m_data, this->m_size );
    }

    ColorType get_pixel_filtered( const vector2f& lookupPoint ) const {
        vector2f normalizedCoord( ( lookupPoint.x / this->size().xsize - 0.5f ) * 2,
                                  ( lookupPoint.y / this->size().ysize - 0.5f ) * 2 );
        return ( this->*m_getPixelFunction )( normalizedCoord );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Pixel func's specific to color data (not float)

    void blend_under( vector2 pixel, const ColorType& color ) {
        if( this->size().is_coord_valid( pixel ) ) {
            int index = this->size().get_index( pixel );
            color_blend_under( this->data()[index], color );
        }
    }

    void blend_over( vector2 pixel, const ColorType& color ) {
        if( this->size().is_coord_valid( pixel ) ) {
            int index = this->size().get_index( pixel );
            color_blend_over( this->data()[index], color );
        } else
            throw std::runtime_error( "frantic::graphics2d::framebuffer.blend_over: Requested pixel index, " +
                                      pixel.str() + ", out of bounds of image size, " + this->size().str() );
    }

    void blend_under( int ix, int iy, const ColorType& color ) { blend_under( vector2( ix, iy ), color ); }

    void blend_over( int ix, int iy, const ColorType& color ) { blend_over( vector2( ix, iy ), color ); }

    void blend_under( const framebuffer<ColorType>& bitmapOnTop ) {
        if( this->size() != bitmapOnTop.size() )
            throw std::runtime_error( "framebuffer.blend_under: Tried to blend bitmaps of incompatible sizes " +
                                      this->size().str() + " and " + bitmapOnTop.size().str() );

        // TODO: Deal with the additional arbitrary channels
        for( unsigned i = 0; i < this->data().size(); ++i ) {
            color_blend_under( this->data()[i], bitmapOnTop.m_data[i] );
        }
    }

    void blend_over( const framebuffer<ColorType>& bitmapOnTop ) {
        if( this->size() != bitmapOnTop.size() )
            throw std::runtime_error( "framebuffer.blend_over: Tried to blend bitmaps of incompatible sizes " +
                                      this->size().str() + " and " + bitmapOnTop.size().str() );

        // TODO: Deal with the additional arbitrary channels
        for( unsigned i = 0; i < this->data().size(); ++i ) {
            color_blend_over( this->data()[i], bitmapOnTop.m_data[i] );
        }
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Channel Functions

    void apply_gamma( float gammaCorrect ) {
        for( unsigned int index = 0; index < this->data().size(); ++index ) {
            this->data()[index].apply_gamma( gammaCorrect );
        }
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------------------------------

    // The area of a voxel on screen at a distance d from the camera is (w/(2*d*tan(fov/2))^2.  This scaling factor is
    // everything except for the d. This value is calculated at the center of the cube faces.
    /*float draw_point_scaling_constant( float xFov ) const {
      // TODO: This needs to be changed for orthographic and spherical views
      float factor = m_size.width / (2 * tan(xFov * 0.5f));
      return factor * factor;
    }*/

    void draw_point( vector2f pixelLocation, const ColorType& color ) {
        return ( this->*m_drawPointFunction )( pixelLocation, color );
    }

    void draw_point_nearest_neighbor( vector2f pixelLocation, const ColorType& color ) {
        int ix = (int)floorf( pixelLocation.x );
        int iy = (int)floorf( pixelLocation.y );
        blend_under( ix, iy, color );
    }

    // If we consider each pixel as being a small cube, the pixel 0,0 is at location [0.5,0.5],
    // and the pixel width-1,height-1 is at [width-0.5,height-0.5]
    void draw_point_bilinear( vector2f pixelLocation, const ColorType& color ) {
        int ix = (int)floorf( pixelLocation.x - 0.5f );
        int iy = (int)floorf( pixelLocation.y - 0.5f );
        float xDelta = pixelLocation.x - 0.5f - ix;
        float yDelta = pixelLocation.y - 0.5f - iy;
        blend_under( ix, iy, ( ( 1 - xDelta ) * ( 1 - yDelta ) ) * color );
        blend_under( ix, iy + 1, ( ( 1 - xDelta ) * yDelta ) * color );
        blend_under( ix + 1, iy, ( xDelta * ( 1 - yDelta ) ) * color );
        blend_under( ix + 1, iy + 1, ( xDelta * yDelta ) * color );
    }

    void draw_point_bicubic( vector2f pixelLocation, const ColorType& color ) {
        // batty 2005/05/13 -I believe this is 3x3 (Wu?) antialiasing, although
        // I don't know where Ben got it from. Spreads the value
        // over the surrounding 9 pixels, so that the total intensity
        // still sums to 1.

        // (ix,iy) is the nearest neighbor index of the pixelLocation
        int ix = (int)floorf( pixelLocation.x );
        if( -2 <= ix && ix <= ( this->size().xsize + 2 ) ) {

            int iy = (int)floorf( pixelLocation.y );
            if( -2 <= iy && iy <= ( this->size().ysize + 2 ) ) {
                // int rx = (int) frantic::math::round( pixelLocation.x );
                // int ry = (int) frantic::math::round( pixelLocation.y );

                // float xDelta = rx - pixelLocation.x;
                // float yDelta = ry - pixelLocation.y;
                float xDelta = ix - ( pixelLocation.x - 0.5f );
                float yDelta = iy - ( pixelLocation.y - 0.5f );

                float cx0 = ( 0.5f + xDelta ) * ( 0.5f + xDelta ) * 0.5f;
                float cx2 = ( 0.5f - xDelta ) * ( 0.5f - xDelta ) * 0.5f;
                float cx1 = cx0 + cx2 + 2 * ( 0.5f - xDelta ) * ( 0.5f + xDelta );

                float cy0 = ( 0.5f + yDelta ) * ( 0.5f + yDelta ) * 0.5f;
                float cy2 = ( 0.5f - yDelta ) * ( 0.5f - yDelta ) * 0.5f;
                float cy1 = cy0 + cy2 + 2 * ( 0.5f - yDelta ) * ( 0.5f + yDelta );

                blend_under( ix - 1, iy - 1, ( cx0 * cy0 ) * color );
                blend_under( ix - 0, iy - 1, ( cx1 * cy0 ) * color );
                blend_under( ix + 1, iy - 1, ( cx2 * cy0 ) * color );

                blend_under( ix - 1, iy - 0, ( cx0 * cy1 ) * color );
                blend_under( ix - 0, iy - 0, ( cx1 * cy1 ) * color );
                blend_under( ix + 1, iy - 0, ( cx2 * cy1 ) * color );

                blend_under( ix - 1, iy + 1, ( cx0 * cy2 ) * color );
                blend_under( ix - 0, iy + 1, ( cx1 * cy2 ) * color );
                blend_under( ix + 1, iy + 1, ( cx2 * cy2 ) * color );
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------
    void draw_line( const vector2f& screenLocationA, const vector2f& screenLocationB, ColorType color ) {
        float ax = ( screenLocationA.x + 1 ) * 0.5f * this->size().xsize;
        float bx = ( screenLocationB.x + 1 ) * 0.5f * this->size().xsize;
        float ay = ( screenLocationA.y + 1 ) * 0.5f * this->size().ysize;
        float by = ( screenLocationB.y + 1 ) * 0.5f * this->size().ysize;

        draw_line( ax, ay, bx, by, color );
    }

    // Expects positions between 0 and image size
    void draw_line( const vector2& imageLocationA, const vector2& imageLocationB, ColorType color ) {
        draw_line( (float)imageLocationA.x, (float)imageLocationA.y, (float)imageLocationB.x, (float)imageLocationB.y,
                   color );
        // data()[imageLocationA.y * width() + imageLocationA.x] = color3f( 1.0f, 0.0f, 0.0f );
        // data()[imageLocationB.y * width() + imageLocationB.x] = color3f( 0.0f, 0.0f, 1.0f );
    }

    void outline_triangle( const std::vector<vector2f>& vertices, ColorType color ) {
        if( vertices.size() != 3 ) {
            throw std::runtime_error( "outline_triangle: Error, Expecting 3 vertices but received " +
                                      boost::lexical_cast<std::string>( vertices.size() ) + " vertices" );
        }
        draw_line( vertices[0].x, vertices[0].y, vertices[1].x, vertices[1].y, color );
        draw_line( vertices[1].x, vertices[1].y, vertices[2].x, vertices[2].y, color );
        draw_line( vertices[2].x, vertices[2].y, vertices[0].x, vertices[0].y, color );
    }

    void outline_polygon( const std::vector<vector2f>& vertices, ColorType color ) {
        const size_t numVertices = vertices.size();
        for( int i = 0; i < numVertices; ++i ) {
            const vector2f vertex1 = vertices[i];
            const vector2f vertex2 = vertices[( i + 1 ) % numVertices];
            draw_line( vertex1.x, vertex1.y, vertex2.x, vertex2.y, color );
        }
    }

    void fill_triangle( const std::vector<vector2f>& vertices, ColorType color ) {
        if( vertices.size() != 3 ) {
            throw std::runtime_error( "fill_triangle: Error, Expecting 3 vertices but received " +
                                      boost::lexical_cast<std::string>( vertices.size() ) + " vertices" );
        }
        int startingRow;
        std::vector<frantic::graphics2d::simple_run> runs;
        bool boundaryEdges[3] = { true, true, true };
        frantic::graphics2d::rasterize_triangle( vertices[0], vertices[1], vertices[2], this->size().xsize,
                                                 this->size().ysize, boundaryEdges, startingRow, runs );
        for( int rowNum = 0; rowNum < runs.size(); ++rowNum ) {
            const frantic::graphics2d::simple_run& run = runs[rowNum];
            const int ri = startingRow + rowNum;
            for( int ci = run.columnRange.first; ci < run.columnRange.second; ++ci ) {
                this->set_pixel( ci, ri, color );
            }
        }
    }

    void fill_polygon( const std::vector<vector2f>& vertices, ColorType color ) {
        const size_t numVertices = vertices.size();
        if( numVertices == 3 ) {
            fill_triangle( vertices, color );
        } else if( numVertices > 3 ) {
            std::vector<frantic::graphics::vector3> triangleBuffer( numVertices - 2 );
            frantic::geometry::triangulate_polygon( &vertices[0], numVertices, &triangleBuffer[0] );
            for( size_t triangleIndex = 0; triangleIndex < numVertices - 2; ++triangleIndex ) {
                std::vector<vector2f> triangleVerts( 3 );
                const frantic::graphics::vector3& t = triangleBuffer[triangleIndex];
                for( size_t corner = 0; corner < 3; ++corner ) {
                    triangleVerts[corner] = vertices[t[corner]];
                }
                fill_triangle( triangleVerts, color );
            }
        } else {
            throw std::runtime_error( "fill_polygon: Error, Expecting at least 3 vertices but received " +
                                      boost::lexical_cast<std::string>( numVertices ) + " vertices" );
        }
    }

    //---------------------------------------------------------------------------------------------------------
    // Some aggregate image processing algorithms
    //---------------------------------------------------------------------------------------------------------

    // The structure tensor is used in some PDE based image processing algorithms,
    // and is the sum of grad(I)*transpose(grad(I)) over all the color channels.  It
    // is supposed to be a way of measuring the aggregate color variation.
    // Because it is a symmetric 2x2 matrix, we're putting it inside of a 3 component
    // vector.  The resulting vector contains [ sum (dI/dx)^2, sum (dI/dx)*(dI/dy), sum (dI/dy)^2 ].
    void compute_structure_tensor( framebuffer<frantic::graphics::vector3f>& outResult ) {
        outResult.resize( this->size() );
        vector2 coord;
        for( coord.y = 0; coord.y < this->size().ysize; ++coord.y ) {
            for( coord.x = 0; coord.x < this->size().xsize; ++coord.x ) {
                ColorType dIdx2 = this->get_x_partial_derivative_squared( coord );
                ColorType dIdy2 = this->get_y_partial_derivative_squared( coord );
                ColorType dIdx_dIdy = this->get_x_partial_derivative( coord ) * this->get_y_partial_derivative( coord );
                frantic::graphics::vector3f structureTensor;
                structureTensor.x = dIdx2.component_sum();
                structureTensor.y = dIdx_dIdy.component_sum();
                structureTensor.z = dIdy2.component_sum();
                outResult.set_pixel( coord, structureTensor );
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------
    void to_file( const std::string& fileName ) {
        // .cxr files are corona's multichannel exr file format
        if( fileName.substr( fileName.length() - 4, 4 ).compare( ".exr" ) == 0 ||
            fileName.substr( fileName.length() - 4, 4 ).compare( ".cxr" ) == 0 )
            this->to_OpenEXR_file( fileName );
        else if( fileName.substr( fileName.length() - 4, 4 ).compare( ".dpx" ) == 0 )
            this->to_dpx_file( fileName );
    }

    void from_file( const std::string& fileName ) {
        // .cxr files are corona's multichannel exr file format
        if( fileName.substr( fileName.length() - 4, 4 ).compare( ".exr" ) == 0 ||
            fileName.substr( fileName.length() - 4, 4 ).compare( ".cxr" ) == 0 )
            this->from_OpenEXR_file( fileName );
        else if( fileName.substr( fileName.length() - 4, 4 ).compare( ".dpx" ) == 0 )
            this->from_dpx_file( fileName );
    }

    void to_OpenEXR_file( const frantic::tstring& fileName ) const {
#if defined( _WIN32 ) || defined( _WIN64 )
        std::basic_string<TCHAR> tempFile = frantic::win32::ffGetTempFilename();
        image::write_to_OpenEXR( tempFile, this->m_data, this->m_size );

        if( CopyFile( tempFile.c_str(), fileName.c_str(), 0 ) == 0 ) {
            // Give the file move one retry after half a second, in case it was a funky network issue.
            Sleep( 500 );
            if( CopyFile( tempFile.c_str(), fileName.c_str(), 0 ) == 0 ) {
                throw std::runtime_error( "framebuffer.to_OpenEXR_file: Failed to move temporary saved file to \"" +
                                          frantic::strings::to_string( fileName ) +
                                          "\": " + frantic::win32::GetLastErrorMessageA() );
            } else {
                DeleteFile( tempFile.c_str() ); // The copy was a success so delete the temp
            }
        } else {
            DeleteFile( tempFile.c_str() ); // The copy was a success so delete the temp
        }
#else
        image::write_to_OpenEXR( fileName, this->data(), this->size() );
#endif
    }

    void from_OpenEXR_file( const frantic::tstring& fileName ) {
        image::read_from_OpenEXR(
            fileName, this->data(),
            this->m_size ); // passing in the member variable is ok here since its about to be deprecated.
    }

    void to_dpx_file( const std::string& fileName ) const {
        file_io::dpx_file_io fileIo;
        fileIo.copy_from_framebuffer( *this );
#ifdef _WIN32
        std::string tempFile = frantic::win32::ffGetTempFilenameA();
        fileIo.write_file( tempFile );

        if( CopyFile( tempFile.c_str(), fileName.c_str(), 0 ) == 0 ) {
            // Give the file move one retry after half a second, in case it was a funky network issue.
            Sleep( 500 );
            if( CopyFile( tempFile.c_str(), fileName.c_str(), 0 ) == 0 ) {
                throw std::runtime_error( "framebuffer.to_OpenEXR_file: Failed to move temporary saved file to \"" +
                                          fileName + "\": " + frantic::win32::GetLastErrorMessage() );
            } else {
                DeleteFile( tempFile.c_str() ); // The copy was a success so delete the temp
            }
        } else {
            DeleteFile( tempFile.c_str() ); // The copy was a success so delete the temp
        }
#else
        fileIo.write_file( fileName );
#endif
    }

    void from_dpx_file( const std::string& fileName ) {
        file_io::dpx_file_io fileIo;
        fileIo.read_file( fileName );
        fileIo.copy_to_framebuffer( this );
    }

    void from_interpolated( const framebuffer<ColorType>& a, const framebuffer<ColorType>& b, float alpha ) {
        if( a.size() != b.size() )
            throw std::runtime_error(
                "framebuffer.from_interpolated: Tried to interpolate between two buffers of different sizes." );

        set_size( a.size() );

        for( unsigned i = 0; i < this->data().size(); ++i ) {
            this->data()[i] = ( 1 - alpha ) * a.data()[i] + alpha * b.data()[i];
        }
    }

    //---------------------------------------------------------------------------------------------------------

    // If the fusion header file Image.h is included, add this functionality
#ifdef _IMAGE_H_ // Digital Fusion image header
    void from_FusionImage( Image* image, bool addAllChannels = true ) {
        image::from_FusionImage( *this, image, addAllChannels );
    }

    void to_FusionImage( Image* image ) { image::to_FusionImage( image, *this ); }
#endif

#ifdef MAX_VERSION
    void from_3dsMaxBitmap( Bitmap* bm ) {
        // Treat a null pointer as if it were an empty bitmap
        if( bm == 0 ) {
            clear();
            return;
        }
        // Get the bitmap size, and clear this bitmap if it's zero
        size2 bmSize( bm->Width(), bm->Height() );
        if( bmSize.xsize == 0 || bmSize.ysize == 0 ) {
            clear();
            return;
        }

        // Resize this bitmap so we can get the image
        set_size( bmSize );

        image::from_3dsMaxBitmap( bm, m_data, m_size );
    }

    void to_3dsMaxBitmap( Bitmap* bm ) { image::to_3dsMaxBitmap( bm, m_data, m_size ); }
#endif

  private:
    void draw_line( float ax, float ay, float bx, float by, ColorType color );
};

template <typename ColorType>
void framebuffer<ColorType>::draw_line( float x0, float y0, float x1, float y1, ColorType color ) {
    bool steep = std::abs( y1 - y0 ) > std::abs( x1 - x0 );

    if( x0 < 0 )
        x0 = 0;
    if( x0 > this->size().xsize - 1 )
        x0 = this->size().xsize - 1.0f;
    if( x1 < 0 )
        x1 = 0;
    if( x1 > this->size().xsize - 1.0f )
        x1 = this->size().xsize - 1.0f;

    if( y0 < 0 )
        y0 = 0;
    if( y0 > this->size().ysize - 1.0f )
        y0 = this->size().ysize - 1.0f;
    if( y1 < 0 )
        y1 = 0;
    if( y1 > this->size().ysize - 1.0f )
        y1 = this->size().ysize - 1.0f;

    if( steep ) {
        float temp = x0;
        x0 = y0;
        y0 = temp;

        temp = x1;
        x1 = y1;
        y1 = temp;
    }
    if( x0 > x1 ) {
        float temp = x0;
        x0 = x1;
        x1 = temp;

        temp = y0;
        y0 = y1;
        y1 = temp;
    }

    int deltax = (int)( x1 - x0 );
    int deltay = (int)std::abs( y1 - y0 );
    int error = 0;
    int ystep;
    int y = (int)y0;

    if( y0 < y1 )
        ystep = 1;
    else
        ystep = -1;

    for( int x = (int)x0; x <= (int)x1; x++ ) {
        if( steep )
            this->data()[x * this->width() + y] = color;
        else
            this->data()[y * this->width() + x] = color;

        error = error + deltay;

        if( ( 2 * error ) >= deltax ) {
            y = y + ystep;
            error = error - deltax;
        }
    }
}

} // namespace graphics2d
} // namespace frantic
