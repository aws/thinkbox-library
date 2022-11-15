// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <frantic/graphics/color3f.hpp>
#if defined( _MSC_VER ) && _MSC_VER > 1200
#include <frantic/graphics/color3h.hpp>
#endif
#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/color_with_alpha.hpp>

#include <frantic/graphics2d/image/algorithms.hpp>
#include <frantic/graphics2d/image/filter.hpp>
#include <frantic/graphics2d/image/pixel.hpp>
#include <frantic/graphics2d/image_channel_base.hpp>

//#include <frantic/graphics2d/image/convert.hpp>

namespace frantic {
namespace graphics2d {
// This should never be done, importing a namespace into another is a bad idea.
// using namespace std;

// Forward declaration of the framebuffer
template <typename FrameBufferType>
class framebuffer;

template <typename DataType>
class image_channel : public image_channel_base {

    // template< typename FrameBufferType >
    //		friend class framebuffer;
    friend class framebuffer<float>;
    friend class framebuffer<graphics::alpha3f>;
    friend class framebuffer<graphics::color3f>;
#if defined( _MSC_VER ) && _MSC_VER > 1200
    friend class framebuffer<graphics::color3h>;
#endif
    friend class framebuffer<graphics::color4f>;
    friend class framebuffer<graphics::color6f>;

    // image_channel implementation

    // Pixel data is laid out:
    //(0,height-1)
    //  ^
    //	|
    //	|
    //	|
    //	|
    //	|
    //	|
    //	|
    //	|
    //	|
    // (0,0)---------------------->(width-1,0)
  protected:
    std::vector<DataType> m_data;

    image_channel( image_channel_base* parent, const std::string& name = "UNNAMED" )
        : image_channel_base( parent, name )
        , m_data() {
        // base checks for valid parent
        if( m_parent->size().get_area() != 0 )
            set_size_helper( m_parent->size() );
    }

    image_channel( image_channel_base* parent, std::vector<DataType> data )
        : image_channel_base( parent )
        , m_data( data ) {
        // base checks for valid parent
        m_size = m_parent->size();
        if( m_size.get_area() != m_data.size() )
            throw std::runtime_error( std::string( "frantic::graphics2d::image_channel() : data size (" ) +
                                      boost::lexical_cast<std::string>( m_data.size() ) +
                                      std::string( ") not the same as parent size (" ) +
                                      boost::lexical_cast<std::string>( parent->size().area() ) + std::string( ")" ) );
    }

    image_channel( image_channel_base* parent, const std::string& name, std::vector<DataType> data )
        : image_channel_base( parent, name )
        , m_data( data ) {
        // base checks for valid parent
        m_size = m_parent->size();
        if( m_size.get_area() != (int)m_data.size() )
            throw std::runtime_error(
                std::string( "frantic::graphics2d::image_channel() : data size (" ) +
                boost::lexical_cast<std::string>( m_data.size() ) + std::string( ") not the same as parent size (" ) +
                boost::lexical_cast<std::string>( parent->size().get_area() ) + std::string( ")" ) );
    }

    void clear_helper() {
        m_size = size2( 0, 0 );
        std::vector<DataType> newData( 0 );
        m_data.swap( newData );
    }

    void set_size_helper( size2 size ) {
        if( size.xsize <= 0 )
            throw std::runtime_error(
                "frantic::graphics2d::image_channel.set_size: size.width() must be greater than zero.  Got size " +
                size.str() );
        if( size.ysize <= 0 )
            throw std::runtime_error(
                "frantic::graphics2d::image_channel.set_size: size.height() must be greater than zero.  Got size " +
                size.str() );

        if( size.get_area() != (int)m_data.size() ) {
            std::vector<DataType> newData( size.get_area() );
            m_data.swap( newData );
        }

        m_size = size;
    }

    void swap_helper( image_channel& rhs ) {
        std::swap( m_size, rhs.m_size );
        m_data.swap( rhs.m_data );
        // NOTE: does not switch channel name
    }

    // OVERLOAD THIS LATER: Overload for when OtherDataType == DataType
    template <typename OtherDataType>
    void copy_from_helper( const image_channel<OtherDataType>& source ) {
        set_size_helper( source.size() );

        for( unsigned int i = 0; i < m_data.size(); ++i )
            m_data[i] = DataType( source.data()[i] );
    }

    void resize_helper( const size2& size ) {
        image_channel<DataType> resized( size );
        image::resize_up<DataType>( m_data, m_size, resized.m_data, size );
        swap( resized );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
  public:
    image_channel( const std::string& name = "UNNAMED" )
        : image_channel_base( name )
        , m_data( 0 ) {}
    image_channel( size2 size )
        : image_channel_base( size )
        , m_data( 0 ) {
        set_size_helper( size );
    }
    image_channel( const std::string& name, size2 size )
        : image_channel_base( name, size )
        , m_data( 0 ) {
        set_size_helper( size );
    }
    image_channel( std::vector<DataType> data, size2 size )
        : image_channel_base( size )
        , m_data( data ) {
        if( m_size.get_area() != (int)m_data.size() )
            throw std::runtime_error( std::string( "frantic::graphics2d::image_channel() : data size (" ) +
                                      boost::lexical_cast<std::string>( m_data.size() ) +
                                      std::string( ") not the same as channel size (" ) +
                                      boost::lexical_cast<std::string>( m_size.get_area() ) + std::string( ")" ) );
    }
    image_channel( const std::string& name, std::vector<DataType> data, size2 size )
        : image_channel_base( name, size )
        , m_data( data ) {
        if( m_size.get_area() != m_data.size() )
            throw std::runtime_error( std::string( "frantic::graphics2d::image_channel() : data size (" ) +
                                      boost::lexical_cast<std::string>( m_data.size() ) +
                                      std::string( ") not the same as channel size (" ) +
                                      boost::lexical_cast<std::string>( m_size.get_area() ) + std::string( ")" ) );
    }
    virtual ~image_channel() {}

    void clear() {
        if( m_parent )
            throw std::runtime_error(
                "frantic::graphics2d::image_channel.clear: clear can only operate on the head channel" );
        clear_helper();
    }

    void set_size( size2 size ) {
        if( m_parent )
            throw std::runtime_error(
                "frantic::graphics2d::image_channel.set_size: set_size can only operate on the head channel" );

        set_size_helper( size );
    }

    void swap( image_channel<DataType>& rhs ) {
        if( m_parent || rhs.m_parent )
            if( m_size != rhs.m_size )
                throw std::runtime_error(
                    "frantic::graphics2d::image_channel.swap: swap can only swap child channels of the same size" );

        swap_helper( rhs );
    }

    template <typename OtherDataType>
    void copy_from( const image_channel<OtherDataType>& source ) {
        if( m_parent )
            if( m_size != source.size() )
                throw std::runtime_error( "frantic::graphics2d::image_channel.copy_from: copy_from can only copy child "
                                          "channels of the same size" );

        copy_from_helper( source );
    }

    void resize( const size2& size ) {
        if( m_parent )
            throw std::runtime_error( "frantic::graphics2d::image_channel.resize: cannot resize child channels" );
        resize_helper( size );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Accessors

    const size2& size() const { return m_size; }

    frantic::graphics2d::size2f size2f() const { return frantic::graphics2d::size2f( m_size ); }

    int width() const { return m_size.xsize; }

    int height() const { return m_size.ysize; }

    const std::vector<DataType>& data() const { return m_data; }

    // NOTE: If you use this, please don't change the size of array...
    std::vector<DataType>& data() { return m_data; }

    const std::string& name() const { return m_name; }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Pixel access as vector

    DataType& operator[]( vector2 pixelCoord ) {
        if( !m_size.is_coord_valid( pixelCoord ) ) {
            throw std::runtime_error( "frantic::graphics2d::image_channel.operator[] : pixelCoord out of range - " +
                                      m_size.str() + ", image size is " + m_size.str() );
        }
        return m_data[m_size.get_index( pixelCoord )];
    }

    const DataType& operator[]( vector2 pixelCoord ) const {
        if( !m_size.is_coord_valid( pixelCoord ) ) {
            throw std::runtime_error( "frantic::graphics2d::image_channel.operator[] : pixelCoord out of range - " +
                                      m_size.str() + ", image size is " + m_size.str() );
        }
        return m_data[m_size.get_index( pixelCoord )];
    }

    DataType& at( int ix, int iy ) {
        vector2 pixelCoord( ix, iy );
        return operator[]( pixelCoord );
    }

    const DataType& at( int ix, int iy ) const {
        vector2 pixelCoord( ix, iy );
        return operator[]( pixelCoord );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Pixel ops

    void set_pixel( vector2 pixel, const DataType& color ) { image::set_pixel( pixel, color, m_data, m_size ); }

    void set_pixel( int ix, int iy, const DataType& color ) { set_pixel( vector2( ix, iy ), color ); }

    void add_pixel( vector2 pixel, const DataType& color ) {
        if( m_size.is_coord_valid( pixel ) ) {
            int index = m_size.get_index( pixel );
            m_data[index] += color;
        }
    }

    void add_pixel( int ix, int iy, const DataType& color ) { add_pixel( vector2( ix, iy ), color ); }

    const DataType& get_pixel( vector2 pixel ) const { return image::get_pixel( pixel, m_data, m_size ); }

    const DataType& get_pixel( int ix, int iy ) const { return image::get_pixel( ix, iy, m_data, m_size ); }

    const DataType& get_pixel_fast( vector2 pixel ) const { return image::get_pixel_fast( pixel, m_data, m_size ); }

    const DataType& get_pixel_fast( int ix, int iy ) const { return image::get_pixel_fast( ix, iy, m_data, m_size ); }

    // Gets a pixel, wrapping along both x and y.
    const DataType& get_pixel_wrap_wrap( int ix, int iy ) const {
        return image::get_pixel_wrap_wrap( vector2( ix, iy ), m_data, m_size );
    }

    const DataType& get_pixel_wrap_wrap( vector2 coords ) const {
        return image::get_pixel_wrap_wrap( coords, m_data, m_size );
    }

    // Gets a pixel, wrapping along x and duplicating the edge along y.
    const DataType& get_pixel_wrap_duplicate( int ix, int iy ) const {
        return image::get_pixel_wrap_duplicate( vector2( ix, iy ), m_data, m_size );
    }

    const DataType& get_pixel_wrap_duplicate( vector2 coords ) const {
        return image::get_pixel_wrap_duplicate( coords, m_data, m_size );
    }

    // Gets a pixel value from normalized coordinates which are in [-1, 1]
    const DataType& get_pixel_nearest_neighbor( vector2f normalizedCoords ) const {
        return image::get_pixel_nearest_neighbor( normalizedCoords, m_data, m_size );
    }

    // Gets a pixel value from normalized coordinates which are in [-1, 1]
    const DataType& get_pixel_nearest_neighbor_int( int x, int y ) const {
        return image::get_pixel_nearest_neighbor_int( x, y, m_data, m_size );
    }

    const DataType& get_pixel_bilinear( vector2f normalizedCoords ) const {
        return image::get_pixel_bilinear( normalizedCoords, m_data, m_size );
    }

    // Gets a pixel value from normalized coordinates which are in [-1,1] and returns the weight
    // for normalization at corners.
    DataType get_pixel_bilinear_return_weight( vector2f normalizedCoords, float& outWeightAccumulator ) const {
        return image::get_pixel_bilinear( normalizedCoords, m_data, m_size, &outWeightAccumulator );
    }

    // Gets a pixel value from normalized coordinates which are in [-1,1] and returns the weight
    // for normalization at edges or corners.
    DataType get_pixel_bicubic_return_weight( vector2f normalizedCoords, float& outWeightAccumulator ) const {
        return image::get_pixel_bicubic( normalizedCoords, m_data, m_size, &outWeightAccumulator );
    }

    // Gets a pixel value from normalized coordinates which are in [-1, 1]
    // This method wraps along both the x and y directions
    DataType get_pixel_bilinear_wrap_wrap( vector2f normalizedCoords ) const {
        return image::get_pixel_bilinear_wrap_wrap( normalizedCoords, m_data, m_size );
    }

    // Gets a pixel value from normalized coordinates which are in [-1, 1]
    // This wraps along the x direction, and mirrors along the y direction
    DataType get_pixel_bilinear_wrap_mirror( vector2f normalizedCoords ) const {
        return image::get_pixel_bilinear_wrap_mirror( normalizedCoords, m_data, m_size );
    }

    DataType get_x_partial_derivative( const vector2& coord ) {
        return image::get_x_partial_derivative( coord, m_data, m_size );
    }

    DataType get_x_partial_derivative_squared( const vector2& coord ) {
        return image::get_x_partial_derivative_squared( coord, m_data, m_size );
    }

    DataType get_y_partial_derivative( const vector2& coord ) {
        return image::get_y_partial_derivative( coord, m_data, m_size );
    }

    DataType get_y_partial_derivative_squared( const vector2& coord ) {
        return image::get_y_partial_derivative_squared( coord, m_data, m_size );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Channel Functions

    void apply_abs();

    void apply_gain( float gain ) {
        for( unsigned int index = 0; index < m_data.size(); ++index ) {
            m_data[index] *= gain;
        }
    }

    void add_image_data( const image_channel<DataType>& channel ) {
        if( m_size != channel.size() )
            throw std::runtime_error( "frantic::graphics2d::image_channel.add_image: Tried to add an image of size " +
                                      m_size.str() + " to an image of size " + channel.size().str() );
        for( unsigned int index = 0; index < m_data.size(); ++index )
            m_data[index] += channel.m_data[index];
    }

    template <typename F>
    void apply_function( const image_channel<DataType>& channel, F func ) {
        if( m_size != channel.size() )
            throw std::runtime_error(
                "frantic::graphics2d::image_channel.apply_function: Tried to apply a function to an image of size " +
                m_size.str() + " to an image of size " + channel.size().str() );
        for( unsigned int index = 0; index < m_data.size(); ++index ) {
            m_data[index] = func( m_data[index], channel.m_data[index] );
        }
    }

    void fill( const DataType& c ) {
        for( unsigned int index = 0; index < m_data.size(); ++index ) {
            m_data[index] = c;
        }
    }

    void fill_under( const DataType& c ) {
        for( unsigned int index = 0; index < m_data.size(); ++index )
            m_data[index].blend_over( c );
    }

    void fill_boundary( const DataType& c ) {
        for( int x = 0; x < m_size.xsize; ++x ) {
            m_data[x] = c;
            m_data[x + m_size.xsize * ( m_size.ysize - 1 )] = c;
        }
        for( int y = 1; y < m_size.ysize - 1; ++y ) {
            m_data[m_size.xsize * y] = c;
            m_data[m_size.xsize * y + m_size.xsize - 1] = c;
        }
    }

    void fill_gradient( const DataType& c1, const DataType& c2 ) {
        double maxDist = sqrt( double( m_size.xsize * m_size.xsize + m_size.ysize * m_size.ysize ) );
        for( int x = 0; x < m_size.xsize; ++x ) {
            for( int y = 0; y < m_size.ysize; ++y ) {
                double dist = sqrt( double( x * x + y * y ) );
                double c1t = dist / maxDist;
                double c2t = 1 - c1t;
                m_data[m_size.get_index( x, y )] = c1 * c1t + c2 * c2t;
            }
        }
    }

    void flood_fill( const vector2& start, const DataType& fillColor, const DataType& boundary,
                     bool overwriteBoundary = false ) {
        image::flood_fill( start, fillColor, boundary, overwriteBoundary, m_data, m_size );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    // Channel Filters

    void apply_filter_x( const std::vector<float> filter ) { image::filter::apply_filter_x( filter, m_data, m_size ); }

    void apply_filter_y( const std::vector<float> filter ) { image::filter::apply_filter_y( filter, m_data, m_size ); }

    void apply_filter( const std::vector<float> filter, const int size ) {
        image::filter::apply_filter( filter, size, m_data, m_size );
    }

    void apply_laplacian_filter( const int size ) { image::filter::apply_laplacian( size, m_data, m_size ); }

    void gaussian_blur( const float sigma ) { image::filter::apply_gaussian( sigma, m_data, m_size ); }
};

//---------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------
// Inlined overloaded functions

template <typename DataType>
inline void image_channel<DataType>::apply_abs() {
    for( unsigned int index = 0; index < m_data.size(); ++index ) {
        m_data[index] = DataType::abs( m_data[index] );
    }
}

template <>
inline void image_channel<float>::apply_abs() {
    for( unsigned int index = 0; index < m_data.size(); ++index ) {
        m_data[index] = fabsf( m_data[index] );
    }
}
} // namespace graphics2d
} // namespace frantic
