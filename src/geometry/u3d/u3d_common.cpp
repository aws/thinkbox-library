// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/u3d/u3d_common.hpp>
#include <frantic/strings/utf8.hpp>

void frantic::geometry::u3d::write_tstring( frantic::graphics::raw_byte_buffer& buf, const frantic::tstring& s ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    std::string utf8String = frantic::strings::to_utf8( s );
    boost::uint16_t len = static_cast<boost::uint16_t>( utf8String.length() );
    write_value<boost::uint16_t>( buf, len );
    memcpy( buf.add_element( len ), utf8String.data(), len );
#else
    boost::uint16_t len = static_cast<boost::uint16_t>( s.length() );
    write_value<boost::uint16_t>( buf, len );
    memcpy( buf.add_element( len ), s.data(), len );
#endif
}

void frantic::geometry::u3d::write_vector3f( frantic::graphics::raw_byte_buffer& buf,
                                             const frantic::graphics::vector3f& v ) {
    write_value<float>( buf, v.x );
    write_value<float>( buf, v.y );
    write_value<float>( buf, v.z );
}

void frantic::geometry::u3d::write_vector4f( frantic::graphics::raw_byte_buffer& buf,
                                             const frantic::graphics::vector4f& v ) {
    write_value<float>( buf, v.x );
    write_value<float>( buf, v.y );
    write_value<float>( buf, v.z );
    write_value<float>( buf, v.w );
}

boost::uint32_t frantic::geometry::u3d::get_alignment_data_padding( const boost::uint32_t position,
                                                                    const boost::uint32_t alignment ) {
    return ( ( position + ( alignment - 1 ) ) & ~( alignment - 1 ) ) - position;
}

boost::uint32_t frantic::geometry::u3d::get_padding_guess( const frantic::tstring& s ) {
    std::string utf8String = frantic::strings::to_utf8( s );
    return 4 - ( utf8String.length() % 4 );
}

void frantic::geometry::u3d::add_padding( std::ostream& stream, const boost::uint32_t paddingSize ) {
    for( boost::uint32_t i = 0; i < paddingSize; ++i ) {
        stream.put( 0 );
    }
}

void frantic::geometry::u3d::write_color3f( frantic::graphics::raw_byte_buffer& buf,
                                            const frantic::graphics::color3f& c ) {
    write_value<float>( buf, c.r );
    write_value<float>( buf, c.g );
    write_value<float>( buf, c.b );
}

void frantic::geometry::u3d::write_color_rgba_f( frantic::graphics::raw_byte_buffer& buf,
                                                 const frantic::graphics::color_rgba_f& c ) {
    write_value<float>( buf, c.get_r() );
    write_value<float>( buf, c.get_g() );
    write_value<float>( buf, c.get_b() );
    write_value<float>( buf, c.get_a() );
}

void frantic::geometry::u3d::write_transform4f( frantic::graphics::raw_byte_buffer& buf,
                                                const frantic::graphics::transform4f& t ) {
    memcpy( buf.add_element( 64 ), reinterpret_cast<const char*>( &t[0] ), 64u );
}

void frantic::geometry::u3d::write_buffered_block( std::ostream& stream, const frantic::graphics::raw_byte_buffer& buf,
                                                   const std::string& blockName ) {
    boost::uint64_t bufferSize = buf.size();
    stream.write( buf.begin(), bufferSize );

    if( !stream ) {
        std::stringstream ss;
        ss << "u3d::write_buffered_block : Failed to write " << blockName << " block to the .U3D file " << std::endl;
        throw std::runtime_error( ss.str() );
    }
}
