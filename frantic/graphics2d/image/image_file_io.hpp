// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <vector>

#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/color3h.hpp>
#include <frantic/graphics/color_with_alpha.hpp>

#include <frantic/graphics2d/size2.hpp>

namespace frantic {
namespace graphics2d {
namespace image {

// If we have already included anything that includes ImfIO.h, then enable the OpenEXR functionality
#ifdef INCLUDED_IMF_IO_H

////// NOTE: if these functions appear to be missing, add #include <ImfIO.h> in your file

void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::alpha3f>& data,
                       const frantic::graphics2d::size2& size );
void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::color4f>& data,
                       const frantic::graphics2d::size2& size );
void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::color6f>& data,
                       const frantic::graphics2d::size2& size );
void write_to_OpenEXR( const std::string& path, const std::vector<frantic::graphics::color3f>& data,
                       const frantic::graphics2d::size2& size );
void write_to_OpenEXR( const std::string& path, const std::vector<float>& data,
                       const frantic::graphics2d::size2& size );

void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::alpha3f>& data,
                        frantic::graphics2d::size2& size );
void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::color3f>& data,
                        frantic::graphics2d::size2& size );
void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::color3h>& outData,
                        frantic::graphics2d::size2& outSize );
void read_from_OpenEXR( const std::string& path, std::vector<frantic::graphics::color4f>& data,
                        frantic::graphics2d::size2& size );
void read_from_OpenEXR( const std::string& path, std::vector<float>& data, size2& size );

frantic::graphics2d::size2 get_exr_size( const std::string& filename );
frantic::graphics2d::size2 get_exr_display_size( const std::string& filename );

#if defined( _WIN32 ) || defined( _WIN64 )

// Windows *should* use wchar_t as its native type, but OpenEXR only supports paths via char (presumably via CP_ACP), so
// I provide these overloads. This has the unfortunate property that some valid .exr files cannot be opened due to the
// path requiring UNICODE characters..

inline void write_to_OpenEXR( const std::wstring& path, const std::vector<frantic::graphics::alpha3f>& data,
                              const frantic::graphics2d::size2& size ) {
    write_to_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void write_to_OpenEXR( const std::wstring& path, const std::vector<frantic::graphics::color4f>& data,
                              const frantic::graphics2d::size2& size ) {
    write_to_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void write_to_OpenEXR( const std::wstring& path, const std::vector<frantic::graphics::color6f>& data,
                              const frantic::graphics2d::size2& size ) {
    write_to_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void write_to_OpenEXR( const std::wstring& path, const std::vector<frantic::graphics::color3f>& data,
                              const frantic::graphics2d::size2& size ) {
    write_to_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void write_to_OpenEXR( const std::wstring& path, const std::vector<float>& data,
                              const frantic::graphics2d::size2& size ) {
    write_to_OpenEXR( frantic::strings::to_string( path ), data, size );
}

inline void read_from_OpenEXR( const std::wstring& path, std::vector<frantic::graphics::alpha3f>& data,
                               frantic::graphics2d::size2& size ) {
    read_from_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void read_from_OpenEXR( const std::wstring& path, std::vector<frantic::graphics::color3f>& data,
                               frantic::graphics2d::size2& size ) {
    read_from_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void read_from_OpenEXR( const std::wstring& path, std::vector<frantic::graphics::color3h>& data,
                               frantic::graphics2d::size2& size ) {
    read_from_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void read_from_OpenEXR( const std::wstring& path, std::vector<frantic::graphics::color4f>& data,
                               frantic::graphics2d::size2& size ) {
    read_from_OpenEXR( frantic::strings::to_string( path ), data, size );
}
inline void read_from_OpenEXR( const std::wstring& path, std::vector<float>& data, size2& size ) {
    read_from_OpenEXR( frantic::strings::to_string( path ), data, size );
}

inline frantic::graphics2d::size2 get_exr_size( const std::wstring& filename ) {
    return get_exr_size( frantic::strings::to_string( filename ) );
}
inline frantic::graphics2d::size2 get_exr_display_size( const std::wstring& filename ) {
    return get_exr_display_size( frantic::strings::to_string( filename ) );
}

#endif

template <typename T>
void write_to_dpx( const std::string& path, std::vector<T>& data, frantic::graphics2d::size2& size );
template <typename T>
void read_from_dpx( const std::string& path, std::vector<T>& data, frantic::graphics2d::size2& size );

#endif

} // namespace image
} // namespace graphics2d
} // namespace frantic
