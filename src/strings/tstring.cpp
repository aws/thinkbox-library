// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/array.hpp>

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace strings {

#if defined( _WIN32 ) || defined( _WIN64 )

namespace {

// character array with short string optimization
template <class CharType>
class char_array {
    CharType* m_buffer;
    boost::array<CharType, 32> m_staticBuffer;
    boost::scoped_array<CharType> m_dynamicBuffer;
    std::size_t m_size;

  public:
    char_array( std::size_t size )
        : m_size( size ) {
        m_buffer = m_staticBuffer.data();
        if( size > m_staticBuffer.size() ) {
            m_dynamicBuffer.reset( new CharType[size] );
            m_buffer = m_dynamicBuffer.get();
        }
    }
    CharType* get() { return m_buffer; }
    std::size_t size() const { return m_size; }
};

std::wstring convert_multibyte_to_wide_string( const char* s, UINT codePage, const wchar_t* errorString = L"<error>" ) {
    if( s ) {
        int size = MultiByteToWideChar( codePage, 0, s, -1, 0, 0 );
        if( size > 0 ) {
            char_array<wchar_t> out( size + 1 );
            int result = MultiByteToWideChar( codePage, 0, s, -1, out.get(), static_cast<int>( out.size() ) );
            if( result == 0 ) {
                return errorString;
            }
            return std::wstring( out.get() );
        } else {
            if( errorString ) {
                return errorString;
            } else {
                throw std::runtime_error( "convert_multibyte_to_wide_string Error code: " +
                                          boost::lexical_cast<std::string>( GetLastError() ) );
            }
        }
    } else {
        return std::wstring();
    }
}

std::string convert_wide_string_to_multibyte( const wchar_t* s, UINT codePage, const char* errorString = "<error>" ) {
    if( s ) {
        int size = WideCharToMultiByte( codePage, 0, s, -1, 0, 0, NULL, 0 );
        if( size > 0 ) {
            char_array<char> out( size + 1 );
            int result = WideCharToMultiByte( codePage, 0, s, -1, out.get(), static_cast<int>( out.size() ), NULL, 0 );
            if( result == 0 ) {
                return errorString;
            }
            return out.get();
        } else {
            if( errorString ) {
                return errorString;
            } else {
                throw std::runtime_error( "convert_wide_string_to_multibyte Error code: " +
                                          boost::lexical_cast<std::string>( GetLastError() ) );
            }
        }
    } else {
        return std::string();
    }
}

} // anonymous namespace

std::wstring to_wstring( const char* s ) {
    return convert_multibyte_to_wide_string( s, CP_ACP, L"<error:to_wstring>" );
}

std::wstring to_wstring( const std::string& s ) { return to_wstring( s.c_str() ); }

std::string to_string( const wchar_t* s ) { return convert_wide_string_to_multibyte( s, CP_ACP, "<error:to_string>" ); }

std::string to_string( const std::wstring& s ) { return to_string( s.c_str() ); }

#endif // #if defined(_WIN32) || defined(_WIN64)

} // namespace strings
} // namespace frantic
