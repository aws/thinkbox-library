// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <string>

#if defined( _WIN32 ) || defined( _WIN64 )
#include <tchar.h>
#else
#include <boost/scoped_array.hpp>
#include <cstdlib>
#include <stdexcept>
#endif

// Define FRANTIC_USE_WCHAR to indicate that the library should use
// wchar_t-based tstring.
// I am separating this from UNICODE and _UNICODE in case some
// platform needs to use a different character type for unicode builds
// (for example, char for UTF-8).
#if defined( UNICODE ) || defined( _UNICODE )
#define FRANTIC_USE_WCHAR
#endif

// Define _T macro on Linux.
//
// This macro is used on Windows to define string literals using the
// current character type.
//
// Note: As far as I know Linux uses UTF-8, so you will probably not want
// to use wchar_t.
#if !defined( _WIN32 ) && !defined( _WIN64 )

#ifndef _T
#ifdef FRANTIC_USE_WCHAR
#define _T( x ) L##x
#else
#define _T( x ) x
#endif
#endif

#endif // if !defined(_WIN32) && !defined(_WIN64)

namespace frantic {
namespace strings {

// Define tstring.
// This is also placed in the frantic namespace below.
#ifdef FRANTIC_USE_WCHAR
typedef std::wstring tstring;
#else
typedef std::string tstring;
#endif

// Define tchar.
// This is also placed in the frantic namespace below.
#if defined( _WIN32 ) || defined( _WIN64 )
typedef TCHAR tchar;
#else
#ifdef FRANTIC_USE_WCHAR
typedef wchar_t tchar;
#else
typedef char tchar;
#endif
#endif

// to_string() and to_wstring() functions that don't perform
// any character-type conversion.
// These exist for convenience.
inline std::wstring to_wstring( const wchar_t* s ) { return std::wstring( s ); }

inline const std::wstring& to_wstring( const std::wstring& s ) { return s; }

inline std::string to_string( const char* s ) { return std::string( s ); }

inline const std::string& to_string( const std::string& s ) { return s; }

// to_string() and to_wstring() functions for converting between
// character types.
//
// These are intentionally not defined for non-Windows platforms -- I
// assume that we are using UTF-8, and I want to avoid using wstring
// by accident.
//
// TODO: add support for different code pages.
#if defined( _WIN32 ) || defined( _WIN64 )

std::wstring to_wstring( const char* s );
std::wstring to_wstring( const std::string& s );

std::string to_string( const wchar_t* s );
std::string to_string( const std::wstring& s );

#else

inline std::wstring to_wstring( const char* s ) {
    if( s ) {
        const std::size_t size = mbstowcs( 0, s, 0 );
        if( size == (size_t)-1 ) {
            throw std::runtime_error( "to_wstring: invalid multibyte character in input" );
        } else {
            boost::scoped_array<wchar_t> buffer( new wchar_t[size + 1] );
            const std::size_t result = mbstowcs( buffer.get(), s, size + 1 );
            if( result == size ) {
                return std::wstring( buffer.get(), buffer.get() + size );
            } else if( result == (size_t)-1 ) {
                throw std::runtime_error( "to_wstring: invalid multibyte character in input" );
            } else {
                throw std::runtime_error( "to_wstring: mismatch between expected and real string length" );
            }
        }
    } else {
        return std::wstring();
    }
}
inline std::wstring to_wstring( const std::string& s ) { return to_wstring( s.c_str() ); }

// TODO: should we provide error handling options for these conversions, or
// provide throwing and non-throwing versions?
// The Windows implementation of this function can replace invalid characters
// with a placeholder character such as '?'.
inline std::string to_string( const wchar_t* ws ) {
    if( ws ) {
        const std::size_t size = wcstombs( 0, ws, 0 );
        if( size == (size_t)-1 ) {
            throw std::runtime_error( "to_string: input does not correspond to valid character" );
        } else {
            boost::scoped_array<char> buffer( new char[size + 1] );
            const std::size_t result = wcstombs( buffer.get(), ws, size + 1 );
            if( result == size ) {
                return std::string( buffer.get(), buffer.get() + size );
            } else if( result == (size_t)-1 ) {
                throw std::runtime_error( "to_string: input does not correspond to valid character" );
            } else {
                throw std::runtime_error( "to_string: mismatch between expected and real string length" );
            }
        }
    } else {
        return std::string();
    }
}
inline std::string to_string( const std::wstring& ws ) { return to_string( ws.c_str() ); }

#endif // #if defined(_WIN32) || defined(_WIN64)

// to_tstring() functions for converting to tstring
#ifdef FRANTIC_USE_WCHAR

inline std::wstring to_tstring( const wchar_t* s ) { return s; }
inline const std::wstring& to_tstring( const std::wstring& s ) { return s; }
inline std::wstring to_tstring( const char* s ) { return to_wstring( s ); }
inline std::wstring to_tstring( const std::string& s ) { return to_wstring( s ); }

#else

// These are intentionally not defined for non-Windows platforms -- I
// assume that we are using UTF-8, and I want to avoid using wstring
// by accident.
#if defined( _WIN32 ) || defined( _WIN64 )
inline std::string to_tstring( const wchar_t* s ) { return to_string( s ); }
inline std::string to_tstring( const std::wstring& s ) { return to_string( s ); }
#endif // #if defined(_WIN32) || defined(_WIN64)

inline std::string to_tstring( const char* s ) { return s; }
inline const std::string& to_tstring( const std::string& s ) { return s; }

#endif

} // namespace strings

using frantic::strings::tchar;
using frantic::strings::tstring;

} // namespace frantic
