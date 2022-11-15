// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <wchar.h>

#include <boost/static_assert.hpp>

#include <utf8cpp/utf8.h>

#include <frantic/strings/tstring.hpp>
#include <frantic/strings/utf8.hpp>

#if !defined( _WIN32 ) && !defined( _WIN64 )
#include <langinfo.h>
#endif

namespace frantic {
namespace strings {

bool is_valid_utf8( const char* sbegin, const char* send ) { return utf8::is_valid( sbegin, send ); }

bool is_valid_utf8( const char* s ) {
    if( s ) {
        return utf8::is_valid( s, s + strlen( s ) );
    } else {
        return false; // TODO: what?
    }
}

bool is_valid_utf8( const std::string& s ) { return utf8::is_valid( s.begin(), s.end() ); }

namespace {
#if defined( _WIN32 ) || defined( _WIN64 )
inline bool is_locale_utf8() { return false; }
#else
inline bool is_locale_utf8() {
    // This checks the current locale, not the system locale.
    // My understanding is that most every application should use the system locale.
    // To set the current locale to the system locale, call setlocale( LC_ALL, "" )
    return !strcmp( nl_langinfo( CODESET ), "UTF-8" );
}
#endif
} // anonymous namespace

std::string to_utf8( const char* s ) {
    if( s ) {
        if( is_locale_utf8() && is_valid_utf8( s ) ) {
            // TODO: what should we do when the user's locale is UTF-8, but
            // the string is not UTF-8?
            return s;
        } else {
            // For now I'm doing a roundabout conversion from
            // locale string -> wstring -> UTF-8, because I'm not aware of
            // any simple library support for directly converting from
            // locale string -> UTF-8.
            //
            // Newer C++ standards provide support such as codecvt_utf8
            // and codecvt_utf8_utf16.  Consider using these where they are
            // available.
            return to_utf8( to_wstring( s ) );
        }
    } else {
        return std::string();
    }
}
std::string to_utf8( const std::string& s ) { return to_utf8( s.c_str() ); }

namespace {

#if defined( _WIN32 ) || defined( _WIN64 )

std::string to_utf8( const wchar_t* wsBegin, const wchar_t* wsEnd ) {
    BOOST_STATIC_ASSERT( sizeof( wchar_t ) == 2 );
    std::string result;
    utf8::utf16to8( wsBegin, wsEnd, std::back_inserter( result ) );
    return result;
}
std::wstring wstring_from_utf8( const char* sBegin, const char* sEnd ) {
    BOOST_STATIC_ASSERT( sizeof( wchar_t ) == 2 );
    std::wstring result;
    utf8::utf8to16( sBegin, sEnd, std::back_inserter( result ) );
    return result;
}

#else

std::string to_utf8( const wchar_t* wsBegin, const wchar_t* wsEnd ) {
    BOOST_STATIC_ASSERT( sizeof( wchar_t ) == 4 );
    std::string result;
    utf8::utf32to8( wsBegin, wsEnd, std::back_inserter( result ) );
    return result;
}
std::wstring wstring_from_utf8( const char* sBegin, const char* sEnd ) {
    BOOST_STATIC_ASSERT( sizeof( wchar_t ) == 4 );
    std::wstring result;
    utf8::utf8to32( sBegin, sEnd, std::back_inserter( result ) );
    return result;
}

#endif

} // anonymous namespace

std::string to_utf8( const wchar_t* s ) {
    if( s ) {
        return to_utf8( s, s + wcslen( s ) );
    } else {
        return std::string();
    }
}

std::string to_utf8( const std::wstring& s ) { return to_utf8( s.data(), s.data() + s.size() ); }

std::wstring wstring_from_utf8( const char* s ) {
    if( s ) {
        return wstring_from_utf8( s, s + strlen( s ) );
    } else {
        return std::wstring();
    }
}

std::wstring wstring_from_utf8( const std::string& s ) { return wstring_from_utf8( s.data(), s.data() + s.size() ); }

} // namespace strings
} // namespace frantic
