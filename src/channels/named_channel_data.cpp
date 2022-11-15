// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/named_channel_data.hpp>

#include <frantic/locale/locale.hpp>

#include <frantic/logging/logging_level.hpp>

#include <frantic/math/utils.hpp>

#include <frantic/simd/float_v.hpp>

#define BOOST_REGEX_STATIC_LINK
#include <boost/io/ios_state.hpp>
#include <boost/regex.hpp>

#include <iomanip>

using namespace std;
using namespace boost;
using namespace frantic::channels;

using frantic::simd::float_v;

///////////////////////////////
// Here are all the named channel type convertors that are possible
///////////////////////////////

namespace {

// Float conversions where the type size is increasing

void convert_half_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = *reinterpret_cast<const half*>( in );
        out += sizeof( float );
        in += sizeof( half );
    } while( --count );
}

void convert_half_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = *reinterpret_cast<const half*>( in );
        out += sizeof( double );
        in += sizeof( half );
    } while( --count );
}

void convert_float_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = *reinterpret_cast<const float*>( in );
        out += sizeof( double );
        in += sizeof( float );
    } while( --count );
}

// Float conversions where the type size is decreasing

void convert_double_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = float( *reinterpret_cast<const double*>( in ) );
        out += sizeof( float );
        in += sizeof( double );
    } while( --count );
}

void convert_double_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const double*>( in ) );
        out += sizeof( half );
        in += sizeof( double );
    } while( --count );
}

void convert_float_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const float*>( in ) );
        out += sizeof( half );
        in += sizeof( float );
    } while( --count );
}

// Signed int conversions where the type size is increasing

void convert_int8_to_int16( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int16_t*>( out ) = *reinterpret_cast<const int8_t*>( in );
        out += 2;
        in++;
    } while( --count );
}

void convert_int8_to_int32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int32_t*>( out ) = *reinterpret_cast<const int8_t*>( in );
        out += 4;
        in++;
    } while( --count );
}

void convert_int8_to_int64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int64_t*>( out ) = *reinterpret_cast<const int8_t*>( in );
        out += 8;
        in++;
    } while( --count );
}

void convert_int16_to_int32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int32_t*>( out ) = *reinterpret_cast<const int16_t*>( in );
        out += 4;
        in += 2;
    } while( --count );
}

void convert_int16_to_int64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int64_t*>( out ) = *reinterpret_cast<const int16_t*>( in );
        out += 8;
        in += 2;
    } while( --count );
}

void convert_int32_to_int64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int64_t*>( out ) = *reinterpret_cast<const int32_t*>( in );
        out += 8;
        in += 4;
    } while( --count );
}

// Signed int conversions where the type size is decreasing

void convert_int64_to_int32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int32_t*>( out ) = int32_t( *reinterpret_cast<const int64_t*>( in ) );
        out += 4;
        in += 8;
    } while( --count );
}

void convert_int64_to_int16( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int16_t*>( out ) = int16_t( *reinterpret_cast<const int64_t*>( in ) );
        out += 2;
        in += 8;
    } while( --count );
}

void convert_int64_to_int8( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int8_t*>( out ) = int8_t( *reinterpret_cast<const int64_t*>( in ) );
        out++;
        in += 8;
    } while( --count );
}

void convert_int32_to_int16( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int16_t*>( out ) = int16_t( *reinterpret_cast<const int32_t*>( in ) );
        out += 2;
        in += 4;
    } while( --count );
}

void convert_int32_to_int8( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int8_t*>( out ) = int8_t( *reinterpret_cast<const int32_t*>( in ) );
        out++;
        in += 4;
    } while( --count );
}

void convert_int16_to_int8( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int8_t*>( out ) = int8_t( *reinterpret_cast<const int16_t*>( in ) );
        out++;
        in += 2;
    } while( --count );
}

// Unsigned int conversions where the type size is increasing

void convert_uint8_to_uint16( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint16_t*>( out ) = *reinterpret_cast<const uint8_t*>( in );
        out += 2;
        in++;
    } while( --count );
}

void convert_uint8_to_uint32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint32_t*>( out ) = *reinterpret_cast<const uint8_t*>( in );
        out += 4;
        in++;
    } while( --count );
}

void convert_uint8_to_uint64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint64_t*>( out ) = *reinterpret_cast<const uint8_t*>( in );
        out += 8;
        in++;
    } while( --count );
}

void convert_uint16_to_uint32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint32_t*>( out ) = *reinterpret_cast<const uint16_t*>( in );
        out += 4;
        in += 2;
    } while( --count );
}

void convert_uint16_to_uint64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint64_t*>( out ) = *reinterpret_cast<const uint16_t*>( in );
        out += 8;
        in += 2;
    } while( --count );
}

void convert_uint32_to_uint64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint64_t*>( out ) = *reinterpret_cast<const uint32_t*>( in );
        out += 8;
        in += 4;
    } while( --count );
}

// Unsigned int conversions where the type size is decreasing

void convert_uint64_to_uint32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint32_t*>( out ) = uint32_t( *reinterpret_cast<const uint64_t*>( in ) );
        out += 4;
        in += 8;
    } while( --count );
}

void convert_uint64_to_uint16( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint16_t*>( out ) = uint16_t( *reinterpret_cast<const uint64_t*>( in ) );
        out += 2;
        in += 8;
    } while( --count );
}

void convert_uint64_to_uint8( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint8_t*>( out ) = uint8_t( *reinterpret_cast<const uint64_t*>( in ) );
        out++;
        in += 8;
    } while( --count );
}

void convert_uint32_to_uint16( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint16_t*>( out ) = uint16_t( *reinterpret_cast<const uint32_t*>( in ) );
        out += 2;
        in += 4;
    } while( --count );
}

void convert_uint32_to_uint8( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint8_t*>( out ) = uint8_t( *reinterpret_cast<const uint32_t*>( in ) );
        out++;
        in += 4;
    } while( --count );
}

void convert_uint16_to_uint8( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<uint8_t*>( out ) = uint8_t( *reinterpret_cast<const uint16_t*>( in ) );
        out++;
        in += 2;
    } while( --count );
}

// Unsigned int to signed int conversions where the type size is increasing

void convert_uint8_to_int16( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int16_t*>( out ) = *reinterpret_cast<const uint8_t*>( in );
        out += 2;
        in++;
    } while( --count );
}

void convert_uint8_to_int32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int32_t*>( out ) = *reinterpret_cast<const uint8_t*>( in );
        out += 4;
        in++;
    } while( --count );
}

void convert_uint8_to_int64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int64_t*>( out ) = *reinterpret_cast<const uint8_t*>( in );
        out += 8;
        in++;
    } while( --count );
}

void convert_uint16_to_int32( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int32_t*>( out ) = *reinterpret_cast<const uint16_t*>( in );
        out += 4;
        in += 2;
    } while( --count );
}

void convert_uint16_to_int64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int64_t*>( out ) = *reinterpret_cast<const uint16_t*>( in );
        out += 8;
        in += 2;
    } while( --count );
}

void convert_uint32_to_int64( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<int64_t*>( out ) = *reinterpret_cast<const uint32_t*>( in );
        out += 8;
        in += 4;
    } while( --count );
}

// Unsigned int to float conversions where the type size is decreasing

void convert_uint16_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const uint16_t*>( in ) );
        out += sizeof( half );
        in += sizeof( uint16_t );
    } while( --count );
}

void convert_uint32_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const uint32_t*>( in ) );
        out += sizeof( half );
        in += sizeof( uint32_t );
    } while( --count );
}

void convert_uint32_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = float( *reinterpret_cast<const uint32_t*>( in ) );
        out += sizeof( float );
        in += sizeof( uint32_t );
    } while( --count );
}

void convert_uint64_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const uint64_t*>( in ) );
        out += sizeof( half );
        in += sizeof( uint64_t );
    } while( --count );
}

void convert_uint64_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = float( *reinterpret_cast<const uint64_t*>( in ) );
        out += sizeof( float );
        in += sizeof( uint64_t );
    } while( --count );
}

void convert_uint64_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = double( *reinterpret_cast<const uint64_t*>( in ) );
        out += sizeof( double );
        in += sizeof( uint64_t );
    } while( --count );
}

// Unsigned int to float conversions where the type size is increasing

void convert_uint8_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = ( *reinterpret_cast<const uint8_t*>( in ) );
        out += sizeof( half );
        in += sizeof( uint8_t );
    } while( --count );
}

void convert_uint8_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = ( *reinterpret_cast<const uint8_t*>( in ) );
        out += sizeof( float );
        in += sizeof( uint8_t );
    } while( --count );
}

void convert_uint8_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = ( *reinterpret_cast<const uint8_t*>( in ) );
        out += sizeof( double );
        in += sizeof( uint8_t );
    } while( --count );
}

void convert_uint16_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = ( *reinterpret_cast<const uint16_t*>( in ) );
        out += sizeof( float );
        in += sizeof( uint16_t );
    } while( --count );
}

void convert_uint16_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = ( *reinterpret_cast<const uint16_t*>( in ) );
        out += sizeof( double );
        in += sizeof( uint16_t );
    } while( --count );
}

void convert_uint32_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = ( *reinterpret_cast<const uint32_t*>( in ) );
        out += sizeof( double );
        in += sizeof( uint32_t );
    } while( --count );
}

// Signed int to float conversions where the type size is decreasing

void convert_int16_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const int16_t*>( in ) );
        out += sizeof( half );
        in += sizeof( int16_t );
    } while( --count );
}

void convert_int32_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const int32_t*>( in ) );
        out += sizeof( half );
        in += sizeof( int32_t );
    } while( --count );
}

void convert_int32_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = float( *reinterpret_cast<const int32_t*>( in ) );
        out += sizeof( float );
        in += sizeof( int32_t );
    } while( --count );
}

void convert_int64_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = float( *reinterpret_cast<const int64_t*>( in ) );
        out += sizeof( half );
        in += sizeof( int64_t );
    } while( --count );
}

void convert_int64_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = float( *reinterpret_cast<const int64_t*>( in ) );
        out += sizeof( float );
        in += sizeof( int64_t );
    } while( --count );
}

void convert_int64_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = double( *reinterpret_cast<const int64_t*>( in ) );
        out += sizeof( double );
        in += sizeof( int64_t );
    } while( --count );
}

// Signed int to float conversions where the type size is increasing

void convert_int8_to_half( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<half*>( out ) = ( *reinterpret_cast<const int8_t*>( in ) );
        out += sizeof( half );
        in += sizeof( int8_t );
    } while( --count );
}

void convert_int8_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = ( *reinterpret_cast<const int8_t*>( in ) );
        out += sizeof( float );
        in += sizeof( int8_t );
    } while( --count );
}

void convert_int8_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = ( *reinterpret_cast<const int8_t*>( in ) );
        out += sizeof( double );
        in += sizeof( int8_t );
    } while( --count );
}

void convert_int16_to_float( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<float*>( out ) = ( *reinterpret_cast<const int16_t*>( in ) );
        out += sizeof( float );
        in += sizeof( int16_t );
    } while( --count );
}

void convert_int16_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = ( *reinterpret_cast<const int16_t*>( in ) );
        out += sizeof( double );
        in += sizeof( int16_t );
    } while( --count );
}

void convert_int32_to_double( char* out, const char* in, std::size_t count ) {
    do {
        *reinterpret_cast<double*>( out ) = ( *reinterpret_cast<const int32_t*>( in ) );
        out += sizeof( double );
        in += sizeof( int32_t );
    } while( --count );
}

// Copy from the same type

void copy_1byte( char* out, const char* in, std::size_t count ) { memcpy( out, in, count ); }

void copy_2bytes( char* out, const char* in, std::size_t count ) { memcpy( out, in, 2 * count ); }

void copy_4bytes( char* out, const char* in, std::size_t count ) { memcpy( out, in, 4 * count ); }

void copy_8bytes( char* out, const char* in, std::size_t count ) { memcpy( out, in, 8 * count ); }

void copy_16bytes( char* out, const char* in, std::size_t count ) { memcpy( out, in, 16 * count ); }

} // anonymous namespace

namespace frantic {
namespace channels {

///////////////////////////////
// CHANNEL STRING FUNCTIONS
// See the .hpp file for more details.
///////////////////////////////

void construct_channel_string( char* data, const frantic::tstring& init ) {
    char** dataInternal = reinterpret_cast<char**>( data );
    if( init.empty() ) {
        // Set the channel string pointer to NULL, indicating an empty string.
        *dataInternal = 0;
    } else {
        char* dataStr =
            reinterpret_cast<char*>( malloc( 8 + sizeof( frantic::tstring::value_type ) * ( init.size() + 1 ) ) );
        if( dataStr == 0 )
            // throw bad_alloc("construct_channel_string() - Failed to allocate memory for channel_string.");
            throw bad_alloc();
        // Set the channel string pointer to the allocated memory
        *dataInternal = dataStr;

        // Set reference count to 1
        *reinterpret_cast<uint32_t*>( dataStr ) = 1;
        // Set string length to the input string size
        *reinterpret_cast<uint32_t*>( dataStr + 4 ) = (uint32_t)init.size();
        // Copy the string data, including the null terminator, to the destination
        memcpy( dataStr + 8, init.c_str(), sizeof( frantic::tstring::value_type ) * ( init.size() + 1 ) );
    }
}

void construct_channel_string( char* data ) {
    char** dataInternal = reinterpret_cast<char**>( data );
    // Set the channel string pointer to NULL, indicating an empty string.
    *dataInternal = 0;
}

void destruct_channel_string( char* data ) {
    char** dataInternal = reinterpret_cast<char**>( data );
    char* dataStr = *dataInternal;
    if( dataStr ) {
        // NOTE: The following sequence of operations is NOT THREAD SAFE.  For the channel string to be thread
        //       safe, this must be changed.
        // Decrease the reference count
        --( *reinterpret_cast<uint32_t*>( dataStr ) );
        // Free the memory if the reference count has reached zero
        if( *reinterpret_cast<uint32_t*>( dataStr ) == 0 ) {
            // Change the string so that dangling pointers to unchanged memory can't slip through easily
            *reinterpret_cast<uint32_t*>( dataStr + 4 ) = 0;
            free( dataStr );
        }
        // Set the channel string pointer to NULL, so the pointer data doesn't accidentally hang around.
        // Could potentially set it to 0xdeadbeef or something like that in debug mode.
        *dataInternal = 0;
    }
}

void assign_channel_string( char* data, const frantic::tstring& init ) {
    destruct_channel_string( data );
    construct_channel_string( data, init );
}

void copy_channel_string( char* destData, const char* sourceData, std::size_t arity ) {
    const char* const* sourceDataInternal = reinterpret_cast<const char* const*>( sourceData );
    char** destDataInternal = reinterpret_cast<char**>( destData );

    for( ; arity != 0; --arity, ++sourceDataInternal, ++destDataInternal ) {
        destruct_channel_string( reinterpret_cast<char*>( destDataInternal ) );

        if( *sourceDataInternal == 0 ) {
            // An empty string is represented by a NULL pointer.
            *destDataInternal = 0;
        } else {
            const char* dataStr = *sourceDataInternal;
            // Increment the reference count. (Do a const_cast so we can use a const void* for sourceData in the
            // interface.)
            ++( *reinterpret_cast<uint32_t*>( const_cast<char*>( dataStr ) ) );
            // Copy the pointer to the destination string.
            *destDataInternal = const_cast<char*>( dataStr );
        }
    }
}

const frantic::tchar* cstring_from_channel_string( const char* data ) {
    const char* const* dataInternal = reinterpret_cast<const char* const*>( data );
    const char* dataStr = *dataInternal;
    if( dataStr ) {
        // The string data starts 8 bytes in (after a 4 byte reference count and a 4 byte length), and is
        // null-terminated
        return reinterpret_cast<const frantic::tchar*>( dataStr + 8 );
    } else {
        // Since the pointer is NULL, we can return a pointer to this pointer
        // as a C-style string, because it starts with a 0 byte and is hence
        // a null-terminated empty string.
        return reinterpret_cast<const frantic::tchar*>( data );
    }
}

frantic::tstring string_from_channel_string( const char* data ) {
    const char* const* dataInternal = reinterpret_cast<const char* const*>( data );
    const char* dataStr = *dataInternal;
    if( dataStr ) {
        // The string data starts 8 bytes in (after a 4 byte reference count and a 4 byte length), and is
        // null-terminated
        return frantic::tstring( reinterpret_cast<const frantic::tchar*>( dataStr + 8 ),
                                 *reinterpret_cast<const uint32_t*>( dataStr + 4 ) );
    } else {
        // a NULL pointer represents an empty string
        return frantic::tstring();
    }
}

std::size_t length_from_channel_string( const char* data ) {
    const char* const* dataInternal = reinterpret_cast<const char* const*>( data );
    const char* dataStr = *dataInternal;

    if( dataStr ) {
        // The string length is a uint32 quantity 4 bytes into the memory.
        return *reinterpret_cast<const uint32_t*>( dataStr + 4 );
    } else {
        // The empty string is represented by a NULL pointer
        return 0;
    }
}

///////////////////////////////
// Some basic functions dealing with named channel data types
///////////////////////////////

const frantic::tchar* channel_data_type_str( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return _T("int8");
    case data_type_int16:
        return _T("int16");
    case data_type_int32:
        return _T("int32");
    case data_type_int64:
        return _T("int64");
    case data_type_uint8:
        return _T("uint8");
    case data_type_uint16:
        return _T("uint16");
    case data_type_uint32:
        return _T("uint32");
    case data_type_uint64:
        return _T("uint64");
    case data_type_float16:
        return _T("float16");
    case data_type_float32:
        return _T("float32");
    case data_type_float64:
        return _T("float64");
    case data_type_string:
        return _T("string");
    default:
        throw std::runtime_error(
            "channel_data_type_str: Attempted to retrieve the string version of an invalid data type enum " +
            boost::lexical_cast<std::string>( type ) );
    }
}

data_type_t channel_data_type_from_string( const std::string& dataTypeStr ) {
    if( dataTypeStr == "int8" )
        return data_type_int8;
    if( dataTypeStr == "int16" )
        return data_type_int16;
    if( dataTypeStr == "int32" || dataTypeStr == "int" )
        return data_type_int32;
    if( dataTypeStr == "int64" || dataTypeStr == "long" )
        return data_type_int64;
    if( dataTypeStr == "uint8" )
        return data_type_uint8;
    if( dataTypeStr == "uint16" )
        return data_type_uint16;
    if( dataTypeStr == "uint32" || dataTypeStr == "uint" )
        return data_type_uint32;
    if( dataTypeStr == "uint64" )
        return data_type_uint64;
    if( dataTypeStr == "float16" || dataTypeStr == "half" )
        return data_type_float16;
    if( dataTypeStr == "float32" || dataTypeStr == "float" )
        return data_type_float32;
    if( dataTypeStr == "float64" || dataTypeStr == "double" )
        return data_type_float64;
    if( dataTypeStr == "string" )
        return data_type_string;

    return data_type_invalid;
}

std::pair<data_type_t, std::size_t> channel_data_type_and_arity_from_string( const std::string& dataTypeStr ) {
    static boost::regex splitDataTypeStr( "^([a-z0-9]+)\\[(\\d+)\\]$" );
    boost::smatch what;
    if( boost::regex_search( dataTypeStr, what, splitDataTypeStr ) ) {
        size_t arity = atoi( what[2].str().c_str() );
        if( arity > 0 )
            return pair<data_type_t, std::size_t>( channel_data_type_from_string( what[1].str() ), arity );
        else
            return pair<data_type_t, std::size_t>( data_type_invalid, 1 );
    } else {
        return pair<data_type_t, std::size_t>( channel_data_type_from_string( dataTypeStr ), 1 );
    }
}

frantic::tstring channel_data_type_str( std::size_t arity, data_type_t type ) {
    if( arity == 1 )
        return channel_data_type_str( type );
    else
        return frantic::tstring() + channel_data_type_str( type ) + _T("[") +
               boost::lexical_cast<frantic::tstring>( arity ) + _T("]");
}

size_t sizeof_channel_data_type( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
    case data_type_uint8:
        return 1;

    case data_type_int16:
    case data_type_uint16:
    case data_type_float16:
        return 2;

    case data_type_int32:
    case data_type_uint32:
    case data_type_float32:
        return 4;

    case data_type_int64:
    case data_type_uint64:
    case data_type_float64:
        return 8;

    // Strings are represented by a separately managed block of memory, thus the size is pointer size.  Note
    // that this means there will need to be special serialization code.
    case data_type_string:
        return sizeof( char* );

    default:
        throw std::runtime_error(
            "sizeof_channel_data_type: Attempted to determine the size of an invalid data type enum " +
            boost::lexical_cast<std::string>( type ) );
    }
}

// This parses an input string, saving the resulting value to the provided memory buffer.  It's the responsibility
// of the caller to ensure that the output buffer has enough room for the data (up to 8 bytes, currently).
void parse_channel_value_from_string( data_type_t type, const std::string& data, char* outValueBuffer ) {
    switch( type ) {
    case data_type_int8:
    case data_type_int16:
    case data_type_int32:
    case data_type_int64: {
        // TODO: This is MS-specific, need to also do cross platform versions.
        const char* str = data.c_str();
        char* end = 0;
#ifdef _WIN32
        __int64 value = _strtoi64( str, &end, 10 );
#else
        long long value = strtoll( str, &end, 10 );
#endif
        // Make sure that there's only whitespace to the end of the string
        while( *end ) {
            if( !isspace( *end ) )
                throw std::runtime_error( "parse_channel_value_from_string: Could not parse input string \"" + data +
                                          "\" as a " + frantic::strings::to_string( channel_data_type_str( type ) ) +
                                          "." );
            ++end;
        }

        switch( type ) {
        case data_type_int8:
            *reinterpret_cast<boost::int8_t*>( outValueBuffer ) = static_cast<boost::int8_t>( value );
            break;
        case data_type_int16:
            *reinterpret_cast<boost::int16_t*>( outValueBuffer ) = static_cast<boost::int16_t>( value );
            break;
        case data_type_int32:
            *reinterpret_cast<boost::int32_t*>( outValueBuffer ) = static_cast<boost::int32_t>( value );
            break;
        case data_type_int64:
            *reinterpret_cast<boost::int64_t*>( outValueBuffer ) = static_cast<boost::int64_t>( value );
            break;
        default:
            break;
        }

        break;
    }
    case data_type_uint8:
    case data_type_uint16:
    case data_type_uint32:
    case data_type_uint64: {
        // TODO: This is MS-specific, need to also do cross platform versions.
        const char* str = data.c_str();
        char* end;
#ifdef _WIN32
        unsigned __int64 value = _strtoui64( str, &end, 10 );
#else
        unsigned long long value = strtoll( str, &end, 10 );
#endif

        // Make sure that there's only whitespace to the end of the string
        while( *end ) {
            if( !isspace( *end ) )
                throw std::runtime_error( "parse_channel_value_from_string: Could not parse input string \"" + data +
                                          "\" as a " + frantic::strings::to_string( channel_data_type_str( type ) ) +
                                          "." );
            ++end;
        }

        switch( type ) {
        case data_type_uint8:
            *reinterpret_cast<boost::uint8_t*>( outValueBuffer ) = static_cast<boost::uint8_t>( value );
            break;
        case data_type_uint16:
            *reinterpret_cast<boost::uint16_t*>( outValueBuffer ) = static_cast<boost::uint16_t>( value );
            break;
        case data_type_uint32:
            *reinterpret_cast<boost::uint32_t*>( outValueBuffer ) = static_cast<boost::uint32_t>( value );
            break;
        case data_type_uint64:
            *reinterpret_cast<boost::uint64_t*>( outValueBuffer ) = static_cast<boost::uint64_t>( value );
            break;
        default:
            break;
        }

        break;
    }

    case data_type_float16:
    case data_type_float32:
    case data_type_float64: {
        const char* str = data.c_str();
        char* end;
        double value = frantic::locale::strtod_c( str, &end );
        // Make sure that there's only whitespace to the end of the string
        while( *end ) {
            if( !isspace( *end ) )
                throw std::runtime_error( "parse_channel_value_from_string: Could not parse input string \"" + data +
                                          "\" as a " + frantic::strings::to_string( channel_data_type_str( type ) ) +
                                          "." );
            ++end;
        }

        switch( type ) {
        case data_type_float16:
            *reinterpret_cast<half*>( outValueBuffer ) = static_cast<float>( value );
            break;
        case data_type_float32:
            *reinterpret_cast<float*>( outValueBuffer ) = static_cast<float>( value );
            break;
        case data_type_float64:
            *reinterpret_cast<double*>( outValueBuffer ) = static_cast<double>( value );
            break;
        default:
            break;
        }

        break;
    }

    case data_type_string:
        assign_channel_string( outValueBuffer, frantic::strings::to_tstring( data ) );
        break;

    default:
        // TODO: Throw exception?
        break;
    }
}

// This prints a tuple of the given type to a stream
template <class CharType>
void channel_data_type_print_impl( std::basic_ostream<CharType>& out, const std::basic_string<CharType>& separator,
                                   std::size_t arity, data_type_t dataType, const char* rawData ) {
    boost::io::basic_ios_all_saver<CharType> oldState( out );

    // TODO: determine float16 required precision
    if( dataType == data_type_float32 ) {
        out << std::setprecision( std::numeric_limits<float>::digits10 + 1 );
    } else if( dataType == data_type_float64 ) {
        out << std::setprecision( std::numeric_limits<double>::digits10 + 1 );
    }

    std::size_t increment = sizeof_channel_data_type( dataType );
    for( unsigned i = 0; i < arity; ++i ) {
        if( i != 0 )
            out << separator;

        switch( dataType ) {
        case data_type_int8:
            out << int( *reinterpret_cast<const boost::int8_t*>( rawData ) );
            break;
        case data_type_int16:
            out << *reinterpret_cast<const boost::int16_t*>( rawData );
            break;
        case data_type_int32:
            out << *reinterpret_cast<const boost::int32_t*>( rawData );
            break;
        case data_type_int64:
            out << *reinterpret_cast<const boost::int64_t*>( rawData );
            break;
        case data_type_uint8:
            out << unsigned( *reinterpret_cast<const boost::uint8_t*>( rawData ) );
            break;
        case data_type_uint16:
            out << *reinterpret_cast<const boost::uint16_t*>( rawData );
            break;
        case data_type_uint32:
            out << *reinterpret_cast<const boost::uint32_t*>( rawData );
            break;
        case data_type_uint64:
            out << *reinterpret_cast<const boost::uint64_t*>( rawData );
            break;
        case data_type_float16:
            out << *reinterpret_cast<const half*>( rawData );
            break;
        case data_type_float32:
            out << *reinterpret_cast<const float*>( rawData );
            break;
        case data_type_float64:
            out << *reinterpret_cast<const double*>( rawData );
            break;
        case data_type_string:
            out << cstring_from_channel_string( rawData );
            break;
        default:
            out << "unknown";
            break;
        }
        // Cast to a char pointer so the pointer increment is in bytes.
        rawData = reinterpret_cast<const char*>( rawData ) + increment;
    }
}

// This prints a tuple of the given type to a stream
void channel_data_type_print( std::ostream& out, const std::string& separator, std::size_t arity, data_type_t dataType,
                              const char* rawData ) {
    channel_data_type_print_impl( out, separator, arity, dataType, rawData );
}

void channel_data_type_print( std::wostream& out, const std::wstring& separator, std::size_t arity,
                              data_type_t dataType, const char* rawData ) {
    channel_data_type_print_impl( out, separator, arity, dataType, rawData );
}

// This returns a function which does a weighted sum combination of data
channel_weighted_sum_combine_function_t channel_weighted_sum_combine_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return channel_data_type_traits<boost::int8_t>::weighted_sum_combine_general;
    case data_type_int16:
        return channel_data_type_traits<boost::int16_t>::weighted_sum_combine_general;
    case data_type_int32:
        return channel_data_type_traits<boost::int32_t>::weighted_sum_combine_general;
    case data_type_int64:
        return channel_data_type_traits<boost::int64_t>::weighted_sum_combine_general;
    case data_type_uint8:
        return channel_data_type_traits<boost::uint8_t>::weighted_sum_combine_general;
    case data_type_uint16:
        return channel_data_type_traits<boost::uint16_t>::weighted_sum_combine_general;
    case data_type_uint32:
        return channel_data_type_traits<boost::uint32_t>::weighted_sum_combine_general;
    case data_type_uint64:
        return channel_data_type_traits<boost::uint64_t>::weighted_sum_combine_general;
    case data_type_float16:
        return channel_data_type_traits<half>::weighted_sum_combine_general;
    case data_type_float32:
        return channel_data_type_traits<float>::weighted_sum_combine_general;
    case data_type_float64:
        return channel_data_type_traits<double>::weighted_sum_combine_general;
    case data_type_string:
        return channel_data_type_traits<frantic::tstring>::weighted_sum_combine_general;
    // TODO: Should we maybe return NULL, so that the caller can create a more contextual error message?
    default:
        throw std::runtime_error( "channel_barycentric_combine_function: Attempted to get a weighted sum combining "
                                  "function for an invalid data type enum " +
                                  boost::lexical_cast<std::string>( type ) );
    }
}

namespace {

template <class T>
inline void offset_input_weighted_sum_integer_combine_general( const float* weights, std::size_t inputOffset,
                                                               const char* const* data, std::size_t weightCount,
                                                               std::size_t arity, char* out ) {
    // Find the index of the biggest weight
    unsigned biggestWeightIndex = 0;
    for( unsigned i = 1; i < weightCount; ++i ) {
        if( weights[i] > weights[biggestWeightIndex] )
            biggestWeightIndex = i;
    }

    if( out != data[biggestWeightIndex] ) {
        memcpy( out, data[biggestWeightIndex] + inputOffset, arity * sizeof( T ) );
    }
}

template <class T>
inline void offset_input_weighted_sum_float_combine_general( const float* weights, std::size_t inputOffset,
                                                             const char* const* data, std::size_t weightCount,
                                                             std::size_t arity, char* out ) {
    while( arity-- > 0 ) {
        T result = 0;

        for( unsigned i = 0; i < weightCount; ++i )
            result += weights[i] * ( *reinterpret_cast<const T*>( data[i] + inputOffset ) );

        *reinterpret_cast<T*>( out ) = result;

        inputOffset += sizeof( T );
        out += sizeof( T );
    }
}

#ifdef FRANTIC_HAS_SSE2

template <typename T, std::size_t Arity>
struct float_v_storage;

template <typename T>
struct float_v_storage<T, 1> {
    static float_v load( const T* data ) { return float_v::load1( data ); }

    static void store( T* dst, const float_v& src ) { src.store1( dst ); }
};

template <typename T>
struct float_v_storage<T, 2> {
    static float_v load( const T* data ) { return float_v::load2( data ); }

    static void store( T* dst, const float_v& src ) { src.store2( dst ); }
};

template <typename T>
struct float_v_storage<T, 3> {
    static float_v load( const T* data ) { return float_v::load3( data ); }

    static void store( T* dst, const float_v& src ) { src.store3( dst ); }
};

template <typename T>
struct float_v_storage<T, 4> {
    static float_v load( const T* data ) { return float_v::load4( data ); }

    static void store( T* dst, const float_v& src ) { src.store4( dst ); }
};

template <class T, std::size_t Arity>
inline void weighted_sum_float_v( const float* weights, const char* const* data, std::size_t weightCount,
                                  std::size_t offset, char* out ) {
    float_v sum0( _mm_setzero_ps() );

    std::size_t i = 0;
    // If there is an odd number of weights, add one of them now, so we can
    // process them in pairs in the loop below.
    if( weightCount & 1 ) {
        float_v w( weights[i] );
        float_v x = float_v_storage<T, Arity>::load( reinterpret_cast<const T*>( data[i] + offset ) );
        sum0 += w * x;
        ++i;
    }
    if( i < weightCount ) {
        // Process items in pairs, so we can interleave operations for the
        // sake of performance.
        float_v sum1( _mm_setzero_ps() );
        for( ; i < weightCount; i += 2 ) {
            float_v w0( weights[i] );
            float_v w1( weights[i + 1] );
            float_v x0 = float_v_storage<T, Arity>::load( reinterpret_cast<const T*>( data[i] + offset ) );
            float_v x1 = float_v_storage<T, Arity>::load( reinterpret_cast<const T*>( data[i + 1] + offset ) );
            sum0 += w0 * x0;
            sum1 += w1 * x1;
        }
        sum0 += sum1;
    }

    float_v_storage<T, Arity>::store( reinterpret_cast<T*>( out ), sum0 );
}

template <typename T>
inline void offset_input_weighted_sum_float_combine_float_v( const float* weights, std::size_t inputOffset,
                                                             const char* const* data, std::size_t weightCount,
                                                             std::size_t arity, char* out ) {
    BOOST_STATIC_ASSERT( float_v::static_size == 4 );

    while( arity >= 4 ) {
        weighted_sum_float_v<T, 4>( weights, data, weightCount, inputOffset, out );
        inputOffset += 4 * sizeof( T );
        out += 4 * sizeof( T );
        arity -= 4;
    }
    switch( arity ) {
    case 3:
        weighted_sum_float_v<T, 3>( weights, data, weightCount, inputOffset, out );
        break;
    case 2:
        weighted_sum_float_v<T, 2>( weights, data, weightCount, inputOffset, out );
        break;
    case 1:
        weighted_sum_float_v<T, 1>( weights, data, weightCount, inputOffset, out );
        break;
    }
}

#endif // #ifdef FRANTIC_HAS_SSE2

inline void offset_input_weighted_sum_float_combine_half( const float* weights, std::size_t inputOffset,
                                                          const char* const* data, std::size_t weightCount,
                                                          std::size_t arity, char* out ) {
#ifdef FRANTIC_HAS_SSE2
    offset_input_weighted_sum_float_combine_float_v<half>( weights, inputOffset, data, weightCount, arity, out );
#else
    while( arity-- > 0 ) {
        float result = 0;

        for( unsigned i = 0; i < weightCount; ++i )
            result += weights[i] * ( *reinterpret_cast<const half*>( data[i] + inputOffset ) );

        *reinterpret_cast<half*>( out ) = result;

        inputOffset += sizeof( half );
        out += sizeof( half );
    }
#endif
}

inline void offset_input_weighted_sum_float_combine_float( const float* weights, std::size_t inputOffset,
                                                           const char* const* data, std::size_t weightCount,
                                                           std::size_t arity, char* out ) {
#ifdef FRANTIC_HAS_SSE2
    offset_input_weighted_sum_float_combine_float_v<float>( weights, inputOffset, data, weightCount, arity, out );
#else
    offset_input_weighted_sum_float_combine_general<float>( weights, inputOffset, data, weightCount, arity, out );
#endif
}

inline void offset_input_weighted_sum_float_combine_double( const float* weights, std::size_t inputOffset,
                                                            const char* const* data, std::size_t weightCount,
                                                            std::size_t arity, char* out ) {
    offset_input_weighted_sum_float_combine_general<double>( weights, inputOffset, data, weightCount, arity, out );
}

} // anonymous namespace

offset_input_channel_weighted_sum_combine_function_t
offset_input_channel_weighted_sum_combine_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return offset_input_weighted_sum_integer_combine_general<boost::int8_t>;
    case data_type_int16:
        return offset_input_weighted_sum_integer_combine_general<boost::int16_t>;
    case data_type_int32:
        return offset_input_weighted_sum_integer_combine_general<boost::int32_t>;
    case data_type_int64:
        return offset_input_weighted_sum_integer_combine_general<boost::int64_t>;
    case data_type_uint8:
        return offset_input_weighted_sum_integer_combine_general<boost::uint8_t>;
    case data_type_uint16:
        return offset_input_weighted_sum_integer_combine_general<boost::uint16_t>;
    case data_type_uint32:
        return offset_input_weighted_sum_integer_combine_general<boost::uint32_t>;
    case data_type_uint64:
        return offset_input_weighted_sum_integer_combine_general<boost::uint64_t>;
    case data_type_float16:
        return offset_input_weighted_sum_float_combine_half;
    case data_type_float32:
        return offset_input_weighted_sum_float_combine_float;
    case data_type_float64:
        return offset_input_weighted_sum_float_combine_double;
    case data_type_string:
        throw std::runtime_error( "offset_input_channel_weighted_sum_combine_function: "
                                  "Attempted to get a weighted sum combining function for "
                                  "string data.  This is not currently supported." );
    default:
        throw std::runtime_error( "offset_input_channel_weighted_sum_combine_function: "
                                  "Attempted to get a weighted sum combining function for an "
                                  "invalid data type enum " +
                                  boost::lexical_cast<std::string>( type ) );
    }
}

frantic::channels::channel_weighted_sum_combine_function_t
channel_weighted_sum_combine_and_convert_function( frantic::channels::data_type_t sourceType,
                                                   frantic::channels::data_type_t destType,
                                                   const frantic::tstring& channelNameForErrorMessage ) {
    switch( sourceType ) {

    case frantic::channels::data_type_int8:
        switch( destType ) {
        case frantic::channels::data_type_int8:
            return &frantic::channels::channel_data_type_traits<boost::int8_t>::weighted_sum_combine_general;
        case frantic::channels::data_type_int16:
            return &weighted_sum_integer_combine_and_convert<boost::int16_t, boost::int8_t>;
        case frantic::channels::data_type_int32:
            return &weighted_sum_integer_combine_and_convert<boost::int32_t, boost::int8_t>;
        case frantic::channels::data_type_int64:
            return &weighted_sum_integer_combine_and_convert<boost::int64_t, boost::int8_t>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_int16:
        switch( destType ) {
        case frantic::channels::data_type_int8:
            return &weighted_sum_integer_combine_and_convert<boost::int8_t, boost::int16_t>;
        case frantic::channels::data_type_int16:
            return &frantic::channels::channel_data_type_traits<boost::int16_t>::weighted_sum_combine_general;
        case frantic::channels::data_type_int32:
            return &weighted_sum_integer_combine_and_convert<boost::int32_t, boost::int16_t>;
        case frantic::channels::data_type_int64:
            return &weighted_sum_integer_combine_and_convert<boost::int64_t, boost::int16_t>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_int32:
        switch( destType ) {
        case frantic::channels::data_type_int8:
            return &weighted_sum_integer_combine_and_convert<boost::int8_t, boost::int32_t>;
        case frantic::channels::data_type_int16:
            return &weighted_sum_integer_combine_and_convert<boost::int16_t, boost::int32_t>;
        case frantic::channels::data_type_int32:
            return &frantic::channels::channel_data_type_traits<boost::int32_t>::weighted_sum_combine_general;
        case frantic::channels::data_type_int64:
            return &weighted_sum_integer_combine_and_convert<boost::int64_t, boost::int32_t>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_int64:
        switch( destType ) {
        case frantic::channels::data_type_int8:
            return &weighted_sum_integer_combine_and_convert<boost::int8_t, boost::int64_t>;
        case frantic::channels::data_type_int16:
            return &weighted_sum_integer_combine_and_convert<boost::int16_t, boost::int64_t>;
        case frantic::channels::data_type_int32:
            return &weighted_sum_integer_combine_and_convert<boost::int32_t, boost::int64_t>;
        case frantic::channels::data_type_int64:
            return &frantic::channels::channel_data_type_traits<boost::int64_t>::weighted_sum_combine_general;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_uint8:
        switch( destType ) {
        case frantic::channels::data_type_uint8:
            return &frantic::channels::channel_data_type_traits<boost::uint8_t>::weighted_sum_combine_general;
        case frantic::channels::data_type_uint16:
            return &weighted_sum_integer_combine_and_convert<boost::uint16_t, boost::uint8_t>;
        case frantic::channels::data_type_uint32:
            return &weighted_sum_integer_combine_and_convert<boost::uint32_t, boost::uint8_t>;
        case frantic::channels::data_type_uint64:
            return &weighted_sum_integer_combine_and_convert<boost::uint64_t, boost::uint8_t>;
        case frantic::channels::data_type_int16:
            return &weighted_sum_integer_combine_and_convert<boost::int16_t, boost::uint8_t>;
        case frantic::channels::data_type_int32:
            return &weighted_sum_integer_combine_and_convert<boost::int32_t, boost::uint8_t>;
        case frantic::channels::data_type_int64:
            return &weighted_sum_integer_combine_and_convert<boost::int64_t, boost::uint8_t>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_uint16:
        switch( destType ) {
        case frantic::channels::data_type_uint8:
            return &weighted_sum_integer_combine_and_convert<boost::uint8_t, boost::uint16_t>;
        case frantic::channels::data_type_uint16:
            return &frantic::channels::channel_data_type_traits<boost::uint16_t>::weighted_sum_combine_general;
        case frantic::channels::data_type_uint32:
            return &weighted_sum_integer_combine_and_convert<boost::uint32_t, boost::uint16_t>;
        case frantic::channels::data_type_uint64:
            return &weighted_sum_integer_combine_and_convert<boost::uint64_t, boost::uint16_t>;
        case frantic::channels::data_type_int32:
            return &weighted_sum_integer_combine_and_convert<boost::int32_t, boost::uint16_t>;
        case frantic::channels::data_type_int64:
            return &weighted_sum_integer_combine_and_convert<boost::int64_t, boost::uint16_t>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_uint32:
        switch( destType ) {
        case frantic::channels::data_type_uint8:
            return &weighted_sum_integer_combine_and_convert<boost::uint8_t, boost::uint32_t>;
        case frantic::channels::data_type_uint16:
            return &weighted_sum_integer_combine_and_convert<boost::uint16_t, boost::uint32_t>;
        case frantic::channels::data_type_uint32:
            return &frantic::channels::channel_data_type_traits<boost::uint32_t>::weighted_sum_combine_general;
        case frantic::channels::data_type_uint64:
            return &weighted_sum_integer_combine_and_convert<boost::uint64_t, boost::uint32_t>;
        case frantic::channels::data_type_int64:
            return &weighted_sum_integer_combine_and_convert<boost::int64_t, boost::uint32_t>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_uint64:
        switch( destType ) {
        case frantic::channels::data_type_uint8:
            return &weighted_sum_integer_combine_and_convert<boost::uint8_t, boost::uint64_t>;
        case frantic::channels::data_type_uint16:
            return &weighted_sum_integer_combine_and_convert<boost::uint16_t, boost::uint64_t>;
        case frantic::channels::data_type_uint32:
            return &weighted_sum_integer_combine_and_convert<boost::uint32_t, boost::uint64_t>;
        case frantic::channels::data_type_uint64:
            return &frantic::channels::channel_data_type_traits<boost::uint64_t>::weighted_sum_combine_general;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_float16:
        switch( destType ) {
        case frantic::channels::data_type_float16:
            return &frantic::channels::channel_data_type_traits<half>::weighted_sum_combine_general;
        case frantic::channels::data_type_float32:
            return &weighted_sum_float_combine_and_convert<float, float, half>;
        case frantic::channels::data_type_float64:
            return &weighted_sum_float_combine_and_convert<double, double, half>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_float32:
        switch( destType ) {
        case frantic::channels::data_type_float16:
            return &weighted_sum_float_combine_and_convert<half, float, float>;
        case frantic::channels::data_type_float32:
            return &frantic::channels::channel_data_type_traits<float>::weighted_sum_combine_general;
        case frantic::channels::data_type_float64:
            return &weighted_sum_float_combine_and_convert<double, double, float>;
        default:
            break;
        }
        break;

    case frantic::channels::data_type_float64:
        switch( destType ) {
        case frantic::channels::data_type_float16:
            return &weighted_sum_float_combine_and_convert<half, float, double>;
        case frantic::channels::data_type_float32:
            return &weighted_sum_float_combine_and_convert<float, double, double>;
        case frantic::channels::data_type_float64:
            return &frantic::channels::channel_data_type_traits<double>::weighted_sum_combine_general;
        default:
            break;
        }
        break;

    default:
        break;
    }

    throw std::runtime_error( "channel_weighted_sum_combine_and_convert_function() - Conversion of channel \"" +
                              frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type " +
                              frantic::strings::to_string( channel_data_type_str( sourceType ) ) + " to type " +
                              frantic::strings::to_string( channel_data_type_str( destType ) ) +
                              " is not an allowed conversion." );
}

// This returns a function which does a weighted increment of data
channel_weighted_increment_function_t channel_weighted_increment_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return weighted_increment_general<boost::int8_t>;
    case data_type_int16:
        return weighted_increment_general<boost::int16_t>;
    case data_type_int32:
        return weighted_increment_general<boost::int32_t>;
    case data_type_int64:
        return weighted_increment_general<boost::int64_t>;
    case data_type_uint8:
        return weighted_increment_general<boost::uint8_t>;
    case data_type_uint16:
        return weighted_increment_general<boost::uint16_t>;
    case data_type_uint32:
        return weighted_increment_general<boost::uint32_t>;
    case data_type_uint64:
        return weighted_increment_general<boost::uint64_t>;
    case data_type_float16:
        return weighted_increment_general<half>;
    case data_type_float32:
        return weighted_increment_general<float>;
    case data_type_float64:
        return weighted_increment_general<double>;
    case data_type_string:
        return channel_data_type_traits<frantic::tstring>::weighted_increment_general;
    // TODO: Should we maybe return NULL, so that the caller can create a more contextual error message?
    default:
        throw std::runtime_error( "channel_barycentric_combine_function: Attempted to get a weighted sum combining "
                                  "function for an invalid data type enum " +
                                  boost::lexical_cast<std::string>( type ) );
    }
}

channel_scale_function_t channel_scale_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return channel_scale_general<boost::int8_t>;
    case data_type_int16:
        return channel_scale_general<boost::int16_t>;
    case data_type_int32:
        return channel_scale_general<boost::int32_t>;
    case data_type_int64:
        return channel_scale_general<boost::int64_t>;
    case data_type_uint8:
        return channel_scale_general<boost::uint8_t>;
    case data_type_uint16:
        return channel_scale_general<boost::uint16_t>;
    case data_type_uint32:
        return channel_scale_general<boost::uint32_t>;
    case data_type_uint64:
        return channel_scale_general<boost::uint64_t>;
    case data_type_float16:
        return channel_scale_general<half>;
    case data_type_float32:
        return channel_scale_general<float>;
    case data_type_float64:
        return channel_scale_general<double>;
    case data_type_string:
        throw std::runtime_error( "channel_scale_function: Attempted to get a channel scaling function for the string "
                                  "type.  This is not currently supported." );
    // TODO: Should we maybe return NULL, so that the caller can create a more contextual error message?
    default:
        throw std::runtime_error(
            "channel_scale_function: Attempted to get a channel scaling function for an invalid data type enum " +
            boost::lexical_cast<std::string>( type ) );
    }
}

namespace channel_range_map_general_detail {
template <class DataType>
struct intermediate {
    typedef double type;
};
// avoid a warning by casting half to float instead of to double
template <>
struct intermediate<half> {
    typedef float type;
};
} // namespace channel_range_map_general_detail
template <class T>
inline void channel_range_map_general( double fromLB, double fromUB, double toLB, double toUB, const char* data,
                                       std::size_t arity, char* out ) {
    typedef typename channel_scale_general_detail::intermediate<T>::type intermediate_type;

    std::size_t offset = 0;
    while( arity-- > 0 ) {
        T& outPrimitive = reinterpret_cast<T*>( out )[offset];
        outPrimitive = frantic::math::linearConvertRange(
            reinterpret_cast<const T*>( data )[offset], T( intermediate_type( fromLB ) ),
            T( intermediate_type( fromUB ) ), T( intermediate_type( toLB ) ), T( intermediate_type( toUB ) ) );
        ++offset;
    }
}

channel_range_map_function_t channel_range_map_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return channel_range_map_general<boost::int8_t>;
    case data_type_int16:
        return channel_range_map_general<boost::int16_t>;
    case data_type_int32:
        return channel_range_map_general<boost::int32_t>;
    case data_type_int64:
        return channel_range_map_general<boost::int64_t>;
    case data_type_uint8:
        return channel_range_map_general<boost::uint8_t>;
    case data_type_uint16:
        return channel_range_map_general<boost::uint16_t>;
    case data_type_uint32:
        return channel_range_map_general<boost::uint32_t>;
    case data_type_uint64:
        return channel_range_map_general<boost::uint64_t>;
    case data_type_float16:
        return channel_range_map_general<half>;
    case data_type_float32:
        return channel_range_map_general<float>;
    case data_type_float64:
        return channel_range_map_general<double>;
    case data_type_string:
        throw std::runtime_error( "channel_scale_function: Attempted to get a channel scaling function for the string "
                                  "type.  This is not currently supported." );
    // TODO: Should we maybe return NULL, so that the caller can create a more contextual error message?
    default:
        throw std::runtime_error(
            "channel_scale_function: Attempted to get a channel scaling function for an invalid data type enum " +
            boost::lexical_cast<std::string>( type ) );
    }
}

namespace detail {
template <class T>
struct channel_op {
    static void add( const char* srcLHS, const char* srcRHS, std::size_t arity, char* dest ) {
        for( std::size_t i = 0; i < arity; ++i )
            reinterpret_cast<T*>( dest )[i] =
                reinterpret_cast<const T*>( srcLHS )[i] + reinterpret_cast<const T*>( srcRHS )[i];
    }
    static void sub( const char* srcLHS, const char* srcRHS, std::size_t arity, char* dest ) {
        for( std::size_t i = 0; i < arity; ++i )
            reinterpret_cast<T*>( dest )[i] =
                reinterpret_cast<const T*>( srcLHS )[i] - reinterpret_cast<const T*>( srcRHS )[i];
    }
    static void mul( const char* srcLHS, const char* srcRHS, std::size_t arity, char* dest ) {
        for( std::size_t i = 0; i < arity; ++i )
            reinterpret_cast<T*>( dest )[i] =
                reinterpret_cast<const T*>( srcLHS )[i] * reinterpret_cast<const T*>( srcRHS )[i];
    }
    static void max( const char* srcLHS, const char* srcRHS, std::size_t arity, char* dest ) {
        for( std::size_t i = 0; i < arity; ++i )
            reinterpret_cast<T*>( dest )[i] =
                std::max<T>( reinterpret_cast<const T*>( srcLHS )[i], reinterpret_cast<const T*>( srcRHS )[i] );
    }
    static void min( const char* srcLHS, const char* srcRHS, std::size_t arity, char* dest ) {
        for( std::size_t i = 0; i < arity; ++i )
            reinterpret_cast<T*>( dest )[i] =
                std::min<T>( reinterpret_cast<const T*>( srcLHS )[i], reinterpret_cast<const T*>( srcRHS )[i] );
    }
};
} // namespace detail

channel_binary_function_t channel_addition_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return detail::channel_op<boost::int8_t>::add;
    case data_type_int16:
        return detail::channel_op<boost::int16_t>::add;
    case data_type_int32:
        return detail::channel_op<boost::int32_t>::add;
    case data_type_int64:
        return detail::channel_op<boost::int64_t>::add;
    case data_type_uint8:
        return detail::channel_op<boost::uint8_t>::add;
    case data_type_uint16:
        return detail::channel_op<boost::uint16_t>::add;
    case data_type_uint32:
        return detail::channel_op<boost::uint32_t>::add;
    case data_type_uint64:
        return detail::channel_op<boost::uint64_t>::add;
    case data_type_float16:
        return detail::channel_op<half>::add;
    case data_type_float32:
        return detail::channel_op<float>::add;
    case data_type_float64:
        return detail::channel_op<double>::add;
    default:
        throw std::runtime_error( "There is no addition function for the data type" );
    }
}

channel_binary_function_t channel_subtraction_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return detail::channel_op<boost::int8_t>::sub;
    case data_type_int16:
        return detail::channel_op<boost::int16_t>::sub;
    case data_type_int32:
        return detail::channel_op<boost::int32_t>::sub;
    case data_type_int64:
        return detail::channel_op<boost::int64_t>::sub;
    case data_type_uint8:
        return detail::channel_op<boost::uint8_t>::sub;
    case data_type_uint16:
        return detail::channel_op<boost::uint16_t>::sub;
    case data_type_uint32:
        return detail::channel_op<boost::uint32_t>::sub;
    case data_type_uint64:
        return detail::channel_op<boost::uint64_t>::sub;
    case data_type_float16:
        return detail::channel_op<half>::sub;
    case data_type_float32:
        return detail::channel_op<float>::sub;
    case data_type_float64:
        return detail::channel_op<double>::sub;
    default:
        throw std::runtime_error( "There is no subtraction function for the data type" );
    }
}

channel_binary_function_t channel_multiplication_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return detail::channel_op<boost::int8_t>::mul;
    case data_type_int16:
        return detail::channel_op<boost::int16_t>::mul;
    case data_type_int32:
        return detail::channel_op<boost::int32_t>::mul;
    case data_type_int64:
        return detail::channel_op<boost::int64_t>::mul;
    case data_type_uint8:
        return detail::channel_op<boost::uint8_t>::mul;
    case data_type_uint16:
        return detail::channel_op<boost::uint16_t>::mul;
    case data_type_uint32:
        return detail::channel_op<boost::uint32_t>::mul;
    case data_type_uint64:
        return detail::channel_op<boost::uint64_t>::mul;
    case data_type_float16:
        return detail::channel_op<half>::mul;
    case data_type_float32:
        return detail::channel_op<float>::mul;
    case data_type_float64:
        return detail::channel_op<double>::mul;
    default:
        throw std::runtime_error( "There is no multiplication function for the data type" );
    }
}

channel_binary_function_t channel_maximum_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return detail::channel_op<boost::int8_t>::max;
    case data_type_int16:
        return detail::channel_op<boost::int16_t>::max;
    case data_type_int32:
        return detail::channel_op<boost::int32_t>::max;
    case data_type_int64:
        return detail::channel_op<boost::int64_t>::max;
    case data_type_uint8:
        return detail::channel_op<boost::uint8_t>::max;
    case data_type_uint16:
        return detail::channel_op<boost::uint16_t>::max;
    case data_type_uint32:
        return detail::channel_op<boost::uint32_t>::max;
    case data_type_uint64:
        return detail::channel_op<boost::uint64_t>::max;
    case data_type_float16:
        return detail::channel_op<half>::max;
    case data_type_float32:
        return detail::channel_op<float>::max;
    case data_type_float64:
        return detail::channel_op<double>::max;
    default:
        throw std::runtime_error( "There is no maximum function for the data type" );
    }
}

channel_binary_function_t channel_minimum_function( data_type_t type ) {
    switch( type ) {
    case data_type_int8:
        return detail::channel_op<boost::int8_t>::min;
    case data_type_int16:
        return detail::channel_op<boost::int16_t>::min;
    case data_type_int32:
        return detail::channel_op<boost::int32_t>::min;
    case data_type_int64:
        return detail::channel_op<boost::int64_t>::min;
    case data_type_uint8:
        return detail::channel_op<boost::uint8_t>::min;
    case data_type_uint16:
        return detail::channel_op<boost::uint16_t>::min;
    case data_type_uint32:
        return detail::channel_op<boost::uint32_t>::min;
    case data_type_uint64:
        return detail::channel_op<boost::uint64_t>::min;
    case data_type_float16:
        return detail::channel_op<half>::min;
    case data_type_float32:
        return detail::channel_op<float>::min;
    case data_type_float64:
        return detail::channel_op<double>::min;
    default:
        throw std::runtime_error( "There is no minimum function for that data type(" +
                                  boost::lexical_cast<std::string>( type ) + ")." );
    }
}

///////////////////////////////
// Some functions dealing with conversions of named channel data types
///////////////////////////////

// This function returns a convertor function which can convert an array of the source type
// into an array of the destination type.
channel_type_convertor_function_t
get_channel_type_convertor_function( data_type_t sourceType, data_type_t destType,
                                     const frantic::tstring& channelNameForErrorMessage ) {
    // Special case string copying, because it can't be done with memcpy
    if( sourceType == data_type_string && destType == data_type_string )
        return channels::copy_channel_string;

    // Deal with matching types first
    if( sourceType == destType ) {
        switch( sizeof_channel_data_type( sourceType ) ) {
        case 1:
            return &copy_1byte;
        case 2:
            return &copy_2bytes;
        case 4:
            return &copy_4bytes;
        case 8:
            return &copy_8bytes;
        case 16:
            return &copy_16bytes;
        default:
            throw std::runtime_error( "get_channel_type_convertor_function: The conversion type for channel \"" +
                                      frantic::strings::to_string( channelNameForErrorMessage ) +
                                      "\" had an unexpected size." );
        }
    }

    switch( sourceType ) {

    case data_type_int8:
        switch( destType ) {
        case data_type_int16:
            return &convert_int8_to_int16;
        case data_type_int32:
            return &convert_int8_to_int32;
        case data_type_int64:
            return &convert_int8_to_int64;
        case data_type_float16:
            return &convert_int8_to_half;
        case data_type_float32:
            return &convert_int8_to_float;
        case data_type_float64:
            return &convert_int8_to_double;
        default:
            break;
        }
        break;
    case data_type_int16:
        switch( destType ) {
        case data_type_int8:
            return &convert_int16_to_int8;
        case data_type_int32:
            return &convert_int16_to_int32;
        case data_type_int64:
            return &convert_int16_to_int64;
        case data_type_float16:
            return &convert_int16_to_half;
        case data_type_float32:
            return &convert_int16_to_float;
        case data_type_float64:
            return &convert_int16_to_double;
        default:
            break;
        }
        break;
    case data_type_int32:
        switch( destType ) {
        case data_type_int8:
            return &convert_int32_to_int8;
        case data_type_int16:
            return &convert_int32_to_int16;
        case data_type_int64:
            return &convert_int32_to_int64;
        case data_type_float16:
            return &convert_int32_to_half;
        case data_type_float32:
            return &convert_int32_to_float;
        case data_type_float64:
            return &convert_int32_to_double;
        default:
            break;
        }
        break;
    case data_type_int64:
        switch( destType ) {
        case data_type_int8:
            return &convert_int64_to_int8;
        case data_type_int16:
            return &convert_int64_to_int16;
        case data_type_int32:
            return &convert_int64_to_int32;
        case data_type_float16:
            return &convert_int64_to_half;
        case data_type_float32:
            return &convert_int64_to_float;
        case data_type_float64:
            return &convert_int64_to_double;
        default:
            break;
        }
        break;
    case data_type_uint8:
        switch( destType ) {
        case data_type_uint16:
            return &convert_uint8_to_uint16;
        case data_type_uint32:
            return &convert_uint8_to_uint32;
        case data_type_uint64:
            return &convert_uint8_to_uint64;
        case data_type_int16:
            return &convert_uint8_to_int16;
        case data_type_int32:
            return &convert_uint8_to_int32;
        case data_type_int64:
            return &convert_uint8_to_int64;
        case data_type_float16:
            return &convert_uint8_to_half;
        case data_type_float32:
            return &convert_uint8_to_float;
        case data_type_float64:
            return &convert_uint8_to_double;
        default:
            break;
        }
        break;
    case data_type_uint16:
        switch( destType ) {
        case data_type_uint8:
            return &convert_uint16_to_uint8;
        case data_type_uint32:
            return &convert_uint16_to_uint32;
        case data_type_uint64:
            return &convert_uint16_to_uint64;
        case data_type_int32:
            return &convert_uint16_to_int32;
        case data_type_int64:
            return &convert_uint16_to_int64;
        case data_type_float16:
            return &convert_uint16_to_half;
        case data_type_float32:
            return &convert_uint16_to_float;
        case data_type_float64:
            return &convert_uint16_to_double;
        default:
            break;
        }
        break;
    case data_type_uint32:
        switch( destType ) {
        case data_type_uint8:
            return &convert_uint32_to_uint8;
        case data_type_uint16:
            return &convert_uint32_to_uint16;
        case data_type_uint64:
            return &convert_uint32_to_uint64;
        case data_type_int64:
            return &convert_uint32_to_int64;
        case data_type_float16:
            return &convert_uint32_to_half;
        case data_type_float32:
            return &convert_uint32_to_float;
        case data_type_float64:
            return &convert_uint32_to_double;
        default:
            break;
        }
        break;
    case data_type_uint64:
        switch( destType ) {
        case data_type_uint8:
            return &convert_uint64_to_uint8;
        case data_type_uint16:
            return &convert_uint64_to_uint16;
        case data_type_uint32:
            return &convert_uint64_to_uint32;
        case data_type_float16:
            return &convert_uint64_to_half;
        case data_type_float32:
            return &convert_uint64_to_float;
        case data_type_float64:
            return &convert_uint64_to_double;
        default:
            break;
        }
        break;
    case data_type_float16:
        switch( destType ) {
        case data_type_float32:
            return &convert_half_to_float;
        case data_type_float64:
            return &convert_half_to_double;
        default:
            break;
        }
        break;
    case data_type_float32:
        switch( destType ) {
        case data_type_float16:
            return &convert_float_to_half;
        case data_type_float64:
            return &convert_float_to_double;
        default:
            break;
        }
        break;
    case data_type_float64:
        switch( destType ) {
        case data_type_float16:
            return &convert_double_to_half;
        case data_type_float32:
            return &convert_double_to_float;
        default:
            break;
        }
        break;
    default:
        break;
    } // switch( sourceType )

    throw runtime_error( "get_channel_type_convertor_function() - Conversion of channel \"" +
                         frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type " +
                         frantic::strings::to_string( channel_data_type_str( sourceType ) ) + " to type " +
                         frantic::strings::to_string( channel_data_type_str( destType ) ) +
                         " is not an allowed conversion." );
}

data_type_t promote_types( data_type_t t1, data_type_t t2 ) {
    std::size_t t1Size = sizeof_channel_data_type( t1 );
    std::size_t t2Size = sizeof_channel_data_type( t2 );

    if( t1 == t2 )
        return t1;
    else {
        if( t1Size == t2Size ) {
            throw std::runtime_error( std::string( "promote_types() - The types: " ) +
                                      frantic::strings::to_string( channel_data_type_str( t1 ) ) + " and " +
                                      frantic::strings::to_string( channel_data_type_str( t2 ) ) +
                                      " are incompatible." );
        } else if( t1Size < t2Size ) {
            if( is_channel_data_type_float( t1 ) && !is_channel_data_type_float( t2 ) )
                throw std::runtime_error( std::string( "promote_types() - The types: " ) +
                                          frantic::strings::to_string( channel_data_type_str( t1 ) ) + " and " +
                                          frantic::strings::to_string( channel_data_type_str( t2 ) ) +
                                          " are incompatible." );
            if( !is_channel_data_type_unsigned( t1 ) && is_channel_data_type_unsigned( t2 ) )
                throw std::runtime_error( std::string( "promote_types() - The types: " ) +
                                          frantic::strings::to_string( channel_data_type_str( t1 ) ) + " and " +
                                          frantic::strings::to_string( channel_data_type_str( t2 ) ) +
                                          " are incompatible." );
            return t2;
        } else {
            if( is_channel_data_type_float( t2 ) && !is_channel_data_type_float( t1 ) )
                throw std::runtime_error( std::string( "promote_types() - The types: " ) +
                                          frantic::strings::to_string( channel_data_type_str( t1 ) ) + " and " +
                                          frantic::strings::to_string( channel_data_type_str( t2 ) ) +
                                          " are incompatible." );
            if( !is_channel_data_type_unsigned( t2 ) && is_channel_data_type_unsigned( t1 ) )
                throw std::runtime_error( std::string( "promote_types() - The types: " ) +
                                          frantic::strings::to_string( channel_data_type_str( t1 ) ) + " and " +
                                          frantic::strings::to_string( channel_data_type_str( t2 ) ) +
                                          " are incompatible." );
            return t1;
        }
    }
}

// NOTE: The MS compiler collapses together identical functions.  The functions convert_uint32_to_uint64 and
// convert_uint32_to_int64
//       produce identical assembly language, so the function pointer for them ends up being the same.  Keep this in
//       mind when looking
//       at debug dumps of a channel_map_adaptor.
const char* get_channel_type_convertor_debug_string( channel_type_convertor_function_t ptc ) {
    if( ptc == &convert_half_to_float )
        return "convert_half_to_float";
    if( ptc == &convert_half_to_double )
        return "convert_half_to_double";
    if( ptc == &convert_float_to_double )
        return "convert_float_to_double";
    if( ptc == &convert_double_to_float )
        return "convert_double_to_float";
    if( ptc == &convert_double_to_half )
        return "convert_double_to_half";
    if( ptc == &convert_float_to_half )
        return "convert_float_to_half";
    if( ptc == &convert_int8_to_int16 )
        return "convert_int8_to_int16";
    if( ptc == &convert_int8_to_int32 )
        return "convert_int8_to_int32";
    if( ptc == &convert_int8_to_int64 )
        return "convert_int8_to_int64";
    if( ptc == &convert_int16_to_int8 )
        return "convert_int16_to_int8";
    if( ptc == &convert_int16_to_int32 )
        return "convert_int16_to_int32";
    if( ptc == &convert_int16_to_int64 )
        return "convert_int16_to_int64";
    if( ptc == &convert_int32_to_int8 )
        return "convert_int32_to_int8";
    if( ptc == &convert_int32_to_int16 )
        return "convert_int32_to_int16";
    if( ptc == &convert_int32_to_int64 )
        return "convert_int32_to_int64";
    if( ptc == &convert_int64_to_int8 )
        return "convert_int64_to_int8";
    if( ptc == &convert_int64_to_int16 )
        return "convert_int64_to_int16";
    if( ptc == &convert_int64_to_int32 )
        return "convert_int64_to_int32";
    if( ptc == &convert_uint8_to_uint16 )
        return "convert_uint8_to_uint16";
    if( ptc == &convert_uint8_to_uint32 )
        return "convert_uint8_to_uint32";
    if( ptc == &convert_uint8_to_uint64 )
        return "convert_uint8_to_uint64";
    if( ptc == &convert_uint16_to_uint8 )
        return "convert_uint16_to_uint8";
    if( ptc == &convert_uint16_to_uint32 )
        return "convert_uint16_to_uint32";
    if( ptc == &convert_uint16_to_uint64 )
        return "convert_uint16_to_uint64";
    if( ptc == &convert_uint32_to_uint8 )
        return "convert_uint32_to_uint8";
    if( ptc == &convert_uint32_to_uint16 )
        return "convert_uint32_to_uint16";
    if( ptc == &convert_uint32_to_uint64 )
        return "convert_uint32_to_uint64";
    if( ptc == &convert_uint64_to_uint8 )
        return "convert_uint64_to_uint8";
    if( ptc == &convert_uint64_to_uint16 )
        return "convert_uint64_to_uint16";
    if( ptc == &convert_uint64_to_uint32 )
        return "convert_uint64_to_uint32";
    if( ptc == &convert_uint8_to_int16 )
        return "convert_uint8_to_int16";
    if( ptc == &convert_uint8_to_int32 )
        return "convert_uint8_to_int32";
    if( ptc == &convert_uint8_to_int64 )
        return "convert_uint8_to_int64";
    if( ptc == &convert_uint16_to_int32 )
        return "convert_uint16_to_int32";
    if( ptc == &convert_uint16_to_int64 )
        return "convert_uint16_to_int64";
    if( ptc == &convert_uint32_to_int64 )
        return "convert_uint32_to_int64";
    if( ptc == &copy_1byte )
        return "copy_1byte";
    if( ptc == &copy_2bytes )
        return "copy_2bytes";
    if( ptc == &copy_4bytes )
        return "copy_4bytes";
    if( ptc == &copy_8bytes )
        return "copy_8bytes";
    if( ptc == &copy_16bytes )
        return "copy_16bytes";
    if( ptc == &copy_channel_string )
        return "copy_channel_sring";
    return "unknown";
}

namespace detail {
// Test to confirm whether a particular name is valid for a channel.  Basically, it's limited to simple identifiers,
// with '.' allowed.
template <class CharType>
bool is_valid_channel_name( const std::basic_string<CharType>& name ) {
    if( name.empty() )
        return false;

    for( unsigned i = 0; i < name.size(); ++i ) {
        CharType c = name[i];
        if( i == 0 ) {
            if( !isalpha( c ) && c != '_' )
                return false;
        } else {
            if( !isalnum( c ) && c != '_' )
                return false;
        }
    }
    return true;
}
}; // namespace detail

// Test to confirm whether a particular name is valid for a channel.  Basically, it's limited to simple identifiers.
bool is_valid_channel_name( const std::string& name ) { return detail::is_valid_channel_name( name ); }

// Test to confirm whether a particular name is valid for a channel.  Basically, it's limited to simple identifiers.
bool is_valid_channel_name( const std::wstring& name ) { return detail::is_valid_channel_name( name ); }
} // namespace channels
} // namespace frantic
