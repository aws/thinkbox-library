// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/misc/exception_stream.hpp>
#include <frantic/prtfile/prt2_common.hpp>
#include <frantic/strings/utf8.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::strings;
using namespace frantic::prtfile;
using namespace frantic::channels;
using frantic::invalid_particle_file_exception;
using frantic::graphics::raw_byte_buffer;

namespace {
static const char* prt2_compression_strings[prt2_compression_count] = { "uncompressed", "zlib", "transpose-zlib",
                                                                        "transpose",    "lz4",  "transpose-lz4" };
}

const char* prtfile::get_prt2_compression_string( prt2_compression_t compressionScheme ) {
    if( compressionScheme >= 0 && compressionScheme < prt2_compression_count ) {
        return prt2_compression_strings[compressionScheme];
    } else {
        stringstream ss;
        ss << "Invalid PRT2 compression scheme enum value " << compressionScheme;
        throw invalid_argument( ss.str() );
    }
}

prt2_compression_t prtfile::get_prt2_compression_scheme_from_string( const frantic::tstring& s ) {
    string compressionScheme = frantic::strings::to_string( s );
    for( int i = 0; i < prt2_compression_count; ++i ) {
        if( compressionScheme == get_prt2_compression_string( static_cast<prt2_compression_t>( i ) ) ) {
            return static_cast<prt2_compression_t>( i );
        }
    }

    stringstream ss;
    ss << "Unrecognized PRT2 compression scheme \"" << compressionScheme << "\"";
    throw invalid_argument( ss.str() );
}

int prtfile::sniff_prt_spec_revision( const frantic::tstring& filename ) {
    ifstream fin( filename.c_str(), ios::in | ios::binary );
    boost::uint64_t magic;
    fin.read( reinterpret_cast<char*>( &magic ), 8u );
    if( !fin ) {
        return 0;
    } else if( magic == PRT2_MAGIC_NUMBER ) {
        boost::int32_t version;
        fin.read( reinterpret_cast<char*>( &version ), 4u );
        return fin ? version : 0;
    } else if( magic == PRT1_MAGIC_NUMBER ) {
        boost::int32_t version;
        // The version is at offset 44 in PRT1 (8 bytes magic number, 4 bytes header length, 32 bytes signature string,
        // followed by version number)
        fin.seekg( 44u );
        fin.read( reinterpret_cast<char*>( &version ), 4u );
        return fin ? version : 0;
    } else {
        return 0;
    }
}

void prtfile::transpose_bytes_forward( size_t particleSize, size_t particleCount, const char* __restrict input,
                                       char* __restrict output ) {
    // Do the blocking along the particleCount axis
    size_t countRemaining = particleCount;
    while( countRemaining >= 128 ) {
        // Use a blocking of 128, to try and fit reasonably in L1 cache which may typically be 32KB data.
        // For example, 100 byte particles * 128 is 12KB, double for input and output is 24KB, which should be ok.
        // Also, by hardcoding the blocking, the compiler might be able to try some trickier unrolling, maybe.
        for( size_t j = 0; j != 128; ++j ) {
            for( size_t i = 0; i != particleSize; ++i ) {
                output[j + i * particleCount] = input[j * particleSize + i];
            }
        }
        input += 128 * particleSize;
        output += 128;
        countRemaining -= 128;
    }
    for( size_t j = 0; j != countRemaining; ++j ) {
        for( size_t i = 0; i != particleSize; ++i ) {
            output[j + i * particleCount] = input[j * particleSize + i];
        }
    }
}

void prtfile::transpose_bytes_reverse( size_t particleSize, size_t particleCount, const char* __restrict input,
                                       char* __restrict output ) {
    // Do the blocking along the particleCount axis
    size_t countRemaining = particleCount;
    while( countRemaining >= 128 ) {
        // Use a blocking of 128, to try and fit reasonably in L1 cache which may typically be 32KB data.
        // For example, 100 byte particles * 128 is 12KB, double for input and output is 24KB, which should be ok.
        // Also, by hardcoding the blocking, the compiler might be able to try some trickier unrolling, maybe.
        for( size_t j = 0; j != 128; ++j ) {
            for( size_t i = 0; i != particleSize; ++i ) {
                output[j * particleSize + i] = input[j + i * particleCount];
            }
        }
        input += 128;
        output += 128 * particleSize;
        countRemaining -= 128;
    }
    for( size_t j = 0; j != countRemaining; ++j ) {
        for( size_t i = 0; i != particleSize; ++i ) {
            output[j * particleSize + i] = input[j + i * particleCount];
        }
    }
}

static void skip_whitespace( const char*& begin, const char* end ) {
    while( begin < end && isspace( *begin ) )
        ++begin;
}

/**
 * (Taken from libdynd's BSD-licensed parsing routines)
 *
 * Skips whitespace, then matches the provided literal string token. On success,
 * returns true and modifies `rbegin` to point after the token. If the token is a
 * single character, use the other `parse_token` function which accepts
 * a char.
 *
 * Example:
 *     // Match the token "while"
 *     if (parse_token(begin, end, "while")) {
 *         // Handle while statement
 *     } else {
 *         // No while token found
 *     }
 */
template <int N>
static bool parse_token( const char*& rbegin, const char* end, const char ( &token )[N] ) {
    const char* begin = rbegin;
    skip_whitespace( begin, end );
    if( N - 1 <= end - begin && memcmp( begin, token, N - 1 ) == 0 ) {
        rbegin = begin + N - 1;
        return true;
    } else {
        return false;
    }
}

/**
 * (Taken from libdynd's BSD-licensed parsing routines)
 *
 * Parses an unsigned integer.
 *
 * Example:
 *     // Match a two digit month
 *     const char *match_begin, *match_end;
 *     if (parse_unsigned_int(begin, end, match_begin, match_end) {
 *         // Convert to int, process
 *     } else {
 *         // Couldn't match unsigned integer
 *     }
 */
static bool parse_unsigned_int( const char*& rbegin, const char* end, const char*& out_strbegin,
                                const char*& out_strend ) {
    const char* begin = rbegin;
    skip_whitespace( begin, end );
    if( begin < end ) {
        if( '1' <= *begin && *begin <= '9' ) {
            ++begin;
            while( begin < end && ( '0' <= *begin && *begin <= '9' ) ) {
                ++begin;
            }
            out_strbegin = rbegin;
            out_strend = begin;
            rbegin = begin;
            return true;
        } else if( *begin == '0' ) {
            if( begin + 1 < end && ( '0' <= *( begin + 1 ) && *( begin + 1 ) <= '9' ) ) {
                // Don't match leading zeros
                return false;
            } else {
                out_strbegin = begin;
                out_strend = begin + 1;
                rbegin = begin + 1;
                return true;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
}

/**
 * Parses a type id from the raw UTF-8 string buffer. Returns true on success, false on failure.
 */
static bool parse_type_id( const char* begin, const char* end, boost::int32_t& out_arity,
                           frantic::channels::data_type_t& out_dtype ) {
    out_arity = 1;
    out_dtype = data_type_invalid;
    const char *tok_begin, *tok_end;
    // Parse the "<arity> *" part if it exists
    if( parse_unsigned_int( begin, end, tok_begin, tok_end ) ) {
        out_arity = atoi( string( tok_begin, tok_end ).c_str() );
        if( !parse_token( begin, end, "*" ) ) {
            return false;
        }
    }
    // Parse the "<dtype>" part (Try common dtypes first for speed)
    if( parse_token( begin, end, "float32" ) ) {
        out_dtype = data_type_float32;
    } else if( parse_token( begin, end, "float16" ) ) {
        out_dtype = data_type_float16;
    } else if( parse_token( begin, end, "int32" ) ) {
        out_dtype = data_type_int32;
    } else if( parse_token( begin, end, "int64" ) ) {
        out_dtype = data_type_int64;
    } else if( parse_token( begin, end, "float64" ) ) {
        out_dtype = data_type_float64;
    } else if( parse_token( begin, end, "int8" ) ) {
        out_dtype = data_type_int8;
    } else if( parse_token( begin, end, "int16" ) ) {
        out_dtype = data_type_int16;
    } else if( parse_token( begin, end, "uint8" ) ) {
        out_dtype = data_type_uint8;
    } else if( parse_token( begin, end, "uint16" ) ) {
        out_dtype = data_type_uint16;
    } else if( parse_token( begin, end, "uint32" ) ) {
        out_dtype = data_type_uint32;
    } else if( parse_token( begin, end, "uint64" ) ) {
        out_dtype = data_type_uint64;
    } else if( parse_token( begin, end, "string" ) ) {
        out_dtype = data_type_string;
    } else {
        return false;
    }

    // Make sure we consumed the whole string
    skip_whitespace( begin, end );
    if( begin == end ) {
        return true;
    } else {
        return false;
    }
}

boost::uint64_t serialize::read_varint( const char*& begin, const char* end, const frantic::tstring& streamName ) {
    boost::uint64_t result = 0;
    boost::uint32_t shift = 0;
    // A varint cannot be more than 10 bytes
    for( int i = 0; i < 10 && begin < end; ++i, shift += 7 ) {
        boost::uint64_t base128digit = *reinterpret_cast<const boost::uint8_t*>( begin );
        ++begin;
        result |= ( base128digit & 0x7f ) << shift;
        if( ( base128digit & 0x80 ) == 0 ) {
            return result;
        }
    }
    throw invalid_particle_file_exception()
        << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName ) << "\" is corrupted.";
}

boost::uint64_t serialize::read_varint( std::istream& is, const frantic::tstring& streamName ) {
    boost::uint64_t result = 0;
    boost::uint32_t shift = 0;
    // A varint cannot be more than 10 bytes
    for( int i = 0; i < 10 && is; ++i, shift += 7 ) {
        boost::uint64_t base128digit = static_cast<boost::uint8_t>( is.get() );
        result |= ( base128digit & 0x7f ) << shift;
        if( ( base128digit & 0x80 ) == 0 ) {
            return result;
        }
    }
    throw invalid_particle_file_exception()
        << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName ) << "\" is corrupted.";
}

static size_t encode_varint( unsigned char* tmp, boost::uint64_t value ) {
    size_t size = 0;
    while( ( value & ~0x7fULL ) != 0 ) {
        tmp[size++] = static_cast<unsigned char>( ( ( value | 0x80U ) & 0xFF ) );
        value >>= 7;
    }
    tmp[size++] = static_cast<unsigned char>( value );
    return size;
}

void serialize::write_varint( raw_byte_buffer& buf, boost::uint64_t value ) {
    // It will take at most 10 bytes
    unsigned char tmp[10];
    size_t size = encode_varint( tmp, value );
    // Append the constructed varint to the buffer
    memcpy( buf.add_element( size ), tmp, size );
}

void serialize::write_varint( std::ostream& os, boost::uint64_t value ) {
    // It will take at most 10 bytes
    unsigned char tmp[10];
    size_t size = encode_varint( tmp, value );
    // Write the constructed varint to the buffer
    os.write( reinterpret_cast<const char*>( tmp ), size );
}

frantic::tstring serialize::read_varstring( const char*& begin, const char* end, const frantic::tstring& streamName ) {
    boost::uint64_t len = serialize::read_varint( begin, end, streamName );
    if( ( len > boost::uint64_t( end - begin ) ) || !is_valid_utf8( begin, begin + len ) ) {
        throw invalid_particle_file_exception() << "prt2_file_reader: The input stream \""
                                                << frantic::strings::to_string( streamName ) << "\" is corrupted.";
    }
    const char* saved_begin = begin;
    begin += len;
#ifdef _WIN32
    // On Windows, have to convert to tstring
    return to_tstring( frantic::strings::wstring_from_utf8( std::string( saved_begin, begin ) ) );
#else
    // on non-Windows platforms, the string stays UTF-8
    return frantic::tstring( saved_begin, begin );
#endif
}

frantic::tstring serialize::read_varstring( std::istream& is, const frantic::tstring& streamName ) {
    boost::uint64_t len = serialize::read_varint( is, streamName );
    std::string s;
    s.resize( len );
    if( len > 0 ) {
        is.read( &s[0], len );
    }
    if( !is || !is_valid_utf8( s ) ) {
        throw invalid_particle_file_exception() << "prt2_file_reader: The input stream \""
                                                << frantic::strings::to_string( streamName ) << "\" is corrupted.";
    }
#ifdef _WIN32
    // On Windows, have to convert to tstring
    return to_tstring( frantic::strings::wstring_from_utf8( s ) );
#else
    // on non-Windows platforms, the string stays UTF-8
    return s;
#endif
}

std::string serialize::read_varstring_utf8( std::istream& is, const frantic::tstring& streamName ) {
    boost::uint64_t len = serialize::read_varint( is, streamName );
    std::string s;
    s.resize( len );
    if( len > 0 ) {
        is.read( &s[0], len );
    }
    if( !is || !is_valid_utf8( s ) ) {
        throw invalid_particle_file_exception() << "prt2_file_reader: The input stream \""
                                                << frantic::strings::to_string( streamName ) << "\" is corrupted.";
    }
    return s;
}

void serialize::write_varstring( raw_byte_buffer& buf, const frantic::tstring& value ) {
#ifdef _WIN32
    string s = strings::to_utf8( value );
#else
    // On non-Windows platforms, already using utf-8
    const string& s = value;
#endif
    write_varint( buf, s.size() );
    memcpy( buf.add_element( s.size() ), s.data(), s.size() );
}

void serialize::write_varstring( std::ostream& os, const frantic::tstring& value ) {
#ifdef _WIN32
    string s = strings::to_utf8( value );
#else
    // On non-Windows platforms, already using utf-8
    const string& s = value;
#endif
    write_varint( os, s.size() );
    if( !s.empty() ) {
        os.write( s.data(), s.size() );
    }
}

void serialize::write_varstring_utf8( std::ostream& os, const std::string& value ) {
    write_varint( os, value.size() );
    if( !value.empty() ) {
        os.write( value.data(), value.size() );
    }
}

void serialize::read_type_id( const char*& begin, const char* end, const frantic::tstring& streamName,
                              boost::int32_t& out_arity, frantic::channels::data_type_t& out_dtype ) {
    boost::uint64_t len = serialize::read_varint( begin, end, streamName );
    if( ( len > boost::uint64_t( end - begin ) ) || !is_valid_utf8( begin, begin + len ) ) {
        throw invalid_particle_file_exception() << "prt2_file_reader: The input stream \""
                                                << frantic::strings::to_string( streamName ) << "\" is corrupted.";
    }
    if( !parse_type_id( begin, begin + len, out_arity, out_dtype ) ) {
        // TODO: Could more gracefully handle unrecognized type IDs
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" contains an unrecognized type id.";
    }
    begin += len;
}

void serialize::write_type_id( raw_byte_buffer& buf, boost::int32_t arity, frantic::channels::data_type_t dtype ) {
    stringstream ss;
    if( arity != 1 ) {
        ss << arity << " * ";
    }
    switch( dtype ) {
    case data_type_int8:
        ss << "int8";
        break;
    case data_type_int16:
        ss << "int16";
        break;
    case data_type_int32:
        ss << "int32";
        break;
    case data_type_int64:
        ss << "int64";
        break;
    case data_type_uint8:
        ss << "uint8";
        break;
    case data_type_uint16:
        ss << "uint16";
        break;
    case data_type_uint32:
        ss << "uint32";
        break;
    case data_type_uint64:
        ss << "uint64";
        break;
    case data_type_float16:
        ss << "float16";
        break;
    case data_type_float32:
        ss << "float32";
        break;
    case data_type_float64:
        ss << "float64";
        break;
    case data_type_string:
        ss << "string";
        break;
    default:
        throw invalid_argument( "prt2_file_writer: Encountered an invalid data type when writing a PRT type ID" );
    }
    string s = ss.str();
    write_varint( buf, s.size() );
    memcpy( buf.add_element( s.size() ), s.data(), s.size() );
}

prt2_compression_t serialize::read_compression_scheme( std::istream& is, const frantic::tstring& streamName ) {
    string compressionScheme = read_varstring_utf8( is, streamName );
    for( int i = 0; i < prt2_compression_count; ++i ) {
        if( compressionScheme == get_prt2_compression_string( static_cast<prt2_compression_t>( i ) ) ) {
            return static_cast<prt2_compression_t>( i );
        }
    }

    throw invalid_particle_file_exception()
        << "prt2_file_reader: Unrecognized PRT2 compression scheme \"" << compressionScheme << "\"";
}

void serialize::write_compression_scheme( std::ostream& os, prt2_compression_t compressionScheme ) {
    write_varstring_utf8( os, get_prt2_compression_string( compressionScheme ) );
}
