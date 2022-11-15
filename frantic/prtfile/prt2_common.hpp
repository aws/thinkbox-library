// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/misc/exception_stream.hpp>

// NOTE: This here is assuming a little-endian machine, as do the PRT and PRT2 I/O routines.

// {0xC0, 'P', 'R', 'T', '\r', '\n', 0x1A, '\n'}
#define PRT1_MAGIC_NUMBER ( 0x0A1A0A0D545250C0ULL )

// {0xC0, 'P', 'R', 'T', '2', '\r', '\n', 0x1A}
#define PRT2_MAGIC_NUMBER ( 0x1A0A0D32545250C0ULL )

/** FileChunk specifying channels */
#define PRT2_FILECHUNK_NAME_Chan ( 0x6e616843U )
/** FileChunk containing particles */
#define PRT2_FILECHUNK_NAME_Part ( 0x74726150U )
/** FileChunk containing particles with a per-chunk Position offset */
#define PRT2_FILECHUNK_NAME_PrtO ( 0x4f747250U )
/** FileChunk with a bit of metadata */
#define PRT2_FILECHUNK_NAME_Meta ( 0x6174654DU )
/** FileChunk with an index for the 'Part' chunk(s) */
#define PRT2_FILECHUNK_NAME_PIdx ( 0x78644950U )
/** For SPRT files, FileChunk with an octree index for the 'Part' chunk(s) */
#define PRT2_FILECHUNK_NAME_SIdx ( 0x78644953U )

namespace frantic {
namespace prtfile {

enum prt2_compression_t {
    prt2_compression_uncompressed,
    prt2_compression_zlib,
    prt2_compression_transpose_zlib,
    prt2_compression_transpose, // No compression, just a temporary hack to experiment with
    prt2_compression_lz4,
    prt2_compression_transpose_lz4,

    // Number of compression modes
    prt2_compression_count,

    // Value for default parameters to use, so easily switchable here
    prt2_compression_default = prt2_compression_transpose_zlib,
};

const char* get_prt2_compression_string( prt2_compression_t compressionScheme );

prt2_compression_t get_prt2_compression_scheme_from_string( const frantic::tstring& s );

/**
 * Peeks at the header of the file, returns 0 if it's neither PRT1 nor PRT2, and the specification revision number
 * otherwise. This is 1 for PRT1.0, 2 for PRT1.1, and 3 for PRT2.0.
 */
int sniff_prt_spec_revision( const frantic::tstring& filename );

/** Transposition for putting all the channels together/splitting them back apart */
void transpose_bytes_forward( size_t particleSize, size_t particleCount, const char* __restrict input,
                              char* __restrict output );
void transpose_bytes_reverse( size_t particleSize, size_t particleCount, const char* __restrict input,
                              char* __restrict output );

namespace serialize {
/**
 * Reads a POD value of type T from the buffer, incrementing 'begin' as necessary.
 *
 * Throws an exception if there's not enough data.
 */
template <typename T>
inline T read_value( const char*& begin, const char* end, const frantic::tstring& streamName ) {
    if( sizeof( T ) <= end - begin ) {
        T result;
        memcpy( &result, begin, sizeof( T ) );
        begin += sizeof( T );
        return result;
    } else {
        throw invalid_particle_file_exception() << "prt2_file_reader: The input stream \""
                                                << frantic::strings::to_string( streamName ) << "\" is corrupted.";
    }
}

/** Writes a POD value of type T to the buffer. */
template <typename T>
inline void write_value( graphics::raw_byte_buffer& buf, const T& value ) {
    memcpy( buf.add_element( sizeof( T ) ), reinterpret_cast<const char*>( &value ), sizeof( T ) );
}

/**
 * Reads a uint64_t stored as varint from the buffer, incrementing 'begin' as necessary.
 *
 * Throws an exception if there's not enough data or the data is invalid.
 */
boost::uint64_t read_varint( const char*& begin, const char* end, const frantic::tstring& streamName );

/** Reads a varint from a std::istream */
boost::uint64_t read_varint( std::istream& is, const frantic::tstring& streamName );

/** Writes a uint64_t as a varint into the buffer. */
void write_varint( graphics::raw_byte_buffer& buf, boost::uint64_t value );

/** Writes a uint64_t as a varint into the std::ostream. */
void write_varint( std::ostream& os, boost::uint64_t value );

/**
 * Extracts a string stored as varstring from the buffer, incrementing 'begin' as necessary.
 *
 * Throws an exception if there's not enough data or the data is invalid.
 */
frantic::tstring read_varstring( const char*& begin, const char* end, const frantic::tstring& streamName );

/** Reads a varstring from a std::istream */
frantic::tstring read_varstring( std::istream& is, const frantic::tstring& streamName );

/** Reads a varstring from a std::istream as utf8 */
std::string read_varstring_utf8( std::istream& is, const frantic::tstring& streamName );

/** Writes a string as a varstring into the buffer. */
void write_varstring( graphics::raw_byte_buffer& buf, const frantic::tstring& value );

/** Writes a string as a varstring into the std::ostream. */
void write_varstring( std::ostream& os, const frantic::tstring& value );

/** Writes a utf8 string as a varstring into the std::ostream. */
void write_varstring_utf8( std::ostream& os, const std::string& value );

/**
 * Read the type ID, stored as a varstring, from the binary buffer. incrementing 'begin' as necessary. Throws an
 * exception if there's a problem.
 */
void read_type_id( const char*& begin, const char* end, const frantic::tstring& streamName, boost::int32_t& out_arity,
                   frantic::channels::data_type_t& out_dtype );

/** Writes a type ID as a varstring into the buffer. */
void write_type_id( graphics::raw_byte_buffer& buf, boost::int32_t arity, frantic::channels::data_type_t dtype );

/** Reads the compression scheme, stored as a varstring */
prt2_compression_t read_compression_scheme( std::istream& is, const frantic::tstring& streamName );

/** Writes the compression scheme as a varstring */
void write_compression_scheme( std::ostream& os, prt2_compression_t compressionScheme );

namespace detail {
// A helper class for writing metadata, use specialization to distinguish the string case
template <class T>
struct metadata_writer {
    static void write( frantic::graphics::raw_byte_buffer& buf, const frantic::tstring& metaName, const T& value ) {
        write_varstring( buf, metaName );
        write_type_id( buf, (boost::int32_t)frantic::channels::channel_data_type_traits<T>::arity(),
                       frantic::channels::channel_data_type_traits<T>::data_type() );
        write_value( buf, value );
    }
};

template <>
struct metadata_writer<frantic::tstring> {
    static void write( frantic::graphics::raw_byte_buffer& buf, const frantic::tstring& metaName,
                       const frantic::tstring& value ) {
        write_varstring( buf, metaName );
        write_type_id( buf, 1, frantic::channels::data_type_string );
        write_varstring( buf, value );
    }
};
} // namespace detail

/** Writes the bytes needed for a 'GMet' chunk of the specified metadata name and value into the buffer */
template <class T>
void write_metadata_filechunk( frantic::graphics::raw_byte_buffer& buf, const frantic::tstring& metaName,
                               const T& value ) {
    detail::metadata_writer<T>::write( buf, metaName, value );
}

} // namespace serialize

} // namespace prtfile
} // namespace frantic
