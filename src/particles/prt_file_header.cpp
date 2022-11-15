// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include <errno.h>

#include <boost/predef.h>

#include <frantic/channels/channel_map_const_iterator.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/exception_stream.hpp>
#include <frantic/particles/prt_file_header.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/strings/utf8.hpp>

#define HYBRID

// We want to enable
#if defined( NDEBUG ) && defined( HYBRID )
#define POPPED_NDEBUG
#undef NDEBUG
#endif

// <cassert> and <boost/assert.hpp> are designed to not use include guards so that we can explicitly change NDEBUG to
// enable their behavior.
#include <boost/assert.hpp>

#ifdef POPPED_NDEBUG
#define NDEBUG
#endif

using namespace frantic::channels;

#include <zlib.h>

#if defined( __APPLE__ ) && defined( BOOST_COMP_CLANG ) && defined( BOOST_LIB_STD_CXX ) && ( BOOST_LIB_STD_CXX )
#define USE_BOOST_FILE_STREAM_WRAPPER
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#elif defined( __GNUC__ )
#include <ext/stdio_filebuf.h>
#endif

// This stream is used for debug logging sometimes.
// extern std::ofstream fout;

// Global variable for setting PRT compression. TODO: Make it not a global.
int g_prtZlibCompressionFactor = Z_DEFAULT_COMPRESSION;

// Isolate the implementation details of the .prt format from the code at large.
namespace {

// The PRT channel types were made differently than the current channel_maps.  I would like
// to migrate to a text-based description of the particle channel map embedded in the file, instead of using
// an integer enumeration.
enum prt_channel_type {
    prt_ct_int16,
    prt_ct_int32,
    prt_ct_int64,

    prt_ct_float16,
    prt_ct_float32,
    prt_ct_float64,

    prt_ct_uint16,
    prt_ct_uint32,
    prt_ct_uint64,

    prt_ct_int8,
    prt_ct_uint8,

    prt_ct_invalid,

    prt_ct_utf8string = -1
};

data_type_t prt_channel_type_to_data_type( prt_channel_type pct ) {
    switch( pct ) {
    case prt_ct_int8:
        return data_type_int8;
    case prt_ct_int16:
        return data_type_int16;
    case prt_ct_int32:
        return data_type_int32;
    case prt_ct_int64:
        return data_type_int64;

    case prt_ct_uint8:
        return data_type_uint8;
    case prt_ct_uint16:
        return data_type_uint16;
    case prt_ct_uint32:
        return data_type_uint32;
    case prt_ct_uint64:
        return data_type_uint64;

    case prt_ct_float16:
        return data_type_float16;
    case prt_ct_float32:
        return data_type_float32;
    case prt_ct_float64:
        return data_type_float64;
    default:
        return data_type_invalid;
    }
}

prt_channel_type data_type_to_prt_channel_type( data_type_t dt ) {
    switch( dt ) {
    case data_type_int8:
        return prt_ct_int8;
    case data_type_int16:
        return prt_ct_int16;
    case data_type_int32:
        return prt_ct_int32;
    case data_type_int64:
        return prt_ct_int64;

    case data_type_uint8:
        return prt_ct_uint8;
    case data_type_uint16:
        return prt_ct_uint16;
    case data_type_uint32:
        return prt_ct_uint32;
    case data_type_uint64:
        return prt_ct_uint64;

    case data_type_float16:
        return prt_ct_float16;
    case data_type_float32:
        return prt_ct_float32;
    case data_type_float64:
        return prt_ct_float64;
    default:
        throw std::runtime_error( "data_type_to_prt_channel_type: Unknown data type: " +
                                  boost::lexical_cast<std::string>( dt ) );
    }
}

// The types supported as metadata extend the base types with a UTF-8 encoded string.
data_type_t prt_metadata_type_to_data_type( prt_channel_type pct ) {
    if( pct == prt_ct_utf8string )
        return data_type_string;
    return prt_channel_type_to_data_type( pct );
}

// The types supported as metadata extend the base types with a UTF-8 encoded string.
prt_channel_type data_type_to_prt_metadata_type( data_type_t dt ) {
    if( dt == data_type_string )
        return prt_ct_utf8string;
    return data_type_to_prt_channel_type( dt );
}

// This is the old header, for reference only
struct prt_header_v1 {
    boost::int64_t magicNumber;
    boost::int32_t headerLength;
    char fmtIdentStr[32];
    boost::int32_t version;
    boost::int64_t particleCount;
};

struct prt_channel_header_v1 {
    char channelName[32];
    boost::int32_t channelType;
    boost::int32_t channelArity;
    boost::int32_t channelOffset;
};

// Returns the 8 byte magic number that indicates this file format
boost::int64_t prt_magic_number() {
    static const unsigned char magic[] = { 192, 'P', 'R', 'T', '\r', '\n', 26, '\n' };
    return *(boost::int64_t*)magic;
}

// Returns the human readable signature string to embed in the file
const char* prt_signature_string() { return "Extensible Particle Format"; }

boost::int32_t meta_chunk_type() {
    static const unsigned char type[] = { 'M', 'e', 't', 'a' };
    return *reinterpret_cast<const boost::int32_t*>( type );
}

boost::int32_t stop_chunk_type() {
    static const unsigned char type[] = { 'S', 't', 'o', 'p' };
    return *reinterpret_cast<const boost::int32_t*>( type );
}

std::string generate_error_msg( const std::string& context, const frantic::tstring& streamName ) {
    std::stringstream ss;
    ss << context << " Failure to write header for file \"" << frantic::strings::to_string( streamName ) << "\"\n";
    ss << "\tError number: " << errno << "\n";
#ifdef _WIN32
    ss << "\tOS Error number: " << _doserrno << "\n";
#endif
    ss << "\tError message: " << strerror( errno ) << std::endl;

    return ss.str();
}

int fseek64( FILE* stream, boost::int64_t offset, int origin ) {
#ifdef _WIN32
    return _fseeki64( stream, offset, origin );
#elif defined( __APPLE__ )
    return fseeko( stream, offset, origin );
#else
    return fseeko64( stream, offset, origin );
#endif
}

boost::int64_t ftell64( FILE* stream ) {
#ifdef _WIN32
    return _ftelli64( stream );
#elif defined( __APPLE__ )
    return ftello( stream );
#else
    return ftello64( stream );
#endif
}

} // anonymous namespace

namespace frantic {
namespace particles {

namespace {
void add_required_metadata( property_map& metadata ) {
    // We manually merge in a BoundBox property so that it is always available.
    channel_map newMap;
    newMap.define_channel<frantic::graphics::boundbox3f>( _T("BoundBox") );
    newMap.end_channel_definition();

    property_map boundboxProp( newMap );
    boundboxProp.set_cvt<frantic::graphics::boundbox3f>(
        _T("BoundBox"),
        frantic::graphics::boundbox3f( frantic::graphics::vector3f( std::numeric_limits<float>::quiet_NaN() ),
                                       frantic::graphics::vector3f( std::numeric_limits<float>::quiet_NaN() ) ) );

    metadata.merge_property_map( boundboxProp );
}
} // anonymous namespace

void prt_file_header::set_general_metadata( const property_map& metadata ) {
    property_map newProps = metadata;

    add_required_metadata( newProps );

    m_metadata.get_general_metadata().swap( newProps );
}

void prt_file_header::set_all_metadata( const particle_file_metadata& data ) {
    m_metadata.clear();

    set_general_metadata( data.get_general_metadata() );

    const std::map<frantic::tstring, property_map>& channelMetadata = m_metadata.get_all_channel_metadata();
    for( std::map<frantic::tstring, property_map>::const_iterator it = channelMetadata.begin(),
                                                                  itEnd = channelMetadata.end();
         it != itEnd; ) {
        set_channel_metadata( it->first, it->second );
    }
}

namespace {
void read_string( std::istream& in, char ( &outString )[32] ) {
    in.getline( outString, 32, '\0' );

    if( !in )
        throw std::istream::failure(
            "prt_file_header.read_header: The input stream had an invalid string longer than 32 characters." );
}

#ifdef _WIN32
// On Windows, we need to convert from UTF8 to UTF16 to tstring.
inline frantic::tstring from_utf8( const std::string& utf8String ) {
    return frantic::strings::to_tstring( frantic::strings::wstring_from_utf8( utf8String ) );
}

inline frantic::tstring from_utf8( const char* szUtf8String ) {
    return frantic::strings::to_tstring( frantic::strings::wstring_from_utf8( szUtf8String ) );
}
#else
// On non-Windows std::string is already UTF8 so do nothing.
inline const std::string& from_utf8( const std::string& utf8String ) { return utf8String; }

inline const char* from_utf8( const char* szUtf8String ) { // TODO: Might make sense to return std::string.
    return szUtf8String;
}
#endif
void read_meta_chunk( std::istream& in, std::size_t chunkLength, property_map& inoutGlobalProps,
                      std::map<frantic::tstring, property_map>& inoutChannelProps ) {
    char szChannelName[32], szValueName[32];
    boost::int32_t metadataType = prt_ct_invalid;

    boost::int32_t valueSize = static_cast<boost::int32_t>( chunkLength );

    in.getline( szChannelName, 32, '\0' );
    if( !in )
        throw std::istream::failure(
            "prt_file_header.read_header: The input stream had an invalid channel name longer than 32 characters." );

    valueSize -= static_cast<boost::int32_t>( in.gcount() );

    in.getline( szValueName, 32, '\0' );
    if( !in )
        throw std::istream::failure(
            "prt_file_header.read_header: The input stream had an invalid value name longer than 32 characters." );

    valueSize -= static_cast<boost::int32_t>( in.gcount() );

    in.read( reinterpret_cast<char*>( &metadataType ), 4u );

    valueSize -= static_cast<boost::int32_t>( in.gcount() ); // Should always be 4.

    if( valueSize < 0 )
        throw std::runtime_error( "prt_file_header.read_header: The 'Meta' chunk had invalid length" );

    frantic::tstring channelName = from_utf8( szChannelName );
    frantic::tstring valueName = from_utf8( szValueName );

    // Interpret non-standard metadata with propery name
    // "channel.property" as a per-channel property for the "channel"
    // channel named "property".
    // Such non-standard files were written by Sequoia.
    if( channelName.empty() ) {
        const std::size_t dotPos = valueName.find( '.' );
        if( dotPos != frantic::tstring::npos ) {
            channelName.assign( valueName, 0, dotPos );
            valueName.assign( valueName, dotPos + 1, frantic::tstring::npos );
        }
    }

    property_map* pProps = &inoutGlobalProps;
    if( !channelName.empty() )
        pProps = &inoutChannelProps[channelName];

    data_type_t type = prt_channel_type_to_data_type( static_cast<prt_channel_type>( metadataType ) );
    if( type == data_type_invalid ) {
        in.seekg( static_cast<std::size_t>( valueSize ), std::ios::cur );
    } else if( type == data_type_string ) {
        std::vector<char> tempStringBuffer( static_cast<std::size_t>( valueSize ), '\0' );

        in.read( &tempStringBuffer.front(), tempStringBuffer.size() );

        channel_map propMap;
        propMap.define_channel( valueName, 1, data_type_string );
        propMap.union_channel_map( pProps->get_channel_map() );
        propMap.end_channel_definition();

        pProps->set_channel_map_with_swap(
            propMap ); // This copies existing items in the property_map if they are present in the new map.
        pProps->set_cvt<frantic::tstring>( valueName, from_utf8( &tempStringBuffer.front() ) );
    } else {
        if( ( static_cast<std::size_t>( valueSize ) % sizeof_channel_data_type( type ) ) != 0 ) {
            std::stringstream ss;
            ss << "Chunk size " << valueSize << " not divisible by size of chunk data type.";
            throw std::runtime_error( ss.str() );
        }

        std::size_t arity = static_cast<std::size_t>( valueSize ) / sizeof_channel_data_type( type );

        channel_map propMap;
        propMap.define_channel( valueName, arity, type );
        propMap.union_channel_map( pProps->get_channel_map() );
        propMap.end_channel_definition();

        pProps->set_channel_map_with_swap( propMap );

        const channel& ch = pProps->get_channel_map()[valueName];
        char* pDest = ch.get_channel_data_pointer( const_cast<char*>( pProps->get_raw_buffer() ) );

        in.read( pDest, ch.primitive_size() );
    }
}
} // anonymous namespace

// Function to read in the PRT header
void prt_file_header::read_header( std::istream& in, const frantic::tstring& streamName ) {
    // The boolean operator for std::istream will return true if the stream position is set beyond the end of the file
    // Even if failbit or badbit are set.
    // Setting these in the exeption mask for the stream causes an exception to be thrown if they are ever set,
    // preventing an infinite loop from occuring since our error detection was conditional on the boolean operator
    // returning false.
    std::ios_base::iostate exceptionMask = in.exceptions();
    exceptionMask |= std::ios_base::badbit;
    exceptionMask |= std::ios_base::failbit;
    in.exceptions( exceptionMask );

    try {

        if( !in )
            throw std::runtime_error( "prt_file_header.read_header: The input stream \"" +
                                      frantic::strings::to_string( streamName ) + "\" failed to open for reading." );

        // Note where the header begins, in case its not the beginning of the file.
        std::streampos headerStart = in.tellg();

        // Read in the main header data.
        // NOTE: When we have the v2 header, we'll first have to read up to the size of the v1 header before finishing
        // off the rest of that header, in order to maintain compatibility with the v1 version
        prt_header_v1 header;

        memset( &header, 0, sizeof( prt_header_v1 ) );

        // This code is used in a cross-platform manner, and I'm concerned that expecting the layout of a C-style struct
        // to be consistent might cause trouble. So far, we have chosen data types such that there are no packing bytes
        // introduced but that is a sketchy proposition since alignment and padding bytes are compiler implementation
        // defined.
        in.read( reinterpret_cast<char*>( &header.magicNumber ), 8u );
        in.read( reinterpret_cast<char*>( &header.headerLength ), 4u );
        in.read( reinterpret_cast<char*>( &header.fmtIdentStr[0] ), 32u );
        in.read( reinterpret_cast<char*>( &header.version ), 4u );
        in.read( reinterpret_cast<char*>( &header.particleCount ), 8u );

        // Throw the invalid_particle_file_exception, which indicates that this is not a prt file (as opposed to a
        // corrupt prt file);
        if( header.magicNumber != prt_magic_number() )
            throw invalid_particle_file_exception()
                << "prt_file_header.read_header: The input stream \"" << frantic::strings::to_string( streamName )
                << "\" did not contain the .prt file magic number: " << prt_magic_number() << ".";

        // Throw the invalid_particle_file_exception, which indicates that this is not a prt file (as opposed to a
        // corrupt prt file);
        if( strncmp( prt_signature_string(), header.fmtIdentStr, 32 ) != 0 )
            throw invalid_particle_file_exception()
                << "prt_file_header.read_header: The input stream \"" << frantic::strings::to_string( streamName )
                << "\" did not contain the signature string '" << prt_signature_string() << "'.";

        m_particleCount = header.particleCount;
        if( m_particleCount < 0 && !m_allowNegativeCount )
            throw invalid_particle_file_exception()
                << "prt_file_header.read_header: The input stream \"" << frantic::strings::to_string( streamName )
                << "\" was not closed correctly and reported " << m_particleCount << " particles within.";

        if( header.version > 1 ) {
            boost::int32_t chunkType = 0;
            boost::int32_t chunkLength = 0;

            in.read( reinterpret_cast<char*>( &chunkType ), 4u );
            in.read( reinterpret_cast<char*>( &chunkLength ), 4u );

            if( chunkLength < 0 ) {
                throw std::runtime_error( "prt_file_header.read_header: The input stream \"" +
                                          frantic::strings::to_string( streamName ) +
                                          "\" contained a chunk with a negative size." );
            }

            while( chunkType != stop_chunk_type() ) {
                std::streampos chunkStart = in.tellg();

                if( chunkType == meta_chunk_type() ) {
                    read_meta_chunk( in, static_cast<std::size_t>( chunkLength ), m_metadata.get_general_metadata(),
                                     m_metadata.get_all_channel_metadata() );
                } else {
                    // This is an unknown type, so skip it using the chunk length.
                    in.seekg( chunkLength, std::ios::cur );
                }

                BOOST_ASSERT( ( in.tellg() - chunkStart ) == chunkLength );

                // Read the next chunk
                in.read( reinterpret_cast<char*>( &chunkType ), 4u );
                in.read( reinterpret_cast<char*>( &chunkLength ), 4u );

                if( chunkLength < 0 ) {
                    throw std::runtime_error( +"prt_file_header.read_header: The input stream \"" +
                                              frantic::strings::to_string( streamName ) +
                                              "\" contained a chunk with a negative size." );
                }

                if( in.eof() ) {
                    throw std::istream::failure( "prt_file_header.read_header: Unexpected end of file while reading "
                                                 "header chunks in the input stream \"" +
                                                 frantic::strings::to_string( streamName ) + "\"." );
                }
            }

            BOOST_ASSERT( chunkLength == 0 );
            BOOST_ASSERT( ( in.tellg() - headerStart ) == header.headerLength );
        } else {
            // Skip parts of the file header which may have been added since the first version of the .prt format
            // We just read 56 bytes (ie. the expected size of the v1 header)
            if( header.headerLength > 56 )
                in.seekg( header.headerLength - 56, std::ios::cur );
        }

        boost::int32_t attrLength;
        in.read( reinterpret_cast<char*>( &attrLength ), 4u );

        if( attrLength != 4 ) {
            throw std::runtime_error( "prt_file_header.read_header: The reserved int value is not set to 4." );
        }

        // Read in the particle channel map
        m_particleChannelMap.reset();

        boost::int32_t channelCount, perChannelLength;
        in.read( reinterpret_cast<char*>( &channelCount ), 4u );
        in.read( reinterpret_cast<char*>( &perChannelLength ), 4u );

        // Calculate the known size of the per-channel header.
        boost::int32_t expectedSize = 44;

        // If we have less that out version of the spec requires, this is an invalid PRT file.
        if( perChannelLength < expectedSize )
            throw invalid_particle_file_exception()
                << "prt_file_header.read_header: The input stream \"" << frantic::strings::to_string( streamName )
                << "\" had invalid per-channel length.";

        for( int i = 0; i < channelCount; ++i ) {
            prt_channel_header_v1 channel;
            // in.read(reinterpret_cast<char*>(&channel), sizeof(prt_channel_header_v1));
            in.read( reinterpret_cast<char*>( &channel.channelName[0] ), 32u );
            in.read( reinterpret_cast<char*>( &channel.channelType ), 4u );
            in.read( reinterpret_cast<char*>( &channel.channelArity ), 4u );
            in.read( reinterpret_cast<char*>( &channel.channelOffset ), 4u );

            // Make sure the channel name is null terminated
            channel.channelName[31] = '\0';

            data_type_t channelDataType = prt_channel_type_to_data_type( prt_channel_type( channel.channelType ) );
            if( channelDataType == data_type_invalid )
                throw std::runtime_error( std::string() +
                                          "prt_file_header.read_header: The data type specified in channel \"" +
                                          channel.channelName + "\" in the input stream \"" +
                                          frantic::strings::to_string( streamName ) + "\" is not valid." );

            frantic::tstring channelName = frantic::strings::to_tstring( channel.channelName );

            m_particleChannelMap.define_channel( channelName, channel.channelArity, channelDataType,
                                                 channel.channelOffset );

            const property_map* channelData = m_metadata.get_channel_metadata( channelName );
            if( channelData != NULL && channelData->has_property( _T("Interpretation") ) ) {
                prt::channel_interpretation::option interp = prt::get_channel_interpretation( *channelData );

                if( !prt::channel_interpretation::is_compatible( interp, channelDataType,
                                                                 static_cast<std::size_t>( channel.channelArity ) ) ) {
                    throw std::runtime_error(
                        std::string() + "prt_file_header.write_header: The channel \"" + channel.channelName +
                        "\" is tagged as \"" +
                        frantic::strings::to_string( prt::channel_interpretation::to_string( interp ) ) +
                        "\" which is not compatible with " +
                        frantic::strings::to_string( channel_data_type_str(
                            static_cast<std::size_t>( channel.channelArity ), channelDataType ) ) );
                }
            }

            if( perChannelLength > expectedSize )
                in.seekg( perChannelLength - expectedSize, std::ios::cur ); // Skip unknown parts of the channel header
        }

        m_particleChannelMap.end_channel_definition( 1, true, false );

        // Discard any metadata for channels we don't actually have.
        std::map<frantic::tstring, property_map>& channelMetadata = m_metadata.get_all_channel_metadata();
        for( std::map<frantic::tstring, property_map>::iterator it = channelMetadata.begin(),
                                                                itEnd = channelMetadata.end();
             it != itEnd; ) {
            if( m_particleChannelMap.has_channel( it->first ) ) {
                ++it;
            } else {
                FF_LOG( warning ) << _T("Discarding metadata for channel \"") << it->first << _T("\"") << std::endl;

                // erase() does not invalidate other iterators in a std::map
                channelMetadata.erase( it++ );
            }
        }
    } catch( std::ios::failure& e ) {
        throw std::runtime_error( "prt_file_header::read_header(): error reading input stream \"" +
                                  frantic::strings::to_string( streamName ) + "\": " + e.what() );
    }
}

namespace {
/**
 * Returns the size in bytes of a 'Meta' chunk's data, for the given value name & type. This does not include the chunk
 * type and length fields.
 */
std::size_t count_meta_chunk_bytes( const frantic::tstring& channelName, const channel& ch, const char* pData ) {
    std::string utf8ChannelName = frantic::strings::to_utf8( channelName );
    if( !utf8ChannelName.empty() && !is_valid_channel_name( utf8ChannelName ) )
        throw std::runtime_error( "prt_file_header.write_header: Invalid channel name \"" + utf8ChannelName +
                                  "\" in metadata chunk" );

    std::string utf8PropName = frantic::strings::to_utf8( ch.name() );
    if( !is_valid_channel_name( utf8PropName ) )
        throw std::runtime_error( "prt_file_header.write_header: Invalid metadata name \"" + utf8PropName + "\" in \"" +
                                  ( utf8ChannelName.empty() ? "Global" : utf8ChannelName.c_str() ) +
                                  "\" metadata chunk" );

    std::size_t result = 0u;
    result += utf8ChannelName.size() + 1u; // Channel name
    result += utf8PropName.size() + 1u;    // Value name
    result += 4u;                          // Value type

    if( ch.data_type() == data_type_string ) {
        std::string utfPropValue = frantic::strings::to_utf8( channels::cstring_from_channel_string( pData ) );

        result += utfPropValue.size() + 1u;
    } else {
        result += ch.primitive_size();
    }

    return result;
}

/**
 * Returns the size in bytes of the 'Meta' chunks for the given channel_map. The result includes the size of the chunks
 * type & length fields.
 */
std::size_t count_meta_chunk_bytes( const frantic::tstring& channelName, const property_map& props ) {
    std::size_t result = 0u;

    for( channel_map_const_iterator it = begin( props.get_channel_map() ), itEnd = end( props.get_channel_map() );
         it != itEnd; ++it ) {
        // Chunk type, length & data section
        result += 4u + 4u +
                  count_meta_chunk_bytes( channelName, *it, it->get_channel_data_pointer( props.get_raw_buffer() ) );
    }

    return result;
}

/**
 * Writes a series 'Meta' chunks to the stream.
 * \param out The stream to write to.
 * \param channelName The channel these 'Meta' chunks are associated with. Can be an
 * empty string for the global 'Meta' chunk.
 * \param props A 'Meta' chunk is emitted for each property in the property_map.
 * \return Returns the stream position of the global chunk's 'BoundBox' metadata value. The result
 * is meaningless for non-global chunks.
 */
std::streampos write_meta_chunks( std::ostream& out, const frantic::tstring& channelName, const property_map& props ) {
    std::streampos result = static_cast<std::streampos>( 0 );
    std::streampos chunkStart;

    std::string utf8ChannelName = frantic::strings::to_utf8( channelName );
    if( !utf8ChannelName.empty() && !is_valid_channel_name( utf8ChannelName ) )
        throw std::runtime_error( "prt_file_header.write_header: Invalid channel name \"" + utf8ChannelName +
                                  "\" in metadata chunk" );

    for( channel_map_const_iterator it = begin( props.get_channel_map() ), itEnd = end( props.get_channel_map() );
         it != itEnd; ++it ) {
        const char* pData = it->get_channel_data_pointer( props.get_raw_buffer() );

        boost::int32_t chunkType = meta_chunk_type();
        boost::int32_t chunkLength =
            static_cast<boost::int32_t>( count_meta_chunk_bytes( channelName, *it, pData ) ); // Inefficient.

        out.write( reinterpret_cast<const char*>( &chunkType ), 4u );
        out.write( reinterpret_cast<const char*>( &chunkLength ), 4u );

        std::string utf8PropName = frantic::strings::to_utf8( it->name() );
        if( !is_valid_channel_name( utf8PropName ) )
            throw std::runtime_error( "prt_file_header.write_header: Invalid metadata name \"" + utf8PropName +
                                      "\" in \"" + ( utf8ChannelName.empty() ? "Global" : utf8ChannelName.c_str() ) +
                                      "\" metadata chunk" );

        boost::int32_t type = static_cast<boost::int32_t>( data_type_to_prt_channel_type( it->data_type() ) );

        chunkStart = out.tellp();

        out.write( utf8ChannelName.c_str(), utf8ChannelName.size() + 1u );
        out.write( utf8PropName.c_str(), utf8PropName.size() + 1u );
        out.write( reinterpret_cast<const char*>( &type ), 4u );

        if( type == prt_ct_utf8string ) {
            std::string utfPropValue = frantic::strings::to_utf8( channels::cstring_from_channel_string( pData ) );

            out.write( utfPropValue.c_str(), utfPropValue.size() + 1u );
        } else {
            if( channelName.empty() && it->name() == _T("BoundBox") )
                result = out.tellp();

            out.write( pData, it->primitive_size() );
        }

        // Make sure our length calculation was consistent with reality.
        BOOST_ASSERT( ( out.tellp() - chunkStart ) == chunkLength );
    }

    return result;
}

void write_stop_chunk( std::ostream& out ) {
    boost::int32_t chunkType = stop_chunk_type();
    boost::int32_t chunkLength = 0;

    out.write( reinterpret_cast<const char*>( &chunkType ), 4u );
    out.write( reinterpret_cast<const char*>( &chunkLength ), 4u );
}
} // namespace

// Function to write out the PRT header.  The streamName parameter should be the filename for
// embedding in error messages.
void prt_file_header::write_header( std::ostream& out, const frantic::tstring& streamName ) {
    if( !out )
        throw std::runtime_error( "prt_file_header.write_header: The output stream \"" +
                                  frantic::strings::to_string( streamName ) + "\" failed to open for writing." );

    add_required_metadata( m_metadata.get_general_metadata() );

    std::streampos fileStart = out.tellp();

    m_boundboxPos = std::streampos( 0 );

    // Calculate the total size of the metadata section.
    std::streamsize totalChunkSize = 0u;

    // Add a 'Meta' chunks for the global metadata.
    totalChunkSize += count_meta_chunk_bytes( _T(""), m_metadata.get_general_metadata() );

    // Add 'Meta' chunks for each channel's metadata.
    const std::map<frantic::tstring, property_map>& channelMetadata = m_metadata.get_all_channel_metadata();
    for( std::map<frantic::tstring, property_map>::const_iterator it = channelMetadata.begin(),
                                                                  itEnd = channelMetadata.end();
         it != itEnd; ++it )
        totalChunkSize += count_meta_chunk_bytes( it->first, it->second );

    // Add the required 'Stop' chunk.
    totalChunkSize += 4u + 4u;

    // Write the main header data
    prt_header_v1 header;
    memset( &header, 0, sizeof( prt_header_v1 ) );

    header.magicNumber = prt_magic_number();
    header.headerLength =
        56 + static_cast<int>( totalChunkSize ); // 56 for v1 header, then the various chunks written after it.
    strncpy( header.fmtIdentStr, prt_signature_string(), 32 );
    header.version = 2;
    header.particleCount = m_particleCount;

    // This code is used in a cross-platform manner, and I'm concerned that expecting the layout of a C-style struct to
    // be consistent might cause trouble. So far, we have chosen data types such that there are no packing bytes
    // introduced but that is a sketchy proposition since alignment and padding bytes are compiler implementation
    // defined. out.write(reinterpret_cast<const char*>(&header), sizeof(prt_header_v1));
    out.write( reinterpret_cast<const char*>( &header.magicNumber ), 8u );
    out.write( reinterpret_cast<const char*>( &header.headerLength ), 4u );
    out.write( reinterpret_cast<const char*>( &header.fmtIdentStr[0] ), 32u );
    out.write( reinterpret_cast<const char*>( &header.version ), 4u );

    m_particleCountPos = out.tellp();
    out.write( reinterpret_cast<const char*>( &header.particleCount ), 8u );

    std::streampos chunkStart = out.tellp();

    m_boundboxPos = write_meta_chunks( out, _T(""), m_metadata.get_general_metadata() );

    for( std::map<frantic::tstring, property_map>::const_iterator it = channelMetadata.begin(),
                                                                  itEnd = channelMetadata.end();
         it != itEnd; ++it )
        write_meta_chunks( out, it->first, it->second );

    write_stop_chunk( out );

    // Make sure our length calculations were correct.
    BOOST_ASSERT( ( out.tellp() - chunkStart ) == totalChunkSize );
    BOOST_ASSERT( ( out.tellp() - fileStart ) == 56 + totalChunkSize );

    //	fout << "written main header" << std::endl;

    // Write the reserved bytes
    boost::int32_t reservedInt = 4;
    out.write( reinterpret_cast<const char*>( &reservedInt ), 4 );

    // fout << "written reserved header info" << std::endl;

    // Write the channel map

    prt_channel_header_v1 prtChannel;

    boost::int32_t channelCount = (boost::int32_t)m_particleChannelMap.channel_count();
    boost::int32_t channelHeaderItemSize = 44;

    out.write( reinterpret_cast<const char*>( &channelCount ), 4 );
    out.write( reinterpret_cast<const char*>( &channelHeaderItemSize ), 4 );
    // fout << "Channel Count: " << channelCount << std::endl;
    // fout << "Channel Header Size: " << channelHeaderItemSize << std::endl;

    for( int i = 0; i < channelCount; ++i ) {
        // fout << "Reset channel buffer" << std::endl;
        memset( &prtChannel, 0, sizeof( prt_channel_header_v1 ) );
        strncpy( prtChannel.channelName, frantic::strings::to_string( m_particleChannelMap[i].name() ).c_str(), 32 );
        // fout << i << ") " << "Channel name " << m_particleChannelMap[i].name() << " copied " << std::endl;
        prtChannel.channelArity = (boost::int32_t)m_particleChannelMap[i].arity();
        prtChannel.channelType = data_type_to_prt_channel_type( m_particleChannelMap[i].data_type() );
        prtChannel.channelOffset = (boost::int32_t)m_particleChannelMap[i].offset();

        // fout << "Checking channel type " << m_particleChannelMap[i].type_str() << std::endl;
        if( prtChannel.channelType == prt_ct_invalid )
            throw std::runtime_error(
                std::string() + "prt_file_header.write_header: The input particle channel map contained a data type " +
                frantic::strings::to_string( channel_data_type_str( m_particleChannelMap[i].data_type() ) ) +
                ", which isn't supported by the .prt format currently." );

        // fout << m_particleChannelMap[i].name() << " channel type OK" << std::endl;

        // out.write(reinterpret_cast<const char*>(&prtChannel), sizeof(prt_channel_header_v1));
        out.write( reinterpret_cast<const char*>( &prtChannel.channelName[0] ), 32u );
        out.write( reinterpret_cast<const char*>( &prtChannel.channelType ), 4u );
        out.write( reinterpret_cast<const char*>( &prtChannel.channelArity ), 4u );
        out.write( reinterpret_cast<const char*>( &prtChannel.channelOffset ), 4u );
        // fout << "Channel written to stream" << std::endl;
    }
    // fout << "written channel map header info" << std::endl;

    if( !out )
        throw std::runtime_error( "prt_file_header.write_header: Failed to write the prt file header to stream \"" +
                                  frantic::strings::to_string( streamName ) + "\"." );
}

// Function to rewrite the particle count once the whole PRT file is written.
// The name parameter should be the filename for error messages.
void prt_file_header::rewrite_particle_count( std::ostream& out, const frantic::tstring& streamName,
                                              boost::int64_t particleCount ) {
    m_particleCount = particleCount;

    out.seekp( m_particleCountPos );

    if( !out )
        throw std::runtime_error(
            "prt_file_header.rewrite_particle_count: Failed to seek in output particle stream \"" +
            frantic::strings::to_string( streamName ) + "\" to save the particle count in the header." );

    // Write the new value
    out.write( reinterpret_cast<const char*>( &particleCount ), sizeof( boost::int64_t ) );

    if( !out )
        throw std::runtime_error( "prt_file_header.rewrite_particle_count: Failed to write the particle count to the "
                                  "header of the stream \"" +
                                  frantic::strings::to_string( streamName ) + "\"." );
}

void prt_file_header::rewrite_particle_boundbox( std::ostream& out, const frantic::tstring& streamName,
                                                 const frantic::graphics::boundbox3f& bounds ) {
    // Only write the boundbox if it has all non-NaN values.
    if( bounds.minimum().is_nan() || bounds.maximum().is_nan() )
        return;

    out.seekp( m_boundboxPos );

    if( !out )
        throw std::runtime_error(
            "prt_file_header.rewrite_particle_count: Failed to seek in output particle stream \"" +
            frantic::strings::to_string( streamName ) + "\" to save the boundbox in the header." );

    float tempBounds[] = { bounds.minimum().x, bounds.minimum().y, bounds.minimum().z,
                           bounds.maximum().x, bounds.maximum().y, bounds.maximum().z };

    // Write the new value
    out.write( reinterpret_cast<const char*>( &tempBounds[0] ), sizeof( float ) * 6u );

    if( !out )
        throw std::runtime_error(
            "prt_file_header.rewrite_particle_count: Failed to write the boundbox to the header of the stream \"" +
            frantic::strings::to_string( streamName ) + "\"." );
}

namespace {
class auto_istream : public std::istream {
    std::unique_ptr<std::streambuf> m_streamBufImpl;

  public:
    explicit auto_istream( std::streambuf* pImpl )
        : std::istream( pImpl )
        , m_streamBufImpl( pImpl ) {}

    static std::streambuf* create_wrapper_streambuf( FILE* f );
};

class auto_ostream : public std::ostream {
    std::unique_ptr<std::streambuf> m_streamBufImpl;

  public:
    explicit auto_ostream( std::streambuf* pImpl )
        : m_streamBufImpl( pImpl )
        , std::ostream( pImpl ) {}

    static std::streambuf* create_wrapper_streambuf( FILE* f );
};

std::streambuf* auto_istream::create_wrapper_streambuf( FILE* f ) {
#if defined( USE_BOOST_FILE_STREAM_WRAPPER )
    using namespace boost::iostreams;
    std::unique_ptr<stream_buffer<file_descriptor>> streamBuffer(
        new stream_buffer<file_descriptor>( fileno( f ), never_close_handle ) );
    return streamBuffer.release();
#elif defined( _WIN32 )
    // This is a M$ extension allowing a C FILE* to be used in a std::istream. It is not automatically closed.
    return new std::filebuf( f );
#elif defined( __GNUC__ )
    // This is a GNU extension allowing a C FILE* to be used in a std::istream. See
    // http://gcc.gnu.org/onlinedocs/gcc-4.6.3/libstdc++/api/a00069.html
    return new __gnu_cxx::stdio_filebuf<char>( f, std::ios::in | std::ios::binary );
#else
#error std::streambuf* auto_istream::create_wrapper_streambuf( FILE* f ) not supported on this platform.
#endif
}

std::streambuf* auto_ostream::create_wrapper_streambuf( FILE* f ) {
#if defined( USE_BOOST_FILE_STREAM_WRAPPER )
    using namespace boost::iostreams;
    std::unique_ptr<stream_buffer<file_descriptor>> streamBuffer(
        new stream_buffer<file_descriptor>( fileno( f ), never_close_handle ) );
    return streamBuffer.release();
#elif defined( _WIN32 )
    // This is a M$ extension allowing a C FILE* to be used in a std::istream. It is not automatically closed.
    return new std::filebuf( f );
#elif defined( __GNUC__ )
    // This is a GNU extension allowing a C FILE* to be used in a std::istream. See
    // http://gcc.gnu.org/onlinedocs/gcc-4.6.3/libstdc++/api/a00069.html
    return new __gnu_cxx::stdio_filebuf<char>( f, std::ios::out | std::ios::binary );
#else
#error std::streambuf* auto_ostream::create_wrapper_streambuf( FILE* f ) not supported on this platform.
#endif
}
} // namespace

void prt_file_header::read_header( FILE* fin, const frantic::tstring& streamName ) {
    if( !fin )
        throw std::runtime_error( "prt_file_header.read_header: The input stream \"" +
                                  frantic::strings::to_string( streamName ) + "\" failed to open for reading." );

    auto_istream fileWrapper( auto_istream::create_wrapper_streambuf( fin ) );

#if defined( USE_BOOST_FILE_STREAM_WRAPPER )
    const boost::int64_t offset = ftell64( fin );
    fileWrapper.seekg( offset );
#endif

    this->read_header( fileWrapper, streamName );

#if defined( USE_BOOST_FILE_STREAM_WRAPPER )
    // The Boost stream_buffer reads extra data from fin to fill its internal buffer.
    // Reset fin to where we actually finished reading.
    fseek64( fin, fileWrapper.tellg(), 0 );
#endif
}

// Function to write out the PRT header.  The streamName parameter should be the filename for
// embedding in error messages.
void prt_file_header::write_header( FILE* out, const frantic::tstring& streamName ) {
    if( !out )
        throw std::runtime_error( "prt_file_header.write_header: The output stream \"" +
                                  frantic::strings::to_string( streamName ) + "\" failed to open for writing." );

    auto_ostream fileWrapper( auto_ostream::create_wrapper_streambuf( out ) );

#if defined( USE_BOOST_FILE_STREAM_WRAPPER )
    const boost::int64_t offset = ftell64( out );
    fileWrapper.seekp( offset );
#endif

    this->write_header( fileWrapper, streamName );

#if defined( USE_BOOST_FILE_STREAM_WRAPPER )
    // The Boost stream_buffer reads extra data from fin to fill its internal buffer.
    // Reset fin to where we actually finished reading.
    fseek64( out, fileWrapper.tellp(), 0 );
#endif
}

// Function to rewrite the particle count once the whole PRT file is written.
// The name parameter should be the filename for error messages.
void prt_file_header::rewrite_particle_count( FILE* out, const frantic::tstring& streamName,
                                              boost::int64_t particleCount ) {
    auto_ostream fileWrapper( auto_ostream::create_wrapper_streambuf( out ) );

    this->rewrite_particle_count( fileWrapper, streamName, particleCount );
}

void prt_file_header::rewrite_particle_boundbox( FILE* out, const frantic::tstring& streamName,
                                                 const frantic::graphics::boundbox3f& bounds ) {
    auto_ostream fileWrapper( auto_ostream::create_wrapper_streambuf( out ) );

    this->rewrite_particle_boundbox( fileWrapper, streamName, bounds );
}

} // namespace particles
} // namespace frantic
