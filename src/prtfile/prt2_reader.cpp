// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/exception_stream.hpp>
#include <frantic/prtfile/prt2_common.hpp>
#include <frantic/prtfile/prt2_reader.hpp>
#include <frantic/strings/utf8.hpp>

#if defined( LZ4_AVAILABLE )
#include <lz4.h>
#endif
#include <zlib.h>

using namespace std;
using namespace frantic;
using namespace frantic::strings;
using namespace frantic::channels;
using namespace frantic::prtfile;
using frantic::invalid_particle_file_exception;
using frantic::particles::particle_file_metadata;

bool prt2_reader::is_position_offset() const {
    // It's Position+Offset mode if any of the particle filechunks were 'PrtO'
    if( m_defaultPRTChunk.isPrtO )
        return true;
    std::map<frantic::tstring, particle_chunk_named>::const_iterator it, itend;
    for( it = m_namedPRTChunks.begin(), itend = m_namedPRTChunks.end(); it != itend; ++it ) {
        if( it->second.isPrtO )
            return true;
    }
    return false;
}

const particles::particle_file_metadata& prt2_reader::get_metadata() const { return m_metadata; }

const property_map& prt2_reader::get_general_metadata() const { return m_metadata.get_general_metadata(); }

const property_map& prt2_reader::get_channel_metadata( const frantic::tstring& channelName ) const {
    const property_map* result = m_metadata.get_channel_metadata( channelName );
    if( !result ) {
        throw std::runtime_error( "prt2_reader.get_channel_metadata: Metadata for "
                                  "channel \"" +
                                  to_string( channelName ) +
                                  "\" is NULL in "
                                  "stream \"" +
                                  to_string( get_stream_name() ) + "\"" );
    }
    return *result;
}

/** Reads the header for a FileChunk, returns whether it succeeded */
static bool read_prt2_chunk_info( istream& in, boost::uint32_t& chunkType, boost::uint64_t& chunkSize ) {
    in.read( reinterpret_cast<char*>( &chunkType ), 4u );
    in.read( reinterpret_cast<char*>( &chunkSize ), 8u );
    if( in )
        return true;
    else
        return false;
}

/** Reads the PRT2 header from the input stream */
static void read_prt2_header( istream& in, const frantic::tstring& streamName ) {
    boost::uint64_t magic = 0;
    in.read( reinterpret_cast<char*>( &magic ), 8u );
    if( magic != PRT2_MAGIC_NUMBER ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" did not begin with the .prt file magic number: " << PRT2_MAGIC_NUMBER << ".";
    }

    boost::uint32_t version = 0;
    in.read( reinterpret_cast<char*>( &version ), 4u );
    if( version != 3 ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" contained an unsupported PRT file version: " << version << ", only version 3 is supported.";
    }
}

/** Reads the PRT2 'Chan' chunk from the input stream, populating the channel_map */
static void read_prt2_channels( istream& in, const frantic::tstring& streamName, channel_map& cm ) {
    boost::uint32_t channelName;
    in.read( reinterpret_cast<char*>( &channelName ), 4u );
    if( !in || channelName != PRT2_FILECHUNK_NAME_Chan ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" did not begin with a 'Chan' chunk.";
    }
    boost::uint64_t chunkSize;
    in.read( reinterpret_cast<char*>( &chunkSize ), 8u );
    if( !in ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" is corrupted, failed to read 'Chan' chunk.";
    }
    if( chunkSize > 50000 ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" had an unexpectedly large 'Chan' chunk of size " << chunkSize << ".";
    }

    // Read the whole chunk into a buffer, then process from there
    vector<char> buffer( chunkSize );
    if( chunkSize > 0 ) {
        in.read( &buffer[0], chunkSize );
    }
    if( !in ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" is corrupted, failed to read 'Chan' chunk.";
    }
    const char* begin = ( buffer.size() > 0 ? &buffer[0] : 0 );
    const char* end = begin + chunkSize;

    size_t channelOffset = 0;
    boost::uint64_t channelCount = serialize::read_varint( begin, end, streamName );
    for( boost::uint64_t i = 0; i != channelCount; ++i ) {
        frantic::tstring channelName = serialize::read_varstring( begin, end, streamName );
        boost::int32_t channelArity;
        data_type_t channelDataType;
        serialize::read_type_id( begin, end, streamName, channelArity, channelDataType );
        if( channelDataType == data_type_string ) {
            // Can't do strings (at least for now...)
            throw invalid_particle_file_exception()
                << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
                << "\" specifies an invalid 'string' channel type.";
        }
        boost::uint64_t sizeBytes = serialize::read_varint( begin, end, streamName );
        boost::uint64_t expectedSizeBytes = channelArity * sizeof_channel_data_type( channelDataType );
        if( sizeBytes != expectedSizeBytes ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
                << "\" had an unexpectedly channel bytes size for channel \"" << to_string( channelName ) << "\", "
                << sizeBytes << " instead of " << expectedSizeBytes << ".";
        }
        cm.define_channel( channelName, channelArity, channelDataType, channelOffset );
        channelOffset += sizeBytes;
    }

    // Make sure the whole chunk was consumed
    if( begin != end ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" has a malformed 'Chan' chunk.";
    }

    cm.end_channel_definition( 1, true, false );
}

static void split_metadata_name( const tstring& metadataName, tstring& outChannelName, tstring& outPropertyName ) {
    outChannelName.clear();
    outPropertyName.clear();

    const tstring::size_type dotPos = metadataName.find( '.' );
    if( dotPos == tstring::npos ) {
        outPropertyName.assign( metadataName );
    } else {
        outChannelName.assign( metadataName, 0, dotPos );
        outPropertyName.assign( metadataName, dotPos + 1, tstring::npos );
    }
}

static property_map* try_get_property_map( particle_file_metadata& metadata, const tstring& channelName ) {
    if( channelName.empty() ) {
        return &metadata.get_general_metadata();
    } else {
        return metadata.get_channel_metadata( channelName );
    }
}

static void read_prt2_meta( const char* begin, const char* end, const tstring& streamName,
                            particle_file_metadata& outMetadata ) {
    // Get the type information
    tstring metadataName = serialize::read_varstring( begin, end, streamName );
    boost::int32_t arity;
    data_type_t dataType;
    serialize::read_type_id( begin, end, streamName, arity, dataType );

    FF_LOG( debug ) << _T("PRT2: 'Meta' chunk has name \"") << metadataName << _T("\" and type ")
                    << channel_data_type_str( arity, dataType ) << endl;

    if( dataType == data_type_string && arity != 1 ) {
        // Arrays of strings not presently supported
        return;
    }

    tstring channelName, propertyName;
    split_metadata_name( metadataName, channelName, propertyName );

    property_map* pm = try_get_property_map( outMetadata, channelName );
    if( !pm ) {
        return;
    }

    // Add the channel to the property map
    channel_map propMap;
    propMap.define_channel( propertyName, arity, dataType );
    propMap.union_channel_map( pm->get_channel_map() );
    propMap.end_channel_definition();

    // This copies existing items in the property_map if they are present in the new map.
    pm->set_channel_map_with_swap( propMap );

    if( dataType == data_type_string ) {
        pm->set_cvt<frantic::tstring>( propertyName, serialize::read_varstring( begin, end, streamName ) );
    } else {
        const channel& ch = pm->get_channel_map()[propertyName];
        char* pDest = ch.get_channel_data_pointer( const_cast<char*>( pm->get_raw_buffer() ) );

        if( ch.primitive_size() != size_t( end - begin ) ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
                << "\" has corrupted metadata \"" << to_string( metadataName ) << "\", "
                << "primitive size " << ch.primitive_size() << " does not match the data size " << ( end - begin )
                << ".";
        }

        memcpy( pDest, begin, ch.primitive_size() );
    }
}

static void read_prt2_meta( istream& in, boost::uint64_t chunkSize, const tstring& streamName,
                            particle_file_metadata& outMetadata ) {
    // Read the whole chunk into a buffer, then process from there
    vector<char> buffer( chunkSize );
    if( chunkSize > 0 ) {
        in.read( &buffer[0], chunkSize );
    }
    if( !in ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" is corrupted, failed to read 'Meta' chunk.";
    }
    const char* begin = ( buffer.size() > 0 ? &buffer[0] : 0 );
    const char* end = begin + chunkSize;

    read_prt2_meta( begin, end, streamName, outMetadata );
}

static void read_prt2_prtchunk_index( istream& in, boost::uint64_t chunkSize, const frantic::tstring& streamName,
                                      vector<prt2_reader::particle_chunk_index_info>& particleChunkIndex ) {
    // Read the whole chunk into a buffer, then process from there
    vector<char> buffer( chunkSize );
    if( chunkSize > 0 ) {
        in.read( &buffer[0], chunkSize );
    }
    if( !in ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" is corrupted, error reading 'PIdx' chunk.";
    }
    const char* begin = ( buffer.size() > 0 ? &buffer[0] : 0 );
    const char* end = begin + chunkSize;

    boost::uint64_t particleChunkCount = serialize::read_value<boost::uint64_t>( begin, end, streamName );
    particleChunkIndex.resize( particleChunkCount + 1 );
    particleChunkIndex.front().chunkSeek = 0;
    particleChunkIndex.front().particleIndex = 0;
    for( boost::uint64_t i = 0; i != particleChunkCount; ++i ) {
        prt2_reader::particle_chunk_index_info& pcii0 = particleChunkIndex[i];
        prt2_reader::particle_chunk_index_info& pcii1 = particleChunkIndex[i + 1];
        boost::uint64_t chunkSize = serialize::read_varint( begin, end, streamName );
        boost::uint64_t chunkParticleCount = serialize::read_varint( begin, end, streamName );
        pcii1.chunkSeek = pcii0.chunkSeek + chunkSize;
        pcii1.particleIndex = pcii0.particleIndex + chunkParticleCount;
    }

    if( begin != end ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( streamName )
            << "\" is corrupted, invalid 'PIdx' chunk.";
    }
}

void prt2_reader::read_file_structure() {
    if( m_file.get() == NULL || !( *m_file ) ) {
        throw std::runtime_error( "prt2_reader.read_file_structure: The input stream \"" + to_string( m_streamname ) +
                                  "\" failed to open for reading." );
    }

    m_defaultPRTChunk.clear();

    // Read the header (it only contains information to validate)
    read_prt2_header( *m_file, m_streamname );

    // The first chunk is required to be the 'Chan' chunk
    read_prt2_channels( *m_file, m_streamname, m_particleChannelMap );

    // Initialize empty channel metadata based on the channel map
    for( std::size_t i = 0, ie = m_particleChannelMap.channel_count(); i < ie; ++i ) {
        m_metadata.get_channel_metadata( m_particleChannelMap[i].name(), true );
    }

    // Now process the rest of the chunks, saving any metadata and
    // saving the location of the 'Part' chunk
    boost::uint32_t chunkType;
    boost::uint64_t chunkSize;
    while( read_prt2_chunk_info( *m_file, chunkType, chunkSize ) ) {
        string chunkName( reinterpret_cast<const char*>( &chunkType ),
                          reinterpret_cast<const char*>( &chunkType ) + 4 );
        FF_LOG( debug ) << _T("PRT2: Reading chunk of type '") << frantic::strings::to_tstring( chunkName )
                        << _T( "' with size " << chunkSize << " at seek position " ) << m_file->tellg() << endl;
        if( chunkType == PRT2_FILECHUNK_NAME_Part || chunkType == PRT2_FILECHUNK_NAME_PrtO ) {
            // The particles, this FileChunk contains simpler PRTChunks within it
            size_t chunkStart = m_file->tellg();
            frantic::tstring particleStreamName = serialize::read_varstring( *m_file, m_streamname );
            particle_chunk_named* pcn;
            if( particleStreamName.empty() ) {
                pcn = &m_defaultPRTChunk;
            } else {
                pcn = &m_namedPRTChunks[particleStreamName];
            }
            pcn->compressionScheme = serialize::read_compression_scheme( *m_file, m_streamname );
            m_file->read( reinterpret_cast<char*>( &pcn->particleCount ), 8u );
            m_file->read( reinterpret_cast<char*>( &pcn->particleChunkCount ), 8u );
            if( !( *m_file ) || pcn->particleCount < 0 || pcn->particleChunkCount < 0 ) {
                throw invalid_particle_file_exception()
                    << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
                    << "\" is corrupted, possibly due to incomplete write.";
            }
            FF_LOG( debug ) << _T("PRT2: 'Part' chunk has name \"") << particleStreamName << _T("\", particle count ")
                            << pcn->particleCount << _T(", and particle chunk count ") << pcn->particleChunkCount
                            << endl;
            pcn->chunkSeek = m_file->tellg();
            pcn->chunkSize = chunkSize - ( pcn->chunkSeek - chunkStart );
            if( chunkType == PRT2_FILECHUNK_NAME_PrtO ) {
                // The 'PrtO' chunk includes an offset to be added to every position
                pcn->isPrtO = true;
                pcn->posAccessor = m_particleChannelMap.get_accessor<graphics::vector3f>( _T("Position") );
            } else {
                pcn->isPrtO = false;
            }
            // Seek past all the particle data to finish processing all the chunks
            m_file->seekg( pcn->chunkSize, ios::cur );
        } else if( chunkType == PRT2_FILECHUNK_NAME_Meta ) {
            read_prt2_meta( *m_file, chunkSize, m_streamname, m_metadata );
        } else if( chunkType == PRT2_FILECHUNK_NAME_PIdx ) {
            // Particle chunk index
            size_t chunkStart = m_file->tellg();
            frantic::tstring particleStreamName = serialize::read_varstring( *m_file, m_streamname );
            FF_LOG( debug ) << _T("PRT2: 'PIdx' chunk has name \"") << particleStreamName << _T("\"") << endl;
            particle_chunk_named* pcn;
            if( particleStreamName.empty() ) {
                pcn = &m_defaultPRTChunk;
            } else {
                pcn = &m_namedPRTChunks[particleStreamName];
            }
            read_prt2_prtchunk_index( *m_file, chunkSize - ( (size_t)m_file->tellg() - chunkStart ), m_streamname,
                                      pcn->particleChunkIndex );
        } else {
            /* TODO: Reintroduce way to warn, maybe a boolean parameter?
            // Unrecognized chunk, log it and skip past
            string chunkName( reinterpret_cast<const char*>( &chunkType ), reinterpret_cast<const char*>( &chunkType ) +
            4 ); FF_LOG( warning ) << _T("PRT2: Unrecognized PRT File chunk \"") << to_tstring( chunkName ) << _T("\" in
            stream
            \"")
                              << m_streamname << _T("\"") << std::endl;
            */
            // Save the unrecognized chunk for the user of the prt2_reader to inspect
            m_unrecognizedFileChunks.push_back( unrecognized_chunk_index_info() );
            m_unrecognizedFileChunks.back().chunkType = chunkType;
            m_unrecognizedFileChunks.back().chunkSeek = m_file->tellg();
            m_unrecognizedFileChunks.back().chunkSize = chunkSize;
            m_file->seekg( chunkSize, ios::cur );
        }
    }

    if( m_defaultPRTChunk.particleChunkIndex.empty() ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
            << "\" is corrupted, contained no 'PIdx' chunk for the default 'Part'/'PrtO' particles filechunk.";
    }
    if( boost::uint64_t( m_defaultPRTChunk.particleChunkCount ) + 1 != m_defaultPRTChunk.particleChunkIndex.size() ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
            << "\" is corrupted, contained an index inconsistent with the particles.";
    }
    for( std::map<frantic::tstring, particle_chunk_named>::const_iterator i = m_namedPRTChunks.begin(),
                                                                          iend = m_namedPRTChunks.end();
         i != iend; ++i ) {
        if( i->second.particleChunkIndex.empty() ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
                << "\" is corrupted, contained no 'PIdx' chunk for the 'Part'/'PrtO' particles filechunk named \""
                << to_string( i->first ) << "\".";
        }
        if( ( boost::uint64_t( i->second.particleChunkCount ) + 1 != i->second.particleChunkIndex.size() ) ||
            ( i->second.particleCount != i->second.particleChunkIndex.back().particleIndex ) ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
                << "\" is corrupted, contained an inconsistent 'PIdx' chunk for the 'Part'/'PrtO' particles filechunk "
                   "named "
                   "\""
                << to_string( i->first ) << "\".";
        }
    }

    if( m_defaultPRTChunk.particleCount < 0 ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
            << "\" is corrupted, did not contain a default 'Part' particles chunk.";
    }

    // Clear the failure bit of the istream, caused by trying to read a chunk at the end of the file
    m_file->clear();

    // After reading the metadata, if this is a file-openable stream, close it to save on resources
    if( m_isFileStream ) {
        close_stream();
    }
}

void prt2_reader::open( const frantic::tstring& filename ) {
    m_file.reset( new fstream( filename.c_str(), ios::in | ios::binary ) );
    m_streamname = filename;
    m_isFileStream = true;
    read_file_structure();
}

void prt2_reader::open( std::istream* is, const frantic::tstring& filename ) {
    m_file.reset( is );
    m_streamname = filename;
    m_isFileStream = false;
    read_file_structure();
}

/*
// This would be a way to do it with C++11
void prt2_reader::open( std::unique_ptr<std::istream>&& is, const frantic::tstring& filename ) {
  m_file = std::move( is );
  m_streamname = filename;
  read_file_structure();
}
*/

void prt2_reader::close() {
    close_stream();
    m_streamname = _T("");
    m_particleChannelMap.reset();
    m_metadata.clear();
    m_defaultPRTChunk.clear();
    m_namedPRTChunks.clear();
}

void prt2_reader::close_stream() { m_file.reset(); }

void prt2_reader::auto_open_stream() {
    if( !m_file.get() ) {
        reopen_stream();
    }
}

void prt2_reader::reopen_stream() {
    if( !m_isFileStream || m_streamname.empty() ) {
        throw std::runtime_error(
            "prt2_reader::reopen_stream : cannot re-open, there was no file associated with the stream." );
    }

    m_file.reset( new fstream( m_streamname.c_str(), ios::in | ios::binary ) );

    // check that the header is still valid
    read_prt2_header( *m_file, m_streamname );

    // TODO: also verify that none of the file information has changed? (e.g. same PCM, same metadata)
}

namespace {
class not_matches_chunk_type {
    boost::uint32_t m_chunkType;

  public:
    not_matches_chunk_type( boost::uint32_t chunkType )
        : m_chunkType( chunkType ) {}

    bool operator()( const prt2_reader::unrecognized_chunk_index_info& ucii ) { return ucii.chunkType != m_chunkType; }
};
} // anonymous namespace

void prt2_reader::remove_unrecognized_filechunks( boost::uint32_t chunkType,
                                                  std::vector<unrecognized_chunk_index_info>& outChunkInfo ) {
    vector<unrecognized_chunk_index_info>::iterator i = stable_partition(
        m_unrecognizedFileChunks.begin(), m_unrecognizedFileChunks.end(), not_matches_chunk_type( chunkType ) );
    copy( i, m_unrecognizedFileChunks.end(), back_inserter<vector<unrecognized_chunk_index_info>>( outChunkInfo ) );
    m_unrecognizedFileChunks.erase( i, m_unrecognizedFileChunks.end() );
}

void prt2_reader::read_unrecognized_filechunk( const unrecognized_chunk_index_info& unrecognizedChunk,
                                               char* outBuffer ) {
    auto_open_stream();

    m_file->seekg( unrecognizedChunk.chunkSeek );
    m_file->read( outBuffer, unrecognizedChunk.chunkSize );
    if( !( *m_file ) ) {
        string chunkName( reinterpret_cast<const char*>( &unrecognizedChunk.chunkType ),
                          reinterpret_cast<const char*>( &unrecognizedChunk.chunkType ) + 4 );
        throw invalid_particle_file_exception() << "prt2_file_reader: Failed to read data for PRT File chunk \""
                                                << chunkName << "\" in stream \"" << to_string( m_streamname ) << "\"";
    }
}

static void zlib_uncompress( size_t bytesIn, size_t bytesOut, const char* input, char* output,
                             const tstring& streamName ) {
    uLongf rawBufferSize = static_cast<uLongf>( bytesOut );
    int result = uncompress( reinterpret_cast<Bytef*>( output ), &rawBufferSize,
                             reinterpret_cast<const Bytef*>( input ), static_cast<uLongf>( bytesIn ) );
    if( result != Z_OK ) {
        stringstream ss;
        ss << "prt2_file_reader: ZLib uncompression failure reading from input stream \"" << to_string( streamName )
           << "\"";
        switch( result ) {
        case Z_MEM_ERROR:
            ss << ", not enough memory";
            break;
        case Z_BUF_ERROR:
            ss << ", not enough room in the buffer";
            break;
        case Z_DATA_ERROR:
            ss << ", corrupted or incomplete input data";
            break;
        }
        throw runtime_error( ss.str() );
    }
    if( rawBufferSize != bytesOut ) {
        stringstream ss;
        ss << "prt2_file_reader: ZLib uncompression failure reading from input stream \"" << to_string( streamName )
           << "\", the uncompressed data was too small";
        throw runtime_error( ss.str() );
    }
}

#if defined( LZ4_AVAILABLE )
static void lz4_uncompress( size_t bytesIn, size_t bytesOut, const char* input, char* output,
                            const tstring& streamName ) {
    int decompressResult =
        LZ4_decompress_safe( input, output, static_cast<int>( bytesIn ), static_cast<int>( bytesOut ) );
    if( decompressResult != (int)bytesOut ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: Failed to read particle chunk from the input stream \""
            << frantic::strings::to_string( streamName ) << "\", LZ4 reported error: "
            << ( decompressResult < 0 ? "invalid data" : "wrong number of bytes decompressed" );
    }
}
#endif

void prt2_reader::read_particle_chunk( const particle_chunk_named& pcn, size_t i, char* rawBuffer,
                                       graphics::vector3f* outPositionOffset ) {
    if( pcn.particleChunkIndex.empty() ) {
        throw invalid_particle_file_exception()
            << "prt2_file_reader: Cannot read particle chunk " << i << " from the input stream \""
            << frantic::strings::to_string( m_streamname ) << "\" because the file had no index";
    }

    auto_open_stream();

    const particle_chunk_index_info* pcii = &pcn.particleChunkIndex[i];
    if( (boost::int64_t)m_file->tellg() == pcn.chunkSeek + pcii[0].chunkSeek ) {
        // If the file pointer is at the beginning of this chunk, because the last thing read was the previous
        // chunk, just do a little read instead of a seek to stream more efficiently.
        char throwaway[8];
        m_file->read( throwaway, 8 );
    } else {
        boost::int64_t seekPosition = pcn.chunkSeek + pcii[0].chunkSeek + 8;
        m_file->seekg( seekPosition, ios::beg );
    }
    boost::uint64_t chunkSize = pcii[1].chunkSeek - pcii[0].chunkSeek - 8;
    // Get the position offset (0 for 'Part' filechunk, stored in the file for 'PrtO' filechunk)
    graphics::vector3f positionOffset( 0.f, 0.f, 0.f );
    if( pcn.isPrtO ) {
        m_file->read( reinterpret_cast<char*>( &positionOffset[0] ), 12u );
        chunkSize -= 12u;
    }
    if( chunkSize == 0 ) {
        if( pcii[1].particleIndex == pcii[0].particleIndex ) {
            return;
        } else {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: Zero-sized data encountered in particle chunk " << i
                << " from the input stream \"" << frantic::strings::to_string( m_streamname ) << "\"";
        }
    }
    boost::int64_t particleCount = pcii[1].particleIndex - pcii[0].particleIndex;
    boost::int64_t expectedOutputByteSize = particleCount * m_particleChannelMap.structure_size();
    switch( pcn.compressionScheme ) {
    case prt2_compression_uncompressed: {
        // Uncompressed, so the bytes should be particle data as is
        if( expectedOutputByteSize != (boost::int64_t)chunkSize ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
                << "\" has an invalid PRT chunk at index " << i << ", has size " << chunkSize << " but expected size "
                << expectedOutputByteSize;
        }
        // Read the particle data directly into the output buffer
        if( *m_file ) {
            m_file->read( rawBuffer, expectedOutputByteSize );
        }
        if( !( *m_file ) ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: Failed to read particle chunk " << i << " from the input stream \""
                << frantic::strings::to_string( m_streamname ) << "\"";
        }
    } break;
    case prt2_compression_transpose: {
        // Uncompressed, so the bytes should be particle data as is
        if( expectedOutputByteSize != (boost::int64_t)chunkSize ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
                << "\" has an invalid PRT chunk at index " << i << ", has size " << chunkSize << " but expected size "
                << expectedOutputByteSize;
        }
        // Read the particle data directly into the output buffer
        vector<char> buffer( expectedOutputByteSize );
        if( *m_file && expectedOutputByteSize > 0 ) {
            m_file->read( &buffer[0], expectedOutputByteSize );
        }
        if( !( *m_file ) ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: Failed to read particle chunk " << i << " from the input stream \""
                << frantic::strings::to_string( m_streamname ) << "\"";
        }
        transpose_bytes_reverse( m_particleChannelMap.structure_size(), particleCount,
                                 ( buffer.size() > 0 ? &buffer[0] : 0 ), rawBuffer );
    } break;
    case prt2_compression_zlib:
    case prt2_compression_transpose_zlib: {
        // Read the particle data into a temporary buffer
        vector<char> buffer( chunkSize );
        if( *m_file ) {
            m_file->read( &buffer[0], chunkSize );
        }
        if( !( *m_file ) ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: Failed to read particle chunk " << i << " from the input stream \""
                << frantic::strings::to_string( m_streamname ) << "\"";
        }
        if( pcn.compressionScheme == prt2_compression_zlib ) {
            zlib_uncompress( chunkSize, expectedOutputByteSize, &buffer[0], rawBuffer, m_streamname );
        } else {
            // Decompress, then reverse transpose the bytes
            vector<char> buffer2( expectedOutputByteSize );
            zlib_uncompress( chunkSize, expectedOutputByteSize, &buffer[0], ( buffer2.size() > 0 ? &buffer2[0] : 0 ),
                             m_streamname );
            transpose_bytes_reverse( m_particleChannelMap.structure_size(), particleCount,
                                     ( buffer2.size() > 0 ? &buffer2[0] : 0 ), rawBuffer );
        }
    } break;
#if defined( LZ4_AVAILABLE )
    case prt2_compression_lz4:
    case prt2_compression_transpose_lz4: {
        // Read the particle data into a temporary buffer
        vector<char> buffer( chunkSize );
        if( *m_file ) {
            m_file->read( &buffer[0], chunkSize );
        }
        if( !( *m_file ) ) {
            throw invalid_particle_file_exception()
                << "prt2_file_reader: Failed to read particle chunk " << i << " from the input stream \""
                << frantic::strings::to_string( m_streamname ) << "\"";
        }
        if( pcn.compressionScheme == prt2_compression_lz4 ) {
            lz4_uncompress( chunkSize, expectedOutputByteSize, &buffer[0], rawBuffer, m_streamname );
        } else {
            // Decompress, then reverse transpose the bytes
            vector<char> buffer2( expectedOutputByteSize );
            lz4_uncompress( chunkSize, expectedOutputByteSize, &buffer[0], ( buffer2.size() > 0 ? &buffer2[0] : 0 ),
                            m_streamname );
            transpose_bytes_reverse( m_particleChannelMap.structure_size(), particleCount,
                                     ( buffer2.size() > 0 ? &buffer2[0] : 0 ), rawBuffer );
        }
    } break;
#endif
    default:
        throw invalid_particle_file_exception()
            << "prt2_file_reader: The input stream \"" << frantic::strings::to_string( m_streamname )
            << "\" has an unrecognized compression scheme " << pcn.compressionScheme << ".";
    }

    if( outPositionOffset ) {
        // The caller is responsible for the position offset
        *outPositionOffset = positionOffset;
    } else if( !positionOffset.is_zero() ) {
        // Apply the position offset to every particle in the chunk
        size_t structureSize = m_particleChannelMap.structure_size();
        for( intptr_t i = 0; i < particleCount; ++i, rawBuffer += structureSize ) {
            pcn.posAccessor( rawBuffer ) += positionOffset;
        }
    }
}

void prt2_reader::read_particle_chunk( const particle_chunk_named& pcn, size_t i, std::vector<char>& rawBuffer,
                                       graphics::vector3f* outPositionOffset ) {
    const particle_chunk_index_info* pcii = &pcn.particleChunkIndex[i];
    boost::int64_t chunkByteSize =
        ( pcii[1].particleIndex - pcii[0].particleIndex ) * m_particleChannelMap.structure_size();
    rawBuffer.resize( chunkByteSize );
    if( chunkByteSize > 0 ) {
        read_particle_chunk( pcn, i, &rawBuffer[0], outPositionOffset );
    } else if( outPositionOffset ) {
        outPositionOffset->set( 0.f );
    }
}

void prt2_reader::read_particle_chunk( const particle_chunk_named& pcn, size_t i,
                                       frantic::particles::particle_array& parray,
                                       graphics::vector3f* outPositionOffset ) {
    if( parray.get_channel_map() != m_particleChannelMap ) {
        parray.resize( 0 );
        parray.set_channel_map( m_particleChannelMap );
    }
    parray.resize( get_particle_chunk_particle_count( i ) );
    read_particle_chunk( pcn, i, parray.at( 0 ), outPositionOffset );
}

void prt2_reader::read_particles( const particle_chunk_named& pcn, frantic::particles::particle_array& parray ) {
    // Create a channel map using  default settings, so as to be aligned instead of packed like the file
    channel_map cm;
    frantic::tstring channelName;
    size_t channelArity;
    data_type_t channelDataType;
    for( size_t i = 0; i < get_channel_count(); ++i ) {
        m_particleChannelMap.get_channel_definition( i, channelName, channelDataType, channelArity );
        cm.define_channel( channelName, channelArity, channelDataType );
    }
    cm.end_channel_definition();

    parray.resize( 0 );
    parray.set_channel_map( cm );
    parray.resize( pcn.particleCount );

    channel_map_adaptor cma( cm, m_particleChannelMap );
    frantic::particles::particle_array chunkBuffer( m_particleChannelMap );
    size_t offset = 0;
    for( boost::int64_t i = 0; i < pcn.particleChunkCount; ++i ) {
        // Read the chunk
        boost::int64_t chunkSize = get_particle_chunk_particle_count( pcn, i );
        chunkBuffer.resize( chunkSize );
        read_particle_chunk( pcn, i, chunkBuffer.at( 0 ) );
        // Copy it into the output parray
        for( boost::int64_t j = 0; j < chunkSize; ++j ) {
            cma.copy_structure( parray.at( offset + j ), chunkBuffer.at( j ) );
        }
        offset += chunkSize;
    }
}
