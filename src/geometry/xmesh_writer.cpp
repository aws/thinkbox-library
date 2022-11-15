// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/files/compression_stream.hpp>
#include <frantic/files/files.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/geometry/xmesh_writer.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/strings/to_string_classic.hpp>
#include <frantic/strings/utf8.hpp>

#include <tinyxml2.h>

#include <boost/filesystem/path.hpp>
#include <boost/numeric/conversion/cast.hpp>


namespace frantic {
namespace geometry {

namespace {

// TODO: It's not obvious how char string channel names will be encoded.
// I think it should be UTF-8, but we have code that can produce different
// encodings.  For now, I assume that it's UTF-8 if it's valid UTF-8.

inline std::wstring wstring_from_channel_name( const std::string& s ) {
    if( frantic::strings::is_valid_utf8( s ) ) {
        return frantic::strings::wstring_from_utf8( s );
    } else {
        return frantic::strings::wstring_from_utf8( frantic::strings::to_utf8( s ) );
    }
}

inline const std::wstring& wstring_from_channel_name( const std::wstring& s ) { return s; }

inline std::string utf8_from_channel_name( const std::string& s ) {
    if( frantic::strings::is_valid_utf8( s ) ) {
        return s;
    } else {
        return frantic::strings::to_utf8( s );
    }
}

inline std::string utf8_from_channel_name( const std::wstring& s ) { return frantic::strings::to_utf8( s ); }

} // anonymous namespace

// Saves a single array of data.
// Format:
//   "xmeshdat"     8 bytes
//   "<datatag>"    12 bytes (based on the named channel data type, so all 12 bytes may not actually be used)
//   count          4 bytes
//   dataSize       4 bytes
//   <array data>   (count * dataSize) bytes
void write_xmesh_array_file( const boost::filesystem::path& path, const void* data, const std::string& dataType,
                             std::size_t dataSize, std::size_t numData, int compressionLevel ) {
    frantic::files::file_ptr fout( frantic::files::fopen( path, "wb" ) );

    if( !fout.get() )
        throw std::runtime_error( "write_xmesh_array_file() - Failed to open file \"" + path.string() +
                                  "\" for output." );

    frantic::files::zlib_gzip_ostream_cstdio zlibOut;
    zlibOut.open( fout, frantic::files::to_tstring( path ), compressionLevel );

    const char* header = "xmeshdat"; // header string

    // the data type is variable length but we'll put it in a 12byte char array and zero what isnt used
    char charDataType[12];
    memset( charDataType, 0, 12 );
    memcpy( charDataType, dataType.c_str(), dataType.length() );

    typedef boost::int32_t num_data_t;
    if( numData > (size_t)std::numeric_limits<num_data_t>::max() ) {
        throw std::runtime_error(
            "write_xmesh_array_file() - the data count " + frantic::strings::to_string_classic<std::string>( numData ) +
            " exceeds the maximum data count (" +
            frantic::strings::to_string_classic<std::string>( std::numeric_limits<num_data_t>::max() ) + ")." );
    }
    const num_data_t numData32 = boost::numeric_cast<num_data_t>( numData );

    typedef boost::int32_t data_size_t;
    if( dataSize > (size_t)std::numeric_limits<data_size_t>::max() ) {
        throw std::runtime_error(
            "write_xmesh_array_file() - the data size " + frantic::strings::to_string_classic<std::string>( dataSize ) +
            " exceeds the maximum data size (" +
            frantic::strings::to_string_classic<std::string>( std::numeric_limits<data_size_t>::max() ) + ")." );
    }
    const data_size_t dataSize32 = boost::numeric_cast<data_size_t>( dataSize );

    // write the header data
    zlibOut.write( header, 8 );                                       // file tag
    zlibOut.write( charDataType, 12 );                                // write data type
    zlibOut.write( reinterpret_cast<const char*>( &numData32 ), 4 );  // write data amount
    zlibOut.write( reinterpret_cast<const char*>( &dataSize32 ), 4 ); // write data size

    // Write the data
    if( numData > 0 ) {
        zlibOut.write( reinterpret_cast<const char*>( data ), numData * dataSize );
    }
}

xmesh_writer::xmesh_writer( const frantic::tstring& path )
    : m_alreadyClosed( false )
    , m_compressionLevel( Z_DEFAULT_COMPRESSION ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    m_rootFilename = frantic::strings::to_wstring( path );
#else
    m_rootFilename = path;
#endif
}

xmesh_writer::~xmesh_writer() {
    // if file was not closed ( e.g. the user canceled ) then remove all the files we wrote
    if( !m_alreadyClosed ) {
        boost::system::error_code errcode;

        boost::filesystem::remove( m_rootFilename, errcode );
        if( errcode )
            FF_LOG( warning ) << "xmesh_writer failed to remove " << m_rootFilename << " : "
                              << frantic::files::to_tstring( errcode.message() );

        for( std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.begin(),
                                                                              itEnd = m_vertexChannels.end();
             it != itEnd; ++it ) {
            boost::filesystem::remove( it->second.vertexPath, errcode );
            if( errcode )
                FF_LOG( warning ) << "xmesh_writer failed to remove " << it->second.vertexPath << " : "
                                  << frantic::files::to_tstring( errcode.message() );
            if( !it->second.facePath.empty() )
                boost::filesystem::remove( it->second.facePath, errcode );
            if( errcode )
                FF_LOG( warning ) << "xmesh_writer failed to remove " << it->second.facePath << " : "
                                  << frantic::files::to_tstring( errcode.message() );
        }

        for( std::map<frantic::tstring, xmesh_face_channel>::const_iterator it = m_faceChannels.begin(),
                                                                            itEnd = m_faceChannels.end();
             it != itEnd; ++it ) {
            boost::filesystem::remove( it->second.facePath, errcode );
            if( errcode )
                FF_LOG( warning ) << "xmesh_writer failed to remove " << it->second.facePath << " : "
                                  << frantic::files::to_tstring( errcode.message() );
        }
    }
}

void xmesh_writer::write_vertex_channel( const frantic::tstring& name, const char* pData, data_type_t type,
                                         std::size_t arity, std::size_t count ) {
    const std::wstring wname = wstring_from_channel_name( name );
    boost::filesystem::path dataPath( frantic::files::replace_extension(
        frantic::files::add_before_sequence_number( m_rootFilename.wstring(),
                                                    name != _T("verts") ? L"channel_" + wname + L"_verts" : L"verts" ),
        L".xmdata" ) );

    // HACK: Do not use frantic::strings::to_string(
    std::string dataType = frantic::strings::to_string( frantic::channels::channel_data_type_str( arity, type ) );
    std::size_t dataSize = frantic::channels::sizeof_channel_data_type( type ) * arity;
    write_xmesh_array_file( dataPath, pData, dataType, dataSize, count, m_compressionLevel );

    { // scope for lock
        mutex_t::scoped_lock lock( m_vertexChannelsMutex );

        xmesh_vertex_channel& vertChannel = m_vertexChannels[name];
        vertChannel.vertexType.first = type;
        vertChannel.vertexType.second = arity;
        vertChannel.vertexCount = count;
        vertChannel.vertexPath = dataPath.native();
    }
}

void xmesh_writer::write_vertex_channel_faces( const frantic::tstring& name, const int* pData, std::size_t count ) {
    const std::wstring wname = wstring_from_channel_name( name );
    boost::filesystem::path dataFacePath( frantic::files::replace_extension(
        frantic::files::add_before_sequence_number( m_rootFilename.wstring(),
                                                    name != _T("verts") ? L"channel_" + wname + L"_faces" : L"faces" ),
        L".xmdata" ) );
    std::string dataFaceType = "int32[1]";
    std::size_t dataFaceSize = sizeof( int );
    write_xmesh_array_file( dataFacePath, (const char*)pData, dataFaceType, dataFaceSize, count, m_compressionLevel );

    { // scope for lock
        mutex_t::scoped_lock lock( m_vertexChannelsMutex );

        xmesh_vertex_channel& vertChannel = m_vertexChannels[name];
        vertChannel.facePath = dataFacePath.native();
        vertChannel.faceCount = count;
    }
}

void xmesh_writer::write_face_channel( const frantic::tstring& name, const char* pData, data_type_t type,
                                       std::size_t arity, std::size_t count ) {
    const std::wstring wname = wstring_from_channel_name( name );
    boost::filesystem::path dataPath( frantic::files::replace_extension(
        frantic::files::add_before_sequence_number( m_rootFilename.wstring(), L"channel_" + wname + L"_facedata" ),
        L".xmdata" ) );
    std::string dataType = frantic::strings::to_string( frantic::channels::channel_data_type_str( arity, type ) );
    std::size_t dataSize = frantic::channels::sizeof_channel_data_type( type ) * arity;
    write_xmesh_array_file( dataPath, pData, dataType, dataSize, count, m_compressionLevel );

    { // scope for lock
        mutex_t::scoped_lock lock( m_faceChannelsMutex );

        xmesh_face_channel& faceChannel = m_faceChannels[name];
        faceChannel.faceType.first = type;
        faceChannel.faceType.second = arity;
        faceChannel.faceCount = count;
        faceChannel.facePath = dataPath.native();
    }
}

void xmesh_writer::set_compression_level( int compressionLevel ) { m_compressionLevel = compressionLevel; }

void xmesh_writer::set_metadata( const frantic::geometry::xmesh_metadata& metadata ) { m_metadata = metadata; }

void xmesh_writer::close() {
    // Make sure we don'y write the file multiple times.
    // TODO: It should probably throw an exception if we try to add channels once the file is written.
    if( m_alreadyClosed )
        return;
    m_alreadyClosed = true;

    tinyxml2::XMLDocument doc( m_rootFilename.c_str() );
    doc.LinkEndChild( doc.NewDeclaration() );

    // Make the single document level, root tag.
    tinyxml2::XMLElement* pXMesh = doc.LinkEndChild( doc.NewElement( "xmesh" ) )->ToElement();

    // Record the version number
    tinyxml2::XMLElement* pVersion = pXMesh->LinkEndChild( doc.NewElement( "version" ) )->ToElement();
    pVersion->LinkEndChild( doc.NewText( "1" ) );

    // Record the geometric vertex count
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.find( _T("verts") );
    if( it == m_vertexChannels.end() )
        throw std::runtime_error(
            "xmesh_writer::close() There must be a vertex channel called 'verts' storing the geometry channel." );

    tinyxml2::XMLElement* pVertCount = pXMesh->LinkEndChild( doc.NewElement( "vertCount" ) )->ToElement();
    pVertCount->LinkEndChild(
        doc.NewText( frantic::strings::to_string_classic<std::string>( it->second.get_vertex_count() ).c_str() ) );

    if( it->second.get_vertex_file_path().empty() ) {
        throw std::runtime_error( "xmesh_writer::close() Error: While writing \"" + m_rootFilename.string() +
                                  "\": No vertex data was provided for the geometry channel." );
    }

    const std::string vertexPathUTF8(
        frantic::strings::to_utf8( it->second.get_vertex_file_path().filename().wstring() ) );
    if( !frantic::strings::is_valid_utf8( vertexPathUTF8 ) ) {
        throw std::runtime_error( "xmesh_writer::close() Internal Error: While writing \"" + m_rootFilename.string() +
                                  "\": Vertex path was not converted to valid UTF-8: \"" +
                                  it->second.get_vertex_file_path().string() + "\"." );
    }
    tinyxml2::XMLElement* pVertFile = pXMesh->LinkEndChild( doc.NewElement( "verts" ) )->ToElement();
    pVertFile->LinkEndChild( doc.NewText( vertexPathUTF8.c_str() ) );

    tinyxml2::XMLElement* pFaceCount = pXMesh->LinkEndChild( doc.NewElement( "faceCount" ) )->ToElement();
    pFaceCount->LinkEndChild(
        doc.NewText( frantic::strings::to_string_classic<std::string>( it->second.get_face_count() ).c_str() ) );

    const std::string facePathUTF8( frantic::strings::to_utf8( it->second.get_face_file_path().filename().wstring() ) );
    if( !frantic::strings::is_valid_utf8( facePathUTF8 ) ) {
        throw std::runtime_error( "xmesh_writer::close() Internal Error: While writing \"" + m_rootFilename.string() +
                                  "\": Face path was not converted to valid UTF-8: \"" +
                                  it->second.get_vertex_file_path().string() + "\"." );
    }
    tinyxml2::XMLElement* pFaceFile = pXMesh->LinkEndChild( doc.NewElement( "faces" ) )->ToElement();
    pFaceFile->LinkEndChild( doc.NewText( facePathUTF8.c_str() ) );

    for( std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.begin(),
                                                                          itEnd = m_vertexChannels.end();
         it != itEnd; ++it ) {
        if( it->first == _T("verts") )
            continue;

        tinyxml2::XMLElement* pChannel = pXMesh->LinkEndChild( doc.NewElement( "vertexChannel" ) )->ToElement();

        if( it->second.get_vertex_file_path().empty() ) {
            throw std::runtime_error( "xmesh_writer::close() Error: While writing \"" + m_rootFilename.string() +
                                      "\": No vertex data was provided for channel \"" +
                                      frantic::strings::to_string( it->first ) + "\"." );
        }

        const std::string channelNameUTF8( utf8_from_channel_name( it->first ) );
        if( !frantic::strings::is_valid_utf8( channelNameUTF8 ) ) {
            throw std::runtime_error( "xmesh_writer::close() Internal Error: While writing \"" +
                                      m_rootFilename.string() +
                                      "\": Vertex channel name was not converted to valid UTF-8: \"" +
                                      frantic::strings::to_string( it->first ) + "\"." );
        }
        tinyxml2::XMLElement* pChannelName = pChannel->LinkEndChild( doc.NewElement( "name" ) )->ToElement();
        pChannelName->LinkEndChild( doc.NewText( channelNameUTF8.c_str() ) );

        tinyxml2::XMLElement* pChannelType = pChannel->LinkEndChild( doc.NewElement( "type" ) )->ToElement();
        pChannelType->LinkEndChild( doc.NewText( frantic::strings::to_string( channel_data_type_str(
            it->second.get_vertex_primitive_arity(), it->second.get_vertex_primitive_type() ) ).c_str() ) );

        tinyxml2::XMLElement* pChannelCount = pChannel->LinkEndChild( doc.NewElement( "vertCount" ) )->ToElement();
        pChannelCount->LinkEndChild(
            doc.NewText( frantic::strings::to_string_classic<std::string>( it->second.get_vertex_count() ).c_str() ) );

        const std::string channelPathUTF8(
            frantic::strings::to_utf8( it->second.get_vertex_file_path().filename().wstring() ) );
        if( !frantic::strings::is_valid_utf8( channelPathUTF8 ) ) {
            throw std::runtime_error( "xmesh_writer::close() Internal Error: While writing \"" +
                                      m_rootFilename.string() +
                                      "\": Vertex channel path was not converted to valid UTF-8: \"" +
                                      it->second.get_vertex_file_path().string() + "\"." );
        }
        tinyxml2::XMLElement* pChannelPath = pChannel->LinkEndChild( doc.NewElement( "vertData" ) )->ToElement();
        pChannelPath->LinkEndChild( doc.NewText( channelPathUTF8.c_str() ) );

        if( it->second.get_face_count() > 0 && !it->second.get_face_file_path().empty() ) {
            tinyxml2::XMLElement* pChannelFaceCount = pChannel->LinkEndChild( doc.NewElement( "faceCount" ) )->ToElement();
            pChannelFaceCount->LinkEndChild(
                doc.NewText( frantic::strings::to_string_classic<std::string>( it->second.get_face_count() ).c_str() ) );

            const std::string channelFacePathUTF8(
                frantic::strings::to_utf8( it->second.get_face_file_path().filename().wstring() ) );
            if( !frantic::strings::is_valid_utf8( channelFacePathUTF8 ) ) {
                throw std::runtime_error( "xmesh_writer::close() Internal Error: While writing \"" +
                                          m_rootFilename.string() +
                                          "\": Vertex channel face path was not converted to valid UTF-8: \"" +
                                          it->second.get_face_file_path().string() + "\"." );
            }
            tinyxml2::XMLElement* pChannelFacePath = pChannel->LinkEndChild( doc.NewElement( "faceData" ) )->ToElement();
            pChannelFacePath->LinkEndChild( doc.NewText( channelFacePathUTF8.c_str() ) );
        }
    }

    for( std::map<frantic::tstring, xmesh_face_channel>::const_iterator it = m_faceChannels.begin(),
                                                                        itEnd = m_faceChannels.end();
         it != itEnd; ++it ) {
        tinyxml2::XMLElement* pChannel = pXMesh->LinkEndChild( doc.NewElement( "faceChannel" ) )->ToElement();

        const std::string channelNameUTF8( utf8_from_channel_name( it->first ) );
        if( !frantic::strings::is_valid_utf8( channelNameUTF8 ) ) {
            throw std::runtime_error( "xmesh_writer::close() Internal Error: While writing \"" +
                                      m_rootFilename.string() +
                                      "\": Face channel name was not converted to valid UTF-8: \"" +
                                      frantic::strings::to_string( it->first ) + "\"." );
        }
        tinyxml2::XMLElement* pChannelName = pChannel->LinkEndChild( doc.NewElement( "name" ) )->ToElement();
        pChannelName->LinkEndChild( doc.NewText( utf8_from_channel_name( it->first ).c_str() ) );

        tinyxml2::XMLElement* pChannelType = pChannel->LinkEndChild( doc.NewElement( "type" ) )->ToElement();
        pChannelType->LinkEndChild( doc.NewText( frantic::strings::to_string(
            channel_data_type_str( it->second.get_face_primitive_arity(), it->second.get_face_primitive_type() ) ).c_str() ) );

        tinyxml2::XMLElement* pChannelCount = pChannel->LinkEndChild( doc.NewElement( "faceCount" ) )->ToElement();
        pChannelCount->LinkEndChild(
            doc.NewText( frantic::strings::to_string_classic<std::string>( it->second.get_face_count() ).c_str() ) );

        const std::string channelPathUTF8(
            frantic::strings::to_utf8( it->second.get_face_file_path().filename().wstring() ) );
        if( !frantic::strings::is_valid_utf8( channelPathUTF8 ) ) {
            throw std::runtime_error( "xmesh_writer::close() Internal Error: While writing \"" +
                                      m_rootFilename.string() +
                                      "\": Face channel face path was not converted to valid UTF-8: \"" +
                                      it->second.get_face_file_path().string() + "\"." );
        }
        tinyxml2::XMLElement* pChannelPath = pChannel->LinkEndChild( doc.NewElement( "faceData" ) )->ToElement();
        pChannelPath->LinkEndChild( doc.NewText( channelPathUTF8.c_str() ) );
    }

    write_xmesh_metadata( doc, m_metadata );

    frantic::files::file_ptr fout( frantic::files::fopen( m_rootFilename, "w" ) );
    const tinyxml2::XMLError result = doc.SaveFile( fout );
    if( result != tinyxml2::XMLError::XML_SUCCESS )
        throw std::runtime_error( "xmesh_writer::close() Error writing XML document to \"" + m_rootFilename.string() +
                                  "\"" );
}

} // namespace geometry
} // namespace frantic
