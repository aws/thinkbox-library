// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/files/files.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/files/zlib_reader.hpp>
#include <frantic/geometry/xmesh_reader.hpp>
#include <frantic/strings/utf8.hpp>
#include <frantic/tinyxml/frantic_tinyxml_utility.hpp>

#include <boost/algorithm/string/case_conv.hpp>

#include <tinyxml2.h>

using namespace frantic::channels;

namespace frantic {
namespace geometry {

// NOTE: Copied from xmesh_sequence_saver.cpp, this implementation should override that one.
// Loads a single array of data.
// Format:
//   "xmeshdat"     8 bytes
//   "<datatag>"    12 bytes (based on the named channel data type, so all 12 bytes may not actually be used)
//   count          4 bytes
//   dataSize       4 bytes
//   <array data>   (count * dataSize) bytes
void load_xmesh_array_file( const boost::filesystem::path& path, char* data, const std::string& correctDataType,
                            size_t correctDataSize, size_t correctDataCount ) {
    frantic::files::zlib_gzip_read_interface stream;

    stream.open( path );

    // Read the header info (a string and the array size)
    char header[9];
    char dataTag[13];

    stream.read( header, 8 );
    header[8] = 0;

    stream.read( dataTag, 12 );
    dataTag[12] = 0;

    std::pair<data_type_t, std::size_t> expectedType =
        frantic::channels::channel_data_type_and_arity_from_string( correctDataType );
    std::pair<data_type_t, std::size_t> actualType =
        frantic::channels::channel_data_type_and_arity_from_string( dataTag );

    if( expectedType.first != actualType.first || expectedType.second != actualType.second ) {
        throw std::runtime_error( "load_xmesh_array_file: Input file \"" + path.string() + "\" has wrong data tag, \"" +
                                  dataTag + "\".  The expected data tag is \"" + correctDataType + "\"." );
    }

    boost::int32_t dataSize = 0;
    boost::int32_t dataCount = 0;
    stream.read( reinterpret_cast<char*>( &dataCount ), 4 );
    stream.read( reinterpret_cast<char*>( &dataSize ), 4 );

    // Perform the appropriate checks
    if( strcmp( header, "xmeshdat" ) != 0 ) {
        throw std::runtime_error( "load_xmesh_array_file: File \"" + path.string() + "\" had invalid header." );
    }
    if( dataSize != (int)correctDataSize ) {
        throw std::runtime_error( "load_xmesh_array_file: File \"" + path.string() +
                                  "\" contains the wrong data size.  Expected " +
                                  boost::lexical_cast<std::string>( correctDataSize ) + ", but got " +
                                  boost::lexical_cast<std::string>( dataSize ) + "." );
    }
    if( dataCount != (int)correctDataCount ) {
        throw std::runtime_error( "load_xmesh_array_file: File \"" + path.string() +
                                  "\" contains the wrong data count.  Expected " +
                                  boost::lexical_cast<std::string>( correctDataCount ) + ", but got " +
                                  boost::lexical_cast<std::string>( dataCount ) + "." );
    }

    // Read the data
    if( dataCount > 0 ) {
        stream.read( data, dataCount * dataSize );
    }
}

/**
 * This function will take an array of int32[3] and convert it to an encoded array of int32[1]
 * @param pOldStyleFaces pointer to the array holding triangle face indices
 * @param numElements the number of int32[1] elements in the array.
 */
void fix_old_style_custom_faces_channel( char* pOldStyleFaces, std::size_t numElements ) {
    int* pData = (int*)pOldStyleFaces;
    for( std::size_t i = 2; i < numElements; i += 3 )
        pData[i] = -pData[i] - 1;
}

namespace {

// Helper functions for interpreting strings in XML files.
//
// We use the utf8 functions if we find encoding="UTF-8" in the xmesh file,
// and the legacy functions otherwise.

// Convert from XML string to native path string

boost::filesystem::path::string_type convert_utf8_string_to_path_string( const char* s ) {
    if( s ) {
        return boost::filesystem::path( frantic::strings::wstring_from_utf8( s ) ).native();
    } else {
        return boost::filesystem::path::string_type();
    }
}

boost::filesystem::path::string_type convert_legacy_string_to_path_string( const char* s ) {
    if( s ) {
#if defined( _WIN32 ) || defined( _WIN64 )
        return boost::filesystem::path( frantic::strings::to_wstring( s ) ).native();
#else
        return boost::filesystem::path( s ).native();
#endif
    } else {
        return boost::filesystem::path::string_type();
    }
}

// Convert from XML string to channel name

// TODO: How are std::string channel names encoded?
// I think they should be UTF-8, but for now we're just passing through the
// string as it was found in the XML file.
frantic::tstring convert_utf8_string_to_channel_name( const char* s ) {
    if( s ) {
#ifdef FRANTIC_USE_WCHAR
        return frantic::strings::wstring_from_utf8( s );
#else
        return s;
#endif
    } else {
        return frantic::tstring();
    }
}

frantic::tstring convert_legacy_string_to_channel_name( const char* s ) {
    if( s ) {
#ifdef FRANTIC_USE_WCHAR
        return frantic::strings::to_wstring( s );
#else
        return s;
#endif
    } else {
        return frantic::tstring();
    }
}

// All "legacy" strings are considered valid.
bool is_valid_legacy_string( const char* ) { return true; }

} // anonymous namespace

xmesh_reader::xmesh_reader( const boost::filesystem::path& path ) {
    m_path = path;
    m_rootPath = path.parent_path();

    // Check that there's an xml file at the given location
    tinyxml2::XMLDocument doc( path.c_str() );
    frantic::files::file_ptr fin( frantic::files::fopen( path, "rb" ) );
    const tinyxml2::XMLError result = doc.LoadFile( fin );
    if( result != tinyxml2::XMLError::XML_SUCCESS ) {
        if( doc.ErrorID() == tinyxml2::XMLError::XML_ERROR_FILE_COULD_NOT_BE_OPENED ) {
            throw frantic::files::file_open_error( "xmesh_reader::xmesh_reader() Failed to open xmesh file \"" +
                                                   path.string() + "\"." );
        } else {
            throw std::runtime_error( "xmesh_reader::xmesh_reader() Failed to load or parse xmesh file \"" +
                                      path.string() + "\"." );
        }
    }

    tinyxml2::XMLElement* pXMesh = doc.RootElement();
    if( !pXMesh )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" did not have a root element." );

    // Interpretation of channel name and file name strings depends on
    // whether the encoding is "UTF-8".
    bool useUTF8 = false;
    const std::string encoding = boost::algorithm::to_upper_copy( frantic::tinyxml::get_encoding_declaration( doc ) );
    if( encoding == "UTF-8" ) {
        useUTF8 = true;
    }
    bool ( *is_valid_external_string )( const char* ) =
        useUTF8 ? static_cast<bool ( * )( const char* )>( &frantic::strings::is_valid_utf8 ) : is_valid_legacy_string;
    boost::filesystem::path::string_type ( *to_internal_file_name )( const char* ) =
        useUTF8 ? convert_utf8_string_to_path_string : convert_legacy_string_to_path_string;
    frantic::tstring ( *to_internal_channel_name )( const char* ) =
        useUTF8 ? convert_utf8_string_to_channel_name : convert_legacy_string_to_channel_name;

    tinyxml2::XMLHandle docHandle( pXMesh );

    m_version = 0;

    tinyxml2::XMLText* pVersion = docHandle.FirstChildElement( "version" ).FirstChild().ToText();
    if( pVersion )
        m_version = boost::lexical_cast<int>( pVersion->Value() );

    const int maxSupportedVersion = 1;
    if( m_version > maxSupportedVersion ) {
        throw std::runtime_error( "xmesh_reader::xmesh_reader() Unsupported version in file: \"" + path.string() +
                                  "\".  "
                                  "File version: " +
                                  boost::lexical_cast<std::string>( m_version ) +
                                  ".  "
                                  "Maximum supported version: " +
                                  boost::lexical_cast<std::string>( maxSupportedVersion ) +
                                  ".  "
                                  "Please contact support@thinkboxsoftware.com for information on newer versions." );
    }

    tinyxml2::XMLText* pVertCount = docHandle.FirstChildElement( "vertCount" ).FirstChild().ToText();
    if( !pVertCount )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" has a missing or invalid <vertCount> tag." );

    tinyxml2::XMLText* pVertFile = docHandle.FirstChildElement( "verts" ).FirstChild().ToText();
    if( !pVertFile )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" has a missing or invalid <verts> tag." );
    if( !is_valid_external_string( pVertFile->Value() ) )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" has an invalid string in its <verts> tag." );

    tinyxml2::XMLText* pFaceCount = docHandle.FirstChildElement( "faceCount" ).FirstChild().ToText();
    if( !pFaceCount )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" has a missing or invalid <faceCount> tag." );

    // It is possible for the faces tag to be empty, so we need to specifically handle that but TinyXML will
    // not create a text node with no text. Instead it just barfs.
    tinyxml2::XMLHandle faceFileHandle = docHandle.FirstChildElement( "faces" );
    if( !faceFileHandle.ToNode() )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" has a missing or invalid <faces> tag." );
    tinyxml2::XMLText* pFaceFile = faceFileHandle.FirstChild().ToText();
    if( pFaceFile && !is_valid_external_string( pFaceFile->Value() ) )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" has an invalid string in its <faces> tag." );

    xmesh_vertex_channel& geomChannel = m_vertexChannels[_T("verts")];
    geomChannel.vertexPath = m_rootPath / to_internal_file_name( pVertFile->Value() );
    geomChannel.vertexType.first = frantic::channels::data_type_float32;
    geomChannel.vertexType.second = 3;
    geomChannel.vertexCount = boost::lexical_cast<std::size_t>( pVertCount->Value() );
    geomChannel.faceCount = boost::lexical_cast<std::size_t>( pFaceCount->Value() );
    if( pFaceFile )
        geomChannel.facePath = m_rootPath / to_internal_file_name( pFaceFile->Value() );
    else if( geomChannel.faceCount != 0 )
        throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                  "\" has a missing or invalid <faces> tag." );

    // Handle the old version way of storing the number of triangles.
    if( m_version < 1 )
        geomChannel.faceCount *= 3;

    std::size_t curVertexChannel = 0;
    std::size_t curFaceChannel = 0;

    // Iterate over all 'vertexChannel' tags, getting information about the per-vertex and the
    // custom face channels.
    tinyxml2::XMLElement* pVertexChannel = docHandle.FirstChildElement( "vertexChannel" ).ToElement();
    while( pVertexChannel ) {
        tinyxml2::XMLHandle channelHandle( pVertexChannel );

        tinyxml2::XMLText* pChannelType = channelHandle.FirstChildElement( "type" ).FirstChild().ToText();
        if( !pChannelType )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " had a missing or invalid <type> tag." );

        tinyxml2::XMLText* pChannelName = channelHandle.FirstChildElement( "name" ).FirstChild().ToText();
        if( !pChannelName )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " had a missing or invalid <name> tag." );
        if( !is_valid_external_string( pChannelName->Value() ) )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " has an invalid string in its <name> tag." );

        tinyxml2::XMLText* pChannelVertCount = channelHandle.FirstChildElement( "vertCount" ).FirstChild().ToText();
        if( !pChannelVertCount )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " had a missing or invalid <vertCount> tag." );

        tinyxml2::XMLText* pChannelVertFile = channelHandle.FirstChildElement( "vertData" ).FirstChild().ToText();
        if( !pChannelVertFile )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " had a missing or invalid <vertData> tag." );
        if( !is_valid_external_string( pChannelVertFile->Value() ) )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " has an invalid string in its <vertData> tag." );

        // We need to either have both or neither of these tags, which is the case for vertex channels
        // and custom face channels.
        tinyxml2::XMLText* pChannelFaceCount = channelHandle.FirstChildElement( "faceCount" ).FirstChild().ToText();
        tinyxml2::XMLText* pChannelFaceFile = channelHandle.FirstChildElement( "faceData" ).FirstChild().ToText();
        if( ( pChannelFaceCount == NULL ) ^ ( pChannelFaceFile == NULL ) )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " had invalid <faceCount> and <faceData> tags." );
        if( pChannelFaceFile && !is_valid_external_string( pChannelFaceFile->Value() ) )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on vertexChannel " + boost::lexical_cast<std::string>( curVertexChannel ) +
                                      " has an invalid string in its <faceData> tag." );

        frantic::tstring name = to_internal_channel_name( pChannelName->Value() );
        std::map<frantic::tstring, xmesh_vertex_channel>::iterator it = m_vertexChannels.find( name );
        if( it != m_vertexChannels.end() )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The vertex channel \"" +
                                      frantic::strings::to_string( name ) + "\" was defined more than once in \"" +
                                      path.string() + "\"" );

        std::pair<data_type_t, std::size_t> typeInfo = channel_data_type_and_arity_from_string( pChannelType->Value() );

        xmesh_vertex_channel& vertChannel = m_vertexChannels[name];
        vertChannel.vertexPath = m_rootPath / to_internal_file_name( pChannelVertFile->Value() );
        vertChannel.vertexType = typeInfo;
        vertChannel.vertexCount = boost::lexical_cast<std::size_t>( pChannelVertCount->Value() );
        if( pChannelFaceFile ) {
            vertChannel.facePath = m_rootPath / to_internal_file_name( pChannelFaceFile->Value() );
            vertChannel.faceCount = boost::lexical_cast<std::size_t>( pChannelFaceCount->Value() );

            // Handle the old version way of storing the number of triangles.
            if( m_version < 1 )
                vertChannel.faceCount *= 3;
        } else {
            vertChannel.facePath.clear();
            vertChannel.faceCount = 0;
        }

        ++curVertexChannel;
        pVertexChannel = pVertexChannel->NextSiblingElement( "vertexChannel" );
    }

    // Iterate over all 'faceChannel' tags, getting information about the per-face channels.
    tinyxml2::XMLElement* pFaceChannel = docHandle.FirstChildElement( "faceChannel" ).ToElement();
    while( pFaceChannel ) {
        tinyxml2::XMLHandle channelHandle( pFaceChannel );

        tinyxml2::XMLText* pChannelType = channelHandle.FirstChildElement( "type" ).FirstChild().ToText();
        if( !pChannelType )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on faceChannel " + boost::lexical_cast<std::string>( curFaceChannel ) +
                                      " had a missing or invalid <type> tag." );

        tinyxml2::XMLText* pChannelName = channelHandle.FirstChildElement( "name" ).FirstChild().ToText();
        if( !pChannelName )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on faceChannel " + boost::lexical_cast<std::string>( curFaceChannel ) +
                                      " had a missing or invalid <name> tag." );
        if( !is_valid_external_string( pChannelName->Value() ) )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on faceChannel " + boost::lexical_cast<std::string>( curFaceChannel ) +
                                      " has an invalid string in its <name> tag." );

        tinyxml2::XMLText* pChannelFaceCount = channelHandle.FirstChildElement( "faceCount" ).FirstChild().ToText();
        if( !pChannelFaceCount )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on faceChannel " + boost::lexical_cast<std::string>( curFaceChannel ) +
                                      " had a missing or invalid <faceCount> tag." );

        tinyxml2::XMLText* pChannelFaceFile = channelHandle.FirstChildElement( "faceData" ).FirstChild().ToText();
        if( !pChannelFaceFile )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on faceChannel " + boost::lexical_cast<std::string>( curFaceChannel ) +
                                      " had a missing or invalid <faceData> tag." );
        if( !is_valid_external_string( pChannelFaceFile->Value() ) )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The xmesh file \"" + path.string() +
                                      "\" on faceChannel " + boost::lexical_cast<std::string>( curFaceChannel ) +
                                      " has an invalid string in its <faceData> tag." );

        frantic::tstring name = to_internal_channel_name( pChannelName->Value() );
        std::map<frantic::tstring, xmesh_face_channel>::iterator it = m_faceChannels.find( name );
        if( it != m_faceChannels.end() )
            throw std::runtime_error( "xmesh_reader::xmesh_reader() The per-face channel \"" +
                                      frantic::strings::to_string( name ) + "\" was defined more than once in \"" +
                                      path.string() + "\"" );

        std::pair<data_type_t, std::size_t> typeInfo = channel_data_type_and_arity_from_string( pChannelType->Value() );

        xmesh_face_channel& faceChannel = m_faceChannels[name];
        faceChannel.facePath = m_rootPath / to_internal_file_name( pChannelFaceFile->Value() );
        faceChannel.faceType = typeInfo;
        faceChannel.faceCount = boost::lexical_cast<std::size_t>( pChannelFaceCount->Value() );

        ++curFaceChannel;
        pFaceChannel = pFaceChannel->NextSiblingElement( "faceChannel" );
    }

    read_xmesh_metadata( doc, m_metadata );
}

void xmesh_reader::get_vertex_channel_names( std::vector<frantic::tstring>& outNames ) const {
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.begin();
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator itEnd = m_vertexChannels.end();
    for( ; it != itEnd; ++it ) {
        if( it->first != _T("verts") )
            outNames.push_back( it->first );
    }
}

void xmesh_reader::get_face_channel_names( std::vector<frantic::tstring>& outNames ) const {
    std::map<frantic::tstring, xmesh_face_channel>::const_iterator it = m_faceChannels.begin();
    std::map<frantic::tstring, xmesh_face_channel>::const_iterator itEnd = m_faceChannels.end();
    for( ; it != itEnd; ++it )
        outNames.push_back( it->first );
}

bool xmesh_reader::has_vertex_channel( const frantic::tstring& name ) const {
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.find( name );
    return ( it != m_vertexChannels.end() );
}

bool xmesh_reader::has_custom_faces( const frantic::tstring& vertexChannelName ) const {
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.find( vertexChannelName );
    if( it == m_vertexChannels.end() )
        throw std::runtime_error( "xmesh_reader::has_custom_faces() The vertex channel \"" +
                                  frantic::strings::to_string( vertexChannelName ) + "\" does not exist" );
    return ( it->second.faceCount > 0 );
}

bool xmesh_reader::has_face_channel( const frantic::tstring& name ) const {
    std::map<frantic::tstring, xmesh_face_channel>::const_iterator it = m_faceChannels.find( name );
    return ( it != m_faceChannels.end() );
}

const xmesh_vertex_channel& xmesh_reader::get_vertex_channel( const frantic::tstring& vertexChannelName ) const {
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.find( vertexChannelName );
    if( it == m_vertexChannels.end() )
        throw std::runtime_error( "xmesh_reader::get_vertex_channel() The vertex channel \"" +
                                  frantic::strings::to_string( vertexChannelName ) + "\" does not exist" );
    return it->second;
}

const xmesh_face_channel& xmesh_reader::get_face_channel( const frantic::tstring& faceChannelName ) const {
    std::map<frantic::tstring, xmesh_face_channel>::const_iterator it = m_faceChannels.find( faceChannelName );
    if( it == m_faceChannels.end() )
        throw std::runtime_error( "xmesh_reader::get_face_channel() The face channel \"" +
                                  frantic::strings::to_string( faceChannelName ) + "\" does not exist" );
    return it->second;
}

void xmesh_reader::load_vertex_channel( const frantic::tstring& name, char* pData, data_type_t expectedType,
                                        std::size_t expectedArity, std::size_t expectedCount ) const {
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.find( name );
    if( it == m_vertexChannels.end() )
        throw std::runtime_error( "xmesh_reader::load_vertex_channel() The vertex channel \"" +
                                  frantic::strings::to_string( name ) + "\" does not exist" );

    // HACK: Use tstring here!
    std::string dataType =
        frantic::strings::to_string( frantic::channels::channel_data_type_str( expectedArity, expectedType ) );
    std::size_t dataSize = frantic::channels::sizeof_channel_data_type( expectedType ) * expectedArity;
    load_xmesh_array_file( it->second.vertexPath, pData, dataType, dataSize, expectedCount );
}

void xmesh_reader::load_vertex_channel_faces( const frantic::tstring& name, char* pData,
                                              std::size_t expectedCount ) const {
    std::map<frantic::tstring, xmesh_vertex_channel>::const_iterator it = m_vertexChannels.find( name );
    if( it == m_vertexChannels.end() )
        throw std::runtime_error( "xmesh_reader::load_vertex_channel_faces() The vertex channel \"" +
                                  frantic::strings::to_string( name ) + "\" does not exist" );

    if( it->second.faceCount == 0 || it->second.facePath.empty() )
        throw std::runtime_error( "xmesh_reader::load_vertex_channel_faces() The vertex channel \"" +
                                  frantic::strings::to_string( name ) + "\" does not have custom faces" );

    if( m_version == 0 ) {
        // Support legacy face encodings that assume triangles.
        std::string dataType = "int32[3]";
        std::size_t dataSize = sizeof( int[3] );
        load_xmesh_array_file( it->second.facePath, pData, dataType, dataSize, expectedCount / 3 );
        fix_old_style_custom_faces_channel( pData, expectedCount );
    } else {
        std::string dataType = "int32[1]";
        std::size_t dataSize = sizeof( int );
        load_xmesh_array_file( it->second.facePath, pData, dataType, dataSize, expectedCount );
    }
}

void xmesh_reader::load_face_channel( const frantic::tstring& name, char* pData, data_type_t expectedType,
                                      std::size_t expectedArity, std::size_t expectedCount ) const {
    std::map<frantic::tstring, xmesh_face_channel>::const_iterator it = m_faceChannels.find( name );
    if( it == m_faceChannels.end() )
        throw std::runtime_error( "xmesh_reader::load_face_channel() The face channel \"" +
                                  frantic::strings::to_string( name ) + "\" does not exist" );

    // HACK: Use tstring here!
    std::string dataType =
        frantic::strings::to_string( frantic::channels::channel_data_type_str( expectedArity, expectedType ) );
    std::size_t dataSize = frantic::channels::sizeof_channel_data_type( expectedType ) * expectedArity;
    load_xmesh_array_file( it->second.facePath, pData, dataType, dataSize, expectedCount );
}

} // namespace geometry
} // namespace frantic
