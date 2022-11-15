// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/ply_reader.hpp>

#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/rply.hpp>
#include <frantic/logging/logging_level.hpp>

using frantic::graphics::vector3f;

static void ply_error_cb_ff_log( p_ply ply, const char* message ) {
    void* pData;
    ply_get_ply_user_data( ply, &pData, 0 );
    FF_LOG( error ) << "ply_reader_file error: " << message << " in " << *(const frantic::tstring*)pData << std::endl;
}

namespace frantic {
namespace geometry {

namespace {

// RAII wrapper for p_ply
class ply_reader_file : boost::noncopyable {
  public:
    ply_reader_file( const frantic::tstring& filename )
        : m_filename( filename ) {
        p_ply_error_cb errorCallback = &ply_error_cb_ff_log;
        m_file = ply_open( filename.c_str(), errorCallback, 0, &m_filename );
        if( !m_file ) {
            throw std::runtime_error( "ply_reader_file Error: unable to open file "
                                      "\"" +
                                      frantic::strings::to_string( filename ) + "\"" );
        }
        if( !ply_read_header( m_file ) ) {
            throw std::runtime_error( "ply_reader_file: There is an error with the header of file "
                                      "\"" +
                                      frantic::strings::to_string( filename ) + "\"" );
        }
    }

    ~ply_reader_file() { close(); }

    void close() {
        if( m_file ) {
            ply_close( m_file );
            m_file = 0;
        }
    }

    operator p_ply() { return m_file; }

    operator const p_ply() const { return m_file; }

    const frantic::tstring& get_filename() const { return m_filename; }

  private:
    p_ply m_file;

    frantic::tstring m_filename;
};

void get_property_names( std::set<std::string>& outPropertyNames, const frantic::tstring& filename,
                         p_ply_element element ) {
    outPropertyNames.clear();

    for( p_ply_property i = ply_get_next_property( element, NULL ); i; i = ply_get_next_property( element, i ) ) {
        const char* cname;
        int success = ply_get_property_info( i, &cname, 0, 0, 0 );
        if( !success ) {
            throw std::runtime_error( "get_property_names: Error getting property name "
                                      "in file: \"" +
                                      frantic::strings::to_string( filename ) );
        }
        if( !cname ) {
            throw std::runtime_error( "read_header: Error property name is NULL "
                                      "in file: \"" +
                                      frantic::strings::to_string( filename ) + "\"" );
        }
        outPropertyNames.insert( cname );
    }
}

// PLY callback function for vertex position, texture, and normal properties
int ply_vertex_callback( p_ply_argument arg ) {
    long vertexIndex, index;
    void* pdata;
    ply_get_argument_element( arg, NULL, &vertexIndex );
    ply_get_argument_user_data( arg, &pdata, &index );

    std::vector<vector3f>& plyVerts = *reinterpret_cast<std::vector<vector3f>*>( pdata );

    if( vertexIndex >= static_cast<long>( plyVerts.size() ) ) {
        plyVerts.resize( vertexIndex + 1 );
    }

    plyVerts[vertexIndex][index] = (float)ply_get_argument_value( arg );

    return 1;
}

// PLY callback function for vertex color
int ply_color_callback( p_ply_argument arg ) {
    long vertexIndex, index;
    void* pdata;
    ply_get_argument_element( arg, NULL, &vertexIndex );
    ply_get_argument_user_data( arg, &pdata, &index );

    std::vector<vector3f>& plyColor = *reinterpret_cast<std::vector<vector3f>*>( pdata );

    if( vertexIndex >= static_cast<long>( plyColor.size() ) ) {
        plyColor.resize( vertexIndex + 1 );
    }

    plyColor[vertexIndex][index] = (float)ply_get_argument_value( arg ) / 255;

    return 1;
}

struct face_callback_data : boost::noncopyable {
    face_callback_data( std::vector<int>& faceIndices, std::vector<int>& faceSizes )
        : faceIndices( faceIndices )
        , faceSizes( faceSizes ) {}

    std::vector<int>& faceIndices;
    std::vector<int>& faceSizes;
};

// PLY callback function for the faces
int ply_face_callback( p_ply_argument arg ) {
    long length, index;
    void* pdata;
    ply_get_argument_property( arg, NULL, &length, &index );

    ply_get_argument_user_data( arg, &pdata, NULL );

    face_callback_data* data = reinterpret_cast<face_callback_data*>( pdata );
    assert( data );

    if( index == 0 ) {
        data->faceSizes.push_back( length );
    }

    if( index >= 0 && index < length ) {
        data->faceIndices.push_back( static_cast<int>( ply_get_argument_value( arg ) ) );
    }

    return 1;
}

template <typename T>
class raw_byte_buffer_helper {
    frantic::graphics::raw_byte_buffer m_buffer;
    std::size_t m_size;

  public:
    raw_byte_buffer_helper()
        : m_size( 0 ) {}
    std::size_t size() const { return m_size; }
    void push_back( const T& val ) {
        m_buffer.add_element( sizeof( T ) );
        memcpy( m_buffer.ptr_at( m_size * sizeof( T ) ), &val, sizeof( val ) );
        ++m_size;
    }
    void steal_raw_byte_buffer( frantic::graphics::raw_byte_buffer& outBuffer ) {
        outBuffer.swap( m_buffer );
        m_buffer.clear();
        m_size = 0;
    }
    void clear() {
        m_buffer.clear();
        m_size = 0;
    }
};

void add_vertex_channel( frantic::geometry::polymesh3_ptr mesh, const frantic::tstring& channelName,
                         const std::vector<frantic::graphics::vector3f>& data ) {
    if( !mesh ) {
        throw std::runtime_error( "add_vertex_channel Error: mesh is NULL" );
    }
    if( mesh->vertex_count() != data.size() ) {
        throw std::runtime_error( "add_vertex_channel Error: "
                                  "mismatch between mesh vertex count and data count" );
    }

    raw_byte_buffer_helper<vector3f> bufferHelper;
    for( std::size_t i = 0, ie = data.size(); i < ie; ++i ) {
        bufferHelper.push_back( data[i] );
    }
    frantic::graphics::raw_byte_buffer buffer;
    bufferHelper.steal_raw_byte_buffer( buffer );
    mesh->add_vertex_channel( channelName, frantic::channels::data_type_float32, 3, buffer );
}

} // anonymous namespace

class ply_reader_header {
  public:
    ply_reader_header( ply_reader_file& ply )
        : vertexCount( 0 )
        , faceCount( 0 ) {
        const std::string filename( frantic::strings::to_string( ply.get_filename() ) );

        p_ply_element element = ply_get_next_element( ply, NULL );
        for( ; element; element = ply_get_next_element( ply, element ) ) {
            const char* cname = 0;
            long count = 0;
            int success = ply_get_element_info( element, &cname, &count );
            if( !success ) {
                throw std::runtime_error( "ply_reader_header: Error getting element info "
                                          "in file: \"" +
                                          filename + "\"" );
            }
            if( !cname ) {
                throw std::runtime_error( "ply_reader_header: Error element name is NULL "
                                          "in file: \"" +
                                          filename + "\"" );
            }
            const std::string name( cname );
            if( name == "vertex" ) {
                vertexCount = count;
                std::set<std::string> propertyNames;
                get_property_names( propertyNames, ply.get_filename(), element );

                if( propertyNames.count( "x" ) ) {
                    add_vertex_channel( _T( "Position" ) );
                }
                if( propertyNames.count( "red" ) ) {
                    add_vertex_channel( _T( "Color" ) );
                }
                if( propertyNames.count( "nx" ) ) {
                    add_vertex_channel( _T( "Normal" ) );
                }
                if( propertyNames.count( "u" ) ) {
                    add_vertex_channel( _T( "TextureCoord" ) );
                }
            } else if( name == "face" ) {
                faceCount = count;
            }
        }
    }

    bool has_vertex_channel( const frantic::tstring& channelName ) const {
        return vertexChannels.count( channelName ) != 0;
    }

    bool operator!=( const ply_reader_header& other ) const {
        return vertexCount != other.vertexCount || faceCount != other.faceCount ||
               vertexChannels != other.vertexChannels;
    }

    typedef std::map<frantic::tstring, ply_reader::channel_type_t> vertex_channels_t;

    std::size_t vertexCount;
    std::size_t faceCount;
    vertex_channels_t vertexChannels;

  private:
    void add_vertex_channel( const frantic::tstring& channelName ) {
        vertexChannels[channelName] = ply_reader::channel_type_t( frantic::channels::data_type_float32, 3 );
    }
};

ply_reader::ply_reader( const frantic::tstring& filename )
    : m_filename( filename ) {
    ply_reader_file ply( filename );

    m_header.reset( new ply_reader_header( ply ) );

    if( !m_header->has_vertex_channel( _T( "Position" ) ) ) {
        throw std::runtime_error( "ply_reader Error reading file "
                                  "\"" +
                                  frantic::strings::to_string( filename ) +
                                  "\": "
                                  "file is missing x, y, or z property." );
    }
}

const frantic::tstring& ply_reader::get_filename() const { return m_filename; }

std::size_t ply_reader::get_vertex_count() const { return get_header().vertexCount; }

std::size_t ply_reader::get_face_count() const { return get_header().faceCount; }

void ply_reader::get_vertex_channel_names( std::vector<frantic::tstring>& outNames ) const {
    typedef ply_reader_header::vertex_channels_t vertex_channels_t;
    const vertex_channels_t& vertexChannels = get_header().vertexChannels;

    outNames.clear();

    for( vertex_channels_t::const_iterator i = vertexChannels.begin(), ie = vertexChannels.end(); i != ie; ++i ) {
        const frantic::tstring& channelName = i->first;
        if( channelName != _T( "Position" ) ) {
            outNames.push_back( channelName );
        }
    }
}

ply_reader::channel_type_t ply_reader::get_vertex_channel_type( const frantic::tstring& name ) const {
    typedef ply_reader_header::vertex_channels_t vertex_channels_t;
    const vertex_channels_t& vertexChannels = get_header().vertexChannels;

    vertex_channels_t::const_iterator i = vertexChannels.find( name );
    if( i == vertexChannels.end() ) {
        throw std::runtime_error( "ply_reader::get_vertex_channel_type Error: no channel named "
                                  "\"" +
                                  frantic::strings::to_string( name ) + "\"" );
    }

    return i->second;
}

frantic::geometry::polymesh3_ptr ply_reader::read_polymesh3() const {
    const std::string filename( frantic::strings::to_string( get_filename() ) );

    ply_reader_file ply( get_filename() );

    ply_reader_header rereadHeader( ply );

    if( get_header() != rereadHeader ) {
        throw std::runtime_error( "ply_reader::read_polymesh3 Error: "
                                  "mismatch between original header and re-read header "
                                  "in file: \"" +
                                  filename + "\"" );
    }

    // Initialize lists
    std::vector<vector3f> plyVerts;
    std::vector<vector3f> plyColor;
    std::vector<vector3f> plyTexture;
    std::vector<vector3f> plyNormal;
    std::vector<int> faceIndices;
    std::vector<int> faceCounts;
    face_callback_data faceCallbackData( faceIndices, faceCounts );

    // Setup the read callbacks
    long vertexCount = ply_set_read_cb( ply, "vertex", "x", ply_vertex_callback, &plyVerts, 0 );
    ply_set_read_cb( ply, "vertex", "y", ply_vertex_callback, &plyVerts, 1 );
    ply_set_read_cb( ply, "vertex", "z", ply_vertex_callback, &plyVerts, 2 );

    plyVerts.reserve( get_header().vertexCount );

    if( (size_t)vertexCount != get_header().vertexCount ) {
        throw std::runtime_error( "ply_reader::read_polymesh3 Internal Error: mismatch between "
                                  "header vertex count "
                                  "(" +
                                  boost::lexical_cast<std::string>( get_header().vertexCount ) +
                                  ") and "
                                  "callback vertex count "
                                  "(" +
                                  boost::lexical_cast<std::string>( vertexCount ) +
                                  ") "
                                  "in file: \"" +
                                  filename + "\"" );
    }

    if( get_header().has_vertex_channel( _T( "Color" ) ) ) {
        plyColor.reserve( vertexCount );
        ply_set_read_cb( ply, "vertex", "red", ply_color_callback, &plyColor, 0 );
        ply_set_read_cb( ply, "vertex", "green", ply_color_callback, &plyColor, 1 );
        ply_set_read_cb( ply, "vertex", "blue", ply_color_callback, &plyColor, 2 );
    }

    if( get_header().has_vertex_channel( _T( "TextureCoord" ) ) ) {
        plyTexture.reserve( vertexCount );
        ply_set_read_cb( ply, "vertex", "u", ply_vertex_callback, &plyTexture, 0 );
        ply_set_read_cb( ply, "vertex", "v", ply_vertex_callback, &plyTexture, 1 );
        ply_set_read_cb( ply, "vertex", "t", ply_vertex_callback, &plyTexture, 2 );
    }

    if( get_header().has_vertex_channel( _T( "Normal" ) ) ) {
        plyNormal.reserve( vertexCount );
        ply_set_read_cb( ply, "vertex", "nx", ply_vertex_callback, &plyNormal, 0 );
        ply_set_read_cb( ply, "vertex", "ny", ply_vertex_callback, &plyNormal, 1 );
        ply_set_read_cb( ply, "vertex", "nz", ply_vertex_callback, &plyNormal, 2 );
    }

    long faceCount = ply_set_read_cb( ply, "face", "vertex_indices", ply_face_callback, &faceCallbackData, 0 );
    faceIndices.reserve( 3 * faceCount );
    faceCounts.reserve( faceCount );

    if( (size_t)faceCount != get_header().faceCount ) {
        throw std::runtime_error( "ply_reader::read_polymesh3 Internal Error: mismatch between "
                                  "header face count "
                                  "(" +
                                  boost::lexical_cast<std::string>( get_header().faceCount ) +
                                  ") and "
                                  "callback face count "
                                  "(" +
                                  boost::lexical_cast<std::string>( faceCount ) +
                                  ") "
                                  "in file: \"" +
                                  filename + "\"" );
    }

    if( !ply_read( ply ) ) {
        throw std::runtime_error( "load_ply_polymesh_file: Could not read file \"" + filename + "\"" );
    }

    polymesh3_builder builder;

    // Add faces to the builder
    std::size_t offset = 0;
    for( std::size_t i = 0, ie = faceCounts.size(); i < ie; ++i ) {
        const std::size_t faceCornerCount = faceCounts[i];
        if( offset + faceCornerCount > faceIndices.size() ) {
            throw std::runtime_error( "load_ply_polymesh_file: Error reading file \"" + filename +
                                      "\": "
                                      "mismatch between face index array and face size array" );
        }
        builder.add_polygon( &faceIndices[offset], faceCornerCount );
        offset += faceCornerCount;
    }
    // Add vertices
    for( std::vector<vector3f>::const_iterator i = plyVerts.begin(); i != plyVerts.end(); ++i ) {
        builder.add_vertex( *i );
    }

    // Build mesh
    polymesh3_ptr mesh = builder.finalize();

    // Add 'Color' channel if needed
    if( plyColor.size() == plyVerts.size() ) {
        add_vertex_channel( mesh, _T( "Color" ), plyColor );
    }

    // Add 'TextureCoord' channel if needed
    if( plyTexture.size() == plyVerts.size() ) {
        add_vertex_channel( mesh, _T( "TextureCoord" ), plyTexture );
    }

    // Add 'Normal' channel if needed
    if( plyNormal.size() == plyVerts.size() ) {
        add_vertex_channel( mesh, _T( "Normal" ), plyNormal );
    }

    return mesh;
}

ply_reader_header& ply_reader::get_header() const {
    if( !m_header ) {
        throw std::runtime_error( "ply_reader::get_header Error: header is NULL" );
    }

    return *m_header;
}

} // namespace geometry
} // namespace frantic
