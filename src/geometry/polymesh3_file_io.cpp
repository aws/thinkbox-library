// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/files/files.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/ply_reader.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_file_io.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/xmesh_reader.hpp>
#include <frantic/geometry/xmesh_writer.hpp>
#include <frantic/locale/locale.hpp>

using frantic::graphics::vector3f;

namespace frantic {
namespace geometry {

polymesh3_ptr load_polymesh_file( const frantic::tstring& path ) {
    frantic::tstring type = frantic::strings::to_lower( frantic::files::extension_from_path( path ) );
    if( type == _T(".obj") ) {
        return load_obj_polymesh_file( path );
    } else if( type == _T(".xmesh") ) {
        return load_xmesh_polymesh_file( path );
    } else {
        throw std::runtime_error( "load_polymesh_file: Didn't recognize the file format of the input mesh file \"" +
                                  frantic::strings::to_string( path ) + "\"" );
    }
}

polymesh3_ptr load_xmesh_polymesh_file( const frantic::tstring& filename ) {
#if defined( _WIN32 ) || defined( _WIN64 )
    boost::filesystem::path path( frantic::strings::to_wstring( filename ) );
#else
    boost::filesystem::path path( filename );
#endif

    xmesh_reader fileReader( path );

    frantic::channels::channel_propagation_policy allChannels;

    return load_xmesh_polymesh_file( fileReader, allChannels, true );
}

polymesh3_ptr load_xmesh_polymesh_file( const xmesh_reader& fileReader,
                                        const frantic::channels::channel_propagation_policy& cpp, bool loadFaces ) {
    using namespace frantic::channels;

    const xmesh_vertex_channel& geomCh = fileReader.get_vertex_channel( _T("verts") );
    std::size_t numVerts = geomCh.get_vertex_count();
    std::size_t numFaceElements = loadFaces ? geomCh.get_face_count() : 0;

    frantic::graphics::raw_byte_buffer vertBuffer;
    vertBuffer.resize( sizeof( float[3] ) * numVerts );

    std::vector<int> faceIndexBuffer, faceEndIndexBuffer;
    faceIndexBuffer.resize( numFaceElements );
    faceEndIndexBuffer.reserve( 5000 );

    if( numVerts > 0 )
        fileReader.load_vertex_channel( _T("verts"), vertBuffer.ptr_at( 0 ), data_type_float32, 3, numVerts );
    if( numFaceElements > 0 )
        fileReader.load_vertex_channel_faces( _T("verts"), (char*)&faceIndexBuffer[0], numFaceElements );

    std::size_t numPolygons = 0;

    // Run through each face index and record the end of each polygon. Remove the sentinel values as well.
    for( std::size_t i = 0; i < numFaceElements; ++i ) {
        if( faceIndexBuffer[i] < 0 ) {
            faceIndexBuffer[i] = -faceIndexBuffer[i] - 1;
            faceEndIndexBuffer.push_back( (int)( i + 1 ) );
            ++numPolygons;
        }
    }

    polymesh3_ptr pResult = new polymesh3( vertBuffer, faceIndexBuffer, faceEndIndexBuffer );

    polymesh3_const_vertex_accessor<frantic::graphics::vector3f> geomAccessor =
        pResult->get_const_vertex_accessor<frantic::graphics::vector3f>( _T("verts" ) );

    // Grab the vertex channel names. This does not include the 'verts' channel so we don't need to check for
    // it.
    std::vector<frantic::tstring> channelNames;
    fileReader.get_vertex_channel_names( channelNames );
    cpp.filter_channel_vector( channelNames );

    for( std::size_t i = 0, iEnd = channelNames.size(); i < iEnd; ++i ) {
        const xmesh_vertex_channel& ch = fileReader.get_vertex_channel( channelNames[i] );
        std::size_t numChannelVerts = ch.get_vertex_count();
        std::size_t numChannelFaceElements = loadFaces ? ch.get_face_count() : 0;
        std::pair<frantic::channels::data_type_t, std::size_t> typeInfo = ch.get_vertex_type();

        if( numChannelFaceElements > 0 ) {
            if( numChannelFaceElements != numFaceElements )
                throw std::runtime_error(
                    "load_xmesh_polymesh_file() The channel \"" + frantic::strings::to_string( channelNames[i] ) +
                    "\" did not have the same custom face layout as the geometry channel in file \"" +
                    fileReader.get_path().string() + "\"" );
        } else {
            if( numChannelVerts != numVerts )
                throw std::runtime_error( "load_xmesh_polymesh_file() The channel \"" +
                                          frantic::strings::to_string( channelNames[i] ) +
                                          "\" did not have the same vertex layout as the geometry channel in file \"" +
                                          fileReader.get_path().string() + "\"" );
        }

        // Resize the vertex buffer to hold the new channel's vertex data, and load it from disk.
        vertBuffer.resize( frantic::channels::sizeof_channel_data_type( typeInfo.first ) * typeInfo.second *
                           numChannelVerts );
        fileReader.load_vertex_channel( channelNames[i], vertBuffer.ptr_at( 0 ), typeInfo.first, typeInfo.second,
                                        numChannelVerts );

        if( numChannelFaceElements > 0 ) {
            // Resize the face buffer to hold the new channel's face data, and load it from disk.
            faceIndexBuffer.resize( numChannelFaceElements );
            fileReader.load_vertex_channel_faces( channelNames[i], (char*)&faceIndexBuffer[0], numChannelFaceElements );

            // Remove the encoded polygon sentinel values, and check that this face's size matches the geometry
            // channel.
            std::size_t curPoly = 0;
            std::size_t curPolySize = 0;
            for( std::size_t j = 0; j < numChannelFaceElements; ++j ) {
                ++curPolySize;

                if( faceIndexBuffer[j] < 0 ) {
                    faceIndexBuffer[j] = -faceIndexBuffer[j] - 1;
                    if( curPolySize != geomAccessor.get_face_degree( curPoly ) )
                        throw std::runtime_error(
                            "load_xmesh_polymesh_file() The channel \"" +
                            frantic::strings::to_string( channelNames[i] ) +
                            "\" had mismatched polygon sizes compared to the geometry channel in file \"" +
                            fileReader.get_path().string() + "\"" );
                    curPolySize = 0;
                    curPoly++;
                }
            }

            pResult->add_vertex_channel( channelNames[i], typeInfo.first, typeInfo.second, vertBuffer,
                                         &faceIndexBuffer );
        } else
            pResult->add_vertex_channel( channelNames[i], typeInfo.first, typeInfo.second, vertBuffer );
    }

    if( loadFaces ) {
        // Grab the face channel names (Separate from vertex channels in the xmesh format)
        channelNames.clear();

        fileReader.get_face_channel_names( channelNames );
        cpp.filter_channel_vector( channelNames );

        for( std::size_t i = 0, iEnd = channelNames.size(); i < iEnd; ++i ) {
            const xmesh_face_channel& ch = fileReader.get_face_channel( channelNames[i] );
            std::size_t numChannelFaces = ch.get_face_count();
            std::pair<frantic::channels::data_type_t, std::size_t> typeInfo = ch.get_face_type();

            if( numChannelFaces != numPolygons )
                throw std::runtime_error(
                    "load_xmesh_polymesh_file() The face channel \"" + frantic::strings::to_string( channelNames[i] ) +
                    "\" did not have one entry per-face in file \"" + fileReader.get_path().string() + "\"" );

            vertBuffer.resize( frantic::channels::sizeof_channel_data_type( typeInfo.first ) * typeInfo.second *
                               numChannelFaces );
            fileReader.load_face_channel( channelNames[i], vertBuffer.ptr_at( 0 ), typeInfo.first, typeInfo.second,
                                          numChannelFaces );

            pResult->add_face_channel( channelNames[i], typeInfo.first, typeInfo.second, vertBuffer );
        }
    }

    return pResult;
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

namespace {

float to_float( const std::string& s ) { return static_cast<float>( frantic::locale::strtod_c( s.c_str(), 0 ) ); }

} // anonymous namespace

// copied from load_obj_mesh_file
// TODO : refactor to share logic with load_obj_mesh_file
polymesh3_ptr load_obj_polymesh_file( const frantic::tstring& objFile ) {
    std::ifstream fin( objFile.c_str() );

    if( !fin )
        throw std::runtime_error( "load_obj_polymesh_file: Could not open file \"" +
                                  frantic::strings::to_string( objFile ) + "\" for reading." );

    polymesh3_builder builder;
    std::vector<vector3f> verts;
    raw_byte_buffer_helper<vector3f> tverts, normals;
    std::vector<int> tvertsFaces, normalsFaces;
    bool hasCustomTVertFaces = false;
    bool hasCustomNormalFaces = false;

    // this will be set to non-negative if face definitions are inconsistent (ie. there are four ways to define face
    // data. they are: "face", "face/tvface", "face/tvface/normal", "face//normal")
    int customFacesErrorOnLine = -1;

    std::string line;
    int lineNumber = 0;
    std::vector<std::string> tokens;
    std::vector<std::string> faceTokens;
    bool isFirstFace = true;
    while( getline( fin, line ) ) {
        lineNumber++;
        tokens.clear();
        strings::split( line, tokens );

        // Process only if this line had any tokens
        if( tokens.size() > 0 ) {
            // Skip lines starting with '#'
            if( tokens[0][0] != '#' ) {
                if( tokens[0] == "v" ) {
                    // Add a vertex
                    if( tokens.size() != 4 )
                        throw std::runtime_error(
                            "load_obj_polymesh_file: Line number " + boost::lexical_cast<std::string>( lineNumber ) +
                            "  of \"" + frantic::strings::to_string( objFile ) + "\" had an invalid vertex." );
                    vector3f v( to_float( tokens[1] ), to_float( tokens[2] ), to_float( tokens[3] ) );
                    verts.push_back( v );
                } else if( tokens[0] == "vt" ) {
                    // Add a texture vertex
                    vector3f vt;
                    if( tokens.size() == 3 ) {
                        vt.set( to_float( tokens[1] ), to_float( tokens[2] ), 0 );
                    } else if( tokens.size() == 4 ) {
                        vt.set( to_float( tokens[1] ), to_float( tokens[2] ), to_float( tokens[3] ) );
                    } else {
                        throw std::runtime_error(
                            "load_obj_polymesh_file: Line number " + boost::lexical_cast<std::string>( lineNumber ) +
                            "  of \"" + frantic::strings::to_string( objFile ) + "\" had an invalid texture vertex." );
                    }
                    tverts.push_back( vt );
                } else if( tokens[0] == "vn" ) {
                    // Add a normal
                    if( tokens.size() != 4 )
                        throw std::runtime_error(
                            "load_obj_polymesh_file: Line number " + boost::lexical_cast<std::string>( lineNumber ) +
                            "  of \"" + frantic::strings::to_string( objFile ) + "\" had an invalid normal." );
                    vector3f vn( to_float( tokens[1] ), to_float( tokens[2] ), to_float( tokens[3] ) );
                    normals.push_back( vn );
                } else if( tokens[0] == "f" ) {
                    // Add a face.  Vertices are 1-based in obj files, so we subtract 1 to compensate.
                    if( tokens.size() < 4 )
                        throw std::runtime_error( "load_obj_polymesh_file: Line number " +
                                                  boost::lexical_cast<std::string>( lineNumber ) + "  of \"" +
                                                  frantic::strings::to_string( objFile ) + "\" had an invalid face." );
                    if( isFirstFace ) {
                        // If this is the first face we're adding, it will determine whether vt or vn have custom faces
                        // or share their faces with v
                        faceTokens.clear();
                        strings::split( tokens[1], faceTokens, '/' );
                        if( faceTokens.size() > 1 && faceTokens[1] != "" ) {
                            if( tverts.size() ) {
                                hasCustomTVertFaces = true;
                            } else {
                                throw std::runtime_error( "load_obj_polymesh_file: Line number " +
                                                          boost::lexical_cast<std::string>( lineNumber ) + "  of \"" +
                                                          frantic::strings::to_string( objFile ) +
                                                          "\" had contains a index for the texture coordinate data, "
                                                          "but no \"vt\" data has been defined." );
                            }
                        }
                        if( faceTokens.size() > 2 && faceTokens[2] != "" ) {
                            if( normals.size() ) {
                                hasCustomNormalFaces = true;
                            } else {
                                throw std::runtime_error( "load_obj_polymesh_file: Line number " +
                                                          boost::lexical_cast<std::string>( lineNumber ) + "  of \"" +
                                                          frantic::strings::to_string( objFile ) +
                                                          "\" had contains a index for the normal data, but no \"vn\" "
                                                          "data has been defined." );
                            }
                        }
                        isFirstFace = false;
                    }
                    std::vector<int> face, tvface, tnface;
                    for( unsigned i = 1; i < tokens.size(); ++i ) {
                        if( !hasCustomTVertFaces && !hasCustomNormalFaces ) {
                            if( tokens[i].find_first_of( '/' ) == std::string::npos ) {
                                face.push_back( atoi( tokens[i].c_str() ) - 1 );
                            } else {
                                // get the face from the first token of the string. note, this is an error case since we
                                // encounted a tvert/normal in a face definition when we didn't expect to.
                                faceTokens.clear();
                                strings::split( tokens[i], faceTokens, '/' );
                                face.push_back( atoi( faceTokens[0].c_str() ) - 1 );
                                if( customFacesErrorOnLine == -1 )
                                    customFacesErrorOnLine = lineNumber;
                            }
                        } else {
                            faceTokens.clear();
                            strings::split( tokens[i], faceTokens, '/' );

                            // add vert index
                            face.push_back( atoi( faceTokens[0].c_str() ) - 1 );

                            // add texture coord index
                            bool faceHasTVertIndex = ( faceTokens.size() > 1 && !faceTokens[1].empty() );
                            if( faceHasTVertIndex && hasCustomTVertFaces ) {
                                tvface.push_back( atoi( faceTokens[1].c_str() ) - 1 );
                            } else if( ( faceHasTVertIndex && !hasCustomTVertFaces ) ||
                                       ( !faceHasTVertIndex && hasCustomTVertFaces ) ) { // error cases
                                if( customFacesErrorOnLine == -1 )
                                    customFacesErrorOnLine = lineNumber;
                            }

                            // add normal index
                            bool faceHasNormalIndex = ( faceTokens.size() > 2 );
                            if( faceHasNormalIndex && hasCustomNormalFaces ) {
                                tnface.push_back( atoi( faceTokens[2].c_str() ) - 1 );
                            } else if( ( faceHasNormalIndex && !hasCustomNormalFaces ) ||
                                       ( !faceHasNormalIndex && hasCustomNormalFaces ) ) { // error cases
                                if( customFacesErrorOnLine == -1 )
                                    customFacesErrorOnLine = lineNumber;
                            }
                        }
                    }
                    // Make sure the face has valid indexes, then add it to the mesh
                    for( size_t i = 0; i < face.size(); ++i ) {
                        if( face[i] < -1 )
                            face[i] += ( 1 + static_cast<int>( verts.size() ) );

                        if( face[i] < 0 || face[i] >= (int)verts.size() )
                            throw std::runtime_error( "load_obj_polymesh_file: Line number " +
                                                      boost::lexical_cast<std::string>( lineNumber ) + "  of \"" +
                                                      frantic::strings::to_string( objFile ) +
                                                      "\" had an invalid face, because a v index is out of bounds." );
                    }
                    builder.add_polygon( face );
                    if( hasCustomTVertFaces ) {
                        for( unsigned i = 0; i < tvface.size(); ++i ) {
                            if( tvface[i] < -1 )
                                tvface[i] += ( 1 + static_cast<int>( tverts.size() ) );

                            if( tvface[i] < 0 || tvface[i] >= (int)tverts.size() )
                                throw std::runtime_error(
                                    "load_obj_polymesh_file: Line number " +
                                    boost::lexical_cast<std::string>( lineNumber ) + "  of \"" +
                                    frantic::strings::to_string( objFile ) +
                                    "\" had an invalid face, because a vt index is out of bounds." );
                        }
                        tvertsFaces.insert( tvertsFaces.end(), tvface.begin(), tvface.end() );
                    }
                    if( hasCustomNormalFaces ) {
                        for( unsigned i = 0; i < tnface.size(); ++i ) {
                            if( tnface[i] < -1 )
                                tnface[i] += ( 1 + static_cast<int>( normals.size() ) );

                            if( tnface[i] < 0 || tnface[i] >= (int)normals.size() )
                                throw std::runtime_error(
                                    "load_obj_polymesh_file: Line number " +
                                    boost::lexical_cast<std::string>( lineNumber ) + "  of \"" +
                                    frantic::strings::to_string( objFile ) +
                                    "\" had an invalid face, because a vn index is out of bounds." );
                        }
                        normalsFaces.insert( normalsFaces.end(), tnface.begin(), tnface.end() );
                    }
                }
            }
        }
    }

    // handle warning the user of invalid faces
    if( customFacesErrorOnLine != -1 ) {
        FF_LOG( warning )
            << "load_obj_polymesh_file: Line number " << customFacesErrorOnLine << "  of \"" << objFile
            << "\" defines a face that inconsistent with previous face definitions. Texture coordinates and "
               "custom normal values will not be loaded."
            << std::endl;
        tverts.clear();
        normals.clear();
    }

    // handle warning case where there are not enough vertex channel values (for non-custom faces)
    if( tverts.size() && !hasCustomTVertFaces ) {
        if( tverts.size() != verts.size() ) {
            FF_LOG( warning )
                << "load_obj_polymesh_file: The number of TextureCoordinate \"vt\" values did not match the "
                   "number of vertex \"v\" values. Texture coordinates will not be loaded."
                << std::endl;
            tverts.clear();
        }
    }

    // handle warning case where there are not enough normal values (for non-custom faces)
    if( normals.size() && !hasCustomNormalFaces ) {
        if( normals.size() != verts.size() ) {
            FF_LOG( warning )
                << "load_obj_polymesh_file: The number of Normal \"vn\" values did not match the number of "
                   "vertex \"v\" values. Normal values will not be loaded."
                << std::endl;
            normals.clear();
        }
    }

    for( std::vector<vector3f>::const_iterator i = verts.begin(); i != verts.end(); ++i ) {
        builder.add_vertex( *i );
    }
    polymesh3_ptr mesh = builder.finalize();
    if( tverts.size() ) {
        frantic::graphics::raw_byte_buffer buffer;
        tverts.steal_raw_byte_buffer( buffer );
        std::vector<int>* pInFaceBuffer = hasCustomTVertFaces ? &tvertsFaces : 0;
        mesh->add_vertex_channel( _T("TextureCoord"), frantic::channels::data_type_float32, 3, buffer, pInFaceBuffer );
    }
    if( normals.size() ) {
        frantic::graphics::raw_byte_buffer buffer;
        normals.steal_raw_byte_buffer( buffer );
        std::vector<int>* pInFaceBuffer = hasCustomNormalFaces ? &normalsFaces : 0;
        mesh->add_vertex_channel( _T("Normal"), frantic::channels::data_type_float32, 3, buffer, pInFaceBuffer );
    }
    return mesh;
}

/**
 * Load a STL file from the specified path.
 *
 * The resulting mesh is guaranteed to have a 'verts' and 'Normal' channel.
 *
 * @param path The filesystem location to load the STL file from.
 * @return A new polymesh3 object with all the channels loaded.
 */
polymesh3_ptr load_stl_polymesh_file( const frantic::tstring& stlFile ) {
    std::ifstream fin( stlFile.c_str() );

    if( !fin )
        throw std::runtime_error( "load_stl_polymesh_file: Could not open file \"" +
                                  frantic::strings::to_string( stlFile ) + "\" for reading." );

    // check first 5 bytes to see if ascii
    std::vector<char> temp( 5 );
    bool ascii = false;
    bool hasColor = false;
    fin.read( &temp[0], 5 );
    std::string tempStr( temp.begin(), temp.end() );
    if( tempStr == "solid" )
        ascii = true;

    polymesh3_builder builder;
    std::vector<vector3f> verts;
    raw_byte_buffer_helper<vector3f> normals;
    raw_byte_buffer_helper<vector3f> colors;

    if( ascii ) {

        int lineNumber = 0;
        std::string line;
        std::vector<std::string> tokens;

        while( getline( fin, line ) ) {
            lineNumber++;
            tokens.clear();
            strings::split( line, tokens );

            // Process only if this line had any tokens
            if( tokens.size() > 0 ) {
                if( tokens[0] == "vertex" ) {
                    // Add a vertex
                    if( tokens.size() != 4 )
                        throw std::runtime_error(
                            "load_stl_polymesh_file: Line number " + boost::lexical_cast<std::string>( lineNumber ) +
                            "  of \"" + frantic::strings::to_string( stlFile ) + "\" had an invalid vertex." );
                    if( verts.size() >= normals.size() * 3 )
                        throw std::runtime_error(
                            "load_stl_polymesh_file: ascii file \"" + frantic::strings::to_string( stlFile ) +
                            "\" is not correctly formatted around Line number " +
                            boost::lexical_cast<std::string>( lineNumber ) + ", missing a normal" );
                    vector3f v( to_float( tokens[1] ), to_float( tokens[2] ), to_float( tokens[3] ) );
                    verts.push_back( v );

                    if( verts.size() % 3 == 0 ) {
                        std::vector<int> face;
                        face.push_back( (int)verts.size() - 3 );
                        face.push_back( (int)verts.size() - 2 );
                        face.push_back( (int)verts.size() - 1 );

                        builder.add_polygon( face );
                    }
                } else if( tokens[0] == "facet" ) {
                    if( tokens.size() != 5 )
                        throw std::runtime_error(
                            "load_stl_polymesh_file: Line number " + boost::lexical_cast<std::string>( lineNumber ) +
                            "  of \"" + frantic::strings::to_string( stlFile ) + "\" had an invalid normal." );
                    if( normals.size() * 3 != verts.size() )
                        throw std::runtime_error( "load_stl_polymesh_file: ascii file \"" +
                                                  frantic::strings::to_string( stlFile ) +
                                                  "\" is not correctly formatted around Line number " +
                                                  boost::lexical_cast<std::string>( lineNumber ) +
                                                  ", not the amount of vertices for the amount of normals" );
                    vector3f n( to_float( tokens[2] ), to_float( tokens[3] ), to_float( tokens[4] ) );
                    normals.push_back( n );
                } else if( tokens[0] == "outer" ) {
                    if( ( normals.size() - 1 ) * 3 != verts.size() )
                        throw std::runtime_error( "load_stl_polymesh_file: ascii file \"" +
                                                  frantic::strings::to_string( stlFile ) +
                                                  "\" is not correctly formatted aaround Line number " +
                                                  boost::lexical_cast<std::string>( lineNumber ) +
                                                  ", not the amount of vertices for the amount of normals" );
                } else if( tokens[0] == "endloop" ) {
                    if( verts.size() % 3 != 0 )
                        throw std::runtime_error(
                            "load_stl_polymesh_file: ascii file \"" + frantic::strings::to_string( stlFile ) +
                            "\" is not correctly formatted around Line number " +
                            boost::lexical_cast<std::string>( lineNumber ) + ", invalid amount of vertices" );
                } else if( tokens[0] == "endfacet" ) {
                    if( normals.size() * 3 != verts.size() )
                        throw std::runtime_error( "load_stl_polymesh_file: ascii file \"" +
                                                  frantic::strings::to_string( stlFile ) +
                                                  "\" is not correctly formatted around Line number " +
                                                  boost::lexical_cast<std::string>( lineNumber ) +
                                                  ", not the amount of vertices for the amount of normals" );
                }
            }
        }
    } else {
        fin.close(); // need to open with binary flag
        // open file and get size
        std::ifstream fin( stlFile.c_str(), std::ios::binary | std::ios::ate );

        if( !fin )
            throw std::runtime_error( "load_stl_polymesh_file: Could not open file \"" +
                                      frantic::strings::to_string( stlFile ) + "\" for reading." );
        // make sure it is larger than 84 bytes
        boost::uint32_t size = (boost::uint32_t)fin.tellg();
        boost::uint32_t facets;
        if( size < 84 )
            throw std::runtime_error( "load_stl_polymesh_file: binary file \"" +
                                      frantic::strings::to_string( stlFile ) +
                                      "\" is not correctly formatted, not large enough for header" );
        fin.seekg( 0, std::ios::beg );
        std::vector<char> header( 80 );
        fin.read( &header[0], 80 );
        // find if header contains the word COLOR= then set to magics
        std::string head( header.begin(), header.end() );
        bool isSolidView = true;
        hasColor = true;
        vector3f magicsDefaultColor;
        size_t pos = head.find( "COLOR=" );
        if( pos != std::string::npos ) {
            isSolidView = false; // it is magics
            char red = head[pos + 6];
            float r = (float)red / 255;
            char green = head[pos + 7];
            ;
            float g = (float)green / 255;
            char blue = head[pos + 8];
            ;
            float b = (float)blue / 255;
            magicsDefaultColor.set( r, g, b );
        }
        fin.read( (char*)&facets, 4 );
        // then check number of facets and make sure numfacets*50bytes are available from the file
        std::size_t facetDataSize = size - 84;
        if( facetDataSize != facets * 50 )
            throw std::runtime_error(
                "load_stl_polymesh_file: binary file \"" + frantic::strings::to_string( stlFile ) +
                "\" is not correctly formatted, not the correct number of bytes for amount of facets" );
        // loop through each facet
        for( boost::uint32_t facetIndex = 0; facetIndex < facets; ++facetIndex ) {
            // grab the first 12 bytes into the normal
            vector3f n( 0, 0, 0 );
            fin.read( (char*)&n.x, 4 );
            fin.read( (char*)&n.y, 4 );
            fin.read( (char*)&n.z, 4 );
            normals.push_back( n );
            // grab the 36 next bytes into vectors
            vector3f v1( 0, 0, 0 );
            fin.read( (char*)&v1.x, 4 );
            fin.read( (char*)&v1.y, 4 );
            fin.read( (char*)&v1.z, 4 );
            verts.push_back( v1 );

            vector3f v2( 0, 0, 0 );
            fin.read( (char*)&v2.x, 4 );
            fin.read( (char*)&v2.y, 4 );
            fin.read( (char*)&v2.z, 4 );
            verts.push_back( v2 );

            vector3f v3( 0, 0, 0 );
            fin.read( (char*)&v3.x, 4 );
            fin.read( (char*)&v3.y, 4 );
            fin.read( (char*)&v3.z, 4 );
            verts.push_back( v3 );

            std::vector<int> face;
            face.push_back( (int)verts.size() - 3 );
            face.push_back( (int)verts.size() - 2 );
            face.push_back( (int)verts.size() - 1 );

            builder.add_polygon( face );

            // get color from last 2 bytes or discard if it doesn't have any
            boost::uint16_t attribute;
            fin.read( (char*)&attribute, 2 );
            boost::uint8_t red;
            boost::uint8_t green;
            boost::uint8_t blue;
            bool valid;
            // STL are assumed little endian so use LSB bit numbering
            valid = ( boost::uint8_t )( attribute >> 15 ) & '\1';
            if( hasColor && isSolidView && valid ) {
                red = ( boost::uint8_t )( attribute >> 10 ) & 31;
                green = ( boost::uint8_t )( attribute >> 5 ) & 31;
                blue = ( boost::uint8_t )(attribute)&31;
                // need them as floats on scale of 0-1
                float r = boost::numeric_cast<float, boost::uint8_t>( red ) / 31;
                float g = boost::numeric_cast<float, boost::uint8_t>( green ) / 31;
                float b = boost::numeric_cast<float, boost::uint8_t>( blue ) / 31;
                vector3f color( r, g, b );
                // added 3 times since this is color for the 3 vertices
                colors.push_back( color );
                colors.push_back( color );
                colors.push_back( color );
            } else if( hasColor && !isSolidView && !valid ) {
                red = ( boost::uint8_t )(attribute)&31;
                green = ( boost::uint8_t )( attribute >> 5 ) & 31;
                blue = ( boost::uint8_t )( attribute >> 10 ) & 31;
                // need them as floats on scale of 0-1
                float r = boost::numeric_cast<float, boost::uint8_t>( red ) / 31;
                float g = boost::numeric_cast<float, boost::uint8_t>( green ) / 31;
                float b = boost::numeric_cast<float, boost::uint8_t>( blue ) / 31;
                vector3f color( r, g, b );
                // added 3 times since this is color for the 3 vertices
                colors.push_back( color );
                colors.push_back( color );
                colors.push_back( color );
            } else if( hasColor && !isSolidView && valid ) {
                // added 3 times since this is color for the 3 vertices
                colors.push_back( magicsDefaultColor );
                colors.push_back( magicsDefaultColor );
                colors.push_back( magicsDefaultColor );
            } else {
                hasColor = false;
            }
        }
    }

    for( std::vector<vector3f>::const_iterator i = verts.begin(); i != verts.end(); ++i ) {
        builder.add_vertex( *i );
    }
    polymesh3_ptr mesh = builder.finalize();

    // STL file requires normals channel so it is guaranteed to have it
    frantic::graphics::raw_byte_buffer normalBuf;
    normals.steal_raw_byte_buffer( normalBuf );
    mesh->add_face_channel( _T("Normal"), frantic::channels::data_type_float32, 3, normalBuf );

    // If the file has color write the channel
    if( hasColor ) {
        frantic::graphics::raw_byte_buffer colorBuf;
        colors.steal_raw_byte_buffer( colorBuf );
        mesh->add_vertex_channel( _T("Color"), frantic::channels::data_type_float32, 3, colorBuf );
    }

    return mesh;
}

polymesh3_ptr load_ply_polymesh_file( const frantic::tstring& plyFile ) {
    ply_reader reader( plyFile );
    return reader.read_polymesh3();
}

void write_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh ) {
    xmesh_metadata metadata;
    write_polymesh_file( path, polymesh, metadata );
}

void write_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh,
                          const xmesh_metadata& metadata ) {
    const frantic::tstring type = frantic::strings::to_lower( frantic::files::extension_from_path( path ) );
    if( type == _T(".obj") ) {
        write_obj_polymesh_file( path, polymesh );
    } else if( type == _T(".xmesh") ) {
        write_xmesh_polymesh_file( path, polymesh, metadata );
    } else {
        throw std::runtime_error( "write_polymesh_file: Didn't recognize the file format of the output mesh file \"" +
                                  frantic::strings::to_string( path ) + "\"" );
    }
}

void write_xmesh_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh,
                                const xmesh_metadata* metadata ) {
    xmesh_writer fileWriter( path );

    if( metadata ) {
        fileWriter.set_metadata( *metadata );
    }

    polymesh3::iterator it = polymesh->begin();
    polymesh3::iterator itEnd = polymesh->end();
    for( ; it != itEnd; ++it ) {
        if( it->second.is_vertex_channel() ) {
            // TODO: A general accessor would help out here

            std::size_t numVertices = it->second.get_data().size() / it->second.element_size();
            fileWriter.write_vertex_channel( it->first, it->second.get_data().ptr_at( 0 ), it->second.type(),
                                             it->second.arity(), numVertices );

            if( it->second.get_faces().size() > 0 ) {
                // There are custom faces. We need to translate them into the appropriate encoding.
                std::vector<int> faceBuffer( it->second.get_faces().size() );

                std::size_t cur = 0;
                std::size_t numPolygons = polymesh->m_pFaceEndOffsets->size();
                for( std::size_t i = 0; i < numPolygons; ++i ) {
                    std::size_t next = polymesh->m_pFaceEndOffsets->operator[]( i );
                    for( ; cur < next; ++cur )
                        faceBuffer[cur] = it->second.get_faces()[cur];
                    faceBuffer[cur - 1] = -faceBuffer[cur - 1] - 1;
                }

                fileWriter.write_vertex_channel_faces( it->first, &faceBuffer[0], faceBuffer.size() );
            }
        } else {
            std::size_t numFaces = it->second.get_data().size() / it->second.element_size();
            fileWriter.write_face_channel( it->first, it->second.get_data().ptr_at( 0 ), it->second.type(),
                                           it->second.arity(), numFaces );
        }
    }

    fileWriter.close();
}

void write_xmesh_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh ) {
    write_xmesh_polymesh_file( path, polymesh, (xmesh_metadata*)NULL );
}

void write_xmesh_polymesh_file( const frantic::tstring& path, const polymesh3_ptr& polymesh,
                                const xmesh_metadata& metadata ) {
    write_xmesh_polymesh_file( path, polymesh, &metadata );
}

void write_obj_polymesh_file( const frantic::tstring& destFile, const polymesh3_ptr& mesh ) {
    if( !mesh )
        throw std::runtime_error( "write_obj_polymesh_file: mesh is NULL" );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );

    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( destFile, meshInterface, progress );
}

/**
 * Writes a given polymesh3_ptr to the given filename
 * using the ascii stl format.
 * <p>
 * this function saves the mesh to an ASCII stl file
 * @param  mesh		the polymesh3_ptr to be saved out
 * @param  filename	the location/name of the destination stl file
 */
void write_ascii_stl_polymesh_file( const frantic::tstring& filename, const polymesh3_ptr& mesh ) {

    if( !mesh )
        throw std::runtime_error( "write_ascii_stl_polymesh_file: mesh is NULL" );

    std::ofstream fout( filename.c_str() );

    if( !fout )
        throw std::runtime_error( "write_ascii_stl_polymesh_file: Could not open file \"" +
                                  frantic::strings::to_string( filename ) + "\" for writing." );

    fout.precision( 6 );

    fout << "solid \n";

    polymesh3_const_vertex_accessor<vector3f> geomAcc = mesh->get_const_vertex_accessor<vector3f>( _T("verts") );

    polymesh3_const_cvt_face_accessor<vector3f> normalAcc;
    if( mesh->has_face_channel( _T("Normal") ) ) {
        normalAcc = mesh->get_const_cvt_face_accessor<vector3f>( _T("Normal") );
    }

    for( std::size_t i = 0, ie = geomAcc.face_count(); i != ie; ++i ) {
        // get 3 vertices based on the face
        polymesh3_const_face_range geomFace = geomAcc.get_face( i );
        int cornerCount = 0;
        vector3f faceVertices[3];
        for( polymesh3_const_face_iterator igeom = geomFace.first; igeom != geomFace.second; ++igeom ) {
            if( cornerCount >= 3 )
                throw std::runtime_error( "write_ascii_stl_polymesh_file: polymesh has an invalid face(more than 3 "
                                          "vertices) for STL file format" );
            faceVertices[cornerCount] = geomAcc.get_vertex( *igeom );
            ++cornerCount;
        }
        if( cornerCount < 3 )
            throw std::runtime_error( "write_ascii_stl_polymesh_file: polymesh has an invalid face(less than 3 "
                                      "vertices) for STL file format" );

        vector3f normal;

        if( normalAcc.is_valid() ) {
            normal = normalAcc.get_face( i ); // in for loop for faces already
        } else {
            normal = frantic::graphics::triangle_normal( faceVertices[0], faceVertices[1], faceVertices[2] );
        }

        fout << "\tfacet normal " << std::scientific << normal.x << " " << normal.y << " " << normal.z << "\n";
        fout << "\t\touter loop\n";
        for( int corner = 0; corner < 3; ++corner ) {
            fout << "\t\t\tvertex ";
            fout << std::scientific << faceVertices[corner].x << " " << faceVertices[corner].y << " "
                 << faceVertices[corner].z << "\n";
        }
        fout << "\t\tendloop\n";
        fout << "\tendfacet\n";
    }

    fout << "endsolid";
}

namespace {

class binary_stl_writer {
  public:
    binary_stl_writer( const frantic::tstring& filename )
        : m_file( frantic::files::tfopen( filename.c_str(), _T("wb") ) )
        , m_filename( filename ) {
        if( !m_file )
            throw std::runtime_error( "binary_stl_writer: Could not open file \"" +
                                      frantic::strings::to_string( filename ) + "\" for writing." );
    }

    void write( char* data, std::size_t size ) {
        if( !size ) {
            return;
        }

        assert( data );
        assert( m_file );

        const std::size_t writeCount = fwrite( data, size, 1, m_file );
        if( writeCount != 1 ) {
            throw std::runtime_error( "binary_stl_writer: Error writing to file \"" +
                                      frantic::strings::to_string( m_filename ) + "\"." );
        }
    }

  private:
    binary_stl_writer(); // not implemented

    frantic::files::file_ptr m_file;
    frantic::tstring m_filename;
};

} // anonymous namespace

/**
 * Writes a given polymesh3_ptr to the given filename
 * using the binary stl format.
 * <p>
 * this function saves the mesh to a binary stl file
 * @param  mesh		the polymesh3_ptr to be saved out
 * @param  filename	the location/name of the destination stl file
 * @param  isSolidView   determines if the color channel should be written in SolidView format or in Magics format
 */
void write_binary_stl_polymesh_file( const frantic::tstring& filename, const polymesh3_ptr& mesh, bool isSolidView ) {
    if( !mesh )
        throw std::runtime_error( "write_binary_stl_polymesh_file: mesh is NULL" );

    binary_stl_writer writer( filename );

    std::vector<char> header( 80 );

    polymesh3_const_cvt_vertex_accessor<vector3f> colorAcc;
    if( mesh->has_vertex_channel( _T("Color") ) ) {
        colorAcc = mesh->get_const_cvt_vertex_accessor<vector3f>( _T("Color") );

        // magics header string
        if( !isSolidView ) {
            std::string colorString =
                "COLOR=\0\0\0\0,MATERIAL=\0\0\0\0\0\0\0\0\0\0\0\0"; // whole object doesnt have color so just zero it
            header.assign( colorString.begin(), colorString.end() );
        }
    }
    writer.write( &header[0], 80 );

    polymesh3_const_vertex_accessor<vector3f> geomAcc = mesh->get_const_vertex_accessor<vector3f>( _T("verts") );

    polymesh3_const_cvt_face_accessor<vector3f> normalAcc;
    if( mesh->has_face_channel( _T("Normal") ) ) {
        normalAcc = mesh->get_const_cvt_face_accessor<vector3f>( _T("Normal") );
    }

    const std::size_t vertexCount = geomAcc.vertex_count();
    const std::size_t faceCount = geomAcc.face_count();

    const boost::int32_t facets = static_cast<boost::int32_t>( faceCount );

    writer.write( (char*)&facets, 4 );

    std::vector<vector3f> faceVertexColors;

    for( std::size_t i = 0, ie = geomAcc.face_count(); i != ie; ++i ) {
        // get 3 vertices based on the face
        polymesh3_const_face_range geomFace = geomAcc.get_face( i );
        size_t cornerCount = 0;
        vector3f faceVertices[3];
        faceVertexColors.clear();
        for( polymesh3_const_face_iterator igeom = geomFace.first; igeom != geomFace.second; ++igeom ) {
            if( cornerCount >= 3 )
                throw std::runtime_error( "write_binary_stl_polymesh_file: polymesh has an invalid face(more than 3 "
                                          "vertices) for STL file format" );
            faceVertices[cornerCount] = geomAcc.get_vertex( *igeom );
            if( colorAcc.is_valid() ) {
                if( !colorAcc.has_custom_faces() && colorAcc.vertex_count() == vertexCount ) {
                    faceVertexColors.push_back( colorAcc.get_vertex( *igeom ) );
                } else if( colorAcc.has_custom_faces() && colorAcc.face_count() == faceCount ) {
                    polymesh3_const_face_range colorFace = colorAcc.get_face( i );
                    for( polymesh3_const_face_iterator colorIterator = colorFace.first;
                         colorIterator != colorFace.second; ++colorIterator ) {
                        faceVertexColors.push_back( colorAcc.get_vertex( *colorIterator ) );
                    }
                }
            }
            ++cornerCount;
        }
        if( cornerCount < 3 )
            throw std::runtime_error( "write_binary_stl_polymesh_file: polymesh has an invalid face(less than 3 "
                                      "vertices) for STL file format" );

        vector3f normal;

        if( normalAcc.is_valid() ) {
            normal = normalAcc.get_face( i ); // in for loop for faces already
        } else {
            normal = frantic::graphics::triangle_normal( faceVertices[0], faceVertices[1], faceVertices[2] );
        }

        writer.write( (char*)&normal.x, 4 );
        writer.write( (char*)&normal.y, 4 );
        writer.write( (char*)&normal.z, 4 );

        // write the 3 vertices
        for( int corner = 0; corner < 3; ++corner ) {
            writer.write( (char*)&faceVertices[corner].x, 4 );
            writer.write( (char*)&faceVertices[corner].y, 4 );
            writer.write( (char*)&faceVertices[corner].z, 4 );
        }

        boost::uint16_t attribute = 0;
        // check if want to save colors
        // if yes check whether solidview or magics
        // if no color channel set not valid or per-object color
        if( colorAcc.is_valid() && faceVertexColors.size() == cornerCount ) {
            // average out vertices colors
            vector3f faceColor;
            for( size_t corner = 0; corner < cornerCount; ++corner ) {
                faceColor += faceVertexColors[corner];
            }
            faceColor /= static_cast<float>( cornerCount );
            boost::uint8_t red = (boost::uint8_t)math::round( faceColor.x * 31 );
            boost::uint8_t green = (boost::uint8_t)math::round( faceColor.y * 31 );
            boost::uint8_t blue = (boost::uint8_t)math::round( faceColor.z * 31 );
            // STL is assumed little endian so it is best to assume lsb bit numbering
            // Solidview
            if( isSolidView ) {
                attribute = 1;
                attribute <<= 5; // left shift the bits of attribute
                attribute |= red;
                attribute <<= 5;    // left shift the bits of attribute
                attribute |= green; // or attribute with green
                attribute <<= 5;
                attribute |= blue;
            } else {
                // Magics
                attribute = blue;
                attribute <<= 5;    // left shift the bits of attribute
                attribute |= green; // or attribute with green
                attribute <<= 5;
                attribute |= red;
            }
        } else {
            attribute = 0; // no color channel just set it as 0
        }

        writer.write( (char*)&attribute, 2 );
    }
}

/**
 * Writes a given polymesh3_ptr to the given filename using the ply format.
 * <p>
 * this function saves the mesh to an ply file
 * @param  mesh		the polymesh3_ptr to be saved out
 * @param  filename	the location/name of the destination ply file
 */
void write_ply_polymesh_file( const frantic::tstring& filename, const polymesh3_ptr& mesh ) {
    if( !mesh )
        throw std::runtime_error( "write_ply_polymesh_file: mesh is NULL" );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );

    frantic::logging::null_progress_logger progress;

    write_ply_mesh_file( filename, meshInterface, progress );
}

} // namespace geometry
} // namespace frantic
