// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/algorithm/string/trim.hpp>

#include <frantic/files/paths.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/trimesh3_file_io.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>
#include <frantic/geometry/xmesh_sequence_saver.hpp>
#include <frantic/locale/locale.hpp>

#include <frantic/tinyxml/frantic_tinyxml_graphics.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::graphics;

namespace frantic {
namespace geometry {

namespace {

float to_float( const std::string& s ) { return static_cast<float>( frantic::locale::strtod_c( s.c_str(), 0 ) ); }

} // anonymous namespace

// This loads an .obj file into a trimesh3
void load_obj_mesh_file( const frantic::tstring& objFile, trimesh3& mesh ) {
    mesh.clear();

    ifstream fin( objFile.c_str() );

    if( !fin )
        throw std::runtime_error( "load_obj_mesh_file: Could not open file \"" +
                                  frantic::strings::to_string( objFile ) + "\" for reading." );

    // These will become the texture coord and normal channel accessors if need be
    trimesh3_vertex_channel_accessor<vector3f> tverts, normals;
    bool hasCustomTVertFaces = false;
    bool hasCustomNormalFaces = false;

    // this will be set to non-negative if face definitions are inconsistent (ie. there are four ways to define face
    // data. they are: "face", "face/tvface", "face/tvface/normal", "face//normal")
    int customFacesErrorOnLine = -1;

    string line;
    int lineNumber = 0;
    vector<string> tokens;
    vector<string> faceTokens;
    while( getline( fin, line ) ) {
        boost::algorithm::trim_right_if( line, boost::is_any_of( "\r" ) );

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
                            "load_obj_mesh_file: Line number " + boost::lexical_cast<string>( lineNumber ) + "  of \"" +
                            frantic::strings::to_string( objFile ) + "\" had an invalid vertex." );
                    vector3f v( to_float( tokens[1] ), to_float( tokens[2] ), to_float( tokens[3] ) );
                    mesh.add_vertex( v );
                } else if( tokens[0] == "vt" ) {
                    // Add a texture vertex
                    frantic::graphics::vector3f vt;
                    if( tokens.size() == 3 ) {
                        vt.set( to_float( tokens[1] ), to_float( tokens[2] ), 0 );
                    } else if( tokens.size() == 4 ) {
                        vt.set( to_float( tokens[1] ), to_float( tokens[2] ), to_float( tokens[3] ) );
                    } else {
                        throw std::runtime_error(
                            "load_obj_mesh_file: Line number " + boost::lexical_cast<string>( lineNumber ) + "  of \"" +
                            frantic::strings::to_string( objFile ) + "\" had an invalid texture vertex." );
                    }
                    if( !tverts ) {
                        mesh.add_vertex_channel<vector3f>( _T("TextureCoord"), 0 );
                        tverts = mesh.get_vertex_channel_accessor<vector3f>( _T("TextureCoord") );
                    }
                    tverts.add_vertex( vt );
                } else if( tokens[0] == "vn" ) {
                    // Add a normal
                    if( tokens.size() != 4 )
                        throw std::runtime_error(
                            "load_obj_mesh_file: Line number " + boost::lexical_cast<string>( lineNumber ) + "  of \"" +
                            frantic::strings::to_string( objFile ) + "\" had an invalid normal." );
                    vector3f vn( to_float( tokens[1] ), to_float( tokens[2] ), to_float( tokens[3] ) );
                    if( !normals ) {
                        mesh.add_vertex_channel<vector3f>( _T("Normal"), 0 );
                        normals = mesh.get_vertex_channel_accessor<vector3f>( _T("Normal") );
                    }
                    normals.add_vertex( vn );
                } else if( tokens[0] == "f" ) {
                    // Add a face.  Vertices are 1-based in obj files, so we subtract 1 to compensate.
                    if( tokens.size() < 4 )
                        throw std::runtime_error( "load_obj_mesh_file: Line number " +
                                                  boost::lexical_cast<string>( lineNumber ) + "  of \"" +
                                                  frantic::strings::to_string( objFile ) + "\" had an invalid face." );
                    if( mesh.face_count() == 0 ) {
                        // If this is the first face we're adding, it will determine whether vt or vn have custom faces
                        // or share their faces with v
                        faceTokens.clear();
                        strings::split( tokens[1], faceTokens, '/' );
                        if( faceTokens.size() > 1 && faceTokens[1] != "" ) {
                            if( tverts ) {
                                hasCustomTVertFaces = true;
                                mesh.set_vertex_channel_custom_faces( _T("TextureCoord"), true );
                                // Changing the custom faces flag invalidates the existing accessor, so we need to get a
                                // new one
                                tverts = mesh.get_vertex_channel_accessor<vector3f>( _T("TextureCoord") );
                            } else {
                                throw std::runtime_error( "load_obj_mesh_file: Line number " +
                                                          boost::lexical_cast<string>( lineNumber ) + "  of \"" +
                                                          frantic::strings::to_string( objFile ) +
                                                          "\" had contains a index for the texture coordinate data, "
                                                          "but no \"vt\" data has been defined." );
                            }
                        }
                        if( faceTokens.size() > 2 && faceTokens[2] != "" ) {
                            if( normals ) {
                                hasCustomNormalFaces = true;
                                mesh.set_vertex_channel_custom_faces( _T("Normal"), true );
                                // Changing the custom faces flag invalidates the existing accessor, so we need to get a
                                // new one
                                normals = mesh.get_vertex_channel_accessor<vector3f>( _T("Normal") );
                            } else {
                                throw std::runtime_error( "load_obj_mesh_file: Line number " +
                                                          boost::lexical_cast<string>( lineNumber ) + "  of \"" +
                                                          frantic::strings::to_string( objFile ) +
                                                          "\" had contains a index for the normal data, but no \"vn\" "
                                                          "data has been defined." );
                            }
                        }
                    }
                    vector<int> face, tvface, tnface;
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
                            face[i] += ( 1 + static_cast<int>( mesh.vertex_count() ) );

                        if( face[i] < 0 || face[i] >= (int)mesh.vertex_count() )
                            throw std::runtime_error( "load_obj_mesh_file: Line number " +
                                                      boost::lexical_cast<string>( lineNumber ) + "  of \"" +
                                                      frantic::strings::to_string( objFile ) +
                                                      "\" had an invalid face, because a v index is out of bounds." );
                    }
                    mesh.add_face( face );

                    // add the vertex coord face to the mesh
                    if( hasCustomTVertFaces ) {
                        for( unsigned i = 0; i < tvface.size(); ++i ) {
                            if( tvface[i] < -1 )
                                tvface[i] += ( 1 + static_cast<int>( tverts.size() ) );

                            if( tvface[i] < 0 || tvface[i] >= (int)tverts.size() )
                                throw std::runtime_error(
                                    "load_obj_mesh_file: Line number " + boost::lexical_cast<string>( lineNumber ) +
                                    "  of \"" + frantic::strings::to_string( objFile ) +
                                    "\" had an invalid face, because a vt index is out of bounds." );
                        }
                        tverts.add_face( tvface );
                    }

                    // add the normal face to the mesh
                    if( hasCustomNormalFaces ) {
                        for( unsigned i = 0; i < tnface.size(); ++i ) {
                            if( tnface[i] < -1 )
                                tnface[i] += ( 1 + static_cast<int>( normals.size() ) );

                            if( tnface[i] < 0 || tnface[i] >= (int)normals.size() )
                                throw std::runtime_error(
                                    "load_obj_mesh_file: Line number " + boost::lexical_cast<string>( lineNumber ) +
                                    "  of \"" + frantic::strings::to_string( objFile ) +
                                    "\" had an invalid face, because a vn index is out of bounds." );
                        }
                        normals.add_face( tnface );
                    }
                }
            }
        }
    }

    // handle warning the user of invalid faces
    if( customFacesErrorOnLine != -1 ) {
        FF_LOG( warning )
            << "load_obj_mesh_file: Line number " << customFacesErrorOnLine << "  of \"" << objFile
            << "\" defines a face that inconsistent with previous face definitions. Texture coordinates and "
               "custom normal values will not be loaded."
            << std::endl;
        if( tverts )
            mesh.erase_vertex_channel( _T("TextureCoord") ); // deleting the channels because they contain invalid data.
        if( normals )
            mesh.erase_vertex_channel( _T("Normal") );
    }

    // handle warning case where there are not enough vertex channel values (for non-custom faces)
    if( tverts && !hasCustomTVertFaces ) {
        if( tverts.size() != mesh.vertex_count() ) {
            FF_LOG( warning )
                << "load_obj_mesh_file: The number of TextureCoordinate \"vt\" values did not match the number "
                   "of vertex \"v\" values. Texture coordinates will not be loaded."
                << std::endl;
            mesh.erase_vertex_channel( _T("TextureCoord") );
        }
    }

    // handle warning case where there are not enough normal values (for non-custom faces)
    if( normals && !hasCustomNormalFaces ) {
        if( normals.size() != mesh.vertex_count() ) {
            FF_LOG( warning )
                << "load_obj_mesh_file: The number of Normal \"vn\" values did not match the number of vertex "
                   "\"v\" values. Normal values will not be loaded."
                << std::endl;
            mesh.erase_vertex_channel( _T("Normal") );
        }
    }
}

void write_obj_mesh_file( const frantic::tstring& destFile, const trimesh3& mesh ) {
    mesh_interface_ptr meshInterface( trimesh3_interface::create_instance( &mesh ).release() );

    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( destFile, meshInterface, progress );
}

// Detect extension/format
void load_mesh_file( const frantic::tstring& srcFile, trimesh3& mesh ) {
    frantic::tstring type = frantic::strings::to_lower( frantic::files::extension_from_path( srcFile ) );
    if( type == _T(".obj") ) {
        load_obj_mesh_file( srcFile, mesh );
    } else if( type == _T(".xmesh") ) {
        load_xmesh_file( srcFile, mesh );
    } else {
        throw std::runtime_error( "load_mesh_file: Didn't recognize the file format of the input mesh file \"" +
                                  frantic::strings::to_string( srcFile ) + "\"" );
    }

    stringstream sout;
    sout << "load_mesh_file: Input mesh \"" << frantic::strings::to_string( srcFile )
         << "\" failed consistency check:\n";
    if( !mesh.check_consistency( sout ) ) {
        throw runtime_error( sout.str() );
    }
}

void write_mesh_file( const frantic::tstring& destFile, const trimesh3& mesh ) {
    frantic::tstring type = frantic::strings::to_lower( frantic::files::extension_from_path( destFile ) );
    if( type == _T(".obj") ) {
        write_obj_mesh_file( destFile, mesh );
    } else if( type == _T(".xmesh") ) {
        // Save the xmesh file!
        xmesh_sequence_saver xss;
        xss.write_xmesh( mesh, destFile );
    } else {
        throw std::runtime_error( "write_mesh_file: Didn't recognize the file format of the output mesh file \"" +
                                  frantic::strings::to_string( destFile ) + "\"" );
    }
}

bool try_load_boundbox_metadata( const frantic::tstring& srcFile, boundbox3f& outBoundbox, const bool showWarnings ) {
    try {
        outBoundbox.set_to_empty();

        vector3f minCorner;
        vector3f maxCorner;

        frantic::tstring type = frantic::strings::to_lower( frantic::files::extension_from_path( srcFile ) );
        if( type == _T(".xmesh") ) {
            tinyxml2::XMLDocument doc( frantic::strings::to_string( srcFile ).c_str() );
            frantic::files::file_ptr f( frantic::files::tfopen( srcFile.c_str(), _T("rb") ) );
            tinyxml2::XMLError result = doc.LoadFile( f );
            if( result != tinyxml2::XMLError::XML_SUCCESS ) {
                FF_LOG( warning ) << "try_load_boundbox_metadata Warning: Unable to open the requested XMesh file \""
                                  << srcFile << "\".\n";
                return false;
            }

            tinyxml2::XMLElement* pXMesh = doc.FirstChildElement( "xmesh" );

            if( pXMesh == 0 ) {
                FF_LOG( warning ) << "try_load_boundbox_metadata Warning: The XML in the XMesh file \"" << srcFile
                                  << "\" didn't contain an <xmesh> node.\n";
                return false;
            }

            tinyxml2::XMLHandle docHandle( &doc );
            tinyxml2::XMLHandle boundboxHandle = docHandle.FirstChildElement( "xmesh" ).FirstChildElement( "boundbox" );

            if( boundboxHandle.ToElement() ) {
                minCorner = frantic::tinyxml::get_vector3f( "minimum", boundboxHandle );
                maxCorner = frantic::tinyxml::get_vector3f( "maximum", boundboxHandle );
                outBoundbox = frantic::graphics::boundbox3f( minCorner, maxCorner );
                return true;
            }
        } else {
            if( showWarnings ) {
                FF_LOG( warning )
                    << "try_load_boundbox_metadata Warning: Did not recognize the file format of the input mesh file \""
                    << srcFile << "\".  Only .xmesh files are supported.\n";
            }
        }
    } catch( const std::exception& e ) {
        if( showWarnings ) {
            FF_LOG( warning ) << "try_load_boundbox_metadata Warning: " << e.what() << "\n";
        }
    }
    return false;
}

} // namespace geometry
} // namespace frantic
