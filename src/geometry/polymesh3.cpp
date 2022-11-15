// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/foreach.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <frantic/files/files.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/mesh_merging.hpp>
#include <frantic/geometry/ply_reader.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/xmesh_reader.hpp>
#include <frantic/geometry/xmesh_writer.hpp>
#include <frantic/graphics/plane3f.hpp>
#include <frantic/locale/locale.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/exception_stream.hpp>

#include <xxhash.h>

using frantic::graphics::raw_byte_buffer;
using frantic::graphics::vector3f;

namespace frantic {
namespace geometry {

// polymesh3_channel_data

polymesh3_channel_data::polymesh3_channel_data()
    : m_arity( 0 )
    , m_elementSize( 0 )
    , m_type( frantic::channels::data_type_invalid ) {}

polymesh3_channel_data::polymesh3_channel_data( frantic::graphics::raw_byte_buffer& data,
                                                frantic::channels::data_type_t type, std::size_t arity )
    : m_arity( arity )
    , m_elementSize( frantic::channels::sizeof_channel_data_type( type ) * arity )
    , m_type( type ) {
    if( data.size() % m_elementSize != 0 ) {
        throw frantic::exception_stream() << "polymesh3_channel_data() Data size (" << data.size() << ") "
                                          << "is not a multiple of the element size (" << m_elementSize << ")";
    }

    m_pData = boost::make_shared<raw_byte_buffer>();
    m_pData->swap( data );
}

std::size_t polymesh3_channel_data::arity() const { return m_arity; }

std::size_t polymesh3_channel_data::element_size() const { return m_elementSize; }

frantic::channels::data_type_t polymesh3_channel_data::type() const { return m_type; }

const frantic::graphics::raw_byte_buffer& polymesh3_channel_data::get() const {
    assert_data();
    return *m_pData;
}

frantic::graphics::raw_byte_buffer& polymesh3_channel_data::get_writable() {
    assert_data();
    if( is_shared() ) {
        // m_pData = boost::make_shared<raw_byte_buffer>( *m_pData );
        throw std::runtime_error( "polymesh3_channel_data::get_writable() Cannot get write access to shared data" );
    }
    return *m_pData;
}

bool polymesh3_channel_data::is_valid() const { return m_pData != 0; }

bool polymesh3_channel_data::is_shared() const {
    assert_data();
    return !m_pData.unique();
}

void polymesh3_channel_data::assert_data() const {
    if( !m_pData ) {
        throw std::runtime_error( "polymesh3_channel_data::assert_data() Data is NULL" );
    }
}

// polymesh3_channel_faces

polymesh3_channel_faces::polymesh3_channel_faces() { m_pFaces = boost::make_shared<std::vector<int>>(); }

polymesh3_channel_faces::polymesh3_channel_faces( std::vector<int>& faces ) {
    m_pFaces = boost::make_shared<std::vector<int>>();
    m_pFaces->swap( faces );
}

const std::vector<int>& polymesh3_channel_faces::get() const {
    assert_data();
    return *m_pFaces;
}

std::vector<int>& polymesh3_channel_faces::get_writable() {
    assert_data();
    if( is_shared() ) {
        // m_pFaces = boost::make_shared<std::vector<int> >( *m_pFaces );
        throw std::runtime_error( "polymesh3_channel_faces::get_writable() Cannot get write access to shared data" );
    }
    return *m_pFaces;
}

bool polymesh3_channel_faces::is_valid() const { return m_pFaces != 0; }

bool polymesh3_channel_faces::is_shared() const {
    assert_data();
    return !m_pFaces.unique();
}

void polymesh3_channel_faces::assert_data() const {
    if( !m_pFaces ) {
        throw std::runtime_error( "polymesh3_channel_faces::assert_data() Faces are NULL" );
    }
}

// polymesh3_channel

const frantic::graphics::raw_byte_buffer& polymesh3_channel::get_data() const { return data.get(); }

frantic::graphics::raw_byte_buffer& polymesh3_channel::get_writable_data() { return data.get_writable(); }

const std::vector<int>& polymesh3_channel::get_faces() const { return faces.get(); }

std::vector<int>& polymesh3_channel::get_writable_faces() { return faces.get_writable(); }

bool polymesh3_channel::is_vertex_channel() const { return m_isVertexChannel; }

std::size_t polymesh3_channel::arity() const { return data.arity(); }

std::size_t polymesh3_channel::element_size() const { return data.element_size(); }

frantic::channels::data_type_t polymesh3_channel::type() const { return data.type(); }

// polymesh3

polymesh3::polymesh3( polymesh3_channel_data vertData, polymesh3_channel_faces geomPolyIndices,
                      polymesh3_channel_faces geomPolyEndIndices ) {
    init( vertData, geomPolyIndices, geomPolyEndIndices );
}

polymesh3::polymesh3( frantic::graphics::raw_byte_buffer& vertData, std::vector<int>& faceIndices,
                      std::vector<int>& faceEndIndices ) {
    polymesh3_channel_data holdVertData( vertData, frantic::channels::data_type_float32, 3 );
    polymesh3_channel_faces holdGeomPolyIndices( faceIndices );
    polymesh3_channel_faces holdGeomPolyEndIndices( faceEndIndices );

    init( holdVertData, holdGeomPolyIndices, holdGeomPolyEndIndices );
}

void polymesh3::acquire() const { ++m_refCount; }

void polymesh3::release() const {
    if( --m_refCount <= 0 )
        delete this;
}

void polymesh3::init( polymesh3_channel_data vertData, polymesh3_channel_faces geomPolyIndices,
                      polymesh3_channel_faces geomPolyEndIndices ) {
    m_refCount = 0;

    if( !vertData.is_valid() ) {
        throw std::runtime_error( "polymesh3::init() The vertex data is uninitialized" );
    }

    if( vertData.type() != frantic::channels::data_type_float32 ) {
        throw std::runtime_error( "polymesh3::init() The vertex data has an unexpected type (" +
                                  boost::lexical_cast<std::string>( vertData.type() ) + ")" );
    }

    if( vertData.arity() != 3 ) {
        throw std::runtime_error( "polymesh3::init() The vertex data has an unexpected arity (" +
                                  boost::lexical_cast<std::string>( vertData.arity() ) + ")" );
    }

    if( !geomPolyIndices.is_valid() ) {
        throw std::runtime_error( "polymesh3::init() The polygon indices are uninitialized" );
    }

    if( !geomPolyEndIndices.is_valid() ) {
        throw std::runtime_error( "polymesh3::init() The polygon end indices are uninitialized" );
    }

    const std::size_t faceVertexCountFromIndices = geomPolyIndices.get().size();
    const std::size_t faceVertexCountFromEndIndices =
        geomPolyEndIndices.get().size() > 0 ? geomPolyEndIndices.get().back() : 0;
    if( faceVertexCountFromIndices != faceVertexCountFromEndIndices ) {
        throw frantic::exception_stream() << "polymesh3::init() Mismatch between number of face vertex indices from "
                                          << "face indices (" << faceVertexCountFromIndices << ") "
                                          << "vs face end indices (" << faceVertexCountFromEndIndices << ")";
    }

    polymesh3_channel& geomChannel = m_channels[_T("verts")];
    geomChannel.m_isVertexChannel = true;
    geomChannel.data = vertData;
    geomChannel.faces = geomPolyIndices;

    m_faceEndOffsetsHoldRef = geomPolyEndIndices;
    m_pFaceEndOffsets = &geomPolyEndIndices.get();
    m_pGeomChannel = &geomChannel;
}

void polymesh3::add_vertex_channel( const frantic::tstring& channel, polymesh3_channel_data& inVertexBuffer,
                                    polymesh3_channel_faces* pInFaceBuffer ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it != m_channels.end() )
        throw std::runtime_error( "polymesh3::add_vertex_channel() A channel by the name: \"" +
                                  frantic::strings::to_string( channel ) + "\" already exists" );
    if( !inVertexBuffer.is_valid() )
        throw std::runtime_error( "polymesh3::add_vertex_channel() The channel data is uninitialzed for channel \"" +
                                  frantic::strings::to_string( channel ) + "\"." );

    if( pInFaceBuffer ) {
        if( !pInFaceBuffer->is_valid() )
            throw std::runtime_error(
                "polymesh3::add_vertex_channel() The channel faces are uninitialzed for channel \"" +
                frantic::strings::to_string( channel ) + "\"." );
        if( pInFaceBuffer->get().size() != face_vertex_count() )
            throw std::runtime_error( "polymesh3::add_vertex_channel() The vertex channel \"" +
                                      frantic::strings::to_string( channel ) +
                                      "\" did not match the custom face mapping fo the 'verts' channel" );
    } else {
        if( inVertexBuffer.get().size() / inVertexBuffer.element_size() != vertex_count() )
            throw std::runtime_error(
                "polymesh3::add_vertex_channel() The vertex channel \"" + frantic::strings::to_string( channel ) +
                "\" did not match the length of the 'verts' channel, and also did not have a custom face mapping" );
    }

    polymesh3_channel& ch = m_channels[channel];
    ch.m_isVertexChannel = true;
    ch.data = inVertexBuffer;

    if( pInFaceBuffer ) {
        ch.faces = *pInFaceBuffer;
    }
}

void polymesh3::add_face_channel( const frantic::tstring& channel, polymesh3_channel_data& inFaceBuffer ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it != m_channels.end() )
        throw std::runtime_error( "polymesh3::add_face_channel() A channel by the name: \"" +
                                  frantic::strings::to_string( channel ) + "\" already exists" );
    if( !inFaceBuffer.is_valid() )
        throw std::runtime_error( "polymesh3::add_face_channel() The channel data is uninitialzed for channel \"" +
                                  frantic::strings::to_string( channel ) + "\"." );

    if( inFaceBuffer.get().size() / inFaceBuffer.element_size() != face_count() )
        throw std::runtime_error(
            "polymesh3::add_face_channel() The face channel \"" + frantic::strings::to_string( channel ) +
            "\" did not have one entry per-polygon.  The channel has " +
            boost::lexical_cast<std::string>( inFaceBuffer.get().size() / inFaceBuffer.element_size() ) +
            " entries, while the mesh has " + boost::lexical_cast<std::string>( face_count() ) + " polygons." );

    polymesh3_channel& ch = m_channels[channel];
    ch.m_isVertexChannel = false;
    ch.data = inFaceBuffer;
}

std::size_t polymesh3::face_vertex_count() const {
    assert( m_pGeomChannel );

    return m_pGeomChannel->get_faces().size();
}

void polymesh3::get_vertex_channel_names( std::vector<frantic::tstring>& outNames ) const {
    outNames.clear();
    for( iterator it = m_channels.begin(); it != m_channels.end(); ++it ) {
        if( it->second.is_vertex_channel() && it->first != _T( "verts" ) )
            outNames.push_back( it->first );
    }
}

void polymesh3::get_face_channel_names( std::vector<frantic::tstring>& outNames ) const {
    outNames.clear();
    for( iterator it = m_channels.begin(); it != m_channels.end(); ++it ) {
        if( !it->second.is_vertex_channel() )
            outNames.push_back( it->first );
    }
}

void polymesh3::get_face_end_offsets( std::vector<int>& outOffsets ) const { outOffsets = *m_pFaceEndOffsets; }

polymesh3_channel& polymesh3::get_vertex_channel( const frantic::tstring& channel ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it == m_channels.end() )
        throw std::runtime_error( "polymesh3::get_vertex_channel() There is no channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\"" );
    if( !it->second.is_vertex_channel() )
        throw std::runtime_error( "polymesh3::get_vertex_channel() The channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\" is a face channel." );
    return it->second;
}

const polymesh3_channel& polymesh3::get_const_vertex_channel( const frantic::tstring& channel ) const {
    std::map<frantic::tstring, polymesh3_channel>::const_iterator it = m_channels.find( channel );
    if( it == m_channels.end() )
        throw std::runtime_error( "polymesh3::get_const_vertex_channel() There is no channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\"" );
    if( !it->second.is_vertex_channel() )
        throw std::runtime_error( "polymesh3::get_const_vertex_channel() The channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\" is a face channel." );
    return it->second;
}

polymesh3_channel& polymesh3::get_face_channel( const frantic::tstring& channel ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it == m_channels.end() )
        throw std::runtime_error( "polymesh3::get_face_channel() There is no channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\"" );
    if( it->second.is_vertex_channel() )
        throw std::runtime_error( "polymesh3::get_face_channel() The channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\" is a vertex channel." );
    return it->second;
}

const polymesh3_channel& polymesh3::get_const_face_channel( const frantic::tstring& channel ) const {
    std::map<frantic::tstring, polymesh3_channel>::const_iterator it = m_channels.find( channel );
    if( it == m_channels.end() )
        throw std::runtime_error( "polymesh3::get_const_face_channel() There is no channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\"" );
    if( it->second.is_vertex_channel() )
        throw std::runtime_error( "polymesh3::get_const_face_channel() The channel named: \"" +
                                  frantic::strings::to_string( channel ) + "\" is a vertex channel." );
    return it->second;
}

bool polymesh3::has_vertex_channel( const frantic::tstring& name ) const {
    std::map<frantic::tstring, polymesh3_channel>::const_iterator i = m_channels.find( name );
    return ( i != m_channels.end() ) && i->second.is_vertex_channel();
}

bool polymesh3::has_face_channel( const frantic::tstring& name ) const {
    std::map<frantic::tstring, polymesh3_channel>::const_iterator i = m_channels.find( name );
    return ( i != m_channels.end() ) && ( !i->second.is_vertex_channel() );
}

bool polymesh3::is_triangle_mesh() const {
    polymesh3_const_vertex_accessor<frantic::graphics::vector3f> acc =
        get_const_vertex_accessor<frantic::graphics::vector3f>( _T("verts") );
    std::size_t numFaces = acc.face_count();
    for( std::size_t i = 0; i < numFaces; ++i ) {
        if( acc.get_face_degree( i ) != 3 )
            return false;
    }
    return true;
}

void polymesh3::add_empty_vertex_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                          std::size_t arity ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it != m_channels.end() )
        throw std::runtime_error( "polymesh3::add_empty_vertex_channel() A channel by the name: \"" +
                                  frantic::strings::to_string( channel ) + "\" already exists" );

    const std::size_t elementSize = arity * frantic::channels::sizeof_channel_data_type( type );

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( vertex_count() * elementSize );

    add_vertex_channel( channel, type, arity, buffer );
}

void polymesh3::add_empty_vertex_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                          std::size_t arity, std::size_t numVertices ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it != m_channels.end() )
        throw std::runtime_error( "polymesh3::add_empty_vertex_channel() A channel by the name: \"" +
                                  frantic::strings::to_string( channel ) + "\" already exists" );

    const std::size_t elementSize = arity * frantic::channels::sizeof_channel_data_type( type );

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( numVertices * elementSize );

    polymesh3_channel& ch = m_channels[channel];
    ch.m_isVertexChannel = true;
    ch.data = polymesh3_channel_data( buffer, type, arity );
    ch.get_writable_faces().resize( face_vertex_count() );
}

void polymesh3::add_vertex_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                    std::size_t arity, frantic::graphics::raw_byte_buffer& inVertexBuffer,
                                    std::vector<int>* pInFaceBuffer ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it != m_channels.end() )
        throw std::runtime_error( "polymesh3::add_vertex_channel() A channel by the name: \"" +
                                  frantic::strings::to_string( channel ) + "\" already exists" );

    const std::size_t elementSize = arity * frantic::channels::sizeof_channel_data_type( type );

    if( pInFaceBuffer ) {
        if( pInFaceBuffer->size() != face_vertex_count() )
            throw std::runtime_error( "polymesh3::add_vertex_channel() The vertex channel \"" +
                                      frantic::strings::to_string( channel ) +
                                      "\" did not match the custom face mapping fo the 'verts' channel" );
    } else {
        if( inVertexBuffer.size() != vertex_count() * elementSize )
            throw std::runtime_error(
                "polymesh3::add_vertex_channel() The vertex channel \"" + frantic::strings::to_string( channel ) +
                "\" did not match the length of the 'verts' channel, and also did not have a custom face mapping" );
    }

    polymesh3_channel_data holdData( inVertexBuffer, type, arity );
    if( pInFaceBuffer ) {
        polymesh3_channel_faces holdFaces( *pInFaceBuffer );
        add_vertex_channel( channel, holdData, &holdFaces );
    } else {
        add_vertex_channel( channel, holdData );
    }
}

void polymesh3::add_empty_face_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                        std::size_t arity ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it != m_channels.end() )
        throw std::runtime_error( "polymesh3::add_empty_face_channel() A channel by the name: \"" +
                                  frantic::strings::to_string( channel ) + "\" already exists" );

    const std::size_t elementSize = arity * frantic::channels::sizeof_channel_data_type( type );

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( face_count() * elementSize );

    add_face_channel( channel, type, arity, buffer );
}

void polymesh3::add_face_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                  std::size_t arity, frantic::graphics::raw_byte_buffer& inFaceBuffer ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );
    if( it != m_channels.end() )
        throw std::runtime_error( "polymesh3::add_face_channel() A channel by the name: \"" +
                                  frantic::strings::to_string( channel ) + "\" already exists" );

    const std::size_t elementSize = arity * frantic::channels::sizeof_channel_data_type( type );

    if( inFaceBuffer.size() / elementSize != face_count() )
        throw std::runtime_error(
            "polymesh3::add_face_channel() The face channel \"" + frantic::strings::to_string( channel ) +
            "\" did not have one entry per-polygon.  The channel has " +
            boost::lexical_cast<std::string>( inFaceBuffer.size() / elementSize ) + " entries, while the mesh has " +
            boost::lexical_cast<std::string>( face_count() ) + " polygons." );

    polymesh3_channel_data holdData( inFaceBuffer, type, arity );
    add_face_channel( channel, holdData );
}

void polymesh3::erase_vertex_channel( const frantic::tstring& channel ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );

    if( it == m_channels.end() ) {
        throw std::runtime_error( "polymesh3::erase_vertex_channel -- No such channel \"" +
                                  frantic::strings::to_string( channel ) + "\" exists." );
    } else if( !it->second.is_vertex_channel() ) {
        throw std::runtime_error( "polymesh3::erase_vertex_channel -- Channel \"" +
                                  frantic::strings::to_string( channel ) + "\" is not a vertex channel." );
    }

    m_channels.erase( it );
}

void polymesh3::erase_face_channel( const frantic::tstring& channel ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( channel );

    if( it == m_channels.end() ) {
        throw std::runtime_error( "polymesh3::erase_vertex_channel -- No such channel \"" +
                                  frantic::strings::to_string( channel ) + "\" exists." );
    } else if( it->second.is_vertex_channel() ) {
        throw std::runtime_error( "polymesh3::erase_vertex_channel -- Channel \"" +
                                  frantic::strings::to_string( channel ) + "\" is not a face channel." );
    }

    m_channels.erase( it );
}

void polymesh3::copy_channel( const frantic::tstring& destChannel, const frantic::tstring& srcChannel ) {
    std::map<frantic::tstring, polymesh3_channel>::iterator it = m_channels.find( srcChannel );
    if( it == m_channels.end() )
        throw std::runtime_error( "polymesh3::copy_channel() There is no source channel named: \"" +
                                  frantic::strings::to_string( srcChannel ) + "\"" );

    m_channels[destChannel] = it->second;
}

void polymesh3::build_vertex_normals( const frantic::tstring& vertexChannelName, bool normalizeNormals ) {
    const frantic::channels::data_type_t dataType = frantic::channels::data_type_float16;
    const std::size_t dataArity = 3;
    const std::size_t vertexCount = vertex_count();

    // Make sure the channel is there
    if( !has_vertex_channel( vertexChannelName ) )
        add_empty_vertex_channel( vertexChannelName, dataType, dataArity );

    polymesh3_cvt_vertex_accessor<vector3f> normals = get_cvt_vertex_accessor<vector3f>( vertexChannelName );
    polymesh3_vertex_accessor<void> faces = get_vertex_accessor( _T("verts") );

    // Make sure all the normals start as 0
    for( unsigned i = 0; i < vertexCount; ++i )
        normals.set_vertex( i, vector3f( 0.0f ) );

    std::vector<std::size_t> indexVector;
    std::vector<float> weightVector;

    // Add the contributions of all the faces to the normals array
    for( size_t faceIndex = 0; faceIndex < face_count(); ++faceIndex ) {
        const std::size_t cornerCount = faces.get_face_degree( faceIndex );
        frantic::geometry::polymesh3_const_face_range face = faces.get_face( faceIndex );

        if( cornerCount > indexVector.size() ) {
            indexVector.resize( cornerCount );
            weightVector.resize( cornerCount );
        }

        if( cornerCount == 3 ) {
            // Get the triangle
            {
                frantic::geometry::polymesh3_const_face_iterator i = face.first;
                indexVector[0] = *i;
                ++i;
                indexVector[1] = *i;
                ++i;
                indexVector[2] = *i;
                ++i;
            }

            const vector3f& a = get_vertex( indexVector[0] );
            const vector3f& b = get_vertex( indexVector[1] );
            const vector3f& c = get_vertex( indexVector[2] );

            // Get the normal, and use the triangle angles as the weights
            vector3f geoNormal = triangle_normal( a, b, c );
            vector3f weights = get_triangle_angles( a, b, c );
            normals.set_vertex( indexVector[0], normals.get_vertex( indexVector[0] ) + weights.x * geoNormal );
            normals.set_vertex( indexVector[1], normals.get_vertex( indexVector[1] ) + weights.y * geoNormal );
            normals.set_vertex( indexVector[2], normals.get_vertex( indexVector[2] ) + weights.z * geoNormal );
        } else {
            // For now, using Newell's Method for polygons with more than three vertices, in
            // case they have collinear vertices.
            std::size_t index = 0;
            for( frantic::geometry::polymesh3_const_face_iterator i = face.first, ie = face.second; i != ie; ++i ) {
                indexVector[index] = *i;
            }

            frantic::graphics::vector3f first, last;
            first = get_vertex( indexVector.front() );
            last = get_vertex( indexVector.back() );

            vector3f normal;
            vector3f previous;
            vector3f current = last;
            vector3f next = first;
            for( std::size_t cornerIndex = 0, cornerIndexEnd = cornerCount - 1; cornerIndex < cornerIndexEnd;
                 ++cornerIndex ) {
                previous = current;
                current = next;
                next = get_vertex( indexVector[cornerIndex + 1] );

                normal.x += ( current.y - next.y ) * ( current.z + next.z );
                normal.y += ( current.z - next.z ) * ( current.x + next.x );
                normal.z += ( current.x - next.x ) * ( current.y + next.y );
                weightVector[cornerIndex] = get_triangle_angles( previous, current, next ).x;
            }

            previous = current;
            current = next;
            next = first;

            normal.x += ( current.y - next.y ) * ( current.z + next.z );
            normal.y += ( current.z - next.z ) * ( current.x + next.x );
            normal.z += ( current.x - next.x ) * ( current.y + next.y );
            weightVector[cornerCount - 1] = get_triangle_angles( previous, current, next ).x;

            const float length = normal.get_magnitude();
            if( length != 0 ) {
                normal /= length;
            }

            for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
                normals.set_vertex( indexVector[cornerIndex], normals.get_vertex( indexVector[cornerIndex] ) +
                                                                  weightVector[cornerIndex] * normal );
            }
        }
    }

    // Normalize all the resulting vectors
    if( normalizeNormals ) {
        for( unsigned i = 0; i < vertexCount; ++i ) {
            vector3f norm = normals.get_vertex( i );
            norm.normalize();
            normals.set_vertex( i, norm );
        }
    }
}

void polymesh3::build_face_normals( const frantic::tstring& faceChannelName ) {

    using namespace frantic::channels;
    using namespace frantic::graphics;

    if( !has_face_channel( faceChannelName ) ) {
        const data_type_t normalsType = data_type_float16;
        const size_t normalsArity = 3;
        add_empty_face_channel( faceChannelName, normalsType, normalsArity );
    }

    polymesh3_cvt_face_accessor<vector3f> normals = get_cvt_face_accessor<vector3f>( faceChannelName );
    polymesh3_vertex_accessor<void> vertices = get_vertex_accessor( _T("verts") );

    for( size_t i = 0; i < normals.face_count(); ++i ) {

        typedef polymesh3_vertex_accessor<void>::const_face_range face_range;
        face_range faceRange = vertices.get_face( i );

        size_t vertexCount = 0;
        vector3f foundVertices[3];

        // Search for 3 non-co-indicent vertices on the current face
        for( polymesh3_vertex_accessor<void>::const_face_iterator it = faceRange.first;
             it != faceRange.second && vertexCount < 3; ++it ) {

            vector3f vertex = get_vertex( *it );

            bool match = false;
            for( size_t j = 0; j < vertexCount && !match; ++j ) {
                if( ( foundVertices[j] - vertex ).get_magnitude() < 0.0001f ) {
                    match = true;
                    break;
                }
            }

            if( !match ) {
                foundVertices[vertexCount] = vertex;
                ++vertexCount;
            }
        }

        // Apply the face normal, if 3 suitable vertices were found
        if( vertexCount >= 3 ) {
            normals.set_face( i,
                              plane3f::from_triangle( foundVertices[0], foundVertices[1], foundVertices[2] ).normal() );
        } else {
            normals.set_face( i, vector3f( 0.0f ) );
        }
    }
}

namespace {

std::size_t get_data_count( polymesh3_vertex_accessor<void>& acc ) { return acc.vertex_count(); }

char* get_data( polymesh3_vertex_accessor<void>& acc, std::size_t i ) { return acc.get_vertex( i ); }

std::size_t get_data_count( polymesh3_face_accessor<void>& acc ) { return acc.face_count(); }

char* get_data( polymesh3_face_accessor<void>& acc, std::size_t i ) { return acc.get_face( i ); }

template <class AccessorType>
void scale_channel( AccessorType& accessor, double scale ) {
    using namespace frantic::channels;

    const data_type_t dataType = accessor.get_type();
    const std::size_t arity = accessor.get_arity();

    channel_scale_function_t scaleFunction = channel_scale_function( dataType );

    for( std::size_t i = 0, ie = get_data_count( accessor ); i < ie; ++i ) {
        char* data = get_data( accessor, i );
        scaleFunction( scale, data, arity, data );
    }
}

} // anonymous namespace

void scale_channel( polymesh3_ptr mesh, const frantic::tstring& channelName, double scale ) {
    if( !mesh ) {
        throw std::runtime_error( "scale_channel Error: the mesh is NULL" );
    }

    const polymesh3_channel& channel = mesh->get_channel_info( channelName );

    if( channel.is_vertex_channel() ) {
        polymesh3_vertex_accessor<void> acc = mesh->get_vertex_accessor( channelName );
        scale_channel( acc, scale );
    } else {
        polymesh3_face_accessor<void> acc = mesh->get_face_accessor( channelName );
        scale_channel( acc, scale );
    }
}

void transform( polymesh3_ptr mesh, const frantic::graphics::transform4f& xform,
                const frantic::graphics::transform4f& xformTimeDerivative ) {
    if( !mesh ) {
        throw std::runtime_error( "transform() - mesh is NULL" );
    }

    if( !xform.is_identity() ) {
        polymesh3_vertex_accessor<vector3f> geomAcc = mesh->get_vertex_accessor<vector3f>( _T("verts") );

        // Get the Velocity channel
        bool hasVelocityChannel = mesh->has_vertex_channel( _T("Velocity") );
        polymesh3_vertex_accessor<void> velAcc;
        if( hasVelocityChannel ) {
            velAcc = mesh->get_vertex_accessor( _T("Velocity") );
            if( velAcc.has_custom_faces() ) {
                throw std::runtime_error(
                    "transform() - The Velocity channel in the provided mesh has custom faces.  The "
                    "velocities should correspond 1-1 to the vertices, which means they should share the "
                    "primary geometry faces." );
            }
            if( velAcc.vertex_count() != mesh->vertex_count() ) {
                throw std::runtime_error(
                    "transform() - The Velocity channel of the mesh has a different number of entries "
                    "than the vertex count.  The velocities should correspond 1-1 to the vertices." );
            }
            if( velAcc.get_arity() != 3 ) {
                throw std::runtime_error( "transform() - The Velocity channel of the mesh has arity " +
                                          boost::lexical_cast<std::string>( velAcc.get_arity() ) +
                                          ", but it should have arity 3." );
            }
        }

        // Get the Normal channel
        bool hasNormalChannel = mesh->has_vertex_channel( _T("Normal") );
        polymesh3_vertex_accessor<void> normalAcc;
        if( hasNormalChannel ) {
            normalAcc = mesh->get_vertex_accessor( _T("Normal") );

            if( normalAcc.get_arity() != 3 ) {
                throw std::runtime_error( "transform() - The Normal channel of the mesh has arity " +
                                          boost::lexical_cast<std::string>( normalAcc.get_arity() ) +
                                          ", but it should have arity 3." );
            }
        }

        // Apply the transform to the positions and Velocities
        if( hasVelocityChannel ) {
            frantic::channels::channel_type_convertor_function_t convertFromChannel =
                get_channel_type_convertor_function( velAcc.get_type(), frantic::channels::data_type_float32,
                                                     _T("Velocity") );
            frantic::channels::channel_type_convertor_function_t convertToChannel = get_channel_type_convertor_function(
                frantic::channels::data_type_float32, velAcc.get_type(), _T("Velocity") );
            vector3f velocity, transformedVelocity;
            for( std::size_t i = 0, ie = geomAcc.vertex_count(); i != ie; ++i ) {
                convertFromChannel( reinterpret_cast<char*>( &velocity ), velAcc.get_vertex( i ), 3 );
                transformedVelocity =
                    xform.transform_no_translation( velocity ) + xformTimeDerivative * geomAcc.get_vertex( i );
                convertToChannel( velAcc.get_vertex( i ), reinterpret_cast<char*>( &transformedVelocity ), 3 );
                geomAcc.get_vertex( i ) = xform * geomAcc.get_vertex( i );
            }
        } else {
            for( std::size_t i = 0, ie = geomAcc.vertex_count(); i != ie; ++i ) {
                geomAcc.get_vertex( i ) = xform * geomAcc.get_vertex( i );
            }
        }

        // Apply the transform to the Normals
        if( hasNormalChannel ) {
            frantic::channels::channel_type_convertor_function_t convertFromChannel =
                get_channel_type_convertor_function( normalAcc.get_type(), frantic::channels::data_type_float32,
                                                     _T("Normal") );
            frantic::channels::channel_type_convertor_function_t convertToChannel = get_channel_type_convertor_function(
                frantic::channels::data_type_float32, normalAcc.get_type(), _T("Normal") );
            frantic::graphics::transform4f xformInverse = xform.to_inverse();
            vector3f normal, transformedNormal;
            for( std::size_t i = 0, ie = normalAcc.vertex_count(); i != ie; ++i ) {
                convertFromChannel( reinterpret_cast<char*>( &normal ), normalAcc.get_vertex( i ), 3 );
                transformedNormal = xformInverse.transpose_transform_no_translation( normal );
                convertToChannel( normalAcc.get_vertex( i ), reinterpret_cast<char*>( &transformedNormal ), 3 );
            }
        }

        if( xform.is_orientation_inverting() ) {
            reverse_face_winding( mesh );
        }
    } else if( !xformTimeDerivative.is_zero() && mesh->has_vertex_channel( _T("Velocity") ) ) {
        polymesh3_vertex_accessor<vector3f> geomAcc = mesh->get_vertex_accessor<vector3f>( _T("verts") );

        // Get the Velocity channel
        polymesh3_vertex_accessor<void> velAcc = mesh->get_vertex_accessor( _T("Velocity") );
        if( velAcc.has_custom_faces() ) {
            throw std::runtime_error(
                "transform() - The Velocity channel in the provided mesh has custom faces.  The "
                "velocities should correspond 1-1 to the vertices, which means they should share the "
                "primary geometry faces." );
        }
        if( velAcc.vertex_count() != mesh->vertex_count() ) {
            throw std::runtime_error(
                "transform() - The Velocity channel of the mesh has a different number of entries than "
                "the vertex count.  The velocities should correspond 1-1 to the vertices." );
        }
        if( velAcc.get_arity() != 3 ) {
            throw std::runtime_error( "transform() - The Velocity channel of the mesh has arity " +
                                      boost::lexical_cast<std::string>( velAcc.get_arity() ) +
                                      ", but it should have arity 3." );
        }
        frantic::channels::channel_type_convertor_function_t convertFromChannel = get_channel_type_convertor_function(
            velAcc.get_type(), frantic::channels::data_type_float32, _T("Velocity") );
        frantic::channels::channel_type_convertor_function_t convertToChannel = get_channel_type_convertor_function(
            frantic::channels::data_type_float32, velAcc.get_type(), _T("Velocity") );
        vector3f velocity, transformedVelocity;
        for( std::size_t i = 0, ie = geomAcc.vertex_count(); i != ie; ++i ) {
            convertFromChannel( reinterpret_cast<char*>( &velocity ), velAcc.get_vertex( i ), 3 );
            transformedVelocity = velocity + xformTimeDerivative * geomAcc.get_vertex( i );
            convertToChannel( velAcc.get_vertex( i ), reinterpret_cast<char*>( &transformedVelocity ), 3 );
        }
    }
}

void scale( polymesh3_ptr mesh, float scale ) {
    transform( mesh, frantic::graphics::transform4f::from_scale( scale, scale, scale ) );
}

void reverse( const polymesh3_face_range& faceRange ) {
    const std::size_t cornerCount = faceRange.second - faceRange.first;
    const std::size_t halfCornerCount = cornerCount / 2;

    for( std::size_t i = 0; i < halfCornerCount; ++i ) {
        std::swap( *( faceRange.first + i ), *( faceRange.second - 1 - i ) );
    }
}

void reverse_face_winding( polymesh3_ptr mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "reverse_face_winding() - mesh is NULL" );
    }

    frantic::geometry::polymesh3_vertex_accessor<vector3f> geomAcc = mesh->get_vertex_accessor<vector3f>( _T("verts") );

    for( std::size_t faceIndex = 0, faceIndexEnd = geomAcc.face_count(); faceIndex != faceIndexEnd; ++faceIndex ) {
        polymesh3_face_range faceRange = geomAcc.get_face( faceIndex );
        reverse( faceRange );
    }
}

frantic::geometry::polymesh3_ptr linear_interpolate( const frantic::geometry::const_polymesh3_ptr mesh1,
                                                     const frantic::geometry::const_polymesh3_ptr mesh2, float alpha ) {
    if( !mesh1 ) {
        throw std::runtime_error( "linear_interpolate() - The first mesh is NULL" );
    }
    if( !mesh2 ) {
        throw std::runtime_error( "linear_interpolate() - The second mesh is NULL" );
    }
    if( alpha < 0 || alpha > 1 ) {
        throw std::runtime_error(
            "linear_interpolate() - Cannot interpolate.  Alpha value must be in the range [0,1], but instead it is " +
            boost::lexical_cast<std::string>( alpha ) + "." );
    }

    if( mesh1->face_count() != mesh2->face_count() ) {
        throw std::runtime_error(
            "linear_interpolate() - Cannot interpolate, inconsistent face count.  The first mesh has " +
            boost::lexical_cast<std::string>( mesh1->face_count() ) + " faces, while the second mesh has " +
            boost::lexical_cast<std::string>( mesh2->face_count() ) + " faces." );
    }
    if( mesh1->vertex_count() != mesh2->vertex_count() ) {
        throw std::runtime_error(
            "linear_interpolate() - Cannot interpolate, inconsistent vertex count.  The first mesh has " +
            boost::lexical_cast<std::string>( mesh1->vertex_count() ) + " vertices, while the second mesh has " +
            boost::lexical_cast<std::string>( mesh2->vertex_count() ) + " vertices." );
    }

    frantic::geometry::polymesh3_builder builder;

    frantic::geometry::polymesh3_const_vertex_accessor<vector3f> geom1 =
        mesh1->get_const_vertex_accessor<vector3f>( _T("verts") );
    frantic::geometry::polymesh3_const_vertex_accessor<vector3f> geom2 =
        mesh2->get_const_vertex_accessor<vector3f>( _T("verts") );

    std::size_t faceCornerCountSum = 0;

    for( std::size_t i = 0, ie = mesh1->face_count(); i != ie; ++i ) {
        frantic::geometry::polymesh3_const_face_range f1 = geom1.get_face( i );
        frantic::geometry::polymesh3_const_face_range f2 = geom2.get_face( i );

        const std::size_t cornerCount = f1.second - f1.first;
        const std::size_t cornerCount2 = f2.second - f2.first;

        if( cornerCount != cornerCount2 ) {
            throw std::runtime_error(
                "linear_interpolate() - Cannot interpolate, inconsistent number of vertices indices in face " +
                boost::lexical_cast<std::string>( i ) + "." );
        }

        const std::size_t numBytes =
            reinterpret_cast<const char*>( f1.second ) - reinterpret_cast<const char*>( f1.first );
        if( memcmp( f1.first, f2.first, numBytes ) != 0 ) {
            throw std::runtime_error( "linear_interpolate() - Cannot interpolate, inconsistent indices in face " +
                                      boost::lexical_cast<std::string>( i ) + "." );
        }

        faceCornerCountSum += cornerCount;

        builder.add_polygon( f1.first, cornerCount );
    }

    for( std::size_t i = 0, ie = mesh1->vertex_count(); i != ie; ++i ) {
        builder.add_vertex( geom1.get_vertex( i ) + alpha * ( geom2.get_vertex( i ) - geom1.get_vertex( i ) ) );
    }

    polymesh3_ptr result = builder.finalize();

    for( polymesh3::iterator i = mesh1->begin(), ie = mesh1->end(); i != ie; ++i ) {
        const frantic::tstring& channelName = i->first;
        if( i->first != _T("verts") && mesh2->has_channel( i->first ) ) {
            const polymesh3_channel& ch1 = i->second;
            const polymesh3_channel& ch2 = mesh2->get_channel_info( i->first );

            if( ch1.is_vertex_channel() == ch2.is_vertex_channel() ) {
                if( ch1.arity() != ch2.arity() ) {
                    throw std::runtime_error(
                        std::string( "linear_interpolate() - Cannot interpolate, inconsistent arity in " ) +
                        ( ch1.is_vertex_channel() ? "vertex" : "face" ) + " channel " +
                        frantic::strings::to_string( channelName ) + ".  First mesh has arity " +
                        boost::lexical_cast<std::string>( ch1.arity() ) + ", while second mesh has arity " +
                        boost::lexical_cast<std::string>( ch2.arity() ) + "." );
                }
                if( ch1.type() != ch2.type() ) {
                    throw std::runtime_error(
                        std::string( "linear_interpolate() - Cannot interpolate, inconsistent data type in " ) +
                        ( ch1.is_vertex_channel() ? "vertex" : "face" ) + " channel " +
                        frantic::strings::to_string( channelName ) + ".  First mesh has data type " +
                        frantic::strings::to_string( frantic::channels::channel_data_type_str( ch1.type() ) ) +
                        ", while second mesh has data type " +
                        frantic::strings::to_string( frantic::channels::channel_data_type_str( ch2.type() ) ) + "." );
                }
                const std::size_t arity = ch1.arity();
                const frantic::channels::data_type_t dataType = ch1.type();

                if( ch1.is_vertex_channel() ) {
                    polymesh3_const_vertex_accessor<void> acc1 = mesh1->get_const_vertex_accessor( i->first );
                    polymesh3_const_vertex_accessor<void> acc2 = mesh2->get_const_vertex_accessor( i->first );

                    // TODO: I'm not sure if this should be an error
                    if( acc1.has_custom_faces() != acc2.has_custom_faces() ) {
                        throw std::runtime_error( "linear_interpolate() - Cannot interpolate vertex channel " +
                                                  frantic::strings::to_string( channelName ) + ".  The first mesh " +
                                                  ( acc1.has_custom_faces() ? "has" : "does not have" ) +
                                                  " custom faces, while the second mesh " +
                                                  ( acc1.has_custom_faces() ? "has" : "does not have" ) +
                                                  " custom faces." );
                    }

                    // TODO: I don't think this should be an error.
                    if( acc1.vertex_count() != acc2.vertex_count() ) {
                        throw std::runtime_error(
                            "linear_interpolate() - Mismatch in vertex count for vertex channel " +
                            frantic::strings::to_string( channelName ) + ".  First mesh has " +
                            boost::lexical_cast<std::string>( acc1.vertex_count() ) +
                            " vertices while the second mesh has " +
                            boost::lexical_cast<std::string>( acc2.vertex_count() ) + " vertices." );
                    }

                    const std::size_t vertexCount = acc1.vertex_count();

                    std::vector<int> inFaceBuffer;
                    std::vector<int>* pInFaceBuffer = 0;

                    if( acc1.has_custom_faces() ) {
                        inFaceBuffer.reserve( faceCornerCountSum );

                        if( acc1.face_count() != acc2.face_count() ) {
                            throw std::runtime_error(
                                "linear_interpolate() - Mismatch in face count for vertex channel " +
                                frantic::strings::to_string( channelName ) + ".  First mesh has " +
                                boost::lexical_cast<std::string>( acc1.face_count() ) +
                                " faces while the second mesh has " +
                                boost::lexical_cast<std::string>( acc2.face_count() ) + " faces." );
                        }

                        for( std::size_t faceIndex = 0, faceIndexEnd = acc1.face_count(); faceIndex != faceIndexEnd;
                             ++faceIndex ) {
                            // copy custom faces
                            frantic::geometry::polymesh3_const_face_range f1 = acc1.get_face( faceIndex );
                            frantic::geometry::polymesh3_const_face_range f2 = acc2.get_face( faceIndex );

                            const std::size_t cornerCount = f1.second - f1.first;
                            const std::size_t cornerCount2 = f2.second - f2.first;

                            if( cornerCount != cornerCount2 ) {
                                throw std::runtime_error( "linear_interpolate() - Cannot interpolate, inconsistent "
                                                          "number of vertices in face " +
                                                          boost::lexical_cast<std::string>( faceIndex ) + "." );
                            }

                            const std::size_t numBytes =
                                reinterpret_cast<const char*>( f1.second ) - reinterpret_cast<const char*>( f1.first );
                            if( memcmp( f1.first, f2.first, numBytes ) != 0 ) {
                                throw std::runtime_error(
                                    "linear_interpolate() - Cannot interpolate, inconsistent indices in face " +
                                    boost::lexical_cast<std::string>( faceIndex ) + "." );
                            }

                            for( polymesh3_const_face_iterator i = f1.first, ie = f1.second; i != ie; ++i ) {
                                inFaceBuffer.push_back( *i );
                            }
                        }

                        pInFaceBuffer = &inFaceBuffer;
                    }

                    // copy verts
                    frantic::graphics::raw_byte_buffer dataBuffer;
                    dataBuffer.resize( vertexCount * arity * frantic::channels::sizeof_channel_data_type( dataType ) );

                    if( vertexCount > 0 ) {
                        frantic::channels::channel_weighted_sum_combine_function_t ws =
                            channel_weighted_sum_combine_function( dataType );
                        const char* data[2] = { acc1.get_vertex( 0 ), acc2.get_vertex( 0 ) };
                        float weights[2] = { 1.f - alpha, alpha };
                        ws( weights, data, 2, vertexCount * arity, dataBuffer.ptr_at( 0 ) );
                    }

                    result->add_vertex_channel( channelName, dataType, arity, dataBuffer, pInFaceBuffer );
                } else {
                    polymesh3_const_face_accessor<void> acc1 = mesh1->get_const_face_accessor( channelName );
                    polymesh3_const_face_accessor<void> acc2 = mesh2->get_const_face_accessor( channelName );

                    if( acc1.face_count() != acc2.face_count() ) {
                        throw std::runtime_error( "linear_interpolate() - Mismatch in face count for face channel " +
                                                  frantic::strings::to_string( channelName ) + ".  First mesh has " +
                                                  boost::lexical_cast<std::string>( acc1.face_count() ) +
                                                  " faces while the second mesh has " +
                                                  boost::lexical_cast<std::string>( acc2.face_count() ) + " faces." );
                    }

                    const std::size_t faceCount = acc1.face_count();

                    frantic::graphics::raw_byte_buffer dataBuffer;
                    dataBuffer.resize( faceCount * arity * frantic::channels::sizeof_channel_data_type( dataType ) );

                    if( faceCount > 0 ) {
                        frantic::channels::channel_weighted_sum_combine_function_t ws =
                            channel_weighted_sum_combine_function( dataType );
                        const char* data[2] = { acc1.get_face( 0 ), acc2.get_face( 0 ) };
                        float weights[2] = { 1.f - alpha, alpha };
                        ws( weights, data, 2, faceCount * arity, dataBuffer.ptr_at( 0 ) );
                    }

                    result->add_face_channel( channelName, dataType, arity, dataBuffer );
                }
            }
        }
    }

    return result;
}

// Culling

void copy_face_channel_culled( polymesh3_ptr& outMesh, const polymesh3_ptr& inMesh, const frantic::tstring& channelName,
                               const std::vector<bool>& usedFace, frantic::logging::progress_logger& logger ) {
    if( !outMesh ) {
        throw std::runtime_error( "copy_face_channel_culled Error: output mesh is NULL" );
    }
    if( !inMesh ) {
        throw std::runtime_error( "copy_face_channel_culled Error: input mesh is NULL" );
    }

    if( inMesh->face_count() != usedFace.size() ) {
        throw std::runtime_error( "copy_face_channel_culled Error: mismatch between number of faces in input mesh and "
                                  "number of entries in face tag array" );
    }

    const size_t updateInterval = 1024;

    const std::size_t outFaceCount = outMesh->face_count();
    const std::size_t inFaceCount = inMesh->face_count();

    const frantic::geometry::polymesh3_channel& inChannelInfo = inMesh->get_channel_info( channelName );
    const std::size_t elementSize = inChannelInfo.element_size();

    frantic::graphics::raw_byte_buffer data;
    data.resize( outFaceCount * elementSize );

    polymesh3_const_face_accessor<void> acc = inMesh->get_const_face_accessor( channelName );
    if( acc.face_count() != inFaceCount ) {
        throw std::runtime_error(
            "combine_face_channel_culled Error: mismatch between geometry and channel face count" );
    }

    std::size_t outIndex = 0;
    for( std::size_t inIndex = 0; inIndex < inFaceCount; ++inIndex ) {
        if( usedFace[inIndex] ) {
            if( outIndex >= outFaceCount ) {
                throw std::runtime_error(
                    "copy_face_channel_culled Error: more tagged faces than faces in output mesh" );
            }
            memcpy( data.ptr_at( outIndex * elementSize ), acc.get_face( inIndex ), elementSize );
            ++outIndex;
        }

        if( inIndex % updateInterval == 0 ) {
            logger.update_progress( inIndex, inFaceCount );
        }
    }

    if( outIndex != outFaceCount ) {
        throw std::runtime_error( "copy_face_channel_culled Error: mismatch between number of faces in output mesh and "
                                  "number of tagged faces" );
    }

    outMesh->add_face_channel( channelName, inChannelInfo.type(), inChannelInfo.arity(), data );

    logger.update_progress( 100.0f );
}

void copy_face_channel_culled( polymesh3_ptr& outMesh, const polymesh3_ptr& inMesh, const frantic::tstring& channelName,
                               const std::vector<bool>& usedFace ) {
    frantic::logging::null_progress_logger nullLogger;
    copy_face_channel_culled( outMesh, inMesh, channelName, usedFace, nullLogger );
}

void copy_custom_vertex_channel_culled( polymesh3_ptr& outMesh, const polymesh3_ptr& inMesh,
                                        const frantic::tstring& channelName, const std::vector<bool>& keepFaces,
                                        ::frantic::logging::progress_logger& logger ) {
    if( !outMesh ) {
        throw( "copy_custom_vertex_channel_culled Error: output mesh is NULL" );
    }

    if( !inMesh ) {
        throw( "copy_custom_vertex_channel_culled Error: input mesh is NULL" );
    }

    if( inMesh->face_count() != keepFaces.size() ) {
        throw( "copy_custom_vertex_channel_culled Error: mismatch between number of faces in input mesh, "
               "and size of face tag array" );
    }

    frantic::geometry::polymesh3_const_vertex_accessor<void> acc = inMesh->get_const_vertex_accessor( channelName );

    if( acc.face_count() != keepFaces.size() ) {
        throw( "copy_custom_vertex_channel_culled Error: mismatch between number of faces in channel "
               "\"" +
               frantic::strings::to_string( channelName ) + "\", and size of face tag array" );
    }

    const std::size_t elementSize = acc.get_arity() * frantic::channels::sizeof_channel_data_type( acc.get_type() );

    const size_t updateInterval = 1024;
    frantic::logging::progress_logger_subinterval_tracker subinterval( logger, 0.0f, 20.0f );

    std::vector<bool> keepVertices( acc.vertex_count(), false );
    std::size_t faceCornerCount = 0;
    for( std::size_t faceIndex = 0, faceCount = acc.face_count(); faceIndex < faceCount; ++faceIndex ) {
        if( keepFaces[faceIndex] ) {
            frantic::geometry::polymesh3_const_face_range face = acc.get_face( faceIndex );
            for( const int* corner = face.first; corner != face.second; ++corner ) {
                keepVertices[*corner] = true;
                ++faceCornerCount;
            }
        }

        if( faceIndex % updateInterval == 0 ) {
            logger.update_progress( faceIndex, faceCount );
        }
    }

    subinterval.reset( 20.0f, 40.0f );

    std::vector<int> vertexRemap( acc.vertex_count() );
    int usedVertexCount = 0;
    for( std::size_t vertexIndex = 0, vertexIndexEnd = keepVertices.size(); vertexIndex < vertexIndexEnd;
         ++vertexIndex ) {
        if( keepVertices[vertexIndex] ) {
            vertexRemap[vertexIndex] = usedVertexCount;
            ++usedVertexCount;
        }

        if( vertexIndex % updateInterval == 0 ) {
            logger.update_progress( vertexIndex, vertexIndexEnd );
        }
    }

    frantic::graphics::raw_byte_buffer data;
    data.resize( usedVertexCount * elementSize );

    subinterval.reset( 40.0f, 70.0f );

    for( std::size_t vertexIndex = 0, vertexIndexEnd = acc.vertex_count(); vertexIndex < vertexIndexEnd;
         ++vertexIndex ) {
        if( keepVertices[vertexIndex] ) {
            memcpy( data.ptr_at( vertexRemap[vertexIndex] * elementSize ), acc.get_vertex( vertexIndex ), elementSize );
        }

        if( vertexIndex % updateInterval == 0 ) {
            logger.update_progress( vertexIndex, vertexIndexEnd );
        }
    }

    subinterval.reset( 70.0f, 100.0f );

    std::vector<int> customFaces;
    customFaces.reserve( faceCornerCount );
    for( std::size_t faceIndex = 0, faceCount = acc.face_count(); faceIndex < faceCount; ++faceIndex ) {
        if( keepFaces[faceIndex] ) {
            frantic::geometry::polymesh3_const_face_range face = acc.get_face( faceIndex );
            for( const int* corner = face.first; corner != face.second; ++corner ) {
                customFaces.push_back( vertexRemap[*corner] );
            }
        }

        if( faceIndex % updateInterval == 0 ) {
            logger.update_progress( faceIndex, faceCount );
        }
    }

    outMesh->add_vertex_channel( channelName, acc.get_type(), acc.get_arity(), data, &customFaces );

    logger.update_progress( 100.0f );
}

void copy_custom_vertex_channel_culled( polymesh3_ptr& outMesh, const polymesh3_ptr& inMesh,
                                        const frantic::tstring& channelName, const std::vector<bool>& keepFaces ) {
    frantic::logging::null_progress_logger nullLogger;
    copy_custom_vertex_channel_culled( outMesh, inMesh, channelName, keepFaces, nullLogger );
}

void copy_simple_vertex_channel_culled( polymesh3_ptr& outMesh, const polymesh3_ptr& inMesh,
                                        const frantic::tstring& channelName, const std::vector<bool>& keepVertices,
                                        frantic::logging::progress_logger& logger ) {
    if( !outMesh ) {
        throw std::runtime_error( "copy_simple_vertex_channel_culled Error: output mesh is NULL" );
    }

    if( !inMesh ) {
        throw std::runtime_error( "copy_simple_vertex_channel_culled Error: input mesh is NULL" );
    }

    const size_t updateInterval = 1024;

    polymesh3_const_vertex_accessor<void> acc = inMesh->get_const_vertex_accessor( channelName );

    const std::size_t elementSize = acc.get_arity() * frantic::channels::sizeof_channel_data_type( acc.get_type() );

    if( acc.vertex_count() != inMesh->vertex_count() ) {
        throw std::runtime_error( "copy_simple_vertex_channel_culled Error: mismatch between "
                                  "number of channel data entries and number of mesh vertices" );
    }
    if( acc.vertex_count() != keepVertices.size() ) {
        throw std::runtime_error( "copy_simple_vertex_channel_culled Error: mismatch between number of "
                                  "channel data entries and number of entries in keep vertices array" );
    }

    frantic::graphics::raw_byte_buffer data;
    data.resize( elementSize * outMesh->vertex_count() );

    std::size_t outIndex = 0;
    for( std::size_t i = 0, ie = acc.vertex_count(); i < ie; ++i ) {
        if( keepVertices[i] ) {
            memcpy( data.ptr_at( outIndex * elementSize ), acc.get_vertex( i ), elementSize );
            ++outIndex;
        }

        if( i % updateInterval == 0 ) {
            logger.update_progress( i, acc.vertex_count() );
        }
    }

    if( outIndex != outMesh->vertex_count() ) {
        throw std::runtime_error( "copy_simple_vertex_channel_culled Error: mismatch between number of "
                                  "tagged vertices and the number of vertices in output mesh" );
    }

    outMesh->add_vertex_channel( channelName, acc.get_type(), acc.get_arity(), data );

    logger.update_progress( 100.0f );
}

void copy_simple_vertex_channel_culled( polymesh3_ptr& outMesh, const polymesh3_ptr& inMesh,
                                        const frantic::tstring& channelName, const std::vector<bool>& keepVertices ) {
    frantic::logging::null_progress_logger nullLogger;
    copy_simple_vertex_channel_culled( outMesh, inMesh, channelName, keepVertices, nullLogger );
}

// Alternative bounding box contains checking to be used with Sequoia's hacksaw functionality
// Required as a "tie breaker" for what happens when a point lies at the edge (hacksaw boundary).
//   A naive check would either cause the point to fall in both hacksaw cells or neither
inline bool bounding_box_contains_impl( const frantic::graphics::boundbox3f& bbox,
                                        const frantic::graphics::transform4f& bboxInverseTransform,
                                        const frantic::graphics::vector3f& testPosition ) {
    if( bbox.is_empty() )
        return true;

    const frantic::graphics::vector3f& minPosition = bbox.minimum();
    const frantic::graphics::vector3f& maxPosition = bbox.maximum();
    const frantic::graphics::vector3f checkPosition( bboxInverseTransform * testPosition );

    return minPosition.x <= checkPosition.x && checkPosition.x < maxPosition.x && minPosition.y <= checkPosition.y &&
           checkPosition.y < maxPosition.y && minPosition.z <= checkPosition.z && checkPosition.z < maxPosition.z;
}

bool check_face_region_impl( const frantic::graphics::boundbox3f& triangleBBox,
                             const frantic::graphics::boundbox3f& bbox,
                             const frantic::graphics::transform4f& bboxInverseTransform ) {
    return !bounding_box_contains_impl( bbox, bboxInverseTransform, triangleBBox.center() );
}

bool check_face_region_impl( const polymesh3_const_face_range& inFace,
                             const frantic::geometry::polymesh3_const_vertex_accessor<vector3f>& geomAcc,
                             const frantic::graphics::boundbox3f& bbox,
                             const frantic::graphics::transform4f& bboxInverseTransform, bool& outIsOnBoundary ) {

    bool containsInside = false;
    bool containsOutside = false;

    frantic::graphics::boundbox3f triangleBBox;
    for( polymesh3_const_face_iterator iter = inFace.first; iter < inFace.second; ++iter ) {
        const frantic::graphics::vector3f& current = geomAcc.get_vertex( *iter );
        triangleBBox += current;

        const bool insideBounds = bounding_box_contains_impl( bbox, bboxInverseTransform, current );
        if( insideBounds )
            containsInside = true;
        else
            containsOutside = true;
    }

    outIsOnBoundary = containsInside && containsOutside;
    return check_face_region_impl( triangleBBox, bbox, bboxInverseTransform );
}

bool check_face_region_impl( const frantic::graphics::vector3f& v1, const frantic::graphics::vector3f& v2,
                             const frantic::graphics::vector3f& v3, const frantic::graphics::boundbox3f& bbox,
                             const frantic::graphics::transform4f& bboxInverseTransform, bool& outIsOnBoundary ) {

    bool containsInside = false;
    bool containsOutside = false;
    bool insideBounds;

    frantic::graphics::boundbox3f triangleBBox;
    triangleBBox += v1;
    triangleBBox += v2;
    triangleBBox += v3;

    insideBounds = bounding_box_contains_impl( bbox, bboxInverseTransform, v1 );
    if( insideBounds )
        containsInside = true;
    else
        containsOutside = true;
    insideBounds = bounding_box_contains_impl( bbox, bboxInverseTransform, v2 );
    if( insideBounds )
        containsInside = true;
    else
        containsOutside = true;
    insideBounds = bounding_box_contains_impl( bbox, bboxInverseTransform, v3 );
    if( insideBounds )
        containsInside = true;
    else
        containsOutside = true;

    outIsOnBoundary = containsInside && containsOutside;
    return check_face_region_impl( triangleBBox, bbox, bboxInverseTransform );
}

polymesh3_ptr cull_geometry_impl( polymesh3_ptr mesh, const frantic::graphics::boundbox3f& bbox,
                                  const frantic::graphics::transform4f& bboxInverseTransform ) {
    frantic::geometry::polymesh3_builder builder;

    frantic::geometry::polymesh3_const_vertex_accessor<vector3f> geomAcc =
        mesh->get_const_vertex_accessor<vector3f>( _T("verts") );
    const std::size_t faceCount = geomAcc.face_count();

    // Determine which faces are used
    std::vector<bool> usedFaces( faceCount, false );
    for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        polymesh3_const_face_range inFace = geomAcc.get_face( faceIndex );

        bool onBounds;
        const bool isUsed = !check_face_region_impl( inFace, geomAcc, bbox, bboxInverseTransform, onBounds );
        if( isUsed ) {
            usedFaces[faceIndex] = true;
        }
    }

    return cull_faces( mesh, usedFaces );
}

polymesh3_ptr combine( const std::vector<polymesh3_ptr>& meshes ) {
    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces;
    meshInterfaces.reserve( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        if( !meshes[meshIndex] ) {
            continue;
        }
        meshInterfaces.push_back(
            mesh_interface_ptr( polymesh3_interface::create_instance( meshes[meshIndex] ).release() ) );
    }

    return combine( meshInterfaces );
}

frantic::geometry::polymesh3_ptr weld( const std::vector<polymesh3_ptr>& meshes, float errorTolerance ) {
    frantic::logging::null_progress_logger progress;
    return weld( meshes, progress, errorTolerance );
}

frantic::geometry::polymesh3_ptr weld( const std::vector<polymesh3_ptr>& meshes,
                                       ::frantic::logging::progress_logger& logger, float errorTolerance ) {
    const std::size_t meshCount = meshes.size();
    std::vector<mesh_interface_ptr> meshInterfaces( meshCount );
    for( std::size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        if( !meshes[meshIndex] ) {
            continue;
        }
        meshInterfaces.push_back(
            mesh_interface_ptr( polymesh3_interface::create_instance( meshes[meshIndex] ).release() ) );
    }

    return weld_vertices( meshInterfaces, errorTolerance, logger );
}

polymesh3_ptr cull_faces( polymesh3_ptr mesh, const std::vector<bool>& keepFaces ) {
    frantic::logging::null_progress_logger nullLogger;
    return cull_faces( mesh, keepFaces, nullLogger );
}

polymesh3_ptr cull_faces( polymesh3_ptr mesh, const std::vector<bool>& keepFaces,
                          ::frantic::logging::progress_logger& logger ) {
    if( !mesh ) {
        throw std::runtime_error( "cull_faces Error: mesh is NULL" );
    }
    if( mesh->face_count() != keepFaces.size() ) {
        throw std::runtime_error(
            "cull_faces Error: mismatch between number of faces in mesh, and size of face tag array" );
    }

    frantic::geometry::polymesh3_builder builder;

    frantic::geometry::polymesh3_const_vertex_accessor<vector3f> geomAcc =
        mesh->get_const_vertex_accessor<vector3f>( _T("verts") );
    const std::size_t vertexCount = geomAcc.vertex_count();
    const std::size_t faceCount = geomAcc.face_count();

    const size_t updateInterval = 1024;
    const size_t iterationIntervals = mesh->channel_count() + 3;
    const float iteratorIntervalSize = 100.0f / float( iterationIntervals );
    frantic::logging::progress_logger_subinterval_tracker subinterval( logger, 0.0f, iteratorIntervalSize );

    // Determine which vertex/face is used
    std::vector<bool> usedVertices( vertexCount, false );
    for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        polymesh3_const_face_range inFace = geomAcc.get_face( faceIndex );

        const bool isUsed = keepFaces[faceIndex];
        if( isUsed ) {
            for( polymesh3_const_face_iterator iter = inFace.first; iter < inFace.second; ++iter ) {
                usedVertices[*iter] = true;
            }
        }

        if( faceIndex % updateInterval == 0 ) {
            logger.update_progress( faceIndex, faceCount );
        }
    }

    subinterval.reset( iteratorIntervalSize, iteratorIntervalSize * 2 );

    // Add the vertices
    std::vector<std::size_t> remapVertex( vertexCount );
    std::size_t usedVerticesCount = 0;
    for( std::size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
        if( usedVertices[vertexIndex] ) {
            builder.add_vertex( geomAcc.get_vertex( vertexIndex ) );
            remapVertex[vertexIndex] = usedVerticesCount;
            ++usedVerticesCount;
        }

        if( vertexIndex % updateInterval == 0 ) {
            logger.update_progress( vertexIndex, vertexCount );
        }
    }

    subinterval.reset( iteratorIntervalSize * 2, iteratorIntervalSize * 3 );

    // Add the faces
    {
        std::vector<int> outFace;
        for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
            if( keepFaces[faceIndex] ) {
                polymesh3_const_face_range inFace = geomAcc.get_face( faceIndex );
                outFace.clear();
                for( polymesh3_const_face_iterator iter = inFace.first; iter < inFace.second; ++iter ) {
                    outFace.push_back( int( remapVertex[*iter] ) );
                }
                builder.add_polygon( outFace );
            }

            if( faceIndex % updateInterval == 0 ) {
                logger.update_progress( faceIndex, faceCount );
            }
        }
    }

    // Create the base geometry
    frantic::geometry::polymesh3_ptr result = builder.finalize();

    logger.update_progress( iteratorIntervalSize * 4 );

    size_t currentInterval = 4;

    // Merge channels
    for( polymesh3::iterator i = mesh->begin(), ie = mesh->end(); i != ie; ++i ) {
        if( i->first == _T("verts") ) {
            continue;
        }

        subinterval.reset( iteratorIntervalSize * currentInterval, iteratorIntervalSize * ( currentInterval + 1 ) );

        if( i->second.is_vertex_channel() ) {
            if( i->second.get_faces().empty() ) {
                copy_simple_vertex_channel_culled( result, mesh, i->first, usedVertices, logger );
            } else {
                copy_custom_vertex_channel_culled( result, mesh, i->first, keepFaces, logger );
            }
        } else {
            copy_face_channel_culled( result, mesh, i->first, keepFaces, logger );
        }

        ++currentInterval;
    }

    return result;
}

polymesh3_ptr cull_geometry( polymesh3_ptr mesh, const frantic::graphics::boundbox3f& bbox,
                             const frantic::graphics::transform4f& bboxTransform ) {
    return cull_geometry_impl( mesh, bbox, bboxTransform.to_inverse() );
}

bool check_face_region( const frantic::graphics::vector3f& v1, const frantic::graphics::vector3f& v2,
                        const frantic::graphics::vector3f& v3, const frantic::graphics::boundbox3f& bbox,
                        const frantic::graphics::transform4f& bboxInverseTransform, bool& outIsOnBoundary ) {
    return check_face_region_impl( v1, v2, v3, bbox, bboxInverseTransform, outIsOnBoundary );
}

namespace {

void copy_vertices( frantic::geometry::polymesh3_builder& builder, const frantic::geometry::mesh_interface* mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "copy_vertices Error: mesh is NULL" );
    }

    const std::size_t vertexCount = mesh->get_num_verts();

    for( std::size_t i = 0; i < vertexCount; ++i ) {
        float buffer[3];
        mesh->get_vert( i, buffer );
        builder.add_vertex( buffer[0], buffer[1], buffer[2] );
    }
}

void copy_faces( frantic::geometry::polymesh3_builder& builder, const frantic::geometry::mesh_interface* mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "copy_faces Error: mesh is NULL" );
    }

    const std::size_t faceCount = mesh->get_num_faces();

    std::vector<std::size_t> indexBuffer;
    std::vector<int> intIndexBuffer;
    for( std::size_t i = 0; i < faceCount; ++i ) {
        const std::size_t cornerCount = mesh->get_num_face_verts( i );
        if( cornerCount == 0 ) {
            throw std::runtime_error( "copy_faces Error: cannot create face with zero corners" );
        }
        if( indexBuffer.size() < cornerCount ) {
            indexBuffer.resize( cornerCount );
            intIndexBuffer.resize( cornerCount );
        }

        mesh->get_face_vert_indices( i, &indexBuffer[0] );
        for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
            intIndexBuffer[cornerIndex] = static_cast<int>( indexBuffer[cornerIndex] );
        }

        builder.add_polygon( &intIndexBuffer[0], cornerCount );
    }
}

void copy_channel_data( frantic::graphics::raw_byte_buffer& outBuffer,
                        const frantic::geometry::mesh_channel* channel ) {
    if( !channel ) {
        throw std::runtime_error( "copy_channel_data Error: channel is NULL" );
    }

    const std::size_t elementSize = channel->get_element_size();
    const std::size_t elementCount = channel->get_num_elements();

    if( elementSize == 0 ) {
        throw std::runtime_error( "copy_channel_data Error: element size is 0" );
    }

    outBuffer.resize( elementSize * elementCount );

    boost::scoped_array<char> elementBuffer( new char[elementSize] );
    for( std::size_t i = 0; i < elementCount; ++i ) {
        channel->get_value( i, elementBuffer.get() );
        memcpy( outBuffer.ptr_at( elementSize * i ), elementBuffer.get(), elementSize );
    }
}

void copy_channel_indices( std::vector<int>& outBuffer, const frantic::geometry::mesh_channel* channel ) {
    if( !channel ) {
        throw std::runtime_error( "copy_channel_indices Error: channel is NULL" );
    }

    const std::size_t faceCount = channel->get_num_faces();

    outBuffer.resize( 0 );
    outBuffer.reserve( 3 * faceCount );

    for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        const std::size_t cornerCount = channel->get_num_face_verts( faceIndex );

        for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
            const std::size_t vertexIndex = channel->get_fv_index( faceIndex, cornerIndex );
            outBuffer.push_back( static_cast<int>( vertexIndex ) );
        }
    }
}

bool has_custom_faces( frantic::geometry::polymesh3_ptr mesh, const frantic::geometry::mesh_channel* channel ) {
    if( !mesh ) {
        throw std::runtime_error( "has_custom_faces Error: mesh is NULL" );
    }
    if( !channel ) {
        throw std::runtime_error( "has_custom_faces Error: channel is NULL" );
    }

    bool hasCustomFaces = false;

    const std::size_t faceCount = mesh->face_count();

    if( faceCount != channel->get_num_faces() ) {
        throw std::runtime_error( "has_custom_faces Error: mismatch between number of mesh faces (" +
                                  boost::lexical_cast<std::string>( mesh->face_count() ) +
                                  ") and number of channel faces (" +
                                  boost::lexical_cast<std::string>( channel->get_num_faces() ) + ")" );
    }

    if( channel->get_channel_type() == mesh_channel::vertex ) {
        if( mesh->vertex_count() != channel->get_num_elements() ) {
            throw frantic::exception_stream()
                << "has_custom_faces Error: channel is a vertex channel, but there is a mismatch "
                << "between the number of mesh vertices (" << mesh->vertex_count() << ") and the "
                << "number of data elements (" << channel->get_num_elements() << ")";
        }

        return false;
    }

    frantic::geometry::polymesh3_const_vertex_accessor<vector3f> geomAcc =
        mesh->get_const_vertex_accessor<vector3f>( _T("verts") );

    for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        frantic::geometry::polymesh3_const_vertex_accessor<vector3f>::const_face_range face =
            geomAcc.get_face( faceIndex );

        const std::size_t cornerCount = face.second - face.first;
        if( cornerCount != channel->get_num_face_verts( faceIndex ) ) {
            throw std::runtime_error( "has_custom_faces Error: mismatch in degree of mesh face (" +
                                      boost::lexical_cast<std::string>( cornerCount ) + ") and channel face (" +
                                      boost::lexical_cast<std::string>( channel->get_num_face_verts( faceIndex ) ) +
                                      ")" );
        }

        for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
            if( static_cast<size_t>( face.first[cornerIndex] ) != channel->get_fv_index( faceIndex, cornerIndex ) ) {
                hasCustomFaces = true;
            }
        }
    }

    return hasCustomFaces;
}

void copy_vertex_channels( frantic::geometry::polymesh3_ptr outMesh, const frantic::geometry::mesh_interface* inMesh ) {
    if( !outMesh ) {
        throw std::runtime_error( "copy_vertex_channels Error: outMesh is NULL" );
    }
    if( !inMesh ) {
        throw std::runtime_error( "copy_vertex_channels Error: inMesh is NULL" );
    }

    using frantic::geometry::mesh_interface;

    const mesh_interface::mesh_channel_map& vertexChannels = inMesh->get_vertex_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = vertexChannels.begin(), ie = vertexChannels.end();
         i != ie; ++i ) {
        const frantic::geometry::mesh_channel* channel = i->second;

        frantic::graphics::raw_byte_buffer dataBuffer;
        copy_channel_data( dataBuffer, channel );

        std::vector<int>* pIndexBuffer = 0; // we change this to point to indexBuffer if the channel has custom faces
        std::vector<int> indexBuffer;
        if( has_custom_faces( outMesh, channel ) ) {
            copy_channel_indices( indexBuffer, channel );
            pIndexBuffer = &indexBuffer;
        }

        outMesh->add_vertex_channel( channel->get_name(), channel->get_data_type(), channel->get_data_arity(),
                                     dataBuffer, pIndexBuffer );
    }
}

void copy_face_channels(
    frantic::geometry::polymesh3_ptr outMesh, const frantic::geometry::mesh_interface* inMesh,
    const frantic::channels::channel_propagation_policy& cpp = frantic::channels::channel_propagation_policy() ) {
    if( !outMesh ) {
        throw std::runtime_error( "copy_face_channels Error: outMesh is NULL" );
    }
    if( !inMesh ) {
        throw std::runtime_error( "copy_face_channels Error: inMesh is NULL" );
    }

    using frantic::geometry::mesh_interface;

    const mesh_interface::mesh_channel_map& faceChannels = inMesh->get_face_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = faceChannels.begin(), ie = faceChannels.end(); i != ie;
         ++i ) {
        const frantic::geometry::mesh_channel* channel = i->second;

        if( !cpp.is_channel_included( i->first ) ) {
            continue;
        }

        frantic::graphics::raw_byte_buffer buffer;
        if( channel->get_num_elements() != outMesh->face_count() ) {
            throw std::runtime_error( "copy_face_channels Error: mismatch between number of faces in geometry (" +
                                      boost::lexical_cast<std::string>( outMesh->face_count() ) +
                                      ") and elements in face channel \"" +
                                      frantic::strings::to_string( channel->get_name() ) + "\" (" +
                                      boost::lexical_cast<std::string>( channel->get_num_elements() ) + ")" );
        }
        copy_channel_data( buffer, channel );
        outMesh->add_face_channel( channel->get_name(), channel->get_data_type(), channel->get_data_arity(), buffer );
    }
}

} // anonymous namespace

frantic::geometry::polymesh3_ptr create_polymesh3( const frantic::geometry::mesh_interface* mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "create_polymesh3 Error: mesh is NULL" );
    }

    frantic::geometry::polymesh3_builder builder;

    copy_vertices( builder, mesh );
    copy_faces( builder, mesh );

    frantic::geometry::polymesh3_ptr result = builder.finalize();

    copy_vertex_channels( result, mesh );
    copy_face_channels( result, mesh );

    return result;
}

namespace {

/**
 *  Map arrays of channel indices to a unique id.  This is intended for finding
 * all of the combinations of custom face indices in a mesh.
 */
class index_array_to_id {
  public:
    typedef int index_t;

    /**
     * @param vertexCount one past the greatest vertex index that will be
     *    inserted into this object.
     * @param channelCount the number of indices in each data entry, including
     *    one for the vertex index.
     */
    index_array_to_id( std::size_t vertexCount, std::size_t channelCount )
        : m_entrySize( channelCount + 1 )
        , m_dataSizeInBytes( channelCount * sizeof( index_t ) )
        , m_hashTableEntryCount( 2 )
        , m_hashTableOccupancy( 0 )
        , m_nextId( 0 )
        , m_channelCount( channelCount )
        , m_vertexCount( vertexCount ) {
        if( channelCount < 1 ) {
            throw std::runtime_error( "index_array_to_id Error: must have one or more channels" );
        }

        assert( is_power_of_two( m_hashTableEntryCount ) );

        m_vertexTable.resize( m_vertexCount * m_entrySize, -1 );
        m_hashTable.resize( m_hashTableEntryCount * m_entrySize, -1 );
    }

    /**
     * @param data an array of indices.  The array size is set by the
     *    constructor's channelCount parameter.  The first entry in the array
     *    is the vertex index, which must less then the vertexCount
     *    specified in the constructor.
     * @return an id corresponding to the data.  This id is a serial number,
     *    which starts at zero.
     */
    index_t get_id( const index_t* data ) {
        assert( data );

        index_t* vertexEntry = get_vertex_entry( data[0] );
        if( is_occupied( vertexEntry ) ) {
            if( matches( vertexEntry, data ) ) {
                return vertexEntry[0];
            } else {
                index_t* hashEntry = find_hash_entry( data );
                if( is_occupied( hashEntry ) ) {
                    return hashEntry[0];
                } else {
                    return insert_hash_entry( hashEntry, data );
                }
            }
        } else {
            return insert_vertex_entry( vertexEntry, data );
        }
    }

    /**
     * @return The number of id's returned so far.  This is one greater than
     *    the value of the greatest id returned.
     */
    index_t get_id_count() const { return m_nextId; }

    /**
     *  Produce an index table, which indicates which index corresponds to a
     * given id and channel.
     *
     *  The table is indexed as outIndexTable[channel][id].
     */
    void get_index_table( std::vector<std::vector<index_t>>& outIndexTable ) {
        outIndexTable.resize( m_channelCount );
        BOOST_FOREACH( std::vector<index_t>& channelIndexTable, outIndexTable ) {
            channelIndexTable.resize( m_nextId );
        }
        copy_to_entry_table( m_vertexTable, outIndexTable );
        copy_to_entry_table( m_hashTable, outIndexTable );
    }

  private:
    bool is_occupied( const index_t* entry ) { return entry[0] != -1; }

    bool matches( const index_t* entry, const index_t* data ) {
        return memcmp( entry + 1, data, m_dataSizeInBytes ) == 0;
    }

    /**
     *  Copy data into entry.  Give entry a new id, and return the new id.
     */
    index_t populate( index_t* entry, const index_t* data ) {
        if( m_nextId == std::numeric_limits<index_t>::max() ) {
            throw std::runtime_error( "index_array_to_id::insert Error: reached maximum id value" );
        }
        index_t result = m_nextId++;
        entry[0] = result;
        memcpy( &entry[1], data, m_dataSizeInBytes );
        return result;
    }

    index_t insert_vertex_entry( index_t* vertexEntry, const index_t* data ) { return populate( vertexEntry, data ); }

    index_t insert_hash_entry( index_t* hashEntry, const index_t* data ) {
        const index_t result = populate( hashEntry, data );

        ++m_hashTableOccupancy;
        maybe_rehash();

        return result;
    }

    index_t* get_vertex_entry( const index_t vertexIndex ) {
        const std::size_t i = vertexIndex * m_entrySize;
        if( i >= m_vertexTable.size() ) {
            throw std::runtime_error( "index_array_to_id::find_vertex_entry Error: vertex index "
                                      "(" +
                                      boost::lexical_cast<std::string>( vertexIndex ) +
                                      ")"
                                      " is out of bounds" );
        }
        return &m_vertexTable[i];
    }

    boost::uint64_t hash( const index_t* entry ) const { return XXH64( entry, m_dataSizeInBytes, 0 ); }

    /**
     *  Attempt to find a hash table entry that matches data.
     *  If we find such an entry, return it.
     *  Otherwise, return an empty entry that is suitable for storing data.
     */
    index_t* find_hash_entry( std::vector<index_t>& table, const std::size_t tableEntryCount, const index_t* data ) {
        if( table.empty() ) {
            throw std::runtime_error( "index_array_to_id::find_hash_entry Error: table is empty" );
        }
        const std::size_t entrySize = m_entrySize;
        const std::size_t startEntryIndex = hash( data ) & ( tableEntryCount - 1 );

        // Start searching at the entry given by the hash function.
        // Search until the end of the table, and then loop back and start searching from index 0.
        typedef std::pair<std::size_t, std::size_t> range_t;
        const boost::array<range_t, 2> ranges = { range_t( startEntryIndex * entrySize, table.size() ),
                                                  range_t( 0, startEntryIndex * entrySize ) };

        BOOST_FOREACH( const range_t& range, ranges ) {
            for( std::size_t i = range.first, ie = range.second; i < ie; i += entrySize ) {
                index_t* entry = &table[i];

                if( !is_occupied( entry ) || matches( entry, data ) ) {
                    return entry;
                }
            }
        }

        throw std::runtime_error( "find_hash_entry Error: unable to find entry" );
    }

    index_t* find_hash_entry( const index_t* data ) {
        return find_hash_entry( m_hashTable, m_hashTableEntryCount, data );
    }

    static bool is_power_of_two( std::size_t n ) {
        if( n > 0 ) {
            return ( n & ( n - 1 ) ) == 0;
        }
        return false;
    }

    void maybe_rehash() {
        const std::size_t threshold = static_cast<std::size_t>( 0.75 * m_hashTableEntryCount );
        const bool needRehash = ( m_hashTableOccupancy >= threshold );
        if( needRehash ) {
            const std::size_t newHashTableEntryCount = 2 * m_hashTableEntryCount;
            assert( is_power_of_two( newHashTableEntryCount ) );

            const std::size_t newHashTableSize = newHashTableEntryCount * m_entrySize;
            assert( newHashTableSize > m_hashTable.size() );

            std::vector<index_t> newHashTable( newHashTableSize, -1 );

            rehash( m_hashTable, newHashTable );

            std::swap( m_hashTable, newHashTable );
            m_hashTableEntryCount = newHashTableEntryCount;
        }
    }

    void rehash( const std::vector<index_t>& in, std::vector<index_t>& out ) {
        const std::size_t outEntryCount = out.size() / m_entrySize;
        const std::size_t entrySizeInBytes = m_entrySize * sizeof( index_t );
        for( std::size_t i = 0, ie = in.size(); i < ie; i += m_entrySize ) {
            const index_t* inEntry = &in[i];
            if( is_occupied( inEntry ) ) {
                index_t* outEntry = find_hash_entry( out, outEntryCount, inEntry + 1 );
                memcpy( outEntry, inEntry, entrySizeInBytes );
            }
        }
    }

    void copy_to_entry_table( const std::vector<index_t>& in, std::vector<std::vector<index_t>>& out ) const {
        const std::size_t entrySize = m_entrySize;
        const std::size_t channelCount = m_channelCount;
        for( std::size_t i = 0, ie = in.size(); i < ie; i += entrySize ) {
            const index_t* entry = &in[i];
            const index_t id = entry[0];
            if( id != -1 ) {
                for( std::size_t channelIndex = 0; channelIndex < channelCount; ++channelIndex ) {
                    out[channelIndex][id] = entry[channelIndex + 1];
                }
            }
        }
    }

    // Size: ( m_channelCount + 1 ) * m_vertexCount
    // The entry for vertex index i begins at i * ( m_channelCount + 1 ).
    // An entry is composed of one index_t indicating the assigned id
    // (or -1 to indicate an empty entry), followed by m_channelCount
    // indices.
    std::vector<index_t> m_vertexTable;
    std::size_t m_entrySize;
    std::size_t m_dataSizeInBytes;
    // Similar to m_vertexTable.  Each hash table entry is of size
    // ( m_channelCount + 1 ), and
    std::vector<index_t> m_hashTable;
    std::size_t m_hashTableEntryCount;
    std::size_t m_hashTableOccupancy;
    index_t m_nextId;
    std::size_t m_channelCount;
    std::size_t m_vertexCount;
};

void copy_exploded_vertex_channel( frantic::geometry::polymesh3_ptr outMesh, const mesh_channel* inChannel,
                                   const boost::int32_t* outputToInputIndexMap ) {
    if( !outMesh ) {
        throw std::runtime_error( "copy_exploded_vertex_channel Error: output mesh is NULL" );
    }
    if( !inChannel ) {
        throw std::runtime_error( "copy_exploded_vertex_channel Error: input channel is NULL" );
    }

    const std::size_t outVertexCount = outMesh->vertex_count();
    const std::size_t inElementCount = inChannel->get_num_elements();
    const std::size_t elementSize = inChannel->get_element_size();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( outVertexCount * elementSize );
    char* out = buffer.begin();

    if( outVertexCount > 0 ) {
        if( !outputToInputIndexMap ) {
            throw std::runtime_error( "copy_exploded_vertex_channel Error: index map is NULL" );
        }

        for( std::size_t outVertexIndex = 0; outVertexIndex < outVertexCount; ++outVertexIndex ) {
            const boost::int32_t inIndex = outputToInputIndexMap[outVertexIndex];
            if( inIndex < 0 || static_cast<size_t>( inIndex ) >= inElementCount ) {
                throw std::runtime_error( "copy_exploded_vertex_channel Error: index " +
                                          boost::lexical_cast<std::string>( inIndex ) +
                                          " is out of range "
                                          "[0, " +
                                          boost::lexical_cast<std::string>( inElementCount ) + ")" );
            }
            inChannel->get_value( inIndex, out );
            out += elementSize;
        }
    }

    outMesh->add_vertex_channel( inChannel->get_name(), inChannel->get_data_type(), inChannel->get_data_arity(),
                                 buffer );
}

} // anonymous namespace

frantic::geometry::polymesh3_ptr explode_custom_faces( const frantic::geometry::mesh_interface* mesh,
                                                       const frantic::channels::channel_propagation_policy& cpp ) {
    typedef index_array_to_id::index_t index_t;
    typedef mesh_interface::mesh_channel_map channel_map_t;

    if( !mesh ) {
        throw std::runtime_error( "explode_custom_faces Error: mesh is NULL" );
    }
    if( !mesh->is_valid() ) {
        throw std::runtime_error( "explode_custom_faces Error: mesh is not valid" );
    }

    // vertex channels that are indexed using the mesh's faces
    std::vector<const mesh_channel*> simpleChannels;

    // vertex channels that use custom face indexing
    std::vector<const mesh_channel*> customChannels;

    const channel_map_t& vertexChannelMap = mesh->get_vertex_channels();
    for( channel_map_t::const_iterator i = vertexChannelMap.begin(), ie = vertexChannelMap.end(); i != ie; ++i ) {
        if( cpp.is_channel_included( i->first ) ) {
            const mesh_channel* channel = i->second;
            if( !channel ) {
                throw std::runtime_error( "explode_custom_faces Error: channel "
                                          "\"" +
                                          frantic::strings::to_string( i->first ) + "\" is NULL" );
            }
            const mesh_channel::channel_type channelType = channel->get_channel_type();
            if( channelType == mesh_channel::vertex ) {
                simpleChannels.push_back( channel );
            } else if( channelType == mesh_channel::face_vertex ) {
                customChannels.push_back( channel );
            } else {
                throw std::runtime_error( "explode_custom_faces Error: channel "
                                          "\"" +
                                          frantic::strings::to_string( i->first ) + "\" has unexpected type: " +
                                          boost::lexical_cast<std::string>( channelType ) );
            }
        }
    }

    frantic::geometry::polymesh3_ptr result;
    std::vector<std::vector<index_t>> idToIndexTable;
    std::size_t outVertexCount = 0;
    { // scope for faceVertexIds
        std::vector<index_t> faceVertexIds( get_face_vertex_count( mesh ) );
        std::size_t faceVertexIndex = 0;

        { // scope for table
            index_array_to_id table( mesh->get_num_verts(), 1 + customChannels.size() );

            std::vector<index_array_to_id::index_t> vertexChannelIndices( 1 + customChannels.size() );

            std::vector<std::size_t> meshFace( 1 );

            for( std::size_t faceIndex = 0, faceIndexEnd = mesh->get_num_faces(); faceIndex < faceIndexEnd;
                 ++faceIndex ) {
                const std::size_t cornerCount = mesh->get_num_face_verts( faceIndex );
                if( meshFace.size() < cornerCount ) {
                    meshFace.resize( cornerCount );
                }
                mesh->get_face_vert_indices( faceIndex, &meshFace[0] );

                for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
                    std::size_t channelIndex = 0;
                    vertexChannelIndices[channelIndex++] = boost::numeric_cast<index_t>( meshFace[cornerIndex] );
                    BOOST_FOREACH( const mesh_channel* channel, customChannels ) {
                        const std::size_t fvIndex = channel->get_fv_index( faceIndex, cornerIndex );
                        vertexChannelIndices[channelIndex++] = boost::numeric_cast<index_t>( fvIndex );
                    }
                    if( faceVertexIndex >= faceVertexIds.size() ) {
                        throw std::runtime_error( "explode_custom_faces Error: face vertex index is out of bounds" );
                    }
                    faceVertexIds[faceVertexIndex++] = table.get_id( &vertexChannelIndices[0] );
                }
            }

            table.get_index_table( idToIndexTable );
            outVertexCount = table.get_id_count();
        }

        // create the output polymesh

        polymesh3_builder builder;

        for( std::size_t vertexIndex = 0; vertexIndex < outVertexCount; ++vertexIndex ) {
            builder.add_vertex( mesh->get_vert( idToIndexTable[0][vertexIndex] ) );
        }

        std::vector<int> face;
        faceVertexIndex = 0;
        for( std::size_t faceIndex = 0, faceIndexEnd = mesh->get_num_faces(); faceIndex < faceIndexEnd; ++faceIndex ) {
            const std::size_t cornerCount = mesh->get_num_face_verts( faceIndex );
            face.resize( cornerCount );
            for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
                face[cornerIndex] = faceVertexIds[faceVertexIndex++];
            }
            builder.add_polygon( face );
        }

        result = builder.finalize();
    }

    assert( result );

    copy_face_channels( result, mesh, cpp );

    // copy vertex channels

    BOOST_FOREACH( const mesh_channel* channel, simpleChannels ) {
        const index_array_to_id::index_t* inputIndices = outVertexCount > 0 ? &idToIndexTable[0][0] : 0;

        copy_exploded_vertex_channel( result, channel, inputIndices );
    }

    for( std::size_t channelIndex = 0, channelCount = customChannels.size(); channelIndex < channelCount;
         ++channelIndex ) {
        const index_array_to_id::index_t* inputIndices = outVertexCount > 0 ? &idToIndexTable[1 + channelIndex][0] : 0;
        const mesh_channel* channel = customChannels[channelIndex];
        copy_exploded_vertex_channel( result, channel, inputIndices );
    }

    return result;
}

bool is_consistent_topology( polymesh3_ptr mesh1, polymesh3_ptr mesh2 ) {
    if( !mesh1->has_vertex_channel( _T("verts") ) || !mesh2->has_vertex_channel( _T("verts") ) ) {
        throw std::runtime_error(
            "is_consistent_topology Error: One or both of the provided meshes lacked a \"verts\" channel." );
    }

    if( mesh1->face_count() != mesh2->face_count() ) {
        return false;
    }

    if( mesh1->vertex_count() != mesh2->vertex_count() ) {
        return false;
    }

    frantic::geometry::polymesh3_const_vertex_accessor<vector3f> accMain =
        mesh1->get_const_vertex_accessor<vector3f>( _T("verts") );
    frantic::geometry::polymesh3_const_vertex_accessor<vector3f> accCheck =
        mesh2->get_const_vertex_accessor<vector3f>( _T("verts") );
    for( std::size_t i = 0, ie = mesh1->face_count(); i < ie; ++i ) {
        frantic::geometry::polymesh3_const_face_range rMain = accMain.get_face( i );
        frantic::geometry::polymesh3_const_face_range rCheck = accCheck.get_face( i );

        if( ( rMain.second - rMain.first ) != ( rCheck.second - rCheck.first ) ) {
            return false;
        }

        frantic::geometry::polymesh3_const_face_iterator itMain = rMain.first;
        frantic::geometry::polymesh3_const_face_iterator itCheck = rCheck.first;

        for( ; itMain != rMain.second; ++itMain, ++itCheck ) {
            if( *itMain != *itCheck ) {
                return false;
            }
        }
    }

    return true;
}

bool is_equal( const const_polymesh3_ptr& mesh1, const const_polymesh3_ptr& mesh2 ) {
    mesh_interface::ptr_type a( polymesh3_interface::create_const_instance( mesh1 ).release() );
    mesh_interface::ptr_type b( polymesh3_interface::create_const_instance( mesh2 ).release() );

    return is_equal( a, b );
}

frantic::graphics::boundbox3f compute_boundbox( frantic::geometry::const_polymesh3_ptr mesh ) {
    frantic::graphics::boundbox3f result( vector3f( 1.f ), vector3f( -1.f ) );

    if( mesh ) {
        polymesh3_const_vertex_accessor<vector3f> geomAcc = mesh->get_const_vertex_accessor<vector3f>( _T("verts") );

        for( std::size_t i = 0, iEnd = geomAcc.vertex_count(); i < iEnd; ++i ) {
            result += geomAcc.get_vertex( i );
        }
    }

    return result;
}

bool is_valid_maya_face( const frantic::geometry::polymesh3_const_face_range& faceRange ) {
    const int degree = (int)( faceRange.second - faceRange.first );
    frantic::geometry::polymesh3_const_face_iterator f = faceRange.first;

    // can't repeat f[0], f[1]
    {
        const int tabu0 = f[0];
        const int tabu1 = f[1];
        for( int i = 1, ie = degree - 1; i < ie; ++i ) {
            if( f[i] == tabu0 && f[i + 1] == tabu1 ) {
                return false;
            }
        }
    }
    // can't repeat f[i], f[i-1] unless f[i] appeared earlier in the sequence
    for( int i = 1, ie = degree - 1; i < ie; ++i ) {
        const int tabu0 = f[i];
        const int tabu1 = f[i - 1];

        bool found = false;
        for( int j = 0; j < i; ++j ) {
            if( f[j] == tabu0 ) {
                found = true;
                break;
            }
        }

        if( !found ) {
            for( int j = i, je = degree - 1; j < je; ++j ) {
                if( f[j] == tabu0 && f[j + 1] == tabu1 ) {
                    return false;
                }
            }
        }
    }
    // can't repeat f[i], f[i] unless f[i] is f[0]
    for( int i = 1, ie = degree - 1; i < ie; ++i ) {
        const int tabu = f[i];

        if( tabu == f[0] ) {
            continue;
        }

        for( int j = i + 1, je = degree - 1; j < je; ++j ) {
            if( f[j] == tabu && f[j + 1] == tabu ) {
                return false;
            }
        }
    }
    return true;
}
} // namespace geometry
} // namespace frantic
