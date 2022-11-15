// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/xmesh_standard_mesh_interface.hpp>

#include <boost/optional.hpp>
#include <boost/regex.hpp>

#include <frantic/geometry/delegated_const_mesh_interface.hpp>

using namespace frantic::channels;

namespace frantic {
namespace geometry {

namespace {

boost::optional<std::pair<data_type_t, std::size_t>>
try_get_xmesh_standard_vertex_channel_type( const frantic::tstring& channel ) {
    typedef std::pair<data_type_t, std::size_t> result_t;

    if( channel == _T("Color") ) {
        return result_t( data_type_float32, 3 );
    } else if( channel == _T("TextureCoord") ) {
        return result_t( data_type_float32, 3 );
    } else if( channel == _T("Normal") ) {
        return result_t( data_type_float32, 3 );
    } else if( channel == _T("Selection") ) {
        return result_t( data_type_float32, 1 );
    } else if( channel == _T("Velocity") ) {
        return result_t( data_type_float32, 3 );
    } else {
        // Defines the set of strings from 'Mapping2' through 'Mapping99'.
        boost::basic_regex<frantic::tchar> expression( _T("Mapping([2-9]|[1-9][0-9])") );
        if( boost::regex_match( channel, expression ) ) {
            return result_t( data_type_float32, 3 );
        } else {
            return boost::optional<result_t>();
        }
    }
}

boost::optional<std::pair<data_type_t, std::size_t>>
try_get_xmesh_standard_face_channel_type( const frantic::tstring& channel ) {
    typedef std::pair<data_type_t, std::size_t> result_t;

    if( channel == _T("MaterialID") ) {
        return result_t( data_type_uint16, 1 );
    } else if( channel == _T("SmoothingGroup") ) {
        return result_t( data_type_int32, 1 );
    } else if( channel == _T("FaceSelection") ) {
        return result_t( data_type_int32, 1 );
    } else {
        return boost::optional<result_t>();
    }
}

boost::optional<std::pair<data_type_t, std::size_t>>
try_get_xmesh_standard_channel_type( const frantic::tstring& channel ) {
    typedef std::pair<data_type_t, std::size_t> result_t;

    const boost::optional<result_t> vertexChannelType = try_get_xmesh_standard_vertex_channel_type( channel );

    if( vertexChannelType ) {
        return vertexChannelType;
    } else {
        return try_get_xmesh_standard_face_channel_type( channel );
    }
}

std::pair<data_type_t, std::size_t> get_xmesh_standard_channel_type( const frantic::tstring& channel ) {
    typedef std::pair<data_type_t, std::size_t> result_t;

    const boost::optional<result_t> standardType = try_get_xmesh_standard_channel_type( channel );

    if( !standardType ) {
        throw std::runtime_error( "get_xmesh_standard_channel_type Error: "
                                  "channel is not a standard xmesh channel: "
                                  "\"" +
                                  frantic::strings::to_string( channel ) + "\"" );
    }

    return *standardType;
}

frantic::channels::data_type_t get_xmesh_standard_data_type( const frantic::tstring& channel ) {
    return get_xmesh_standard_channel_type( channel ).first;
}

std::size_t get_xmesh_standard_arity( const frantic::tstring& channel ) {
    return get_xmesh_standard_channel_type( channel ).second;
}

bool is_xmesh_standard_vertex_channel_name( const frantic::tstring& channel ) {
    return (bool)try_get_xmesh_standard_vertex_channel_type( channel );
}

bool is_xmesh_standard_face_channel_name( const frantic::tstring& channel ) {
    return (bool)try_get_xmesh_standard_face_channel_type( channel );
}

class xmesh_standard_mesh_channel : public mesh_channel {
  public:
    xmesh_standard_mesh_channel( const mesh_channel* channel );
    void get_value( std::size_t index, void* outValue ) const;
    void set_value( std::size_t index, const void* value ) const;
    std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const;
    std::size_t get_num_face_verts( std::size_t faceIndex ) const;

  private:
    const mesh_channel* m_channel;
    channel_type_convertor_function_t m_convert;
    std::size_t m_convertArity;
};

xmesh_standard_mesh_channel::xmesh_standard_mesh_channel( const mesh_channel* channel )
    : mesh_channel(
          channel->get_name(), channel->get_channel_type(), get_xmesh_standard_data_type( channel->get_name() ),
          get_xmesh_standard_arity( channel->get_name() ), channel->get_num_elements(), channel->get_num_faces() )
    , m_channel( channel ) {
    if( !m_channel ) {
        throw std::runtime_error( "xmesh_standard_mesh_channel Error: channel is NULL" );
    }

    if( is_stored_at_vertex( get_channel_type() ) ) {
        if( !is_xmesh_standard_vertex_channel_name( get_name() ) ) {
            throw std::runtime_error( "xmesh_standard_mesh_channel Internal Error: "
                                      "got vertex channel "
                                      "\"" +
                                      frantic::strings::to_string( get_name() ) +
                                      "\", "
                                      "but this is not a standard vertex channel" );
        }
    } else if( get_channel_type() == frantic::geometry::mesh_channel::face ) {
        if( !is_xmesh_standard_face_channel_name( get_name() ) ) {
            throw std::runtime_error( "xmesh_standard_mesh_channel Internal Error: "
                                      "got face channel "
                                      "\"" +
                                      frantic::strings::to_string( get_name() ) +
                                      "\", "
                                      "but this is not a standard face channel" );
        }
    } else {
        throw std::runtime_error( "xmesh_standard_mesh_channel Internal Error: "
                                  "got channel "
                                  "\"" +
                                  frantic::strings::to_string( get_name() ) +
                                  "\" "
                                  "with unhandled channel type: " +
                                  boost::lexical_cast<std::string>( get_channel_type() ) );
    }

    m_convert =
        frantic::channels::get_channel_type_convertor_function( channel->get_data_type(), get_data_type(), get_name() );
    m_convertArity = std::min( get_data_arity(), channel->get_data_arity() );
}

void xmesh_standard_mesh_channel::get_value( std::size_t index, void* outValue ) const {
    assert( m_channel );
    assert( outValue );

    memset( outValue, 0, get_element_size() );

    char* temp = reinterpret_cast<char*>( alloca( m_channel->get_element_size() ) );

    m_channel->get_value( index, temp );
    m_convert( reinterpret_cast<char*>( outValue ), temp, m_convertArity );
}

void xmesh_standard_mesh_channel::set_value( std::size_t /*index*/, const void* /*value*/ ) const {
    // pass
}

std::size_t xmesh_standard_mesh_channel::get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
    assert( m_channel );
    return m_channel->get_fv_index( faceIndex, fvertIndex );
}

std::size_t xmesh_standard_mesh_channel::get_num_face_verts( std::size_t faceIndex ) const {
    assert( m_channel );
    return m_channel->get_num_face_verts( faceIndex );
}

class xmesh_standard_mesh_interface : public delegated_const_mesh_interface {
  public:
    xmesh_standard_mesh_interface( const mesh_interface* meshInterface )
        : delegated_const_mesh_interface( meshInterface ) {
        if( !meshInterface ) {
            throw std::runtime_error( "xmesh_standard_mesh_interface Error: mesh interface is NULL" );
        }

        append_all_xmesh_standard_delegate_channels();
    }

  private:
    void append_all_xmesh_standard_delegate_channels() {
        append_xmesh_standard_delegate_channels( m_mesh->get_vertex_channels() );
        append_xmesh_standard_delegate_channels( m_mesh->get_face_channels() );
    }

    void
    append_xmesh_standard_delegate_channels( const frantic::geometry::mesh_interface::mesh_channel_map& channels ) {
        typedef frantic::geometry::mesh_interface::mesh_channel_map channel_map_t;

        for( channel_map_t::const_iterator i = channels.begin(), ie = channels.end(); i != ie; ++i ) {
            const mesh_channel::channel_type channelType = i->second->get_channel_type();
            if( mesh_channel::is_stored_at_vertex( channelType ) ) {
                if( is_xmesh_standard_vertex_channel_name( i->first ) ) {
                    std::unique_ptr<mesh_channel> channel( new xmesh_standard_mesh_channel( i->second ) );
                    append_vertex_channel( std::move( channel ) );
                }
            } else if( channelType == mesh_channel::face ) {
                if( is_xmesh_standard_face_channel_name( i->first ) ) {
                    std::unique_ptr<mesh_channel> channel( new xmesh_standard_mesh_channel( i->second ) );
                    append_face_channel( std::move( channel ) );
                }
            }
        }
    }
};

} // anonymous namespace

std::unique_ptr<mesh_interface> create_xmesh_standard_mesh_interface( const mesh_interface* meshInterface ) {
    return std::unique_ptr<mesh_interface>( new xmesh_standard_mesh_interface( meshInterface ) );
}

} // namespace geometry
} // namespace frantic
