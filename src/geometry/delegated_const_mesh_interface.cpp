// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/delegated_const_mesh_interface.hpp>

#include <frantic/channels/channel_propagation_policy.hpp>

namespace frantic {
namespace geometry {

namespace {

class delegated_const_mesh_channel : public mesh_channel {
  public:
    delegated_const_mesh_channel( const mesh_channel* channel );
    void get_value( std::size_t index, void* outValue ) const;
    void set_value( std::size_t index, const void* value ) const;
    std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const;
    std::size_t get_num_face_verts( std::size_t faceIndex ) const;

  private:
    delegated_const_mesh_channel(); // not implemented

    const mesh_channel* m_channel;
};

} // anonymous namespace

//
// delegated_const_mesh_interface
//

delegated_const_mesh_interface::delegated_const_mesh_interface( const mesh_interface* meshInterface )
    : m_mesh( meshInterface ) {
    if( !m_mesh ) {
        throw std::runtime_error( "delegated_const_mesh_interface() Mesh is NULL" );
    }
}

bool delegated_const_mesh_interface::is_valid() const { return m_mesh != NULL; }

bool delegated_const_mesh_interface::request_channel( const frantic::tstring& channelName, bool vertexChannel,
                                                      bool forOutput, bool throwOnError ) {
    if( forOutput ) {
        if( throwOnError ) {
            throw std::runtime_error( "delegated_const_mesh_interface::request_channel() "
                                      "Cannot provide writable channel: "
                                      "\"" +
                                      frantic::strings::to_string( channelName ) + "\"" );
        }
        return false;
    }

    const mesh_channel* ch = 0;

    if( vertexChannel ) {
        ch = this->get_vertex_channels().get_channel( channelName );
    } else {
        ch = this->get_face_channels().get_channel( channelName );
    }

    if( ch ) {
        return true;
    } else {
        if( throwOnError ) {
            throw std::runtime_error( "delegated_const_mesh_interface::request_channel() "
                                      "Cannot add channel: "
                                      "\"" +
                                      frantic::strings::to_string( channelName ) + "\"" );
        }
        return false;
    }
}

std::size_t delegated_const_mesh_interface::get_num_verts() const {
    assert( m_mesh );
    return m_mesh->get_num_verts();
}

void delegated_const_mesh_interface::get_vert( std::size_t index, float ( &outValues )[3] ) const {
    assert( m_mesh );
    m_mesh->get_vert( index, outValues );
}

std::size_t delegated_const_mesh_interface::get_num_faces() const {
    assert( m_mesh );
    return m_mesh->get_num_faces();
}

std::size_t delegated_const_mesh_interface::get_num_face_verts( std::size_t faceIndex ) const {
    assert( m_mesh );
    return m_mesh->get_num_face_verts( faceIndex );
}

std::size_t delegated_const_mesh_interface::get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
    assert( m_mesh );
    return m_mesh->get_face_vert_index( faceIndex, fvertIndex );
}

void delegated_const_mesh_interface::get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const {
    assert( m_mesh );
    m_mesh->get_face_vert_indices( faceIndex, outValues );
}

void delegated_const_mesh_interface::get_face_verts( std::size_t faceIndex, float outValues[][3] ) const {
    assert( m_mesh );
    m_mesh->get_face_verts( faceIndex, outValues );
}

std::size_t delegated_const_mesh_interface::get_num_elements() const {
    assert( m_mesh );
    return m_mesh->get_num_elements();
}

std::size_t delegated_const_mesh_interface::get_face_element_index( std::size_t faceIndex ) const {
    assert( m_mesh );
    return m_mesh->get_face_element_index( faceIndex );
}

void delegated_const_mesh_interface::init_adjacency() {
    assert( m_mesh );
    if( !m_mesh->has_adjacency() ) {
        throw std::runtime_error( "delegated_const_mesh_interface::init_adjacency() Not implemented" );
    }
}

bool delegated_const_mesh_interface::has_adjacency() const {
    assert( m_mesh );
    return m_mesh->has_adjacency();
}

bool delegated_const_mesh_interface::init_vertex_iterator( frantic::geometry::vertex_iterator& vIt,
                                                           std::size_t vertexIndex ) const {
    assert( m_mesh );
    return m_mesh->init_vertex_iterator( vIt, vertexIndex );
}

bool delegated_const_mesh_interface::advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const {
    assert( m_mesh );
    return m_mesh->advance_vertex_iterator( vIt );
}

std::size_t delegated_const_mesh_interface::get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const {
    assert( m_mesh );
    return m_mesh->get_edge_endpoint( vIt );
}

std::size_t delegated_const_mesh_interface::get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const {
    assert( m_mesh );
    return m_mesh->get_edge_left_face( vIt );
}

std::size_t delegated_const_mesh_interface::get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const {
    assert( m_mesh );
    return m_mesh->get_edge_right_face( vIt );
}

bool delegated_const_mesh_interface::is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const {
    assert( m_mesh );
    return m_mesh->is_edge_visible( vIt );
}

bool delegated_const_mesh_interface::is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const {
    assert( m_mesh );
    return m_mesh->is_edge_boundary( vIt );
}

void delegated_const_mesh_interface::init_face_iterator( frantic::geometry::face_iterator& fIt,
                                                         std::size_t faceIndex ) const {
    assert( m_mesh );
    m_mesh->init_face_iterator( fIt, faceIndex );
}

bool delegated_const_mesh_interface::advance_face_iterator( frantic::geometry::face_iterator& fIt ) const {
    assert( m_mesh );
    return m_mesh->advance_face_iterator( fIt );
}

std::size_t delegated_const_mesh_interface::get_face_neighbor( frantic::geometry::face_iterator& fIt ) const {
    assert( m_mesh );
    return m_mesh->get_face_neighbor( fIt );
}

std::size_t delegated_const_mesh_interface::get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const {
    assert( m_mesh );
    return m_mesh->get_face_next_vertex( fIt );
}

std::size_t delegated_const_mesh_interface::get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const {
    assert( m_mesh );
    return m_mesh->get_face_prev_vertex( fIt );
}

void delegated_const_mesh_interface::append_delegate_vertex_channel( const frantic::tstring& channelName ) {
    const mesh_channel* delegateChannel = m_mesh->get_vertex_channels().get_channel( channelName );
    append_delegate_channel( delegateChannel, channelName );
}

void delegated_const_mesh_interface::append_delegate_face_channel( const frantic::tstring& channelName ) {
    const mesh_channel* delegateChannel = m_mesh->get_face_channels().get_channel( channelName );
    append_delegate_channel( delegateChannel, channelName );
}

void delegated_const_mesh_interface::append_delegate_channel( const frantic::geometry::mesh_channel* delegateChannel,
                                                              const frantic::tstring& channelNameForErrorMessage ) {
    if( !delegateChannel ) {
        throw std::runtime_error( "delegated_const_mesh_interface::append_delegate_channel() Channel "
                                  "\"" +
                                  frantic::strings::to_string( channelNameForErrorMessage ) +
                                  "\" "
                                  " is NULL" );
    }

    std::unique_ptr<delegated_const_mesh_channel> myChannel( new delegated_const_mesh_channel( delegateChannel ) );

    const mesh_channel::channel_type channelType = myChannel->get_channel_type();
    if( mesh_channel::is_stored_at_vertex( channelType ) ) {
        append_vertex_channel( std::move( myChannel ) );
    } else if( channelType == mesh_channel::face ) {
        append_face_channel( std::move( myChannel ) );
    } else {
        throw std::runtime_error( "delegated_const_mesh_interface::add_delegate_channel() Unknown channel type: " +
                                  boost::lexical_cast<std::string>( channelType ) );
    }
}

void delegated_const_mesh_interface::append_delegate_face_channels(
    const frantic::channels::channel_propagation_policy& cpp ) {
    append_delegate_channels( m_mesh->get_face_channels(), cpp );
}

void delegated_const_mesh_interface::append_delegate_vertex_channels(
    const frantic::channels::channel_propagation_policy& cpp ) {
    append_delegate_channels( m_mesh->get_vertex_channels(), cpp );
}

void delegated_const_mesh_interface::append_delegate_channels(
    const frantic::channels::channel_propagation_policy& cpp ) {
    append_delegate_vertex_channels( cpp );
    append_delegate_face_channels( cpp );
}

void delegated_const_mesh_interface::append_delegate_channels(
    const frantic::geometry::mesh_interface::mesh_channel_map& channels,
    const frantic::channels::channel_propagation_policy& cpp ) {
    typedef frantic::geometry::mesh_interface::mesh_channel_map channel_map_t;

    for( channel_map_t::const_iterator i = channels.begin(), ie = channels.end(); i != ie; ++i ) {
        if( cpp.is_channel_included( i->first ) ) {
            const mesh_channel::channel_type channelType = i->second->get_channel_type();
            if( mesh_channel::is_stored_at_vertex( channelType ) ) {
                append_delegate_vertex_channel( i->first );
            } else if( channelType == mesh_channel::face ) {
                append_delegate_face_channel( i->first );
            } else {
                throw std::runtime_error( "create_vertex_normal_channel_mesh_interface_impl::add_delegate_channels "
                                          "Error: Unknown channel type: " +
                                          boost::lexical_cast<std::string>( channelType ) );
            }
        }
    }
}

//
// delegated_const_mesh_channel
//

delegated_const_mesh_channel::delegated_const_mesh_channel( const mesh_channel* ch )
    : mesh_channel( ch->get_name(), ch->get_channel_type(), ch->get_data_type(), ch->get_data_arity(),
                    ch->get_num_elements(), ch->get_num_faces(), true )
    , m_channel( ch ) {}

void delegated_const_mesh_channel::get_value( std::size_t index, void* outValue ) const {
    assert( m_channel );
    m_channel->get_value( index, outValue );
}

void delegated_const_mesh_channel::set_value( std::size_t /*index*/, const void* /*value*/ ) const {
    // pass
}

std::size_t delegated_const_mesh_channel::get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
    assert( m_channel );
    return m_channel->get_fv_index( faceIndex, fvertIndex );
}

std::size_t delegated_const_mesh_channel::get_num_face_verts( std::size_t faceIndex ) const {
    assert( m_channel );
    return m_channel->get_num_face_verts( faceIndex );
}

} // namespace geometry
} // namespace frantic
