// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/dcel_mesh_interface.hpp>

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/dcel_iterators.hpp>

#include <boost/config.hpp>
#include <boost/foreach.hpp>
#include <frantic/misc/iterator.hpp>

namespace frantic {
namespace geometry {

namespace {

BOOST_STATIC_ASSERT( sizeof( dcel::const_halfedge_handle ) <= frantic::geometry::detail::iterator_storage_size );
BOOST_STATIC_ASSERT( sizeof( dcel::const_halfedge_handle ) <= frantic::geometry::detail::iterator_storage_size );

dcel::const_halfedge_handle& to_dcel_vertex_iterator( vertex_iterator& vIt ) {
    return *static_cast<dcel::const_halfedge_handle*>( vIt.m_data.address() );
}

dcel::const_halfedge_handle& to_dcel_face_iterator( face_iterator& vIt ) {
    return *static_cast<dcel::const_halfedge_handle*>( vIt.m_data.address() );
}

} // namespace

#if( defined( __linux__ ) || defined( __APPLE__ ) ) && !defined( BOOST_NO_CXX11_HDR_SMART_PTR )
dcel_mesh_interface::dcel_mesh_interface( const dcel* impl, bool release )
    : m_impl( impl )
    , m_hasDcelOwnership( release )
    , m_vertices( NULL )
    , m_hasVerticesOwnership( false ) {}
#endif

dcel_mesh_interface::dcel_mesh_interface( const dcel* impl, const std::vector<frantic::graphics::vector3f>* vertices )
    : m_impl( impl )
    , m_hasDcelOwnership( false )
    , m_vertices( vertices )
    , m_hasVerticesOwnership( false ) {}

dcel_mesh_interface::~dcel_mesh_interface() {
    if( m_hasDcelOwnership ) {
        delete m_impl;
    }

    if( m_hasVerticesOwnership ) {
        delete m_vertices;
    }
}

bool dcel_mesh_interface::is_valid() const { return m_impl != NULL; }

bool dcel_mesh_interface::request_channel( const frantic::tstring& /*channelName*/, bool /*vertexChannel*/,
                                           bool /*forOutput*/, bool throwOnError ) {
    if( throwOnError ) {
        throw std::runtime_error(
            "dcel_mesh_interface::request_channel : No channels are supported on dcel_mesh_interface." );
    }
    return false;
}

std::size_t dcel_mesh_interface::get_num_verts() const { return m_impl->vertex_count(); }

void dcel_mesh_interface::get_vert( std::size_t index, float ( &outValues )[3] ) const {
    for( size_t i = 0; i < 3; ++i ) {
        outValues[i] = m_vertices->at( index )[int( i )];
    }
}

std::size_t dcel_mesh_interface::get_num_faces() const { return m_impl->face_count(); }

std::size_t dcel_mesh_interface::get_num_face_verts( std::size_t faceIndex ) const {
    if( m_impl->is_triangle_mesh() ) {
        return 3;
    } else {
        return std::distance( dcel_face_cycle_begin( m_impl->get_face_halfedge( faceIndex ) ),
                              dcel_face_cycle_end( m_impl->get_face_halfedge( faceIndex ) ) );
    }
}

std::size_t dcel_mesh_interface::get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
    size_t current = 0;

    BOOST_FOREACH( dcel::const_halfedge_handle handle,
                   dcel_face_cycle_range( m_impl->get_face_halfedge( faceIndex ) ) ) {
        if( current == fvertIndex ) {
            return handle.source_vertex();
        }
        ++current;
    }

    throw std::runtime_error( "dcel_mesh_interface::get_face_vert_index : face index " +
                              boost::lexical_cast<std::string>( fvertIndex ) + " out of range." );
}

void dcel_mesh_interface::get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const {
    size_t current = 0;

    BOOST_FOREACH( dcel::const_halfedge_handle handle,
                   dcel_face_cycle_range( m_impl->get_face_halfedge( faceIndex ) ) ) {
        outValues[current] = handle.source_vertex();
        ++current;
    }
}

void dcel_mesh_interface::get_face_verts( std::size_t faceIndex, float outValues[][3] ) const {
    size_t current = 0;

    BOOST_FOREACH( dcel::const_halfedge_handle handle,
                   dcel_face_cycle_range( m_impl->get_face_halfedge( faceIndex ) ) ) {
        get_vert( handle.source_vertex(), outValues[current] );
        ++current;
    }
}

// TODO: maybe change this to allow caching a face element set
// (& then possibly generate a mesh_segmentation from it?)
std::size_t dcel_mesh_interface::get_num_elements() const { return 1; }

std::size_t dcel_mesh_interface::get_face_element_index( std::size_t /*faceIndex*/ ) const { return 0; }

// Can we drop these methods from the interface?
vertex_adjacency_interface& dcel_mesh_interface::get_vertex_adjacency() {
    throw std::runtime_error( "dcel_mesh_interface::get_vertex_adjacency : This method is not usable." );
}

face_adjacency_interface& dcel_mesh_interface::get_face_adjacency() {
    throw std::runtime_error( "dcel_mesh_interface::get_face_adjacency : This method is not usable." );
}

void dcel_mesh_interface::init_adjacency() {
    // the dcel is an adjacency structure already
}

bool dcel_mesh_interface::has_adjacency() const { return true; }

bool dcel_mesh_interface::init_vertex_iterator( frantic::geometry::vertex_iterator& vIt,
                                                std::size_t vertexIndex ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_vertex_iterator( vIt );
    edgeHandle = m_impl->get_vertex_halfedge( vertexIndex );
    return m_impl->has_vertex_halfedge( vertexIndex );
}

bool dcel_mesh_interface::advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_vertex_iterator( vIt );
    edgeHandle = edgeHandle.vertex_next();
    return edgeHandle != m_impl->get_vertex_halfedge( edgeHandle.target_vertex() );
}

std::size_t dcel_mesh_interface::get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_vertex_iterator( vIt );
    return edgeHandle.source_vertex();
}

std::size_t dcel_mesh_interface::get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_vertex_iterator( vIt );
    return get_as_interface_face( edgeHandle.opposite_face() );
}

std::size_t dcel_mesh_interface::get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_vertex_iterator( vIt );
    return get_as_interface_face( edgeHandle.current_face() );
}

bool dcel_mesh_interface::is_edge_visible( frantic::geometry::vertex_iterator& /*vIt*/ ) const { return true; }

bool dcel_mesh_interface::is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_vertex_iterator( vIt );
    return edgeHandle.is_boundary_edge();
}

void dcel_mesh_interface::init_face_iterator( frantic::geometry::face_iterator& fIt, std::size_t faceIndex ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_face_iterator( fIt );
    edgeHandle = m_impl->get_face_halfedge( faceIndex );
}

bool dcel_mesh_interface::advance_face_iterator( frantic::geometry::face_iterator& fIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_face_iterator( fIt );
    edgeHandle = edgeHandle.face_next();
    return edgeHandle != m_impl->get_face_halfedge( edgeHandle.current_face() );
}

std::size_t dcel_mesh_interface::get_face_neighbor( frantic::geometry::face_iterator& fIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_face_iterator( fIt );
    return get_as_interface_face( edgeHandle.opposite_face() );
}

std::size_t dcel_mesh_interface::get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_face_iterator( fIt );
    return edgeHandle.source_vertex();
}

std::size_t dcel_mesh_interface::get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const {
    dcel::const_halfedge_handle& edgeHandle = to_dcel_face_iterator( fIt );
    return edgeHandle.target_vertex();
}

const dcel& dcel_mesh_interface::get_dcel() const { return *m_impl; }

const std::vector<frantic::graphics::vector3f>& dcel_mesh_interface::vertices_ref() const { return *m_vertices; }

bool dcel_mesh_interface::has_dcel_ownership() const { return m_hasDcelOwnership; }

bool dcel_mesh_interface::has_vertices_ownership() const { return m_hasVerticesOwnership; }

size_t dcel_mesh_interface::get_as_interface_face( size_t faceId ) const {
    if( faceId >= m_impl->face_count() ) {
        return HOLE_INDEX;
    } else {
        return size_t( faceId );
    }
}

} // namespace geometry
} // namespace frantic
