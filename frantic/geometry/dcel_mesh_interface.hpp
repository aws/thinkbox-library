// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/mesh_interface.hpp>

#include <boost/config.hpp>
#include <boost/move/move.hpp>

namespace frantic {
namespace geometry {

class dcel_mesh_interface : public mesh_interface {
  public:
#if( defined( __linux__ ) || defined( __APPLE__ ) ) && !defined( BOOST_NO_CXX11_HDR_SMART_PTR )
    // For the version of Boost we use, Boost.Move doesn't seem to work correctly in C++11
    dcel_mesh_interface( const dcel* impl, bool release );
#endif

    dcel_mesh_interface( const dcel* impl, BOOST_RV_REF( std::vector<frantic::graphics::vector3f> ) vertices )
        : m_impl( impl )
        , m_hasDcelOwnership( false )
        , m_vertices( NULL )
        , m_hasVerticesOwnership( true ) {
        m_vertices = new std::vector<frantic::graphics::vector3f>( boost::move( vertices ) );
    }

    dcel_mesh_interface( const dcel* impl, const std::vector<frantic::graphics::vector3f>* vertices = NULL );

    dcel_mesh_interface( BOOST_RV_REF( dcel ) impl, const std::vector<frantic::graphics::vector3f>* vertices = NULL )
        : m_impl( NULL )
        , m_hasDcelOwnership( true )
        , m_vertices( vertices )
        , m_hasVerticesOwnership( false ) {
        m_impl = new dcel( boost::move( impl ) );
    }

    dcel_mesh_interface( BOOST_RV_REF( dcel ) impl, BOOST_RV_REF( std::vector<frantic::graphics::vector3f> ) vertices )
        : m_impl( NULL )
        , m_hasDcelOwnership( true )
        , m_vertices( NULL )
        , m_hasVerticesOwnership( true ) {
        m_impl = new dcel( boost::move( impl ) );
        m_vertices = new std::vector<frantic::graphics::vector3f>( boost::move( vertices ) );
    }

    ~dcel_mesh_interface();

    virtual bool is_valid() const;

    virtual bool request_channel( const frantic::tstring& channelName, bool vertexChannel, bool forOutput,
                                  bool throwOnError = true );

    virtual std::size_t get_num_verts() const;

    virtual void get_vert( std::size_t index, float ( &outValues )[3] ) const;

    virtual std::size_t get_num_faces() const;

    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const;

    virtual std::size_t get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const;

    virtual void get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const;

    virtual void get_face_verts( std::size_t faceIndex, float outValues[][3] ) const;

    virtual std::size_t get_num_elements() const;

    virtual std::size_t get_face_element_index( std::size_t faceIndex ) const;

    virtual vertex_adjacency_interface& get_vertex_adjacency();

    virtual face_adjacency_interface& get_face_adjacency();

    virtual void init_adjacency();

    virtual bool has_adjacency() const;

    virtual bool init_vertex_iterator( frantic::geometry::vertex_iterator& vIt, std::size_t vertexIndex ) const;

    virtual bool advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const;

    virtual std::size_t get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const;

    virtual std::size_t get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const;

    virtual std::size_t get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const;

    virtual bool is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const;

    virtual bool is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const;

    virtual void init_face_iterator( frantic::geometry::face_iterator& fIt, std::size_t faceIndex ) const;

    virtual bool advance_face_iterator( frantic::geometry::face_iterator& fIt ) const;

    virtual std::size_t get_face_neighbor( frantic::geometry::face_iterator& fIt ) const;

    virtual std::size_t get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const;

    virtual std::size_t get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const;

    virtual const dcel& get_dcel() const;
    virtual const std::vector<frantic::graphics::vector3f>& vertices_ref() const;
    virtual bool has_dcel_ownership() const;
    virtual bool has_vertices_ownership() const;

  private:
    size_t get_as_interface_face( size_t faceId ) const;

    const dcel* m_impl;
    bool m_hasDcelOwnership;

    const std::vector<frantic::graphics::vector3f>* m_vertices;
    bool m_hasVerticesOwnership;
};

} // namespace geometry
} // namespace frantic
