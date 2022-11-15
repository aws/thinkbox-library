// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/trimesh3_interface.hpp>

#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/dcel_mesh_interface.hpp>

#include <boost/cstdint.hpp>

#include <frantic/misc/hybrid_assert.hpp>

namespace frantic {
namespace geometry {

class trimesh3_interface_impl : public trimesh3_interface {
  protected:
    virtual mesh_channel* create_vertex_channel( const frantic::tstring& channelName,
                                                 frantic::channels::data_type_t dataType, size_t arity );

    virtual mesh_channel* create_vertex_channel_custom_faces( const frantic::tstring& channelName,
                                                              frantic::channels::data_type_t dataType, size_t arity,
                                                              size_t vertexCount );

    virtual void destroy_vertex_channel( mesh_channel* channel );

    virtual mesh_channel* create_face_channel( const frantic::tstring& channelName,
                                               frantic::channels::data_type_t dataType, size_t arity );

    virtual void destroy_face_channel( mesh_channel* channel );

  public:
    trimesh3_interface_impl();

    virtual ~trimesh3_interface_impl() {}

    void set_mesh_pointer( const trimesh3* mesh );

    void set_mesh( BOOST_RV_REF( trimesh3 ) mesh );

    virtual const trimesh3& get_trimesh() const;

    virtual bool is_valid() const;

    virtual bool request_channel( const frantic::tstring& channelName, bool vertexChannel, bool forOutput,
                                  bool throwOnError = true );

    virtual std::size_t get_num_verts() const;

    virtual void get_vert( std::size_t index, float ( &outValues )[3] ) const;

    virtual bool is_read_only() const;

    virtual void set_vert( std::size_t index, const float* v );

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

    virtual std::size_t get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const;

    virtual std::size_t get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const;

  private:
    frantic::geometry::trimesh3* m_mesh;
    frantic::geometry::trimesh3 m_ownedMesh;

    std::unique_ptr<dcel_mesh_interface> m_adjacencyDelegate;

  private:
    void apply_underlying_channels();

    void init_iterator_data();
};

trimesh3_interface_impl::trimesh3_interface_impl()
    : m_mesh( NULL ) {}

std::unique_ptr<trimesh3_interface> trimesh3_interface::create_instance( const trimesh3* theMesh ) {
    std::unique_ptr<trimesh3_interface_impl> result( new trimesh3_interface_impl );
    result->set_mesh_pointer( theMesh );
    return std::unique_ptr<trimesh3_interface>( result.release() );
}

std::unique_ptr<trimesh3_interface> trimesh3_interface::create_instance( BOOST_RV_REF( trimesh3 ) theMesh ) {
    std::unique_ptr<trimesh3_interface_impl> result( new trimesh3_interface_impl );
    result->set_mesh( boost::move( theMesh ) );
    return std::unique_ptr<trimesh3_interface>( result.release() );
}

class trimesh_vertex_channel : public mesh_channel {
    trimesh3_vertex_channel_general_accessor m_acc;

  public:
    trimesh_vertex_channel( const trimesh3_vertex_channel_general_accessor& acc, const frantic::tstring& name )
        : mesh_channel( name, acc.has_custom_faces() ? face_vertex : vertex, acc.data_type(), acc.arity(), acc.size(),
                        acc.face_count() )
        , m_acc( acc ) {}

    virtual void get_value( std::size_t index, void* outValue ) const {
        memcpy( outValue, const_cast<trimesh3_vertex_channel_general_accessor&>( m_acc ).data( index ),
                m_acc.primitive_size() );
    }

    virtual void set_value( std::size_t index, const void* value ) const {
        memcpy( const_cast<trimesh3_vertex_channel_general_accessor&>( m_acc ).data( index ), value,
                m_acc.primitive_size() );
    }

    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        return const_cast<trimesh3_vertex_channel_general_accessor&>( m_acc ).face( faceIndex )[fvertIndex];
    }

    virtual void set_fv_index( std::size_t faceIndex, std::size_t fvertIndex, std::size_t value ) const {
        const_cast<trimesh3_vertex_channel_general_accessor&>( m_acc ).face( faceIndex )[fvertIndex] =
            vector3::value_type( value );
    }

    virtual std::size_t get_num_face_verts( std::size_t /*faceIndex*/ ) const { return 3; }
};

class trimesh_face_channel : public mesh_channel {
    trimesh3_face_channel_general_accessor m_acc;

  public:
    trimesh_face_channel( const trimesh3_face_channel_general_accessor& acc, const frantic::tstring& name )
        : mesh_channel( name, face, acc.data_type(), acc.arity(), acc.size(), acc.size() )
        , m_acc( acc ) {}

    virtual void get_value( std::size_t index, void* outValue ) const {
        memcpy( outValue, const_cast<trimesh3_face_channel_general_accessor&>( m_acc ).data( index ),
                m_acc.primitive_size() );
    }

    virtual void set_value( std::size_t index, const void* value ) const {
        memcpy( const_cast<trimesh3_face_channel_general_accessor&>( m_acc ).data( index ), value,
                m_acc.primitive_size() );
    }

    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t /*fvertIndex*/ ) const { return faceIndex; }

    virtual std::size_t get_num_face_verts( std::size_t /*faceIndex*/ ) const { return 3; }
};

mesh_channel* trimesh3_interface_impl::create_vertex_channel( const frantic::tstring& channelName,
                                                              frantic::channels::data_type_t dataType, size_t arity ) {
    m_mesh->add_vertex_channel_raw( channelName, arity, dataType );
    return new trimesh_vertex_channel( m_mesh->get_vertex_channel_general_accessor( channelName ), channelName );
}

mesh_channel* trimesh3_interface_impl::create_vertex_channel_custom_faces( const frantic::tstring& channelName,
                                                                           frantic::channels::data_type_t dataType,
                                                                           size_t arity, size_t vertexCount ) {
    m_mesh->add_vertex_channel_raw( channelName, arity, dataType, vertexCount, true );
    return new trimesh_vertex_channel( m_mesh->get_vertex_channel_general_accessor( channelName ), channelName );
}

void trimesh3_interface_impl::destroy_vertex_channel( mesh_channel* channel ) {
    m_mesh->erase_vertex_channel( channel->get_name() );
    delete channel;
}

mesh_channel* trimesh3_interface_impl::create_face_channel( const frantic::tstring& channelName,
                                                            frantic::channels::data_type_t dataType, size_t arity ) {
    m_mesh->add_face_channel_raw( channelName, arity, dataType );
    return new trimesh_face_channel( m_mesh->get_face_channel_general_accessor( channelName ), channelName );
}

void trimesh3_interface_impl::destroy_face_channel( mesh_channel* channel ) {
    m_mesh->erase_face_channel( channel->get_name() );
    delete channel;
}

void trimesh3_interface_impl::set_mesh( BOOST_RV_REF( trimesh3 ) mesh ) {
    m_ownedMesh.swap( mesh );

    m_mesh = &m_ownedMesh;

    apply_underlying_channels();
}

void trimesh3_interface_impl::set_mesh_pointer( const trimesh3* mesh ) {
    m_ownedMesh.clear_and_deallocate();

    // TODO: why const_cast?  this seems dangerous
    m_mesh = const_cast<trimesh3*>( mesh );

    apply_underlying_channels();
}

void trimesh3_interface_impl::apply_underlying_channels() {
    std::vector<frantic::tstring> channels;

    m_mesh->get_vertex_channel_names( channels );
    for( std::vector<frantic::tstring>::const_iterator it = channels.begin(), itEnd = channels.end(); it != itEnd;
         ++it ) {
        std::unique_ptr<trimesh_vertex_channel> ch(
            new trimesh_vertex_channel( m_mesh->get_vertex_channel_general_accessor( *it ), *it ) );

        this->append_vertex_channel( std::move( ch ) );
    }

    channels.clear();

    m_mesh->get_face_channel_names( channels );
    for( std::vector<frantic::tstring>::const_iterator it = channels.begin(), itEnd = channels.end(); it != itEnd;
         ++it ) {
        std::unique_ptr<trimesh_face_channel> ch(
            new trimesh_face_channel( m_mesh->get_face_channel_general_accessor( *it ), *it ) );

        this->append_face_channel( std::move( ch ) );
    }
}

const trimesh3& trimesh3_interface_impl::get_trimesh() const { return *m_mesh; }

bool trimesh3_interface_impl::is_valid() const { return true; }

bool trimesh3_interface_impl::request_channel( const frantic::tstring& channelName, bool vertexChannel,
                                               bool /*forOutput*/, bool throwOnError ) {
    if( vertexChannel ) {
        if( !m_mesh->has_vertex_channel( channelName ) ) {
            if( throwOnError )
                throw std::runtime_error( "trimesh3_interface::request_channel() Failed to add vertex channel: \"" +
                                          frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
    } else {
        if( !m_mesh->has_face_channel( channelName ) ) {
            if( throwOnError )
                throw std::runtime_error( "trimesh3_interface::request_channel() Failed to add face channel: \"" +
                                          frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
    }
    return true;
}

std::size_t trimesh3_interface_impl::get_num_verts() const { return m_mesh->vertex_count(); }

void trimesh3_interface_impl::get_vert( std::size_t index, float ( &outValues )[3] ) const {
    const frantic::graphics::vector3f& v = m_mesh->get_vertex( index );
    outValues[0] = v.x;
    outValues[1] = v.y;
    outValues[2] = v.z;
}

bool trimesh3_interface_impl::is_read_only() const { return false; }

void trimesh3_interface_impl::set_vert( std::size_t index, const float* v ) {
    m_mesh->get_vertex( index ) = frantic::graphics::vector3f( v[0], v[1], v[2] );
}

std::size_t trimesh3_interface_impl::get_num_faces() const { return m_mesh->face_count(); }

std::size_t trimesh3_interface_impl::get_num_face_verts( std::size_t /*faceIndex*/ ) const { return 3; }

std::size_t trimesh3_interface_impl::get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
    return m_mesh->get_face( faceIndex )[fvertIndex];
}

void trimesh3_interface_impl::get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const {
    const frantic::graphics::vector3& f = m_mesh->get_face( faceIndex );
    outValues[0] = f.x;
    outValues[1] = f.y;
    outValues[2] = f.z;
}

void trimesh3_interface_impl::get_face_verts( std::size_t faceIndex, float outValues[][3] ) const {
    const frantic::graphics::vector3& f = m_mesh->get_face( faceIndex );
    const frantic::graphics::vector3f& v0 = m_mesh->get_vertex( f.x );
    const frantic::graphics::vector3f& v1 = m_mesh->get_vertex( f.y );
    const frantic::graphics::vector3f& v2 = m_mesh->get_vertex( f.z );
    outValues[0][0] = v0.x;
    outValues[0][1] = v0.y;
    outValues[0][2] = v0.z;
    outValues[1][0] = v1.x;
    outValues[1][1] = v1.y;
    outValues[1][2] = v1.z;
    outValues[2][0] = v2.x;
    outValues[2][1] = v2.y;
    outValues[2][2] = v2.z;
}

// Trimeshes don't support the concept of elements. TODO: Add that if enabled via some sort of option.
std::size_t trimesh3_interface_impl::get_num_elements() const { return 1; }

std::size_t trimesh3_interface_impl::get_face_element_index( std::size_t /*faceIndex*/ ) const { return 0; }

void trimesh3_interface_impl::init_adjacency() {
    if( !has_adjacency() ) {
        this->init_iterator_data();
    }
}

bool trimesh3_interface_impl::has_adjacency() const { return m_adjacencyDelegate.get() != NULL; }

bool trimesh3_interface_impl::init_vertex_iterator( vertex_iterator& vIt, std::size_t vertex ) const {
    assert( has_adjacency() && "init_iterator_data() must be called before using init_vertex_iterator()" );
    return m_adjacencyDelegate->init_vertex_iterator( vIt, vertex );
}

bool trimesh3_interface_impl::advance_vertex_iterator( vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->advance_vertex_iterator( vIt );
}

std::size_t trimesh3_interface_impl::get_edge_endpoint( vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_endpoint( vIt );
}

std::size_t trimesh3_interface_impl::get_edge_left_face( vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_left_face( vIt );
}

std::size_t trimesh3_interface_impl::get_edge_right_face( vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_right_face( vIt );
}

bool trimesh3_interface_impl::is_edge_visible( vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->is_edge_visible( vIt );
}

bool trimesh3_interface_impl::is_edge_boundary( vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->is_edge_boundary( vIt );
}

void trimesh3_interface_impl::init_face_iterator( face_iterator& fIt, std::size_t faceIndex ) const {
    assert( has_adjacency() );
    m_adjacencyDelegate->init_face_iterator( fIt, faceIndex );
}

bool trimesh3_interface_impl::advance_face_iterator( face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->advance_face_iterator( fIt );
}

std::size_t trimesh3_interface_impl::get_face_neighbor( face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_neighbor( fIt );
}

std::size_t trimesh3_interface_impl::get_face_next_vertex( face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_next_vertex( fIt );
}

std::size_t trimesh3_interface_impl::get_face_prev_vertex( face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_prev_vertex( fIt );
}

vertex_adjacency_interface& trimesh3_interface_impl::get_vertex_adjacency() {
    return *static_cast<vertex_adjacency_interface*>( NULL );
}

face_adjacency_interface& trimesh3_interface_impl::get_face_adjacency() {
    return *static_cast<face_adjacency_interface*>( NULL );
}

void trimesh3_interface_impl::init_iterator_data() {
    dcel adjacency;
    trimesh3_to_dcel( *m_mesh, adjacency );
    m_adjacencyDelegate.reset( new dcel_mesh_interface( boost::move( adjacency ) ) );
}

} // namespace geometry
} // namespace frantic
