// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/polymesh3_interface.hpp>

#include <boost/foreach.hpp>

#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/dcel_mesh_interface.hpp>
#include <frantic/geometry/polymesh3_accessors.hpp>

namespace frantic {
namespace geometry {

template <class Polymesh3Ptr>
struct polymesh3_ptr_traits;

template <>
struct polymesh3_ptr_traits<const_polymesh3_ptr> {
    typedef polymesh3_const_vertex_accessor<void> vertex_accessor_t;
    typedef polymesh3_const_face_accessor<void> face_accessor_t;

    static bool is_read_only() { return true; }

    static vertex_accessor_t get_vertex_accessor( const_polymesh3_ptr mesh, const frantic::tstring& channelName ) {
        return mesh->get_const_vertex_accessor( channelName );
    }

    static face_accessor_t get_face_accessor( const_polymesh3_ptr mesh, const frantic::tstring& channelName ) {
        return mesh->get_const_face_accessor( channelName );
    }
};

template <>
struct polymesh3_ptr_traits<polymesh3_ptr> {
    typedef polymesh3_vertex_accessor<void> vertex_accessor_t;
    typedef polymesh3_face_accessor<void> face_accessor_t;

    static bool is_read_only() { return false; }

    static vertex_accessor_t get_vertex_accessor( polymesh3_ptr mesh, const frantic::tstring& channelName ) {
        return mesh->get_vertex_accessor( channelName );
    }

    static face_accessor_t get_face_accessor( polymesh3_ptr mesh, const frantic::tstring& channelName ) {
        return mesh->get_face_accessor( channelName );
    }
};

template <class Polymesh3Ptr>
class polymesh3_interface_impl : public polymesh3_interface {
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
    polymesh3_interface_impl( Polymesh3Ptr mesh );

    void set_mesh( Polymesh3Ptr mesh );

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

  private:
    polymesh3_interface_impl(); // not implemented

    Polymesh3Ptr m_mesh;
    frantic::geometry::polymesh3_const_vertex_accessor<frantic::graphics::vector3f> m_geomAcc;

    std::unique_ptr<dcel_mesh_interface> m_adjacencyDelegate;
};

std::unique_ptr<polymesh3_interface> polymesh3_interface::create_instance( frantic::geometry::polymesh3_ptr mesh ) {
    std::unique_ptr<polymesh3_interface> result(
        new polymesh3_interface_impl<frantic::geometry::polymesh3_ptr>( mesh ) );
    return result;
}

std::unique_ptr<polymesh3_interface>
polymesh3_interface::create_const_instance( frantic::geometry::const_polymesh3_ptr mesh ) {
    std::unique_ptr<polymesh3_interface> result(
        new polymesh3_interface_impl<frantic::geometry::const_polymesh3_ptr>( mesh ) );
    return result;
}

namespace {

frantic::geometry::mesh_channel::channel_type get_object_type( frantic::geometry::const_polymesh3_ptr mesh,
                                                               const frantic::tstring& channelName ) {
    const frantic::geometry::polymesh3_channel& channel = mesh->get_channel_info( channelName );
    if( channel.is_vertex_channel() ) {
        if( channel.get_faces().empty() ) {
            return frantic::geometry::mesh_channel::vertex;
        } else {
            return frantic::geometry::mesh_channel::face_vertex;
        }
    } else {
        return frantic::geometry::mesh_channel::face;
    }
}

frantic::channels::data_type_t get_type( frantic::geometry::const_polymesh3_ptr mesh,
                                         const frantic::tstring& channelName ) {
    const frantic::geometry::polymesh3_channel& channel = mesh->get_channel_info( channelName );
    return channel.type();
}

std::size_t get_arity( frantic::geometry::const_polymesh3_ptr mesh, const frantic::tstring& channelName ) {
    const frantic::geometry::polymesh3_channel& channel = mesh->get_channel_info( channelName );
    return channel.arity();
}

std::size_t get_num_data( frantic::geometry::const_polymesh3_ptr mesh, const frantic::tstring& channelName ) {
    const frantic::geometry::polymesh3_channel& channel = mesh->get_channel_info( channelName );
    if( channel.is_vertex_channel() ) {
        polymesh3_const_vertex_accessor<void> acc = mesh->get_const_vertex_accessor( channelName );
        return acc.vertex_count();
    } else {
        polymesh3_const_face_accessor<void> acc = mesh->get_const_face_accessor( channelName );
        return acc.face_count();
    }
}

} // anonymous namespace

template <class Polymesh3Ptr>
class polymesh_vertex_channel : public mesh_channel {
  public:
    polymesh_vertex_channel( Polymesh3Ptr mesh, const frantic::tstring& name )
        : mesh_channel( name, get_object_type( mesh, name ), get_type( mesh, name ), get_arity( mesh, name ),
                        get_num_data( mesh, name ), mesh->face_count(),
                        polymesh3_ptr_traits<Polymesh3Ptr>::is_read_only() )
        , m_acc( polymesh3_ptr_traits<Polymesh3Ptr>::get_vertex_accessor( mesh, name ) )
        , m_geomAcc( mesh->get_const_vertex_accessor( _T("verts") ) ) {}

    virtual void get_value( std::size_t index, void* outValue ) const {
        memcpy( outValue, m_acc.get_vertex( index ), get_element_size() );
    }

    virtual void set_value( std::size_t index, const void* value ) const;

    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        if( m_acc.face_count() ) {
            return m_acc.get_face( faceIndex ).first[fvertIndex];
        } else {
            return m_geomAcc.get_face( faceIndex ).first[fvertIndex];
        }
    }

    virtual void set_fv_index( std::size_t faceIndex, std::size_t fvertIndex, std::size_t value ) const;

    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const {
        return m_geomAcc.get_face_degree( faceIndex );
    }

  private:
    mutable typename polymesh3_ptr_traits<Polymesh3Ptr>::vertex_accessor_t
        m_acc; // mutable because set_value() is const in interface
    polymesh3_const_vertex_accessor<void> m_geomAcc;
};

template <>
void polymesh_vertex_channel<const_polymesh3_ptr>::set_value( std::size_t /*index*/, const void* /*value*/ ) const {
    // From docs for is_writeable(): "If false, any calls to this->set_value() are discarded."
}

template <>
void polymesh_vertex_channel<polymesh3_ptr>::set_value( std::size_t index, const void* value ) const {
    memcpy( static_cast<void*>( m_acc.get_vertex( index ) ), value, get_element_size() );
}

template <>
void polymesh_vertex_channel<const_polymesh3_ptr>::set_fv_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/,
                                                                 std::size_t /*value*/ ) const {
    throw std::runtime_error( "polymesh_vertex_channel::set_fv_index() Cannot set index in read-only mesh" );
}

template <>
void polymesh_vertex_channel<polymesh3_ptr>::set_fv_index( std::size_t faceIndex, std::size_t fvertIndex,
                                                           std::size_t value ) const {
    if( m_acc.face_count() ) {
        m_acc.get_face( faceIndex ).first[fvertIndex] = (int)value;
    } else {
        throw std::runtime_error(
            "polymesh_vertex_channel::set_fv_index -- This channel does not have custom face indices" );
    }
}

template <class Polymesh3Ptr>
class polymesh_face_channel : public mesh_channel {
  public:
    polymesh_face_channel( Polymesh3Ptr mesh, const frantic::tstring& name )
        : mesh_channel( name, get_object_type( mesh, name ), get_type( mesh, name ), get_arity( mesh, name ),
                        get_num_data( mesh, name ), mesh->face_count(),
                        polymesh3_ptr_traits<Polymesh3Ptr>::is_read_only() )
        , m_acc( polymesh3_ptr_traits<Polymesh3Ptr>::get_face_accessor( mesh, name ) )
        , m_geomAcc( mesh->get_const_vertex_accessor( _T("verts") ) ) {}

    virtual void get_value( std::size_t index, void* outValue ) const {
        memcpy( outValue, m_acc.get_face( index ), get_element_size() );
    }

    virtual void set_value( std::size_t index, const void* value ) const {
        memcpy( static_cast<void*>( m_acc.get_face( index ) ), value, get_element_size() );
    }

    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        return m_geomAcc.get_face( faceIndex ).first[fvertIndex];
    }

    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const {
        return m_geomAcc.get_face_degree( faceIndex );
    }

  private:
    mutable typename polymesh3_ptr_traits<Polymesh3Ptr>::face_accessor_t
        m_acc; // mutable because set_value() is const in interface
    polymesh3_const_vertex_accessor<void> m_geomAcc;
};

template <>
void polymesh_face_channel<const_polymesh3_ptr>::set_value( std::size_t /*index*/, const void* /*value*/ ) const {
    // From docs for is_writeable(): "If false, any calls to this->set_value() are discarded."
}

template <>
void polymesh_face_channel<polymesh3_ptr>::set_value( std::size_t index, const void* value ) const {
    memcpy( static_cast<void*>( m_acc.get_face( index ) ), value, get_element_size() );
}

template <>
mesh_channel* polymesh3_interface_impl<const_polymesh3_ptr>::create_vertex_channel(
    const frantic::tstring& channelName, frantic::channels::data_type_t /*dataType*/, size_t /*arity*/ ) {
    throw std::runtime_error( "polymesh3_interface::create_vertex_channel() Cannot add vertex channel \"" +
                              frantic::strings::to_string( channelName ) + "\" to read-only mesh" );
}

template <>
mesh_channel* polymesh3_interface_impl<polymesh3_ptr>::create_vertex_channel( const frantic::tstring& channelName,
                                                                              frantic::channels::data_type_t dataType,
                                                                              size_t arity ) {
    m_mesh->add_empty_vertex_channel( channelName, dataType, arity );
    return new polymesh_vertex_channel<polymesh3_ptr>( m_mesh, channelName );
}

template <>
mesh_channel* polymesh3_interface_impl<const_polymesh3_ptr>::create_vertex_channel_custom_faces(
    const frantic::tstring& channelName, frantic::channels::data_type_t /*dataType*/, size_t /*arity*/,
    size_t /*vertexCount*/ ) {
    throw std::runtime_error( "polymesh3_interface::create_vertex_channel_custom_faces() Cannot add vertex channel \"" +
                              frantic::strings::to_string( channelName ) + "\" to read-only mesh" );
}

template <>
mesh_channel* polymesh3_interface_impl<polymesh3_ptr>::create_vertex_channel_custom_faces(
    const frantic::tstring& channelName, frantic::channels::data_type_t dataType, size_t arity, size_t vertexCount ) {
    m_mesh->add_empty_vertex_channel( channelName, dataType, arity, vertexCount );
    return new polymesh_vertex_channel<polymesh3_ptr>( m_mesh, channelName );
}

template <>
void polymesh3_interface_impl<const_polymesh3_ptr>::destroy_vertex_channel( mesh_channel* channel ) {
    throw std::runtime_error( "polymesh3_interface::destroy_vertex_channel() Cannot destroy vertex channel \"" +
                              frantic::strings::to_string( channel->get_name() ) + "\" in read-only mesh" );
}

template <>
void polymesh3_interface_impl<polymesh3_ptr>::destroy_vertex_channel( mesh_channel* channel ) {
    m_mesh->erase_vertex_channel( channel->get_name() );
    delete channel;
}

template <>
mesh_channel* polymesh3_interface_impl<const_polymesh3_ptr>::create_face_channel(
    const frantic::tstring& channelName, frantic::channels::data_type_t /*dataType*/, size_t /*arity*/ ) {
    throw std::runtime_error( "polymesh3_interface::create_face_channel() Cannot add face channel \"" +
                              frantic::strings::to_string( channelName ) + "\" to read-only mesh" );
}

template <>
mesh_channel* polymesh3_interface_impl<polymesh3_ptr>::create_face_channel( const frantic::tstring& channelName,
                                                                            frantic::channels::data_type_t dataType,
                                                                            size_t arity ) {
    m_mesh->add_empty_face_channel( channelName, dataType, arity );
    return new polymesh_face_channel<polymesh3_ptr>( m_mesh, channelName );
}

template <>
void polymesh3_interface_impl<const_polymesh3_ptr>::destroy_face_channel( mesh_channel* channel ) {
    throw std::runtime_error( "polymesh3_interface::destroy_face_channel() Cannot destroy face channel \"" +
                              frantic::strings::to_string( channel->get_name() ) + "\" in read-only mesh" );
}

template <>
void polymesh3_interface_impl<polymesh3_ptr>::destroy_face_channel( mesh_channel* channel ) {
    m_mesh->erase_face_channel( channel->get_name() );
    delete channel;
}

template <class Polymesh3Ptr>
polymesh3_interface_impl<Polymesh3Ptr>::polymesh3_interface_impl( Polymesh3Ptr mesh ) {
    set_mesh( mesh );
}

template <class Polymesh3Ptr>
void polymesh3_interface_impl<Polymesh3Ptr>::set_mesh( Polymesh3Ptr mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "polymesh3_interface::set_mesh Error: mesh is NULL" );
    }

    m_geomAcc = mesh->template get_const_vertex_accessor<frantic::graphics::vector3f>( _T("verts") );
    m_mesh = mesh;

    std::vector<frantic::tstring> channelNames;

    m_mesh->get_vertex_channel_names( channelNames );
    BOOST_FOREACH( const frantic::tstring& channelName, channelNames ) {
        std::unique_ptr<polymesh_vertex_channel<Polymesh3Ptr>> ch(
            new polymesh_vertex_channel<Polymesh3Ptr>( m_mesh, channelName ) );

        this->append_vertex_channel( std::move( ch ) );
    }

    channelNames.clear();

    m_mesh->get_face_channel_names( channelNames );
    BOOST_FOREACH( const frantic::tstring& channelName, channelNames ) {
        std::unique_ptr<polymesh_face_channel<Polymesh3Ptr>> ch(
            new polymesh_face_channel<Polymesh3Ptr>( m_mesh, channelName ) );

        this->append_face_channel( std::move( ch ) );
    }
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::is_valid() const {
    return m_mesh != 0;
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::request_channel( const frantic::tstring& channelName, bool vertexChannel,
                                                              bool forOutput, bool throwOnError ) {
    if( !m_mesh ) {
        throw std::runtime_error( "polymesh3_interface::request_channel() mesh is NULL" );
    }

    if( vertexChannel ) {
        if( !m_mesh->has_vertex_channel( channelName ) ) {
            if( throwOnError )
                throw std::runtime_error( "polymesh3_interface::request_channel() Failed to add vertex channel: \"" +
                                          frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
        if( forOutput && polymesh3_ptr_traits<Polymesh3Ptr>::is_read_only() ) {
            if( throwOnError )
                throw std::runtime_error( "polymesh3_interface::request_channel() Cannot get writable vertex channel "
                                          "from read-only mesh: \"" +
                                          frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
    } else {
        if( !m_mesh->has_face_channel( channelName ) ) {
            if( throwOnError )
                throw std::runtime_error( "polymesh3_interface::request_channel() Failed to add face channel: \"" +
                                          frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
        if( forOutput && polymesh3_ptr_traits<Polymesh3Ptr>::is_read_only() ) {
            if( throwOnError )
                throw std::runtime_error(
                    "polymesh3_interface::request_channel() Cannot get writable face channel from read-only mesh: \"" +
                    frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
    }
    return true;
}

template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_num_verts() const {
    if( !m_mesh ) {
        throw std::runtime_error( "polymesh3_interface::get_num_verts() mesh is NULL" );
    }

    return m_mesh->vertex_count();
}

template <class Polymesh3Ptr>
void polymesh3_interface_impl<Polymesh3Ptr>::get_vert( std::size_t index, float ( &outValues )[3] ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "polymesh3_interface::get_vert() mesh is NULL" );
    }

    const frantic::graphics::vector3f vert = m_mesh->get_vertex( index );

    for( int axis = 0; axis < 3; ++axis ) {
        outValues[axis] = vert[axis];
    }
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::is_read_only() const {
    return polymesh3_ptr_traits<Polymesh3Ptr>::is_read_only();
}

template <>
void polymesh3_interface_impl<const_polymesh3_ptr>::set_vert( std::size_t /*index*/, const float* /*v*/ ) {
    throw std::runtime_error( "polymesh3_interface::set_vert() cannot operate on read-only mesh" );
}

template <>
void polymesh3_interface_impl<polymesh3_ptr>::set_vert( std::size_t index, const float* v ) {
    m_mesh->get_cvt_vertex_accessor<frantic::graphics::vector3f>( _T( "verts" ) )
        .set_vertex( index, frantic::graphics::vector3f( v[0], v[1], v[2] ) );
}

template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_num_faces() const {
    if( !m_mesh ) {
        throw std::runtime_error( "polymesh3_interface::get_num_faces() mesh is NULL" );
    }

    return m_mesh->face_count();
}

template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_num_face_verts( std::size_t faceIndex ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "polymesh3_interface::get_num_face_verts() mesh is NULL" );
    }

    return m_geomAcc.get_face_degree( faceIndex );
}

template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_face_vert_index( std::size_t faceIndex,
                                                                         std::size_t fvertIndex ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "polymesh3_interface::get_face_vert_index() mesh is NULL" );
    }

    return m_geomAcc.get_face( faceIndex ).first[fvertIndex];
}

template <class Polymesh3Ptr>
void polymesh3_interface_impl<Polymesh3Ptr>::get_face_vert_indices( std::size_t faceIndex,
                                                                    std::size_t outValues[] ) const {
    frantic::geometry::polymesh3_const_face_range face = m_geomAcc.get_face( faceIndex );
    for( frantic::geometry::polymesh3_const_face_iterator i = face.first, ie = face.second; i != ie; ++i ) {
        outValues[0] = *i;
        ++outValues;
    }
}

template <class Polymesh3Ptr>
void polymesh3_interface_impl<Polymesh3Ptr>::get_face_verts( std::size_t faceIndex, float outValues[][3] ) const {
    frantic::geometry::polymesh3_const_face_range face = m_geomAcc.get_face( faceIndex );
    for( frantic::geometry::polymesh3_const_face_iterator i = face.first, ie = face.second; i != ie; ++i ) {
        const frantic::graphics::vector3f vert = m_geomAcc.get_vertex( *i );
        for( int axis = 0; axis < 3; ++axis ) {
            outValues[0][axis] = vert[axis];
        }
        ++outValues;
    }
}

// Polymeshes don't support the concept of elements. TODO: Add that if enabled via some sort of option.
template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_num_elements() const {
    return 1;
}

template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_face_element_index( std::size_t /*faceIndex*/ ) const {
    return 0;
}

template <class Polymesh3Ptr>
void polymesh3_interface_impl<Polymesh3Ptr>::init_adjacency() {
    if( !has_adjacency() ) {
        dcel adjacency;
        polymesh3_to_dcel( m_mesh, adjacency );
        m_adjacencyDelegate.reset( new dcel_mesh_interface( boost::move( adjacency ) ) );
    }
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::has_adjacency() const {
    return m_adjacencyDelegate.get() != NULL;
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::init_vertex_iterator( frantic::geometry::vertex_iterator& vIt,
                                                                   std::size_t vertexIndex ) const {
    assert( has_adjacency() && "init_adjacency() must be called before using init_vertex_iterator()" );
    return m_adjacencyDelegate->init_vertex_iterator( vIt, vertexIndex );
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->advance_vertex_iterator( vIt );
}

template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_endpoint( vIt );
}

template <class Polymesh3Ptr>
std::size_t
polymesh3_interface_impl<Polymesh3Ptr>::get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_left_face( vIt );
}

template <class Polymesh3Ptr>
std::size_t
polymesh3_interface_impl<Polymesh3Ptr>::get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_right_face( vIt );
}

// As far as I know, no information regarding edge visibility is stored with polymesh3
template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->is_edge_visible( vIt );
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->is_edge_boundary( vIt );
}

template <class Polymesh3Ptr>
void polymesh3_interface_impl<Polymesh3Ptr>::init_face_iterator( frantic::geometry::face_iterator& fIt,
                                                                 std::size_t faceIndex ) const {
    assert( has_adjacency() && "init_adjacency() must be called before using init_face_iterator()" );
    return m_adjacencyDelegate->init_face_iterator( fIt, faceIndex );
}

template <class Polymesh3Ptr>
bool polymesh3_interface_impl<Polymesh3Ptr>::advance_face_iterator( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->advance_face_iterator( fIt );
}

template <class Polymesh3Ptr>
std::size_t polymesh3_interface_impl<Polymesh3Ptr>::get_face_neighbor( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_neighbor( fIt );
}

template <class Polymesh3Ptr>
std::size_t
polymesh3_interface_impl<Polymesh3Ptr>::get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_prev_vertex( fIt );
}

template <class Polymesh3Ptr>
std::size_t
polymesh3_interface_impl<Polymesh3Ptr>::get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_next_vertex( fIt );
}

} // namespace geometry
} // namespace frantic
