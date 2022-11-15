// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/transformed_mesh_interface.hpp>
#include <frantic/math/eigen.hpp>

namespace {
// we have to create our own face iterator if m_reverseFaceWinding
struct face_half_edge {
    std::size_t next;
    std::size_t face;
    std::size_t twinFace;
};
struct face_iterator_impl {
    const face_half_edge *cur, *end;
};
BOOST_STATIC_ASSERT( sizeof( face_iterator_impl ) <= frantic::geometry::detail::iterator_storage_size );
} // anonymous namespace

namespace frantic {
namespace geometry {

class transformed_mesh_interface_impl : public transformed_mesh_interface {
  private:
    frantic::geometry::mesh_interface_ptr m_mesh;
    const frantic::graphics::transform4f m_xform;
    const bool m_reverseFaceWinding;
    std::vector<face_half_edge> m_faceHalfEdgeData;
    std::vector<size_t> m_faceHalfEdges;

    transformed_mesh_interface_impl( const transformed_mesh_interface_impl& );            // not implemented
    transformed_mesh_interface_impl& operator=( const transformed_mesh_interface_impl& ); // not implemented

  protected:
    // does not support creating new or deleting channels because it does not know if the mesh_interface it is
    // transforming supports them
    //    and we can't access a mesh_interface's protected functions.

  public:
    transformed_mesh_interface_impl( frantic::geometry::mesh_interface_ptr mesh,
                                     const frantic::graphics::transform4f& xform );

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
};

std::unique_ptr<transformed_mesh_interface>
transformed_mesh_interface::create_instance( frantic::geometry::mesh_interface_ptr mesh,
                                             const frantic::graphics::transform4f& xform ) {
    std::unique_ptr<transformed_mesh_interface> result( new transformed_mesh_interface_impl( mesh, xform ) );
    return result;
}

/*
 * a transformed mesh channel
 * this mesh_channel is read only
 * \note If Normal or Velocity channel and transform_type::none class will
 *    default them to transform_type::normal or transform_type::vector respectively
 */
class transformed_mesh_channel : public mesh_channel {
    const mesh_channel* m_channel;
    const frantic::graphics::transform4f m_xform;
    const frantic::graphics::transform4f m_xformInverse;
    const bool m_reverseFaceWinding;
    float m_transformScaleMagnitude;

    transformed_mesh_channel( const transformed_mesh_channel& );            // not implemented
    transformed_mesh_channel& operator=( const transformed_mesh_channel& ); // not implemented

    void transform_channel_value( void* value ) const {
        if( this->get_transform_type() != transform_type::none && !m_xform.is_identity() ) {
            frantic::channels::channel_type_convertor_function_t convertFromChannel =
                get_channel_type_convertor_function( this->get_data_type(), frantic::channels::data_type_float32,
                                                     this->get_name() );
            frantic::channels::channel_type_convertor_function_t convertToChannel = get_channel_type_convertor_function(
                frantic::channels::data_type_float32, this->get_data_type(), this->get_name() );

            if( this->get_transform_type() == transform_type::point ||
                this->get_transform_type() == transform_type::vector ||
                this->get_transform_type() == transform_type::normal ) {
                frantic::graphics::vector3f convertedValue, transformedValue;
                convertFromChannel( reinterpret_cast<char*>( &convertedValue ), reinterpret_cast<char*>( value ), 3 );

                if( this->get_transform_type() == transform_type::point ) {
                    transformedValue = m_xform * convertedValue;

                } else if( this->get_transform_type() == transform_type::vector ) {
                    transformedValue = m_xform.transform_no_translation( convertedValue );

                } else if( this->get_transform_type() == transform_type::normal ) {
                    transformedValue = m_xformInverse.transpose_transform_no_translation( convertedValue );
                }

                convertToChannel( reinterpret_cast<char*>( value ), reinterpret_cast<char*>( &transformedValue ), 3 );

            } else if( this->get_transform_type() == transform_type::scale ) {
                float convertedValue, transformedValue;
                convertFromChannel( reinterpret_cast<char*>( &convertedValue ), reinterpret_cast<char*>( value ), 1 );

                transformedValue = m_transformScaleMagnitude * convertedValue;

                convertToChannel( reinterpret_cast<char*>( value ), reinterpret_cast<char*>( &transformedValue ), 1 );
            }
        }
    }

  public:
    transformed_mesh_channel( const mesh_channel* ch, const frantic::graphics::transform4f& xform )
        : mesh_channel( ch->get_name(), ch->get_channel_type(), ch->get_data_type(), ch->get_data_arity(),
                        ch->get_num_elements(), ch->get_num_faces(), true )
        , m_channel( ch )
        , m_xform( xform )
        , m_xformInverse( xform.to_inverse() )
        , m_reverseFaceWinding( xform.is_orientation_inverting() )
        , m_transformScaleMagnitude( 1 ) {

        if( m_channel->get_transform_type() != transform_type::none ) {
            // copy the channels transfrom
            this->set_transform_type( m_channel->get_transform_type() );
        } else {
            // we want these channels to be transformed by default
            if( this->get_channel_type() == mesh_channel::vertex ) {
                if( this->get_name() == _T("Normal") ) {
                    this->set_transform_type( transform_type::normal );
                } else if( this->get_name() == _T("Velocity") ) {
                    this->set_transform_type( transform_type::vector );
                }
            }
        }
        // this m_transformScaleMagnitude code copied from transformed_particle_istream
        frantic::graphics::vector3f translation;
        frantic::graphics::transform4f perspective, rotation, stretch;
        m_xform.decompose( perspective, translation, rotation, stretch );

        // Extract the scale from the stretch matrix (its eigenvalues).
        // Use maximum scale component as m_transformScaleMagnitude,
        // which is used to scale the particle radius.
        float stretchUpperTriangle[6] = { stretch.get( 0, 0 ), stretch.get( 1, 0 ), stretch.get( 2, 0 ),
                                          stretch.get( 1, 1 ), stretch.get( 2, 1 ), stretch.get( 2, 2 ) };
        frantic::graphics::vector3f eigenvalues;
        frantic::math::linearalgebra::get_eigenvalues_symmetric_3x3( stretchUpperTriangle, eigenvalues[0],
                                                                     eigenvalues[1], eigenvalues[2] );
        m_transformScaleMagnitude = eigenvalues.max_abs_component();
    }

    virtual void get_value( std::size_t index, void* outValue ) const {
        m_channel->get_value( index, outValue );
        transform_channel_value( outValue );
    }
    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        if( m_reverseFaceWinding ) {
            fvertIndex = this->get_num_face_verts( faceIndex ) - 1 - fvertIndex;
        }
        return m_channel->get_fv_index( faceIndex, fvertIndex );
    }

    virtual void set_value( std::size_t /*index*/, const void* /*value*/ ) const {
        // not writable
    }
    virtual void set_fv_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/, std::size_t /*value*/ ) const {
        // not writable
    }
    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const {
        return m_channel->get_num_face_verts( faceIndex );
    }
};

transformed_mesh_interface_impl::transformed_mesh_interface_impl( frantic::geometry::mesh_interface_ptr mesh,
                                                                  const frantic::graphics::transform4f& xform )
    : m_mesh( mesh )
    , m_xform( xform )
    , m_reverseFaceWinding( xform.is_orientation_inverting() ) {
    if( !m_mesh ) {
        throw std::runtime_error( "transformed_mesh_interface::set_mesh Error: mesh is NULL" );
    }
    // init the channel maps
    mesh_interface::mesh_channel_map::const_iterator it = m_mesh->get_vertex_channels().begin();
    mesh_interface::mesh_channel_map::const_iterator itEnd = m_mesh->get_vertex_channels().end();
    for( ; it != itEnd; ++it ) {
        this->append_vertex_channel(
            std::unique_ptr<transformed_mesh_channel>( new transformed_mesh_channel( it->second, m_xform ) ) );
    }

    it = m_mesh->get_face_channels().begin();
    itEnd = m_mesh->get_face_channels().end();
    for( ; it != itEnd; ++it ) {
        this->append_face_channel(
            std::unique_ptr<transformed_mesh_channel>( new transformed_mesh_channel( it->second, m_xform ) ) );
    }
}

bool transformed_mesh_interface_impl::is_valid() const { return m_mesh != 0; }
bool transformed_mesh_interface_impl::request_channel( const frantic::tstring& channelName, bool vertexChannel,
                                                       bool /*forOutput*/, bool throwOnError ) {
    if( vertexChannel ) {
        if( !m_mesh->has_vertex_channel( channelName ) ) {
            if( throwOnError )
                throw std::runtime_error(
                    "transformed_mesh_interface::request_channel() Failed to add vertex channel: \"" +
                    frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
    } else {
        if( !m_mesh->has_face_channel( channelName ) ) {
            if( throwOnError )
                throw std::runtime_error(
                    "transformed_mesh_interface::request_channel() Failed to add face channel: \"" +
                    frantic::strings::to_string( channelName ) + "\"" );
            return false;
        }
    }
    return true;
}
void transformed_mesh_interface_impl::get_vert( std::size_t index, float ( &outValues )[3] ) const {
    m_mesh->get_vert( index, outValues );
    frantic::graphics::vector3f* vert = reinterpret_cast<frantic::graphics::vector3f*>( &outValues[0] );
    ( *vert ) = m_xform * ( *vert );
}

std::size_t transformed_mesh_interface_impl::get_face_vert_index( std::size_t faceIndex,
                                                                  std::size_t fvertIndex ) const {
    if( m_reverseFaceWinding ) {
        fvertIndex = m_mesh->get_num_face_verts( faceIndex ) - 1 - fvertIndex;
    }
    return m_mesh->get_face_vert_index( faceIndex, fvertIndex );
}
void transformed_mesh_interface_impl::get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const {
    m_mesh->get_face_vert_indices( faceIndex, outValues );
    if( m_reverseFaceWinding ) {
        for( std::size_t i = 0; i < get_num_face_verts( faceIndex ) / 2; ++i ) {
            std::swap( outValues[i], outValues[get_num_face_verts( faceIndex ) - 1 - i] );
        }
    }
}
void transformed_mesh_interface_impl::get_face_verts( std::size_t faceIndex, float outValues[][3] ) const {
    m_mesh->get_face_verts( faceIndex, outValues );
    if( m_reverseFaceWinding ) {
        for( std::size_t i = 0; i < get_num_face_verts( faceIndex ) / 2; ++i ) {
            for( std::size_t j = 0; j < 3; ++j )
                std::swap( outValues[i][j], outValues[get_num_face_verts( faceIndex ) - 1 - i][j] );
        }
    }
}

std::size_t transformed_mesh_interface_impl::get_num_verts() const { return m_mesh->get_num_verts(); }
std::size_t transformed_mesh_interface_impl::get_num_faces() const { return m_mesh->get_num_faces(); }
std::size_t transformed_mesh_interface_impl::get_num_face_verts( std::size_t faceIndex ) const {
    return m_mesh->get_num_face_verts( faceIndex );
}
std::size_t transformed_mesh_interface_impl::get_num_elements() const { return m_mesh->get_num_elements(); }

std::size_t transformed_mesh_interface_impl::get_face_element_index( std::size_t faceIndex ) const {
    return m_mesh->get_face_element_index( faceIndex );
}

void transformed_mesh_interface_impl::init_adjacency() {
    m_mesh->init_adjacency();
    // we have to create our own face iterator if m_reverseFaceWinding
    if( m_reverseFaceWinding ) {
        if( m_faceHalfEdgeData.empty() ) {
            m_faceHalfEdgeData.reserve( m_mesh->get_num_faces() * 3 ); // we have at least this many half edges

            for( std::size_t i = 0; i < m_mesh->get_num_faces(); ++i ) {
                frantic::geometry::face_iterator fIt;
                m_mesh->init_face_iterator( fIt, i );

                face_half_edge faceHalfEdge;
                faceHalfEdge.face = i;
                // point to the last half edge in this face
                faceHalfEdge.next = m_faceHalfEdgeData.size() - 1 + m_mesh->get_num_face_verts( i );
                faceHalfEdge.twinFace = m_mesh->get_face_neighbor( fIt );
                m_faceHalfEdgeData.push_back( faceHalfEdge );

                while( m_mesh->advance_face_iterator( fIt ) ) {
                    face_half_edge faceHalfEdge;
                    faceHalfEdge.face = i;
                    // point to the previous half edge in this face
                    faceHalfEdge.next = m_faceHalfEdgeData.size() - 1;
                    faceHalfEdge.twinFace = m_mesh->get_face_neighbor( fIt );
                    m_faceHalfEdgeData.push_back( faceHalfEdge );
                }

                m_faceHalfEdges.push_back( m_faceHalfEdgeData.size() - 1 );
            }
        }
    }
}

bool transformed_mesh_interface_impl::has_adjacency() const { return !m_faceHalfEdgeData.empty(); }

std::size_t transformed_mesh_interface_impl::get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const {
    if( !m_reverseFaceWinding ) {
        return m_mesh->get_edge_left_face( vIt );
    } else {
        return m_mesh->get_edge_right_face( vIt );
    }
}
std::size_t transformed_mesh_interface_impl::get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const {
    if( !m_reverseFaceWinding ) {
        return m_mesh->get_edge_right_face( vIt );
    } else {
        return m_mesh->get_edge_left_face( vIt );
    }
}

void transformed_mesh_interface_impl::init_face_iterator( frantic::geometry::face_iterator& fIt,
                                                          std::size_t faceIndex ) const {
    if( !m_reverseFaceWinding ) {
        m_mesh->init_face_iterator( fIt, faceIndex );
    } else {
        assert( ( !m_faceHalfEdgeData.empty() || m_mesh->get_num_faces() == 0 ) &&
                "init_adjacency() must be called before using init_face_iterator()" );
        face_iterator_impl& it = *static_cast<face_iterator_impl*>( fIt.m_data.address() );
        it.cur = it.end = &m_faceHalfEdgeData[m_faceHalfEdges[faceIndex]];
    }
}
bool transformed_mesh_interface_impl::advance_face_iterator( frantic::geometry::face_iterator& fIt ) const {
    if( !m_reverseFaceWinding ) {
        return m_mesh->advance_face_iterator( fIt );
    } else {
        face_iterator_impl& it = *static_cast<face_iterator_impl*>( fIt.m_data.address() );
        it.cur = &m_faceHalfEdgeData[it.cur->next];
        return it.cur != it.end;
    }
}
std::size_t transformed_mesh_interface_impl::get_face_neighbor( frantic::geometry::face_iterator& fIt ) const {
    if( !m_reverseFaceWinding ) {
        return m_mesh->get_face_neighbor( fIt );
    } else {
        return static_cast<face_iterator_impl*>( fIt.m_data.address() )->cur->twinFace;
    }
}

bool transformed_mesh_interface_impl::init_vertex_iterator( frantic::geometry::vertex_iterator& vIt,
                                                            std::size_t vertexIndex ) const {
    return m_mesh->init_vertex_iterator( vIt, vertexIndex );
}
bool transformed_mesh_interface_impl::advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const {
    return m_mesh->advance_vertex_iterator( vIt );
}
std::size_t transformed_mesh_interface_impl::get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const {
    return m_mesh->get_edge_endpoint( vIt );
}
bool transformed_mesh_interface_impl::is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const {
    return m_mesh->is_edge_visible( vIt );
}
bool transformed_mesh_interface_impl::is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const {
    return m_mesh->is_edge_boundary( vIt );
}

std::size_t transformed_mesh_interface_impl::get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const {
    if( !m_reverseFaceWinding ) {
        return m_mesh->get_face_prev_vertex( fIt );
    } else {
        return m_mesh->get_face_next_vertex( fIt );
    }
}

std::size_t transformed_mesh_interface_impl::get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const {
    if( !m_reverseFaceWinding ) {
        return m_mesh->get_face_next_vertex( fIt );
    } else {
        return m_mesh->get_face_prev_vertex( fIt );
    }
}
} // namespace geometry
} // namespace frantic
