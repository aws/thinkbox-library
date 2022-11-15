// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/fractional_mesh_interface.hpp>
#include <frantic/math/fractions.hpp>
#include <frantic/math/uint128.hpp>

namespace frantic {
namespace geometry {

/**
 * This class provides access to a list of evenly distributed indicies between 0 and totalCount based upon a given
 * fraction. If limit is less than totalCount*fraction, the list will be limited to that amount.
 */
class fractional_index_map {
  private:
    std::vector<boost::int64_t> m_inputs;
    boost::int64_t m_totalCount;

  public:
    fractional_index_map( boost::int64_t totalCount, double fraction,
                          boost::int64_t limit = std::numeric_limits<boost::int64_t>::max() ) {
        boost::int64_t lastOutput = std::numeric_limits<boost::int64_t>::max();
        boost::int64_t lastInput = std::numeric_limits<boost::int64_t>::max();

        fraction = frantic::math::clamp( fraction, 0.0, 1.0 );

        std::pair<boost::int64_t, boost::int64_t> rational = frantic::math::get_rational_representation( fraction );
        boost::int64_t numerator = rational.first;
        boost::int64_t denominator = rational.second;
        boost::int64_t accumulator = 0;

        using frantic::math::uint128;

        m_totalCount = std::min(
            limit, ( uint128( totalCount ) * uint128( numerator ) / denominator ).to_integral<boost::int64_t>() );
        m_inputs.resize( m_totalCount );

        for( boost::int64_t output = 0; output < m_totalCount; ++output ) {
            boost::int64_t input;
            if( output > lastOutput ) {
                input = lastInput;
                bool done = false;

                while( !done ) {
                    ++input;
                    accumulator += numerator;
                    if( accumulator >= denominator ) {
                        accumulator -= denominator;
                        done = true;
                    }
                }

                lastInput = input;
                lastOutput = output;

            } else {
                accumulator = 0;
                boost::int64_t outputAcc = -1;
                input = -1;

                while( outputAcc != output ) {
                    ++input;
                    accumulator += numerator;
                    if( accumulator >= denominator ) {
                        accumulator -= denominator;
                        ++outputAcc;
                    }
                }

                lastInput = input;
                lastOutput = output;
            }
            m_inputs[output] = input;
        }
    }

    /**
     * Gets the actual index for the nth element in the list
     */
    boost::int64_t get_input_for_output( boost::int64_t output ) {
        if( output >= m_totalCount )
            throw std::out_of_range(
                ( boost::format( "fractional_index_map::get_input_for_output Value %d outside of bounds [ 0, %d ]" ) %
                  output % ( m_totalCount - 1 ) )
                    .str() );

        return m_inputs[output];
    }
    std::size_t size() { return static_cast<std::size_t>( m_totalCount ); }
};

/**
 * This buffer compresses the list of vertices to only include those used by the faces defined by a
 * fractional_index_map.
 */
class vertex_index_map {
  private:
    std::vector<size_t> m_newToOld;
    std::vector<boost::int64_t> m_oldToNew;

  public:
    vertex_index_map( const mesh_channel* channel, fractional_index_map* filter ) {
        std::vector<bool> used;
        used.resize( channel->get_num_elements(), false );
        for( size_t i = 0; i < filter->size(); i++ ) {
            for( size_t j = 0; j < channel->get_num_face_verts( filter->get_input_for_output( i ) ); j++ ) {
                size_t index = channel->get_fv_index( filter->get_input_for_output( i ), j );
                if( index < used.size() )
                    used[index] = true;
            }
        }

        m_oldToNew.resize( channel->get_num_elements(), -1 );
        size_t count = 0;
        for( size_t i = 0; i < used.size(); i++ ) {
            if( used[i] ) {
                m_oldToNew[i] = count;
                m_newToOld.push_back( i );
                count++;
            }
        }
    }

    vertex_index_map( const mesh_interface* mesh, fractional_index_map* filter ) {
        std::vector<bool> used;
        used.resize( mesh->get_num_verts(), false );
        for( size_t i = 0; i < filter->size(); i++ ) {
            for( size_t j = 0; j < mesh->get_num_face_verts( filter->get_input_for_output( i ) ); j++ ) {
                size_t index = mesh->get_face_vert_index( filter->get_input_for_output( i ), j );
                if( index < used.size() )
                    used[index] = true;
            }
        }
        m_oldToNew.resize( mesh->get_num_verts(), -1 );
        int count = 0;
        for( size_t i = 0; i < used.size(); i++ ) {
            if( used[i] ) {
                m_oldToNew[i] = count;
                m_newToOld.push_back( i );
                count++;
            }
        }
    }

    /**
     * Gets the actual index of the vertex
     */
    size_t get_orig_index( size_t index ) {
        if( index >= m_newToOld.size() )
            throw std::out_of_range(
                ( boost::format( "vertex_index_map::get_orig_index Value %d outside of bounds [ 0, %d ]" ) % index %
                  ( m_newToOld.size() - 1 ) )
                    .str() );
        return m_newToOld[index];
    }

    size_t get_num_verts() { return m_newToOld.size(); }

    /**
     * Gets the index of the vertex in the list of used vertices
     */
    size_t get_new_index( size_t index ) {
        if( index >= m_oldToNew.size() )
            throw std::out_of_range(
                ( boost::format( "vertex_index_map::get_new_index Value %d outside of bounds [ 0, %d ]" ) % index %
                  ( m_oldToNew.size() - 1 ) )
                    .str() );
        boost::int64_t result = m_oldToNew[index];
        if( result < 0 )
            throw std::runtime_error(
                ( boost::format(
                      "vertex_index_map::get_new_index Input value %d does not have a corresponding output value." ) %
                  index )
                    .str() );
        return m_oldToNew[index];
    }
};

class fractional_vertex_channel : public mesh_channel {
  public:
    fractional_vertex_channel( const mesh_channel* channel, boost::shared_ptr<fractional_index_map> vertexIndexMap )
        : mesh_channel( channel->get_name(), channel->get_channel_type(), channel->get_data_type(),
                        channel->get_data_arity(), vertexIndexMap->size(), 0, true )
        , m_channel( channel )
        , m_vertexIndexMap( vertexIndexMap ) {}

    virtual void get_value( std::size_t index, void* outValue ) const {
        m_channel->get_value( m_vertexIndexMap->get_input_for_output( index ), outValue );
    }

    virtual void set_value( std::size_t /*index*/, const void* /*value*/ ) const {
        // should not cause error
    }

    virtual std::size_t get_fv_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/ ) const {
        throw std::runtime_error( "fractional_vertex_channel::get_fv_index() FV index does not exist" );
    }

    virtual std::size_t get_num_face_verts( std::size_t /*faceIndex*/ ) const { return 0; }

  private:
    const mesh_channel* m_channel;
    boost::shared_ptr<fractional_index_map> m_vertexIndexMap;
};

class fractional_face_vertex_channel : public mesh_channel {
  public:
    fractional_face_vertex_channel( const mesh_channel* channel, boost::shared_ptr<fractional_index_map> faceIndexMap,
                                    boost::shared_ptr<vertex_index_map> vertexIndexMap )
        : mesh_channel( channel->get_name(), channel->get_channel_type(), channel->get_data_type(),
                        channel->get_data_arity(), vertexIndexMap->get_num_verts(), faceIndexMap->size(), true )
        , m_channel( channel )
        , m_faceIndexMap( faceIndexMap )
        , m_vertexIndexMap( vertexIndexMap ) {}

    virtual void get_value( std::size_t index, void* outValue ) const {
        m_channel->get_value( m_vertexIndexMap->get_orig_index( index ), outValue );
    }

    virtual void set_value( std::size_t /*index*/, const void* /*value*/ ) const {}

    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        return m_vertexIndexMap->get_new_index(
            m_channel->get_fv_index( m_faceIndexMap->get_input_for_output( faceIndex ), fvertIndex ) );
    }

    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const {
        return m_channel->get_num_face_verts( m_faceIndexMap->get_input_for_output( faceIndex ) );
    }

  private:
    const mesh_channel* m_channel;
    boost::shared_ptr<fractional_index_map> m_faceIndexMap;
    boost::shared_ptr<vertex_index_map> m_vertexIndexMap;
};

class fractional_face_channel : public mesh_channel {
  public:
    fractional_face_channel( const mesh_channel* channel, boost::shared_ptr<fractional_index_map> faceIndexMap )
        : mesh_channel( channel->get_name(), channel->get_channel_type(), channel->get_data_type(),
                        channel->get_data_arity(), faceIndexMap->size(), faceIndexMap->size(), true )
        , m_channel( channel )
        , m_faceIndexMap( faceIndexMap ) {}

    virtual void get_value( std::size_t index, void* outValue ) const {
        m_channel->get_value( m_faceIndexMap->get_input_for_output( index ), outValue );
    }

    virtual void set_value( std::size_t /*index*/, const void* /*value*/ ) const {}

    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        return m_channel->get_fv_index( m_faceIndexMap->get_input_for_output( faceIndex ), fvertIndex );
    }

    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const {
        return m_channel->get_num_face_verts( m_faceIndexMap->get_input_for_output( faceIndex ) );
    }

  private:
    const mesh_channel* m_channel;
    boost::shared_ptr<fractional_index_map> m_faceIndexMap;
};

std::size_t fractional_mesh_interface_base::get_num_elements() const { return m_mesh->get_num_elements(); }

std::size_t fractional_mesh_interface_base::get_face_element_index( std::size_t faceIndex ) const {
    return m_mesh->get_face_element_index( faceIndex );
}

bool fractional_mesh_interface_base::has_adjacency() const { return m_adjacencyDelegate.get() != NULL; }

bool fractional_mesh_interface_base::request_channel( const frantic::tstring& channelName, bool vertexChannel,
                                                      bool forOutput, bool throwOnError ) {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_mesh_interface_base::request_channel() Error: mesh is NULL" );
    }
    if( forOutput ) {
        if( throwOnError ) {
            throw std::runtime_error( "fractional_mesh_interface_base::request_channel() "
                                      "Cannot provide writable channel: "
                                      "\"" +
                                      frantic::strings::to_string( channelName ) + "\"" );
        }
        return false;
    }

    return m_mesh->request_channel( channelName, vertexChannel, forOutput, throwOnError );
}

void fractional_mesh_interface_base::init_adjacency() {
    dcel adjacency;
    mesh_interface_to_dcel( this, adjacency );
    m_adjacencyDelegate.reset( new dcel_mesh_interface( boost::move( adjacency ) ) );
}

bool fractional_mesh_interface_base::is_valid() const { return m_mesh != 0; }

bool fractional_mesh_interface_base::init_vertex_iterator( frantic::geometry::vertex_iterator& vIt,
                                                           std::size_t vertexIndex ) const {
    assert( has_adjacency() && "init_adjacency() must be called before using init_vertex_iterator()" );
    return m_adjacencyDelegate->init_vertex_iterator( vIt, vertexIndex );
}

bool fractional_mesh_interface_base::advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->advance_vertex_iterator( vIt );
}

std::size_t fractional_mesh_interface_base::get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_endpoint( vIt );
}

std::size_t fractional_mesh_interface_base::get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_left_face( vIt );
}

std::size_t fractional_mesh_interface_base::get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_edge_right_face( vIt );
}

bool fractional_mesh_interface_base::is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->is_edge_visible( vIt );
}

bool fractional_mesh_interface_base::is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->is_edge_boundary( vIt );
}

void fractional_mesh_interface_base::init_face_iterator( frantic::geometry::face_iterator& fIt,
                                                         std::size_t faceIndex ) const {
    assert( has_adjacency() && "init_adjacency() must be called before using init_face_iterator()" );
    return m_adjacencyDelegate->init_face_iterator( fIt, faceIndex );
}

bool fractional_mesh_interface_base::advance_face_iterator( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->advance_face_iterator( fIt );
}

std::size_t fractional_mesh_interface_base::get_face_neighbor( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_neighbor( fIt );
}

std::size_t fractional_mesh_interface_base::get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_prev_vertex( fIt );
}

std::size_t fractional_mesh_interface_base::get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const {
    assert( has_adjacency() );
    return m_adjacencyDelegate->get_face_next_vertex( fIt );
}

fractional_face_interface::fractional_face_interface( mesh_interface* mesh, double fraction, boost::int64_t limit ) {
    set_mesh( mesh, fraction, limit );
}

void fractional_face_interface::set_mesh( mesh_interface* mesh, double fraction, boost::int64_t limit ) {
    if( !mesh ) {
        throw std::runtime_error( "fractional_face_interface::set_mesh Error: mesh is NULL" );
    }
    m_faceIndexMap =
        boost::shared_ptr<fractional_index_map>( new fractional_index_map( mesh->get_num_faces(), fraction, limit ) );
    m_mesh = mesh;
    m_fraction = fraction;
    m_limit = limit;
    m_vertexIndexMap = boost::shared_ptr<vertex_index_map>( new vertex_index_map( m_mesh, m_faceIndexMap.get() ) );

    const mesh_channel_map& vertexChannels = m_mesh->get_vertex_channels();
    for( mesh_channel_map::const_iterator channel( vertexChannels.begin() ); channel != vertexChannels.end();
         channel++ ) {
        boost::shared_ptr<vertex_index_map> chBuffer =
            boost::shared_ptr<vertex_index_map>( new vertex_index_map( channel->second, m_faceIndexMap.get() ) );
        std::unique_ptr<fractional_face_vertex_channel> ch(
            new fractional_face_vertex_channel( channel->second, m_faceIndexMap, chBuffer ) );

        this->append_vertex_channel( std::move( ch ) );
    }

    const mesh_channel_map& faceChannels = m_mesh->get_face_channels();
    for( mesh_channel_map::const_iterator channel( faceChannels.begin() ); channel != faceChannels.end(); channel++ ) {
        std::unique_ptr<fractional_face_channel> ch( new fractional_face_channel( channel->second, m_faceIndexMap ) );

        this->append_face_channel( std::move( ch ) );
    }
}

void fractional_face_interface::reset_mesh() {
    reset();
    set_mesh( m_mesh, m_fraction, m_limit );
}

std::size_t fractional_face_interface::get_num_verts() const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_face_interface::get_num_verts() mesh is NULL" );
    }
    return m_vertexIndexMap->get_num_verts();
}

void fractional_face_interface::get_vert( std::size_t index, float ( &outValues )[3] ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_face_interface::get_vert() mesh is NULL" );
    }
    m_mesh->get_vert( m_vertexIndexMap->get_orig_index( index ), outValues );
}

std::size_t fractional_face_interface::get_num_faces() const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_face_interface::get_num_faces() mesh is NULL" );
    }
    return m_faceIndexMap->size();
}

std::size_t fractional_face_interface::get_num_face_verts( std::size_t faceIndex ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_face_interface::get_num_face_verts() mesh is NULL" );
    }
    return m_mesh->get_num_face_verts( m_faceIndexMap->get_input_for_output( faceIndex ) );
}

std::size_t fractional_face_interface::get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_face_interface::get_face_vert_index() mesh is NULL" );
    }
    return m_vertexIndexMap->get_new_index(
        m_mesh->get_face_vert_index( m_faceIndexMap->get_input_for_output( faceIndex ), fvertIndex ) );
}

void fractional_face_interface::get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_face_interface::get_face_vert_indicies() mesh is NULL" );
    }
    m_mesh->get_face_vert_indices( m_faceIndexMap->get_input_for_output( faceIndex ), outValues );
    for( size_t i = 0; i < m_mesh->get_num_face_verts( m_faceIndexMap->get_input_for_output( faceIndex ) ); i++ ) {
        outValues[i] = m_vertexIndexMap->get_new_index( outValues[i] );
    }
}

void fractional_face_interface::get_face_verts( std::size_t faceIndex, float outValues[][3] ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_face_interface::get_face_verts() mesh is NULL" );
    }
    m_mesh->get_face_verts( m_faceIndexMap->get_input_for_output( faceIndex ), outValues );
}

fractional_vertex_interface::fractional_vertex_interface( mesh_interface* mesh, double fraction,
                                                          boost::int64_t limit ) {
    set_mesh( mesh, fraction, limit );
}

void fractional_vertex_interface::set_mesh( mesh_interface* mesh, double fraction, boost::int64_t limit ) {
    if( !mesh ) {
        throw std::runtime_error( "fractional_vertex_interface::set_mesh Error: mesh is NULL" );
    }
    m_vertexIndexMap =
        boost::shared_ptr<fractional_index_map>( new fractional_index_map( mesh->get_num_verts(), fraction, limit ) );
    m_mesh = mesh;
    m_fraction = fraction;
    m_limit = limit;
    const mesh_channel_map& channels = m_mesh->get_vertex_channels();
    for( mesh_channel_map::const_iterator channel = channels.begin(); channel != channels.end(); channel++ ) {
        if( channel->second->get_channel_type() == mesh_channel::vertex ) {
            std::unique_ptr<fractional_vertex_channel> ch(
                new fractional_vertex_channel( channel->second, m_vertexIndexMap ) );

            this->append_vertex_channel( std::move( ch ) );
        }
    }
}

void fractional_vertex_interface::reset_mesh() {
    reset();
    set_mesh( m_mesh, m_fraction, m_limit );
}

bool fractional_vertex_interface::request_channel( const frantic::tstring& channelName, bool vertexChannel,
                                                   bool forOutput, bool throwOnError ) {
    if( !vertexChannel ) {
        if( throwOnError )
            throw std::runtime_error(
                "fractional_vertex_interface::request_channel() Failed to add vertex channel: \"" +
                frantic::strings::to_string( channelName ) + "\"" );

        return false;
    }

    return fractional_mesh_interface_base::request_channel( channelName, vertexChannel, forOutput, throwOnError );
}

std::size_t fractional_vertex_interface::get_num_verts() const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_vertex_interface::get_num_verts() mesh is NULL" );
    }

    return m_vertexIndexMap->size();
}

void fractional_vertex_interface::get_vert( std::size_t index, float ( &outValues )[3] ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_vertex_interface::get_vert() mesh is NULL" );
    }

    m_mesh->get_vert( m_vertexIndexMap->get_input_for_output( index ), outValues );
}

std::size_t fractional_vertex_interface::get_num_faces() const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_vertex_interface::get_num_faces() mesh is NULL" );
    }

    return 0;
}

std::size_t fractional_vertex_interface::get_num_face_verts( std::size_t /*faceIndex*/ ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_vertex_interface::get_num_face_verts() mesh is NULL" );
    }

    return 0;
}

std::size_t fractional_vertex_interface::get_face_vert_index( std::size_t /*faceIndex*/,
                                                              std::size_t /*fvertIndex*/ ) const {
    if( !m_mesh ) {
        throw std::runtime_error( "fractional_vertex_interface::get_face_vert_index() mesh is NULL" );
    }
    throw std::runtime_error( "fractional_vertex_interface::get_face_vert_index() face does not exist" );
}

void fractional_vertex_interface::get_face_vert_indices( std::size_t /*faceIndex*/,
                                                         std::size_t /*outValues*/[] ) const {
    throw std::runtime_error( "fractional_vertex_interface::get_face_vert_indices() face does not exist" );
}

void fractional_vertex_interface::get_face_verts( std::size_t /*faceIndex*/, float /*outValues*/[][3] ) const {
    throw std::runtime_error( "fractional_vertex_interface::get_face_verts() face does not exist" );
}
} // namespace geometry
} // namespace frantic
