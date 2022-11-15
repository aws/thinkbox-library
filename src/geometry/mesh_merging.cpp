// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/foreach.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <frantic/geometry/mesh_merging.hpp>

#include <frantic/geometry/edge_topology.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/logging/progress_logger.hpp>

using std::map;
using std::pair;
using std::size_t;
using std::vector;

using namespace frantic;
using namespace geometry::topology;
using geometry::mesh_channel;
using geometry::mesh_interface;
using geometry::mesh_interface_ptr;
using geometry::polymesh3_builder;
using geometry::polymesh3_ptr;
using graphics::raw_byte_buffer;
using graphics::vector3f;

namespace {

/**
 * This arbiter is a replacement for the old ClosestVertexCalculator from polymesh3.cpp. Its implementation is taken
 * directly from ClosestVertexCalculator, this could be improved in the future to use a cleaner algorithm without
 * multi-level map's.
 */
class duplicate_vertex_arbiter {
  private:
    vector3f m_offset;
    vector3f m_scale;
    float m_tolerance;

    typedef vector<size_t> vert_list_t;
    typedef map<size_t, vert_list_t> map_z_to_vert_list_t;
    typedef map<size_t, map_z_to_vert_list_t> map_y_to_map_z_t;
    typedef map<size_t, map_y_to_map_z_t> map_x_to_map_y_t;
    typedef map_x_to_map_y_t vertex_by_position_t;

    vector<size_t> m_vertexRemap;
    vector<bool> m_vertexMask;
    vertex_by_position_t m_verticesByPosition;
    vector<vector3f> m_vertices;

  public:
    duplicate_vertex_arbiter( vector3f offset, vector3f scale, float tolerance )
        : m_offset( offset )
        , m_scale( scale )
        , m_tolerance( tolerance ) {}

    void begin_mesh( size_t /*mesh*/ ) {}

    bool add_vertex( const vector3f& vertex ) {
        size_t x, y, z;
        to_index( vertex, x, y, z );

        // Get the relevant subgroups to check
        vector<vector<size_t>*> toCheck;
        vector<size_t>* currentGroup;
        vector<size_t>& centerGroup = get_vertex_group_at( x, y, z );

        for( size_t i = 0; i < 3; ++i ) {
            if( x + i < 1 )
                continue;
            for( size_t j = 0; j < 3; ++j ) {
                if( y + j < 1 )
                    continue;
                for( size_t k = 0; k < 3; ++k ) {
                    if( z + k < 1 )
                        continue;

                    if( i == 1 && j == 1 && k == 1 )
                        toCheck.push_back( &centerGroup );
                    else if( try_get_vertex_group_at( x + i - 1, y + j - 1, z + k - 1, currentGroup ) )
                        toCheck.push_back( currentGroup );
                }
            }
        }

        // Check for closest match
        const size_t nextIndex = m_vertices.size();
        size_t bestIndex = nextIndex;
        float bestError = std::numeric_limits<float>::max();
        for( vector<vector<size_t>*>::const_iterator iter = toCheck.begin(); iter != toCheck.end(); ++iter ) {
            for( vector<size_t>::const_iterator iter2 = ( *iter )->begin(); iter2 != ( *iter )->end(); ++iter2 ) {
                const vector3f& checkPoint = m_vertices[*iter2];
                const float error = vector3f::distance_squared( checkPoint, vertex );
                if( error < bestError ) {
                    bestError = error;
                    bestIndex = *iter2;
                }
            }
        }

        if( bestIndex >= m_vertices.size() || bestError > m_tolerance * m_tolerance ) {
            // No closest match, add as new vertex
            m_vertices.push_back( vertex );
            centerGroup.push_back( nextIndex );
            m_vertexRemap.push_back( nextIndex );
            m_vertexMask.push_back( true );
            return true;
        } else {
            // Closest match found, remap the vertex
            m_vertexRemap.push_back( bestIndex );
            m_vertexMask.push_back( false );
            return false;
        }
    }

    size_t vertex_remap( size_t index ) const { return m_vertexRemap[index]; }

    bool vertex_mask( size_t index ) const { return m_vertexMask[index]; }

  private:
    void to_index( const vector3f& vertex, size_t& outX, size_t& outY, size_t& outZ ) const {
        outX = static_cast<size_t>( 0.5f + ( vertex.x - m_offset.x ) / m_scale.x );
        outY = static_cast<size_t>( 0.5f + ( vertex.y - m_offset.y ) / m_scale.y );
        outZ = static_cast<size_t>( 0.5f + ( vertex.z - m_offset.z ) / m_scale.z );
    }

    vector<size_t>& get_vertex_group_at( size_t x, size_t y, size_t z ) { return m_verticesByPosition[x][y][z]; }

    bool try_get_vertex_group_at( size_t x, size_t y, size_t z, vector<size_t>*& result ) {
        map_x_to_map_y_t::iterator iterX = m_verticesByPosition.find( x );
        if( iterX != m_verticesByPosition.end() ) {
            map_y_to_map_z_t::iterator iterY = iterX->second.find( y );
            if( iterY != iterX->second.end() ) {
                map_z_to_vert_list_t::iterator iterZ = iterY->second.find( z );
                if( iterZ != iterY->second.end() ) {
                    result = &iterZ->second;
                    return true;
                }
            }
        }
        return false;
    }
};

class vertex_cluster_arbiter {
  private:
    vector<size_t> m_vertexRemap;
    vector<bool> m_vertexMask;

    vector<pair<bool, size_t>> m_clusterReps;
    vector<vector<pair<size_t, size_t>>> m_vertexToClusterMaps;
    vector<pair<size_t, size_t>>::iterator m_nextClusteredVertex;

    size_t m_currentMesh;
    size_t m_vertexIndex;
    size_t m_nextVertex;

  public:
    vertex_cluster_arbiter( vector<vector<mesh_vertex>> vertexClusters, size_t meshCount )
        : m_clusterReps( vertexClusters.size() )
        , m_vertexToClusterMaps( meshCount )
        , m_nextVertex( 0 ) {
        for( size_t clusterIndex = 0, clusterCount = vertexClusters.size(); clusterIndex < clusterCount;
             ++clusterIndex ) {
            BOOST_FOREACH( mesh_vertex vertex, vertexClusters[clusterIndex] ) {
                m_vertexToClusterMaps[vertex.m_mesh].push_back( pair<size_t, size_t>( vertex.m_index, clusterIndex ) );
            }
        }

        for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
            std::sort( m_vertexToClusterMaps[meshIndex].begin(), m_vertexToClusterMaps[meshIndex].end() );
        }
    }

    void begin_mesh( size_t mesh ) {
        m_currentMesh = mesh;
        m_nextClusteredVertex = m_vertexToClusterMaps[mesh].begin();
        m_vertexIndex = 0;
    }

    bool add_vertex( const vector3f& /*vertex*/ ) {
        bool include = true;
        size_t clusterIndex;

        if( m_nextClusteredVertex != m_vertexToClusterMaps[m_currentMesh].end() &&
            m_nextClusteredVertex->first == m_vertexIndex ) {
            if( m_clusterReps[m_nextClusteredVertex->second].first ) {
                // A vertex in the cluster has already been added, reference it
                clusterIndex = m_clusterReps[m_nextClusteredVertex->second].second;
                include = false;
            } else {
                // Add yourself as the first vertex in the cluster
                m_clusterReps[m_nextClusteredVertex->second].first = true;
                m_clusterReps[m_nextClusteredVertex->second].second = m_nextVertex;
            }
            ++m_nextClusteredVertex;
        }

        if( include ) {
            m_vertexRemap.push_back( m_nextVertex++ );
            m_vertexMask.push_back( true );
        } else {
            m_vertexRemap.push_back( clusterIndex );
            m_vertexMask.push_back( false );
        }

        ++m_vertexIndex;
        return include;
    }

    size_t vertex_remap( size_t index ) const { return m_vertexRemap[index]; }

    bool vertex_mask( size_t index ) const { return m_vertexMask[index]; }
};

/**
 * NOTE: The following data structures and static helper functions have been moved from their original location in
 * polymesh3.cpp, where they supported the original weld and combine functions. Their parameters and function calls have
 * been touched up to use mesh_interface_ptr instead of polymesh3_ptr but their structure has not been redesigned to
 * take advantage of the mesh_interface design. They should be re-evaluated at some point to take advantage of any
 * mesh_interface functions or utilities.
 */

struct combine_channel_info {
  private:
    size_t m_hitCount;
    size_t m_customFacesHitCount;
    bool m_isVertexChannel;
    bool m_isFaceChannel;
    size_t m_arity;
    channels::data_type_t m_dataType;
    bool m_isEmpty;

    void reset() {
        m_isEmpty = true;
        m_hitCount = 0;
        m_customFacesHitCount = 0;
        m_isVertexChannel = false;
        m_isFaceChannel = false;
        m_arity = 0;
        m_dataType = channels::data_type_invalid;
    }

  public:
    combine_channel_info() { reset(); }

    void add_vertex_channel( const tstring& channelName, channels::data_type_t dataType, size_t arity,
                             bool hasCustomFaces = false ) {
        if( m_isFaceChannel ) {
            throw std::runtime_error( "combine_channel_info.add_vertex_channel Error: channel \'" +
                                      strings::to_string( channelName ) + "\' was already added as a face channel." );
        }
        m_isVertexChannel = true;
        ++m_hitCount;
        if( hasCustomFaces ) {
            ++m_customFacesHitCount;
        }
        if( m_isEmpty ) {
            m_dataType = dataType;
            m_arity = arity;
        } else {
            m_dataType = channels::promote_types( m_dataType, dataType );
            if( m_arity != arity ) {
                throw std::runtime_error(
                    "combine_channel_info.add_vertex_channel Error: arity mismatch in channel \'" +
                    strings::to_string( channelName ) + "\'.  (" + boost::lexical_cast<std::string>( m_arity ) +
                    " vs " + boost::lexical_cast<std::string>( arity ) + ")" );
            }
        }
        m_isEmpty = false;
    }

    void add_face_channel( const tstring& channelName, channels::data_type_t dataType, size_t arity ) {
        if( m_isVertexChannel ) {
            throw std::runtime_error( "combine_channel_info.add_face_channel Error: channel \'" +
                                      strings::to_string( channelName ) + "\' was already added as a vertex channel." );
        }
        m_isFaceChannel = true;
        ++m_hitCount;
        if( m_isEmpty ) {
            m_dataType = dataType;
            m_arity = arity;
        } else {
            m_dataType = channels::promote_types( m_dataType, dataType );
            if( m_arity != arity ) {
                throw std::runtime_error( "combine_channel_info.add_face_channel Error: arity mismatch in channel \'" +
                                          strings::to_string( channelName ) + "\'.  (" +
                                          boost::lexical_cast<std::string>( m_arity ) + " vs " +
                                          boost::lexical_cast<std::string>( arity ) + ")" );
            }
        }
        m_isEmpty = false;
    }

    bool is_face_channel() const { return m_isFaceChannel; }

    channels::data_type_t get_data_type() const { return m_dataType; }

    size_t get_arity() const { return m_arity; }

    bool has_custom_faces() const { return m_customFacesHitCount > 0; }

    size_t get_hit_count() const { return m_hitCount; }
};
} // anonymous namespace

static size_t get_vertex_channel_element_count_sum( const mesh_interface* mesh, const tstring& channelName ) {
    if( mesh ) {
        const mesh_channel* channel = mesh->get_vertex_channel( channelName );
        if( channel ) {
            return channel->get_num_elements();
        }
    }
    return 0;
}

static size_t get_vertex_channel_element_count_sum( const vector<const mesh_interface*>& meshes,
                                                    const tstring& channelName ) {
    size_t count = 0;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        count += get_vertex_channel_element_count_sum( mesh, channelName );
    }
    return count;
}

static size_t get_face_count_sum( const vector<const mesh_interface*>& meshes ) {
    size_t count = 0;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( mesh ) {
            count += mesh->get_num_faces();
        }
    }
    return count;
}

static polymesh3_ptr combine_geometry( const vector<const mesh_interface*>& meshes ) {
    geometry::polymesh3_builder builder;

    size_t vertexIndexOffset = 0;

    // Create geometry
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }

        const size_t vertexCount = mesh->get_num_verts();
        const size_t faceCount = mesh->get_num_faces();

        for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
            builder.add_vertex( mesh->get_vert( vertexIndex ) );
        }

        vector<size_t> inFace;
        vector<int> outFace;
        for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
            const size_t fvCount = mesh->get_num_face_verts( faceIndex );
            inFace.resize( fvCount );
            outFace.resize( fvCount );

            mesh->get_face_vert_indices( faceIndex, &inFace[0] );
            for( size_t fvIndex = 0; fvIndex < fvCount; ++fvIndex ) {
                outFace[fvIndex] = static_cast<int>( vertexIndexOffset + inFace[fvIndex] );
            }

            builder.add_polygon( outFace );
        }

        vertexIndexOffset += vertexCount;
    }

    geometry::polymesh3_ptr result = builder.finalize();

    return result;
}

static bool is_non_empty_mesh( const mesh_interface* mesh ) {
    if( mesh ) {
        if( mesh->get_num_verts() > 0 ) {
            return true;
        }
    }
    return false;
}

static size_t get_non_empty_mesh_count( const vector<const mesh_interface*>& meshes ) {
    size_t count = 0;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( is_non_empty_mesh( mesh ) ) {
            ++count;
        }
    }
    return count;
}

static void get_combined_channel_info( map<tstring, combine_channel_info>& outInfo, const mesh_interface* mesh ) {
    if( mesh->get_num_verts() ) {
        const mesh_interface::mesh_channel_map& vertexChannels = mesh->get_vertex_channels();
        for( mesh_interface::mesh_channel_map::const_iterator channel = vertexChannels.begin(),
                                                              channelEnd = vertexChannels.end();
             channel != channelEnd; ++channel ) {
            combine_channel_info& info = outInfo[channel->first];
            info.add_vertex_channel( channel->first, channel->second->get_data_type(),
                                     channel->second->get_data_arity(),
                                     channel->second->get_channel_type() == mesh_channel::face_vertex );
        }
        const mesh_interface::mesh_channel_map& faceChannels = mesh->get_face_channels();
        for( mesh_interface::mesh_channel_map::const_iterator channel = faceChannels.begin(),
                                                              channelEnd = faceChannels.end();
             channel != channelEnd; ++channel ) {
            combine_channel_info& info = outInfo[channel->first];
            info.add_face_channel( channel->first, channel->second->get_data_type(),
                                   channel->second->get_data_arity() );
        }
    }
}

static void get_combined_channel_info( map<tstring, combine_channel_info>& outInfo,
                                       const vector<const mesh_interface*>& meshes ) {
    outInfo.clear();

    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }
        get_combined_channel_info( outInfo, mesh );
    }
}

static void combine_face_channel( polymesh3_ptr& result, const vector<const mesh_interface*>& meshes,
                                  const tstring& channelName, const combine_channel_info& channelInfo ) {
    const size_t faceCountSum = get_face_count_sum( meshes );
    const size_t elementSize =
        channelInfo.get_arity() * channels::sizeof_channel_data_type( channelInfo.get_data_type() );

    raw_byte_buffer data;
    data.resize( faceCountSum * elementSize );

    size_t dataIndexOffset = 0;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }
        const size_t faceCount = mesh->get_num_faces();
        if( faceCount == 0 ) {
            continue;
        }

        const mesh_channel* channel = mesh->get_face_channel( channelName );
        if( channel ) {
            if( channel->get_data_arity() != channelInfo.get_arity() ) {
                throw std::runtime_error( "combine_face_channel Error: arity mismatch" );
            }

            if( channel->get_data_type() != channelInfo.get_data_type() ) {
                const size_t inputElementSize = channel->get_element_size();
                raw_byte_buffer elementBuffer;
                elementBuffer.resize( faceCount * inputElementSize );

                for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                    channel->get_value( faceIndex, elementBuffer.ptr_at( faceIndex * inputElementSize ) );
                }

                channels::channel_type_convertor_function_t cvt = channels::get_channel_type_convertor_function(
                    channel->get_data_type(), channelInfo.get_data_type(), channelName );
                cvt( data.ptr_at( dataIndexOffset * elementSize ), elementBuffer.ptr_at( 0 ),
                     faceCount * channelInfo.get_arity() );
            } else {
                for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                    channel->get_value( faceIndex, data.ptr_at( ( dataIndexOffset + faceIndex ) * elementSize ) );
                }
            }
        } else {
            memset( data.ptr_at( dataIndexOffset * elementSize ), 0, faceCount * elementSize );
        }

        dataIndexOffset += faceCount;
    }

    result->add_face_channel( channelName, channelInfo.get_data_type(), channelInfo.get_arity(), data );
}

static void combine_custom_vertex_channel( polymesh3_ptr& result, const vector<const mesh_interface*>& meshes,
                                           const tstring& channelName, const combine_channel_info& channelInfo ) {
    const size_t elementSize =
        channelInfo.get_arity() * channels::sizeof_channel_data_type( channelInfo.get_data_type() );
    size_t dataCountSum = get_vertex_channel_element_count_sum( meshes, channelName );

    const bool createDefaultData = ( channelInfo.get_hit_count() < get_non_empty_mesh_count( meshes ) );
    if( createDefaultData ) {
        ++dataCountSum;
    }

    size_t dataIndexOffset = 0;
    raw_byte_buffer data;
    data.resize( dataCountSum * elementSize );

    if( createDefaultData ) {
        memset( data.ptr_at( 0 ), 0, elementSize );
        ++dataIndexOffset;
    }

    vector<int> customFaces;

    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }
        if( mesh->get_num_verts() == 0 ) {
            continue;
        }

        const mesh_channel* channel = mesh->get_vertex_channel( channelName );
        if( channel ) {
            // copy data
            const size_t elementCount = channel->get_num_elements();
            if( channel->get_data_type() != channelInfo.get_data_type() ) {
                const size_t inputElementSize = channel->get_element_size();
                raw_byte_buffer elementBuffer;
                elementBuffer.resize( elementCount * inputElementSize );

                for( size_t elementIndex = 0; elementIndex < elementCount; ++elementIndex ) {
                    channel->get_value( elementIndex, elementBuffer.ptr_at( elementIndex * inputElementSize ) );
                }

                channels::channel_type_convertor_function_t cvt = channels::get_channel_type_convertor_function(
                    channel->get_data_type(), channelInfo.get_data_type(), channelName );
                cvt( data.ptr_at( dataIndexOffset * elementSize ), elementBuffer.ptr_at( 0 ),
                     elementCount * channelInfo.get_arity() );
            } else {
                for( size_t elementIndex = 0; elementIndex < elementCount; ++elementIndex ) {
                    channel->get_value( elementIndex, data.ptr_at( ( dataIndexOffset + elementIndex ) * elementSize ) );
                }
            }

            // copy faces with offset
            if( channel->get_channel_type() == mesh_channel::face_vertex ) {
                const size_t faceCount = channel->get_num_faces();
                for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                    for( size_t fvIndex = 0, fvCount = channel->get_num_face_verts( faceIndex ); fvIndex < fvCount;
                         ++fvIndex ) {
                        const size_t dataIndex = channel->get_fv_index( faceIndex, fvIndex ) + dataIndexOffset;
                        try {
                            customFaces.push_back( boost::numeric_cast<int>( dataIndex ) );
                        } catch( boost::bad_numeric_cast& ) {
                            throw std::runtime_error( "combine_custom_vertex_channel Error: face index overflow" );
                        }
                    }
                }
            } else {
                const size_t faceCount = mesh->get_num_faces();
                for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                    for( size_t fvIndex = 0, fvCount = mesh->get_num_face_verts( faceIndex ); fvIndex < fvCount;
                         ++fvIndex ) {
                        const size_t dataIndex = mesh->get_face_vert_index( faceIndex, fvIndex ) + dataIndexOffset;
                        try {
                            customFaces.push_back( boost::numeric_cast<int>( dataIndex ) );
                        } catch( boost::bad_numeric_cast& ) {
                            throw std::runtime_error( "combine_custom_vertex_channel Error: face index overflow" );
                        }
                    }
                }
            }

            dataIndexOffset += elementCount;
        } else {
            if( !createDefaultData ) {
                throw std::runtime_error(
                    "combine_vertex_channel Internal Error: unprepared for mesh without channel" );
            }

            // fill custom faces with vertex index 0
            // the 0th data entry was already set to the default value (0) earlier
            for( size_t faceIndex = 0, faceCount = mesh->get_num_faces(); faceIndex < faceCount; ++faceIndex ) {
                for( size_t corner = 0, cornerCount = mesh->get_num_face_verts( faceIndex ); corner < cornerCount;
                     ++corner ) {
                    customFaces.push_back( 0 );
                }
            }
        }
    }

    result->add_vertex_channel( channelName, channelInfo.get_data_type(), channelInfo.get_arity(), data, &customFaces );
}

static void combine_simple_vertex_channel( polymesh3_ptr& result, const vector<const mesh_interface*>& meshes,
                                           const tstring& channelName, const combine_channel_info& channelInfo ) {
    const size_t elementSize =
        channelInfo.get_arity() * channels::sizeof_channel_data_type( channelInfo.get_data_type() );
    const size_t dataCountSum = result->vertex_count();

    raw_byte_buffer data;
    data.resize( elementSize * dataCountSum );

    size_t dataIndexOffset = 0;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }
        const size_t vertexCount = mesh->get_num_verts();
        if( vertexCount == 0 ) {
            continue;
        }

        const mesh_channel* channel = mesh->get_vertex_channel( channelName );
        if( channel ) {
            if( channel->get_data_arity() != channelInfo.get_arity() ) {
                throw std::runtime_error( "combine_simple_vertex_channel Error: arity mismatch" );
            }

            if( channel->get_data_type() != channelInfo.get_data_type() ) {
                const size_t inputElementSize = channel->get_element_size();
                raw_byte_buffer elementBuffer;
                elementBuffer.resize( vertexCount * inputElementSize );

                for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
                    channel->get_value( vertexIndex, elementBuffer.ptr_at( vertexIndex * inputElementSize ) );
                }

                channels::channel_type_convertor_function_t cvt = channels::get_channel_type_convertor_function(
                    channel->get_data_type(), channelInfo.get_data_type(), channelName );
                cvt( data.ptr_at( dataIndexOffset * elementSize ), elementBuffer.ptr_at( 0 ),
                     vertexCount * channelInfo.get_arity() );
            } else {
                for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
                    channel->get_value( vertexIndex, data.ptr_at( ( dataIndexOffset + vertexIndex ) * elementSize ) );
                }
            }
        } else {
            memset( data.ptr_at( dataIndexOffset * elementSize ), 0, vertexCount * elementSize );
        }

        dataIndexOffset += vertexCount;
    }

    result->add_vertex_channel( channelName, channelInfo.get_data_type(), channelInfo.get_arity(), data );
}

static void combine_vertex_channel( polymesh3_ptr& result, const vector<const mesh_interface*>& meshes,
                                    const tstring& channelName, const combine_channel_info& channelInfo ) {
    if( channelInfo.has_custom_faces() ) {
        combine_custom_vertex_channel( result, meshes, channelName, channelInfo );
    } else {
        combine_simple_vertex_channel( result, meshes, channelName, channelInfo );
    }
}

static void combine_channels( polymesh3_ptr& result, const vector<const mesh_interface*>& meshes ) {
    typedef map<tstring, combine_channel_info> combined_channel_info_t;
    combined_channel_info_t combinedChannelInfo;
    get_combined_channel_info( combinedChannelInfo, meshes );

    for( combined_channel_info_t::iterator i = combinedChannelInfo.begin(); i != combinedChannelInfo.end(); ++i ) {
        if( i->second.is_face_channel() ) {
            combine_face_channel( result, meshes, i->first, i->second );
        } else {
            combine_vertex_channel( result, meshes, i->first, i->second );
        }
    }
}

template <class WeldingArbiter>
static void combine_simple_vertex_channel_welded( polymesh3_ptr& result, const vector<const mesh_interface*>& meshes,
                                                  const tstring& channelName, const combine_channel_info& channelInfo,
                                                  WeldingArbiter arbiter ) {
    const size_t elementSize =
        channelInfo.get_arity() * channels::sizeof_channel_data_type( channelInfo.get_data_type() );
    const size_t dataCountSum = result->vertex_count();

    raw_byte_buffer data;
    data.resize( elementSize * dataCountSum );

    size_t dataIndexOffset = 0;
    size_t vertexIndexOffset = 0;

    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        // TODO: blend channels if vertex merged
        if( !mesh ) {
            continue;
        }
        const size_t vertexCount = mesh->get_num_verts();
        if( vertexCount == 0 ) {
            continue;
        }

        size_t currentIndex = 0;
        size_t currentCount = 0;
        size_t currentOffset = 0;
        size_t currentTotal = 0;
        const mesh_channel* channel = mesh->get_vertex_channel( channelName );
        if( channel ) {
            channels::channel_type_convertor_function_t cvt = channels::get_channel_type_convertor_function(
                channel->get_data_type(), channelInfo.get_data_type(), channelName );

            const size_t sourceVertexCountSum = get_vertex_channel_element_count_sum( mesh, channelName );
            for( size_t i = 0; i < sourceVertexCountSum; ++i ) {
                const size_t index = i + vertexIndexOffset;
                if( !arbiter.vertex_mask( index ) ) {
                    if( currentCount > 0 ) {
                        if( channel->get_data_type() != channelInfo.get_data_type() ) {
                            const size_t inputElementSize = channel->get_element_size();
                            raw_byte_buffer elementBuffer;
                            elementBuffer.resize( currentCount * inputElementSize );

                            for( size_t subsetIndex = 0; subsetIndex < currentCount; ++subsetIndex ) {
                                channel->get_value( currentIndex + subsetIndex,
                                                    elementBuffer.ptr_at( subsetIndex * inputElementSize ) );
                            }

                            cvt( data.ptr_at( ( dataIndexOffset + currentOffset ) * elementSize ),
                                 elementBuffer.ptr_at( 0 ), currentCount * channelInfo.get_arity() );
                        } else {
                            for( size_t subsetIndex = 0; subsetIndex < currentCount; ++subsetIndex ) {
                                channel->get_value(
                                    currentIndex + subsetIndex,
                                    data.ptr_at( ( dataIndexOffset + currentOffset + subsetIndex ) * elementSize ) );
                            }
                        }
                    }
                    currentTotal += currentCount;
                    currentOffset += currentCount;
                    currentIndex = i + 1;
                    currentCount = 0;
                } else {
                    ++currentCount;
                }
            }
            if( currentCount > 0 ) {
                currentTotal += currentCount;
                if( channel->get_data_type() != channelInfo.get_data_type() ) {
                    const size_t inputElementSize = channel->get_element_size();
                    raw_byte_buffer elementBuffer;
                    elementBuffer.resize( currentCount * inputElementSize );

                    for( size_t subsetIndex = 0; subsetIndex < currentCount; ++subsetIndex ) {
                        channel->get_value( currentIndex + subsetIndex,
                                            elementBuffer.ptr_at( subsetIndex * inputElementSize ) );
                    }

                    cvt( data.ptr_at( ( dataIndexOffset + currentOffset ) * elementSize ), elementBuffer.ptr_at( 0 ),
                         currentCount * channelInfo.get_arity() );
                } else {
                    for( size_t subsetIndex = 0; subsetIndex < currentCount; ++subsetIndex ) {
                        channel->get_value(
                            currentIndex + subsetIndex,
                            data.ptr_at( ( dataIndexOffset + currentOffset + subsetIndex ) * elementSize ) );
                    }
                }
            }
        } else {
            const size_t sourceVertexCountSum = vertexCount;
            for( size_t i = 0; i < sourceVertexCountSum; ++i ) {
                const size_t index = i + vertexIndexOffset;
                if( !arbiter.vertex_mask( index ) ) {
                    if( currentCount > 0 ) {
                        memset( data.ptr_at( ( dataIndexOffset + currentOffset ) * elementSize ), 0,
                                currentCount * elementSize );
                    }
                    currentTotal += currentCount;
                    currentOffset += currentCount;
                    currentIndex = i + 1;
                    currentCount = 0;
                } else {
                    ++currentCount;
                }
            }
            if( currentCount > 0 ) {
                currentTotal += currentCount;
                memset( data.ptr_at( ( dataIndexOffset + currentOffset ) * elementSize ), 0,
                        currentCount * elementSize );
            }
        }
        dataIndexOffset += currentTotal;
        vertexIndexOffset += vertexCount;
    }
    result->add_vertex_channel( channelName, channelInfo.get_data_type(), channelInfo.get_arity(), data );
}

template <class WeldingArbiter>
static void combine_vertex_channel_welded( polymesh3_ptr& result, const vector<const mesh_interface*>& meshes,
                                           const tstring& channelName, const combine_channel_info& channelInfo,
                                           WeldingArbiter& arbiter ) {
    if( channelInfo.has_custom_faces() ) {
        combine_custom_vertex_channel( result, meshes, channelName, channelInfo );
    } else {
        combine_simple_vertex_channel_welded( result, meshes, channelName, channelInfo, arbiter );
    }
}

template <class WeldingArbiter>
static polymesh3_ptr weld_meshes( const vector<const mesh_interface*>& meshes, WeldingArbiter& arbiter,
                                  logging::progress_logger& progress ) {
    polymesh3_builder builder;

    size_t totalVertices = 0;
    size_t totalFaces = 0;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }
        totalVertices += mesh->get_num_verts();
        totalFaces += mesh->get_num_faces();
    }

    size_t totalProgress = 0;
    size_t progressAccumulator = 0;
    size_t progressThreshold = std::min<size_t>( 10000, std::min( totalVertices / 10, totalFaces / 10 ) );
    logging::progress_logger_subinterval_tracker progressInterval( progress, 0, 100 );

    // Add the vertices
    progress.set_title( _T("Merging Vertices") );
    progressInterval.reset( 0.0f, 75.0f );

    const size_t meshCount = meshes.size();
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        if( !meshes[meshIndex] ) {
            continue;
        }

        arbiter.begin_mesh( meshIndex );

        const size_t vertexCount = meshes[meshIndex]->get_num_verts();
        for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
            const vector3f checkPosition = meshes[meshIndex]->get_vert( vertexIndex );
            if( arbiter.add_vertex( checkPosition ) ) {
                builder.add_vertex( checkPosition );
            }

            ++totalProgress;
            ++progressAccumulator;
            if( progressAccumulator >= progressThreshold ) {
                progress.update_progress( totalProgress, totalVertices );
                progressAccumulator = 0;
                progress.check_for_abort();
            }
        }
        progress.check_for_abort();
    }
    progress.check_for_abort();

    // Add the faces
    progress.set_title( _T("Merging Faces") );
    progressInterval.reset( 75.0f, 87.5f );
    totalProgress = 0;
    progressAccumulator = 0;

    size_t vertexIndexOffset = 0;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }

        const size_t vertexCount = mesh->get_num_verts();
        const size_t faceCount = mesh->get_num_faces();

        vector<size_t> inFace;
        vector<int> outFace;
        for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
            const size_t fvCount = mesh->get_num_face_verts( faceIndex );
            inFace.resize( fvCount );
            outFace.resize( fvCount );

            mesh->get_face_vert_indices( faceIndex, &inFace[0] );
            for( size_t fvIndex = 0; fvIndex < fvCount; ++fvIndex ) {
                outFace[fvIndex] = static_cast<int>( arbiter.vertex_remap( vertexIndexOffset + inFace[fvIndex] ) );
            }

            builder.add_polygon( outFace );

            ++totalProgress;
            ++progressAccumulator;
            if( progressAccumulator >= progressThreshold ) {
                progress.update_progress( totalProgress, totalFaces );
                progressAccumulator = 0;
                progress.check_for_abort();
            }
        }

        vertexIndexOffset += vertexCount;
        progress.check_for_abort();
    }
    progress.check_for_abort();

    // Finalize base geometry
    polymesh3_ptr result = builder.finalize();

    // Merge channels
    progress.set_title( _T("Generating Geometry") );
    progressInterval.reset( 87.5f, 100.0f );
    totalProgress = 0;

    typedef map<frantic::tstring, combine_channel_info> combined_channel_info_t;
    combined_channel_info_t combinedChannelInfo;
    combinedChannelInfo.clear();
    get_combined_channel_info( combinedChannelInfo, meshes );

    const size_t totalChannels = combinedChannelInfo.size();
    for( combined_channel_info_t::iterator i = combinedChannelInfo.begin(); i != combinedChannelInfo.end(); ++i ) {
        if( i->second.is_face_channel() ) {
            combine_face_channel( result, meshes, i->first, i->second );
        } else {
            combine_vertex_channel_welded( result, meshes, i->first, i->second, arbiter );
        }

        ++totalProgress;
        progress.update_progress( totalProgress, totalChannels );
        progress.check_for_abort();
    }

    return result;
}

polymesh3_ptr geometry::weld_vertices( const vector<mesh_interface_ptr>& meshes, float tolerance ) {
    logging::null_progress_logger progress;

    const size_t meshCount = meshes.size();
    vector<const mesh_interface*> interfaceRefs( meshCount );
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        interfaceRefs[meshIndex] = meshes[meshIndex].get();
    }

    return weld_vertices( interfaceRefs, tolerance, progress );
}

polymesh3_ptr geometry::weld_vertices( const vector<const mesh_interface*>& meshes, float tolerance ) {
    logging::null_progress_logger progress;
    return weld_vertices( meshes, tolerance, progress );
}

polymesh3_ptr geometry::weld_vertices( const vector<mesh_interface_ptr>& meshes, float tolerance,
                                       logging::progress_logger& progress ) {
    const size_t meshCount = meshes.size();
    vector<const mesh_interface*> interfaceRefs( meshCount );
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        interfaceRefs[meshIndex] = meshes[meshIndex].get();
    }

    return weld_vertices( interfaceRefs, tolerance, progress );
}

polymesh3_ptr geometry::weld_vertices( const vector<const mesh_interface*>& meshes, float tolerance,
                                       logging::progress_logger& progress ) {
    // Calculate bounding box for all meshes
    graphics::boundbox3f bbox;
    BOOST_FOREACH( const mesh_interface* mesh, meshes ) {
        if( !mesh ) {
            continue;
        }

        const size_t vertexCount = mesh->get_num_verts();
        for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
            bbox += mesh->get_vert( vertexIndex );
        }
    }

    // Compute adequate space division
    vector3f offset = bbox.minimum();
    vector3f scale = bbox.maximum() - bbox.minimum();
    if( scale.x <= 0 )
        scale.x = 1;
    if( scale.y <= 0 )
        scale.y = 1;
    if( scale.z <= 0 )
        scale.z = 1;
    if( tolerance > 0 ) {
        if( scale.x > tolerance )
            scale.x = tolerance;
        if( scale.y > tolerance )
            scale.y = tolerance;
        if( scale.z > tolerance )
            scale.z = tolerance;
    } else {
        scale.x /= 32;
        scale.y /= 32;
        scale.z /= 32;
    }

    // Create a welding policy
    duplicate_vertex_arbiter arbiter( offset, scale, tolerance );

    return weld_meshes( meshes, arbiter, progress );
}

typedef pair<mesh_half_edge, mesh_half_edge> mesh_half_edge_pair;

static void compute_clusters( const vector<mesh_half_edge_pair>& edgePairs, const vector<vector<half_edge>>& meshEdges,
                              vector<vector<mesh_vertex>>& outClusters ) {
    // Create set of mesh_vertices from set of mesh_edge pairs, giving each a reference to it's vertex pair
    vector<pair<mesh_vertex, size_t>> meshVertices;
    meshVertices.reserve( edgePairs.size() * 2 );

    // Each edge pair produces two vertex pairs which are stored as 4 adjacent entries in this vector
    vector<size_t> vertexPairs( edgePairs.size() * 4 );

    BOOST_FOREACH( mesh_half_edge_pair edgePair, edgePairs ) {
        mesh_vertex firstHead =
            mesh_vertex( edgePair.first.m_mesh, meshEdges[edgePair.first.m_mesh][edgePair.first.m_edge].m_head );
        mesh_vertex firstTail =
            mesh_vertex( edgePair.first.m_mesh, meshEdges[edgePair.first.m_mesh][edgePair.first.m_edge].m_tail );
        mesh_vertex secondHead =
            mesh_vertex( edgePair.second.m_mesh, meshEdges[edgePair.second.m_mesh][edgePair.second.m_edge].m_head );
        mesh_vertex secondTail =
            mesh_vertex( edgePair.second.m_mesh, meshEdges[edgePair.second.m_mesh][edgePair.second.m_edge].m_tail );

        // Flatten vertex pairs into the vector
        const size_t index = meshVertices.size();
        meshVertices.push_back( pair<mesh_vertex, size_t>( firstHead, index ) );
        meshVertices.push_back( pair<mesh_vertex, size_t>( secondTail, index + 1 ) );
        meshVertices.push_back( pair<mesh_vertex, size_t>( firstTail, index + 2 ) );
        meshVertices.push_back( pair<mesh_vertex, size_t>( secondHead, index + 3 ) );
    }

    // Sort set of mesh_vertices
    std::sort( meshVertices.begin(), meshVertices.end() );

    // Remove duplicates from the set of mesh_vertices, redirecting vertex pair indicies to the new vertex positions
    vertexPairs[meshVertices[0].second] = 0;
    size_t uniqueIndex = 0;
    for( size_t mvIndex = 1, mvCount = meshVertices.size(); mvIndex < mvCount; ++mvIndex ) {
        if( meshVertices[mvIndex].first == meshVertices[mvIndex - 1].first ) {
            vertexPairs[meshVertices[mvIndex].second] = uniqueIndex;
        } else {
            ++uniqueIndex;
            vertexPairs[meshVertices[mvIndex].second] = uniqueIndex;
            meshVertices[uniqueIndex] = meshVertices[mvIndex];
        }
    }
    const size_t vertexCount = uniqueIndex + 1;
    meshVertices.resize( vertexCount );

    // Perform union-find on the set of mesh_vertices, using the vertex pairs to create connections
    vector<size_t> ranks( vertexCount );
    vector<size_t> parents( vertexCount );
    boost::disjoint_sets<size_t*, size_t*> disjointSets( &ranks[0], &parents[0] );

    for( size_t i = 0; i < vertexCount; ++i ) {
        disjointSets.make_set( i );
    }

    for( size_t pairIndex = 0, pairCount = vertexPairs.size(); pairIndex < pairCount; pairIndex += 2 ) {
        disjointSets.union_set( vertexPairs[pairIndex], vertexPairs[pairIndex + 1] );
    }

    // Construct vertex clusters from the union-find structure
    outClusters.clear();
    outClusters.reserve( vertexCount );

    const size_t uninitialized = std::numeric_limits<size_t>::max();

    vector<size_t> clusterIndices( vertexCount, uninitialized );
    for( size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
        size_t rep = disjointSets.find_set( vertexIndex );
        if( clusterIndices[rep] != uninitialized ) {
            outClusters[clusterIndices[rep]].push_back( meshVertices[vertexIndex].first );
        } else {
            clusterIndices[rep] = outClusters.size();
            outClusters.push_back( vector<mesh_vertex>() );
            outClusters.back().push_back( meshVertices[vertexIndex].first );
        }
    }
}

polymesh3_ptr geometry::weld_boundary_edges( const vector<mesh_interface_ptr>& meshes, float tolerance ) {
    logging::null_progress_logger progress;

    const size_t meshCount = meshes.size();
    vector<const mesh_interface*> interfaceRefs( meshCount );
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        interfaceRefs[meshIndex] = meshes[meshIndex].get();
    }

    return weld_boundary_edges( interfaceRefs, tolerance, progress );
}

polymesh3_ptr geometry::weld_boundary_edges( const vector<const mesh_interface*>& meshes, float tolerance ) {
    logging::null_progress_logger progress;
    return weld_boundary_edges( meshes, tolerance, progress );
}

polymesh3_ptr geometry::weld_boundary_edges( const vector<mesh_interface_ptr>& meshes, float tolerance,
                                             logging::progress_logger& progress ) {
    const size_t meshCount = meshes.size();
    vector<const mesh_interface*> interfaceRefs( meshCount );
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        interfaceRefs[meshIndex] = meshes[meshIndex].get();
    }

    return weld_boundary_edges( interfaceRefs, tolerance, progress );
}

polymesh3_ptr geometry::weld_boundary_edges( const vector<const mesh_interface*>& meshes, float tolerance,
                                             logging::progress_logger& progress ) {
    logging::progress_logger_subinterval_tracker progressInterval( progress, 0, 100 );

    // Collect boundary edges for each mesh
    progress.set_title( _T("Finding Boundary Edges") );
    progressInterval.reset( 0.0f, 35.0f );

    const size_t meshCount = meshes.size();
    vector<vector<half_edge>> meshBoundaries( meshCount );
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        if( !meshes[meshIndex] ) {
            continue;
        }

        get_boundary_edges( meshes[meshIndex], meshBoundaries[meshIndex] );

        progress.update_progress( meshIndex + 1, meshCount );
        progress.check_for_abort();
    }
    progress.check_for_abort();

    // Find boundary edges which complement one another
    progress.set_title( _T("Computing Boundary Edge Pairs") );
    progressInterval.reset( 35.0f, 50.0f );

    vector<mesh_half_edge_pair> edgePairs;
    find_complement_edges( meshes, meshBoundaries, tolerance, edgePairs );
    progress.check_for_abort();

    // Find boundary edges which complement one another
    progress.set_title( _T("Clustering Boundary Vertices") );
    progressInterval.reset( 50.0f, 60.0f );

    if( edgePairs.size() > 0 ) {
        vector<vector<mesh_vertex>> clusters;
        compute_clusters( edgePairs, meshBoundaries, clusters );
        progress.check_for_abort();

        // Create a welding policy
        vertex_cluster_arbiter arbiter( clusters, meshCount );
        progress.check_for_abort();

        // Weld the meshes
        progress.set_title( _T("Welding Meshes") );
        progressInterval.reset( 60.0f, 100.0f );

        return weld_meshes( meshes, arbiter, progress );
    } else {
        progress.set_title( _T("Combining Meshes") );
        progressInterval.reset( 60.0f, 100.0f );

        return combine( meshes, progress );
    }
}

polymesh3_ptr geometry::combine( const vector<mesh_interface_ptr>& meshes ) {
    logging::null_progress_logger progress;

    const size_t meshCount = meshes.size();
    vector<const mesh_interface*> interfaceRefs( meshCount );
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        interfaceRefs[meshIndex] = meshes[meshIndex].get();
    }

    return combine( interfaceRefs, progress );
}

polymesh3_ptr geometry::combine( const vector<const mesh_interface*>& meshes ) {
    logging::null_progress_logger progress;
    return combine( meshes, progress );
}

polymesh3_ptr geometry::combine( const vector<mesh_interface_ptr>& meshes, logging::progress_logger& progress ) {
    const size_t meshCount = meshes.size();
    vector<const mesh_interface*> interfaceRefs( meshCount );
    for( size_t meshIndex = 0; meshIndex < meshCount; ++meshIndex ) {
        interfaceRefs[meshIndex] = meshes[meshIndex].get();
    }

    return combine( interfaceRefs, progress );
}

polymesh3_ptr geometry::combine( const vector<const mesh_interface*>& meshes, logging::progress_logger& /*progress*/ ) {
    polymesh3_ptr result = combine_geometry( meshes );

    combine_channels( result, meshes );

    return result;
}
