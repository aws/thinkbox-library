// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/container/small_vector.hpp>
#include <boost/foreach.hpp>
#if defined( FRANTIC_TBB_AVAILABLE )
#include <tbb/parallel_sort.h>
#else
#include <algorithm>
#endif

#include <frantic/geometry/edge_topology.hpp>

using std::deque;
using std::pair;
using std::size_t;
using std::vector;

using namespace frantic;
using geometry::topology::half_edge;
using geometry::topology::mesh_half_edge;
using graphics::vector3f;

namespace {

/**
 * Representation of a directed connection in an adjacency list with a flag marking the connection as boundary or
 * non-boundary.
 */
struct connection {
  private:
    size_t m_head;

  public:
    connection() {}
    connection( size_t head, bool boundary )
        : m_head( head ) {
        if( boundary ) {
            set_boundary();
        } else {
            set_non_boundary();
        }
    }

    size_t get_head() { return m_head & ( std::numeric_limits<size_t>::max() >> 1 ); }

    void set_head( size_t head ) {
        m_head = ( m_head & ~( std::numeric_limits<size_t>::max() >> 1 ) ) |
                 ( head & ( std::numeric_limits<size_t>::max() >> 1 ) );
    }

    bool is_boundary() { return ( m_head & ~( std::numeric_limits<size_t>::max() >> 1 ) ) != 0; }

    void set_boundary() { m_head |= ~( std::numeric_limits<size_t>::max() >> 1 ); }

    void set_non_boundary() { m_head &= ( std::numeric_limits<size_t>::max() >> 1 ); }
};
} // anonymous namespace

void geometry::get_boundary_edges( const mesh_interface* mesh, vector<half_edge>& outEdges ) {
    using boost::container::small_vector;

    if( !mesh ) {
        throw std::runtime_error( "get_boundary_edges Error: mesh is NULL" );
    }

    const size_t faceCount = mesh->get_num_faces();
    const size_t vertexCount = mesh->get_num_verts();
    size_t boundaryCount = 0;

    // An adjacency list to record edges from each vertex
    // Note: All connections are one-sided (if A->B is in the list, B->A will not be)
    vector<small_vector<connection, 6>> connections( vertexCount );

    vector<size_t> faceVerts;
    for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        const size_t fvCount = mesh->get_num_face_verts( faceIndex );

        faceVerts.resize( fvCount );
        mesh->get_face_vert_indices( faceIndex, &faceVerts[0] );

        // Record each edge on the face, in proper winding order
        for( size_t fvIndex = 0; fvIndex < fvCount; ++fvIndex ) {
            const size_t tail = faceVerts[fvIndex];
            const size_t head = ( fvIndex != fvCount - 1 ) ? faceVerts[fvIndex + 1] : faceVerts[0];

            small_vector<connection, 6>::iterator edge, edgeEnd;
            for( edge = connections[tail].begin(), edgeEnd = connections[tail].end(); edge < edgeEnd; ++edge ) {
                if( edge->get_head() == head ) {
                    // This half edge has been seen before, so it is not a boundary
                    if( edge->is_boundary() ) {
                        --boundaryCount;
                        edge->set_non_boundary();
                    }
                    break;
                }
            }
            if( edge == edgeEnd ) {
                for( edge = connections[head].begin(), edgeEnd = connections[head].end(); edge < edgeEnd; ++edge ) {
                    if( edge->get_head() == tail ) {
                        // The opposite half edge has been seen before, so it is not a boundary
                        if( edge->is_boundary() ) {
                            --boundaryCount;
                            edge->set_non_boundary();
                        }
                        break;
                    }
                }
                if( edge == edgeEnd ) {
                    // The edge has not been seen before, record it and assume it is a boundary for now
                    connections[tail].push_back( connection( head, true ) );
                    ++boundaryCount;
                }
            }
        }
    }

    // Output only the boundary edges, with their original winding order
    outEdges.clear();
    outEdges.reserve( boundaryCount );
    for( size_t tail = 0, headEnd = connections.size(); tail < headEnd; ++tail ) {
        for( small_vector<connection, 6>::iterator edge = connections[tail].begin(), edgeEnd = connections[tail].end();
             edge < edgeEnd; ++edge ) {
            if( edge->is_boundary() ) {
                outEdges.push_back( half_edge( edge->get_head(), tail ) );
            }
        }
    }
}

namespace {

struct projection_bound {
    mesh_half_edge m_meshEdge;
    float m_value;

    projection_bound() {}
    projection_bound( mesh_half_edge meshEdge, float value )
        : m_meshEdge( meshEdge )
        , m_value( value ) {}

    bool operator<( const projection_bound& other ) const { return m_value < other.m_value; }

    bool operator==( const projection_bound& other ) const { return m_value == other.m_value; }
};

struct bound_data {
    size_t m_bound;
    vector3f m_head;
    vector3f m_tail;

    bound_data() {}
    bound_data( size_t bound, vector3f head, vector3f tail )
        : m_bound( bound )
        , m_head( head )
        , m_tail( tail ) {}
};
} // anonymous namespace

void geometry::find_complement_edges( const vector<const mesh_interface*>& meshes,
                                      const vector<vector<half_edge>>& meshEdges, float tolerance,
                                      vector<pair<mesh_half_edge, mesh_half_edge>>& outComplements,
                                      vector3f scanAxis ) {
    // Edges must be within tolerance on both ends, so we perform a broad-phase collision check with a sort-and-sweep
    // algorithm to find edges whose scan-axis projections are within tolerance on their lower bound.

    if( meshes.size() != meshEdges.size() ) {
        throw std::runtime_error( "find_complement_edges Error: mesh and mesh-edge lists are not the same size" );
    }
    vector<size_t> vertexCounts( meshes.size() );
    for( size_t meshIndex = 0, meshCount = meshes.size(); meshIndex < meshCount; ++meshIndex ) {
        if( !meshes[meshIndex] ) {
            continue;
        }
        vertexCounts[meshIndex] = meshes[meshIndex]->get_num_verts();
    }

    // Collect lower bounds for all edges
    float toleranceSquared = tolerance * tolerance;
    scanAxis.normalize();

    vector<projection_bound> bounds;
    for( size_t meshIndex = 0, meshCount = meshes.size(); meshIndex < meshCount; ++meshIndex ) {
        if( !meshes[meshIndex] ) {
            continue;
        }

        for( size_t edgeIndex = 0, edgeCount = meshEdges[meshIndex].size(); edgeIndex < edgeCount; ++edgeIndex ) {
            const half_edge& edge = meshEdges[meshIndex][edgeIndex];
            if( edge.m_head >= vertexCounts[meshIndex] || edge.m_tail >= vertexCounts[meshIndex] ) {
                throw std::runtime_error( "find_complement_edges Error: mesh edge referenced out-of-range vertex" );
            }

            vector3f head = meshes[meshIndex]->get_vert( edge.m_head );
            vector3f tail = meshes[meshIndex]->get_vert( edge.m_tail );
            float minimumBound = std::min( vector3f::dot( head, scanAxis ), vector3f::dot( tail, scanAxis ) );

            bounds.push_back( projection_bound( mesh_half_edge( meshIndex, edgeIndex ), minimumBound ) );
        }
    }

    // Use Sort-and-Sweep to find collisions
    outComplements.clear();
#if defined( FRANTIC_TBB_AVAILABLE )
    tbb::parallel_sort( bounds.begin(), bounds.end() );
#else
    std::sort( bounds.begin(), bounds.end() );
#endif

    deque<bound_data> workingSet;
    for( size_t boundIndex = 0, boundCount = bounds.size(); boundIndex < boundCount; ++boundIndex ) {
        const projection_bound& nextBound = bounds[boundIndex];
        const half_edge& nextEdge = meshEdges[nextBound.m_meshEdge.m_mesh][nextBound.m_meshEdge.m_edge];

        while( workingSet.size() > 0 && bounds[workingSet.front().m_bound].m_value < nextBound.m_value - tolerance ) {
            workingSet.pop_front();
        }

        bound_data nextBoundData( boundIndex, meshes[nextBound.m_meshEdge.m_mesh]->get_vert( nextEdge.m_head ),
                                  meshes[nextBound.m_meshEdge.m_mesh]->get_vert( nextEdge.m_tail ) );

        // Check edge criteria against all other edges in the working set
        bool matchFound = false;
        for( deque<bound_data>::reverse_iterator boundData = workingSet.rbegin(), boundDataEnd = workingSet.rend();
             boundData != boundDataEnd; ++boundData ) {
            if( bounds[boundData->m_bound].m_meshEdge.m_mesh == nextBound.m_meshEdge.m_mesh ) {
                continue;
            }

            if( vector3f::distance_squared( boundData->m_head, nextBoundData.m_tail ) <= toleranceSquared &&
                vector3f::distance_squared( boundData->m_tail, nextBoundData.m_head ) <= toleranceSquared ) {
                outComplements.push_back( pair<mesh_half_edge, mesh_half_edge>( bounds[boundData->m_bound].m_meshEdge,
                                                                                nextBound.m_meshEdge ) );
                workingSet.erase( ( boundData + 1 ).base() );
                matchFound = true;
                break;
            }
        }

        if( !matchFound ) {
            workingSet.push_back( nextBoundData );
        }
    }
}
