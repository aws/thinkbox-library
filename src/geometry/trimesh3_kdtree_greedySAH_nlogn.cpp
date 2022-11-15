// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_kdtree_node.hpp>
#include <frantic/graphics/boundbox3f.hpp>

using namespace frantic::graphics;

namespace frantic {
namespace geometry {

namespace {
enum kdtreeEvent_type { END = 0, PARALLEL, BEGIN };

const int COST_TRAVERSAL = 5;
const int COST_INTERSECT = 1;
} // namespace

namespace detail {
struct kdtreeEvent {
    kdtreeEvent_type type;
    float location;
    int index;
    int axis;

    kdtreeEvent() {}
    kdtreeEvent( int i, kdtreeEvent_type t, int a, float l )
        : type( t )
        , location( l )
        , index( i )
        , axis( a ) {}
    bool operator<( const kdtreeEvent& rhs ) const {
        return location < rhs.location ||
               ( location == rhs.location && ( axis < rhs.axis || ( axis == rhs.axis && type < rhs.type ) ) );
    }
};
} // namespace detail

inline float SAH_cost( float voxelSA, float leftSA, float rightSA, int nLeft, int nRight ) {
    return ( ( nLeft == 0 || nRight == 0 ) ? 0.8f : 1.f ) *
           ( COST_TRAVERSAL + COST_INTERSECT * ( leftSA / voxelSA * nLeft + rightSA / voxelSA * nRight ) );
}

inline void extract_events( std::vector<detail::kdtreeEvent>& e, int tri, const boundbox3f& v ) {
    for( int i = 0; i < 3; ++i ) {
        if( v.size( i ) > 0.00001f ) {
            e.push_back( detail::kdtreeEvent( tri, BEGIN, i, v.minimum()[i] ) );
            e.push_back( detail::kdtreeEvent( tri, END, i, v.maximum()[i] ) );
        } else
            e.push_back( detail::kdtreeEvent( tri, PARALLEL, i, v.minimum()[i] ) );
    }
}

void build_kdtree_greedy_SAH_nlogn( trimesh3_kdtree_node& node, const trimesh3& mesh, const boundbox3f& bounds,
                                    std::vector<detail::kdtreeEvent>& events, std::vector<int>& indices,
                                    std::vector<int>& triFlags ) {
    if( indices.size() <= 8 ) {
        node.initialize( indices );
        return;
    }

    float voxelSA = bounds.get_surface_area();

    int bestAxis = -1;
    float bestCost = 0.9f * COST_INTERSECT * indices.size();
    float bestSplit = std::numeric_limits<float>::max();

    int nLeft[3], nRight[3];
    nLeft[0] = nLeft[1] = nLeft[2] = 0;
    nRight[0] = nRight[1] = nRight[2] = static_cast<int>( indices.size() );

    for( std::size_t i = 0; i < events.size(); /*do nothing*/ ) {
        int axis = events[i].axis;
        float split = events[i].location;

        int counters[] = { 0, 0, 0 };
        for( ; i < events.size() && events[i].axis == axis && events[i].location == split; ++i )
            ++counters[events[i].type];

        nRight[axis] -= ( counters[END] + counters[PARALLEL] );

        if( split > bounds.minimum()[axis] && split < bounds.maximum()[axis] ) {
            boundbox3f left( bounds ), right( bounds );
            left.maximum()[axis] = right.minimum()[axis] = split;

            float cost = SAH_cost( voxelSA, left.get_surface_area(), right.get_surface_area(),
                                   nLeft[axis] + counters[PARALLEL], nRight[axis] + counters[PARALLEL] );
            if( cost < bestCost ) {
                bestCost = cost;
                bestSplit = split;
                bestAxis = axis;
            }
        }

        nLeft[axis] += ( counters[BEGIN] + counters[PARALLEL] );
    }

    if( bestAxis != -1 ) {
        // We need to classify triangles as being on either side of the split, or crossing
        for( std::vector<int>::const_iterator index = indices.begin(), end = indices.end(); index != end; ++index )
            triFlags[*index] = 0;

        for( std::vector<detail::kdtreeEvent>::const_iterator it = events.begin(), end = events.end(); it != end;
             ++it ) {
            if( it->axis != bestAxis )
                continue;

            if( it->location < bestSplit ) {
                if( it->type != BEGIN )
                    triFlags[it->index] = -1;
            } else if( it->location > bestSplit ) {
                if( it->type != END )
                    triFlags[it->index] = 1;
            } else {
                if( it->type == END )
                    triFlags[it->index] = -1;
                else if( it->type == BEGIN )
                    triFlags[it->index] = 1;
            }
        }

        std::vector<detail::kdtreeEvent> leftEvents[3], rightEvents[3];

        // Most events should be split into one side xor the other
        leftEvents[0].reserve( events.size() / 2 );
        rightEvents[0].reserve( events.size() / 2 );

        for( std::vector<detail::kdtreeEvent>::const_iterator it = events.begin(), end = events.end(); it != end;
             ++it ) {
            if( triFlags[it->index] < 0 )
                leftEvents[0].push_back( *it );
            else if( triFlags[it->index] > 0 )
                rightEvents[0].push_back( *it );
        }

        frantic::clear_with_swap( events );

        std::vector<int> leftIndices, rightIndices;
        leftIndices.reserve( indices.size() / 2 );
        rightIndices.reserve( indices.size() / 2 );

        boundbox3f leftBounds( bounds ), rightBounds( bounds );
        leftBounds.maximum()[bestAxis] = rightBounds.minimum()[bestAxis] = bestSplit;

        for( std::vector<int>::const_iterator it = indices.begin(), end = indices.end(); it != end; ++it ) {
            if( triFlags[*it] < 0 )
                leftIndices.push_back( *it );
            else if( triFlags[*it] > 0 )
                rightIndices.push_back( *it );
            else {
                const vector3& f = mesh.get_face( *it );

                boundbox3f clippedBounds;
                if( leftBounds.intersect_with_triangle( mesh.get_vertex( f.x ), mesh.get_vertex( f.y ),
                                                        mesh.get_vertex( f.z ), clippedBounds ) ) {
                    leftIndices.push_back( *it );
                    extract_events( leftEvents[1], *it, clippedBounds );
                }

                if( rightBounds.intersect_with_triangle( mesh.get_vertex( f.x ), mesh.get_vertex( f.y ),
                                                         mesh.get_vertex( f.z ), clippedBounds ) ) {
                    rightIndices.push_back( *it );
                    extract_events( rightEvents[1], *it, clippedBounds );
                }
            } // triFlags[*it] == 0
        }

        frantic::clear_with_swap( indices );

        node.initialize( new trimesh3_kdtree_node[2], bestAxis, bestSplit );

        leftEvents[2].resize( leftEvents[0].size() + leftEvents[1].size() );
        std::sort( leftEvents[1].begin(), leftEvents[1].end() );
        std::merge( leftEvents[0].begin(), leftEvents[0].end(), leftEvents[1].begin(), leftEvents[1].end(),
                    leftEvents[2].begin() );

        frantic::clear_with_swap( leftEvents[0] );
        frantic::clear_with_swap( leftEvents[1] );

        rightEvents[2].resize( rightEvents[0].size() + rightEvents[1].size() );
        std::sort( rightEvents[1].begin(), rightEvents[1].end() );
        std::merge( rightEvents[0].begin(), rightEvents[0].end(), rightEvents[1].begin(), rightEvents[1].end(),
                    rightEvents[2].begin() );

        frantic::clear_with_swap( rightEvents[0] );
        frantic::clear_with_swap( rightEvents[1] );

        build_kdtree_greedy_SAH_nlogn( *node.left_child(), mesh, leftBounds, leftEvents[2], leftIndices, triFlags );
        build_kdtree_greedy_SAH_nlogn( *node.right_child(), mesh, rightBounds, rightEvents[2], rightIndices, triFlags );
    } else
        node.initialize( indices );
}

void build_kdtree_greedy_SAH_nlogn( trimesh3_kdtree_node& node, const trimesh3& mesh, const boundbox3f& bounds ) {
    std::vector<int> indices, triFlags( mesh.face_count(), 0 );
    std::vector<detail::kdtreeEvent> events;

    indices.reserve( mesh.face_count() );
    events.reserve( 6 * mesh.face_count() );
    for( std::size_t i = 0; i < mesh.face_count(); ++i ) {
        const vector3& f = mesh.get_face( i );

        boundbox3f clippedBounds;
        if( bounds.intersect_with_triangle( mesh.get_vertex( f.x ), mesh.get_vertex( f.y ), mesh.get_vertex( f.z ),
                                            clippedBounds ) ) {
            indices.push_back( static_cast<int>( i ) );
            extract_events( events, static_cast<int>( i ), clippedBounds );
        }
    }

    std::sort( events.begin(), events.end() );
    build_kdtree_greedy_SAH_nlogn( node, mesh, bounds, events, indices, triFlags );
}
} // namespace geometry
} // namespace frantic
