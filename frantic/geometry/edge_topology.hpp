// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace geometry {

namespace topology {

/**
 * Given a list of vertices for a mesh, this struct denotes an edge from one vertex to another by index.
 */
struct half_edge {
    std::size_t m_head;
    std::size_t m_tail;

    half_edge() {}
    half_edge( std::size_t head, std::size_t tail )
        : m_head( head )
        , m_tail( tail ) {}

    bool operator<( const half_edge& other ) const {
        return ( m_head < other.m_head || ( m_head == other.m_head && m_tail < other.m_tail ) );
    }

    bool operator==( const half_edge& other ) const { return ( m_head == other.m_head && m_tail == other.m_tail ); }
};

/**
 * Given a list of meshes, and a list of edges for each mesh, this struct is used to index a specific edge from one of
 * the meshes. The internal values index the list of meshes and the list of edges for the selected mesh.
 */
struct mesh_half_edge {
    std::size_t m_mesh;
    std::size_t m_edge;

    mesh_half_edge() {}
    mesh_half_edge( std::size_t mesh, std::size_t edge )
        : m_mesh( mesh )
        , m_edge( edge ) {}

    bool operator<( const mesh_half_edge& other ) const {
        return ( m_mesh < other.m_mesh || ( m_mesh == other.m_mesh && m_edge < other.m_edge ) );
    }

    bool operator==( const mesh_half_edge& other ) const {
        return ( m_mesh == other.m_mesh && m_edge == other.m_edge );
    }
};

/**
 * Given a list of meshes, this struct is used to index a specific vertex from one of the meshes. The internal values
 * index the list of meshes and the list of vertices for the selected mesh.
 */
struct mesh_vertex {
    std::size_t m_mesh;
    std::size_t m_index;

    mesh_vertex() {}
    mesh_vertex( std::size_t mesh, std::size_t index )
        : m_mesh( mesh )
        , m_index( index ) {}

    bool operator<( const mesh_vertex& other ) const {
        return ( m_mesh < other.m_mesh || ( m_mesh == other.m_mesh && m_index < other.m_index ) );
    }

    bool operator==( const mesh_vertex& other ) const { return ( m_mesh == other.m_mesh && m_index == other.m_index ); }
};
} // namespace topology

/**
 * Return a list of all boundary edges in the mesh. Boundary edges are edges with exactly one incident face.
 *
 * The winding order of the output edges is the same as the winding order of its single incident face.
 *
 * @param mesh The mesh to check.
 * @param[out] outEdges The list of boundary edges.
 */
void get_boundary_edges( const mesh_interface* mesh, std::vector<topology::half_edge>& outEdges );

/**
 * Return a list of complementary edges between a set of meshes.
 *
 * Two edges are considered complementary if they meet the following criteria.
 *  - They are from different input meshes
 *  - The head of one edge is within a given tolerance of the tail of the other edge, and vice versa
 *  - Neither edge is already being matched with another edge
 * The scan axis is used to determine broad-phase collisions by projecting each edge onto the axis. When providing a
 * custom axis, it should be chosen so that edges project onto the axis in an evenly distributed manner.
 *
 * @param meshes A list of meshes.
 * @param meshEdges A list of edges for each mesh.
 * @param tolerance The maximum distance between edge endpoints.
 * @param[out] outComplements A list of complementary edge pairs.
 * @param scanAxis The axis to project the edge bounding boxes onto.
 */
void find_complement_edges( const std::vector<const mesh_interface*>& meshes,
                            const std::vector<std::vector<topology::half_edge>>& meshEdges, float tolerance,
                            std::vector<std::pair<topology::mesh_half_edge, topology::mesh_half_edge>>& outComplements,
                            graphics::vector3f scanAxis = graphics::vector3f( 1, 1, 1 ) );
} // namespace geometry
} // namespace frantic
