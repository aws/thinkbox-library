// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file
 * @author Stephen Kiazyk
 *
 * A place for methods for removing bad parts of meshes.
 */

#pragma once

#include <frantic/geometry/trimesh3.hpp>

#include <frantic/strings/tstring.hpp>

#include <frantic/logging/progress_logger.hpp>

#include <boost/foreach.hpp>

namespace frantic {
namespace geometry {

template <class FaceFunctor>
bool remove_faces( frantic::geometry::trimesh3& mesh, FaceFunctor aliveCheck,
                   frantic::logging::progress_logger* logger = NULL ) {
    frantic::logging::null_progress_logger nullLogger;

    if( !logger ) {
        logger = &nullLogger;
    }

    logger->update_progress( 0.0f );

    std::vector<frantic::tstring> vertexChannelNames;
    mesh.get_vertex_channel_names( vertexChannelNames );

    std::vector<trimesh3_vertex_channel_general_accessor> customVertexChannelAccessors;
    customVertexChannelAccessors.reserve( vertexChannelNames.size() );
    BOOST_FOREACH( const frantic::tstring& channelName, vertexChannelNames ) {
        trimesh3_vertex_channel_general_accessor acc = mesh.get_vertex_channel_general_accessor( channelName );
        if( acc.has_custom_faces() ) {
            customVertexChannelAccessors.push_back( acc );
        }
    }

    std::vector<frantic::tstring> faceChannelNames;
    mesh.get_face_channel_names( faceChannelNames );

    std::vector<trimesh3_face_channel_general_accessor> faceChannelAccessors;
    faceChannelAccessors.reserve( faceChannelNames.size() );
    BOOST_FOREACH( const frantic::tstring& channelName, faceChannelNames ) {
        faceChannelAccessors.push_back( mesh.get_face_channel_general_accessor( channelName ) );
    }

    logger->update_progress( 5.0f );

    size_t insertionPoint = 0;

    bool anyMods = false;

    for( size_t faceId = 0; faceId < mesh.face_count(); ++faceId ) {
        if( !aliveCheck( faceId ) ) {
            anyMods = true;
        } else {
            if( insertionPoint != faceId ) {
                mesh.get_face( insertionPoint ) = mesh.get_face( faceId );

                BOOST_FOREACH( trimesh3_face_channel_general_accessor& acc, faceChannelAccessors ) {
                    memcpy( acc.data( insertionPoint ), acc.data( faceId ), acc.primitive_size() );
                }

                BOOST_FOREACH( trimesh3_vertex_channel_general_accessor& acc, customVertexChannelAccessors ) {
                    acc.face( insertionPoint ) = acc.face( faceId );
                }
            }

            ++insertionPoint;
        }
    }

    logger->update_progress( 70.0f );

    // shift down the sizes of the face arrays
    mesh.faces_ref().resize( insertionPoint );

    logger->update_progress( 80.0f );

    BOOST_FOREACH( trimesh3_face_channel_general_accessor& acc, faceChannelAccessors ) {
        acc.set_face_count( insertionPoint );
    }

    logger->update_progress( 90.0f );

    BOOST_FOREACH( trimesh3_vertex_channel_general_accessor& acc, customVertexChannelAccessors ) {
        acc.set_face_count( insertionPoint );
    }

    logger->update_progress( 100.0f );

    return anyMods;
}

/**
 * Removes any vertices unreferenced by the mesh.
 *
 * @param mesh the mesh to be modified
 * @return true if the mesh was modified by this method
 */
bool remove_dead_vertices( frantic::geometry::trimesh3& mesh, frantic::logging::progress_logger* logger = NULL );

/**
 * Removes any faces in the mesh which share the same three vertices,
 * though not neccessarily in the same order. These are almost always
 * unwanted, and furthermore causes adjacency information generation
 * to break.
 *
 * @param mesh the mesh to be modified
 * @return true if the mesh was modified by this method
 */
bool remove_doubled_faces( frantic::geometry::trimesh3& mesh, frantic::logging::progress_logger* logger = NULL );

/**
 * Removes any face in the mesh which has more than one of the
 * same vertex index on it
 *
 * @param mesh the mesh to be modified
 * @return true if the mesh was modified by this method
 */
bool remove_duplicate_face_indices( frantic::geometry::trimesh3& mesh,
                                    frantic::logging::progress_logger* logger = NULL );

/**
 * Looks for overlapping geometry (specifically faces with area zero) and
 * attempts to repair them.
 * Note that existing face and vertex references might be invalided by this
 * call if any modifications are made.
 *
 * @remark This method will fail if there are any vertex channels with custom face storage,
 *         or vertex/face channels of different size than the position/geometry channels
 * @remark This method will only detect co-incident vertices if they share at least one face
 *         a more involved check involving a kd-tree nearest neighbour query would be required
 *         to detect and join arbitrary vertices.
 *
 * @param epsilon threshold for detecting two co-incident vertices
 * @return true if any modifications were made to the mesh
 * @throws std::runtime_error if any custom storage channels were detected
 */
bool remove_coincident_geometry( frantic::geometry::trimesh3& mesh, double epsilon = 1.0e-10f,
                                 frantic::logging::progress_logger* logger = NULL );

struct edge_with_face {
    std::pair<boost::uint32_t, boost::uint32_t> m_edge;
    boost::uint32_t m_face;

    edge_with_face() {}

    edge_with_face( boost::uint32_t edge0, boost::uint32_t edge1, boost::uint32_t face )
        : m_edge( edge0, edge1 )
        , m_face( face ) {}

    bool operator<( const edge_with_face& other ) const {
        if( m_edge < other.m_edge ) {
            return true;
        } else if( m_edge == other.m_edge ) {
            return m_face < other.m_face;
        }
        return false;
    }

    bool operator==( const edge_with_face& other ) const { return m_edge == other.m_edge && m_face == other.m_face; }
};

/**
 * Returns a list of all halfedges with faces in the mesh
 *
 * @param mesh the mesh to read
 * @param outEdges the resulting list
 * @param sorted (optional) if true, outEdges will be sorted by edge indices
 */
void collect_edges_with_faces( const frantic::geometry::trimesh3& mesh, std::vector<edge_with_face>& outEdges,
                               bool sorted = false );

/**
 * Detects and removes faces such that the resulting mesh has no duplicated halfedges.
 *
 * @param mesh The mesh to fix
 * @param logger (optional) progress logger
 * @return true if the mesh was changed in any way
 */
bool remove_duplicate_halfedge_faces( frantic::geometry::trimesh3& mesh, frantic::logging::progress_logger& logger,
                                      bool fullyDegenerateOnly = true );
bool remove_duplicate_halfedge_faces( frantic::geometry::trimesh3& mesh, bool fullyDegenerateOnly = true );

/**
 * Identifies vertices which are causing zero-area faces.
 * Currently handles:
 * - epsilon-length edges
 * - T-vertices (vertices which are incident to a non-endpoint part of the edge)
 * @todo It may make sense to split this out into multiple degeneracy detection methods
 *
 * @param mesh the mesh to test
 * @param outSlivers the list of sliver vertices located
 * @param epsilon (optional) the distance tolerance for identifying features as 'small'
 */
void locate_sliver_vertices( const frantic::geometry::trimesh3& mesh, std::vector<size_t>& outSlivers,
                             double epsilon = 1.0e-5 );

/**
 * Runs localized relaxation on vertices which are causing zero-area faces.
 * The vertices are those identified by `locate_sliver_vertices`
 *
 * @param mesh the mesh to operate on
 * @param iterations the number of iterations
 * @param scale the scaling factor to use in relaxation
 * @param epsilon (optional) the distance tolerance for identifying features as 'small'
 */
bool relax_sliver_vertices( frantic::geometry::trimesh3& mesh, size_t iterations, float scale, float epsilon = 1.0e-5 );

} // namespace geometry
} // namespace frantic
