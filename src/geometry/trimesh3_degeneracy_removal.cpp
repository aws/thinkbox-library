// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/trimesh3_degeneracy_removal.hpp>

#include <frantic/logging/progress_logger.hpp>

#include <frantic/geometry/trimesh3.hpp>

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/dcel_iterators.hpp>

#include <frantic/geometry/relaxation.hpp>

#include <frantic/misc/functor.hpp>

#include <boost/foreach.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <tbb/atomic.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>

#include <algorithm>
#include <utility>
#include <vector>

namespace frantic {
namespace geometry {

namespace {

typedef std::pair<frantic::graphics::vector3, size_t> face_ref;

struct face_ref_comparator {
    bool operator()( const face_ref& lhs, const face_ref& rhs ) const {
        for( size_t i = 0; i < 3; ++i ) {
            if( lhs.first[i] < rhs.first[i] ) {
                return true;
            } else if( lhs.first[i] > rhs.first[i] ) {
                return false;
            }
        }
        return false;
    }
};

} // namespace

// TODO:
// - Add parallelism (everything except the last part can be easily parallelized)
bool remove_doubled_faces( frantic::geometry::trimesh3& mesh, frantic::logging::progress_logger* logger ) {

    frantic::logging::null_progress_logger nullLogger;

    if( !logger ) {
        logger = &nullLogger;
    }

    logger->update_progress( 0.0f );

    std::vector<frantic::tstring> vertexChannelNames;
    mesh.get_vertex_channel_names( vertexChannelNames );

    std::vector<frantic::tstring> faceChannelNames;
    mesh.get_face_channel_names( faceChannelNames );

    std::vector<face_ref> facesRef( mesh.face_count() );

    for( size_t i = 0; i < mesh.face_count(); ++i ) {
        facesRef[i] = std::make_pair( mesh.get_face( i ).to_sorted(), i );
    }

    logger->update_progress( 10.0f );

    std::sort( facesRef.begin(), facesRef.end(), face_ref_comparator() );

    logger->update_progress( 50.0f );

    std::vector<char> living( mesh.face_count(), true );

    bool anyMods = false;

    for( size_t i = 0; i < facesRef.size() - 1; ++i ) {
        if( facesRef[i].first == facesRef[i + 1].first ) {
            living[facesRef[i].second] = false;
            living[facesRef[i + 1].second] = false;
            anyMods = true;
        }
    }

    logger->update_progress( 60.0f );

    {
        std::vector<face_ref> swapAndKill;
        facesRef.swap( swapAndKill );
    }

    logger->update_progress( 70.0f );

    {
        frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, 70.0f, 100.0f );
        remove_faces( mesh, frantic::make_bool_array_functor( living ), logger );
    }

    return anyMods;
}

// TODO:
// - Add parallelism (though might not gain all that much in this case)
bool remove_dead_vertices( frantic::geometry::trimesh3& mesh, frantic::logging::progress_logger* logger ) {
    frantic::logging::null_progress_logger nullLogger;

    if( !logger ) {
        logger = &nullLogger;
    }

    logger->update_progress( 0.0f );

    size_t totalUpdates = mesh.count_named_vertex_channels_with_custom_faces() + 1;
    float updateStride = 100.0f / totalUpdates;
    bool anyMods = false;

    std::vector<frantic::tstring> vertexChannelNames;
    mesh.get_vertex_channel_names( vertexChannelNames );

    typedef trimesh3_vertex_channel_general_accessor vertex_channel_t;
    std::vector<vertex_channel_t> vertexChannels;
    BOOST_FOREACH( const frantic::tstring& channelName, vertexChannelNames ) {
        vertexChannels.push_back( mesh.get_vertex_channel_general_accessor( channelName ) );
    }

    std::vector<char> living( mesh.vertex_count(), false );
    std::vector<trimesh3::index_t> vertexRemap( mesh.vertex_count() );

    {
        frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, 0.0, updateStride );

        logger->update_progress( 10.0f );

        for( size_t i = 0; i < mesh.face_count(); ++i ) {
            for( size_t j = 0; j < 3; ++j ) {
                living[mesh.get_face( i )[j]] = true;
            }
        }

        logger->update_progress( 30.0f );

        size_t insertionPoint = 0;

        for( size_t i = 0; i < mesh.vertex_count(); ++i ) {
            if( living[i] ) {
                vertexRemap[i] = trimesh3::index_t( insertionPoint );
                mesh.get_vertex( insertionPoint ) = mesh.get_vertex( i );
                BOOST_FOREACH( vertex_channel_t& acc, vertexChannels ) {
                    if( !acc.has_custom_faces() && insertionPoint != i ) {
                        std::memcpy( acc.data( insertionPoint ), acc.data( i ), acc.primitive_size() );
                    }
                }
                ++insertionPoint;
            } else {
                anyMods = true;
            }
        }

        mesh.vertices_ref().resize( insertionPoint );

        BOOST_FOREACH( vertex_channel_t& acc, vertexChannels ) {
            if( !acc.has_custom_faces() ) {
                acc.set_vertex_count( insertionPoint );
            }
        }

        logger->update_progress( 70.0f );

        for( size_t i = 0; i < mesh.face_count(); ++i ) {
            vector3& face = mesh.get_face( i );
            for( size_t j = 0; j < 3; ++j ) {
                if( living[face[j]] ) {
                    face[j] = vertexRemap[face[j]];
                }
            }
        }

        logger->update_progress( 100.0f );
    }

    size_t currentChannel = 1;

    BOOST_FOREACH( vertex_channel_t& acc, vertexChannels ) {
        if( acc.has_custom_faces() ) {
            frantic::logging::progress_logger_subinterval_tracker subinterval( *logger, updateStride * currentChannel,
                                                                               updateStride * ( currentChannel + 1 ) );
            ++currentChannel;

            logger->update_progress( 0.0f );

            living.resize( acc.size() );
            living.assign( acc.size(), false );

            vertexRemap.resize( acc.size() );

            for( size_t i = 0; i < acc.face_count(); ++i ) {
                for( size_t j = 0; j < 3; ++j ) {
                    living[acc.face( i )[j]] = true;
                }
            }

            logger->update_progress( 30.0f );

            size_t insertionPoint = 0;

            for( size_t i = 0, ie = acc.size(); i < ie; ++i ) {
                if( living[i] ) {
                    vertexRemap[i] = trimesh3::index_t( insertionPoint );
                    if( insertionPoint != i ) {
                        std::memcpy( acc.data( insertionPoint ), acc.data( i ), acc.primitive_size() );
                    }
                    ++insertionPoint;
                } else {
                    anyMods = true;
                }
            }

            acc.set_vertex_count( insertionPoint );

            logger->update_progress( 70.0f );

            for( size_t i = 0; i < acc.face_count(); ++i ) {
                vector3& face = acc.face( i );
                for( size_t j = 0; j < 3; ++j ) {
                    if( living[face[j]] ) {
                        face[j] = vertexRemap[face[j]];
                    }
                }
            }

            logger->update_progress( 100.0f );
        }
    }

    return anyMods;
}

namespace {

struct duplicate_face_indices_predicate {
  private:
    duplicate_face_indices_predicate& operator=( const duplicate_face_indices_predicate& ); // prevent warning
    const frantic::geometry::trimesh3& m_mesh;

  public:
    duplicate_face_indices_predicate( const frantic::geometry::trimesh3& mesh )
        : m_mesh( mesh ) {}
    bool operator()( size_t faceId ) const {
        const frantic::graphics::vector3& face( m_mesh.get_face( faceId ) );
        for( size_t faceVertex = 0; faceVertex < 3; ++faceVertex ) {
            if( face[faceVertex] == face[( faceVertex + 1 ) % 3] ) {
                return false;
            }
        }
        return true;
    }
};

} // namespace

// TODO: add parallelism
bool remove_duplicate_face_indices( trimesh3& mesh, frantic::logging::progress_logger* logger ) {
    bool removed = remove_faces( mesh, duplicate_face_indices_predicate( mesh ), logger );
    if( removed ) {
        remove_dead_vertices( mesh, logger );
    }
    return removed;
}

namespace {

void get_face_channel_accessors( trimesh3& mesh,
                                 std::vector<trimesh3_face_channel_general_accessor>& outFaceChannels ) {
    outFaceChannels.clear();

    std::vector<frantic::tstring> channelNames;
    mesh.get_face_channel_names( channelNames );
    BOOST_FOREACH( const frantic::tstring& channelName, channelNames ) {
        outFaceChannels.push_back( mesh.get_face_channel_general_accessor( channelName ) );
    }
}

void get_custom_vertex_channel_accessors(
    trimesh3& mesh, std::vector<trimesh3_vertex_channel_general_accessor>& outCustomVertexChannels ) {
    outCustomVertexChannels.clear();

    std::vector<frantic::tstring> channelNames;
    mesh.get_vertex_channel_names( channelNames );
    BOOST_FOREACH( const frantic::tstring& channelName, channelNames ) {
        trimesh3_vertex_channel_general_accessor acc = mesh.get_vertex_channel_general_accessor( channelName );
        if( acc.has_custom_faces() ) {
            outCustomVertexChannels.push_back( acc );
        }
    }
}

} // anonymous namespace

bool remove_coincident_geometry( trimesh3& mesh, double epsilon, frantic::logging::progress_logger* logger ) {
    using namespace frantic::graphics;
    using namespace frantic::logging;

    frantic::logging::null_progress_logger nullLogger;
    if( !logger ) {
        logger = &nullLogger;
    }

    logger->update_progress( 0.f );

    bool anyMods = false;

    dcel edgeStructure;
    trimesh3_to_dcel( mesh, edgeStructure );

    logger->update_progress( 40.f );

    size_t vertexRemovalCount = 0;
    size_t faceRemovalCount = 0;

    {
        progress_logger_subinterval_tracker subinterval( *logger, 40.f, 50.f );

        const size_t halfedgeUpdateInterval = edgeStructure.halfedge_count() / 30 + 1;

        for( size_t i = 0; i < edgeStructure.halfedge_count(); ++i ) {
            dcel::halfedge_handle handle = edgeStructure.get_halfedge( i );

            if( handle.is_valid() ) {
                size_t incidentFaces = 2;

                if( handle.current_face() == dcel::INVALID_FACE_INDEX ) {
                    --incidentFaces;
                }
                if( handle.opposite_face() == dcel::INVALID_FACE_INDEX ) {
                    --incidentFaces;
                }

                if( incidentFaces > 0 ) {
                    dcel::index_t sourceVertex = handle.source_vertex();
                    dcel::index_t targetVertex = handle.target_vertex();

                    const vector3f& sourceVertexPos = mesh.get_vertex( sourceVertex );
                    const vector3f& targetVertexPos = mesh.get_vertex( targetVertex );

                    double diff = vector3f::distance_squared_double( sourceVertexPos, targetVertexPos );

                    if( diff <= epsilon ) {
                        bool collapse = edgeStructure.try_collapse_edge( handle );
                        if( collapse ) {
                            ++vertexRemovalCount;
                            faceRemovalCount += incidentFaces;
                            anyMods = true;
                        }
                    }
                }
            }

            if( i % halfedgeUpdateInterval == 0 ) {
                logger->update_progress( ( 100.f * i ) / edgeStructure.halfedge_count() );
            }
        }
    }

    logger->update_progress( 50.f );

    if( anyMods ) {
        std::vector<trimesh3_face_channel_general_accessor> faceChannels;
        get_face_channel_accessors( mesh, faceChannels );

        std::vector<trimesh3_vertex_channel_general_accessor> customVertexChannels;
        get_custom_vertex_channel_accessors( mesh, customVertexChannels );

        const size_t faceUpdateInterval = edgeStructure.face_count() / 20 + 1;

        size_t faceWritePosition = 0;

        progress_logger_subinterval_tracker subinterval( *logger, 50.f, 75.f );

        for( size_t i = 0; i < edgeStructure.face_count(); ++i ) {
            if( edgeStructure.has_face_halfedge( i ) ) {
                vector3& face = mesh.get_face( faceWritePosition );
                size_t faceVertex = 0;
                BOOST_FOREACH( dcel::const_halfedge_handle handle,
                               dcel_face_cycle_range( edgeStructure.get_face_halfedge( i ) ) ) {
                    face[faceVertex] = trimesh3::index_t( handle.source_vertex() );
                    ++faceVertex;
                }

                if( faceWritePosition != i ) {
                    BOOST_FOREACH( trimesh3_face_channel_general_accessor& channel, faceChannels ) {
                        memcpy( channel.data( faceWritePosition ), channel.data( i ), channel.primitive_size() );
                    }

                    BOOST_FOREACH( trimesh3_vertex_channel_general_accessor& channel, customVertexChannels ) {
                        channel.face( faceWritePosition ) = channel.face( i );
                    }
                }

                ++faceWritePosition;
            }

            if( i % faceUpdateInterval == 0 ) {
                logger->update_progress( ( 100.f * i ) / edgeStructure.face_count() );
            }
        }

        mesh.set_face_count( faceWritePosition );

        BOOST_FOREACH( trimesh3_face_channel_general_accessor& channel, faceChannels ) {
            channel.set_face_count( faceWritePosition );
        }

        BOOST_FOREACH( trimesh3_vertex_channel_general_accessor& channel, customVertexChannels ) {
            channel.set_face_count( faceWritePosition );
        }

        subinterval.reset( 75.f, 100.f );
        remove_dead_vertices( mesh, logger );
    }

    logger->update_progress( 100.f );

    return anyMods;
}

namespace {

class parallel_collect_edges_with_face {
  private:
    const std::vector<frantic::graphics::vector3>& m_faces;
    std::vector<edge_with_face>& m_edgeSet;

    parallel_collect_edges_with_face& operator=( const parallel_collect_edges_with_face& );

  public:
    parallel_collect_edges_with_face( const std::vector<frantic::graphics::vector3>& faces,
                                      std::vector<edge_with_face>& edgeSet )
        : m_faces( faces )
        , m_edgeSet( edgeSet ) {}

    void operator()( const tbb::blocked_range<size_t>& range ) const {
        for( size_t i = range.begin(); i != range.end(); ++i ) {
            const frantic::graphics::vector3& face = m_faces[i];
            for( size_t j = 0; j < 3; ++j ) {
                boost::uint32_t i0 = face[j];
                boost::uint32_t i1 = face[( j + 1 ) % 3];
                m_edgeSet[i * 3 + j] = edge_with_face( i0, i1, boost::uint32_t( i ) );
            }
        }
    }

    static void invoke( const std::vector<frantic::graphics::vector3>& faces, std::vector<edge_with_face>& edgeSet ) {
        edgeSet.resize( faces.size() * 3 );
        parallel_collect_edges_with_face op( faces, edgeSet );
        tbb::parallel_for( tbb::blocked_range<size_t>( 0, faces.size() ), op );
    }
};

} // namespace

void collect_edges_with_faces( const frantic::geometry::trimesh3& mesh, std::vector<edge_with_face>& outEdges,
                               bool sorted ) {
    using namespace frantic::logging;

    parallel_collect_edges_with_face::invoke( mesh.faces_ref(), outEdges );

    if( sorted ) {
        tbb::parallel_sort( outEdges.begin(), outEdges.end() );
    }
}

namespace {

bool collect_manifold_vertices( const std::vector<edge_with_face>& edges, size_t numVertices,
                                std::vector<char>& isManifoldVertex, frantic::logging::progress_logger& logger ) {
    isManifoldVertex.clear();
    isManifoldVertex.resize( numVertices, char( true ) );

    const size_t progressUpdateInterval = 1 << 25;
    bool anyFound = false;

    for( size_t i = 1; i < edges.size(); ++i ) {
        if( edges[i].m_edge == edges[i - 1].m_edge ) {
            FF_LOG( debug ) << "Duplicate halfedge (" << edges[i].m_edge.first << "," << edges[i].m_edge.second
                            << ") detected on face " << edges[i].m_face << std::endl;
            isManifoldVertex[edges[i].m_edge.first] = char( false );
            isManifoldVertex[edges[i].m_edge.second] = char( false );
            anyFound = true;
        }

        if( i % progressUpdateInterval == 0 ) {
            logger.update_progress( i, edges.size() );
            logger.check_for_abort();
        }
    }

    return anyFound;
}

// kill any face where every incident vertex is non-manifold
bool find_non_fully_degenerate_faces( const frantic::geometry::trimesh3& mesh, const std::vector<edge_with_face>& edges,
                                      std::vector<char>& keepFaces, frantic::logging::progress_logger& logger ) {
    using namespace frantic::logging;

    const size_t faceCount = mesh.face_count();
    keepFaces.clear();
    keepFaces.resize( faceCount, char( true ) );

    size_t progressUpdateInterval = ( faceCount / 100 ) + 1;
    float splitPoint = float( edges.size() ) / ( float( edges.size() + faceCount ) );
    std::vector<char> isManifoldVertex;

    bool anyVerticesFound = false;

    {
        progress_logger_subinterval_tracker subinterval( logger, 0.f, splitPoint );
        anyVerticesFound = collect_manifold_vertices( edges, mesh.vertex_count(), isManifoldVertex, logger );
    }

    bool anyFacesFound = false;

    if( anyVerticesFound ) {
        progress_logger_subinterval_tracker subinterval( logger, splitPoint, 100.f );

        for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
            const frantic::graphics::vector3& face = mesh.get_face( faceIndex );
            bool hasManifoldVertex = false;

            for( int corner = 0; corner < 3; ++corner ) {
                int vertexIndex = face[corner];
                if( isManifoldVertex[vertexIndex] ) {
                    hasManifoldVertex = true;
                    break;
                }
            }

            if( !hasManifoldVertex ) {
                keepFaces[faceIndex] = char( false );
                anyFacesFound = true;
            }

            if( faceIndex % progressUpdateInterval == 0 ) {
                logger.update_progress( faceIndex, faceCount );
                logger.check_for_abort();
            }
        }
    }

    return anyFacesFound;
}

// kill any face that has at least one non-manifold vertex incident to it
bool find_non_partially_degenerate_faces( const frantic::geometry::trimesh3& mesh,
                                          const std::vector<edge_with_face>& edges, std::vector<char>& aliveFaces,
                                          frantic::logging::progress_logger& logger ) {
    bool anyMods = false;

    aliveFaces.clear();
    aliveFaces.resize( mesh.face_count(), char( true ) );

    size_t progressUpdateInterval = ( edges.size() / 100 ) + 1;

    for( size_t i = 1; i < edges.size(); ++i ) {
        if( edges[i].m_edge == edges[i - 1].m_edge ) {
            FF_LOG( debug ) << "Duplicate halfedge (" << edges[i].m_edge.first << "," << edges[i].m_edge.second
                            << ") detected on face " << edges[i].m_face << std::endl;
            aliveFaces[edges[i].m_face] = char( false );
            anyMods = true;
        }

        if( i % progressUpdateInterval == 0 ) {
            logger.update_progress( i, edges.size() );
            logger.check_for_abort();
        }
    }

    return anyMods;
}

} // namespace

bool remove_duplicate_halfedge_faces( frantic::geometry::trimesh3& mesh, frantic::logging::progress_logger& logger,
                                      bool fullyDegenerateOnly ) {
    using namespace frantic::logging;

    logger.update_progress( 0.f );

    bool anyMods = false;

    std::vector<edge_with_face> edges;

    {
        progress_logger_subinterval_tracker subinterval( logger, 0.f, 30.f );
        collect_edges_with_faces( mesh, edges, true );
    }

    std::vector<char> aliveFaces( mesh.face_count(), char( true ) );

    {
        progress_logger_subinterval_tracker subinterval( logger, 30.f, 60.f );

        if( fullyDegenerateOnly ) {
            anyMods = find_non_fully_degenerate_faces( mesh, edges, aliveFaces, logger );
        } else {
            anyMods = find_non_partially_degenerate_faces( mesh, edges, aliveFaces, logger );
        }
    }

    if( anyMods ) {
        {
            progress_logger_subinterval_tracker subinterval( logger, 60.f, 80.f );
            remove_faces( mesh, frantic::make_bool_array_functor( aliveFaces ), &logger );
        }

        {
            progress_logger_subinterval_tracker subinterval( logger, 80.f, 100.f );
            remove_dead_vertices( mesh, &logger );
        }
    }

    return anyMods;
}

bool remove_duplicate_halfedge_faces( frantic::geometry::trimesh3& mesh, bool fullyDegenerateOnly ) {
    frantic::logging::null_progress_logger logger;
    return remove_duplicate_halfedge_faces( mesh, logger, fullyDegenerateOnly );
}

void locate_sliver_vertices( const frantic::geometry::trimesh3& mesh, std::vector<size_t>& outSlivers,
                             double epsilon ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    const double epsilonSquared = frantic::math::square( epsilon );

    std::vector<char> vertexState( mesh.vertex_count(), char( false ) );
    outSlivers.clear();

    bool anyMods = false;

    for( size_t faceId = 0; faceId < mesh.face_count(); ++faceId ) {
        const vector3& face = mesh.get_face( faceId );

        for( size_t vertex = 0; vertex < 3; ++vertex ) {
            const vector3::value_type sVertIndex = face[vertex];
            const vector3::value_type tVertIndex = face[( vertex + 1 ) % 3];
            const vector3::value_type pVertIndex = face[( vertex + 2 ) % 3];

            const vector3f& sVert = mesh.get_vertex( sVertIndex );
            const vector3f& tVert = mesh.get_vertex( tVertIndex );
            const vector3f& pVert = mesh.get_vertex( pVertIndex );

            if( vector3f::distance_squared_double( sVert, tVert ) < epsilonSquared ) {
                // epsilon-sized edges
                vertexState[sVertIndex] = char( true );
                vertexState[tVertIndex] = char( true );
                anyMods = true;
            } else {
                // 'T'-vertices
                ray3f oppositeEdge = ray3f::from_line_segment( tVert, pVert );

                const float offsetValue = oppositeEdge.nearest_point_parameter( sVert );

                if( offsetValue > 0.0f && offsetValue < 1.0f &&
                    vector3f::distance_squared_double( oppositeEdge.at( offsetValue ), sVert ) < epsilonSquared ) {
                    vertexState[sVertIndex] = char( true );
                    anyMods = true;
                }
            }
        }
    }

    if( anyMods ) {
        for( size_t i = 0; i < vertexState.size(); ++i ) {
            if( vertexState[i] ) {
                outSlivers.push_back( i );
            }
        }
    }
}

bool relax_sliver_vertices( frantic::geometry::trimesh3& mesh, size_t iterations, float scale, float epsilon ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<size_t> sliverVertices;
    locate_sliver_vertices( mesh, sliverVertices, epsilon );

    if( sliverVertices.size() > 0 ) {
        relaxation::laplacian_smooth( mesh, iterations, scale, sliverVertices );
        return true;
    } else {
        return false;
    }
}

} // namespace geometry
} // namespace frantic
