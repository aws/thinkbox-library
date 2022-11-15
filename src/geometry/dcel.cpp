// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/dcel_iterators.hpp>
#include <frantic/geometry/dcel_surface_iterators.hpp>

#include <frantic/misc/exception_stream.hpp>
#include <frantic/misc/iterator.hpp>

#include <boost/algorithm/cxx11/copy_if.hpp>
#include <boost/foreach.hpp>
#include <boost/scoped_array.hpp>

#include <sstream>

namespace frantic {
namespace geometry {

const dcel::index_t dcel::INVALID_VERTEX_INDEX = index_t( -1 );
const dcel::index_t dcel::INVALID_HALFEDGE_INDEX = index_t( -1 );
const dcel::index_t dcel::INVALID_FACE_INDEX = index_t( -1 );
const dcel::halfedge_handle dcel::INVALID_HALFEDGE_HANDLE;

namespace {

dcel::halfedge_handle find_existing_previous( dcel::halfedge_handle base, dcel::halfedge_handle next ) {
    BOOST_FOREACH( dcel::halfedge_handle he, dcel_vertex_cycle_range( base ) ) {
        if( he.face_next() == next ) {
            return he;
        }
    }

    return dcel::INVALID_HALFEDGE_HANDLE;
}

void perform_insertion_fixup( dcel::halfedge_handle& addedEdge, dcel::halfedge_handle& expectedNext,
                              bool alreadyExisted ) {
    if( expectedNext != dcel::INVALID_HALFEDGE_HANDLE && addedEdge != expectedNext ) {
        dcel::halfedge_handle nextBound = next_boundary_edge( addedEdge.twin() );

        if( alreadyExisted ) {
            dcel::halfedge_handle addedPrev = find_existing_previous( addedEdge.twin(), addedEdge );
            dcel::halfedge_handle::link_edges( addedPrev, nextBound.face_next() );
        }

        dcel::halfedge_handle::link_edges( nextBound, expectedNext );
    }
}

} // namespace

class dcel::private_face_construction_info {
  public:
    void initialize() {
        expectedNextHalfedge = dcel::INVALID_HALFEDGE_HANDLE;
        previousHalfedge = dcel::INVALID_HALFEDGE_HANDLE;
        previousExisted = false;
        firstExisted = false;
        firstHalfedge = dcel::INVALID_HALFEDGE_HANDLE;
        count = 0;
        lastVertexId = dcel::INVALID_VERTEX_INDEX;
        firstVertexId = dcel::INVALID_VERTEX_INDEX;
        faceId = dcel::INVALID_FACE_INDEX;
    }

    index_t faceId;
    index_t lastVertexId;
    index_t firstVertexId;
    index_t count;
    dcel::halfedge_handle expectedNextHalfedge;
    dcel::halfedge_handle previousHalfedge;
    dcel::halfedge_handle firstHalfedge;
    bool previousExisted;
    bool firstExisted;
};

void dcel::halfedge_handle::link_edges( dcel::halfedge_handle prev, dcel::halfedge_handle next ) {
    if( prev.m_owner != next.m_owner ) {
        throw std::runtime_error(
            "dcel::halfedge_handle::link_edges : Error, cannot link edges from different structures" );
    }
    prev.m_owner->set_face_next( prev.m_index, next.m_index );
}

void dcel::halfedge_handle::insert_after( dcel::halfedge_handle handle ) {
    link_edges( handle, face_next() );
    link_edges( *this, handle.twin() );
}

dcel::dcel()
    : m_faceConstructionInfo( new private_face_construction_info ) {
    clear();
}

dcel::dcel( const dcel& other )
    : m_vertexHalfedges( other.m_vertexHalfedges )
    , m_faceHalfedges( other.m_faceHalfedges )
    , m_boundaryHalfedges( other.m_boundaryHalfedges )
    , m_halfedges( other.m_halfedges )
    , m_isAddingFace( other.m_isAddingFace )
    , m_isTriangleMesh( other.m_isTriangleMesh )
    , m_isBoundaryInfoCached( other.m_isBoundaryInfoCached )
    , m_faceConstructionInfo( new private_face_construction_info ) {}

dcel::dcel( size_t vertexCount, size_t faceCount )
    : m_faceConstructionInfo( new private_face_construction_info ) {
    initialize( vertexCount, faceCount );
}

dcel::~dcel() { delete m_faceConstructionInfo; }

void dcel::swap( dcel& other ) {
    if( this != &other ) {
        std::swap( m_vertexHalfedges, other.m_vertexHalfedges );
        std::swap( m_faceHalfedges, other.m_faceHalfedges );
        std::swap( m_boundaryHalfedges, other.m_boundaryHalfedges );
        std::swap( m_halfedges, other.m_halfedges );
        std::swap( m_isAddingFace, other.m_isAddingFace );
        std::swap( m_isTriangleMesh, other.m_isTriangleMesh );
        std::swap( m_isBoundaryInfoCached, other.m_isBoundaryInfoCached );
    }
}

void dcel::clear() {
    m_vertexHalfedges.clear();
    m_faceHalfedges.clear();
    m_boundaryHalfedges.clear();
    m_halfedges.clear();
    m_isAddingFace = false;
    m_isTriangleMesh = true;
    m_isBoundaryInfoCached = false;
}

void dcel::initialize( size_t vertexCount, size_t faceCount ) {
    m_vertexHalfedges.resize( vertexCount );
    m_vertexHalfedges.assign( vertexCount, INVALID_HALFEDGE_INDEX );
    m_faceHalfedges.resize( faceCount );
    m_faceHalfedges.assign( faceCount, INVALID_HALFEDGE_INDEX );
    m_boundaryHalfedges.clear();
    m_halfedges.clear();

    // rough approximation of (one of) Euler's Formula
    size_t edgeEstimate = vertexCount + faceCount;
    m_halfedges.reserve( edgeEstimate * 2 );
    m_isBoundaryInfoCached = false;
}

size_t dcel::vertex_count() const { return m_vertexHalfedges.size(); }

size_t dcel::face_count() const { return m_faceHalfedges.size(); }

size_t dcel::halfedge_count() const { return m_halfedges.size(); }

size_t dcel::boundary_count() const { return m_boundaryHalfedges.size(); }

bool dcel::is_triangle_mesh() const { return m_isTriangleMesh; }

bool dcel::is_boundary_vertex( size_t vertexId ) const {
    if( !has_vertex_halfedge( vertexId ) ) {
        return true;
    }

    BOOST_FOREACH( dcel::const_halfedge_handle handle, dcel_vertex_cycle_range( get_vertex_halfedge( vertexId ) ) ) {
        if( handle.is_boundary_edge() ) {
            return true;
        }
    }

    return false;
}

bool dcel::has_vertex_halfedge( size_t vertexId ) const {
    return m_vertexHalfedges.at( vertexId ) != INVALID_HALFEDGE_INDEX;
}

dcel::const_halfedge_handle dcel::get_vertex_halfedge( size_t vertexId ) const {
    return get_halfedge( m_vertexHalfedges.at( vertexId ) );
}

dcel::halfedge_handle dcel::get_vertex_halfedge( size_t vertexId ) {
    return get_halfedge( m_vertexHalfedges.at( vertexId ) );
}

void dcel::set_vertex_halfedge( size_t vertexId, dcel::const_halfedge_handle handle ) {
    m_vertexHalfedges.at( vertexId ) = handle.get_index();
}

void dcel::insert_vertex_halfedge( dcel::halfedge_handle handle ) {
    index_t newVertexId = index_t( m_vertexHalfedges.size() );
    m_vertexHalfedges.push_back( handle.get_index() );

    BOOST_FOREACH( halfedge_handle vertexIt, dcel_vertex_cycle_range( handle ) ) {
        vertexIt.set_target_vertex( newVertexId );
    }
}

void dcel::adjust_vertex_halfedge( size_t vertexId ) {
    halfedge_handle vertexHandle = get_vertex_halfedge( vertexId );

    if( vertexHandle != INVALID_HALFEDGE_HANDLE ) {
        BOOST_FOREACH( halfedge_handle vertexIt, dcel_vertex_cycle_range( vertexHandle ) ) {
            if( vertexIt.is_valid() ) {
                set_vertex_halfedge( vertexId, vertexIt );
                return;
            }
        }

        // if no valid halfedge was found, set the representative to the invalid index
        set_vertex_halfedge( vertexId, INVALID_HALFEDGE_HANDLE );
    }
}

bool dcel::has_face_halfedge( size_t faceId ) const { return m_faceHalfedges.at( faceId ) != INVALID_HALFEDGE_INDEX; }

dcel::const_halfedge_handle dcel::get_face_halfedge( size_t faceId ) const {
    return get_halfedge( m_faceHalfedges.at( faceId ) );
}

dcel::halfedge_handle dcel::get_face_halfedge( size_t faceId ) { return get_halfedge( m_faceHalfedges.at( faceId ) ); }

void dcel::set_face_halfedge( size_t faceId, dcel::const_halfedge_handle handle ) {
    m_faceHalfedges.at( faceId ) = handle.get_index();
}

void dcel::adjust_face_halfedge( size_t faceId ) {
    halfedge_handle faceHandle = get_face_halfedge( faceId );

    if( faceHandle != INVALID_HALFEDGE_HANDLE ) {
        BOOST_FOREACH( halfedge_handle faceIt, dcel_face_cycle_range( faceHandle ) ) {
            if( faceIt.is_valid() ) {
                set_face_halfedge( faceId, faceIt );
                return;
            }
        }

        // if no valid halfedge was found, set the representative to the invalid index
        set_face_halfedge( faceId, INVALID_HALFEDGE_HANDLE );
    }
}

bool dcel::has_boundary_halfedge( size_t boundaryId ) const {
    return m_boundaryHalfedges.at( boundaryId ) != INVALID_HALFEDGE_INDEX;
}

dcel::const_halfedge_handle dcel::get_boundary_halfedge( size_t boundaryId ) const {
    return get_halfedge( m_boundaryHalfedges.at( boundaryId ) );
}

dcel::halfedge_handle dcel::get_boundary_halfedge( size_t boundaryId ) {
    return get_halfedge( m_boundaryHalfedges.at( boundaryId ) );
}

void dcel::set_boundary_halfedge( size_t boundaryId, dcel::const_halfedge_handle handle ) {
    m_boundaryHalfedges.at( boundaryId ) = handle.get_index();
}

void dcel::insert_boundary_halfedge( dcel::halfedge_handle handle ) {
    m_boundaryHalfedges.push_back( handle.get_index() );
}

dcel::const_halfedge_handle dcel::get_halfedge( size_t halfedgeId ) const {
    if( halfedgeId == INVALID_HALFEDGE_INDEX ) {
        return INVALID_HALFEDGE_HANDLE;
    }
    return const_halfedge_handle( this, index_t( halfedgeId ) );
}

dcel::halfedge_handle dcel::get_halfedge( size_t halfedgeId ) {
    if( halfedgeId == INVALID_HALFEDGE_INDEX ) {
        return INVALID_HALFEDGE_HANDLE;
    }
    return halfedge_handle( this, index_t( halfedgeId ) );
}

std::pair<dcel::halfedge_handle, dcel::halfedge_handle> dcel::insert_edge() {
    index_t newIndex = index_t( m_halfedges.size() );
    m_halfedges.resize( newIndex + 2 );
    return std::make_pair( get_halfedge( newIndex ), get_halfedge( newIndex + 1 ) );
}

dcel::const_halfedge_handle dcel::find_vertex_halfedge( size_t sourceVertexId, size_t destVertexId ) const {
    dcel::const_halfedge_handle vertexStart = get_vertex_halfedge( destVertexId );
    if( vertexStart != INVALID_HALFEDGE_HANDLE ) {
        BOOST_FOREACH( dcel::const_halfedge_handle handle, dcel_vertex_cycle_range( vertexStart ) ) {
            if( handle.source_vertex() == sourceVertexId ) {
                return handle;
            }
        }
    }
    return INVALID_HALFEDGE_HANDLE;
}

dcel::halfedge_handle dcel::find_vertex_halfedge( size_t sourceVertexId, size_t destVertexId ) {
    const_halfedge_handle handle =
        const_cast<const dcel*>( this )->find_vertex_halfedge( sourceVertexId, destVertexId );

    if( handle != INVALID_HALFEDGE_HANDLE ) {
        return halfedge_handle( this, handle.get_index() );
    } else {
        return INVALID_HALFEDGE_HANDLE;
    }
}

dcel::const_halfedge_handle dcel::find_face_halfedge( size_t face, size_t oppositeFace ) const {
    const_halfedge_handle faceStart = get_face_halfedge( face );
    if( faceStart != INVALID_HALFEDGE_HANDLE ) {
        BOOST_FOREACH( const_halfedge_handle handle, dcel_face_cycle_range( faceStart ) ) {
            if( handle.opposite_face() == oppositeFace ) {
                return handle;
            }
        }
    }
    return INVALID_HALFEDGE_HANDLE;
}

dcel::halfedge_handle dcel::find_face_halfedge( size_t face, size_t oppositeFace ) {
    const_halfedge_handle handle = const_cast<const dcel*>( this )->find_face_halfedge( face, oppositeFace );

    if( handle != INVALID_HALFEDGE_HANDLE ) {
        return halfedge_handle( this, handle.get_index() );
    } else {
        return INVALID_HALFEDGE_HANDLE;
    }
}

size_t dcel::add_face3( size_t v0, size_t v1, size_t v2 ) {
    index_t newFaceId = insert_new_face_id();
    add_face3( newFaceId, v0, v1, v2 );
    return newFaceId;
}

void dcel::add_face3( size_t faceId, size_t v0, size_t v1, size_t v2 ) {
    index_t list[3] = { index_t( v0 ), index_t( v1 ), index_t( v2 ) };
    add_face( faceId, list, list + 3 );
}

size_t dcel::add_face4( size_t v0, size_t v1, size_t v2, size_t v3 ) {
    size_t newFaceId = insert_new_face_id();
    add_face4( newFaceId, v0, v1, v2, v3 );
    return newFaceId;
}

void dcel::add_face4( size_t faceId, size_t v0, size_t v1, size_t v2, size_t v3 ) {
    index_t list[4] = { index_t( v0 ), index_t( v1 ), index_t( v2 ), index_t( v3 ) };
    add_face( faceId, list, list + 4 );
}

void dcel::begin_face( size_t faceId ) {
    if( m_isAddingFace ) {
        throw std::runtime_error( "dcel::begin_face: already adding a face." );
    }
    m_isAddingFace = true;

    m_faceConstructionInfo->initialize();
    m_faceConstructionInfo->faceId = index_t( faceId );
}

size_t dcel::begin_face() {
    size_t newFaceId = insert_new_face_id();
    begin_face( newFaceId );
    return newFaceId;
}

void dcel::add_face_vertex( size_t vertexId ) {
    if( m_faceConstructionInfo->count > 0 ) {
        add_edge_internal( m_faceConstructionInfo->faceId, m_faceConstructionInfo->lastVertexId, index_t( vertexId ) );
    } else {
        m_faceConstructionInfo->firstVertexId = index_t( vertexId );
    }
    m_faceConstructionInfo->lastVertexId = index_t( vertexId );
    ++m_faceConstructionInfo->count;
}

void dcel::end_face() {
    if( !m_isAddingFace ) {
        throw std::runtime_error( "dcel::end_face: face was not completed." );
    }

    if( m_faceConstructionInfo->count < 3 ) {
        throw std::runtime_error( "dcel::end_face : face has fewer than 3 vertices." );
    } else if( m_faceConstructionInfo->count > 3 ) {
        m_isTriangleMesh = false;
    }

    add_edge_internal( m_faceConstructionInfo->faceId, m_faceConstructionInfo->lastVertexId,
                       m_faceConstructionInfo->firstVertexId );

    perform_insertion_fixup( m_faceConstructionInfo->firstHalfedge, m_faceConstructionInfo->expectedNextHalfedge,
                             m_faceConstructionInfo->firstExisted );
    halfedge_handle::link_edges( m_faceConstructionInfo->previousHalfedge, m_faceConstructionInfo->firstHalfedge );

    m_isAddingFace = false;
}

dcel::index_t dcel::get_face_prev( dcel::index_t halfedgeId ) const {

    const_halfedge_handle currentHandle = get_halfedge( halfedgeId );

    if( currentHandle.is_boundary_face() ) {
        const_halfedge_handle lastEdgeBefore;
        BOOST_FOREACH( const_halfedge_handle edge, dcel_vertex_cycle_range( currentHandle.twin() ) ) {
            lastEdgeBefore = edge;
        }
        return lastEdgeBefore.get_index();
    } else {
        const_halfedge_handle lastEdgeBefore;
        BOOST_FOREACH( const_halfedge_handle edge, dcel_face_cycle_range( currentHandle ) ) {
            lastEdgeBefore = edge;
        }
        return lastEdgeBefore.get_index();
    }
}

dcel::index_t dcel::get_vertex_prev( dcel::index_t halfedgeId ) const {
    const_halfedge_handle lastEdgeBefore;

    BOOST_FOREACH( const_halfedge_handle edge, dcel_vertex_cycle_range( get_halfedge( halfedgeId ) ) ) {
        lastEdgeBefore = edge;
    }

    return lastEdgeBefore.get_index();
}

dcel::index_t dcel::insert_new_face_id() {
    index_t faceId = index_t( m_faceHalfedges.size() );
    m_faceHalfedges.push_back( INVALID_HALFEDGE_INDEX );
    return faceId;
}

dcel::halfedge_handle dcel::create_halfedge_internal( index_t face, index_t startVertex, index_t endVertex,
                                                      bool& alreadyExisted ) {
    // check if this edge's opposite has already been created
    dcel::halfedge_handle existingOpposite = find_vertex_halfedge( endVertex, startVertex );

    dcel::halfedge_handle added;

    if( existingOpposite.is_valid() ) {
        added = existingOpposite.twin();

        if( added.current_face() != INVALID_FACE_INDEX ) {
            std::ostringstream oss;
            oss << "halfedge_structure::add_halfedge: Not a manifold mesh, halfedge (" << startVertex << ","
                << endVertex << ") appeared twice (faces " << added.current_face() << " and " << face << ").";
            throw std::runtime_error( oss.str() );
        }

        alreadyExisted = true;
    } else {
        size_t addedHalfedge = m_halfedges.size();
        m_halfedges.resize( addedHalfedge + 2 );
        added = get_halfedge( index_t( addedHalfedge ) );
        added.set_target_vertex( endVertex );
        dcel::halfedge_handle twin = get_halfedge( index_t( addedHalfedge + 1 ) );
        twin.set_target_vertex( startVertex );
        halfedge_handle::link_edges( added, twin );
        halfedge_handle::link_edges( twin, added );
        alreadyExisted = false;
    }

    added.set_current_face( face );
    return added;
}

void dcel::add_edge_internal( index_t faceId, index_t startVertex, index_t endVertex ) {
    bool alreadyExisted;
    halfedge_handle addedEdge = create_halfedge_internal( faceId, startVertex, endVertex, alreadyExisted );

    if( !m_faceConstructionInfo->firstHalfedge.is_valid() ) {
        m_faceConstructionInfo->firstHalfedge = addedEdge;
        m_faceConstructionInfo->firstExisted = alreadyExisted;
    }

    // link the tip of the current node
    if( m_vertexHalfedges.at( addedEdge.target_vertex() ) == INVALID_HALFEDGE_INDEX ) {
        m_vertexHalfedges.at( addedEdge.target_vertex() ) = addedEdge.get_index();
    } else if( !alreadyExisted ) {
        halfedge_handle nextBoundary =
            next_boundary_edge( get_halfedge( m_vertexHalfedges.at( addedEdge.target_vertex() ) ) );
        if( !nextBoundary.is_valid() ) {
            std::ostringstream oss;
            oss << "halfedge_structure::add_face : Not a manifold mesh. Vertex " << addedEdge.target_vertex()
                << " has two or more indicent local components.";
            throw std::runtime_error( oss.str() );
        }
        nextBoundary.insert_after( addedEdge );
    }

    perform_insertion_fixup( addedEdge, m_faceConstructionInfo->expectedNextHalfedge, alreadyExisted );

    if( m_faceConstructionInfo->previousHalfedge.is_valid() ) {
        halfedge_handle::link_edges( m_faceConstructionInfo->previousHalfedge, addedEdge );
    }

    if( m_faceHalfedges.at( addedEdge.current_face() ) == INVALID_HALFEDGE_INDEX ) {
        m_faceHalfedges.at( addedEdge.current_face() ) = addedEdge.get_index();
    }

    m_faceConstructionInfo->expectedNextHalfedge = addedEdge.face_next();
    m_faceConstructionInfo->previousHalfedge = addedEdge;
    m_faceConstructionInfo->previousExisted = alreadyExisted;
}

void dcel::collapse_edge( dcel::halfedge_handle handle, bool checkManifold ) {
    if( !collapse_edge_internal( handle, checkManifold ) ) {
        throw std::runtime_error( "dcel::collapse_edge : Collapsing this edge would result in a non-manifold mesh." );
    }
}

bool dcel::try_collapse_edge( dcel::halfedge_handle handle ) { return collapse_edge_internal( handle, true ); }

namespace {
struct apex_filter {
    apex_filter()
        : filterVertexCount( 0 ) {}

    dcel::index_t vertices[2];
    size_t filterVertexCount;
    bool operator()( dcel::index_t value ) const {
        for( size_t i = 0; i < filterVertexCount; ++i ) {
            if( value == vertices[i] ) {
                return false;
            }
        }
        return true;
    }
};
} // namespace

// The goal here is to make this method general (i.e. handle any possible case),
// while keeping it efficient in the general case (where the vertex incidence can
// generally be seen as constant).
bool dcel::check_manifold_collapse( dcel::const_halfedge_handle edgeHandle ) const {
    // The expected incidence of a given vertex is around 6
    const size_t expectedBufferSize = 128;

    if( !edgeHandle.is_valid() ) {
        throw std::runtime_error( "dcel::check_manifold_collapse : halfedge is not valid" );
    }

    if( ( is_boundary_vertex( edgeHandle.target_vertex() ) && is_boundary_vertex( edgeHandle.source_vertex() ) ) &&
        !edgeHandle.is_boundary_edge() ) {
        return false;
    }

    dcel::const_halfedge_handle twinEdgeHandle = edgeHandle.twin();

    dcel::index_t constBufferTarget[expectedBufferSize];
    dcel::index_t constBufferSource[expectedBufferSize];

    dcel::index_t* targetAdjacencies = constBufferTarget;
    dcel::index_t* sourceAdjacencies = constBufferSource;

    boost::scoped_array<dcel::index_t> dynamicBufferTarget;
    boost::scoped_array<dcel::index_t> dynamicBufferSource;

    // if either side of the current edge is a triangular face, exclude the triangle
    // apex from the set of adjacent vertices to check
    apex_filter filter;
    filter.filterVertexCount = 0;

    if( edgeHandle.face_next().face_next().face_next() == edgeHandle ) {
        filter.vertices[filter.filterVertexCount] = edgeHandle.face_next().target_vertex();
        ++filter.filterVertexCount;
    }

    if( twinEdgeHandle.face_next().face_next().face_next() == twinEdgeHandle ) {
        filter.vertices[filter.filterVertexCount] = twinEdgeHandle.face_next().target_vertex();
        ++filter.filterVertexCount;
    }

    const dcel::index_t targetVertex = edgeHandle.target_vertex();
    const dcel::index_t sourceVertex = twinEdgeHandle.target_vertex();

    limited_output_iterator<dcel::index_t*> targetAdjacenciesEnd = boost::algorithm::copy_if(
        dcel_vertex_adjacency_cycle_begin( edgeHandle ), dcel_vertex_adjacency_cycle_end( edgeHandle ),
        make_limited_output_iterator( targetAdjacencies, expectedBufferSize ), filter );

    if( targetAdjacenciesEnd.count() > expectedBufferSize ) {
        dynamicBufferTarget.reset( new dcel::index_t[targetAdjacenciesEnd.count()] );
        targetAdjacencies = dynamicBufferTarget.get();
        boost::algorithm::copy_if( dcel_vertex_adjacency_cycle_begin( edgeHandle ),
                                   dcel_vertex_adjacency_cycle_end( edgeHandle ), targetAdjacencies, filter );
    }

    dcel::const_halfedge_handle edgeTwinHandle = edgeHandle.twin();

    limited_output_iterator<dcel::index_t*> sourceAdjacenciesEnd = boost::algorithm::copy_if(
        dcel_vertex_adjacency_cycle_begin( edgeTwinHandle ), dcel_vertex_adjacency_cycle_end( edgeTwinHandle ),
        make_limited_output_iterator( sourceAdjacencies, expectedBufferSize ), filter );

    if( sourceAdjacenciesEnd.count() > expectedBufferSize ) {
        dynamicBufferSource.reset( new dcel::index_t[sourceAdjacenciesEnd.count()] );
        sourceAdjacencies = dynamicBufferSource.get();
        boost::algorithm::copy_if( dcel_vertex_adjacency_cycle_begin( edgeTwinHandle ),
                                   dcel_vertex_adjacency_cycle_end( edgeTwinHandle ), sourceAdjacencies, filter );
    }

    if( ( targetAdjacenciesEnd.count() == 1 && targetAdjacencies[0] == sourceVertex ) &&
        ( sourceAdjacenciesEnd.count() == 1 && sourceAdjacencies[0] == targetVertex ) ) {
        return false;
    }

    std::sort( targetAdjacencies, targetAdjacencies + targetAdjacenciesEnd.count() );
    std::sort( sourceAdjacencies, sourceAdjacencies + sourceAdjacenciesEnd.count() );

    count_output_iterator counter =
        std::set_intersection( targetAdjacencies, targetAdjacencies + targetAdjacenciesEnd.count(), sourceAdjacencies,
                               sourceAdjacencies + sourceAdjacenciesEnd.count(), count_output_iterator() );

    if( counter.count() == 0 ) {
        return true;
    } else {
        return false;
    }
}

// TODO: this isn't terribly efficient
// Part of the reason why is that we are assuming we are starting from scratch in terms of information
// (i.e. no existing boundary information is known, other than all halfedges with face references out-
// side fo the face count are boundary edges)
// One way to speed this up at the expense of storage (and code complexity) would be to use disjoint
// sets at all times to ensure all boundary information is up to date
void dcel::cache_boundary_info( bool force ) {
    if( m_isBoundaryInfoCached && !force ) {
        return;
    }

    // first pass, normalize all out of range faces to be the invalid index
    for( size_t i = 0; i < halfedge_count(); ++i ) {
        if( is_boundary_face( m_halfedges.at( i ).m_face ) ) {
            m_halfedges.at( i ).m_face = INVALID_FACE_INDEX;
        }
    }

    m_boundaryHalfedges.clear();

    // second pass, loop through each boundary, assigning it a unique id
    for( size_t i = 0; i < halfedge_count(); ++i ) {
        if( m_halfedges.at( i ).m_face == INVALID_FACE_INDEX ) {
            size_t loopCounter = 0;

            index_t currentBoundaryId = index_t( m_boundaryHalfedges.size() );
            m_boundaryHalfedges.push_back( index_t( i ) );

            BOOST_FOREACH( halfedge_handle handle, dcel_face_cycle_range( get_halfedge( i ) ) ) {
                if( handle.current_boundary() != INVALID_FACE_INDEX ) {
                    throw std::runtime_error(
                        "dcel::cache_boundary_info : Error, non-boundary edge on boundary detected." );
                }

                handle.set_current_boundary( currentBoundaryId );

                // simple infinite loop checker
                ++loopCounter;
                if( loopCounter > halfedge_count() ) {
                    throw std::runtime_error( "dcel::cache_boundary_info : Error, invalid boundary loop detected." );
                }
            }
        }
    }

    m_isBoundaryInfoCached = true;
}

dcel::private_face_construction_info* dcel::create_private_face_construction_info() {
    return new private_face_construction_info;
}

namespace {

void remove_spike( dcel::halfedge_handle curr ) {
    dcel::halfedge_handle twinPrev = curr.face_prev();
    dcel::halfedge_handle currNext = curr.twin().face_next();
    dcel::halfedge_handle::link_edges( twinPrev, currNext );
}

void remove_edge( dcel::halfedge_handle curr, dcel& owner ) {
    dcel::index_t targetVertex = curr.target_vertex();
    dcel::index_t sourceVertex = curr.source_vertex();

    dcel::halfedge_handle twin = curr.twin();

    dcel::index_t currentFace = curr.current_face();
    dcel::index_t twinFace = twin.current_face();

    const bool floatingTarget = curr.face_next() == twin;
    const bool floatingSource = twin.face_next() == curr;

    if( !floatingSource ) {
        remove_spike( curr );
    }

    if( !floatingTarget ) {
        remove_spike( twin );
    }

    curr.invalidate();
    twin.invalidate();

    // Re-adjust the vertex references, unless the endpoint was floating
    if( !floatingTarget ) {
        owner.adjust_vertex_halfedge( targetVertex );
    } else {
        owner.set_vertex_halfedge( targetVertex, dcel::INVALID_HALFEDGE_HANDLE );
    }

    if( !floatingSource ) {
        owner.adjust_vertex_halfedge( sourceVertex );
    } else {
        owner.set_vertex_halfedge( sourceVertex, dcel::INVALID_HALFEDGE_HANDLE );
    }

    // if neither edge is a boundary, this is now a unified face
    if( !curr.is_boundary_edge() && !twin.is_boundary_edge() ) {

        if( currentFace != twinFace ) {
            BOOST_FOREACH( dcel::halfedge_handle handle, dcel_face_cycle_range( curr ) ) {
                handle.set_current_face( currentFace );
            }
            owner.set_face_halfedge( twinFace, dcel::INVALID_HALFEDGE_HANDLE );
        }

        owner.adjust_face_halfedge( currentFace );
    } else {
        // Otherwise, this is now part of the boundary

        if( !floatingSource && !floatingTarget ) {
            dcel::index_t boundaryFace = curr.is_boundary_edge() ? currentFace : twinFace;

            BOOST_FOREACH( dcel::halfedge_handle handle, dcel_face_cycle_range( curr.face_next() ) ) {
                handle.set_current_face( boundaryFace );
            }
        }

        if( !curr.is_boundary_edge() ) {
            owner.set_face_halfedge( currentFace, dcel::INVALID_HALFEDGE_HANDLE );
        }

        if( !twin.is_boundary_edge() ) {
            owner.set_face_halfedge( twinFace, dcel::INVALID_HALFEDGE_HANDLE );
        }
    }
}

void remove_loop( dcel::halfedge_handle current, dcel& owner ) {
    dcel::index_t targetVertex = current.target_vertex();
    dcel::index_t sourceVertex = current.source_vertex();

    dcel::halfedge_handle next = current.face_next();

    if( next.face_next() != current ) {
        throw std::runtime_error( "remove_loop : Error, was not a loop." );
    }

    dcel::halfedge_handle currTwin = current.twin();
    dcel::halfedge_handle nextTwin = next.twin();

    dcel::halfedge_handle nextTwinNext = nextTwin.face_next();
    dcel::halfedge_handle nextTwinPrev = nextTwin.face_prev();

    dcel::halfedge_handle::link_edges( current, nextTwinNext );
    dcel::halfedge_handle::link_edges( nextTwinPrev, current );

    if( !current.is_boundary_face() ) {
        owner.set_face_halfedge( current.current_face(), dcel::INVALID_HALFEDGE_HANDLE );
    }

    current.set_current_face( nextTwin.current_face() );

    next.set_current_face( dcel::INVALID_HALFEDGE_HANDLE );
    nextTwin.set_current_face( dcel::INVALID_HALFEDGE_HANDLE );
    next.invalidate();
    nextTwin.invalidate();

    if( !current.is_boundary_face() ) {
        owner.adjust_face_halfedge( current.current_face() );
    }

    const bool removeEdge = ( currTwin.is_boundary_face() && current.is_boundary_face() ) ||
                            currTwin.face_next() == current || current.face_next() == currTwin;

    if( removeEdge ) {
        remove_edge( current, owner );
    } else {
        owner.adjust_vertex_halfedge( targetVertex );
        owner.adjust_vertex_halfedge( sourceVertex );
    }
}

} // namespace

bool dcel::collapse_edge_internal( dcel::halfedge_handle edgeHandle, bool checkManifold ) {
    if( !edgeHandle.is_valid() ) {
        throw std::runtime_error( "dcel::collapse_edge_internal : halfedge is not valid" );
    }

    if( checkManifold && !check_manifold_collapse( edgeHandle ) ) {
        return false;
    }

    dcel::halfedge_handle edgeHandleTwin = edgeHandle.twin();

    dcel::index_t targetVertex = edgeHandle.target_vertex();
    dcel::index_t sourceVertex = edgeHandleTwin.target_vertex();

    BOOST_FOREACH( dcel::halfedge_handle it, dcel_vertex_cycle_range( edgeHandleTwin ) ) {
        it.set_target_vertex( targetVertex );
    }

    halfedge_handle edgeNext = edgeHandle.face_next();
    halfedge_handle twinNext = edgeHandleTwin.face_next();

    dcel::halfedge_handle handles[2] = { edgeHandle, edgeHandle.twin() };

    for( size_t i = 0; i < 2; ++i ) {
        // find the edge that is previous to the collapsed edge on the current face
        // N.B. this is an example of where having a previous pointer on the
        //   halfedges would be faster, though I'm not sure if the extra memory is worth it?
        dcel::halfedge_handle prevEdge = handles[i].face_prev();
        dcel::halfedge_handle nextEdge = handles[i].face_next();

        if( !prevEdge.is_valid() || !nextEdge.is_valid() ) {
            throw std::runtime_error(
                "perform_half_collapse : Error, incident edges to collapse edge were not valid." );
        }

        dcel::halfedge_handle::link_edges( prevEdge, nextEdge );

        handles[i].invalidate();

        if( !handles[i].is_boundary_face() ) {
            adjust_face_halfedge( handles[i].current_face() );
        }
    }

    adjust_vertex_halfedge( targetVertex );
    set_vertex_halfedge( sourceVertex, INVALID_HALFEDGE_HANDLE );

    if( edgeNext.face_next().face_next() == edgeNext ) {
        remove_loop( edgeNext, *this );
    }

    if( twinNext.face_next().face_next() == twinNext ) {
        remove_loop( twinNext, *this );
    }

    return true;
}

void dcel::flip_edge( halfedge_handle edge ) {
    if( edge.is_boundary_edge() ) {
        throw std::runtime_error( "dcel::flip_edge : Error, cannot flip, the edge (" +
                                  boost::lexical_cast<std::string>( edge.source_vertex() ) + "," +
                                  boost::lexical_cast<std::string>( edge.target_vertex() ) + ") was a boundary edge." );
    }

    halfedge_handle edgeTwin = edge.twin();

    index_t f0 = edge.current_face();
    halfedge_handle f0Next = edge.face_next();
    halfedge_handle f0Prev = f0Next.face_next();

    if( f0Prev.face_next() != edge ) {
        throw std::runtime_error( "dcel::flip_edge : Error, cannot flip, incident face " +
                                  boost::lexical_cast<std::string>( f0 ) + " was not a triangle." );
    }

    index_t f1 = edgeTwin.current_face();
    halfedge_handle f1Next = edgeTwin.face_next();
    halfedge_handle f1Prev = f1Next.face_next();

    if( f1Prev.face_next() != edgeTwin ) {
        throw std::runtime_error( "dcel::flip_edge : Error, cannot flip, incident face " +
                                  boost::lexical_cast<std::string>( f1 ) + " was not a triangle." );
    }

    if( find_vertex_halfedge( f0Next.target_vertex(), f1Next.target_vertex() ) != INVALID_HALFEDGE_HANDLE ) {
        throw std::runtime_error( "dcel::flip_edge : Error, cannot flip, the opposite edge (" +
                                  boost::lexical_cast<std::string>( f0Next.target_vertex() ) + "," +
                                  boost::lexical_cast<std::string>( f1Next.target_vertex() ) + ") already exists." );
    }

    // update the vertex refs for the two halfedges which will be moved
    if( m_vertexHalfedges.at( edge.target_vertex() ) == edge.get_index() ) {
        m_vertexHalfedges.at( edge.target_vertex() ) = edge.vertex_next().get_index();
    }

    edge.set_target_vertex( f0Next.target_vertex() );

    if( m_vertexHalfedges.at( edgeTwin.target_vertex() ) == edgeTwin.get_index() ) {
        m_vertexHalfedges.at( edgeTwin.target_vertex() ) = edgeTwin.vertex_next().get_index();
    }

    edgeTwin.set_target_vertex( f1Next.target_vertex() );

    // update the face refs for the two face-shifted edges
    if( m_faceHalfedges.at( f0Next.current_face() ) == f0Next.get_index() ) {
        m_faceHalfedges.at( f0Next.current_face() ) = edge.get_index();
    }

    f0Next.set_current_face( f1 );

    if( m_faceHalfedges.at( f1Next.current_face() ) == f1Next.get_index() ) {
        m_faceHalfedges.at( f1Next.current_face() ) = edgeTwin.get_index();
    }

    f1Next.set_current_face( f0 );

    // finally, re-link the corresponding edges
    halfedge_handle::link_edges( edge, f0Prev );
    halfedge_handle::link_edges( f0Prev, f1Next );
    halfedge_handle::link_edges( f1Next, edge );

    halfedge_handle::link_edges( edgeTwin, f1Prev );
    halfedge_handle::link_edges( f1Prev, f0Next );
    halfedge_handle::link_edges( f0Next, edgeTwin );
}

} // namespace geometry
} // namespace frantic
