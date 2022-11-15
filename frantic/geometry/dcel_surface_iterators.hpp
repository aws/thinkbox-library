// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <frantic/geometry/dcel_iterators.hpp>

#include <boost/foreach.hpp>

// These methods are a little more complex, this involve finding (if present) a surface (non-boundary) range around a
// vertex Note that the assumption is that there is at most one boundary range at each vertex (the alternative is not
// disallowed, but is really wonky and shouldn't happen in general)

namespace frantic {
namespace geometry {

/**
 * Returns the next boundary edge after the current vertex reference halfedge
 *
 * @param vertexRef a halfedge handle pointing to the current vertex
 * @return a halfedge handle to the first boundary handle pointing to the current vertex, or an invalid halfedge if
 * there is no boundary on the current vertex
 */
template <class Halfedge_t>
Halfedge_t next_boundary_edge( const Halfedge_t& vertexRef ) {
    BOOST_FOREACH( dcel::halfedge_handle he, dcel_vertex_cycle_range( vertexRef ) ) {
        if( he.is_boundary_face() ) {
            return he;
        }
    }
    return dcel::INVALID_HALFEDGE_HANDLE;
}

/**
 * Return the start of a range around the vertex which iterates through a contiguous set of surface faces
 */
template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t
dcel_vertex_surface_begin( const Halfedge_t& vertexRef ) {
    BOOST_FOREACH( dcel::halfedge_handle he, dcel_vertex_cycle_range( vertexRef ) ) {
        if( he.twin().is_boundary_face() ) {
            return dcel_vertex_begin( he );
        }
    }

    return dcel_vertex_begin( vertexRef );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t
dcel_vertex_surface_end( const Halfedge_t& vertexRef ) {
    Halfedge_t nextBoundary = next_boundary_edge( vertexRef );

    if( nextBoundary != dcel::INVALID_HALFEDGE_HANDLE ) {
        return dcel_vertex_end( nextBoundary );
    } else {
        return dcel_vertex_cycle_end( vertexRef );
    }
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_range_t
dcel_vertex_surface_range( const Halfedge_t& vertexRef ) {
    Halfedge_t nextBoundary = next_boundary_edge( vertexRef );

    if( nextBoundary != dcel::INVALID_HALFEDGE_HANDLE ) {
        return boost::make_iterator_range( dcel_vertex_begin( nextBoundary.vertex_next() ),
                                           dcel_vertex_end( nextBoundary ) );
    } else {
        return dcel_vertex_cycle_range( vertexRef );
    }
}

} // namespace geometry
} // namespace frantic
