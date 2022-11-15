// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <boost/move/move.hpp>

#include <vector>

namespace frantic {
namespace geometry {

/**
 * Implementation of the Doubly-Connected Edge List data structure for
 * representing a planar surface mesh. The objective of this structure
 * is to provide constant (or near-constant) lookup for any local
 * adjacency property in the mesh, given a reference to a directed
 * edge (aka halfedge).
 *
 * Note that for space efficiency reasons, this implementation is not
 * doubly-connected (saves about 1/6th memory usage), at a slight cost
 * to reverse-iteration through face adjacency (which is a very uncommon
 * operation).
 * TODO: consider adding this as an option to the structure
 *
 * To save space (another 1/6th) and also allow consistent edge indexing,
 * all halfedges are located adjacent in memory (such that edgeId == halfedgeId / 2).
 *
 * All accesses to the internal structure are done through const_halfedge_handle
 * or halfedge_handle, which are a wrapper around a reference to the dcel, and
 * the index of the halfedge. This was to allow the structure to store all data
 * as a single contiguous array (which makes allocation and de-allocation more efficient),
 * while still allowing incremental construction. Handle objects become invalidated
 * when the dcel is destroyed.
 */
class dcel {
  public:
    // Currently it is explicitly using a 32-bit index to save on space.
    // Ideally, this could be a template parameter to allow a larger index space.
    typedef boost::uint32_t index_t;

  private:
    BOOST_COPYABLE_AND_MOVABLE( dcel )

    struct halfedge {
      public:
        halfedge()
            : m_next( INVALID_HALFEDGE_INDEX )
            , m_vertex( INVALID_VERTEX_INDEX )
            , m_face( INVALID_FACE_INDEX ) {}

        index_t m_next;
        index_t m_vertex;
        index_t m_face;
    };

  public:
    class const_halfedge_handle;
    class halfedge_handle;
    friend class const_halfedge_handle;
    friend class halfedge_handle;

    static const index_t INVALID_HALFEDGE_INDEX;
    static const index_t INVALID_FACE_INDEX;
    static const index_t INVALID_VERTEX_INDEX;

    class const_halfedge_handle {
      public:
        const_halfedge_handle()
            : m_owner( NULL )
            , m_index( INVALID_HALFEDGE_INDEX ) {}

        const_halfedge_handle( const dcel* owner, index_t index )
            : m_owner( owner )
            , m_index( index ) {}

        const_halfedge_handle( const const_halfedge_handle& other )
            : m_owner( other.m_owner )
            , m_index( other.m_index ) {}

        const_halfedge_handle( const halfedge_handle& other );

        const_halfedge_handle& operator=( const halfedge_handle& other );

        const_halfedge_handle& operator=( const const_halfedge_handle& other ) {
            m_owner = other.m_owner;
            m_index = other.m_index;
            return *this;
        }

        operator bool() const { return m_owner != NULL && m_index != INVALID_HALFEDGE_INDEX; }

        bool operator==( const const_halfedge_handle& other ) const {
            return m_owner == other.m_owner && m_index == other.m_index;
        }

        bool operator!=( const const_halfedge_handle& other ) const { return !( *this == other ); }

        bool operator==( const halfedge_handle& other ) const;

        bool operator!=( const halfedge_handle& other ) const { return !( *this == other ); }

        index_t get_index() const { return m_index; }

        index_t get_edge_index() const { return m_owner->get_edge( m_index ); }

        const dcel* get_owner() const { return m_owner; }

        const_halfedge_handle face_next() const {
            return const_halfedge_handle( m_owner, m_owner->get_face_next( m_index ) );
        }

        const_halfedge_handle face_prev() const {
            return const_halfedge_handle( m_owner, m_owner->get_face_prev( m_index ) );
        }

        const_halfedge_handle vertex_next() const {
            return const_halfedge_handle( m_owner, m_owner->get_vertex_next( m_index ) );
        }

        const_halfedge_handle vertex_prev() const {
            return const_halfedge_handle( m_owner, m_owner->get_vertex_prev( m_index ) );
        }

        const_halfedge_handle twin() const { return const_halfedge_handle( m_owner, m_owner->get_twin( m_index ) ); }

        index_t target_vertex() const { return m_owner->get_target_vertex( m_index ); }

        index_t source_vertex() const { return m_owner->get_source_vertex( m_index ); }

        index_t current_face() const { return m_owner->get_current_face( m_index ); }

        index_t opposite_face() const { return m_owner->get_opposite_face( m_index ); }

        const_halfedge_handle canonical_edge() const {
            if( is_canonical_edge() ) {
                return *this;
            } else {
                return twin();
            }
        }

        bool is_canonical_edge() const { return source_vertex() > target_vertex(); }

        bool is_valid() const { return m_owner && m_owner->is_halfedge_valid( m_index ); }

        bool is_boundary_face() const { return m_owner->is_boundary_edge_face( m_index ); }

        bool is_boundary_edge() const { return m_owner->is_boundary_edge( m_index ); }

        index_t current_boundary() const { return m_owner->get_current_boundary( m_index ); }

      private:
        const dcel* m_owner;
        index_t m_index;
    };

    class halfedge_handle {
      public:
        halfedge_handle()
            : m_owner( NULL )
            , m_index( INVALID_HALFEDGE_INDEX ) {}

        halfedge_handle( dcel* owner, index_t index )
            : m_owner( owner )
            , m_index( index ) {}

        halfedge_handle( const halfedge_handle& other )
            : m_owner( other.m_owner )
            , m_index( other.m_index ) {}

        const dcel* get_owner() const { return m_owner; }

        index_t get_index() const { return m_index; }

        index_t get_edge_index() const { return m_owner->get_edge( m_index ); }

        operator bool() const { return m_owner != NULL && m_index != INVALID_HALFEDGE_INDEX; }

        halfedge_handle& operator=( const halfedge_handle& other ) {
            m_owner = other.m_owner;
            m_index = other.m_index;
            return *this;
        }

        bool operator==( const const_halfedge_handle& other ) const {
            return m_owner == other.get_owner() && m_index == other.get_index();
        }

        bool operator==( const halfedge_handle& other ) const {
            return m_owner == other.m_owner && m_index == other.m_index;
        }

        bool operator!=( const halfedge_handle& other ) const { return !( *this == other ); }

        bool operator!=( const const_halfedge_handle& other ) const { return !( *this == other ); }

        halfedge_handle face_next() const { return halfedge_handle( m_owner, m_owner->get_face_next( m_index ) ); }

        halfedge_handle face_prev() const { return halfedge_handle( m_owner, m_owner->get_face_prev( m_index ) ); }

        halfedge_handle vertex_next() const { return halfedge_handle( m_owner, m_owner->get_vertex_next( m_index ) ); }

        halfedge_handle vertex_prev() const { return halfedge_handle( m_owner, m_owner->get_vertex_prev( m_index ) ); }

        halfedge_handle twin() const { return halfedge_handle( m_owner, m_owner->get_twin( m_index ) ); }

        index_t target_vertex() const { return m_owner->get_target_vertex( m_index ); }

        void set_target_vertex( index_t targetId ) { m_owner->set_target_vertex( m_index, targetId ); }

        index_t source_vertex() const { return m_owner->get_source_vertex( m_index ); }

        index_t current_face() const { return m_owner->get_current_face( m_index ); }

        void set_current_face( index_t faceId ) { m_owner->set_current_face( m_index, faceId ); }

        index_t opposite_face() const { return m_owner->get_opposite_face( m_index ); }

        halfedge_handle canonical_edge() const {
            if( is_canonical_edge() ) {
                return *this;
            } else {
                return twin();
            }
        }

        bool is_canonical_edge() const { return source_vertex() > target_vertex(); }

        bool is_valid() const { return m_owner && m_owner->is_halfedge_valid( m_index ); }

        void invalidate() {
            if( is_valid() ) {
                set_target_vertex( INVALID_VERTEX_INDEX );
            }
        }

        bool is_boundary_face() const { return m_owner->is_boundary_edge_face( m_index ); }

        bool is_boundary_edge() const { return m_owner->is_boundary_edge( m_index ); }

        index_t current_boundary() const { return m_owner->get_current_boundary( m_index ); }

        void set_current_boundary( index_t boundaryId ) const {
            return m_owner->set_current_boundary( m_index, boundaryId );
        }

        static void link_edges( halfedge_handle prev, halfedge_handle next );
        void insert_after( halfedge_handle edge );

      private:
        dcel* m_owner;
        index_t m_index;
    };

    static const halfedge_handle INVALID_HALFEDGE_HANDLE;

  public:
    dcel();

    dcel( const dcel& other );

    dcel( BOOST_RV_REF( dcel ) other )
        : m_faceConstructionInfo( create_private_face_construction_info() ) {
        swap( other );
    }

    dcel( size_t vertexCount, size_t faceCount = 0 );
    ~dcel();

    dcel& operator=( BOOST_COPY_ASSIGN_REF( dcel ) other ) {
        dcel temp( other );
        swap( temp );
        return *this;
    }

    dcel& operator=( BOOST_RV_REF( dcel ) other ) {
        swap( other );
        return *this;
    }

    void swap( dcel& other );
    void initialize( size_t vertexCount, size_t faceCount = 0 );
    void clear();

    size_t vertex_count() const;
    size_t face_count() const;
    size_t halfedge_count() const;
    size_t edge_count() const { return halfedge_count() / 2; }
    size_t boundary_count() const;

    bool is_triangle_mesh() const;

    bool is_boundary_vertex( size_t vertexId ) const;
    bool has_vertex_halfedge( size_t vertexId ) const;
    const_halfedge_handle get_vertex_halfedge( size_t vertexId ) const;
    halfedge_handle get_vertex_halfedge( size_t vertexId );
    void set_vertex_halfedge( size_t vertexId, const_halfedge_handle handle );
    void insert_vertex_halfedge( halfedge_handle handle );
    void adjust_vertex_halfedge( size_t vertexId );

    bool has_face_halfedge( size_t faceId ) const;
    const_halfedge_handle get_face_halfedge( size_t faceId ) const;
    halfedge_handle get_face_halfedge( size_t faceId );
    void set_face_halfedge( size_t faceId, const_halfedge_handle handle );
    void adjust_face_halfedge( size_t faceId );

    bool has_boundary_halfedge( size_t boundaryId ) const;
    const_halfedge_handle get_boundary_halfedge( size_t boundaryId ) const;
    halfedge_handle get_boundary_halfedge( size_t boundaryId );
    void set_boundary_halfedge( size_t boundaryId, const_halfedge_handle handle );
    void insert_boundary_halfedge( halfedge_handle handle );

    bool is_halfedge_valid( index_t halfedgeId ) const {
        return halfedgeId < m_halfedges.size() && m_halfedges[halfedgeId].m_vertex != INVALID_VERTEX_INDEX;
    }

    bool is_edge_valid( index_t edgeId ) const { return is_halfedge_valid( edgeId * 2 ); }

    const_halfedge_handle get_halfedge( size_t halfedgeId ) const;
    halfedge_handle get_halfedge( size_t halfedgeId );
    const_halfedge_handle get_edge_halfedge( size_t edgeId ) const { return get_halfedge( edgeId * 2 ); }
    halfedge_handle get_edge_halfedge( size_t edgeId ) { return get_halfedge( edgeId * 2 ); }

    const_halfedge_handle find_vertex_halfedge( size_t sourceVertexId, size_t destVertexId ) const;
    halfedge_handle find_vertex_halfedge( size_t sourceVertexId, size_t destVertexId );

    const_halfedge_handle find_face_halfedge( size_t face, size_t oppositeFace ) const;
    halfedge_handle find_face_halfedge( size_t face, size_t oppositeFace );

    std::pair<halfedge_handle, halfedge_handle> insert_edge();

    size_t add_face3( size_t v0, size_t v1, size_t v2 );
    void add_face3( size_t faceId, size_t v0, size_t v1, size_t v2 );
    size_t add_face4( size_t v0, size_t v1, size_t v2, size_t v3 );
    void add_face4( size_t faceId, size_t v0, size_t v1, size_t v2, size_t v3 );

    template <class InputIterator>
    size_t add_face( InputIterator begin, InputIterator end ) {
        index_t newFaceId = insert_new_face_id();
        add_face( newFaceId, begin, end );
        return newFaceId;
    }

    template <class InputIterator>
    void add_face( size_t faceId, InputIterator begin, InputIterator end ) {
        begin_face( faceId );
        for( InputIterator it = begin; it != end; ++it ) {
            add_face_vertex( *it );
        }
        end_face();
    }

    // face construction interface
    void begin_face( size_t faceId );
    size_t begin_face();
    void add_face_vertex( size_t vertexId );
    void end_face();

    /**
     * Check if the mesh will retain its manifold topology if the given edge is
     * collapsed
     */
    bool check_manifold_collapse( const_halfedge_handle handle ) const;

    /**
     * Collapses the `source` of the specified edge into the `target`
     *
     * @remark Note that the direction of the collapse is different depending
     *   on which side of the halfedge you are using.
     *
     * @param handle the edge to collapse
     * @param checkManifold (optional) if true, will check if the collapse is manifold before executing it.
     *   This option is provided since the check may have already been done separately.
     */
    void collapse_edge( halfedge_handle handle, bool checkManifold = true );
    bool try_collapse_edge( halfedge_handle handle );

    /**
     * Performs an edge-flip operation on the two triangles incident to this edge
     *
     * @remark The edge must be adjacent to two triangular faces
     * @remark The edge to be flipped must not be a boundary edge
     *
     * @param edge handle to the edge to be flipped
     * @throw std::runtime_error if the handle is to a boundary edge, or is not adjacent to two triangles
     */
    void flip_edge( halfedge_handle edge );

    /**
     * Finds all boundary loops in the surface, and caches them
     * @param force forces the caching (this may be neeccessary if external changes were
     *   made to the dcel without updating the boundary information)
     */
    void cache_boundary_info( bool force = false );

    bool boundary_info_cached() const { return m_isBoundaryInfoCached; }

    void set_boundary_info_cached( bool setting ) { m_isBoundaryInfoCached = setting; }

    bool is_boundary_face( index_t faceId ) const { return faceId > face_count(); }

    static index_t face_boundary_id( index_t faceId ) { return ( INVALID_FACE_INDEX - faceId ) - 1; }

    static index_t boundary_face_id( index_t boundaryId ) { return ( INVALID_FACE_INDEX - boundaryId ) - 1; }

    // Raw halfedge interface methods. This allows clients to perform operations on the structure
    // without the potential overhead of the halfedge_handle/const_halfedge_handle
  public:
    index_t get_face_next( index_t halfedgeId ) const { return m_halfedges[halfedgeId].m_next; }

    void set_face_next( index_t halfedgeId, index_t nextId ) { m_halfedges[halfedgeId].m_next = nextId; }

    // currently a non-trivial operation
    index_t get_face_prev( index_t halfedgeId ) const;

    index_t get_vertex_next( index_t halfedgeId ) const { return get_twin( get_face_next( halfedgeId ) ); }

    index_t get_vertex_prev( index_t halfedgeId ) const;

    index_t get_twin( index_t halfedgeId ) const {
        if( halfedgeId % 2 == 0 ) {
            return halfedgeId + 1;
        } else {
            return halfedgeId - 1;
        }
    }

    // The new method of edge management ensures this property
    index_t get_edge( index_t halfedgeId ) const { return halfedgeId / 2; }

    index_t get_target_vertex( index_t halfedgeId ) const { return m_halfedges[halfedgeId].m_vertex; }

    void set_target_vertex( index_t halfedgeId, index_t targetId ) { m_halfedges[halfedgeId].m_vertex = targetId; }

    index_t get_source_vertex( index_t halfedgeId ) const { return get_target_vertex( get_twin( halfedgeId ) ); }

    index_t get_current_face( index_t halfedgeId ) const { return m_halfedges[halfedgeId].m_face; }

    void set_current_face( index_t halfedgeId, index_t faceId ) { m_halfedges[halfedgeId].m_face = faceId; }

    index_t get_opposite_face( index_t halfedgeId ) const { return get_current_face( get_twin( halfedgeId ) ); }

    bool is_boundary_edge_face( index_t halfedgeId ) const {
        return is_boundary_face( get_current_face( halfedgeId ) );
    }

    bool is_boundary_edge( index_t halfedgeId ) const {
        return is_boundary_face( get_current_face( halfedgeId ) ) ||
               is_boundary_face( get_opposite_face( halfedgeId ) );
    }

    index_t get_current_boundary( index_t halfedgeId ) const {
        return face_boundary_id( get_current_face( halfedgeId ) );
    }

    void set_current_boundary( index_t halfedgeId, index_t boundaryId ) {
        set_current_face( halfedgeId, boundary_face_id( boundaryId ) );
    }

  private:
    class private_face_construction_info;

    static private_face_construction_info* create_private_face_construction_info();

    bool collapse_edge_internal( halfedge_handle handle, bool checkManifold );

    index_t insert_new_face_id();

    void add_edge_internal( index_t faceId, index_t startVertex, index_t endVertex );
    halfedge_handle create_halfedge_internal( index_t faceId, index_t startVertex, index_t endVertex,
                                              bool& alreadyExisted );

    std::vector<index_t> m_vertexHalfedges;
    std::vector<index_t> m_faceHalfedges;
    std::vector<index_t> m_boundaryHalfedges;
    std::vector<halfedge> m_halfedges;

    bool m_isAddingFace;
    bool m_isTriangleMesh;
    bool m_isBoundaryInfoCached;

    private_face_construction_info* m_faceConstructionInfo;
};

inline dcel::const_halfedge_handle::const_halfedge_handle( const dcel::halfedge_handle& other )
    : m_owner( other.get_owner() )
    , m_index( other.get_index() ) {}

inline dcel::const_halfedge_handle& dcel::const_halfedge_handle::operator=( const dcel::halfedge_handle& other ) {
    m_owner = other.get_owner();
    m_index = other.get_index();
    return *this;
}

inline bool dcel::const_halfedge_handle::operator==( const dcel::halfedge_handle& other ) const {
    return m_owner == other.get_owner() && m_index == other.get_index();
}

} // namespace geometry
} // namespace frantic

namespace std {

inline void swap( frantic::geometry::dcel& lhs, frantic::geometry::dcel& rhs ) { lhs.swap( rhs ); }

} // namespace std
