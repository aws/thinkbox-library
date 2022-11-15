// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/dcel.hpp>

#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>

#include <exception>
#include <iterator>
#include <utility>

namespace frantic {
namespace geometry {

template <class Halfedge_handle_t>
class dcel_halfedge_iterator_base {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef Halfedge_handle_t value_type;
    typedef std::ptrdiff_t difference_type;
    typedef Halfedge_handle_t* pointer;
    typedef Halfedge_handle_t reference;
    enum tag_end {
        END,
    };

  public:
    explicit dcel_halfedge_iterator_base()
        : m_isBegin( false ) {}

    explicit dcel_halfedge_iterator_base( const Halfedge_handle_t& start )
        : m_current( start )
        , m_isBegin( start != dcel::INVALID_HALFEDGE_HANDLE ) {}

    explicit dcel_halfedge_iterator_base( const Halfedge_handle_t& start, tag_end )
        : m_current( start )
        , m_isBegin( false ) {}

    dcel_halfedge_iterator_base& operator=( const dcel_halfedge_iterator_base& other ) {
        m_current = other.m_current;
        m_isBegin = other.m_isBegin;
        return *this;
    }

    bool operator==( const dcel_halfedge_iterator_base& other ) const {
        return m_current == other.m_current && ( m_isBegin == other.m_isBegin );
    }

    bool operator!=( const dcel_halfedge_iterator_base& other ) const { return !( *this == other ); }

    const Halfedge_handle_t operator*() const { return get(); }

    const Halfedge_handle_t* operator->() const { return &get(); }

    const Halfedge_handle_t get() const { return m_current; }

    bool is_start_iterator() const { return m_isBegin; }

  protected:
    Halfedge_handle_t m_current;
    bool m_isBegin;
};

// iteration rotation direction is the same as the face winding of the underlying dcel
// for most sane applications, that means counter-clockwise
// TODO: might it make sense to allow reverse iterators as well?
// it would not be at all difficult to support, even as a
// parameter on the iterators themselves
template <class Halfedge_handle_t>
class dcel_vertex_iterator_template : public dcel_halfedge_iterator_base<Halfedge_handle_t> {
  public:
    explicit dcel_vertex_iterator_template()
        : dcel_halfedge_iterator_base<Halfedge_handle_t>() {}
    explicit dcel_vertex_iterator_template( const Halfedge_handle_t& start )
        : dcel_halfedge_iterator_base<Halfedge_handle_t>( start ) {}
    explicit dcel_vertex_iterator_template( const Halfedge_handle_t& start,
                                            typename dcel_halfedge_iterator_base<Halfedge_handle_t>::tag_end )
        : dcel_halfedge_iterator_base<Halfedge_handle_t>( start, dcel_halfedge_iterator_base<Halfedge_handle_t>::END ) {
    }

    dcel_vertex_iterator_template& operator++() {
        if( this->m_isBegin ) {
            this->m_isBegin = false;
        }
        this->m_current = this->m_current.vertex_next();
        return *this;
    }

    dcel_vertex_iterator_template operator++( int ) {
        dcel_vertex_iterator_template copy = *this;
        ++( *this );
        return copy;
    }
};

template <class Halfedge_handle_t>
class dcel_face_iterator_template : public dcel_halfedge_iterator_base<Halfedge_handle_t> {
  public:
    explicit dcel_face_iterator_template()
        : dcel_halfedge_iterator_base<Halfedge_handle_t>() {}
    explicit dcel_face_iterator_template( const Halfedge_handle_t& start )
        : dcel_halfedge_iterator_base<Halfedge_handle_t>( start ) {}
    explicit dcel_face_iterator_template( const Halfedge_handle_t& start,
                                          typename dcel_halfedge_iterator_base<Halfedge_handle_t>::tag_end )
        : dcel_halfedge_iterator_base<Halfedge_handle_t>( start, dcel_halfedge_iterator_base<Halfedge_handle_t>::END ) {
    }

    dcel_face_iterator_template& operator++() {
        if( this->m_isBegin ) {
            this->m_isBegin = false;
        }
        this->m_current = this->m_current.face_next();
        return *this;
    }

    dcel_face_iterator_template operator++( int ) {
        dcel_face_iterator_template copy = *this;
        ++( *this );
        return copy;
    }
};

typedef dcel_vertex_iterator_template<dcel::halfedge_handle> dcel_vertex_iterator;
typedef dcel_vertex_iterator_template<dcel::const_halfedge_handle> const_dcel_vertex_iterator;
typedef dcel_face_iterator_template<dcel::halfedge_handle> dcel_face_iterator;
typedef dcel_face_iterator_template<dcel::const_halfedge_handle> const_dcel_face_iterator;

typedef boost::iterator_range<dcel_vertex_iterator> dcel_vertex_range_t;
typedef boost::iterator_range<const_dcel_vertex_iterator> const_dcel_vertex_range_t;
typedef boost::iterator_range<dcel_face_iterator> dcel_face_range_t;
typedef boost::iterator_range<const_dcel_face_iterator> const_dcel_face_range_t;

template <class T>
struct halfedge_iterator_traits;

template <>
struct halfedge_iterator_traits<dcel::halfedge_handle> {
    typedef dcel_vertex_iterator vertex_iterator_t;
    typedef dcel_vertex_range_t vertex_range_t;
    typedef dcel_face_iterator face_iterator_t;
    typedef dcel_face_range_t face_range_t;
};

template <>
struct halfedge_iterator_traits<dcel::const_halfedge_handle> {
    typedef const_dcel_vertex_iterator vertex_iterator_t;
    typedef const_dcel_vertex_range_t vertex_range_t;
    typedef const_dcel_face_iterator face_iterator_t;
    typedef const_dcel_face_range_t face_range_t;
};

namespace detail {
template <class Halfedge_t>
inline void ensure_same_vertex( const Halfedge_t& start, const Halfedge_t& end ) {
    if( start.target_vertex() != end.target_vertex() ) {
        throw std::runtime_error(
            "ensure_same_vertex : the halfedges provided for the vertex range did not match the same vertex" );
    }
}

template <class Halfedge_t>
inline void ensure_same_face( const Halfedge_t& start, const Halfedge_t& end ) {
    if( start.current_face() != end.current_face() ) {
        throw std::runtime_error(
            "ensure_same_face : the halfedges provided for the face range did not match the same face" );
    }
}
} // namespace detail

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_range_t dcel_vertex_range( const Halfedge_t& start,
                                                                                        const Halfedge_t& end ) {
    detail::ensure_same_vertex( start, end );
    return boost::make_iterator_range( dcel_vertex_begin( start ), dcel_vertex_end( end ) );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t dcel_vertex_begin( const Halfedge_t& start ) {
    return typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t( start );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t dcel_vertex_end( const Halfedge_t& end ) {
    return typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t(
        end, halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t::END );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_range_t
dcel_vertex_inclusive_range( const Halfedge_t& start, const Halfedge_t& end ) {
    detail::ensure_same_vertex( start, end );
    return boost::make_iterator_range( dcel_vertex_inclusive_begin( start ), dcel_vertex_inclusive_end( end ) );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t
dcel_vertex_inclusive_begin( const Halfedge_t& start ) {
    return dcel_vertex_begin( start );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t
dcel_vertex_inclusive_end( const Halfedge_t& end ) {
    return dcel_vertex_end( end.vertex_next() );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_range_t
dcel_vertex_cycle_range( const Halfedge_t& start ) {
    return boost::make_iterator_range( dcel_vertex_cycle_begin( start ), dcel_vertex_cycle_end( start ) );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t
dcel_vertex_cycle_begin( const Halfedge_t& start ) {
    return dcel_vertex_begin( start );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::vertex_iterator_t
dcel_vertex_cycle_end( const Halfedge_t& start ) {
    return dcel_vertex_end( start );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_range_t dcel_face_range( const Halfedge_t& start,
                                                                                    const Halfedge_t& end ) {
    detail::ensure_same_face( start, end );
    return boost::make_iterator_range( dcel_face_begin( start ), dcel_face_end( end ) );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t dcel_face_begin( const Halfedge_t& start ) {
    return typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t( start );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t dcel_face_end( const Halfedge_t& end ) {
    return typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t(
        end, halfedge_iterator_traits<Halfedge_t>::face_iterator_t::END );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_range_t dcel_face_inclusive_range( const Halfedge_t& start,
                                                                                              const Halfedge_t& end ) {
    detail::ensure_same_face( start, end );
    return boost::make_iterator_range( dcel_face_inclusive_begin( start ), dcel_face_inclusive_end( end ) );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t
dcel_face_inclusive_begin( const Halfedge_t& start ) {
    return dcel_face_begin( start );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t dcel_face_inclusive_end( const Halfedge_t& end ) {
    return dcel_face_end( end.face_next() );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_range_t dcel_face_cycle_range( const Halfedge_t& start ) {
    return boost::make_iterator_range( dcel_face_cycle_begin( start ), dcel_face_cycle_end( start ) );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t dcel_face_cycle_begin( const Halfedge_t& start ) {
    return dcel_face_begin( start );
}

template <class Halfedge_t>
inline typename halfedge_iterator_traits<Halfedge_t>::face_iterator_t dcel_face_cycle_end( const Halfedge_t& end ) {
    return dcel_face_end( end );
}

class dcel_vertex_adjacency_iterator {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef dcel::index_t value_type;
    typedef ptrdiff_t difference_type;
    typedef dcel::index_t* pointer;
    typedef dcel::index_t reference;

  public:
    dcel_vertex_adjacency_iterator() {}
    explicit dcel_vertex_adjacency_iterator( const const_dcel_vertex_iterator& it )
        : m_impl( it ) {}
    dcel_vertex_adjacency_iterator( const dcel_vertex_adjacency_iterator& other )
        : m_impl( other.m_impl ) {}

    bool operator==( const dcel_vertex_adjacency_iterator& other ) { return m_impl == other.m_impl; }

    bool operator!=( const dcel_vertex_adjacency_iterator& other ) { return m_impl != other.m_impl; }

    dcel_vertex_adjacency_iterator& operator++() {
        ++m_impl;
        return *this;
    }

    dcel_vertex_adjacency_iterator operator++( int ) {
        dcel_vertex_adjacency_iterator copy = *this;
        ++m_impl;
        return copy;
    }

    dcel::index_t operator*() const { return ( *m_impl ).source_vertex(); }

  private:
    const_dcel_vertex_iterator m_impl;
};

class dcel_face_adjacency_iterator {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef dcel::index_t value_type;
    typedef ptrdiff_t difference_type;
    typedef dcel::index_t* pointer;
    typedef dcel::index_t reference;

    dcel_face_adjacency_iterator() {}
    dcel_face_adjacency_iterator( const const_dcel_face_iterator& it )
        : m_impl( it ) {}
    dcel_face_adjacency_iterator( const dcel_face_adjacency_iterator& other )
        : m_impl( other.m_impl ) {}

    bool operator==( const dcel_face_adjacency_iterator& other ) { return m_impl == other.m_impl; }

    bool operator!=( const dcel_face_adjacency_iterator& other ) { return m_impl != other.m_impl; }

    dcel_face_adjacency_iterator& operator++() {
        ++m_impl;
        return *this;
    }

    dcel_face_adjacency_iterator operator++( int ) {
        dcel_face_adjacency_iterator copy = *this;
        ++m_impl;
        return copy;
    }

    dcel::index_t operator*() { return ( *m_impl ).opposite_face(); }

  private:
    const_dcel_face_iterator m_impl;
};

typedef boost::iterator_range<dcel_vertex_adjacency_iterator> dcel_vertex_adjacency_range_t;
typedef boost::iterator_range<dcel_face_adjacency_iterator> dcel_face_adjacency_range_t;

inline dcel_vertex_adjacency_iterator dcel_vertex_adjacency_begin( const dcel::const_halfedge_handle& start ) {
    return dcel_vertex_adjacency_iterator( dcel_vertex_begin( start ) );
}

inline dcel_vertex_adjacency_iterator dcel_vertex_adjacency_end( const dcel::const_halfedge_handle& start ) {
    return dcel_vertex_adjacency_iterator( dcel_vertex_end( start ) );
}

inline dcel_vertex_adjacency_range_t dcel_vertex_adjacency_range( const dcel::const_halfedge_handle& start,
                                                                  const dcel::const_halfedge_handle& end ) {
    return dcel_vertex_adjacency_range_t( dcel_vertex_adjacency_begin( start ), dcel_vertex_adjacency_end( end ) );
}

inline dcel_vertex_adjacency_iterator
dcel_vertex_adjacency_inclusive_begin( const dcel::const_halfedge_handle& start ) {
    return dcel_vertex_adjacency_iterator( dcel_vertex_inclusive_begin( start ) );
}

inline dcel_vertex_adjacency_iterator dcel_vertex_adjacency_inclusive_end( const dcel::const_halfedge_handle& start ) {
    return dcel_vertex_adjacency_iterator( dcel_vertex_inclusive_end( start ) );
}

inline dcel_vertex_adjacency_range_t dcel_vertex_adjacency_inclusive_range( const dcel::const_halfedge_handle& start,
                                                                            const dcel::const_halfedge_handle& end ) {
    return dcel_vertex_adjacency_range_t( dcel_vertex_adjacency_inclusive_begin( start ),
                                          dcel_vertex_adjacency_inclusive_end( end ) );
}

inline dcel_vertex_adjacency_iterator dcel_vertex_adjacency_cycle_begin( const dcel::const_halfedge_handle& start ) {
    return dcel_vertex_adjacency_iterator( dcel_vertex_cycle_begin( start ) );
}

inline dcel_vertex_adjacency_iterator dcel_vertex_adjacency_cycle_end( const dcel::const_halfedge_handle& start ) {
    return dcel_vertex_adjacency_iterator( dcel_vertex_cycle_end( start ) );
}

inline dcel_vertex_adjacency_range_t dcel_vertex_adjacency_cycle_range( const dcel::const_halfedge_handle& start ) {
    return dcel_vertex_adjacency_range_t( dcel_vertex_adjacency_cycle_begin( start ),
                                          dcel_vertex_adjacency_cycle_end( start ) );
}

inline dcel_face_adjacency_iterator dcel_face_adjacency_begin( const dcel::const_halfedge_handle& start ) {
    return dcel_face_adjacency_iterator( dcel_face_begin( start ) );
}

inline dcel_face_adjacency_iterator dcel_face_adjacency_end( const dcel::const_halfedge_handle& start ) {
    return dcel_face_adjacency_iterator( dcel_face_end( start ) );
}

inline dcel_face_adjacency_range_t dcel_face_adjacency_range( const dcel::const_halfedge_handle& start,
                                                              const dcel::const_halfedge_handle& end ) {
    return dcel_face_adjacency_range_t( dcel_face_adjacency_begin( start ), dcel_face_adjacency_end( end ) );
}

inline dcel_face_adjacency_iterator dcel_face_adjacency_inclusive_begin( const dcel::const_halfedge_handle& start ) {
    return dcel_face_adjacency_iterator( dcel_face_inclusive_begin( start ) );
}

inline dcel_face_adjacency_iterator dcel_face_adjacency_inclusive_end( const dcel::const_halfedge_handle& start ) {
    return dcel_face_adjacency_iterator( dcel_face_inclusive_end( start ) );
}

inline dcel_face_adjacency_range_t dcel_face_adjacency_inclusive_range( const dcel::const_halfedge_handle& start,
                                                                        const dcel::const_halfedge_handle& end ) {
    return dcel_face_adjacency_range_t( dcel_face_adjacency_inclusive_begin( start ),
                                        dcel_face_adjacency_inclusive_end( end ) );
}

inline dcel_face_adjacency_iterator dcel_face_adjacency_cycle_begin( const dcel::const_halfedge_handle& start ) {
    return dcel_face_adjacency_iterator( dcel_face_cycle_begin( start ) );
}

inline dcel_face_adjacency_iterator dcel_face_adjacency_cycle_end( const dcel::const_halfedge_handle& start ) {
    return dcel_face_adjacency_iterator( dcel_face_cycle_end( start ) );
}

inline dcel_face_adjacency_range_t dcel_face_adjacency_cycle_range( const dcel::const_halfedge_handle& start ) {
    return dcel_face_adjacency_range_t( dcel_face_adjacency_cycle_begin( start ),
                                        dcel_face_adjacency_cycle_end( start ) );
}

} // namespace geometry
} // namespace frantic
