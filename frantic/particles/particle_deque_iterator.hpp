// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frantic {
namespace particles {

class particle_deque;

class particle_deque_iterator : std::iterator<std::random_access_iterator_tag, char*, std::ptrdiff_t, char*, char*> {
  private:
    friend class particle_deque;
    friend class particle_deque_const_iterator;

    particle_deque* m_container;
    std::size_t m_offset;

    particle_deque_iterator( particle_deque& owner, std::size_t offset )
        : m_container( &owner )
        , m_offset( offset ) {}

  public:
    particle_deque_iterator()
        : m_container( NULL )
        , m_offset( 0 ) {}

    inline std::size_t structure_size();                 // Implemented in particle_deque.hpp
    inline value_type operator*();                       // Implemented in particle_deque.hpp
    inline value_type operator[]( const std::size_t i ); // Implemented in particle_deque.hpp

    particle_deque_iterator& operator--() {
        --m_offset;
        return *this;
    }
    particle_deque_iterator& operator++() {
        ++m_offset;
        return *this;
    }
    particle_deque_iterator& operator--( const int i ) {
        m_offset -= i;
        return *this;
    }
    particle_deque_iterator& operator++( const int i ) {
        m_offset += i;
        return *this;
    }
    particle_deque_iterator& operator+=( const difference_type i ) {
        m_offset += i;
        return *this;
    }
    particle_deque_iterator& operator-=( const difference_type i ) {
        m_offset -= i;
        return *this;
    }

    friend difference_type operator-( const particle_deque_iterator& lhs, const particle_deque_iterator& rhs );
    friend particle_deque_iterator operator+( const particle_deque_iterator& lhs, const difference_type rhs );
    friend particle_deque_iterator operator-( const particle_deque_iterator& lhs, const difference_type rhs );

    bool operator<( const particle_deque_iterator& rhs ) const { return m_offset < rhs.m_offset; }
    bool operator>( const particle_deque_iterator& rhs ) const { return m_offset > rhs.m_offset; }
    bool operator<=( const particle_deque_iterator& rhs ) const { return m_offset <= rhs.m_offset; }
    bool operator>=( const particle_deque_iterator& rhs ) const { return m_offset >= rhs.m_offset; }
    bool operator!=( const particle_deque_iterator& rhs ) const { return m_offset != rhs.m_offset; }
    bool operator==( const particle_deque_iterator& rhs ) const { return m_offset == rhs.m_offset; }
};

inline particle_deque_iterator operator+( const particle_deque_iterator& lhs,
                                          const particle_deque_iterator::difference_type rhs ) {
    return particle_deque_iterator( *lhs.m_container, lhs.m_offset + rhs );
}

inline particle_deque_iterator operator-( const particle_deque_iterator& lhs,
                                          const particle_deque_iterator::difference_type rhs ) {
    return particle_deque_iterator( *lhs.m_container, lhs.m_offset - rhs );
}

inline particle_deque_iterator::difference_type operator-( const particle_deque_iterator& lhs,
                                                           const particle_deque_iterator& rhs ) {
    return lhs.m_offset - rhs.m_offset;
}

class particle_deque_const_iterator
    : std::iterator<std::random_access_iterator_tag, const char*, std::ptrdiff_t, const char*, const char*> {
  private:
    friend class particle_deque;

    const particle_deque* m_container;
    std::size_t m_offset;

    particle_deque_const_iterator( const particle_deque& owner, std::size_t offset )
        : m_container( &owner )
        , m_offset( offset ) {}

  public:
    particle_deque_const_iterator()
        : m_container( NULL )
        , m_offset( 0 ) {}
    particle_deque_const_iterator( const particle_deque_iterator& it )
        : m_container( it.m_container )
        , m_offset( it.m_offset ) {}

    inline std::size_t structure_size();                 // Implemented in particle_deque.hpp
    inline value_type operator*();                       // Implemented in particle_deque.hpp
    inline value_type operator[]( const std::size_t i ); // Implemented in particle_deque.hpp

    particle_deque_const_iterator& operator--() {
        --m_offset;
        return *this;
    }
    particle_deque_const_iterator& operator++() {
        ++m_offset;
        return *this;
    }
    particle_deque_const_iterator& operator--( const int i ) {
        m_offset -= i;
        return *this;
    }
    particle_deque_const_iterator& operator++( const int i ) {
        m_offset += i;
        return *this;
    }
    particle_deque_const_iterator& operator+=( const difference_type i ) {
        m_offset += i;
        return *this;
    }
    particle_deque_const_iterator& operator-=( const difference_type i ) {
        m_offset -= i;
        return *this;
    }

    friend difference_type operator-( const particle_deque_const_iterator& lhs,
                                      const particle_deque_const_iterator& rhs );
    friend particle_deque_const_iterator operator+( const particle_deque_const_iterator& lhs,
                                                    const difference_type rhs );
    friend particle_deque_const_iterator operator-( const particle_deque_const_iterator& lhs,
                                                    const difference_type rhs );

    bool operator<( const particle_deque_const_iterator& rhs ) const { return m_offset < rhs.m_offset; }
    bool operator>( const particle_deque_const_iterator& rhs ) const { return m_offset > rhs.m_offset; }
    bool operator<=( const particle_deque_const_iterator& rhs ) const { return m_offset <= rhs.m_offset; }
    bool operator>=( const particle_deque_const_iterator& rhs ) const { return m_offset >= rhs.m_offset; }
    bool operator!=( const particle_deque_const_iterator& rhs ) const { return m_offset != rhs.m_offset; }
    bool operator==( const particle_deque_const_iterator& rhs ) const { return m_offset == rhs.m_offset; }
};

inline particle_deque_const_iterator operator+( const particle_deque_const_iterator& lhs,
                                                const particle_deque_const_iterator::difference_type rhs ) {
    return particle_deque_const_iterator( *lhs.m_container, lhs.m_offset + rhs );
}

inline particle_deque_const_iterator operator-( const particle_deque_const_iterator& lhs,
                                                const particle_deque_const_iterator::difference_type rhs ) {
    return particle_deque_const_iterator( *lhs.m_container, lhs.m_offset - rhs );
}

inline particle_deque_const_iterator::difference_type operator-( const particle_deque_const_iterator& lhs,
                                                                 const particle_deque_const_iterator& rhs ) {
    return lhs.m_offset - rhs.m_offset;
}

} // namespace particles
} // namespace frantic
