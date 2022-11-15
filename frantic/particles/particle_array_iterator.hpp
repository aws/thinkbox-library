// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iostream>
#include <iterator>

namespace frantic {
namespace particles {
class particle_array;
class const_particle_array_iterator;

class particle_array_iterator {
    char* m_current;
    size_t m_particleSize;

    friend class const_particle_array_iterator;

  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef std::ptrdiff_t difference_type;
    typedef void value_type;
    typedef void pointer;
    typedef void reference;

  public:
    particle_array_iterator( char* begin, std::size_t particleSize ) {
        m_current = begin;
        m_particleSize = particleSize;
    }

    particle_array_iterator()
        : m_current( NULL )
        , m_particleSize( 0 ) {}

    void dump( std::ostream& out ) {
        out << "Begin const_particle_grid_tree_iterator Dump:" << std::endl;
        out << "m_current: " << (void*)m_current << std::endl;
        out << "m_particleSize: " << m_particleSize << std::endl;
        out << "End const_particle_grid_tree_iterator Dump:" << std::endl;
    }

    std::size_t structure_size() const { return m_particleSize; }

    char* operator*() const { return m_current; }

    char* operator[]( difference_type i ) const { return m_current + i * m_particleSize; }

    // preincrement operator
    particle_array_iterator& operator++() {
        m_current += m_particleSize;
        return *this;
    }

    // predecrement operator
    particle_array_iterator& operator--() {
        m_current -= m_particleSize;
        return *this;
    }

    particle_array_iterator operator++( int offset ) {
        char* current = m_current;

        m_current += offset * m_particleSize;

        return particle_array_iterator( current, m_particleSize );
    }

    particle_array_iterator operator+( std::ptrdiff_t offset ) const {
        return particle_array_iterator( m_current + offset * m_particleSize, m_particleSize );
    }

    particle_array_iterator operator-( std::ptrdiff_t offset ) const {
        return particle_array_iterator( m_current - offset * m_particleSize, m_particleSize );
    }

    std::ptrdiff_t operator-( const particle_array_iterator& rhs ) const {
        return ( m_current - rhs.m_current ) / m_particleSize;
    }

    bool operator==( const particle_array_iterator& rhs ) const { return m_current == rhs.m_current; }

    bool operator!=( const particle_array_iterator& rhs ) const { return m_current != rhs.m_current; }

    bool operator<( const particle_array_iterator& rhs ) const { return m_current < rhs.m_current; }

    bool operator<=( const particle_array_iterator& rhs ) const { return m_current <= rhs.m_current; }

    bool operator>( const particle_array_iterator& rhs ) const { return m_current > rhs.m_current; }

    bool operator>=( const particle_array_iterator& rhs ) const { return m_current >= rhs.m_current; }

    particle_array_iterator& operator+=( ptrdiff_t offs ) {
        m_current += m_particleSize * offs;
        return *this;
    }
};

class const_particle_array_iterator {
    const char* m_current;
    size_t m_particleSize;

  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef std::ptrdiff_t difference_type;
    typedef void value_type;
    typedef void pointer;
    typedef void reference;

    const_particle_array_iterator()
        : m_current( NULL )
        , m_particleSize( 0 ) {}

    const_particle_array_iterator( const char* begin, std::size_t particleSize ) {
        m_current = begin;
        m_particleSize = particleSize;
    }

    const_particle_array_iterator( const particle_array_iterator& rhs )
        : m_current( rhs.m_current )
        , m_particleSize( rhs.m_particleSize ) {}

    const_particle_array_iterator& operator=( const particle_array_iterator& rhs ) {
        m_current = rhs.m_current;
        m_particleSize = rhs.m_particleSize;
        return *this;
    }

    void dump( std::ostream& out ) {
        out << "Begin const_particle_grid_tree_iterator Dump:" << std::endl;
        out << "m_current: " << (void*)m_current << std::endl;
        out << "m_particleSize: " << m_particleSize << std::endl;
        out << "End const_particle_grid_tree_iterator Dump:" << std::endl;
    }

    std::size_t structure_size() { return m_particleSize; }

    const char* operator*() const { return m_current; }

    // preincrement operator
    const_particle_array_iterator& operator++() {
        m_current += m_particleSize;
        return *this;
    }

    // predecrement operator
    const_particle_array_iterator& operator--() {
        m_current -= m_particleSize;
        return *this;
    }

    // TODO: postincrement operator
    // particle_grid_tree_iterator operator++(int);

    const_particle_array_iterator operator+( std::ptrdiff_t offset ) const {
        return const_particle_array_iterator( m_current + offset * m_particleSize, m_particleSize );
    }

    const_particle_array_iterator operator-( std::ptrdiff_t offset ) const {
        return const_particle_array_iterator( m_current - offset * m_particleSize, m_particleSize );
    }

    std::ptrdiff_t operator-( const const_particle_array_iterator& rhs ) const {
        return ( m_current - rhs.m_current ) / m_particleSize;
    }

    bool operator==( const const_particle_array_iterator& rhs ) const { return m_current == rhs.m_current; }

    bool operator!=( const const_particle_array_iterator& rhs ) const { return m_current != rhs.m_current; }

    bool operator<( const const_particle_array_iterator& rhs ) const { return m_current < rhs.m_current; }

    bool operator<=( const const_particle_array_iterator& rhs ) const { return m_current <= rhs.m_current; }

    bool operator>( const const_particle_array_iterator& rhs ) const { return m_current > rhs.m_current; }

    bool operator>=( const const_particle_array_iterator& rhs ) const { return m_current >= rhs.m_current; }

    const_particle_array_iterator& operator+=( ptrdiff_t offs ) {
        m_current += m_particleSize * offs;
        return *this;
    }
};
} // namespace particles
} // namespace frantic
