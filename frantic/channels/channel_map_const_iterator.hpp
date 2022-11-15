// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>

namespace frantic {
namespace channels {

class channel_map_const_iterator;

channel_map_const_iterator begin( const channel_map& map );
channel_map_const_iterator end( const channel_map& map );

/**
 * An iterator for visiting each channel in a channel_map. This is external from channel_map in order to avoid changing
 * channel_map which is part of the external DLL interface of various plugins we are shipping. That means we cannot
 * change the channel_map class at all.
 *
 * Use begin() and end() free functions on a channel_map to get the begin/end iterators.
 */
class channel_map_const_iterator {
  public:
    channel_map_const_iterator( const channel_map& map, std::size_t index = 0u );

    channel_map_const_iterator& operator++();
    channel_map_const_iterator operator++( int ) const;
    channel_map_const_iterator& operator+=( int step );

    channel_map_const_iterator& operator--();
    channel_map_const_iterator operator--( int ) const;
    channel_map_const_iterator& operator-=( int step );

    const channel& operator*() const;
    const channel& operator[]( std::size_t i ) const;
    const channel* operator->() const;

    bool operator==( const channel_map_const_iterator& rhs ) const;

  private:
    const channel_map* m_map;
    std::size_t m_chIndex;
};

inline channel_map_const_iterator operator+( channel_map_const_iterator lhs, int step ) {
    lhs += step;
    return lhs;
}

inline channel_map_const_iterator operator-( channel_map_const_iterator lhs, int step ) {
    lhs -= step;
    return lhs;
}

inline bool operator!=( const channel_map_const_iterator& lhs, const channel_map_const_iterator& rhs ) {
    return !lhs.operator==( rhs );
}

inline channel_map_const_iterator::channel_map_const_iterator( const channel_map& map, std::size_t index )
    : m_map( &map )
    , m_chIndex( index ) {}

inline channel_map_const_iterator& channel_map_const_iterator::operator++() {
    ++m_chIndex;
    return *this;
}

inline channel_map_const_iterator channel_map_const_iterator::operator++( int ) const {
    return channel_map_const_iterator( *m_map, m_chIndex + 1 );
}

inline channel_map_const_iterator& channel_map_const_iterator::operator+=( int step ) {
    m_chIndex += step;
    return *this;
}

inline channel_map_const_iterator& channel_map_const_iterator::operator--() {
    --m_chIndex;
    return *this;
}

inline channel_map_const_iterator channel_map_const_iterator::operator--( int ) const {
    return channel_map_const_iterator( *m_map, m_chIndex - 1 );
}

inline channel_map_const_iterator& channel_map_const_iterator::operator-=( int step ) {
    m_chIndex -= step;
    return *this;
}

inline const channel& channel_map_const_iterator::operator*() const { return ( *m_map )[m_chIndex]; }

inline const channel& channel_map_const_iterator::operator[]( std::size_t i ) const {
    return ( *m_map )[m_chIndex + i];
}

inline const channel* channel_map_const_iterator::operator->() const { return &( *m_map )[m_chIndex]; }

inline bool channel_map_const_iterator::operator==( const channel_map_const_iterator& rhs ) const {
    return m_map == rhs.m_map && m_chIndex == rhs.m_chIndex;
}

inline channel_map_const_iterator begin( const channel_map& map ) { return channel_map_const_iterator( map, 0 ); }

inline channel_map_const_iterator end( const channel_map& map ) {
    return channel_map_const_iterator( map, map.channel_count() );
}

} // namespace channels
} // namespace frantic
