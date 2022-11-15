// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>

#include <frantic/misc/exception_stream.hpp>

#include <utility>
#include <vector>

namespace frantic {
namespace channels {

class channel_column_map {

  private:
    channel_map m_channelMap;
    std::vector<std::pair<size_t, size_t>> m_columnMap;
    std::vector<bool> m_filledColumns;
    bool m_isEmpty;

  public:
    channel_column_map()
        : m_isEmpty( true ) {}

    channel_column_map( const channel_map& channelMap, size_t columnCount )
        : m_channelMap( channelMap )
        , m_columnMap( columnCount )
        , m_filledColumns( columnCount, false )
        , m_isEmpty( true ) {}

    channel_column_map( const channel_column_map& other )
        : m_channelMap( other.m_channelMap )
        , m_columnMap( other.m_columnMap )
        , m_filledColumns( other.m_filledColumns )
        , m_isEmpty( other.m_isEmpty ) {}

    ~channel_column_map() { m_columnMap.clear(); }

    channel_column_map& operator=( const channel_column_map& other ) {
        m_channelMap = other.m_channelMap;
        m_columnMap = other.m_columnMap;
        m_filledColumns = other.m_filledColumns;
        m_isEmpty = other.m_isEmpty;
        return *this;
    }

    size_t column_count() const { return m_columnMap.size(); }

    void remove_trailing_columns();

    void reset( const channel_map& channelMap, size_t columnCount ) {
        m_channelMap = channelMap;
        reset( columnCount );
    }

    void reset( size_t columnCount ) {
        m_columnMap.resize( columnCount );
        m_filledColumns.resize( columnCount );
        m_filledColumns.assign( columnCount, false );
        m_isEmpty = true;
    }

    void reset() {
        m_channelMap.reset();
        m_columnMap.clear();
        m_filledColumns.clear();
        m_isEmpty = true;
    }

    void resize( std::size_t newSize ) {
        m_columnMap.resize( newSize );
        std::size_t originalSize = m_filledColumns.size();
        m_filledColumns.resize( newSize );
        if( newSize > originalSize ) {
            for( std::size_t i = originalSize; i < m_filledColumns.size(); ++i ) {
                m_filledColumns[i] = false;
            }
        }
    }

    const channel_map& get_channel_map() const { return m_channelMap; }

    size_t find_mapping( const frantic::tstring& channelName, size_t channelElement ) const;
    size_t find_mapping( size_t channelIndex, size_t channelElement ) const;

    void set_column_mapping( size_t column, const frantic::tstring& channelName, size_t channelElement = 0 );
    void set_column_mapping( size_t column, size_t channelIndex, size_t channelElement = 0 );

    bool has_column_mapping( size_t column ) const { return m_filledColumns[column]; }

    std::pair<size_t, size_t> column_mapping( size_t column ) const {
        if( !has_column_mapping( column ) ) {
            throw frantic::exception_stream()
                << "channel_column_map::get_column_mapping -- No mapping for column " << column << ".";
        }

        return m_columnMap[column];
    }

    bool is_empty() const { return m_isEmpty; }

    void copy_structure_from_strings( std::vector<char>& buffer, const std::vector<std::string>& srcData ) const;
    void copy_structure_from_strings( char* buffer, const std::vector<std::string>& srcData ) const;

    /**
     * Writes out the given column mapping as a set of csv headers
     *
     * \param columnMap the column mapping to dump
     * \param columns output array to store the dump
     */
    void dump_as_headers( std::vector<std::string>& columns );
};

} // namespace channels
} // namespace frantic
