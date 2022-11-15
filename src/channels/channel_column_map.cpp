// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/channel_column_map.hpp>

namespace frantic {
namespace channels {

void channel_column_map::remove_trailing_columns() {
    while( column_count() > 0 && !m_filledColumns.back() ) {
        m_filledColumns.pop_back();
        m_columnMap.pop_back();
    }
}

size_t channel_column_map::find_mapping( const frantic::tstring& channelName, size_t channelElement ) const {
    if( !m_channelMap.has_channel( channelName ) ) {
        throw frantic::exception_stream() << "channel_column_map::set_column_mapping -- No such channel "
                                          << frantic::strings::to_string( channelName ) << ".";
    }

    return find_mapping( m_channelMap.channel_index( channelName ), channelElement );
}

size_t channel_column_map::find_mapping( size_t channelIndex, size_t channelElement ) const {
    for( size_t i = 0; i < m_columnMap.size(); ++i ) {
        if( m_filledColumns[i] && m_columnMap[i].first == channelIndex && m_columnMap[i].second == channelElement ) {
            return i;
        }
    }
    return static_cast<size_t>( -1 );
}

void channel_column_map::set_column_mapping( size_t column, const frantic::tstring& channelName,
                                             size_t channelElement ) {
    if( !m_channelMap.has_channel( channelName ) ) {
        throw frantic::exception_stream() << "channel_column_map::set_column_mapping -- No such channel "
                                          << frantic::strings::to_string( channelName ) << ".";
    }

    set_column_mapping( column, m_channelMap.channel_index( channelName ), channelElement );
}

void channel_column_map::set_column_mapping( size_t column, size_t channelIndex, size_t channelElement ) {

    if( channelIndex >= m_channelMap.channel_count() ) {
        throw frantic::exception_stream() << "channel_column_map::set_column_mapping -- Invalid channel "
                                          << channelIndex << " (max = " << m_channelMap.channel_count() << ").";
    } else if( channelElement >= m_channelMap[channelIndex].arity() ) {
        throw frantic::exception_stream()
            << "channel_column_map::set_column_mapping -- Invalid channel element " << channelElement << " in channel "
            << channelIndex << " (arity = " << m_channelMap[channelIndex].arity() << ").";
    } else if( column >= m_columnMap.size() ) {
        throw frantic::exception_stream() << "channel_column_map::set_column_mapping -- Invalid column " << column
                                          << " in column map (max = " << m_columnMap.size() << ").";
    } else if( column >= m_filledColumns.size() ) {
        throw frantic::exception_stream() << "channel_column_map::set_column_mapping -- Invalid column " << column
                                          << " in filled columns (max = " << m_filledColumns.size() << ").";
    }

    size_t existingMapping = find_mapping( channelIndex, channelElement );

    if( existingMapping < m_columnMap.size() && existingMapping != column ) {
        throw frantic::exception_stream()
            << "channel_column_map::set_column_mapping -- Element ("
            << frantic::strings::to_string( m_channelMap[channelIndex].name() ) << "," << channelElement
            << ") is already mapped from column " << existingMapping << ".";
    }

    m_columnMap[column] = std::make_pair( channelIndex, channelElement );
    m_filledColumns[column] = true;
    m_isEmpty = false;
}

void channel_column_map::copy_structure_from_strings( std::vector<char>& buffer,
                                                      const std::vector<std::string>& srcData ) const {
    copy_structure_from_strings( &buffer[0], srcData );
}

void channel_column_map::copy_structure_from_strings( char* buffer, const std::vector<std::string>& srcData ) const {

    if( m_isEmpty ) {
        m_channelMap.copy_structure_from_strings( buffer, srcData );
    } else if( srcData.size() < m_columnMap.size() ) {
        throw frantic::exception_stream()
            << "channel_column_map::copy_structure_from_strings -- Insufficient columns mapped to copy data.";
    } else if( m_channelMap.needs_scope_management() ) {
        throw frantic::exception_stream()
            << "channel_column_map::copy_structure_from_strings -- Scope management types not implemented yet.";
    } else {
        for( size_t i = 0; i < m_columnMap.size(); ++i ) {
            if( m_filledColumns[i] ) { // Skip columns that have not been assigned a mapping
                const size_t channelIndex = m_columnMap[i].first;
                const size_t channelElement = m_columnMap[i].second;
                // Parse the data in column i into channel mapped to column i at index mapped to column i
                parse_channel_value_from_string(
                    m_channelMap[channelIndex].data_type(), srcData[i],
                    m_channelMap[channelIndex].get_channel_data_pointer( buffer, (int)channelElement ) );
            }
        }
    }
}

void channel_column_map::dump_as_headers( std::vector<std::string>& columns ) {
    for( size_t i = 0; i < m_columnMap.size(); ++i ) {
        if( has_column_mapping( i ) ) {
            size_t channel = m_columnMap[i].first;
            size_t index = m_columnMap[i].second;

            frantic::tstring typestr = frantic::tstring( channel_data_type_str( m_channelMap[channel].data_type() ) ) +
                                       _T(" ") + m_channelMap[channel].name();

            if( m_channelMap[channel].arity() > 1 ) {
                typestr += _T("[") + boost::lexical_cast<frantic::tstring>( index ) + _T("]");
            }

            columns.push_back( frantic::strings::to_string( typestr ) );
        } else {
            // an empty column signifies a gap/unmapped data
            columns.push_back( "" );
        }
    }
}

} // namespace channels
} // namespace frantic
