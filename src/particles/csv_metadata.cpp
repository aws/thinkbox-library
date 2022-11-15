// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/csv_metadata.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/type.hpp>

namespace frantic {
namespace particles {
namespace csv {

namespace {

frantic::tstring serialize_single_column_map_to_metadata_format( const frantic::tstring& channel, size_t element ) {
    return channel + _T( "[" ) + boost::lexical_cast<frantic::tstring>( element ) + _T( "]" );
}

std::pair<frantic::tstring, size_t> parse_single_column_map_from_metadata_format( const frantic::tstring& metadata ) {
    size_t bracketPos = metadata.find_first_of( _T( "[" ) );

    // alternatively, throw an exception
    if( bracketPos == frantic::tstring::npos ) {
        return std::make_pair( metadata, 0 );
    }

    size_t closeBracketPos = metadata.find_first_of( _T( "]" ) );

    std::pair<frantic::tstring, size_t> result;
    result.first = metadata.substr( 0, bracketPos );
    boost::trim( result.first );
    frantic::tstring c = metadata.substr( bracketPos + 1, closeBracketPos - bracketPos - 1 );
    result.second = boost::lexical_cast<int>( c );
    return result;
}

} // namespace

frantic::tstring serialize_channel_map_to_metadata_format( const frantic::channels::channel_map& cMap ) {
    frantic::tstring channelRep;

    for( size_t i = 0; i < cMap.channel_count(); ++i ) {
        channelRep += channel_data_type_str( cMap[i].arity(), cMap[i].data_type() ) + _T( " " ) + cMap[i].name();
        if( i < cMap.channel_count() - 1 ) {
            channelRep += _T( "," );
        }
    }

    return channelRep;
}

frantic::tstring serialize_column_map_to_metadata_format( const frantic::channels::channel_column_map& columnMap ) {

    frantic::tstring result = _T( "" );

    const frantic::channels::channel_map& channelMap = columnMap.get_channel_map();

    for( size_t i = 0; i < columnMap.column_count(); ++i ) {
        if( columnMap.has_column_mapping( i ) ) {
            std::pair<size_t, size_t> pairing = columnMap.column_mapping( i );
            result +=
                serialize_single_column_map_to_metadata_format( channelMap[pairing.first].name(), pairing.second );
        }

        if( i < columnMap.column_count() - 1 ) {
            result += _T( "," );
        }
    }

    return result;
}

frantic::channels::channel_map parse_channel_map_from_metadata_format( const frantic::tstring& s ) {
    using namespace frantic::channels;

    std::vector<frantic::tstring> splitVec;
    boost::split( splitVec, s, boost::is_any_of( _T( "," ) ) );

    channel_map channelMap;

    for( size_t i = 0; i < splitVec.size(); ++i ) {
        size_t splitPoint = splitVec[i].find( _T( " " ) );

        if( splitPoint == frantic::tstring::npos ) {
            throw std::runtime_error( "parse_channel_map_from_metadata_format -- Invalid channel map metadata format" );
        }

        frantic::tstring channelDefStr = splitVec[i].substr( 0, splitPoint );
        boost::trim( channelDefStr );
        std::pair<data_type_t, std::size_t> channelDef = channel_data_type_and_arity_from_string( channelDefStr );
        frantic::tstring channelName = splitVec[i].substr( splitPoint + 1 );

        channelMap.define_channel( channelName, channelDef.second, channelDef.first );
    }

    channelMap.end_channel_definition();

    return channelMap;
}

frantic::channels::channel_column_map
parse_column_map_from_metadata_format( const frantic::channels::channel_map& channelMap, const frantic::tstring& s ) {
    std::vector<frantic::tstring> splitVec;
    boost::split( splitVec, s, boost::is_any_of( _T( "," ) ) );

    frantic::channels::channel_column_map result( channelMap, splitVec.size() );

    for( size_t i = 0; i < splitVec.size(); ++i ) {
        boost::trim( splitVec[i] );

        if( splitVec[i].length() > 0 ) {
            std::pair<frantic::tstring, size_t> mapping( parse_single_column_map_from_metadata_format( splitVec[i] ) );
            result.set_column_mapping( i, mapping.first, mapping.second );
        }
    }

    return result;
}

namespace {

static frantic::tstring g_columnMappingProperty( _T( "ColumnMapping" ) );
static frantic::tstring g_columnMappingChannelMapProperty( _T( "ColumnMappingChannelMap" ) );
static frantic::tstring g_textDelimiterProperty( _T( "TextDelimiter" ) );
static frantic::tstring g_headerRowCountProperty( _T( "HeaderRowCount" ) );

void set_channel_map_property( frantic::channels::property_map& pMap, const frantic::tstring& channelName,
                               const frantic::channels::channel_map& cMap ) {
    using namespace frantic::channels;

    if( pMap.has_property( channelName ) ) {
        channel_map tempMap = pMap.get_channel_map();
        if( tempMap[channelName].data_type() != frantic::channels::data_type_string ) {
            throw std::runtime_error( "set_channel_map_property -- native channel map type must be a string" );
        }
    } else {
        throw std::runtime_error( "set_channel_map_property -- no channel column map property found" );
    }

    pMap.set_cvt( channelName, serialize_channel_map_to_metadata_format( cMap ) );
}

frantic::channels::channel_map get_channel_map_property( const frantic::channels::property_map& pMap,
                                                         const frantic::tstring& channelName ) {
    if( !pMap.has_property( channelName ) ) {
        throw std::runtime_error( "get_channel_map_property -- channel map property not found" );
    }

    return parse_channel_map_from_metadata_format( pMap.get_cvt<frantic::tstring>( channelName ) );
}

void set_channel_column_map_property( frantic::channels::property_map& pMap,
                                      const frantic::tstring& columnMapChannelName,
                                      const frantic::tstring& channelMapChannelName,
                                      const frantic::channels::channel_column_map& columnMap ) {
    using namespace frantic::channels;

    set_channel_map_property( pMap, channelMapChannelName, columnMap.get_channel_map() );

    if( pMap.has_property( columnMapChannelName ) ) {
        channel_map tempMap = pMap.get_channel_map();
        if( tempMap[columnMapChannelName].data_type() != frantic::channels::data_type_string ) {
            throw std::runtime_error( "set_channel_column_map_property -- channel column map type must be a string" );
        }
    } else {
        throw std::runtime_error( "set_channel_column_map_property -- no channel column map property found" );
    }

    pMap.set_cvt( columnMapChannelName, serialize_column_map_to_metadata_format( columnMap ) );
}

frantic::channels::channel_column_map get_channel_column_map_property( const frantic::channels::property_map& pMap,
                                                                       const frantic::tstring& columnMapChannelName,
                                                                       const frantic::tstring& channelMapChannelName ) {
    if( !pMap.has_property( columnMapChannelName ) ) {
        throw std::runtime_error( "get_channel_column_map_property -- channel column map property not found" );
    }

    const frantic::channels::channel_map channelMap = get_channel_map_property( pMap, channelMapChannelName );

    return parse_column_map_from_metadata_format( channelMap, pMap.get_cvt<frantic::tstring>( columnMapChannelName ) );
}

} // namespace

bool has_column_mapping( const frantic::channels::property_map& props ) {
    try {
        get_column_mapping( props );
        return true;
    } catch( std::exception& ) {
        return false;
    }
}

frantic::channels::channel_column_map get_column_mapping( const frantic::channels::property_map& props ) {
    return get_channel_column_map_property( props, g_columnMappingProperty, g_columnMappingChannelMapProperty );
}

bool get_column_mapping( const frantic::channels::property_map& props,
                         frantic::channels::channel_column_map& outColumnMap ) {
    try {
        outColumnMap = get_column_mapping( props );
        return true;
    } catch( std::exception& ) {
        return false;
    }
}

void add_column_mapping_property( frantic::channels::property_map& props ) {
    using namespace frantic::channels;

    channel_map newMap = props.get_channel_map();
    bool added = false;

    if( !props.has_property( g_columnMappingProperty ) ) {
        newMap.append_channel( g_columnMappingProperty, 1, frantic::channels::data_type_string );
        added = true;
    } else if( props.get_channel_map()[g_columnMappingProperty].data_type() != data_type_string ) {
        throw new std::runtime_error( "prt::add_column_mapping_property -- Column map property must be a string." );
    }

    if( !props.has_property( g_columnMappingChannelMapProperty ) ) {
        newMap.append_channel( g_columnMappingChannelMapProperty, 1, frantic::channels::data_type_string );
        added = true;
    } else if( props.get_channel_map()[g_columnMappingChannelMapProperty].data_type() != data_type_string ) {
        throw new std::runtime_error(
            "prt::add_column_mapping_property -- Column map channel map property must be a string." );
    }

    if( added ) {
        props.set_channel_map_with_swap( newMap );
    }
}

void add_column_mapping_property( frantic::channels::property_map& props,
                                  const frantic::channels::channel_column_map& columnMap ) {
    add_column_mapping_property( props );
    set_column_mapping( props, columnMap );
}

void set_column_mapping( frantic::channels::property_map& props,
                         const frantic::channels::channel_column_map& columnMap ) {
    set_channel_column_map_property( props, g_columnMappingProperty, g_columnMappingChannelMapProperty, columnMap );
}

bool has_text_delimiter( const frantic::channels::property_map& props ) {
    return props.has_property( g_textDelimiterProperty );
}

char get_text_delimiter( const frantic::channels::property_map& props ) {
    if( props.has_property( g_textDelimiterProperty ) ) {
        return static_cast<char>( props.get_cvt<boost::uint8_t>( g_textDelimiterProperty ) );
    } else {
        return _T( '\0' );
    }
}

void add_text_delimiter_property( frantic::channels::property_map& props, char initial ) {
    using namespace frantic::channels;

    if( !has_text_delimiter( props ) ) {
        channel_map newMap;
        newMap.union_channel_map( props.get_channel_map() );
        newMap.define_channel( g_textDelimiterProperty, 1u, frantic::channels::data_type_uint8 );
        newMap.end_channel_definition();
        props.set_channel_map_with_swap( newMap );
    }

    set_text_delimiter( props, initial );
}

void set_text_delimiter( frantic::channels::property_map& props, char delimiter ) {
    props.set_cvt( g_textDelimiterProperty, static_cast<boost::uint8_t>( delimiter ) );
}

bool has_header_row_count( const frantic::channels::property_map& props ) {
    return props.has_property( g_headerRowCountProperty );
}

size_t get_header_row_count( const frantic::channels::property_map& props ) {
    if( has_header_row_count( props ) ) {
        return props.get_cvt<boost::uint32_t>( g_headerRowCountProperty );
    } else {
        return 0;
    }
}

void add_header_row_count_property( frantic::channels::property_map& props, size_t initial ) {
    using namespace frantic::channels;

    if( !has_header_row_count( props ) ) {
        channel_map newMap;
        newMap.union_channel_map( props.get_channel_map() );
        newMap.define_channel( g_headerRowCountProperty, 1u, frantic::channels::data_type_uint8 );
        newMap.end_channel_definition();
        props.set_channel_map_with_swap( newMap );
    }

    set_header_row_count( props, initial );
}

void set_header_row_count( frantic::channels::property_map& props, size_t count ) {
    props.set_cvt( g_headerRowCountProperty, static_cast<boost::uint32_t>( count ) );
}

} // namespace csv
} // namespace particles
} // namespace frantic
