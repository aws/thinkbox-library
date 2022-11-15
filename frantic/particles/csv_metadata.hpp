// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @author Stephen Kiazyk
 *
 * Contains a set of property definitions for text-based particle files
 */

#pragma once

#include <frantic/channels/channel_column_map.hpp>
#include <frantic/channels/property_map.hpp>

namespace frantic {
namespace particles {
namespace csv {

/**
 * Reads a channel column map from a string format
 */
frantic::channels::channel_map parse_channel_map_from_metadata_format( const frantic::tstring& s );

/**
 * Reads a channel map from a string format
 */
frantic::channels::channel_column_map
parse_column_map_from_metadata_format( const frantic::channels::channel_map& channelMap, const frantic::tstring& s );

/**
 * Serializes a channel column map into a string format
 */
frantic::tstring serialize_column_map_to_metadata_format( const frantic::channels::channel_column_map& columnMap );

/**
 * Serializes a channel map into a string format
 */
frantic::tstring serialize_channel_map_to_metadata_format( const frantic::channels::channel_map& cMap );

/**
 * Check if the given metadata has a valid channel column map property
 *
 * @param props the property map to query
 * @return true if a valid mapping could be found, false otherwise
 */
bool has_column_mapping( const frantic::channels::property_map& props );

/**
 * Return the column mapping propety in the given metadata, if present
 *
 * @param props the property map to query
 * @return the column map, if found
 * @throws a std::runtime_error if no property was found or an error occurred when parsing the property
 */
frantic::channels::channel_column_map get_column_mapping( const frantic::channels::property_map& props );

/**
 * Return the column mapping propety in the given metadata, if present
 *
 * @param props the property map to query
 * @param outColumnMap stores the result of the query
 * @return true if a column map could be found and parsed, false otherwise
 */
bool get_column_mapping( const frantic::channels::property_map& props,
                         frantic::channels::channel_column_map& outColumnMap );

/**
 * Adds a column mapping property to the given property map
 *
 * @param props the property map to modify
 */
void add_column_mapping_property( frantic::channels::property_map& props );

/**
 * Adds a column mapping property to the given property map, and sets its value to the given column map
 *
 * @param props the property map to modify
 * @param columnMap the initial value for the specified property
 */
void add_column_mapping_property( frantic::channels::property_map& props,
                                  const frantic::channels::channel_column_map& columnMap );

/**
 * Sets the column mapping propety in the given metadata
 *
 * @param props the property map to modify
 * @param columnMap the column map to apply to the metadata
 */
void set_column_mapping( frantic::channels::property_map& props,
                         const frantic::channels::channel_column_map& columnMap );

bool has_text_delimiter( const frantic::channels::property_map& props );

char get_text_delimiter( const frantic::channels::property_map& props );

void add_text_delimiter_property( frantic::channels::property_map& props, char initial = ',' );

void set_text_delimiter( frantic::channels::property_map& props, char delimiter );

bool has_header_row_count( const frantic::channels::property_map& props );

size_t get_header_row_count( const frantic::channels::property_map& props );

void add_header_row_count_property( frantic::channels::property_map& props, size_t initial = 0 );

void set_header_row_count( frantic::channels::property_map& props, size_t count );

} // namespace csv
} // namespace particles
} // namespace frantic
