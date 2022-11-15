// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

// These functions are for reading and writing the legacy file formats used in the C# version of Flood.
// They provide a way to emulate the functionality of the BinaryReader2 and BinaryWriter2 classes.

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <iostream>

namespace frantic {
namespace files {
namespace legacyflood {

///////////
// Read Functions
///////////

void read_header( std::istream& in );

boost::int16_t read_int16( std::istream& in );
boost::int32_t read_int32( std::istream& in );
float read_float32( std::istream& in );
std::string read_string( std::istream& in );
graphics::vector3f read_vector3f( std::istream& in );
graphics::vector3 read_vector3( std::istream& in );
graphics::size3f read_size3f( std::istream& in );
graphics::size3 read_size3( std::istream& in );
graphics::boundbox3 read_boundbox3( std::istream& in );

void read_int16_vector( std::istream& in, std::vector<boost::int16_t>& outResult );
void read_int32_vector( std::istream& in, std::vector<boost::int32_t>& outResult );
void read_float32_vector( std::istream& in, std::vector<float>& outResult );
void read_vector3f_vector( std::istream& in, std::vector<graphics::vector3f>& outResult );

// The multiple begin_section functions are for when there have been class name changes in the C# code, renaming the
// sections of files that we want to be able to read.
void begin_section( const std::string& streamName, std::istream& in, const std::string& sectionName );
void begin_section( const std::string& streamName, std::istream& in, const std::string& sectionNameA,
                    const std::string& sectionNameB );
void begin_section( const std::string& streamName, std::istream& in, const std::string& sectionNameA,
                    const std::string& sectionNameB, const std::string& sectionNameC );

///////////
// Write Functions
///////////

void write_header( std::ostream& out );

void write_int16( std::ostream& out, boost::int16_t value );
void write_int32( std::ostream& out, boost::int32_t value );
void write_float32( std::ostream& out, float value );
void write_string( std::ostream& out, const std::string& value );
void write_vector3f( std::ostream& out, const graphics::vector3f& value );
void write_vector3( std::ostream& out, const graphics::vector3& value );
void write_boundbox3( std::ostream& out, const graphics::boundbox3& value );

void write_int16_vector( std::ostream& out, const std::vector<boost::int16_t>& value );
void write_int32_vector( std::ostream& out, const std::vector<boost::int32_t>& value );
void write_float32_vector( std::ostream& out, const std::vector<float>& value );
void write_vector3f_vector( std::ostream& out, const std::vector<graphics::vector3f>& value );

void begin_section( std::ostream& out, const std::string& sectionName );

} // namespace legacyflood
} // namespace files
} // namespace frantic
