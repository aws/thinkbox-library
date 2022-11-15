// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/files/legacy_flood_file_io.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::graphics;

namespace frantic {
namespace files {
namespace legacyflood {

///////////
// Read Functions
///////////

void read_header( std::istream& in ) {
    if( in.get() != 0 )
        throw runtime_error( "legacyflood::read_header: This legacy Flood file either uses the unsupported checksums "
                             "feature, or wasn't created by Flood." );
}

template <class T>
inline T read_raw( std::istream& in ) {
    T result = T();
    in.read( (char*)&result, sizeof( result ) );
    return result;
}

boost::int16_t read_int16( std::istream& in ) { return read_raw<boost::int16_t>( in ); }

boost::int32_t read_int32( std::istream& in ) { return read_raw<boost::int32_t>( in ); }

float read_float32( std::istream& in ) { return read_raw<float>( in ); }

// NOTE: This is implementing the functionality of the Read7BitEncodedInt function defined in the .NET BinaryReader
// class.
boost::int32_t read_7bit_encoded_int( std::istream& in ) {
    boost::int32_t result = 0;
    int shift = 0, value = 0;

    do {
        value = in.get();
        if( in.rdstate() & ifstream::failbit )
            throw std::runtime_error( "frantic::files::legacy_flood::read_7bit_encoded_int() - Read operation failed.  "
                                      "Perhaps one of your legacy flood files is corrupt?" );
        result |= ( value & 0x7f ) << shift;
        shift += 7;
    } while( ( value & 0x80 ) != 0 );

    return result;
}

// NOTE: This is implementing the functionality of the ReadString function defined in the .NET BinaryReader class.
string read_string( std::istream& in ) {
    // Get the number of bytes in the string
    boost::int32_t byteCount = read_7bit_encoded_int( in );

    // Handle errors and the empty string
    if( byteCount < 0 )
        throw std::runtime_error(
            "legacyflood::read_string: The input stream gave an invalid length for a string input." );
    if( byteCount == 0 )
        return "";

    // Read the string data
    vector<char> buffer( byteCount );
    in.read( &buffer[0], byteCount );

    if( !in )
        throw std::runtime_error(
            "legacyflood::read_string: The input stream failed to read the character data for a string input." );

    return string( &buffer[0], byteCount );
}

graphics::vector3f read_vector3f( std::istream& in ) {
    float x = read_float32( in );
    float y = read_float32( in );
    float z = read_float32( in );
    return vector3f( x, y, z );
}

graphics::vector3 read_vector3( std::istream& in ) {
    int x = read_int32( in );
    int y = read_int32( in );
    int z = read_int32( in );
    return vector3( x, y, z );
}

graphics::size3f read_size3f( std::istream& in ) {
    float xsize = read_float32( in );
    float ysize = read_float32( in );
    float zsize = read_float32( in );
    return size3f( xsize, ysize, zsize );
}

graphics::size3 read_size3( std::istream& in ) {
    int xsize = read_int32( in );
    int ysize = read_int32( in );
    int zsize = read_int32( in );
    return size3( xsize, ysize, zsize );
}

graphics::boundbox3 read_boundbox3( std::istream& in ) {
    vector3 minimum = read_vector3( in );
    vector3 maximum = read_vector3( in );
    return boundbox3( minimum, maximum );
}

template <class T>
inline void read_vector_raw( std::istream& in, std::vector<T>& outResult ) {
    int length = read_int32( in );
    if( length < 0 )
        throw std::runtime_error(
            "legacyflood::read_vector: The input stream gave a negative length for an input vector." );
    outResult.resize( length );
    if( length > 0 )
        in.read( (char*)&outResult[0], length * sizeof( T ) );
}

void read_int16_vector( std::istream& in, std::vector<boost::int16_t>& outResult ) { read_vector_raw( in, outResult ); }

void read_int32_vector( std::istream& in, std::vector<boost::int32_t>& outResult ) { read_vector_raw( in, outResult ); }

void read_float32_vector( std::istream& in, std::vector<float>& outResult ) { read_vector_raw( in, outResult ); }

void read_vector3f_vector( std::istream& in, std::vector<graphics::vector3f>& outResult ) {
    read_vector_raw( in, outResult );
}

// The multiple begin_section functions are for when there have been class name changes in the C# code, renaming the
// sections of files that we want to be able to read.
void begin_section( const std::string& streamName, std::istream& in, const std::string& sectionName ) {
    string fileSectionName = read_string( in );
    if( fileSectionName != sectionName )
        throw runtime_error( "legacyFlood::begin_section: The input stream \"" + streamName +
                             "\" didn't contain the section \"" + sectionName +
                             "\" as expected, instead it was named \"" + fileSectionName + "\"." );
}

void begin_section( const std::string& streamName, std::istream& in, const std::string& sectionNameA,
                    const std::string& sectionNameB ) {
    string fileSectionName = read_string( in );
    if( fileSectionName != sectionNameA && fileSectionName != sectionNameB )
        throw runtime_error( "legacyFlood::begin_section: The input stream \"" + streamName +
                             "\" didn't contain the section \"" + sectionNameA +
                             "\" as expected, instead it was named \"" + fileSectionName + "\"." );
}

void begin_section( const std::string& streamName, std::istream& in, const std::string& sectionNameA,
                    const std::string& sectionNameB, const std::string& sectionNameC ) {
    string fileSectionName = read_string( in );
    if( fileSectionName != sectionNameA && fileSectionName != sectionNameB && fileSectionName != sectionNameC )
        throw runtime_error( "legacyFlood::begin_section: The input stream \"" + streamName +
                             "\" didn't contain the section \"" + sectionNameA +
                             "\" as expected, instead it was named \"" + fileSectionName + "\"." );
}

///////////
// Write Functions
///////////

void write_header( std::ostream& out ) {
    // Indicate no checksums
    out.put( 0 );
}

template <class T>
inline void write_raw( std::ostream& out, const T& value ) {
    out.write( (char*)&value, sizeof( T ) );
}

void write_int16( std::ostream& out, boost::int16_t value ) { write_raw( out, value ); }

void write_int32( std::ostream& out, boost::int32_t value ) { write_raw( out, value ); }

void write_float32( std::ostream& out, float value ) { write_raw( out, value ); }

void write_7bit_encoded_int( std::ostream& out, boost::int32_t value ) {
    do {
        char byte = (char)( value & 0x7f );
        value = ( value >> 7 ) & 0x01ffffff;

        if( value != 0 )
            byte |= 0x80;

        out.put( byte );
    } while( value != 0 );
}

void write_string( std::ostream& out, const std::string& value ) {
    // Write the string data
    write_7bit_encoded_int( out, (boost::uint32_t)value.size() );
    // out.write( (char*)&wideCharConverted[0], byteCount );
    out.write( value.data(), value.size() );
}

void write_vector3f( std::ostream& out, const graphics::vector3f& value ) {
    write_float32( out, value.x );
    write_float32( out, value.y );
    write_float32( out, value.z );
}

void write_vector3( std::ostream& out, const graphics::vector3& value ) {
    write_int32( out, value.x );
    write_int32( out, value.y );
    write_int32( out, value.z );
}

void write_boundbox3( std::ostream& out, const graphics::boundbox3& value ) {
    write_vector3( out, value.minimum() );
    write_vector3( out, value.maximum() );
}

template <class T>
void write_vector_raw( std::ostream& out, const std::vector<T>& value ) {
    write_int32( out, (boost::int32_t)value.size() );
    if( value.size() > 0 )
        out.write( (const char*)&value[0], value.size() * sizeof( T ) );
}

void write_int16_vector( std::ostream& out, const std::vector<boost::int16_t>& value ) {
    write_vector_raw( out, value );
}

void write_int32_vector( std::ostream& out, const std::vector<boost::int32_t>& value ) {
    write_vector_raw( out, value );
}

void write_float32_vector( std::ostream& out, const std::vector<float>& value ) { write_vector_raw( out, value ); }

void write_vector3f_vector( std::ostream& out, const std::vector<graphics::vector3f>& value ) {
    write_vector_raw( out, value );
}

void begin_section( std::ostream& out, const std::string& sectionName ) { write_string( out, sectionName ); }

} // namespace legacyflood
} // namespace files
} // namespace frantic
