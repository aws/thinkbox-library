// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <cstdio>

#include <boost/algorithm/string/join.hpp>

#include <frantic/graphics2d/vector2.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/csv_particle_istream.hpp>

#include <frantic/files/scoped_file_cleanup.hpp>

#include "gtest-helper.h"

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::files;
using frantic::graphics::vector3f;
using frantic::graphics::vector3fd;
namespace fs = boost::filesystem;

TEST( CSV, Tokenizer ) {
    fs::path csvPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( csvPath );

    frantic::tstring tmpName = frantic::files::to_tstring( csvPath );
    files::file_ptr fin;
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );

    std::fputs( "abc,\"def\",  \"gh i\",\"jkl\",m no", fin );
    std::rewind( fin );

    vector<string> tokens;
    files::read_csv_file_line( fin, _T("test"), tokens );
    ASSERT_EQ( tokens.size(), 5 );
    ASSERT_EQ( tokens[0], "abc" );
    ASSERT_EQ( tokens[1], "def" );
    ASSERT_EQ( tokens[2], "gh i" );
    ASSERT_EQ( tokens[3], "jkl" );
    ASSERT_EQ( tokens[4], "m no" );

    fin.close();
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );
    std::fputs( "a,b\nc,d\n", fin );
    std::rewind( fin );

    tokens.clear();
    files::read_csv_file_line( fin, _T("test"), tokens );
    ASSERT_EQ( tokens.size(), 2 );
    ASSERT_EQ( tokens[0], "a" );
    ASSERT_EQ( tokens[1], "b" );
    tokens.clear();
    files::read_csv_file_line( fin, _T("test"), tokens );
    ASSERT_EQ( tokens.size(), 2 );
    ASSERT_EQ( tokens[0], "c" );
    ASSERT_EQ( tokens[1], "d" );

    fin.close();
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );
    std::fputs( "\"a\nb\nc\",d,e, f ", fin );
    std::rewind( fin );

    tokens.clear();
    files::read_csv_file_line( fin, _T("test"), tokens );
    ASSERT_EQ( tokens.size(), 4 );
    ASSERT_EQ( tokens[0], "a\nb\nc" );
    ASSERT_EQ( tokens[1], "d" );
    ASSERT_EQ( tokens[2], "e" );
    ASSERT_EQ( tokens[3], " f " );

    fin.close();
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );
    std::fputs( "\"a\", \"\"\"a\"\"\", \"a\"\"b\"", fin );
    std::rewind( fin );
    tokens.clear();
    files::read_csv_file_line( fin, _T("test"), tokens );
    ASSERT_EQ( tokens.size(), 3 );
    ASSERT_EQ( tokens[0], "a" );
    ASSERT_EQ( tokens[1], "\"a\"" );
    ASSERT_EQ( tokens[2], "a\"b" );

    fin.close();
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );
    std::fputs( "a,\nb,", fin );
    std::rewind( fin );

    tokens.clear();
    files::read_csv_file_line( fin, _T("test"), tokens );
    ASSERT_EQ( tokens.size(), 2 );
    ASSERT_EQ( tokens[0], "a" );
    ASSERT_EQ( tokens[1], "" );
    tokens.clear();
    files::read_csv_file_line( fin, _T("test"), tokens );
    ASSERT_EQ( tokens.size(), 2 );
    ASSERT_EQ( tokens[0], "b" );
    ASSERT_EQ( tokens[1], "" );

    fin.close();
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );
    std::fputs( "a,b\"c,d", fin );
    std::rewind( fin );

    tokens.clear();
    EXPECT_THROW( files::read_csv_file_line( fin, _T("test"), tokens ), std::exception )
        << "An invalid quote inside of a csv token should have triggered an std::exception.";

    fin.close();
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );
    std::fputs( "a,\"b\"c,d", fin );
    std::rewind( fin );

    tokens.clear();
    EXPECT_THROW( files::read_csv_file_line( fin, _T("test"), tokens ), std::exception )
        << "An invalid quote inside of a csv token should have triggered an std::exception.";

    fin.close();
    fin.reset( frantic::files::tfopen( tmpName.c_str(), _T("w+") ) );
    std::fputs( "a,bc,\"d", fin );
    std::rewind( fin );

    tokens.clear();
    EXPECT_THROW( files::read_csv_file_line( fin, _T("test"), tokens ), std::exception )
        << "An unterminated quoted string at the end of a csv file should have triggered an std::exception.";

    fin.close();
}

TEST( CSVIStream, GermanLocale ) {
#ifdef _WIN32
#if defined( _MSC_VER ) && _MSC_VER >= 1700
    const char* german = "de-DE";
#else
    const char* german = "German";
#endif
#else
    const char* german = "de_DE.UTF-8";
#endif
    frantic::locale::set_locale_in_scope setLocale( german );

    fs::path csvPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( csvPath );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( csvPath ).c_str(), _T("w+") ) );
    std::fputs( "1.2,3.4,5.6", f );
    f.close();

    frantic::particles::particle_istream_ptr pin(
        new frantic::particles::streams::csv_particle_istream( frantic::files::to_tstring( csvPath ) ) );

    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAcc(
        pin->get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") ) );
    std::vector<char> buffer( pin->get_channel_map().structure_size() );

    ASSERT_TRUE( pin->get_particle( buffer ) );

    EXPECT_EQ( vector3f( 1.2f, 3.4f, 5.6f ), positionAcc.get( buffer ) );
}

TEST( CSVIStream, OneColumn ) {
    fs::path csvPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( csvPath );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( csvPath ).c_str(), _T("w+") ) );
    std::fputs( "1.2", f );
    f.close();

    frantic::particles::particle_istream_ptr pin;

    EXPECT_THROW(
        pin.reset( new frantic::particles::streams::csv_particle_istream( frantic::files::to_tstring( csvPath ) ) ),
        std::exception );
}

TEST( CSVIStream, TwoColumns ) {
    fs::path csvPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( csvPath );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( csvPath ).c_str(), _T("w+") ) );
    std::fputs( "1.2,3.4", f );
    f.close();

    frantic::particles::particle_istream_ptr pin;

    EXPECT_THROW(
        pin.reset( new frantic::particles::streams::csv_particle_istream( frantic::files::to_tstring( csvPath ) ) ),
        std::exception );
}

TEST( CSVIStream, SixColumns ) {
    fs::path csvPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( csvPath );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( csvPath ).c_str(), _T("w+") ) );
    std::fputs( "1 2 3 4 5 6", f );
    f.close();

    frantic::particles::particle_istream_ptr pin(
        new frantic::particles::streams::csv_particle_istream( frantic::files::to_tstring( csvPath ) ) );

    EXPECT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );
    EXPECT_EQ( pin->get_channel_map()[_T("Position")].arity(), 3 );
    EXPECT_TRUE( pin->get_channel_map().has_channel( _T("Color") ) );
    EXPECT_EQ( pin->get_channel_map()[_T("Color")].arity(), 3 );
}

TEST( CSVIStream, DataColumns ) {
    for( int columnCount = 4; columnCount <= 12; ++columnCount ) {
        if( columnCount == 6 ) {
            // 6 columns are interpreted as Position and Color
            continue;
        }

        const int dataColumnCount = columnCount - 3;

        fs::path csvPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

        scoped_file_cleanup fileCleanup;
        fileCleanup.add( csvPath );

        std::vector<std::string> columns;
        for( int column = 1; column <= columnCount; ++column ) {
            columns.push_back( boost::lexical_cast<std::string>( column ) );
        }

        frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( csvPath ).c_str(), _T("w+") ) );
        std::fputs( boost::algorithm::join( columns, " " ).c_str(), f );
        f.close();

        frantic::particles::particle_istream_ptr pin(
            new frantic::particles::streams::csv_particle_istream( frantic::files::to_tstring( csvPath ) ) );

        ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );
        ASSERT_EQ( 3, pin->get_channel_map()[_T("Position")].arity() );
        for( int i = 1; i <= dataColumnCount; ++i ) {
            const frantic::tstring channelName = _T("Data") + boost::lexical_cast<frantic::tstring>( i );
            EXPECT_TRUE( pin->get_channel_map().has_channel( channelName ) );
            EXPECT_EQ( pin->get_channel_map()[channelName].arity(), 1 );
        }
    }
}

// This is specifically just meant to check if it can read a 'normal' header row
TEST( CSVIStream, GuessFormat ) {
    using namespace frantic::channels;
    using namespace frantic::particles::streams;

    std::vector<std::string> columns;
    columns.push_back( "float32 Viscosity[0]" );
    columns.push_back( "float32 Viscosity[1]" );
    columns.push_back( "float32 Viscosity[2]" );
    columns.push_back( "int8 Luster" );
    columns.push_back( "uint16 Color[0]" );
    columns.push_back( "uint16 Color[1]" );

    channel_column_map columnMap;
    csv_particle_istream::guess_channel_column_map( columns, _T(""), false, frantic::channels::data_type_float32,
                                                    columnMap );

    channel_map underlyingMap( columnMap.get_channel_map() );

    ASSERT_EQ( 3, underlyingMap.channel_count() );
    ASSERT_TRUE( underlyingMap.has_channel( _T("Viscosity") ) );
    ASSERT_EQ( 3, underlyingMap[_T("Viscosity")].arity() );
    size_t viscosityIndex = underlyingMap.channel_index( _T("Viscosity") );
    size_t currentColumn = 0;
    for( size_t i = 0; i < 3; ++i ) {
        EXPECT_EQ( viscosityIndex, columnMap.column_mapping( currentColumn ).first );
        EXPECT_EQ( i, columnMap.column_mapping( currentColumn ).second );
        ++currentColumn;
    }

    ASSERT_TRUE( underlyingMap.has_channel( _T("Luster") ) );
    ASSERT_EQ( 1, underlyingMap[_T("Luster")].arity() );
    size_t lusterIndex = underlyingMap.channel_index( _T("Luster") );
    EXPECT_EQ( lusterIndex, columnMap.column_mapping( currentColumn ).first );
    EXPECT_EQ( 0, columnMap.column_mapping( currentColumn ).second );
    ++currentColumn;

    ASSERT_TRUE( underlyingMap.has_channel( _T("Color") ) );
    ASSERT_EQ( 2, underlyingMap[_T("Color")].arity() );
    size_t colorIndex = underlyingMap.channel_index( _T("Color") );
    for( size_t i = 0; i < 2; ++i ) {
        EXPECT_EQ( colorIndex, columnMap.column_mapping( currentColumn ).first );
        EXPECT_EQ( i, columnMap.column_mapping( currentColumn ).second );
        ++currentColumn;
    }
}

TEST( CSVIStream, WithColumnMap ) {
    using namespace frantic::channels;
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;

    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "0 1 2 3 4 5\n", f );
    f.close();

    channel_map channelMap;
    channelMap.define_channel( _T("P"), 3, data_type_float32 );
    channelMap.define_channel( _T("V"), 1, data_type_uint32 );
    channelMap.define_channel( _T("I"), 2, data_type_int32 );
    channelMap.end_channel_definition( 1, true );

    channel_column_map columnMap( channelMap, 6 );
    columnMap.set_column_mapping( 0, _T("V"), 0 );
    columnMap.set_column_mapping( 1, _T("P"), 2 );
    columnMap.set_column_mapping( 2, _T("I"), 0 );
    columnMap.set_column_mapping( 3, _T("P"), 1 );
    columnMap.set_column_mapping( 4, _T("P"), 0 );
    columnMap.set_column_mapping( 5, _T("I"), 1 );

    frantic::particles::particle_istream_ptr pin(
        new frantic::particles::streams::csv_particle_istream( frantic::files::to_tstring( path ), columnMap ) );

    std::vector<char> buffer( pin->get_channel_map().structure_size() );

    pin->get_particle( buffer );

    frantic::channels::channel_accessor<vector3f> pAcc( pin->get_channel_map().get_accessor<vector3f>( _T("P") ) );
    frantic::channels::channel_accessor<vector2> iAcc( pin->get_channel_map().get_accessor<vector2>( _T("I") ) );
    frantic::channels::channel_accessor<uint32_t> vAcc( pin->get_channel_map().get_accessor<uint32_t>( _T("V") ) );

    vector3f P = pAcc.get( buffer );
    EXPECT_FLOAT_EQ( 4.f, P[0] );
    EXPECT_FLOAT_EQ( 3.f, P[1] );
    EXPECT_FLOAT_EQ( 1.f, P[2] );
    vector2 I = iAcc.get( buffer );
    EXPECT_EQ( 2, I[0] );
    EXPECT_EQ( 5, I[1] );
    uint32_t V = vAcc.get( buffer );
    EXPECT_EQ( 0, V );
}

TEST( CSVIStream, Float64Position ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    // Test the code path without a header
    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "1.234567890123456 3.14159265358979 1.1e200 3 4 5\n", f );
    f.close();

    // Test this through the factory object, set to request float64 Position data
    particles::particle_file_stream_factory_object pfactory;
    frantic::particles::particle_istream_ptr pin;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> positionAcc;
    vector<char> buffer;
    vector3fd v;

    pfactory.set_position_type_hint( channels::data_type_float64 );
    pin = pfactory.create_istream( frantic::files::to_tstring( path ) );
    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );
    ASSERT_EQ( 3, pin->get_channel_map()[_T("Position")].arity() );
    ASSERT_EQ( channels::data_type_float64, pin->get_channel_map()[_T("Position")].data_type() );
    positionAcc = pin->get_channel_map().get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
    buffer.resize( pin->get_channel_map().structure_size() );
    ASSERT_TRUE( pin->get_particle( buffer ) );
    v = positionAcc.get( buffer );
    EXPECT_DOUBLE_EQ( 1.234567890123456, v.x );
    EXPECT_DOUBLE_EQ( 3.14159265358979, v.y );
    EXPECT_DOUBLE_EQ( 1.1e200, v.z );

    // Test the code path with a header
    f.reset( frantic::files::tfopen( frantic::files::to_tstring( path ).c_str(), _T("w+") ) );
    std::fputs( "float32 Density, float16 Color, float16 Color, float16 Color, "
                "float32 Position, float32 Position, float32 Position, int32 Data\n",
                f );
    std::fputs( "0.5, 1, 1, 0.3, 1.111111111111111, -2.2, 1e-100, 12345\n", f );
    f.close();

    pfactory.set_position_type_hint( channels::data_type_float64 );
    pin = pfactory.create_istream( frantic::files::to_tstring( path ) );
    ASSERT_TRUE( pin->get_channel_map().has_channel( _T("Position") ) );
    ASSERT_EQ( 3, pin->get_channel_map()[_T("Position")].arity() );
    ASSERT_EQ( channels::data_type_float64, pin->get_channel_map()[_T("Position")].data_type() );
    positionAcc = pin->get_channel_map().get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
    buffer.resize( pin->get_channel_map().structure_size() );
    ASSERT_TRUE( pin->get_particle( buffer ) );
    v = positionAcc.get( buffer );
    EXPECT_DOUBLE_EQ( 1.111111111111111, v.x );
    EXPECT_DOUBLE_EQ( -2.2, v.y );
    EXPECT_DOUBLE_EQ( 1e-100, v.z );
}

namespace {

// These functions are for getting macro functions and templates to work together cleanly below
inline void expect_vector_eq( const vector3f& expectedPoint, const vector3f& actualPoint ) {
    EXPECT_VECTOR3F_EQ( expectedPoint, actualPoint );
}

inline void expect_vector_eq( const vector3fd& expectedPoint, const vector3fd& actualPoint ) {
    EXPECT_VECTOR3FD_EQ( expectedPoint, actualPoint );
}

template <typename FloatType>
void testPrecisionOStream() {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.csv" );

    typedef graphics::vector3t<FloatType> vector_t;

    channels::channel_map cm;
    cm.define_channel<vector_t>( _T("Position") );
    cm.end_channel_definition();

    std::vector<vector_t> testPoints;
    testPoints.push_back( vector_t( 3, 4, 5 ) );
    testPoints.push_back(
        vector_t( FloatType( 1.234567890123456 ), FloatType( 3.14159265358979 ), FloatType( 1.1e20 ) ) );
    testPoints.push_back(
        vector_t( FloatType( -6543210987654321 ), FloatType( 2718281828.459 ), FloatType( 1.01e-20 ) ) );
    testPoints.push_back( vector_t( std::numeric_limits<FloatType>::max(), std::numeric_limits<FloatType>::epsilon(),
                                    std::numeric_limits<FloatType>::min() ) );

    channels::channel_accessor<vector_t> positionAccessor = cm.get_accessor<vector_t>( _T("Position") );
    std::vector<char> buffer( cm.structure_size() );

    particles::particle_file_stream_factory_object pfactory;

    {
        particles::particle_ostream_ptr ostream = pfactory.create_ostream( frantic::files::to_tstring( path ), cm, cm );

        for( typename std::vector<vector_t>::const_iterator iter = testPoints.begin(); iter != testPoints.end();
             ++iter ) {
            positionAccessor( &buffer[0] ) = *iter;
            ostream->put_particle( &buffer[0] );
        }

        ostream->close();
    }

    {
        particles::particle_istream_ptr istream = pfactory.create_istream( frantic::files::to_tstring( path ), cm );

        for( typename std::vector<vector_t>::const_iterator iter = testPoints.begin(); iter != testPoints.end();
             ++iter ) {
            const bool ok = istream->get_particle( &buffer[0] );
            ASSERT_TRUE( ok );

            const vector_t& expectedPoint = *iter;
            const vector_t& actualPoint = positionAccessor( &buffer[0] );
            expect_vector_eq( expectedPoint, actualPoint );
        }
        const bool endOfStreamCheck = istream->get_particle( &buffer[0] );
        EXPECT_FALSE( endOfStreamCheck );
    }
}

} // anonymous namespace

TEST( CSVIStream, Float32PositionOStream ) { testPrecisionOStream<float>(); }

TEST( CSVIStream, Float64PositionOStream ) { testPrecisionOStream<double>(); }

// Tests a simple file with headers and HDR colors
TEST( CSV, HDRColors ) {
    particles::particle_file_stream_factory_object pfactory;
    frantic::particles::particle_istream_ptr pin;
    pin = pfactory.create_istream( _T("TestInputs/hdri_colors.csv") );
    ASSERT_EQ( 2u, pin->get_channel_map().channel_count() );
    // Both Position and Color are 32-bit floats
    ASSERT_EQ( channels::data_type_float32, pin->get_channel_map()[0].data_type() );
    ASSERT_EQ( 3u, pin->get_channel_map()[0].arity() );
    EXPECT_EQ( _T("Position"), pin->get_channel_map()[0].name() );
    ASSERT_EQ( channels::data_type_float32, pin->get_channel_map()[1].data_type() );
    ASSERT_EQ( 3u, pin->get_channel_map()[1].arity() );
    EXPECT_EQ( _T("Color"), pin->get_channel_map()[1].name() );

    struct {
        vector3f pos;
        vector3f col;
    } points[4];
    size_t numPoints = 4;

    ASSERT_TRUE( pin->get_particles( (char*)&points, numPoints ) );
    ASSERT_EQ( 4u, numPoints );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 0 ), points[0].pos );
    EXPECT_VECTOR3F_EQ( vector3f( 5, 3, 1 ), points[0].col );
    EXPECT_VECTOR3F_EQ( vector3f( 1, 0, 0 ), points[1].pos );
    EXPECT_VECTOR3F_EQ( vector3f( 0.001f, 0.002f, 0.0005f ), points[1].col );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 1, 0 ), points[2].pos );
    EXPECT_VECTOR3F_EQ( vector3f( 0.5f, 1.0f, 0.4f ), points[2].col );
    EXPECT_VECTOR3F_EQ( vector3f( 0, 0, 1 ), points[3].pos );
    EXPECT_VECTOR3F_EQ( vector3f( 0.8f, 1.2f, 0.8f ), points[3].col );
}

TEST( CSV, InfiniteLoop ) {
    frantic::files::file_ptr maliciousFile( std::fopen( "TestInputs/malicious_whitespace_after_quote.csv", "r" ) );
    ASSERT_TRUE( maliciousFile );

    csv_reader reader( maliciousFile );
    std::vector<std::string> output;

    // If the test fails, this call will result in an infinite loop.
    // If run in the context of CI, it will timeout, otherwise, it will
    // need to be stopped manually.
    reader.read_line( output );

    ASSERT_EQ( output.size(), 1 );
    ASSERT_EQ( output[0], "" );
}
