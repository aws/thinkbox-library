// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"
#include "utilities/scoped_temp_file.hpp"

#include <frantic/files/files.hpp>
#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_file_io.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>

#include <frantic/locale/locale.hpp>

TEST( STL, AsciiBasic ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "stl" );
    polymesh3_ptr output = make_triangle_polymesh();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );

    frantic::logging::null_progress_logger progress;

    write_ascii_stl_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_stl_polymesh_file( path.get_path() );

    EXPECT_TRUE( is_consistent_topology( output, input ) );
    EXPECT_TRUE( input->has_face_channel( _T( "Normal" ) ) );
}

TEST( STL, AsciiNormals ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "stl" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 0, 0.5, 1 );

    output->add_face_channel( _T( "Normal" ), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_ascii_stl_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_stl_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( STL, AsciiGermanLocale ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

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

    scoped_temp_file path( "stl" );

    // make sure we can read numbers that include a decimal point
    frantic::files::file_ptr f( frantic::files::tfopen( path.get_path().c_str(), _T("w+") ) );
    std::fputs( "solid\n", f );
    std::fputs( "  facet normal 0.0 0.0 0.0\n", f );
    std::fputs( "    outer loop\n", f );
    std::fputs( "      vertex 1.2 0.0 0.0\n", f );
    std::fputs( "      vertex 0.0 3.4 0.0\n", f );
    std::fputs( "      vertex 0.0 0.0 5.6\n", f );
    std::fputs( "    endloop\n", f );
    std::fputs( "  endfacet\n", f );
    std::fputs( "endsolid\n", f );
    f.close();

    polymesh3_ptr mesh = load_stl_polymesh_file( path.get_path() );

    EXPECT_TRUE( bool( mesh ) );

    EXPECT_EQ( mesh->face_count(), 1 );
    EXPECT_EQ( mesh->vertex_count(), 3 );
    EXPECT_EQ( mesh->get_vertex( 0 ), vector3f( 1.2f, 0, 0 ) );
    EXPECT_EQ( mesh->get_vertex( 1 ), vector3f( 0, 3.4f, 0 ) );
    EXPECT_EQ( mesh->get_vertex( 2 ), vector3f( 0, 0, 5.6f ) );
}

TEST( STL, BinaryBasic ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "stl" );
    polymesh3_ptr output = make_triangle_polymesh();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );

    frantic::logging::null_progress_logger progress;

    write_binary_stl_mesh_file( path.get_path(), meshInterface, true, progress );
    polymesh3_ptr input = load_stl_polymesh_file( path.get_path() );

    EXPECT_TRUE( is_consistent_topology( output, input ) );
    EXPECT_TRUE( input->has_face_channel( _T( "Normal" ) ) );
}

TEST( STL, BinaryNormals ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "stl" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 0, 0.5, 1 );

    output->add_face_channel( _T( "Normal" ), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_binary_stl_mesh_file( path.get_path(), meshInterface, true, progress );
    polymesh3_ptr input = load_stl_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( STL, BinaryColorsSolid ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "stl" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 0, 0.5, 0.75 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 0, 0.5, 0.75 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 0, 0.5, 0.75 );

    output->add_vertex_channel( _T( "Color" ), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_binary_stl_mesh_file( path.get_path(), meshInterface, true, progress );
    polymesh3_ptr input = load_stl_polymesh_file( path.get_path() );

    polymesh3_const_vertex_accessor<vector3f> acc = input->get_const_vertex_accessor<vector3f>( _T("Color") );

    EXPECT_EQ( acc.vertex_count(), 3 );

    vector3f color;
    for( std::size_t i = 0; i < 3; ++i ) {
        color = acc.get_vertex( i );
        EXPECT_NEAR( color.x, 0.0, 1.0f / 32.0f );
        EXPECT_NEAR( color.y, 0.5, 1.0f / 32.0f );
        EXPECT_NEAR( color.z, 0.75, 1.0f / 32.0f );
    }
}

TEST( STL, BinaryColorsMagics ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "stl" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 0, 0.5, 0.75 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 0, 0.5, 0.75 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 0, 0.5, 0.75 );

    output->add_vertex_channel( _T( "Color" ), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_binary_stl_mesh_file( path.get_path(), meshInterface, false, progress );
    polymesh3_ptr input = load_stl_polymesh_file( path.get_path() );

    polymesh3_const_vertex_accessor<vector3f> acc = input->get_const_vertex_accessor<vector3f>( _T("Color") );

    EXPECT_EQ( acc.vertex_count(), 3 );

    vector3f color;
    for( std::size_t i = 0; i < 3; ++i ) {
        color = acc.get_vertex( i );
        EXPECT_NEAR( color.x, 0, 1.0f / 32.0f );
        EXPECT_NEAR( color.y, 0.5, 1.0f / 32.0f );
        EXPECT_NEAR( color.z, 0.75, 1.0f / 32.0f );
    }
}

TEST( STL, InvalidFace ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "stl" );

    polymesh3_builder builder;

    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );
    builder.add_vertex( 0, 1, 0 );

    std::vector<int> polygon;
    polygon.push_back( 0 );
    polygon.push_back( 1 );
    polygon.push_back( 2 );
    polygon.push_back( 3 );
    builder.add_polygon( polygon );

    polymesh3_ptr output = builder.finalize();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );

    frantic::logging::null_progress_logger progress;

    ASSERT_THROW( write_ascii_stl_mesh_file( path.get_path(), meshInterface, progress ), std::runtime_error );
}
