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

TEST( OBJ, Basic ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "obj" );
    polymesh3_ptr output = make_triangle_polymesh();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );

    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_obj_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( OBJ, Normals ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "obj" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 1, 0, 0 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 0, 1, 0 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 0, 0, 1 );

    output->add_vertex_channel( _T( "Normal" ), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_obj_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( OBJ, TextureCoords ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "obj" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 1 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 2 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 4 );

    std::vector<int> polygon;
    polygon.push_back( 2 );
    polygon.push_back( 1 );
    polygon.push_back( 0 );

    output->add_vertex_channel( _T( "TextureCoord" ), frantic::channels::data_type_float32, 3, buffer, &polygon );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_obj_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( OBJ, TextureCoordsAndNormals ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "obj" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 1 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 2 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 4 );

    output->add_vertex_channel( _T( "TextureCoord" ), frantic::channels::data_type_float32, 3, buffer );

    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 1, 0, 0 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 0, 1, 0 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 0, 0, 1 );

    std::vector<int> polygon;
    polygon.push_back( 2 );
    polygon.push_back( 1 );
    polygon.push_back( 0 );

    output->add_vertex_channel( _T( "Normal" ), frantic::channels::data_type_float32, 3, buffer, &polygon );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_obj_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( OBJ, Polygon ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "obj" );
    polymesh3_ptr output = make_quad_polymesh();

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_obj_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( OBJ, VerticesOnly ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "obj" );

    polymesh3_builder builder;

    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );
    builder.add_vertex( 0, 1, 0 );

    polymesh3_ptr output = builder.finalize();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );

    frantic::logging::null_progress_logger progress;

    write_obj_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_obj_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( OBJ, GermanLocale ) {
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

    scoped_temp_file path( "obj" );

    // make sure we can read numbers that include a decimal point
    frantic::files::file_ptr f( frantic::files::tfopen( path.get_path().c_str(), _T("w+") ) );
    std::fputs( "v 0 0 0\n", f );
    std::fputs( "v 1.2 0 0\n", f );
    std::fputs( "v 1.2 1.2 0\n", f );
    std::fputs( "f 1 2 3\n", f );
    f.close();

    polymesh3_ptr mesh = load_obj_polymesh_file( path.get_path() );

    EXPECT_TRUE( bool( mesh ) );

    EXPECT_EQ( mesh->face_count(), 1 );
    EXPECT_EQ( mesh->vertex_count(), 3 );
    EXPECT_EQ( mesh->get_vertex( 0 ), vector3f( 0, 0, 0 ) );
    EXPECT_EQ( mesh->get_vertex( 1 ), vector3f( 1.2f, 0, 0 ) );
    EXPECT_EQ( mesh->get_vertex( 2 ), vector3f( 1.2f, 1.2f, 0 ) );
}
