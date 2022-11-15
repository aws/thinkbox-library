// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"
#include "utilities/scoped_temp_file.hpp"

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_file_io.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>

TEST( PLY, Basic ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "ply" );
    polymesh3_ptr output = make_triangle_polymesh();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );

    frantic::logging::null_progress_logger progress;

    write_ply_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_ply_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( PLY, Normals ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "ply" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 1, 0, 0 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 0, 1, 0 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 0, 0, 1 );

    output->add_vertex_channel( _T( "Normal" ), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_ply_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_ply_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( PLY, TextureCoords ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "ply" );
    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 1 );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 2 );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 4 );

    output->add_vertex_channel( _T( "TextureCoord" ), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    write_ply_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_ply_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( PLY, Polygon ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "ply" );

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

    write_ply_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_ply_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( PLY, VerticesOnly ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "ply" );

    polymesh3_builder builder;

    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );
    builder.add_vertex( 0, 1, 0 );

    polymesh3_ptr output = builder.finalize();
    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );

    frantic::logging::null_progress_logger progress;

    write_ply_mesh_file( path.get_path(), meshInterface, progress );
    polymesh3_ptr input = load_ply_polymesh_file( path.get_path() );

    ASSERT_TRUE( is_equal( output, input ) );
}

TEST( PLY, ColorOutOfRange ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "ply" );

    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( sizeof( float[3] ) * 3 );

    reinterpret_cast<vector3f*>( buffer.begin() )[0] = vector3f( 1.1f );
    reinterpret_cast<vector3f*>( buffer.begin() )[1] = vector3f( 1.1f );
    reinterpret_cast<vector3f*>( buffer.begin() )[2] = vector3f( 1.1f );

    output->add_vertex_channel( _T("Color"), frantic::channels::data_type_float32, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    ASSERT_NO_THROW( write_ply_mesh_file( path.get_path(), meshInterface, progress ) );
    polymesh3_ptr input = load_ply_polymesh_file( path.get_path() );

    polymesh3_const_vertex_accessor<vector3f> acc = input->get_const_vertex_accessor<vector3f>( _T("Color") );

    ASSERT_EQ( acc.vertex_count(), 3 );

    for( std::size_t i = 0; i < 3; ++i ) {
        ASSERT_EQ( vector3f( 1 ), acc.get_vertex( i ) );
    }
}

TEST( PLY, Float16Color ) {
    using namespace frantic::geometry;
    using frantic::graphics::vector3f;

    scoped_temp_file path( "ply" );

    polymesh3_ptr output = make_triangle_polymesh();

    frantic::graphics::raw_byte_buffer buffer;

    const std::size_t bufferSize =
        3 * 3 * frantic::channels::sizeof_channel_data_type( frantic::channels::data_type_float16 );
    buffer.resize( bufferSize );
    memset( buffer.begin(), 0, bufferSize );

    output->add_vertex_channel( _T("Color"), frantic::channels::data_type_float16, 3, buffer );

    mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( output ).release() );
    frantic::logging::null_progress_logger progress;

    ASSERT_NO_THROW( write_ply_mesh_file( path.get_path(), meshInterface, progress ) );
}
