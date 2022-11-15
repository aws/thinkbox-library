// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"

#include <frantic/geometry/const_shared_polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_accessors.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_file_io.hpp>

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>

#include <frantic/graphics/vector3f.hpp>

#include <frantic/files/files.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/filesystem/path.hpp>

TEST( Polymesh3, BasicTest ) {
    using namespace frantic::geometry;

    polymesh3_builder builder;

    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );

    std::vector<int> polygon;
    polygon.push_back( 0 );
    polygon.push_back( 1 );
    polygon.push_back( 2 );
    builder.add_polygon( polygon );

    polymesh3_ptr polymesh = builder.finalize();

    frantic::graphics::raw_byte_buffer buffer;

    buffer.resize( 3 * sizeof( frantic::graphics::vector3f ) );
    reinterpret_cast<frantic::graphics::vector3f*>( buffer.ptr_at( 0 ) )[0].set( 1, 0, 0 );
    reinterpret_cast<frantic::graphics::vector3f*>( buffer.ptr_at( 0 ) )[1].set( 0, 2, 0 );
    reinterpret_cast<frantic::graphics::vector3f*>( buffer.ptr_at( 0 ) )[2].set( 0, 0, 3 );
    polymesh->add_vertex_channel( _T("VertexDataChannel"), frantic::channels::data_type_float32, 3, buffer );

    buffer.resize( sizeof( float ) );
    reinterpret_cast<float*>( buffer.ptr_at( 0 ) )[0] = 1.f;
    polymesh->add_face_channel( _T("FaceDataChannel"), frantic::channels::data_type_float32, 1, buffer );

    EXPECT_TRUE( polymesh->has_channel( _T("VertexDataChannel") ) );
    EXPECT_TRUE( polymesh->has_vertex_channel( _T("VertexDataChannel") ) );
    EXPECT_TRUE( !polymesh->has_face_channel( _T("VertexDataChannel") ) );

    EXPECT_TRUE( polymesh->has_channel( _T("FaceDataChannel") ) );
    EXPECT_TRUE( polymesh->has_face_channel( _T("FaceDataChannel") ) );
    EXPECT_TRUE( !polymesh->has_vertex_channel( _T("FaceDataChannel") ) );

    // Shouldn't be able to get a face accessor for a vertex channel
    try {
        EXPECT_ANY_THROW( polymesh->get_face_accessor( _T("VertexDataChannel") ) );
    } catch( std::exception& ) {
    }
    // Shouldn't be able to get an accessor of the wrong type
    try {
        EXPECT_ANY_THROW( polymesh->get_vertex_accessor<float>( _T("VertexDataChannel") ) );
    } catch( std::exception& ) {
    }
    try {
        EXPECT_ANY_THROW( polymesh->get_vertex_accessor<frantic::graphics::vector3>( _T("VertexDataChannel") ) );
    } catch( std::exception& ) {
    }
    // Test typed accessor
    {
        polymesh3_vertex_accessor<frantic::graphics::vector3f> acc =
            polymesh->get_vertex_accessor<frantic::graphics::vector3f>( _T("VertexDataChannel") );
        EXPECT_EQ( acc.get_vertex( 0 ), frantic::graphics::vector3f( 1, 0, 0 ) );
    }
    // Test general accessor
    {
        polymesh3_vertex_accessor<void> acc = polymesh->get_vertex_accessor( _T("VertexDataChannel") );
        EXPECT_EQ( *reinterpret_cast<frantic::graphics::vector3f*>( acc.get_vertex( 1 ) ),
                   frantic::graphics::vector3f( 0, 2, 0 ) );
    }

    // Shouldn't be able to get a vertex accessor for a face channel
    try {
        EXPECT_ANY_THROW( polymesh->get_vertex_accessor( _T("FaceDataChannel") ) );
    } catch( std::exception& ) {
    }
    // Shouldn't be able to get an accessor of the wrong type
    try {
        EXPECT_ANY_THROW( polymesh->get_face_accessor<frantic::graphics::vector3f>( _T("FaceDataChannel") ) );
    } catch( std::exception& ) {
    }
    try {
        EXPECT_ANY_THROW( polymesh->get_face_accessor<int>( _T("FaceDataChannel") ) );
    } catch( std::exception& ) {
    }
    // Test typed face accessorf
    {
        polymesh3_face_accessor<float> acc = polymesh->get_face_accessor<float>( _T("FaceDataChannel") );
        EXPECT_EQ( acc.get_face( 0 ), 1.f );
    }
    // Test general face accessor
    {
        polymesh3_face_accessor<void> acc = polymesh->get_face_accessor( _T("FaceDataChannel") );
        EXPECT_EQ( *reinterpret_cast<float*>( acc.get_face( 0 ) ), 1.f );
    }

    polymesh->add_empty_face_channel( _T("EmptyFaceChannel"), frantic::channels::data_type_float32, 3 );
    {
        polymesh3_const_face_accessor<frantic::graphics::vector3f> acc =
            polymesh->get_const_face_accessor<frantic::graphics::vector3f>( _T("EmptyFaceChannel") );
        EXPECT_EQ( acc.face_count(), 1 );
    }
}

frantic::geometry::polymesh3_ptr build_polymesh3() {
    using namespace frantic::geometry;

    polymesh3_builder builder;

    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );

    std::vector<int> polygon;
    polygon.push_back( 0 );
    polygon.push_back( 1 );
    polygon.push_back( 2 );
    builder.add_polygon( polygon );

    polymesh3_ptr polymesh = builder.finalize();
    return polymesh;
}

frantic::geometry::polymesh3_ptr build_polymesh3_with_channels() {
    using namespace frantic::geometry;

    polymesh3_ptr polymesh = build_polymesh3();

    frantic::graphics::raw_byte_buffer buffer;

    buffer.resize( 3 * sizeof( frantic::graphics::vector3f ) );
    reinterpret_cast<frantic::graphics::vector3f*>( buffer.ptr_at( 0 ) )[0].set( 1, 0, 0 );
    reinterpret_cast<frantic::graphics::vector3f*>( buffer.ptr_at( 0 ) )[1].set( 0, 2, 0 );
    reinterpret_cast<frantic::graphics::vector3f*>( buffer.ptr_at( 0 ) )[2].set( 0, 0, 3 );
    polymesh->add_vertex_channel( _T("VertexDataChannel"), frantic::channels::data_type_float32, 3, buffer );

    buffer.resize( sizeof( float ) );
    reinterpret_cast<float*>( buffer.ptr_at( 0 ) )[0] = 1.f;
    polymesh->add_face_channel( _T("FaceDataChannel"), frantic::channels::data_type_float32, 1, buffer );

    return polymesh;
}

TEST( Polymesh3, Accessors ) {
    using namespace frantic::geometry;

    // vertex accessors
    {
        polymesh3_vertex_accessor<frantic::graphics::vector3f> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_vertex_accessor<frantic::graphics::vector3f> acc =
            polymesh->get_vertex_accessor<frantic::graphics::vector3f>( _T("VertexDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( acc.get_vertex( 0 ), frantic::graphics::vector3f( 1, 0, 0 ) );
        acc.get_vertex( 0 ) = frantic::graphics::vector3f( 2, 0, 0 );
        EXPECT_EQ( acc.get_vertex( 0 ), frantic::graphics::vector3f( 2, 0, 0 ) );
    }
    // const
    {
        polymesh3_const_vertex_accessor<frantic::graphics::vector3f> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        const polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_const_vertex_accessor<frantic::graphics::vector3f> acc =
            polymesh->get_const_vertex_accessor<frantic::graphics::vector3f>( _T("VertexDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( acc.get_vertex( 0 ), frantic::graphics::vector3f( 1, 0, 0 ) );
    }
    // void
    {
        polymesh3_vertex_accessor<void> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_vertex_accessor<void> acc = polymesh->get_vertex_accessor( _T("VertexDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( *reinterpret_cast<frantic::graphics::vector3f*>( acc.get_vertex( 1 ) ),
                   frantic::graphics::vector3f( 0, 2, 0 ) );
    }
    // const void
    {
        polymesh3_const_vertex_accessor<void> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        const polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_const_vertex_accessor<void> acc = polymesh->get_const_vertex_accessor( _T("VertexDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( *reinterpret_cast<const frantic::graphics::vector3f*>( acc.get_vertex( 1 ) ),
                   frantic::graphics::vector3f( 0, 2, 0 ) );
    }
    // cvt
    {
        polymesh3_cvt_vertex_accessor<frantic::graphics::vector3fd> cvtAcc;
        EXPECT_TRUE( !cvtAcc.is_valid() );
    }
    {
        polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_cvt_vertex_accessor<frantic::graphics::vector3fd> cvtAcc =
            polymesh->get_cvt_vertex_accessor<frantic::graphics::vector3fd>( _T("VertexDataChannel") );
        EXPECT_TRUE( cvtAcc.is_valid() );
        EXPECT_EQ( cvtAcc.get_vertex( 0 ), frantic::graphics::vector3fd( 1, 0, 0 ) );
        cvtAcc.set_vertex( 0, frantic::graphics::vector3fd( 2, 0, 0 ) );
        EXPECT_EQ( cvtAcc.get_vertex( 0 ), frantic::graphics::vector3fd( 2, 0, 0 ) );
        EXPECT_EQ( frantic::graphics::vector3fd( 1, 0, 0 ), frantic::graphics::vector3fd( 1, 0, 0 ) );
    }
    // const cvt
    {
        polymesh3_const_cvt_vertex_accessor<frantic::graphics::vector3fd> cvtAcc;
        EXPECT_TRUE( !cvtAcc.is_valid() );
    }
    {
        const polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_const_cvt_vertex_accessor<frantic::graphics::vector3fd> cvtAcc =
            polymesh->get_const_cvt_vertex_accessor<frantic::graphics::vector3fd>( _T("VertexDataChannel") );
        EXPECT_TRUE( cvtAcc.is_valid() );
        EXPECT_EQ( cvtAcc.get_vertex( 0 ), frantic::graphics::vector3fd( 1, 0, 0 ) );
    }

    // face accessors
    {
        polymesh3_face_accessor<float> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_face_accessor<float> acc = polymesh->get_face_accessor<float>( _T("FaceDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( acc.get_face( 0 ), 1.f );
        acc.get_face( 0 ) = 2.f;
        EXPECT_EQ( acc.get_face( 0 ), 2.f );
    }
    // const
    {
        polymesh3_const_face_accessor<float> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        const polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_const_face_accessor<float> acc = polymesh->get_const_face_accessor<float>( _T("FaceDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( acc.get_face( 0 ), 1.f );
    }
    // cvt
    {
        polymesh3_cvt_face_accessor<double> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_cvt_face_accessor<double> acc = polymesh->get_cvt_face_accessor<double>( _T("FaceDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( acc.get_face( 0 ), 1.0 );
        acc.set_face( 0, 2.0 );
        EXPECT_EQ( acc.get_face( 0 ), 2.0 );
    }
    // const cvt
    {
        polymesh3_const_cvt_face_accessor<double> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        const polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_const_cvt_face_accessor<double> acc =
            polymesh->get_const_cvt_face_accessor<double>( _T("FaceDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( acc.get_face( 0 ), 1.0 );
    }
    // void
    {
        polymesh3_face_accessor<void> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_face_accessor<void> acc = polymesh->get_face_accessor( _T("FaceDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( *reinterpret_cast<float*>( acc.get_face( 0 ) ), 1.f );
    }
    // const void
    {
        polymesh3_const_face_accessor<void> acc;
        EXPECT_TRUE( !acc.is_valid() );
    }
    {
        const polymesh3_ptr polymesh = build_polymesh3_with_channels();
        polymesh3_const_face_accessor<void> acc = polymesh->get_const_face_accessor( _T("FaceDataChannel") );
        EXPECT_TRUE( acc.is_valid() );
        EXPECT_EQ( *reinterpret_cast<const float*>( acc.get_face( 0 ) ), 1.f );
    }
}

TEST( Polymesh3, ScaleChannel ) {
    using namespace frantic::geometry;

    polymesh3_ptr inMesh = build_polymesh3_with_channels();
    polymesh3_ptr outMesh = build_polymesh3_with_channels();

    scale_channel( outMesh, _T("VertexDataChannel"), 2 );
    scale_channel( outMesh, _T("FaceDataChannel"), 4 );

    { // Scope for VertexDataChannel accessor
        typedef frantic::graphics::vector3f data_t;
        const frantic::tstring name = _T("VertexDataChannel");

        polymesh3_vertex_accessor<data_t> inAcc = inMesh->get_vertex_accessor<data_t>( name );
        polymesh3_vertex_accessor<data_t> outAcc = outMesh->get_vertex_accessor<data_t>( name );

        ASSERT_EQ( inAcc.vertex_count(), outAcc.vertex_count() );
        EXPECT_TRUE( inAcc.vertex_count() > 0 );

        for( int i = 0; i < inAcc.vertex_count(); ++i ) {
            EXPECT_EQ( 2 * inAcc.get_vertex( i ), outAcc.get_vertex( i ) );
        }
    }

    { // Scope for FaceDataChannel accessor
        typedef float data_t;
        const frantic::tstring name = _T("FaceDataChannel");

        polymesh3_face_accessor<data_t> inAcc = inMesh->get_face_accessor<data_t>( name );
        polymesh3_face_accessor<data_t> outAcc = outMesh->get_face_accessor<data_t>( name );

        ASSERT_EQ( inAcc.face_count(), outAcc.face_count() );
        EXPECT_TRUE( inAcc.face_count() > 0 );

        for( int i = 0; i < inAcc.face_count(); ++i ) {
            EXPECT_EQ( 4 * inAcc.get_face( i ), outAcc.get_face( i ) );
        }
    }
}

TEST( Polymesh3, IsConsistentTopology ) {
    using namespace frantic::geometry;

    { // same mesh
        polymesh3_ptr a = build_polymesh3();
        polymesh3_ptr b = build_polymesh3();
        EXPECT_TRUE( is_consistent_topology( a, b ) );
    }

    { // different corner count
        polymesh3_ptr a = build_polymesh3();

        polymesh3_builder builder;

        builder.add_vertex( 0, 0, 0 );
        builder.add_vertex( 1, 0, 0 );
        builder.add_vertex( 1, 1, 0 );

        std::vector<int> polygon;
        polygon.push_back( 0 );
        polygon.push_back( 1 );
        builder.add_polygon( polygon );

        polymesh3_ptr b = builder.finalize();

        EXPECT_TRUE( !is_consistent_topology( a, b ) );
        EXPECT_TRUE( !is_consistent_topology( b, a ) );
    }

    { // different face count
        polymesh3_ptr a = build_polymesh3();

        polymesh3_builder builder;

        builder.add_vertex( 0, 0, 0 );
        builder.add_vertex( 1, 0, 0 );
        builder.add_vertex( 1, 1, 0 );

        std::vector<int> polygon;
        polygon.push_back( 0 );
        polygon.push_back( 1 );
        polygon.push_back( 2 );
        builder.add_polygon( polygon );
        builder.add_polygon( polygon );

        polymesh3_ptr b = builder.finalize();

        EXPECT_TRUE( !is_consistent_topology( a, b ) );
        EXPECT_TRUE( !is_consistent_topology( b, a ) );
    }

    { // different face indices
        polymesh3_ptr a, b;
        {
            polymesh3_builder builder;

            builder.add_vertex( 0, 0, 0 );
            builder.add_vertex( 1, 0, 0 );
            builder.add_vertex( 1, 1, 0 );
            builder.add_vertex( 0, 1, 0 );

            std::vector<int> polygon;
            polygon.push_back( 0 );
            polygon.push_back( 1 );
            polygon.push_back( 2 );
            builder.add_polygon( polygon );
            builder.add_polygon( polygon );

            a = builder.finalize();
        }

        {
            polymesh3_builder builder;

            builder.add_vertex( 0, 0, 0 );
            builder.add_vertex( 1, 0, 0 );
            builder.add_vertex( 1, 1, 0 );
            builder.add_vertex( 0, 1, 0 );

            std::vector<int> polygon;
            polygon.push_back( 0 );
            polygon.push_back( 1 );
            polygon.push_back( 2 );
            builder.add_polygon( polygon );

            polygon[2] = 3;
            builder.add_polygon( polygon );

            b = builder.finalize();
        }

        EXPECT_TRUE( a->face_count() == b->face_count() );

        EXPECT_TRUE( !is_consistent_topology( a, b ) );
        EXPECT_TRUE( !is_consistent_topology( b, a ) );
    }
}

TEST( Polymesh3, CreateFromMeshInterface ) {
    frantic::geometry::trimesh3 trimesh;

    { // scope for populating trimesh
        trimesh.add_vertex( 0, 0, 0 );
        trimesh.add_vertex( 1, 0, 0 );
        trimesh.add_vertex( 1, 1, 0 );

        trimesh.add_face( 0, 1, 2 );

        trimesh.add_face_channel<boost::uint16_t>( _T("Face") );
        frantic::geometry::trimesh3_face_channel_accessor<boost::uint16_t> materialIdAcc =
            trimesh.get_face_channel_accessor<boost::uint16_t>( _T("Face") );
        materialIdAcc[0] = 1;

        trimesh.add_vertex_channel<float>( _T("SimpleVertex") );
        frantic::geometry::trimesh3_vertex_channel_accessor<float> simpleAcc =
            trimesh.get_vertex_channel_accessor<float>( _T("SimpleVertex") );
        simpleAcc[0] = 1;
        simpleAcc[1] = 2;
        simpleAcc[2] = 3;

        trimesh.add_vertex_channel<int>( _T("CustomVertex"), 2, true );
        frantic::geometry::trimesh3_vertex_channel_accessor<int> customAcc =
            trimesh.get_vertex_channel_accessor<int>( _T("CustomVertex") );
        customAcc[0] = 4;
        customAcc[1] = 5;

        customAcc.face( 0 ).set( 1, 0, 1 );
    }

    std::unique_ptr<frantic::geometry::trimesh3_interface> meshInterface(
        frantic::geometry::trimesh3_interface::create_instance( boost::move( trimesh ) ) );

    frantic::geometry::polymesh3_ptr polymesh = frantic::geometry::create_polymesh3( meshInterface.get() );
    EXPECT_TRUE( bool( polymesh ) );
    EXPECT_EQ( polymesh->vertex_count(), 3 );
    EXPECT_EQ( polymesh->face_count(), 1 );

    { // scope for geomAcc
        frantic::geometry::polymesh3_const_vertex_accessor<frantic::graphics::vector3f> geomAcc =
            polymesh->get_const_vertex_accessor<frantic::graphics::vector3f>( _T("verts") );
        EXPECT_EQ( geomAcc.get_vertex( 0 ), frantic::graphics::vector3f( 0 ) );
        EXPECT_EQ( geomAcc.get_vertex( 1 ), frantic::graphics::vector3f( 1, 0, 0 ) );
        EXPECT_EQ( geomAcc.get_vertex( 2 ), frantic::graphics::vector3f( 1, 1, 0 ) );

        frantic::geometry::polymesh3_const_vertex_accessor<frantic::graphics::vector3f>::const_face_range face =
            geomAcc.get_face( 0 );
        EXPECT_EQ( face.second - face.first, 3 );
        EXPECT_EQ( face.first[0], 0 );
        EXPECT_EQ( face.first[1], 1 );
        EXPECT_EQ( face.first[2], 2 );
    }

    { // scope for faceAcc
        EXPECT_TRUE( polymesh->has_face_channel( _T("Face") ) );
        frantic::geometry::polymesh3_const_face_accessor<boost::uint16_t> faceAcc =
            polymesh->get_const_face_accessor<boost::uint16_t>( _T("Face") );
        EXPECT_EQ( faceAcc.face_count(), 1 );
        EXPECT_EQ( faceAcc.get_face( 0 ), 1 );
    }

    { // scope for simpleAcc
        EXPECT_TRUE( polymesh->has_vertex_channel( _T("SimpleVertex") ) );
        frantic::geometry::polymesh3_const_vertex_accessor<float> simpleAcc =
            polymesh->get_const_vertex_accessor<float>( _T("SimpleVertex") );
        EXPECT_EQ( simpleAcc.vertex_count(), 3 );
        EXPECT_EQ( simpleAcc.get_vertex( 0 ), 1 );
        EXPECT_EQ( simpleAcc.get_vertex( 1 ), 2 );
        EXPECT_EQ( simpleAcc.get_vertex( 2 ), 3 );

        EXPECT_EQ( simpleAcc.face_count(), 0 );
    }

    { // scope for customAcc
        EXPECT_TRUE( polymesh->has_vertex_channel( _T("CustomVertex") ) );
        frantic::geometry::polymesh3_const_vertex_accessor<int> customAcc =
            polymesh->get_const_vertex_accessor<int>( _T("CustomVertex") );
        EXPECT_EQ( customAcc.vertex_count(), 2 );
        EXPECT_EQ( customAcc.get_vertex( 0 ), 4 );
        EXPECT_EQ( customAcc.get_vertex( 1 ), 5 );

        EXPECT_EQ( customAcc.face_count(), 1 );
        frantic::geometry::polymesh3_const_vertex_accessor<frantic::graphics::vector3f>::const_face_range customFace =
            customAcc.get_face( 0 );
        EXPECT_EQ( customFace.second - customFace.first, 3 );
        EXPECT_EQ( customFace.first[0], 1 );
        EXPECT_EQ( customFace.first[1], 0 );
        EXPECT_EQ( customFace.first[2], 1 );
    }
}

TEST( Polymesh3, EraseVertexChannel ) {
    using namespace frantic::geometry;

    polymesh3_ptr simpleMesh = make_cube_polymesh();

    frantic::tstring myChannelName = _T( "MyNormalChannel" );

    EXPECT_ANY_THROW( simpleMesh->erase_vertex_channel( myChannelName ) )
        << "Calling erase on non-existent channel did not throw.";

    simpleMesh->build_vertex_normals( myChannelName, true );

    EXPECT_TRUE( simpleMesh->has_vertex_channel( myChannelName ) ) << "Channel not created.";
    EXPECT_TRUE( !simpleMesh->has_face_channel( myChannelName ) )
        << "Misidentified vertex normals channel as face channel.";

    simpleMesh->erase_vertex_channel( myChannelName );

    EXPECT_TRUE( !simpleMesh->has_vertex_channel( myChannelName ) ) << "Channel not destroyed.";

    EXPECT_ANY_THROW( simpleMesh->erase_vertex_channel( myChannelName ) )
        << "Calling erase on non-existent channel did not throw.";
}

TEST( Polymesh3, BuildFaceNormalChannel ) {
    using namespace frantic::geometry;

    polymesh3_ptr simpleMesh = make_cube_polymesh();

    frantic::tstring myChannelName = _T( "MyNormalChannel" );

    simpleMesh->build_face_normals( myChannelName );

    EXPECT_TRUE( simpleMesh->has_face_channel( myChannelName ) ) << "Face normals channel was not created.";

    polymesh3_cvt_face_accessor<vector3f> normalsChannel = simpleMesh->get_cvt_face_accessor<vector3f>( myChannelName );

    for( size_t i = 0; i < simpleMesh->face_count(); ++i ) {
        vector3f storedNormal = normalsChannel.get_face( i );

        EXPECT_EQ( 1.0f, storedNormal.get_magnitude() ) << "Normal was not normalized.";

        // TODO: more tests might be nice
    }
}

TEST( Polymesh3, EraseFaceChannel ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_ptr simpleMesh = make_cube_polymesh();

    frantic::tstring myChannelName = _T( "MyNormalChannel" );

    EXPECT_ANY_THROW( simpleMesh->erase_face_channel( myChannelName ) )
        << "Calling erase on non-existent channel did not throw.";

    simpleMesh->build_face_normals( myChannelName );

    EXPECT_TRUE( simpleMesh->has_face_channel( myChannelName ) ) << "Channel not created.";
    EXPECT_TRUE( !simpleMesh->has_vertex_channel( myChannelName ) )
        << "Misidentified face normals channel as vertex channel.";

    simpleMesh->erase_face_channel( myChannelName );

    EXPECT_TRUE( !simpleMesh->has_face_channel( myChannelName ) ) << "Channel not destroyed.";

    EXPECT_ANY_THROW( simpleMesh->erase_face_channel( myChannelName ) )
        << "Calling erase on non-existent channel did not throw.";
}

TEST( Polymesh3, CullFaces ) {
    frantic::geometry::polymesh3_builder builder;

    builder.add_vertex( 0, 0, 0 );
    builder.add_vertex( 1, 0, 0 );
    builder.add_vertex( 1, 1, 0 );
    builder.add_vertex( 0, 1, 0 );

    std::vector<int> face;
    face = boost::assign::list_of( 0 )( 1 )( 2 ).convert_to_container<std::vector<int>>();
    builder.add_polygon( face );
    face = boost::assign::list_of( 0 )( 2 )( 3 ).convert_to_container<std::vector<int>>();
    builder.add_polygon( face );

    frantic::geometry::polymesh3_ptr original = builder.finalize();

    {
        frantic::graphics::raw_byte_buffer buffer;
        buffer.resize( 2 * sizeof( boost::int32_t ) );
        boost::int32_t* data = reinterpret_cast<boost::int32_t*>( buffer.ptr_at( 0 ) );
        data[0] = 0;
        data[1] = 1;
        original->add_face_channel( _T("Face"), frantic::channels::data_type_int32, 1, buffer );
    }

    {
        frantic::graphics::raw_byte_buffer buffer;
        buffer.resize( 4 * sizeof( float ) );
        float* data = reinterpret_cast<float*>( buffer.ptr_at( 0 ) );
        for( int i = 0; i < 4; ++i ) {
            data[i] = static_cast<float>( i );
        }
        original->add_vertex_channel( _T("SimpleVertex"), frantic::channels::data_type_float32, 1, buffer );
    }

    {
        frantic::graphics::raw_byte_buffer buffer;
        buffer.resize( 6 * sizeof( float ) );
        float* data = reinterpret_cast<float*>( buffer.ptr_at( 0 ) );
        for( int i = 0; i < 6; ++i ) {
            data[i] = static_cast<float>( i );
        }
        std::vector<int> faces;
        for( int i = 0; i < 6; ++i ) {
            faces.push_back( i );
        }
        original->add_vertex_channel( _T("CustomVertex"), frantic::channels::data_type_float32, 1, buffer, &faces );
    }

    std::vector<bool> keepFaces = boost::assign::list_of( false )( true );

    frantic::geometry::polymesh3_ptr culled = frantic::geometry::cull_faces( original, keepFaces );

    ASSERT_TRUE( culled != 0 );

    ASSERT_EQ( 1, culled->face_count() );
    ASSERT_EQ( 3, culled->vertex_count() );

    frantic::geometry::polymesh3_const_vertex_accessor<frantic::graphics::vector3f> geom =
        culled->get_const_vertex_accessor<frantic::graphics::vector3f>( _T("verts") );

    ASSERT_EQ( 3, geom.get_face_degree( 0 ) );
    EXPECT_EQ( 0, geom.get_face( 0 ).first[0] );
    EXPECT_EQ( 1, geom.get_face( 0 ).first[1] );
    EXPECT_EQ( 2, geom.get_face( 0 ).first[2] );

    EXPECT_EQ( frantic::graphics::vector3f( 0, 0, 0 ), culled->get_vertex( 0 ) );
    EXPECT_EQ( frantic::graphics::vector3f( 1, 1, 0 ), culled->get_vertex( 1 ) );
    EXPECT_EQ( frantic::graphics::vector3f( 0, 1, 0 ), culled->get_vertex( 2 ) );

    {
        frantic::geometry::polymesh3_const_face_accessor<boost::int32_t> acc =
            culled->get_const_face_accessor<boost::int32_t>( _T("Face") );
        ASSERT_EQ( 1, acc.face_count() );
        EXPECT_EQ( 1, acc.get_face( 0 ) );
    }

    {
        frantic::geometry::polymesh3_const_vertex_accessor<float> acc =
            culled->get_const_vertex_accessor<float>( _T("SimpleVertex") );
        EXPECT_FALSE( acc.has_custom_faces() );
        ASSERT_EQ( 3, acc.vertex_count() );
        EXPECT_EQ( 0, acc.get_vertex( 0 ) );
        EXPECT_EQ( 2, acc.get_vertex( 1 ) );
        EXPECT_EQ( 3, acc.get_vertex( 2 ) );
    }

    {
        frantic::geometry::polymesh3_const_vertex_accessor<float> acc =
            culled->get_const_vertex_accessor<float>( _T("CustomVertex") );
        ASSERT_TRUE( acc.has_custom_faces() );
        ASSERT_EQ( 1, acc.face_count() );
        ASSERT_EQ( 3, acc.get_face_degree( 0 ) );
        const int* face = acc.get_face( 0 ).first;
        EXPECT_EQ( 3, acc.get_vertex( face[0] ) );
        EXPECT_EQ( 4, acc.get_vertex( face[1] ) );
        EXPECT_EQ( 5, acc.get_vertex( face[2] ) );
    }
}

TEST( Polymesh3, ConstSharedPolymesh3Builder ) {
    using namespace frantic::channels;
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    // Make simple channel data that we can use to build a mesh.
    polymesh3_channel_data verts;
    {
        vector3f data[] = { vector3f( 0, 0, 0 ), vector3f( 1, 0, 0 ), vector3f( 0, 1, 0 ) };
        raw_byte_buffer geomBuffer( data, sizeof( data ) );
        verts = polymesh3_channel_data( geomBuffer, data_type_float32, 3 );
    }

    polymesh3_channel_faces faceIndices, faceEndIndices;
    {
        std::vector<int> polyIndices;
        polyIndices.push_back( 0 );
        polyIndices.push_back( 1 );
        polyIndices.push_back( 2 );

        std::vector<int> polyEndIndices;
        polyEndIndices.push_back( static_cast<int>( polyIndices.size() ) );

        faceIndices = polymesh3_channel_faces( polyIndices );
        faceEndIndices = polymesh3_channel_faces( polyEndIndices );
    }

    polymesh3_channel_data vertexChannel;
    {
        float data[] = { 0, 1, 2 };
        raw_byte_buffer vertexChannelBuffer( data, sizeof( data ) );
        vertexChannel = polymesh3_channel_data( vertexChannelBuffer, data_type_float32, 1 );
    }

    polymesh3_channel_data customVertexChannel;
    {
        float data[] = { 4, 8 };
        raw_byte_buffer customVertexChannelBuffer( data, sizeof( data ) );
        customVertexChannel = polymesh3_channel_data( customVertexChannelBuffer, data_type_float32, 1 );
    }

    polymesh3_channel_faces customVertexFaces;
    {
        std::vector<int> customVertexChannelIndices;
        customVertexChannelIndices.push_back( 1 );
        customVertexChannelIndices.push_back( 0 );
        customVertexChannelIndices.push_back( 1 );
        customVertexFaces = polymesh3_channel_faces( customVertexChannelIndices );
    }

    polymesh3_channel_data faceChannel;
    {
        float data[] = { 16 };
        raw_byte_buffer faceChannelBuffer( data, sizeof( data ) );
        faceChannel = polymesh3_channel_data( faceChannelBuffer, data_type_float32, 1 );
    }

    // None of the channel data should be shared yet.
    EXPECT_FALSE( verts.is_shared() );
    EXPECT_FALSE( faceIndices.is_shared() );
    EXPECT_FALSE( faceEndIndices.is_shared() );

    EXPECT_FALSE( vertexChannel.is_shared() );
    EXPECT_FALSE( customVertexChannel.is_shared() );
    EXPECT_FALSE( customVertexFaces.is_shared() );
    EXPECT_FALSE( faceChannel.is_shared() );

    // Build a mesh from the channel data.
    const_polymesh3_ptr meshA;
    {
        const_shared_polymesh3_builder builder( verts, faceIndices, faceEndIndices );

        builder.add_vertex_channel( _T("Vertex"), vertexChannel );
        builder.add_vertex_channel( _T("CustomVertex"), customVertexChannel, &customVertexFaces );
        builder.add_face_channel( _T("Face"), faceChannel );

        meshA = builder.finalize();
    }

    // Now all of the channel data should be shared.
    EXPECT_TRUE( verts.is_shared() );
    EXPECT_TRUE( faceIndices.is_shared() );
    EXPECT_TRUE( faceEndIndices.is_shared() );

    EXPECT_TRUE( vertexChannel.is_shared() );
    EXPECT_TRUE( customVertexChannel.is_shared() );
    EXPECT_TRUE( customVertexFaces.is_shared() );
    EXPECT_TRUE( faceChannel.is_shared() );

    // Build a second mesh from the channel data.
    const_polymesh3_ptr meshB;
    {
        const_shared_polymesh3_builder builder( verts, faceIndices, faceEndIndices );

        builder.add_vertex_channel( _T("Vertex"), vertexChannel );
        builder.add_vertex_channel( _T("CustomVertex"), customVertexChannel, &customVertexFaces );
        builder.add_face_channel( _T("Face"), faceChannel );

        meshB = builder.finalize();
    }

    // Make sure that meshA contains the expected data.
    EXPECT_EQ( meshA->vertex_count(), 3 );
    EXPECT_EQ( meshA->face_count(), 1 );

    {
        polymesh3_const_vertex_accessor<vector3f> geomAcc = meshA->get_const_vertex_accessor<vector3f>( _T("verts") );
        EXPECT_TRUE( geomAcc.has_custom_faces() );
        EXPECT_EQ( 3, geomAcc.vertex_count() );
        EXPECT_EQ( 1, geomAcc.face_count() );

        EXPECT_EQ( vector3f( 0, 0, 0 ), geomAcc.get_vertex( 0 ) );
        EXPECT_EQ( vector3f( 1, 0, 0 ), geomAcc.get_vertex( 1 ) );
        EXPECT_EQ( vector3f( 0, 1, 0 ), geomAcc.get_vertex( 2 ) );

        polymesh3_const_face_range face = geomAcc.get_face( 0 );
        EXPECT_EQ( 3, face.second - face.first );
        EXPECT_EQ( 0, face.first[0] );
        EXPECT_EQ( 1, face.first[1] );
        EXPECT_EQ( 2, face.first[2] );
    }

    {
        polymesh3_const_vertex_accessor<float> vertexAcc = meshA->get_const_vertex_accessor<float>( _T("Vertex") );
        EXPECT_FALSE( vertexAcc.has_custom_faces() );
        EXPECT_EQ( 3, vertexAcc.vertex_count() );
        EXPECT_EQ( 0, vertexAcc.get_vertex( 0 ) );
        EXPECT_EQ( 1, vertexAcc.get_vertex( 1 ) );
        EXPECT_EQ( 2, vertexAcc.get_vertex( 2 ) );
    }

    {
        polymesh3_const_vertex_accessor<float> customVertexAcc =
            meshA->get_const_vertex_accessor<float>( _T("CustomVertex") );
        EXPECT_TRUE( customVertexAcc.has_custom_faces() );
        EXPECT_EQ( 2, customVertexAcc.vertex_count() );
        ASSERT_EQ( 1, customVertexAcc.face_count() );
        EXPECT_EQ( 4, customVertexAcc.get_vertex( 0 ) );
        EXPECT_EQ( 8, customVertexAcc.get_vertex( 1 ) );
        polymesh3_const_face_range face = customVertexAcc.get_face( 0 );
        EXPECT_EQ( 3, face.second - face.first );
        EXPECT_EQ( 1, face.first[0] );
        EXPECT_EQ( 0, face.first[1] );
        EXPECT_EQ( 1, face.first[2] );
    }

    {
        polymesh3_const_face_accessor<float> faceAcc = meshA->get_const_face_accessor<float>( _T("Face") );
        EXPECT_EQ( 1, faceAcc.face_count() );
        EXPECT_EQ( 16, faceAcc.get_face( 0 ) );
    }

    // Make sure that the meshes and channel data buffers all use the exact
    // same memory.
    EXPECT_EQ( &meshA->get_channel_info( _T("verts" ) ).get_data(), &verts.get() );
    EXPECT_EQ( &meshA->get_channel_info( _T("verts" ) ).get_data(),
               &meshB->get_channel_info( _T("verts" ) ).get_data() );

    EXPECT_EQ( &meshA->get_channel_info( _T("verts" ) ).get_faces(), &faceIndices.get() );
    EXPECT_EQ( &meshA->get_channel_info( _T("verts" ) ).get_faces(),
               &meshB->get_channel_info( _T("verts" ) ).get_faces() );

    EXPECT_EQ( &meshA->get_channel_info( _T("Vertex" ) ).get_data(), &vertexChannel.get() );
    EXPECT_EQ( &meshA->get_channel_info( _T("Vertex" ) ).get_data(),
               &meshB->get_channel_info( _T("Vertex" ) ).get_data() );

    EXPECT_EQ( &meshA->get_channel_info( _T("CustomVertex" ) ).get_data(), &customVertexChannel.get() );
    EXPECT_EQ( &meshA->get_channel_info( _T("CustomVertex" ) ).get_data(),
               &meshB->get_channel_info( _T("CustomVertex" ) ).get_data() );

    EXPECT_EQ( &meshA->get_channel_info( _T("CustomVertex" ) ).get_faces(), &customVertexFaces.get() );
    EXPECT_EQ( &meshA->get_channel_info( _T("CustomVertex" ) ).get_faces(),
               &meshB->get_channel_info( _T("CustomVertex" ) ).get_faces() );

    EXPECT_EQ( &meshA->get_channel_info( _T("Face" ) ).get_data(), &faceChannel.get() );
    EXPECT_EQ( &meshA->get_channel_info( _T("Face" ) ).get_data(), &meshB->get_channel_info( _T("Face" ) ).get_data() );

    { // shouldn't be able to get write access to the shared data
        polymesh3* mesh = const_cast<polymesh3*>( meshA.get() );

        EXPECT_ANY_THROW( mesh->get_vertex_accessor<vector3f>( _T("verts") ) );
        EXPECT_ANY_THROW( mesh->get_vertex_accessor<float>( _T("Vertex") ) );
        EXPECT_ANY_THROW( mesh->get_vertex_accessor<float>( _T("CustomVertex") ) );
        EXPECT_ANY_THROW( mesh->get_face_accessor<float>( _T("Face") ) );
    }

    // After clearing the meshes, we expect the data to be unshared again.
    meshA.reset();
    meshB.reset();

    EXPECT_FALSE( verts.is_shared() );
    EXPECT_FALSE( faceIndices.is_shared() );
    EXPECT_FALSE( faceEndIndices.is_shared() );

    EXPECT_FALSE( vertexChannel.is_shared() );
    EXPECT_FALSE( customVertexChannel.is_shared() );
    EXPECT_FALSE( customVertexFaces.is_shared() );
    EXPECT_FALSE( faceChannel.is_shared() );
}

TEST( Polymesh3, FaceVertexCount ) {
    using namespace frantic::geometry;

    polymesh3_ptr mesh = make_triangle_polymesh();
    EXPECT_EQ( 3, mesh->face_vertex_count() );

    mesh = make_quad_polymesh();
    EXPECT_EQ( 4, mesh->face_vertex_count() );

    mesh = make_split_quad_polymesh();
    EXPECT_EQ( 6, mesh->face_vertex_count() );
}

TEST( Polymesh3, AddChannelsToEmptyMesh ) {
    using namespace frantic::channels;
    using namespace frantic::graphics;
    using namespace frantic::geometry;

    // Ensure that the add_*_channel() functions don't crash on an empty mesh.

    polymesh3_builder builder;
    polymesh3_ptr mesh = builder.finalize();

    raw_byte_buffer buffer;
    mesh->add_vertex_channel( _T("Velocity"), data_type_float32, 3, buffer );

    std::vector<int> faces;
    mesh->add_vertex_channel( _T("TextureCoord"), data_type_float32, 3, buffer, &faces );

    mesh->add_face_channel( _T("MaterialID"), data_type_uint16, 1, buffer );
}
