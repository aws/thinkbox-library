// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"
#include <frantic/geometry/particle_collision_detector.hpp>
#include <frantic/geometry/polygon_utils.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_file_io.hpp>
#include <frantic/geometry/xmesh_sequence_saver.hpp>
#include <frantic/volumetrics/levelset/rle_level_set_file_io.hpp>

using namespace std;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic;
using namespace frantic::diagnostics;
using namespace frantic::graphics2d;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::volumetrics::implicitsurface;

// Some utility methods for the xmesh sequence saver
static const frantic::tstring& test_temp_dir() {
    static frantic::tstring _static_test_temp_dir = _T( "./XMeshSequenceSaverTestScratchDirectory/" );
    return _static_test_temp_dir;
}

static void clean_temp_directory() {
    std::vector<frantic::tstring> fileListing;
    frantic::files::get_filenames_in_directory( test_temp_dir(), fileListing );

    for( size_t i = 0; i < fileListing.size(); ++i ) {
        frantic::files::delete_file( test_temp_dir() + fileListing[i] );
    }
}

static frantic::tstring get_verts_file( const frantic::tstring& baseName ) {
    return frantic::files::replace_extension( frantic::files::add_before_sequence_number( baseName, _T( "_verts" ) ),
                                              _T( ".xmdat" ) );
}

static frantic::tstring get_faces_file( const frantic::tstring& baseName ) {
    return frantic::files::replace_extension( frantic::files::add_before_sequence_number( baseName, _T( "_faces" ) ),
                                              _T( ".xmdat" ) );
}

static frantic::tstring get_channel_verts_file( const frantic::tstring& baseName,
                                                const frantic::tstring& channelName ) {
    return get_verts_file( frantic::files::add_before_sequence_number( baseName, _T( "_channel_" ) + channelName ) );
}

static frantic::tstring get_channel_faces_file( const frantic::tstring& baseName,
                                                const frantic::tstring& channelName ) {
    return get_faces_file( frantic::files::add_before_sequence_number( baseName, _T( "_channel_" ) + channelName ) );
}

static frantic::tstring get_facedata_file( const frantic::tstring& baseName, const frantic::tstring& channelName ) {
    return frantic::files::replace_extension(
        frantic::files::add_before_sequence_number( baseName, _T( "_channel_" ) + channelName + _T( "_facedata" ) ),
        _T( ".xmdat" ) );
}

static frantic::geometry::polymesh3_ptr make_polymesh( const std::vector<frantic::graphics::vector3f>& verts,
                                                       const std::vector<int>& indices,
                                                       const std::vector<int>& faceEnds ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    polymesh3_builder polymeshBuilder;

    for( size_t i = 0; i < verts.size(); ++i ) {
        polymeshBuilder.add_vertex( verts[i] );
    }

    size_t currentIndex = 0;
    for( size_t i = 0; i < faceEnds.size(); ++i ) {
        polymeshBuilder.add_polygon( &indices[currentIndex], faceEnds[i] - currentIndex );
        currentIndex = faceEnds[i];
    }

    return polymeshBuilder.finalize();
}

static bool compare_polymesh( frantic::geometry::polymesh3_ptr meshA, frantic::geometry::polymesh3_ptr meshB ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    std::vector<frantic::tstring> meshAVertexChannels;
    meshA->get_vertex_channel_names( meshAVertexChannels );
    meshAVertexChannels.push_back( _T( "verts" ) );

    std::vector<frantic::tstring> meshBVertexChannels;
    meshB->get_vertex_channel_names( meshBVertexChannels );
    meshBVertexChannels.push_back( _T( "verts" ) );

    if( meshAVertexChannels.size() != meshBVertexChannels.size() )
        return false;

    std::vector<frantic::tstring> meshAFaceChannels;
    meshA->get_face_channel_names( meshAFaceChannels );

    std::vector<frantic::tstring> meshBFaceChannels;
    meshB->get_face_channel_names( meshBFaceChannels );

    if( meshBFaceChannels.size() != meshAFaceChannels.size() )
        return false;

    std::vector<int> meshAFaceEndOffsets;
    meshA->get_face_end_offsets( meshAFaceEndOffsets );
    std::vector<int> meshBFaceEndOffsets;
    meshB->get_face_end_offsets( meshBFaceEndOffsets );

    if( meshAFaceEndOffsets != meshBFaceEndOffsets )
        return false;

    for( size_t i = 0; i < meshAVertexChannels.size(); ++i ) {

        if( !meshB->has_vertex_channel( meshAVertexChannels[i] ) )
            return false;

        const polymesh3_channel& meshAChannel = meshA->get_channel_info( meshAVertexChannels[i] );
        const polymesh3_channel& meshBChannel = meshB->get_channel_info( meshAVertexChannels[i] );

        if( meshAChannel.arity() != meshBChannel.arity() || meshAChannel.type() != meshBChannel.type() )
            return false;

        if( meshAChannel.get_data() != meshBChannel.get_data() )
            return false;

        if( meshAChannel.get_faces() != meshBChannel.get_faces() )
            return false;
    }

    for( size_t i = 0; i < meshAFaceChannels.size(); ++i ) {

        if( !meshB->has_face_channel( meshAFaceChannels[i] ) )
            return false;

        const polymesh3_channel& meshAChannel = meshA->get_channel_info( meshAFaceChannels[i] );
        const polymesh3_channel& meshBChannel = meshB->get_channel_info( meshAFaceChannels[i] );

        if( meshAChannel.arity() != meshBChannel.arity() || meshAChannel.type() != meshBChannel.type() )
            return false;

        if( meshAChannel.get_data() != meshBChannel.get_data() )
            return false;
    }

    return true;
}

// Google fixture for all tests
class GeometryTests : public ::testing::Test {
  public:
    virtual void SetUp() {
        if( frantic::files::directory_exists( test_temp_dir() ) ) {
            clean_temp_directory();
        } else {
            frantic::files::make_directory( test_temp_dir() );
            cout << "directory made" << endl;
        }
    }

    virtual void TearDown() { clean_temp_directory(); }
};

TEST_F( GeometryTests, ChangingTrimeshOutputToMeshSequenceSaver ) {
    trimesh3 simpleTrimesh;
    trimesh3 testTrimesh;

    make_quad_trimesh( simpleTrimesh );

    xmesh_sequence_saver testSaver;

    frantic::tstring rootFilename = _T( "atestxmeshfile" );

    frantic::files::filename_pattern filePattern( test_temp_dir() + rootFilename + _T( "0000.xmesh" ) );

    for( int i = 0; i < 2; ++i ) {
        // cout << frantic::strings::to_string(filePattern[i]) << endl;
        testSaver.write_xmesh( simpleTrimesh, filePattern[i] );

        load_xmesh_file( filePattern[i], testTrimesh );
        EXPECT_EQ( testTrimesh, simpleTrimesh );
    }

    EXPECT_TRUE( files::file_exists( get_verts_file( filePattern[0] ) ) )
        << "Vertex data file for frame 0 should exist";
    EXPECT_TRUE( files::file_exists( get_faces_file( filePattern[0] ) ) ) << "Face data file for frame 0 should exist";

    EXPECT_TRUE( !files::file_exists( get_verts_file( filePattern[1] ) ) )
        << "Vertex data file for frame 1 should not exist";
    EXPECT_TRUE( !files::file_exists( get_faces_file( filePattern[1] ) ) )
        << "Face data file for frame 1 should not exist";

    // change a vertex
    simpleTrimesh.get_vertex( 3 ) = vector3f( 4, 5, 6 );

    testSaver.write_xmesh( simpleTrimesh, filePattern[2] );

    load_xmesh_file( filePattern[2], testTrimesh );
    EXPECT_TRUE( testTrimesh == simpleTrimesh ) << "Loaded trimesh was not the same as the saved trimesh.";

    EXPECT_TRUE( files::file_exists( get_verts_file( filePattern[2] ) ) )
        << "Vertex data file for frame 2 should exist";
    EXPECT_TRUE( !files::file_exists( get_faces_file( filePattern[2] ) ) )
        << "Face data file for frame 2 should not exist";

    // now, add another face to cause a re-write of the faces file, but not the
    // vertices (this technically is no longer a value mesh, but w/e)
    simpleTrimesh.add_face( vector3( 1, 2, 3 ) );

    testSaver.write_xmesh( simpleTrimesh, filePattern[3] );

    load_xmesh_file( filePattern[3], testTrimesh );
    EXPECT_TRUE( testTrimesh == simpleTrimesh ) << "Loaded trimesh was not the same as the saved trimesh.";

    EXPECT_TRUE( !files::file_exists( get_verts_file( filePattern[3] ) ) )
        << "Vertex data file for frame 3 should exist";
    EXPECT_TRUE( files::file_exists( get_faces_file( filePattern[3] ) ) ) << "Face data file for frame 3 should exist";

    const frantic::tstring customChannelName = _T("SuperChannel");

    simpleTrimesh.add_vertex_channel<vector3f>( customChannelName, simpleTrimesh.vertex_count(), true );

    trimesh3_vertex_channel_accessor<vector3f> customChannel =
        simpleTrimesh.get_vertex_channel_accessor<vector3f>( customChannelName );

    for( size_t i = 0; i < customChannel.size(); ++i ) {
        customChannel[i] = vector3f( static_cast<float>( i ) );
    }

    testSaver.write_xmesh( simpleTrimesh, filePattern[4] );

    load_xmesh_file( filePattern[4], testTrimesh );
    EXPECT_TRUE( testTrimesh == simpleTrimesh ) << "Loaded trimesh was not the same as the saved trimesh.";

    EXPECT_TRUE( !files::file_exists( get_verts_file( filePattern[4] ) ) )
        << "Vertex data file for frame 4 should not exist";
    EXPECT_TRUE( !files::file_exists( get_faces_file( filePattern[4] ) ) )
        << "Face data file for frame 4 should not exist";
    EXPECT_TRUE( files::file_exists( get_channel_verts_file( filePattern[4], customChannelName ) ) )
        << "Vertex custom channel data file for frame 4 should exist";
    EXPECT_TRUE( files::file_exists( get_channel_faces_file( filePattern[4], customChannelName ) ) )
        << "Face custom channel data file for frame 4 should exist";

    testSaver.write_xmesh( simpleTrimesh, filePattern[5] );

    load_xmesh_file( filePattern[5], testTrimesh );
    EXPECT_TRUE( testTrimesh == simpleTrimesh ) << "Loaded trimesh was not the same as the saved trimesh.";

    EXPECT_TRUE( !files::file_exists( get_channel_verts_file( filePattern[5], customChannelName ) ) )
        << "Vertex custom channel data file for frame 5 should not exist";
    EXPECT_TRUE( !files::file_exists( get_channel_faces_file( filePattern[5], customChannelName ) ) )
        << "Face custom channel data file for frame 5 should not exist";
}

TEST_F( GeometryTests, ChangingPolymeshOutputToXMeshSequenceSaver ) {

    xmesh_sequence_saver testSaver;

    frantic::tstring rootFilename = _T("atestxmeshfile");

    frantic::files::filename_pattern filePattern( test_temp_dir() + rootFilename + _T("0000.xmesh") );

    std::vector<vector3f> vertices;
    std::vector<int> indices;
    std::vector<int> faceEnds;

    vertices.push_back( vector3f( 0.0, 0.0, 0.0 ) );
    vertices.push_back( vector3f( 0.0, 1.0, 0.0 ) );
    vertices.push_back( vector3f( 1.0, 1.0, 0.0 ) );
    vertices.push_back( vector3f( 1.0, 0.0, 0.0 ) );

    indices.push_back( 0 );
    indices.push_back( 1 );
    indices.push_back( 2 );
    indices.push_back( 3 );

    faceEnds.push_back( static_cast<int>( indices.size() ) );

    for( int i = 0; i < 2; ++i ) {
        polymesh3_ptr polymesh = make_polymesh( vertices, indices, faceEnds );
        testSaver.write_xmesh( polymesh, filePattern[i] );
        polymesh3_ptr testPolymesh = load_xmesh_polymesh_file( filePattern[i] );
        ASSERT_TRUE( compare_polymesh( testPolymesh, polymesh ) )
            << "Loaded polymesh was not the same as the saved polymesh.";
    }

    ASSERT_TRUE( files::file_exists( get_verts_file( filePattern[0] ) ) )
        << "Vertex data file for frame 0 should exist";
    ASSERT_TRUE( files::file_exists( get_faces_file( filePattern[0] ) ) ) << "Face data file for frame 0 should exist";
    ASSERT_TRUE( !files::file_exists( get_verts_file( filePattern[1] ) ) )
        << "Vertex data file for frame 1 should not exist";
    ASSERT_TRUE( !files::file_exists( get_faces_file( filePattern[1] ) ) )
        << "Face data file for frame 1 should not exist";

    vertices.push_back( vector3f( 4, 5, 6 ) );

    polymesh3_ptr frame2 = make_polymesh( vertices, indices, faceEnds );

    testSaver.write_xmesh( frame2, filePattern[2] );

    polymesh3_ptr frame2OutAgain = load_xmesh_polymesh_file( filePattern[2] );

    ASSERT_TRUE( compare_polymesh( frame2, frame2OutAgain ) )
        << "Loaded trimesh was not the same as the saved trimesh.";
    ASSERT_TRUE( files::file_exists( get_verts_file( filePattern[2] ) ) )
        << "Vertex data file for frame 2 should exist";
    ASSERT_TRUE( !files::file_exists( get_faces_file( filePattern[2] ) ) )
        << "Face data file for frame 2 should not exist";

    indices.push_back( 2 );
    indices.push_back( 3 );
    indices.push_back( 4 );
    faceEnds.push_back( static_cast<int>( indices.size() ) );

    polymesh3_ptr frame3 = make_polymesh( vertices, indices, faceEnds );

    testSaver.write_xmesh( frame3, filePattern[3] );

    polymesh3_ptr frame3OutAgain = load_xmesh_polymesh_file( filePattern[3] );

    ASSERT_TRUE( compare_polymesh( frame3, frame3OutAgain ) )
        << "Loaded trimesh was not the same as the saved trimesh.";

    ASSERT_TRUE( !files::file_exists( get_verts_file( filePattern[3] ) ) )
        << "Vertex data file for frame 3 should exist";
    ASSERT_TRUE( files::file_exists( get_faces_file( filePattern[3] ) ) ) << "Face data file for frame 3 should exist";

    const frantic::tstring customChannelName = _T("SuperChannel");

    polymesh3_ptr frame4 = make_polymesh( vertices, indices, faceEnds );

    frame4->add_empty_vertex_channel( customChannelName, frantic::channels::data_type_float32, 3, vertices.size() );

    polymesh3_vertex_accessor<vector3f> customChannel = frame4->get_vertex_accessor<vector3f>( customChannelName );

    for( size_t i = 0; i < customChannel.vertex_count(); ++i ) {
        customChannel.get_vertex( i ) = vector3f( static_cast<float>( i ) );
    }

    testSaver.write_xmesh( frame4, filePattern[4] );

    polymesh3_ptr frame4OutAgain = load_xmesh_polymesh_file( filePattern[4] );

    ASSERT_TRUE( compare_polymesh( frame4, frame4OutAgain ) )
        << "Loaded trimesh was not the same as the saved trimesh.";
    ASSERT_TRUE( !files::file_exists( get_verts_file( filePattern[4] ) ) )
        << "Vertex data file for frame 4 should not exist";
    ASSERT_TRUE( !files::file_exists( get_faces_file( filePattern[4] ) ) )
        << "Face data file for frame 4 should not exist";
    ASSERT_TRUE( files::file_exists( get_channel_verts_file( filePattern[4], customChannelName ) ) )
        << "Vertex custom channel data file for frame 4 should exist";
    ASSERT_TRUE( files::file_exists( get_channel_faces_file( filePattern[4], customChannelName ) ) )
        << "Face custom channel data file for frame 4 should exist";

    polymesh3_ptr frame5 = frame4;

    testSaver.write_xmesh( frame5, filePattern[5] );

    polymesh3_ptr frame5OutAgain = load_xmesh_polymesh_file( filePattern[5] );
    ASSERT_TRUE( compare_polymesh( frame5, frame5OutAgain ) )
        << "Loaded trimesh was not the same as the saved trimesh.";
    ASSERT_TRUE( !files::file_exists( get_channel_verts_file( filePattern[5], customChannelName ) ) )
        << "Vertex custom channel data file for frame 5 should not exist";

    ASSERT_TRUE( !files::file_exists( get_channel_faces_file( filePattern[5], customChannelName ) ) )
        << "Face custom channel data file for frame 5 should not exist";
}
