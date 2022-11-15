// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <tbb/task_scheduler_init.h>

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_degeneracy_removal.hpp>
#include <frantic/geometry/trimesh3_file_io.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>
#include <frantic/geometry/xmesh_sequence_saver.hpp>

#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/stream.hpp>

#include <frantic/files/filename_sequence.hpp>
#include <frantic/files/files.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>

#include <frantic/locale/locale.hpp>

#include <boost/cstdint.hpp>

#include "utilities/mesh_generators.hpp"

using frantic::files::scoped_file_cleanup;

TEST( Trimesh3, RemoveCoincidentGeometry ) {

    using namespace frantic::geometry;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;

    make_quad_trimesh( mesh );
    EXPECT_FALSE( remove_coincident_geometry( mesh ) );

    make_regular_tetrahedron( mesh );
    EXPECT_FALSE( remove_coincident_geometry( mesh ) );

    make_tetrahedron_beveled_edge_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( remove_coincident_geometry( mesh ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_EQ( 4, mesh.vertex_count() );
    EXPECT_EQ( 4, mesh.face_count() );

    EXPECT_FALSE( remove_coincident_geometry( mesh ) );

    make_triangle_beveled_vertex_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( remove_coincident_geometry( mesh ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_EQ( 4, mesh.vertex_count() );
    EXPECT_EQ( 3, mesh.face_count() );
    EXPECT_FALSE( remove_coincident_geometry( mesh ) );

    make_square_beveled_vertex_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( remove_coincident_geometry( mesh ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_EQ( 5, mesh.vertex_count() );
    EXPECT_EQ( 4, mesh.face_count() );
    EXPECT_FALSE( remove_coincident_geometry( mesh ) );

    make_chained_combined_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( remove_coincident_geometry( mesh ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_EQ( 5, mesh.vertex_count() );
    EXPECT_EQ( 4, mesh.face_count() );
    EXPECT_FALSE( remove_coincident_geometry( mesh ) );
}

TEST( Trimesh3, RemoveCoincidentGeometryWithChannels ) {
    using namespace frantic::geometry;

    trimesh3 mesh;

    mesh.add_vertex( 0, 0, 0 );
    mesh.add_vertex( 1, 0, 0 );
    mesh.add_vertex( 1, 1, 0 );
    mesh.add_vertex( 0, 1, 0 );
    mesh.add_vertex( 0, 0, 0 );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 2, 4, 0 ); // zero area face
    mesh.add_face( 2, 3, 4 );

    mesh.add_face_channel<boost::uint16_t>( _T("MaterialID") );
    {
        trimesh3_face_channel_accessor<boost::uint16_t> acc(
            mesh.get_face_channel_accessor<boost::uint16_t>( _T("MaterialID") ) );
        for( boost::uint16_t i = 0; i < 3; ++i ) {
            acc[i] = i;
        }
    }

    mesh.add_vertex_channel<float>( _T("Selection") );
    {
        trimesh3_vertex_channel_accessor<float> acc( mesh.get_vertex_channel_accessor<float>( _T("Selection") ) );
        for( int i = 0; i < 5; ++i ) {
            // % so the first and last vertices, which are both at the
            // origin, have the same value
            acc[i] = static_cast<float>( i % 4 );
        }
    }

    mesh.add_vertex_channel<frantic::graphics::vector3f>( _T("Color"), 9, true );
    {
        trimesh3_vertex_channel_accessor<frantic::graphics::vector3f> acc(
            mesh.get_vertex_channel_accessor<frantic::graphics::vector3f>( _T("Color") ) );
        // different color per face corner
        for( int i = 0; i < 9; ++i ) {
            acc[i].set( static_cast<float>( i ) );
        }
        acc.face( 0 ).set( 0, 1, 2 );
        acc.face( 1 ).set( 3, 4, 5 );
        acc.face( 2 ).set( 6, 7, 8 );
    }

    EXPECT_TRUE( remove_coincident_geometry( mesh ) );

    ASSERT_EQ( mesh.vertex_count(), 4 );
    ASSERT_EQ( mesh.face_count(), 2 );

    {
        trimesh3_face_channel_accessor<boost::uint16_t> acc(
            mesh.get_face_channel_accessor<boost::uint16_t>( _T("MaterialID") ) );
        ASSERT_EQ( acc.size(), 2 );
        EXPECT_EQ( acc[0], 0 );
        EXPECT_EQ( acc[1], 2 );
    }

    {
        trimesh3_vertex_channel_accessor<float> acc( mesh.get_vertex_channel_accessor<float>( _T("Selection") ) );
        ASSERT_EQ( acc.size(), 4 );
        EXPECT_EQ( acc[0], 0 );
        EXPECT_EQ( acc[1], 1 );
        EXPECT_EQ( acc[2], 2 );
        EXPECT_EQ( acc[3], 3 );
    }

    {
        trimesh3_vertex_channel_accessor<frantic::graphics::vector3f> acc(
            mesh.get_vertex_channel_accessor<frantic::graphics::vector3f>( _T("Color") ) );
        ASSERT_EQ( acc.size(), 6 );
        ASSERT_EQ( acc.face_count(), 2 );
        EXPECT_EQ( acc[acc.face( 0 )[0]], frantic::geometry::vector3f( 0 ) );
        EXPECT_EQ( acc[acc.face( 0 )[1]], frantic::geometry::vector3f( 1 ) );
        EXPECT_EQ( acc[acc.face( 0 )[2]], frantic::geometry::vector3f( 2 ) );
        EXPECT_EQ( acc[acc.face( 1 )[0]], frantic::geometry::vector3f( 6 ) );
        EXPECT_EQ( acc[acc.face( 1 )[1]], frantic::geometry::vector3f( 7 ) );
        EXPECT_EQ( acc[acc.face( 1 )[2]], frantic::geometry::vector3f( 8 ) );
    }
}

TEST( Trimesh3, CheckDuplicateFaceIndices ) {
    using namespace frantic::geometry;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    EXPECT_TRUE( mesh.check_duplicate_face_indices( nullStream ) );

    mesh.add_face( 0, 0, 1 );

    EXPECT_FALSE( mesh.check_duplicate_face_indices( nullStream ) );

    remove_duplicate_face_indices( mesh );

    EXPECT_TRUE( mesh.check_duplicate_face_indices( nullStream ) );
}

TEST( Trimesh3, CheckZeroAreaFaces ) {
    using namespace frantic::geometry;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    EXPECT_TRUE( mesh.check_duplicate_face_indices( nullStream ) );

    mesh.clear();

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    mesh.add_vertex( 0.0f, 0.0f, 1.0f );
    mesh.add_vertex( 0.0f, 1.0f, 0.0f );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 3, 4 );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
}

TEST( Trimesh3, CheckFiniteVertices ) {
    using namespace frantic::geometry;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    EXPECT_TRUE( mesh.check_finite_vertices( nullStream ) );

    mesh.add_vertex( 0.0f, 0.0f, std::numeric_limits<float>::infinity() );

    mesh.add_face( 0, 1, int( mesh.vertex_count() - 1 ) );

    EXPECT_FALSE( mesh.check_finite_vertices( nullStream ) );
}

TEST( Trimesh3, CheckIndexRanges ) {
    using namespace frantic::geometry;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    EXPECT_TRUE( mesh.check_index_ranges( nullStream ) );

    mesh.add_face( 0, 1, int( mesh.vertex_count() ) );

    EXPECT_FALSE( mesh.check_index_ranges( nullStream ) );
}

TEST( Trimesh3, LoadObjMeshFileWithGermanLocale ) {
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

    boost::filesystem::path objPath = boost::filesystem::unique_path(
        ( boost::filesystem::temp_directory_path() / L"%%%%-%%%%-%%%%-%%%%.csv" ).wstring() );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( objPath );

    // make sure we can read numbers that include a decimal point
    frantic::files::file_ptr f( frantic::files::tfopen( frantic::files::to_tstring( objPath ).c_str(), _T("w+") ) );
    std::fputs( "v 0 0 0\n", f );
    std::fputs( "v 1.2 0 0\n", f );
    std::fputs( "v 1.2 1.2 0\n", f );
    std::fputs( "f 1 2 3\n", f );
    f.close();

    frantic::geometry::trimesh3 mesh;
    frantic::geometry::load_obj_mesh_file( frantic::files::to_tstring( objPath ), mesh );

    EXPECT_EQ( 1, mesh.face_count() );
    EXPECT_EQ( 3, mesh.vertex_count() );
    EXPECT_EQ( frantic::graphics::vector3f( 0, 0, 0 ), mesh.get_vertex( 0 ) );
    EXPECT_EQ( frantic::graphics::vector3f( 1.2f, 0, 0 ), mesh.get_vertex( 1 ) );
    EXPECT_EQ( frantic::graphics::vector3f( 1.2f, 1.2f, 0 ), mesh.get_vertex( 2 ) );
}

void testTrimesh3ResizeChannels() {

    using frantic::geometry::trimesh3;
    using frantic::graphics::vector3f;

    trimesh3 mesh;

    //
    // vertex channels
    //
    mesh.add_vertex_channel<boost::uint8_t>( _T("byte") );
    mesh.add_vertex_channel<vector3f>( _T("vector3f") );

    EXPECT_EQ( 0, mesh.get_vertex_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 0, mesh.get_vertex_channel_general_accessor( _T("vector3f") ).size() );

    // add no vertices
    mesh.get_vertex_channel_general_accessor( _T("byte") ).add_vertices( 0 );
    mesh.get_vertex_channel_general_accessor( _T("vector3f") ).add_vertices( 0 );
    EXPECT_EQ( 0, mesh.get_vertex_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 0, mesh.get_vertex_channel_general_accessor( _T("vector3f") ).size() );

    // add 1 vertex
    mesh.add_vertex( 0, 0, 0 );
    mesh.get_vertex_channel_general_accessor( _T("byte") ).add_vertices( 1 );
    mesh.get_vertex_channel_general_accessor( _T("vector3f") ).add_vertices( 1 );
    EXPECT_EQ( 1, mesh.vertex_count() );
    EXPECT_EQ( 1, mesh.get_vertex_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 1, mesh.get_vertex_channel_general_accessor( _T("vector3f") ).size() );

    // add 2 vertices ( 1 + 2 = 3 )
    mesh.add_vertices( 2 );
    mesh.get_vertex_channel_general_accessor( _T("byte") ).add_vertices( 2 );
    mesh.get_vertex_channel_general_accessor( _T("vector3f") ).add_vertices( 2 );
    EXPECT_EQ( 3, mesh.vertex_count() );
    EXPECT_EQ( 3, mesh.get_vertex_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 3, mesh.get_vertex_channel_general_accessor( _T("vector3f") ).size() );

    //
    // face channels
    //
    mesh.add_face_channel<boost::uint8_t>( _T("byte") );
    mesh.add_face_channel<vector3f>( _T("vector3f") );

    EXPECT_EQ( 0, mesh.get_face_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 0, mesh.get_face_channel_general_accessor( _T("vector3f") ).size() );

    // add no faces
    mesh.get_face_channel_general_accessor( _T("byte") ).add_faces( 0 );
    mesh.get_face_channel_general_accessor( _T("vector3f") ).add_faces( 0 );
    EXPECT_EQ( 0, mesh.get_face_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 0, mesh.get_face_channel_general_accessor( _T("vector3f") ).size() );

    // add one face
    mesh.add_face( 0, 0, 0 );
    mesh.get_face_channel_general_accessor( _T("byte") ).add_faces( 1 );
    mesh.get_face_channel_general_accessor( _T("vector3f") ).add_faces( 1 );
    EXPECT_EQ( 1, mesh.get_face_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 1, mesh.get_face_channel_general_accessor( _T("vector3f") ).size() );

    // add two faces
    mesh.add_face( 0, 0, 1 );
    mesh.add_face( 0, 1, 1 );
    mesh.get_face_channel_general_accessor( _T("byte") ).add_faces( 2 );
    mesh.get_face_channel_general_accessor( _T("vector3f") ).add_faces( 2 );
    EXPECT_EQ( 3, mesh.get_face_channel_general_accessor( _T("byte") ).size() );
    EXPECT_EQ( 3, mesh.get_face_channel_general_accessor( _T("vector3f") ).size() );
}

TEST( Trimesh3, CopyToMeshInterface ) {
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

    frantic::geometry::trimesh3 outMesh;
    frantic::geometry::copy_to_trimesh3( meshInterface.get(), outMesh );

    EXPECT_EQ( outMesh, meshInterface->get_trimesh() );
}

// TODO: Perhaps these belong in a separate file for the kdtree only?
TEST( Trimesh3, Trimesh3KDTree ) {

    using namespace std;
    using namespace frantic::graphics;
    using namespace frantic::geometry;

    trimesh3 mesh;
    boundbox3f box( vector3f( 0 ), vector3f( 1 ) );
    mesh.set_to_box( box );
    trimesh3_kdtree kd( mesh );

    boundbox3f biggerBox( vector3f( -1 ), vector3f( 2 ) );
    // Test the intersects_ray_segment function
    for( int test = 0; test < 300; ++test ) {
        vector3f a = biggerBox.random_vector(), b = biggerBox.random_vector();
        ray3f ray( a, ( b - a ).to_normalized() );
        float distance = vector3f::distance( a, b );

        bool boundBoxIsect = ray.is_intersecting_box_surface( box, 0, distance );
        bool kdIsect = kd.intersects_ray_segment( a, b );

        if( boundBoxIsect != kdIsect || ( ray.is_intersecting_box_volume( box, 0, distance ) !=
                                          box.is_volume_intersecting_line_segment( a, b ) ) ) {
            cerr << endl << "ray: " << ray << ", distance: " << distance << endl;
            double i1, i2;
            vector3f n1, n2;
            if( ray.intersect_with_box( box, i1, i2, &n1, &n2 ) ) {
                cerr << "intersected boundbox at " << i1 << " and " << i2 << ", normals " << n1 << " and " << n2
                     << endl;
                cerr << "positions are " << ray.at( i1 ) << " and " << ray.at( i2 ) << endl;
            }
        }

        EXPECT_EQ( boundBoxIsect, kdIsect );
        EXPECT_EQ( ray.is_intersecting_box_volume( box, 0, distance ),
                   box.is_volume_intersecting_line_segment( a, b ) );
    }

    // Test the intersect_ray function
    for( int test = 0; test < 300; ++test ) {
        vector3f a = biggerBox.random_vector(), b = biggerBox.random_vector();
        ray3f ray( a, ( b - a ).to_normalized() );
        float distance = vector3f::distance( a, b );

        double i1, i2;
        vector3f n1, n2;
        raytrace_intersection isect;
        bool boundBoxIsect = ray.intersect_with_box( box, i1, i2, &n1, &n2 );
        if( i1 < 0 ) {
            i1 = i2;
            i2 = 0;
            n1 = n2;
            n2 = vector3f( 0 );
            if( i1 < 0 ) {
                i1 = 0;
                n1 = vector3f( 0 );
                boundBoxIsect = false;
            }
        }
        if( i2 > distance ) {
            i2 = 0;
            n2 = vector3f( 0 );
            if( i1 > distance ) {
                i1 = 0;
                n1 = vector3f( 0 );
                boundBoxIsect = false;
            }
        }
        if( i1 > distance ) {
            i1 = i2;
            n1 = n2;
            i2 = 0;
            n2 = vector3f( 0 );
            if( i1 == 0 )
                boundBoxIsect = false;
        }
        bool kdIsect = kd.intersect_ray( ray, 0, distance, isect );

        if( boundBoxIsect != kdIsect ||
            ( kdIsect && boundBoxIsect &&
              ( vector3f::distance_squared( ray.at( i1 ), isect.position ) > 0.000000001f ||
                vector3f::distance_squared( n1, isect.geometricNormal ) > 0.000000001f ) ) ) {
            cerr << endl << "ray: " << ray << ", distance: " << distance << endl;
            if( boundBoxIsect ) {
                cerr << "intersected boundbox at " << i1 << " and " << i2 << ", normals " << n1 << " and " << n2
                     << endl;
                cerr << "positions are " << ray.at( i1 ) << " and " << ray.at( i2 ) << endl;
                cerr << "normals are " << n1 << " and " << n2 << endl;
            }
            if( kdIsect ) {
                cerr << "intersected kdtree at " << isect.distance << endl;
                cerr << "position is " << isect.position << endl;
                cerr << "normal is " << isect.geometricNormal << endl;
            }
        }

        EXPECT_EQ( boundBoxIsect, kdIsect );
        if( kdIsect && boundBoxIsect ) {
            EXPECT_TRUE( vector3f::distance_squared( ray.at( i1 ), isect.position ) <= 0.000000001f );
            EXPECT_TRUE( vector3f::distance_squared( n1, isect.geometricNormal ) <= 0.000000001f );
        }
    }

    // Test the intersect_ray_all function
    for( int test = 0; test < 300; ++test ) {
        vector3f a = biggerBox.random_vector(), b = biggerBox.random_vector();
        ray3f ray( a, ( b - a ).to_normalized() );
        float distance = vector3f::distance( a, b );

        double i1, i2;
        vector3f n1, n2;
        int boundBoxIsectCount = 0;
        if( ray.intersect_with_box( box, i1, i2, &n1, &n2 ) ) {
            boundBoxIsectCount = 2;
            // Test i2 > distance before i1 < 0, so we are sure i2 is not zeroed
            if( i2 > distance ) {
                i2 = -1;
                n2 = vector3f( 0 );
                boundBoxIsectCount = 1;
                if( i1 > distance ) {
                    i1 = -1;
                    n1 = vector3f( 0 );
                    boundBoxIsectCount = 0;
                }
            }
            if( i1 < 0 ) {
                i1 = i2;
                i2 = -1;
                n1 = n2;
                n2 = vector3f( 0 );
                boundBoxIsectCount = 1;
                if( i1 < 0 ) {
                    i1 = -1;
                    n1 = vector3f( 0 );
                    boundBoxIsectCount = 0;
                }
            }
        }
        std::vector<raytrace_intersection> isects;
        kd.intersect_ray_all( ray, distance, isects );

        if( (int)isects.size() != boundBoxIsectCount ) {
            cerr << endl << "ray: " << ray << ", distance: " << distance << endl;
            if( boundBoxIsectCount > 0 ) {
                cerr << "intersected boundbox at " << i1 << " and " << i2 << ", normals " << n1 << " and " << n2
                     << endl;
                cerr << "positions are " << ray.at( i1 ) << " and " << ray.at( i2 ) << endl;
                cerr << "normals are " << n1 << " and " << n2 << endl;
            }
            for( unsigned i = 0; i < isects.size(); ++i ) {
                cerr << "intersected kdtree at " << isects[i].distance << endl;
                cerr << "position is " << isects[i].position << endl;
                cerr << "normal is " << isects[i].geometricNormal << endl;
            }
        } else {
            if( boundBoxIsectCount > 0 ) {
                if( vector3f::distance_squared( ray.at( i1 ), isects[0].position ) > 0.000000001f ||
                    vector3f::distance_squared( n1, isects[0].geometricNormal ) > 0.000000001f ) {
                    cerr << "intersected boundbox at " << i1 << " and " << i2 << ", normals " << n1 << " and " << n2
                         << endl;
                    cerr << "positions are " << ray.at( i1 ) << " and " << ray.at( i2 ) << endl;
                    cerr << "normals are " << n1 << " and " << n2 << endl;
                    cerr << "(first) intersected kdtree at " << isects[0].distance << endl;
                    cerr << "position is " << isects[0].position << endl;
                    cerr << "normal is " << isects[0].geometricNormal << endl;
                }
                EXPECT_TRUE( vector3f::distance_squared( ray.at( i1 ), isects[0].position ) <= 0.000000001f );
                EXPECT_TRUE( vector3f::distance_squared( n1, isects[0].geometricNormal ) <= 0.000000001f );
            }
            if( boundBoxIsectCount > 1 ) {
                if( vector3f::distance_squared( ray.at( i2 ), isects[1].position ) > 0.000000001f ||
                    vector3f::distance_squared( n2, isects[1].geometricNormal ) > 0.000000001f ) {
                    cerr << "intersected boundbox at " << i1 << " and " << i2 << ", normals " << n1 << " and " << n2
                         << endl;
                    cerr << "positions are " << ray.at( i1 ) << " and " << ray.at( i2 ) << endl;
                    cerr << "normals are " << n1 << " and " << n2 << endl;
                    cerr << "(second) intersected kdtree at " << isects[0].distance << endl;
                    cerr << "position is " << isects[1].position << endl;
                    cerr << "normal is " << isects[1].geometricNormal << endl;
                }
                EXPECT_TRUE( vector3f::distance_squared( ray.at( i2 ), isects[1].position ) <= 0.000000001f );
                EXPECT_TRUE( vector3f::distance_squared( n2, isects[1].geometricNormal ) <= 0.000000001f );
            }
        }

        EXPECT_EQ( int( isects.size() ), boundBoxIsectCount );
    }

    // Test the find_nearest_point function

    for( int test = 0; test < 300; ++test ) {
        vector3f a = biggerBox.random_vector();

        vector3f boundBoxP = box.get_nearest_surface_point( a );
        nearest_point_search_result np;
        EXPECT_TRUE( kd.find_nearest_point( a, 50, np ) );

        if( vector3f::distance_squared( boundBoxP, np.position ) > 0.0000000001f ) {
            cerr << vector3f::distance_squared( boundBoxP, np.position ) << endl;
            cerr << endl << "test point: " << a << endl;
            cerr << "bounding box point: " << boundBoxP << endl;
            cerr << "kd point: " << np.position << endl;
        }

        EXPECT_TRUE( vector3f::distance_squared( boundBoxP, np.position ) <= 0.0000000001f );
    }

    { // Test triangle clipping for a specific case that was killing the kdtree
        unsigned int min[] = { 0xc1827e50, 0x359929f5, 0x4208ee4b };
        unsigned int max[] = { 0xc0324bff, 0x41200069, 0x422c3e09 };
        boundbox3f box( *reinterpret_cast<vector3f*>( min ), *reinterpret_cast<vector3f*>( max ) );

        unsigned int v0[] = { 0xc0dd4ea8, 0x00000000, 0x422c3def };
        unsigned int v1[] = { 0xc0324da2, 0x35eaac7a, 0x41ec4743 };
        unsigned int v2[] = { 0xc139d1fb, 0x41200000, 0x4208ee4b };
        vector3f vec0 = *reinterpret_cast<vector3f*>( v0 );
        vector3f vec1 = *reinterpret_cast<vector3f*>( v1 );
        vector3f vec2 = *reinterpret_cast<vector3f*>( v2 );

        cerr << endl;
        cerr << "Boundbox: " << box << endl;
        cerr << "Triangle: " << endl << vec0 << endl << vec1 << endl << vec2 << endl;

        boundbox3f boxOut;
        EXPECT_TRUE( box.intersect_with_triangle( vec0, vec1, vec2, boxOut ) );
        EXPECT_TRUE( boxOut.maximum().x > vec0.x );
    }
}

TEST( Trimesh3, KDTreeCollectTriangles ) {
    using frantic::geometry::nearest_point_on_triangle;
    using frantic::geometry::trimesh3;
    using frantic::geometry::trimesh3_kdtree;
    using frantic::graphics::vector3;
    using frantic::graphics::vector3f;

    using namespace frantic::logging;

    // Temporarily logging debug to help track down an occasional test failure
    set_logging_level_in_scope logDebug( level::debug );

    //{	//A few problem cases to try specifically
    vector3f v[5];
    v[0] = vector3f( 6.90481f, 49.6185f, -45.5535f );
    v[1] = vector3f( -1.11545f, 30.4804f, -41.8394f );
    v[2] = vector3f( -0.642414f, -33.4986f, -3.2487f );
    v[3] = vector3f( 2.77184f, 68.0006f, -73.0178f );
    v[4] = nearest_point_on_triangle( v[3], v[0], v[1], v[2] );

    EXPECT_TRUE( vector3f::distance( v[3], v[4] ) > 26.87f );
    //}

    trimesh3 mesh;

    unsigned seed = 12345678u;

    FF_LOG( debug ) << "Random seed: " << seed << std::endl;

    srand( seed );
    for( int i = 0; i < 1000; ++i ) {
        mesh.add_face( 3 * i, 3 * i + 1, 3 * i + 2 );
        mesh.add_vertex( 100.f * vector3f::from_random() - vector3f( 50.f ) );
        mesh.add_vertex( 100.f * vector3f::from_random() - vector3f( 50.f ) );
        mesh.add_vertex( 100.f * vector3f::from_random() - vector3f( 50.f ) );
    }

    // write_mesh_file("c:\\dharrison\\temp_0000.xmesh", mesh);
    // std::ofstream fout( "c:\\dharrison\\temp_selection.ms" );
    // std::cout << "\nSeed: " << seed << "\n";

    trimesh3_kdtree tree( mesh );
    for( int i = 0; i < 1000; ++i ) {
        vector3f testPt = 150.f * vector3f::from_random() - vector3f( 75.f );
        float testRadius = 50.f * ( rand() / float( RAND_MAX ) );

        std::vector<int> outFaces;
        std::vector<int> outNaiveFaces;

        // Naively trying each triangle, find the faces with a point within testRadius of testPoint
        for( size_t j = 0; j < mesh.face_count(); ++j ) {
            vector3 f = mesh.get_face( j );
            vector3f a = mesh.get_vertex( f.x );
            vector3f b = mesh.get_vertex( f.y );
            vector3f c = mesh.get_vertex( f.z );
            if( vector3f::distance( nearest_point_on_triangle( testPt, a, b, c ), testPt ) < testRadius )
                outNaiveFaces.push_back( (int)j );
        }

        // Using the kdtree, collect the faces with a point within testRadius of testPoint
        tree.collect_faces_within_range( testPt, testRadius, outFaces );

        if( outFaces.size() != outNaiveFaces.size() ||
            !std::equal( outFaces.begin(), outFaces.end(), outNaiveFaces.begin() ) ) {
            // This will generate a script that will illuminate the issue if run within Max.

            EXPECT_EQ( outFaces.size(), outNaiveFaces.size() );
            EXPECT_TRUE( std::equal( outFaces.begin(), outFaces.end(), outNaiveFaces.begin() ) );
            return;
        }
        // std::cout << "\r" << i;
    }
}

struct nsprSorter {
    bool operator()( const frantic::geometry::nearest_point_search_result& lhs,
                     const frantic::geometry::nearest_point_search_result& rhs ) const {
        return lhs.distance == rhs.distance ? lhs.faceIndex < rhs.faceIndex : lhs.distance < rhs.distance;
    }
};

struct nsprEquals {
    bool operator()( const frantic::geometry::nearest_point_search_result& lhs,
                     const frantic::geometry::nearest_point_search_result& rhs ) const {
        return lhs.faceIndex == rhs.faceIndex;
    }
};

TEST( Trimesh3, KDTreeClosestTriangles ) {
    using namespace frantic::graphics;
    using namespace frantic::geometry;

    unsigned int randomSeed = 1458282101u;
    std::cerr << "Seed: " << randomSeed << std::endl;
    srand( randomSeed );

    trimesh3 mesh;

    for( int i = 0; i < 1000; ++i ) {
        mesh.add_vertex( vector3f::from_random() * 100.f - vector3f( 50.f ) );
        mesh.add_vertex( vector3f::from_random() * 100.f - vector3f( 50.f ) );
        mesh.add_vertex( vector3f::from_random() * 100.f - vector3f( 50.f ) );
        mesh.add_face( 3 * i, 3 * i + 1, 3 * i + 2 );
    }

    trimesh3_kdtree tree( mesh );
    for( std::size_t i = 0; i < 1000; ++i ) {
        vector3f testPt = vector3f::from_random() * 100.f - vector3f( 50.f );

        std::vector<nearest_point_search_result> results;
        std::vector<nearest_point_search_result> naiveResults;

        tree.collect_nearest_faces( testPt, 10, results );

        for( std::size_t j = 0; j < mesh.face_count(); ++j ) {
            vector3 face = mesh.get_face( j );
            vector3f v0 = mesh.get_vertex( face.x );
            vector3f v1 = mesh.get_vertex( face.y );
            vector3f v2 = mesh.get_vertex( face.z );

            nearest_point_search_result nspr;
            frantic::geometry::nearest_point_on_triangle( testPt, v0, v1, v2, nspr );
            nspr.faceIndex = (int)j;

            naiveResults.push_back( nspr );
            std::sort( naiveResults.begin(), naiveResults.end(), nsprSorter() );
            if( naiveResults.size() > 10 )
                naiveResults.pop_back();
        }

        if( results.size() != naiveResults.size() ||
            !std::equal( results.begin(), results.end(), naiveResults.begin(), nsprEquals() ) ) {
            std::cerr << "On pass: " << i << " of 1000" << std::endl;
            std::cerr << "KDtree found: " << results.size() << " faces" << std::endl;
            std::cerr << "Naive found: " << naiveResults.size() << " faces" << std::endl;

            for( std::vector<nearest_point_search_result>::iterator it = results.begin(), itEnd = results.end();
                 it != itEnd; ++it )
                std::cerr << it->faceIndex << " ";
            std::cerr << std::endl;
            for( std::vector<nearest_point_search_result>::iterator it = naiveResults.begin(),
                                                                    itEnd = naiveResults.end();
                 it != itEnd; ++it )
                std::cerr << it->faceIndex << " ";
            std::cerr << std::endl;

            std::sort( results.begin(), results.end() );
            std::sort( naiveResults.begin(), naiveResults.end() );

            EXPECT_TRUE( results.size() == naiveResults.size() );
            EXPECT_TRUE( std::equal( results.begin(), results.end(), naiveResults.begin(), nsprEquals() ) );
            return;
        }
    } // for(std::size_t i = 0; i < 1000; ++i)

    std::clog << "testKDTreeClosestTriangles() completed." << std::endl;
}

TEST( Trimesh3, RemoveDoubledFaces ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;

    make_regular_tetrahedron( mesh );
    mesh.add_face( mesh.get_face( 0 ) );

    EXPECT_FALSE( mesh.check_consistency( nullStream, trimesh3::duplicate_edges ) );

    remove_doubled_faces( mesh );

    EXPECT_EQ( 4, mesh.vertex_count() );
    EXPECT_EQ( 3, mesh.face_count() );
    EXPECT_TRUE( mesh.check_consistency( nullStream, trimesh3::duplicate_edges ) );

    make_regular_tetrahedron( mesh );
    vector3 face( mesh.get_face( 0 ) );
    std::swap( face[0], face[2] );
    mesh.add_face( face );

    EXPECT_FALSE( mesh.check_consistency( nullStream, trimesh3::duplicate_edges ) );

    remove_doubled_faces( mesh );

    EXPECT_EQ( 4, mesh.vertex_count() );
    EXPECT_EQ( 3, mesh.face_count() );
    EXPECT_TRUE( mesh.check_consistency( nullStream, trimesh3::duplicate_edges ) );

    // TODO: add tests for named channels
}

TEST( Trimesh3, RemoveDeadVertices ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    trimesh3 mesh;

    mesh.add_vertex( 0.0f, 0.0f, 0.0f ); // 0
    mesh.add_vertex( 0.0f, 0.0f, 1.0f ); // 1
    mesh.add_vertex( 0.0f, 1.0f, 0.0f ); // 2
    mesh.add_vertex( 0.0f, 1.0f, 1.0f ); // 3
    mesh.add_vertex( 1.0f, 0.0f, 0.0f ); // 4
    mesh.add_vertex( 1.0f, 0.0f, 1.0f ); // 5
    mesh.add_vertex( 1.0f, 1.0f, 0.0f ); // 6
    mesh.add_vertex( 1.0f, 1.0f, 1.0f ); // 7

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 5, 6, 7 );

    remove_dead_vertices( mesh );

    EXPECT_EQ( 6, mesh.vertex_count() );
    // This check is a little brittle, since it relies on vertex removal
    // being stable w.r.t. the original order, a general 'vertex-in-face' check would be better
    EXPECT_EQ( 3, mesh.get_face( 1 )[0] );
    EXPECT_EQ( 4, mesh.get_face( 1 )[1] );
    EXPECT_EQ( 5, mesh.get_face( 1 )[2] );
}

TEST( Trimesh3, RemoveDuplicateHalfedges ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    tbb::task_scheduler_init taskScheduleInit;

    trimesh3 mesh;

    make_sharkfin_mesh( mesh );

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    EXPECT_FALSE( mesh.check_duplicate_edges( nullStream ) );

    EXPECT_FALSE( remove_duplicate_halfedge_faces( mesh ) );

    EXPECT_TRUE( remove_duplicate_halfedge_faces( mesh, false ) );

    EXPECT_TRUE( mesh.check_duplicate_edges( nullStream ) );
}

TEST( Trimesh3, RemoveFrostMarchingCubeDegeneracy ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    tbb::task_scheduler_init taskScheduleInit;

    trimesh3 mesh;

    make_bubble_straw_mesh( mesh );

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    EXPECT_FALSE( mesh.check_duplicate_edges( nullStream ) );

    EXPECT_TRUE( remove_duplicate_halfedge_faces( mesh ) );

    EXPECT_TRUE( mesh.check_duplicate_edges( nullStream ) );
}

TEST( Trimesh3, XMeshSequenceSaver ) {
    using namespace std;
    using namespace frantic::geometry;
    using namespace frantic::files;
    namespace fs = boost::filesystem;

    trimesh3 mesh;

    mesh.add_vertex( vector3f( -1.0, -1.0, 0 ) );
    mesh.add_vertex( vector3f( 1.0, -1.0, 0 ) );
    mesh.add_vertex( vector3f( -1.0, 1.0, 0 ) );
    mesh.add_vertex( vector3f( 1.0, 1.0, 0 ) );
    mesh.add_vertex( vector3f( -1.0, -1.0, 1.0 ) );
    mesh.add_vertex( vector3f( 1.0, -1.0, 1.0 ) );
    mesh.add_vertex( vector3f( -1.0, 1.0, 1.0 ) );
    mesh.add_vertex( vector3f( 1.0, 1.0, 1.0 ) );

    mesh.add_face( vector3( 0, 2, 3 ) );
    mesh.add_face( vector3( 3, 1, 0 ) );
    mesh.add_face( vector3( 4, 5, 7 ) );
    mesh.add_face( vector3( 7, 6, 4 ) );
    mesh.add_face( vector3( 0, 1, 5 ) );
    mesh.add_face( vector3( 5, 4, 0 ) );
    mesh.add_face( vector3( 1, 3, 7 ) );
    mesh.add_face( vector3( 7, 5, 1 ) );
    mesh.add_face( vector3( 3, 2, 6 ) );
    mesh.add_face( vector3( 6, 7, 3 ) );
    mesh.add_face( vector3( 2, 0, 4 ) );
    mesh.add_face( vector3( 4, 6, 2 ) );

    xmesh_sequence_saver xss;
    trimesh3 tempMesh;

    fs::path tempPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%" );
    fs::create_directory( tempPath );
    frantic::tstring filename = to_tstring( tempPath / "test.xmesh" );

    scoped_file_cleanup cleanup;
    cleanup.add( tempPath );

    filename_sequence fseq( filename );

    // make sure i'm using 4 digits for the sequence number
    fseq.get_filename_pattern().set_num_digits( 4 );

    // build the xmesh file names for existence checks
    int fileCount = 0;
    frantic::tstring vertsFile, facesFile, channelVertsFile, channelFacesFile, channelName;
    channelName = _T("someChannel");
    vertsFile = fseq.get_filename_pattern().get_directory( true ) + fseq.get_filename_pattern().get_prefix() +
                _T("_verts0000.xmdat");
    facesFile = fseq.get_filename_pattern().get_directory( true ) + fseq.get_filename_pattern().get_prefix() +
                _T("_faces0000.xmdat");
    channelVertsFile = fseq.get_filename_pattern().get_directory( true ) + fseq.get_filename_pattern().get_prefix() +
                       _T("_channel_") + channelName + _T("_verts0000.xmdat");
    channelFacesFile = fseq.get_filename_pattern().get_directory( true ) + fseq.get_filename_pattern().get_prefix() +
                       _T("_channel_") + channelName + _T("_faces0000.xmdat");

    // write a mesh
    xss.clear();
    xss.write_xmesh( mesh, fseq[fileCount] );

    // load it back in a check for equality
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    fileCount++;

    // write another mesh that hasnt changed
    xss.write_xmesh( mesh, fseq[fileCount] );

    // load it back in and check for equality
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );

    // check that no extra files were created
    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_FALSE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_FALSE( file_exists( facesFile ) );
    fileCount++;

    // alter a vertex
    mesh.vertices_ref()[0] = vector3f( -0.5, -0.5, 0.0 );
    xss.write_xmesh( mesh, fseq[fileCount] );

    // check for equality and new file
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_TRUE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_FALSE( file_exists( facesFile ) );

    fileCount++;

    // add a vertex
    mesh.add_vertex( vector3f( 0.0, 0.0, 2.0 ) );
    xss.write_xmesh( mesh, fseq[fileCount] );
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_TRUE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_FALSE( file_exists( facesFile ) );
    fileCount++;

    // alter a face
    vector3 tempVec3 = mesh.faces_ref()[0];
    mesh.faces_ref()[0] = mesh.faces_ref()[1];
    mesh.faces_ref()[1] = mesh.faces_ref()[0];
    xss.write_xmesh( mesh, fseq[fileCount] );

    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_FALSE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_TRUE( file_exists( facesFile ) );
    fileCount++;

    // add a face
    mesh.add_face( vector3( 7, 6, 8 ) );
    xss.write_xmesh( mesh, fseq[fileCount] );
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_FALSE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_TRUE( file_exists( facesFile ) );
    fileCount++;

    // add both
    mesh.add_vertex( vector3f( 0.0, 0.0, -2.0 ) );
    mesh.add_face( vector3( 0, 2, 9 ) );
    xss.write_xmesh( mesh, fseq[fileCount] );
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_TRUE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_TRUE( file_exists( facesFile ) );
    fileCount++;

    // add a channel without custom facedata
    vector<vector3f> someChannel;
    for( unsigned i = 0; i < mesh.vertices_ref().size(); i++ )
        someChannel.push_back( vector3f( 1.0, 1.0, 1.0 ) );
    mesh.add_vertex_channel_raw( channelName, 3, frantic::channels::data_type_uint32, someChannel.size() );
    memcpy( mesh.get_vertex_channel_general_accessor( channelName ).data( 0 ), &( someChannel[0] ),
            someChannel.size() * sizeof( vector3f ) );
    mesh.set_vertex_channel_custom_faces( channelName, false );

    xss.write_xmesh( mesh, fseq[fileCount] );
    load_xmesh_file( fseq[fileCount], tempMesh );

    vector<frantic::tstring> meshChannels, tempMeshChannels;
    mesh.get_vertex_channel_names( meshChannels );
    mesh.get_vertex_channel_names( tempMeshChannels );

    EXPECT_EQ( 1, meshChannels.size() );
    EXPECT_EQ( meshChannels, tempMeshChannels );
    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ( mesh.get_vertex_channel_general_accessor( channelName ).size(),
               tempMesh.get_vertex_channel_general_accessor( channelName ).size() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( mesh.get_vertex_channel_general_accessor( channelName ).data( 0 ),
                          tempMesh.get_vertex_channel_general_accessor( channelName ).data( 0 ),
                          mesh.get_vertex_channel_general_accessor( channelName ).size() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_FALSE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_FALSE( file_exists( facesFile ) );
    channelVertsFile = replace_sequence_number( channelVertsFile, fileCount );
    EXPECT_TRUE( file_exists( channelVertsFile ) );
    channelFacesFile = replace_sequence_number( channelFacesFile, fileCount );
    EXPECT_FALSE( file_exists( channelFacesFile ) );
    fileCount++;

    // add custom face data
    mesh.set_vertex_channel_custom_faces( channelName, true );
    xss.write_xmesh( mesh, fseq[fileCount] );
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ( mesh.get_vertex_channel_general_accessor( channelName ).size(),
               tempMesh.get_vertex_channel_general_accessor( channelName ).size() );
    EXPECT_EQ( mesh.get_vertex_channel_general_accessor( channelName ).face_count(),
               tempMesh.get_vertex_channel_general_accessor( channelName ).face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( mesh.get_vertex_channel_general_accessor( channelName ).data( 0 ),
                          tempMesh.get_vertex_channel_general_accessor( channelName ).data( 0 ),
                          mesh.get_vertex_channel_general_accessor( channelName ).size() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.get_vertex_channel_general_accessor( channelName ).face( 0 ),
                          &tempMesh.get_vertex_channel_general_accessor( channelName ).face( 0 ),
                          mesh.get_vertex_channel_general_accessor( channelName ).face_count() * sizeof( vector3 ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_FALSE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_FALSE( file_exists( facesFile ) );
    channelVertsFile = replace_sequence_number( channelVertsFile, fileCount );
    EXPECT_FALSE( file_exists( channelVertsFile ) );
    channelFacesFile = replace_sequence_number( channelFacesFile, fileCount );
    EXPECT_TRUE( file_exists( channelFacesFile ) );
    fileCount++;

    // remove custom face data
    mesh.set_vertex_channel_custom_faces( channelName, false );
    xss.write_xmesh( mesh, fseq[fileCount] );
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ( mesh.get_vertex_channel_general_accessor( channelName ).size(),
               tempMesh.get_vertex_channel_general_accessor( channelName ).size() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( mesh.get_vertex_channel_general_accessor( channelName ).data( 0 ),
                          tempMesh.get_vertex_channel_general_accessor( channelName ).data( 0 ),
                          mesh.get_vertex_channel_general_accessor( channelName ).size() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_FALSE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_FALSE( file_exists( facesFile ) );
    channelVertsFile = replace_sequence_number( channelVertsFile, fileCount );
    EXPECT_FALSE( file_exists( channelVertsFile ) );
    channelFacesFile = replace_sequence_number( channelFacesFile, fileCount );
    EXPECT_FALSE( file_exists( channelFacesFile ) );
    fileCount++;

    // remove the channel
    mesh.erase_vertex_channel( channelName );
    xss.write_xmesh( mesh, fseq[fileCount] );
    load_xmesh_file( fseq[fileCount], tempMesh );

    EXPECT_EQ( mesh.vertex_count(), tempMesh.vertex_count() );
    EXPECT_EQ( mesh.face_count(), tempMesh.face_count() );
    EXPECT_EQ(
        0, memcmp( &mesh.vertices_ref()[0], &tempMesh.vertices_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );
    EXPECT_EQ( 0, memcmp( &mesh.faces_ref()[0], &tempMesh.faces_ref()[0], mesh.vertex_count() * sizeof( vector3f ) ) );

    vertsFile = replace_sequence_number( vertsFile, fileCount );
    EXPECT_FALSE( file_exists( vertsFile ) );
    facesFile = replace_sequence_number( facesFile, fileCount );
    EXPECT_FALSE( file_exists( facesFile ) );
    channelVertsFile = replace_sequence_number( channelVertsFile, fileCount );
    EXPECT_FALSE( file_exists( channelVertsFile ) );
    channelFacesFile = replace_sequence_number( channelFacesFile, fileCount );
    EXPECT_FALSE( file_exists( channelFacesFile ) );
    fileCount++;
}

TEST( Trimesh3, RelaxSliverFaces ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;
    make_tetrahedron_beveled_edge_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( relax_sliver_vertices( mesh, 1, 0.2f ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_FALSE( relax_sliver_vertices( mesh, 1, 0.2f ) );

    make_triangle_beveled_vertex_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( relax_sliver_vertices( mesh, 1, 0.2f ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_FALSE( relax_sliver_vertices( mesh, 1, 0.2f ) );

    make_square_beveled_vertex_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( relax_sliver_vertices( mesh, 1, 0.2f ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_FALSE( relax_sliver_vertices( mesh, 1, 0.2f ) );

    make_chained_combined_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( relax_sliver_vertices( mesh, 1, 0.2f ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_FALSE( relax_sliver_vertices( mesh, 1, 0.2f ) );
}

TEST( Trimesh3, RelaxTVertices ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    boost::iostreams::stream<boost::iostreams::null_sink> nullStream( ( boost::iostreams::null_sink() ) );

    trimesh3 mesh;

    make_t_vertex_mesh( mesh );

    EXPECT_FALSE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_TRUE( relax_sliver_vertices( mesh, 1, 0.2f ) );
    EXPECT_TRUE( mesh.check_zero_area_faces( nullStream ) );
    EXPECT_FALSE( relax_sliver_vertices( mesh, 1, 0.2f ) );
}
