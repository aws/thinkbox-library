// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>

namespace {

unsigned char get_cube_case( const float cubeCorners[8] ) {
    unsigned char result = 0;

    for( int i = 0; i < 8; ++i ) {
        result >>= 1;

        const float density = cubeCorners[i];

        result |= ( density < 0 ? 0x80 : 0 );
    }

    return result;
}

} // anonymous namespace

TEST( MarchingCubes, Configuration13WithZeroAndNegativeOne ) {
    const float cubeCorners[8] = { -1, 0, 0, -1, 0, -1, -1, 0 };
    const unsigned char cubeCase = get_cube_case( cubeCorners );

    frantic::volumetrics::marching_cubes_table mct;

    std::vector<frantic::graphics::vector3> faces;
    ASSERT_NO_THROW( mct.get_cubecase_faces( cubeCase, cubeCorners, faces ) );
}

TEST( MarchingCubes, Configuration13EpsilonRegressionTest ) {
    // This case produced an error in marching_cubes_table.get_cubecase_faces(),
    // because marching_cubes_table::test_face() erroneously returned true
    // for two extra faces due to a large epsilon.
    const float cubeCorners[8] = { -0.001757f, 0.000576f,  0.000104f,  -0.000239f,
                                   0.000921f,  -0.000404f, -0.000136f, 0.000078f };
    const unsigned char cubeCase = get_cube_case( cubeCorners );

    frantic::volumetrics::marching_cubes_table mct;

    std::vector<frantic::graphics::vector3> faces;
    ASSERT_NO_THROW( mct.get_cubecase_faces( cubeCase, cubeCorners, faces ) );
}

TEST( MarchingCubes, Topology ) {
    using namespace std;
    using namespace frantic::graphics;
    using namespace frantic::geometry;
    using namespace frantic::volumetrics;
    using namespace frantic::volumetrics::levelset;
    using namespace frantic::volumetrics::implicitsurface;

    /*
    cubeCaseCounter["1"] = 0;
    cubeCaseCounter["2"] = 0;
    cubeCaseCounter["3.1"] = 0;
    cubeCaseCounter["3.2"] = 0;
    cubeCaseCounter["4.1"] = 0;
    cubeCaseCounter["4.2"] = 0;
    cubeCaseCounter["5"] = 0;
    cubeCaseCounter["6.1.1"] = 0;
    cubeCaseCounter["6.1.2"] = 0;
    cubeCaseCounter["6.2"] = 0;
    cubeCaseCounter["7.1"] = 0;
    cubeCaseCounter["7.2"] = 0;
    cubeCaseCounter["7.3"] = 0;
    cubeCaseCounter["7.4.1"] = 0;
    cubeCaseCounter["7.4.2"] = 0;
    cubeCaseCounter["8"] = 0;
    cubeCaseCounter["9"] = 0;
    cubeCaseCounter["10.1.1"] = 0;
    cubeCaseCounter["10.1.2"] = 0;
    cubeCaseCounter["10.2"] = 0;
    cubeCaseCounter["11"] = 0;
    cubeCaseCounter["12.1.1"] = 0;
    cubeCaseCounter["12.1.2"] = 0;
    cubeCaseCounter["12.2"] = 0;
    cubeCaseCounter["13.1"] = 0;
    cubeCaseCounter["13.2"] = 0;
    cubeCaseCounter["13.3"] = 0;
    cubeCaseCounter["13.4"] = 0;
    cubeCaseCounter["13.5.1"] = 0;
    cubeCaseCounter["13.5.2"] = 0;
    */

    // Generate a 20x20x20 chunk of data that is a closed surface.
    size3 dim( 10, 10, 10 );

    int numTests = 10;
    unsigned seed = (unsigned)time( NULL );
    std::cout << "seed: " << seed << std::endl;
    for( int n = 0; n < numTests; ++n ) {
        ++seed;
        srand( seed );
        // std::cout << "test: " << n << "\t\tseed: " << seed << std::endl;

        std::vector<frantic::graphics::vector3> voxelCoords;
        vector<float> distanceData;

        for( int x = 0; x < dim.xsize(); ++x ) {
            for( int y = 0; y < dim.ysize(); ++y ) {
                for( int z = 0; z < dim.zsize(); ++z ) {

                    // fill randomly with outside values only
                    // voxelCoords.push_back(vector3(x,y,z));
                    // distanceData.push_back( 10.f *((float)rand()/RAND_MAX));

                    // fill randomly with inside/outside values
                    if( x == 0 || x == dim.xsize() - 1 || y == 0 || y == dim.ysize() - 1 || z == 0 ||
                        z == dim.zsize() - 1 ) {
                        voxelCoords.push_back( vector3( x, y, z ) );
                        distanceData.push_back( 10.f * ( (float)rand() / RAND_MAX ) );
                    } else {
                        voxelCoords.push_back( vector3( x, y, z ) );
                        distanceData.push_back( 20.f * ( (float)rand() / RAND_MAX ) - 10.f );
                    }
                }
            }
        }

        rle_index_spec ris;
        ris.build_from_voxel_array( voxelCoords );
        // ris.dump(std::cout);

        voxel_coord_system vcs( vector3f( 0.f, 0.f, 0.f ), 1.f );
        rle_level_set rls( vcs, ris, distanceData, 1.f, 1.f );
        rls.add_channel<vector3f>( _T("Velocity") );

        frantic::geometry::trimesh3 mesh;
        frantic::volumetrics::implicitsurface::convert_levelset_to_trimesh3( rls, mesh, 0 );
        // trimesh3_face_channel_accessor<int> caseAcc = mesh.get_face_channel_accessor<int>("CubeCase");
        // trimesh3_face_channel_accessor<vector3> configAcc = mesh.get_face_channel_accessor<vector3>("CubeConfig");

        // test to make sure that each edge in the resulting mesh is shared by exactly two triangles.
        // it's just a test, we're not really concerned with speed/efficiency.
        std::vector<vector3>& faces = mesh.faces_ref();
        for( int i = 0; i < (int)faces.size(); ++i ) {
            for( int j = 0; j < 3; ++j ) {
                std::vector<vector3> sharedFaces;
                int v0 = faces[i][j], v1 = faces[i][( j + 1 ) % 3];
                for( int a = 0; a < (int)faces.size(); ++a ) {

                    for( int b = 0; b < 3; ++b ) {
                        int u0 = faces[a][b], u1 = faces[a][( b + 1 ) % 3];
                        if( v0 == u0 && v1 == u1 || v0 == u1 && v1 == u0 ) {
                            sharedFaces.push_back( faces[a] );
                        }
                    }
                }

                if( sharedFaces.size() != 2 ) {

                    // ignore any duplicated faces, which can result from certain cube case combinations
                    int sameFaceCount = 0;
                    for( int a = 0; a < (int)sharedFaces.size(); ++a ) {
                        std::set<int> face_a;
                        face_a.insert( sharedFaces[a].x );
                        face_a.insert( sharedFaces[a].y );
                        face_a.insert( sharedFaces[a].z );

                        for( int b = 0; b < (int)sharedFaces.size(); ++b ) {
                            std::set<int> face_b;
                            if( a != b ) {
                                face_b.insert( sharedFaces[b].x );
                                face_b.insert( sharedFaces[b].y );
                                face_b.insert( sharedFaces[b].z );
                                if( face_a == face_b ) {
                                    // std::cout << a << ": " << sharedFaces[a] << "   " << b << ": " << sharedFaces[b]
                                    // << std::endl;
                                    ++sameFaceCount;
                                }
                            }
                        }
                    }
                    size_t uniqueFaceCount = sharedFaces.size() - sameFaceCount;
                    if( uniqueFaceCount == 2 || uniqueFaceCount == 0 )
                        continue;

                    // write_rls_rle_level_set_file("C:\\temp\\marchingcubes\\levelset.rls", rls);
                    // write_mesh_file("C:\\temp\\marchingcubes\\mesh.xmesh", mesh);
                    // int total = 0;
                    // for (std::map<std::string,int>::iterator iter = cubeCaseCounter.begin(); iter !=
                    // cubeCaseCounter.end();
                    // ++iter) { 	std::cout << iter->first << "\t\t" << iter->second << std::endl; 	total +=
                    // iter->second;
                    // }
                    // std::cout << "total\t\t" << total << std::endl;

                    // std::cout << std::string("face " + boost::lexical_cast<std::string>(i) + " edge " +
                    // boost::lexical_cast<std::string>(j) + " used " +
                    // boost::lexical_cast<std::string>(sharedFaces.size()) + " times.") << std::endl; std::cout <<
                    // "shared by faces: " << std::endl; for ( int q = 0; q < sharedFaces.size(); ++q) 	std::cout <<
                    // "\t" << sharedFaces[q] << std::endl; std::cout << std::endl; std::cout << "face generated by cube
                    // case: " << caseAcc[i] << "  config: " << configAcc[i] << std::endl; std::cout << std::endl;
                    ADD_FAILURE() << "testMarchingCubesTopology() - Topology error:  face " +
                                         boost::lexical_cast<std::string>( i ) + " edge " +
                                         boost::lexical_cast<std::string>( j ) + " used " +
                                         boost::lexical_cast<std::string>( sharedFaces.size() ) + " times."
                                  << std::endl;
                }
            }
        }
        // make sure that the data channel has the same amount of vertex info as the mesh
        trimesh3_vertex_channel_accessor<vector3f> velAcc =
            mesh.get_vertex_channel_accessor<vector3f>( _T("Velocity") );
        EXPECT_FALSE( velAcc.has_custom_faces() )
            << "The vertex channel generated during meshing has custom faces." << std::endl;
        EXPECT_EQ( velAcc.size(), mesh.vertex_count() )
            << "The vertex channel generated during meshing does not have the same vertex count as the mesh."
            << std::endl;
    }
    /*
    int total = 0;
    for (std::map<std::string,int>::iterator iter = cubeCaseCounter.begin(); iter != cubeCaseCounter.end(); ++iter) {
      std::cout << iter->first << "\t\t" << iter->second << std::endl;
      total += iter->second;
    }
    std::cout << "total\t\t" << total << std::endl;
    */
}

TEST( MarchingCubes, FixMarchingCubesTopology ) {
    tbb::task_scheduler_init taskScheduleInit;
    using namespace std;
    using namespace frantic::graphics;
    using namespace frantic::geometry;
    using namespace frantic::volumetrics;
    using namespace frantic::volumetrics::levelset;
    using namespace frantic::volumetrics::implicitsurface;

    boost::random::mt19937 gen;
    boost::random::uniform_real_distribution<float> dist( -1, 1 );

    // Generate a 10x10x10 chunk of data that is a closed surface.
    size3 dim( 10, 10, 10 );

    int numTests = 20;

    const unsigned badSeed = 1441218265; // known bad case
    unsigned seed = (unsigned)time( NULL );

    bool foundAnyNonClosedManifold = false;

    for( int n = 0; n < numTests; ++n ) {
        ++seed;

        gen.seed( n == 0 ? badSeed : seed );
        dist.reset();

        std::vector<frantic::graphics::vector3> voxelCoords;
        vector<float> distanceData;

        for( int x = 0; x < dim.xsize(); ++x ) {
            for( int y = 0; y < dim.ysize(); ++y ) {
                for( int z = 0; z < dim.zsize(); ++z ) {
                    // fill randomly with inside/outside values
                    float value = dist( gen );

                    // fill edges with outside (positive) values only
                    if( x == 0 || x == dim.xsize() - 1 || y == 0 || y == dim.ysize() - 1 || z == 0 ||
                        z == dim.zsize() - 1 ) {
                        value = std::abs( value ) + std::numeric_limits<float>::epsilon();
                    }

                    voxelCoords.push_back( vector3( x, y, z ) );
                    distanceData.push_back( value );
                }
            }
        }

        rle_index_spec ris;
        ris.build_from_voxel_array( voxelCoords );

        voxel_coord_system vcs( vector3f( 0.f, 0.f, 0.f ), 1.f );
        rle_level_set rls( vcs, ris, distanceData, 1.f, 1.f );

        frantic::geometry::trimesh3 mesh;
        frantic::volumetrics::implicitsurface::convert_levelset_to_trimesh3( rls, mesh, 0 );

        boost::shared_ptr<frantic::geometry::mesh_interface> meshInterface;
        meshInterface.reset( frantic::geometry::trimesh3_interface::create_instance( &mesh ).release() );

        if( !frantic::geometry::is_closed_manifold( meshInterface ) ) {
            foundAnyNonClosedManifold = true;

            frantic::logging::null_progress_logger progress;

            frantic::volumetrics::implicitsurface::fix_marching_cubes_topology( mesh, progress );
        }

        meshInterface.reset( frantic::geometry::trimesh3_interface::create_instance( &mesh ).release() );
        EXPECT_TRUE( frantic::geometry::is_closed_manifold( meshInterface ) )
            << "The mesh is not a closed manifold after fix.  Seed: " << seed << std::endl;
    }

    EXPECT_TRUE( foundAnyNonClosedManifold )
        << "All input meshes were closed manifold.\n\n"
        << "If you've fixed the mesh generation so it always outputs a closed manifold,\n"
        << "great!  You can remove this test.  Otherwise, please fix this test so that it\n"
        << "always operates on at least one mesh that is not a closed manifold." << std::endl;
}
