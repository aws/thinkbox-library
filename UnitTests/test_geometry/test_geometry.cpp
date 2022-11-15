// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/geometry/raytraced_geometry_collection.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>
#include <frantic/geometry/trimesh3_scan_conversion.hpp>

#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>
#include <frantic/volumetrics/levelset/geometry_to_levelset.hpp>
#include <frantic/volumetrics/levelset/rle_level_set_file_io.hpp>

using namespace std;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::diagnostics;
using namespace frantic::graphics2d;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::volumetrics::implicitsurface;

TEST( GeometryTest, RaytracedGeometryCollection ) {

    raytraced_geometry_collection rgc;
    // Fill in the pcd with a static box and a moving box.
    trimesh3 box;
    //		box.set_to_box( boundbox3f(vector3f(0),vector3f(1)) );
    //		rgc.add_static_object( transform4f::from_translation(3,3,3), box
    //);

    box.set_to_box( boundbox3f( vector3f( 3 ), vector3f( 4 ) ) );
    rgc.add_static_object( transform4f::identity(), box );

    //		motion_blurred_transform mbt(
    // transform4f::from_translation(0,0,0),
    // transform4f::from_translation(1,0,0) );
    //		pcd.add_rigid_object( mbt, box );

    rgc.prepare_kdtrees();
    // pcd.print_statistics(std::cout);

    raytrace_intersection isect;
    ASSERT_TRUE(
        !rgc.intersect_ray( ray3f::from_line_segment( vector3f( 10, 4, 8 ), vector3f( 11, 12, 4 ) ), 0, 1, isect ) );

    ray3f ray = ray3f::from_line_segment( vector3f( 3.5f, 3.5f, 2 ), vector3f( 3.5f, 3.5f, 4 ) );
    double entry, exit;
    ASSERT_TRUE( ray.intersect_with_box( boundbox3f( vector3f( 3 ), vector3f( 4 ) ), entry, exit, 0, 0 ) );
    ASSERT_TRUE( rgc.intersect_ray( ray, 0, 1, isect ) );
    ASSERT_LT( fabs( isect.distance - 0.5 ), 0.000001f );
    ASSERT_LT( vector3f::distance( isect.position, vector3f( 3.5f, 3.5f, 3 ) ), 0.000001f );
}

TEST( GeometryTest, Trimesh3ScanConversion ) {

    trimesh3 mesh;

    /*
    mesh.add_vertex( vector3f(13,4,6) );
    mesh.add_vertex( vector3f(2,18,3) );
    mesh.add_vertex( vector3f(21.6f,14.5f,2) );
    mesh.add_vertex( vector3f(26,10,8) );
    mesh.add_face( vector3(0,1,2) );
    mesh.add_face( vector3(2,3,0) );
    */

    boundbox3f box( 3, 5, 4, 13, 2, 25 );
    mesh.set_to_box( box );

    size2 dim( 30, 20 );

    //			transform4f xform = transform4f::from_translation(15,10,4) *
    // transform4f::from_euler_xyz_rotation_degrees(vector3f(10,20,15)) * transform4f::from_translation(-15,-10,-4);

    //			vector< vector<scan_conversion_intersection> > samples(dim.get_area());
    //			trimesh3_scan_convert( xform, mesh, dim, samples );

    /*
    cout << "\n";
    for( int y = 0; y < dim.ysize; ++y ) {
      for( int x = 0; x < dim.xsize; ++x ) {
        vector<scan_conversion_intersection>& sci = samples[dim.get_index(x,y)];
        //cout << samples[dim.get_index(x,y)].size();
        if( sci.empty() )
          cout << "XX";
        else
          //cout << sci[0].barycentricCoord.get_largest_axis() << " ";
          //cout << (int)sci.back().z << "";
          //cout << sci.size() << "";
          cout << sci.front().normalFacingZPositive << sci.back().normalFacingZPositive;
      }
      cout << "\n";
    }
    cout << endl;
    */

    rle_level_set rls( voxel_coord_system( vector3f(), 1 ) );
    // rle_level_set rls( voxel_coord_system(vector3f(76,12,-32),2.0f) );
    // rle_level_set rls( voxel_coord_system(vector3f(),0.1143f) );
    // rle_level_set rls( voxel_coord_system(vector3f(-4.3688f, -17.4752f, 0.1016f),0.2032f) );

    //			load_mesh_file( "C:\\temp\\mtol.obj", mesh );
    //			load_mesh_file( "C:\\temp\\teapot_small.obj", mesh );
    //			load_mesh_file( "C:\\temp\\testparticles0005.xmesh", mesh );
    //			load_mesh_file( "C:\\temp\\vertplane.obj", mesh );
    //			load_mesh_file( "C:\\temp\\sphere.obj", mesh );
    //			load_mesh_file( "C:\\temp\\plane.obj", mesh );
    //			load_mesh_file( "C:\\temp\\gtoltest.xmesh", mesh );
    //			load_mesh_file( "C:\\temp\\gtols\\gtoltest0044.xmesh", mesh );

    convert_geometry_to_levelset( mesh, -3, 3, rls );
    // convert_geometry_to_levelset( mesh, -3, 2, boundbox3(-72,168,-184,-10,-11,57), rls );

    //			write_rls_rle_level_set_file( "c:\\temp\\testing.rls", rls );

    //			rle_level_set rls2;
    //			read_rle_level_set_file( "c:\\temp\\testing.rls", rls2 );
    //			write_rls_rle_level_set_file( "c:\\temp\\testingRewrite.rls", rls2 );

    //			convert_levelset_to_trimesh3( rls, mesh );
    //			write_mesh_file( "c:\\temp\\test0000.xmesh", mesh );

    /*
    const rle_index_spec& ris = rls.get_rle_index_spec();
    ris_adjacency adj;
    adj.compute(ris);
    for( rle_index_spec_defined_iterator i = ris.begin_x(), ie = ris.end_x(); i != ie; ++i ) {
      //cout << "Testing data index " << i.get_data_index() << endl;
      TS_ASSERT_EQUALS( adj[i.get_data_index()].x_neg, ris.XYZtoDataIndex( i.get_coord() - vector3(1,0,0) ) );
      TS_ASSERT_EQUALS( adj[i.get_data_index()].x_pos, ris.XYZtoDataIndex( i.get_coord() + vector3(1,0,0) ) );
      TS_ASSERT_EQUALS( adj[i.get_data_index()].y_neg, ris.XYZtoDataIndex( i.get_coord() - vector3(0,1,0) ) );
      TS_ASSERT_EQUALS( adj[i.get_data_index()].y_pos, ris.XYZtoDataIndex( i.get_coord() + vector3(0,1,0) ) );
      TS_ASSERT_EQUALS( adj[i.get_data_index()].z_neg, ris.XYZtoDataIndex( i.get_coord() - vector3(0,0,1) ) );
      TS_ASSERT_EQUALS( adj[i.get_data_index()].z_pos, ris.XYZtoDataIndex( i.get_coord() + vector3(0,0,1) ) );
    }
    //*/

    /*
    profiling_section psKD("KD Tree Building"), psDF("Build Distance Function");
    psKD.enter();
    trimesh3_kdtree kd(mesh);
    psKD.exit();
    psDF.enter();
    rls.build_using_distance_function( kd, 3, 3 );
    psDF.exit();
    cout << psKD << endl;
    cout << psDF << endl;
    cout << endl;
    //*/

    cout << "Resultant voxel count: " << rls.get_rle_index_spec().data_size() << endl;

    //			frantic::volumetrics::levelset::detail::convert_intersections_to_level_set( samples, mesh, 4, 4, false,
    //-1,
    // boundbox3(0,29,0,dim.xsize-1,0,dim.ysize-1), rls );
    //			TS_ASSERT( rls.get_rle_index_spec().check_consistency(cout) );
    // legacyflood::write_legacy_rle_level_set( "c:\\temp\\singleconversion.cssf", rls );
    //			frantic::volumetrics::levelset::detail::finalize_geometry_to_levelset_conversion( 4, 4, -1, rls
    //); 			TS_ASSERT( rls.get_rle_index_spec().check_consistency(cout) );
    // legacyflood::write_legacy_rle_level_set( "c:\\temp\\finalized.cssf", rls );
}

TEST( GeometryTest, TricubicMTOL ) {
    rle_level_set rls;

    // legacyflood::read_legacy_rle_level_set( "c:\\temp\\occ_levelset.cssf", rls );

    /*
    legacyflood::read_legacy_rle_level_set(
    "\\\\sim-02\\array\\0082_Voyage\\Seq_RS\\VO_RS_270\\VO_RS_270_P2_plesio_Test_ML_v008\\VO_RS_270_P2_plesio_Test_ML_v008_06_07_2007_2_33_03_PM\\main\\fluidsparseLevelSet_194400.cssf",
    rls );

    trimesh3 mesh;
    convert_levelset_to_trimesh3( rls, rls.get_voxel_coord_system(), 3, mesh );
    write_mesh_file( "c:\\temp\\test0000.xmesh", mesh );
    */

    /*
    const rle_index_spec& ris = rls.get_rle_index_spec();

    boundbox3 randomTestBounds = ris.outer_bounds();
    randomTestBounds.expand(4);
    for( int randomTest = 0; randomTest < 1000000; ++randomTest ) {
      // Create a random box within the randomTestBounds
      boundbox3 testBox;
      testBox += randomTestBounds.random_vector();
      testBox += randomTestBounds.random_vector();

      vector<int> dataIndices(testBox.get_volume());

      // Get the indexes
      ris.fill_data_index_box( testBox, &dataIndices[0] );

      for( int z = testBox.minimum().z; z < testBox.maximum().z; ++z ) {
        for( int y = testBox.minimum().y; y < testBox.maximum().y; ++y ) {
          for( int x = testBox.minimum().x; x < testBox.maximum().x; ++x ) {
            // Compare the values retrieved with the block function versus using the XYZtoDataIndex function
            int rawDataIndex = ris.XYZtoDataIndex( vector3(x,y,z) );
            int blockDataIndex = dataIndices[testBox.size().get_index(vector3(x,y,z) - testBox.minimum())];
            if( rawDataIndex != blockDataIndex ) {
              cerr << "outer bounds: " << ris.outer_bounds() << "\n";
              cerr << "data index box: " << testBox << "\n";
              cerr << "outer position: (" << x << "," << y << "," << z << ")\n";
              cerr << "raw data index: " << rawDataIndex << "\n";
              cerr << "general block data index: " << blockDataIndex << "\n";
            }
            TS_ASSERT_EQUALS( rawDataIndex, blockDataIndex );
          }
        }
      }
    }
    */

    //			voxel_coord_system vcs( vector3f(), rls.get_voxel_coord_system().voxel_length()/1.5f );
    /*
    voxel_coord_system vcs = rls.get_voxel_coord_system();

    trimesh3 mesh;

    boundrect2 extents( rls.get_rle_index_spec().outer_bounds().xminimum() - 2,
      rls.get_rle_index_spec().outer_bounds().xmaximum() + 2,
      rls.get_rle_index_spec().outer_bounds().yminimum() - 2,
      rls.get_rle_index_spec().outer_bounds().ymaximum() + 2 );
    int middleZ = (int)rls.get_rle_index_spec().outer_bounds().center().z+1;
    //*/

    /*
    vector<float> values(extents.get_area());
    {
      rls.fill_plane( extents, middleZ, values );
      ofstream fout("c:\\temp\\raw.txt");
      for( int y = extents.minimum().y; y <= extents.maximum().y; ++y ) {
        for( int x = extents.minimum().x; x <= extents.maximum().x; ++x ) {
          fout << values[(x-extents.minimum().x) + (y-extents.minimum().y)*extents.xsize()] << " ";
        }
        fout << endl;
      }
    }
    //*/

    /*
    {
      ofstream fout("c:\\temp\\ind.txt");
      //for( int z = rls.get_rle_index_spec().outer_bounds().yminimum()-2; z <=
    rls.get_rle_index_spec().outer_bounds().zmaximum()+2; ++z ) {
      int z = middleZ;
      {
        fout << "Z: " << z << endl;
        for( int y = extents.minimum().y; y <= extents.maximum().y; ++y ) {
          levelset::rle_index_spec_block_iterator_x iter( rls.get_rle_index_spec(), extents.minimum().x, y, y, z, z );
          for( int x = extents.minimum().x; x <= extents.maximum().x; ++x ) {
            int actualIndex = rls.get_rle_index_spec().XYZtoDataIndex(vector3(x,y,z));
            fout << actualIndex << "/";
            int index = 0;
            iter.get_indexes(&index);
            fout << index << " ";
            iter.increment_x();
            if( actualIndex != index ) {
              fout << endl;
              throw std::runtime_error( "Mismatch!" );
            }
          }
          fout << endl;
        }
      }
    }
    //*/

    /*
    {
      rls.fill_plane( vcs, extents, middleZ, values );
      ofstream fout("c:\\temp\\tric.txt");
      for( int y = extents.minimum().y; y <= extents.maximum().y; ++y ) {
        for( int x = extents.minimum().x; x <= extents.maximum().x; ++x ) {
          fout << values[(x-extents.minimum().x) + (y-extents.minimum().y)*extents.xsize()] << " ";
        }
        fout << endl;
      }
    }
    //*/

    //			convert_levelset_to_trimesh3( rls, vcs, 2, mesh );
    //			write_mesh_file( "c:\\temp\\test0000.xmesh", mesh );
}
