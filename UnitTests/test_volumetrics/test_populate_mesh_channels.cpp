// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "tbb/task_scheduler_init.h"

#include <frantic/volumetrics/implicitsurface/level_set_implicit_surface_policies.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

using namespace frantic::graphics;
using namespace frantic::geometry;
using namespace frantic::channels;
using namespace frantic::particles;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::volumetrics::implicitsurface;
using namespace frantic::volumetrics::implicitsurface::detail;

TEST( PopulateMeshChannels, DirectLinearRleLevelSetISpolicy ) {
    voxel_coord_system vcs( vector3f( 0.f, 0.f, 0.f ), 1 );
    std::vector<vector3> coordinateArray;
    std::vector<float> distanceData;

    for( boost::int32_t z = 0; z < 2; z++ ) {
        for( boost::int32_t y = 0; y < 2; y++ ) {
            for( boost::int32_t x = 0; x < 2; x++ ) {
                coordinateArray.push_back( vector3( x, y, z ) );
                distanceData.push_back( 1.f );
            }
        }
    }

    rle_index_spec ris;
    ris.build_from_voxel_array( coordinateArray );

    rle_level_set rls( vcs, ris, distanceData, 1, 1 );

    rls.add_channel<vector3f>( _T("Velocity") );
    rle_channel_accessor<vector3f> velAcc = rls.get_channel_accessor<vector3f>( _T("Velocity") );

    vector3f vec;
    for( size_t z = 0; z < 2; z++ ) {
        for( size_t y = 0; y < 2; y++ ) {
            for( size_t x = 0; x < 2; x++ ) {
                vec = vector3f( float( x ), float( y ), float( z ) );
                memcpy( velAcc.data( 4 * z + 2 * y + x ), &vec, sizeof( vector3f ) );
            }
        }
    }

    trimesh3 mesh;

    // check for points landed on
    mesh.add_vertex( 0.5f, 0.5f, 0.5f );
    mesh.add_vertex( 1.5f, 0.5f, 0.5f );
    mesh.add_vertex( 0.5f, 1.5f, 0.5f );
    mesh.add_vertex( 1.5f, 1.5f, 0.5f );
    mesh.add_vertex( 0.5f, 0.5f, 1.5f );
    mesh.add_vertex( 1.5f, 0.5f, 1.5f );
    mesh.add_vertex( 0.5f, 1.5f, 1.5f );
    mesh.add_vertex( 1.5f, 1.5f, 1.5f );

    // check for point in middle
    mesh.add_vertex( 1.f, 1.f, 1.f );

    // check for obscure points inside
    mesh.add_vertex( 0.7f, 0.2f, 0.3f );

    // check for points outside box, but inside range
    mesh.add_vertex( 2.f, 2.f, 2.f );

    // check for points outside range
    mesh.add_vertex( 3.f, 3.f, 3.f );

    direct_linear_rle_level_set_is_policy isp( rls );

    channel_propagation_policy cpp( true );
    cpp.add_channel( _T("Velocity") );
    isp.populate_mesh_channels( mesh, cpp );

    // COMFIRM VERTEX CHANNELS
    trimesh3_vertex_channel_accessor<vector3f> acc = mesh.get_vertex_channel_accessor<vector3f>( _T("Velocity") );

    EXPECT_TRUE( ( *(vector3f*)acc.data( 0 ) ).is_equal( vector3f( 0.f, 0.f, 0.f ), 0.00001f ) ); // pos = 0.5,0.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 1 ) ).is_equal( vector3f( 1.f, 0.f, 0.f ), 0.00001f ) ); // pos = 1.5,0.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 2 ) ).is_equal( vector3f( 0.f, 1.f, 0.f ), 0.00001f ) ); // pos = 0.5,1.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 3 ) ).is_equal( vector3f( 1.f, 1.f, 0.f ), 0.00001f ) ); // pos = 1.5,1.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 4 ) ).is_equal( vector3f( 0.f, 0.f, 1.f ), 0.00001f ) ); // pos = 0.5,0.5,1.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 5 ) ).is_equal( vector3f( 1.f, 0.f, 1.f ), 0.00001f ) ); // pos = 1.5,0.5,1.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 6 ) ).is_equal( vector3f( 0.f, 1.f, 1.f ), 0.00001f ) ); // pos = 0.5,1.5,1.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 7 ) ).is_equal( vector3f( 1.f, 1.f, 1.f ), 0.00001f ) ); // pos = 1.5,1.5,1.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 8 ) ).is_equal( vector3f( 0.5f, 0.5f, 0.5f ), 0.00001f ) ); // pos = 0.5,0.5,0.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 9 ) ).is_equal( vector3f( 0.112f, 0.f, 0.f ), 0.00001f ) ); // pos = 0.7,0.2,0.3
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 10 ) ).is_equal( vector3f( 0.125f, 0.125f, 0.125f ), 0.00001f ) ); // pos = 2,2,2
    EXPECT_TRUE( ( *(vector3f*)acc.data( 11 ) ).is_equal( vector3f( 0.f, 0.f, 0.f ), 0.00001f ) ); // pos = 3,3,3
}

TEST( PopulateMeshChannels, MeatballISPolicy ) {
    tbb::task_scheduler_init taskScheduleInit;

    trimesh3 mesh;

    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<vector3f>( _T("Velocity") );
    pcm.define_channel<float>( _T("Radius") );
    pcm.end_channel_definition();
    pcm.channel_definition_complete();

    particle_grid_tree pgt( pcm );

    channel_accessor<vector3f> posA = pcm.get_accessor<vector3f>( _T("Position") );
    channel_accessor<vector3f> velA = pcm.get_accessor<vector3f>( _T("Velocity") );
    channel_accessor<float> radA = pcm.get_accessor<float>( _T("Radius") );
    std::vector<char> particle( pcm.structure_size() );
    for( size_t z = 0; z < 2; z++ ) {
        for( size_t y = 0; y < 2; y++ ) {
            for( size_t x = 0; x < 2; x++ ) {
                posA( &particle[0] ) = vector3f( (float)x, (float)y, (float)z );
                velA( &particle[0] ) = vector3f( (float)( x + 1 ), (float)( y + 1 ), (float)( z + 1 ) );
                radA( &particle[0] ) = 1.f;

                pgt.insert( &particle[0] );
            }
        }
    }
    // check for points landed on
    mesh.add_vertex( 0.f, 0.f, 0.f );
    mesh.add_vertex( 1.f, 0.f, 0.f );
    mesh.add_vertex( 0.f, 1.f, 0.f );
    mesh.add_vertex( 1.f, 1.f, 0.f );
    mesh.add_vertex( 0.f, 0.f, 1.f );
    mesh.add_vertex( 1.f, 0.f, 1.f );
    mesh.add_vertex( 0.f, 1.f, 1.f );
    mesh.add_vertex( 1.f, 1.f, 1.f );

    // check for point in middle
    mesh.add_vertex( 0.5f, 0.5f, 0.5f );

    // check for obscure points inside
    mesh.add_vertex( 0.7f, 0.2f, 0.3f );

    // check for points outside box, but inside range
    mesh.add_vertex( -1.1f, -1.1f, -1.1f );

    // check for points outside range
    mesh.add_vertex( 3.f, 3.f, 3.f );

    voxel_coord_system vcs( vector3f(), 1.f );
    particle_metaball_is_policy isp( pgt, 1.f, 2.f, 1.f, vcs, 1 );

    channel_propagation_policy cpp( true );
    cpp.add_channel( _T("Velocity") );
    isp.populate_mesh_channels( mesh, cpp );

    // Slowest way of calculating metaball channels.
    // Used to check results.
    std::vector<vector3f> result( mesh.vertex_count() );
    std::vector<float> d( 8 );
    float distance;
    float sum;
    float particleEffectRadius = 2.f;
    vector3f pos;
    for( size_t i = 0; i < mesh.vertex_count(); i++ ) {
        sum = 0;
        pos = mesh.get_vertex( i );
        for( size_t z = 0; z < 2; z++ ) {
            for( size_t y = 0; y < 2; y++ ) {
                for( size_t x = 0; x < 2; x++ ) {
                    // embedded metaball function, b/c real one is not accessible
                    distance = vector3f::distance( vector3f( (float)x, (float)y, (float)z ), pos );
                    if( distance < ( 0.33333333f ) * particleEffectRadius )
                        d[4 * z + 2 * y + x] = 1 - 3 * frantic::math::square( distance / particleEffectRadius );
                    else if( distance < particleEffectRadius )
                        d[4 * z + 2 * y + x] =
                            ( 3.f / 2.f ) * frantic::math::square( 1 - distance / particleEffectRadius );
                    else
                        d[4 * z + 2 * y + x] = 0;
                }
            }
        }

        for( int w = 0; w < 8; w++ )
            sum += d[w];
        for( int w = 0; w < 8 && sum != 0; w++ )
            d[w] /= sum;
        result[i] = d[0] * vector3f( 1, 1, 1 ) + d[1] * vector3f( 2, 1, 1 ) + d[2] * vector3f( 1, 2, 1 ) +
                    d[3] * vector3f( 2, 2, 1 ) + d[4] * vector3f( 1, 1, 2 ) + d[5] * vector3f( 2, 1, 2 ) +
                    d[6] * vector3f( 1, 2, 2 ) + d[7] * vector3f( 2, 2, 2 );
    }

    // COMFIRM VERTEX CHANNELS
    trimesh3_vertex_channel_accessor<frantic::graphics::vector3f> acc =
        mesh.get_vertex_channel_accessor<frantic::graphics::vector3f>( _T("Velocity") );

    EXPECT_TRUE( ( *(vector3f*)acc.data( 0 ) ).is_equal( result[0], 0.00001f ) );   // pos = 0,0,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 1 ) ).is_equal( result[1], 0.00001f ) );   // pos = 1,0,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 2 ) ).is_equal( result[2], 0.00001f ) );   // pos = 0,1,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 3 ) ).is_equal( result[3], 0.00001f ) );   // pos = 1,1,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 4 ) ).is_equal( result[4], 0.00001f ) );   // pos = 0,0,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 5 ) ).is_equal( result[5], 0.00001f ) );   // pos = 1,0,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 6 ) ).is_equal( result[6], 0.00001f ) );   // pos = 0,1,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 7 ) ).is_equal( result[7], 0.00001f ) );   // pos = 1,1,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 8 ) ).is_equal( result[8], 0.00001f ) );   // pos = 0.5,0.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 9 ) ).is_equal( result[9], 0.00001f ) );   // pos = 0.7,0.2,0.3
    EXPECT_TRUE( ( *(vector3f*)acc.data( 10 ) ).is_equal( result[10], 0.00001f ) ); // pos = -1.1,-1.1,-1.1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 11 ) ).is_equal( result[11], 0.00001f ) ); // pos = 3,3,3
}

TEST( PopulateMeshChannels, ParticleZhuBridsonISPolicy ) {
    tbb::task_scheduler_init taskScheduleInit;

    trimesh3 mesh;

    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<vector3f>( _T("Velocity") );
    pcm.define_channel<float>( _T("Radius") );
    pcm.end_channel_definition();
    pcm.channel_definition_complete();

    particle_grid_tree pgt( pcm );

    channel_accessor<vector3f> posA = pcm.get_accessor<vector3f>( _T("Position") );
    channel_accessor<vector3f> velA = pcm.get_accessor<vector3f>( _T("Velocity") );
    channel_accessor<float> radA = pcm.get_accessor<float>( _T("Radius") );
    std::vector<char> particle( pcm.structure_size() );
    for( size_t z = 0; z < 2; z++ ) {
        for( size_t y = 0; y < 2; y++ ) {
            for( size_t x = 0; x < 2; x++ ) {
                posA( &particle[0] ) = vector3f( (float)x, (float)y, (float)z );
                velA( &particle[0] ) = vector3f( (float)( x + 1 ), (float)( y + 1 ), (float)( z + 1 ) );
                radA( &particle[0] ) = 1.f;

                pgt.insert( &particle[0] );
            }
        }
    }
    // check for points landed on
    mesh.add_vertex( 0.f, 0.f, 0.f );
    mesh.add_vertex( 1.f, 0.f, 0.f );
    mesh.add_vertex( 0.f, 1.f, 0.f );
    mesh.add_vertex( 1.f, 1.f, 0.f );
    mesh.add_vertex( 0.f, 0.f, 1.f );
    mesh.add_vertex( 1.f, 0.f, 1.f );
    mesh.add_vertex( 0.f, 1.f, 1.f );
    mesh.add_vertex( 1.f, 1.f, 1.f );

    // check for point in middle
    mesh.add_vertex( 0.5f, 0.5f, 0.5f );

    // check for obscure points inside
    mesh.add_vertex( 0.7f, 0.2f, 0.3f );

    // check for points outside box, but inside range
    mesh.add_vertex( -1.1f, -1.1f, -1.1f );

    // check for points outside range
    mesh.add_vertex( 3.f, 3.f, 3.f );

    voxel_coord_system vcs( vector3f(), 1.f );
    particle_zhu_bridson_is_policy isp( pgt, 1.f, 2.f, 1.f, 1.f, vcs, 1 );

    channel_propagation_policy cpp( true );
    cpp.add_channel( _T("Velocity") );
    isp.populate_mesh_channels( mesh, cpp );

    std::vector<vector3f> result( mesh.vertex_count() );
    std::vector<float> d( 8 );
    float kernelCompactSupportSquared = 4.f;
    float distanceSquared;
    float sum;
    vector3f pos;
    for( size_t i = 0; i < mesh.vertex_count(); i++ ) {
        sum = 0;
        pos = mesh.get_vertex( i );
        for( size_t z = 0; z < 2; z++ ) {
            for( size_t y = 0; y < 2; y++ ) {
                for( size_t x = 0; x < 2; x++ ) {
                    distanceSquared = vector3f::distance_squared( vector3f( float( x ), float( y ), float( z ) ), pos );
                    if( distanceSquared < kernelCompactSupportSquared ) {
                        float k = 1 - distanceSquared / kernelCompactSupportSquared;
                        d[4 * z + 2 * y + x] = k * k * k;
                    } else {
                        d[4 * z + 2 * y + x] = 0;
                    }
                }
            }
        }

        for( int w = 0; w < 8; w++ )
            sum += d[w];
        for( int w = 0; w < 8 && sum != 0; w++ )
            d[w] /= sum;
        result[i] = d[0] * vector3f( 1, 1, 1 ) + d[1] * vector3f( 2, 1, 1 ) + d[2] * vector3f( 1, 2, 1 ) +
                    d[3] * vector3f( 2, 2, 1 ) + d[4] * vector3f( 1, 1, 2 ) + d[5] * vector3f( 2, 1, 2 ) +
                    d[6] * vector3f( 1, 2, 2 ) + d[7] * vector3f( 2, 2, 2 );
    }

    // COMFIRM VERTEX CHANNELS
    trimesh3_vertex_channel_accessor<frantic::graphics::vector3f> acc =
        mesh.get_vertex_channel_accessor<frantic::graphics::vector3f>( _T("Velocity") );

    EXPECT_TRUE( ( *(vector3f*)acc.data( 0 ) ).is_equal( result[0], 0.00001f ) );   // pos = 0,0,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 1 ) ).is_equal( result[1], 0.00001f ) );   // pos = 1,0,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 2 ) ).is_equal( result[2], 0.00001f ) );   // pos = 0,1,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 3 ) ).is_equal( result[3], 0.00001f ) );   // pos = 1,1,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 4 ) ).is_equal( result[4], 0.00001f ) );   // pos = 0,0,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 5 ) ).is_equal( result[5], 0.00001f ) );   // pos = 1,0,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 6 ) ).is_equal( result[6], 0.00001f ) );   // pos = 0,1,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 7 ) ).is_equal( result[7], 0.00001f ) );   // pos = 1,1,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 8 ) ).is_equal( result[8], 0.00001f ) );   // pos = 0.5,0.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 9 ) ).is_equal( result[9], 0.00001f ) );   // pos = 0.7,0.2,0.3
    EXPECT_TRUE( ( *(vector3f*)acc.data( 10 ) ).is_equal( result[10], 0.00001f ) ); // pos = -1.1,-1.1,-1.1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 11 ) ).is_equal( result[11], 0.00001f ) ); // pos = 3,3,3
}

TEST( PopulateMeshChannels, ReconsturctionFilteredRleLevelSetISPolicy ) {
    voxel_coord_system vcs( vector3f( 0.f, 0.f, 0.f ), 1 );
    std::vector<vector3> coordinateArray;
    std::vector<float> distanceData;

    for( boost::int32_t z = 0; z < 2; z++ ) {
        for( boost::int32_t y = 0; y < 2; y++ ) {
            for( boost::int32_t x = 0; x < 2; x++ ) {
                coordinateArray.push_back( vector3( x, y, z ) );
                distanceData.push_back( 1.f );
            }
        }
    }

    rle_index_spec ris;
    ris.build_from_voxel_array( coordinateArray );

    rle_level_set rls( vcs, ris, distanceData, 1, 1 );

    rls.add_channel<vector3f>( _T("Velocity") );
    rle_channel_accessor<vector3f> velAcc = rls.get_channel_accessor<vector3f>( _T("Velocity") );

    vector3f vec;
    for( size_t z = 0; z < 2; z++ )
        for( size_t y = 0; y < 2; y++ )
            for( size_t x = 0; x < 2; x++ ) {
                vec = vector3f( float( x ), float( y ), float( z ) );
                memcpy( velAcc.data( 4 * z + 2 * y + x ), &vec, sizeof( vector3f ) );
            }

    trimesh3 mesh;

    // check for points landed on
    mesh.add_vertex( 0.5f, 0.5f, 0.5f );
    mesh.add_vertex( 1.5f, 0.5f, 0.5f );
    mesh.add_vertex( 0.5f, 1.5f, 0.5f );
    mesh.add_vertex( 1.5f, 1.5f, 0.5f );
    mesh.add_vertex( 0.5f, 0.5f, 1.5f );
    mesh.add_vertex( 1.5f, 0.5f, 1.5f );
    mesh.add_vertex( 0.5f, 1.5f, 1.5f );
    mesh.add_vertex( 1.5f, 1.5f, 1.5f );

    // check for point in middle
    mesh.add_vertex( 1.f, 1.f, 1.f );

    // check for obscure points inside
    mesh.add_vertex( 0.7f, 0.2f, 0.3f );

    // check for points outside box, but inside range
    mesh.add_vertex( 2.f, 2.f, 2.f );

    // check for points that only hit one point
    mesh.add_vertex( 3.f, 3.f, 3.f );

    // check for points outside range
    mesh.add_vertex( 5.f, 5.f, 5.f );

    frantic::math::b_spline_filter reconFilter;
    reconstruction_filtered_rle_level_set_is_policy<frantic::math::b_spline_filter> isp( rls, vcs, 1, reconFilter );

    channel_propagation_policy cpp( true );
    cpp.add_channel( _T("Velocity") );
    isp.populate_mesh_channels( mesh, cpp );

    // COMFIRM VERTEX CHANNELS
    trimesh3_vertex_channel_accessor<vector3f> acc = mesh.get_vertex_channel_accessor<vector3f>( _T("Velocity") );

    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 0 ) ).is_equal( vector3f( 0.2f, 0.2f, 0.2f ), 0.00001f ) ); // pos = 0.5,0.5,0.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 1 ) ).is_equal( vector3f( 0.8f, 0.2f, 0.2f ), 0.00001f ) ); // pos = 1.5,0.5,0.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 2 ) ).is_equal( vector3f( 0.2f, 0.8f, 0.2f ), 0.00001f ) ); // pos = 0.5,1.5,0.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 3 ) ).is_equal( vector3f( 0.8f, 0.8f, 0.2f ), 0.00001f ) ); // pos = 1.5,1.5,0.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 4 ) ).is_equal( vector3f( 0.2f, 0.2f, 0.8f ), 0.00001f ) ); // pos = 0.5,0.5,1.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 5 ) ).is_equal( vector3f( 0.8f, 0.2f, 0.8f ), 0.00001f ) ); // pos = 1.5,0.5,1.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 6 ) ).is_equal( vector3f( 0.2f, 0.8f, 0.8f ), 0.00001f ) ); // pos = 0.5,1.5,1.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 7 ) ).is_equal( vector3f( 0.8f, 0.8f, 0.8f ), 0.00001f ) ); // pos = 1.5,1.5,1.5
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 8 ) ).is_equal( vector3f( 0.5f, 0.5f, 0.5f ), 0.00001f ) ); // pos = 0.5,0.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 9 ) )
                     .is_equal( vector3f( 0.309489f, 0.088311f, 0.119181f ), 0.00001f ) ); // pos = 0.7,0.2,0.3
    EXPECT_TRUE( ( *(vector3f*)acc.data( 10 ) )
                     .is_equal( vector3f( 0.958333f, 0.958333f, 0.958333f ), 0.00001f ) );         // pos = 2,2,2
    EXPECT_TRUE( ( *(vector3f*)acc.data( 11 ) ).is_equal( vector3f( 1.f, 1.f, 1.f ), 0.00001f ) ); // pos = 3,3,3
    EXPECT_TRUE( ( *(vector3f*)acc.data( 12 ) ).is_equal( vector3f( 0.f, 0.f, 0.f ), 0.00001f ) ); // pos = 5,5,5
}

TEST( PopulateMeshChannels, UnionOfSphereISPolicy ) {
    tbb::task_scheduler_init taskScheduleInit;

    trimesh3 mesh;

    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<vector3f>( _T("Velocity") );
    pcm.define_channel<float>( _T("Radius") );
    pcm.end_channel_definition();
    pcm.channel_definition_complete();

    particle_grid_tree pgt( pcm );

    channel_accessor<vector3f> posA = pcm.get_accessor<vector3f>( _T("Position") );
    channel_accessor<vector3f> velA = pcm.get_accessor<vector3f>( _T("Velocity") );
    channel_accessor<float> radA = pcm.get_accessor<float>( _T("Radius") );
    std::vector<char> particle( pcm.structure_size() );
    for( size_t z = 0; z < 2; z++ )
        for( size_t y = 0; y < 2; y++ )
            for( size_t x = 0; x < 2; x++ ) {
                posA( &particle[0] ) = vector3f( (float)x, (float)y, (float)z );
                velA( &particle[0] ) = vector3f( (float)( x + 1 ), (float)( y + 1 ), (float)( z + 1 ) );
                radA( &particle[0] ) = 1.f;

                pgt.insert( &particle[0] );
            }
    // check for points landed on
    mesh.add_vertex( 0.f, 0.f, 0.f );
    mesh.add_vertex( 1.f, 0.f, 0.f );
    mesh.add_vertex( 0.f, 1.f, 0.f );
    mesh.add_vertex( 1.f, 1.f, 0.f );
    mesh.add_vertex( 0.f, 0.f, 1.f );
    mesh.add_vertex( 1.f, 0.f, 1.f );
    mesh.add_vertex( 0.f, 1.f, 1.f );
    mesh.add_vertex( 1.f, 1.f, 1.f );

    // check for point in middle
    mesh.add_vertex( 0.5f, 0.5f, 0.5f );

    // check for obscure points inside
    mesh.add_vertex( 0.7f, 0.2f, 0.3f );

    // check for points outside box, but inside range
    mesh.add_vertex( -1.1f, -1.1f, -1.1f );

    // check for points outside range
    mesh.add_vertex( 3.f, 3.f, 3.f );

    voxel_coord_system vcs( vector3f(), 1.f );
    particle_union_of_spheres_is_policy isp( pgt, 1.f, 2.f, 1.f, vcs, 1 );

    channel_propagation_policy cpp( true );
    cpp.add_channel( _T("Velocity") );
    isp.populate_mesh_channels( mesh, cpp );

    // COMFIRM VERTEX CHANNELS
    trimesh3_vertex_channel_accessor<frantic::graphics::vector3f> acc =
        mesh.get_vertex_channel_accessor<frantic::graphics::vector3f>( _T("Velocity") );

    EXPECT_TRUE( ( *(vector3f*)acc.data( 0 ) ).is_equal( vector3f( 1.f, 1.f, 1.f ), 0.00001f ) ); // pos = 0,0,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 1 ) ).is_equal( vector3f( 2.f, 1.f, 1.f ), 0.00001f ) ); // pos = 1,0,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 2 ) ).is_equal( vector3f( 1.f, 2.f, 1.f ), 0.00001f ) ); // pos = 0,1,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 3 ) ).is_equal( vector3f( 2.f, 2.f, 1.f ), 0.00001f ) ); // pos = 1,1,0
    EXPECT_TRUE( ( *(vector3f*)acc.data( 4 ) ).is_equal( vector3f( 1.f, 1.f, 2.f ), 0.00001f ) ); // pos = 0,0,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 5 ) ).is_equal( vector3f( 2.f, 1.f, 2.f ), 0.00001f ) ); // pos = 1,0,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 6 ) ).is_equal( vector3f( 1.f, 2.f, 2.f ), 0.00001f ) ); // pos = 0,1,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 7 ) ).is_equal( vector3f( 2.f, 2.f, 2.f ), 0.00001f ) ); // pos = 1,1,1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 8 ) ).is_equal( vector3f( 1.f, 1.f, 1.f ), 0.00001f ) ); // pos = 0.5,0.5,0.5
    EXPECT_TRUE( ( *(vector3f*)acc.data( 9 ) ).is_equal( vector3f( 2.f, 1.f, 1.f ), 0.00001f ) ); // pos = 0.7,0.2,0.3
    EXPECT_TRUE(
        ( *(vector3f*)acc.data( 10 ) ).is_equal( vector3f( 1.f, 1.f, 1.f ), 0.00001f ) ); // pos = -1.1,-1.1,-1.1
    EXPECT_TRUE( ( *(vector3f*)acc.data( 11 ) ).is_equal( vector3f( 0.f, 0.f, 0.f ), 0.00001f ) ); // pos = 3,3,3
}
