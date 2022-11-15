// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <vector>

#include "gtest/gtest.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>

#include <tbb/task_scheduler_init.h>

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/channel_map_lerp.hpp>
#include <frantic/graphics/dart_thrower.hpp>
#include <frantic/particles/particle_classes.hpp>
#include <frantic/particles/particle_cursor.hpp>
#include <frantic/particles/particle_kdtree.hpp>
#include <frantic/particles/particle_kmeans.hpp>
#include <frantic/particles/proxy_particle_cursor.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/geometry/particle_collision_detector.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>

#include <frantic/graphics/dart_thrower.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_grid_tree.hpp>

#include <frantic/files/paths.hpp>

#include <frantic/particles/streams/prt_particle_istream.hpp>
#include <frantic/particles/streams/prt_particle_ostream.hpp>
#include <frantic/particles/streams/pts_particle_istream.hpp>

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/particles/particle_grid_tree.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

TEST( Particles, KMeans ) {

    std::cout << "\nK-means clustering test" << std::endl;

    frantic::channels::channel_map pcm;
    pcm.define_channel<double>( _T("X") );
    pcm.define_channel<double>( _T("Y") );
    pcm.define_channel<double>( _T("Z") );
    pcm.end_channel_definition();

    frantic::particles::particle_array pa( pcm );

    pa.resize( 4 );
    *(double*)pa[0] = -6.f;
    *( (double*)pa[0] + 1 ) = -6.f;
    *( (double*)pa[0] + 2 ) = -6.f;
    *(double*)pa[1] = -4.f;
    *( (double*)pa[1] + 1 ) = -4.f;
    *( (double*)pa[1] + 2 ) = -4.f;
    *(double*)pa[2] = 4.f;
    *( (double*)pa[2] + 1 ) = 4.f;
    *( (double*)pa[2] + 2 ) = 4.f;
    *(double*)pa[3] = 6.f;
    *( (double*)pa[3] + 1 ) = 6.f;
    *( (double*)pa[3] + 2 ) = 6.f;

    //;std::cout << "particles:" << std::endl;
    // for ( int i = 0; i < 4; ++i )
    //	std::cout << *(double*)pa[i] << "," << *((double*)pa[i]+1) << "," << *((double*)pa[i]+1) << std::endl;

    frantic::particles::particle_array ma( pcm );

    ma.resize( 2 );
    *(double*)ma[0] = -6.f;
    *( (double*)ma[0] + 1 ) = -6.f;
    *( (double*)ma[0] + 2 ) = -6.f;
    *(double*)ma[1] = 6.f;
    *( (double*)ma[1] + 1 ) = 6.f;
    *( (double*)ma[1] + 2 ) = 6.f;

    std::vector<int> clusters;

    int K = 2;
    int N = 3;
    frantic::particles::particle_kmeans( pa, N, K, clusters, ma );

    std::cout << "k-means cluster centres" << std::endl;
    for( int i = 0; i < K; ++i ) {
        std::cout << *(double*)ma[i] << ", " << *( (double*)ma[i] + 1 ) << ", " << *( (double*)ma[i] + 1 ) << std::endl;
    }

    float tol = 1e-5f;
    for( int i = 0; i < 3; ++i )
        // Test to see if it finds cluster centre 1
        EXPECT_LT( fabs( *( (double*)ma[0] + i ) + 5.f ), tol );

    for( int i = 0; i < 3; ++i )
        // Test to see if it finds cluster centre 2
        EXPECT_LT( fabs( *( (double*)ma[1] + i ) - 5.f ), tol );
}

TEST( Particle, Classes ) {
    using namespace frantic::particles;
    // Check that the particle classes are actually the size we expect
    EXPECT_EQ( sizeof( basic_particle ), 26 );
    EXPECT_EQ( sizeof( renderable_particle ), 20 );
    EXPECT_EQ( sizeof( motion_blurred_particle ), 26 );
    EXPECT_EQ( sizeof( colored_particle ), 32 );
}

TEST( Particle, ContainerClasses ) {
    using namespace std;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::channels;

    std::vector<char> particleData;
    std::vector<boost::uint32_t> proxyParticleIDs;

    channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<vector3f>( _T("Velocity") );
    channelMap.define_channel<vector3f>( _T("Acceleration") );
    channelMap.define_channel<half>( _T("U") );
    channelMap.define_channel<half>( _T("V") );
    channelMap.define_channel<boost::uint32_t>( _T("ID") );
    channelMap.end_channel_definition();

    // Make some data
    int n_particles = 10;
    for( unsigned i = 0; i < channelMap.structure_size() * n_particles; ++i ) {
        particleData.push_back( 0 );
    }

    particle_cursor pCursor( &particleData[0], &particleData[particleData.size() - 1], channelMap.structure_size() );

    channel_accessor<vector3f> pos = channelMap.get_accessor<vector3f>( _T("Position") );
    channel_accessor<vector3f> vel = channelMap.get_accessor<vector3f>( _T("Velocity") );
    channel_accessor<vector3f> acc = channelMap.get_accessor<vector3f>( _T("Acceleration") );
    channel_accessor<half> u = channelMap.get_accessor<half>( _T("U") );
    channel_accessor<half> v = channelMap.get_accessor<half>( _T("V") );
    channel_accessor<boost::uint32_t> id = channelMap.get_accessor<boost::uint32_t>( _T("ID") );

    //
    EXPECT_THROW( channel_accessor<half> temp = channelMap.get_accessor<half>( _T("Position") ), std::runtime_error )
        << "Trying to access a channel with the incorrect type";

    // Type Checking
    for( int i = 0; i < n_particles && pCursor.next_particle(); ++i ) {
        vector3f& p = pos.get( pCursor.raw_particle_buffer() );
        p.x = p.y = p.z = float( 50 * i );
        u.get( pCursor.raw_particle_buffer() ) = 5.f;
        ;
        v.get( pCursor.raw_particle_buffer() ) = 6.f;

        //*u1 = (float) 5.0;
        //*v1 = (float) 6.0;
        id.get( pCursor.raw_particle_buffer() ) = (boost::uint32_t)i;

        if( i % 3 == 0 ) {
            proxyParticleIDs.push_back( i );
        }

        if( i % 2 == 0 ) {
            vector3f& p1 = vel.get( pCursor.raw_particle_buffer() );
            p1.x = 10;
            p1.y = 10;
            p1.z = 10;
        }
        if( i % 3 == 0 ) {
            vector3f& p2 = acc.get( pCursor.raw_particle_buffer() );
            p2.x = -1;
            p2.y = -1;
            p2.z = -1;
        }
    }

    const channel_map& constChannelMap = channelMap;
    // output
    channel_accessor<vector3f> c_pos = constChannelMap.get_accessor<vector3f>( _T("Position") );
    channel_accessor<vector3f> c_vel = constChannelMap.get_accessor<vector3f>( _T("Velocity") );
    channel_accessor<vector3f> c_acc = constChannelMap.get_accessor<vector3f>( _T("Acceleration") );
    channel_accessor<half> c_u = constChannelMap.get_accessor<half>( _T("U") );
    channel_accessor<half> c_v = constChannelMap.get_accessor<half>( _T("V") );
    channel_accessor<boost::uint32_t> c_id = constChannelMap.get_accessor<boost::uint32_t>( _T("ID") );

    pCursor.reset();

    for( int i = 0; i < n_particles && pCursor.next_particle(); ++i ) {
        const vector3f& p = c_pos.get( pCursor.raw_particle_buffer() );
        const vector3f& velc = c_vel.get( pCursor.raw_particle_buffer() );
        const half& u1 = c_u.get( pCursor.raw_particle_buffer() );
        const half& v1 = c_v.get( pCursor.raw_particle_buffer() );
        const vector3f& a = c_acc.get( pCursor.raw_particle_buffer() );

        EXPECT_EQ( p.x, 50 * i );
        EXPECT_EQ( p.y, 50 * i );
        EXPECT_EQ( p.z, 50 * i );

        EXPECT_EQ( c_id.get( pCursor.raw_particle_buffer() ), i );

        float testVel = 0.f;
        if( i % 2 == 0 )
            testVel = 10.f;

        EXPECT_EQ( velc.x, testVel );
        EXPECT_EQ( velc.y, testVel );
        EXPECT_EQ( velc.z, testVel );

        float testAcc = 0.0f;
        if( i % 3 == 0 )
            testAcc = -1.f;

        EXPECT_EQ( a.x, testAcc );
        EXPECT_EQ( a.y, testAcc );
        EXPECT_EQ( a.z, testAcc );

        EXPECT_EQ( u1, 5 );
        EXPECT_EQ( v1, 6 );
    }
    proxy_particle_cursor proxyCursor( &particleData[0], &particleData[0] + particleData.size(), proxyParticleIDs,
                                       channelMap );

    //		cout << "\nParticle IDs of the Proxy Particles:\n";
    //		for(size_t i=0; i < proxyParticleIDs.size(); ++i){
    //			cout << proxyParticleIDs[i] << " ";
    //		}
    //		cout << "\n\nExtracted Proxy Particle IDs:\n";

    while( proxyCursor.next_particle() ) {
        vector3f& p = pos.get( proxyCursor.raw_particle_buffer() );
        vector3f& velc = vel.get( proxyCursor.raw_particle_buffer() );
        half& u1 = u.get( proxyCursor.raw_particle_buffer() );
        half& v1 = v.get( proxyCursor.raw_particle_buffer() );
        vector3f& a = acc.get( proxyCursor.raw_particle_buffer() );

        boost::uint32_t theID = id.get( proxyCursor.raw_particle_buffer() );
        //			cout << theID << " " ;

        // Access the variables to get rid of the warning messages
        p;
        velc;
        u1;
        v1;
        a;
        theID;
    }

    while( proxyCursor.prev_particle() ) {
        vector3f& p = pos.get( proxyCursor.raw_particle_buffer() );
        vector3f& velc = vel.get( proxyCursor.raw_particle_buffer() );
        half& u1 = u.get( proxyCursor.raw_particle_buffer() );
        half& v1 = v.get( proxyCursor.raw_particle_buffer() );
        vector3f& a = acc.get( proxyCursor.raw_particle_buffer() );

        // Access the variables to get rid of the warning messages
        p;
        velc;
        u1;
        v1;
        a;

        id.get( proxyCursor.raw_particle_buffer() ) = id.get( proxyCursor.raw_particle_buffer() ) | PRT_FLG_MODIFIED;
    }

    const_proxy_particle_cursor c_proxyCursor(
        &particleData[0], &particleData[0] + n_particles * channelMap.structure_size(), proxyParticleIDs, channelMap );

    cout << "\n\nConst Extracted Proxy Particle IDs:\n";
    {
        //			unsigned i = 0;
        while( c_proxyCursor.next_particle() ) {
            const vector3f& p = c_pos.get( c_proxyCursor.raw_particle_buffer() );
            const vector3f& velc = c_vel.get( c_proxyCursor.raw_particle_buffer() );
            const half& u1 = c_u.get( c_proxyCursor.raw_particle_buffer() );
            const half& v1 = c_v.get( c_proxyCursor.raw_particle_buffer() );
            const vector3f& a = c_acc.get( c_proxyCursor.raw_particle_buffer() );

            const boost::uint32_t theID = c_id.get( c_proxyCursor.raw_particle_buffer() ) & ~PRT_FLG_MODIFIED;
            //				if( i < proxyParticleIDs.size() && masked_uint32<PRT_ID_MASK >(c_id.get(c_proxyCursor))
            //==
            // masked_uint32<PRT_ID_MASK>( proxyParticleIDs[i])) {
            //					cout << theID << " " ;
            //					++i;
            //				}

            // Access the variables to get rid of the warning messages
            p;
            velc;
            u1;
            v1;
            a;
            theID;
        }
    }

    // particle_cursor p_cursor(&particleData[0],&particleData[particleData.size()-1], channelMap.structure_size() );

    // cout << boost::lexical_cast<std::string>(pa.get(p_cursor));
    //		cout << "Test Complete\n";
}
TEST( Particle, Array ) {
    using namespace std;
    using namespace boost;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<boost::uint32_t>( _T("ID") );
    pcm.end_channel_definition();

    particle_array parray( pcm );

    channel_accessor<boost::uint32_t> id = pcm.get_accessor<boost::uint32_t>( _T("ID") );
    channel_accessor<vector3f> pos = pcm.get_accessor<vector3f>( _T("Position") );

    vector<char> particle( pcm.structure_size() );
    size_t count = 100000;
    boost::uint32_t i;

    for( i = 0; i < count; ++i ) {
        vector3f p = vector3f::from_unit_random();
        p = p * (float)( rand() % 10 );

        pos.get( &particle[0] ) = p;
        id.get( &particle[0] ) = i;

        parray.push_back( &particle[0] );
    }

    parray.remove_channel( _T("Position") );

    EXPECT_TRUE( parray.get_channel_map().channel_count() == 1 );

    parray.set_channel_map( pcm );

    EXPECT_TRUE( parray.get_channel_map().channel_count() == 2 );

    particle_array::iterator iter = parray.begin();

    for( i = 0; iter != parray.end(); ++iter, ++i ) {
        EXPECT_TRUE( id.get( *iter ) == i );
    }

    particle_array::const_iterator const_iter = parray.begin();

    for( i = 0; const_iter != parray.end(); ++const_iter, ++i ) {
        EXPECT_TRUE( id.get( *const_iter ) == i );
    }

    EXPECT_TRUE( parray.particle_count() == count );

    channel_map pcm2;
    pcm2.define_channel<boost::uint32_t>( _T("ID") );
    pcm2.end_channel_definition();

    particle_array parray2( pcm2 );

    for( i = 0; i < count; ++i ) {
        id.get( &particle[0] ) = (int)count + i;
        parray2.push_back( &particle[0] );
    }

    parray.copy_particles_from( parray2 );

    EXPECT_TRUE( parray.particle_count() == 2 * count );
}

TEST( Particle, GeometryCollisions ) {
    using namespace std;
    using namespace frantic::geometry;
    using namespace frantic::particles;

    particle_collision_detector pcd;
    // Fill in the pcd with a static box and a moving box
    trimesh3 box;
    box.set_to_box( boundbox3f( vector3f( 0 ), vector3f( 1 ) ) );
    pcd.add_static_object( transform4f::from_translation( 3, 3, 3 ), box );
    motion_blurred_transform<float> mbt( transform4f::from_translation( 0, 0, 0 ),
                                         transform4f::from_translation( 1, 0, 0 ) );
    pcd.add_rigid_object( mbt, box );

    pcd.prepare_kdtrees();
    // pcd.print_statistics(std::cout);

    raytrace_intersection isect;

    // Miss everything by a mile
    EXPECT_TRUE( !pcd.collides_with_particle( vector3f( 10, 4, 8 ), vector3f( 11, 12, 4 ), isect ) );

    // Find a collision with the static box
    EXPECT_TRUE( pcd.collides_with_particle( vector3f( 3.5f, 3.5f, 2 ), vector3f( 3.5f, 3.5f, 4 ), isect ) );
    EXPECT_LT( fabs( isect.distance - 0.5 ), 0.000001f );
    EXPECT_LT( vector3f::distance( isect.position, vector3f( 3.5f, 3.5f, 3 ) ), 0.000001f );

    // Find a collision with the moving box
    EXPECT_TRUE( pcd.collides_with_particle( vector3f( 1, 0.5f, 1.5f ), vector3f( 1, 0.5f, 0.5f ), isect ) );
    EXPECT_LT( fabs( isect.distance - 0.5 ), 0.000001f );
    EXPECT_LT( vector3f::distance( isect.position, vector3f( 1, 0.5f, 1 ) ), 0.000001f );

    // Find a collision with a moving box
    EXPECT_TRUE( pcd.collides_with_particle( vector3f( 2, 0, 1 ), vector3f( 1, 1, 0 ), isect ) );
    EXPECT_LT( fabs( isect.distance - 0.5 ), 0.000001f );
    EXPECT_LT( vector3f::distance( isect.position, vector3f( 1.5f, 0.5f, 0.5f ) ), 0.000001f );

    // Narrowly miss the moving box
    EXPECT_FALSE( pcd.collides_with_particle( vector3f( -0.1f, 0.5f, 0.5f ), vector3f( 0.9f, 0.95f, 0.5f ), isect ) );
    EXPECT_FALSE( pcd.collides_with_particle( vector3f( 1.1f, 0.5f, 0.5f ), vector3f( 2.1f, 0.5f, 0.5f ), isect ) );
    EXPECT_FALSE( pcd.collides_with_particle( vector3f( 2, 0.5f, 0.5f ), vector3f( 1, 1.6f, 0.5f ), isect ) );
}
