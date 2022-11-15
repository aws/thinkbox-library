// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <vector>

#include "gtest/gtest.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scope_exit.hpp>

#include <tbb/task_scheduler_init.h>

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/channel_map_lerp.hpp>
#include <frantic/graphics/dart_thrower.hpp>
#include <frantic/particles/particle_classes.hpp>
#include <frantic/particles/particle_cursor.hpp>
#include <frantic/particles/particle_kdtree.hpp>
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

TEST( Particle_KDTree, Creation ) {
    using namespace std;
    using namespace frantic::graphics;

    frantic::particles::particle_kdtree<vector3f> kd;

    kd.add_particle( vector3f( 0, 0, 0 ) );
    kd.add_particle( vector3f( 0, 0, 1 ) );
    kd.add_particle( vector3f( 0, 1, 0 ) );
    kd.add_particle( vector3f( 0, 1, 1 ) );
    kd.add_particle( vector3f( 1, 0, 0 ) );
    kd.add_particle( vector3f( 1, 0, 1 ) );
    kd.add_particle( vector3f( 1, 1, 0 ) );
    kd.add_particle( vector3f( 1, 1, 1 ) );

    kd.balance_kdtree();

    //		cout << "Balanced kd-tree:" << endl;
    //		for( unsigned i = 0; i < kd.size(); ++i )
    //			cout << "particle " << i << ": " << kd[i] << endl;

    EXPECT_EQ( (int)kd.size(), 8 );

    std::vector<std::pair<float, vector3f>> nearestParticles;

    // Should be able to find each particle individually
    for( unsigned i = 0; i < kd.size(); ++i ) {
        kd.locate_particles_range( kd[i], 0.5f * 0.5f, nearestParticles );
        EXPECT_EQ( (int)nearestParticles.size(), 1 );
        if( nearestParticles.size() > 0 )
            EXPECT_EQ( kd[i], nearestParticles[0].second );
        nearestParticles.clear();
    }

    // Should find 4 particles
    kd.locate_particles_range( vector3f( 0, 1, 0 ), 1.1f * 1.1f, nearestParticles );
    EXPECT_EQ( (int)nearestParticles.size(), 4 );
    nearestParticles.clear();

    // Should find no particles
    kd.locate_particles_range( vector3f( 0.5f, 0.5f, 0.3f ), 0.3f * 0.3f, nearestParticles );
    EXPECT_EQ( (int)nearestParticles.size(), 0 );
    nearestParticles.clear();

    // Should find all particles
    kd.locate_particles_range( vector3f( 0.5f, 0.5f, 0.5f ), 2.5f * 2.5f, nearestParticles );
    EXPECT_EQ( (int)nearestParticles.size(), 8 );
    //		for( unsigned i = 0; i < nearestParticles.size(); ++i )
    //			cout << "particle " << i << ": " << nearestParticles[i].second << endl;
}

// This tests the support function in the particle kd-tree
TEST( Particle_KDTree, Support ) {
    using namespace std;
    using namespace frantic::graphics;

    frantic::particles::particle_kdtree<vector3f> kd;
    int pointCount = 5000, testCount = 100;

    // First create a random set of points to test
    for( int i = 0; i < pointCount; ++i ) {
        kd.add_particle( vector3f::from_random_gaussian() );
    }

    frantic::diagnostics::profiling_section pf_kdSupport( _T("kd-tree support") ),
        pf_manSupport( _T("manual support") ), pf_balance( _T("kd-tree balance") );

    pf_balance.enter();
    kd.balance_kdtree();
    pf_balance.exit();

    cerr << pf_balance << endl;

    // Now test the support function with a bunch of different directions.
    for( int test = 0; test < testCount; ++test ) {
        vector3f d = vector3f::from_random_gaussian();
        pf_kdSupport.enter();
        vector3f kdSupport = kd.support( d );
        pf_kdSupport.exit();
        pf_manSupport.enter();
        vector3f manSupport = kd[0];
        float manDotProd = vector3f::dot( manSupport, d );
        for( unsigned i = 1; i < kd.size(); ++i ) {
            vector3f testSupport = kd[i];
            float testDotProd = vector3f::dot( testSupport, d );
            if( testDotProd > manDotProd ) {
                manSupport = testSupport;
                manDotProd = testDotProd;
            }
        }
        pf_manSupport.exit();
        if( kdSupport != manSupport ) {
            cerr << "kd-tree support: " << kdSupport << endl;
            cerr << "naive support: " << manSupport << endl;
        }
        EXPECT_EQ( kdSupport, manSupport );
    }
    cerr << pf_kdSupport << endl;
    cerr << pf_manSupport << endl;
}