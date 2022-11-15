// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/files/paths.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/graphics/dart_thrower.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

TEST( DartThrow, Rejection ) {
    using namespace std;
    using namespace frantic;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    cout << endl;
    boundbox3f boundBox( vector3f( -1.f ), vector3f( 1.f ) );
    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<boost::uint32_t>( _T("ID") );

    pcm.end_channel_definition();

    boost::mt19937 baseGen;

    voxel_coord_system vcs( vector3f( 0.f, 0.f, 0.f ), 0.25f );
    rejection_dart_thrower3d<boost::mt19937> rejectThrower( pcm, boundBox, 0.5f, baseGen );
    basic_dart_thrower3d<boost::mt19937> basicThrower( pcm, boundBox, 0.5f, baseGen );

    vector3f point;
    vector3f p;
    //		int nIterations = 10000;
    int fail = 0, trueCount = 0, falseCount = 0;

    float radius = 0.125f, targetDensity = 0.7f;
    std::size_t theoreticalCount = rejectThrower.estimate_particle_count( radius, targetDensity );
    size_t maxIterations = 30;

    cout << "3D Theoretical Particle Count: " << theoreticalCount << endl;

    //*
    diagnostics::profiling_section psThrowing;

    psThrowing.enter();
    particle_grid_tree rejectResults;
    rejectThrower.throw_darts( radius, targetDensity, maxIterations, rejectResults );
    psThrowing.exit();

    cout << "3D Total particle count = " << rejectResults.particle_count() << endl;
    cout << "# of Iterations Taken = " << maxIterations << endl;

    particles::particle_array particles( rejectResults.get_channel_map() );

    for( particles::particle_grid_tree::const_iterator i = rejectResults.begin(), endIter = rejectResults.end();
         i != endIter; ++i )
        particles.push_back( *i );

    channel_accessor<boost::uint32_t> id = pcm.get_accessor<boost::uint32_t>( _T("ID") );
    channel_accessor<vector3f> pos = pcm.get_accessor<vector3f>( _T("Position") );

    int result = 0;
    for( particle_array::iterator i = particles.begin(); i != particles.end(); ++i ) {
        for( particle_array::iterator j = i + 1; j != particles.end(); ++j ) {

            float d = vector3f::distance( pos.get( *i ), pos.get( *j ) );

            if( d < 2 * radius ) {
                std::cout << "BAD V " << pos.get( *i ) << " is " << d << " from:  " << pos.get( *j ) << std::endl;

                ++result;
                break;
            }
        }
    }

    ASSERT_TRUE( result == 0 );
    cout << "Test Radius: " << radius << endl;
    std::cout << "Counted Failures: " << result << endl;
    std::cout << "Execution Timing: " << psThrowing << endl;

    psThrowing.reset();
    maxIterations = 1000000;

    psThrowing.enter();
    particle_grid_tree basicResults;
    basicThrower.throw_darts( radius, targetDensity, maxIterations, basicResults );
    psThrowing.exit();

    cout << "3D Total particle count = " << basicResults.particle_count() << endl;
    cout << "# of Iterations Taken = " << maxIterations << endl;

    fail = 0;
    trueCount = 0;
    falseCount = 0;

    particles.clear();

    for( particles::particle_grid_tree::const_iterator i = basicResults.begin(), endIter = basicResults.end();
         i != endIter; ++i )
        particles.push_back( *i );

    /*
    for( int i =0; i < particles.size(); ++i){
      cout << "Particle " << i << " " << particles[i].second << endl;
    }//*/

    result = 0;
    for( particle_array::iterator i = particles.begin(); i != particles.end(); ++i ) {
        for( particle_array::iterator j = i + 1; j != particles.end(); ++j ) {

            float d = vector3f::distance( pos.get( *i ), pos.get( *j ) );

            if( d < 2 * radius ) {
                std::cout << "BAD V " << pos.get( *i ) << " is " << d << " from:  " << pos.get( *j ) << std::endl;

                ++result;
                break;
            }
        }
    }

    cout << "Test Radius: " << radius << endl;

    std::cout << "Counted Failures: " << result << endl;
    std::cout << "Execution Timing: " << psThrowing << endl;
    ASSERT_TRUE( result == 0 );
}