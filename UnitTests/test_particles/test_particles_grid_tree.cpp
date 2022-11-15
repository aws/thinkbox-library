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

// A function used by the testParticleGridTreeParticleParticleInteractions test to count the number of interactions
// The channel_map used with this function has the interaction count as its first member, hence the simple
// cast and increment is valid.
static void particle_particle_interaction_increment_count( void* /*userData*/, char* firstParticle,
                                                           char* secondParticle ) {
    ( *reinterpret_cast<boost::int32_t*>( firstParticle ) )++;
    ( *reinterpret_cast<boost::int32_t*>( secondParticle ) )++;
    //	std::cout << "Interacting particle " << *(reinterpret_cast<boost::int32_t*>(firstParticle)+1) << " with particle
    //"
    //<< *(reinterpret_cast<boost::int32_t*>(secondParticle)+1) << std::endl;
}

// This does an asymmetric increment, so that we can check that the interactions between two different particle grid
// trees call this
// function with the arguments in the right order.
static void particle_particle_interaction_asymmetric_increment_count( void* /*userData*/, char* firstParticle,
                                                                      char* secondParticle ) {
    ( *reinterpret_cast<boost::int32_t*>( firstParticle ) )++;
    ( *reinterpret_cast<boost::int32_t*>( secondParticle ) ) += 2;
    //	std::cout << "Interacting particle " << *(reinterpret_cast<boost::int32_t*>(firstParticle)+1) << " with particle
    //"
    //<< *(reinterpret_cast<boost::int32_t*>(secondParticle)+1) << std::endl;
}

// This does an asymmetric increment, so that we can check that the interactions between two different particle grid
// trees call this
// function with the arguments in the right order.
static void
particle_grid_interaction_asymmetric_increment_count( void* /*userData*/, char* firstParticle,
                                                      const frantic::graphics::vector3f& /*gridParticlePosition*/,
                                                      char* gridParticle ) {
    ( *reinterpret_cast<boost::int32_t*>( firstParticle ) )++;
    ( *reinterpret_cast<boost::int32_t*>( gridParticle ) ) += 2;
    //	std::cout << "Interacting particle " << *(reinterpret_cast<boost::int32_t*>(firstParticle)+1) << " with particle
    //"
    //<< *(reinterpret_cast<boost::int32_t*>(secondParticle)+1) << std::endl;
}

struct particle_eval_data {
    frantic::channels::channel_accessor<frantic::graphics::vector3f> posAcc;
    frantic::graphics::boundbox3f bounds;
    int numParticles;
};

static void particle_eval_func( void* userData, char* particle ) {
    particle_eval_data* data = reinterpret_cast<particle_eval_data*>( userData );
    data->numParticles++;
    EXPECT_TRUE( data->bounds.contains( data->posAcc( particle ) ) );
}

TEST( Particles_GridTree, Creation ) {
    using namespace std;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    boundbox3f bounds( vector3f( 0.f ), vector3f( 1.f ) );
    float voxelSize = 1.0f;

    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<boost::uint32_t>( _T("ID") );
    pcm.end_channel_definition();

    std::vector<char> rawParticle;
    rawParticle.resize( pcm.structure_size() );

    channel_accessor<boost::uint32_t> id = pcm.get_accessor<boost::uint32_t>( _T("ID") );
    channel_accessor<vector3f> pos = pcm.get_accessor<vector3f>( _T("Position") );

    cout << endl;

    particle_grid_tree cursorTree( pcm, 0.25f );

    particle_grid_tree::const_iterator begin = cursorTree.begin();

    EXPECT_TRUE( begin == cursorTree.end() );

    for( int i = 0; i < 4; ++i ) {
        for( int j = 0; j < 4; ++j ) {
            for( int k = 0; k < 4; ++k ) {
                id.get( &rawParticle[0] ) = i;
                pos.get( &rawParticle[0] ) = vector3f( 0.125f + i * 0.25f, 0.125f + j * 0.25f, 0.125f + k * 0.25f );
                cursorTree.insert( &rawParticle[0] );
            }
        }
    }

    pos.get( &rawParticle[0] ) = vector3f( -0.125f, 0.125f, 0.15f );
    cursorTree.insert( &rawParticle[0] );

    pos.get( &rawParticle[0] ) = vector3f( -0.375f, -0.125f, -0.15f );
    cursorTree.insert( &rawParticle[0] );

    particle_grid_tree::const_iterator iter = cursorTree.begin();
    particle_grid_tree::const_iterator end = cursorTree.end();

    EXPECT_TRUE( iter != cursorTree.end() );
    size_t count = 0;
    for( ; iter != cursorTree.end(); ++iter ) {
        ++count;
    }
    cout << "1) GridTree: Particle Count: " << (size_t)cursorTree.particle_count() << endl;
    cout << "1) Grid Tree Cursor Count: " << count << endl;

    EXPECT_TRUE( count == (size_t)cursorTree.particle_count() );

    bounds = boundbox3f( vector3f( -1.f ), vector3f( 1.f ) );
    particle_grid_tree theTree( pcm, voxelSize );

    size_t addedCount = 0;
    for( int i = 0; i < 1000; ++i ) {
        id.get( &rawParticle[0] ) = i;
        pos.get( &rawParticle[0] ) = vector3f::from_unit_random();
        theTree.insert( &rawParticle[0] );
        ++addedCount;
    }

    cout << "Random Test, Added Count: " << addedCount << endl;

    EXPECT_TRUE( (size_t)theTree.particle_count() == addedCount );

    iter = theTree.begin();

    EXPECT_TRUE( iter != theTree.end() );

    count = 0;
    for( ; iter != theTree.end(); ++iter ) {
        // cout << "Particle: " << pos.get(*iter) << endl;
        ++count;
    }

    cout << "2) GridTree: Particle Count: " << (size_t)theTree.particle_count() << endl;
    cout << "2) Grid Tree Cursor Count: " << count << endl;

    EXPECT_TRUE( count == (size_t)theTree.particle_count() );

    particle_grid_tree::node_iterator nodes_start = theTree.nodes_begin();
    particle_grid_tree::node_iterator nodes_end = theTree.nodes_end();

    size_t leafCount = 0;
    size_t sumCount = 0;
    for( ; nodes_start != nodes_end; ++nodes_start ) {
        ++leafCount;
        particle_grid_tree::node_iterator::particle_iterator_pair iterPair = nodes_start.particle_iterators();

        // iterPair

        for( ; iterPair.first != iterPair.second; ++iterPair.first ) {
            ++sumCount;
        }
    }
    cout << "leaf count: " << leafCount << endl;
    EXPECT_TRUE( sumCount == (size_t)theTree.particle_count() );

    particle_grid_tree theTree2( pcm, voxelSize );

    pos.get( &rawParticle[0] ) = vector3f( 0.385f, 0.385f, 0.f );
    theTree2.insert( &rawParticle[0] );

    EXPECT_TRUE( theTree2.has_particle_near( vector3f( 0.386f, 0.375f, 0.f ), 0.15f ) );

    iter = theTree2.begin();

    for( count = 0; iter != theTree2.end(); ++iter ) {
        // cout << "Particle: " << pos.get(*iter) << endl;
        ++count;
    }

    EXPECT_TRUE( count == 1 );

    // test the particle_grid_tree node swap
    theTree2.swap( theTree );

    cout << "Tree2 Should have " << addedCount << " particles, it really has: " << (size_t)theTree2.particle_count()
         << endl;
    cout << "Tree1 Should have 1 particle, it really has: " << (size_t)theTree.particle_count() << endl;
    EXPECT_TRUE( (size_t)theTree2.particle_count() == addedCount );
    EXPECT_TRUE( (size_t)theTree.particle_count() == 1 );

    std::vector<char> tempData;
    tempData.resize( pcm.structure_size() );
    char* tempParticle = &tempData[0];

    for( particle_grid_tree::const_iterator i = theTree2.begin(); i != theTree2.end(); ++i ) {
        pcm.copy_structure( tempParticle, *i );

        theTree.insert( tempParticle );
    }

    cout << "Tree Should have " << addedCount + 1 << " particles, it really has: " << (size_t)theTree.particle_count()
         << endl;
}

TEST( Particles_GridTree, Insert_Particles ) {
    using namespace std;
    using namespace boost;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    channel_map pcmA, pcmB;
    pcmA.define_channel<vector3f>( _T("Position") );
    pcmA.define_channel<boost::uint32_t>( _T("ID") );
    pcmA.end_channel_definition();

    pcmB.define_channel<vector3f>( _T("Position") );
    pcmB.define_channel<float>( _T("Radius") );
    pcmB.define_channel<boost::uint32_t>( _T("ID") );
    pcmB.define_channel<vector3f>( _T("Velocity") );

    pcmB.end_channel_definition();

    voxel_coord_system vcs( vector3f(), 1.f );

    particle_grid_tree treeA( pcmA, vcs );
    particle_grid_tree treeB( pcmB, vcs );

    channel_accessor<boost::uint32_t> idA = pcmA.get_accessor<boost::uint32_t>( _T("ID") );
    channel_accessor<vector3f> posA = pcmA.get_accessor<vector3f>( _T("Position") );

    channel_accessor<boost::uint32_t> idB = pcmB.get_accessor<boost::uint32_t>( _T("ID") );
    channel_accessor<vector3f> posB = pcmB.get_accessor<vector3f>( _T("Position") );
    channel_accessor<float> radB = pcmB.get_accessor<float>( _T("Radius") );
    channel_accessor<vector3f> velB = pcmB.get_accessor<vector3f>( _T("Velocity") );

    vector<char> tempA( pcmA.structure_size() );
    vector<char> tempB( pcmB.structure_size() );
    int nA = 15000;
    int nB = nA / 10;

    for( int i = 0; i < nA; ++i ) {
        idA( &tempA[0] ) = i;
        posA( &tempA[0] ) = vector3f::from_random();

        treeA.insert( &tempA[0] );
    }

    for( int i = 0; i < nB; ++i ) {
        idB( &tempB[0] ) = i;
        posB( &tempB[0] ) = vector3f::from_random();
        velB( &tempB[0] ) = vector3f::from_random();
        radB( &tempB[0] ) = rand() / (float)RAND_MAX;

        treeB.insert( &tempB[0] );
    }

    particle_grid_tree treeAintoB( pcmB, vcs );

    treeAintoB.insert_particles( treeB );
    particle_grid_tree::const_iterator iter = treeAintoB.begin(), iterEnd = treeAintoB.end();
    for( ; iter != iterEnd; ++iter ) {
        bool found = false;
        particle_grid_tree::const_iterator treeBIter = treeB.begin(), treeBIterEnd = treeB.end();
        for( ; treeBIter != treeBIterEnd && !found; ++treeBIter ) {
            if( idB( *iter ) == idB( *treeBIter ) && posB( *iter ) == posB( *treeBIter ) &&
                radB( *iter ) == radB( *treeBIter ) && velB( *iter ) == velB( *treeBIter ) ) {
                found = true;
            }
        }

        if( !found ) {
            cout << "_";

            // TS_FAIL( "Tree B was not properly inserted into the comparisiton grid tree" );
        }
    }

    cout << "A particle count: " << treeA.particle_count() << endl;
    cout << "B particle count: " << treeAintoB.particle_count() << endl;
    treeAintoB.insert_particles( treeA );
    cout << "TOTAL particle count: " << treeAintoB.particle_count() << endl;

    iter = treeAintoB.begin();
    iterEnd = treeAintoB.end();

    for( ; iter != iterEnd; ++iter ) {
        bool found = false;

        particle_grid_tree::const_iterator treeBIter = treeB.begin(), treeBIterEnd = treeB.end();
        for( ; treeBIter != treeBIterEnd && !found; ++treeBIter ) {
            if( idB( *iter ) == idB( *treeBIter ) && posB( *iter ) == posB( *treeBIter ) &&
                radB( *iter ) == radB( *treeBIter ) && velB( *iter ) == velB( *treeBIter ) ) {
                found = true;
            }
        }

        if( found ) {
            continue;
        }

        bool foundB = false;
        particle_grid_tree::const_iterator treeAIter = treeA.begin(), treeAIterEnd = treeA.end();
        for( ; treeAIter != treeAIterEnd && !foundB; ++treeAIter ) {
            if( idB( *iter ) == idA( *treeAIter ) && posB( *iter ) == posA( *treeAIter ) && radB( *iter ) == 0.f &&
                velB( *iter ) == vector3f() ) {
                foundB = true;
            }
        }
        // Particle in Tree A was not properly inserted into the comparisiton grid tree
        EXPECT_TRUE( foundB );
    }

    particle_grid_tree treeBintoA( pcmA, vcs );

    treeBintoA.insert_particles( treeA );
    cout << "A particle count: " << treeA.particle_count() << endl;
    cout << "A particle count: " << treeBintoA.particle_count() << endl;
    cout << "B particle count: " << treeB.particle_count() << endl;
    treeBintoA.insert_particles( treeB );
    cout << "TOTAL particle count: " << treeBintoA.particle_count() << endl;

    iter = treeBintoA.begin();
    iterEnd = treeBintoA.end();

    for( ; iter != iterEnd; ++iter ) {
        bool found = false;
        particle_grid_tree::const_iterator treeAIter = treeA.begin(), treeAIterEnd = treeA.end();

        /*	if( treeAIterEnd == treeA.end())
            cout << "@" << endl; */

        // cout << "i_ptr" << *treeAIter << endl;
        // cout << "ie_ptr" << *treeAIterEnd << endl;
        // cout << "\n";
        for( int c = 0; treeAIter != treeAIterEnd && !found; ++treeAIter, ++c ) {
            // cout << c << " ";
            if( idA( *iter ) == idA( *treeAIter ) && posA( *iter ) == posA( *treeAIter ) ) {
                found = true;
            }
        }

        //*
        if( found ) {
            continue;
        } //*/

        bool foundB = false;
        particle_grid_tree::const_iterator treeBIter = treeB.begin(), treeBIterEnd = treeB.end();
        for( ; treeBIter != treeBIterEnd && !foundB; ++treeBIter ) {
            if( idA( *iter ) == idB( *treeBIter ) && posA( *iter ) == posB( *treeBIter ) ) {
                foundB = true;
            }
        }
        // Particle in Tree B was not properly inserted into the comparisiton grid tree
        EXPECT_TRUE( foundB );
    }
}

TEST( Particles_GridTree, Particle_Particle_Interaction ) {
    using namespace std;
    using namespace boost;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    tbb::task_scheduler_init taskScheduleInit;

    // Put the interaction count first, so that our increment and decrement functions don't have to know a non-zero
    // offset.
    channel_map pcm;
    pcm.define_channel<boost::int32_t>( _T("InteractionCount") );
    pcm.define_channel<boost::int32_t>( _T("ID") );
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.end_channel_definition( 4, true );

    channel_accessor<boost::int32_t> countAcc = pcm.get_accessor<boost::int32_t>( _T("InteractionCount") );
    channel_accessor<boost::int32_t> idAcc = pcm.get_accessor<boost::int32_t>( _T("ID") );
    channel_accessor<vector3f> posAcc = pcm.get_accessor<vector3f>( _T("Position") );

    boundbox3f bounds( vector3f(), vector3f( 11 ) );
    int particleCount = 3000;

    particle_grid_tree pgt( pcm, voxel_coord_system() );

    struct {
        boost::int32_t interactionCount;
        boost::int32_t id;
        vector3f position;
    } tempParticle;

    tempParticle.interactionCount = 0;
    tempParticle.id = 0;

    // Generate 27 particles, one in each of the 27 voxels of the bounding box, all of them adjacent to the center
    // particle
    /*
    for( int z = 0; z < 3; ++z ) {
      for( int y = 0; y < 3; ++y ) {
        for( int x = 0; x < 3; ++x ) {
          //if( x < 2 && y > 0 && z < 2 ) {
            tempParticle.position.x = (x-1) * 1.1f / 2 + 1.5f;
            tempParticle.position.y = (y-1) * 1.1f / 2 + 1.5f;
            tempParticle.position.z = (z-1) * 1.1f / 2 + 1.5f;
            cout << "Particle " << tempParticle.id << ": " << tempParticle.position << endl;
            pgt.insert( &tempParticle );
            ++tempParticle.id;
          //}
        }
      }
    }
    //*/
    //*
    // Generate a bunch of randomly positioned particles with interaction count 0 to start
    for( int i = 0; i < particleCount; ++i ) {
        tempParticle.position = bounds.random_vector();
        pgt.insert( (char*)&tempParticle );
        ++tempParticle.id;
    }
    //*/
    EXPECT_EQ( pgt.size(), particleCount );

    cout << "Generated a particle_grid_tree with " << particleCount << " particles, and volume " << bounds.volume()
         << " voxels" << endl;

    //////////////////////////
    // Test particle-particle interactions within one particle grid tree.
    //////////////////////////

    // This function will count the number of interactions each particle partakes in
    pgt.particle_particle_interactions( particle_particle_interaction_increment_count, 0, 1 );

    int maxCount = 0, minCount = ( numeric_limits<int>::max )(), totalCount = 0;
    for( particle_grid_tree::const_iterator i = pgt.begin(), iterEnd = pgt.end(); i != iterEnd; ++i ) {
        int count = countAcc.get( *i );
        if( count > maxCount )
            maxCount = count;
        if( count < minCount )
            minCount = count;
        totalCount += count;
    }

    // The additions were always in pairs, so the total interaction count should be even
    EXPECT_EQ( totalCount % 2, 0 );
    cout << "With an interaction radius of 1" << endl;
    cout << " Minimum interaction count: " << minCount << endl;
    cout << " Maximum interaction count: " << maxCount << endl;
    cout << " Twice the total interaction count: " << totalCount << endl;
    cout << " Average interactions per particle: " << ( float( totalCount ) / pgt.size() ) << endl;

    // Now, let's subtract out all those interactions with the brute force method, to make sure it matches up
    for( particle_grid_tree::iterator i = pgt.begin(), iterEnd = pgt.end(); i != iterEnd; ++i ) {
        particle_grid_tree::iterator j = i;
        ++j;
        for( ; j != iterEnd; ++j ) {

            if( vector3f::distance_squared( posAcc.get( *i ), posAcc.get( *j ) ) < 1 ) {
                // NOTE: It's really important to not be modifying the position here, we're just modifying the counts.
                countAcc.get( *i )--;
                countAcc.get( *j )--;
            }
        }
    }

    int interactionFailureCount = 0;
    // Check that all the interaction counts are back at zero now
    for( particle_grid_tree::const_iterator i = pgt.begin(), iterEnd = pgt.end(); i != iterEnd; ++i ) {
        // cout << "Position: " << posAcc.get(*i) << endl;
        if( countAcc.get( *i ) != 0 )
            interactionFailureCount++;
        // TS_ASSERT_EQUALS( countAcc.get(*i), 0 );
    }
    EXPECT_EQ( interactionFailureCount, 0 );
}

TEST( Particles_GridTree, Particles_In_Range ) {
    using namespace std;
    using namespace boost;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.end_channel_definition( 1, true );

    channel_accessor<vector3f> posAcc = pcm.get_accessor<vector3f>( _T("Position") );

    boundbox3f bounds( vector3f(), vector3f( 11 ) );
    int particleCount = 3000;

    particle_grid_tree pgt( pcm, voxel_coord_system() );

    struct {
        vector3f position;
    } tempParticle;

    // Generate a bunch of randomly positioned particles
    for( int i = 0; i < particleCount; ++i ) {
        tempParticle.position = bounds.random_vector();
        pgt.insert( (char*)&tempParticle );
    }

    vector<char*> nearbyParticles;
    // Test get_particles_in_range from a bunch of different positions
    for( int i = 0; i < 20; ++i ) {
        nearbyParticles.clear();
        vector3f testPos = bounds.random_vector();
        pgt.get_particles_in_range( testPos, 1, nearbyParticles );
        // Make sure all the particles we got were within range
        for( unsigned j = 0; j < nearbyParticles.size(); ++j )
            EXPECT_LT( vector3f::distance_squared( testPos, posAcc( nearbyParticles[j] ) ), 1 );

        // Count how many particles are nearby the slow way
        int count = 0;
        for( particle_grid_tree::iterator j = pgt.begin(); j != pgt.end(); ++j ) {
            if( vector3f::distance_squared( testPos, posAcc( *j ) ) < 1 )
                ++count;
        }

        // Make sure that the counts match
        EXPECT_EQ( count, (int)nearbyParticles.size() );
    }
}

TEST( Particles_GridTree, GetParticlesInRangeBox ) {
    using namespace std;
    using namespace boost;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    channel_map pcm;
    // Offset the Position to make sure that I handle it
    pcm.define_channel( _T("Position"), 3, data_type_float32, 4 );
    pcm.end_channel_definition( 1, false, false );

    channel_accessor<vector3f> posAcc = pcm.get_accessor<vector3f>( _T("Position") );

    const boundbox3f bounds( vector3f(), vector3f( 11 ) );
    const int particleCount = 300;

    particle_grid_tree pgt( pcm, voxel_coord_system() );

    std::vector<char> tempParticle( pcm.structure_size() );

    // Generate a bunch of randomly positioned particles
    for( int i = 0; i < particleCount; ++i ) {
        posAcc( tempParticle ) = bounds.random_vector();
        pgt.insert( &tempParticle[0] );
    }

    const float queryRadius = 1;

    vector<char*> nearbyParticles;

    // Test get_particles_in_range with different positions and box sizes.
    for( int i = 0; i < 100; ++i ) {
        nearbyParticles.clear();

        const boundbox3f queryBounds = boundbox3f::from_random( bounds );

        pgt.get_particles_in_range( queryBounds, queryRadius, nearbyParticles );

        // Make sure all the particles we got were within range
        for( std::size_t j = 0; j < nearbyParticles.size(); ++j ) {
            const vector3f position = posAcc( nearbyParticles[j] );
            if( !queryBounds.contains( position ) ) {
                const vector3f nearestPoint = queryBounds.get_nearest_surface_point( position );
                const float distanceToBounds = vector3f::distance_squared( position, nearestPoint );
                EXPECT_LE( distanceToBounds, queryRadius );
            }
        }

        // Count how many particles are nearby the slow way
        std::size_t count = 0;
        for( particle_grid_tree::iterator j = pgt.begin(); j != pgt.end(); ++j ) {
            const vector3f position = posAcc( *j );
            if( queryBounds.contains( position ) ) {
                ++count;
            } else {
                const vector3f nearestPoint = queryBounds.get_nearest_surface_point( position );
                const float distanceToBounds = vector3f::distance_squared( position, nearestPoint );
                if( distanceToBounds <= queryRadius ) {
                    ++count;
                }
            }
        }

        // Make sure that the counts match
        EXPECT_EQ( count, nearbyParticles.size() );
    }
}

TEST( Particles_GridTree, Process_Particles_In_Bounds ) {
    using namespace std;
    using namespace boost;
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::volumetrics;
    using namespace frantic::channels;

    particle_eval_data data;
    srand( (unsigned int)time( NULL ) );

    // Put the interaction count first, so that our increment and decrement functions don't have to know a non-zero
    // offset.
    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.end_channel_definition( 1, true );

    channel_accessor<vector3f> posAcc = pcm.get_accessor<vector3f>( _T("Position") );

    // Bounds for generating random particles and then random non empty bounding boxes
    boundbox3f bounds( vector3f(), vector3f( 11 ) );
    boundbox3f lowerBounds( vector3f(), vector3f( 5 ) );
    boundbox3f upperBounds( vector3f( 6 ), vector3f( 11 ) );

    int particleCount = 3000;

    particle_grid_tree pgt( pcm, voxel_coord_system() );

    struct {
        vector3f position;
    } tempParticle;

    // Generate a bunch of randomly positioned particles
    for( int i = 0; i < particleCount; ++i ) {
        tempParticle.position = bounds.random_vector();
        pgt.insert( (char*)&tempParticle );
    }

    vector<char*> nearbyParticles;
    // Test get_particles_in_range for a bunch of different bounding boxes
    for( int i = 0; i < 1; ++i ) {

        // get some bounding box
        boundbox3f searchBounds( lowerBounds.random_vector(), upperBounds.random_vector() );

        data.posAcc = pcm.get_accessor<vector3f>( _T("Position") );
        data.bounds = searchBounds;
        data.numParticles = 0;

        // process the particles in the bounds, the eval function takes
        // care of the test
        pgt.process_particles_in_bounds( &data, searchBounds, &particle_eval_func );

        // Count how many particles are in the bounds the slow way
        int count = 0;
        for( particle_grid_tree::iterator j = pgt.begin(); j != pgt.end(); ++j ) {
            if( posAcc( *j ).x <= searchBounds.maximum().x && posAcc( *j ).x >= searchBounds.minimum().x &&
                posAcc( *j ).y <= searchBounds.maximum().y && posAcc( *j ).y >= searchBounds.minimum().y &&
                posAcc( *j ).z <= searchBounds.maximum().z && posAcc( *j ).z >= searchBounds.minimum().z )
                ++count;
        }

        // Make sure that the counts match
        EXPECT_EQ( count, data.numParticles );
    }
}
