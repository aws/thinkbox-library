// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics2d/boundrect2f.hpp>
#include <frantic/particles/particle_cursor.hpp>
#include <frantic/particles/particle_grid_tree.hpp>
#include <frantic/particles/streams/particle_ostream.hpp>

#include <boost/random.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

// extern std::ofstream fout;

#pragma warning( push )
#pragma warning( disable : 4127 )

namespace frantic {
namespace graphics {

class default_rejection_policy {
    const boundbox3f* m_bounds;

  public:
    default_rejection_policy() {}

    default_rejection_policy( const boundbox3f* bounds )
        : m_bounds( bounds ) {}

    bool accept_particle( const vector3f& particlePosition ) const { return m_bounds->contains( particlePosition ); }

    std::size_t estimate_particle_count( float particlesPerUnitVolume ) const {
        return ( std::size_t )( m_bounds->volume() * particlesPerUnitVolume );
    }
};

class level_set_rejection_policy {
    const volumetrics::levelset::rle_level_set* m_clipLevelSet;
    const level_set_rejection_policy* m_parentPolicy;
    float m_upperBound, m_lowerBound;
    boundbox3f m_worldBounds;

  public:
    frantic::diagnostics::profiling_section psRejectTree;
    frantic::diagnostics::profiling_section psRejectLevelSet;

    level_set_rejection_policy() {}

    /**
     * Creates a level set based rejection policy for a dart thrower
     *
     * @param parentPolicy - This policy will be checked before the created policy. Use this to nest policies
     */
    level_set_rejection_policy( const volumetrics::levelset::rle_level_set* clipLevelSet, float upperBound = 0,
                                float lowerBound = -std::numeric_limits<float>::max(),
                                const level_set_rejection_policy* parentPolicy = 0 )
        : m_clipLevelSet( clipLevelSet ) {
        if( upperBound < lowerBound )
            throw std::runtime_error( "level_set_rejection_policy: upper bound " +
                                      boost::lexical_cast<std::string>( upperBound ) + " is less than lower bound " +
                                      boost::lexical_cast<std::string>( lowerBound ) );

        m_upperBound = upperBound;
        m_lowerBound = lowerBound;

        if( clipLevelSet ) {
            m_worldBounds = m_clipLevelSet->get_voxel_coord_system().get_world_bounds(
                m_clipLevelSet->get_rle_index_spec().outer_bounds() );
            FF_LOG( debug ) << "\t\tWorld Bounds: " << m_worldBounds << std::endl;
        }

        if( parentPolicy ) {
            FF_LOG( debug ) << "\t\tParent Policy provided on construction" << std::endl;
            m_parentPolicy = parentPolicy;
        } else {
            FF_LOG( debug ) << "\t\tNo Parent Policy provided on construction" << std::endl;
            m_parentPolicy = 0;
        }
    }

    std::size_t estimate_particle_count( float particlesPerUnitVolume ) const {
        if( m_clipLevelSet )
            return ( std::size_t )( m_clipLevelSet->compute_volume() * particlesPerUnitVolume );
        else
            return 0;
    }

    bool accept_particle( const vector3f& particlePosition ) const {
        if( m_parentPolicy ) {
            if( !m_parentPolicy->accept_particle( particlePosition ) )
                return false;
        }

        if( m_clipLevelSet ) {
            if( m_worldBounds.contains( particlePosition ) ) {
                float phi = m_clipLevelSet->trilerp_signed_distance( particlePosition );
                // FF_LOG(debug) << "Testing phi=" << phi << " against upper: " << m_upperBound << " and lower: " <<
                // m_lowerBound << std::endl;
                return phi < m_upperBound && phi > m_lowerBound;
            } else {
                // FF_LOG(debug) << "\t\t" << particlePosition << " is outside of the world bounds " << m_worldBounds <<
                // std::endl;
                //  if the exterior code is inside (-1), then return true,
                return m_clipLevelSet->get_rle_index_spec().get_exterior_region_code() < -1;
            }
        } else
            return false;
    }

    level_set_rejection_policy& operator=( const level_set_rejection_policy& rhs ) {
        m_clipLevelSet = rhs.m_clipLevelSet;
        m_upperBound = rhs.m_upperBound;
        m_lowerBound = rhs.m_lowerBound;
        m_parentPolicy = rhs.m_parentPolicy;
        m_worldBounds = rhs.m_worldBounds;

        return *this;
    }
};

class particle_grid_tree_rejection_policy {
    const particles::particle_grid_tree* m_particleTree;
    float m_testRadius;

  public:
    frantic::diagnostics::profiling_section psRejectTree;
    frantic::diagnostics::profiling_section psRejectLevelSet;

    particle_grid_tree_rejection_policy() {}

    particle_grid_tree_rejection_policy( const particles::particle_grid_tree* particleTree, float radius )
        : m_particleTree( particleTree ) {
        // m_particleTree = particleTree;
        m_testRadius = radius;
        // std::cout << "set radius: " << m_testRadius << std::endl;
    }

    bool accept_particle( const vector3f& particlePosition ) const {
        return !m_particleTree->has_particle_near( particlePosition, 2 * m_testRadius );
    }

    std::size_t estimate_particle_count( float particlesPerUnitVolume ) const {
        boundbox3f bounds = m_particleTree->compute_particle_bounds();
        return ( std::size_t )( bounds.volume() * particlesPerUnitVolume ) - m_particleTree->particle_count();
    }

    particle_grid_tree_rejection_policy& operator=( const particle_grid_tree_rejection_policy& rhs ) {
        m_particleTree = rhs.m_particleTree;
        m_testRadius = rhs.m_testRadius;
        return *this;
    }
};

class combined_rejection_policy {
    level_set_rejection_policy m_levelSetPolicy;
    particle_grid_tree_rejection_policy m_particlePolicy;

  public:
    frantic::diagnostics::profiling_section psRejectTree;
    frantic::diagnostics::profiling_section psRejectLevelSet;

    combined_rejection_policy() {
        psRejectTree.set_name( _T("Rejection:Tree") );
        psRejectLevelSet.set_name( _T("Rejection:LevelSet") );
    }

    combined_rejection_policy( particle_grid_tree_rejection_policy particlePolicy,
                               level_set_rejection_policy m_particlePolicy )
        : m_levelSetPolicy( m_particlePolicy )
        , m_particlePolicy( particlePolicy ) {}

    bool
    accept_particle( const vector3f& particlePosition ) /*const*/ { // disabling const for timing tests, should be const
        psRejectLevelSet.enter();
        bool levelSetResult = m_levelSetPolicy.accept_particle( particlePosition );
        psRejectLevelSet.exit();

        bool treeResult = false;

        if( levelSetResult ) {
            psRejectTree.enter();
            treeResult = m_particlePolicy.accept_particle( particlePosition );
            psRejectTree.exit();
        }

        return treeResult;
    }

    std::size_t estimate_particle_count( float particlesPerUnitVolume ) const {
        size_t lsEstimate = m_levelSetPolicy.estimate_particle_count( particlesPerUnitVolume );
        size_t prtEstimate = m_particlePolicy.estimate_particle_count( particlesPerUnitVolume );

        return ( lsEstimate < prtEstimate ) ? lsEstimate : prtEstimate;
    }

    combined_rejection_policy& operator=( const combined_rejection_policy& rhs ) {
        m_levelSetPolicy = rhs.m_levelSetPolicy;
        m_particlePolicy = rhs.m_particlePolicy;
        return *this;
    }
};

// NOTE: THIS CLASS IS NOT CURRENTLY USED ANYWHERE
template <class Generator>
class basic_dart_thrower3d {
    frantic::graphics::boundbox3f m_boundBox;
    float m_voxelLength;
    frantic::channels::channel_map m_channelMap;

    Generator m_baseGen;
    boost::uniform_real<float> m_disType;
    boost::variate_generator<Generator, boost::uniform_real<float>> m_rnd;

    boost::shared_ptr<frantic::particles::streams::particle_istream> m_preload;

  public:
    basic_dart_thrower3d( const frantic::channels::channel_map& pcm, frantic::graphics::boundbox3f bounds,
                          float voxelLength, Generator& theGen )
        : m_baseGen( theGen )
        , m_disType( 0.f, 1.f )
        , m_rnd( m_baseGen, m_disType )
        , m_channelMap( pcm ) {
        m_voxelLength = voxelLength;
        m_boundBox = bounds;
    }

    void set_preload_stream( boost::shared_ptr<frantic::particles::streams::particle_istream> pin ) { m_preload = pin; }

    void set_bound_box( graphics::boundbox3f bounds ) { m_boundBox = bounds; }

    graphics::boundbox3f get_bound_box() { return m_boundBox; }

    size_t estimate_particle_count( float radius, float relativePacking ) {
        double fourSqrtTwo = 4 * std::sqrt( 2.0 );
        // because
        boundbox3f realBox( m_boundBox.minimum() - vector3f( radius ), m_boundBox.maximum() + vector3f( radius ) );

        // float boxVolume = ( (m_boundBox.xsize() +  2*radius)*(m_boundBox.ysize() +  2*radius)*(m_boundBox.zsize() +
        // 2*radius) );
        return (size_t)( realBox.volume() * relativePacking * relativePacking * relativePacking /
                         ( fourSqrtTwo * radius * radius * radius ) );
    }

    void throw_darts( float radius, float relativePacking, std::size_t& maxIterations,
                      particles::particle_grid_tree& outDarts, bool debugPretend = false ) {
        // ensure that the particle grid tree is setup properly for the dart thrower
        outDarts.reset( m_channelMap, m_voxelLength );
        boost::uint32_t insertedCount = 0;

        if( m_preload ) {
            outDarts.insert_particles( m_preload );
        }

        boost::uint32_t count = insertedCount;
        // we offset the theoretical count to reflect that we are starting the count the insertedCount instead of at 0
        size_t theoreticalCount = estimate_particle_count( radius, relativePacking );

        unsigned i;
        std::vector<char> rawParticle( m_channelMap.structure_size() );

        // channels::channel_accessor<boost::uint32_t> id = m_channelMap.get_accessor<boost::uint32_t>("ID");
        channels::channel_accessor<vector3f> pos = m_channelMap.get_accessor<vector3f>( _T("Position") );

        vector3f min = m_boundBox.minimum();
        vector3f diff = m_boundBox.maximum() - min;

        diagnostics::profiling_section psNear( _T("PoissonSearch") );
        for( i = 0; i < maxIterations; ++i ) {
            // vector3f dart( m_rnd()*m_boundBox.maximum().x , m_rnd() *m_boundBox.maximum().y,
            // m_rnd()*m_boundBox.maximum().z
            // );
            vector3f dart( m_rnd(), m_rnd(), m_rnd() );
            dart.x *= diff.x;
            dart.y *= diff.y;
            dart.z *= diff.z;
            dart = min + dart;

            psNear.enter();
            if( !outDarts.has_particle_near( dart, 2 * radius ) ) {

                if( !debugPretend ) {
                    // id.get(&rawParticle[0]) = count;
                    pos.get( &rawParticle[0] ) = dart;

                    outDarts.insert( &rawParticle[0] );
                }

                count++;

                if( (int)count >= theoreticalCount ) {
                    break;
                }
            }
            psNear.exit();
        }

        std::cout << psNear << "\n";

        maxIterations = i; // return the number of iterations actually taken

        float particleVolume = 4 * (float)frantic::math::pi_value * radius * radius * radius / 3.f;
        float totalVolume = ( count * particleVolume );
        float boxVolume = ( ( m_boundBox.xsize() + 2 * radius ) * ( m_boundBox.ysize() + 2 * radius ) *
                            ( m_boundBox.zsize() + 2 * radius ) );
        float actualPacking = totalVolume / ( 0.7405f * boxVolume );

        std::cout << "Particle Vol    : " << particleVolume << "  Total Volume  : " << totalVolume
                  << "  Box Volume : " << boxVolume << std::endl;
        std::cout << "Maximum Possible Packing Count: " << estimate_particle_count( radius, 1.f ) << std::endl;
        std::cout << "Expected Packing^3: " << relativePacking * relativePacking * relativePacking
                  << " Actual Packing^3: " << actualPacking << std::endl;

        std::cout << "Count: " << count << " expected: " << theoreticalCount << std::endl;
        std::cout << "Internal Particle Count: " << outDarts.particle_count() << std::endl;

        std::cout << "Number of Iterations: " << i << std::endl;
        std::cout << "Max Iterations: " << maxIterations << std::endl;

        maxIterations = i; // return the number of iterations actually taken

        // TODO: Need to return the particles somehow
        // return darts;
    } //*/
};

// the default rejection policy tests only against the given throwing bounds
// other policies could test against levelsets, clipping planes, etc...
template <class Generator, class RejectionPolicy = default_rejection_policy>
class rejection_dart_thrower3d {
    graphics::boundbox3f m_boundBox;
    float m_voxelLength;
    frantic::channels::channel_map m_channelMap;

    Generator m_baseGen;
    boost::uniform_real<float> m_disType;
    boost::variate_generator<Generator, boost::uniform_real<float>> m_rand;

    boost::shared_ptr<frantic::particles::streams::particle_istream> m_preload;

    RejectionPolicy m_rejectPolicy;

    struct available_dart {
        vector3f point;
        std::size_t tryCount;
        std::vector<vector3f> placedNeighbors;
    };

    int m_genCount;
    std::vector<std::vector<available_dart>>
        availableDarts; // a vector of avialable darts, where each base index is a generation

  public:
    rejection_dart_thrower3d( const frantic::channels::channel_map& pcm, frantic::graphics::boundbox3f bounds,
                              float voxelLength, RejectionPolicy reject, Generator& theGen )
        : m_baseGen( theGen )
        , m_disType( 0.f, 1.f )
        , m_rand( m_baseGen, m_disType )
        , m_channelMap( pcm ) {
        m_voxelLength = voxelLength;
        m_boundBox = bounds;
        m_rejectPolicy = reject;

        m_genCount = 2;

        for( int generation = 0; generation < m_genCount; ++generation ) {
            availableDarts.push_back( std::vector<available_dart>() ); // push back generation one.
        }
    }

    rejection_dart_thrower3d( const frantic::channels::channel_map& pcm, frantic::graphics::boundbox3f bounds,
                              float voxelLength, Generator& theGen )
        : m_baseGen( theGen )
        , m_disType( 0.f, 1.f )
        , m_rand( m_baseGen, m_disType )
        , m_channelMap( pcm ) {
        m_voxelLength = voxelLength;
        m_boundBox = bounds;
        m_rejectPolicy = RejectionPolicy( &m_boundBox );

        m_genCount = 2;

        for( int generation = 0; generation < m_genCount; ++generation ) {
            availableDarts.push_back( std::vector<available_dart>() ); // push back generation one.
        }
    }

    void set_preload_particles( boost::shared_ptr<frantic::particles::streams::particle_istream> pin ) {
        m_preload = pin;
    }

    void add_available_particle( const vector3f& pos ) {
        available_dart temp;
        temp.point = pos;
        temp.tryCount = 0;
        availableDarts[0].push_back( temp );
    }

    // sets the bound box. This is important as the thrower only throws within the bound box
    void set_bound_box( graphics::boundbox3f bounds ) { m_boundBox = bounds; }

    // returns the bound box of the thrower
    graphics::boundbox3f get_bound_box() const { return m_boundBox; }

    // estimates the particle count based on the radius and the requested relative packing
    std::size_t estimate_particle_count( float radius, float relativePacking ) const {
        double fourSqrtTwo = 4 * std::sqrt( 2.0 );
        // because
        // boundbox3f realBox( m_boundBox.minimum() - vector3f(radius), m_boundBox.maximum() + vector3f(radius) );

        // float boxVolume = ( (m_boundBox.xsize() +  2*radius)*(m_boundBox.ysize() +  2*radius)*(m_boundBox.zsize() +
        // 2*radius) );
        float particlesPerUnitVol =
            (float)( relativePacking * relativePacking * relativePacking / ( fourSqrtTwo * radius * radius * radius ) );
        return m_rejectPolicy.estimate_particle_count( particlesPerUnitVol );
    }

    // attempts to throw darts of the given radius in the boundbox and achieve the requested relative packing
    void throw_darts( float radius, float relativePacking, std::size_t& maxIterations,
                      particles::particle_grid_tree& outDarts ) {

        frantic::diagnostics::profiling_section psTreeReset( _T("Particle Grid Tree Reset") );
        psTreeReset.enter();
        outDarts.reset( m_channelMap, m_voxelLength );
        psTreeReset.exit();

        std::size_t insertedCount = 0;
        size_t maxAttemptCount = maxIterations;
        float innerRadius = radius * 2.f;
        float outerRadius = innerRadius * 2.f;
        float innerRadiusSqrd = innerRadius * innerRadius;

        // availablepoint, iteration count, placed neighbors

        vector3f dart;

        std::vector<char> rawParticle( m_channelMap.structure_size() );

        // channels::channel_accessor<boost::uint32_t> id = m_channelMap.get_accessor<boost::uint32_t>("ID");
        channels::channel_accessor<vector3f> pos = m_channelMap.get_accessor<vector3f>( _T("Position") );

        frantic::diagnostics::profiling_section psPreLoad( _T("Preload Particles") );

        if( m_preload ) {
            psPreLoad.enter();
            outDarts.insert_particles( m_preload ); // TODO This is a problem, we aren'\t using the preload particles as
                                                    // seed available particles
            psPreLoad.exit();
        }

        /*
        if( m_preloadParticles.size() > 0 ) {
          for( unsigned i = 0; i < m_preloadParticles.size();  ++i) {
            id.get(&rawParticle[0]) = (boost::uint32_t) insertedCount;
            pos.get(&rawParticle[0]) = m_preloadParticles[i];

            outDarts.insert( &rawParticle[0] );

            availablePoints.push_back(std::make_pair(m_preloadParticles[i], 0) );
            ++insertedCount;
          }
        }//*/

        boost::variate_generator<Generator, boost::uniform_on_sphere<float>> sphereRand(
            m_baseGen, boost::uniform_on_sphere<float>( 3 ) );

        vector3f boundsMin = m_boundBox.minimum();
        vector3f boundsMax = m_boundBox.maximum();

        // FF_LOG(debug) << "\t\tAvailable generations: " << availableDarts.size() << std::endl;
        // for(int generation=0; generation < m_genCount; ++generation) {
        //	availableDarts[generation].reserve(100000);
        // } // try to do some pre-allocation

        std::size_t theoreticalCount = estimate_particle_count( radius, relativePacking );

        std::size_t seedCount = theoreticalCount / 10;
        if( seedCount < 50 )
            seedCount = 50;

        // FF_LOG(debug) << "\t\tThrower:Target Bound Box: " << m_boundBox << "\n";
        // FF_LOG(debug) << "\t\tThrower:Initial Seed Count: " << seedCount << endl;

        frantic::diagnostics::profiling_section psInitSeeding( _T("InitialSeeding") );
        frantic::diagnostics::profiling_section psInsert( _T("Insertion") );
        frantic::diagnostics::profiling_section psThrow( _T("Throwing") );
        frantic::diagnostics::profiling_section psSphereRand( _T("Spherical Random Number") );

        frantic::diagnostics::profiling_section psRejection( _T("RejectionPolicy") );
        frantic::diagnostics::profiling_section psNear( _T("PoissonSearch") );
        frantic::diagnostics::profiling_section psNeighbor( _T("NeighbourCache") );
        frantic::diagnostics::profiling_section psRemoveAttempt( _T("Removing Failed Attempt") );

        psInitSeeding.enter();

        if( availableDarts[0].size() == 0 ) {
            for( unsigned c = 0; c < seedCount; ++c ) {
                // add the first 100 completely random seed points to the available list
                dart = m_boundBox.random_vector( m_rand );

                // FF_LOG(debug) << "Trying dart " << c << ": " << dart << endl;

                if( m_rejectPolicy.accept_particle( dart ) ) {
                    // FF_LOG(debug) << "\t\t\tPassed policy test! " << endl;

                    // if( !outDarts.has_particle_near( dart, innerRadius ) ){
                    pos.get( &rawParticle[0] ) = dart;
                    if( outDarts.insert_with_proximity_constraint( &rawParticle[0], innerRadius ) ) {
                        //	if( outDarts.insert( &rawParticle[0] ) ){
                        available_dart temp;
                        temp.point = dart;
                        temp.tryCount = 0;
                        temp.placedNeighbors.reserve( 2 );
                        availableDarts[0].push_back( temp );
                        ++insertedCount;

                        //	}
                    }
                }
            }
        }

        psInitSeeding.exit();

        // FF_LOG(debug) << "\t\tThrower:Available Count Gen 0: " << availableDarts[0].size() << "\n";
        // FF_LOG(debug) << "\t\tThrower:Particle Count: " << outDarts.particle_count() << "\n";
        // FF_LOG(debug) << "\t\tThrower:Max attempts per particle: " << maxAttemptCount << std::endl;

        // return;

        vector3f offset( outerRadius );
        std::size_t count = insertedCount;

        std::vector<char> rejectedParticles;

        size_t promotedCount = 0;
        size_t iterationCount = 0;
        available_dart temp;

        temp.tryCount = 0;
        temp.placedNeighbors.reserve( 2 );

        for( int generation = 0; generation < m_genCount /*&& count < theoreticalCount*/; ++generation ) {

            size_t generationParticleLimit = generation ? ( generation * 3 ) : 1;
            // FF_LOG(debug) << "\tSeeding Generation " << generation << std::endl;
            // FF_LOG(debug) << "\tAvailable Seeding " << availableDarts[generation].size() << std::endl;

            while( !availableDarts[generation].empty() ) {
                psThrow.enter();

                // we use the std integer rand() function since we are grabbing a random index rather than a floating
                // point value
                size_t p = rand() % availableDarts[generation].size();

                psSphereRand.enter();
                // TODO change the from unit random spherical to use a [0,1] distribution instead of [-1,1]
                // vector3f random = innerRadius * vector3f::from_unit_random_spherical( sphereRand );
                std::vector<float> pointOnSphere = sphereRand();
                vector3f random = innerRadius * vector3f( pointOnSphere[0], pointOnSphere[1], pointOnSphere[2] );

                psSphereRand.exit();

                // In order to save the pointer additions we keep pointer to the struct.
                // There is a DANGER here in that the the available darts list if changed could invalidate the pointer
                available_dart* centerPoint = &availableDarts[generation][p];

                // dart = availableDarts[generation][p].point + random*1.001f ;
                dart = centerPoint->point + random * 1.001f;

                bool succeededThrow = false;
                if( m_boundBox.contains( dart ) ) {

                    //						FF_LOG(debug) << "Generation: " << generation << "\n";
                    //						FF_LOG(debug) << availableDarts[generation][p].point << "\t" <<  dart <<
                    //"\t"
                    //<< availableDarts[generation][p].tryCount << "\t" << outDarts.particle_count() << endl;

                    psRejection.enter();
                    if( m_rejectPolicy.accept_particle( dart ) ) {
                        psRejection.exit();
                        //*
                        psNeighbor.enter();
                        // check the other darts thrown from the same seed and see if they can give an early rejection
                        std::size_t neighborSize = centerPoint->placedNeighbors.size();
                        for( std::size_t i = 0; i < neighborSize; ++i ) {
                            vector3f testDart = centerPoint->placedNeighbors[i];

                            float d = vector3f::distance_squared( dart, testDart );

                            if( d < innerRadiusSqrd ) {
                                ++centerPoint->tryCount;
                                psNeighbor.exit();
                                continue;
                            }
                        }
                        psNeighbor.exit();
                        //*/

                        pos.get( &rawParticle[0] ) = dart;
                        psNear.enter();
                        if( outDarts.insert_with_proximity_constraint( &rawParticle[0], innerRadius ) ) {
                            psNear.exit();
                            psInsert.enter();
                            succeededThrow = true;

                            centerPoint->placedNeighbors.push_back( dart );

                            // promote the seed dart if needed
                            if( centerPoint->placedNeighbors.size() > generationParticleLimit ) {
                                // move the dart to the back to make it easier to remove it
                                std::swap( *centerPoint,
                                           availableDarts[generation][availableDarts[generation].size() - 1] );
                                // add it to the next oldest generation
                                if( generation + 1 < m_genCount )
                                    availableDarts[generation + 1].push_back( availableDarts[generation].back() );
                                // remove it from the younger generation
                                availableDarts[generation].pop_back();
                                // std::cout << "promoting particle...\n";
                                ++promotedCount;
                            }

                            temp.point = dart;

                            availableDarts[generation].push_back( temp );

                            count++;
                            psInsert.exit();
                        } //*/
                        else {
                            psNear.exit();
                        }
                    } else {
                        psRejection.exit();
                    }
                }

                if( !succeededThrow ) {
                    if( centerPoint->tryCount >= maxAttemptCount ) {
                        // FF_LOG(debug) << "Removing particle " << availableDarts[generation][p].point << endl;

                        psRemoveAttempt.enter();

                        // remove the particle p from the available list
                        if( p != availableDarts[generation].size() - 1 )
                            std::swap( *centerPoint,
                                       availableDarts[generation][availableDarts[generation].size() - 1] );

                        availableDarts[generation].pop_back();

                        psRemoveAttempt.exit();

                        // pos.get(&rawParticle[0]) = dart;
                        // std::copy(rawParticle.begin(), rawParticle.end(), std::back_inserter( rejectedParticles ) );

                    } else { // The particle was rejected, but not yet at its attempt limit
                        ++centerPoint->tryCount;

                        // pos.get(&rawParticle[0]) = dart;
                        // std::copy(rawParticle.begin(), rawParticle.end(), std::back_inserter( rejectedParticles ) );
                    }
                }

                ++iterationCount;
                psThrow.exit();
                // if(  count >= theoreticalCount ) {
                //	break;
                // }
            }
        }

        // FF_LOG(debug) << "\t\tGen 1 Promotion Count: " << promotedCount << "\n";

        // std::ofstream fout("c:\\temp\\darts_log.txt");
        // FF_LOG(debug) << "Iterations: " << iterationCount << "\n";

        // FF_LOG(debug) << "Remaining Test Points: " << youngAvailableDarts.size() << "\n";
        FF_LOG( debug ) << "\t*** " << psTreeReset << "\n";
        FF_LOG( debug ) << "\t*** " << psPreLoad << "\n";
        FF_LOG( debug ) << "\t*** " << psInitSeeding << "\n";
        FF_LOG( debug ) << "\t*** " << psInsert << "\n";
        FF_LOG( debug ) << "\t*** " << psRejection << "\n";
        // FF_LOG(debug) << "\t\t*** " << m_rejectPolicy.psRejectTree << std::endl;
        // FF_LOG(debug) << "\t\t*** " << m_rejectPolicy.psRejectLevelSet << std::endl;
        FF_LOG( debug ) << "\t*** " << psSphereRand << "\n";
        FF_LOG( debug ) << "\t*** " << psNear << "\n";
        FF_LOG( debug ) << "\t*** " << psRemoveAttempt << "\n";
        FF_LOG( debug ) << "\t*** " << psNeighbor << "\n";
        FF_LOG( debug ) << "\t*** " << psThrow << std::endl;

        /* verbose outputting...
        float particleVolume = 4 * (float) frantic::math::pi_value * radius * radius * radius / 3.f;
        float totalVolume = (count * particleVolume ) ;
        float boxVolume = ( (m_boundBox.xsize() +  2*radius)*(m_boundBox.ysize() +  2*radius)*(m_boundBox.zsize() +
        2*radius) ); float actualPacking = totalVolume/(0.7405f * boxVolume);


        std::cout << "Particle Vol    : " << particleVolume<< "  Total Volume  : " << totalVolume << "  Box Volume : "
        << boxVolume << std::endl; std::cout << "Maximum Possible Packing Count: " <<
        estimate_particle_count(radius,1.f) << std::endl; std::cout << "Expected Packing^3: " <<
        relativePacking*relativePacking*relativePacking << " Actual Packing^3: " << actualPacking << std::endl;

        std::cout << "Count: " << count << std::endl;
        std::cout << "AvailablePoint Size: " << availablePoints.size() << std::endl;
        std::cout << "Internal Particle Count: " << darts->particle_count() << std::endl;

        std::cout << "Number of Iterations: " << iterationCount << std::endl;
        std::cout << "Max Iterations: " << maxIterations << std::endl;
        //*/

        for( int generation = 0; generation < m_genCount; ++generation )
            availableDarts[generation].clear();
    }
};

// the default rejection policy tests only against the given throwing bounds
// other policies could test against levelsets, clipping planes, etc...
// NOTE: THIS CLASS IS NOT CURRENTLY USED ANYWHERE
template <class Generator, class RejectionPolicy = default_rejection_policy>
class generation_limit_dart_thrower3d {
    graphics::boundbox3f m_boundBox;
    float m_voxelLength;
    frantic::channels::channel_map m_channelMap;

    Generator m_baseGen;
    boost::uniform_real<float> m_disType;
    boost::variate_generator<Generator, boost::uniform_real<float>> m_rand;

    boost::shared_ptr<frantic::particles::streams::particle_istream> m_preload;

    RejectionPolicy m_rejectPolicy;

    struct available_dart {
        vector3f point;
        std::size_t tryCount;
        std::vector<vector3f> placedNeighbors;
    };

    int m_genCount;
    std::vector<std::vector<available_dart>>
        availableDarts; // a vector of avialable darts, where each base index is a generation

  public:
    generation_limit_dart_thrower3d( const frantic::channels::channel_map& pcm, graphics::boundbox3f bounds,
                                     float voxelLength, RejectionPolicy reject, Generator& theGen )
        : m_baseGen( theGen )
        , m_disType( 0.f, 1.f )
        , m_rand( m_baseGen, m_disType )
        , m_channelMap( pcm ) {
        m_voxelLength = voxelLength;
        m_boundBox = bounds;
        m_rejectPolicy = reject;

        m_genCount = 2;

        for( int generation = 0; generation < m_genCount; ++generation ) {
            availableDarts.push_back( std::vector<available_dart>() ); // push back generation one.
        }
    }

    generation_limit_dart_thrower3d( const frantic::channels::channel_map& pcm, graphics::boundbox3f bounds,
                                     float voxelLength, Generator& theGen )
        : m_baseGen( theGen )
        , m_disType( 0.f, 1.f )
        , m_rand( m_baseGen, m_disType )
        , m_channelMap( pcm ) {
        m_voxelLength = voxelLength;
        m_boundBox = bounds;
        m_rejectPolicy = RejectionPolicy( &m_boundBox );

        m_genCount = 2;

        for( int generation = 0; generation < m_genCount; ++generation ) {
            availableDarts.push_back( std::vector<available_dart>() ); // push back generation one.
        }
    }

    void set_preload_particles( boost::shared_ptr<frantic::particles::streams::particle_istream> pin ) {
        m_preload = pin;
    }

    void add_available_particle( const vector3f& pos ) {
        available_dart temp;
        temp.point = pos;
        temp.tryCount = 0;
        availableDarts[0].push_back( temp );
        // std::cout <<  "Adding available Particle: " << pos << std::endl;
    }

    // sets the bound box. This is important as the thrower only throws within the bound box
    void set_bound_box( graphics::boundbox3f bounds ) { m_boundBox = bounds; }

    void set_generation_limit( int genLimit ) {
        m_genCount = genLimit;

        while( (int)availableDarts.size() < m_genCount ) {
            availableDarts.push_back( std::vector<available_dart>() ); // push back generation one.
        }
    }

    // returns the bound box of the thrower
    graphics::boundbox3f get_bound_box() const { return m_boundBox; }

    // estimates the particle count based on the radius and the requested relative packing
    std::size_t estimate_particle_count( float radius, float relativePacking ) const {
        double fourSqrtTwo = 4 * std::sqrt( 2.0 );
        // because
        // boundbox3f realBox( m_boundBox.minimum() - vector3f(radius), m_boundBox.maximum() + vector3f(radius) );

        // float boxVolume = ( (m_boundBox.xsize() +  2*radius)*(m_boundBox.ysize() +  2*radius)*(m_boundBox.zsize() +
        // 2*radius) );
        float particlesPerUnitVol =
            (float)( relativePacking * relativePacking * relativePacking / ( fourSqrtTwo * radius * radius * radius ) );
        return m_rejectPolicy.estimate_particle_count( particlesPerUnitVol );
    }

    // attempts to throw darts of the given radius in the boundbox and achieve the requested relative packing
    void throw_darts( float radius, float relativePacking, std::size_t& maxIterations,
                      particles::particle_grid_tree& outDarts ) {

        outDarts.reset( m_channelMap, m_voxelLength );
        std::size_t insertedCount = 0;
        size_t maxAttemptCount = maxIterations;
        float innerRadius = radius * 2.f;
        float outerRadius = innerRadius * 2.f;
        float innerRadiusSqrd = innerRadius * innerRadius;

        // availablepoint, iteration count, placed neighbors

        // no space to throw into
        if( m_boundBox.is_empty() ) {
            FF_LOG( debug ) << "No space in bound box to throw darts: " << m_boundBox << std::endl;
            return;
        }

        vector3f dart;

        std::vector<char> rawParticle( m_channelMap.structure_size() );

        // channels::channel_accessor<boost::uint32_t> id = m_channelMap.get_accessor<boost::uint32_t>("ID");
        channels::channel_accessor<vector3f> pos = m_channelMap.get_accessor<vector3f>( "Position" );
        channels::channel_accessor<color3f> col = m_channelMap.get_accessor<color3f>( "Color" );

        if( m_preload ) {
            outDarts.insert_particles( m_preload ); // TODO This is a problem, we aren'\t using the preload particles as
                                                    // seed available particles
        }

        /*
        if( m_preloadParticles.size() > 0 ) {
          for( unsigned i = 0; i < m_preloadParticles.size();  ++i) {
            id.get(&rawParticle[0]) = (boost::uint32_t) insertedCount;
            pos.get(&rawParticle[0]) = m_preloadParticles[i];

            outDarts.insert( &rawParticle[0] );

            availablePoints.push_back(std::make_pair(m_preloadParticles[i], 0) );
            ++insertedCount;
          }
        }//*/

        boost::uniform_on_sphere<float> sphereBoost( 3 );
        // boost::variate_generator<boost::mt19937,boost::uniform_real<float> > sphereRand(m_baseGen,
        // boost::uniform_real<float>(-1.f,1.f) );

        boost::variate_generator<Generator, boost::uniform_on_sphere<float>> sphereRand( m_baseGen, sphereBoost );

        vector3f boundsMin = m_boundBox.minimum();
        vector3f boundsMax = m_boundBox.maximum();

        // std::ofstream fout ("c:\\temp\\darts.log");

        // FF_LOG(debug) << "Available darts size: " << availableDarts.size() << std::endl;
        // for(int generation=0; generation < m_genCount; ++generation) {
        //	availableDarts[generation].reserve(10000);
        // } // try to do some pre-allocation

        std::size_t theoreticalCount = estimate_particle_count( radius, relativePacking );

        available_dart temp;

        temp.tryCount = 0;
        temp.placedNeighbors.reserve( 2 );

        color3f genColor;

        std::size_t seedCount = theoreticalCount / 10;

        if( availableDarts[0].size() == 0 ) {
            FF_LOG( debug ) << "seeding available darts" << std::endl;
            for( unsigned c = 0; c < seedCount; ++c ) {
                // add the first 100 completely random seed points to the available list
                dart = m_boundBox.random_vector( m_rand );

                //					FF_LOG(debug) << "Trying dart " << c << ": " << dart << endl;

                if( m_rejectPolicy.accept_particle( dart ) ) {
                    pos.get( &rawParticle[0] ) = dart;
                    if( outDarts.insert_with_proximity_constraint( &rawParticle[0], innerRadius ) ) {
                        temp.point = dart;
                        availableDarts[0].push_back( temp );
                        ++insertedCount;

                        //	}
                    }
                }
            }
        }

        FF_LOG( debug ) << "\tAvailable Count: " << availableDarts[0].size() << std::endl;
        FF_LOG( debug ) << "\tParticle Count: " << outDarts.particle_count() << std::endl;
        FF_LOG( debug ) << "\tMax attempts per particle: " << maxAttemptCount << std::endl;

        // return;

        vector3f offset( outerRadius );
        std::size_t count = insertedCount;

        frantic::diagnostics::profiling_section psInsert( _T("Insertion") );
        frantic::diagnostics::profiling_section psThrow( _T("Throwing") );
        frantic::diagnostics::profiling_section psBounds( _T("BoundsTesting") );
        frantic::diagnostics::profiling_section psNear( _T("PoissonSearch") );
        frantic::diagnostics::profiling_section psNeighbor( _T("NeighbourCache") );

        std::vector<char> rejectedParticles;

        size_t promotedCount = 0;
        size_t iterationCount = 0;

        genColor.g = 0.15f;

        for( int generation = 0; generation < m_genCount && count < theoreticalCount; ++generation ) {

            // size_t generationParticleLimit = generation ? (generation * 3) : 1;
            if( generation == 0 )
                genColor = color3f( 0.f, 1.f, 0.f );
            else if( generation == 1 )
                genColor = color3f( 0.f, 0.f, 1.f );
            else if( generation == 2 )
                genColor = color3f( 1.f, 1.f, 1.f );

            col.get( &rawParticle[0] ) = genColor;

            while( !availableDarts[generation].empty() ) {
                psThrow.enter();

                // we use the std integer rand() function since we are grabbing a random index rather than a floating
                // point value
                size_t p = rand() % availableDarts[generation].size();

                psBounds.enter();

                // vector3f random = innerRadius * vector3f::from_unit_random_spherical( sphereRand );
                const std::vector<float>& randomResult = sphereRand();

                vector3f random( innerRadius * randomResult[0], innerRadius * randomResult[1],
                                 innerRadius * randomResult[2] );

                // In order to save the pointer additions we keep pointer to the struct.
                // There is a DANGER here in that the the available darts list if changed could invalidate the pointer
                available_dart* centerPoint = &availableDarts[generation][p];

                // dart = availableDarts[generation][p].point + random*1.001f ;
                dart = centerPoint->point + random * 1.001f;

                bool succeededThrow = false;
                if( m_boundBox.contains( dart ) ) {

                    //						FF_LOG(debug) << "Generation: " << generation << "\n";
                    //						FF_LOG(debug) << availableDarts[generation][p].point << "\t" <<  dart <<
                    //"\t"
                    //<< availableDarts[generation][p].tryCount << "\t" << outDarts.particle_count() << endl;

                    if( m_rejectPolicy.accept_particle( dart ) ) {

                        psNeighbor.enter();
                        // check the other darts thrown from the same seed and see if they can give an early rejection
                        std::size_t neighborSize = centerPoint->placedNeighbors.size();
                        for( std::size_t i = 0; i < neighborSize; ++i ) {
                            vector3f testDart = centerPoint->placedNeighbors[i];

                            float d = vector3f::distance_squared( dart, testDart );

                            if( d < innerRadiusSqrd ) {
                                ++centerPoint->tryCount;
                                psNeighbor.exit();
                                continue;
                            }
                        }
                        psNeighbor.exit();

                        pos.get( &rawParticle[0] ) = dart;

                        psNear.enter();
                        if( outDarts.insert_with_proximity_constraint( &rawParticle[0], innerRadius ) ) {
                            psNear.exit();
                            psInsert.enter();
                            succeededThrow = true;

                            centerPoint->placedNeighbors.push_back( dart );

                            // remove the seed dart if needed
                            // if( centerPoint->placedNeighbors.size() > generationParticleLimit ) {
                            //	// move the dart to the back to make it easier to remove it
                            //	std::swap( *centerPoint,
                            // availableDarts[generation][availableDarts[generation].size()-1]);

                            //	availableDarts[generation].pop_back();
                            //}

                            temp.point = dart;
                            // we insert the new dart in the next generation
                            if( generation + 1 < m_genCount )
                                availableDarts[generation + 1].push_back( temp );

                            count++;
                            psInsert.exit();
                            psBounds.exit();

                        } //*/
                        psNear.exit();
                    }
                }
                if( !succeededThrow ) {
                    if( centerPoint->tryCount >= maxAttemptCount ) {
                        // FF_LOG(debug) << "Removing particle " << availableDarts[generation][p].point << endl;

                        // remove the particle p from the available list
                        if( p != availableDarts[generation].size() - 1 )
                            std::swap( *centerPoint,
                                       availableDarts[generation][availableDarts[generation].size() - 1] );

                        availableDarts[generation].pop_back();

                    } else { // The particle was rejected, but not yet at its attempt limit
                        ++centerPoint->tryCount;
                    }
                }

                psThrow.exit();

                ++iterationCount;

                // if(  count >= theoreticalCount ) {
                //	break;
                // }
            }
        }

        FF_LOG( debug ) << "\tGen 1 Promotion Count: " << promotedCount << std::endl;

        // std::ofstream fout("c:\\temp\\darts_log.txt");
        // FF_LOG(debug) << "Iterations: " << iterationCount << "\n";

        // FF_LOG(debug) << "Remaining Test Points: " << youngAvailableDarts.size() << "\n";

        FF_LOG( debug ) << "\t" << psInsert << std::endl;
        FF_LOG( debug ) << "\t" << psBounds << std::endl;
        FF_LOG( debug ) << "\t" << psThrow << std::endl;
        FF_LOG( debug ) << "\t" << psNear << std::endl;
        FF_LOG( debug ) << "\t" << psNeighbor << std::endl;

        /* verbose outputting...
        float particleVolume = 4 * (float) frantic::math::pi_value * radius * radius * radius / 3.f;
        float totalVolume = (count * particleVolume ) ;
        float boxVolume = ( (m_boundBox.xsize() +  2*radius)*(m_boundBox.ysize() +  2*radius)*(m_boundBox.zsize() +
        2*radius) ); float actualPacking = totalVolume/(0.7405f * boxVolume);


        std::cout << "Particle Vol    : " << particleVolume<< "  Total Volume  : " << totalVolume << "  Box Volume : "
        << boxVolume << std::endl; std::cout << "Maximum Possible Packing Count: " <<
        estimate_particle_count(radius,1.f) << std::endl; std::cout << "Expected Packing^3: " <<
        relativePacking*relativePacking*relativePacking << " Actual Packing^3: " << actualPacking << std::endl;

        std::cout << "Count: " << count << std::endl;
        std::cout << "AvailablePoint Size: " << availablePoints.size() << std::endl;
        std::cout << "Internal Particle Count: " << darts->particle_count() << std::endl;

        std::cout << "Number of Iterations: " << iterationCount << std::endl;
        std::cout << "Max Iterations: " << maxIterations << std::endl;
        //*/

        // maxIterations = iterationCount; // return the number of iterations actually taken

        for( int generation = 0; generation < m_genCount; ++generation )
            availableDarts[generation].clear();

        // It hardly makes sense to error out here, because the halting condition is not the number of particles
        // provided by the generational limit
        /*if( count != theoreticalCount ) {
          throw std::runtime_error("Failed to seed particle to the requested packing level " +
        boost::lexical_cast<std::string>(relativePacking)
            + " with the maximum attempt count set at: " + boost::lexical_cast<std::string>(maxAttemptCount) );
        }//*/
    }
};

} // namespace graphics
} // namespace frantic

#pragma warning( pop )
