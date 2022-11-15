// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/particle_grid_tree.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/simd/float_v.hpp>

#include <tbb/parallel_for.h>

using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::simd;
using namespace frantic::volumetrics;

namespace frantic {
namespace particles {

//
//
// Iterator implementations
//
//

particle_grid_tree::iterator::iterator( pgt_node_type* currentNode, size_t particleSize ) {
    // creates an iterator object starting at the first particle, or an "end" iterator if currentNode is NULL.
    m_currentNode = currentNode;
    m_particleSize = particleSize;
    m_currentParticle = NULL;
    m_particleEndPtr = NULL;
    m_gridIndex = 0;
    if( m_currentNode ) {
        particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
        m_currentParticle = leafData->grid[m_gridIndex].begin();
        m_particleEndPtr = leafData->grid[m_gridIndex].end();
        // find first non-empty grid index
        while( m_currentParticle == NULL ) {
            ++m_gridIndex;
            m_currentParticle = leafData->grid[m_gridIndex].begin();
            m_particleEndPtr = leafData->grid[m_gridIndex].end();
        }
    }
}

particle_grid_tree::iterator& particle_grid_tree::iterator::operator++() {
    // check to see if we have already reached the end
    if( m_currentNode ) {
        // increment within this grid
        m_currentParticle += m_particleSize;
        // check to see if we've gone over our current grid's buffer
        if( m_currentParticle == m_particleEndPtr ) {

            // move to the next populated grid
            // the loop is here because sometimes grids don't have particles. in this case we move to the next grid,
            // etc.
            particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
            do {
                // move to the next grid
                ++m_gridIndex;

                // check to see if we've overrun our grid array for this leaf
                if( m_gridIndex >= PGT_SL * PGT_SL * PGT_SL ) {
                    // move to the next leaf node
                    m_gridIndex = 0;
                    m_currentNode = leafData->nextNode;
                    if( m_currentNode ) {
                        // we have another leaf node, move to it.
                        leafData = m_currentNode->leafData;
                    } else {
                        // we've reached the end of the leaf nodes, we're done all iteration
                        m_currentParticle = NULL;
                        m_particleEndPtr = NULL;
                        break;
                    }
                }

                // set our pointers to the first particle in this grid
                m_currentParticle = leafData->grid[m_gridIndex].begin();
                m_particleEndPtr = leafData->grid[m_gridIndex].end();

            } while( !m_currentParticle );
        }
    }
    return *this;
}

particle_grid_tree::const_iterator::const_iterator( const pgt_node_type* currentNode, size_t particleSize ) {
    // SAME AS THE NON-CONST VERSION
    // creates an iterator object starting at the first particle, or an "end" iterator if currentNode is NULL.
    m_currentNode = currentNode;
    m_particleSize = particleSize;
    m_currentParticle = NULL;
    m_particleEndPtr = NULL;
    m_gridIndex = 0;
    if( m_currentNode ) {
        particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
        m_currentParticle = leafData->grid[m_gridIndex].begin();
        m_particleEndPtr = leafData->grid[m_gridIndex].end();
        // find first non-empty grid index
        while( m_currentParticle == NULL ) {
            ++m_gridIndex;
            m_currentParticle = leafData->grid[m_gridIndex].begin();
            m_particleEndPtr = leafData->grid[m_gridIndex].end();
        }
    }
}

particle_grid_tree::const_iterator::const_iterator( const iterator& rhs ) {
    m_currentNode = rhs.m_currentNode;
    m_gridIndex = rhs.m_gridIndex;
    m_currentParticle = rhs.m_currentParticle;
    m_particleEndPtr = rhs.m_particleEndPtr;
    m_particleSize = rhs.m_particleSize;
}

particle_grid_tree::const_iterator& particle_grid_tree::const_iterator::operator++() {
    // SAME AS THE NON-CONST VERSION
    // check to see if we haven't already reached the end
    if( m_currentNode ) {
        // increment within this grid
        m_currentParticle += m_particleSize;
        // check to see if we've gone over our current grid's buffer
        if( m_currentParticle == m_particleEndPtr ) {
            // move to the next populated grid
            // the loop is here because sometimes grids don't have particles. in this case we move to the next grid,
            // etc.
            particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
            do {
                // move to the next grid
                ++m_gridIndex;

                // check to see if we've overrun our grid array for this leaf
                if( m_gridIndex >= PGT_SL * PGT_SL * PGT_SL ) {
                    // move to the next leaf node
                    m_gridIndex = 0;
                    m_currentNode = leafData->nextNode;
                    if( m_currentNode ) {
                        // we have another leaf node, move to it.
                        leafData = m_currentNode->leafData;
                    } else {
                        // we've reached the end of the leaf nodes, we're done all iteration
                        m_currentParticle = NULL;
                        m_particleEndPtr = NULL;
                        break;
                    }
                }

                // set our pointers to the first particle in this grid
                m_currentParticle = leafData->grid[m_gridIndex].begin();
                m_particleEndPtr = leafData->grid[m_gridIndex].end();

            } while( !m_currentParticle );
        }
    }
    return *this;
}

particle_grid_tree::node_iterator::node_iterator(
    pgt_node_type* currentNode,
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
    const frantic::volumetrics::voxel_coord_system& vcs, size_t particleSize ) {
    // creates an iterator object starting at the first grid, or an "end" iterator if currentNode is NULL.
    m_currentNode = currentNode;
    m_particleSize = particleSize;
    m_positionAccessor = positionAccessor;
    m_vcs = vcs;
    m_gridIndex = 0;
    if( m_currentNode ) {
        particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
        char* gridStart = leafData->grid[m_gridIndex].begin();
        // find first non-empty grid index
        while( gridStart == NULL ) {
            ++m_gridIndex;
            gridStart = leafData->grid[m_gridIndex].begin();
        }
    }
}

particle_grid_tree::node_iterator::particle_iterator_pair
particle_grid_tree::node_iterator::particle_iterators() const {

    if( m_currentNode )
        return particle_iterator_pair(
            particle_array_iterator( m_currentNode->leafData->grid[m_gridIndex].begin(), m_particleSize ),
            particle_array_iterator( m_currentNode->leafData->grid[m_gridIndex].end(), m_particleSize ) );
    else
        return particle_iterator_pair( particle_array_iterator( NULL, m_particleSize ),
                                       particle_array_iterator( NULL, m_particleSize ) );
}

frantic::graphics::boundbox3f particle_grid_tree::node_iterator::world_bounds() const {
    if( m_currentNode ) {

        // get our leaf node bounds in world coordinates
        vector3f worldCoord = m_vcs.get_world_coord( m_currentNode->coord );

        // determine grid x,y,z from our m_gridIndex
        size_t gridZ = m_gridIndex / ( PGT_SL * PGT_SL );
        size_t gridY = ( m_gridIndex / PGT_SL ) % PGT_SL;
        size_t gridX = m_gridIndex % PGT_SL;

        // include grid offset
        float gridLength = m_vcs.voxel_length() / PGT_SL;
        worldCoord.x += gridX * gridLength;
        worldCoord.y += gridY * gridLength;
        worldCoord.z += gridZ * gridLength;

        // create boundbox
        return boundbox3f( worldCoord, worldCoord + vector3f( gridLength ) );
    }

    // end iterators will return empty bounds
    return frantic::graphics::boundbox3f();
}

particle_grid_tree::node_iterator& particle_grid_tree::node_iterator::operator++() {
    // check to see if we have already reached the end
    if( m_currentNode ) {
        // move to the next populated grid
        // the loop is here because sometimes grids don't have particles. in this case we move to the next grid, etc.
        particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
        char* gridStart = NULL;
        do {
            // move to the next grid
            ++m_gridIndex;

            // check to see if we've overrun our grid array for this leaf
            if( m_gridIndex >= PGT_SL * PGT_SL * PGT_SL ) {
                // move to the next leaf node
                m_gridIndex = 0;
                m_currentNode = leafData->nextNode;
                if( m_currentNode ) {
                    // we have another leaf node, move to it.
                    leafData = m_currentNode->leafData;
                } else {
                    // we've reached the end of the leaf nodes, we're done all iteration
                    m_gridIndex = 0;
                    break;
                }
            }

            // set the grid start to see if we have any particles in this grid
            gridStart = leafData->grid[m_gridIndex].begin();

        } while( !gridStart );
    }
    return *this;
}

particle_grid_tree::const_node_iterator::const_node_iterator(
    const pgt_node_type* currentNode,
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
    const frantic::volumetrics::voxel_coord_system& vcs, size_t particleSize ) {
    // SAME AS NON-CONST VERSION
    // creates an iterator object starting at the first grid, or an "end" iterator if currentNode is NULL.
    m_currentNode = currentNode;
    m_particleSize = particleSize;
    m_positionAccessor = positionAccessor;
    m_vcs = vcs;
    m_gridIndex = 0;
    if( m_currentNode ) {
        particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
        char* gridStart = leafData->grid[m_gridIndex].begin();
        // find first non-empty grid index
        while( gridStart == NULL ) {
            ++m_gridIndex;
            gridStart = leafData->grid[m_gridIndex].begin();
        }
    }
}

particle_grid_tree::const_node_iterator::const_node_iterator( const node_iterator& rhs ) {
    m_currentNode = rhs.m_currentNode;
    m_gridIndex = rhs.m_gridIndex;
    m_particleSize = rhs.m_particleSize;
    m_positionAccessor = rhs.m_positionAccessor;
    m_vcs = rhs.m_vcs;
}

particle_grid_tree::const_node_iterator::const_particle_iterator_pair
particle_grid_tree::const_node_iterator::const_particle_iterators() const {
    if( m_currentNode )
        return const_particle_iterator_pair(
            const_particle_array_iterator( m_currentNode->leafData->grid[m_gridIndex].begin(), m_particleSize ),
            const_particle_array_iterator( m_currentNode->leafData->grid[m_gridIndex].end(), m_particleSize ) );
    else
        return const_particle_iterator_pair( const_particle_array_iterator( NULL, m_particleSize ),
                                             const_particle_array_iterator( NULL, m_particleSize ) );
}

frantic::graphics::boundbox3f particle_grid_tree::const_node_iterator::world_bounds() const {
    if( m_currentNode ) {

        // get our leaf node bounds in world coordinates
        vector3f worldCoord = m_vcs.get_world_coord( m_currentNode->coord );

        // determine grid x,y,z from our m_gridIndex
        size_t gridZ = m_gridIndex / ( PGT_SL * PGT_SL );
        size_t gridY = ( m_gridIndex / PGT_SL ) % PGT_SL;
        size_t gridX = m_gridIndex % PGT_SL;

        // include grid offset
        float gridLength = m_vcs.voxel_length() / PGT_SL;
        worldCoord.x += gridX * gridLength;
        worldCoord.y += gridY * gridLength;
        worldCoord.z += gridZ * gridLength;

        // create boundbox
        return boundbox3f( worldCoord, worldCoord + vector3f( gridLength ) );
    }

    // end iterators will return empty bounds
    return frantic::graphics::boundbox3f();
}

particle_grid_tree::const_node_iterator& particle_grid_tree::const_node_iterator::operator++() {
    // SAME AS NON-CONST VERSION
    // check to see if we have already reached the end
    if( m_currentNode ) {
        // move to the next populated grid
        // the loop is here because sometimes grids don't have particles. in this case we move to the next grid, etc.
        particle_grid_tree_leaf_data* leafData = m_currentNode->leafData;
        char* gridStart = NULL;
        do {
            // move to the next grid
            ++m_gridIndex;

            // check to see if we've overrun our grid array for this leaf
            if( m_gridIndex >= PGT_SL * PGT_SL * PGT_SL ) {
                // move to the next leaf node
                m_gridIndex = 0;
                m_currentNode = leafData->nextNode;
                if( m_currentNode ) {
                    // we have another leaf node, move to it.
                    leafData = m_currentNode->leafData;
                } else {
                    // we've reached the end of the leaf nodes, we're done all iteration
                    m_gridIndex = 0;
                    break;
                }
            }

            // set the grid start to see if we have any particles in this grid
            gridStart = leafData->grid[m_gridIndex].begin();

        } while( !gridStart );
    }
    return *this;
}

//
//
// Tree construction
//
//

particle_grid_tree::particle_grid_tree() { reset( channel_map(), voxel_coord_system() ); }

particle_grid_tree::particle_grid_tree( const frantic::channels::channel_map& pcm ) {
    reset( pcm, voxel_coord_system() );
}

particle_grid_tree::particle_grid_tree( const frantic::channels::channel_map& pcm, float voxelLength ) {
    reset( pcm, voxel_coord_system( vector3f(), voxelLength ) );
}

particle_grid_tree::particle_grid_tree( const frantic::channels::channel_map& pcm,
                                        const frantic::volumetrics::voxel_coord_system& vcs ) {
    reset( pcm, vcs );
}

particle_grid_tree::~particle_grid_tree() {
    // the parent class's destructor will take care of freeing memory
}

void particle_grid_tree::clear() { reset( channel_map(), voxel_coord_system() ); }

void particle_grid_tree::reset( const frantic::channels::channel_map& pcm, float voxelLength ) {
    reset( pcm, voxel_coord_system( vector3f(), voxelLength ) );
}

void particle_grid_tree::reset( const frantic::channels::channel_map& pcm,
                                const frantic::volumetrics::voxel_coord_system& vcs ) {
    // called from all the constructors, clears, and resets
    m_channelMap = pcm;

    // set the "Position" channel accessor
    if( m_channelMap.channel_definition_complete() ) {
        if( pcm.has_channel( _T("Position") ) )
            m_positionAccessor = pcm.get_accessor<vector3f>( _T("Position") );
        else
            throw std::runtime_error(
                "particle_grid_tree.reset: Provided channel_map object does not have a \"Position\" channel." );
    } else {
        // default constructor allows this to be defined later
        m_positionAccessor = frantic::channels::channel_accessor<frantic::graphics::vector3f>();
    }
    m_leafNodeLinkedListStart = NULL;
    m_leafNodeListArray.clear();
    m_leafNodeListArray.resize( 8 ); // ensure the size of the array is always 8
    m_vcs = voxel_coord_system( vcs.world_origin(),
                                vcs.voxel_length() *
                                    PGT_SL ); // voxels store a PGT_SL^3 grid of particles, so increase the voxel length
    m_gcs = voxel_coord_system( vcs.world_origin(), vcs.voxel_length() );
    m_particleCount = 0;
    grid_tree_base<particle_grid_tree_leaf_data, PGT_SL>::clear();
}

void particle_grid_tree::swap( particle_grid_tree& rhs ) {
    // swap all our member variables
    m_channelMap.swap( rhs.m_channelMap );
    std::swap( m_positionAccessor, rhs.m_positionAccessor );
    m_vcs.swap( rhs.m_vcs );
    m_gcs.swap( rhs.m_gcs );
    std::swap( m_particleCount, rhs.m_particleCount );
    std::swap( m_leafNodeLinkedListStart, rhs.m_leafNodeLinkedListStart );
    m_leafNodeListArray.swap( rhs.m_leafNodeListArray );
    grid_tree_base<particle_grid_tree_leaf_data, PGT_SL>::swap(
        static_cast<grid_tree_base<particle_grid_tree_leaf_data, PGT_SL>&>( rhs ) ); // call our superclass's swap
}

frantic::graphics::boundbox3f particle_grid_tree::compute_particle_bounds() const {
    // this call is expensive and should not be called in performance situations
    // also, this function could be optimized

    vector3f minimum( std::numeric_limits<float>::max() );
    vector3f maximum( -std::numeric_limits<float>::max() );

    // iterate through the particles and get the boundbox
    const_iterator iter = begin();
    const_iterator iterEnd = end();
    for( ; iter != iterEnd; ++iter ) {
        const vector3f& pos = m_positionAccessor( *iter );
        minimum.x = std::min( minimum.x, pos.x );
        minimum.y = std::min( minimum.y, pos.y );
        minimum.z = std::min( minimum.z, pos.z );
        maximum.x = std::max( maximum.x, pos.x );
        maximum.y = std::max( maximum.y, pos.y );
        maximum.z = std::max( maximum.z, pos.z );
    }
    return boundbox3f( minimum, maximum );
}

//
//
// Particle access
//
//

bool particle_grid_tree::has_particle_near( const frantic::graphics::vector3f& worldLocation,
                                            float searchRadius ) const {
    float searchRadiusSquared = searchRadius * searchRadius;
    size_t particleSize = m_channelMap.structure_size();

    std::vector<frantic::graphics::raw_byte_buffer*> byteBuffers;
    get_grid_neighborhood( worldLocation, searchRadius, byteBuffers );
    for( size_t i = 0; i < byteBuffers.size(); ++i ) {
        char* iter = byteBuffers[i]->begin();
        char* iterEnd = byteBuffers[i]->end();
        for( ; iter != iterEnd; iter = iter + particleSize ) {
            if( vector3f::distance_squared( m_positionAccessor( iter ), worldLocation ) < searchRadiusSquared )
                return true;
        }
    }
    return false;
}

bool particle_grid_tree::get_closest_particle( const frantic::graphics::vector3f& worldLocation, float searchRadius,
                                               char* outParticle, float& outDistance ) const {
    float searchRadiusSquared = searchRadius * searchRadius;
    size_t particleSize = m_channelMap.structure_size();

    bool foundParticle = false;
    outDistance = std::numeric_limits<float>::max();

    std::vector<frantic::graphics::raw_byte_buffer*> byteBuffers;
    get_grid_neighborhood( worldLocation, searchRadius, byteBuffers );
    for( size_t i = 0; i < byteBuffers.size(); ++i ) {
        char* iter = byteBuffers[i]->begin();
        char* iterEnd = byteBuffers[i]->end();
        for( ; iter != iterEnd; iter = iter + particleSize ) {
            float distanceSquared = vector3f::distance_squared( m_positionAccessor( iter ), worldLocation );
            if( distanceSquared < searchRadiusSquared ) {
                if( distanceSquared < outDistance ) {
                    outDistance = distanceSquared;
                    memcpy( outParticle, iter, particleSize );
                    foundParticle = true;
                }
            }
        }
    }
    if( foundParticle )
        outDistance = sqrtf( outDistance );
    return foundParticle;
}

namespace {
struct get_particles_in_point_range_body {
    const size_t m_particleSize;
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& m_positionAccessor;
    const frantic::graphics::vector3f& m_worldLocation;
    const float m_searchRadiusSquared;
    std::vector<char*>& m_outParticles;

#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4822 ) // local class member function does not have a body
#endif
    get_particles_in_point_range_body& operator=( get_particles_in_point_range_body& ); // not implemented
#if defined( _MSC_VER )
#pragma warning( pop )
#endif

    get_particles_in_point_range_body(
        const std::size_t particleSize,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::graphics::vector3f& worldLocation, const float searchRadius, std::vector<char*>& outParticles )
        : m_particleSize( particleSize )
        , m_positionAccessor( positionAccessor )
        , m_worldLocation( worldLocation )
        , m_searchRadiusSquared( searchRadius * searchRadius )
        , m_outParticles( outParticles ) {}

    bool operator()( frantic::graphics::raw_byte_buffer* byteBuffer ) {
        vector3f nearestInBox;
        char* iter = byteBuffer->begin();
        char* iterEnd = byteBuffer->end();
        for( ; iter != iterEnd; iter = iter + m_particleSize ) {
            if( vector3f::distance_squared( m_positionAccessor( iter ), m_worldLocation ) < m_searchRadiusSquared )
                m_outParticles.push_back( iter );
        }
        return true;
    }
};
} // namespace

void particle_grid_tree::get_particles_in_range( const frantic::graphics::vector3f& worldLocation, float searchRadius,
                                                 std::vector<char*>& outParticles ) const {
    for_each_grid_in_neighborhood( worldLocation, searchRadius,
                                   get_particles_in_point_range_body( m_channelMap.structure_size(), m_positionAccessor,
                                                                      worldLocation, searchRadius, outParticles ) );
}

namespace {

#ifdef FRANTIC_HAS_SSE2
// Minimum features of boundbox3f required for clamping,
// implemented using SSE.
class sse_clamp_box {
  public:
    typedef float_v vector3f_type;

    sse_clamp_box( const boundbox3f& bounds )
        : m_min( bounds.minimum() )
        , m_max( bounds.maximum() ) {}

    float_v clamp_nothrow( const float_v& x ) const { return std::min( std::max( x, m_min ), m_max ); }

    float xminimum() const { return m_min[0]; }
    float yminimum() const { return m_min[1]; }
    float zminimum() const { return m_min[2]; }

    float xmaximum() const { return m_max[0]; }
    float ymaximum() const { return m_max[1]; }
    float zmaximum() const { return m_max[2]; }

  private:
    float_v m_min;
    float_v m_max;
};

// Workalike for the distance_squared( vector3f, vector3f ) function.
float distance_squared( const float_v& a, const float_v& b ) {
    float_v x2 = frantic::math::square( a - b );
    return x2.sum3();
}

typedef sse_clamp_box clamp_box_t;
#else
typedef boundbox3f clamp_box_t;
#endif

class get_particles_in_box_range_body {
  public:
    get_particles_in_box_range_body( const std::size_t particleSize, const channel_accessor<vector3f>& positionAccessor,
                                     const boundbox3f& worldBox, const float searchRadius,
                                     std::vector<char*>& outParticles )
        : m_particleSize( particleSize )
        , m_particleSizeX4( particleSize * 4 )
        , m_positionAccessor( positionAccessor )
        , m_searchRadiusSquared( searchRadius * searchRadius )
        , m_worldBox( worldBox )
        , m_outParticles( outParticles ) {}

    bool operator()( raw_byte_buffer* byteBuffer ) const {
        typedef clamp_box_t::vector3f_type vector3f_t;

        using frantic::math::square;
        char* iter = byteBuffer->begin();
        char* iterEnd = byteBuffer->end();

        const std::size_t particleSize = m_particleSize;

#ifdef FRANTIC_HAS_SSE2
        const std::size_t particleSizeX4 = m_particleSizeX4;

        if( static_cast<std::size_t>( iterEnd - iter ) > particleSizeX4 ) {
            BOOST_STATIC_ASSERT( float_v::static_size == 4 );

            const char* iterX4End = iter + ( iterEnd - iter - 1 ) / particleSizeX4 * particleSizeX4;

            const float_v xMin( m_worldBox.xminimum() );
            const float_v yMin( m_worldBox.yminimum() );
            const float_v zMin( m_worldBox.zminimum() );

            const float_v xMax( m_worldBox.xmaximum() );
            const float_v yMax( m_worldBox.ymaximum() );
            const float_v zMax( m_worldBox.zmaximum() );

            for( ; iter != iterX4End; iter += particleSizeX4 ) {
                // Loading 4 floats for the Position, while only the first 3 are valid,
                // because it's faster.  This could go out of bounds on the last
                // particle, so we set the iteration bounds to always process the last
                // particle using the single-particle code path.
                float_v x = float_v::load4( &m_positionAccessor( iter ).x );
                float_v y = float_v::load4( &m_positionAccessor( iter + particleSize ).x );
                float_v z = float_v::load4( &m_positionAccessor( iter + 2 * particleSize ).x );
                float_v w = float_v::load4( &m_positionAccessor( iter + 3 * particleSize ).x );

                // Transpose the particle data so that x holds the x coordinate for all
                // four particles, etc.
                _MM_TRANSPOSE4_PS( x.native(), y.native(), z.native(), w.native() );

                const float_v xClamped = std::min( std::max( x, xMin ), xMax );
                const float_v yClamped = std::min( std::max( y, yMin ), yMax );
                const float_v zClamped = std::min( std::max( z, zMin ), zMax );

                const float_v dx = x - xClamped;
                const float_v dy = y - yClamped;
                const float_v dz = z - zClamped;

                const float_v distanceSquared = square( dx ) + square( dy ) + square( dz );

                const float_v inRange = distanceSquared < m_searchRadiusSquared;

                const int inRangeMask = _mm_movemask_ps( inRange.native() );
                if( inRangeMask & 1 ) {
                    m_outParticles.push_back( iter );
                }
                if( inRangeMask & 2 ) {
                    m_outParticles.push_back( iter + particleSize );
                }
                if( inRangeMask & 4 ) {
                    m_outParticles.push_back( iter + 2 * particleSize );
                }
                if( inRangeMask & 8 ) {
                    m_outParticles.push_back( iter + 3 * particleSize );
                }
            }
        }
#endif

        for( ; iter != iterEnd; iter += particleSize ) {
            const vector3f_t particlePosition( m_positionAccessor( iter ) );
            const vector3f_t nearestInBox = m_worldBox.clamp_nothrow( particlePosition );
            if( distance_squared( particlePosition, nearestInBox ) < m_searchRadiusSquared ) {
                m_outParticles.push_back( iter );
            }
        }
        return true;
    }

  private:
    get_particles_in_box_range_body& operator=( const get_particles_in_box_range_body& ); // not implemented

    const size_t m_particleSize;
    const size_t m_particleSizeX4;
    const channel_accessor<vector3f> m_positionAccessor;
    const float m_searchRadiusSquared;
    const clamp_box_t m_worldBox;
    std::vector<char*>& m_outParticles;
};

} // anonymous namespace

void particle_grid_tree::get_particles_in_range( const frantic::graphics::boundbox3f& worldBox, float searchRadius,
                                                 std::vector<char*>& outParticles ) const {
    boundbox3f worldBoxExpanded( worldBox );
    worldBoxExpanded.expand( searchRadius );

    if( worldBox.is_empty() ) {
        throw std::runtime_error( "particle_grid_tree.get_particles_in_range: The world box is empty" );
    }

    for_each_grid_in_neighborhood_box(
        worldBoxExpanded, get_particles_in_box_range_body( m_channelMap.structure_size(), m_positionAccessor, worldBox,
                                                           searchRadius, outParticles ) );
}

void particle_grid_tree::particle_particle_interactions( particle_particle_interaction_function_t interactionFunction,
                                                         void* userData, float searchRadius ) {
    if( searchRadius > m_gcs.voxel_length() )
        throw std::runtime_error( "particle_grid_tree.particle_particle_interactions: Particle-particle interactions "
                                  "require that the search radius " +
                                  boost::lexical_cast<std::string>( searchRadius ) +
                                  " is at most the length of one voxel of the particle_grid_tree: " +
                                  boost::lexical_cast<std::string>( m_gcs.voxel_length() ) );

    // create the task schedular object. this is tied to some sort of global tbb object that is reference counted.
    // all we need to do is create this and hope for the best. this is already done in python library main, but this can
    // be used from anywhere.
    // tbb::task_scheduler_init taskScheduleInit; //please move out.

    // for all 8 non-interacting node sets (see wiki page)
    for( int i = 0; i < 8; ++i ) {

        size_t numLeaves = m_leafNodeListArray[i].size();
        if( numLeaves > 0 ) {

            // create a tbb body object that has all the data needed for particle-particle interactions
            detail::particle_grid_tree_particle_particle_interactions_tbb_body tbbObject(
                &m_leafNodeListArray[i], interactionFunction, userData, searchRadius, m_positionAccessor,
                m_channelMap.structure_size() );

            // 1/16th of the total number of leaf nodes seems a good starting point for ttb's grain size
            size_t tbbGrainSize = ( numLeaves < 16 ) ? 1 : ( numLeaves / 16 );

// call parallel_for on this object
#ifndef FRANTIC_DISABLE_THREADS
            tbb::parallel_for( tbb::blocked_range<size_t>( 0, m_leafNodeListArray[i].size(), tbbGrainSize ),
                               tbbObject );
#else
#pragma message( "Threads are disabled" )
            tbbObject( tbb::blocked_range<size_t>( 0, m_leafNodeListArray[i].size(), tbbGrainSize ) );
#endif

            // NOTE: to do the non-parallel version, call this function instead:
            // tbbObject( tbb::blocked_range<size_t>(0,m_leafNodeListArray[i].size()) );
        }
    }
}

void particle_grid_tree::process_particles_in_bounds( void* userData, const frantic::graphics::boundbox3f& worldBox,
                                                      particle_evaluation_function_t evaluationFunction ) {
    size_t particleSize = m_channelMap.structure_size();

    std::vector<frantic::graphics::raw_byte_buffer*> byteBuffers;
    get_grid_neighborhood_box( worldBox, byteBuffers );
    for( size_t i = 0; i < byteBuffers.size(); ++i ) {
        char* iter = byteBuffers[i]->begin();
        char* iterEnd = byteBuffers[i]->end();
        for( ; iter != iterEnd; iter = iter + particleSize ) {
            if( worldBox.contains( m_positionAccessor( iter ) ) )
                evaluationFunction( userData, iter );
        }
    }
}

namespace {
struct process_particles_in_bounds_body {
    const std::vector<frantic::graphics::raw_byte_buffer*>& m_byteBuffers;
    const size_t m_particleSize;
    const frantic::graphics::boundbox3f& m_worldBox;
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& m_positionAccessor;
    void* m_userData;
    particle_grid_tree::particle_evaluation_function_t m_evaluationFunction;

#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4822 ) // local class member function does not have a body
#endif
    process_particles_in_bounds_body& operator=( process_particles_in_bounds_body& ); // not implemented
#if defined( _MSC_VER )
#pragma warning( pop )
#endif

    process_particles_in_bounds_body(
        const std::vector<frantic::graphics::raw_byte_buffer*>& byteBuffers, const size_t particleSize,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::graphics::boundbox3f& worldBox, void* userData,
        particle_grid_tree::particle_evaluation_function_t evaluationFunction )
        : m_byteBuffers( byteBuffers )
        , m_particleSize( particleSize )
        , m_worldBox( worldBox )
        , m_positionAccessor( positionAccessor )
        , m_userData( userData )
        , m_evaluationFunction( evaluationFunction ) {}
    bool operator()( const tbb::blocked_range<size_t>& r ) const {
        for( size_t i = r.begin(); i < r.end(); ++i ) {
            char* iter = m_byteBuffers[i]->begin();
            char* iterEnd = m_byteBuffers[i]->end();
            for( ; iter != iterEnd; iter = iter + m_particleSize ) {
                m_evaluationFunction( m_userData, iter );
            }
        }
        return true;
    }
};
} // namespace

void particle_grid_tree::process_particles_in_bounds_mt( void* userData, const frantic::graphics::boundbox3f& worldBox,
                                                         particle_evaluation_function_t evaluationFunction ) {
    const size_t particleSize = m_channelMap.structure_size();

    std::vector<frantic::graphics::raw_byte_buffer*> byteBuffers;
    get_grid_neighborhood_box( worldBox, byteBuffers );

    tbb::parallel_for( tbb::blocked_range<size_t>( 0, byteBuffers.size() ),
                       process_particles_in_bounds_body( byteBuffers, particleSize, m_positionAccessor, worldBox,
                                                         userData, evaluationFunction ),
                       tbb::auto_partitioner() );
}

namespace {
struct process_particles_in_range_body {
    const size_t m_particleSize;
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& m_positionAccessor;
    const frantic::graphics::vector3f& m_worldLocation;
    const float m_searchRadiusSquared;
    void* m_userData;
    particle_grid_tree::particle_evaluation_function_t m_evaluationFunction;

#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4822 ) // local class member function does not have a body
#endif
    process_particles_in_range_body& operator=( process_particles_in_range_body& ); // not implemented
#if defined( _MSC_VER )
#pragma warning( pop )
#endif

    process_particles_in_range_body(
        const std::size_t particleSize,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::graphics::vector3f& worldLocation, const float searchRadius, void* userData,
        particle_grid_tree::particle_evaluation_function_t evaluationFunction )
        : m_particleSize( particleSize )
        , m_positionAccessor( positionAccessor )
        , m_worldLocation( worldLocation )
        , m_searchRadiusSquared( searchRadius * searchRadius )
        , m_userData( userData )
        , m_evaluationFunction( evaluationFunction ) {}

    bool operator()( frantic::graphics::raw_byte_buffer* byteBuffer ) {
        char* iter = byteBuffer->begin();
        char* iterEnd = byteBuffer->end();
        for( ; iter != iterEnd; iter += m_particleSize ) {
            if( vector3f::distance_squared( m_positionAccessor( iter ), m_worldLocation ) < m_searchRadiusSquared ) {
                m_evaluationFunction( m_userData, iter );
            }
        }
        return true;
    }
};
} // namespace

void particle_grid_tree::process_particles_in_range( void* userData, const frantic::graphics::vector3f& worldPosition,
                                                     float searchRadius,
                                                     particle_evaluation_function_t evaluationFunction ) {
    const size_t particleSize = m_channelMap.structure_size();

    for_each_grid_in_neighborhood( worldPosition, searchRadius,
                                   process_particles_in_range_body( particleSize, m_positionAccessor, worldPosition,
                                                                    searchRadius, userData, evaluationFunction ) );
}

void particle_grid_tree::particle_grid_interactions( const frantic::graphics::boundbox3& /*voxelGrid*/,
                                                     const frantic::volumetrics::voxel_coord_system& /*gridVCS*/,
                                                     std::size_t /*gridParticleSize*/, void* /*gridParticleArray*/,
                                                     particle_grid_interaction_function_t /*interactionFunction*/,
                                                     void* /*userData*/, float /*searchRadius*/ ) {
    // THIS FUNCTION IS TO BE REMOVED ONCE THE UNUSED PARTICLE IMPLICIT POLICY THAT USED IT IS CLEANED OUT!
    throw std::runtime_error( "particle_grid_interactions isn't implemented." );
}

//
//
// Tree insertion functions
//
//

void particle_grid_tree::insert( const char* rawParticleData ) {
    const vector3f& position = m_positionAccessor( rawParticleData );
    if( !position.is_finite() )
        throw std::runtime_error(
            "particle_grid_tree::insert() - Attempted to insert a particle with non-finite Position channel value (" +
            position.str() + ")." );

    vector3 voxelCoord = vector3::from_floor( m_vcs.get_voxel_coord( position ) );

    // go to the leaf node. if this is a newly created leaf node, create leaf data for it
    pgt_node_type* leafNode = navigate_to_leaf_with_create( voxelCoord );
    if( !leafNode->leafData )
        add_new_leaf_data( voxelCoord, leafNode );

    // determine the correct grid index in our leaf node
    vector3 gridCoord = vector3::from_floor( m_gcs.get_voxel_coord( position ) );
    vector3 gridStart = voxelCoord * PGT_SL;

    // because of floating point errors for particles that sit right on grid borders, we need to do the clamp
    int gridIndexX = frantic::math::clamp<int>( gridCoord.x - gridStart.x, 0, PGT_SL - 1 );
    int gridIndexY = frantic::math::clamp<int>( gridCoord.y - gridStart.y, 0, PGT_SL - 1 );
    int gridIndexZ = frantic::math::clamp<int>( gridCoord.z - gridStart.z, 0, PGT_SL - 1 );

    // add this particle to the correct grid in our leaf node
    int gridIndex = gridIndexX + ( gridIndexY + gridIndexZ * PGT_SL ) * PGT_SL;
    leafNode->leafData->grid[gridIndex].add_element( rawParticleData, m_channelMap.structure_size() );

    // increment particle count
    ++m_particleCount;
}

bool particle_grid_tree::insert_with_proximity_constraint( const char* rawParticleData, float searchRadius ) {
    if( !has_particle_near( m_positionAccessor( rawParticleData ), searchRadius ) ) {
        insert( rawParticleData );
        return true;
    }
    return false;
}

bool particle_grid_tree::insert_with_proximity_constraint( const char* rawParticleData, float searchRadius,
                                                           const frantic::graphics::boundbox3f& cyclicBounds ) {
    // make a copy of the particle because we may need to modify it.
    // an alternative is to change in place which would provide the mapped particle.
    std::vector<char> tempParticle( m_channelMap.structure_size() );
    m_channelMap.copy_structure( &tempParticle[0], rawParticleData );

    const vector3f& worldPosition = m_positionAccessor.get( rawParticleData );

    // map the point into the bound box and collect any permutations of the point that intersect the boundary
    std::vector<vector3f> testPoints;
    cyclicBounds.compute_mod_boundary_points( worldPosition, searchRadius, testPoints );

    size_t testCount = testPoints.size();
    bool insertFailed = false;

    // conduct the insertion test of the mapped points using the usual functions
    for( size_t i = 0; i < testCount && !insertFailed; ++i ) {
        // test the boundary permutations first
        if( i + 1 < testCount ) {
            insertFailed = has_particle_near( testPoints[i], searchRadius );
        } else {
            // once all the boundary permutations are computed we can try and do the final insert
            m_positionAccessor.get( &tempParticle[0] ) = testPoints[i];
            // the insert function returns true if it's successful
            insertFailed = !insert_with_proximity_constraint( &tempParticle[0], searchRadius );
        }
    }

    return !insertFailed;
}

void particle_grid_tree::insert_particles( const particle_grid_tree& theParticles ) {
    particle_grid_tree::const_iterator iter = theParticles.begin(), iterEnd = theParticles.end();
    const channel_map& source = theParticles.get_channel_map();

    if( source != m_channelMap ) {
        // std::cout << "source size:" << source.structure_size() << std::endl;
        // std::cout << "dest size: " << m_channelMap.structure_size() << std::endl;

        channel_map_adaptor cma( m_channelMap, source );

        std::vector<char> temp( m_channelMap.structure_size() );
        for( ; iter != iterEnd; ++iter ) {
            cma.copy_structure( &temp[0], *iter );

            insert( &temp[0] );
            // reset things to zero
            memset( &temp[0], 0, m_channelMap.structure_size() );
        }

    } else {
        for( ; iter != iterEnd; ++iter )
            insert( *iter );
    }
}

void particle_grid_tree::insert_particles_with_proximity_constraint( const particle_grid_tree& theParticles,
                                                                     float searchRadius ) {
    particle_grid_tree::const_iterator iter = theParticles.begin(), iterEnd = theParticles.end();
    const channel_map& source = theParticles.get_channel_map();

    if( source != m_channelMap ) {
        channels::channel_map_adaptor cma( m_channelMap, source );

        std::vector<char> temp( m_channelMap.structure_size() );
        for( ; iter != iterEnd; ++iter ) {
            cma.copy_structure( &temp[0], *iter );
            insert_with_proximity_constraint( &temp[0], searchRadius );
            // reset things to zero
            memset( &temp[0], 0, m_channelMap.structure_size() );
        }
    } else {
        for( ; iter != iterEnd; ++iter ) {
            insert_with_proximity_constraint( *iter, searchRadius );
        }
    }
}

void particle_grid_tree::insert_particles( boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                           frantic::logging::progress_logger& progress ) {
    // set the channel_map on the stream
    pin->set_channel_map( m_channelMap );

    // can use a single particleData array as data storage with this particle channel map
    std::vector<char> rawParticle( m_channelMap.structure_size() );

    boost::int64_t loadedParticles = 0;
    boost::int64_t totalParticlesToLoad = pin->particle_count_left();

    while( pin->get_particle( &rawParticle[0] ) ) {
        insert( &rawParticle[0] );

        // If totalParticlesToLoad is less than zero, there's an unknown number to load and the progress is meaningless.
        if( totalParticlesToLoad > 0 ) {
            if( ( loadedParticles & 0xfff ) == 0 ) {
                progress.update_progress( loadedParticles, totalParticlesToLoad );
            }
            ++loadedParticles;
        }
    }
}

void particle_grid_tree::write_particles(
    boost::shared_ptr<frantic::particles::streams::particle_ostream> pout ) const {
    pout->set_channel_map( m_channelMap );
    const_iterator iter = begin();
    const_iterator iterEnd = end();
    for( ; iter != iterEnd; ++iter )
        pout->put_particle( *iter );
}

void particle_grid_tree::add_new_leaf_data( const frantic::graphics::vector3& voxelCoordinate,
                                            pgt_node_type* leafNode ) {
    // create a new leaf data object for this leaf
    particle_grid_tree_leaf_data* leafData = new particle_grid_tree_leaf_data;

    // hook this leaf node into our leaf-node linked list (for fast iteration)
    leafData->nextNode = m_leafNodeLinkedListStart;
    m_leafNodeLinkedListStart = leafNode;

    // add this to our leaf node list (for iteration purposes)
    // this is added to one of 8 lists. they are separated into buckets that cannot interact with each other
    // during particle-to-particle interactions. this is done for fast multi-threading of the p-to-p algorithm.
    int arrayIndex = ( ( leafNode->coord.x % 2 ) ? 0 : 1 ) + ( ( leafNode->coord.y % 2 ) ? 0 : 2 ) +
                     ( ( leafNode->coord.z % 2 ) ? 0 : 4 );
    m_leafNodeListArray[arrayIndex].push_back( leafNode );

    // link us to our neighbours, and link our neighbours to us
    for( int z = 0; z < 3; ++z ) {
        for( int y = 0; y < 3; ++y ) {
            for( int x = 0; x < 3; ++x ) {
                vector3 neighborCoord( voxelCoordinate.x + x - 1, voxelCoordinate.y + y - 1,
                                       voxelCoordinate.z + z - 1 );
                if( x != 1 || y != 1 || z != 1 ) {
                    pgt_node_type* neighbor = navigate_to_leaf( neighborCoord );
                    leafData->neighborhood[x + ( y + z * 3 ) * 3] = neighbor;
                    if( neighbor )
                        neighbor->leafData->neighborhood[( -x + 2 ) + ( ( -y + 2 ) + ( -z + 2 ) * 3 ) * 3] = leafNode;
                } else {
                    leafData->neighborhood[x + ( y + z * 3 ) * 3] = leafNode; // us (stored for convenience)
                }
            }
        }
    }

    // attach this leaf data to our node
    leafNode->leafData = leafData;
}

//
//
// private particle grid retrieval functions
//
//

/*
namespace detail {

struct append_to_byte_buffers {
  std::vector<frantic::graphics::raw_byte_buffer*>& m_byteBuffers;

#	pragma warning( push )
#	pragma warning( disable : 4822 ) // local class member function does not have a body
  append_to_byte_buffers& operator=( append_to_byte_buffers& ); // not implemented
#	pragma warning( pop )

  append_to_byte_buffers( std::vector<frantic::graphics::raw_byte_buffer*>& byteBuffers )
    : m_byteBuffers( byteBuffers )
  {
  }

  bool operator()( frantic::graphics::raw_byte_buffer* buffer ) {
    m_byteBuffers.push_back( buffer );
    return true;
  }
};

} // namespace detail
*/

void particle_grid_tree::get_grid_neighborhood(
    const frantic::graphics::vector3f& worldCoord, float radius,
    std::vector<frantic::graphics::raw_byte_buffer*>& outByteBuffers ) const {
    // bounds for this lookup
    boundbox3f worldBounds( worldCoord - vector3f( radius ), worldCoord + vector3f( radius ) );

    if( radius > m_vcs.voxel_length() ) {
        // if our lookup radius is larger than the voxel length, we can't use the simple
        // neighbour pointer lookup that our current function provides. instead we must call this function.
        get_grid_neighborhood_box( worldBounds, outByteBuffers );
        return;
    }

    // get center voxel (the one from which to do neighbour lookups)
    vector3 voxelCoord = vector3::from_floor( m_vcs.get_voxel_coord( worldCoord ) );
    pgt_node_type* leafNode = navigate_to_leaf( voxelCoord );

    if( !leafNode ) {
        // if we don't have a leaf directly in the lookup location, we still need to check the
        // neighbours. except we can't use the handy neighbour pointers, so we call this funciton.
        get_grid_neighborhood_box( worldBounds, outByteBuffers );
        return;
    }

    // get the start and end voxels, and the start and end grids
    vector3f voxelCoordMinFloat = m_vcs.get_voxel_coord( worldBounds.minimum() );
    vector3f voxelCoordMaxFloat = m_vcs.get_voxel_coord( worldBounds.maximum() );
    vector3 voxelCoordMin = vector3::from_floor( voxelCoordMinFloat );
    vector3 voxelCoordMax = vector3::from_ceil( voxelCoordMaxFloat );
    vector3 gridCoordMin = vector3::from_floor( voxelCoordMinFloat * PGT_SL );
    vector3 gridCoordMax = vector3::from_ceil( voxelCoordMaxFloat * PGT_SL );

    // iteration variables
    vector3 currVoxelCoord;
    vector3 currGridCoord;

    // iterate through all our voxels
    for( currVoxelCoord.z = voxelCoordMin.z; currVoxelCoord.z < voxelCoordMax.z; ++currVoxelCoord.z ) {
        currGridCoord.z = currVoxelCoord.z * PGT_SL;
        for( currVoxelCoord.y = voxelCoordMin.y; currVoxelCoord.y < voxelCoordMax.y; ++currVoxelCoord.y ) {
            currGridCoord.y = currVoxelCoord.y * PGT_SL;
            for( currVoxelCoord.x = voxelCoordMin.x; currVoxelCoord.x < voxelCoordMax.x; ++currVoxelCoord.x ) {
                currGridCoord.x = currVoxelCoord.x * PGT_SL;

                // get the appropriate neighbour from our base voxel for our current voxel
                int neighborIndex =
                    ( currVoxelCoord.x - voxelCoord.x + 1 ) +
                    ( ( currVoxelCoord.y - voxelCoord.y + 1 ) + ( currVoxelCoord.z - voxelCoord.z + 1 ) * 3 ) * 3;
                pgt_node_type* node = leafNode->leafData->neighborhood[neighborIndex];
                if( node ) {

                    // compute which grid to start and end on for this leaf node
                    vector3 start( gridCoordMin - currGridCoord );
                    vector3 end( gridCoordMax - currGridCoord );
                    vector3 startGrid( std::max<boost::int32_t>( start.x, 0 ), std::max<boost::int32_t>( start.y, 0 ),
                                       std::max<boost::int32_t>( start.z, 0 ) );
                    vector3 endGrid( std::min<boost::int32_t>( end.x, PGT_SL ),
                                     std::min<boost::int32_t>( end.y, PGT_SL ),
                                     std::min<boost::int32_t>( end.z, PGT_SL ) );

                    // go through each of the byte buffers in the grids and add them to our output vector
                    for( int gz = startGrid.z; gz < endGrid.z; ++gz ) {
                        for( int gy = startGrid.y; gy < endGrid.y; ++gy ) {
                            for( int gx = startGrid.x; gx < endGrid.x; ++gx ) {
                                int gridIndex = gx + ( gy + gz * PGT_SL ) * PGT_SL;
                                raw_byte_buffer* buffer = &node->leafData->grid[gridIndex];
                                if( buffer->ptr_at( 0 ) ) // only add non-empty buffers
                                    outByteBuffers.push_back( buffer );
                            }
                        }
                    } // grid iteration loop
                }
            }
        }
    } // voxel iteration loop
}

void particle_grid_tree::get_grid_neighborhood_box(
    const frantic::graphics::boundbox3f& worldBounds,
    std::vector<frantic::graphics::raw_byte_buffer*>& outByteBuffers ) const {
    // get the start and end voxels, and the start and end grids
    vector3f voxelCoordMinFloat = m_vcs.get_voxel_coord( worldBounds.minimum() );
    vector3f voxelCoordMaxFloat = m_vcs.get_voxel_coord( worldBounds.maximum() );
    vector3 voxelCoordMin = vector3::from_floor( voxelCoordMinFloat );
    vector3 voxelCoordMax = vector3::from_ceil( voxelCoordMaxFloat );
    vector3 gridCoordMin = vector3::from_floor( voxelCoordMinFloat * PGT_SL );
    vector3 gridCoordMax = vector3::from_ceil( voxelCoordMaxFloat * PGT_SL );

    // iteration variables
    vector3 currVoxelCoord;
    vector3 currGridCoord;

    // iterate through all our voxels
    for( currVoxelCoord.z = voxelCoordMin.z; currVoxelCoord.z < voxelCoordMax.z; ++currVoxelCoord.z ) {
        currGridCoord.z = currVoxelCoord.z * PGT_SL;
        for( currVoxelCoord.y = voxelCoordMin.y; currVoxelCoord.y < voxelCoordMax.y; ++currVoxelCoord.y ) {
            currGridCoord.y = currVoxelCoord.y * PGT_SL;
            for( currVoxelCoord.x = voxelCoordMin.x; currVoxelCoord.x < voxelCoordMax.x; ++currVoxelCoord.x ) {
                currGridCoord.x = currVoxelCoord.x * PGT_SL;

                // NOTE: it could be written so that we use the neighbor node pointers if they're available
                // however, it didn't seem to be a performance issue during profiling.
                pgt_node_type* node = navigate_to_leaf( currVoxelCoord );

                if( node ) {
                    // compute which grid to start and end on for this leaf node
                    vector3 start( gridCoordMin - currGridCoord );
                    vector3 end( gridCoordMax - currGridCoord );
                    vector3 startGrid( start.x > 0 ? start.x : 0, start.y > 0 ? start.y : 0,
                                       start.z > 0 ? start.z : 0 );
                    vector3 endGrid( end.x < PGT_SL ? end.x : PGT_SL, end.y < PGT_SL ? end.y : PGT_SL,
                                     end.z < PGT_SL ? end.z : PGT_SL );

                    // go through each of the byte buffers in the grids and add them to our output vector
                    for( int gz = startGrid.z; gz < endGrid.z; ++gz ) {
                        for( int gy = startGrid.y; gy < endGrid.y; ++gy ) {
                            for( int gx = startGrid.x; gx < endGrid.x; ++gx ) {
                                int gridIndex = gx + ( gy + gz * PGT_SL ) * PGT_SL;
                                raw_byte_buffer* buffer = &node->leafData->grid[gridIndex];
                                if( buffer->ptr_at( 0 ) ) // only add non-empty buffers
                                    outByteBuffers.push_back( buffer );
                            }
                        }
                    } // grid iteration loop
                }
            }
        }
    } // voxel iteration loop
}

void particle_grid_tree::get_next_populated_grid( const std::vector<std::deque<pgt_node_type*>>& leafNodeList,
                                                  size_t& outNextLeafNodeListIndex, size_t& outNextLeafIndex,
                                                  size_t& outNextGridIndex, char*& outNextCurrentParticle,
                                                  char*& outNextParticleEndPtr ) {
    // static function provided strictly for the iterator classes

    // start at the first populated leaf-node list
    while( outNextLeafNodeListIndex < 8 && leafNodeList[outNextLeafNodeListIndex].size() == 0 ) {
        ++outNextLeafNodeListIndex;
        outNextLeafIndex = 0;
        outNextGridIndex = 0;
    }

    // check if this iterator's at the end already
    if( outNextLeafNodeListIndex >= 8 ) {
        outNextCurrentParticle = NULL;
        outNextParticleEndPtr = NULL;
        return;
    }

    // loop until we've found our first populated grid, or there are no more grids to look at
    bool endReached = false;
    while( !endReached ) {

        // first: find the next populated grid
        for( ; outNextGridIndex < PGT_SL * PGT_SL * PGT_SL; ++outNextGridIndex ) {
            particle_grid_tree_leaf_data* leafData = leafNodeList[outNextLeafNodeListIndex][outNextLeafIndex]->leafData;
            outNextCurrentParticle = leafData->grid[outNextGridIndex].begin();
            outNextParticleEndPtr = leafData->grid[outNextGridIndex].end();
            if( leafData->grid[outNextGridIndex].begin() )
                return; // all done, we've found a grid
        }

        // second: if we've overrun our last grids, move to the next leaf node and continue
        outNextGridIndex = 0;
        ++outNextLeafIndex;

        // check if we've overrun our current leaf array
        if( outNextLeafIndex >= leafNodeList[outNextLeafNodeListIndex].size() ) {

            // find the next non-empty leafNodeList index
            ++outNextLeafNodeListIndex;
            outNextLeafIndex = 0;
            while( outNextLeafNodeListIndex < 8 && leafNodeList[outNextLeafNodeListIndex].size() == 0 )
                ++outNextLeafNodeListIndex;

            // third: check to see if we've overrun our last leaf-array
            if( outNextLeafNodeListIndex >= 8 ) {
                // done everything. end is reached
                outNextCurrentParticle = NULL;
                outNextParticleEndPtr = NULL;
                endReached = true; // exit our main loop
            }
        }
    }
}

//
//
// tbb body object for our particle to particle interactions
//
//

detail::particle_grid_tree_particle_particle_interactions_tbb_body::
    particle_grid_tree_particle_particle_interactions_tbb_body(
        std::deque<pgt_node_type*>* leafNodeList,
        frantic::particles::particle_grid_tree::particle_particle_interaction_function_t interactionFunction,
        void* userData, float searchRadius,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor, size_t particleSize )
    : m_leafNodeList( leafNodeList )
    , m_interactionFunction( interactionFunction )
    , m_userData( userData )
    , m_searchRadiusSquared( searchRadius * searchRadius )
    , m_positionAccessor( positionAccessor )
    , m_particleSize( particleSize ) {
    // constructor
}

void detail::particle_grid_tree_particle_particle_interactions_tbb_body::operator()(
    const tbb::blocked_range<size_t>& r ) const {
    // iterate through the range of leaf nodes that tbb specifies
    for( size_t i = r.begin(); i < r.end(); ++i ) {
        pgt_node_type* node = ( *m_leafNodeList )[i];

        // iterate through all the grids for this leaf
        for( int gridZ = 0; gridZ < PGT_SL; ++gridZ ) {
            for( int gridY = 0; gridY < PGT_SL; ++gridY ) {
                for( int gridX = 0; gridX < PGT_SL; ++gridX ) {
                    int gridIndex = gridX + ( gridY + gridZ * PGT_SL ) * PGT_SL;

                    // iterate through all the particles in our current grid
                    char* currParticle = node->leafData->grid[gridIndex].begin();
                    char* gridIterEnd = node->leafData->grid[gridIndex].end();

                    // continue if this grid is not empty
                    if( currParticle ) {

                        // find all the interacting grids with our current grid
                        frantic::graphics::raw_byte_buffer*
                            neighborGrids[13]; // there are 13 neighbouring grids that we can interact with (see wiki)
                        get_grid_particle_particle_neighborhood( node, gridX, gridY, gridZ, neighborGrids );
                        for( ; currParticle != gridIterEnd; currParticle = currParticle + m_particleSize ) {
                            const vector3f& currParticlePos = m_positionAccessor( currParticle );

                            // interact with all the particles in OUR grid (every particle after our current particle)
                            char* interactionParticle = currParticle + m_particleSize; // the next particle in our grid
                            for( ; interactionParticle != gridIterEnd;
                                 interactionParticle = interactionParticle + m_particleSize ) {
                                if( vector3f::distance_squared( currParticlePos,
                                                                m_positionAccessor( interactionParticle ) ) <
                                    m_searchRadiusSquared )
                                    m_interactionFunction( m_userData, currParticle, interactionParticle );
                            }

                            // interact with all the particle from NEIGHBOURING grids
                            for( int i = 0; i < 13; ++i ) {
                                // continue if this grid is populated
                                if( neighborGrids[i] ) {
                                    char* interactionParticle = neighborGrids[i]->begin();
                                    char* interactionGridIterEnd = neighborGrids[i]->end();
                                    for( ; interactionParticle != interactionGridIterEnd;
                                         interactionParticle = interactionParticle + m_particleSize ) {
                                        if( vector3f::distance_squared( currParticlePos,
                                                                        m_positionAccessor( interactionParticle ) ) <
                                            m_searchRadiusSquared )
                                            m_interactionFunction( m_userData, currParticle, interactionParticle );
                                    }
                                }
                            }
                        } // for: every particle in our current grid
                    }
                } // for: X grid
            }     // for: Y grid
        }         // for: Z grid (every grid in our current leaf node)

    } // for: every leaf node
}

frantic::graphics::raw_byte_buffer* detail::particle_grid_tree_particle_particle_interactions_tbb_body::get_grid_buffer(
    particle_grid_tree_leaf_data* leafData, int gridXIndex, int gridYIndex, int gridZIndex ) const {
    // this function handles getting the particle grid for a given leafData object and grid index
    // the real trick here is that if the index is greater than the side length, or less than zero,
    // it will use its neighbour pointer to navigate and get its correct particle grid.

    int neighborXIndex = 1;
    int neighborYIndex = 1;
    int neighborZIndex = 1;
    if( gridXIndex < 0 ) {
        gridXIndex += PGT_SL;
        neighborXIndex = 0;
    } else if( gridXIndex >= PGT_SL ) {
        gridXIndex -= PGT_SL;
        neighborXIndex = 2;
    }
    if( gridYIndex < 0 ) {
        gridYIndex += PGT_SL;
        neighborYIndex = 0;
    } else if( gridYIndex >= PGT_SL ) {
        gridYIndex -= PGT_SL;
        neighborYIndex = 2;
    }
    if( gridZIndex < 0 ) {
        gridZIndex += PGT_SL;
        neighborZIndex = 0;
    } else if( gridZIndex >= PGT_SL ) {
        gridZIndex -= PGT_SL;
        neighborZIndex = 2;
    }

    pgt_node_type* node = leafData->neighborhood[neighborXIndex + ( neighborYIndex + neighborZIndex * 3 ) * 3];
    if( node )
        return &node->leafData->grid[gridXIndex + ( gridYIndex + gridZIndex * PGT_SL ) * PGT_SL];
    return NULL;
}

void detail::particle_grid_tree_particle_particle_interactions_tbb_body::get_grid_particle_particle_neighborhood(
    pgt_node_type* node, int gridXIndex, int gridYIndex, int gridZIndex,
    frantic::graphics::raw_byte_buffer** outByteBuffers ) const {
    // gets the grid buffers of all 13 grids that we need to interact with in our particle-to-particle interactions.
    // see the wiki page for more information.
    particle_grid_tree_leaf_data* leafData = node->leafData;

    // 9 grids below us
    outByteBuffers[0] = get_grid_buffer( leafData, gridXIndex - 1, gridYIndex - 1, gridZIndex - 1 );
    outByteBuffers[1] = get_grid_buffer( leafData, gridXIndex, gridYIndex - 1, gridZIndex - 1 );
    outByteBuffers[2] = get_grid_buffer( leafData, gridXIndex + 1, gridYIndex - 1, gridZIndex - 1 );
    outByteBuffers[3] = get_grid_buffer( leafData, gridXIndex - 1, gridYIndex, gridZIndex - 1 );
    outByteBuffers[4] = get_grid_buffer( leafData, gridXIndex, gridYIndex, gridZIndex - 1 );
    outByteBuffers[5] = get_grid_buffer( leafData, gridXIndex + 1, gridYIndex, gridZIndex - 1 );
    outByteBuffers[6] = get_grid_buffer( leafData, gridXIndex - 1, gridYIndex + 1, gridZIndex - 1 );
    outByteBuffers[7] = get_grid_buffer( leafData, gridXIndex, gridYIndex + 1, gridZIndex - 1 );
    outByteBuffers[8] = get_grid_buffer( leafData, gridXIndex + 1, gridYIndex + 1, gridZIndex - 1 );

    // 4 grids on the same z axis as us
    outByteBuffers[9] = get_grid_buffer( leafData, gridXIndex - 1, gridYIndex - 1, gridZIndex );
    outByteBuffers[10] = get_grid_buffer( leafData, gridXIndex, gridYIndex - 1, gridZIndex );
    outByteBuffers[11] = get_grid_buffer( leafData, gridXIndex + 1, gridYIndex - 1, gridZIndex );
    outByteBuffers[12] = get_grid_buffer( leafData, gridXIndex - 1, gridYIndex, gridZIndex );
}
} // namespace particles
} // namespace frantic
