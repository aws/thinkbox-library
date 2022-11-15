// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <deque>

#include <frantic/channels/channel_map.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/particles/particle_array_iterator.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/particles/streams/particle_ostream.hpp>
#include <frantic/volumetrics/grid_tree_base.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

#include <tbb/blocked_range.h>

// PGT_SL: this definition is the side length for the particle grid tree.
// it is a fixed value to avoid extra pointer indirection when creating nodes.
// setting it to 4 seemed to be the reasonable choice for our test cases.
// it cannot be defined within the class because the class definition depends on it.
#define PGT_SL 4

namespace frantic {
namespace particles {

// forward declarations
class particle_grid_tree;
struct particle_grid_tree_leaf_data;

// handy typedef for our grid tree node type
typedef frantic::volumetrics::grid_tree_base_node<particle_grid_tree_leaf_data, PGT_SL> pgt_node_type;

/**
 * particle_grid_tree_leaf_data
 * This is a struct for leaf data in our tree. Each leaf node in
 * the tree has one of these objects. The particle_grid_tree's parent object (grid_tree_base)
 * is templated on this type.
 */
struct particle_grid_tree_leaf_data {
    // particle "grids". this is the actual particle data held by this leaf
    // our leaf data is a grid of length SIDELENGTH^3. essentially, the
    // tree is always 1 less in depth and the extra data is stored at the leaf.
    // this is done to reduce the tree depth to decrease the traversal time.
    frantic::graphics::raw_byte_buffer grid[PGT_SL * PGT_SL * PGT_SL];

    // next node in our linked list of leaf nodes
    // these form a linked list which is used for fast node iteration
    pgt_node_type* nextNode;

    // for accessing neighbours quickly
    // all 3x3x3 neighbours are stored (including ourself) to reduce tree lookups
    pgt_node_type* neighborhood[27];

    ~particle_grid_tree_leaf_data() {
        for( size_t i = 0; i < PGT_SL * PGT_SL * PGT_SL; ++i )
            grid[i].clear();
    }
};

/**
 * particle_grid_tree
 * This class contains all the logic and data for creating a particle grid tree.
 * It extends a grid_tree_base using a custom leaf data class and compile-time side length
 */
class particle_grid_tree : public frantic::volumetrics::grid_tree_base<particle_grid_tree_leaf_data, PGT_SL> {
  private:
    // the channel map for the particles in the grids
    frantic::channels::channel_map m_channelMap;

    // the position channel accessor for the particles in the grids
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;

    // the voxel coordinate system for the tree nodes
    frantic::volumetrics::voxel_coord_system m_vcs;

    // the voxel coordinate system for the particle grids
    // note that this is different from the tree node's vcs since the tree's leaf nodes
    // store SIDELENGTH^3 grids of particles. the difference is m_vcs = m_gcs * PGT_SL.
    frantic::volumetrics::voxel_coord_system m_gcs;

    // the number of particles. returned by size().
    std::size_t m_particleCount;

    // pointer to the beginning of our leaf-node linked list (for fast iteration purposes)
    pgt_node_type* m_leafNodeLinkedListStart;

    // collection of deques used for multi-threaded particle-to-particle interactions.
    // vector is always of size 8. each vector item contains a deque of leaf nodes that don't
    // interact with leaf nodes in the same deque. this basically amounts to storing odd and
    // even axis leaf nodes in different deques (see wiki page on the multi-threaded algorithm).
    std::vector<std::deque<pgt_node_type*>> m_leafNodeListArray;

  public:
    // This function prototype is used for doing particle-particle interactions.
    typedef void ( *particle_particle_interaction_function_t )( void* userData, char* firstParticle,
                                                                char* secondParticle );

    // This function prototype is used for doing particle-grid interactions.  The grid particles don't necessarily have
    // their positions embedded within them, so the position must be specified separately.
    typedef void ( *particle_grid_interaction_function_t )( void* userData, char* firstParticle,
                                                            const frantic::graphics::vector3f& gridParticlePosition,
                                                            char* gridParticle );

    // This function prototype is used for doing evaluating particles.
    typedef void ( *particle_evaluation_function_t )( void* userData, char* firstParticle );

    // forward declare iterator classes
    class const_iterator;
    class const_node_iterator;

    /**
     * iterator inner class for the particle grid tree.
     */
    class iterator {
      private:
        friend class particle_grid_tree; // friends with the particle grid tree
        friend class const_iterator;     // friends with the const version for implicit conversion
        pgt_node_type* m_currentNode;    // when this is NULL, its the end of the iterator
        size_t m_gridIndex;
        char* m_currentParticle;
        char* m_particleEndPtr;
        size_t m_particleSize;
        iterator( pgt_node_type* currentNode, size_t particleSize ); // private constructor
      public:
        iterator& operator++(); // pre-increment
        bool operator==( const iterator& rhs ) const { return m_currentParticle == rhs.m_currentParticle; }
        bool operator!=( const iterator& rhs ) const { return m_currentParticle != rhs.m_currentParticle; }
        char* operator*() const { return m_currentParticle; }
    };

    /**
     * const_iterator inner class for the particle grid tree. the implementation is identical to the non-const iterator.
     */
    class const_iterator {
      private:
        friend class particle_grid_tree;    // friends with the particle grid tree
        const pgt_node_type* m_currentNode; // when this is NULL, its the end of the iterator
        size_t m_gridIndex;
        const char* m_currentParticle;
        const char* m_particleEndPtr;
        size_t m_particleSize;
        const_iterator( const pgt_node_type* currentNode, size_t particleSize ); // private constructor
      public:
        const_iterator( const iterator& rhs ); // convert from non-const iterator
        const_iterator& operator++();          // pre-increment
        bool operator==( const const_iterator& rhs ) const { return m_currentParticle == rhs.m_currentParticle; }
        bool operator!=( const const_iterator& rhs ) const { return m_currentParticle != rhs.m_currentParticle; }
        const char* operator*() const { return m_currentParticle; }
    };

    /**
     * node_iterator inner class for the particle grid tree.
     */
    class node_iterator {
      private:
        friend class particle_grid_tree;  // friends with the particle grid tree
        friend class const_node_iterator; // friends with the const version for implicit conversion
        pgt_node_type* m_currentNode;     // when this is NULL, its the end of the iterator
        size_t m_gridIndex;
        size_t m_particleSize;
        frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
        frantic::volumetrics::voxel_coord_system m_vcs;
        node_iterator( pgt_node_type* currentNode,
                       const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
                       const frantic::volumetrics::voxel_coord_system& vcs,
                       size_t particleSize );                              // private constructor
        node_iterator& operator=( const node_iterator& ) { return *this; } // disabled assignment operator
      public:
        typedef std::pair<frantic::particles::particle_array_iterator, frantic::particles::particle_array_iterator>
            particle_iterator_pair; // pair of starting and ending particle iterators for this node
        particle_iterator_pair particle_iterators() const; // accessor for our particle iterators
        frantic::graphics::boundbox3f world_bounds() const;
        node_iterator& operator++(); // pre-increment
        bool operator==( const node_iterator& rhs ) const {
            return ( m_currentNode == rhs.m_currentNode ) && ( m_gridIndex == rhs.m_gridIndex );
        }
        bool operator!=( const node_iterator& rhs ) const {
            return ( m_currentNode != rhs.m_currentNode ) || ( m_gridIndex != rhs.m_gridIndex );
        }
    };

    /**
     * const_node_iterator inner class for the particle grid tree. the implementation is identical to the non-const
     * iterator.
     */
    class const_node_iterator {
      private:
        friend class particle_grid_tree;    // friends with the particle grid tree
        const pgt_node_type* m_currentNode; // when this is NULL, its the end of the iterator
        size_t m_gridIndex;
        size_t m_particleSize;
        frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
        frantic::volumetrics::voxel_coord_system m_vcs;
        const_node_iterator( const pgt_node_type* currentNode,
                             const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
                             const frantic::volumetrics::voxel_coord_system& vcs,
                             size_t particleSize );                              // private constructor
        const_node_iterator& operator=( const node_iterator& ) { return *this; } // disabled assignment operator
      public:
        const_node_iterator( const node_iterator& rhs ); // convert from non-const version
        typedef std::pair<frantic::particles::const_particle_array_iterator,
                          frantic::particles::const_particle_array_iterator>
            const_particle_iterator_pair; // pair of starting and ending particle iterators for this node
        const_particle_iterator_pair const_particle_iterators() const; // accessor for our particle iterators
        frantic::graphics::boundbox3f world_bounds() const;
        const_node_iterator& operator++(); // pre-increment
        bool operator==( const const_node_iterator& rhs ) const {
            return ( m_currentNode == rhs.m_currentNode ) && ( m_gridIndex == rhs.m_gridIndex );
        }
        bool operator!=( const const_node_iterator& rhs ) const {
            return ( m_currentNode != rhs.m_currentNode ) || ( m_gridIndex != rhs.m_gridIndex );
        }
    };

  public:
    //
    //
    // Tree construction
    //
    //

    /*
     * default constructor
     */
    particle_grid_tree();

    /**
     * constructor
     * @param  pcm  the particle channel map to use.
     */
    particle_grid_tree( const frantic::channels::channel_map& pcm );

    /**
     * constructor
     * @param  pcm  the particle channel map to use.
     * @param  voxelLength  the voxel length of the particle grids.
     */
    particle_grid_tree( const frantic::channels::channel_map& pcm, float voxelLength );

    /**
     * constructor with a particular particle channel map, and voxel coordinate system.
     * @param  pcm  the particle channel map to use.
     * @param  vcs  the voxel coordinate system of the particle grids.
     */
    particle_grid_tree( const frantic::channels::channel_map& pcm,
                        const frantic::volumetrics::voxel_coord_system& vcs );

    /**
     * destructor
     */
    ~particle_grid_tree();

    /**
     * clears all the particles from the tree
     */
    void clear();

    /**
     * clears the particles and starts a new tree
     * @param  pcm  the new particle channel map to use.
     * @param  voxelLength  the new voxel length of the particle grids.
     */
    void reset( const frantic::channels::channel_map& pcm, float voxelLength );

    /**
     * clears the particles and starts a new tree
     * @param  pcm  the particle channel map to use.
     * @param  vcs  the voxel coordinate system of the particle grids.
     */
    void reset( const frantic::channels::channel_map& pcm, const frantic::volumetrics::voxel_coord_system& vcs );

    /**
     * swaps this particle grid with another
     * @param  rhs  the tree to swap with.
     */
    void swap( particle_grid_tree& rhs );

    //
    //
    // Tree iteration
    //
    //

    /**
     * access to the begining particle iterator for our tree
     * @note Writable iterators MUST NOT be used to modify the position of the particles in the particle_grid_tree.
     *       Doing so will produce incorrect results
     * @return  an iterator object set to the first particle
     */
    iterator begin() {
        return particle_grid_tree::iterator( m_leafNodeLinkedListStart, m_channelMap.structure_size() );
    }

    /**
     * access to the one-past-end particle iterator for our tree
     * @return  an iterator object set to the one-past-end particle
     */
    iterator end() { return particle_grid_tree::iterator( NULL, m_channelMap.structure_size() ); }

    /**
     * access to the begining particle const iterator for our tree
     * @return  a const_iterator object set to the first particle
     */
    const_iterator begin() const {
        return particle_grid_tree::const_iterator( m_leafNodeLinkedListStart, m_channelMap.structure_size() );
    }

    /**
     * access to the one-past-end particle const iterator for our tree
     * @return  a const_iterator object set to the one-past-end particle
     */
    const_iterator end() const { return particle_grid_tree::const_iterator( NULL, m_channelMap.structure_size() ); }

    /**
     * access to the begining node iterator for our tree
     * @note Writable iterators MUST NOT be used to modify the position of the particles in the particle_grid_tree.
     *       Doing so will produce incorrect results
     * @return  a node_iterator object set to the first grid
     */
    node_iterator nodes_begin() {
        return particle_grid_tree::node_iterator( m_leafNodeLinkedListStart, m_positionAccessor, m_vcs,
                                                  m_channelMap.structure_size() );
    }

    /**
     * access to the one-past-end node iterator for our tree
     * @return  an node_iterator object set to the one-past-end node
     */
    node_iterator nodes_end() {
        return particle_grid_tree::node_iterator( NULL, m_positionAccessor, m_vcs, m_channelMap.structure_size() );
    }

    /**
     * access to the begining node const iterator for our tree
     * @return  a const_node_iterator object set to the first grid
     */
    const_node_iterator const_nodes_begin() const {
        return particle_grid_tree::const_node_iterator( m_leafNodeLinkedListStart, m_positionAccessor, m_vcs,
                                                        m_channelMap.structure_size() );
    }

    /**
     * access to the one-past-end node const iterator for our tree
     * @return  an const_node_iterator object set to the one-past-end node
     */
    const_node_iterator const_nodes_end() const {
        return particle_grid_tree::const_node_iterator( NULL, m_positionAccessor, m_vcs,
                                                        m_channelMap.structure_size() );
    }

    //
    //
    // Accessors for tree data
    //
    //

    /**
     * vcs accessor
     * @return  the voxel coordinate system of our particle grids.
     */
    const frantic::volumetrics::voxel_coord_system& get_voxel_coord_system() const {
        return m_gcs; // not the vcs, that's for the leaf nodes in the tree
    }

    /**
     * channel map accessor
     * @return  the channel map for our particles
     */
    const frantic::channels::channel_map& get_channel_map() const { return m_channelMap; }

    /**
     * size accessor
     * @return  the number of particles in the tree
     */
    std::size_t particle_count() const { return m_particleCount; }

    /**
     * a synonym for particle_count (to be STL container-like)
     * @return  the number of particles in the tree
     */
    std::size_t size() const { return m_particleCount; }

    /**
     * computes the bounds of the particles in the tree.
     * this call is expensive and should not be called in performance situations.
     * @return  voxel bounding box
     */
    frantic::graphics::boundbox3f compute_particle_bounds() const;

    //
    //
    // Particle access
    //
    //

    /**
     * Searches to see if there is one or more particles near a given location
     * @param  worldLocation  The location to look for particles
     * @param  searchRadius  The radius around the world location to look
     * @return  boolean result
     */
    bool has_particle_near( const frantic::graphics::vector3f& worldLocation, float searchRadius ) const;

    /**
     * Finds and returns the closest particle to a given location
     * @param  worldLocation  The location to look for the nearest particle.
     * @param  searchRadius  The maximum radius to look for a particle.
     * @param  outParticle  An output parameter containing the data of the nearest particle. This must be allocated by
     * the user.
     * @param  outDistance  The distance at which the nearest particle was found.
     * @return  whether or not a particle was found within the given search radius.
     */
    bool get_closest_particle( const frantic::graphics::vector3f& worldLocation, float searchRadius, char* outParticle,
                               float& outDistance ) const;

    /**
     * This returns pointers to all the particles which are within searchRadius of worldLocation.
     * @param  worldLocation  the world location around which to search.
     * @param  searchRadius  the search radius.
     * @param  outParticles  the output parameter into which the particles will be placed. The pointers point directly
     * into the tree's data.
     */
    void get_particles_in_range( const frantic::graphics::vector3f& worldLocation, float searchRadius,
                                 std::vector<char*>& outParticles ) const;

    /**
     * This returns pointers to all the particles which are within searchRadius of the worldBox bounding box.
     * @param  worldBox  the bounding box in world space around which to search.
     * @param  searchRadius  the search radius.
     * @param  outParticles  the output parameter into which the particles will be placed.
     */
    void get_particles_in_range( const frantic::graphics::boundbox3f& worldBox, float searchRadius,
                                 std::vector<char*>& outParticles ) const;

    /**
     * Particle interation with the particle grid tree:
     * This function calls the provided function pointer on all the pairs of particles that are closer together than the
     * provided search radius.  It does not call the function pointer with a particle repeated twice, interactions of a
     * particle with itself need to be handled separately.
     * IMPORTANT NOTE: The following constraints MUST be followed by the provided function for correctness:
     *                  * The function MUST NOT modify the "Position" channel.
     *                  * The function MUST ONLY add or subtract from destination channels.
     *                  * The function MUST NOT write to any channel it is also reading from.
     *
     * @param  interactionFunction  the function that is called on each particle pair.
     * @param  userData  user data provided for the interaction function.
     * @param  searchRadius  the radius of interaction between particles. must be less than the particle grid's length.
     */
    void particle_particle_interactions( particle_particle_interaction_function_t interactionFunction, void* userData,
                                         float searchRadius );

    /**
     * Calls the user function provided on every particle within a defined bounds.
     * @param  userData  user data provided for the evaluation function.
     * @param  worldBox  the bounding box for the particles which we need to call the evaluation function on.
     * @param  evaluationFunction  the function which will be evaluated for each particle in the given bounds.
     */
    void process_particles_in_bounds( void* userData, const frantic::graphics::boundbox3f& worldBox,
                                      particle_evaluation_function_t evaluationFunction );

    /**
     *  Calls the user function provided on every particle within a defined bounds
     * using multiple worker threads.
     *
     * @note multiple threads may access userData simultaneously.  You are responsible
     *		for ensuring thread safe access to the userData.
     *
     * @param  userData  user data provided for the evaluation function.  You are responsible for ensuring thread safe
     *access to this data.
     * @param  worldBox  the bounding box in world space in which to search for particles.
     * @param  evaluationFunction  the function which will be evaluated for each particle in the given bounds.
     */
    void process_particles_in_bounds_mt( void* userData, const frantic::graphics::boundbox3f& worldBox,
                                         particle_evaluation_function_t evaluationFunction );

    /**
     * Calls the user function provided on every particle within searchRadius of worldPosition.
     * @param  userData  user data provided for the evaluation function.
     * @param  worldPosition  the world position around which to search.
     * @param  searchRadius  the search radius.
     * @param  evaluationFunction  the function which will be evaluated for each particle in the given bounds.
     */
    void process_particles_in_range( void* userData, const frantic::graphics::vector3f& worldPosition,
                                     float searchRadius, particle_evaluation_function_t evaluationFunction );

    /**
     * Particle interaction with a regular grid of particles:
     * This function interacts with an array of particles defined on a regular grid. The interaction function is
     * called with first parameter from the particle grid tree, and second parameter from the regular grid of particles.
     * The function is intended for use in interacting with dense grid structures, such as a dense grid level set.
     * IMPORTANT NOTE: The following constraints MUST be followed by the provided function for correctness:
     *                  * The function MUST NOT modify the "Position" channel of the particle from the
     * particle_grid_tree.
     *                  * The function MUST ONLY add or subtract from destination channels.
     *                  * The function MUST NOT write to any channel it is also reading from.
     *
     * NOTE: THIS FUNCTION IS TO BE REMOVED ONCE THE UNUSED PARTICLE IMPLICIT POLICY THAT USED IT IS CLEANED OUT!
     *
     * @param  voxelGrid  the bounds of the grid that we are interacting with.
     * @param  gridVCS  the voxel coordinate system of the grid.
     * @param  gridParticleSize  the number of bytes allocated for each particle in our regular grid.
     * @param  gridParticleArray  the buffer containing the particle data for our regular grid.
     * @param  interactionFunction  the function to is called on each particle pair.
     * @param  userData  user data provided for the interaction function.
     * @param  searchRadius  the radius of interaction between particles. must be less than the particle grid's length.
     */
    void particle_grid_interactions( const frantic::graphics::boundbox3& voxelGrid,
                                     const frantic::volumetrics::voxel_coord_system& gridVCS,
                                     std::size_t gridParticleSize, void* gridParticleArray,
                                     particle_grid_interaction_function_t interactionFunction, void* userData,
                                     float searchRadius );

    //
    //
    // Particle insertion functions
    //
    //

    /**
     * Adds the given particle to the tree.
     * @param  rawParticleData  A pointer to the particle to insert.
     */
    void insert( const char* rawParticleData );

    /**
     * Adds the given particle to the tree, subject to the constraint that no particle already in the tree
     * is within searchRadius of it.
     * @param  rawParticleData  A pointer to the particle to insert.
     * @param  searchRadius  The radius around the particle within which no other particle must exist.
     */
    bool insert_with_proximity_constraint( const char* rawParticleData, float searchRadius );

    /**
     * Insert the particle with cyclic boundary conditions, as well as a constraint that no particle already in
     * the tree is within searchRadius of it.
     * @param  rawParticleData  A pointer to the particle to insert.
     * @param  searchRadius  The radius around the particle within which no other particle must exist.
     * @param  cyclicBounds  A bounding box defining the cyclic bounds.
     */
    bool insert_with_proximity_constraint( const char* rawParticleData, float searchRadius,
                                           const frantic::graphics::boundbox3f& cyclicBounds );

    /**
     * Inserts particles from a second grid tree into our tree.
     * @param  theParticles  tree containing the particles to insert.
     */
    void insert_particles( const particle_grid_tree& theParticles );

    /**
     * Inserts particles from a second grid tree into our tree, subject to the
     * constraint that no particle already in the tree is within searchRadius of it.
     * @param  theParticles  tree containing the particles to insert.
     * @param  searchRadius  The radius around the particle within which no other particle must exist.
     */
    void insert_particles_with_proximity_constraint( const particle_grid_tree& theParticles, float searchRadius );

    /**
     * Inserts particles from a particle_istream into our tree.
     * @param  pin  The stream containing the particles to insert.
     */
    inline void insert_particles( boost::shared_ptr<frantic::particles::streams::particle_istream> pin ) {
        frantic::logging::null_progress_logger nullProgress;
        insert_particles( pin, nullProgress );
    }

    /**
     * Inserts particles from a particle_istream into our tree.
     * @param  pin  The stream containing the particles to insert.
     * @param  progress  A progress_logger for logging the progress of the particle insertions.
     */
    void insert_particles( boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                           frantic::logging::progress_logger& progress );

    /**
     * Write out all the particles in the tree to the given particle_ostream.
     * @param  pout  The stream to write the particles to.
     */
    void write_particles( boost::shared_ptr<frantic::particles::streams::particle_ostream> pout ) const;

  private:
    /// Internal function to add a leaf data object to a newly created leaf node.
    void add_new_leaf_data( const frantic::graphics::vector3& voxelCoordinate, pgt_node_type* leafNode );

    /// Evaluate the UnaryFunction for all grids within the radius of a given point.
    /// Note that this function is optimized to use neighbour leaf node pointers if it is possible.
    ///
    /// @param UnaryFunction a bool( raw_byte_buffer* ) function to evaluate for each
    ///			grid within the radius of the worldCoord.  This function must return
    ///			true to continue iteration, or false to stop iteration.
    template <class UnaryFunction>
    void for_each_grid_in_neighborhood( const frantic::graphics::vector3f& worldCoord, float radius,
                                        UnaryFunction f ) const {
        using namespace frantic::graphics;

        // bounds for this lookup
        boundbox3f worldBounds( worldCoord - vector3f( radius ), worldCoord + vector3f( radius ) );

        if( radius > m_vcs.voxel_length() ) {
            // if our lookup radius is larger than the voxel length, we can't use the simple
            // neighbour pointer lookup that our current function provides. instead we must call this function.
            for_each_grid_in_neighborhood_box( worldBounds, f );
            return;
        }

        // get center voxel (the one from which to do neighbour lookups)
        vector3 voxelCoord = vector3::from_floor( m_vcs.get_voxel_coord( worldCoord ) );
        pgt_node_type* leafNode = navigate_to_leaf( voxelCoord );

        if( !leafNode ) {
            // if we don't have a leaf directly in the lookup location, we still need to check the
            // neighbours. except we can't use the handy neighbour pointers, so we call this funciton.
            for_each_grid_in_neighborhood_box( worldBounds, f );
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
                        vector3 startGrid( std::max<boost::int32_t>( start.x, 0 ),
                                           std::max<boost::int32_t>( start.y, 0 ),
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
                                    if( buffer->ptr_at( 0 ) ) { // only process non-empty buffers
                                        if( !f( buffer ) )
                                            return;
                                    }
                                }
                            }
                        } // grid iteration loop
                    }
                }
            }
        } // voxel iteration loop
    }

    /// Evaluate the UnaryFunction for all grids for the given bounding box.
    /// This function is not optimized with neighbour pointer lookups like the one above.
    ///
    /// @param UnaryFunction a bool( raw_byte_buffer* ) function to evaluate for each
    ///			grid on or inside the bounding box.  This function must return
    ///			true to continue iteration, or false to stop iteration.
    template <class UnaryFunction>
    void for_each_grid_in_neighborhood_box( const frantic::graphics::boundbox3f& worldBounds, UnaryFunction f ) const {
        using namespace frantic::graphics;

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
                                    if( buffer->ptr_at( 0 ) ) { // only process non-empty buffers
                                        if( !f( buffer ) )
                                            return;
                                    }
                                }
                            }
                        } // grid iteration loop
                    }
                }
            }
        } // voxel iteration loop
    }

    /// Internal function to get a neighbourhood of grids within the radius of a given point.
    /// Note that this function is optimized to use neighbour leaf node pointers if it is possible.
    void get_grid_neighborhood( const frantic::graphics::vector3f& worldCoord, float radius,
                                std::vector<frantic::graphics::raw_byte_buffer*>& outByteBuffers ) const;

    /// Internal function to get a neighbourhood of grids for a given bounding box.
    /// This function is not optimized with neighbour pointer lookups like the one above.
    void get_grid_neighborhood_box( const frantic::graphics::boundbox3f& worldBounds,
                                    std::vector<frantic::graphics::raw_byte_buffer*>& outByteBuffers ) const;

    /// This static function is provided strictly for use by the iterators. it allows us to get at the next populated
    /// particle grid.
    static void get_next_populated_grid( const std::vector<std::deque<pgt_node_type*>>& leafNodeList,
                                         size_t& outNextLeafNodeListIndex, size_t& outNextLeafIndex,
                                         size_t& outNextGridIndex, char*& outNextCurrentParticle,
                                         char*& outNextParticleEndPtr );
};

//
//
// tbb body object for our particle to particle interactions
//
//

namespace detail {
/**
 * This class is a "body" class for tbb. its purpose is to do particle to particle interactions.
 * One of these objects is created in the pgt's particle_particle_interaction function and executed
 * with tbb's parallel_for.
 */
class particle_grid_tree_particle_particle_interactions_tbb_body {
  private:
    std::deque<pgt_node_type*>* m_leafNodeList;
    frantic::particles::particle_grid_tree::particle_particle_interaction_function_t m_interactionFunction;
    void* m_userData;
    float m_searchRadiusSquared;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    size_t m_particleSize;

  public:
    /**
     * constructor- just takes all the parameters for the particle to particle interaction function and which data to
     * apply it to
     */
    particle_grid_tree_particle_particle_interactions_tbb_body(
        std::deque<pgt_node_type*>* leafNodeList,
        frantic::particles::particle_grid_tree::particle_particle_interaction_function_t interactionFunction,
        void* userData, float searchRadius,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor, size_t particleSize );

    /// called by tbb to start processing: this function does particle particle interactions for a range decided by tbb
    void operator()( const tbb::blocked_range<size_t>& r ) const;

  private:
    /// Internal function used in get_grid_particle_particle_neighborhood to get a byte buffer given a leaf and grid
    /// index.
    frantic::graphics::raw_byte_buffer* get_grid_buffer( particle_grid_tree_leaf_data* leafData, int gridXIndex,
                                                         int gridYIndex, int gridZIndex ) const;

    /// Internal function to get the grids that interact with a defined grid for particle-particle interactions.
    void get_grid_particle_particle_neighborhood( pgt_node_type* node, int gridXIndex, int gridYIndex, int gridZIndex,
                                                  frantic::graphics::raw_byte_buffer** outByteBuffers ) const;
};
} // namespace detail

} // namespace particles
} // namespace frantic
