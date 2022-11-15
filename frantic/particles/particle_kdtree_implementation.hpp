// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

template <class P>
inline frantic::particles::particle_kdtree<P>::particle_kdtree() {
    clear();
    m_KDTreeBuiltFlag = false;
}

template <class P>
inline frantic::particles::particle_kdtree<P>::~particle_kdtree() {}

template <class P>
inline void frantic::particles::particle_kdtree<P>::clear() {
    m_particles.clear();
    m_axes.clear();
    // m_elems.clear();
    m_bounds.set_to_empty();
    m_KDTreeBuiltFlag = false;
}

// This adds a single particle to the kd-tree.  Once all the particles are added, balance_kdtree should be called.
// Adding additional particles resets the balanced flag of the structure.
template <class P>
inline void frantic::particles::particle_kdtree<P>::add_particle( const P& particle ) {
    m_particles.push_back( particle );
    // m_elems.push_back(std::make_pair(particle,0));
    m_KDTreeBuiltFlag = false;
}

/*
 *	Balance: Must be called before particles can be located. Creates a kd-tree.
 */

template <class P>
inline void frantic::particles::particle_kdtree<P>::balance_kdtree() {
    //	frantic::diagnostics::profiling_section pf_bbox("compute bounding box"),
    //		pf_allocate_temp("allocate temp vectors"), pf_balance("balance"), pf_copy("copy to result");

    //	pf_bbox.enter();
    // Compute the bounding box, in case the get_particle_ref function was used
    m_bounds.set_to_empty();
    for( unsigned i = 0; i < m_particles.size(); ++i )
        m_bounds += m_particles[i];
    //	pf_bbox.exit();

    m_axes.resize( m_particles.size() );

    if( m_particles.size() > 1 ) {
        //		pf_allocate_temp.enter();
        std::vector<P> kdTree( m_particles.size() );
        //		pf_allocate_temp.exit();

        //		pf_balance.enter();
        // TODO: Make a second balance_segment function which uses pointers, and use either the by-value one
        //       or the by-pointer one depending on sizeof(P).  For vector3f's, with 12 bytes, the by-value
        //       version is significantly faster.  For larger particle classes, the by-pointer one will likely
        //       be faster.
        // Balance the tree, filling the allocated array with the kd-tree as we go.
        balance_segment( kdTree, 0, 0, (int)size() - 1 );
        //		pf_balance.exit();

        //		pf_copy.enter();
        // Swap the finished kd-tree into the primary particles array
        m_particles.swap( kdTree );
        //		pf_copy.exit();
    }

    m_KDTreeBuiltFlag = true;

    //	if( frantic::logging::is_logging_stats() ) {
    //		std::cerr << pf_bbox << std::endl;
    //		std::cerr << pf_allocate_temp << std::endl;
    //		std::cerr << pf_balance << std::endl;
    //		std::cerr << pf_copy << std::endl;
    //	}
}

template <class P>
inline void frantic::particles::particle_kdtree<P>::balance_segment( std::vector<P>& outKDTree, int index, int start,
                                                                     int end ) {
    int median = 1;

    // Find the desired median index for a left-balanced kd-tree
    while( ( 4 * median ) <= ( end - start + 1 ) )
        median += median;

    if( ( 3 * median ) <= ( end - start + 1 ) ) {
        median += median;
        median += start - 1;
    } else {
        median = end - median + 1;
    }

    // Find the axis to split on.  During creation, the global bounding box
    // is updated to contain the bounds of the region being processed.
    int axis = m_bounds.get_max_dimension_axis();

    // Partition the segment on the median
    median_split( m_particles, start, end, median, axis );

    //	std::cerr << "\n\n" << start << ", " << median << ", " << end << " axis " << axis << ", bounds " << m_bounds;

    //	std::cerr << "\n";
    //	for( int i = start; i <= end; ++i )
    //		std::cerr << m_particles[i] << " ";

    // Fill in this element in the output kd-tree
    outKDTree[index] = m_particles[median];
    m_axes[index] = (char)axis;

    // Recursively balance the left and right block
    if( median > start ) // Left
    {
        if( start < median - 1 ) {
            const float temp = m_bounds.maximum()[axis];
            m_bounds.maximum()[axis] = outKDTree[index][axis];
            balance_segment( outKDTree, 2 * index + 1, start, median - 1 );
            m_bounds.maximum()[axis] = temp;
        } else {
            outKDTree[2 * index + 1] = m_particles[start];
            m_axes[2 * index + 1] = 0;
        }
    }

    if( median < end ) // Right
    {
        if( median + 1 < end ) {
            const float temp = m_bounds.minimum()[axis];
            m_bounds.minimum()[axis] = outKDTree[index][axis];
            balance_segment( outKDTree, 2 * index + 2, median + 1, end );
            m_bounds.minimum()[axis] = temp;
        } else {
            outKDTree[2 * index + 2] = m_particles[end];
            m_axes[2 * index + 2] = 0;
        }
    }
}

// This function pushes elements leftward and rightward so that the value at pv[median] is greater or equal to
// all the values pv[i] for i < median, and less than or equal to all the values pv[i] for i > median.  This algorithm
// has an expected linear time.
//
// TODO: Switch to calling the more generic median_split implementation in threaded_sort.hpp
template <class P>
inline void frantic::particles::particle_kdtree<P>::median_split( std::vector<P>& particles, int start, int end,
                                                                  int median, int axis ) {
    // Insertion sort when fewer elements
    if( end <= start + 2 ) {
        for( int i = start; i < end; ++i ) {
            float v = particles[i][axis];
            int vi = i;
            for( int j = i + 1; j <= end; ++j ) {
                if( particles[j][axis] < v ) {
                    v = particles[j][axis];
                    vi = j;
                }
            }
            if( vi != i )
                std::swap( particles[i], particles[vi] );
        }

    } else {
        int left = start;
        int right = end;

        while( right > left ) {
            int pivotIndex = rand() % ( right - left + 1 ) + left;

            const float v = particles[pivotIndex][axis];
            int i = left - 1;
            int j = right;

            // Randomizing the pivot improved the performance by about 30% in some tests.
            if( pivotIndex != right )
                std::swap( particles[pivotIndex], particles[right] );
            // std::swap( pv[pivotIndex], pv[right] );

            for( ;; ) {
                ++i;
                while( particles[i][axis] < v )
                    ++i;

                --j;
                while( particles[j][axis] > v && j > left )
                    --j;

                if( i >= j )
                    break;

                // std::swap( pv[i], pv[j] );
                std::swap( particles[i], particles[j] );
            }

            // std::swap( pv[i], pv[right] );
            std::swap( particles[i], particles[right] );

            if( i >= median )
                right = i - 1;

            if( i <= median )
                left = i + 1;
        }
    }
}

// Locates all the particles within the specified range (squared), and adds them to the provided vector
template <class P>
inline void frantic::particles::particle_kdtree<P>::locate_particles_range(
    const frantic::graphics::vector3f& position, float rangeSquared,
    std::vector<std::pair<float, P>>& outLocatedParticles ) const {
    if( !m_KDTreeBuiltFlag )
        throw std::runtime_error( "Tried to locate particles in a kdtree that hadn't been balanced yet." );

    if( m_particles.size() > 0 )
        recursive_locate_particles_range( position, rangeSquared, outLocatedParticles, 0 );
}

// Recursively locates all the particles within the specified range (squared), and adds them to the provided vector
template <class P>
inline void frantic::particles::particle_kdtree<P>::recursive_locate_particles_range(
    const frantic::graphics::vector3f& position, float rangeSquared,
    std::vector<std::pair<float, P>>& outLocatedParticles, int index ) const {
    const P& p = m_particles[index];
    int axis = m_axes[index];

    if( 2 * index + 1 < (int)m_particles.size() ) {
        float axisDistance = position[axis] - p[axis];

        if( axisDistance > 0 ) // search right plane
        {
            if( 2 * index + 2 < (int)m_particles.size() )
                recursive_locate_particles_range( position, rangeSquared, outLocatedParticles, 2 * index + 2 );
            if( axisDistance * axisDistance < rangeSquared )
                recursive_locate_particles_range( position, rangeSquared, outLocatedParticles, 2 * index + 1 );

        } else // dist1 is negative, search left first
        {
            recursive_locate_particles_range( position, rangeSquared, outLocatedParticles, 2 * index + 1 );
            if( 2 * index + 2 < (int)m_particles.size() && axisDistance * axisDistance < rangeSquared )
                recursive_locate_particles_range( position, rangeSquared, outLocatedParticles, 2 * index + 2 );
        }
    }

    float particleDistanceSquared = frantic::graphics::vector3f::distance_squared( p, position );

    if( particleDistanceSquared < rangeSquared ) // It's a particle within range
        outLocatedParticles.push_back( std::make_pair( particleDistanceSquared, p ) );
}

// Locates all the particles within the specified range (squared), and adds them to the provided vector
template <class P>
inline void frantic::particles::particle_kdtree<P>::locate_particles_count(
    const frantic::graphics::vector3f& position, int particleCount,
    std::vector<std::pair<float, P>>& outLocatedParticles ) const {

    if( !m_KDTreeBuiltFlag )
        throw std::runtime_error( "Tried to locate particles in a kdtree that hadn't been balanced yet." );

    typename detail::locate_particles_count_sorter<P>::priority_queue particles;

    if( m_particles.size() > 0 )
        recursive_locate_particles_count( position, particleCount, particles, 0 );

    // Allocate space for the closest particles we found
    outLocatedParticles.resize( outLocatedParticles.size() + particles.size() );
    // Pop all the particles, and add them to the outLocatedParticles vector in reverse order so the closest one ends up
    // first
    unsigned setPosition = (unsigned)outLocatedParticles.size() - 1;
    while( particles.size() > 0 ) {
        outLocatedParticles[setPosition--] = particles.top();
        particles.pop();
    }
}

// Recursively locates all the particles, maintaining the particleCount closest in the priority queue
template <class P>
void frantic::particles::particle_kdtree<P>::recursive_locate_particles_count(
    const frantic::graphics::vector3f& position, int particleCount,
    typename detail::locate_particles_count_sorter<P>::priority_queue& outParticles, int index ) const {
    const P& p = m_particles[index];
    int axis = m_axes[index];

    if( 2 * index + 1 < (int)m_particles.size() ) {
        float axisDistance = position[axis] - p[axis];

        if( axisDistance > 0 ) {
            // First search the right branch
            if( 2 * index + 2 < (int)m_particles.size() )
                recursive_locate_particles_count( position, particleCount, outParticles, 2 * index + 2 );
            // search the left branch if we need more particles or if it's close enough that it could replace some
            // particles
            if( (int)outParticles.size() < particleCount || axisDistance * axisDistance < outParticles.top().first )
                recursive_locate_particles_count( position, particleCount, outParticles, 2 * index + 1 );

        } else {
            // First search the left branch
            recursive_locate_particles_count( position, particleCount, outParticles, 2 * index + 1 );
            // search the right branch if we need more particles or if it's close enough that it could replace some
            // particles
            if( 2 * index + 2 < (int)m_particles.size() &&
                ( (int)outParticles.size() < particleCount || axisDistance * axisDistance < outParticles.top().first ) )
                recursive_locate_particles_count( position, particleCount, outParticles, 2 * index + 2 );
        }
    }

    float particleDistanceSquared = frantic::graphics::vector3f::distance_squared( p, position );

    if( (int)outParticles.size() < particleCount ) {
        outParticles.push( std::make_pair( particleDistanceSquared, p ) );
    } else if( particleDistanceSquared < outParticles.top().first ) {
        // remove the most distant particle
        outParticles.pop();
        // and replace it with the one we found
        outParticles.push( std::make_pair( particleDistanceSquared, p ) );
    }
}

// Locates the first particle within the specified range (squared), and returns true as soon as one is found. otherwise,
// it will return false.
template <class P>
inline bool frantic::particles::particle_kdtree<P>::proximity_test( const frantic::graphics::vector3f& position,
                                                                    float rangeSquared ) const {
    if( !m_KDTreeBuiltFlag )
        throw std::runtime_error( "Tried to locate particles in a kdtree that hadn't been balanced yet." );

    if( m_particles.size() > 0 )
        return recursive_proximity_test( position, rangeSquared, 0 );
    return false;
}

// Recursively locates the first particle within the specified range (squared), and returns true as soon as one is
// found. otherwise, it will return false.
template <class P>
inline bool
frantic::particles::particle_kdtree<P>::recursive_proximity_test( const frantic::graphics::vector3f& position,
                                                                  float rangeSquared, int index ) const {
    const P& p = m_particles[index];
    int axis = m_axes[index];

    if( 2 * index + 1 < (int)m_particles.size() ) {
        float axisDistance = position[axis] - p[axis];

        if( axisDistance > 0 ) { // search right plane
            if( 2 * index + 2 < (int)m_particles.size() )
                if( recursive_proximity_test( position, rangeSquared, 2 * index + 2 ) )
                    return true;
            if( axisDistance * axisDistance < rangeSquared )
                if( recursive_proximity_test( position, rangeSquared, 2 * index + 1 ) )
                    return true;
        } else { // dist1 is negative, search left first
            if( recursive_proximity_test( position, rangeSquared, 2 * index + 1 ) )
                return true;
            if( 2 * index + 2 < (int)m_particles.size() && axisDistance * axisDistance < rangeSquared )
                if( recursive_proximity_test( position, rangeSquared, 2 * index + 2 ) )
                    return true;
        }
    }

    float particleDistanceSquared = frantic::graphics::vector3f::distance_squared( p, position );
    if( particleDistanceSquared < rangeSquared ) // It's a particle within range
        return true;

    return false;
}

template <class Type>
struct AlwaysTrue {
    bool operator()( const Type& ) const { return true; }
};

// Finds the nearest particle to the specified position
// Benchmarks were done to compare this speed, BUT, we think these benchmarks were invalid, so retesting
// this versus an implementation with no predicate may be worthwhile.
template <class P>
inline const P&
frantic::particles::particle_kdtree<P>::locate_closest_particle( const frantic::graphics::vector3f& position ) const {
    return locate_closest_particle_if( position, AlwaysTrue<P>() );
}

// Finds the single nearest particle to the specified position. Using the predicate version may return a particle
// that does not actually satisfy the predicate if there are no particles that do. Always check the returned particle
// against the predicate again. Or change this to use pointers. Your call.
// BARF: locate_closest_particle_if moved into the particle_kdtree.hpp file because of Visual Studio 6.

// Recursively finds the nearest particle to the specified position
// BARF: recursive_locate_closest_particle_if moved into the particle_kdtree.hpp file because of Visual Studio 6.

// This uses the kd-tree to accelerate the support function computation.  I did a benchmark
// against a procedure which simply ran through all the points to find the support.  See
// the support function unit test in the graphics test suite for the code.  The points being
// used were based on a random gaussian distribution.
// Here are the results:
//  Count |  KD-Tree | Naive
//  ------------------------
//    100 | 0.0066ms | 0.0069ms
//     1K |   0.02ms |   0.07ms
//    10K |   0.07ms |   0.65ms
//   100K |   0.19ms |    7.3ms
//     1M |   0.31ms |   73.9ms
//    10M |   0.78ms | not tested
template <class P>
inline const P& frantic::particles::particle_kdtree<P>::support( const frantic::graphics::vector3f& d ) const {
    using namespace frantic::graphics;

    if( m_particles.size() == 0 ) {
        throw std::runtime_error( "Tried to call particle_kdtree::support on an empty kd-tree." );
    } else if( !m_KDTreeBuiltFlag ) {
        throw std::runtime_error( "Tried to call particle_kdtree::support on a kd-tree before it was balanced." );
    } else {
        if( m_particles.size() < 50 ) {
            const P* resultSupport = &m_particles[0];
            float resultDotProd = vector3f::dot( *resultSupport, d );
            for( unsigned i = 1; i < m_particles.size(); ++i ) {
                const P* testSupport = &m_particles[i];
                float testDotProd = vector3f::dot( *testSupport, d );
                if( testDotProd > resultDotProd ) {
                    resultSupport = testSupport;
                    resultDotProd = testDotProd;
                }
            }
            return *resultSupport;
        } else {
            boundbox3f nodeBounds = m_bounds;
            const P* result = 0;
            recursive_support( d, nodeBounds, 0, &result );
            if( result != 0 )
                return *result;
            else
                throw std::runtime_error(
                    "particle_kdtree::support: Unexpected null value returned from recursive_support." );
        }
    }
}

template <class P>
inline void frantic::particles::particle_kdtree<P>::recursive_support( const frantic::graphics::vector3f& d,
                                                                       frantic::graphics::boundbox3f& nodeBounds,
                                                                       int index, const P** outSupport ) const {
    using namespace frantic::graphics;

    const P& p = m_particles[index];
    int axis = m_axes[index];

    // If we've already got a candidate point, we could throw away this whole node
    // based on the support of the bounding box.
    if( *outSupport != 0 ) {
        vector3f boundsSupport = nodeBounds.support( d );
        // When the support of bounding box isn't in the direction of d from our current
        // candidate support point, there's no point in checking the points in that bounding box.
        // This test would provide the main speed boost over just testing all the points.
        if( vector3f::dot( d, boundsSupport - **outSupport ) <= 0 ) {
            return;
        }
    }

    if( 2 * index + 1 < (int)m_particles.size() ) // Check that there's a left child
    {
        float axisDistance = d[axis];

        if( axisDistance > 0 ) {
            // search right child first
            if( 2 * index + 2 < (int)m_particles.size() ) {
                const float temp = nodeBounds.minimum()[axis];
                nodeBounds.minimum()[axis] = p[axis];
                recursive_support( d, nodeBounds, 2 * index + 2, outSupport );
                nodeBounds.minimum()[axis] = temp;
            }
            const float temp = nodeBounds.maximum()[axis];
            nodeBounds.maximum()[axis] = p[axis];
            recursive_support( d, nodeBounds, 2 * index + 1, outSupport );
            nodeBounds.maximum()[axis] = temp;

        } else {
            // search left child first
            const float temp = nodeBounds.maximum()[axis];
            nodeBounds.maximum()[axis] = p[axis];
            recursive_support( d, nodeBounds, 2 * index + 1, outSupport );
            nodeBounds.maximum()[axis] = temp;

            if( 2 * index + 2 < (int)m_particles.size() ) {
                const float temp = nodeBounds.minimum()[axis];
                nodeBounds.minimum()[axis] = p[axis];
                recursive_support( d, nodeBounds, 2 * index + 2, outSupport );
                nodeBounds.minimum()[axis] = temp;
            }
        }
    }

    // Process the current particle last.  This is because the nodes around the edge
    // are more likely to be leaf nodes, so skipping the processing of parent nodes
    // until later might save some time.
    if( *outSupport == 0 ) {
        // If no point selected yet, use this one
        *outSupport = &p;
    } else {
        // If a candidate point was already found, check whether this point beats it
        if( vector3f::dot( d, p - **outSupport ) > 0 )
            *outSupport = &p;
    }
}
