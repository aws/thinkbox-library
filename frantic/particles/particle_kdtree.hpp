// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#ifndef FRANTIC_KDTREE_H
#define FRANTIC_KDTREE_H

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <queue>

namespace frantic {
namespace particles {

namespace detail {
template <class P>
class locate_particles_count_sorter {
  public:
    locate_particles_count_sorter() {}

    bool operator()( const std::pair<float, P>& lhs, const std::pair<float, P>& rhs ) { return lhs.first < rhs.first; }
    typedef std::priority_queue<std::pair<float, P>, std::vector<std::pair<float, P>>, locate_particles_count_sorter>
        priority_queue;
};
} // namespace detail

template <class P>
class particle_kdtree {
    typedef size_t size_type;

    bool m_KDTreeBuiltFlag;

    frantic::graphics::boundbox3f m_bounds;

    std::vector<P> m_particles;
    std::vector<char> m_axes;

  public:
    particle_kdtree();
    ~particle_kdtree();

    size_type size() const { return m_particles.size(); }
    void clear(); // Clears elems, and resets all BBox values

    const P& operator[]( size_type index ) const { return m_particles[index]; }
    // Resets the kd-tree built flag, and returns a writable reference
    P& get_particle_ref( size_type index ) {
        m_KDTreeBuiltFlag = false;
        return m_particles[index];
    }

    void set_particles( const std::vector<P>& particles ) {
        m_KDTreeBuiltFlag = false;
        m_particles = particles;
    }

    void set_particles_with_swap( std::vector<P>& particles ) {
        m_KDTreeBuiltFlag = false;
        m_particles.swap( particles );
    }

    // If this flag is false, both the kd-tree and the bounding box values aren't to be trusted.
    // If this flag is true, both are correct.
    bool is_kdtree_built() const { return m_KDTreeBuiltFlag; }

    // Adds a particle to the kd-tree.  This resets the 'built' flag to false.
    void add_particle( const P& particle );

    // Creates the KD tree structure.  This sets the 'built' flag to true, so
    // the algorithms that use it will work.
    void balance_kdtree();

    // Finds all the particles within the specified range (squared).
    // The pairs returned contain the distance squared and the particle found.
    void locate_particles_range( const frantic::graphics::vector3f& position, float rangeSquared,
                                 std::vector<std::pair<float, P>>& outLocatedParticles ) const;

    // Finds the nearest particleCount particles.
    // The pairs returned contain the distance squared and the particle found.
    // The closest particle will be the first in the resulting list.
    void locate_particles_count( const frantic::graphics::vector3f& position, int particleCount,
                                 std::vector<std::pair<float, P>>& outLocatedParticles ) const;

    // Locates the first particle within the specified range (squared).
    // Returns true as soon as a particle is found. If none are found, it will return false.
    bool proximity_test( const frantic::graphics::vector3f& position, float rangeSquared ) const;

    // Finds the single nearest particle to the specified position.
    const P& locate_closest_particle( const frantic::graphics::vector3f& position ) const;
    // BARF: This was moved here from particle_kdtree_implementation.hpp because of Visual Studio 6 limitations.
    template <class Predicate>
    const P& locate_closest_particle_if( const frantic::graphics::vector3f& position, const Predicate& func ) const {
        if( m_particles.size() > 0 ) {
            const P* closestParticle = &m_particles[0];
            float closestParticleDistance = 1e38f;
            recursive_locate_closest_particle_if( position, 0, closestParticle, closestParticleDistance, func );
            return *closestParticle;
        } else {
            throw std::runtime_error(
                "Tried to call particle_kdtree::locate_closest_particle_if on an empty kd-tree." );
        }
    }

    // This computes the support of the convex hull of the points in the direction d.  The
    // support is the (not necessarily unique) point in the convex hull which
    // is furthest along in the direction d.  That is,
    // (p in this) such that { dot(p,d) >= p' for all (p' in this) }.
    // Mathematically, the support would be the set of all such points, but we just
    // return one.
    const P& support( const frantic::graphics::vector3f& d ) const;

  private:
    void balance_segment( std::vector<P>& outKDTree, int index, int start, int end );
    void median_split( std::vector<P>& particles, int start, int end, int median, int axis );

    void recursive_locate_particles_range( const frantic::graphics::vector3f& position, float rangeSquared,
                                           std::vector<std::pair<float, P>>& particles, int index ) const;
    void recursive_locate_particles_count( const frantic::graphics::vector3f& position, int particleCount,
                                           typename detail::locate_particles_count_sorter<P>::priority_queue& particles,
                                           int index ) const;
    bool recursive_proximity_test( const frantic::graphics::vector3f& position, float rangeSquared, int index ) const;

    // BARF: This was moved here from particle_kdtree_implementation.hpp because of Visual Studio 6 limitations.
    template <class Predicate>
    void recursive_locate_closest_particle_if( const frantic::graphics::vector3f& position, int index,
                                               const P*& closestParticle, float& closestParticleDistanceSquared,
                                               const Predicate& func ) const {
        const P& p = m_particles[index];
        int axis = m_axes[index];

        float axisDistance = position[axis] - p[axis];

        if( 2 * index + 1 < (int)m_particles.size() ) // Check that there's a left child
        {
            if( axisDistance > 0 ) {
                // search right child first
                if( 2 * index + 2 < (int)m_particles.size() )
                    recursive_locate_closest_particle_if( position, 2 * index + 2, closestParticle,
                                                          closestParticleDistanceSquared, func );
                if( axisDistance * axisDistance < closestParticleDistanceSquared )
                    recursive_locate_closest_particle_if( position, 2 * index + 1, closestParticle,
                                                          closestParticleDistanceSquared, func );

            } else {
                // search left child first
                recursive_locate_closest_particle_if( position, 2 * index + 1, closestParticle,
                                                      closestParticleDistanceSquared, func );
                if( 2 * index + 2 < (int)m_particles.size() &&
                    axisDistance * axisDistance < closestParticleDistanceSquared )
                    recursive_locate_closest_particle_if( position, 2 * index + 2, closestParticle,
                                                          closestParticleDistanceSquared, func );
            }
        }

        // We do this test afterwards because the recursive procedure above will descend in the
        // tree to an approximation of the closest point, thus narrowing the search sooner.  Also
        // test the axis distance (computed earlier for choosing which subtree to descend first)
        // before getting the distance squared as a possible early culling.

        if( axisDistance * axisDistance < closestParticleDistanceSquared ) {
            // Compute squared distance between current particle and ParticlesInRange
            float distanceSquaredToParticle = frantic::graphics::vector3f::distance_squared( p, position );

            // If this particle is closer than our closest so far and it meets the predicate value, set this to our
            // closest
            if( distanceSquaredToParticle < closestParticleDistanceSquared && func( p ) ) {
                closestParticleDistanceSquared = distanceSquaredToParticle;
                closestParticle = &p;
            }
        }
    }

    void recursive_support( const frantic::graphics::vector3f& d, frantic::graphics::boundbox3f& nodeBounds, int index,
                            const P** outSupport ) const;
};

} // namespace particles
} // namespace frantic

#include <frantic/particles/particle_kdtree_implementation.hpp>

#endif // FRANTIC_KDTREE_H
