// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

// THIS HEADER SHOULD ONLY BE INCLUDED AT THE END OF generic_mesh_kdtree.hpp

#include <boost/foreach.hpp>

#include <frantic/geometry/mixed_kdtree.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/utility.hpp>

namespace frantic {
namespace geometry {

struct generic_mesh_kdtree_node {
    enum { flag_leaf = ( 1 << 31 ), mask_axis = 0x03 };

    union {
        // Interior node property
        generic_mesh_kdtree_node*
            m_children; // The children are stored as an array of two nodes, so delete[] must be used.
        // Leaf node property
        unsigned* m_primitiveIndices; // First entry is num primitives, the rest are primitive indices.
    };

    // Bottom 2 bits determine axis, sign bit determines whether the node is a leaf or not
    unsigned m_axisAndLeafFlag;
    float m_split;

    inline bool is_leaf() const { return ( m_axisAndLeafFlag & flag_leaf ) != 0; }
    inline bool is_internal() const { return ( m_axisAndLeafFlag & flag_leaf ) == 0; }

    inline unsigned get_axis() const { return ( m_axisAndLeafFlag & mask_axis ); }
    inline float get_split() const { return m_split; }
    inline const generic_mesh_kdtree_node& get_left_child() const { return m_children[0]; }
    inline const generic_mesh_kdtree_node& get_right_child() const { return m_children[1]; }
    inline generic_mesh_kdtree_node& get_left_child() { return m_children[0]; }
    inline generic_mesh_kdtree_node& get_right_child() { return m_children[1]; }

    // Only call this if we are sure there are actually primitives!
    inline bool is_empty() const { return ( m_primitiveIndices == NULL ); }
    inline unsigned* begin_primitives() const { return m_primitiveIndices + 1; }
    inline unsigned* end_primitives() const { return m_primitiveIndices + m_primitiveIndices[0] + 1; }

    inline void initialize( int axis, float split ) {
        m_split = split;
        m_axisAndLeafFlag = axis;
        m_children = new generic_mesh_kdtree_node[2];
    }

    template <class IndexContainer>
    inline void initialize( const IndexContainer& indices ) {
        m_axisAndLeafFlag = 0u | flag_leaf;
        if( indices.empty() ) {
            m_primitiveIndices = NULL;
        } else {
            m_primitiveIndices = new unsigned[indices.size() + 1];
            m_primitiveIndices[0] = (unsigned)indices.size();

            std::copy( indices.begin(), indices.end(), m_primitiveIndices + 1 );
        }
    }

    ~generic_mesh_kdtree_node() {
        if( is_leaf() )
            delete[] m_primitiveIndices;
        else
            delete[] m_children;
    }
};

namespace detail {
template <class NodeType, class MeshTraits>
void build_kdtree_greedy_SAH_nlogn( NodeType& rootNode, const typename MeshTraits::mesh_type& mesh,
                                    const frantic::graphics::boundbox3f& bounds );
}

template <class MeshTraits>
generic_mesh_kdtree<MeshTraits>::generic_mesh_kdtree()
    : m_mesh( NULL ) {}

template <class MeshTraits>
generic_mesh_kdtree<MeshTraits>::generic_mesh_kdtree( typename MeshTraits::mesh_type_const_ptr mesh, bool final ) {
    set_mesh( mesh );
    if( final )
        finalize();
}

template <class MeshTraits>
generic_mesh_kdtree<MeshTraits>::~generic_mesh_kdtree() {}

template <class MeshTraits>
void generic_mesh_kdtree<MeshTraits>::set_mesh( typename MeshTraits::mesh_type_const_ptr mesh ) {
    m_mesh = mesh;
}

template <class MeshTraits>
void generic_mesh_kdtree<MeshTraits>::finalize() {
    m_rootBounds = MeshTraits::get_bounds( *m_mesh );

    // Expand the box somewhat to avoid clipping any triangles lying on the outer edges.
    m_rootBounds.expand_fractional( 0.00001f );

    m_rootNode.reset( new node_type );

    detail::build_kdtree_greedy_SAH_nlogn<node_type, MeshTraits>( *m_rootNode, *m_mesh, m_rootBounds );
}

template <class MeshTraits>
bool generic_mesh_kdtree<MeshTraits>::intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                                     raytrace_result& outIntersection, bool ignoreBackfaces ) const {
    if( !ray.clamp_to_box( m_rootBounds, tMin, tMax ) )
        return false;

    outIntersection.distance = tMax;

    return intersect_ray_recursive( *m_rootNode, ray, tMin, tMax, outIntersection, ignoreBackfaces );
}

template <class MeshTraits>
bool generic_mesh_kdtree<MeshTraits>::intersect_ray_recursive( const node_type& curNode,
                                                               const frantic::graphics::ray3f& ray, double tMin,
                                                               double tMax, raytrace_result& outIntersection,
                                                               bool ignoreBackfaces ) const {
    if( curNode.is_leaf() ) {
        if( curNode.is_empty() )
            return false;

        bool result = false;

        for( unsigned *it = curNode.begin_primitives(), *itEnd = curNode.end_primitives(); it != itEnd; ++it ) {
            if( MeshTraits::intersect_ray( *m_mesh, *it, ray, tMin, tMax, outIntersection, ignoreBackfaces ) ) {
                result = true;
                tMax = outIntersection.distance;
            }
        }

        return result;
    } else {
        unsigned axis = curNode.get_axis();
        double rayOriginAlongAxis = (double)ray.origin()[axis];
        double rayDirectionAlongAxis = (double)ray.direction()[axis];

        if( rayDirectionAlongAxis == 0 ) {
            // When the ray is parallel, we only need to check one of the sides
            if( rayOriginAlongAxis <= curNode.get_split() ) {
                return intersect_ray_recursive( curNode.get_left_child(), ray, tMin, tMax, outIntersection,
                                                ignoreBackfaces );
            } else {
                return intersect_ray_recursive( curNode.get_right_child(), ray, tMin, tMax, outIntersection,
                                                ignoreBackfaces );
            }
        } else {
            double distance = ( (double)curNode.get_split() - rayOriginAlongAxis ) / rayDirectionAlongAxis;

            // HACK: Due to rounding error in the intersection code, we need to adjust the splitting plane location to
            // compensate for triangles
            //       that lie directly on the plane.

            // HACK: This epsilon is somehwat arbitrary, but it should capture a relative error measurement that can be
            // used
            //       to compensate for inaccuracy when calculating the ray distance to triangles who lie directly on
            //       the splitting plane. That calculation is more error prone than this

            // HACK: Shouldn't this actually be in the intersection code, where it might accept an answer slightly
            // outside
            //       of the range passed in? We might still want an error measure here too...
            double eps = distance * 16.0 * std::numeric_limits<float>::epsilon();

            std::pair<double, double> interval( distance - eps, distance + eps );

            // Traverse the children in the order so the closest child is tried first
            if( rayDirectionAlongAxis > 0 ) {
                if( interval.first > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the left side
                    return intersect_ray_recursive( curNode.get_left_child(), ray, tMin, tMax, outIntersection,
                                                    ignoreBackfaces );
                } else if( interval.second < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the right side
                    return intersect_ray_recursive( curNode.get_right_child(), ray, tMin, tMax, outIntersection,
                                                    ignoreBackfaces );
                } else {
                    bool result = intersect_ray_recursive( curNode.get_left_child(), ray, tMin, interval.second,
                                                           outIntersection, ignoreBackfaces );
                    if( !result )
                        result = intersect_ray_recursive( curNode.get_right_child(), ray, interval.first, tMax,
                                                          outIntersection, ignoreBackfaces );
                    else if( outIntersection.distance > interval.first )
                        // If the intersection was found in the region of uncertainty, search the other side too. We may
                        // have found a triangle that overlaps both regions and is not actually the closest.
                        intersect_ray_recursive( curNode.get_right_child(), ray, interval.first,
                                                 outIntersection.distance, outIntersection );
                    return result;
                }
            } else {
                if( interval.second < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the left side
                    return intersect_ray_recursive( curNode.get_left_child(), ray, tMin, tMax, outIntersection,
                                                    ignoreBackfaces );
                } else if( interval.first > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the right side
                    return intersect_ray_recursive( curNode.get_right_child(), ray, tMin, tMax, outIntersection,
                                                    ignoreBackfaces );
                } else {
                    bool result = intersect_ray_recursive( curNode.get_right_child(), ray, tMin, interval.second,
                                                           outIntersection, ignoreBackfaces );
                    if( !result )
                        result = intersect_ray_recursive( curNode.get_left_child(), ray, interval.first, tMax,
                                                          outIntersection, ignoreBackfaces );
                    else if( outIntersection.distance > interval.first )
                        // If the intersection was found in the region of uncertainty, search the other side too. We may
                        // have found a triangle that overlaps both regions and is not actually the closest.
                        intersect_ray_recursive( curNode.get_left_child(), ray, interval.first,
                                                 outIntersection.distance, outIntersection );
                    return result;
                }
            }
        }
    }
}

template <class MeshTraits>
bool generic_mesh_kdtree<MeshTraits>::find_nearest_point( const frantic::graphics::vector3f& p, double tMax,
                                                          nearest_point_result& outResult,
                                                          bool ignoreBackfaces ) const {
    if( !m_rootNode )
        return false;

    ////Compute axis-wise squared distances to one corner (I arbitrarily chose the maximum one).
    // double cornerDist[] = {
    //	frantic::math::square( static_cast<double>( m_rootBounds.maximum().x - p.x ) ),
    //	frantic::math::square( static_cast<double>( m_rootBounds.maximum().y - p.y ) ),
    //	frantic::math::square( static_cast<double>( m_rootBounds.maximum().z - p.z ) )
    // };
    //
    // double oppDist;
    //
    ////Determine if the opposite side of the box is further, and grab it if so.
    // oppDist = frantic::math::square( static_cast<double>( p.x - m_rootBounds.minimum().x ) );
    // if( oppDist > cornerDist[0] )
    //	cornerDist[0] = oppDist;

    // oppDist = frantic::math::square( static_cast<double>( p.y - m_rootBounds.minimum().y ) );
    // if( oppDist > cornerDist[1] )
    //	cornerDist[1] = oppDist;

    // oppDist = frantic::math::square( static_cast<double>( p.z - m_rootBounds.minimum().z ) );
    // if( oppDist > cornerDist[2] )
    //	cornerDist[2] = oppDist;
    //
    ////Compute the distance to the furthest corner of the bounding box. Use that as tMax if closer.
    // double dist = std::sqrt( cornerDist[0] + cornerDist[1] + cornerDist[2] );
    // if( dist < tMax )
    //	tMax = dist;

    frantic::graphics::boundbox3f bounds = m_rootBounds;

    outResult.distance = tMax;

    return find_nearest_point_recursive( *m_rootNode.get(), p, tMax, bounds, outResult, ignoreBackfaces );
}

template <class MeshTraits>
bool generic_mesh_kdtree<MeshTraits>::find_nearest_point_recursive( const node_type& curNode,
                                                                    const frantic::graphics::vector3f& p, double tMax,
                                                                    frantic::graphics::boundbox3f& currentBounds,
                                                                    nearest_point_result& outResult,
                                                                    bool ignoreBackfaces ) const {
    if( curNode.is_leaf() ) {
        if( curNode.is_empty() || vector3f::distance_squared( currentBounds.clamp_nothrow( p ), p ) >= tMax * tMax )
            return false;

        bool result = false;

        for( unsigned *it = curNode.begin_primitives(), *itEnd = curNode.end_primitives(); it != itEnd; ++it ) {
            if( MeshTraits::find_nearest_point( *m_mesh, *it, p, tMax, outResult, ignoreBackfaces ) ) {
                result = true;
                tMax = outResult.distance;
            }
        }

        return result;
    } else {
        unsigned axis = curNode.get_axis();
        // double axisDistance = curNode.get_split() - p[axis];

        // If we're too far away from the bounding box, don't even bother looking.  Use clamp_nothrow to avoid
        // the is_empty check, for efficiency.
        if( frantic::graphics::vector3f::distance_squared( currentBounds.clamp_nothrow( p ), p ) >= tMax * tMax )
            return false;

        // Traverse the children so the closest child is tried first
        if( p[axis] < curNode.get_split() ) {
            float savedBoundsValue = currentBounds.maximum()[axis];
            currentBounds.maximum()[axis] = curNode.get_split();
            // First find the nearest point in the left child
            bool found = find_nearest_point_recursive( curNode.get_left_child(), p, tMax, currentBounds, outResult,
                                                       ignoreBackfaces );
            currentBounds.maximum()[axis] = savedBoundsValue;

            // If found, constrict the search radius for the search in the right child
            if( found )
                tMax = outResult.distance;

            savedBoundsValue = currentBounds.minimum()[axis];
            currentBounds.minimum()[axis] = curNode.get_split();
            // Then find the nearest point in the right child, using this new constrained radius
            if( find_nearest_point_recursive( curNode.get_right_child(), p, tMax, currentBounds, outResult,
                                              ignoreBackfaces ) )
                found = true;
            currentBounds.minimum()[axis] = savedBoundsValue;

            return found;
        } else {
            float savedBoundsValue = currentBounds.minimum()[axis];
            currentBounds.minimum()[axis] = curNode.get_split();
            // First find the nearest point in the right child
            bool found = find_nearest_point_recursive( curNode.get_right_child(), p, tMax, currentBounds, outResult,
                                                       ignoreBackfaces );
            currentBounds.minimum()[axis] = savedBoundsValue;

            // If found, constrict the search radius for the search in the left child
            if( found )
                tMax = outResult.distance;

            savedBoundsValue = currentBounds.maximum()[axis];
            currentBounds.maximum()[axis] = curNode.get_split();
            // Then find the nearest point in the left child, using this new constrained radius
            if( find_nearest_point_recursive( curNode.get_left_child(), p, tMax, currentBounds, outResult,
                                              ignoreBackfaces ) )
                found = true;
            currentBounds.maximum()[axis] = savedBoundsValue;

            return found;
        }

        //// Traverse the children so the closest child is tried first
        // if( axisDistance > 0 ) {
        //	float savedBoundsValue = currentBounds.maximum()[axis];
        //
        //	currentBounds.maximum()[axis] = curNode.get_split();

        //	// First find the nearest point in the left child
        //	bool result = find_nearest_point_recursive( curNode.get_left_child(), p, tMax, currentBounds, outResult
        //);

        //	//We can modify the max distance if we found something on the left side.
        //	if( result )
        //		tMax = outResult.distance;

        //	// Reset the modified boundbox.
        //	currentBounds.maximum()[axis] = savedBoundsValue;

        //	//If the right side could be closer, check it too.
        //	if( axisDistance < tMax ){
        //		savedBoundsValue = currentBounds.minimum()[axis];

        //		currentBounds.minimum()[axis] = curNode.get_split();

        //		result = find_nearest_point_recursive( curNode.get_right_child(), p, tMax, currentBounds,
        // outResult ) || result;

        //		currentBounds.minimum()[axis] = savedBoundsValue;
        //	}
        //
        //	return result;
        //} else {
        //	float savedBoundsValue = currentBounds.minimum()[axis];
        //
        //	currentBounds.minimum()[axis] = curNode.get_split();

        //	// First find the nearest point in the left child
        //	bool result = find_nearest_point_recursive( curNode.get_right_child(), p, tMax, currentBounds, outResult
        //);

        //	//We can modify the max distance if we found something on the left side.
        //	if( result )
        //		tMax = outResult.distance;

        //	// Reset the modified boundbox.
        //	currentBounds.minimum()[axis] = savedBoundsValue;

        //	//If the right side could be closer, check it too.
        //	if( -axisDistance < tMax ){
        //		savedBoundsValue = currentBounds.maximum()[axis];

        //		currentBounds.maximum()[axis] = curNode.get_split();

        //		result = find_nearest_point_recursive( curNode.get_right_child(), p, tMax, currentBounds,
        // outResult ) || result;

        //		currentBounds.maximum()[axis] = savedBoundsValue;
        //	}
        //
        //	return result;
        //}
    }
}

// template <class MeshTraits>
// bool generic_mesh_kdtree<MeshTraits>::find_nearest_point_recursive( const node_type& curNode, const
// frantic::graphics::vector3f& p, double tMax, nearest_point_result& outResult ) const { 	if( curNode.is_leaf() ){
// if(
// curNode.is_empty() ) 			return false;
//
//		bool result = false;
//
//		for( unsigned *it = curNode.begin_primitives(), *itEnd = curNode.end_primitives(); it != itEnd; ++it ){
//			if( MeshTraits::find_nearest_point( *m_mesh, *it, p, tMax, outResult ) ){
//				result = true;
//				tMax = outIntersection.distance;
//			}
//		}
//
//		return result;
//	}else{
//		unsigned axis = curNode.get_axis();
//		double axisDistance = curNode.get_split() - p[axis];
//
//		if( axisDistance >= 0 ){
//			//We're in the left sub-space, so search there first. If we find something closer than tMax we get
//a 'true' result
//			//and update the max distance. Then we optionally search the other space if it could have
//something
// closer. 			bool result = find_nearest_point_recursive( curNode.get_left_child(), p, tMax, outResult );
// if( result ) 				tMax
//= outResult.distance; 			if( axisDistance < tMax ) 				result =
//find_nearest_point_recursive( curNode.get_right_child(), p,
// tMax, outResult ) || result; 			return result; 		}else{ 			bool result =
// find_nearest_point_recursive( curNode.get_right_child(), p, tMax, outResult ); 			if( result )
// tMax = outResult.distance; 			if( -axisDistance < tMax ) 				result = find_nearest_point_recursive(
//curNode.get_left_child(), p, tMax, outResult ) || result; 			return result;
//		}
//	}
// }

namespace detail {

// This simultaneously checks for an empty bounding box, and also confirms that the bounds are not set to indeterminant
// values which would otherwise cause the kdtree creation to enter an infinite loop.
inline bool is_boundbox_empty_or_invalid( const frantic::graphics::boundbox3f& bounds ) {
    return !( bounds.maximum().x >= bounds.minimum().x ) || !( bounds.maximum().y >= bounds.minimum().y ) ||
           !( bounds.maximum().y >= bounds.minimum().y );
}

// I totally stole this from Paul's implementation in mixed_kdtree.cpp
template <class NodeType, class MeshTraits>
void build_kdtree_greedy_SAH_nlogn( NodeType& theNode, const typename MeshTraits::mesh_type& mesh,
                                    const frantic::graphics::boundbox3f& bounds,
                                    std::vector<mixed_kdtree_detail::kdtreeEvent>& events,
                                    std::vector<mixed_kdtree_detail::index_t>& indices,
                                    std::vector<boost::int8_t>& objFlags, const int maximumDepth, const int depth ) {
    // Could use Wald's stopping condition, instead of stopping at 8
    // triangles.
    //
    // Wald's stopping condition is part of the initial bestCost below:
    // if we don't find anything better, then the recursion stops.
    // from Wald:
    // Terminate( triangles T, voxel V) =
    // {
    //   true if argmin_p C_V(p) > K_I|T|
    //   false otherwise
    // }
    // here we have T : indices( need anything else? ), V: bounds
    // that is: true if (bestCost > COST_INTERSECT * indices.size())
    // or set initial bestCost = COST_INTERSECT * indices.size()

    // if(indices.size() <= 8){
    // node.initialize(indices);
    // return;
    //}
    bool doneSplitting = false;

    const float voxelSA = bounds.get_surface_area();

    // TODO: I think this used to be done indirectly by the
    // event detection code, which counted an AABB as parallel
    // once it got smaller than a threshold.
    //
    // Since this is gone now, I'm instead checking whether
    // the clipped bounds are greater than some fraction
    // of their original extents.  Maybe this should be
    // depend on the average size of features in the scene,
    // or a function of the size of features in the current
    // box ?
    // if( bounds.get_max_dimension() < minimumVoxelLength ) {

    // I've switched to limiting the depth for now

    const std::size_t nMin = 1;

    if( depth > maximumDepth || indices.size() <= nMin ) {
        doneSplitting = true;
    }

    int bestAxis = -1;
    // float bestCost = 0.9f * static_cast<float>( mixed_kdtree_detail::COST_INTERSECT ) * indices.size();
    const float intersectionCost = static_cast<float>( mixed_kdtree_detail::COST_INTERSECT ) * indices.size();
    float bestCost = intersectionCost;
    float bestSplit = std::numeric_limits<float>::max();

    // The side that triangles on the splitting plane should go to.
    // This is \hat{p}_{side} in Wald's paper.
    mixed_kdtree_detail::kdtreeSide_type bestPlanarPartitionSide = mixed_kdtree_detail::LEFT;

    std::size_t nLeft[3], nRight[3];
    nLeft[0] = nLeft[1] = nLeft[2] = 0;
    nRight[0] = nRight[1] = nRight[2] = indices.size();

    if( !doneSplitting ) {
        for( std::size_t i = 0; i < events.size(); /*do nothing*/ ) {
            int axis = events[i].axis;
            float split = events[i].location;

            std::size_t counters[] = { 0, 0, 0 };
            for( ; i < events.size() && events[i].location == split && events[i].axis == axis; ++i )
                ++counters[events[i].type];

            nRight[axis] -= ( counters[mixed_kdtree_detail::END_] + counters[mixed_kdtree_detail::PARALLEL] );

            boundbox3f left( bounds ), right( bounds );
            left.maximum()[axis] = right.minimum()[axis] = split;
            const float leftSA = left.get_surface_area();
            const float rightSA = right.get_surface_area();

            // todo: I'm not sure if this is a reasonable criteria for
            // accepting splits:  each half must have either a primitive
            // in it, or it must have some volume.
            if( ( ( nLeft[axis] + counters[mixed_kdtree_detail::PARALLEL] ) > 0 || left.size( axis ) > 0 ) &&
                ( nRight[axis] > 0 || right.size( axis ) > 0 ) ) {
                float cost = mixed_kdtree_detail::SAH_cost(
                    voxelSA, leftSA, rightSA, nLeft[axis] + counters[mixed_kdtree_detail::PARALLEL], nRight[axis] );
                if( cost < bestCost ) {
                    bestCost = cost;
                    bestSplit = split;
                    bestAxis = axis;
                    bestPlanarPartitionSide = mixed_kdtree_detail::LEFT;
                }
            }
            // If there are no parallel triangles, then there is no
            // need to check for their cost on the right-hand side.
            if( counters[mixed_kdtree_detail::PARALLEL] > 0 ) {
                if( ( nLeft[axis] > 0 || left.size( axis ) > 0 ) &&
                    ( ( nRight[axis] + counters[mixed_kdtree_detail::PARALLEL] ) > 0 || right.size( axis ) > 0 ) ) {
                    float cost = mixed_kdtree_detail::SAH_cost(
                        voxelSA, leftSA, rightSA, nLeft[axis], nRight[axis] + counters[mixed_kdtree_detail::PARALLEL] );
                    if( cost < bestCost ) {
                        bestCost = cost;
                        bestSplit = split;
                        bestAxis = axis;
                        bestPlanarPartitionSide = mixed_kdtree_detail::RIGHT;
                    }
                }
            }

            nLeft[axis] += ( counters[mixed_kdtree_detail::BEGIN] + counters[mixed_kdtree_detail::PARALLEL] );
        }
    }

    // this count may be unnecessary
    // I added it to avoid reallocations.  This seemed beneficial
    // in test cases but it may not be in general.
    std::size_t nBestLeft = 0;
    std::size_t nBestRight = 0;

    // if( bestAxis != -1 && indices.size() > 2 && ! doneSplitting ){
    if( bestAxis != -1 && !doneSplitting ) {
        // std::cout << "splitting axis " << bestAxis << " at " << bestSplit << "\n";
        // We need to classify triangles as being on either side of the split, or crossing
        for( std::vector<mixed_kdtree_detail::index_t>::const_iterator index = indices.begin(), end = indices.end();
             index != end; ++index )
            objFlags[*index] = 0;

        for( std::vector<mixed_kdtree_detail::kdtreeEvent>::const_iterator it = events.begin(), end = events.end();
             it != end; ++it ) {
            if( it->axis != bestAxis )
                continue;

            if( it->location < bestSplit ) {
                if( it->type != mixed_kdtree_detail::BEGIN ) {
                    ++nBestLeft;
                    objFlags[it->index] = -1;
                }
            } else if( it->location > bestSplit ) {
                if( it->type != mixed_kdtree_detail::END_ ) {
                    ++nBestRight;
                    objFlags[it->index] = 1;
                }
            } else {
                if( it->type == mixed_kdtree_detail::END_ ) {
                    ++nBestLeft;
                    objFlags[it->index] = -1;
                } else if( it->type == mixed_kdtree_detail::BEGIN ) {
                    ++nBestRight;
                    objFlags[it->index] = 1;
                } else if( it->type == mixed_kdtree_detail::PARALLEL && it->location == bestSplit ) {
                    if( bestPlanarPartitionSide == mixed_kdtree_detail::LEFT ) {
                        ++nBestLeft;
                        objFlags[it->index] = -1;
                    } else {
                        ++nBestRight;
                        objFlags[it->index] = +1;
                    }
                }
            }
        }

        const std::size_t nBestIntersecting = indices.size() - nBestLeft - nBestRight;

        const std::size_t reserveLeftIndices = nBestIntersecting + nBestLeft;
        const std::size_t reserveRightIndices = nBestIntersecting + nBestRight;
        const std::size_t reserveLeftEvents0 = 6 * nBestLeft;
        const std::size_t reserveRightEvents0 = 6 * nBestRight;
        const std::size_t reserveLeftEvents1 = 6 * nBestIntersecting;
        const std::size_t reserveRightEvents1 = 6 * nBestIntersecting;

        std::vector<mixed_kdtree_detail::kdtreeEvent> leftEvents[3], rightEvents[3];

        // Most events should be split into one side xor the other
        // leftEvents[0].reserve(events.size()/2);
        // rightEvents[0].reserve(events.size()/2);
        leftEvents[0].reserve( reserveLeftEvents0 );
        rightEvents[0].reserve( reserveRightEvents0 );

        for( std::vector<mixed_kdtree_detail::kdtreeEvent>::const_iterator it = events.begin(), end = events.end();
             it != end; ++it ) {
            if( objFlags[it->index] < 0 )
                leftEvents[0].push_back( *it );
            else if( objFlags[it->index] > 0 )
                rightEvents[0].push_back( *it );
        }

        frantic::clear_with_swap( events );

        std::vector<mixed_kdtree_detail::index_t> leftIndices, rightIndices;
        // reserve more than this?
        // leftIndices.reserve(indices.size()/2);
        // rightIndices.reserve(indices.size()/2);
        leftIndices.reserve( reserveLeftIndices );
        rightIndices.reserve( reserveRightIndices );

        boundbox3f leftBounds( bounds ), rightBounds( bounds );
        leftBounds.maximum()[bestAxis] = rightBounds.minimum()[bestAxis] = bestSplit;

        leftEvents[1].reserve( reserveLeftEvents1 );
        rightEvents[1].reserve( reserveRightEvents1 );

        BOOST_FOREACH( mixed_kdtree_detail::index_t i, indices ) {
            if( objFlags[i] < 0 )
                leftIndices.push_back( i );
            else if( objFlags[i] > 0 )
                rightIndices.push_back( i );
            else {
                // const vector3& f = mesh.get_face(i);
                // const boundbox3f unclippedBounds( primitives[i]->get_bounds() );
                // const float unclippedVolume = unclippedBounds.get_volume();

                boundbox3f clippedBounds;
                // clippedBounds = primitives[i]->intersect_with( leftBounds );
                clippedBounds = MeshTraits::get_clipped_bounds( mesh, i, leftBounds );

                // if( leftBounds.intersect_with_triangle(mesh.get_vertex(f.x), mesh.get_vertex(f.y),
                // mesh.get_vertex(f.z), clippedBounds) ){
                if( !is_boundbox_empty_or_invalid( clippedBounds ) ) {
                    // const int maxAxis = clippedBounds.get_max_dimension_axis();
                    // if ( unclippedBounds.size( maxAxis ) < 16.f * clippedBounds.size( maxAxis ) ) {
                    leftIndices.push_back( i );
                    extract_events( leftEvents[1], i, clippedBounds );
                }

                // clippedBounds = primitives[i]->intersect_with( rightBounds );
                clippedBounds = MeshTraits::get_clipped_bounds( mesh, i, rightBounds );

                // if( rightBounds.intersect_with_triangle(mesh.get_vertex(f.x), mesh.get_vertex(f.y),
                // mesh.get_vertex(f.z), clippedBounds) ){
                if( !is_boundbox_empty_or_invalid( clippedBounds ) ) {
                    // const int maxAxis = clippedBounds.get_max_dimension_axis();
                    // if ( unclippedBounds.size( maxAxis ) < 16.f * clippedBounds.size( maxAxis ) ) {
                    rightIndices.push_back( i );
                    extract_events( rightEvents[1], i, clippedBounds );
                }
            }
        }

        // std::cout << "0: " << leftEvents[0].size() << "/" << leftEvents[0].capacity() << " " << rightEvents[0].size()
        // <<
        // "/" << rightEvents[0].capacity() << std::endl; std::cout << "1: " << leftEvents[1].size() << "/" <<
        // leftEvents[1].capacity() << " " << rightEvents[1].size() << "/" << rightEvents[1].capacity() << std::endl;
        // std::cout << "i: " << leftIndices.size() << "/" << leftIndices.capacity() << " " << rightIndices.size() <<
        // "/" << rightIndices.capacity() << std::endl;

        frantic::clear_with_swap( indices );

        theNode.initialize( bestAxis, bestSplit );
        // node.initialize( new mixed_kdtree_node[2], bestAxis, bestSplit );

        leftEvents[2].resize( leftEvents[0].size() + leftEvents[1].size() );
        std::sort( leftEvents[1].begin(), leftEvents[1].end() );
        std::merge( leftEvents[0].begin(), leftEvents[0].end(), leftEvents[1].begin(), leftEvents[1].end(),
                    leftEvents[2].begin() );

        // leftEvents[2] and leftIndices
        // will be deallocated inside the recursive function calls.
        // Clear the others now.
        for( int i = 0; i < 2; ++i ) {
            frantic::clear_with_swap( leftEvents[i] );
        }

        rightEvents[2].resize( rightEvents[0].size() + rightEvents[1].size() );
        std::sort( rightEvents[1].begin(), rightEvents[1].end() );
        std::merge( rightEvents[0].begin(), rightEvents[0].end(), rightEvents[1].begin(), rightEvents[1].end(),
                    rightEvents[2].begin() );

        // rightEvents[2] and rightIndices
        // will be deallocated inside the recursive function calls.
        // Clear the others now.
        for( int i = 0; i < 2; ++i ) {
            frantic::clear_with_swap( rightEvents[i] );
        }

        build_kdtree_greedy_SAH_nlogn<NodeType, MeshTraits>( theNode.get_left_child(), mesh, leftBounds, leftEvents[2],
                                                             leftIndices, objFlags, maximumDepth, depth + 1 );
        build_kdtree_greedy_SAH_nlogn<NodeType, MeshTraits>( theNode.get_right_child(), mesh, rightBounds,
                                                             rightEvents[2], rightIndices, objFlags, maximumDepth,
                                                             depth + 1 );
        // build_kdtree_greedy_SAH_nlogn( *node.left_child(), primitives, leftBounds, leftEvents[2], leftIndices,
        // objFlags, maximumDepth, depth + 1 ); build_kdtree_greedy_SAH_nlogn( *node.right_child(), primitives,
        // rightBounds, rightEvents[2], rightIndices, objFlags, maximumDepth, depth + 1 );
    } else {
        // node.initialize( indices, primitives );
        theNode.initialize( indices );
    }
}

// template<class mixed_kdtree_primitive>
// void build_kdtree_greedy_SAH_nlogn(mixed_kdtree_node<mixed_kdtree_primitive>& node,
// std::vector<boost::shared_ptr<mixed_kdtree_primitive> >& primitives, const boundbox3f& bounds) {
template <class NodeType, class MeshTraits>
void build_kdtree_greedy_SAH_nlogn( NodeType& rootNode, const typename MeshTraits::mesh_type& mesh,
                                    const frantic::graphics::boundbox3f& bounds ) {
    // const std::size_t primitiveCount = primitives.size();
    const unsigned primitiveCount = MeshTraits::get_count( mesh );

    if( primitiveCount > std::numeric_limits<mixed_kdtree_detail::index_t>::max() ) {
        throw std::runtime_error(
            "build_kdtree_greedy_SAH_nlogn Internal Error: the number of primitives (" +
            boost::lexical_cast<std::string>( primitiveCount ) + ") exceeds the maximum primitive index (" +
            boost::lexical_cast<std::string>( std::numeric_limits<mixed_kdtree_detail::index_t>::max() ) + ")." );
    }

    std::vector<mixed_kdtree_detail::index_t> indices;
    std::vector<boost::int8_t> objFlags( primitiveCount, 0 );
    std::vector<mixed_kdtree_detail::kdtreeEvent> events;

    indices.reserve( primitiveCount );
    events.reserve( 6 * primitiveCount );

    size_t invalidFaceCount = 0;

    for( unsigned i = 0; i < primitiveCount; ++i ) {
        boundbox3f clippedBounds;
        // clippedBounds = primitives[i]->intersect_with( bounds );
        clippedBounds = MeshTraits::get_clipped_bounds( mesh, i, bounds );

        if( !is_boundbox_empty_or_invalid( clippedBounds ) ) {
            indices.push_back( static_cast<mixed_kdtree_detail::index_t>( i ) );
            extract_events( events, static_cast<mixed_kdtree_detail::index_t>( i ), clippedBounds );
        } else {
            ++invalidFaceCount;
        }
    }

    if( invalidFaceCount > 0 )
        FF_LOG( debug ) << _T("There are a total of ") << invalidFaceCount << _T(" invalid faces in the mesh!");

    // not sure if this is a good idea
    // This is an attempt to avoid degenerately small voxels.
    // Maybe I should instead adapt this to what's in the current voxel,
    // or set relative to the size of features in the scene ?
    // that I need to do this may also be a symptom of a problem with
    // the splitting planes that I accept.
    // const float minimumVoxelLength = std::min<float>( 0.00001f, 10.f * std::numeric_limits<float>::epsilon() *
    // bounds.get_max_dimension() );

    // heuristic suggested in Havran's PhD thesis, "Heuristic Ray Shooting Algorithms", 2001
    const int maximumDepth = static_cast<int>( 1.2 * log( static_cast<double>( primitiveCount ) ) / log( 2.0 ) + 2.0 );

    std::sort( events.begin(), events.end() );

    build_kdtree_greedy_SAH_nlogn<NodeType, MeshTraits>( rootNode, mesh, bounds, events, indices, objFlags,
                                                         maximumDepth, 0 );
}
} // namespace detail

} // namespace geometry
} // namespace frantic
