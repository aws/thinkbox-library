// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <boost/shared_ptr.hpp>

#include <frantic/geometry/raytracing.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/motion_blurred_transform.hpp>

namespace frantic {
namespace geometry {
using frantic::graphics::ray3f;

template <class PrimitiveType>
class primitive_kdtree;

// This is a node class for primitive_kdtree.  This shouldn't be created independently of that class
template <class PrimitiveType>
class primitive_kdtree_node {
    // Everything is private, because only the declared friend classes need access.
    friend class primitive_kdtree<PrimitiveType>;

    bool find_nearest_point( const vector3f& point, float maxDistance, frantic::graphics::boundbox3f& bounds,
                             const std::vector<boost::shared_ptr<PrimitiveType>> primitives,
                             nearest_point_search_result& outNearestPoint ) const {
        if( m_axisAndLeafFlag < 0 ) {
            nearest_point_search_result testNearestPoint;
            double bestDistance = std::numeric_limits<double>::infinity();

            bool found = false;

            for( unsigned i = 0; i < primitives.size(); ++i ) {
                primitives[i]->find_nearest_point( point, maxDistance, testNearestPoint );

                if( testNearestPoint.distance < bestDistance && testNearestPoint.distance <= maxDistance ) {
                    bestDistance = testNearestPoint.distance;

                    outNearestPoint = testNearestPoint;

                    outNearestPoint.position = testNearestPoint.position;

                    outNearestPoint.primitiveIndex = i;

                    found = true;
                }
            }

            return found;
        } else {
            if( vector3f::distance_squared( bounds.clamp_nothrow( point ), point ) >= maxDistance * maxDistance )
                return false;

            int axis = m_axisAndLeafFlag & 0x3;

            if( point[axis] < m_split ) {
                float savedBoundsValue = bounds.maximum()[axis];
                bounds.maximum()[axis] = m_split;

                // First find the nearest point in the left child
                bool found =
                    left_child()->find_nearest_point( point, maxDistance, bounds, primitives, outNearestPoint );
                bounds.maximum()[axis] = savedBoundsValue;

                // If found, constrict the search radius for the search in the right child
                if( found )
                    maxDistance = (float)outNearestPoint.distance;

                savedBoundsValue = bounds.minimum()[axis];
                bounds.minimum()[axis] = m_split;
                // Then find the nearest point in the right child, using this new constrained radius
                if( right_child()->find_nearest_point( point, maxDistance, bounds, primitives, outNearestPoint ) )
                    found = true;
                bounds.minimum()[axis] = savedBoundsValue;

                return found;
            } else if( point[axis] > m_split ) {
                float savedBoundsValue = bounds.minimum()[axis];
                bounds.minimum()[axis] = m_split;
                // First find the nearest point in the right child
                bool found =
                    right_child()->find_nearest_point( point, maxDistance, bounds, primitives, outNearestPoint );
                bounds.minimum()[axis] = savedBoundsValue;

                // If found, constrict the search radius for the search in the left child
                if( found )
                    maxDistance = (float)outNearestPoint.distance;

                savedBoundsValue = bounds.maximum()[axis];
                bounds.maximum()[axis] = m_split;
                // Then find the nearest point in the left child, using this new constrained radius
                if( left_child()->find_nearest_point( point, maxDistance, bounds, primitives, outNearestPoint ) )
                    found = true;
                bounds.maximum()[axis] = savedBoundsValue;

                return found;
            }
        }

        return false;
    }

    bool is_point_in_volume( const vector3f& point,
                             const std::vector<boost::shared_ptr<PrimitiveType>> primitives ) const {
        if( m_axisAndLeafFlag < 0 ) {
            for( unsigned i = 0; i < primitives.size(); ++i ) {
                if( primitives[i]->is_point_in_volume( point ) ) {
                    return true;
                }
            }

            return false;
        } else {
            int axis = m_axisAndLeafFlag & 0x3;

            // double distance = m_split - point[axis];

            if( point[axis] < m_split ) {
                return right_child()->is_point_in_volume( point, primitives );
            } else if( point[axis] > m_split ) {
                return left_child()->is_point_in_volume( point, primitives );
            } else {
                // try both sides if the point is in the splitting plane
                return left_child()->is_point_in_volume( point, primitives ) ||
                       right_child()->is_point_in_volume( point, primitives );
            }
        }
    }

    union {
        // Interior node property
        primitive_kdtree_node*
            m_children; // The children are stored as an array of two nodes, so delete[] must be used.
        // Leaf node property
        std::vector<int>* m_primitiveIndices;
    };
    // Bottom 2 bits determine axis, sign bit determines whether the node is a leaf or not
    int m_axisAndLeafFlag;
    float m_split;

    primitive_kdtree_node() {
        // Initialize it indicating it's not a leaf node.
        m_axisAndLeafFlag = 0;
        m_children = 0;
    }

    // Constructor for making an interior node
    // NOTE: Only one initialize function should be called once on a node.  No checks are done to make sure things are
    // valid.
    void initialize( primitive_kdtree_node* children, int axis, float split ) {
        m_split = split;
        m_children = children;
        // Set the axis and indicate it's not a leaf
        m_axisAndLeafFlag = axis | 0x00000000;
    }

    // Constructor for making a leaf node
    // This destructively swaps out the primitive indices from the vector that gets passed in, to avoid the extra copy.
    // NOTE: Only one initialize function should be called once on a node.  No checks are done to make sure things are
    // valid.
    void initialize( std::vector<int>& primitiveIndices ) {
        if( primitiveIndices.size() > 0 ) {
            m_primitiveIndices = new std::vector<int>();
            m_primitiveIndices->swap( primitiveIndices );
        } else {
            m_primitiveIndices = 0;
        }

        // Indicate it's a leaf
        m_axisAndLeafFlag = 0x80000000;
    }

    ~primitive_kdtree_node() {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices != 0 ) {
                delete m_primitiveIndices;
                m_primitiveIndices = 0;
            }
        } else {
            if( m_children != 0 ) {
                delete[] m_children;
                m_children = 0;
            }
        }
    }

    primitive_kdtree_node* left_child() { return m_children; }
    primitive_kdtree_node* right_child() { return m_children + 1; }

    const primitive_kdtree_node* left_child() const { return m_children; }
    const primitive_kdtree_node* right_child() const { return m_children + 1; }

    int get_node_count() const {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            return 1;
        } else {
            return 1 + left_child()->get_node_count() + right_child()->get_node_count();
        }
    }

    int get_maximum_depth() const {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            return 0;
        } else {
            return 1 + ( std::max )( left_child()->get_maximum_depth(), right_child()->get_maximum_depth() );
        }
    }

    int get_largest_leaf_size() const {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices != 0 )
                return (int)m_primitiveIndices->size();
            else
                return 0;
        } else {
            return ( std::max )( left_child()->get_largest_leaf_size(), right_child()->get_largest_leaf_size() );
        }
    }

    void dump_tree( std::ostream& out, int depth,
                    const std::vector<boost::shared_ptr<PrimitiveType>>& primitives ) const {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices != 0 ) {
                out << depth << ": child node with " << m_primitiveIndices->size() << " children\n";
                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    out << depth << ":  primitive " << ( *m_primitiveIndices )[i] << ", box "
                        << primitives[( *m_primitiveIndices )[i]]->get_bounds() << "\n";
                    // boost::shared_ptr<PrimitiveType> primitive = primitives[(*m_faceIndices)[i]];
                    // out << depth << ": face " << mesh.get_vertex(face.x) << " " << mesh.get_vertex(face.y) << " " <<
                    // mesh.get_vertex(face.z) << "\n";
                }
            } else {
                out << depth << ": empty leaf node\n";
            }
        } else {
            char dimensions[] = { 'X', 'Y', 'Z' };
            out << depth << ": interior node " << dimensions[m_axisAndLeafFlag & 0x00000003] << " " << m_split << "\n";
            out << depth << ": left child\n";
            left_child()->dump_tree( out, depth + 1, primitives );
            out << depth << ": right child\n";
            right_child()->dump_tree( out, depth + 1, primitives );
        }
    }

    // This returns the first intersection between the ray and the mesh, or false if none is found.
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                        const std::vector<boost::shared_ptr<PrimitiveType>>& primitives,
                        frantic::geometry::raytrace_intersection& outIntersection ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.

            // Empty nodes are optimized by not allocating the primitive indices pointer
            if( m_primitiveIndices != 0 ) {
                bool found = false; // m_primitiveIndices->nosuchthing().12 = 47;

                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    int primitiveIndex = ( *m_primitiveIndices )[i];

                    if( primitives[primitiveIndex]->intersect_ray( ray, tMin, tMax, outIntersection ) ) {
                        outIntersection.primitiveIndex = primitiveIndex;

                        tMax = outIntersection.distance;
                        found = true;
                    }
                }
                return found;
            } else {
                return false;
            }
        } else {
            // It's an interior node.
            int axis = m_axisAndLeafFlag & 0x00000003;
            double distance = ( (double)m_split - ray.origin()[axis] ) / ray.direction()[axis];

            // if the ray is parallel to the split axis, we only test one child
            if( ray.direction()[axis] == 0 ) {
                if( ray.origin()[axis] <= m_split )
                    return left_child()->intersect_ray( ray, tMin, tMax, primitives, outIntersection );
                else
                    return right_child()->intersect_ray( ray, tMin, tMax, primitives, outIntersection );
            } else {

                // Traverse the children in the order so the closest child is tried first
                bool result = false;
                if( ray.direction()[axis] > 0 ) {

                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the right side
                        return right_child()->intersect_ray( ray, tMin, tMax, primitives, outIntersection );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the left side
                        return left_child()->intersect_ray( ray, tMin, tMax, primitives, outIntersection );
                    } else {
                        result = left_child()->intersect_ray( ray, tMin, distance, primitives, outIntersection );
                        if( !result )
                            result = right_child()->intersect_ray( ray, distance, tMax, primitives, outIntersection );
                    }
                } else {
                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the left side
                        return left_child()->intersect_ray( ray, tMin, tMax, primitives, outIntersection );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the right side
                        return right_child()->intersect_ray( ray, tMin, tMax, primitives, outIntersection );
                    } else {
                        result = right_child()->intersect_ray( ray, tMin, distance, primitives, outIntersection );
                        if( !result )
                            result = left_child()->intersect_ray( ray, distance, tMax, primitives, outIntersection );
                    }
                }

                return result;
            }
        }
    }

    // This returns an array of all the intersections between the ray and the mesh.
    void intersect_ray_all( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                            const std::vector<boost::shared_ptr<PrimitiveType>>& primitives,
                            std::vector<frantic::geometry::raytrace_intersection>& outIntersections ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.

            // Empty nodes are optimized by not allocating the primitive indices pointer
            if( m_primitiveIndices != 0 ) {
                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    int primitiveIndex = ( *m_primitiveIndices )[i];

                    primitives[primitiveIndex]->intersect_ray_all( ray, tMin, tMax, outIntersections );
                }
            }
        } else {
            // It's an interior node
            int axis = m_axisAndLeafFlag & 0x00000003;

            double distance = ( (double)m_split - ray.origin()[axis] ) / ray.direction()[axis];

            // Traverse the children in the order so the closest child is processed first,
            // thus putting all the intersections in ascending order
            if( ray.direction()[axis] > 0 ) {
                if( distance < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the right side
                    right_child()->intersect_ray_all( ray, tMin, tMax, primitives, outIntersections );
                } else if( distance > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the left side
                    left_child()->intersect_ray_all( ray, tMin, tMax, primitives, outIntersections );
                } else {
                    left_child()->intersect_ray_all( ray, tMin, distance, primitives, outIntersections );
                    right_child()->intersect_ray_all( ray, distance, tMax, primitives, outIntersections );
                }
            } else {
                if( distance < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the right side
                    left_child()->intersect_ray_all( ray, tMin, tMax, primitives, outIntersections );
                } else if( distance > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the left side
                    right_child()->intersect_ray_all( ray, tMin, tMax, primitives, outIntersections );
                } else {
                    right_child()->intersect_ray_all( ray, tMin, distance, primitives, outIntersections );
                    left_child()->intersect_ray_all( ray, distance, tMax, primitives, outIntersections );
                }
            }
        }
    }

    // This returns whether or not the ray intersects with the mesh within the given distance.
    bool intersects_ray_segment( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                 const std::vector<boost::shared_ptr<PrimitiveType>>& primitives ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.

            if( m_primitiveIndices != 0 ) {
                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    int primitiveIndex = ( *m_primitiveIndices )[i];

                    if( primitives[primitiveIndex]->intersects_ray_segment( ray, tMin, tMax ) )
                        return true;
                }
            }

            return false;
        } else {
            // It's an interior node
            int axis = m_axisAndLeafFlag & 0x00000003;

            double distance = ( (double)m_split - ray.origin()[axis] ) / ray.direction()[axis];

            // Traverse the children in the order so the closest child is tried first
            bool result = false;
            if( ray.direction()[axis] > 0 ) {
                if( distance < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the right side
                    return right_child()->intersects_ray_segment( ray, tMin, tMax, primitives );
                } else if( distance > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the left side
                    return left_child()->intersects_ray_segment( ray, tMin, tMax, primitives );
                } else {
                    result = left_child()->intersects_ray_segment( ray, tMin, distance, primitives );
                    if( !result )
                        result = right_child()->intersects_ray_segment( ray, distance, tMax, primitives );
                }
            } else {
                if( distance < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the left side
                    return left_child()->intersects_ray_segment( ray, tMin, tMax, primitives );
                } else if( distance > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the right side
                    return right_child()->intersects_ray_segment( ray, tMin, tMax, primitives );
                } else {
                    result = right_child()->intersects_ray_segment( ray, tMin, distance, primitives );
                    if( !result )
                        result = left_child()->intersects_ray_segment( ray, distance, tMax, primitives );
                }
            }
            return result;
        }
    }
};

} // namespace geometry
} // namespace frantic
