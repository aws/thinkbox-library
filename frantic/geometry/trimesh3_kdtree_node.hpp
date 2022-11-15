// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <boost/tuple/tuple.hpp>
#include <frantic/geometry/raytracing.hpp>
#include <frantic/geometry/triangle_utils.hpp>

namespace frantic {
namespace geometry {
using frantic::graphics::ray3f;

class trimesh3_kdtree;
class trimesh3_kdtree_ext;

namespace detail {
struct kdtreeEvent;
}

struct traced_ray_info { // Debug output for following a ray_trace
    ray3f ray;
    float tMin, tMax;           // Ray in this space
    float intersect;            // Valid if internal, distance along ray to plane
    std::vector<int> triangles; // Valid if leaf
    bool isInternal;            // indicates either leaf or interal
    int axis;                   // split axis for internal
    float split;                // split value for internal
    boundbox3f voxel;           // space we are considering
};

// This is the node class for trimesh3_kdtree.  This shouldn't be created independently of that class
class trimesh3_kdtree_node {
    // Everything is private, because only the declared friend classes need access.
    friend class trimesh3_kdtree;
    friend void build_kdtree_greedy_SAH_nlog2n( trimesh3_kdtree_node& node, const trimesh3& mesh,
                                                std::vector<int>& indices, const boundbox3f& bounds );
    friend void build_kdtree_greedy_SAH_nlogn( trimesh3_kdtree_node& node, const trimesh3& mesh,
                                               const boundbox3f& bounds );
    friend void build_kdtree_greedy_SAH_nlogn( trimesh3_kdtree_node& node, const trimesh3& mesh,
                                               const boundbox3f& bounds, std::vector<detail::kdtreeEvent>& events,
                                               std::vector<int>& indices, std::vector<int>& triFlags );

    friend class trimesh3_kdtree_ext;

    union {
        // Interior node property
        trimesh3_kdtree_node* m_children; // The children are stored as an array of two nodes, so delete[] must be used.
        // Leaf node property
        std::vector<int>* m_primitiveIndices;
    };
    // Bottom 2 bits determine axis, sign bit determines whether the node is a leaf or not
    int m_axisAndLeafFlag;
    float m_split;

    trimesh3_kdtree_node() {
        // Initialize it indicating it's not a leaf node.
        m_axisAndLeafFlag = 0;
        m_children = 0;
    }

    // Constructor for making an interior node
    // NOTE: Only one initialize function should be called once on a node.  No checks are done to make sure things are
    // valid.
    void initialize( trimesh3_kdtree_node* children, int axis, float split ) {
        m_split = split;
        m_children = children;
        // Set the axis and indicate it's not a leaf
        m_axisAndLeafFlag = axis | 0x00000000;
    }

    // Constructor for making a leaf node
    // This destructively swaps out the face indices from the vector that gets passed in, to avoid the extra copy.
    // NOTE: Only one initialize function should be called once on a node.  No checks are done to make sure things are
    // valid.
    void initialize( std::vector<int>& faceIndices ) {
        if( faceIndices.size() > 0 ) {
            m_primitiveIndices = new std::vector<int>();
            m_primitiveIndices->swap( faceIndices );
        } else {
            m_primitiveIndices = 0;
        }

        // Indicate it's a leaf
        m_axisAndLeafFlag = 0x80000000;
    }

    ~trimesh3_kdtree_node() {
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

    trimesh3_kdtree_node* left_child() { return m_children; }
    trimesh3_kdtree_node* right_child() { return m_children + 1; }

    const trimesh3_kdtree_node* left_child() const { return m_children; }
    const trimesh3_kdtree_node* right_child() const { return m_children + 1; }

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

    void dump_tree( std::ostream& out, int depth, const trimesh3& mesh ) const {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices != 0 ) {
                out << depth << ": child node with " << m_primitiveIndices->size() << " children\n";
                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    vector3 face = mesh.get_face( ( *m_primitiveIndices )[i] );
                    out << depth << ": face " << ( *m_primitiveIndices )[i] << " / " << face << " / "
                        << mesh.get_vertex( face.x ) << " " << mesh.get_vertex( face.y ) << " "
                        << mesh.get_vertex( face.z ) << "\n";
                }
            } else {
                out << depth << ": empty leaf node\n";
            }
        } else {
            char dimensions[] = { 'X', 'Y', 'Z' };
            out << depth << ": interior node " << dimensions[m_axisAndLeafFlag & 0x00000003] << " " << m_split << "\n";
            out << depth << ": left child\n";
            left_child()->dump_tree( out, depth + 1, mesh );
            out << depth << ": right child\n";
            right_child()->dump_tree( out, depth + 1, mesh );
        }
    }

    // This returns the first intersection between the ray and the mesh, or false if none is found.
    bool intersect_ray( const ray3f& ray, double tMin, double tMax, const trimesh3& mesh,
                        raytrace_intersection& outIntersection ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.

            if( m_primitiveIndices != 0 ) {
                bool found = false;

                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    int faceIndex = ( *m_primitiveIndices )[i];

                    const plane3f& plane = mesh.get_face_plane( faceIndex );

                    // Get the intersection with the plane
                    double distance = plane.get_distance_to_intersection( ray );

                    if( distance >= tMin && distance <= tMax ) {
                        // Check that the intersection is inside the bounding box (we're intersecting the ray with the
                        // intersection of the triangle and the bounding box)
                        vector3f pt = ray.at( distance );

                        const boundbox3f& faceBounds = mesh.get_face_bound_box( faceIndex );

                        if( faceBounds.contains( pt ) ) {
                            vector3f barycentricCoord = mesh.compute_barycentric_coordinates( faceIndex, pt );

                            // If all the coordinates are positive, it's inside the triangle
                            if( barycentricCoord.x >= -0.0000001f && barycentricCoord.y >= -0.0000001f &&
                                barycentricCoord.z >= -0.0000001f ) {
                                found = true;

                                outIntersection.ray = ray;
                                outIntersection.position = pt;
                                outIntersection.geometricNormal = plane.normal();
                                outIntersection.distance = distance;
                                outIntersection.faceIndex = faceIndex;
                                outIntersection.barycentricCoords = barycentricCoord;

                                // distance < maxDistance from the if before, so constrict maxDistance
                                // for the rest of the intersection tests.
                                tMax = (float)distance;
                            }
                        }
                    }
                }

                return found;
            } else {
                return false;
            }
        } else {
            // It's an interior node.
            int axis = m_axisAndLeafFlag & 0x00000003;

            // If the ray is parallel to the plane, the distance computed would be invalid, so we special case it
            float rayDirectionAlongAxis = ray.direction()[axis];
            if( rayDirectionAlongAxis == 0 ) {
                // When the ray is parallel, we only need to check one of the sides
                if( ray.origin()[axis] <= m_split ) {
                    return left_child()->intersect_ray( ray, tMin, tMax, mesh, outIntersection );
                } else {
                    return right_child()->intersect_ray( ray, tMin, tMax, mesh, outIntersection );
                }
            } else {
                double distance = ( (double)m_split - ray.origin()[axis] ) / rayDirectionAlongAxis;

                // Traverse the children in the order so the closest child is tried first
                bool result = false;
                if( ray.direction()[axis] > 0 ) {
                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the right side
                        return right_child()->intersect_ray( ray, tMin, tMax, mesh, outIntersection );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the left side
                        return left_child()->intersect_ray( ray, tMin, tMax, mesh, outIntersection );
                    } else {
                        result = left_child()->intersect_ray( ray, tMin, distance, mesh, outIntersection );
                        if( !result )
                            result = right_child()->intersect_ray( ray, distance, tMax, mesh, outIntersection );
                    }
                } else {
                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the left side
                        return left_child()->intersect_ray( ray, tMin, tMax, mesh, outIntersection );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the right side
                        return right_child()->intersect_ray( ray, tMin, tMax, mesh, outIntersection );
                    } else {
                        result = right_child()->intersect_ray( ray, tMin, distance, mesh, outIntersection );
                        if( !result )
                            result = left_child()->intersect_ray( ray, distance, tMax, mesh, outIntersection );
                    }
                }

                return result;
            }
        }
    }

    // This returns an array of all the intersections between the ray and the mesh.
    void intersect_ray_all( const ray3f& ray, double tMin, double tMax, const trimesh3& mesh,
                            std::vector<raytrace_intersection>& outIntersections ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices != 0 ) {
                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    int faceIndex = ( *m_primitiveIndices )[i];

                    plane3f plane = mesh.get_face_plane( faceIndex );

                    // Get the intersection with the plane
                    double distance = plane.get_distance_to_intersection( ray );
                    if( distance >= tMin && distance <= tMax ) {

                        // Check that the intersection is inside the bounding box (we're intersecting the ray with the
                        // intersection of the triangle and the bounding box)
                        vector3f pt = ray.at( distance );

                        boundbox3f faceBounds = mesh.get_face_bound_box( faceIndex );

                        if( faceBounds.contains( pt ) ) {

                            vector3f barycentricCoord = mesh.compute_barycentric_coordinates( faceIndex, pt );

                            // If all the coordinates are positive, it's inside the triangle
                            if( barycentricCoord.x >= -0.0000001f && barycentricCoord.y >= -0.0000001f &&
                                barycentricCoord.z >= -0.0000001f ) {
                                raytrace_intersection result;

                                result.ray = ray;
                                result.position = pt;
                                result.geometricNormal = plane.normal();
                                result.distance = distance;
                                result.faceIndex = faceIndex;
                                result.barycentricCoords = barycentricCoord;

                                outIntersections.push_back( result );
                            }
                        }
                    }
                }
            }
        } else {
            // It's an interor node
            int axis = m_axisAndLeafFlag & 0x00000003;

            // If the ray is parallel to the plane, the distance computed would be invalid, so we special case it
            float rayDirectionAlongAxis = ray.direction()[axis];
            if( rayDirectionAlongAxis == 0 ) {
                // When the ray is parallel, we only need to check one of the sides
                if( ray.origin()[axis] <= m_split ) {
                    left_child()->intersect_ray_all( ray, tMin, tMax, mesh, outIntersections );
                } else {
                    right_child()->intersect_ray_all( ray, tMin, tMax, mesh, outIntersections );
                }
            } else {
                double distance = ( (double)m_split - ray.origin()[axis] ) / rayDirectionAlongAxis;

                // Traverse the children in the order so the closest child is processed first,
                // thus putting all the intersections in ascending order
                if( ray.direction()[axis] > 0 ) {
                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the right side
                        right_child()->intersect_ray_all( ray, tMin, tMax, mesh, outIntersections );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the left side
                        left_child()->intersect_ray_all( ray, tMin, tMax, mesh, outIntersections );
                    } else {
                        left_child()->intersect_ray_all( ray, tMin, distance, mesh, outIntersections );
                        right_child()->intersect_ray_all( ray, distance, tMax, mesh, outIntersections );
                    }
                } else {
                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the right side
                        left_child()->intersect_ray_all( ray, tMin, tMax, mesh, outIntersections );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the left side
                        right_child()->intersect_ray_all( ray, tMin, tMax, mesh, outIntersections );
                    } else {
                        right_child()->intersect_ray_all( ray, tMin, distance, mesh, outIntersections );
                        left_child()->intersect_ray_all( ray, distance, tMax, mesh, outIntersections );
                    }
                }
            }
        }
    }

    // This returns whether or not the ray intersects with the mesh within the given distance.
    bool intersects_ray_segment( const ray3f& ray, double tMin, double tMax, const trimesh3& mesh ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices != 0 ) {
                for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                    int faceIndex = ( *m_primitiveIndices )[i];

                    const plane3f& plane = mesh.get_face_plane( faceIndex );

                    // Get the intersection with the plane
                    double distance = plane.get_distance_to_intersection( ray );
                    if( distance >= tMin && distance <= tMax ) {
                        // Check that the intersection is inside the bounding box (we're intersecting the ray with the
                        // intersection of the triangle and the bounding box)

                        vector3f pt = ray.at( distance );

                        const boundbox3f& faceBounds = mesh.get_face_bound_box( faceIndex );

                        if( faceBounds.contains( pt ) ) {
                            vector3f barycentricCoord = mesh.compute_barycentric_coordinates( faceIndex, pt );

                            // If all the coordinates are positive, it's inside the triangle
                            if( barycentricCoord.x >= -0.0000001f && barycentricCoord.y >= -0.0000001f &&
                                barycentricCoord.z >= -0.0000001f ) {
                                return true;
                            }
                        }
                    }
                }
            }

            return false;
        } else {
            // It's an interior node
            int axis = m_axisAndLeafFlag & 0x00000003;

            // If the ray is parallel to the plane, the distance computed would be invalid, so we special case it
            float rayDirectionAlongAxis = ray.direction()[axis];
            if( rayDirectionAlongAxis == 0 ) {
                // When the ray is parallel, we only need to check one of the sides
                if( ray.origin()[axis] <= m_split ) {
                    return left_child()->intersects_ray_segment( ray, tMin, tMax, mesh );
                } else {
                    return right_child()->intersects_ray_segment( ray, tMin, tMax, mesh );
                }
            } else {
                double distance = ( (double)m_split - ray.origin()[axis] ) / rayDirectionAlongAxis;

                // Traverse the children in the order so the closest child is tried first
                bool result = false;
                if( ray.direction()[axis] > 0 ) {
                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the right side
                        return right_child()->intersects_ray_segment( ray, tMin, tMax, mesh );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the left side
                        return left_child()->intersects_ray_segment( ray, tMin, tMax, mesh );
                    } else {
                        result = left_child()->intersects_ray_segment( ray, tMin, distance, mesh );
                        if( !result )
                            result = right_child()->intersects_ray_segment( ray, distance, tMax, mesh );
                    }
                } else {
                    if( distance < tMin ) {
                        // The plane is to the left of the ray segment, so only consider the left side
                        return left_child()->intersects_ray_segment( ray, tMin, tMax, mesh );
                    } else if( distance > tMax ) {
                        // The plane is to the right of the ray segment, so only consider the right side
                        return right_child()->intersects_ray_segment( ray, tMin, tMax, mesh );
                    } else {
                        result = right_child()->intersects_ray_segment( ray, tMin, distance, mesh );
                        if( !result )
                            result = left_child()->intersects_ray_segment( ray, distance, tMax, mesh );
                    }
                }
                return result;
            }
        }
    }

    // This returns the nearest point on the mesh.  If no point on the mesh is within maxDistance,
    // returns false.
    bool find_nearest_point( const vector3f& point, boundbox3f& nodeBounds, float maxDistance, const trimesh3& mesh,
                             nearest_point_search_result& outNearestPoint ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices == 0 ||
                vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >= maxDistance * maxDistance )
                return false;

            bool found = false;

            for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                int faceIndex = ( *m_primitiveIndices )[i];

                plane3f plane = mesh.get_face_plane( faceIndex );

                // Project the point onto the triangle's plane
                float distance = plane.get_signed_distance_to_plane( point );

                if( fabsf( distance ) < maxDistance ) {

                    // Check that the intersection is inside the bounding box (we're intersecting the ray with the
                    // intersection of the triangle and the bounding box)
                    vector3f projected = point - distance * plane.normal();

                    // Now get the barycentric coordinates of this projection
                    vector3f barycentricCoord = mesh.compute_barycentric_coordinates( faceIndex, projected );
                    if( barycentricCoord.is_inf() )
                        continue;

                    vector3 face = mesh.get_face( faceIndex );
                    vector3f A = mesh.get_vertex( face.x );
                    vector3f B = mesh.get_vertex( face.y );
                    vector3f C = mesh.get_vertex( face.z );

                    if( barycentricCoord.x >= 0 && barycentricCoord.y >= 0 && barycentricCoord.z >= 0 ) {
                        distance = std::abs( distance ); // Projected and barycentricCoord are both accurate
                    } else {
                        detail::nearest_point_on_triangle( projected, barycentricCoord, A, B, C );
                        distance = vector3f::distance( point, projected );
                    }

                    if( distance < maxDistance ) {
                        outNearestPoint.distance = distance;
                        outNearestPoint.position = projected;
                        outNearestPoint.geometricNormal = plane.normal();
                        outNearestPoint.faceIndex = faceIndex;
                        outNearestPoint.barycentricCoords = barycentricCoord;

                        maxDistance = distance;
                        found = true;
                    }
                }
            }
            return found;
        } else {
            // It's an interor node
            int axis = m_axisAndLeafFlag & 0x00000003;

            // If we're too far away from the bounding box, don't even bother looking.  Use clamp_nothrow to avoid
            // the is_empty check, for efficiency.
            if( vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >= maxDistance * maxDistance )
                return false;

            // Traverse the children so the closest child is tried first
            if( point[axis] < m_split ) {
                float savedBoundsValue = nodeBounds.maximum()[axis];
                nodeBounds.maximum()[axis] = m_split;
                // First find the nearest point in the left child
                bool found = left_child()->find_nearest_point( point, nodeBounds, maxDistance, mesh, outNearestPoint );
                nodeBounds.maximum()[axis] = savedBoundsValue;

                // If found, constrict the search radius for the search in the right child
                if( found )
                    maxDistance = (float)outNearestPoint.distance;

                savedBoundsValue = nodeBounds.minimum()[axis];
                nodeBounds.minimum()[axis] = m_split;
                // Then find the nearest point in the right child, using this new constrained radius
                if( right_child()->find_nearest_point( point, nodeBounds, maxDistance, mesh, outNearestPoint ) )
                    found = true;
                nodeBounds.minimum()[axis] = savedBoundsValue;

                return found;
            } else {
                float savedBoundsValue = nodeBounds.minimum()[axis];
                nodeBounds.minimum()[axis] = m_split;
                // First find the nearest point in the right child
                bool found = right_child()->find_nearest_point( point, nodeBounds, maxDistance, mesh, outNearestPoint );
                nodeBounds.minimum()[axis] = savedBoundsValue;

                // If found, constrict the search radius for the search in the left child
                if( found )
                    maxDistance = (float)outNearestPoint.distance;

                savedBoundsValue = nodeBounds.maximum()[axis];
                nodeBounds.maximum()[axis] = m_split;
                // Then find the nearest point in the left child, using this new constrained radius
                if( left_child()->find_nearest_point( point, nodeBounds, maxDistance, mesh, outNearestPoint ) )
                    found = true;
                nodeBounds.maximum()[axis] = savedBoundsValue;

                return found;
            }
        }
    }

    void collect_nearest_faces( const vector3f& point, boundbox3f& nodeBounds, int nFaces, const trimesh3& mesh,
                                std::vector<nearest_point_search_result>& outNearestPoints ) {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            // If this leaf is empty, return
            if( m_primitiveIndices == 0 )
                return;

            // If there cannot be faces inside this region, return.
            if( (int)( outNearestPoints.size() ) == nFaces &&
                math::square( outNearestPoints.back().distance ) <
                    vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) )
                return;

            for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                int faceIndex = ( *m_primitiveIndices )[i];

                plane3f plane = mesh.get_face_plane( faceIndex );

                // Project the point onto the triangle's plane
                float distance = plane.get_signed_distance_to_plane( point );

                int n = (int)outNearestPoints.size();
                if( n < nFaces || fabsf( distance ) < outNearestPoints[n - 1].distance ) {

                    // Check that the intersection is inside the bounding box (we're intersecting the ray with the
                    // intersection of the triangle and the bounding box)
                    vector3f projected = point - distance * plane.normal();

                    // Now get the barycentric coordinates of this projection
                    vector3f barycentricCoord = mesh.compute_barycentric_coordinates( faceIndex, projected );

                    vector3 face = mesh.get_face( faceIndex );
                    vector3f A = mesh.get_vertex( face.x );
                    vector3f B = mesh.get_vertex( face.y );
                    vector3f C = mesh.get_vertex( face.z );

                    if( barycentricCoord.x >= 0 && barycentricCoord.y >= 0 && barycentricCoord.z >= 0 ) {
                        distance = std::abs( distance ); // Projected and barycentricCoord are both correct already
                    } else {
                        detail::nearest_point_on_triangle( projected, barycentricCoord, A, B, C );
                        distance = vector3f::distance( point, projected );
                    }

                    // This is mildly complicated.
                    // 1. If there are no faces, just add this one.
                    // 2. If there are less than nFaces,
                    //   a. If this face is further than any other, just add this one.
                    //   b. Otherwise, find the correct location for this face.
                    //       i. If this face is already in the list (They are not unique in the kd-tree), discard it.
                    //       ii. Otherwise, push back the other faces and insert this face.
                    // 3. Otherwise,
                    //   a. if this face is further than the last face, discard this face
                    //   b. Otherwise, find the correct location for this face.
                    //       i. If this face is already in the list discard it.
                    //       ii. Otherwise, insert this face in the correct location, discarding the last face.
                    //
                    // NOTE: I am assuming that if the same face is found twice, that the distance computed will
                    //       be exactly equal. I think this is safe since the exact same code-path and data will be
                    //       used to calculate the distances.
                    // TODO: There is a subtle bug where if two distinct faces x and y are added, then x is re-added
                    //       it won't catch the duplication because the first triangle at equal distance is not x.
                    if( n > 0 && distance < outNearestPoints[n - 1].distance ) {
                        std::size_t i = n - 1;
                        while( i > 0 && distance < outNearestPoints[i - 1].distance )
                            --i;
                        if( i == 0 || faceIndex != outNearestPoints[i - 1].faceIndex ) {
                            if( n < nFaces )
                                outNearestPoints.resize( n + 1 );
                            for( std::size_t j = outNearestPoints.size() - 1; j > i; --j )
                                outNearestPoints[j] = outNearestPoints[j - 1];
                            outNearestPoints[i].distance = distance;
                            outNearestPoints[i].position = projected;
                            outNearestPoints[i].geometricNormal = plane.normal();
                            outNearestPoints[i].faceIndex = faceIndex;
                            outNearestPoints[i].barycentricCoords = barycentricCoord;
                        }
                    } else if( n < nFaces && ( n == 0 || faceIndex != outNearestPoints.back().faceIndex ) ) {
                        outNearestPoints.resize( n + 1 );
                        outNearestPoints.back().distance = distance;
                        outNearestPoints.back().position = projected;
                        outNearestPoints.back().geometricNormal = plane.normal();
                        outNearestPoints.back().faceIndex = faceIndex;
                        outNearestPoints.back().barycentricCoords = barycentricCoord;
                    }
                }
            } // for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
        } else {
            // It's an interor node
            int axis = m_axisAndLeafFlag & 0x00000003;

            // If we're too far away from the bounding box, don't even bother looking.  Use clamp_nothrow to avoid
            // the is_empty check, for efficiency.
            if( (int)( outNearestPoints.size() ) == nFaces &&
                vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >=
                    math::square( outNearestPoints.back().distance ) )
                return;

            // Traverse the children so the closest child is tried first
            if( point[axis] < m_split ) {
                float savedBoundsValue = nodeBounds.maximum()[axis];
                nodeBounds.maximum()[axis] = m_split;
                // First find the nearest point in the left child
                left_child()->collect_nearest_faces( point, nodeBounds, nFaces, mesh, outNearestPoints );
                nodeBounds.maximum()[axis] = savedBoundsValue;

                savedBoundsValue = nodeBounds.minimum()[axis];
                nodeBounds.minimum()[axis] = m_split;
                // Then find the nearest point in the right child, using this new constrained radius
                right_child()->collect_nearest_faces( point, nodeBounds, nFaces, mesh, outNearestPoints );
                nodeBounds.minimum()[axis] = savedBoundsValue;
            } else {
                float savedBoundsValue = nodeBounds.minimum()[axis];
                nodeBounds.minimum()[axis] = m_split;
                // First find the nearest point in the right child
                right_child()->collect_nearest_faces( point, nodeBounds, nFaces, mesh, outNearestPoints );
                nodeBounds.minimum()[axis] = savedBoundsValue;

                savedBoundsValue = nodeBounds.maximum()[axis];
                nodeBounds.maximum()[axis] = m_split;
                // Then find the nearest point in the left child, using this new constrained radius
                left_child()->collect_nearest_faces( point, nodeBounds, nFaces, mesh, outNearestPoints );
                nodeBounds.maximum()[axis] = savedBoundsValue;
            }
        }
    }

    void collect_nearest_points_within_range( const vector3f& point, boundbox3f& nodeBounds, float maxDistance,
                                              const trimesh3& mesh,
                                              std::vector<nearest_point_search_result>& outResults ) const {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices == 0 ||
                vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >= maxDistance * maxDistance )
                return;

            for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                int faceIndex = ( *m_primitiveIndices )[i];

                plane3f plane = mesh.get_face_plane( faceIndex );

                // Project the point onto the triangle's plane
                float planeDistance = plane.get_signed_distance_to_plane( point );
                if( planeDistance < maxDistance ) {
                    // Check that the intersection is inside the bounding box (we're intersecting the ray with the
                    // intersection of the triangle and the bounding box)
                    vector3f projected = point - planeDistance * plane.normal();

                    // Now get the barycentric coordinates of this intersection
                    vector3f barycentricCoord = mesh.compute_barycentric_coordinates( faceIndex, projected );

                    vector3 face = mesh.get_face( faceIndex );
                    vector3f A = mesh.get_vertex( face.x );
                    vector3f B = mesh.get_vertex( face.y );
                    vector3f C = mesh.get_vertex( face.z );

                    detail::nearest_point_on_triangle( projected, barycentricCoord, A, B, C );

                    if( vector3f::distance_squared( point, projected ) < maxDistance * maxDistance ) {
                        outResults.resize( outResults.size() + 1 );

                        outResults.back().barycentricCoords = barycentricCoord;
                        outResults.back().distance = vector3f::distance( point, projected );
                        outResults.back().faceIndex = faceIndex;
                        outResults.back().geometricNormal = plane.normal();
                        outResults.back().position = projected;
                    }
                }
            }
            return;
        } else {
            // It's an interor node
            int axis = m_axisAndLeafFlag & 0x00000003;

            // If we're too far away from the bounding box, don't even bother looking.  Use clamp_nothrow to avoid
            // the is_empty check, for efficiency.
            if( vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >= maxDistance * maxDistance )
                return;

            // Traverse the children
            float savedBoundsValue = nodeBounds.maximum()[axis];
            nodeBounds.maximum()[axis] = m_split;
            left_child()->collect_nearest_points_within_range( point, nodeBounds, maxDistance, mesh, outResults );
            nodeBounds.maximum()[axis] = savedBoundsValue;

            savedBoundsValue = nodeBounds.minimum()[axis];
            nodeBounds.minimum()[axis] = m_split;
            right_child()->collect_nearest_points_within_range( point, nodeBounds, maxDistance, mesh, outResults );
            nodeBounds.minimum()[axis] = savedBoundsValue;

            return;
        }
    }

    // This collects all triangles indices within the sphere with center "point" and radius "maxDistance"
    void collect_faces_within_range( const vector3f& point, boundbox3f& nodeBounds, float maxDistance,
                                     const trimesh3& mesh, std::vector<int>& outFaces ) const {
        if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
            if( m_primitiveIndices == 0 ||
                vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >= maxDistance * maxDistance )
                return;

            for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
                int faceIndex = ( *m_primitiveIndices )[i];

                plane3f plane = mesh.get_face_plane( faceIndex );

                // Project the point onto the triangle's plane
                float planeDistance = plane.get_signed_distance_to_plane( point );
                if( planeDistance < maxDistance ) {
                    // Check that the intersection is inside the bounding box (we're intersecting the ray with the
                    // intersection of the triangle and the bounding box)
                    vector3f projected = point - planeDistance * plane.normal();

                    // Now get the barycentric coordinates of this intersection
                    vector3f barycentricCoord = mesh.compute_barycentric_coordinates( faceIndex, projected );

                    vector3 face = mesh.get_face( faceIndex );
                    vector3f A = mesh.get_vertex( face.x );
                    vector3f B = mesh.get_vertex( face.y );
                    vector3f C = mesh.get_vertex( face.z );

                    detail::nearest_point_on_triangle( projected, barycentricCoord, A, B, C );

                    if( vector3f::distance_squared( point, projected ) < maxDistance * maxDistance )
                        outFaces.push_back( faceIndex );
                }
            }
            return;
        } else {
            // It's an interor node
            int axis = m_axisAndLeafFlag & 0x00000003;

            // If we're too far away from the bounding box, don't even bother looking.  Use clamp_nothrow to avoid
            // the is_empty check, for efficiency.
            if( vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >= maxDistance * maxDistance )
                return;

            // Traverse the children
            float savedBoundsValue = nodeBounds.maximum()[axis];
            nodeBounds.maximum()[axis] = m_split;
            left_child()->collect_faces_within_range( point, nodeBounds, maxDistance, mesh, outFaces );
            nodeBounds.maximum()[axis] = savedBoundsValue;

            savedBoundsValue = nodeBounds.minimum()[axis];
            nodeBounds.minimum()[axis] = m_split;
            right_child()->collect_faces_within_range( point, nodeBounds, maxDistance, mesh, outFaces );
            nodeBounds.minimum()[axis] = savedBoundsValue;

            return;
        }
    }
};

} // namespace geometry
} // namespace frantic
