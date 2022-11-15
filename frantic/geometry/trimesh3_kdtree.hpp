// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_kdtree_node.hpp>

namespace frantic {
namespace geometry {

class trimesh3_kdtree {
    // Pointer to the mesh object.  The mesh is owned by someone else, not by this kd-tree.  We may want to change this
    // to use a shared_ptr<> at some point.
    trimesh3* m_mesh;
    // Pointer to the kd tree.  This member is owned by this class, so is both allocated and destroyed by it.
    trimesh3_kdtree_node* m_rootNode;

    boundbox3f m_rootBounds;

  public:
    friend class trimesh3_kdtree_ext;

    trimesh3_kdtree( trimesh3& mesh );

    ~trimesh3_kdtree() {
        if( m_rootNode != 0 ) {
            delete m_rootNode;
            m_rootNode = 0;
        }
    }

    const trimesh3& get_mesh() const { return *m_mesh; }

    // Returns the total number of nodes in the tree.
    int get_node_count() const { return m_rootNode->get_node_count(); }

    // Returns the depth of the deepest leaf node.
    int get_maximum_depth() const { return m_rootNode->get_maximum_depth(); }

    // Returns the largest count of primitives in a leaf node.
    int get_largest_leaf_size() const { return m_rootNode->get_largest_leaf_size(); }

    void dump_tree( std::ostream& out ) const {
        out << "----- Dumping KD-Tree\n";
        out << " boundbox: " << m_rootBounds << "\n";
        m_rootNode->dump_tree( out, 0, *m_mesh );
        out << "----- Finished dumping KD-Tree" << std::endl;
    }

    ///// Queries/////

    /**
     * This returns the first intersection between the ray and the mesh, or false if none is found.
     *
     * @returns ray intersection if found
     * @returns true if intersection found
     */
    bool intersect_ray( const ray3f& ray, double tMin, double tMax, raytrace_intersection& outIntersection ) const {
        /*
        std::cout << "Casting ray " << ray << "\n";
        //*/

        if( !ray.clamp_to_box( m_rootBounds, tMin, tMax ) )
            return false;

        return m_rootNode->intersect_ray( ray, tMin, tMax, *m_mesh, outIntersection );
    }

    // This adds the intersections found to the array, and maintains its sorted order
    void intersect_ray_all( const ray3f& ray, float segmentLength,
                            std::vector<raytrace_intersection>& outIntersections ) const {
        m_rootNode->intersect_ray_all( ray, 0, segmentLength, *m_mesh, outIntersections );
        std::sort( outIntersections.begin(), outIntersections.end() );
    }

    /**
     *  This returns whether or not the ray intersects with the mesh within the given distance.
     */
    bool intersects_ray_segment( const vector3f& start, const vector3f& end ) const {
        ray3f ray( start, end - start );
        return m_rootNode->intersects_ray_segment( ray, 0, 1, *m_mesh );
    }

    /**
     * This returns whether or not the ray intersects with the mesh within the given distance.
     */
    bool intersects_ray_segment( const ray3f& ray, double tMin, double tMax ) const {
        return m_rootNode->intersects_ray_segment( ray, tMin, tMax, *m_mesh );
    }

    // This returns the nearest point on the mesh.  If no point on the mesh is within maxDistance,
    // returns false.
    bool find_nearest_point( const vector3f& point, float maxDistance,
                             nearest_point_search_result& outNearestPoint ) const {
        boundbox3f bounds = m_rootBounds;
        return m_rootNode->find_nearest_point( point, bounds, maxDistance, *m_mesh, outNearestPoint );
    }

    void collect_nearest_faces( const vector3f& point, int nFaces,
                                std::vector<nearest_point_search_result>& outResults ) const {
        outResults.clear();

        if( nFaces <= 0 )
            return;

        boundbox3f bounds = m_rootBounds;
        return m_rootNode->collect_nearest_faces( point, bounds, nFaces, *m_mesh, outResults );
    }

    // This collects all triangles indices within the sphere with center "point" and radius "maxDistance"
    void collect_faces_within_range( const vector3f& point, float maxDistance, std::vector<int>& outFaces ) const {
        boundbox3f bounds = m_rootBounds;
        m_rootNode->collect_faces_within_range( point, bounds, maxDistance, *m_mesh, outFaces );

        if( outFaces.size() == 0 )
            return;

        // Sort the indices, then remove duplicates. This is O(nLogn) instead of O(n^2) when doing the naive
        // removal. Perhaps the naive approach is still better ... this could be tested.
        std::sort( outFaces.begin(), outFaces.end() );
        for( std::size_t i = outFaces.size() - 1; i > 0; --i ) {
            if( outFaces[i - 1] == outFaces[i] ) {
                outFaces[i] = outFaces.back();
                outFaces.pop_back();
            }
        }
    }

    void collect_nearest_points_within_range( const vector3f point, float maxDistance,
                                              std::vector<nearest_point_search_result>& outResults ) const {
        outResults.clear();
        boundbox3f bounds = m_rootBounds;
        m_rootNode->collect_nearest_points_within_range( point, bounds, maxDistance, *m_mesh, outResults );
    }

    /**
     * Returns the distance to the nearest point on the surface, and optionally returns the normal
     * at that point.
     *
     * @param  p                The point to test.
     * @param  pOutNearestPoint If non-null, this is populated with the point nearest the surface
     * @param  pOutNormal       If non-null, this is populated with the geometric normal.
     */
    float get_distance_to_surface( const vector3f& p, vector3f* pOutNearestPoint, vector3f* pOutNormal ) const {
        nearest_point_search_result npsr;
        if( !find_nearest_point( p, get_bounds().get_diagonal().get_magnitude(), npsr ) )
            return ( std::numeric_limits<float>::max )();

        if( pOutNearestPoint )
            *pOutNearestPoint = npsr.position;
        if( pOutNormal )
            *pOutNormal = npsr.geometricNormal;
        return npsr.distance;
    }

    /**
     * Returns true if the point is in the volume of the mesh, false otherwise.
     * The implementation first searches for a face on which the point in question does not lie.
     * Then two sample rays are generated and tested against the mesh. If both intersections face towards
     * or away from the point (using the face normal) then the result is returned. A third tie breaking
     * ray will be generated if there is disagreement btw. the first two rays. Finally if the third ray
     * is inconclusive (ie. it misses) then a new triangle is chosen and the algorithm repeats.
     *
     * @param  p  The point to test.
     */
    bool is_point_in_volume( const vector3f& p ) const {
        const trimesh3& mesh = get_mesh();
        const std::size_t numFaces = mesh.face_count();

        int result = 0;
        for( std::size_t i = 0; i < numFaces; ++i ) {
            const vector3 f = mesh.get_face( i );
            const vector3f v0 = mesh.get_vertex( f.x );
            const vector3f v1 = mesh.get_vertex( f.y );
            const vector3f v2 = mesh.get_vertex( f.z );

            // Verify this point doesn't lie in the triangle's plane.
            // TODO: This should be a relative error
            vector3f norm = vector3f::normalize( vector3f::cross( v1 - v0, v2 - v0 ) );
            if( std::abs( vector3f::dot( norm, ( p - v0 ) ) ) < 0.001f )
                continue;

            raytrace_intersection ri;

            vector3f faceCenter = ( v0 + v1 + v2 ) / 3.f;
            if( intersect_ray( ray3f( p, faceCenter - p ), 0.001, 1.001, ri ) ) {
                // We hit the triangle ... so now add a vote in either the positive or negative direction
                if( vector3f::dot( ri.geometricNormal, ri.ray.direction() ) > 0 ) {
                    if( ++result > 1 )
                        return true;
                } else {
                    if( --result < -1 )
                        return false;
                }
            }
        }

        // If all the faces were skipped, then the point is co-planar with the
        // entire mesh (or the mesh is empty).
        return false;
    }

    boundbox3f get_bounds() const { return m_rootBounds; }

  private:
    // Recursively generates the kd-tree from the mesh.  For efficiency, this destructively swaps out the passed
    // in faceIndices array if it is needed.
    void construct_kdtree_n_log2_n( std::vector<int>& faceIndices, const boundbox3f& bounds, int depth,
                                    trimesh3_kdtree_node& outNode ) const;
};

} // namespace geometry
} // namespace frantic
