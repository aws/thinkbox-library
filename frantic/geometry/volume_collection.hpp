// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/trimesh3_kdtree.hpp>
//#include <frantic/geometry/trimesh3_file_io.hpp>

#include <frantic/volumetrics/levelset/rle_level_set.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/thread/tss.hpp>

namespace frantic {
namespace geometry {

class volume_collection {
  public:
    virtual float get_segment_distance_in_volume( const vector3f& p1, const vector3f& p2 ) const = 0;
    virtual float get_distance_to_surface( const vector3f& p, vector3f* pOutNormal ) const = 0;
    virtual bool is_point_in_volume( const vector3f& p ) const = 0;
};

class trimesh3_kdtree_volume_collection : public volume_collection {
  private:
    trimesh3 m_mesh;
    std::unique_ptr<trimesh3_kdtree> m_pTree; // unique_ptr used because this is owned, not shared

    // I added this cached vector object because profiling showed that creating a std::vector
    // for each call to get_segment_distance_in_volume() was very expensive. Using a boost::thread_specific_ptr
    // allowed me to cache the allocation of the vector between calls while still being thread-safe.
    mutable boost::thread_specific_ptr<std::vector<raytrace_intersection>> m_pIntersections;

  public:
    trimesh3_kdtree_volume_collection( trimesh3& mesh ) {
        m_mesh.swap( mesh );
        m_pTree.reset( new trimesh3_kdtree( m_mesh ) );
    }

    float get_segment_distance_in_volume( const vector3f& p1, const vector3f& p2 ) const {
        vector3f rayDir( p2 - p1 );
        float distance = rayDir.get_magnitude();

        // Converted to a thread_specific_ptr because allocating this vector was expensive.
        // std::vector<raytrace_intersection> intersections;
        if( !m_pIntersections.get() ) {
            m_pIntersections.reset( new std::vector<raytrace_intersection> );
            m_pIntersections->reserve( 10 );
        }

        std::vector<raytrace_intersection>& intersections = *m_pIntersections;
        intersections.clear();

        m_pTree->intersect_ray_all( ray3f( p1, rayDir ), 1, intersections );

        // Determine if the ray is entirely inside or outside the volume
        if( intersections.size() == 0 )
            return is_point_in_volume( p1 ) ? distance : 0;

        double result = 0;
        std::vector<raytrace_intersection>::iterator it = intersections.begin();

        // If the ray starts inside
        if( vector3f::dot( it->geometricNormal, rayDir ) > 0 )
            result = ( it++ )->distance;

        while( it != intersections.end() ) {
            result -= ( it++ )->distance;
            if( it != intersections.end() )
                result += ( it++ )->distance;
        }

        // If the ray ends inside
        if( vector3f::dot( intersections.back().geometricNormal, rayDir ) < 0 )
            result += 1.0;

        result = math::clamp( result, 0.0, 1.0 ); // TODO: remove this hack
        /*if(result > 1.0 || result < 0){
          bool p1Inside = vector3f::dot(intersections.front().geometricNormal, rayDir) > 0;
          bool p2Inside = vector3f::dot(intersections.back().geometricNormal, rayDir) < 0;

          mprintf("p1 @ %s was %sinside.\n", p1.str().c_str(), p1Inside ? "" : "not ");
          mprintf("p2 @ %s was %sinside.\n", p2.str().c_str(), p2Inside ? "" : "not ");
          mprintf("Total distance was: %f\n", distance);
          mprintf("There were %u intersections\n", intersections.size());
          for(std::size_t i = 0; i < intersections.size(); ++i){
            float d = vector3f::dot(intersections[i].geometricNormal, rayDir) / distance;
            mprintf("\tIntersection[%u] @ %f w/ normal: %s and inner product %f\n", i, (float)intersections[i].distance,
        intersections[i].geometricNormal.str().c_str(), d);
          }
        }*/

        return static_cast<float>( result * distance );
    }

    float get_distance_to_surface( const vector3f& p, vector3f* pOutNormal ) const {
        return m_pTree->get_distance_to_surface( p, 0, pOutNormal );
    }

    bool is_point_in_volume( const vector3f& p ) const { return m_pTree->is_point_in_volume( p ); }
};

} // namespace geometry
} // namespace frantic
