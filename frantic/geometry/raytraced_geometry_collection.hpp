// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/geometry/primitive_kdtree.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>
#include <frantic/graphics/motion_blurred_transform.hpp>

namespace frantic {
namespace geometry {
using frantic::graphics::motion_blurred_transform;
using frantic::graphics::transform4f;

// This provides an interface for a single raytraced object which is undergoing a rigid transformation.
class animated_rigid_raytraced_mesh {
    trimesh3 m_staticGeometry;
    motion_blurred_transform<float> m_transform;
    boost::shared_ptr<trimesh3_kdtree> m_staticGeometryKdtree;
    float m_cachedMotionSegmentTime;
    transform4f m_cachedTransform, m_cachedInverse;

  public:
    animated_rigid_raytraced_mesh() {
        // Start it at an unreasonable value
        m_cachedMotionSegmentTime = ( std::numeric_limits<float>::max )();
    }

    ///// Setters /////

    void clear() {
        m_staticGeometryKdtree.reset();
        m_staticGeometry.clear();
        m_transform.set_to_identity();
        m_cachedMotionSegmentTime = 0;
        m_cachedTransform.set_to_identity();
        m_cachedInverse.set_to_identity();
    }

    void set( const motion_blurred_transform<float>& xform, const trimesh3& mesh ) {
        m_staticGeometryKdtree.reset();
        m_staticGeometry = mesh;
        m_transform = xform;
    }

    void set_with_swap( const motion_blurred_transform<float>& xform, trimesh3& mesh ) {
        m_staticGeometryKdtree.reset();
        m_staticGeometry.swap( mesh );
        m_transform = xform;
    }

    ///// Getters /////

    boundbox3f get_bounds() const { return m_staticGeometry.compute_bound_box( m_cachedTransform ); }

    const motion_blurred_transform<float>& get_motion_blurred_transform() const { return m_transform; }

    void print_statistics( std::ostream& out, int objectNumber ) const {
        using namespace std;

        out << " Rigid geometry " << objectNumber << " has " << m_staticGeometry.face_count() << " faces and "
            << m_staticGeometry.vertex_count() << " vertices" << endl;
        if( m_staticGeometryKdtree ) {
            out << " Rigid geometry " << objectNumber << " kd-tree has " << m_staticGeometryKdtree->get_node_count()
                << " nodes, with a maximum depth of " << m_staticGeometryKdtree->get_maximum_depth() << endl;
        } else {
            out << " Rigid geometry " << objectNumber << " kd-tree hasn't been generated." << endl;
        }
    }

    ///// Operaters /////

    void add_object( const transform4f& xform, const trimesh3& mesh ) {
        if( mesh.face_count() > 0 ) {
            // invalidate the kd tree
            m_staticGeometryKdtree.reset();
            // add the new mesh
            m_staticGeometry.combine( xform, mesh );
        }
    }

    void prepare_kdtrees( float motionSegmentTime = 0.5f, bool forceRefresh = false ) {
        if( m_cachedMotionSegmentTime != motionSegmentTime || forceRefresh ) {
            m_cachedMotionSegmentTime = motionSegmentTime;
            m_cachedTransform = m_transform.get_transform( m_cachedMotionSegmentTime );
            m_cachedInverse = m_transform.get_inverse_transform( m_cachedMotionSegmentTime );
            if( m_staticGeometryKdtree.get() == 0 || forceRefresh )
                m_staticGeometryKdtree.reset( new trimesh3_kdtree( m_staticGeometry ) );
        }

        //		std::cout << "prepared kdtree at time " << motionSegmentTime << std::endl;
    }

    ///// Queries /////

    bool intersects_ray_segment( const vector3f& start, const vector3f& end, float motionSegmentTime ) const {
        if( m_staticGeometryKdtree.get() == 0 )
            throw std::runtime_error(
                "animated_rigid_raytraced_mesh.intersects_ray_segment: Tried to check whether a mesh "
                "intersects a ray segment without first initializing the kdtree." );
        vector3f transformedStart, transformedEnd;
        // If the requested time is the cached time, use the cached matrices
        if( fabs( motionSegmentTime - m_cachedMotionSegmentTime ) < 0.0001f ) {
            transformedStart = m_cachedInverse * start;
            transformedEnd = m_cachedInverse * end;
        } else {
            transformedStart = m_transform.inverse_transform_point( start, motionSegmentTime );
            transformedEnd = m_transform.inverse_transform_point( end, motionSegmentTime );
        }
        return m_staticGeometryKdtree->intersects_ray_segment( transformedStart, transformedEnd );
    }

    bool intersects_ray_segment( const ray3f& ray, double tMin, double tMax ) const {
        if( m_staticGeometryKdtree.get() == 0 )
            throw std::runtime_error(
                "animated_rigid_raytraced_mesh.intersects_ray_segment: Tried to check whether a mesh "
                "intersects a ray segment without first initializing the kdtree." );

        return m_staticGeometryKdtree->intersects_ray_segment( m_cachedInverse * ray, tMin, tMax );
    }

    bool intersect_ray( const ray3f& ray, double tMin, double tMax, float motionSegmentTime,
                        frantic::geometry::raytrace_intersection& outIntersection ) const {
        if( m_staticGeometryKdtree.get() == 0 )
            throw std::runtime_error(
                "animated_rigid_raytraced_mesh.intersect_ray: Tried to do a ray intersection without "
                "first initializing the kdtree." );
        ray3f transformedRay;
        // If the requested time is the cached time, use the cached matrices
        if( fabs( motionSegmentTime - m_cachedMotionSegmentTime ) < 0.0001f ) {
            transformedRay = m_cachedInverse * ray;
        } else {
            transformedRay = m_transform.get_inverse_transform( motionSegmentTime ) * ray;
        }
        if( m_staticGeometryKdtree->intersect_ray( transformedRay, tMin, tMax, outIntersection ) ) {
            outIntersection.position = m_cachedTransform * outIntersection.position;
            // Translate the normal using the correct method (transpose of the inverse)
            outIntersection.geometricNormal =
                m_cachedInverse.transpose_transform_no_translation( outIntersection.geometricNormal );
            outIntersection.geometricNormal.normalize();
            return true;
        } else {
            return false;
        }
    }

    bool intersect_ray( const ray3f& ray, double tMin, double tMax,
                        frantic::geometry::raytrace_intersection& outIntersection ) const {
        if( m_staticGeometryKdtree.get() == 0 )
            throw std::runtime_error(
                "animated_rigid_raytraced_mesh.intersect_ray: Tried to do a ray intersection without "
                "first initializing the kdtree." );
        if( m_staticGeometryKdtree->intersect_ray( m_cachedInverse * ray, tMin, tMax, outIntersection ) ) {
            outIntersection.position = m_cachedTransform * outIntersection.position;
            // Translate the normal using the correct method (transpose of the inverse)
            outIntersection.geometricNormal =
                m_cachedInverse.transpose_transform_no_translation( outIntersection.geometricNormal );
            outIntersection.geometricNormal.normalize();
            return true;
        } else {
            return false;
        }
    }

    void intersect_ray_all( const ray3f& ray, float segmentLength, float motionSegmentTime,
                            std::vector<frantic::geometry::raytrace_intersection>& outIntersections ) const {
        if( m_staticGeometryKdtree.get() == 0 )
            throw std::runtime_error(
                "animated_rigid_raytraced_mesh.intersect_ray: Tried to do a ray intersection without "
                "first initializing the kdtree." );
        ray3f transformedRay;
        // If the requested time is the cached time, use the cached matrices
        if( fabs( motionSegmentTime - m_cachedMotionSegmentTime ) < 0.0001f ) {
            transformedRay = m_cachedInverse * ray;
        } else {
            transformedRay = m_transform.get_inverse_transform( motionSegmentTime ) * ray;
        }
        std::vector<frantic::geometry::raytrace_intersection> tempIntersections;
        m_staticGeometryKdtree->intersect_ray_all( transformedRay, segmentLength, tempIntersections );
        for( unsigned i = 0; i < tempIntersections.size(); ++i ) {
            tempIntersections[i].position = m_transform.transform_point( outIntersections[i].position );
            tempIntersections[i].geometricNormal = m_transform.transform_normal( outIntersections[i].geometricNormal );
            tempIntersections[i].geometricNormal.normalize();
        }
        std::copy( tempIntersections.begin(), tempIntersections.end(), std::back_inserter( outIntersections ) );
        std::sort( outIntersections.begin(), outIntersections.end() );
    }
};

// This class represents a number of raytraced objects with possible rigid motion.
// For speed, all the static objects are combined into a single triangle mesh, while all the
// moving objects are kept separate.
class raytraced_geometry_collection {
    trimesh3 m_staticGeometry;
    boost::shared_ptr<trimesh3_kdtree> m_staticGeometryKdtree;
    std::vector<boost::shared_ptr<animated_rigid_raytraced_mesh>> m_rigidGeometryObjects;

    boost::shared_ptr<primitive_kdtree<animated_rigid_raytraced_mesh>> m_masterKdtree;

    float m_cachedMotionSegmentTime;

  public:
    raytraced_geometry_collection() {
        // Start it at an unreasonable value
        m_cachedMotionSegmentTime = ( std::numeric_limits<float>::max )();
    }

    void clear() {
        m_staticGeometryKdtree.reset();
        m_staticGeometry.clear();
        m_rigidGeometryObjects.clear();
        m_masterKdtree.reset();
    }

    void swap( raytraced_geometry_collection& rhs ) {
        m_staticGeometry.swap( rhs.m_staticGeometry );
        m_staticGeometryKdtree.swap( rhs.m_staticGeometryKdtree );
        m_rigidGeometryObjects.swap( rhs.m_rigidGeometryObjects );
        m_masterKdtree.swap( rhs.m_masterKdtree );
        std::swap( m_cachedMotionSegmentTime, rhs.m_cachedMotionSegmentTime );
    }

    void add_static_object( const transform4f& xform, const trimesh3& mesh ) {
        if( mesh.face_count() > 0 ) {
            // invalidate the kd tree
            m_staticGeometryKdtree.reset();
            // add the new mesh
            m_staticGeometry.combine( xform, mesh );
        }
    }

    void add_rigid_object( const motion_blurred_transform<float>& xform, const trimesh3& mesh ) {
        // Try to add it to an existing rigid object, but if not add a new one
        if( mesh.face_count() > 0 ) {
            if( !add_rigid_object_to_existing( xform, mesh ) ) {
                m_rigidGeometryObjects.push_back(
                    boost::shared_ptr<animated_rigid_raytraced_mesh>( new animated_rigid_raytraced_mesh() ) );
                m_rigidGeometryObjects.back()->set( xform, mesh );
            }
        }
    }

    void add_rigid_object_with_swap( const motion_blurred_transform<float>& xform, trimesh3& mesh ) {
        // Try to add it to an existing rigid object, but if not add a new one
        if( mesh.face_count() > 0 ) {
            if( !add_rigid_object_to_existing( xform, mesh ) ) {
                m_rigidGeometryObjects.push_back(
                    boost::shared_ptr<animated_rigid_raytraced_mesh>( new animated_rigid_raytraced_mesh() ) );
                m_rigidGeometryObjects.back()->set_with_swap( xform, mesh );
            }
        }
    }

  private:
    bool add_rigid_object_to_existing( const motion_blurred_transform<float>& xform, const trimesh3& mesh ) {
        if( xform.is_static() ) {
            // If the transform is static, add it to the collection of static objects
            add_static_object( xform.get_transform(), mesh );
            return true;
        } else {
            // Otherwise look through all the rigid objects for a compatible transform animation, and if one is found
            // add it to that
            /*
            for( unsigned i = 0; i < m_rigidGeometryObjects.size(); ++i ) {
              if( xform.is_rigid_transformation_of( m_rigidGeometryObjects[i]->get_motion_blurred_transform() ) ) {
                transform4f relativeTransform =
            m_rigidGeometryObjects[i]->get_motion_blurred_transform().get_inverse_transform() * xform.get_transform();
                m_rigidGeometryObjects[i]->add_object( relativeTransform, mesh );
                return true;
              }
            }
            //*/
            return false;
        }
    }

  public:
    void prepare_kdtrees( float motionSegmentTime = 0.5f, bool forceRefresh = false ) {
        if( m_cachedMotionSegmentTime != motionSegmentTime || forceRefresh ) {
            // The static geometry
            if( ( !m_staticGeometryKdtree || forceRefresh ) && m_staticGeometry.face_count() > 0 )
                m_staticGeometryKdtree.reset( new trimesh3_kdtree( m_staticGeometry ) );

            // m_staticGeometryKdtree->dump_tree(std::cout);

            // Destroy the old master kd-tree
            m_masterKdtree.reset();

            for( std::vector<boost::shared_ptr<animated_rigid_raytraced_mesh>>::iterator i =
                     m_rigidGeometryObjects.begin();
                 i != m_rigidGeometryObjects.end(); ++i ) {
                ( *i )->prepare_kdtrees( motionSegmentTime, forceRefresh );
            }

            // Build a new master kd-tree
            m_masterKdtree.reset( new primitive_kdtree<animated_rigid_raytraced_mesh>( m_rigidGeometryObjects ) );

            m_cachedMotionSegmentTime = motionSegmentTime;
        }
    }

    boundbox3f get_bounds() const {
        boundbox3f box = m_staticGeometry.compute_bound_box();
        for( unsigned i = 0; i < m_rigidGeometryObjects.size(); ++i )
            box += m_rigidGeometryObjects[i]->get_bounds();
        return box;
    }

    bool empty() const { return m_staticGeometry.face_count() == 0 && m_rigidGeometryObjects.size() == 0; }

    void print_statistics( std::ostream& out ) const {
        using namespace std;

        out << " ===== Raytraced Geometry Collection Statistics =====" << endl;
        out << " Static geometry has " << m_staticGeometry.face_count() << " faces and "
            << m_staticGeometry.vertex_count() << " vertices" << endl;
        if( m_staticGeometryKdtree ) {
            out << " Static geometry kd-tree has " << m_staticGeometryKdtree->get_node_count()
                << " nodes, with a maximum depth of " << m_staticGeometryKdtree->get_maximum_depth() << endl;
            // m_staticGeometryKdtree->dump_tree(std::cout);
        } else {
            out << " Static geometry kd-tree hasn't been generated." << endl;
        }
        for( unsigned i = 0; i < m_rigidGeometryObjects.size(); ++i ) {
            m_rigidGeometryObjects[i]->print_statistics( out, i );
        }
    }

    bool intersects_ray_segment( const ray3f& ray, double tMin, double tMax ) const {
        if( m_staticGeometry.face_count() > 0 ) {
            if( m_staticGeometryKdtree.get() == 0 )
                throw std::runtime_error(
                    "raytraced_geometry_collection.intersects_ray_segment: Tried to check whether a mesh "
                    "intersects a ray segment without first initializing the kdtrees." );

            if( m_staticGeometryKdtree->intersects_ray_segment( ray, tMin, tMax ) )
                return true;
        }

        if( m_masterKdtree->intersects_ray_segment( ray, tMin, tMax ) )
            return true;

        return false;
    }

    bool intersects_ray_segment( const vector3f& start, const vector3f& end ) const {
        // Don't let the NaNs and Infs slip into here, they cause problems!
        if( !( start.is_finite() && end.is_finite() ) )
            return false;

        ray3f ray( start, end - start );
        return intersects_ray_segment( ray, 0, 1 );
    }

    bool intersect_ray( const ray3f& ray, double tMin, double tMax,
                        frantic::geometry::raytrace_intersection& outIntersection ) const {
        bool intersected = false;

        if( m_staticGeometry.face_count() > 0 ) {
            if( m_staticGeometryKdtree.get() == 0 )
                throw std::runtime_error(
                    "raytraced_geometry_collection.intersect_ray: Tried to intersect a ray with a mesh "
                    "without first initializing the kdtrees." );

            if( m_staticGeometryKdtree->intersect_ray( ray, tMin, tMax, outIntersection ) ) {
                tMax = outIntersection.distance;
                intersected = true;
            }
        }

        if( m_masterKdtree.get() != 0 ) {
            //			std::cout << "intersecting obj with ray " << ray << ", " << tMin << ", " << tMax <<
            // std::endl;
            if( m_masterKdtree->intersect_ray( ray, tMin, tMax, outIntersection ) ) {
                //				std::cout << "intersected!" << std::endl;
                tMax = outIntersection.distance;
                intersected = true;
            }
        }

        return intersected;
    }

    void intersect_ray_all( const ray3f& ray, float segmentLength, float motionSegmentTime,
                            std::vector<frantic::geometry::raytrace_intersection>& outIntersections ) const {
        if( m_staticGeometry.face_count() > 0 ) {
            if( m_staticGeometryKdtree.get() == 0 )
                throw std::runtime_error(
                    "raytraced_geometry_collection.intersect_ray: Tried to intersect a ray with a mesh "
                    "without first initializing the kdtrees." );

            m_staticGeometryKdtree->intersect_ray_all( ray, segmentLength, outIntersections );
        }

        for( std::vector<boost::shared_ptr<animated_rigid_raytraced_mesh>>::const_iterator i =
                 m_rigidGeometryObjects.begin();
             i != m_rigidGeometryObjects.end(); ++i ) {
            ( *i )->intersect_ray_all( ray, segmentLength, motionSegmentTime, outIntersections );
        }
    }
};

} // namespace geometry
} // namespace frantic
