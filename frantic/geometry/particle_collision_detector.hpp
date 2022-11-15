// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <frantic/geometry/mixed_kdtree.hpp>
#include <frantic/geometry/primitive_kdtree.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>
#include <frantic/particles/particle_kdtree.hpp>

#include <frantic/particles/particle_step.hpp>

namespace frantic {
namespace geometry {
using frantic::graphics::motion_blurred_transform;

// Avoid destructing the trimesh3 when the kdtree is done with it
struct null_deleter {
    void operator()( void const* ) const {}
};

class animated_particle_collision_mesh {
    trimesh3 m_staticGeometry;
    boost::shared_ptr<mixed_kdtree> m_staticGeometryKdtree;
    transform4f m_firstTransform, m_firstTransformInverse, m_secondTransform, m_secondTransformInverse;

  public:
    animated_particle_collision_mesh() {}

    ///// Setters /////

    void clear() {
        m_staticGeometryKdtree.reset();
        m_staticGeometry.clear();
        m_firstTransform.set_to_identity();
        m_firstTransformInverse.set_to_identity();
        m_secondTransform.set_to_identity();
        m_secondTransformInverse.set_to_identity();
    }

    void set( const motion_blurred_transform<float>& xform, const trimesh3& mesh ) {
        m_staticGeometryKdtree.reset();
        m_staticGeometry = mesh;
        m_firstTransform = xform.get_transform( 0 );
        m_firstTransformInverse = xform.get_inverse_transform( 0 );
        m_secondTransform = xform.get_transform( 1 );
        m_secondTransformInverse = xform.get_inverse_transform( 1 );
    }

    void set_with_swap( const motion_blurred_transform<float>& xform, trimesh3& mesh ) {
        m_staticGeometryKdtree.reset();
        m_staticGeometry.swap( mesh );
        m_firstTransform = xform.get_transform( 0 );
        m_firstTransformInverse = xform.get_inverse_transform( 0 );
        m_secondTransform = xform.get_transform( 1 );
        m_secondTransformInverse = xform.get_inverse_transform( 1 );
    }

    void set_transform( const motion_blurred_transform<float>& xform ) {
        m_firstTransform = xform.get_transform( 0 );
        m_firstTransformInverse = xform.get_inverse_transform( 0 );
        m_secondTransform = xform.get_transform( 1 );
        m_secondTransformInverse = xform.get_inverse_transform( 1 );
    }

    ///// Getters /////

    // Use the bounding box of the convex hull of the motion.
    boundbox3f get_bounds() const {
        boundbox3f result = m_staticGeometry.compute_bound_box( m_firstTransform );
        result += m_staticGeometry.compute_bound_box( m_secondTransform );
        return result;
    }

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

    ///// Operators /////

    void prepare_kdtrees( bool forceRefresh = false ) {
        if( m_staticGeometryKdtree.get() == 0 || forceRefresh ) {
            boost::shared_ptr<trimesh3> staticGeometry( &m_staticGeometry, null_deleter() );
            m_staticGeometryKdtree.reset( new mixed_kdtree( staticGeometry ) );
        }

        // std::cout << "prepared kdtree at time " << motionSegmentTime << std::endl;
    }

    ///// Queries /////

    // Detects a collision between the ray (defined as a particle moving from time 0 to time 1), and this moving mesh.
    bool intersect_ray( const ray3f& ray, double tMin, double tMax,
                        frantic::geometry::raytrace_intersection& outIntersection ) const {
        if( m_staticGeometryKdtree.get() == 0 )
            throw std::runtime_error(
                "animated_particle_collision_mesh.intersect_ray: Tried to do a particle-mesh collision "
                "check without first initializing the kdtree." );

        // Transform the particle path via the inverse of the object's path
        // vector3f transformedOrigin = m_firstTransformInverse * ray.origin();
        // ray3f transformedRay( transformedOrigin, m_secondTransformInverse * (ray.origin() + ray.direction()) -
        // transformedOrigin );
        vector3f transformedStart = ( ( 1.f - static_cast<float>( tMin ) ) * m_firstTransformInverse +
                                      static_cast<float>( tMin ) * m_secondTransformInverse ) *
                                    ray.at( tMin );
        vector3f transformedEnd = m_secondTransformInverse * ray.at( 1.0 );
        ray3f transformedRay( ray3f::from_line_segment(
            transformedStart, transformedEnd ) ); // transformedOrigin, m_secondTransformInverse * (ray.origin() +
                                                  // ray.direction()) - transformedOrigin );

        // this is a hack.  the kdtree intersection code needs to be changed because it doesn't currently
        // support negative t values.  changing that behaviour at this time would cause problems, and since
        // we're in the middle of a kdtree rewrite, this hack will suffice for now.
        /*
        double scale = tMax - tMin, offset = tMin;
        if ( offset < 0 ) {
          transformedRay.set_origin( transformedRay.origin() + transformedRay.direction()*float(offset));
          transformedRay.set_direction( transformedRay.direction()*float(1+tMin) );
          tMin = 0.0;
          tMax = 1.0;
        }
        */

        // if( m_staticGeometryKdtree->intersect_ray( transformedRay, tMin, tMax, outIntersection ) ) {
        if( m_staticGeometryKdtree->intersect_ray( transformedRay, 0.0, 1.0, outIntersection ) ) {
            /*
            if ( offset < 0 ) {
              outIntersection.distance = outIntersection.distance*scale + offset;
            }
            */

            outIntersection.distance = tMin + outIntersection.distance * ( tMax - tMin );

            // Compute the transform matrix and its inverse at the intersection time, and use them to compute the
            // position and normal of the intersection
            transform4f transform = transform4f::linear_interpolate( m_firstTransform, m_secondTransform,
                                                                     static_cast<float>( outIntersection.distance ) );
            transform4f transformInverse = transform4f::linear_interpolate(
                m_firstTransformInverse, m_secondTransformInverse, static_cast<float>( outIntersection.distance ) );

            // Compute the motion of the intersected point
            outIntersection.motionDuringTimeStep = ( m_secondTransform - m_firstTransform ) * outIntersection.position;

            outIntersection.position = transform * outIntersection.position;
            outIntersection.ray = transform * outIntersection.ray;

            // Translate the normal using the correct method (transpose of the inverse)
            outIntersection.geometricNormal =
                transformInverse.transpose_transform_no_translation( outIntersection.geometricNormal );
            outIntersection.geometricNormal.normalize();

            return true;

        } else {
            return false;
        }
    }
};

// This class represents a number of triangle mesh objects with possible rigid motion.
// For speed, all the static objects are combined into a single triangle mesh, while all the
// moving objects are kept separate.
class particle_collision_detector {
    trimesh3 m_staticGeometry;
    boost::shared_ptr<trimesh3_kdtree> m_staticGeometryKdtree;
    std::vector<boost::shared_ptr<animated_particle_collision_mesh>> m_rigidGeometryObjects;

    boost::shared_ptr<primitive_kdtree<animated_particle_collision_mesh>> m_masterKdtree;

  public:
    particle_collision_detector() {}

    void clear() {
        m_staticGeometryKdtree.reset();
        m_staticGeometry.clear();
        m_rigidGeometryObjects.clear();
        m_masterKdtree.reset();
    }

    void swap( particle_collision_detector& rhs ) {
        m_staticGeometry.swap( rhs.m_staticGeometry );
        m_staticGeometryKdtree.swap( rhs.m_staticGeometryKdtree );
        m_rigidGeometryObjects.swap( rhs.m_rigidGeometryObjects );
        m_masterKdtree.swap( rhs.m_masterKdtree );
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
            // if( !add_rigid_object_to_existing( xform, mesh ) ) {
            m_rigidGeometryObjects.push_back(
                boost::shared_ptr<animated_particle_collision_mesh>( new animated_particle_collision_mesh() ) );
            m_rigidGeometryObjects.back()->set( xform, mesh );
            //}
        }
    }

    void add_rigid_object_with_swap( const motion_blurred_transform<float>& xform, trimesh3& mesh ) {
        // Try to add it to an existing rigid object, but if not add a new one
        if( mesh.face_count() > 0 ) {
            if( !add_rigid_object_to_existing( xform, mesh ) ) {
                m_rigidGeometryObjects.push_back(
                    boost::shared_ptr<animated_particle_collision_mesh>( new animated_particle_collision_mesh() ) );
                m_rigidGeometryObjects.back()->set_with_swap( xform, mesh );
            }
        }
    }

    void change_transform( const motion_blurred_transform<float>& xform, int i ) {
        if( i < 0 || i >= (int)m_rigidGeometryObjects.size() ) {
            throw std::runtime_error(
                "frantic::geometry::particle_collision_detector::change_transform() - The object index (" +
                boost::lexical_cast<std::string>( i ) + ") provided is out of bounds of the geometry objects array." );
        }
        m_rigidGeometryObjects[i]->set_transform( xform );
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
    void prepare_kdtrees( bool forceRefresh = false ) {
        if( forceRefresh || m_staticGeometryKdtree == 0 || m_masterKdtree == 0 ) {
            // The static geometry
            m_staticGeometryKdtree.reset( new trimesh3_kdtree( m_staticGeometry ) );

            // Destroy the old master kd-tree
            m_masterKdtree.reset();

            for( std::vector<boost::shared_ptr<animated_particle_collision_mesh>>::iterator i =
                     m_rigidGeometryObjects.begin();
                 i != m_rigidGeometryObjects.end(); ++i ) {
                ( *i )->prepare_kdtrees( forceRefresh );
            }

            // Build a new master kd-tree
            m_masterKdtree.reset( new primitive_kdtree<animated_particle_collision_mesh>( m_rigidGeometryObjects ) );
        }
    }

    void prepare_master_kdtree() {
        m_masterKdtree.reset( new primitive_kdtree<animated_particle_collision_mesh>( m_rigidGeometryObjects ) );
    }

    bool empty() const { return m_staticGeometry.face_count() == 0 && m_rigidGeometryObjects.size() == 0; }

    bool collides_with_particle( const vector3f& start, const vector3f& end,
                                 frantic::geometry::raytrace_intersection& outIntersection ) const {
        ray3f ray( start, end - start );
        return collides_with_particle( ray, 0, 1, outIntersection );
    }

    bool collides_with_particle( const ray3f& ray, double tMin, double tMax,
                                 frantic::geometry::raytrace_intersection& outIntersection ) const {
        bool intersected = false;

        if( m_staticGeometry.face_count() > 0 ) {
            if( m_staticGeometryKdtree.get() == 0 )
                throw std::runtime_error(
                    "particle_collision_detector.collides_with_particle: Tried to intersect a ray with a "
                    "mesh without first initializing the kdtrees." );

            if( m_staticGeometryKdtree->intersect_ray( ray, tMin, tMax, outIntersection ) ) {
                tMax = outIntersection.distance;
                outIntersection.motionDuringTimeStep = vector3f();
                intersected = true;
            }
        }

        if( m_masterKdtree.get() != 0 ) {
            if( m_masterKdtree->intersect_ray( ray, tMin, tMax, outIntersection ) ) {
                intersected = true;
            }
        }

        return intersected;
    }

    void update_particle( frantic::particles::particle_step& ps ) {
        frantic::geometry::raytrace_intersection intersect;

        // handle the collision by updating the end position with the appropriate position
        if( collides_with_particle( ps.oldPosition(), ps.newPosition(), intersect ) ) {
            float dt = ( 1 - (float)intersect.distance ) * ps.time_step();

            // reconstruct the velocity vector of the intersect
            vector3f objMotionAfterCollision = intersect.motionDuringTimeStep * ( 1 - (float)intersect.distance );

            // retrieve the expected end position and the velocity
            vector3f& newVelocity = ps.newVelocity();
            vector3f& newPos = ps.newPosition();

            // find the penetrating vector
            vector3f v = ( newPos - ( objMotionAfterCollision + intersect.position ) );

            // find the end position from the reflected unit vector and expected end position
            newPos = newPos - 2 * vector3f::dot( v, intersect.geometricNormal ) * intersect.geometricNormal;

            // find the adjusted velocity
            newVelocity = ( newPos - intersect.position ) / dt;

            // finally set the new position and the new velocity of the particle
            ps.set( newPos, newVelocity );
        }
    }

    void update_particle( frantic::particles::particle_step& ps, frantic::geometry::raytrace_intersection& intersect ) {

        // handle the collision by updating the end position with the appropriate position
        if( collides_with_particle( ps.oldPosition(), ps.newPosition(), intersect ) ) {
            float dt = ( 1 - (float)intersect.distance ) * ps.time_step();

            // reconstruct the velocity vector of the intersect
            vector3f objMotionAfterCollision = intersect.motionDuringTimeStep * ( 1 - (float)intersect.distance );

            // retrieve the expected end position and the velocity
            vector3f& newVelocity = ps.newVelocity();
            vector3f& newPos = ps.newPosition();

            // find the penetrating vector
            vector3f v = ( newPos - ( objMotionAfterCollision + intersect.position ) );

            // find the end position from the reflected unit vector and expected end position
            newPos = newPos - 2 * vector3f::dot( v, intersect.geometricNormal ) * intersect.geometricNormal;

            // find the adjusted velocity
            newVelocity = ( newPos - intersect.position ) / dt;

            // finally set the new position and the new velocity of the particle
            ps.set( newPos, newVelocity );
        }
    }

    void print_statistics( std::ostream& out ) const {
        using namespace std;

        out << " ===== Particle-Geometry Collision Detector Statistics =====" << endl;
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
};

} // namespace geometry
} // namespace frantic
