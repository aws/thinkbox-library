// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

//#undef NDEBUG
//#include <assert.h>

#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/ray3f.hpp>

namespace frantic {
namespace volumetrics {
using frantic::graphics::boundbox3f;
using frantic::graphics::ray3f;
using frantic::graphics::size3;
using frantic::graphics::vector3;
using frantic::graphics::vector3f;

// This class can be used to
// 1) step through every intersection of any axis' faces in a voxel grid
// 2) step through every intersection of a specified axis' faces  in a voxel grid
// These two things cannot be done simultaneously.  If you call step_x, step_y, step_z, or step_axis,
// do not call step.  The class will end up in an inconsistent state if you do that.

class voxel_edge_stepper {
    vector3f m_direction;

    // This is the scale of the voxels versus world space (in which the distances are specified)
    float m_scale;

    // The current intersection point
    float m_isectDist;
    vector3f m_isectPoint;
    int m_isectAxis;

    // The values necessary to step along the x=c intersections
    int m_xX, m_xXStep;
    float m_xY, m_xZ, m_xYStep, m_xZStep;
    float m_xDist, m_xDistStep;

    // The values necessary to step along the y=c intersections
    int m_yY, m_yYStep;
    float m_yX, m_yZ, m_yXStep, m_yZStep;
    float m_yDist, m_yDistStep;

    // The values necessary to step along the z=c intersections
    int m_zZ, m_zZStep;
    float m_zX, m_zY, m_zXStep, m_zYStep;
    float m_zDist, m_zDistStep;

    int m_debug_lastOperation;

    void initialize_step_sizes() {
        // along the x=c intersections
        if( m_direction.x == 0 ) {
            m_xXStep = 0;
            m_xYStep = 0;
            m_xZStep = 0;
            m_xDistStep = 0;
        } else {
            if( m_direction.x > 0 )
                m_xXStep = 1;
            else
                m_xXStep = -1;
            m_xYStep = m_xXStep * m_direction.y / m_direction.x;
            m_xZStep = m_xXStep * m_direction.z / m_direction.x;
            m_xDistStep = sqrt( 1 + m_xYStep * m_xYStep + m_xZStep * m_xZStep );
        }

        // along the y=c intersections
        if( m_direction.y == 0 ) {
            m_yXStep = 0;
            m_yYStep = 0;
            m_yZStep = 0;
            m_yDistStep = 0;
        } else {
            if( m_direction.y > 0 )
                m_yYStep = 1;
            else
                m_yYStep = -1;
            m_yXStep = m_yYStep * m_direction.x / m_direction.y;
            m_yZStep = m_yYStep * m_direction.z / m_direction.y;
            m_yDistStep = sqrt( 1 + m_yXStep * m_yXStep + m_yZStep * m_yZStep );
        }

        // along the z=c intersections
        if( m_direction.z == 0 ) {
            m_zXStep = 0;
            m_zYStep = 0;
            m_zZStep = 0;
            m_zDistStep = 0;
        } else {
            if( m_direction.z > 0 )
                m_zZStep = 1;
            else
                m_zZStep = -1;
            m_zXStep = m_zZStep * m_direction.x / m_direction.z;
            m_zYStep = m_zZStep * m_direction.y / m_direction.z;
            m_zDistStep = sqrt( 1 + m_zXStep * m_zXStep + m_zYStep * m_zYStep );
        }

        m_debug_lastOperation = 0;
    }

  public:
    voxel_edge_stepper() {
        m_debug_lastOperation = -1;
        m_scale = 1;
    }

    voxel_edge_stepper( vector3f direction ) {
        m_debug_lastOperation = -1;
        initialize_direction( direction );
    }

    voxel_edge_stepper( const vector3f& position, const vector3f& direction ) {
        m_debug_lastOperation = -1;
        initialize_stepper( position, direction );
    }

    voxel_edge_stepper( const ray3f& ray ) {
        m_debug_lastOperation = -1;
        initialize_stepper( ray );
    }

    void initialize_direction( const vector3f& direction ) {
        if( direction.get_magnitude_squared() == 0 )
            throw std::runtime_error(
                "voxel_edge_stepper.initialize_direction: Cannot create a stepper with direction 0" );
        m_direction = direction;
        // Normalize the direction, and clamp any extreme axis directions to zero.
        m_direction.normalize();
        if( fabs( m_direction.x ) < 0.0001f )
            m_direction.x = 0;
        if( fabs( m_direction.y ) < 0.0001f )
            m_direction.y = 0;
        if( fabs( m_direction.z ) < 0.0001f )
            m_direction.z = 0;
        initialize_step_sizes();
        m_scale = 1;
    }

    void initialize_stepper( const vector3f& position, const vector3f& direction ) {
        initialize_direction( direction );
        intersect_from( position );
    }

    void initialize_stepper( ray3f ray ) {
        initialize_direction( ray.direction() );
        intersect_from( ray.origin() );
    }

    void intersect_from( vector3f pos, float distance = 0 ) {
        m_isectAxis = 0;

        // If it's starting on  one of the axes, then detect and indicate that
        if( fabs( pos.x - floor( pos.x + 0.5f ) ) < 0.00001f * fabs( pos.x ) || fabs( pos.x ) < 0.00001f ) {
            pos.x = floor( pos.x + 0.5f );
            m_isectAxis |= 1;
        }
        if( fabs( pos.y - floor( pos.y + 0.5f ) ) < 0.00001f * fabs( pos.y ) || fabs( pos.y ) < 0.00001f ) {
            pos.y = floor( pos.y + 0.5f );
            m_isectAxis |= 2;
        }
        if( fabs( pos.z - floor( pos.z + 0.5f ) ) < 0.00001f * fabs( pos.z ) || fabs( pos.z ) < 0.00001f ) {
            pos.z = floor( pos.z + 0.5f );
            m_isectAxis |= 4;
        }

        // Starting from an arbitrary point
        m_isectPoint = pos;
        m_isectDist = distance;

        // along the x=c intersections
        if( m_xXStep == -1 ) {
            m_xX = (int)floor( pos.x );
            if( m_isectAxis & 1 )
                m_xX += m_xXStep;
            float delta = pos.x - m_xX;
            m_xY = pos.y + delta * m_xYStep;
            m_xZ = pos.z + delta * m_xZStep;
            m_xDist = distance + delta * m_xDistStep;
        } else if( m_xXStep == 1 ) {
            m_xX = (int)ceil( pos.x );
            if( m_isectAxis & 1 )
                m_xX += m_xXStep;
            float delta = m_xX - pos.x;
            m_xY = pos.y + delta * m_xYStep;
            m_xZ = pos.z + delta * m_xZStep;
            m_xDist = distance + delta * m_xDistStep;
        } else {
            m_xX = ( std::numeric_limits<int>::max )();
            m_xY = m_xZ = m_xDist = ( std::numeric_limits<float>::max )();
        }

        // along the y=c intersections
        if( m_yYStep == -1 ) {
            m_yY = (int)floor( pos.y );
            if( m_isectAxis & 2 )
                m_yY += m_yYStep;
            float delta = ( pos.y - m_yY );
            m_yX = pos.x + delta * m_yXStep;
            m_yZ = pos.z + delta * m_yZStep;
            m_yDist = distance + delta * m_yDistStep;
        } else if( m_yYStep == 1 ) {
            m_yY = (int)ceil( pos.y );
            if( m_isectAxis & 2 )
                m_yY += m_yYStep;
            float delta = ( m_yY - pos.y );
            m_yX = pos.x + delta * m_yXStep;
            m_yZ = pos.z + delta * m_yZStep;
            m_yDist = distance + delta * m_yDistStep;
        } else {
            m_yY = ( std::numeric_limits<int>::max )();
            m_yX = m_yZ = m_yDist = ( std::numeric_limits<float>::max )();
        }

        // along the z=c intersections
        if( m_zZStep == -1 ) {
            m_zZ = (int)floor( pos.z );
            if( m_isectAxis & 4 )
                m_zZ += m_zZStep;
            float delta = ( pos.z - m_zZ );
            m_zX = pos.x + delta * m_zXStep;
            m_zY = pos.y + delta * m_zYStep;
            m_zDist = distance + delta * m_zDistStep;
        } else if( m_zZStep == 1 ) {
            m_zZ = (int)ceil( pos.z );
            if( m_isectAxis & 4 )
                m_zZ += m_zZStep;
            float delta = ( m_zZ - pos.z );
            m_zX = pos.x + delta * m_zXStep;
            m_zY = pos.y + delta * m_zYStep;
            m_zDist = distance + delta * m_zDistStep;
        } else {
            m_zZ = ( std::numeric_limits<int>::max )();
            m_zX = m_zY = m_zDist = ( std::numeric_limits<float>::max )();
        }
        m_debug_lastOperation = 1;
    }

    // Sets up a ves in a child voxel space, with the following properties:
    //  - The distance() values are in the same space as the current (parent) ves
    void initialize_child_ves( voxel_edge_stepper& ves, size3 childVoxels ) {
        if( !childVoxels.is_cube() ) {
            std::string message =
                "voxel_edge_stepper.initialize_child_ves: The child voxel space must be a cube. " + childVoxels.str();
            std::cerr << message << std::endl;
            throw std::runtime_error( message );
        }
        // Copy the direction and step sizes, but scale the distances
        ves.m_direction = m_direction;

        ves.m_xXStep = m_xXStep;
        ves.m_xYStep = m_xYStep;
        ves.m_xZStep = m_xZStep;
        ves.m_xDistStep = m_xDistStep / childVoxels.xsize();

        ves.m_yXStep = m_yXStep;
        ves.m_yYStep = m_yYStep;
        ves.m_yZStep = m_yZStep;
        ves.m_yDistStep = m_yDistStep / childVoxels.xsize();

        ves.m_zXStep = m_zXStep;
        ves.m_zYStep = m_zYStep;
        ves.m_zZStep = m_zZStep;
        ves.m_zDistStep = m_zDistStep / childVoxels.xsize();

        // Copy the current intersection
        vector3 destVoxel = isect_voxel_being_entered();
        boundbox3f currentVoxel = boundbox3f( destVoxel, destVoxel + vector3( 1 ) );
        ves.intersect_from( currentVoxel.get_voxel_coord( m_isectPoint, childVoxels ), m_isectDist );
        // Start the child hitting the same axis direction (the detection in instersect_from is good enough)
        //		ves.m_isectAxis = m_isectAxis;

        ves.m_scale = m_scale / childVoxels.xsize();

        m_debug_lastOperation = 2;
    }

    // For a child ves which has already been initialized with initialize_child_ves,
    // this moves the child ves to start from the current position of the parent ves
    void move_child_ves( voxel_edge_stepper& ves, size3 childVoxels ) const {
        // Copy the current intersection
        vector3 destVoxel = isect_voxel_being_entered();
        boundbox3f currentVoxelBox = boundbox3f( destVoxel, destVoxel + vector3( 1 ) );
        vector3f childVoxelCoord = currentVoxelBox.get_voxel_coord( m_isectPoint, childVoxels );
        ves.intersect_from( childVoxelCoord, m_isectDist );
        // Start the child hitting the same axis direction (the detection in instersect_from is good enough)
        //		ves.m_isectAxis = m_isectAxis;
    }

    void apply_offset( vector3 offset ) {
        m_xX += offset.x;
        m_xY += offset.y;
        m_xZ += offset.z;

        m_yX += offset.x;
        m_yY += offset.y;
        m_yZ += offset.z;

        m_zX += offset.x;
        m_zY += offset.y;
        m_zZ += offset.z;

        m_isectPoint.x += offset.x;
        m_isectPoint.y += offset.y;
        m_isectPoint.z += offset.z;

        m_debug_lastOperation = 4;
    }

  private:
    // Steps to the next x plane intersection
    void step_x( bool unionIntersections ) {
        if( unionIntersections ) {
            m_isectAxis |= 1;
            // NOTE: the following assignments cause the distances to deviate from the actual
            // distance a lot.  The increase the stability of the integration, though.
            m_isectPoint.x = (float)m_xX;
            // To keep things in synch, copy the intersection data back over the x data
            m_xY = m_isectPoint.y;
            m_xZ = m_isectPoint.z;
            m_xDist = m_isectDist;
        } else {
            m_isectAxis = 1;
            m_isectPoint.set( (float)m_xX, m_xY, m_xZ );
            m_isectDist = m_xDist;
        }

        m_xX += m_xXStep;
        m_xY += m_xYStep;
        m_xZ += m_xZStep;
        m_xDist += m_xDistStep;
    }

    // Steps to the next y plane intersection
    void step_y( bool unionIntersections ) {
        if( unionIntersections ) {
            m_isectAxis |= 2;
            // NOTE: the following assignments cause the distances to deviate from the actual
            // distance a lot.  The increase the stability of the integration, though.
            m_isectPoint.y = (float)m_yY;
            // To keep things in synch, copy the intersection data back over the y data
            m_yX = m_isectPoint.x;
            m_yZ = m_isectPoint.z;
            m_yDist = m_isectDist;
        } else {
            m_isectAxis = 2;
            m_isectPoint.set( m_yX, (float)m_yY, m_yZ );
            m_isectDist = m_yDist;
        }

        m_yX += m_yXStep;
        m_yY += m_yYStep;
        m_yZ += m_yZStep;
        m_yDist += m_yDistStep;
    }

    // Steps to the next z plane intersection
    void step_z( bool unionIntersections ) {
        if( unionIntersections ) {
            m_isectAxis |= 4;
            // NOTE: the following assignments cause the distances to deviate from the actual
            // distance a lot.  The increase the stability of the integration, though.
            m_isectPoint.z = (float)m_zZ;
            // To keep things in synch, copy the intersection data back over the z data
            m_zX = m_isectPoint.x;
            m_zY = m_isectPoint.y;
            m_zDist = m_isectDist;
        } else {
            m_isectAxis = 4;
            m_isectPoint.set( m_zX, m_zY, (float)m_zZ );
            m_isectDist = m_zDist;
        }

        m_zX += m_zXStep;
        m_zY += m_zYStep;
        m_zZ += m_zZStep;
        m_zDist += m_zDistStep;
    }

  public:
    void step_axis( int axis, bool unionIntersections = false ) {
        // Save the current distance
        float savedDist = m_isectDist;

        switch( axis ) {
        case 0:
            step_x( unionIntersections );
            break;
        case 1:
            step_y( unionIntersections );
            break;
        case 2:
            step_z( unionIntersections );
            break;
        default: {
            std::string message = "voxel_edge_stepper.step_axis: Attempted to step along invalid axis, " +
                                  boost::lexical_cast<std::string>( axis );
            std::cerr << message << std::endl;
            throw std::runtime_error( message );
        }
        }
        // Within a small tolerance, multiple intersections happening at the same position
        // are skipped
        if( ( m_isectDist - savedDist ) <= 0.00001f * m_isectDist ) {
            step_axis( axis, true );
        }
    }

    // Take a step
    void step() {
        // assert( m_direction.get_magnitude_squared() > 0 );

        // Find the closest next intersection, assign it to the intersection point,
        // and increment that one.  If multiple intersections happen within a threshold,
        // we collapse those together.
        // NOTE: Collapsing multiple intersections used to be done by comparing distances.
        //       This turned out to be inaccurate in a number of cases, so has been converted
        //       to compare the actual relevant axis component.
        if( m_xDist <= m_yDist && m_xDist <= m_zDist ) {
            step_x( false );
            // Determine the threshold within which to merge multiple intersections
            // float stepUnionThreshold = std::max( 0.0005f * m_isectDist, 0.000001f );
            // if( m_yDist - m_isectDist < stepUnionThreshold )
            if( fabs( m_yY - m_isectPoint.y ) < ( std::max )( 0.00001f, (float)fabs( 0.00001f * m_isectPoint.y ) ) )
                step_y( true );
            // if( m_zDist - m_isectDist < stepUnionThreshold )
            if( fabs( m_zZ - m_isectPoint.z ) < ( std::max )( 0.00001f, (float)fabs( 0.00001f * m_isectPoint.z ) ) )
                step_z( true );
        } else if( m_yDist <= m_zDist ) {
            step_y( false );
            // Determine the threshold within which to merge multiple intersections
            // float stepUnionThreshold = std::max( 0.0005f * m_isectDist, 0.000001f );
            // if( m_xDist - m_isectDist < stepUnionThreshold )
            if( fabs( m_xX - m_isectPoint.x ) < ( std::max )( 0.00001f, (float)fabs( 0.00001f * m_isectPoint.x ) ) )
                step_x( true );
            // if( m_zDist - m_isectDist < stepUnionThreshold )
            if( fabs( m_zZ - m_isectPoint.z ) < ( std::max )( 0.00001f, (float)fabs( 0.00001f * m_isectPoint.z ) ) )
                step_z( true );
        } else {
            step_z( false );
            // Determine the threshold within which to merge multiple intersections
            // float stepUnionThreshold = std::max( 0.0005f * m_isectDist, 0.000001f );
            // if( m_xDist - m_isectDist < stepUnionThreshold )
            if( fabs( m_xX - m_isectPoint.x ) < ( std::max )( 0.00001f, (float)fabs( 0.00001f * m_isectPoint.x ) ) )
                step_x( true );
            // if( m_yDist - m_isectDist < stepUnionThreshold )
            if( fabs( m_yY - m_isectPoint.y ) < ( std::max )( 0.00001f, (float)fabs( 0.00001f * m_isectPoint.y ) ) )
                step_y( true );
        }

        m_debug_lastOperation = 5;
    }

    vector3f direction() const { return m_direction; }

    // Returns the square of the distance from the given point to the line
    // segment being traversed
    float distance_squared_from_point( vector3f p ) const {
        vector3f delta = p - m_isectPoint;
        delta -= vector3f::dot( delta, m_direction ) * m_direction;
        //		std::cerr << "p: " << p << " m_isectPoint: " << m_isectPoint << " vec: " << delta << std::endl;
        return delta.get_magnitude_squared();
    }

    // Returns which voxel is being entered.  The voxel that is returned
    // has centered at result + vector3f(0.5f)
    vector3 isect_voxel_being_entered() const {
        vector3 result( (int)floor( m_isectPoint.x ), (int)floor( m_isectPoint.y ), (int)floor( m_isectPoint.z ) );

        if( m_isectAxis & 1 && m_xXStep == -1 )
            --result.x;
        if( m_isectAxis & 2 && m_yYStep == -1 )
            --result.y;
        if( m_isectAxis & 4 && m_zZStep == -1 )
            --result.z;

        return result;
    }

    vector3 isect_face_corner() const {
        return vector3( (int)floor( m_isectPoint.x ), (int)floor( m_isectPoint.y ), (int)floor( m_isectPoint.z ) );
    }

    vector3f isect_position() const { return m_isectPoint; }

    int isect_axis() const { return m_isectAxis; }

    float isect_distance() const { return m_isectDist; }

    vector3f position_at( float distance ) { return m_isectPoint + m_scale * ( distance - m_isectDist ) * m_direction; }

    void debug_print( std::ostream& out ) const {
        using namespace std;

        out << "begin voxel_edge_stepper" << std::endl;
        out << "Direction: " << m_direction << std::endl;
        out << "X Next Isect: " << m_xX << ", " << m_xY << ", " << m_xZ << std::endl;
        out << "X Step: " << m_xXStep << ", " << m_xYStep << ", " << m_xZStep << std::endl;
        out << "X Dist: " << m_xDist << ", X Dist Step: " << m_xDistStep << std::endl;
        out << "Y Next Isect: " << m_yX << ", " << m_yY << ", " << m_yZ << std::endl;
        out << "Y Step: " << m_yXStep << ", " << m_yYStep << ", " << m_yZStep << std::endl;
        out << "Y Dist: " << m_yDist << ", Y Dist Step: " << m_yDistStep << std::endl;
        out << "Z Next Isect: " << m_zX << ", " << m_zY << ", " << m_zZ << std::endl;
        out << "Z Step: " << m_zXStep << ", " << m_zYStep << ", " << m_zZStep << std::endl;
        out << "Z Dist: " << m_zDist << ", Z Dist Step: " << m_zDistStep << std::endl;
        out << "Current Intersection Point: " << m_isectPoint << std::endl;
        out << "Current Intersection Distance: " << m_isectDist << std::endl;
        out << "Current Intersection Axis: " << m_isectAxis << std::endl;
        out << "Voxel Being Entered: " << isect_voxel_being_entered() << std::endl;
        out << "Last operation: ";
        switch( m_debug_lastOperation ) {
        case 0:
            out << "initialize_step_sizes";
            break;
        case 1:
            out << "intersect_from";
            break;
        case 2:
            out << "initialize_child_ves";
            break;
        case 3:
            out << "move_child_ves";
            break;
        case 4:
            out << "apply_offset";
            break;
        case 5:
            out << "step";
            break;
        default:
            out << "unknown";
            break;
        }
        out << std::endl;

        out << "end voxel_edge_stepper" << std::endl;
    }
};

} // namespace volumetrics
} // namespace frantic
