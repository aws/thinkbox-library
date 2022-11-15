// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/logic/tribool.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <frantic/geometry/generic_mesh_kdtree.hpp>
#include <frantic/geometry/generic_mesh_kdtree_impl.hpp>
#include <frantic/geometry/trimesh3.hpp>

#include "utilities/mesh_generators.hpp"

#include "gtest-helper.h"

struct reporting_mesh_traits {
    struct fake_result {
        double distance;
    };

    typedef frantic::geometry::trimesh3 mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_type_const_ptr;
    typedef fake_result raytrace_result;
    typedef fake_result nearest_point_result;

    static boost::tribool s_intersectRayIgnoreBackfaces;
    static boost::tribool s_nearestPointIgnoreBackfaces;

    inline static void reset_results() {
        s_intersectRayIgnoreBackfaces = boost::indeterminate;
        s_nearestPointIgnoreBackfaces = boost::indeterminate;
    }

    // Extracts the number of primitives( faces ) in the mesh.
    inline static unsigned get_count( const mesh_type& inst ) { return static_cast<unsigned>( inst.face_count() ); }

    // Extracts the bounding box of all primitives in the mesh.
    inline static frantic::graphics::boundbox3f get_bounds( const mesh_type& inst ) { return inst.compute_bound_box(); }

    // Extracts the bounds of a single primitive( ie.face ) after clipping the primitive to the supplied bounding box.
    inline static frantic::graphics::boundbox3f get_clipped_bounds( const mesh_type& inst, unsigned index,
                                                                    const frantic::graphics::boundbox3f& bounds ) {
        const frantic::graphics::vector3 face = inst.get_face( index );
        frantic::graphics::boundbox3f result;
        bounds.intersect_with_triangle( inst.get_vertex( face.x ), inst.get_vertex( face.y ), inst.get_vertex( face.z ),
                                        result );
        return result;
    }

    inline static bool intersect_ray( const mesh_type& /*inst*/, unsigned /*index*/,
                                      const frantic::graphics::ray3f& /*ray*/, double /*tMin*/, double /*tMax*/,
                                      raytrace_result& /*outIntersection*/, bool ignoreBackfaces ) {
        s_intersectRayIgnoreBackfaces = ignoreBackfaces;
        return true;
    }

    // Interestects the given ray with the specified primitive.Returns true if an intersection occurred in the( tMin,
    // tMax ) range.If returning true, then 'outIntersection' is updated with the relevant intersection data.
    inline static bool find_nearest_point( const mesh_type& /*inst*/, unsigned /*index*/,
                                           const frantic::graphics::vector3f& /*point*/, double /*maxDistance*/,
                                           nearest_point_result& /*outNearestPoint*/, bool ignoreBackfaces ) {
        s_nearestPointIgnoreBackfaces = ignoreBackfaces;
        return true;
    }
};

boost::tribool reporting_mesh_traits::s_intersectRayIgnoreBackfaces = boost::indeterminate;
boost::tribool reporting_mesh_traits::s_nearestPointIgnoreBackfaces = boost::indeterminate;

TEST( GenericMeshKDTree, BackfaceIgnorancePropagationIntersectRay ) {
    boost::shared_ptr<frantic::geometry::trimesh3> testMesh = boost::make_shared<frantic::geometry::trimesh3>();
    make_cube_mesh( *testMesh );

    frantic::geometry::generic_mesh_kdtree<reporting_mesh_traits> tree( testMesh );

    frantic::graphics::ray3f ray( frantic::graphics::vector3f( 0.0f, 0.0f, 2.0f ),
                                  frantic::graphics::vector3f( 0.0f, 0.0f, -1.0f ) );
    reporting_mesh_traits::fake_result raytraceResult;

    reporting_mesh_traits::reset_results();

    tree.intersect_ray( ray, 0.0, std::numeric_limits<double>::max(), raytraceResult, true );

    EXPECT_FALSE( boost::indeterminate( reporting_mesh_traits::s_intersectRayIgnoreBackfaces ) );
    EXPECT_EQ( true, reporting_mesh_traits::s_intersectRayIgnoreBackfaces );

    reporting_mesh_traits::reset_results();

    tree.intersect_ray( ray, 0.0, std::numeric_limits<double>::max(), raytraceResult, false );

    EXPECT_FALSE( boost::indeterminate( reporting_mesh_traits::s_intersectRayIgnoreBackfaces ) );
    EXPECT_EQ( false, reporting_mesh_traits::s_intersectRayIgnoreBackfaces );
}

TEST( GenericMeshKDTree, BackfaceIgnorancePropagationNearestPoint ) {
    boost::shared_ptr<frantic::geometry::trimesh3> testMesh = boost::make_shared<frantic::geometry::trimesh3>();
    make_cube_mesh( *testMesh );

    frantic::geometry::generic_mesh_kdtree<reporting_mesh_traits> tree( testMesh );
    reporting_mesh_traits::fake_result nearestPointResult;

    reporting_mesh_traits::reset_results();

    tree.find_nearest_point( frantic::graphics::vector3f( 0.0f, 0.0f, 2.0f ), std::numeric_limits<double>::max(),
                             nearestPointResult, true );

    EXPECT_FALSE( boost::indeterminate( reporting_mesh_traits::s_nearestPointIgnoreBackfaces ) );
    EXPECT_EQ( true, reporting_mesh_traits::s_nearestPointIgnoreBackfaces );

    reporting_mesh_traits::reset_results();

    tree.find_nearest_point( frantic::graphics::vector3f( 0.0f, 0.0f, 2.0f ), std::numeric_limits<double>::max(),
                             nearestPointResult, false );

    EXPECT_FALSE( boost::indeterminate( reporting_mesh_traits::s_nearestPointIgnoreBackfaces ) );
    EXPECT_EQ( false, reporting_mesh_traits::s_nearestPointIgnoreBackfaces );
}
