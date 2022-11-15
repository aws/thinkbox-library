// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// boundbox3f.cpp
//
// 3D AABB functions.

// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/graphics_utils.hpp>
#include <frantic/graphics/plane3f.hpp>
#include <frantic/graphics/ray3f.hpp>

namespace frantic {
namespace graphics {

// TODO: Maybe move this class somewhere else? This isn't being used anywhere?
#if 0
template<class FloatType>
class triangle_intersector {
	typedef FloatType float_type;
	typedef vector3f<FloatType> vector3f_type;
	plane3t<FloatType> m_plane;
	float_type m_edge0axis0, m_edge0axis1, m_edge1axis0, m_edge1axis1;
	vector3f_type m_vert0; //, m_edge0, m_edge1;
	char m_barycentric0Axis, m_barycentric1Axis;

public:
	triangle_intersector( const vector3f_type& vert0, const vector3f_type& vert1, const vector3f_type& vert2 )
		: m_vert0(vert0), m_plane( plane3t<FloatType>::from_triangle(vert0,vert1,vert2) )
	{
		m_plane = plane3t<FloatType>::from_triangle(vert0, vert1, vert2);
		float barycentricInverseDeterminant;
		compute_barycentric_helpers(vert0, vert1, vert2, m_barycentric0Axis, m_barycentric1Axis, barycentricInverseDeterminant);
		vector3f_type edge0 = vert1-vert0, edge1 = vert2-vert0;
		m_edge0axis0 = edge0[m_barycentric0Axis] * barycentricInverseDeterminant;
		m_edge0axis1 = edge0[m_barycentric1Axis] * barycentricInverseDeterminant;
		m_edge1axis0 = edge1[m_barycentric0Axis] * barycentricInverseDeterminant;
		m_edge1axis1 = edge1[m_barycentric1Axis] * barycentricInverseDeterminant;
	}

	bool intersect_with_line_segment( const vector3f_type& start, const vector3f_type& end, vector3f_type& outIntersection ) const {
		if( m_plane.get_intersection( start, end, outIntersection ) ) {
			// This is a tweaked version of vector3f::compute_barycentric_coordinates_with_helpers
			// In spite of this being smaller and simpler, it did not give a big performance boost over the original function call.
			vector3f_type relativePosition = outIntersection - m_vert0;

			// Solve the resulting 2x2 linear equation
			float_type barycentricB = m_edge1axis1*relativePosition[m_barycentric0Axis] - m_edge1axis0*relativePosition[m_barycentric1Axis];

			// Giving a little bit of leeway on the barycentric values
			// is really important for some applications of this function.
			// Most notably, the function which computes the intersection of
			// a box with a triangle, returning the bounding box of that
			// intersection, will fail if 0.f and 1.f are used instead of
			// these slightly expanded values.
			if( barycentricB < -0.00001f ) // TODO: Different tolerence for float vs double?
				return false;

			float barycentricC = m_edge0axis0*relativePosition[m_barycentric1Axis] - m_edge0axis1*relativePosition[m_barycentric0Axis];

			return (barycentricC >= -0.00001f && (barycentricB + barycentricC) <= 1.00001f);
		} else {
			return false;
		}
	}

};
#endif

// Triangle Box intersection test as described by Tomas Akenine-Moller in his 2001 paper, "Fast 3D Triangle-Box Overlap
// Testing." This uses the separating axis theorem to determine whether the triangle and box are intersecting or not.  I
// believe the 13 directions being tested are all the face normals of the Minkowski difference between the box and the
// triangle, but I'm not sure about that.
template <class FloatType>
inline bool is_intersecting_triangle_templ( const boundbox3t<FloatType>& box, vector3t<FloatType> vert0,
                                            vector3t<FloatType> vert1, vector3t<FloatType> vert2 ) {
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;

    vector3f_type boxCenter = box.center(), boxHalfWidth( box.xsize() * 0.5f, box.ysize() * 0.5f, box.zsize() * 0.5f );
    // Adjust the triangle so the box is centered at the origin.
    vert0 -= boxCenter;
    vert1 -= boxCenter;
    vert2 -= boxCenter;

    vector3f_type edge0 = vert1 - vert0, edge1 = vert2 - vert1, edge2 = vert0 - vert2;
    float_type p0, p1, p2, minProject, maxProject, boxProject;

    // First do the 9 tests of cross(axis[i], edge[i]) (Bullet 3 in Akenine-Moller's paper)
    // NOTE: In Akenine's example code, these nine cases are collapsed using macro magic.  I implemented this in
    //       expanded form to be able to actually understand what's going on.  We may want to do something like this
    //       here as well.
    // a_00 = cross([1,0,0], edge0) = [0, -edge0.z, edge0.y]
    //  3 triangle vert projections:
    //      p0 = dot(a_00, vert0) = vert1.y * vert0.z - vert1.z * vert0.y
    //      p1 = dot(a_00, vert1) = p0
    //      p2 = dot(a_00, vert2) = edge0.y * vert2.z - edge0.z * vert2.y
    p0 = -edge0.z * vert0.y + edge0.y * vert0.z;
    p2 = -edge0.z * vert2.y + edge0.y * vert2.z;
    if( p0 < p2 ) {
        minProject = p0;
        maxProject = p2;
    } else {
        minProject = p2;
        maxProject = p0;
    }
    boxProject = fabs( edge0.z ) * boxHalfWidth.y + fabs( edge0.y ) * boxHalfWidth.z;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_01 = cross([1,0,0], edge1) = [0, -edge1.z, edge1.y]
    //  3 triangle vert projections:
    //      p0 = dot(a_01, vert0) = edge1.y * vert0.z - edge1.z * vert0.y
    //      p1 = dot(a_01, vert1) = vert1.z * vert2.y - vert1.y * vert2.z
    //      p2 = dot(a_01, vert2) = p1
    p0 = -edge1.z * vert0.y + edge1.y * vert0.z;
    p1 = -edge1.z * vert1.y + edge1.y * vert1.z;
    if( p0 < p1 ) {
        minProject = p0;
        maxProject = p1;
    } else {
        minProject = p1;
        maxProject = p0;
    }
    boxProject = fabs( edge1.z ) * boxHalfWidth.y + fabs( edge1.y ) * boxHalfWidth.z;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_02 = cross([1,0,0], edge2) = [0, -edge2.z, edge2.y]
    //   Here, p2 = p0
    p0 = -edge2.z * vert0.y + edge2.y * vert0.z;
    p1 = -edge2.z * vert1.y + edge2.y * vert1.z;
    if( p0 < p1 ) {
        minProject = p0;
        maxProject = p1;
    } else {
        minProject = p1;
        maxProject = p0;
    }
    boxProject = fabs( edge2.z ) * boxHalfWidth.y + fabs( edge2.y ) * boxHalfWidth.z;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_10 = cross([0,1,0], edge0) = [edge0.z, 0, -edge0.x]
    //   Here, p1 = p0
    p0 = edge0.z * vert0.x - edge0.x * vert0.z;
    p2 = edge0.z * vert2.x - edge0.x * vert2.z;
    if( p0 < p2 ) {
        minProject = p0;
        maxProject = p2;
    } else {
        minProject = p2;
        maxProject = p0;
    }
    boxProject = fabs( edge0.z ) * boxHalfWidth.x + fabs( edge0.x ) * boxHalfWidth.z;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_11 = cross([0,1,0], edge1) = [edge1.z, 0, -edge1.x]
    //   Here, p2 = p1
    p0 = edge1.z * vert0.x - edge1.x * vert0.z;
    p1 = edge1.z * vert1.x - edge1.x * vert1.z;
    if( p0 < p1 ) {
        minProject = p0;
        maxProject = p1;
    } else {
        minProject = p1;
        maxProject = p0;
    }
    boxProject = fabs( edge1.z ) * boxHalfWidth.x + fabs( edge1.x ) * boxHalfWidth.z;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_12 = cross([0,1,0], edge2) = [edge2.z, 0, -edge2.x]
    //   Here, p2 = p0
    p0 = edge2.z * vert0.x - edge2.x * vert0.z;
    p1 = edge2.z * vert1.x - edge2.x * vert1.z;
    if( p0 < p1 ) {
        minProject = p0;
        maxProject = p1;
    } else {
        minProject = p1;
        maxProject = p0;
    }
    boxProject = fabs( edge2.z ) * boxHalfWidth.x + fabs( edge2.x ) * boxHalfWidth.z;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_20 = cross([0,0,1], edge0) = [-edge0.y, edge0.x, 0]
    //   Here, p1 = p0
    p0 = -edge0.y * vert0.x + edge0.x * vert0.y;
    p2 = -edge0.y * vert2.x + edge0.x * vert2.y;
    if( p0 < p2 ) {
        minProject = p0;
        maxProject = p2;
    } else {
        minProject = p2;
        maxProject = p0;
    }
    boxProject = fabs( edge0.y ) * boxHalfWidth.x + fabs( edge0.x ) * boxHalfWidth.y;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_21 = cross([0,0,1], edge1) = [-edge1.y, edge1.x, 0]
    //   Here, p2 = p1
    p0 = -edge1.y * vert0.x + edge1.x * vert0.y;
    p1 = -edge1.y * vert1.x + edge1.x * vert1.y;
    if( p0 < p1 ) {
        minProject = p0;
        maxProject = p1;
    } else {
        minProject = p1;
        maxProject = p0;
    }
    boxProject = fabs( edge1.y ) * boxHalfWidth.x + fabs( edge1.x ) * boxHalfWidth.y;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;
    // a_22 = cross([0,0,1], edge2) = [-edge2.y, edge2.x, 0]
    //   Here, p2 = p0
    p0 = -edge2.y * vert0.x + edge2.x * vert0.y;
    p1 = -edge2.y * vert1.x + edge2.x * vert1.y;
    if( p0 < p1 ) {
        minProject = p0;
        maxProject = p1;
    } else {
        minProject = p1;
        maxProject = p0;
    }
    boxProject = fabs( edge2.y ) * boxHalfWidth.x + fabs( edge2.x ) * boxHalfWidth.y;
    if( minProject > boxProject || maxProject < -boxProject )
        return false;

    // Second do the 3 tests of along the axis directions. (Bullet 1 in Akenine-Moller's paper)
    // This is the same as testing whether the AABB around the triangle intersects this bound box
    // e_0 = [1,0,0]
    minProject = maxProject = vert0.x;
    if( vert1.x < minProject )
        minProject = vert1.x;
    if( vert1.x > maxProject )
        maxProject = vert1.x;
    if( vert2.x < minProject )
        minProject = vert2.x;
    if( vert2.x > maxProject )
        maxProject = vert2.x;
    if( minProject > boxHalfWidth.x || maxProject < -boxHalfWidth.x )
        return false;
    // e_1 = [0,1,0]
    minProject = maxProject = vert0.y;
    if( vert1.y < minProject )
        minProject = vert1.y;
    if( vert1.y > maxProject )
        maxProject = vert1.y;
    if( vert2.y < minProject )
        minProject = vert2.y;
    if( vert2.y > maxProject )
        maxProject = vert2.y;
    if( minProject > boxHalfWidth.y || maxProject < -boxHalfWidth.y )
        return false;
    // e_2 = [0,0,1]
    minProject = maxProject = vert0.z;
    if( vert1.z < minProject )
        minProject = vert1.z;
    if( vert1.z > maxProject )
        maxProject = vert1.z;
    if( vert2.z < minProject )
        minProject = vert2.z;
    if( vert2.z > maxProject )
        maxProject = vert2.z;
    if( minProject > boxHalfWidth.z || maxProject < -boxHalfWidth.z )
        return false;

    // Finally do the test along the triangle normal. (Bullet 2 in Akenine-Moller's paper)
    // The normal calculation here is the only potential robustness problem in the algorithm.
    // NOTE: The example code used a more complex plane-box intersection.  I don't know why the same
    //       technique as is used in the 9 cases of bullet 2 shouldn't be used, so I'm using it.
    //       Sent Tomas an email about the idea.
    vector3f_type triNormal = vector3f_type::cross( edge0, edge1 );
    // Determine the interval of the triangle along the normal (the interval is only a single point)
    minProject = vector3f_type::dot( triNormal, vert0 );
    boxProject = fabs( triNormal.x ) * boxHalfWidth.x + fabs( triNormal.y ) * boxHalfWidth.y +
                 fabs( triNormal.z ) * boxHalfWidth.z;
    if( minProject > boxProject || minProject < -boxProject )
        return false;

    return true;
}

bool detail::is_intersecting_triangle( const boundbox3f& box, const vector3f& vert0, const vector3f& vert1,
                                       const vector3f& vert2 ) {
    return is_intersecting_triangle_templ<float>( box, vert0, vert1, vert2 );
}

bool detail::is_intersecting_triangle( const boundbox3fd& box, const vector3fd& vert0, const vector3fd& vert1,
                                       const vector3fd& vert2 ) {
    return is_intersecting_triangle_templ<double>( box, vert0, vert1, vert2 );
}

// This function finds the intersection of the triangle and the bound box, returning the result in another bound box.
// It is used by the kd-tree construction code to compute perfect splits.  An example of where a more naive algorithm
// would erroneously include a triangle when it shouldn't is as follows:
// |-|   /
// | |  /
// |-| /
//    /
//   /
//  /
//
template <class FloatType>
inline bool intersect_with_triangle_templ( const boundbox3t<FloatType>& box, const vector3t<FloatType>& vert0,
                                           const vector3t<FloatType>& vert1, const vector3t<FloatType>& vert2,
                                           boundbox3t<FloatType>& outIntersectedBounds ) {
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;

    int count[2];
    // A planar-convex polygon can produce at most 1 more vertex per plane intersection. 6 intersections w/ a triangle
    // means at most 9 verts
    vector3f_type storage[2][9];

    count[0] = 3;
    storage[0][0] = vert0;
    storage[0][1] = vert1;
    storage[0][2] = vert2;

    for( int dim = 0; dim < 3; ++dim ) {
        count[1] = 0;
        for( int i = 0; i < count[0]; ++i ) {
            const vector3f_type& v0 = storage[0][i];
            const vector3f_type& v1 = storage[0][( i + 1 ) % count[0]];

            if( v0[dim] == box.minimum()[dim] ) {
                storage[1][count[1]++] = v0;
            } else if( v0[dim] > box.minimum()[dim] ) {
                storage[1][count[1]++] = v0;
                if( v1[dim] < box.minimum()[dim] ) {
                    float_type t = ( box.minimum()[dim] - v0[dim] ) / ( v1[dim] - v0[dim] );
                    storage[1][count[1]++] = ( 1 - t ) * v0 + t * v1;
                }
            } else if( v1[dim] > box.minimum()[dim] ) {
                float_type t = ( box.minimum()[dim] - v0[dim] ) / ( v1[dim] - v0[dim] );
                storage[1][count[1]++] = ( 1 - t ) * v0 + t * v1;
            }
        }

        if( count[1] == 0 )
            return false;

        count[0] = 0;
        for( int i = 0; i < count[1]; ++i ) {
            const vector3f_type& v0 = storage[1][i];
            const vector3f_type& v1 = storage[1][( i + 1 ) % count[1]];

            if( v0[dim] == box.maximum()[dim] ) {
                storage[0][count[0]++] = v0;
            } else if( v0[dim] < box.maximum()[dim] ) {
                storage[0][count[0]++] = v0;
                if( v1[dim] > box.maximum()[dim] ) {
                    float_type t = ( box.maximum()[dim] - v0[dim] ) / ( v1[dim] - v0[dim] );
                    storage[0][count[0]++] = ( 1 - t ) * v0 + t * v1;
                }
            } else if( v1[dim] < box.maximum()[dim] ) {
                float_type t = ( box.maximum()[dim] - v0[dim] ) / ( v1[dim] - v0[dim] );
                storage[0][count[0]++] = ( 1 - t ) * v0 + t * v1;
            }
        }

        if( count[0] == 0 )
            return false;
    }

    outIntersectedBounds.set_to_empty();
    for( int i = 0; i < count[0]; ++i )
        outIntersectedBounds += storage[0][i];

    return true;
}

bool detail::intersect_with_triangle( const boundbox3f& box, const vector3f& vert0, const vector3f& vert1,
                                      const vector3f& vert2, boundbox3f& outIntersectedBounds ) {
    return intersect_with_triangle_templ<float>( box, vert0, vert1, vert2, outIntersectedBounds );
}
bool detail::intersect_with_triangle( const boundbox3fd& box, const vector3fd& vert0, const vector3fd& vert1,
                                      const vector3fd& vert2, boundbox3fd& outIntersectedBounds ) {
    return intersect_with_triangle_templ<double>( box, vert0, vert1, vert2, outIntersectedBounds );
}

// This intersects a line segment with the bounding box, and returns both the intersection points and the number of
// intersections (0, 1, or 2)
template <class FloatType>
inline int intersect_with_line_segment_templ( const boundbox3t<FloatType>& box, const vector3t<FloatType>& start,
                                              const vector3t<FloatType>& end, vector3t<FloatType>& outFirstIntersection,
                                              vector3t<FloatType>& outSecondIntersection ) {
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;

    ray3t<float_type> ray( start, end - start );
    int intersectionCount = 0;
    double t0 = 0, t1 = 0;
    vector3f_type entryNormal, exitNormal;
    if( ray.intersect_with_box( box.minimum(), box.maximum(), t0, t1, &entryNormal, &exitNormal ) ) {
        if( t0 >= 0 && t0 <= 1 ) {
            outFirstIntersection = ray.at( t0 );
            ++intersectionCount;
        }
        if( t1 >= 0 && t1 <= 1 ) {
            if( intersectionCount == 0 )
                outFirstIntersection = ray.at( t1 );
            else
                outSecondIntersection = ray.at( t1 );
            ++intersectionCount;
        }
    }
    return intersectionCount;
}

int detail::intersect_with_line_segment( const boundbox3f& box, const vector3f& start, const vector3f& end,
                                         vector3f& outFirstIntersection, vector3f& outSecondIntersection ) {
    return intersect_with_line_segment_templ<float>( box, start, end, outFirstIntersection, outSecondIntersection );
}
int detail::intersect_with_line_segment( const boundbox3fd& box, const vector3fd& start, const vector3fd& end,
                                         vector3fd& outFirstIntersection, vector3fd& outSecondIntersection ) {
    return intersect_with_line_segment_templ<double>( box, start, end, outFirstIntersection, outSecondIntersection );
}

template <class FloatType>
inline bool is_volume_intersecting_line_segment_templ( const boundbox3t<FloatType>& box,
                                                       const vector3t<FloatType>& start,
                                                       const vector3t<FloatType>& end ) {
    typedef FloatType float_type;
    typedef vector3t<FloatType> vector3f_type;

    if( box.is_empty() ) {
        return false;
    } else if( box.contains( start ) || box.contains( end ) ) {
        return true;
    } else if( ( start.x < box.minimum().x && end.x < box.minimum().x ) ||
               ( start.y < box.minimum().y && end.y < box.minimum().y ) ||
               ( start.z < box.minimum().z && end.z < box.minimum().z ) ||
               ( start.x > box.maximum().x && end.x > box.maximum().x ) ||
               ( start.y > box.maximum().y && end.y > box.maximum().y ) ||
               ( start.z > box.maximum().z && end.z > box.maximum().z ) ) {
        // A number of special cases of guaranteed non-intersection, which were very numerically
        // unstable in the test below.
        return false;
    } else {
        vector3f_type delta = end - start;
        float_type maximumOfNearDistances = 0, minimumOfFarDistances = 1;

        if( delta.x != 0 ) {
            if( delta.x > 0 ) {
                maximumOfNearDistances = ( box.minimum().x - start.x ) / delta.x;
                minimumOfFarDistances = ( box.maximum().x - start.x ) / delta.x;
            } else {
                maximumOfNearDistances = ( box.maximum().x - start.x ) / delta.x;
                minimumOfFarDistances = ( box.minimum().x - start.x ) / delta.x;
            }
        }

        if( delta.y != 0 ) {
            if( delta.y > 0 ) {
                maximumOfNearDistances =
                    ( std::max )( ( box.minimum().y - start.y ) / delta.y, maximumOfNearDistances );
                minimumOfFarDistances = ( std::min )( ( box.maximum().y - start.y ) / delta.y, minimumOfFarDistances );
            } else {
                maximumOfNearDistances =
                    ( std::max )( ( box.maximum().y - start.y ) / delta.y, maximumOfNearDistances );
                minimumOfFarDistances = ( std::min )( ( box.minimum().y - start.y ) / delta.y, minimumOfFarDistances );
            }
        }

        if( delta.z != 0 ) {
            if( delta.z > 0 ) {
                maximumOfNearDistances =
                    ( std::max )( ( box.minimum().z - start.z ) / delta.z, maximumOfNearDistances );
                minimumOfFarDistances = ( std::min )( ( box.maximum().z - start.z ) / delta.z, minimumOfFarDistances );
            } else {
                maximumOfNearDistances =
                    ( std::max )( ( box.maximum().z - start.z ) / delta.z, maximumOfNearDistances );
                minimumOfFarDistances = ( std::min )( ( box.minimum().z - start.z ) / delta.z, minimumOfFarDistances );
            }
        }

        if( maximumOfNearDistances >= 1 || minimumOfFarDistances <= 0 ) {
            return false;
        }

        return ( maximumOfNearDistances <= minimumOfFarDistances );
    }
}

bool detail::is_volume_intersecting_line_segment( const boundbox3f& box, const vector3f& start, const vector3f& end ) {
    return is_volume_intersecting_line_segment_templ<float>( box, start, end );
}
bool detail::is_volume_intersecting_line_segment( const boundbox3fd& box, const vector3fd& start,
                                                  const vector3fd& end ) {
    return is_volume_intersecting_line_segment_templ<double>( box, start, end );
}

// takes a point and performs a modulus operation to place the point within the boundbox
// the return container will then have any points that are modulo permutations of the original point
// that are within the given threshold of the faces of the boundbox
// The modded original point is also returned and will be the last point in the container
// This could be used to create cyclic boundary conditions
template <class FloatType>
inline void compute_mod_boundary_points_templ( const boundbox3t<FloatType>& box, const vector3t<FloatType>& point,
                                               FloatType boundaryThreshold,
                                               std::vector<vector3t<FloatType>>& outBoundaryPoints ) {
    typedef vector3t<FloatType> vector3f_type;

    bool testPointFlag[6];
    int aindex, bindex, cindex;
    vector3f_type aTest, bTest, cTest;
    vector3f_type difference = box.maximum() - box.minimum();

    // mod the point into the boundbox
    vector3f_type finalPosition = box.mod_point( point );

    // std::cout << "\nTest Point: " << finalPosition << std::endl;
    //  find the distances to faces and set if the flag if the point is close to the face
    testPointFlag[0] = fabs( finalPosition.x - box.minimum().x ) < boundaryThreshold;
    testPointFlag[1] = fabs( finalPosition.x - box.maximum().x ) < boundaryThreshold;

    testPointFlag[2] = fabs( finalPosition.y - box.minimum().y ) < boundaryThreshold;
    testPointFlag[3] = fabs( finalPosition.y - box.maximum().y ) < boundaryThreshold;

    testPointFlag[4] = fabs( finalPosition.z - box.minimum().z ) < boundaryThreshold;
    testPointFlag[5] = fabs( finalPosition.z - box.maximum().z ) < boundaryThreshold;

    // find all permutations of the locations that close to each of the six faces
    for( int a = 0; a < 6; ++a ) {

        if( testPointFlag[a] ) {
            // std::cout << "a: " << a << ", " << testPointFlag[a] << std::endl;

            aTest = finalPosition;

            aindex = a / 2;
            aTest[aindex] += ( a % 2 == 0 ) ? difference[aindex] : -difference[aindex];
            // std::cout << "increment value: " << ( (a%2==0) ? difference[aindex] : -difference[aindex] ) << std::endl;
            // std::cout << " a: insert " << aTest << std::endl;
            outBoundaryPoints.push_back( aTest );

            for( int b = ( a % 2 == 0 ) ? a + 2 : a + 1; b < 6; ++b ) {
                if( testPointFlag[b] ) {
                    bTest = aTest;
                    bindex = b / 2;
                    bTest[bindex] += ( b % 2 == 0 ) ? difference[bindex] : -difference[bindex];
                    // std::cout << " b: insert " << bTest << std::endl;
                    outBoundaryPoints.push_back( bTest );

                    for( int c = ( b % 2 == 0 ) ? b + 2 : b + 1; c < 6; ++c ) {
                        if( testPointFlag[c] ) {
                            cTest = bTest;
                            cindex = c / 2;
                            cTest[cindex] += ( c % 2 == 0 ) ? difference[cindex] : -difference[cindex];
                            // cout << "(a,b,c) (" << a << ", " << b << ", " << c << ")\n";
                            // std::cout << " c: insert " << cTest << std::endl;
                            outBoundaryPoints.push_back( cTest );
                        }
                    }
                }
            }
        }
    }

    // add the actually true position
    outBoundaryPoints.push_back( finalPosition );
}

void detail::compute_mod_boundary_points( const boundbox3f& box, const vector3f& point, float boundaryThreshold,
                                          std::vector<vector3f>& outBoundaryPoints ) {
    compute_mod_boundary_points_templ<float>( box, point, boundaryThreshold, outBoundaryPoints );
}
void detail::compute_mod_boundary_points( const boundbox3fd& box, const vector3fd& point, double boundaryThreshold,
                                          std::vector<vector3fd>& outBoundaryPoints ) {
    compute_mod_boundary_points_templ<double>( box, point, boundaryThreshold, outBoundaryPoints );
}
} // namespace graphics
} // namespace frantic
