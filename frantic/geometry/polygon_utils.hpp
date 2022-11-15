// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <frantic/graphics/graphics_utils.hpp>
#include <frantic/graphics/plane3f.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <boost/range/irange.hpp>

#include <list>
#include <utility>
#include <vector>

namespace frantic {
namespace geometry {

namespace boundary_relation {
enum boundary_relation {
    outside = 0,
    boundary = 1,
    inside = 2,
};
}

/**
 * Point in polygon test, using horizontal scan technique (integer version)
 *
 * \param polygon list of polygon vertices
 * \param point location to test
 * \return true if point in polygon, false otherwise
 */
boundary_relation::boundary_relation point_in_polygon( const std::vector<frantic::graphics2d::vector2>& polygon,
                                                       const frantic::graphics2d::vector2& point );

/**
 * Point in polygon test, using horizontal scan technique (float version)
 *
 * \param polygon list of polygon vertices
 * \param point location to test
 * \return true if point in polygon, false otherwise
 */
boundary_relation::boundary_relation point_in_polygon( const std::vector<frantic::graphics2d::vector2f>& polygon,
                                                       const frantic::graphics2d::vector2f& point,
                                                       float tolerance = 1.0e-5f );

/**
 * Find the minimum area rectangle which encloses all of the input points. Note that the output will not (in general) be
 * axis-aligned
 *
 * @param points the point set to enclose
 * @param outRect store the resulting 4 points found
 */
void minimum_enclosing_rectange( const std::vector<frantic::graphics2d::vector2f>& points,
                                 std::vector<frantic::graphics2d::vector2f>& outRect );

/**
 * Return the area of the given polygon
 *
 * @tparam InputIterator A read-only forward iterator which dereferences to a vector2f
 * @param begin iterator to the start of the polygon sequence
 * @param end iterator to the one-past-the-end of the polygon sequence
 * @param isSigned (optional) switch if this algorithm should report the 'signed' area or not.
 *   The signed area is positive if the polygon appears counter-clockwise, negative otherwise.
 * @return the area of the polygon
 */
template <class InputIterator>
inline double polygon2_area( InputIterator begin, InputIterator end, bool isSigned = false ) {
    typedef typename std::iterator_traits<InputIterator>::value_type vector_type;

    double sumCross = 0.0;
    InputIterator it = begin;

    while( it != end ) {
        const vector_type p0 = *it;

        ++it;

        const vector_type p1 = ( it == end ) ? *begin : *it;

        sumCross += double( p0[0] ) * double( p1[1] ) - double( p1[0] ) * double( p0[1] );
    }

    sumCross /= 2.0;

    return isSigned ? sumCross : frantic::math::get_absolute( sumCross );
}

template <class VectorType>
inline double triangle2_area( const VectorType& p0, const VectorType& p1, const VectorType& p2,
                              bool isSigned = false ) {
    const VectorType p[3] = { p0, p1, p2 };
    return polygon2_area( p, p + 3, isSigned );
}

template <class VectorType>
inline frantic::graphics::vector3fd get_barycentric_coordinates( const VectorType& p, const VectorType& p0,
                                                                 const VectorType& p1, const VectorType& p2 ) {
    using namespace frantic::graphics2d;

    VectorType edge0 = p1 - p0;
    VectorType edge1 = p2 - p0;
    VectorType diff = p - p0;

    double d00 = VectorType::dot_double( edge0, edge0 );
    double d01 = VectorType::dot_double( edge0, edge1 );
    double d11 = VectorType::dot_double( edge1, edge1 );
    double di0 = VectorType::dot_double( diff, edge0 );
    double di1 = VectorType::dot_double( diff, edge1 );
    double denominator = d00 * d11 - d01 * d01;
    double s = ( d11 * di0 - d01 * di1 ) / denominator;
    double t = ( d00 * di1 - d01 * di0 ) / denominator;

    return frantic::graphics::vector3fd( 1.0f - s - t, s, t );
}

/**
 * Determine if point 'p' is in Triangle(p0,p1,p2)
 * It is assumed that you have already determined that 'p' lies in the same plane (if in 3d)
 * @param p the query point
 * @param p0 the first point of the triangle
 * @param p1 the second point of the triangle
 * @param p2 the third point of the triangle
 * @return true only if the point lies inside the triangle
 */
template <class VectorType>
inline bool point_in_triangle( const VectorType& p, const VectorType& p0, const VectorType& p1, const VectorType& p2 ) {
    frantic::graphics::vector3fd coords = get_barycentric_coordinates( p, p0, p1, p2 );
    return coords[0] >= 0.0f && coords[1] >= 0.0f && coords[2] >= 0.0f;
}

/**
 * Computes the area of a polygon in 3d space
 *
 * @tparam InputIterator A read-only forward iterator which dereferences to a vector3f
 * @param begin iterator to first element
 * @param end iterator to one-past-end element
 * @return the computed area
 */
template <class InputIterator>
inline double polygon3_area( InputIterator begin, InputIterator end ) {
    using namespace frantic::graphics;

    InputIterator it = begin;

    vector3fd planePoints[3];
    size_t planeIndex = 0;

    vector3fd sumCross( 0.0f, 0.0f, 0.0f );

    while( it != end ) {
        vector3fd current = *it;
        if( planeIndex < 3 ) {
            planePoints[planeIndex] = current;
            ++planeIndex;
        }
        ++it;
        vector3fd next = ( it == end ) ? *begin : *it;
        sumCross += vector3fd::cross( current, next );
    }

    const vector3fd planeNormal = plane3fd::from_triangle( planePoints[0], planePoints[1], planePoints[2] ).normal();
    const double doubleArea = vector3fd::dot_double( sumCross, planeNormal );
    return std::abs( doubleArea ) / 2.0f;
}

/**
 * Computes the area of a triangle in 3d space
 *
 * @param p0 first point in the triangle
 * @param p1 second point in the triangle
 * @param p2 third point in the triangle
 * @return the computed area
 */
inline double triangle3_area( const frantic::graphics::vector3fd& p0, const frantic::graphics::vector3fd& p1,
                              const frantic::graphics::vector3fd& p2 ) {
    const frantic::graphics::vector3fd p[3] = { p0, p1, p2 };
    return polygon3_area( p, p + 3 );
}

/**
 * Collects the sequence of disjointed edges (the output of frantic::geometry::quickhull)
 * into an ordered polygon.
 *
 * @param edges An unordered list of edges, assumed to form a polygon if joined end to end
 * @param polygon The resulting ordered sequence from edges
 */
void collect_hull_polygon( const std::vector<frantic::graphics2d::vector2>& edges, std::vector<size_t>& polygon );

namespace polygon_winding {
enum polygon_winding {
    invalid = 0,
    clockwise,
    counter_clockwise,

    num_types,
};
}

inline bool counter_clockwise_convexity_functor2( const frantic::graphics2d::vector2f& vPrev,
                                                  const frantic::graphics2d::vector2f& vCurr,
                                                  const frantic::graphics2d::vector2f& vNext ) {
    return frantic::graphics2d::vector2f::cross( vNext - vCurr, vPrev - vCurr ) > 0.0f;
}

/**
 * Determines the 'convexity' of the given vertex sequence. It is assumed the vertices are presented
 * in the order they appear in the polygon (i.e. this will compute the 'convexity' of vCurr. This
 * assumes the winding of the polygon is clockwise when viewed from (0,0,1)
 *
 * @param vPrev The vertex before the vertex whose convexity is being tested
 * @param vCurr The current vertex in the sequence
 * @param vNext The vertex after the vertex whose convexity is being tested
 * @return true if vCurr is convex, false otherwise
 */
inline bool clockwise_convexity_functor2( const frantic::graphics2d::vector2f& vPrev,
                                          const frantic::graphics2d::vector2f& vCurr,
                                          const frantic::graphics2d::vector2f& vNext ) {
    return frantic::graphics2d::vector2f::cross( vNext - vCurr, vPrev - vCurr ) < 0.0f;
}

/**
 * Evaluates the 'winding' of the given polygon by finding an extreme vertex, and determining its cross product's sign\
 *
 * @param polygon list of vertices in the polygon
 * @param N number of vertices in the polygon
 * @return a value indicated the computed winding of the polygon
 */
inline polygon_winding::polygon_winding determine_polygon_winding2( const frantic::graphics2d::vector2f* polygon,
                                                                    const size_t N, const float epsilon = 1.0e-6f ) {

    // Find a max/min vertex in the x direction, and then compute the cross product
    frantic::graphics2d::vector2f minVert = polygon[0];
    size_t minIndex = 0;

    for( size_t i = 1; i < N; ++i ) {
        if( polygon[i].x < minVert.x || ( polygon[i].x == minVert.x && polygon[i].y < minVert.y ) ) {
            minIndex = i;
            minVert = polygon[i];
        }
    }

    const frantic::graphics2d::vector2f a = polygon[( minIndex + 1 ) % N] - polygon[minIndex];
    const frantic::graphics2d::vector2f b = polygon[( minIndex + N - 1 ) % N] - polygon[minIndex];

    const float cross = frantic::graphics2d::vector2f::cross( a, b );

    if( cross > std::abs( epsilon ) ) {
        return polygon_winding::counter_clockwise;
    } else if( cross < std::abs( -epsilon ) ) {
        return polygon_winding::clockwise;
    }

    return polygon_winding::invalid;
}

template <class InputIterator>
typename std::iterator_traits<InputIterator>::value_type get_polygon_normal3( InputIterator polyBegin,
                                                                              InputIterator polyEnd ) {
    using namespace frantic::graphics;

    typedef typename std::iterator_traits<InputIterator>::value_type vector_type;

    vector_type normal( 0.0, 0.0, 0.0 );

    vector_type first = *polyBegin;
    InputIterator it = polyBegin;

    while( it != polyEnd ) {
        vector_type current = *it;
        ++it;
        vector_type next = ( it == polyEnd ) ? first : *it;

        normal[0] += ( current[1] - next[1] ) * ( current[2] + next[2] );
        normal[1] += ( current[2] - next[2] ) * ( current[0] + next[0] );
        normal[2] += ( current[0] - next[0] ) * ( current[1] + next[1] );
    }

    return normal.to_normalized();
}

/**
 * Returns the normal of the given polygon, using Newell's method.
 * The orientation returned assumes the polygon winding order is CCW
 *
 * @param polygon the list of vertices in the polygon
 * @param N the number of vertices in the polygon
 * @return the surface normal of the polygon, in the orientation that makes it CCW
 */
template <class VectorType>
inline VectorType get_polygon_normal3( const VectorType* polygon, const size_t N ) {
    return get_polygon_normal3( polygon, polygon + N );
}

namespace detail {

// TODO: change the input to be a RandomAccessIterator>
template <class VectorType>
inline bool is_ear( const VectorType* polygon, size_t i0, size_t i1, size_t i2, const std::list<size_t>& indices ) {
    for( std::list<size_t>::const_iterator it = indices.begin(); it != indices.end(); ++it ) {
        if( *it != i0 && *it != i1 && *it != i2 &&
            frantic::geometry::point_in_triangle( polygon[*it], polygon[i0], polygon[i1], polygon[i2] ) ) {
            return false;
        }
    }

    return true;
}

// TODO: make this work with a RandomAccessIterator range rather than an array
// @return true if triangulation is successful.  If unsuccessful, returns false and sets outError to the error message
template <class VectorType, class ConvexityFunctor, class OutputIterator>
inline bool impl_triangulate_polygon_ear_clip( const VectorType* polygon, std::list<size_t>& indices,
                                               ConvexityFunctor convexTest, OutputIterator outTriangles,
                                               const char*& outError ) {
    using frantic::graphics::vector3;

    size_t outputCount = 0;
    size_t initialSize = indices.size();

    std::list<size_t>::iterator itPrev = indices.begin();
    std::list<size_t>::iterator it = itPrev;
    ++it;
    std::list<size_t>::iterator itNext = it;
    ++itNext;

    std::list<size_t>::iterator cycleEnd = itPrev;

    while( outputCount < initialSize - 3 ) {
        if( convexTest( polygon[*itPrev], polygon[*it], polygon[*itNext] ) &&
            is_ear( polygon, *itPrev, *it, *itNext, indices ) ) {
            *outTriangles = vector3( boost::int32_t( *itPrev ), boost::int32_t( *it ), boost::int32_t( *itNext ) );
            ++outTriangles;
            indices.erase( it );
            ++outputCount;
            cycleEnd = itPrev;
        } else {
            if( it == cycleEnd ) {
                outError = "triangulate_polygon_ear_clip -- infinite loop detected";
                return false;
            }
            itPrev = it;
        }

        it = itNext;
        ++itNext;

        if( itNext == indices.end() ) {
            itNext = indices.begin();
        }
    }

    if( indices.size() > 3 ) {
        outError = "triangulate_polygon_ear_clip -- failed to reduce the polygon to a single triangle";
        return false;
    }

    *outTriangles = vector3( boost::int32_t( *itPrev ), boost::int32_t( *it ), boost::int32_t( *itNext ) );
    return true;
}

template <class VectorType>
struct convexity_functor3 {
    VectorType m_normal;

    convexity_functor3( const VectorType& normal )
        : m_normal( normal ) {}

    bool operator()( const VectorType& vPrev, const VectorType& vCurr, const VectorType& vNext ) const {
        const VectorType crossProduct = VectorType::cross( vNext - vCurr, vPrev - vCurr );
        return VectorType::dot( m_normal, crossProduct ) > 0.0f;
    }
};

} // namespace detail

template <class VectorType, class ConvexityFunctor, class OutputIterator>
inline bool triangulate_polygon_ear_clip( const VectorType* polygon, const size_t N, ConvexityFunctor convexTest,
                                          OutputIterator outTriangles, const char*& outError ) {
    boost::integer_range<size_t>::type range = boost::irange( size_t( 0 ), N );
    std::list<size_t> indices( boost::begin( range ), boost::end( range ) );
    return detail::impl_triangulate_polygon_ear_clip( polygon, indices, convexTest, outTriangles, outError );
}

/**
 * Attempts to triangulate the given 2D polygon.
 *
 * This version returns true if triangulation was successful.  Otherwise, outError is set to an error message describing
 * what went wrong and returns false.
 *
 * @param polygon the list of vertices in the polygon
 * @param N the number of vertices in the polygon
 * @param outTriangles the list or iterator to insert resulting triangulation in
 * @param outError error message to set if triangulation is unsuccessful
 * @return true if triangulation is successful, false otherwise.
 */
template <class OutputIterator>
inline bool triangulate_polygon_nothrow( const frantic::graphics2d::vector2f* polygon, const size_t N,
                                         OutputIterator outTriangles, const char*& outError ) {
    polygon_winding::polygon_winding winding = determine_polygon_winding2( polygon, N );

    if( winding == polygon_winding::counter_clockwise ) {
        return triangulate_polygon_ear_clip( polygon, N, counter_clockwise_convexity_functor2, outTriangles, outError );
    } else {
        return triangulate_polygon_ear_clip( polygon, N, clockwise_convexity_functor2, outTriangles, outError );
    }
}

/**
 * Attempts to triangulate the given 3D polygon.
 *
 * This version returns true if triangulation was successful.  Otherwise, outError is set to an error message describing
 * what went wrong and returns false.
 *
 * @param polygon the list of vertices in the polygon
 * @param N the number of vertices in the polygon
 * @param outTriangles the list or iterator to insert resulting triangulation in
 * @param outError error message to set if triangulation is unsuccessful
 * @return true if triangulation is successful, false otherwise.
 */
template <class VectorType, class OutputIterator>
inline bool triangulate_polygon_nothrow( const VectorType* polygon, const size_t N, OutputIterator outTriangles,
                                         const char*& outError ) {
    VectorType normal = get_polygon_normal3( polygon, N );
    if( normal.get_magnitude_squared() < 0.0001f ) {
        outError = "triangulate_polygon -- Could not determine polygon normal";
        return false;
    }

    return triangulate_polygon_ear_clip( polygon, N, detail::convexity_functor3<VectorType>( normal ), outTriangles,
                                         outError );
}

/**
 * Attempts to triangulate the given 2D or 3D polygon.
 *
 * This version throws a std::logic_error exception if triangulation fails
 *
 * @param polygon the list of vertices in the polygon
 * @param N the number of vertices in the polygon
 * @param outTriangles the list or iterator to insert resulting triangulation in
 */
template <class VectorType, class OutputIterator>
inline void triangulate_polygon( const VectorType* polygon, const size_t N, OutputIterator outTriangles ) {
    const char* error = NULL;
    if( !triangulate_polygon_nothrow( polygon, N, outTriangles, error ) ) {
        throw std::logic_error( error );
    }
}

/**
 * Triangulate the given 2D or 3D polygon.
 *
 * This version always forces a triangulation even if the polygon itself is degenerate.  outTriangles is guaranteed to
 * always have N - 2 entries describing the triangles making up the polygon upon returning from this function call.
 *
 * If the polygon was degenerate in some way, it returns false and sets outError to the error message describing what
 * went wrong. Otherwise it returns true.
 *
 * @param polygon the list of vertices in the polygon
 * @param N the number of vertices in the polygon
 * @param outTriangles the list to insert resulting triangulation in
 * @param outError error message to set if triangulation is unsuccessful
 * @return true if triangulation is successful and the polygon was not degenerate in some way, false otherwise.
 */
template <class VectorType, class OutputVector>
inline bool triangulate_polygon_robust( const VectorType* polygon, const size_t N, OutputVector& outTriangles,
                                        const char*& outError ) {

    // Check the vertices to see if they all overlap.  This is to prevent calling triangulate_polygon and
    // having it fail with an "Infinite Loop Detected" error.  This doesn't catch all the cases though.
    bool canTriangulate = false;
    {
        const VectorType& checkVertex = polygon[0];
        for( std::size_t i = 1; i < N; ++i ) {
            if( checkVertex != polygon[i] ) {
                canTriangulate = true;
                break;
            }
        }
    }

    if( canTriangulate ) {
        canTriangulate = triangulate_polygon_nothrow( polygon, N, std::back_inserter( outTriangles ), outError );
    } else {
        outError = "triangulate_polygon -- polygon is degenerate";
    }

    if( !canTriangulate ) {
        outTriangles.clear();
        for( int i = 0; i < int( N - 2 ); ++i ) {
            outTriangles.push_back( frantic::graphics::vector3( 0, i + 1, i + 2 ) );
        }
    }

    return canTriangulate;
}

} // namespace geometry
} // namespace frantic
