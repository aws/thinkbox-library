// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/polygon_utils.hpp>

#include <frantic/geometry/quickhull2.hpp>

#include <frantic/math/utils.hpp>

#include <Eigen/Geometry>

#include <map>
#include <vector>

namespace {

typedef Eigen::Hyperplane<float, 2> Line2f;

frantic::graphics2d::vector2f from_eigen_t( const Eigen::Vector2f& v ) {
    return frantic::graphics2d::vector2f( v[0], v[1] );
}

Eigen::Vector2f to_eigen_t( const frantic::graphics2d::vector2f& v ) { return Eigen::Vector2f( v.x, v.y ); }

Eigen::Vector2f to_eigen_t( const frantic::graphics2d::vector2& v ) {
    return Eigen::Vector2f( float( v.x ), float( v.y ) );
}

} // namespace

namespace frantic {
namespace geometry {

namespace {

template <typename VectorType>
bool point_on_polygon_bound_raycast_impl( const std::vector<VectorType>& polygon, const VectorType& point,
                                          const typename VectorType::value_type tolerance ) {
    size_t j = polygon.size() - 1;

    for( size_t i = 0; i < polygon.size(); ++i ) {
        if( ( polygon[i].y <= point.y && polygon[j].y >= point.y ) ||
            ( polygon[j].y <= point.y && polygon[i].y >= point.y ) ) {
            Line2f line( Line2f::Through( to_eigen_t( polygon[i] ), to_eigen_t( polygon[j] ) ) );

            if( line.absDistance( to_eigen_t( point ) ) <= tolerance ) {
                return true;
            }
        }

        j = i;
    }

    return false;
}

template <>
bool point_on_polygon_bound_raycast_impl<frantic::graphics2d::vector2>(
    const std::vector<frantic::graphics2d::vector2>& polygon, const frantic::graphics2d::vector2& point,
    int tolerance ) {
    size_t j = polygon.size() - 1;

    for( size_t i = 0; i < polygon.size(); ++i ) {

        if( ( polygon[i].y <= point.y && polygon[j].y >= point.y ) ||
            ( polygon[j].y <= point.y && polygon[i].y >= point.y ) ) {
            const int rise = polygon[i].y - polygon[j].y;
            const int run = polygon[i].x - polygon[j].x;

            if( frantic::math::get_absolute( rise ) > tolerance ) {
                const int intersectionPoint = polygon[i].x + ( ( run * ( point.y - polygon[i].y ) ) / rise );
                if( intersectionPoint == point.x ) {
                    return true;
                }
            } else if( frantic::math::get_absolute( polygon[i].y - polygon[j].y ) <= tolerance &&
                       frantic::math::get_absolute( polygon[i].y - point.y ) <= tolerance &&
                       ( ( polygon[i].x <= point.x && polygon[j].x >= point.x ) ||
                         ( polygon[j].x <= point.x && polygon[i].x >= point.x ) ) ) {
                return true;
            }
        }

        j = i;
    }

    return false;
}

template <typename VectorType>
boundary_relation::boundary_relation point_in_polygon_raycast_impl( const std::vector<VectorType>& polygon,
                                                                    const VectorType& point,
                                                                    const typename VectorType::value_type tolerance ) {

    bool bound = point_on_polygon_bound_raycast_impl<VectorType>( polygon, point, tolerance );

    if( bound ) {
        return boundary_relation::boundary;
    }

    size_t j = polygon.size() - 1;

    bool parity = false;

    for( size_t i = 0; i < polygon.size(); ++i ) {
        if( ( polygon[i].y < point.y && polygon[j].y >= point.y ) ||
            ( polygon[j].y < point.y && polygon[i].y >= point.y ) ) {
            const typename VectorType::value_type rise = polygon[i].y - polygon[j].y;
            const typename VectorType::value_type run = polygon[i].x - polygon[j].x;

            if( frantic::math::get_absolute( rise ) > tolerance ) {
                const typename VectorType::value_type intersectionPoint =
                    polygon[i].x + ( ( run * ( point.y - polygon[i].y ) ) / rise );
                if( intersectionPoint < point.x ) {
                    parity = !parity;
                }
            }
        }

        j = i;
    }

    return parity ? boundary_relation::inside : boundary_relation::outside;
}

/**
 * Returns the relative 'quadrant' of q, assuming 'base' is the origin
 * 1 | 0
 * --+--
 * 2 | 3
 *
 * @param q the query point
 * @param base the point acting as the origin
 * @return the associated quadrant
 * @remark I'm not making an enum for this, just deal with it
 * @remark points that lie exactly on axis boundaries are (arbitrarily)
 *         assigned to the most positive quadrant that axis is adjacent to.
 *         So, for example, a point of (-3,0) would be assigned to Q1, and
 *         a point of (0, -4) would be assigned to Q3. (0,0) is assigned to Q0
 *         I'm fairly certain that this works out, as long as the resolution
 *         used is consistent for the entire algorithm
 */
template <typename VectorType>
int get_quadrant( const VectorType& q, const VectorType& base ) {
    if( q.x < base.x ) {
        if( q.y < base.y ) {
            return 2;
        } else {
            return 1;
        }
    } else {
        if( q.y < base.y ) {
            return 3;
        } else {
            return 0;
        }
    }
}

template <typename VectorType>
boundary_relation::boundary_relation point_in_polygon_winding_impl( const std::vector<VectorType>& polygon,
                                                                    const VectorType& point,
                                                                    const typename VectorType::value_type tolerance ) {
    VectorType lastPoint = polygon[polygon.size() - 1];
    int lastQuadrant = get_quadrant( polygon[polygon.size() - 1], point );

    int summation = 0;

    for( size_t i = 0; i < polygon.size(); ++i ) {
        VectorType currentPoint = polygon[i];
        int quadrant = get_quadrant( polygon[i], point );

        if( VectorType::distance_squared( point, currentPoint ) <= tolerance ) {
            return boundary_relation::boundary;
        }

        int transition = 0;

        // compute (start - end) mod 4
        const int transitionMod4 = ( ( lastQuadrant + 4 ) - quadrant ) % 4;

        // We single out the check for co-boundary and double-quadrant crossing
        // so we only ever compute the line equation once per iteration
        if( point.x >= std::min( lastPoint.x, currentPoint.x ) && point.y >= std::min( lastPoint.y, currentPoint.y ) &&
            point.x <= std::max( lastPoint.x, currentPoint.x ) && point.y <= std::max( lastPoint.y, currentPoint.y ) ) {
            Line2f line( Line2f::Through( to_eigen_t( lastPoint ), to_eigen_t( currentPoint ) ) );
            float dist = line.signedDistance( to_eigen_t( point ) );
            if( frantic::math::get_absolute( dist ) <= tolerance ) {
                return boundary_relation::boundary;
            } else if( transitionMod4 == 2 ) {
                if( dist > 0.0f ) {
                    transition = -2;
                } else {
                    transition = 2;
                }
            }
        }

        switch( transitionMod4 ) {
        // no change in quadrant
        case 0:
            transition = 0;
            break;
        // clockwise change in quadrant
        case 1:
            transition = 1;
            break;
        case 2:
            // this case is handled by the line check above
#if !defined( NDEBUG ) // quick sanity check
            if( frantic::math::get_absolute( transition ) != 2 ) {
                throw std::logic_error(
                    "point_in_polygon_winding_impl -- did not correctly assign transition value above for "
                    "double quadrant crossing" );
            }
#endif
            break;
        // counter-clockwise change in quadrant
        case 3:
            transition = -1;
            break;
        default:
            throw std::logic_error( "point_in_polygon_winding_impl -- invalid transition state" );
        }

        summation += transition;

        lastPoint = currentPoint;
        lastQuadrant = quadrant;
    }

    return summation == 0 ? boundary_relation::outside : boundary_relation::inside;
}

struct vertex_indexer {
    size_t index[2];
    size_t coord[2];
    bool filled[2];

    void add_index( size_t i, size_t c ) {
        if( filled[0] ) {
            if( filled[1] ) {
                throw std::logic_error( "collect_quickhull_polygon -- invalid input edge set" );
            } else {
                index[1] = i;
                coord[1] = c;
                filled[1] = true;
            }
        } else {
            index[0] = i;
            coord[0] = c;
            filled[0] = true;
        }
    }
};

} // namespace

void collect_hull_polygon( const std::vector<frantic::graphics2d::vector2>& edges, std::vector<size_t>& polygon ) {
    std::map<size_t, vertex_indexer> indices;

    // Fill out an index structure, indicating the location of each vertex in the edge set
    // though gaps will exist for non-convex-hull vertices
    for( size_t i = 0; i < edges.size(); ++i ) {
        indices[edges[i].x].add_index( i, 0 );
        indices[edges[i].y].add_index( i, 1 );
    }

    size_t j = 0;
    size_t currentVert = edges[j].x;
    size_t currentEntry = 0;

    polygon.clear();
    polygon.push_back( currentVert );

    while( polygon.size() < edges.size() ) {
        size_t endIndex = indices[currentVert].index[( currentEntry + 1 ) % 2];
        size_t endCoord = indices[currentVert].coord[( currentEntry + 1 ) % 2];

        size_t nextVert = edges[endIndex][( endCoord + 1 ) % 2];

        polygon.push_back( nextVert );

        if( indices[nextVert].index[0] == endIndex ) {
            currentEntry = 0;
        } else {
            currentEntry = 1;
        }

        currentVert = nextVert;
    }
}

namespace {

inline double distance_to_supporting_line( const frantic::graphics2d::vector2f& segStart,
                                           const frantic::graphics2d::vector2f& segEnd,
                                           const frantic::graphics2d::vector2f& point ) {
    Line2f line( Line2f::Through( to_eigen_t( segStart ), to_eigen_t( segEnd ) ) );
    return line.absDistance( to_eigen_t( point ) );
}

} // namespace

boundary_relation::boundary_relation point_in_polygon( const std::vector<frantic::graphics2d::vector2>& polygon,
                                                       const frantic::graphics2d::vector2& point ) {
    return point_in_polygon_winding_impl<frantic::graphics2d::vector2>( polygon, point, 0 );
}

boundary_relation::boundary_relation point_in_polygon( const std::vector<frantic::graphics2d::vector2f>& polygon,
                                                       const frantic::graphics2d::vector2f& point, float tolerance ) {
    return point_in_polygon_winding_impl<frantic::graphics2d::vector2f>( polygon, point, tolerance );
}

namespace {

void make_enclosing_rectangle( const Eigen::Vector2f& normal, const Eigen::Vector2f& parMin,
                               const Eigen::Vector2f& parMax, const Eigen::Vector2f& perpMin,
                               const Eigen::Vector2f& perpMax, Eigen::Vector2f* outPoints ) {

    Line2f lMin( normal, parMin );
    Line2f lMax( normal, parMax );
    Line2f pMin( Eigen::Vector2f( -normal[1], normal[0] ), perpMin );
    Line2f pMax( Eigen::Vector2f( -normal[1], normal[0] ), perpMax );

    outPoints[0] = lMin.intersection( pMin );
    outPoints[1] = lMin.intersection( pMax );
    outPoints[2] = lMax.intersection( pMax );
    outPoints[3] = lMax.intersection( pMin );
}

float rectangle_area_squared( const Eigen::Vector2f* rect ) {
    Eigen::Vector2f length = rect[1] - rect[0];
    Eigen::Vector2f height = rect[2] - rect[1];
    return length.dot( length ) * height.dot( height );
}

} // namespace

void minimum_enclosing_rectange( const std::vector<frantic::graphics2d::vector2f>& points,
                                 std::vector<frantic::graphics2d::vector2f>& outRect ) {
    using namespace frantic::graphics2d;

    std::vector<vector2> convexHull;
    frantic::geometry::quickhull( points, convexHull );

    std::vector<Eigen::Vector2f> ePolygon( convexHull.size() );
    for( size_t i = 0; i < convexHull.size(); ++i ) {
        ePolygon[i] = to_eigen_t( points[convexHull[i].x] );
    }

    // simple quadratic algorithm for finding minimum enclosing rectangle
    // A linear algorithm exists, but it is more complex, and I am assuming
    // nearly all polygons will be 3 or 4 verts only

    Eigen::Vector2f minRect[4];
    float minArea = std::numeric_limits<float>::max();

    for( size_t i = 0; i < ePolygon.size(); ++i ) {
        Line2f lines[2];
        lines[0] = Line2f::Through( ePolygon[i], ePolygon[( i + 1 ) % ePolygon.size()] );

        lines[1] = Line2f::Through( lines[0].normal() + ePolygon[i], ePolygon[i] );

        bool setMin[2] = { false, false };
        float minDist[2] = { std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
        size_t minIndex[2];
        bool setMax[2] = { false, false };
        float maxDist[2] = { -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max() };
        size_t maxIndex[2];

        for( size_t j = 0; j < ePolygon.size(); ++j ) {
            for( size_t k = 0; k < 2; ++k ) {
                float dist = lines[k].signedDistance( ePolygon[j] );

                if( !setMin[k] || dist < minDist[k] ) {
                    minDist[k] = dist;
                    minIndex[k] = j;
                    setMin[k] = true;
                }
                if( !setMax[k] || dist > maxDist[k] ) {
                    maxDist[k] = dist;
                    maxIndex[k] = j;
                    setMax[k] = true;
                }
            }
        }

        Eigen::Vector2f rect[4];
        make_enclosing_rectangle( lines[0].normal(), ePolygon[minIndex[0]], ePolygon[maxIndex[0]],
                                  ePolygon[minIndex[1]], ePolygon[maxIndex[1]], rect );
        float squaredArea = rectangle_area_squared( rect );

        if( squaredArea < minArea ) {
            minArea = squaredArea;
            for( size_t j = 0; j < 4; ++j ) {
                minRect[j] = rect[j];
            }
        }
    }

    // copy the result to our output structure
    outRect.resize( 4 );
    std::transform( minRect, minRect + 4, outRect.begin(), from_eigen_t );
}

} // namespace geometry
} // namespace frantic
