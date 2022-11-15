// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

#include <frantic/geometry/raytracing.hpp>
#include <frantic/graphics/graphics_utils.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace geometry {

namespace detail {
/**
 * This version assumes that pt is in the triangle plane but not inside the triangle, and baryC represents it.
 * Useful for meshes which have pre-computed the face planes and barycentric helpers.
 *
 * "pt" and "baryC" are updated to reflect the location in the triangle closest to the input "pt" and "baryC".
 */
inline void nearest_point_on_triangle( vector3f& pt, vector3f& baryC, const vector3f& A, const vector3f& B,
                                       const vector3f& C ) {
    if( baryC.x <= 0 ) {
        if( baryC.y <= 0 ) {
            // Check against BC and AC
            vector3f AC = ( A - C ), BC = ( B - C ), PTC = ( pt - C );
            float a = vector3f::dot( AC, PTC ) / AC.get_magnitude_squared();
            float b = vector3f::dot( BC, PTC ) / BC.get_magnitude_squared();
            if( a > 0 ) {
                if( a >= 1 ) {
                    pt = A;
                    baryC = vector3f( 1, 0, 0 );
                } else {
                    pt = a * A + ( 1 - a ) * C;
                    baryC = vector3f( a, 0, 1 - a );
                }
            } else if( b > 0 ) {
                if( b >= 1 ) {
                    pt = B;
                    baryC = vector3f( 0, 1, 0 );
                } else {
                    pt = b * B + ( 1 - b ) * C;
                    baryC = vector3f( 0, b, 1 - b );
                }
            } else {
                pt = C;
                baryC = vector3f( 0, 0, 1 );
            }
        } else if( baryC.z <= 0 ) {
            // Check against AB and CB
            vector3f AB = ( A - B ), CB = ( C - B ), PB = ( pt - B );
            float a = vector3f::dot( AB, PB ) / AB.get_magnitude_squared();
            float c = vector3f::dot( CB, PB ) / CB.get_magnitude_squared();
            if( a > 0 ) {
                if( a >= 1 ) {
                    pt = A;
                    baryC = vector3f( 1, 0, 0 );
                } else {
                    pt = a * A + ( 1 - a ) * B;
                    baryC = vector3f( a, 1 - a, 0 );
                }
            } else if( c > 0 ) {
                if( c >= 1 ) {
                    pt = C;
                    baryC = vector3f( 0, 0, 1 );
                } else {
                    pt = c * C + ( 1 - c ) * B;
                    baryC = vector3f( 0, 1 - c, c );
                }
            } else {
                pt = B;
                baryC = vector3f( 0, 1, 0 );
            }
        } else {
            // Check against BC or alternatively CB
            vector3f BC = ( B - C );
            float b = vector3f::dot( BC, ( pt - C ) ) / BC.get_magnitude_squared();
            if( b <= 0 ) {
                pt = C;
                baryC = vector3f( 0, 0, 1 );
            } else if( b >= 1 ) {
                pt = B;
                baryC = vector3f( 0, 1, 0 );
            } else {
                pt = b * B + ( 1 - b ) * C;
                baryC = vector3f( 0, b, 1 - b );
            }
        }
    } else if( baryC.y <= 0 ) {
        if( baryC.z <= 0 ) {
            // Check against CA and BA
            vector3f BA = ( B - A ), CA = ( C - A ), PA = ( pt - A );
            float b = vector3f::dot( BA, PA ) / BA.get_magnitude_squared();
            float c = vector3f::dot( CA, PA ) / CA.get_magnitude_squared();
            if( b > 0 ) {
                if( b >= 1 ) {
                    pt = B;
                    baryC = vector3f( 0, 1, 0 );
                } else {
                    pt = b * B + ( 1 - b ) * A;
                    baryC = vector3f( 1 - b, b, 0 );
                }
            } else if( c > 0 ) {
                if( c >= 1 ) {
                    pt = C;
                    baryC = vector3f( 0, 0, 1 );
                } else {
                    pt = c * C + ( 1 - c ) * A;
                    baryC = vector3f( 1 - c, 0, c );
                }
            } else {
                pt = A;
                baryC = vector3f( 1, 0, 0 );
            }
        } else {
            // Check against AC or alternatively CA
            vector3f AC = ( A - C );
            float a = vector3f::dot( AC, ( pt - C ) ) / AC.get_magnitude_squared();
            if( a <= 0 ) {
                pt = C;
                baryC = vector3f( 0, 0, 1 );
            } else if( a >= 1 ) {
                pt = A;
                baryC = vector3f( 1, 0, 0 );
            } else {
                pt = a * A + ( 1 - a ) * C;
                baryC = vector3f( a, 0, 1 - a );
            }
        }
    } else if( baryC.z <= 0 ) {
        // Check against AB or alternatively BA
        vector3f AB = ( A - B );
        float a = vector3f::dot( AB, ( pt - B ) ) / AB.get_magnitude_squared();
        if( a <= 0 ) {
            pt = B;
            baryC = vector3f( 0, 1, 0 );
        } else if( a >= 1 ) {
            pt = A;
            baryC = vector3f( 1, 0, 0 );
        } else {
            pt = a * A + ( 1 - a ) * B;
            baryC = vector3f( a, 1 - a, 0 );
        }
    }
}
} // End of namespace detail

/**
 * Calculates the point closest to p on the triangle tri(A,B,C)
 */
inline vector3f nearest_point_on_triangle( const vector3f& p, const vector3f& A, const vector3f& B,
                                           const vector3f& C ) {
    using frantic::graphics::vector3f;

    vector3f norm = vector3f::normalize( vector3f::cross( B - A, C - A ) );
    vector3f proj = p - vector3f::dot( ( p - A ), norm ) * norm;
    vector3f baryC = compute_barycentric_coordinates( proj, A, B, C );

    detail::nearest_point_on_triangle( proj, baryC, A, B, C );
    return proj;
}

/**
 * Calculates the point closest to p on the triangle tri(A,B,C) and returns the barycentric coordinate of that point.
 */
inline vector3f nearest_barycoord_on_triangle( const vector3f& p, const vector3f& A, const vector3f& B,
                                               const vector3f& C ) {
    using frantic::graphics::vector3f;

    vector3f norm = vector3f::normalize( vector3f::cross( B - A, C - A ) );
    vector3f proj = p - vector3f::dot( ( p - A ), norm ) * norm;
    vector3f baryC = compute_barycentric_coordinates( proj, A, B, C );

    detail::nearest_point_on_triangle( proj, baryC, A, B, C );
    return baryC;
}

inline void nearest_point_on_triangle( const vector3f& p, const vector3f& A, const vector3f& B, const vector3f& C,
                                       nearest_point_search_result& npsr ) {
    using frantic::graphics::vector3f;

    vector3f norm = vector3f::normalize( vector3f::cross( B - A, C - A ) );
    vector3f proj = p - vector3f::dot( ( p - A ), norm ) * norm;
    vector3f baryC = compute_barycentric_coordinates( proj, A, B, C );

    detail::nearest_point_on_triangle( proj, baryC, A, B, C );

    npsr.distance = vector3f::distance( p, proj );
    npsr.position = proj;
    npsr.geometricNormal = norm;
    npsr.faceIndex = -1;
    npsr.barycentricCoords = baryC;
}

inline bool compute_dPduv( const vector3f ( &tri )[3], const float ( &uvs )[3][2], vector3f& dPdu, vector3f& dPdv ) {
    // We assume a surface parameterization of a triangle with UV values at the corners (only UVs, no W, because a
    // surface has 2 degrees of freedom). We can define the triangle as all points P + dPdu(u) + dPdv(v). Therefore,
    // Pi-P0 = dPDu(ui-u0) + dPdv(vi-v0) which can be rearranged to solve for dP/du & dP/dv.
    vector3f e1 = tri[1] - tri[0];
    vector3f e2 = tri[2] - tri[0];

    float u1 = uvs[1][0] - uvs[0][0];
    float u2 = uvs[2][0] - uvs[0][0];
    float v1 = uvs[1][1] - uvs[0][1];
    float v2 = uvs[2][1] - uvs[0][1];

    float determinant = u1 * v2 - u2 * v1;
    if( std::abs( determinant ) > 1e-5f ) { // Cannot divide by 0 determinant.
        dPdu = ( v2 * e1 - v1 * e2 ) / determinant;
        dPdv = ( u1 * e2 - u2 * e1 ) / determinant;
        return true;
    } else {
        return false;
    }
}

template <class RandomEngine>
inline void random_barycentric_coordinate( float ( &outBaryCoords )[3], RandomEngine& eng ) {
    boost::variate_generator<RandomEngine&, boost::random::uniform_01<float>> rnd( eng,
                                                                                   boost::random::uniform_01<float>() );

    float ba = rnd(), bb = rnd();
    float bc = ba + bb;

    if( bc > 1.f ) {
        ba = 1.f - ba;
        bb = 1.f - bb;
        bc = bc - 1.f;
    } else {
        bc = 1.f - bc;
    }

    outBaryCoords[0] = ba;
    outBaryCoords[1] = bb;
    outBaryCoords[2] = bc;
}

template <class RandomEngine>
inline frantic::graphics::vector3f random_point_in_triangle( const vector3f& a, const vector3f& b, const vector3f& c,
                                                             RandomEngine& eng ) {
    float baryCoords[3];

    random_barycentric_coordinate( baryCoords, eng );

    return baryCoords[0] * a + baryCoords[1] * b + baryCoords[2] * c;
}

} // namespace geometry
} // namespace frantic
