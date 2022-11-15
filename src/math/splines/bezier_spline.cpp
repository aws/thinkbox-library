// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/math/polynomial_roots.hpp>
#include <frantic/math/splines/bezier_spline.hpp>
#include <frantic/math/utils.hpp>

using namespace std;
using namespace frantic::graphics2d;

namespace frantic {
namespace math {

float bezier_curve_x_to_y( const frantic::graphics2d::vector2f& position0,
                           const frantic::graphics2d::vector2f& control0, const frantic::graphics2d::vector2f& control1,
                           const frantic::graphics2d::vector2f& position1, float x ) {
    double a = position1.x - 3 * control1.x + 3 * control0.x - position0.x,
           b = 3 * control1.x - 6 * control0.x + 3 * position0.x, c = 3 * control0.x - 3 * position0.x,
           d = position0.x - x;
    double root1 = 0, root2 = 0, root3 = 0;
    int rootCount = 0;
    if( fabs( a ) < 0.00000001 * ( fabs( b ) + fabs( c ) + fabs( d ) ) ) {
        if( fabs( b ) < 0.00000001 * ( fabs( c ) + fabs( d ) ) ) {
            rootCount = 1;
            if( c != 0 )
                root1 = -d / c;
            else
                throw std::runtime_error(
                    "frantic::math::bezier_curve_x_to_y: Determined to be linear case, invalid parameter: d = 0" );
        } else {
            rootCount = frantic::math::polynomial_roots::get_quadratic_roots( b, c, d, root1, root2 );
        }
    } else {
        rootCount = frantic::math::polynomial_roots::get_cubic_roots( b / a, c / a, d / a, root1, root2, root3 );
    }

    if( rootCount == 0 ) {
        throw std::runtime_error( "bezier_curve_x_to_y: Bezier curve provided, " + position0.str() + " " +
                                  control0.str() + " " + control1.str() + " " + position1.str() +
                                  ", doesn't intersect x coordinate " + boost::lexical_cast<std::string>( x ) );
    } else {
        float t = (float)root1;
        // Find a root in the [0,1] interval
        if( root1 < -0.000001 || root1 > 1.000001 ) {
            if( rootCount > 1 && root2 >= -0.000001 && root2 <= 1.000001 ) {
                t = (float)root2;
            } else {
                if( rootCount > 2 && root3 >= -0.000001 && root3 <= 1.000001 ) {
                    t = (float)root3;
                } else {
                    std::string roots;
                    if( rootCount >= 1 )
                        roots += boost::lexical_cast<std::string>( root1 );
                    if( rootCount >= 2 )
                        roots += " " + boost::lexical_cast<std::string>( root2 );
                    if( rootCount >= 3 )
                        roots += " " + boost::lexical_cast<std::string>( root3 );
                    throw std::runtime_error(
                        "bezier_curve_x_to_y: Bezier curve provided, " + position0.str() + " " + control0.str() + " " +
                        control1.str() + " " + position1.str() + ", doesn't intersect x coordinate " +
                        boost::lexical_cast<std::string>( x ) + ".  Roots found: " + roots +
                        ", Cubic equation coefficients: " + boost::lexical_cast<std::string>( a ) + " " +
                        boost::lexical_cast<std::string>( b ) + " " + boost::lexical_cast<std::string>( c ) + " " +
                        boost::lexical_cast<std::string>( d ) );
                }
            }
        }

        t = math::clamp( t, 0.f, 1.f );

        return ( 1 - t ) * ( 1 - t ) * ( 1 - t ) * position0.y + ( 3 * t * ( 1 - t ) * ( 1 - t ) ) * control0.y +
               ( 3 * t * t * ( 1 - t ) ) * control1.y + ( t * t * t ) * position1.y;
    }
}

// Uses a curvature-based adaptive approach to convert the bezier curve into a linear polyline.
void convert_bezier_curve_to_samples( const frantic::graphics2d::vector2f& position0,
                                      const frantic::graphics2d::vector2f& control0,
                                      const frantic::graphics2d::vector2f& control1,
                                      const frantic::graphics2d::vector2f& position1, bool includeFirstPoint,
                                      float dPhi, std::vector<std::pair<float, float>>& outSamples ) {
    if( includeFirstPoint )
        outSamples.push_back( std::make_pair( position0.x, position0.y ) );

    float t = 0;
    t += get_bezier_dt_from_dPhi( position0, control0, control1, position1, dPhi, t );
    while( t < 1 ) {
        vector2f position = get_bezier_curve_position( position0, control0, control1, position1, t );
        outSamples.push_back( std::make_pair( position.x, position.y ) );

        t += get_bezier_dt_from_dPhi( position0, control0, control1, position1, dPhi, t );
    }

    outSamples.push_back( std::make_pair( position1.x, position1.y ) );
}

// Finds the nearest point on the bezier curve
float find_nearest_point_on_bezier_curve( const frantic::graphics2d::vector2f& position0,
                                          const frantic::graphics2d::vector2f& control0,
                                          const frantic::graphics2d::vector2f& control1,
                                          const frantic::graphics2d::vector2f& position1,
                                          const frantic::graphics2d::vector2f& p ) {
    // First off, one of the endpoints could be the closest point
    float tBest = 0, distanceBest = vector2f::distance( p, position0 );
    float tCandidate = 1, distanceCandidate = vector2f::distance( p, position1 );
    if( distanceCandidate < distanceBest ) {
        tBest = tCandidate;
        distanceBest = distanceCandidate;
    }

    // TODO: get the distance to the convex hull of the control points, and if that == the distance to position0 or
    // position1, we've already found the nearest point.

    // Next, points where the tangent of the curve is perpendicular to the vector from p to the position on the curve
    // could also be the closest point That is, solutions to the equation dot(curve'(t), curve(t) - p) = 0, which is a
    // fifth degree polynomial.  We attempt to find these points by breaking down the curve into twenty segments, and
    // solving for a numeric root wherever the function switches signs. Unproven hypotheses: If there's a solution where
    // the root doesn't switch signs, it can't be a minimum, because along one of the directions along t the curve will
    // be approach.
    float tAtLastSegment = 0;
    float valueAtLastSegment =
        vector2f::dot( get_bezier_curve_derivative( position0, control0, control1, position0, tAtLastSegment ),
                       get_bezier_curve_position( position0, control0, control1, position0, tAtLastSegment ) - p );
    int testSegments = 20;
    for( int segment = 0; segment < testSegments; ++segment ) {
        float tAtCurrentSegment = ( segment + 1.f ) / testSegments;
        float valueAtCurrentSegment = vector2f::dot(
            get_bezier_curve_derivative( position0, control0, control1, position0, tAtCurrentSegment ),
            get_bezier_curve_position( position0, control0, control1, position0, tAtCurrentSegment ) - p );
        // If there's a 0 crossing, do a couple of iterations of a binary search to approximate the root
        if( ( valueAtLastSegment >= 0 && valueAtCurrentSegment <= 0 ) ||
            ( valueAtLastSegment <= 0 && valueAtCurrentSegment >= 0 ) ) {
            // Initilize so that value0 < value1
            float t0, value0, t1, value1;
            if( valueAtLastSegment < valueAtCurrentSegment ) {
                t0 = tAtLastSegment;
                value0 = valueAtLastSegment;
                t1 = tAtCurrentSegment;
                value1 = valueAtCurrentSegment;
            } else {
                t1 = tAtLastSegment;
                value1 = valueAtLastSegment;
                t0 = tAtCurrentSegment;
                value0 = valueAtCurrentSegment;
            }
            // Do the binary search with a maximum number of iterations and a threshold deviation
            int maximumIterations = 10;
            float acceptableDeviation = 0.05f;
            while( ( std::min )( fabsf( value0 ), fabsf( value1 ) ) > acceptableDeviation && maximumIterations-- > 0 ) {
                float tMiddle = t0 + ( t1 - t0 ) * ( 0 - value0 ) / ( value1 - value0 );
                float valueMiddle =
                    vector2f::dot( get_bezier_curve_derivative( position0, control0, control1, position0, tMiddle ),
                                   get_bezier_curve_position( position0, control0, control1, position0, tMiddle ) - p );
                if( tMiddle < 0 ) {
                    t0 = tMiddle;
                    value0 = valueMiddle;
                } else {
                    t1 = tMiddle;
                    value1 = valueMiddle;
                }
            }
            // Check whether the candidates are better than the best one so far
            tCandidate = t0;
            distanceCandidate = vector2f::distance(
                p, get_bezier_curve_position( position0, control0, control1, position0, tCandidate ) );
            if( distanceCandidate < distanceBest ) {
                tBest = tCandidate;
                distanceBest = distanceCandidate;
            }
            tCandidate = t1;
            distanceCandidate = vector2f::distance(
                p, get_bezier_curve_position( position0, control0, control1, position0, tCandidate ) );
            if( distanceCandidate < distanceBest ) {
                tBest = tCandidate;
                distanceBest = distanceCandidate;
            }
        }
    }
    // Return the best minimum we found
    return tBest;
}

// Given the curve endpoints and two sample points on the curve, solves for the control points
void solve_for_bezier_curve_control_points( const frantic::graphics2d::vector2f& position0,
                                            const frantic::graphics2d::vector2f& position1, float t0,
                                            const frantic::graphics2d::vector2f& t0Position, float t1,
                                            const frantic::graphics2d::vector2f& t1Position,
                                            frantic::graphics2d::vector2f& outControl0,
                                            frantic::graphics2d::vector2f& outControl1 ) {
    /*
    cout << "solve_for_bezier_curve_control_points:" << endl;
    cout << "pos0: " << position0 << endl;
    cout << "t0Pos: " << t0Position << ", t0: " << t0 << endl;
    cout << "t1Pos: " << t1Position << ", t1: " << t1 << endl;
    cout << "pos1: " << position1 << endl;
    */

    // Basic formula is
    //      curve(t) = (t-1)*(t-1)*(t-1)*position0 + 3*t*(t-1)*(t-1)*control0 + 3*t*t*(t-1)*control1 + t*t*t*position1
    // We are given
    //      curve(t0) = t0Position
    //      curve(t1) = t1Position
    // So we get two pairs of linear equations, one for the x coordinates and one for the y coordinats, both looking
    // like this:
    //      3*t0*(1-t0)*(1-t0)*outControl0 + 3*t0*t0*(1-t0)*outControl1 = t0Position - (1-t0)*(1-t0)*(1-t0)*position0 -
    //      t0*t0*t0*position1 3*t1*(1-t1)*(1-t1)*outControl0 + 3*t1*t1*(1-t1)*outControl1 = t1Position -
    //      (1-t1)*(1-t1)*(1-t1)*position0 - t1*t1*t1*position1
    // Represented in matrix form
    //     [ 3*t0*(1-t0)*(1-t0)  3*t0*t0*(1-t0) ][ outControl0 ]   [ t0Position - (1-t0)*(1-t0)*(1-t0)*position0 -
    //     t0*t0*t0*position1 ] [ 3*t1*(1-t1)*(1-t1)  3*t1*t1*(1-t1) ][ outControl1 ] = [ t1Position -
    //     (1-t1)*(1-t1)*(1-t1)*position0 - t1*t1*t1*position1 ]
    // Or
    //     [ a b ] [ outControl0 ]   [ e ]
    //     [ c d ] [ outControl1 ] = [ f ]
    // Inverting the 2x2 matrix:
    //     [ outControl0 ]   [ ai bi ] [ e ]
    //     [ outControl1 ] = [ ci di ] [ f ]

    // Construct the 2x2 matrix
    double a = 3.0 * t0 * ( 1 - t0 ) * ( 1 - t0 ), b = 3.0 * t0 * t0 * ( 1 - t0 ),
           c = 3.0 * t1 * ( 1 - t1 ) * ( 1 - t1 ), d = 3.0 * t1 * t1 * ( 1 - t1 );
    // Invert the 2x2 matrix
    double determinant = a * d - b * c;
    double ai = d / determinant, bi = -b / determinant, ci = -c / determinant, di = a / determinant;

    /*
    cout << "abcd: " << a << " " << b << " " << c << " " << d << endl;
    cout << "abcd inverse: " << ai << " " << bi << " " << ci << " " << di << endl;
    */

    // Find the e, f coordinates
    vector2f e = t0Position - ( 1 - t0 ) * ( 1 - t0 ) * ( 1 - t0 ) * position0 - t0 * t0 * t0 * position1;
    vector2f f = t1Position - ( 1 - t1 ) * ( 1 - t1 ) * ( 1 - t1 ) * position0 - t1 * t1 * t1 * position1;

    // Multiply the inverse matrix by the e, f coordinates to get the control points
    outControl0.x = float( ai * e.x + bi * f.x );
    outControl0.y = float( ai * e.y + bi * f.y );
    outControl1.x = float( ci * e.x + di * f.x );
    outControl1.y = float( ci * e.y + di * f.y );

    /*
    cout << "outControl0: " << outControl0 << endl;
    cout << "outControl1: " << outControl1 << endl;
    cout << "curve at point t0: " << get_bezier_curve_position( position0, outControl0, outControl1, position1, t0 ) <<
    endl; cout << "curve at point t1: " << get_bezier_curve_position( position0, outControl0, outControl1, position1, t1
    )
    << endl;
    */
}

// Assuming polyLine.front() and polyLine.back() are the start and end positions,
// Computes the least squares control points corresponding to the t values provided
// for the poly line vertices.
bool best_fit_bezier_center_control_points( const std::vector<frantic::graphics2d::vector2f>& polyLine,
                                            const std::vector<float>& polyLineTValues,
                                            frantic::graphics2d::vector2f& outControl0,
                                            frantic::graphics2d::vector2f& outControl1 ) {
    if( polyLine.size() < 2 || polyLine.size() != polyLineTValues.size() )
        return false;

    if( polyLine.size() == 2 ) {
        // Straight line is an exact fit
        outControl0 = 0.667f * polyLine.front() + 0.333f * polyLine.back();
        outControl1 = 0.333f * polyLine.front() + 0.667f * polyLine.back();
        return true;
    }

    // Come up with something for this!

    // fit it to a quadradic bezier, then just convert it to a cubic directly (essentially what we're currently doing
    // for 2 sample points to a linear bezier)
    if( polyLine.size() == 3 ) {
        float t = polyLineTValues[1];
        float tinv = 1.0f - t;
        float t2 = t * t;
        float tinv2 = tinv * tinv;

        // solve for the middle control point of the quadradic bezier that fits these 3 points
        vector2f B = ( polyLine[1] - ( polyLine[0] * tinv2 ) - ( polyLine[2] * t2 ) ) / ( 2.0f * tinv * t );

        // then set the cubic center control points to be 2/3 along the lines to the quadradic center control point
        //(I believe this is equivalent, please correct me if I'm wrong)
        outControl0 = 0.667f * B + 0.333f * polyLine.front();
        outControl1 = 0.667f * B + 0.333f * polyLine.back();

        return true;
    }

    // Let's say the four control points are A, B, C, and D.  We know A and D, and have
    // estimates for t for each point X. For each point X with parameter estimate t, this
    // gives the system:
    //
    // X = (1-t)^3 * A + 3*t*(1-t)^2 * B + 3*t^2*(1-t) * C + t^3 * D
    //
    // Or in matrix form, isolating B and C,
    //
    // ( 3*t*(1-t)^2  3*t^2*(1-t) )( B ) = ( X - (1-t)^3 * A - t^3 * D )
    //                             ( C )
    //
    // Combining all the points produces 2 (usually) overconstrained linear systems
    // looking like this (substituting x and y for #):
    //
    // ( 3*t1*(1-t1)^2  3*t1^2*(1-t1) )( B.# )   ( X1.# - (1-t1)^3 * A.# - t1^3 * D.# )
    // ( 3*t2*(1-t2)^2  3*t2^2*(1-t2) )( C.# ) = ( X2.# - (1-t2)^3 * A.# - t2^3 * D.# )
    // ( 3*t3*(1-t3)^2  3*t3^2*(1-t3) )          ( X3.# - (1-t3)^3 * A.# - t3^3 * D.# )
    // (             ...              )          (                ...                 )
    //
    // We solve this by finding the pseudo inverse, A^(+) = (A^(T)A)^(-1)A^(T).

    // To minimize memory allocations, we calculate two parts separately.
    // A^(T)A is a 2x2 matrix, and A^(T)R where R is the right hand side vector is a 2-vector.

    // Atranspose * A is symmetric, so we only need three values for it
    float Atranspose_A[3] = { 0, 0, 0 };
    vector2f Atranspose_b[2];
    for( unsigned i = 1; i < polyLine.size() - 1; ++i ) {
        float t = polyLineTValues[i];
        float temp = 3 * t * ( 1 - t );
        float x = temp * ( 1 - t ), y = temp * t;
        Atranspose_A[0] += x * x;
        Atranspose_A[1] += x * y;
        Atranspose_A[2] += y * y;

        vector2f r = polyLine[i] - ( 1 - t ) * ( 1 - t ) * ( 1 - t ) * polyLine.front() - t * t * t * polyLine.back();
        Atranspose_b[0] += x * r;
        Atranspose_b[1] += y * r;
    }

    float determinant = Atranspose_A[0] * Atranspose_A[2] - Atranspose_A[1] * Atranspose_A[1];
    if( determinant != 0 ) {
        // Compute the inverse of Atranspose * A
        float inverse_Atranspose_A[3];
        inverse_Atranspose_A[0] = Atranspose_A[2] / determinant;
        inverse_Atranspose_A[1] = -Atranspose_A[1] / determinant;
        inverse_Atranspose_A[2] = Atranspose_A[0] / determinant;

        // And complete the multiplication by the pseudo-inverse to get the desired vectors
        outControl0 = inverse_Atranspose_A[0] * Atranspose_b[0] + inverse_Atranspose_A[1] * Atranspose_b[1];
        outControl1 = inverse_Atranspose_A[1] * Atranspose_b[0] + inverse_Atranspose_A[2] * Atranspose_b[1];

        return true;
    } else {
        // A typical case of a singular matrix will be a straight line, so set it to a straight
        // line
        outControl0 = 0.667f * polyLine.front() + 0.333f * polyLine.back();
        outControl1 = 0.333f * polyLine.front() + 0.667f * polyLine.back();
        return true;
    }
}

// Does a least squares style best-fit approximation using a Bezier cubic spline
// Returns true if it is returning an answer, false otherwise.
bool best_fit_bezier_cubic_spline( const std::vector<frantic::graphics2d::vector2f>& polyLine,
                                   frantic::graphics2d::vector2f& outPosition0,
                                   frantic::graphics2d::vector2f& outPosition1,
                                   frantic::graphics2d::vector2f& outControl0,
                                   frantic::graphics2d::vector2f& outControl1 ) {
    if( polyLine.size() == 2 ) {
        // Straight line is an exact fit
        outPosition0 = polyLine.front();
        outPosition1 = polyLine.back();
        outControl0 = 0.667f * outPosition0 + 0.333f * outPosition1;
        outControl1 = 0.333f * outPosition0 + 0.667f * outPosition1;
        return true;
    }

    // TODO: Come up with a strategy for the best fit with 3 points.
    // there now is a strategy for best fit with 3 points

    // The initial and end positions are the first and last points
    outPosition0 = polyLine.front();
    outPosition1 = polyLine.back();

    ////////////////
    // Compute an initial t value for each point in the polyline
    ////////////////
    float polyLineLength = 0;
    std::vector<float> polyLineTValues;
    polyLineTValues.resize( polyLine.size() );
    // Set the initial t estimates to be based on the distance along the polyline
    polyLineTValues[0] = 0;
    for( unsigned i = 1; i < polyLineTValues.size(); ++i ) {
        polyLineLength += vector2f::distance( polyLine[i - 1], polyLine[i] );
        polyLineTValues[i] = polyLineLength;
    }
    // Normalize all the T values based on the total length
    for( unsigned i = 1; i < polyLineTValues.size(); ++i ) {
        polyLineTValues[i] /= polyLineLength;
    }

    if( !best_fit_bezier_center_control_points( polyLine, polyLineTValues, outControl0, outControl1 ) )
        return false;

    // TODO: Produce a better estimate for the T values, and then do the best fit again.
    //       Iterate until reasonable convergence

    return true;
}

float bezier_curve_x_to_y( const std::vector<bezier_curve_point<frantic::graphics2d::vector2f>>& points, float x ) {
    size_t index = 0;
    // find the two points x is inbetween
    while( index < points.size() && x > points[index].position.x )
        index++;

    if( index == 0 ) {
        if( points.size() > 0 )
            return points[0].position.y;
        else
            return 0;
    } else if( index >= points.size() )
        return points[points.size() - 1].position.y;
    else {
        // check if it is right on the point
        // the bezier_curve_x_to_y does not seem to handle this
        if( points[index].position.x == x )
            return points[index].position.y;

        vector2f control0 = points[index - 1].position + points[index - 1].outTan;
        vector2f control1 = points[index].position + points[index].inTan;

        return bezier_curve_x_to_y( points[index - 1].position, control0, control1, points[index].position, x );
    }
}

} // namespace math
} // namespace frantic
