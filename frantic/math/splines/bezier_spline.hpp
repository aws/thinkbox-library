// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics2d/vector2f.hpp>
#include <frantic/math/iterative_root_solve.hpp>
#include <utility>
#include <vector>

namespace frantic {
namespace math {

/**
 * Evaluates a cubic, bezier spline at the specified parametric value.
 * @tparam VectorType The vector type for the spline. Typically vector2f or vector3f.
 * @param position0 The start of the spline
 * @param control0 The first control point of the spline. This point is generally not ON the spline.
 * @param control1 The second control point of the spline. This point is generally not ON the spline.
 * @param position1 The end of the spline
 * @param t The parametric value to evaluate the spline at. Must be in the interval [0, 1].
 * @return The spline position at the supplied parametric value.
 */
template <typename VectorType, typename TimeType>
inline VectorType get_bezier_curve_position( const VectorType& position0, const VectorType& control0,
                                             const VectorType& control1, const VectorType& position1,
                                             const TimeType& t ) {
    const TimeType t2 = t * t;
    const TimeType t3 = t * t2;
    const TimeType tInv = 1.f - t;
    const TimeType tInv2 = tInv * tInv;
    const TimeType tInv3 = tInv * tInv2;

    return position0 * tInv3 + control0 * ( 3 * t * tInv2 ) + control1 * ( 3 * t2 * tInv ) + position1 * t3;
}

/**
 * Evaluates the derivative of a cubic, bezier spline at the specified parametric value.
 * @tparam VectorType The vector type for the spline. Typically vector2f or vector3f.
 * @param position0 The start of the spline
 * @param control0 The first control point of the spline. This point is generally not ON the spline.
 * @param control1 The second control point of the spline. This point is generally not ON the spline.
 * @param position1 The end of the spline
 * @param t The parametric value to evaluate the spline at. Must be in the interval [0, 1].
 * @return The derivative of the cubic bezier spline at the supplied parametric value.
 */
template <typename VectorType>
inline VectorType get_bezier_curve_derivative( const VectorType& position0, const VectorType& control0,
                                               const VectorType& control1, const VectorType& position1, float t ) {
    float tInv = 1.f - t;
    float three_t2 = 3 * t * t;
    float six_t_tInv = 6 * t * tInv;
    float three_tInv2 = 3 * tInv * tInv;

    // This takes 9 mults and 15 adds (for vector3f)
    return ( control0 - position0 ) * three_tInv2 + ( control1 - control0 ) * six_t_tInv +
           ( position1 - control1 ) * three_t2;

    // This takes 12 mults and 9 + 2 adds (for vector3f)
    //	position0 *  -three_tInv2 +
    //	control0  * (-six_t_tInv  + three_tInv2) +
    //	control1  * (-three_t2    + six_t_tInv)  +
    //	position1 *   three_t2;
}

/**
 * Evaluates the second derivative of a cubic, bezier spline at the specified parametric value.
 * @tparam VectorType The vector type for the spline. Typically vector2f or vector3f.
 * @param position0 The start of the spline
 * @param control0 The first control point of the spline. This point is generally not ON the spline.
 * @param control1 The second control point of the spline. This point is generally not ON the spline.
 * @param position1 The end of the spline
 * @param t The parametric value to evaluate the spline at. Must be in the interval [0, 1].
 * @return The second derivative of the cubic bezier spline at the supplied parametric value.
 */
template <typename VectorType>
inline VectorType get_bezier_curve_second_derivative( const VectorType& position0, const VectorType& control0,
                                                      const VectorType& control1, const VectorType& position1,
                                                      float t ) {
    return
        // This is 9 + 3 mults and 6 + 2 adds (for vector3f)
        ( control0 - position0 ) *
            ( 6 * t -
              6 ) + // BUG: This used to be: '(control0 - position0) * (6 -  6 * t) +' which was wrong. Fixed 06/03/2012
        ( control1 - control0 ) * ( 6 - 12 * t ) +
        ( position1 - control1 ) * ( 6 * t );

    // This is 12 + 4 mults and 9 + 3 adds (for vector3f)
    //	position0 * (6.f  * (1.f - t)) +
    //	control0  * (18.f * t - 12.f) +
    //	control1  * (6.f      - 18.f * t)  +
    //	position1 * (6.f  * t);
}

/**
 * Evaluates the derivative of a cubic, bezier spline at the start of the spline.
 * @tparam VectorType The vector type for the spline. Typically vector2f or vector3f.
 * @param position0 The start of the spline
 * @param control0 The first control point of the spline. This point is generally not ON the spline.
 * @param control1 The second control point of the spline. This point is generally not ON the spline.
 * @param position1 The end of the spline
 * @return The derivative of the cubic bezier spline at the start of the spline.
 */
template <typename VectorType>
inline VectorType get_bezier_curve_start_derivative( const VectorType& position0, const VectorType& control0,
                                                     const VectorType& /*control1*/, const VectorType& /*position1*/ ) {
    return 3 * ( control0 - position0 );
}

template <typename VectorType>
inline VectorType get_bezier_curve_mid_derivative( const VectorType& position0, const VectorType& control0,
                                                   const VectorType& control1, const VectorType& position1 ) {
    return 0.75f * ( position1 + control1 - control0 - position0 );
}

template <typename VectorType>
inline VectorType get_bezier_curve_end_derivative( const VectorType& position0, const VectorType& control0,
                                                   const VectorType& control1, const VectorType& position1 ) {
    return 3 * ( position1 - control1 );
}

float bezier_curve_x_to_y( const frantic::graphics2d::vector2f& position0,
                           const frantic::graphics2d::vector2f& control0, const frantic::graphics2d::vector2f& control1,
                           const frantic::graphics2d::vector2f& position1, float x );

// The curvature of a bezier curve at parameter t
template <typename VectorType>
inline float get_bezier_curve_curvature( const VectorType& position0, const VectorType& control0,
                                         const VectorType& control1, const VectorType& position1, float t ) {
    VectorType derivative = get_bezier_curve_derivative( position0, control0, control1, position1, t );
    VectorType secondDerivative = get_bezier_curve_second_derivative( position0, control0, control1, position1, t );

    float derivMag = (float)derivative.get_magnitude_squared();
    return (float)( ( derivative.x * secondDerivative.y - derivative.y * secondDerivative.x ) /
                    ( derivMag * sqrtf( derivMag ) ) );
}

/**
 * Evaluates an esitmate of the arc length for the given curve interval using Gauss-Legendre quadrature
 * @param position0 The start of the curve
 * @param control0 The first control point of the curve. This point is generally not ON the spline.
 * @param control1 The second control point of the curve. This point is generally not ON the spline.
 * @param position1 The end of the curve
 * @param t0 the start of the interval to calculate
 * @param t1 the end of the interval to calculate
 * @return the arc length of the interval
 */
template <typename VectorType>
float get_bezier_curve_arc_length( const VectorType& position0, const VectorType& control0, const VectorType& control1,
                                   const VectorType& position1, float t0, float t1 ) {
    float tDiff = t1 - t0;
    float f0 = float( get_bezier_curve_derivative( position0, control0, control1, position1, t0 ).get_magnitude() );
    float f1 = float(
        get_bezier_curve_derivative( position0, control0, control1, position1, t0 + tDiff * 0.25f ).get_magnitude() );
    float f2 = float(
        get_bezier_curve_derivative( position0, control0, control1, position1, t0 + tDiff * 0.5f ).get_magnitude() );
    float f3 = float(
        get_bezier_curve_derivative( position0, control0, control1, position1, t0 + tDiff * 0.75f ).get_magnitude() );
    float f4 = float( get_bezier_curve_derivative( position0, control0, control1, position1, t1 ).get_magnitude() );
    return ( ( f0 + 4 * ( f1 + f3 ) + 2 * f2 + f4 ) * tDiff ) / 12.f;
}

/**
 * Evaluates an esitmate of the arc length as above, but also returns the derivative of the endpoint
 * this method is used for integrating over a bezier curve's length
 * @param position0 The start of the curve
 * @param control0 The first control point of the curve. This point is generally not ON the spline.
 * @param control1 The second control point of the curve. This point is generally not ON the spline.
 * @param position1 The end of the curve
 * @param t0 the start of the interval to calculate
 * @param t1 the end of the interval to calculate
 * @return the arc length of the interval and the derivative at the endpoint
 */
template <typename VectorType>
std::pair<float, float> get_bezier_curve_arc_length_deriv( const VectorType& position0, const VectorType& control0,
                                                           const VectorType& control1, const VectorType& position1,
                                                           float t0, float t1 ) {
    float tDiff = t1 - t0;
    float f0 = float( get_bezier_curve_derivative( position0, control0, control1, position1, t0 ).get_magnitude() );
    float f1 = float(
        get_bezier_curve_derivative( position0, control0, control1, position1, t0 + tDiff * 0.25f ).get_magnitude() );
    float f2 = float(
        get_bezier_curve_derivative( position0, control0, control1, position1, t0 + tDiff * 0.5f ).get_magnitude() );
    float f3 = float(
        get_bezier_curve_derivative( position0, control0, control1, position1, t0 + tDiff * 0.75f ).get_magnitude() );
    float f4 = float( get_bezier_curve_derivative( position0, control0, control1, position1, t1 ).get_magnitude() );
    return std::make_pair( ( ( f0 + 4 * ( f1 + f3 ) + 2 * f2 + f4 ) * tDiff ) / 12.f, f4 );
}

/**
 * Evaluates an esitmate of the arc length of the entire bezier curve.  Slightly more efficient than
 * specifying an explicit interval
 * @param position0 The start of the curve
 * @param control0 The first control point of the curve. This point is generally not ON the spline.
 * @param control1 The second control point of the curve. This point is generally not ON the spline.
 * @param position1 The end of the curve
 * @return the arc length of the curve
 */
template <typename VectorType>
float get_bezier_curve_arc_length( const VectorType& position0, const VectorType& control0, const VectorType& control1,
                                   const VectorType& position1 ) {
    float f0 = float( get_bezier_curve_start_derivative( position0, control0, control1, position1 ).get_magnitude() );
    float f1 = float( get_bezier_curve_derivative( position0, control0, control1, position1, 0.25f ).get_magnitude() );
    float f2 = float( get_bezier_curve_mid_derivative( position0, control0, control1, position1 ).get_magnitude() );
    float f3 = float( get_bezier_curve_derivative( position0, control0, control1, position1, 0.75f ).get_magnitude() );
    float f4 = float( get_bezier_curve_end_derivative( position0, control0, control1, position1 ).get_magnitude() );
    return ( f0 + 4 * ( f1 + f3 ) + 2 * f2 + f4 ) / 12.f;
}

// The rate of change of the angle with respect to t, at a given parameter t
template <typename VectorType>
inline float get_bezier_curve_dPhidt( const VectorType& position0, const VectorType& control0,
                                      const VectorType& control1, const VectorType& position1, float t ) {
    VectorType derivative = get_bezier_curve_derivative( position0, control0, control1, position1, t );
    VectorType secondDerivative = get_bezier_curve_second_derivative( position0, control0, control1, position1, t );

    return ( derivative.x * secondDerivative.y - derivative.y * secondDerivative.x ) /
           derivative.get_magnitude_squared();
}

// Solves for the t value > than the provided t which is dPhi away from the provided t
template <typename VectorType>
inline float get_bezier_dt_from_dPhi( const VectorType& position0, const VectorType& control0,
                                      const VectorType& control1, const VectorType& position1, float dPhi, float t ) {
    float abs_dPhidt = fabsf( get_bezier_curve_dPhidt( position0, control0, control1, position1, t ) );
    // Limit the step to at most 0.1f
    float dt = 0.1f;
    if( abs_dPhidt > 0.00001f ) {
        dt = dPhi / abs_dPhidt;
        if( dt > 0.1f )
            dt = 0.1f;
    }
    // Check out what the step size will be at the next step
    float abs_dPhidtSecond = fabsf( get_bezier_curve_dPhidt( position0, control0, control1, position1, t + dt ) );
    if( abs_dPhidt > 0.00001f ) {
        float dtSecond = dPhi / abs_dPhidtSecond;
        // If the step where we went would be smaller, take a smaller step.
        if( dtSecond < dt )
            dt = dtSecond;
    }
    // std::cout << "t: " << t << ", dPhi: " << dPhi << ", dt: " << dt << " abs(dPhidt): " << abs_dPhidt << ",
    // curvature: " << get_bezier_curve_curvature( position0, control0, control1, position1, t ) << std::endl;
    return dt;
}

// Uses a curvature-based adaptive approach to convert the bezier curve into a linear polyline.
void convert_bezier_curve_to_samples( const frantic::graphics2d::vector2f& position0,
                                      const frantic::graphics2d::vector2f& control0,
                                      const frantic::graphics2d::vector2f& control1,
                                      const frantic::graphics2d::vector2f& position1, bool includeFirstPoint,
                                      float dPhi, std::vector<std::pair<float, float>>& outSamples );

// Finds the nearest point on the bezier curve
float find_nearest_point_on_bezier_curve( const frantic::graphics2d::vector2f& position0,
                                          const frantic::graphics2d::vector2f& control0,
                                          const frantic::graphics2d::vector2f& control1,
                                          const frantic::graphics2d::vector2f& position1,
                                          const frantic::graphics2d::vector2f& p );

// Given the curve endpoints and two sample points on the curve, solves for the control points
void solve_for_bezier_curve_control_points( const frantic::graphics2d::vector2f& position0,
                                            const frantic::graphics2d::vector2f& position1, float t0,
                                            const frantic::graphics2d::vector2f& t0Position, float t1,
                                            const frantic::graphics2d::vector2f& t1Position,
                                            frantic::graphics2d::vector2f& outControl0,
                                            frantic::graphics2d::vector2f& outControl1 );

// Assuming polyLine.front() and polyLine.back() are the start and end positions,
// Computes the least squares control points corresponding to the t values provided
// for the poly line vertices.
bool best_fit_bezier_center_control_points( const std::vector<frantic::graphics2d::vector2f>& polyLine,
                                            const std::vector<float>& polyLineTValues,
                                            frantic::graphics2d::vector2f& outControl0,
                                            frantic::graphics2d::vector2f& outControl1 );

// Does a least squares style best-fit approximation using a Bezier cubic spline
// Returns true if it is returning an answer, false otherwise.
bool best_fit_bezier_cubic_spline( const std::vector<frantic::graphics2d::vector2f>& polyLine,
                                   frantic::graphics2d::vector2f& outPosition0,
                                   frantic::graphics2d::vector2f& outPosition1,
                                   frantic::graphics2d::vector2f& outControl0,
                                   frantic::graphics2d::vector2f& outControl1 );

/**
A simple template structure for storing a knot on a bezier spline.
To evaluate one of the curve functions for a given pair of knots 'a' and 'b'
position0 = a.position
control0 = a.position + a.outTan
control1 = b.position + a.inTan
position1 = b.position
*/
template <typename VectorType>
struct bezier_curve_point {
    VectorType position;
    VectorType inTan;
    VectorType outTan;

    bezier_curve_point() {}

    bezier_curve_point( const VectorType& _position, const VectorType& _inTan, const VectorType& _outTan )
        : position( _position )
        , inTan( _inTan )
        , outTan( _outTan ) {}
};

// template function object used in 'bezier_curve_integrate_length'
template <typename VectorType>
class integrate_length_solver {
  private:
    integrate_length_solver& operator=( const integrate_length_solver& ) { return *this; }

  public:
    const VectorType& position0;
    const VectorType& control0;
    const VectorType& control1;
    const VectorType& position1;
    const float lowerT;
    const float targetLength;

    integrate_length_solver( const VectorType& inPosition0, const VectorType& inControl0, const VectorType& inControl1,
                             const VectorType& inPosition1, const float inLowerT, const float inTargetLength )
        : position0( inPosition0 )
        , control0( inControl0 )
        , control1( inControl1 )
        , position1( inPosition1 )
        , lowerT( inLowerT )
        , targetLength( inTargetLength ) {}

    float operator()( float t, float& outDerivative ) const {
        std::pair<float, float> result =
            get_bezier_curve_arc_length_deriv( position0, control0, control1, position1, lowerT, t );
        outDerivative = result.second;
        return result.first - targetLength;
    }
};

/**
 * Returns the 't' value along the given bezier curve such that the arc length between
 * it and 'startT' is equal to 'targetLength'.  Repeated calls to this function can
 * generate a set of equidistant sample points along the curve
 * @param position0 The start of the curve
 * @param control0 The first control point of the curve. This point is generally not ON the spline.
 * @param control1 The second control point of the curve. This point is generally not ON the spline.
 * @param position1 The end of the curve
 * @param startT the t value to begin integrating the length
 * @param targetLength the desired arc length to cover
 */
template <typename VectorType>
float bezier_curve_integrate_length( const VectorType& position0, const VectorType& control0,
                                     const VectorType& control1, const VectorType& position1, float startT,
                                     float targetLength ) {
    const float arcLength = get_bezier_curve_arc_length( position0, control0, control1, position1, startT, 1.0f );

    if( targetLength < arcLength ) {
        integrate_length_solver<VectorType> fn( position0, control0, control1, position1, startT, targetLength );

        float guess = startT + ( ( 1.0f - startT ) * targetLength / arcLength );
        return frantic::math::find_root_newton_raphson( fn, startT, 1.0f, guess );
    }

    return 1.0f;
}

/**
 * Generates a set of equidistant samples over the entire bezier spline (assuming that the spline does not 'loop' from
 * the last point to the first) Each knot in the spline will be included as one of the sample points, regardless of the
 * actual distance
 * @param spline the spline to operate on
 * @param targetLength the arclength between each sample
 * @param outSamplePoints after the method will contain a set of pairs, containing the 'segment' and t value of each of
 * the sample points, the segement index represents the curve that spans a given knot and the next knot after it
 */
template <typename VectorType>
void bezier_spline_equidistant_samples( const std::vector<bezier_curve_point<VectorType>>& spline, float targetLength,
                                        bool closed, std::vector<std::pair<size_t, float>>& outSamplePoints ) {
    if( !spline.empty() ) {
        outSamplePoints.clear();

        outSamplePoints.push_back( std::make_pair( 0, 0.0f ) );

        const size_t numIters = ( closed ? spline.size() : spline.size() - 1 );

        for( size_t i = 0; i < numIters; ++i ) {
            const VectorType position0 = spline[i].position;
            const VectorType control0 = spline[i].position + spline[i].outTan;

            // ensure that the points wrap around to the first one for closed splines
            const size_t nextIndex = ( i + 1 ) % spline.size();

            const VectorType control1 = spline[nextIndex].position + spline[nextIndex].inTan;
            const VectorType position1 = spline[nextIndex].position;

            float lastT = 0.0f;

            while( lastT < 1.0f ) {
                lastT = bezier_curve_integrate_length( position0, control0, control1, position1, lastT, targetLength );
                outSamplePoints.push_back( std::make_pair( i, lastT ) );
            }
        }
    }
}

float bezier_curve_x_to_y( const std::vector<bezier_curve_point<frantic::graphics2d::vector2f>>& points, float x );

} // namespace math
} // namespace frantic
