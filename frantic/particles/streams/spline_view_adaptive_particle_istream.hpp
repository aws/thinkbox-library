// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/graphics/camera.hpp>
#include <frantic/graphics/quat4f.hpp>
#include <frantic/math/iterative_root_solve.hpp>
#include <frantic/math/splines/bezier_spline.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <boost/random.hpp>

// extern frantic::diagnostics::profiling_section
//	psSplineTotal,
//	psComputeLength,
//	psBezierSolve,
//	psRootSolve;

namespace frantic {
namespace particles {
namespace streams {

class bezier_cubic_spline {
    std::vector<frantic::graphics::vector3f> m_points;

  public:
    bezier_cubic_spline() {}

    template <typename Iterator>
    bezier_cubic_spline( Iterator begin, Iterator end ) {
        m_points.assign( begin, end );
    }

    template <typename Iterator>
    void assign( Iterator begin, Iterator end ) {
        m_points.assign( begin, end );
    }

    frantic::graphics::vector3f get_vertex( std::size_t vert ) const { return m_points[3 * vert]; }

    std::size_t num_segments() const { return ( m_points.size() - 1 ) / 3; }

    frantic::graphics::vector3f bezier_interp( std::size_t seg, float t ) const {
        float t2 = t * t;
        float t3 = t * t2;
        float tInv = 1.f - t;
        float tInv2 = tInv * tInv;
        float tInv3 = tInv * tInv2;

        std::size_t index = seg * 3;

        return m_points[index] * tInv3 + m_points[index + 1] * ( 3 * t * tInv2 ) +
               m_points[index + 2] * ( 3 * t2 * tInv ) + m_points[index + 3] * t3;
    }

    frantic::graphics::vector3f bezier_tangent( std::size_t seg, float t ) const {
        float tInv = 1.f - t;
        float three_t2 = 3 * t * t;
        float six_t_tInv = 6 * t * tInv;
        float three_tInv2 = 3 * tInv * tInv;

        std::size_t index = seg * 3;

        return m_points[index] * -three_tInv2 + m_points[index + 1] * ( -six_t_tInv + three_tInv2 ) +
               m_points[index + 2] * ( -three_t2 + six_t_tInv ) + m_points[index + 3] * three_t2;
    }

    // Special case when t = 0
    frantic::graphics::vector3f bezier_tangent_start( std::size_t seg ) const {
        std::size_t index = seg * 3;
        return 3 * ( m_points[index + 1] - m_points[index] );
    }

    // Special case when t = 0.5
    frantic::graphics::vector3f bezier_tangent_mid( std::size_t seg ) const {
        std::size_t index = seg * 3;
        return 0.75f * ( m_points[index + 3] + m_points[index + 2] - m_points[index + 1] - m_points[index] );
    }

    // Special case when t = 1
    frantic::graphics::vector3f bezier_tangent_end( std::size_t seg ) const {
        std::size_t index = seg * 3;
        return 3 * ( m_points[index + 3] - m_points[index + 2] );
    }

    float bezier_arc_length( std::size_t seg ) const {
        float f0 = bezier_tangent_start( seg ).get_magnitude();
        float f1 = bezier_tangent( seg, 0.25f ).get_magnitude();
        float f2 = bezier_tangent_mid( seg ).get_magnitude();
        float f3 = bezier_tangent( seg, 0.75f ).get_magnitude();
        float f4 = bezier_tangent_end( seg ).get_magnitude();
        return ( f0 + 4 * ( f1 + f3 ) + 2 * f2 + f4 ) / 12.f;
    }

    std::pair<float, float> bezier_arc_length( std::size_t seg, float t ) const {
        float f0 = bezier_tangent_start( seg ).get_magnitude();
        float f1 = bezier_tangent( seg, 0.25f * t ).get_magnitude();
        float f2 = bezier_tangent( seg, 0.5f * t ).get_magnitude();
        float f3 = bezier_tangent( seg, 0.75f * t ).get_magnitude();
        float f4 = bezier_tangent( seg, t ).get_magnitude();
        return std::pair<float, float>( ( f0 + 4 * ( f1 + f3 ) + 2 * f2 + f4 ) * t / 12.f, f4 );
    }

    std::pair<float, float> bezier_arc_length( std::size_t seg, float tStart, float tEnd ) const {
        float tDiff = tEnd - tStart;
        float f0 = bezier_tangent( seg, tStart ).get_magnitude();
        float f1 = bezier_tangent( seg, tStart + 0.25f * tDiff ).get_magnitude();
        float f2 = bezier_tangent( seg, tStart + 0.5f * tDiff ).get_magnitude();
        float f3 = bezier_tangent( seg, tStart + 0.75f * tDiff ).get_magnitude();
        float f4 = bezier_tangent( seg, tEnd ).get_magnitude();
        return std::pair<float, float>( ( f0 + 4 * ( f1 + f3 ) + 2 * f2 + f4 ) * tDiff / 12.f, f4 );
    }

    float compute_length() const {
        float curLength = 0.f;
        for( std::size_t i = 0, iEnd = ( m_points.size() - 1 ) / 3; i < iEnd; ++i )
            curLength += bezier_arc_length( i );
        return curLength;
    }

    float compute_lengths( std::vector<float>& outLengths ) const {
        outLengths.resize( ( m_points.size() - 1 ) / 3 );

        float curLength = 0.f;
        for( std::size_t i = 0, iEnd = outLengths.size(); i < iEnd; ++i )
            curLength = outLengths[i] = curLength + bezier_arc_length( i );
        return curLength;
    }

    struct length_fn {
        const bezier_cubic_spline* m_pSpline;
        std::size_t m_splineSeg;
        float m_targetLength;

        length_fn( const bezier_cubic_spline& spline, std::size_t seg, float targetLength )
            : m_pSpline( &spline )
            , m_splineSeg( seg )
            , m_targetLength( targetLength ) {}

        float operator()( float t ) const {
            return m_pSpline->bezier_arc_length( m_splineSeg, t ).first - m_targetLength;
        }

        float operator()( float t, float& outDerivative ) const {
            std::pair<float, float> result = m_pSpline->bezier_arc_length( m_splineSeg, t );
            outDerivative = result.second;
            return result.first - m_targetLength;
        }
    };

    struct integrate_length_solver {
        const bezier_cubic_spline* pSpline;
        std::size_t curSeg;
        float lb;
        float k;

        float operator()( float t, float& outDerivative ) const {
            std::pair<float, float> result = pSpline->bezier_arc_length( curSeg, lb, t );
            outDerivative = result.second;
            return result.first - k;
        }
    };

    std::pair<std::size_t, float> integrate_length( std::pair<std::size_t, float> lb, float desiredLength ) const {
        float arcLength = bezier_arc_length( lb.first, lb.second, 1.f ).first;
        if( desiredLength <= arcLength ) {
            integrate_length_solver fn;
            fn.pSpline = this;
            fn.curSeg = lb.first;
            fn.lb = lb.second;
            fn.k = desiredLength;

            float guess = lb.second + ( 1.f - lb.second ) * desiredLength / arcLength;
            return std::pair<std::size_t, float>(
                lb.first, frantic::math::find_root_newton_raphson( fn, lb.second, 1.f, guess ) );
        } else {
            // We need to go to a different spline in order to get the correct length.
            desiredLength -= arcLength;

            std::size_t curSeg = lb.first + 1;
            std::size_t numSegs = ( m_points.size() - 1 ) / 3;
            for( ; curSeg < numSegs; ++curSeg ) {
                arcLength = bezier_arc_length( curSeg );
                if( desiredLength < arcLength )
                    break;
                desiredLength -= arcLength;
            }
            if( curSeg == numSegs )
                return std::pair<std::size_t, float>( numSegs - 1, 1.f );

            integrate_length_solver fn;
            fn.pSpline = this;
            fn.curSeg = curSeg;
            fn.lb = 0.f;
            fn.k = desiredLength;

            float guess = desiredLength / arcLength;
            return std::pair<std::size_t, float>( curSeg,
                                                  frantic::math::find_root_newton_raphson( fn, 0.f, 1.f, guess ) );
        }
    }

    std::pair<std::size_t, float> compute_length_parameter( const std::vector<float>& segLengths, float desiredLength,
                                                            std::size_t segGuess = 0 ) const {
        typedef std::pair<std::size_t, float> result_t;

        if( desiredLength <= 0 )
            return result_t( 0, 0.f );

        // Find the segment that this point falls into.
        std::size_t segID, iEnd;
        for( segID = segGuess, iEnd = segLengths.size(); segID < iEnd; ++segID ) {
            if( segLengths[segID] > desiredLength )
                break;
        }

        if( segID == iEnd )
            return result_t( segLengths.size() - 1, 1.f );

        float thisLength = segLengths[segID];
        if( segID > 0 ) {
            thisLength -= segLengths[segID - 1];
            desiredLength -= segLengths[segID - 1];
        }

        // frantic::diagnostics::scoped_profile spRootSolve(psRootSolve);

        // Make a guess by assuming this is a straight line.
        float guess = desiredLength / thisLength;

        length_fn fn( *this, segID, desiredLength );
        return std::make_pair( segID, frantic::math::find_root_newton_raphson( fn, 0.f, 1.f, guess, 1e-3f ) );
    }
};

namespace detail {
struct find_arc_distance_fn {
    const bezier_cubic_spline* pSpline;
    std::size_t curSeg;
    float lb;
    float k;

    float operator()( float t, float& outDerivative ) const {
        std::pair<float, float> result = pSpline->bezier_arc_length( curSeg, lb, t );
        outDerivative = result.second;
        return result.first - k;
    }
};
} // namespace detail

class spline_particle_istream : public particle_istream {
    std::vector<bezier_cubic_spline> m_curvePoints;
    std::vector<bezier_cubic_spline> m_curvePointsOffsetTime; // For motion blur
    std::vector<bezier_cubic_spline>
        m_curvePointsReferenceTime; // For providing a time-consistent XYZ coord. Probably for camera mapping.

    // The time difference in seconds between the curves in m_curvePoints, and m_curvePointsOffsetTime
    float m_timeStep;

    // Spacing in pixels between particles.
    float m_pixelStep;

    // The world-space distance along the curve between particles.
    float m_distanceStep;

    std::size_t m_curSpline; // Spline index of the current spline being seeded on.
    float m_curSplineLength; // Length of current spline.

    struct {
        std::size_t segment; // Segment index of previous particle
        float t;             // Parametric value of previous particle
        float length;        // Distance along curve of previous particle

        frantic::graphics::quat4f orientation; // Orientation of previous particle

        bool pixelPosValid;                     // Pixel position of previous particle is valid
        frantic::graphics2d::vector2f pixelPos; // Pixel position of previous particle
    } m_prevParticle;

    // Camera for view-adaptive generation
    boost::shared_ptr<const frantic::graphics::camera<float>> m_pCamera;

    // Random generator for randomizing base orientation
    boost::mt19937 m_rngGen;
    boost::uniform_real<> m_rngRange;
    boost::variate_generator<boost::mt19937&, boost::uniform_real<>> m_rng;

    boost::int64_t m_particleIndex;
    frantic::channels::channel_map m_outMap;
    frantic::channels::channel_map m_nativeMap;
    boost::scoped_array<char> m_defaultParticle;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_velocityAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_refPosAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_tangentAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_normalAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_hairRootPositionAccessor;
    frantic::channels::channel_cvt_accessor<float> m_hairLengthAccessor;
    frantic::channels::channel_cvt_accessor<float> m_distanceAccessor;
    frantic::channels::channel_cvt_accessor<float> m_densityAccessor;
    frantic::channels::channel_cvt_accessor<int> m_splineIDAccessor;

  private:
    void init_spline() {
        using frantic::graphics::vector3f;
        using frantic::graphics2d::vector2f;

        m_prevParticle.segment = 0;
        m_prevParticle.t = 0.f;
        m_prevParticle.length = 0.f;
        m_prevParticle.orientation =
            frantic::graphics::quat4f::from_angle_axis( (float)m_rng(), vector3f::from_zaxis() );

        vector3f tangent = vector3f::normalize( m_curvePoints[m_curSpline].bezier_tangent_start( 0 ) );
        vector3f crossUp = vector3f::normalize( vector3f::cross( vector3f::from_zaxis(), tangent ) );

        float cosAngle = tangent.z;
        if( fabsf( cosAngle ) < 0.9999f )
            m_prevParticle.orientation =
                frantic::graphics::quat4f::from_angle_axis( acos( cosAngle ), crossUp ) * m_prevParticle.orientation;

        if( m_pCamera.get() ) {
            vector3f basePos = m_curvePoints[m_curSpline].get_vertex( 0 );

            m_prevParticle.pixelPosValid = true;
            m_prevParticle.pixelPos = m_pCamera->from_worldspace_position( basePos, m_prevParticle.pixelPosValid );
        } else {
            m_prevParticle.pixelPosValid = false;
        }

        // This is approximate length only.
        // TODO: a more accurate computation
        m_curSplineLength = m_curvePoints[m_curSpline].compute_length();
    }

    struct find_pixel_dist_solver {
        const frantic::graphics::camera<float>* pCamera;
        frantic::graphics2d::vector2f refPoint;

        float k;
        std::size_t curSeg;
        const bezier_cubic_spline* pSpline;

        float operator()( float t ) const {
            bool isValid = true;
            frantic::graphics::vector3f p = pSpline->bezier_interp( curSeg, t );
            frantic::graphics2d::vector2f sp = pCamera->from_worldspace_position( p, isValid );

            return frantic::graphics2d::vector2f::distance( sp, refPoint ) - k;
        }
    };

    static float find_arc_distance_in_range( const bezier_cubic_spline& spline, std::size_t segment, float lb, float ub,
                                             float guess, float target, float epsilon ) {
        detail::find_arc_distance_fn fn;

        fn.pSpline = &spline;
        fn.curSeg = segment;
        fn.lb = lb;
        fn.k = target;

        return frantic::math::find_root_newton_raphson( fn, lb, ub, guess, epsilon );
    }

    void apply_minimum_step_constraint( std::size_t& inoutSegment, float& inoutT, float& inoutDistance ) const {
        static const float EPSILON = 1e-4f;

        if( m_prevParticle.segment == inoutSegment && inoutT - m_prevParticle.t < EPSILON ) {
            inoutT = ( m_prevParticle.t + EPSILON );
            if( inoutT >= 1.f ) {
                if( inoutSegment < m_curvePoints[m_curSpline].num_segments() - 1 ) {
                    inoutT = 0.f;
                    inoutSegment += 1;
                } else
                    inoutT = 1.f;
                inoutDistance =
                    m_curvePoints[m_curSpline].bezier_arc_length( m_prevParticle.segment, m_prevParticle.t, 1.f ).first;
            } else
                inoutDistance = m_curvePoints[m_curSpline]
                                    .bezier_arc_length( m_prevParticle.segment, m_prevParticle.t, inoutT )
                                    .first;
        }
    }

    void apply_distance_constraint( std::size_t& inoutSegment, float& inoutT, float& outDistance ) const {
        using frantic::graphics::vector3f;
        using frantic::graphics2d::vector2f;

        static const float EPSILON = 1e-2f;

        float resultT = inoutT;
        std::size_t seg = inoutSegment;

        // We have chosen a point which is m_pixelSpacing away. Now make sure that the world distance constraint is met
        // as well.
        float distance;
        if( seg == m_prevParticle.segment ) {
            // The upper bound on the distance search lies in the same spline segment as the lower bound, so we only
            // need to inspect this single segment.
            distance =
                m_curvePoints[m_curSpline].bezier_arc_length( m_prevParticle.segment, m_prevParticle.t, resultT ).first;
            if( distance > m_distanceStep ) {
                // Guess that this is a straight line to start the root-solve.
                float guess = m_prevParticle.t + ( resultT - m_prevParticle.t ) * m_distanceStep / distance;
                resultT = find_arc_distance_in_range( m_curvePoints[m_curSpline], m_prevParticle.segment,
                                                      m_prevParticle.t, resultT, guess, m_distanceStep, EPSILON );

                inoutT = resultT;
                inoutSegment = m_prevParticle.segment;
                outDistance = m_distanceStep;
            } else
                outDistance = distance;
        } else {
            // The upper bound on the distance search lies in a different segment from the lower bound so we have three
            // parts. The region [lb,1] on the first segment, the regions [0,1] on in-between segments, and [0,ub] on
            // the final segment.
            distance =
                m_curvePoints[m_curSpline].bezier_arc_length( m_prevParticle.segment, m_prevParticle.t, 1.f ).first;
            if( distance > m_distanceStep ) {
                // Guess that this is a straight line to start the root-solve.
                float guess = m_prevParticle.t + ( 1.f - m_prevParticle.t ) * m_distanceStep / distance;
                resultT = find_arc_distance_in_range( m_curvePoints[m_curSpline], m_prevParticle.segment,
                                                      m_prevParticle.t, 1.f, guess, m_distanceStep, EPSILON );

                inoutT = resultT;
                inoutSegment = m_prevParticle.segment;
                outDistance = m_distanceStep;
            } else {
                bool done = false;
                for( std::size_t i = m_prevParticle.segment + 1; i < seg && !done; ++i ) {
                    float thisSegmentTarget = m_distanceStep - distance;
                    float thisSegmentDistance = m_curvePoints[m_curSpline].bezier_arc_length( i );

                    distance += thisSegmentDistance;
                    if( distance > m_distanceStep ) {
                        // Guess that this is a straight line to start the root-solve.
                        float guess = thisSegmentTarget / thisSegmentDistance;
                        resultT = find_arc_distance_in_range( m_curvePoints[m_curSpline], i, 0.f, 1.f, guess,
                                                              thisSegmentTarget, EPSILON );

                        inoutT = resultT;
                        inoutSegment = i;
                        outDistance = m_distanceStep;
                        done = true;
                    }
                }
                if( !done ) {
                    float thisSegmentTarget = m_distanceStep - distance;
                    float thisSegmentDistance = m_curvePoints[m_curSpline].bezier_arc_length( seg, 0.f, resultT ).first;

                    distance += thisSegmentDistance;
                    if( distance > m_distanceStep ) {
                        // Guess that this is a straight line to start the root-solve.
                        float guess = resultT * thisSegmentTarget / thisSegmentDistance;
                        resultT = find_arc_distance_in_range( m_curvePoints[m_curSpline], seg, 0.f, resultT, guess,
                                                              thisSegmentTarget, EPSILON );

                        inoutT = resultT;
                        inoutSegment = seg;
                        outDistance = m_distanceStep;
                    } else
                        outDistance = distance;
                } // if( !done ){
            }     // if( distance > m_distanceStep ){
        }         // if( seg == m_prevParticle.segment ){}else{
    }

    void find_next_point( std::size_t& outSegment, float& outT, float& outDistance ) const {
        using frantic::graphics::vector3f;
        using frantic::graphics2d::vector2f;

        // If the previous particle can not be placed on the screen correctly, just apply the distance constraint.
        if( !m_pCamera.get() || !m_prevParticle.pixelPosValid ) {
            outT = 1.f;
            outSegment = m_prevParticle.segment;

            apply_distance_constraint( outSegment, outT, outDistance );
            apply_minimum_step_constraint( outSegment, outT, outDistance );
            return;
        }

        vector3f nextVert;
        vector2f nextVertPixel;
        float pixelDistance;

        std::size_t seg = m_prevParticle.segment;
        float lb = m_prevParticle.t;
        float ub = 1.f;

        do {
            bool isValid = true;
            nextVert = m_curvePoints[m_curSpline].get_vertex( seg + 1 );
            nextVertPixel = m_pCamera->from_worldspace_position( nextVert, isValid );

            // If this next particle can not be placed on the screen correctly, just apply the distance constraint.
            if( !isValid ) {
                outT = 1.f;
                outSegment = seg;

                apply_distance_constraint( outSegment, outT, outDistance );
                apply_minimum_step_constraint( outSegment, outT, outDistance );
                return;
            }

            pixelDistance = vector2f::distance( m_prevParticle.pixelPos, nextVertPixel );

            if( pixelDistance >= m_pixelStep ) {
                find_pixel_dist_solver fn;
                fn.pCamera = m_pCamera.get();
                fn.refPoint = m_prevParticle.pixelPos;
                fn.k = m_pixelStep;
                fn.curSeg = seg;
                fn.pSpline = &m_curvePoints[m_curSpline];

                float resultT = frantic::math::find_root_bisection( fn, lb, ub, 1e-1f );

                outSegment = seg;
                outT = resultT;

                apply_distance_constraint( outSegment, outT, outDistance );
                apply_minimum_step_constraint( outSegment, outT, outDistance );
                return;
            }

            ++seg;
            lb = 0.f;
            ub = 1.f;
        } while( seg < m_curvePoints[m_curSpline].num_segments() );

        // By reaching here, this means we have gone to the end of the spline w/o choosing a particle.
        outSegment = m_curvePoints[m_curSpline].num_segments() - 1;
        outT = 1.f;

        apply_distance_constraint( outSegment, outT, outDistance );
        apply_minimum_step_constraint( outSegment, outT, outDistance );
    }

  public:
    /**
     * constructor input:
     *	the constructor call's first 3 params are vectors of "splines" where each "spline" is a list of the cubic bezier
     *control points -param1 is the actual hair spline at the current time -param2 is used for doing finite differences
     *to compute the velocity -param3 is the spline at some reference time. used for camera mapping -it fills in a
     *float32[3] channel called "ReferencePosition" that is the equivalent point on the spline, except from the
     *reference splines you would use that for shading that uses the position of the particle that way it doesn't chane
     *color if the hair is animated. -timeStepSeconds is the time difference in seconds between the spline and its
     *time-offset version for finite differencing -distanceStep is the maximum distance in object space units between
     *each particle seeded on the spline -pCamera is a camera object (can be null if not doing camera dependent seeding)
     *-The camera is used so that particles can be seeded at a constant screen resolution for scenes with a lot of
     *camera motion, its difficult to pick a world unit spacing that will make it look good and be fast and not have
     *gaps -spacing is a maximum spacing (in pixels) between particles if pCamera is non-NULL -this parameter is ignored
     *if pCamera is NULL.
     *
     * output state of spline_particle_istream:
     *	-the channels exposed by the stream are: Position, Velocity, ReferencePosition, Tangent, Normal, HairRoot,
     *HairLength, Distance, and Density -Density is Krakatoa's hacky version of "thickness". in 3dsmax there is a
     *spinner that affects how quickly the hair tapered to 0 density -using Density = 1 * (Distance / HairLength) ^ X
     *where X is a user parameter defaulting to 10 -HairLength is the length of the spline -Distance is the length along
     *the spline of the current hair
     */
    spline_particle_istream( const std::vector<std::vector<frantic::graphics::vector3f>>& splines,
                             const std::vector<std::vector<frantic::graphics::vector3f>>& splinesOffsetTime,
                             const std::vector<std::vector<frantic::graphics::vector3f>>& splinesReferenceTime,
                             float timeStepSeconds, float distanceStep,
                             boost::shared_ptr<const frantic::graphics::camera<float>> pCamera, float spacing )
        : m_rngRange( 0, 2 * M_PI )
        , m_rng( m_rngGen, m_rngRange ) {
        using namespace frantic::graphics;
        using namespace frantic::graphics2d;

        m_curvePoints.resize( splines.size() );
        m_curvePointsOffsetTime.resize( splinesOffsetTime.size() );
        m_curvePointsReferenceTime.resize( splinesReferenceTime.size() );

        if( m_curvePointsOffsetTime.size() != 0 && m_curvePointsOffsetTime.size() != m_curvePoints.size() )
            throw std::runtime_error(
                "spline_particle_istream() - The spline counts at different time steps did not match" );

        if( m_curvePointsReferenceTime.size() != 0 && m_curvePointsReferenceTime.size() != m_curvePoints.size() )
            throw std::runtime_error(
                "spline_particle_istream() - The spline count at the reference time did not match the current time" );

        std::vector<std::vector<frantic::graphics::vector3f>>::const_iterator it, itEnd;
        std::vector<bezier_cubic_spline>::iterator outIt;

        it = splines.begin();
        itEnd = splines.end();
        outIt = m_curvePoints.begin();
        for( ; it != itEnd; ++it, ++outIt )
            outIt->assign( it->begin(), it->end() );

        it = splinesOffsetTime.begin();
        itEnd = splinesOffsetTime.end();
        outIt = m_curvePointsOffsetTime.begin();
        for( ; it != itEnd; ++it, ++outIt )
            outIt->assign( it->begin(), it->end() );

        it = splinesReferenceTime.begin();
        itEnd = splinesReferenceTime.end();
        outIt = m_curvePointsReferenceTime.begin();
        for( ; it != itEnd; ++it, ++outIt )
            outIt->assign( it->begin(), it->end() );

        m_pCamera = pCamera;
        m_pixelStep = spacing;
        m_distanceStep = distanceStep;
        m_timeStep = timeStepSeconds;
        m_curSpline = 0;

        // TODO: Un-hardcode this.
        m_rngGen.seed( 42 );

        // Skip any empty splines
        while( m_curSpline < m_curvePoints.size() && m_curvePoints[m_curSpline].num_segments() == 0 )
            ++m_curSpline;

        if( m_curSpline < m_curvePoints.size() )
            init_spline();

        m_particleIndex = -1;

        m_nativeMap.define_channel<frantic::graphics::vector3f>( _T("Position") );
        m_nativeMap.define_channel<frantic::graphics::vector3f>( _T("Velocity") );
        m_nativeMap.define_channel<frantic::graphics::vector3f>( _T("ReferencePosition") );
        m_nativeMap.define_channel<frantic::graphics::vector3f>( _T("Tangent") );
        m_nativeMap.define_channel<frantic::graphics::vector3f>( _T("Normal") );
        m_nativeMap.define_channel<frantic::graphics::vector3f>( _T("HairRoot") );
        m_nativeMap.define_channel<float>( _T("HairLength") );
        m_nativeMap.define_channel<float>( _T("Distance") );
        m_nativeMap.define_channel<float>( _T("Density") );
        m_nativeMap.define_channel<int>( _T("SplineID") );
        m_nativeMap.end_channel_definition();

        set_channel_map( m_nativeMap );
    }

    // Virtual destructor so that we can use allocated pointers (generally with boost::shared_ptr)
    virtual ~spline_particle_istream() { close(); }

    virtual void close() {}

    // The stream can return its filename or other identifier for better error messages.
    virtual frantic::tstring name() const { return _T("spline_stream"); }

    // This is the size of the particle structure which will be loaded, in bytes.
    virtual std::size_t particle_size() const { return m_outMap.structure_size(); }

    // Returns the number of particles, or -1 if unknown
    virtual boost::int64_t particle_count() const { return -1; }
    virtual boost::int64_t particle_index() const { return m_particleIndex; }
    virtual boost::int64_t particle_count_left() const { return -1; }

    virtual boost::int64_t particle_progress_count() const { return m_curvePoints.size(); }
    virtual boost::int64_t particle_progress_index() const { return m_curSpline; }

    virtual boost::int64_t particle_count_guess() const {
        // Make a guess based on the distance in a curve and spacing per particle. This isn't really accurate
        // when doing view-adaptive particles, but that would be pretty expensive to figure out.
        boost::int64_t result = 0;
        std::vector<bezier_cubic_spline>::const_iterator it, itEnd;
        for( it = m_curvePoints.begin(), itEnd = m_curvePoints.end(); it != itEnd; ++it )
            result += static_cast<boost::int64_t>( it->compute_length() / m_distanceStep + 0.5f ) + 1;
        return result;
    }

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    virtual void set_channel_map( const frantic::channels::channel_map& pcm ) {
        { // Swap in a new default particle.
            boost::scoped_array<char> newDefault( new char[pcm.structure_size()] );
            memset( newDefault.get(), 0, pcm.structure_size() );
            if( m_defaultParticle ) {
                frantic::channels::channel_map_adaptor tempAdaptor( pcm, m_outMap );
                tempAdaptor.copy_structure( newDefault.get(), m_defaultParticle.get() );
            }
            m_defaultParticle.swap( newDefault );
        }

        m_outMap = pcm;
        m_posAccessor = m_outMap.get_accessor<frantic::graphics::vector3f>( _T("Position") );

        m_velocityAccessor.reset();
        m_refPosAccessor.reset();
        m_tangentAccessor.reset();
        m_normalAccessor.reset();
        m_hairLengthAccessor.reset();
        m_hairRootPositionAccessor.reset();
        m_distanceAccessor.reset();
        m_densityAccessor.reset();
        m_splineIDAccessor.reset();

        if( m_outMap.has_channel( _T("Velocity") ) && m_curvePointsOffsetTime.size() > 0 )
            m_velocityAccessor = m_outMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Velocity") );

        if( m_outMap.has_channel( _T("ReferencePosition") ) && m_curvePointsReferenceTime.size() > 0 )
            m_refPosAccessor = m_outMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("ReferencePosition") );

        if( m_outMap.has_channel( _T("Tangent") ) )
            m_tangentAccessor = m_outMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Tangent") );

        if( m_outMap.has_channel( _T("Normal") ) )
            m_normalAccessor = m_outMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Normal") );

        if( m_outMap.has_channel( _T("HairRoot") ) )
            m_hairRootPositionAccessor = m_outMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("HairRoot") );

        if( m_outMap.has_channel( _T("Distance") ) )
            m_distanceAccessor = m_outMap.get_cvt_accessor<float>( _T("Distance") );

        if( m_outMap.has_channel( _T("HairLength") ) )
            m_hairLengthAccessor = m_outMap.get_cvt_accessor<float>( _T("HairLength") );

        if( m_outMap.has_channel( _T("Density") ) )
            m_densityAccessor = m_outMap.get_cvt_accessor<float>( _T("Density") );

        if( m_outMap.has_channel( _T("SplineID") ) )
            m_splineIDAccessor = m_outMap.get_cvt_accessor<int>( _T("SplineID") );
    }

    // This is the particle channel map which specifies the byte layout of the particle structure that is being used.
    virtual const frantic::channels::channel_map& get_channel_map() const { return m_outMap; }

    // This is the particle channel map which specifies the byte layout of the input to this stream.
    // NOTE: This value is allowed to change after the following conditions:
    //    * set_channel_map() is called (for example, the empty_particle_istream equates the native map with the
    //    external map)
    //    * get_particle() is called (for example, a concatenated particle stream will switch to a different stream for
    //    the next particle)
    virtual const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    /** This provides a default particle which should be used to fill in channels of the requested channel map
     *	which are not supplied by the native channel map.
     *	IMPORTANT: Make sure the buffer you pass in is at least as big as particle_size() bytes.
     */
    virtual void set_default_particle( char* rawParticleBuffer ) {
        boost::scoped_array<char> newDefault( new char[m_outMap.structure_size()] );
        m_outMap.copy_structure( newDefault.get(), rawParticleBuffer );
        m_defaultParticle.swap( newDefault );
    }

    // This reads a particle into a buffer matching the channel_map.
    // It returns true if a particle was read, false otherwise.
    // IMPORTANT: Make sure the buffer you pass in is at least as big as particle_size() bytes.
    virtual bool get_particle( char* rawParticleBuffer ) {
        using namespace frantic::graphics;
        using namespace frantic::graphics2d;

        if( m_curSpline == m_curvePoints.size() )
            return false;

        // frantic::diagnostics::scoped_profile spSplineTotal(psSplineTotal);

        ++m_particleIndex;

        // Copy the default values into the particle
        m_outMap.copy_structure( rawParticleBuffer, m_defaultParticle.get() );

        // psBezierSolve.enter();

        std::size_t nextSegment;
        float nextT;
        float distance;

        find_next_point( nextSegment, nextT, distance );

        // psBezierSolve.exit();

        vector3f pos = m_curvePoints[m_curSpline].bezier_interp( nextSegment, nextT );

        // Update the orientation
        transform4f tm;
        m_prevParticle.orientation.as_transform4f( tm );

        vector3f tangent = vector3f::normalize( m_curvePoints[m_curSpline].bezier_tangent( nextSegment, nextT ) );
        vector3f oldTangent = tm.get_column( 2 );
        vector3f parallelTransportAxis = vector3f::normalize( vector3f::cross( oldTangent, tangent ) );

        float cosAngle = vector3f::dot( oldTangent, tangent );
        if( fabsf( cosAngle ) < ( 0.9999f ) ) {
            m_prevParticle.orientation =
                quat4f::from_angle_axis( acos( cosAngle ), parallelTransportAxis ) * m_prevParticle.orientation;
            m_prevParticle.orientation.as_transform4f( tm );
        }

        m_prevParticle.length += distance;
        m_prevParticle.segment = nextSegment;
        m_prevParticle.t = nextT;

        if( m_pCamera.get() ) {
            bool isValid = true;
            vector2f posScreen = m_pCamera->from_worldspace_position( pos, isValid );

            m_prevParticle.pixelPos = posScreen;

            // Check the view frustrum as well, since we don't need to be adaptive outside the viewing area.
            if( isValid &&
                ( posScreen.x < 0.f || posScreen.y < 0.f || posScreen.x >= m_pCamera->get_output_size().xsize ||
                  posScreen.y >= m_pCamera->get_output_size().ysize ) )
                isValid = false;
            m_prevParticle.pixelPosValid = isValid;
        }

        m_posAccessor.get( rawParticleBuffer ) = pos;

        if( m_velocityAccessor.is_valid() ) {
            vector3f timeOffsetPos = m_curvePointsOffsetTime[m_curSpline].bezier_interp( nextSegment, nextT );
            vector3f velocity = ( pos - timeOffsetPos ) / m_timeStep;
            m_velocityAccessor.set( rawParticleBuffer, velocity );
        }

        if( m_refPosAccessor.is_valid() ) {
            vector3f refPos = m_curvePointsReferenceTime[m_curSpline].bezier_interp( nextSegment, nextT );
            m_refPosAccessor.set( rawParticleBuffer, refPos );
        }

        if( m_normalAccessor.is_valid() )
            m_normalAccessor.set( rawParticleBuffer, tm.get_column( 0 ) );

        if( m_tangentAccessor.is_valid() )
            m_tangentAccessor.set( rawParticleBuffer, tm.get_column( 2 ) );

        if( m_hairRootPositionAccessor.is_valid() )
            m_hairRootPositionAccessor.set( rawParticleBuffer, m_curvePoints[m_curSpline].get_vertex( 0 ) );

        if( m_distanceAccessor.is_valid() )
            m_distanceAccessor.set( rawParticleBuffer, std::min( m_prevParticle.length, m_curSplineLength ) );

        if( m_hairLengthAccessor.is_valid() )
            m_hairLengthAccessor.set( rawParticleBuffer, m_curSplineLength );

        if( m_densityAccessor.is_valid() )
            m_densityAccessor.set( rawParticleBuffer, distance );

        int splineID = (int)m_curSpline;
        if( m_splineIDAccessor.is_valid() )
            m_splineIDAccessor.set( rawParticleBuffer, splineID );

        if( nextSegment == m_curvePoints[m_curSpline].num_segments() - 1 && ( 1.f - nextT ) < 1e-4f ) {
            do {
                ++m_curSpline;
            } while( m_curSpline < m_curvePoints.size() && m_curvePoints[m_curSpline].num_segments() == 0 );

            if( m_curSpline < m_curvePoints.size() )
                init_spline();
        }

        return true;
    }

    // This reads a group of particles. Returns false if the end of the source
    // was reached during the read.
    virtual bool get_particles( char* rawParticleBuffer, std::size_t& numParticles ) {
        std::size_t particleSize = particle_size();
        for( std::size_t i = 0; i < numParticles; ++i ) {
            if( !get_particle( rawParticleBuffer ) ) {
                numParticles = i;
                return false;
            }
            rawParticleBuffer += particleSize;
        }
        return true;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
