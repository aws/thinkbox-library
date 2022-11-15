// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/graphics/quat4f.hpp>
#include <frantic/math/splines/bezier_spline.hpp>
#include <frantic/math/utils.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/age_culled_particle_istream.hpp>
#include <frantic/particles/streams/concatenated_particle_istream.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/particles/streams/shared_particle_container_particle_istream.hpp>
#include <frantic/sort/sort.hpp>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace frantic {
namespace particles {
namespace streams {

namespace detail {
void make_other_particles_map( frantic::channels::channel_map& outMap, const frantic::channels::channel_map& map );
}

/**
 * This class is responsible for using two samples of a particle system in order to smoothly interpolate the particle's
 * motion at times between the two samples. It can do cubic or linear interpolation, or velocity extrapolation if there
 * is no matching particle in the second sample.
 *
 * @note There must be an ID channel in both streams to do proper matching of particles. Maybe the index could be
 * substituted in certain situations.
 */
class time_interpolation_particle_istream : public delegated_particle_istream {
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_streamPosAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_streamVelAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_streamNormAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_streamTangentAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_streamColorAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_streamTextureCoordAccessor;
    frantic::channels::channel_cvt_accessor<float> m_streamAgeAccessor;
    frantic::channels::channel_cvt_accessor<float> m_streamDensityAccessor;
    frantic::channels::channel_static_cast_const_accessor<boost::int64_t> m_streamIDAccessor;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_otherPosAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_otherVelAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_otherNormAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_otherTangentAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_otherColorAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_otherTextureCoordAccessor;
    frantic::channels::channel_cvt_accessor<float> m_otherAgeAccessor;
    frantic::channels::channel_cvt_accessor<float> m_otherDensityAccessor;
    frantic::channels::channel_static_cast_const_accessor<boost::int64_t> m_otherIDAccessor;

    float m_interpVal;       // Interpolation parameter in [0,1].
    float m_timeStepSeconds; // Amount of time between samples
    float m_offsetSeconds;   // = (m_interpVal * m_timeStepSeconds)

    // If neither the channel map or native channel map has an "ID" channel we want to force extrapolation based on
    // velocity.
    bool m_forceExtrapolation;

    frantic::channels::channel_map m_outMap;
    frantic::channels::channel_map_adaptor m_outAdaptor;

    frantic::particles::particle_array m_otherParticles;

    void set_channel_map_internal( const frantic::channels::channel_map& pcm );
    void set_other_accessors();
    const char* find_particle_match( boost::int64_t id, bool& outSuccess ) const;
    void process_particle( char* particle ) const;
    void process_particles_in_place( char* pBuffer, const tbb::blocked_range<std::size_t>& range ) const;
    void process_particles( char* pDestBuffer, char* pSrcBuffer, const tbb::blocked_range<std::size_t>& range ) const;

    static void operator_delete( void* ptr ) { operator delete( ptr ); }

    template <class T>
    class particle_compare_id {
        frantic::channels::channel_accessor<T> m_idAccessor;

      public:
        particle_compare_id( frantic::channels::channel_map map ) { m_idAccessor = map.get_accessor<T>( _T("ID") ); }

        bool operator()( const char* lhs, const char* rhs ) const {
            return m_idAccessor.get( lhs ) < m_idAccessor.get( rhs );
        }

        bool operator()( const char* lhs, const T& rhs ) const { return m_idAccessor.get( lhs ) < rhs; }

        bool operator()( const T& lhs, const char* rhs ) const { return lhs < m_idAccessor.get( rhs ); }

        T operator()( const char* p ) const { return m_idAccessor.get( p ); }
    };

    void sort_other_particles();

  public:
    /**
     * Constructor
     * @param mainPin The stream containing the particles that are "closest" to the current time
     * @param otherTimePin The stream containing the same particles, sampled at a different time
     * @param timeStepSeconds The time distance between samples in seconds. Should be negative if otherTimePin is
     * "before" mainPin.
     * @param interpVal The [0,1] interpolation position between the two samples.
     */
    time_interpolation_particle_istream( boost::shared_ptr<particle_istream> mainPin,
                                         boost::shared_ptr<particle_istream> otherTimePin, float timeStepSeconds,
                                         float interpVal );

    // called by apply_interpolation_particle_istream_with_lifespan_culling
    time_interpolation_particle_istream( boost::shared_ptr<particle_istream> mainPin,
                                         frantic::particles::particle_array interpParticles, float timeStepSeconds,
                                         float interpVal );

    virtual void set_channel_map( const frantic::channels::channel_map& pcm );

    virtual const frantic::channels::channel_map& get_channel_map() const;

    virtual void set_default_particle( char* defaultParticle );

    virtual bool get_particle( char* outParticle );

    virtual bool get_particles( char* buffer, std::size_t& numParticles );
};

/**
 * Given the closest sample (ie. 'pin') and closest sample from the other side (ie. 'otherPin') this function will
 * create a stream that interpolates the particles to a given time. It will use smooth, cubic interpolation if the
 * particles have a Velocity channel. Particle birth & death times are also interpolated such that there is no
 * objectionable popping in or out of particles. \param pin A stream of particles closest in time to the requested
 * interpolation. \param otherPin The stream of particles closest on the other side of the requested time. This is the
 * second closest sample when using uniform sample times. \param timeStepSeconds The amount of time between samples.
 * Will  e negative if 'otherPin' comes before 'pin', or positive otherwise. \param interpParam The [0, 1] interpolation
 * parameter between the two samples. Must be positive regardless of the ordering in time of 'pin' and 'otherPin'.
 * \return A stream that cubically interpolates the two samples.
 */
boost::shared_ptr<particle_istream>
apply_interpolation_particle_istream_with_lifespan_culling( boost::shared_ptr<particle_istream> pin,
                                                            boost::shared_ptr<particle_istream> otherPin,
                                                            double timeStepSeconds, double interpParam );

} // namespace streams
} // namespace particles
} // namespace frantic
