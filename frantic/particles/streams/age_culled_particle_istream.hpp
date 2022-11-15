// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/streams/particle_istream.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>

#include <tbb/blocked_range.h>
#include <tbb/tbb_allocator.h>

namespace frantic {
namespace particles {
namespace streams {

class age_culled_particle_istream : public frantic::particles::streams::delegated_particle_istream {
    boost::int64_t m_particleIndex;

  public:
    age_culled_particle_istream( particle_istream_ptr pDelegate, float minAge );

    static particle_istream_ptr apply_to_stream( particle_istream_ptr pin, float minAge );

    virtual std::size_t particle_size() const { return this->get_channel_map().structure_size(); }
    virtual boost::int64_t particle_count() const { return -1; }
    virtual boost::int64_t particle_index() const { return m_particleIndex; }
    virtual boost::int64_t particle_count_left() const { return -1; }
    virtual boost::int64_t particle_progress_count() const { return m_delegate->particle_progress_count(); }
    virtual boost::int64_t particle_progress_index() const { return m_delegate->particle_progress_index(); }
    virtual boost::int64_t particle_count_guess() const { return m_delegate->particle_count_guess(); }

    virtual void set_channel_map( const frantic::channels::channel_map& particleChannelMap );
    virtual void set_default_particle( char* rawParticleBuffer );
    virtual const frantic::channels::channel_map& get_channel_map() const;

    virtual bool get_particle( char* pBuffer );
    virtual bool get_particles( char* pBuffer, std::size_t& numParticles );

  private:
    class culling_body;

    typedef std::vector<tbb::blocked_range<std::size_t>, tbb::tbb_allocator<tbb::blocked_range<std::size_t>>>
        vector_type;

    void cull_particles1( char* pBuffer, const tbb::blocked_range<std::size_t>& range,
                          vector_type& outValidRanges ) const;
    void cull_particles2( char* pOutBuffer, const char* pBuffer, const vector_type& outValidRanges,
                          std::size_t& outCount ) const;

    bool cull( char* pParticle ) const;

    void set_channel_map_impl( const frantic::channels::channel_map& channelMap );

  private:
    frantic::channels::channel_map_adaptor m_adaptor;
    frantic::channels::channel_map m_outMap;
    frantic::channels::channel_cvt_accessor<float> m_ageAccessor;

    float m_minAge;
};

} // namespace streams
} // namespace particles
} // namespace frantic
