// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

// An empty particle_istream. Useful as a return value for a function that has had an internal error,
// but still needs to produce a particle_istream.
class empty_particle_istream : public particle_istream {
    frantic::channels::channel_map m_particleChannelMap;
    frantic::channels::channel_map m_nativeMap;

  public:
    // NOTE: We don't provide a default constructor, because that's too error-prone.
    empty_particle_istream( const frantic::channels::channel_map& particleChannelMap )
        : m_particleChannelMap( particleChannelMap )
        , m_nativeMap( particleChannelMap ) {}

    // This is the preferred usage of empty_particle_istream
    empty_particle_istream( const frantic::channels::channel_map& particleChannelMap,
                            const frantic::channels::channel_map& nativeMap )
        : m_particleChannelMap( particleChannelMap )
        , m_nativeMap( nativeMap ) {}

    virtual ~empty_particle_istream() {}

    void close() {}

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return _T("Empty Particle Input Stream"); }
    boost::int64_t particle_count() const { return 0; }
    boost::int64_t particle_index() const { return -1; }
    boost::int64_t particle_count_left() const { return 0; }
    boost::int64_t particle_progress_count() const { return 0; }
    boost::int64_t particle_progress_index() const { return -1; }

    const frantic::channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    void set_default_particle( char* /*buffer*/ ) {}

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        m_particleChannelMap = particleChannelMap;
    }

    bool get_particle( char* /*rawParticleBuffer*/ ) { return false; }

    bool get_particles( char* /*particleBuffer*/, std::size_t& numParticles ) {
        numParticles = 0;
        return false;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
