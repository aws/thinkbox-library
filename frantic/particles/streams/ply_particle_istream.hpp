// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/streams/particle_istream.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/geometry/ply_reader.hpp>
#include <frantic/particles/particle_array.hpp>

namespace frantic {
namespace particles {
namespace streams {

class ply_particle_istream : public particle_istream {
  public:
    ply_particle_istream( const frantic::tstring& filename );

    ~ply_particle_istream();

    void close();

    std::size_t particle_size() const;

    frantic::tstring name() const;

    boost::int64_t particle_count() const;
    boost::int64_t particle_index() const;
    boost::int64_t particle_count_left() const;

    boost::int64_t particle_progress_count() const;
    boost::int64_t particle_progress_index() const;

    void set_channel_map( const frantic::channels::channel_map& channelMap );

    const frantic::channels::channel_map& get_channel_map() const;

    const frantic::channels::channel_map& get_native_channel_map() const;

    void set_default_particle( char* rawParticleBuffer );

    bool get_particle( char* rawParticleBuffer );

    bool get_particles( char* rawParticleBuffer, std::size_t& numParticles );

  private:
    bool m_isOpen;
    frantic::particles::particle_array m_particles;
    boost::int64_t m_particleIndex;
    frantic::channels::channel_map_adaptor m_channelMapAdaptor;
    std::vector<char> m_defaultParticle;
    frantic::channels::channel_map m_channelMap;
    frantic::geometry::ply_reader m_reader;
};

} // namespace streams
} // namespace particles
} // namespace frantic
