// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <boost/cstdint.hpp>

#include <frantic/channels/channel_map.hpp>

namespace frantic {
namespace particles {
namespace streams {

using frantic::graphics::color3f;
using frantic::graphics::vector3f;

using frantic::channels::channel_map;

//////////////////////
// PARTICLE OUTPUT STREAM BASE CLASS
//////////////////////

class particle_ostream {
  public:
    particle_ostream() {}

    // Virtual destructor so that we can use allocated pointers (generally with boost::shared_ptr)
    virtual ~particle_ostream() {}

    // This is the particle channel map which specifies the byte layout of the particle structure.
    virtual const channel_map& get_channel_map() const = 0;

    // This allows you to change the particle layout that's being saved on the fly, in case it couldn't
    // be set correctly at creation time.
    virtual void set_channel_map( const channel_map& particleChannelMap ) = 0;

    // This is the size of the particle structure which should be passed to put_particle.  It corresponds to
    // get_channel_map().structure_size()
    virtual std::size_t particle_size() const = 0;

    virtual void close() = 0;

    virtual void put_particle( const char* rawParticleData ) = 0;

    // Convenience function for std::vector
    void put_particle( const std::vector<char>& rawParticleData ) {
        // TODO: In a debug mode, we could confirm that the vector is big enough to hold the particle.
        put_particle( &rawParticleData[0] );
    }
};

typedef boost::shared_ptr<particle_ostream> particle_ostream_ptr;

} // namespace streams

using streams::particle_ostream_ptr;

} // namespace particles
} // namespace frantic
