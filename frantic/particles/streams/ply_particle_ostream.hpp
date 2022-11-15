// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/streams/particle_ostream.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/particle_array.hpp>

namespace frantic {
namespace particles {
namespace streams {

class ply_particle_ostream : public particle_ostream {
  public:
    /**
     * @param filename .ply file to write to.
     * @param channelMap channel layout of the buffer that will be passed to
     *    put_particle().
     * @param diskChannelMap channels to write to disk.  Only the Position,
     *    Color, TextureCoord, and Normal channels are saved.  These channels
     *    are always saved as float32[3]; their data type in the channel map
     *    is ignored.  All other channels are ignored.
     */
    ply_particle_ostream( const frantic::tstring& filename, const frantic::channels::channel_map& channelMap,
                          const frantic::channels::channel_map& diskChannelMap );

    ~ply_particle_ostream();

    const channel_map& get_channel_map() const;

    void set_channel_map( const channel_map& channelMap );

    std::size_t particle_size() const;

    /**
     *  Close the stream, and write all of the particles to disk.
     *
     *  Unlike most others, this stream does not keep a file open while you
     * insert particles.  Instead, the particles are added to a buffer in
     * memory, and they are written to disk by this close() function.  I do
     * this because the PLY IO library we're using requires the particle count
     * before we begin to write the particle data.  (As an alternative, I
     * think we could write our own IO code, reserve space for the particle
     * count using a comment, and write the particle count when we close the
     * file.)
     *
     * @note No particles will be written to disk unless you call close().
     */
    void close();

    void put_particle( const char* rawParticleData );

  private:
    bool m_isClosed;

    frantic::particles::particle_array m_particles;

    // convert from m_channelMap to the m_particles channel map
    frantic::channels::channel_map_adaptor m_channelMapAdaptor;

    frantic::channels::channel_map m_channelMap;

    frantic::tstring m_filename;
};

} // namespace streams
} // namespace particles
} // namespace frantic
