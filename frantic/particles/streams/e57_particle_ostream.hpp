// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#if defined( E57_AVAILABLE )

#include <E57Format/E57Format.h>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/particle_ostream.hpp>
#include <frantic/strings/tstring.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/thread/mutex.hpp>

namespace frantic {
namespace particles {
namespace streams {

/**
 * E57 File Format Writer
 *
 * @see http://www.libe57.org/
 * @see http://www.ri.cmu.edu/publication_view.html?pub_id=6767
 */
class e57_particle_ostream : public particle_ostream {
    /// Maximum size of the buffer before flushing
    static const size_t m_ACCEPTABLE_BUFFER_SIZE = 262144;
    /// Channel layout.
    channel_map m_particleChannelMap;
    /// Channel layout for file.
    channel_map m_particleChannelMapForFile;
    /// Particle channel map adaptor.
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    /// The file being written to.
    const frantic::tstring m_file;
    /// Number of particles processed so far.
    boost::int64_t m_currentParticleIndex;
    /// Number of particles to stream. If this is -1, then we don't know.
    boost::int64_t m_expectedParticleCount;
    /// Buffer for storing particles (fixed size).
    frantic::particles::particle_array m_particleBuffer;
    /// Number of particles in m_particleBuffer (note that m_particleBuffer.size() is NOT accurate).
    size_t m_particleBufferSize;
    /// File being written to.
    boost::scoped_ptr<e57::ImageFile> m_imf;
    /// Pointer to writer for writing particles. Owner pointer. Should be non-NULL always by end of constructor.
    boost::scoped_ptr<e57::CompressedVectorWriter> m_writer;

    // Private copy constructor and assignment operator to disable copying.
    e57_particle_ostream( const e57_particle_ostream& );            // not implemented
    e57_particle_ostream& operator=( const e57_particle_ostream& ); // not implemented

    /**
     * Initialize buffers and set up file-level metadata.
     */
    void initialize_stream( const frantic::channels::property_map* generalMetadata,
                            const std::map<frantic::tstring, frantic::channels::property_map>* channelMetadata );

    /**
     * Flush the particle buffers and schedule them for writing.
     */
    void flush();

  public:
    /**
     * Create an `e57_particle_ostream`.
     * @param file  The file to write to
     * @param particleChannelMap  Layout of the particle channel
     * @param expectedParticleCount  Number of particles to write
     * @param generalMetadata  Global metadata
     * @param channelMetadata  Channel metadata
     */
    e57_particle_ostream( const frantic::tstring& file, const channel_map& particleChannelMap,
                          boost::int64_t expectedParticleCount,
                          const frantic::channels::property_map* generalMetadata = NULL,
                          const std::map<frantic::tstring, frantic::channels::property_map>* channelMetadata = NULL );

    /** Destructor */
    virtual ~e57_particle_ostream();

    /**
     * Get particle channel map which specifies the byte layout of the particle structure.
     */
    const channel_map& get_channel_map() const;

    /**
     * Change the particle byte layout that's being saved on the fly.
     */
    void set_channel_map( const channel_map& particleChannelMap );

    /**
     * Close the file.
     */
    void close();

    /**
     * Return how big one particle is.
     */
    std::size_t particle_size() const;

    /**
     * Put a particle into the file.
     */
    void put_particle( const char* rawParticleData );
};
} // namespace streams
} // namespace particles
} // namespace frantic

#endif
