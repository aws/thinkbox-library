// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <fstream>

#include <boost/scoped_array.hpp>

#include <frantic/files/zlib_reader.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/prt_file_header.hpp>
#include <frantic/particles/prt_metadata.hpp>

namespace frantic {
namespace particles {
namespace streams {

using frantic::channels::channel_map;
using frantic::channels::channel_map_adaptor;

namespace detail {

template <class ReadInterface>
class prt_particle_istream : public particle_istream {
    static const int BUFFER_SIZE = ( 32 << 10 ); // Use a 32Kb internal buffer when transforming the particle layout

    prt_file_header m_header;
    channel_map m_particleChannelMap;
    channel_map_adaptor m_pcmAdaptor;

    // These store the metadata in PRT2 format.  The original unaltered metadata are in m_header
    particle_file_metadata m_metadata;

    const frantic::tstring m_name;
    std::ios::off_type m_finOffset; // Offset to seek to when re-opening the file.
    bool m_finOpen;                 // Indicates whether the file has been opened for the first particle read yet.

    ReadInterface m_internalStream;

    boost::scoped_array<char> m_tempBuffer;
    boost::scoped_array<char> m_defaultParticle;
    boost::int64_t m_currentParticleIndex;

    void initialize_stream( const frantic::tstring& file ) {
        std::ifstream fin( file.c_str(), std::ios::binary );

        if( !fin )
            throw std::runtime_error( "prt_particle_istream: Failed to open file \"" +
                                      frantic::strings::to_string( file ) + "\" for reading." );

        // Load in the header
        m_header.read_header( fin, name() );

        // Allocate the temporary particle buffer
        m_tempBuffer.reset( new char[BUFFER_SIZE] );

        // Start the particle index one before the first valid index
        m_currentParticleIndex = -1;

        m_finOpen = false;
        m_finOffset = fin.tellg();

        // Setup the metadata so that they are in PRT2 format
        m_metadata = prt::convert_metadata_prt1_to_prt2( m_header.get_all_metadata() );
    }

    prt_particle_istream( const prt_particle_istream& rhs );
    prt_particle_istream& operator=( const prt_particle_istream& rhs );

  public:
    // This loads a particle file, and feeds back the results in the particle's native particle storage layout.
    prt_particle_istream( const frantic::tstring& file )
        : m_name( file ) {
        initialize_stream( file );

        // Now that the header is loaded, this function can be called to set the particle channel map
        set_channel_map( m_header.get_channel_map() );
    }

    // This loads a particle file, and feeds back the results in the requested particle storage layout.
    prt_particle_istream( const frantic::tstring& file, const channel_map& particleChannelMap )
        : m_name( file ) {
        initialize_stream( file );

        // Now that the header is loaded, this function can be called to set the particle channel map
        set_channel_map( particleChannelMap );
    }

    virtual ~prt_particle_istream() { close(); }

    /**
     * Get the general metadata field for the particle stream.
     * The metadata will be in PRT2 format
     */
    const frantic::channels::property_map& get_general_metadata() const { return m_metadata.get_general_metadata(); }

    /**
     * Get the given channel's metadata field for the particle stream
     * The metadata will be in PRT2 format
     */
    const frantic::channels::property_map* get_channel_metadata( const frantic::tstring& channelName ) const {
        return m_metadata.get_channel_metadata( channelName );
    }

    void close() {
        // Have to set 'm_finOpen = true' since it isn't actually opened until the first particle is read.
        m_finOpen = true;
        m_internalStream.close();
    }

    //////////////////////
    // Information functions
    //////////////////////

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_name; }

    boost::int64_t particle_count() const { return m_header.particle_count(); }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return m_header.particle_count() - m_currentParticleIndex - 1; }

    boost::int64_t particle_progress_count() const { return m_header.particle_count(); }

    boost::int64_t particle_progress_index() const { return m_currentParticleIndex; }

    const channel_map& get_channel_map() const { return m_particleChannelMap; }

    // Access to the channel_map in the file
    const channel_map& get_native_channel_map() const { return m_header.get_channel_map(); }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const channel_map& particleChannelMap ) {
        boost::scoped_array<char> newDefaultParticle( new char[particleChannelMap.structure_size()] );

        particleChannelMap.construct_structure( newDefaultParticle.get() );

        // If we previously had a default particle set, adapt its channels to the new layout.
        if( m_defaultParticle ) {
            frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
            defaultAdaptor.copy_structure( newDefaultParticle.get(), m_defaultParticle.get() );
        }

        m_defaultParticle.swap( newDefaultParticle );

        // Set the map and the adaptor
        m_particleChannelMap = particleChannelMap;
        m_pcmAdaptor.set( m_particleChannelMap, m_header.get_channel_map() );
    }

    void set_default_particle( char* rawParticleBuffer ) {
        m_particleChannelMap.copy_structure( m_defaultParticle.get(), rawParticleBuffer );
    }

    // This reads a particle into a buffer matching the channel_map.  It returns true if a particle was read, false
    // otherwise.
    // IMPORTANT: Make sure the buffer you pass in is at least as big as particle_size() bytes.
    bool get_particle( char* rawParticleBuffer ) {
        if( !m_finOpen ) {
            m_finOpen = true;
            m_internalStream.open( m_name );
            m_internalStream.seekg( m_finOffset );
        }

        if( m_currentParticleIndex + 1 < m_header.particle_count() ) {
            if( m_pcmAdaptor.is_identity() ) {
                // If we can, just read the particle data straight into the output buffer
                m_internalStream.read( rawParticleBuffer, m_particleChannelMap.structure_size() );
            } else {
                // Otherwise, read it into the temp buffer, then use the adaptor to copy it into the output buffer
                m_internalStream.read( m_tempBuffer.get(), m_header.get_channel_map().structure_size() );
                m_pcmAdaptor.copy_structure( rawParticleBuffer, m_tempBuffer.get(), m_defaultParticle.get() );
            }

            // Go to the next particle
            ++m_currentParticleIndex;

            return true;
        } else {
            return false;
        }
    }

    bool get_particles( char* particleBuffer, std::size_t& numParticles ) {
        if( !m_finOpen ) {
            m_finOpen = true;
            m_internalStream.open( m_name );
            m_internalStream.seekg( m_finOffset );
        }

        boost::int64_t particlesLeft = m_header.particle_count() - m_currentParticleIndex - 1;

        if( particlesLeft == 0 ) {
            numParticles = 0;
            return false;
        }

        bool result = true;

        if( (std::size_t)particlesLeft <= numParticles ) {
            numParticles = (std::size_t)particlesLeft;
            result = false;
        }

        m_currentParticleIndex += numParticles;

        if( m_pcmAdaptor.is_identity() ) {
            m_internalStream.read( particleBuffer, numParticles * m_particleChannelMap.structure_size() );
        } else {
            std::size_t numLeftToRead = numParticles;

            std::size_t numPerChunk = ( BUFFER_SIZE / m_pcmAdaptor.source_size() );
            std::size_t chunkSize = numPerChunk * m_pcmAdaptor.source_size();

            while( numLeftToRead > 0 ) {
                if( numLeftToRead < numPerChunk ) {
                    numPerChunk = numLeftToRead;
                    chunkSize = numPerChunk * m_pcmAdaptor.source_size();
                }

                char* curSrcParticle = m_tempBuffer.get();

                m_internalStream.read( curSrcParticle, chunkSize );

                for( std::size_t i = 0; i < numPerChunk;
                     ++i, particleBuffer += m_pcmAdaptor.dest_size(), curSrcParticle += m_pcmAdaptor.source_size() )
                    m_pcmAdaptor.copy_structure( particleBuffer, curSrcParticle, m_defaultParticle.get() );

                numLeftToRead -= numPerChunk;
            }
        }

        return result;
    }
};

} // End of namespace detail

// Replace the old definition with this typedef pointing to the new one.
typedef detail::prt_particle_istream<frantic::files::zlib_read_interface> prt_particle_istream;

} // namespace streams
} // namespace particles
} // namespace frantic
