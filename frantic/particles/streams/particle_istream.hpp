// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// Input stream for particles
#pragma once

#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>

#include <half.h>

#include <frantic/channels/channel_map.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/exception_stream.hpp>

#include <tbb/blocked_range.h>
#include <tbb/tbb_stddef.h>

//*************************************************************************************************************************************************************
// This 64 bit number is used to uniquely define the particle_istream interface's version. Any changes made to
// particle_istream MUST be paired with
// a change to this number! Millions of lives hang in the balance.
//
// As an arbitrary decision, I suggest we increment this number by one each time
// particle_istream is changed. If you don't follow this rule, you may choose a number that has been used already.
//*************************************************************************************************************************************************************
#define PARTICLE_ISTREAM_INTERFACE_VERSION 0x0000000000000001

namespace frantic {
namespace particles {
namespace streams {

// Implements the TBB Range concept
class particle_range {
    std::size_t m_begin, m_end;
    std::size_t m_grainSize;
    std::size_t m_particleSize;
    char* m_buffer;

  public:
    particle_range( std::size_t begin, std::size_t end, char* buffer, std::size_t particleSize,
                    std::size_t grainSize = 1 )
        : m_begin( begin )
        , m_end( end )
        , m_grainSize( grainSize )
        , m_particleSize( particleSize )
        , m_buffer( buffer ) {}

    particle_range( particle_range& lhs, tbb::split )
        : m_end( lhs.m_end )
        , m_grainSize( lhs.m_grainSize )
        , m_particleSize( lhs.m_particleSize )
        , m_buffer( lhs.m_buffer ) {
        m_begin = lhs.m_end = lhs.m_begin + ( lhs.m_end - lhs.m_begin ) / 2u;
    }

    bool empty() const { return m_begin >= m_end; }
    bool is_divisible() const { return m_end - m_begin > m_grainSize; }

    std::size_t begin_index() const { return m_begin; }
    std::size_t end_index() const { return m_end; }
    char* begin_ptr() const { return m_buffer + m_particleSize * m_begin; }
    char* end_ptr() const { return m_buffer + m_particleSize * m_end; }

    std::size_t structure_size() const { return m_particleSize; }
};

//////////////////////
// PARTICLE INPUT STREAM BASE CLASS
//////////////////////
class particle_istream {
  public:
    particle_istream() {}

    // Virtual destructor so that we can use allocated pointers (generally with boost::shared_ptr)
    virtual ~particle_istream() {}

    virtual void close() = 0;

    // The stream can return its filename or other identifier for better error messages.
    virtual frantic::tstring name() const = 0;

    // TODO: We should add a verbose_name function, which all wrapping streams are required to mark up in some way

    // This is the size of the particle structure which will be loaded, in bytes.
    virtual std::size_t particle_size() const = 0;

    // Returns the number of particles, or -1 if unknown
    virtual boost::int64_t particle_count() const = 0;

    // Returns the index of the last particle that was read.
    virtual boost::int64_t particle_index() const = 0;

    // Returns the number of particles left to be read, or -1 if unknown
    virtual boost::int64_t particle_count_left() const = 0;

    // Returns the total number of particles. In cases like that of csv, it approximates the number of particles
    // by using the file size and the position in the file. This is done in order to provide some progress feedback.
    virtual boost::int64_t particle_progress_count() const = 0;

    // Returns the index of the last particle that was read. In cases like that of csv, it approximates the index of the
    // particle by using the position in the file. This is done in order to provide some progress feedback.
    virtual boost::int64_t particle_progress_index() const = 0;

    // If a stream does not know how many particles it has, it can optionally override this function
    // to produce a guess of how many there will be. This guess will be used to pre-allocate storage
    // for this many particles if the user is concerned about memory performance.
    virtual boost::int64_t particle_count_guess() const { return particle_count(); }

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    virtual void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) = 0;

    // This is the particle channel map which specifies the byte layout of the particle structure that is being used.
    virtual const frantic::channels::channel_map& get_channel_map() const = 0;

    // This is the particle channel map which specifies the byte layout of the input to this stream.
    // NOTE: This value is allowed to change after the following conditions:
    //    * set_channel_map() is called (for example, the empty_particle_istream equates
    //      the native map with the external map)
    virtual const frantic::channels::channel_map& get_native_channel_map() const = 0;

    /**
     * This provides a default particle which should be used to fill in channels of the requested channel map
     * which are not supplied by the native channel map.
     * IMPORTANT: Make sure the buffer you pass in is at least as big as particle_size() bytes.
     */
    virtual void set_default_particle( char* rawParticleBuffer ) = 0;

    // This reads a particle into a buffer matching the channel_map.
    // It returns true if a particle was read, false otherwise.
    // IMPORTANT: Make sure the buffer you pass in is at least as big as particle_size() bytes.
    virtual bool get_particle( char* rawParticleBuffer ) = 0;

    // This reads a group of particles. Returns false if the end of the source
    // was reached during the read.
    // IMPORTANT: If the number of particles read is smaller than the provided `inoutNumParticles`,
    //            that value must be modified to reflect the actual number of particles read.
    virtual bool get_particles( char* rawParticleBuffer, std::size_t& inoutNumParticles ) = 0;

    // Convenience function for std::vector
    bool get_particle( std::vector<char>& rawParticleBuffer ) {
        // TODO: In a debug mode, we could confirm that the vector is big enough to hold the particle.
        return get_particle( &rawParticleBuffer[0] );
    }

    void set_default_particle( std::vector<char>& rawParticleBuffer ) {
        // TODO: In a debug mode, we could confirm that the vector is big enough to hold the particle.
        set_default_particle( &rawParticleBuffer[0] );
    }
};

// class delegated_particle_istream
//
// This class is a simple wrapper that delegates common functions to a member particle_istream
// The intended usage is for "post-processing" streams that modify data members without modifying the structure or
// number of particles in the stream. An example is transformed_particle_istream.
class delegated_particle_istream : public particle_istream {
  protected:
    boost::shared_ptr<particle_istream> m_delegate;

  public:
    delegated_particle_istream( boost::shared_ptr<particle_istream> pin )
        : m_delegate( pin ) {}
    virtual ~delegated_particle_istream() {}

    void close() { m_delegate->close(); }
    std::size_t particle_size() const { return get_channel_map().structure_size(); }
    boost::int64_t particle_count() const { return m_delegate->particle_count(); }
    boost::int64_t particle_index() const { return m_delegate->particle_index(); }
    boost::int64_t particle_count_left() const { return m_delegate->particle_count_left(); }
    boost::int64_t particle_progress_count() const { return m_delegate->particle_progress_count(); }
    boost::int64_t particle_progress_index() const { return m_delegate->particle_progress_index(); }
    boost::int64_t particle_count_guess() const { return m_delegate->particle_count_guess(); }
    frantic::tstring name() const { return m_delegate->name(); }

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        m_delegate->set_channel_map( particleChannelMap );
    }
    void set_default_particle( char* rawParticleBuffer ) { m_delegate->set_default_particle( rawParticleBuffer ); }
    const frantic::channels::channel_map& get_channel_map() const { return m_delegate->get_channel_map(); }
    const frantic::channels::channel_map& get_native_channel_map() const {
        return m_delegate->get_native_channel_map();
    }

    // These functions must still be implemented
    /*
     virtual bool get_particle( char* rawParticleBuffer ) = 0;
     virtual bool get_particles( char* buffer, std::size_t& numParticles ) = 0;
    */
};

typedef boost::shared_ptr<particle_istream> particle_istream_ptr;

} // namespace streams

using streams::particle_istream_ptr;
} // namespace particles
} // namespace frantic
