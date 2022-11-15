// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/streams/particle_istream.hpp>

#include <tbb/parallel_reduce.h>

namespace frantic {
namespace particles {
namespace streams {

/**	This is the base policy class for culling policies used by the culling_particle_istream.
 *	It provides the necessary supporting code to parallelize culling operations.
 *  The typical use of culling_particle_istream is to instantiate it with a culling_policy
 *	derived from culling_policy_base. Ex:
 *
 *	new culling_particle_istream<volume_culling_policy>(delegateStream, boost::make_tuple(... volume culling args
 *...)); new culling_particle_istream<surface_culling_policy>(delegateStream, boost::make_tuple(... surface culling args
 *...));
 *
 *  The template arguement is for the Curiously Recursive template pattern. Set it to the derived class.
 */
template <class T>
class culling_policy_base {
  protected:
    // Used for TBB parallel_reduce purposes.
    char* m_outBuffer;   // Buffer storing output particles
    char* m_writeBuffer; // Ptr to end of m_outBuffer

  public:
    culling_policy_base()
        : m_outBuffer( NULL )
        , m_writeBuffer( NULL ) {}

    void reset_tbb_storage() { m_outBuffer = m_writeBuffer = NULL; }

    std::size_t result_size() { return ( m_writeBuffer - m_outBuffer ); }

    void operator()( const particle_range& range ) {
        if( !m_outBuffer ) {
            m_writeBuffer = m_outBuffer = range.begin_ptr();
        }

        char* p = range.begin_ptr();
        char* pEnd = range.end_ptr();
        for( ; p != pEnd; p += range.structure_size() ) {
            if( !static_cast<T*>( this )->cull( p ) ) {
                if( p != m_writeBuffer )
                    memcpy( m_writeBuffer, p, range.structure_size() );
                m_writeBuffer += range.structure_size();
            }
        }
    }

    void join( const culling_policy_base& rhs ) {
        std::size_t size = rhs.m_writeBuffer - rhs.m_outBuffer;
        memmove( m_writeBuffer, rhs.m_outBuffer, size );
        m_writeBuffer += size;
    }
};

template <class CullingPolicy>
class culling_particle_istream : public delegated_particle_istream, private CullingPolicy {
    boost::int64_t m_actualIndex;

  public:
    template <class T>
    culling_particle_istream( boost::shared_ptr<particle_istream> delegateStream, const T& cullingArgs )
        : delegated_particle_istream( delegateStream )
        , CullingPolicy( cullingArgs, m_delegate->get_channel_map() ) {
        m_actualIndex = -1;
    }

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        delegated_particle_istream::set_channel_map( pcm );
        static_cast<CullingPolicy*>( this )->set_channel_map( pcm );
    }

    boost::int64_t particle_index() const { return m_actualIndex; }
    boost::int64_t particle_count() const { return -1; }
    boost::int64_t particle_count_left() const { return -1; }

    bool get_particle( char* buffer ) {
        do {
            if( !m_delegate->get_particle( buffer ) )
                return false;
        } while( static_cast<CullingPolicy*>( this )->cull( buffer ) );

        ++m_actualIndex;
        return true;
    }

    bool get_particles( char* writeBuffer, std::size_t& numParticles ) {
        std::size_t totalOut = 0, particleSize = particle_size();
        bool notEos;

        do {
            std::size_t totalLeft = numParticles - totalOut;

            notEos = m_delegate->get_particles( writeBuffer, totalLeft );

#ifndef FRANTIC_DISABLE_THREADS
            CullingPolicy::reset_tbb_storage();
            tbb::parallel_reduce( particle_range( 0, totalLeft, writeBuffer, particleSize, 2500 ),
                                  *static_cast<CullingPolicy*>( this ) );
#else
#pragma message( "Threads are disabled" )
            ( *static_cast<CullingPolicy*>( this ) )( particle_range( 0, totalLeft, writeBuffer, particleSize, 2500 ) );
#endif

            writeBuffer += CullingPolicy::result_size();
            totalOut += CullingPolicy::result_size() / particleSize;
        } while( notEos && 2 * totalOut < numParticles );

        numParticles = totalOut;
        m_actualIndex += totalOut;
        return notEos;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
