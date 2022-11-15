// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <frantic/strings/tstring.hpp>

#include <boost/bind.hpp>

#ifndef FRANTIC_DISABLE_THREADS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#endif

namespace frantic {
namespace particles {
namespace streams {

/**
 * This class acts as a decorator that applies a functor to each particle moving through the stream
 * @tparam Functor A class acting as a function object. Must implement the following functions:
 *          void set_channel_map( channel_map& inoutChannelMap, const channel_map& nativeMap );
 *          OutType operator()( const char* particle ) const;
 *
 * @tparam OutType The type returned by the Functor object
 */
template <class Functor, class OutType>
class functor_particle_istream : public delegated_particle_istream {
  protected:
    frantic::channels::channel_map m_outMap, m_nativeMap;
    frantic::channels::channel_map_adaptor m_adaptor;

    Functor m_functor;
    frantic::tstring m_channelName;

    frantic::channels::channel_cvt_accessor<OutType> m_outAccessor;

  private:
    /**
     * Applies the functor to a range of particle data
     * @param dest Pointer to the raw particle data (described by 'm_outMap')
     * @param range The range of particles to process
     */
    void apply( char* dest, const tbb::blocked_range<std::size_t>& range ) const {
        char* it = dest + m_outMap.structure_size() * range.begin();
        char* itEnd = dest + m_outMap.structure_size() * range.end();

        for( ; it != itEnd; it += m_outMap.structure_size() )
            m_outAccessor.set( it, m_functor( it ) );
    }

    /**
     * Applies the functor to a range of particle data, and copies from on memory location to another using a
     * channel_map_adaptor
     * @param dest Pointer to the destination particle data (described by 'm_outMap')
     * @param src Pointer to the source particle data (described by 'm_delegate->get_channel_map()')
     * @param range The range of particles to process
     */
    void apply_with_adaptor( char* dest, const char* src, const tbb::blocked_range<std::size_t>& range ) const {
        const char* itSrc = src + m_adaptor.source_size() * range.begin();
        char* itDest = dest + m_adaptor.dest_size() * range.begin();
        char* itDestEnd = dest + m_adaptor.dest_size() * range.end();

        for( ; itDest != itDestEnd; itDest += m_adaptor.dest_size(), itSrc += m_adaptor.source_size() ) {
            m_adaptor.copy_structure( itDest, itSrc );
            m_outAccessor.set( itDest, m_functor( itSrc ) );
        }
    }

  public:
    /**
     * Constructor.
     * @param delegatePin The particle_istream to decorate.
     * @param channelName The name of the channel to write to.
     * @param fn The functor to apply to each particle. A copy of the argument will be made.
     */
    functor_particle_istream( particle_istream_ptr delegatePin, const frantic::tstring& channelName, const Functor& fn )
        : delegated_particle_istream( delegatePin )
        , m_functor( fn )
        , m_channelName( channelName ) {
        m_outMap = m_delegate->get_channel_map();
        m_nativeMap = m_delegate->get_native_channel_map();

        if( !m_nativeMap.has_channel( m_channelName ) )
            m_nativeMap.append_channel<OutType>( m_channelName );

        if( m_outMap.has_channel( m_channelName ) ) {
            m_outAccessor = m_outMap.get_cvt_accessor<OutType>( m_channelName );

            frantic::channels::channel_map delegateMap = m_outMap;
            m_functor.set_channel_map( delegateMap, m_nativeMap );

            m_adaptor.set( m_outMap, delegateMap );

            if( !m_adaptor.is_identity() )
                m_delegate->set_channel_map( delegateMap );
        } else {
            m_adaptor.set( m_outMap, m_outMap );
        }
    }

    virtual ~functor_particle_istream() {}

    virtual void set_channel_map( const frantic::channels::channel_map& pcm ) {
        if( m_outMap != pcm ) {
            m_outMap = pcm;

            if( m_outMap.has_channel( m_channelName ) ) {
                m_outAccessor = m_outMap.get_cvt_accessor<OutType>( m_channelName );

                frantic::channels::channel_map delegateMap = pcm;
                m_functor.set_channel_map( delegateMap, m_nativeMap );

                m_adaptor.set( m_outMap, delegateMap );

                m_delegate->set_channel_map( delegateMap );
            } else {
                m_outAccessor.reset();
                m_adaptor.set( m_outMap, m_outMap );
                m_delegate->set_channel_map( m_outMap );
            }
        }
    }

    virtual const frantic::channels::channel_map& get_channel_map() const { return m_outMap; }

    virtual const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    virtual void set_default_particle( char* rawParticleBuffer ) {
        if( !m_adaptor.is_identity() ) {
            const frantic::channels::channel_map& delegateMap = m_delegate->get_channel_map();

            frantic::channels::channel_map_adaptor tempAdaptor( delegateMap, m_outMap );

            char* tempDefault = (char*)alloca( delegateMap.structure_size() );
            delegateMap.construct_structure( tempDefault );
            tempAdaptor.copy_structure( tempDefault, rawParticleBuffer );

            m_delegate->set_default_particle( tempDefault );
        } else {
            m_delegate->set_default_particle( rawParticleBuffer );
        }
    }

    virtual bool get_particle( char* particle ) {
        if( m_adaptor.is_identity() ) {
            if( !m_delegate->get_particle( particle ) )
                return false;

            if( m_outAccessor.is_valid() )
                m_outAccessor.set( particle, m_functor( particle ) );
        } else {
            char* temp = (char*)alloca( m_adaptor.source_size() );
            if( !m_delegate->get_particle( temp ) )
                return false;

            m_adaptor.copy_structure( particle, temp );

            if( m_outAccessor.is_valid() )
                m_outAccessor.set( particle, m_functor( temp ) );
        }

        return true;
    }

    virtual bool get_particles( char* particles, std::size_t& numParticles ) {
        bool result;

        if( m_adaptor.is_identity() ) {
            result = m_delegate->get_particles( particles, numParticles );

            // Only apply the functor if the output channel is actually requested.
            if( m_outAccessor.is_valid() ) {
#ifndef FRANTIC_DISABLE_THREADS
                tbb::parallel_for( tbb::blocked_range<std::size_t>( 0, numParticles, 2000 ),
                                   boost::bind( &functor_particle_istream::apply, this, particles, _1 ) );
#else
                this->apply( particles, tbb::blocked_range<std::size_t>( 0, numParticles, 2000 ) );
#endif
            }
        } else {
            // TODO: This could be allocated once and reused across calls to get_particles().
            boost::scoped_array<char> temp( new char[numParticles * m_adaptor.source_size()] );

            result = m_delegate->get_particles( temp.get(), numParticles );

            if( m_outAccessor.is_valid() ) {
#ifndef FRANTIC_DISABLE_THREADS
                tbb::parallel_for(
                    tbb::blocked_range<std::size_t>( 0, numParticles, 2000 ),
                    boost::bind( &functor_particle_istream::apply_with_adaptor, this, particles, temp.get(), _1 ) );
#else
                this->apply_with_adaptor( particles, temp.get(),
                                          tbb::blocked_range<std::size_t>( 0, numParticles, 2000 ) );
#endif
            }
        }

        return result;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
