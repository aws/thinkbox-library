// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once


// This particle input stream adaptor can be used to set the value of a single individual channel of arbitrary supported
// type.

#include <boost/call_traits.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <tbb/parallel_for.h>

namespace frantic {
namespace particles {
namespace streams {

template <typename T>
class set_channel_functor {
    frantic::channels::channel_cvt_accessor<T> m_accessor;
    T m_value;

  public:
    typedef typename boost::call_traits<T>::param_type param_type;

  public:
    set_channel_functor( const frantic::channels::channel_cvt_accessor<T>& accessor, param_type value )
        : m_accessor( accessor )
        , m_value( value ) {}

    void operator()( const particle_range& range ) const {
        for( char *p = range.begin_ptr(), *pEnd = range.end_ptr(); p != pEnd; p += range.structure_size() )
            m_accessor.set( p, m_value );
    }
};

template <typename T>
class set_channel_particle_istream : public delegated_particle_istream {
    frantic::channels::channel_cvt_accessor<T> m_accessor;
    frantic::channels::channel_map m_nativeMap; // May have to add to the native map;
    frantic::tstring m_channelName;
    T m_value;

    void init( const frantic::channels::channel_map& pcm ) {
        // I made this handle situations where the overriden channel isn't actually
        // requested by the user.
        m_accessor.reset( T() );
        if( pcm.has_channel( m_channelName ) )
            m_accessor = pcm.get_cvt_accessor<T>( m_channelName );
    }

  public:
    typedef typename boost::call_traits<T>::param_type param_type;

  public:
    set_channel_particle_istream( boost::shared_ptr<particle_istream> pin, const frantic::tstring& channelName,
                                  param_type channelValue )
        : delegated_particle_istream( pin )
        , m_channelName( channelName )
        , m_value( channelValue ) {
        m_nativeMap = m_delegate->get_native_channel_map();
        if( !m_nativeMap.has_channel( channelName ) )
            m_nativeMap.append_channel<T>( channelName );
        init( m_delegate->get_channel_map() );
    }

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_delegate->set_channel_map( pcm );
        init( pcm );
    }

    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    bool get_particle( char* p ) {
        if( !m_delegate->get_particle( p ) )
            return false;

        m_accessor.set( p, m_value );
        return true;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        bool notEos = m_delegate->get_particles( buffer, numParticles );

#ifndef FRANTIC_DISABLE_THREADS
        tbb::parallel_for(
            particle_range( 0, numParticles, buffer, m_delegate->get_channel_map().structure_size(), 1000 ),
            set_channel_functor<T>( m_accessor, m_value ) );
#else
#pragma message( "Threads are disabled" )
        set_channel_functor<T> f( m_accessor, m_value );
        f( particle_range( 0, numParticles, buffer, m_delegate->get_channel_map().structure_size(), 1000 ) );
#endif

        return notEos;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
