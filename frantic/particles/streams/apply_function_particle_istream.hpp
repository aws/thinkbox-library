// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/array.hpp>
#include <boost/function.hpp>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <frantic/channels/channel_accessor.hpp>
#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace frantic {
namespace particles {
namespace streams {

template <class T>
struct remove_const_ref {
    typedef typename boost::remove_const<typename boost::remove_reference<T>::type>::type value_type;
};

/**
 * This template and its relatives are responsible for taking an output accessor
 * and a tuple of input accessors, then hooking the accessors i/o into the supplied function
 */
template <class F>
struct ParticleDispatcher;
template <class R>
struct ParticleDispatcher<R( void )> {
    static const int NUM_ARGS = 0;
    typedef R result_type;
    typedef boost::function<R( void )> function_type;
    // It doesn't matter what type args_type is in this case
    typedef int args_type;

    static void dispatch( const function_type& f, const frantic::channels::channel_cvt_accessor<R>& dest,
                          args_type /*src*/, char* buffer ) {
        dest.set( buffer, f() );
    }
};

template <class R, class A1>
struct ParticleDispatcher<R( A1 )> {
    static const int NUM_ARGS = 1;
    typedef R result_type;
    typedef boost::function<R( A1 )> function_type;
    typedef boost::tuple<channels::channel_const_cvt_accessor<typename remove_const_ref<A1>::value_type>> args_type;

    static void dispatch( const function_type& f, const frantic::channels::channel_cvt_accessor<R>& dest,
                          const args_type& src, char* buffer ) {
        dest.set( buffer, f( src.template get<0>().get( buffer ) ) );
    }
};

template <class R, class A1, class A2>
struct ParticleDispatcher<R( A1, A2 )> {
    static const int NUM_ARGS = 2;
    typedef R result_type;
    typedef boost::function<R( A1, A2 )> function_type;
    typedef boost::tuple<channels::channel_const_cvt_accessor<typename remove_const_ref<A1>::value_type>,
                         channels::channel_const_cvt_accessor<typename remove_const_ref<A2>::value_type>>
        args_type;

    static void dispatch( const function_type& f, const frantic::channels::channel_cvt_accessor<R>& dest,
                          const args_type& src, char* buffer ) {
        dest.set( buffer, f( src.template get<0>().get( buffer ), src.template get<1>().get( buffer ) ) );
    }
};

template <class R, class A1, class A2, class A3>
struct ParticleDispatcher<R( A1, A2, A3 )> {
    static const int NUM_ARGS = 3;
    typedef R result_type;
    typedef boost::function<R( A1, A2, A3 )> function_type;
    typedef boost::tuple<channels::channel_const_cvt_accessor<typename remove_const_ref<A1>::value_type>,
                         channels::channel_const_cvt_accessor<typename remove_const_ref<A2>::value_type>,
                         channels::channel_const_cvt_accessor<typename remove_const_ref<A3>::value_type>>
        args_type;

    static void dispatch( const function_type& f, const frantic::channels::channel_cvt_accessor<R>& dest,
                          const args_type& src, char* buffer ) {
        dest.set( buffer, f( src.template get<0>().get( buffer ), src.template get<1>().get( buffer ),
                             src.template get<2>().get( buffer ) ) );
    }
};

template <class R, class A1, class A2, class A3, class A4>
struct ParticleDispatcher<R( A1, A2, A3, A4 )> {
    static const int NUM_ARGS = 4;
    typedef R result_type;
    typedef boost::function<R( A1, A2, A3, A4 )> function_type;
    typedef boost::tuple<channels::channel_const_cvt_accessor<typename remove_const_ref<A1>::value_type>,
                         channels::channel_const_cvt_accessor<typename remove_const_ref<A2>::value_type>,
                         channels::channel_const_cvt_accessor<typename remove_const_ref<A3>::value_type>,
                         channels::channel_const_cvt_accessor<typename remove_const_ref<A4>::value_type>>
        args_type;

    static void dispatch( const function_type& f, const frantic::channels::channel_cvt_accessor<R>& dest,
                          const args_type& src, char* buffer ) {
        dest.set( buffer, f( src.template get<0>().get( buffer ), src.template get<1>().get( buffer ),
                             src.template get<2>().get( buffer ), src.template get<3>().get( buffer ) ) );
    }
};

template <class F>
class apply_function_impl {
    typedef typename ParticleDispatcher<F>::result_type destination_type;
    typedef typename ParticleDispatcher<F>::function_type function_type;
    typedef typename ParticleDispatcher<F>::args_type args_type;
    enum { NUM_ARGS = ParticleDispatcher<F>::NUM_ARGS };

  private:
    frantic::channels::channel_cvt_accessor<destination_type> m_destAccessor;
    args_type m_srcAccessors;
    function_type m_function;
    const channels::channel_map_adaptor& m_adaptor;
    char* m_destBuffer;

    apply_function_impl& operator=( const apply_function_impl& ) {} // disabled assignment operator

  public:
    apply_function_impl( const frantic::channels::channel_cvt_accessor<destination_type>& dest, const args_type& srcs,
                         const function_type& func, const channels::channel_map_adaptor& adaptor, char* destBuffer )
        : m_destAccessor( dest )
        , m_srcAccessors( srcs )
        , m_function( func )
        , m_adaptor( adaptor )
        , m_destBuffer( destBuffer ) {}

    void operator()( const particle_range& range ) const {
        char* p = range.begin_ptr();
        char* pEnd = range.end_ptr();

        if( m_adaptor.is_identity() ) {
            for( ; p != pEnd; p += range.structure_size() )
                ParticleDispatcher<F>::dispatch( m_function, m_destAccessor, m_srcAccessors, p );
        } else {
            char* pDest = m_destBuffer + range.begin_index() * m_adaptor.dest_size();
            for( ; p != pEnd; p += range.structure_size(), pDest += m_adaptor.dest_size() ) {
                ParticleDispatcher<F>::dispatch( m_function, m_destAccessor, m_srcAccessors, p );
                m_adaptor.copy_structure( pDest, p );
            }
        }
    }
};

/**
 * This class allows an arbitrary function to be applied to every particle in a stream
 * by taking an output channel and boost::array of input channels. The supplied channels must
 * must match the signature of the function. It will also handle modifying its channel map to request
 * additional channels that are not in its output.
 */
template <class F>
class apply_function_particle_istream : public delegated_particle_istream {
  private:
    typedef typename ParticleDispatcher<F>::result_type destination_type;
    typedef typename ParticleDispatcher<F>::function_type function_type;
    typedef typename ParticleDispatcher<F>::args_type args_type;
    enum { NUM_ARGS = ParticleDispatcher<F>::NUM_ARGS };

  private:
    frantic::tstring m_destChannel;
    frantic::channels::channel_cvt_accessor<destination_type> m_destAccessor;

    args_type m_srcAccessors;
    function_type m_function;
    boost::array<frantic::tstring, NUM_ARGS> m_srcChannels;

    // When an argument references a channel not specified in the user-supplied channel
    // map, it will request a modified channel map from the delegate with the extra channels.
    channels::channel_map m_outMap, m_nativeMap;
    channels::channel_map_adaptor m_adaptor;
    boost::scoped_array<char> m_tempBuffer;

  public:
    apply_function_particle_istream( boost::shared_ptr<particle_istream> pin, const function_type& function,
                                     const frantic::tstring& destChannel,
                                     const boost::array<frantic::tstring, NUM_ARGS>& srcChannels )
        : delegated_particle_istream( pin )
        , m_function( function )
        , m_destChannel( destChannel )
        , m_srcChannels( srcChannels ) {
        set_channel_map( m_delegate->get_channel_map() );

        // Add the output of this operation to the native channel map if it isn't already there
        m_nativeMap = m_delegate->get_native_channel_map();
        if( !m_nativeMap.has_channel( m_destChannel ) )
            m_nativeMap.append_channel<destination_type>( m_destChannel );
    }

    virtual ~apply_function_particle_istream() {}

    inline void set_channel_map( const channels::channel_map& pcm );

    std::size_t particle_size() const { return m_outMap.structure_size(); }

    const channels::channel_map& get_channel_map() const { return m_outMap; }
    const channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    void set_default_particle( char* buffer ) {
        if( m_adaptor.is_identity() )
            m_delegate->set_default_particle( buffer );
        else {
            channels::channel_map_adaptor tempAdaptor( m_delegate->get_channel_map(), m_outMap );

            boost::scoped_ptr<char> pDefault( new char[tempAdaptor.dest_size()] );
            memset( pDefault.get(), 0, tempAdaptor.dest_size() );
            tempAdaptor.copy_structure( pDefault.get(), buffer );

            m_delegate->set_default_particle( pDefault.get() );
        }
    }

    bool get_particle( char* outBuffer ) {
        char* inBuffer = ( m_adaptor.is_identity() ) ? outBuffer : m_tempBuffer.get();

        if( !m_delegate->get_particle( inBuffer ) )
            return false;

        // It is imperative that the function be dispatched on the input buffer, because there are
        // no guarantees that the required channels exist in the output buffer
        ParticleDispatcher<F>::dispatch( m_function, m_destAccessor, m_srcAccessors, inBuffer );

        if( inBuffer != outBuffer )
            m_adaptor.copy_structure( outBuffer, inBuffer );

        return true;
    }

    bool get_particles( char* outBuffer, std::size_t& numParticles ) {
        boost::scoped_array<char> tempBuffer;
        char* target;

        if( !m_adaptor.is_identity() ) {
            tempBuffer.reset( new char[numParticles * m_adaptor.source_size()] );
            target = tempBuffer.get();
        } else
            target = outBuffer;

        bool notEos = m_delegate->get_particles( target, numParticles );

#ifndef FRANTIC_DISABLE_THREADS
        tbb::parallel_for( particle_range( 0, numParticles, target, m_adaptor.source_size(), 1000 ),
                           apply_function_impl<F>( m_destAccessor, m_srcAccessors, m_function, m_adaptor, outBuffer ) );
#else
#pragma message( "Threads are disabled" )
        apply_function_impl<F> f( m_destAccessor, m_srcAccessors, m_function, m_adaptor, outBuffer );
        f( particle_range( 0, numParticles, target, m_adaptor.source_size(), 1000 ) );
#endif
        return notEos;
    }
};

namespace detail {

/**
 * This class template will traverse the typelist of the arguments to the template function F
 * For each arg, it will get the relevant cvt_accessor from the supplied map and store it in
 * srcAccessors[Index]. The template allows for a compile-time loop over the arguments.
 */
template <class Tuple, int Index>
struct init_arg_accessors {
    inline static void apply( const channels::channel_map& pcm, const frantic::tstring* channelArray,
                              Tuple& srcAccessors ) {
        typedef typename boost::tuples::element<Index, Tuple>::type::value_type T;

        if( !pcm.has_channel( channelArray[Index] ) )
            srcAccessors.template get<Index>().reset( T() );
        else
            srcAccessors.template get<Index>() = pcm.get_const_cvt_accessor<T>( channelArray[Index] );
        init_arg_accessors<Tuple, Index - 1>::apply( pcm, channelArray, srcAccessors );
    }
};

template <class Tuple>
struct init_arg_accessors<Tuple, -1> {
    inline static void apply( const channels::channel_map&, const frantic::tstring*, Tuple& ) {}
};

/**
 * This class template will traverse the typelist of the arguments to the template function F
 * For each arg, it will ensure the supplied map has a channel with the given name, otherwise
 * it will create a new channel with the relevant type.
 */
template <class Tuple, int Index>
struct ensure_arg_channels {
    inline static void apply( channels::channel_map& pcm, const channels::channel_map& nativeMap,
                              const frantic::tstring* channelArray ) {
        typedef typename boost::tuples::element<Index, Tuple>::type::value_type T;

        if( !pcm.has_channel( channelArray[Index] ) && nativeMap.has_channel( channelArray[Index] ) )
            pcm.append_channel<T>( channelArray[Index] );
        ensure_arg_channels<Tuple, Index - 1>::apply( pcm, nativeMap, channelArray );
    }
};

template <class Tuple>
struct ensure_arg_channels<Tuple, -1> {
    inline static void apply( channels::channel_map&, const channels::channel_map&, const frantic::tstring* ) {}
};

} // namespace detail

template <class F>
void apply_function_particle_istream<F>::set_channel_map( const channels::channel_map& pcm ) {
    channels::channel_map requested = m_outMap = pcm;
    const channels::channel_map& nativeMap = m_delegate->get_native_channel_map();

    // Append any missing channels to the delegate's channel map
    if( !requested.has_channel( m_destChannel ) )
        requested.append_channel<destination_type>( m_destChannel );

    // Loop through all the channels and add them to the channel map if not already there.
    detail::ensure_arg_channels<args_type, boost::tuples::length<args_type>::value - 1>::apply( requested, nativeMap,
                                                                                                m_srcChannels.data() );

    // Initialize the accessors.
    m_destAccessor = requested.get_cvt_accessor<destination_type>( m_destChannel );

    // Loop through all channels and initialize the source accessors
    detail::init_arg_accessors<args_type, boost::tuples::length<args_type>::value - 1>::apply(
        requested, m_srcChannels.data(), m_srcAccessors );

    m_delegate->set_channel_map( requested );
    m_adaptor.set( m_outMap, requested );
    m_tempBuffer.reset( m_adaptor.is_identity() ? NULL : new char[requested.structure_size()] );
}

} // namespace streams
} // namespace particles
} // namespace frantic
