// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/channel_operation_compiler.hpp>
#include <frantic/particles/particle_array_iterator.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace frantic {
namespace particles {
namespace streams {

namespace detail {
template <class Value>
class indexed_blocked_range : public tbb::blocked_range<Value> {
    std::size_t m_startIndex;

  public:
    indexed_blocked_range()
        : m_startIndex( 0 )
        , tbb::blocked_range<Value>() {}

    indexed_blocked_range( Value begin, Value end, std::size_t grainsize, int startIndex )
        : m_startIndex( startIndex )
        , tbb::blocked_range<Value>( begin, end, grainsize ) {}

    indexed_blocked_range( indexed_blocked_range& r, tbb::split s )
        : tbb::blocked_range<Value>( r, s ) {
        m_startIndex = ( r.m_startIndex + r.size() );
    }

    std::size_t get_start_index() const { return m_startIndex; }
};

class parallel_eval {
    const frantic::channels::channel_operation_compiler* m_compData;

  public:
    typedef indexed_blocked_range<frantic::particles::particle_array_iterator> range_type;

  public:
    parallel_eval( const frantic::channels::channel_operation_compiler& compData )
        : m_compData( &compData ) {}
    void operator()( const range_type& range ) const {
        std::size_t index = range.get_start_index();
        range_type::const_iterator it = range.begin();
        range_type::const_iterator itEnd = range.end();

        for( ; it != itEnd; ++it, ++index )
            m_compData->eval( *it, index );
    }
};
} // end of namespace detail

class channel_operation_particle_istream : public delegated_particle_istream {
    typedef frantic::channels::channel_map channel_map;

    frantic::channels::channel_map m_outMap;
    frantic::channels::channel_map_adaptor m_adaptor;

    frantic::channels::channel_operation_compiler m_compiledData;

    boost::function<bool( int, const char* )> m_errorCallback;
    std::vector<frantic::channels::channel_op_node*> m_inputNodes; // Store these in case we need to re-compile.

    void init( const frantic::channels::channel_map& pcm ) {
        m_outMap = pcm;
        m_compiledData.reset( pcm, m_delegate->get_native_channel_map() );

        if( m_inputNodes.size() > 0 ) {
            try {
                m_inputNodes[0]->compile( m_inputNodes, m_compiledData );
            } catch( const frantic::channels::channel_compiler_error& e ) {
                // if( frantic::logging::is_logging_debug() )
                //	m_compiledData.debug_print(frantic::logging::debug);
                if( m_errorCallback ) {
                    bool showError = m_errorCallback( e.which_node(), e.what() );
                    if( !showError )
                        throw std::runtime_error( "" ); // Use an empty exception so it won't show up in the log window.
                }
                throw;
            }
        }

        m_delegate->set_channel_map( m_compiledData.get_channel_map() );
        m_adaptor.set( m_outMap, m_compiledData.get_channel_map() );
    }

  public:
    channel_operation_particle_istream( boost::shared_ptr<particle_istream> delegatePin,
                                        const std::vector<frantic::channels::channel_op_node*>& inNodes,
                                        boost::function<bool( int, const char* )> errorCallback )
        : delegated_particle_istream( delegatePin )
        , m_inputNodes( inNodes )
        , m_errorCallback( errorCallback ) {
        if( m_inputNodes.size() == 0 )
            throw std::runtime_error( "channel_operation_particle_istream() - The expression tree was empty" );
        init( m_delegate->get_channel_map() );
    }

    virtual ~channel_operation_particle_istream() {
        for( std::size_t i = 0; i < m_inputNodes.size(); ++i )
            delete m_inputNodes[i];
    }

    void set_channel_map( const frantic::channels::channel_map& pcm ) { init( pcm ); }

    std::size_t particle_size() const { return m_outMap.structure_size(); }

    const frantic::channels::channel_map& get_channel_map() const { return m_outMap; }

    const frantic::channels::channel_map& get_native_channel_map() const {
        return m_compiledData.get_native_channel_map();
    }

    bool get_particle( char* outBuffer ) {
        if( m_adaptor.is_identity() ) {
            if( !m_delegate->get_particle( outBuffer ) )
                return false;
            m_compiledData.eval( outBuffer, (int)m_delegate->particle_index() );
        } else {
            char* internalBuffer = (char*)_alloca( m_adaptor.source_size() );
            if( !m_delegate->get_particle( internalBuffer ) )
                return false;
            m_compiledData.eval( internalBuffer, (int)m_delegate->particle_index() );
            m_adaptor.copy_structure( outBuffer, internalBuffer );
        }

        return true;
    }

    bool get_particles( char* outBuffer, std::size_t& numParticles ) {
        bool result;
        int startIndex = (int)m_delegate->particle_index() + 1;

        if( m_adaptor.is_identity() ) {
            result = m_delegate->get_particles( outBuffer, numParticles );
            frantic::particles::particle_array_iterator itStart( outBuffer, m_adaptor.source_size() );
            frantic::particles::particle_array_iterator itEnd = itStart + numParticles;

#ifndef FRANTIC_DISABLE_THREADS
            tbb::parallel_for( detail::parallel_eval::range_type( itStart, itEnd, 2000, startIndex ),
                               detail::parallel_eval( m_compiledData ) );
#else
#pragma message( "Threads are disabled" )
            detail::parallel_eval func( m_compiledData );
            func( detail::parallel_eval::range_type( itStart, itEnd, 2000, startIndex ) );
#endif
        } else {
            boost::scoped_array<char> tempBuffer( new char[m_adaptor.source_size() * numParticles] );
            result = m_delegate->get_particles( tempBuffer.get(), numParticles );
            frantic::particles::particle_array_iterator itStart( tempBuffer.get(), m_adaptor.source_size() );
            frantic::particles::particle_array_iterator itEnd = itStart + numParticles;

#ifndef FRANTIC_DISABLE_THREADS
            tbb::parallel_for( detail::parallel_eval::range_type( itStart, itEnd, 2000, startIndex ),
                               detail::parallel_eval( m_compiledData ), tbb::simple_partitioner() );
#else
#pragma message( "Threads are disabled" )
            detail::parallel_eval func( m_compiledData );
            func( detail::parallel_eval::range_type( itStart, itEnd, 2000, startIndex ) );
#endif
            for( frantic::particles::particle_array_iterator it = itStart; it != itEnd;
                 ++it, outBuffer += m_adaptor.dest_size() )
                m_adaptor.copy_structure( outBuffer, *it );
        }

        return result;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
