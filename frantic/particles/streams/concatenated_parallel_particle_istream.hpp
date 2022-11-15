// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/threads/buffered_producer_consumer.hpp>

#include <tbb/spin_mutex.h>

#include <memory>

namespace frantic {
namespace particles {
namespace streams {

class concatenated_parallel_particle_istream : public particle_istream {
    friend class frantic::threads::buffered_producer_consumer<concatenated_parallel_particle_istream>;

    typedef frantic::particles::particle_array buffer_type;
    typedef boost::shared_ptr<frantic::particles::streams::particle_istream> particle_istream_ptr;

    tbb::spin_mutex m_pinMutex;
    std::vector<boost::shared_ptr<particle_istream>> m_delegates;
    std::size_t m_currentDelegate;

    /**
     * This class is instantiated in worker threads to handle the filling of particle buffers using the collection
     * of particle_istream objects.
     */
    class producer_instance : boost::noncopyable {
        concatenated_parallel_particle_istream& m_owner;
        particle_istream_ptr m_currentPin;

        std::size_t m_currentSize;

      public:
        producer_instance( concatenated_parallel_particle_istream& owner )
            : m_owner( owner )
            , m_currentSize( 0 ) {}

        bool can_produce_more() {
            if( !m_currentPin ) {
                tbb::spin_mutex::scoped_lock lock( m_owner.m_pinMutex );
                if( m_owner.m_currentDelegate >= m_owner.m_delegates.size() )
                    return false;

                m_currentPin = m_owner.m_delegates[m_owner.m_currentDelegate++];
            }

            return true;
        }

        bool is_buffer_empty( buffer_type* ) { return m_currentSize == 0; }

        bool is_buffer_full( buffer_type* ) { return m_currentSize > 50000; }

        void init_buffer( frantic::particles::particle_array* buffer ) {
            m_currentSize = 0;
            buffer->resize( 100000u );
        }

        void fill_buffer( frantic::particles::particle_array* buffer ) {
            boost::int64_t startProgress = m_currentPin->particle_progress_index();

            do {
                std::size_t numParticles = 100000u - m_currentSize;
                if( !m_currentPin->get_particles( buffer->at( m_currentSize ), numParticles ) ) {
                    m_owner.m_currentProgressParticle += ( m_currentPin->particle_progress_index() - startProgress );
                    m_currentPin.reset();
                }

                m_currentSize += numParticles;
            } while( m_currentSize < 50000u && m_currentPin );

            if( m_currentPin )
                m_owner.m_currentProgressParticle += ( m_currentPin->particle_progress_index() - startProgress );
        }

        void finish_buffer( frantic::particles::particle_array* buffer ) { buffer->resize( m_currentSize ); }
    };

    frantic::particles::particle_array* create_buffer() {
        return new frantic::particles::particle_array( m_channelMap );
    }

  private:
    boost::int64_t m_totalParticles;
    boost::int64_t m_currentParticle;
    boost::int64_t m_totalProgressParticles;

    tbb::atomic<boost::int64_t> m_currentProgressParticle;

    frantic::channels::channel_map m_nativeMap, m_channelMap;

    frantic::threads::buffered_producer_consumer<concatenated_parallel_particle_istream> m_pcImpl;

    std::unique_ptr<buffer_type> m_outputBuffer;
    std::size_t m_outBufferLeft;

    enum states { kInitNeeded, kReady };

    int m_state;

    void init_streams() {
        m_state = kInitNeeded;
        m_outBufferLeft = 0;

        m_currentDelegate = 0;
        m_currentParticle = -1;
        m_currentProgressParticle = -1;

        m_totalParticles = 0;
        for( std::size_t i = 0; i < m_delegates.size(); ++i ) {
            boost::int64_t count = m_delegates[i]->particle_count();
            if( count < 0 ) {
                m_totalParticles = -1;
                break;
            }
            m_totalParticles += count;
        }

        m_totalProgressParticles = 0; // sum of GetProgressCount() from all sub streams
        for( std::size_t i = 0; i < m_delegates.size(); ++i ) {
            boost::int64_t count = m_delegates[i]->particle_progress_count();
            if( count < 0 ) {
                m_totalProgressParticles = -1;
                break;
            }
            m_totalProgressParticles += count;
        }

        channel_map internalMap;
        m_nativeMap.reset();
        for( std::size_t i = 0; i < m_delegates.size(); ++i ) {
            internalMap.union_channel_map( m_delegates[i]->get_channel_map() );
            m_nativeMap.union_channel_map( m_delegates[i]->get_native_channel_map() );
        }
        m_nativeMap.end_channel_definition();
        internalMap.end_channel_definition( 4, true );

        set_channel_map( internalMap );
    }

  public:
    concatenated_parallel_particle_istream( const std::vector<boost::shared_ptr<particle_istream>>& pins )
        : m_delegates( pins ) {
        if( m_delegates.empty() )
            throw std::runtime_error(
                "concatenated_particle_istream() - The provided particle_istream array was empty.  It "
                "should contain at least one stream." );

        init_streams();
    }

    virtual ~concatenated_parallel_particle_istream() {}

    void close() {
        for( std::vector<boost::shared_ptr<particle_istream>>::iterator it = m_delegates.begin();
             it != m_delegates.end(); ++it )
            ( *it )->close();
    }

    void set_channel_map( const channel_map& particleChannelMap ) {
        m_channelMap = particleChannelMap;
        for( std::vector<boost::shared_ptr<particle_istream>>::iterator it = m_delegates.begin();
             it != m_delegates.end(); ++it )
            ( *it )->set_channel_map( particleChannelMap );
    }

    const channel_map& get_channel_map() const { return m_channelMap; }

    const channel_map& get_native_channel_map() const { return m_nativeMap; }

    void set_default_particle( char* rawParticleBuffer ) {
        for( std::vector<boost::shared_ptr<particle_istream>>::iterator it = m_delegates.begin();
             it != m_delegates.end(); ++it )
            ( *it )->set_default_particle( rawParticleBuffer );
    }

    bool get_particle( char* rawParticleBuffer ) {
        if( m_state == kInitNeeded ) {
            m_pcImpl.reset( *this, 0 );
            m_state = kReady;
        }

        if( !m_outputBuffer.get() ) {
            m_outputBuffer = m_pcImpl.steal_finished_item( true );
            if( !m_outputBuffer.get() )
                return false;

            m_outBufferLeft = m_outputBuffer->size();
        }

        m_channelMap.copy_structure( rawParticleBuffer, m_outputBuffer->at( --m_outBufferLeft ) );

        if( m_outBufferLeft == 0 )
            m_pcImpl.return_finished_item( std::move( m_outputBuffer ) );

        ++m_currentParticle;

        return true;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        if( m_state == kInitNeeded ) {
            m_pcImpl.reset( *this, 0 );
            m_state = kReady;
        }

        if( !m_outputBuffer.get() ) {
            m_outputBuffer = m_pcImpl.steal_finished_item( true );
            if( !m_outputBuffer.get() ) {
                numParticles = 0;
                return false;
            }

            m_outBufferLeft = m_outputBuffer->size();
        }

        numParticles = std::min( numParticles, m_outBufferLeft );

        for( std::size_t i = 0; i < numParticles; ++i, buffer += m_channelMap.structure_size() )
            m_channelMap.copy_structure( buffer, m_outputBuffer->at( --m_outBufferLeft ) );

        if( m_outBufferLeft == 0 )
            m_pcImpl.return_finished_item( std::move( m_outputBuffer ) );

        m_currentParticle += numParticles;

        return true;
    }

    frantic::tstring name() const {
        std::basic_stringstream<frantic::tchar> ss;
        ss << _T("Concatenated streams: ");
        for( std::size_t i = 0; i < m_delegates.size(); ++i )
            ss << m_delegates[i]->name() << _T(", ");
        return ss.str();
    }

    std::size_t particle_size() const { return m_delegates[0]->particle_size(); }

    boost::int64_t particle_count() const { return m_totalParticles; }

    boost::int64_t particle_count_guess() const {
        boost::int64_t result = 0;
        for( std::size_t i = 0; i < m_delegates.size(); ++i )
            result += m_delegates[i]->particle_count_guess();
        return result;
    }

    boost::int64_t particle_index() const { return m_currentParticle; }

    boost::int64_t particle_count_left() const {
        return ( m_totalParticles < 0 ) ? -1 : m_totalParticles - m_currentParticle - 1;
    }

    boost::int64_t particle_progress_count() const { return m_totalProgressParticles; }

    boost::int64_t particle_progress_index() const { return m_currentProgressParticle; }
};

} // namespace streams
} // namespace particles
} // namespace frantic
