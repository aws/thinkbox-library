// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/make_shared.hpp>

#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

class concatenated_particle_istream : public particle_istream {
  private:
    std::vector<boost::shared_ptr<particle_istream>> m_delegates;
    std::size_t m_currentDelegate;

    boost::int64_t m_totalParticles;
    boost::int64_t m_currentParticle;
    boost::int64_t m_currentProgressParticle;
    boost::int64_t m_totalProgressParticles;

    frantic::channels::channel_map m_nativeMap;

    void init_streams() {
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

        frantic::channels::channel_map internalMap;
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
    concatenated_particle_istream( boost::shared_ptr<particle_istream> pin1,
                                   boost::shared_ptr<particle_istream> pin2 ) {
        m_delegates.push_back( pin1 );
        m_delegates.push_back( pin2 );
        init_streams();
    }

    concatenated_particle_istream( const std::vector<boost::shared_ptr<particle_istream>>& pins );

    virtual ~concatenated_particle_istream() {}

    void close() {
        for( std::vector<boost::shared_ptr<particle_istream>>::iterator it = m_delegates.begin();
             it != m_delegates.end(); ++it )
            ( *it )->close();
    }

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        for( std::vector<boost::shared_ptr<particle_istream>>::iterator it = m_delegates.begin();
             it != m_delegates.end(); ++it )
            ( *it )->set_channel_map( particleChannelMap );
    }

    const frantic::channels::channel_map& get_channel_map() const { return m_delegates[0]->get_channel_map(); }

    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    void set_default_particle( char* rawParticleBuffer ) {
        for( std::vector<boost::shared_ptr<particle_istream>>::iterator it = m_delegates.begin();
             it != m_delegates.end(); ++it )
            ( *it )->set_default_particle( rawParticleBuffer );
    }

    bool get_particle( char* rawParticleBuffer ) {
        while( m_currentDelegate < m_delegates.size() ) {
            if( m_delegates[m_currentDelegate]->get_particle( rawParticleBuffer ) ) {
                ++m_currentParticle;
                return true;
            } else {
                m_currentProgressParticle += m_delegates[m_currentDelegate]->particle_progress_index();
                m_delegates[m_currentDelegate++]->close();
            }
        }

        return false;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        if( m_currentDelegate >= m_delegates.size() ) {
            numParticles = 0;
            return false;
        }

        std::size_t totalOut = 0;
        std::size_t particleSize = this->particle_size();

        do {
            std::size_t desiredOut = numParticles - totalOut;
            bool notEos = m_delegates[m_currentDelegate]->get_particles( buffer, desiredOut );

            totalOut += desiredOut;
            buffer += desiredOut * particleSize;

            if( !notEos ) {
                FF_LOG( debug ) << "Exhausted stream " << m_currentDelegate << " named "
                                << m_delegates[m_currentDelegate]->name() << std::endl;
                m_currentProgressParticle += m_delegates[m_currentDelegate]->particle_progress_index();
                m_delegates[m_currentDelegate++]->close();
            }
        } while( m_currentDelegate < m_delegates.size() && totalOut < numParticles );

        numParticles = totalOut;
        m_currentParticle += totalOut;
        return m_currentDelegate < m_delegates.size();
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

    boost::int64_t particle_progress_index() const {
        if( m_currentDelegate < m_delegates.size() )
            return m_currentProgressParticle + m_delegates[m_currentDelegate]->particle_progress_index();
        else
            return m_currentProgressParticle;
    }
};

inline particle_istream_ptr create_concatenated_particle_istream( const std::vector<particle_istream_ptr>& pins ) {
    if( pins.size() == 1 ) {
        return pins[0];
    } else {
        return boost::make_shared<concatenated_particle_istream>( pins );
    }
}

} // namespace streams
} // namespace particles
} // namespace frantic
