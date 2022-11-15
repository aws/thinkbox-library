// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/particles/particle_cursor.hpp>

namespace frantic {
namespace particles {

class const_proxy_particle_cursor;

class proxy_particle_cursor : public particle_cursor {
    std::vector<boost::uint32_t>& m_proxyID;
    std::size_t offset;

    int m_proxyCursor;

    friend class const_proxy_particle_cursor;

    proxy_particle_cursor& operator=( const proxy_particle_cursor& ) {}

  public:
    proxy_particle_cursor( char* begin, char* end, std::vector<boost::uint32_t>& proxyParticleIDs,
                           channels::channel_map& channelMap )
        : particle_cursor( begin, end, channelMap.structure_size() )
        , m_proxyID( proxyParticleIDs ) {
        offset = channelMap.channel_offset( _T("ID") ); // the offset for the id channel in the channel map

        m_proxyCursor = -1;
    }

    void reset() {
        m_current = m_begin - m_particleSize;
        m_proxyCursor = -1;
    }

    // Resest the current value of the cursor to point to an invalid particle in front of the actual particles.
    bool next_particle() {
        if( m_proxyID.size() == 0 ) {
            return particle_cursor::next_particle();
        }

        ++m_proxyCursor;

        // checks for bounds conditions
        if( m_proxyCursor >= (int)m_proxyID.size() ) {
            m_current = m_end;
            return false;
        }

        masked_uint32<PRT_ID_MASK> checkid = m_proxyID[m_proxyCursor]; // the id to search for

        return binary_array_search( checkid, m_current + m_particleSize, m_end, offset, m_particleSize, m_current );
    }

    bool prev_particle() {
        if( m_proxyID.size() == 0 ) {
            return particle_cursor::prev_particle();
        }

        --m_proxyCursor;

        if( m_proxyCursor < 0 ) {
            m_current = m_begin - m_particleSize;
            return false;
        }

        masked_uint32<PRT_ID_MASK> checkid = m_proxyID[m_proxyCursor];

        return binary_array_search( checkid, m_begin, m_current, offset, m_particleSize, m_current );
    }

    bool move_to( int particleNumber ) {
        m_proxyCursor = particleNumber;

        if( m_proxyID.size() == 0 ) {
            return particle_cursor::move_to( particleNumber );
        }

        if( m_proxyCursor < 0 ) {
            m_current = m_begin - m_particleSize;
            return false;
        }

        if( m_proxyCursor >= (int)m_proxyID.size() ) {
            m_current = m_end;
            return false;
        }

        masked_uint32<PRT_ID_MASK> checkid = m_proxyID[m_proxyCursor];

        return binary_array_search( checkid, m_begin, m_end, offset, m_particleSize, m_current );
    }
};

class const_proxy_particle_cursor : public const_particle_cursor {
    const std::vector<boost::uint32_t>& m_proxyID;
    std::size_t offset;
    int m_proxyCursor;

    const_proxy_particle_cursor& operator=( const const_proxy_particle_cursor& ) {}

  public:
    const_proxy_particle_cursor( const char* begin, const char* end,
                                 const std::vector<boost::uint32_t>& proxyParticleIDs,
                                 const channels::channel_map& channelMap )
        : const_particle_cursor( begin, end, channelMap.structure_size() )
        , m_proxyID( proxyParticleIDs ) {

        offset = channelMap.channel_offset( _T("ID") ); // the offset for the id channel in the channel map

        m_proxyCursor = -1;
    }

    const_proxy_particle_cursor( const proxy_particle_cursor& rhs )
        : const_particle_cursor( rhs.m_begin, rhs.m_end, rhs.m_particleSize )
        , m_proxyID( rhs.m_proxyID ) {
        offset = rhs.offset;
        m_proxyCursor = -1;
    }

    void reset() {
        m_current = m_begin - m_particleSize;
        m_proxyCursor = -1;
    }

    // Resest the current value of the cursor to point to an invalid particle in front of the actual particles.
    bool next_particle() {
        if( m_proxyID.size() == 0 ) {
            return const_particle_cursor::next_particle();
        }

        ++m_proxyCursor;

        // checks for bounds conditions
        if( m_proxyCursor >= (int)m_proxyID.size() ) {
            m_current = m_end;
            return false;
        }

        masked_uint32<PRT_ID_MASK> checkid = m_proxyID[m_proxyCursor]; // the id to search for

        return binary_array_search( checkid, m_current + m_particleSize, m_end, offset, m_particleSize, m_current );
    }

    bool prev_particle() {
        if( m_proxyID.size() == 0 ) {
            return const_particle_cursor::prev_particle();
        }

        --m_proxyCursor;

        if( m_proxyCursor < 0 ) {
            m_current = m_begin - m_particleSize;
            return false;
        }

        masked_uint32<PRT_ID_MASK> checkid = m_proxyID[m_proxyCursor];

        return binary_array_search( checkid, m_begin, m_current, offset, m_particleSize, m_current );
    }

    bool move_to( int particleNumber ) {
        if( m_proxyID.size() == 0 ) {
            return const_particle_cursor::move_to( particleNumber );
        }

        m_proxyCursor = particleNumber;

        if( m_proxyCursor < 0 ) {
            m_current = m_begin - m_particleSize;
            return false;
        }

        if( m_proxyCursor >= (int)m_proxyID.size() ) {
            m_current = m_end;
            return false;
        }

        masked_uint32<PRT_ID_MASK> checkid = m_proxyID[m_proxyCursor];

        return binary_array_search( checkid, m_begin, m_end, offset, m_particleSize, m_current );
    }
};

} // namespace particles
} // namespace frantic
