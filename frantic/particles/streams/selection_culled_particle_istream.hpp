// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

namespace detail {

inline boost::uint32_t hash_uint32( boost::uint32_t a ) {
    a = ( a ^ 61 ) ^ ( a >> 16 );
    a = a + ( a << 3 );
    a = a ^ ( a >> 4 );
    a = a * 0x27d4eb2d;
    a = a ^ ( a >> 15 );
    return a;
}
} // namespace detail

class selection_culled_particle_istream : public delegated_particle_istream {
    frantic::channels::channel_map m_outMap;
    frantic::channels::channel_map_adaptor m_adaptor;

    bool m_useSoftSelection;
    bool m_resetSelection;

    boost::int64_t m_particleIndex;
    frantic::channels::channel_cvt_accessor<float> m_selAccessor;
    frantic::channels::channel_cvt_accessor<int> m_idAccessor;

  private:
    bool cull_particle( char* particle, boost::int64_t particleIndex ) {
        float selWeight = m_selAccessor.get( particle );
        if( selWeight > 0 ) {
            if( selWeight >= 1 || !m_useSoftSelection )
                return true;
            boost::uint32_t particleID;
            if( m_idAccessor.is_valid() )
                particleID = static_cast<boost::uint32_t>( m_idAccessor.get( particle ) );
            else
                particleID = static_cast<boost::uint32_t>( particleIndex );

            // For soft selections we will use a mixing function on the particle ID and see
            // if it falls in the lower fraction as indicated by the selection weight.
            boost::uint32_t maxVal = std::numeric_limits<boost::uint32_t>::max();
            boost::uint32_t cullThresh = static_cast<boost::uint32_t>( double( selWeight ) * double( maxVal ) );
            if( detail::hash_uint32( particleID ) < cullThresh )
                return true;
        }
        return false;
    }

  public:
    selection_culled_particle_istream( boost::shared_ptr<particle_istream> pDelegate, bool useSoftSelection,
                                       bool resetSelection = false )
        : delegated_particle_istream( pDelegate )
        , m_particleIndex( -1 )
        , m_useSoftSelection( useSoftSelection )
        , m_resetSelection( resetSelection ) {
        set_channel_map( pDelegate->get_channel_map() );
    }

    virtual ~selection_culled_particle_istream() {}

    std::size_t particle_size() const { return m_outMap.structure_size(); }
    boost::int64_t particle_count() const { return -1; }
    boost::int64_t particle_index() const { return m_particleIndex; }
    boost::int64_t particle_count_left() const { return -1; }

    const frantic::channels::channel_map& get_channel_map() const { return m_outMap; }

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        const frantic::channels::channel_map& nativePcm = m_delegate->get_native_channel_map();

        m_outMap = particleChannelMap;

        m_idAccessor.reset();
        m_selAccessor.reset( 0.f );

        frantic::channels::channel_map modifiedPcm = particleChannelMap;

        // We want an ID channel, iff we are using a soft selection, and there is a
        // valid, filled, ID channel from the delegate stream. We check the the native
        // channel map first, since an upstream particle_istream may be requesting an ID
        // but it won't neccessarily be meaningful at this point unless someone downstream
        // is filling it.
        if( m_useSoftSelection && nativePcm.has_channel( _T("ID") ) ) {
            if( !modifiedPcm.has_channel( _T("ID") ) )
                modifiedPcm.append_channel<int>( _T("ID") );
            m_idAccessor = modifiedPcm.get_cvt_accessor<int>( _T("ID") );
        }

        if( nativePcm.has_channel( _T("Selection") ) ) {
            if( !modifiedPcm.has_channel( _T("Selection") ) )
                modifiedPcm.append_channel<float>( _T("Selection") );
            m_selAccessor = modifiedPcm.get_cvt_accessor<float>( _T("Selection") );
        }

        m_adaptor.set( m_outMap, modifiedPcm );
        m_delegate->set_channel_map( modifiedPcm );
    }

    void set_default_particle( char* rawParticleBuffer ) {
        if( m_adaptor.is_identity() )
            m_delegate->set_default_particle( rawParticleBuffer );
        else {
            // The assumption is that I have added some memory to the end of the structure, so this
            // zero pads the extra space.
            char* temp = (char*)_alloca( m_adaptor.source_size() );
            memcpy( temp, rawParticleBuffer, m_adaptor.dest_size() );
            memset( temp + m_adaptor.dest_size(), 0, m_adaptor.source_size() - m_adaptor.dest_size() );
            m_delegate->set_default_particle( temp );
        }
    }

    virtual bool get_particle( char* rawParticleBuffer ) {
        char* targ = m_adaptor.is_identity() ? rawParticleBuffer : (char*)_alloca( m_adaptor.source_size() );

        do {
            if( !m_delegate->get_particle( targ ) )
                return false;
        } while( cull_particle( targ, m_delegate->particle_index() ) );

        if( m_resetSelection )
            m_selAccessor.set( targ, 1.f );

        ++m_particleIndex;
        if( !m_adaptor.is_identity() )
            m_adaptor.copy_structure( rawParticleBuffer, targ );
        return true;
    }

    virtual bool get_particles( char* buffer, std::size_t& numParticles ) {
        bool result = false;
        std::size_t count = 0;
        boost::int64_t particleBaseIndex = m_delegate->particle_index();

        if( m_adaptor.is_identity() ) {
            char* pIn = buffer;

            result = m_delegate->get_particles( pIn, numParticles );
            for( std::size_t i = 0; i < numParticles; ++i, pIn += m_adaptor.source_size() ) {
                if( !cull_particle( pIn, ++particleBaseIndex ) ) {
                    if( m_resetSelection )
                        m_selAccessor.set( pIn, 1.f );

                    if( pIn != buffer )
                        memcpy( buffer, pIn, m_adaptor.dest_size() );
                    count += 1;
                    buffer += m_adaptor.dest_size();
                }
            }
        } else {
            // TODO: This is probably not good, but it would require a very large architectural
            //        change to get rid of.
            boost::scoped_array<char> temp( new char[m_adaptor.source_size() * numParticles] );
            char* pIn = temp.get();

            result = m_delegate->get_particles( pIn, numParticles );
            for( std::size_t i = 0; i < numParticles; ++i, pIn += m_adaptor.source_size() ) {
                if( !cull_particle( pIn, ++particleBaseIndex ) ) {
                    if( m_resetSelection )
                        m_selAccessor.set( pIn, 1.f );

                    m_adaptor.copy_structure( buffer, pIn );
                    count += 1;
                    buffer += m_adaptor.dest_size();
                }
            }
        }

        numParticles = count;
        m_particleIndex += count;
        return result;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
