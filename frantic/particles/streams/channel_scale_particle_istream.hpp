// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace particles {
namespace streams {

//	A particle stream that scales the density by a provided constant.
class channel_scale_particle_istream : public delegated_particle_istream {
  private:
    float m_channelScale;

    frantic::tstring m_channelName;
    bool m_hasChannel;
    frantic::channels::channel_general_accessor m_channel;

    void prepare_accessors() {
        const frantic::channels::channel_map& pcm = m_delegate->get_channel_map();

        m_hasChannel = pcm.has_channel( m_channelName );
        if( m_hasChannel )
            m_channel = pcm.get_general_accessor( m_channelName );
    }

  public:
    channel_scale_particle_istream( const boost::shared_ptr<particle_istream>& istream,
                                    const frantic::tstring& channelName, float channelScale )
        : delegated_particle_istream( istream )
        , m_channelName( channelName )
        , m_channelScale( channelScale ) {
        prepare_accessors();
    }

    virtual ~channel_scale_particle_istream() {}

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        m_delegate->set_channel_map( particleChannelMap );
        prepare_accessors();
    }

    bool get_particle( char* rawParticleBuffer ) {
        bool result = m_delegate->get_particle( rawParticleBuffer );
        if( result ) {
            if( m_hasChannel ) {

                char* data = m_channel.get_channel_data_pointer( rawParticleBuffer );
                m_channel.weighted_sum( &m_channelScale, &data, 1, data );
            }
        }
        return result;
    }

    bool get_particles( char* buffer, std::size_t& outSize ) {
        std::size_t particleSize = m_delegate->get_channel_map().structure_size();
        for( std::size_t i = 0; i < outSize; ++i, buffer += particleSize ) {
            if( !get_particle( buffer ) ) {
                outSize = i;
                return false;
            }
        }

        return true;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
