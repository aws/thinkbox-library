// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 *
 *
 */

#pragma once

#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace particles {
namespace streams {

struct channel_range_map {
    double channelMinimum;
    double channelMaximum;
    double outputMinimum;
    double outputMaximum;
};

class channel_range_map_particle_istream : public frantic::particles::streams::delegated_particle_istream {
  private:
    frantic::tstring m_channelName;
    channel_range_map m_mapping;
    bool m_hasChannel;
    frantic::channels::channel_general_accessor m_channel;
    frantic::channels::channel_range_map_function_t m_rangeMap;

    void prepare_accessors() {
        using namespace frantic::particles;
        using namespace frantic::channels;

        const channel_map& pcm = m_delegate->get_channel_map();

        m_hasChannel = pcm.has_channel( m_channelName );
        if( m_hasChannel ) {
            m_channel = pcm.get_general_accessor( m_channelName );
            m_rangeMap = channel_range_map_function( m_channel.data_type() );
        }
    }

  public:
    channel_range_map_particle_istream( const frantic::particles::particle_istream_ptr& istream,
                                        const frantic::tstring& channel, const channel_range_map& mapping )
        : delegated_particle_istream( istream )
        , m_channelName( channel )
        , m_mapping( mapping ) {
        prepare_accessors();
    }

    virtual ~channel_range_map_particle_istream() {}

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        m_delegate->set_channel_map( particleChannelMap );
        prepare_accessors();
    }

    bool get_particle( char* rawParticleBuffer ) {
        bool result = m_delegate->get_particle( rawParticleBuffer );
        if( result ) {
            if( m_hasChannel ) {
                char* data = m_channel.get_channel_data_pointer( rawParticleBuffer );
                m_rangeMap( m_mapping.channelMinimum, m_mapping.channelMaximum, m_mapping.outputMinimum,
                            m_mapping.outputMaximum, data, m_channel.arity(), data );
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
