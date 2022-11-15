// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/particle_file_metadata.hpp>

namespace frantic {
namespace particles {

namespace {

void merge_property_maps_no_repeats( frantic::channels::property_map& mergedInto,
                                     const frantic::channels::property_map& mergedFrom ) {

    for( size_t i = 0; i < mergedInto.get_channel_map().channel_count(); ++i ) {
        const frantic::channels::channel& ch = mergedInto.get_channel_map()[i];
        if( mergedFrom.get_channel_map().has_channel( ch.name() ) ) {
            throw std::runtime_error( "merge_property_maps_no_repeats -- a property with the name " +
                                      frantic::strings::to_string( ch.name() ) + " was already defined." );
        }
    }

    mergedInto.merge_property_map( mergedFrom );
}

} // namespace

particle_file_metadata::particle_file_metadata() {}

particle_file_metadata::particle_file_metadata( const frantic::channels::property_map& generalMetadata )
    : m_generalMetadata( generalMetadata ) {}

particle_file_metadata::particle_file_metadata( const frantic::channels::property_map& generalMetadata,
                                                const channel_metadata_t& channelMetadata )
    : m_generalMetadata( generalMetadata )
    , m_channelMetadata( channelMetadata ) {}

void particle_file_metadata::clear() {
    m_generalMetadata.clear();
    m_channelMetadata.clear();
}

void particle_file_metadata::get_channels_with_metadata( std::vector<frantic::tstring>& outNames ) const {
    outNames.clear();

    for( channel_metadata_t::const_iterator it = m_channelMetadata.begin(); it != m_channelMetadata.end(); ++it ) {
        outNames.push_back( it->first );
    }
}

const frantic::channels::property_map& particle_file_metadata::get_general_metadata() const {
    return m_generalMetadata;
}

frantic::channels::property_map& particle_file_metadata::get_general_metadata() { return m_generalMetadata; }

const frantic::channels::property_map*
particle_file_metadata::get_channel_metadata( const frantic::tstring& channelName ) const {
    channel_metadata_t::const_iterator it = m_channelMetadata.find( channelName );

    if( it != m_channelMetadata.end() ) {
        return &it->second;
    } else {
        return NULL;
    }
}

frantic::channels::property_map* particle_file_metadata::get_channel_metadata( const frantic::tstring& channelName,
                                                                               bool autoCreate ) {
    channel_metadata_t::iterator it = m_channelMetadata.find( channelName );

    if( it != m_channelMetadata.end() ) {
        return &it->second;
    } else if( autoCreate ) {
        return &m_channelMetadata[channelName];
    } else {
        return NULL;
    }
}

const particle_file_metadata::channel_metadata_t& particle_file_metadata::get_all_channel_metadata() const {
    return m_channelMetadata;
}

particle_file_metadata::channel_metadata_t& particle_file_metadata::get_all_channel_metadata() {
    return m_channelMetadata;
}

void particle_file_metadata::set_general_metadata( const frantic::channels::property_map& metadata ) {
    m_generalMetadata = metadata;
}

void particle_file_metadata::set_channel_metadata( const frantic::tstring& channelName,
                                                   const frantic::channels::property_map& metadata ) {
    if( metadata.empty() ) {
        channel_metadata_t::iterator it = m_channelMetadata.find( channelName );
        if( it != m_channelMetadata.end() ) {
            m_channelMetadata.erase( it );
        }
    } else {
        m_channelMetadata[channelName] = metadata;
    }
}

void particle_file_metadata::append( const particle_file_metadata& other, bool noRepeats ) {
    append_general_metadata( other.m_generalMetadata, noRepeats );

    for( channel_metadata_t::const_iterator it = other.m_channelMetadata.begin(); it != other.m_channelMetadata.end();
         ++it ) {
        append_channel_metadata( it->first, it->second, noRepeats );
    }
}

void particle_file_metadata::append_general_metadata( const frantic::channels::property_map& metadata,
                                                      bool noRepeats ) {
    if( !metadata.empty() ) {
        if( noRepeats ) {
            merge_property_maps_no_repeats( m_generalMetadata, metadata );
        } else {
            m_generalMetadata.merge_property_map( metadata );
        }
    }
}

void particle_file_metadata::append_channel_metadata( const frantic::tstring& channelName,
                                                      const frantic::channels::property_map& metadata,
                                                      bool noRepeats ) {
    if( !metadata.empty() ) {
        if( noRepeats ) {
            merge_property_maps_no_repeats( m_channelMetadata[channelName], metadata );
        } else {
            m_channelMetadata[channelName].merge_property_map( metadata );
        }
    }
}

bool particle_file_metadata::empty() const {
    if( !m_generalMetadata.empty() ) {
        return false;
    }

    if( !m_channelMetadata.empty() ) {
        for( channel_metadata_t::const_iterator it = m_channelMetadata.begin(); it != m_channelMetadata.end(); ++it ) {
            if( !it->second.empty() ) {
                return false;
            }
        }
    }

    return true;
}

} // namespace particles
} // namespace frantic
