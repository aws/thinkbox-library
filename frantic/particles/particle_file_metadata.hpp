// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/property_map.hpp>
#include <frantic/strings/tstring.hpp>

#include <map>

namespace frantic {
namespace particles {

class particle_file_metadata {

  private:
    frantic::channels::property_map m_generalMetadata;
    typedef std::map<frantic::tstring, frantic::channels::property_map> channel_metadata_t;
    channel_metadata_t m_channelMetadata;

  public:
    particle_file_metadata();
    particle_file_metadata( const frantic::channels::property_map& generalMetadata );
    particle_file_metadata( const frantic::channels::property_map& generalMetadata,
                            const channel_metadata_t& channelMetadata );

    void clear();

    bool empty() const;

    void get_channels_with_metadata( std::vector<frantic::tstring>& outNames ) const;

    const frantic::channels::property_map& get_general_metadata() const;
    frantic::channels::property_map& get_general_metadata();
    const frantic::channels::property_map* get_channel_metadata( const frantic::tstring& channelName ) const;
    frantic::channels::property_map* get_channel_metadata( const frantic::tstring& channelName,
                                                           bool autoCreate = false );
    const channel_metadata_t& get_all_channel_metadata() const;
    channel_metadata_t& get_all_channel_metadata();

    void set_general_metadata( const frantic::channels::property_map& metadata );
    void set_channel_metadata( const frantic::tstring& channelName, const frantic::channels::property_map& metadata );

    void append( const particle_file_metadata& other, bool noRepeats = true );
    void append_general_metadata( const frantic::channels::property_map& metadata, bool noRepeats = true );
    void append_channel_metadata( const frantic::tstring& channelName, const frantic::channels::property_map& metadata,
                                  bool noRepeats = true );
};

} // namespace particles
} // namespace frantic
