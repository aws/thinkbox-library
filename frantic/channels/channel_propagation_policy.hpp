// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/strings/tstring.hpp>
#include <set>
#include <string>
#include <vector>

namespace frantic {
namespace channels {

/**
 * A class which defines which channels are to propagate
 */
class channel_propagation_policy {
  private:
    std::set<frantic::tstring> m_channelNames;
    bool m_isIncludePolicy;

  public:
    /**
     * Default constructor that creates a policy to exclude the enumerated channels
     */
    channel_propagation_policy()
        : m_isIncludePolicy( false ) {}

    /**
     * Constructor that allows the user to specify if the policy is to include or exclude the enumerated channels
     *
     * @param  isIncludePolicy   true will only include channels listed, false will only include channels that are not
     * listed
     */
    channel_propagation_policy( bool isIncludePolicy )
        : m_isIncludePolicy( isIncludePolicy ) {}

    /**
     * Constructor that allows the user to specify if the policy is to include or exclude the enumerated channels
     * and a set of channel names to include or exclude
     *
     * @param  isIncludePolicy   true will only include channels listed, false will only include channels that are not
     * listed
     * @param  channelNames       a set of channel name strings to include/exclude
     */
    channel_propagation_policy( bool isIncludePolicy, std::set<frantic::tstring> channelNames )
        : m_channelNames( channelNames )
        , m_isIncludePolicy( isIncludePolicy ) {}

    /**
     * Adds a channel name to our include/exclude list
     *
     * @param  channelName  the name of the channel to add
     */
    void add_channel( const frantic::tstring& channelName ) { m_channelNames.insert( channelName ); }

    /**
     * Removes a channel name from our include/exclude list
     *
     * @param  channelName  the name of the channel to remove
     */
    void remove_channel( const frantic::tstring& channelName ) {
        std::set<frantic::tstring>::iterator iter = m_channelNames.find( channelName );
        if( iter != m_channelNames.end() )
            m_channelNames.erase( iter );
    }

    /**
     * Sets our list of channel name to be used as the include/exclude list
     *
     * @param  channelNames  the set of channel names to be used
     */
    void set_channels( const std::set<frantic::tstring>& channelNames ) { m_channelNames = channelNames; }

    /**
     * Sets our policy to only include the specified channels
     */
    void set_to_include_policy() { m_isIncludePolicy = true; }

    /**
     * Sets our policy to only exclude the specified channels
     */
    void set_to_exclude_policy() { m_isIncludePolicy = false; }

    /**
     * Query function to see if a channel should be included in propagation
     *
     * @param  channelName  the name of the channel to query
     */
    bool is_channel_included( const frantic::tstring& channelName ) const {
        std::set<frantic::tstring>::const_iterator iter = m_channelNames.find( channelName );
        if( iter != m_channelNames.end() )
            return m_isIncludePolicy;
        return !m_isIncludePolicy;
    }

    /**
     * Takes an input vector of channel names and removes the channels
     * that are not to be propagated based on the current policy.
     *
     * @param  channelNames  input vector of channel names. vector will be changed.
     */
    void filter_channel_vector( std::vector<frantic::tstring>& channelNames ) const {
        std::vector<frantic::tstring> outChannelNames;
        for( std::size_t i = 0; i < channelNames.size(); ++i )
            if( is_channel_included( channelNames[i] ) )
                outChannelNames.push_back( channelNames[i] );
        channelNames.swap( outChannelNames );
    }

    /**
     * Is it a include list?
     *
     * @return retures true if the policy is an include list
     */
    bool is_include_list() const { return m_isIncludePolicy; }

    /**
     * Is it a exclude list?
     *
     * @return retures true if the policy is an exclude list
     */
    bool is_exclude_list() const { return !m_isIncludePolicy; }
};

} // namespace channels
} // namespace frantic
