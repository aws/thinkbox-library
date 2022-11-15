// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/particles/particle_file_metadata.hpp>

namespace frantic {
namespace particles {

// This class represents a PRT file header, and is able to read or write it to an istream or ostream, respectively, on
// disk.
class prt_file_header {
    // This is the particle channel map of the file representation
    channels::channel_map m_particleChannelMap;

    boost::int64_t m_particleCount;

    bool m_allowNegativeCount;

    particle_file_metadata m_metadata;

    std::streampos m_particleCountPos; // Stores the stream position where the particle count must be written.

    std::streampos m_boundboxPos; // Stores the stream position where the boundbox value must be written.

  public:
    prt_file_header() {
        m_particleCount = -1;
        m_allowNegativeCount = false;
    }

    prt_file_header( const channels::channel_map& particleChannelMap )
        : m_particleChannelMap( particleChannelMap ) {
        if( particleChannelMap.needs_scope_management() )
            throw std::runtime_error(
                "prt_file_header() - The channel map particle layout requested includes non-raw data, "
                "such as a string.  The particle layout must be raw numeric data." );
        m_particleCount = -1;
        m_allowNegativeCount = false;
    }

    void set_allow_negative_count( bool allow ) { m_allowNegativeCount = allow; }

    /**
     * \return A property_map storing the general metadata from the prt file.
     */
    const frantic::channels::property_map& get_general_metadata() const { return m_metadata.get_general_metadata(); }

    /**
     * Use this function to assign metadata that will be written by the next call to write_header().
     * \param metadata The collection of names & values to write into the prt file.
     */
    void set_general_metadata( const frantic::channels::property_map& metadata );

    const frantic::channels::property_map* get_channel_metadata( const frantic::tstring& channelName ) const {
        return m_metadata.get_channel_metadata( channelName );
    }

    void set_channel_metadata( const frantic::tstring& channelName, const frantic::channels::property_map& metadata ) {
        if( !m_particleChannelMap.has_channel( channelName ) )
            throw std::invalid_argument(
                "prt_file_header.get_channel_metadata() Cannot add metadata to non-existing channel named \"" +
                frantic::strings::to_string( channelName ) + "\"" );

        m_metadata.set_channel_metadata( channelName, metadata );
    }

    const particle_file_metadata& get_all_metadata() const { return m_metadata; }

    void set_all_metadata( const particle_file_metadata& data );

    //////////////////////
    // Functions for dealing with the particle channel map
    //////////////////////

    void set_channel_map( const channels::channel_map& particleChannelMap ) {
        if( particleChannelMap.needs_scope_management() )
            throw std::runtime_error( "prt_file_header.set_channel_map() - The channel map particle layout requested "
                                      "includes non-raw data, such "
                                      "as a string.  The particle layout must be raw numeric data." );

        m_particleChannelMap = particleChannelMap;
        // We don't want to have any padding in the files we save
        // (i.e. if the particle uses 26 bytes, save 26 bytes and
        // not 28 bytes as it would be with 4 byte padding)
        m_particleChannelMap.remove_padding();
    }

    const channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    //////////////////////
    // Functions for dealing with header properties
    //////////////////////

    boost::int64_t particle_count() const { return m_particleCount; }

    void set_particle_count( boost::int64_t particleCount ) { m_particleCount = particleCount; }

    //////////////////////
    // Functions for dealing with file i/o of the header
    //////////////////////

    // Function to read in the PRT header.  The name parameter should be the filename for error messages.
    void read_header( std::istream& in, const frantic::tstring& streamName );

    // Function to read in the PRT header.  The name parameter should be the filename for error messages.
    void read_header( FILE* fin, const frantic::tstring& streamName );

    // Function to write out the PRT header.  The name parameter should be the filename for error messages.
    void write_header( std::ostream& out, const frantic::tstring& streamName );

    // Function to write out the PRT header.  The name parameter should be the filename for error messages.
    void write_header( FILE* out, const frantic::tstring& streamName );

    // Function to rewrite the particle count once the whole PRT file is written.
    // The name parameter should be the filename for error messages.
    // NOTE: This function seeks in the stream, but doesn't return the seek position to where it was before.
    void rewrite_particle_count( std::ostream& out, const frantic::tstring& streamName, boost::int64_t particleCount );

    // Function to rewrite the particle count once the whole PRT file is written.
    // The name parameter should be the filename for error messages.
    // NOTE: This function seeks in the stream, but doesn't return the seek position to where it was before.
    void rewrite_particle_count( FILE* out, const frantic::tstring& streamName, boost::int64_t particleCount );

    /**
     * Seeks the stream back into the header section of the file and updates the bounding box.
     * \note This function seeks in the stream, but doesn't return the seek position to where it was before.
     * \param out The stream that just finished writing the PRT file.
     * \param streamName Provided for error messages
     * \param bounds The boundbox of the particles.
     */
    void rewrite_particle_boundbox( std::ostream& out, const frantic::tstring& streamName,
                                    const frantic::graphics::boundbox3f& bounds );

    /**
     * \overload This overload accepts a C style file pointer for file I/O
     */
    void rewrite_particle_boundbox( FILE* out, const frantic::tstring& streamName,
                                    const frantic::graphics::boundbox3f& bounds );
};

} // namespace particles
} // namespace frantic
