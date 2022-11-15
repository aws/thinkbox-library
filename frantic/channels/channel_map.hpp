// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <map>
#include <string>

#include <frantic/channels/channel_accessor.hpp>
#include <frantic/channels/named_channel_data.hpp>

namespace frantic {
namespace channels {

/**
 * Stores the information necessary for a particle data channel to be added to
 * a particle data collection
 */
class channel {
    friend class channel_map;

    /// string name of the channel
    frantic::tstring m_name;
    /// # of data values for this channel per particle
    std::size_t m_arity;
    /// primitive data type
    data_type_t m_dataType;
    /// total channel size = arity* the size of m_dataType
    std::size_t m_primitiveSize;
    /// The offset for this channel, relative to the start of a particle
    std::size_t m_offset;

    channel() {}

  public:
    channel( const frantic::tstring& n, std::size_t off, std::size_t arity, data_type_t dataType ) {
        m_name = n;
        m_offset = off;
        m_arity = arity;
        m_dataType = dataType;
        m_primitiveSize = m_arity * sizeof_channel_data_type( m_dataType );
    }

    void swap( channel& rhs ) {
        m_name.swap( rhs.m_name );
        std::swap( m_arity, rhs.m_arity );
        std::swap( m_dataType, rhs.m_dataType );
        std::swap( m_primitiveSize, rhs.m_primitiveSize );
        std::swap( m_offset, rhs.m_offset );
    }

    const frantic::tstring& name() const { return m_name; }

    std::size_t offset() const { return m_offset; }

    void set_offset( std::size_t offset ) { m_offset = offset; }

    std::size_t arity() const { return m_arity; }
    data_type_t data_type() const { return m_dataType; }

    /**
     * This returns the size of the primitive of this channel.  That's equivalent to arity() * (data type size).
     */
    std::size_t primitive_size() const { return m_primitiveSize; }

    /**
     * Given a pointer to a particle buffer, this returns a pointer to the
     * channel within that buffer
     */
    const char* get_channel_data_pointer( const void* rawParticleData ) const {
        return reinterpret_cast<const char*>( rawParticleData ) + m_offset;
    }

    /**
     * Given a pointer to a particle buffer, this returns a pointer to the
     * channel within that buffer
     */
    char* get_channel_data_pointer( void* rawParticleData ) const {
        return reinterpret_cast<char*>( rawParticleData ) + m_offset;
    }

    /**
     * Given a pointer to a particle buffer, this returns a pointer to the channel
     * within that buffer, indexed within the channel.
     */
    char* get_channel_data_pointer( char* rawParticleData, int index ) const {
        return rawParticleData + m_offset + index * sizeof_channel_data_type( m_dataType );
    }

    /**
     * Given a pointer to a particle buffer, this returns a pointer to the channel
     * within that buffer, indexed within the channel.
     */
    const char* get_channel_data_pointer( const char* rawParticleData, int index ) const {
        return rawParticleData + m_offset + index * sizeof_channel_data_type( m_dataType );
    }

    frantic::tstring type_str() const { return channel_data_type_str( m_arity, m_dataType ); }

    bool operator==( const channel& rhs ) const {
        return m_arity == rhs.m_arity && m_dataType == rhs.m_dataType && m_offset == rhs.m_offset &&
               m_name == rhs.m_name;
    }

    bool operator!=( const channel& rhs ) const { return !( operator==( rhs ) ); }

    void dump( std::ostream& out ) const {
        out << "Named Channel:\n";
        out << "Name:     " << frantic::strings::to_string( m_name ) << "\n";
        out << "Offset:   " << (unsigned int)m_offset << "\n";
        out << "Data Type: " << frantic::strings::to_string( channel_data_type_str( m_arity, m_dataType ) ) << "\n";
    }

    void dump_value( std::ostream& out, const char* rawParticleData ) const {
        out << frantic::strings::to_string( m_name ) << ": ";
        channels::channel_data_type_print( out, ",", m_arity, m_dataType, rawParticleData + m_offset );
        out << "\n";
    }
};

class channel_map {
    /// This maps channel names to the index in the m_channels vector
    std::map<frantic::tstring, int> m_nameMapping;
    /// This contains all the channels, in order of their appearance within the structure.
    std::vector<channel> m_channels;
    /// This contains the indices of all the particular channels that need scope management.
    std::vector<std::size_t> m_scopedChannels;
    /// This is the total size of one particle
    std::size_t m_structureSize;
    /// Total structure alignment
    std::size_t m_structureAlignment;
    bool m_channelDefinitionComplete;
    /// This boolean flag indicates whether a memory block represented by this channel_map needs scope management.
    bool m_needsScopeManagement;

  public:
    channel_map() {
        m_channelDefinitionComplete = false;
        m_needsScopeManagement = true;
        m_structureSize =
            0; // It was undefined before, which is ok but makes for some really irritating bugs when it works in debug
               // mode but not release. Maybe this could be UINT_MAX so it causes *defined* stupid behavior.
        m_structureAlignment = 1; // It was undefined before, which is ok but makes for some really irritating bugs when
                                  // it works in debug mode but not release.
    }

    // NOTE: It is extremely important that any *data members* added to channel_map are swapped correctly!!!!
    void swap( channel_map& rhs ) {
        m_nameMapping.swap( rhs.m_nameMapping );
        m_channels.swap( rhs.m_channels );
        m_scopedChannels.swap( rhs.m_scopedChannels );
        std::swap( m_structureSize, rhs.m_structureSize );
        std::swap( m_structureAlignment, rhs.m_structureAlignment );
        std::swap( m_channelDefinitionComplete, rhs.m_channelDefinitionComplete );
        std::swap( m_needsScopeManagement, rhs.m_needsScopeManagement );
    }

    ///////////////
    // Access to information about channels
    ///////////////

    /**
     * This initializes an instance of this structure.  If needs_scope_management() is true,
     * you must use this and destruct_structure around the lifetime of an instance.
     *
     * @param  data  A pointer to the memory comprising the structure.
     */
    void construct_structure( char* data ) const;

    /**
     * @overload
     */
    void construct_structure( std::vector<char>& data ) const { construct_structure( &data[0] ); }

    /**
     * This destroys an instance of this structure.  If needs_scope_management() is true,
     * you must use this and destruct_structure around the lifetime of an instance.  It
     * deallocates the memory for strings and other elements which also contain memory
     * from the heap.
     *
     * @param  data  A pointer to the memory comprising the structure.
     */
    void destruct_structure( char* data ) const;

    /**
     * @overload
     */
    void destruct_structure( std::vector<char>& data ) const { destruct_structure( &data[0] ); }

    void copy_structure( char* destData, const char* sourceData ) const {
        // If the structure needs scope management, we have to do a more complicated copy.
        if( needs_scope_management() ) {
            std::size_t prevOffset = 0;
            for( std::vector<std::size_t>::const_iterator i = m_scopedChannels.begin(), ie = m_scopedChannels.end();
                 i != ie; ++i ) {
                std::size_t index = *i, curOffset = m_channels[index].offset();

                // Copy everything up to this scoped channel
                if( prevOffset < curOffset )
                    memcpy( destData + prevOffset, sourceData + prevOffset, curOffset - prevOffset );

                // Copy the data for the scoped channel
                switch( m_channels[index].data_type() ) {
                case data_type_string:
                    copy_channel_string( destData + curOffset, sourceData + curOffset, m_channels[index].arity() );
                    break;
                default:
                    throw std::runtime_error(
                        "channel_map::copy_structure() - Channel \"" +
                        frantic::strings::to_string( m_channels[index].name() ) +
                        "\" was flagged as needing scope management, but it is of type " +
                        frantic::strings::to_string( channel_data_type_str( m_channels[index].data_type() ) ) +
                        ", which doesn't have a copy mechanism implemented." );
                }

                // This is where to continue copying later
                prevOffset = curOffset + m_channels[index].primitive_size();
            }
            // Copy the last bit
            if( prevOffset < m_structureSize )
                memcpy( destData + prevOffset, sourceData + prevOffset, m_structureSize - prevOffset );
        } else {
            memcpy( destData, sourceData, m_structureSize );
        }
    }

    void copy_structure( std::vector<char>& destData, const char* sourceData ) const {
        // TODO: Can do checking in a DEBUG mode.
        copy_structure( &destData[0], sourceData );
    }

    void copy_structure( char* destData, const std::vector<char>& sourceData ) const {
        // TODO: Can do checking in a DEBUG mode.
        copy_structure( destData, &sourceData[0] );
    }

    void copy_structure( std::vector<char>& destData, const std::vector<char>& sourceData ) const {
        // TODO: Can do checking in a DEBUG mode.
        copy_structure( &destData[0], &sourceData[0] );
    }

    /**
     * This fills in the particle with data from the input string array.  An example is
     * parsing the data in a line of a .csv file, then using this function to convert it
     * into a binary particle structure.
     */
    void copy_structure_from_strings( char* destParticle, const std::vector<std::string>& srcParticleData ) const;

    size_t channel_count() const { return m_channels.size(); }

    bool has_channel( const frantic::tstring& name ) const {
        if( m_channelDefinitionComplete ) {
            return m_nameMapping.find( name ) != m_nameMapping.end();
        } else {
            std::vector<channel>::const_iterator i;
            for( i = m_channels.begin(); i != m_channels.end() && i->name() != name; ++i )
                ;
            return i != m_channels.end();
        }
    }

    // allows the user to query what a channel is, based on a channel index.
    // this function allows the user to iterate a channel_map, and discover the details of each channel.
    void get_channel_definition( const std::size_t& channelIndex, frantic::tstring& outName, data_type_t& outDataType,
                                 size_t& outArity ) const {
        if( channelIndex >= m_channels.size() )
            throw std::runtime_error( "channel_map::get_channel_definition(): Index is out of bounds." );
        outName = m_channels[channelIndex].name();
        outDataType = m_channels[channelIndex].data_type();
        outArity = m_channels[channelIndex].arity();
    }

    // allows the use to query what a channel is, based on a channel name.
    // this function allows users to find out channel details without using a channel_general_accessor (for example, if
    // the channel_map definition is not yet complete).
    void get_channel_definition( const frantic::tstring& channelName, data_type_t& outDataType, size_t& outArity,
                                 size_t& outOffset ) const {
        if( m_channelDefinitionComplete ) {
            std::map<frantic::tstring, int>::const_iterator iter;
            iter = m_nameMapping.find( channelName );
            if( iter == m_nameMapping.end() )
                throw std::runtime_error( "channel_map::get_channel_definition(): \"" +
                                          frantic::strings::to_string( channelName ) +
                                          "\" was not found in channel map." );
            const frantic::channels::channel& c = m_channels[iter->second];
            outDataType = c.data_type();
            outArity = c.arity();
            outOffset = c.offset();
        } else {
            std::vector<channel>::const_iterator iter;
            for( iter = m_channels.begin(); iter != m_channels.end() && iter->name() != channelName; ++iter )
                ;
            if( iter == m_channels.end() )
                throw std::runtime_error( "channel_map::get_channel_definition(): \"" +
                                          frantic::strings::to_string( channelName ) +
                                          "\" was not found in channel map." );
            outDataType = iter->data_type();
            outArity = iter->arity();
            outOffset = iter->offset();
        }
    }

    void get_channel_definition( const frantic::tstring& channelName, data_type_t& outDataType,
                                 size_t& outArity ) const {
        size_t outOffset;
        get_channel_definition( channelName, outDataType, outArity, outOffset );
    }

    std::size_t structure_size() const { return m_structureSize; }

    // This method returns the size of the structure without the inter-structure padding. Ie. If there are a junk 3bytes
    // at the end of each struture to ensure proper alignment of an adjacent structure, this will return
    // structure_size()
    // - 3.
    //  @note This is only valid if channel_definition_complete() == true.
    //  @note This is similar to remove_padding() except it does not change the structure, and it takes advantage of the
    //  sorted nature of m_channels.
    std::size_t structure_size_without_padding() const {
        if( !m_channels.empty() ) {
            return m_channels.back().offset() + m_channels.back().primitive_size();
        } else {
            return 0;
        }
    }

    /**
     * Low level function which gets the offset of the named channel within the particle structure.
     */
    size_t channel_offset( const frantic::tstring& name ) const { return m_channels[channel_index( name )].offset(); }

    /**
     * Return the index of the given named channel in the map
     *
     * @param name the channel to find the index for
     *
     */
    size_t channel_index( const frantic::tstring& name ) const {
        if( !m_channelDefinitionComplete )
            throw std::runtime_error( "channel_map.channel_offset: Cannot return an accessor for channel \"" +
                                      frantic::strings::to_string( name ) +
                                      "\" until end_channel_definition() is called" );

        std::map<frantic::tstring, int>::const_iterator i = m_nameMapping.find( name );

        if( i == m_nameMapping.end() )
            throw std::runtime_error( "channel_map.channel_offset: This particle channel map has no channel named \"" +
                                      frantic::strings::to_string( name ) + "\"" );

        return i->second;
    }

    /**
     * This returns true if the channel_map structure has been fully defined.  It returns false if
     * the structure is in the process of being defined.
     */
    bool channel_definition_complete() const { return m_channelDefinitionComplete; }

    /**
     * This returns true is a memory block following the layout of this channel_map will need its scope managed.
     * If this returns true, then the construct_structure should be called when memory is allocated, and
     * destruct_structure should be called when it is deallocated.
     */
    bool needs_scope_management() const { return m_needsScopeManagement; }

    const channel& operator[]( std::size_t index ) const { return m_channels[index]; }

    const channel& operator[]( const frantic::tstring& name ) const {
        std::map<frantic::tstring, int>::const_iterator i = m_nameMapping.find( name );

        if( i == m_nameMapping.end() )
            throw std::runtime_error( "channel_map.operator[]: This particle channel map has no channel named \"" +
                                      frantic::strings::to_string( name ) + "\"" );

        return m_channels[i->second];
    }

    bool operator==( const channel_map& rhs ) const {
        return m_structureSize == rhs.m_structureSize && m_channels == rhs.m_channels &&
               m_nameMapping == rhs.m_nameMapping;
    }

    bool operator!=( const channel_map& rhs ) const {
        return m_structureSize != rhs.m_structureSize || m_channels != rhs.m_channels ||
               m_nameMapping != rhs.m_nameMapping;
    }

    ///////////////
    // Functions for creating the channel map
    ///////////////

    /**
     * This resets the channel_map to nothing, allowing it to be recreated
     */
    void reset() {
        m_channelDefinitionComplete = false;
        m_needsScopeManagement = true;
        m_structureSize = 0;
        m_nameMapping.clear();
        m_channels.clear();
    }

    /**
     * Convenience function which translates a compile time type into the
     * appropriate enum values.
     */
    template <class T>
    void define_channel( const frantic::tstring& name ) {
        define_channel( name, channel_data_type_traits<T>::arity(), channel_data_type_traits<T>::data_type() );
    }

    /**
     * @note Normally you don't want to bother specifying the offset, that's just for when you're exactly specifying the
     * structure.
     */
    void define_channel( const frantic::tstring& name, size_t arity, data_type_t dataType, std::size_t offset = 0 );

    /**
     * For in-memory structures, the padding should generally be 4 or 8, while for writing to disk structures it should
     * generally be 1. Set preserveChannelOrder to true if you want to define a data structure which is packed, without
     * manually computing the offsets yourself. Set assignOffsets to false only if you've provided the offest value for
     * all the define_channel calls, for instance when loading a particle file from disk.
     */
    void end_channel_definition( std::size_t alignmentPadding = 4, bool preserveChannelOrder = false,
                                 bool assignOffsets = true );

    /**
     * This function allows a channel to be modified after calling define_channel()
     * @note Calling this on a channel map for which you have specified explicit offsets
     *        is a very BAD idea.
     * @param channelName The name of the channel to modify.
     * @param type The new type of the channel to modify.
     * @param arity The new arity of the channel to modify. If 0, then the old channel's arity will be used.
     */
    void set_channel_data_type( const frantic::tstring& channelName, data_type_t type, std::size_t arity = 0 ) {
        if( channel_definition_complete() )
            throw std::runtime_error(
                "channel_map::set_channel_data_type() - Cannot change the channel type after it has been finalized." );

        std::vector<channel>::iterator it = m_channels.begin(), itEnd = m_channels.end();
        for( ; it != itEnd; ++it ) {
            if( it->name() == channelName ) {
                ( *it ) = channel( it->name(), it->offset(), arity == 0 ? it->arity() : arity, type );
                return;
            }
        }

        throw std::runtime_error( "channel_map::set_channel_data_type() - The channel \"" +
                                  frantic::strings::to_string( channelName ) + "\" was not found" );
    }

    /**
     * This function modifies the channel map, shrinking the size to remove any alignment padding.
     */
    void remove_padding() {
        std::size_t particleSize = 0;
        for( std::size_t i = 0; i < m_channels.size(); ++i ) {
            std::size_t sizeCandidate = m_channels[i].offset() + m_channels[i].primitive_size();
            if( sizeCandidate > particleSize )
                particleSize = sizeCandidate;
        }
        m_structureSize = particleSize;
    }

    ///////////////
    // Functions for modifying a finalized channel map
    ///////////////

    /**
     * Appends a channel to the end of the structure. Works on finalized channels
     */
    void append_channel( const frantic::tstring& name, size_t arity, data_type_t dataType ) {
        if( !channel_definition_complete() )
            throw std::runtime_error( "channel_map::append_channel() - Called on a non-finalized channel_map" );

        if( m_nameMapping.find( name ) != m_nameMapping.end() )
            throw std::runtime_error( "channel_map::append_channel() - Map already has channel: " +
                                      frantic::strings::to_string( name ) );

        m_nameMapping[name] = (int)m_channels.size(); // Cast for x64

        if( m_channels.size() == 0 ) {
            m_channels.push_back( channel( name, 0, arity, dataType ) );
        } else {
            const channel& ch = m_channels.back();
            m_channels.push_back( channel( name, ch.offset() + ch.primitive_size(), arity, dataType ) );
        }

        m_structureSize = m_channels.back().offset() + m_channels.back().primitive_size();
        if( m_structureSize % m_structureAlignment != 0 )
            m_structureSize += ( m_structureAlignment - m_structureSize % m_structureAlignment );

        // Determine whether the structure needs scope management after the change.
        m_needsScopeManagement = false;
        m_scopedChannels.clear();
        for( std::size_t i = 0; i < m_channels.size(); ++i ) {
            // If this data type isn't POD (plain old data), then the channel_map structures need scope management.
            if( !is_channel_data_type_pod( m_channels[i].data_type() ) ) {
                m_needsScopeManagement = true;
                m_scopedChannels.push_back( i );
            }
        }
    }

    template <class T>
    void append_channel( const frantic::tstring& name ) {
        append_channel( name, channel_data_type_traits<T>::arity(), channel_data_type_traits<T>::data_type() );
    }

    /**
     * Rename a channel in a map.
     */
    void rename_channel( const frantic::tstring& from, const frantic::tstring& to );

    /**
     * Deletes a channel from a map. Can leave a gap is leaveGaps is set to true.
     */
    void delete_channel( const frantic::tstring& channel, bool leaveGaps = false );

    /**
     * Sets this channel map to be the union of its current channels and the channels of the supplied map. Will promote
     *data types data types if they don't match and one can be promoted to the other. This may break the supplied
     *offsets if you plan to call end_channel_definition() with assignOffsets=false.
     */
    void union_channel_map( const channel_map& other );

    ///////////////
    // Functions for retrieving channel accessors
    ///////////////

    template <class T>
    channel_accessor<T> get_accessor( const frantic::tstring& name ) const {
        if( !m_channelDefinitionComplete )
            throw std::runtime_error( "channel_map.get_accessor: Cannot return an accessor for channel \"" +
                                      frantic::strings::to_string( name ) +
                                      "\" until end_channel_definition() is called" );

        std::map<frantic::tstring, int>::const_iterator i = m_nameMapping.find( name );

        if( i == m_nameMapping.end() )
            throw std::runtime_error( "channel_map.get_accessor: There is no channel named \"" +
                                      frantic::strings::to_string( name ) + "\"." );

        if( channel_data_type_traits<T>::data_type() != m_channels[i->second].data_type() ||
            channel_data_type_traits<T>::arity() != m_channels[i->second].arity() ) {
            throw std::runtime_error(
                "channel_map.get_accessor: The type requested is not a valid type for the channel \"" +
                frantic::strings::to_string( name ) +
                "\"."
                " The requested type is " +
                frantic::strings::to_string( channel_data_type_traits<T>::type_str() ) +
                ","
                " the stored type is " +
                frantic::strings::to_string(
                    channel_data_type_str( m_channels[i->second].arity(), m_channels[i->second].data_type() ) ) +
                "." );
        }

        return channel_accessor<T>( m_channels[i->second].offset() );
    }

    template <class T>
    channel_cvt_accessor<T> get_cvt_accessor( const frantic::tstring& name ) const {
        if( !m_channelDefinitionComplete )
            throw std::runtime_error( "channel_map.get_cvt_accessor: Cannot return an accessor for channel \"" +
                                      frantic::strings::to_string( name ) +
                                      "\" until end_channel_definition() is called" );

        std::map<frantic::tstring, int>::const_iterator i = m_nameMapping.find( name );

        if( i == m_nameMapping.end() )
            throw std::runtime_error( "channel_map.get_cvt_accessor: There is no channel named \"" +
                                      frantic::strings::to_string( name ) + "\"" );

        const channel& ch = m_channels[i->second];
        return channel_cvt_accessor<T>( ch.offset(), ch.arity(), ch.data_type(), ch.name() );
    }

    template <class T>
    channel_const_cvt_accessor<T> get_const_cvt_accessor( const frantic::tstring& name ) const {
        if( !m_channelDefinitionComplete )
            throw std::runtime_error( "channel_map.get_const_cvt_accessor: Cannot return an accessor for channel \"" +
                                      frantic::strings::to_string( name ) +
                                      "\" until end_channel_definition() is called" );

        std::map<frantic::tstring, int>::const_iterator i = m_nameMapping.find( name );

        if( i == m_nameMapping.end() )
            throw std::runtime_error( "channel_map.get_const_cvt_accessor: There is no channel named \"" +
                                      frantic::strings::to_string( name ) + "\"" );

        const channel& ch = m_channels[i->second];
        return channel_const_cvt_accessor<T>( ch.offset(), ch.arity(), ch.data_type(), ch.name() );
    }

    channel_general_accessor get_general_accessor( const frantic::tstring& name ) const {
        if( !m_channelDefinitionComplete )
            throw std::runtime_error( "channel_map.get_general_accessor: Cannot return an accessor for channel \"" +
                                      frantic::strings::to_string( name ) +
                                      "\" until end_channel_definition() is called" );

        std::map<frantic::tstring, int>::const_iterator i = m_nameMapping.find( name );

        if( i == m_nameMapping.end() )
            throw std::runtime_error( "channel_map.get_general_accessor: There is no channel named \"" +
                                      frantic::strings::to_string( name ) + "\"" );

        return channel_general_accessor( m_channels[i->second].offset(), m_channels[i->second].arity(),
                                         m_channels[i->second].data_type() );
    }

    channel_general_accessor get_general_accessor( const std::size_t& channelIndex ) const {
        if( !m_channelDefinitionComplete )
            throw std::runtime_error( "channel_map.get_general_accessor: Cannot return an accessor for channel " +
                                      boost::lexical_cast<std::string>( channelIndex ) +
                                      " until end_channel_definition() is called" );
        return channel_general_accessor( m_channels[channelIndex].offset(), m_channels[channelIndex].arity(),
                                         m_channels[channelIndex].data_type() );
    }

    ///////////////
    // I/O functions
    ///////////////

    std::string str() const;

    ///////////////
    // Debugging functions
    ///////////////

    void dump( std::ostream& out ) const;

    void dump_particle( std::ostream& out, const char* rawParticleBuffer ) const;
};

std::basic_ostream<frantic::tchar>& operator<<( std::basic_ostream<frantic::tchar>& out, const channel_map& cm );

} // namespace channels
} // namespace frantic
