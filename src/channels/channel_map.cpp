// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/channel_map.hpp>

using namespace std;

namespace frantic {
namespace channels {

void channel_map::construct_structure( char* data ) const {
    // Currently, a NULL pointer indicates an empty string, so we can get by with a single memset.
    // Note that this may change in the future, so you shouldn't rely on this behavior in other code.
    //
    // IMPORTANT NOTE: If construct_structure is made such that it could throw an exception, it MUST
    //                 clean up the partially constructed part before propagating that exception!
    memset( data, 0, m_structureSize );
}

void channel_map::destruct_structure( char* data ) const {
    if( needs_scope_management() ) {
        for( std::vector<channel>::const_iterator i = m_channels.begin(), ie = m_channels.end(); i != ie; ++i ) {
            // String channels need to be destructed specially
            if( i->data_type() == data_type_string ) {
                for( std::size_t j = 0, je = i->arity(); j != je; ++j ) {
                    destruct_channel_string( reinterpret_cast<char*>( data ) + i->offset() + j * i->primitive_size() );
                }
            }
        }
    }
}

// This fills in the particle with data from the input string array.  An example is parsing the data in a line of a .csv
// file, then using this function to convert it into a binary particle structure.
void channel_map::copy_structure_from_strings( char* destParticle,
                                               const std::vector<std::string>& srcParticleData ) const {
    if( needs_scope_management() )
        throw runtime_error( "channel_map::copy_structure_from_strings() - This function is not yet supported for "
                             "channel_map structures that require scope management." );

    unsigned srcIndex = 0;
    for( unsigned i = 0; i < m_channels.size(); ++i ) {
        for( unsigned j = 0; j < m_channels[i].arity(); ++j ) {
            if( srcIndex >= srcParticleData.size() )
                throw std::runtime_error(
                    "channel_map.copy_stucture_from_strings: Could not convert the input string array "
                    "into particle data, because it did not contain enough entries." );
            parse_channel_value_from_string( m_channels[i].data_type(), srcParticleData[srcIndex],
                                             m_channels[i].get_channel_data_pointer( destParticle, j ) );
            ++srcIndex;
        }
    }
}

// NOTE: Normally you don't want to bother specifying the offset, that's just for when you're exactly specifying the
// structure.
void channel_map::define_channel( const frantic::tstring& name, size_t arity, data_type_t dataType,
                                  std::size_t offset ) {
    if( m_channelDefinitionComplete )
        throw std::runtime_error( "channel_map.define_channel: Cannot add the channel definition for \"" +
                                  frantic::strings::to_string( name ) +
                                  "\" to the channel map after the end_channel_definition() is called." );

    if( !is_valid_channel_name( name ) )
        throw std::runtime_error( "channel_map.define_channel: Cannot add a channel named \"" +
                                  frantic::strings::to_string( name ) +
                                  "\" to the channel map, because it is an invalid identifier name." );

    for( std::vector<channel>::iterator i = m_channels.begin(); i != m_channels.end(); ++i ) {
        if( i->name() == name )
            throw std::runtime_error(
                "channel_map.define_channel: A channel definition for \"" + frantic::strings::to_string( name ) +
                "\" already exists in this channel map.  Each channel should be defined exactly once." );
    }

    m_channels.push_back( channel( name, offset, arity, dataType ) );
}

struct decreasingChannelSize {
    bool operator()( const channel& p1, const channel& p2 ) {
        // A particle channel should be sorted based on the size of the primitive type first
        // and then the arity. The following example demonstrates this ordering.
        // e.g.) int64 id, vector3f pos, vector3f vel, color3h
        //        8 bytes, 3*4 bytes   , 3*4 bytes   , 3*2 bytes
        //
        // This makes it easier to pad the particle channels so that they are aligned to 4 or 8 bytes.

        size_t size1 = sizeof_channel_data_type( p1.data_type() );
        size_t size2 = sizeof_channel_data_type( p2.data_type() );

        if( size1 == size2 ) {
            return p1.arity() > p2.arity();
        } else {
            return size1 > size2;
        }
    }
};

struct channel_name_equals {
  private:
    const channel* ch;

  public:
    channel_name_equals( const channel& c )
        : ch( &c ) {}
    bool operator()( const channel& c ) const { return c.name() == ch->name(); }
};

void channel_map::union_channel_map( const channel_map& otherChannelMap ) {
    if( channel_definition_complete() )
        throw std::runtime_error( "channel_map::union_channel_map() - Cannot union onto an already defined map" );

    bool needToRedoChannelMap = false;

    // go though each channel in the otherChannelMap channel map, and add them to our own.
    for( std::size_t i = 0; i < otherChannelMap.channel_count(); ++i ) {
        const channel& otherChannel = otherChannelMap[i];

        // find this channel in our channel map (linear search... bad?)
        int index = -1;
        for( std::size_t j = 0; j < channel_count(); ++j ) {
            if( m_channels[j].name() == otherChannel.name() ) {
                index = (int)j;
                break;
            }
        }

        if( index == -1 ) {
            // found a channel in "otherChannelMap" that is not in our channel map. insert it into ours.
            define_channel( otherChannelMap[i].name(), otherChannelMap[i].arity(), otherChannelMap[i].data_type() );
        } else {
            // found a channel in "otherChannelMap" this *is* in our channel map.
            const channel& ourChannel = m_channels[index];

            // check the arity is the same for both channels
            if( ourChannel.arity() != otherChannel.arity() )
                throw std::runtime_error(
                    "channel_map::union_channel_map() - Both channel maps have a channel named \"" +
                    frantic::strings::to_string( ourChannel.name() ) +
                    "\". These two channels can not be merged because they are of different arity." );

            // check the data type
            if( ourChannel.data_type() != otherChannel.data_type() ) {

                // check if the data type is the same "class" of data type.
                bool dataTypesMatch = false;
                if( frantic::channels::is_channel_data_type_float( ourChannel.data_type() ) ) {
                    if( frantic::channels::is_channel_data_type_float( otherChannel.data_type() ) )
                        dataTypesMatch = true;
                } else if( frantic::channels::is_channel_data_type_signed( ourChannel.data_type() ) ) {
                    if( frantic::channels::is_channel_data_type_signed( otherChannel.data_type() ) )
                        dataTypesMatch = true;
                } else if( frantic::channels::is_channel_data_type_unsigned( ourChannel.data_type() ) ) {
                    if( frantic::channels::is_channel_data_type_unsigned( otherChannel.data_type() ) )
                        dataTypesMatch = true;
                }

                // if the data types work together, check to see if our channel needs its type upgraded to a higher
                // bitsize. a higher bitsize channel means we need to recreate our channel map.
                if( dataTypesMatch ) {
                    if( otherChannel.data_type() > ourChannel.data_type() )
                        needToRedoChannelMap = true;
                } else {
                    throw std::runtime_error(
                        "channel_map::union_channel_map() - Both channel maps have a channel named \"" +
                        frantic::strings::to_string( ourChannel.name() ) +
                        "\". These two channels can not be merged because they have incompatible data "
                        "types. Our channel has data type " +
                        frantic::strings::to_string( channel_data_type_str( ourChannel.data_type() ) ) +
                        ". Other channel has data type " +
                        frantic::strings::to_string( channel_data_type_str( otherChannel.data_type() ) ) );
                }
            }
        }
    }

    // at the end of the above, we have a fully unioned channel map.
    // however, if there was one or more channels that were flagged as "needing bit upgrades", we will need to redo the
    // whole channel map.
    if( needToRedoChannelMap ) {
        channel_map newThis;

        // go though each channel in the otherChannelMap channel map, and add them to our own.
        for( std::size_t i = 0; i < channel_count(); ++i ) {
            const channel& ourChannel = m_channels[i];

            // find this channel in the other channel map (linear search... bad?)
            int index = -1;
            for( std::size_t j = 0; j < otherChannelMap.channel_count(); ++j ) {
                if( ourChannel.name() == otherChannelMap.m_channels[j].name() ) {
                    index = (int)j;
                    break;
                }
            }

            // define this new channel
            if( index == -1 ) {
                newThis.define_channel( ourChannel.name(), ourChannel.arity(), ourChannel.data_type() );
            } else {
                const channel& otherChannel = otherChannelMap.m_channels[index];
                newThis.define_channel( ourChannel.name(), ourChannel.arity(),
                                        std::max( ourChannel.data_type(),
                                                  otherChannel.data_type() ) ); // we are ensured that arity is correct,
                                                                                // and data type is promoteable.
            }
        }

        // assign "us" to be this newly created channel map
        this->swap( newThis );
    }
}

struct increasingChannelOffset {
    bool operator()( const channel& p1, const channel& p2 ) {
        // When explicitly constructed, the channels need to be sorted in the
        // order such that the offset values are increasing

        return p1.offset() < p2.offset();
    }
};

// For in-memory structures, the padding should generally be 4 or 8, while for writing to disk structures it should
// generally be 1. Set preserveChannelOrder to true if you want to define a data structure which is packed, without
// manually computing the offsets yourself. Set assignOffsets to false only if you've provided the offest value for all
// the define_channel calls, for instance when loading a particle file from disk.
void channel_map::end_channel_definition( std::size_t alignmentPadding, bool preserveChannelOrder,
                                          bool assignOffsets ) {
    m_structureAlignment = alignmentPadding;
    if( assignOffsets ) {
        if( !preserveChannelOrder )
            stable_sort( m_channels.begin(), m_channels.end(), decreasingChannelSize() );

        // for each channel compute offset
        std::size_t offset = 0;
        for( unsigned i = 0; i < m_channels.size(); ++i ) {
            m_channels[i].set_offset( offset );
            offset += m_channels[i].primitive_size();
        }

        m_structureSize = offset;
    } else {
        stable_sort( m_channels.begin(), m_channels.end(), increasingChannelOffset() );

        // Check to make sure the offsets are sane (ie. no overlapping channels)
        for( std::size_t i = 1; i < m_channels.size(); ++i ) {
            if( m_channels[i].offset() < m_channels[i - 1].offset() + m_channels[i - 1].primitive_size() )
                throw std::runtime_error( "channel_map::end_channel_definition() - The channel " +
                                          frantic::strings::to_string( m_channels[i].name() ) +
                                          " is overlapping with channel " +
                                          frantic::strings::to_string( m_channels[i - 1].name() ) );
        }

        if( m_channels.empty() ) {
            m_structureSize = 0;
        } else {
            m_structureSize = m_channels.back().offset() + m_channels.back().primitive_size();
        }
    }

    // Add the required padding to the structure size
    if( m_structureSize % alignmentPadding != 0 )
        m_structureSize += alignmentPadding - m_structureSize % alignmentPadding;

    // Create the name mapping, and determine whether the structure needs scope management.
    m_nameMapping.clear();
    m_needsScopeManagement = false;
    m_scopedChannels.clear();
    for( std::size_t i = 0; i < m_channels.size(); ++i ) {
        // Add this channel to the name mapping
        m_nameMapping[m_channels[i].name()] = (int)i;
        // If this data type isn't POD (plain old data), then the channel_map structures need scope management.
        if( !is_channel_data_type_pod( m_channels[i].data_type() ) ) {
            m_needsScopeManagement = true;
            m_scopedChannels.push_back( i );
        }
    }

    m_channelDefinitionComplete = true;
}

struct compare_channel_name {
    frantic::tstring m_name;
    compare_channel_name( const frantic::tstring& name )
        : m_name( name ) {}
    bool operator()( const channel& c ) const { return c.name() == m_name; }
};

void channel_map::rename_channel( const frantic::tstring& from, const frantic::tstring& to ) {
    if( !channel_definition_complete() )
        throw std::runtime_error( "channel_map::rename_channel() - Called on a non-finalized channel_map" );

    std::map<frantic::tstring, int>::iterator it = m_nameMapping.find( from );
    if( it == m_nameMapping.end() )
        throw std::runtime_error( "channel_map::rename_channel() - Map does not have channel \"" +
                                  frantic::strings::to_string( from ) + "\" to rename" );

    int index = it->second;
    m_nameMapping.erase( it );
    m_nameMapping.insert( std::make_pair( to, index ) );
    m_channels[index].m_name = to;
}

void channel_map::delete_channel( const frantic::tstring& channel, bool leaveGaps ) {
    if( !channel_definition_complete() )
        throw std::runtime_error( "channel_map::delete_channel() - Called on a non-final channel_map" );

    std::map<frantic::tstring, int>::iterator it = m_nameMapping.find( channel );
    if( it == m_nameMapping.end() )
        throw std::runtime_error( "channel_map::delete_channel() - Map does not have channel \"" +
                                  frantic::strings::to_string( channel ) + "\" to delete" );

    std::size_t size = m_channels[it->second].primitive_size();

    for( std::size_t i = it->second + 1; i != m_channels.size(); ++i ) {
        if( !leaveGaps )
            m_channels[i].set_offset( m_channels[i].offset() - size ); // Fix the offset

        m_nameMapping[m_channels[i].name()] -= 1; // Repoint the nameMapping to the correct index
        m_channels[i - 1] = m_channels[i];        // Push the channel backwards one
    }

    // Fix up the structure size.
    m_structureSize = 0;
    if( !m_channels.empty() ) {
        m_structureSize = m_channels.back().offset() +
                          m_channels.back().primitive_size(); // The structure ends after the last channel.

        // We might need to re-add padding to fit the alignment requirements.
        if( m_structureSize % m_structureAlignment != 0 )
            m_structureSize += m_structureAlignment - m_structureSize % m_structureAlignment;
    }

    m_channels.pop_back();
    m_nameMapping.erase( it );

    // Determine whether the structure still needs scope management.
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

std::basic_ostream<frantic::tchar>& operator<<( std::basic_ostream<frantic::tchar>& out, const channel_map& cm ) {
    out << _T("channel_map {\n");
    for( unsigned i = 0; i < cm.channel_count(); ++i ) {
        out << _T( ' ' ) << cm[i].type_str() << _T( ' ' ) << cm[i].name() << _T(";\n");
    }
    out << _T("};\n");

    return out;
}

std::string channel_map::str() const {
    std::basic_stringstream<frantic::tchar> ss;
    ss << *this;
    return frantic::strings::to_string( ss.str() );
}

void channel_map::dump( std::ostream& out ) const {
    out << "Dumping particle channel map:\n";
    for( std::vector<channel>::const_iterator i = m_channels.begin(); i != m_channels.end(); ++i ) {
        i->dump( out );
    }
    out << "Particle Size: " << m_structureSize << std::endl;
    out << "Done particle channel map dump" << std::endl;
}

void channel_map::dump_particle( std::ostream& out, const char* rawParticleBuffer ) const {
    out << "<<<Dumping particle:\n";
    for( std::vector<channel>::const_iterator i = m_channels.begin(); i != m_channels.end(); ++i ) {
        i->dump_value( out, rawParticleBuffer );
    }
    out << ">>>Done particle dump" << std::endl;
}

} // namespace channels
} // namespace frantic
