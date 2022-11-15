// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/property_map.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>

using namespace std;
using namespace boost;

namespace frantic {
namespace channels {

namespace detail {
/**
 * This is a small RAII helper class, which is used by the copy constructor, for example, to avoid a potential
 * memory leak if channel_map::construct_structure or channel_map::copy_structure were to throw an exception.
 */
class property_map_mem_holder {
    char* m_buffer;

  public:
    explicit property_map_mem_holder( char* buffer )
        : m_buffer( buffer ) {}

    ~property_map_mem_holder() {
        if( m_buffer ) {
            delete[] m_buffer;
            m_buffer = 0;
        }
    }

    char* get() { return m_buffer; }

    char* release() {
        char* temp = m_buffer;
        m_buffer = 0;
        return temp;
    }
};
} // namespace detail

property_map::property_map() {
    m_channelMap.end_channel_definition();
    m_buffer = 0;
}

// Copy constructor
property_map::property_map( const property_map& rhs ) {
    // Initialize the property_map to empty.
    m_channelMap.end_channel_definition();
    m_buffer = 0;
    // If the rhs property_map is not empty, then we have to also copy the data.
    if( !rhs.empty() ) {
        // First construct a property_map using the channel map of the rhs.  This provides scope
        // for calling construct_structure on allocation, and will call destruct_structure on that memory
        // if rhs.copy_structure throws an exception for any reason.
        property_map temp( rhs.m_channelMap );
        // Copy the structure from rhs to the temp buffer.
        rhs.m_channelMap.copy_structure( temp.m_buffer, rhs.m_buffer );
        // Finally, swap the temp structure into *this
        temp.swap( *this );
    }
}

// Construct a 0-valued property map using the given layout
property_map::property_map( const channel_map& channelMap )
    : m_channelMap( channelMap ) {
    if( m_channelMap.structure_size() == 0 ) {
        m_buffer = 0;
    } else {
        // Use the property_map_mem_holder to protect this allocated memory in case construct_structure throws.
        detail::property_map_mem_holder mem( new char[m_channelMap.structure_size()] );
        m_channelMap.construct_structure( mem.get() );
        m_buffer = mem.release();
    }
}

property_map::~property_map() {
    if( m_buffer ) {
        m_channelMap.destruct_structure( m_buffer );
        delete[] m_buffer;
        m_buffer = 0;
    }
}

void property_map::clear() {
    if( m_buffer ) {
        m_channelMap.destruct_structure( m_buffer );
        delete[] m_buffer;
        m_buffer = 0;
    }
    m_channelMap.reset();
    m_channelMap.end_channel_definition();
}

void property_map::swap( property_map& rhs ) {
    m_channelMap.swap( rhs.m_channelMap );
    std::swap( m_buffer, rhs.m_buffer );
}

property_map& property_map::operator=( const property_map& rhs ) {
    // Create a new property_map, initializing it to equal the rhs
    property_map newPropertyMap( rhs );
    // Swap that property map into *this
    newPropertyMap.swap( *this );
    return *this;
}

void property_map::dump( std::ostream& out ) const {
    for( size_t i = 0; i < m_channelMap.channel_count(); ++i ) {
        out << frantic::strings::to_string( m_channelMap[i].name() ) << " -> ";
        channels::channel_data_type_print( out, ", ", m_channelMap[i].arity(), m_channelMap[i].data_type(),
                                           m_channelMap[i].get_channel_data_pointer( m_buffer ) );
        out << "\n";
    }
}

void property_map::dump( std::wostream& out ) const {
    for( size_t i = 0; i < m_channelMap.channel_count(); ++i ) {
        out << frantic::strings::to_wstring( m_channelMap[i].name() ) << L" -> ";
        channels::channel_data_type_print( out, L", ", m_channelMap[i].arity(), m_channelMap[i].data_type(),
                                           m_channelMap[i].get_channel_data_pointer( m_buffer ) );
        out << L"\n";
    }
}

const channel_map& property_map::get_channel_map() const { return m_channelMap; }

void property_map::set_channel_map( const channel_map& cm ) {
    if( !cm.channel_definition_complete() )
        throw runtime_error( "property_map.set_channel_map() - The channel map provided has not been completed." );

    if( m_buffer ) {
        if( cm.structure_size() != 0 ) {
            // This initializes a new empty property_map, which will have m_buffer == NULL so will follow the other code
            // path inside of set_channel_map
            property_map newPropertyMap;
            newPropertyMap.set_channel_map( cm );

            // In this case, first we create the adaptor to remap the parameters (may throw an exception),
            // then we allocate a new buffer and copy the existing data into the new data layout.
            channel_map_adaptor cma( cm, m_channelMap );

            cma.copy_structure( newPropertyMap.m_buffer, m_buffer );
            newPropertyMap.swap( *this );
        } else {
            // If the channel map being set is empty, then delete the buffer
            m_channelMap.destruct_structure( m_buffer );
            delete[] m_buffer;
            m_buffer = 0;
            m_channelMap = cm;
        }
    } else {
        // If no buffer was previously allocated, copy the channel map, and allocate the buffer if necessary.
        m_channelMap = cm;
        if( m_channelMap.structure_size() != 0 ) {
            // Use the property_map_mem_holder to protect this allocated memory in case construct_structure throws.
            detail::property_map_mem_holder mem( new char[m_channelMap.structure_size()] );
            m_channelMap.construct_structure( mem.get() );
            m_buffer = mem.release();
        }
    }
}

void property_map::set_channel_map_with_swap( channel_map& cm ) {
    if( !cm.channel_definition_complete() )
        throw runtime_error(
            "property_map.set_channel_map_with_swap() - The channel map provided has not been completed." );

    if( m_buffer ) {
        if( cm.structure_size() != 0 ) {
            // This initializes a new empty property_map, which will have m_buffer == NULL so will follow the other code
            // path inside of set_channel_map
            property_map newPropertyMap;
            newPropertyMap.set_channel_map_with_swap( cm );

            // In this case, first we create the adaptor to remap the parameters (may throw an exception),
            // then we allocate a new buffer and copy the existing data into the new data layout.
            channel_map_adaptor cma( newPropertyMap.get_channel_map(), m_channelMap );

            cma.copy_structure( newPropertyMap.m_buffer, m_buffer );
            newPropertyMap.swap( *this );
        } else {
            // If the channel map being set is empty, then delete the buffer
            m_channelMap.destruct_structure( m_buffer );
            delete[] m_buffer;
            m_buffer = 0;
            m_channelMap.swap( cm );
        }
    } else {
        // If no buffer was previously allocated, copy the channel map, and allocate the buffer if necessary.
        m_channelMap.swap( cm );
        if( m_channelMap.structure_size() != 0 ) {
            // Use the property_map_mem_holder to protect this allocated memory in case construct_structure throws.
            detail::property_map_mem_holder mem( new char[m_channelMap.structure_size()] );
            m_channelMap.construct_structure( mem.get() );
            m_buffer = mem.release();
        }
    }
}

bool property_map::has_property( const frantic::tstring& propertyName ) const {
    return m_channelMap.has_channel( propertyName );
}

void property_map::merge_property_map( const property_map& rhs ) {
    // Build the unioned channel map
    channel_map cm;
    cm.union_channel_map( m_channelMap );
    cm.union_channel_map( rhs.m_channelMap );
    cm.end_channel_definition();

    // Create a temporary property_map, and set its channel map, using
    // swap instead of copy to avoid an extra copy of the channel map.
    property_map newPropertyMap;
    newPropertyMap.set_channel_map_with_swap( cm );

    // First copy the parameters from *this
    channel_map_adaptor cma( newPropertyMap.m_channelMap, m_channelMap );
    cma.copy_structure( newPropertyMap.m_buffer, m_buffer );

    // Then copy the parameters from rhs
    cma.set( newPropertyMap.m_channelMap, rhs.m_channelMap );
    cma.copy_structure( newPropertyMap.m_buffer, rhs.m_buffer );

    // Swap the completed property_map into *this
    newPropertyMap.swap( *this );
}

void property_map::delete_property( const frantic::tstring& propertyName ) {
    if( !has_property( propertyName ) ) {
        return;
    }
    channel_map cm = get_channel_map();
    cm.delete_channel( propertyName );
    set_channel_map_with_swap( cm );
}

} // namespace channels
} // namespace frantic
