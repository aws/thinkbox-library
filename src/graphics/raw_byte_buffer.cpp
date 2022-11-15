// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics/raw_byte_buffer.hpp>

using namespace std;
using namespace boost;

namespace frantic {
namespace graphics {

// Copy constructor
raw_byte_buffer::raw_byte_buffer( const raw_byte_buffer& rhs ) {
    if( rhs.m_data ) {
        // QUESTION: Maybe it would be better to replicate the same capacity as well?

        // Allocate the buffer, trimming the capacity down to the size actually used
        size_type dataSize = *reinterpret_cast<const size_type*>( rhs.m_data );
        m_data = (char*)malloc( 2 * sizeof( size_type ) + dataSize );
        if( m_data == 0 )
            throw std::runtime_error( "raw_byte_buffer: Failed to malloc " +
                                      lexical_cast<string>( 2 * sizeof( size_type ) + dataSize ) +
                                      " bytes of memory." );
        // With this trimming, the capacity equals the datasize, so set both to the same thing
        ( *reinterpret_cast<size_type*>( m_data ) ) = dataSize;
        ( *reinterpret_cast<size_type*>( m_data + sizeof( size_type ) ) ) = dataSize;
        memcpy( m_data + 2 * sizeof( size_type ), rhs.m_data + 2 * sizeof( size_type ), dataSize );
    } else {
        m_data = 0;
    }
}

raw_byte_buffer::raw_byte_buffer( const void* data, std::size_t size ) {
    m_data = 0;

    if( size > 0 ) {
        if( !data ) {
            throw std::runtime_error( "raw_byte_buffer: data is NULL" );
        }

        resize( size );
        memcpy( begin(), data, size );
    }
}

// Assignment operator
raw_byte_buffer& raw_byte_buffer::operator=( const raw_byte_buffer& rhs ) {
    if( rhs.m_data ) {
        size_type dataSize = *reinterpret_cast<const size_type*>( rhs.m_data );
        if( m_data ) {
            // Check whether we already have enough capacity for the data
            size_type currentCapacity = *reinterpret_cast<const size_type*>( m_data + sizeof( size_type ) );
            if( dataSize <= currentCapacity ) {
                // QUESTION: Maybe we should reallocate if this assignment would shrink the buffer significantly?
                // Set the data
                ( *reinterpret_cast<size_type*>( m_data ) ) = dataSize;
                memcpy( m_data + 2 * sizeof( size_type ), rhs.m_data + 2 * sizeof( size_type ), dataSize );
            } else {
                // Delete the buffer, set the pointer to zero, and fall through to the null m_data version
                free( m_data );
                m_data = 0;
            }
        }
        // This isn't in an 'else', to allow the case where we deleted the buffer to fall through to this code.
        if( !m_data ) {
            // QUESTION: Maybe it would be better to replicate the same capacity as well?

            // Allocate the buffer, trimming the capacity down to the size actually used
            m_data = (char*)malloc( 2 * sizeof( size_type ) + dataSize );
            if( m_data == 0 )
                throw std::runtime_error( "raw_byte_buffer.operator=: Failed to malloc " +
                                          lexical_cast<string>( 2 * sizeof( size_type ) + dataSize ) +
                                          " bytes of memory." );
            // With this trimming, the capacity equals the datasize, so set both to the same thing
            ( *reinterpret_cast<size_type*>( m_data ) ) = dataSize;
            ( *reinterpret_cast<size_type*>( m_data + sizeof( size_type ) ) ) = dataSize;
            memcpy( m_data + 2 * sizeof( size_type ), rhs.m_data + 2 * sizeof( size_type ), dataSize );
        }
    } else if( m_data ) {
        free( m_data );
        m_data = 0;
    }
    return *this;
}

// Equality comparison
bool raw_byte_buffer::operator==( const raw_byte_buffer& rhs ) const {
    if( m_data ) {
        if( rhs.m_data ) {
            size_type dataSize = *reinterpret_cast<const size_type*>( m_data );
            if( dataSize == *reinterpret_cast<const size_type*>( rhs.m_data ) ) {
                // memcpy returns 0 if the two buffers match
                return memcmp( m_data + 2 * sizeof( size_type ), rhs.m_data + 2 * sizeof( size_type ), dataSize ) == 0;
            } else {
                return false;
            }
        } else {
            // Need to return true if the data size is zero.
            return *reinterpret_cast<const size_type*>( m_data ) == 0;
        }
    } else {
        // Need to return true if either rhs's container is NULL, or it has data size set to zero.
        return rhs.empty();
    }
}

void raw_byte_buffer::resize( size_type requestedSize ) {
    if( m_data ) {
        size_type currentCapacity = *reinterpret_cast<const size_type*>( m_data + sizeof( size_type ) );
        if( requestedSize <= currentCapacity ) {
            // If the request shrunk the size, just set the size member
            *reinterpret_cast<size_type*>( m_data ) = requestedSize;
        } else {
            // In this case we need to allocate a bigger buffer
            char* newData = (char*)realloc( m_data, 2 * sizeof( size_type ) + requestedSize );
            if( newData == 0 )
                throw std::runtime_error( "raw_byte_buffer.resize: Failed to realloc " +
                                          lexical_cast<string>( 2 * sizeof( size_type ) + requestedSize ) +
                                          " bytes of memory." );
            m_data = newData;

            // Set the size member
            *reinterpret_cast<size_type*>( m_data ) = requestedSize;
            // Set the capacity member
            *reinterpret_cast<size_type*>( m_data + sizeof( size_type ) ) = requestedSize;
        }
    } else {
        // Allocate exactly the requested size, use add_element if you want the doubling allocation behavior
        m_data = (char*)malloc( 2 * sizeof( size_type ) + requestedSize );
        if( m_data == 0 )
            throw std::runtime_error( "raw_byte_buffer.resize: Failed to malloc " +
                                      lexical_cast<string>( 2 * sizeof( size_type ) + requestedSize ) +
                                      " bytes of memory." );
        // Set the data size
        *reinterpret_cast<std::size_t*>( m_data ) = requestedSize;
        // Set the capacity to the value requested
        *reinterpret_cast<std::size_t*>( m_data + sizeof( size_type ) ) = requestedSize;
    }
}

void raw_byte_buffer::resize_with_exponential_growth( size_type requestedSize ) {
    const size_type currentSize = size();

    if( currentSize < requestedSize ) {
        // The requested size will grow the buffer.  Call add_element to increase the allocated capacity exponentially.
        add_element( requestedSize - currentSize, 1 );
    } else if( currentSize > requestedSize ) {
        // The requested size will shrink the buffer.
        resize( requestedSize );
    }
}

// Make sure that the buffer has at least the requested capacity
void raw_byte_buffer::reserve( size_type requestedCapacity ) {
    if( m_data ) {
        size_type currentCapacity = *reinterpret_cast<const size_type*>( m_data + sizeof( size_type ) );
        if( currentCapacity < requestedCapacity ) {
            char* newData = (char*)realloc( m_data, 2 * sizeof( size_type ) + requestedCapacity );
            if( newData == 0 )
                throw std::runtime_error( "raw_byte_buffer.reserve: Failed to realloc " +
                                          lexical_cast<string>( 2 * sizeof( size_type ) + requestedCapacity ) +
                                          " bytes of memory." );
            m_data = newData;

            // Set the capacity member
            *reinterpret_cast<size_type*>( m_data + sizeof( size_type ) ) = requestedCapacity;
        }
    } else {
        m_data = (char*)malloc( 2 * sizeof( size_type ) + requestedCapacity );
        if( m_data == 0 )
            throw std::runtime_error( "raw_byte_buffer.reserve: Failed to malloc " +
                                      lexical_cast<string>( 2 * sizeof( size_type ) + requestedCapacity ) +
                                      " bytes of memory." );
        // Set the data size to 0
        *reinterpret_cast<std::size_t*>( m_data ) = 0;
        // Set the capacity to the value requested
        *reinterpret_cast<std::size_t*>( m_data + sizeof( size_type ) ) = requestedCapacity;
    }
}

void raw_byte_buffer::trim() {
    if( m_data ) {
        size_type currentSize = *reinterpret_cast<const size_type*>( m_data );
        size_type currentCapacity = *reinterpret_cast<const size_type*>( m_data + sizeof( size_type ) );

        if( currentCapacity > currentSize ) {
            char* newData = (char*)realloc( m_data, 2 * sizeof( size_type ) + currentSize );
            if( !newData )
                throw std::runtime_error( "raw_byte_buffer.trim: Failed to realloc " +
                                          lexical_cast<string>( 2 * sizeof( size_type ) + currentSize ) +
                                          " bytes of memory." );
            m_data = newData;

            *reinterpret_cast<size_type*>( m_data + sizeof( size_type ) ) = currentSize;
        }
    }
}

// Add an element of the requested size, returning a pointer to the newly allocated memory.
// The initialElementAllocation parameter specifies how much memory to allocate when adding the very first element.
raw_byte_buffer::iterator raw_byte_buffer::add_element( size_type elementSize, int initialElementAllocation ) {
    if( m_data ) {
        std::size_t dataSize = *reinterpret_cast<const size_type*>( m_data );
        std::size_t currentCapacity = *reinterpret_cast<const size_type*>( m_data + sizeof( size_type ) );
        std::size_t unusedSpace = currentCapacity - dataSize;
        // Check whether the new element fits within the existing buffer
        if( unusedSpace >= elementSize ) {
            // Increment the data size
            ( *reinterpret_cast<size_type*>( m_data ) ) += elementSize;
            // Return a pointer to what used to be the end of the buffer
            return m_data + 2 * sizeof( size_type ) + dataSize;
        }
        // It didn't fit, so we need to allocate a new buffer
        else {
            // Figure out how big a buffer to allocate by doubling the buffer size repeatedly
            std::size_t newCapacity = 2 * currentCapacity + 1;
            while( newCapacity < dataSize + elementSize )
                newCapacity *= 2;
            // TODO: If this allocation fails (and throws a bad alloc), should we try again with a smaller requested
            // size, to make things work closer to memory limits?
            char* newData = (char*)realloc( m_data, 2 * sizeof( size_type ) + newCapacity );
            if( newData == 0 )
                throw std::runtime_error( "raw_byte_buffer.add_element: Failed to realloc " +
                                          lexical_cast<string>( 2 * sizeof( size_type ) + newCapacity ) +
                                          " bytes of memory." );
            m_data = newData;

            // Set the size member
            ( *reinterpret_cast<size_type*>( m_data ) ) = dataSize + elementSize;
            // Set the capacity member
            ( *reinterpret_cast<size_type*>( m_data + sizeof( size_type ) ) ) = newCapacity;
            // Return a pointer to the newly added element
            return m_data + 2 * sizeof( size_type ) + dataSize;
        }
    } else {
        // Make the first allocation big enough for 'initialElementAllocation' structures of the requested size
        std::size_t newCapacity = initialElementAllocation * elementSize;
        m_data = (char*)malloc( 2 * sizeof( size_type ) + newCapacity );
        if( m_data == 0 )
            throw std::runtime_error( "raw_byte_buffer.add_element: Failed to malloc " +
                                      lexical_cast<string>( 2 * sizeof( size_type ) + newCapacity ) +
                                      " bytes of memory." );
        // Set the size member
        ( *reinterpret_cast<size_type*>( m_data ) ) = elementSize;
        // Set the capacity member
        ( *reinterpret_cast<size_type*>( m_data + sizeof( size_type ) ) ) = newCapacity;
        // Return a pointer to the newly added element
        return m_data + 2 * sizeof( size_type );
    }
}

} // namespace graphics
} // namespace frantic
