// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <limits>

namespace frantic {
namespace graphics {

/**
 * This class stores a raw buffer of bytes.  It may be desirable to move it outside of the graphics namespace
 * into a more general namespace at some point.
 *
 * The raw_byte_buffer is designed to be as low-impact as possible when it is zero-sized (it's just a null pointer),
 * and when it's allocated, the idea is to have just one buffer allocated and store everything in that buffer.  Having
 * two std::size_t elements at the beginning should keep alignment fine, because that's 8 or 16 bytes.  The buffer
 * also doesn't initialize its values to anything, so newly created values always need to be set to something.
 *
 */
class raw_byte_buffer {
    // m_data is structured as follows:
    //   std::size_t size;
    //   std::size_t capacity;
    //   char data[]
    char* m_data;

  public:
    // Typedefs required by the STL container concept
    typedef char* iterator;
    typedef const char* const_iterator;

    typedef char& reference;
    typedef const char& const_reference;

    typedef char* pointer;
    typedef std::ptrdiff_t difference_type;
    typedef std::size_t size_type;

    raw_byte_buffer()
        : m_data( 0 ) {}

    raw_byte_buffer( const raw_byte_buffer& rhs );

    /**
     *  Create a new buffer containing size bytes copied from data.
     *
     * @param data the data to copy.
     * @param size number of bytes to copy from data.
     */
    raw_byte_buffer( const void* data, std::size_t size );

    ~raw_byte_buffer() {
        if( m_data )
            free( m_data );
    }

    raw_byte_buffer& operator=( const raw_byte_buffer& rhs );

    /////////////
    // Basic container information
    /////////////

    size_type size() const {
        if( m_data )
            return *reinterpret_cast<const size_type*>( m_data );
        else
            return 0;
    }

    size_type capacity() const {
        if( m_data )
            return *reinterpret_cast<const size_type*>( m_data + sizeof( size_type ) );
        else
            return 0;
    }

    size_type max_size() const { return ( std::numeric_limits<size_type>::max )() - 2 * sizeof( size_type ); }

    bool empty() const {
        if( m_data )
            return *reinterpret_cast<const size_type*>( m_data ) == 0;
        else
            return true;
    }

    bool operator==( const raw_byte_buffer& rhs ) const;

    bool operator!=( const raw_byte_buffer& rhs ) const { return !operator==( rhs ); }

    /////////////
    // Iterator access through char* pointers
    /////////////

    iterator begin() {
        if( m_data )
            return m_data + 2 * sizeof( size_type );
        else
            return 0;
    }

    iterator end() {
        if( m_data )
            return m_data + 2 * sizeof( size_type ) + *reinterpret_cast<const size_type*>( m_data );
        else
            return 0;
    }

    const_iterator begin() const {
        if( m_data )
            return m_data + 2 * sizeof( size_type );
        else
            return 0;
    }

    const_iterator end() const {
        if( m_data )
            return m_data + 2 * sizeof( size_type ) + *reinterpret_cast<const size_type*>( m_data );
        else
            return 0;
    }

    /////////////
    // Access to pointers within the buffer
    /////////////

    iterator ptr_at( std::size_t offset ) {
        if( m_data )
            return m_data + 2 * sizeof( size_type ) + offset;
        else
            return 0;
    }

    const_iterator ptr_at( std::size_t offset ) const {
        if( m_data )
            return m_data + 2 * sizeof( size_type ) + offset;
        else
            return 0;
    }

    /////////////
    // Functions to modify the buffer
    /////////////

    // Clears the array by deallocating the buffer.  This returns the buffer to 0-sized, using only the space of one
    // pointer.
    void clear() {
        if( m_data ) {
            free( m_data );
            m_data = 0;
        }
    }

    char* release() {
        char* ret = m_data;
        m_data = 0;
        return ret;
    }

    void swap( raw_byte_buffer& rhs ) { std::swap( m_data, rhs.m_data ); }

    // Resize the buffer to the requested size
    void resize( size_type requestedSize );

    /**
     *  Resize the buffer to the requested size.
     *
     *  If the capacity must increase to fit the requested size,
     * then the capacity will increase exponentially.
     */
    void resize_with_exponential_growth( size_type requestedSize );

    // Make sure that the buffer has at least the requested capacity
    void reserve( size_type requestedCapacity );

    // Will release any memory that has been reserved, but is not currently used.
    void trim();

    /**
     * Add an element of the requested size, returning a pointer to the newly allocated memory.
     * The initialElementAllocation parameter specifies how much memory (in multiples of elementSize)
     * to allocate when adding the very first element.
     *
     * @param  elementSize  The size of the element to add.
     * @param  initialElementAllocation  When starting from an empty raw_byte_buffer, the number of elements to allocate
     * initially.
     */
    iterator add_element( size_type elementSize, int initialElementAllocation = 6 );

    // Add an element of the requested size, copying the provided data
    // The initialElementAllocation parameter specifies how much memory to allocate when adding the very first element.
    void add_element( const void* rawElementData, size_type elementSize, int initialElementAllocation = 6 ) {
        memcpy( add_element( elementSize, initialElementAllocation ), rawElementData, elementSize );
    }
};

} // namespace graphics
} // namespace frantic
