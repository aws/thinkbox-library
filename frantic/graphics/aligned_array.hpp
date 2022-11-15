// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdlib>
#include <tbb/tbb_stddef.h>
//#include <stdexcept>

namespace frantic {
namespace graphics {

template <class T>
class aligned_blocked_range {
  private:
    T m_begin;
    T m_end;
    size_t m_grainsize;
    size_t m_tokenSize;

  public:
    typedef size_t size_type;
    typedef T const_iterator;

    aligned_blocked_range( T begin, T end, size_t grainsize, size_t minTokenSize = 2 ) {
        m_begin = begin;
        m_end = end;
        m_grainsize = grainsize;
        m_tokenSize = minTokenSize;

        // std::cout << m_begin << ", " << m_end << std::endl;
    }

    aligned_blocked_range( aligned_blocked_range& rhs, tbb::split ) {
        // Partitions range into two subranges. The newly constructed blocked_range is
        // approximately the second half of the original range, and range is updated to
        // remainder. Each subrange has the same grainsize as the original range.

        // std::cout << "." <<  endl;

        if( !rhs.is_divisible() )
            throw std::runtime_error( "uh oh tbb...." );

        size_t start = rhs.size() / 2;

        m_begin = rhs.m_begin + ( start - start % rhs.m_tokenSize );
        m_end = rhs.m_end;

        m_grainsize = rhs.m_grainsize;
        m_tokenSize = rhs.m_tokenSize;

        rhs.m_end = m_begin;
    }

    size_type size() const { return m_end - m_begin; }

    bool empty() const { return size() <= 0; }

    size_type grainsize() const { return m_grainsize; }

    bool is_divisible() const {
        // std::cout << ": " <<  (size() > m_grainsize)  << endl;
        return ( size() > m_grainsize );
    }

    // iterators
    const_iterator begin() const { return m_begin; }
    const_iterator end() const { return m_end; }
};

template <class T>
class aligned_array {
  private:
    T* m_dataPtr;
    size_t m_size;

    aligned_array& operator=( const aligned_array& rhs );

    static const std::size_t ALIGNMENT = 16;

    inline static T* alloc( std::size_t count ) {
#ifdef _WIN32
        return reinterpret_cast<T*>( _aligned_malloc( size * sizeof( T ), ALIGNMENT ) );
#else
        void* result = NULL;
        if( 0 != posix_memalign( &result, ALIGNMENT, count * sizeof( T ) ) )
            throw std::bad_alloc();
        return reinterpret_cast<T*>( result );
#endif
    }

    inline static void alloc_free( T* ptr ) {
#ifdef _WIN32
        _aligned_free( ptr );
#else
        free( ptr );
#endif
    }

  public:
    aligned_array()
        : m_size( 0 )
        , m_dataPtr( 0 ) {}

    aligned_array( size_t size )
        : m_size( size ) {
        m_dataPtr = alloc( size );

        if( !m_dataPtr )
            throw std::runtime_error( "Aligned Memory allocation failed when trying to allocate " +
                                      boost::lexical_cast<std::string>( size * sizeof( T ) ) + " bytes " );

        m_size = size;
    }

    aligned_array( const aligned_array& rhs ) {
        if( rhs.size() != 0 ) {
            m_dataPtr = alloc( rhs.m_size );

            if( !m_dataPtr )
                throw std::runtime_error( "Aligned Memory allocation failed when trying to allocate " +
                                          boost::lexical_cast<std::string>( rhs.m_size * sizeof( T ) ) + " bytes " );

            m_size = rhs.size();

            memcpy( m_dataPtr, rhs.m_dataPtr, rhs.m_size * sizeof( T ) );
        } else {
            m_dataPtr = rhs.m_dataPtr;
            m_size = rhs.m_size;
        }
    }

    ~aligned_array() {
        //_aligned_free( m_dataPtr );
        alloc_free( m_dataPtr );
    }

    T& operator[]( size_t i ) { return *( m_dataPtr + i ); }

    const T& operator[]( size_t i ) const { return *( m_dataPtr + i ); }

    void resize( size_t size ) {
        // do nothing if the array is already big enough
        if( m_size == size )
            return;

        T* temp = alloc( size );

        if( !temp )
            throw std::runtime_error( "Aligned Memory allocation failed when trying to allocate " +
                                      boost::lexical_cast<std::string>( size * sizeof( T ) ) + " bytes " );

        if( m_dataPtr ) {
            // copy any existing data
            memcpy( temp, m_dataPtr, m_size * sizeof( T ) );
            // free the previous block
            //_aligned_free( m_dataPtr );
            alloc_free( m_dataPtr );
        }

        m_dataPtr = temp;
        m_size = size;
    }

    void clear() {
        if( m_dataPtr ) {
            // free the previous block
            //_aligned_free( m_dataPtr );
            alloc_free( m_dataPtr );
            m_dataPtr = 0;
            m_size = 0;
        }
    }

    void swap( aligned_array<T>& rhs ) {
        T* temp = rhs.m_dataPtr;
        size_t size = rhs.m_size;

        rhs.m_dataPtr = m_dataPtr;
        rhs.m_size = m_size;

        m_dataPtr = temp;
        m_size = size;
    }

    size_t size() const { return m_size; }

    const T* get() const { return m_dataPtr; }

    T* get() { return m_dataPtr; }
};

} // namespace graphics
} // namespace frantic
