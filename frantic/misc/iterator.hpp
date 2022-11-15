// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file
 * @author Stephen Kiazyk
 *
 * Some fun iterator objects I've found useful sometimes.
 */

#pragma once

#include <iterator>

namespace frantic {

/**
 * An iterator for when you need to know the _size_ of an output operation, but dealing
 * with buffer allocation is unnessary/inefficient.
 *
 * Typically, the use of this object with standard algorithm methods would be:
 *
 * count_output_iterator counter = std::some_algorithm( obj.begin(), obj.end(), count_output_iterator() );
 */
class count_output_iterator {
  public:
    typedef std::output_iterator_tag iterator_category;
    typedef void value_type;
    typedef void difference_type;
    typedef void pointer;
    typedef void reference;

  private:
    size_t m_count;

  public:
    explicit count_output_iterator()
        : m_count( 0 ) {}
    template <class T>
    count_output_iterator& operator=( const T& /*v*/ ) {
        ++m_count;
        return *this;
    }
    count_output_iterator& operator*() { return *this; }
    count_output_iterator& operator++() { return *this; }
    count_output_iterator& operator++( int ) { return *this; }
    size_t count() const { return m_count; }
    void reset() { m_count = 0; }
};

/**
 * An iterator for when the expected output size of an operation is bounded and/or
 * ignoring extra output values beyond the first few is acceptable (especially if it
 * avoids a dynamic allocation).
 *
 * Typically one would expect to use this with a static array, for example.
 * float array[16];
 * limited_output_iterator<float*> output = some_generator( limited_output_iterator<float*>( array, 16 ) );
 */
template <class OutputIterator>
class limited_output_iterator {
  public:
    typedef std::output_iterator_tag iterator_category;
    typedef void value_type;
    typedef void difference_type;
    typedef void pointer;
    typedef void reference;

  private:
    OutputIterator m_current;
    size_t m_count;
    size_t m_limit;

  public:
    explicit limited_output_iterator( OutputIterator current, size_t limit )
        : m_current( current )
        , m_count( 0 )
        , m_limit( limit ) {}

    limited_output_iterator& operator=( const typename std::iterator_traits<OutputIterator>::value_type& v ) {
        if( m_count < m_limit ) {
            *m_current = v;
            ++m_current;
        }
        ++m_count;
        return *this;
    }
    limited_output_iterator& operator*() { return *this; }
    limited_output_iterator& operator++() { return *this; }
    limited_output_iterator& operator++( int ) { return *this; }
    size_t count() const { return m_count; }
    size_t limit() const { return m_limit; }
    OutputIterator current() const { return m_current; }
};

/**
 * A method to avoid needing explicit template instanciation when creating a limited output iterator
 */
template <class OutputIterator>
limited_output_iterator<OutputIterator> make_limited_output_iterator( OutputIterator it, size_t limit ) {
    return limited_output_iterator<OutputIterator>( it, limit );
}

} // namespace frantic
