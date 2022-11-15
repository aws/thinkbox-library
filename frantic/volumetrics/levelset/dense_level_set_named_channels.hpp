// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

using frantic::channels::data_type_t;
using frantic::graphics::raw_byte_buffer;

// TODO: Much of this stuff is shared with the trimesh3 named channels.  We'll probably want to refactor the shared code
// into a common
//       base class once this has settled down a bit.

// Forward declaration so we can be friends with this class.
class dense_level_set;

///////////////////////
// CHANNEL GENERAL ACCESSOR BASE
///////////////////////

// This class is used to provide faster access within a named channel.
template <class DataVector>
class dense_level_set_channel_general_accessor_base {
  protected:
    DataVector* m_data;

    std::size_t m_arity;         // # of data values for this channel per defined voxel
    data_type_t m_dataType;      // primitive data type
    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    // Function which knows how to do weighted sum combinations of this data type
    frantic::channels::channel_weighted_sum_combine_function_t m_weightedSumFunction;

    dense_level_set_channel_general_accessor_base( DataVector* data, std::size_t arity, data_type_t dataType )
        : m_data( data )
        , m_arity( arity )
        , m_dataType( dataType ) {
        m_weightedSumFunction = channel_weighted_sum_combine_function( dataType );
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    friend class dense_level_set;

  public:
    std::size_t size() const { return m_data ? ( m_data->size() / m_primitiveSize ) : 0; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    std::string type_str() const { return channel_data_type_str( m_arity, m_dataType ); }

    bool valid() const { return m_data != 0; }

    operator bool() const { return m_data != 0; }

    /**
     * This returns the size of the primitive of this channel.  That's equivalent to arity() * (data type size).
     */
    std::size_t primitive_size() const { return m_primitiveSize; }

    ///////////
    // Access to the data
    ///////////

    const char* data( std::size_t i ) const { return m_data->ptr_at( i * m_primitiveSize ); }

    // m_inputChannels[i].get_linear_interpolated( dataIndex0, dataIndex1, alpha, outputChannels[i].add_vertex() );
    void get_linear_interpolated( int i0, int i1, float alpha, char* outResult ) const {
        float weightsArray[2];
        weightsArray[0] = 1 - alpha;
        weightsArray[1] = alpha;
        const char* dataArray[2];
        dataArray[0] = m_data->ptr_at( i0 * m_primitiveSize );
        dataArray[1] = m_data->ptr_at( i1 * m_primitiveSize );

        m_weightedSumFunction( weightsArray, dataArray, 2, m_arity, outResult );
    }
};

///////////////////////
// CHANNEL GENERAL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
class const_dense_level_set_channel_general_accessor;

class dense_level_set_channel_general_accessor : public dense_level_set_channel_general_accessor_base<raw_byte_buffer> {
  protected:
    dense_level_set_channel_general_accessor()
        : dense_level_set_channel_general_accessor_base<raw_byte_buffer>( 0, 0, frantic::channels::data_type_int8 ) {}
    dense_level_set_channel_general_accessor( raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : dense_level_set_channel_general_accessor_base<raw_byte_buffer>( data, arity, dataType ) {}

    // This function adds a new element to the end of the array, and returns a pointer to its memory.
    // NOTE: Ordinary code shouldn't access this, but the dense_level_set itself should.
    char* add_element() { return m_data->add_element( m_primitiveSize ); }

    friend class const_dense_level_set_channel_general_accessor;
    friend class dense_level_set_channel;
    friend class dense_level_set;

  public:
    ///////////
    // Access to the data
    ///////////

    char* data( std::size_t i ) { return m_data->ptr_at( i * m_primitiveSize ); }
};

///////////////////////
// CONST CHANNEL GENERAL ACCESSOR
///////////////////////

class const_dense_level_set_channel_general_accessor
    : public dense_level_set_channel_general_accessor_base<const raw_byte_buffer> {
  protected:
    const_dense_level_set_channel_general_accessor()
        : dense_level_set_channel_general_accessor_base<const raw_byte_buffer>( 0, 0,
                                                                                frantic::channels::data_type_int8 ) {}
    const_dense_level_set_channel_general_accessor( const raw_byte_buffer* data, std::size_t arity,
                                                    data_type_t dataType )
        : dense_level_set_channel_general_accessor_base<const raw_byte_buffer>( data, arity, dataType ) {}

    friend class dense_level_set_channel;

  public:
    const_dense_level_set_channel_general_accessor( const dense_level_set_channel_general_accessor& rhs )
        : dense_level_set_channel_general_accessor_base<const raw_byte_buffer>( rhs.m_data, rhs.m_arity,
                                                                                rhs.m_dataType ) {}
};

///////////////////////
// CHANNEL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
template <class DataType>
class const_dense_level_set_channel_accessor;

// The non-const accessor adds additional methods for writing to the channel
template <class DataType>
class dense_level_set_channel_accessor : public dense_level_set_channel_general_accessor {

    dense_level_set_channel_accessor( raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : dense_level_set_channel_general_accessor( data, arity, dataType ) {}

    friend class const_dense_level_set_channel_accessor<DataType>;
    friend class dense_level_set_channel;
    friend class dense_level_set;

  public:
    dense_level_set_channel_accessor()
        : dense_level_set_channel_general_accessor( 0, 0, frantic::channels::data_type_int8 ) {}

    ///////////
    // Access to the data
    ///////////

    const DataType& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }

    DataType& operator[]( std::size_t i ) {
        return *reinterpret_cast<DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }

    void add_element( const DataType& element ) {
        *reinterpret_cast<DataType*>( dense_level_set_channel_general_accessor::add_element() ) = element;
    }
};

///////////////////////
// CONST CHANNEL ACCESSOR
///////////////////////

// The const accessor just uses all the methods in the base accessor class
template <class DataType>
class const_dense_level_set_channel_accessor : public const_dense_level_set_channel_general_accessor {

    const_dense_level_set_channel_accessor( const raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : const_dense_level_set_channel_general_accessor( data, arity, dataType ) {}

    friend class dense_level_set_channel;

  public:
    const_dense_level_set_channel_accessor()
        : const_dense_level_set_channel_general_accessor( 0, 0, frantic::channels::data_type_int8 ) {}

    const_dense_level_set_channel_accessor( const dense_level_set_channel_accessor<DataType>& rhs )
        : const_dense_level_set_channel_general_accessor( rhs.m_data, rhs.m_arity, rhs.m_dataType ) {}
    ///////////
    // Access to the data
    ///////////

    const DataType& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }
};

class dense_level_set_channel {
    frantic::tstring m_name; // string name of the channel
    std::size_t m_arity;     // # of data values for this channel per defined voxel
    data_type_t m_dataType;  // primitive data type

    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    raw_byte_buffer m_data;

    dense_level_set_channel( const frantic::tstring& name, std::size_t arity, data_type_t dataType )
        : m_name( name )
        , m_arity( arity )
        , m_dataType( dataType ) {
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    void set( const frantic::tstring& name, std::size_t arity, data_type_t dataType ) {
        m_name = name;
        m_arity = arity;
        m_dataType = dataType;
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    friend class dense_level_set;

  public:
    const frantic::tstring& name() const { return m_name; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    std::string type_str() const {
        return frantic::strings::to_string( channel_data_type_str( m_arity, m_dataType ) ); // HACK: Remove to_string
    }

    /**
     * This returns the size of the primitive of this channel.  That's equivalent to arity() * (data type size).
     */
    std::size_t primitive_size() const { return m_primitiveSize; }

    std::size_t size() const { return m_data.size() / m_primitiveSize; }

    template <typename T>
    bool supports_type() const {
        return channels::channel_data_type_traits<T>::arity() == m_arity &&
               channels::channel_data_type_traits<T>::data_type == m_dataType;
    }

    ///////////
    // Access to the data
    ///////////

    void* data() { return m_data.begin(); }

    const void* data() const { return m_data.begin(); }

    dense_level_set_channel_general_accessor get_general_accessor() {
        return dense_level_set_channel_general_accessor( &m_data, m_arity, m_dataType );
    }

    const_dense_level_set_channel_general_accessor get_general_accessor() const {
        return const_dense_level_set_channel_general_accessor( &m_data, m_arity, m_dataType );
    }

    // This returns an accessor to this channel, viewed as the requested type.
    template <class T>
    dense_level_set_channel_accessor<T> get_accessor() {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            return dense_level_set_channel_accessor<T>( &m_data, m_arity, m_dataType );
        } else {
            throw std::runtime_error(
                "dense_level_set_channel.get_accessor: Could not access channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                " is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }

    // This returns a const accessor to this channel, viewed as the requested type.
    template <class T>
    const_dense_level_set_channel_accessor<T> get_accessor() const {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            return const_dense_level_set_channel_accessor<T>( &m_data, m_arity, m_dataType );
        } else {
            throw std::runtime_error(
                "dense_level_set_channel.get_accessor: Could not access channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                " is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
