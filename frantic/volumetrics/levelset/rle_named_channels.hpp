// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>

namespace frantic {
namespace fluids {
// Forward declaration so we can be friends with this class.
class rle_voxel_field;
class rle_staggered_vel_accessor;
} // namespace fluids
} // namespace frantic

namespace frantic {
namespace volumetrics {
namespace levelset {

// This should never be done, importing a namespace into another is a bad idea.
using frantic::channels::data_type_t;
using frantic::graphics::raw_byte_buffer;

// TODO: Much of this stuff is shared with the trimesh3 named channels.  We'll probably want to refactor the shared code
// into a common
//       base class once this has settled down a bit.

// Forward declaration so we can be friends with this class.
class rle_level_set;

///////////////////////
// CHANNEL GENERAL ACCESSOR BASE
///////////////////////

// This class is used to provide faster access within a named channel.
template <class DataVector>
class rle_channel_general_accessor_base {
  protected:
    DataVector* m_data;

    std::size_t m_arity;         // # of data values for this channel per defined voxel
    data_type_t m_dataType;      // primitive data type
    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    // Function which knows how to do weighted sum combinations of this data type
    frantic::channels::channel_weighted_sum_combine_function_t m_weightedSumFunction;

    rle_channel_general_accessor_base( DataVector* data, std::size_t arity, data_type_t dataType )
        : m_data( data )
        , m_arity( arity )
        , m_dataType( dataType ) {
        m_weightedSumFunction = channel_weighted_sum_combine_function( dataType );
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    friend class rle_level_set;
    friend class frantic::fluids::rle_voxel_field;

  public:
    std::size_t size() const { return m_data ? ( m_data->size() / m_primitiveSize ) : 0; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    std::string type_str() const {
        return frantic::strings::to_string(
            channel_data_type_str( m_arity, m_dataType ) ); // HACK: Remove frantic::strings::to_string(
    }

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

    frantic::channels::channel_weighted_sum_combine_function_t get_weighted_sum_combine_function() const {
        return m_weightedSumFunction;
    }

    // trilerp by precalculated values.
    // TODO: consider calculating weights and passing them in.
    //		 This would save calculations when using the same weights multiple times
    void get_trilinear_interpolated( const float* multipliers, const boost::int32_t* trilerpDataIndices,
                                     char* out ) const {
        const char* dataArray[8];
        char* undefined;

        // Set undefined to all zeros.
        // TODO: Renormalize weights instead of using all zero's or
        // find a better way to choose this value;
        undefined = (char*)alloca( m_primitiveSize );
        memset( undefined, 0, m_primitiveSize );

        // Fill data array, setting any out of bound point to 0
        for( int i = 0; i < 8; ++i ) {
            if( trilerpDataIndices[i] < 0 )
                dataArray[i] = undefined;
            else
                dataArray[i] = data( trilerpDataIndices[i] );
        }

        m_weightedSumFunction( multipliers, dataArray, 8, m_arity, out );
    }
};

///////////////////////
// CHANNEL GENERAL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
class const_rle_channel_general_accessor;

class rle_channel_general_accessor : public rle_channel_general_accessor_base<raw_byte_buffer> {
  protected:
    rle_channel_general_accessor( raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : rle_channel_general_accessor_base<raw_byte_buffer>( data, arity, dataType ) {}

    // This function adds a new element to the end of the array, and returns a pointer to its memory.
    // NOTE: Ordinary code shouldn't access this, but the rle_level_set itself should.
    char* add_element() { return m_data->add_element( m_primitiveSize ); }

    friend class const_rle_channel_general_accessor;
    friend class rle_channel;
    friend class rle_level_set;
    friend class frantic::fluids::rle_voxel_field;

    // Some friend functions
    friend void detail::convert_intersections_to_level_set(
        const std::vector<std::vector<geometry::scan_conversion_intersection>>& intersectionDepths,
        const geometry::trimesh3& geometry, float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside,
        bool isHalfOpenSurface, int exteriorRegionCodeNeg, int exteriorRegionCodePos,
        const frantic::graphics::boundbox3& resultVoxelBounds, rle_level_set& outResult );
    friend void detail::finalize_geometry_to_levelset_conversion( int exteriorRegionCode, rle_level_set& levelSet );
    friend void detail::intermediate_mesh_to_level_set_union( rle_level_set& levelSet,
                                                              const rle_level_set& inputLevelSet );
    friend void detail::mesh_to_level_set_region_x_dilation( rle_level_set& levelSet );

  public:
    rle_channel_general_accessor()
        : rle_channel_general_accessor_base<raw_byte_buffer>( 0, 0, frantic::channels::data_type_int8 ) {}

    ///////////
    // Access to the data
    ///////////

    char* data( std::size_t i ) { return m_data->ptr_at( i * m_primitiveSize ); }

    const char* data( std::size_t i ) const { return m_data->ptr_at( i * m_primitiveSize ); }
};

///////////////////////
// CONST CHANNEL GENERAL ACCESSOR
///////////////////////

class const_rle_channel_general_accessor : public rle_channel_general_accessor_base<const raw_byte_buffer> {
  protected:
    const_rle_channel_general_accessor( const raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : rle_channel_general_accessor_base<const raw_byte_buffer>( data, arity, dataType ) {}

    friend class rle_channel;

  public:
    const_rle_channel_general_accessor()
        : rle_channel_general_accessor_base<const raw_byte_buffer>( 0, 0, frantic::channels::data_type_int8 ) {}
    const_rle_channel_general_accessor( const rle_channel_general_accessor& rhs )
        : rle_channel_general_accessor_base<const raw_byte_buffer>( rhs.m_data, rhs.m_arity, rhs.m_dataType ) {}
};

///////////////////////
// CHANNEL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
template <class DataType>
class const_rle_channel_accessor;

// The non-const accessor adds additional methods for writing to the channel
template <class DataType>
class rle_channel_accessor : public rle_channel_general_accessor {

    rle_channel_accessor( raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : rle_channel_general_accessor( data, arity, dataType ) {}

    friend class const_rle_channel_accessor<DataType>;
    friend class rle_channel;
    friend class rle_level_set;
    friend class frantic::fluids::rle_voxel_field;
    friend class frantic::fluids::rle_staggered_vel_accessor;

  public:
    rle_channel_accessor()
        : rle_channel_general_accessor( 0, 0, frantic::channels::data_type_int8 ) {}

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
        *reinterpret_cast<DataType*>( rle_channel_general_accessor::add_element() ) = element;
    }
};

///////////////////////
// CONST CHANNEL ACCESSOR
///////////////////////

// The const accessor just uses all the methods in the base accessor class
template <class DataType>
class const_rle_channel_accessor : public const_rle_channel_general_accessor {

    const_rle_channel_accessor( const raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : const_rle_channel_general_accessor( data, arity, dataType ) {}

    friend class rle_channel;

  public:
    const_rle_channel_accessor()
        : const_rle_channel_general_accessor( 0, 0, frantic::channels::data_type_int8 ) {}

    const_rle_channel_accessor( const rle_channel_accessor<DataType>& rhs )
        : const_rle_channel_general_accessor( rhs.m_data, rhs.m_arity, rhs.m_dataType ) {}

    ///////////
    // Access to the data
    ///////////

    const DataType& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }
};

class rle_channel {
    frantic::tstring m_name; // string name of the channel
    std::size_t m_arity;     // # of data values for this channel per defined voxel
    data_type_t m_dataType;  // primitive data type

    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    raw_byte_buffer m_data;

    rle_channel( const frantic::tstring& name, std::size_t arity, data_type_t dataType )
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

    friend class rle_level_set;
    friend class frantic::fluids::rle_voxel_field;

  public:
    const frantic::tstring& name() const { return m_name; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    frantic::tstring type_str() const { return channel_data_type_str( m_arity, m_dataType ); }

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

    char* data() { return m_data.begin(); }

    const char* data() const { return m_data.begin(); }

    rle_channel_general_accessor get_general_accessor() {
        return rle_channel_general_accessor( &m_data, m_arity, m_dataType );
    }

    const_rle_channel_general_accessor get_general_accessor() const {
        return const_rle_channel_general_accessor( &m_data, m_arity, m_dataType );
    }

    // This returns an accessor to this channel, viewed as the requested type.
    template <class T>
    rle_channel_accessor<T> get_accessor() {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            return rle_channel_accessor<T>( &m_data, m_arity, m_dataType );
        } else {
            throw std::runtime_error(
                "rle_channel.get_accessor: Could not access channel \"" + frantic::strings::to_string( m_name ) +
                "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                " is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }

    // This returns a const accessor to this channel, viewed as the requested type.
    template <class T>
    const_rle_channel_accessor<T> get_accessor() const {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            return const_rle_channel_accessor<T>( &m_data, m_arity, m_dataType );
        } else {
            throw std::runtime_error(
                "rle_channel.get_accessor: Could not access channel \"" + frantic::strings::to_string( m_name ) +
                "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                " is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }

    bool check_consistency( std::ostream& out ) const {
        if( m_primitiveSize != channels::sizeof_channel_data_type( m_dataType ) * m_arity ) {
            out << "Channel " << frantic::strings::to_string( m_name ) << " has primitive size " << m_primitiveSize
                << " which does not match its data type " << m_dataType << " and arity " << m_arity << std::endl;
            return false;
        }

        if( !channels::is_valid_channel_name( m_name ) ) {
            out << "Channel name is not valid: \"" << frantic::strings::to_string( m_name ) << "\"" << std::endl;
            return false;
        }

        if( m_data.capacity() < m_data.size() ) {
            out << "Channel " << frantic::strings::to_string( m_name ) << " has data capacity " << m_data.capacity()
                << " that is less then the channel data size " << m_data.size() << std::endl;
            return false;
        }

        return true;
    }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
