// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

#include <boost/call_traits.hpp>
#include <boost/shared_ptr.hpp>

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace channels {

///////////////////////
// Provides native access to a channel in a channel_map structure. Only use this if
// you are working with channel_maps you have defined and can be sure of the types.
///////////////////////
template <class T>
class channel_accessor {
    std::size_t m_offset;

  public:
    // Default constructor makes it point to an invalid position.  Putting in an
    // "invalid" state would be more overhead than we want.
    channel_accessor()
        : m_offset( ( std::numeric_limits<std::size_t>::max )() ) {}

    channel_accessor( const channel_accessor& accessor )
        : m_offset( accessor.m_offset ) {}

    channel_accessor( std::size_t offset )
        : m_offset( offset ) {}

    bool is_valid() const { return m_offset != ( std::numeric_limits<std::size_t>::max )(); }
    void reset() { m_offset = ( std::numeric_limits<std::size_t>::max )(); }

    /**
     * Retrieves a read/write reference to the data within the provided particle for the accessor's channel.
     */
    T& get( char* rawParticleData ) const {
        return *reinterpret_cast<T*>( reinterpret_cast<char*>( rawParticleData ) + m_offset );
    }

    /**
     * Retrieves a read-only reference to the data within the provided particle for the accessor's channel.
     */
    const T& get( const char* rawParticleData ) const {
        return *reinterpret_cast<const T*>( reinterpret_cast<const char*>( rawParticleData ) + m_offset );
    }

    /**
     * Retrieves a read/write reference to the data within the provided particle for the accessor's channel.
     */
    T& get( std::vector<char>& rawParticleData ) const {
        return *reinterpret_cast<T*>( &rawParticleData[0] + m_offset );
    }

    /**
     * Retrieves a read-only reference to the data within the provided particle for the accessor's channel.
     */
    const T& get( const std::vector<char>& rawParticleData ) const {
        return *reinterpret_cast<const T*>( &rawParticleData[0] + m_offset );
    }

    T& operator()( char* rawParticleData ) const {
        return *reinterpret_cast<T*>( reinterpret_cast<char*>( rawParticleData ) + m_offset );
    }

    const T& operator()( const char* rawParticleData ) const {
        return *reinterpret_cast<const T*>( reinterpret_cast<const char*>( rawParticleData ) + m_offset );
    }

    // Convenience function provided for vector<char>
    T& operator()( std::vector<char>& rawParticleData ) const {
        return *reinterpret_cast<T*>( &rawParticleData[0] + m_offset );
    }

    // Convenience function provided for vector<char>
    const T& operator()( const std::vector<char>& rawParticleData ) const {
        return *reinterpret_cast<const T*>( &rawParticleData[0] + m_offset );
    }
};

namespace detail {
class string_accessor_helper {
    char* m_rawChannelStringData;

    friend class channel_accessor<frantic::tstring>;

    string_accessor_helper( char* rawChannelStringData, int )
        : m_rawChannelStringData( rawChannelStringData ) {}

    // non-copyable (for now, until it becomes necessary to copy for some reason)
    string_accessor_helper& operator=( const string_accessor_helper& rhs );

  public:
    string_accessor_helper( const string_accessor_helper& rhs ) { m_rawChannelStringData = rhs.m_rawChannelStringData; }

    string_accessor_helper& operator=( const frantic::tstring& rhs ) {
        assign_channel_string( m_rawChannelStringData, rhs );
        return *this;
    }

    operator frantic::tstring() const { return string_from_channel_string( m_rawChannelStringData ); }
};

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator==( const string_accessor_helper& lhs, const frantic::tstring& rhs ) {
    return frantic::tstring( lhs ) == rhs;
}

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator==( const frantic::tstring& lhs, const string_accessor_helper& rhs ) {
    return lhs == frantic::tstring( rhs );
}

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator!=( const string_accessor_helper& lhs, const frantic::tstring& rhs ) {
    return frantic::tstring( lhs ) != rhs;
}

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator!=( const frantic::tstring& lhs, const string_accessor_helper& rhs ) {
    return lhs != frantic::tstring( rhs );
}

class const_string_accessor_helper {
    const char* m_rawChannelStringData;

    friend class channel_accessor<frantic::tstring>;

    const_string_accessor_helper( const char* rawChannelStringData, int )
        : m_rawChannelStringData( rawChannelStringData ) {}

    // non-copyable (for now, until it becomes necessary to copy for some reason)
    const_string_accessor_helper& operator=( const const_string_accessor_helper& rhs );
    // No assignment from a string, because it's const.
    const_string_accessor_helper& operator=( const frantic::tstring& rhs );

  public:
    const_string_accessor_helper( const const_string_accessor_helper& rhs ) {
        m_rawChannelStringData = rhs.m_rawChannelStringData;
    }

    operator frantic::tstring() const { return string_from_channel_string( m_rawChannelStringData ); }
};

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator==( const const_string_accessor_helper& lhs, const frantic::tstring& rhs ) {
    return frantic::tstring( lhs ) == rhs;
}

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator==( const frantic::tstring& lhs, const const_string_accessor_helper& rhs ) {
    return lhs == frantic::tstring( rhs );
}

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator!=( const const_string_accessor_helper& lhs, const frantic::tstring& rhs ) {
    return frantic::tstring( lhs ) != rhs;
}

// Implement string comparison with C-style string for the string_accessor_helper
inline bool operator!=( const frantic::tstring& lhs, const const_string_accessor_helper& rhs ) {
    return lhs != frantic::tstring( rhs );
}

} // namespace detail

/**
 * Specialized channel_accessor for the data_type_string type.  It uses a special string_accessor_helper
 * to make assignment and usage of array access act similar to a std::string.
 */
template <>
class channel_accessor<frantic::tstring> {
    std::size_t m_offset;

  public:
    // Default constructor makes it point to and invalid position.  Putting in an
    // "invalid" state would be more overhead than we want.
    channel_accessor()
        : m_offset( ( std::numeric_limits<std::size_t>::max )() ) {}

    channel_accessor( std::size_t offset )
        : m_offset( offset ) {}

    bool is_valid() const { return m_offset != ( std::numeric_limits<std::size_t>::max )(); }
    void reset() { m_offset = ( std::numeric_limits<std::size_t>::max )(); }

    // Functions to use the accessor with a raw memory particle
    detail::string_accessor_helper get( char* rawParticleData ) {
        return detail::string_accessor_helper( rawParticleData + m_offset, 0 );
    }

    const detail::const_string_accessor_helper get( const char* rawParticleData ) const {
        return detail::const_string_accessor_helper( rawParticleData + m_offset, 0 );
    }

    // Convenience function provided for vector<char>
    detail::string_accessor_helper get( std::vector<char>& rawParticleData ) const {
        return detail::string_accessor_helper( &rawParticleData[0] + m_offset, 0 );
    }

    // Convenience function provided for vector<char>
    const detail::const_string_accessor_helper get( const std::vector<char>& rawParticleData ) const {
        return detail::const_string_accessor_helper( &rawParticleData[0] + m_offset, 0 );
    }

    detail::string_accessor_helper operator()( char* rawParticleData ) {
        return detail::string_accessor_helper( &rawParticleData[0] + m_offset, 0 );
    }

    const detail::const_string_accessor_helper operator()( const char* rawParticleData ) const {
        return detail::const_string_accessor_helper( rawParticleData + m_offset, 0 );
    }

    /**
     * Convenience function provided for vector<char>
     */
    detail::string_accessor_helper operator()( std::vector<char>& rawParticleData ) const {
        return detail::string_accessor_helper( &rawParticleData[0] + m_offset, 0 );
    }

    /**
     * Convenience function provided for vector<char>
     */
    const detail::const_string_accessor_helper operator()( const std::vector<char>& rawParticleData ) const {
        return detail::const_string_accessor_helper( &rawParticleData[0] + m_offset, 0 );
    }
};

///////////////////////
// Provides access to a particle's data member with automatic type conversion. For use
// when the actual data type is not known at compile time.
///////////////////////
template <class T>
class channel_cvt_accessor {
  private:
    T m_default;
    std::size_t m_offset;
    channels::channel_type_convertor_function_t m_convertGet;
    channels::channel_type_convertor_function_t m_convertSet;

  public:
    typedef T value_type;

  public:
    channel_cvt_accessor()
        : m_offset( ( std::numeric_limits<std::size_t>::max )() )
        , m_convertGet( NULL )
        , m_convertSet( NULL ) {}

    channel_cvt_accessor( typename boost::call_traits<T>::param_type defaultVal )
        : m_default( defaultVal )
        , m_offset( 0 )
        , m_convertGet( NULL )
        , m_convertSet( NULL ) {}

    channel_cvt_accessor( std::size_t offset, std::size_t sourceDataArity, data_type_t sourceDataType,
                          const frantic::tstring& channelNameForErrorMessage )
        : m_offset( offset ) {
        data_type_t targetDataType = channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity )
            throw std::runtime_error(
                std::string() + "cvt_channel_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );

        m_convertGet =
            get_channel_type_convertor_function( sourceDataType, targetDataType, channelNameForErrorMessage );
        m_convertSet =
            get_channel_type_convertor_function( targetDataType, sourceDataType, channelNameForErrorMessage );
    }

    bool is_valid() const { return m_offset != ( std::numeric_limits<std::size_t>::max )(); }

    bool is_default() const { return m_offset == 0 && m_convertSet == NULL; }

    void reset() {
        m_offset = ( std::numeric_limits<std::size_t>::max )();
        m_convertGet = m_convertSet = NULL;
    }

    void reset( typename boost::call_traits<T>::param_type defaultVal ) {
        m_offset = 0;
        m_convertGet = m_convertSet = NULL;
        m_default = defaultVal;
    }

    const T get( const char* rawParticleData ) const {
        if( !m_convertGet )
            return m_default;

        T result;
        m_convertGet( (char*)&result, rawParticleData + m_offset, channel_data_type_traits<T>::arity() );

        return result;
    }

    void set( char* rawParticleData, const T& value ) const {
        if( !m_convertSet )
            return;
        m_convertSet( rawParticleData + m_offset, (const char*)&value, channel_data_type_traits<T>::arity() );
    }

    // Convenience function provided for vector<char>
    const T get( const std::vector<char>& rawParticleData ) const { return get( &rawParticleData[0] ); }

    // Convenience function provided for alternate syntax
    const T operator()( const char* rawParticleData ) const { return get( rawParticleData ); }

    // Convenience function provided for alternate syntax and vector<char>
    const T operator()( const std::vector<char>& rawParticleData ) const { return get( &rawParticleData[0] ); }

    // Convenience function provided for vector<char>
    void set( std::vector<char>& rawParticleData, const T& value ) const { set( &rawParticleData[0], value ); }
};

///////////////////////
// Provides access to a particle's data member with automatic type conversion. For use
// when the actual data type is not known at compile time.  Const version so it only
// converts one way.
///////////////////////
template <class T>
class channel_const_cvt_accessor {
  private:
    T m_default;
    std::size_t m_offset;
    channels::channel_type_convertor_function_t m_convertGet;

  public:
    typedef T value_type;

  public:
    channel_const_cvt_accessor()
        : m_offset( ( std::numeric_limits<std::size_t>::max )() )
        , m_convertGet( NULL ) {}

    channel_const_cvt_accessor( typename boost::call_traits<T>::param_type defaultVal )
        : m_default( defaultVal )
        , m_offset( 0 )
        , m_convertGet( NULL ) {}

    channel_const_cvt_accessor( std::size_t offset, std::size_t sourceDataArity, data_type_t sourceDataType,
                                const frantic::tstring& channelNameForErrorMessage )
        : m_offset( offset ) {
        data_type_t targetDataType = channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity )
            throw std::runtime_error(
                std::string() + "cvt_const_channel_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );

        m_convertGet =
            get_channel_type_convertor_function( sourceDataType, targetDataType, channelNameForErrorMessage );
    }

    bool is_valid() const { return m_offset != ( std::numeric_limits<std::size_t>::max )(); }

    bool is_default() const { return is_valid() && m_convertGet == NULL; }

    void reset() {
        m_offset = ( std::numeric_limits<std::size_t>::max )();
        m_convertGet = NULL;
    }

    void reset( typename boost::call_traits<T>::param_type defaultVal ) {
        m_offset = 0;
        m_convertGet = NULL;
        m_default = defaultVal;
    }

    const T get( const char* rawParticleData ) const {
        if( !m_convertGet )
            return m_default;

        T result;
        m_convertGet( (char*)&result, rawParticleData + m_offset, channel_data_type_traits<T>::arity() );

        return result;
    }

    // Convenience function provided for vector<char>
    const T get( const std::vector<char>& rawParticleData ) const { return get( &rawParticleData[0] ); }

    // Convenience function provided for alternate syntax
    const T operator()( const char* rawParticleData ) const { return get( rawParticleData ); }

    // Convenience function provided for alternate syntax and vector<char>
    const T operator()( const std::vector<char>& rawParticleData ) const { return get( &rawParticleData[0] ); }
};

/**
 * Specialized channel_cvt_accessor for the data_type_string type.  It uses a special string_accessor_helper
 * to make assignment and usage of array access act similar to a std::string.
 *
 * @note  Because the string data type can't be implicitly converted to any other type, this is effectively
 *        the same as the channel_accessor<std::string>.
 */
template <>
class channel_cvt_accessor<frantic::tstring> {
    frantic::tstring m_default;
    std::size_t m_offset;

  public:
    // Default constructor makes it point to and invalid position.  Putting in an
    // "invalid" state would be more overhead than we want.
    channel_cvt_accessor()
        : m_offset( ( std::numeric_limits<std::size_t>::max )() ) {}

    channel_cvt_accessor( const frantic::tstring& defaultVal )
        : m_default( defaultVal )
        , m_offset( ( std::numeric_limits<std::size_t>::max )() ) {}

    channel_cvt_accessor( std::size_t offset, std::size_t sourceDataArity, data_type_t sourceDataType,
                          const frantic::tstring& channelNameForErrorMessage )
        : m_offset( offset ) {
        data_type_t targetDataType = channel_data_type_traits<frantic::tstring>::data_type();
        std::size_t targetDataArity = channel_data_type_traits<frantic::tstring>::arity();

        if( sourceDataArity != targetDataArity )
            throw std::runtime_error(
                std::string() + "cvt_channel_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );
    }

    bool is_valid() const { return m_offset != ( std::numeric_limits<std::size_t>::max )(); }
    void reset() { m_offset = ( std::numeric_limits<std::size_t>::max )(); }
    void reset( const frantic::tstring& defaultVal ) {
        m_offset = ( std::numeric_limits<std::size_t>::max )();
        m_default = defaultVal;
    }

    // Functions to use the accessor with a raw memory particle
    frantic::tstring get( const char* rawParticleData ) {
        return string_from_channel_string( rawParticleData + m_offset );
    }

    // Convenience function provided for vector<char>
    frantic::tstring get( const std::vector<char>& rawParticleData ) const {
        return string_from_channel_string( &rawParticleData[0] + m_offset );
    }

    frantic::tstring operator()( const char* rawParticleData ) const {
        return string_from_channel_string( rawParticleData + m_offset );
    }

    /**
     * Convenience function provided for vector<char>
     */
    frantic::tstring operator()( const std::vector<char>& rawParticleData ) const {
        return string_from_channel_string( &rawParticleData[0] + m_offset );
    }

    void set( char* rawParticleData, const frantic::tstring& value ) const {
        if( !is_valid() )
            return;
        assign_channel_string( rawParticleData + m_offset, value );
    }

    // Convenience function provided for vector<char>
    void set( std::vector<char>& rawParticleData, const frantic::tstring& value ) const {
        set( &rawParticleData[0], value );
    }
};

///////////////////////
// Provides access to a particle's data member in a runtime general fashion.
///////////////////////
class channel_general_accessor {
    std::size_t m_arity;         // # of data values for this channel per particle
    data_type_t m_dataType;      // primitive data type
    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    std::size_t m_offset;

    // Function which knows how to do weighted sum combinations of this data type
    frantic::channels::channel_weighted_sum_combine_function_t m_weightedSumFunction;
    frantic::channels::channel_weighted_increment_function_t m_weightedIncrementFunction;

  public:
    // Default constructor makes it point to an invalid position.  Putting in an
    // "invalid" state would be more overhead than we want.
    channel_general_accessor()
        : m_arity( 1 )
        , m_dataType( data_type_uint8 )
        , m_offset( ( std::numeric_limits<std::size_t>::max )() ) {
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
        m_weightedSumFunction = channel_weighted_sum_combine_function( m_dataType );
        m_weightedIncrementFunction = channel_weighted_increment_function( m_dataType );
    }

    channel_general_accessor( std::size_t offset, std::size_t arity, data_type_t dataType )
        : m_arity( arity )
        , m_dataType( dataType )
        , m_offset( offset ) {
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
        m_weightedSumFunction = channel_weighted_sum_combine_function( m_dataType );
        m_weightedIncrementFunction = channel_weighted_increment_function( m_dataType );
    }

    bool is_valid() const { return m_offset != ( std::numeric_limits<std::size_t>::max )(); }
    void reset() { m_offset = ( std::numeric_limits<std::size_t>::max )(); }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    frantic::tstring type_str() const { return channel_data_type_str( m_arity, m_dataType ); }

    std::size_t primitive_size() const { return m_primitiveSize; }

    // Given a pointer to a particle buffer, this returns a pointer to the channel within that buffer
    const char* get_channel_data_pointer( const char* rawParticleData ) const { return rawParticleData + m_offset; }

    // Given a pointer to a particle buffer, this returns a pointer to the channel within that buffer
    char* get_channel_data_pointer( char* rawParticleData ) const { return rawParticleData + m_offset; }

    // Given a pointer to a particle buffer, this returns a pointer to the channel within that buffer
    const char* get_channel_data_pointer( const std::vector<char>& rawParticleData ) const {
        // TODO: Can do extra checking in a DEBUG mode
        return &rawParticleData[0] + m_offset;
    }

    // Given a pointer to a particle buffer, this returns a pointer to the channel within that buffer
    char* get_channel_data_pointer( std::vector<char>& rawParticleData ) const {
        // TODO: Can do extra checking in a DEBUG mode
        return &rawParticleData[0] + m_offset;
    }
    /*
      // This copies this primitive from a source particle to a destination primitive.
      void copy_primitive( char* destPrimitiveBuffer, const char* sourceParticleBuffer ) {
        memcpy( destPrimitiveBuffer, sourceParticleBuffer + m_offset, m_primitiveSize );
      }
    */
    // Copy the channel's data from the source particle buffer into
    // a destination primitive buffer.  The destination buffer must have
    // been previously allocated with a size of at least m_primitiveSize.
    void copy_primitive_from_channel( char* destPrimitiveBuffer, const char* sourceParticleBuffer ) const {
        memcpy( destPrimitiveBuffer, sourceParticleBuffer + m_offset, m_primitiveSize );
    }

    // This copies the data from the source primitive buffer into the
    // destination particle buffer's channel.
    void set_channel_from_primitive( char* destParticleBuffer, const char* sourcePrimitiveBuffer ) const {
        memcpy( destParticleBuffer + m_offset, sourcePrimitiveBuffer, m_primitiveSize );
    }

    // This computes a weighted sum of the input type, saving it to the output
    void weighted_sum( const float* weights, const char* const* data, std::size_t weightCount, char* out ) const {
        m_weightedSumFunction( weights, data, weightCount, m_arity, out );
    }

    // This computes a weighted sum of the input type, saving it to the output
    void weighted_sum( const std::vector<float>& weights, const std::vector<char*>& data, char* out ) const {
        m_weightedSumFunction( &weights[0], &data[0], weights.size(), m_arity, out );
    }

    // This increments input type
    void weighted_increment( float weight, const char* data, char* out ) const {
        m_weightedIncrementFunction( weight, data, m_arity, out );
    }
};

namespace detail {
template <class TDest, class TSrc>
struct static_converter {
    static TDest apply( const void* pSrc ) { return static_cast<TDest>( *static_cast<const TSrc*>( pSrc ) ); }
};

template <class TDest>
struct null_converter {
    static TDest apply( const void* ) { return static_cast<TDest>( 0 ); }
};
} // namespace detail

template <class OutType>
class channel_static_cast_const_accessor {
  private:
    channel_general_accessor m_impl;
    OutType ( *m_pConvertFn )( const void* );

  public:
    void reset() {
        m_impl.reset();
        m_pConvertFn = &detail::null_converter<OutType>::apply;
    }

    void reset( const channel_general_accessor& impl ) {
        m_impl = impl;
        m_pConvertFn = &detail::null_converter<OutType>::apply;

        // Switch on all the static_cast applicable types.
        switch( impl.data_type() ) {
        case frantic::channels::data_type_float16:
            m_pConvertFn = &detail::static_converter<OutType, half>::apply;
            break;
        case frantic::channels::data_type_float32:
            m_pConvertFn = &detail::static_converter<OutType, float>::apply;
            break;
        case frantic::channels::data_type_float64:
            m_pConvertFn = &detail::static_converter<OutType, double>::apply;
            break;
        case frantic::channels::data_type_int8:
            m_pConvertFn = &detail::static_converter<OutType, boost::int8_t>::apply;
            break;
        case frantic::channels::data_type_int16:
            m_pConvertFn = &detail::static_converter<OutType, boost::int16_t>::apply;
            break;
        case frantic::channels::data_type_int32:
            m_pConvertFn = &detail::static_converter<OutType, boost::int32_t>::apply;
            break;
        case frantic::channels::data_type_int64:
            m_pConvertFn = &detail::static_converter<OutType, boost::int64_t>::apply;
            break;
        case frantic::channels::data_type_uint8:
            m_pConvertFn = &detail::static_converter<OutType, boost::uint8_t>::apply;
            break;
        case frantic::channels::data_type_uint16:
            m_pConvertFn = &detail::static_converter<OutType, boost::uint16_t>::apply;
            break;
        case frantic::channels::data_type_uint32:
            m_pConvertFn = &detail::static_converter<OutType, boost::uint32_t>::apply;
            break;
        case frantic::channels::data_type_uint64:
            m_pConvertFn = &detail::static_converter<OutType, boost::uint64_t>::apply;
            break;
        default:
            // TODO: Do we want an exception here?
            break;
        }
    }

    channel_static_cast_const_accessor()
        : m_pConvertFn( &detail::null_converter<OutType>::apply ) {
        // BOOST_STATIC_ASSERT( channel_data_type_traits<OutType>::arity() == 1 );
        if( channel_data_type_traits<OutType>::arity() != 1 )
            throw std::runtime_error( "Invalid output type for channel_static_cast_const_accessor<T>" );
    }

    channel_static_cast_const_accessor( const channel_general_accessor& impl ) {
        // BOOST_STATIC_ASSERT( channel_data_type_traits<OutType>::arity() == 1 );
        if( channel_data_type_traits<OutType>::arity() != 1 )
            throw std::runtime_error( "Invalid output type for channel_static_cast_const_accessor<T>" );

        this->reset( impl );
    }

    bool is_valid() const { return m_impl.is_valid(); }

    OutType get( const void* pSrc ) const {
        return m_pConvertFn( m_impl.get_channel_data_pointer( static_cast<const char*>( pSrc ) ) );
    }

    OutType get( const void* pSrc, size_t offset ) const {
        return m_pConvertFn( m_impl.get_channel_data_pointer( static_cast<const char*>( pSrc ) ) +
                             sizeof_channel_data_type( m_impl.data_type() ) * offset );
    }

    OutType operator()( const void* pSrc ) const { return get( pSrc ); }

    OutType operator()( const void* pSrc, size_t offset ) const { return get( pSrc, offset ); }
};

} // namespace channels
} // namespace frantic
