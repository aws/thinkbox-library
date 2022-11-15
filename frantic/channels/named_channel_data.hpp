// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <half.h>

#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <frantic/graphics/quat4f.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics/vector4f.hpp>

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/color3h.hpp>

namespace frantic {
namespace channels {

/**
 * @note IF YOU CHANGE THIS, be VERY careful to update dependent functions. Ex. is_channel_data_type_float(),
 * is_channel_data_type_pod().
 *
 * Currently, all the data types here are POD (plain old data) except for the string type.  The string type
 * must be dealt with specially, as follows:
 *
 * <h1>Named Channel String Type</h1>
 *
 * First, an empty string is represented by a NULL pointer.  Thus, the semantics of initializing a channel_map
 * block of memory by memsetting it to 0 can be preserved.  Still, for future expansion to other types which
 * may not match this behavior, a construction/destruction pair of functions in classes like channel_map should be
 * used.
 *
 * Second, an allocated string is immutable.  The main reasoning behind this is to allow for the use of reference
 * counting so that if a particle is copied sequentially through a number of memory blocks, this doesn't incur a
 * memory allocation per copy.  This choice will also, I believe, positively impact multithreading, but be
 * aware that the initial implementation is NOT thread safe.
 *
 * Finally, this string will also be used as the string representation for the Frantic/Flood Language compiler which
 * is planned.  This is of course forcing some design choices on that language, but if we later decide we want different
 * string semantics in that language, we should change the representation there and here simultaneously.
 *
 * All string access should be through the set of string helper functions that can be used to construct, assign,
 * copy, and destruct strings.
 *
 * A string is allocated with malloc/free, and the string layout is as follows:
 *
 *  - 4 bytes -- Integer reference count.
 *  - 4 bytes -- Number of bytes of string data (N).
 *  - N bytes -- The string data.
 *  - 1 byte  -- A null terminator, so the string data can be treated directly as a C-style string.
 */
enum data_type_t {
    data_type_invalid,
    data_type_int8,
    data_type_int16,
    data_type_int32,
    data_type_int64,
    data_type_uint8,
    data_type_uint16,
    data_type_uint32,
    data_type_uint64,
    data_type_float16,
    data_type_float32,
    data_type_float64,
    data_type_string // The string data type is special - it is stored as a pointer.  Extra work for managing lifetime
                     // is necessary.
};

/**
 * This function returns the string form of the given data type.
 *
 * @param  type  The data type to convert into a string representation.
 */
const frantic::tchar* channel_data_type_str( data_type_t type );

/**
 * Convert a data type and arity into a string representation, like "int8", "float32[3]", etc.
 *
 * @param  arity  The number of data elements in this type.
 * @param  type   The basic data type.
 */
frantic::tstring channel_data_type_str( std::size_t arity, data_type_t type );

/**
 * Convert a string into a data type enum.  Input is like "int8", "float32", etc.
 *
 * If the input string is not a valid data type, it returns data_type_invalid.
 *
 * @param  dataTypeStr  The string representation of the data type.
 */
data_type_t channel_data_type_from_string( const std::string& dataTypeStr );

// HACK: Remove this when fixing channel_data_type_from_string() to correctly use frantic::tstring
inline data_type_t channel_data_type_from_string( const std::wstring& dataTypeStr ) {
    return channel_data_type_from_string( frantic::strings::to_string( dataTypeStr ) );
}

/**
 * Returns a data type and arity pair for the given data type string.  Input is like "int8", "float32[3]", etc.
 *
 * If the input string is not a valid data type, the returned data type in the pair will be data_type_invalid.
 *
 * @param  dataTypeStr  The string representation of the data type.
 */
std::pair<data_type_t, std::size_t> channel_data_type_and_arity_from_string( const std::string& dataTypeStr );

// HACK: Remove this when fixing channel_data_type_and_arity_from_string() to correctly use frantic::tstring
inline std::pair<data_type_t, std::size_t> channel_data_type_and_arity_from_string( const std::wstring& dataTypeStr ) {
    return channel_data_type_and_arity_from_string( frantic::strings::to_string( dataTypeStr ) );
}

/**
 * Returns the number of bytes in an element of the given data type.
 *
 * @param  type  The data type.
 */
size_t sizeof_channel_data_type( data_type_t type );

/**
 * This parses an input string, saving the resulting value to the provided memory buffer.  It's the responsibility
 * of the caller to ensure that the output buffer has enough room for the data (up to 8 bytes, currently).
 *
 * @param  type            The data type of the value to parse
 * @param  data            The string containing the value to parse.
 * @param  outValueBuffer  A pointer to the memory where the parsed value should go.
 */
void parse_channel_value_from_string( data_type_t type, const std::string& data, char* outValueBuffer );

/**
 * This prints a tuple of the given type to a stream, separated by the given separator.
 * For example, you could use a "," to produce part of a .csv file output.
 *
 * @param  out        The target ostream for printing.
 * @param  separator  A string to output between each element in the data.
 * @param  arity      The number of basic data elements in the data.
 * @param  dataType   The basic data type.
 * @param  data       A pointer to the memory containing the data to print.
 */
void channel_data_type_print( std::ostream& out, const std::string& separator, std::size_t arity, data_type_t dataType,
                              const char* data );

/**
 * @overload For wide character streams
 */
void channel_data_type_print( std::wostream& out, const std::wstring& separator, std::size_t arity,
                              data_type_t dataType, const char* data );

/**
 * This prints a tuple of the given type to a stream, separating each element by a space.
 *
 * @param  out        The target ostream for printing.
 * @param  arity      The number of basic data elements in the data.
 * @param  dataType   The basic data type.
 * @param  data       A pointer to the memory containing the data to print.
 */
inline void channel_data_type_print( std::ostream& out, std::size_t arity, data_type_t dataType, const char* data ) {
    channel_data_type_print( out, " ", arity, dataType, data );
}

/**
 * @overload For wide character streams
 */
inline void channel_data_type_print( std::wostream& out, std::size_t arity, data_type_t dataType, const char* data ) {
    channel_data_type_print( out, L" ", arity, dataType, data );
}

/**
 * Returns true if the given data type is a floating point type, false otherwise.
 *
 * @param  type  The data type to check.
 */
inline bool is_channel_data_type_float( data_type_t type ) {
    return type >= data_type_float16 && type <= data_type_float64;
}

/**
 * Returns true if the given data type is an integer type (signed or unsigned), false otherwise.
 *
 * @param  type  The data type to check.
 */
inline bool is_channel_data_type_int( data_type_t type ) { return type >= data_type_int8 && type <= data_type_uint64; }

/**
 * Returns true if the data type is unsigned, false otherwise.
 *
 * @param  type  The data type to check.
 */
inline bool is_channel_data_type_unsigned( data_type_t type ) {
    return type >= data_type_uint8 && type <= data_type_uint64;
}

/**
 * Returns true if the data type is unsigned, false otherwise.
 *
 * @param  type  The data type to check.
 */
inline bool is_channel_data_type_signed( data_type_t type ) {
    return type >= data_type_int8 && type <= data_type_int64;
}

/**
 * Returns true if the data type is plain old data, false otherwise.  Currently
 * only the channel string type is not plain old data.
 *
 * @param  type  The data type to check.
 */
inline bool is_channel_data_type_pod( data_type_t type ) { return type != data_type_string; }

/////////////////////
// CHANNEL STRING HELPER FUNCTIONS
// The channel string type is a special string layout designed to fit within channel data.  See the comments
// with the data_type_t enumeration.
/////////////////////

/**
 * This function constructs a new channel string object.  It assumes the memory it is writing to is uninitialized,
 * so should only be used in that context.  If a string has already been constructed there, use assign_channel_string
 * instead.
 *
 * @param  data  A pointer to the data (generally within a channel_map managed memory block) for the string data.
 *               In the current implementation, this is pointing to a char* pointer under the hood.
 * @param  init  The newly constructed string will have a copy of the value from this string.
 */
void construct_channel_string( char* data, const frantic::tstring& init );

/**
 * This function constructs a new channel string object and initializes it to the empty string.  It assumes the memory
 * it is writing to is uninitialized, so should only be used in that context.  If a string has already been
 * constructed there, use assign_channel_string instead.
 *
 * @param  data  A pointer to the data (generally within a channel_map managed memory block) for the string data.
 *               In the current implementation, this is pointing to a char* pointer under the hood.
 */
void construct_channel_string( char* data );

/**
 * This function destructs an existing channel string.  This function must be called when a block of memory is
 * going out of its lifetime scope as a channel_map memory block or other mechanism which uses the channel string.
 *
 * @param  data  A pointer to the data (generally within a channel_map managed memory block) for the string data.
 *               In the current implementation, this is pointing to a char* pointer under the hood.
 */
void destruct_channel_string( char* data );

/**
 * This function assigns a value to an existing channel string object.  It assumes the memory it is writing to already
 * contains a valid channel string.
 *
 * @param  data  A pointer to the data (generally within a channel_map managed memory block) for the string data.
 *               In the current implementation, this is pointing to a char* pointer under the hood.
 * @param  init  The output string will have a copy of the value from this string.
 */
void assign_channel_string( char* data, const frantic::tstring& init );

/**
 * This function copies the string value pointed to by sourceData into the string value pointed to by destData.
 * It assumes that both destData and sourceData point at valid channel strings.
 *
 * @param  destData    A pointer to the data (generally within a channel_map managed memory block) for the destination
 * string data. In the current implementation, this is pointing to a char* pointer under the hood.
 * @param  sourceData  A pointer to the data (generally within a channel_map managed memory block) for the source string
 * data. In the current implementation, this is pointing to a char* pointer under the hood.
 * @param  arity       The number of string variables in a row to copy.
 */
void copy_channel_string( char* destData, const char* sourceData, std::size_t arity = 1 );

/**
 * This function returns a C-style string from an existing channel string.  This operation does not allocate any
 * new memory.
 *
 * @param  data  A pointer to the data (generally within a channel_map managed memory block) for the string data.
 *               In the current implementation, this is pointing to a char* pointer under the hood.
 */
const frantic::tchar* cstring_from_channel_string( const char* data );

/**
 * This function returns a C++ std::string from an existing channel string.  Because this function will allocate
 * new memory to construct the string class, it is recommended that you use cstring_from_channel_string where
 * beneficial.
 *
 * @param  data  A pointer to the data (generally within a channel_map managed memory block) for the string data.
 *               In the current implementation, this is pointing to a char* pointer under the hood.
 */
frantic::tstring string_from_channel_string( const char* data );

/**
 * This function returns the length of the channel string pointed to by the data pointer.
 *
 * @param  data  A pointer to the data (generally within a channel_map managed memory block) for the string data.
 *               In the current implementation, this is pointing to a char* pointer under the hood.
 */
std::size_t length_from_channel_string( const char* data );

/////////////////////
// NAMED CHANNEL DATA TYPE TRAITS
// The channel_data_type_traits provide a way to map from various data types into [arity,data_type_t] pairs.  This way,
// if code requests a named channel of an invalid type, it will error out because the type isn't specialized in these
// type traits, and for valid types, this provides a general way to do the mapping.
/////////////////////

// Default type traits don't define any values, so if you try to use them, it produces a compile error.
// Every valid type has a specialization of these type traits.
template <class T>
struct channel_data_type_traits;

// For integer types, the barycentric combine gets the nearest value instead of a weighted sum.
template <class T>
inline T barycentric_integer_combine( const frantic::graphics::vector3f& barycentricCoords, T a, T b, T c ) {
    if( barycentricCoords.x >= barycentricCoords.y ) {
        if( barycentricCoords.x >= barycentricCoords.z ) {
            return a;
        } else { // x < z
            return c;
        }
    } else { // x < y
        if( barycentricCoords.y >= barycentricCoords.z ) {
            return b;
        } else { // y < z
            return c;
        }
    }
}

// For integer types, the combine function gets the value corresponding to the largest weight
template <class T>
inline void weighted_sum_integer_combine_general( const float* weights, const char* const* data,
                                                  std::size_t weightCount, std::size_t arity, char* out ) {
    // Find the index of the biggest weight
    unsigned biggestWeightIndex = 0;
    for( unsigned i = 1; i < weightCount; ++i ) {
        if( weights[i] > weights[biggestWeightIndex] )
            biggestWeightIndex = i;
    }

    if( out != data[biggestWeightIndex] ) {
        memcpy( out, data[biggestWeightIndex], arity * sizeof( T ) );
    }
}

// This general function is used to do weighted sum queries in accessor classes for named channels
// NOTE: Currently we're implementing linear interpolation and barycentric combination using this, but later
//       we may want to special case those for performance.
template <class T>
inline void weighted_sum_float_combine_general( const float* weights, const char* const* data, std::size_t weightCount,
                                                std::size_t arity, char* out ) {
    std::size_t offset = 0;

    while( arity-- > 0 ) {
        T result = 0;
        for( unsigned i = 0; i < weightCount; ++i )
            result += weights[i] * ( *reinterpret_cast<const T*>( data[i] + offset ) );
        *reinterpret_cast<T*>( out + offset ) = result;
        offset += sizeof( T );
    }
}

// Special case the combination of half values so the math is done as floats for all intermediate values
inline void weighted_sum_float_combine_half( const float* weights, const char* const* data, std::size_t weightCount,
                                             std::size_t arity, char* out ) {
    std::size_t offset = 0;
    while( arity-- > 0 ) {
        float result = 0;

        for( unsigned i = 0; i < weightCount; ++i )
            result += weights[i] * ( *reinterpret_cast<const half*>( data[i] + offset ) );

        *reinterpret_cast<half*>( out + offset ) = result;

        offset += sizeof( half );
    }
}

// This general function is used to do weighted sum queries in accessor classes for named channels
template <class T>
inline void weighted_increment_general( float weight, const char* data, std::size_t arity, char* out ) {
    std::size_t offset = 0;
    while( arity-- > 0 ) {
        T& inc = reinterpret_cast<T*>( out )[offset];
        inc = inc + (T)( weight * (float)( reinterpret_cast<const T*>( data )[offset] ) );
        ++offset;
    }
}

// Scale data by the specified weight, and store the result in out.
namespace channel_scale_general_detail {
template <class DataType>
struct intermediate {
    typedef double type;
};
// avoid a warning by casting half to float instead of to double
template <>
struct intermediate<half> {
    typedef float type;
};
} // namespace channel_scale_general_detail
template <class T>
inline void channel_scale_general( double scale, const char* data, std::size_t arity, char* out ) {
    typedef typename channel_scale_general_detail::intermediate<T>::type intermediate_type;

    std::size_t offset = 0;
    while( arity-- > 0 ) {
        T& outPrimitive = reinterpret_cast<T*>( out )[offset];
        outPrimitive = (T)( static_cast<intermediate_type>( scale ) *
                            static_cast<intermediate_type>( reinterpret_cast<const T*>( data )[offset] ) );
        ++offset;
    }
}

// Signed integer types
template <>
struct channel_data_type_traits<boost::int8_t> {
    typedef boost::int8_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_int8; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<boost::int16_t> {
    typedef boost::int16_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_int16; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<boost::int32_t> {
    typedef boost::int32_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_int32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<boost::int64_t> {
    typedef boost::int64_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_int64; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

// Unsigned integer types
template <>
struct channel_data_type_traits<boost::uint8_t> {
    typedef boost::uint8_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_uint8; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<boost::uint16_t> {
    typedef boost::uint16_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_uint16; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<boost::uint32_t> {
    typedef boost::uint32_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_uint32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<boost::uint64_t> {
    typedef boost::uint64_t value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_uint64; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_integer_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

// Floating point types
template <>
struct channel_data_type_traits<half> {
    typedef half value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_float16; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        // Use a special function for half, which uses floats for accumulation to avoid clipping problems
        weighted_sum_float_combine_half( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<float> {
    typedef float value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_float32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_float_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

template <>
struct channel_data_type_traits<double> {
    typedef double value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_float64; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
    inline static void weighted_sum_combine_general( const float* weights, const char* const* data,
                                                     std::size_t weightCount, std::size_t arity, char* out ) {
        weighted_sum_float_combine_general<value_type>( weights, data, weightCount, arity, out );
    }
};

// Vector types
template <>
struct channel_data_type_traits<frantic::graphics2d::vector2f> {
    typedef frantic::graphics2d::vector2f value_type;
    inline static size_t arity() { return 2; }
    inline static data_type_t data_type() { return data_type_float32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
};

template <>
struct channel_data_type_traits<frantic::graphics2d::vector2> {
    typedef frantic::graphics2d::vector2 value_type;
    inline static size_t arity() { return 2; }
    inline static data_type_t data_type() { return data_type_int32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
};

template <class FloatType>
struct channel_data_type_traits<frantic::graphics::vector3t<FloatType>> {
    typedef frantic::graphics::vector3t<FloatType> value_type;
    inline static size_t arity() { return 3; }
    inline static data_type_t data_type() { return channel_data_type_traits<FloatType>::data_type(); }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
};

template <class FloatType>
struct channel_data_type_traits<frantic::graphics::vector4t<FloatType>> {
    typedef frantic::graphics::vector4t<FloatType> value_type;
    inline static size_t arity() { return 4; }
    inline static data_type_t data_type() { return channel_data_type_traits<FloatType>::data_type(); }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
};

template <class FloatType>
struct channel_data_type_traits<frantic::graphics::quat4t<FloatType>> {
    typedef frantic::graphics::quat4t<FloatType> value_type;
    inline static size_t arity() { return 4; }
    inline static data_type_t data_type() { return channel_data_type_traits<FloatType>::data_type(); }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
};

template <>
struct channel_data_type_traits<frantic::graphics::vector3> {
    typedef frantic::graphics::vector3 value_type;
    inline static size_t arity() { return 3; }
    inline static data_type_t data_type() { return data_type_int32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentric_integer_combine( barycentricCoords, a, b, c );
    }
};

// Matrix types
template <class FloatType>
struct channel_data_type_traits<frantic::graphics::transform4t<FloatType>> {
    typedef frantic::graphics::transform4t<FloatType> value_type;
    inline static size_t arity() { return 16; }
    inline static data_type_t data_type() { return channel_data_type_traits<FloatType>::data_type(); }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords,
                                                  const value_type& a, const value_type& b, const value_type& c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
};

// Boundbox types
template <class FloatType>
struct channel_data_type_traits<frantic::graphics::boundbox3t<FloatType>> {
    typedef frantic::graphics::boundbox3t<FloatType> value_type;
    inline static size_t arity() { return 6; }
    inline static data_type_t data_type() { return channel_data_type_traits<FloatType>::data_type(); }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
};

// Color types
template <>
struct channel_data_type_traits<frantic::graphics::color3f> {
    typedef frantic::graphics::color3f value_type;
    inline static size_t arity() { return 3; }
    inline static data_type_t data_type() { return data_type_float32; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
};

template <>
struct channel_data_type_traits<frantic::graphics::color3h> {
    typedef frantic::graphics::color3h value_type;
    inline static size_t arity() { return 3; }
    inline static data_type_t data_type() { return data_type_float16; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static value_type barycentric_combine( const frantic::graphics::vector3f& barycentricCoords, value_type a,
                                                  value_type b, value_type c ) {
        return barycentricCoords.x * a + barycentricCoords.y * b + barycentricCoords.z * c;
    }
};

// String types
template <>
struct channel_data_type_traits<frantic::tstring> {
    typedef std::string value_type;
    inline static size_t arity() { return 1; }
    inline static data_type_t data_type() { return data_type_string; }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
    inline static void weighted_sum_combine_general( const float* /*weights*/, const char* const* /*data*/,
                                                     std::size_t /*weightCount*/, std::size_t /*arity*/,
                                                     char* /*out*/ ) {
        throw std::runtime_error( "channel_weighted_sum_combine_function: Attempted to apply a weighted sum combining "
                                  "function for the string type.  This is not currently supported." );
    }
    inline static void weighted_increment_general( float /*weight*/, const char* /*data*/, std::size_t /*arity*/,
                                                   char* /*out*/ ) {
        throw std::runtime_error(
            "channel_weighted_increment_function: Attempted to apply a weighted increment function "
            "for the string type.  This is not currently supported." );
    }
};

// std::pair<T, T>
template <typename T>
struct channel_data_type_traits<std::pair<T, T>> {
    typedef std::pair<T, T> value_type;
    inline static size_t arity() { return 2; }
    inline static data_type_t data_type() { return channel_data_type_traits<T>::data_type(); }
    inline static frantic::tstring type_str() { return channel_data_type_str( arity(), data_type() ); }
};

///////////////////////
// Helpers for doing named channel type conversions
///////////////////////

/**
 * This function type converts an input array of a given number of elements to
 * an output array of a given number of elements.  The total number of input bytes read
 * will be (input data type size)*arity, and the total number of output bytes written
 * will be (output data type size)*arity.
 *
 * @param  dest    The destination memory.
 * @param  source  The source memory.
 * @param  arity   The number of basic data types (i.e. float, int, etc) to copy.
 */
typedef void ( *channel_type_convertor_function_t )( char* dest, const char* source, std::size_t arity );

/**
 * This returns a function which does a type conversion of an array of data
 *
 * @param  sourceType  The input type for the converter function.
 * @param  destType    The output type for the converter function.
 * @param  channelNameForErrorMessage  The channel name to use when producing an error message.
 */
channel_type_convertor_function_t
get_channel_type_convertor_function( data_type_t sourceType, data_type_t destType,
                                     const frantic::tstring& channelNameForErrorMessage );

/**
 * This function returns the promoted data type from the two inputs.  For instance, a float and a double promotes
 * to a double, or an int16 and an int32 promotes to an int32.
 *
 * @param  t1  The first data type.
 * @param  t2  The second data type.
 */
data_type_t promote_types( data_type_t t1, data_type_t t2 );

/**
 * Given a type conversion function, this returns a string indicating which function it is.
 *
 * @note The MS compiler collapses together identical functions.  The functions convert_uint32_to_uint64 and
 * convert_uint32_to_int64 produce identical assembly language, so the function pointer for them ends up being the same.
 * Keep this in mind when looking at debug dumps of a channel_map_adaptor.
 *
 * @param  ptc  The channel type converter function to convert to a string for debug purposes.
 */
const char* get_channel_type_convertor_debug_string( channel_type_convertor_function_t ptc );

/**
 * Test to confirm whether a particular name is valid for a channel.  Basically, it's limited to simple identifiers.
 *
 * @note this function's behaviour is assumed to match the PRT spec for
 *       channel names: 'Must match the regex "[a-zA-Z_][0-9a-zA-Z_]*"'.
 *       Changing this behaviour will break things, so please be very
 *       careful!
 *
 * @param  name  The channel name to test for validity.
 */
bool is_valid_channel_name( const std::string& name );
bool is_valid_channel_name( const std::wstring& name );

///////////////////////
// Helpers for doing weighted sum combinations
///////////////////////

/**
 * This is the type of functions used to do barycentric combining in the trimesh3 general named vertex channel accessor
 * class.
 *
 * For example, to do a linear interpolation between two sets of data, you would have weightCount == 2 and arity == # of
 * elements to interpolate.  In such a case you can often make arity == dataType.arity * # array elements to interpolate
 * whole arrays at once.
 *
 * To do a reconstruction of a value in an rle level set, you will make one weight and one data index for each input.
 * See the implicit surface policy for filtered reconstruction of the level set to see how this can work.
 *
 * @param  weights      An array of floating point weights.
 * @param  data         An array of pointers to data elements, each element corresponding to a weight.
 * @param  weightCount  The number of weights and data elements.
 * @param  arity        The number of basic data types (i.e. float, int, etc) in each data elemment.
 * @param  out          A pointer to a data element for output.
 */
typedef void ( *channel_weighted_sum_combine_function_t )( const float* weights, const char* const* data,
                                                           std::size_t weightCount, std::size_t arity, char* out );

/**
 *
 * Weighted sum function intended to operate on one variable in an array of structures data set,
 * for example, to blend one channel from many input particles.
 *
 * @param  weights        An array of floating point weights.
 * @param  channelOffset  An offset, measured in bytes, that will be added to the particles
 *                        to get the input data element pointers.
 * @param  particles      An array of pointers, each element corresponding to a weight.  The
 *                        channelOffset will be added to these pointers to get the data elements.
 * @param  weightCount    The number of weights and data elements.
 * @param  arity          The number of basic data types (i.e. float, int, etc) in each data elemment.
 * @param  out            A pointer to a data element for output.
 */
typedef void ( *offset_input_channel_weighted_sum_combine_function_t )( const float* weights, std::size_t channelOffset,
                                                                        const char* const* particles,
                                                                        std::size_t weightCount, std::size_t arity,
                                                                        char* out );

/**
 * This is the type of functions used to do flexible channel incrementing.  It allows you to weight the
 * incremental value.
 *
 * @note This is equivalent to the += operator.
 *
 * @param  weight      This is the float weight that the incrementing data will be mulitplied by.
 * @param  data        A pointer to the data to increment by.
 * @param  arity       The number of basic data types (i.e. float, int, etc) in each data elemment.
 * @param  out         A pointer to a data element to be incremented.
 */
typedef void ( *channel_weighted_increment_function_t )( float weight, const char* data, std::size_t arity, char* out );

/**
 * This is the type of functions used to do channel scaling.
 *
 * @note This is equivalent to the *= operator.
 *
 * @param  weight  This is the double weight that the data will be multiplied by.
 * @param  data    A pointer to the data to scale.
 * @param  arity   The number of basic data types (i.e. float, int, etc) in each data elemment.
 * @param  out     A pointer to the output data.
 */
typedef void ( *channel_scale_function_t )( double weight, const char* data, std::size_t arity, char* out );

/**
 * This is a type of function that performs a binary operation (such as addition or multiplication)
 * on the channels passed in, storing the result in a third channel.
 *
 * @param  srcLHS      The left-hand value to the binary op
 * @param  srcRHS      The right-hand value to the binary op
 * @param  arity       The number of basic data types (i.e. float, int, etc) in each value.
 * @param  out         A pointer to a value where to store the result
 */
typedef void ( *channel_binary_function_t )( const char* srcLHS, const char* srcRHS, std::size_t arity, char* out );

/**
 * This function will perform a range re-mapping on a channel
 * This is used for mapping between data formats with dissimilar ranges
 */
typedef void ( *channel_range_map_function_t )( double fromLB, double fromUB, double toLB, double toUB,
                                                const char* srcRHS, std::size_t arity, char* out );

/**
 * This returns a function which does a weighted sum combination of data.
 *
 * @param  type  The data type that the returned function should use for weighted sums.
 */
channel_weighted_sum_combine_function_t channel_weighted_sum_combine_function( data_type_t type );

/**
 * @param  type  The data type that the returned function should use for weighted sums.
 * @return a function which does a weighted sum combine of the data, with a byte offset applied
 *         to the input data pointers.
 */
offset_input_channel_weighted_sum_combine_function_t
offset_input_channel_weighted_sum_combine_function( data_type_t type );

/**
 *  This function performs a weighted sum of integer data and converts
 * the result into the out_type.
 *
 * @param  weights      An array of floating point weights.
 * @param  data         An array of pointers to data elements of type input_type, each element corresponding to a
 * weight.
 * @param  weightCount  The number of weights and data elements.
 * @param  arity        The number of basic data types (i.e. float, int, etc) in each data element.
 * @param  out          A pointer to a data element of type output_type for output.
 */
template <class out_type, class in_type>
inline void weighted_sum_integer_combine_and_convert( const float* weights, const char* const* data,
                                                      std::size_t weightCount, std::size_t arity, char* out ) {
    // Find the index of the biggest weight
    unsigned biggestWeightIndex = 0;
    for( unsigned i = 1; i < weightCount; ++i ) {
        if( weights[i] > weights[biggestWeightIndex] )
            biggestWeightIndex = i;
    }

    frantic::channels::channel_type_convertor_function_t convertor =
        frantic::channels::get_channel_type_convertor_function(
            frantic::channels::channel_data_type_traits<in_type>::data_type(),
            frantic::channels::channel_data_type_traits<out_type>::data_type(), _T("") );
    convertor( out, data[biggestWeightIndex], arity );
}

/**
 *  This function performs a weighted sum of floating-point data and
 * converts the result into the out_type.
 *
 * @param  weights      An array of floating point weights.
 * @param  data         An array of pointers to data elements of type input_type, each element corresponding to a
 * weight.
 * @param  weightCount  The number of weights and data elements.
 * @param  arity        The number of basic data types (i.e. float, int, etc) in each data element.
 * @param  out          A pointer to a data element of type output_type for output.
 */
template <class output_type, class intermediate_type, class input_type>
inline void weighted_sum_float_combine_and_convert( const float* weights, const char* const* data,
                                                    const std::size_t weightCount, std::size_t arity, char* out ) {
    std::size_t inputOffset = 0;
    std::size_t outputOffset = 0;

    while( arity-- > 0 ) {
        intermediate_type result = 0;
        for( std::size_t i = 0; i < weightCount; ++i ) {
            result += static_cast<intermediate_type>( weights[i] ) *
                      static_cast<intermediate_type>( *reinterpret_cast<const input_type*>( data[i] + inputOffset ) );
        }
        *reinterpret_cast<output_type*>( out + outputOffset ) = static_cast<output_type>( result );
        inputOffset += sizeof( input_type );
        outputOffset += sizeof( output_type );
    }
}

/**
 *  Get a function which performs a weighted sum of data of type
 * sourceType and converts it into an output of type destType.
 *
 *  An exception is thrown if the sourceType cannot be converted
 * into the destType
 *
 * @param sourceType the input type for the desired function.
 * @param destType the output type for the desired function.
 * @param channelNameForErrorMessage the channel name to use in
 *		error messages, in case no appropriate conversion
 *		function could be found.
 * @return a function which can perform a weighted sum of data of
 *		type sourceType, and convert it into an output of type
 *		destType.
 */
frantic::channels::channel_weighted_sum_combine_function_t
channel_weighted_sum_combine_and_convert_function( frantic::channels::data_type_t sourceType,
                                                   frantic::channels::data_type_t destType,
                                                   const frantic::tstring& channelNameForErrorMessage );

/**
 * This returns a function which does a weighted increment of data.
 *
 * @note This is equivalent to the += operator.
 *
 * @param  type  The data type that the returned function should use for incrementing.
 */
channel_weighted_increment_function_t channel_weighted_increment_function( data_type_t type );

/**
 *  This returns a function that scales data.
 *
 * @note This is equivalent to the *= operator.
 *
 * @param  type  The data type that the returned function should scale.
 */
channel_scale_function_t channel_scale_function( data_type_t type );

/**
 * This returns a function that does range re-mapping on data.
 *
 * @param type Data type this function maps
 */
channel_range_map_function_t channel_range_map_function( data_type_t type );

/**
 * This returns a function which does an addition of two inputs, storing the result in a third location.
 *
 * @param  type  The data type that the returned function should use for addition.
 */
channel_binary_function_t channel_addition_function( data_type_t type );

/**
 * This returns a function which does a subtraction of two inputs, storing the result in a third location.
 *
 * @param  type  The data type that the returned function should use for subtraction.
 */
channel_binary_function_t channel_subtraction_function( data_type_t type );

/**
 * This returns a function which does a multiplication of two inputs, storing the result in a third location.
 *
 * @param  type  The data type that the returned function should use for multiplication.
 */
channel_binary_function_t channel_multiplication_function( data_type_t type );

/**
 *  This returns a function which finds the maximum of two inputs, storing the maximum in a third location.
 *
 * @note the returned function finds the per-element maximum, which may not
 *		be what you expect if the channel arity is > 1.
 *
 * @param  type  The data type that the returned function should operate on.
 */
channel_binary_function_t channel_maximum_function( data_type_t type );

/**
 * This returns a function which finds the minimum of two inputs, storeing the maximum in a third location
 *
 * @note the returned function finds the per-element minimum, which may not
 *		be what you expect if the channel arity > 1.
 *
 * @param	type The data type that the returned function should operate on.
 */
channel_binary_function_t channel_minimum_function( data_type_t type );

} // namespace channels
} // namespace frantic
