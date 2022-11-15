// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <map>
#include <string>

#include <frantic/channels/channel_map.hpp>

namespace frantic {
namespace channels {

namespace detail {
// For read/write access of a property, we need to specialize the type so that the std::string can propagate
// the string_accessor_helper.
template <class T>
struct property_map_access_type {
    typedef T& type;
};
template <>
struct property_map_access_type<frantic::tstring> {
    typedef frantic::channels::detail::string_accessor_helper type;
};

template <class T>
struct const_property_map_access_type {
    typedef const T& type;
};
template <>
struct const_property_map_access_type<frantic::tstring> {
    typedef frantic::channels::detail::const_string_accessor_helper type;
};
} // namespace detail

/**
 * The property_map class is a container for holding properties.  It's intended use is for things like
 * holding metadata in PRT files and passing properties to simulation clients from a Flood scene source.
 * It's current design is focused less on high performance of all operations, and more on convenience of
 * use for a variety of contexts.  For instance, it isn't optimized for quick addition/removal of properties.
 *
 * For string properties, the memory pointed to by the channel_map is no longer POD (plain old data), so
 * the scope management becomes really important in that case.  The property_map class uses the channel_map
 * functions construct_structure, destruct_structure, as well as the channel_map_adaptor when the channel_map
 * changes.
 */
class property_map {
    // m_channelMap and m_buffer are used for all the numeric properties
    channel_map m_channelMap;
    char* m_buffer;

  public:
    /**
     * Constructs an empty property_map.
     */
    property_map();

    /**
     * Copy constructor for property_map.
     */
    property_map( const property_map& rhs );

    /**
     * Constructs a property_map with the given channel_map for data layout.
     */
    explicit property_map( const channel_map& channelMap );

    /**
     * Destructor for property_map.
     */
    ~property_map();

    /**
     * This returns true if the property_map is empty.
     */
    bool empty() const { return m_channelMap.structure_size() == 0; }

    /**
     * This function clears the property_map so it is empty.
     */
    void clear();

    /**
     * This function swaps two property_map instances.
     */
    void swap( property_map& rhs );

    /**
     * Assignment operator for property_map.
     */
    property_map& operator=( const property_map& rhs );

    void dump( std::ostream& out ) const;
    void dump( std::wostream& out ) const;

    ////////////////////
    // Lower-level access to the internals
    ////////////////////

    /**
     * Returns the channel map being used to store the numeric property data.  Note that
     * there are also string properties not represented in the channel map.
     */
    const channel_map& get_channel_map() const;

    /**
     * This function allows you to change the channel map being used to store the property
     * data.  This function will attempt to preserve the existing data as well as it can.
     * Note that there are also string properties not represented in the channel map.
     *
     * Under some circumstances, in particular when a parameter type is changed in an incompatible
     * way, this function will throw an exception.
     *
     * @param  cm  The channel map to use.
     */
    void set_channel_map( const channel_map& cm );

    /**
     * This function allows you to change the channel map being used to store the property
     * data.  This function will attempt to preserve the existing data as well as it can.
     * Note that there are also string properties not represented in the channel map.
     *
     * By using a swap instead of a copy, it can make some code more efficient, where you don't
     * need to keep the channel_map you're passing into 'cm' around any longer.
     *
     * Under some circumstances, in particular when a parameter type is changed in an incompatible
     * way, this function will throw an exception.
     *
     * @param  cm  The channel map to use.  NOTE that it will be swapped with the channel map internal to the class.
     */
    void set_channel_map_with_swap( channel_map& cm );

    /**
     * Returns a const version of the raw pointer to the numeric buffer.  Note that this pointer
     * may be invalidated by other operations on the property_map.
     */
    const char* get_raw_buffer() const { return m_buffer; }

    /**
     * Returns a non-const version of the raw pointer to the numeric buffer.  Note that this pointer
     * may be invalidated by other operations on the property_map.
     */
    char* get_raw_buffer() { return m_buffer; }

    /**
     * Returns a pointer to the data of the specific given channel
     *
     * @param channel the name of the channel to get
     * @return pointer to the buffer
     */
    char* get_channel_buffer( const frantic::tstring& channel ) {
        channel_general_accessor acc = m_channelMap.get_general_accessor( channel );
        return acc.get_channel_data_pointer( m_buffer );
    }

    /**
     * Returns a pointer to the data of the specific given channel
     *
     * @param channel the name of the channel to get
     * @return pointer to the buffer
     */
    const char* get_channel_buffer( const frantic::tstring& channel ) const {
        channel_general_accessor acc = m_channelMap.get_general_accessor( channel );
        return acc.get_channel_data_pointer( m_buffer );
    }

    ////////////////////
    // General functions.
    ////////////////////

    /**
     * Returns true if the property_map has the given property.
     *
     * @param  propertyName  The name of the property.
     */
    bool has_property( const frantic::tstring& propertyName ) const;

    /**
     * Deletes a property from the property_map. If the property doesn't exist, does nothing.
     *
     * @note This function is VERY SLOW! The property_map is not optimized for dynamic modification.
     *
     * @param propertyName The name of the property.
     */
    void delete_property( const frantic::tstring& propertyName );

    /**
     * Adds a property to the property_map. If the property exists, then it sets the value if the arity and data_type
     * match.
     *
     * @note This function is VERY SLOW! The property_map is not optimized for dynamic modification.
     *
     * @param propertyName The name of the property.
     * @param value The value of the property.
     *
     * @throws runtime_error Throws a runtime error if the data type and arity of the value dont match with the existing
     * property
     */
    template <class T>
    void set_property( const frantic::tstring& propertyName, const T& value ) {
        if( !has_property( propertyName ) ) {
            channel_map newChannelMap;
            for( size_t i = 0; i < m_channelMap.channel_count(); ++i ) {
                newChannelMap.define_channel( m_channelMap[i].name(), m_channelMap[i].arity(),
                                              m_channelMap[i].data_type() );
            }
            newChannelMap.define_channel( propertyName, channel_data_type_traits<T>::arity(),
                                          channel_data_type_traits<T>::data_type() );
            newChannelMap.end_channel_definition();
            set_channel_map_with_swap( newChannelMap );
            set_cvt<T>( propertyName, value );
        } else {
            size_t i = m_channelMap.channel_index( propertyName );
            if( channel_data_type_traits<T>::arity() == m_channelMap[i].arity() &&
                channel_data_type_traits<T>::data_type() == m_channelMap[i].data_type() ) {
                set_cvt<T>( propertyName, value );
            } else {
                throw std::runtime_error(
                    "property_map add_property: Property" + frantic::strings::to_string( propertyName ) +
                    "already exists. The data type and arity of the "
                    "value do not match the existing property. The data type and arity of the existing property are " +
                    frantic::strings::to_string(
                        channel_data_type_str( m_channelMap[i].arity(), m_channelMap[i].data_type() ) ) +
                    ". The data type and arity of the passed value were " +
                    frantic::strings::to_string( channel_data_type_str( channel_data_type_traits<T>::arity(),
                                                                        channel_data_type_traits<T>::data_type() ) ) );
            }
        }
    }

    /**
     * Merges the two property maps into the current one.
     *
     * @param  rhs  The properties to merge into *this.
     */
    void merge_property_map( const property_map& rhs );

    ////////////////////
    // Get functions for ease of use.
    ////////////////////

    /**
     * Const access to numeric property data.
     *
     * @param  propertyName  The name of the property.
     */
    template <class T>
    typename detail::const_property_map_access_type<T>::type get( const frantic::tstring& propertyName ) const {
        return m_channelMap.get_accessor<T>( propertyName )( m_buffer );
    }

    /**
     * Mutable access to numeric property data.
     *
     * @param  propertyName  The name of the property.
     */
    template <class T>
    typename detail::property_map_access_type<T>::type get( const frantic::tstring& propertyName ) {
        return m_channelMap.get_accessor<T>( propertyName )( m_buffer );
    }

    /**
     * Access to numeric property data, returning a boolean indicating whether the property exists.  Note
     * that an exception could still be thrown if the property exists but the types are incompatible.
     *
     * @param  propertyName  The name of the property.
     * @param  outValue  A reference to the location to place the value.
     *
     * @return True if the property exists and was populated, false if the property doesn't exist.
     */
    template <class T>
    bool get( const frantic::tstring& propertyName, T& outValue ) const {
        if( m_channelMap.has_channel( propertyName ) ) {
            outValue = m_channelMap.get_cvt_accessor<T>( propertyName )( m_buffer );
            return true;
        } else {
            return false;
        }
    }

    /**
     * Const access to numeric property data, with type conversion if necessary.
     *
     *
     * @param  propertyName  The name of the property.
     */
    template <class T>
    const T get_cvt( const frantic::tstring& propertyName ) const {
        return m_channelMap.get_cvt_accessor<T>( propertyName )( m_buffer );
    }

    /**
     * Set the numeric property value, with type conversion if necessary.
     *
     *
     * @param  propertyName  The name of the property.
     * @param  value  The new value for the property.
     */
    template <class T>
    void set_cvt( const frantic::tstring& propertyName, const T& value ) const {
        m_channelMap.get_cvt_accessor<T>( propertyName ).set( m_buffer, value );
    }
};

} // namespace channels
} // namespace frantic
