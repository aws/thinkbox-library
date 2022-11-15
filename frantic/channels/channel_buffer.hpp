// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "frantic/channels/channel_buffer_iterator.hpp"
#include "frantic/channels/channel_map.hpp"
#include "frantic/graphics/raw_byte_buffer.hpp"
#include "frantic/graphics/size3.hpp"
#include "frantic/graphics2d/size2.hpp"
#include "frantic/misc/exception_stream.hpp"

namespace frantic {
namespace channels {

using frantic::exception_stream;
using frantic::graphics::raw_byte_buffer;
using frantic::graphics::size3;
using frantic::graphics2d::size2;

class fast_channel_buffer_iterator_exception : public exception_stream_base<fast_channel_buffer_iterator_exception> {
  public:
    fast_channel_buffer_iterator_exception( const std::string& msg )
        : exception_stream_base<fast_channel_buffer_iterator_exception>( msg ) {}
    fast_channel_buffer_iterator_exception() {}
    virtual ~fast_channel_buffer_iterator_exception() throw() {}
};

class generic_channel_buffer_iterator_exception
    : public exception_stream_base<generic_channel_buffer_iterator_exception> {
  public:
    generic_channel_buffer_iterator_exception( const std::string& msg )
        : exception_stream_base<generic_channel_buffer_iterator_exception>( msg ) {}
    generic_channel_buffer_iterator_exception() {}
    virtual ~generic_channel_buffer_iterator_exception() throw() {}
};

/**
 * The core class for the storage of a buffer of channel map data. This
 * class supports 1, 2, and 3 dimensional data.
 *
 * @author   Brian McKinnon
 * @since    Apr 24, 2007
 */
class channel_buffer {
  public:
    channel_buffer();
    channel_buffer( const channel_map& map );
    channel_buffer( const channel_map& map, int size );
    channel_buffer( const channel_map& map, size2 size );
    channel_buffer( const channel_map& map, size3 size );

    void clear();

    void resize( int size );
    void resize( size2 size );
    void resize( size3 size );

    int get_size() const { return m_size.xsize(); }
    size2 get_size2() const { return size2( m_size.xsize(), m_size.ysize() ); }
    size3 get_size3() const { return m_size; }

    virtual void set_channel_map( const channel_map& map );
    const channel_map& get_channel_map() const { return m_channelMap; }

    template <class NamedType>
    void copy_channel( channel_buffer& src );
    template <typename T>
    void copy_channel( channel_buffer& src, const frantic::tstring& channelName, data_type_t channelType );

    template <class NamedType>
    fast_channel_buffer_iterator<NamedType> get_fast_iterator();
    template <typename T>
    fast_channel_buffer_iterator<T> get_fast_iterator( const frantic::tstring& channelName, data_type_t channelType );
    template <class NamedType, class Converter>
    generic_channel_buffer_iterator<NamedType> get_generic_iterator();
    template <class T, class Converter>
    generic_channel_buffer_iterator<T> get_generic_iterator( const frantic::tstring& channelName,
                                                             data_type_t channelType );

    template <class NamedType>
    const_fast_channel_buffer_iterator<NamedType> get_const_fast_iterator() const;
    template <typename T>
    const_fast_channel_buffer_iterator<T> get_const_fast_iterator( const frantic::tstring& channelName,
                                                                   data_type_t channelType ) const;
    template <class NamedType, class Converter>
    const_generic_channel_buffer_iterator<NamedType> get_const_generic_iterator() const;
    template <class NamedType, class Converter>
    const_generic_channel_buffer_iterator<NamedType> get_const_generic_iterator( const frantic::tstring& channelName,
                                                                                 data_type_t channelType ) const;

    raw_byte_buffer& get_raw_buffer() { return m_channelBuffer; }
    const raw_byte_buffer& get_raw_buffer() const { return m_channelBuffer; }

  protected:
    /// Stores the format of the image data
    channel_map m_channelMap;
    /// The buffer that holds the image data
    raw_byte_buffer m_channelBuffer;
    /// The dimensions of the image
    size3 m_size;
};

template <class NamedType>
void channel_buffer::copy_channel( channel_buffer& src ) {
    return copy_channel<NamedType>( src, NamedType::get_name(), NamedType::get_data_type() );
}

template <typename T>
void channel_buffer::copy_channel( channel_buffer& src, const frantic::tstring& channelName, data_type_t channelType ) {
    try {
        fast_channel_buffer_iterator<T> destIter = this->get_fast_iterator<T>( channelName, channelType );
        fast_channel_buffer_iterator<T> srcIter = src.get_fast_iterator<T>( channelName, channelType );

        if( this->m_size != src.m_size ) {
            std::stringstream strstm;
            strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot copy buffer. "
                   << "The source " << src.m_size << " is not the same size as the destination " << this->m_size << ".";
            throw std::runtime_error( strstm.str() );
        }

        for( int z = 0; z < this->m_size.zsize(); z++ ) {
            for( int y = 0; y < this->m_size.ysize(); y++ ) {
                for( int x = 0; x < this->m_size.xsize(); x++ ) {
                    *destIter = *srcIter;

                    srcIter.next_x();
                    destIter.next_x();
                }
                srcIter.go_x( 0 );
                destIter.go_x( 0 );
                srcIter.next_y();
                destIter.next_y();
            }
            srcIter.go_y( 0 );
            destIter.go_y( 0 );
            srcIter.next_z();
            destIter.next_z();
        }
    } catch( const fast_channel_buffer_iterator_exception& e ) {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot copy buffer. " << channelName
               << " is not a defined in one of the channel maps.";
        throw std::runtime_error( strstm.str() );
    }
}

/** The fast_iterator provides a fast method of accessing
 * data from a channel buffer.  If the NamedType is not
 * the type of data stored in the channel_map then an
 * exception is thrown, and no iterator is returned.
 *
 * \return	iter		A templated fast iterator
 * \throws	runtime_error	Throws a runtime_error if the requested channel
 *							does not exist
 */
template <class NamedType>
fast_channel_buffer_iterator<NamedType> channel_buffer::get_fast_iterator() {
    return get_fast_iterator<NamedType>( NamedType::get_name(), NamedType::get_data_type() );
}

/** The fast_iterator provides a fast method of accessing
 * data from a channel buffer.
 *
 * \return	iter		A templated fast iterator
 * \throws	runtime_error	Throws a runtime_error if the requested channel
 *							does not exist
 */
template <typename T>
fast_channel_buffer_iterator<T> channel_buffer::get_fast_iterator( const frantic::tstring& channelName,
                                                                   data_type_t channelType ) {
    if( m_channelMap.channel_definition_complete() ) {
        for( size_t i = 0; i < m_channelMap.channel_count(); i++ ) {
            if( channelName.compare( m_channelMap[i].name() ) == 0 ) {
                if( channelType == m_channelMap[i].data_type() )
                    return fast_channel_buffer_iterator<T>(
                        m_channelBuffer.begin() + m_channelMap[i].offset(), (long)( m_channelMap.structure_size() ),
                        (long)( m_channelMap.structure_size() * this->m_size.xsize() ),
                        (long)( m_channelMap.structure_size() * this->m_size.xsize() * this->m_size.ysize() ) );
                else {
                    std::stringstream strstm;
                    strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot get fast iterator "
                           << frantic::strings::to_string( channelName ) << " is defined but data type " << channelType
                           << " is incorrect.";
                    throw fast_channel_buffer_iterator_exception( strstm.str() );
                }
            }
        }

        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot get fast iterator "
               << frantic::strings::to_string( channelName ) << " is not a defined channel." << std::endl;
        m_channelMap.dump( strstm );
        throw fast_channel_buffer_iterator_exception( strstm.str() );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Cannot get fast iterator without a valid channel map.";
        throw fast_channel_buffer_iterator_exception( strstm.str() );
    }
    return fast_channel_buffer_iterator<T>( 0, 0, 0, 0 );
}

/** The generic_iterator provides a generic method of accessing
 * data from a channel buffer.  By using a magic conversion the
 * NamedType is returned regardless of the internal data storage.
 *
 * \return	iter		a templated generic iterator
 * \throws	runtime_error	Throws a runtime_error if none of the channels
 *							can be converted to the requested type
 */
template <class NamedType, class Converter>
generic_channel_buffer_iterator<NamedType> channel_buffer::get_generic_iterator() {
    void ( *convertToExtern )( void*, void* ) = 0;
    void ( *convertToIntern )( void*, void* ) = 0;

    if( m_channelMap.channel_definition_complete() ) {
        for( size_t i = 0; i < m_channelMap.channel_count(); i++ ) {
            convertToExtern = Converter::get_conversion_function( m_channelMap[i].name(), m_channelMap[i].data_type(),
                                                                  NamedType::get_name(), NamedType::get_data_type() );
            convertToIntern = Converter::get_conversion_function( NamedType::get_name(), NamedType::get_data_type(),
                                                                  m_channelMap[i].name(), m_channelMap[i].data_type() );
            if( ( convertToExtern != 0 ) && ( convertToIntern != 0 ) ) {
                return generic_channel_buffer_iterator<NamedType>(
                    m_channelBuffer.begin() + m_channelMap[i].offset(), (long)m_channelMap.structure_size(),
                    (long)( m_channelMap.structure_size() * this->m_size.xsize() ),
                    (long)( m_channelMap.structure_size() * this->m_size.xsize() * this->m_size.ysize() ),
                    convertToExtern, convertToIntern );
            }
        }

        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot get generic iterator "
               << frantic::strings::to_string( NamedType::get_name() ) << " no conversion funtion is defined.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Cannot get fast iterator without a valid channel map.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    }
    return generic_channel_buffer_iterator<NamedType>( 0, 0, 0, 0, 0, 0 );
}

/** The generic_iterator provides a generic method of accessing
 * data from a channel buffer.  By using a magic conversion the
 * type T is returned regardless of the internal data storage.
 *
 * \return	iter			a templated generic iterator
 * \throws	runtime_error	Throws a runtime_error if none of the channels
 *							can be converted to the requested type
 */
template <class T, class Converter>
generic_channel_buffer_iterator<T> channel_buffer::get_generic_iterator( const frantic::tstring& channelName,
                                                                         data_type_t channelType ) {
    void ( *convertToExtern )( void*, void* ) = 0;
    void ( *convertToIntern )( void*, void* ) = 0;

    if( m_channelMap.channel_definition_complete() ) {
        for( size_t i = 0; i < m_channelMap.channel_count(); i++ ) {
            if( m_channelMap[i].name().compare( channelName ) == 0 ) {

                convertToExtern = Converter::get_conversion_function(
                    m_channelMap[i].name(), m_channelMap[i].data_type(), channelName, channelType );
                convertToIntern = Converter::get_conversion_function( channelName, channelType, m_channelMap[i].name(),
                                                                      m_channelMap[i].data_type() );
                if( ( convertToExtern != 0 ) && ( convertToIntern != 0 ) ) {
                    return generic_channel_buffer_iterator<T>(
                        m_channelBuffer.begin() + m_channelMap[i].offset(), m_channelMap.structure_size(),
                        m_channelMap.structure_size() * this->m_size.xsize(),
                        m_channelMap.structure_size() * this->m_size.xsize() * this->m_size.ysize(), convertToExtern,
                        convertToIntern );
                }
            }
        }

        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot get generic iterator " << channelName
               << " no conversion funtion is defined.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Cannot get fast iterator without a valid channel map.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    }
    return generic_channel_buffer_iterator<T>( 0, 0, 0, 0, 0, 0 );
}

/** The const_fast_iterator provides a fast method of accessing
 * data from a constant channel buffer.  If the NamedType is not
 * the type of data stored in the channel_map then an
 * exception is thrown, and no iterator is returned.
 *
 * \return	iter		A templated fast iterator
 * \throws	runtime_error	Throws a runtime_error if the requested channel
 *							does not exist
 */
template <class NamedType>
const_fast_channel_buffer_iterator<NamedType> channel_buffer::get_const_fast_iterator() const {
    return get_const_fast_iterator<NamedType>( NamedType::get_name(), NamedType::get_data_type() );
}

/** The fast_iterator provides a fast method of accessing
 * data from a channel buffer.
 *
 * \return	iter		A templated fast iterator
 * \throws	runtime_error	Throws a runtime_error if the requested channel
 *							does not exist
 */
template <typename T>
const_fast_channel_buffer_iterator<T> channel_buffer::get_const_fast_iterator( const frantic::tstring& channelName,
                                                                               data_type_t channelType ) const {
    if( m_channelMap.channel_definition_complete() ) {
        for( size_t i = 0; i < m_channelMap.channel_count(); i++ ) {
            if( ( channelName.compare( m_channelMap[i].name() ) == 0 ) &&
                ( channelType == m_channelMap[i].data_type() ) ) {
                return const_fast_channel_buffer_iterator<T>(
                    m_channelBuffer.begin() + m_channelMap[i].offset(), (long)m_channelMap.structure_size(),
                    (long)m_channelMap.structure_size() * this->m_size.xsize(),
                    (long)m_channelMap.structure_size() * this->m_size.xsize() * this->m_size.ysize() );
            }
        }

        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot get const fast iterator "
               << frantic::strings::to_string( channelName ) << " is not a defined channel.";
        throw fast_channel_buffer_iterator_exception( strstm.str() );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Cannot get const fast iterator without a valid channel map.";
        throw fast_channel_buffer_iterator_exception( strstm.str() );
    }
    return const_fast_channel_buffer_iterator<T>( 0, 0, 0, 0 );
}

/** The generic_iterator provides a generic method of accessing
 * data from a channel buffer.  By using a magic conversion the
 * NamedType is returned regardless of the internal data storage.
 *
 * \return	iter		a templated generic iterator
 * \throws	runtime_error	Throws a runtime_error if none of the channels
 *							can be converted to the requested type
 */
template <class NamedType, class Converter>
const_generic_channel_buffer_iterator<NamedType> channel_buffer::get_const_generic_iterator() const {
    void ( *convertToExtern )( void*, void* ) = 0;

    if( m_channelMap.channel_definition_complete() ) {
        for( size_t i = 0; i < m_channelMap.channel_count(); i++ ) {
            convertToExtern = Converter::get_conversion_function( m_channelMap[i].name(), m_channelMap[i].data_type(),
                                                                  NamedType::get_name(), NamedType::get_data_type() );
            if( convertToExtern != 0 ) {
                return const_generic_channel_buffer_iterator<NamedType>(
                    m_channelBuffer.begin() + m_channelMap[i].offset(), m_channelMap.structure_size(),
                    m_channelMap.structure_size() * this->m_size.xsize(),
                    m_channelMap.structure_size() * this->m_size.xsize() * this->m_size.ysize(), convertToExtern );
            }
        }

        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot get const generic iterator "
               << frantic::strings::to_string( NamedType::get_name() ) << " no conversion funtion is defined.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Cannot get const generic iterator without a valid channel map.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    }
}

/** The generic_iterator provides a generic method of accessing
 * data from a channel buffer.  By using a magic conversion the
 * NamedType is returned regardless of the internal data storage.
 *
 * \return	iter		a templated generic iterator
 * \throws	runtime_error	Throws a runtime_error if none of the channels
 *							can be converted to the requested type
 */
template <class T, class Converter>
const_generic_channel_buffer_iterator<T>
channel_buffer::get_const_generic_iterator( const frantic::tstring& channelName, data_type_t channelType ) const {
    void ( *convertToExtern )( void*, void* ) = 0;

    if( m_channelMap.channel_definition_complete() ) {
        for( size_t i = 0; i < m_channelMap.channel_count(); i++ ) {
            if( m_channelMap[i].name().compare( channelName ) == 0 ) {
                convertToExtern = Converter::get_conversion_function(
                    m_channelMap[i].name(), m_channelMap[i].data_type(), channelName, channelType );
                if( ( convertToExtern != 0 ) ) {
                    return const_generic_channel_buffer_iterator<T>(
                        m_channelBuffer.begin() + m_channelMap[i].offset(), m_channelMap.structure_size(),
                        m_channelMap.structure_size() * this->m_size.xsize(),
                        m_channelMap.structure_size() * this->m_size.xsize() * this->m_size.ysize(), convertToExtern );
                }
            }
        }

        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Cannot get generic iterator " << channelName
               << " no conversion funtion is defined.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Cannot get fast iterator without a valid channel map.";
        throw generic_channel_buffer_iterator_exception( strstm.str() );
    }
    return const_generic_channel_buffer_iterator<T>( 0, 0, 0, 0, 0, 0 );
}

} // namespace channels
} // namespace frantic
