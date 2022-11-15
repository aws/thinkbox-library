// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace frantic {
namespace channels {

/**
 * A base class for handling moves over a 3D buffer.  It is only used for
 * storing step sizes in the X/Y/Z direction, as well as, the current offset
 * in the X/Y/Z direction.  No bounds checking is done, this will remain the
 * users responsability.
 */
template <class IterType>
class iterator_3d_base {
  public:
    IterType& operator++() { return next_x(); }

    IterType& next_x() {
        m_xOffset += m_channelMapSizeX;
        return *(IterType*)this;
    }

    IterType& prev_x() {
        m_xOffset -= m_channelMapSizeX;
        return *(IterType*)this;
    }

    IterType& next_y() {
        m_yOffset += m_channelMapSizeY;
        return *(IterType*)this;
    }

    IterType& prev_y() {
        m_yOffset -= m_channelMapSizeY;
        return *(IterType*)this;
    }

    IterType& next_z() {
        m_zOffset += m_channelMapSizeZ;
        return *(IterType*)this;
    }

    IterType& prev_z() {
        m_zOffset -= m_channelMapSizeZ;
        return *(IterType*)this;
    }

    IterType& next_x( long step ) {
        m_xOffset += step * m_channelMapSizeX;
        return *(IterType*)this;
    }

    IterType& prev_x( long step ) {
        m_xOffset -= step * m_channelMapSizeX;
        return *(IterType*)this;
    }

    IterType& next_y( long step ) {
        m_yOffset += step * m_channelMapSizeY;
        return *(IterType*)this;
    }

    IterType& prev_y( long step ) {
        m_yOffset -= step * m_channelMapSizeY;
        return *(IterType*)this;
    }

    IterType& next_z( long step ) {
        m_zOffset += step * m_channelMapSizeZ;
        return *(IterType*)this;
    }

    IterType& prev_z( long step ) {
        m_zOffset -= step * m_channelMapSizeZ;
        return *(IterType*)this;
    }

    IterType& go_x( long x ) {
        m_xOffset = x * m_channelMapSizeX;
        return *(IterType*)this;
    }

    IterType& go_y( long y ) {
        m_yOffset = y * m_channelMapSizeY;
        return *(IterType*)this;
    }

    IterType& go_z( long z ) {
        m_zOffset = z * m_channelMapSizeZ;
        return *(IterType*)this;
    }

    IterType& go( long x = 0, long y = 0, long z = 0 ) {
        m_xOffset = x * m_channelMapSizeX;
        m_yOffset = y * m_channelMapSizeY;
        m_zOffset = z * m_channelMapSizeZ;
        return *(IterType*)this;
    }

    long get_offset() const { return m_xOffset + m_yOffset + m_zOffset; }

  protected:
    iterator_3d_base( long channelMapSizeX, long channelMapSizeY, long channelMapSizeZ, long xOffset = 0,
                      long yOffset = 0, long zOffset = 0 )
        : m_channelMapSizeX( channelMapSizeX )
        , m_channelMapSizeY( channelMapSizeY )
        , m_channelMapSizeZ( channelMapSizeZ )
        , m_xOffset( xOffset )
        , m_yOffset( yOffset )
        , m_zOffset( zOffset ) {}

    long m_channelMapSizeX;
    long m_channelMapSizeY;
    long m_channelMapSizeZ;

    long m_xOffset;
    long m_yOffset;
    long m_zOffset;
};

/**
 * A fast access iterator, used when the type of data in the buffer is known.
 */
template <class DataType>
class fast_channel_buffer_iterator : public iterator_3d_base<fast_channel_buffer_iterator<DataType>> {

    friend class channel_buffer;
    char* m_begin;

    fast_channel_buffer_iterator( char* begin, long channelMapSizeX, long channelMapSizeY, long channelMapSizeZ )
        : iterator_3d_base<fast_channel_buffer_iterator<DataType>>( channelMapSizeX, channelMapSizeY,
                                                                    channelMapSizeZ ) {
        m_begin = begin;
    }

  public:
    typedef DataType Type;

    DataType& operator*() { return *(DataType*)( m_begin + this->get_offset() ); }

    DataType get() const { return *(DataType*)( m_begin + this->get_offset() ); }

    DataType get( long x ) const { return *(DataType*)( m_begin + x * this->m_channelMapSizeX ); }

    DataType get( long x, long y ) const {
        return *(DataType*)( m_begin + y * this->m_channelMapSizeY + x * this->m_channelMapSizeX );
    }

    DataType get( long x, long y, long z ) const {
        return *(DataType*)( m_begin + z * this->m_channelMapSizeZ + y * this->m_channelMapSizeY +
                             x * this->m_channelMapSizeX );
    }

    void set( long x, const DataType& data ) { *(DataType*)( m_begin + x * this->m_channelMapSizeX ) = data; }

    void set( long x, long y, const DataType& data ) {
        *(DataType*)( m_begin + y * this->m_channelMapSizeY + x * this->m_channelMapSizeX ) = data;
    }

    void set( long x, long y, long z, const DataType& data ) {
        *(DataType*)( m_begin + z * this->m_channelMapSizeZ + y * this->m_channelMapSizeY +
                      x * this->m_channelMapSizeX ) = data;
    }

    void set( const DataType& data ) { *(DataType*)( m_begin + this->get_offset() ) = data; }

    bool operator==( const fast_channel_buffer_iterator& rhs ) const {
        return ( m_begin + this->get_offset() ) == ( rhs.m_begin + rhs.get_offset() );
    }

    bool operator!=( const fast_channel_buffer_iterator& rhs ) const {
        return ( m_begin + this->get_offset() ) != ( rhs.m_begin + rhs.get_offset() );
    }
};

/**
 * A fast access iterator to a const buffer, used when the type of data in the buffer is known.
 */
template <class DataType>
class const_fast_channel_buffer_iterator : public iterator_3d_base<const_fast_channel_buffer_iterator<DataType>> {

    friend class channel_buffer;
    const char* m_begin;

    const_fast_channel_buffer_iterator( const char* begin, long channelMapSizeX, long channelMapSizeY,
                                        long channelMapSizeZ )
        : iterator_3d_base<const_fast_channel_buffer_iterator<DataType>>( channelMapSizeX, channelMapSizeY,
                                                                          channelMapSizeZ ) {
        m_begin = begin;
    }

  public:
    typedef DataType Type;

    const DataType& operator*() const { return *(DataType*)( m_begin + this->get_offset() ); }

    DataType get( long x ) const { return *(DataType*)( m_begin + x * this->m_channelMapSizeX ); }

    DataType get( long x, long y ) const {
        return *(DataType*)( m_begin + y * this->m_channelMapSizeY + x * this->m_channelMapSizeX );
    }

    DataType get( long x, long y, long z ) const {
        return *(DataType*)( m_begin + z * this->m_channelMapSizeZ + y * this->m_channelMapSizeY +
                             x * this->m_channelMapSizeX );
    }

    DataType get() const { return *(DataType*)( m_begin + this->get_offset() ); }

    bool operator==( const const_fast_channel_buffer_iterator& rhs ) const {
        return ( m_begin + this->get_offset() ) == ( rhs.m_begin + rhs.get_offset() );
    }

    bool operator!=( const const_fast_channel_buffer_iterator& rhs ) const {
        return ( m_begin + this->get_offset() ) != ( rhs.m_begin + rhs.get_offset() );
    }
};

/**
 * A generic access iterator, used when the type of data in the buffer not specifically known,
 * but the user still wants access to a specific data type.
 */
template <class DataType>
class generic_channel_buffer_iterator : public iterator_3d_base<generic_channel_buffer_iterator<DataType>> {

    friend class channel_buffer;
    char* m_begin;
    void ( *m_convertToExtern )( void*, void* );
    void ( *m_convertToIntern )( void*, void* );

    generic_channel_buffer_iterator( char* begin, long channelMapSizeX, long channelMapSizeY, long channelMapSizeZ,
                                     void ( *convertToExtern )( void*, void* ),
                                     void ( *convertToIntern )( void*, void* ) )
        : iterator_3d_base<generic_channel_buffer_iterator<DataType>>( channelMapSizeX, channelMapSizeY,
                                                                       channelMapSizeZ ) {
        m_begin = begin;
        m_convertToExtern = convertToExtern;
        m_convertToIntern = convertToIntern;
    }

  public:
    typedef DataType Type;

    DataType get( long x ) const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + x * this->m_channelMapSizeX ), (void*)&data );
        return data;
    }

    DataType get( long x, long y ) const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + y * this->m_channelMapSizeY + x * this->m_channelMapSizeX ),
                                (void*)&data );
        return data;
    }

    DataType get( long x, long y, long z ) const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + z * this->m_channelMapSizeZ + y * this->m_channelMapSizeY +
                                         x * this->m_channelMapSizeX ),
                                (void*)&data );
        return data;
    }

    DataType get() const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + this->get_offset() ), (void*)&data );
        return data;
    }

    void set( long x, const DataType& data ) {
        ( *m_convertToIntern )( (void*)&data, (void*)( m_begin + x * this->m_channelMapSizeX ) );
    }

    void set( long x, long y, const DataType& data ) {
        ( *m_convertToIntern )( (void*)&data,
                                (void*)( m_begin + y * this->m_channelMapSizeY + x * this->m_channelMapSizeX ) );
    }

    void set( long x, long y, long z, const DataType& data ) {
        ( *m_convertToIntern )( (void*)&data, (void*)( m_begin + z * this->m_channelMapSizeZ +
                                                       y * this->m_channelMapSizeY + x * this->m_channelMapSizeX ) );
    }

    void set( const DataType& data ) {
        ( *m_convertToIntern )( (void*)&data, (void*)( m_begin + this->get_offset() ) );
    }

    bool operator==( const generic_channel_buffer_iterator& rhs ) const {
        return ( m_begin + this->get_offset() ) == ( rhs.m_begin + rhs.get_offset() );
    }

    bool operator!=( const generic_channel_buffer_iterator& rhs ) const {
        return ( m_begin + this->get_offset() ) != ( rhs.m_begin + rhs.get_offset() );
    }
};

/**
 * A generic access iterator to a const buffer, used when the type of data in the buffer not specifically known,
 * but the user still wants access to a specific data type.
 */
template <class DataType>
class const_generic_channel_buffer_iterator : public iterator_3d_base<const_generic_channel_buffer_iterator<DataType>> {

    friend class channel_buffer;
    const char* m_begin;
    void ( *m_convertToExtern )( void*, void* );

    const_generic_channel_buffer_iterator( const char* begin, long channelMapSizeX, long channelMapSizeY,
                                           long channelMapSizeZ, void ( *convertToExtern )( void*, void* ) )
        : iterator_3d_base<const_generic_channel_buffer_iterator<DataType>>( channelMapSizeX, channelMapSizeY,
                                                                             channelMapSizeZ ) {
        m_begin = begin;
        m_convertToExtern = convertToExtern;
    }

  public:
    typedef DataType Type;

    DataType get( long x ) const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + x * this->m_channelMapSizeX ), (void*)&data );
        return data;
    }

    DataType get( long x, long y ) const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + y * this->m_channelMapSizeY + x * this->m_channelMapSizeX ),
                                (void*)&data );
        return data;
    }

    DataType get( long x, long y, long z ) const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + z * this->m_channelMapSizeZ + y * this->m_channelMapSizeY +
                                         x * this->m_channelMapSizeX ),
                                (void*)&data );
        return data;
    }

    DataType get() const {
        DataType data;
        ( *m_convertToExtern )( (void*)( m_begin + this->get_offset() ), (void*)&data );
        return data;
    }

    void set( const DataType& data ) { ( m_convertToExtern )( (void*)&data, (void*)( m_begin + this->get_offset() ) ); }

    bool operator==( const generic_channel_buffer_iterator<DataType>& rhs ) const {
        return ( m_begin + this->get_offset() ) == ( rhs.m_begin + rhs.get_offset() );
    }

    bool operator!=( const generic_channel_buffer_iterator<DataType>& rhs ) const {
        return ( m_begin + this->get_offset() ) != ( rhs.m_begin + rhs.get_offset() );
    }
};

} // namespace channels
} // namespace frantic
