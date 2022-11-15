// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <zlib.h>

#include <boost/tuple/tuple.hpp>
#include <frantic/files/pipelined_ifstream.hpp>
#include <frantic/files/pipelined_ofstream.hpp>

#include <frantic/misc/exception_stream.hpp>

namespace frantic {
namespace files {

/**
 * This stream wraps an std::istream with zlib compression.
 */
class zlib_deflate_ostream {
  private:
    int m_bufferSize;
    char* m_buffer;
    std::ostream* m_out;
    frantic::tstring m_streamName;
    z_stream m_zstream;

    bool m_bOpen;

    void flush() {
        if( m_zstream.next_out > reinterpret_cast<unsigned char*>( m_buffer ) ) {
            m_out->write( m_buffer, m_zstream.next_out - reinterpret_cast<unsigned char*>( m_buffer ) );
            m_zstream.avail_out = m_bufferSize;
            m_zstream.next_out = reinterpret_cast<unsigned char*>( m_buffer );
            if( !m_out && !std::uncaught_exceptions() )
                throw std::runtime_error(
                    "zlib_deflate_ostream.flush: Failed to write compressed data to the output stream \"" +
                    frantic::strings::to_string( m_streamName ) + "\"." );
        }
    }

    // Disable copy construction
    zlib_deflate_ostream( const zlib_deflate_ostream& ); // not implemented

    // Disable assignment
    zlib_deflate_ostream& operator=( const zlib_deflate_ostream& ); // not implemented

  public:
    zlib_deflate_ostream( int bufferSize = ( 1 << 19 ) )
        : m_bufferSize( bufferSize )
        , m_bOpen( false ) {
        m_buffer = new char[bufferSize];
        m_streamName = _T("<unknown>");
    }
    zlib_deflate_ostream( std::ostream& out, const frantic::tstring& streamName,
                          int compressionLevel = Z_DEFAULT_COMPRESSION, int bufferSize = 32768 )
        : m_bufferSize( bufferSize )
        , m_bOpen( false ) {
        m_buffer = new char[bufferSize];
        open( out, streamName, compressionLevel );
    }

    ~zlib_deflate_ostream() {
        close();
        delete[] m_buffer;
    }

    void open( std::ostream& out, const frantic::tstring& streamName, int compressionLevel = Z_DEFAULT_COMPRESSION ) {
        m_streamName = streamName;

        if( !out )
            throw std::runtime_error( "zlib_deflate_ostream.open: Received an invalid std::ostream + \"" +
                                      frantic::strings::to_string( m_streamName ) + "\"." );

        m_out = &out;

        memset( &m_zstream, 0, sizeof( z_stream ) );
        if( Z_OK != deflateInit( &m_zstream, compressionLevel ) )
            throw std::runtime_error(
                "zlib_deflate_ostream.open: Unable to initialize a zlib deflate stream for the output stream \"" +
                frantic::strings::to_string( m_streamName ) + "\"." );

        m_zstream.avail_out = m_bufferSize;
        m_zstream.next_out = reinterpret_cast<unsigned char*>( m_buffer );
        m_bOpen = true;
    }

    void close() {
        if( m_bOpen ) {
            // Write out all the rest of the stream data
            while( Z_STREAM_END != deflate( &m_zstream, Z_FINISH ) )
                flush();
            flush();

            m_bOpen = false;
            m_out = NULL;

            // If this function fails, there was probably some corrupted state somewhere, so we should throw.
            int ret = deflateEnd( &m_zstream );
            if( !std::uncaught_exceptions() ) {
                if( ret == Z_STREAM_ERROR )
                    throw std::runtime_error(
                        "zlib_deflate_ostream.close: An error occured when shutting down a zlib deflate "
                        "stream for output stream \"" +
                        frantic::strings::to_string( m_streamName ) + "\", got error code Z_STREAM_ERROR." );
                if( ret == Z_DATA_ERROR )
                    throw std::runtime_error(
                        "zlib_deflate_ostream.close: An error occured when shutting down a zlib deflate "
                        "stream for output stream \"" +
                        frantic::strings::to_string( m_streamName ) + "\", got error code Z_DATA_ERROR" );
                if( ret != Z_OK )
                    throw std::runtime_error(
                        "zlib_deflate_ostream.close: An error occured when shutting down a zlib deflate "
                        "stream for output stream \"" +
                        frantic::strings::to_string( m_streamName ) + "\", got error code " +
                        boost::lexical_cast<std::string>( ret ) );
            }
        }
    }

    void write( const char* data, std::size_t dataSize ) {
        if( !m_bOpen )
            throw std::runtime_error( "zlib_deflate_ostream.write: Tried to write data to output stream \"" +
                                      frantic::strings::to_string( m_streamName ) + "\" after the stream was closed." );

        m_zstream.avail_in = (unsigned)dataSize;
        // NOTE: It should be ok to const_cast here, because this is what zlib is reading from...  If zlib modifies
        //       the values temporarily, though, we should make a copy of the input first.
        m_zstream.next_in = reinterpret_cast<unsigned char*>( const_cast<char*>( data ) );

        for( ;; ) {
            int ret = deflate( &m_zstream, Z_NO_FLUSH );
            if( ret == Z_STREAM_ERROR )
                throw std::runtime_error(
                    "zlib_deflate_ostream.write: Failed to write data to the zlib deflate stream, got "
                    "error code Z_STREAM_ERROR." );

            if( m_zstream.avail_out == 0 ) { // Did we fill the output buffer completely? If so flush and loop.
                flush();
            } else {
                break;
            }
        }

        return;
    }
};

/**
 * This stream wraps an std::istream with zlib decompression.
 */
class zlib_inflate_istream {
    int m_bufferSize;
    char* m_buffer;
    std::istream* m_in;
    frantic::tstring m_streamName;
    z_stream m_zstream;
    bool m_bOpen;

    // Disable copy construction
    zlib_inflate_istream( const zlib_inflate_istream& ); // not implemented

    // Disable assignment
    zlib_inflate_istream& operator=( const zlib_inflate_istream& ); // not implemented

  public:
    zlib_inflate_istream( int bufferSize = 32768 )
        : m_bufferSize( bufferSize )
        , m_streamName( _T("<unknown>") )
        , m_bOpen( false ) {
        m_buffer = new char[bufferSize];
    }

    zlib_inflate_istream( std::istream& in, const frantic::tstring& streamName, int bufferSize = 32768 )
        : m_bufferSize( bufferSize )
        , m_streamName( _T("<unknown>") )
        , m_bOpen( false ) {
        m_buffer = new char[bufferSize];
        open( in, streamName );
    }

    ~zlib_inflate_istream() {
        close();
        delete[] m_buffer;
    }

    void close() {
        if( m_bOpen ) {
            inflateEnd( &m_zstream );
            m_bOpen = false;
        }
    }

    void open( std::istream& in, const frantic::tstring& streamName ) {
        m_streamName = streamName;
        if( !in )
            throw std::runtime_error( "zlib_inflate_istream.open: Received a bad std::istream for input stream \"" +
                                      frantic::strings::to_string( m_streamName ) + "\"." );

        m_in = &in;

        memset( &m_zstream, 0, sizeof( m_zstream ) );
        if( Z_OK != inflateInit( &m_zstream ) )
            throw std::runtime_error(
                "zlib_inflate_istream.open: Unable to initialize a zlib inflate stream for input stream \"" +
                frantic::strings::to_string( m_streamName ) + "\"." );

        m_bOpen = true;
    }

    void read( char* data, int nBytes ) {
        if( !m_bOpen )
            throw std::runtime_error( "zlib_inflate_istream.read: Attempted to read data from stream \"" +
                                      frantic::strings::to_string( m_streamName ) + "\" after it was closed." );

        m_zstream.avail_out = nBytes;
        m_zstream.next_out = reinterpret_cast<unsigned char*>( data );

        do {
            if( m_zstream.avail_in == 0 ) {
                m_in->read( m_buffer, m_bufferSize );
                m_zstream.avail_in = static_cast<uInt>( m_in->gcount() );
                m_zstream.next_in = reinterpret_cast<unsigned char*>( m_buffer );
            }

            int ret = inflate( &m_zstream, Z_SYNC_FLUSH );
            if( Z_OK != ret && Z_STREAM_END != ret )
                throw exception_stream() << "zlib_inflate_istream.read: ZLib inflate call for stream \"" +
                                                frantic::strings::to_string( m_streamName ) +
                                                "\" failed with error code: "
                                         << ret;

        } while( m_zstream.avail_out != 0 );
    }
};

struct zlib_deflate_cstdio_traits {
    static std::string type_str() { return std::string( "zlib deflate" ); }
    static int init_stream( z_stream* strm, int compressionLevel ) { return deflateInit( strm, compressionLevel ); }
};

struct zlib_gzip_cstdio_traits {
    static std::string type_str() { return std::string( "gzip" ); }
    static int init_stream( z_stream* strm, int compressionLevel ) {
        return deflateInit2( strm, compressionLevel, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY );
    }
};

template <class ZLibCompressionType>
class zlib_ostream_cstdio_base {
  private:
    int m_bufferSize;
    pipelined_ofstream m_out;
    frantic::tstring m_streamName;
    z_stream m_zstream;

    bool m_bOpen;

    void flush() {
        std::size_t numOut = ( m_out.get_write_buffer_size() - m_zstream.avail_out );
        if( numOut > 0 ) {
            m_out.write( numOut );
            m_zstream.avail_out = static_cast<unsigned int>( m_out.get_write_buffer_size() );
            m_zstream.next_out = reinterpret_cast<unsigned char*>( m_out.get_write_buffer() );
        }
    }

    // Disable copy construction
    zlib_ostream_cstdio_base( const zlib_ostream_cstdio_base& ) {}

    // Disable assignment
    zlib_ostream_cstdio_base& operator=( const zlib_ostream_cstdio_base& ) {}

  public:
    zlib_ostream_cstdio_base()
        : m_streamName( _T("<unknown>") )
        , m_bOpen( false ) {}

    zlib_ostream_cstdio_base( FILE* out, const frantic::tstring& streamName,
                              int compressionLevel = Z_DEFAULT_COMPRESSION, int bufferSize = ( 1 << 20 ) )
        : m_bOpen( false ) {
        open( out, streamName, compressionLevel, bufferSize );
    }

    ~zlib_ostream_cstdio_base() { close(); }

    void open( FILE* out, const frantic::tstring& streamName, int compressionLevel = Z_DEFAULT_COMPRESSION,
               int bufferSize = ( 1 << 20 ) ) {
        m_streamName = streamName;

        if( !out )
            throw std::runtime_error( "zlib_ostream_cstdio_base.open: Received an invalid FILE* \"" +
                                      frantic::strings::to_string( m_streamName ) + "\"." );

        m_out.open( out, bufferSize );

        memset( &m_zstream, 0, sizeof( z_stream ) );
        if( Z_OK != ZLibCompressionType::init_stream( &m_zstream, compressionLevel ) )
            throw std::runtime_error( "zlib_ostream_cstdio_base.open: Unable to initialize a " +
                                      ZLibCompressionType::type_str() + " stream for the output stream \"" +
                                      frantic::strings::to_string( m_streamName ) + "\"." );

        m_zstream.avail_out = static_cast<unsigned int>( m_out.get_write_buffer_size() );
        m_zstream.next_out = reinterpret_cast<unsigned char*>( m_out.get_write_buffer() );
        m_bOpen = true;
    }

    void close() {
        if( m_bOpen ) {
            // Write out all the rest of the stream data
            while( Z_STREAM_END != deflate( &m_zstream, Z_FINISH ) )
                flush();
            flush();

            m_bOpen = false;
            m_out.close();

            // If this function fails, there was probably some corrupted state somewhere, so we should throw.
            int ret = deflateEnd( &m_zstream );
            if( !std::uncaught_exceptions() ) {
                if( ret == Z_STREAM_ERROR )
                    throw std::runtime_error( "zlib_ostream_cstdio_base.close: An error occured when shutting down a " +
                                              ZLibCompressionType::type_str() + " stream for output stream \"" +
                                              frantic::strings::to_string( m_streamName ) +
                                              "\", got error code Z_STREAM_ERROR." );
                if( ret == Z_DATA_ERROR )
                    throw std::runtime_error( "zlib_ostream_cstdio.close: An error occured when shutting down a " +
                                              ZLibCompressionType::type_str() + " stream for output stream \"" +
                                              frantic::strings::to_string( m_streamName ) +
                                              "\", got error code Z_DATA_ERROR" );
                if( ret != Z_OK )
                    throw std::runtime_error( "zlib_ostream_cstdio.close: An error occured when shutting down a " +
                                              ZLibCompressionType::type_str() + " stream for output stream \"" +
                                              frantic::strings::to_string( m_streamName ) + "\", got error code " +
                                              boost::lexical_cast<std::string>( ret ) );
            }
        }
    }

    void write( const char* data, std::size_t dataSize ) {
        if( !m_bOpen )
            throw std::runtime_error( "zlib_ostream_cstdio_base.write: Tried to write data to output stream \"" +
                                      frantic::strings::to_string( m_streamName ) + "\" after the stream was closed." );

        m_zstream.avail_in = (unsigned)dataSize;
        // NOTE: It should be ok to const_cast here, because this is what zlib is reading from...  If zlib modifies
        //       the values temporarily, though, we should make a copy of the input first.
        m_zstream.next_in = reinterpret_cast<unsigned char*>( const_cast<char*>( data ) );

        for( ;; ) {
            int ret = deflate( &m_zstream, Z_NO_FLUSH );
            if( ret == Z_STREAM_ERROR )
                throw std::runtime_error( "zlib_ostream_cstdio_base.write: Failed to write data to the " +
                                          ZLibCompressionType::type_str() + " stream, got error code Z_STREAM_ERROR." );

            if( m_zstream.avail_out == 0 ) { // Did we fill the output buffer completely? If so flush and loop.
                flush();
            } else {
                break;
            }
        }

        return;
    }
};

typedef zlib_ostream_cstdio_base<zlib_deflate_cstdio_traits> zlib_deflate_ostream_cstdio;
typedef zlib_ostream_cstdio_base<zlib_gzip_cstdio_traits> zlib_gzip_ostream_cstdio;

struct zlib_inflate_cstdio_traits {
    static std::string type_str() { return std::string( "zlib inflate" ); }
    static int init_stream( z_stream& strm ) { return inflateInit( &strm ); }
};

struct zlib_gunzip_cstdio_traits {
    static std::string type_str() { return std::string( "zlib" ); }
    static int init_stream( z_stream& strm ) { return inflateInit2( &strm, 32 | MAX_WBITS ); }
};

template <class ZLibDecompressionType>
class zlib_istream_cstdio_base {
  private:
    int m_bufferSize;
    pipelined_ifstream m_pipe;
    frantic::tstring m_streamName;
    z_stream m_zstream;
    bool m_bOpen;

    // Disable copy construction
    zlib_istream_cstdio_base( const zlib_istream_cstdio_base& ) {}

    // Disable assignment
    zlib_istream_cstdio_base& operator=( const zlib_istream_cstdio_base& ) {}

  public:
    zlib_istream_cstdio_base( int bufferSize = ( 1 << 16 ) )
        : m_bOpen( false )
        , m_bufferSize( bufferSize )
        , m_streamName( _T("<unknown>") ) {}

    zlib_istream_cstdio_base( FILE* in, const frantic::tstring& streamName, int bufferSize = ( 1 << 16 ) )
        : m_bOpen( false )
        , m_bufferSize( bufferSize )
        , m_streamName( streamName ) {
        open( in, streamName );
    }

    ~zlib_istream_cstdio_base() { close(); }

    void close() {
        if( m_bOpen ) {
            inflateEnd( &m_zstream );
            m_pipe.close();
            m_bOpen = false;
        }
    }

    void open( FILE* fin, const frantic::tstring& streamName ) {
        m_streamName = streamName;
        if( !fin )
            throw std::runtime_error( "zlib_istream_cstdio_base.open: Unable to open file: \"" +
                                      frantic::strings::to_string( m_streamName ) + "\"." );

        m_pipe.open( fin, m_bufferSize );

        memset( &m_zstream, 0, sizeof( m_zstream ) );
        if( Z_OK != ZLibDecompressionType::init_stream( m_zstream ) )
            throw std::runtime_error( "zlib_istream_cstdio_base.open: Unable to initialize a " +
                                      ZLibDecompressionType::type_str() + " stream for input stream \"" +
                                      frantic::strings::to_string( m_streamName ) + "\"." );

        m_bOpen = true;
    }

    void read( char* data, int nBytes ) {
        if( !m_bOpen )
            throw std::runtime_error( "zlib_istream_cstdio_base.read: Attempted to read data from stream \"" +
                                      frantic::strings::to_string( m_streamName ) + "\" after it was closed." );

        m_zstream.avail_out = nBytes;
        m_zstream.next_out = reinterpret_cast<unsigned char*>( data );

        do {
            if( m_zstream.avail_in == 0 ) {
                char* pBuffer;
                std::size_t bufferSize;

                boost::tuples::tie( pBuffer, bufferSize ) = m_pipe.read();

                m_zstream.avail_in = static_cast<uInt>( bufferSize );
                m_zstream.next_in = reinterpret_cast<Byte*>( pBuffer );
            }

            int ret = inflate( &m_zstream, Z_SYNC_FLUSH );
            if( Z_STREAM_END == ret ) {
                if( m_zstream.avail_out != 0 )
                    throw std::runtime_error(
                        "zlib_istream_cstdio_base.read: End of stream reached during read from stream\"" +
                        frantic::strings::to_string( m_streamName ) + "\"" );
            } else if( Z_OK != ret )
                throw exception_stream() << "zlib_istream_cstdio_base.read: ZLib inflate call for stream \"" +
                                                frantic::strings::to_string( m_streamName ) +
                                                "\" failed with error code: "
                                         << ret;

        } while( m_zstream.avail_out != 0 );
    }
};

typedef zlib_istream_cstdio_base<zlib_inflate_cstdio_traits> zlib_inflate_istream_cstdio;
typedef zlib_istream_cstdio_base<zlib_gunzip_cstdio_traits> zlib_gzip_istream_cstdio;

} // namespace files
} // namespace frantic
