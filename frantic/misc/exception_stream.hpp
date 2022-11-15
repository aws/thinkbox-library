// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <boost/shared_ptr.hpp>
#include <sstream>
#include <streambuf>

namespace frantic {

// This specialization of std::basic_stringbuf exposes the underlying char* of the buffer and adds methods
// for manipulating the ending character of the buffer. This is convenient for maintaining the terminating NULL.
template <class CharType, class Traits>
class basic_charbuf : public std::basic_stringbuf<CharType, Traits> {
  private:
    // To stop compiler complaints
    basic_charbuf( const basic_charbuf& buf );
    basic_charbuf& operator=( const basic_charbuf& buf );

  public:
    basic_charbuf() {}

    // Returns ptr to the out buffer
    CharType* c_str() { return std::basic_stringbuf<CharType, Traits>::gptr(); }

    void pop_back() {
        if( std::basic_stringbuf<CharType, Traits>::pptr() > std::basic_stringbuf<CharType, Traits>::pbase() )
            std::basic_stringbuf<CharType, Traits>::pbump( -1 );
    }

    void push_back( CharType c ) { std::basic_stringbuf<CharType, Traits>::sputc( c ); }
};
typedef basic_charbuf<char, std::char_traits<char>> char_buf;

// class throw_exception_stream
//
//  Example Usage:
//    int value = 3;
//    throw exception_stream() << "testing " << value;
//
// This class is a derived exception that provides stream operators in order
// to build the message it contains.
template <class ExceptionClass>
class exception_stream_base : public std::exception {
  private:
    // Stores the buffer we want the stream to use. Make sure not to
    // create multiple buffers when copying exceptions because it is a bad idea
    // to allocate memory after an exception has been thrown.
    boost::shared_ptr<char_buf> m_pBuffer;

    // The stream for this object.
    std::ostream m_stream;

    exception_stream_base& operator=( const exception_stream_base& /*e*/ ) throw() {}

  public:
    exception_stream_base()
        : m_pBuffer( new char_buf )
        , m_stream( m_pBuffer.get() ) // This relies on constructor order. Maybe set this to NULL and use the code below
    {
        // m_stream.rdbuf( m_pBuffer.get() );

        // Mark the end of the buffer w/ a NULL chasracter
        m_pBuffer->push_back( '\0' );
    }

    exception_stream_base( const std::string& msg )
        : m_pBuffer( new char_buf )
        , m_stream( m_pBuffer.get() ) // This relies on constructor order. Maybe set this to NULL and use the code below
    {
        // m_stream.rdbuf( m_pBuffer.get() );

        m_stream << msg;

        // Mark the end of the buffer w/ a NULL chasracter
        m_pBuffer->push_back( '\0' );
    }

    exception_stream_base( const exception_stream_base& e )
        : m_pBuffer( e.m_pBuffer )
        , m_stream( e.m_pBuffer.get() ) {}

    virtual ~exception_stream_base() throw() {}

    const char* what() const throw() { return m_pBuffer->c_str(); }

    template <class Type>
    ExceptionClass& operator<<( const Type& t ) {
        m_pBuffer->pop_back(); // Pop the trailing NULL char out of the stream in order to correctly concatenate
        m_stream << t;
        m_pBuffer->push_back( '\0' ); // Push a trailing NULL onto the stream in order to mark the end of the c_str()

        return *static_cast<ExceptionClass*>( this );
    }
};

class exception_stream : public exception_stream_base<exception_stream> {
  public:
    exception_stream() {}
    virtual ~exception_stream() throw() {}
};

class invalid_particle_file_exception : public exception_stream_base<invalid_particle_file_exception> {
  public:
    invalid_particle_file_exception( const std::string& msg )
        : exception_stream_base<invalid_particle_file_exception>( msg ) {}
    invalid_particle_file_exception() {}
    virtual ~invalid_particle_file_exception() throw() {}
};

} // namespace frantic
