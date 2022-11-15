// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/locale/locale.hpp>

#if defined( __APPLE__ )
#include <xlocale.h>
#endif

#include <boost/noncopyable.hpp>
#include <boost/thread/once.hpp>

#include <boost/math/special_functions/nonfinite_num_facets.hpp>

namespace frantic {
namespace locale {

// Windows-only, because we're currently using an inline function instead
// on Linux and Mac OS X
#ifdef _WIN32

namespace {

#ifdef _WIN32
typedef _locale_t locale_type;
#else
typedef locale_t locale_type;
#endif

/**
 * RAII wrapper for a locale.
 */
class scoped_locale : boost::noncopyable {
  public:
    scoped_locale()
        : m_locale( 0 ) {}

    scoped_locale( locale_type loc )
        : m_locale( loc ) {}

    virtual ~scoped_locale() { close(); }

    inline operator locale_type() { return m_locale; }

  private:
    locale_type m_locale;

    void close() {
        if( m_locale ) {
#ifdef _WIN32
            _free_locale( m_locale );
#else
            freelocale( m_locale );
#endif

            m_locale = 0;
        }
    }
};

/**
 * Return a POSIX / C standard locale.
 *
 * The "C" locale is useful for us because it always uses a '.' decimal
 * separator.
 */
inline locale_type create_c_locale() {
#if defined( _WIN32 )
    _locale_t loc = _create_locale( LC_ALL, "C" );
#elif defined( __APPLE__ )
    locale_t loc = newlocale( LC_ALL_MASK, NULL, NULL );
#else
    locale_t loc = newlocale( LC_ALL, "C", NULL );
#endif

    if( !loc ) {
        throw std::runtime_error( "create_c_locale Error: unable to create C locale" );
    }

    return loc;
}

/**
 * RAII wrapper for the "C" locale.
 */
class scoped_c_locale : public scoped_locale {
  public:
    scoped_c_locale()
        : scoped_locale( create_c_locale() ) {}
};

class c_locale_singleton {
  public:
    static inline locale_type get() {
        // serialize the first call to init
        boost::call_once( s_flag, &init );
        return init();
    }

  private:
    static boost::once_flag s_flag;

    static inline scoped_c_locale& init() {
        static scoped_c_locale s_loc;
        return s_loc;
    }
};

boost::once_flag c_locale_singleton::s_flag;

} // anonymous namespace

double strtod_c( const char* str, char** endptr ) {
#if defined( _WIN32 )
    return _strtod_l( str, endptr, c_locale_singleton::get() );
#elif defined( __APPLE__ )
    // From the xlocale(3) manpage:
    // "If a NULL locale_t is passed, the C locale will be used."
    return strtod_l( str, endptr, NULL );
#else
    return strtod_l( str, endptr, c_locale_singleton::get() );
#endif
}

int fprintf_c( FILE* stream, const char* format, ... ) {
    locale_type loc( c_locale_singleton::get() );

    va_list varargs;
    va_start( varargs, format );

    const int result = _vfprintf_l( stream, format, loc, varargs );

    va_end( varargs );

    return result;
}

#endif // #ifdef _WIN32

namespace {

template <class CharType, class FloatType>
std::basic_string<CharType> to_string_c_impl( FloatType value ) {
    std::basic_ostringstream<CharType> ss;
    // nonfinite_num_put for a standard representation of non-finite numbers,
    // for example, "inf" instead of the "1.#INF" in Visual Studio
    ss.imbue( std::locale( std::locale::classic(), new boost::math::nonfinite_num_put<CharType> ) );
    ss << value;
    return ss.str();
}

} // anonymous namespace

std::string to_string_c( double value ) { return to_string_c_impl<char>( value ); }

set_locale_in_scope::set_locale_in_scope( const char* locale ) {
    const char* oldLocale = setlocale( LC_ALL, NULL );
    if( !oldLocale ) {
        throw std::runtime_error( "set_locale_in_scope Error: unable to get current locale" );
    }
    m_oldLocale = oldLocale;

    const char* newLocale = setlocale( LC_ALL, locale );
    if( !newLocale ) {
        throw std::runtime_error( "set_locale_in_scope Error: unable to set locale to "
                                  "\"" +
                                  std::string( locale ) + "\"" );
    }
}

set_locale_in_scope::~set_locale_in_scope() { setlocale( LC_ALL, m_oldLocale.c_str() ); }

} // namespace locale
} // namespace frantic
