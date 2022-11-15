// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef _WIN32
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#endif

namespace frantic {
namespace locale {

/**
 * Convert string to double, using the "C" locale.
 *
 * The "C" locale is useful here because it always uses a '.' decimal
 * separator.  Other locales may use a different decimal separator.
 * For example, the German locale uses a ',', which causes problems when
 * parsing floating point numbers in a comma-separated value (CSV) file.
 */
#if defined( _WIN32 )
double strtod_c( const char* str, char** endptr );
#else
inline double strtod_c( const char* str, char** endptr ) {
#if defined( __APPLE__ )
    // From the xlocale(3) manpage:
    // "If a NULL locale_t is passed, the C locale will be used."
    return strtod_l( str, endptr, NULL );
#else
    locale_t loc = newlocale( LC_ALL, "C", NULL );
    const double result = strtod_l( str, endptr, loc );
    freelocale( loc );
    return result;
#endif
}

#endif

/**
 * Write formatted string to a file, using the "C" locale.
 */
#if defined( _WIN32 )
int fprintf_c( FILE* stream, const char* format, ... );
#else
inline int fprintf_c( FILE* stream, const char* format, ... ) {
#if defined( __APPLE__ )
    va_list varargs;
    va_start( varargs, format );

    const int result = vfprintf_l( stream, NULL, format, varargs );

    va_end( varargs );

    return result;
#else
    locale_t loc = newlocale( LC_ALL, "C", NULL );
    locale_t oldLoc = uselocale( loc );

    va_list varargs;
    va_start( varargs, format );

    const int result = vfprintf( stream, format, varargs );

    va_end( varargs );

    uselocale( oldLoc );
    freelocale( loc );

    return result;
#endif
}
#endif

/**
 *  Convert the value to a string, using the "C" locale.
 */
std::string to_string_c( double value );

/**
 * Set the program locale by calling setlocale().
 * This operation is typically per-process, and not thread safe.
 *
 * @note This is intended for use in test suites only.
 */
class set_locale_in_scope {
  public:
    set_locale_in_scope( const char* locale );
    ~set_locale_in_scope();

  private:
    std::string m_oldLocale;
};

} // namespace locale
} // namespace frantic
