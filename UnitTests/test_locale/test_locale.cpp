// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/lexical_cast.hpp>

#include <frantic/locale/locale.hpp>

TEST( Locale, SetLocaleInScope ) {
#ifdef _WIN32
#if defined( _MSC_VER ) && _MSC_VER >= 1700
    const char* german = "de-DE";
#else
    const char* german = "German";
#endif
#else
    const char* german = "de_DE.UTF-8";
#endif

    { // test in German locale
        frantic::locale::set_locale_in_scope setLocale( german );

        EXPECT_EQ( "0,5", boost::lexical_cast<std::string>( 0.5 ) );
    }

    { // back to "C" locale
        EXPECT_EQ( "0.5", boost::lexical_cast<std::string>( 0.5 ) );
    }
}
