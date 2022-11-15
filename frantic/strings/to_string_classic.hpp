// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <sstream>

namespace frantic {
namespace strings {

/**
 * boost::lexical_cast() uses the current C++ locale.
 * This is a problem for certain locales that cause integers >999 to
 * have commas inserted into them (eg 1000 -> 1,000).
 * We do not want to set the locale ourselves, since this is not thread safe.
 * This should be used as a safer alternative to boost::lexical_cast()
 * when converting integers to strings.
 */
template <typename T, typename I>
T to_string_classic( I input ) {
    std::basic_ostringstream<typename T::value_type> ss;
    ss.imbue( std::locale::classic() );
    ss << input;
    return ss.str();
}

} // namespace strings
} // namespace frantic
