// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <iosfwd>
#include <sstream>

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace tsstream {
#ifdef FRANTIC_USE_WCHAR
typedef std::wostringstream tostringstream;
#else
typedef std::ostringstream tostringstream;
#endif
} // namespace tsstream

using frantic::tsstream::tostringstream;
} // namespace frantic