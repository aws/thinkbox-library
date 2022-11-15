// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#if defined( _MSC_VER )
// Disables deprecated function warnings
#pragma warning( disable : 4996 )

// Disables unreachable code warnings
#pragma warning( disable : 4702 )
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0501
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <winsock2.h>
#endif

// All the standard C++ headers we need
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined( __GNUC__ ) && __GNUC__ < 3
#include <limits.h>
#else
#include <limits>
#endif
//#include <hash_set>

// Because 3ds max irresponsibly defines lower-case preprocessor macros,
// we need to do this here instead of elsewhere.
#include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/integer_fwd.hpp>

#include <boost/lexical_cast.hpp>

#include <boost/smart_ptr.hpp>

// Eigen SIMD is ICEing on MSVC 2008
#if defined( _MSC_VER ) && _MSC_VER < 1600
#define EIGEN_DONT_VECTORIZE
#endif

// Avoid error C2719 (parameter with __declspec(align('16')) won't be aligned)
// when building src/geometry/quadtree.cpp
#if defined( _MSC_VER ) && _MSC_VER < 1700 && defined( _M_IX86 )
#define EIGEN_DONT_ALIGN
#endif

#include <Eigen/Core>
