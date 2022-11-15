// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#define _VARIADIC_MAX 10
#define GTEST_HAS_TR1_TUPLE 0

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "targetver.h"

#ifdef _WIN32
#include <tchar.h>

#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#include <algorithm>
#include <limits>
#include <memory>
#include <stdio.h>
#include <vector>

#include "gtest/gtest.h"

// Eigen SIMD is ICEing on MSVC 2008
#if defined( _MSC_VER ) && _MSC_VER < 1600
#define EIGEN_DONT_VECTORIZE
#endif
