// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

// <cassert> (and therefire <assert.h>) are unique in that they are designed to be included multiple times and can
// change the meaning of the assert macro depending on the state of the NDEBUG macro.
//
// This file will temporarily disable NDEBUG so that asserts are enabled (if FRANTIC_HYBRID is defined at the time of
// inclusion)

#ifdef FRANTIC_HYBRID
#ifdef NDEBUG
#undef NDEBUG
#define FRANTIC_HYBRID_NDEBUG_TRUE
#endif

#include <cassert>

#ifdef FRANTIC_HYBRID_NDEBUG_TRUE
#undef FRANTIC_HYBRID_NDEBUG_TRUE
#define NDEBUG
#endif
#endif
