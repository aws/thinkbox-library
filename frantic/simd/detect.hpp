// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#ifndef FRANTIC_DISABLE_SIMD
#if defined( __SSE2__ ) || defined( _M_X64 ) || ( defined( _M_IX86_FP ) && _M_IX86_FP >= 2 )
#define FRANTIC_HAS_SSE2
#endif
#endif

#ifdef FRANTIC_HAS_SSE2
#define FRANTIC_SIMD_NAMESPACE sse2
#else
#define FRANTIC_SIMD_NAMESPACE scalar
#endif
