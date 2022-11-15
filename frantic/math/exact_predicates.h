// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#define FRANTIC_EXACT_PREDICATES

namespace frantic {
namespace math {

extern "C" {
// Any program seeking to use the exact arithmetic predicates must call this function to initialize some
// eror bounds information.
extern void exactinit();

extern float insphere( float* pa, float* pb, float* pc, float* pd, float* pe );
extern float incircle( float* pa, float* pb, float* pc, float* pd );
extern float orient3d( float* pa, float* pb, float* pc, float* pd );
extern float orient2d( float* pa, float* pb, float* pc );

/*extern double insphere_double(double* pa, double* pb, double* pc, double* pd, double* pe);
extern double incircle_double(double* pa, double* pb, double* pc, double* pd);
extern double orient3d_double(double* pa, double* pb, double* pc, double* pd);
extern double orient2d_double(double* pa, double* pb, double* pc);*/
}

} // namespace math
} // namespace frantic
