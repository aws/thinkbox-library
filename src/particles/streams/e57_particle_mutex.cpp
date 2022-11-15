// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/thread/mutex.hpp>
#include <frantic/particles/streams/e57_particle_mutex.hpp>

boost::mutex frantic::particles::streams::g_e57ImageFileConstructorMutex;
