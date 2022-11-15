// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/math/hash.hpp>

#if defined( _MSC_VER )
#pragma warning( disable : 4244 4800 4127 4101 4245 )
#endif

using namespace std;
using namespace boost;

// Include Bob Jenkins' hash function
extern "C" {
#include <frantic/math/bob_jenkins_hash.h>
}
