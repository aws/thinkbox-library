// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/thread/mutex.hpp>

namespace frantic {
namespace particles {
namespace streams {
/// Mutex to ensure that only one thread calls the e57::ImageFile constructor at a time. The ImageFile constructor
/// calls xercesc::XMLPlatforumUtils::Initialize() and Terminate(), which are not thread safe, and can cause a crash
/// if multiple threads create an e57_particle_ostream simultaneously.
extern boost::mutex g_e57ImageFileConstructorMutex;
} // namespace streams
} // namespace particles
} // namespace frantic
