// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/concatenated_particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

concatenated_particle_istream::concatenated_particle_istream(
    const std::vector<boost::shared_ptr<particle_istream>>& pins ) {
    // Using a member initializer list here seems to induce a crash on GCC 5.2.1 on Ubuntu 15.10 under -O2
    m_delegates = pins;
    if( m_delegates.empty() )
        throw std::runtime_error(
            "concatenated_particle_istream() - The provided particle_istream array was empty.  It "
            "should contain at least one stream." );

    init_streams();
}
} // namespace streams
} // namespace particles
} // namespace frantic
