// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics2d/framebuffer.hpp>
#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace logging {

// This adds framebuffer update functions to the progress logger.
class render_progress_logger : public progress_logger {
  public:
    render_progress_logger( float progressStart = 0, float progressEnd = 100 )
        : progress_logger( progressStart, progressEnd ) {}
    virtual ~render_progress_logger() {}

    // Frame buffer update functions
    virtual void update_frame_buffer( frantic::graphics2d::framebuffer<frantic::graphics::color6f>& buffer ) = 0;
    // TODO: Add variants for updating individual rectangle regions, buckets, etc.
};

// A null logger for situations where logging is just not required.
class null_render_progress_logger : public render_progress_logger {
    void set_title( const frantic::tstring& /*title*/ ) {}

    void update_progress( long long /*completed*/, long long /*maximum*/ ) {}
    void update_progress( float /*percent*/ ) {}
    void update_frame_buffer( frantic::graphics2d::framebuffer<frantic::graphics::color6f>& /*buffer*/ ) {}
};

} // namespace logging
} // namespace frantic
