// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <boost/shared_ptr.hpp>

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <frantic/particles/particle_utilities.hpp>

#include <frantic/channels/channel_map.hpp>

#include <frantic/particles/particle_cursor.hpp>

namespace frantic {
namespace particles {

class particle_step {
    // constant old values
    const frantic::graphics::vector3f& m_oldVelocity;
    const frantic::graphics::vector3f& m_oldPosition;

    int m_startTime;
    int m_endTime;
    float m_timeStep;

    // new values that can be set
    boost::uint32_t& m_flag;
    frantic::graphics::vector3f& m_newVelocity;
    frantic::graphics::vector3f& m_newPosition;

    // Make the assignment operator private so no one can use it.
    particle_step& operator=( const particle_step& ) {}

  public:
    particle_step( const frantic::graphics::vector3f& oldVelocity, const frantic::graphics::vector3f& oldPosition,
                   frantic::graphics::vector3f& newVelocity, frantic::graphics::vector3f& newPosition,
                   boost::uint32_t& flag, int startTime, int endTime, float timeStep )
        : m_oldVelocity( oldVelocity )
        , m_oldPosition( oldPosition )
        , m_newVelocity( newVelocity )
        , m_newPosition( newPosition )
        , m_flag( flag ) {
        m_startTime = startTime;
        m_endTime = endTime;
        m_timeStep = timeStep;
    }

    void set( const frantic::graphics::vector3f& position, const frantic::graphics::vector3f& velocity ) {
        m_newPosition = position;
        m_newVelocity = velocity;
        m_flag = m_flag | frantic::particles::PRT_FLG_MODIFIED; // set the modified flag within the ID number
    }

    const frantic::graphics::vector3f& oldPosition() const { return m_oldPosition; }
    const frantic::graphics::vector3f& oldVelocity() const { return m_oldVelocity; }

    frantic::graphics::vector3f& newPosition() const { return m_newPosition; }
    frantic::graphics::vector3f& newVelocity() const { return m_newVelocity; }

    int start_time() const { return m_startTime; }
    int end_time() const { return m_endTime; }
    float time_step() const { return m_timeStep; }
    boost::uint32_t flag() const { return m_flag; }
};

} // namespace particles
} // namespace frantic
