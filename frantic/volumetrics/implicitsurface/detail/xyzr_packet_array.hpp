// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_accessor.hpp>

#include <frantic/simd/float_v.hpp>

namespace frantic {
namespace volumetrics {
namespace implicitsurface {
namespace detail {

struct xyzr_packet {
    frantic::graphics::vector3t<frantic::simd::float_v> position;
    frantic::simd::float_v radius;
};

/**
 *  A class to help working with particle data in an "array of structures of
 * arrays" layout.
 */
class xyzr_packet_array {
  public:
    xyzr_packet_array();

    void clear();

    void load( const char* const* particlesBegin, const char* const* particlesEnd,
               const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAcc,
               const frantic::channels::channel_accessor<float>& radiusAcc );

    const std::vector<xyzr_packet>& get_particles() const { return m_particles; };

    std::size_t get_particle_count() const { return m_particleCount; }

    /**
     * @return the number of packets in which all entries are populated.
     */
    std::size_t get_filled_packet_count() const { return m_filledPacketCount; }

    /**
     * @return the total number of particles contained in filled packets.
     */
    std::size_t get_filled_particle_count() const { return m_filledParticleCount; }

    /**
     * @return true if at least one entry in the last packet is not populated.
     */
    bool has_remainder() const { return m_remainderParticleCount > 0; }

    /**
     * @return true if at least one entry in the last packet is not populated.
     */
    std::size_t get_remainder_particle_count() const { return m_remainderParticleCount; }

    /**
     * @return a mask that is set for all populated entries in the remainder
     *         packet, and cleared for all unpopulated entries.
     */
    const frantic::simd::float_v& get_remainder_mask() const { return m_remainderParticleMask; }

  private:
    void resize( std::size_t particleCount );

#ifdef FRANTIC_HAS_SSE2
    void load_adjacent_internal( const char* const* particlesBegin,
                                 const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAcc );
#endif

    std::vector<xyzr_packet> m_particles;

    std::size_t m_particleCount;

    std::size_t m_filledPacketCount;

    std::size_t m_filledParticleCount;

    std::size_t m_remainderParticleCount;

    frantic::simd::float_v m_remainderParticleMask;
};

} // namespace detail
} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic
