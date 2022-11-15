// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/math/fractions.hpp>
#include <frantic/math/uint128.hpp>
#include <frantic/math/utils.hpp>
#include <frantic/particles/streams/empty_particle_istream.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <boost/make_shared.hpp>

namespace frantic {
namespace particles {
namespace streams {

class fractional_particle_istream : public delegated_particle_istream {
    boost::int64_t m_numerator, m_denominator, m_accumulator;
    boost::int64_t m_curIndex, m_totalCount;
    bool m_evenly;

    bool is_next_in_set() {
        m_accumulator += m_numerator;
        if( m_accumulator >= m_denominator ) {
            m_accumulator -= m_denominator;
            return true;
        }
        return false;
    }

  public:
    fractional_particle_istream( boost::shared_ptr<particle_istream> stream, double fraction,
                                 boost::int64_t limit = 2000000000, bool evenlyDistribute = true )
        : delegated_particle_istream( stream )
        , m_evenly( evenlyDistribute ) {
        m_curIndex = m_delegate->particle_index();
        m_totalCount = m_delegate->particle_count();

        fraction = frantic::math::clamp( fraction, 0.0, 1.0 );

        std::pair<boost::int64_t, boost::int64_t> rational = frantic::math::get_rational_representation( fraction );
        m_numerator = rational.first;
        m_denominator = rational.second;
        m_accumulator = 0;

        if( m_totalCount < 0 ) {
            m_evenly = true;
            m_totalCount = limit;
        } else if( m_totalCount > 0 ) {
            using frantic::math::uint128;
            m_totalCount = std::min(
                limit,
                ( uint128( m_totalCount ) * uint128( m_numerator ) / m_denominator ).to_integral<boost::int64_t>() );

            // Force a single particle instead of rounding down to zero.
            if( m_totalCount == 0 )
                ++m_totalCount;
        }

        // Set up the accumulator such that it always accepts the first particle so we can always grab at least a single
        // particle.
        m_accumulator = std::max<boost::int64_t>( m_denominator - m_numerator, 0 );
    }

    virtual ~fractional_particle_istream() {}

    boost::int64_t particle_count() const { return m_delegate->particle_count() >= 0 ? m_totalCount : -1; }
    boost::int64_t particle_count_left() const {
        return m_delegate->particle_count() >= 0 ? m_totalCount - m_curIndex - 1 : -1;
    }
    boost::int64_t particle_index() const { return m_curIndex; }

    virtual boost::int64_t particle_count_guess() const {
        using frantic::math::uint128;
        return ( uint128( std::max( static_cast<boost::int64_t>( 0 ), m_delegate->particle_count_guess() ) ) *
                 uint128( m_numerator ) / m_denominator )
            .to_integral<boost::int64_t>();
    }

    bool get_particle( char* buffer ) {
        if( ++m_curIndex >= m_totalCount )
            return false;

        if( !m_evenly ) {
            if( !m_delegate->get_particle( buffer ) )
                return false;
        } else {
            do {
                if( !m_delegate->get_particle( buffer ) )
                    return false;
            } while( !is_next_in_set() );
        }

        return true;
    }

    bool get_particles( char* particleBuffer, std::size_t& numParticles ) {
        if( m_curIndex >= m_totalCount ) {
            numParticles = 0;
            return false;
        }

        if( !m_evenly ) {
            // Request whichever is less, the size of the buffer, or the amount which will bring us to the desired
            // stream total.
            numParticles = std::min( numParticles, static_cast<std::size_t>( m_totalCount - m_curIndex - 1 ) );
            bool notEos = m_delegate->get_particles( particleBuffer, numParticles );
            m_curIndex += numParticles;

            return notEos && m_curIndex + 1 < m_totalCount;
        } else {
            std::size_t desiredParticleCount =
                std::min( numParticles, static_cast<std::size_t>( m_totalCount - m_curIndex - 1 ) );
            std::size_t totalParticleCount = 0;
            std::size_t particleSize = m_delegate->particle_size(); // HACK: This should be stored somewhere else.

            bool notEos;

            do {
                std::size_t requestedCount = desiredParticleCount - totalParticleCount;
                notEos = m_delegate->get_particles( particleBuffer, requestedCount );

                for( char *readBuffer = particleBuffer, *readEnd = particleBuffer + requestedCount * particleSize;
                     readBuffer != readEnd; readBuffer += particleSize ) {
                    if( is_next_in_set() ) {
                        if( readBuffer != particleBuffer )
                            memcpy( particleBuffer, readBuffer, particleSize );
                        particleBuffer += particleSize;
                        ++totalParticleCount;
                    }
                }
            } while( notEos && totalParticleCount < desiredParticleCount );

            m_curIndex += totalParticleCount;
            numParticles = totalParticleCount;

            return notEos && m_curIndex + 1 < m_totalCount;
        }
    }
};

// This helper function will handle
inline boost::shared_ptr<particle_istream>
apply_fractional_particle_istream( boost::shared_ptr<particle_istream> stream, double fraction,
                                   boost::int64_t limit = 2000000000, bool evenlyDistribute = true ) {
    fraction = frantic::math::clamp( fraction, 0.0, 1.0 );

    boost::int64_t streamCount = stream->particle_count();
    if( fraction == 0.0 || streamCount == 0 )
        return boost::make_shared<empty_particle_istream>( stream->get_channel_map() );

    if( fraction == 1.0 && streamCount >= 0 && limit >= streamCount )
        return stream;

    return boost::make_shared<fractional_particle_istream>( stream, fraction, limit, evenlyDistribute );
}

} // namespace streams
} // namespace particles
} // namespace frantic
