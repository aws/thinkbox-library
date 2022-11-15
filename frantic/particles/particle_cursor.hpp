// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <boost/cstdint.hpp>

namespace frantic {
namespace particles {

const boost::uint32_t PRT_FLG_MODIFIED = 0x80000000;
const boost::uint32_t PRT_FLG_PROXY = 0x40000000;
const boost::uint32_t PRT_ID_MASK = 0x3FFFFFFF;

const boost::uint32_t MAX_PARTICLE_ID = 0xFFFFFFFF - PRT_FLG_PROXY - PRT_FLG_MODIFIED;

// eg. masked_uint32<0x3fffffff>
template <boost::uint32_t MaskValue>
class masked_uint32 {
    boost::uint32_t m_value;

  public:
    masked_uint32() {}
    masked_uint32( const boost::uint32_t& value ) { m_value = value; }
    bool operator<( const masked_uint32& rhs ) const { return ( m_value & MaskValue ) < ( rhs.m_value & MaskValue ); }
    bool operator>( const masked_uint32& rhs ) const { return ( m_value & MaskValue ) > ( rhs.m_value & MaskValue ); }
    bool operator==( const masked_uint32& rhs ) const { return ( m_value & MaskValue ) == ( rhs.m_value & MaskValue ); }
};

template <class T, class SearchValue>
static bool binary_array_search( SearchValue requestedValue, T begin, T end, std::size_t offset,
                                 std::size_t particleSize, T& outLocation ) {
    int headParticle = 0;
    int tailParticle = (int)( ( end - begin ) / (int)particleSize );

    while( headParticle <= tailParticle ) {
        int lookupParticle = ( headParticle + tailParticle ) / 2; // the head particle id is always '0'

        T current = begin + lookupParticle * particleSize;

        // extract the value to test against
        SearchValue currentValue = *reinterpret_cast<const SearchValue*>( current + offset );

        if( currentValue > requestedValue ) { // jump backwards
            tailParticle = lookupParticle - 1;
        } else if( currentValue < requestedValue ) { // jump forwards
            headParticle = lookupParticle + 1;
        } else {
            // we have found the value, move the actual cursor to the found value
            outLocation = current;
            return true;
        }
    }
    return false;
}

class const_particle_cursor;

//** Particle Cursor
template <class CharPtr>
class template_particle_cursor {
    // std::map<std::string,particle_channel> m_channels; // map between channel names and channel markers

    friend class const_particle_cursor;

  protected:
    CharPtr m_begin;
    CharPtr m_end;
    CharPtr m_current;
    std::size_t m_particleSize;

  public:
    template_particle_cursor( CharPtr begin, CharPtr end, std::size_t particleSize ) {
        m_begin = begin;
        m_end = end;
        // the current value of the cursor initially points to an invalid particle in front of the actual particles.
        // This is because the first call to next_particle() moves the cursor to the first particle in the array.
        // It is also the same behaviour as the particle_istream.
        m_current = m_begin - particleSize;

        m_particleSize = particleSize;
    }

    std::size_t size_of_particle() { return m_particleSize; }

    // Resest the current value of the cursor to point to an invalid particle in front of the actual particles.
    void reset() { m_current = m_begin - m_particleSize; }

    bool valid_particle() const { return m_begin <= m_current && m_current < m_end; }

    bool next_particle() {
        m_current += m_particleSize;
        return m_begin <= m_current && m_current < m_end;
    }
    bool prev_particle() {
        m_current -= m_particleSize;
        return m_begin <= m_current && m_current < m_end;
    }

    bool move_to( int particleNumber ) {
        m_current = m_begin + particleNumber * m_particleSize;
        return m_begin <= m_current && m_current < m_end;
    }

    CharPtr raw_particle_buffer() const { return m_current; }
};

class particle_cursor : public template_particle_cursor<char*> {
  public:
    particle_cursor( char* begin, char* end, std::size_t particleSize )
        : template_particle_cursor<char*>( begin, end, particleSize ) {}
};

class const_particle_cursor : public template_particle_cursor<const char*> {
  public:
    const_particle_cursor( const char* begin, const char* end, std::size_t particleSize )
        : template_particle_cursor<const char*>( begin, end, particleSize ) {}
    // Need to be able to convert a particle_cursor into a const_particle_cursor, but don't want to be able to do the
    // reverse conversion.
    const_particle_cursor( const particle_cursor& rhs )
        : template_particle_cursor<const char*>( rhs.m_begin, rhs.m_end, rhs.m_particleSize ) {
        m_current = rhs.m_current;
    }
};

} // namespace particles
} // namespace frantic
