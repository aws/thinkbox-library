// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/make_shared.hpp>
#include <boost/scoped_array.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/math/fractions.hpp>
#include <frantic/math/utils.hpp>
#include <frantic/particles/streams/empty_particle_istream.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

class fractional_by_id_particle_istream : public delegated_particle_istream {
    boost::int64_t m_curIndex, m_totalCount, m_guessCount;
    boost::uint64_t m_threshold;
    frantic::tstring m_idChannel;

    frantic::channels::channel_map_adaptor m_adaptor;
    frantic::channels::channel_map m_outPcm;
    frantic::channels::channel_const_cvt_accessor<boost::int64_t> m_idAccessor;

    boost::uint64_t mix( boost::uint64_t c ) {
        boost::uint64_t a = 0, b = 1;
        a = a - b;
        a = a - c;
        a = a ^ ( c >> 43 );
        b = b - c;
        b = b - a;
        b = b ^ ( a << 9 );
        c = c - a;
        c = c - b;
        c = c ^ ( b >> 8 );
        a = a - b;
        a = a - c;
        a = a ^ ( c >> 38 );
        b = b - c;
        b = b - a;
        b = b ^ ( a << 23 );
        c = c - a;
        c = c - b;
        c = c ^ ( b >> 5 );
        a = a - b;
        a = a - c;
        a = a ^ ( c >> 35 );
        b = b - c;
        b = b - a;
        b = b ^ ( a << 49 );
        c = c - a;
        c = c - b;
        c = c ^ ( b >> 11 );
        a = a - b;
        a = a - c;
        a = a ^ ( c >> 12 );
        b = b - c;
        b = b - a;
        b = b ^ ( a << 18 );
        c = c - a;
        c = c - b;
        c = c ^ ( b >> 22 );
        return c;
    }

    bool is_next_in_set( boost::uint64_t id ) { return mix( id ) < m_threshold; }

  public:
    fractional_by_id_particle_istream( boost::shared_ptr<particle_istream> stream, double fraction,
                                       const frantic::tstring& idChannel, boost::int64_t limit = 2000000000 )
        : delegated_particle_istream( stream ) {
        m_idChannel = idChannel;
        m_curIndex = m_delegate->particle_index();
        m_totalCount = limit;

        fraction = frantic::math::clamp( fraction, 0.0, 1.0 );
        // m_threshold = (boost::uint64_t)(fraction*std::numeric_limits<boost::uint64_t>::max());

        std::pair<boost::int64_t, boost::int64_t> rational = frantic::math::get_rational_representation( fraction );

        boost::uint64_t a = static_cast<boost::uint64_t>( std::max<boost::int64_t>( 0, rational.first ) );
        boost::uint64_t b = static_cast<boost::uint64_t>( std::max<boost::int64_t>( 0, rational.second ) );
        boost::uint64_t c = std::numeric_limits<boost::uint64_t>::max();

        assert( 0 <= a && a <= b );

        // Given (a/b)*c where 0 <= a <= b, and b << c we can calculate (c/b)*a + ((c%b)*a)/b which will not overflow if
        // a*b does not overflow since (c%b)*a < b*a.
        m_threshold = ( c / b ) * a;
        m_threshold += ( a * ( c % b ) ) / b;

        boost::uint64_t gc = std::max( static_cast<boost::uint64_t>( 0 ),
                                       static_cast<boost::uint64_t>( m_delegate->particle_count_guess() ) );

        m_guessCount = static_cast<boost::int64_t>( ( gc / b ) * a );
        m_guessCount += static_cast<boost::int64_t>( ( a * ( gc % b ) ) / b );

        set_channel_map( m_delegate->get_channel_map() );
    }

    virtual ~fractional_by_id_particle_istream() {}

    boost::int64_t particle_count() const { return -1; }
    boost::int64_t particle_count_left() const { return -1; }
    boost::int64_t particle_index() const { return m_curIndex; }

    virtual boost::int64_t particle_count_guess() const { return m_guessCount; }

    std::size_t particle_size() const { return m_outPcm.structure_size(); }

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        const frantic::channels::channel_map& delegateNativePcm = m_delegate->get_native_channel_map();

        frantic::channels::channel_map newPcm = m_outPcm = pcm;

        // Append in the ID channel of the native map if it fits constraints (integer, arity one)
        // and there isn't already one.
        if( newPcm.has_channel( m_idChannel ) ) {
            frantic::channels::channel_general_accessor gca = pcm.get_general_accessor( m_idChannel );
            if( gca.arity() != 1 )
                throw std::runtime_error(
                    "frantic::channels::streams::fractional_by_id_particle_istream() - The requested ID channel (\"" +
                    frantic::strings::to_string( m_idChannel ) + "\") has arity greater than 1 (it has arity " +
                    boost::lexical_cast<std::string>( gca.arity() ) + ")." );
            if( !frantic::channels::is_channel_data_type_int( gca.data_type() ) )
                throw std::runtime_error(
                    "frantic::channels::streams::fractional_by_id_particle_istream() - The requested ID channel (\"" +
                    frantic::strings::to_string( m_idChannel ) + "\") is not of an integer type (it is of type " +
                    frantic::strings::to_string( frantic::channels::channel_data_type_str( gca.data_type() ) ) + ")." );
        } else {
            if( !delegateNativePcm.has_channel( m_idChannel ) )
                throw std::runtime_error(
                    "frantic::channels::streams::fractional_by_id_particle_istream() - The requested ID channel (\"" +
                    frantic::strings::to_string( m_idChannel ) + "\") does not exist in the input particle stream." );
            frantic::channels::channel_general_accessor gca = delegateNativePcm.get_general_accessor( m_idChannel );
            if( gca.arity() != 1 )
                throw std::runtime_error(
                    "frantic::channels::streams::fractional_by_id_particle_istream() - The requested ID channel (\"" +
                    frantic::strings::to_string( m_idChannel ) + "\") has arity greater than 1 (it has arity " +
                    boost::lexical_cast<std::string>( gca.arity() ) + ")." );
            if( !frantic::channels::is_channel_data_type_int( gca.data_type() ) )
                throw std::runtime_error(
                    "frantic::channels::streams::fractional_by_id_particle_istream() - The requested ID channel (\"" +
                    frantic::strings::to_string( m_idChannel ) + "\") is not of an integer type (it is of type " +
                    frantic::strings::to_string( frantic::channels::channel_data_type_str( gca.data_type() ) ) + ")." );
            newPcm.append_channel( m_idChannel, gca.arity(), gca.data_type() );

            /*
            // This wont work.  What if the ID is an unsigned type?
            // Or if it is a signed type, why can't you have negative IDs?
            // Also, without knowing the type, how do I set it to -1 with an accessor that has to only go 1 way?

            // set the default particle value to -1 if
            std::vector<char> p(newPcm.structure_size());
            memset(&p[0], 0, p.size());
            m_idAccessor = newPcm.get_const_cvt_accessor<boost::uint64_t>(m_idChannel);
            m_idAccessor(p) = -1;
            m_delegate->set_default_particle(&p[0]);
            */
        }

        m_idAccessor = newPcm.get_const_cvt_accessor<boost::int64_t>( m_idChannel );
        m_adaptor.set( m_outPcm, newPcm );
        m_delegate->set_channel_map( newPcm );
    }

    void set_default_particle( char* particle ) {
        if( m_adaptor.is_identity() ) {
            m_delegate->set_default_particle( particle );
        } else {
            frantic::channels::channel_map_adaptor tempAdaptor( m_delegate->get_channel_map(), m_outPcm );
            char* temp = (char*)alloca( tempAdaptor.dest_size() );
            m_delegate->get_channel_map().construct_structure( temp );
            tempAdaptor.copy_structure( temp, particle );
            m_delegate->set_default_particle( temp );
        }
    }

    bool get_particle( char* buffer ) {

        if( ++m_curIndex >= m_totalCount )
            return false;

        if( m_adaptor.is_identity() ) {
            do {
                if( !m_delegate->get_particle( buffer ) )
                    return false;
            } while( !is_next_in_set( static_cast<boost::uint64_t>( m_idAccessor( buffer ) ) ) );
        } else {
            char* tempBuffer = (char*)alloca( m_adaptor.source_size() );
            do {
                if( !m_delegate->get_particle( tempBuffer ) )
                    return false;
            } while( !is_next_in_set( static_cast<boost::uint64_t>( m_idAccessor( tempBuffer ) ) ) );
            m_adaptor.copy_structure( buffer, tempBuffer );
        }

        return true;
    }

    bool get_particles( char* particleBuffer, std::size_t& numParticles ) {
        std::size_t desiredParticleCount =
            std::min( numParticles, static_cast<std::size_t>( m_totalCount - m_curIndex - 1 ) );
        std::size_t totalParticleCount = 0;

        bool notEos;

        if( m_adaptor.is_identity() ) {
            do {
                std::size_t requestCount = desiredParticleCount - totalParticleCount;

                notEos = m_delegate->get_particles( particleBuffer, requestCount );

                for( char *srcBuffer = particleBuffer,
                          *srcBufferEnd = particleBuffer + requestCount * m_outPcm.structure_size();
                     srcBuffer != srcBufferEnd; srcBuffer += m_outPcm.structure_size() ) {
                    if( is_next_in_set( static_cast<boost::uint64_t>( m_idAccessor( srcBuffer ) ) ) ) {
                        if( particleBuffer != srcBuffer )
                            m_outPcm.copy_structure( particleBuffer, srcBuffer );
                        particleBuffer += m_outPcm.structure_size();
                        ++totalParticleCount;
                    }
                }

                // Keep asking for more particles until we at least half-fill the output buffer.
            } while( notEos && 2 * totalParticleCount < desiredParticleCount );
        } else {
            boost::scoped_array<char> tempBuffer( new char[desiredParticleCount * m_adaptor.source_size()] );
            do {
                std::size_t requestCount = desiredParticleCount - totalParticleCount;

                notEos = m_delegate->get_particles( tempBuffer.get(), requestCount );

                for( char *srcBuffer = tempBuffer.get(),
                          *srcBufferEnd = tempBuffer.get() + requestCount * m_adaptor.source_size();
                     srcBuffer != srcBufferEnd; srcBuffer += m_adaptor.source_size() ) {
                    if( is_next_in_set( static_cast<boost::uint64_t>( m_idAccessor( srcBuffer ) ) ) ) {
                        m_adaptor.copy_structure( particleBuffer, srcBuffer );
                        particleBuffer += m_adaptor.dest_size();
                        ++totalParticleCount;
                    }
                }
            } while( notEos && 2 * totalParticleCount < desiredParticleCount );
        }

        m_curIndex += totalParticleCount;
        numParticles = totalParticleCount;

        return notEos && m_curIndex + 1 < m_totalCount;
    }
};

inline boost::shared_ptr<particle_istream>
apply_fractional_by_id_particle_istream( boost::shared_ptr<particle_istream> stream, double fraction,
                                         const frantic::tstring& idChannel, boost::int64_t limit = 2000000000 ) {
    fraction = frantic::math::clamp( fraction, 0.0, 1.0 );

    if( fraction == 0.0 || stream->particle_count() == 0 )
        return boost::make_shared<empty_particle_istream>( stream->get_channel_map() );

    // We don't need to remove any particles if the limit is not set (ie. < 0 or the max value) or if we know the
    // particle count and it is less than the limit.
    if( fraction == 1.0 && ( ( limit < 0 || limit == std::numeric_limits<boost::int64_t>::max() ) ||
                             ( stream->particle_count() >= 0 && limit >= stream->particle_count() ) ) )
        return stream;

    return boost::make_shared<fractional_by_id_particle_istream>( stream, fraction, idChannel, limit );
}

} // namespace streams
} // namespace particles
} // namespace frantic
