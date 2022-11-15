// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/age_culled_particle_istream.hpp>

#pragma warning( push, 3 )
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#pragma warning( pop )

namespace frantic {
namespace particles {
namespace streams {

particle_istream_ptr age_culled_particle_istream::apply_to_stream( particle_istream_ptr pin, float minAge ) {
    if( minAge > 1e-5f )
        pin.reset( new age_culled_particle_istream( pin, minAge ) );
    return pin;
}

age_culled_particle_istream::age_culled_particle_istream( particle_istream_ptr pDelegate, float minAge )
    : delegated_particle_istream( pDelegate )
    , m_minAge( minAge )
    , m_particleIndex( -1 ) {
    this->set_channel_map_impl( pDelegate->get_channel_map() );
}

void age_culled_particle_istream::set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
    this->set_channel_map_impl( particleChannelMap );
}

void age_culled_particle_istream::set_default_particle( char* rawParticleBuffer ) {
    frantic::channels::channel_map_adaptor tempAdaptor( m_delegate->get_channel_map(), m_outMap );

    char* tempBuffer = reinterpret_cast<char*>( alloca( tempAdaptor.dest_size() ) );

    tempAdaptor.copy_structure( tempBuffer, rawParticleBuffer );

    m_delegate->set_default_particle( tempBuffer );
}

const frantic::channels::channel_map& age_culled_particle_istream::get_channel_map() const { return m_outMap; }

inline bool age_culled_particle_istream::cull( char* pParticle ) const {
    return !m_ageAccessor.is_valid() || m_ageAccessor.get( pParticle ) < m_minAge;
}

bool age_culled_particle_istream::get_particle( char* pBuffer ) {
    if( !m_adaptor.is_identity() ) {
        char* tempBuffer = reinterpret_cast<char*>( alloca( m_adaptor.source_size() ) );

        do {
            if( !m_delegate->get_particle( tempBuffer ) )
                return false;
        } while( this->cull( tempBuffer ) );

        m_adaptor.copy_structure( pBuffer, tempBuffer );
    } else {
        do {
            if( !m_delegate->get_particle( pBuffer ) )
                return false;
        } while( this->cull( pBuffer ) );
    }

    ++m_particleIndex;

    return true;
}

namespace {
void operator_delete( void* ptr ) { operator delete( ptr ); }

template <class T>
struct blocked_range_less {
    bool operator()( const tbb::blocked_range<T>& lhs, const tbb::blocked_range<T>& rhs ) const {
        return lhs.begin() < rhs.begin();
    }
};
} // namespace

class age_culled_particle_istream::culling_body {
  private:
    age_culled_particle_istream* m_pOwner;

    char* m_pBuffer;

    typedef age_culled_particle_istream::vector_type vector_type;

    vector_type m_resultRanges;

  public:
    culling_body( age_culled_particle_istream& owner, char* pBuffer )
        : m_pOwner( &owner )
        , m_pBuffer( pBuffer ) {
        m_resultRanges.reserve( 2 * tbb::task_scheduler_init::default_num_threads() );
    }

    culling_body( culling_body& lhs, tbb::split )
        : m_pOwner( lhs.m_pOwner )
        , m_pBuffer( lhs.m_pBuffer ) {}

    void operator()( const tbb::blocked_range<std::size_t>& range ) {
        m_pOwner->cull_particles1( m_pBuffer, range, m_resultRanges );
    }

    void join( culling_body& rhs ) {
        vector_type::const_iterator itRhs = rhs.m_resultRanges.begin(), itRhsEnd = rhs.m_resultRanges.end();

        if( !m_resultRanges.empty() && itRhs != itRhsEnd && m_resultRanges.back().end() == itRhs->begin() ) {
            // Prefer to merge blocked_range objects instead of just merging the lists.
            m_resultRanges.back() =
                tbb::blocked_range<std::size_t>( m_resultRanges.back().begin(),
                                                 ( itRhs++ )->end() ); // NOTE: Post-increment used here.
        }

        std::size_t oldSize = m_resultRanges.size();

        // Append the elements onto the end of the list, then do an inplace_merge to get them resorted.
        m_resultRanges.insert( m_resultRanges.end(), itRhs, itRhsEnd );

        std::inplace_merge( m_resultRanges.begin(), m_resultRanges.begin() + oldSize, m_resultRanges.end(),
                            blocked_range_less<std::size_t>() );
    }

    const vector_type& get_result() const { return m_resultRanges; }
};

bool age_culled_particle_istream::get_particles( char* pBuffer, std::size_t& numParticles ) {
    bool eos;

    if( !m_adaptor.is_identity() ) {
        // We need to extract particles to a separate buffer first.
        boost::shared_ptr<void> pTempBuffer( operator new( numParticles* m_adaptor.source_size() ), &operator_delete );

        eos = !m_delegate->get_particles( reinterpret_cast<char*>( pTempBuffer.get() ), numParticles );

        culling_body theBody( *this, reinterpret_cast<char*>( pTempBuffer.get() ) );

        tbb::parallel_reduce( tbb::blocked_range<std::size_t>( 0, numParticles, 1000 ), theBody,
                              tbb::auto_partitioner() );

        this->cull_particles2( pBuffer, reinterpret_cast<char*>( pTempBuffer.get() ), theBody.get_result(),
                               numParticles );
    } else {
        eos = !m_delegate->get_particles( pBuffer, numParticles );

        culling_body theBody( *this, pBuffer );

        tbb::parallel_reduce( tbb::blocked_range<std::size_t>( 0, numParticles, 1000 ), theBody,
                              tbb::auto_partitioner() );

        this->cull_particles2( pBuffer, pBuffer, theBody.get_result(), numParticles );
    }

    m_particleIndex += numParticles;

    return !eos;
}

void age_culled_particle_istream::cull_particles1( char* pBuffer, const tbb::blocked_range<std::size_t>& range,
                                                   vector_type& outValidRanges ) const {
    pBuffer = pBuffer + m_adaptor.source_size() * range.begin();

    char* pOut = pBuffer;

    std::size_t numLeft = 0;

    for( std::size_t i = range.begin(); i < range.end(); ++i, pBuffer += m_adaptor.source_size() ) {
        if( !this->cull( pBuffer ) ) {
            ++numLeft;

            if( pOut != pBuffer )
                memcpy( pOut, pBuffer, m_adaptor.source_size() );

            pOut += m_adaptor.source_size();
        }
    }

    if( numLeft > 0 ) {
        if( !outValidRanges.empty() && outValidRanges.back().end() == range.begin() ) {
            outValidRanges.back() = tbb::blocked_range<std::size_t>(
                outValidRanges.back().begin(),
                range.begin() + numLeft ); // Prefer to extend the range than to make a new one.
        } else {
            outValidRanges.push_back( tbb::blocked_range<std::size_t>( range.begin(), range.begin() + numLeft ) );
        }
    }
}

void age_culled_particle_istream::cull_particles2( char* pOutBuffer, const char* pBuffer,
                                                   const vector_type& validRanges, std::size_t& outCount ) const {
    outCount = 0;

    if( validRanges.empty() )
        return;

    if( pBuffer != pOutBuffer ) {
        vector_type::const_iterator it = validRanges.begin(), itEnd = validRanges.end();

        for( ; it != itEnd; ++it ) {
            const char* pSrc = pBuffer + it->begin() * m_adaptor.source_size();

            for( std::size_t i = it->begin(), iEnd = it->end(); i < iEnd;
                 ++i, pSrc += m_adaptor.source_size(), pOutBuffer += m_adaptor.dest_size() )
                m_adaptor.copy_structure( pOutBuffer, pSrc );

            outCount += it->size();
        }
    } else {
        // We can optimize by copying out of order
        vector_type::const_iterator it = validRanges.begin(), itEnd = validRanges.end();
        vector_type::const_iterator itSrc = itEnd - 1;

        std::size_t dest = 0, src = itSrc->end();

        // Search for a gap to fill.
        while( it != itEnd && dest < src && dest == it->begin() )
            dest = ( it++ )->end();

        while( it != itEnd && dest < src ) {
            char* pDest = pOutBuffer + dest * m_adaptor.dest_size();
            const char* pSrc;

            std::size_t gapSize = it->begin() - dest;
            std::size_t srcSize = src - itSrc->begin();

            if( gapSize >= srcSize ) {
                pSrc = pBuffer + itSrc->begin() * m_adaptor.source_size();

                // Our destination region can hold this entire source range so copy all of them.
                for( std::size_t destEnd = dest + srcSize; dest < destEnd;
                     ++dest, pDest += m_adaptor.dest_size(), pSrc += m_adaptor.source_size() )
                    m_adaptor.copy_structure( pDest, pSrc );

                // Move the source block iterator unless it has already reached the dest block in which case we are
                // done.
                src = ( itSrc != it ) ? ( --itSrc )->end() : dest;
            } else {
                src -= gapSize;
                pSrc = pBuffer + src * m_adaptor.source_size();

                // Our destination can only copy part of the source range, so copy what we can. We copy from the range
                // end so that we only move the items that fill in gaps without producing new gaps.
                for( std::size_t destEnd = it->begin(); dest < destEnd;
                     ++dest, pDest += m_adaptor.dest_size(), pSrc += m_adaptor.source_size() )
                    m_adaptor.copy_structure( pDest, pSrc );

                dest = ( it++ )->end();
            }

            // Search for a gap to fill if the current gap was filled.
            while( it != itEnd && dest < src && dest == it->begin() )
                dest = ( it++ )->end();
        }

        outCount = src;
    }
}

void age_culled_particle_istream::set_channel_map_impl( const frantic::channels::channel_map& channelMap ) {
    m_outMap = channelMap;

    frantic::channels::channel_map delegateMap = channelMap;
    if( !delegateMap.has_channel( _T("Age") ) && m_delegate->get_native_channel_map().has_channel( _T("Age") ) )
        delegateMap.append_channel<float>( _T("Age") );

    if( delegateMap != m_delegate->get_channel_map() )
        m_delegate->set_channel_map( delegateMap );

    m_adaptor.set( m_outMap, delegateMap );

    m_ageAccessor.reset();

    if( delegateMap.has_channel( _T("Age") ) )
        m_ageAccessor = delegateMap.get_cvt_accessor<float>( _T("Age") );
}

} // namespace streams
} // namespace particles
} // namespace frantic
