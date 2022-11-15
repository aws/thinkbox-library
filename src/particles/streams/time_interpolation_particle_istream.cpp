// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/time_interpolation_particle_istream.hpp>

using namespace frantic::particles::streams;

void frantic::particles::streams::detail::make_other_particles_map( frantic::channels::channel_map& outMap,
                                                                    const frantic::channels::channel_map& map ) {
    outMap.reset();
    outMap.define_channel( _T("Position"), 3, frantic::channels::data_type_float32 );
    // handle both the signed and unsigned cases.
    bool isUint = ( map.has_channel( _T("ID") ) &&
                    frantic::channels::is_channel_data_type_unsigned( map[_T("ID")].data_type() ) );
    outMap.define_channel( _T("ID"), 1,
                           isUint ? frantic::channels::data_type_uint64 : frantic::channels::data_type_int64 );
    if( map.has_channel( _T("Velocity") ) )
        outMap.define_channel<frantic::graphics::vector3f>( _T("Velocity") );

    if( map.has_channel( _T("Normal") ) )
        outMap.define_channel<frantic::graphics::vector3f>( _T("Normal") );

    if( map.has_channel( _T("Tangent") ) )
        outMap.define_channel<frantic::graphics::vector3f>( _T("Tangent") );

    if( map.has_channel( _T("Color") ) )
        outMap.define_channel<frantic::graphics::vector3f>( _T("Color") );

    if( map.has_channel( _T("TextureCoord") ) )
        outMap.define_channel<frantic::graphics::vector3f>( _T("TextureCoord") );

    if( map.has_channel( _T("Age") ) )
        outMap.define_channel<float>( _T("Age") );

    if( map.has_channel( _T("Density") ) )
        outMap.define_channel<float>( _T("Density") );

    outMap.end_channel_definition();
}

void frantic::particles::streams::time_interpolation_particle_istream::set_channel_map_internal(
    const frantic::channels::channel_map& pcm ) {
    m_outMap = pcm;

    const frantic::channels::channel_map& nativeMap = m_delegate->get_native_channel_map();

    frantic::channels::channel_map delegateMap = m_outMap;
    if( !delegateMap.has_channel( _T("ID") ) ) {
        if( nativeMap.has_channel( _T("ID") ) ) {
            const frantic::channels::channel& ch = nativeMap[_T("ID")];

            delegateMap.append_channel( _T("ID"), ch.arity(), ch.data_type() );
            m_forceExtrapolation = false;
        } else {
            delegateMap.append_channel<boost::int64_t>( _T("ID") );
            m_forceExtrapolation = true;
        }
    } else {
        m_forceExtrapolation = false;
    }

    if( !delegateMap.has_channel( _T("Velocity") ) &&
        m_delegate->get_native_channel_map().has_channel( _T("Velocity") ) )
        delegateMap.append_channel<frantic::graphics::vector3f>( _T("Velocity") );

    if( !delegateMap.has_channel( _T("Normal") ) && m_delegate->get_native_channel_map().has_channel( _T("Normal") ) )
        delegateMap.append_channel<frantic::graphics::vector3f>( _T("Normal") );

    if( !delegateMap.has_channel( _T("Tangent") ) && m_delegate->get_native_channel_map().has_channel( _T("Tangent") ) )
        delegateMap.append_channel<frantic::graphics::vector3f>( _T("Tangent") );

    if( !delegateMap.has_channel( _T("Color") ) && m_delegate->get_native_channel_map().has_channel( _T("Color") ) )
        delegateMap.append_channel<frantic::graphics::vector3f>( _T("Color") );

    if( !delegateMap.has_channel( _T("TextureCoord") ) &&
        m_delegate->get_native_channel_map().has_channel( _T("TextureCoord") ) )
        delegateMap.append_channel<frantic::graphics::vector3f>( _T("TextureCoord") );

    if( !delegateMap.has_channel( _T("Age") ) && m_delegate->get_native_channel_map().has_channel( _T("Age") ) )
        delegateMap.append_channel<float>( _T("Age") );

    if( !delegateMap.has_channel( _T("Density") ) && m_delegate->get_native_channel_map().has_channel( _T("Density") ) )
        delegateMap.append_channel<float>( _T("Density") );

    m_outAdaptor.set( m_outMap, delegateMap );

    m_streamPosAccessor = delegateMap.get_accessor<frantic::graphics::vector3f>( _T("Position") );
    m_streamIDAccessor.reset( delegateMap.get_general_accessor( _T("ID") ) );

    m_streamVelAccessor.reset();
    if( delegateMap.has_channel( _T("Velocity") ) )
        m_streamVelAccessor = delegateMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Velocity") );

    m_streamNormAccessor.reset();
    if( delegateMap.has_channel( _T("Normal") ) )
        m_streamNormAccessor = delegateMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Normal") );

    m_streamTangentAccessor.reset();
    if( delegateMap.has_channel( _T("Tangent") ) )
        m_streamTangentAccessor = delegateMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Tangent") );

    m_streamColorAccessor.reset();
    if( delegateMap.has_channel( _T("Color") ) )
        m_streamColorAccessor = delegateMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Color") );

    m_streamTextureCoordAccessor.reset();
    if( delegateMap.has_channel( _T("TextureCoord") ) )
        m_streamTextureCoordAccessor = delegateMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("TextureCoord") );

    m_streamAgeAccessor.reset();
    if( delegateMap.has_channel( _T("Age") ) )
        m_streamAgeAccessor = delegateMap.get_cvt_accessor<float>( _T("Age") );

    m_streamDensityAccessor.reset();
    if( delegateMap.has_channel( _T("Density") ) )
        m_streamDensityAccessor = delegateMap.get_cvt_accessor<float>( _T("Density") );

    m_delegate->set_channel_map( delegateMap );
}

void frantic::particles::streams::time_interpolation_particle_istream::set_other_accessors() {
    const frantic::channels::channel_map& map = m_otherParticles.get_channel_map();
    m_otherPosAccessor = map.get_accessor<frantic::graphics::vector3f>( _T("Position") );
    m_otherIDAccessor = map.get_general_accessor( _T("ID") );
    if( map.has_channel( _T("Velocity") ) )
        m_otherVelAccessor = map.get_cvt_accessor<frantic::graphics::vector3f>( _T("Velocity") );
    if( map.has_channel( _T("Normal") ) )
        m_otherNormAccessor = map.get_cvt_accessor<frantic::graphics::vector3f>( _T("Normal") );
    if( map.has_channel( _T("Tangent") ) )
        m_otherTangentAccessor = map.get_cvt_accessor<frantic::graphics::vector3f>( _T("Tangent") );
    if( map.has_channel( _T("Color") ) )
        m_otherColorAccessor = map.get_cvt_accessor<frantic::graphics::vector3f>( _T("Color") );
    if( map.has_channel( _T("TextureCoord") ) )
        m_otherTextureCoordAccessor = map.get_cvt_accessor<frantic::graphics::vector3f>( _T("TextureCoord") );
    if( map.has_channel( _T("Age") ) )
        m_otherAgeAccessor = map.get_cvt_accessor<float>( _T("Age") );
    if( map.has_channel( _T("Density") ) )
        m_otherDensityAccessor = map.get_cvt_accessor<float>( _T("Density") );
}

const char*
frantic::particles::streams::time_interpolation_particle_istream::find_particle_match( boost::int64_t id,
                                                                                       bool& outSuccess ) const {
    const bool isUint =
        frantic::channels::is_channel_data_type_unsigned( m_otherParticles.get_channel_map()[_T( "ID" )].data_type() );
    if( isUint ) {
        particle_compare_id<boost::uint64_t> comp( m_otherParticles.get_channel_map() );
        frantic::particles::particle_array::const_iterator result =
            std::lower_bound( m_otherParticles.begin(), m_otherParticles.end(), id, comp );
        if( result != m_otherParticles.end() ) {
            outSuccess = true;
            return *result;
        }
    } else {
        particle_compare_id<boost::int64_t> comp( m_otherParticles.get_channel_map() );
        frantic::particles::particle_array::const_iterator result =
            std::lower_bound( m_otherParticles.begin(), m_otherParticles.end(), id, comp );
        if( result != m_otherParticles.end() ) {
            outSuccess = true;
            return *result;
        }
    }

    outSuccess = false;
    return NULL;
}

void frantic::particles::streams::time_interpolation_particle_istream::process_particle( char* particle ) const {
    using namespace frantic::graphics;

    boost::int64_t id = m_streamIDAccessor.get( particle );

    const char* particleMatch = 0;
    bool foundMatch = false;
    if( m_otherParticles.size() > 0 && m_otherIDAccessor.is_valid() )
        particleMatch = find_particle_match( id, foundMatch );

    // change the particle
    if( !m_forceExtrapolation && foundMatch && m_otherParticles.size() > 0 &&
        m_otherIDAccessor( particleMatch ) == id ) {

        vector3f lhs = m_streamPosAccessor.get( particle );
        vector3f rhs = m_otherPosAccessor.get( particleMatch );

        if( m_streamVelAccessor.is_valid() && m_otherVelAccessor.is_valid() ) {
            // Do cubic Hermite interpolation using the provided positions and derivatives.
            vector3f lhsVel = m_streamVelAccessor.get( particle );
            vector3f rhsVel = m_otherVelAccessor.get( particleMatch );
            vector3f lhsCont = ( lhs + ( m_timeStepSeconds / 3.f ) * lhsVel );
            vector3f rhsCont = ( rhs - ( m_timeStepSeconds / 3.f ) * rhsVel );

            vector3f outPos = frantic::math::get_bezier_curve_position( lhs, lhsCont, rhsCont, rhs, m_interpVal );
            vector3f outVel = frantic::math::get_bezier_curve_derivative( lhs, lhsCont, rhsCont, rhs, m_interpVal ) /
                              m_timeStepSeconds;

            m_streamPosAccessor.get( particle ) = outPos;
            m_streamVelAccessor.set( particle, outVel );
        } else {
            // There was no velocity information so we do straight lerping.
            m_streamPosAccessor.get( particle ) = frantic::math::lerp( lhs, rhs, m_interpVal );
        }

        if( m_streamNormAccessor.is_valid() && m_otherNormAccessor.is_valid() ) {
            quat4f lhsQuat( 0, m_streamNormAccessor.get( particle ) );
            quat4f rhsQuat( 0, m_otherNormAccessor.get( particleMatch ) );
            lhsQuat.normalize();
            rhsQuat.normalize();
            m_streamNormAccessor.set( particle,
                                      frantic::graphics::quat4f::slerp( lhsQuat, rhsQuat, m_interpVal ).vector_part() );
        }

        if( m_streamTangentAccessor.is_valid() && m_otherTangentAccessor.is_valid() ) {
            quat4f lhsQuat( 0, m_streamTangentAccessor.get( particle ) );
            quat4f rhsQuat( 0, m_otherTangentAccessor.get( particleMatch ) );
            lhsQuat.normalize();
            rhsQuat.normalize();
            m_streamTangentAccessor.set(
                particle, frantic::graphics::quat4f::slerp( lhsQuat, rhsQuat, m_interpVal ).vector_part() );
        }

        if( m_streamColorAccessor.is_valid() && m_otherColorAccessor.is_valid() )
            m_streamColorAccessor.set( particle,
                                       frantic::math::lerp( m_streamColorAccessor.get( particle ),
                                                            m_otherColorAccessor.get( particleMatch ), m_interpVal ) );

        if( m_streamTextureCoordAccessor.is_valid() && m_otherTextureCoordAccessor.is_valid() )
            m_streamTextureCoordAccessor.set(
                particle, frantic::math::lerp( m_streamTextureCoordAccessor.get( particle ),
                                               m_otherTextureCoordAccessor.get( particleMatch ), m_interpVal ) );

        if( m_streamAgeAccessor.is_valid() ) {
            if( m_otherAgeAccessor.is_valid() ) {
                m_streamAgeAccessor.set( particle,
                                         frantic::math::lerp( m_streamAgeAccessor.get( particle ),
                                                              m_otherAgeAccessor.get( particleMatch ), m_interpVal ) );
            } else {
                m_streamAgeAccessor.set( particle, m_streamAgeAccessor.get( particle ) + m_offsetSeconds );
            }
        }

        if( m_streamDensityAccessor.is_valid() && m_otherDensityAccessor.is_valid() )
            m_streamDensityAccessor.set( particle, frantic::math::lerp( m_streamDensityAccessor.get( particle ),
                                                                        m_otherDensityAccessor.get( particleMatch ),
                                                                        m_interpVal ) );

    } else {
        // There was no matching particle, so we just use the velocity to offset.
        if( m_streamVelAccessor.is_valid() )
            m_streamPosAccessor.get( particle ) += m_offsetSeconds * m_streamVelAccessor.get( particle );
        if( m_streamAgeAccessor.is_valid() )
            m_streamAgeAccessor.set( particle, m_streamAgeAccessor.get( particle ) + m_offsetSeconds );
    }
}

void frantic::particles::streams::time_interpolation_particle_istream::process_particles_in_place(
    char* pBuffer, const tbb::blocked_range<std::size_t>& range ) const {
    pBuffer += range.begin() * m_outMap.structure_size();

    for( std::size_t i = range.begin(); i < range.end(); ++i, pBuffer += m_outMap.structure_size() )
        this->process_particle( pBuffer );
}

void frantic::particles::streams::time_interpolation_particle_istream::process_particles(
    char* pDestBuffer, char* pSrcBuffer, const tbb::blocked_range<std::size_t>& range ) const {
    pSrcBuffer += range.begin() * m_outAdaptor.source_size();
    pDestBuffer += range.begin() * m_outAdaptor.dest_size();

    for( std::size_t i = range.begin(); i < range.end();
         ++i, pSrcBuffer += m_outAdaptor.source_size(), pDestBuffer += m_outAdaptor.dest_size() ) {
        this->process_particle( pSrcBuffer );

        m_outAdaptor.copy_structure( pDestBuffer, pSrcBuffer );
    }
}

void frantic::particles::streams::time_interpolation_particle_istream::sort_other_particles() {
    bool isUint =
        frantic::channels::is_channel_data_type_unsigned( m_otherParticles.get_channel_map()[_T("ID")].data_type() );
    if( isUint ) {
        frantic::sort::parallel_sort( m_otherParticles.begin(), m_otherParticles.end(),
                                      particle_compare_id<boost::uint64_t>( m_otherParticles.get_channel_map() ) );
    } else {
        frantic::sort::parallel_sort( m_otherParticles.begin(), m_otherParticles.end(),
                                      particle_compare_id<boost::int64_t>( m_otherParticles.get_channel_map() ) );
    }
}

frantic::particles::streams::time_interpolation_particle_istream::time_interpolation_particle_istream(
    boost::shared_ptr<particle_istream> mainPin, boost::shared_ptr<particle_istream> otherTimePin,
    float timeStepSeconds, float interpVal )
    : m_interpVal( interpVal )
    , m_timeStepSeconds( timeStepSeconds )
    , m_forceExtrapolation( false )
    , delegated_particle_istream( mainPin ) {

    m_offsetSeconds = ( interpVal * timeStepSeconds );

    set_channel_map_internal( m_delegate->get_channel_map() );

    frantic::channels::channel_map otherMap;
    detail::make_other_particles_map( otherMap, m_delegate->get_native_channel_map() );

    m_otherParticles.reset( otherMap );

    m_otherParticles.insert_particles( otherTimePin );

    set_other_accessors();

    sort_other_particles();
}

// called by apply_interpolation_particle_istream_with_lifespan_culling
frantic::particles::streams::time_interpolation_particle_istream::time_interpolation_particle_istream(
    boost::shared_ptr<particle_istream> mainPin, frantic::particles::particle_array interpParticles,
    float timeStepSeconds, float interpVal )
    : m_interpVal( interpVal )
    , m_timeStepSeconds( timeStepSeconds )
    , m_forceExtrapolation( false )
    , delegated_particle_istream( mainPin ) {

    m_offsetSeconds = ( interpVal * timeStepSeconds );

    set_channel_map_internal( m_delegate->get_channel_map() );

    m_otherParticles.swap( interpParticles );

    set_other_accessors();

    sort_other_particles();
}

void frantic::particles::streams::time_interpolation_particle_istream::set_channel_map(
    const frantic::channels::channel_map& pcm ) {
    set_channel_map_internal( pcm );
}

const frantic::channels::channel_map&
frantic::particles::streams::time_interpolation_particle_istream::get_channel_map() const {
    return m_outMap;
}

void frantic::particles::streams::time_interpolation_particle_istream::set_default_particle( char* defaultParticle ) {
    if( m_outAdaptor.is_identity() ) {
        m_delegate->set_default_particle( defaultParticle );
    } else {
        frantic::channels::channel_map_adaptor tempAdaptor( m_delegate->get_channel_map(), m_outMap );

        char* buffer = (char*)alloca( m_delegate->get_channel_map().structure_size() );

        m_delegate->get_channel_map().construct_structure( buffer );

        tempAdaptor.copy_structure( buffer, defaultParticle );

        m_delegate->set_default_particle( buffer );
    }
}

bool frantic::particles::streams::time_interpolation_particle_istream::get_particle( char* outParticle ) {
    char* buffer = m_outAdaptor.is_identity() ? outParticle : (char*)alloca( m_outAdaptor.source_size() );

    if( !m_delegate->get_particle( buffer ) )
        return false;

    process_particle( buffer );

    if( !m_outAdaptor.is_identity() )
        m_outAdaptor.copy_structure( outParticle, buffer );

    return true;
}

bool frantic::particles::streams::time_interpolation_particle_istream::get_particles( char* buffer,
                                                                                      std::size_t& numParticles ) {
    bool eos;

    if( m_outAdaptor.is_identity() ) {
        eos = !m_delegate->get_particles( buffer, numParticles );

        tbb::parallel_for(
            tbb::blocked_range<std::size_t>( 0, numParticles, 100 ),
            boost::bind( &time_interpolation_particle_istream::process_particles_in_place, this, buffer, _1 ),
            tbb::auto_partitioner() );
    } else {
        // We need to extract particles to a separate buffer first.
        boost::shared_ptr<void> pTempBuffer( operator new( numParticles* m_outAdaptor.source_size() ),
                                             &operator_delete );

        eos = !m_delegate->get_particles( reinterpret_cast<char*>( pTempBuffer.get() ), numParticles );

        tbb::parallel_for( tbb::blocked_range<std::size_t>( 0, numParticles, 100 ),
                           boost::bind( &time_interpolation_particle_istream::process_particles, this, buffer,
                                        reinterpret_cast<char*>( pTempBuffer.get() ), _1 ),
                           tbb::auto_partitioner() );
    }
    return !eos;
}

namespace frantic {
namespace particles {
namespace streams {

boost::shared_ptr<particle_istream>
apply_interpolation_particle_istream_with_lifespan_culling( boost::shared_ptr<particle_istream> pin,
                                                            boost::shared_ptr<particle_istream> otherPin,
                                                            double timeStepSeconds, double interpParam ) {

    // If interpolation will do basically nothing (due to being almost exactly at the endpoint) we just return the
    // 'main' stream.
    if( std::abs( interpParam ) < 1e-5 )
        return pin;

    // TODO: We are screwed if we don't have an ID channel.

    // If we don't have age and lifespan, then we can't do anything fancy to smoothly create or delete particles. We can
    // only interpolate from the nearest sample
    if( !pin->get_native_channel_map().has_channel( _T("Age") ) ||
        !pin->get_native_channel_map().has_channel( _T("LifeSpan") ) ) {
        pin = boost::make_shared<time_interpolation_particle_istream>(
            pin, otherPin, static_cast<float>( timeStepSeconds ), static_cast<float>( interpParam ) );
    } else {
        // We want to always interpolate backward from furthest in time sample, so we might swap the 'main' and 'other'
        // streams. This will make sure timeStepSeconds is always negative.
        if( timeStepSeconds > 0 ) {
            timeStepSeconds = -timeStepSeconds;
            interpParam = 1.0 - interpParam;
            pin.swap( otherPin );
        }

        // The timeOffset is calculated relative to the earliest stream in time (ie. otherPin)
        float otherTimeStep = static_cast<float>( -timeStepSeconds );
        float otherTimeOffset = static_cast<float>( otherTimeStep * ( 1.0 - interpParam ) );

        // Collect interpolation data into this std::vector<>. We need position, velocity, and ID information.
        frantic::particles::particle_array interpParticles;

        boost::shared_ptr<frantic::particles::particle_array> aliveParticles;

        const frantic::channels::channel_map& map = otherPin->get_native_channel_map();

        otherPin->set_channel_map( map );

        frantic::channels::channel_map otherMap;
        detail::make_other_particles_map( otherMap, map );

        interpParticles.reset( otherMap );

        frantic::channels::channel_map_adaptor adaptor( otherMap, map );

        frantic::channels::channel_accessor<frantic::graphics::vector3f> posAccessor =
            map.get_accessor<frantic::graphics::vector3f>( _T("Position") );
        frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> velAccessor(
            frantic::graphics::vector3f( 0.f ) );
        frantic::channels::channel_cvt_accessor<float> ageAccessor( 0.f ), lifespanAccessor( 0.f );

        if( map.has_channel( _T("Velocity") ) )
            velAccessor = map.get_cvt_accessor<frantic::graphics::vector3f>( _T("Velocity") );

        if( map.has_channel( _T("Age") ) )
            ageAccessor = map.get_cvt_accessor<float>( _T("Age") );

        if( map.has_channel( _T("LifeSpan") ) )
            lifespanAccessor = map.get_cvt_accessor<float>( _T("LifeSpan") );

        // Assume that no particles died between the two samples, and reserve space for them.
        boost::int64_t expectedCount = pin->particle_count();
        if( expectedCount > 0 )
            interpParticles.reserve( expectedCount );

        // We need to store the entire particle if it is missing from 'pin' but should still be alive.
        aliveParticles = boost::make_shared<frantic::particles::particle_array>( map );

        bool done;

        char* tempParticle = (char*)alloca( otherMap.structure_size() );

        boost::scoped_array<char> buffer( new char[10000 * map.structure_size()] );

        do {
            std::size_t numParticles = 10000;

            done = !otherPin->get_particles( buffer.get(), numParticles );

            for( std::size_t i = 0; i < numParticles; ++i ) {
                char* particle = buffer.get() + map.structure_size() * i;

                float age = ageAccessor.get( particle );
                float lifespan = lifespanAccessor.get( particle );
                float lifeLeft = lifespan - age;

                if( lifeLeft >= otherTimeStep ) {
                    // This particle will be alive in the main sample. We need it for interpolation.
                    // Note: This might be violated if the framerate was changed compared to when the particle set was
                    // generated.
                    adaptor.copy_structure( tempParticle, particle );

                    interpParticles.push_back( tempParticle );
                } else if( lifeLeft >= otherTimeOffset ) {
                    // This particle will be dead in the main sample, but should still be alive at the requested time!
                    posAccessor.get( particle ) += velAccessor.get( particle ) * otherTimeOffset;
                    ageAccessor.set( particle, age + otherTimeOffset );

                    aliveParticles->push_back( particle );
                } else {
                    // This particle should be dead at the requested time.
                    // Do nothing.
                }
            }
        } while( !done );

        // Filter out any particles in the "main" stream that shouldn't be born yet.
        pin = age_culled_particle_istream::apply_to_stream( pin, static_cast<float>( -timeStepSeconds * interpParam ) );

        // Interpolate all particles that were alive in both streams.
        pin.reset( new time_interpolation_particle_istream( pin, interpParticles, static_cast<float>( timeStepSeconds ),
                                                            static_cast<float>( interpParam ) ) );

        // Append the particles that are dead in the "main" stream but are still alive at the interpolation time.
        std::size_t numExtra = aliveParticles->size();
        if( numExtra > 0 ) {
            boost::shared_ptr<shared_particle_container_particle_istream<frantic::particles::particle_array>> pExtra =
                boost::make_shared<shared_particle_container_particle_istream<frantic::particles::particle_array>>(
                    aliveParticles );

            pExtra->set_channel_map( pin->get_channel_map() );

            pin = boost::make_shared<concatenated_particle_istream>( pin, pExtra );
        }
    }
    return pin;
}

} // namespace streams
} // namespace particles
} // namespace frantic
