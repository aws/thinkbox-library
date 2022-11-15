// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/particles/streams/rle_levelset_particle_istream.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

using namespace frantic::graphics;
using namespace frantic::volumetrics::levelset;

namespace frantic {
namespace particles {
namespace streams {

rle_levelset_particle_istream::rle_levelset_particle_istream(
    const frantic::channels::channel_map& pcm,
    boost::shared_ptr<frantic::volumetrics::levelset::rle_level_set> pLevelset,
    frantic::volumetrics::voxel_sampler_interface_ptr pParticleGenerator, float innerDistance, float outerDistance,
    bool compensateDensity )
    : m_pLevelset( pLevelset )
    , m_pParticleGenerator( pParticleGenerator )
    , m_outerDistance( outerDistance )
    , m_innerDistance( innerDistance )
    , m_iter( pLevelset->get_rle_index_spec(), false )
    , m_endIter( pLevelset->get_rle_index_spec(), true )
    , m_compensateDensity( compensateDensity ) {
    m_particleIndex = -1;
    m_particleProgressIndex = -1;
    m_particleProgressCount = (boost::int64_t)pLevelset->size();

    m_compensationFactor = pLevelset->get_voxel_coord_system().voxel_length();
    m_compensationFactor *= m_compensationFactor * m_compensationFactor;

    set_channel_map( pcm );

    std::vector<frantic::tstring> levelsetChannels;
    pLevelset->get_channel_names( levelsetChannels );

    m_nativeMap.define_channel<vector3f>( _T("Position") );
    m_nativeMap.define_channel<vector3f>( _T("Normal") );
    m_nativeMap.define_channel<vector3f>( _T("SignedDistanceGradient") );
    m_nativeMap.define_channel<float>( _T("SignedDistance") );
    m_nativeMap.define_channel<float>( _T("Density") );
    for( std::vector<frantic::tstring>::iterator it = levelsetChannels.begin(), itEnd = levelsetChannels.end();
         it != itEnd; ++it ) {
        rle_channel_general_accessor acc = pLevelset->get_channel_general_accessor( *it );
        if( !m_nativeMap.has_channel( *it ) )
            m_nativeMap.define_channel( *it, acc.arity(), acc.data_type() );
    }
    m_nativeMap.end_channel_definition();
}

rle_levelset_particle_istream::~rle_levelset_particle_istream() {}

void rle_levelset_particle_istream::copy_channel_data( const frantic::graphics::vector3f& localCoord,
                                                       float localDistance, boost::int32_t trilerpIndices[],
                                                       float trilerpWeights[], char* pParticle ) {
    // Copy the default values.
    m_outMap.copy_structure( pParticle, m_defaultParticle.get() );

    // Set the world-space position
    const frantic::volumetrics::voxel_coord_system& vcs = m_pLevelset->get_voxel_coord_system();
    m_posAccessor.get( pParticle ) = vcs.get_world_coord( localCoord );

    // Set the signed distance channel ... TODO: Should this be positive instead?
    if( m_distAccessor.is_valid() )
        m_distAccessor.set( pParticle, localDistance );

    if( m_signedDistanceGradientAccessor.is_valid() ) {
        vector3f signedDistGrad;
        trilerp_staggered_centered_signed_distance_gradient( *m_pLevelset, localCoord, signedDistGrad );

        m_signedDistanceGradientAccessor.set( pParticle, signedDistGrad );

        if( m_normalAccessor.is_valid() )
            m_normalAccessor.set( pParticle, vector3f::normalize( signedDistGrad ) );
    } else if( m_normalAccessor.is_valid() ) {
        vector3f signedDistGrad;
        trilerp_staggered_centered_signed_distance_gradient( *m_pLevelset, localCoord, signedDistGrad );

        m_normalAccessor.set( pParticle, vector3f::normalize( signedDistGrad ) );
    }

    // TODO: currently this will barf on channels with more than 128 bits. The correct version
    //        would allocate storage for the largest channel that needs to get copied.
    float tempSpace[4];
    for( std::vector<channel_connection>::iterator it = m_channels.begin(), itEnd = m_channels.end(); it != itEnd;
         ++it ) {
        char* particleChannel = it->particleAccessor.get_channel_data_pointer( pParticle );
        char* levelsetChannel = reinterpret_cast<char*>( tempSpace );

        std::size_t arity = it->particleAccessor.arity();
        it->levelsetAccessor.get_trilinear_interpolated( trilerpWeights, trilerpIndices, levelsetChannel );
        it->conversionFunction( particleChannel, levelsetChannel, arity );
    }
}

bool rle_levelset_particle_istream::advance_iterator() {
    ++m_particleProgressIndex;
    if( ++m_iter == m_endIter )
        return false;

    const frantic::volumetrics::voxel_coord_system& vcs = m_pLevelset->get_voxel_coord_system();

    static const float SQRT3 = 1.7320508075688772935274463415059f;
    float minDistance = m_innerDistance - vcs.voxel_length() * SQRT3;
    float maxDistance = m_outerDistance + vcs.voxel_length() * SQRT3;

    // Iterate through the defined voxels until the cell defined by (x,y,z)->(x+1,y+1,z+1)
    // is overlapping the seeding region.
    int voxelIndex = m_iter.get_data_index();
    float distance = ( *m_pLevelset )[voxelIndex];
    while( distance > maxDistance || distance < minDistance ) {
        ++m_particleProgressIndex;
        if( ++m_iter == m_endIter )
            return false;
        voxelIndex = m_iter.get_data_index();
        distance = ( *m_pLevelset )[voxelIndex];
    }

    m_pParticleGenerator->update_for_voxel( m_iter.get_coord() );

    return true;
}

void rle_levelset_particle_istream::set_channel_map( const frantic::channels::channel_map& pcm ) {
    { // Swap in a new default particle.
        boost::scoped_array<char> newDefault( new char[pcm.structure_size()] );
        pcm.construct_structure( newDefault.get() );

        if( m_defaultParticle ) {
            frantic::channels::channel_map_adaptor tempAdaptor( pcm, m_outMap );
            tempAdaptor.copy_structure( newDefault.get(), m_defaultParticle.get() );
        }

        m_defaultParticle.swap( newDefault );
    }

    m_outMap = pcm;
    m_posAccessor = m_outMap.get_accessor<vector3f>( _T("Position") );

    m_normalAccessor.reset();
    m_signedDistanceGradientAccessor.reset();
    m_distAccessor.reset();
    m_densityAccessor.reset();
    m_channels.clear();

    if( m_outMap.has_channel( _T("Normal") ) )
        m_normalAccessor = m_outMap.get_cvt_accessor<vector3f>( _T("Normal") );
    if( m_outMap.has_channel( _T("SignedDistanceGradient") ) )
        m_signedDistanceGradientAccessor = m_outMap.get_cvt_accessor<vector3f>( _T("SignedDistanceGradient") );
    if( m_outMap.has_channel( _T("SignedDistance") ) )
        m_distAccessor = m_outMap.get_cvt_accessor<float>( _T("SignedDistance") );
    if( m_outMap.has_channel( _T("Density") ) && m_compensateDensity )
        m_densityAccessor = m_outMap.get_cvt_accessor<float>( _T("Density") );

    for( std::size_t i = 0, iEnd = m_outMap.channel_count(); i != iEnd; ++i ) {
        const frantic::channels::channel& ch = m_outMap[i];
        // TODO: Should I try to handle this?
        // if( ch.name() == "Position" || ch.name() == "Normal" || ch.name() == "SignedDistance" )
        //	continue;
        if( m_pLevelset->has_channel( ch.name() ) ) {
            channel_connection con;
            con.levelsetAccessor = m_pLevelset->get_channel_general_accessor( ch.name() );
            con.particleAccessor = m_outMap.get_general_accessor( ch.name() );

            if( con.levelsetAccessor.arity() != con.particleAccessor.arity() )
                throw std::runtime_error( "rle_levelset_particle_istream() - Could not match the volume's channel \"" +
                                          frantic::strings::to_string( ch.name() ) + "\" of arity " +
                                          boost::lexical_cast<std::string>( con.levelsetAccessor.arity() ) +
                                          " to the requested particle channel's arity of " +
                                          boost::lexical_cast<std::string>( con.particleAccessor.arity() ) );
            if( con.levelsetAccessor.primitive_size() > 16 )
                throw std::runtime_error(
                    "rle_levelset_particle_istream() - Channel \"" + frantic::strings::to_string( ch.name() ) + "\" " +
                    frantic::strings::to_string( ch.type_str() ) + " is too large for the current implementation" );

            con.conversionFunction = frantic::channels::get_channel_type_convertor_function(
                con.levelsetAccessor.data_type(), con.particleAccessor.data_type(), ch.name() );
            m_channels.push_back( con );
        }
    }
}

void rle_levelset_particle_istream::set_default_particle( char* rawParticleBuffer ) {
    m_outMap.copy_structure( m_defaultParticle.get(), rawParticleBuffer );
}

bool rle_levelset_particle_istream::get_particle( char* rawParticleBuffer ) {
    float distance;
    float cellWeights[8];
    float compensationFactor;
    vector3f localOffset;

    do {
        if( !m_pParticleGenerator->get_next_position( localOffset, compensationFactor ) ) {
            do {
                if( !advance_iterator() )
                    return false;
                ++m_particleProgressIndex;
            } while( !m_pParticleGenerator->get_next_position( localOffset, compensationFactor ) );

            m_pLevelset->get_rle_index_spec().fill_2x2x2_data_index_box( m_iter.get_coord(), m_cachedIndices );
        }

        get_trilerp_weights( &localOffset.x, cellWeights );

        distance = 0;
        for( std::size_t i = 0; i < 8; ++i )
            distance += cellWeights[i] * m_pLevelset->get_using_data_index( m_cachedIndices[i] );
    } while( distance > m_outerDistance || distance < m_innerDistance );

    vector3f localCoord = m_iter.get_coord() + vector3f( 0.5f ) + localOffset;
    copy_channel_data( localCoord, distance, m_cachedIndices, cellWeights, rawParticleBuffer );

    // Compensate the density if we are seeding more than one particle per voxel
    // BUG: Why is this squared?!?
    if( m_densityAccessor.is_valid() )
        m_densityAccessor.set( rawParticleBuffer, m_compensationFactor * compensationFactor );

    ++m_particleIndex;
    return true;
}

bool rle_levelset_particle_istream::get_particles( char* rawParticleBuffer, std::size_t& inoutNumParticles ) {
    std::size_t particleSize = particle_size();
    for( std::size_t i = 0; i < inoutNumParticles; ++i ) {
        if( !get_particle( rawParticleBuffer ) ) {
            inoutNumParticles = i;
            return false;
        }
        rawParticleBuffer += particleSize;
    }
    return true;
}

} // namespace streams
} // namespace particles
} // namespace frantic
