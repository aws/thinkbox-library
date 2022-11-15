// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/ply_particle_istream.hpp>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include <frantic/channels/channel_propagation_policy.hpp>

namespace frantic {
namespace particles {
namespace streams {

namespace {

frantic::channels::channel_map get_file_channel_map( const frantic::geometry::ply_reader& reader ) {
    frantic::channels::channel_propagation_policy cpp( true );
    cpp.add_channel( _T("Color") );
    cpp.add_channel( _T("TextureCoord") );
    cpp.add_channel( _T("Normal") );

    std::vector<frantic::tstring> vertexChannelNames;
    reader.get_vertex_channel_names( vertexChannelNames );

    frantic::channels::channel_map channelMap;
    channelMap.define_channel<frantic::graphics::vector3f>( _T("Position") );
    BOOST_FOREACH( const frantic::tstring& channelName, vertexChannelNames ) {
        if( cpp.is_channel_included( channelName ) ) {
            frantic::channels::data_type_t dataType;
            std::size_t arity;
            boost::tie( dataType, arity ) = reader.get_vertex_channel_type( channelName );
            channelMap.define_channel( channelName, arity, dataType );
        }
    }
    channelMap.end_channel_definition();
    return channelMap;
}

std::string get_data_type_string( const frantic::channels::channel_general_accessor& acc ) {
    return frantic::strings::to_string( frantic::channels::channel_data_type_str( acc.arity(), acc.data_type() ) );
}

std::string get_data_type_string( const frantic::geometry::polymesh3_const_vertex_accessor<void>& acc ) {
    return frantic::strings::to_string( frantic::channels::channel_data_type_str( acc.get_arity(), acc.get_type() ) );
}

/**
 *  Copy mesh vertex positions and channels to the outParticles.
 *
 *  outParticles must have its channel map set before calling this function,
 * and mesh must have a matching vertex channel for each particle channel.
 */
void copy_vertices_to_particles( const frantic::geometry::polymesh3_ptr mesh,
                                 frantic::particles::particle_array& outParticles ) {
    using frantic::graphics::vector3f;

    typedef frantic::channels::channel_general_accessor particle_accessor_t;
    typedef frantic::geometry::polymesh3_const_vertex_accessor<void> mesh_accessor_t;
    typedef std::pair<particle_accessor_t, mesh_accessor_t> accessor_pair_t;

    const std::size_t vertexCount = mesh->vertex_count();

    const frantic::channels::channel_map& channelMap = outParticles.get_channel_map();

    frantic::channels::channel_accessor<vector3f> positionAcc = channelMap.get_accessor<vector3f>( _T( "Position" ) );

    std::vector<accessor_pair_t> accessors;

    for( std::size_t i = 0, ie = channelMap.channel_count(); i < ie; ++i ) {
        const frantic::channels::channel& ch = channelMap[i];
        if( ch.name() == _T( "Position" ) ) {
            continue;
        }

        if( !mesh->has_vertex_channel( ch.name() ) ) {
            throw std::runtime_error( "copy_vertices_to_particles Internal Error: mesh is missing a "
                                      "\"" +
                                      frantic::strings::to_string( ch.name() ) + "\" channel" );
        }

        particle_accessor_t particleAcc = channelMap.get_general_accessor( ch.name() );
        mesh_accessor_t meshAcc = mesh->get_const_vertex_accessor( ch.name() );

        if( particleAcc.data_type() != meshAcc.get_type() || particleAcc.arity() != meshAcc.get_arity() ) {
            const std::string msg = "copy_vertices_to_particles Internal Error: mismatch between "
                                    "particle type (" +
                                    get_data_type_string( particleAcc ) +
                                    ") and "
                                    "mesh type (" +
                                    get_data_type_string( meshAcc ) +
                                    ") "
                                    "in channel \"" +
                                    frantic::strings::to_string( ch.name() ) + "\"";
            throw std::runtime_error( msg.c_str() );
        }

        accessors.push_back( accessor_pair_t( particleAcc, meshAcc ) );
    }

    outParticles.resize( vertexCount );

    for( std::size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {
        char* outParticle = outParticles.at( vertexIndex );

        positionAcc( outParticles.at( vertexIndex ) ) = mesh->get_vertex( vertexIndex );

        BOOST_FOREACH( accessor_pair_t& accessorPair, accessors ) {
            particle_accessor_t& particleAcc = accessorPair.first;
            mesh_accessor_t& meshAcc = accessorPair.second;

            char* dest = particleAcc.get_channel_data_pointer( outParticle );
            const char* src = meshAcc.get_vertex( vertexIndex );

            memcpy( dest, src, particleAcc.primitive_size() );
        }
    }
}

} // anonymous namespace

ply_particle_istream::ply_particle_istream( const frantic::tstring& filename )
    : m_isOpen( false )
    , m_particleIndex( -1 )
    , m_reader( filename ) {
    m_particles.reset( get_file_channel_map( m_reader ) );

    set_channel_map( m_particles.get_channel_map() );
}

ply_particle_istream::~ply_particle_istream() {}

void ply_particle_istream::close() {}

std::size_t ply_particle_istream::particle_size() const { return m_channelMap.structure_size(); }

frantic::tstring ply_particle_istream::name() const { return m_reader.get_filename(); }

boost::int64_t ply_particle_istream::particle_count() const {
    return static_cast<boost::int64_t>( m_reader.get_vertex_count() );
}

boost::int64_t ply_particle_istream::particle_index() const { return m_particleIndex; }

boost::int64_t ply_particle_istream::particle_count_left() const { return particle_count() - m_particleIndex - 1; }

boost::int64_t ply_particle_istream::particle_progress_count() const { return particle_count(); }

boost::int64_t ply_particle_istream::particle_progress_index() const { return m_particleIndex; }

void ply_particle_istream::set_channel_map( const frantic::channels::channel_map& channelMap ) {
    std::vector<char> newDefaultParticle( channelMap.structure_size() );
    if( newDefaultParticle.size() > 0 && m_defaultParticle.size() > 0 ) {
        frantic::channels::channel_map_adaptor adaptor( channelMap, m_channelMap );
        adaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticle[0] );
    }
    m_defaultParticle.swap( newDefaultParticle );

    m_channelMap = channelMap;
    m_channelMapAdaptor.set( m_channelMap, get_native_channel_map() );
}

const frantic::channels::channel_map& ply_particle_istream::get_channel_map() const { return m_channelMap; }

const frantic::channels::channel_map& ply_particle_istream::get_native_channel_map() const {
    return m_particles.get_channel_map();
}

void ply_particle_istream::set_default_particle( char* rawParticleBuffer ) {
    if( !m_defaultParticle.empty() ) {
        m_channelMap.copy_structure( &m_defaultParticle[0], rawParticleBuffer );
    }
}

bool ply_particle_istream::get_particle( char* rawParticleBuffer ) {
    assert( rawParticleBuffer );

    if( !m_isOpen ) {
        m_isOpen = true;

        frantic::geometry::polymesh3_ptr mesh = m_reader.read_polymesh3();
        if( !mesh ) {
            throw std::runtime_error( "ply_particle_istream::get_particle Error: mesh is NULL" );
        }

        copy_vertices_to_particles( mesh, m_particles );
    }

    const boost::int64_t currentIndex = m_particleIndex + 1;
    if( currentIndex < particle_count() ) {
        m_channelMapAdaptor.copy_structure( rawParticleBuffer, m_particles.at( currentIndex ), &m_defaultParticle[0] );

        ++m_particleIndex;

        return true;
    } else {
        return false;
    }
}

bool ply_particle_istream::get_particles( char* rawParticleBuffer, std::size_t& numParticles ) {
    const std::size_t particleSize = m_channelMap.structure_size();
    for( std::size_t i = 0; i < numParticles; ++i ) {
        if( !get_particle( rawParticleBuffer + i * particleSize ) ) {
            numParticles = i;
            return false;
        }
    }

    return true;
}

} // namespace streams
} // namespace particles
} // namespace frantic
