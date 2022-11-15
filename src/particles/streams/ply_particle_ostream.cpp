// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/ply_particle_ostream.hpp>

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>

#include <boost/move/move.hpp>

namespace {

/**
 *  Create a new channel map for the file, so that we propagate only known
 * channels, and with a known type (float32[3]).
 */
frantic::channels::channel_map create_filtered_disk_channel_map( const frantic::channels::channel_map& channelMap ) {
    frantic::channels::channel_propagation_policy cpp( true );
    cpp.add_channel( _T( "Position" ) );
    cpp.add_channel( _T( "Normal" ) );
    cpp.add_channel( _T( "Color" ) );
    cpp.add_channel( _T( "TextureCoord" ) );

    frantic::channels::channel_map result;

    for( std::size_t i = 0, ie = channelMap.channel_count(); i < ie; ++i ) {
        const frantic::channels::channel& ch = channelMap[i];
        if( cpp.is_channel_included( ch.name() ) ) {
            result.define_channel( ch.name(), 3, frantic::channels::data_type_float32 );
        }
    }

    result.end_channel_definition();

    return result;
}

} // anonymous namespace

namespace frantic {
namespace particles {
namespace streams {

ply_particle_ostream::ply_particle_ostream( const frantic::tstring& filename,
                                            const frantic::channels::channel_map& channelMap,
                                            const frantic::channels::channel_map& diskChannelMap )
    : m_isClosed( false )
    , m_filename( filename ) {
    m_particles.set_channel_map( create_filtered_disk_channel_map( diskChannelMap ) );

    set_channel_map( channelMap );
}

ply_particle_ostream::~ply_particle_ostream() {}

const channel_map& ply_particle_ostream::get_channel_map() const { return m_channelMap; }

void ply_particle_ostream::set_channel_map( const channel_map& channelMap ) {
    frantic::channels::channel_map_adaptor channelMapAdaptor( m_particles.get_channel_map(), channelMap );

    m_channelMap = channelMap;
    m_channelMapAdaptor = channelMapAdaptor;
}

std::size_t ply_particle_ostream::particle_size() const { return get_channel_map().structure_size(); }

void ply_particle_ostream::close() {
    using namespace frantic::geometry;

    if( m_isClosed ) {
        throw std::runtime_error( "ply_particle_ostream::close Error: The file \"" +
                                  frantic::strings::to_string( m_filename ) + "\" was already closed." );
    }
    m_isClosed = true;

    frantic::logging::null_progress_logger progress;

    mesh_interface_ptr meshInterface( create_particle_array_mesh_interface( boost::move( m_particles ) ) );

    write_ply_mesh_file( m_filename, meshInterface, progress );
}

void ply_particle_ostream::put_particle( const char* rawParticleData ) {
    if( m_isClosed ) {
        throw std::runtime_error( "ply_particle_ostream::put_particle Error: Tried to write particle to file \"" +
                                  frantic::strings::to_string( m_filename ) + "\" after it was closed." );
    }

    m_particles.push_back( rawParticleData, m_channelMapAdaptor );
}

} // namespace streams
} // namespace particles
} // namespace frantic
