// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/geometry/volume_collection.hpp>
#include <frantic/particles/streams/culling_particle_istream.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

class volume_culling_policy : public culling_policy_base<volume_culling_policy> {
    boost::shared_ptr<frantic::geometry::volume_collection> m_volume;
    frantic::channels::channel_accessor<vector3f> m_posAccessor;
    bool m_cullOutside;

  public:
    typedef boost::tuples::tuple<boost::shared_ptr<frantic::geometry::volume_collection>, bool> args_type;

    volume_culling_policy( const args_type& args, const frantic::channels::channel_map& pcm )
        : m_volume( args.get<0>() )
        , m_cullOutside( args.get<1>() ) {
        set_channel_map( pcm );
    }

    // Need a splitting constructor to support TBB
    volume_culling_policy( const volume_culling_policy& lhs, tbb::split )
        : m_volume( lhs.m_volume )
        , m_cullOutside( lhs.m_cullOutside )
        , m_posAccessor( lhs.m_posAccessor ) {}

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_posAccessor = pcm.get_accessor<vector3f>( _T("Position") );
    }

    bool cull( char* particle ) const {
        return m_volume->is_point_in_volume( m_posAccessor.get( particle ) ) != m_cullOutside;
    }
};

typedef culling_particle_istream<volume_culling_policy> volume_culled_particle_istream;

} // namespace streams
} // namespace particles
} // namespace frantic
