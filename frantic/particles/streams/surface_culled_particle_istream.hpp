// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/geometry/volume_collection.hpp>
#include <frantic/particles/streams/culling_particle_istream.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

class surface_culling_policy : public culling_policy_base<surface_culling_policy> {
    boost::shared_ptr<frantic::geometry::volume_collection> m_volume;
    channels::channel_accessor<vector3f> m_posAccessor;
    channels::channel_cvt_accessor<vector3f> m_normalAccessor;
    float m_cullDistance;
    bool m_surfaceCull;
    bool m_setNormal;

  public:
    // Arguments are:
    // 1. The volume collection to surface cull with.
    // 2. Whether to cull particles not "near" the surface
    // 3. Whether to assign a normal to particles "near" the surface
    // 4. The distance that is considered "near" the surface.
    typedef boost::tuples::tuple<boost::shared_ptr<frantic::geometry::volume_collection>, bool, bool, float> args_type;

    surface_culling_policy( const args_type& args, const frantic::channels::channel_map& pcm )
        : m_volume( args.get<0>() )
        , m_surfaceCull( args.get<1>() )
        , m_setNormal( args.get<2>() )
        , m_cullDistance( args.get<3>() ) {
        set_channel_map( pcm );
    }

    // Need a splitting constructor to support TBB
    surface_culling_policy( const surface_culling_policy& lhs, tbb::split )
        : m_volume( lhs.m_volume )
        , m_surfaceCull( lhs.m_surfaceCull )
        , m_setNormal( lhs.m_setNormal )
        , m_cullDistance( lhs.m_cullDistance )
        , m_posAccessor( lhs.m_posAccessor )
        , m_normalAccessor( lhs.m_normalAccessor ) {}

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_posAccessor = pcm.get_accessor<vector3f>( _T("Position") );
        if( m_setNormal && pcm.has_channel( _T("Normal") ) )
            m_normalAccessor = pcm.get_cvt_accessor<vector3f>( _T("Normal") );
        else
            m_normalAccessor.reset();
    }

    bool cull( char* particle ) const {
        vector3f normal;
        float dist = m_volume->get_distance_to_surface( m_posAccessor.get( particle ), &normal );

        if( !m_surfaceCull || dist <= m_cullDistance ) {
            if( m_normalAccessor.is_valid() ) {
                normal = frantic::math::lerp( normal, m_normalAccessor.get( particle ),
                                              frantic::math::smoothstep( dist, 0.f, m_cullDistance ) );
                normal.normalize();
                m_normalAccessor.set( particle, normal );
            }
            return false;
        }
        return true;
    }
};

typedef culling_particle_istream<surface_culling_policy> surface_culled_particle_istream;

} // namespace streams
} // namespace particles
} // namespace frantic
