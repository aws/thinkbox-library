// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/channel_map_lerp.hpp>
#include <frantic/graphics/dart_thrower.hpp>
#include <frantic/particles/particle_classes.hpp>
#include <frantic/particles/particle_cursor.hpp>
#include <frantic/particles/particle_kdtree.hpp>
#include <frantic/particles/proxy_particle_cursor.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/geometry/particle_collision_detector.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>

#include <frantic/graphics/dart_thrower.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_grid_tree.hpp>

TEST( Particle_SPH, Evaluation ) {
    using namespace frantic::particles;
    using namespace frantic::graphics;
    using namespace frantic::volumetrics;
    using namespace frantic::volumetrics::levelset;

    channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.define_channel<vector3f>( _T("Velocity") );
    pcm.define_channel<vector3f>( _T("Acceleration") );
    pcm.define_channel<float>( _T("Density") );
    pcm.define_channel<vector3f>( _T("GradDensity") );
    pcm.define_channel<float>( _T("InelasticCollisionWeightedMass") );
    pcm.define_channel<vector3f>( _T("InelasticCollisionCOMPos") );
    pcm.define_channel<vector3f>( _T("InelasticCollisionCOMVel") );
    pcm.define_channel<vector3f>( _T("InelasticCollisionMomentumTransfer") );
    pcm.define_channel<vector3f>( _T("AccGravityTerm") );
    pcm.define_channel<vector3f>( _T("AccPressureTerm") );
    pcm.define_channel<vector3f>( _T("GhostGradDensityTerm") );
    pcm.define_channel<vector3f>( _T("VelInelasticPrtTerm") );
    pcm.define_channel<vector3f>( _T("VelInelasticOccTerm") );
    pcm.define_channel<vector3f>( _T("VelEulerIntTerm") );
    pcm.end_channel_definition();

    boundbox3f bounds( vector3f(), vector3f( 20 ) );

    particle_grid_tree pgt( pcm, frantic::volumetrics::voxel_coord_system() );
}