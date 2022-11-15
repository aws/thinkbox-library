// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/trimesh3.hpp>

namespace frantic {
namespace geometry {
float get_mesh_max_velocity_magnitude( const trimesh3& mesh ) {

    if( !mesh.has_vertex_channel( _T("Velocity") ) )
        return 0.f;

    float maxMagnitudeSquared = 0, tempMagnitudeSquared = 0;

    const_trimesh3_vertex_channel_accessor<vector3f> velAcc =
        mesh.get_vertex_channel_accessor<vector3f>( _T("Velocity") );

    for( unsigned int vertexCount = 0; vertexCount < mesh.vertex_count(); vertexCount++ ) {
        tempMagnitudeSquared = velAcc[vertexCount].get_magnitude_squared();
        if( tempMagnitudeSquared > maxMagnitudeSquared )
            maxMagnitudeSquared = tempMagnitudeSquared;
    }

    return sqrt( maxMagnitudeSquared );
}
} // namespace geometry
} // namespace frantic
