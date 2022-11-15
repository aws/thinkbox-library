// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/mesh_measurement.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/surface_particle_istream.hpp>

#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/range/irange.hpp>

namespace frantic {
namespace geometry {

namespace detail {

// This is the primary building block of the volume computation forumla
float signed_origin_tetrahedron_volume( const frantic::graphics::vector3f& a, const frantic::graphics::vector3f& b,
                                        const frantic::graphics::vector3f& c ) {
    return vector3f::dot( a, vector3f::cross( b, c ) ) / 6.0f;
}

} // namespace detail

float tetrahedron_volume( const frantic::graphics::vector3f& a, const frantic::graphics::vector3f& b,
                          const frantic::graphics::vector3f& c, const frantic::graphics::vector3f& d ) {
    return std::abs( detail::signed_origin_tetrahedron_volume( a - d, b - d, c - d ) );
}

float get_mesh_volume( const mesh_interface* mesh, frantic::logging::progress_logger& logger ) {
    const size_t numFaces = mesh->get_num_faces();
    boost::integer_range<size_t> range( 0, numFaces );
    return get_mesh_component_volume( mesh, range.begin(), range.end(), logger, numFaces );
}

float get_mesh_volume( const mesh_interface* mesh ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_volume( mesh, nullLogger );
}

float get_mesh_volume( const trimesh3& mesh, frantic::logging::progress_logger& logger ) {
    std::unique_ptr<trimesh3_interface> iMesh = trimesh3_interface::create_instance( &mesh );
    return get_mesh_volume( iMesh.get(), logger );
}

float get_mesh_volume( const trimesh3& mesh ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_volume( mesh, nullLogger );
}

float get_mesh_volume( const polymesh3_ptr mesh, frantic::logging::progress_logger& logger ) {
    std::unique_ptr<polymesh3_interface> iMesh = polymesh3_interface::create_instance( mesh );
    return get_mesh_volume( iMesh.get(), logger );
}

float get_mesh_volume( const polymesh3_ptr mesh ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_volume( mesh, nullLogger );
}

float get_mesh_surface_area( const mesh_interface* mesh, frantic::logging::progress_logger& logger ) {
    const size_t numFaces = mesh->get_num_faces();
    boost::integer_range<size_t> range( 0, numFaces );
    return get_mesh_component_surface_area( mesh, range.begin(), range.end(), logger, numFaces );
}

float get_mesh_surface_area( const mesh_interface* mesh ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_surface_area( mesh, nullLogger );
}

float get_mesh_surface_area( const trimesh3& mesh, frantic::logging::progress_logger& logger ) {
    std::unique_ptr<trimesh3_interface> iMesh = trimesh3_interface::create_instance( &mesh );
    return get_mesh_surface_area( iMesh.get(), logger );
}

float get_mesh_surface_area( const trimesh3& mesh ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_surface_area( mesh, nullLogger );
}

float get_mesh_surface_area( const polymesh3_ptr mesh, frantic::logging::progress_logger& logger ) {
    std::unique_ptr<polymesh3_interface> iMesh = polymesh3_interface::create_instance( mesh );
    return get_mesh_surface_area( iMesh.get(), logger );
}

float get_mesh_surface_area( const polymesh3_ptr mesh ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_surface_area( mesh, nullLogger );
}

double hausdorff_distance_two_sided( const trimesh3& meshA, const trimesh3& meshB, const bool sampleVertices,
                                     const bool sampleFaces, const size_t numSamples,
                                     boost::shared_ptr<frantic::particles::particle_array> outParticlesA,
                                     boost::shared_ptr<frantic::particles::particle_array> outParticlesB ) {

    return std::max(
        hausdorff_distance_one_sided( meshA, meshB, sampleVertices, sampleFaces, numSamples, outParticlesB ),
        hausdorff_distance_one_sided( meshB, meshA, sampleVertices, sampleFaces, numSamples, outParticlesA ) );
}

double hausdorff_distance_one_sided( const trimesh3& meshA, const trimesh3& meshB, const bool sampleVertices,
                                     const bool sampleFaces, const size_t numSamples,
                                     boost::shared_ptr<frantic::particles::particle_array> outParticles ) {
    using namespace frantic::particles;
    using namespace frantic::particles::streams;
    using namespace frantic::channels;

    double distance = 0.0;

    // if one of the meshes is empty, return 0 for the distance.
    if( meshA.is_empty() || meshB.is_empty() )
        return distance;

    // set up the kdtree for finding nearest point.
    trimesh3_kdtree aKdtree( const_cast<trimesh3&>( meshA ) );
    nearest_point_search_result nearestPoint;
    size_t outIndex = 0;

    const float maxDistance = FLT_MAX;

    // Only used if particle array is passed in.
    channel_accessor<frantic::graphics::vector3f> position;
    channel_accessor<double> meshDistance;
    particle_array_iterator particleBuffer;

    if( outParticles != NULL ) {
        // set up the particle array to store all the positions sampled on meshB and the corresponding closest distances
        // to meshA.
        channel_map map;
        map.define_channel<double>( _T( "Distance" ) );
        map.define_channel( _T( "Position" ), 3, data_type_float32 );
        map.end_channel_definition();
        outParticles->set_channel_map( map );
        outParticles->resize( numSamples * 2 );

        // Get the necessary channels and buffer (particle_array).
        position = map.get_accessor<frantic::graphics::vector3f>( _T( "Position" ) );
        meshDistance = map.get_accessor<double>( _T( "Distance" ) );
        particleBuffer = outParticles->begin();
    }

    // take the first numSamples vertices and put them in the particle_array.
    if( sampleVertices ) {
        // set up a vector of all indices for vertices, shuffle the numbers.
        std::vector<size_t> indices;
        indices.reserve( meshB.vertex_count() );

        // Randomize the indices taken for sampling.
        for( size_t i = 0; i < meshB.vertex_count(); ++i ) {
            indices.push_back( i );
        }
        boost::random_shuffle( indices );

        for( size_t i = 0; i < numSamples && i < meshB.vertex_count(); ++i ) {
            if( aKdtree.find_nearest_point( meshB.get_vertex( indices[i] ), maxDistance, nearestPoint ) ) {
                if( outParticles != NULL ) {
                    position( particleBuffer[outIndex] ) = meshB.get_vertex( indices[i] );
                    meshDistance( particleBuffer[outIndex++] ) = nearestPoint.distance;
                }

                if( nearestPoint.distance > distance ) {
                    distance = nearestPoint.distance;
                }
            }
        }
    }

    // take numSamples points from the surface of the mesh and put them in the particle_array
    if( sampleFaces ) {
        // Create a stream of particles from the faces of the mesh.
        particle_istream_ptr surfaceStream = create_surface_particle_istream_using_count( meshB, numSamples );

        // Get the necessary channels and buffer (surface_particle_istream).
        channel_map faceChannelMap = surfaceStream->get_channel_map();
        std::vector<char> faceBuffer( faceChannelMap.structure_size() );
        channel_accessor<vector3f> facePosition = faceChannelMap.get_accessor<vector3f>( _T( "Position" ) );

        while( surfaceStream->get_particle( &faceBuffer[0] ) ) {
            vector3f faceSample = facePosition.get( faceBuffer );
            if( aKdtree.find_nearest_point( faceSample, maxDistance, nearestPoint ) ) {
                if( outParticles != NULL ) {
                    position( particleBuffer[outIndex] ) = faceSample;
                    meshDistance( particleBuffer[outIndex++] ) = nearestPoint.distance;
                }

                if( nearestPoint.distance > distance ) {
                    distance = nearestPoint.distance;
                }
            }
        }
    }

    // set the array to the actual size used.
    if( outParticles != NULL ) {
        outParticles->resize( outIndex );
    }

    return distance;
}

} // namespace geometry
} // namespace frantic
