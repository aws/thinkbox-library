// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file
 * @author Stephen Kiazyk
 * @brief Methods for computing measurements on meshes
 */
#pragma once

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>
#include <frantic/geometry/polygon_utils.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/trimesh3_interface.hpp>

#include <frantic/graphics/vector3f.hpp>

#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace geometry {

namespace detail {

/**
 * An optimized version that assumes one of the points is at the origin, and returns a 'signed' volume.
 * Used to compute more complex mesh volumes.
 */
float signed_origin_tetrahedron_volume( const frantic::graphics::vector3f& a, const frantic::graphics::vector3f& b,
                                        const frantic::graphics::vector3f& c );

} // namespace detail

/**
 * Compute the volume of a raw tetrahedron.
 *
 * @param a first vertex of the tetrahedron
 * @param b second vertex of the tetrahedron
 * @param c third vertex of the tetrahedron
 * @param d fourth vertex of the tetrahedron
 * @return the total area of the tetrahedron
 */
float tetrahedron_volume( const frantic::graphics::vector3f& a, const frantic::graphics::vector3f& b,
                          const frantic::graphics::vector3f& c, const frantic::graphics::vector3f& d );

/**
 * Get the volume of a single mesh component
 * @remark the set of faces must specific a watertight component, or set of watertight components
 *
 * @tparam FaceIterator a ForwardIterator to the face indices
 *
 * @param mesh the mesh to evaluate
 * @param facesBegin start of the face range
 * @param facesEnd end of the face range
 * @param logger progress logger
 * @param sizeHint the size of the range (0 if you don't know)
 * @return the total area of the component
 */
template <class FaceIterator>
float get_mesh_component_volume( const mesh_interface* mesh, FaceIterator facesBegin, FaceIterator facesEnd,
                                 frantic::logging::progress_logger& logger, size_t sizeHint ) {
    using namespace frantic::graphics;

    const size_t updateInterval = 1024;

    float result = 0.0f;

    std::vector<vector3f> vertices( 3 );
    std::vector<vector3> triangles( 1 );

    size_t currentIndex = 0;

    for( FaceIterator it = facesBegin; it != facesEnd; ++it ) {

        const size_t faceSize = mesh->get_num_face_verts( *it );
        vertices.resize( faceSize );
        vertices.assign( face_vertex_iterator( mesh, *it, face_vertex_iterator::begin ),
                         face_vertex_iterator( mesh, *it, face_vertex_iterator::end ) );

        if( vertices.size() > 3 ) {
            triangles.resize( vertices.size() - 2 );
            triangulate_polygon( &vertices[0], vertices.size(), &triangles[0] );
        } else {
            triangles[0] = vector3( 0, 1, 2 );
        }

        for( size_t triangleId = 0; triangleId < triangles.size(); ++triangleId ) {
            const vector3& triangle = triangles[triangleId];
            const float faceVolume = detail::signed_origin_tetrahedron_volume(
                vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]] );
            result += faceVolume;
        }

        if( currentIndex % updateInterval == 0 ) {
            if( sizeHint > 0 ) {
                logger.update_progress( currentIndex, sizeHint );
            } else {
                logger.check_for_abort();
            }
        }

        ++currentIndex;
    }

    return std::abs( result );
}

template <class FaceIterator>
float get_mesh_component_volume( const trimesh3& mesh, FaceIterator facesBegin, FaceIterator facesEnd,
                                 frantic::logging::progress_logger& logger, size_t sizeHint ) {
    std::unique_ptr<trimesh3_interface> iMesh = trimesh3_interface::create_instance( &mesh );
    return get_mesh_component_volume( iMesh.get(), facesBegin, facesEnd, logger, sizeHint );
}

template <class FaceIterator>
float get_mesh_component_volume( const trimesh3& mesh, FaceIterator facesBegin, FaceIterator facesEnd ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_component_volume( mesh, facesBegin, facesEnd, nullLogger, 0 );
}

template <class FaceIterator>
float get_mesh_component_volume( const polymesh3_ptr& mesh, FaceIterator facesBegin, FaceIterator facesEnd,
                                 frantic::logging::progress_logger& logger, size_t sizeHint ) {
    std::unique_ptr<polymesh3_interface> iMesh = polymesh3_interface::create_instance( mesh );
    return get_mesh_component_volume( iMesh.get(), facesBegin, facesEnd, logger, sizeHint );
}

template <class FaceIterator>
float get_mesh_component_volume( const polymesh3_ptr& mesh, FaceIterator facesBegin, FaceIterator facesEnd ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_component_volume( mesh, facesBegin, facesEnd, nullLogger, 0 );
}

/**
 * Get the volume of a mesh component
 * @remark the mesh must be watertight
 *
 * @param mesh the mesh to evaluate
 * @return the total area of the mesh
 */
float get_mesh_volume( const mesh_interface* mesh, frantic::logging::progress_logger& logger );
float get_mesh_volume( const mesh_interface* mesh );

float get_mesh_volume( const trimesh3& mesh, frantic::logging::progress_logger& logger );
float get_mesh_volume( const trimesh3& mesh );

float get_mesh_volume( const polymesh3_ptr mesh, frantic::logging::progress_logger& logger );
float get_mesh_volume( const polymesh3_ptr mesh );

/**
 * Get the surface area of all the faces in a mesh component
 *
 * @tparam FaceIterator a forward iterator to the face indices
 *
 * @param mesh
 * @param facesBegin iterator to the first face in the component
 * @param facesEnd iterator to the end face in the component
 * @param logger progress logger
 * @param sizeHint the size of the range (0 if you don't know)
 * @return the sum area of all the faces
 */
template <class FaceIterator>
float get_mesh_component_surface_area( const mesh_interface* mesh, FaceIterator facesBegin, FaceIterator facesEnd,
                                       frantic::logging::progress_logger& logger, size_t sizeHint ) {
    const size_t updateInterval = 1024;

    size_t currentIndex = 0;

    double result = 0;

    for( FaceIterator it = facesBegin; it != facesEnd; ++it ) {
        result += polygon3_area( face_vertex_iterator( mesh, *it, face_vertex_iterator::begin ),
                                 face_vertex_iterator( mesh, *it, face_vertex_iterator::end ) );

        if( currentIndex % updateInterval == 0 ) {
            if( sizeHint > 0 ) {
                logger.update_progress( currentIndex, sizeHint );
            } else {
                logger.check_for_abort();
            }
        }

        ++currentIndex;
    }

    return static_cast<float>( result );
}

template <class FaceIterator>
float get_mesh_component_surface_area( const mesh_interface* mesh, FaceIterator facesBegin, FaceIterator facesEnd ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_component_surface_area( mesh, facesBegin, facesEnd, nullLogger, 0 );
}

template <class FaceIterator>
float get_mesh_component_surface_area( const trimesh3& mesh, FaceIterator facesBegin, FaceIterator facesEnd,
                                       frantic::logging::progress_logger& logger, size_t sizeHint ) {
    std::unique_ptr<trimesh3_interface> iMesh = trimesh3_interface::create_instance( &mesh );
    return get_mesh_component_surface_area( iMesh.get(), facesBegin, facesEnd, logger, sizeHint );
}

template <class FaceIterator>
float get_mesh_component_surface_area( const trimesh3& mesh, FaceIterator facesBegin, FaceIterator facesEnd ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_component_surface_area( mesh, facesBegin, facesEnd, nullLogger, 0 );
}

template <class FaceIterator>
float get_mesh_component_surface_area( const polymesh3_ptr mesh, FaceIterator facesBegin, FaceIterator facesEnd,
                                       frantic::logging::progress_logger& logger, size_t sizeHint ) {
    std::unique_ptr<polymesh3_interface> iMesh = polymesh3_interface::create_instance( mesh );
    return get_mesh_component_surface_area( iMesh.get(), facesBegin, facesEnd, logger, sizeHint );
}

template <class FaceIterator>
float get_mesh_component_surface_area( const polymesh3_ptr mesh, FaceIterator facesBegin, FaceIterator facesEnd ) {
    frantic::logging::null_progress_logger nullLogger;
    return get_mesh_component_surface_area( mesh, facesBegin, facesEnd, nullLogger, 0 );
}

/**
 * Get the surface area of all the faces in a mesh
 *
 * @param mesh
 * @return the sum area of all the faces
 */
float get_mesh_surface_area( const mesh_interface* mesh, frantic::logging::progress_logger& logger );
float get_mesh_surface_area( const mesh_interface* mesh );

float get_mesh_surface_area( const trimesh3& mesh, frantic::logging::progress_logger& logger );
float get_mesh_surface_area( const trimesh3& mesh );

float get_mesh_surface_area( const polymesh3_ptr mesh, frantic::logging::progress_logger& logger );
float get_mesh_surface_area( const polymesh3_ptr mesh );

/**
 *	Find the symmetric hausdorff distance between two meshes. (Maximum of the maximum minimal distances)
 *
 *	@param meshA The mesh that closest points will be found on.
 *	@param meshB The mesh to check each vertex of for a closest point on meshA.
 *	@param sampleVertices Do we want to sample points from vertices?
 *	@param sampleFaces Do we want to sample points from faces?
 *	@param numSamples The number of samples to take from each of vertices and faces.
 *	@param[out] outParticlesA The particle array to store the info in (optional).
 *	@param[out] outParticlesB The particle array to store the info in (optional).
 *
 *	@return The Hausdorff distance between the two meshes (The maximum of the maximum minimal distances)
 */
double hausdorff_distance_two_sided( const trimesh3& meshA, const trimesh3& meshB, const bool sampleVertices,
                                     const bool sampleFaces, const size_t numSamples,
                                     boost::shared_ptr<frantic::particles::particle_array> outParticlesA =
                                         boost::shared_ptr<frantic::particles::particle_array>(),
                                     boost::shared_ptr<frantic::particles::particle_array> outParticlesB =
                                         boost::shared_ptr<frantic::particles::particle_array>() );

/**
 *	Calculates the closest point on MeshA to each sample point on MeshB, stores the points from meshB and the
 *distance into a particle array.
 *
 *	@param meshA The mesh that closest points will be found on.
 *	@param meshB The mesh to sample points from.
 *	@param sampleVertices Do we want to sample points from vertices?
 *	@param sampleFaces Do we want to sample points from faces?
 *	@param numSamples The number of samples to take from each of vertices and faces.
 *   @param[out] outParticles The particle array to store the info in.
 *
 *	@return The maximum of the minimal distances between the two meshes.
 */
double hausdorff_distance_one_sided( const trimesh3& meshA, const trimesh3& meshB, const bool sampleVertices,
                                     const bool sampleFaces, const size_t numSamples,
                                     boost::shared_ptr<frantic::particles::particle_array> outParticles );

} // namespace geometry
} // namespace frantic
