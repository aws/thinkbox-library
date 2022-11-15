// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/*
Equations for smoothing taken from

[1] A. Belyaev, I. Bogaevski, and Y. Ohtake. Polyhedral Surface Smoothing
with Simultaneous Mesh Regularization. Geometric Modeling and Processing,
pages 229-237, 2000.

[2] M. Desbrun, M. Meyer, P. Schr�oder, and A. H. Barr. Implicit Fairing of
Irregular Meshes using Diffusion and Curvature Flow.
Computer Graphics (SIGGRAPH 99 Proceedings), pages 317�324, 1999.

[3] A. Belyaev, Y. Ohtake, and A. Pasko. Dynamic Mesh Optimization for Polygonized
Implicit Surfaces with Sharp Features. The Visual Computer, 2003, Volume 19,
Numbers 2-3, Pages 115-126.
*/

#pragma once

#include <vector>

#include <frantic/geometry/mesh_interface.hpp>

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/logging/progress_logger.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace frantic {
namespace geometry {
namespace relaxation {

namespace detail {

/**
 * Sums the area of all triangles in the faces vector.
 */
inline float sum_adjacent_areas( const std::vector<int>& adjFaces,
                                 const std::vector<frantic::graphics::vector3f>& vertices,
                                 const std::vector<frantic::graphics::vector3>& faces ) {
    float sum = 0;
    for( std::size_t i = 0; i < adjFaces.size(); ++i ) {
        sum += triangle_area( vertices[faces[adjFaces[i]].x], vertices[faces[adjFaces[i]].y],
                              vertices[faces[adjFaces[i]].z] );
    }

    return sum;
}

/**
 * The Z operator moves vertices towards an implicit surface define with the gradient of the surface. Equation (5) from
 * [3].
 */
template <class ImplicitSurfacePolicy>
frantic::graphics::vector3f z_force( std::vector<frantic::graphics::vector3f>& vertices,
                                     std::vector<frantic::graphics::vector3>& faces, std::vector<int>& adjFaces,
                                     int pos, /*float t, */ ImplicitSurfacePolicy& isp, float h ) {
    // Z(P, v(P, f) ) = -2 t A(P) f(P) gradf(P)
    // t = 1 / 100 * max( A |f grad f| )

    float A = sum_adjacent_areas( adjFaces, vertices, faces );
    float f = isp.get_density( vertices[pos] );
    frantic::graphics::vector3f gradient = isp.get_gradient( vertices[pos], h );

    float t = 1 / ( 100 * std::max( A, ( f * gradient ).get_magnitude() ) );
    return ( -2 * t * A * f * gradient );
}

/********************** Classes for the tbb parallel loops **********************/
template <class ImplicitSurfacePolicy>
class apply_z_force {
  private:
    std::vector<frantic::graphics::vector3f>& m_currVertices;
    std::vector<frantic::graphics::vector3f>& m_prevVertices;
    std::vector<frantic::graphics::vector3>& m_faces;
    std::vector<std::vector<int>>& m_adjFaces;
    ImplicitSurfacePolicy& m_isp;
    float h;

    apply_z_force& operator=( const apply_z_force& ) {}

  public:
    apply_z_force( std::vector<frantic::graphics::vector3f>& currVertices,
                   std::vector<frantic::graphics::vector3f>& prevVertices,
                   std::vector<frantic::graphics::vector3>& faces, std::vector<std::vector<int>>& adjFaces,
                   ImplicitSurfacePolicy& isp, float h )
        : m_currVertices( currVertices )
        , m_prevVertices( prevVertices )
        , m_faces( faces )
        , m_adjFaces( adjFaces )
        , m_isp( isp )
        , h( h ) {}

    void operator()( const tbb::blocked_range<size_t>& r ) const {
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            m_currVertices[i] = m_prevVertices[i] + z_force( m_prevVertices, m_faces, m_adjFaces[i], (int)i, m_isp, h );
        }
    }
};

class apply_n_force {
  private:
    std::vector<frantic::graphics::vector3f>& m_currVertices;
    std::vector<frantic::graphics::vector3f>& m_prevVertices;
    std::vector<frantic::graphics::vector3>& m_faces;
    std::vector<std::vector<int>>& m_adjFaces;
    std::vector<frantic::graphics::vector3f>& m_gradients;

    apply_n_force& operator=( const apply_n_force& ); // not implemented

  public:
    apply_n_force( std::vector<frantic::graphics::vector3f>& currVertices,
                   std::vector<frantic::graphics::vector3f>& prevVertices,
                   std::vector<frantic::graphics::vector3>& faces, std::vector<std::vector<int>>& adjFaces,
                   std::vector<frantic::graphics::vector3f>& gradients )
        : m_currVertices( currVertices )
        , m_prevVertices( prevVertices )
        , m_faces( faces )
        , m_adjFaces( adjFaces )
        , m_gradients( gradients ) {}

    // equation between 5 and 6 in [3]
    void operator()( const tbb::blocked_range<size_t>& r ) const;
};

template <class ImplicitSurfacePolicy>
class apply_find_gradients {
  private:
    std::vector<frantic::geometry::vector3f>& m_gradients;
    std::vector<frantic::geometry::vector3>& m_faces;
    std::vector<frantic::graphics::vector3f>& m_vertices;
    ImplicitSurfacePolicy& m_isp;
    float h;

    apply_find_gradients& operator=( const apply_find_gradients& ) {}

  public:
    apply_find_gradients( std::vector<frantic::geometry::vector3f>& gradients,
                          std::vector<frantic::geometry::vector3>& faces,
                          std::vector<frantic::graphics::vector3f>& vertices, ImplicitSurfacePolicy& isp, float h )
        : m_gradients( gradients )
        , m_faces( faces )
        , m_vertices( vertices )
        , m_isp( isp )
        , h( h ) {}

    void operator()( const tbb::blocked_range<size_t>& r ) const {
        frantic::graphics::vector3f centroid;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            centroid = ( m_vertices[m_faces[i].x] + m_vertices[m_faces[i].y] + m_vertices[m_faces[i].z] ) / 3;
            m_gradients[i] = m_isp.get_gradient( centroid, h );
        }
    }
};

/**
 * For every vertex in the mesh, find all its neighbouring faces.
 */
void find_adjacent_faces( std::vector<std::vector<int>>& adjFaces, const frantic::geometry::trimesh3& mesh );

/**
 * For every vertex in the mesh, find all its neighbouring vertices.
 */
void find_adjacent_vertices( std::vector<std::vector<int>>& connections, frantic::geometry::trimesh3& mesh );

/**
 * Finds the gradient at the center point of each face and stores it in 'gradients'.
 */
template <class ImplicitSurfacePolicy>
void find_face_gradients( std::vector<frantic::geometry::vector3f>& gradients,
                          std::vector<frantic::geometry::vector3>& faces,
                          std::vector<frantic::graphics::vector3f>& vertices, ImplicitSurfacePolicy& isp, float h ) {
    tbb::task_scheduler_init taskScheduleInit;
    size_t numFaces = faces.size();

    gradients.clear();
    gradients.resize( numFaces );

    tbb::blocked_range<size_t> range( 0, numFaces );
    apply_find_gradients<ImplicitSurfacePolicy> body( gradients, faces, vertices, isp, h );
#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body );
#else
    body( range );
#endif
}

/**
 * The N operator fits the normals of the mesh to the normals of the implicit surface. Equation 5.5 (not labeled, in
 * between 5 and 6) in [3].
 */
frantic::graphics::vector3f n_force( std::vector<frantic::graphics::vector3f>& vertices,
                                     std::vector<frantic::graphics::vector3>& faces, std::vector<int>& adjFaces,
                                     int pos, std::vector<frantic::graphics::vector3f>& gradients );

/**
 * The R operator equalizes the mesh sampling rate with equation (6) from [3]. Used for the mesh optimization to an
 * implicit surface.
 */
frantic::graphics::vector3f r_force( std::vector<frantic::graphics::vector3f>& vertices, std::vector<int>& adjVertices,
                                     frantic::graphics::vector3f& normal, int pos, float C = 0.1f );

} // namespace detail

/**
 * Laplacian smoothing. Equation 2 from [1]. Points move towards the median point of all neighbours multiplied by a
 * weighted amount. Reduces high frequency surface information and tends to flatten the surface. Works on the mesh with
 * holes.
 */
void laplacian_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount, float weight );
void laplacian_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount, float weight,
                       const std::vector<size_t>& subset );

void laplacian_smooth( frantic::geometry::mesh_interface* mesh, std::size_t iterationCount, double scale );
void laplacian_smooth( frantic::geometry::mesh_interface* mesh, const frantic::geometry::mesh_channel* channel,
                       std::size_t iterationCount, double scale );

/**
 * Uses the discrete mean curvature flow for smoothing and the Laplacian flow for improving mesh sampling rate. Equation
 * 11 from [1]. Scale / step should also be a small positive number. Gives undesired results for mesh with 'bubbles'
 * within the surface.
 */
void laplacian_mean_curvature_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount, float scale,
                                      float constant );

/**
 * Smoothing that examines the angle between the mean curvature vector and normal speed. Smoothing as suggested in
 * equation 13 from [1]. Epsilon should be small, .1 is good. Scale / step should also be a small positive number. Gives
 * undesired results for mesh with 'bubbles' within the surface.
 */
void adaptive_regularizing_mean_curvature_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount,
                                                  float scale, float elipson );

/**
 * Accurately evolves mesh toward an implicit surface. Performs three subsequent subsection vector functions, n-flow
 * corrects the mesh normals to the implicit surface, z-flow pushes mesh vertices toward the implicit surface, and
 * r-flow equalizes the mesh sampling rate. Defined in equations (1-4) in [3].
 */
template <class ImplicitSurfacePolicy>
void evolve_mesh_to_implicit_surface( frantic::geometry::trimesh3& mesh, ImplicitSurfacePolicy& isp, int iterations,
                                      float h, float C, frantic::logging::progress_logger* progressLogger ) {
    tbb::task_scheduler_init taskScheduleInit;

    std::vector<std::vector<int>> adjFaces;
    std::vector<std::vector<int>> adjVertices;
    detail::find_adjacent_faces( adjFaces, mesh );
    detail::find_adjacent_vertices( adjVertices, mesh );

    std::vector<frantic::graphics::vector3f>& currVertices = mesh.vertices_ref();
    std::vector<frantic::graphics::vector3f> prevVertices = mesh.vertices_ref();
    size_t numVertices = mesh.vertices_ref().size();

    std::vector<frantic::geometry::vector3>& faces = mesh.faces_ref();
    std::vector<frantic::graphics::vector3f> faceGradients;

    int progress = 0;
    int finalProgress = iterations * ( C ? 3 : 2 );

    for( int i = 0; i < iterations; ++i ) {
        { // scope for range and body
            tbb::blocked_range<size_t> range( 0, numVertices );
            detail::apply_z_force<ImplicitSurfacePolicy> body( currVertices, prevVertices, faces, adjFaces, isp, h );
#ifndef FRANTIC_DISABLE_THREADS
            tbb::parallel_for( range, body );
#else
            body( range );
#endif
        }
        prevVertices = currVertices;
        progressLogger->update_progress( 100.f * ++progress / finalProgress );

        if( C != 0 ) {
            mesh.build_vertex_normals();
            trimesh3_vertex_channel_accessor<vector3f> normals =
                mesh.get_vertex_channel_accessor<vector3f>( _T("Normal") );

            for( int v = 0; v < numVertices; ++v ) {
                currVertices[v] = prevVertices[v] + detail::r_force( prevVertices, adjVertices[v], normals[v], v, C );
            }
            prevVertices = currVertices;
            progressLogger->update_progress( 100.f * ++progress / finalProgress );
        }

        detail::find_face_gradients<ImplicitSurfacePolicy>( faceGradients, faces, currVertices, isp, h );

        { // scope for range and body
            tbb::blocked_range<size_t> range( 0, numVertices );
            detail::apply_n_force body( currVertices, prevVertices, faces, adjFaces, faceGradients );
#ifndef FRANTIC_DISABLE_THREADS
            tbb::parallel_for( range, body );
#else
            body( range );
#endif
        }
        prevVertices = currVertices;
        progressLogger->update_progress( 100.f * ++progress / finalProgress );
    }
}

} // namespace relaxation
} // namespace geometry
} // namespace frantic
