// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/dcel_construction.hpp>

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/dcel_iterators.hpp>

#include <frantic/graphics/vector3.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>

#include <vector>

namespace frantic {
namespace geometry {

// TODO: make a parallel version
void trimesh3_to_dcel( const trimesh3& mesh, dcel& out ) {
    out.initialize( mesh.vertex_count(), mesh.face_count() );

    for( size_t i = 0; i < mesh.face_count(); ++i ) {
        out.add_face3( i, mesh.get_face( i )[0], mesh.get_face( i )[1], mesh.get_face( i )[2] );
    }
}

void dcel_to_trimesh3( const dcel& input, const std::vector<frantic::graphics::vector3f>& vertices,
                       trimesh3& outMesh ) {
    using namespace frantic::graphics;

    if( !input.is_triangle_mesh() ) {
        throw std::runtime_error( "dcel_to_trimesh3 : Error, mesh was not triangulated." );
    }

    if( &outMesh.vertices_ref() != &vertices ) {
        outMesh.vertices_ref() = vertices;
    }

    outMesh.faces_ref().resize( input.face_count() );

    for( size_t i = 0; i < input.face_count(); ++i ) {
        size_t j = 0;
        BOOST_FOREACH( dcel::const_halfedge_handle handle, dcel_face_cycle_range( input.get_face_halfedge( i ) ) ) {
            if( j >= 3 ) {
                throw std::runtime_error( "dcel_to_trimesh3 : Error, mesh was not triangulated." );
            }

            outMesh.get_face( i )[j] = handle.source_vertex();
            ++j;
        }
    }
}

void polymesh3_to_dcel( const_polymesh3_ptr polymesh, dcel& out ) {
    out.initialize( polymesh->vertex_count(), polymesh->face_count() );

    polymesh3_const_vertex_accessor_base acc = polymesh->get_const_vertex_accessor( _T( "verts" ) );

    for( size_t i = 0; i < polymesh->face_count(); ++i ) {
        out.add_face( i, acc.get_face( i ).first, acc.get_face( i ).second );
    }
}

void mesh_interface_to_dcel( const mesh_interface_ptr meshInterface, dcel& out ) {
    mesh_interface_to_dcel( meshInterface.get(), out );
}

void mesh_interface_to_dcel( const mesh_interface* meshInterface, dcel& out ) {
    out.initialize( meshInterface->get_num_verts(), meshInterface->get_num_faces() );

    for( size_t i = 0; i < meshInterface->get_num_faces(); ++i ) {
        out.begin_face( i );
        for( size_t j = 0; j < meshInterface->get_num_face_verts( i ); ++j ) {
            out.add_face_vertex( meshInterface->get_face_vert_index( i, j ) );
        }
        out.end_face();
    }
}

void mesh_channel_to_dcel( const mesh_channel* meshChannel, dcel& out ) {
    if( meshChannel->get_channel_type() != mesh_channel::face_vertex ) {
        throw std::runtime_error(
            "mesh_interface_channel_to_dcel: Error, cannot get a dcel for a channel without custom faces." );
    }

    out.initialize( meshChannel->get_num_elements(), meshChannel->get_num_faces() );

    for( size_t i = 0; i < meshChannel->get_num_faces(); ++i ) {
        out.begin_face( i );
        for( size_t j = 0; j < meshChannel->get_num_face_verts( i ); ++j ) {
            out.add_face_vertex( meshChannel->get_fv_index( i, j ) );
        }
        out.end_face();
    }
}

void parameterization_chart_to_dcel( const parameterization::parameterization_chart& chart, dcel& out ) {
    out.initialize( chart.get_num_vertices(), chart.get_num_faces() );

    for( size_t i = 0; i < chart.get_num_faces(); ++i ) {
        out.begin_face( i );
        for( size_t j = 0; j < chart.get_num_face_verts( i ); ++j ) {
            out.add_face_vertex( chart.get_face_vert_index( i, j ) );
        }
        out.end_face();
    }
}

} // namespace geometry
} // namespace frantic
