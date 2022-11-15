// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include "utilities/mesh_generators.hpp"

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/dcel_iterators.hpp>
#include <frantic/geometry/dcel_surface_iterators.hpp>
#include <frantic/geometry/trimesh3.hpp>

#include <frantic/graphics/vector3.hpp>

#include <frantic/misc/iterator.hpp>

#include <boost/foreach.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>

namespace {

void check_face_consistency( const frantic::geometry::dcel& structure ) {
    using namespace frantic::geometry;

    for( size_t i = 0; i < structure.face_count(); ++i ) {
        dcel::const_halfedge_handle handle = structure.get_face_halfedge( i );
        if( handle.is_valid() ) {
            EXPECT_EQ( i, handle.current_face() );
        }
    }
}

void check_vertex_consistency( const frantic::geometry::dcel& structure ) {
    using namespace frantic::geometry;

    for( size_t i = 0; i < structure.vertex_count(); ++i ) {
        dcel::const_halfedge_handle handle = structure.get_vertex_halfedge( i );
        if( handle.is_valid() ) {
            EXPECT_EQ( i, handle.target_vertex() );
        }
    }
}

// This is assuming our new convention of ensuring
//   edgeId == floor(halfedgeId / 2)
void check_halfedge_index_consistency( const frantic::geometry::dcel& structure ) {

    for( size_t i = 0; i < structure.halfedge_count(); ++i ) {
        EXPECT_EQ( i + ( i % 2 == 0 ? 1 : -1 ), structure.get_halfedge( i ).twin().get_index() );
    }
}

template <class InputIterator>
void face_consistency_check( const frantic::geometry::dcel& structure, InputIterator begin, InputIterator end ) {
    using namespace frantic::geometry;

    std::vector<size_t> faceVerts( begin, end );

    dcel::const_halfedge_handle startFaceEdge = structure.find_vertex_halfedge( faceVerts.back(), faceVerts.front() );

    ASSERT_TRUE( startFaceEdge.is_valid() );

    size_t currentFaceVertex = 0;

    // Check the circular order of the newly added face
    BOOST_FOREACH( dcel::const_halfedge_handle it, dcel_face_range( startFaceEdge, startFaceEdge ) ) {
        EXPECT_EQ( faceVerts[currentFaceVertex], it.target_vertex() );
        ++currentFaceVertex;
    }

    EXPECT_EQ( currentFaceVertex, faceVerts.size() );

    // Check that the vertex range is correct over every corner of the new face
    for( size_t i = 0; i < faceVerts.size(); ++i ) {
        const size_t current = faceVerts[i];
        const size_t before = faceVerts[( ( i + faceVerts.size() ) - 1 ) % faceVerts.size()];
        const size_t after = faceVerts[( i + 1 ) % faceVerts.size()];

        dcel::const_halfedge_handle beforeEdge = structure.find_vertex_halfedge( before, current );
        dcel::const_halfedge_handle afterEdge = structure.find_vertex_halfedge( after, current );

        ASSERT_TRUE( beforeEdge.is_valid() );
        ASSERT_TRUE( afterEdge.is_valid() );
        EXPECT_EQ( afterEdge, beforeEdge.vertex_next() );
        // EXPECT_EQ( beforeEdge, afterEdge.vertex_prev() );
    }
}

void rotate_face( frantic::graphics::vector3& face, size_t rot ) { std::rotate( &face[0], &face[rot % 3], &face[3] ); }

size_t trimesh3_permutation_count( const frantic::geometry::trimesh3& mesh ) {
    if( mesh.face_count() > 18 ) {
        throw std::runtime_error(
            "trimesh3_face_permutation_count -- This mesh has more than 18 faces, are you sure you "
            "want to compute 3^18 iterations?" );
    }

    size_t count = 1;
    for( size_t i = 0; i < mesh.face_count(); ++i ) {
        count *= 3;
    }
    return count;
}

void permute_trimesh3( frantic::geometry::trimesh3& dest, size_t permutation ) {
    for( size_t i = 0; i < dest.face_count(); ++i ) {
        rotate_face( dest.get_face( i ), permutation % 3 );
        permutation /= 3;
        if( permutation == 0 ) {
            break;
        }
    }
}

/// A method to continuously check consistency as we incrementally build a dcel from a trimesh
void create_halfedge_structure_from_trimesh3( const frantic::geometry::trimesh3& mesh,
                                              frantic::geometry::dcel& outStructure ) {
    outStructure.initialize( mesh.vertex_count(), mesh.face_count() );

    for( size_t i = 0; i < mesh.face_count(); ++i ) {
        size_t indices[3] = { (size_t)mesh.get_face( i )[0], (size_t)mesh.get_face( i )[1],
                              (size_t)mesh.get_face( i )[2] };
        outStructure.add_face( i, indices, indices + 3 );

        // check this and all previous faces for consistency
        for( size_t j = 0; j <= i; ++j ) {
            size_t subIndices[3] = { (size_t)mesh.get_face( j )[0], (size_t)mesh.get_face( j )[1],
                                     (size_t)mesh.get_face( j )[2] };
            face_consistency_check( outStructure, subIndices, subIndices + 3 );
        }
    }

    // Check that the transfer was successful
    EXPECT_EQ( mesh.vertex_count(), outStructure.vertex_count() );
    EXPECT_EQ( mesh.face_count(), outStructure.face_count() );
    EXPECT_LE( mesh.face_count() * 3, outStructure.halfedge_count() );
    EXPECT_TRUE( outStructure.is_triangle_mesh() );
}

} // namespace

TEST( DCEL, MakeFromTetrahedron ) {
    using namespace frantic::geometry;

    trimesh3 tetrahedronMesh;

    make_regular_tetrahedron( tetrahedronMesh );

    const size_t numPerms = trimesh3_permutation_count( tetrahedronMesh );

    for( size_t perm = 0; perm < numPerms; ++perm ) {
        dcel halfedgeStructure;
        trimesh3 copyMesh( tetrahedronMesh );
        permute_trimesh3( copyMesh, perm );

        create_halfedge_structure_from_trimesh3( copyMesh, halfedgeStructure );

        // tetrahedron specific: check that the circular order around all vertices is the set of all other vertices
        for( size_t i = 0; i < copyMesh.vertex_count(); ++i ) {
            std::set<size_t> vertexSet;
            for( size_t j = 0; j < 4; ++j ) {
                if( j != i ) {
                    vertexSet.insert( j );
                }
            }

            std::set<size_t> halfedgeVertexSet;
            BOOST_FOREACH( dcel::halfedge_handle he,
                           dcel_vertex_cycle_range( halfedgeStructure.get_vertex_halfedge( i ) ) ) {
                halfedgeVertexSet.insert( he.source_vertex() );
            }

            EXPECT_EQ( vertexSet, halfedgeVertexSet );
        }
    }
}

TEST( DCEL, IncrementalFanConstructionSimple ) {
    using namespace frantic::geometry;

    frantic::geometry::trimesh3 mesh;

    // actual geometry is irrelevant
    for( size_t i = 0; i < 8; ++i ) {
        mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    }

    int faces[5][3] = {
        { 0, 1, 2 }, { 3, 2, 4 }, { 5, 2, 6 }, { 1, 6, 2 }, { 5, 7, 2 },
    };

    for( size_t i = 0; i < 5; ++i ) {
        mesh.add_face( faces[i][0], faces[i][1], faces[i][2] );
    }

    const size_t numPerms = trimesh3_permutation_count( mesh );

    for( size_t perm = 0; perm < numPerms; ++perm ) {
        dcel halfedgeStructure;
        trimesh3 copyMesh( mesh );
        permute_trimesh3( copyMesh, perm );
        create_halfedge_structure_from_trimesh3( copyMesh, halfedgeStructure );

        size_t expectedIncidence[7] = { 0, 3, 4, 7, 5, 6, 1 };

        const dcel::halfedge_handle startPoint = halfedgeStructure.find_vertex_halfedge( 0, 2 );

        ASSERT_TRUE( startPoint.is_valid() );

        size_t currentIndex = 0;

        BOOST_FOREACH( dcel::halfedge_handle edge, dcel_vertex_range( startPoint, startPoint ) ) {
            ASSERT_LT( currentIndex, 7 );
            EXPECT_EQ( expectedIncidence[currentIndex], edge.source_vertex() );
            ++currentIndex;
        }

        EXPECT_EQ( 7, currentIndex );
    }
}

TEST( DCEL, OctofanPermutations ) {
    using namespace frantic::geometry;

    trimesh3 octofanMesh;
    make_octofan_mesh( octofanMesh );

    const size_t numTests = 13;

    // ideally we could run every permutation (as I believe that would cover every possible
    // problem input sequence), but that would take waaaay to long.
    size_t testPermutations[numTests][8] = {
        {
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
        }, // baseline test
        {
            0,
            4,
            2,
            6,
            7,
            5,
            3,
            1,
        }, // test 'leading' insertions with correction
        {
            0,
            4,
            6,
            2,
            7,
            5,
            3,
            1,
        },
        {
            0,
            4,
            2,
            6,
            5,
            7,
            1,
            3,
        }, // test 'following' insertions with correction
        {
            0,
            4,
            6,
            2,
            5,
            7,
            1,
            3,
        },
        {
            0,
            1,
            2,
            6,
            4,
            7,
            3,
            5,
        }, // test simple (1-entry) cycle correction
        {
            0,
            1,
            2,
            6,
            4,
            7,
            3,
            5,
        },
        {
            0,
            2,
            3,
            5,
            6,
            1,
            7,
            4,
        }, // test complex cycle correction
        {
            0,
            2,
            3,
            5,
            6,
            7,
            1,
            4,
        },
        {
            0,
            5,
            6,
            2,
            3,
            1,
            7,
            4,
        },
        {
            0,
            5,
            6,
            2,
            3,
            7,
            1,
            4,
        },
        {
            0,
            2,
            3,
            5,
            6,
            4,
            7,
            1,
        }, // test continued cycle correction
        {
            0,
            2,
            3,
            5,
            6,
            4,
            1,
            7,
        },
    };

    for( size_t currentTest = 0; currentTest < numTests; ++currentTest ) {
        trimesh3 octofanCopy( octofanMesh );
        for( size_t i = 0; i < 8; ++i ) {
            octofanCopy.get_face( i ) = octofanMesh.get_face( testPermutations[currentTest][i] );
        }

        // Adding the permutations is a really good test, but it makes the running time a little too slow
        // for( size_t perm = 0; perm < numVertexPerms; ++perm ) {
        trimesh3 copyMesh( octofanCopy );
        // permute_trimesh3( octofanCopy, perm );
        frantic::geometry::dcel halfedgeStructure;
        create_halfedge_structure_from_trimesh3( copyMesh, halfedgeStructure );
        //}
    }
}

TEST( DCEL, ThrowOnSharedApexNonManifold ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    dcel halfedgeStructure;

    make_hourglass_mesh( mesh );

    // The mesh is non-manifold, and therefore should generate an error when trying to construct a dcel from it
    EXPECT_ANY_THROW( create_halfedge_structure_from_trimesh3( mesh, halfedgeStructure ) );
}

TEST( DCEL, CollapseEdge ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    dcel halfedgeStructure;

    make_sphere_mesh( 50, mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );

    dcel::halfedge_handle handle = halfedgeStructure.get_halfedge( 0 );
    dcel::halfedge_handle oldTwin = handle.twin();

    dcel::index_t leftFace = handle.current_face();
    dcel::index_t rightFace = handle.opposite_face();
    dcel::index_t removedVertex = handle.source_vertex();
    dcel::index_t remainingVertex = handle.target_vertex();

    EXPECT_TRUE( handle.is_valid() );
    EXPECT_TRUE( oldTwin.is_valid() );
    EXPECT_TRUE( halfedgeStructure.try_collapse_edge( handle ) );
    check_face_consistency( halfedgeStructure );
    check_vertex_consistency( halfedgeStructure );
    check_halfedge_index_consistency( halfedgeStructure );
    EXPECT_FALSE( handle.is_valid() );
    EXPECT_FALSE( oldTwin.is_valid() );

    EXPECT_EQ( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_vertex_halfedge( removedVertex ) );
    EXPECT_NE( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_vertex_halfedge( remainingVertex ) );
    EXPECT_EQ( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_face_halfedge( leftFace ) );
    EXPECT_EQ( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_face_halfedge( rightFace ) );

    size_t i = 0;
    while( i < halfedgeStructure.halfedge_count() && !handle.is_valid() ) {
        handle = halfedgeStructure.get_halfedge( i );
        ++i;
    }

    EXPECT_TRUE( handle.is_valid() );
    EXPECT_TRUE( handle.twin().is_valid() );
}

TEST( DCEL, CollapseBoundaryEdge ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    dcel halfedgeStructure;

    make_tetrahedral_cap( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );

    dcel::halfedge_handle boundaryEdge = halfedgeStructure.find_vertex_halfedge( 0, 1 );

    EXPECT_TRUE( halfedgeStructure.check_manifold_collapse( boundaryEdge ) );
    EXPECT_TRUE( halfedgeStructure.try_collapse_edge( boundaryEdge ) );
    check_face_consistency( halfedgeStructure );
    check_vertex_consistency( halfedgeStructure );
    check_halfedge_index_consistency( halfedgeStructure );

    // This collapse will actually remove the boundary itself and close the manifold
    for( size_t i = 0; i < halfedgeStructure.halfedge_count(); ++i ) {
        dcel::halfedge_handle handle = halfedgeStructure.get_halfedge( i );
        if( handle.is_valid() ) {
            EXPECT_NE( dcel::INVALID_FACE_INDEX, handle.current_face() );
        }
    }
}

TEST( DCEL, CollapseToNonManifold ) {

    using namespace frantic::geometry;

    trimesh3 mesh;

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    mesh.add_vertex( 1.0f, 0.0f, 0.0f );
    mesh.add_vertex( 1.0f, 1.0f, 0.0f );
    mesh.add_vertex( 0.0f, 1.0f, 0.0f );

    mesh.add_face( 0, 1, 2 );
    mesh.add_face( 0, 2, 3 );

    add_hole_face( mesh, 0 );
    add_hole_face( mesh, 1 );

    dcel halfedgeStructure;
    trimesh3_to_dcel( mesh, halfedgeStructure );

    dcel::halfedge_handle crossEdge = halfedgeStructure.find_vertex_halfedge( 0, 2 );

    EXPECT_FALSE( halfedgeStructure.check_manifold_collapse( crossEdge ) );

    halfedgeStructure.flip_edge( crossEdge );

    EXPECT_FALSE( halfedgeStructure.check_manifold_collapse( crossEdge ) );
}

TEST( DCEL, BoundaryEdgeCollapseTrimesh ) {
    using namespace frantic::geometry;

    trimesh3 mesh;

    make_quad_trimesh( mesh );

    dcel halfedgeStructure;
    trimesh3_to_dcel( mesh, halfedgeStructure );

    dcel::halfedge_handle boundaryEdge = halfedgeStructure.find_vertex_halfedge( 0, 1 );

    bool isManifoldCollapse = halfedgeStructure.check_manifold_collapse( boundaryEdge );

    EXPECT_TRUE( isManifoldCollapse );

    halfedgeStructure.collapse_edge( boundaryEdge );
    check_face_consistency( halfedgeStructure );
    check_vertex_consistency( halfedgeStructure );
    check_halfedge_index_consistency( halfedgeStructure );

    dcel::halfedge_handle otherEdge = halfedgeStructure.find_vertex_halfedge( 2, 3 );

    EXPECT_EQ( 3, std::distance( dcel_face_begin( otherEdge ), dcel_face_end( otherEdge ) ) );

    EXPECT_EQ( 3, std::distance( dcel_face_begin( otherEdge.twin() ), dcel_face_end( otherEdge.twin() ) ) );
}

TEST( DCEL, BoundaryEdgeCollapsePolymesh ) {
    using namespace frantic::geometry;

    polymesh3_ptr mesh = make_quad_polymesh();

    dcel halfedgeStructure;
    polymesh3_to_dcel( mesh, halfedgeStructure );

    dcel::halfedge_handle boundaryEdge = halfedgeStructure.find_vertex_halfedge( 0, 1 );

    EXPECT_TRUE( halfedgeStructure.check_manifold_collapse( boundaryEdge ) );

    halfedgeStructure.collapse_edge( boundaryEdge );
    check_face_consistency( halfedgeStructure );
    check_vertex_consistency( halfedgeStructure );
    check_halfedge_index_consistency( halfedgeStructure );

    dcel::halfedge_handle otherEdge = halfedgeStructure.find_vertex_halfedge( 2, 3 );

    EXPECT_EQ( 3, std::distance( dcel_face_begin( otherEdge ), dcel_face_end( otherEdge ) ) );

    EXPECT_EQ( 3, std::distance( dcel_face_begin( otherEdge.twin() ), dcel_face_end( otherEdge.twin() ) ) );
}

TEST( DCEL, InternalBoundaryCollapse ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    trimesh3 mesh;

    mesh.add_vertex( 0.0f, 0.0f, 0.0f );
    mesh.add_vertex( 1.0f, 0.0f, 0.0f );
    mesh.add_vertex( 0.0f, 1.0f, 0.0f );

    mesh.add_face( 0, 1, 2 );

    vector3 hole = add_hole_face( mesh, 0 );

    dcel halfedgeStructure;
    trimesh3_to_dcel( mesh, halfedgeStructure );

    dcel::halfedge_handle boundaryEdge = halfedgeStructure.find_vertex_halfedge( hole[0], hole[1] );

    const dcel::index_t foldVertex1 = boundaryEdge.target_vertex();

    EXPECT_TRUE( halfedgeStructure.check_manifold_collapse( boundaryEdge ) );

    halfedgeStructure.collapse_edge( boundaryEdge );
    check_face_consistency( halfedgeStructure );
    check_vertex_consistency( halfedgeStructure );
    check_halfedge_index_consistency( halfedgeStructure );

    dcel::halfedge_handle remantHalfedge = halfedgeStructure.find_vertex_halfedge( foldVertex1, hole[2] );

    const dcel::index_t foldVertex2 = remantHalfedge.target_vertex();

    EXPECT_TRUE( halfedgeStructure.check_manifold_collapse( remantHalfedge ) );

    halfedgeStructure.collapse_edge( remantHalfedge );
    check_face_consistency( halfedgeStructure );
    check_vertex_consistency( halfedgeStructure );
    check_halfedge_index_consistency( halfedgeStructure );

    std::vector<dcel::index_t> surroundingVerts(
        dcel_vertex_adjacency_begin( halfedgeStructure.get_vertex_halfedge( foldVertex2 ) ),
        dcel_vertex_adjacency_end( halfedgeStructure.get_vertex_halfedge( foldVertex2 ) ) );

    std::sort( surroundingVerts.begin(), surroundingVerts.end() );

    EXPECT_EQ( 3, surroundingVerts.size() );

    for( size_t i = 0; i < 3; ++i ) {
        EXPECT_EQ( i, surroundingVerts[i] );
    }
}

TEST( DCEL, TetrahedralCapCollapseEdge ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    dcel halfedgeStructure;

    make_tetrahedral_cap( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );

    // the apex is 6, surrounded by 3, 4, 5

    dcel::halfedge_handle invalidCollapseEdge = halfedgeStructure.find_vertex_halfedge( 3, 4 );

    EXPECT_FALSE( halfedgeStructure.check_manifold_collapse( invalidCollapseEdge ) );
    EXPECT_FALSE( halfedgeStructure.try_collapse_edge( invalidCollapseEdge ) );

    dcel::halfedge_handle validCollapseEdge = halfedgeStructure.find_vertex_halfedge( 6, 4 );

    dcel::index_t leftFace = validCollapseEdge.current_face();
    dcel::index_t rightFace = validCollapseEdge.opposite_face();
    dcel::index_t removedVertex = validCollapseEdge.source_vertex();
    dcel::index_t remainingVertex = validCollapseEdge.target_vertex();

    EXPECT_TRUE( halfedgeStructure.check_manifold_collapse( validCollapseEdge ) );
    EXPECT_TRUE( halfedgeStructure.try_collapse_edge( validCollapseEdge ) );
    check_face_consistency( halfedgeStructure );
    check_vertex_consistency( halfedgeStructure );
    check_halfedge_index_consistency( halfedgeStructure );

    EXPECT_EQ( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_vertex_halfedge( removedVertex ) );
    EXPECT_NE( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_vertex_halfedge( remainingVertex ) );
    EXPECT_EQ( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_face_halfedge( leftFace ) );
    EXPECT_EQ( dcel::INVALID_HALFEDGE_HANDLE, halfedgeStructure.get_face_halfedge( rightFace ) );
}

TEST( DCEL, CheckBoundaryVertex ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    dcel halfedgeStructure;

    make_regular_tetrahedron( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );

    for( dcel::index_t i = 0; i < halfedgeStructure.vertex_count(); ++i ) {
        EXPECT_FALSE( halfedgeStructure.is_boundary_vertex( i ) );
    }

    make_quad_trimesh( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );

    for( dcel::index_t i = 0; i < 4; ++i ) {
        EXPECT_TRUE( halfedgeStructure.is_boundary_vertex( i ) );
    }

    for( size_t i = 0; i < 2; ++i ) {
        add_face_vertex( mesh, i );
    }

    trimesh3_to_dcel( mesh, halfedgeStructure );

    for( dcel::index_t i = 0; i < 4; ++i ) {
        EXPECT_TRUE( halfedgeStructure.is_boundary_vertex( i ) );
    }

    for( dcel::index_t i = 5; i < 6; ++i ) {
        EXPECT_FALSE( halfedgeStructure.is_boundary_vertex( i ) );
    }
}

TEST( DCEL, SimpleBoundaryInfo ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    dcel halfedgeStructure;

    make_regular_tetrahedron( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );
    halfedgeStructure.cache_boundary_info();

    EXPECT_TRUE( halfedgeStructure.boundary_info_cached() );
    EXPECT_EQ( 0, halfedgeStructure.boundary_count() );

    make_octofan_mesh( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );
    halfedgeStructure.cache_boundary_info();

    EXPECT_TRUE( halfedgeStructure.boundary_info_cached() );

    EXPECT_EQ( 1, halfedgeStructure.boundary_count() );

    dcel::const_halfedge_handle boudaryHandle( halfedgeStructure.get_boundary_halfedge( 0 ) );

    EXPECT_EQ( 8, std::distance( dcel_face_cycle_begin( boudaryHandle ), dcel_face_cycle_end( boudaryHandle ) ) );

    make_holed_tetrahedron( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );

    halfedgeStructure.cache_boundary_info();

    EXPECT_TRUE( halfedgeStructure.boundary_info_cached() );

    EXPECT_EQ( 4, halfedgeStructure.boundary_count() );

    for( size_t i = 0; i < halfedgeStructure.boundary_count(); ++i ) {
        dcel::const_halfedge_handle boudaryHandle( halfedgeStructure.get_boundary_halfedge( i ) );
        EXPECT_EQ( 3, std::distance( dcel_face_cycle_begin( boudaryHandle ), dcel_face_cycle_end( boudaryHandle ) ) );
    }
}

TEST( DCEL, VertexSurfaceRange ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    dcel halfedgeStructure;

    make_regular_tetrahedron( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );
    halfedgeStructure.cache_boundary_info();

    for( size_t i = 0; i < 4; ++i ) {
        dcel::halfedge_handle boundary = next_boundary_edge( halfedgeStructure.get_vertex_halfedge( i ) );
        EXPECT_EQ( dcel::INVALID_HALFEDGE_HANDLE, boundary );

        size_t count = 0;
        BOOST_FOREACH( dcel::halfedge_handle handle,
                       dcel_vertex_surface_range( halfedgeStructure.get_vertex_halfedge( i ) ) ) {
            EXPECT_NE( i, handle.source_vertex() );
            EXPECT_GT( 4u, handle.source_vertex() );
            ++count;
        }

        EXPECT_EQ( 3, count );
    }

    make_holed_tetrahedron( mesh );
    trimesh3_to_dcel( mesh, halfedgeStructure );
    halfedgeStructure.cache_boundary_info();

    for( size_t i = 4; i < mesh.vertex_count(); ++i ) {
        dcel::halfedge_handle boundary = next_boundary_edge( halfedgeStructure.get_vertex_halfedge( i ) );
        EXPECT_NE( dcel::INVALID_HALFEDGE_HANDLE, boundary );

        size_t count = 0;
        BOOST_FOREACH( dcel::halfedge_handle handle,
                       dcel_vertex_surface_range( halfedgeStructure.get_vertex_halfedge( i ) ) ) {
            EXPECT_NE( i, handle.source_vertex() );
            EXPECT_FALSE( handle.is_boundary_face() );
            ++count;
        }

        EXPECT_EQ( 3, count );
    }

    for( size_t i = 0; i < halfedgeStructure.boundary_count(); ++i ) {
        BOOST_FOREACH( dcel::halfedge_handle handle,
                       dcel_face_cycle_range( halfedgeStructure.get_boundary_halfedge( i ) ) ) {
            dcel::halfedge_handle boundary = next_boundary_edge( handle );
            EXPECT_NE( dcel::INVALID_HALFEDGE_HANDLE, boundary );

            size_t count = 0;
            BOOST_FOREACH( dcel::halfedge_handle vertexHandle, dcel_vertex_surface_range( handle ) ) {
                EXPECT_FALSE( vertexHandle.is_boundary_face() );
                ++count;
            }

            EXPECT_EQ( 3, count );
        }
    }
}

TEST( DCEL, FlipSingleEdge ) {
    using namespace frantic::geometry;
    using namespace frantic::graphics;

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    dcel edgeStructure;
    trimesh3_to_dcel( mesh, edgeStructure );

    dcel::halfedge_handle flipEdge = edgeStructure.find_vertex_halfedge( 0, 2 );

    edgeStructure.flip_edge( flipEdge );

    EXPECT_EQ( 1, std::min( flipEdge.source_vertex(), flipEdge.target_vertex() ) );
    EXPECT_EQ( 3, std::max( flipEdge.source_vertex(), flipEdge.target_vertex() ) );

    for( size_t i = 0; i < 2; ++i ) {

        vector3 face;

        if( flipEdge.source_vertex() == 1 ) {
            EXPECT_EQ( 0, flipEdge.face_next().target_vertex() );
            face = vector3( 3, 0, 1 );
        } else {
            EXPECT_EQ( 2, flipEdge.face_next().target_vertex() );
            face = vector3( 1, 2, 3 );
        }

        std::set<dcel::index_t> faceSet;

        BOOST_FOREACH( dcel::halfedge_handle handle, dcel_face_cycle_range( flipEdge ) ) {
            EXPECT_EQ( flipEdge.current_face(), handle.current_face() );
            faceSet.insert( handle.target_vertex() );
        }

        EXPECT_EQ( 3, faceSet.size() );

        for( size_t i = 0; i < 3; ++i ) {
            EXPECT_TRUE( faceSet.find( face[i] ) != faceSet.end() );
        }

        flipEdge = flipEdge.twin();
    }
}

TEST( DCEL, CannotFlipBoundaryEdge ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    make_quad_trimesh( mesh );

    dcel edgeStructure;
    trimesh3_to_dcel( mesh, edgeStructure );

    for( dcel::index_t halfedgeId = 0; halfedgeId < edgeStructure.halfedge_count(); ++halfedgeId ) {
        dcel::halfedge_handle flipEdge = edgeStructure.get_halfedge( halfedgeId );

        if( flipEdge.is_boundary_edge() ) {
            ASSERT_ANY_THROW( edgeStructure.flip_edge( flipEdge ) );
        }
    }
}

TEST( DCEL, CannotFlipTetrahedronEdge ) {
    using namespace frantic::geometry;

    trimesh3 mesh;
    make_regular_tetrahedron( mesh );

    dcel edgeStructure;
    trimesh3_to_dcel( mesh, edgeStructure );

    for( dcel::index_t halfedgeId = 0; halfedgeId < edgeStructure.halfedge_count(); ++halfedgeId ) {
        dcel::halfedge_handle flipEdge = edgeStructure.get_halfedge( halfedgeId );
        ASSERT_ANY_THROW( edgeStructure.flip_edge( flipEdge ) );
    }
}
