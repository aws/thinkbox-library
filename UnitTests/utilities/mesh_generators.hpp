// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file
 * @author Stephen Kiazyk
 *
 * The goal is to create a set of procedural mesh methods for use in testing.
 * Ideally this would be a separate library that could be included in live
 * projects as well so users could generate test data
 */
#pragma once

#include <vector>

#include <frantic/geometry/trimesh3.hpp>

#include <frantic/geometry/polymesh3.hpp>

#include <frantic/misc/indexer.hpp>

/**
 * A single triangle
 */
void make_single_triangle_mesh( frantic::geometry::trimesh3& mesh );

/**
 * Trimesh with two adjacent faces.
 *
 * 3---2
 * |1 /|
 * | / |
 * |/ 0|
 * 0---1
 *
 */
void make_quad_trimesh( frantic::geometry::trimesh3& outMesh );

/**
 * Unit tetrahedron centered at the origin
 *
 * The face and vertex indexing is such that given face `i`, the only
 * vertex (of the set {0,1,2,3}) _not_ on that face is vertex `i`.
 * Likewise, the only face that vertex `i` is not on is face `i`.
 */
void make_regular_tetrahedron( frantic::geometry::trimesh3& mesh );

frantic::geometry::polymesh3_ptr make_regular_tetrahedron();

/**
 * A mesh containing 3 sliver triangles and 1 zero area face, where the outer
 * edges is a triangle
 */
void make_triangle_beveled_vertex_mesh( frantic::geometry::trimesh3& mesh );

/**
 * A mesh containing 3 sliver triangles and 1 zero area face, where the outer
 * edges is a square
 */
void make_square_beveled_vertex_mesh( frantic::geometry::trimesh3& mesh );

/**
 * A tetrahedron mesh, except that each vertex is doubled, and each
 * edge has been replaced by a pair of faces forming a quadrilateral.
 */
void make_tetrahedron_beveled_edge_mesh( frantic::geometry::trimesh3& mesh );

/**
 * A mesh where there is a chain of 3 or more co-incident vertices
 * which are not all on the same face.
 */
void make_chained_combined_mesh( frantic::geometry::trimesh3& mesh );

/**
 * A mesh that forms an N-point 'wheel', such that all faces are incident to a single
 * vertex in the middle. The middle vertex is guaranteed to have index 0.
 *
 * @param wheelCount the number of exterior points, must be at least 3
 */
void make_wheel_mesh( frantic::geometry::trimesh3& mesh, size_t wheelCount );

/**
 * A mesh that forms an 8-point 'wheel', such that all faces are incident to a single
 * vertex in the middle. Equivalent to calling `make_wheel_mesh` with a parameter of 8.
 * The middle vertex is guaranteed to have index 0.
 */
void make_octofan_mesh( frantic::geometry::trimesh3& mesh );

/**
 * Just a pair of triangles, back-to-back. The simplest possible
 * manifold mesh, and often a cause of problems
 */
void make_doubly_covered_triangle( frantic::geometry::trimesh3& mesh );

/**
 * Simplest possible mesh to contain a saddle vertex (> 2pi surface angle)
 */
void make_saddle_vertex_mesh( frantic::geometry::trimesh3& mesh );

/**
 * A simple one-sided mesh with a peak in the middle
 */
void make_tetrahedral_cap( frantic::geometry::trimesh3& mesh );

/**
 * Creates a simple sphere mesh, centered at the origin and with radius 1.0
 *
 * @remark By simple, I mean it will not attempt to make a well-proportioned sphere, but
 *   rather will generate it as a level-set mesh using spherical coordinates
 *
 * @param tesselationFactor a parameter to control the number of vertices to generate. Specifically
 *   the complexity of the sphere will be quadratic in the tesselationFactor. Any value from 0 upwards
 *   will produce a valid mesh (subject of course to the upper constraints of the integer types involved).
 */
void make_sphere_mesh( size_t tesselationFactor, frantic::geometry::trimesh3& mesh );

/**
 * An hourglass is two otherwise independent pyramid shapes which share an apex vertex
 * It is a very simple example of a non-manifold mesh.
 * @param resulting mesh
 */
void make_hourglass_mesh( frantic::geometry::trimesh3& mesh );

/**
 * A sharkfin is a face that juts out from an edge of an otherwise manifold surface
 *
 * @param resulting mesh
 */
void make_sharkfin_mesh( frantic::geometry::trimesh3& mesh );

/**
 * Makes a 2x2 cube centered at the origin
 *
 *    6----------7
 *   /|         /|
 *  / |        / |
 * 2----------3  |
 * |  4-------|--5
 * | /        | /
 * |/         |/
 * 0----------1
 *
 */
void make_cube_mesh( frantic::geometry::trimesh3& mesh );

/**
 * According to Paul: The feature is a double-sided quad connecting two surfaces together -- imagine a flattened
 * drinking straw connecting two inflated balloons. In this case, I've made it about as minimal as possible, two
 * tetrahedra, with the flat quad connecting two of their edges
 *
 * @param mesh I drink your milkshake! I drink it up!
 */
void make_bubble_straw_mesh( frantic::geometry::trimesh3& mesh );

/**
 * Make a (closed) mesh with the specified genus.
 *
 * @param genus the number of 'donut holes' in the resulting mesh
 * @param mesh (output) the resulting mesh
 */
void make_genus_mesh( size_t genus, frantic::geometry::trimesh3& mesh );

/**
 * Adds a 'spike' to the middle of the specified face
 * @remark replaces the target face with 9 other faces (hence, increases the number of faces by 8)
 *
 * @param mesh the mesh to be modified
 * @param faceId the face to be modified. The original face will no longer be in the mesh, its id will be occupied by
 * one of the 9 new faces
 * @param offsetFactor a value in the range [0.0,1.0] which controls how far/close to the center the base vertices will
 * be
 * @param faceHeight the height at which to place the center peak vertex away from the face
 */
void add_spikey_face( frantic::geometry::trimesh3& mesh, size_t faceId, float offsetFactor = 0.5f,
                      float faceHeight = 0.5f );

/**
 * Adds a 'hole' to the middle of the specified face
 * @remark replaces the target face with 6 other faces (increases the number of faces by 5)
 *
 * @param mesh the mesh to be modified
 * @param faceId the face to be modified.
 * @param offsetFactor a value in the range [0.0,1.0] which controls how far/close to the center the hole will be
 * @return the 3 vertices surrounding the newly created hole
 */
frantic::graphics::vector3 add_hole_face( frantic::geometry::trimesh3& mesh, size_t faceId, float offsetFactor = 0.5f );

/**
 * Adds a single vertex to the specified barycentric location on the specific face, then connects all of its internal
 * vertices to it.
 *
 * @param mesh the mesh to be modified
 * @param faceId the face to be modified
 * @param baryCoords the location in the face to insert the vertex
 * @return the index of the newly created vertex
 */
frantic::geometry::trimesh3::index_t
add_face_vertex( frantic::geometry::trimesh3& mesh, size_t faceId,
                 const frantic::graphics::vector3f& baryCoords = frantic::graphics::vector3f( 1.0f / 3.0f, 1.0f / 3.0f,
                                                                                              1.0f / 3.0f ) );

/**
 * A unit tetrahedron with a 'spike' emanating from each face
 */
void make_spikey_tetrahedron( frantic::geometry::trimesh3& mesh, float offsetFactor = 0.5f, float faceHeight = 0.5f );

/**
 * A unit tetrahedron with a 'hole' on each face
 */
void make_holed_tetrahedron( frantic::geometry::trimesh3& mesh, float offsetFactor = 0.5f );

/**
 * Creates a simple mesh with a t-vertex (that is, a vertex which intersects a
 * point on an edge other than its endpoints.
 *
 * 0 --3---2
 *  \  |  /
 *   \ | /
 *     1
 *
 * (Where triangle (0,3,2) has zero area)
 *
 * @param mesh the generated mesh
 */
void make_t_vertex_mesh( frantic::geometry::trimesh3& mesh );

/**
 * Creates a mesh that lies on the XY plane, in the +/+ quadrant.
 * The mesh will be a grid of (xSize+1)*(ySize+1) vertices, with 2*xSize*ySize faces.
 *
 * The vertex and face numbering will follow row-major ordering, note that this
 * means indexing must treat Y as the most-signifigant axis
 *
 * Example: make_plane_mesh( 3, 3, mesh )
 *
 *    12--13--14--15
 *    |\13|\15|\17|
 *    | \ | \ | \ |
 *    |12\|14\|16\|
 *    8---9---10--11
 *    |\ 7|\ 9|\11|
 * ^  | \ | \ | \ |
 * |  |6 \|8 \|10\|
 * |  4---5---6---7
 * |  |\ 1|\ 3|\ 5|
 * |  | \ | \ | \ |
 * |  |0 \|2 \|4 \|
 * Y  0---1---2---3
 *
 *    X -------->
 *
 * @param xSize the number of faces in the x direction
 * @param ySize the number of faces in the y direction
 * @param xDelta the distance between vertices in the x direction
 * @param yDelta the distance between vertices in the y direction
 */
void make_plane_mesh( frantic::geometry::trimesh3& mesh, int xSize, int ySize, float xDelta = 1.0f,
                      float yDelta = 1.0f );

/**
 * Create a plane mesh, that also has a varying Z component. The Z values will vary within
 * `zDelta` of 0.0f, and will be selected psuedo-randomly using `zSeed`
 */
void make_terrain_mesh( frantic::geometry::trimesh3& mesh, int xSize, int ySize, float xDelta = 1.0f,
                        float yDelta = 1.0f, float zDelta = 1.0f, boost::uint32_t zSeed = 0 );

/**
 * Specific problem case for UV unwrapping, a slightly elongated tetrahedron
 */
void make_elongated_tetrahedron_mesh( frantic::geometry::trimesh3& mesh );

/**
 * Polymesh composed of a single slivered triangle which references the same vertex twice.
 */
frantic::geometry::polymesh3_ptr make_sliver_polymesh();

/**
 * Polymesh composed of a single, triangle face.
 */
frantic::geometry::polymesh3_ptr make_triangle_polymesh();

/**
 * Polymesh composed of a single, square face.
 */
frantic::geometry::polymesh3_ptr make_quad_polymesh();

/**
 * Polymesh composed of a square split into two triangle faces.
 */
frantic::geometry::polymesh3_ptr make_split_quad_polymesh();

/**
 * Polymesh of a 2x2 cube, centered at the origin, composed of six quads
 *
 * @param removedFaces (optional) Specify a subset of faces in [0..5] that you do not want to appear on the generated
 * mesh
 */
frantic::geometry::polymesh3_ptr make_cube_polymesh( const std::set<size_t>& removedFaces = std::set<size_t>() );

/**
 * Creates a simple sphere mesh, centered at the origin and with radius 1.0
 */
frantic::geometry::polymesh3_ptr make_sphere_polymesh( size_t tesselationFactor );

/**
 * Create a quad split down the middle, resulting in two meshes with matching boundaries.
 *
 * @param separation (optional) How far to separate the boundaries of the half-cubes
 */
std::vector<frantic::geometry::polymesh3_ptr> make_separated_quad( float separation );

/**
 * Create a quad split down the middle, resulting in one mesh with an edge down the center.
 */
frantic::geometry::polymesh3_ptr make_welded_quad();

/**
 * Create a quad split down the middle, resulting in one mesh with an two adjacent quads.
 */
frantic::geometry::polymesh3_ptr make_combined_quad();
