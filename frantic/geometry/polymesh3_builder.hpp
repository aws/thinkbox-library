// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/polymesh3.hpp>

namespace frantic {
namespace geometry {

/**
 * This class is an intermediate object used for constructing the geometry of a polymesh3.
 * It is intended to build the initial geometry and polygon connectivity in this object,
 * which will produce a polymesh3 when finished.
 */
class polymesh3_builder {
    frantic::graphics::raw_byte_buffer m_vertexBuffer;
    std::vector<int> m_faceBuffer;
    std::vector<int> m_faceEndOffsets;

    /**
     * This function will remove any polygons which are non-planar. Ie. All verts don't lie on a single plane.
     */
    void removeNonPlanarPolygons();

    /**
     * This function will replace any vertices which share the same location with a single copy. The polygons
     * will be re-indexed accordingly.
     */
    void removeCoincidentVertices();

    /**
     * This function will remove any parts of a polygon that repeat vertex indices.
     * Ex. Polygon with indices 123234 will become 1234
     */
    void removePolygonLoops();

  public:
    polymesh3_builder();

    /**
     * This function will create the polymesh3 with the data supplied to the polymesh3_builder.
     * It will optionally run several checks to verify cetains aspects of correctness  for the
     * polymesh3.
     * @param removeNonPlanar Will remove any non-planar polygons in the mesh.
     * @param removeCoincident Will cause all coincident vertices to be replaced by just one
     *                          vertex. The polygon indices will be re-ordered as required.
     * @param removePolyRepeats Will remove any regions in a polygon's perimeter that loop back on
     *                           itself. Ex. Polygon 123234 will become 1234
     * @todo Add more cleanup options
     * @return A new polymesh3 object initialized with the data here.
     */
    polymesh3_ptr finalize( bool removeNonPlanar = false, bool removeCoincident = false,
                            bool removePolyRepeats = false );

    /**
     * Will add the vertex [x,y,z] to the polymesh
     * @param x The x-coordinate of the vertex to add.
     * @param y The y-coordinate of the vertex to add.
     * @param z The z-coordinate of the vertex to add.
     */
    void add_vertex( float x, float y, float z );

    /**
     * Will add the vertex [p.x, p.y, p.z] to the polymesh
     * @param p The coordinates of the vertex to add.
     */
    void add_vertex( const frantic::graphics::vector3f& p );

    /**
     * Will add new vertices to the polymesh and and a new polygon that indexes these vertices
     * @param pts A list of points to add as vertices, and to make into a new polygon.
     */
    void add_polygon( const std::vector<frantic::graphics::vector3f>& pts );

    /**
     * @overload
     * @param pPts A list of floating point numbers, where every group of three is an individual
     *             new vertex.
     * @param count The number of points in the list.
     */
    void add_polygon( const frantic::graphics::vector3f* pPts, std::size_t count );

    /**
     * @overload
     * @param pPts A list of floating point numbers, where every group of three is an individual
     *             new vertex.
     * @param count The number of points in the list.
     */
    void add_polygon( const float ( *pPts )[3], std::size_t count );

    /**
     * @overload
     * This function will add an new polygon to the polymesh, indexing vertices aleady existing.
     * @param ptIndices A list of indices into the vertex array.
     */
    void add_polygon( const std::vector<int>& ptIndices );

    /**
     * @overload
     * This function will add an new polygon to the polymesh, indexing vertices aleady existing.
     * @param pIndices A list of indices into the vertex array.
     * @param count The number of indices in the list.
     */
    void add_polygon( const int pIndices[], std::size_t count );
};

} // namespace geometry
} // namespace frantic
