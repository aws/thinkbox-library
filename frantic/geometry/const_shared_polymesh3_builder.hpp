// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/polymesh3.hpp>

namespace frantic {
namespace geometry {

/**
 *  Build a const_polymesh3_ptr from data that can be shared with other
 * meshes.
 *
 * @note This class is intended for internal use.  Users should rarely
 *       interact with this class.
 */
class const_shared_polymesh3_builder {
  public:
    /**
     *  Initialize a polymesh3 using vertex and face information that can be
     * shared with other meshes.
     *
     * @note You must not modify any of the channel data or faces after
     *       passing them into this function.
     *
     * @param vertData the vertex positions.  This should be an array of
     *                 vector3f.
     * @param geomPolyIndices the vertex indices for all faces in the mesh.
     * @param geomPolyEndIndices for every face in the mesh, this is one past
     *                           the last index in geomPolyIndices used for
     *                           the face.
     */
    const_shared_polymesh3_builder( polymesh3_channel_data& vertData, polymesh3_channel_faces& geomPolyIndices,
                                    polymesh3_channel_faces& geomPolyEndIndices );

    /**
     *  Add a vertex channel using data that can be shared with other meshes.
     *
     * @note You must not modify any of the channel data or faces after
     *       passing them into this function.
     *
     * @param channel name of the channel to add.
     * @param type primitive data type of the channel.
     * @param arity number of primitive data items per vertex.
     * @param inVertexBuffer the channel data.
     * @param pInFaceBuffer optional custom face indices.
     */
    void add_vertex_channel( const frantic::tstring& channel, polymesh3_channel_data& inVertexBuffer,
                             polymesh3_channel_faces* pInFaceBuffer = (polymesh3_channel_faces*)NULL );

    /**
     *  Add a face channel using data that can be shared with other meshes.
     *
     * @note You must not modify any of the channel data after passing it into
     *       this function.
     *
     * @param channel name of the channel to add.
     * @param type primitive data type of the channel.
     * @param arity number of primitive data items per vertex.
     * @param inFaceBuffer the channel data.
     */
    void add_face_channel( const frantic::tstring& channel, polymesh3_channel_data& inFaceBuffer );

    /**
     * @return the mesh built from data passed into this class.
     */
    const_polymesh3_ptr finalize();

  private:
    polymesh3_ptr m_mesh;
};

} // namespace geometry
} // namespace frantic
