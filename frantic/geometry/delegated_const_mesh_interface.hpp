// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/mesh_interface.hpp>

namespace frantic {
namespace channels {

class channel_propagation_policy;

}
} // namespace frantic

namespace frantic {
namespace geometry {

/**
 * mesh_interface implementation that delegates its function calls to a delegate mesh_interface.
 * This class is intended to be overloaded by a "post-processing" mesh interface.
 *
 * @note The caller is responsible for maintaining the lifetime of the delegate mesh.
 *
 * @note You must not modify the delegate mesh or its channels during the lifetime of this object.
 *
 * @todo create a better system for this?
 */
class delegated_const_mesh_interface : public mesh_interface {
  public:
    /**
     *  Create a class instance.
     *
     *  Initially, the instance does not have any vertex or face channels.
     * You can propagate channels from the delegate mesh by calling
     * add_delegate_vertex_channel() and add_delegate_face_channel().
     *
     * @param meshInterface the mesh_interface to use as a delegate.
     *    The caller is responsible for maintaining the lifetime of the
     *    meshInterface and its channels.
     */
    delegated_const_mesh_interface( const mesh_interface* meshInterface );

    bool is_valid() const;
    bool request_channel( const frantic::tstring& channelName, bool vertexChannel, bool forOutput,
                          bool throwOnError = true );
    std::size_t get_num_verts() const;
    void get_vert( std::size_t index, float ( &outValues )[3] ) const;
    std::size_t get_num_faces() const;
    std::size_t get_num_face_verts( std::size_t faceIndex ) const;
    std::size_t get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const;
    void get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const;
    void get_face_verts( std::size_t faceIndex, float outValues[][3] ) const;
    std::size_t get_num_elements() const;
    std::size_t get_face_element_index( std::size_t faceIndex ) const;
    /**
     * Currently throws an exception to avoid modifying the delegate mesh.
     * @todo implement this?
     */
    void init_adjacency();
    bool has_adjacency() const;
    bool init_vertex_iterator( frantic::geometry::vertex_iterator& vIt, std::size_t vertexIndex ) const;
    bool advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const;
    std::size_t get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const;
    std::size_t get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const;
    std::size_t get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const;
    bool is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const;
    bool is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const;
    void init_face_iterator( frantic::geometry::face_iterator& fIt, std::size_t faceIndex ) const;
    bool advance_face_iterator( frantic::geometry::face_iterator& fIt ) const;
    std::size_t get_face_neighbor( frantic::geometry::face_iterator& fIt ) const;
    std::size_t get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const;
    std::size_t get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const;

  protected:
    /**
     *  Add a vertex channel from the delegate mesh to this object's
     * vertex channels.
     *
     * @param channelName the name of the vertex channel to add.
     */
    void append_delegate_vertex_channel( const frantic::tstring& channelName );

    /**
     *  Add a face channel from the delegate mesh to this object's face
     * channels.
     *
     * @param channelName the name of the face channel to add.
     */
    void append_delegate_face_channel( const frantic::tstring& channelName );

    /**
     * Add all vertex channels allowed by the given channel_propagation_policy
     */
    void append_delegate_vertex_channels( const frantic::channels::channel_propagation_policy& cpp );

    /**
     * Add all vertex channels allowed by the given channel_propagation_policy
     */
    void append_delegate_face_channels( const frantic::channels::channel_propagation_policy& cpp );

    /**
     * Add all vertex and face channels allowed by the given channel_propagation_policy
     */
    void append_delegate_channels( const frantic::channels::channel_propagation_policy& cpp );

    const mesh_interface* m_mesh;

  private:
    delegated_const_mesh_interface(); // not implemented

    void append_delegate_channel( const frantic::geometry::mesh_channel* delegateChannel,
                                  const frantic::tstring& channelNameForErrorMessage );

    void append_delegate_channels( const frantic::geometry::mesh_interface::mesh_channel_map& channels,
                                   const frantic::channels::channel_propagation_policy& cpp );
};

} // namespace geometry
} // namespace frantic
