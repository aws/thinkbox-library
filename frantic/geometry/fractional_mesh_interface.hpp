// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <frantic/geometry/dcel_mesh_interface.hpp>
#include <frantic/geometry/mesh_interface.hpp>

namespace frantic {
namespace geometry {

class fractional_index_map;
class vertex_index_map;

/**
 * Abstract class to implement any mesh_interface functions common to fractional_mesh_interfaces.
 */
class fractional_mesh_interface_base : public mesh_interface {
  protected:
    std::unique_ptr<dcel_mesh_interface> m_adjacencyDelegate;
    mesh_interface* m_mesh;

  public:
    virtual std::size_t get_num_elements() const;

    virtual std::size_t get_face_element_index( std::size_t faceIndex ) const;

    virtual bool has_adjacency() const;

    virtual bool is_valid() const;

    virtual bool request_channel( const frantic::tstring& channelName, bool vertexChannel, bool forOutput,
                                  bool throwOnError = true );

    virtual void init_adjacency();

    virtual bool init_vertex_iterator( frantic::geometry::vertex_iterator& vIt, std::size_t vertexIndex ) const;

    virtual bool advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const;

    virtual std::size_t get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const;

    virtual std::size_t get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const;

    virtual std::size_t get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const;

    virtual bool is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const;

    virtual bool is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const;

    virtual void init_face_iterator( frantic::geometry::face_iterator& fIt, std::size_t faceIndex ) const;

    virtual bool advance_face_iterator( frantic::geometry::face_iterator& fIt ) const;

    virtual std::size_t get_face_neighbor( frantic::geometry::face_iterator& fIt ) const;

    virtual std::size_t get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const;

    virtual std::size_t get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const;

    virtual ~fractional_mesh_interface_base() {}
};

/**
 * This class accesses a fraction of the faces (and their corresponding vertices) from a mesh_interface.
 * \warning It is meant as a read only accessor and any changes to the underlying mesh_interface should be followed by
 * calling mesh_reset() to avoid undefined behavior.
 */
class fractional_face_interface : public fractional_mesh_interface_base {
  private:
    boost::shared_ptr<fractional_index_map> m_faceIndexMap;
    boost::shared_ptr<vertex_index_map> m_vertexIndexMap;
    double m_fraction;
    boost::int64_t m_limit;

  public:
    fractional_face_interface( mesh_interface* mesh, double fraction,
                               boost::int64_t limit = std::numeric_limits<boost::int64_t>::max() );

    void set_mesh( mesh_interface* mesh, double fraction, boost::int64_t limit );

    void reset_mesh();

    virtual std::size_t get_num_verts() const;

    virtual void get_vert( std::size_t index, float ( &outValues )[3] ) const;

    virtual std::size_t get_num_faces() const;

    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const;

    virtual std::size_t get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const;

    virtual void get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const;

    virtual void get_face_verts( std::size_t faceIndex, float outValues[][3] ) const;

    virtual ~fractional_face_interface() {}
};

/**
 * This class accesses a fraction of the vertices from a mesh_interface and will exclude any faces or face channels.
 * \warning It is meant as a read only accessor and any changes to the underlying mesh_interface should be followed by
 * calling mesh_reset() to avoid undefined  behavior.
 */
class fractional_vertex_interface : public fractional_mesh_interface_base {
  private:
    boost::shared_ptr<fractional_index_map> m_vertexIndexMap;
    double m_fraction;
    boost::int64_t m_limit;

  public:
    fractional_vertex_interface( mesh_interface* mesh, double fraction,
                                 boost::int64_t limit = std::numeric_limits<boost::int64_t>::max() );

    void set_mesh( mesh_interface* mesh, double fraction, boost::int64_t limit );

    void reset_mesh();

    virtual bool request_channel( const frantic::tstring& channelName, bool vertexChannel, bool forOutput,
                                  bool throwOnError = true );

    virtual std::size_t get_num_verts() const;

    virtual void get_vert( std::size_t index, float ( &outValues )[3] ) const;

    virtual std::size_t get_num_faces() const;

    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const;

    virtual std::size_t get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const;

    virtual void get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const;

    virtual void get_face_verts( std::size_t faceIndex, float outValues[][3] ) const;

    virtual ~fractional_vertex_interface() {}
};
} // namespace geometry
} // namespace frantic
