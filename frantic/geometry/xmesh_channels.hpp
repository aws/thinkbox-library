// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file xmesh_channels.hpp
 * @author Darcy Harrison
 * This file contains objects used by xmesh_reader and xmesh_writer for encapsulating channel information.
 */
#pragma once

#include <boost/filesystem/path.hpp>

#include <frantic/channels/named_channel_data.hpp>

namespace frantic {
namespace geometry {

/**
 * This class encapsulates the useful information about a vertex channel in an .xmesh file.
 * It provides read-only access to the data stored within, so that a user may conveniently inspect
 * the channel information.
 */
class xmesh_vertex_channel {
    // Stores the path to the vertex .xmdata file.
    boost::filesystem::path vertexPath;

    // Stores the path to the face .xmdata file, or an empty string if there are no faces.
    boost::filesystem::path facePath;

    // Stores the number of vertices expected in the vertex .xmdata file.
    std::size_t vertexCount;

    // Stores the number of faces expected in the custom faces .xmdata file
    std::size_t faceCount;

    // Stores the type and arity of each vertex in .xmdata file
    std::pair<frantic::channels::data_type_t, std::size_t> vertexType;

    friend class xmesh_reader;
    friend class xmesh_writer;

  public:
    xmesh_vertex_channel()
        : vertexCount( 0 )
        , faceCount( 0 )
        , vertexType( frantic::channels::data_type_invalid, 0 ) {}

    /**
     * @return The filesystem path to the .xmdata file for this channel.
     */
    const boost::filesystem::path& get_vertex_file_path() const { return vertexPath; }

    /**
     * @return The filesystem path to the .xmdata file for ther custom faces of this channel. Will be empty
     *          if this channel does not have custom faces.
     */
    const boost::filesystem::path& get_face_file_path() const { return facePath; }

    /**
     * @return The number of vertices in this channel.
     */
    std::size_t get_vertex_count() const { return vertexCount; }

    /**
     * @return The number of custom faces in this channel. Will be 0 if there are no custom faces.
     */
    std::size_t get_face_count() const { return faceCount; }

    /**
     * @return The data type stored in each vertex of this channel.
     */
    std::pair<frantic::channels::data_type_t, std::size_t> get_vertex_type() const { return vertexType; }

    /**
     * @return The numeric type that is stored in each vertex of this channel.
     */
    frantic::channels::data_type_t get_vertex_primitive_type() const { return vertexType.first; }

    /**
     * @return The number of primitive numeric types stored per-vertex in this channel.
     */
    std::size_t get_vertex_primitive_arity() const { return vertexType.second; }

    /**
     * @return The size of the primitives being stored in this channel, ie sizof(primitive type)*arity
     */
    std::size_t get_vertex_primitive_size() const {
        return frantic::channels::sizeof_channel_data_type( vertexType.first ) * vertexType.second;
    }
};

/**
 * This class encapsulates the usefule information about a face channel in an .xmesh file.
 * It provides read-only access to the data stored within, so that a user may conveniently inspect
 * the channel information.
 */
class xmesh_face_channel {
    // Stores the path to the .xmdata file.
    boost::filesystem::path facePath;

    // Stores the number of elements in the .xmdata file
    std::size_t faceCount;

    // Stores the type of data stored per-face in the .xmdata file
    std::pair<frantic::channels::data_type_t, std::size_t> faceType;

    friend class xmesh_reader;
    friend class xmesh_writer;

  public:
    xmesh_face_channel()
        : faceCount( 0 )
        , faceType( frantic::channels::data_type_invalid, 0 ) {}

    /**
     * @return The filesystem path to the .xmdata file for this channel.
     */
    const boost::filesystem::path& get_face_file_path() const { return facePath; }

    /**
     * @return The number of per-face elements in this channel.
     */
    std::size_t get_face_count() const { return faceCount; }

    /**
     * @return The data type stored in each per-face element of this channel.
     */
    std::pair<frantic::channels::data_type_t, std::size_t> get_face_type() const { return faceType; }

    /**
     * @return The numeric type that is stored in each per-face element of this channel.
     */
    frantic::channels::data_type_t get_face_primitive_type() const { return faceType.first; }

    /**
     * @return The number of primitive numeric types stored per-face in this channel.
     */
    std::size_t get_face_primitive_arity() const { return faceType.second; }

    /**
     * @return The size of the primitives being stored in this channel, ie sizof(primitive type)*arity
     */
    std::size_t get_face_primitive_size() const {
        return frantic::channels::sizeof_channel_data_type( faceType.first ) * faceType.second;
    }
};

} // namespace geometry
} // namespace frantic
