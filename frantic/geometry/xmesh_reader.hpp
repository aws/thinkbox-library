// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file xmesh_reader.hpp
 * @author Darcy Harrison
 * This file contains an object capable of parsing .xmesh files and providing access to the data contained within.
 */
#pragma once

#include <boost/filesystem/path.hpp>

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/geometry/xmesh_channels.hpp>
#include <frantic/geometry/xmesh_metadata.hpp>

namespace frantic {
namespace geometry {

/**
 * This class provides facilities for parsing an .xmesh file on disk. It allows enumeration of available channels,
 * access to channel properties, and provides a mechanism for loading the data associated with a channel.
 */
class xmesh_reader {
    // Stores the path of the .xmesh file we are parsing.
    boost::filesystem::path m_path;

    // Stores the path to the directory containing the .xmesh file we are parsing.
    boost::filesystem::path m_rootPath;

    // Stores all the vertex channels in this .xmesh, indexed by name. The 'verts' channel is a reserved
    // channel guaranteed to exist in all .xmesh files.
    std::map<frantic::tstring, xmesh_vertex_channel> m_vertexChannels;

    // Stores all the face channels in this .xmesh, indexed by name.
    std::map<frantic::tstring, xmesh_face_channel> m_faceChannels;

    // Records the version string found in the .xmesh file (or 0 if no version was found).
    int m_version;

    xmesh_metadata m_metadata;

  public:
    typedef frantic::channels::data_type_t data_type_t;
    typedef std::pair<data_type_t, std::size_t> channel_type_t;

  public:
    /**
     * This function will construct a new xmesh_reader, initializing it with the file at location 'path'.
     * It will throw an exception if the .xmesh file at 'path' is not a valid .xmesh file.
     * @param path The filesystem path to the .xmesh file to read.
     */
    xmesh_reader( const boost::filesystem::path& path );

    virtual ~xmesh_reader() {}

    /**
     *	This function returns the version number of the xmesh file.
     *	@return int version number
     */
    int get_version_number() const { return m_version; }

    /**
     * @return the metadata stored in the xmesh file.
     */
    const xmesh_metadata& get_metadata() const { return m_metadata; }

    /**
     * @return the path of the xmesh file we are parsing.
     */
    const boost::filesystem::path& get_path() const { return m_path; }
    /**
     * This will fill a vector with the names of all vertex channels. It will not include the special
     * vertex channel 'verts', which is guaranteed to exist in all .xmesh files.
     * @param outNames The vector to store the channel names in.
     */
    void get_vertex_channel_names( std::vector<frantic::tstring>& outNames ) const;

    /**
     * This will fill a vector with the names of all face channels.
     * @param outNames The vector to store the channel names in.
     */
    void get_face_channel_names( std::vector<frantic::tstring>& outNames ) const;

    /**
     * @param name The name of the channel
     * @return true iff a vertex channel with name 'name' exists, otherwise false.
     */
    bool has_vertex_channel( const frantic::tstring& name ) const;

    /**
     * @return true iff a vertex channel with name 'vertexChannelName' has custom faces, otherwise false.
     */
    bool has_custom_faces( const frantic::tstring& vertexChannelName ) const;

    /**
     * @return true iff a face channel with name 'name' exists, otherwise false.
     */
    bool has_face_channel( const frantic::tstring& name ) const;

    /**
     * @return The descriptor object for the vertex channel with name 'vertexChannelName'
     */
    const xmesh_vertex_channel& get_vertex_channel( const frantic::tstring& vertexChannelName ) const;

    /**
     * @return The descriptor object for the face channel with name 'faceChannelName'
     */
    const xmesh_face_channel& get_face_channel( const frantic::tstring& faceChannelName ) const;

    /**
     * This function will load the vertex data for the given channel into a user-supplied buffer.
     * The function will make sure that the channel data actually matches the type, arity, and count expected.
     *
     * @note It is the user's responsibility to ensure that the buffer 'pData' is of the correct size to receive
     *        'count' elements of size 'expectedArity' * sizeof('expectedType').
     * @param name The name of the vertex channel
     * @param pData A pointer to a data buffer to store the vertex channel data in.
     * @param expectedType The primitive type expected per-vertex
     * @param expectedArity The number of primitive type elements expected per-vertex.
     * @param expectedCount The number of vertices expected in this channel.
     */
    void load_vertex_channel( const frantic::tstring& name, char* pData, data_type_t expectedType,
                              std::size_t expectedArity, std::size_t expectedCount ) const;

    /**
     * This function will load the custom face indices for the given vertex channel into a user-supplied buffer.
     * The function will make sure that the channel data is int32[1], with 'expectedCount' elements. The
     * data is returned encoded as runs of indices where each polygon is terminated with a negative index.
     * The negative indices should be interpreted as (-value - 1) where 'value' is the value stored in the array.
     *
     * @note It is the user's responsibility to ensure that the buffer 'pData' is of the correct size to receive
     *        'count' elements of size: sizeof(int32).
     * @param name The name of the vertex channel
     * @param pData A pointer to a data buffer to store the custom faces.
     * @param expectedCount The number of custom face indices expected in this channel.
     */
    void load_vertex_channel_faces( const frantic::tstring& name, char* pData, std::size_t expectedCount ) const;

    /**
     * This function will load the per-face data for the given face channel into a user-supplied buffer.
     * The function will make sure that the channel data actually matches the type, arity, and count expected.
     *
     * @note It is the user's responsibility to ensure that the buffer 'pData' is of the correct size to receive
     *        'count' elements of size 'expectedArity' * sizeof('expectedType').
     * @param name The name of the face channel to load
     * @param pData A pointer to a data buffer to store the per-face channel data in.
     * @param expectedType The primitive type expected per-face
     * @param expectedArity The number of primitive type elements expected per-face.
     * @param expectedCount The number of per-face elements expected in this channel.
     */
    void load_face_channel( const frantic::tstring& name, char* pData, data_type_t expectedType,
                            std::size_t expectedArity, std::size_t expectedCount ) const;
};

/**
 *  Loads a single array of data.
 *  This is intended for internal use by xmesh_reader and load_xmesh_file.
 *
 *  Format:
 *   "xmeshdat"     8 bytes
 *   "<datatag>"    12 bytes (based on the named channel data type, so all 12 bytes may not actually be used)
 *   count          4 bytes
 *   dataSize       4 bytes
 *   <array data>   (count * dataSize) bytes
 *
 * @todo put this somewhere more appropriate
 */
void load_xmesh_array_file( const boost::filesystem::path& path, char* data, const std::string& correctDataType,
                            size_t correctDataSize, size_t correctDataCount );

} // namespace geometry
} // namespace frantic
