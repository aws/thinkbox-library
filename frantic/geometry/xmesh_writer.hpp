// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file xmesh_writer.hpp
 * @author Darcy Harrison
 * This file contains an object capable of creating .xmesh files and the .xmdata files that store
 * individual channel information.
 */
#pragma once

#include <boost/thread/mutex.hpp>

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/geometry/xmesh_channels.hpp>
#include <frantic/geometry/xmesh_metadata.hpp>

namespace frantic {
namespace geometry {

/**
 * This class provides facilities for creating new .xmesh files, and the associated .xmdata files that store
 * the .xmesh channel data.
 */
class xmesh_writer {
    typedef boost::mutex mutex_t;

    boost::filesystem::path m_rootFilename;

    mutex_t m_vertexChannelsMutex;
    std::map<frantic::tstring, xmesh_vertex_channel> m_vertexChannels;
    mutex_t m_faceChannelsMutex;
    std::map<frantic::tstring, xmesh_face_channel> m_faceChannels;

    bool m_alreadyClosed;
    int m_compressionLevel;

    frantic::geometry::xmesh_metadata m_metadata;

    typedef frantic::channels::data_type_t data_type_t;

  public:
    /**
     * Will create a new xmesh_writer, which will create an xmesh file at the specified location once
     * the xmesh_writer is destroyed or when close() is called.
     * @param path The filesystem location at which to create the new xmesh file.
     */
    xmesh_writer( const frantic::tstring& path );

    /**
     * Will destroy this xmesh_writer.
     * If close() was not called before this it will remove all the files it wrote.
     */
    virtual ~xmesh_writer();

    /**
     * This function will create a new channel with name 'name' and the specified type. The data file will be
     * immediately written to disk, and the channel's description will be stored for writing to the
     * .xmesh file later.
     * @param name The name of the new channel
     * @param pData Pointer to the data to write to disk
     * @param type The type of the data to write
     * @param arity The arity of the data to write
     * @param count The number of type[arity] elements to in pData
     */
    void write_vertex_channel( const frantic::tstring& name, const char* pData, data_type_t type, std::size_t arity,
                               std::size_t count );

    /*
     * This function will create a new channel with the given name, and point it at pre-existing
     * data files.
     * @param name The name of the new channel
     * @param vertexDataPath The path to the .xmdata of the vertex channel
     * @param faceDataPath The path to the .xmdata of the custom faces for this vertex channel
     */
    /*void add_existing_vertex_channel(
      const std::string& name,
      const std::string& vertexDataPath
    );*/

    /**
     * This function will write the custom faces for an existing vertex channel. This channel must already have been
     * created by a call to write_vertex_channel().
     * @param name The name of the channel to write the custom faces for
     * @param pData Pointer to the face data
     * @param count The number of face indices to write
     */
    void write_vertex_channel_faces( const frantic::tstring& name, const int* pData, std::size_t count );

    /**
     * This function will create a new per-face channel with name 'name' and the specified type. The data file will be
     * immediately written to disk, and the channel's description will be stored for writing to the
     * .xmesh file later.
     * @param name The name of the new face channel
     * @param pData Pointer to the data to write to disk
     * @param type The type of the data to write
     * @param arity The arity of the data to write
     * @param count The number of type[arity] elements to in pData
     */
    void write_face_channel( const frantic::tstring& name, const char* pData, data_type_t type, std::size_t arity,
                             std::size_t count );

    /**
     *  Set the zlib compression level for channel data files.
     *
     * @param compressionLevel the zlib compression level to use for
     *        compressing channel data files.
     */
    void set_compression_level( int compressionLevel );

    /**
     *  Set the metadata to include in the output xmesh file.
     *
     * @param metadata the metadata to include in the output xmesh file.
     */
    void set_metadata( const frantic::geometry::xmesh_metadata& metadata );

    /**
     * This function will end the creation of the .xmesh file and write the master .xmesh file to disk with the
     * information about the channels created through previous calls to write_channel, and write_channel_faces.
     * The .xmesh file will be written to the path as specified in the constructor call.
     */
    void close();
};

/**
 *  Write a single array of data.
 *  This is intended for internal use by xmesh_writer and xmesh_sequence_saver.
 *
 *  Format:
 *    "xmeshdat"     8 bytes
 *    "<datatag>"    12 bytes (based on the named channel data type, so all 12 bytes may not actually be used)
 *    count          4 bytes
 *    dataSize       4 bytes
 *    <array data>   (count * dataSize) bytes
 *
 * @todo put this somewhere more appropriate.
 */
void write_xmesh_array_file( const boost::filesystem::path& path, const void* data, const std::string& dataType,
                             std::size_t dataSize, std::size_t numData, int zlibCompressionLevel );

} // namespace geometry
} // namespace frantic
