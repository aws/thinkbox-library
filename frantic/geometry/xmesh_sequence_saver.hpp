// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>
#include <boost/filesystem/path.hpp>
#include <frantic/channels/named_channel_data.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/xmesh_metadata.hpp>

namespace frantic {
namespace geometry {

/**
 * This function loads an .xmesh file from disk.
 *
 * @param  srcFile  The filename of the input .xmesh file to load.
 * @param  outMesh  The trimesh3 into which the .xmesh is loaded.
 */
void load_xmesh_file( const frantic::tstring& srcFile, frantic::geometry::trimesh3& outMesh );

/**
 * This is the class used to save xmesh sequences.  Usage requires simply the instantiation of
 * an xmesh_sequence_saver object and then repeated calls to the member function write_xmesh.
 */
class xmesh_sequence_saver {
  public:
    // These are some utility structs used with the xmesh sequence saver
    struct xmesh_cache_channel {
        frantic::channels::data_type_t m_dataType;
        size_t m_arity;
        frantic::graphics::raw_byte_buffer m_data;

        xmesh_cache_channel()
            : m_dataType( frantic::channels::data_type_invalid )
            , m_arity( 0 ) {}
    };

    struct xmesh_cache_vertex_channel : xmesh_cache_channel {
        std::vector<boost::int32_t> m_customFaces;
        std::vector<boost::int32_t> m_serializationCustomFaces;
    };

    struct xmesh_cache_face_channel : xmesh_cache_channel {};

  private:
    typedef std::map<frantic::tstring, xmesh_cache_vertex_channel> vertex_channels_cache_t;
    typedef std::map<frantic::tstring, xmesh_cache_face_channel> face_channels_cache_t;

    // This used to operate with separate member variables for each component in a mesh.
    // I'm just going to use a single trimesh object, which I can poke new information
    // into as required and save out the various parts of it when necessary.  I'll still
    // require the file name variables though.

    xmesh_cache_vertex_channel m_cachedGeometryChannel;
    bool m_isCachedTriangleMesh;
    std::vector<int> m_cachedFaceOffsets;

    vertex_channels_cache_t m_cachedVertexChannels;
    face_channels_cache_t m_cachedFaceChannels;

    boost::filesystem::path::string_type m_savedVertsFile;
    boost::filesystem::path::string_type m_savedFacesFile;

    int m_compressionLevel;
    std::size_t m_threadCount;

    // This is just a map from channel names to filenames, also used to indicate when a vertex
    // or face data has changed by an empty filename
    std::map<frantic::tstring, boost::filesystem::path::string_type> m_vertChannelVertDataChanged;
    std::map<frantic::tstring, boost::filesystem::path::string_type> m_vertChannelFaceDataChanged;
    std::map<frantic::tstring, boost::filesystem::path::string_type> m_faceChannelDataChanged;

    void retrieve_and_compare_mesh( const frantic::geometry::trimesh3& mesh );
    void retrieve_and_compare_vertex_channel( const frantic::geometry::trimesh3& mesh,
                                              const frantic::tstring& channelName );
    void retrieve_and_compare_face_channel( const frantic::geometry::trimesh3& mesh,
                                            const frantic::tstring& channelName );

    void retrieve_and_compare_mesh( const frantic::geometry::polymesh3_ptr mesh );
    void retrieve_and_compare_vertex_channel( const frantic::geometry::polymesh3_ptr mesh,
                                              const frantic::tstring& channelName, bool invalidateFaces = false );
    void retrieve_and_compare_face_channel( const frantic::geometry::polymesh3_ptr mesh,
                                            const frantic::tstring& channelName, bool invalidateFaces = false );

    void write_new_channels( const boost::filesystem::path& filename );
    void write_xmesh_xml_file( const boost::filesystem::path& filename,
                               const frantic::geometry::xmesh_metadata* metadata );

    void write_xmesh( const frantic::geometry::trimesh3& mesh, const frantic::geometry::xmesh_metadata* metadata,
                      const frantic::tstring& filename );
    void write_xmesh( const frantic::geometry::polymesh3_ptr mesh, const frantic::geometry::xmesh_metadata* metadata,
                      const frantic::tstring& filename );

  public:
    xmesh_sequence_saver();
    ~xmesh_sequence_saver();

    void clear();
    void write_xmesh( const frantic::geometry::trimesh3& mesh, const frantic::tstring& filename );
    void write_xmesh( const frantic::geometry::trimesh3& mesh, const frantic::geometry::xmesh_metadata& metadata,
                      const frantic::tstring& filename );

    void write_xmesh( const frantic::geometry::polymesh3_ptr polymesh, const frantic::tstring& filename );
    void write_xmesh( const frantic::geometry::polymesh3_ptr polymesh,
                      const frantic::geometry::xmesh_metadata& metadata, const frantic::tstring& filename );

    /**
     *  Set the zlib compression level for channel data files.
     *
     * @param compressionLevel the zlib compression level to use for
     *        compressing channel data files.
     */
    void set_compression_level( int compressionLevel );

    /**
     *  Set the number of threads to use for saving channel data files.
     *
     * @param threadCount the number of threads to use for saving channel
     *        data files.
     */
    void set_thread_count( std::size_t threadCount );
};

} // namespace geometry
}; // namespace frantic
