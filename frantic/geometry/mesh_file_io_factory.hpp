// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/filesystem/path.hpp>

#include <frantic/geometry/mesh_interface.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/logging/progress_logger.hpp>

#ifdef OIIO_LIB_AVAILABLE
#include <OpenImageIO/imagebuf.h>
#endif

namespace frantic {
namespace geometry {

class mesh_file_io_factory {
  public:
    mesh_file_io_factory();

    /**
     * @return the coordinate system used for meshes returned by and passed into this class.
     */
    frantic::graphics::coordinate_system::option get_coordinate_system();

    /**
     * @param coordinateSystem the coordinate system to use for the meshes returned by and
     *    passed into this class.
     */
    void set_coordinate_system( frantic::graphics::coordinate_system::option coordinateSystem );

    /**
     * @return the coordinate system used for mesh files that don't have a specified coordinate
     *    system.
     */
    frantic::graphics::coordinate_system::option get_default_file_coordinate_system();

    /**
     * @param  the coordinate system to use for mesh files that don't have a specified coordinate
     *    system.
     */
    void set_default_file_coordinate_system( frantic::graphics::coordinate_system::option defaultFileCoordinateSystem );

    /**
     * @return the unit length used for meshes returned by and passed into this class, measured in meters.
     *    Default is 0. If 0 the default unit length or "no unit" will be used instead
     */
    double get_unit_length_in_meters();

    /**
     * @param unitLengthInMeters the unit length to use for meshes returned by and passed into this class, measured in
     * meters. Default is 0. If 0 the default unit length or "no unit" will be used instead
     */
    void set_unit_length_in_meters( double unitLengthInMeters );

    /**
     * @return the number of threads to use for saving channel data files.
     */
    std::size_t get_thread_count() const;

    /**
     * @param threadCount the number of threads to use for saving channel
     *    data files.
     */
    void set_thread_count( std::size_t threadCount );

    /**
     * @return true if only XMesh standard channels will be written to
     *    .xmesh files.
     */
    bool get_use_xmesh_standard() const;

    /**
     * @param useXMeshStandard true to write only standard XMesh channels
     *    to .xmesh files, or false to write all channels.
     */
    void set_use_xmesh_standard( bool useXMeshStandard );

    /**
     *  Read a polymesh from the specified path.
     *
     * @param path the path to read a mesh from.
     * @return a new polymesh containing data read from the path.
     */
    frantic::geometry::polymesh3_ptr read_polymesh3( const boost::filesystem::path& path );

    /**
     *  Write a polymesh to the specified file path.
     *
     * @param path the path to write a mesh to.
     * @param mesh the mesh to write.
     */
    void write( const boost::filesystem::path& path, frantic::geometry::polymesh3_ptr mesh );

    /**
     *  Write a polymesh to the specified file path.
     *
     * @param path the path to write a mesh to.
     * @param mesh the mesh to write.
     * @param progress The progress logger to update while writing the file.
     */
    void write( const boost::filesystem::path& path, frantic::geometry::mesh_interface_ptr mesh,
                frantic::logging::progress_logger& progress );

#ifdef OIIO_LIB_AVAILABLE
    /**
     * Sets the texture data that will be applied to the mesh
     *
     * @param texture The texture to be applied to the mesh
     */
    void set_texture( const OIIO::ImageBuf& texture );
#endif

  private:
    mesh_file_io_factory( const mesh_file_io_factory& );            // not implemented
    mesh_file_io_factory& operator=( const mesh_file_io_factory& ); // not implemented

    frantic::graphics::coordinate_system::option m_coordinateSystem;
    frantic::graphics::coordinate_system::option m_defaultFileCoordinateSystem;
    double m_unitLengthInMeters;
    double m_defaultFileUnitLengthInMeters;
    std::size_t m_threadCount;
    bool m_useXMeshStandard;

#ifdef OIIO_LIB_AVAILABLE
    boost::shared_ptr<OIIO::ImageBuf> m_texture;
#endif
};

} // namespace geometry
} // namespace frantic
