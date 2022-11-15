// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/mesh_file_io_factory.hpp>

#include <frantic/files/paths.hpp>
#include <frantic/geometry/mesh_interface_file_io.hpp>
#include <frantic/geometry/polymesh3_file_io.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/geometry/transformed_mesh_interface.hpp>
#include <frantic/geometry/u3d_writer.hpp>
#include <frantic/geometry/xmesh_sequence_saver.hpp>
#include <frantic/geometry/xmesh_standard_mesh_interface.hpp>

using frantic::graphics::transform4f;

namespace {

frantic::geometry::polymesh3_ptr clone( frantic::geometry::polymesh3_ptr mesh ) {
    // TODO: less roundabout method to clone a polymesh3 ?
    std::unique_ptr<frantic::geometry::mesh_interface> meshInterface(
        frantic::geometry::polymesh3_interface::create_instance( mesh ) );
    return frantic::geometry::create_polymesh3( meshInterface.get() );
}

} // anonymous namespace

namespace frantic {
namespace geometry {

mesh_file_io_factory::mesh_file_io_factory()
    : m_coordinateSystem( frantic::graphics::coordinate_system::unspecified )
    , m_defaultFileCoordinateSystem( frantic::graphics::coordinate_system::unspecified )
    , m_unitLengthInMeters( 0 )
    , m_defaultFileUnitLengthInMeters( 0 )
    , m_threadCount( 0 )
    , m_useXMeshStandard( true ) {}

frantic::graphics::coordinate_system::option mesh_file_io_factory::get_coordinate_system() {
    return m_coordinateSystem;
}

void mesh_file_io_factory::set_coordinate_system( frantic::graphics::coordinate_system::option coordinateSystem ) {
    m_coordinateSystem = coordinateSystem;
}

frantic::graphics::coordinate_system::option mesh_file_io_factory::get_default_file_coordinate_system() {
    return m_defaultFileCoordinateSystem;
}

void mesh_file_io_factory::set_default_file_coordinate_system(
    frantic::graphics::coordinate_system::option defaultFileCoordinateSystem ) {
    m_defaultFileCoordinateSystem = defaultFileCoordinateSystem;
}

double mesh_file_io_factory::get_unit_length_in_meters() { return m_unitLengthInMeters; }

void mesh_file_io_factory::set_unit_length_in_meters( double unitLength ) { m_unitLengthInMeters = unitLength; }

std::size_t mesh_file_io_factory::get_thread_count() const { return m_threadCount; }

void mesh_file_io_factory::set_thread_count( std::size_t threadCount ) { m_threadCount = threadCount; }

bool mesh_file_io_factory::get_use_xmesh_standard() const { return m_useXMeshStandard; }

void mesh_file_io_factory::set_use_xmesh_standard( bool useXMeshStandard ) { m_useXMeshStandard = useXMeshStandard; }

#ifdef OIIO_LIB_AVAILABLE
void mesh_file_io_factory::set_texture( const OIIO::ImageBuf& texture ) {
    m_texture.reset( new OIIO::ImageBuf( texture ) );
}
#endif

frantic::geometry::polymesh3_ptr mesh_file_io_factory::read_polymesh3( const boost::filesystem::path& path ) {
    const frantic::tstring type = frantic::strings::to_lower( frantic::files::to_tstring( path.extension() ) );

    frantic::geometry::polymesh3_ptr mesh;
    frantic::graphics::coordinate_system::option fileCoordinateSystem = m_defaultFileCoordinateSystem;
    double fileLengthUnitInMeters = m_defaultFileUnitLengthInMeters;

    if( type == _T(".obj") ) {
        mesh = load_obj_polymesh_file( frantic::files::to_tstring( path ) );
    } else if( type == _T(".xmesh") ) {
        mesh = load_xmesh_polymesh_file( frantic::files::to_tstring( path ) );
        fileCoordinateSystem = frantic::graphics::coordinate_system::right_handed_zup;
        // get the length unit
        frantic::geometry::xmesh_metadata fileMetadata;
        frantic::geometry::read_xmesh_metadata( path, fileMetadata );
        if( fileMetadata.has_length_unit() ) {
            fileLengthUnitInMeters = frantic::geometry::get_meters_from_xmesh_length_unit(
                fileMetadata.get_length_unit_scale(), fileMetadata.get_length_unit() );
        }
    } else if( type == _T(".stl") ) {
        mesh = load_stl_polymesh_file( frantic::files::to_tstring( path ) );
    } else if( type == _T(".ply") ) {
        mesh = load_ply_polymesh_file( frantic::files::to_tstring( path ) );
    } else {
        throw std::runtime_error(
            "mesh_file_io_factory::read_polymesh3 Error: Didn't recognize the file format of the input mesh file \"" +
            path.string() + "\"" );
    }

    if( !mesh ) {
        throw std::runtime_error( "mesh_file_io_factory::read_polymesh3 Error: mesh is NULL" );
    }

    frantic::graphics::transform4f xform;
    if( frantic::graphics::coordinate_system::create_transform( xform, fileCoordinateSystem, m_coordinateSystem ) ) {
        frantic::geometry::transform( mesh, xform );
    }

    if( fileLengthUnitInMeters != 0 && m_unitLengthInMeters != 0 ) {
        frantic::geometry::scale( mesh, static_cast<float>( fileLengthUnitInMeters / m_unitLengthInMeters ) );
    }

    return mesh;
}

void mesh_file_io_factory::write( const boost::filesystem::path& path, frantic::geometry::polymesh3_ptr mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "mesh_file_io_factory::write Error: mesh is NULL" );
    }

    const frantic::tstring type = frantic::strings::to_lower( frantic::files::to_tstring( path.extension() ) );

    frantic::graphics::coordinate_system::option fileCoordinateSystem = m_defaultFileCoordinateSystem;
    double fileLengthUnitInMeters = 0;

    if( type == _T(".xmesh") ) {
        fileCoordinateSystem = frantic::graphics::coordinate_system::right_handed_zup;
    } else if( type == _T(".fbx") ) {
        fileCoordinateSystem = m_coordinateSystem;
        fileLengthUnitInMeters = m_unitLengthInMeters;
    }

    transform4f xform;
    frantic::graphics::coordinate_system::create_transform( xform, m_coordinateSystem, fileCoordinateSystem );
    if( fileLengthUnitInMeters != 0 && m_unitLengthInMeters != 0 ) {
        xform = xform * transform4f::from_scale( float( m_unitLengthInMeters / fileLengthUnitInMeters ) );
    }
    if( !xform.is_identity() ) {
        // TODO: I'm cloning because I think it would be surprising if this function modified the
        // input mesh.  Should we do something else instead ?
        mesh = clone( mesh );
        frantic::geometry::transform( mesh, xform );
    }

    if( type == _T(".obj") ) {
        write_obj_polymesh_file( frantic::files::to_tstring( path ), mesh );
    } else if( type == _T(".xmesh") ) {
        frantic::logging::null_progress_logger progress;

        frantic::geometry::xmesh_metadata metaData;
        if( m_unitLengthInMeters != 0 ) {
            std::pair<double, frantic::geometry::xmesh_metadata::length_unit_t> lengthUnit =
                frantic::geometry::get_xmesh_length_unit_from_meters( m_unitLengthInMeters );
            metaData.set_length_unit( lengthUnit.first, lengthUnit.second );
        }
        mesh_interface_ptr meshInterface( polymesh3_interface::create_instance( mesh ).release() );
        mesh_interface_ptr outputMesh = meshInterface;
        if( m_useXMeshStandard ) {
            std::unique_ptr<mesh_interface> xmeshStandardMesh =
                create_xmesh_standard_mesh_interface( meshInterface.get() );
            outputMesh.reset( xmeshStandardMesh.release() );
        }
        if( m_threadCount > 1 ) {
            write_xmesh_mesh_file( frantic::files::to_tstring( path ), outputMesh, metaData, progress, m_threadCount );
        } else {
            write_xmesh_mesh_file( frantic::files::to_tstring( path ), outputMesh, metaData, progress );
        }
    } else if( type == _T(".stl") ) {
        // write_ascii_stl_polymesh_file( frantic::files::to_tstring( path ), mesh );

        // Using false (Magics, not SolidView) convention for the Color
        // channel, because the color doesn't seem to appear in
        // Autodesk Inventor otherwise.
        const bool isSolidView = false;
        write_binary_stl_polymesh_file( frantic::files::to_tstring( path ), mesh, isSolidView );
    } else if( type == _T(".ply") ) {
        write_ply_polymesh_file( frantic::files::to_tstring( path ), mesh );
    } else if( type == _T(".u3d") ) {
        u3d_writer writer( frantic::files::to_tstring( path ) );
        mesh_interface_ptr meshPtr( polymesh3_interface::create_instance( mesh ).release() );
#ifdef OIIO_LIB_AVAILABLE
        if( !m_texture.get() ) {
            writer.write_mesh( meshPtr.get() );
        } else {
            writer.write_mesh_with_texture( meshPtr.get(), m_texture.get() );
        }
#else
        writer.write_mesh( meshPtr.get() );
#endif
    } else {
        throw std::runtime_error(
            "mesh_file_io_factory::write Error: Didn't recognize the file format of the output mesh file \"" +
            path.string() + "\"" );
    }
}

void mesh_file_io_factory::write( const boost::filesystem::path& path, frantic::geometry::mesh_interface_ptr mesh,
                                  frantic::logging::progress_logger& progress ) {
    if( !mesh ) {
        throw std::runtime_error( "mesh_file_io_factory::write Error: mesh is NULL" );
    }

    const frantic::tstring type = frantic::strings::to_lower( frantic::files::to_tstring( path.extension() ) );

    frantic::graphics::coordinate_system::option fileCoordinateSystem = m_defaultFileCoordinateSystem;
    double fileLengthUnitInMeters = 0;

    if( type == _T(".xmesh") ) {
        fileCoordinateSystem = frantic::graphics::coordinate_system::right_handed_zup;
    } else if( type == _T(".fbx") ) {
        fileCoordinateSystem = m_coordinateSystem;
        fileLengthUnitInMeters = m_unitLengthInMeters;
    }

    transform4f xform;
    frantic::graphics::coordinate_system::create_transform( xform, m_coordinateSystem, fileCoordinateSystem );
    if( fileLengthUnitInMeters != 0 && m_unitLengthInMeters != 0 ) {
        xform = xform * transform4f::from_scale( float( m_unitLengthInMeters / fileLengthUnitInMeters ) );
    }
    if( !xform.is_identity() ) {
        mesh = transformed_mesh_interface::create_instance( mesh, xform );
    }

    if( type == _T(".obj") ) {
        write_obj_mesh_file( frantic::files::to_tstring( path ), mesh, progress );
    } else if( type == _T(".xmesh") ) {
        frantic::geometry::xmesh_metadata metaData;
        if( m_unitLengthInMeters != 0 ) {
            std::pair<double, frantic::geometry::xmesh_metadata::length_unit_t> lengthUnit =
                frantic::geometry::get_xmesh_length_unit_from_meters( m_unitLengthInMeters );
            metaData.set_length_unit( lengthUnit.first, lengthUnit.second );
        }
        frantic::geometry::mesh_interface_ptr outputMesh = mesh;
        if( m_useXMeshStandard ) {
            std::unique_ptr<mesh_interface> xmeshStandardMesh = create_xmesh_standard_mesh_interface( mesh.get() );
            outputMesh.reset( xmeshStandardMesh.release() );
        }
        if( m_threadCount > 1 ) {
            write_xmesh_mesh_file( frantic::files::to_tstring( path ), outputMesh, metaData, progress, m_threadCount );
        } else {
            write_xmesh_mesh_file( frantic::files::to_tstring( path ), outputMesh, metaData, progress );
        }
    } else if( type == _T(".stl") ) {
        // write_ascii_stl_polymesh_file( frantic::files::to_tstring( path ), mesh );

        // Using false (Magics, not SolidView) convention for the Color
        // channel, because the color doesn't seem to appear in
        // Autodesk Inventor otherwise.
        const bool isSolidView = false;
        write_binary_stl_mesh_file( frantic::files::to_tstring( path ), mesh, isSolidView, progress );
    } else if( type == _T(".ply") ) {
        write_ply_mesh_file( frantic::files::to_tstring( path ), mesh, progress );
    } else if( type == _T(".u3d") ) {
        u3d_writer writer( frantic::files::to_tstring( path ) );
#ifdef OIIO_LIB_AVAILABLE
        if( !m_texture.get() ) {
            writer.write_mesh( mesh.get() );
        } else {
            writer.write_mesh_with_texture( mesh.get(), m_texture.get() );
        }
#else
        writer.write_mesh( mesh.get() );
#endif
    } else {
        throw std::runtime_error(
            "mesh_file_io_factory::write Error: Didn't recognize the file format of the output mesh file \"" +
            path.string() + "\"" );
    }
}
} // namespace geometry
} // namespace frantic
