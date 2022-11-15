// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( E57_AVAILABLE )

#include <frantic/files/files.hpp>
#include <frantic/misc/exception_stream.hpp>
#include <frantic/particles/streams/e57_particle_mutex.hpp>
#include <frantic/particles/streams/e57_particle_ostream.hpp>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace frantic {
namespace particles {
namespace streams {

using namespace frantic;
using namespace frantic::particles;
using namespace e57;
using frantic::channels::data_type_float32;
using frantic::channels::data_type_float64;
using std::runtime_error;

namespace {
template <typename T>
void add_source_buffer_triple( ImageFile& imf, std::vector<SourceDestBuffer>& sourceBuffers, const char* xname,
                               const char* yname, const char* zname, const particle_array& particleBuffer,
                               const size_t offset, const size_t bufferSize, const bool doConversion,
                               const bool doScaling, const size_t stride ) {
    sourceBuffers.push_back( SourceDestBuffer( imf, xname, (T*)( particleBuffer[0] + offset ) + 0, bufferSize,
                                               doConversion, doScaling, stride ) );
    sourceBuffers.push_back( SourceDestBuffer( imf, yname, (T*)( particleBuffer[0] + offset ) + 1, bufferSize,
                                               doConversion, doScaling, stride ) );
    sourceBuffers.push_back( SourceDestBuffer( imf, zname, (T*)( particleBuffer[0] + offset ) + 2, bufferSize,
                                               doConversion, doScaling, stride ) );
}
} // namespace

// Private methods

void e57_particle_ostream::initialize_stream( const channels::property_map* generalMetadata,
                                              const std::map<tstring, channels::property_map>* channelMetadata ) {
    try {
        // Get the initialized root node (a Structure).
        StructureNode root = m_imf->root();

        // Register surface normals extension in case the file has normals
        m_imf->extensionsAdd( "nor", "http://www.libe57.org/E57_NOR_surface_normals.txt" );

        // Set per-file properties.
        root.set( "formatName", StringNode( *m_imf, "ASTM E57 3D Imaging Data File" ) ); // Path name: "/formatName"
        root.set( "guid", StringNode( *m_imf, "3F2504E0-4F89-11D3-9A0C-0305E82C3300" ) );
        int astmMajor, astmMinor;
        ustring libraryId;
        e57::Utilities::getVersions( astmMajor, astmMinor, libraryId );
        root.set( "versionMajor", IntegerNode( *m_imf, astmMajor ) ); // Path name: "/versionMajor"
        root.set( "versionMinor", IntegerNode( *m_imf, astmMinor ) ); // Path name: "/versionMinor"

        if( generalMetadata ) {
            // Save a WKT string identifying the coordinate reference system (CRS).
            bool hasCoordSys = generalMetadata->has_property( _T("CoordSys") );
            bool hasLengthUnitInMicro = generalMetadata->has_property( _T("LengthUnitInMicrometers") );
            std::string lengthUnitStr, axesStr;
            if( hasCoordSys ) {
                int coordSysEnum = generalMetadata->get<int>( _T("CoordSys") );
                if( 0 < coordSysEnum && coordSysEnum < 5 ) {
                    const std::string nameY[] = { "", "up", "north", "up", "south" };
                    const std::string nameZ[] = { "", "south", "up", "north", "up" };
                    axesStr = boost::str( boost::format( "  AXIS[\"(X)\",east],\n"
                                                         "  AXIS[\"(Y)\",%1%],\n"
                                                         "  AXIS[\"(Z)\",%2%]" ) %
                                          nameY[coordSysEnum] % nameZ[coordSysEnum] );
                }
                if( hasLengthUnitInMicro )
                    axesStr += ",\n";
            }
            if( hasLengthUnitInMicro ) {
                double lengthUnitInMicro = generalMetadata->get<double>( _T("LengthUnitInMicrometers") );
                lengthUnitStr =
                    boost::str( boost::format( "  LENGTHUNIT[\"metre\",%1%]" ) % ( lengthUnitInMicro * 10e-6 ) );
            }
            if( hasCoordSys || hasLengthUnitInMicro ) {
                std::string wkt = boost::str( boost::format( "CS[Cartesian,3],\n"
                                                             "%1%"
                                                             "%2%" ) %
                                              axesStr % lengthUnitStr );
                root.set( "coordinateMetadata", StringNode( *m_imf, wkt ) ); // Path name: "/coordinateMetadata"
            }
        }

        // Create 3D data area.
        VectorNode data3D = VectorNode( *m_imf, true );
        root.set( "data3D", data3D ); // Path name: "/data3D"

        // Add first scan
        StructureNode scan = StructureNode( *m_imf );
        data3D.append( scan ); // Path name: "/data3D/0"
        // Path name: "/data3D/0/guid"
        scan.set( "guid",
                  StringNode( *m_imf, boost::lexical_cast<std::string>( boost::uuids::random_generator()() ) ) );

        // Make a prototype of data types that will be stored in points record.
        StructureNode proto = StructureNode( *m_imf );
        if( m_particleBuffer.has_channel( _T("Position") ) ) {
            FloatPrecision precision =
                ( m_particleChannelMapForFile[_T("Position")].data_type() == data_type_float64 ? E57_DOUBLE
                                                                                               : E57_SINGLE );
            proto.set( "cartesianX",
                       FloatNode( *m_imf, 0.0, precision ) ); // Path name: "/data3D/0/points/0/cartesianX"
            proto.set( "cartesianY", FloatNode( *m_imf, 0.0, precision ) ); // ...
            proto.set( "cartesianZ", FloatNode( *m_imf, 0.0, precision ) );
        }
        if( m_particleBuffer.has_channel( _T("Intensity") ) ) {
            FloatPrecision precision =
                ( m_particleChannelMapForFile[_T("Intensity")].data_type() == data_type_float64 ? E57_DOUBLE
                                                                                                : E57_SINGLE );
            proto.set( "intensity", FloatNode( *m_imf, 0.0, precision ) );
        }
        if( m_particleBuffer.has_channel( _T("Color") ) ) {
            FloatPrecision precision =
                ( m_particleChannelMapForFile[_T("Color")].data_type() == data_type_float64 ? E57_DOUBLE : E57_SINGLE );
            proto.set( "colorRed", FloatNode( *m_imf, 0.0, precision ) );
            proto.set( "colorGreen", FloatNode( *m_imf, 0.0, precision ) );
            proto.set( "colorBlue", FloatNode( *m_imf, 0.0, precision ) );
        }
        if( m_particleBuffer.has_channel( _T("Normal") ) ) {
            FloatPrecision precision =
                ( m_particleChannelMapForFile[_T("Normal")].data_type() == data_type_float64 ? E57_DOUBLE
                                                                                             : E57_SINGLE );
            proto.set( "nor:normalX", FloatNode( *m_imf, 0.0, precision ) );
            proto.set( "nor:normalY", FloatNode( *m_imf, 0.0, precision ) );
            proto.set( "nor:normalZ", FloatNode( *m_imf, 0.0, precision ) );
        }

        // Make empty codecs vector for use in creating `m_points`
        // Since codecs parameter is empty, FoundationAPI assumes all fields will use the BitPack codec.
        VectorNode codecs = VectorNode( *m_imf, true );
        CompressedVectorNode points = CompressedVectorNode( *m_imf, proto, codecs );
        scan.set( "points", points ); // Path Name: "/data3D/0/points".

        // Add Cartesian bounding box to scan.
        if( channelMetadata ) {
            std::map<tstring, channels::property_map>::const_iterator p = channelMetadata->find( _T("Position") );
            if( p != channelMetadata->end() ) {
                const channels::property_map pm = p->second;
                graphics::boundbox3fd bounds;
                if( pm.get( _T("Extents"), bounds ) ) {
                    // Set the bounding box values
                    // Path names: "/data3D/0/cartesianBounds/xMinimum", etc.
                    StructureNode bbox = StructureNode( *m_imf );
                    bbox.set( "xMinimum", FloatNode( *m_imf, bounds.xminimum() ) );
                    bbox.set( "yMinimum", FloatNode( *m_imf, bounds.yminimum() ) );
                    bbox.set( "zMinimum", FloatNode( *m_imf, bounds.zminimum() ) );
                    bbox.set( "xMaximum", FloatNode( *m_imf, bounds.xmaximum() ) );
                    bbox.set( "yMaximum", FloatNode( *m_imf, bounds.ymaximum() ) );
                    bbox.set( "zMaximum", FloatNode( *m_imf, bounds.zmaximum() ) );
                    scan.set( "cartesianBounds", bbox );
                }
            }
        }

        // Set up buffers to read and write from
        const bool doConversion = true, doScaling = true;
        const size_t stride = particle_size();
        std::vector<SourceDestBuffer> sourceBuffers;

        if( m_particleBuffer.has_channel( _T("Position") ) ) {
            size_t offset = m_particleChannelMapForFile.channel_offset( _T("Position") );
            switch( m_particleChannelMapForFile[_T("Position")].data_type() ) {
            case data_type_float32:
                add_source_buffer_triple<float>( *m_imf, sourceBuffers, "cartesianX", "cartesianY", "cartesianZ",
                                                 m_particleBuffer, offset, m_ACCEPTABLE_BUFFER_SIZE, doConversion,
                                                 doScaling, stride );
                break;
            case data_type_float64:
                add_source_buffer_triple<double>( *m_imf, sourceBuffers, "cartesianX", "cartesianY", "cartesianZ",
                                                  m_particleBuffer, offset, m_ACCEPTABLE_BUFFER_SIZE, doConversion,
                                                  doScaling, stride );
                break;
            default:
                throw exception_stream() << "e57_particle_ostream.initialize_stream: Particle position expected to be "
                                            "either float32 or float64, "
                                            "got frantic::channels::data_type_t "
                                         << m_particleChannelMapForFile[_T("Position")].data_type();
            }
        }
        if( m_particleChannelMap.has_channel( _T("Intensity") ) ) {
            size_t offset = m_particleChannelMapForFile.channel_offset( _T("Intensity") );
            switch( m_particleChannelMapForFile[_T("Intensity")].data_type() ) {
            case data_type_float32:
                sourceBuffers.push_back(
                    SourceDestBuffer( *m_imf, "intensity", (float*)( m_particleBuffer[0] + offset ),
                                      m_ACCEPTABLE_BUFFER_SIZE, doConversion, doScaling, stride ) );
                break;
            case data_type_float64:
                sourceBuffers.push_back(
                    SourceDestBuffer( *m_imf, "intensity", (double*)( m_particleBuffer[0] + offset ),
                                      m_ACCEPTABLE_BUFFER_SIZE, doConversion, doScaling, stride ) );
                break;
            default:
                throw exception_stream() << "e57_particle_ostream.initialize_stream: Particle intensity expected to be "
                                            "either float32 or float64, "
                                            "got frantic::channels::data_type_t "
                                         << m_particleChannelMapForFile[_T("Intensity")].data_type();
            }
        }
        if( m_particleChannelMap.has_channel( _T("Color") ) ) {
            size_t offset = m_particleChannelMapForFile.channel_offset( _T("Color") );
            switch( m_particleChannelMapForFile[_T("Color")].data_type() ) {
            case data_type_float32:
                add_source_buffer_triple<float>( *m_imf, sourceBuffers, "colorRed", "colorGreen", "colorBlue",
                                                 m_particleBuffer, offset, m_ACCEPTABLE_BUFFER_SIZE, doConversion,
                                                 doScaling, stride );
                break;
            case data_type_float64:
                add_source_buffer_triple<double>( *m_imf, sourceBuffers, "colorRed", "colorGreen", "colorBlue",
                                                  m_particleBuffer, offset, m_ACCEPTABLE_BUFFER_SIZE, doConversion,
                                                  doScaling, stride );
                break;
            default:
                throw exception_stream() << "e57_particle_ostream.initialize_stream: Particle color expected to be "
                                            "either float32 or float64, "
                                            "got frantic::channels::data_type_t "
                                         << m_particleChannelMapForFile[_T("Color")].data_type();
            }
        }
        if( m_particleChannelMap.has_channel( _T("Normal") ) ) {
            size_t offset = m_particleChannelMapForFile.channel_offset( _T("Normal") );
            switch( m_particleChannelMapForFile[_T("Normal")].data_type() ) {
            case data_type_float32:
                add_source_buffer_triple<float>( *m_imf, sourceBuffers, "nor:normalX", "nor:normalY", "nor:normalZ",
                                                 m_particleBuffer, offset, m_ACCEPTABLE_BUFFER_SIZE, doConversion,
                                                 doScaling, stride );
                break;
            case data_type_float64:
                add_source_buffer_triple<double>( *m_imf, sourceBuffers, "nor:normalX", "nor:normalY", "nor:normalZ",
                                                  m_particleBuffer, offset, m_ACCEPTABLE_BUFFER_SIZE, doConversion,
                                                  doScaling, stride );
                break;
            default:
                throw exception_stream() << "e57_particle_ostream.initialize_stream: Particle normal expected to be "
                                            "either float32 or float64, "
                                            "got frantic::channels::data_type_t "
                                         << m_particleChannelMapForFile[_T("Normal")].data_type();
            }
        }

        // Create the writer
        m_writer.reset( new CompressedVectorWriter( points.writer( sourceBuffers ) ) );
    } catch( const E57Exception& e ) {
        throw exception_stream() << "e57_particle_ostream.initialize_stream: " << e.what() << " Error code "
                                 << e.errorCode() << " Context " << e.context().c_str() << " in "
                                 << e.sourceFunctionName() << " line " << e.sourceLineNumber();
    }
}

void e57_particle_ostream::flush() {
    if( m_particleBufferSize == 0 )
        return;

    try {
        m_writer->write( m_particleBufferSize );
    } catch( const E57Exception& e ) {
        throw exception_stream() << "e57_particle_ostream.flush: " << e.what() << " Error code " << e.errorCode()
                                 << " Context " << e.context().c_str() << " in " << e.sourceFunctionName() << " line "
                                 << e.sourceLineNumber();
    }
    m_particleBufferSize = 0;
}

// Public methods

e57_particle_ostream::e57_particle_ostream( const tstring& file, const channel_map& particleChannelMap,
                                            boost::int64_t expectedParticleCount = -1,
                                            const channels::property_map* generalMetadata,
                                            const std::map<tstring, channels::property_map>* channelMetadata )
    : m_file( file )
    , m_particleChannelMap( particleChannelMap )
    , m_currentParticleIndex( 0 )
    , m_expectedParticleCount( expectedParticleCount )
    // `*m_writer` is not properly initialized right now and does not have a default
    // constructor, but we will ensure it is set up properly in `initialize_stream()`
    , m_writer( NULL )
    , m_imf( NULL )
    , m_particleBufferSize( 0 ) {

    boost::mutex::scoped_lock lock( g_e57ImageFileConstructorMutex );
    m_imf.reset( new ImageFile( strings::to_string( file ), "w" ) );
    lock.unlock();

    m_particleChannelMapForFile = channels::channel_map();
    for( std::size_t i = 0, iEnd = particleChannelMap.channel_count(); i < iEnd; ++i ) {
        tstring name;
        channels::data_type_t dataType;
        size_t arity;
        particleChannelMap.get_channel_definition( i, name, dataType, arity );
        // E57 doesn't support float16
        if( dataType == channels::data_type_float16 )
            dataType = channels::data_type_float32;
        m_particleChannelMapForFile.define_channel( name, arity, dataType );
    }
    m_particleChannelMapForFile.end_channel_definition();
    m_pcmAdaptor.set( m_particleChannelMapForFile, m_particleChannelMap );
    m_particleBuffer = particle_array( m_particleChannelMapForFile );
    m_particleBuffer.reserve( m_ACCEPTABLE_BUFFER_SIZE );

    initialize_stream( generalMetadata, channelMetadata );
};

e57_particle_ostream::~e57_particle_ostream() {
    try {
        close();
    } catch( const exception_stream& e ) {
        FF_LOG( error ) << e.what();
    }
};

const channel_map& e57_particle_ostream::get_channel_map() const { return m_particleChannelMap; };

void e57_particle_ostream::set_channel_map( const channel_map& particleChannelMap ) {
    m_particleChannelMap = particleChannelMap;

    // Initialize the adaptor for converting the particle format to the one in the file
    m_pcmAdaptor.set( m_particleChannelMapForFile, m_particleChannelMap );
};

void e57_particle_ostream::close() {
    try {
        if( m_imf->isOpen() ) {
            flush();

            if( m_expectedParticleCount >= 0 && m_expectedParticleCount != m_currentParticleIndex ) {
                throw exception_stream()
                    << "e57_particle_ostream.close: Closed a file without writing the specified number of particles. "
                    << m_expectedParticleCount << " expected vs. " << m_currentParticleIndex << " written.";
            }

            if( m_writer.get() && m_writer->isOpen() )
                m_writer->close();

            m_imf->close();
        }
    } catch( const E57Exception& e ) {
        throw exception_stream() << "e57_particle_ostream.close: Could not close m_imf. " << e.what() << " Error code "
                                 << e.errorCode();
    }
};

std::size_t e57_particle_ostream::particle_size() const { return m_particleChannelMapForFile.structure_size(); };

void e57_particle_ostream::put_particle( const char* rawParticleData ) {
    if( m_expectedParticleCount >= 0 && m_currentParticleIndex >= m_expectedParticleCount ) {
        throw runtime_error( boost::str(
            boost::format(
                "e57_particle_ostream.put_particle: Tried to write more particles than were specified to the "
                "file \"%1%\". The specified number of particles is %2%" ) %
            strings::to_string( m_file ) % ( boost::lexical_cast<std::string>( m_expectedParticleCount ) ) ) );
    }

    if( m_particleBufferSize >= m_ACCEPTABLE_BUFFER_SIZE )
        flush();
    m_pcmAdaptor.copy_structure( m_particleBuffer[m_particleBufferSize], rawParticleData ); // add to end
    m_particleBufferSize++;

    m_currentParticleIndex++;
};

} // namespace streams
} // namespace particles
} // namespace frantic

#endif
