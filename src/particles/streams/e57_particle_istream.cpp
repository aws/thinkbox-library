// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( E57_AVAILABLE )

#include <E57Format/E57Export.h>

#include <frantic/particles/streams/e57_particle_istream.hpp>
#include <frantic/particles/streams/e57_particle_mutex.hpp>

using namespace frantic;
using namespace frantic::particles;
using boost::int64_t;
using frantic::particles::streams::e57_particle_istream;

static graphics::transform4fd get_transform_from_node( const e57::StructureNode& poseNode ) {
    if( !poseNode.isDefined( "rotation" ) && !poseNode.isDefined( "translation" ) ) {
        FF_LOG( warning ) << _T("E57 pose:\n scan pose is missing rotation and translation.") << std::endl;
    }
    frantic::graphics::quat4fd quat;
    frantic::graphics::vector3fd translationVector;
    if( poseNode.isDefined( "rotation" ) ) {
        e57::StructureNode rotation_node( poseNode.get( "rotation" ) );

        const double w = e57::FloatNode( rotation_node.get( "w" ) ).value();
        const double x = e57::FloatNode( rotation_node.get( "x" ) ).value();
        const double y = e57::FloatNode( rotation_node.get( "y" ) ).value();
        const double z = e57::FloatNode( rotation_node.get( "z" ) ).value();

        quat = frantic::graphics::quat4fd( w, x, y, z );
    }
    if( poseNode.isDefined( "translation" ) ) {
        e57::StructureNode translation_node( poseNode.get( "translation" ) );

        const double x = e57::FloatNode( translation_node.get( "x" ) ).value();
        const double y = e57::FloatNode( translation_node.get( "y" ) ).value();
        const double z = e57::FloatNode( translation_node.get( "z" ) ).value();

        translationVector.set( x, y, z );
    }

    graphics::transform4fd transform;
    quat.as_transform4f( transform );         // set the rotation quaternion to transform4f
    transform.translate( translationVector ); // combine the transform matrix with translation_vector
    return transform;
}

/**
 * @param bbox A "cartesianBounds" node.
 * @param transform The scan transform.
 * @param[out] hasBounds Does the file have bounds?
 * @return The boundbox transformed by transform.
 * @note This result is suboptimal and we might want to improve on it in the future.
 */
static graphics::boundbox3fd get_boundbox_from_node( const e57::StructureNode& bbox,
                                                     const graphics::transform4fd& transform, bool& hasBounds ) {
    graphics::vector3fd min, max;
    graphics::boundbox3fd transformedBbox;
    e57::NodeType type = bbox.get( "xMinimum" ).type();
    if( type == e57::E57_SCALED_INTEGER ) {
        min.x = e57::ScaledIntegerNode( bbox.get( "xMinimum" ) ).scaledValue();
        max.x = e57::ScaledIntegerNode( bbox.get( "xMaximum" ) ).scaledValue();
        min.y = e57::ScaledIntegerNode( bbox.get( "yMinimum" ) ).scaledValue();
        max.y = e57::ScaledIntegerNode( bbox.get( "yMaximum" ) ).scaledValue();
        min.z = e57::ScaledIntegerNode( bbox.get( "zMinimum" ) ).scaledValue();
        max.z = e57::ScaledIntegerNode( bbox.get( "zMaximum" ) ).scaledValue();
    } else if( type == e57::E57_INTEGER ) {
        min.x = double( e57::IntegerNode( bbox.get( "xMinimum" ) ).value() );
        max.x = double( e57::IntegerNode( bbox.get( "xMaximum" ) ).value() );
        min.y = double( e57::IntegerNode( bbox.get( "yMinimum" ) ).value() );
        max.y = double( e57::IntegerNode( bbox.get( "yMaximum" ) ).value() );
        min.z = double( e57::IntegerNode( bbox.get( "zMinimum" ) ).value() );
        max.z = double( e57::IntegerNode( bbox.get( "zMaximum" ) ).value() );
    } else if( type == e57::E57_FLOAT ) {
        min.x = e57::FloatNode( bbox.get( "xMinimum" ) ).value();
        max.x = e57::FloatNode( bbox.get( "xMaximum" ) ).value();
        min.y = e57::FloatNode( bbox.get( "yMinimum" ) ).value();
        max.y = e57::FloatNode( bbox.get( "yMaximum" ) ).value();
        min.z = e57::FloatNode( bbox.get( "zMinimum" ) ).value();
        max.z = e57::FloatNode( bbox.get( "zMaximum" ) ).value();
    }
    hasBounds = ( type == e57::E57_SCALED_INTEGER || type == e57::E57_INTEGER || type == e57::E57_FLOAT );
    if( hasBounds ) {
        frantic::graphics::boundbox3fd scannerBounds( min.x, max.x, min.y, max.y, min.z, max.z );
        transformedBbox = transform * scannerBounds;
    }
    return transformedBbox;
}

/**
 * Checks if the given int can be represented in 32 bit float
 */
static bool is_representable_in_float( int64_t value ) {
    // Check if there's enough precision for float32 based on the number of bits in the mantissa.
    // If bigger than the mantissa, we won't be able to represent the least significant digit(s).
    const boost::int64_t absValue = value >= 0 ? value : -value;
    if( absValue >= ( 1 << 23 ) )
        return false;
    // In case the minimum is the absolute smallest int value (the above check would be ill-defined since abs( MIN_INT )
    // is ill-defined)
    if( value == std::numeric_limits<boost::int64_t>::min() )
        return false;
    return true;
}

/**
 * Checks if the given e57 channel is represented or should be represented in doubles
 */
static bool is_double( const e57::Node& channel ) {
    switch( channel.type() ) {
    case e57::E57_INTEGER: {
        const e57::IntegerNode intNode( channel );
        const boost::int64_t maxValue = intNode.maximum();
        const boost::int64_t minValue = intNode.minimum();
        return !is_representable_in_float( minValue ) || !is_representable_in_float( maxValue );
    }
    case e57::E57_FLOAT:
        return e57::FloatNode( channel ).precision() == e57::E57_DOUBLE;
    default:
        return true;
    }
}

/**
 * Checks if the input precision data is or should be represented in doubles
 */
static bool is_position_double( const e57::StructureNode& prototype ) {
    if( prototype.isDefined( "cartesianX" ) ) {
        return is_double( prototype.get( "cartesianX" ) ) || is_double( prototype.get( "cartesianY" ) ) ||
               is_double( prototype.get( "cartesianZ" ) );
    } else {
        return true;
    }
}

e57_particle_istream::e57_particle_istream( const tstring& file, channels::data_type_t positionTypeHint )
    : m_filename( file )
    , m_transformData( true )
    , m_boundsValid( false ) {
    initialize_stream( positionTypeHint );
    set_channel_map( m_onDiskParticleChannelMap );
    m_tempParticleBuffer.resize( m_onDiskParticleChannelMap.structure_size() );
}

e57_particle_istream::e57_particle_istream( const tstring& file, const tstring& guid, bool transformData,
                                            channels::data_type_t positionTypeHint )
    : m_filename( file )
    , m_scanGuid( guid )
    , m_transformData( transformData )
    , m_boundsValid( false ) {
    initialize_stream( positionTypeHint );
    set_channel_map( m_onDiskParticleChannelMap );
    m_tempParticleBuffer.resize( m_onDiskParticleChannelMap.structure_size() );
}

void e57_particle_istream::initialize_stream( channels::data_type_t positionTypeHint ) {
    m_currentParticleIndex = -1;
    m_currentScanIndex = 0;
    m_bufferParticleIndex = 0;
    m_scannerTransformsIndex = 0;
    m_maxColor = 1;
    open_file( positionTypeHint );
    setup_metadata();
}

void e57_particle_istream::setup_metadata() {
    frantic::channels::channel_map newMap;
    prt::add_coordinate_system( newMap );
    prt::length_unit_in_micrometers::add_channel( newMap );
    newMap.end_channel_definition();
    frantic::channels::property_map generalMetadata;
    generalMetadata.set_channel_map( newMap );

    prt::set_coordinate_system( generalMetadata, frantic::graphics::coordinate_system::right_handed_zup );
    prt::length_unit_in_micrometers::set_value( generalMetadata, 1e6 );

    if( m_scannerTransforms.size() != 0 ) {
        prt::set_scanner_transforms( generalMetadata, m_scannerTransforms );
    }

    m_metadata.set_general_metadata( generalMetadata );

    if( m_boundsValid ) {
        frantic::channels::property_map positionMetadata;
        prt::add_channel_extents( positionMetadata, m_onDiskParticleChannelMap[_T("Position")].data_type() );
        prt::set_extents( positionMetadata, m_bounds );
        m_metadata.set_channel_metadata( _T("Position"), positionMetadata );
    }
}

void e57_particle_istream::open_file( channels::data_type_t positionTypeHint ) {
    try {
        boost::mutex::scoped_lock lock( g_e57ImageFileConstructorMutex );
        m_imageFile.reset( new e57::ImageFile( frantic::strings::to_utf8( m_filename ), "r" ) );
        lock.unlock();

        e57::StructureNode readRoot = m_imageFile->root();
        // Confirm root has data3D inside
        if( !readRoot.isDefined( "/data3D" ) ) {
            throw std::runtime_error( "e57_particle_istream: File \"" + frantic::strings::to_string( m_filename ) +
                                      "\" doesnt contain 3D data." );
        }
        e57::Node n = readRoot.get( "/data3D" );
        if( n.type() != e57::E57_VECTOR ) {
            throw std::runtime_error( "e57_particle_istream: File \"" + frantic::strings::to_string( m_filename ) +
                                      "\" has a data3D root that isn't of type E57_VECTOR." );
        }
        e57::VectorNode data3D( n );
        m_scanCount = data3D.childCount();

        m_onDiskParticleChannelMap.reset();

        // Go through all the scans to see if there exists other channels & get m_particleCount
        m_particleCount = 0;

        bool someBoundsDefined = false;
        bool noBoundsUndefined = true;
        m_bounds.set_to_empty();
        m_boundsValid = false;

        float minIntensity = std::numeric_limits<float>::max();
        float maxIntensity = boost::numeric::bounds<float>::lowest();

        bool noneHasIntensityLimits = true;

        m_scannerTransforms.clear();

        bool any64BitPositions = false;

        for( int i = 0; i < m_scanCount; ++i ) {
            e57::StructureNode scan( data3D.get( i ) );

            frantic::tstring guid;
            try {
                guid = frantic::strings::to_tstring( e57::StringNode( scan.get( "guid" ) ).value() );
            } catch( e57::E57Exception e ) {
                if( e.errorCode() != e57::E57_ERROR_PATH_UNDEFINED )
                    throw;
                // Otherwise, there was no guid. That's ok.
            }
            if( m_scanGuid == _T("") || m_scanGuid == guid ) {
                e57::CompressedVectorNode points( scan.get( "points" ) );
                m_particleCount += points.childCount();

                frantic::graphics::transform4fd transform;
                // Get pose data
                if( scan.isDefined( "pose" ) && m_transformData ) {
                    transform = get_transform_from_node( e57::StructureNode( scan.get( "pose" ) ) );
                    m_scannerTransforms.push_back( transform );
                } else {
                    FF_LOG( debug ) << "E57 has no pose" << std::endl;
                }

                e57::StructureNode proto( points.prototype() );

                if( scan.isDefined( "intensityLimits" ) ) {
                    noneHasIntensityLimits = false;
                    e57::StructureNode intensityRange( scan.get( "intensityLimits" ) );
                    if( intensityRange.get( "intensityMaximum" ).type() == e57::E57_SCALED_INTEGER ) {
                        const float tempMax =
                            (float)e57::ScaledIntegerNode( intensityRange.get( "intensityMaximum" ) ).scaledValue();
                        const float tempMin =
                            (float)e57::ScaledIntegerNode( intensityRange.get( "intensityMinimum" ) ).scaledValue();
                        maxIntensity = std::max( maxIntensity, tempMax );
                        minIntensity = std::min( minIntensity, tempMin );
                    } else if( intensityRange.get( "intensityMaximum" ).type() == e57::E57_FLOAT ) {
                        const float tempMax = (float)e57::FloatNode( intensityRange.get( "intensityMaximum" ) ).value();
                        const float tempMin = (float)e57::FloatNode( intensityRange.get( "intensityMinimum" ) ).value();
                        maxIntensity = std::max( maxIntensity, tempMax );
                        minIntensity = std::min( minIntensity, tempMin );
                    } else if( intensityRange.get( "intensityMaximum" ).type() == e57::E57_INTEGER ) {
                        const float tempMax =
                            (float)e57::IntegerNode( intensityRange.get( "intensityMaximum" ) ).value();
                        const float tempMin =
                            (float)e57::IntegerNode( intensityRange.get( "intensityMinimum" ) ).value();
                        maxIntensity = std::max( maxIntensity, tempMax );
                        minIntensity = std::min( minIntensity, tempMin );
                    }
                }
                if( is_position_double( proto ) || !transform.is_identity() ) {
                    any64BitPositions = true;
                }
                if( !m_onDiskParticleChannelMap.has_channel( _T("Intensity") ) && proto.isDefined( "intensity" ) ) {
                    m_onDiskParticleChannelMap.define_channel( _T("Intensity"), 1, channels::data_type_float32 );
                }
                if( !m_onDiskParticleChannelMap.has_channel( _T("Color") ) && proto.isDefined( "colorRed" ) &&
                    proto.isDefined( "colorGreen" ) && proto.isDefined( "colorBlue" ) ) {
                    m_onDiskParticleChannelMap.define_channel( _T("Color"), 3, channels::data_type_float32 );
                }

                e57::ustring norUri;
                if( !m_onDiskParticleChannelMap.has_channel( _T("Normal") ) &&
                    m_imageFile->extensionsLookupPrefix( "nor", norUri ) && proto.isDefined( "nor:normalX" ) &&
                    proto.isDefined( "nor:normalY" ) && proto.isDefined( "nor:normalZ" ) ) {
                    m_onDiskParticleChannelMap.define_channel( _T("Normal"), 3, channels::data_type_float32 );
                }

                // Provide an overestimate bounding box by applying the scanner transform to the scanner
                // cartesianBounds.
                bool hasBounds = false;
                frantic::graphics::boundbox3fd transformedBbox;
                if( scan.isDefined( "cartesianBounds" ) && noBoundsUndefined ) {
                    const e57::StructureNode bbox( scan.get( "cartesianBounds" ) );
                    transformedBbox = get_boundbox_from_node( bbox, transform, hasBounds );
                }

                if( hasBounds ) {
                    if( someBoundsDefined ) {
                        m_bounds += transformedBbox;
                    } else {
                        someBoundsDefined = true;
                        m_bounds = transformedBbox;
                    }
                } else {
                    noBoundsUndefined = false;
                }
            }
        }

        m_onDiskParticleChannelMap.define_channel(
            _T("Position"), 3,
            ( positionTypeHint == channels::data_type_float64 && any64BitPositions ) ? channels::data_type_float64
                                                                                     : channels::data_type_float32 );

        if( noneHasIntensityLimits ) {
            m_intensityCoefficient = 1.0f;
            m_intensityOffset = 0.0f;
        } else {
            m_intensityOffset = minIntensity;
            m_intensityCoefficient = 1.0f / ( maxIntensity - m_intensityOffset );
        }

        m_boundsValid = someBoundsDefined && noBoundsUndefined;

        if( m_scannerTransforms.size() != 0 ) {
            m_onDiskParticleChannelMap.define_channel( _T("ScannerIndex"), 1, channels::data_type_uint32 );
        }

        m_onDiskParticleChannelMap.end_channel_definition( 1, true, true );

        // Add channel accessors for each channel
        m_posAccessor = m_onDiskParticleChannelMap.get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
        if( m_onDiskParticleChannelMap.has_channel( _T("Intensity") ) ) {
            m_intensityAccessor = m_onDiskParticleChannelMap.get_cvt_accessor<float>( _T("Intensity") );
        }
        if( m_onDiskParticleChannelMap.has_channel( _T("Color") ) ) {
            m_colorAccessor = m_onDiskParticleChannelMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Color") );
        }
        if( m_onDiskParticleChannelMap.has_channel( _T("Normal") ) ) {
            m_normalAccessor = m_onDiskParticleChannelMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Normal") );
        }
        if( m_scannerTransforms.size() != 0 ) {
            m_scannerIndexAccessor = m_onDiskParticleChannelMap.get_accessor<boost::uint32_t>( _T("ScannerIndex") );
        }

        // Get a scan in the file loaded
        if( m_scanGuid == _T("") )
            load_new_scan();
        else
            load_requested_scan();

    } catch( e57::E57Exception& ex ) { // Can happen if there are errors in the file
        std::string what = ex.what();
        std::string description( e57::Utilities::errorCodeToString( ex.errorCode() ) );
        throw std::runtime_error( "e57_particle_istream: Problem occured with file \"" +
                                  frantic::strings::to_string( m_filename ) + "\": " + what +
                                  " with error: " + description );
    } catch( std::exception& ex ) {
        std::string what = ex.what(); // Can't add the const char* returned from ex.what() in the runtime error throw
        throw std::runtime_error( "e57_particle_istream: Problem occured with file \"" +
                                  frantic::strings::to_string( m_filename ) + "\": " + what );
    } catch( ... ) {
        throw std::runtime_error( "e57_particle_istream: Got an unknown exception while trying to open File \"" +
                                  frantic::strings::to_string( m_filename ) + "\"" );
    }
}

void e57_particle_istream::load_new_scan() {
    e57::StructureNode readRoot = m_imageFile->root();
    if( !readRoot.isDefined( "/data3D" ) ) {
        throw std::runtime_error( "e57_particle_istream: File \"" + frantic::strings::to_string( m_filename ) +
                                  "\" doesnt contain 3D data." );
    }
    e57::Node n = readRoot.get( "/data3D" );
    if( n.type() != e57::E57_VECTOR ) {
        throw std::runtime_error( "e57_particle_istream: File \"" + frantic::strings::to_string( m_filename ) +
                                  "\" has a data3D root that isn't of type E57_VECTOR." );
    }

    if( m_currentScanIndex < m_scanCount ) {
        load_buffers();
    }
}

void e57_particle_istream::load_requested_scan() {
    e57::StructureNode readRoot = m_imageFile->root();
    if( !readRoot.isDefined( "/data3D" ) ) {
        throw std::runtime_error( "e57_particle_istream: File \"" + strings::to_string( m_filename ) +
                                  "\" doesnt contain 3D data." );
    }
    e57::Node n = readRoot.get( "/data3D" );
    if( n.type() != e57::E57_VECTOR ) {
        throw std::runtime_error( "e57_particle_istream: File \"" + strings::to_string( m_filename ) +
                                  "\" has a data3D root that isn't of type E57_VECTOR." );
    }
    e57::VectorNode data3D( n );

    bool scanFound = false;
    for( int i = 0; i < m_scanCount; ++i ) {
        e57::StructureNode scan( data3D.get( i ) );
        if( scan.isDefined( "guid" ) ) {
            if( e57::StringNode( scan.get( "guid" ) ).value() == strings::to_string( m_scanGuid ) ) {
                m_currentScanIndex = i;
                scanFound = true;
                break;
            }
        }
    }
    if( !scanFound ) {
        throw std::runtime_error( "e57_particle_istream: Scan guid \"" + strings::to_string( m_scanGuid ) +
                                  "\" not found in file." );
    }
    load_buffers();
}

void e57_particle_istream::load_buffers() {
    // We already checked that this is OK in load_new_scan or load_requested_scan
    e57::VectorNode data3D( m_imageFile->root().get( "/data3D" ) );
    // Get scan from "/data3D", assume its a Structure (else get exception)
    e57::StructureNode scan( data3D.get( m_currentScanIndex ) );

    if( scan.isDefined( "pose" ) && m_transformData ) {
        m_hasTransform = true;
        ++m_scannerTransformsIndex;
        assert( m_scannerTransformsIndex <= m_scannerTransforms.size() &&
                "m_scannerTransformsIndex should never be greater than m_scannerTransforms.size()" );
    } else {
        m_hasTransform = false;
    }

    // Get "points" field in scan.  Should be a CompressedVectorNode.
    e57::CompressedVectorNode points( scan.get( "points" ) );

    e57::StructureNode proto( points.prototype() );

    if( scan.isDefined( "colorLimits" ) ) {
        e57::StructureNode colorLimits( scan.get( "colorLimits" ) );
        e57::IntegerNode colorRedMaximum( colorLimits.get( "colorRedMaximum" ) );
        // Only check the Red color maximum because they should all be stored same way
        m_maxColor = (float)colorRedMaximum.value();
    } else if( proto.isDefined( "colorRed" ) && proto.isDefined( "colorGreen" ) && proto.isDefined( "colorBlue" ) &&
               proto.get( "colorRed" ).type() == e57::E57_INTEGER ) {
        // Colors are scaled by 1 / m_maxColor
        m_maxColor = 255;
    }

    m_destBuffers.clear();

    // Make sure there are position points
    if( proto.isDefined( "cartesianX" ) && proto.isDefined( "cartesianY" ) && proto.isDefined( "cartesianZ" ) ) {
        m_coordSys = cartesian_coord;
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "cartesianX", m_x, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "cartesianY", m_y, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "cartesianZ", m_z, N, true, true ) );
    } else if( proto.isDefined( "sphericalRange" ) && proto.isDefined( "sphericalAzimuth" ) &&
               proto.isDefined( "sphericalElevation" ) ) {
        m_coordSys = spherical_coord;
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "sphericalRange", m_x, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "sphericalAzimuth", m_y, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "sphericalElevation", m_z, N, true, true ) );
    } else {
        throw std::runtime_error( "e57_particle_istream: A scan in File \"" +
                                  frantic::strings::to_string( m_filename ) +
                                  "\" doesn't contain Cartesian or Spherical points." );
    }

    if( proto.isDefined( "intensity" ) ) {
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "intensity", m_intensity, N, true, true ) );
        m_hasIntensity = true;
    } else {
        m_hasIntensity = false;
    }

    if( proto.isDefined( "colorRed" ) && proto.isDefined( "colorGreen" ) && proto.isDefined( "colorBlue" ) ) {
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "colorRed", m_red, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "colorGreen", m_green, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "colorBlue", m_blue, N, true, true ) );
        m_hasColor = true;
    } else {
        m_hasColor = false;
    }

    /// NOTE: we call `extensionsLookupPrefix` before `isDefined` because
    /// `isDefined` checks XML namespace extensions and will throw an exception
    /// if the `nor` namespace wasn't enabled when the file was saved.
    e57::ustring norUri;
    if( m_imageFile->extensionsLookupPrefix( "nor", norUri ) && proto.isDefined( "nor:normalX" ) &&
        proto.isDefined( "nor:normalY" ) && proto.isDefined( "nor:normalZ" ) ) {
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "nor:normalX", m_normalX, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "nor:normalY", m_normalY, N, true, true ) );
        m_destBuffers.push_back( e57::SourceDestBuffer( *m_imageFile, "nor:normalZ", m_normalZ, N, true, true ) );
        m_hasNormal = true;
    } else {
        m_hasNormal = false;
    }

    // Create a reader of the points CompressedVector. The read calls all occur in get_particle
    m_reader.reset();
    m_reader.reset( new e57::CompressedVectorReader( points.reader( m_destBuffers ) ) );
    m_bufferParticleCount = 0;
    m_bufferParticleIndex = 0;
}

void e57_particle_istream::close() {
    if( m_imageFile->isOpen() ) {
        m_imageFile->close();
    }
    m_particleCount = 0;
    m_scanCount = 0;
}

void e57_particle_istream::set_channel_map( const channels::channel_map& particleChannelMap ) {
    std::vector<char> newDefaultParticle( particleChannelMap.structure_size() );
    if( m_defaultParticleBuffer.size() > 0 ) {
        frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
        defaultAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticleBuffer[0] );
    } else {
        memset( &newDefaultParticle[0], 0, particleChannelMap.structure_size() );
    }
    m_defaultParticleBuffer.swap( newDefaultParticle );

    // Set the map and the adaptor
    m_particleChannelMap = particleChannelMap;
    m_pcmAdaptor.set( m_particleChannelMap, m_onDiskParticleChannelMap );
}

void e57_particle_istream::set_default_particle( char* buffer ) {
    m_particleChannelMap.copy_structure( &m_defaultParticleBuffer[0], buffer );
}

bool e57_particle_istream::get_particle( char* rawParticleBuffer ) {
    if( m_currentScanIndex >= m_scanCount ) {
        if( m_imageFile->isOpen() ) {
            m_imageFile->close();
        }
        return false;
    }

    if( m_bufferParticleIndex >= m_bufferParticleCount ) {
        m_bufferParticleCount = m_reader->read();
        m_bufferParticleIndex = 0;
        // See http://www.libe57.org/bestCoordinates.html for conversion formulas
        switch( m_coordSys ) {
        case cartesian_coord:
            break;
        case spherical_coord:
            for( unsigned i = 0; i < m_bufferParticleCount; ++i ) {
                double sinAzimuth = sin( m_y[i] ), cosAzimuth = cos( m_y[i] );
                double sinElevation = sin( m_z[i] ), cosElevation = cos( m_z[i] );
                double x = m_x[i] * cosElevation * cosAzimuth;
                double y = m_x[i] * cosElevation * sinAzimuth;
                double z = m_x[i] * sinElevation;
                m_x[i] = x;
                m_y[i] = y;
                m_z[i] = z;
            }
            break;
        }
    }

    if( m_bufferParticleCount <= 0 ) {
        if( m_scanGuid == _T("") ) {
            ++m_currentScanIndex;
            m_reader->close();
            load_new_scan();
            return get_particle( rawParticleBuffer );
        } else {
            // Current scan exhausted, but we only wanted the one
            if( m_imageFile->isOpen() )
                m_imageFile->close();
            return false;
        }
    }

    memset( &m_tempParticleBuffer[0], 0, m_onDiskParticleChannelMap.structure_size() );

    frantic::graphics::vector3fd position( m_x[m_bufferParticleIndex], m_y[m_bufferParticleIndex],
                                           m_z[m_bufferParticleIndex] );
    // Apply the transform matrix to the particle
    if( m_hasTransform && m_transformData ) {
        // [m_scannerTransformsIndex - 1] be m_scannerTransformsIndex is 1 baced and vectors are 0 baced.
        position = m_scannerTransforms[m_scannerTransformsIndex - 1].projection_transform( position );
    }

    // Add particle position to the temp buffer
    m_posAccessor.set( m_tempParticleBuffer, position );

    // If intensity is a channel add to the temp buffer
    if( m_hasIntensity ) {
        float intensity = m_intensity[m_bufferParticleIndex];
        intensity = ( intensity - m_intensityOffset ) * m_intensityCoefficient;
        m_intensityAccessor.set( m_tempParticleBuffer, intensity );
    }

    // If color is a channel add to the temp buffer
    if( m_hasColor ) {
        frantic::graphics::vector3f color( m_red[m_bufferParticleIndex], m_green[m_bufferParticleIndex],
                                           m_blue[m_bufferParticleIndex] );
        color /= m_maxColor;
        m_colorAccessor.set( m_tempParticleBuffer, color );
    }

    // If normal is a channel add to the temp buffer
    if( m_hasNormal ) {
        frantic::graphics::vector3f normal( m_normalX[m_bufferParticleIndex], m_normalY[m_bufferParticleIndex],
                                            m_normalZ[m_bufferParticleIndex] );
        m_normalAccessor.set( m_tempParticleBuffer, normal );
    }

    // Add scanner index to the temp buffer
    if( m_scannerTransforms.size() != 0 ) {
        boost::uint32_t currScannerIndex;
        if( m_hasTransform ) {
            currScannerIndex = m_scannerTransformsIndex;
        } else {
            currScannerIndex = 0; // no transform
        }

        m_scannerIndexAccessor.get( m_tempParticleBuffer ) = currScannerIndex;
    }

    // Copy the temp buffer over to rawParticleBuffer
    m_pcmAdaptor.copy_structure( rawParticleBuffer, &m_tempParticleBuffer[0], &m_defaultParticleBuffer[0] );

    ++m_bufferParticleIndex;
    ++m_currentParticleIndex;

    return true;
}

bool e57_particle_istream::get_particles( char* particleBuffer, std::size_t& numParticles ) {
    std::size_t particleSize = m_particleChannelMap.structure_size();
    for( std::size_t i = 0; i < numParticles; ++i ) {
        if( !get_particle( particleBuffer + i * particleSize ) ) {
            numParticles = i;
            return false;
        }
    }

    return true;
}

const particles::particle_file_metadata& e57_particle_istream::get_metadata() const { return m_metadata; }

#endif
