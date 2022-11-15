// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/prt_metadata.hpp>

#include <frantic/channels/channel_column_map.hpp>
#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/named_channel_data.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/particles/particle_file_metadata.hpp>

#include <boost/foreach.hpp>

namespace frantic {
namespace particles {
namespace prt {

using namespace std;
using namespace frantic;
using namespace frantic::particles;

void set_channel_interpretation( frantic::channels::property_map& props, channel_interpretation::option val ) {
    set_channel_interpretation( props, channel_interpretation::to_string( val ) );
}

void set_channel_interpretation( frantic::channels::property_map& props, const frantic::tstring& val ) {
    props.set_cvt<tstring>( _T("Interpretation"), val );
}

namespace channel_interpretation {

option parse( const frantic::tstring& interp ) {
    for( int i = 0; i < channel_interpretation::invalid; ++i ) {
        if( interp == channel_interpretation::to_string( static_cast<channel_interpretation::option>( i ) ) ) {
            return static_cast<channel_interpretation::option>( i );
        }
    }

    return interp.empty() ? channel_interpretation::unspecified : channel_interpretation::invalid;
}

} // namespace channel_interpretation

channel_interpretation::option get_channel_interpretation( const frantic::channels::property_map& props ) {
    if( !props.has_property( _T("Interpretation") ) )
        return prt::channel_interpretation::unspecified;

    channels::data_type_t dataType;
    size_t arity;
    props.get_channel_map().get_channel_definition( _T("Interpretation"), dataType, arity );
    if( dataType == channels::data_type_string ) {
        tstring interp = props.get_cvt<tstring>( _T("Interpretation") );
        return channel_interpretation::parse( interp );
    } else {
        return static_cast<channel_interpretation::option>( props.get<boost::int32_t>( _T("Interpretation") ) );
    }
}

bool has_channel_range( const frantic::channels::property_map& props ) {
    return props.has_property( _T( "ChannelRange" ) );
}

void add_channel_range_property( frantic::channels::property_map& props, frantic::channels::data_type_t type ) {
    if( !has_channel_range( props ) ) {
        frantic::channels::channel_map newMap;
        newMap.union_channel_map( props.get_channel_map() );
        newMap.define_channel( _T( "ChannelRange" ), 2u, type );
        newMap.end_channel_definition();
        props.set_channel_map_with_swap( newMap );
    }
}

void get_scanner_transforms( const frantic::channels::property_map& propertyMap,
                             std::vector<frantic::graphics::transform4fd>& outScannerTransforms ) {
    const frantic::tstring propertyName( _T("ScannerTransforms") );

    outScannerTransforms.clear();

    const frantic::channels::channel_map& channelMap = propertyMap.get_channel_map();

    if( channelMap.has_channel( propertyName ) ) {
        const std::size_t arity = channelMap[propertyName].arity();
        if( arity % 16 ) {
            throw std::runtime_error( "get_scanner_transforms Error: unexpected arity for ScannerTransforms.  Expected "
                                      "a multiple of 16, but got " +
                                      boost::lexical_cast<std::string>( arity ) + " instead." );
        }
        if( arity == 0 ) {
            return;
        }

        // convert input property data to double precision
        frantic::channels::channel_type_convertor_function_t cvt =
            frantic::channels::get_channel_type_convertor_function(
                channelMap[propertyName].data_type(), frantic::channels::data_type_float64, propertyName );

        std::vector<double> doubleBuffer;
        doubleBuffer.resize( arity );

        cvt( reinterpret_cast<char*>( &doubleBuffer[0] ),
             propertyMap.get_raw_buffer() + channelMap[propertyName].offset(), arity );

        const std::size_t scannerTransformCount = arity / 16;
        for( std::size_t scannerTransformIndex = 0; scannerTransformIndex < scannerTransformCount;
             ++scannerTransformIndex ) {
            outScannerTransforms.push_back(
                frantic::graphics::transform4fd( &doubleBuffer[16 * scannerTransformIndex] ) );
        }
    }
}

void set_scanner_transforms( frantic::channels::property_map& outPropertyMap,
                             const std::vector<frantic::graphics::transform4fd>& scannerTransforms ) {
    const frantic::tstring propertyName( _T("ScannerTransforms") );

    // Remove ScannerTransforms property.
    // We'll recreate it with the correct size and data type if it is needed.
    if( outPropertyMap.has_property( propertyName ) ) {
        frantic::channels::channel_map channelMap = outPropertyMap.get_channel_map();
        channelMap.delete_channel( propertyName );
        outPropertyMap.set_channel_map( channelMap );
    }

    if( scannerTransforms.size() > 0 ) {
        frantic::channels::channel_map newChannelMap = outPropertyMap.get_channel_map();

        const std::size_t arity = 16 * scannerTransforms.size();
        const frantic::channels::data_type_t dataType = frantic::channels::data_type_float64;
        newChannelMap.append_channel( propertyName, arity, dataType );

        outPropertyMap.set_channel_map( newChannelMap );

        double* outBuffer =
            reinterpret_cast<double*>( outPropertyMap.get_raw_buffer() + newChannelMap[propertyName].offset() );

        for( std::size_t scannerIndex = 0, scannerIndexEnd = scannerTransforms.size(); scannerIndex != scannerIndexEnd;
             ++scannerIndex ) {
            for( int i = 0; i < 16; ++i ) {
                outBuffer[16 * scannerIndex + i] = scannerTransforms[scannerIndex][i];
            }
        }
    }
}

frantic::graphics::transform4f
get_unit_length_coordinate_transform( const frantic::channels::property_map& props, double targetScaleToMeters,
                                      frantic::graphics::coordinate_system::option targetCoordSys ) {
    const double scaleToMeters = length_unit_in_meters::get_value( props );
    const frantic::graphics::coordinate_system::option fromCoordSys = get_coordinate_system( props );

    // Compute the matrix that accomodates units & coordinate system changes. Order doesn't matter due to the nature of
    // respective matrices.
    frantic::graphics::transform4f tm;
    frantic::graphics::coordinate_system::create_transform( tm, fromCoordSys, targetCoordSys );
    prt::length_unit_in_micrometers::create_transform( tm, scaleToMeters, targetScaleToMeters );
    return tm;
}

frantic::tstring convert_channel_interpretation_prt1_to_prt2( int interpretationValue ) {
    if( interpretationValue == channel_interpretation::unspecified )
        return _T("");
    if( interpretationValue < 0 || interpretationValue >= channel_interpretation::invalid )
        return channel_interpretation::to_string( channel_interpretation::invalid );
    return channel_interpretation::to_string( static_cast<channel_interpretation::option>( interpretationValue ) );
}

frantic::channels::property_map
convert_channel_metadata_prt1_to_prt2( const frantic::channels::property_map& inMetadata ) {
    frantic::channels::property_map pm = inMetadata;
    if( pm.has_property( _T("Interpretation") ) ) {
        const boost::int32_t interpretationValue = pm.get_cvt<boost::int32_t>( _T("Interpretation") );
        pm.delete_property( _T("Interpretation") );
        pm.set_property( _T("Interpretation"), convert_channel_interpretation_prt1_to_prt2( interpretationValue ) );
    }

    return pm;
}

frantic::channels::property_map
convert_channel_metadata_prt2_to_prt1( const frantic::channels::property_map& inMetadata ) {
    frantic::channels::property_map pm = inMetadata;
    if( pm.has_property( _T("Interpretation") ) ) {
        const tstring interpretationValue = pm.get_cvt<tstring>( _T("Interpretation") );
        pm.delete_property( _T("Interpretation") );
        pm.set_property( _T("Interpretation"), boost::int32_t( channel_interpretation::parse( interpretationValue ) ) );
    }

    return pm;
}

frantic::particles::particle_file_metadata
convert_metadata_prt1_to_prt2( const frantic::particles::particle_file_metadata& inMetadata ) {
    particle_file_metadata outMetadata;

    // Set up the per channel metadata
    std::vector<frantic::tstring> channelNames;
    inMetadata.get_channels_with_metadata( channelNames );
    BOOST_FOREACH( const frantic::tstring& ch, channelNames ) {
        frantic::channels::property_map* outChannelMetadata = outMetadata.get_channel_metadata( ch, true );
        if( !outChannelMetadata ) {
            throw std::runtime_error( "convert_metadata_prt1_to_prt2 Internal Error: output property map is NULL" );
        }

        const frantic::channels::property_map* inChannelMetadata = inMetadata.get_channel_metadata( ch );
        if( !inChannelMetadata ) {
            throw std::runtime_error( "convert_metadata_prt1_to_prt2 Internal Error: input property map is NULL" );
        }

        outChannelMetadata->merge_property_map( convert_channel_metadata_prt1_to_prt2( *inChannelMetadata ) );
    }

    // Set up the global available metadata
    const frantic::channels::property_map& pm = inMetadata.get_general_metadata();
    const frantic::channels::channel_map& cm = pm.get_channel_map();
    frantic::channels::channel_map generalChannelMap;
    for( size_t i = 0; i < cm.channel_count(); ++i ) {
        if( cm[i].name() == _T("BoundBox") ) {
            // BoundBox is Position.Extents in PRT2
            frantic::channels::channel_map extentsChannelMap;
            extentsChannelMap.define_channel( _T("Extents"), cm[i].arity(), cm[i].data_type() );
            extentsChannelMap.end_channel_definition();
            frantic::channels::property_map extentsPropertyMap( extentsChannelMap );

            const char* src = pm.get_channel_buffer( _T("BoundBox") );
            char* dst = extentsPropertyMap.get_channel_buffer( _T("Extents") );
            memcpy( dst, src, cm[i].primitive_size() );

            outMetadata.append_channel_metadata( _T("Position"), extentsPropertyMap, false );

        } else if( cm[i].name() == _T("LengthUnitInMeters") ) {
            // LengthUnitInMeters is replaced by LengthUnitInMicrometers in PRT2
            if( !pm.has_property( _T("LengthUnitInMicrometers") ) ) {
                generalChannelMap.define_channel( _T("LengthUnitInMicrometers"), cm[i].arity(), cm[i].data_type() );
            }

        } else {
            // Allow all other properties to pass through
            // (in particular: CoordSys is unchanged)
            generalChannelMap.define_channel( cm[i].name(), cm[i].arity(), cm[i].data_type() );
        }
    }
    generalChannelMap.end_channel_definition();

    // Set up the actual metadata values
    frantic::channels::property_map generalMetadata( generalChannelMap );
    for( size_t i = 0; i < cm.channel_count(); ++i ) {
        if( cm[i].name() == _T("BoundBox") ) {
            // Bounding box has been dealt with already above, ignore it

        } else if( cm[i].name() == _T("LengthUnitInMeters") ) {
            // Translate unit length
            generalMetadata.set_cvt<double>( _T("LengthUnitInMicrometers"),
                                             1.0e6 * pm.get_cvt<double>( _T("LengthUnitInMeters") ) );

        } else {
            // Pass through others
            const char* src = pm.get_channel_buffer( cm[i].name() );
            char* dst = generalMetadata.get_channel_buffer( cm[i].name() );
            memcpy( dst, src, cm[i].primitive_size() );
        }
    }
    outMetadata.set_general_metadata( generalMetadata );
    return outMetadata;
}

frantic::particles::particle_file_metadata
convert_metadata_prt2_to_prt1( const frantic::particles::particle_file_metadata& inMetadata ) {
    particle_file_metadata outMetadata;

    // Set up the per channel metadata
    std::vector<frantic::tstring> channelNames;
    inMetadata.get_channels_with_metadata( channelNames );

    bool hasExtents = false;
    std::size_t extentsArity = 0;
    frantic::channels::data_type_t extentsType = frantic::channels::data_type_invalid;
    std::vector<char> bboxData;

    BOOST_FOREACH( const frantic::tstring& ch, channelNames ) {
        frantic::channels::property_map* outChannelMetadata = outMetadata.get_channel_metadata( ch, true );
        if( !outChannelMetadata ) {
            throw std::runtime_error(
                "convert_channel_metadata_prt2_to_prt1 Internal Error: output property map is NULL" );
        }

        const frantic::channels::property_map* inChannelMetadata = inMetadata.get_channel_metadata( ch );
        if( !inChannelMetadata ) {
            throw std::runtime_error(
                "convert_channel_metadata_prt2_to_prt1 Internal Error: input property map is NULL" );
        }

        if( ch == _T("Position") && inChannelMetadata->has_property( _T("Extents") ) ) {
            // Position.Extents is BoundBox in PRT1
            hasExtents = true;
            const channels::channel_map& sourceMap = inChannelMetadata->get_channel_map();
            bboxData.resize( sourceMap[_T("Extents")].primitive_size() );
            sourceMap.get_channel_definition( _T("Extents"), extentsType, extentsArity );

            const char* src = inChannelMetadata->get_channel_buffer( _T("Extents") );
            char* dst = &bboxData[0];
            memcpy( dst, src, bboxData.size() );
        }

        outChannelMetadata->merge_property_map( convert_channel_metadata_prt2_to_prt1( *inChannelMetadata ) );
    }

    // Set up the global available metadata
    const frantic::channels::property_map& pm = inMetadata.get_general_metadata();
    const frantic::channels::channel_map& cm = pm.get_channel_map();
    frantic::channels::channel_map generalChannelMap, extentsChannelMap;
    for( size_t i = 0; i < cm.channel_count(); ++i ) {
        if( cm[i].name() == _T("LengthUnitInMicrometers") ) {
            // LengthUnitInMicrometers is replaced by LengthUnitInMeters in PRT1
            if( !pm.has_property( _T("LengthUnitInMeters") ) ) {
                generalChannelMap.define_channel( _T("LengthUnitInMeters"), cm[i].arity(), cm[i].data_type() );
            }

        } else {
            // Allow all other properties to pass through
            // (in particular: CoordSys is unchanged)
            generalChannelMap.define_channel( cm[i].name(), cm[i].arity(), cm[i].data_type() );
        }
    }
    if( hasExtents ) {
        // BoundBox data from Extents
        if( !pm.has_property( _T("BoundBox") ) ) {
            generalChannelMap.define_channel( _T("BoundBox"), extentsArity, extentsType );
        }
    }

    generalChannelMap.end_channel_definition();

    // Set up the actual metadata values
    frantic::channels::property_map generalMetadata( generalChannelMap );
    for( size_t i = 0; i < cm.channel_count(); ++i ) {
        if( cm[i].name() == _T("LengthUnitInMicrometers") ) {
            // Translate unit length
            generalMetadata.set_cvt<double>( _T("LengthUnitInMeters"),
                                             1.0e-6 * pm.get_cvt<double>( _T("LengthUnitInMicrometers") ) );

        } else {
            // Pass through others
            const char* src = pm.get_channel_buffer( cm[i].name() );
            char* dst = generalMetadata.get_channel_buffer( cm[i].name() );
            memcpy( dst, src, cm[i].primitive_size() );
        }
    }
    if( hasExtents ) {
        // Set bounding box data
        const char* src = &bboxData[0];
        char* dst = generalMetadata.get_channel_buffer( _T("BoundBox") );
        memcpy( dst, src, bboxData.size() );
    }
    outMetadata.set_general_metadata( generalMetadata );
    return outMetadata;
}

frantic::tstring get_default_prt2_channel_interpretation( const frantic::tstring& channelName ) {
    if( channelName == _T("Position") || channelName == _T("BirthPosition") ) {
        return _T("Point");
    } else if( channelName == _T("Velocity") || channelName == _T("Acceleration") ) {
        return _T("Vector");
    } else if( channelName == _T("Normal") || channelName == _T("Tangent") || channelName == _T("Binormal") ) {
        return _T("Normal");
    } else if( channelName == _T("Orientation") ) {
        return _T("Orientation");
    } else if( channelName == _T("Spin") ) {
        return _T("Rotation");
    } else if( channelName == _T("Radius") ) {
        return _T("Scalar");
    }
    return _T("");
}

void add_default_prt2_channel_interpretation_data( frantic::particles::particle_file_metadata& metadata,
                                                   const frantic::channels::channel_map& channelMap ) {
    for( size_t i = 0; i < channelMap.channel_count(); ++i ) {
        frantic::tstring value( get_default_prt2_channel_interpretation( channelMap[i].name() ) );
        if( value != _T("") ) {
            frantic::channels::property_map* channelMetadata =
                metadata.get_channel_metadata( channelMap[i].name(), true );
            if( !channelMetadata->has_property( _T("Interpretation") ) ) {
                channelMetadata->set_property( _T("Interpretation"), value );
            }
        }
    }
}

} // namespace prt
} // namespace particles
} // namespace frantic
