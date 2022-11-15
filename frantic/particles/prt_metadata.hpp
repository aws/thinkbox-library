// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_column_map.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/particles/particle_file_metadata.hpp>
#include <frantic/strings/tstring.hpp>

#include <boost/math/special_functions/next.hpp>

namespace frantic {
namespace particles {
namespace prt {

namespace length_unit_in_meters {
inline bool is_equivalent( double lhsScaleToMeters, double rhsScaleToMeters ) {
    // FIXME: Magic number (absolute, not relative)
    return std::abs( lhsScaleToMeters - rhsScaleToMeters ) < 1e-5;
}

inline const frantic::tchar* to_string( double scaleToMeters ) {
    // Search for a known unit exactly.
    for( std::size_t i = 0, iEnd = frantic::graphics::length_unit::invalid; i < iEnd; ++i ) {
        frantic::graphics::length_unit::option unit = static_cast<frantic::graphics::length_unit::option>( i );
        if( is_equivalent( frantic::graphics::length_unit::to_meters( unit ), scaleToMeters ) )
            return frantic::graphics::length_unit::to_string( unit );
    }

    // TODO: Could also find multiples of known units (ie. divide them and see if its a whole number result)

    return _T("");
}

template <typename FloatType>
inline bool create_transform( frantic::graphics::transform4t<FloatType>& outTM, double srcScaleToMeters,
                              double destScaleToMeters ) {
    // If either scale is unknown (indicated by 0.0 scale) we cannot convert to meters so do nothing.
    if( srcScaleToMeters == 0.0 || destScaleToMeters == 0.0 )
        return false;

    if( !is_equivalent( srcScaleToMeters, destScaleToMeters ) )
        outTM.scale( static_cast<FloatType>( srcScaleToMeters / destScaleToMeters ) );
    return true;
}

inline void add_channel( frantic::channels::channel_map& inoutMap ) {
    inoutMap.define_channel( _T("LengthUnitInMeters"), 1u, frantic::channels::data_type_float64 );
}

inline void add_channel( frantic::channels::property_map& props ) {
    if( !props.has_property( _T("LengthUnitInMeters") ) ) {
        frantic::channels::channel_map newMap;
        newMap.union_channel_map( props.get_channel_map() );
        newMap.define_channel( _T("LengthUnitInMeters"), 1u, frantic::channels::data_type_float64 );
        newMap.end_channel_definition();

        props.set_channel_map_with_swap( newMap );
    }
}

inline void set_value( frantic::channels::property_map& props, double scaleToMeters ) {
    props.get<double>( _T("LengthUnitInMeters") ) = scaleToMeters;
}

inline void set_value( frantic::channels::property_map& props, frantic::graphics::length_unit::option unit ) {
    set_value( props, frantic::graphics::length_unit::to_meters( unit ) );
}

inline double get_value( const frantic::channels::property_map& props ) {
    if( props.has_property( _T("LengthUnitInMicrometers") ) ) {
        return props.get<double>( _T("LengthUnitInMicrometers") ) / 1.0e6;
    } else if( props.has_property( _T("LengthUnitInMeters") ) ) {
        return props.get<double>( _T("LengthUnitInMeters") );
    } else {
        return 0.f;
    }
}
} // namespace length_unit_in_meters

// For PRT2
namespace length_unit_in_micrometers {
inline bool is_equivalent( double lhsScaleToMicrometers, double rhsScaleToMicrometers ) {
    // TODO: double(0.001f) and 0.001 are 219043332 ULP apart, should we be that lax, or add special float routines?
    return std::abs( boost::math::float_distance( lhsScaleToMicrometers, rhsScaleToMicrometers ) ) < 100;
}

inline const frantic::tchar* to_string( double scaleToMicrometers ) {
    // Search for a known unit exactly.
    for( std::size_t i = 0, iEnd = frantic::graphics::length_unit::invalid; i < iEnd; ++i ) {
        frantic::graphics::length_unit::option unit = static_cast<frantic::graphics::length_unit::option>( i );
        if( is_equivalent( frantic::graphics::length_unit::to_micrometers( unit ), scaleToMicrometers ) )
            return frantic::graphics::length_unit::to_string( unit );
    }

    // TODO: Could also find multiples of known units (ie. divide them and see if its a whole number result)

    return _T("");
}

template <typename FloatType>
inline bool create_transform( frantic::graphics::transform4t<FloatType>& outTM, double srcScaleToMicrometers,
                              double destScaleToMicrometers ) {
    // If either scale is unknown (indicated by 0.0 scale) we cannot convert to meters so do nothing.
    if( srcScaleToMicrometers == 0.0 || destScaleToMicrometers == 0.0 )
        return false;

    if( !is_equivalent( srcScaleToMicrometers, destScaleToMicrometers ) )
        outTM.scale( static_cast<FloatType>( srcScaleToMicrometers / destScaleToMicrometers ) );
    return true;
}

inline void add_channel( frantic::channels::channel_map& inoutMap ) {
    inoutMap.define_channel( _T("LengthUnitInMicrometers"), 1u, frantic::channels::data_type_float64 );
}

inline void add_channel( frantic::channels::property_map& props ) {
    if( !props.has_property( _T("LengthUnitInMicrometers") ) ) {
        frantic::channels::channel_map newMap;
        newMap.union_channel_map( props.get_channel_map() );
        newMap.define_channel( _T("LengthUnitInMicrometers"), 1u, frantic::channels::data_type_float64 );
        newMap.end_channel_definition();

        props.set_channel_map_with_swap( newMap );
    }
}

inline void set_value( frantic::channels::property_map& props, double scaleToMicrometers ) {
    props.get<double>( _T("LengthUnitInMicrometers") ) = scaleToMicrometers;
}

inline void set_value( frantic::channels::property_map& props, frantic::graphics::length_unit::option unit ) {
    set_value( props, frantic::graphics::length_unit::to_micrometers( unit ) );
}

inline double get_value( const frantic::channels::property_map& props ) {
    if( props.has_property( _T("LengthUnitInMicrometers") ) ) {
        return props.get<double>( _T("LengthUnitInMicrometers") );
    } else if( props.has_property( _T("LengthUnitInMeters") ) ) {
        return props.get<double>( _T("LengthUnitInMeters") ) * 1e6;
    } else {
        return 0.f;
    }
}
} // namespace length_unit_in_micrometers

namespace channel_interpretation {
/**
 * Channels in PRT files can be tagged with this enumeration to provide more information about how to interpret the
 * values associated with a channel. This allows us to determine how a transformation applied to a particle set will
 * affect the data within it. Ex. point channels are modified by the full rotation, scale, skew, translation and
 * projection aspects of the transformation.
 */
enum option {
    unspecified,
    point,       // Normal 3D transformation
    vector,      // Not affected by translation
    normal,      // Not affected by translation or skew (maybe not scale either...)
    orientation, // Only affected by rotation
    rotation,    // Only affected by rotation
    scalar,      // Only affected by scale (somehow combining the 3 scale parts into 1 value)
    invalid
};

/**
 * \return A human-readable string, which is also the value used in PRT2 metadata, for the specified option.
 */
inline const frantic::tchar* to_string( option val ) {
    static const frantic::tchar* names[] = { _T("Unspecified"), _T("Point"),    _T("Vector"), _T("Normal"),
                                             _T("Orientation"), _T("Rotation"), _T("Scalar"), _T("Invalid") };

    return names[val];
}

/**
 * \return The option corresponding to the given string, or invalid if no match was found
 */
option parse( const frantic::tstring& s );

/**
 * Determines if a data_type_t & arity pair is compatible with a given channel_interpretation::option. For example a
 * float32[1] channel or an int32[3] channel cannot be interpreted as a point because of incorrect arity and data type
 * respectively. \param val The channel_interpretation::option to check for compatibility with. \param dt The data type
 * of the channel. \param arity The arity of the channel. \return True if the data type and arity are appropriate for
 * the given channel_interpretation::option.
 */
inline bool is_compatible( option val, frantic::channels::data_type_t dt, std::size_t arity ) {
    switch( val ) {
    case point:
    case vector:
    case normal:
        // Intentional fall-through
        return frantic::channels::is_channel_data_type_float( dt ) && arity == 3u;
    case orientation:
    case rotation:
        // Intentional fall-through
        return frantic::channels::is_channel_data_type_float( dt ) && arity == 4u;
    case scalar:
        return frantic::channels::is_channel_data_type_float( dt ) && arity == 1u;
    default:
        return true;
    }
}
} // namespace channel_interpretation

/**
 * Defines an appropriately named channel in the provided channel map to hold a coordinate_system::option value.
 * \param ch The map to modify.
 */
inline void add_coordinate_system( frantic::channels::channel_map& ch ) {
    ch.define_channel( _T("CoordSys"), 1, frantic::channels::data_type_int32 );
}

/**
 * Assigns the coordinate_system::option to the appropriate channel in the property_map. The property_map's channel_map
 * should have been intialized by add_coordinate_system(). \param props The property_map to modify \param val The value
 * to set in the property_map
 */
inline void set_coordinate_system( frantic::channels::property_map& props,
                                   frantic::graphics::coordinate_system::option val ) {
    props.get<boost::int32_t>( _T("CoordSys") ) = val;
}

/**
 * Retrieves a coordinate_system::option from a property_map.
 * \param props The property_map presumed to contain a coordinate_system::option value.
 * \return coordinate_system::unspecified if the property_map did not contain a value, otherwise the contained value.
 */
inline frantic::graphics::coordinate_system::option
get_coordinate_system( const frantic::channels::property_map& props ) {
    if( !props.has_property( _T("CoordSys") ) )
        return frantic::graphics::coordinate_system::unspecified;

    return static_cast<frantic::graphics::coordinate_system::option>( props.get<boost::int32_t>( _T("CoordSys") ) );
}

inline void add_framerate( frantic::channels::channel_map& ch ) {
    ch.define_channel( _T("FrameRate"), 2, frantic::channels::data_type_uint32 );
}

inline void set_framerate( frantic::channels::property_map& props, boost::uint32_t numerator,
                           boost::uint32_t denominator ) {
    props.get<frantic::graphics2d::vector2>( _T("FrameRate") )
        .set( static_cast<int>( numerator ), static_cast<int>( denominator ) );
}

inline bool get_framerate( const frantic::channels::property_map& props, boost::uint32_t& outNumerator,
                           boost::uint32_t& outDenominator ) {
    if( !props.has_property( _T("FrameRate") ) )
        return false;

    const std::pair<boost::uint32_t, boost::uint32_t>& fps =
        props.get<std::pair<boost::uint32_t, boost::uint32_t>>( _T("FrameRate") );

    outNumerator = fps.first;
    outDenominator = fps.second;

    return true;
}

inline void add_channel_interpretation( frantic::channels::channel_map& ch ) {
    // In PRT2, Interpretation is a string
    // If we want PRT1, we'll need to call a converter method (see convert_channel_metadata_prt2_to_prt1)
    ch.define_channel( _T("Interpretation"), 1, frantic::channels::data_type_string );
}

inline void add_channel_interpretation( frantic::channels::property_map& props ) {
    if( !props.get_channel_map().has_channel( _T( "Interpretation" ) ) ) {
        frantic::channels::channel_map newMap;
        newMap.union_channel_map( props.get_channel_map() );
        add_channel_interpretation( newMap );
        newMap.end_channel_definition();
        props.set_channel_map_with_swap( newMap );
    }
}

inline void add_channel_extents( frantic::channels::channel_map& ch, channels::data_type_t type ) {
    ch.define_channel( _T("Extents"), 6, type );
}

inline void add_channel_extents( frantic::channels::property_map& props, channels::data_type_t type ) {
    if( !props.get_channel_map().has_channel( _T( "Extents" ) ) ) {
        frantic::channels::channel_map newMap;
        newMap.union_channel_map( props.get_channel_map() );
        add_channel_extents( newMap, type );
        newMap.end_channel_definition();
        props.set_channel_map_with_swap( newMap );
    }
}

template <typename FloatType>
inline void set_extents( frantic::channels::property_map& props, const graphics::boundbox3t<FloatType>& bounds ) {
    props.set_cvt( _T("Extents"), bounds );
}

/**
 * Sets the channel interpretation metadata value.  Auto converts to an argument appropriate for PRT2 format.
 */
void set_channel_interpretation( frantic::channels::property_map& props, channel_interpretation::option val );

/**
 * Sets the channel interpretation metadata value.  Assumes PRT2 format.
 */
void set_channel_interpretation( frantic::channels::property_map& props, const frantic::tstring& val );

/**
 * Gets the channel interpretation metadata value.  Auto converts to int if the property is a string type (PRT2 Format)
 */
channel_interpretation::option get_channel_interpretation( const frantic::channels::property_map& props );

bool has_channel_range( const frantic::channels::property_map& props );

template <typename T>
std::pair<T, T> get_channel_range( const frantic::channels::property_map& props ) {
    using namespace frantic::channels;

    if( !has_channel_range( props ) ) {
        throw std::runtime_error( "get_channel_range -- no channel range on this channel" );
    }

    channel_static_cast_const_accessor<T> acc( props.get_channel_map().get_general_accessor( _T( "ChannelRange" ) ) );
    const char* buffer = props.get_channel_buffer( _T( "ChannelRange" ) );
    return std::make_pair( acc.get( buffer, 0 ), acc.get( buffer, 1 ) );
}

void add_channel_range_property( frantic::channels::property_map& props, frantic::channels::data_type_t type );

template <typename T>
void set_channel_range( frantic::channels::property_map& props, const std::pair<T, T>& range ) {
    if( !has_channel_range( props ) ) {
        throw std::runtime_error( "get_channel_range -- no channel range on this channel" );
    }
    props.set_cvt<std::pair<T, T>>( _T( "ChannelRange" ), range );
}

template <typename T>
void add_channel_range_property( frantic::channels::property_map& props, const std::pair<T, T>& initial ) {
    add_channel_range_property( props, frantic::channels::channel_data_type_traits<T>::data_type() );
    set_channel_range( props, initial );
}

void get_scanner_transforms( const frantic::channels::property_map& props,
                             std::vector<frantic::graphics::transform4fd>& outScannerTransforms );

void set_scanner_transforms( frantic::channels::property_map& props,
                             const std::vector<frantic::graphics::transform4fd>& scannerTransforms );

frantic::graphics::transform4f
get_unit_length_coordinate_transform( const frantic::channels::property_map& props, double targetScaleToMeters,
                                      frantic::graphics::coordinate_system::option targetCoordSys );

/**
 * Returns the value of the "interpretation" property of the channel following prt2 specification.
 *
 * @param interpretationValue The value of the "interpretation" property following prt1 specification
 */
frantic::tstring convert_channel_interpretation_prt1_to_prt2( boost::int32_t interpretationValue );

/**
 * Returns the property_map with metadata following prt2 specification for a channel.
 *
 * @param inMetadata property_map of a particular channel following prt1 specification
 */
frantic::channels::property_map
convert_channel_metadata_prt1_to_prt2( const frantic::channels::property_map& inMetadata );

/**
 * Returns the property_map with metadata following prt1 specification for a channel.
 *
 * @param inMetadata property_map of a particular channel following prt2 specification
 */
frantic::channels::property_map
convert_channel_metadata_prt2_to_prt1( const frantic::channels::property_map& inMetadata );

/**
 * Returns the particle_file_metadata following prt2 specification
 *
 * @note This function only converts standard prt1 metadata channels.
 *
 * @param inMetadata particle_file_meatdata from a file following prt1 specification
 */
frantic::particles::particle_file_metadata
convert_metadata_prt1_to_prt2( const frantic::particles::particle_file_metadata& inMetadata );

/**
 * Returns the particle_file_metadata following prt1 specification
 *
 * @note This function only converts standard prt2 metadata channels.
 *
 * @param inMetadata particle_file_meatdata from a file following prt2 specification
 */
frantic::particles::particle_file_metadata
convert_metadata_prt2_to_prt1( const frantic::particles::particle_file_metadata& inMetadata );

/**
 * Returns the default value of the "interpretation" property for a given channel.
 *
 * @param channelName Name of the channel for which we want the value of the "interpretation" proeprty.
 */
frantic::tstring get_default_prt2_channel_interpretation( const frantic::tstring& channelName );

/**
 * Adds interpretation property to the channels if the channels already don't have it.
 *
 * @param[in,out] metadata particle_file_metadata which may or may not have the "interpreatation" property for its
 * channels.
 * @param channelMap channel_map defining all the channels that will be written in the prt2 file.
 */
void add_default_prt2_channel_interpretation_data( frantic::particles::particle_file_metadata& metadata,
                                                   const frantic::channels::channel_map& channelMap );

} // namespace prt
} // namespace particles
} // namespace frantic
