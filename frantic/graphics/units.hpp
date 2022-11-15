// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/transform4f.hpp>
#include <frantic/strings/tstring.hpp>

// TODO: A different namespace?
namespace frantic {
namespace graphics {

namespace length_unit {
enum option { inches, feet, yards, miles, millimeters, centimeters, meters, kilometers, invalid };

struct unit_data {
    const frantic::tchar* name;
    const frantic::tchar* shortName;
    double scaleToMeters;
    double scaleToMicrometers;
};

inline const unit_data& get_unit_data( option unit ) {
    static const unit_data g_unitData[] = { { _T("inches"), _T("in"), 0.0254, 25400.0 },
                                            { _T("feet"), _T("ft"), 0.3048, 304800.0 },
                                            { _T("yards"), _T("yd"), 0.9144, 914400.0 },
                                            { _T("miles"), _T("mi"), 1609.34, 1609344000.0 },
                                            { _T("millimeters"), _T("mm"), 0.001, 1.0e+3 },
                                            { _T("centimeters"), _T("cm"), 0.01, 1.0e+4 },
                                            { _T("meters"), _T("m"), 1.0, 1.0e+6 },
                                            { _T("kilometers"), _T("km"), 1000.0, 1.0e+9 },
                                            { _T("invalid"), _T(""), 0.0, 0.0 } };

    return g_unitData[unit];
}

inline const frantic::tchar* to_string( option unit ) { return get_unit_data( unit ).name; }

inline const frantic::tchar* to_abbreviation( option unit ) { return get_unit_data( unit ).shortName; }

inline double to_meters( option unit ) { return get_unit_data( unit ).scaleToMeters; }

inline double to_micrometers( option unit ) { return get_unit_data( unit ).scaleToMicrometers; }

/**
 * Formats a numeric length of the given units as a string, using the provided number of precision digits to guide the
 * number of decimal places. The units for autoscaling is determined by whether the length unit is SI or imperial.
 *
 * \param length  The numeric value of the length to format.
 * \param lengthUnit  What unit the length represents. If this is an SI unit, the formatting will use SI units for
 *                    autoscaling, and if it is an imperial unit, inches/feet/miles will be used.
 * \param precisionDigits  How many digits of precision to use in formatting the length.
 */
frantic::strings::tstring format_autoscale( double length, option lengthUnit, int precisionDigits );
} // namespace length_unit

namespace coordinate_system {
/**
 * The various coordinate systems that PRT data could be saved in. Knowing this allows us to maintain the expected
 * orienation of 3D data when using PRT files to transfer data between host programs that use differernt coordinate
 * systems.
 */
enum option { unspecified, right_handed_yup, right_handed_zup, left_handed_yup, left_handed_zup, invalid };

/**
 * Some helper names that translate to the actual coordinate system used by default in the named product.
 */
enum {
    max = right_handed_zup,
    maya = right_handed_yup,
    xsi = right_handed_yup,
    houdini = right_handed_yup,
    cinema4d = left_handed_yup,
    realflow = left_handed_yup,
};

/**
 * \return A human readable string corresponding to the specified coordinate system
 */
inline const frantic::tchar* to_string( option sysType ) {
    static const frantic::tchar* names[] = { _T("unspecified"),     _T("right_handed_yup"), _T("right_handed_zup"),
                                             _T("left_handed_yup"), _T("left_handed_zup"),  _T("invalid") };

    return names[sysType];
}

/**
 * Creates a transformation that converts from one coordinate system to another.
 * \param outTM This matrix will be modified to apply the transformation between coordinate systems. Should be identity
 * beforehand. \param from The coordinate system to convert from \param to The coordiante system to convert to \return
 * True if a valid transformation was possible between the coordinate systems. If not, we return false and do not affect
 * outTM.
 */
template <typename FloatType>
inline bool create_transform( frantic::graphics::transform4t<FloatType>& outTM, option from, option to ) {
    if( from >= invalid || to >= invalid )
        throw std::logic_error( "Invalid length_unit type for length_unit::create_transform()" );
    if( from <= unspecified || to <= unspecified )
        return false;

    // This table contains a compressed representation of the transform required. It is indexed like a 2 dimensional
    // matrix with first index being "from". Legend: bit 0 is set if we need to swap y/z axis. bit 1 set if we then
    // negate y. bit 2 set if we then negate z.
    static const int bits[][4] = { { 0, ( 1 << 0 ) | ( 1 << 2 ), ( 1 << 2 ), ( 1 << 0 ) },
                                   { ( 1 << 0 ) | ( 1 << 1 ), 0, ( 1 << 0 ), ( 1 << 1 ) },
                                   { ( 1 << 2 ), ( 1 << 0 ), 0, ( 1 << 0 ) | ( 1 << 2 ) },
                                   { ( 1 << 0 ), ( 1 << 1 ), ( 1 << 0 ) | ( 1 << 1 ), 0 } };

    int curBits = bits[from - right_handed_yup][to - right_handed_yup];

    if( curBits & ( 1 << 0 ) )
        std::swap( outTM.get_column( 1 ), outTM.get_column( 2 ) );
    if( curBits & ( 1 << 1 ) )
        outTM.get_column( 1 ) = -outTM.get_column( 1 );
    if( curBits & ( 1 << 2 ) )
        outTM.get_column( 2 ) = -outTM.get_column( 2 );

    return true;
}
} // namespace coordinate_system

} // namespace graphics
} // namespace frantic
