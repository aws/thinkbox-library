// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics/units.hpp>
#include <frantic/math/utils.hpp>

#include <iomanip>

using namespace frantic::graphics;
using frantic::strings::tstring;

tstring length_unit::format_autoscale( double length, length_unit::option lengthUnit, int precisionDigits ) {
    // Normalize the Infinity and NaN cases, units aren't relevant for them.
    if( frantic::math::is_nan( length ) ) {
        return _T("nan");
    } else if( length == std::numeric_limits<double>::infinity() ) {
        return _T("inf");
    } else if( length == -std::numeric_limits<double>::infinity() ) {
        return _T("-inf");
    }

    // If no valid unit is provided, just format the number
    if( lengthUnit >= length_unit::invalid ) {
        std::stringstream ss;
        ss << std::setprecision( precisionDigits ) << length;
        return frantic::strings::to_tstring( ss.str() );
    }

    // Zero is also a special case, because of autoscaling
    if( length == 0 ) {
        std::stringstream ss;
        // We print 0.0 instead of `length` here to throw away -0
        ss << std::fixed << std::setprecision( precisionDigits - 1 ) << 0.0 << " ";
        return frantic::strings::to_tstring( ss.str() ) + length_unit::get_unit_data( lengthUnit ).shortName;
    }

    double absLength = fabs( length );

    // For values within a reasonable range, stick with the base unit instead of autoscaling it. This has the effect of
    // keeping 'yards' when that's the unit, even though we don't use 'yards' for the autoscaling bit.
    if( 0.1 < absLength && absLength < 1000.0 ) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision( precisionDigits - 1 ) << length << " ";
        return frantic::strings::to_tstring( ss.str() ) + length_unit::get_unit_data( lengthUnit ).shortName;
    }

    std::vector<length_unit::option> autoScales;

    if( lengthUnit == length_unit::inches || lengthUnit == length_unit::feet || lengthUnit == length_unit::yards ||
        lengthUnit == length_unit::miles ) {
        // It's an imperial unit, autoscale between inches/feet/miles
        autoScales.push_back( length_unit::inches );
        autoScales.push_back( length_unit::feet );
        autoScales.push_back( length_unit::miles );
    } else {
        // It's an SI unit, autoscale between mm/m/km
        autoScales.push_back( length_unit::millimeters );
        autoScales.push_back( length_unit::meters );
        autoScales.push_back( length_unit::kilometers );
    }

    // Pick the format unit that best fits the number
    length_unit::option formatUnit = autoScales.back();
    for( size_t i = 1; i < autoScales.size(); ++i ) {
        double asUnit = length * ( length_unit::get_unit_data( lengthUnit ).scaleToMicrometers /
                                   length_unit::get_unit_data( autoScales[i] ).scaleToMicrometers );
        if( fabs( asUnit ) < 1 ) {
            formatUnit = autoScales[i - 1];
            break;
        }
    }

    std::stringstream ss;
    double asUnit = length * ( length_unit::get_unit_data( lengthUnit ).scaleToMicrometers /
                               length_unit::get_unit_data( formatUnit ).scaleToMicrometers );
    ss << std::setprecision( precisionDigits ) << asUnit << " ";
    return frantic::strings::to_tstring( ss.str() ) + length_unit::get_unit_data( formatUnit ).shortName;
}
