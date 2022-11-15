// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/xmesh_metadata.hpp>

#include <frantic/locale/locale.hpp>
#include <frantic/strings/utf8.hpp>
#include <frantic/tinyxml/frantic_tinyxml_graphics.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>

namespace frantic {
namespace geometry {

namespace {

void parse_frames_per_second( const std::string& s, boost::int64_t& outNumerator, boost::int64_t& outDenominator ) {
    outNumerator = 1;
    outDenominator = 1;

    std::vector<std::string> parts;
    boost::algorithm::split( parts, s, boost::algorithm::is_any_of( "/" ) );

    if( parts.size() > 2 ) {
        throw std::runtime_error( "parse_frames_per_second: Unexpected number of \'/\' separators in \'" + s +
                                  "\'.  Expected 0 or 1." );
    }

    boost::int64_t numerator = 1;
    boost::int64_t denominator = 1;

    if( parts.size() < 1 ) {
        throw std::runtime_error( "parse_frames_per_second: Internal error: Could not find string part in \'" + s +
                                  "\'." );
    }
    numerator = boost::lexical_cast<boost::int64_t>( parts[0] );

    if( parts.size() > 1 ) {
        denominator = boost::lexical_cast<boost::int64_t>( parts[1] );
    }

    outNumerator = numerator;
    outDenominator = denominator;
}

std::string get_frames_per_second_string( const boost::rational<boost::int64_t>& framesPerSecond ) {
    std::stringstream ss;
    ss << boost::lexical_cast<std::string>( framesPerSecond.numerator() );
    if( framesPerSecond.denominator() != 1 ) {
        ss << "/" << boost::lexical_cast<std::string>( framesPerSecond.denominator() );
    }
    return ss.str();
}

const char* get_transform_type_string( xmesh_metadata::transform_type_t transformType ) {
    switch( transformType ) {
    case xmesh_metadata::transform_type_none:
        return "none";
    case xmesh_metadata::transform_type_point:
        return "point";
    case xmesh_metadata::transform_type_vector:
        return "vector";
    case xmesh_metadata::transform_type_normal:
        return "normal";
    default:
        throw std::runtime_error( "get_transform_type_string Error: unknown transform type " +
                                  boost::lexical_cast<std::string>( transformType ) );
    }
}

xmesh_metadata::transform_type_t get_transform_type_from_string( const std::string& s ) {
    if( s.empty() || s == "none" ) {
        return xmesh_metadata::transform_type_none;
    } else if( s == "point" ) {
        return xmesh_metadata::transform_type_point;
    } else if( s == "vector" ) {
        return xmesh_metadata::transform_type_vector;
    } else if( s == "normal" ) {
        return xmesh_metadata::transform_type_normal;
    }
    return xmesh_metadata::transform_type_invalid;
}

xmesh_metadata::length_unit_t get_length_unit_from_string( const std::string& s ) {
    if( s == "" || s == "none" ) {
        return xmesh_metadata::length_unit_unitless;
    } else if( s == "in" ) {
        return xmesh_metadata::length_unit_inches;
    } else if( s == "ft" ) {
        return xmesh_metadata::length_unit_feet;
    } else if( s == "mi" ) {
        return xmesh_metadata::length_unit_miles;
    } else if( s == "mm" ) {
        return xmesh_metadata::length_unit_millimeters;
    } else if( s == "cm" ) {
        return xmesh_metadata::length_unit_centimeters;
    } else if( s == "m" ) {
        return xmesh_metadata::length_unit_meters;
    } else if( s == "km" ) {
        return xmesh_metadata::length_unit_kilometers;
    }
    return xmesh_metadata::length_unit_invalid;
}

std::string get_length_unit_string( xmesh_metadata::length_unit_t lengthUnit ) {
    switch( lengthUnit ) {
    case xmesh_metadata::length_unit_unitless:
        return "none";
    case xmesh_metadata::length_unit_inches:
        return "in";
    case xmesh_metadata::length_unit_feet:
        return "ft";
    case xmesh_metadata::length_unit_miles:
        return "mi";
    case xmesh_metadata::length_unit_millimeters:
        return "mm";
    case xmesh_metadata::length_unit_centimeters:
        return "cm";
    case xmesh_metadata::length_unit_meters:
        return "m";
    case xmesh_metadata::length_unit_kilometers:
        return "km";
    default:
        throw std::runtime_error( "get_length_unit_string: invalid length unit: " +
                                  boost::lexical_cast<std::string>( lengthUnit ) );
    }
}

bool not_alpha_char( const char& c ) { return !isalpha( c, std::locale::classic() ); }

void parse_length_unit( std::string s, double& outLength, frantic::geometry::xmesh_metadata::length_unit_t& outUnit ) {
    outLength = 1.0;
    outUnit = xmesh_metadata::length_unit_invalid;

    boost::algorithm::trim( s );

    std::string::iterator unitSubstringBegin( std::find_if( s.rbegin(), s.rend(), not_alpha_char ).base() );

    const std::string lengthString( s.begin(), unitSubstringBegin );
    const std::string unitString( unitSubstringBegin, s.end() );

    if( lengthString.size() > 0 ) {
        double d = boost::lexical_cast<double>( lengthString );
        outLength = d;
    }

    xmesh_metadata::length_unit_t lengthUnit = get_length_unit_from_string( unitString );
    if( lengthUnit == frantic::geometry::xmesh_metadata::length_unit_invalid ) {
        throw std::runtime_error( "parse_length_unit_string: Unknown unit string: \'" + unitString + "\' in string \'" +
                                  s + "\'." );
    } else {
        outUnit = lengthUnit;
    }
}

double get_to_meter_conversion_factor( xmesh_metadata::length_unit_t from ) {
    switch( from ) {
    case xmesh_metadata::length_unit_unitless:
        return 1.0;
    case xmesh_metadata::length_unit_inches:
        return 0.0254;
    case xmesh_metadata::length_unit_feet:
        return 12.0 * 0.0254;
    case xmesh_metadata::length_unit_miles:
        return 5280.0 * 12.0 * 0.0254;
    case xmesh_metadata::length_unit_millimeters:
        return 0.001;
    case xmesh_metadata::length_unit_centimeters:
        return 0.01;
    case xmesh_metadata::length_unit_meters:
        return 1.0;
    case xmesh_metadata::length_unit_kilometers:
        return 1000.0;
    default:
        throw std::runtime_error( "get_to_meter_conversion_factor: Unknown unit: " +
                                  boost::lexical_cast<std::string>( from ) );
    }
}

boost::rational<boost::int64_t> get_to_meter_conversion_fraction( xmesh_metadata::length_unit_t from ) {
    switch( from ) {
    case xmesh_metadata::length_unit_unitless:
        return boost::rational<boost::int64_t>( 1, 1 );
    case xmesh_metadata::length_unit_inches:
        return boost::rational<boost::int64_t>( 254, 10000 );
    case xmesh_metadata::length_unit_feet:
        return boost::rational<boost::int64_t>( 12 * 254, 10000 );
    case xmesh_metadata::length_unit_miles:
        return boost::rational<boost::int64_t>( 5280 * 12 * 254, 10000 );
    case xmesh_metadata::length_unit_millimeters:
        return boost::rational<boost::int64_t>( 1, 1000 );
    case xmesh_metadata::length_unit_centimeters:
        return boost::rational<boost::int64_t>( 1, 100 );
    case xmesh_metadata::length_unit_meters:
        return boost::rational<boost::int64_t>( 1, 1 );
    case xmesh_metadata::length_unit_kilometers:
        return boost::rational<boost::int64_t>( 1000, 1 );
    default:
        throw std::runtime_error( "get_to_meter_conversion_factor: Unknown unit: " +
                                  boost::lexical_cast<std::string>( from ) );
    }
}

double get_from_meter_conversion_factor( xmesh_metadata::length_unit_t to ) {
    return 1.0 / get_to_meter_conversion_factor( to );
}

boost::rational<boost::int64_t> get_from_meter_conversion_fraction( xmesh_metadata::length_unit_t to ) {
    return 1 / get_to_meter_conversion_fraction( to );
}

double get_length_unit_conversion_factor( xmesh_metadata::length_unit_t from, xmesh_metadata::length_unit_t to ) {
    return boost::rational_cast<double>( get_to_meter_conversion_fraction( from ) *
                                         get_from_meter_conversion_fraction( to ) );
}

// TODO: these are shared with xmesh_writer.cpp.  Move them somewhere in common.
inline std::string utf8_from_channel_name( const std::string& s ) {
    if( frantic::strings::is_valid_utf8( s ) ) {
        return s;
    } else {
        return frantic::strings::to_utf8( s );
    }
}

inline std::string utf8_from_channel_name( const std::wstring& s ) { return frantic::strings::to_utf8( s ); }

} // anonymous namespace

void insert_channel_transform_type_tags( const xmesh_metadata& metadata, tinyxml2::XMLHandle xmeshHandle,
                                         const std::string& elementName, bool enableCreateChannels ) {
    if( metadata.channelTransformType.size() == 0 ) {
        return;
    }
    tinyxml2::XMLNode* xmeshHandleNode = xmeshHandle.ToNode();
    if (!xmeshHandleNode) {
        return;
    }
    tinyxml2::XMLDocument* doc = xmeshHandleNode->GetDocument();

    // This function assumes that the xmesh file is saved using UTF-8 string encoding.

    // "untagged" channels are channels that have transforms but are not present in the
    // xmesh file.  If enableCreateChannels is true, we create a channel for all such
    // untagged channels.  This is meant for internal use: currently it's used to transport
    // metadata onto Deadline slaves.
    std::set<std::string> untaggedChannelsUTF8;
    std::map<std::string, frantic::tstring> utf8ToChannelName;
    for( xmesh_metadata::channel_transform_type_collection_t::const_iterator i = metadata.channelTransformType.begin();
         i != metadata.channelTransformType.end(); ++i ) {
        const std::string channelNameUTF8 = utf8_from_channel_name( i->first );
        if( !frantic::strings::is_valid_utf8( channelNameUTF8 ) ) {
            throw std::runtime_error( "Internal Error: channel name \'" + frantic::strings::to_string( i->first ) +
                                      "\' was not converted to valid UTF-8" );
        }
        utf8ToChannelName[channelNameUTF8] = i->first;
        untaggedChannelsUTF8.insert( channelNameUTF8 );
    }

    tinyxml2::XMLElement* channelElement = xmeshHandle.FirstChildElement( elementName.c_str() ).ToElement();
    while( channelElement ) {
        tinyxml2::XMLElement* nameElement = channelElement->FirstChildElement( "name" );
        if( nameElement ) {
            tinyxml2::XMLHandle nameHandle( nameElement );
            tinyxml2::XMLText* nameText = nameHandle.FirstChild().ToText();
            if( nameText ) {
                const std::string channelNameUTF8( nameText->Value() );
                std::map<std::string, frantic::tstring>::const_iterator i = utf8ToChannelName.find( channelNameUTF8 );
                if( i != utf8ToChannelName.end() ) {
                    const frantic::tstring& channelName = i->second;
                    xmesh_metadata::channel_transform_type_collection_t::const_iterator j =
                        metadata.channelTransformType.find( channelName );
                    if( j != metadata.channelTransformType.end() ) {
                        untaggedChannelsUTF8.erase( channelNameUTF8 );
                        if( j->second != xmesh_metadata::transform_type_none &&
                            j->second != xmesh_metadata::transform_type_invalid ) {
                            tinyxml2::XMLElement* dataClassElement = doc->NewElement( "transformType" );
                            dataClassElement->LinkEndChild( doc->NewText( get_transform_type_string( j->second ) ) );
                            channelElement->LinkEndChild( dataClassElement );
                        }
                    }
                }
            }
        }
        channelElement = channelElement->NextSiblingElement( elementName.c_str() );
    }
    if( enableCreateChannels ) {
        tinyxml2::XMLElement* pXMesh = xmeshHandle.ToElement();
        if( pXMesh ) {
            BOOST_FOREACH( const std::string& channelNameUTF8, untaggedChannelsUTF8 ) {
                xmesh_metadata::transform_type_t transformType =
                    metadata.get_channel_transform_type( utf8ToChannelName[channelNameUTF8] );
                if( transformType != xmesh_metadata::transform_type_none &&
                    transformType != xmesh_metadata::transform_type_invalid ) {
                    tinyxml2::XMLElement* channelElement = doc->NewElement( elementName.c_str() );
                    pXMesh->LinkEndChild( channelElement );

                    tinyxml2::XMLElement* nameElement = doc->NewElement( "name" );
                    nameElement->LinkEndChild( doc->NewText( channelNameUTF8.c_str() ) );
                    channelElement->LinkEndChild( nameElement );

                    tinyxml2::XMLElement* dataClassElement = doc->NewElement( "transformType" );
                    dataClassElement->LinkEndChild( doc->NewText( get_transform_type_string( transformType ) ) );
                    channelElement->LinkEndChild( dataClassElement );
                }
            }
        }
    }
}

xmesh_metadata::xmesh_metadata() { clear(); }

// xmesh_metadata::ptr_type xmesh_metadata::create() {
//	xmesh_metadata::ptr_type result( new xmesh_metadata() );
//	return result;
// }

// xmesh_metadata::ptr_type xmesh_metadata::from_xml( tinyxml2::XMLDocument& doc, const frantic::tstring& pathForErrorMessage )
// { 	xmesh_metadata::ptr_type result( new xmesh_metadata() ); 	result->from_xml( doc, pathForErrorMessage );
// return
// result;
// }

void xmesh_metadata::clear( void ) {
    clear_length_unit();
    clear_frames_per_second();
    clear_boundbox();
    channelTransformType.clear();
    userData.clear();
}

// Frames per second
void xmesh_metadata::clear_frames_per_second() {
    hasFramesPerSecond = false;
    framesPerSecond.assign( 1, 1 );
}

bool xmesh_metadata::has_frames_per_second() const { return hasFramesPerSecond; }

void xmesh_metadata::set_frames_per_second( boost::int64_t numerator, boost::int64_t denominator ) {
    framesPerSecond.assign( numerator, denominator );
    hasFramesPerSecond = true;
}

void xmesh_metadata::set_frames_per_second( const boost::rational<boost::int64_t>& framesPerSecond ) {
    this->framesPerSecond = framesPerSecond;
    hasFramesPerSecond = true;
}

boost::rational<boost::int64_t> xmesh_metadata::get_frames_per_second() const { return framesPerSecond; }

// Boundbox
void xmesh_metadata::clear_boundbox() {
    hasBoundbox = false;
    boundbox.set_to_empty();
}

void xmesh_metadata::set_boundbox( const frantic::graphics::boundbox3f& boundbox ) {
    this->boundbox = boundbox;
    hasBoundbox = true;
}

bool xmesh_metadata::has_boundbox() const { return hasBoundbox; }

frantic::graphics::boundbox3f xmesh_metadata::get_boundbox() const { return boundbox; }

// Length unit
void xmesh_metadata::clear_length_unit() {
    hasLengthUnit = false;
    lengthUnitScale = 1.0;
    lengthUnit = length_unit_unitless;
}
void xmesh_metadata::set_length_unit( double lengthUnitScale, xmesh_metadata::length_unit_t lengthUnit ) {
    this->lengthUnitScale = lengthUnitScale;
    this->lengthUnit = lengthUnit;
    hasLengthUnit = true;
}
bool xmesh_metadata::has_length_unit() const { return hasLengthUnit; }
double xmesh_metadata::get_length_unit_scale() const { return lengthUnitScale; }
xmesh_metadata::length_unit_t xmesh_metadata::get_length_unit() const { return lengthUnit; }

// Channel transform types
void xmesh_metadata::clear_channel_transform_type() { channelTransformType.clear(); }

void xmesh_metadata::set_channel_transform_type( const frantic::tstring& channelName, transform_type_t transformType ) {
    if( transformType == transform_type_invalid ) {
        throw std::runtime_error( "set_channel_transform_type Error: cannot set channel \'" +
                                  frantic::strings::to_string( channelName ) + "\' to invalid transform type." );
    }
    if( transformType == transform_type_none ) {
        channel_transform_type_collection_t::iterator i = channelTransformType.find( channelName );
        if( i != channelTransformType.end() ) {
            channelTransformType.erase( i );
        }
    } else {
        channelTransformType[channelName] = transformType;
    }
}
xmesh_metadata::transform_type_t
xmesh_metadata::get_channel_transform_type( const frantic::tstring& channelName ) const {
    channel_transform_type_collection_t::const_iterator i = channelTransformType.find( channelName );
    if( i == channelTransformType.end() ) {
        return xmesh_metadata::transform_type_none;
    } else {
        return i->second;
    }
}

// User data
bool xmesh_metadata::has_user_data( const frantic::tstring& key ) const {
    user_data_collection_t::const_iterator i = userData.find( key );
    return i != userData.end();
}
frantic::tstring xmesh_metadata::get_user_data( const frantic::tstring& key ) const {
    user_data_collection_t::const_iterator i = userData.find( key );
    if( i == userData.end() ) {
        throw std::runtime_error( "No user data with key: " + frantic::strings::to_string( key ) );
    } else {
        return i->second;
    }
}
void xmesh_metadata::set_user_data( const frantic::tstring& key, const frantic::tstring& value ) {
    userData[key] = value;
}
void xmesh_metadata::clear_user_data() { userData.clear(); }
void xmesh_metadata::erase_user_data( const frantic::tstring& key ) {
    user_data_collection_t::iterator i = userData.find( key );
    if( i != userData.end() ) {
        userData.erase( i );
    }
}
void xmesh_metadata::get_user_data_keys( std::vector<frantic::tstring>& out ) const {
    out.clear();
    for( user_data_collection_t::const_iterator i = userData.begin(); i != userData.end(); ++i ) {
        out.push_back( i->first );
    }
}

namespace {
frantic::tstring convert_utf8_string_to_user_data( const char* s ) {
    if( s ) {
#ifdef FRANTIC_USE_WCHAR
        return frantic::strings::wstring_from_utf8( s );
#else
        return frantic::strings::to_string( frantic::strings::wstring_from_utf8( s ) );
#endif
    } else {
        return frantic::tstring();
    }
}

frantic::tstring convert_legacy_string_to_user_data( const char* s ) {
    if( s ) {
        return frantic::strings::to_tstring( s );
    } else {
        return frantic::tstring();
    }
}

// Convert from XML string to channel name
// Copied from xmesh_writer.cpp
// TODO: These functions should be moved somewhere shared with xmesh_writer.cpp

// TODO: how are std::string channel names encoded?  I think they
// should be UTF-8, but for now we're just passing through the string
// as it was found in the XML file.
frantic::tstring convert_utf8_string_to_channel_name( const char* s ) {
    if( s ) {
#ifdef FRANTIC_USE_WCHAR
        return frantic::strings::wstring_from_utf8( s );
#else
        return s;
#endif
    } else {
        return frantic::tstring();
    }
}

frantic::tstring convert_legacy_string_to_channel_name( const char* s ) {
    if( s ) {
#ifdef FRANTIC_USE_WCHAR
        return frantic::strings::to_wstring( s );
#else
        return s;
#endif
    } else {
        return frantic::tstring();
    }
}

bool is_valid_legacy_string( const char* /*s*/ ) { return true; }

} // anonymous namespace

void read_xmesh_metadata( tinyxml2::XMLDocument& xmeshDocument, xmesh_metadata& outMetadata ) {
    const char* pPathForErrorMessage = xmeshDocument.Value();
    const std::string pathForErrorMessage( pPathForErrorMessage ? pPathForErrorMessage : "" );
    tinyxml2::XMLElement* pXMesh = xmeshDocument.RootElement();
    tinyxml2::XMLHandle xmeshHandle( pXMesh );
    outMetadata.clear();
    // units
    // hasLengthUnit = false;
    tinyxml2::XMLElement* units = xmeshHandle.FirstChildElement( "units" ).ToElement();
    if( units ) {
        tinyxml2::XMLHandle unitsHandle( units );
        tinyxml2::XMLElement* lengthUnitElement = unitsHandle.FirstChildElement( "length" ).ToElement();
        if( lengthUnitElement ) {
            tinyxml2::XMLHandle lengthUnitHandle( lengthUnitElement );
            tinyxml2::XMLText* pText = lengthUnitHandle.FirstChild().ToText();
            if( pText ) {
                const std::string s = pText->Value();
                double length;
                xmesh_metadata::length_unit_t unit;
                parse_length_unit( s, length, unit );
                if( unit == xmesh_metadata::length_unit_invalid ) {
                    throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage +
                                              "\" has an invalid value in <units><length>: " + s );
                } else {
                    outMetadata.set_length_unit( length, unit );
                }
            } else {
                throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage +
                                          "\" is missing a value in <units><length>." );
            }
        }
    }

    // frames per second
    tinyxml2::XMLElement* framesPerSecondElement = xmeshHandle.FirstChildElement( "framesPerSecond" ).ToElement();
    if( framesPerSecondElement ) {
        tinyxml2::XMLHandle framesPerSecondHandle( framesPerSecondElement );
        tinyxml2::XMLText* pText = framesPerSecondHandle.FirstChild().ToText();
        if( pText ) {
            const std::string s = pText->Value();
            boost::int64_t numerator, denominator;
            parse_frames_per_second( s, numerator, denominator );
            outMetadata.set_frames_per_second( numerator, denominator );
        } else {
            throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage +
                                      "\" is missing a value in <framesPerSecond>." );
        }
    }

    // boundbox
    tinyxml2::XMLHandle boundboxHandle = xmeshHandle.FirstChildElement( "boundbox" );
    if( boundboxHandle.ToElement() ) {
        const frantic::graphics::vector3f minCorner = frantic::tinyxml::get_vector3f( "minimum", boundboxHandle );
        const frantic::graphics::vector3f maxCorner = frantic::tinyxml::get_vector3f( "maximum", boundboxHandle );
        outMetadata.set_boundbox( frantic::graphics::boundbox3f( minCorner, maxCorner ) );
    }

    bool useUTF8 = false;
    const std::string encoding =
        boost::algorithm::to_upper_copy( frantic::tinyxml::get_encoding_declaration( xmeshDocument ) );
    if( encoding == "UTF-8" ) {
        useUTF8 = true;
    }

    frantic::tstring ( *to_internal_user_data )( const char* ) =
        useUTF8 ? convert_utf8_string_to_user_data : convert_legacy_string_to_user_data;
    frantic::tstring ( *to_internal_channel_name )( const char* ) =
        useUTF8 ? convert_utf8_string_to_channel_name : convert_legacy_string_to_channel_name;
    bool ( *is_valid_external_string )( const char* ) =
        useUTF8 ? static_cast<bool ( * )( const char* )>( &frantic::strings::is_valid_utf8 ) : is_valid_legacy_string;

    // userData
    { // scope for processing user data
        std::set<std::string> inputNames;

        std::size_t curUserData = 0;
        tinyxml2::XMLElement* pUserData = xmeshHandle.FirstChildElement( "userData" ).ToElement();
        while( pUserData ) {
            tinyxml2::XMLHandle dataHandle( pUserData );

            tinyxml2::XMLText* pDataName = dataHandle.FirstChildElement( "name" ).FirstChild().ToText();
            if( !pDataName ) {
                throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage + "\" in userData " +
                                          boost::lexical_cast<std::string>( curUserData ) +
                                          " has a missing or invalid <name> tag." );
            }

            const char* inputName = pDataName->Value();
            if( !is_valid_external_string( inputName ) ) {
                throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage + "\" in userData " +
                                          boost::lexical_cast<std::string>( curUserData ) +
                                          " has an invalid <name> string: \"" +
                                          ( inputName ? std::string( inputName ) : std::string() ) + "\"" );
            }
            const frantic::tstring name = to_internal_user_data( inputName );
            // Using separate checks for previous inputName and name,
            // because different inputNames may map to the same name
            // depending on the user's locale.
            if( inputNames.count( inputName ) == 0 ) {
                if( outMetadata.has_user_data( name ) ) {
                    // TODO: warn?
                }
            } else {
                throw std::runtime_error( "The userData \"" + frantic::strings::to_string( name ) +
                                          "\" was defined more than once in \"" + pathForErrorMessage + "\"." );
            }

            tinyxml2::XMLElement* pValue = dataHandle.FirstChildElement( "value" ).ToElement();
            if( pValue ) {
                tinyxml2::XMLHandle valueHandle( pValue );
                tinyxml2::XMLText* pValueText = valueHandle.FirstChild().ToText();
                if( pValueText ) {
                    const char* inputValue = pValueText->Value();
                    if( !is_valid_external_string( inputValue ) ) {
                        throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage + "\" in userData " +
                                                  boost::lexical_cast<std::string>( curUserData ) +
                                                  " has an invalid <value> string: \"" +
                                                  ( inputValue ? std::string( inputValue ) : std::string() ) + "\"" );
                    }
                    const frantic::tstring s = to_internal_user_data( inputValue );
                    outMetadata.set_user_data( name, s );
                } else {
                    throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage + "\" in userData \"" +
                                              frantic::strings::to_string( name ) + "\" is missing a value." );
                }
            } else {
                throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage + "\" in userData \"" +
                                          frantic::strings::to_string( name ) + "\" is missing a <value> tag." );
            }

            pUserData = pUserData->NextSiblingElement( "userData" );
            ++curUserData;
        }
    }

    // channel transform type
    const std::vector<std::string> channelElementNames =
        boost::assign::list_of( "vertexChannel" )( "faceChannel" )( "channel" );
    BOOST_FOREACH( const std::string& elementName, channelElementNames ) {
        tinyxml2::XMLElement* channelElement = xmeshHandle.FirstChildElement( elementName.c_str() ).ToElement();
        while( channelElement ) {
            frantic::tstring channelName;
            xmesh_metadata::transform_type_t transformType = xmesh_metadata::transform_type_invalid;

            tinyxml2::XMLElement* nameElement = channelElement->FirstChildElement( "name" );
            if( nameElement ) {
                tinyxml2::XMLHandle nameHandle( nameElement );
                tinyxml2::XMLText* nameText = nameHandle.FirstChild().ToText();
                if( nameText ) {
                    const char* inputChannelName( nameText->Value() );
                    if( !is_valid_external_string( inputChannelName ) ) {
                        throw std::runtime_error( "The xmesh file \"" + pathForErrorMessage +
                                                  "\" has an invalid channel name string: \"" + inputChannelName +
                                                  "\"" );
                    }
                    channelName = to_internal_channel_name( inputChannelName );
                }
            }
            tinyxml2::XMLElement* transformTypeElement = channelElement->FirstChildElement( "transformType" );
            if( transformTypeElement ) {
                tinyxml2::XMLHandle transformTypeHandle( transformTypeElement );
                tinyxml2::XMLText* transformTypeText = transformTypeHandle.FirstChild().ToText();
                if( transformTypeText ) {
                    transformType = get_transform_type_from_string( transformTypeText->Value() );
                }
            }
            if( !channelName.empty() && transformType != xmesh_metadata::transform_type_invalid &&
                transformType != xmesh_metadata::transform_type_none ) {
                outMetadata.set_channel_transform_type( channelName, transformType );
            }
            channelElement = channelElement->NextSiblingElement( elementName.c_str() );
        }
    }
}

void read_xmesh_metadata( const boost::filesystem::path& path, xmesh_metadata& outMetadata ) {
    outMetadata.clear();

    frantic::files::file_ptr f( frantic::files::fopen( path, "rb" ) );
    tinyxml2::XMLDocument doc( path.c_str() );
    tinyxml2::XMLError result = doc.LoadFile( f );
    if( result != tinyxml2::XMLError::XML_SUCCESS ) {
        if( doc.ErrorID() == tinyxml2::XMLError::XML_ERROR_FILE_COULD_NOT_BE_OPENED ) {
            throw frantic::files::file_open_error( "Failed to open xmesh file \"" + path.string() + "\"." );
        } else {
            throw std::runtime_error( "Failed to load or parse xmesh file \"" + path.string() + "\"." );
        }
    }

    tinyxml2::XMLElement* pXMesh = doc.RootElement();
    if( !pXMesh )
        throw std::runtime_error( "The xmesh file \"" + path.string() + "\" did not have a root element." );

    tinyxml2::XMLHandle xmeshHandle( pXMesh );

    read_xmesh_metadata( doc, outMetadata );
}

void write_xmesh_metadata( tinyxml2::XMLDocument& xmeshDocument, const xmesh_metadata& metadata, bool standalone ) {
    tinyxml2::XMLElement* pXMesh = xmeshDocument.RootElement();
    tinyxml2::XMLHandle xmeshHandle( pXMesh );
    if( pXMesh ) {
        tinyxml2::XMLDocument* doc = pXMesh->GetDocument();
        if( metadata.has_frames_per_second() ) {
            tinyxml2::XMLElement* framesPerSecondElement = doc->NewElement( "framesPerSecond" );
            pXMesh->LinkEndChild( framesPerSecondElement );

            framesPerSecondElement->LinkEndChild(
                doc->NewText( get_frames_per_second_string( metadata.get_frames_per_second() ).c_str() ) );
        }
        if( metadata.has_length_unit() &&
            ( metadata.get_length_unit_scale() != 1.0 ||
              metadata.get_length_unit() != xmesh_metadata::length_unit_unitless ) &&
            metadata.get_length_unit() != xmesh_metadata::length_unit_invalid ) {
            tinyxml2::XMLElement* units = doc->NewElement( "units" );
            pXMesh->LinkEndChild( units );

            tinyxml2::XMLElement* unitsLength = doc->NewElement( "length" );
            units->LinkEndChild( unitsLength );

            std::stringstream ss;
            if( metadata.get_length_unit_scale() != 1.0 ) {
                ss << frantic::locale::to_string_c( metadata.get_length_unit_scale() );
            }
            if( metadata.get_length_unit() != xmesh_metadata::length_unit_unitless ) {
                ss << get_length_unit_string( metadata.get_length_unit() );
            }
            unitsLength->LinkEndChild( doc->NewText( ss.str().c_str() ) );
        }
        if( metadata.has_boundbox() ) {
            tinyxml2::XMLHandle boundboxHandle = xmeshHandle.FirstChildElement( "boundbox" );
            if( !boundboxHandle.ToElement() ) {
                const frantic::graphics::boundbox3f boundbox = metadata.get_boundbox();
                tinyxml2::XMLElement* boundboxElement = doc->NewElement( "boundbox" );
                pXMesh->LinkEndChild( boundboxElement );
                tinyxml2::XMLElement* boundboxMin = doc->NewElement( "minimum" );
                boundboxElement->LinkEndChild( boundboxMin );
                boundboxMin->LinkEndChild( doc->NewText( boundbox.minimum().str().c_str() ) );
                tinyxml2::XMLElement* boundboxMax = doc->NewElement( "maximum" );
                boundboxElement->LinkEndChild( boundboxMax );
                boundboxMax->LinkEndChild( doc->NewText( boundbox.maximum().str().c_str() ) );
            }
        }
        std::vector<frantic::tstring> userDataKeys;
        metadata.get_user_data_keys( userDataKeys );
        BOOST_FOREACH( const frantic::tstring& key, userDataKeys ) {
            tinyxml2::XMLElement* userDataElement = doc->NewElement( "userData" );
            pXMesh->LinkEndChild( userDataElement );

            tinyxml2::XMLElement* userDataName = doc->NewElement( "name" );
            const std::string keyUTF8 = frantic::strings::to_utf8( key );
            if( !frantic::strings::is_valid_utf8( keyUTF8 ) ) {
                throw std::runtime_error( "write_xmesh_metadata Internal Error: user data name \'" +
                                          frantic::strings::to_string( key ) + "\' was not converted to valid UTF-8" );
            }
            userDataName->LinkEndChild( doc->NewText( keyUTF8.c_str() ) );
            userDataElement->LinkEndChild( userDataName );

            tinyxml2::XMLElement* userDataValue = doc->NewElement( "value" );
            const frantic::tstring value = metadata.get_user_data( key );
            const std::string valueUTF8 = frantic::strings::to_utf8( value );
            if( !frantic::strings::is_valid_utf8( valueUTF8 ) ) {
                throw std::runtime_error( "write_xmesh_metadata Internal Error: user data value \'" +
                                          frantic::strings::to_string( value ) +
                                          "\' was not converted to valid UTF-8" );
            }
            userDataElement->LinkEndChild( userDataValue );
            userDataValue->LinkEndChild( doc->NewText( valueUTF8.c_str() ) );
        }
        if( standalone ) {
            insert_channel_transform_type_tags( metadata, xmeshHandle, "channel", standalone );
        } else {
            insert_channel_transform_type_tags( metadata, xmeshHandle, "vertexChannel", standalone );
            insert_channel_transform_type_tags( metadata, xmeshHandle, "faceChannel", standalone );
        }
    }
}

void write_xmesh_metadata( const boost::filesystem::path& path, const xmesh_metadata& metadata ) {
    tinyxml2::XMLDocument doc( path.c_str() );

    doc.LinkEndChild( doc.NewDeclaration() );

    tinyxml2::XMLElement* pXMesh = doc.LinkEndChild( doc.NewElement( "xmesh" ) )->ToElement();

    tinyxml2::XMLElement* pVersion = pXMesh->LinkEndChild( doc.NewElement( "version" ) )->ToElement();
    pVersion->LinkEndChild( doc.NewText( "1" ) );

    write_xmesh_metadata( doc, metadata, true );

    frantic::files::file_ptr f( frantic::files::fopen( path, "w" ) );
    tinyxml2::XMLError result = doc.SaveFile( f );
    if( result != tinyxml2::XMLError::XML_SUCCESS )
        throw std::runtime_error( "write_xmesh_metadata: Error writing XML document to \"" + path.string() + "\"" );
}

std::pair<double, xmesh_metadata::length_unit_t> get_xmesh_length_unit_from_meters( double metersScale ) {
    double precision = 1e-5; // arbitrary fractional error that we accept
    boost::array<xmesh_metadata::length_unit_t, 7> unitCase = {
        xmesh_metadata::length_unit_inches,      xmesh_metadata::length_unit_feet,
        xmesh_metadata::length_unit_miles,       xmesh_metadata::length_unit_millimeters,
        xmesh_metadata::length_unit_centimeters, xmesh_metadata::length_unit_meters,
        xmesh_metadata::length_unit_kilometers };
    // check if it is any of the units
    for( size_t i = 0; i < unitCase.size(); ++i ) {
        if( std::abs( ( get_to_meter_conversion_factor( unitCase[i] ) - metersScale ) / metersScale ) < precision ) {
            return std::pair<double, xmesh_metadata::length_unit_t>( 1, unitCase[i] );
        }
    }
    return std::pair<double, xmesh_metadata::length_unit_t>( metersScale, xmesh_metadata::length_unit_meters );
}

double get_meters_from_xmesh_length_unit( double scale, frantic::geometry::xmesh_metadata::length_unit_t unit ) {
    return scale * get_to_meter_conversion_factor( unit );
}
} // namespace geometry
} // namespace frantic
