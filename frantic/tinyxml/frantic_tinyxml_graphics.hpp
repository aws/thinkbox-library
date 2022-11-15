// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/animated_trimesh3_cache.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/motion_blurred_transform.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/tinyxml/frantic_tinyxml_utility.hpp>

#include <tinyxml2.h>

namespace frantic {
namespace tinyxml {

inline frantic::graphics::transform4f get_transform4f( const std::string& name,
                                                       const frantic::graphics::transform4f& defaultValue,
                                                       tinyxml2::XMLHandle handle, bool warnOnFail = true ) {

    tinyxml2::XMLHandle transformHandle = get_tinyxml_handle( handle, name );
    if( transformHandle.ToNode() == 0 ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        return defaultValue;
    }
    tinyxml2::XMLText* text = transformHandle.FirstChild().ToText();
    if( text != 0 ) {
        return frantic::graphics::transform4f::parse( text->Value() );
    } else {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't parse parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        // TODO: check whether there's some other form of matrix specified in the xml
        return defaultValue;
    }
}

inline frantic::graphics::transform4f get_transform4f( const std::string& name, tinyxml2::XMLHandle handle ) {
    tinyxml2::XMLHandle transformHandle = get_tinyxml_handle( handle, name );
    if( transformHandle.ToNode() == 0 )
        throw std::runtime_error(
            "XML variable " + get_xml_handle_path( handle, name ) +
            ", intended to be parsed as a transform matrix, could not be found in the .xml file" );
    tinyxml2::XMLText* text = transformHandle.FirstChild().ToText();
    if( text != 0 ) {
        return frantic::graphics::transform4f::parse( text->Value() );
    } else {
        // TODO: check whether there's some other form of matrix specified in the xml
        throw std::runtime_error(
            "XML variable " + get_xml_handle_path( handle, name ) +
            ", intended to be parsed as a transform matrix, could not be found in the .xml file" );
    }
}

inline void get_motion_blurred_transform( const std::string& name, tinyxml2::XMLHandle handle,
                                          frantic::graphics::motion_blurred_transform<float>& outXform ) {
    tinyxml2::XMLHandle transformHandle = get_tinyxml_handle( handle, name );

    frantic::graphics::transform4f xform = get_transform4f( "transform", transformHandle );
    std::vector<frantic::graphics::transform4f> xformArray;
    for( tinyxml2::XMLNode* node = transformHandle.FirstChildElement( "animatedTransform" ).FirstChildElement( "transform" ).ToNode(); node;
         node = node->NextSiblingElement( "transform" ) ) {
        xformArray.push_back( tinyxml::get_transform4f( ".", tinyxml2::XMLHandle( node ) ) );
    }

    outXform.set( xform, xformArray );
}

inline frantic::graphics::boundbox3f get_boundbox3f( const std::string& name,
                                                     const frantic::graphics::boundbox3f& defaultValue,
                                                     tinyxml2::XMLHandle handle, bool warnOnFail = true ) {
    tinyxml2::XMLHandle boundboxHandle = get_tinyxml_handle( handle, name );
    if( boundboxHandle.ToNode() == 0 ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        return defaultValue;
    }
    tinyxml2::XMLText* text = boundboxHandle.FirstChild().ToText();
    if( text != 0 ) {
        return frantic::graphics::boundbox3f::parse( text->Value() );
    } else {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        // TODO: check whether there's some other form of matrix specified in the xml
        return defaultValue;
    }
}

inline frantic::graphics::boundbox3f get_boundbox3f( const std::string& name, tinyxml2::XMLHandle handle ) {
    tinyxml2::XMLHandle transformHandle = get_tinyxml_handle( handle, name );
    if( transformHandle.ToNode() == 0 )
        throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a bounding box, could not be found in the .xml file" );
    tinyxml2::XMLText* text = transformHandle.FirstChild().ToText();
    if( text != 0 ) {
        return frantic::graphics::boundbox3f::parse( text->Value() );
    } else {
        // TODO: check whether there's some other form of matrix specified in the xml
        throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a bounding box, could not be found in the .xml file" );
    }
}

inline frantic::graphics2d::size2 get_size2( const std::string& name, const frantic::graphics2d::size2& defaultValue,
                                             tinyxml2::XMLHandle handle, bool warnOnFail = true ) {
    // return frantic::graphics::size2( tinyxml::get_xml_value<int>("width", xml),
    // tinyxml::get_xml_value<int>("height",xml) );

    tinyxml2::XMLHandle sizeHandle = get_tinyxml_handle( handle, name );
    if( sizeHandle.ToNode() == 0 || sizeHandle.FirstChildElement( "width" ).ToNode() == 0 ||
        sizeHandle.FirstChildElement( "height" ).ToNode() == 0 ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        return defaultValue;
    }

    try {
        return frantic::graphics2d::size2( tinyxml::get_xml_value<int>( "width", sizeHandle ),
                                           tinyxml::get_xml_value<int>( "height", sizeHandle ) );
    } catch( std::exception& e ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read parameter \"" << get_xml_handle_path( handle, name )
                      << "\" because of error \"" << e.what() << "\", assuming default of " << defaultValue
                      << std::endl;
        // TODO: check whether there's some other form of size2 specified in the xml
        return defaultValue;
    }
}

inline frantic::graphics2d::size2 get_size2( const std::string& name, tinyxml2::XMLHandle handle ) {
    // return frantic::graphics::size2( tinyxml::get_xml_value<int>("width", xml),
    // tinyxml::get_xml_value<int>("height",xml) );

    tinyxml2::XMLHandle sizeHandle = get_tinyxml_handle( handle, name );
    if( sizeHandle.ToNode() == 0 || sizeHandle.FirstChildElement( "width" ).ToNode() == 0 ||
        sizeHandle.FirstChildElement( "height" ).ToNode() == 0 )
        throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a size2, could not be found in the .xml file" );

    try {
        return frantic::graphics2d::size2( tinyxml::get_xml_value<int>( "width", sizeHandle ),
                                           tinyxml::get_xml_value<int>( "height", sizeHandle ) );
    } catch( std::exception& ) {
        throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a size2, was not valid" );
    }
}

inline frantic::graphics2d::vector2 get_vector2( const std::string& name,
                                                 const frantic::graphics2d::vector2& defaultValue, tinyxml2::XMLHandle handle,
                                                 bool warnOnFail = true ) {
    // return frantic::graphics::size2( tinyxml::get_xml_value<int>("width", xml),
    // tinyxml::get_xml_value<int>("height",xml) );

    tinyxml2::XMLHandle sizeHandle = get_tinyxml_handle( handle, name );
    if( sizeHandle.ToNode() == 0 || sizeHandle.FirstChildElement( "x" ).ToNode() == 0 ||
        sizeHandle.FirstChildElement( "x" ).ToNode() == 0 ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        return defaultValue;
    }

    try {
        return frantic::graphics2d::vector2( tinyxml::get_xml_value<int>( "x", sizeHandle ),
                                             tinyxml::get_xml_value<int>( "y", sizeHandle ) );
    } catch( std::exception& e ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't parse parameter \"" << get_xml_handle_path( handle, name )
                      << "\" because of error \"" << e.what() << "\", assuming default of " << defaultValue
                      << std::endl;
        // TODO: check whether there's some other form of vector2 specified in the xml
        return defaultValue;
    }
}

inline frantic::graphics::color3f get_color3f( const std::string& name, const frantic::graphics::color3f& defaultValue,
                                               tinyxml2::XMLHandle handle, bool warnOnFail = true ) {
    // return frantic::graphics::size2( tinyxml::get_xml_value<int>("width", xml),
    // tinyxml::get_xml_value<int>("height",xml) );

    tinyxml2::XMLHandle colorHandle = get_tinyxml_handle( handle, name );
    if( colorHandle.ToNode() == 0 ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        return defaultValue;
    }

    try {
        if( colorHandle.FirstChildElement( "r" ).ToNode() != 0 ) {
            return frantic::graphics::color3f( tinyxml::get_xml_value<float>( "r", colorHandle ),
                                               tinyxml::get_xml_value<float>( "g", colorHandle ),
                                               tinyxml::get_xml_value<float>( "b", colorHandle ) );
        } else if( colorHandle.FirstChildElement( "red" ).ToNode() != 0 ) {
            return frantic::graphics::color3f( tinyxml::get_xml_value<float>( "red", colorHandle ),
                                               tinyxml::get_xml_value<float>( "green", colorHandle ),
                                               tinyxml::get_xml_value<float>( "blue", colorHandle ) );
        } else {
            std::string colorString = tinyxml::get_xml_value<std::string>( "", colorHandle );
            return frantic::graphics::color3f::parse( frantic::strings::to_tstring( colorString ) );
        }
    } catch( std::exception& e ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't parse parameter \"" << get_xml_handle_path( handle, name )
                      << "\" because of error \"" << e.what() << "\", assuming default of " << defaultValue
                      << std::endl;
        // TODO: check whether there's some other form of color specified in the xml
        return defaultValue;
    }
    return defaultValue;
}

inline frantic::graphics::color3f get_color3f( const std::string& name, tinyxml2::XMLHandle handle ) {
    // return frantic::graphics::size2( tinyxml::get_xml_value<int>("width", xml),
    // tinyxml::get_xml_value<int>("height",xml) );

    tinyxml2::XMLHandle colorHandle = get_tinyxml_handle( handle, name );

    try {
        if( colorHandle.FirstChildElement( "r" ).ToNode() != 0 ) {
            return frantic::graphics::color3f( tinyxml::get_xml_value<float>( "r", colorHandle ),
                                               tinyxml::get_xml_value<float>( "g", colorHandle ),
                                               tinyxml::get_xml_value<float>( "b", colorHandle ) );
        } else if( colorHandle.FirstChildElement( "red" ).ToNode() != 0 ) {
            return frantic::graphics::color3f( tinyxml::get_xml_value<float>( "red", colorHandle ),
                                               tinyxml::get_xml_value<float>( "green", colorHandle ),
                                               tinyxml::get_xml_value<float>( "blue", colorHandle ) );
        } else {
            throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                      ", intended to be parsed as a color3f, could not be found in the .xml file" );
        }
    } catch( const std::exception& ) {
        throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a color3f, could not be found in the .xml file" );
    }
}

/**
 *  Attempt to get the value of attributes x, y, and z from the xml element
 * and store these values in the output vector.
 *
 * @param h the handle of an xml element to get the x, y and z attribute
 *		values from.
 * @param[out] outVector if the function returns true, then this vector
 *		holds the x, y and z attribute values in the given xml handle.
 * @return true if the outVector could be set using the attributes of the
 *		given xml handle.
 */
inline bool try_get_vector3f_from_attributes( tinyxml2::XMLHandle h, frantic::graphics::vector3f& outVector ) {
    const tinyxml2::XMLElement* e = h.ToElement();

    if( !e ) {
        return false;
    }

    int status;
    double d;
    const std::string axisLetters( "xyz" );

    for( int axis = 0; axis < 3; ++axis ) {
        // In the current tinyxml there is a QueryFloatAttribute, which
        // is that I want to use here.  In the version we have there is
        // instead an overloaded QueryDoubleAttribute for floats.
        status = e->QueryDoubleAttribute( axisLetters.substr( axis, 1 ).c_str(), &d );
        if( status != tinyxml2::XMLError::XML_SUCCESS )
            return false;
        outVector[axis] = static_cast<float>( d );
    }
    return true;
}

inline frantic::graphics::vector3f get_vector3f( const std::string& name,
                                                 const frantic::graphics::vector3f& defaultValue, tinyxml2::XMLHandle handle,
                                                 bool warnOnFail = true ) {

    tinyxml2::XMLHandle vectorHandle = get_tinyxml_handle( handle, name );
    try {
        frantic::graphics::vector3f out;
        if( try_get_vector3f_from_attributes( vectorHandle, out ) ) {
            return out;
        }
        if( vectorHandle.FirstChildElement( "x" ).ToNode() != 0 ) {
            return frantic::graphics::vector3f( tinyxml::get_xml_value<float>( "x", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "y", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "z", vectorHandle ) );
        } else if( vectorHandle.FirstChildElement( "X" ).ToNode() != 0 ) {
            return frantic::graphics::vector3f( tinyxml::get_xml_value<float>( "X", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "Y", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "Z", vectorHandle ) );
        } else {
            return frantic::graphics::vector3f::parse( tinyxml::get_string( name, handle ) );
        }
    } catch( std::exception& ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't parse parameter \"" << get_xml_handle_path( handle, name )
                      << "\", assuming default of " << defaultValue << std::endl;
        return defaultValue;
    }
}

inline frantic::graphics::vector3f get_vector3f( const std::string& name, tinyxml2::XMLHandle handle ) {

    tinyxml2::XMLHandle vectorHandle = get_tinyxml_handle( handle, name );
    try {
        frantic::graphics::vector3f out;
        if( try_get_vector3f_from_attributes( vectorHandle, out ) ) {
            return out;
        }
        if( vectorHandle.FirstChildElement( "x" ).ToNode() != 0 ) {
            return frantic::graphics::vector3f( tinyxml::get_xml_value<float>( "x", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "y", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "z", vectorHandle ) );
        } else if( vectorHandle.FirstChildElement( "X" ).ToNode() != 0 ) {
            return frantic::graphics::vector3f( tinyxml::get_xml_value<float>( "X", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "Y", vectorHandle ),
                                                tinyxml::get_xml_value<float>( "Z", vectorHandle ) );
        } else {
            return frantic::graphics::vector3f::parse( tinyxml::get_string( name, handle ) );
        }
    } catch( std::exception& ) {
        throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a vector3f, could not be found in the .xml file" );
    }
}

inline frantic::geometry::animated_trimesh3_cache get_animated_trimesh3_cache( const std::string& name,
                                                                               tinyxml2::XMLHandle handle ) {
    tinyxml2::XMLHandle cacheHandle = get_tinyxml_handle( handle, name );
    if( cacheHandle.ToNode() == 0 )
        throw std::runtime_error(
            "XML variable " + get_xml_handle_path( handle, name ) +
            ", intended to be parsed as an animated_trimesh3_cache, could not be found in the .xml file" );

    try {
        frantic::geometry::animated_trimesh3_cache result;
        result.set( frantic::strings::to_tstring( tinyxml::get_string( "filePrefix", cacheHandle ) ),
                    tinyxml::get_int( "framePadding", cacheHandle ),
                    frantic::strings::to_tstring( tinyxml::get_string( "filePostfix", cacheHandle ) ),
                    tinyxml::get_int( "countsPerFrame", cacheHandle ), tinyxml::get_int( "countOffset", cacheHandle ),
                    tinyxml::get_bool( "animated", cacheHandle ), tinyxml::get_float( "timeSpeedFactor", cacheHandle ),
                    tinyxml::get_int( "countIntervalStart", cacheHandle ),
                    tinyxml::get_int( "countIntervalEnd", cacheHandle ) );
        return result;
    } catch( std::exception& ) {
        throw std::runtime_error( "XML variable " + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a size2, could not be found in the .xml file" );
    }
}

// This adds a node called <name> to the handle, containing the data of meshCache
inline void add_animated_trimesh3_cache( const std::string& name, tinyxml2::XMLHandle handle,
                                         const frantic::geometry::animated_trimesh3_cache& meshCache ) {
    tinyxml2::XMLNode* handleNode = handle.ToNode();
    if( !handleNode )
        throw std::runtime_error( "add_animated_trimesh3_cache: tried to add data to a null xml node." );
    tinyxml2::XMLDocument* doc = handleNode->GetDocument();

    tinyxml2::XMLElement* node = doc->NewElement( name.c_str() );
    tinyxml2::XMLElement* filePrefixNode = doc->NewElement( "filePrefix" );
    tinyxml2::XMLElement* framePaddingNode = doc->NewElement( "framePadding" );
    tinyxml2::XMLElement* filePostfixNode = doc->NewElement( "filePostfix" );
    tinyxml2::XMLElement* countsPerFrameNode = doc->NewElement( "countsPerFrame" );
    tinyxml2::XMLElement* countOffsetNode = doc->NewElement( "countOffset" );
    tinyxml2::XMLElement* animatedNode = doc->NewElement( "animated" );
    tinyxml2::XMLElement* countIntervalStartNode = doc->NewElement( "countIntervalStart" );
    tinyxml2::XMLElement* countIntervalEndNode = doc->NewElement( "countIntervalEnd" );
    tinyxml2::XMLElement* timeSpeedFactorNode = doc->NewElement( "timeSpeedFactor" );

    filePrefixNode->InsertEndChild( doc->NewText( frantic::strings::to_string( meshCache.get_filePrefix() ).c_str() ) );
    framePaddingNode->InsertEndChild( doc->NewText( boost::lexical_cast<std::string>( meshCache.get_framePadding() ).c_str() ) );
    filePostfixNode->InsertEndChild( doc->NewText( frantic::strings::to_string( meshCache.get_filePostfix() ).c_str() ) );
    countsPerFrameNode->InsertEndChild(
        doc->NewText( boost::lexical_cast<std::string>( meshCache.get_countsPerFrame() ).c_str() ) );
    countOffsetNode->InsertEndChild( doc->NewText( boost::lexical_cast<std::string>( meshCache.get_countOffset() ).c_str() ) );
    animatedNode->InsertEndChild( doc->NewText( meshCache.get_animated() ? "true" : "false" ) );
    countIntervalStartNode->InsertEndChild(
        doc->NewText( boost::lexical_cast<std::string>( meshCache.get_countInterval().first ).c_str() ) );
    countIntervalEndNode->InsertEndChild(
        doc->NewText( boost::lexical_cast<std::string>( meshCache.get_countInterval().second ).c_str() ) );
    timeSpeedFactorNode->InsertEndChild(
        doc->NewText( boost::lexical_cast<std::string>( meshCache.get_timeSpeedFactor() ).c_str() ) );

    node->InsertEndChild( filePrefixNode );
    node->InsertEndChild( framePaddingNode );
    node->InsertEndChild( filePostfixNode );
    node->InsertEndChild( countsPerFrameNode );
    node->InsertEndChild( countOffsetNode );
    node->InsertEndChild( animatedNode );
    node->InsertEndChild( countIntervalStartNode );
    node->InsertEndChild( countIntervalEndNode );
    node->InsertEndChild( timeSpeedFactorNode );

    handle.ToNode()->InsertEndChild( node );
}

/**
 *  Set the x, y and z attributes of the given xml element from the values
 * in the given vector.
 *
 * @param e an xml element whose x, y, and z attributes should be set.
 * @param v a vector whose values will be used to set the x, y and z
 *		attributes.
 */
inline void set_attributes_from_vector3f( tinyxml2::XMLElement* e, const frantic::graphics::vector3f& v ) {
    if( !e ) {
        throw std::runtime_error( "set_attributes_from_vector3f Error: the specified XmlElement is NULL." );
    }
    e->SetAttribute( "x", boost::lexical_cast<std::string>( v.x ).c_str() );
    e->SetAttribute( "y", boost::lexical_cast<std::string>( v.y ).c_str() );
    e->SetAttribute( "z", boost::lexical_cast<std::string>( v.z ).c_str() );
}

} // namespace tinyxml
} // namespace frantic
