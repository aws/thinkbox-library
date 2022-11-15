// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/string_functions.hpp>

#include <tinyxml2.h>

#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include <fstream>
#include <sstream>
#include <typeinfo>

// File: frantic_tinyxml_utility.hpp
//
// A header file containing utility functions when using TinyXml.  It
// is typical to get an integer or float value specified in an XML file.
//
// This further provides the ability to have command line parameters of equivalent names override
// the xml if desired.
namespace frantic {
namespace tinyxml {

inline std::string get_xml_handle_path( tinyxml2::XMLHandle handle ) {
    if( handle.ToNode() == 0 )
        return "<nul>";
    else if( handle.ToNode()->ToDocument() != 0 )
        return "/";
    else if( handle.ToNode()->ToElement() != 0 ) {
        std::string root = get_xml_handle_path( tinyxml2::XMLHandle( handle.ToNode()->Parent() ) );
        if( root.size() == 0 || root[root.size() - 1] != '/' )
            root += "/";
        return root + handle.ToNode()->Value();
    } else
        return "unknown";
}

inline std::string get_xml_handle_path( tinyxml2::XMLHandle handle, const std::string& nodeName ) {
    if( handle.ToNode() == 0 )
        return nodeName == "." ? "<nul>" : "<nul>/" + nodeName;
    else if( handle.ToNode()->ToDocument() != 0 ) {
        if( nodeName == "." || nodeName.size() == 0 )
            return "/";
        if( nodeName[0] == '/' )
            return nodeName;
        else
            return "/" + nodeName;
    } else if( handle.ToNode()->ToElement() != 0 ) {
        std::string root = get_xml_handle_path( tinyxml2::XMLHandle( handle.ToNode()->Parent() ) );
        if( root.size() == 0 || root[root.size() - 1] != '/' )
            root += "/";
        return root + handle.ToNode()->Value() + ( nodeName == "." ? "" : ( "/" + nodeName ) );
    } else
        return "unknown" + ( nodeName == "." ? "" : ( "/" + nodeName ) );
}

namespace detail {
inline tinyxml2::XMLHandle get_tinyxml_handle_helper( tinyxml2::XMLHandle handle, const std::vector<std::string>& nodePath,
                                              int start ) {
    if( start >= static_cast<int>( nodePath.size() ) )
        return tinyxml2::XMLHandle( nullptr );

    tinyxml2::XMLNode* childNode = nullptr;
    for( childNode = handle.FirstChildElement( nodePath[start].c_str() ).ToNode(); childNode;
         childNode = childNode->NextSiblingElement( nodePath[start].c_str() ) ) {
        if( start == (int)nodePath.size() - 1 )
            return tinyxml2::XMLHandle( childNode );
        tinyxml2::XMLHandle foundNode = get_tinyxml_handle_helper( tinyxml2::XMLHandle( childNode ), nodePath, start + 1 );
        if( foundNode.ToNode() != nullptr ) {
            return foundNode;
        }
    }

    return tinyxml2::XMLHandle( nullptr );
}

} // namespace detail

inline tinyxml2::XMLHandle get_tinyxml_handle( tinyxml2::XMLHandle handle, const std::string& xmlPath ) {
    if( xmlPath == "" || xmlPath == "." || handle.ToNode() == 0 )
        return handle;

    std::vector<std::string> nodePath;
    frantic::strings::split( xmlPath, nodePath, '/' );

    return detail::get_tinyxml_handle_helper( handle, nodePath, 0 );
}

inline bool node_exists( tinyxml2::XMLHandle handle, const std::string& xmlPath ) {
    return get_tinyxml_handle( handle, xmlPath ).ToNode() != 0;
}

//
// Function: get_xml_value
//
// Gets a value from XML.
//
// Parameters:
//	name - the name of the field in XML that contains the float.
//  defaultValue - the value this function returns if "name" is not specified in XML.
//  handle - a tinyxml2::XMLHandle to the node that contains the "name" field.
//  warnOnFail - a flag indicating that a warning should be made if the value of "name" could not be cast to a float.
//
// Returns:
//  The float specified by "name" or "defaultValue" if "name" was not found.
//
template <class T>
inline T get_xml_value( const std::string& name, const T& defaultValue, tinyxml2::XMLHandle handle, bool warnOnFail = true,
                        boost::program_options::variables_map* cmd_line_vars = 0 ) {
    // Check whether the value was provided at the command line
    if( cmd_line_vars != 0 ) {
        if( cmd_line_vars->count( name ) > 0 ) {
            // my copy of gcc chokes on the commented-out line below.
            // The following two lines should be equivalent -PTF
            const boost::any& val = ( *cmd_line_vars )[name].value();
            return boost::any_cast<T>( val );
            // return (*cmd_line_vars)[name].as<T>();
        }
    }

    // Now check whether the value is provided in the XML file
    tinyxml2::XMLHandle childHandle = get_tinyxml_handle( handle, name );
    tinyxml2::XMLText* text = childHandle.FirstChild().ToText();
    if( text ) {
        try {
            return boost::lexical_cast<T>( text->Value() );
        } catch( const boost::bad_lexical_cast& ) {
            throw std::runtime_error( "XML variable \"" + get_xml_handle_path( handle, name ) + "\", with value \"" +
                                      text->Value() + "\" could not be parsed as a " + typeid( T ).name() );
        }
    } else if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() ) {
        std::cerr << "Warning: Couldn't read xml parameter \"" << get_xml_handle_path( handle, name )
                  << "\", assuming default of " << defaultValue << std::endl;
    }

    // Otherwise return the default value
    return defaultValue;
}

// TODO: The variables map thing here hasn't ended up being useful, we should excise this from all the code everywhere.
template <class T>
inline T get_xml_value( const std::string& name, tinyxml2::XMLHandle handle,
                        boost::program_options::variables_map* cmd_line_vars = 0 ) {
    // Check whether the value was provided at the command line
    if( cmd_line_vars != 0 ) {
        if( cmd_line_vars->count( name ) > 0 ) {
            // my copy of gcc chokes on the commented-out line below.
            // The following two lines should be equivalent -PTF
            const boost::any& val = ( *cmd_line_vars )[name].value();
            return boost::any_cast<T>( val );
            //
            // return (*cmd_line_vars)[name].as<T>();
        }
    }

    // Now check whether the value is provided in the XML file
    tinyxml2::XMLHandle childHandle = get_tinyxml_handle( handle, name );
    tinyxml2::XMLText* text = childHandle.FirstChild().ToText();
    if( text ) {
        try {
            return boost::lexical_cast<T>( text->Value() );
        } catch( const boost::bad_lexical_cast& ) {
            throw std::runtime_error( "XML variable \"" + get_xml_handle_path( handle, name ) + "\", with value \"" +
                                      text->Value() + "\" could not be parsed as a " + typeid( T ).name() );
        }
    }
    // Otherwise throw an exception
    throw std::runtime_error( "XML variable \"" + get_xml_handle_path( handle, name ) +
                              "\", intended to be parsed as a " + typeid( T ).name() +
                              ", could not be found in the .xml file" );
}

//
// Function: get_int
//
// Gets an int from XML.
//
// Parameters:
//	name - the name of the field in XML that contains the int.
//  defaultValue - the value this function returns if "name" is not specified in XML.
//  handle - a tinyxml2::XMLHandle to the node that contains the "name" field.
//  warnOnFail - a flag indicating that a warning should be made if the value of "name" could not be cast to an int.
//
// Returns:
//  The int specified by "name" or "defaultValue" if "name" was not found.
//
inline int get_int( const std::string& name, int defaultValue, tinyxml2::XMLHandle handle, bool warnOnFail = true,
                    boost::program_options::variables_map* cmd_line_vars = 0 ) {
    return get_xml_value<int>( name, defaultValue, handle, warnOnFail, cmd_line_vars );
}

inline int get_int( const std::string& name, tinyxml2::XMLHandle handle,
                    boost::program_options::variables_map* cmd_line_vars = 0 ) {
    return get_xml_value<int>( name, handle, cmd_line_vars );
}

inline float get_float( const std::string& name, float defaultValue, tinyxml2::XMLHandle handle, bool warnOnFail = true,
                        boost::program_options::variables_map* cmd_line_vars = 0 ) {
    return get_xml_value<float>( name, defaultValue, handle, warnOnFail, cmd_line_vars );
}

inline float get_float( const std::string& name, tinyxml2::XMLHandle handle,
                        boost::program_options::variables_map* cmd_line_vars = 0 ) {
    return get_xml_value<float>( name, handle, cmd_line_vars );
}

inline bool get_bool( const std::string& name, bool defaultValue, tinyxml2::XMLHandle handle, bool warnOnFail = true,
                      boost::program_options::variables_map* cmd_line_vars = 0 ) {
    std::string v;
    try {
        v = get_xml_value<std::string>( name, handle, cmd_line_vars );
        v = frantic::strings::to_lower( v );
        v = frantic::strings::trim( v );

        if( v == "true" || v == "on" || v == "1" || v == "yes" || v == "oui" )
            return true;

        if( v == "false" || v == "off" || v == "0" || v == "no" || v == "non" )
            return false;

    } catch( std::exception& ) {
        if( ( warnOnFail && frantic::logging::is_logging_warnings() ) || frantic::logging::is_logging_debug() )
            std::cerr << "Warning: Couldn't read xml parameter \"" + get_xml_handle_path( handle, name ) +
                             "\", assuming default of "
                      << defaultValue << std::endl;
        return defaultValue;
    }

    throw std::runtime_error( "XML variable \"" + get_xml_handle_path( handle, name ) + "\" with value \"" + v +
                              "\", intended to be parsed as a bool, could not be parsed properly" );
}

inline bool get_bool( const std::string& name, tinyxml2::XMLHandle handle ) {
    std::string v;
    try {
        v = get_xml_value<std::string>( name, handle );
        v = frantic::strings::to_lower( v );
        v = frantic::strings::trim( v );

        if( v == "true" || v == "on" || v == "1" || v == "yes" || v == "oui" )
            return true;

        if( v == "false" || v == "off" || v == "0" || v == "no" || v == "non" )
            return false;

        throw std::runtime_error( "XML variable \"" + get_xml_handle_path( handle, name ) + " with value \"" + v +
                                  "\", intended to be parsed as a bool, could not be parsed properly" );
    } catch( std::exception& ) {
        throw std::runtime_error( "XML variable \"" + get_xml_handle_path( handle, name ) +
                                  ", intended to be parsed as a bool, could not be parsed properly" );
    }
}

inline std::string get_string( const std::string& name, const std::string& defaultValue, tinyxml2::XMLHandle handle,
                               bool warnOnFail = true, boost::program_options::variables_map* cmd_line_vars = 0 ) {
    return get_xml_value<std::string>( name, defaultValue, handle, warnOnFail, cmd_line_vars );
}

inline std::string get_string( const std::string& name, tinyxml2::XMLHandle handle,
                               boost::program_options::variables_map* cmd_line_vars = 0 ) {
    return get_xml_value<std::string>( name, handle, cmd_line_vars );
}

// Function: set_node_value
// Sets the value of the element that is pointed to by the specified handle.
//  If the element already has a tinyxml2::XMLText child, it simply sets the value of the text.
//  If it doesn't exists, it creates a tinyxml2::XMLText with the specified value and inserts it
//  as a child of element in handle.  The tinyxml2::XMLText will always be the first child of the
//  element in handle.
// Parameters:
//	handle - the element who's value we are setting
//  value - the new value of the element.
// Returns:
//  true if the set value was successful, false otherwise.
inline bool set_node_value( tinyxml2::XMLHandle handle, std::string value ) {
    bool success = false;
    tinyxml2::XMLText* text = handle.FirstChild().ToText();
    if( text != nullptr ) {
        text->SetValue( value.c_str() );
        success = true;
    } else {
		tinyxml2::XMLNode* currentNode = handle.ToNode();
		if( !currentNode ) {
			return false;
		}
		if( handle.ToNode()->InsertFirstChild( currentNode->GetDocument()->NewText( value.c_str() ) ) ) {
			success = true;
		}
    }
    return success;
}

// Function: attribute_exists
//  Determines if the attribute exists for the given handle.
// Parameters:
//		handle - the element to check if it contains the attribute
//		attributeName - the name of the attribute we are checking for existence.
// Returns:
//		true if an attribute with the name "attributeName" was found in the handle, false otherwise.
inline bool attribute_exists( tinyxml2::XMLHandle handle, const std::string& attributeName ) {
    tinyxml2::XMLElement* element = handle.ToElement();
    if( element && element->Attribute( attributeName.c_str() ) )
        return true;
    return false;
}

// Function: get_all_attribute_names
//  Gets all the attribute names of the given handle.
// Parameters:
//	handle - the handle containing all the attributes.
// Returns:
//	A vector of string which are all the attribute names of the given handle.
inline std::vector<std::string> get_all_attribute_names( tinyxml2::XMLHandle handle ) {
    std::vector<std::string> attributeNames;

    tinyxml2::XMLElement* element = handle.ToElement();
    if( element ) {
        const tinyxml2::XMLAttribute* attribute = element->FirstAttribute();
        while( attribute ) {
            std::string attributeName = attribute->Name();
            attributeNames.push_back( attributeName );
            attribute = attribute->Next();
        }
    }

    return attributeNames;
}

// Function: get_attribute_value
//		Gets the value of the attribute with the given attribute name.
// Parameters:
//	handle - the handle to the node that contains the attribute.
//  attributeName - the name of the attribute of who's value we're interested in.
//  attributeValue - a field to store the value of the attribute.
// Returns: the true if it was successful in retrieving the value of the attribute, false otherwise.
inline bool get_attribute_value( tinyxml2::XMLHandle handle, const std::string& attributeName, std::string& attributeValue ) {
    tinyxml2::XMLElement* element = handle.ToElement();
    if( element ) {
        if( element->Attribute( attributeName.c_str() ) ) {
            attributeValue = *element->Attribute( attributeName.c_str() );
            return true;
        }
    }
    return false;
}

// Function: set_attribute_value
//  Sets the value of the attribute specified by "attributeName" to "attributeValue".
//  This method does nothing if the attribute does not exists.
// Parameters:
//  handle - a handle to the node which contains the attribute.
//	attributeName - the name of the attribute to set
//  attribuetValue- the new value of the attribute (as a string)
// Returns: true if the attribute's value was set, false if the attribute doesn't exist.
inline bool set_attribute_value( tinyxml2::XMLHandle handle, const std::string& attributeName,
                                 const std::string& attributeValue ) {
    tinyxml2::XMLElement* element = handle.ToElement();
    if( element ) {
        if( element->Attribute( attributeName.c_str() ) ) {
            element->SetAttribute( attributeName.c_str(), attributeValue.c_str() );
            return true;
        }
    }
    return false;
}

// Function: add_attribute
//  Adds the attribute (with the specified name and value) to the current handle.
//  This doesn't do anything if an attribute with a name of "attributeName" already exists or
//   if the current handle doesn't point to an element.
// Parameters:
//		handle - a handle to the node which the attribute should be added to.
//		attributeName - The name of the attribute to add.
//		attributeValue- The value of the attribute to add.
// Returns: true if the attribute was successfully added, false otherwise.
inline bool add_attribute( tinyxml2::XMLHandle handle, const std::string& attributeName, const std::string& attributeValue ) {
    tinyxml2::XMLElement* element = handle.ToElement();
    if( element ) {
        if( element->Attribute( attributeName.c_str() ) ) 
        {
			// We only add the attribute if one with the same name doesn't already exist.
            element->SetAttribute( attributeName.c_str(), attributeValue.c_str() );
            return true;
        }
    }
    return false;
}

// Function: remove_attribute
//  Removes the attribute with the name "attributeName".
// Parameters:
//		handle - a handle to the node which contains the attribute to be removed
//		attributeName - the name of the attribute to remove
// Returns:
//      true if the attribute existed and was removed, false otherwise.
inline bool remove_attribute( tinyxml2::XMLHandle handle, const std::string& attributeName ) {
    tinyxml2::XMLElement* element = handle.ToElement();
    if( element ) {
        if( element->Attribute( attributeName.c_str() ) ) // We only remove the attribute if one exists.
        {
            element->DeleteAttribute( attributeName.c_str() );
            return true;
        }
    }
    return false;
}

template <typename T>
inline T safeCmdLineGet( boost::program_options::variables_map& cmd_line_vars, std::string name ) {
    try {
        // my copy of gcc chokes on the commented-out line below.
        // The following two lines should be equivalent -PTF
        const boost::any& val = cmd_line_vars[name].value();
        return boost::any_cast<T>( val );
        // return cmd_line_vars[name].as<T>();
    } catch( boost::bad_any_cast const& e ) {
        throw std::runtime_error( "Unable to convert command line variable \"" + name + "\" to a " +
                                  typeid( T() ).name() + " " + e.what() );
    }
}

/**
 *  Return the XML document's encoding declaration, if it exists.  Otherwise
 * return an empty string.
 *
 * @param doc the XML document to process.
 * @return the XML document's encoding declaration, if it exists.  Returns an empty
 *        string otherwise.
 */
inline std::string get_encoding_declaration( const tinyxml2::XMLDocument& doc ) {
    return std::string( "UTF-8" );
}

/**
 * Parses the xml file into a string which only contains the elements supported by TinyXml
 * Removes elements like <?xml... or <!DOCTYPE... or [<!ENTITY...
 *
 * @param filename The file from which the elements must be removed
 *
 * @return Returns the xmlString of the file which only contains supported elements
 */
inline std::string remove_unsupported_xml_elements( const boost::filesystem::path& filename ) {
    std::ifstream stream( filename.c_str(), std::ios::binary | std::ios::ate );
    size_t fileSize = stream.tellg();
    stream.seekg( 0, std::ios::beg );
    std::vector<char> bytes( fileSize );
    stream.read( &bytes[0], fileSize );

    std::string xmlString( &bytes[0], fileSize );

    size_t pos = 0;
    while( true ) {
        pos = xmlString.find_first_of( "<", pos );
        if( pos == -1 ) {
            // We couldn't find a "<". If the first character in the file was a "?" or "!"
            // we would enter an infinite loop because xmlString[pos + 1] would always be a "?" or "!"
            // and the next condition would always evaluate to true.
            throw std::runtime_error( "remove_unsupported_xml_elements(): Input file " + filename.string() +
                                      " contains no XML elements." );
        }
        if( xmlString[pos + 1] == '?' ||  // <?xml...
            xmlString[pos + 1] == '!' ) { // <!DOCTYPE... or [<!ENTITY...
            // Skip this line
            pos = xmlString.find_first_of( "\n", pos );
        } else
            break;
    }
    xmlString = xmlString.substr( pos );
    return xmlString;
}

} // namespace tinyxml
} // namespace frantic
