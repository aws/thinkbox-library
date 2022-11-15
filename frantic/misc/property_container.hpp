// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

//#pragma message( "frantic/misc/property_container.hpp is deprecated!" )

#include <frantic/strings/tstring.hpp>
#include <map>
#include <set>
#include <string>

namespace frantic {
namespace properties {

// This class holds properties in the form of a string->string std::map
class property_container {
    std::map<frantic::tstring, frantic::tstring> m_properties;
    std::set<frantic::tstring> m_validNames;
    bool m_nameValidation;

  public:
    property_container( bool nameValidation = false )
        : m_nameValidation( nameValidation ) {}

    std::map<frantic::tstring, frantic::tstring>& storage() { return m_properties; }

    void erase( const frantic::tstring& name ) { m_properties.erase( name ); }

    bool exists( const frantic::tstring& name ) const { return m_properties.find( name ) != m_properties.end(); }

    bool valid_property_name( const frantic::tstring& name ) const {
        if( m_nameValidation )
            return m_validNames.find( name ) != m_validNames.end();
        else
            return true;
    }

    // When name validation is on, use this function to add valid names of properties
    void add( const frantic::tstring& name ) { m_validNames.insert( name ); }

    // If a property exists in the array, return it even if it's not valid.  Otherwise, throw an exception.
    frantic::tstring get( const frantic::tstring& name ) const {
        std::map<frantic::tstring, frantic::tstring>::const_iterator i = m_properties.find( name );
        if( i != m_properties.end() ) {
            return i->second;
        } else {
            if( m_nameValidation && m_validNames.find( name ) == m_validNames.end() )
                throw std::runtime_error( "Attempted to get property \"" + frantic::strings::to_string( name ) +
                                          "\", which is not a valid property name." );
            else
                throw std::runtime_error( "Property \"" + frantic::strings::to_string( name ) +
                                          "\" does not have a value stored in the properties array." );
        }
    }
    void set( const frantic::tstring& name, const frantic::tstring& value ) {
        if( m_nameValidation ) {
            if( m_validNames.find( name ) == m_validNames.end() )
                throw std::runtime_error( "Attempted to set property \"" + frantic::strings::to_string( name ) +
                                          "\", which is not a valid property name." );
        }
        m_properties[name] = value;
    }
    void set( const frantic::tstring& name, const frantic::tchar* value ) {
        if( m_nameValidation ) {
            if( m_validNames.find( name ) == m_validNames.end() )
                throw std::runtime_error( "Attempted to set property \"" + frantic::strings::to_string( name ) +
                                          "\", which is not a valid property name." );
        }
        m_properties[name] = value;
    }

    bool get_bool( const frantic::tstring& name ) const {
        frantic::tstring value = frantic::strings::trim( frantic::strings::to_lower( get( name ) ) );
        if( value == _T("true") || value == _T("1") || value == _T("yes") || value == _T("on") )
            return true;
        if( value == _T("false") || value == _T("off") || value == _T("0") || value == _T("no") )
            return false;
        throw std::runtime_error( "Bool Property \"" + frantic::strings::to_string( name ) + "\" had invalid value \"" +
                                  frantic::strings::to_string( value ) + "\"" );
    }
    void set( const frantic::tstring& name, bool value ) {
        if( value )
            set( name, _T("true") );
        else
            set( name, _T("false") );
    }

    int get_int( const frantic::tstring& name ) const {
        frantic::tstring value = get( name );
        try {
            return boost::lexical_cast<int>( value );
        } catch( const boost::bad_lexical_cast& ) {
            throw std::runtime_error( "Integer Property \"" + frantic::strings::to_string( name ) +
                                      "\" had invalid value \"" + frantic::strings::to_string( value ) + "\"" );
        }
    }
    void set( const frantic::tstring& name, int value ) { set( name, boost::lexical_cast<frantic::tstring>( value ) ); }

    float get_float( const frantic::tstring& name ) const {
        frantic::tstring value = get( name );
        try {
            return boost::lexical_cast<float>( value );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "Float Property \"" + frantic::strings::to_string( name ) +
                                      "\" had invalid value \"" + frantic::strings::to_string( value ) + "\"" );
        }
    }
    void set( const frantic::tstring& name, float value ) {
        set( name, boost::lexical_cast<frantic::tstring>( value ) );
    }
};

} // namespace properties
} // namespace frantic
