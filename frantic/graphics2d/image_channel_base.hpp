// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace frantic {
namespace graphics2d {
// This should never be done, importing a namespace into another is a bad idea.
// using namespace std;

class image_channel_base {

  protected:
    image_channel_base* m_parent;
    size2 m_size;
    std::string m_name;

    image_channel_base( image_channel_base* parent, const std::string& name = "UNNAMED" )
        : m_parent( parent )
        , m_size()
        , m_name( name ) {
        if( !m_parent )
            throw std::runtime_error( "frantic::graphics2d::image_channel_base() : invalid parent" );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
  public:
    image_channel_base( const std::string& name = "UNNAMED" )
        : m_parent( 0 )
        , m_size( 0, 0 )
        , m_name( name ) {}
    image_channel_base( size2 size )
        : m_parent( 0 )
        , m_size( size )
        , m_name( "UNNAMED" ) {}
    image_channel_base( const std::string& name, size2 size )
        : m_parent( 0 )
        , m_size( size )
        , m_name( name ) {}

    virtual ~image_channel_base() {}

    const size2& size() const { return m_size; }

    frantic::graphics2d::size2f size2f() const { return frantic::graphics2d::size2f( m_size ); }

    int xsize() const { return m_size.xsize; }

    int ysize() const { return m_size.ysize; }

    const std::string& name() const { return m_name; }
};
} // namespace graphics2d
} // namespace frantic
