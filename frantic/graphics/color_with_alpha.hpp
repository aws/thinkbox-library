// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/alpha1f.hpp>
#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/color3f.hpp>

namespace frantic {
namespace graphics {

// In general this represents a color with alpha in premultiplied form.  That means it is in the form
// (alpha*color,alpha).
template <class ColorType, class AlphaType>
class color_with_alpha {
  public:
    ColorType c;
    AlphaType a;

    enum { Channels = ColorType::Channels + AlphaType::Channels };

    color_with_alpha() {}

    explicit color_with_alpha( const ColorType& color )
        : c( color ) {}

    // Note: this assumes you have already premultiplied the color with the alpha
    color_with_alpha( const ColorType& color, const AlphaType& alpha )
        : c( color )
        , a( alpha ) {}

    explicit color_with_alpha( float intensity )
        : c( intensity ) {}

    color_with_alpha( float red, float green, float blue )
        : c( red, green, blue ) {}

    color_with_alpha( float red, float green, float blue, float alpha )
        : c( red, green, blue )
        , a( alpha ) {}

    static color_with_alpha<ColorType, AlphaType> from_nonpremultiplied( const ColorType& color,
                                                                         const AlphaType& alpha ) {
        return color_with_alpha<ColorType, AlphaType>( alpha.premultiply( color ), alpha );
    }

    // Accessors to return pure color and alpha
    ColorType& color() { return c; }

    const ColorType& color() const { return c; }

    AlphaType& alpha() { return a; }

    const AlphaType& alpha() const { return a; }

    //////////////
    // Queries
    //////////////

    float component_sum() const { return c.component_sum() + a.component_sum(); }

    static color_with_alpha<ColorType, AlphaType> abs( const color_with_alpha<ColorType, AlphaType>& color ) {
        return color_with_alpha<ColorType, AlphaType>( ColorType::abs( color.c ), AlphaType::abs( color.a ) );
    }

    //////////////
    // Operators
    //////////////

    // Assumes premultiplied alpha.  Suppose we want (alpha1*color1, alpha1) over (alpha2*color2, alpha2).  Then,
    // alpha3 = 1 - (1-alpha1)*(1-alpha2) = alpha1.blend_over(alpha2)
    // alpha3 * color3 = alpha1*color1 + (1-alpha1)*alpha2*color2 = alpha1*color1 + alpha1.occlude(alpha2*color2)
    void blend_over( const color_with_alpha<ColorType, AlphaType>& color ) {
        c += a.occlude( color.c );
        a.blend_over( color.a );
    }

    void blend_under( const color_with_alpha<ColorType, AlphaType>& color ) {
        c = color.c + color.a.occlude( c );
        a.blend_under( color.a );
    }

    color_with_alpha<color3f, alpha1f> to_color4f() const {
        return color_with_alpha<color3f, alpha1f>( c, a.to_float() );
    }

    color_with_alpha& operator+=( const color_with_alpha& cwa ) {
        c += cwa.c;
        a += cwa.a;
        return *this;
    }

    color_with_alpha& operator*=( float x ) {
        c *= x;
        a *= x;
        return *this;
    }

    color_with_alpha operator/( float x ) { return color_with_alpha( c / x, a / x ); }

    color_with_alpha operator/=( float x ) {
        c /= x;
        a /= x;
        return *this;
    }

    bool operator!=( const color_with_alpha& other ) { return c != other.c || a != other.a; }

    bool operator==( const color_with_alpha& other ) { return c == other.c && a == other.a; }

    // A parse function, which is specialized for the various types we use
    static color_with_alpha parse( const frantic::tstring& input );

    static color_with_alpha parse_maxscript( const frantic::tstring& input );
};

template <class ColorType, class AlphaType>
color_with_alpha<ColorType, AlphaType> operator*( float x, const color_with_alpha<ColorType, AlphaType>& color ) {
    return color_with_alpha<ColorType, AlphaType>( x * color.c, x * color.a );
}

template <class ColorType, class AlphaType>
color_with_alpha<ColorType, AlphaType> operator*( const color_with_alpha<ColorType, AlphaType>& color, float x ) {
    return color_with_alpha<ColorType, AlphaType>( color.c * x, color.a * x );
}

typedef color_with_alpha<color3f, alpha1f> color4f;
typedef color_with_alpha<color3f, alpha3f> color6f;

// g++ gives errors saying the following two methods have invalid declaration
#if defined( _WIN32 ) || defined( _WIN64 )
template <>
inline color_with_alpha<color3f, alpha3f> color_with_alpha<color3f, alpha3f>::parse( const frantic::tstring& input ) {
    std::vector<frantic::tstring> values;
    frantic::strings::split( frantic::strings::string_replace( input, _T("color "), _T("") ), values, _T("()[], ") );
    if( values.size() == 1 ) {
        float v;
        try {
            v = boost::lexical_cast<float>( values[0] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                      frantic::strings::to_string( input ) + "\"" );
        }
        return color6f( color3f( v ), alpha3f( 1.f ) );
    } else if( values.size() == 3 ) {
        float r, g, b;
        try {
            r = boost::lexical_cast<float>( values[0] );
            g = boost::lexical_cast<float>( values[1] );
            b = boost::lexical_cast<float>( values[2] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                      frantic::strings::to_string( input ) + "\"" );
        }

        return color6f( color3f( r, g, b ), alpha3f( 1.f ) );
    } else if( values.size() == 4 ) {
        float r, g, b, a;
        try {
            r = boost::lexical_cast<float>( values[0] );
            g = boost::lexical_cast<float>( values[1] );
            b = boost::lexical_cast<float>( values[2] );
            a = boost::lexical_cast<float>( values[3] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                      frantic::strings::to_string( input ) + "\"" );
        }

        return color6f( color3f( r, g, b ), alpha3f( a ) );
    } else if( values.size() == 6 ) {
        float r, g, b, ar, ag, ab;
        try {
            r = boost::lexical_cast<float>( values[0] );
            g = boost::lexical_cast<float>( values[1] );
            b = boost::lexical_cast<float>( values[2] );
            ar = boost::lexical_cast<float>( values[3] );
            ag = boost::lexical_cast<float>( values[4] );
            ab = boost::lexical_cast<float>( values[5] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                      frantic::strings::to_string( input ) + "\"" );
        }

        return color6f( color3f( r, g, b ), alpha3f( ar, ag, ab ) );
    } else {
        throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                  frantic::strings::to_string( input ) + "\"" );
    }
}

template <>
inline color_with_alpha<color3f, alpha3f>
color_with_alpha<color3f, alpha3f>::parse_maxscript( const frantic::tstring& input ) {
    std::vector<frantic::tstring> values;
    frantic::strings::split( frantic::strings::string_replace( input, _T("color "), _T("") ), values, _T("()[], ") );
    if( values.size() == 3 ) {
        float r, g, b;
        try {
            r = boost::lexical_cast<float>( values[0] );
            g = boost::lexical_cast<float>( values[1] );
            b = boost::lexical_cast<float>( values[2] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                      frantic::strings::to_string( input ) + "\"" );
        }

        return color6f( color3f( r / 255.f, g / 255.f, b / 255.f ), alpha3f( 1.f ) );
    } else if( values.size() == 4 ) {
        float r, g, b, a;
        try {
            r = boost::lexical_cast<float>( values[0] );
            g = boost::lexical_cast<float>( values[1] );
            b = boost::lexical_cast<float>( values[2] );
            a = boost::lexical_cast<float>( values[3] );
        } catch( boost::bad_lexical_cast& ) {
            throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                      frantic::strings::to_string( input ) + "\"" );
        }

        return color6f( color3f( r / 255.f, g / 255.f, b / 255.f ), alpha3f( a / 255.f ) );
    } else {
        throw std::runtime_error( "color3f::parse: Couldn't parse input value \"" +
                                  frantic::strings::to_string( input ) + "\"" );
    }
}
#endif

template <class Color, class Alpha>
inline color_with_alpha<Color, Alpha> operator+( const color_with_alpha<Color, Alpha>& c1,
                                                 const color_with_alpha<Color, Alpha>& c2 ) {
    return color_with_alpha<Color, Alpha>( c1.c + c2.c, c1.a + c2.a );
}

#ifndef FRANTIC_DISABLE_IOSTREAM
inline std::ostream& operator<<( std::ostream& out, const color4f& a ) {
    out << "(color " << a.c.r << ", " << a.c.g << ", " << a.c.b << ", " << a.a.a << " )";
    return out;
}

inline std::ostream& operator<<( std::ostream& out, const color6f& a ) {
    out << "(color " << a.c.r << ", " << a.c.g << ", " << a.c.b << ", " << a.a.ar << ", " << a.a.ag << ", " << a.a.ab
        << " )";
    return out;
}
#endif

} // namespace graphics
} // namespace frantic
