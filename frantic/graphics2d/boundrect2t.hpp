// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <boost/lexical_cast.hpp>

#include <frantic/graphics2d/vector2f.hpp>
#include <frantic/graphics2d/vector2t.hpp>

namespace frantic {
namespace graphics2d {

template <class VectorType, class SizeType, class BoundRectType>
class boundrect2t {
  protected:
    typedef typename VectorType::value_type value_type;

    VectorType m_minimum, m_maximum;

  public:
    //////////////
    // Constructors
    //////////////
    boundrect2t( const VectorType& position )
        : m_minimum( position )
        , m_maximum( position ) {}

    boundrect2t( const VectorType& minimum, const VectorType& maximum )
        : m_minimum( minimum )
        , m_maximum( maximum ) {}

    boundrect2t( value_type xmin, value_type xmax, value_type ymin, value_type ymax )
        : m_minimum( xmin, ymin )
        , m_maximum( xmax, ymax ) {}

    //////////////
    // Setters
    //////////////

    void set_minimum( const VectorType& minimum ) { m_minimum = minimum; }

    void set_maximum( const VectorType& maximum ) { m_maximum = maximum; }

    void set( const VectorType& minimum, const VectorType& maximum ) {
        m_minimum = minimum;
        m_maximum = maximum;
    }

    void set( const VectorType& minimum, const SizeType& size ) {
        m_minimum = minimum;
        m_maximum = minimum + size;
    }

    //////////////
    // Member properties
    //////////////

    VectorType minimum() const { return m_minimum; }

    VectorType maximum() const { return m_maximum; }

    VectorType& minimum() { return m_minimum; }

    VectorType& maximum() { return m_maximum; }

    bool is_empty() const {
        return static_cast<const BoundRectType*>( this )->xsize() <= 0 ||
               static_cast<const BoundRectType*>( this )->ysize() <= 0;
    }

    vector2f center() const {
        return vector2f( m_minimum.get_x() + static_cast<const BoundRectType*>( this )->xsize() / 2.0f,
                         m_minimum.get_y() + static_cast<const BoundRectType*>( this )->ysize() / 2.0f );
    }

    SizeType size() const {
        return SizeType( static_cast<const BoundRectType*>( this )->xsize(),
                         static_cast<const BoundRectType*>( this )->ysize() );
    }

    value_type get_area() const {
        if( static_cast<const BoundRectType*>( this )->is_empty() )
            return 0;
        else
            return static_cast<const BoundRectType*>( this )->xsize() *
                   static_cast<const BoundRectType*>( this )->ysize();
    }

    VectorType get_corner( int i ) const {
        switch( i ) {
        case 0:
            return VectorType( m_minimum.x, m_minimum.y );
        case 1:
            return VectorType( m_minimum.x, m_maximum.y );
        case 2:
            return VectorType( m_maximum.x, m_minimum.y );
        case 3:
            return VectorType( m_maximum.x, m_maximum.y );
        default:
            throw std::runtime_error(
                "boundrect2t.get_corner: Tried to get an invalid corner (must be between 0 and 3)" );
        }
    }

    bool contains( const VectorType& position ) const {
        return m_minimum.x <= position.x && position.x <= m_maximum.x && m_minimum.y <= position.y &&
               position.y <= m_maximum.y;
    }

    bool contains_as_open_set( const VectorType& position ) const {
        return m_minimum.x < position.x && position.x < m_maximum.x && m_minimum.y < position.y &&
               position.y < m_maximum.y;
    }

    bool contains( const boundrect2t<VectorType, SizeType, BoundRectType>& bounds ) const {
        return contains( bounds.minimum() ) && contains( bounds.maximum() );
    }

    VectorType clamp( VectorType v ) const {
        if( is_empty() ) {
            throw std::runtime_error( "boundrect2t: Cannot clamp vector " + boost::lexical_cast<std::string>( v ) +
                                      " because rectangle is empty" );
        }
        if( v.x < m_minimum.x )
            v.x = m_minimum.x;
        if( v.y < m_minimum.y )
            v.y = m_minimum.y;
        if( v.x > m_maximum.x )
            v.x = m_maximum.x;
        if( v.y > m_maximum.y )
            v.y = m_maximum.y;

        return v;
    }

    //////////////
    // Operators
    //////////////

    bool expand( value_type distance ) { return expand( VectorType( distance ) ); }

    bool expand( const VectorType& distance ) {
        if( !is_empty() ) {
            m_minimum -= distance;
            m_maximum += distance;

            return true;
        } else
            return false;
    }

    BoundRectType& operator+=( const boundrect2t<VectorType, SizeType, BoundRectType>& box ) {
        if( is_empty() ) {
            m_minimum = box.minimum();
            m_maximum = box.maximum();
        } else if( !box.is_empty() ) {
            if( m_minimum.x > box.m_minimum.x )
                m_minimum.x = box.m_minimum.x;
            if( m_maximum.x < box.m_maximum.x )
                m_maximum.x = box.m_maximum.x;
            if( m_minimum.y > box.m_minimum.y )
                m_minimum.y = box.m_minimum.y;
            if( m_maximum.y < box.m_maximum.y )
                m_maximum.y = box.m_maximum.y;
        }
        return *static_cast<BoundRectType*>( this );
    }

    BoundRectType& operator+=( const VectorType& position ) {
        if( is_empty() ) {
            m_minimum = position;
            m_maximum = position;
        } else {
            if( m_minimum.x > position.x )
                m_minimum.x = position.x;
            if( m_maximum.x < position.x )
                m_maximum.x = position.x;
            if( m_minimum.y > position.y )
                m_minimum.y = position.y;
            if( m_maximum.y < position.y )
                m_maximum.y = position.y;
        }
        return *static_cast<BoundRectType*>( this );
    }

    std::string str() const;
};

template <class VectorType, class SizeType, class BoundRectType>
bool operator==( const boundrect2t<VectorType, SizeType, BoundRectType>& a,
                 const boundrect2t<VectorType, SizeType, BoundRectType>& b ) {
    return a.minimum() == b.minimum() && a.maximum() == b.maximum();
}

template <class VectorType, class SizeType, class BoundRectType>
bool operator!=( const boundrect2t<VectorType, SizeType, BoundRectType>& a,
                 const boundrect2t<VectorType, SizeType, BoundRectType>& b ) {
    return a.minimum() != b.minimum() || a.maximum() != b.maximum();
}

template <class VectorType, class SizeType, class BoundRectType>
inline std::ostream& operator<<( std::ostream& out, const boundrect2t<VectorType, SizeType, BoundRectType>& box ) {
    out << "(boundrect2t " << box.minimum() << " " << box.maximum() << " )";
    return out;
}

template <class VectorType, class SizeType, class BoundRectType>
inline std::string boundrect2t<VectorType, SizeType, BoundRectType>::str() const {
    std::stringstream ss;
    ss << *static_cast<const BoundRectType*>( this );
    return ss.str();
}

} // namespace graphics2d
} // namespace frantic
