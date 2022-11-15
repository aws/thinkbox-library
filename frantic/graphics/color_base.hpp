// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "frantic/channels/channel_map.hpp"
#include "frantic/math/math_array.hpp"
#include <boost/static_assert.hpp>

namespace frantic {
namespace graphics {

using frantic::channels::channel_map;
using frantic::math::math_array;

template <class T, int Arity>
struct color_array : public math_array<T, Arity> {};

template <class T, int Arity>
struct alpha_array : public math_array<T, Arity> {};

/**
 * The core class for the storage of color and alpha information.
 *
 * @note Any class that inherits this class should add a conversion
 *	function to @ref color_converter.hpp
 *
 * @author   Brian McKinnon
 * @since    Apr 25, 2007
 */
template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
class color_base : private color_array<DataType, ColorCount>, private alpha_array<DataType, AlphaCount> {
    BOOST_STATIC_ASSERT( ColorCount > 0 );
    BOOST_STATIC_ASSERT( ( AlphaCount == 0 ) || ( AlphaCount == 1 ) || ( AlphaCount == ColorCount ) );

    typedef color_array<DataType, ColorCount> color_array_t;
    typedef alpha_array<DataType, AlphaCount> alpha_array_t;

  public:
    typedef DataType Type;

    static int get_color_count() { return ColorCount; }
    static int get_alpha_count() { return AlphaCount; }

    DataType get_color( int i ) const { return color_array_t::get_data( i ); }
    DataType get_alpha( int i ) const { return alpha_array_t::get_data( i ); }

    void set_color( int i, DataType d ) { return get_color().set_data( i, d ); }
    void set_alpha( int i, DataType d ) { return get_alpha().set_data( i, d ); }

    math_array<DataType, ColorCount>& get_color() { return *static_cast<color_array_t*>( this ); }
    math_array<DataType, ColorCount> get_color() const { return *static_cast<const color_array_t*>( this ); }
    math_array<DataType, AlphaCount>& get_alpha() { return *static_cast<alpha_array_t*>( this ); }
    math_array<DataType, AlphaCount> get_alpha() const { return *static_cast<const alpha_array_t*>( this ); }

    PixelType operator-() const;
    PixelType operator+( const PixelType& p ) const;
    PixelType operator-( const PixelType& p ) const;
    PixelType operator*( const PixelType& p ) const;
    PixelType operator/( const PixelType& p ) const;

    PixelType& operator+=( const PixelType& p );
    PixelType& operator-=( const PixelType& p );
    PixelType& operator*=( const PixelType& p );
    PixelType& operator/=( const PixelType& p );

    PixelType operator+( DataType d ) const;
    PixelType operator-( DataType d ) const;
    PixelType operator*( DataType d ) const;
    PixelType operator/( DataType d ) const;

    PixelType& operator+=( DataType d );
    PixelType& operator-=( DataType d );
    PixelType& operator*=( DataType d );
    PixelType& operator/=( DataType d );

    PixelType get_abs() const;
    PixelType get_sqrt() const;
    DataType get_sum() const {
        return color_array<DataType, ColorCount>::get_sum() + alpha_array<DataType, AlphaCount>::get_sum();
    }
    DataType get_mean() const { return get_sum() / ( ColorCount + AlphaCount ); }

    bool operator==( const PixelType& p ) const;
    bool operator!=( const PixelType& p ) const;

    /** Applies the alpha channel values to the colors channels.
     * This generates a premultiplied alpha pixel.
     *
     * \return A premultiplied alpha pixel
     */
    PixelType apply_alpha_to_color() const;

    /** Removes the alpha channel values from the colors channels.
     * This removes the premultiplied alpha.
     * \note If the alpha is zero this function will generate a
     *		  divide by zero error
     *
     * \return A color pixel without alpha applied
     */
    PixelType remove_alpha_from_color() const;

    /** Uses the alpha values to blend the color and alpha
     * channels together.
     *
     * \arg	p	The background pixel that this pixel will merge onto.
     * \return	The blended pixel
     */
    PixelType blend_over( const PixelType& p ) const;

    /** Uses the alpha values to blend the color and alpha
     * channels together.
     *
     * \arg	p	The background pixel that this pixel will merge onto.
     */
    void apply_blend_over( const PixelType& p );

    /** Uses the alpha values to blend the color and alpha
     * channels together.
     *
     * \arg	p	The foreground pixel that will merge onto this pixel.
     * \return	The blended pixel
     */
    PixelType blend_under( const PixelType& p ) const;

    /** Uses the alpha values to blend the color and alpha
     * channels together.
     *
     * \arg	p	The foreground pixel that will merge onto this pixel.
     */
    void apply_blend_under( const PixelType& p );
};

#if defined( _MSC_VER )
#pragma warning( push, 3 )
#pragma warning( disable : 4127 )
#endif

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator-() const {
    PixelType result;

    result.get_color() = -this->get_color();

    if( AlphaCount > 0 ) {
        result.get_alpha() = -this->get_alpha();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator+( const PixelType& p ) const {
    PixelType result;

    result.get_color() = this->get_color() + p.get_color();

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() + p.get_alpha();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator-( const PixelType& p ) const {
    PixelType result;

    result.get_color() = this->get_color() - p.get_color();

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() - p.get_alpha();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator*( const PixelType& p ) const {
    PixelType result;

    result.get_color() = this->get_color() * p.get_color();

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() * p.get_alpha();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator/( const PixelType& p ) const {
    PixelType result;

    result.get_color() = this->get_color() / p.get_color();

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() / p.get_alpha();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator+=( const PixelType& p ) {
    this->get_color() += p.get_color();

    if( AlphaCount > 0 ) {
        this->get_alpha() += p.get_alpha();
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator-=( const PixelType& p ) {
    this->get_color() -= p.get_color();

    if( AlphaCount > 0 ) {
        this->get_alpha() -= p.get_alpha();
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator*=( const PixelType& p ) {
    this->get_color() *= p.get_color();

    if( AlphaCount > 0 ) {
        this->get_alpha() *= p.get_alpha();
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator/=( const PixelType& p ) {
    this->get_color() /= p.get_color();

    if( AlphaCount > 0 ) {
        this->get_alpha() /= p.get_alpha();
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator+( DataType d ) const {
    PixelType result;

    result.get_color() = this->get_color() + d;

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() + d;
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator-( DataType d ) const {
    PixelType result;

    result.get_color() = this->get_color() - d;

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() - d;
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator*( DataType d ) const {
    PixelType result;

    result.get_color() = this->get_color() * d;

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() * d;
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::operator/( DataType d ) const {
    PixelType result;

    result.get_color() = this->get_color() / d;

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha() / d;
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator+=( DataType d ) {
    this->get_color() += d;

    if( AlphaCount > 0 ) {
        this->get_alpha() += d;
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator-=( DataType d ) {
    this->get_color() -= d;

    if( AlphaCount > 0 ) {
        this->get_alpha() -= d;
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator*=( DataType d ) {
    this->get_color() *= d;

    if( AlphaCount > 0 ) {
        this->get_alpha() *= d;
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType& color_base<PixelType, DataType, ColorCount, AlphaCount>::operator/=( DataType d ) {
    this->get_color() /= d;

    if( AlphaCount > 0 ) {
        this->get_alpha() /= d;
    }

    return *(PixelType*)this;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::get_abs() const {
    PixelType result;

    result.get_color() = this->get_color().get_abs();

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha().get_abs();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::get_sqrt() const {
    PixelType result;

    result.get_color() = this->get_color().get_sqrt();

    if( AlphaCount > 0 ) {
        result.get_alpha() = this->get_alpha().get_sqrt();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
bool color_base<PixelType, DataType, ColorCount, AlphaCount>::operator==( const PixelType& p ) const {
    bool resultC = true;
    bool resultA = true;

    if( ColorCount > 0 ) {
        resultC = this->get_color() == p.get_color();
    }
    if( AlphaCount > 0 ) {
        resultA = this->get_alpha() == p.get_alpha();
    }

    return resultC && resultA;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
bool color_base<PixelType, DataType, ColorCount, AlphaCount>::operator!=( const PixelType& p ) const {
    return !( ( *this ) == p );
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::apply_alpha_to_color() const {
    PixelType result;

    if( AlphaCount == 0 ) {
        result.get_color() = this->get_color();
    } else {
        if( AlphaCount == 1 ) {
            result.get_color() = this->get_color() * this->get_alpha().get_data( 0 );
        } else if( ColorCount == AlphaCount ) {
            for( int i = 0; i < ColorCount; i++ ) {
                result.get_color().set_data( i, this->get_color().get_data( i ) * this->get_alpha().get_data( i ) );
            }
        }
        result.get_alpha() = this->get_alpha();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::remove_alpha_from_color() const {
    PixelType result;

    if( AlphaCount == 0 ) {
        result.get_color() = this->get_color();
    } else {
        if( AlphaCount == 1 ) {
            result.get_color() = this->get_color() / this->get_alpha().get_data( 0 );
        } else if( ColorCount == AlphaCount ) {
            for( int i = 0; i < ColorCount; i++ ) {
                result.get_color().set_data( i, this->get_color().get_data( i ) / this->get_alpha().get_data( i ) );
            }
        }
        result.get_alpha() = this->get_alpha();
    }

    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::blend_over( const PixelType& p ) const {
    PixelType result = *(PixelType*)this;
    result.apply_blend_over( p );
    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
void color_base<PixelType, DataType, ColorCount, AlphaCount>::apply_blend_over( const PixelType& p ) {
    if( AlphaCount == 0 ) {
        get_color() = this->get_color();
    } else {
        if( AlphaCount == 1 ) {
            this->get_color() = this->get_color() + p.get_color() - p.get_color() * this->get_alpha().get_data( 0 );
        } else if( ColorCount == AlphaCount ) {
            for( int i = 0; i < ColorCount; i++ ) {
                this->get_color().set_data( i, this->get_color().get_data( i ) + p.get_color().get_data( i ) -
                                                   p.get_color().get_data( i ) * this->get_alpha().get_data( i ) );
            }
        }
        this->get_alpha() = this->get_alpha() + p.get_alpha() - this->get_alpha() * p.get_alpha();
    }
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
PixelType color_base<PixelType, DataType, ColorCount, AlphaCount>::blend_under( const PixelType& p ) const {
    PixelType result = *(PixelType*)this;
    result.apply_blend_under( p );
    return result;
}

template <class PixelType, typename DataType, int ColorCount, int AlphaCount>
void color_base<PixelType, DataType, ColorCount, AlphaCount>::apply_blend_under( const PixelType& p ) {
    if( AlphaCount == 0 ) {
        this->get_color() = p.get_color();
    } else {
        if( AlphaCount == 1 ) {
            this->get_color() = this->get_color() + p.get_color() - this->get_color() * p.get_alpha().get_data( 0 );
        } else if( ColorCount == AlphaCount ) {
            for( int i = 0; i < ColorCount; i++ ) {
                this->get_color().set_data( i, this->get_color().get_data( i ) + p.get_color().get_data( i ) -
                                                   this->get_color().get_data( i ) * p.get_alpha().get_data( i ) );
            }
        }
        this->get_alpha() = this->get_alpha() + p.get_alpha() - this->get_alpha() * p.get_alpha();
    }
}

#if defined( _MSC_VER )
#pragma warning( pop )
#endif
} // namespace graphics
} // namespace frantic
