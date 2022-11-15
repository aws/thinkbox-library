// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "frantic/graphics/color_converter.hpp"

#include <iostream>

#include "frantic/graphics/color_rgb3a_f.hpp"
#include "frantic/graphics/color_rgb_f.hpp"
#include "frantic/graphics/color_rgba_f.hpp"
#include "frantic/graphics/color_rgba_h.hpp"

using namespace std;
using namespace frantic::graphics;
using namespace frantic::channels;

namespace frantic {
namespace graphics {
namespace details {

/** The one step pixel conversion is the simplest conversion that
 * can take place with two pixel objects.  This class should include
 * functions for converting between all pixel formats using a FLOAT
 * data type, as well as conversion functions from non-FLOAT to FLOAT
 * data types.  This will help ensure that all pixel types can be handled
 * without an exponential explosion of conversion functions.
 */
class convert_color_1_step {
  public:
    /*-------------------------------------------------------------------------------------
     * Functions that operate on an RGB_F source pixel type
     *-------------------------------------------------------------------------------------*/
    /** A null conversion from FLOAT RGB to FLOAT RGB.
     */
    static void convert_rgb_f_to_rgb_f( void* s, void* d ) {
        // cout<< "RGB_F to RGB_F" << endl;
        *(color_rgb_f*)d = *(color_rgb_f*)s;
    }

    /** A conversion from FLOAT RGB to FLOAT RGBA. The alpha channel is untouched.
     */
    static void convert_rgb_f_to_rgba_f( void* s, void* d ) {
        color_rgb_f& pix = *(color_rgb_f*)s;
        ( (color_rgba_f*)d )->set_r( pix.get_r() );
        ( (color_rgba_f*)d )->set_g( pix.get_g() );
        ( (color_rgba_f*)d )->set_b( pix.get_b() );
        ( (color_rgba_f*)d )->set_a( 1.0f );
    }

    /** A conversion from FLOAT RGB to HALF RGBA. The alpha channel is untouched.
     */
    static void convert_rgb_f_to_rgba_h( void* s, void* d ) {
        color_rgb_f& pix = *(color_rgb_f*)s;
        ( (color_rgba_h*)d )->set_r( (half)pix.get_r() );
        ( (color_rgba_h*)d )->set_g( (half)pix.get_g() );
        ( (color_rgba_h*)d )->set_b( (half)pix.get_b() );
        ( (color_rgba_h*)d )->set_a( (half)1.0f );
    }

    /** A conversion from FLOAT RGB to FLOAT RGB3A. The alpha channel is untouched.
     */
    static void convert_rgb_f_to_rgb3a_f( void* s, void* d ) {
        color_rgb_f& pix = *(color_rgb_f*)s;
        ( (color_rgb3a_f*)d )->set_r( pix.get_r() );
        ( (color_rgb3a_f*)d )->set_g( pix.get_g() );
        ( (color_rgb3a_f*)d )->set_b( pix.get_b() );
        ( (color_rgb3a_f*)d )->set_ar( 1.0f );
        ( (color_rgb3a_f*)d )->set_ag( 1.0f );
        ( (color_rgb3a_f*)d )->set_ab( 1.0f );
    }

    /*-------------------------------------------------------------------------------------
     * Functions that operate on an RGBA_F source pixel type
     *-------------------------------------------------------------------------------------*/
    /** A conversion from FLOAT RGBA to FLOAT RGB.  The alpha channel is dropped.
     */
    static void convert_rgba_f_to_rgb_f( void* s, void* d ) {
        color_rgba_f pix = ( (color_rgba_f*)s )->remove_alpha_from_color();
        *(color_rgb_f*)d = color_rgb_f( pix.get_r(), pix.get_g(), pix.get_b() );
    }

    /** A null conversion from FLOAT RGBA to FLOAT RGBA.
     */
    static void convert_rgba_f_to_rgba_f( void* s, void* d ) { *(color_rgba_f*)d = *(color_rgba_f*)s; }

    /** A datatype conversion from FLOAT RGBA to HALF RGBA. The 32-bit float is converted to a 16-bit float.
     */
    static void convert_rgba_f_to_rgba_h( void* s, void* d ) {
        color_rgba_f& pix = *(color_rgba_f*)s;
        *(color_rgba_h*)d = color_rgba_h( (half)pix.get_r(), (half)pix.get_g(), (half)pix.get_b(), (half)pix.get_a() );
    }

    /** A conversion from FLOAT RGBA to FLOAT RGB3A. All alpha channels are set to value of the single source alpha
     * channel.
     */
    static void convert_rgba_f_to_rgb3a_f( void* s, void* d ) {
        color_rgba_f& pix = *(color_rgba_f*)s;
        *(color_rgb3a_f*)d =
            color_rgb3a_f( pix.get_r(), pix.get_g(), pix.get_b(), pix.get_a(), pix.get_a(), pix.get_a() );
    }

    /*-------------------------------------------------------------------------------------
     * Functions that operate on an RGBA_H source pixel type
     *-------------------------------------------------------------------------------------*/
    /** A datatype conversion from HALF RGBA to FLOAT RGBA. The 16-bit float is converted to 32-bit float.
     */
    static void convert_rgba_h_to_rgba_f( void* s, void* d ) {
        color_rgba_h& pix = *(color_rgba_h*)s;
        *(color_rgba_f*)d =
            color_rgba_f( (float)pix.get_r(), (float)pix.get_g(), (float)pix.get_b(), (float)pix.get_a() );
    }

    /** A null conversion from HALF RGBA to HALF RGBA.
     */
    static void convert_rgba_h_to_rgba_h( void* s, void* d ) { *(color_rgba_h*)d = *(color_rgba_h*)s; }

    /*-------------------------------------------------------------------------------------
     * Functions that operate on an RGB3A_F source pixel type
     *-------------------------------------------------------------------------------------*/
    /** A conversion from FLOAT RGB3A to FLOAT RGB.  All alpha channels are dropped.
     */
    static void convert_rgb3a_f_to_rgb_f( void* s, void* d ) {
        color_rgb3a_f pix = ( (color_rgb3a_f*)s )->remove_alpha_from_color();
        *(color_rgb_f*)d = color_rgb_f( pix.get_r(), pix.get_g(), pix.get_b() );
    }

    /** A conversion from FLOAT RGB3A to FLOAT RGBA. The individual alphas are removed from the pixel and
     * the mean of the 3 alpha channels is applied and stored.
     */
    static void convert_rgb3a_f_to_rgba_f( void* s, void* d ) {
        color_rgb3a_f& pix = *(color_rgb3a_f*)s;
        float alpha = ( pix.get_ar() + pix.get_ag() + pix.get_ab() ) / 3.0f;
        *(color_rgba_f*)d = color_rgba_f( alpha * pix.get_r() / pix.get_ar(), alpha * pix.get_g() / pix.get_ag(),
                                          alpha * pix.get_b() / pix.get_ab(), alpha );
    }

    /** A null conversion from FLOAT RGB3A to FLOAT RGB3A.
     */
    static void convert_rgb3a_f_to_rgb3a_f( void* s, void* d ) { *(color_rgb3a_f*)d = *(color_rgb3a_f*)s; }
};

/** This class makes it possible to convert from one pixel format to another
 * with an intermediate step in the middle.  This is used when either the
 * source or destionation pixel is a FLOAT.  It contains only one static
 * template function.
 */
template <class PixelTypeA, void ( *StoA )( void*, void* ), void ( *AtoD )( void*, void* )>
class convert_color_2_step {
  public:
    static void convert_s_to_d( void* s, void* d ) {
        PixelTypeA a;
        StoA( s, &a );
        AtoD( &a, d );
    }
};

/** This class makes it possible to convert from one pixel format to another
 * with an intermediate step in the middle.  This is used when either the
 * source or destionation pixel is a FLOAT.  It contains only one static
 * template function.
 */
template <class PixelTypeA, class PixelTypeB, void ( *StoA )( void*, void* ), void ( *AtoB )( void*, void* ),
          void ( *BtoD )( void*, void* )>
class convert_color_3_step {
  public:
    static void convert_s_to_d( void* s, void* d ) {
        PixelTypeA a;
        PixelTypeA b;
        StoA( s, &a );
        AtoB( &a, &b );
        BtoD( &b, d );
    }
};

} // namespace details
} // namespace graphics
} // namespace frantic

using namespace frantic::graphics::details;

/** Selects the correct conversion function to transform one pixel format to another.  The can involve a 1, 2, or 3
 * step conversion process depending on the source and destination pixel.  However, this entire process is hidden
 * from the user inside of the frantic::channels::generic_channel_buffer_iteratator.
 *
 * @note		Any pixel format that is added to the library must be handled by this function or the
 *				frantic::channels::generic_channel_buffer_iteratator will throw an exception
 *
 * @arg	s	The name of the source pixel format
 * @arg	d	The name of the destination pixel format
 *	@return		Returns a void function that accepts two void pointers, and will convert the source pixel to the
 *				pixel destination, or 0 if no conversion function is defined.
 */
void ( *frantic::graphics::color_converter::get_conversion_function( const frantic::tstring& s, data_type_t sType,
                                                                     const frantic::tstring& d,
                                                                     data_type_t dType ) )( void*, void* ) {
    // Source type RGB Float
    if( s.compare( color_rgb_f::get_name() ) == 0 && sType == data_type_float32 ) {
        // Destination type RGB Float
        if( d.compare( color_rgb_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgb_f_to_rgb_f;
        }
        // Destination type RGBA Float
        else if( d.compare( color_rgba_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgb_f_to_rgba_f;
        }
        // Destination type RGBA Half
        else if( d.compare( color_rgba_h::get_name() ) == 0 && dType == data_type_float16 ) {
            return convert_color_1_step::convert_rgb_f_to_rgba_h;
        }
        // Destination type RGB3A Float
        else if( d.compare( color_rgb3a_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgb_f_to_rgb3a_f;
        }
        // Source type RGBA Float
    } else if( s.compare( color_rgba_f::get_name() ) == 0 && sType == data_type_float32 ) {
        // Destination type RGB Float
        if( d.compare( color_rgb_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgba_f_to_rgb_f;
        }
        // Destination type RGBA Float
        else if( d.compare( color_rgba_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgba_f_to_rgba_f;
        }
        // Destination type RGBA Half
        else if( d.compare( color_rgba_h::get_name() ) == 0 && dType == data_type_float16 ) {
            return convert_color_1_step::convert_rgba_f_to_rgba_h;
        }
        // Destination type RGB3A Float
        else if( d.compare( color_rgb3a_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgba_f_to_rgb3a_f;
        }
        // Source type RGBA Half
    } else if( s.compare( color_rgba_h::get_name() ) == 0 && sType == data_type_float16 ) {
        // Destination type RGB Float
        if( d.compare( color_rgb_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_2_step<color_rgba_f, convert_color_1_step::convert_rgba_h_to_rgba_f,
                                        convert_color_1_step::convert_rgba_f_to_rgb_f>::convert_s_to_d;
        }
        // Destination type RGBA Float
        else if( d.compare( color_rgba_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgba_h_to_rgba_f;
        }
        // Destination type RGBA Half
        else if( d.compare( color_rgba_h::get_name() ) == 0 && dType == data_type_float16 ) {
            return convert_color_1_step::convert_rgba_h_to_rgba_h;
        }
        // Destination type RGB3A Float
        else if( d.compare( color_rgb3a_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_2_step<color_rgba_f, convert_color_1_step::convert_rgba_h_to_rgba_f,
                                        convert_color_1_step::convert_rgba_f_to_rgb3a_f>::convert_s_to_d;
        }
        // Source type RGB3A Float
    } else if( s.compare( color_rgb3a_f::get_name() ) == 0 && sType == data_type_float32 ) {
        // Destination type RGB Float
        if( d.compare( color_rgb_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgb3a_f_to_rgb_f;
        }
        // Destination type RGBA Float
        else if( d.compare( color_rgba_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgb3a_f_to_rgba_f;
        }
        // Destination type RGBA Half
        else if( d.compare( color_rgba_h::get_name() ) == 0 && dType == data_type_float16 ) {
            return convert_color_2_step<color_rgba_f, convert_color_1_step::convert_rgb3a_f_to_rgba_f,
                                        convert_color_1_step::convert_rgba_f_to_rgba_h>::convert_s_to_d;
        }
        // Destination type RGB3A Float
        else if( d.compare( color_rgb3a_f::get_name() ) == 0 && dType == data_type_float32 ) {
            return convert_color_1_step::convert_rgb3a_f_to_rgb3a_f;
        }
    }

    return 0;
}
