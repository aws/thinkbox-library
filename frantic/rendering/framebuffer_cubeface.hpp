// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/files/files.hpp>
#include <frantic/graphics/spherical_coords.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics2d/framebuffer.hpp>

#include <frantic/graphics/alpha3f.hpp>
#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/color_with_alpha.hpp>

namespace frantic {
namespace rendering {

inline frantic::graphics::alpha3f cubeface_compensation( const frantic::graphics::alpha3f& a, const float factor ) {
    return frantic::graphics::alpha3f( 1 - std::pow( 1 - a.ar, factor ), 1 - std::pow( 1 - a.ag, factor ),
                                       1 - std::pow( 1 - a.ab, factor ) );
}

inline frantic::graphics::alpha1f cubeface_compensation( const frantic::graphics::alpha1f a, const float factor ) {
    return frantic::graphics::alpha1f( 1 - std::pow( 1 - a.a, factor ) );
}

inline frantic::graphics::color3f cubeface_compensation( const frantic::graphics::color3f& c, const float factor ) {
    return c * factor;
}

inline frantic::graphics::color4f cubeface_compensation( const frantic::graphics::color4f& c, const float factor ) {
    frantic::graphics::alpha1f newAlpha = cubeface_compensation( c.a, factor );
    return frantic::graphics::color4f( ( newAlpha.a / c.a.a ) * c.c, newAlpha );
}

inline frantic::graphics::color6f cubeface_compensation( const frantic::graphics::color6f& c, const float factor ) {
    frantic::graphics::alpha3f newAlpha = cubeface_compensation( c.a, factor );
    frantic::graphics::color3f newColor( ( newAlpha.ar / c.a.ar ) * c.c.r, ( newAlpha.ag / c.a.ag ) * c.c.g,
                                         ( newAlpha.ab / c.a.ab ) * c.c.b );
    return frantic::graphics::color6f( newColor, newAlpha );
}

template <class ColorType>
class framebuffer_cubeface /* : public environment_map_provider<ColorType>*/ {
  public:
    // A member-function typedef for the get pixel function
    typedef ColorType ( framebuffer_cubeface<ColorType>::*get_pixel_function_type )(
        const frantic::graphics::vector3f& ) const;

  private:
    frantic::graphics2d::framebuffer<ColorType> m_cubefaces[6];
    int m_size;

    // This is the default mode for pixel lookup
    get_pixel_function_type m_getPixelFunction;

  public:
    framebuffer_cubeface() {
        set_size( 1 );
        m_getPixelFunction = &framebuffer_cubeface<ColorType>::get_pixel_bilinear;
    }

    framebuffer_cubeface( int size ) {
        set_size( size );
        m_getPixelFunction = &framebuffer_cubeface<ColorType>::get_pixel_bilinear;
    }

    virtual ~framebuffer_cubeface() {}

    void set_size( int size ) {
        m_size = size;
        for( int i = 0; i < 6; ++i )
            m_cubefaces[i].set_size( frantic::graphics2d::size2( m_size, m_size ) );
    }

    void set_size( frantic::graphics2d::size2 size ) {
        if( size.xsize != size.ysize )
            throw std::runtime_error( "framebuffer_cubeface<>.set_size() The provided size was not square." );
        set_size( size.xsize );
    }

    int size() const { return m_size; }

    int width() const { return m_size; }

    int height() const { return m_size; }

    void set_pixel_lookup_filter( frantic::graphics2d::pixel_lookup_filter::pixel_lookup_filter_enum lookupFilter ) {
        switch( lookupFilter ) {
        case frantic::graphics2d::pixel_lookup_filter::nearest_neighbor_filter:
            m_getPixelFunction = &framebuffer_cubeface<ColorType>::get_pixel_nearest_neighbor;
            break;
        case frantic::graphics2d::pixel_lookup_filter::bilinear_filter:
            m_getPixelFunction = &framebuffer_cubeface<ColorType>::get_pixel_bilinear;
            break;
        case frantic::graphics2d::pixel_lookup_filter::bicubic_filter:
            m_getPixelFunction = &framebuffer_cubeface<ColorType>::get_pixel_bicubic;
            break;
        default:
            throw std::runtime_error(
                "framebuffer_cubeface.set_pixel_lookup_filter: Invalid pixel lookup filter specified." );
        }
    }

    frantic::graphics2d::pixel_lookup_filter::pixel_lookup_filter_enum get_pixel_lookup_filter() const {
        if( m_getPixelFunction == &framebuffer_cubeface<ColorType>::get_pixel_nearest_neighbor )
            return frantic::graphics2d::pixel_lookup_filter::nearest_neighbor_filter;
        else if( m_getPixelFunction == &framebuffer_cubeface<ColorType>::get_pixel_bilinear )
            return frantic::graphics2d::pixel_lookup_filter::bilinear_filter;
        else if( m_getPixelFunction == &framebuffer_cubeface<ColorType>::get_pixel_bicubic )
            return frantic::graphics2d::pixel_lookup_filter::bicubic_filter;
        else
            throw std::runtime_error( "framebuffer_cubeface.get_pixel_lookup_filter: The pixel lookup filter was in an "
                                      "invalid state, this is likely a bug in the software somewhere." );
    }

    void swap( framebuffer_cubeface<ColorType>& rhs ) {
        for( int i = 0; i < 6; ++i )
            m_cubefaces[i].swap( rhs.m_cubefaces[i] );
        std::swap( m_size, rhs.m_size );
    }

    void clear() {
        for( int i = 0; i < 6; ++i )
            m_cubefaces[i].clear();
        m_size = 0;
    }

    frantic::graphics2d::framebuffer<ColorType>& get_cubeface_framebuffer( int cubeFace ) {
        return m_cubefaces[cubeFace];
    }

    const frantic::graphics2d::framebuffer<ColorType>& get_cubeface_framebuffer( int cubeFace ) const {
        return m_cubefaces[cubeFace];
    }

    ColorType& at( int cubeFace, frantic::graphics2d::vector2 pixelCoord ) {
        if( cubeFace < 0 || cubeFace >= 6 ) {
            throw std::runtime_error( "framebuffer_cube.at: Tried to access out of bounds cube face " +
                                      boost::lexical_cast<std::string>( cubeFace ) );
        }
        return m_cubefaces[cubeFace][pixelCoord];
    }

    // This get pixel function calls a member function, which will be one of the other get pixel functions in the class
    // depending on the pixel lookup filter setting.
    ColorType get_pixel( const frantic::graphics::vector3f& direction ) const {
        return ( this->*m_getPixelFunction )( direction );
    }

    // Copy of above function to match framebuffer<> function name
    ColorType get_pixel_filtered( const frantic::graphics::vector3f& direction ) const {
        return ( this->*m_getPixelFunction )( direction );
    }

    ColorType get_pixel_nearest_neighbor( const frantic::graphics::vector3f& direction ) const {
        frantic::graphics::cube_face::default_cube_face cubeFace = get_cube_face( direction );
        frantic::graphics2d::vector2f coord = get_cube_face_coordinate( direction, cubeFace );
        return m_cubefaces[cubeFace].get_pixel_nearest_neighbor( coord );
    }

    ColorType get_pixel_bilinear( const frantic::graphics::vector3f& direction ) const {
        // Get the cube faces within a tolerance of 0.5 of a pixel.
        float ratioTolerance = ( m_size - 1.f ) / m_size;
        frantic::graphics::cube_face::default_cube_face cubeFaces[6];
        // This gets all the cube faces that are within 2 pixels of the direction
        int cubeFaceCount = get_cube_faces( direction, ratioTolerance, cubeFaces );

        ColorType result = ColorType();
        float weightAccumulator = 0;
        for( int i = 0; i < cubeFaceCount; ++i ) {
            result += m_cubefaces[cubeFaces[i]].get_pixel_bilinear_return_weight(
                get_cube_face_coordinate( direction, cubeFaces[i] ), weightAccumulator );
        }
        if( weightAccumulator != 1.f ) {
            result /= weightAccumulator;
        }
        return result;
    }

    ColorType get_pixel_bicubic( const frantic::graphics::vector3f& direction ) const {
        // Get the cube faces within a tolerance of 2 of a pixel.
        float ratioTolerance = ( m_size - 4.f ) / m_size;
        frantic::graphics::cube_face::default_cube_face cubeFaces[6];

        int cubeFaceCount = get_cube_faces( direction, ratioTolerance, cubeFaces );

        ColorType result = ColorType();
        float weightAccumulator = 0;
        for( int i = 0; i < cubeFaceCount; ++i ) {
            result += m_cubefaces[cubeFaces[i]].get_pixel_bicubic_return_weight(
                get_cube_face_coordinate( direction, cubeFaces[i] ), weightAccumulator );
        }
        if( weightAccumulator != 1.f ) {
            result /= weightAccumulator;
        }

        return result;
    }

    ColorType lookup_environment( const frantic::graphics::vector3f& direction ) const {
        return get_pixel_bilinear( direction );
    }

    frantic::graphics2d::draw_point_filter::draw_point_filter_enum get_draw_point_filter() const {
        return m_cubefaces[0].get_draw_point_filter();
    }

    void set_draw_point_filter( frantic::graphics2d::draw_point_filter::draw_point_filter_enum filter ) {
        for( int i = 0; i < 6; ++i ) {
            m_cubefaces[i].set_draw_point_filter( filter );
        }
    }

    // The area of a voxel on screen at a distance d from the camera is (w/(2*d*tan(fov/2))^2.  In the case of a 90
    // degree field of view, the tangent term is 1.  This value is calculated at the center of the cube faces.
    float draw_point_scaling_constant() const { return m_size * m_size / 4.f; }

    // The function draw_point is guaranteed to add exactly the amount of energy specified in the color to the image,
    // normalized to behave as if it is drawn on the center of one of the cube faces.
    void draw_point( frantic::graphics::vector3f direction, const ColorType& color ) {
        using namespace std;

        //		cerr << "Drawing point " << direction << endl;
        // The radius of the framebuffer::draw_point is 2, so this is the maximum radius of effect
        float ratioTolerance = ( m_size - 4.f ) / m_size;
        frantic::graphics::cube_face::default_cube_face cubeFaces[6];
        // This gets all the cube faces that are within 2 pixels of the direction
        int cubeFacesCount = get_cube_faces( direction, ratioTolerance, cubeFaces );

        // Can safely assume the size of the array is at least 1
        frantic::graphics2d::vector2f coordinate0 = get_cube_face_coordinate( direction, cubeFaces[0] );

        // The compensation factor is cos^3(theta).  This calculation depends on the fact that the face has a 90 degree
        // field of view.
        //
        // The reason for the cos^3(theta) is twofold.  First, consider the ratio of the distance to the point being
        // drawn versus the perpendicular distance. If you draw out the little triangle, it becomes clear that this
        // distance is cos(theta).  Because the inverse square falloff behavior of point sources, this means that the
        // ray density at that point is decreased by 1/cos^2(theta). Next, consider the angle at which the rays are
        // hitting the surface. If you compute the ratio of area of the stretched surface versus the area of he
        // perpendicular surface, you get that the surface area is increased by cos(theta).  Multiplying these two
        // factors together as a compensation for the decrease in density gives us cos^3(theta).
        //
        // This compensation works correctly for both additive drawing of particles and alpha blending of particles.
        float compensationFactor = ( ( coordinate0.x * coordinate0.x + coordinate0.y * coordinate0.y + 1 ) );
        compensationFactor *= sqrt( compensationFactor );

        ColorType compensatedColor = cubeface_compensation( color, compensationFactor );

        // draw on the first face
        // m_cubefaces[cubeFaces[0]].draw_point( (0.5f * m_size) * (coordinate0 + frantic::graphics2d::vector2f(1,1)),
        // compensationFactor * color );
        m_cubefaces[cubeFaces[0]].draw_point(
            ( 0.5f * m_size ) * ( coordinate0 + frantic::graphics2d::vector2f( 1, 1 ) ), compensatedColor );

        // draw on the rest of the faces
        for( int i = 1; i < cubeFacesCount; ++i ) {
            //			cerr << "Drawing on cube face " << cubeFaces[i] << " with coordinates " <<
            // direction.get_cube_face_coordinate( cubeFaces[i] ) << endl; m_cubefaces[cubeFaces[i]].draw_point( (0.5f *
            // m_size) * (direction.get_cube_face_coordinate( cubeFaces[i] ) + frantic::graphics2d::vector2f(1,1)),
            // compensationFactor * color );
            m_cubefaces[cubeFaces[i]].draw_point(
                ( 0.5f * m_size ) *
                    ( get_cube_face_coordinate( direction, cubeFaces[i] ) + frantic::graphics2d::vector2f( 1, 1 ) ),
                compensatedColor );
        }
    }

    void set_pixel( int cubeFace, frantic::graphics2d::vector2 pixel, const ColorType& color ) {
        m_cubefaces[cubeFace].set_pixel( pixel, color );
    }

    void add_image_data( const framebuffer_cubeface& other ) {
        for( unsigned i = 0; i < 6; ++i )
            m_cubefaces[i].add_image_data( other.m_cubefaces[i] );
    }

    void apply_gain( float gain ) {
        for( unsigned i = 0; i < 6; ++i )
            m_cubefaces[i].apply_gain( gain );
    }

    void fill( const ColorType& color ) {
        for( unsigned i = 0; i < 6; ++i ) {
            m_cubefaces[i].fill( color );
        }
    }

    void fill_face( int cubeFace, const ColorType& color ) {
        if( cubeFace < 0 || cubeFace >= 6 )
            throw std::runtime_error( "framebuffer_cubeface.fill_face: Tried to access out of bounds cube face " +
                                      boost::lexical_cast<std::string>( cubeFace ) );
        m_cubefaces[cubeFace].fill( color );
    }

    void fill_face_boundary( int cubeFace, const ColorType& color ) {
        if( cubeFace < 0 || cubeFace >= 6 )
            throw std::runtime_error( "framebuffer_cubeface.fill_face: Tried to access out of bounds cube face " +
                                      boost::lexical_cast<std::string>( cubeFace ) );
        m_cubefaces[cubeFace].fill_boundary( color );
    }

    void fill_under( const ColorType& color ) {
        for( unsigned i = 0; i < 6; ++i )
            m_cubefaces[i].fill_under( color );
    }

    void blend_under( framebuffer_cubeface& over ) {
        for( int i = 0; i < 6; ++i )
            m_cubefaces[i].blend_under( over.m_cubefaces[i] );
    }

    void blend_under( int cubeFace, int x, int y, const ColorType& color ) {
        m_cubefaces[cubeFace].blend_over( x, y, color );
    }

    void blend_over( framebuffer_cubeface& under ) {
        for( int i = 0; i < 6; ++i )
            m_cubefaces[i].blend_over( under.m_cubefaces[i] );
    }

    void blend_over( int cubeFace, int x, int y, const ColorType& color ) {
        m_cubefaces[cubeFace].blend_over( x, y, color );
    }

    void create_borders( const ColorType& borderColor ) {
        for( int i = 0; i < 6; i++ ) {
            m_cubefaces[i].fill_boundary( borderColor );
        }
    }

    // Fills the target framebuffer with a longlat image of the cubeface framebuffer
    // In the target framebuffer, the top of the image will correspond to positive Z
    void
    to_longlat_zup( frantic::graphics2d::framebuffer<ColorType>& fb, int superSampling = 2,
                    const frantic::graphics::transform4f& xform = frantic::graphics::transform4f::identity() ) const {
        frantic::graphics2d::vector2 p;
        frantic::graphics2d::size2 sz = fb.size();
        frantic::graphics2d::size2 szSuper = superSampling * sz;
        for( p.y = 0; p.y < sz.ysize; ++p.y ) {
            for( p.x = 0; p.x < sz.xsize; ++p.x ) {
                ColorType color = ColorType();
                for( int sy = 0; sy < superSampling; ++sy ) {
                    for( int sx = 0; sx < superSampling; ++sx ) {
                        frantic::graphics::spherical_coords sc = frantic::graphics::spherical_coords::from_longlat(
                            1.f - ( superSampling * p.x + sx + 0.5f ) / szSuper.xsize,
                            1.f - ( superSampling * p.y + sy + 0.5f ) / szSuper.ysize );
                        color += get_pixel( xform.transform_no_translation( sc.to_vector3f() ) );
                    }
                }
                color *= 1.f / ( superSampling * superSampling );
                fb.set_pixel( p, color );
            }
        }
    }

    // Fills the target framebuffer with a longlat image of the cubeface framebuffer
    // In the target framebuffer, the top of the image will correspond to positive Y
    void to_longlat_yup( frantic::graphics2d::framebuffer<ColorType>& fb, int superSampling = 2 ) const {
        // Use the zup function with an appropriate transform matrix
        to_longlat_zup( fb, superSampling,
                        frantic::graphics::transform4f( frantic::graphics::vector3f( 0, 0, 1 ),
                                                        frantic::graphics::vector3f( 1, 0, 0 ),
                                                        frantic::graphics::vector3f( 0, 1, 0 ) ) );
    }

    void from_longlat_zup( const frantic::graphics2d::framebuffer<ColorType>& fb, float latitudeBottom = 0,
                           float latitudeTop = 1, int superSampling = 2 ) {
        latitudeBottom = frantic::math::clamp( latitudeBottom, 0.f, 1.f );
        latitudeTop = frantic::math::clamp( latitudeTop, latitudeBottom, 1.f );

        fb.to_OpenEXR_file( "from_longlat_zup.exr" );

        frantic::graphics2d::vector2 p;
        int sz = m_size;
        int szSuper = superSampling * sz;
        for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
            frantic::graphics2d::framebuffer<ColorType>& cubeFaceBuffer = m_cubefaces[cubeFace];
            for( p.y = 0; p.y < sz; ++p.y ) {
                for( p.x = 0; p.x < sz; ++p.x ) {
                    ColorType color;
                    for( int sy = 0; sy < superSampling; ++sy ) {
                        for( int sx = 0; sx < superSampling; ++sx ) {
                            frantic::graphics2d::vector2f cubeCoord( ( superSampling * p.x + sx + 0.5f ) / szSuper,
                                                                     ( superSampling * p.y + sy + 0.5f ) / szSuper );
                            frantic::graphics::vector3f direction = frantic::graphics::from_cube_face_coordinate(
                                cubeCoord, (frantic::graphics::cube_face::default_cube_face)cubeFace );

                            frantic::graphics::spherical_coords sc( direction );
                            frantic::graphics2d::vector2f longLatCoord = sc.to_longlat();

                            // std::cout << "Face " << cubeFace << ", coordinate " << p << ", direction " << direction
                            // << ", longLatCoord " << longLatCoord << std::endl;

                            longLatCoord.y -= latitudeBottom;
                            longLatCoord.y /= ( latitudeTop - latitudeBottom );
                            if( longLatCoord.y >= 0 && longLatCoord.y <= 1 ) {
                                longLatCoord.x = 1 - longLatCoord.x * 2;
                                longLatCoord.y = 1 - longLatCoord.y * 2;
                                // std::cout << "Getting color at " << longLatCoord << std::endl;
                                color += fb.get_pixel_bilinear_wrap( longLatCoord );
                            }
                        }
                    }
                    color *= 1.f / ( superSampling * superSampling );
                    cubeFaceBuffer.set_pixel( p, color );
                }
            }
        }
    }

    void to_longlat_OpenEXR_file_zup(
        const frantic::tstring& filename, int width = -1, int height = -1, int superSampling = 2,
        const frantic::graphics::transform4f& xform = frantic::graphics::transform4f::identity() ) {
        if( width < 0 )
            width = m_size * 2;
        if( height < 0 )
            height = m_size;
        frantic::graphics2d::framebuffer<ColorType> fbLongLat( frantic::graphics2d::size2( width, height ) );
        to_longlat_zup( fbLongLat, superSampling, xform );
        fbLongLat.to_OpenEXR_file( filename );
    }

    void to_longlat_OpenEXR_file_yup( const frantic::tstring& filename, int width = -1, int height = -1,
                                      int superSampling = 2 ) {
        to_longlat_OpenEXR_file_zup( filename, width, height, superSampling,
                                     frantic::graphics::transform4f( frantic::graphics::vector3f( 0, 0, 1 ),
                                                                     frantic::graphics::vector3f( 1, 0, 0 ),
                                                                     frantic::graphics::vector3f( 0, 1, 0 ) ) );
    }

    void from_cubeface_files( const std::string& filename ) {
        std::string basename = frantic::files::basename_from_path( filename );
        if( basename.size() < 5 )
            throw std::runtime_error( "framebuffer_cubeface.from_cubeface_files:  The filename \"" + filename +
                                      "\" doesn't have a cube face specifier at the end (for example _xpos)." );

        static const char* specifiers[6] = { "_xpos", "_xneg", "_ypos", "_yneg", "_zpos", "_zneg" };
        // TODO: On unix we shouldn't do the to_lower call.
        std::string specifier = frantic::strings::to_lower( basename.substr( basename.size() - 5 ) );
        bool foundValidSpecifier = false;
        for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
            if( specifier == specifiers[cubeFace] )
                foundValidSpecifier = true;
        }
        if( !foundValidSpecifier )
            throw std::runtime_error( "framebuffer_cubeface.from_cubeface_files:  The filename \"" + filename +
                                      "\" doesn't have a cube face specifier at the end (for example _xpos)." );

        std::string filePath =
            frantic::files::ensure_trailing_pathseparator( frantic::files::directory_from_path( filename ) );
        std::string prefix = basename.substr( 0, basename.size() - 5 );
        std::string extension = frantic::files::extension_from_path( filename );
        for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
            try {
                std::string cubeFaceFile = filePath + prefix + specifiers[cubeFace] + extension;

                m_cubefaces[cubeFace].from_OpenEXR_file( cubeFaceFile );
            } catch( const std::exception& ) {
                clear();
                throw;
            }
        }

        // Make sure that the loaded images are square and all of the same size
        m_size = m_cubefaces[0].xsize();
        for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
            if( m_cubefaces[cubeFace].xsize() != m_size || m_cubefaces[cubeFace].ysize() != m_size ) {
                clear();
                throw std::runtime_error( "framebuffer_cubeface.from_cubeface_files:  The filename \"" + filename +
                                          "\" doesn't have a cube face specifier at the end (for example _xpos)." );
            }
        }
    }

    // Writes the cube faces to the framebuffer, in a 2x3 pattern (vertical or horizontal depending on the 'outImage'
    // dimensions).
    void to_framebuffer( frantic::graphics2d::framebuffer<ColorType>& outImage ) const {
        bool tall = false;
        if( outImage.width() < outImage.height() ) {
            tall = true;
            if( outImage.height() < 3 * m_size )
                throw std::runtime_error( "framebuffer_cubeface.to_framebuffer:  The supplied image was too small" );
        } else if( outImage.width() < 3 * m_size ) {
            throw std::runtime_error( "framebuffer_cubeface.to_framebuffer:  The supplied image was too small" );
        }

        if( tall ) {
            for( int i = 0; i < 3; ++i ) {
                int startY = m_size * i;

                for( int k = 0; k < m_size; ++k ) {
                    for( int j = 0; j < m_size; ++j )
                        outImage.set_pixel( j, startY + k, m_cubefaces[2 * i].get_pixel( j, k ) );
                    for( int j = 0; j < m_size; ++j )
                        outImage.set_pixel( j + m_size, startY + k, m_cubefaces[2 * i + 1].get_pixel( j, k ) );
                }
            }
        } else {
            for( int i = 0; i < 3; ++i ) {
                int startX = m_size * i;

                for( int k = 0; k < m_size; ++k ) {
                    for( int j = 0; j < m_size; ++j )
                        outImage.set_pixel( startX + j, k, m_cubefaces[2 * i].get_pixel( j, k ) );
                }
                for( int k = 0; k < m_size; ++k ) {
                    for( int j = 0; j < m_size; ++j )
                        outImage.set_pixel( startX + j, k + m_size, m_cubefaces[2 * i + 1].get_pixel( j, k ) );
                }
            }
        }
    }

    void to_cubeface_files( const std::string& filename ) const {
        std::string basename = frantic::files::basename_from_path( filename );

        static const char* specifiers[6] = { "_xpos", "_xneg", "_ypos", "_yneg", "_zpos", "_zneg" };

        bool appendSpecifier = true;
        if( basename.size() >= 5 ) {

            // TODO: On unix we shouldn't do the to_lower call.
            std::string specifier = frantic::strings::to_lower( basename.substr( basename.size() - 5 ) );
            for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
                if( specifier == specifiers[cubeFace] )
                    appendSpecifier = false;
            }
        }

        std::string filePath =
            frantic::files::ensure_trailing_pathseparator( frantic::files::directory_from_path( filename ) );
        std::string prefix = basename.substr( 0, basename.size() - ( appendSpecifier ? 0 : 5 ) );
        std::string extension = frantic::files::extension_from_path( filename );
        for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
            std::string cubeFaceFile = filePath + prefix + specifiers[cubeFace] + extension;

            m_cubefaces[cubeFace].to_OpenEXR_file( cubeFaceFile );
        }
    }

    // ==================================================================================
    // The following functions are used only for debugging purposes.
    // They fill each cubeface with predetermined images which can be
    // used as test patterns
    // ==================================================================================

    // Assigns a fill and boundary color to each cube face.
    void color_code_faces() {

        // Fills

        m_cubefaces[0].fill( ColorType( 1, 0, 0 ) ); // cube_face::CF_RIGHT		->	RED
        m_cubefaces[1].fill( ColorType( 0, 1, 1 ) ); // cube_face::CF_LEFT		->	TURQUOISE
        m_cubefaces[2].fill( ColorType( 0, 1, 0 ) ); // cube_face::CF_TOP		->	GREEN
        m_cubefaces[3].fill( ColorType( 1, 0, 1 ) ); // cube_face::CF_BOTTOM		->	PURPLE
        m_cubefaces[4].fill( ColorType( 0, 0, 1 ) ); // cube_face::CF_REAR		->	BLUE
        m_cubefaces[5].fill( ColorType( 1, 1, 0 ) ); // cube_face::CF_FRONT		->	YELLOW

        // Boundaries

        /*
        m_cubefaces[0].fill_boundary( ColorType(1,.5f,.5f) );			// cube_face::CF_RIGHT
        ->	LIGHT RED m_cubefaces[1].fill_boundary( ColorType(.5f,1,1) );		// cube_face::CF_LEFT
        ->	LIGHT TURQUOISE m_cubefaces[2].fill_boundary( ColorType(.5f,1,.5f) );			//
        cube_face::CF_TOP		->	LIGHT GREEN m_cubefaces[3].fill_boundary( ColorType(1,.5f,1) );
        // cube_face::CF_BOTTOM		->	LIGHT PURPLE m_cubefaces[4].fill_boundary( ColorType(.5f,.5f,1) );
        // cube_face::CF_REAR		->	LIGHT BLUE m_cubefaces[5].fill_boundary( ColorType(1,1,.5f) );
        // cube_face::CF_FRONT		->	LIGHT YELLOW
        */
    }

    // Fills each face with a gradient which sweeps from top left to bottom right
    // The top left corner of each face will be solid green
    // The bottom right corner of each face will be solid red
    void color_code_faces_gradient() {
        using namespace frantic::graphics::cube_face;

        // If the boolean is true, the face will be filled with the gradient
        // If the boolean is false, the face will be filled with black

        const bool FILL_RIGHT = true;
        const bool FILL_LEFT = true;
        const bool FILL_TOP = true;
        const bool FILL_BOTTOM = true;
        const bool FILL_REAR = true;
        const bool FILL_FRONT = true;

        const ColorType c1( 0, 1, 0 ); // Green
        const ColorType c2( 1, 0, 0 ); // Red
        const ColorType emptyColor;    // Black

        FILL_RIGHT ? m_cubefaces[CF_X_POS].fill_gradient( c1, c2 ) : m_cubefaces[CF_X_POS].fill( emptyColor );
        FILL_LEFT ? m_cubefaces[CF_X_NEG].fill_gradient( c1, c2 ) : m_cubefaces[CF_X_NEG].fill( emptyColor );
        FILL_TOP ? m_cubefaces[CF_Y_POS].fill_gradient( c1, c2 ) : m_cubefaces[CF_Y_POS].fill( emptyColor );
        FILL_BOTTOM ? m_cubefaces[CF_Y_NEG].fill_gradient( c1, c2 ) : m_cubefaces[CF_Y_NEG].fill( emptyColor );
        FILL_REAR ? m_cubefaces[CF_Z_POS].fill_gradient( c1, c2 ) : m_cubefaces[CF_Z_POS].fill( emptyColor );
        FILL_FRONT ? m_cubefaces[CF_Z_NEG].fill_gradient( c1, c2 ) : m_cubefaces[CF_Z_NEG].fill( emptyColor );
    }

    // This loads textures onto each cubeface's framebuffer with a corresponding
    // label for the cubeface.  It is only needed for debugging purposes and will
    // only be compiled in debug configuration
    void picture_code_faces() {
        using namespace frantic::graphics;

        frantic::graphics2d::framebuffer<ColorType> tempBuffer;
        frantic::graphics2d::size2 finalSize( m_size, m_size );

        tempBuffer.from_OpenEXR_file( "..\\TestFaces\\right.exr" );
        m_cubefaces[cube_face::CF_X_POS] = frantic::graphics2d::framebuffer<ColorType>( tempBuffer, finalSize );

        tempBuffer.from_OpenEXR_file( "..\\TestFaces\\left.exr" );
        m_cubefaces[cube_face::CF_X_NEG] = frantic::graphics2d::framebuffer<ColorType>( tempBuffer, finalSize );

        tempBuffer.from_OpenEXR_file( "..\\TestFaces\\top.exr" );
        m_cubefaces[cube_face::CF_Y_POS] = frantic::graphics2d::framebuffer<ColorType>( tempBuffer, finalSize );

        tempBuffer.from_OpenEXR_file( "..\\TestFaces\\bottom.exr" );
        m_cubefaces[cube_face::CF_Y_NEG] = frantic::graphics2d::framebuffer<ColorType>( tempBuffer, finalSize );

        tempBuffer.from_OpenEXR_file( "..\\TestFaces\\rear.exr" );
        m_cubefaces[cube_face::CF_Z_POS] = frantic::graphics2d::framebuffer<ColorType>( tempBuffer, finalSize );

        tempBuffer.from_OpenEXR_file( "..\\TestFaces\\front.exr" );
        m_cubefaces[cube_face::CF_Z_NEG] = frantic::graphics2d::framebuffer<ColorType>( tempBuffer, finalSize );
    }
};

} // namespace rendering
} // namespace frantic
