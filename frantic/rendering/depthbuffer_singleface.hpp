// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
//#include <frantic/rendering/framebuffer_cubeface.hpp>
#include <frantic/rendering/geometry_renderer.hpp>

#include <frantic/graphics2d/framebufferiterator.hpp>

#include <frantic/logging/console_progress_logger.hpp>

struct pick_min_float {
    float operator()( float lhs, float rhs ) { return lhs <= rhs ? lhs : rhs; }
};

namespace frantic {
namespace rendering {

class depthbuffer_singleface {
    transform4f m_transform, m_transformInverse;
    frantic::graphics2d::framebuffer<float> m_depthbuffer;

  public:
    depthbuffer_singleface() {
        m_depthbuffer.set_draw_point_filter( frantic::graphics2d::draw_point_filter::bilinear );
        m_depthbuffer.set_pixel_lookup_filter( frantic::graphics2d::pixel_lookup_filter::bilinear_filter );
    }

    depthbuffer_singleface( int mapSize )
        : m_depthbuffer( frantic::graphics2d::size2( mapSize, mapSize ) ) {
        m_depthbuffer.set_draw_point_filter( frantic::graphics2d::draw_point_filter::bilinear );
        m_depthbuffer.set_pixel_lookup_filter( frantic::graphics2d::pixel_lookup_filter::bilinear_filter );
    }

    depthbuffer_singleface( const frantic::graphics2d::size2& mapSize )
        : m_depthbuffer( mapSize ) {
        m_depthbuffer.set_draw_point_filter( frantic::graphics2d::draw_point_filter::bilinear );
        m_depthbuffer.set_pixel_lookup_filter( frantic::graphics2d::pixel_lookup_filter::bilinear_filter );
    }

    // Loads a depth image from file
    depthbuffer_singleface( const frantic::tstring& file ) {
        m_depthbuffer.set_draw_point_filter( frantic::graphics2d::draw_point_filter::bilinear );
        m_depthbuffer.set_pixel_lookup_filter( frantic::graphics2d::pixel_lookup_filter::nearest_neighbor_filter );
        m_depthbuffer.from_OpenEXR_file( file );
    }

    void clear() { m_depthbuffer.fill( std::numeric_limits<float>::max() ); }

    depthbuffer_singleface& operator=( const depthbuffer_singleface& rhs ) {
        m_depthbuffer = rhs.m_depthbuffer;
        return *this;
    }

    void set( const frantic::graphics2d::framebuffer<float>& depthBuffer ) { m_depthbuffer = depthBuffer; }

    void set_with_swap( frantic::graphics2d::framebuffer<float>& depthBuffer ) { m_depthbuffer.swap( depthBuffer ); }

    void swap( depthbuffer_singleface& rhs ) { m_depthbuffer.swap( rhs.m_depthbuffer ); }

    // Resizes the depth buffer
    void set_size( int mapSize ) { m_depthbuffer.set_size( frantic::graphics2d::size2( mapSize, mapSize ) ); }

    void set_size( const frantic::graphics2d::size2& mapSize ) { m_depthbuffer.set_size( mapSize ); }

    void set_transform( const frantic::graphics::transform4f& tm ) {
        m_transform = tm;
        m_transformInverse = m_transform.to_inverse();
    }

    frantic::graphics2d::size2 size() const { return m_depthbuffer.size(); }

    int width() const { return m_depthbuffer.width(); }

    int height() const { return m_depthbuffer.height(); }

    void combine( const frantic::graphics2d::framebuffer<float>& rhs ) {
        if( rhs.width() != width() || rhs.height() != height() )
            throw std::runtime_error(
                "depthbuffer_singleface::combine() - Input size: " + boost::lexical_cast<std::string>( rhs.size() ) +
                " does not match: " + boost::lexical_cast<std::string>( size() ) );

        std::transform( m_depthbuffer.data().begin(), m_depthbuffer.data().end(), rhs.data().begin(),
                        m_depthbuffer.data().begin(), pick_min_float() );
    }

    void combine( const depthbuffer_singleface& rhs ) { combine( rhs.m_depthbuffer ); }

    void resample_combine( const depthbuffer_singleface& rhs ) {
        // If the incoming image does not need to be resampled then just forward to the normal combine function.
        if( rhs.size() == this->size() )
            return this->combine( rhs );

        // We are going to use a nearest-neighbor magnification filter, and a non-overlapping 'minimum' minification
        // filter on rhs in order to combine it with the existing depth image. NOTE The magnification filter is used
        // when a dimension is increasing ( ex. rhs.width() < this->width() ) and the minification filter is applied
        // when
        //      the dimension is decreasing ( ex. rhs.width() > this->width() ). A min filter chooses the minimum value
        //      from its inputs.

        if( rhs.height() < this->height() ) {
            // We are expanding the image vertically
            if( rhs.width() < this->width() ) {
                // We are expanding the image horizontally too. Use the nearest-neighbor in rhs to get the value.
                for( int y = 0; y < this->height(); ++y ) {
                    int yRhs = static_cast<int>(
                        std::floor( static_cast<float>( rhs.height() ) *
                                    ( ( static_cast<float>( y ) + 0.5f ) / static_cast<float>( this->height() ) ) ) );
                    for( int x = 0; x < this->width(); ++x ) {
                        int xRhs = static_cast<int>(
                            std::floor( static_cast<float>( rhs.width() ) * ( ( static_cast<float>( x ) + 0.5f ) /
                                                                              static_cast<float>( this->width() ) ) ) );

                        float v = m_depthbuffer.get_pixel( x, y );
                        float vRhs = rhs.m_depthbuffer.get_pixel( xRhs, yRhs );
                        if( vRhs < v )
                            m_depthbuffer.set_pixel( x, y, vRhs );
                    }
                }
            } else {
                // We are expanding vertically, but shrinking (or maintaining the same size) horizontally. Use
                // nearest-neighbor on y, but a min filter on a range of x.

                boost::int64_t xRound =
                    rhs.width() +
                    this->width() / 2; // Store a portion of the xNext calculation that isn't affect by the loop.

                for( int y = 0; y < this->height(); ++y ) {
                    int xRhs = 0;
                    int yRhs = static_cast<int>(
                        std::floor( static_cast<float>( rhs.height() ) *
                                    ( ( static_cast<float>( y ) + 0.5f ) / static_cast<float>( this->height() ) ) ) );

                    for( int x = 0; x < this->width(); ++x ) {
                        // Calculating this way will handle non-even multiples and include the remainder appropriately.
                        // Rounds in a manner consistent with nearest neighbor. Some parts of the calculation can be
                        // pre-computed. xNext = static_cast<int>( ( static_cast<boost::int64_t>(x + 1) *
                        // static_cast<boost::int64_t>( rhs.width() ) + this.width() / 2 ) /
                        // static_cast<boost::int64_t>( this->width() ) );
                        int xNext = static_cast<int>(
                            ( static_cast<boost::int64_t>( x ) * static_cast<boost::int64_t>( rhs.width() ) + xRound ) /
                            static_cast<boost::int64_t>( this->width() ) );

                        assert( xNext > 0 && xNext <= rhs.width() );

                        float v = m_depthbuffer.get_pixel( x, y );
                        for( ; xRhs < xNext; ++xRhs ) {
                            float vRhs = rhs.m_depthbuffer.get_pixel( xRhs, yRhs );
                            if( vRhs < v )
                                v = vRhs;
                        }

                        m_depthbuffer.set_pixel( x, y, v );

                        assert( xRhs == xNext );
                    }
                }
            }
        } else if( rhs.width() < this->width() ) {
            // We are expanding horizontally, but shrinking (or maintaining the same size) vertically. Use
            // nearest-neighbor on x, but a min filter on a range of y.

            boost::int64_t yRound =
                rhs.height() +
                this->height() / 2; // Store a portion of the yNext calculation that isn't affect by the loop.

            int yRhs = 0;
            for( int y = 0; y < this->height(); ++y ) {
                // Calculating this way will handle non-even multiples and include the remainder appropriately. Rounds
                // in a manner consistent with nearest neighbor. Some parts of the calculation can be pre-computed.
                // yNext = static_cast<int>( ( static_cast<boost::int64_t>(y + 1) * static_cast<boost::int64_t>(
                // rhs.height() )
                // + this.width() / 2 ) / static_cast<boost::int64_t>( this->height() ) );
                int yNext = static_cast<int>(
                    ( static_cast<boost::int64_t>( y ) * static_cast<boost::int64_t>( rhs.height() ) + yRound ) /
                    static_cast<boost::int64_t>( this->height() ) );

                assert( yNext > 0 && yNext <= rhs.height() );

                int yCur = yRhs;
                for( int x = 0; x < this->width(); ++x ) {
                    int xRhs = static_cast<int>(
                        std::floor( static_cast<float>( rhs.width() ) *
                                    ( ( static_cast<float>( x ) + 0.5f ) / static_cast<float>( this->width() ) ) ) );

                    float v = m_depthbuffer.get_pixel( x, y );
                    for( yRhs = yCur; yRhs < yNext; ++yRhs ) {
                        float vRhs = rhs.m_depthbuffer.get_pixel( xRhs, yRhs );
                        if( vRhs < v )
                            v = vRhs;
                    }

                    m_depthbuffer.set_pixel( x, y, v );
                }

                assert( yRhs == yNext );
            }
        } else {
            // We are shrinking (or maintaining the same size) both vertically and horizontally.

            boost::int64_t xRound =
                rhs.width() +
                this->width() / 2; // Store a portion of the xNext calculation that isn't affect by the loop.
            boost::int64_t yRound =
                rhs.height() +
                this->height() / 2; // Store a portion of the yNext calculation that isn't affect by the loop.

            int yRhs = 0;
            for( int y = 0; y < this->height(); ++y ) {
                // Calculating this way will handle non-even multiples and include the remainder appropriately. Rounds
                // in a manner consistent with nearest neighbor. Some parts of the calculation can be pre-computed.
                // yNext = static_cast<int>( ( static_cast<boost::int64_t>(y + 1) * static_cast<boost::int64_t>(
                // rhs.height() )
                // + this.width() / 2 ) / static_cast<boost::int64_t>( this->height() ) );
                int yNext = static_cast<int>(
                    ( static_cast<boost::int64_t>( y ) * static_cast<boost::int64_t>( rhs.height() ) + yRound ) /
                    static_cast<boost::int64_t>( this->height() ) );

                assert( yNext > 0 && yNext <= rhs.height() );

                int xRhs = 0;
                int yCur = yRhs;
                for( int x = 0; x < this->width(); ++x ) {
                    // Calculating this way will handle non-even multiples and include the remainder appropriately.
                    // Rounds in a manner consistent with nearest neighbor. Some parts of the calculation can be
                    // pre-computed. xNext = static_cast<int>( ( static_cast<boost::int64_t>(x + 1) *
                    // static_cast<boost::int64_t>( rhs.width() )
                    // + this.width() / 2 ) / static_cast<boost::int64_t>( this->width() ) );
                    int xNext = static_cast<int>(
                        ( static_cast<boost::int64_t>( x ) * static_cast<boost::int64_t>( rhs.width() ) + xRound ) /
                        static_cast<boost::int64_t>( this->width() ) );
                    int xCur = xRhs;

                    assert( xNext > 0 && xNext <= rhs.width() );

                    float v = m_depthbuffer.get_pixel( x, y );
                    for( yRhs = yCur; yRhs < yNext; ++yRhs ) {
                        for( xRhs = xCur; xRhs < xNext; ++xRhs ) {
                            float vRhs = rhs.m_depthbuffer.get_pixel( xRhs, yRhs );
                            if( vRhs < v )
                                v = vRhs;
                        }
                    }

                    m_depthbuffer.set_pixel( x, y, v );

                    assert( xRhs == xNext );
                }

                assert( yRhs == yNext );
            }
        }
    }

    void render( const transform4f& xform,
                 const boost::shared_ptr<frantic::geometry::raytraced_geometry_collection>& geometry,
                 const frantic::graphics::camera<float>& cam, int superSampling = 1 ) {
        frantic::logging::console_progress_logger progress;
        render( xform, geometry, cam, superSampling, progress );
    }

    // Renders the given geometry into a depth buffer, from the point of view of the passed camera.
    void render( const transform4f& xform,
                 const boost::shared_ptr<frantic::geometry::raytraced_geometry_collection>& geometry,
                 const frantic::graphics::camera<float>& cam, int superSampling,
                 frantic::logging::progress_logger& progress ) {
        geometry_renderer gr;

        // Set the geometry of the renderer
        gr.set_geometry( geometry );

        gr.render_depth_image( m_depthbuffer, cam, superSampling, progress );
        m_transform = xform;
        m_transformInverse = m_transform.to_inverse();
    }

    // Returns true if the depth value at pixel(x,y) is > distance.
    bool is_visible( int x, int y, float distance ) const { return m_depthbuffer.get_pixel( x, y ) > distance; }

    // Does a visibility check of the world-space coordinate using the depth map
    bool is_visible( const vector3f& worldPosition, const frantic::graphics::camera<float>& cam ) {
        vector3f localPosition = m_transformInverse * worldPosition;
        return is_visible_cameraspace( localPosition, cam );
    }

    // Does a visibility check of the camera-space coordinate using the depth map
    // Currently this works using camera distance, NOT camera Z-Depth
    bool is_visible_cameraspace( const vector3f& cameraSpacePosition, const frantic::graphics::camera<float>& cam ) {
        bool outSuccess = true;

        frantic::graphics2d::vector2f lookupPixel = cam.from_cameraspace_position( cameraSpacePosition, outSuccess );
        if( outSuccess ) {
            float depth = m_depthbuffer.get_pixel_filtered( lookupPixel );
            return depth * depth > cameraSpacePosition.get_magnitude_squared();
        }

        return false;
    }

    bool is_visible_cameraspace_z( const vector3f& cameraSpacePosition, const frantic::graphics::camera<float>& cam ) {
        bool outSuccess = true;

        frantic::graphics2d::vector2f lookupPixel = cam.from_cameraspace_position( cameraSpacePosition, outSuccess );
        if( outSuccess ) {
            float depth = m_depthbuffer.get_pixel_filtered( lookupPixel );
            return depth > -cameraSpacePosition.z;
        }

        return false;
    }

    /**
     * This function will provide a visibility of a cameraspace position, using percentage closer filtering of the
     * nearest 4 pixels.
     *
     * @param cameraSpacePosition a location relative to cam
     * @param cam the camera that the depthbuffer_singleface was rendered from
     * @result a floating point value [0,1] that is the fraction of rays which fall upon
     */
    float get_pcf_visibility_cameraspace_z( const vector3f& cameraSpacePosition,
                                            const frantic::graphics::camera<float>& cam ) const {
        bool outSuccess = true;

        frantic::graphics2d::vector2f lookupPixel = cam.from_cameraspace_position( cameraSpacePosition, outSuccess );
        if( !outSuccess )
            return 1.f;

        float floorX = floor( lookupPixel.x - 0.5f );
        float floorY = floor( lookupPixel.y - 0.5f );
        float weightX = lookupPixel.x - floorX - 0.5f;
        float weightY = lookupPixel.y - floorY - 0.5f;

        int x = (int)floorX;
        int y = (int)floorY;

        // Start fully visible, removing visibility as we go. This makes areas out of bounds be fully visible.
        float vis = 1.f;
        if( (unsigned)y < (unsigned)cam.get_output_size().ysize ) {
            if( (unsigned)x < (unsigned)cam.get_output_size().xsize && !is_visible( x, y, -cameraSpacePosition.z ) )
                vis -= ( 1.f - weightX ) * ( 1.f - weightY );
            if( (unsigned)( x + 1 ) < (unsigned)cam.get_output_size().xsize &&
                !is_visible( x + 1, y, -cameraSpacePosition.z ) )
                vis -= weightX * ( 1.f - weightY );
        }
        if( (unsigned)( y + 1 ) < (unsigned)cam.get_output_size().ysize ) {
            if( (unsigned)x < (unsigned)cam.get_output_size().xsize && !is_visible( x, y + 1, -cameraSpacePosition.z ) )
                vis -= ( 1.f - weightX ) * weightY;
            if( (unsigned)( x + 1 ) < (unsigned)cam.get_output_size().xsize &&
                !is_visible( x + 1, y + 1, -cameraSpacePosition.z ) )
                vis -= weightX * weightY;
        }

        // Clamp to 0 in case FP arithmetic caused funkiness.
        if( vis < 0.f )
            return 0.f;
        return vis;
    }

    /**
     * This function will provide a visibility of a cameraspace position, using percentage closer filtering of the
     * nearest 4 pixels.
     *
     * @param lookupPixel a location in the depthbuffer to lookup
     * @param zDepth the z depth of the query point
     * @result a floating point value [0,1] that is the fraction of rays which fall upon
     */
    float get_pcf_visibility_z( const frantic::graphics2d::vector2f& lookupPixel, float zDepth ) const {
        float floorX = floor( lookupPixel.x - 0.5f );
        float floorY = floor( lookupPixel.y - 0.5f );
        float weightX = lookupPixel.x - floorX - 0.5f;
        float weightY = lookupPixel.y - floorY - 0.5f;

        int x = (int)floorX;
        int y = (int)floorY;

        float weightXY = weightX * weightY;

        // Start fully visible, removing visibility as we go. This makes areas out of bounds be fully visible.
        float vis = 1.f;
        if( (unsigned)y < (unsigned)height() ) {
            if( (unsigned)x < (unsigned)width() && !is_visible( x, y, zDepth ) )
                vis -= 1.f - weightX - weightY + weightXY;
            if( (unsigned)( x + 1 ) < (unsigned)width() && !is_visible( x + 1, y, zDepth ) )
                vis -= weightX - weightXY;
        }
        if( (unsigned)( y + 1 ) < (unsigned)height() ) {
            if( (unsigned)x < (unsigned)width() && !is_visible( x, y + 1, zDepth ) )
                vis -= weightY - weightXY;
            if( (unsigned)( x + 1 ) < (unsigned)width() && !is_visible( x + 1, y + 1, zDepth ) )
                vis -= weightXY;
        }

        // Clamp to 0 in case FP arithmetic caused funkiness.
        if( vis < 0.f )
            return 0.f;
        return vis;
    }

    float get_depth_value( int x, int y ) const { return m_depthbuffer.get_pixel( x, y ); }

    float get_depth_value_fast( int x, int y ) const { return m_depthbuffer.get_pixel_fast( x, y ); }

    void set_depth_value( int x, int y, float depth ) { m_depthbuffer.set_pixel( x, y, depth ); }

    void denormalize_values( float nearPlane, float farPlane, bool bInvert = false ) {
        float range = farPlane - nearPlane;
        float temp;

        frantic::graphics2d::framebufferiterator<float> it( &m_depthbuffer );
        while( it.isRowValid() ) {
            while( it.isColValid() ) {
                it.getData( &temp );

                if( bInvert )
                    it.setData( ( 1.f - temp ) * range + nearPlane );
                else
                    it.setData( temp * range + nearPlane );

                it.nextCol();
            }

            it.nextRowStart();
        }
    }

    void normalize_values( float nearPlane, float farPlane ) {
        float range = farPlane - nearPlane;
        float temp;

        frantic::graphics2d::framebufferiterator<float> it( &m_depthbuffer );
        while( it.isRowValid() ) {
            while( it.isColValid() ) {
                it.getData( &temp );

                if( temp < nearPlane )
                    temp = nearPlane;
                else if( temp > farPlane )
                    temp = farPlane;

                it.setData( ( temp - nearPlane ) / range );
                it.nextCol();
            }

            it.nextRowStart();
        }
    }

    void to_OpenEXR_file( const frantic::tstring& filename ) const { m_depthbuffer.to_OpenEXR_file( filename ); }

    void from_OpenEXR_file( const frantic::tstring& filename ) { m_depthbuffer.from_OpenEXR_file( filename ); }

    // Use this for file I/O only please
    const std::vector<float>& data() const { return m_depthbuffer.data(); }

    frantic::graphics2d::framebuffer<float>& as_framebuffer() { return m_depthbuffer; }
};

} // namespace rendering
} // namespace frantic
