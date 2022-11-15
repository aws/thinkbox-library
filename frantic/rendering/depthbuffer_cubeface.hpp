// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/rendering/framebuffer_cubeface.hpp>
#include <frantic/rendering/geometry_renderer.hpp>

#include <frantic/logging/console_progress_logger.hpp>

namespace frantic {
namespace rendering {

class depthbuffer_cubeface {
    transform4f m_transform, m_transformInverse;
    framebuffer_cubeface<float> m_depthbuffer;

  public:
    depthbuffer_cubeface() { m_depthbuffer.clear(); }

    depthbuffer_cubeface( int mapSize )
        : m_depthbuffer( mapSize ) {}

    // Resizes the depth buffer
    void set_size( int mapSize ) { m_depthbuffer.set_size( mapSize ); }

    void set( const framebuffer_cubeface<float>& depthBuffer ) { m_depthbuffer = depthBuffer; }

    void set_with_swap( framebuffer_cubeface<float>& depthBuffer ) { m_depthbuffer.swap( depthBuffer ); }

    void set_transform( const transform4f& xform ) {
        m_transform = xform;
        m_transformInverse = m_transform.to_inverse();
    }

    void render( const transform4f& xform,
                 const boost::shared_ptr<frantic::geometry::raytraced_geometry_collection>& geometry,
                 int superSampling = 1 ) {
        frantic::logging::console_progress_logger progress;
        render( xform, geometry, superSampling, progress );
    }

    // Renders the given geometry into a depth buffer, from the
    void render( const transform4f& xform,
                 const boost::shared_ptr<frantic::geometry::raytraced_geometry_collection>& geometry, int superSampling,
                 frantic::logging::progress_logger& progress ) {
        geometry_renderer gr;

        // Set the geometry of the renderer
        gr.set_geometry( geometry );

        gr.render_depth_image( m_depthbuffer, xform, superSampling, progress );
        m_transform = xform;
        m_transformInverse = m_transform.to_inverse();
    }

    // Does a visibility check of the world-space coordinate using the depth map
    bool is_visible( const vector3f& worldPosition ) {
        vector3f localPosition = m_transformInverse * worldPosition;

        float depth = m_depthbuffer.get_pixel_bilinear( localPosition );

        return depth * depth > localPosition.get_magnitude_squared();
    }

    bool is_visible_zdepth( const vector3f& worldPosition ) {
        float zdepth;
        vector3f localPosition = m_transformInverse * worldPosition;

        frantic::graphics::cube_face::default_cube_face cubeFace = get_cube_face( localPosition );
        frantic::graphics2d::vector2f coord = get_cube_face_coordinate_and_zdepth( localPosition, cubeFace, zdepth );

        float weight = 0.f;
        float imgDepth =
            m_depthbuffer.get_cubeface_framebuffer( cubeFace ).get_pixel_bilinear_return_weight( coord, weight );

        // TODO: This isn't the best treatment of this issue, however I'm not going to dwell on it particularly long.
        //        It may be better to do something similar to framebuffer_cubeface<>::get_pixel_bilinear().
        return zdepth < ( imgDepth / weight );
    }

    float get_pcf_visibility_cameraspace_z( const frantic::graphics::vector3f& camSpacePos ) const {
        frantic::graphics::cube_face::default_cube_face cubeFaces[6];
        int numFaces = frantic::graphics::get_cube_faces(
            camSpacePos, ( (float)m_depthbuffer.size() - 2.f ) / (float)m_depthbuffer.size(), cubeFaces );

        float totalBlocked = 0.f;
        float totalWeight = 0.f;

        float z;
        frantic::graphics2d::vector2f coord;
        for( int i = 0; i < numFaces; ++i ) {
            const frantic::graphics2d::framebuffer<float>& cubefaceBuffer =
                m_depthbuffer.get_cubeface_framebuffer( cubeFaces[i] );
            coord = frantic::graphics::get_cube_face_coordinate_and_zdepth( camSpacePos, cubeFaces[i], z );
            coord += frantic::graphics2d::vector2f( 1, 1 );
            coord *= ( 0.5f * (float)m_depthbuffer.size() );

            float floorX = floor( coord.x - 0.5f );
            float alphaX = coord.x - floorX - 0.5f;
            int pixelX = (int)floorX;

            float floorY = floor( coord.y - 0.5f );
            float alphaY = coord.y - floorY - 0.5f;
            int pixelY = (int)floorY;

            if( (unsigned)pixelY < (unsigned)m_depthbuffer.size() ) {
                if( (unsigned)pixelX < (unsigned)m_depthbuffer.size() ) {
                    float imgZ = cubefaceBuffer.get_pixel( pixelX, pixelY );
                    float weight = ( 1.f - alphaX ) * ( 1.f - alphaY );
                    if( imgZ < z )
                        totalBlocked += weight;
                    totalWeight += weight;
                }
                if( (unsigned)( pixelX + 1 ) < (unsigned)m_depthbuffer.size() ) {
                    float imgZ = cubefaceBuffer.get_pixel( pixelX + 1, pixelY );
                    float weight = alphaX * ( 1.f - alphaY );
                    if( imgZ < z )
                        totalBlocked += weight;
                    totalWeight += weight;
                }
            }
            if( (unsigned)( pixelY + 1 ) < (unsigned)m_depthbuffer.size() ) {
                if( (unsigned)pixelX < (unsigned)m_depthbuffer.size() ) {
                    float imgZ = cubefaceBuffer.get_pixel( pixelX, pixelY + 1 );
                    float weight = ( 1.f - alphaX ) * alphaY;
                    if( imgZ < z )
                        totalBlocked += weight;
                    totalWeight += weight;
                }
                if( (unsigned)( pixelX + 1 ) < (unsigned)m_depthbuffer.size() ) {
                    float imgZ = cubefaceBuffer.get_pixel( pixelX + 1, pixelY + 1 );
                    float weight = alphaX * alphaY;
                    if( imgZ < z )
                        totalBlocked += weight;
                    totalWeight += weight;
                }
            }
        }

        if( totalWeight == 0.f )
            return 1.f;
        return 1.f - ( totalBlocked / totalWeight );
    }

    void to_longlat_OpenEXR_file_zup( const frantic::tstring& filename, int width = -1, int height = -1,
                                      int superSampling = 2 ) {
        m_depthbuffer.to_longlat_OpenEXR_file_zup( filename, width, height, superSampling );
    }

    void to_longlat_OpenEXR_file_yup( const frantic::tstring& filename, int width = -1, int height = -1,
                                      int superSampling = 2 ) {
        m_depthbuffer.to_longlat_OpenEXR_file_yup( filename, width, height, superSampling );
    }

    frantic::rendering::framebuffer_cubeface<float>& as_framebuffer() { return m_depthbuffer; }
};

} // namespace rendering
} // namespace frantic
