// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#pragma warning( push, 3 )
#pragma warning( disable : 4512 )
#include <boost/random.hpp>
#pragma warning( pop )

#include <frantic/graphics2d/framebuffer.hpp>
#include <frantic/rendering/framebuffer_cubeface.hpp>
#include <frantic/rendering/lights/light_list.hpp>
#include <frantic/volumetrics/volumetric_scattering_functions.hpp>

#include <frantic/geometry/raytraced_geometry_collection.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>
#include <frantic/graphics/color_with_alpha.hpp>

#include <frantic/logging/console_progress_logger.hpp>

namespace frantic {
namespace rendering {
using frantic::geometry::raytrace_intersection;
using frantic::geometry::trimesh3;
using frantic::graphics::camera;
using frantic::graphics::color4f;
using frantic::graphics::ray3f;
using frantic::graphics::transform4f;
using frantic::graphics::vector3f;

class geometry_renderer {
    // This is the object which encapsulates all the geometry in the scene.
    boost::shared_ptr<frantic::geometry::raytraced_geometry_collection> m_geometry;
    // This is how many motion blur segments to use.  If this is less than 1, motion blur is disabled.
    int m_motionBlurSegments;
    // When motion blur is enabled, this indicates whether or not to use random jittering.
    bool m_jitteredMotionBlur;

  public:
    geometry_renderer() {
        m_motionBlurSegments = 0;
        m_jitteredMotionBlur = false;
    }

    void set_geometry( const boost::shared_ptr<frantic::geometry::raytraced_geometry_collection>& geometry ) {
        m_geometry = geometry;
    }

    void set_motionBlurSegments( int segments ) { m_motionBlurSegments = segments; }

    void set_jitteredMotionBlur( bool jitteredMotionBlur ) { m_jitteredMotionBlur = jitteredMotionBlur; }

    void render_z_depth_image( graphics2d::framebuffer<float>& renderFramebuffer,
                               const frantic::graphics::camera<float>& renderCamera, int superSampling = 1 ) {
        frantic::logging::console_progress_logger progress;
        render_z_depth_image( renderFramebuffer, renderCamera, superSampling, progress );
    }

    void render_z_depth_image( graphics2d::framebuffer<float>& renderFramebuffer,
                               const frantic::graphics::camera<float>& renderCamera, int superSampling,
                               frantic::logging::progress_logger& progress ) {
        using namespace frantic::graphics2d;

        if( renderCamera.output_size() != renderFramebuffer.size() )
            throw std::runtime_error( "geometry_renderer.render_z_depth_image: The camera output size, " +
                                      renderCamera.output_size().str() + ", is different from the framebuffer size, " +
                                      renderFramebuffer.size().str() );

        // Initialize the geometry object acceleration
        if( m_geometry )
            m_geometry->prepare_kdtrees( 0.5f );

        raytrace_intersection isect;

        vector2 p;
        size2 sz = renderFramebuffer.size();
        size2 szSuper = superSampling * sz;
        bool isValid;
        for( p.y = 0; p.y < sz.ysize; ++p.y ) {
            //			cerr.flush();
            for( p.x = 0; p.x < sz.xsize; ++p.x ) {
                float depth = ( std::numeric_limits<float>::max )();
                for( int sy = 0; sy < superSampling; ++sy ) {
                    for( int sx = 0; sx < superSampling; ++sx ) {
                        vector2f fPixel( ( superSampling * p.x + sx + 0.5f ) / superSampling,
                                         ( superSampling * p.y + sy + 0.5f ) / superSampling );
                        isValid = true;
                        ray3f ray = renderCamera.get_worldspace_ray( fPixel, isValid );

                        if( isValid ) {
                            bool rayIntersected =
                                m_geometry->intersect_ray( ray, 0, ( std::numeric_limits<float>::max )(), isect );

                            if( rayIntersected ) {
                                // Apply the min filter
                                float testDepth = -( renderCamera.world_transform_inverse() * isect.position ).z;
                                if( testDepth < depth )
                                    depth = testDepth;
                            }
                        }
                    }
                }
                renderFramebuffer.set_pixel( p, depth );
            }
            progress.update_progress( p.y + 1, sz.ysize );
        }
        if( frantic::logging::is_logging_progress() )
            std::cout << std::endl;
    }

    void render_depth_image( graphics2d::framebuffer<float>& renderFramebuffer,
                             const frantic::graphics::camera<float>& renderCamera, int superSampling = 1 ) {
        frantic::logging::console_progress_logger progress;
        render_depth_image( renderFramebuffer, renderCamera, superSampling, progress );
    }

    void render_depth_image( graphics2d::framebuffer<float>& renderFramebuffer,
                             const frantic::graphics::camera<float>& renderCamera, int superSampling,
                             frantic::logging::progress_logger& progress ) {
        using namespace frantic::graphics2d;

        if( renderCamera.output_size() != renderFramebuffer.size() )
            throw std::runtime_error( "geometry_renderer.render_depth_image: The camera output size, " +
                                      renderCamera.output_size().str() + ", is different from the framebuffer size, " +
                                      renderFramebuffer.size().str() );

        // Initialize the geometry object acceleration
        m_geometry->prepare_kdtrees( 0.5f );

        raytrace_intersection isect;

        vector2 p;
        size2 sz = renderFramebuffer.size();
        size2 szSuper = superSampling * sz;
        bool isValid;
        for( p.y = 0; p.y < sz.ysize; ++p.y ) {
            //			cerr.flush();
            for( p.x = 0; p.x < sz.xsize; ++p.x ) {
                float depth = ( std::numeric_limits<float>::max )();
                for( int sy = 0; sy < superSampling; ++sy ) {
                    for( int sx = 0; sx < superSampling; ++sx ) {
                        vector2f fPixel( ( superSampling * p.x + sx + 0.5f ) / superSampling,
                                         ( superSampling * p.y + sy + 0.5f ) / superSampling );
                        isValid = true;
                        ray3f ray = renderCamera.get_worldspace_ray( fPixel, isValid );

                        if( isValid ) {
                            bool rayIntersected =
                                m_geometry->intersect_ray( ray, 0, ( std::numeric_limits<float>::max )(), isect );

                            if( rayIntersected ) {
                                // Apply the min filter
                                float testDepth =
                                    ( renderCamera.world_transform_inverse() * isect.position ).get_magnitude();
                                if( testDepth < depth )
                                    depth = testDepth;
                            }
                        }
                    }
                }
                renderFramebuffer.set_pixel( p, depth );
            }
            progress.update_progress( p.y + 1, sz.ysize );
        }
        if( frantic::logging::is_logging_progress() )
            std::cout << std::endl;
    }

    void render( graphics2d::framebuffer<frantic::graphics::color4f>& renderFramebuffer,
                 const frantic::graphics::camera<float>& renderCamera, int superSampling = 1 ) {
        using namespace frantic::graphics;

        if( m_motionBlurSegments <= 0 ) {
            render_subinterval( renderFramebuffer, renderCamera, superSampling, 0.5f, 0.5f, 42 );
        } else if( m_motionBlurSegments == 1 ) {
            render_subinterval( renderFramebuffer, renderCamera, superSampling, 0.f, 1.f, 42 );
        } else {
            renderFramebuffer.fill( color4f() );
            graphics2d::framebuffer<frantic::graphics::color4f> tempBuffer( renderFramebuffer.size() );
            for( int segment = 0; segment < m_motionBlurSegments; ++segment ) {
                float startTime = float( segment ) / m_motionBlurSegments,
                      endTime = float( segment + 1 ) / m_motionBlurSegments;
                std::cout << "Processing motion blur pass " << ( segment + 1 ) << " of " << m_motionBlurSegments
                          << " (time " << startTime << " to " << endTime << ")" << std::endl;
                render_subinterval( tempBuffer, renderCamera, superSampling, startTime, endTime, 67 * segment + 52,
                                    startTime * 100, endTime * 100 );
                tempBuffer.apply_gain( 1.f / m_motionBlurSegments );
                renderFramebuffer.add_image_data( tempBuffer );
            }
        }
    }

    void render( framebuffer_cubeface<frantic::graphics::color4f>& renderFramebuffer,
                 const frantic::graphics::camera<float>& renderCamera ) {
        for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
            frantic::graphics::camera<float> cubeFaceCamera = frantic::graphics::camera<float>::from_cube_face(
                renderCamera.world_transform(), (frantic::graphics::cube_face::default_cube_face)cubeFace );
            cubeFaceCamera.set_output_size(
                frantic::graphics2d::size2( renderFramebuffer.size(), renderFramebuffer.size() ) );
            render( renderFramebuffer.get_cubeface_framebuffer( cubeFace ), cubeFaceCamera, 1 );
        }
    }

    // Renders the distance from the camera
    void render_depth_image( rendering::framebuffer_cubeface<float>& renderFramebuffer, const transform4f& xform,
                             int superSampling = 1 ) {
        frantic::logging::console_progress_logger progress;
        render_depth_image( renderFramebuffer, xform, superSampling, progress );
    }

    // Renders the distance from the camera
    void render_depth_image( rendering::framebuffer_cubeface<float>& renderFramebuffer, const transform4f& xform,
                             int superSampling, frantic::logging::progress_logger& progress ) {
        if( frantic::logging::is_logging_stats() )
            std::cerr << "Rendering depth image at resolution " << renderFramebuffer.size() << std::endl;

        for( int cubeFace = 0; cubeFace < 6; ++cubeFace ) {
            frantic::graphics::camera<float> cubeFaceCamera = frantic::graphics::camera<float>::from_cube_face(
                xform, (frantic::graphics::cube_face::default_cube_face)cubeFace );
            cubeFaceCamera.set_output_size(
                frantic::graphics2d::size2( renderFramebuffer.size(), renderFramebuffer.size() ) );
            frantic::logging::progress_logger_subinterval_tracker plst( progress, 100 * cubeFace / 6.f,
                                                                        100 * ( cubeFace + 1 ) / 6.f );
            render_depth_image( renderFramebuffer.get_cubeface_framebuffer( cubeFace ), cubeFaceCamera, superSampling,
                                progress );
        }

        /*
        // Initialize the geometry object acceleration
        m_geometry->prepare_kdtrees( 0.5f );

        raytrace_intersection isect;

        vector3f pos = renderCamera.camera_position();
        for( int sample = 0; sample < 500000; ++sample ) {
          vector3f dir = vector3f::from_unit_random();

          bool rayIntersected = m_geometry->intersect_ray( ray3f(pos,dir), 0, (numeric_limits<float>::max)(), isect );

          if( rayIntersected ) {
            vector3f isectPoint = renderCamera.world_transform_inverse() * isect.position;
            renderFramebuffer.draw_point( renderCamera.world_transform_inverse().transform_no_translation(dir), color4f(
        color3f(isectPoint.get_magnitude()), 1) );
          }
        }
        */
    }

    void render_subinterval( graphics2d::framebuffer<frantic::graphics::color4f>& renderFramebuffer,
                             const frantic::graphics::camera<float>& renderCamera, int superSampling, float startTime,
                             float endTime, boost::uint32_t seed, float progressStart = 0, float progressEnd = 100 ) {
        using namespace std;
        using namespace frantic::graphics2d;

        if( renderCamera.output_size() != renderFramebuffer.size() )
            throw std::runtime_error( "geometry_renderer.render_subinterval: The camera output size, " +
                                      renderCamera.output_size().str() + ", is different from the framebuffer size, " +
                                      renderFramebuffer.size().str() );

        float centerTime = 0.5f * ( startTime + endTime );

        // Create a random number generator
        boost::mt19937 generator( seed );
        boost::uniform_real<float> uni_dist( 0, 1 );
        boost::variate_generator<boost::mt19937&, boost::uniform_real<float>> rng( generator, uni_dist );

        // Initialize the geometry object acceleration at the interval center
        m_geometry->prepare_kdtrees( centerTime );

        // Determine whether or not to use jittering
        bool jittered = m_jitteredMotionBlur;
        if( startTime >= endTime )
            jittered = false;

        raytrace_intersection isect;

        vector2 pixel;
        size2 sz = renderFramebuffer.size();
        size2 szSuper = superSampling * sz;
        bool isValid;
        for( pixel.y = 0; pixel.y < sz.ysize; ++pixel.y ) {
            cerr.flush();
            for( pixel.x = 0; pixel.x < sz.xsize; ++pixel.x ) {
                /*
                std::cout << "\nRendering pixel " << pixel << "\n";
                //*/
                color4f pixelColor;
                std::vector<raytrace_intersection> intersections;
                for( int sy = 0; sy < superSampling; ++sy ) {
                    for( int sx = 0; sx < superSampling; ++sx ) {
                        vector2f fPixel( ( superSampling * pixel.x + sx + 0.5f ) / superSampling,
                                         ( superSampling * pixel.y + sy + 0.5f ) / superSampling );

                        float sampleTime = centerTime;
                        if( jittered ) {
                            sampleTime = rng() * ( endTime - startTime ) + startTime;
                        }

                        isValid = true;
                        ray3f ray = renderCamera.get_worldspace_ray( fPixel, sampleTime, isValid );

                        if( isValid ) {
                            double probeTime = diagnostics::get_milliseconds();
                            bool rayIntersected = false;
                            // if( true ) {
                            rayIntersected =
                                m_geometry->intersect_ray( ray, 0, ( numeric_limits<float>::max )(), isect );
                            //} else {
                            //	intersections.clear();
                            //	m_geometry->intersect_ray_all( ray, (numeric_limits<float>::max)(), sampleTime,
                            // intersections ); 	rayIntersected = intersections.size() > 0; 	if( rayIntersected )
                            // isect = intersections[0];
                            //}
                            probeTime = diagnostics::get_milliseconds() - probeTime;

                            float bChannel;
                            // bChannel = (float)probeTime;
                            bChannel = (float)isect.distance;

                            if( rayIntersected ) {
                                float IdotN = vector3f::dot( ray.direction(), isect.geometricNormal );
                                if( IdotN < 0 )
                                    pixelColor += color4f( -IdotN, -IdotN, bChannel, 1 );
                                else {
                                    pixelColor += color4f( IdotN, 0, bChannel, 1 );
                                    // pixelColor += color4f( IdotN, IdotN, IdotN, 1 );
                                }
                                // std::cout << pixelColor.c << std::endl;
                            }
                        }
                    }
                }
                pixelColor *= 1.f / ( superSampling * superSampling );
                renderFramebuffer.set_pixel( pixel, pixelColor );
            }
            cout << "\rProgress: " << pixel.y + 1 << "/" << sz.ysize << " ("
                 << ( ( (float)pixel.y / ( sz.ysize - 1 ) ) * ( progressEnd - progressStart ) + progressStart )
                 << "%)                     ";
            cout.flush();
        }
        cout << endl;
    }
};

} // namespace rendering
} // namespace frantic
