// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/motion_blurred_transform.hpp>
#include <frantic/graphics/plane3f.hpp>
#include <frantic/graphics/quat4f.hpp>
#include <frantic/graphics/ray3f.hpp>
#include <frantic/logging/logging_level.hpp>

#include <boost/io/ios_state.hpp>
#include <iomanip>

// Only do these 3ds Max things if 3ds Max is already included.  That way, the camera class can integrate nicely in 3ds
// Max, but still be used outside of 3ds Max.
#if defined( MAX_VERSION )
#include <frantic/max3d/convert.hpp>
#include <frantic/max3d/geopipe/xref_utils.hpp>
#include <frantic/max3d/parameter_extraction.hpp>
#include <frantic/max3d/units.hpp>
#endif

namespace frantic {
namespace graphics {

namespace projection_mode {
enum projection_mode {
    perspective,
    orthographic,
    spherical, // aka fisheye
    panoramic,
    vray_quadratic, // attempts to match the vray quadratic spherical distortion.
    vray_cubic      // attempts to match the vray cubic spherical distortion.
};
}

template <class T>
inline void from_panoramic_yup( const T* v2, const T* size, T hfov, T pixelAspect, T* v3Out );

template <class T>
inline void to_panoramic_yup( const T* v3, const T* size, T hfov, T pixelAspect, T* v2Out );

template <typename FloatType>
inline frantic::graphics::vector3t<FloatType> from_panoramic_yup( const frantic::graphics2d::vector2f& v,
                                                                  const frantic::graphics2d::size2f& size,
                                                                  FloatType hfov, FloatType pixelAspect = 1.0f ) {
    frantic::graphics::vector3t<FloatType> result;
    FloatType v2[2] = { v.x, v.y };
    FloatType s2[2] = { size.xsize, size.ysize };
    from_panoramic_yup<FloatType>( v2, s2, hfov, pixelAspect, &result[0] );
    return result;
}

template <typename FloatType>
inline frantic::graphics::vector3t<FloatType> from_panoramic_yup( const frantic::graphics2d::vector2& v,
                                                                  const frantic::graphics2d::size2& size,
                                                                  FloatType hfov, FloatType pixelAspect = 1.0f ) {
    return from_panoramic_yup( frantic::graphics2d::vector2f( v.x + 0.5f, v.y + 0.5f ),
                               frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ), hfov,
                               pixelAspect );
}

template <typename FloatType>
inline frantic::graphics2d::vector2f to_panoramic_yup( const frantic::graphics::vector3t<FloatType>& v,
                                                       const frantic::graphics2d::size2f& size, FloatType hfov,
                                                       FloatType pixelAspect = 1.0f ) {
    FloatType v2[2] = { 0, 0 };
    FloatType s2[2] = { size.xsize, size.ysize };
    to_panoramic_yup( &v[0], s2, hfov, pixelAspect, v2 );
    return frantic::graphics2d::vector2f( (float)v2[0], (float)v2[1] );
}

template <typename FloatType>
inline frantic::graphics2d::vector2 to_panoramic_yup( const frantic::graphics::vector3t<FloatType>& v,
                                                      const frantic::graphics2d::size2& size, FloatType hfov,
                                                      FloatType pixelAspect = 1.0f ) {
    frantic::graphics2d::vector2f result = to_panoramic_yup(
        v, frantic::graphics2d::size2f( float( size.xsize ), float( size.ysize ) ), hfov, pixelAspect );
    return frantic::graphics2d::vector2( (int)floorf( result.x ), (int)floorf( result.y ) );
}

/**
 * Notes about transform spaces:
 *   Unless specified in the method comments, there are three spaces that we are working with:
 *     [2D image space] <-> [2D camera space] <-> [3D world space]
 *   Where:
 *     Image Space is given in pixels with the x domain given in [0, m_outputSize.xsize) and y domain given in [0,
 * m_outputSize.ysize).  m_subpixelOffset is applied automatically for each pixel when going to and from camera space.
 *     Camera Space x, y, and optionally z positions are given in [-1, 1] and is pixel aspect corrected.
 *     World Space is, well, world space.  Internally we store the transform going from camera to world space
 */
template <class FloatType>
class camera {
    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
  protected:
    frantic::graphics2d::size2 m_outputSize;

    // The type of projection the camera is using.  Currently only perspective is implemented
    projection_mode::projection_mode m_projectionMode;

    // The field of view, and the tangent of 1/2 the field of view.
    // The direction depends on the image size ratio after applying the pixel aspect
    // If the image is wider than it is tall, it specifies the vertical field of view.  Otherwise, it specifies the
    // horizontal field of view
    FloatType m_fov, m_tanHalfFov;

    // For an orthographic camera, this is the width/height of the view in camera space units.
    // The direction depends on the image size ratio after applying the pixel aspect.
    // If the image is wider than it is tall, it specifies the height.  Otherwise, it specifies the width
    FloatType m_orthographicSize;

    // pixel aspect ratio given as x over y
    FloatType m_pixelAspect;

    // Near and far clipping planes
    FloatType m_near, m_far;

    // A subpixel offset applied after the camera transform
    frantic::graphics2d::vector2f m_subpixelOffset;

    // Transform matrix, includes camera animation over the motion blur segment
    // This transforms a camera space point into a world space point.
    motion_blurred_transform<FloatType> m_cameraToWorldTransform;

    /////////////
    // Depth of field information
    //  fStop:            Number of f/stops of the lens setup.
    //  focalLength:      The focal length of the lens setup, in camera-space units.
    //  focalDistance:    The distance from the camera (lens) that appears sharp, in camera-space units.
    //                    (Often equivalent to the target distance in a 3ds Max camera.)
    FloatType m_fStop, m_focalLength, m_focalDistance;
    //  filmDistance:     The distance from the film to the lens.
    //  apertureDiameter: The diameter of the lens letting in light.  For CG rendering we generally just
    //                    have this affect the depth of field blurring, not the light gathering power.
    //
    // Some depth of field identities:
    // 1/focalLength = 1/focalDistance + 1/filmDistance
    // fStop = focalLength/apertureDiameter
    ////////////

    // Kappa1 from a cubic radial distortion.
    FloatType m_distortionKappa;

    //---------------------------------------------------------------------------------------------------------

  public:
    ///// Constructors /////

    camera() {
        using namespace std;

        m_projectionMode = projection_mode::perspective;

        m_near = 0.001f;
        m_far = 1e+10;

        m_outputSize = frantic::graphics2d::size2( 640, 480 );

        // Default to no motion blur and no depth of field
        m_fStop = 1e30f;
        m_focalLength = 30;
        m_focalDistance = 100;

        // Default to pixel aspect of 1.0
        m_pixelAspect = 1.0f;
        m_subpixelOffset = frantic::graphics2d::vector2f( 0, 0 );

        set_horizontal_fov_degrees( 45.f );
        set_orthographic_width( 400 );
    }

    camera( const transform4t<FloatType>& transform, FloatType horizontalFovDegrees )
        : m_cameraToWorldTransform( transform ) {
        using namespace std;

        m_cameraToWorldTransform.remove_scale();

        m_projectionMode = projection_mode::perspective;

        m_near = 0.001f;
        m_far = 1e+10;

        m_outputSize = frantic::graphics2d::size2( 640, 480 );

        // Default to no depth of field
        m_fStop = 1e30f;
        m_focalLength = 30;
        m_focalDistance = 100;

        // Default to pixel aspect of 1.0
        m_pixelAspect = 1.0f;

        set_horizontal_fov_degrees( horizontalFovDegrees );
        set_orthographic_width( 400 );
    }

    // The point of this constructor is to allow the camera to be constructed for
    // any projection without any DOF by default.
    camera( projection_mode::projection_mode projectionMode, const motion_blurred_transform<FloatType>& transform,
            FloatType FOV_or_Width, frantic::graphics2d::size2 outputSize = frantic::graphics2d::size2( 640, 480 ),
            FloatType nearClip = 0.001f, FloatType farClip = 1e+10, FloatType pixelAspect = 1.0f )
        : m_projectionMode( projectionMode )
        , m_cameraToWorldTransform( transform )
        , m_outputSize( outputSize )
        , m_near( nearClip )
        , m_far( farClip )

        , m_fStop( 1e30f ) // Default to no depth of field
        , m_focalLength( 30 )
        , m_focalDistance( 100 )

        , m_pixelAspect( pixelAspect ) {
        using namespace std;

        set_horizontal_fov( FOV_or_Width );
        set_orthographic_width( FOV_or_Width );

        m_cameraToWorldTransform.remove_scale();
    }

    camera( projection_mode::projection_mode projectionMode, FloatType horizontalFovDegrees, FloatType nearDist,
            FloatType farDist, FloatType fStop, FloatType focalLength, FloatType focalDistance,
            transform4t<FloatType> worldTransform,
            const std::vector<transform4t<FloatType>>& animatedWorldTransform = std::vector<transform4t<FloatType>>(),
            frantic::graphics2d::size2 outputSize = frantic::graphics2d::size2( 640, 480 ),
            FloatType pixelAspect = 1.0f, FloatType motionBlurInterval = 0.5f )
        : m_cameraToWorldTransform( worldTransform, animatedWorldTransform, motionBlurInterval ) {
        using namespace std;

        m_cameraToWorldTransform.remove_scale();

        m_projectionMode = projectionMode;

        m_near = nearDist;
        m_far = farDist;

        m_outputSize = outputSize;

        m_fStop = fStop;
        m_focalLength = focalLength;
        m_focalDistance = focalDistance;

        m_pixelAspect = pixelAspect;

        set_horizontal_fov_degrees( horizontalFovDegrees );
        set_orthographic_width( 400 );
    }

// If the 3ds Max header object.h was included, add a constructor from a 3ds Max camera state.
// According to
// https://knowledge.autodesk.com/support/3ds-max/learn-explore/caas/CloudHelp/cloudhelp/2017/ENU/3DSMax/files/GUID-BFED400B-685A-47D1-A52C-33213FFAE356-htm.html,
// 3ds Max field of view is the horizontal field of view
#ifdef MAX_VERSION
    camera( const transform4f& cameraTransform, const CameraState& cameraState )
        : m_cameraToWorldTransform( cameraTransform ) {
        m_cameraToWorldTransform.remove_scale();

        if( cameraState.isOrtho )
            m_projectionMode = projection_mode::orthographic;
        else
            m_projectionMode = projection_mode::perspective;

        m_near = cameraState.hither;
        m_far = cameraState.yon;

        m_outputSize =
            frantic::graphics2d::size2( GetCOREInterface()->GetRendWidth(), GetCOREInterface()->GetRendHeight() );

        // Default to no depth of field
        m_fStop = 1e30f;
        m_focalLength = 30;
        m_focalDistance = cameraState.tdist;

        m_pixelAspect = GetCOREInterface()->GetRendApect();

        set_horizontal_fov( cameraState.fov );
        // set_orthographic_width( 400 * viewParams.zoom );
        set_orthographic_width( 2 * ( cameraState.tdist * horizontal_tan_half_fov() ) );
    }

    camera( const ViewParams& viewParams, frantic::graphics2d::size2 outputSize ) {
        // The affine transform in the viewparams sends world coordinates to object coordinates,
        // so we need to invert it.
        m_cameraToWorldTransform.set( viewParams.affineTM );
        m_cameraToWorldTransform.invert();
        m_cameraToWorldTransform.remove_scale_not_skew();

        // TODO: Use the prevAffineTM transform to derive camera motion blur.

        switch( viewParams.projType ) {
        case PROJ_PERSPECTIVE:
            m_projectionMode = projection_mode::perspective;
            break;
        case PROJ_PARALLEL:
            m_projectionMode = projection_mode::orthographic;
            break;
        default:
            throw std::runtime_error(
                "camera::camera: Cannot construct camera from a ViewParams, it has an invalid projection type " +
                boost::lexical_cast<std::string>( viewParams.projType ) );
        }

        m_near = viewParams.hither;
        m_far = viewParams.yon;

        m_outputSize = outputSize;
        // m_outputSize = size2( GetCOREInterface()->GetRendWidth(), GetCOREInterface()->GetRendHeight() );

        // Default to no depth of field
        m_fStop = 1e30f;
        m_focalLength = 30;
        m_focalDistance = 100;

        m_pixelAspect = GetCOREInterface()->GetRendApect();

        set_horizontal_fov( viewParams.fov );
        set_orthographic_width( 400 * viewParams.zoom );
    }

  private:
    // A helper function that will extract camera information from a light source. The cameras built from
    // direct lights have significant differences compared to those shown in the max viewports but it seems like
    // ours are prettier.
    void set_camera_from_light( LightObject* obj, TimeValue t, frantic::graphics2d::size2 outputSize ) {
        Class_ID objClassId = obj->ClassID();

        RefResult r;
        LightState ls;
        Interval ivl = FOREVER;
        r = obj->EvalLightState( t, ivl, &ls );
        if( r == REF_FAIL )
            throw std::runtime_error( "Failed to evaluate light in order to extract camera information." );

        if( ls.type == SPOT_LGT || ls.type == DIRECT_LGT ) {
            float fallsizeRads = frantic::math::degrees_to_radians( ls.fallsize );
            float sqrtAspect = sqrt( ls.aspect ); // Used to keep visible area constant over aspect changes.

            m_projectionMode = ( ls.type == SPOT_LGT ) ? projection_mode::perspective : projection_mode::orthographic;

            m_near = 0.001f;
            m_far = 1e+10;

            m_outputSize = outputSize;

            m_pixelAspect = ls.aspect * m_outputSize.ysize / m_outputSize.xsize;

            // Default to no depth of field
            m_fStop = 1e30f;
            m_focalLength = 30;
            m_focalDistance = 100;

            set_horizontal_fov( ( ls.shape == RECT_LIGHT ) ? 2.0f * ::atanf( ::tanf( fallsizeRads / 2 ) * sqrtAspect )
                                                           : fallsizeRads );
            set_orthographic_width( ls.fallsize );
        } else
            throw std::runtime_error( "Cannot extract camera info from lights that aren't spot or direct types." );

        return;
    }

  public:
    camera( INode* cameraNode, TimeValue t, frantic::graphics2d::size2 outputSize, float motionBlurInterval = 0.5f,
            float shutterBias = 0.f, Interval* ivalid = 0, bool disableCameraMotionBlur = false ) {
        m_cameraToWorldTransform = frantic::graphics::motion_blurred_transform<FloatType>::from_objtmafterwsm(
            cameraNode, t, motionBlurInterval, shutterBias, ivalid, disableCameraMotionBlur );
        m_cameraToWorldTransform.remove_scale_not_skew();

        FF_LOG( debug ) << "Processing camera " << cameraNode->GetName() << std::endl;

        Object* cameraObjectPtr = frantic::max3d::geopipe::get_object_ref( cameraNode )->FindBaseObject();
        if( cameraObjectPtr->SuperClassID() == LIGHT_CLASS_ID ) {
            set_camera_from_light( static_cast<LightObject*>( cameraObjectPtr ), t, outputSize );
            return;
        }

        if( cameraObjectPtr->SuperClassID() != CAMERA_CLASS_ID ) {
            throw std::runtime_error( std::string() + "camera::camera: Cannot construct camera from INode \"" +
                                      frantic::strings::to_string( cameraNode->GetName() ) +
                                      "\", because it does not have superclass of CAMERA_CLASS_ID" );
        }

        CameraObject* cameraObject = static_cast<CameraObject*>( cameraObjectPtr );

        CameraState cameraState;
        Interval cameraValid = FOREVER;
        if( cameraObject->EvalCameraState( t, cameraValid, &cameraState ) != REF_SUCCEED ) {
            throw std::runtime_error( std::string() + "camera::camera: Couldn't evaluate the camera state of node \"" +
                                      frantic::strings::to_string( cameraNode->GetName() ) + "\"." );
        }

        // Intersect the camera validity interval
        if( ivalid != 0 )
            *ivalid &= cameraValid;

        m_projectionMode = ( cameraState.isOrtho ) ? projection_mode::orthographic : projection_mode::perspective;

        if( cameraState.manualClip ) {
            m_near = cameraState.hither;
            m_far = cameraState.yon;
        } else {
            m_near = 0.001f;
            m_far = 1e+10f;
        }

        m_outputSize = outputSize;

        m_pixelAspect = GetCOREInterface()->GetRendApect();

        set_horizontal_fov( cameraState.fov );

        FloatType tanHalfFov = horizontal_tan_half_fov();

        set_orthographic_width( 2 * ( cameraState.tdist * tanHalfFov ) );

        // Default to no depth of field
        m_fStop = 1e30f;

        bool useTargetDistance = true;
        m_focalLength = GetCOREInterface()->GetRendApertureWidth() / ( 2 * tanHalfFov );

        // The focal length is provided in millimeters, so we need to convert it to camera-space units.
        m_focalLength /= (float)frantic::max3d::get_scale_to_millimeters();

        // Brazil Camera has DOF settings we can grab
        if( cameraObject->ClassID() == Class_ID( 966714866, 1481273768 ) ) {
            FF_LOG( debug ) << "It's a Brazil 1 camera!" << std::endl;

            bool dofOn = max3d::get_parameter<bool>( cameraObject, t, _T("dof_on") );
            if( dofOn ) {
                useTargetDistance = false;
                m_focalDistance = max3d::get_parameter<float>( cameraObject, t, _T("dof_focus_distance") );
                m_fStop = max3d::get_parameter<float>( cameraObject, t, _T("dof_fstop") );
                m_focalLength = max3d::get_parameter<float>( cameraObject, t, _T("focal_length") );

                // The focal length is provided in millimeters, so we need to convert it to camera-space units.
                m_focalLength /= (float)frantic::max3d::get_scale_to_millimeters();
            }

            set_horizontal_fov( max3d::get_parameter<float>( cameraObject, t, _T("fov_internal") ) );

            m_projectionMode = static_cast<projection_mode::projection_mode>(
                max3d::get_parameter<int>( cameraObject, t, _T("lens_type") ) );
            if( m_projectionMode > projection_mode::spherical )
                throw std::runtime_error( std::string() + "camera::camera: Brazil camera \"" +
                                          frantic::strings::to_string( cameraNode->GetName() ) +
                                          "\" has an unsupported lens type." );

        } else if( cameraObject->ClassID() == Class_ID( 1532704034, 1337596672 ) ) {
            FF_LOG( debug ) << "It's a Brazil 2 camera!" << std::endl;

            bool dofOn = max3d::get_parameter<bool>( cameraObject, t, _T("Depth_of_Field_PB_Holder.dof_on") );
            if( dofOn ) {
                useTargetDistance = false;
                m_focalDistance =
                    max3d::get_parameter<float>( cameraObject, t, _T("Depth_of_Field_PB_Holder.dof_focus_distance") );
                m_fStop = max3d::get_parameter<float>( cameraObject, t, _T("Depth_of_Field_PB_Holder.dof_fstop") );
                m_focalLength = max3d::get_parameter<float>( cameraObject, t, _T("Lens.focal_length") );

                // The focal length is provided in millimeters, so we need to convert it to camera-space units.
                m_focalLength /= (float)frantic::max3d::get_scale_to_millimeters();
            }

            set_horizontal_fov( frantic::math::degrees_to_radians(
                max3d::get_parameter<float>( cameraObject, t, _T("Lens.fov_internal") ) ) );

            ReferenceTarget* theLens = max3d::get_parameter<ReferenceTarget*>( cameraObject, t, _T("Lens") );

            TSTR lensClass;
            theLens->GetClassName( lensClass );

            std::string lensType = frantic::strings::to_lower( frantic::strings::to_string( lensClass.data() ) );

            if( lensType == "perspective lens" )
                m_projectionMode = projection_mode::perspective;
            else if( lensType == "ortho lens" )
                m_projectionMode = projection_mode::orthographic;
            else if( lensType == "spherical lens" )
                m_projectionMode = projection_mode::spherical;
            else
                throw std::runtime_error( std::string() + "camera::camera: Brazil camera \"" +
                                          frantic::strings::to_string( cameraNode->GetName() ) +
                                          "\" has an unsupported lens type \"" + lensType + "\"." );

        } else if( cameraObject->ClassID() == Class_ID( 1079918574, 1332024254 ) ) {
            FF_LOG( debug ) << "Its a VRay Physical Camera!" << std::endl;

            useTargetDistance = false;

            m_distortionKappa = frantic::max3d::get_parameter<float>( cameraObject, t, _T("distortion") );
            m_fStop = frantic::max3d::get_parameter<float>( cameraObject, t, _T("f_number") );
            m_focalLength = frantic::max3d::get_parameter<float>( cameraObject, t, _T("focal_length") );
            m_focalLength /= (float)frantic::max3d::get_scale_to_millimeters();

            if( frantic::max3d::get_parameter<bool>( cameraObject, t, _T("specify_focus") ) )
                m_focalDistance = frantic::max3d::get_parameter<float>( cameraObject, t, _T("focus_distance") );
            else if( !frantic::max3d::get_parameter<bool>( cameraObject, t, _T("targeted") ) )
                m_focalDistance = frantic::max3d::get_parameter<float>( cameraObject, t, _T("target_distance") );
            else {
                INode* cameraTargetNode = cameraNode->GetTarget();
                if( cameraTargetNode ) {
                    Point3 p = cameraTargetNode->GetNodeTM( t ).GetTrans();
                    m_focalDistance = vector3f::distance( camera_position(), vector3f( p.x, p.y, p.z ) );
                } else
                    m_focalDistance = cameraState.tdist;
            }

            float filmWidth = frantic::max3d::get_parameter<float>( cameraObject, t, _T("film_width") );
            filmWidth /= (float)frantic::max3d::get_scale_to_millimeters();

            // I'm not sure what this is exactly. It seems to be the (3,4) element of a perspective transform
            // matrix with near clipping plane of focalLength and far clipping plane focalDistance.
            float aperatureDistance = ( m_focalDistance * m_focalLength / ( m_focalDistance - m_focalLength ) );

            float tanHalfFov = filmWidth / ( 2.f * aperatureDistance );
            set_horizontal_fov( 2.f * std::atan( tanHalfFov ) );

            int distortionType = frantic::max3d::get_parameter<int>( cameraObject, t, _T("distortion_type") );
            if( distortionType == 0 )
                m_projectionMode = projection_mode::vray_quadratic;
            else if( distortionType == 1 )
                m_projectionMode = projection_mode::vray_cubic;
            else
                throw std::runtime_error( std::string() + "camera::camera: VRay Physical camera \"" +
                                          frantic::strings::to_string( cameraNode->GetName() ) +
                                          "\" has an unsupported distortion type" );
        }

        // If there was no other focal distance, use the target distance
        if( useTargetDistance ) {
            INode* cameraTargetNode = cameraNode->GetTarget();
            if( cameraTargetNode == 0 ) {
                // If this is a free camera, the tdist entry in the CameraState has the target distance.
                m_focalDistance = cameraState.tdist;
            } else {
                // If this camera has a target object, the tdist entry isn't guaranteed to be right,
                // so we should take the distance from the camera object to the target object
                Point3 p = cameraTargetNode->GetNodeTM( t ).GetTrans();
                m_focalDistance = vector3f::distance( camera_position(), vector3f( p.x, p.y, p.z ) );
            }
        }
    }
#endif //_OBJECT_

    static camera from_cube_face( const transform4t<FloatType>& worldTransform,
                                  cube_face::default_cube_face cubeFace ) {
        return camera( worldTransform * transform4t<FloatType>::from_cubeface( cubeFace ), 90 );
    }

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------

    ///////////
    // Getters
    ///////////
  private:
    /**
     * Although camera has been changed so that m_fov's direction depends on the image size, there are still some
     * calculations
     * that requires us to get the tan( m_fov / 2 ) for the horizontal direction.
     */
    FloatType horizontal_tan_half_fov() const {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize <= m_outputSize.ysize * m_pixelAspect ) {
#endif
            return m_tanHalfFov;
        } else {
            return m_tanHalfFov * image_aspect_x_over_y();
        }
    }

  public:
    FloatType fov() const { return m_fov; }

    FloatType fov_degrees() const { return math::radians_to_degrees( fov() ); }

    FloatType horizontal_fov_degrees() const { return math::radians_to_degrees( horizontal_fov() ); }

    FloatType horizontal_fov() const {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize <= m_outputSize.ysize * m_pixelAspect ) {
#endif
            return fov();
        } else {
            using namespace std;

            switch( m_projectionMode ) {
            case projection_mode::panoramic:
            case projection_mode::spherical:
                return m_fov * image_aspect_x_over_y();
            case projection_mode::perspective:
            default: // By default, treat it like a perspective projection
                return (FloatType)( 2.0f * atan( m_tanHalfFov * image_aspect_x_over_y() ) );
            }
        }
    }

    FloatType vertical_fov_degrees() const { return math::radians_to_degrees( vertical_fov() ); }

    FloatType vertical_fov() const {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize < m_outputSize.ysize * m_pixelAspect ) {
#endif
            using namespace std;

            switch( m_projectionMode ) {
            case projection_mode::panoramic:
            case projection_mode::spherical:
                return m_fov * image_aspect_y_over_x();
            case projection_mode::perspective:
            default: // By default, treat it like a perspective projection
                return (FloatType)( 2.0f * atan( m_tanHalfFov * image_aspect_y_over_x() ) );
            }
        } else {
            return fov();
        }
    }

    FloatType orthographic_size() const { return m_orthographicSize; }

    FloatType orthographic_height() const {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize < m_outputSize.ysize * m_pixelAspect ) {
#endif
            return orthographic_size() * image_aspect_y_over_x();
        } else {
            return orthographic_size();
        }
    }

    FloatType orthographic_width() const {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize <= m_outputSize.ysize * m_pixelAspect ) {
#endif
            return orthographic_size();
        } else {
            return orthographic_size() * image_aspect_x_over_y();
        }
    }

    FloatType near_distance() const { return m_near; }
    FloatType far_distance() const { return m_far; }

    FloatType fstop() const { return m_fStop; }

    bool is_static() const { return m_cameraToWorldTransform.is_static(); }

    const vector3t<FloatType>& camera_position() const {
        return m_cameraToWorldTransform.get_transform().get_column( 3 );
    }

    /**
     * This function returns the position of the camera, given a world space position.  For
     * some camera projections, most notably orthographic, the camera is defined by a plane
     * rather than a single point.  In this case, the position on the camera might be different
     * depending on what world position is being considered.  This function takes that into
     * account.
     *
     * @param  worldSpacePos  The world space position under consideration for use when generating the camera position.
     */
    vector3t<FloatType> camera_position( const vector3t<FloatType>& worldSpacePos ) const {
        if( m_projectionMode == projection_mode::orthographic ) {
            vector3t<FloatType> diff = worldSpacePos - camera_position();
            vector3t<FloatType> view_direction = this->view_direction();

            // Compute the distance from the eye-plane, then calculate the intersection
            // of a ray cast from worldSpacePos to the eye-plane.
            FloatType parallelDist = vector3t<FloatType>::dot( diff, view_direction );
            return worldSpacePos + ( -parallelDist * view_direction );
        } else {
            return camera_position();
        }
    }

    // motionSegmentTime should be in the range [0,1] to cover the full camera animation
    vector3t<FloatType> camera_position( FloatType motionSegmentTime ) const {
        return m_cameraToWorldTransform.get_translation( motionSegmentTime );
    }

    vector3t<FloatType> camera_position( FloatType motionSegmentTime, const vector3t<FloatType>& worldSpacePos ) const {
        if( m_projectionMode == projection_mode::orthographic ) {
            vector3t<FloatType> diff = worldSpacePos - camera_position( motionSegmentTime );
            vector3t<FloatType> view_direction = this->view_direction( motionSegmentTime );

            // Compute the distance from the eye-plane, then calculate the intersection
            // of a ray cast from worldSpacePos to the eye-plane.
            FloatType parallelDist = vector3t<FloatType>::dot( diff, view_direction );
            return worldSpacePos + ( -parallelDist * view_direction );
        } else {
            return camera_position( motionSegmentTime );
        }
    }

    // The orientation of the camera as a quaternion
    quat4t<FloatType> camera_orientation() const {
        return quat4t<FloatType>::from_transform4f( m_cameraToWorldTransform.get_transform() );
    }

    quat4t<FloatType> camera_orientation( FloatType motionSegmentTime ) const {
        return quat4t<FloatType>::from_transform4f( m_cameraToWorldTransform.get_transform( motionSegmentTime ) );
    }

    // The direction the camera is looking
    vector3t<FloatType> view_direction() const { return -m_cameraToWorldTransform.get_transform().get_column( 2 ); }

    vector3t<FloatType> view_direction( FloatType motionSegmentTime ) const {
        return -m_cameraToWorldTransform.get_transform_column( motionSegmentTime, 2 );
    }

    // The "up" direction of the camera
    vector3t<FloatType> view_up() const { return m_cameraToWorldTransform.get_transform().get_column( 1 ); }

    vector3t<FloatType> view_up( FloatType motionSegmentTime ) const {
        return m_cameraToWorldTransform.get_transform_column( motionSegmentTime, 1 );
    }

    // The "right" direction of the camera
    vector3t<FloatType> view_right() const { return m_cameraToWorldTransform.get_transform().get_column( 0 ); }

    vector3t<FloatType> view_right( FloatType motionSegmentTime ) const {
        return m_cameraToWorldTransform.get_transform_column( motionSegmentTime, 0 );
    }

    // Returns the camera to world transform
    const transform4t<FloatType>& world_transform() const { return m_cameraToWorldTransform.get_transform(); }

    // Returns the camera to world transform
    transform4t<FloatType> world_transform( FloatType motionSegmentTime ) const {
        return m_cameraToWorldTransform.get_transform( motionSegmentTime );
    }

    // Returns the world to camera transform
    const transform4t<FloatType>& world_transform_inverse() const {
        return m_cameraToWorldTransform.get_inverse_transform();
    }

    // Returns the world to camera transform
    transform4t<FloatType> world_transform_inverse( FloatType motionSegmentTime ) const {
        return m_cameraToWorldTransform.get_inverse_transform( motionSegmentTime );
    }

    // Returns the projection transform (only valid for orthographic and perspective)
    //
    // Assumes screen space is a cube with corners (-1,-1,-1) and (1,1,1).  positions outside are considered not
    // visible.
    //
    // This returns a transform that is suitable for use in OpenGL which assumes z is in [-1,1].  However, it is not
    // numerically  precise if the clipping range is wide and depending on the application, developers should use
    // alternative methods such as get_worldspace_ray instead (or at least modify the z range to be [0,1] instead.
    //
    transform4t<FloatType> get_projection_transform() const {
        switch( m_projectionMode ) {
        case projection_mode::orthographic: {
            return transform4t<FloatType>( 2 / orthographic_width(), 0, 0, 0, 0, 2 / orthographic_height(), 0, 0, 0, 0,
                                           -2 / ( m_far - m_near ), 0, 0, 0, -2 * m_near / ( m_far - m_near ) - 1, 1 );
        }
        case projection_mode::perspective: {
            const FloatType cotFovX = math::cot( horizontal_fov() / 2 );
            const FloatType cotFovY = math::cot( vertical_fov() / 2 );

            return transform4t<FloatType>( cotFovX, 0, 0, 0, 0, cotFovY, 0, 0, 0, 0,
                                           ( -m_near - m_far ) / ( m_far - m_near ), -1, 0, 0,
                                           -2 * m_far * m_near / ( m_far - m_near ), 0 );
        }
        default:
            break;
        }
        return transform4t<FloatType>();
    }

    bool has_dof() const {
        // TODO: probably need a better criteria?
        return m_fStop < 100.f;
    }

    FloatType pixel_aspect() const { return m_pixelAspect; }

    frantic::graphics2d::size2 get_output_size() const { return m_outputSize; }

    frantic::graphics2d::vector2f get_subpixel_offset() const { return m_subpixelOffset; }

    // Returns a set of 6 cameras, facing in all 6 directions, parallel to an axis, based on the transform matrix of the
    // calling camera
    // order is: Z+, Z-, Y-, Y+, X-, X+
    std::vector<camera<FloatType>> get_cubic_cameras( int outputSizeSide ) const {
        using namespace frantic::graphics::cube_face;
        std::vector<camera<FloatType>> result;

        result.push_back( from_cube_face( world_transform(), CF_Z_POS ) );
        result.push_back( from_cube_face( world_transform(), CF_Z_NEG ) );
        result.push_back( from_cube_face( world_transform(), CF_Y_NEG ) );
        result.push_back( from_cube_face( world_transform(), CF_Y_POS ) );
        result.push_back( from_cube_face( world_transform(), CF_X_NEG ) );
        result.push_back( from_cube_face( world_transform(), CF_X_POS ) );

        for( typename std::vector<camera<FloatType>>::iterator it = result.begin(); it != result.end(); ++it ) {
            it->set_output_size( frantic::graphics2d::size2( outputSizeSide ) );
            it->set_fov_degrees( 90 );
            it->set_pixel_aspect( 1.0f );

            it->set_near( m_near );
            it->set_far( m_far );
            it->set_fstop( m_fStop );
            it->set_focal_length( m_focalLength );
            it->set_focal_distance( m_focalDistance );
            it->set_subpixel_offset( m_subpixelOffset );
        }

        return result;
    }

    ///////////
    // Setters
    ///////////
  public:
    void set_transform( const motion_blurred_transform<FloatType>& trans ) { m_cameraToWorldTransform = trans; }

    /**
     * Set the camera's transform.
     *
     * @note Overrides camera animation
     */
    void set_transform( const transform4t<FloatType>& t ) { m_cameraToWorldTransform = t; }

    /**
     * Set the camera's position.
     *
     * @note Overrides camera animation
     */
    void set_position( const frantic::graphics::vector3t<FloatType>& newPosition ) {
        frantic::graphics::transform4t<FloatType> cameraTrans = world_transform();
        cameraTrans.set_translation( newPosition );
        m_cameraToWorldTransform = cameraTrans;
    }

    /**
     * Set the camera's orientation.
     *
     * @note Overrides camera animation
     */
    void set_orientation( const frantic::graphics::quat4t<FloatType>& newRotation ) {
        frantic::graphics::vector3t<FloatType> oldPosition;
        frantic::graphics::quat4t<FloatType> oldRotation;
        world_transform().decompose( oldPosition, oldRotation );

        const frantic::graphics::transform4t<FloatType> cameraTrans =
            frantic::graphics::transform4t<FloatType>::from_translation( oldPosition ) * newRotation.to_transform4f();
        set_transform( cameraTrans );
    }

    void set_fov( FloatType fovRad ) {
        m_fov = fovRad;
        m_tanHalfFov = tan( m_fov / 2 );
    }

    void set_fov_degrees( FloatType fovDeg ) { set_fov( math::degrees_to_radians( fovDeg ) ); }

    void set_horizontal_fov( FloatType horizontalFov ) {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize <= m_outputSize.ysize * m_pixelAspect ) {
#endif
            set_fov( horizontalFov );
        } else {
            using namespace std;

            FloatType verticalFov;
            switch( m_projectionMode ) {
            case projection_mode::panoramic:
            case projection_mode::spherical:
                verticalFov = horizontalFov * image_aspect_y_over_x();
                break;
            case projection_mode::perspective:
            default: // By default, treat it like a perspective projection
                verticalFov = (FloatType)( 2.0f * atan( tan( horizontalFov * 0.5f ) * image_aspect_y_over_x() ) );
            }
            set_fov( verticalFov );
        }
    }

    void set_horizontal_fov_degrees( FloatType horizontalFovDeg ) {
        set_horizontal_fov( math::degrees_to_radians( horizontalFovDeg ) );
    }

    void set_vertical_fov( FloatType verticalFov ) {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize < m_outputSize.ysize * m_pixelAspect ) {
#endif
            using namespace std;

            FloatType horizontalFov;
            switch( m_projectionMode ) {
            case projection_mode::panoramic:
            case projection_mode::spherical:
                horizontalFov = verticalFov * image_aspect_x_over_y();
            case projection_mode::perspective:
            default: // By default, treat it like a perspective projection
                horizontalFov = (FloatType)( 2.0f * atan( tan( verticalFov * 0.5f ) * image_aspect_x_over_y() ) );
            }

            set_fov( horizontalFov );
        } else {
            set_fov( verticalFov );
        }
    }

    void set_vertical_fov_degrees( FloatType verticalFovDeg ) {
        set_vertical_fov( math::degrees_to_radians( verticalFovDeg ) );
    }

    void set_orthographic_size( FloatType orthographicSize ) { m_orthographicSize = orthographicSize; }

    void set_orthographic_width( FloatType orthographicWidth ) {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize <= m_outputSize.ysize * m_pixelAspect ) {
#endif
            set_orthographic_size( orthographicWidth );
        } else {
            set_orthographic_size( orthographicWidth * image_aspect_y_over_x() );
        }
    }

    void set_orthographic_height( FloatType orthographicHeight ) {
#if defined( KVRY_MODE )
        if( true ) {
#else
        if( m_outputSize.xsize < m_outputSize.ysize * m_pixelAspect ) {
#endif
            set_orthographic_size( orthographicHeight * image_aspect_x_over_y() );
        } else {
            set_orthographic_size( orthographicHeight );
        }
    }

    void set_output_size( const frantic::graphics2d::size2& outputSize ) { m_outputSize = outputSize; }

    void set_pixel_aspect( FloatType aspect ) { m_pixelAspect = aspect; }

    void set_near( FloatType nearDist ) { m_near = nearDist; }

    void set_far( FloatType farDist ) { m_far = farDist; }

    void set_projection_mode( projection_mode::projection_mode projectionMode ) { m_projectionMode = projectionMode; }

    projection_mode::projection_mode projection_mode() const { return m_projectionMode; }

    std::string projection_mode_string() const {
        switch( m_projectionMode ) {
        case projection_mode::perspective:
            return "perspective";
        case projection_mode::orthographic:
            return "orthographic";
        case projection_mode::spherical:
            return "spherical";
        case projection_mode::panoramic:
            return "panoramic";
        case projection_mode::vray_quadratic:
            return "VRay quadratic distortion";
        case projection_mode::vray_cubic:
            return "VRay cubic distortion";
        default:
            return "unknown";
        }
    }

    void set_fstop( FloatType fStop ) { m_fStop = fStop; }

    void set_focal_distance( FloatType dist ) { m_focalDistance = dist; }

    void set_focal_length( FloatType len ) { m_focalLength = len; }

    void set_subpixel_offset( frantic::graphics2d::vector2f offset ) { m_subpixelOffset = offset; }

    /**
     * Rotate the camera by `upAmount` radians around the camera's `view_right()`, then `rightAmount` radians around the
     * world up. The position of the camera is preserved.
     *
     * Assuming the axis is facing the viewer, the rotation around `view_right()` is clockwise, and the rotation around
     * the world up is counterclockwise as long as the world up points above the camera's perceived horizon, otherwise
     * it rotates clockwise.
     *
     * @returns The rotation applied.
     */
    quat4t<FloatType> rotate_camera( FloatType upAmount, FloatType rightAmount, const vector3t<FloatType>& worldUp ) {
        const quat4t<FloatType> startRotation = camera_orientation();

        // If the view is upside-down, reverse horizontal rotation to keep things feeling right
        double angle = std::acos( vector3t<FloatType>::dot( worldUp, view_up() ) /
                                  ( worldUp.get_magnitude() * view_up().get_magnitude() ) );
        if( angle > FloatType( M_PI / 2 ) ) {
            rightAmount = -rightAmount;
        }

        const quat4t<FloatType> rotation = quat4t<FloatType>::from_angle_axis( rightAmount, worldUp ) *
                                           quat4t<FloatType>::from_angle_axis( upAmount, view_right() );

        set_orientation( rotation * startRotation );

        return camera_orientation() * startRotation.get_inverse();
    }

    /**
     * Similar to `rotate_camera` except that position is changed as if the camera orbited around `center`.
     *
     * @returns The rotation applied and the new camera position.
     */
    std::pair<frantic::graphics::quat4t<FloatType>, frantic::graphics::vector3t<FloatType>>
    orbit_camera( const frantic::graphics::vector3t<FloatType>& center, FloatType upAmount, FloatType rightAmount,
                  const vector3t<FloatType>& worldUp ) {
        const frantic::graphics::quat4t<FloatType> rotation = rotate_camera( upAmount, rightAmount, worldUp );

        const frantic::graphics::vector3t<FloatType> cameraPosition = camera_position();
        const frantic::graphics::vector3t<FloatType> position =
            rotation.to_transform4f() * ( cameraPosition - center ) + center;
        set_position( position );

        return std::make_pair( rotation, position );
    }

    ///////////
    // Projection functions
    ///////////

    // Given a pixel coordinate, returns the direction in camera space
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    vector3t<FloatType> to_cameraspace_direction( const frantic::graphics2d::vector2f& p,
                                                  bool& outIsValidInput ) const {
        switch( m_projectionMode ) {
        case projection_mode::spherical:
            return vector3t<FloatType>::from_fisheye( p - m_subpixelOffset, frantic::graphics2d::size2f( m_outputSize ),
                                                      horizontal_fov(), outIsValidInput, m_pixelAspect );
        // TODO: implement this stuff
        case projection_mode::panoramic:
            outIsValidInput = true;
            return from_panoramic_yup( p - m_subpixelOffset, frantic::graphics2d::size2f( m_outputSize ),
                                       horizontal_fov(), m_pixelAspect );
        case projection_mode::orthographic:
            outIsValidInput = true;
            return vector3t<FloatType>( 0, 0, -1 );
        case projection_mode::vray_quadratic: {
            FloatType xfactor = 2.f * m_focalDistance * horizontal_tan_half_fov();
            FloatType yfactor = xfactor * image_aspect_y_over_x();

            FloatType nx = xfactor * ( p.x / m_outputSize.xsize - 0.5f );
            FloatType ny = yfactor * ( p.y / m_outputSize.ysize - 0.5f );
            FloatType nz = -std::sqrt( m_focalDistance * m_focalDistance - m_distortionKappa * ( nx * nx + ny * ny ) );

            return vector3t<FloatType>( nx, ny, nz );
        }
        case projection_mode::perspective:
        default:
            return vector3t<FloatType>::from_perspective_projection( p - m_subpixelOffset,
                                                                     frantic::graphics2d::size2f( m_outputSize ),
                                                                     horizontal_tan_half_fov(), m_pixelAspect );
        }
    }

    // Given a direction in camera space, returns a pixel coordinate
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2f from_cameraspace_direction( vector3t<FloatType> v, bool& outIsValidInput ) const {
        switch( m_projectionMode ) {
        case projection_mode::spherical:
            outIsValidInput = true;
            return m_subpixelOffset + v.to_fisheye( frantic::graphics2d::size2f( m_outputSize ), horizontal_fov() );
        case projection_mode::panoramic:
            outIsValidInput = true;
            return m_subpixelOffset +
                   to_panoramic_yup( v, frantic::graphics2d::size2f( m_outputSize ), horizontal_fov(), m_pixelAspect );
        case projection_mode::orthographic:
            outIsValidInput = false;
            return frantic::graphics2d::vector2f( 0, 0 );
        case projection_mode::perspective:
        default:
            return m_subpixelOffset + v.to_perspective_projection( frantic::graphics2d::size2f( m_outputSize ),
                                                                   horizontal_tan_half_fov(), outIsValidInput,
                                                                   m_pixelAspect );
        }
    }

    // Given a position in camera space, returns a pixel coordinate
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2f from_cameraspace_position( const vector3t<FloatType>& v,
                                                             bool& outIsValidInput ) const {
        using namespace std;

        switch( m_projectionMode ) {
        case projection_mode::spherical:
            return m_subpixelOffset + v.to_fisheye( frantic::graphics2d::size2f( m_outputSize ), horizontal_fov() );
        case projection_mode::panoramic:
            outIsValidInput = true;
            return m_subpixelOffset +
                   to_panoramic_yup( v, frantic::graphics2d::size2f( m_outputSize ), horizontal_fov(), m_pixelAspect );
        case projection_mode::orthographic: {
            FloatType xFactorInv, yFactorInv;
#if defined( KVRY_MODE )
            if( true ) {
#else
            if( m_outputSize.xsize < m_outputSize.ysize * m_pixelAspect ) {
#endif
                xFactorInv = m_outputSize.xsize / m_orthographicSize;
                yFactorInv = m_outputSize.xsize / ( m_orthographicSize * m_pixelAspect );
            } else {
                xFactorInv = m_pixelAspect * m_outputSize.ysize / m_orthographicSize;
                yFactorInv = m_outputSize.ysize / m_orthographicSize;
            }

            FloatType rx = m_subpixelOffset.x + v.x * xFactorInv + m_outputSize.xsize * 0.5f;
            FloatType ry = m_subpixelOffset.y + v.y * yFactorInv + m_outputSize.ysize * 0.5f;
            return frantic::graphics2d::vector2f( static_cast<float>( rx ), static_cast<float>( ry ) );
        }
        case projection_mode::vray_quadratic:
            if( v.z < 0 ) {
                FloatType t = m_focalDistance / sqrt( ( v.x * v.x + v.y * v.y ) * m_distortionKappa + v.z * v.z );
                FloatType focalDistanceTimesTanHalfFov = m_focalDistance * m_tanHalfFov;

                FloatType pX, pY;
#if defined( KVRY_MODE )
                if( true ) {
#else
                if( m_outputSize.xsize <= m_outputSize.ysize * m_pixelAspect ) {
#endif
                    pX = v.x * t / focalDistanceTimesTanHalfFov;
                    pY = v.y * t / ( focalDistanceTimesTanHalfFov * image_aspect_y_over_x() );
                } else {
                    pX = v.x * t * image_aspect_y_over_x() / focalDistanceTimesTanHalfFov;
                    pY = v.y * t / focalDistanceTimesTanHalfFov;
                }

                return frantic::graphics2d::vector2f(
                    float( m_subpixelOffset.x + 0.5f * ( 1.f + pX ) * m_outputSize.xsize ),
                    float( m_subpixelOffset.y + 0.5f * ( 1.f + pY ) * m_outputSize.ysize ) );
            }

            outIsValidInput = false;
            return frantic::graphics2d::vector2f( -1, -1 );
        case projection_mode::perspective:
        default:
            return m_subpixelOffset + v.to_perspective_projection( frantic::graphics2d::size2f( m_outputSize ),
                                                                   static_cast<float>( horizontal_tan_half_fov() ),
                                                                   outIsValidInput,
                                                                   static_cast<float>( m_pixelAspect ) );
        }
    }

    // Given a position in camera space, returns a pixel coordinate
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    template <class RandomNumberGenerator>
    frantic::graphics2d::vector2f from_cameraspace_position_with_dof_jitter( vector3t<FloatType> v,
                                                                             RandomNumberGenerator& rng,
                                                                             bool& outIsValidInput ) const {
        using frantic::graphics2d::vector2f;

        vector2f pixel = from_cameraspace_position( v, outIsValidInput );
        if( has_dof() ) {

            // this was for square pixels
            // float pixelRadius = dof_circle_of_confusion_pixel_radius( -v.z );
            // return pixel + pixelRadius * vector2f::from_unit_disk_random( rng );

            vector2f pixelRadii;
            dof_ellipse_of_confusion_pixel_radii( -v.z, pixelRadii.x, pixelRadii.y );
            return pixel + pixelRadii * vector2f::from_unit_disk_random( rng );

        } else {
            return pixel;
        }
    }

    // Given a pixel coordinate, returns the direction in world space
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    vector3t<FloatType> to_worldspace_direction( const frantic::graphics2d::vector2f& p, bool& outIsValidInput ) const {
        return world_transform().transform_no_translation( to_cameraspace_direction( p, outIsValidInput ) );
    }

    // Given a direction in world space, returns a pixel coordinate
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2f from_worldspace_direction( const vector3t<FloatType>& v,
                                                             bool& outIsValidInput ) const {
        vector3t<FloatType> cameraSpaceV = world_transform_inverse().transform_no_translation( v );
        return from_cameraspace_direction( cameraSpaceV, outIsValidInput );
    }

    // Given a direction in world space, returns a pixel coordinate
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2f from_worldspace_direction( const vector3t<FloatType>& v, FloatType motionSegmentTime,
                                                             bool& outIsValidInput ) const {
        vector3t<FloatType> cameraSpaceV = world_transform_inverse( motionSegmentTime ).transform_no_translation( v );
        return from_cameraspace_direction( cameraSpaceV, outIsValidInput );
    }

    // Given a position in world space, returns a pixel coordinate
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2f from_worldspace_position( const vector3t<FloatType>& p,
                                                            bool& outIsValidInput ) const {
        vector3t<FloatType> cameraSpaceP = world_transform_inverse() * p;
        return from_cameraspace_position( cameraSpaceP, outIsValidInput );
    }

    // Given a position in world space, returns a pixel coordinate
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics2d::vector2f from_worldspace_position( const vector3t<FloatType>& p, FloatType motionSegmentTime,
                                                            bool& outIsValidInput ) const {
        vector3t<FloatType> cameraSpaceP = world_transform_inverse( motionSegmentTime ) * p;
        return from_cameraspace_position( cameraSpaceP, outIsValidInput );
    }

    // Given a position in image space, returns the position/direction in world space
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics::ray3t<FloatType> get_cameraspace_ray( const frantic::graphics2d::vector2f& p,
                                                             bool& outIsValidInput ) const {
        if( m_projectionMode == projection_mode::orthographic ) {
            FloatType xFactor, yFactor;
#if defined( KVRY_MODE )
            if( true ) {
#else
            if( m_outputSize.xsize < m_outputSize.ysize * m_pixelAspect ) {
#endif
                xFactor = m_orthographicSize / m_outputSize.xsize;
                yFactor = m_orthographicSize * m_pixelAspect / m_outputSize.xsize;
            } else {
                xFactor = m_orthographicSize / ( m_pixelAspect * m_outputSize.ysize );
                yFactor = m_orthographicSize / m_outputSize.ysize;
            }

            return ray3t<FloatType>(
                vector3t<FloatType>( ( p.x - m_subpixelOffset.x - m_outputSize.xsize * 0.5f ) * xFactor,
                                     ( p.y - m_subpixelOffset.y - m_outputSize.ysize * 0.5f ) * yFactor, 0 ),
                vector3f( 0, 0, -1 ) );
        } else {
            return ray3t<FloatType>( vector3t<FloatType>(), to_cameraspace_direction( p, outIsValidInput ) );
        }
    }

    // Given a position in image space, returns the position/direction in world space
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics::ray3t<FloatType> get_worldspace_ray( const frantic::graphics2d::vector2f& p,
                                                            bool& outIsValidInput ) const {
        return world_transform() * get_cameraspace_ray( p, outIsValidInput );
    }

    // Given a position in image space, returns the position/direction in world space
    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics::ray3t<FloatType> get_worldspace_ray( const frantic::graphics2d::vector2f& p,
                                                            FloatType motionSegmentTime, bool& outIsValidInput ) const {
        return world_transform( motionSegmentTime ) * get_cameraspace_ray( p, outIsValidInput );
    }

    // Imagine, if you will, a camera-facing unit cube. This function will return the
    // area, in pixels, of said cube, at position p, projected onto the screen.
    FloatType area_differential_from_cameraspace_position( const frantic::graphics::vector3t<FloatType>& p ) const {
        using namespace std;

        FloatType area = 0.f;
        switch( m_projectionMode ) {
        case projection_mode::orthographic: {
            FloatType dx = m_outputSize.xsize / orthographic_width();
            FloatType dy = dx * m_pixelAspect;
            return dx * dy;
        }
        case projection_mode::spherical: {
            // This is grossly ugly
            // TODO: find a better method for computing the spherical projection area differential
            vector3t<FloatType> localAxisX = vector3t<FloatType>::cross( p, vector3t<FloatType>( 0, 1, 0 ) );
            vector3t<FloatType> localAxisY = vector3t<FloatType>::cross( p, localAxisX );

            localAxisX.normalize();
            localAxisY.normalize();

            bool bTest = true;
            frantic::graphics2d::vector2f p1 = from_cameraspace_position( p + localAxisX, bTest );
            frantic::graphics2d::vector2f p2 = from_cameraspace_position( p + localAxisY, bTest );
            frantic::graphics2d::vector2f p3 = from_cameraspace_position( p, bTest );

            return sqrt( frantic::graphics2d::vector2f::distance_squared( p1, p3 ) *
                         frantic::graphics2d::vector2f::distance_squared( p2, p3 ) );
        }
        case projection_mode::perspective: // Intentional fall-through
        default: {
            FloatType tanHalfFov = horizontal_tan_half_fov();
            /*area  = ( m_outputSize.xsize / tanHalfFov ); // *( m_pixelAspect * m_outputSize.xsize / tanHalfFov );
            area = area * area * m_pixelAspect;
            area *= (0.25f / p.get_magnitude_squared());
            return area;*/

            FloatType numer = (FloatType)( m_outputSize.xsize * m_outputSize.xsize ) * m_pixelAspect;
            FloatType denom = 4.f * p.get_magnitude_squared() * ( tanHalfFov * tanHalfFov );
            return numer / denom;
        }
        }
    }

    // NOTE: outIsValidInput should be true when it's passed in, and is set to false if it fails.
    //       It remains false if it was false when passed in.
    frantic::graphics::vector3t<FloatType>
    get_point_projected_to_plane( frantic::graphics2d::vector2f p,
                                  const frantic::graphics::plane3t<FloatType>& projectionPlane,
                                  bool& outIsValidInput ) const {
        frantic::graphics::ray3t<FloatType> projectionRay = get_worldspace_ray( p, outIsValidInput );
        if( outIsValidInput ) {
            double d = projectionPlane.get_distance_to_intersection( projectionRay );
            if( d >= 0 ) {
                return projectionRay.at( d );
            } else {
                outIsValidInput = false;
                return vector3t<FloatType>();
            }
        } else {
            outIsValidInput = false;
            return vector3t<FloatType>();
        }
    }

    // This returns the distance that the camera path travels in terms of rotation, in radians
    FloatType get_pixel_motion_distance( frantic::graphics2d::vector2f p, frantic::graphics2d::size2f s ) const {
        return motion_blurred_transform<FloatType>::get_pixel_motion_distance( m_cameraToWorldTransform, p, s,
                                                                               horizontal_fov() );
    }

    /**
     * Helper method to get the appropriate value such that when the camera is placed at
     *
     *     camera_position() - <return value> * view_direction()
     *
     * the point at `objectToView` will be visible.
     *
     * Takes into account FOV but not clipping.
     *
     * @param objectToView The point that should appear within the view frustum after the above transformation.
     * @param epsilon The acceptable margin of error. Currently this only affects the case where the camera is already
     *                really close to the objectToView, in which 0 is returned.
     * @throws std::runtime_error If this method is called on a non-perspective camera.
     */
    FloatType get_suggested_depth( const frantic::graphics::vector3t<FloatType>& objectToView,
                                   FloatType epsilon = std::numeric_limits<FloatType>::epsilon() ) const {
        if( m_projectionMode != projection_mode::perspective )
            throw std::runtime_error(
                "frantic::graphics::camera.get_suggested_depth was called on a non-perspective camera" );

        const frantic::graphics::vector3t<FloatType> centerToObject = objectToView - camera_position();
        const FloatType centerToObjectDistance = centerToObject.get_magnitude();

        if( centerToObjectDistance <= epsilon )
            return 0;

        const FloatType cosA =
            frantic::graphics::vector3t<FloatType>::dot( centerToObject, -view_direction() ) / centerToObjectDistance;
        const FloatType centerAngle = std::acos( cosA );
        const FloatType sinA = std::sin( centerAngle );

        // `centerToObjectDistance * cosA` is the amount we need to move the camera back by so that both the
        // `objectToView` and the camera position lie on a plane with a normal of `view_direction()`.
        //
        // `centerToObjectDistance * sinA / m_tanHalfFov` is the additional amount we need to move the camera back by so
        // that the point falls within the viewing frustum.
        return centerToObjectDistance * ( cosA + sinA / m_tanHalfFov );
    }

    /**
     * Helper method to get the appropriate value such that when the camera is placed at
     *
     *     camera_position() - <return value> * view_direction()
     *
     * the bounding box `objectToView` will be entirely visible.
     *
     * @throws std::runtime_error If this method is called on a non-perspective camera.
     */
    FloatType get_suggested_depth( const frantic::graphics::boundbox3t<FloatType>& objectToView ) const {
        if( m_projectionMode != projection_mode::perspective )
            throw std::runtime_error(
                "frantic::graphics::camera.get_suggested_depth was called on a non-perspective camera" );

        double maxDistance = -std::numeric_limits<FloatType>::infinity();

        for( int i = 0; i < 8; ++i ) {
            const frantic::graphics::vector3t<FloatType> testPosition = objectToView.get_corner( i );
            const FloatType testDistance = get_suggested_depth( testPosition );
            maxDistance = std::max( maxDistance, testDistance );
        }

        return maxDistance;
    }

    // The maximum distance will be at one of the corners.
    FloatType get_maximum_pixel_motion_distance( frantic::graphics2d::size2 s ) {
        FloatType result =
            get_pixel_motion_distance( frantic::graphics2d::vector2f( 0, 0 ), frantic::graphics2d::size2f( s ) );
        FloatType test = get_pixel_motion_distance( frantic::graphics2d::vector2f( (float)s.xsize - 1, 0 ),
                                                    frantic::graphics2d::size2f( s ) );
        if( test > result )
            result = test;
        test = get_pixel_motion_distance( frantic::graphics2d::vector2f( 0, (float)s.ysize - 1 ),
                                          frantic::graphics2d::size2f( s ) );
        if( test > result )
            result = test;
        test = get_pixel_motion_distance( frantic::graphics2d::vector2f( (float)s.xsize - 1, (float)s.ysize - 1 ),
                                          frantic::graphics2d::size2f( s ) );
        if( test > result )
            result = test;
        return result;
    }

    // This approximates the number of radians taken up by the pixel distance given
    FloatType pixels_to_radians( FloatType pixels ) const {
        using namespace std;
        // TODO: different value for spherical projection
        return atan( 2 * pixels / output_size().xsize * horizontal_tan_half_fov() );
    }

    FloatType vertical_pixels_to_radians( FloatType pixels ) const {
        using namespace std;
        return atan( 2 * pixels / output_size().ysize * tan( vertical_fov() / 2 ) );
    }

    FloatType get_horizontal_angle_from_pixel( FloatType pixels ) const {
        return ( 2 * pixels / output_size().xsize ) * horizontal_fov() / 2;
    }

    FloatType get_vertical_angle_from_pixel( FloatType pixels ) const {
        return ( 2 * pixels / output_size().ysize ) * vertical_fov() / 2;
    }

    const frantic::graphics2d::size2& output_size() const { return m_outputSize; }

    FloatType focal_distance() const { return m_focalDistance; }

    FloatType focal_length() const { return m_focalLength; }

    frantic::graphics2d::size2f output_size2f() const {
        return frantic::graphics2d::size2f( (float)m_outputSize.xsize, (float)m_outputSize.ysize );
    }

    // Returns the image aspect ratio for converting from horizontal to vertical
    // For example, to get the vertical field of view given the horizontal field of view for spherical and panoramic as
    // well as orthographic height from width use:
    //   vfov = hfov * image_aspect_y_over_x();
    // For perspective, use:
    //   vfov = 2 * atan( tan( hfov / 2 ) * image_aspect_y_over_x() );
    FloatType image_aspect_y_over_x() const { return m_pixelAspect * m_outputSize.ysize / m_outputSize.xsize; }

    // Returns the image aspect ratio for converting from vertical to horizontal
    // For example, to get the horizontal field of view given the vertical field of view for spherical and panoramic, as
    // well as orthographic width from height use:
    //   hfov = vfov * image_aspect_x_over_y();
    // For perspective, use:
    //   hfov = 2 * atan( tan( vfov / 2 ) * image_aspect_x_over_y() );
    FloatType image_aspect_x_over_y() const { return m_outputSize.xsize / ( m_outputSize.ysize * m_pixelAspect ); }

    //////////////////////////////
    // Depth of Field Functions
    //////////////////////////////

    // Returns the radius of the circle of confusion at given distance along Z, measured
    // in pixels on the final image.
    FloatType dof_circle_of_confusion_pixel_radius( FloatType zDistance ) const {
        if( has_dof() ) {
            // This computes the radius assuming a perspective projection.  Except for really
            // wide projections, this is probably fine for spherical lenses too.
            FloatType lensRadius = 0.5f * m_focalLength / m_fStop;
            FloatType imageDistance = 1 / ( 1 / m_focalLength - 1 / m_focalDistance );
            FloatType objectImageDistance = 1 / ( 1 / m_focalLength - 1 / zDistance );
            FloatType circleOfConfusionWorldRadius =
                lensRadius * fabsf( imageDistance - objectImageDistance ) / objectImageDistance;
            return circleOfConfusionWorldRadius * m_outputSize.xsize /
                   ( 2 * imageDistance * horizontal_tan_half_fov() );
        } else {
            return 0;
        }
    }

    // Returns the radii of the ellipse of confusion at given distance along Z, measured
    // in pixels on the final image.
    void dof_ellipse_of_confusion_pixel_radii( FloatType zDistance, float& outXradius, float& outYradius ) const {
        using namespace std;

        if( has_dof() ) {
            // This computes the radii of the ellipse assuming a perspective projection.  Except for really
            // wide projections, this is probably fine for spherical lenses too.
            FloatType lensRadius = 0.5f * m_focalLength / m_fStop;
            FloatType imageDistance = 1 / ( 1 / m_focalLength - 1 / m_focalDistance );
            FloatType objectImageDistance = 1 / ( 1 / m_focalLength - 1 / zDistance );
            FloatType circleOfConfusionWorldRadius =
                lensRadius * fabsf( imageDistance - objectImageDistance ) / objectImageDistance;
            FloatType Xradius =
                circleOfConfusionWorldRadius * m_outputSize.xsize / ( 2 * imageDistance * horizontal_tan_half_fov() );
            Xradius = m_focalLength * m_focalLength * fabs( m_focalDistance - zDistance ) * m_outputSize.xsize;
            Xradius /= 2 * m_fStop * ( m_focalDistance - m_focalLength ) * zDistance *
                       ( 2 * imageDistance * horizontal_tan_half_fov() );
            outXradius = float( Xradius );
            outYradius = float( outXradius * m_pixelAspect );
        } else {
            outXradius = 0;
            outYradius = 0;
        }
    }

    // XML Output function
    void write_xml( std::ostream& out, const std::string& prefix = "" ) const {
        out << prefix << "<xfovdegrees>" << horizontal_fov_degrees() << "</xfovdegrees>\n";
        out << prefix << "<orthographicwidth>" << orthographic_width() << "</orthographicwidth>\n";
        out << prefix << "<pixelaspect>" << m_pixelAspect << "</pixelaspect>\n";
        out << prefix << "<imagewidth>" << output_size().xsize << "</imagewidth>\n";
        out << prefix << "<imageheight>" << output_size().ysize << "</imageheight>\n";
        out << prefix << "<nearclip>" << near_distance() << "</nearclip>\n";
        out << prefix << "<farclip>" << far_distance() << "</farclip>\n";
        out << prefix << "<projectionmode>" << projection_mode_string() << "</projectionmode>\n";
        if( has_dof() ) {
            out << prefix << "<dof>\n";
            out << prefix << " <fstop>" << m_fStop << "</fstop>\n";
            out << prefix << " <focallength>" << m_focalLength << "</focallength>\n";
            out << prefix << " <focaldistance>" << m_focalDistance << "</focaldistance>\n";
            out << prefix << "</dof>\n";
        }
        m_cameraToWorldTransform.write_xml( out, prefix );
    }
};

template <class T>
inline void from_panoramic_yup( const T* v2, const T* size, T hfov, T pixelAspect, T* v3Out ) {
    T xNorm = v2[0] / size[0];
    T yNorm = v2[1] / size[1];

    T xHalfRange = hfov / T( 2.0 );
    T yHalfRange = ( hfov * ( size[1] / ( pixelAspect * size[0] ) ) ) / T( 2.0 );

    const T tPi = T( M_PI );

    T xLoRange = tPi + xHalfRange;
    T xHiRange = tPi - xHalfRange;
    T yLoRange = ( tPi / T( 2.0 ) ) + yHalfRange;
    T yHiRange = ( tPi / T( 2.0 ) ) - yHalfRange;

    T phi = ( ( T( 1.0 ) - xNorm ) * xLoRange ) + ( xNorm * xHiRange );
    T theta = ( ( T( 1.0 ) - yNorm ) * yLoRange ) + ( yNorm * yHiRange );

    v3Out[0] = sin( theta ) * sin( phi );
    v3Out[1] = cos( theta );
    v3Out[2] = sin( theta ) * cos( phi );
}

template <class T>
inline void to_panoramic_yup( const T* v3, const T* size, T hfov, T pixelAspect, T* v2Out ) {

    T r = sqrt( frantic::math::square( v3[0] ) + frantic::math::square( v3[1] ) + frantic::math::square( v3[2] ) );
    T phi = T( 0.0 );

    if( !( v3[0] == T( 0.0 ) && v3[2] == T( 0.0 ) ) ) {
        phi = atan2( v3[0], -v3[2] );
    }

    T theta = -acos( v3[1] / r );

    T vfov = ( hfov * ( size[1] / ( pixelAspect * size[0] ) ) );

    v2Out[0] = ( ( phi / hfov ) + T( 0.5 ) ) * size[0];
    v2Out[1] = ( ( ( theta + T( M_PI / 2.0 ) ) / vfov ) + T( 0.5 ) ) * size[1];
}
} // namespace graphics
} // namespace frantic
