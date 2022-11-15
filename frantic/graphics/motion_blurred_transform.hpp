// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics2d/vector2f.hpp>

namespace frantic {
namespace graphics {

template <typename FloatType>
class motion_blurred_transform {
    frantic::graphics::transform4t<FloatType> m_transform, m_inverse;
    std::vector<frantic::graphics::transform4t<FloatType>> m_transformArray, m_inverseArray;

    // This value represents the fraction of a 360 degree shutter that the motion path represents across its [0,1]
    // range. For instance, a 180 degree shutter would have a value of 0.5 here.
    FloatType m_motionBlurInterval;

  public:
    //////////////
    // Constructors
    //////////////

    motion_blurred_transform() {
        // Everything defaults to the identity transform
        m_motionBlurInterval = 0;
    }

    motion_blurred_transform( const frantic::graphics::transform4t<FloatType>& staticTransform ) {
        set( staticTransform );
    }

    motion_blurred_transform( const frantic::graphics::transform4t<FloatType>& startTransform,
                              const transform4t<FloatType>& endTransform ) {
        set( startTransform, endTransform );
    }

    motion_blurred_transform( const std::vector<frantic::graphics::transform4t<FloatType>>& animatedTransform,
                              FloatType motionBlurInterval = 0.5f ) {
        set( animatedTransform, motionBlurInterval );
    }

    motion_blurred_transform( const frantic::graphics::transform4t<FloatType>& staticTransform,
                              const std::vector<frantic::graphics::transform4t<FloatType>>& animatedTransform,
                              FloatType motionBlurInterval = 0.5f ) {
        set( staticTransform, animatedTransform, motionBlurInterval );
    }

#ifdef MAX_VERSION
    motion_blurred_transform( INode* inode, TimeValue t, float motionBlurInterval = 0.5f, float shutterBias = 0.f,
                              Interval* ivalid = 0 ) {
        m_transformArray.clear();

        Interval transformValid = FOREVER;
        m_transform = inode->GetNodeTM( t, &transformValid );

        if( ivalid != 0 )
            *ivalid &= transformValid;

        if( motionBlurInterval > 0 ) {
            Interval motionInterval( int( t - 0.5f * ( 1 - shutterBias ) * motionBlurInterval * GetTicksPerFrame() ),
                                     int( t + 0.5f * ( 1 + shutterBias ) * motionBlurInterval * GetTicksPerFrame() ) );

            if( transformValid.Start() > motionInterval.Start() || transformValid.End() < motionInterval.End() ) {
                m_transformArray.reserve( 10 );
                // There's motion blur, so we should get the animated transform.
                for( int motionSample = 0; motionSample < 10; ++motionSample ) {
                    float alpha = motionSample / 9.f;
                    TimeValue sampleTime =
                        TimeValue( motionInterval.Start() * ( 1 - alpha ) + motionInterval.End() * alpha );
                    m_transformArray.push_back( inode->GetNodeTM( sampleTime ) );
                }
            }
        }
        compute_inverses();
        m_motionBlurInterval = motionBlurInterval;
    }

    static motion_blurred_transform from_objtmafterwsm( INode* inode, TimeValue t, float motionBlurInterval = 0.5f,
                                                        float shutterBias = 0.f, Interval* ivalid = 0,
                                                        bool disableCameraMotionBlur = false ) {
        motion_blurred_transform result;

        Interval transformValid = FOREVER;
        // The main transform is always at the time
        result.m_transform = inode->GetObjTMAfterWSM( t, &transformValid );

        if( ivalid != 0 )
            *ivalid &= transformValid;

        if( motionBlurInterval > 0 ) {
            int mBlurStart = int( t - 0.5f * ( 1 - shutterBias ) * motionBlurInterval * GetTicksPerFrame() );
            int mBlurEnd = int( t + 0.5f * ( 1 + shutterBias ) * motionBlurInterval * GetTicksPerFrame() );

            if( disableCameraMotionBlur ) {
                int mBlurMiddle = mBlurStart + mBlurEnd / 2;
                mBlurStart = mBlurMiddle;
                mBlurEnd = mBlurMiddle;
            }

            Interval motionInterval( mBlurStart, mBlurEnd );

            // mprintf( "Time %d, interval %d %d\n", t, motionInterval.Start(), motionInterval.End() );

            if( transformValid.Start() > motionInterval.Start() || transformValid.End() < motionInterval.End() ) {
                if( !disableCameraMotionBlur ) {
                    result.m_transformArray.reserve( 10 );
                    // There's motion blur, so we should get the animated transform.
                    for( int motionSample = 0; motionSample < 10; ++motionSample ) {
                        float alpha = motionSample / 9.f;
                        TimeValue sampleTime =
                            TimeValue( motionInterval.Start() * ( 1 - alpha ) + motionInterval.End() * alpha );
                        result.m_transformArray.push_back( inode->GetObjTMAfterWSM( sampleTime ) );
                    }
                } else {
                    TimeValue sampleTime = TimeValue( motionInterval.Start() );
                    result.m_transformArray.push_back( inode->GetObjTMAfterWSM( sampleTime ) );
                }
            }
        }
        result.compute_inverses();
        result.m_motionBlurInterval = motionBlurInterval;

        return result;
    }
#endif

    //////////////
    // Setters
    //////////////

    void set( const frantic::graphics::transform4t<FloatType>& staticTransform ) {
        m_transform = staticTransform;
        m_transformArray.clear();
        compute_inverses();
        m_motionBlurInterval = 0;
    }

    void set( const frantic::graphics::transform4t<FloatType>& startTransform,
              const frantic::graphics::transform4t<FloatType>& endTransform ) {
        // Set the transform to the average
        m_transform = startTransform;
        m_transform += endTransform;
        m_transform *= 0.5f;

        // Set the array to hold the two values
        m_transformArray.clear();
        m_transformArray.push_back( startTransform );
        m_transformArray.push_back( endTransform );

        compute_inverses();
        m_motionBlurInterval = 0;
    }

    void set( const std::vector<frantic::graphics::transform4t<FloatType>>& animatedTransform,
              FloatType motionBlurInterval = 0.5f ) {
        m_transformArray = animatedTransform;
        if( m_transformArray.size() > 0 ) {
            if( m_transformArray.size() % 2 == 1 ) {
                m_transform = m_transformArray[m_transformArray.size() / 2];
            } else {
                m_transform = 0.5f * ( m_transformArray[m_transformArray.size() / 2 - 1] +
                                       m_transformArray[m_transformArray.size() / 2] );
            }
            compute_inverses();
        } else {
            m_transform.set_to_identity();
            m_inverse.set_to_identity();
            m_inverseArray.clear();
        }
        m_motionBlurInterval = motionBlurInterval;
    }

    void set( const frantic::graphics::transform4t<FloatType>& staticTransform,
              const std::vector<frantic::graphics::transform4t<FloatType>>& animatedTransform,
              FloatType motionBlurInterval = 0.5f ) {
        m_transform = staticTransform;
        m_transformArray = animatedTransform;
        compute_inverses();
        m_motionBlurInterval = motionBlurInterval;
    }

    void set_to_identity() {
        m_transform.set_to_identity();
        m_inverse.set_to_identity();
        m_transformArray.clear();
        m_inverseArray.clear();
        m_motionBlurInterval = 0;
    }

    //////////////
    // Operators
    //////////////

    void invert() {
        std::swap( m_transform, m_inverse );
        m_transformArray.swap( m_inverseArray );
    }

    void remove_scale() {
        m_transform = m_transform.to_scale_free();
        for( unsigned i = 0; i < m_transformArray.size(); ++i ) {
            m_transformArray[i] = m_transformArray[i].to_scale_free();
        }
        compute_inverses();
    }

    // TODO: The remove_scale function uses polar decomposition which discards skew in a matrix. This is not desireable
    //        certain situations, so this method will remove scale without affecting skew. If that explanation is wrong,
    //        note that this function allows us to match the 3dsMax camera.
    void remove_scale_not_skew() {
        m_transform = m_transform.to_scale_free_not_skew();
        for( unsigned i = 0; i < m_transformArray.size(); ++i )
            m_transformArray[i] = m_transformArray[i].to_scale_free_not_skew();
        compute_inverses();
    }

    //////////////
    // Transformation methods
    //////////////

    frantic::graphics::vector3t<FloatType> transform_point( const frantic::graphics::vector3t<FloatType>& p ) const {
        return m_transform * p;
    }

    frantic::graphics::vector3t<FloatType> transform_point( const frantic::graphics::vector3t<FloatType>& p,
                                                            FloatType motionSegmentTime ) const {
        if( m_transformArray.size() > 1 ) {
            // Scale the position into the animation range
            motionSegmentTime *= ( m_transformArray.size() - 1 );
            // Clip to the range
            if( motionSegmentTime <= 0 ) {
                return m_transformArray.front() * p;
            }
            if( motionSegmentTime >= m_transformArray.size() - 1 ) {
                return m_transformArray.back() * p;
            }
            // Linear interpolate the actual camera position
            int lowerIndex = (int)motionSegmentTime;
            FloatType alpha = motionSegmentTime - lowerIndex;
            return frantic::graphics::linear_interpolate( m_transformArray[lowerIndex] * p,
                                                          m_transformArray[lowerIndex + 1] * p, alpha );
        } else {
            return m_transform * p;
        }
    }

    frantic::graphics::vector3t<FloatType>
    inverse_transform_point( const frantic::graphics::vector3t<FloatType>& p ) const {
        return m_inverse * p;
    }

    frantic::graphics::vector3t<FloatType> inverse_transform_point( const frantic::graphics::vector3t<FloatType>& p,
                                                                    FloatType motionSegmentTime ) const {
        if( m_inverseArray.size() > 1 ) {
            // Scale the position into the animation range
            motionSegmentTime *= ( m_inverseArray.size() - 1 );
            // Clip to the range
            if( motionSegmentTime <= 0 ) {
                return m_inverseArray.front() * p;
            }
            if( motionSegmentTime >= m_inverseArray.size() ) {
                return m_inverseArray.back() * p;
            }
            // Linear interpolate the actual camera position
            int lowerIndex = (int)motionSegmentTime;
            FloatType alpha = motionSegmentTime - lowerIndex;
            return frantic::graphics::linear_interpolate( m_inverseArray[lowerIndex] * p,
                                                          m_inverseArray[lowerIndex + 1] * p, alpha );
        } else {
            return m_inverse * p;
        }
    }

    frantic::graphics::vector3t<FloatType> transform_normal( const frantic::graphics::vector3t<FloatType>& p ) const {
        return m_inverse.transpose_transform_no_translation( p );
    }

    //////////////
    // Functions to retrieve and compute information
    //////////////

    bool is_static() const { return m_transformArray.size() < 2; }

    bool is_rigid_transformation_of( const motion_blurred_transform& xform ) const {
        // If both are static, this is trivially true
        if( is_static() && xform.is_static() )
            return true;

        // Only do the comparison if both have the same number of samples
        if( m_transformArray.size() != xform.m_transformArray.size() )
            return false;

        frantic::graphics::transform4t<FloatType> rigidTransform = xform.get_transform() * get_inverse_transform();
        for( unsigned i = 0; i < m_transformArray.size(); ++i ) {
            frantic::graphics::transform4t<FloatType> testRigidTransform =
                xform.m_transformArray[i] * m_inverseArray[i];
            // std::cout << "comparing matrix " << rigidTransform << std::endl;
            // std::cout << "with matrix " << testRigidTransform << std::endl;
            //  NOTE: if something is going wrong, it's probably here...
            if( !testRigidTransform.equals_relative( rigidTransform, 0.00001f ) )
                return false;
        }
        return true;
    }

    FloatType get_motion_blur_interval() const { return m_motionBlurInterval; }

    frantic::graphics::vector3t<FloatType> get_translation() const { return m_transform.translation(); }

    frantic::graphics::vector3t<FloatType> get_translation( FloatType motionSegmentTime ) const {
        if( m_transformArray.size() > 1 ) {
            // Scale the position into the animation range
            motionSegmentTime *= ( m_transformArray.size() - 1 );
            // Clip to the range
            if( motionSegmentTime <= 0 ) {
                return m_transformArray.front().translation();
            }
            if( motionSegmentTime >= m_transformArray.size() - 1 ) {
                return m_transformArray.back().translation();
            }
            // Linear interpolate the actual camera position
            int lowerIndex = (int)motionSegmentTime;
            FloatType alpha = motionSegmentTime - lowerIndex;
            return linear_interpolate( m_transformArray[lowerIndex].translation(),
                                       m_transformArray[lowerIndex + 1].translation(), alpha );
        } else {
            return m_transform.translation();
        }
    }

    const frantic::graphics::transform4t<FloatType>& get_transform() const { return m_transform; }

    frantic::graphics::transform4t<FloatType> get_transform( FloatType motionSegmentTime ) const {
        if( m_transformArray.size() > 1 ) {
            // Scale the position into the animation range
            motionSegmentTime *= ( m_transformArray.size() - 1 );
            // Clip to the range
            if( motionSegmentTime <= 0 ) {
                return m_transformArray.front();
            }
            if( motionSegmentTime >= m_transformArray.size() - 1 ) {
                return m_transformArray.back();
            }
            // Linear interpolate the actual camera position
            int lowerIndex = (int)motionSegmentTime;
            FloatType alpha = motionSegmentTime - lowerIndex;
            return transform4t<FloatType>::linear_interpolate( m_transformArray[lowerIndex],
                                                               m_transformArray[lowerIndex + 1], alpha );
        } else {
            return m_transform;
        }
    }

    /**
     * Returns a specific column of the linearly interpolated transform at the specified time. Equivalent
     * to get_transform( motionSegmentTime ).get_column( column );
     *
     * @param motionSegmentTime The [0,1] "time" value for evaluating the animated transform at.
     * @param column Which column of the 4x4 transform matrix to return. Discards the 4th element of that column.
     * @return The specified column of the transformation interpolated at the specified time.
     *
     * @note
     * This is more or less only going to get used for camera::view_direction() since we want just the 3rd column
     * interpolated. Its worth having this special cased though. Ex. Sorting particles along the viewing direction with
     * jittered mblur evaluates this in a tight loop.
     */
    inline frantic::graphics::vector3t<FloatType> get_transform_column( FloatType motionSegmentTime,
                                                                        int column ) const {
        if( m_transformArray.size() > 1 ) {
            // Scale the position into the animation range
            motionSegmentTime *= ( m_transformArray.size() - 1 );

            // Clip to the range
            if( motionSegmentTime <= 0 )
                return m_transformArray.front().get_column( column );

            if( motionSegmentTime >= m_transformArray.size() - 1 )
                return m_transformArray.back().get_column( column );

            // Linear interpolate the actual camera position
            int lowerIndex = (int)motionSegmentTime;
            FloatType alpha = motionSegmentTime - lowerIndex;

            return frantic::math::lerp( m_transformArray[lowerIndex].get_column( column ),
                                        m_transformArray[lowerIndex + 1].get_column( column ), alpha );
        } else {
            return m_transform.get_column( column );
        }
    }

    const frantic::graphics::transform4t<FloatType>& get_inverse_transform() const { return m_inverse; }

    frantic::graphics::transform4t<FloatType> get_inverse_transform( FloatType motionSegmentTime ) const {
        if( m_inverseArray.size() > 1 ) {
            // Scale the position into the animation range
            motionSegmentTime *= ( m_inverseArray.size() - 1 );
            // Clip to the range
            if( motionSegmentTime <= 0 ) {
                return m_inverseArray.front();
            }
            if( motionSegmentTime >= m_inverseArray.size() - 1 ) {
                return m_inverseArray.back();
            }
            // Linear interpolate the actual camera position
            int lowerIndex = (int)motionSegmentTime;
            FloatType alpha = motionSegmentTime - lowerIndex;
            return transform4t<FloatType>::linear_interpolate( m_inverseArray[lowerIndex],
                                                               m_inverseArray[lowerIndex + 1], alpha );
        } else {
            return m_inverse;
        }
    }

    // This function computes how much motion a given pixel travels through.  It
    // is the screen-space distance travelled by a point at infinity in the direction of
    // the given pixel.
    static FloatType get_pixel_motion_distance( const motion_blurred_transform& cameraTransform,
                                                frantic::graphics2d::vector2f pixel,
                                                frantic::graphics2d::size2f imageSize, FloatType xFieldOfView ) {
        FloatType tanHalfFov = tan( xFieldOfView / 2 );
        if( cameraTransform.m_transformArray.size() > 1 ) {
            FloatType result = 0;
            // Convert the screen space coordinates into a direction
            vector3t<FloatType> direction =
                vector3t<FloatType>::from_perspective_projection( pixel, imageSize, tanHalfFov );
            for( unsigned i = 0; i < cameraTransform.m_transformArray.size() - 1; ++i ) {
                // Transform the camera space direction into world space at time i
                vector3t<FloatType> worldDirection =
                    cameraTransform.m_transformArray[i].transform_no_translation( direction );
                // Transform the world space direction back to camera space at time i+1
                vector3t<FloatType> directionMoved =
                    cameraTransform.m_inverseArray[i + 1].transform_no_translation( worldDirection );
                // Convert the moved direction into screen space coordinates
                bool isValid = false; // TODO: refactor this a bit, with the spherical projection, etc.
                frantic::graphics2d::vector2f pixelMoved =
                    directionMoved.to_perspective_projection( imageSize, tanHalfFov, isValid );
                // Add the distance moved by this transformation to the total
                result += frantic::graphics2d::vector2f::distance( pixel, pixelMoved );
                //				std::cerr << "DEBUG: segment " << i << " has distance " <<
                //vector2f::distance( p, pMoved ) << ", and it went from " << p << " to " << pMoved << std::endl;
            }
            //			std::cerr << "DEBUG: got pixel distance of " << result << std::endl;
            return result;
        } else {
            //			std::cerr << "DEBUG: No animation in camera" << std::endl;
            return 0;
        }
    }

    void write_xml( std::ostream& out, const std::string& prefix ) const {
        out << prefix << "<transform>" << m_transform << "</transform>\n";
        if( m_transformArray.size() > 0 ) {
            out << prefix << "<animatedTransform>\n";
            out << prefix << " <shutter>" << 360 * get_motion_blur_interval() << "</shutter>\n";
            for( unsigned i = 0; i < m_transformArray.size(); ++i ) {
                out << prefix << " <transform>" << m_transformArray[i] << "</transform>\n";
            }
            out << prefix << "</animatedTransform>\n";
        }
    }

  private:
    static bool is_transform_array_static( const std::vector<frantic::graphics::transform4t<FloatType>>& transforms ) {
        // If the transform array only has one element, don't use it at all
        if( transforms.size() <= 1 )
            return true;

        // Check whether all the transforms are equal.  Tried to use a relative error
        // bound but even with a small epsilon it allowed non-static cases to slip through
        // as static.
        for( unsigned i = 0; i < transforms.size() - 1; ++i ) {
            if( !transforms[i].equals( transforms.back() ) )
                return false;
        }
        return true;
    }

    void compute_inverses() {
        m_inverse = m_transform.to_inverse();
        // If the transform array is static, get rid of it
        if( is_transform_array_static( m_transformArray ) )
            m_transformArray.clear();

        m_inverseArray.resize( m_transformArray.size() );
        for( unsigned i = 0; i < m_transformArray.size(); ++i ) {
            m_inverseArray[i] = m_transformArray[i].to_inverse();
        }
    }
};

} // namespace graphics
} // namespace frantic
