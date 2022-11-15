// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <tbb/parallel_for.h>

#include <frantic/graphics/transform4f.hpp>
#include <frantic/math/eigen.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

template <typename FloatType>
class transform_impl {
    typedef FloatType float_type;
    typedef frantic::graphics::vector3t<FloatType> vector3f_type;
    typedef frantic::graphics::vector4t<FloatType> vector4f_type;
    typedef frantic::graphics::quat4t<FloatType> quat4f_type;
    typedef frantic::graphics::transform4t<FloatType> transform4f_type;

    frantic::channels::channel_cvt_accessor<vector3f_type> m_posAccessor;
    frantic::channels::channel_cvt_accessor<vector3f_type>
        m_velAccessor; // Special treatment to handle derivative effect

    std::vector<frantic::channels::channel_cvt_accessor<vector3f_type>> m_pointChannels;
    std::vector<frantic::channels::channel_cvt_accessor<vector3f_type>> m_vectorChannels;
    std::vector<frantic::channels::channel_cvt_accessor<vector3f_type>> m_normalChannels;
    std::vector<frantic::channels::channel_cvt_accessor<quat4f_type>> m_orientationChannels;
    std::vector<frantic::channels::channel_cvt_accessor<vector4f_type>> m_rotationChannels;
    std::vector<frantic::channels::channel_cvt_accessor<FloatType>> m_scalarChannels;

    transform4f_type m_transform, m_transformInverse, m_transformDerivative;
    quat4f_type m_rotationPart;
    float_type m_transformScaleMagnitude;

    std::size_t m_particleSize;
    char* m_particles;

  public:
    transform_impl( const transform4f_type& tm, const transform4f_type& tmDerivative,
                    const frantic::channels::channel_map& pcm,
                    const std::map<frantic::tstring, prt::channel_interpretation::option>& channelInterpretations )
        : m_transform( tm )
        , m_transformInverse( tm.to_inverse() )
        , m_transformDerivative( tmDerivative )
        , m_particles( NULL ) {
        vector3f_type translation;
        transform4f_type perspective, rotation, stretch;
        m_transform.decompose( perspective, translation, rotation, stretch );

        // TODO: Convert the rotation matrix into a quaternion.

        // Extract the scale from the stretch matrix (its eigenvalues).
        // Use maximum scale component as m_transformScaleMagnitude,
        // which is used to scale the particle radius.
        float_type stretchUpperTriangle[6] = { stretch.get( 0, 0 ), stretch.get( 1, 0 ), stretch.get( 2, 0 ),
                                               stretch.get( 1, 1 ), stretch.get( 2, 1 ), stretch.get( 2, 2 ) };
        vector3f_type eigenvalues;
        frantic::math::linearalgebra::get_eigenvalues_symmetric_3x3( stretchUpperTriangle, eigenvalues[0],
                                                                     eigenvalues[1], eigenvalues[2] );
        m_transformScaleMagnitude = eigenvalues.max_abs_component();

        set_channel_map( pcm, channelInterpretations );
    }

    void
    set_channel_map( const frantic::channels::channel_map& pcm,
                     const std::map<frantic::tstring, prt::channel_interpretation::option>& channelInterpretations ) {
        m_posAccessor.reset();
        m_velAccessor.reset();

        if( pcm.has_channel( _T("Position") ) )
            m_posAccessor = pcm.get_cvt_accessor<vector3f_type>( _T("Position") );
        if( pcm.has_channel( _T("Velocity") ) )
            m_velAccessor = pcm.get_cvt_accessor<vector3f_type>( _T("Velocity") );

        m_pointChannels.clear();
        m_vectorChannels.clear();
        m_normalChannels.clear();
        m_orientationChannels.clear();
        m_rotationChannels.clear();
        m_scalarChannels.clear();

        // If we haven't specified anything specifically, use some assumptions based on known channel names.
        if( channelInterpretations.empty() ) {
            if( pcm.has_channel( _T("Acceleration") ) )
                m_vectorChannels.push_back( pcm.get_cvt_accessor<vector3f_type>( _T("Acceleration") ) );
            if( pcm.has_channel( _T("Normal") ) )
                m_normalChannels.push_back( pcm.get_cvt_accessor<vector3f_type>( _T("Normal") ) );
            if( pcm.has_channel( _T("Tangent") ) )
                m_normalChannels.push_back( pcm.get_cvt_accessor<vector3f_type>( _T("Tangent") ) );
            if( pcm.has_channel( _T("Binormal") ) )
                m_normalChannels.push_back( pcm.get_cvt_accessor<vector3f_type>( _T("Binormal") ) );
            if( pcm.has_channel( _T("Orientation") ) )
                m_orientationChannels.push_back( pcm.get_cvt_accessor<quat4f_type>( _T("Orientation") ) );
            if( pcm.has_channel( _T("Spin") ) )
                m_rotationChannels.push_back( pcm.get_cvt_accessor<vector4f_type>( _T("Spin") ) );
            if( pcm.has_channel( _T("Radius") ) )
                m_scalarChannels.push_back( pcm.get_cvt_accessor<float_type>( _T("Radius") ) );
        } else {
            for( std::size_t i = 0, iEnd = pcm.channel_count(); i < iEnd; ++i ) {
                const frantic::channels::channel& ch = pcm[i];

                // We've manually dealt w/ Position & Velocity so don't do so here.
                if( ch.name() == _T("Position") || ch.name() == _T("Velocity") )
                    continue;

                std::map<frantic::tstring, prt::channel_interpretation::option>::const_iterator it =
                    channelInterpretations.find( ch.name() );
                if( it != channelInterpretations.end() && it->second ) {
                    switch( it->second ) {
                    case prt::channel_interpretation::point:
                        m_pointChannels.push_back( pcm.get_cvt_accessor<vector3f_type>( ch.name() ) );
                        break;
                    case prt::channel_interpretation::vector:
                        m_vectorChannels.push_back( pcm.get_cvt_accessor<vector3f_type>( ch.name() ) );
                        break;
                    case prt::channel_interpretation::normal:
                        m_normalChannels.push_back( pcm.get_cvt_accessor<vector3f_type>( ch.name() ) );
                        break;
                    case prt::channel_interpretation::orientation:
                        m_orientationChannels.push_back( pcm.get_cvt_accessor<quat4f_type>( ch.name() ) );
                        break;
                    case prt::channel_interpretation::rotation:
                        m_rotationChannels.push_back( pcm.get_cvt_accessor<vector4f_type>( ch.name() ) );
                        break;
                    case prt::channel_interpretation::scalar:
                        m_scalarChannels.push_back( pcm.get_cvt_accessor<float_type>( ch.name() ) );
                        break;
                    default:
                        break;
                    }
                }
            } // For each channel
        }

        m_particleSize = pcm.structure_size();
    }

    // TODO: Evaluate a better way to pass the range around. Probably a custom Range object.
    void set_buffer( char* buffer ) { m_particles = buffer; }

    void operator()( char* p ) const {
        if( m_posAccessor.is_valid() ) {
            if( m_velAccessor.is_valid() ) // By chain rule of derivatives
                m_velAccessor.set( p, m_transform.transform_no_translation( m_velAccessor.get( p ) ) +
                                          m_transformDerivative * m_posAccessor.get( p ) );
            m_posAccessor.set( p, m_transform * m_posAccessor.get( p ) );
        } else if( m_velAccessor.is_valid() )
            m_velAccessor.set( p, m_transform.transform_no_translation( m_velAccessor.get( p ) ) );

        for( typename std::vector<frantic::channels::channel_cvt_accessor<vector3f_type>>::const_iterator
                 it = m_pointChannels.begin(),
                 itEnd = m_pointChannels.end();
             it != itEnd; ++it )
            it->set( p, m_transform * it->get( p ) );

        for( typename std::vector<frantic::channels::channel_cvt_accessor<vector3f_type>>::const_iterator
                 it = m_vectorChannels.begin(),
                 itEnd = m_vectorChannels.end();
             it != itEnd; ++it )
            it->set( p, m_transform.transform_no_translation( it->get( p ) ) );

        for( typename std::vector<frantic::channels::channel_cvt_accessor<vector3f_type>>::const_iterator
                 it = m_normalChannels.begin(),
                 itEnd = m_normalChannels.end();
             it != itEnd; ++it )
            it->set( p, m_transformInverse.transpose_transform_no_translation( it->get( p ) ) );

        for( typename std::vector<frantic::channels::channel_cvt_accessor<quat4f_type>>::const_iterator
                 it = m_orientationChannels.begin(),
                 itEnd = m_orientationChannels.end();
             it != itEnd; ++it )
            it->set( p, m_rotationPart * it->get( p ) );

        for( typename std::vector<frantic::channels::channel_cvt_accessor<vector4f_type>>::const_iterator
                 it = m_rotationChannels.begin(),
                 itEnd = m_rotationChannels.end();
             it != itEnd; ++it ) {
            vector4f_type angAxis = it->get( p );
            vector3f_type axis =
                frantic::graphics::rotate_point( m_rotationPart, vector3f_type( angAxis.x, angAxis.y, angAxis.z ) );

            it->set( p, vector4f_type( axis.x, axis.y, axis.z, angAxis.w ) );
        }

        for( typename std::vector<frantic::channels::channel_cvt_accessor<float_type>>::const_iterator
                 it = m_scalarChannels.begin(),
                 itEnd = m_scalarChannels.end();
             it != itEnd; ++it )
            it->set( p, m_transformScaleMagnitude * it->get( p ) );
    }

    void operator()( const tbb::blocked_range<std::size_t>& range ) const {
        char* p = m_particles + m_particleSize * range.begin();
        char* pEnd = m_particles + m_particleSize * range.end();
        for( ; p != pEnd; p += m_particleSize )
            ( *this )( p );
    }
};

// Transforms a particle stream using the given matrices.
//   transform: The position transform, to transform the positions
//   transformDerivative: The time derivative of the position transform, to add the transform motion into the
//   velocities.
template <typename FloatType>
class transformed_particle_istream : public delegated_particle_istream {
    transformed_particle_istream& operator=( const transformed_particle_istream& /*rhs*/ ); // not implemented
    // transformed_particle_istream(const transformed_particle_istream& rhs); // not implemented

    std::map<frantic::tstring, prt::channel_interpretation::option> m_transformTypes;

    transform_impl<FloatType> m_impl; // Must be listed after m_transformTypes.

  public:
    typedef std::map<frantic::tstring, prt::channel_interpretation::option> transform_types_t;

    transformed_particle_istream(
        boost::shared_ptr<particle_istream> stream, const frantic::graphics::transform4t<FloatType>& tm,
        const frantic::graphics::transform4t<FloatType>& tmDeriv = frantic::graphics::transform4t<FloatType>::zero(),
        const std::map<frantic::tstring, prt::channel_interpretation::option>& transformTypes = transform_types_t() )
        : delegated_particle_istream( stream )
        , m_transformTypes( transformTypes )
        , m_impl( tm, tmDeriv, stream->get_channel_map(), m_transformTypes ) {}

    virtual ~transformed_particle_istream() {}

    virtual void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_delegate->set_channel_map( pcm );
        m_impl.set_channel_map( pcm, m_transformTypes );
    }

    virtual bool get_particle( char* particleBuffer ) {
        if( !m_delegate->get_particle( particleBuffer ) )
            return false;

        m_impl( particleBuffer );
        return true;
    }

    virtual bool get_particles( char* particleBuffer, std::size_t& bufferSize ) {
        bool result = m_delegate->get_particles( particleBuffer, bufferSize );
        if( bufferSize > 0 ) {
            m_impl.set_buffer( particleBuffer );

#ifndef FRANTIC_DISABLE_THREADS
            tbb::parallel_for( tbb::blocked_range<std::size_t>( 0, bufferSize, 1000 ), m_impl );
#else
#pragma message( "Threads are disabled" )
            char* p = particleBuffer;
            std::size_t pSize = m_delegate->particle_size();
            for( std::size_t i = 0; i < bufferSize; ++i, p += pSize )
                m_impl( p );
#endif
        }

        return result;
    }
};

// This helper function
template <typename FloatType>
inline boost::shared_ptr<particle_istream>
apply_transform_to_particle_istream( const boost::shared_ptr<particle_istream>& pin,
                                     const frantic::graphics::transform4t<FloatType>& tm,
                                     const frantic::graphics::transform4t<FloatType>& tmDeriv ) {
    if( !( tm.is_identity() && tmDeriv.is_zero() ) )
        return boost::shared_ptr<particle_istream>( new transformed_particle_istream<FloatType>( pin, tm, tmDeriv ) );
    return pin;
}

template <typename FloatType>
inline boost::shared_ptr<particle_istream> apply_transform_to_particle_istream(
    const boost::shared_ptr<particle_istream>& pin, const frantic::graphics::transform4t<FloatType>& tm,
    const frantic::graphics::transform4t<FloatType>& tmDeriv,
    const std::map<frantic::tstring, prt::channel_interpretation::option>& transformTypes ) {
    if( !( tm.is_identity() && tmDeriv.is_zero() ) )
        return boost::shared_ptr<particle_istream>(
            new transformed_particle_istream<FloatType>( pin, tm, tmDeriv, transformTypes ) );
    return pin;
}

} // namespace streams
} // namespace particles
} // namespace frantic
