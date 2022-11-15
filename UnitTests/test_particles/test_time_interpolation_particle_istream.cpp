// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/cstdint.hpp>
#include <boost/scoped_array.hpp>

#include <frantic/math/utils.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/particle_array_particle_istream.hpp>
#include <frantic/particles/streams/time_interpolation_particle_istream.hpp>

using frantic::graphics::vector3f;

TEST( TimeInterpolationParticleIstream, NoIDChannel ) {
    frantic::channels::channel_map pcm;
    pcm.define_channel<vector3f>( _T( "Position" ) );
    pcm.define_channel<vector3f>( _T( "Velocity" ) );
    pcm.define_channel<float>( _T( "Age" ) );
    pcm.end_channel_definition();

    frantic::particles::particle_array leftParticles( pcm );
    frantic::particles::particle_array rightParticles( pcm );

    boost::scoped_array<char> particleBuffer( new char[pcm.structure_size()] );

    const frantic::channels::channel_cvt_accessor<vector3f> posAccessor =
        pcm.get_cvt_accessor<vector3f>( _T( "Position" ) );
    const frantic::channels::channel_cvt_accessor<vector3f> velAccessor =
        pcm.get_cvt_accessor<vector3f>( _T( "Velocity" ) );
    const frantic::channels::channel_cvt_accessor<float> ageAccessor = pcm.get_cvt_accessor<float>( _T( "Age" ) );

    posAccessor.set( particleBuffer.get(), vector3f( 0.0f ) );
    velAccessor.set( particleBuffer.get(), vector3f( 1.0f, 0.0f, 0.0f ) );
    ageAccessor.set( particleBuffer.get(), 0.0f );
    leftParticles.push_back( particleBuffer.get() );

    posAccessor.set( particleBuffer.get(), vector3f( -10.0f, 0.0f, 0.0f ) );
    velAccessor.set( particleBuffer.get(), vector3f( -10.0f, 0.0f, 0.0f ) );
    ageAccessor.set( particleBuffer.get(), 10.0f );
    rightParticles.push_back( particleBuffer.get() );

    // If there's no "ID" channel we just want the position to be extrapolated using velocity.
    // The right particles should be ignored.

    frantic::particles::streams::particle_istream_ptr leftStream =
        boost::make_shared<frantic::particles::streams::particle_array_particle_istream>( leftParticles );
    frantic::particles::streams::particle_istream_ptr rightStream =
        boost::make_shared<frantic::particles::streams::particle_array_particle_istream>( rightParticles );

    frantic::particles::streams::time_interpolation_particle_istream interpolationStream( leftStream, rightStream, 1.0f,
                                                                                          0.5f );

    interpolationStream.get_particle( particleBuffer.get() );

    const vector3f position = posAccessor.get( particleBuffer.get() );
    const vector3f velocity = velAccessor.get( particleBuffer.get() );
    const float age = ageAccessor.get( particleBuffer.get() );

    const vector3f expectedPosition( 0.5f, 0.0f, 0.0f );
    const vector3f expectedVelocity( 1.0f, 0.0f, 0.0f );

    EXPECT_EQ( expectedPosition, position );
    EXPECT_EQ( expectedVelocity, velocity );

    EXPECT_EQ( 0.5f, age );
}

TEST( TimeInterpolationParticleIstream, NoMatchingParticle ) {
    frantic::channels::channel_map pcm;
    pcm.define_channel<vector3f>( _T( "Position" ) );
    pcm.define_channel<vector3f>( _T( "Velocity" ) );
    pcm.define_channel<float>( _T( "Age" ) );
    pcm.define_channel<boost::int64_t>( _T( "ID" ) );
    pcm.end_channel_definition();

    frantic::particles::particle_array leftParticles( pcm );
    frantic::particles::particle_array rightParticles( pcm );

    boost::scoped_array<char> particleBuffer( new char[pcm.structure_size()] );

    const frantic::channels::channel_cvt_accessor<vector3f> posAccessor =
        pcm.get_cvt_accessor<vector3f>( _T( "Position" ) );
    const frantic::channels::channel_cvt_accessor<vector3f> velAccessor =
        pcm.get_cvt_accessor<vector3f>( _T( "Velocity" ) );
    const frantic::channels::channel_cvt_accessor<float> ageAccessor = pcm.get_cvt_accessor<float>( _T( "Age" ) );
    const frantic::channels::channel_cvt_accessor<boost::int64_t> idAccessor =
        pcm.get_cvt_accessor<boost::int64_t>( _T( "ID" ) );

    posAccessor.set( particleBuffer.get(), vector3f( 512.0f ) );
    velAccessor.set( particleBuffer.get(), vector3f( 512.0f ) );
    ageAccessor.set( particleBuffer.get(), 256.0f );
    idAccessor.set( particleBuffer.get(), 0 );
    leftParticles.push_back( particleBuffer.get() );

    posAccessor.set( particleBuffer.get(), vector3f( 0.0f ) );
    velAccessor.set( particleBuffer.get(), vector3f( 1.0f, 0.0f, 0.0f ) );
    ageAccessor.set( particleBuffer.get(), 0.0f );
    idAccessor.set( particleBuffer.get(), 1 );
    leftParticles.push_back( particleBuffer.get() );

    posAccessor.set( particleBuffer.get(), vector3f( 512.0f ) );
    velAccessor.set( particleBuffer.get(), vector3f( 512.0f ) );
    ageAccessor.set( particleBuffer.get(), 257.0f );
    idAccessor.set( particleBuffer.get(), 0 );
    rightParticles.push_back( particleBuffer.get() );

    posAccessor.set( particleBuffer.get(), vector3f( -10.0f, 0.0f, 0.0f ) );
    velAccessor.set( particleBuffer.get(), vector3f( -10.0f, 0.0f, 0.0f ) );
    ageAccessor.set( particleBuffer.get(), 10.0f );
    idAccessor.set( particleBuffer.get(), 2 );
    rightParticles.push_back( particleBuffer.get() );

    // If there is no matching particle by ID we just want the position to be extrapolated using velocity.
    // The right particle should be ignored.
    // The "512" particle only exist to make sure the old behaviour of defaulting to the first particle no longer
    // occurs. We don't care what the stream outputs for it.

    frantic::particles::streams::particle_istream_ptr leftStream =
        boost::make_shared<frantic::particles::streams::particle_array_particle_istream>( leftParticles );
    frantic::particles::streams::particle_istream_ptr rightStream =
        boost::make_shared<frantic::particles::streams::particle_array_particle_istream>( rightParticles );

    frantic::particles::streams::time_interpolation_particle_istream interpolationStream( leftStream, rightStream, 1.0f,
                                                                                          0.5f );

    interpolationStream.get_particle( particleBuffer.get() );
    interpolationStream.get_particle( particleBuffer.get() );

    const vector3f position = posAccessor.get( particleBuffer.get() );
    const vector3f velocity = velAccessor.get( particleBuffer.get() );
    const float age = ageAccessor.get( particleBuffer.get() );

    const vector3f expectedPosition( 0.5f, 0.0f, 0.0f );
    const vector3f expectedVelocity( 1.0f, 0.0f, 0.0f );

    EXPECT_EQ( expectedPosition, position );
    EXPECT_EQ( expectedVelocity, velocity );

    EXPECT_EQ( 0.5f, age );
}

TEST( TimeInterpolationParticleIStream, MatchingParticle ) {
    frantic::channels::channel_map pcm;
    pcm.define_channel<vector3f>( _T( "Position" ) );
    pcm.define_channel<float>( _T( "Age" ) );
    pcm.define_channel<boost::int64_t>( _T( "ID" ) );
    pcm.define_channel<vector3f>( _T( "Normal" ) );
    pcm.define_channel<vector3f>( _T( "Tangent" ) );
    pcm.end_channel_definition();

    frantic::particles::particle_array leftParticles( pcm );
    frantic::particles::particle_array rightParticles( pcm );

    boost::scoped_array<char> particleBuffer( new char[pcm.structure_size()] );

    const frantic::channels::channel_cvt_accessor<vector3f> posAccessor =
        pcm.get_cvt_accessor<vector3f>( _T( "Position" ) );
    const frantic::channels::channel_cvt_accessor<float> ageAccessor = pcm.get_cvt_accessor<float>( _T( "Age" ) );
    const frantic::channels::channel_cvt_accessor<boost::int64_t> idAccessor =
        pcm.get_cvt_accessor<boost::int64_t>( _T( "ID" ) );
    const frantic::channels::channel_cvt_accessor<vector3f> normalAccessor =
        pcm.get_cvt_accessor<vector3f>( _T( "Normal" ) );
    const frantic::channels::channel_cvt_accessor<vector3f> tangentAccessor =
        pcm.get_cvt_accessor<vector3f>( _T( "Tangent" ) );

    const vector3f leftPosition( 0.0f );
    const vector3f rightPosition( -10.0f, 0.0f, 0.0f );

    const vector3f leftNormal( 1.f, 1.f, 0.f );
    const vector3f rightNormal( 1.f, 0.f, 0.f );

    posAccessor.set( particleBuffer.get(), vector3f( 512.0f ) );
    ageAccessor.set( particleBuffer.get(), 256.0f );
    idAccessor.set( particleBuffer.get(), 0 );
    normalAccessor.set( particleBuffer.get(), leftNormal );
    tangentAccessor.set( particleBuffer.get(), leftNormal );
    leftParticles.push_back( particleBuffer.get() );

    posAccessor.set( particleBuffer.get(), leftPosition );
    ageAccessor.set( particleBuffer.get(), 0.0f );
    idAccessor.set( particleBuffer.get(), 1 );
    normalAccessor.set( particleBuffer.get(), leftNormal );
    tangentAccessor.set( particleBuffer.get(), leftNormal );
    leftParticles.push_back( particleBuffer.get() );

    posAccessor.set( particleBuffer.get(), vector3f( 512.0f ) );
    ageAccessor.set( particleBuffer.get(), 257.0f );
    idAccessor.set( particleBuffer.get(), 0 );
    normalAccessor.set( particleBuffer.get(), leftNormal );
    tangentAccessor.set( particleBuffer.get(), leftNormal );
    rightParticles.push_back( particleBuffer.get() );

    posAccessor.set( particleBuffer.get(), rightPosition );
    ageAccessor.set( particleBuffer.get(), 1.0f );
    idAccessor.set( particleBuffer.get(), 1 );
    normalAccessor.set( particleBuffer.get(), rightNormal );
    tangentAccessor.set( particleBuffer.get(), rightNormal );
    rightParticles.push_back( particleBuffer.get() );

    // If there is no matching particle by ID we just want the position to be extrapolated using velocity.
    // The right particle should be ignored.
    // The "512" particle only exist to make sure the old behaviour of defaulting to the first particle no longer
    // occurs. We don't care what the stream outputs for it.

    frantic::particles::streams::particle_istream_ptr leftStream =
        boost::make_shared<frantic::particles::streams::particle_array_particle_istream>( leftParticles );
    frantic::particles::streams::particle_istream_ptr rightStream =
        boost::make_shared<frantic::particles::streams::particle_array_particle_istream>( rightParticles );

    frantic::particles::streams::time_interpolation_particle_istream interpolationStream( leftStream, rightStream, 1.0f,
                                                                                          0.5f );

    interpolationStream.get_particle( particleBuffer.get() );
    interpolationStream.get_particle( particleBuffer.get() );

    const vector3f position = posAccessor.get( particleBuffer.get() );
    const float age = ageAccessor.get( particleBuffer.get() );

    const vector3f expectedPosition = frantic::math::lerp( leftPosition, rightPosition, 0.5f );

    frantic::graphics::quat4f normalQuatLeft( 0, leftNormal );
    frantic::graphics::quat4f normalQuatRight( 0, rightNormal );
    normalQuatLeft.normalize();
    normalQuatRight.normalize();
    const vector3f expectedNormal =
        frantic::graphics::quat4f::slerp( normalQuatLeft, normalQuatRight, 0.5f ).vector_part();

    const vector3f normal = normalAccessor.get( particleBuffer.get() );
    const vector3f tangent = tangentAccessor.get( particleBuffer.get() );

    EXPECT_EQ( normal, expectedNormal );
    EXPECT_EQ( tangent, expectedNormal );
    EXPECT_EQ( expectedPosition, position );
    EXPECT_EQ( 0.5f, age );
}
