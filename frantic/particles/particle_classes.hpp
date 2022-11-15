// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/color3h.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <frantic/channels/channel_map.hpp>

namespace frantic {
namespace particles {

// Align on 2 byte boundaries so half values are fully packed
#pragma pack( push, 2 )

// Basic particle class, matching the structure of the Flood v1 particle file format.
//
// Switched from using floats to halfs for the velocity and density values, to conserve memory.
// The particle was 12 + 12 + 4 = 28 bytes before, now it's 12 + 6 + 2 + 6 = 26 bytes.
class basic_particle {
    frantic::graphics::vector3f m_position;
    half m_velocity[3];
    half m_density;
    frantic::graphics::color3h m_color;

  public:
    basic_particle() { m_density = 1; }

    basic_particle( const frantic::graphics::vector3f& position, float density = 1 )
        : m_position( position )
        , m_density( density ) {}

    basic_particle( const frantic::graphics::vector3f& position, const frantic::graphics::vector3f& velocity,
                    float density = 1, frantic::graphics::color3f color = frantic::graphics::color3f( 1.f ) )
        : m_position( position )
        , m_density( density )
        , m_color( color ) {
        m_velocity[0] = velocity.x;
        m_velocity[1] = velocity.y;
        m_velocity[2] = velocity.z;
    }

    frantic::graphics::vector3f get_position() const { return m_position; }

    void set_position( const frantic::graphics::vector3f& position ) { m_position = position; }

    float get_density() const { return m_density; }

    void set_density( float density ) { m_density = density; }

    frantic::graphics::vector3f get_velocity() const {
        return frantic::graphics::vector3f( m_velocity[0], m_velocity[1], m_velocity[2] );
    }

    void set_velocity( const frantic::graphics::vector3f& velocity ) {
        m_velocity[0] = velocity.x;
        m_velocity[1] = velocity.y;
        m_velocity[2] = velocity.z;
    }

    void set_color( const frantic::graphics::color3f& color ) { m_color = color; }

    frantic::graphics::color3f get_color() const { return m_color; }

    operator frantic::graphics::vector3f() const { return m_position; }

    static void define_channel_map( frantic::channels::channel_map& pcm ) {
        pcm.define_channel( _T("Position"), 3, frantic::channels::data_type_float32 );
        pcm.define_channel( _T("Velocity"), 3, frantic::channels::data_type_float16 );
        pcm.define_channel( _T("Density"), 1, frantic::channels::data_type_float16 );
        pcm.define_channel( _T("Color"), 3, frantic::channels::data_type_float16 );
    }
};

// This is a very basic particle class, used in the KrakatoaPRTLoader as the cached particle type.
class viewport_particle {
    frantic::graphics::vector3f m_position;
    frantic::graphics::vector3f m_velocity;
    boost::uint32_t m_color;

  public:
    viewport_particle()
        : m_position( 0 )
        , m_color( 0 ) {}
    viewport_particle( const frantic::graphics::vector3f& pos )
        : m_position( pos ) {}
    viewport_particle( const frantic::graphics::vector3f& pos, const frantic::graphics::color3f& col )
        : m_position( pos )
        , m_color( col.to_RGBA() ) {}

    frantic::graphics::vector3f get_position() const { return m_position; }
    frantic::graphics::vector3f get_velocity() const { return m_velocity; }
    frantic::graphics::color3f get_color() const { return frantic::graphics::color3f::from_RGBA( m_color ); }

    void set_position( const frantic::graphics::vector3f& pos ) { m_position = pos; }
    void set_velocity( const frantic::graphics::vector3f& vel ) { m_velocity = vel; }
    void set_color( const frantic::graphics::color3f& col ) { m_color = col.to_RGBA(); }

    operator frantic::graphics::vector3f() const { return m_position; }

    static void define_channel_map( frantic::channels::channel_map& pcm ) {
        pcm.define_channel( _T("Position"), 3, frantic::channels::data_type_float32 );
        pcm.define_channel( _T("Velocity"), 3, frantic::channels::data_type_float32 );
        pcm.define_channel( _T("UINTColor"), 1, frantic::channels::data_type_uint32 );
    }
};

// This class is for rendering particles with the volumetric particle renderer.
// It has a position, and is convertible to a vector3f, which is the way the renderer gets the position.
// The lighting color is used to accumulate the lighting each particle receives, or just to communicate what
// color the particle should be rendered as if no lighting is being done.
//
// Switched from using floats to halfs for the lighting and density values, to conserve memory.
// The particle was 12 + 12 + 4 = 28 bytes before, now it's 12 + 6 + 2 = 20 bytes.
class renderable_particle {
    // Position of the particle
    frantic::graphics::vector3f m_position;
    // Value used to compute the lighting of the render
    frantic::graphics::color3h m_lighting;
    // Density of this particle.  This is used for adaptive seeding.
    half m_density;

  public:
    renderable_particle() { m_density = 1; }

    renderable_particle( const frantic::graphics::vector3f& position, float density = 1 )
        : m_position( position )
        , m_density( density ) {}

    // No color support
    void set_color( const frantic::graphics::color3f& ) {}

    frantic::graphics::color3f get_color() const { return frantic::graphics::color3f( 1 ); }

    void set_lighting( const frantic::graphics::color3f& lighting ) { m_lighting = lighting; }

    void add_lighting( const frantic::graphics::color3f& lighting ) { m_lighting += lighting; }

    frantic::graphics::color3f get_raw_lighting() const { return m_lighting; }

    frantic::graphics::color3f get_lighting() const { return m_lighting; }

    frantic::graphics::vector3f get_position() const { return m_position; }

    // Assuming a motion interval of [0,1], this returns the particle position at the specified time
    frantic::graphics::vector3f get_position_at_time( float /*time*/ ) const {
        // No velocities stored in this particle type
        return m_position;
    }

    void set_position( const frantic::graphics::vector3f& position ) { m_position = position; }

    float get_density() const { return m_density; }

    void set_density( float density ) { m_density = density; }

    // velocity is 0
    frantic::graphics::vector3f get_normalized_velocity() const { return frantic::graphics::vector3f( 0.f, 0.f, 0.f ); }

    // for interface compatibility, just discard the velocities
    void set_normalized_velocity( const frantic::graphics::vector3f& /*normalizedVelocity*/ ) {}

    // for interface compatibility, just discard the velocities
    void set_velocity( const frantic::graphics::vector3f& /*velocity*/ ) {}

    // for interface compatibility, just discard the velocities
    void set_velocity( const frantic::graphics::vector3f& /*velocity*/, int /*frameRate*/,
                       float /*shutterAngleInDegrees*/ ) {}

    operator frantic::graphics::vector3f() const { return m_position; }

    static void define_channel_map( frantic::channels::channel_map& pcm ) {
        pcm.define_channel( _T("Position"), 3, frantic::channels::data_type_float32 );
        pcm.define_channel( _T("Lighting"), 3, frantic::channels::data_type_float16 );
        pcm.define_channel( _T("Density"), 1, frantic::channels::data_type_float16 );
    }
};

// This particle class is the same as the renderable particle, but adds a motion vector.  The
// normalized velocity is scaled such that (position + (time - 0.5f) * normalizedVelocity) where
// time is in [0, 1] covers the entire centered motion interval.
//
// Switched from using floats to halfs for the velocity value, to conserve memory.
// The particle was 32 + 12 = 44 bytes before, now it's 20 + 6 = 26 bytes.
class motion_blurred_particle : public renderable_particle {
    // TODO: explore various ways to quantize this velocity.  Could 24 bits be enough for decent
    // motion blur, if we do appropriate randomization within the quantization interval?
    half m_normalizedVelocity[3];

  public:
    motion_blurred_particle() {
        m_normalizedVelocity[0] = 0;
        m_normalizedVelocity[1] = 0;
        m_normalizedVelocity[2] = 0;
    }

    motion_blurred_particle( const frantic::graphics::vector3f& position )
        : renderable_particle( position ) {
        m_normalizedVelocity[0] = 0;
        m_normalizedVelocity[1] = 0;
        m_normalizedVelocity[2] = 0;
    }

    motion_blurred_particle( const frantic::graphics::vector3f& position, float density )
        : renderable_particle( position, density ) {}

    frantic::graphics::vector3f get_normalized_velocity() const {
        return frantic::graphics::vector3f( m_normalizedVelocity[0], m_normalizedVelocity[1], m_normalizedVelocity[2] );
    }

    void set_normalized_velocity( const frantic::graphics::vector3f& normalizedVelocity ) {
        m_normalizedVelocity[0] = normalizedVelocity.x;
        m_normalizedVelocity[1] = normalizedVelocity.y;
        m_normalizedVelocity[2] = normalizedVelocity.z;
    }

    void set_velocity( const frantic::graphics::vector3f& velocity, int frameRate = 24,
                       float shutterAngleInDegrees = 180.f ) {
        m_normalizedVelocity[0] = velocity.x * ( shutterAngleInDegrees / ( 360.f * frameRate ) );
        m_normalizedVelocity[1] = velocity.y * ( shutterAngleInDegrees / ( 360.f * frameRate ) );
        m_normalizedVelocity[2] = velocity.z * ( shutterAngleInDegrees / ( 360.f * frameRate ) );
    }

    // Assuming a motion interval of [0,1], this returns the particle position at the specified time
    frantic::graphics::vector3f get_position_at_time( float time ) const {
        // No velocities stored in this particle type
        return get_position() + ( time - 0.5f ) * get_normalized_velocity();
    }

    static void define_channel_map( frantic::channels::channel_map& pcm ) {
        renderable_particle::define_channel_map( pcm );
        pcm.define_channel( _T("Velocity"), 3, frantic::channels::data_type_float16 );
    }
};

// Same as the motion blurred particle, but has a color.
//
// Switched from using floats to halfs for the color value, to conserve memory.
// The particle was 44 + 12 = 56 bytes before, now it's 26 + 6 = 32 bytes.
class colored_particle : public motion_blurred_particle {
    frantic::graphics::color3h m_color;

  public:
    colored_particle()
        : m_color( 1.f ) {}

    colored_particle( const frantic::graphics::vector3f& position )
        : motion_blurred_particle( position )
        , m_color( 1.f ) {}

    colored_particle( const frantic::graphics::vector3f& position, float density )
        : motion_blurred_particle( position, density )
        , m_color( 1.f ) {}

    void set_color( const frantic::graphics::color3f& color ) { m_color = color; }

    frantic::graphics::color3f get_color() const { return m_color; }

    // The color modulates the lighting
    frantic::graphics::color3f get_lighting() const { return m_color * get_raw_lighting(); }

    static void define_channel_map( frantic::channels::channel_map& pcm ) {
        motion_blurred_particle::define_channel_map( pcm );
        pcm.define_channel( _T("Color"), 3, frantic::channels::data_type_float16 );
    }
};

// This particle class is the same as the renderable particle, but adds a normal.  It will probably
// gain a scale at some point as well.  This is currently used for loading particle flow systems which contain
// orientations.
class particle_with_normal : public colored_particle {
    // vector3f m_orientation;
    half m_normal[3];

  public:
    particle_with_normal() {}

    particle_with_normal( const frantic::graphics::vector3f& position )
        : colored_particle( position ) {}

    particle_with_normal( const frantic::graphics::vector3f& position, float density )
        : colored_particle( position, density ) {}

    particle_with_normal( const frantic::graphics::vector3f& position, const frantic::graphics::vector3f& orientation )
        : colored_particle( position ) {
        set_normal( orientation );
    }

    frantic::graphics::vector3f get_normal() const {
        return frantic::graphics::vector3f( m_normal[0], m_normal[1], m_normal[2] );
    }

    void set_normal( const frantic::graphics::vector3f& normal ) {
        m_normal[0] = normal.x;
        m_normal[1] = normal.y;
        m_normal[2] = normal.z;
    }

    static void define_channel_map( frantic::channels::channel_map& pcm ) {
        colored_particle::define_channel_map( pcm );
        pcm.define_channel( _T("Normal"), 3, frantic::channels::data_type_float16 );
    }
};

#pragma pack( pop )

} // namespace particles
} // namespace frantic
