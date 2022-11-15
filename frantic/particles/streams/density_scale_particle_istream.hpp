// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/particles/streams/particle_istream.hpp>

#ifndef FRANTIC_DISABLE_THREADS
#include <tbb/parallel_for.h>
#endif

namespace frantic {
namespace particles {
namespace streams {

class density_scale_impl {
    frantic::channels::channel_cvt_accessor<float> m_densityAccessor;
    float m_densityScale;

  public:
    density_scale_impl( float densityScale, const frantic::channels::channel_cvt_accessor<float>& densityChannel )
        : m_densityScale( densityScale )
        , m_densityAccessor( densityChannel ) {}

    void operator()( const particle_range& range ) const {
        for( char *p = range.begin_ptr(), *pEnd = range.end_ptr(); p != pEnd; p += range.structure_size() )
            m_densityAccessor.set( p, m_densityScale * m_densityAccessor.get( p ) );
    }
};

class density_scale_particle_istream : public delegated_particle_istream {
    frantic::channels::channel_cvt_accessor<float> m_densityAccessor;
    float m_densityScale;

  public:
    density_scale_particle_istream( const boost::shared_ptr<particle_istream>& istream, float densityScale )
        : delegated_particle_istream( istream )
        , m_densityScale( densityScale ) {
        const frantic::channels::channel_map& pcm = delegated_particle_istream::get_channel_map();

        m_densityAccessor.reset( 1.f );
        if( pcm.has_channel( _T("Density") ) )
            m_densityAccessor = pcm.get_cvt_accessor<float>( _T("Density") );
    }

    virtual ~density_scale_particle_istream() {}

    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_densityAccessor.reset( 1.f );
        if( pcm.has_channel( _T("Density") ) )
            m_densityAccessor = pcm.get_cvt_accessor<float>( _T("Density") );
        delegated_particle_istream::set_channel_map( pcm );
    }

    bool get_particle( char* p ) {
        if( !m_delegate->get_particle( p ) )
            return false;

        m_densityAccessor.set( p, m_densityScale * m_densityAccessor.get( p ) );

        return true;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        bool notEos = m_delegate->get_particles( buffer, numParticles );

#ifndef FRANTIC_DISABLE_THREADS
        tbb::parallel_for( particle_range( 0, numParticles, buffer, m_delegate->particle_size(), 1000 ),
                           density_scale_impl( m_densityScale, m_densityAccessor ) );
#else
#pragma message( "Threads are disabled" )
        density_scale_impl f( m_densityScale, m_densityAccessor );
        f( particle_range( 0, numParticles, buffer, m_delegate->particle_size(), 1000 ) );
#endif

        return notEos;
    }
};

/*
//A particle stream that scales the density by a provided constant.
//
class density_scale_particle_istream : public delegated_particle_istream{
private:
  float m_densityScale;

  bool m_hasDensity;
  frantic::channels::channel_cvt_accessor<float> m_density;

  void prepare_accessors() {
    const channel_map& pcm = m_delegate->get_channel_map();

    m_hasDensity = pcm.has_channel( "Density" );
    if( m_hasDensity )
      m_density = pcm.get_cvt_accessor<float>("Density");
  }

public:
  density_scale_particle_istream( const boost::shared_ptr<particle_istream>& istream, float densityScale )
    :	delegated_particle_istream(istream), m_densityScale( densityScale )
  {
    prepare_accessors();
  }

  virtual ~density_scale_particle_istream() {}

  void set_channel_map( const channel_map& particleChannelMap ) {
    m_delegate->set_channel_map( particleChannelMap );
    prepare_accessors();
  }

  bool get_particle( char* rawParticleBuffer ) {
    bool result = m_delegate->get_particle(rawParticleBuffer);
    if( result ) {
      if( m_hasDensity )
        m_density.set( rawParticleBuffer, m_density.get(rawParticleBuffer) * m_densityScale );
    }
    return result;
  }
};
*/

} // namespace streams
} // namespace particles
} // namespace frantic
