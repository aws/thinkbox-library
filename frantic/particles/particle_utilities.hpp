// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <zlib.h>

#include <boost/random.hpp>

#include <boost/unordered_map.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/logging/console_progress_logger.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/math/morton_code.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_classes.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/concatenated_particle_istream.hpp>
#include <frantic/particles/streams/fractional_particle_istream.hpp>
#include <frantic/particles/streams/transformed_particle_istream.hpp>
#include <frantic/sort/sort.hpp>
#include <frantic/tinyxml/frantic_tinyxml_graphics.hpp>

// TODO: the mesh loading and saving code in here needs to be refactored a lot.

namespace frantic {
namespace particles {

namespace detail {

// Helper function to allow reserve to be called on vector but ignored on deque
template <class ContainerType>
struct container_helper {
    static void reserve( ContainerType&, std::size_t ) {}
};
template <class T>
struct container_helper<std::vector<T>> {
    static void reserve( std::vector<T>& container, std::size_t size ) { container.reserve( size ); }
};

template <>
struct container_helper<particle_array> {
    static void reserve( particle_array& container, std::size_t size ) { container.reserve( size ); }
};

class float_array_predicate {
    const std::vector<float>* m_floats;

  public:
    // NOTE: the default constructor puts the class in an invalid state which will crash if the predicate is used.
    float_array_predicate() { m_floats = 0; }
    float_array_predicate( const std::vector<float>& floats )
        : m_floats( &floats ) {}

    float operator()( std::size_t i ) const { return ( *m_floats )[i]; }

    bool operator()( std::size_t i, std::size_t j ) const { return ( *m_floats )[i] < ( *m_floats )[j]; }
};
} // namespace detail
////////// Predicates to sort particles indirectly through containers ////////////

// A predicate function object to be used for generating a permutation that indicates the order
// of the particles in ascending value based on the supplied predicate.
template <class Container, class Predicate>
class container_predicate {
  private:
    const Container* m_data;
    Predicate m_predicate;

  public:
    container_predicate()
        : m_data( NULL )
        , m_predicate() {}
    container_predicate( const Container& data, const Predicate& predicate )
        : m_data( &data )
        , m_predicate( predicate ) {}

    bool operator()( std::size_t i, std::size_t j ) const { return m_predicate( m_data->at( i ), m_data->at( j ) ); }

    float operator()( std::size_t i ) const { return m_predicate( m_data->at( i ) ); }
};

template <class Container, class Predicate>
inline container_predicate<Container, Predicate> make_container_predicate( const Container& c, const Predicate& p ) {
    return container_predicate<Container, Predicate>( c, p );
}

////////// Predicates to sort particles directly ////////////

// A predicate function object to be used for sorting in ascending radius from m_point.
template <typename VectorType>
struct point_distance {
    VectorType m_point;
    frantic::channels::channel_accessor<VectorType> m_posAccess;

  public:
    point_distance() {}
    point_distance( const frantic::channels::channel_map& pcm, const VectorType& pt )
        : m_point( pt )
        , m_posAccess( pcm.get_accessor<VectorType>( _T("Position") ) ) {}
    point_distance( const frantic::channels::channel_accessor<VectorType>& posAcc, const VectorType& pt )
        : m_point( pt )
        , m_posAccess( posAcc ) {}

    /** Returns true if p0 is closer to m_point than p1 is */
    bool operator()( const char* p0, const char* p1 ) const {
        return VectorType::distance_squared( m_posAccess.get( p0 ), m_point ) <
               VectorType::distance_squared( m_posAccess.get( p1 ), m_point );
    }

    /** Returns true if p0 is closer to m_point than sqrt(distanceSquared1), enables using std::lower_bound and friends
     */
    bool operator()( const char* p0, typename VectorType::float_type distanceSquared1 ) const {
        return VectorType::distance_squared( m_posAccess.get( p0 ), m_point ) < distanceSquared1;
    }
    /**
     * Returns true if sqrt(distanceSquared0) is less than the distance from p1 to m_point, enables using
     * std::lower_bound and friends
     */
    bool operator()( typename VectorType::float_type distanceSquared0, const char* p1 ) const {
        return distanceSquared0 < VectorType::distance_squared( m_posAccess.get( p1 ), m_point );
    }

    // A pre-computation function for caching predicate results
    typename VectorType::float_type operator()( const char* p ) const {
        return VectorType::distance_squared( m_posAccess.get( p ), m_point );
    }
};

// A predicate function object to be used for sorting in ascending radius from m_point.
// This functor applies motion blur at the particle MBlurTime time.
struct point_distance_at_time {
  private:
    frantic::graphics::vector3f m_point;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccess;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_velAccess;
    frantic::channels::channel_cvt_accessor<float> m_timeAccess;

  public:
    point_distance_at_time() {}
    point_distance_at_time( const frantic::channels::channel_map& pcm, const frantic::graphics::vector3f& pt,
                            float defaultTime )
        : m_point( pt )
        , m_posAccess( pcm.get_accessor<frantic::graphics::vector3f>( _T("Position") ) )
        , m_velAccess( pcm.get_cvt_accessor<frantic::graphics::vector3f>( _T("Velocity") ) )
        , m_timeAccess( defaultTime ) {
        if( pcm.has_channel( _T("MBlurTime") ) )
            m_timeAccess = pcm.get_cvt_accessor<float>( _T("MBlurTime") );
    }

    // Returns true if p0 is closer to m_point than p1 is.
    bool operator()( const char* p0, const char* p1 ) const {
        return frantic::graphics::vector3f::distance_squared(
                   m_posAccess.get( p0 ) + ( m_timeAccess.get( p0 ) - 0.5f ) * m_velAccess.get( p0 ), m_point ) <
               frantic::graphics::vector3f::distance_squared(
                   m_posAccess.get( p1 ) + ( m_timeAccess.get( p1 ) - 0.5f ) * m_velAccess.get( p1 ), m_point );
    }

    // A pre-computation function for caching predicate results
    float operator()( const char* p ) const {
        return frantic::graphics::vector3f::distance_squared(
            m_posAccess.get( p ) + ( m_timeAccess.get( p ) - 0.5f ) * m_velAccess.get( p ), m_point );
    }
};

// A predicate function object to be used for sorting in ascending position along the direction vector m_direction.
struct directed_distance {
  private:
    frantic::graphics::vector3f m_direction;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccess;

  public:
    directed_distance() {}
    directed_distance( const frantic::channels::channel_map& pcm, const frantic::graphics::vector3f& dir )
        : m_direction( dir )
        , m_posAccess( pcm.get_accessor<frantic::graphics::vector3f>( _T("Position") ) ) {}

    // Returns true if p0 projects closer to negative infinity on the ray at the origin
    // point in m_direction than p1 does.
    bool operator()( const char* p0, const char* p1 ) const {
        return frantic::graphics::vector3f::dot( m_posAccess.get( p0 ), m_direction ) <
               frantic::graphics::vector3f::dot( m_posAccess.get( p1 ), m_direction );
    }

    // A pre-computation function for caching predicate results
    float operator()( const char* p ) const {
        return frantic::graphics::vector3f::dot( m_posAccess.get( p ), m_direction );
    }
};

// A predicate function object to be used for sorting in ascending position along the direction vector m_direction.
// This functor applies motion blur at the particle MBlurTime time.
struct directed_distance_at_time {
  private:
    frantic::graphics::vector3f m_direction;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccess;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_velAccess;
    frantic::channels::channel_cvt_accessor<float> m_timeAccess;
    float m_mblurTimeScale;

  public:
    directed_distance_at_time() {}
    directed_distance_at_time( const frantic::channels::channel_map& pcm, const frantic::graphics::vector3f& dir,
                               float defaultTime, float mblurTimeScale )
        : m_direction( dir )
        , m_posAccess( pcm.get_accessor<frantic::graphics::vector3f>( _T("Position") ) )
        , m_velAccess( frantic::graphics::vector3f( 0 ) )
        , // Made this not required.
        m_timeAccess( defaultTime )
        , m_mblurTimeScale( mblurTimeScale ) {
        if( pcm.has_channel( _T("Velocity") ) )
            m_velAccess = pcm.get_cvt_accessor<frantic::graphics::vector3f>( _T("Velocity") );
        if( pcm.has_channel( _T("MBlurTime") ) )
            m_timeAccess = pcm.get_cvt_accessor<float>( _T("MBlurTime") );
    }

    // Returns true if p0 projects closer to negative infinity on the ray at the origin
    // point in m_direction than p1 does.
    bool operator()( const char* p0, const char* p1 ) const {
        // m_timeAccess is in [0,1] motion blur interval. This needs to be scaled into real time
        // in seconds when using it with the velocity.
        float t0 = m_mblurTimeScale * ( m_timeAccess.get( p0 ) - 0.5f );
        float t1 = m_mblurTimeScale * ( m_timeAccess.get( p1 ) - 0.5f );
        return frantic::graphics::vector3f::dot( m_posAccess.get( p0 ) + t0 * m_velAccess.get( p0 ), m_direction ) <
               frantic::graphics::vector3f::dot( m_posAccess.get( p1 ) + t1 * m_velAccess.get( p1 ), m_direction );
    }

    // A pre-computation function for caching predicate results (also for radix sorts ...)
    float operator()( const char* p ) const {
        // m_timeAccess is in [0,1] motion blur interval. This needs to be scaled into real time
        // in seconds when using it with the velocity.
        float t = m_mblurTimeScale * ( m_timeAccess.get( p ) - 0.5f );
        return frantic::graphics::vector3f::dot( m_posAccess.get( p ) + t * m_velAccess.get( p ), m_direction );
    }
};

// A predicate function object to be used for sorting in ascending position along the direction vector m_direction.
// This functor applies motion blur at the particle MBlurTime time. This version uses frame relative times,
// instead of centerTime relative time.
struct directed_distance_at_abs_time {
  private:
    frantic::graphics::vector3f m_direction;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccess;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_velAccess;
    frantic::channels::channel_cvt_accessor<float> m_timeAccess;

  public:
    directed_distance_at_abs_time() {}
    directed_distance_at_abs_time( const frantic::channels::channel_map& pcm, const frantic::graphics::vector3f& dir,
                                   float defaultTime )
        : m_direction( dir )
        , m_posAccess( pcm.get_accessor<frantic::graphics::vector3f>( _T("Position") ) )
        , m_velAccess( frantic::graphics::vector3f( 0 ) )
        , // Made this not required.
        m_timeAccess( defaultTime ) {
        if( pcm.has_channel( _T("Velocity") ) )
            m_velAccess = pcm.get_cvt_accessor<frantic::graphics::vector3f>( _T("Velocity") );
        if( pcm.has_channel( _T("MBlurTime") ) )
            m_timeAccess = pcm.get_cvt_accessor<float>( _T("MBlurTime") );
    }

    // Returns true if p0 projects closer to negative infinity on the ray at the origin
    // point in m_direction than p1 does.
    bool operator()( const char* p0, const char* p1 ) const {
        return frantic::graphics::vector3f::dot( m_posAccess.get( p0 ) + m_timeAccess.get( p0 ) * m_velAccess.get( p0 ),
                                                 m_direction ) <
               frantic::graphics::vector3f::dot( m_posAccess.get( p1 ) + m_timeAccess.get( p1 ) * m_velAccess.get( p1 ),
                                                 m_direction );
    }

    // A pre-computation function for caching predicate results (also for radix sorts ...)
    float operator()( const char* p ) const {
        return frantic::graphics::vector3f::dot( m_posAccess.get( p ) + m_timeAccess.get( p ) * m_velAccess.get( p ),
                                                 m_direction );
    }
};

////////// Functions to load particles from different file formats ///////////

// Loads per-particle color information from an external file. An ugly methodology,
// but it is required (for now) since most particle formats do not use colors.
template <class InputParticleType>
void load_particle_color_file( const std::string& filename, std::vector<InputParticleType>& particlesLoaded ) {
    if( !frantic::files::file_exists( filename ) )
        throw std::runtime_error( "load_particle_color_file() : filename invalid." );

    gzFile fin = gzopen( filename.c_str(), "rb" );
    if( !fin )
        throw std::runtime_error( "load_particle_color_file() : failed to load file: " + filename );

    char headerBuffer[sizeof( char ) + sizeof( int )];
    int nRead = gzread( fin, headerBuffer, sizeof( headerBuffer ) / sizeof( char ) );

    if( headerBuffer[0] != char( 1 ) )
        throw std::runtime_error( "load_particle_color_file() : color file version unsupported." );

    int nParticles = *(int*)&headerBuffer[1];
    if( nParticles != (int)particlesLoaded.size() )
        throw std::runtime_error( "load_particle_color_file() : particle count mismatch." );

    for( typename std::vector<InputParticleType>::iterator it = particlesLoaded.begin(); it != particlesLoaded.end();
         ++it ) {
        float colorBuffer[3];

        nRead = gzread( fin, colorBuffer, sizeof( float ) * 3 );
        if( nRead != sizeof( float ) * 3 )
            throw std::runtime_error( "load_particle_color_file() : zlib compressed read failed." );

        ( *it ).set_color( frantic::graphics::color3f( colorBuffer[0], colorBuffer[1], colorBuffer[2] ) );
    }

    gzclose( fin );
}

/// TODO: Remove this function
template <class ParticleType>
void load_particleflow_file( const std::string& filename, std::vector<ParticleType>& particles, int frameRate,
                             float shutterAngleInDegrees ) {
    std::fstream file( filename.c_str(), std::ios::in | std::ios::binary );
    if( !file.is_open() )
        throw std::runtime_error( "Error: " + filename + " could not be opened for reading." );

    frantic::graphics::vector3f position, rotation, velocity;

    // Original pf version does not have an identifier so we need a string to represent the version
    char version[14];
    file.read( version, sizeof( char ) * 14 );

    if( strncmp( version, "PFlowExportOne", 14 ) == 0 ) {
        float particlePtr[9];

        // get all 9 floats for a particle
        while( file.read( (char*)&particlePtr, sizeof( float ) * 9 ) ) {
            position.x = particlePtr[0];
            position.y = particlePtr[1];
            position.z = particlePtr[2];
            rotation.x = particlePtr[3];
            rotation.y = particlePtr[4];
            rotation.z = particlePtr[5];
            velocity.x = particlePtr[6];
            velocity.y = particlePtr[7];
            velocity.z = particlePtr[8];

            // add this particle to our vector
            ParticleType p( position, rotation );
            p.set_velocity( velocity, frameRate, shutterAngleInDegrees );

            particles.push_back( p );
        }
    } else {
        float particlePtr[6];
        file.seekg( std::ios::beg );

        // get all 6 floats for a particle
        while( file.read( (char*)&particlePtr, sizeof( float ) * 6 ) ) {
            position.x = particlePtr[0];
            position.y = particlePtr[1];
            position.z = particlePtr[2];
            rotation.x = particlePtr[3];
            rotation.y = particlePtr[4];
            rotation.z = particlePtr[5];
            // add this particle to our vector
            particles.push_back( ParticleType( position, rotation ) );
        }
    }

    file.close();

} // load_particleflow_file

/// Transforms a set of particles by a given matrix
template <class ParticleType>
void transform_particles( std::vector<ParticleType>& particles, const frantic::graphics::transform4f& xform ) {
    for( unsigned int i = 0; i < particles.size(); ++i ) // transform each point
        particles[i].set_position( xform * particles[i] );
}

// This adds one particle system to another, applying a specified transform and color
// to it.
template <class InputParticleType, class OutputParticleType>
void add_transformed_particles( const std::vector<InputParticleType>& particlesNew,
                                const frantic::graphics::transform4f& particleXform,
                                const frantic::graphics::vector3f& velocity, const frantic::graphics::color3f& color,
                                std::vector<OutputParticleType>& particlesOut, bool bPerParticleColor = false ) {
    for( unsigned i = 0; i < particlesNew.size(); ++i ) {
        // Get the particle
        OutputParticleType particle = particlesNew[i];
        // Transform its position
        particle.set_position( particleXform * particle );
        // Set the velocity
        particle.set_normalized_velocity(
            velocity + particleXform.transform_no_translation( particle.get_normalized_velocity() ) );
        // Set the color
        if( !bPerParticleColor )
            particle.set_color( color );
        else {
            particle.set_color( particlesNew[i].get_color() );
        }

        // Add it to the output vector
        particlesOut.push_back( particle );
    }
}

template <class ParticleContainerType>
void load_particles_from_stream( streams::particle_istream& pin, ParticleContainerType& outParticles,
                                 frantic::logging::progress_logger& progress ) {
    typedef typename ParticleContainerType::value_type ParticleType;

    // Set the reading structure to match that of this particle
    frantic::channels::channel_map pcm;
    ParticleType::define_channel_map( pcm );
    pcm.end_channel_definition( 1, true );

    pin.set_channel_map( pcm );

    boost::int64_t particleCount;
    particleCount = pin.particle_count();

    if( particleCount > 0 ) {
        // detail::container_helper<ParticleContainerType>::reserve( outParticles, (std::size_t)particleCount +
        // outParticles.size() );
        std::size_t startIndex = 0;
        outParticles.resize( ( std::size_t )( particleCount + outParticles.size() ) );
        for( boost::int64_t i = 0; i < particleCount; ++i ) {
            if( !pin.get_particle( &outParticles[startIndex + (std::size_t)i] ) ) {
                throw std::runtime_error( "load_particles_from_stream: Failed to load particle number " +
                                          boost::lexical_cast<std::string>( i ) + " from particle stream \"" +
                                          pin.name() + "\"" );
            }

            progress.update_progress( i + 1, particleCount );
        }
    } else if( particleCount == -1 ) {
        ParticleType particle;
        while( pin.get_particle( &particle ) ) {
            outParticles.push_back( particle );
        }

        progress.update_progress( 1, 1 );
    }
}

inline boost::shared_ptr<streams::particle_istream> get_flood_reflow_istream( const std::string& particleFile ) {
    using namespace frantic::particles::streams;

    if( particleFile.find( '@' ) == std::string::npos ) {
        return boost::shared_ptr<particle_istream>(
            particle_file_istream_factory( frantic::strings::to_tstring( particleFile ) ) );
    } else {
        // Frame "test0035.prt@0.34" should be (1-0.34) * "test0035.prt" + 0.34 * "test0036.prt".
        std::vector<std::string> parts;
        frantic::strings::split( particleFile, parts, '@' );
        if( parts.size() != 2 )
            throw std::runtime_error( "load_flood_reflow_particles: Filename \"" + particleFile +
                                      "\" could not be parsed as a proper '@' frame-interpolated filename (like "
                                      "\"test0035.prt@0.34\" for example)" );

        std::string particleFileA = parts[0];
        std::string particleFileB = frantic::files::increment_sequence_number( particleFileA );

        float alpha = (float)atof( parts[1].c_str() );

        if( alpha < 0 || alpha > 1 )
            throw std::runtime_error( "get_flood_reflow_istream: Filename \"" + particleFile +
                                      "\" could not be parsed as a proper '@' frame-interpolated filename (like "
                                      "\"test0035.prt@0.34\" for example)" );

        if( frantic::logging::is_logging_stats() )
            std::cout << "Blending particles with fraction " << alpha << " between files \"" << particleFileA
                      << "\" and \"" << particleFileB << "\"" << std::endl;

        if( alpha == 0 )
            return boost::shared_ptr<particle_istream>(
                particle_file_istream_factory( frantic::strings::to_tstring( particleFileA ) ) );
        else if( alpha == 1 )
            return boost::shared_ptr<particle_istream>(
                particle_file_istream_factory( frantic::strings::to_tstring( particleFileB ) ) );
        else {
            /*boost::shared_ptr<particle_istream> finA( particle_file_istream_factory(particleFileA) );
            boost::shared_ptr<particle_istream> finB( particle_file_istream_factory(particleFileB) );
            return boost::shared_ptr<particle_istream>( new interpolating_particle_istream( finA, finB, alpha ) );*/
            throw std::runtime_error( "get_flood_reflow_istream() has been deprecated and destroyed." );
        }
    }
}

inline boost::shared_ptr<streams::particle_istream> get_flood_reflow_istream( tinyxml2::XMLHandle xml ) {
    using namespace frantic::particles::streams;
    using namespace frantic::graphics;

    if( xml.FirstChildElement( "reflowParticles" ).ToNode() == 0 ) {
        transform4f transform = frantic::tinyxml::get_transform4f( "transform", transform4f::identity(), xml, false );
        transform4f transformDerivative =
            tinyxml::get_transform4f( "transformDerivative", transform4f::zero(), xml, false );
        std::string fileName = tinyxml::get_xml_value<std::string>( "datafile", "", xml );
        if( fileName == "" ) {
            throw std::runtime_error( "get_flood_reflow_istream: No particle file was specified in xml node \"" +
                                      tinyxml::get_xml_handle_path( xml, "datafile" ) + "\"." );
        }

        boost::shared_ptr<particle_istream> result = get_flood_reflow_istream( fileName );
        if( transform.is_identity() && transformDerivative.is_zero() )
            return result;
        else
            return boost::shared_ptr<particle_istream>(
                new transformed_particle_istream<float>( result, transform, transformDerivative ) );
    } else {
        int count = 0;
        tinyxml2::XMLHandle reflowHandle = xml.FirstChildElement( "reflowParticles" );
        while( reflowHandle.ToNode() != 0 ) {
            ++count;
            reflowHandle = tinyxml2::XMLHandle( reflowHandle.ToNode()->NextSiblingElement( "reflowParticles" ) );
        }
        if( count == 1 ) {
            // If there's only one reflow particle reference then just return the one stream object
            return get_flood_reflow_istream( xml.FirstChildElement( "reflowParticles" ) );
        } else {
            // Otherwise create a concatenated stream
            std::vector<boost::shared_ptr<particle_istream>> streams;
            reflowHandle = xml.FirstChildElement( "reflowParticles" );
            while( reflowHandle.ToNode() != 0 ) {
                streams.push_back( get_flood_reflow_istream( reflowHandle ) );
                reflowHandle = tinyxml2::XMLHandle( reflowHandle.ToNode()->NextSiblingElement( "reflowParticles" ) );
            }
            return boost::shared_ptr<particle_istream>( new concatenated_particle_istream( streams ) );
        }
    }
}

// For testing, this creates a cube of constant density particles
template <class ParticleContainerType>
void load_test_particles( tinyxml2::XMLHandle xml, ParticleContainerType& outParticles, float particleFraction = 1.f ) {
    typedef typename ParticleContainerType::value_type ParticleType;

    frantic::graphics::transform4f particleXform =
        tinyxml::get_transform4f( "transform", frantic::graphics::transform4f(), xml, false );
    frantic::graphics::boundbox3f bounds = tinyxml::get_boundbox3f( "boundbox", xml );
    int particleCount = tinyxml::get_int( "particleCount", xml );
    boost::uint32_t seed = tinyxml::get_int( "randomSeed", 42, xml );
    frantic::graphics::color3f color = tinyxml::get_color3f( "color", xml );
    float densityPerVolume = tinyxml::get_float( "densityPerVolume", 1, xml );

    if( particleFraction != 1 )
        particleCount = (int)( particleCount * particleFraction );

    float densityPerParticle = densityPerVolume * bounds.volume() / particleCount;

    boost::mt19937 generator( seed );
    boost::uniform_real<float> uni_dist( 0, 1 );
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float>> rng( generator, uni_dist );

    for( int i = 0; i < particleCount; ++i ) {
        ParticleType particle( particleXform * bounds.random_vector( rng ), densityPerParticle );
        particle.set_color( color );
        outParticles.push_back( particle );
    }
}

// This parses the scene xml, and loads whatever type of particles are specified in the
// scene.  The particle fraction specifies what fraction of the total particles to actually load, and can be used
// for viewport previews or low quality test renders.
template <class ParticleContainerType>
void load_particles_from_scene( tinyxml2::XMLHandle sceneXml, ParticleContainerType& outParticles, float particleFraction,
                                boost::int64_t particleLimit, frantic::logging::progress_logger& progress ) {
    typedef typename ParticleContainerType::value_type ParticleType;

    /////// Load any particle reflow particle simulations ///////

    // First count how many particles there are
    boost::shared_ptr<streams::particle_istream> pin = get_flood_reflow_istream( sceneXml );
    // Load a fraction of the particles if requested
    if( particleFraction < 1 || particleLimit < pin->particle_count() ) {
        pin = boost::shared_ptr<streams::particle_istream>(
            new streams::fractional_particle_istream( pin, particleFraction, particleLimit ) );
    }

    load_particles_from_stream( *pin, outParticles, progress );
    tinyxml2::XMLHandle testHandle = sceneXml.FirstChildElement( "testParticles" );
    while( testHandle.ToNode() != nullptr ) {
        load_test_particles( testHandle, outParticles, particleFraction );

        testHandle = testHandle.ToNode()->NextSiblingElement( "testParticles" );
    }
}

/**
 * Takes a vector of splines and a time offset vector of splines and writes particles to a stream that represent the
 * splines. The stream sets the "SplineIndex", "KnotIndex" channels to identify the particle as part of a spline. The
 * output stream must have these channels.
 */
inline void
save_splines_to_stream( boost::shared_ptr<streams::particle_ostream> pout,
                        const std::vector<std::vector<frantic::graphics::vector3f>>& cachedSplines,
                        const std::vector<std::vector<frantic::graphics::vector3f>>& cachedSplinesOffsetTime,
                        float timeStep ) {
    using namespace frantic::channels;
    using namespace frantic::graphics;

    const channel_map& cm = pout->get_channel_map();

    std::vector<char> particle( cm.structure_size() );
    char* rawParticle = &particle[0];

    // stream MUST have these four channels.
    channel_accessor<vector3f> knotPosAcc = cm.get_accessor<vector3f>( _T("Position") );
    channel_cvt_accessor<vector3f> knotVelAcc = cm.get_cvt_accessor<vector3f>( _T("Velocity") );
    channel_accessor<boost::uint32_t> splineIndexAcc = cm.get_accessor<boost::uint32_t>( _T("SplineIndex") );
    channel_accessor<boost::uint32_t> knotIndexAcc = cm.get_accessor<boost::uint32_t>( _T("KnotIndex") );
    size_t numSplines = cachedSplines.size();
    for( size_t s = 0; s < numSplines; ++s ) {
        const std::vector<vector3f>& spline = cachedSplines[s];
        const std::vector<vector3f>& splineOffset = cachedSplinesOffsetTime[s];
        size_t numKnots = spline.size();
        if( numKnots != splineOffset.size() )
            throw std::runtime_error( "save_splines_to_stream error: Spline " + boost::lexical_cast<std::string>( s ) +
                                      " has different number of knots than its time offset spline." );
        for( size_t k = 0; k < numKnots; ++k ) {

            // set rawParticle with this particle's data
            const vector3f& pos = spline[k];
            knotPosAcc( rawParticle ) = pos;
            knotVelAcc.set( rawParticle, ( splineOffset[k] - pos ) * timeStep );
            splineIndexAcc( rawParticle ) = static_cast<boost::uint32_t>( s );
            knotIndexAcc( rawParticle ) = static_cast<boost::uint32_t>( k );

            // add this particle to the IO stream
            pout->put_particle( rawParticle );
        }
    }
}

namespace detail {
// struct for data from knots. add more members if needed.
// this struct used to be defined within load_spline_from_stream, but gcc did not like that.
struct knot_data {
    frantic::graphics::vector3f position;
    frantic::graphics::vector3f velocity;
    void set( const frantic::graphics::vector3f& p, const frantic::graphics::vector3f& v ) {
        position = p;
        velocity = v;
    }
};
} // namespace detail

inline void load_splines_from_stream( boost::shared_ptr<streams::particle_istream> pin,
                                      std::vector<std::vector<frantic::graphics::vector3f>>& outSplines,
                                      std::vector<std::vector<frantic::graphics::vector3f>>& outSplinesTimeOffset,
                                      float timeStep ) {

    using namespace frantic::graphics;
    using namespace frantic::channels;

    const channel_map& cm = pin->get_channel_map();
    if( !cm.has_channel( _T("SplineIndex") ) || !cm.has_channel( _T("KnotIndex") ) )
        throw std::runtime_error(
            "load_splines_from_stream error: Particles must have a \"SplineIndex\" and \"KnotIndex\" "
            "channel in order to be processed as splines." );

    // this temporarily holds the splines as they are being constructed. we do hash (unordered_map) because it performs
    // better than a tree (map) with our common case (inserting in order by spline index)
    boost::unordered_map<int, std::map<int, detail::knot_data>> inSplines;

    channel_accessor<vector3f> knotPosAcc = cm.get_accessor<vector3f>( _T("Position") );
    channel_cvt_accessor<vector3f> knotVelAcc = cm.get_cvt_accessor<vector3f>( _T("Velocity") );
    channel_accessor<boost::uint32_t> splineIndexAcc = cm.get_accessor<boost::uint32_t>( _T("SplineIndex") );
    channel_accessor<boost::uint32_t> knotIndexAcc = cm.get_accessor<boost::uint32_t>( _T("KnotIndex") );

    std::vector<char> buffer( cm.structure_size() );
    char* rawBuffer = &buffer[0];

    // populate inSplines
    while( pin->get_particle( rawBuffer ) ) {
        std::map<int, detail::knot_data>& splineData = inSplines[splineIndexAcc( rawBuffer )];
        splineData[knotIndexAcc( rawBuffer )].set( knotPosAcc( rawBuffer ), knotVelAcc( rawBuffer ) );
    }

    // allocate the right number of output splines.
    int numSplines = static_cast<int>( inSplines.size() );
    outSplines.resize( numSplines );
    outSplinesTimeOffset.resize( numSplines );

    // convert our collection format (inSplines) to the output format of vectors (outSplines and outSplinesTimeOffset)
    std::vector<std::vector<vector3f>>::iterator outSplineIter0 = outSplines.begin();
    std::vector<std::vector<vector3f>>::iterator outSplineIter1 = outSplinesTimeOffset.begin();
    boost::unordered_map<int, std::map<int, detail::knot_data>>::const_iterator inSplineIter,
        inSplineIterEnd = inSplines.end();
#pragma warning( push )
#pragma warning( disable : 4913 )
    for( inSplineIter = inSplines.begin(); inSplineIter != inSplineIterEnd;
         ++inSplineIter, ++outSplineIter0, ++outSplineIter1 ) {
#pragma warning( pop )
        const std::map<int, detail::knot_data>& inSpline = inSplineIter->second;
        int numKnots = static_cast<int>( inSpline.size() );
        outSplineIter0->reserve( numKnots );
        outSplineIter1->reserve( numKnots );
        std::map<int, detail::knot_data>::const_iterator inKnotIter, inKnotIterEnd = inSpline.end();
        for( inKnotIter = inSpline.begin(); inKnotIter != inKnotIterEnd; ++inKnotIter ) {
            const detail::knot_data& inKnot = inKnotIter->second;
            outSplineIter0->push_back( inKnot.position );
            outSplineIter1->push_back( inKnot.position + inKnot.velocity * timeStep );
        }
    }
}

////////// Functions to query information about particle sets ///////////

template <class ParticleType>
frantic::graphics::vector3f get_largest_normalized_velocity( const std::vector<ParticleType>& particles ) {
    frantic::graphics::vector3f result;
    float resultMagnitudeSquared = 0;
    for( unsigned i = 0; i < particles.size(); ++i ) {
        frantic::graphics::vector3f normalizedVelocity = particles[i].get_normalized_velocity();
        float magnitudeSquared = normalizedVelocity.get_magnitude_squared();
        if( magnitudeSquared > resultMagnitudeSquared ) {
            resultMagnitudeSquared = magnitudeSquared;
            result = normalizedVelocity;
        }
    }
    return result;
}

template <class ParticleType>
frantic::graphics::boundbox3f get_particles_bounding_box( const std::vector<ParticleType>& particles,
                                                          float time = 0.5f ) {
    frantic::graphics::boundbox3f result;
    for( unsigned i = 0; i < particles.size(); ++i ) {
        result += particles[i].get_position_at_time( time );
    }
    return result;
}

/////////// Function for sorting particles /////////////

namespace detail {
// used in sort_particles_by_voxel
class sorter_fcn {
  public:
    frantic::channels::channel_accessor<frantic::graphics::vector3f> posAcc;
    frantic::volumetrics::voxel_coord_system vcs;
    bool operator()( const char* first, const char* second ) const {
        frantic::graphics::vector3 firstVoxel = frantic::graphics::vector3::from_floor(
            vcs.get_voxel_coord( posAcc( first ) ) - frantic::graphics::vector3f( 0.5f ) );
        frantic::graphics::vector3 secondVoxel = frantic::graphics::vector3::from_floor(
            vcs.get_voxel_coord( posAcc( second ) ) - frantic::graphics::vector3f( 0.5f ) );
        return firstVoxel < secondVoxel;
    }
};

/*
 * Finds the bounds of the supplied inserted particles and some additional statistics
 * Variance and mean is calculated using Knuth's incremental algorithm at
 * http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
 */
class running_position_statistics {
  public:
    running_position_statistics()
        : m_count( 0 ) {}

    boost::uint64_t get_count() const { return m_count; }

    frantic::graphics::vector3f get_mean() const { return m_mean; }

    frantic::graphics::vector3f get_variance() const {
        if( m_count > 1 ) {
            return m_q / static_cast<float>( m_count - 1 );
        } else {
            return frantic::graphics::vector3f( 0 );
        }
    }

    frantic::graphics::boundbox3f get_bounds() const { return m_bounds; }

    void insert( const frantic::graphics::vector3f& position ) {
        const frantic::graphics::vector3f delta = position - m_mean;

        ++m_count;
        m_mean += delta / static_cast<float>( m_count );
        m_q += frantic::graphics::vector3f::component_multiply( delta, position - m_mean );
        m_bounds += position;
    }

  private:
    boost::uint64_t m_count;
    frantic::graphics::vector3f m_mean;
    frantic::graphics::vector3f m_q;
    frantic::graphics::boundbox3f m_bounds;
};

// used for sort_particles_by_morton_code
class morton_sorter {
    frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor;
    frantic::graphics::vector3f boundingBoxCorner;
    float voxelLength;

    morton_sorter& operator=( const morton_sorter& ); // Disable assignment.

  public:
    morton_sorter( frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
                   const frantic::graphics::vector3f& boundingBoxCorner, float voxelLength )
        : positionAccessor( positionAccessor )
        , boundingBoxCorner( boundingBoxCorner )
        , voxelLength( voxelLength ) {}

    // compare morton codes
    bool operator()( const char* first, const char* second ) const {
        boost::uint64_t code1 = frantic::math::morton_code( boundingBoxCorner, voxelLength, positionAccessor( first ) );
        boost::uint64_t code2 =
            frantic::math::morton_code( boundingBoxCorner, voxelLength, positionAccessor( second ) );

        return code1 < code2;
    }

    /*
     * Finds the bounds of the box of particles.
     */
    static frantic::graphics::boundbox3f
    compute_particle_bounds( const frantic::particles::particle_array& particles ) {
        frantic::particles::particle_array::const_iterator it = particles.begin();
        frantic::particles::particle_array::const_iterator end = particles.end();
        frantic::channels::channel_accessor<frantic::graphics::vector3f> posAccessor =
            particles.get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") );
        frantic::graphics::boundbox3f bounds;
        for( ; it != end; ++it )
            bounds += posAccessor( *it );
        return bounds;
    }

    static inline frantic::graphics::vector3f multiply( const frantic::graphics::vector3f& lfs,
                                                        const frantic::graphics::vector3f& rhs ) {
        return frantic::graphics::vector3f( lfs.x * rhs.x, lfs.y * rhs.y, lfs.z * rhs.z );
    }

    /*
     * Finds the bounds of the box of particles and some additional statistics
     * Variance and mean is calculated using Knuth's incremental algorithm at
     * http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
     */
    static frantic::graphics::boundbox3f
    compute_particle_bounds( boost::shared_ptr<frantic::particles::streams::particle_istream> particleStream,
                             frantic::graphics::vector3f& outMean, frantic::graphics::vector3f& outVariance,
                             std::size_t& outCount, frantic::logging::progress_logger& logger ) {
        frantic::channels::channel_accessor<frantic::graphics::vector3f> posAccessor =
            particleStream->get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") );

        running_position_statistics stats;

        const std::size_t UpdatePauseTime = 1;
        const std::size_t BlockSize = 500000;
        std::size_t updateCounter = 0;

        const std::size_t structureSize = particleStream->get_channel_map().structure_size();
        std::vector<char> buffer( BlockSize * structureSize );
        std::size_t toRead = BlockSize;
        bool ok = particleStream->get_particles( &buffer[0], toRead );
        while( toRead > 0 ) {
            for( std::size_t n = 0; n < toRead; ++n ) {
                const frantic::graphics::vector3f position = posAccessor( &buffer[n * structureSize] );
                stats.insert( position );
            }

            ++updateCounter;
            if( updateCounter >= UpdatePauseTime ) {
                logger.update_progress( particleStream->particle_progress_index(),
                                        particleStream->particle_progress_count() );
                updateCounter = 0;
            }

            if( !ok )
                break;
            ok = particleStream->get_particles( &buffer[0], toRead );
        }

        outVariance = stats.get_variance();
        outMean = stats.get_mean();
        outCount = stats.get_count();
        return stats.get_bounds();
    }
};

} // namespace detail

inline void sort_particles_by_voxel( particle_array& particles, frantic::volumetrics::voxel_coord_system vcs ) {
    // use special sorting function to sort
    detail::sorter_fcn sorterFcn;
    sorterFcn.posAcc = particles.get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") );
    sorterFcn.vcs = vcs;
    frantic::sort::threaded_sort( particles.begin(), particles.end(), sorterFcn, 4 ); // num threads?
}

/**
 * Convenience function designed for computing a good vcs for the sort_particles_by_morton_code function.
 * This function computes the boundbox of the particles, and a reasonable voxel length.
 */
inline frantic::volumetrics::voxel_coord_system
get_default_vcs_for_morton_sort( const frantic::graphics::boundbox3f& bounds,
                                 const frantic::graphics::vector3f& /*mean*/,
                                 const frantic::graphics::vector3f& variance, std::size_t particleCount ) {
    const float maxLength = bounds.get_max_dimension();
    const float edgeLength = maxLength > 0 ? maxLength : 1;

    // Approximate Bounding Sphere heuristic
    const float sizeMultiplier = 1024;
    const frantic::graphics::vector3f& v = variance;
    const float newRadius =
        particleCount > 0
            ? (float)( sizeMultiplier * std::pow( math::pi_value * 4.0 / 3.0 / particleCount, 1.0 / 3.0 ) *
                       std::sqrt( v.x + v.y + v.z ) )
            : 1;

    float voxelLength = std::min(
        newRadius,
        edgeLength ); // The voxel length shouldn't be bigger than the bounding box size.  That's unnecessarily large
    voxelLength =
        std::max( voxelLength, edgeLength / std::pow( 2.f, 20.f ) ); // Give it a nonzero absolute minimum size
    return frantic::volumetrics::voxel_coord_system( bounds.minimum(), voxelLength );
}

/**
 * Convenience function designed for computing a good vcs for the sort_particles_by_morton_code function.
 * This function computes the boundbox of the particles, and a reasonable voxel length.
 */
inline void
get_default_vcs_for_morton_sort( boost::shared_ptr<frantic::particles::streams::particle_istream> particleStream,
                                 frantic::volumetrics::voxel_coord_system& outVcs,
                                 frantic::graphics::boundbox3f& outBoundingBox,
                                 frantic::logging::progress_logger& logger ) {
    frantic::graphics::vector3f mean;
    frantic::graphics::vector3f variance;
    std::size_t pCount;
    const frantic::graphics::boundbox3f particleBoundingBox =
        detail::morton_sorter::compute_particle_bounds( particleStream, mean, variance, pCount, logger );

    outVcs = get_default_vcs_for_morton_sort( particleBoundingBox, mean, variance, pCount );

    const float maxLength = particleBoundingBox.get_max_dimension();
    const float useLength = maxLength > 0 ? maxLength : 1;
    const frantic::graphics::vector3f minBox( particleBoundingBox.minimum() );
    const frantic::graphics::vector3f lengthBox( useLength );
    outBoundingBox = frantic::graphics::boundbox3f( minBox, minBox + lengthBox );
}

/*
 * Sorts the particles by the morton code for each particle.
 */
inline void sort_particles_by_morton_code( frantic::particles::particle_array::iterator& start,
                                           frantic::particles::particle_array::iterator& end,
                                           frantic::channels::channel_accessor<frantic::graphics::vector3f> posAcc,
                                           const frantic::volumetrics::voxel_coord_system& vcs,
                                           frantic::logging::progress_logger& logger ) {
    detail::morton_sorter mortonSorter( posAcc, vcs.world_origin(), vcs.voxel_length() );
    frantic::sort::parallel_sort( start, end, mortonSorter, logger );
}

/*
 * Sorts the particles by the morton code for each particle.
 */
inline void sort_particles_by_morton_code( frantic::particles::particle_array& particles,
                                           const frantic::volumetrics::voxel_coord_system& vcs,
                                           frantic::logging::progress_logger& logger ) {
    if( !particles.has_channel( _T("Position") ) )
        throw std::runtime_error( "The particle array does not contain a 'Position' channel.\n" );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> posAcc =
        particles.get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") );
    frantic::particles::particle_array::iterator iterStart = particles.begin();
    frantic::particles::particle_array::iterator iterEnd = particles.end();
    sort_particles_by_morton_code( iterStart, iterEnd, posAcc, vcs, logger );
}

/**
 * Attempts to determine if the given filename is associated with a csv particle file object
 * Note: this does more than just check the file extension, it will also attempt to load
 * the file to see if it matches any of the other known formats first, since we consider
 * csv to be a default last-resort load
 */
bool is_csv_particle_file( const frantic::tstring& filename );

/**
 * Just checks to see if the extension is ".pts"
 */
bool is_pts_particle_file( const frantic::tstring& filename );

/**
 * Checks if the extension matches any of our known binary formats
 * Please update this if anything new comes up
 */
bool is_binary_particle_file( const frantic::tstring& filename );

} // namespace particles
} // namespace frantic
