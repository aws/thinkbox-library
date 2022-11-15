// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/files/files.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/particles/particle_file_metadata.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <map>

namespace frantic {
namespace particles {
namespace streams {

// A particle_istream that loads RealFlow .bin files. Currently works to version 9.
class realflow_bin_particle_istream : public particle_istream {
    frantic::tstring m_filename;
    frantic::channels::channel_map m_onDiskParticleChannelMap, m_particleChannelMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    boost::int64_t m_particleCount, m_currentParticleIndex;
    float m_scale;
    int m_fileVersion;
    frantic::particles::particle_file_metadata m_metadata;
    // frantic::graphics::coordinate_system::option m_coordinateSystem;

    // std::ifstream	m_fin;
    files::file_ptr m_fin;
    bool m_finOpen;
    long m_finOffset;

    std::vector<char> m_tempParticleBuffer;
    std::vector<char> m_defaultParticleBuffer;

    static const int BIN_MAGIC_SIGNATURE = 0xFABADA;
    static const int BIN_FLUID_NAME_LENGTH = 250;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccessor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_velAccessor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_forceAccessor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_vorticityAccessor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_normalAccessor;

    // Private copy constructor to disable copying
    realflow_bin_particle_istream( const realflow_bin_particle_istream& ); // not implemented

    // Private assignment operator to disable assignment
    realflow_bin_particle_istream& operator=( const realflow_bin_particle_istream& ); // not implemented

    template <typename Type>
    static Type file_read( FILE* in, const frantic::tstring& filename ) {
        Type t;
        if( 1 != std::fread( &t, sizeof( Type ), 1, in ) )
            throw std::runtime_error( "realflow_bin_particle_istream: Failed to read from file \"" +
                                      frantic::strings::to_string( filename ) + "\"." );

        return t;
    }

    /**
     * A helper functoion to get a vector into the desired coordinate system.
     * @param v The vector to swap values for.
     */
    /*void convert_vector_to_current_space( frantic::graphics::vector3f& v ) {
      switch( m_coordinateSystem ) {
        case frantic::graphics::coordinate_system::unspecified: //KSR began by using right-handed-yup as the default for
    unspecified, so I have left that in. case frantic::graphics::coordinate_system::right_handed_yup: v.z = -v.z; break;
        case frantic::graphics::coordinate_system::right_handed_zup:
          std::swap( v.y, v.z );
          break;
        case frantic::graphics::coordinate_system::left_handed_zup:
          std::swap( v.y, v.z );
          v.y = -v.y;
          break;
        default:
          //left_handed_yup, do nothing
          break;
      }
    }*/

    /**
     * This function will convert a particle that is stored in a left-handed Y-up coordinate system to a right-handed
     * Z-up coordinate system.
     * @param p A pointer to the particle to convert. This particle is assumed to described by
     * 'm_onDiskParticleChannelMap'.
     */
    /*void convert_particle_to_current_space( char* p ){
      //Since the matrix form of the conversion is a symmetric, 3x3 row swap; points, vectors and normals get the same
    treatment. convert_vector_to_current_space( m_posAccessor.get( p ) ); convert_vector_to_current_space(
    m_velAccessor.get( p ) ); convert_vector_to_current_space( m_forceAccessor.get( p ) ); if(
    m_vorticityAccessor.is_valid() ) convert_vector_to_current_space( m_vorticityAccessor.get( p ) ); if(
    m_normalAccessor.is_valid() ) convert_vector_to_current_space( m_normalAccessor.get( p ) );
    }*/

    void initialize_stream() {
        if( !m_fin )
            throw std::runtime_error( "realflow_bin_particle_istream: Failed to open file \"" +
                                      frantic::strings::to_string( m_filename ) + "\" for reading." );

        // First check the signature to verify this is actually a .bin file
        boost::int32_t sig = 0;
        std::fread( &sig, sizeof( sig ), 1, m_fin );
        if( sig != BIN_MAGIC_SIGNATURE )
            throw invalid_particle_file_exception()
                << "realflow_bin_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" does not contain the RealFlow .bin file identifier. Got " << sig << " instead.";

        // Skip the fluid name(250 chars) then get the file version
        std::fseek( m_fin, BIN_FLUID_NAME_LENGTH, SEEK_CUR );
        m_fileVersion = static_cast<int>( file_read<boost::int16_t>( m_fin, m_filename ) );

        // get the scene scale (expressed as a ratio of meters)
        m_scale = static_cast<float>( file_read<float>( m_fin, m_filename ) );

        // Skip fluid type, sim time, frame number, fps (4 bytes each, 4 * 4 = 16)
        std::fseek( m_fin, 16, SEEK_CUR );
        m_particleCount = static_cast<boost::int64_t>( file_read<boost::int32_t>( m_fin, m_filename ) );

        // Skip radius, pressure, speed, temp (4 + 3 * 12 = 40)
        std::fseek( m_fin, 40, SEEK_CUR );

        // For versions later than 7 there are a few more fields to skip
        // emitter pos, emitter rotation, emitter scale (3 * 12 = 36)
        if( m_fileVersion >= 7 )
            std::fseek( m_fin, 36, SEEK_CUR );

        m_finOpen = false;
        m_finOffset = std::ftell( m_fin );
        m_fin.close();

        setup_metadata();

        // Create a particle channel map matching the realflow on-disk data structure
        m_onDiskParticleChannelMap.reset();
        m_onDiskParticleChannelMap.define_channel( _T("Position"), 3, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("Velocity"), 3, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("Force"), 3, channels::data_type_float32 );

        if( m_fileVersion >= 9 )
            m_onDiskParticleChannelMap.define_channel( _T("Vorticity"), 3, channels::data_type_float32 );
        if( m_fileVersion >= 3 )
            m_onDiskParticleChannelMap.define_channel( _T("Normal"), 3, channels::data_type_float32 );
        if( m_fileVersion >= 4 )
            m_onDiskParticleChannelMap.define_channel( _T("NeighborCount"), 1, channels::data_type_int32 );
        if( m_fileVersion >= 5 ) {
            m_onDiskParticleChannelMap.define_channel( _T("TextureCoord"), 3, channels::data_type_float32 );
            m_onDiskParticleChannelMap.define_channel( _T("RealflowInfoBits"), 1, channels::data_type_int16 );
        }

        m_onDiskParticleChannelMap.define_channel( _T("Age"), 1, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("IsolationTime"), 1, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("Viscosity"), 1, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("RFDensity"), 1, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("Pressure"), 1, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("Mass"), 1, channels::data_type_float32 );
        m_onDiskParticleChannelMap.define_channel( _T("Temperature"), 1, channels::data_type_float32 );
        if( m_fileVersion >= 12 )
            m_onDiskParticleChannelMap.define_channel( _T("ID"), 1, channels::data_type_int64 );
        else
            m_onDiskParticleChannelMap.define_channel( _T("ID"), 1, channels::data_type_int32 );

        m_onDiskParticleChannelMap.end_channel_definition( 1, true, true );

        m_posAccessor = m_onDiskParticleChannelMap.get_accessor<frantic::graphics::vector3f>( _T("Position") );
        m_velAccessor = m_onDiskParticleChannelMap.get_accessor<frantic::graphics::vector3f>( _T("Velocity") );
        m_forceAccessor = m_onDiskParticleChannelMap.get_accessor<frantic::graphics::vector3f>( _T("Force") );

        if( m_onDiskParticleChannelMap.has_channel( _T("Vorticity") ) )
            m_vorticityAccessor =
                m_onDiskParticleChannelMap.get_accessor<frantic::graphics::vector3f>( _T("Vorticity") );
        if( m_onDiskParticleChannelMap.has_channel( _T("Normal") ) )
            m_normalAccessor = m_onDiskParticleChannelMap.get_accessor<frantic::graphics::vector3f>( _T("Normal") );
    }

    void setup_metadata() {
        using namespace frantic::channels;

        // set file-level metadata
        channel_map newMap;
        prt::add_coordinate_system( newMap );
        prt::length_unit_in_micrometers::add_channel( newMap );
        frantic::channels::property_map generalMetadata;
        newMap.end_channel_definition();
        generalMetadata.set_channel_map_with_swap( newMap );

        prt::set_coordinate_system( generalMetadata, frantic::graphics::coordinate_system::left_handed_yup );
        prt::length_unit_in_micrometers::set_value( generalMetadata, m_scale * 1e6 );

        m_metadata.set_general_metadata( generalMetadata );

        // set channel-specific metadata
        channel_map channelInterpretationMap;
        prt::add_channel_interpretation( channelInterpretationMap );
        channelInterpretationMap.end_channel_definition();

        frantic::channels::property_map positionMetadata;
        frantic::channels::property_map velocityMetadata;
        frantic::channels::property_map forceMetadata;
        frantic::channels::property_map vorticityMetadata;
        frantic::channels::property_map normalMetadata;

        positionMetadata.set_channel_map( channelInterpretationMap );
        velocityMetadata.set_channel_map( channelInterpretationMap );
        forceMetadata.set_channel_map( channelInterpretationMap );

        if( m_onDiskParticleChannelMap.has_channel( _T("Vorticity") ) ) {
            vorticityMetadata.set_channel_map( channelInterpretationMap );
        }

        if( m_onDiskParticleChannelMap.has_channel( _T("Normal") ) ) {
            normalMetadata.set_channel_map( channelInterpretationMap );
        }

        prt::set_channel_interpretation( positionMetadata, prt::channel_interpretation::point );
        prt::set_channel_interpretation( velocityMetadata, prt::channel_interpretation::vector );
        prt::set_channel_interpretation( forceMetadata, prt::channel_interpretation::vector );

        m_metadata.set_channel_metadata( _T("Position"), positionMetadata );
        m_metadata.set_channel_metadata( _T("Velocity"), velocityMetadata );
        m_metadata.set_channel_metadata( _T("Force"), forceMetadata );

        if( m_onDiskParticleChannelMap.has_channel( _T("Vorticity") ) ) {
            prt::set_channel_interpretation( vorticityMetadata, prt::channel_interpretation::vector );
            m_metadata.set_channel_metadata( _T("Vorticity"), vorticityMetadata );
        }

        if( m_onDiskParticleChannelMap.has_channel( _T("Normal") ) ) {
            prt::set_channel_interpretation( normalMetadata, prt::channel_interpretation::vector );
            m_metadata.set_channel_metadata( _T("Normal"), normalMetadata );
        }
    }

  public:
    realflow_bin_particle_istream(
        const frantic::tstring& file ) //, frantic::graphics::coordinate_system::option coordinateSystem =
                                       // frantic::graphics::coordinate_system::right_handed_zup )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) ) //, m_coordinateSystem( coordinateSystem )
    {
        initialize_stream();
        m_currentParticleIndex = -1;
        set_channel_map( m_onDiskParticleChannelMap );
        m_tempParticleBuffer.resize( m_onDiskParticleChannelMap.structure_size() );
    }

    realflow_bin_particle_istream(
        const frantic::tstring& file,
        const frantic::channels::channel_map&
            particleChannelMap ) //, frantic::graphics::coordinate_system::option coordinateSystem =
                                 // frantic::graphics::coordinate_system::right_handed_zup )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) ) //, m_coordinateSystem( coordinateSystem )
    {
        initialize_stream();
        m_currentParticleIndex = -1;
        set_channel_map( particleChannelMap );
        m_tempParticleBuffer.resize( m_onDiskParticleChannelMap.structure_size() );
    }

    virtual ~realflow_bin_particle_istream() { close(); }

    void close() {
        m_fin.close();
        m_particleCount = 0;
    }

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_filename; }

    boost::int64_t particle_count() const { return m_particleCount; }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return m_particleCount - m_currentParticleIndex - 1; }

    boost::int64_t particle_progress_count() const { return m_particleCount; }

    boost::int64_t particle_progress_index() const { return m_currentParticleIndex; }

    const frantic::channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    // Access to the channel_map in the file
    const frantic::channels::channel_map& get_native_channel_map() const { return m_onDiskParticleChannelMap; }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        std::vector<char> newDefaultParticle( particleChannelMap.structure_size() );
        if( m_defaultParticleBuffer.size() > 0 ) {
            frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
            defaultAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticleBuffer[0] );
        } else
            memset( &newDefaultParticle[0], 0, particleChannelMap.structure_size() );
        m_defaultParticleBuffer.swap( newDefaultParticle );

        // Set the map and the adaptor
        m_particleChannelMap = particleChannelMap;
        m_pcmAdaptor.set( m_particleChannelMap, m_onDiskParticleChannelMap );
    }

    void set_default_particle( char* buffer ) {
        m_particleChannelMap.copy_structure( &m_defaultParticleBuffer[0], buffer );
    }

    bool get_particle( char* rawParticleBuffer ) {
        if( !m_finOpen ) {
            m_finOpen = true;
            m_fin.reset( frantic::files::tfopen( m_filename.c_str(), _T("rb") ) );
            if( !m_fin )
                throw std::runtime_error(
                    "realflow_bin_particle_istream.get_particle: Failed to re-open the particle file \"" +
                    frantic::strings::to_string( m_filename ) + "\"" );
            std::fseek( m_fin, m_finOffset, SEEK_SET );
        } else if( !m_fin )
            throw std::runtime_error(
                "realflow_bin_particle_istream.get_particle: Tried to read from particle file \"" +
                frantic::strings::to_string( m_filename ) + "\" after it was already closed." );

        if( particle_index() + 1 < particle_count() ) {
            if( m_pcmAdaptor.is_identity() ) {
                // If we can, just read the particle data straight into the output buffer
                if( 1 != std::fread( rawParticleBuffer, m_onDiskParticleChannelMap.structure_size(), 1, m_fin ) )
                    throw std::runtime_error( "realflow_bin_particle_istream.next_particle: Error reading particle " +
                                              boost::lexical_cast<std::string>( m_currentParticleIndex + 1 ) + " of " +
                                              boost::lexical_cast<std::string>( m_particleCount ) +
                                              " from the file \"" + frantic::strings::to_string( m_filename ) + "\"" );
                // convert_particle_to_current_space( rawParticleBuffer );
            } else {
                // Otherwise, read it into the temp buffer, then use the adaptor to copy it into the output buffer
                if( 1 != std::fread( &m_tempParticleBuffer[0], m_onDiskParticleChannelMap.structure_size(), 1, m_fin ) )
                    throw std::runtime_error( "realflow_bin_particle_istream.next_particle: Error reading particle " +
                                              boost::lexical_cast<std::string>( m_currentParticleIndex + 1 ) + " of " +
                                              boost::lexical_cast<std::string>( m_particleCount ) +
                                              " from the file \"" + frantic::strings::to_string( m_filename ) + "\"" );
                // convert_particle_to_current_space( &m_tempParticleBuffer[0] );
                m_pcmAdaptor.copy_structure( rawParticleBuffer, &m_tempParticleBuffer[0], &m_defaultParticleBuffer[0] );
            }

            // Go to the next particle
            ++m_currentParticleIndex;

            return true;
        } else {
            return false;
        }
    }

    bool get_particles( char* particleBuffer, std::size_t& numParticles ) {
        std::size_t particleSize = m_particleChannelMap.structure_size();
        for( std::size_t i = 0; i < numParticles; ++i ) {
            if( !get_particle( particleBuffer + i * particleSize ) ) {
                numParticles = i;
                return false;
            }
        }

        return true;
    }

    const frantic::particles::particle_file_metadata& get_metadata() const { return m_metadata; }
};

} // namespace streams
} // namespace particles
} // namespace frantic
