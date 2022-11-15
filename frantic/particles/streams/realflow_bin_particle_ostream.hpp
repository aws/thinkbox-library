// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/particles/streams/particle_ostream.hpp>

namespace frantic {
namespace particles {
namespace streams {

class realflow_bin_particle_ostream : public particle_ostream {
  private:
    frantic::tstring m_file;
    std::ofstream m_fout;
    channel_map m_pcmIn, m_pcmOut;
    channel_map_adaptor m_pcmAdaptor;
    std::vector<char> m_outBuffer;
    boost::int64_t m_index, m_expectedParticleCount;
    frantic::graphics::coordinate_system::option m_coordinateSystem;

    // We need to store the file offset if we are going back to write the particle count
    // at the end. This could be calculated manually, but this way it is resistant to header
    // changes.
    std::ostream::pos_type m_countOffset;

    channels::channel_accessor<boost::int32_t> m_idAccessor;
    channels::channel_accessor<float> m_massAccessor, m_temperatureAccessor, m_densityAccessor;
    channels::channel_accessor<frantic::graphics::vector3f> m_posAccessor, m_velAccessor, m_forceAccessor,
        m_vorticityAccessor, m_normalAccessor;

    static void convert_vector_to_LHYup( frantic::graphics::vector3f& v ) { std::swap( v.y, v.z ); }

    // Disallow copy construction
    realflow_bin_particle_ostream( const realflow_bin_particle_ostream& ) {}

    // Disallow assignment
    realflow_bin_particle_ostream& operator=( const realflow_bin_particle_ostream& ) { return *this; }

  private:
    template <class T>
    void write( std::ofstream& o, T t ) {
        o.write( (char*)&t, sizeof( T ) );
    }

    void write_realflow_bin_header() {
        char fluidName[250];
        memset( fluidName, 0, sizeof( fluidName ) );
        strncpy( fluidName, frantic::strings::to_string( m_file ).c_str(), sizeof( fluidName ) );

        // Extract the frame number from the file
        frantic::tstring prefix, postfix;
        int frameNum = 0, digitCount = 0;
        files::split_sequence_path( m_file, prefix, digitCount, frameNum, postfix );

        write<boost::uint32_t>( m_fout, 0xFABADA );
        m_fout.write( fluidName, sizeof( fluidName ) );
        write<boost::uint16_t>( m_fout, 9 );        // Version
        write<float>( m_fout, 1.f );                // Scale
        write<boost::uint32_t>( m_fout, 8 );        // Fluid type
        write<float>( m_fout, 0.f );                // Elapsed simulation time
        write<boost::uint32_t>( m_fout, frameNum ); // Frame number
        write<boost::uint32_t>( m_fout, 24 );       // FPS
        m_countOffset = m_fout.tellp();
        write<boost::uint32_t>( m_fout,
                                static_cast<boost::uint32_t>( m_expectedParticleCount ) ); // Number of particles
        write<float>( m_fout, 1.f );                                                       // Radius
        write<vector3f>( m_fout, vector3f( 1.f, 0.f, 0.5f ) );                             // Pressure
        write<vector3f>( m_fout, vector3f( 1.f, 0.f, 0.5f ) );                             // Speed
        write<vector3f>( m_fout, vector3f( 1.f, 0.f, 0.5f ) );                             // Temperature
        write<vector3f>( m_fout, vector3f( 0.f, 0.f, 0.f ) );                              // Emitter pos
        write<vector3f>( m_fout, vector3f( 0.f, 0.f, 1.f ) );                              // Emitter rot
        write<vector3f>( m_fout, vector3f( 1.f, 1.f, 1.f ) );                              // Emitter scale
    }

    void create_realflow_bin_channel_map( channel_map& outPcm ) {
        outPcm.reset();
        outPcm.define_channel( _T("Position"), 3, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Velocity"), 3, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Force"), 3, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Vorticity"), 3, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Normal"), 3, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("NeighbourCount"), 1, frantic::channels::data_type_int32 );
        outPcm.define_channel( _T("TextureCoord"), 3, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("RealflowInfoBits"), 1, frantic::channels::data_type_int16 );
        outPcm.define_channel( _T("Age"), 1, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("IsolationTime"), 1, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Viscosity"), 1, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Density"), 1, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Pressure"), 1, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Mass"), 1, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("Temperature"), 1, frantic::channels::data_type_float32 );
        outPcm.define_channel( _T("ID"), 1, frantic::channels::data_type_int32 );
        outPcm.end_channel_definition( 1, true, true );
    }

    void convert_vector_from_current_space( frantic::graphics::vector3f& v ) {
        switch( m_coordinateSystem ) {
        case frantic::graphics::coordinate_system::unspecified: // KSR began by using right-handed-yup as the default
                                                                // for unspecified, so I have left that in.
        case frantic::graphics::coordinate_system::right_handed_yup:
            v.z = -v.z;
            break;
        case frantic::graphics::coordinate_system::right_handed_zup:
            std::swap( v.y, v.z );
            break;
        case frantic::graphics::coordinate_system::left_handed_zup:
            v.y = -v.y;
            std::swap( v.y, v.z );
            break;
        default:
            // left_handed_yup, do nothing
            break;
        }
    }

  public:
    realflow_bin_particle_ostream( const frantic::tstring& file, const channel_map& pcm,
                                   boost::int64_t expectedParticleCount = -1,
                                   frantic::graphics::coordinate_system::option coordinateSystem =
                                       frantic::graphics::coordinate_system::right_handed_zup )
        : m_file( file )
        , m_fout( ( file + _T(".tmp") ).c_str(), std::ios::out | std::ios::binary )
        , m_expectedParticleCount( expectedParticleCount )
        , m_index( 0 )
        , m_countOffset( 0 )
        , m_coordinateSystem( coordinateSystem ) {
        if( !m_fout )
            throw std::runtime_error( "realflow_bin_particle_ostream: Unable to open file: " +
                                      frantic::strings::to_string( file ) + " for output." );

        create_realflow_bin_channel_map( m_pcmOut );
        m_idAccessor = m_pcmOut.get_accessor<boost::int32_t>( _T("ID") );
        m_massAccessor = m_pcmOut.get_accessor<float>( _T("Mass") );
        m_temperatureAccessor = m_pcmOut.get_accessor<float>( _T("Temperature") );
        m_densityAccessor = m_pcmOut.get_accessor<float>( _T("Density") );
        m_posAccessor = m_pcmOut.get_accessor<frantic::graphics::vector3f>( _T("Position") );
        m_velAccessor = m_pcmOut.get_accessor<frantic::graphics::vector3f>( _T("Velocity") );
        m_forceAccessor = m_pcmOut.get_accessor<frantic::graphics::vector3f>( _T("Force") );
        m_vorticityAccessor = m_pcmOut.get_accessor<frantic::graphics::vector3f>( _T("Vorticity") );
        m_normalAccessor = m_pcmOut.get_accessor<frantic::graphics::vector3f>( _T("Normal") );
        m_outBuffer.resize( m_pcmOut.structure_size() );

        write_realflow_bin_header();

        set_channel_map( pcm );
    }

    virtual ~realflow_bin_particle_ostream() {
        if( m_fout.is_open() ) {
            m_fout.close();
            frantic::files::remove_file( ( m_file + _T(".tmp") ).c_str() );
            FF_LOG( warning ) << "~realflow_bin_particle_ostream() - Destroyed before the stream was closed"
                              << std::endl;
        }
    }

    const channel_map& get_channel_map() const { return m_pcmIn; }

    void set_channel_map( const channel_map& pcm ) {
        m_pcmIn = pcm;
        m_pcmAdaptor.set( m_pcmOut, m_pcmIn );
        memset( &m_outBuffer[0], 0, m_outBuffer.size() );
    }

    std::size_t particle_size() const { return m_pcmIn.structure_size(); }

    void close() {
        if( m_fout.is_open() ) {
            // The .bin files created by RealFlow have 5 additional 0 bytes at the end.
            // RealFlow hangs if you don't put them there.
            char finalBytes[] = { 0, 0, 0, 0, 0 };
            m_fout.write( finalBytes, 5 );

            if( m_expectedParticleCount < 0 ) {
                boost::uint32_t actualCount = static_cast<boost::uint32_t>( m_index );
                m_fout.seekp( m_countOffset, std::ios::beg );
                write<boost::uint32_t>( m_fout, actualCount );
            } else if( m_expectedParticleCount != m_index )
                throw std::runtime_error(
                    "realflow_bin_particle_ostream: Did not write the expected number of particles to file \"" +
                    frantic::strings::to_string( m_file ) +
                    "\". Wrote: " + boost::lexical_cast<std::string>( m_index ) + " instead of " +
                    boost::lexical_cast<std::string>( m_expectedParticleCount ) + " particles." );

            if( !m_fout )
                throw std::runtime_error(
                    "realflow_bin_particle_ostream: Failed to write final bytes to output stream \"" +
                    frantic::strings::to_string( m_file ) + "\"." );

            m_fout.close();
            if( m_fout.bad() ) {
                frantic::files::remove_file( ( m_file + _T(".tmp") ).c_str() );
                throw std::runtime_error( "realflow_bin_particle_ostream.close: Failed to close the output stream \"" +
                                          frantic::strings::to_string( m_file ) + "\"." );
            }

            frantic::files::remove_file( m_file.c_str() );
            // if( !MoveFile((m_file+".tmp").c_str(), m_file.c_str()) )
            if( 0 != frantic::files::rename_file( ( m_file + _T(".tmp") ).c_str(), m_file.c_str() ) )
                throw std::runtime_error(
                    "realflow_bin_particle_ostream::close() - Failed to move the temporary file \"" +
                    frantic::strings::to_string( m_file + _T(".tmp") ) + "\" to destination \"" +
                    frantic::strings::to_string( m_file ) + "\": "
#ifdef _WIN32
                    + win32::GetLastErrorMessageA() );
#else
                    + strerror( errno ) );
#endif
        }
    }

    void put_particle( const char* rawParticleData ) {
        if( m_index == m_expectedParticleCount )
            throw std::runtime_error( "realflow_bin_particle_ostream: wrote too many particles to file \"" +
                                      frantic::strings::to_string( m_file ) + "\"." );

        // Set the ID before copying the particle, so that if the input has no ID, it gets set to incremental
        // values, but if it does have an ID that gets propagated.
        m_idAccessor( m_outBuffer ) = ( boost::int32_t )( m_index + 1 );
        // Also set the some other properties to default values
        m_massAccessor( m_outBuffer ) = 1;
        m_temperatureAccessor( m_outBuffer ) = 300;
        m_densityAccessor( m_outBuffer ) = 1;

        m_pcmAdaptor.copy_structure( &m_outBuffer[0], rawParticleData );

        // Convert from right-handed, Z-up coordinates to left-handed, Y-up coordinates.
        convert_vector_from_current_space( m_posAccessor.get( m_outBuffer ) );
        convert_vector_from_current_space( m_velAccessor.get( m_outBuffer ) );
        convert_vector_from_current_space( m_forceAccessor.get( m_outBuffer ) );
        convert_vector_from_current_space( m_vorticityAccessor.get( m_outBuffer ) );
        convert_vector_from_current_space( m_normalAccessor.get( m_outBuffer ) );

        m_fout.write( &m_outBuffer[0], m_outBuffer.size() );
        if( !m_fout )
            throw std::runtime_error( "realflow_bin_particle_ostream: Failed to write particle to output stream \"" +
                                      frantic::strings::to_string( m_file ) + "\"." );

        m_index++;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
