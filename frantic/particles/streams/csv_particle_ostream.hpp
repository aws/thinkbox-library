// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <frantic/particles/streams/particle_ostream.hpp>

#include <frantic/strings/tstring.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>

namespace frantic {
namespace particles {
namespace streams {

//////////////////////
// FRANTIC GENERIC PARTICLES V1 OUTPUT STREAM
//////////////////////
class csv_particle_ostream : public particle_ostream {
    channel_map m_particleChannelMapForFile, m_particleChannelMap;
    channel_map_adaptor m_pcmAdaptor;

    const frantic::tstring m_file;
    std::basic_ofstream<frantic::tchar> m_fout;

    std::vector<char> m_tempParticleBuffer;
    boost::int64_t m_currentParticleIndex, m_expectedParticleCount;

    // Private copy constructor to disable copying
    csv_particle_ostream( const csv_particle_ostream& ); // not implemented

    // Private assignment operator to disable copying
    csv_particle_ostream& operator=( const csv_particle_ostream& ); // not implemented

    void initialize_stream() {
        // fout << "trying to open file stream" << std::endl;
        if( !m_fout )
            throw std::runtime_error( "csv_particle_ostream: Failed to open file \"" +
                                      frantic::strings::to_string( m_file ) + "\" for writing." );

        // Write out the column headers
        // This will look something like "float32 Position[0], float32 Position[1], float32 Position[2], float32
        // Density"
        for( std::size_t i = 0; i < m_particleChannelMapForFile.channel_count(); ++i ) {
            if( i != 0 )
                m_fout << _T(",");
            const channels::channel& pc = m_particleChannelMapForFile[i];
            if( pc.arity() == 1 ) {
                m_fout << channel_data_type_str( pc.data_type() ) << _T( ' ' ) << pc.name();
            } else {
                for( std::size_t j = 0; j < pc.arity(); ++j ) {
                    if( j != 0 )
                        m_fout << _T( ',' );
                    m_fout << channel_data_type_str( pc.data_type() ) << _T( ' ' ) << pc.name() << _T( '[' ) << j
                           << _T( ']' );
                }
            }
        }
        m_fout << std::endl;

        // Resize our temp particle buffer, note that the header removes any padding, so it might be
        // smaller than the channel_map that was passed in for initialization.
        m_tempParticleBuffer.resize( m_particleChannelMapForFile.structure_size() );

        // Start the particle index one before the first valid index
        m_currentParticleIndex = -1;
    }

    void write_csv_particle( const char* rawFileParticleData ) {
        for( std::size_t i = 0; i < m_particleChannelMapForFile.channel_count(); ++i ) {
            if( i != 0 )
                m_fout << _T(",");
            const channels::channel& pc = m_particleChannelMapForFile[i];
            // Get a pointer to the data for this particular channel
            const char* channelData = pc.get_channel_data_pointer( rawFileParticleData );
            // Print out the channel, using a comma as the separator
            channels::channel_data_type_print( m_fout, _T(","), pc.arity(), pc.data_type(), channelData );
        }
        m_fout << _T("\n");
    }

  public:
    csv_particle_ostream( const frantic::tstring& file, const channel_map& particleChannelMap,
                          const channel_map& particleChannelMapForFile, boost::int64_t expectedParticleCount = -1 )
        : m_file( file )
        , m_particleChannelMap( particleChannelMap )
        , m_particleChannelMapForFile( particleChannelMapForFile )
        , m_expectedParticleCount( expectedParticleCount ) {
        // TODO: We need a better strategy for generating a temporary file.
        m_fout.open( ( file + _T(".tmp") ).c_str() );

        // fout << "started making stream" << std::endl;
        initialize_stream();

        // Initialize the adaptor for converting the particle format to the one in the file
        m_pcmAdaptor.set( m_particleChannelMapForFile, m_particleChannelMap );
        // Reset the temporary buffer values to zero
        memset( &m_tempParticleBuffer[0], 0, m_tempParticleBuffer.size() );
    }

    virtual ~csv_particle_ostream() {
        // close();
        if( m_fout.is_open() ) {
            m_fout.close();
            frantic::files::remove_file( ( m_file + _T(".tmp") ).c_str() );
            FF_LOG( warning ) << "~csv_particle_ostream() - The output particle stream \""
                              << frantic::strings::to_tstring( m_file )
                              << "\" was destroyed without being properly closed." << std::endl;
        }
    }

    // This is the particle channel map which specifies the byte layout of the particle structure.
    const channel_map& get_channel_map() const { return m_particleChannelMap; }

    // This allows you to change the particle layout that's being saved on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const channel_map& particleChannelMap ) {
        m_particleChannelMap = particleChannelMap;
        // Initialize the adaptor for converting the particle format to the one in the file
        m_pcmAdaptor.set( m_particleChannelMapForFile, m_particleChannelMap );
        // Reset the temporary buffer values to zero
        memset( &m_tempParticleBuffer[0], 0, m_tempParticleBuffer.size() );

        // fout << "set the particle channel map" << std::endl;
    }

    void close() {
        if( m_fout.is_open() ) {
            // If a valid expected amount was given during construction, throw an exception if the
            // expected amount doesn't match the actual amount.
            if( m_expectedParticleCount >= 0 && m_expectedParticleCount != ( m_currentParticleIndex + 1 ) )
                throw exception_stream()
                    << "csv_particle_ostream.close: Closed a file without writing the specified number of particles. "
                    << m_expectedParticleCount << " expected vs. " << ( m_currentParticleIndex + 1 ) << " written.";

            m_fout.close();
            if( m_fout.bad() ) {
                frantic::files::remove_file( ( m_file + _T(".tmp") ).c_str() );
                throw std::runtime_error( "csv_particle_ostream::close() - Failed to close file: " +
                                          frantic::strings::to_string( m_file ) );
            }

            frantic::files::remove_file( m_file.c_str() ); // Delete the destination so we ensure we can write to it.
            if( 0 != frantic::files::rename_file( m_file + _T(".tmp"), m_file ) )
                throw std::runtime_error(
                    "csv_particle_ostream::close() - Failed to move temporary particle file to its destination, \"" +
                    frantic::strings::to_string( m_file ) + "\": "
#ifdef _WIN32
                    + win32::GetLastErrorMessageA() );
#else
                    + strerror( errno ) );
#endif
        }
    }

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    void put_particle( const char* rawParticleData ) {
        if( !m_fout.is_open() )
            throw std::runtime_error( "csv_particle_ostream.put_particle: Tried to write to particle file \"" +
                                      frantic::strings::to_string( m_file ) + "\" after it was closed." );

        // Make sure that we haven't written more than the expected count
        if( m_expectedParticleCount > 0 && ( m_currentParticleIndex + 1 ) >= m_expectedParticleCount )
            throw std::runtime_error(
                "csv_particle_ostream.put_particle: Tried to write more particles than were specified to the file \"" +
                frantic::strings::to_string( m_file ) + "\".  The specified number of particles is " +
                boost::lexical_cast<std::string>( m_expectedParticleCount ) );

        if( m_pcmAdaptor.is_identity() ) {
            write_csv_particle( rawParticleData );
        } else {
            m_pcmAdaptor.copy_structure( &m_tempParticleBuffer[0], rawParticleData );
            write_csv_particle( &m_tempParticleBuffer[0] );
        }
        ++m_currentParticleIndex;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
