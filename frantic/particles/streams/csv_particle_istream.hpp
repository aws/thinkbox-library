// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/streams/particle_istream.hpp>

#include <frantic/files/csv_files.hpp>
#include <frantic/files/files.hpp>

#include <frantic/strings/tstring.hpp>

#include <frantic/channels/channel_column_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/property_map.hpp>

#include <frantic/locale/locale.hpp>

namespace frantic {
namespace particles {
namespace streams {

using frantic::channels::channel_column_map;
using frantic::channels::channel_map;
using frantic::channels::channel_map_adaptor;

// CSV Format for particle file
// Two modes, detected by checking whether all of the values in the column header are floating point numbers:
//  - Without a column header.
//     * Must have 3 or 6 columns
//     * First 3 columns are the Position, next 3 columns are the Color.
//  - With a column header.
//     * General format for column heading is "type TypeName[index]"
//     * If [index] is excluded, then having multiple repeated values will determine the arity.
//       The indices can go from 0 to N-1 or from 1 to N, the code will automatically detect this.
//     * If the type is left out, it will be defaulted to float32.
//     * The individual elements of a multi-value channel must be adjacent and in increasing order.

// A particle_istream that loads .csv files.
class csv_particle_istream : public particle_istream {
  public:
    typedef std::map<frantic::tstring, frantic::channels::property_map> channel_metadata_t;

    enum csv_format_guess {
        CSV_PARTICLE_GUESS_UNKNOWN = 0,
        CSV_PARTICLE_GUESS_NUMBERS,
        CSV_PARTICLE_GUESS_HEADERS,
        CSV_PARTICLE_GUESS_MIXED,
    };

    /**
     * Attempt to guess the underlying format of the columns in the given csv file from its first row.
     */
    static csv_format_guess guess_channel_column_map( const std::vector<std::string>& headers,
                                                      const frantic::tstring& file, bool noError,
                                                      frantic::channels::data_type_t positionTypeHint,
                                                      channel_column_map& mapping );

  private:
    frantic::tstring m_name;
    channel_column_map m_channelColumnMapping;
    channel_map m_particleChannelMap;
    channel_map_adaptor m_pcmAdaptor;
    boost::int64_t m_currentParticleIndex;
    std::vector<char> m_tempParticleBuffer;
    std::vector<char> m_defaultParticleBuffer;
    std::vector<std::string> m_particleData;
    frantic::files::csv_reader m_csvReader;
    char m_delimiter;

    frantic::channels::property_map m_metadata;
    channel_metadata_t m_channelMetadata;

    files::file_ptr m_fin;
    bool m_finReopened;
    int m_headerRowCount;

    // Private copy constructor to disable copying
    csv_particle_istream( const csv_particle_istream& ); // not implemented

    // Private assignment operator to disable assignment
    csv_particle_istream& operator=( const csv_particle_istream& ); // not implemented

    static frantic::files::csv_reader get_csv_reader( FILE* file, const frantic::tstring& streamName, bool useBuffer,
                                                      char delimiter ) {
        frantic::files::csv_reader csvReader( file, streamName, useBuffer );
        csvReader.set_delimiter( delimiter );
        if( delimiter == ' ' || delimiter == '\t' ) {
            csvReader.set_emit_empty_columns( false );
        } else {
            csvReader.set_emit_empty_columns( true );
        }
        return csvReader;
    }

    void initialize_stream( const frantic::tstring& file, const channel_column_map& expectedFormat, char delimiter,
                            int headerRowCount, frantic::channels::data_type_t positionTypeHint, bool noError );
    void initialize_stream( const frantic::tstring& file, char delimiter,
                            frantic::channels::data_type_t positionTypeHint, bool noError );

  public:
    // This loads a particle file, and feeds back the results in the particle's native particle storage layout.
    csv_particle_istream( const frantic::tstring& file, char delimiter = '\0',
                          frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid,
                          bool noError = false ) {
        initialize_stream( file, delimiter, positionTypeHint, noError );
        // Start the particle index one before the first valid index
        m_currentParticleIndex = -1;
        // Initialize the particle channel map that the user of this class sees to match the native one provided by the
        // .csv
        set_channel_map( m_channelColumnMapping.get_channel_map() );
        // Allocate the temporary particle buffer
        m_tempParticleBuffer.resize( m_channelColumnMapping.get_channel_map().structure_size() );
    }

    csv_particle_istream( const frantic::tstring& file, const channel_map& particleChannelMap, char delimiter = '\0',
                          frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid,
                          bool noError = false ) {
        if( positionTypeHint == channels::data_type_invalid && particleChannelMap.has_channel( _T("Position") ) ) {
            size_t positionArity = 0;
            particleChannelMap.get_channel_definition( _T("Position"), positionTypeHint, positionArity );
        }
        initialize_stream( file, delimiter, positionTypeHint, noError );
        m_currentParticleIndex = -1;
        set_channel_map( particleChannelMap );
        m_tempParticleBuffer.resize( m_channelColumnMapping.get_channel_map().structure_size() );
    }

    csv_particle_istream( const frantic::tstring& file, const channel_column_map& columnMap, char delimiter = '\0',
                          int headerRowCount = 0,
                          frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid,
                          bool noError = false ) {
        if( positionTypeHint == channels::data_type_invalid &&
            columnMap.get_channel_map().has_channel( _T("Position") ) ) {
            size_t positionArity = 0;
            columnMap.get_channel_map().get_channel_definition( _T("Position"), positionTypeHint, positionArity );
        }
        initialize_stream( file, columnMap, delimiter, headerRowCount, positionTypeHint, noError );
        m_currentParticleIndex = -1;
        set_channel_map( m_channelColumnMapping.get_channel_map() );
        m_tempParticleBuffer.resize( m_channelColumnMapping.get_channel_map().structure_size() );
    }

    csv_particle_istream( const frantic::tstring& file, const channel_map& particleChannelMap,
                          const channel_column_map& columnMap, char delimiter = '\0', int headerRowCount = 0,
                          frantic::channels::data_type_t positionTypeHint = frantic::channels::data_type_invalid,
                          bool noError = false ) {
        if( positionTypeHint == channels::data_type_invalid && particleChannelMap.has_channel( _T("Position") ) ) {
            size_t positionArity = 0;
            particleChannelMap.get_channel_definition( _T("Position"), positionTypeHint, positionArity );
        }
        initialize_stream( file, columnMap, delimiter, headerRowCount, positionTypeHint, noError );
        m_currentParticleIndex = -1;
        set_channel_map( particleChannelMap );
        m_tempParticleBuffer.resize( m_channelColumnMapping.get_channel_map().structure_size() );
    }

    virtual ~csv_particle_istream() { close(); }

    void close() { m_fin.close(); }

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_name; }

    boost::int64_t particle_count() const {
        // .csv files don't have a header with the particle count embedded in it.
        return -1;
    }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return -1; }

    boost::int64_t particle_progress_count() const {
        // .csv files don't have a header with the particle count embedded in it.
        return m_csvReader.file_progress_count();
    }

    boost::int64_t particle_progress_index() const { return m_csvReader.file_progress_index(); }

    const channel_map& get_channel_map() const { return m_particleChannelMap; }

    // Access to the channel_map in the file
    const channel_map& get_native_channel_map() const { return m_channelColumnMapping.get_channel_map(); }

    void set_default_particle( char* buffer ) {
        m_particleChannelMap.copy_structure( &m_defaultParticleBuffer[0], buffer );
    }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const channel_map& particleChannelMap ) {
        std::vector<char> newDefaultParticle( particleChannelMap.structure_size() );
        if( m_defaultParticleBuffer.size() > 0 ) {
            frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
            defaultAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticleBuffer[0] );
        } else
            memset( &newDefaultParticle[0], 0, particleChannelMap.structure_size() );
        m_defaultParticleBuffer.swap( newDefaultParticle );

        // Set the map and the adaptor
        m_particleChannelMap = particleChannelMap;
        m_pcmAdaptor.set( m_particleChannelMap, m_channelColumnMapping.get_channel_map() );
    }

    bool get_particle( char* rawParticleBuffer ) {
        if( !m_finReopened ) {
            m_finReopened = true;
            m_fin.reset( frantic::files::tfopen( m_name.c_str(), _T("rb") ) );
            if( !m_fin )
                throw std::runtime_error( "csv_particle_istream.get_particle: Could not reopen the file \"" +
                                          frantic::strings::to_string( name() ) + '"' );

            m_csvReader = get_csv_reader( m_fin, name(), true, m_delimiter );
            for( int row = 0; row < m_headerRowCount; ++row ) {
                std::vector<std::string> temp;
                m_csvReader.read_line( temp );
            }
        } else if( !m_fin )
            throw std::runtime_error( "csv_particle_istream.get_particle: Tried to read from particle file \"" +
                                      frantic::strings::to_string( name() ) + "\" after it was already closed." );

        // Load in the .csv line
        bool gotLine = m_csvReader.read_line( m_particleData );

        // If the stream reached the end, then we don't get a particle
        if( !gotLine )
            return false;

        try {
            if( m_particleData.size() < m_channelColumnMapping.column_count() ) {
                throw std::runtime_error( "Unexpected number of columns in particle.  "
                                          "Expected " +
                                          boost::lexical_cast<std::string>( m_channelColumnMapping.column_count() ) +
                                          " columns, "
                                          "but got " +
                                          boost::lexical_cast<std::string>( m_particleData.size() ) + " instead." );
            }

            if( m_pcmAdaptor.is_identity() ) {
                // If we can, dump the data straight to the output
                m_channelColumnMapping.copy_structure_from_strings( rawParticleBuffer, m_particleData );
            } else {
                // Otherwise, read it into the temp buffer, then use the adaptor to copy it into the output buffer
                m_channelColumnMapping.copy_structure_from_strings( &m_tempParticleBuffer[0], m_particleData );
                m_pcmAdaptor.copy_structure( rawParticleBuffer, &m_tempParticleBuffer[0], &m_defaultParticleBuffer[0] );
            }
        } catch( const std::exception& e ) {
            throw std::runtime_error(
                "csv_particle_istream.get_particle: Error reading file \"" + frantic::strings::to_string( name() ) +
                "\", line " + boost::lexical_cast<std::string>( m_csvReader.get_line_number() + 1 ) + ": " + e.what() );
        }

        // Go to the next particle
        ++m_currentParticleIndex;

        return true;
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
};

} // namespace streams
} // namespace particles
} // namespace frantic
