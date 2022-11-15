// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/algorithm/string/trim.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/files/csv_files.hpp>
#include <frantic/files/files.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace particles {
namespace streams {

/**
 * PTX File Format Reader
 */
class ptx_particle_istream : public particle_istream {
    frantic::tstring m_name;
    frantic::files::file_ptr m_fin;
    bool m_finReopened;
    frantic::files::csv_reader m_csvReader;
    std::vector<std::string> m_particleData;
    std::size_t m_expectedParticleColumnCount;

    frantic::channels::channel_map m_particleChannelMap;
    frantic::channels::channel_map m_onDiskParticleChannelMap;
    frantic::channels::channel_map m_nativeParticleChannelMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> m_positionAcc;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_colorAcc;
    frantic::channels::channel_accessor<boost::uint32_t> m_scannerIndexAccessor;

    boost::int64_t m_currentParticleIndex;

    boost::int64_t m_blockParticleIndex;
    boost::int64_t m_blockParticleCount;

    std::vector<char> m_tempParticleBuffer;
    std::vector<char> m_defaultParticleBuffer;

    bool m_applyTransform;
    frantic::graphics::transform4fd m_transform;
    std::vector<frantic::graphics::transform4fd> m_scannerTransforms;

    boost::int64_t m_particleProgressCount;

    frantic::particles::particle_file_metadata m_metadata;

    static void remove_comments( std::vector<std::string>& data ) {
        std::vector<std::string>::iterator beginErase = data.end();

        for( std::vector<std::string>::iterator col = data.begin(); col != data.end(); ++col ) {
            // look for '#'
            std::string::size_type i = col->find( '#' );
            if( i != std::string::npos ) {
                col->erase( i );
                boost::algorithm::trim( *col, std::locale::classic() );
                if( col->size() > 0 ) {
                    beginErase = col + 1;
                } else {
                    beginErase = col;
                }
            }
        }

        data.erase( beginErase, data.end() );
    }

    static bool try_read_non_empty_line( frantic::files::csv_reader& csvReader, std::vector<std::string>& out ) {
        while( csvReader.read_line( out ) ) {
            remove_comments( out );
            if( out.size() > 0 ) {
                return true;
            }
        }
        return false;
    }

    static std::string get_file_error_message( const frantic::tstring& filename, const std::string& message ) {
        return "ptx_particle_istream: In file '" + frantic::strings::to_string( filename ) + "': " + message;
    }

    static std::string get_line_number_string( const frantic::files::csv_reader& csvReader ) {
        return boost::lexical_cast<std::string>( csvReader.get_line_number() + 1 );
    }

    static std::string get_line_error_message( const frantic::tstring& filename,
                                               const frantic::files::csv_reader& csvReader,
                                               const std::string& message ) {
        return "ptx_particle_istream: In file '" + frantic::strings::to_string( filename ) + "', line " +
               get_line_number_string( csvReader ) + ": " + message;
    }

    static bool try_read_header( const frantic::tstring& filename, frantic::files::csv_reader& csvReader,
                                 boost::int64_t& outParticleCount, frantic::graphics::transform4fd& outTransform ) {
        try {
            std::vector<std::string> line;
            if( !try_read_non_empty_line( csvReader, line ) ) {
                return false;
            }
            if( line.size() != 1 ) {
                throw frantic::invalid_particle_file_exception(
                    "Unexpected characters in first dimension of scan size" );
            }
            const boost::int64_t rowCount = boost::lexical_cast<boost::int64_t>( line[0] );
            if( rowCount < 0 ) {
                throw frantic::invalid_particle_file_exception( "First dimension of scan size is a negative number" );
            }
            if( !try_read_non_empty_line( csvReader, line ) ) {
                throw frantic::invalid_particle_file_exception( "Unable to read second dimension of scan size" );
            }
            if( line.size() != 1 ) {
                throw frantic::invalid_particle_file_exception(
                    "Unexpected characters in second dimension of scan size" );
            }
            const boost::int64_t columnCount = boost::lexical_cast<boost::int64_t>( line[0] );
            if( columnCount < 0 ) {
                throw frantic::invalid_particle_file_exception( "Second dimension of scan size is a negative number" );
            }
            outParticleCount = rowCount * columnCount;

            // xyz offset
            if( try_read_non_empty_line( csvReader, line ) ) {
                if( line.size() != 3 ) {
                    throw frantic::invalid_particle_file_exception( "Unexpected number of columns in XYZ offset line" );
                }
            } else {
                throw frantic::invalid_particle_file_exception( "Unable to read XYZ offset" );
            }

            // transformation matrix (3x3)
            // validate and discard for now
            for( std::size_t i = 0; i < 3; ++i ) {
                if( try_read_non_empty_line( csvReader, line ) ) {
                    if( line.size() != 3 ) {
                        throw frantic::invalid_particle_file_exception(
                            "Unexpected number of columns in 3x3 transform line " +
                            boost::lexical_cast<std::string>( i + 1 ) );
                    }
                } else {
                    throw frantic::invalid_particle_file_exception( "Unable to read line " +
                                                                    boost::lexical_cast<std::string>( i + 1 ) +
                                                                    " of 3x3 transform matrix" );
                }
            }

            // angular transformation matrix (4x4)
            for( int i = 0; i < 4; ++i ) {
                if( try_read_non_empty_line( csvReader, line ) ) {
                    if( line.size() == 4 ) {
                        for( int j = 0; j < 4; ++j ) {
                            bool success = true;
                            double val = 0;
                            try {
                                val = boost::lexical_cast<double>( line[j] );
                            } catch( boost::bad_lexical_cast& ) {
                                success = false;
                            }
                            if( !success ) {
                                throw frantic::invalid_particle_file_exception( "Unable to convert " + line[j] +
                                                                                " to float in transform matrix" );
                            }
                            outTransform.set( i, j, val );
                        }
                    } else {
                        throw frantic::invalid_particle_file_exception(
                            "Unexpected number of columns in 4x4 transform line " +
                            boost::lexical_cast<std::string>( i + 1 ) );
                    }
                } else {
                    throw frantic::invalid_particle_file_exception( "Unable to read line " +
                                                                    boost::lexical_cast<std::string>( i + 1 ) +
                                                                    " of 4x4 transform matrix" );
                }
            }

            return true;
        } catch( const frantic::invalid_particle_file_exception& e ) {
            throw frantic::invalid_particle_file_exception( get_line_error_message( filename, csvReader, e.what() ) );
        } catch( const std::exception& e ) {
            throw std::runtime_error( get_line_error_message( filename, csvReader, e.what() ) );
        }
    }

    static void read_header( const frantic::tstring& filename, frantic::files::csv_reader& csvReader,
                             boost::int64_t& outParticleCount, frantic::graphics::transform4fd& outTransform ) {
        const bool done = try_read_header( filename, csvReader, outParticleCount, outTransform );
        if( !done ) {
            throw frantic::invalid_particle_file_exception(
                get_file_error_message( filename, "Missing PTX file header" ) );
        }
    }

    void assert_sane_header( const frantic::tstring& filename, FILE* fin ) {
        // First, make sure the file seems valid
        // It looks like PTX may support comments, which begin with '#'
        // We're looking for a line with a single integer
        // And we'll ignore comments
        bool gotHeaderInt = false;
        bool doneHeaderInt = false;
        int col = 0;
        bool inComment = false;
        for( ;; ) {
            const int c = fgetc( fin );
            if( c == EOF ) {
                if( gotHeaderInt ) {
                    doneHeaderInt = true;
                }
                break;
            }

            if( !inComment ) {
                if( isdigit( c ) ) {
                    if( doneHeaderInt ) {
                        throw frantic::invalid_particle_file_exception(
                            get_file_error_message( filename, "Unexpected digits after integer in scan size line." ) );
                    } else {
                        gotHeaderInt = true;
                    }
                } else if( c == '#' ) {
                    if( gotHeaderInt ) {
                        doneHeaderInt = true;
                    }
                    inComment = true;
                } else if( isspace( c ) ) {
                    if( gotHeaderInt ) {
                        doneHeaderInt = true;
                    }
                } else {
                    throw frantic::invalid_particle_file_exception( get_file_error_message(
                        filename,
                        "Unexpected characters in scan size line.  Expected only an integer in the scan size line." ) );
                }
            }

            if( c == '\n' ) {
                if( gotHeaderInt ) {
                    doneHeaderInt = true;
                    break;
                }

                inComment = false;
                col = 0;
            }

            if( col > 1024 ) {
                throw frantic::invalid_particle_file_exception(
                    get_file_error_message( filename, "Found unexpectedly long line with > 1024 characters." ) );
            }
        }
        if( !gotHeaderInt ) {
            throw frantic::invalid_particle_file_exception(
                get_file_error_message( filename, "Missing scan size integer." ) );
        }
        if( gotHeaderInt && !doneHeaderInt ) {
            throw frantic::invalid_particle_file_exception(
                get_file_error_message( filename, "Too many characters in scan size line." ) );
        }
    }

    void initialize_stream( const frantic::tstring& filename, channels::data_type_t positionTypeHint ) {
        m_finReopened = false;
        m_scannerTransforms.clear();

        frantic::files::file_ptr fin( frantic::files::tfopen( filename.c_str(), _T("rb") ) );
        if( !fin ) {
            throw std::runtime_error( get_file_error_message( filename, "Unable to open file for reading." ) );
        }

        assert_sane_header( filename, fin );

        // Reset the file
        std::clearerr( fin );
        std::fseek( fin, 0, SEEK_SET );

        frantic::files::csv_reader csvReader( fin, filename );
        csvReader.set_delimiter( ' ' );
        csvReader.set_emit_empty_columns( false );

        m_particleProgressCount = csvReader.file_progress_count();

        std::vector<std::string> data;
        boost::int64_t blockParticleCount = -1;
        frantic::graphics::transform4fd xform;

        while( blockParticleCount <= 0 ) {
            read_header( name(), csvReader, blockParticleCount, xform );
        }

        const bool gotLine = try_read_non_empty_line( csvReader, data );

        if( gotLine ) {
            positionTypeHint =
                ( positionTypeHint == channels::data_type_invalid ) ? channels::data_type_float32 : positionTypeHint;
            if( data.size() == 3 ) {
                // XYZ
                m_onDiskParticleChannelMap.define_channel( _T("Position"), 3, positionTypeHint );
            } else if( data.size() == 4 ) {
                // XYZI
                m_onDiskParticleChannelMap.define_channel( _T("Position"), 3, positionTypeHint );
                m_onDiskParticleChannelMap.define_channel( _T("Intensity"), 1, channels::data_type_float32 );
            } else if( data.size() == 6 ) {
                // XYZRGB
                m_onDiskParticleChannelMap.define_channel( _T("Position"), 3, positionTypeHint );
                m_onDiskParticleChannelMap.define_channel( _T("Color"), 3, channels::data_type_float16 );
            } else if( data.size() == 7 ) {
                // XYZIRGB
                m_onDiskParticleChannelMap.define_channel( _T("Position"), 3, positionTypeHint );
                m_onDiskParticleChannelMap.define_channel( _T("Intensity"), 1, channels::data_type_float32 );
                m_onDiskParticleChannelMap.define_channel( _T("Color"), 3, channels::data_type_float16 );
            } else {
                throw frantic::invalid_particle_file_exception( get_line_error_message(
                    filename, csvReader,
                    "Unexpected number of columns: " + boost::lexical_cast<std::string>( data.size() ) +
                        ".  Possible columns are: 3 columns (X, Y, Z), 4 (X, Y, Z, Intensity), 6 (X, Y, "
                        "Z, R, G, B) or 7 (X, Y, Z, Intensity, R, G, B)." ) );
            }
        } else {
            // default to XYZ ?
            m_onDiskParticleChannelMap.define_channel( _T("Position"), 3, positionTypeHint );
        }

        m_onDiskParticleChannelMap.end_channel_definition( 1, true );

        m_expectedParticleColumnCount = 0;
        for( std::size_t i = 0; i < m_onDiskParticleChannelMap.channel_count(); ++i ) {
            m_expectedParticleColumnCount += m_onDiskParticleChannelMap[i].arity();
        }

        m_nativeParticleChannelMap = m_onDiskParticleChannelMap;
        m_nativeParticleChannelMap.append_channel( _T("ScannerIndex"), 1, channels::data_type_uint32 );
        m_scannerIndexAccessor = m_nativeParticleChannelMap.get_accessor<boost::uint32_t>( _T("ScannerIndex") );
        m_positionAcc = m_nativeParticleChannelMap.get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
        if( m_nativeParticleChannelMap.has_channel( _T("Color") ) ) {
            m_colorAcc = m_nativeParticleChannelMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Color") );
        } else {
            m_colorAcc.reset();
        }
    }

    void initialize_metadata() {
        frantic::channels::property_map generalMetadata;
        prt::length_unit_in_micrometers::add_channel( generalMetadata );

        // Leica states that "The coordinate unit is always in meters."
        prt::length_unit_in_micrometers::set_value( generalMetadata, 1e6 );

        if( get_channel_map().has_channel( _T( "Intensity" ) ) ) {
            frantic::channels::property_map intensityMetadata;
            prt::add_channel_range_property( intensityMetadata, std::make_pair( 0.0, 1.0 ) );
            m_metadata.set_channel_metadata( _T( "Intensity" ), intensityMetadata );
        }

        m_metadata.set_general_metadata( generalMetadata );
    }

    /** Applies any post-loading scaling or transformation to the point while in the native channel map */
    void postprocess_native_particle( char* buffer ) {
        m_scannerIndexAccessor.get( buffer ) = (unsigned int)m_scannerTransforms.size();
        if( m_applyTransform ) {
            m_positionAcc.set( buffer, m_transform * m_positionAcc.get( buffer ) );
        }
        if( m_colorAcc.is_valid() ) {
            m_colorAcc.set( buffer, m_colorAcc.get( buffer ) / 255.f );
        }
    }

  public:
    ptx_particle_istream( const frantic::tstring& filename, bool applyTransform = true,
                          channels::data_type_t positionTypeHint = channels::data_type_invalid )
        : m_name( filename )
        , m_currentParticleIndex( -1 )
        , m_blockParticleIndex( -1 )
        , m_blockParticleCount( -1 )
        , m_finReopened( false )
        , m_applyTransform( applyTransform )
        , m_expectedParticleColumnCount( 0 ) {
        initialize_metadata();
        initialize_stream( m_name, positionTypeHint );
        set_channel_map( m_nativeParticleChannelMap );
        m_tempParticleBuffer.resize( m_nativeParticleChannelMap.structure_size() );
    }

    virtual ~ptx_particle_istream() { close(); }

    void close() { m_fin.close(); }

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_name; }

    boost::int64_t particle_count() const { return -1; }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return -1; }

    boost::int64_t particle_progress_count() const { return m_particleProgressCount; }

    boost::int64_t particle_progress_index() const { return m_csvReader.file_progress_index(); }

    const frantic::channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeParticleChannelMap; }

    void set_default_particle( char* buffer ) {
        m_particleChannelMap.copy_structure( &m_defaultParticleBuffer[0], buffer );
    }

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) {
        std::vector<char> newDefaultParticle( particleChannelMap.structure_size() );
        if( newDefaultParticle.size() > 0 ) {
            if( m_defaultParticleBuffer.size() > 0 ) {
                frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
                defaultAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticleBuffer[0] );
            } else {
                memset( &newDefaultParticle[0], 0, newDefaultParticle.size() );
            }
        }

        m_defaultParticleBuffer.swap( newDefaultParticle );

        // Set the map and the adaptor
        m_particleChannelMap = particleChannelMap;
        m_pcmAdaptor.set( m_particleChannelMap, m_nativeParticleChannelMap );
    }

    bool get_particle( char* rawParticleBuffer ) {
        if( !m_finReopened ) {
            m_finReopened = true;
            m_fin.reset( frantic::files::tfopen( m_name.c_str(), _T("rb") ) );
            if( !m_fin ) {
                throw std::runtime_error( get_file_error_message( name(), "Could not reopen file." ) );
            }

            m_csvReader = frantic::files::csv_reader( m_fin, name(), true );
            m_csvReader.set_delimiter( ' ' );
            m_csvReader.set_emit_empty_columns( false );

            read_header( m_name, m_csvReader, m_blockParticleCount, m_transform );
            m_scannerTransforms.push_back( m_transform );

        } else if( !m_fin ) {
            throw std::runtime_error(
                get_file_error_message( name(), "Tried to read particle from file after it was already closed." ) );
        }

        while( m_blockParticleCount <= 0 || m_blockParticleIndex + 1 >= m_blockParticleCount ) {
            m_blockParticleCount = -1;
            m_blockParticleIndex = -1;

            const bool gotHeader = try_read_header( m_name, m_csvReader, m_blockParticleCount, m_transform );
            if( !gotHeader ) {
                return false;
            }
            m_scannerTransforms.push_back( m_transform );
        }

        bool gotLine = try_read_non_empty_line( m_csvReader, m_particleData );
        if( gotLine ) {
            try {
                if( m_particleData.size() != m_expectedParticleColumnCount ) {
                    throw std::runtime_error( "Unexpected number of columns in particle.  "
                                              "Expected " +
                                              boost::lexical_cast<std::string>( m_expectedParticleColumnCount ) +
                                              " columns, "
                                              "but got " +
                                              boost::lexical_cast<std::string>( m_particleData.size() ) + " instead." );
                }
                if( m_pcmAdaptor.is_identity() ) {
                    // If we can, dump the data straight to the output
                    m_onDiskParticleChannelMap.copy_structure_from_strings( rawParticleBuffer, m_particleData );
                    postprocess_native_particle( rawParticleBuffer );
                } else {
                    // Otherwise, read it into the temp buffer, then use the adaptor to copy it into the output buffer
                    m_onDiskParticleChannelMap.copy_structure_from_strings( &m_tempParticleBuffer[0], m_particleData );
                    postprocess_native_particle( &m_tempParticleBuffer[0] );
                    m_pcmAdaptor.copy_structure( rawParticleBuffer, &m_tempParticleBuffer[0],
                                                 &m_defaultParticleBuffer[0] );
                }
            } catch( const std::exception& e ) {
                throw std::runtime_error( get_line_error_message( name(), m_csvReader, e.what() ) );
            }

            ++m_currentParticleIndex;
            ++m_blockParticleIndex;
            return true;
        } else {
            return false;
        }
    }

    bool get_particles( char* particleBuffer, std::size_t& numParticles ) {
        const std::size_t particleSize = m_particleChannelMap.structure_size();

        for( std::size_t i = 0; i < numParticles; ++i ) {
            if( !get_particle( particleBuffer + i * particleSize ) ) {
                numParticles = i;
                return false;
            }
        }

        return true;
    }

    const frantic::particles::particle_file_metadata& get_metadata() const { return m_metadata; }

    const std::vector<frantic::graphics::transform4fd>& get_seen_transforms() const { return m_scannerTransforms; }
};
} // namespace streams
} // namespace particles
} // namespace frantic
