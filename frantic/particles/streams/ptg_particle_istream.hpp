// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/files/files.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

/**
 *A particle_istream that loads PTG files. Currently works to version 1.
 *
 * http://www.xdesy.de/freeware/PTG-DLL/PTG-1.0.pdf
 */
class ptg_particle_istream : public particle_istream {
  public:
    enum properties { PTG_POSITION_AS_FLOAT = 1, PTG_POSITION_AS_DOUBLE = 2, PTG_INTENSITY = 4, PTG_COLOR = 8 };

  private:
    frantic::tstring m_filename;
    frantic::channels::channel_map m_onDiskParticleChannelMap, m_particleChannelMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    boost::int64_t m_particleCountGuess, m_currentParticleIndex;
    boost::int32_t m_fileVersion;

    frantic::particles::particle_file_metadata m_metadata;

    files::file_ptr m_fin;
    bool m_finOpen;

    std::vector<char> m_tempParticleBuffer;
    std::vector<char> m_defaultParticleBuffer;

    static const boost::int32_t PTG_MAGIC_NUMBER = 0x928FA3C7;

    boost::int32_t m_numCols;
    boost::int32_t m_numRows;
    bool m_hasTransform;
    frantic::graphics::transform4fd m_transform;
    boost::int32_t m_properties;

    std::vector<boost::int64_t> m_columnOffsets;
    boost::int32_t m_curCol;
    boost::int32_t m_curRowPoint;

    boost::int32_t m_numValidPointsForCurCol;

    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> m_posAccessor;
    frantic::channels::channel_cvt_accessor<float> m_intensityAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_colorAccessor;
    frantic::channels::channel_accessor<boost::uint32_t> m_scannerIndexAccessor;

    // Private copy constructor to disable copying
    ptg_particle_istream( const ptg_particle_istream& ); // not implemented

    // Private assignment operator to disable assignment
    ptg_particle_istream& operator=( const ptg_particle_istream& ); // not implemented

    template <typename Type>
    static Type file_read( FILE* in, const frantic::tstring& filename ) {
        Type t;
        if( 1 != std::fread( &t, sizeof( Type ), 1, in ) )
            throw std::runtime_error( "ptg_particle_istream: Failed to read from file \"" +
                                      frantic::strings::to_string( filename ) + "\"." );

        return t;
    }

    // A helper function to load the next col
    void load_next_column() {
        do {
            ++m_curCol;
            m_curRowPoint = 0;

            if( m_curCol < m_numCols ) {
                // get next row
                int errorFlag = std::fseek( m_fin, (long)m_columnOffsets[m_curCol], SEEK_SET );
                if( errorFlag )
                    throw std::runtime_error(
                        "ptg_particle_istream.next_particle: Error seeking row for particle col " +
                        boost::lexical_cast<std::string>( m_curCol ) + " from the file \"" +
                        frantic::strings::to_string( m_filename ) + "\"" );
                boost::int32_t bitmaskSize = ( m_numRows + 7 ) / 8; // [ceil(NumRows/8) bytes]
                std::vector<char> bitmask( bitmaskSize );
                if( bitmask.size() > 0 )
                    if( 1 != std::fread( &bitmask[0], bitmaskSize, 1, m_fin ) )
                        throw std::runtime_error( "ptg_particle_istream.next_particle: Error reading particle col " +
                                                  boost::lexical_cast<std::string>( m_curCol ) + " from the file \"" +
                                                  frantic::strings::to_string( m_filename ) + "\"" );

                // count the 1 bits
                m_numValidPointsForCurCol = 0;
                for( boost::int32_t i = 0; i < bitmaskSize; ++i ) {
                    for( int k = 0; k < 8; ++k, bitmask[i] >>= 1 ) {
                        if( bitmask[i] & 1 )
                            ++m_numValidPointsForCurCol;
                    }
                }
            } else {
                return;
            }
        } while( m_numValidPointsForCurCol == 0 );
    }
    /**
     * A helper function to apply the transform if it exists
     * @oaram particle: the particle to transform
     */
    void apply_transform( char* particle ) {
        if( m_hasTransform ) {
            m_posAccessor.set( particle, m_transform * m_posAccessor.get( particle ) );
        }
    }

    /**
     * A helper function to get the next string from m_fin
     * @return the next string in the ptg file
     */
    std::string get_next_string( FILE* in ) {
        boost::int32_t stringSize = file_read<boost::int32_t>( in, m_filename );
        ++stringSize; // need the +1 so it also reads the \0
        if( stringSize < 1 )
            throw std::runtime_error( "ptg_particle_istream:: Invalid string length in file \"" +
                                      frantic::strings::to_string( m_filename ) + "\"." );
        std::vector<char> keyBuff( stringSize );
        fgets( &keyBuff[0], stringSize, in );
        return std::string( &keyBuff[0] );
    }

    /**
     * A helper function to get move the file to directly after a header key
     * @param key: the key you are looking for (e.g.  "%%transform" )
     * @param in: the file to move
     * @return True if the key was found and the file is now at directly after the key. If false the file will be at the
     * end of the header.
     */
    bool seek_header_key( std::string key, FILE* in ) {
        std::fseek( in, 0, SEEK_SET );
        std::string currkey = "";
        char c = (char)std::fgetc( in );
        while( currkey.compare( "%%header_end" ) != 0 ) {
            if( c == '%' ) {
                c = (char)std::fgetc( in );
                if( c == '%' ) {
                    std::fseek( in, -6, SEEK_CUR );
                    currkey = get_next_string( in );
                    if( currkey.compare( key ) == 0 ) {
                        return true;
                    }
                }
            }
            c = (char)std::fgetc( in );
        }
        return false;
    }

    void initialize_stream( channels::data_type_t positionTypeHint ) {
        if( !m_fin )
            throw std::runtime_error( "ptg_particle_istream: Failed to open file \"" +
                                      frantic::strings::to_string( m_filename ) + "\" for reading." );

        // check to make sure this is a PTG file
        char tag[4];
        if( 4 != std::fread( tag, sizeof( char ), 4, m_fin ) )
            throw std::runtime_error( "ptg_particle_istream: Failed to read from file \"" +
                                      frantic::strings::to_string( m_filename ) + "\"." );
        if( std::strcmp( tag, "PTG" ) != 0 )
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" does not contain the PTG file identifier. Got " << std::hex << *(const int*)tag << std::dec
                << " instead.";

        // check the magic number to verify this is accually a ptg file
        boost::int32_t magic = file_read<boost::int32_t>( m_fin, m_filename );
        if( magic != PTG_MAGIC_NUMBER )
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" does not contain the PTG file Magic number. Got" << magic << " instead.";

        std::string curKey; // a string to temporarily store current key

        // check for %%header_begin
        curKey = get_next_string( m_fin );
        if( curKey.compare( "%%header_begin" ) != 0 )
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" could not find header.";

        // check version
        curKey = get_next_string( m_fin );
        if( curKey.compare( "%%version" ) == 0 ) {
            m_fileVersion = file_read<boost::int32_t>( m_fin, m_filename );
            // check if it is a supported file version
            if( m_fileVersion != 1 )
                throw invalid_particle_file_exception()
                    << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                    << "\" does not contain a supported file verion number (1). Got " << m_fileVersion << " instead.";
        } else {
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" could not find version.";
        }

        // get numCols
        if( seek_header_key( "%%cols", m_fin ) ) {
            m_numCols = file_read<boost::int32_t>( m_fin, m_filename );
        } else {
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" could not find col.";
        }
        // get numRows
        if( seek_header_key( "%%rows", m_fin ) ) {
            m_numRows = file_read<boost::int32_t>( m_fin, m_filename );
        } else {
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" could not find rows.";
        }
        // get Transform
        if( seek_header_key( "%%transform", m_fin ) ) {
            for( int i = 0; i < 4; ++i ) {
                for( int j = 0; j < 4; ++j ) {
                    m_transform[i * 4 + j] = file_read<double>( m_fin, m_filename );
                }
            }
            m_hasTransform = true;
        } else {
            m_hasTransform = false;
        }
        // get properties
        if( seek_header_key( "%%properties", m_fin ) ) {
            m_properties = file_read<boost::int32_t>( m_fin, m_filename );
        } else {
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" could not find properties.";
        }

        // Seek to the end of the header and determine the size of what has been read.
        // That amount, plus m_numCols * 8 bytes, plus m_numCols * ceil( m_numRows / 8 ) bytes, is the amount of
        // non-particle data in the file.
        seek_header_key( "%%header_end", m_fin );
        boost::int64_t nonParticleBytes = ftell( m_fin.get() );
        nonParticleBytes += ( 8 * m_numCols ) + ( m_numCols * ( ( m_numRows + 7 ) / 8 ) );

        m_finOpen = false;
        m_fin.close();

        // Create a particle channel map matching the ptg on-disk data structure.
        // At the same time, accumulate the on-disk particle size.
        m_onDiskParticleChannelMap.reset();
        boost::int64_t onDiskParticleSize = 0;

        if( m_properties &
            PTG_POSITION_AS_FLOAT ) { // 0x1: if set then xyz is in float[3] format (mutually exclusive with 0x2)
            m_onDiskParticleChannelMap.define_channel( _T("Position"), 3, channels::data_type_float32 );
            onDiskParticleSize += 12;
        } else if( m_properties & PTG_POSITION_AS_DOUBLE ) { // 0x2: if set then xyz is in double[3] format (mutually
                                                             // exclusive with 0x1)
            // TODO support loading double precision instead of converting to a float.
            m_onDiskParticleChannelMap.define_channel(
                _T("Position"), 3,
                ( positionTypeHint == channels::data_type_invalid ) ? channels::data_type_float32 : positionTypeHint );
            onDiskParticleSize += 24;
        } else {
            throw invalid_particle_file_exception()
                << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                << "\" properties has no coordinate type.";
        }
        if( m_properties & PTG_INTENSITY ) { // 0x4: if set then intensity is available as a float (range [0-1])
            m_onDiskParticleChannelMap.define_channel( _T("Intensity"), 1, channels::data_type_float16 );
            onDiskParticleSize += 4;
        }
        if( m_properties &
            PTG_COLOR ) { // 0x8: if set then RGB is available as unsigned char[3] (each with range [0-255])
            m_onDiskParticleChannelMap.define_channel( _T("Color"), 3, channels::data_type_float16 );
            onDiskParticleSize += 3;
        }
        if( m_hasTransform ) {
            m_onDiskParticleChannelMap.define_channel( _T("ScannerIndex"), 1, channels::data_type_uint32 );
        }

        m_onDiskParticleChannelMap.end_channel_definition( 1, true, true );

        // Guess the file's particle count from the number of bytes that are particles.
        boost::int64_t fileSize = frantic::files::file_size( frantic::strings::to_string( m_filename ) );
        m_particleCountGuess = ( fileSize - nonParticleBytes ) / onDiskParticleSize;

        // Add channel accessors for each channel
        if( ( m_properties & PTG_POSITION_AS_FLOAT ) || ( m_properties & PTG_POSITION_AS_DOUBLE ) ) {
            m_posAccessor = m_onDiskParticleChannelMap.get_cvt_accessor<frantic::graphics::vector3fd>( _T("Position") );
        }
        if( m_properties & PTG_INTENSITY ) {
            m_intensityAccessor = m_onDiskParticleChannelMap.get_cvt_accessor<float>( _T("Intensity") );
        }
        if( m_properties & PTG_COLOR ) {
            m_colorAccessor = m_onDiskParticleChannelMap.get_cvt_accessor<frantic::graphics::vector3f>( _T("Color") );
        }
        if( m_hasTransform ) {
            m_scannerIndexAccessor = m_onDiskParticleChannelMap.get_accessor<boost::uint32_t>( _T("ScannerIndex") );
        }
    }

    void setup_metadata() {
        frantic::channels::channel_map newMap;
        prt::add_coordinate_system( newMap );
        prt::length_unit_in_micrometers::add_channel( newMap );
        newMap.end_channel_definition();

        frantic::channels::property_map generalMetadata;

        generalMetadata.set_channel_map( newMap );

        prt::set_coordinate_system( generalMetadata, frantic::graphics::coordinate_system::right_handed_zup );
        prt::length_unit_in_micrometers::set_value( generalMetadata, 1e6 );

        if( m_hasTransform ) {
            std::vector<frantic::graphics::transform4fd> transforms;
            transforms.push_back( m_transform );
            prt::set_scanner_transforms( generalMetadata, transforms );
        }

        m_metadata.set_general_metadata( generalMetadata );
    }

  public:
    ptg_particle_istream( const frantic::tstring& file,
                          channels::data_type_t positionTypeHint = channels::data_type_invalid )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) ) {
        initialize_stream( positionTypeHint );
        m_currentParticleIndex = -1;
        m_curCol = -1;
        m_curRowPoint = -1;
        setup_metadata();
        set_channel_map( m_onDiskParticleChannelMap );
        m_tempParticleBuffer.resize( m_onDiskParticleChannelMap.structure_size() );
    }

    ptg_particle_istream( const frantic::tstring& file, const frantic::channels::channel_map& particleChannelMap,
                          channels::data_type_t positionTypeHint = channels::data_type_invalid )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) ) {
        if( positionTypeHint == channels::data_type_invalid && particleChannelMap.has_channel( _T("Position") ) ) {
            size_t positionArity = 0;
            particleChannelMap.get_channel_definition( _T("Position"), positionTypeHint, positionArity );
        }
        initialize_stream( positionTypeHint );
        m_currentParticleIndex = -1;
        m_curCol = -1;
        m_curRowPoint = -1;
        setup_metadata();
        set_channel_map( particleChannelMap );
        m_tempParticleBuffer.resize( m_onDiskParticleChannelMap.structure_size() );
    }

    virtual ~ptg_particle_istream() { close(); }

    void close() {
        m_fin.close();
        m_particleCountGuess = 0;
    }

    /**
     * Get the metadata field for the particle stream
     */
    const frantic::particles::particle_file_metadata& get_metadata() const { return m_metadata; }

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_filename; }

    boost::int64_t particle_count() const { return -1; }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return -1; }

    boost::int64_t particle_progress_count() const { return m_numCols * m_numRows; }

    boost::int64_t particle_progress_index() const { return m_numRows * m_curCol + m_curRowPoint; }

    boost::int64_t particle_count_guess() const { return m_particleCountGuess; }

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
        if( m_numCols == 0 || m_numRows == 0 ) {
            return false;
        }
        if( !m_finOpen ) {
            m_finOpen = true;
            m_fin.reset( frantic::files::tfopen( m_filename.c_str(), _T("rb") ) );
            if( !m_fin )
                throw std::runtime_error( "ptg_particle_istream.get_particle: Failed to re-open the particle file \"" +
                                          frantic::strings::to_string( m_filename ) + "\"" );

            // go to the end of the header
            if( !seek_header_key( "%%header_end", m_fin ) )
                throw invalid_particle_file_exception()
                    << "ptg_particle_istream: File \"" << frantic::strings::to_string( m_filename )
                    << "\" could not find header end.";

            // get column offsets
            m_columnOffsets.clear();
            for( boost::int32_t i = 0; i < m_numCols; ++i ) {
                boost::int64_t temp = file_read<boost::int64_t>( m_fin, m_filename );
                m_columnOffsets.push_back( temp );
            }
            m_curCol = -1;
            load_next_column();
        } else if( !m_fin )
            throw std::runtime_error( "ptg_particle_istream.get_particle: Tried to read from particle file \"" +
                                      frantic::strings::to_string( m_filename ) + "\" after it was already closed." );

        if( m_curCol < m_numCols ) {
            char* bufferPtr = &m_tempParticleBuffer[0];

            if( m_properties & PTG_POSITION_AS_FLOAT ) {
                float x = file_read<float>( m_fin, m_filename );
                float y = file_read<float>( m_fin, m_filename );
                float z = file_read<float>( m_fin, m_filename );
                m_posAccessor.set( bufferPtr, frantic::graphics::vector3fd( x, y, z ) );
            } else if( m_properties & PTG_POSITION_AS_DOUBLE ) {
                double x = file_read<double>( m_fin, m_filename );
                double y = file_read<double>( m_fin, m_filename );
                double z = file_read<double>( m_fin, m_filename );
                m_posAccessor.set( bufferPtr, frantic::graphics::vector3fd( x, y, z ) );
            }
            if( m_properties & PTG_INTENSITY ) {
                m_intensityAccessor.set( bufferPtr, file_read<float>( m_fin, m_filename ) );
            }
            if( m_properties & PTG_COLOR ) {
                float x = static_cast<float>( file_read<unsigned char>( m_fin, m_filename ) ) / 255.0f;
                float y = static_cast<float>( file_read<unsigned char>( m_fin, m_filename ) ) / 255.0f;
                float z = static_cast<float>( file_read<unsigned char>( m_fin, m_filename ) ) / 255.0f;
                m_colorAccessor.set( bufferPtr, frantic::graphics::vector3f( x, y, z ) );
            }
            if( m_hasTransform ) {
                m_scannerIndexAccessor.get( bufferPtr ) = 1;
            }

            apply_transform( bufferPtr );

            m_pcmAdaptor.copy_structure( rawParticleBuffer, &m_tempParticleBuffer[0], &m_defaultParticleBuffer[0] );

            // Go to the next particle
            ++m_currentParticleIndex;
            ++m_curRowPoint;

            // Go to the next particle
            if( m_curRowPoint >= m_numValidPointsForCurCol ) {
                load_next_column();
            }

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
};

} // namespace streams
} // namespace particles
} // namespace frantic
