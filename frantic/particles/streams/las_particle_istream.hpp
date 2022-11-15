// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/files/files.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <map>

namespace frantic {
namespace particles {
namespace streams {

/**
 * LAS File format reader.
 *
 * http://www.asprs.org/a/society/committees/standards/asprs_las_format_v12.pdf
 */
class las_particle_istream : public particle_istream {
    static const boost::int32_t LAS_MAGIC_SIGNATURE = 0x4653414C;
    static const int FILE_SRC_ID_AND_PROJ_ID_LEN = 20;
    static const int FILE_CREATE_SRC_AND_DATE_LEN = 70;
    static const int WAVE_PACKET_DATA_LEN = 20;
    static const int NUM_POINT_RETURN_LEN = 120;

    frantic::tstring m_filename;
    frantic::channels::channel_map m_onDiskChannelMap, m_nativeChannelMap, m_channelMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    frantic::particles::particle_file_metadata m_metadata;
    frantic::graphics::boundbox3fd m_boundingBox;

    boost::uint32_t m_finOffset;
    files::file_ptr m_fin;
    boost::int64_t m_particleCount, m_currentParticleIndex;
    frantic::graphics::vector3fd m_scale, m_offset;

    float m_colorRange;
    float m_colorScale;
    bool m_doneColorRangeWarning;

    frantic::channels::channel_accessor<frantic::graphics::vector3> m_onDiskPosAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> m_posAccessor;
    frantic::channels::channel_accessor<boost::uint16_t> m_onDiskIntensityAccessor;
    frantic::channels::channel_cvt_accessor<float> m_intensityAccessor;
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector3f> m_onDiskColorAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_colorAccessor;

    std::vector<char> m_defaultParticleBuffer;
    std::vector<char> m_onDiskParticleBuffer;

  public:
    las_particle_istream( const frantic::tstring& file,
                          channels::data_type_t positionTypeHint = channels::data_type_invalid )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) )
        , m_currentParticleIndex( -1 )
        , m_doneColorRangeWarning( false ) {
        initialize( positionTypeHint );
        set_channel_map( get_native_channel_map() );
    }

    las_particle_istream( const frantic::tstring& file, const frantic::channels::channel_map& particleChannelMap,
                          channels::data_type_t positionTypeHint = channels::data_type_invalid )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) )
        , m_currentParticleIndex( -1 )
        , m_doneColorRangeWarning( false ) {
        if( positionTypeHint == channels::data_type_invalid && particleChannelMap.has_channel( _T("Position") ) ) {
            size_t positionArity = 0;
            particleChannelMap.get_channel_definition( _T("Position"), positionTypeHint, positionArity );
        }
        initialize( positionTypeHint );
        set_channel_map( particleChannelMap );
    }

    virtual ~las_particle_istream() { close(); }

    void close() {
        m_fin.close();
        m_particleCount = 0;
    }

    std::size_t particle_size() const { return m_channelMap.structure_size(); }

    frantic::tstring name() const { return m_filename; }

    boost::int64_t particle_count() const { return m_particleCount; }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return m_particleCount - m_currentParticleIndex - 1; }

    boost::int64_t particle_progress_count() const { return m_particleCount; }

    boost::int64_t particle_progress_index() const { return m_currentParticleIndex; }

    const frantic::channels::channel_map& get_channel_map() const { return m_channelMap; }

    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeChannelMap; }

    /**
     * Get the metadata field for the particle stream
     */
    const frantic::particles::particle_file_metadata get_metadata() const { return m_metadata; }

    void set_channel_map( const frantic::channels::channel_map& newChannelMap ) {
        std::vector<char> newDefaultParticle( newChannelMap.structure_size() );

        if( m_defaultParticleBuffer.size() > 0 ) {
            frantic::channels::channel_map_adaptor defaultAdaptor( newChannelMap, m_channelMap );
            defaultAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticleBuffer[0] );
        } else {
            memset( &newDefaultParticle[0], 0, newChannelMap.structure_size() );
        }
        m_defaultParticleBuffer.swap( newDefaultParticle );

        m_channelMap = newChannelMap;
        m_pcmAdaptor.set( m_channelMap, m_onDiskChannelMap );

        setup_channel_accessors();
    }

    void set_default_particle( char* buffer ) { m_channelMap.copy_structure( &m_defaultParticleBuffer[0], buffer ); }

    bool get_particle( char* rawParticleBuffer ) {
        bool gotParticle = load_disk_particle( m_onDiskParticleBuffer );
        if( gotParticle ) {
            // Convert to the publicParticleChannelMap's data formats
            m_pcmAdaptor.copy_structure( rawParticleBuffer, &m_onDiskParticleBuffer[0], &m_defaultParticleBuffer[0] );
            // Scale and offset the particle according to the las file's scale and offset values
            convert_and_copy_position( rawParticleBuffer, &m_onDiskParticleBuffer[0] );
            convert_particle_intensity_to_float_scale( rawParticleBuffer, &m_onDiskParticleBuffer[0] );
            convert_particle_color_to_float_scale( rawParticleBuffer, &m_onDiskParticleBuffer[0] );
        }
        return gotParticle;
    }

    bool get_particles( char* particleBuffer, std::size_t& numParticles ) {
        std::size_t particleSize = m_channelMap.structure_size();
        for( std::size_t i = 0; i < numParticles; ++i ) {
            if( !get_particle( particleBuffer + i * particleSize ) ) {
                numParticles = i;
                return false;
            }
        }

        return true;
    }

  private:
    template <typename Type>
    static Type file_read( FILE* in, const frantic::tstring& filename ) {
        Type t;
        if( 1 != std::fread( &t, sizeof( Type ), 1, in ) )
            throw std::runtime_error( "las_particle_istream.file_read: Failed to read from file \"" +
                                      frantic::strings::to_string( filename ) + "\"." );

        return t;
    }

    void initialize( channels::data_type_t positionTypeHint ) {
        initialize_stream( positionTypeHint );
        initialize_metadata();
    }

    void initialize_metadata() {
        frantic::channels::property_map positionMetadata;
        prt::add_channel_extents( positionMetadata, m_nativeChannelMap[_T("Position")].data_type() );
        prt::set_extents( positionMetadata, m_boundingBox );

        m_metadata.set_channel_metadata( _T("Position"), positionMetadata );
    }

    void initialize_stream( channels::data_type_t positionTypeHint ) {
        const boost::uint8_t SUPPORTED_VERSION_MAJOR = 1;
        const boost::uint8_t SUPPORTED_VERSION_MINOR = 4;
        const boost::uint8_t SUPPORTED_POINT_FORMAT = 10;

        frantic::graphics::vector3fd minPos, maxPos;

        if( !m_fin )
            throw std::runtime_error( "las_particle_istream.initialize_stream: Failed to open file \"" +
                                      frantic::strings::to_string( m_filename ) + "\" for reading." );

        // Read file signature and check its actually .las
        boost::int32_t sig = file_read<boost::int32_t>( m_fin, m_filename );
        if( sig != LAS_MAGIC_SIGNATURE ) {
            throw invalid_particle_file_exception()
                << "las_particle_istream.initialize_stream: File \"" << frantic::strings::to_string( m_filename )
                << "\" does not contain the ASPRS .las file identifier. Got " << sig << " instead.";
        }

        // Skip File Source Id, Global Encoding, and ProjectID fields(20 bytes) then get the file version major and
        // minor numbers
        std::fseek( m_fin, FILE_SRC_ID_AND_PROJ_ID_LEN, SEEK_CUR );

        const boost::uint8_t versionMajor = file_read<boost::uint8_t>( m_fin, m_filename );
        const boost::uint8_t versionMinor = file_read<boost::uint8_t>( m_fin, m_filename );

        if( versionMajor > 1 || versionMinor > 4 ) {
            throw std::runtime_error(
                "las_particle_istream.initialize_stream: Failed to initialize istream because LAS file version " +
                boost::lexical_cast<std::string>( versionMajor ) + "." +
                boost::lexical_cast<std::string>( versionMinor ) + " was encountered. Highest supported version is " +
                boost::lexical_cast<std::string>( SUPPORTED_VERSION_MAJOR ) + "." +
                boost::lexical_cast<std::string>( SUPPORTED_VERSION_MINOR ) + "." );
        }

        // Skip the System Identification, Generating Software, Header Size, File Creation Day and year fields(68 bytes)
        // then get header size, offset from beginning of file to
        // point data,
        std::fseek( m_fin, FILE_CREATE_SRC_AND_DATE_LEN, SEEK_CUR );
        const boost::uint32_t offsetToData = file_read<boost::uint32_t>( m_fin, m_filename );

        // Skip the number of variable length records field
        std::fseek( m_fin, sizeof( boost::int32_t ), SEEK_CUR );
        boost::uint8_t pointDataFormat = file_read<boost::uint8_t>( m_fin, m_filename );
        boost::uint16_t pointDataLen = file_read<boost::uint16_t>( m_fin, m_filename );

        if( pointDataFormat > SUPPORTED_POINT_FORMAT ) {
            throw std::runtime_error(
                "las_particle_istream.initialize_stream: Failed to initialize istream because point record format " +
                boost::lexical_cast<std::string>( pointDataFormat ) +
                " was encountered. Highest supported version is " +
                boost::lexical_cast<std::string>( SUPPORTED_POINT_FORMAT ) + "." );
        }

        if( versionMajor == 1 && versionMinor < 4 ) {
            // Gets the number of point records in < 1.4 file format version
            m_particleCount = static_cast<long>( file_read<boost::int32_t>( m_fin, m_filename ) );

            // Skip the legacy Number of points by return field
            std::fseek( m_fin, 5 * sizeof( boost::int32_t ), SEEK_CUR );
        } else {
            // Skip both legacy fields
            std::fseek( m_fin, 6 * sizeof( boost::int32_t ), SEEK_CUR );
        }

        m_scale.x = file_read<double>( m_fin, m_filename );
        m_scale.y = file_read<double>( m_fin, m_filename );
        m_scale.z = file_read<double>( m_fin, m_filename );
        m_offset.x = file_read<double>( m_fin, m_filename );
        m_offset.y = file_read<double>( m_fin, m_filename );
        m_offset.z = file_read<double>( m_fin, m_filename );
        maxPos.x = file_read<double>( m_fin, m_filename );
        minPos.x = file_read<double>( m_fin, m_filename );
        maxPos.y = file_read<double>( m_fin, m_filename );
        minPos.y = file_read<double>( m_fin, m_filename );
        maxPos.z = file_read<double>( m_fin, m_filename );
        minPos.z = file_read<double>( m_fin, m_filename );

        m_boundingBox.set( minPos, maxPos );

        std::fseek( m_fin, WAVE_PACKET_DATA_LEN, SEEK_CUR );

        if( versionMajor == 1 && versionMinor >= 4 )
            m_particleCount = static_cast<long>( file_read<boost::int64_t>( m_fin, m_filename ) );
        else
            std::fseek( m_fin, sizeof( boost::int64_t ), SEEK_CUR );

        // Skip the rest of the header
        std::fseek( m_fin, NUM_POINT_RETURN_LEN, SEEK_CUR );
        m_finOffset = offsetToData;

        setup_on_disk_channel_map( pointDataFormat );

        // Note: The pointDataLen may be larger than the m_onDiskChannelMap.
        // From the LAS Specification:
        //
        //   'If the specified size is larger than implied by the point format
        //   type (e.g. 32 bytes instead of 28 bytes for type 1) the remaining
        //   bytes are user-specific "extra bytes.'
        if( pointDataLen < m_onDiskChannelMap.structure_size() ) {
            throw frantic::exception_stream()
                << "las_particle_istream.initialize_stream Error: "
                << "Point Data Record Length in header (" << pointDataLen << ") is smaller than "
                << "internal channel map (" << m_onDiskChannelMap.structure_size() << ") "
                << "while reading file \"" << frantic::strings::to_string( m_filename ) << "\".";
        }
        m_onDiskParticleBuffer.resize( pointDataLen );
        setup_native_channel_map( m_nativeChannelMap, m_onDiskChannelMap, positionTypeHint );

        initialize_color_range();

        m_fin.close();
        m_currentParticleIndex = -1;
    }

    void initialize_color_range() {
        m_colorRange = estimate_color_range();
        m_colorScale = 1 / m_colorRange;
    }

    float estimate_color_range() {
        using namespace frantic::channels;
        using namespace frantic::graphics;

        // Full standard range of Color channel values
        const float fullRange = 65535;
        // Smaller range of Color channel values produced by some non-standard writers
        const float smallRange = 255;

        // Number of particles to scan for the Color channel range
        const std::size_t sampleCount = 1000;

        if( !m_fin ) {
            throw frantic::exception_stream()
                << "las_particle_istream::estimate_color_range Error: Stream is not open for file: "
                << frantic::strings::to_string( m_filename );
        }

        if( !m_onDiskChannelMap.has_channel( _T("Color") ) ) {
            return fullRange;
        }

        channel_const_cvt_accessor<vector3f> colorAcc =
            m_onDiskChannelMap.get_const_cvt_accessor<vector3f>( _T("Color") );

        std::fseek( m_fin, m_finOffset, SEEK_SET );

        bool gotNonZeroColor = false;
        bool gotFullRangeColor = false;

        for( std::size_t i = 0; i < sampleCount && !gotFullRangeColor; ++i ) {
            bool gotParticle = load_disk_particle( m_onDiskParticleBuffer );
            if( !gotParticle ) {
                break;
            }
            const vector3f color = colorAcc( m_onDiskParticleBuffer );
            for( int dim = 0; dim < 3; ++dim ) {
                if( color[dim] != 0 ) {
                    gotNonZeroColor = true;
                }
                if( color[dim] > smallRange ) {
                    gotFullRangeColor = true;
                }
            }
        }

        if( gotNonZeroColor && !gotFullRangeColor ) {
            FF_LOG( warning ) << "las_particle_istream: In file " << m_filename << ": "
                              << "In the first " << ( particle_index() + 1 )
                              << " particles, all Color values are in the range [0," << smallRange << "].  "
                              << "Assuming non-standard Color range [0," << smallRange << "]." << std::endl;
            return smallRange;
        } else {
            return fullRange;
        }
    }

    void setup_on_disk_channel_map( int pointDataFormat ) {
        m_onDiskChannelMap.reset();
        // X, Y, Z, and Intensity are the same for all particle record formats
        m_onDiskChannelMap.define_channel( _T("Position"), 3, channels::data_type_int32 );
        m_onDiskChannelMap.define_channel( _T("Intensity"), 1, channels::data_type_uint16 );
        m_onDiskChannelMap.define_channel( _T("Returns"), 1, channels::data_type_int8 );

        if( pointDataFormat > 5 )
            m_onDiskChannelMap.define_channel( _T("Scanner"), 1, channels::data_type_int8 );

        m_onDiskChannelMap.define_channel( _T("Classification"), 1, channels::data_type_uint8 );

        if( pointDataFormat <= 5 )
            m_onDiskChannelMap.define_channel( _T("ScanAngleRank"), 1, channels::data_type_int8 );

        m_onDiskChannelMap.define_channel( _T("UserData"), 1, channels::data_type_uint8 );

        if( pointDataFormat > 5 )
            m_onDiskChannelMap.define_channel( _T("ScanAngle"), 1, channels::data_type_int16 );

        m_onDiskChannelMap.define_channel( _T("PointSourceID"), 1, channels::data_type_uint16 );

        if( pointDataFormat != 0 && pointDataFormat != 2 )
            m_onDiskChannelMap.define_channel( _T("GPSTime"), 1, channels::data_type_float64 );

        if( pointDataFormat != 0 && pointDataFormat != 1 && pointDataFormat != 4 && pointDataFormat != 6 &&
            pointDataFormat != 9 )
            m_onDiskChannelMap.define_channel( _T("Color"), 3, channels::data_type_uint16 );

        if( pointDataFormat == 10 || pointDataFormat == 8 )
            m_onDiskChannelMap.define_channel( _T("NIR"), 1, channels::data_type_int16 );

        if( pointDataFormat == 10 || pointDataFormat == 9 || pointDataFormat == 5 || pointDataFormat == 4 ) {
            m_onDiskChannelMap.define_channel( _T("WavePacketDescriptorIndex"), 1, channels::data_type_int8 );
            m_onDiskChannelMap.define_channel( _T("ByteOffsetToWaveformData"), 1, channels::data_type_int64 );
            m_onDiskChannelMap.define_channel( _T("WaveformPacketSizeinBytes"), 1, channels::data_type_int32 );
            m_onDiskChannelMap.define_channel( _T("ReturnPointWaveformLocation"), 1, channels::data_type_float32 );
            m_onDiskChannelMap.define_channel( _T("PositionT"), 3, channels::data_type_float32 );
        }

        m_onDiskChannelMap.end_channel_definition( 1, true, true );
    }

    static void setup_native_channel_map( frantic::channels::channel_map& nativeChannelMap,
                                          const frantic::channels::channel_map& onDiskChannelMap,
                                          channels::data_type_t positionTypeHint ) {
        nativeChannelMap.reset();
        nativeChannelMap.define_channel( _T("Position"), 3,
                                         positionTypeHint == channels::data_type_invalid ? channels::data_type_float32
                                                                                         : positionTypeHint );
        nativeChannelMap.define_channel( _T("Intensity"), 1, channels::data_type_float32 );
        if( onDiskChannelMap.has_channel( _T("Color") ) ) {
            nativeChannelMap.define_channel( _T("Color"), 3, channels::data_type_float32 );
        }
        nativeChannelMap.end_channel_definition();
    }

    void setup_channel_accessors() {
        const frantic::tstring position( _T("Position") );
        const frantic::tstring intensity( _T("Intensity") );
        const frantic::tstring color( _T("Color") );

        m_onDiskPosAccessor.reset();
        m_posAccessor.reset();
        if( m_channelMap.has_channel( position ) && m_onDiskChannelMap.has_channel( position ) ) {
            m_onDiskPosAccessor = m_onDiskChannelMap.get_accessor<frantic::graphics::vector3>( position );
            m_posAccessor = m_channelMap.get_cvt_accessor<frantic::graphics::vector3fd>( position );
        }

        m_onDiskIntensityAccessor.reset();
        m_intensityAccessor.reset();
        if( m_channelMap.has_channel( intensity ) && m_onDiskChannelMap.has_channel( intensity ) ) {
            m_onDiskIntensityAccessor = m_onDiskChannelMap.get_accessor<boost::uint16_t>( intensity );
            m_intensityAccessor = m_channelMap.get_cvt_accessor<float>( intensity );
        }

        m_onDiskColorAccessor.reset();
        m_colorAccessor.reset();
        if( m_channelMap.has_channel( color ) && m_onDiskChannelMap.has_channel( color ) ) {
            m_onDiskColorAccessor = m_onDiskChannelMap.get_const_cvt_accessor<frantic::graphics::vector3f>( color );
            m_colorAccessor = m_channelMap.get_cvt_accessor<frantic::graphics::vector3f>( color );
        }
    }

    bool load_disk_particle( std::vector<char>& diskParticleBuffer ) {
        if( !m_fin ) {
            m_fin.reset( frantic::files::tfopen( m_filename.c_str(), _T("rb") ) );
            if( !m_fin )
                throw std::runtime_error( "las_particle_istream.get_particle: Failed to re-open the particle file \"" +
                                          frantic::strings::to_string( m_filename ) + "\"" );
            std::fseek( m_fin, m_finOffset, SEEK_SET );
        }

        if( particle_index() + 1 < particle_count() ) {
            // Read particle into the public buffer, then use the adaptor to copy it into an intermediate buffer then to
            // an output buffer
            if( 1 != std::fread( &diskParticleBuffer[0], diskParticleBuffer.size(), 1, m_fin ) )
                throw std::runtime_error( "las_particle_istream.get_particle: Error reading particle " +
                                          boost::lexical_cast<std::string>( m_currentParticleIndex + 1 ) + " of " +
                                          boost::lexical_cast<std::string>( m_particleCount ) + " from the file \"" +
                                          frantic::strings::to_string( m_filename ) + "\"" );

            ++m_currentParticleIndex;

            return true;
        } else {
            return false;
        }
    }

    /**
     * A helper function to copy the Position channel, also applying the scale and
     * offset the LAS file specified.
     */
    void convert_and_copy_position( char* outBuffer, const char* onDiskBuffer ) {
        using frantic::graphics::vector3fd;

        vector3fd pos( m_onDiskPosAccessor.get( onDiskBuffer ) );
        m_posAccessor.set( outBuffer, vector3fd::component_multiply( pos, m_scale ) + m_offset );
    }

    void convert_particle_intensity_to_float_scale( char* outBuffer, const char* onDiskBuffer ) {
        const float INTENSITY_MAX = 65535.f;
        const float scale = 1.f / INTENSITY_MAX;

        if( m_intensityAccessor.is_valid() ) {
            const float intensity = m_onDiskIntensityAccessor.get( onDiskBuffer );
            m_intensityAccessor.set( outBuffer, scale * intensity );
        }
    }

    /**
     * Transforms each of a particle's color components from a color on scale 0-m_colorRange to a color on the scale
     * 0.0-1.0.
     */
    void convert_particle_color_to_float_scale( char* outBuffer, const char* onDiskBuffer ) {
        if( m_colorAccessor.is_valid() ) {
            const frantic::graphics::vector3f currParticleColor = m_onDiskColorAccessor.get( onDiskBuffer );

            for( int dim = 0; dim < 3; ++dim ) {
                if( currParticleColor[dim] > m_colorRange && !m_doneColorRangeWarning ) {
                    FF_LOG( warning )
                        << "las_particle_istream: In file " << m_filename << ": "
                        << "Found Color (" << currParticleColor << ") "
                        << "with value outside of assumed range [0," << (int)m_colorRange << "].  "
                        << "This file's colors will be incorrect.  Please contact support@thinkboxsoftware.com."
                        << std::endl;
                    m_doneColorRangeWarning = true;
                }
            }

            m_colorAccessor.set( outBuffer, m_colorScale * currParticleColor );
        }
    }

    // Private copy constructor to disable copying
    las_particle_istream( const las_particle_istream& ); // not implemented

    // Private assignment operator to disable assignment
    las_particle_istream& operator=( const las_particle_istream& ); // not implemented
};
} // namespace streams
} // namespace particles
} // namespace frantic
