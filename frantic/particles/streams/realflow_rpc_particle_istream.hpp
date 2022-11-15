// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/files/compression_stream.hpp>
#include <frantic/files/files.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

namespace detail {

class rpc_channel {
  private:
    std::ifstream m_fin;
    frantic::tstring m_name;
    std::size_t m_offset;
    std::size_t m_compressedSize, m_uncompressedSize;
    std::size_t m_arity;
    std::size_t m_rpcType;
    std::size_t m_dataSize;
    boost::uint64_t m_byteStart;
    frantic::channels::data_type_t m_channelType;
    frantic::files::zlib_inflate_istream m_inflateStream;

  public:
    ~rpc_channel() {
        m_inflateStream.close();
        m_fin.close();
    };

    rpc_channel( const frantic::tstring& filename, const frantic::tstring& name, const std::size_t offset,
                 const std::size_t compressedSize, const std::size_t arity,
                 const frantic::channels::data_type_t dataType, const std::size_t rpcType )
        : m_name( name )
        , m_offset( offset )
        , m_fin( filename.c_str(), std::ios::binary | std::ios::in )
        , m_compressedSize( compressedSize )
        , m_uncompressedSize( 0 )
        , m_arity( arity )
        , m_channelType( dataType )
        , m_rpcType( rpcType )
        , m_byteStart( 0 )
        , m_dataSize( frantic::channels::sizeof_channel_data_type( dataType ) ) {
        m_fin.seekg( m_offset, std::ios::beg );
        m_inflateStream.open( m_fin, m_name );
    }

    const frantic::tstring& get_name() const { return m_name; }

    const std::size_t get_arity() const { return m_arity; }

    const frantic::channels::data_type_t get_channel_type() const { return m_channelType; }

    const std::size_t get_data_size() const { return m_dataSize; }

    void close() {
        m_inflateStream.close();
        m_fin.close();
    }

    /**
     * Reads in a portion of the compressed channel, uncompresses it and plates the uncompressed data in the
     * uncompressedBuffer.
     */
    void uncompress_channel_data( std::vector<char>& uncompressedBuffer ) {
        m_inflateStream.read( &uncompressedBuffer[0], static_cast<int>( uncompressedBuffer.size() ) );
    }

    rpc_channel( const rpc_channel& ); // not implemented
};

} // namespace detail

class realflow_rpc_particle_istream : public particle_istream {
  private:
    static const boost::int32_t RPC_FILE_SIGNATURE = 0x70FABADA;
    static const boost::uint64_t CHUNK_SIZE = 5000; // The number of particles' data that will be cached

    bool m_particleRead;
    files::file_ptr m_fin;
    frantic::tstring m_filename;
    std::size_t m_channelCount, m_particleLocationInChunk;
    boost::uint32_t m_version, m_particleCount;
    frantic::graphics::boundbox3f m_boundingBox;
    std::size_t m_particleChunkStart;
    long m_particleIndex;

    frantic::channels::channel_map m_onDiskParticleChannelMap, m_particleChannelMap, m_nativeChannelMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_finalTextureCoordAccessor;
    frantic::channels::channel_accessor<frantic::graphics2d::vector2f> m_readTextureCoordAccessorVec2f;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_readTextureCoordAccessorVec3f;
    frantic::channels::channel_accessor<frantic::graphics::vector4f> m_readTextureCoordAccessorVec4f;
    bool m_convertTextureChannels;
    std::vector<boost::shared_ptr<detail::rpc_channel>> m_channels;
    boost::shared_ptr<detail::rpc_channel> m_textureChannel;

    frantic::particles::particle_file_metadata m_metadata;

    std::vector<char> m_tempParticleBuffer;
    std::vector<char> m_defaultParticleBuffer;

    std::vector<std::vector<char>> m_cachedChannels;

  public:
    realflow_rpc_particle_istream( const frantic::tstring& file )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) )
        , m_particleIndex( -1 )
        , m_particleRead( false )
        , m_particleChunkStart( 0 )
        , m_particleLocationInChunk( 0 )
        , m_convertTextureChannels( false ) {
        setup_istream( m_particleChannelMap );
    }

    realflow_rpc_particle_istream( const frantic::tstring& file,
                                   const frantic::channels::channel_map& particleChannelMap )
        : m_filename( file )
        , m_fin( frantic::files::tfopen( file.c_str(), _T("rb") ) )
        , m_particleIndex( -1 )
        , m_particleRead( false )
        , m_particleChunkStart( 0 )
        , m_particleLocationInChunk( 0 )
        , m_convertTextureChannels( false ) {
        setup_istream( particleChannelMap );
    }

    virtual ~realflow_rpc_particle_istream() {
        close();

        for( std::size_t i = 0; i < m_channels.size(); ++i ) {
            m_channels[i]->close();
        }

        m_channels.resize( 0 );
    }

    /**
     * Clears the cached channels and closes the file if it is open
     */
    void close() {
        m_cachedChannels.resize( 0 );
        m_fin.close();
    }

    /**
     * Get the metadata field for the particle stream
     */
    const frantic::particles::particle_file_metadata& get_metadata() const { return m_metadata; }

    /**
     * Get the size of a single particle in bytes
     */
    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    /**
     * Get the name of the file being read
     */
    frantic::tstring name() const { return m_filename; }

    /**
     * Get the number of particles contained in this stream
     */
    boost::int64_t particle_count() const { return m_particleCount; }

    /**
     * Get the index of the next particle to be read from the stream
     */
    boost::int64_t particle_index() const { return m_particleIndex; }

    /**
     * Get the number of particles left to be read from this stream
     */
    boost::int64_t particle_count_left() const { return m_particleCount - m_particleIndex - 1; }

    boost::int64_t particle_progress_count() const { return m_particleCount; }

    boost::int64_t particle_progress_index() const { return m_particleIndex; }

    /**
     * Gets the current channel map that is being used to place particles into the buffer passed to get_particles
     */
    const frantic::channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    /**
     * Gets the channel map of the particles in the RPC file
     */
    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeChannelMap; }

    /**
     * Replaces the current channel map being used to place particles into the buffer passed to get_particles with a new
     * channel map
     */
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

        if( m_particleChannelMap.has_channel( _T( "TextureCoord" ) ) ) {
            m_finalTextureCoordAccessor =
                m_particleChannelMap.get_accessor<frantic::graphics::vector3f>( _T( "TextureCoord" ) );
            m_convertTextureChannels = m_onDiskParticleChannelMap.has_channel( _T( "Texture" ) );
        } else {
            m_convertTextureChannels = false;
        }
    }

    /**
     * Sets the contents of the default particle buffer. Used to populate channels that are not in the rpc file, but
     * have default values in Krakatoa( for example the color channel )
     */
    void set_default_particle( char* buffer ) {
        m_particleChannelMap.copy_structure( &m_defaultParticleBuffer[0], buffer );
    }

    /**
     * Gets a single particle from the cached channels of the rpc file and places it into the rawParticleBuffer such
     * that the order of the data in the rawParticleBuffer matches the channel map specified by the m_particleChannelMap
     * channel map
     */
    bool get_particle( char* rawParticleBuffer ) {
        // Initialize the cached channel values here so that the first chunk of values aren't all 0's
        if( m_particleIndex == -1 )
            cache_channel_data();

        if( !m_particleRead ) {
            m_particleRead = true;
            m_fin.reset( frantic::files::tfopen( m_filename.c_str(), _T("rb") ) );
            if( !m_fin )
                throw std::runtime_error(
                    "realflow_rpc_particle_istream.get_particle: Failed to re-open the particle file \"" +
                    frantic::strings::to_string( m_filename ) + "\"" );
        } else if( !m_fin )
            throw std::runtime_error(
                "realflow_rpc_particle_istream.get_particle: Tried to read from particle file \"" +
                frantic::strings::to_string( m_filename ) + "\" after it was already closed." );

        if( particle_index() + 1 < particle_count() ) {
            char* bufferPtr = &m_tempParticleBuffer[0];
            memset( bufferPtr, 0, m_onDiskParticleChannelMap.structure_size() );

            // Go through each channel and pull the particle's channel data out, and place it in the temporary buffer
            for( std::size_t i = 0; i < m_channelCount; ++i ) {
                // If the size of the cachedChannel is greater than zero copy the current particles data into the buffer
                if( m_cachedChannels[i].size() > 0 ) {
                    // Get a pointer to the particle's data location within the cached channel and copy it to the buffer
                    char* channelPtr = &m_cachedChannels[i][0] + m_particleLocationInChunk *
                                                                     m_channels[i]->get_data_size() *
                                                                     m_channels[i]->get_arity();
                    memcpy( bufferPtr, channelPtr, m_channels[i]->get_data_size() * m_channels[i]->get_arity() );
                }

                // Update the buffer ptr
                bufferPtr += m_channels[i]->get_data_size() * m_channels[i]->get_arity();
            }

            m_pcmAdaptor.copy_structure( rawParticleBuffer, &m_tempParticleBuffer[0], &m_defaultParticleBuffer[0] );

            if( m_convertTextureChannels ) {
                copy_texture_coord_channel( rawParticleBuffer, &m_tempParticleBuffer[0] );
            }

            // Update indices of the particle to be read
            ++m_particleLocationInChunk;
            ++m_particleIndex;

            // Check to see if we need to cache a new chunk
            if( m_particleLocationInChunk >= CHUNK_SIZE ) {
                if( CHUNK_SIZE <= static_cast<std::size_t>( m_particleCount ) - m_particleChunkStart ) {
                    m_particleChunkStart += CHUNK_SIZE;
                    cache_channel_data();
                }
            }

            return true;
        } else {
            return false;
        }
    }

    /**
     * Gets particles equal to numParticles from the file/cached data and places them in the particle buffer
     */
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

  private:
    template <typename Type>
    static Type file_read( FILE* in, const frantic::tstring& filename ) {
        Type t;
        if( 1 != std::fread( &t, sizeof( Type ), 1, in ) )
            throw std::runtime_error( "realflow_rpc_particle_istream.file_read: Failed to read from file \"" +
                                      frantic::strings::to_string( filename ) + "\"." );

        return t;
    }

    /**
     * Sets up the rpc file istream
     */
    void setup_istream( const frantic::channels::channel_map& particleChannelMap ) {
        initialize_stream();
        setup_metadata();
        set_channel_map( particleChannelMap );
        m_tempParticleBuffer.resize( m_onDiskParticleChannelMap.structure_size() );
    }

    /**
     * Sets up the rpc metadata
     */
    void setup_metadata() {
        frantic::channels::channel_map newMap;
        prt::add_coordinate_system( newMap );
        newMap.end_channel_definition();

        frantic::channels::property_map generalMetadata;

        generalMetadata.set_channel_map( newMap );

        prt::set_coordinate_system( generalMetadata, frantic::graphics::coordinate_system::right_handed_yup );

        frantic::channels::property_map positionMetadata;
        prt::add_channel_extents( positionMetadata, m_onDiskParticleChannelMap[_T("Position")].data_type() );
        prt::set_extents( positionMetadata, m_boundingBox );

        m_metadata.set_general_metadata( generalMetadata );
        m_metadata.set_channel_metadata( _T("Position"), positionMetadata );
    }

    /**
     * Reads the RPC file's header data and channel map and use the values found to setup the stream and its channel map
     */
    void initialize_stream() {
        if( !m_fin )
            throw std::runtime_error( "realflow_rpc_particle_istream.initialize_stream: Failed to open file \"" +
                                      frantic::strings::to_string( m_filename ) + "\" for reading." );

        // Read file signature and check its actually .rpc
        boost::int32_t sig = 0;
        std::fread( &sig, sizeof( sig ), 1, m_fin );
        if( sig != RPC_FILE_SIGNATURE ) {
            throw invalid_particle_file_exception()
                << "realflow_rpc_particle_istream.initialize_stream: File \""
                << frantic::strings::to_string( m_filename )
                << "\" does not contain the Realflow .rpc file identifier. Got " << sig << " instead.";
        }

        // Get the header data, and on disk channel data
        read_header_data();
        read_channel_definitions();

        // Close the particle buffer
        close();
    }

    /**
     * Reads and stores all the header infomation from the file
     */
    void read_header_data() {
        if( m_fin ) {
            frantic::graphics::vector3f mins, maxes;

            m_version = file_read<boost::uint32_t>( m_fin, m_filename ); // File version should be >= 3

            if( m_version != 3 ) {
                throw std::runtime_error( "realflow_rpc_particle_istream.read_header_data: Failed to read headers "
                                          "because an invalid version type(" +
                                          boost::lexical_cast<std::string>( m_version ) + ") was found in " +
                                          frantic::strings::to_string( m_filename ) + "." );
            }

            m_particleCount = file_read<boost::uint32_t>( m_fin, m_filename );
            m_channelCount = static_cast<std::size_t>( file_read<boost::uint32_t>( m_fin, m_filename ) );

            for( int i = 0; i < 3; ++i )
                mins[i] = file_read<float>( m_fin, m_filename );

            for( int i = 0; i < 3; ++i )
                maxes[i] = file_read<float>( m_fin, m_filename );

            m_boundingBox = frantic::graphics::boundbox3f( mins, maxes );
        }
    }

    /**
     * Reads the channel definitions from the file, stores relevant channel information, and adds the channel to the
     * onDisk channel map
     */
    void read_channel_definitions() {
        if( m_fin && m_channelCount > 0 ) {
            m_onDiskParticleChannelMap.reset();
            m_particleChannelMap.reset();
            m_nativeChannelMap.reset();

            for( std::size_t i = 0; i < m_channelCount; ++i ) {
                std::string nameString;
                boost::uint32_t nameSize;
                boost::uint32_t chanType;
                long skipSize;
                std::size_t offset;
                std::size_t chanSize;

                nameSize = file_read<boost::uint32_t>( m_fin, m_filename );
                nameString.clear();
                nameString.reserve( nameSize );

                // Read the name byte-by-byte since our file_read function doesn't support reading by string
                for( boost::uint32_t j = 0; j < nameSize; ++j ) {
                    nameString.push_back( file_read<char>( m_fin, m_filename ) );
                }

                const frantic::tstring channelName = frantic::strings::to_tstring( convert_channel_name( nameString ) );

                chanType = file_read<boost::uint32_t>( m_fin, m_filename );

                offset = static_cast<std::size_t>( file_read<boost::uint64_t>( m_fin, m_filename ) );
                chanSize = static_cast<std::size_t>( file_read<boost::uint64_t>( m_fin, m_filename ) );

                // Create a channel object
                create_channel( channelName, chanType, offset, chanSize );

                m_onDiskParticleChannelMap.define_channel( m_channels[i]->get_name(), m_channels[i]->get_arity(),
                                                           m_channels[i]->get_channel_type() );

                if( channelName == _T( "Texture" ) ) {
                    m_textureChannel = m_channels[i];
                    m_particleChannelMap.define_channel( _T( "TextureCoord" ), 3,
                                                         frantic::channels::data_type_float32 );
                    m_nativeChannelMap.define_channel( _T( "TextureCoord" ), 3, frantic::channels::data_type_float32 );
                } else {
                    m_particleChannelMap.define_channel( m_channels[i]->get_name(), m_channels[i]->get_arity(),
                                                         m_channels[i]->get_channel_type() );
                    m_nativeChannelMap.define_channel( m_channels[i]->get_name(), m_channels[i]->get_arity(),
                                                       m_channels[i]->get_channel_type() );
                }

                // Skip value fields
                skipSize = static_cast<long>( ( 3 * m_channels[i]->get_arity() + 2 ) * m_channels[i]->get_data_size() );
                std::fseek( m_fin, skipSize, SEEK_CUR );
            }

            // adjust the channel count to the number of channels created rather, so that there will be no index out of
            // bounds problems
            if( m_channelCount > m_channels.size() ) {
                m_channelCount = m_channels.size();
            }

            m_onDiskParticleChannelMap.end_channel_definition();
            m_particleChannelMap.end_channel_definition();
            m_nativeChannelMap.end_channel_definition();

            cache_texture_channel_accessors();
        }
    }

    /**
     * Creates a rpc channel object that can be used to represent a channel in the rpc file
     */
    void create_channel( const frantic::tstring& name, const boost::uint32_t rpcType, const std::size_t offset,
                         const std::size_t chanSize ) {
        switch( rpcType ) {
        case 1:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 3, frantic::channels::data_type_float32, rpcType ) ) );
            return;
        case 2:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 3, frantic::channels::data_type_float64, rpcType ) ) );
            return;
        case 3:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 2, frantic::channels::data_type_float32, rpcType ) ) );
            return;
        case 4:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_int64, rpcType ) ) );
            return;
        case 5:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_uint16, rpcType ) ) );
            return;
        case 6:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_uint8, rpcType ) ) );
            return;
        case 7:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_int32, rpcType ) ) );
            return;
        case 8:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_int8, rpcType ) ) );
            return;
        case 9:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_float32, rpcType ) ) );
            return;
        case 10:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_float64, rpcType ) ) );
            return;
        case 11:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 4, frantic::channels::data_type_float32, rpcType ) ) );
            return;
        case 12:
            m_channels.push_back( boost::shared_ptr<detail::rpc_channel>( new detail::rpc_channel(
                m_filename, name, offset, chanSize, 1, frantic::channels::data_type_uint32, rpcType ) ) );
            return;
        default:
            throw std::runtime_error( "realflow_rpc_particle_istream.create_channel: Invalid channel type found(" +
                                      boost::lexical_cast<std::string>( rpcType ) + ") cannot continue reading from " +
                                      frantic::strings::to_string( m_filename ) + "." );
        }
    }

    /**
     * Reads the channel data for each channel using the ReadInterface and caches it to be used later for reading
     * particle data
     */
    void cache_channel_data() {
        m_cachedChannels.resize( m_channelCount );

        for( std::size_t i = 0; i < m_channelCount; ++i ) {
            // If the channel is part of m_particleChannelMap, then cache a chunk of the channel otherwise set the
            // cached structure's length to 0
            if( ( m_particleChannelMap.has_channel( _T("TextureCoord") ) &&
                  m_channels[i]->get_name() == _T("Texture") ) ||
                m_particleChannelMap.has_channel( m_channels[i]->get_name() ) ) {
                populate_stream( m_cachedChannels[i], i );
            } else {
                m_cachedChannels[i].resize( 0 );
            }
        }

        if( m_particleRead ) {
            m_particleRead = false;
            m_fin.close();
        }

        m_particleLocationInChunk = 0;
    }

    /**
     * Converts the name of the channel to our name scheme (ie One word with First Letter capitalized, and the rest
     * lower case) and stores it.
     */
    std::string convert_channel_name( const std::string& channelName ) const {
        std::string converted;

        converted.reserve( channelName.size() );

        // As long the substring from the beginning of the string to the current iterator position is only character or
        // digits copy it
        for( std::string::const_iterator it = channelName.begin(); it != channelName.end(); ++it ) {
            if( ( *it >= 'a' && *it <= 'z' ) || ( *it >= 'A' && *it <= 'Z' ) || ( *it >= '0' && *it <= '9' ) ) {
                converted.push_back( *it );
            }
        }

        if( converted.size() > 0 )
            converted[0] = static_cast<char>( toupper( converted[0] ) );

        if( converted.length() >= 2 && converted[converted.length() - 2] == 'I' &&
            converted[converted.length() - 1] == 'd' ) {
            converted[converted.length() - 1] = 'D';
        }

        return converted;
    }

    /**
     * Reads a compressed channel from a file, uncompresses it and places it in a buffer for use
     */
    void populate_stream( std::vector<char>& channelData, const std::size_t channelIndex ) {
        unsigned int bufferSize = CHUNK_SIZE;
        unsigned int particlesLeft =
            static_cast<unsigned int>( particle_count_left() ); // Particle count left should never be negative
        // Resizes the channel data buffer so that they can be used for input and output to the zlib functions

        if( particlesLeft > 0 ) {
            if( bufferSize > particlesLeft ) {
                bufferSize = particlesLeft; // Ensures that last particle will be placed into stream
            }

            // Multiply be the arity and data size to get the total size the buffer needs to be in bytes
            bufferSize *= static_cast<unsigned int>( m_channels[channelIndex]->get_arity() *
                                                     m_channels[channelIndex]->get_data_size() );

            channelData.resize( bufferSize );

            m_channels[channelIndex]->uncompress_channel_data( channelData );
        }
    }

    /**
     * Copies the values found in the texture channel in the file data buffer to the TextureCoord channel in the final
     * particle data by either extending the size of the of the vector found or truncating it.
     */
    void copy_texture_coord_channel( char* finalParticleData, const char* readParticleData ) const {
        if( m_textureChannel->get_channel_type() != frantic::channels::data_type_float32 ) {
            throw std::runtime_error(
                "realflow_rpc_particle_istream.copy__texture_coord_channel: Failed to copy texture "
                "coordinates from file to particle stream because " +
                ( "texture coordinate data type was type " +
                  boost::lexical_cast<std::string>( m_textureChannel->get_channel_type() ) ) +
                " instead of float." );
        }

        frantic::graphics::vector3f finalTextureCoord;

        // Convert data from file to the format necessary for Krakatoa
        if( m_textureChannel->get_arity() == 2 ) {
            frantic::graphics2d::vector2f readTextureCoord = m_readTextureCoordAccessorVec2f.get( readParticleData );

            finalTextureCoord.x = readTextureCoord.x;
            finalTextureCoord.y = readTextureCoord.y;
            finalTextureCoord.z = 0.0f;
        } else if( m_textureChannel->get_arity() == 3 ) {
            finalTextureCoord = m_readTextureCoordAccessorVec3f.get( readParticleData );
        } else if( m_textureChannel->get_arity() == 4 ) {
            frantic::graphics::vector4f readTextureCoord = m_readTextureCoordAccessorVec4f.get( readParticleData );

            for( int i = 0; i < 3; ++i )
                finalTextureCoord[i] = readTextureCoord[i];
        }

        m_finalTextureCoordAccessor.get( finalParticleData ) = finalTextureCoord;
    }

    /**
     * Cache the texture channel accessors of the onDiskChannelMap
     */
    void cache_texture_channel_accessors() {
        if( m_onDiskParticleChannelMap.has_channel( _T( "Texture" ) ) ) {
            if( m_textureChannel->get_arity() == 2 ) {
                m_readTextureCoordAccessorVec2f =
                    m_onDiskParticleChannelMap.get_accessor<frantic::graphics2d::vector2f>( _T( "Texture" ) );
            } else if( m_textureChannel->get_arity() == 3 ) {
                m_readTextureCoordAccessorVec3f =
                    m_onDiskParticleChannelMap.get_accessor<frantic::graphics::vector3f>( _T( "Texture" ) );
            } else if( m_textureChannel->get_arity() == 4 ) {
                m_readTextureCoordAccessorVec4f =
                    m_onDiskParticleChannelMap.get_accessor<frantic::graphics::vector4f>( _T( "Texture" ) );
            }
        }
    }

    // Private copy constructor to disable copying
    realflow_rpc_particle_istream( const realflow_rpc_particle_istream& ); // not implemented

    // Private assignment operator to disable assignment
    realflow_rpc_particle_istream& operator=( const realflow_rpc_particle_istream& ); // not implemented
};

} // namespace streams
} // namespace particles
} // namespace frantic
