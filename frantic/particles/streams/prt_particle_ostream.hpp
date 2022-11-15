// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/streams/particle_ostream.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/prt_file_header.hpp>

#include <frantic/files/compression_stream.hpp>
#include <frantic/files/files.hpp>

#include <frantic/graphics/color3f.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <frantic/logging/logging_level.hpp>

#include <boost/exception/all.hpp>
#include <exception>

#ifdef _WIN32
#include <Shlwapi.h>
#pragma comment( lib, "Shlwapi.lib" )
#endif

#include <cstdio>

#undef ferror

namespace frantic {
namespace particles {
namespace streams {

using frantic::channels::channel_map_adaptor;

class prt_ostream_exception : public virtual boost::exception, public virtual std::exception {
  public:
    typedef boost::error_info<struct tag_target_file, boost::filesystem::path> errinfo_target_file;
    typedef boost::error_info<struct tag_target_file, std::string> errinfo_message;

    virtual const char* what() const throw() { return boost::diagnostic_information_what( *this ); }
};

namespace {
static const std::size_t DEFAULT_BUFFER_SIZE = static_cast<std::size_t>( INT_MAX ) + 1u;
}

class prt_particle_ostream : public particle_ostream {
    prt_file_header m_header;
    channel_map m_particleChannelMap;
    channel_map_adaptor m_pcmAdaptor;

    frantic::tstring m_name;
    frantic::tstring m_targetFile;
    frantic::tstring m_tempLocalFile;
    frantic::tstring m_tempRemoteFile;
    files::file_ptr m_fout;

    boost::filesystem::path m_tempDir; // Directory used to generate the initial temporary file.

    // files::zlib_deflate_ostream m_deflateStream;
    files::zlib_deflate_ostream_cstdio m_deflateStream;

    std::vector<char> m_tempParticleBuffer;
    boost::int64_t m_currentParticleIndex, m_expectedParticleCount;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccessor;
    frantic::graphics::boundbox3f m_boundbox;

    // Private copy constructor to disable copying
    prt_particle_ostream( const prt_particle_ostream& ); // not implemented

    // Private assignment operator to disable copying
    prt_particle_ostream& operator=( const prt_particle_ostream& ); // not implemented

    void initialize_stream( const frantic::tstring& file, const channel_map& particleChannelMapForFile,
                            int zlibCompressionLevel, std::size_t fileBufferSize, std::size_t writeBufferSize,
                            const frantic::channels::property_map* generalMetadata,
                            const std::map<frantic::tstring, frantic::channels::property_map>* channelMetadata ) {
        boost::filesystem::path filePath( file );

        m_targetFile = file;

        bool useTempfolder = true;

#ifdef _WIN32
        if( !PathIsUNCW( filePath.c_str() ) ) {
            UINT dt = GetDriveTypeW( filePath.root_path().c_str() );

            // We want to use the temporary folder for output if we are targetting a remote or removable drive.
            useTempfolder = ( dt == DRIVE_REMOTE || dt == DRIVE_REMOVABLE );

            if( frantic::logging::is_logging_debug() ) {
                frantic::logging::debug << _T("Target file: \"") << filePath.string<frantic::tstring>()
                                        << _T("\" has drive type: ");
                switch( dt ) {
                case DRIVE_UNKNOWN:
                    frantic::logging::debug << _T("DRIVE_UNKNOWN");
                    break;
                case DRIVE_NO_ROOT_DIR:
                    frantic::logging::debug << _T("DRIVE_NO_ROOT_DIR");
                    break;
                case DRIVE_REMOVABLE:
                    frantic::logging::debug << _T("DRIVE_REMOVABLE");
                    break;
                case DRIVE_FIXED:
                    frantic::logging::debug << _T("DRIVE_FIXED");
                    break;
                case DRIVE_REMOTE:
                    frantic::logging::debug << _T("DRIVE_REMOTE");
                    break;
                case DRIVE_CDROM:
                    frantic::logging::debug << _T("DRIVE_CDROM");
                    break;
                case DRIVE_RAMDISK:
                    frantic::logging::debug << _T("DRIVE_RAMDISK");
                    break;
                default:
                    frantic::logging::debug << _T("Unknown");
                    break;
                }

                frantic::logging::debug << std::endl;
            }
        }
#endif

        assert( !useTempfolder ||
                ( boost::filesystem::exists( m_tempDir ) && boost::filesystem::is_directory( m_tempDir ) ) );

        boost::filesystem::path tempLocalPath, tempRemotePath;

        try {
            tempRemotePath = boost::filesystem::unique_path(
                filePath.parent_path() / ( filePath.filename().wstring() + L".%%%%-%%%%-%%%%-%%%%.tmp" ) );
            if( useTempfolder )
                tempLocalPath = boost::filesystem::unique_path(
                    m_tempDir / ( filePath.filename().wstring() + L".%%%%-%%%%-%%%%-%%%%.tmp" ) );
            else
                tempLocalPath = tempRemotePath;
        } catch( boost::filesystem::filesystem_error& e ) {
            throw prt_ostream_exception()
                << prt_ostream_exception::errinfo_target_file( filePath ) << boost::errinfo_errno( e.code().value() )
                << prt_ostream_exception::errinfo_message( e.code().message() )
                << boost::errinfo_file_name( e.path1().string() )
                << boost::errinfo_nested_exception( boost::copy_exception( e ) );
        }

        m_tempLocalFile = tempLocalPath.string<frantic::tstring>();
        m_tempRemoteFile = tempRemotePath.string<frantic::tstring>();

        // Give this stream a name that indicates both the temporary file path and the target file path.
        m_name = m_targetFile + _T(" via temp file: ") + m_tempLocalFile;

        m_fout.reset( frantic::files::tfopen( m_tempLocalFile.c_str(), _T("wb") ) );
        if( !m_fout ) {
            std::stringstream ss;
            ss << "prt_particle_ostream::initialize_stream() Failure to open file \""
               << frantic::strings::to_string( m_tempLocalFile ) << "\" as a temporary for writing to \""
               << frantic::strings::to_string( m_targetFile ) << "\"\n";
            ss << "\tError number: " << errno << "\n";
#ifdef _WIN32
            ss << "\tOS Error number: " << _doserrno << "\n";
            char errorString[256];
            strerror_s( errorString, 256, errno );
            ss << "\tError message: " << errorString << std::endl;
#else
            ss << "\tError message: " << strerror( errno ) << std::endl;
#endif

            throw std::runtime_error( ss.str() );
        }

        if( fileBufferSize == 0 ) {
            FF_LOG( debug ) << _T("Disabled buffering on file: ") << m_tempLocalFile << std::endl;

            if( 0 != setvbuf( m_fout.get(), NULL, _IONBF, 0 ) )
                FF_LOG( warning ) << _T("Failed to disable file buffer for ") << m_tempLocalFile << std::endl;
        } else if( fileBufferSize <= static_cast<std::size_t>( INT_MAX ) ) {
            FF_LOG( debug ) << _T("Set ") << fileBufferSize << _T(" buffer on file: ") << m_tempLocalFile << std::endl;

            fileBufferSize = /*std::min(*/ std::max(
                fileBufferSize, static_cast<std::size_t>( 2u ) ) /*, static_cast<std::size_t>(INT_MAX) )*/;

            if( 0 != setvbuf( m_fout.get(), NULL, _IOFBF, fileBufferSize ) )
                FF_LOG( warning ) << _T("Failed to set file buffer for ") << m_tempLocalFile << std::endl;
        }

        FF_LOG( debug ) << _T("Set ") << writeBufferSize << _T(" write/compression buffer on file: ") << m_tempLocalFile
                        << std::endl;

        // fout << "setting particle channel maps" << std::endl;
        if( particleChannelMapForFile.has_channel( _T("Position") ) ) {
            const frantic::channels::data_type_t type = particleChannelMapForFile[_T("Position")].data_type();
            if( type != frantic::channels::data_type_float32 ) {
                std::stringstream ss;
                ss << "prt_particle_ostream::initialize_stream() Position channel must be of type "
                   << frantic::strings::to_string(
                          frantic::channels::channel_data_type_str( frantic::channels::data_type_float32 ) )
                   << " "
                      "("
                   << frantic::strings::to_string( frantic::channels::channel_data_type_str( type ) ) << " given)\n";
                throw std::runtime_error( ss.str() );
            }
            m_posAccessor = particleChannelMapForFile.get_accessor<frantic::graphics::vector3f>( _T("Position") );
        } else {
            m_posAccessor.reset();
        }

        m_header.set_channel_map( particleChannelMapForFile );
        m_header.set_particle_count( m_expectedParticleCount );

        if( generalMetadata )
            m_header.set_general_metadata( *generalMetadata );

        if( channelMetadata ) {
            for( std::map<frantic::tstring, frantic::channels::property_map>::const_iterator
                     it = channelMetadata->begin(),
                     itEnd = channelMetadata->end();
                 it != itEnd; ++it )
                m_header.set_channel_metadata( it->first, it->second );
        }

        // fout << "writing the header" << std::endl;
        m_header.write_header( m_fout, m_name );

        if( 0 != fflush( m_fout.get() ) )
            FF_LOG( warning ) << _T("Failed to flush header for ") << m_tempLocalFile << std::endl;

        // Resize our temp particle buffer, note that the header removes any padding, so it might be
        // smaller than the channel_map that was passed in for initialization.
        m_tempParticleBuffer.resize( m_header.get_channel_map().structure_size() );

        // Start the compression stream for the rest of the particles
        m_deflateStream.open( m_fout, m_name, zlibCompressionLevel, static_cast<int>( writeBufferSize ) );

        // Start the particle index one before the first valid index
        m_currentParticleIndex = -1;
    }

  public:
    prt_particle_ostream( const frantic::tstring& file, const channel_map& particleChannelMap,
                          const channel_map& particleChannelMapForFile, boost::int64_t expectedParticleCount = -1,
                          int zlibCompressionLevel = Z_DEFAULT_COMPRESSION,
                          const boost::filesystem::path& tempDir = boost::filesystem::path(),
                          std::size_t fileBufferSize = DEFAULT_BUFFER_SIZE, std::size_t writeBufferSize = ( 1u << 20 ),
                          const frantic::channels::property_map* generalMetadata = NULL,
                          const std::map<frantic::tstring, frantic::channels::property_map>* channelMetadata = NULL )
        : m_fout( NULL )
        , m_expectedParticleCount( expectedParticleCount )
        , m_tempDir( tempDir ) {
        if( m_tempDir.empty() )
            m_tempDir = boost::filesystem::temp_directory_path();

        initialize_stream( file, particleChannelMapForFile, zlibCompressionLevel, fileBufferSize, writeBufferSize,
                           generalMetadata, channelMetadata );

        // Now that the header is initialized, this function can be called to set the particle channel map
        set_channel_map( particleChannelMap );
    }

    virtual ~prt_particle_ostream() { close(); }

    /** Get the file path where we are writing. */
    const frantic::tstring& get_target_file() const { return m_targetFile; }

    // This is the particle channel map which specifies the byte layout of the particle structure.
    const channel_map& get_channel_map() const { return m_particleChannelMap; }

    // This allows you to change the particle layout that's being saved on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const channel_map& particleChannelMap ) {
        m_particleChannelMap = particleChannelMap;
        // Initialize the adaptor for converting the particle format to the one in the file
        m_pcmAdaptor.set( m_header.get_channel_map(), m_particleChannelMap );
        // Reset the temporary buffer values to zero
        memset( &m_tempParticleBuffer[0], 0, m_tempParticleBuffer.size() );

        // fout << "set the particle channel map" << std::endl;
    }

    void close() {
        m_deflateStream.close();

        if( 0 != fflush( m_fout.get() ) )
            FF_LOG( warning ) << _T("Failed to flush particle data for ") << m_tempLocalFile << std::endl;

        // Only finish the file if there is no thrown exception.
        if( !std::uncaught_exceptions() ) {
            // How could m_fout be NULL here?
            if( m_fout ) {

                // What situations would make ferror(m_fout) here? Write operations should check for an error after each
                // write.
                if( ferror( m_fout ) ) {
                    m_fout.close();

                    std::stringstream ss;
                    ss << "prt_particle_ostream::close() The output stream \"" << frantic::strings::to_string( m_name )
                       << "\" had a error while writing.\n";

                    int delResult = frantic::files::remove_file( m_tempLocalFile.c_str() );
                    if( delResult != 0 ) {
                        ss << "\n";
                        ss << "\tFailed to remove temporary file \"" << frantic::strings::to_string( m_tempLocalFile )
                           << "\"\n";
                        ss << "\tError number: " << errno << "\n";
#ifdef _WIN32
                        ss << "\tOS Error number: " << _doserrno << "\n";
                        char errorString[256];
                        strerror_s( errorString, 256, errno );
                        ss << "\tError message: " << errorString << std::endl;
#else
                        ss << "\tError message: " << strerror( errno ) << std::endl;
#endif
                    }

                    throw std::runtime_error( ss.str() );
                }

                // Only rewrite the header of the file if an expected count was not supplied at the beginning
                if( m_expectedParticleCount < 0 ) {
                    m_header.rewrite_particle_count( m_fout, m_name, m_currentParticleIndex + 1 );
                } else if( m_expectedParticleCount != ( m_currentParticleIndex + 1 ) ) {
                    m_fout.close();

                    std::stringstream ss;
                    ss << "prt_particle_ostream::close() The output stream \"" << frantic::strings::to_string( m_name )
                       << "\" did not receive the expected number of particles.\n";
                    ss << "\tExpected: " << m_expectedParticleCount << "\n";
                    ss << "\tWritten: " << ( m_currentParticleIndex + 1 ) << "\n";

                    int delResult = frantic::files::remove_file( m_tempLocalFile.c_str() );
                    if( !delResult ) {
                        ss << "\n";
                        ss << "\tFailed to remove temporary file \"" << frantic::strings::to_string( m_tempLocalFile )
                           << "\"\n";
                        ss << "\tError number: " << errno << "\n";
#ifdef _WIN32
                        ss << "\tOS Error number: " << _doserrno << "\n";
                        char errorString[256];
                        strerror_s( errorString, 256, errno );
                        ss << "\tError message: " << errorString << std::endl;
#else
                        ss << "\tError message: " << strerror( errno ) << std::endl;
#endif
                    }

                    throw std::runtime_error( ss.str() );
                }

                if( m_posAccessor.is_valid() ) {
                    m_header.rewrite_particle_boundbox( m_fout, m_name, m_boundbox );
                }

                int closeResult = m_fout.close();
                if( closeResult != 0 ) {
                    std::stringstream ss;
                    ss << "prt_particle_ostream::close() The temporary file \""
                       << frantic::strings::to_string( m_tempLocalFile ) << "\" could not be closed\n";
                    ss << "\tError number: " << errno << "\n";
#ifdef _WIN32
                    ss << "\tOS Error number: " << _doserrno << "\n";
                    char errorString[256];
                    strerror_s( errorString, 256, errno );
                    ss << "\tError message: " << errorString << std::endl;
#else
                    ss << "\tError message: " << strerror( errno ) << std::endl;
#endif
                    throw std::runtime_error( ss.str() );
                }

                // We have to do this carefully:
                //  1. Move the temporary local file to the temporary remote destination. This part might fail, so we
                //  don't want it to clobber the local temp file, and also we
                //      don't want it to clobber any existing file at the target destination.
                //  2. Move the remote local file to the target remote destination. This should be a simple rename, and
                //  hence atomic.

                // TODO: Store the path objects instead of converting here.
                boost::filesystem::path tempLocalFile( m_tempLocalFile ), tempRemoteFile( m_tempRemoteFile ),
                    targetFile( m_targetFile );

                boost::system::error_code errcode;

                boost::filesystem::rename( tempLocalFile, tempRemoteFile, errcode );
                if( errcode ) {
                    if( errcode != boost::system::errc::cross_device_link )
                        BOOST_THROW_EXCEPTION( boost::filesystem::filesystem_error(
                            "PRT ostream failed to rename local temp file to remote location", tempLocalFile,
                            tempRemoteFile, errcode ) );

                    // If the error is a result of a cross-device link error, we must perform a copy-remove instead. We
                    // want to fail if the file exists in the remote destination since we chosen an appropriately random
                    // name that its a problem if that file now exists.
                    boost::filesystem::copy_file( tempLocalFile, tempRemoteFile,
                                                  boost::filesystem::copy_option::fail_if_exists, errcode );
                    if( errcode )
                        BOOST_THROW_EXCEPTION( boost::filesystem::filesystem_error(
                            "PRT ostream fail to copy local temp file to remote location", tempLocalFile,
                            tempRemoteFile, errcode ) );

                    boost::filesystem::remove( tempLocalFile, errcode );
                    if( errcode )
                        FF_LOG( error ) << _T("Failed to remove local temp file: ")
                                        << tempLocalFile.string<frantic::tstring>() << _T(" after copying to ")
                                        << tempRemoteFile.string<frantic::tstring>() << std::endl;
                }

                // Delete the target file (if it exists) so that we can rename our new file to the final name it needs.
                boost::filesystem::remove( targetFile, errcode );
                if( errcode )
                    BOOST_THROW_EXCEPTION( boost::filesystem::filesystem_error(
                        "PRT ostream failed to remove existing file with target filename", targetFile, errcode ) );

                boost::filesystem::rename( tempRemoteFile, targetFile, errcode );
                if( errcode )
                    BOOST_THROW_EXCEPTION( boost::filesystem::filesystem_error(
                        "PRT ostream fail to rename remote temp file to target filename", tempRemoteFile, targetFile,
                        errcode ) );
            }
        } else {
            m_fout.close();

            if( boost::filesystem::exists( m_tempLocalFile ) ) {
                boost::filesystem::remove( m_tempLocalFile );
            }

            if( boost::filesystem::exists( m_tempRemoteFile ) ) {
                boost::filesystem::remove( m_tempRemoteFile );
            }
        }
    }

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    void put_particle( const char* rawParticleData ) {
        if( !m_fout )
            throw std::runtime_error( "prt_particle_ostream.put_particle: Tried to write to particle file \"" +
                                      frantic::strings::to_string( m_name ) + "\" after it was closed." );

        // Make sure that we haven't written more than the expected count
        if( m_expectedParticleCount > 0 && ( m_currentParticleIndex + 1 ) >= m_expectedParticleCount )
            throw std::runtime_error(
                "prt_particle_ostream.put_particle: Tried to write more particles than were specified to the file \"" +
                frantic::strings::to_string( m_name ) + "\".  The specified number of particles is " +
                boost::lexical_cast<std::string>( m_expectedParticleCount ) );

        const char* writeBuffer;
        if( m_pcmAdaptor.is_identity() ) {
            m_deflateStream.write( reinterpret_cast<const char*>( rawParticleData ), particle_size() );
            writeBuffer = rawParticleData;
        } else {
            m_pcmAdaptor.copy_structure( &m_tempParticleBuffer[0], rawParticleData );
            m_deflateStream.write( &m_tempParticleBuffer[0], m_header.get_channel_map().structure_size() );
            writeBuffer = &m_tempParticleBuffer[0];
        }
        if( m_posAccessor.is_valid() ) {
            m_boundbox += m_posAccessor( writeBuffer );
        }

        ++m_currentParticleIndex;

        if( ferror( m_fout ) != 0 )
            throw std::runtime_error( "prt_particle_ostream.put_particle: Error while writing particle to file \"" +
                                      frantic::strings::to_string( m_name ) + "\"" );
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
