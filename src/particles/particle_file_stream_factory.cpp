// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/particle_file_stream_factory.hpp>

#include <boost/predef.h>

#include <frantic/particles/csv_metadata.hpp>
#include <frantic/particles/prt_metadata.hpp>

#include <frantic/particles/streams/csv_particle_istream.hpp>
#include <frantic/particles/streams/csv_particle_ostream.hpp>
#include <frantic/particles/streams/e57_particle_istream.hpp>
#include <frantic/particles/streams/e57_particle_ostream.hpp>
#include <frantic/particles/streams/las_particle_istream.hpp>
#include <frantic/particles/streams/ply_particle_istream.hpp>
#include <frantic/particles/streams/ply_particle_ostream.hpp>
#include <frantic/particles/streams/prt2_particle_istream.hpp>
#include <frantic/particles/streams/prt2_particle_ostream.hpp>
#include <frantic/particles/streams/prt_particle_istream.hpp>
#include <frantic/particles/streams/prt_particle_ostream.hpp>
#include <frantic/particles/streams/ptg_particle_istream.hpp>
#include <frantic/particles/streams/pts_particle_istream.hpp>
#include <frantic/particles/streams/ptx_particle_istream.hpp>
#include <frantic/particles/streams/realflow_bin_particle_istream.hpp>
#include <frantic/particles/streams/realflow_bin_particle_ostream.hpp>
#include <frantic/particles/streams/realflow_rpc_particle_istream.hpp>

#include <frantic/particles/streams/channel_range_map_particle_istream.hpp>
#include <frantic/particles/streams/concatenated_particle_istream.hpp>
#include <frantic/particles/streams/set_channel_particle_istream.hpp>
#include <frantic/particles/streams/transformed_particle_istream.hpp>

#include <frantic/prtfile/sprt_common.hpp>

#include <frantic/files/aln.hpp>
#include <frantic/files/files.hpp>
#include <frantic/logging/progress_logger.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include <boost/system/error_code.hpp>

using namespace std;
using namespace boost;
using namespace frantic::particles::streams;

namespace frantic {
#ifdef _WIN32
frantic::tchar* getenv( const frantic::tchar* envVar ) {
    static frantic::tchar tempStorage[_MAX_ENV];
    if( 0 == GetEnvironmentVariable( envVar, tempStorage, _MAX_ENV ) )
        return NULL;
    return tempStorage;
}
#else
using std::getenv;
#endif
} // namespace frantic

namespace frantic {
namespace particles {

typedef std::map<frantic::tstring, prt::channel_interpretation::option> channel_metadata_type;

class particle_file_stream_factory_object::impl_type {
  public:
    /**
     * \return A path to a directory suitable for writing temporary files. This is determined by:
     *          1. A specified path via particle_file_stream_factory_object::set_temp_directory(), or
     *          2. The PRT_TEMPDIR environment variable, or
     *          3. boost::filesystem::temp_directory_path()
     */
    boost::filesystem::path temp_directory_path() const;

    /**
     * \return The size of the file buffer that particle file streams should use. This is the OS's buffer backing the
     * file. Setting it to 0
     *         will disable buffering.
     */
    std::size_t get_buffer_size() const;

    /**
     * \return The size of data buffer to write to the file stream at a time. This is the amount of data to compress and
     * write as a unit.
     */
    std::size_t get_write_size() const;

    /**
     * \return The frame rate (in frames per second) associated with the particle data. Represented as a rational
     * number.
     */
    const std::pair<unsigned, unsigned> get_frame_rate() const;

    /**
     * \return The scale value to convert length measures into meters.
     */
    double get_length_unit_in_micrometers() const;

    /**
     * \return The type of coordinate system the scene is using.
     */
    frantic::graphics::coordinate_system::option get_coordinate_system() const;

    /**
     * \return Check if a column mapping has been set
     */
    bool has_column_mapping() const;

    /**
     * Applies the appropriate transformation to convert the data in 'pin' into the coordinate system and units
     * configured in the factory. \param pin The stream to modify \param scaleToMeters The scale to apply to the source
     * in order to get meters. \param fromCoordSys The coordinate system type stored in the stream \param
     * channelMetadata An object containing transformation information for the channels affected. This determines which
     * channels are affected, and how they are transformed. \param toTransform List of transforms to apply
     * unit/coordinate transforms to \return A transformed stream if a conversion could be computed, or the original
     * stream otherwise.
     */
    particle_istream_ptr apply_coordinate_and_unit_conversion(
        particle_istream_ptr pin, double scaleToMeters, frantic::graphics::coordinate_system::option fromCoordSys,
        const channel_metadata_type& channelMetadata, std::vector<frantic::graphics::transform4fd>& toTransform ) const;

    particle_istream_ptr apply_range_mapping_conversion( particle_istream_ptr pin,
                                                         const particle_file_metadata& fileMetadata );

    particle_istream_ptr apply_factory_decorators( boost::shared_ptr<particle_istream> pin,
                                                   particle_file_metadata& outMetadata );

    /**
     * Creates a new istream out of the given file, optionally using the supplied channel map, and applying the
     * appropriate decorators to the object
     *
     * @param file The filename to construct the stream from
     * @param channelMap An optional channel map to use
     * @param outMetadata Metadata read from the stream
     * @return the particle stream
     * @throws runtime_error if no stream object could be created
     */
    particle_istream_ptr create_istream( const frantic::tstring& file,
                                         const frantic::channels::channel_map* channelMap = NULL,
                                         particle_file_metadata* outMetadata = NULL );

    /**
     * Creates a new istream out of the given file, optionally using the supplied channel map
     *
     * @param file The filename to construct the stream from
     * @param channelMap An optional channel map to use
     * @param outMetadata Metadata read from the stream
     * @return the file particle stream object
     * @throws runtime_error if no stream object could be created
     */
    boost::shared_ptr<streams::particle_istream> create_raw_istream( const frantic::tstring& file,
                                                                     const frantic::channels::channel_map* channelMap,
                                                                     particle_file_metadata& outMetadata );

  public:
    frantic::graphics::coordinate_system::option m_coordinateSystem;
    boost::filesystem::path m_tempDirectory;
    double m_scaleToMicrometers;
    std::pair<unsigned, unsigned>
        m_frameRate; // TODO: We can use changes in framerate to adjust time dependent parameters like Velocity.
    channels::data_type_t m_positionTypeHint;

    // Options for PRT2 saving
    bool m_enablePRT2Saving;
    prtfile::prt2_compression_t m_compressionScheme;
    intptr_t m_desiredChunkSizeInBytes;

    particle_file_metadata m_overrideMetadata;
};

boost::filesystem::path particle_file_stream_factory_object::impl_type::temp_directory_path() const {
    if( !m_tempDirectory.empty() )
        return m_tempDirectory;

    boost::filesystem::path result;

    frantic::tchar* pTempEnv = frantic::getenv( _T("PRT_TEMPDIR") );
    if( pTempEnv ) {
        result = pTempEnv;

        if( boost::filesystem::exists( result ) && boost::filesystem::is_directory( result ) )
            return result;

        FF_LOG( warning ) << _T("PRT_TEMPDIR specifies an invalid directory: \"") << result.string<frantic::tstring>()
                          << std::endl;
    }

    return boost::filesystem::temp_directory_path();
}

std::size_t particle_file_stream_factory_object::impl_type::get_buffer_size() const {
    frantic::tchar* pBufSizeEnv = frantic::getenv( _T("PRT_BUFFERSIZE") );
    if( pBufSizeEnv )
        return boost::lexical_cast<std::size_t>( pBufSizeEnv );

    // return BUFSIZ;
    return frantic::particles::streams::DEFAULT_BUFFER_SIZE;
}

std::size_t particle_file_stream_factory_object::impl_type::get_write_size() const {
    frantic::tchar* pBufSizeEnv = frantic::getenv( _T("PRT_WRITESIZE") );
    if( pBufSizeEnv )
        return boost::lexical_cast<std::size_t>( pBufSizeEnv );

    return ( 1u << 20 ); // 1MB write size.
}

const std::pair<unsigned, unsigned> particle_file_stream_factory_object::impl_type::get_frame_rate() const {
    return m_frameRate;
}

double particle_file_stream_factory_object::impl_type::get_length_unit_in_micrometers() const {
    return m_scaleToMicrometers;
}

frantic::graphics::coordinate_system::option
particle_file_stream_factory_object::impl_type::get_coordinate_system() const {
    return m_coordinateSystem;
}

bool particle_file_stream_factory_object::impl_type::has_column_mapping() const {
    return csv::has_column_mapping( m_overrideMetadata.get_general_metadata() );
}

void get_channel_interpretation_metadata( const frantic::channels::channel_map& channelMap,
                                          const particle_file_metadata& fileMetadata,
                                          channel_metadata_type& outChannelMetadata ) {
    // Initialize default values.
    outChannelMetadata[_T("Position")] = prt::channel_interpretation::point;
    outChannelMetadata[_T("BirthPosition")] = prt::channel_interpretation::point;
    outChannelMetadata[_T("Velocity")] = prt::channel_interpretation::vector;
    outChannelMetadata[_T("Acceleration")] = prt::channel_interpretation::vector;
    outChannelMetadata[_T("Normal")] = prt::channel_interpretation::normal;
    outChannelMetadata[_T("Tangent")] = prt::channel_interpretation::normal;
    outChannelMetadata[_T("Binormal")] = prt::channel_interpretation::normal;
    outChannelMetadata[_T("Orientation")] = prt::channel_interpretation::orientation;
    outChannelMetadata[_T("Spin")] = prt::channel_interpretation::rotation;
    outChannelMetadata[_T("Radius")] = prt::channel_interpretation::scalar;

    for( std::size_t i = 0, iEnd = channelMap.channel_count(); i < iEnd; ++i ) {
        const frantic::channels::channel& ch = channelMap[i];

        if( const frantic::channels::property_map* channelMetadata = fileMetadata.get_channel_metadata( ch.name() ) ) {
            prt::channel_interpretation::option channelType = prt::get_channel_interpretation( *channelMetadata );
            if( channelType != prt::channel_interpretation::unspecified ) {
                outChannelMetadata[ch.name()] = prt::get_channel_interpretation( *channelMetadata );
            }
        }
    }
}

particle_istream_ptr particle_file_stream_factory_object::impl_type::apply_coordinate_and_unit_conversion(
    particle_istream_ptr pin, double scaleToMeters, frantic::graphics::coordinate_system::option fromCoordSys,
    const channel_metadata_type& channelMetadata, std::vector<frantic::graphics::transform4fd>& toTransform ) const {

    // Compute the matrix that accomodates units & coordinate system changes. Order doesn't matter due to the nature of
    // respective matrices.
    frantic::graphics::transform4fd tm;
    frantic::graphics::coordinate_system::create_transform( tm, fromCoordSys, this->get_coordinate_system() );
    prt::length_unit_in_micrometers::create_transform( tm, scaleToMeters * 1.0e6,
                                                       this->get_length_unit_in_micrometers() );

    if( !tm.is_identity() && !toTransform.empty() ) {
        for( std::vector<frantic::graphics::transform4fd>::iterator transIter = toTransform.begin();
             transIter != toTransform.end(); ++transIter ) {
            *transIter = tm * ( *transIter );
        }
    }

    // Will only create a transform decorator if the tm is non-identity.
    return frantic::particles::streams::apply_transform_to_particle_istream(
        pin, tm, frantic::graphics::transform4fd::zero(), channelMetadata );
}

particle_file_stream_factory_object::particle_file_stream_factory_object()
    : m_pImpl( new impl_type ) {
    set_to_defaults();
}

particle_file_stream_factory_object::~particle_file_stream_factory_object() { delete m_pImpl; }

frantic::graphics::coordinate_system::option particle_file_stream_factory_object::get_coordinate_system() {
    return m_pImpl->m_coordinateSystem;
}

void particle_file_stream_factory_object::set_coordinate_system(
    frantic::graphics::coordinate_system::option coordinateSystem ) {
    m_pImpl->m_coordinateSystem = coordinateSystem;
}

double particle_file_stream_factory_object::get_length_unit_in_meters() const {
    return m_pImpl->get_length_unit_in_micrometers() / 1.0e6;
}

void particle_file_stream_factory_object::set_length_unit_in_meters( double scaleToMeters ) {
    m_pImpl->m_scaleToMicrometers = 1.0e6 * scaleToMeters;
}

double particle_file_stream_factory_object::get_length_unit_in_micrometers() const {
    return m_pImpl->get_length_unit_in_micrometers();
}

void particle_file_stream_factory_object::set_length_unit_in_micrometers( double scaleToMicrometers ) {
    m_pImpl->m_scaleToMicrometers = scaleToMicrometers;
}

channels::data_type_t particle_file_stream_factory_object::get_position_type_hint() const {
    return m_pImpl->m_positionTypeHint;
}

void particle_file_stream_factory_object::set_position_type_hint( channels::data_type_t positionTypeHint ) {
    m_pImpl->m_positionTypeHint = positionTypeHint;
}

std::pair<unsigned, unsigned> particle_file_stream_factory_object::get_frame_rate() const {
    return m_pImpl->get_frame_rate();
}

void particle_file_stream_factory_object::set_frame_rate( unsigned n, unsigned d ) {
    m_pImpl->m_frameRate.first = n;
    m_pImpl->m_frameRate.second = d;
}

const particle_file_metadata& particle_file_stream_factory_object::get_override_metadata() const {
    return m_pImpl->m_overrideMetadata;
}

void particle_file_stream_factory_object::set_override_metadata( const particle_file_metadata& overrideMetadata ) {
    m_pImpl->m_overrideMetadata = overrideMetadata;
}

void particle_file_stream_factory_object::enable_prt2_saving( frantic::prtfile::prt2_compression_t compressionScheme,
                                                              intptr_t desiredChunkSizeInBytes ) {
    m_pImpl->m_enablePRT2Saving = true;
    m_pImpl->m_compressionScheme = compressionScheme;
    m_pImpl->m_desiredChunkSizeInBytes = desiredChunkSizeInBytes;
}

void particle_file_stream_factory_object::set_temp_directory( const boost::filesystem::path& tempDir ) {
    if( !tempDir.empty() && ( !boost::filesystem::exists( tempDir ) || !boost::filesystem::is_directory( tempDir ) ) )
        BOOST_THROW_EXCEPTION( boost::filesystem::filesystem_error(
            "Invalid temp directory", tempDir,
            boost::system::errc::make_error_code( boost::system::errc::not_a_directory ) ) );

    m_pImpl->m_tempDirectory = tempDir;
}

streams::particle_istream_ptr
particle_file_stream_factory_object::create_istream( const frantic::tstring& file,
                                                     particle_file_metadata* outMetadata ) {
    return m_pImpl->create_istream( file, NULL, outMetadata );
}

streams::particle_istream_ptr particle_file_stream_factory_object::create_istream(
    const frantic::tstring& file, const channel_map& particleChannelMap, particle_file_metadata* outMetadata ) {
    return m_pImpl->create_istream( file, &particleChannelMap, outMetadata );
}

streams::particle_istream_ptr particle_file_stream_factory_object::impl_type::create_istream(
    const frantic::tstring& file, const channel_map* particleChannelMap, particle_file_metadata* outMetadata ) {

    particle_file_metadata localMetadata;
    const bool hasMetadata = outMetadata != NULL;

    if( !hasMetadata ) {
        outMetadata = &localMetadata;
    }

    streams::particle_istream_ptr pin = create_raw_istream( file, particleChannelMap, *outMetadata );

    outMetadata->append( m_overrideMetadata );

    return apply_factory_decorators( pin, *outMetadata );
}

particle_istream_ptr
particle_file_stream_factory_object::impl_type::apply_factory_decorators( boost::shared_ptr<particle_istream> rawStream,
                                                                          particle_file_metadata& metadata ) {

    double lengthConversion = prt::length_unit_in_meters::get_value( metadata.get_general_metadata() );
    frantic::graphics::coordinate_system::option coordinateSystem =
        prt::get_coordinate_system( metadata.get_general_metadata() );

    channel_metadata_type channelRoleMetadata;
    get_channel_interpretation_metadata( rawStream->get_native_channel_map(), metadata, channelRoleMetadata );

    std::vector<frantic::graphics::transform4fd> scannerPositions;
    frantic::particles::prt::get_scanner_transforms( metadata.get_general_metadata(), scannerPositions );

    // default values will not apply any additional decorators, so this should not incur a major performance penalty
    rawStream = apply_coordinate_and_unit_conversion( rawStream, lengthConversion, coordinateSystem,
                                                      channelRoleMetadata, scannerPositions );

    if( !scannerPositions.empty() )
        frantic::particles::prt::set_scanner_transforms( metadata.get_general_metadata(), scannerPositions );

    rawStream = apply_range_mapping_conversion( rawStream, metadata );

    return rawStream;
}

particle_istream_ptr particle_file_stream_factory_object::impl_type::apply_range_mapping_conversion(
    particle_istream_ptr pin, const particle_file_metadata& fileMetadata ) {
    static const frantic::tstring g_targetChannels[] = {
        _T( "Color" ),
        _T( "Intensity" ),
    };

    // Look through the channels for ones which we expect a normalized range, and apply them
    static const size_t g_numTargetChannels = sizeof( g_targetChannels ) / sizeof( frantic::tstring );

    for( size_t i = 0; i < g_numTargetChannels; ++i ) {
        const frantic::channels::property_map* channelProps = fileMetadata.get_channel_metadata( g_targetChannels[i] );
        if( channelProps && prt::has_channel_range( *channelProps ) ) {
            std::pair<double, double> range = prt::get_channel_range<double>( *channelProps );

            streams::channel_range_map rangeMap;
            rangeMap.channelMinimum = range.first;
            rangeMap.channelMaximum = range.second;
            rangeMap.outputMinimum = 0.0;
            rangeMap.outputMaximum = 1.0;

            pin = particle_istream_ptr(
                new streams::channel_range_map_particle_istream( pin, g_targetChannels[i], rangeMap ) );
        }
    }

    return pin;
}

particle_istream_ptr particle_file_stream_factory_object::impl_type::create_raw_istream(
    const frantic::tstring& file, const channel_map* particleChannelMap, particle_file_metadata& outMetadata ) {
    outMetadata.clear();

    std::stringstream errorStream;
    const frantic::tstring ext = frantic::strings::to_lower( frantic::files::extension_from_path( file ) );

    if( !files::file_exists( file ) )
        throw runtime_error( "particle_file_istream_factory: The input file \"" + frantic::strings::to_string( file ) +
                             "\" does not exist." );

    if( ext == _T(".sprt") && frantic::prtfile::has_sprt1_magic_number( file ) ) {
        errorStream << "The file \"" << frantic::strings::to_string( file )
                    << "\" is in the original prototype SPRT format, which is now obsolete.\n"
                    << "Please convert your file to the current SPRT format using Sequoia beta version 0.1.11.";
        throw std::runtime_error( errorStream.str() );
    }

    if( ext == _T(".prt") && prtfile::sniff_prt_spec_revision( file ) < 3 ) {
        boost::shared_ptr<prt_particle_istream> pin;

        if( particleChannelMap != NULL ) {
            pin.reset( new prt_particle_istream( file, *particleChannelMap ) );
        } else {
            pin.reset( new prt_particle_istream( file ) );
        }

        if( frantic::logging::is_logging_debug() ) {
            frantic::logging::debug << _T("PRT file \"") << file << _T("\" metadata:\n");
            pin->get_general_metadata().dump( frantic::logging::debug );

            for( std::size_t i = 0, iEnd = pin->get_channel_map().channel_count(); i < iEnd; ++i ) {
                if( const frantic::channels::property_map* chMetadata =
                        pin->get_channel_metadata( pin->get_channel_map()[i].name() ) ) {
                    frantic::logging::debug << _T("Channel \"") << pin->get_channel_map()[i].name()
                                            << _T("\" metadata:\n");
                    chMetadata->dump( frantic::logging::debug );
                }
            }

            frantic::logging::debug << std::endl;
        }

        outMetadata.set_general_metadata( pin->get_general_metadata() );

        for( size_t i = 0; i < pin->get_channel_map().channel_count(); ++i ) {
            const frantic::channels::channel& ch = pin->get_channel_map()[i];
            const frantic::channels::property_map* channelMetadata = pin->get_channel_metadata( ch.name() );
            if( channelMetadata ) {
                outMetadata.set_channel_metadata( ch.name(), *channelMetadata );
            }
        }

        return pin;
    }

    // Support for the PRT2 format
    if( ext == _T(".prt") || ext == _T(".sprt") ) {
        boost::shared_ptr<prt2_particle_istream> pin;

        if( particleChannelMap != NULL ) {
            pin.reset( new prt2_particle_istream( file, *particleChannelMap ) );
        } else {
            pin.reset( new prt2_particle_istream( file, m_positionTypeHint ) );
        }

        outMetadata.set_general_metadata( pin->get_general_metadata() );

        for( size_t i = 0; i < pin->get_channel_map().channel_count(); ++i ) {
            const frantic::channels::channel& ch = pin->get_channel_map()[i];
            const frantic::channels::property_map* channelMetadata = pin->get_channel_metadata( ch.name() );
            if( channelMetadata ) {
                outMetadata.set_channel_metadata( ch.name(), *channelMetadata );
            }
        }

        return pin;
    }

    if( ext == _T(".bin") ) {
        boost::shared_ptr<realflow_bin_particle_istream> pin;

        if( particleChannelMap != NULL ) {
            pin.reset( new realflow_bin_particle_istream( file, *particleChannelMap ) );
        } else {
            pin.reset( new realflow_bin_particle_istream( file ) );
        }

        outMetadata = pin->get_metadata();

        return pin;
    }

    if( ext == _T(".las") ) {
        boost::shared_ptr<las_particle_istream> pin;

        if( particleChannelMap != NULL ) {
            pin.reset( new las_particle_istream( file, *particleChannelMap, m_positionTypeHint ) );
        } else {
            pin.reset( new las_particle_istream( file, m_positionTypeHint ) );
        }

        outMetadata = pin->get_metadata();

        return pin;
    }

    if( ext == _T(".ptg") ) {
        boost::shared_ptr<ptg_particle_istream> pin;

        if( particleChannelMap != NULL ) {
            pin.reset( new ptg_particle_istream( file, *particleChannelMap, m_positionTypeHint ) );
        } else {
            pin.reset( new ptg_particle_istream( file, m_positionTypeHint ) );
        }

        outMetadata = pin->get_metadata();

        return pin;
    }

    if( ext == _T(".pts") ) {
        boost::shared_ptr<pts_particle_istream> pin;

        if( particleChannelMap != NULL ) {
            pin.reset( new pts_particle_istream( file, *particleChannelMap ) );
        } else {
            pin.reset( new pts_particle_istream( file, m_positionTypeHint ) );
        }

        if( particleChannelMap != NULL ) {
            pin->set_channel_map( *particleChannelMap );
        }

        return pin;
    }

    if( ext == _T(".ptx") ) {
        boost::shared_ptr<ptx_particle_istream> pin( new ptx_particle_istream( file, true, m_positionTypeHint ) );

        if( particleChannelMap != NULL ) {
            pin->set_channel_map( *particleChannelMap );
        }

        outMetadata = pin->get_metadata();

        return pin;
    }

    if( ext == _T(".rpc") ) {
        boost::shared_ptr<realflow_rpc_particle_istream> pin;

        if( particleChannelMap != NULL ) {
            pin.reset( new realflow_rpc_particle_istream( file, *particleChannelMap ) );
        } else {
            pin.reset( new realflow_rpc_particle_istream( file ) );
        }

        outMetadata = pin->get_metadata();

        return pin;
    }

#if defined( E57_AVAILABLE )
    if( ext == _T(".e57") ) {
        boost::shared_ptr<e57_particle_istream> pin;

        pin.reset( new e57_particle_istream( file, m_positionTypeHint ) );

        if( particleChannelMap != NULL )
            pin->set_channel_map( *particleChannelMap );

        outMetadata = pin->get_metadata();

        return pin;
    }
#endif

    if( ext == _T(".ply") ) {
        boost::shared_ptr<ply_particle_istream> pin( new ply_particle_istream( file ) );

        if( particleChannelMap != NULL ) {
            pin->set_channel_map( *particleChannelMap );
        }

        return pin;
    }

    if( ext == _T( ".aln" ) ) {
        boost::shared_ptr<particle_istream> pin(
            frantic::files::aln::create_aln_particle_istream( file, outMetadata, m_positionTypeHint ) );
        if( particleChannelMap != NULL ) {
            pin->set_channel_map( *particleChannelMap );
        }

        return pin;
    }

    try {
        boost::shared_ptr<csv_particle_istream> pResult;

        if( has_column_mapping() ) {
            frantic::channels::channel_column_map columnMap(
                csv::get_column_mapping( m_overrideMetadata.get_general_metadata() ) );
            pResult.reset( new csv_particle_istream(
                file, columnMap, csv::get_text_delimiter( m_overrideMetadata.get_general_metadata() ),
                (int)csv::get_header_row_count( m_overrideMetadata.get_general_metadata() ), m_positionTypeHint ) );
        } else {
            pResult.reset( new csv_particle_istream(
                file, csv::get_text_delimiter( m_overrideMetadata.get_general_metadata() ), m_positionTypeHint ) );
        }

        if( particleChannelMap != NULL ) {
            pResult->set_channel_map( *particleChannelMap );
        }

        return pResult;
    } catch( const invalid_particle_file_exception& e ) {
        errorStream << e.what() << "\n";
    }

    errorStream << "particle_istream_factory: The file \"" << frantic::strings::to_string( file )
                << "\" was not a particle file, or could not be loaded by any known particle_istream.";
    throw std::runtime_error( errorStream.str() );
}

boost::shared_ptr<streams::particle_ostream> particle_file_stream_factory_object::create_ostream(
    const frantic::tstring& file, const channel_map& particleChannelMap, const channel_map& particleChannelMapForFile,
    const particle_file_metadata* extraMetadata, boost::int64_t expectedParticleCount, int zlibCompressionLevel ) {
    frantic::tstring ext = strings::to_lower( files::extension_from_path( file ) );
    if( ext == _T( ".csv" ) ) {
        return boost::shared_ptr<streams::particle_ostream>(
            new csv_particle_ostream( file, particleChannelMap, particleChannelMapForFile, expectedParticleCount ) );
    } else if( ext == _T( ".bin" ) ) {
        return boost::shared_ptr<streams::particle_ostream>( new realflow_bin_particle_ostream(
            file, particleChannelMap, expectedParticleCount, m_pImpl->get_coordinate_system() ) );
    } else if( ext == _T( ".ply" ) ) {
        return boost::shared_ptr<streams::particle_ostream>(
            new ply_particle_ostream( file, particleChannelMap, particleChannelMapForFile ) );
#if defined( E57_AVAILABLE )
    } else if( ext == _T(".e57") ) {
        bool hasCoordSystem =
            this->m_pImpl->get_coordinate_system() != frantic::graphics::coordinate_system::unspecified;
        bool hasUnitInMeters = this->m_pImpl->get_length_unit_in_micrometers() != 0.0;

        frantic::channels::channel_map metadataMap;
        if( hasCoordSystem )
            prt::add_coordinate_system( metadataMap );
        if( hasUnitInMeters ) {
            prt::length_unit_in_micrometers::add_channel( metadataMap );
        }
        metadataMap.end_channel_definition();
        frantic::channels::property_map metadata( metadataMap );
        if( hasCoordSystem )
            prt::set_coordinate_system( metadata, this->m_pImpl->get_coordinate_system() );
        if( hasUnitInMeters ) {
            prt::length_unit_in_micrometers::set_value( metadata, this->m_pImpl->get_length_unit_in_micrometers() );
        }

        frantic::channels::channel_map channelMetadataMap;
        prt::add_channel_interpretation( channelMetadataMap );
        channelMetadataMap.end_channel_definition();

        std::map<frantic::tstring, frantic::channels::property_map> channelMetadata;

        for( std::size_t i = 0, iEnd = particleChannelMapForFile.channel_count(); i < iEnd; ++i ) {
            const frantic::channels::channel& ch = particleChannelMapForFile[i];
            if( extraMetadata ) {
                const frantic::channels::property_map* chSourceMetadata =
                    extraMetadata->get_channel_metadata( ch.name() );
                if( chSourceMetadata )
                    channelMetadata[ch.name()] = *chSourceMetadata;
            }
        }
        return boost::shared_ptr<streams::particle_ostream>(
            new e57_particle_ostream( file, particleChannelMap, expectedParticleCount, &metadata, &channelMetadata ) );
#endif
    } else {
        bool enablePRT2Saving = m_pImpl->m_enablePRT2Saving;

        boost::filesystem::path tempDirPath = m_pImpl->temp_directory_path();
        std::size_t fileBufferSize = m_pImpl->get_buffer_size();
        std::size_t writeBufferSize = m_pImpl->get_write_size();

        bool hasCoordSystem =
            this->m_pImpl->get_coordinate_system() != frantic::graphics::coordinate_system::unspecified;
        bool hasUnitInMeters = this->m_pImpl->get_length_unit_in_micrometers() != 0.0;

        frantic::channels::channel_map metadataMap;
        if( hasCoordSystem )
            prt::add_coordinate_system( metadataMap );
        if( hasUnitInMeters ) {
            if( enablePRT2Saving ) {
                prt::length_unit_in_micrometers::add_channel( metadataMap );
            } else {
                prt::length_unit_in_meters::add_channel( metadataMap );
            }
        }
        metadataMap.end_channel_definition();

        frantic::channels::property_map metadata( metadataMap );
        if( extraMetadata )
            metadata.merge_property_map( extraMetadata->get_general_metadata() );
        if( hasCoordSystem )
            prt::set_coordinate_system( metadata, this->m_pImpl->get_coordinate_system() );
        if( hasUnitInMeters ) {
            if( enablePRT2Saving ) {
                prt::length_unit_in_micrometers::set_value( metadata, this->m_pImpl->get_length_unit_in_micrometers() );
            } else {
                prt::length_unit_in_meters::set_value( metadata,
                                                       this->m_pImpl->get_length_unit_in_micrometers() / 1.0e6 );
            }
        }

        // TODO: Merge in user-supplied metadata first.

        frantic::channels::channel_map channelMetadataMap;
        prt::add_channel_interpretation( channelMetadataMap );
        channelMetadataMap.end_channel_definition();

        std::map<frantic::tstring, frantic::channels::property_map> channelMetadata;

        for( std::size_t i = 0, iEnd = particleChannelMapForFile.channel_count(); i < iEnd; ++i ) {
            const frantic::channels::channel& ch = particleChannelMapForFile[i];

            if( extraMetadata ) {
                const frantic::channels::property_map* chSourceMetadata =
                    extraMetadata->get_channel_metadata( ch.name() );
                if( chSourceMetadata ) {
                    channelMetadata[ch.name()] = *chSourceMetadata;
                }
            }

            const frantic::tstring chMeaning = prt::get_default_prt2_channel_interpretation( ch.name() );

            if( !chMeaning.empty() ) {
                frantic::channels::property_map& chMetadata = channelMetadata[ch.name()];
                chMetadata.set_channel_map( channelMetadataMap );
                prt::set_channel_interpretation( chMetadata, chMeaning );
            }
        }

        if( enablePRT2Saving ) {
            boost::shared_ptr<prt2_particle_ostream> pResult;
            pResult.reset( new prt2_particle_ostream( file, particleChannelMap, particleChannelMapForFile,
                                                      m_pImpl->m_compressionScheme, true, tempDirPath, &metadata,
                                                      &channelMetadata, m_pImpl->m_desiredChunkSizeInBytes ) );
            return pResult;
        } else {
            const particle_file_metadata prt1Metadata =
                prt::convert_metadata_prt2_to_prt1( particle_file_metadata( metadata, channelMetadata ) );

            boost::shared_ptr<prt_particle_ostream> pResult;
            pResult.reset( new prt_particle_ostream(
                file, particleChannelMap, particleChannelMapForFile, expectedParticleCount, zlibCompressionLevel,
                tempDirPath, fileBufferSize, writeBufferSize, &prt1Metadata.get_general_metadata(),
                &prt1Metadata.get_all_channel_metadata() ) );
            return pResult;
        }
    }
}

boost::shared_ptr<streams::particle_ostream> particle_file_stream_factory_object::create_ostream(
    const frantic::tstring& file, const channel_map& particleChannelMap, const channel_map& particleChannelMapForFile,
    boost::int64_t expectedParticleCount1, int zlibCompressionLevel ) {

    return create_ostream( file, particleChannelMap, particleChannelMapForFile, NULL, expectedParticleCount1,
                           zlibCompressionLevel );
}

void particle_file_stream_factory_object::set_to_defaults() {
    m_pImpl->m_coordinateSystem = frantic::graphics::coordinate_system::right_handed_zup;
    m_pImpl->m_tempDirectory.clear();
    m_pImpl->m_scaleToMicrometers = 0.0;
    m_pImpl->m_frameRate.first = m_pImpl->m_frameRate.second = 0;
    m_pImpl->m_positionTypeHint = channels::data_type_invalid;
    m_pImpl->m_enablePRT2Saving = false;
    m_pImpl->m_compressionScheme = prtfile::prt2_compression_default;
    m_pImpl->m_desiredChunkSizeInBytes = 2000000;
    m_pImpl->m_overrideMetadata.clear();
}

// Legacy function, use the factory object instead (particle_file_stream_factory_object).
boost::shared_ptr<particle_istream> particle_file_istream_factory( const frantic::tstring& file ) {
    particle_file_stream_factory_object defaultStreamFactoryObject; // The default loading code will use right-hand
                                                                    // z-up. This is for backwards compatibility with
                                                                    // old code, specifically for reading BIN files.
    return defaultStreamFactoryObject.create_istream( file );
}

// Legacy function, use the factory object instead (particle_file_stream_factory_object).
boost::shared_ptr<particle_istream> particle_file_istream_factory( const frantic::tstring& file,
                                                                   const channel_map& particleChannelMap ) {
    particle_file_stream_factory_object defaultStreamFactoryObject; // The default loading code will use right-hand
                                                                    // z-up. This is for backwards compatibility with
                                                                    // old code, specifically for reading BIN files.
    return defaultStreamFactoryObject.create_istream( file, particleChannelMap );
}

// Legacy function, use the factory object instead (particle_file_stream_factory_object).
// This creats an output file stream, using the provided particle channel maps for particle memory and particle disk
// layout, respectively
// @param zlibCompressionLevel Optionally specify the zlib compression level when saving PRT files.
boost::shared_ptr<streams::particle_ostream>
particle_file_ostream_factory( const frantic::tstring& file, const channel_map& particleChannelMap,
                               const channel_map& particleChannelMapForFile, boost::int64_t expectedParticleCount,
                               int zlibCompressionLevel ) {
    particle_file_stream_factory_object defaultStreamFactoryObject; // The default loading code will use right-hand
                                                                    // z-up. This is for backwards compatibility with
                                                                    // old code, specifically for reading BIN files.
    return defaultStreamFactoryObject.create_ostream( file, particleChannelMap, particleChannelMapForFile,
                                                      expectedParticleCount, zlibCompressionLevel );
}

void save_particle_stream( boost::shared_ptr<streams::particle_istream> pin,
                           boost::shared_ptr<streams::particle_ostream> pout, logging::progress_logger& progress ) {
    static const int CHUNK_SIZE = 50000;

    if( pin->get_channel_map() != pout->get_channel_map() )
        throw std::runtime_error(
            "save_particle_stream_to_file() The input and output streams did not have the same channel map" );

    std::size_t particleSize = pin->particle_size();

    boost::scoped_array<char> buffer( new char[particleSize * CHUNK_SIZE] );

    char* p;
    bool notDone;
    std::size_t particleCount;
    do {
        p = buffer.get();
        particleCount = CHUNK_SIZE;

        notDone = pin->get_particles( p, particleCount );
        for( std::size_t i = 0; i < particleCount; ++i, p += particleSize )
            pout->put_particle( p );
        progress.update_progress( pin->particle_progress_index(), pin->particle_progress_count() );
    } while( notDone );

    pout->close();
}

void save_particle_stream_to_file( boost::shared_ptr<streams::particle_istream> pin,
                                   const frantic::tstring& particleFile, logging::progress_logger& progress ) {
    // Open the output stream using the default factory
    boost::shared_ptr<streams::particle_ostream> pout = particle_file_ostream_factory(
        particleFile, pin->get_channel_map(), pin->get_channel_map(), pin->particle_count() );

    save_particle_stream( pin, pout, progress );
}

void save_particle_stream_to_file( boost::shared_ptr<streams::particle_istream> pin,
                                   const frantic::tstring& particleFile ) {
    frantic::logging::null_progress_logger npl;
    save_particle_stream_to_file( pin, particleFile, npl );
}
} // namespace particles
} // namespace frantic
