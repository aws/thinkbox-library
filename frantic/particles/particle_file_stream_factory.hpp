// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <boost/version.hpp>
#include <frantic/channels/channel_column_map.hpp>
#include <frantic/particles/particle_file_metadata.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/particles/streams/particle_ostream.hpp>
#include <frantic/prtfile/prt2_common.hpp>
#include <frantic/strings/tstring.hpp>

// TODO: remove this once we're using newer boost everywhere
#if BOOST_VERSION >= 105000
namespace boost {
namespace filesystem {
class path;
}
} // namespace boost
#else
namespace boost {
namespace filesystem3 {
class path;
}
namespace filesystem {
using filesystem3::path;
}
} // namespace boost
#endif

namespace frantic {
namespace logging {
class progress_logger;
}
} // namespace frantic

namespace frantic {
namespace particles {

using frantic::channels::channel_column_map;
using frantic::channels::channel_map;

/**
 * This class implements a factory object that creates particle_istream and particle_ostream instances when given a file
 * path. This object is  superior to using free functions, since we can configure extra options and functionality that
 * controls the object created.
 */
class particle_file_stream_factory_object {
  public:
    particle_file_stream_factory_object();

    ~particle_file_stream_factory_object();

    /**
     * \return The configured type of coordinate system this factory will be using with the particle_istreams produced.
     */
    frantic::graphics::coordinate_system::option get_coordinate_system();

    /**
     * \param coordinateSystem The new coordinate system that streams produced by this factory will use.
     */
    void set_coordinate_system( frantic::graphics::coordinate_system::option coordinateSystem );

    /**
     * \deprecated As of PRT2, micrometers is used.
     *
     * \return The value to scale length measures by in order to get meters. Will be 1.0 if the data is already in
     * meters.
     */
    double get_length_unit_in_meters() const;

    /**
     * \deprecated As of PRT2, micrometers is used.
     *
     * \param scaleToMeters The value to scale length measures by in order to get meters.
     */
    void set_length_unit_in_meters( double scaleToMeters );

    /**
     * \return The value to scale length measures by in order to get meters. Will be 1.0 if the data is already in
     * meters.
     */
    double get_length_unit_in_micrometers() const;

    /**
     * \param scaleToMeters The value to scale length measures by in order to get meters.
     */
    void set_length_unit_in_micrometers( double scaleToMicrometers );

    /**
     * Returns the desired dtype for the Position channel (data_type_undefined, data_type_float32, or
     * data_type_float64).
     */
    channels::data_type_t get_position_type_hint() const;

    /**
     * Sets the dtype desired for the Position channel (data_type_undefined, data_type_float32, or data_type_float64).
     * When data_type_float64 is selected, the input will typically return float64 if there is extra precision
     * available, and float32 otherwise.
     */
    void set_position_type_hint( channels::data_type_t positionTypeHint );

    /**
     * \return The frame rate associated with the particle data, in frames per second. Encoded as a rational number.
     */
    std::pair<unsigned, unsigned> get_frame_rate() const;

    /**
     * \param n The numerator of the frame rate in frames per second.
     * \param d The denominator of the frame rate in frames per second.
     */
    void set_frame_rate( unsigned n, unsigned d );

    /**
     * \return the current requested 'emulated' channel map for non-metadata streams
     */
    const frantic::channels::channel_column_map get_channel_column_map() const;

    /**
     * Gets the extra metadata that should be appended onto any file loaded using this factory.
     * Under most normal circumstances, this should be left blank.
     *
     * @return the metadata object on this factory
     */
    const particle_file_metadata& get_override_metadata() const;

    /**
     * Sets the extra metadata to be appended onto any file loaded using this factory.
     * This allows you to, for example, set the native units of a file which has no scale information
     *
     * @param overrideMetadata the set of metadata and channel metadata to be applied
     */
    void set_override_metadata( const particle_file_metadata& overrideMetadata );

    /**
     * Enables saving as PRT2-format files instead of the default PRT1.
     */
    void enable_prt2_saving(
        frantic::prtfile::prt2_compression_t compressionScheme = frantic::prtfile::prt2_compression_default,
        intptr_t desiredChunkSizeInBytes = 2000000 );

    /**
     * Creates a new particle_istream instance from the specified file. The returned particle_istream will have its
     * current channel map equal to its native one.
     * \param file The path to the file to create a particle_istream instance from.
     * \return A new instance of a particle_istream subclass that can read the specified file.
     * \throws runtime_error If the specified file does not exist, or could not be interpreted as a known particle file
     * format.
     */
    boost::shared_ptr<streams::particle_istream> create_istream( const frantic::tstring& file,
                                                                 particle_file_metadata* outMetadata = NULL );

    /**
     * Creates a new particle_istream instance from the specified file, using a specified channel map for the result.
     * \param file The path to the file to create a particle_istream instance from.
     * \param particleChannelMap The channel map to assign to the newly created instance.
     * \return A new instance of a particle_istream subclass that can read the specified file.
     * \throws runtime_error If the specified file does not exist, or could not be interpreted as a known particle file
     * format.
     */
    boost::shared_ptr<streams::particle_istream> create_istream( const frantic::tstring& file,
                                                                 const channel_map& particleChannelMap,
                                                                 particle_file_metadata* outMetadata = NULL );

    /**
     * Creates a new particle_ostream instace for writing to the specified file. The type of particle_ostream is
     * determined from the file's extension.
     * \param file The path to the file to create.
     * \param particleChannelMap The particle layout of data provided to the created particle_istream
     * \param particleChannelMapForFile The particle layout to use in the actual file (can be a subset, or modified
     * version of particleChannelMap). This can be ignored by formats with a fixed layout.
     * \param extraMetadata additional metadata to save to the file
     * \param expectedParticleCount The number of particles that will be written to the particle_ostream. Can be -1 to
     * indicate the amount is not known beforehand.
     * \param zlibCompressionLevel For formats that use zlib, this indicates the compression level to use (see zlib.h
     * for more details). Can be -1 to use the default compression level.
     */
    boost::shared_ptr<streams::particle_ostream>
    create_ostream( const frantic::tstring& file, const channel_map& particleChannelMap,
                    const channel_map& particleChannelMapForFile, const particle_file_metadata* extraMetadata,
                    boost::int64_t expectedParticleCount = -1, int zlibCompressionLevel = -1 );

    boost::shared_ptr<streams::particle_ostream> create_ostream( const frantic::tstring& file,
                                                                 const channel_map& particleChannelMap,
                                                                 const channel_map& particleChannelMapForFile,
                                                                 boost::int64_t expectedParticleCount = -1,
                                                                 int zlibCompressionLevel = -1 );

    /**
     * Resets the factory object to a default configuration compatible with the legacy implementations of
     * particle_file_ostream_factory().
     */
    void set_to_defaults();

    /**
     * Sets the folder that ostream objects will use to write particle files before moving them to their final
     * destination. This is initialized from the environment variable
     * 'PRT_TEMPDIR' or an OS-specific temp location if the former doesn't exist.
     * \param tempDirectory The directory to initially write particle data files to before moving them to the target
     * destination.
     * \pre boost::filesystem::exists(tempDirectory) && boost::filesystem::is_directory(tempDirectory)
     */
    void set_temp_directory( const boost::filesystem::path& tempDirectory );

  private:
    // Disabling copy constructor and assignment operator until we determine they are necessary. They need to be
    // implemented to copy the pimpl inside.
    particle_file_stream_factory_object( const particle_file_stream_factory_object& );
    particle_file_stream_factory_object& operator=( const particle_file_stream_factory_object& );

  private:
    class impl_type;

    impl_type* m_pImpl; // Can't use a scoped_pointer since impl_type is incomplete.
};

// This creates an input file stream, using the file's native particle channel map for the particle memory layout
boost::shared_ptr<streams::particle_istream> particle_file_istream_factory( const frantic::tstring& file );

// This creates an input file stream, using the provided particle channel map for the particle memory layout
boost::shared_ptr<streams::particle_istream> particle_file_istream_factory( const frantic::tstring& file,
                                                                            const channel_map& particleChannelMap );

// This creates an output file stream, using the provided particle channel maps for particle memory and particle disk
// layout, respectively
boost::shared_ptr<streams::particle_ostream>
particle_file_ostream_factory( const frantic::tstring& file, const channel_map& particleChannelMap,
                               const channel_map& particleChannelMapForFile, boost::int64_t expectedParticleCount = -1,
                               int zlibCompressionLevel = -1 );

// This will save a particle istream to a file, by using particle_file_ostream_factory to make the requisite ostream
void save_particle_stream_to_file( boost::shared_ptr<streams::particle_istream> pin,
                                   const frantic::tstring& particleFile, logging::progress_logger& progress );

// OVERLOAD
void save_particle_stream_to_file( boost::shared_ptr<streams::particle_istream> pin,
                                   const frantic::tstring& particleFile );

// Save the particles from pin to pout in chunks suitable for taking advantage of concurrency.
void save_particle_stream( boost::shared_ptr<streams::particle_istream> pin,
                           boost::shared_ptr<streams::particle_ostream> pout, logging::progress_logger& progress );
} // namespace particles
} // namespace frantic
