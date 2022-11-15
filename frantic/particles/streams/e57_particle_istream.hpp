// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#if defined( E57_AVAILABLE )

#include <E57Format/E57Format.h>

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread/mutex.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/files/files.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/strings/utf8.hpp>

namespace frantic {
namespace particles {
namespace streams {

/**
 * E57 File Format Reader
 *
 * http://www.libe57.org/
 * http://www.ri.cmu.edu/publication_view.html?pub_id=6767
 */
class e57_particle_istream : public particle_istream {
  private:
    frantic::tstring m_filename;
    frantic::channels::channel_map m_onDiskParticleChannelMap, m_particleChannelMap;
    frantic::channels::channel_map_adaptor m_pcmAdaptor;
    boost::int64_t m_particleCount, m_currentParticleIndex;
    /// The guid of the scan to load. If equal to "", will load the entire file.
    frantic::tstring m_scanGuid;
    /// Should the data be transformed according to its scanner transform before being fetched?
    bool m_transformData;

    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3fd> m_posAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_colorAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_normalAccessor;
    frantic::channels::channel_cvt_accessor<float> m_intensityAccessor;
    frantic::channels::channel_accessor<boost::uint32_t> m_scannerIndexAccessor;

    frantic::particles::particle_file_metadata m_metadata;
    frantic::graphics::boundbox3fd m_bounds;
    bool m_boundsValid;

    std::vector<char> m_defaultParticleBuffer;
    std::vector<char> m_tempParticleBuffer;

    float m_intensityOffset;
    float m_intensityCoefficient;
    bool m_hasIntensity;

    static const int N = 10; // used for setting the amount of points to read in at a time to SourceDestBuffer's

    boost::int64_t m_currentScanIndex, m_scanCount;

    boost::scoped_ptr<e57::ImageFile> m_imageFile;
    std::vector<e57::SourceDestBuffer> m_destBuffers;
    boost::scoped_ptr<e57::CompressedVectorReader> m_reader;
    unsigned m_bufferParticleCount;
    unsigned m_bufferParticleIndex;
    double m_x[N];
    double m_y[N];
    double m_z[N];
    float m_intensity[N];
    float m_red[N];
    float m_green[N];
    float m_blue[N];
    float m_normalX[N];
    float m_normalY[N];
    float m_normalZ[N];

    std::vector<frantic::graphics::transform4fd> m_scannerTransforms;
    boost::uint32_t m_scannerTransformsIndex;

    // The input coordinates may be in any of these coordinate systems
    enum e57_coord_sys { cartesian_coord, spherical_coord };

    e57_coord_sys m_coordSys;
    bool m_hasColor;
    bool m_hasTransform;
    bool m_hasNormal;
    // The max number for color that the file is using for the scan
    float m_maxColor;

    // Private copy constructor to disable copying
    e57_particle_istream( const e57_particle_istream& ); // not implemented

    // Private assignment operator to disable assignment
    e57_particle_istream& operator=( const e57_particle_istream& ); // not implemented

    void initialize_stream( channels::data_type_t positionTypeHint );

    void open_file( channels::data_type_t positionTypeHint );

    void load_new_scan();

    void load_requested_scan();

    void load_buffers();

    void setup_metadata();

  public:
    e57_particle_istream( const frantic::tstring& file,
                          channels::data_type_t positionTypeHint = channels::data_type_invalid );

    /**
     * Create an e57_particle_istream for the purposes of loading a single scan, identified by `guid`.
     * If multiple scans with the specified `guid` exist, the stream will only load the first one.
     * @param file The file name.
     * @param guid The guid identifying the desired scan.
     * @param transformData Should the scan be transformed according to its transform?
     * @param positionTypeHint The precision of the position channel.
     */
    e57_particle_istream( const frantic::tstring& file, const frantic::tstring& guid, bool transformData = false,
                          channels::data_type_t positionTypeHint = channels::data_type_invalid );

    virtual ~e57_particle_istream() { close(); }

    void close();

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_filename; }

    boost::int64_t particle_count() const { return m_particleCount; }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return m_particleCount - m_currentParticleIndex - 1; }

    boost::int64_t particle_progress_count() const { return m_particleCount; }

    boost::int64_t particle_progress_index() const { return m_currentParticleIndex; }

    const channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    const channels::channel_map& get_native_channel_map() const { return m_onDiskParticleChannelMap; }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const frantic::channels::channel_map& particleChannelMap );

    void set_default_particle( char* buffer );

    bool get_particle( char* rawParticleBuffer );

    bool get_particles( char* particleBuffer, std::size_t& numParticles );

    /**
     * Get the metadata field for the particle stream
     */
    const frantic::particles::particle_file_metadata& get_metadata() const;
};
} // namespace streams
} // namespace particles
} // namespace frantic

#endif
