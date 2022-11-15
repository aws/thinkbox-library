// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/smart_ptr.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/voxel_sampler_interface.hpp>

namespace frantic {
namespace particles {
namespace streams {

/**
 * This class is a base class for particle_istream objects that produce particles
 * based on the defined voxels in a rle_levelset object. The base class is responsible for
 * for maintaining and advancing a rle_defined_iterator. The base class is also responsible
 * for data access to and from the levelset. The child classes should only have to implement
 * a new constructor, and the get_particle() functions.
 */
class rle_levelset_particle_istream : public particle_istream {
  protected:
    /**
     * This internal struct holds the data required to copy a rle_levelset channel into a
     * particle channel, with possibly a different data type.
     */
    struct channel_connection {
        frantic::channels::channel_general_accessor particleAccessor;
        frantic::volumetrics::levelset::rle_channel_general_accessor levelsetAccessor;
        frantic::channels::channel_type_convertor_function_t conversionFunction;
    };

    boost::shared_ptr<frantic::volumetrics::levelset::rle_level_set> m_pLevelset;
    frantic::volumetrics::levelset::rle_defined_iterator m_iter, m_endIter;

    boost::int32_t m_cachedIndices[8];

    //'m_innerDistance' and 'm_outerDistance' define a region of the levelset that is actively
    // generating particles.
    float m_outerDistance;
    float m_innerDistance;

    bool m_compensateDensity;
    float m_compensationFactor;

    frantic::channels::channel_map m_outMap;
    frantic::channels::channel_map m_nativeMap;

    boost::int64_t m_particleIndex;
    boost::int64_t m_particleProgressIndex;
    boost::int64_t m_particleProgressCount;
    boost::scoped_array<char> m_defaultParticle;

    // Accessors for the various extra channels.
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_normalAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_signedDistanceGradientAccessor;
    frantic::channels::channel_cvt_accessor<float> m_distAccessor;
    frantic::channels::channel_cvt_accessor<float> m_densityAccessor;
    std::vector<channel_connection> m_channels;

    frantic::volumetrics::voxel_sampler_interface_ptr m_pParticleGenerator;

  protected:
    /**
     * This function will copy the channel data for a given voxel into the provided particle. This
     * version will trilerp the data values using the 8 indices in trilerpIndices, and the 8 weights
     * in trilerpWeights.
     *
     * @param localCoord The voxel coordinate (relative the centered sample location).
     * @param localDistance The signed-distance value at localCoord.
     * @param trilerpIndices The 8 indices of the voxels bracketing localCoord.
     * @param trilerpWeights The 8 trilerp weights of the voxels bracketing localCoord.
     * @param pParticle A pointer to the particle to copy the data channels into.
     */
    void copy_channel_data( const frantic::graphics::vector3f& localCoord, float localDistance,
                            boost::int32_t trilerpIndices[], float trilerpWeights[], char* pParticle );

    bool advance_iterator();

  public:
    rle_levelset_particle_istream( const frantic::channels::channel_map& pcm,
                                   boost::shared_ptr<frantic::volumetrics::levelset::rle_level_set> pLevelset,
                                   frantic::volumetrics::voxel_sampler_interface_ptr pParticleGenerator,
                                   float innderDistance, float outerDistance, bool compensateDensity );

    virtual ~rle_levelset_particle_istream();

    virtual void close() { m_pLevelset.reset(); }
    virtual frantic::tstring name() const { return _T("rle_levelset_particle_istream"); }

    virtual std::size_t particle_size() const { return m_outMap.structure_size(); }

    virtual boost::int64_t particle_count() const { return -1; }
    virtual boost::int64_t particle_index() const { return m_particleIndex; }
    virtual boost::int64_t particle_count_left() const { return -1; }
    virtual boost::int64_t particle_progress_count() const { return m_particleProgressCount; }
    virtual boost::int64_t particle_progress_index() const { return m_particleProgressIndex; }

    virtual const frantic::channels::channel_map& get_channel_map() const { return m_outMap; }
    virtual const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeMap; }

    virtual void set_channel_map( const frantic::channels::channel_map& particleChannelMap );
    virtual void set_default_particle( char* rawParticleBuffer );

    virtual bool get_particle( char* rawParticleBuffer );
    virtual bool get_particles( char* rawParticleBuffer, std::size_t& numParticles );
};

} // namespace streams
} // namespace particles
} // namespace frantic
