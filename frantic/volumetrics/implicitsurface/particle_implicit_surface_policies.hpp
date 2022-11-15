// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
//#define FRANTIC_DISABLE_THREADS

#pragma once

#include <boost/function.hpp>
#include <boost/math/special_functions/pow.hpp>

#include <frantic/channels/channel_map_weighted_sum.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/graphics/aligned_array.hpp>
#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/particles/particle_grid_tree.hpp>
#include <frantic/particles/sph/sph_kernel_functions.hpp>
#include <frantic/simd/float_v.hpp>
#include <frantic/volumetrics/implicitsurface/detail/xyzr_packet_array.hpp>
#include <frantic/volumetrics/rle_plane.hpp>
#include <frantic/volumetrics/run_tree.hpp>

#pragma warning( push )
#pragma warning( disable : 4512 4100 4244 4245 )
#include <tbb/atomic.h>
#include <tbb/spin_mutex.h>
#include <tbb/spin_rw_mutex.h>
#pragma warning( pop )

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

class shared_progress_logger_proxy {
    boost::int32_t m_isCancelled;
    char m_padding[60];
    tbb::atomic<std::size_t> m_progress;
    std::size_t m_progressEnd;

  public:
    shared_progress_logger_proxy()
        : m_progressEnd( 100 )
        , m_isCancelled( 0 ) {
        m_progress = 0;
    }

    void cancel() { m_isCancelled = 1; }

    bool is_cancelled() { return m_isCancelled != 0; }

    void set_progress_end( std::size_t progressEnd ) { m_progressEnd = progressEnd; }

    float get_progress() { return 100.f * ( float( m_progress ) / m_progressEnd ); }

    void add_progress( std::size_t count ) { m_progress += count; }

    void inc_progress() { ++m_progress; }
};

/**
 *  Act as an intermediate between a progress_logger and a shared_progress_logger_proxy.
 *
 *  The progress_logger is handled on the main thread while the user-specified function
 * is run on another thread.
 *
 *  This is intended to help keep the UI responsive while updating a progress_logger.
 */
class shared_progress_logger_adapter {
    frantic::logging::progress_logger& m_progressLogger;
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy m_proxy;

    /**
     *  Run a user specified function, and return a value indicating whether an
     * exception was thrown, and a string holding the error message if applicable.
     *
     *  A return value of 0 indicates that no exception was thrown.
     *  1 indicates a frantic::logging::progress_cancel_exception exception.
     *  2 indicates a std::exception.
     *  3 indicates any other exception.
     */
    static std::pair<int, std::string> progress_logger_adapter_threadproc( boost::function<void( void )>& f );

    shared_progress_logger_adapter& operator=( const shared_progress_logger_adapter& ); // not implemented
  public:
    shared_progress_logger_adapter( frantic::logging::progress_logger& progressLogger );
    virtual ~shared_progress_logger_adapter();

    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy* get_proxy();
    void run( boost::function<void( void )>& f );
};

namespace detail {

struct particle_order {
    frantic::graphics::vector3 voxelCoord;
    boost::int32_t index;

    void set( const frantic::graphics::vector3& voxelCoord, boost::int32_t index ) {
        this->voxelCoord = voxelCoord;
        this->index = index;
    }
};

struct particle_order_comparison {
    bool operator()( const particle_order& a, const particle_order& b ) const { return a.voxelCoord < b.voxelCoord; }
};

struct voxel_coord_hasher {
    size_t operator()( const frantic::graphics::vector3& k ) const {
        return frantic::graphics::vector3_hasher::hash( k );
    }
};

} // namespace detail

template <class Implementation>
class particle_is_policy_base {
  public:
    particle_is_policy_base( frantic::particles::particle_grid_tree& particles )
        : m_particles( particles ) {}

    const frantic::channels::channel_map& get_channel_map() const { return m_particles.get_channel_map(); }

    void get_particles_in_range( const frantic::graphics::vector3f& position, std::vector<char*>& outParticles ) const {
        outParticles.clear();
        const float compactSupport = derived().get_compact_support();
        m_particles.get_particles_in_range( position, compactSupport, outParticles );
    }

    ////////
    // Functions for dealing with the named channels
    ////////

    /**
     *	This retrieves the named channels in the particles that the policy is operating on.
     *	It adds the names to the existing vector, so if you want just the names this function provides,
     *	make sure to initialize the output vector to empty yourself.
     *
     *	@param	outNames	A vector of channel name strings.
     */
    void get_channel_names( std::vector<frantic::tstring>& outNames ) const {
        outNames.reserve( outNames.size() + get_channel_map().channel_count() );
        for( unsigned i = 0; i < get_channel_map().channel_count(); ++i ) {
            const std::string& channelName = get_channel_map()[i].name();
            outNames.push_back( channelName );
        }
    }

    /**
     *	This retrieves the named channels in the particles that the policy is operating on.
     *	It adds the names to the existing vector, so if you want just the names this function provides,
     *	make sure to initialize the output vector to empty yourself.  It uses the channel propagation
     *	policy to determine which channels should be returned.
     *
     *	@param	cpp		A policy that determines which particle channels should be propagated to the output.
     *	@param	outNames	A vector of channel name strings.
     */
    void get_channel_names( const frantic::channels::channel_propagation_policy& cpp,
                            std::vector<frantic::tstring>& outNames ) const {
        outNames.reserve( outNames.size() + get_channel_map().channel_count() );
        for( unsigned i = 0; i < get_channel_map().channel_count(); ++i ) {
            const frantic::tstring& channelName = get_channel_map()[i].name();
            if( cpp.is_channel_included( channelName ) )
                outNames.push_back( channelName );
        }
    }

    /**
     * This retrieves all the channel information associated with the named channels
     * in isp that are included in the channel propagation policy.
     *
     * @param	cpp						A channel propagation policy.
     * @param	outChannelName				A vector of strings into which the channel names are added.
     * @param	outChannelType			A vector of data_type_t for each channel
     * @param	outChannelArity			A vector of size_t arities for each channel
     * @param	outChannelPrimitiveSize	A vector of size_t primitive sizes for each channel.
     */
    void get_channel_info( const frantic::channels::channel_propagation_policy& cpp,
                           std::vector<frantic::tstring>& outChannelName,
                           std::vector<frantic::channels::data_type_t>& outChannelType,
                           std::vector<size_t>& outChannelArity, std::vector<size_t>& outChannelPrimitiveSize ) const {

        size_t channelCount = get_channel_map().channel_count();

        if( cpp.is_channel_included( _T("SignedDistance") ) ) {
            outChannelName.push_back( _T("SignedDistance") );
            outChannelType.push_back( frantic::channels::data_type_float32 );
            outChannelArity.push_back( 1 );
            outChannelPrimitiveSize.push_back( 4 );
        }

        for( size_t i = 0; i < channelCount; ++i ) {
            if( cpp.is_channel_included( get_channel_map()[i].name() ) ) {
                outChannelName.push_back( get_channel_map()[i].name() );
                outChannelType.push_back( get_channel_map()[i].data_type() );
                outChannelArity.push_back( get_channel_map()[i].arity() );
                outChannelPrimitiveSize.push_back( get_channel_map()[i].primitive_size() );
            }
        }
    }

    /**
     *	This function is called so that the channel accessors can be prepared, and to report the data types for the
     *algorithm to
     *	create corresponding channels in the output data structure.
     */
    void
    prepare_channel_accessors( const std::vector<frantic::tstring>& names,
                               std::vector<std::pair<std::size_t, frantic::channels::data_type_t>>& outDataTypes ) {
        outDataTypes.clear();
        outDataTypes.reserve( names.size() );
        m_inputChannels.clear();
        m_inputChannels.reserve( names.size() );
        frantic::channels::channel_propagation_policy cpp( true );
        for( unsigned i = 0; i < names.size(); ++i ) {
            m_inputChannels.push_back( get_channel_map().get_general_accessor( names[i] ) );
            outDataTypes.push_back(
                std::make_pair( m_inputChannels.back().arity(), m_inputChannels.back().data_type() ) );
            cpp.add_channel( names[i] );
        }
        m_channelMapWeightedSum = frantic::channels::channel_map_weighted_sum( m_particles.get_channel_map(), cpp );
    }

  protected:
    std::vector<frantic::channels::channel_general_accessor> m_inputChannels;
    frantic::particles::particle_grid_tree& m_particles;
    frantic::channels::channel_map_weighted_sum m_channelMapWeightedSum;

  private:
    const Implementation& derived() const { return *static_cast<const Implementation*>( this ); }

    particle_is_policy_base& operator=( const particle_is_policy_base& ); // not implemented
};

namespace detail {
struct union_of_spheres_sparse_data {

    void reset_for_dense_plane_evaluation(
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::channels::channel_accessor<float>& radiusAccessor, const float particleRadiusToEffectRadiusScale,
        float implicitThreshold, const frantic::graphics::boundbox3& voxelExtents,
        const frantic::volumetrics::voxel_coord_system& meshingVCS, float* voxelCornerDensities );

    void reset_for_dense_block_evaluation(
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::channels::channel_accessor<float>& radiusAccessor, const float particleRadiusToEffectRadiusScale,
        float implicitThreshold, const frantic::graphics::boundbox3& voxelExtents,
        const frantic::volumetrics::voxel_coord_system& meshingVCS, float* voxelCornerDensities );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> position;
    frantic::channels::channel_accessor<float> radius;
    float particleRadiusToEffectRadiusScale;
    frantic::graphics::boundbox3 voxelExtents;
    frantic::volumetrics::voxel_coord_system meshingVCS; // voxel coord system for meshing

    // output channel data/accessors
    std::vector<char*> channelData;
    std::vector<frantic::channels::channel_general_accessor> channelAccessors;

    // only used when sparse run data is needed (ie, particle to level set conversion)
    bool collectRunData;
    std::vector<frantic::volumetrics::run_tree> runTrees;
    std::vector<std::vector<char*>> extentData;
    std::vector<std::vector<char>> initChannelData;
    int signedDistanceChannel;

    // for multithreading
    bool useMutexes;
    boost::shared_array<tbb::spin_mutex> mutexes;

    std::size_t gridParticleCount;
};
// per-thread storage that is reused between calls to add_edge_isosurface_vertex_to_mesh
// TODO: newer versions of tbb support thread-local storage.  Switch to that instead if it's more appropriate.
struct union_of_spheres_vertex_workspace {
    typedef std::vector<char*> particles_t;
    particles_t particles;
};
} // namespace detail

/**
 * This is a policy class which defines the implicit surface of a simple union of bounding spheres for a
 * particle system.
 */
class particle_union_of_spheres_is_policy : public particle_is_policy_base<particle_union_of_spheres_is_policy> {
  public:
    typedef detail::union_of_spheres_sparse_data sparse_voxel_corner_density_workspace_t;
    typedef detail::union_of_spheres_vertex_workspace vertex_workspace_t;

  private:
    // The union_of_spheres parameters
    float m_maximumParticleRadius, m_particleRadiusToEffectRadiusScale, m_implicitThreshold;
    // How much to refine the vertex positions
    int m_vertexRefinement;
    float m_vertexRefinementEpsilon;

    // mutable tbb::atomic<std::size_t> m_vrEvalCount;

    // The voxel grid used to sample the implicit surface.
    frantic::volumetrics::voxel_coord_system m_meshingVCS;

    // Position accessor
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_accessor<float> m_radiusAccessor;

    // The bounding box of the particles
    frantic::graphics::boundbox3 m_particleVCSBounds;

    // std::vector<frantic::channels::channel_general_accessor> m_inputChannels;

    // When an isosurface location is found, these parameters specify where it is.  In other implementations of the mc
    // policy,
    // the array of filter weights and data indices would be stored to reuse for all the different named channels.
    frantic::graphics::vector3f m_isosurfaceLocationWorld;
    // Also collected is the array of particles that contribute to the density value at that location
    std::vector<char*> m_isosurfaceLocationParticles;

    frantic::graphics::vector3f find_edge_isosurface_location( float density0, float density1,
                                                               const frantic::graphics::vector3& voxelCorner0,
                                                               const frantic::graphics::vector3& voxelCorner1,
                                                               std::vector<char*>& particles ) const;

    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              const std::vector<char*>& particlesInRange ) const;
    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert ) const;

    static void get_density( void* userData, char* p );

    // Disable copy constructor
    particle_union_of_spheres_is_policy& operator=( const particle_union_of_spheres_is_policy& ); // not implemented
  public:
    // function docs are in the cpp file
    particle_union_of_spheres_is_policy( frantic::particles::particle_grid_tree& particles, float maximumParticleRadius,
                                         float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                         const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                         int vertexRefinement );
    ~particle_union_of_spheres_is_policy();
    const frantic::volumetrics::voxel_coord_system& get_voxel_coord_system() const;
    float get_default_outside_distance() const { return m_implicitThreshold; };
    frantic::graphics::boundbox3 get_voxel_bounds() const;
    float get_compact_support() const { return m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale; };
    void get_sample_weights( const frantic::graphics::vector3f& position, const std::vector<char*>& particles,
                             std::vector<float>& outWeights ) const;
    void get_affected_blocks( const frantic::graphics::size3f blockSize,
                              std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
                              frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );
    float get_density( const frantic::graphics::vector3f& worldLocation ) const;
    frantic::graphics::vector3f get_gradient( const frantic::graphics::vector3f& worldLocation, float h );
    void fill_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                      std::vector<float>& outVoxelCornerValues );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, frantic::volumetrics::rle_plane& rlp,
                                             detail::union_of_spheres_sparse_data& data );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, detail::union_of_spheres_sparse_data& data );
    std::size_t fill_voxel_corner_densities( const frantic::graphics::boundbox3& voxelExtents,
                                             float* voxelCornerDensities,
                                             detail::union_of_spheres_sparse_data& data ) const;
    void fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                   std::vector<frantic::tstring>& channelNames, std::vector<char*>& channelData,
                                   frantic::volumetrics::rle_plane& rlp );
    frantic::graphics::vector3f
    corner_sample_coord_to_world( const frantic::graphics::vector3& voxelCornerCoord ) const;
    void find_edge_isosurface_location( float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
                                        const frantic::graphics::vector3& voxelCorner1 );
    const frantic::graphics::vector3f& get_isosurface_location_world() const;
    void set_isosurface_location_world( const frantic::graphics::vector3f& worldLocation );
    void
    add_vertex_to_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels );
    void get_interface_widths( float& interfaceWidthInside, float& interfaceWidthOutside ) const;
    int get_exterior_region_code() const;
    void add_edge_isosurface_vertex_to_mesh(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh, detail::union_of_spheres_vertex_workspace& workspace ) const;
    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 const frantic::channels::channel_propagation_policy& cpp ) const;
    void
    populate_vertex_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              vertex_workspace_t& workspace ) const;
    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              vertex_workspace_t& workspace ) const;
};

namespace detail {
struct metaball_sparse_data {

    void reset_for_dense_plane_evaluation(
        frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
        frantic::channels::channel_accessor<float> radiusAccessor, float particleRadiusToEffectRadiusScale,
        const frantic::graphics::boundbox3& voxelExtents, const frantic::volumetrics::voxel_coord_system& meshingVCS,
        float* voxelCornerDensities );

    void reset_for_dense_block_evaluation(
        frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
        frantic::channels::channel_accessor<float> radiusAccessor, float particleRadiusToEffectRadiusScale,
        const frantic::graphics::boundbox3& voxelExtents, const frantic::volumetrics::voxel_coord_system& meshingVCS,
        float* voxelCornerDensities );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> position;
    frantic::channels::channel_accessor<float> radius;
    float particleRadiusToEffectRadiusScale;
    frantic::graphics::boundbox3 voxelExtents;
    frantic::volumetrics::voxel_coord_system meshingVCS; // voxel coord system for meshing

    // output channel data/accessors
    std::vector<char*> channelData;
    std::vector<frantic::channels::channel_general_accessor> channelAccessors;

    // only used when sparse run data is needed (ie, particle to level set conversion)
    bool collectRunData;
    std::vector<frantic::volumetrics::run_tree> runTrees;
    std::vector<std::vector<char*>> extentData;
    std::vector<std::vector<char>> initChannelData;
    int signedDistanceChannel;

    // for multithreading
    bool useMutexes;
    boost::shared_array<tbb::spin_mutex> mutexes;

    std::size_t gridParticleCount;
};
// per-thread storage that is reused between calls to add_edge_isosurface_vertex_to_mesh
struct metaball_vertex_workspace {
    typedef std::vector<char*> particles_t;
    particles_t particles;
    std::vector<float> sampleWeights;
    std::vector<char> weightedSum;
};
} // namespace detail

/**
 * This is a policy class which defines the implicit surface of metaballs based on a
 * particle system.
 */
class particle_metaball_is_policy : public particle_is_policy_base<particle_metaball_is_policy> {
  public:
    typedef detail::metaball_sparse_data sparse_voxel_corner_density_workspace_t;
    typedef detail::metaball_vertex_workspace vertex_workspace_t;

  private:
    // The metaball parameters
    float m_maximumParticleRadius, m_particleRadiusToEffectRadiusScale, m_implicitThreshold;
    // How much to refine the vertex positions
    int m_vertexRefinement;
    float m_vertexRefinementEpsilon;

    // mutable tbb::atomic<std::size_t> m_vrEvalCount;

    // The voxel grid used to sample the implicit surface.
    frantic::volumetrics::voxel_coord_system m_meshingVCS;

    // Position accessor
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_accessor<float> m_radiusAccessor;

    // The bounding box of the particles
    frantic::graphics::boundbox3 m_particleVCSBounds;

    // std::vector<frantic::channels::channel_general_accessor> m_inputChannels;

    // When an isosurface location is found, these parameters specify where it is.  In other implementations of the mc
    // policy,
    // the array of filter weights and data indices would be stored to reuse for all the different named channels.
    frantic::graphics::vector3f m_isosurfaceLocationWorld;
    // Also collected is the array of particles that contribute to the density value that that location
    std::vector<char*> m_isosurfaceLocationParticles;

    frantic::graphics::vector3f find_edge_isosurface_location( float density0, float density1,
                                                               const frantic::graphics::vector3& voxelCorner0,
                                                               const frantic::graphics::vector3& voxelCorner1,
                                                               std::vector<char*>& particles ) const;

    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              const std::vector<char*>& particlesInRange, vertex_workspace_t& workspace,
                              const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const;
    void
    populate_vertex_channels( std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert ) const;

    static void get_density( void* userData, char* p );

    // Disable copy constructor
    particle_metaball_is_policy& operator=( const particle_metaball_is_policy& );
    particle_metaball_is_policy( const particle_metaball_is_policy& );

  public:
    // function docs are in the cpp file
    particle_metaball_is_policy( frantic::particles::particle_grid_tree& particles, float maximumParticleRadius,
                                 float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                 const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement );
    ~particle_metaball_is_policy();
    const frantic::volumetrics::voxel_coord_system& get_voxel_coord_system() const;
    frantic::graphics::boundbox3 get_voxel_bounds() const;
    float get_default_outside_distance() const { return m_implicitThreshold; };
    float get_compact_support() const { return m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale; };
    void get_sample_weights( const frantic::graphics::vector3f& position, const std::vector<char*>& particles,
                             std::vector<float>& outWeights ) const;
    void get_affected_blocks( const frantic::graphics::size3f blockSize,
                              std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
                              frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );
    float get_density( const frantic::graphics::vector3f& worldLocation ) const;
    frantic::graphics::vector3f get_gradient( const frantic::graphics::vector3f& worldLocation, float h );
    void fill_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                      std::vector<float>& outVoxelCornerValues );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, frantic::volumetrics::rle_plane& rlp,
                                             detail::metaball_sparse_data& data );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, detail::metaball_sparse_data& data );
    std::size_t fill_voxel_corner_densities( const frantic::graphics::boundbox3& xyzExtents,
                                             float* voxelCornerDensities, detail::metaball_sparse_data& data ) const;
    void fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                   std::vector<frantic::tstring>& channelNames, std::vector<char*>& channelData,
                                   frantic::volumetrics::rle_plane& rlp );
    frantic::graphics::vector3f
    corner_sample_coord_to_world( const frantic::graphics::vector3& voxelCornerCoord ) const;
    void find_edge_isosurface_location( float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
                                        const frantic::graphics::vector3& voxelCorner1 );
    const frantic::graphics::vector3f& get_isosurface_location_world() const;
    void set_isosurface_location_world( const frantic::graphics::vector3f& worldLocation );
    void
    add_vertex_to_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels );
    void get_interface_widths( float& interfaceWidthInside, float& interfaceWidthOutside ) const;
    int get_exterior_region_code() const;
    void add_edge_isosurface_vertex_to_mesh(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh, detail::metaball_vertex_workspace& workspace ) const;
    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 const frantic::channels::channel_propagation_policy& cpp ) const;
    void
    populate_vertex_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              vertex_workspace_t& workspace ) const;
    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace,
                              const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const;
};

/**
 * This is a policy class which defines the implicit surface of Zhu/Bridson based on a
 * particle system.
 *
 */
namespace detail {
/**
 *	Data structures used internally
 */
struct zhu_bridson_grid_particle {
    float weight;
    float blendedRadius;
    // This is the blended offset from the grid particle position
    frantic::graphics::vector3f blendedOffset;
};

struct zhu_bridson_grid_particle_simd {
    frantic::simd::float_v weight;
    frantic::simd::float_v blendedRadius;

    // This is the blended offset from the grid particle position
    frantic::graphics::vector3t<frantic::simd::float_v> blendedOffset;

    enum width { width = frantic::simd::float_v::static_size };
};

struct zhu_bridson_sparse_data {
    zhu_bridson_sparse_data();

    void reset_for_dense_plane_evaluation(
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::channels::channel_accessor<float>& radiusAccessor, float kernelCompactSupport,
        const frantic::graphics::boundbox3& voxelExtents, const frantic::volumetrics::voxel_coord_system& meshingVCS );

    void reset_for_dense_block_evaluation(
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::channels::channel_accessor<float>& radiusAccessor, float kernelCompactSupport,
        const frantic::graphics::boundbox3& voxelExtents, const frantic::volumetrics::voxel_coord_system& meshingVCS );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> position;
    frantic::channels::channel_accessor<float> radius;
    float kernelCompactSupport;
    float kernelCompactSupportSquared;
    frantic::graphics::boundbox3 voxelExtents;
    frantic::volumetrics::voxel_coord_system meshingVCS; // voxel coord system for meshing

    // output channel data/accessors
    int gridParticleChannel;
    std::vector<char*> channelData;
    std::vector<frantic::channels::channel_general_accessor> channelAccessors;

    // use gridEvaluationSimd instead of gridEvaluation?
    bool useSimd;

    // only used when sparse run data is needed (ie, particle to level set conversion)
    bool collectRunData;
    std::vector<frantic::volumetrics::run_tree> runTrees;
    std::vector<std::vector<char*>> extentData;
    std::vector<std::vector<char>> initChannelData;

    std::vector<detail::zhu_bridson_grid_particle> gridEvaluation;
    std::vector<detail::zhu_bridson_grid_particle_simd> gridEvaluationSimd;

    // A per-voxel flag that keeps track of whether the m_gridEvaluation
    // has been initialized for a voxel.  0 : the voxel is uninitialized
    // and holds undefined data; 1: the voxel holds valid data.
    std::vector<boost::uint8_t> isInitializedVoxel;
    // true if we are using the isInitializedVoxel flags
    bool useIsInitializedVoxel;

    // for multithreading
    bool useMutexes;
    boost::shared_array<tbb::spin_mutex> mutexes;

    std::size_t gridParticleCount;
};
// per-thread storage that is reused between calls to add_edge_isosurface_vertex_to_mesh
struct zhu_bridson_vertex_workspace {
    typedef std::vector<char*> particles_t;
    particles_t particles;
    std::vector<float> sampleWeights;
    std::vector<char> weightedSum;
    frantic::volumetrics::implicitsurface::detail::xyzr_packet_array simdParticles;
};
} // namespace detail

struct GridEvaluationToDensityBody;

class particle_zhu_bridson_is_policy : public particle_is_policy_base<particle_zhu_bridson_is_policy> {
  public:
    typedef detail::zhu_bridson_sparse_data sparse_voxel_corner_density_workspace_t;
    typedef detail::zhu_bridson_vertex_workspace vertex_workspace_t;

  private:
    friend struct GridEvaluationToDensityBody;

    // The Zhu/Bridson parameters
    float m_kernelCompactSupport, m_maximumParticleRadius;
    // An attempt to compensate for extra volume added in concave areas.
    float m_lowDensityTrimmingDensity, m_lowDensityTrimmingStrength;
    // How much to refine the vertex positions
    int m_vertexRefinement;
    float m_vertexRefinementEpsilon;

    BOOST_STATIC_CONSTANT( std::size_t, SIMD_PARTICLE_THRESHOLD = 4 );

    // Derived parameters used for the low density trimming
    float m_kernelAtZero;

    // mutable tbb::atomic<std::size_t> m_vrEvalCount;

    // The voxel grid used to sample the implicit surface.
    frantic::volumetrics::voxel_coord_system m_meshingVCS;

    // Position accessor
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_accessor<float> m_radiusAccessor;

    // The bounding box of the particles
    frantic::graphics::boundbox3 m_particleVCSBounds;

    // When an isosurface location is found, these parameters specify where it is.  In other implementations of the mc
    // policy,
    // the array of filter weights and data indices would be stored to reuse for all the different named channels.
    frantic::graphics::vector3f m_isosurfaceLocationWorld;
    // Also collected is the array of particles that contribute to the density value that that location
    std::vector<char*> m_isosurfaceLocationParticles;

    // Disable copy constructor
    particle_zhu_bridson_is_policy& operator=( const particle_zhu_bridson_is_policy& );

    void get_density_from_grid_evaluation_mt( const frantic::graphics2d::boundrect2& xyExtents, const int z,
                                              const std::vector<boost::uint8_t>& isInitializedVoxel,
                                              const std::vector<detail::zhu_bridson_grid_particle>& gridEvaluation,
                                              float* outVoxelCornerValues ) const;
    void get_density_from_grid_evaluation( const frantic::graphics::boundbox3& xyzExtents,
                                           const std::vector<boost::uint8_t>& isInitializedVoxel,
                                           const std::vector<detail::zhu_bridson_grid_particle>& gridEvaluation,
                                           float* outVoxelCornerValues ) const;
    void get_density_from_grid_evaluation_simd(
        const frantic::graphics::boundbox3& xyzExtents, const std::vector<boost::uint8_t>& isInitializedVoxel,
        const std::vector<detail::zhu_bridson_grid_particle_simd>& gridEvaluation, float* outVoxelCornerValues ) const;
    void
    get_density_from_grid_evaluation_simd_mt( const frantic::graphics2d::boundrect2& xyExtents, const int z,
                                              const std::vector<boost::uint8_t>& isInitializedVoxel,
                                              const std::vector<detail::zhu_bridson_grid_particle_simd>& gridEvaluation,
                                              float* outVoxelCornerValues ) const;
    static void get_density( void* userData, char* p );

    frantic::graphics::vector3f find_edge_isosurface_location(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, std::vector<char*>& particles,
        frantic::volumetrics::implicitsurface::detail::xyzr_packet_array& simdParticles ) const;

    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              const std::vector<char*>& particlesInRange, vertex_workspace_t& workspace,
                              const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const;
    void
    populate_vertex_channels( std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert ) const;

  public:
    // function docs are in the cpp file
    particle_zhu_bridson_is_policy( frantic::particles::particle_grid_tree& particles, float maxParticleRadius,
                                    float effectRadius, float lowDensityTrimmingDensity,
                                    float lowDensityTrimmingStrength, const voxel_coord_system& meshingVCS,
                                    int vertexRefinement );
    ~particle_zhu_bridson_is_policy();

    const frantic::volumetrics::voxel_coord_system& get_voxel_coord_system() const;
    frantic::graphics::boundbox3 get_voxel_bounds() const;
    float get_default_outside_distance() const { return m_kernelCompactSupport; };
    float get_compact_support() const { return m_kernelCompactSupport; };
    void get_sample_weights( const frantic::graphics::vector3f& position, const std::vector<char*>& particles,
                             std::vector<float>& outWeights ) const;
    void get_sample_weights_simd( const frantic::graphics::vector3f& position,
                                  const frantic::volumetrics::implicitsurface::detail::xyzr_packet_array& particles,
                                  std::vector<float>& outWeights ) const;
    void get_affected_blocks( const frantic::graphics::size3f blockSize,
                              std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
                              frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );
    void fill_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                      std::vector<float>& outVoxelCornerValues );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, frantic::volumetrics::rle_plane& rlp,
                                             detail::zhu_bridson_sparse_data& data );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, detail::zhu_bridson_sparse_data& data );
    std::size_t fill_voxel_corner_densities( const frantic::graphics::boundbox3& xyzExtents,
                                             float* voxelCornerDensities, detail::zhu_bridson_sparse_data& data ) const;
    void fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                   std::vector<frantic::tstring>& channelNames, std::vector<char*>& channelData,
                                   frantic::volumetrics::rle_plane& rlp );
    void find_edge_isosurface_location( float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
                                        const frantic::graphics::vector3& voxelCorner1 );
    float get_density( const frantic::graphics::vector3f& worldLocation ) const;
    frantic::graphics::vector3f get_gradient( const frantic::graphics::vector3f& worldLocation, float h );
    const frantic::graphics::vector3f& get_isosurface_location_world() const;

    void set_isosurface_location_world( const frantic::graphics::vector3f& worldLocation );
    void
    add_vertex_to_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels );
    inline frantic::graphics::vector3f
    corner_sample_coord_to_world( const frantic::graphics::vector3& voxelCornerCoord ) const {
        return m_meshingVCS.get_world_voxel_center( voxelCornerCoord );
    }
    void get_interface_widths( float& interfaceWidthInside, float& interfaceWidthOutside ) const;
    int get_exterior_region_code() const;

    void add_edge_isosurface_vertex_to_mesh(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh, detail::zhu_bridson_vertex_workspace& workspace ) const;
    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 const frantic::channels::channel_propagation_policy& cpp ) const;
    void
    populate_vertex_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              vertex_workspace_t& workspace ) const;
    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace,
                              const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const;

    /**
     * Utiliy function.  This amount gets added to the isosurface function value, based on the weight.
     * @param	weight	The weight.
     * @return	float	The amount.
     */
    inline float isosurface_function_compensation( float weight ) const {
        if( weight < m_lowDensityTrimmingDensity ) {
            float x = 1 - weight / m_lowDensityTrimmingDensity;
            return m_lowDensityTrimmingStrength * boost::math::pow<4>( x );
        } else {
            return 0;
        }
    }

    inline frantic::simd::float_v isosurface_function_compensation( const frantic::simd::float_v& weight ) const {
        using frantic::simd::float_v;
        const float_v mask = weight < m_lowDensityTrimmingDensity;
        const float_v x = 1 - weight / m_lowDensityTrimmingDensity;
        const float_v result = m_lowDensityTrimmingStrength * boost::math::pow<4>( x );
        return mask & result;
    }
};

namespace detail {
struct anisotropic_poly6_kernel {

    anisotropic_poly6_kernel() {}

    // for an isotropic kernel :
    // distanceFraction2 = (r/h)^2
    // invCompactSupportVolume = 3 / ( 4 * pi * h^3 ) = 1 / volume of sphere radius h
    static inline float kernel_distance_fraction2( const float distanceFraction2,
                                                   const float invCompactSupportVolume ) {
        if( distanceFraction2 > 1.f ) {
            return 0;
        } else {
            return 105.f / 16 * invCompactSupportVolume * boost::math::pow<3>( 1 - distanceFraction2 );
        }
    }
};

// mass
// from volume..?  assume constant?
// for now assume it is proportional to the radius:
// m = r^3
inline float get_particle_mass( float radius ) {
    const float PI4_3 = 4.1887902f;
    return PI4_3 * boost::math::pow<3>( radius );
}
inline frantic::graphics::vector3f smv33( const float* m, const frantic::graphics::vector3f& v ) {
    return frantic::graphics::vector3f( m[0] * v.x + m[1] * v.y + m[2] * v.z, m[1] * v.x + m[3] * v.y + m[4] * v.z,
                                        m[2] * v.x + m[4] * v.y + m[5] * v.z );
}

struct anisotropic_sparse_data {
    void
    reset_for_dense_plane_evaluation( frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
                                      frantic::channels::channel_accessor<float> radiusAccessor,
                                      frantic::channels::channel_accessor<float> volumeAccessor,
                                      frantic::channels::channel_accessor<float> maxDistanceAccessor,
                                      frantic::channels::channel_accessor<float> invCompactSupportVolumeAccessor,
                                      const frantic::channels::channel_general_accessor& anisotropyAccessor,
                                      const frantic::graphics::boundbox3& voxelExtents,
                                      const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                      float* voxelCornerDensities );

    void
    reset_for_dense_block_evaluation( frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
                                      frantic::channels::channel_accessor<float> radiusAccessor,
                                      frantic::channels::channel_accessor<float> volumeAccessor,
                                      frantic::channels::channel_accessor<float> maxDistanceAccessor,
                                      frantic::channels::channel_accessor<float> invCompactSupportVolumeAccessor,
                                      const frantic::channels::channel_general_accessor& anisotropyAccessor,
                                      const frantic::graphics::boundbox3& voxelExtents,
                                      const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                      float* voxelCornerDensities );

    frantic::channels::channel_accessor<frantic::graphics::vector3f> position;
    frantic::channels::channel_accessor<float> radius;
    frantic::channels::channel_accessor<float> volume;
    frantic::channels::channel_accessor<float> effectRadius;
    frantic::channels::channel_accessor<float> invCompactSupportVolume;
    frantic::channels::channel_general_accessor anisotropy;
    typedef anisotropic_poly6_kernel kernel_t;
    frantic::graphics::boundbox3 voxelExtents;
    frantic::volumetrics::voxel_coord_system meshingVCS; // voxel coord system for meshing

    // output channel data/accessors
    std::vector<char*> channelData;
    std::vector<frantic::channels::channel_general_accessor> channelAccessors;

    // only used when sparse run data is needed (ie, particle to level set conversion)
    bool collectRunData;
    std::vector<frantic::volumetrics::run_tree> runTrees;
    std::vector<std::vector<char*>> extentData;
    std::vector<std::vector<char>> initChannelData;
    int signedDistanceChannel;

    // for multithreading
    bool useMutexes;
    boost::shared_array<tbb::spin_mutex> mutexes;

    std::size_t gridParticleCount;

    std::size_t insideKernelCount;
    std::size_t outsideKernelCount;
};

struct anisotropic_vertex_workspace {
    typedef std::vector<char*> particles_t;
    particles_t particles;
    std::vector<char> weightedSum;
    std::vector<float> sampleWeights;
};
}; // namespace detail

/**
 * See:
 *
 *   Jihun Yu and Greg Turk.  "Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels".  ACM
 * SIGGRAPH Symposium on Computer Animation (2010).
 *   Available Online: http://www.cc.gatech.edu/~turk/my_papers/sph_surfaces.pdf
 */
// template< class KernelType >
class particle_anisotropic_is_policy : public particle_is_policy_base<particle_anisotropic_is_policy> {
  public:
    typedef detail::anisotropic_sparse_data sparse_voxel_corner_density_workspace_t;
    typedef detail::anisotropic_vertex_workspace vertex_workspace_t;

  private:
    // The anisotropic parameters
    float m_maximumEffectRadius, m_implicitThreshold;
    // How much to refine the vertex positions
    int m_vertexRefinement;
    float m_vertexRefinementEpsilon;

    // The voxel grid used to sample the implicit surface.
    frantic::volumetrics::voxel_coord_system m_meshingVCS;

    // Position accessor
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_accessor<float> m_radiusAccessor;
    frantic::channels::channel_accessor<float> m_volumeAccessor;
    frantic::channels::channel_accessor<float> m_invCompactSupportVolumeAccessor;
    frantic::channels::channel_accessor<float> m_maxDistanceAccessor;
    frantic::channels::channel_general_accessor m_anisotropyAccessor;

    // The bounding box of the particles
    frantic::graphics::boundbox3 m_particleVCSBounds;

    // std::vector<frantic::channels::channel_general_accessor> m_inputChannels;

    // When an isosurface location is found, these parameters specify where it is.  In other implementations of the mc
    // policy,
    // the array of filter weights and data indices would be stored to reuse for all the different named channels.
    // frantic::graphics::vector3f m_isosurfaceLocationWorld;
    // Also collected is the array of particles that contribute to the density value at that location
    // std::vector<char*> m_isosurfaceLocationParticles;

    typedef detail::anisotropic_poly6_kernel kernel_t;

    mutable std::size_t m_insideKernelCount;
    mutable std::size_t m_outsideKernelCount;

    frantic::graphics::vector3f find_edge_isosurface_location( float density0, float density1,
                                                               const frantic::graphics::vector3& voxelCorner0,
                                                               const frantic::graphics::vector3& voxelCorner1,
                                                               std::vector<char*>& particles ) const;

    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              const std::vector<char*>& particlesInRange, vertex_workspace_t& workspace,
                              const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const;
    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert ) const;

    static void get_density( void* userData, char* p );

    // Disable copy constructor
    particle_anisotropic_is_policy& operator=( const particle_anisotropic_is_policy& ); // not implemented
  public:
    // function docs are in the cpp file
    particle_anisotropic_is_policy( frantic::particles::particle_grid_tree& particles, float implicitThreshold,
                                    const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement );
    ~particle_anisotropic_is_policy();
    const frantic::volumetrics::voxel_coord_system& get_voxel_coord_system() const;
    float get_default_outside_distance() const {
        return m_implicitThreshold;
    }; // TODO this seems wrong..  are we using this for other policies?
    frantic::graphics::boundbox3 get_voxel_bounds() const;
    float get_compact_support() const {
        return m_maximumEffectRadius /*m_maximumParticleRadius * m_effectRadiusScale*/;
    };
    void get_sample_weights( const frantic::graphics::vector3f& position, const std::vector<char*>& particles,
                             std::vector<float>& outWeights ) const;
    void get_affected_blocks( const frantic::graphics::size3f blockSize,
                              std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
                              frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger );
    void fill_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                      std::vector<float>& outVoxelCornerValues );
    // void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z, float*
    // outVoxelCornerValues, frantic::volumetrics::rle_plane& rlp, detail::union_of_spheres_sparse_data & data );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues,
                                             sparse_voxel_corner_density_workspace_t& data );
    std::size_t fill_voxel_corner_densities( const frantic::graphics::boundbox3& voxelExtents,
                                             float* voxelCornerDensities,
                                             sparse_voxel_corner_density_workspace_t& data ) const;
    /*
    void fill_sparse_channel_data(
        const frantic::graphics2d::boundrect2& xyExtents,
        int z,
        std::vector< frantic::tstring >& channelNames,
        std::vector< char* >& channelData,
        frantic::volumetrics::rle_plane& rlp );
        */
    frantic::graphics::vector3f
    corner_sample_coord_to_world( const frantic::graphics::vector3& voxelCornerCoord ) const;
    void find_edge_isosurface_location( float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
                                        const frantic::graphics::vector3& voxelCorner1 );
    const frantic::graphics::vector3f& get_isosurface_location_world() const;
    void set_isosurface_location_world( const frantic::graphics::vector3f& worldLocation );
    void
    add_vertex_to_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels );
    void get_interface_widths( float& interfaceWidthInside, float& interfaceWidthOutside ) const;
    int get_exterior_region_code() const { return -1; };
    void add_edge_isosurface_vertex_to_mesh(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh, vertex_workspace_t& workspace ) const;
    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 const frantic::channels::channel_propagation_policy& cpp ) const;
    void
    populate_vertex_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              vertex_workspace_t& workspace ) const;
    void
    populate_vertex_channels( const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
                              std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace,
                              const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const;

    // float evaluate_kernel_density( const float distance, const float radius ) const;
    float get_density( const char* particle, const frantic::graphics::vector3f& testPosition ) const;
    float get_density( char** particlesBegin, char** particlesEnd,
                       const frantic::graphics::vector3f& testPosition ) const;
    float get_density( const frantic::graphics::vector3f& worldLocation ) const;
    frantic::graphics::vector3f get_gradient( const frantic::graphics::vector3f& worldLocation, float h );
};
} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic
