// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/bind.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_set.hpp>
#pragma warning( push, 3 )
#pragma warning( disable : 4913 )
#include <boost/thread.hpp>
#pragma warning( pop )

#include <frantic/simd/float_v.hpp>
#include <frantic/simd/int_v.hpp>

#include <frantic/graphics/vector3f.hpp>
#include <frantic/math/eigen.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>
#include <frantic/volumetrics/run_tree.hpp>
//#include <frantic/diagnostics/profiling_manager.hpp>

#pragma warning( push )
#pragma warning( disable : 4512 4100 )
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#pragma warning( pop )

using namespace std;
using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::geometry;
using namespace frantic::particles;
using namespace frantic::simd;
using namespace frantic::channels;

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

// frantic::diagnostics::profiling_manager sparsePM;

// shared_progress_logger_adapter implementation

std::pair<int, std::string>
shared_progress_logger_adapter::progress_logger_adapter_threadproc( boost::function<void( void )>& f ) {
    int result = 4;
    std::string resultString;

    try {
        f();
        result = 0;
    } catch( frantic::logging::progress_cancel_exception& ) {
        result = 1;
    } catch( std::exception& e ) {
        result = 2;
        resultString = e.what();
    } catch( ... ) {
        result = 3;
        resultString = "Unknown error";
    }

    return std::pair<int, std::string>( result, resultString );
}

shared_progress_logger_adapter::shared_progress_logger_adapter( frantic::logging::progress_logger& progressLogger )
    : m_progressLogger( progressLogger ) {}

shared_progress_logger_adapter::~shared_progress_logger_adapter() {}

frantic::volumetrics::implicitsurface::shared_progress_logger_proxy* shared_progress_logger_adapter::get_proxy() {
    return &m_proxy;
}

void shared_progress_logger_adapter::run( boost::function<void( void )>& f ) {
    boost::function<std::pair<int, std::string>( void )> fWrapped(
        boost::bind( progress_logger_adapter_threadproc, boost::ref( f ) ) );

    boost::packaged_task<std::pair<int, std::string>> task( fWrapped );
    boost::unique_future<std::pair<int, std::string>> future = task.get_future();
#if defined _WIN32 && defined _MSC_VER && _MSC_VER >= 1600
    boost::thread thread( std::move( task ) );
#else
    boost::thread thread( boost::move( task ) );
#endif

    bool done = false;
    while( !done ) {
        if( future.timed_wait( boost::posix_time::milliseconds( 100 ) ) ) {
            done = true;
        }

        if( future.has_exception() ) {
            throw std::runtime_error( "Internal Error: Exception in thread." );
        }

        if( future.has_value() ) {
            done = true;
        }

        bool cancel = false;
        try {
            m_progressLogger.update_progress( m_proxy.get_progress() );
        } catch( frantic::logging::progress_cancel_exception& /*e*/ ) {
            cancel = true;
        }

        if( cancel ) {
            m_proxy.cancel();
        }
    }

    std::pair<int, std::string> result = future.get();
    int exitCode = result.first;
    if( exitCode == 0 ) {
        // yay
    } else if( exitCode == 1 ) {
        throw frantic::logging::progress_cancel_exception( "Progress Cancelled" );
    } else {
        throw std::runtime_error( "Internal Error: Exception in thread: " + result.second + " (" +
                                  boost::lexical_cast<std::string>( exitCode ) + ")" );
    }
}

namespace {

/**
 *  Build a new list of index runs, that includes only voxels inside the signedDistanceThreshold.
 *
 * @param[out] outRuns list of runs with signedDistance inside the signedDistanceThreshold.
 * @param inRuns list of runs with defined data in signedDistance.
 * @param offset offset from run indices to signedDistance indices.
 * @param signedDistance the voxels' signed distance values.
 * @param signedDistanceThreshold voxels with signedDistance below this threshold will be kept.
 */
void get_runs_inside_threshold( std::vector<std::pair<int, int>>& outRuns,
                                const std::vector<std::pair<int, int>>& inRuns, int offset, const float* signedDistance,
                                float signedDistanceThreshold ) {
    outRuns.clear();
    outRuns.reserve( inRuns.size() );
    for( std::size_t runIndex = 0; runIndex < inRuns.size(); ++runIndex ) {
        bool outRunActive = false;
        std::pair<int, int> outRun;
        const std::pair<int, int> inRun = inRuns[runIndex];
        for( int i = inRun.first; i <= inRun.second; ++i ) {
            const float phi = signedDistance[i + offset];
            if( phi < signedDistanceThreshold ) {
                if( !outRunActive ) {
                    outRun.first = i;
                    outRunActive = true;
                }
            } else {
                if( outRunActive ) {
                    outRun.second = i - 1;
                    outRunActive = false;

                    outRuns.push_back( outRun );
                }
            }
        }
        if( outRunActive ) {
            outRun.second = inRun.second;
            outRunActive = false;

            outRuns.push_back( outRun );
        }
    }
}

inline int_v create_int_v_run() {
#ifdef FRANTIC_HAS_SSE2
    BOOST_STATIC_ASSERT( int_v::static_size == 4 );
    return int_v( 0, 1, 2, 3 );
#else
    BOOST_STATIC_ASSERT( int_v::static_size == 1 );
    return int_v( 0 );
#endif
}

/**
 * Control the evaluation of a particle's influence on a regular voxel grid.
 * This class determines which voxels are within the particle's sphere of
 * influence, while a derived class computes the actual effect on the grid.
 *
 * This is called "sparse" because it computes the particle's effect only
 * within a sphere surrounding the particle, and not across the entire grid.
 */
template <class Derived>
class sparse_contribution_evaluator {
  public:
    sparse_contribution_evaluator()
        : m_hasParticleOnGrid( false ) {}

    template <class SparseData>
    void evaluate( const char* particle, vector3f particlePosition, float particleEffectRadius, SparseData& data ) {
        using frantic::math::square;

        const float particleEffectRadiusSquared = square( particleEffectRadius );

        // Now we determine the runs of grid points in the meshing VCS that the tree particle will interact with,
        // based on the particle position and the position of the plane that we're interested in.

        const float halfVoxelScalar = data.meshingVCS.voxel_length() / 2;

        // Y Range of interaction for this Z.  particle interaction is determined by radius
        int zMin, zMax;
        zMin = int( floorf(
            data.meshingVCS.get_voxel_z_coord( particlePosition.z - particleEffectRadius + halfVoxelScalar ) ) );
        zMin = std::max<int>( zMin, data.voxelExtents.minimum().z );
        zMax = int( floorf(
            data.meshingVCS.get_voxel_z_coord( particlePosition.z + particleEffectRadius + halfVoxelScalar ) ) );
        zMax = std::min<int>( zMax, data.voxelExtents.maximum().z );

        const int xyArea = data.voxelExtents.xsize() * data.voxelExtents.ysize();

#ifdef FRANTIC_HAS_SSE2
        const float_v voxelLength = data.meshingVCS.voxel_length();
        const float_v halfVoxel = halfVoxelScalar;

        const vector3f origin = data.meshingVCS.world_origin();

        for( int zStepped = zMin; zStepped <= zMax; zStepped += float_v::static_size ) {
            const int_v zVector = int_v( zStepped ) + create_int_v_run();
            int_v zMask( zVector <= zMax );

            const float_v zPlane = float_v( zVector ) * voxelLength + halfVoxel + origin.z;
            const float_v zDistToPlane = particlePosition.z - zPlane;
            const float_v zDistToPlaneSquared = square( zDistToPlane );

            zMask &= int_v( particleEffectRadiusSquared >= zDistToPlaneSquared );
            const float_v particleRadiusOnZSquared = particleEffectRadiusSquared - zDistToPlaneSquared;

            const float_v particleRadiusOnZ = std::sqrt( particleRadiusOnZSquared );

            int_v yMinVector, yMaxVector;
            yMinVector =
                int_v( std::floor( ( particlePosition.y - particleRadiusOnZ + halfVoxel - origin.y ) / voxelLength ) );
            yMaxVector =
                int_v( std::floor( ( particlePosition.y + particleRadiusOnZ + halfVoxel - origin.y ) / voxelLength ) );

            yMinVector = std::max( yMinVector, data.voxelExtents.minimum().y );
            yMaxVector = std::min( yMaxVector, data.voxelExtents.maximum().y );

            for( int zo = 0; zo < float_v::static_size; ++zo ) {
                if( !zMask[zo] ) {
                    continue;
                }

                const int z = zVector[zo];
                const int yMin = yMinVector[zo];
                const int yMax = yMaxVector[zo];
                const float yParticleRadiusOnZSquared = particleRadiusOnZSquared[zo];

                const float worldVoxelCenterZ = data.meshingVCS.get_world_z_voxel_center( z );
                const float particleOffsetZ = particlePosition.z - worldVoxelCenterZ;
                const float distanceSquaredZ = frantic::math::square( particleOffsetZ );

                for( int yStepped = yMin; yStepped <= yMax; yStepped += float_v::static_size ) {
                    const int_v yVector = int_v( yStepped ) + create_int_v_run();

                    const float_v yDiff = particlePosition.y - ( float_v( yVector ) * voxelLength + halfVoxel );
                    const float_v xRangeSquared = yParticleRadiusOnZSquared - square( yDiff );
                    const float_v xRangeSquaredPositive = xRangeSquared >= float_v( 0 );

                    const float_v xRange = std::sqrt( xRangeSquared );

                    int_v xMinVector, xMaxVector;
                    xMinVector =
                        int_v( std::floor( ( particlePosition.x - xRange + halfVoxel - origin.x ) / voxelLength ) );
                    xMaxVector =
                        int_v( std::floor( ( particlePosition.x + xRange + halfVoxel - origin.x ) / voxelLength ) );

                    xMinVector = std::max( xMinVector, data.voxelExtents.minimum().x );
                    xMaxVector = std::min( xMaxVector, data.voxelExtents.maximum().x );

                    int_v yMask = int_v::reinterpret( xRangeSquaredPositive );
                    yMask.and_not( yVector > yMax );
                    yMask.and_not( xMinVector > xMaxVector );

                    for( int yo = 0; yo < float_v::static_size; ++yo ) {
                        if( !yMask[yo] ) {
                            continue;
                        }
                        const int y = yVector[yo];
                        const int xMin = xMinVector[yo];
                        const int xMax = xMaxVector[yo];

                        const float worldVoxelCenterY = data.meshingVCS.get_world_y_voxel_center( y );
                        const float particleOffsetY = particlePosition.y - worldVoxelCenterY;
                        const float distanceSquaredY = frantic::math::square( particleOffsetY );

                        const float distanceSquaredYZ( distanceSquaredY + distanceSquaredZ );

                        evaluate_x_run( particle, particlePosition, xyArea, xMin, xMax, y, z, particleOffsetY,
                                        particleOffsetZ, distanceSquaredYZ, data );
                    }
                }
            }
        }
#else
        const float halfVoxel = halfVoxelScalar;

        for( int z = zMin; z <= zMax; ++z ) {
            float zDistToPlane = particlePosition.z - data.meshingVCS.get_world_z_voxel_center( z );
            float zDistToPlaneSquared = zDistToPlane * zDistToPlane;

            // If the particle is too far away, then it won't interact with the plane.  This can happen
            // when the particle radius is smaller than the max particle radius used to calculate the
            // kernal support and search radius for particles in range of the plane.
            if( zDistToPlaneSquared > particleEffectRadiusSquared )
                continue;

            // radius of the particle on the z plane
            float particleRadiusOnZ = sqrtf( particleEffectRadiusSquared - zDistToPlaneSquared );

            const float worldVoxelCenterZ = data.meshingVCS.get_world_z_voxel_center( z );
            const float particleOffsetZ = particlePosition.z - worldVoxelCenterZ;
            const float distanceSquaredZ = frantic::math::square( particleOffsetZ );

            int yMin, yMax;
            yMin = int(
                floorf( data.meshingVCS.get_voxel_y_coord( particlePosition.y - particleRadiusOnZ + halfVoxel ) ) );
            yMin = std::max<int>( yMin, data.voxelExtents.minimum().y );
            yMax = int(
                floorf( data.meshingVCS.get_voxel_y_coord( particlePosition.y + particleRadiusOnZ + halfVoxel ) ) );
            yMax = std::min<int>( yMax, data.voxelExtents.maximum().y );

            for( int y = yMin; y <= yMax; ++y ) {

                // X range of interaction for this Y
                float yDiff = particlePosition.y - ( data.meshingVCS.get_world_y_coord( float( y ) ) + halfVoxel );
                float xRangeSquared = particleRadiusOnZ * particleRadiusOnZ - ( yDiff * yDiff );
                if( xRangeSquared < 0 )
                    continue;
                float xRange = std::sqrt( xRangeSquared );

                int xMin, xMax;

                xMin = int( floorf( data.meshingVCS.get_voxel_x_coord( particlePosition.x - xRange + halfVoxel ) ) );
                xMin = std::max<int>( xMin, data.voxelExtents.minimum().x );
                xMax = int( floorf( data.meshingVCS.get_voxel_x_coord( particlePosition.x + xRange + halfVoxel ) ) );
                xMax = std::min<int>( xMax, data.voxelExtents.maximum().x );

                if( xMax < xMin )
                    continue;

                const float worldVoxelCenterY = data.meshingVCS.get_world_y_voxel_center( y );
                const float particleOffsetY = particlePosition.y - worldVoxelCenterY;
                const float distanceSquaredY = frantic::math::square( particleOffsetY );

                const float distanceSquaredYZ( distanceSquaredY + distanceSquaredZ );

                evaluate_x_run( particle, particlePosition, xyArea, xMin, xMax, y, z, particleOffsetY, particleOffsetZ,
                                distanceSquaredYZ, data );
            }
        }
#endif
    }

    bool has_particle_on_grid() const { return m_hasParticleOnGrid; }

  private:
    template <class SparseData>
    BOOST_FORCEINLINE void evaluate_x_run( const char* particle, const vector3f& particlePosition, int xyArea, int xMin,
                                           int xMax, int y, int z, float particleOffsetY, float particleOffsetZ,
                                           float distanceSquaredYZ, SparseData& data ) {
        m_hasParticleOnGrid = true;

        int extentNum = y - data.voxelExtents.minimum().y;

        // Populate the grid data and channel data with the appropriate info for this run
        // int dataOffset = extentNum*data.voxelExtents.xsize() - data.voxelExtents.minimum().x;
        int dataOffset = ( z - data.voxelExtents.minimum().z ) * xyArea + extentNum * data.voxelExtents.xsize() -
                         data.voxelExtents.minimum().x;

        // lock down the mutex for this row
        tbb::spin_mutex::scoped_lock lock;
        if( data.useMutexes )
            lock.acquire( data.mutexes[extentNum] );

        implementation().evaluate_x_run( particle, particlePosition, xMin, xMax, particleOffsetY, particleOffsetZ,
                                         distanceSquaredYZ, extentNum, dataOffset, data );
    }

    Derived& implementation() { return *reinterpret_cast<Derived*>( this ); }

    bool m_hasParticleOnGrid;
};

template <class VertexWorkspace>
inline void populate_vertex_channels_impl( const std::vector<channel_general_accessor>& inputChannels_,
                                           std::vector<trimesh3_vertex_channel_general_accessor>& outputChannels_,
                                           std::size_t vertIndex, const std::vector<float>& weights,
                                           const std::vector<char*>& particles, VertexWorkspace& workspace,
                                           const channel_map_weighted_sum& channelMapWeightedSum ) {
    assert( inputChannels_.size() == outputChannels_.size() );
    assert( weights.size() == particles.size() );

    const std::size_t channelCount = outputChannels_.size();
    const channel_general_accessor* inputChannels = inputChannels_.size() > 0 ? &inputChannels_[0] : 0;
    trimesh3_vertex_channel_general_accessor* outputChannels = outputChannels_.size() > 0 ? &outputChannels_[0] : 0;

    const std::size_t particleCount = particles.size();

    // If there are no particles in the range, simply fill the output data with 0's.
    if( particleCount == 0 ) {
        for( std::size_t i = 0; i < channelCount; ++i ) {
            memset( outputChannels[i].data( vertIndex ), 0, outputChannels[i].primitive_size() );
        }
        return;
    }

    workspace.weightedSum.resize( channelMapWeightedSum.structure_size() );

    channelMapWeightedSum.channel_weighted_sum( weights, particles, workspace.weightedSum );

    // Set the data for each trimesh3_vertex_channel_general_accessor
    // Each will be the memory at the accessors offset.
    for( size_t i = 0; i < channelCount; ++i ) {
        size_t offset = inputChannels[i].get_channel_data_pointer( particles[0] ) - particles[0];
        memcpy( outputChannels[i].data( vertIndex ), &workspace.weightedSum[offset],
                inputChannels[i].primitive_size() );
    }
}

void filter_positive_weights( std::vector<char*>& particles, std::vector<float>& weights ) {
    assert( particles.size() == weights.size() );

    for( unsigned i = 0; i < particles.size(); ) {
        if( weights[i] > 0 ) {
            ++i;
        } else {
            weights[i] = weights.back();
            particles[i] = particles.back();
            weights.pop_back();
            particles.pop_back();
        }
    }
}

void load( frantic::volumetrics::implicitsurface::detail::xyzr_packet_array& out, const std::vector<char*>& in,
           const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAcc,
           const frantic::channels::channel_accessor<float>& radiusAcc ) {
    if( in.size() > 0 ) {
        out.load( &in[0], &in[0] + in.size(), positionAcc, radiusAcc );
    } else {
        out.load( 0, 0, positionAcc, radiusAcc );
    }
}

} // anonymous namespace

// details for both Union of Spheres and Metaball
namespace {
/**
 *	Data structures used internally for Union of Spheres and Metaball get_density functions
 */
struct density_user_data {
    bool getGradient;
    frantic::graphics::vector3f worldLocation;
    frantic::channels::channel_accessor<float> radiusAccessor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor;
    float particleRadiusToEffectRadiusScale;
    float h;
    float densities[6];

    density_user_data( const frantic::graphics::vector3f& worldLocation,
                       const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
                       const frantic::channels::channel_accessor<float>& radiusAccessor, float implicitThreshold,
                       float particleRadiusToEffectRadiusScale = 0, float h = 0, bool gradient = false )
        : worldLocation( worldLocation )
        , positionAccessor( positionAccessor )
        , radiusAccessor( radiusAccessor )
        , h( h )
        , particleRadiusToEffectRadiusScale( particleRadiusToEffectRadiusScale )
        , getGradient( gradient ) {

        densities[0] = implicitThreshold;
        if( gradient )
            for( int i = 1; i < 6; i++ )
                densities[i] = implicitThreshold;
    }
};
} // namespace

/*******
 * Functions for the union_of_spheres policy
 *******/

namespace {

struct particle_out_of_range_squared {
    vector3f m_surfaceLocation;
    channel_accessor<vector3f> m_positionAccessor;
    float m_kernelCompactSupportSquared;

    particle_out_of_range_squared( const vector3f& surfaceLocation, const channel_accessor<vector3f>& positionAccessor,
                                   const float distanceSquared )
        : m_surfaceLocation( surfaceLocation )
        , m_positionAccessor( positionAccessor )
        , m_kernelCompactSupportSquared( distanceSquared ) {}
    bool operator()( const char* particle ) {
        const float distanceSquared = vector3f::distance_squared( m_surfaceLocation, m_positionAccessor( particle ) );
        if( distanceSquared >= m_kernelCompactSupportSquared ) {
            return true;
        }
        return false;
    }
};

void get_affected_blocks_impl( const frantic::particles::particle_grid_tree& pgt, float compactSupport,
                               const frantic::graphics::size3f blockSize,
                               std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
                               shared_progress_logger_proxy& progressLogger ) {
    using boost::int32_t;
    using frantic::particles::particle_grid_tree;

    outAffectedBlockCoordinates.clear();

    boost::unordered_set<frantic::graphics::vector3, frantic::graphics::vector3_hasher> blockSet;

    for( particle_grid_tree::const_node_iterator i = pgt.const_nodes_begin(), ie = pgt.const_nodes_end(); i != ie;
         ++i ) {
        frantic::graphics::boundbox3f bounds( i.world_bounds() );
        bounds.expand( compactSupport );
        const int32_t zbegin = static_cast<int32_t>( floor( bounds.minimum().z / blockSize.zsize() ) );
        const int32_t zend = static_cast<int32_t>( floor( bounds.maximum().z / blockSize.zsize() ) );
        const int32_t ybegin = static_cast<int32_t>( floor( bounds.minimum().y / blockSize.ysize() ) );
        const int32_t yend = static_cast<int32_t>( floor( bounds.maximum().y / blockSize.ysize() ) );
        const int32_t xbegin = static_cast<int32_t>( floor( bounds.minimum().x / blockSize.xsize() ) );
        const int32_t xend = static_cast<int32_t>( floor( bounds.maximum().x / blockSize.xsize() ) );
        for( int32_t z = zbegin; z <= zend; ++z ) {
            for( int32_t y = ybegin; y <= yend; ++y ) {
                for( int32_t x = xbegin; x <= xend; ++x ) {
                    blockSet.insert( vector3( x, y, z ) );
                    if( progressLogger.is_cancelled() ) {
                        return;
                    }
                }
            }
        }
    }

    outAffectedBlockCoordinates.insert( outAffectedBlockCoordinates.begin(), blockSet.begin(), blockSet.end() );
}

template <class UnaryFunction>
float find_root_false_position( UnaryFunction& f, float a, float b, float fa, float fb, float tol, int maxIter,
                                int& outIter ) {

    if( fa == 0 ) {
        return a;
    } else if( fb == 0 ) {
        return b;
    } else if( fa * fb >= 0 ) {
        throw std::runtime_error( "find_root_false_position: f(a) and f(b) have the same sign" );
    }

    if( fa > fb ) {
        std::swap( a, b );
        std::swap( fa, fb );
    }

    const float epsilon = 2 * std::numeric_limits<float>::epsilon() * std::max<float>( fabsf( a ), fabsf( b ) ) + tol;

    for( int i = 0; i < maxIter; ++i ) {
        if( fabsf( b - a ) <= epsilon ) {
            break;
        }
        if( fb == 0 ) {
            break;
        }

        const float alpha = fa / ( fa - fb );
        const float x = ( 1 - alpha ) * a + alpha * b;

        const float y = f( x );
        if( y < 0 ) {
            a = x;
            fa = y;
        } else {
            b = x;
            fb = y;
        }

        outIter = i;
    }

    const float alpha = fa / ( fa - fb );
    return ( 1 - alpha ) * a + alpha * b;
}

template <class UnaryFunction>
float find_root_brent( UnaryFunction& f, float a, float b, float fa, float fb, float tol, int maxIter, int& outIter ) {

    outIter = 0;

    if( fa == 0 ) {
        return a;
    } else if( fb == 0 ) {
        return b;
    } else if( fa * fb >= 0 ) {
        throw std::runtime_error( "brent: f(a) and f(b) have the same sign" );
    }

    if( fabsf( fa ) < fabsf( fb ) ) {
        std::swap( a, b );
        std::swap( fa, fb );
    }

    float c = a;
    float fc = fa;

    bool mflag = true;
    float d = 0; // unused on first iteration because mflag is true
    int iter = 0;

    const float epsilon = 2 * std::numeric_limits<float>::epsilon() * std::max<float>( fabsf( a ), fabsf( b ) ) + tol;

    while( fb != 0 && fabsf( b - a ) > epsilon && iter < maxIter ) {

        float s;
        if( fa != fc && fb != fc ) {
            s = a * fb * fc / ( ( fa - fb ) * ( fa - fc ) ) + b * fa * fc / ( ( fb - fa ) * ( fb - fc ) ) +
                c * fa * fb / ( ( fc - fa ) * ( fc - fb ) );
        } else {
            s = b - fb * ( b - a ) / ( fb - fa );
        }
        if( s < 0.25f * ( 3 * a + b ) || s > b || ( mflag && ( fabsf( s - b ) >= 0.5f * fabsf( b - c ) ) ) ||
            ( !mflag && ( fabsf( s - b ) >= 0.5f * fabsf( c - d ) ) ) || ( mflag && ( fabsf( b - c ) < epsilon ) ) ||
            ( !mflag && ( fabsf( c - d ) < epsilon ) ) ) {
            s = 0.5f * ( a + b );
            mflag = true;
        } else {
            mflag = false;
        }

        const float fs = f( s );
        d = c;
        c = b;
        fc = fb;
        if( fa * fs < 0 ) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if( fabsf( fa ) < fabsf( fb ) ) {
            std::swap( a, b );
            std::swap( fa, fb );
        }

        ++iter;
    }

    outIter = iter;
    return b;
}

template <class UnaryFunction>
float find_root( UnaryFunction& f, float a, float b, float fa, float fb, float /*tol*/, int maxIter, int& outIter ) {
    return find_root_false_position( f, a, b, fa, fb, 0, maxIter, outIter );
}

/**
 *	Data structures used internally
 */
struct union_of_spheres_data {
    channel_accessor<vector3f> position;
    channel_accessor<float> radius;
    float particleRadiusToEffectRadiusScale;
};

/**
 *	The union_of_spheres implicit function is subtracted, so that -ve is inside, +ve is outside.
 *	For the union_of_spheres reconstruction, the grid used consists entirely of floats.
 */
void union_of_spheres_contribution_function( void* userData, char* treeParticle, const vector3f& gridParticlePosition,
                                             char* gridParticle ) {
    const union_of_spheres_data* data = reinterpret_cast<union_of_spheres_data*>( userData );

    float distance = vector3f::distance( gridParticlePosition, data->position( treeParticle ) );
    *reinterpret_cast<float*>( gridParticle ) =
        min( distance - data->radius( treeParticle ), *reinterpret_cast<float*>( gridParticle ) );
}

class union_of_spheres_sparse_contribution_evaluator
    : public sparse_contribution_evaluator<union_of_spheres_sparse_contribution_evaluator> {
  public:
    typedef detail::union_of_spheres_sparse_data sparse_data_t;

    BOOST_FORCEINLINE void evaluate_x_run( const char* particle, const vector3f& particlePosition, int xMin, int xMax,
                                           float /*particleOffsetY*/, float /*particleOffsetZ*/,
                                           float distanceSquaredYZ, int extentNum, int dataOffset,
                                           sparse_data_t& data ) {
        if( data.collectRunData ) {
            // Push the run into the appropriate run_tree.  This will initialize the data in if it isn't already.
            std::pair<int, int> run( xMin - data.voxelExtents.minimum().x, xMax - data.voxelExtents.minimum().x - 1 );
            if( run.second >= run.first )
                data.runTrees[extentNum].insert_run( run, data.extentData[extentNum], data.initChannelData );
        }

        for( int x = xMin; x <= xMax; x++ ) {
            int i = x + dataOffset;

            float worldVoxelCenterX = data.meshingVCS.get_world_x_voxel_center( x );
            float particleOffsetX = particlePosition.x - worldVoxelCenterX;
            float distanceSquaredX = frantic::math::square( particleOffsetX );

            float distance = std::sqrt( distanceSquaredX + distanceSquaredYZ ) - data.radius( particle );
            float& gridDistance = ( (float*)( data.channelData[data.signedDistanceChannel] ) )[i];
            if( distance < gridDistance ) {
                gridDistance = distance;
                // memcpy(data.channelData[data.signedDistanceChannel]+i*sizeof(float), &distance, sizeof(float));

                // replace the channel information with this particle data
                for( size_t channelNum = 0; channelNum < data.channelData.size(); ++channelNum ) {
                    if( channelNum != (size_t)data.signedDistanceChannel ) {
                        memcpy( data.channelData[channelNum] + i * data.channelAccessors[channelNum].primitive_size(),
                                data.channelAccessors[channelNum].get_channel_data_pointer( particle ),
                                data.channelAccessors[channelNum].primitive_size() );
                    }
                }
            }
        }
    }
};

/**
 *	Utility function.  This gets run on every particle that will interact with a given
 *	plane of data in the fill_sparse_channel_data function, to determine which grid points
 *	on the plane that the particle interacts with.	The arrays of data that get passed in the user
 *	data struct will get populated accordingly, and the run indexing data for the arrays will be
 *	saved in the vector of runtrees also provided in the user data struct.
 *
 *	@param	userData	void pointer to a union_of_spheres_sparse_data struct
 *	@param	treeParticle	pointer to the data in a tree particle that this function will interact with
 */
void union_of_spheres_sparse_contribution_function( void* userData, char* treeParticle ) {
    detail::union_of_spheres_sparse_data* data = reinterpret_cast<detail::union_of_spheres_sparse_data*>( userData );

    const vector3f treeParticlePosition = data->position( treeParticle );
    const float treeParticleInteractionRadius = data->radius( treeParticle ) * data->particleRadiusToEffectRadiusScale;

    union_of_spheres_sparse_contribution_evaluator f;
    f.evaluate( treeParticle, treeParticlePosition, treeParticleInteractionRadius, *data );

    if( f.has_particle_on_grid() ) {
        ++( data->gridParticleCount );
    }
}

// necessary for multithreading with tbb
class UnionOfSpheresBody {
    detail::union_of_spheres_sparse_data* const m_data;
    std::vector<char*>* const m_particles;

  public:
    void operator()( const tbb::blocked_range<size_t>& r ) const {
        // std::cout << "operator()" << std::endl;
        detail::union_of_spheres_sparse_data* data = m_data;
        std::vector<char*>* particles = m_particles;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            // std::cout << "processing particle " << i << std::endl;
            union_of_spheres_sparse_contribution_function( data, particles->at( i ) );
        }
    }

    UnionOfSpheresBody( detail::union_of_spheres_sparse_data* data, vector<char*>* particles )
        : m_data( data )
        , m_particles( particles ) {}

    UnionOfSpheresBody& operator=( const UnionOfSpheresBody& ) { return *this; } // unimplemented
};

class union_of_spheres_vert_refine_eval {
    frantic::graphics::vector3f m_voxelPosition0;

    int m_solveAxis;

    char** m_particlesBegin;
    char** m_particlesEnd;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_accessor<float> m_radiusAccessor;

    float m_implicitThreshold;

  public:
    union_of_spheres_vert_refine_eval(
        const frantic::graphics::vector3f& voxelCoord0, const int solveAxis, char** particlesBegin, char** particlesEnd,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::channels::channel_accessor<float>& radiusAccessor, float implicitThreshold )
        : m_voxelPosition0( voxelCoord0 )
        , m_solveAxis( solveAxis )
        , m_particlesBegin( particlesBegin )
        , m_particlesEnd( particlesEnd )
        , m_positionAccessor( positionAccessor )
        , m_radiusAccessor( radiusAccessor )
        , m_implicitThreshold( implicitThreshold ) {}

    float operator()( float x ) const {
        vector3f vertTest( m_voxelPosition0 );
        vertTest[m_solveAxis] = x;

        // evaluate the density at this location
        float density = m_implicitThreshold;
        for( char** i = m_particlesBegin; i != m_particlesEnd; ++i ) {
            density = min( vector3f::distance( vertTest, m_positionAccessor( *i ) ) - m_radiusAccessor( *i ), density );
        }
        return density;
    }
};

class union_of_spheres_populate_vertex_channels {
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels;
    frantic::geometry::trimesh3& outMesh;
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels;
    const particle_union_of_spheres_is_policy* isp;

    union_of_spheres_populate_vertex_channels&
    operator=( const union_of_spheres_populate_vertex_channels& ); // not implemented

  public:
    union_of_spheres_populate_vertex_channels(
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
        frantic::geometry::trimesh3& outMesh,
        const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
        const particle_union_of_spheres_is_policy* isp )
        : outputChannels( outputChannels )
        , outMesh( outMesh )
        , inputChannels( inputChannels )
        , isp( isp ) {}

    void operator()( const tbb::blocked_range<size_t>& r ) const {
        detail::union_of_spheres_vertex_workspace workspace;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            isp->populate_vertex_channels( inputChannels, outputChannels, i, outMesh.get_vertex( i ), workspace );
        }
    }
};

} // anonymous namespace

namespace detail {

void union_of_spheres_sparse_data::reset_for_dense_plane_evaluation(
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
    const frantic::channels::channel_accessor<float>& radiusAccessor, const float particleRadiusToEffectRadiusScale_,
    float implicitThreshold, const frantic::graphics::boundbox3& voxelExtents_,
    const frantic::volumetrics::voxel_coord_system& meshingVCS_, float* voxelCornerDensities ) {

    assert( voxelExtents_.zsize() == 1 );

    reset_for_dense_block_evaluation( positionAccessor, radiusAccessor, particleRadiusToEffectRadiusScale_,
                                      implicitThreshold, voxelExtents_, meshingVCS_, voxelCornerDensities );

    this->useMutexes = true;
    this->mutexes.reset( new tbb::spin_mutex[voxelExtents_.ysize()] );
}

void union_of_spheres_sparse_data::reset_for_dense_block_evaluation(
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
    const frantic::channels::channel_accessor<float>& radiusAccessor, const float particleRadiusToEffectRadiusScale_,
    float implicitThreshold, const frantic::graphics::boundbox3& voxelExtents_,
    const frantic::volumetrics::voxel_coord_system& meshingVCS_, float* voxelCornerDensities ) {

    this->position = positionAccessor;
    this->radius = radiusAccessor;
    this->particleRadiusToEffectRadiusScale = particleRadiusToEffectRadiusScale_;
    this->voxelExtents = voxelExtents_;
    this->meshingVCS = meshingVCS_;

    // set up the signed distance channel
    this->collectRunData = false;
    this->signedDistanceChannel = 0;
    this->channelAccessors.clear();
    this->channelAccessors.push_back(
        frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
    this->channelData.clear();
    this->channelData.push_back( (char*)voxelCornerDensities );
    this->initChannelData.clear();
    this->initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
    ( (float*)( &this->initChannelData[0][0] ) )[0] = implicitThreshold;
    this->gridParticleCount = 0;

    // not threaded
    this->mutexes.reset();
    this->useMutexes = false;
}

} // namespace detail

/**
 *	This constructor initializes the policy using the following data
 *	@param	particles	A particle grid tree containing the particle system the policy is to operate on.
 *	@param	maximumParticleRadius	The maximum particle radius in the system.
 *	@param	particleRadiusToEffectRadiusScale	Multiplier for the interaction radius applied to max particle
 *radius
 *	@param	implicitThreshold	Distance threshold, outside of which distance values are ignored
 *	@param	meshingVCS	The coord system that the policy is to work in.
 *	@param	vertRefinement	The number of refinement steps to take when solving for vertex locations.
 */
particle_union_of_spheres_is_policy::particle_union_of_spheres_is_policy(
    particle_grid_tree& particles, float maximumParticleRadius, float particleRadiusToEffectRadiusScale,
    float implicitThreshold, const voxel_coord_system& meshingVCS, int vertexRefinement )
    : particle_is_policy_base<particle_union_of_spheres_is_policy>( particles )
    , m_maximumParticleRadius( maximumParticleRadius )
    , m_particleRadiusToEffectRadiusScale( particleRadiusToEffectRadiusScale )
    , m_implicitThreshold( implicitThreshold )
    , m_meshingVCS( meshingVCS )
    , m_vertexRefinement( vertexRefinement ) {
    m_positionAccessor = m_particles.get_channel_map().get_accessor<vector3f>( _T("Position") );
    m_radiusAccessor = m_particles.get_channel_map().get_accessor<float>( _T("Radius") );

    // Compute the voxel bounds over which to do the implicit surface conversion
    boundbox3f particleWorldBounds = m_particles.compute_particle_bounds();
    particleWorldBounds.expand( maximumParticleRadius * particleRadiusToEffectRadiusScale );
    m_particleVCSBounds = m_meshingVCS.get_voxel_bounds( particleWorldBounds );
    m_particleVCSBounds.expand( 1 );

    m_vertexRefinementEpsilon = 0.001f * m_meshingVCS.voxel_length();
    // m_vrEvalCount = 0;
}

particle_union_of_spheres_is_policy::~particle_union_of_spheres_is_policy() {
    // FF_LOG( debug ) << m_vrEvalCount << std::endl;
}

const voxel_coord_system& particle_union_of_spheres_is_policy::get_voxel_coord_system() const { return m_meshingVCS; }

/**
 *	This returns the XY bounds within which to operate
 */
boundbox3 particle_union_of_spheres_is_policy::get_voxel_bounds() const { return m_particleVCSBounds; }

/**
 *	This returns the interface widths in voxels for the policy.
 *
 *	@param	interfaceWidthInside	assigned the inside interface width
 *	@param	interfaceWidthOutside	assigned the outside interface width
 */
void particle_union_of_spheres_is_policy::get_interface_widths( float& interfaceWidthInside,
                                                                float& interfaceWidthOutside ) const {
    vector3 dist = m_particleVCSBounds.maximum() - m_particleVCSBounds.minimum();
    vector3f distf( float( dist.x ), float( dist.y ), float( dist.z ) );
    interfaceWidthInside = distf.get_magnitude();
    interfaceWidthOutside = 2;
}

/**
 *	Returns the exterior region code for the particle is policy.  Particle is policies have an exterior
 *	region code of -1;
 */
int particle_union_of_spheres_is_policy::get_exterior_region_code() const { return -1; }

void particle_union_of_spheres_is_policy::get_sample_weights( const vector3f& position,
                                                              const std::vector<char*>& particles,
                                                              std::vector<float>& outWeights ) const {
    const size_t particleCount = particles.size();

    outWeights.resize( particleCount );

    if( particleCount == 0 ) {
        return;
    }

    memset( &outWeights[0], 0, particleCount * sizeof( outWeights[0] ) );

    size_t closestParticleIndex = 0;
    float closestDistance = std::numeric_limits<float>::max();

    for( size_t i = 0; i < particleCount; ++i ) {
        const vector3f& particlePosition = m_positionAccessor( particles[i] );
        const float radius = m_radiusAccessor( particles[i] );
        const float distance = vector3f::distance( position, particlePosition ) - radius;

        if( distance < closestDistance ) {
            closestParticleIndex = i;
            closestDistance = distance;
        }
    }

    outWeights[closestParticleIndex] = 1;
}

void particle_union_of_spheres_is_policy::get_affected_blocks(
    const frantic::graphics::size3f blockSize, std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
    shared_progress_logger_proxy& progressLogger ) {
    get_affected_blocks_impl( m_particles, get_compact_support(), blockSize, outAffectedBlockCoordinates,
                              progressLogger );
}

/**
 *	Calculates the minimum distance for a given point to a worldPoint (or the 6 surround worldPoints for gradient
 *calculation). If no distance is as small
 *	as the original implicitThreshold, that implicitThreshold remains as the density.
 *
 *	@param	userData				Expected density_user_data - contains particle information, point
 *locations, and current densities
 *	@param	p						The particle to use for density calculation
 */
void particle_union_of_spheres_is_policy::get_density( void* userData, char* p ) {
    density_user_data* data = reinterpret_cast<density_user_data*>( userData );

    vector3f position = data->positionAccessor( p );
    float radius = data->radiusAccessor( p );
    float distance;

    if( data->getGradient ) {
#ifdef FRANTIC_HAS_SSE2
        //_mm_sfence();

        __m128 densityXY = _mm_set_ps( data->densities[3], data->densities[2], data->densities[1], data->densities[0] );
        __m128 densityZ = _mm_set_ps( 0, 0, data->densities[5], data->densities[4] );

        __m128 delta = _mm_set_ps( 0, 0, -data->h, +data->h );
        __m128 r4 = _mm_set1_ps( radius );
        __m128 accumXY;
        __m128 accumZ;

        __m128 sample4 = _mm_set1_ps( data->worldLocation.x );
        __m128 particle4 = _mm_set1_ps( position.x );
        __m128 diff = _mm_sub_ps( sample4, particle4 );
        __m128 diffWithDelta = _mm_add_ps( diff, delta );
        __m128 diffWithDelta2 = _mm_mul_ps( diffWithDelta, diffWithDelta );

        accumXY = diffWithDelta2;
        accumZ = _mm_shuffle_ps( diffWithDelta2, diffWithDelta2, _MM_SHUFFLE( 1, 0, 3, 2 ) );

        sample4 = _mm_set1_ps( data->worldLocation.y );
        particle4 = _mm_set1_ps( position.y );
        diff = _mm_sub_ps( sample4, particle4 );
        diffWithDelta = _mm_add_ps( diff, delta );
        diffWithDelta2 = _mm_mul_ps( diffWithDelta, diffWithDelta );
        diffWithDelta2 = _mm_shuffle_ps( diffWithDelta2, diffWithDelta2, _MM_SHUFFLE( 1, 0, 3, 2 ) );

        accumXY = _mm_add_ps( accumXY, diffWithDelta2 );
        accumZ = _mm_add_ps( accumZ, diffWithDelta2 );

        sample4 = _mm_set1_ps( data->worldLocation.z );
        particle4 = _mm_set1_ps( position.z );
        diff = _mm_sub_ps( sample4, particle4 );
        diffWithDelta = _mm_add_ps( diff, delta );
        diffWithDelta2 = _mm_mul_ps( diffWithDelta, diffWithDelta );

        accumXY = _mm_add_ps( accumXY, _mm_shuffle_ps( diffWithDelta2, diffWithDelta2, _MM_SHUFFLE( 2, 2, 2, 2 ) ) );
        accumZ = _mm_add_ps( accumZ, diffWithDelta2 );

        densityXY = _mm_min_ps( densityXY, _mm_sub_ps( _mm_sqrt_ps( accumXY ), r4 ) );
        densityZ = _mm_min_ps( densityZ, _mm_sub_ps( _mm_sqrt_ps( accumZ ), r4 ) );

        // TODO : improve this
        union sse_vec {
            __m128 v;
            float data[4];
        };

        sse_vec dataDensityXY;
        sse_vec dataDensityZ;

        dataDensityXY.v = densityXY;
        dataDensityZ.v = densityZ;

        data->densities[0] = dataDensityXY.data[0];
        data->densities[1] = dataDensityXY.data[1];
        data->densities[2] = dataDensityXY.data[2];
        data->densities[3] = dataDensityXY.data[3];
        data->densities[4] = dataDensityZ.data[0];
        data->densities[5] = dataDensityZ.data[1];

//_mm_sfence();
#else
        distance = vector3f::distance( data->worldLocation + vector3f( data->h, 0, 0 ), position );
        data->densities[0] = std::min<float>( data->densities[0], distance - radius );

        distance = vector3f::distance( data->worldLocation - vector3f( data->h, 0, 0 ), position );
        data->densities[1] = std::min<float>( data->densities[1], distance - radius );

        distance = vector3f::distance( data->worldLocation + vector3f( 0, data->h, 0 ), position );
        data->densities[2] = std::min<float>( data->densities[2], distance - radius );

        distance = vector3f::distance( data->worldLocation - vector3f( 0, data->h, 0 ), position );
        data->densities[3] = std::min<float>( data->densities[3], distance - radius );

        distance = vector3f::distance( data->worldLocation + vector3f( 0, 0, data->h ), position );
        data->densities[4] = std::min<float>( data->densities[4], distance - radius );

        distance = vector3f::distance( data->worldLocation - vector3f( 0, 0, data->h ), position );
        data->densities[5] = std::min<float>( data->densities[5], distance - radius );
#endif
    } else {
        distance = vector3f::distance( data->worldLocation, position );
        data->densities[0] = std::min<float>( data->densities[0], distance - radius );
    }
}

/**
 *	Calculates the minimum distance for all points in the worldBounds to the worldLocation
 *
 *	@param	worldLocation			Find the density at this point.
 */
float particle_union_of_spheres_is_policy::get_density( const frantic::graphics::vector3f& worldLocation ) const {
    const float support = m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale;

    density_user_data data( worldLocation, m_positionAccessor, m_radiusAccessor, m_implicitThreshold );
    m_particles.process_particles_in_range( &data, worldLocation, support, get_density );

    return data.densities[0];
}

/**
 *	Calculates the densities at the six points surrounding worldLocation when h is applied in all directions. Then
 *uses those densities to estimate the gradient with standard central-difference formulas.
 *
 *	@param	worldLocation			Find the density at this point.
 *	@param	h						The small step to use.
 */
frantic::graphics::vector3f
particle_union_of_spheres_is_policy::get_gradient( const frantic::graphics::vector3f& worldLocation, float h ) {
    const float support = m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale + h;

    density_user_data data( worldLocation, m_positionAccessor, m_radiusAccessor, m_implicitThreshold, 0, h, true );
    m_particles.process_particles_in_range( &data, worldLocation, support, get_density );

    float x = data.densities[0] - data.densities[1];
    float y = data.densities[2] - data.densities[3];
    float z = data.densities[4] - data.densities[5];

    return vector3f( x, y, z ) / 2 * h;
}

/**
 *	This function takes xyExtents and a z coord for a plane, and populates it densely with
 *	density information as determined by the policy.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A vector that will be populated with data.
 */
void particle_union_of_spheres_is_policy::fill_voxel_corner_densities( const boundrect2& xyExtents, int z,
                                                                       std::vector<float>& outVoxelCornerValues ) {

    boundbox3 xyzExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y, z,
                          z );
    union_of_spheres_data data;
    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.particleRadiusToEffectRadiusScale = m_particleRadiusToEffectRadiusScale;

    // Initialize the output corner values to +m_implicitThreshold
    outVoxelCornerValues.resize( xyExtents.get_area() );
    for( vector<float>::iterator i = outVoxelCornerValues.begin(), iterEnd = outVoxelCornerValues.end(); i != iterEnd;
         ++i )
        *i = m_implicitThreshold;

    m_particles.particle_grid_interactions( xyzExtents, m_meshingVCS, sizeof( float ), &outVoxelCornerValues[0],
                                            &union_of_spheres_contribution_function, &data,
                                            m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	in the provided array.  An rle plane structure is also populated to facilitate sparse access of the
 *	planar data.  The given coordinates are interpreted as being in the coord system already defined
 *	in the policy.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A shared array that will be populated sparsely with data.
 *	@param	rlp			An rle plane indexing structure to be populated to allow sparse access to the
 *						above planar data.
 */
void particle_union_of_spheres_is_policy::fill_sparse_voxel_corner_densities(
    const boundrect2& xyExtents, int z, float* outVoxelCornerValues, rle_plane& rlp,
    detail::union_of_spheres_sparse_data& /*data*/ ) {

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    // set the sparse data structure params
    detail::union_of_spheres_sparse_data data;
    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.particleRadiusToEffectRadiusScale = m_particleRadiusToEffectRadiusScale;
    data.meshingVCS = m_meshingVCS;
    data.voxelExtents = voxelExtents;

    for( int i = 0; i < xyExtents.ysize(); ++i )
        data.runTrees.push_back( run_tree( xyExtents.xsize() ) );

    // set up the signed distance channel
    data.signedDistanceChannel = 0;
    data.channelAccessors.push_back(
        frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
    data.channelData.push_back( (char*)outVoxelCornerValues );
    data.initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
    ( (float*)( &data.initChannelData[0][0] ) )[0] = m_implicitThreshold;

    // set up pointers to the extent data for each channel
    data.extentData = std::vector<std::vector<char*>>( xyExtents.ysize() );
    for( size_t y = 0; y < data.extentData.size(); ++y ) {
        data.extentData[y] = vector<char*>( data.channelData.size() );
        // create a vector of pointers to the current extents in the channel data and the grid particles
        size_t dataOffset = y * data.voxelExtents.xsize();
        for( size_t channelNum = 0; channelNum < data.channelData.size(); channelNum++ )
            data.extentData[y][channelNum] =
                data.channelData[channelNum] + data.channelAccessors[channelNum].primitive_size() * dataOffset;
    }
    data.collectRunData = true;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );
    // data.planeParticles.clear();
    // m_particles.get_particles_in_range( worldSearchBox, m_maximumParticleRadius *
    // m_particleRadiusToEffectRadiusScale, data.planeParticles );

    // set up mutexes and start threads
    data.useMutexes = true;
    data.mutexes.reset( new tbb::spin_mutex[xyExtents.ysize()] );
#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, union_of_spheres_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, union_of_spheres_sparse_contribution_function );
#endif
    /*
    tbb::parallel_for(tbb::blocked_range<size_t>(0,data.planeParticles.size()),
              UnionOfSpheresBody(&data,&data.planeParticles),
              tbb::auto_partitioner());*/

    // Build the run indexing structure from the run trees
    rlp.reset( xyExtents );
    for( int y = 0; y < xyExtents.ysize(); ++y ) {

        vector<pair<int, int>> runs;

        data.runTrees[y].traverse_runs( runs );
        for( int i = 0; i < (int)runs.size(); ++i ) {
            runs[i].first += y * xyExtents.xsize();
            runs[i].second += y * xyExtents.xsize();
        }
        rlp.append_runs_by_extent( runs, y );
    }
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	in the provided array.  An rle plane structure is also populated to facilitate sparse access of the
 *	planar data.  The given coordinates are interpreted as being in the coord system already defined
 *	in the policy.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A shared array that will be populated sparsely with data.
 *
 */
void particle_union_of_spheres_is_policy::fill_sparse_voxel_corner_densities(
    const boundrect2& xyExtents, int z, float* outVoxelCornerValues, detail::union_of_spheres_sparse_data& data ) {

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_plane_evaluation( m_positionAccessor, m_radiusAccessor, m_particleRadiusToEffectRadiusScale,
                                           m_implicitThreshold, voxelExtents, m_meshingVCS, outVoxelCornerValues );

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );

    // have to initialize the plane data, since we don't have an rlp to keep track
    for( int i = 0, ie = xyExtents.get_area(); i < ie; ++i )
        outVoxelCornerValues[i] = m_implicitThreshold;

#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, union_of_spheres_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, union_of_spheres_sparse_contribution_function );
#endif
}

std::size_t
particle_union_of_spheres_is_policy::fill_voxel_corner_densities( const frantic::graphics::boundbox3& voxelExtents,
                                                                  float* voxelCornerDensities,
                                                                  detail::union_of_spheres_sparse_data& data ) const {
    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    const boundbox3f worldExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_block_evaluation( m_positionAccessor, m_radiusAccessor, m_particleRadiusToEffectRadiusScale,
                                           m_implicitThreshold, voxelExtents, m_meshingVCS, voxelCornerDensities );

    // Compute the interactions
    boundbox3f worldSearchBox( worldExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );

    // have to initialize the plane data, since we don't have an rlp to keep track
    for( int i = 0, ie = voxelExtents.get_volume(); i < ie; ++i )
        voxelCornerDensities[i] = m_implicitThreshold;

    m_particles.process_particles_in_bounds( &data, worldSearchBox, union_of_spheres_sparse_contribution_function );

    return data.gridParticleCount;
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	with the channel wise information placing it in the provided char arrays.  An rle plane structure
 *  is also populated to facilitate sparse access of the planar data.  The given coordinates are
 *  interpreted as being in the coord system already defined in the policy.  Only the channels
 *	described in the channel propogation policy argument are propagated, although a signed distance
 *	channel is always propagated.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	channelName	A vector of names for desired channels to be copied.
 *	@param	channelData	A vector of shared arrays that will be populated sparsely with data.
 *	@param	rlp			An rle plane indexing structure to be populated to allow sparse access to the
 *						above planar data.
 */
void particle_union_of_spheres_is_policy::fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents,
                                                                    int z, std::vector<frantic::tstring>& channelNames,
                                                                    std::vector<char*>& channelData,
                                                                    frantic::volumetrics::rle_plane& rlp ) {

    if( channelNames.empty() )
        throw std::runtime_error(
            "particle_union_of_spheres_is_policy::fill_sparse_channel_data() - No channel specified, "
            "channelNames argument was empty." );

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    // Set up the sparse data structure with all the info we need to determine which grid points a particle
    // will interact with
    detail::union_of_spheres_sparse_data data;

    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.particleRadiusToEffectRadiusScale = m_particleRadiusToEffectRadiusScale;
    data.voxelExtents = voxelExtents;
    data.meshingVCS = m_meshingVCS;

    // set up the run tree data
    for( int i = 0; i < xyExtents.ysize(); ++i )
        data.runTrees.push_back( run_tree( xyExtents.xsize() ) );

    // set up the channel data and accessors
    data.signedDistanceChannel = -1;
    for( size_t i = 0; i < channelNames.size(); ++i ) {
        if( channelNames[i] == _T("SignedDistance") ) {
            data.signedDistanceChannel = (int)i;
            data.channelAccessors.push_back(
                frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
            data.channelData.push_back( channelData[i] );
            data.initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
            ( (float*)( &data.initChannelData[i][0] ) )[0] = m_implicitThreshold;
        } else {
            data.channelAccessors.push_back( m_particles.get_channel_map().get_general_accessor( channelNames[i] ) );
            data.channelData.push_back( channelData[i] );
            data.initChannelData.push_back( std::vector<char>( data.channelAccessors[i].primitive_size(), 0 ) );
        }
    }

    // if there wasn't a signed distance channel requested, we'll need a temporary one
    std::unique_ptr<char> signedDistanceData;
    if( data.signedDistanceChannel < 0 ) {
        data.signedDistanceChannel = (int)data.channelAccessors.size();
        data.channelAccessors.push_back(
            frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
        signedDistanceData.reset( new char[sizeof( float ) * xyExtents.get_area()] );
        data.channelData.push_back( signedDistanceData.get() );
        data.initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
        ( (float*)( &data.initChannelData[data.signedDistanceChannel][0] ) )[0] = m_implicitThreshold;
    }

    // set up pointers to the extent data for each channel
    data.extentData = std::vector<std::vector<char*>>( xyExtents.ysize() );
    for( size_t y = 0; y < data.extentData.size(); ++y ) {
        data.extentData[y] = vector<char*>( data.channelData.size() );
        // create a vector of pointers to the current extents in the channel data and the grid particles
        size_t dataOffset = y * data.voxelExtents.xsize();
        for( size_t channelNum = 0; channelNum < data.channelData.size(); channelNum++ )
            data.extentData[y][channelNum] =
                data.channelData[channelNum] + data.channelAccessors[channelNum].primitive_size() * dataOffset;
    }
    data.collectRunData = true;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );

    data.mutexes;
    data.mutexes.reset( new tbb::spin_mutex[xyExtents.ysize()] );
    data.useMutexes = true;
#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, union_of_spheres_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, union_of_spheres_sparse_contribution_function );
#endif

    // Build the run indexing structure from the run trees
    rlp.reset( xyExtents );

    for( int y = 0; y < xyExtents.ysize(); ++y ) {
        vector<pair<int, int>> runs;
        data.runTrees[y].traverse_runs( runs );
        int offset = y * xyExtents.xsize();

        if( runs.size() == 0 ) {
            rlp.append_runs_by_extent( runs, y, -1 );
            continue;
        }

        // I need to pop in undefined runs here, also.  Due to the nature of the particle IS policy,
        // I will only need outside undefined runs, and only between pairs of defined runs.
        vector<pair<int, int>> allRuns;
        allRuns.reserve( 2 * runs.size() );
        vector<int> allRunCodes;
        allRunCodes.reserve( 2 * runs.size() );
        size_t i;
        for( i = 0; i < runs.size() - 1; ++i ) {
            // defined run/code
            allRuns.push_back( std::pair<int, int>( runs[i].first + offset, runs[i].second + offset ) );
            allRunCodes.push_back( 0 );

            // undefined run/code
            allRuns.push_back( std::pair<int, int>( runs[i].second + 1 + offset, runs[i + 1].first - 1 + offset ) );
            allRunCodes.push_back( -1 );
        }
        allRuns.push_back( std::pair<int, int>( runs[i].first + offset, runs[i].second + offset ) );
        allRunCodes.push_back( 0 );

        // pop in the new runs/codes
        rlp.append_runs_by_extent( allRuns, allRunCodes, y );
    }
}

vector3f particle_union_of_spheres_is_policy::corner_sample_coord_to_world( const vector3& voxelCornerCoord ) const {
    return m_meshingVCS.get_world_voxel_center( voxelCornerCoord );
}

frantic::graphics::vector3f particle_union_of_spheres_is_policy::find_edge_isosurface_location(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::vector<char*>& particles ) const {

    particles.clear();

    vector3f worldCorner0 = corner_sample_coord_to_world( voxelCorner0 );
    vector3f worldCorner1 = corner_sample_coord_to_world( voxelCorner1 );
    const int solveAxis = ( worldCorner0 - worldCorner1 ).get_largest_axis();

    vector3f surfaceLocation;

    // Refine the vertex position the requested number of times, to make it more closely approximate the isosurface.
    if( m_vertexRefinement > 0 ) {
        // Make a bounding box containing the extents of this voxel (this will always be an axis-aligned line)
        // This way, we only need to get the nearby particles once, though as a trade off we will get more particles
        // than we need for any individual evaluation.
        boundbox3f vertRange( worldCorner0 );
        vertRange += worldCorner1;
        m_particles.get_particles_in_range( vertRange, m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale,
                                            particles );

        char** particlesBegin = particles.size() > 0 ? &particles[0] : 0;
        char** particlesEnd = particlesBegin + particles.size();

        union_of_spheres_vert_refine_eval solver( worldCorner0, solveAxis, particlesBegin, particlesEnd,
                                                  m_positionAccessor, m_radiusAccessor, m_implicitThreshold );

        int iterCount = 0;
        float refined = find_root( solver, worldCorner0[solveAxis], worldCorner1[solveAxis], density0, density1,
                                   m_vertexRefinementEpsilon, m_vertexRefinement, iterCount );
        // m_vrEvalCount += iterCount;

        surfaceLocation = worldCorner0;
        surfaceLocation[solveAxis] = refined;

        // Filter the particle array to only contain the nearby particles we need
        const float maximumDistanceSquared =
            frantic::math::square( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );
        for( unsigned i = 0; i < m_isosurfaceLocationParticles.size(); ) {
            float distanceSquared = vector3f::distance_squared( surfaceLocation, m_positionAccessor( particles[i] ) );
            if( distanceSquared >= maximumDistanceSquared ) {
                // Remove the particle at the back
                particles[i] = particles.back();
                particles.pop_back();
            } else {
                ++i;
            }
        }
    } else {
        surfaceLocation = worldCorner0;
        float alpha = fabsf( density0 / ( density1 - density0 ) );
        surfaceLocation[solveAxis] = ( 1 - alpha ) * worldCorner0[solveAxis] + alpha * worldCorner1[solveAxis];

        // make sure this array is empty, so the extra channel code knows it needs to generate it.
        // this is done at the start of the function
        // m_isosurfaceLocationParticles.clear();
    }

    return surfaceLocation;
}

/**
 *	This function finds the location of the isosurface between the two voxel corners provided,
 *	and sets the internal state of the policy.
 *	@param	density0		First voxel corner value
 *	@param	density1		Second voxel corner value
 *	@param	voxelCorner0	First voxel corner position
 *	@param	voxelCorner0	Second voxel corner position
 */
void particle_union_of_spheres_is_policy::find_edge_isosurface_location( float density0, float density1,
                                                                         const vector3& voxelCorner0,
                                                                         const vector3& voxelCorner1 ) {
    m_isosurfaceLocationWorld =
        find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1, m_isosurfaceLocationParticles );
}

const vector3f& particle_union_of_spheres_is_policy::get_isosurface_location_world() const {
    return m_isosurfaceLocationWorld;
}

////////
// Functions for dealing with the named channels
////////

/**
 *  Set the internal state of the policy so that it points to the
 * specified location.  This is intended to be called before
 * add_vertex_to_channels if you have manually added a vertex at
 * the specified location, particularly as in generate_disambiguated_faces_for_plane.
 */
void particle_union_of_spheres_is_policy::set_isosurface_location_world(
    const frantic::graphics::vector3f& worldLocation ) {
    m_isosurfaceLocationWorld = worldLocation;
    m_isosurfaceLocationParticles.clear();
}

/**
 *	Given a set of writable trimesh3 accessors, this should add one data element to each of the channels,
 *	which corresponds to the channel names given in prepare_channel_accessors.  The function uses the internal state
 *of the mc policy for the location at which to look up the data
 */
void particle_union_of_spheres_is_policy::add_vertex_to_channels(
    std::vector<geometry::trimesh3_vertex_channel_general_accessor>& outputChannels ) {
    // If necessary, get the particles that are near the sample location
    if( m_isosurfaceLocationParticles.empty() ) {
        m_particles.get_particles_in_range( m_isosurfaceLocationWorld,
                                            m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale,
                                            m_isosurfaceLocationParticles );
    }

    if( outputChannels.size() > 0 ) {
        const std::size_t vertIndex = outputChannels[0].size();

        for( std::size_t i = 0; i < outputChannels.size(); ++i ) {
            if( outputChannels[i].size() <= vertIndex ) {
                outputChannels[i].set_vertex_count( vertIndex + 1 );
            }
        }

        populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, m_isosurfaceLocationWorld,
                                  m_isosurfaceLocationParticles );
    }
}

/**
 *	This function finds the location of the isosurface between the two voxel corners provided,
 *	and sets the internal state of the policy.
 *	@param	density0		First voxel corner value
 *	@param	density1		Second voxel corner value
 *	@param	voxelCorner0	First voxel corner position
 *	@param	voxelCorner0	Second voxel corner position
 *	@param	vertIndex		Index of the vert in the accompanying mesh
 *	@param	outputChannels	Trimesh3 channel accessors for any channel data also to be added
 *	@param	outMesh			The mesh to add the vert to.
 *	@param	meshMutex		A mutex for write locking if the function is used in a multithreaded context.
 */
void particle_union_of_spheres_is_policy::add_edge_isosurface_vertex_to_mesh(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
    frantic::geometry::trimesh3& outMesh, detail::union_of_spheres_vertex_workspace& workspace ) const {

    typedef detail::union_of_spheres_vertex_workspace::particles_t particles_t;
    particles_t& particles = workspace.particles;

    vector3f surfaceLocation =
        find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1, particles );

    float support = m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale;

    // add a vert to the mesh at the surface location
    { // i'm just paranoid about the lock getting released on a manual call, so i scoped it
        // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
        outMesh.vertices_ref()[vertIndex] = surfaceLocation;
    }

    if( outputChannels.empty() )
        return;

    // If we didnt have to fetch particles for the surface solve (ie, 0 vertex refinement)
    // we need to fetch them now.
    if( particles.empty() )
        m_particles.get_particles_in_range( surfaceLocation, support, particles );

    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, surfaceLocation, particles );
}

/**
 *	This function populates all vertex channels specified in the channel propagation
 *	policy using info from the surface policy.
 *	<p>
 *	All other channel info is preserved.
 *
 *	@param	outMesh		Mesh containing vertex and channel info.
 *	@param	cpp			A channel propagation policy.
 */
void particle_union_of_spheres_is_policy::populate_mesh_channels(
    frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_propagation_policy& cpp ) const {
    std::vector<frantic::tstring> previousNames;                // all channels names in the mesh
    std::vector<frantic::tstring> names;                        // all channels that are going to be populated
    std::vector<std::pair<std::size_t, data_type_t>> dataTypes; // types and arity of each channel to be populated
    std::vector<frantic::channels::channel_general_accessor> inputChannels; // channel accessors to channel_map
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>
        outputChannels; // channel accessors to vertex_channel

    // tbb::task_scheduler_init taskSchedulerInit; //please move out.

    // If channels already exist throw error.
    outMesh.get_vertex_channel_names( previousNames );
    for( size_t i = 0; i < previousNames.size(); ++i )
        if( cpp.is_channel_included( previousNames[i] ) )
            throw runtime_error(
                "particle_union_of_spheres_is_policy.populate_mesh_channels: Attempt to populate a channel "
                "that already exists." );

    // Find the names of all channels found in both cpp and the particles.
    get_channel_names( cpp, names );

    // Prepare the input channel accessors and get the data types which
    // will be used to create new channels in the mesh.
    dataTypes.reserve( names.size() );
    inputChannels.reserve( names.size() );
    for( size_t i = 0; i < names.size(); ++i ) {
        inputChannels.push_back( m_particles.get_channel_map().get_general_accessor( names[i] ) );
        dataTypes.push_back( std::make_pair( inputChannels.back().arity(), inputChannels.back().data_type() ) );
    }

    // Create the new vertex and get accessors to those vertex channels.
    for( size_t i = 0; i < names.size(); ++i ) {
        outMesh.add_vertex_channel_raw( names[i], dataTypes[i].first, dataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( names[i] ) );
    }

    // Populate each vertex channel.
    tbb::parallel_for( tbb::blocked_range<size_t>( 0, outMesh.vertex_count() ),
                       union_of_spheres_populate_vertex_channels( outputChannels, outMesh, inputChannels, this ) );
}

/**
 *	This function populates all vertex channels using data from the implicit
 *	surface policy for the input vertex.
 *
 *	@param	inputChannels	Accessors to particle channels used for population.
 *	@param	outputChannels	Accessors to each channel that should be populated.
 *	@param	vert			Location of the vertex.
 *	@param	vertIndex		Index of the vertex.
 */
void particle_union_of_spheres_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert ) const {
    if( inputChannels.size() != outputChannels.size() )
        throw std::runtime_error( "particle_union_of_spheres_is_policy.populate_vertex_channels: inputChannels and "
                                  "outputChannels are not the same size" );

    // Get the particles that are near the sample location
    std::vector<char*> particlesInRange;
    m_particles.get_particles_in_range( vert, m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale,
                                        particlesInRange );

    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, particlesInRange );
}

void particle_union_of_spheres_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, const std::vector<char*>& particlesInRange ) const {
    // Find the closest radial distance
    float closestParticleDistance = boost::numeric::bounds<float>::highest();
    int closestParticle = -1;
    for( unsigned i = 0; i < particlesInRange.size(); ++i ) {
        const float iDistance = vector3f::distance( vert, m_positionAccessor( particlesInRange[i] ) ) -
                                m_radiusAccessor( particlesInRange[i] );
        if( iDistance < closestParticleDistance ) {
            closestParticleDistance = iDistance;
            closestParticle = i;
        }
    }

    // If there are no particles nearby to interpolate, then add a zero value to all the channels
    if( closestParticle < 0 ) {
        for( unsigned i = 0; i < outputChannels.size(); ++i ) {
            memset( outputChannels[i].data( vertIndex ), 0, outputChannels[i].primitive_size() );
        }
        return;
    }

    // copy that particle's data into the output channels
    // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
    for( size_t channelNum = 0; channelNum < outputChannels.size(); ++channelNum ) {
        memcpy( outputChannels[channelNum].data( vertIndex ),
                inputChannels[channelNum].get_channel_data_pointer( particlesInRange[closestParticle] ),
                inputChannels[channelNum].primitive_size() );
    }
}

void particle_union_of_spheres_is_policy::populate_vertex_channels(
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace ) const {
    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, vert, workspace );
}

void particle_union_of_spheres_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace ) const {
    workspace.particles.clear();
    m_particles.get_particles_in_range( vert, m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale,
                                        workspace.particles );
    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, workspace.particles );
}

/*******
 * Functions for the metaball policy
 *******/

namespace {
/**
 * This is the blending function used to create metaball surfaces from particles.
 *
 * @param  distance               The distance from the current sample point to the particle.
 * @param  particleEffectRadius   The radius defining a sphere around the particle within which
 *                                it has an effect on the implicit surface.
 */
inline float metaball_function( float distance, float particleEffectRadius ) {
    if( distance < ( 0.33333333f ) * particleEffectRadius ) {
        return 1 - 3 * math::square( distance / particleEffectRadius );
    } else if( distance < particleEffectRadius ) {
        return ( 3.f / 2.f ) * math::square( 1 - distance / particleEffectRadius );
    } else {
        return 0;
    }
}

/*
float metaball_function( float distance, float particleEffectRadius ) {
  if( distance < particleEffectRadius ) {
    return 1.f - (4.f/9.f) * powf(distance/particleEffectRadius,6.f) + (17.f/9.f) *
powf(distance/particleEffectRadius,4.f) - (22.f/9.f) * powf(distance/particleEffectRadius,2.f);
  }
  else
    return 0;

}
*/
/*
float metaball_function( float distanceSquared, float particleEffectRadiusSquared ) {
using boost::math::pow;
if( distanceSquared < particleEffectRadiusSquared ) {
  const float frac = distanceSquared / particleEffectRadiusSquared;
  return 1.f + frac * ( frac * ( -0.44444444f * frac + 1.88888889f ) - 2.44444444f );
}
return 0;
}
*/

/**
 *	Data structures used internally
 */
struct metaball_data {
    channel_accessor<vector3f> position;
    channel_accessor<float> radius;
    float particleRadiusToEffectRadiusScale;
};

/**
 *	The metaball implicit function is subtracted, so that -ve is inside, +ve is outside.
 *	For the metaball reconstruction, the grid used consists entirely of floats.
 */
void metaball_contribution_function( void* userData, char* treeParticle, const vector3f& gridParticlePosition,
                                     char* gridParticle ) {
    const metaball_data* data = reinterpret_cast<metaball_data*>( userData );

    float distance = vector3f::distance( gridParticlePosition, data->position( treeParticle ) );
    *reinterpret_cast<float*>( gridParticle ) -=
        metaball_function( distance, data->particleRadiusToEffectRadiusScale * data->radius( treeParticle ) );
}

class metaball_sparse_contribution_evaluator
    : public sparse_contribution_evaluator<metaball_sparse_contribution_evaluator> {
  public:
    typedef detail::metaball_sparse_data sparse_data_t;

    metaball_sparse_contribution_evaluator( const char* particle, sparse_data_t& data )
        : m_signedDistanceChannel( data.signedDistanceChannel )
        , m_particleEffectRadius( data.radius( particle ) * data.particleRadiusToEffectRadiusScale ) {}

    BOOST_FORCEINLINE void evaluate_x_run( const char* particle, const vector3f& particlePosition, int xMin, int xMax,
                                           float /*particleOffsetY*/, float /*particleOffsetZ*/,
                                           float distanceSquaredYZ, int extentNum, int dataOffset,
                                           sparse_data_t& data ) {
        if( data.collectRunData ) {
            // Push the run into the appropriate run_tree.  This will initialize the data in if it isn't already.
            std::pair<int, int> run( xMin - data.voxelExtents.minimum().x, xMax - data.voxelExtents.minimum().x - 1 );
            if( run.second >= run.first )
                data.runTrees[extentNum].insert_run( run, data.extentData[extentNum], data.initChannelData );
        }

        for( int x = xMin; x <= xMax; x++ ) {
            int i = x + dataOffset;

            float worldVoxelCenterX = data.meshingVCS.get_world_x_voxel_center( x );
            float particleOffsetX = particlePosition.x - worldVoxelCenterX;
            float distanceSquaredX = frantic::math::square( particleOffsetX );

            float distance = std::sqrt( distanceSquaredX + distanceSquaredYZ );
            const float weight = metaball_function( distance, m_particleEffectRadius );
            ( (float*)( data.channelData[m_signedDistanceChannel] ) )[i] -= weight;

            if( !data.collectRunData )
                continue;

            // add in the effect of this particle to this location in each channel
            size_t channelCount = data.channelData.size();
            for( size_t channelNum = 0; channelNum < channelCount; ++channelNum ) {

                // skip the grid particle channel
                if( channelNum == (size_t)m_signedDistanceChannel )
                    continue;

                // add in the weighted contribution of this particle for the channel
                frantic::channels::channel_general_accessor& ca = data.channelAccessors[channelNum];
                ca.weighted_increment( weight, ca.get_channel_data_pointer( particle ),
                                       data.channelData[channelNum] + i * ca.primitive_size() );
            }
        }
    }

  private:
    int m_signedDistanceChannel;
    float m_particleEffectRadius;
};

/**
 *	Utility function.  This gets run on every particle that will interact with a given
 *	plane of data in the fill_sparse_channel_data function, to determine which grid points
 *	on the plane that the particle interacts with.	The arrays of data that get passed in the user
 *	data struct will get populated accordingly, and the run indexing data for the arrays will be
 *	saved in the vector of runtrees also provided in the user data struct.
 *
 *	@param	userData	void pointer to a metaball_sparse_data struct
 *	@param	treeParticle	pointer to the data in a tree particle that this function will interact with
 */
void metaball_sparse_contribution_function( void* userData, char* treeParticle ) {
    detail::metaball_sparse_data* data = reinterpret_cast<detail::metaball_sparse_data*>( userData );

    const vector3f treeParticlePosition = data->position( treeParticle );
    const float treeParticleInteractionRadius = data->radius( treeParticle ) * data->particleRadiusToEffectRadiusScale;

    metaball_sparse_contribution_evaluator f( treeParticle, *data );
    f.evaluate( treeParticle, treeParticlePosition, treeParticleInteractionRadius, *data );

    if( f.has_particle_on_grid() ) {
        ++( data->gridParticleCount );
    }
}

// necessary for multithreading with tbb
class MetaballsBody {
    detail::metaball_sparse_data* const m_data;
    std::vector<char*>* const m_particles;

  public:
    void operator()( const tbb::blocked_range<size_t>& r ) const {
        // std::cout << "operator()" << std::endl;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            // std::cout << "processing particle " << i << std::endl;
            metaball_sparse_contribution_function( m_data, m_particles->at( i ) );
        }
    }

    MetaballsBody( detail::metaball_sparse_data* data, vector<char*>* particles )
        : m_data( data )
        , m_particles( particles ) {}

    MetaballsBody& operator=( const MetaballsBody& ) { return *this; } // unimplemented
};

class metaball_vert_refine_eval {
    frantic::graphics::vector3f m_voxelPosition0;

    int m_solveAxis;

    char** m_particlesBegin;
    char** m_particlesEnd;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_accessor<float> m_radiusAccessor;

    float m_implicitThreshold;
    float m_particleRadiusToEffectRadiusScale;

  public:
    metaball_vert_refine_eval( const frantic::graphics::vector3f& voxelCoord0, const int solveAxis,
                               char** particlesBegin, char** particlesEnd,
                               const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
                               const frantic::channels::channel_accessor<float>& radiusAccessor,
                               float implicitThreshold, float particleRadiusToEffectRadiusScale )
        : m_voxelPosition0( voxelCoord0 )
        , m_solveAxis( solveAxis )
        , m_particlesBegin( particlesBegin )
        , m_particlesEnd( particlesEnd )
        , m_positionAccessor( positionAccessor )
        , m_radiusAccessor( radiusAccessor )
        , m_implicitThreshold( implicitThreshold )
        , m_particleRadiusToEffectRadiusScale( particleRadiusToEffectRadiusScale ) {}

    float operator()( float x ) const {
        vector3f vertTest( m_voxelPosition0 );
        vertTest[m_solveAxis] = x;

        // evaluate the density at this location
        float density = m_implicitThreshold;
        for( char** i = m_particlesBegin; i != m_particlesEnd; ++i ) {
            float distance = vector3f::distance( vertTest, m_positionAccessor( *i ) );
            density -= metaball_function( distance, m_radiusAccessor( *i ) * m_particleRadiusToEffectRadiusScale );
        }
        return density;
    }
};

class metaball_populate_vertex_channels {
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels;
    frantic::geometry::trimesh3& outMesh;
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels;
    const particle_metaball_is_policy* isp;
    const channel_map_weighted_sum& channelMapWeightedSum;

    metaball_populate_vertex_channels& operator=( const metaball_populate_vertex_channels& ); // not implemented

  public:
    metaball_populate_vertex_channels(
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
        frantic::geometry::trimesh3& outMesh,
        const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
        const particle_metaball_is_policy* isp, const channel_map_weighted_sum& channelMapWeightedSum )
        : outputChannels( outputChannels )
        , outMesh( outMesh )
        , inputChannels( inputChannels )
        , isp( isp )
        , channelMapWeightedSum( channelMapWeightedSum ) {}

    void operator()( const tbb::blocked_range<size_t>& r ) const {
        detail::metaball_vertex_workspace workspace;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            isp->populate_vertex_channels( inputChannels, outputChannels, i, outMesh.get_vertex( i ), workspace,
                                           channelMapWeightedSum );
        }
    }
};

} // anonymous namespace

namespace detail {

void metaball_sparse_data::reset_for_dense_plane_evaluation(
    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
    frantic::channels::channel_accessor<float> radiusAccessor, float particleRadiusToEffectRadiusScale_,
    const frantic::graphics::boundbox3& voxelExtents_, const frantic::volumetrics::voxel_coord_system& meshingVCS_,
    float* voxelCornerDensities_ ) {

    assert( voxelExtents_.zsize() == 1 );

    reset_for_dense_block_evaluation( positionAccessor, radiusAccessor, particleRadiusToEffectRadiusScale_,
                                      voxelExtents_, meshingVCS_, voxelCornerDensities_ );

    // set up mutexes
    this->useMutexes = true;
    this->mutexes.reset( new tbb::spin_mutex[voxelExtents_.ysize()] );
}

void metaball_sparse_data::reset_for_dense_block_evaluation(
    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
    frantic::channels::channel_accessor<float> radiusAccessor, float particleRadiusToEffectRadiusScale_,
    const frantic::graphics::boundbox3& voxelExtents_, const frantic::volumetrics::voxel_coord_system& meshingVCS_,
    float* voxelCornerDensities_ ) {

    // Set up the sparse data structure with all the info we need to determine which grid points a particle
    // will interact with
    this->position = positionAccessor;
    this->radius = radiusAccessor;
    this->particleRadiusToEffectRadiusScale = particleRadiusToEffectRadiusScale_;
    this->voxelExtents = voxelExtents_;
    this->meshingVCS = meshingVCS_;

    // set up the signed distance channel
    this->signedDistanceChannel = 0;
    this->channelAccessors.clear();
    this->channelAccessors.push_back(
        frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
    this->channelData.clear();
    this->channelData.push_back( (char*)voxelCornerDensities_ );
    this->initChannelData.clear();
    this->initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
    this->gridParticleCount = 0;

    this->collectRunData = false;

    // not multithreaded
    this->mutexes.reset();
    this->useMutexes = false;
}

} // namespace detail

/**
 *	This constructor initializes the policy using the following data
 *	@param	particles	A particle grid tree containing the particle system the policy is to operate on.
 *	@param	maximumParticleRadius	The maximum particle radius in the system.
 *	@param	particleRadiusToEffectRadiusScale	Required paramater for metaballs
 *	@param	implicitThreshold	Required paramater for metaballs
 *	@param	meshingVCS	The coord system that the policy is to work in.
 *	@param	vertRefinement	The number of refinement steps to take when solving for vertex locations.
 */
particle_metaball_is_policy::particle_metaball_is_policy( particle_grid_tree& particles, float maximumParticleRadius,
                                                          float particleRadiusToEffectRadiusScale,
                                                          float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                          int vertexRefinement )
    : particle_is_policy_base<particle_metaball_is_policy>( particles )
    , m_maximumParticleRadius( maximumParticleRadius )
    , m_particleRadiusToEffectRadiusScale( particleRadiusToEffectRadiusScale )
    , m_implicitThreshold( implicitThreshold )
    , m_meshingVCS( meshingVCS )
    , m_vertexRefinement( vertexRefinement ) {
    m_positionAccessor = m_particles.get_channel_map().get_accessor<vector3f>( _T("Position") );
    m_radiusAccessor = m_particles.get_channel_map().get_accessor<float>( _T("Radius") );

    // Compute the voxel bounds over which to do the implicit surface conversion
    boundbox3f particleWorldBounds = m_particles.compute_particle_bounds();
    particleWorldBounds.expand( maximumParticleRadius * particleRadiusToEffectRadiusScale );
    m_particleVCSBounds = m_meshingVCS.get_voxel_bounds( particleWorldBounds );
    m_particleVCSBounds.expand( 1 );

    m_vertexRefinementEpsilon = 0.001f * m_meshingVCS.voxel_length();

    // m_vrEvalCount = 0;
}

particle_metaball_is_policy::~particle_metaball_is_policy() {
    // FF_LOG( debug ) << m_vrEvalCount << std::endl;
}

const voxel_coord_system& particle_metaball_is_policy::get_voxel_coord_system() const { return m_meshingVCS; }

/**
 *	This returns the XY bounds within which to operate
 */
boundbox3 particle_metaball_is_policy::get_voxel_bounds() const { return m_particleVCSBounds; }

/**
 *	This returns the interface widths in voxels for the policy.
 *
 *	@param	interfaceWidthInside	assigned the inside interface width
 *	@param	interfaceWidthOutside	assigned the outside interface width
 */
void particle_metaball_is_policy::get_interface_widths( float& interfaceWidthInside,
                                                        float& interfaceWidthOutside ) const {
    vector3 dist = m_particleVCSBounds.maximum() - m_particleVCSBounds.minimum();
    vector3f distf( float( dist.x ), float( dist.y ), float( dist.z ) );
    interfaceWidthInside = distf.get_magnitude();
    interfaceWidthOutside = 2; // TODO: this is going to depend on the minimal particle radius
}

/**
 *	Returns the exterior region code for the particle is policy.  Particle is policies have an exterior
 *	region code of -1;
 */
int particle_metaball_is_policy::get_exterior_region_code() const { return -1; }

void particle_metaball_is_policy::get_sample_weights( const vector3f& position, const std::vector<char*>& particles,
                                                      std::vector<float>& outWeights ) const {
    const size_t particleCount = particles.size();

    outWeights.resize( particleCount );

    if( particleCount == 0 ) {
        return;
    }

    // Compute the sample weights for blending between the different values
    float weightSum = 0;
    for( size_t i = 0; i < particleCount; ++i ) {
        const vector3f& particlePosition = m_positionAccessor( particles[i] );
        const float distance = vector3f::distance( position, particlePosition );

        const float radius = m_radiusAccessor( particles[i] );
        const float effectRadius = radius * m_particleRadiusToEffectRadiusScale;

        const float weight = metaball_function( distance, effectRadius );

        outWeights[i] = weight;
        weightSum += weight;
    }

    // TODO: this is a hack to avoid infs in the output channel when there is
    // no particle in range.  What should we do in this case?
    if( weightSum == 0 ) {
        size_t nearestParticleIndex = 0;
        float nearestParticleDistance2 = std::numeric_limits<float>::max();
        for( size_t i = 0; i < particleCount; ++i ) {
            const vector3f& particlePosition = m_positionAccessor( particles[i] );
            const float distance2 = vector3f::distance_squared( position, particlePosition );
            if( distance2 < nearestParticleDistance2 ) {
                nearestParticleIndex = i;
                nearestParticleDistance2 = distance2;
            }
            outWeights[i] = 0.f;
        }
        outWeights[nearestParticleIndex] = 1.f;
        weightSum = 1.f;
    }

    const float invWeightSum = 1.f / weightSum;

    // Normalize the sample weights
    for( size_t i = 0; i < particleCount; ++i ) {
        outWeights[i] *= invWeightSum;
    }
}

void particle_metaball_is_policy::get_affected_blocks(
    const frantic::graphics::size3f blockSize, std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
    shared_progress_logger_proxy& progressLogger ) {
    get_affected_blocks_impl( m_particles, get_compact_support(), blockSize, outAffectedBlockCoordinates,
                              progressLogger );
}

/**
 *	This function takes xyExtents and a z coord for a plane, and populates it densely with
 *	density information as determined by the policy.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A vector that will be populated with data.
 */
void particle_metaball_is_policy::fill_voxel_corner_densities( const boundrect2& xyExtents, int z,
                                                               std::vector<float>& outVoxelCornerValues ) {

    boundbox3 xyzExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y, z,
                          z );
    metaball_data data;
    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.particleRadiusToEffectRadiusScale = m_particleRadiusToEffectRadiusScale;

    // Initialize the output corner values to +m_implicitThreshold
    outVoxelCornerValues.resize( xyExtents.get_area() );
    for( vector<float>::iterator i = outVoxelCornerValues.begin(), iterEnd = outVoxelCornerValues.end(); i != iterEnd;
         ++i )
        *i = m_implicitThreshold;

    m_particles.particle_grid_interactions( xyzExtents, m_meshingVCS, sizeof( float ), &outVoxelCornerValues[0],
                                            &metaball_contribution_function, &data,
                                            m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	in the provided array.  An rle plane structure is also populated to facilitate sparse access of the
 *	planar data.  The given coordinates are interpreted as being in the coord system already defined
 *	in the policy.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A shared array that will be populated sparsely with data.
 *	@param	rlp			An rle plane indexing structure to be populated to allow sparse access to the
 *						above planar data.
 */
void particle_metaball_is_policy::fill_sparse_voxel_corner_densities( const boundrect2& xyExtents, int z,
                                                                      float* outVoxelCornerValues, rle_plane& rlp,
                                                                      detail::metaball_sparse_data& /*data*/ ) {

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    detail::metaball_sparse_data data;
    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.particleRadiusToEffectRadiusScale = m_particleRadiusToEffectRadiusScale;
    data.meshingVCS = m_meshingVCS;
    data.voxelExtents = voxelExtents;

    for( int i = 0; i < xyExtents.ysize(); ++i )
        data.runTrees.push_back( run_tree( xyExtents.xsize() ) );

    // set up the signed distance channel
    data.signedDistanceChannel = 0;
    data.channelAccessors.push_back(
        frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
    data.channelData.push_back( (char*)outVoxelCornerValues );
    data.initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
    ( (float*)( &data.initChannelData[0][0] ) )[0] = m_implicitThreshold;

    // set up pointers to the extent data for each channel
    data.extentData = std::vector<std::vector<char*>>( xyExtents.ysize() );
    for( size_t y = 0; y < data.extentData.size(); ++y ) {
        data.extentData[y] = vector<char*>( data.channelData.size() );
        // create a vector of pointers to the current extents in the channel data and the grid particles
        size_t dataOffset = y * data.voxelExtents.xsize();
        for( size_t channelNum = 0; channelNum < data.channelData.size(); channelNum++ )
            data.extentData[y][channelNum] =
                data.channelData[channelNum] + data.channelAccessors[channelNum].primitive_size() * dataOffset;
    }
    data.collectRunData = true;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );

    // set up mutexes and start threads
    data.mutexes.reset( new tbb::spin_mutex[xyExtents.ysize()] );
    data.useMutexes = true;
#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, metaball_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, metaball_sparse_contribution_function );
#endif

    // Build the run indexing structure from the run trees
    rlp.reset( xyExtents );
    for( int y = 0; y < xyExtents.ysize(); ++y ) {

        vector<pair<int, int>> runs;

        data.runTrees[y].traverse_runs( runs );
        for( int i = 0; i < (int)runs.size(); ++i ) {
            runs[i].first += y * xyExtents.xsize();
            runs[i].second += y * xyExtents.xsize();
        }
        rlp.append_runs_by_extent( runs, y );
    }
}

std::size_t particle_metaball_is_policy::fill_voxel_corner_densities( const frantic::graphics::boundbox3& voxelExtents,
                                                                      float* voxelCornerDensities,
                                                                      detail::metaball_sparse_data& data ) const {
    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    const boundbox3f worldExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_block_evaluation( m_positionAccessor, m_radiusAccessor, m_particleRadiusToEffectRadiusScale,
                                           voxelExtents, m_meshingVCS, voxelCornerDensities );

    for( int i = 0, ie = voxelExtents.get_volume(); i != ie; ++i )
        voxelCornerDensities[i] = m_implicitThreshold;

    // Compute the interactions
    boundbox3f worldSearchBox( worldExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );

    m_particles.process_particles_in_bounds( &data, worldSearchBox, metaball_sparse_contribution_function );

    return data.gridParticleCount;
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	in the provided array.  By sparsely, we mean that it calculates only those interactions
 *  as dictated by each particle and its radii, rather than by each dense grid point.
 *  The given coordinates are interpreted as being in the coord system already defined
 *	in the policy.
 *
 *	This function has been multithreaded using tbb.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A shared array that will be populated sparsely with data.
 */
void particle_metaball_is_policy::fill_sparse_voxel_corner_densities( const boundrect2& xyExtents, int z,
                                                                      float* outVoxelCornerValues,
                                                                      detail::metaball_sparse_data& data ) {

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_plane_evaluation( m_positionAccessor, m_radiusAccessor, m_particleRadiusToEffectRadiusScale,
                                           voxelExtents, m_meshingVCS, outVoxelCornerValues );

    for( int i = 0, ie = xyExtents.get_area(); i != ie; ++i )
        outVoxelCornerValues[i] = m_implicitThreshold;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );

#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, metaball_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, metaball_sparse_contribution_function );
#endif
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	with the channel wise information using the in the provided char arrays.  An rle plane structure
 *  is also populated to facilitate sparse access of the planar data.  The given coordinates are
 *  interpreted as being in the coord system already defined in the policy.  Only channels described
 *	in the channel propagation policy argument will be propagated.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	channelData	A vector of shared arrays that will be populated sparsely with data.
 *	@param	rlp			An rle plane indexing structure to be populated to allow sparse access to the
 *						above planar data.
 *	@param	cpp			Channel propagation policy for which channels to propagate from the particles
 *						to the level set.
 */
void particle_metaball_is_policy::fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                                            std::vector<frantic::tstring>& channelNames,
                                                            std::vector<char*>& channelData,
                                                            frantic::volumetrics::rle_plane& rlp ) {

    if( channelNames.empty() )
        throw std::runtime_error( "particle_metaball_is_policy::fill_sparse_channel_data() - No channel specified, "
                                  "channelNames argument was empty." );

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    // Set up the sparse data structure with all the info we need to determine which grid points a particle
    // will interact with
    detail::metaball_sparse_data data;

    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.particleRadiusToEffectRadiusScale = m_particleRadiusToEffectRadiusScale;
    data.voxelExtents = voxelExtents;
    data.meshingVCS = m_meshingVCS;

    // set up the run tree data
    for( int i = 0; i < xyExtents.ysize(); ++i )
        data.runTrees.push_back( run_tree( xyExtents.xsize() ) );

    // set up the channel data and accessors
    data.signedDistanceChannel = -1;
    for( size_t i = 0; i < channelNames.size(); ++i ) {
        if( channelNames[i] == _T("SignedDistance") ) {
            data.signedDistanceChannel = (int)i;
            data.channelAccessors.push_back(
                frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
            data.channelData.push_back( channelData[i] );
            data.initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
            ( (float*)( &data.initChannelData[i][0] ) )[0] = m_implicitThreshold;
        } else {
            data.channelAccessors.push_back( m_particles.get_channel_map().get_general_accessor( channelNames[i] ) );
            data.channelData.push_back( channelData[i] );
            data.initChannelData.push_back( std::vector<char>( data.channelAccessors[i].primitive_size(), 0 ) );
        }
    }

    // if there wasn't a signed distance channel requested, we'll need a temporary one
    std::unique_ptr<char> signedDistanceData;
    if( data.signedDistanceChannel < 0 ) {
        data.signedDistanceChannel = (int)data.channelAccessors.size();
        data.channelAccessors.push_back(
            frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
        signedDistanceData.reset( new char[sizeof( float ) * xyExtents.get_area()] );
        data.channelData.push_back( signedDistanceData.get() );
        data.initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
        ( (float*)( &data.initChannelData[data.signedDistanceChannel][0] ) )[0] = m_implicitThreshold;
    }

    // set up pointers to the extent data for each channel
    data.extentData = std::vector<std::vector<char*>>( xyExtents.ysize() );
    for( size_t y = 0; y < data.extentData.size(); ++y ) {
        data.extentData[y] = vector<char*>( data.channelData.size() );
        // create a vector of pointers to the current extents in the channel data and the grid particles
        size_t dataOffset = y * data.voxelExtents.xsize();
        for( size_t channelNum = 0; channelNum < data.channelData.size(); channelNum++ )
            data.extentData[y][channelNum] =
                data.channelData[channelNum] + data.channelAccessors[channelNum].primitive_size() * dataOffset;
    }
    data.collectRunData = true;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale );

    // set up mutexes and start threads
    data.mutexes.reset( new tbb::spin_mutex[xyExtents.ysize()] );
    data.useMutexes = true;
// tbb::task_scheduler_init taskSchedulerInit; //please move out.
#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, metaball_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, metaball_sparse_contribution_function );
#endif

    // Build the run indexing structure from the run trees
    rlp.reset( xyExtents );

    for( int y = 0; y < xyExtents.ysize(); ++y ) {
        vector<pair<int, int>> dataRuns;
        data.runTrees[y].traverse_runs( dataRuns );
        int offset = y * xyExtents.xsize();

        vector<pair<int, int>> runs;
        get_runs_inside_threshold( runs, dataRuns, offset,
                                   reinterpret_cast<float*>( data.channelData[data.signedDistanceChannel] ),
                                   m_implicitThreshold );

        if( runs.size() == 0 ) {
            rlp.append_runs_by_extent( runs, y, -1 );
            continue;
        }

        // I need to pop in undefined runs here, also.  Due to the nature of the particle IS policy,
        // I will only need outside undefined runs, and only between pairs of defined runs.
        vector<pair<int, int>> allRuns;
        allRuns.reserve( 2 * runs.size() );
        vector<int> allRunCodes;
        allRunCodes.reserve( 2 * runs.size() );
        size_t i;
        for( i = 0; i < runs.size() - 1; ++i ) {
            // defined run/code
            allRuns.push_back( std::pair<int, int>( runs[i].first + offset, runs[i].second + offset ) );
            allRunCodes.push_back( 0 );

            // undefined run/code
            allRuns.push_back( std::pair<int, int>( runs[i].second + 1 + offset, runs[i + 1].first - 1 + offset ) );
            allRunCodes.push_back( -1 );
        }
        allRuns.push_back( std::pair<int, int>( runs[i].first + offset, runs[i].second + offset ) );
        allRunCodes.push_back( 0 );

        // pop in the new runs/codes
        rlp.append_runs_by_extent( allRuns, allRunCodes, y );
    }

    // Use it to traverse the grid point values and calculate the final channel values
    std::vector<float> weights( 1 );
    vector<char*> inData( 1 );

    if( !rlp.is_empty() ) {
        for( int y = 0; y < rlp.get_xy_extents().ysize(); y++ ) {
            int firstRun = rlp.get_y_extent_first_run( y );
            if( firstRun >= 0 ) {
                int lastRun = rlp.get_y_extent_last_run( y );
                for( int run = firstRun; run <= lastRun; ++run ) {
                    // skip undefined runs
                    if( rlp.get_run_code( run ) < 0 )
                        continue;
                    for( int index = rlp[run].first; index <= rlp[run].second; ++index ) {
                        weights[0] = -1.f / ( ( (float*)( data.channelData[data.signedDistanceChannel] ) )[index] -
                                              m_implicitThreshold );
                        for( size_t channelNum = 1; channelNum < channelData.size(); ++channelNum ) {
                            inData[0] = data.channelData[channelNum] +
                                        index * data.channelAccessors[channelNum].primitive_size();
                            data.channelAccessors[channelNum].weighted_sum( weights, inData, inData[0] );
                        }
                    }
                }
            }
        }
    }
}

vector3f particle_metaball_is_policy::corner_sample_coord_to_world( const vector3& voxelCornerCoord ) const {
    return m_meshingVCS.get_world_voxel_center( voxelCornerCoord );
}

/**
 *	Calculates the metaball densities for the userData in relation to p.
 *
 *	@param	userData				Expected struct density_user_data. contains particle information,
 *point locations, and current densities
 *	@param	p						The particle to use for density calculation
 */
void particle_metaball_is_policy::get_density( void* userData, char* p ) {
    density_user_data* data = reinterpret_cast<density_user_data*>( userData );

    const frantic::graphics::vector3f& position = data->positionAccessor( p );
    float radius = data->radiusAccessor( p );

    if( data->getGradient ) {

        float distance = vector3f::distance( data->worldLocation + vector3f( data->h, 0, 0 ), position );
        data->densities[0] -= metaball_function( distance, radius * data->particleRadiusToEffectRadiusScale );

        distance = vector3f::distance( data->worldLocation - vector3f( data->h, 0, 0 ), position );
        data->densities[1] -= metaball_function( distance, radius * data->particleRadiusToEffectRadiusScale );

        distance = vector3f::distance( data->worldLocation + vector3f( 0, data->h, 0 ), position );
        data->densities[2] -= metaball_function( distance, radius * data->particleRadiusToEffectRadiusScale );

        distance = vector3f::distance( data->worldLocation - vector3f( 0, data->h, 0 ), position );
        data->densities[3] -= metaball_function( distance, radius * data->particleRadiusToEffectRadiusScale );

        distance = vector3f::distance( data->worldLocation + vector3f( 0, 0, data->h ), position );
        data->densities[4] -= metaball_function( distance, radius * data->particleRadiusToEffectRadiusScale );

        distance = vector3f::distance( data->worldLocation - vector3f( 0, 0, data->h ), position );
        data->densities[5] -= metaball_function( distance, radius * data->particleRadiusToEffectRadiusScale );

    } else {

        float distance = vector3f::distance( data->worldLocation, position );
        data->densities[0] -= metaball_function( distance, radius * data->particleRadiusToEffectRadiusScale );
    }
}

/**
 *	Calculates the densities for all points in the worldBounds to the worldLocation
 *
 *	@param	worldLocation			Find the density at this point.
 */
float particle_metaball_is_policy::get_density( const frantic::graphics::vector3f& worldLocation ) const {
    const float support = m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale;

    density_user_data data( worldLocation, m_positionAccessor, m_radiusAccessor, m_implicitThreshold,
                            m_particleRadiusToEffectRadiusScale );
    m_particles.process_particles_in_range( &data, worldLocation, support, &particle_metaball_is_policy::get_density );

    return data.densities[0];
}

/**
 *	Calculates the densities at the six points surrounding worldLocation when h is applied in all directions. Then
 *uses those densities to estimate the gradient with standard central-difference formulas.
 *
 *	@param	worldLocation			Find the density at this point.
 *	@param	h						The small step to use.
 */
frantic::graphics::vector3f particle_metaball_is_policy::get_gradient( const frantic::graphics::vector3f& worldLocation,
                                                                       float h ) {
    const float support = m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale + h;

    density_user_data data( worldLocation, m_positionAccessor, m_radiusAccessor, m_implicitThreshold,
                            m_particleRadiusToEffectRadiusScale, h, true );
    m_particles.process_particles_in_range( &data, worldLocation, support, &particle_metaball_is_policy::get_density );

    float x = data.densities[0] - data.densities[1];
    float y = data.densities[2] - data.densities[3];
    float z = data.densities[4] - data.densities[5];

    return vector3f( x, y, z ) / 2 * h;
}

frantic::graphics::vector3f particle_metaball_is_policy::find_edge_isosurface_location(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::vector<char*>& particles ) const {

    particles.clear();

    vector3f worldCorner0 = corner_sample_coord_to_world( voxelCorner0 );
    vector3f worldCorner1 = corner_sample_coord_to_world( voxelCorner1 );
    const int solveAxis = ( worldCorner0 - worldCorner1 ).get_largest_axis();

    vector3f surfaceLocation;

    // Refine the vertex position the requested number of times, to make it more closely approximate the isosurface.
    if( m_vertexRefinement > 0 ) {
        const float support = m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale;

        // Make a bounding box containing the extents of this voxel (this will always be an axis-aligned line)
        // This way, we only need to get the nearby particles once, though as a trade off we will get more particles
        // than we need for any individual evaluation.
        boundbox3f vertRange( worldCorner0 );
        vertRange += worldCorner1;
        m_particles.get_particles_in_range( vertRange, support, particles );

        char** particlesBegin = particles.size() > 0 ? &particles[0] : 0;
        char** particlesEnd = particlesBegin + particles.size();

        metaball_vert_refine_eval solver( worldCorner0, solveAxis, particlesBegin, particlesEnd, m_positionAccessor,
                                          m_radiusAccessor, m_implicitThreshold, m_particleRadiusToEffectRadiusScale );

        int iterCount = 0;
        float refined = find_root( solver, worldCorner0[solveAxis], worldCorner1[solveAxis], density0, density1,
                                   m_vertexRefinementEpsilon, m_vertexRefinement, iterCount );
        // m_vrEvalCount += iterCount;

        surfaceLocation = worldCorner0;
        surfaceLocation[solveAxis] = refined;

        // Filter the particle array to only contain the nearby particles we need
        const float supportSquared = frantic::math::square( support );
        for( unsigned i = 0; i < particles.size(); ) {
            float distanceSquared = vector3f::distance_squared( surfaceLocation, m_positionAccessor( particles[i] ) );
            if( distanceSquared >= supportSquared ) {
                // Remove the particle at the back
                particles[i] = particles.back();
                particles.pop_back();
            } else {
                ++i;
            }
        }
    } else {
        surfaceLocation = worldCorner0;
        float alpha = fabsf( density0 / ( density1 - density0 ) );
        surfaceLocation[solveAxis] = ( 1 - alpha ) * worldCorner0[solveAxis] + alpha * worldCorner1[solveAxis];

        // make sure this array is empty, so the extra channel code knows it needs to generate it.
        // this is done at the start of the function
        // particles.clear();
    }

    return surfaceLocation;
}

/**
 *	This function finds the location of the isosurface between the two voxel corners provided,
 *	and sets the internal state of the policy.
 *	@param	density0		First voxel corner value
 *	@param	density1		Second voxel corner value
 *	@param	voxelCorner0	First voxel corner position
 *	@param	voxelCorner0	Second voxel corner position
 */
void particle_metaball_is_policy::find_edge_isosurface_location( float density0, float density1,
                                                                 const vector3& voxelCorner0,
                                                                 const vector3& voxelCorner1 ) {
    m_isosurfaceLocationWorld =
        find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1, m_isosurfaceLocationParticles );
}

const vector3f& particle_metaball_is_policy::get_isosurface_location_world() const { return m_isosurfaceLocationWorld; }

////////
// Functions for dealing with the named channels
////////

/**
 *  Set the internal state of the policy so that it points to the
 * specified location.  This is intended to be called before
 * add_vertex_to_channels if you have manually added a vertex at
 * the specified location, particularly as in generate_disambiguated_faces_for_plane.
 */
void particle_metaball_is_policy::set_isosurface_location_world( const frantic::graphics::vector3f& worldLocation ) {
    m_isosurfaceLocationWorld = worldLocation;
    m_isosurfaceLocationParticles.clear();
}

/**
 *	Given a set of writable trimesh3 accessors, this should add one data element to each of the channels,
 *	which corresponds the the channel names given in prepare_channel_accessors.  The function uses the internal
 *state of the mc policy for the location at which to look up the data
 */
void particle_metaball_is_policy::add_vertex_to_channels(
    std::vector<geometry::trimesh3_vertex_channel_general_accessor>& outputChannels ) {
    // If necessary, get the particles that are near the sample location
    if( m_isosurfaceLocationParticles.empty() ) {
        m_particles.get_particles_in_range( m_isosurfaceLocationWorld,
                                            m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale,
                                            m_isosurfaceLocationParticles );
    }

    if( outputChannels.size() > 0 ) {
        const std::size_t vertIndex = outputChannels[0].size();
        for( unsigned i = 0; i < outputChannels.size(); ++i ) {
            if( outputChannels[i].size() <= vertIndex ) {
                outputChannels[i].set_vertex_count( vertIndex + 1 );
            }
        }

        detail::metaball_vertex_workspace workspace;
        populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, m_isosurfaceLocationWorld,
                                  m_isosurfaceLocationParticles, workspace, m_channelMapWeightedSum );
    }
}

/**
 *	This function finds the location of the isosurface between the two voxel corners provided,
 *	and sets the internal state of the policy.
 *	@param	density0		First voxel corner value
 *	@param	density1		Second voxel corner value
 *	@param	voxelCorner0	First voxel corner position
 *	@param	voxelCorner0	Second voxel corner position
 *	@param	vertIndex		Index of the vert in the accompanying mesh
 *	@param	outputChannels	Trimesh3 channel accessors for any channel data also to be added
 *	@param	outMesh			The mesh to add the vert to.
 *	@param	meshMutex		A mutex for write locking if the function is used in a multithreaded context.
 */
void particle_metaball_is_policy::add_edge_isosurface_vertex_to_mesh(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
    frantic::geometry::trimesh3& outMesh, detail::metaball_vertex_workspace& workspace ) const {

    typedef detail::metaball_vertex_workspace::particles_t particles_t;
    particles_t& particles = workspace.particles;

    vector3f surfaceLocation =
        find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1, particles );

    const float kernelCompactSupport = m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale;

    {
        // add a vert to the mesh at the surface location
        // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
        outMesh.vertices_ref()[vertIndex] = surfaceLocation;
    }

    if( outputChannels.empty() )
        return;

    // If we didnt have to fetch particles for the surface solve (ie, 0 vertex refinement)
    // we need to fetch them now.
    if( particles.empty() ) {
        m_particles.get_particles_in_range( surfaceLocation, kernelCompactSupport, particles );
    }

    // In some cases with no vertex refinement, there may be no particles
    // in range of the selected vertex position.
    // Expand the search radius to pick up such distant particles.
    if( particles.empty() ) {
        m_particles.get_particles_in_range( surfaceLocation, kernelCompactSupport + m_meshingVCS.voxel_length(),
                                            particles );
    }

    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, surfaceLocation, particles, workspace,
                              m_channelMapWeightedSum );
}

/**
 *	This function populates all vertex channels specified in the channel propagation
 *	policy using info from the surface policy.
 *	<p>
 *	All other channel info is preserved.
 *
 *	@param	outMesh		Mesh containing vertex and channel info.
 *	@param	cpp			A channel propagation policy.
 */
void particle_metaball_is_policy::populate_mesh_channels(
    frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_propagation_policy& cpp ) const {
    std::vector<frantic::tstring> previousNames;                // all channels names in the mesh
    std::vector<frantic::tstring> names;                        // all channels that are going to be populated
    std::vector<std::pair<std::size_t, data_type_t>> dataTypes; // types and arity of each channel to be populated
    std::vector<frantic::channels::channel_general_accessor> inputChannels; // channel accessors to channel_map
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>
        outputChannels; // channel accessors to vertex_channel

    // tbb::task_scheduler_init taskSchedulerInit; //please move out.

    // If channels already exist throw error.
    outMesh.get_vertex_channel_names( previousNames );
    for( unsigned i = 0; i < previousNames.size(); ++i )
        if( cpp.is_channel_included( previousNames[i] ) )
            throw runtime_error( "particle_metaball_is_policy.populate_mesh_channels: Attempt to populate a channel "
                                 "that already exists." );

    // Find the names of all channels found in both cpp and the particles.
    get_channel_names( cpp, names );

    // Prepare the input channel accessors and get the data types which
    // will be used to create new channels in the mesh.
    dataTypes.reserve( names.size() );
    inputChannels.reserve( names.size() );
    for( unsigned i = 0; i < names.size(); ++i ) {
        inputChannels.push_back( m_particles.get_channel_map().get_general_accessor( names[i] ) );
        dataTypes.push_back( std::make_pair( inputChannels.back().arity(), inputChannels.back().data_type() ) );
    }

    // Create the new vertex and get accessors to those vertex channels.
    for( size_t i = 0; i < names.size(); ++i ) {
        outMesh.add_vertex_channel_raw( names[i], dataTypes[i].first, dataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( names[i] ) );
    }

    channel_map_weighted_sum channelMapWeightedSum( m_particles.get_channel_map(), cpp );

    // Populate each vertex channel.
    tbb::parallel_for(
        tbb::blocked_range<size_t>( 0, outMesh.vertex_count() ),
        metaball_populate_vertex_channels( outputChannels, outMesh, inputChannels, this, channelMapWeightedSum ) );
}

/**
 *	This function populates all vertex channels using data from the implicit
 *	surface policy for the input vertex.
 *
 *	@param	inputChannels	Accessors to particle channels used for population.
 *	@param	outputChannels	Accessors to each channel that should be populated.
 *	@param	vert			Location of the vertex.
 *	@param	vertIndex		Index of the vertex.
 */
void particle_metaball_is_policy::populate_vertex_channels(
    std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert ) const {
    if( inputChannels.size() != outputChannels.size() )
        throw std::runtime_error(
            "particle_metaball_is_policy.populate_vertex_channels: inputChannels and outputChannels "
            "are not the same size" );

    // Get the particles that are near the sample location
    std::vector<char*> particlesInRange;
    m_particles.get_particles_in_range( vert, m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale,
                                        particlesInRange );

    vertex_workspace_t vertexWorkspace;
    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, particlesInRange, vertexWorkspace,
                              m_channelMapWeightedSum );
}

void particle_metaball_is_policy::populate_vertex_channels(
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace ) const {
    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, vert, workspace, m_channelMapWeightedSum );
}

void particle_metaball_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace,
    const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const {
    workspace.particles.clear();
    m_particles.get_particles_in_range( vert, m_maximumParticleRadius * m_particleRadiusToEffectRadiusScale,
                                        workspace.particles );
    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, workspace.particles, workspace,
                              channelMapWeightedSum );
}

void particle_metaball_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, const std::vector<char*>& particlesInRange, vertex_workspace_t& workspace,
    const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const {
    typedef vertex_workspace_t::particles_t particles_t;
    vector<float>& sampleWeights = workspace.sampleWeights;

    get_sample_weights( vert, particlesInRange, sampleWeights );

    populate_vertex_channels_impl( inputChannels, outputChannels, vertIndex, sampleWeights, particlesInRange, workspace,
                                   channelMapWeightedSum );
}

/*******
 * Functions for the zhu bridson policy
 *******/

namespace {
/**
 * This is the function used to do interpolation in the Zhu/Bridson implicit surface.
 *
 * @param  distanceSquared               The square of the distance from the current sample point to the particle.
 * @param  kernelCompactSupportSquared   The square of the radius of the compact support of the function.
 */
inline float zhu_bridson_kernel( float distanceSquared, float kernelCompactSupportSquared ) {
    if( distanceSquared < kernelCompactSupportSquared ) {
        float x = 1 - distanceSquared / kernelCompactSupportSquared;
        return x * x * x;
    } else {
        return 0;
    }
}

inline float zhu_bridson_kernel( float distanceSquared, float kernelCompactSupportSquared,
                                 float invKernelCompactSupportSquared ) {
    if( distanceSquared < kernelCompactSupportSquared ) {
        float x = 1 - distanceSquared * invKernelCompactSupportSquared;
        return x * x * x;
    } else {
        return 0;
    }
}

inline float_v zhu_bridson_kernel( const float_v& distanceSquared, float kernelCompactSupportSquared,
                                   float invKernelCompactSupportSquared ) {
    const float_v mask = distanceSquared < float_v( kernelCompactSupportSquared );
    const float_v x = float_v( 1 ) - distanceSquared * float_v( invKernelCompactSupportSquared );
    return mask & boost::math::pow<3>( x );
}

/**
 * Compute the Zhu/Bridson kernel function at the specified distanceSquared.
 * The "inv_support" suffix is there to indicate that the second parameter is inverted.
 *
 * @param  distanceSquared                  The square of the distance from the current sample point to the particle.
 * @param  invKernelCompactSupportSquared   The inverse square of the radius of the compact support of the function.
 */
inline float_v zhu_bridson_kernel_inv_support( float_v distanceSquared, float_v invKernelCompactSupportSquared ) {
    float_v x = float_v( 1 ) - distanceSquared * invKernelCompactSupportSquared;
    float_v c = boost::math::pow<3>( x );
#ifdef FRANTIC_HAS_SSE2
    return std::max( c, _mm_setzero_ps() );
#else
    return std::max( c, 0 );
#endif
}

} // anonymous namespace

namespace detail {

zhu_bridson_sparse_data::zhu_bridson_sparse_data()
    : useIsInitializedVoxel( false )
    , useSimd( false ) {}

void zhu_bridson_sparse_data::reset_for_dense_plane_evaluation(
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
    const frantic::channels::channel_accessor<float>& radiusAccessor, const float kernelCompactSupport_,
    const frantic::graphics::boundbox3& voxelExtents_, const frantic::volumetrics::voxel_coord_system& meshingVCS_ ) {
    assert( voxelExtents_.zsize() == 1 );

    reset_for_dense_block_evaluation( positionAccessor, radiusAccessor, kernelCompactSupport_, voxelExtents_,
                                      meshingVCS_ );

    // set up mutexes
    this->useMutexes = true;
    this->mutexes.reset( new tbb::spin_mutex[voxelExtents_.ysize()] );
}

void zhu_bridson_sparse_data::reset_for_dense_block_evaluation(
    const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
    const frantic::channels::channel_accessor<float>& radiusAccessor, const float kernelCompactSupport_,
    const frantic::graphics::boundbox3& voxelExtents_, const frantic::volumetrics::voxel_coord_system& meshingVCS_ ) {
    this->position = positionAccessor;
    this->radius = radiusAccessor;
    this->kernelCompactSupport = kernelCompactSupport_;
    this->kernelCompactSupportSquared = frantic::math::square( kernelCompactSupport_ );
    this->voxelExtents = voxelExtents_;
    this->meshingVCS = meshingVCS_;

    this->useSimd = float_v::static_size > 1 && voxelExtents_.xsize() % float_v::static_size == 0;

    std::size_t gridParticleCount =
        this->useSimd ? voxelExtents_.get_volume() / float_v::static_size : voxelExtents_.get_volume();

    std::size_t gridParticleSize =
        this->useSimd ? sizeof( zhu_bridson_grid_particle_simd ) : sizeof( zhu_bridson_grid_particle );

    // Set up the sparse data structure with all the info we need to determine which grid points a particle
    // will interact with
    char* pGridEvaluation = 0;
    if( this->useSimd ) {
        this->gridEvaluationSimd.resize( gridParticleCount );
        pGridEvaluation = reinterpret_cast<char*>( &this->gridEvaluationSimd[0] );
    } else {
        this->gridEvaluation.resize( gridParticleCount );
        pGridEvaluation = reinterpret_cast<char*>( &this->gridEvaluation[0] );
    }
    assert( pGridEvaluation );

    this->gridParticleChannel = 0;
    this->channelAccessors.clear();
    this->channelAccessors.push_back(
        frantic::channels::channel_general_accessor( 0, gridParticleSize, frantic::channels::data_type_int8 ) );
    this->channelData.clear();
    this->channelData.push_back( reinterpret_cast<char*>( pGridEvaluation ) );
    this->collectRunData = false;
    this->gridParticleCount = 0;

    this->useIsInitializedVoxel = true;
    this->isInitializedVoxel.resize( gridParticleCount );
    memset( &this->isInitializedVoxel[0], 0, gridParticleCount );

    // not threaded
    this->mutexes.reset();
    this->useMutexes = false;
}

} // namespace detail

namespace {

/**
 *	Data structures used internally
 */
struct zhu_bridson_data {
    channel_accessor<vector3f> position;
    channel_accessor<float> radius;
    float kernelCompactSupportSquared;
    std::vector<detail::zhu_bridson_grid_particle> gridEvaluation;
};

struct zhu_bridson_density_user_data {
    frantic::graphics::vector3f worldLocation;
    frantic::channels::channel_accessor<float> radiusAccessor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor;
    float h;
    float kernelCompactSupportSquared;
    float blendedRadius[6];
    float totalWeight[6];
    frantic::graphics::vector3f blendedPosition[6];
    bool getGradient;

    zhu_bridson_density_user_data(
        const frantic::graphics::vector3f& worldLocation,
        const frantic::channels::channel_accessor<float>& radiusAccessor,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        float kernelCompactSupportSquared, float h = 0, bool getGradient = false )
        : worldLocation( worldLocation )
        , radiusAccessor( radiusAccessor )
        , positionAccessor( positionAccessor )
        , kernelCompactSupportSquared( kernelCompactSupportSquared )
        , h( h )
        , getGradient( getGradient ) {
        for( int i = 0; i < 6; i++ ) {
            blendedRadius[i] = 0;
            totalWeight[i] = 0;
            blendedPosition[i] = vector3f( 0 );
        }
    }
};

/**
 *	Utility function.  This gets run on every particle and grid point pair that interacts.
 *	The interaction data is saved in the grid particle.  The function gets called by the
 *	particle_grid_interactions function of the particle grid tree.
 *
 *	@param	userData	void pointer to a zhu_bridson_sparse_data struct
 *	@param	treeParticle	pointer to the data in a tree particle that this function will interact with
 *	@param	gridParticlePosition	position of the particle on the grid to interact with
 *	@param	gridParticleRaw		the grid particle data to be written to
 */
void zhu_bridson_contribution_function( void* userData, char* treeParticle, const vector3f& gridParticlePosition,
                                        char* gridParticleRaw ) {
    const zhu_bridson_data* data = reinterpret_cast<zhu_bridson_data*>( userData );
    detail::zhu_bridson_grid_particle* gridParticle =
        reinterpret_cast<detail::zhu_bridson_grid_particle*>( gridParticleRaw );

    vector3f treeParticleOffset = data->position( treeParticle ) - gridParticlePosition;

    float distanceSquared = treeParticleOffset.get_magnitude_squared();

    if( distanceSquared > data->kernelCompactSupportSquared )
        return;

    float weight = zhu_bridson_kernel( distanceSquared, data->kernelCompactSupportSquared );

    gridParticle->weight += weight;
    gridParticle->blendedOffset += weight * treeParticleOffset;
    gridParticle->blendedRadius += weight * data->radius( treeParticle );
}

class zhu_bridson_sparse_contribution_evaluator
    : public sparse_contribution_evaluator<zhu_bridson_sparse_contribution_evaluator> {
  public:
    typedef detail::zhu_bridson_sparse_data sparse_data_t;

    typedef detail::zhu_bridson_grid_particle grid_particle_t;

    zhu_bridson_sparse_contribution_evaluator( sparse_data_t& data )
        : m_particleEffectRadiusSquared( data.kernelCompactSupportSquared )
        , m_invParticleEffectRadiusSquared( 1 / data.kernelCompactSupportSquared )
        , m_gridParticleChannel( data.gridParticleChannel )
        , m_gridParticles( (grid_particle_t*)( data.channelData[data.gridParticleChannel] ) ) {}

    BOOST_FORCEINLINE void evaluate_x_run( const char* particle, const vector3f& particlePosition, int xMin, int xMax,
                                           float particleOffsetY, float particleOffsetZ, float distanceSquaredYZ,
                                           int extentNum, int dataOffset, sparse_data_t& data ) {
        if( data.collectRunData ) {
            // Push the run into the appropriate run_tree.  This will initialize the data in if it isn't already.
            std::pair<int, int> run( xMin - data.voxelExtents.minimum().x, xMax - data.voxelExtents.minimum().x - 1 );
            if( run.second >= run.first )
                data.runTrees[extentNum].insert_run( run, data.extentData[extentNum], data.initChannelData );

            for( int x = xMin; x <= xMax; x++ ) {
                int i = x + dataOffset;

                float worldVoxelCenterX = data.meshingVCS.get_world_x_voxel_center( x );
                float particleOffsetX = particlePosition.x - worldVoxelCenterX;
                float distanceSquaredX = frantic::math::square( particleOffsetX );

                vector3f particleOffset( particleOffsetX, particleOffsetY, particleOffsetZ );

                float distanceSquared = distanceSquaredX + distanceSquaredYZ;
                float weight = zhu_bridson_kernel( distanceSquared, m_particleEffectRadiusSquared,
                                                   m_invParticleEffectRadiusSquared );

                m_gridParticles[i].weight += weight;
                m_gridParticles[i].blendedOffset += weight * particleOffset;
                m_gridParticles[i].blendedRadius += weight * data.radius( particle );

                // add in the effect of this particle to this location in each channel
                size_t channelCount = data.channelData.size();
                for( size_t channelNum = 0; channelNum < channelCount; ++channelNum ) {

                    // skip the grid particle channel
                    if( channelNum == (size_t)m_gridParticleChannel )
                        continue;

                    // add in the weighted contribution of this particle for the channel
                    frantic::channels::channel_general_accessor& ca = data.channelAccessors[channelNum];
                    ca.weighted_increment( weight, ca.get_channel_data_pointer( particle ),
                                           data.channelData[channelNum] + i * ca.primitive_size() );
                }
            }

        } else {

            for( int x = xMin; x <= xMax; x++ ) {
                int i = x + dataOffset;

                float worldVoxelCenterX = data.meshingVCS.get_world_x_voxel_center( x );
                float particleOffsetX = particlePosition.x - worldVoxelCenterX;
                float distanceSquaredX = frantic::math::square( particleOffsetX );

                vector3f particleOffset( particleOffsetX, particleOffsetY, particleOffsetZ );

                float distanceSquared = distanceSquaredX + distanceSquaredYZ;
                float weight = zhu_bridson_kernel( distanceSquared, m_particleEffectRadiusSquared,
                                                   m_invParticleEffectRadiusSquared );

                if( data.useIsInitializedVoxel && data.isInitializedVoxel[i] == 0 ) {
                    data.isInitializedVoxel[i] = true;
                    m_gridParticles[i].weight = weight;
                    m_gridParticles[i].blendedOffset = weight * particleOffset;
                    m_gridParticles[i].blendedRadius = weight * data.radius( particle );
                } else {
                    m_gridParticles[i].weight += weight;
                    m_gridParticles[i].blendedOffset += weight * particleOffset;
                    m_gridParticles[i].blendedRadius += weight * data.radius( particle );
                }
            }
        }
    }

  private:
    float m_particleEffectRadiusSquared;
    float m_invParticleEffectRadiusSquared;
    int m_gridParticleChannel;
    grid_particle_t* m_gridParticles;
};

class zhu_bridson_sparse_contribution_simd_evaluator
    : public sparse_contribution_evaluator<zhu_bridson_sparse_contribution_simd_evaluator> {
  public:
    typedef detail::zhu_bridson_sparse_data sparse_data_t;

    zhu_bridson_sparse_contribution_simd_evaluator( sparse_data_t& data )
        : m_particleEffectRadiusSquared( data.kernelCompactSupportSquared )
        , m_invParticleEffectRadiusSquared( 1 / data.kernelCompactSupportSquared )
        , m_gridParticleChannel( data.gridParticleChannel )
        , m_gridParticles( &data.gridEvaluationSimd[0] ) {}

    BOOST_FORCEINLINE void evaluate_x_run( const char* particle, const vector3f& particlePosition, int xMin, int xMax,
                                           float particleOffsetY, float particleOffsetZ, float distanceSquaredYZ_,
                                           int /*extentNum*/, int dataOffset, sparse_data_t& data ) {
        float_v distanceSquaredYZ( distanceSquaredYZ_ );

        float radius = data.radius( particle );

        const int w = detail::zhu_bridson_grid_particle_simd::width;
        const int extentsMinimumX = data.voxelExtents.xminimum();
        // round down to nearest multiple of w
        int xBegin = ( xMin - extentsMinimumX ) / w * w + extentsMinimumX;
        // round up to next multiple of w
        int xEnd = ( xMax - extentsMinimumX + w ) / w * w + extentsMinimumX;
        int i = ( xBegin + dataOffset ) / w;
        for( int xStepped = xBegin; xStepped < xEnd; xStepped += w, ++i ) {
            int_v xInt = int_v( xStepped ) + create_int_v_run();

            float_v x( xInt );

            int_v mask( xInt >= xMin );
            mask.and_not( xInt > xMax );

            detail::zhu_bridson_grid_particle_simd& gridParticle = m_gridParticles[i];

            float_v gridPointX = ( x + 0.5 ) * data.meshingVCS.voxel_length() + data.meshingVCS.world_origin().x;
            float_v particleOffsetX = float_v( particlePosition.x ) - gridPointX;
            float_v distanceSquaredX = frantic::math::square( particleOffsetX );
            float_v distanceSquared = float_v( distanceSquaredX ) + distanceSquaredYZ;

            float_v weight =
                zhu_bridson_kernel( distanceSquared, m_particleEffectRadiusSquared, m_invParticleEffectRadiusSquared );
            weight &= float_v::reinterpret( mask );

            if( data.useIsInitializedVoxel && data.isInitializedVoxel[i] == 0 ) {
                gridParticle.weight = weight;
                gridParticle.blendedRadius = weight * radius;
                gridParticle.blendedOffset =
                    weight * vector3t<float_v>( particleOffsetX, particleOffsetY, particleOffsetZ );
                data.isInitializedVoxel[i] = true;
            } else {
                gridParticle.weight += weight;
                gridParticle.blendedRadius += weight * radius;
                gridParticle.blendedOffset +=
                    weight * vector3t<float_v>( particleOffsetX, particleOffsetY, particleOffsetZ );
            }
        }
    }

  private:
    float m_particleEffectRadiusSquared;
    float m_invParticleEffectRadiusSquared;
    int m_gridParticleChannel;
    detail::zhu_bridson_grid_particle_simd* m_gridParticles;
};

/**
 *	Utility function.  This gets run on every particle that will interact with a given
 *	plane of data in the fill_sparse_channel_data function, to determine which grid points
 *	on the plane that the particle interacts with.	The arrays of data that get passed in the user
 *	data struct will get populated accordingly, and the run indexing data for the arrays will be
 *	saved in the vector of runtrees also provided in the user data struct.
 *
 *	@param	userData	void pointer to a zhu_bridson_sparse_data struct
 *	@param	treeParticle	pointer to the data in a tree particle that this function will interact with
 */
void zhu_bridson_sparse_contribution_function( void* userData, char* treeParticle ) {
    detail::zhu_bridson_sparse_data* data = reinterpret_cast<detail::zhu_bridson_sparse_data*>( userData );

    vector3f treeParticlePosition = data->position( treeParticle );

    bool hasParticle = false;

    if( data->useSimd ) {
        zhu_bridson_sparse_contribution_simd_evaluator f( *data );
        f.evaluate( treeParticle, treeParticlePosition, data->kernelCompactSupport, *data );
        hasParticle = f.has_particle_on_grid();
    } else {
        zhu_bridson_sparse_contribution_evaluator f( *data );
        f.evaluate( treeParticle, treeParticlePosition, data->kernelCompactSupport, *data );
        hasParticle = f.has_particle_on_grid();
    }

    if( hasParticle ) {
        ++( data->gridParticleCount );
    }
}

// necessary for multithreading with tbb
class ZhuBridsonBody {

    detail::zhu_bridson_sparse_data* const m_data;
    std::vector<char*>* const m_particles;

  public:
    void operator()( const tbb::blocked_range<size_t>& r ) const {
        // std::cout << "operator()" << std::endl;
        detail::zhu_bridson_sparse_data* data = m_data;
        std::vector<char*>* particles = m_particles;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            // std::cout << "processing particle " << i << std::endl;
            zhu_bridson_sparse_contribution_function( data, particles->at( i ) );
        }
    }

    ZhuBridsonBody( detail::zhu_bridson_sparse_data* data, vector<char*>* particles )
        : m_data( data )
        , m_particles( particles ) {}

    ZhuBridsonBody& operator=( const ZhuBridsonBody& ) { return *this; } // unimplemented
};

class zhu_bridson_vert_refine_eval {
    frantic::graphics::vector3f m_voxelPosition0;

    int m_solveAxis;

    char** m_particlesBegin;
    char** m_particlesEnd;

    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_accessor<float> m_radiusAccessor;

    float m_kernelCompactSupportSquared;
    float m_invKernelCompactSupportSquared;
    const particle_zhu_bridson_is_policy& m_zhuBridsonPolicy;

    zhu_bridson_vert_refine_eval& operator=( const zhu_bridson_vert_refine_eval& ); // not implemented

  public:
    zhu_bridson_vert_refine_eval(
        const particle_zhu_bridson_is_policy& zhuBridsonPolicy, const frantic::graphics::vector3f& voxelCoord0,
        const int solveAxis, char** particlesBegin, char** particlesEnd,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::channels::channel_accessor<float>& radiusAccessor, float kernelCompactSupport )
        : m_zhuBridsonPolicy( zhuBridsonPolicy )
        , m_voxelPosition0( voxelCoord0 )
        , m_solveAxis( solveAxis )
        , m_particlesBegin( particlesBegin )
        , m_particlesEnd( particlesEnd )
        , m_positionAccessor( positionAccessor )
        , m_radiusAccessor( radiusAccessor ) {
        m_kernelCompactSupportSquared = boost::math::pow<2>( kernelCompactSupport );
        m_invKernelCompactSupportSquared = 1.f / m_kernelCompactSupportSquared;
    }

    float operator()( float x ) const {
        vector3f vertTest( m_voxelPosition0 );
        vertTest[m_solveAxis] = x;

        // evaluate the density at this location
        frantic::graphics::vector3f blendedOffset( 0 );
        float blendedRadius = 0, totalWeight = 0;
        for( char** i = m_particlesBegin; i != m_particlesEnd; ++i ) {
            const vector3f particleOffset = m_positionAccessor( *i ) - vertTest;
            const float radius = m_radiusAccessor( *i );
            float distanceSquared = particleOffset.get_magnitude_squared();
            float weight =
                zhu_bridson_kernel( distanceSquared, m_kernelCompactSupportSquared, m_invKernelCompactSupportSquared );
            totalWeight += weight;
            blendedOffset += weight * particleOffset;
            blendedRadius += weight * radius;
        }

        float distanceFromSurface;
        if( totalWeight > 0 ) {
            blendedOffset *= ( 1.f / totalWeight );
            blendedRadius *= ( 1.f / totalWeight );
            distanceFromSurface = blendedOffset.get_magnitude() - blendedRadius;
        } else {
            totalWeight = 0;
            distanceFromSurface = m_zhuBridsonPolicy.get_default_outside_distance();
        }

        return distanceFromSurface + m_zhuBridsonPolicy.isosurface_function_compensation( totalWeight );
    }
};

class zhu_bridson_vert_refine_eval_simd {
    float_v m_voxelPositionA;
    float_v m_voxelPositionB;

    float_v m_invKernelCompactSupportSquared;

    float m_defaultOutsideDistance;

    int m_solveAxis;

    frantic::volumetrics::implicitsurface::detail::xyzr_packet_array& m_particles;

    const particle_zhu_bridson_is_policy& m_isp;

    zhu_bridson_vert_refine_eval_simd& operator=( const zhu_bridson_vert_refine_eval_simd& ); // not implemented

  public:
    typedef frantic::volumetrics::implicitsurface::detail::xyzr_packet_array particles_t;

    zhu_bridson_vert_refine_eval_simd( const particle_zhu_bridson_is_policy& isp,
                                       const frantic::graphics::vector3f& voxelCoord0, const int solveAxis,
                                       particles_t& particles, float kernelCompactSupport )
        : m_solveAxis( solveAxis )
        , m_particles( particles )
        , m_isp( isp ) {
        using frantic::math::square;

        if( solveAxis == 0 ) {
            m_voxelPositionA = voxelCoord0.y;
            m_voxelPositionB = voxelCoord0.z;
        } else if( solveAxis == 1 ) {
            m_voxelPositionA = voxelCoord0.x;
            m_voxelPositionB = voxelCoord0.z;
        } else {
            m_voxelPositionA = voxelCoord0.x;
            m_voxelPositionB = voxelCoord0.y;
        }

        m_invKernelCompactSupportSquared = 1.f / square( kernelCompactSupport );
    }

    float operator()( float xScalar ) const {
        using frantic::math::square;

        vector3t<float_v> vertTest;
        if( m_solveAxis == 0 ) {
            vertTest.x = xScalar;
            vertTest.y = m_voxelPositionA;
            vertTest.z = m_voxelPositionB;
        } else if( m_solveAxis == 1 ) {
            vertTest.x = m_voxelPositionA;
            vertTest.y = xScalar;
            vertTest.z = m_voxelPositionB;
        } else {
            vertTest.x = m_voxelPositionA;
            vertTest.y = m_voxelPositionB;
            vertTest.z = xScalar;
        }

        // evaluate the density at this location
        vector3t<float_v> blendedOffset( float_v( 0 ) );
        float_v blendedRadius = 0, totalWeight = 0;
        std::size_t inIndex = 0;
        for( std::size_t endIndex = m_particles.get_filled_packet_count(); inIndex < endIndex; ++inIndex ) {
            vector3t<float_v> position = m_particles.get_particles()[inIndex].position;
            float_v radius = m_particles.get_particles()[inIndex].radius;

            vector3t<float_v> particleOffset = position - vertTest;
            float_v distanceSquared = particleOffset.get_magnitude_squared();
            float_v weight = zhu_bridson_kernel_inv_support( distanceSquared, m_invKernelCompactSupportSquared );
            totalWeight += weight;
            blendedOffset += weight * particleOffset;
            blendedRadius += weight * radius;
        }
        if( m_particles.has_remainder() ) {
            vector3t<float_v> position = m_particles.get_particles()[inIndex].position;
            float_v radius = m_particles.get_particles()[inIndex].radius;

            vector3t<float_v> particleOffset = position - vertTest;
            float_v distanceSquared = vector3t<float_v>::distance_squared( position, vertTest );
            float_v weight = zhu_bridson_kernel_inv_support( distanceSquared, m_invKernelCompactSupportSquared );
            weight &= m_particles.get_remainder_mask();
            totalWeight += weight;
            blendedOffset += weight * particleOffset;
            blendedRadius += weight * radius;
        }

        const float totalWeightSum = totalWeight.sum();
        float distanceFromSurface = 0;
        if( totalWeightSum > 0 ) {
#ifdef FRANTIC_HAS_SSE2
            float invTotalWeight = 1 / totalWeightSum;
            _MM_TRANSPOSE4_PS( blendedOffset.x.native(), blendedOffset.y.native(), blendedOffset.z.native(),
                               blendedRadius.native() );
            // Note that these variables were just transposed, so their
            // names are misleading.  The first component of each variable
            // is actually the x, the second component is the y, etc.
            // So the first component of sum is the sum of all x values, etc.
            float_v sum = blendedOffset.x + blendedOffset.y + blendedOffset.z + blendedRadius;
            sum *= invTotalWeight;
            float blendedRadiusScalar = sum[3];
            float_v sumSquared = square( sum );
            float blendedOffsetMagnitude = std::sqrt( sumSquared[0] + sumSquared[1] + sumSquared[2] );
            distanceFromSurface = blendedOffsetMagnitude - blendedRadiusScalar;
#else
            blendedOffset *= ( 1.f / totalWeight );
            blendedRadius *= ( 1.f / totalWeight );
            distanceFromSurface = blendedOffset.get_magnitude().native() - blendedRadius.native();
#endif
        } else {
            totalWeight = 0;
            distanceFromSurface = m_defaultOutsideDistance;
        }

        return distanceFromSurface + m_isp.isosurface_function_compensation( totalWeightSum );
    }
};

class zhu_bridson_populate_vertex_channels {
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels;
    frantic::geometry::trimesh3& outMesh;
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels;
    const particle_zhu_bridson_is_policy* isp;
    const channel_map_weighted_sum& channelMapWeightedSum;

    zhu_bridson_populate_vertex_channels& operator=( const zhu_bridson_populate_vertex_channels& ); // not implemented

  public:
    zhu_bridson_populate_vertex_channels(
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
        frantic::geometry::trimesh3& outMesh,
        const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
        const particle_zhu_bridson_is_policy* isp, const channel_map_weighted_sum& channelMapWeightedSum )
        : outputChannels( outputChannels )
        , outMesh( outMesh )
        , inputChannels( inputChannels )
        , isp( isp )
        , channelMapWeightedSum( channelMapWeightedSum ) {}

    void operator()( const tbb::blocked_range<size_t>& r ) const {
        detail::zhu_bridson_vertex_workspace workspace;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            isp->populate_vertex_channels( inputChannels, outputChannels, i, outMesh.get_vertex( i ), workspace,
                                           channelMapWeightedSum );
        }
    }
};

} // anonymous namespace

/**
 *	This constructor initializes the policy using the following data
 *	@param	particles	A particle grid tree containing the particle system the policy is to operate on.
 *	@param	maxParticleRadius	The compact support for the zhu/bridson kernal function.
 *	@param	kernelCompactSupport	The compact support for the zhu/bridson kernal function.
 *	@param	lowDensityTrimmingStrength	Required paramater for zhu/bridson
 *	@param	lowDensityTrimmingStrength	Required paramater for zhu/bridson
 *	@param	meshingVCS	The coord system that the policy is to work in.
 *	@param	vertRefinement	The number of refinement steps to take when solving for vertex locations.
 */
particle_zhu_bridson_is_policy::particle_zhu_bridson_is_policy( particle_grid_tree& particles, float maxParticleRadius,
                                                                float effectRadius, float lowDensityTrimmingDensity,
                                                                float lowDensityTrimmingStrength,
                                                                const voxel_coord_system& meshingVCS,
                                                                int vertexRefinement )
    : particle_is_policy_base<particle_zhu_bridson_is_policy>( particles )
    , m_maximumParticleRadius( maxParticleRadius )
    , m_kernelCompactSupport( maxParticleRadius * effectRadius )
    , m_lowDensityTrimmingDensity( lowDensityTrimmingDensity )
    , m_lowDensityTrimmingStrength( lowDensityTrimmingStrength )
    , m_meshingVCS( meshingVCS )
    , m_vertexRefinement( vertexRefinement ) {

    m_positionAccessor = m_particles.get_channel_map().get_accessor<vector3f>( _T("Position") );
    m_radiusAccessor = m_particles.get_channel_map().get_accessor<float>( _T("Radius") );

    // Compute the voxel bounds over which to do the implicit surface conversion
    boundbox3f particleWorldBounds = m_particles.compute_particle_bounds();
    particleWorldBounds.expand( maxParticleRadius * effectRadius );
    m_particleVCSBounds = m_meshingVCS.get_voxel_bounds( particleWorldBounds );
    m_particleVCSBounds.expand( 1 );

    m_kernelAtZero = zhu_bridson_kernel( 0, m_kernelCompactSupport * m_kernelCompactSupport );

    m_vertexRefinementEpsilon = 0.001f * m_meshingVCS.voxel_length(); // could also account for max extents?
                                                                      // m_vrEvalCount = 0;
}

particle_zhu_bridson_is_policy::~particle_zhu_bridson_is_policy() {
    // FF_LOG( debug ) << m_vrEvalCount << std::endl;
}

/**
 *	This returns the voxel coordinate system the policy is working in.
 *	@return	voxel_coord_system	The coord system.
 */
const voxel_coord_system& particle_zhu_bridson_is_policy::get_voxel_coord_system() const { return m_meshingVCS; }

/**
 *	Returns the XYZ bounds of the particle system the policy is working on.
 *	@return	boundbox3	The bounds.
 */
boundbox3 particle_zhu_bridson_is_policy::get_voxel_bounds() const { return m_particleVCSBounds; }

/**
 *	This returns the interface widths in voxels for the policy.
 *
 *	@param	interfaceWidthInside	assigned the inside interface width
 *	@param	interfaceWidthOutside	assigned the outside interface width
 */
void particle_zhu_bridson_is_policy::get_interface_widths( float& interfaceWidthInside,
                                                           float& interfaceWidthOutside ) const {
    vector3 dist = m_particleVCSBounds.maximum() - m_particleVCSBounds.minimum();
    vector3f distf( float( dist.x ), float( dist.y ), float( dist.z ) );
    interfaceWidthInside = distf.get_magnitude();
    interfaceWidthOutside = m_kernelCompactSupport - m_maximumParticleRadius;
}

/**
 *	Returns the exterior region code for the particle is policy.  Particle is policies have an exterior
 *	region code of -1;
 */
int particle_zhu_bridson_is_policy::get_exterior_region_code() const { return -1; }

void particle_zhu_bridson_is_policy::get_sample_weights( const vector3f& position, const std::vector<char*>& particles,
                                                         std::vector<float>& outWeights ) const {
    const size_t particleCount = particles.size();

    outWeights.resize( particleCount );

    if( particleCount == 0 ) {
        return;
    }

    // Compute the sample weights for blending between the different values
    float weightSum = 0;
    const float kernelCompactSupport2 = boost::math::pow<2>( m_kernelCompactSupport );
    for( size_t i = 0; i < particleCount; ++i ) {
        const vector3f& particlePosition = m_positionAccessor( particles[i] );
        const float distance = vector3f::distance_squared( position, particlePosition );
        const float weight = zhu_bridson_kernel( distance, kernelCompactSupport2 );
        outWeights[i] = weight;
        weightSum += weight;
    }

    // Normalize the sample weights
    const float invWeightSum = 1.f / weightSum;
    for( size_t i = 0; i < particleCount; ++i ) {
        outWeights[i] *= invWeightSum;
    }
}

void particle_zhu_bridson_is_policy::get_sample_weights_simd(
    const frantic::graphics::vector3f& position_,
    const frantic::volumetrics::implicitsurface::detail::xyzr_packet_array& particles,
    std::vector<float>& outWeights ) const {

    using frantic::math::square;

    const size_t particleCount = particles.get_particle_count();

    if( particleCount == 0 ) {
        outWeights.resize( 0 );
        return;
    }

    const vector3t<float_v> position( position_.x, position_.y, position_.z );

    outWeights.resize( particles.get_particles().size() * float_v::static_size );

    // Compute the sample weights for blending between the different values
    float_v weightSum = 0;
    const float_v invKernelCompactSupport2 = boost::math::pow<2>( 1.f / m_kernelCompactSupport );
    size_t i = 0;
    for( ; i < particles.get_filled_packet_count(); ++i ) {
        const vector3t<float_v> particlePosition = particles.get_particles()[i].position;
        const float_v distanceSquared = vector3t<float_v>::distance_squared( position, particlePosition );
        const float_v weight = zhu_bridson_kernel_inv_support( distanceSquared, invKernelCompactSupport2 );
        weight.store( &outWeights[i * float_v::static_size] );
        weightSum += weight;
    }
    if( particles.has_remainder() ) {
        const vector3t<float_v> particlePosition = particles.get_particles()[i].position;
        const float_v distanceSquared = vector3t<float_v>::distance_squared( position, particlePosition );
        float_v weight = zhu_bridson_kernel_inv_support( distanceSquared, invKernelCompactSupport2 );
        weight &= particles.get_remainder_mask();
        weight.store( &outWeights[i * float_v::static_size] );
        weightSum += weight;
    }

    // Normalize the sample weights
    const float_v invWeightSum = 1 / weightSum.sum();
    for( size_t i = 0, ie = outWeights.size(); i < ie; i += float_v::static_size ) {
        float_v w = float_v::load( &outWeights[i] );
        w *= invWeightSum;
        w.store( &outWeights[i] );
    }

    outWeights.resize( particleCount );
}

void particle_zhu_bridson_is_policy::get_affected_blocks(
    const frantic::graphics::size3f blockSize, std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
    shared_progress_logger_proxy& progressLogger ) {
    get_affected_blocks_impl( m_particles, get_compact_support(), blockSize, outAffectedBlockCoordinates,
                              progressLogger );
}

/**
 *	This function takes xyExtents and a z coord for a plane, and populates it densely with
 *	density information as determined by the policy.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A vector that will be populated with data.
 */
void particle_zhu_bridson_is_policy::fill_voxel_corner_densities( const boundrect2& xyExtents, int z,
                                                                  std::vector<float>& outVoxelCornerValues ) {

    zhu_bridson_data data;
    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.kernelCompactSupportSquared = m_kernelCompactSupport * m_kernelCompactSupport;

    // Initialize the evaluation grid to all zeros
    // TODO: this should not be allocated every call, but we aren't using this function now
    data.gridEvaluation.resize( xyExtents.get_area() );
    for( int i = 0, ie = xyExtents.get_area(); i != ie; ++i ) {
        data.gridEvaluation[i].blendedOffset.set( 0 );
        data.gridEvaluation[i].blendedRadius = 0;
        data.gridEvaluation[i].weight = 0;
    }

    boundbox3 xyzExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y, z,
                          z );
    m_particles.particle_grid_interactions( xyzExtents, m_meshingVCS, sizeof( detail::zhu_bridson_grid_particle ),
                                            &data.gridEvaluation[0], &zhu_bridson_contribution_function, &data,
                                            m_kernelCompactSupport );

    // Evaluate the implicit function
    outVoxelCornerValues.resize( xyExtents.get_area() );
    int index = 0;
    for( int y = xyExtents.minimum().y; y <= xyExtents.maximum().y; ++y ) {
        for( int x = xyExtents.minimum().x; x <= xyExtents.maximum().x; ++x ) {
            float weight = data.gridEvaluation[index].weight;
            if( weight > 0 ) {
                vector3f blendedOffset = data.gridEvaluation[index].blendedOffset / weight;
                float particleRadius = data.gridEvaluation[index].blendedRadius / weight;
                outVoxelCornerValues[index] =
                    blendedOffset.get_magnitude() - particleRadius + isosurface_function_compensation( weight );
            } else {
                outVoxelCornerValues[index] = m_kernelCompactSupport;
            }
            ++index;
        }
    }
}

struct GridEvaluationToDensityBody {
    const particle_zhu_bridson_is_policy& m_isp;
    const frantic::graphics2d::boundrect2& m_xyExtents;
    const std::vector<boost::uint8_t>& m_isInitializedVoxel;
    const std::vector<detail::zhu_bridson_grid_particle>& m_gridEvaluation;
    float* m_outVoxelCornerValues;
    int m_z;

#pragma warning( push )
#pragma warning( disable : 4822 ) // local class member function does not have a body
    GridEvaluationToDensityBody& operator=( const GridEvaluationToDensityBody& ); // not implemented
#pragma warning( pop )

    void operator()( const tbb::blocked_range<int>& r ) const {
        int index = ( r.begin() - m_xyExtents.minimum().y ) * m_xyExtents.xsize();
        const boost::uint8_t* isInitializedVoxel = m_isInitializedVoxel.size() > 0 ? &m_isInitializedVoxel[0] : 0;
        const detail::zhu_bridson_grid_particle* gridEvaluation =
            m_gridEvaluation.size() > 0 ? &m_gridEvaluation[0] : 0;
        for( int y = r.begin(); y != r.end(); ++y ) {
            for( int x = m_xyExtents.minimum().x; x <= m_xyExtents.maximum().x; ++x ) {
                if( isInitializedVoxel[index] ) {
                    const float weight = gridEvaluation[index].weight;
                    if( weight > 0 ) {
                        const vector3f blendedOffset = gridEvaluation[index].blendedOffset / weight;
                        const float particleRadius = gridEvaluation[index].blendedRadius / weight;
                        m_outVoxelCornerValues[index] = blendedOffset.get_magnitude() - particleRadius +
                                                        m_isp.isosurface_function_compensation( weight );
                    } else {
                        m_outVoxelCornerValues[index] = m_isp.m_kernelCompactSupport;
                    }
                } else {
                    m_outVoxelCornerValues[index] = m_isp.m_kernelCompactSupport;
                }

                ++index;
            }
        }
    }

    GridEvaluationToDensityBody( const particle_zhu_bridson_is_policy& isp,
                                 const frantic::graphics2d::boundrect2& xyExtents, const int z,
                                 const std::vector<boost::uint8_t>& isInitializedVoxel,
                                 const std::vector<detail::zhu_bridson_grid_particle>& gridEvaluation,
                                 float* outVoxelCornerValues )
        : m_isp( isp )
        , m_xyExtents( xyExtents )
        , m_z( z )
        , m_isInitializedVoxel( isInitializedVoxel )
        , m_gridEvaluation( gridEvaluation )
        , m_outVoxelCornerValues( outVoxelCornerValues ) {}
};

void particle_zhu_bridson_is_policy::get_density_from_grid_evaluation_mt(
    const frantic::graphics2d::boundrect2& xyExtents, const int z,
    const std::vector<boost::uint8_t>& isInitializedVoxel,
    const std::vector<detail::zhu_bridson_grid_particle>& gridEvaluation, float* outVoxelCornerValues ) const {
    tbb::blocked_range<int> range( xyExtents.minimum().y, xyExtents.maximum().y + 1 );
    GridEvaluationToDensityBody body( *this, xyExtents, z, isInitializedVoxel, gridEvaluation, outVoxelCornerValues );

#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif
}

void particle_zhu_bridson_is_policy::get_density_from_grid_evaluation(
    const frantic::graphics::boundbox3& xyzExtents, const std::vector<boost::uint8_t>& isInitializedVoxel_,
    const std::vector<detail::zhu_bridson_grid_particle>& gridEvaluation_, float* outVoxelCornerValues ) const {

    boost::int32_t zMin = xyzExtents.zminimum(), zMax = xyzExtents.zmaximum();
    boost::int32_t yMin = xyzExtents.yminimum(), yMax = xyzExtents.ymaximum();
    boost::int32_t xMin = xyzExtents.xminimum(), xMax = xyzExtents.xmaximum();
    const boost::uint8_t* isInitializedVoxel = isInitializedVoxel_.size() > 0 ? &isInitializedVoxel_[0] : 0;
    const detail::zhu_bridson_grid_particle* gridEvaluation = gridEvaluation_.size() > 0 ? &gridEvaluation_[0] : 0;
    int index = 0;
    for( int z = zMin; z <= zMax; ++z ) {
        for( int y = yMin; y <= yMax; ++y ) {
            for( int x = xMin; x <= xMax; ++x ) {
                if( isInitializedVoxel[index] ) {
                    const float weight = gridEvaluation[index].weight;
                    if( weight > 0 ) {
                        const vector3f blendedOffset = gridEvaluation[index].blendedOffset / weight;
                        const float particleRadius = gridEvaluation[index].blendedRadius / weight;
                        outVoxelCornerValues[index] =
                            blendedOffset.get_magnitude() - particleRadius + isosurface_function_compensation( weight );
                    } else {
                        outVoxelCornerValues[index] = m_kernelCompactSupport;
                    }
                } else {
                    outVoxelCornerValues[index] = m_kernelCompactSupport;
                }
                ++index;
            }
        }
    }
}

namespace {

void get_density_from_grid_evaluation_simd_impl(
    const particle_zhu_bridson_is_policy& isp, const frantic::graphics::boundbox3& xyzExtents,
    const std::vector<boost::uint8_t>& isInitializedVoxel_,
    const std::vector<detail::zhu_bridson_grid_particle_simd>& gridEvaluation_, float* outVoxelCornerValues,
    std::size_t indexOffset = 0 ) {
    boost::int32_t zMin = xyzExtents.zminimum(), zMax = xyzExtents.zmaximum();
    boost::int32_t yMin = xyzExtents.yminimum(), yMax = xyzExtents.ymaximum();
    boost::int32_t xMin = xyzExtents.xminimum(), xMax = xyzExtents.xmaximum();

    const boost::uint8_t* isInitializedVoxel = isInitializedVoxel_.size() > 0 ? &isInitializedVoxel_[0] : 0;
    const detail::zhu_bridson_grid_particle_simd* gridEvaluation = gridEvaluation_.size() > 0 ? &gridEvaluation_[0] : 0;

    const int xStep = detail::zhu_bridson_grid_particle_simd::width;

    const float_v kernelCompactSupport( isp.get_compact_support() );

    std::size_t index = indexOffset / float_v::static_size;
    float* out = outVoxelCornerValues + indexOffset;
    for( int z = zMin; z <= zMax; ++z ) {
        for( int y = yMin; y <= yMax; ++y ) {
            for( int x = xMin; x <= xMax; x += xStep ) {
                if( isInitializedVoxel[index] ) {
                    const float_v weight = gridEvaluation[index].weight;
                    const float_v mask = weight > 0;

                    vector3t<float_v> offset = gridEvaluation[index].blendedOffset / weight;

                    float_v distance = offset.get_magnitude();

                    float_v radius = gridEvaluation[index].blendedRadius / weight;

                    float_v density = distance - radius + isp.isosurface_function_compensation( weight );

                    float_v result = float_v::select( kernelCompactSupport, density, mask );

                    result.store( out );
                } else {
                    kernelCompactSupport.store( out );
                }
                ++index;
                out += xStep;
            }
        }
    }
}

struct get_density_from_grid_evaluation_simd_mt_body {
    const particle_zhu_bridson_is_policy& m_isp;
    const frantic::graphics2d::boundrect2& m_xyExtents;
    const std::vector<boost::uint8_t>& m_isInitializedVoxel;
    const std::vector<detail::zhu_bridson_grid_particle_simd>& m_gridEvaluation;
    float* m_outVoxelCornerValues;
    int m_z;

    get_density_from_grid_evaluation_simd_mt_body&
    operator=( const get_density_from_grid_evaluation_simd_mt_body& ); // not implemented

    void operator()( const tbb::blocked_range<int>& r ) const {
        boundbox3 xyzExtents( vector3( m_xyExtents.minimum().x, r.begin(), m_z ),
                              vector3( m_xyExtents.maximum().x, r.end() - 1, m_z ) );
        const std::size_t index = ( r.begin() - m_xyExtents.minimum().y ) * m_xyExtents.xsize();
        get_density_from_grid_evaluation_simd_impl( m_isp, xyzExtents, m_isInitializedVoxel, m_gridEvaluation,
                                                    m_outVoxelCornerValues, index );
    }

    get_density_from_grid_evaluation_simd_mt_body(
        const particle_zhu_bridson_is_policy& isp, const frantic::graphics2d::boundrect2& xyExtents, const int z,
        const std::vector<boost::uint8_t>& isInitializedVoxel,
        const std::vector<detail::zhu_bridson_grid_particle_simd>& gridEvaluation, float* outVoxelCornerValues )
        : m_isp( isp )
        , m_xyExtents( xyExtents )
        , m_z( z )
        , m_isInitializedVoxel( isInitializedVoxel )
        , m_gridEvaluation( gridEvaluation )
        , m_outVoxelCornerValues( outVoxelCornerValues ) {}
};

} // anonymous namespace

void particle_zhu_bridson_is_policy::get_density_from_grid_evaluation_simd(
    const frantic::graphics::boundbox3& xyzExtents, const std::vector<boost::uint8_t>& isInitializedVoxel,
    const std::vector<detail::zhu_bridson_grid_particle_simd>& gridEvaluation, float* outVoxelCornerValues ) const {
    get_density_from_grid_evaluation_simd_impl( *this, xyzExtents, isInitializedVoxel, gridEvaluation,
                                                outVoxelCornerValues );
}

void particle_zhu_bridson_is_policy::get_density_from_grid_evaluation_simd_mt(
    const frantic::graphics2d::boundrect2& xyExtents, const int z,
    const std::vector<boost::uint8_t>& isInitializedVoxel,
    const std::vector<detail::zhu_bridson_grid_particle_simd>& gridEvaluation, float* outVoxelCornerValues ) const {
    tbb::blocked_range<int> range( xyExtents.minimum().y, xyExtents.maximum().y + 1 );
    get_density_from_grid_evaluation_simd_mt_body body( *this, xyExtents, z, isInitializedVoxel, gridEvaluation,
                                                        outVoxelCornerValues );

#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif
}

/**
 *	Calculates the Density using the userData at a giving particle.
 *
 *	@param	userData				Expected zhu_bridson_density_user_data - contains particle
 *information, point locations, and variables to holds the density information
 *	@param	p						The particle to use for density calculation
 */
void particle_zhu_bridson_is_policy::get_density( void* userData, char* p ) {
    zhu_bridson_density_user_data* data = reinterpret_cast<zhu_bridson_density_user_data*>( userData );

    vector3f position = data->positionAccessor( p );
    float radius = data->radiusAccessor( p );

    float invKernelCompactSupportSquared = 1.f / data->kernelCompactSupportSquared;
    float distanceSquared, weight;
    frantic::graphics::vector3f hVector;

    if( data->getGradient ) {

        for( int i = 0; i < 6; i++ ) {
            if( i < 2 )
                hVector = vector3f( data->h, 0, 0 );
            else if( i < 4 )
                hVector = vector3f( 0, data->h, 0 );
            else
                hVector = vector3f( 0, 0, data->h );
            if( i % 2 != 0 )
                hVector = -hVector;

            distanceSquared = frantic::graphics::vector3f::distance_squared( data->worldLocation + hVector, position );
            weight = zhu_bridson_kernel( distanceSquared, data->kernelCompactSupportSquared,
                                         invKernelCompactSupportSquared );
            data->blendedPosition[i] += weight * position;
            data->blendedRadius[i] += weight * radius;
            data->totalWeight[i] += weight;
        }

    } else {
        distanceSquared = frantic::graphics::vector3f::distance_squared( data->worldLocation, position );
        weight =
            zhu_bridson_kernel( distanceSquared, data->kernelCompactSupportSquared, invKernelCompactSupportSquared );
        data->blendedPosition[0] += weight * position;
        data->blendedRadius[0] += weight * radius;
        data->totalWeight[0] += weight;
    }
}

/**
 *	Calculates the density for all points in the worldBounds to the worldLocation
 *
 *	@param	worldLocation			Find the density at this point.
 */
float particle_zhu_bridson_is_policy::get_density( const frantic::graphics::vector3f& worldLocation ) const {
    const float support = m_kernelCompactSupport;

    zhu_bridson_density_user_data data( worldLocation, m_radiusAccessor, m_positionAccessor,
                                        boost::math::pow<2>( support ) );
    m_particles.process_particles_in_range( &data, worldLocation, support, get_density );

    float distanceFromSurface;
    if( data.totalWeight[0] > 0 ) {
        data.blendedPosition[0] *= ( 1.f / data.totalWeight[0] );
        data.blendedRadius[0] *= ( 1.f / data.totalWeight[0] );
        distanceFromSurface =
            frantic::graphics::vector3f::distance( data.blendedPosition[0], worldLocation ) - data.blendedRadius[0];
    } else {
        data.totalWeight[0] = 0;
        distanceFromSurface = get_default_outside_distance();
    }

    return distanceFromSurface + isosurface_function_compensation( data.totalWeight[0] );
}

/**
 *	Calculates the densities at the six points surrounding worldLocation when h is applied in all directions. Then
 *uses those densities to estimate the gradient with standard central-difference formulas.
 *
 *	@param	worldLocation			Find the density at this point.
 *	@param	h						The small step to use.
 */
frantic::graphics::vector3f
particle_zhu_bridson_is_policy::get_gradient( const frantic::graphics::vector3f& worldLocation, float h ) {
    const float support = m_kernelCompactSupport + h;

    zhu_bridson_density_user_data data( worldLocation, m_radiusAccessor, m_positionAccessor,
                                        boost::math::pow<2>( support ), h, true );
    m_particles.process_particles_in_range( &data, worldLocation, support, get_density );

    float distanceFromSurface[6];

    // calculate density for all 6 surrounding points
    for( int i = 0; i < 6; i++ ) {
        if( data.totalWeight[i] > 0 ) {

            frantic::graphics::vector3f hVector( 0 );
            if( i < 2 )
                hVector = vector3f( h, 0, 0 );
            else if( i < 4 )
                hVector = vector3f( 0, h, 0 );
            else
                hVector = vector3f( 0, 0, h );

            if( i % 2 != 0 )
                hVector = -hVector;

            data.blendedPosition[i] *= ( 1.f / data.totalWeight[i] );
            data.blendedRadius[i] *= ( 1.f / data.totalWeight[i] );
            distanceFromSurface[i] =
                frantic::graphics::vector3f::distance( data.blendedPosition[i], worldLocation + hVector ) -
                data.blendedRadius[i];
        } else {
            data.totalWeight[i] = 0;
            distanceFromSurface[i] = get_default_outside_distance();
        }
    }

    float x = ( distanceFromSurface[0] + isosurface_function_compensation( data.totalWeight[0] ) ) -
              ( distanceFromSurface[1] + isosurface_function_compensation( data.totalWeight[1] ) );
    float y = ( distanceFromSurface[2] + isosurface_function_compensation( data.totalWeight[2] ) ) -
              ( distanceFromSurface[3] + isosurface_function_compensation( data.totalWeight[3] ) );
    float z = ( distanceFromSurface[4] + isosurface_function_compensation( data.totalWeight[4] ) ) -
              ( distanceFromSurface[5] + isosurface_function_compensation( data.totalWeight[5] ) );

    return vector3f( x, y, z ) / 2 * h;
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	in the provided array.  By sparsely, we mean that it calculates only those interactions
 *  as dictated by each particle and its radii, rather than by each dense grid point.
 *  The given coordinates are interpreted as being in the coord system already defined
 *	in the policy.
 *
 *	This function has been multithreaded using tbb.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A shared array that will be populated sparsely with data.
 */
void particle_zhu_bridson_is_policy::fill_sparse_voxel_corner_densities( const boundrect2& xyExtents, int z,
                                                                         float* outVoxelCornerValues,
                                                                         detail::zhu_bridson_sparse_data& data ) {

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_plane_evaluation( m_positionAccessor, m_radiusAccessor, m_kernelCompactSupport, voxelExtents,
                                           m_meshingVCS );

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_kernelCompactSupport );
    worldSearchBox.expand( 4 * m_meshingVCS.voxel_length() );

// start threads
#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, zhu_bridson_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, zhu_bridson_sparse_contribution_function );
#endif

    // Evaluate the implicit function
    if( data.useSimd ) {
        get_density_from_grid_evaluation_simd_mt( xyExtents, z, data.isInitializedVoxel, data.gridEvaluationSimd,
                                                  outVoxelCornerValues );
    } else {
        get_density_from_grid_evaluation_mt( xyExtents, z, data.isInitializedVoxel, data.gridEvaluation,
                                             outVoxelCornerValues );
    }
}

std::size_t
particle_zhu_bridson_is_policy::fill_voxel_corner_densities( const frantic::graphics::boundbox3& voxelExtents,
                                                             float* outVoxelCornerValues,
                                                             detail::zhu_bridson_sparse_data& data ) const {
    boundbox3f worldExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_block_evaluation( m_positionAccessor, m_radiusAccessor, m_kernelCompactSupport, voxelExtents,
                                           m_meshingVCS );

    // Compute the interactions
    boundbox3f worldSearchBox( worldExtents );
    // need an epsilon too ?
    worldSearchBox.expand( m_kernelCompactSupport );

    m_particles.process_particles_in_bounds( &data, worldSearchBox, zhu_bridson_sparse_contribution_function );

    // Evaluate the implicit function
    if( data.useSimd ) {
        get_density_from_grid_evaluation_simd( voxelExtents, data.isInitializedVoxel, data.gridEvaluationSimd,
                                               outVoxelCornerValues );
    } else {
        get_density_from_grid_evaluation( voxelExtents, data.isInitializedVoxel, data.gridEvaluation,
                                          outVoxelCornerValues );
    }

    return data.gridParticleCount;
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	in the provided array.  By sparsely, we mean that it calculates only those interactions
 *  as dictated by each particle and its radii, rather than by each dense grid point.
 *  The given coordinates are interpreted as being in the coord system already defined
 *	in the policy.
 *
 *	This function has been multithreaded using tbb.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	outVoxelCornerValues	A shared array that will be populated sparsely with data.
 */
void particle_zhu_bridson_is_policy::fill_sparse_voxel_corner_densities( const boundrect2& xyExtents, int z,
                                                                         float* outVoxelCornerValues, rle_plane& rlp,
                                                                         detail::zhu_bridson_sparse_data& data ) {

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    // Set up the sparse data structure with all the info we need to determine which grid points a particle
    // will interact with
    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.kernelCompactSupport = m_kernelCompactSupport;
    data.kernelCompactSupportSquared = m_kernelCompactSupport * m_kernelCompactSupport;
    data.voxelExtents = voxelExtents;
    data.meshingVCS = m_meshingVCS;

    for( int i = 0; i < xyExtents.ysize(); ++i )
        data.runTrees.push_back( run_tree( xyExtents.xsize() ) );

    // Set up the temp grid particle channel
    data.gridEvaluation.resize( xyExtents.get_area() );

    data.gridParticleChannel = 0;
    data.channelAccessors.clear();
    data.channelAccessors.push_back( frantic::channels::channel_general_accessor(
        0, sizeof( detail::zhu_bridson_grid_particle ), frantic::channels::data_type_int8 ) );
    data.channelData.clear();
    data.channelData.push_back( (char*)( &data.gridEvaluation[0] ) );
    data.initChannelData.clear();
    data.initChannelData.push_back( vector<char>( sizeof( detail::zhu_bridson_grid_particle ), 0 ) );

    // set up pointers to the extent data for each channel
    // data.extentData.clear();
    data.extentData = std::vector<std::vector<char*>>( xyExtents.ysize() );
    for( size_t y = 0; y < data.extentData.size(); ++y ) {
        data.extentData[y] = vector<char*>( data.channelData.size() );
        // create a vector of pointers to the current extents in the channel data and the grid particles
        size_t dataOffset = y * data.voxelExtents.xsize();
        for( size_t channelNum = 0; channelNum < data.channelData.size(); channelNum++ )
            data.extentData[y][channelNum] =
                data.channelData[channelNum] + data.channelAccessors[channelNum].primitive_size() * dataOffset;
    }
    data.collectRunData = true;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_kernelCompactSupport );

    // set up mutexes and start threads
    data.useMutexes = true;
    data.mutexes.reset( new tbb::spin_mutex[xyExtents.ysize()] );
#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, zhu_bridson_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, zhu_bridson_sparse_contribution_function );
#endif

    // Build the run indexing structure from the run trees
    rlp.reset( xyExtents );
    for( int y = 0; y < xyExtents.ysize(); ++y ) {

        vector<pair<int, int>> runs;

        data.runTrees[y].traverse_runs( runs );
        for( int i = 0; i < (int)runs.size(); ++i ) {
            runs[i].first += y * xyExtents.xsize();
            runs[i].second += y * xyExtents.xsize();
        }
        rlp.append_runs_by_extent( runs, y );
    }

    // Use it to traverse the grid point values and calculate the final channel values
    if( !rlp.is_empty() ) {
        for( int y = 0; y < rlp.get_xy_extents().ysize(); y++ ) {
            int xOffset = y * rlp.get_xy_extents().xsize();
            int firstRun = rlp.get_y_extent_first_run( y );
            if( firstRun >= 0 ) {
                int lastRun = rlp.get_y_extent_last_run( y );
                for( int run = firstRun; run <= lastRun; ++run ) {
                    // skip undefined runs
                    if( rlp.get_run_code( run ) < 0 )
                        continue;
                    for( int index = rlp[run].first; index <= rlp[run].second; ++index ) {
                        float weight = data.gridEvaluation[index].weight;
                        if( weight > 0.f ) {
                            vector3f blendedOffset = data.gridEvaluation[index].blendedOffset / weight;
                            float particleRadius = data.gridEvaluation[index].blendedRadius / weight;
                            vector3 cornerSample( xyExtents.minimum().x + ( index - xOffset ),
                                                  xyExtents.minimum().y + y, z );
                            outVoxelCornerValues[index] = blendedOffset.get_magnitude() - particleRadius +
                                                          isosurface_function_compensation( weight );
                        } else {
                            // Sometimes a 0 weight works its way in here, im not sure how atm.
                            // Eg, a couple of 0 weights in a 30 mb level set.  For the moment, I will
                            // zero the resulting channel information out, so I don't get any bad numerical
                            // values, but I need to determine how it gets there...
                            outVoxelCornerValues[index] = m_kernelCompactSupport;
                        }
                    }
                }
            }
        }
    }
    // sparsePM.exit_section( 7 );
}

/**
 *	This function takes the extents and z coord of a plane, and populates its cornervalues sparsely
 *	with the channel wise information placing it in the provided char arrays.  An rle plane structure
 *  which includes defined/undefined run information, is also populated to facilitate sparse access
 *  of the planar data.  The given coordinates are interpreted as being in the coord system already
 *  defined in the policy.  Only channels described in the channel propagation policy are propagated.
 *
 *	@param	xyExtents	A boundrect2 of the xyExtents of the plane
 *	@param	z			The z coord of the plane
 *	@param	channelName	A vector of names for desired channels to be copied.
 *	@param	channelData	A vector of shared arrays that will be populated sparsely with data.
 *	@param	rlp			An rle plane indexing structure to be populated to allow sparse access to the
 *						above planar data.
 */
void particle_zhu_bridson_is_policy::fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                                               std::vector<frantic::tstring>& channelNames,
                                                               std::vector<char*>& channelData,
                                                               frantic::volumetrics::rle_plane& rlp ) {

    if( channelNames.empty() )
        throw std::runtime_error( "particle_zhu_bridson_is_policy::fill_sparse_channel_data() - No channel specified, "
                                  "channelNames argument was empty." );

    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    // Set up the sparse data structure with all the info we need to determine which grid points a particle
    // will interact with
    detail::zhu_bridson_sparse_data data;

    data.position = m_positionAccessor;
    data.radius = m_radiusAccessor;
    data.kernelCompactSupport = m_kernelCompactSupport;
    data.kernelCompactSupportSquared = m_kernelCompactSupport * m_kernelCompactSupport;
    data.voxelExtents = voxelExtents;
    data.meshingVCS = m_meshingVCS;

    // set up the run tree data
    for( int i = 0; i < xyExtents.ysize(); ++i )
        data.runTrees.push_back( run_tree( xyExtents.xsize() ) );

    // set up the channel data and accessors
    int signedDistanceChannel = -1;
    for( size_t i = 0; i < channelNames.size(); ++i ) {
        // skip the signed distance channel, it gets special cased
        if( channelNames[i] == _T("SignedDistance") )
            signedDistanceChannel = (int)i;
        else {
            data.channelAccessors.push_back( m_particles.get_channel_map().get_general_accessor( channelNames[i] ) );
            data.channelData.push_back( channelData[i] );
            data.initChannelData.push_back(
                std::vector<char>( data.channelAccessors[data.channelAccessors.size() - 1].primitive_size(), 0 ) );
        }
    }

    // set up a channel for grid particles
    data.gridEvaluation.resize( xyExtents.get_area() );
    data.gridParticleChannel = (int)data.channelData.size();
    ;
    data.channelAccessors.push_back( frantic::channels::channel_general_accessor(
        0, sizeof( detail::zhu_bridson_grid_particle ), frantic::channels::data_type_int8 ) );
    data.channelData.push_back( (char*)( &data.gridEvaluation[0] ) );
    data.initChannelData.push_back( vector<char>( sizeof( detail::zhu_bridson_grid_particle ), 0 ) );

    // set up pointers to the extent data for each channel
    data.extentData = std::vector<std::vector<char*>>( xyExtents.ysize() );
    for( size_t y = 0; y < data.extentData.size(); ++y ) {
        data.extentData[y] = vector<char*>( data.channelData.size() );
        // create a vector of pointers to the current extents in the channel data and the grid particles
        size_t dataOffset = y * data.voxelExtents.xsize();
        for( size_t channelNum = 0; channelNum < data.channelData.size(); channelNum++ )
            data.extentData[y][channelNum] =
                data.channelData[channelNum] + data.channelAccessors[channelNum].primitive_size() * dataOffset;
    }
    data.collectRunData = true;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_kernelCompactSupport );

    // set up mutexes and start threads
    data.useMutexes = true;
    data.mutexes.reset( new tbb::spin_mutex[xyExtents.ysize()] );
// tbb::task_scheduler_init taskSchedulerInit; //please move out.
#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, zhu_bridson_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, zhu_bridson_sparse_contribution_function );
#endif

    // Build the run indexing structure from the run trees
    rlp.reset( xyExtents );

    for( int y = 0; y < xyExtents.ysize(); ++y ) {
        vector<pair<int, int>> runs;
        data.runTrees[y].traverse_runs( runs );
        int offset = y * xyExtents.xsize();

        if( runs.size() == 0 ) {
            rlp.append_runs_by_extent( runs, y, -1 );
            continue;
        }

        // I need to pop in undefined runs here, also.  Due to the nature of the particle IS policy,
        // I will only need outside undefined runs, and only between pairs of defined runs.
        vector<pair<int, int>> allRuns;
        allRuns.reserve( 2 * runs.size() );
        vector<int> allRunCodes;
        allRunCodes.reserve( 2 * runs.size() );
        size_t i;
        for( i = 0; i < runs.size() - 1; ++i ) {
            // defined run/code
            allRuns.push_back( std::pair<int, int>( runs[i].first + offset, runs[i].second + offset ) );
            allRunCodes.push_back( 0 );

            // undefined run/code
            allRuns.push_back( std::pair<int, int>( runs[i].second + 1 + offset, runs[i + 1].first - 1 + offset ) );
            allRunCodes.push_back( -1 );
        }
        allRuns.push_back( std::pair<int, int>( runs[i].first + offset, runs[i].second + offset ) );
        allRunCodes.push_back( 0 );

        // pop in the new runs/codes
        rlp.append_runs_by_extent( allRuns, allRunCodes, y );
    }

    // Use it to traverse the grid point values and calculate the final channel values
    if( !rlp.is_empty() ) {
        for( int y = 0; y < rlp.get_xy_extents().ysize(); y++ ) {
            int xOffset = y * rlp.get_xy_extents().xsize();
            int firstRun = rlp.get_y_extent_first_run( y );
            if( firstRun >= 0 ) {
                int lastRun = rlp.get_y_extent_last_run( y );
                for( int run = firstRun; run <= lastRun; ++run ) {
                    // skip undefined runs
                    if( rlp.get_run_code( run ) < 0 )
                        continue;
                    for( int index = rlp[run].first; index <= rlp[run].second; ++index ) {
                        float weight = data.gridEvaluation[index].weight;
                        if( weight > 0.f ) {
                            if( signedDistanceChannel >= 0 ) {
                                vector3f blendedOffset = data.gridEvaluation[index].blendedOffset / weight;
                                float particleRadius = data.gridEvaluation[index].blendedRadius / weight;
                                vector3 cornerSample( xyExtents.minimum().x + ( index - xOffset ),
                                                      xyExtents.minimum().y + y, z );
                                ( (float*)( channelData[signedDistanceChannel] ) )[index] =
                                    blendedOffset.get_magnitude() - particleRadius +
                                    isosurface_function_compensation( weight );
                            }
                            for( size_t channelNum = 0; channelNum < data.channelData.size(); ++channelNum ) {
                                std::vector<float> weights;
                                weights.push_back( 1.f / weight );
                                vector<char*> inData;
                                inData.push_back( data.channelData[channelNum] +
                                                  index * data.channelAccessors[channelNum].primitive_size() );
                                data.channelAccessors[channelNum].weighted_sum( weights, inData, inData[0] );
                            }
                        } else {
                            // Sometimes a 0 weight works its way in here, im not sure how atm.
                            // Eg, a couple of 0 weights in a 30 mb level set.  For the moment, I will
                            // zero the resulting channel information out, so I don't get any bad numerical
                            // values, but I need to determine how it gets there...
                            if( signedDistanceChannel >= 0 )
                                ( (float*)( channelData[signedDistanceChannel] ) )[index] = m_kernelCompactSupport;
                            for( size_t channelNum = 0; channelNum < data.channelData.size(); ++channelNum )
                                memset( data.channelData[channelNum] +
                                            index * data.channelAccessors[channelNum].primitive_size(),
                                        0, data.channelAccessors[channelNum].primitive_size() );
                        }
                    }
                }
            }
        }
    }
}

frantic::graphics::vector3f particle_zhu_bridson_is_policy::find_edge_isosurface_location(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::vector<char*>& particles,
    frantic::volumetrics::implicitsurface::detail::xyzr_packet_array& simdParticles ) const {

    particles.clear();
    simdParticles.clear();

    // SIMD is slower for smaller numbers of particles
    const std::size_t simdThreshold = SIMD_PARTICLE_THRESHOLD;

    bool useSimd = ( float_v::static_size > 1 );

    vector3f worldCorner0 = corner_sample_coord_to_world( voxelCorner0 );
    vector3f worldCorner1 = corner_sample_coord_to_world( voxelCorner1 );
    const int solveAxis = ( worldCorner0 - worldCorner1 ).get_largest_axis();

    vector3f surfaceLocation;

    // Refine the vertex position the requested number of times, to make it more closely approximate the isosurface.
    if( m_vertexRefinement > 0 ) {
        // Make a bounding box containing the extents of this voxel (this will always be an axis-aligned line)
        // This way, we only need to get the nearby particles once, though as a trade off we will get more particles
        // than we need for any individual evaluation.
        boundbox3f vertRange( worldCorner0 );
        vertRange += worldCorner1;
        m_particles.get_particles_in_range( vertRange, m_kernelCompactSupport, particles );
        if( particles.size() < simdThreshold ) {
            useSimd = false;
        }

        int iterCount = 0;
        float refined;
        if( useSimd ) {
            load( simdParticles, particles, m_positionAccessor, m_radiusAccessor );
            zhu_bridson_vert_refine_eval_simd solver( *this, worldCorner0, solveAxis, simdParticles,
                                                      m_kernelCompactSupport );
            refined = find_root( solver, worldCorner0[solveAxis], worldCorner1[solveAxis], density0, density1,
                                 m_vertexRefinementEpsilon, m_vertexRefinement, iterCount );
        } else {
            char** particlesBegin = particles.size() > 0 ? &particles[0] : 0;
            char** particlesEnd = particlesBegin + particles.size();

            zhu_bridson_vert_refine_eval solver( *this, worldCorner0, solveAxis, particlesBegin, particlesEnd,
                                                 m_positionAccessor, m_radiusAccessor, m_kernelCompactSupport );
            refined = find_root( solver, worldCorner0[solveAxis], worldCorner1[solveAxis], density0, density1,
                                 m_vertexRefinementEpsilon, m_vertexRefinement, iterCount );
        }

        surfaceLocation = worldCorner0;
        surfaceLocation[solveAxis] = refined;

        if( !useSimd ) {
            // Filter the particle array to only contain the nearby particles we need.
            // I'm skipping this in the SIMD case because it doesn't appear to improve
            // performance.
            const float kernelCompactSupportSquared = boost::math::pow<2>( m_kernelCompactSupport );
            for( unsigned i = 0; i < particles.size(); ) {
                float distanceSquared =
                    vector3f::distance_squared( surfaceLocation, m_positionAccessor( particles[i] ) );
                if( distanceSquared >= kernelCompactSupportSquared ) {
                    // Remove the particle at the back
                    particles[i] = particles.back();
                    particles.pop_back();
                } else {
                    ++i;
                }
            }
        }
    } else {
        surfaceLocation = worldCorner0;
        float alpha = fabsf( density0 / ( density1 - density0 ) );
        surfaceLocation[solveAxis] = ( 1 - alpha ) * worldCorner0[solveAxis] + alpha * worldCorner1[solveAxis];
    }

    return surfaceLocation;
}

/**
 *	This function finds the location of the isosurface between the two voxel corners provided,
 *	and sets the internal state of the policy.
 *	@param	density0		First voxel corner value
 *	@param	density1		Second voxel corner value
 *	@param	voxelCorner0	First voxel corner position
 *	@param	voxelCorner0	Second voxel corner position
 */
void particle_zhu_bridson_is_policy::find_edge_isosurface_location( float density0, float density1,
                                                                    const vector3& voxelCorner0,
                                                                    const vector3& voxelCorner1 ) {
    // TODO: Could move this variable into class scope so we can re-use its
    // memory between invocations.  But I don't think this function is in
    // active use so I'm not going to bother yet.
    frantic::volumetrics::implicitsurface::detail::xyzr_packet_array simdParticles;
    m_isosurfaceLocationWorld = find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1,
                                                               m_isosurfaceLocationParticles, simdParticles );
}

const vector3f& particle_zhu_bridson_is_policy::get_isosurface_location_world() const {
    return m_isosurfaceLocationWorld;
}

////////
// Functions for dealing with the named channels
////////

/**
 *  Set the internal state of the policy so that it points to the
 * specified location.  This is intended to be called before
 * add_vertex_to_channels if you have manually added a vertex at
 * the specified location, particularly as in generate_disambiguated_faces_for_plane.
 */
void particle_zhu_bridson_is_policy::set_isosurface_location_world( const frantic::graphics::vector3f& worldLocation ) {
    m_isosurfaceLocationWorld = worldLocation;
    m_isosurfaceLocationParticles.clear();
}

/**
 *	Given a set of writable trimesh3 accessors, this should add one data element to each of the channels,
 *	which corresponds the the channel names given in prepare_channel_accessors.  The function uses the internal
 *state of the mc policy for the location at which to look up the data
 */
void particle_zhu_bridson_is_policy::add_vertex_to_channels(
    std::vector<geometry::trimesh3_vertex_channel_general_accessor>& outputChannels ) {
    // If necessary, get the particles that are near the sample location
    if( m_isosurfaceLocationParticles.empty() ) {
        m_particles.get_particles_in_range( m_isosurfaceLocationWorld, m_kernelCompactSupport,
                                            m_isosurfaceLocationParticles );
    }

    if( outputChannels.size() > 0 ) {
        const std::size_t vertIndex = outputChannels[0].size();
        detail::zhu_bridson_vertex_workspace workspace;

        for( unsigned i = 0; i < outputChannels.size(); ++i ) {
            if( outputChannels[i].size() <= vertIndex ) {
                outputChannels[i].set_vertex_count( vertIndex + 1 );
            }
        }

        populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, m_isosurfaceLocationWorld,
                                  m_isosurfaceLocationParticles, workspace, m_channelMapWeightedSum );
    }
}

/**
 *	This function finds the location of the isosurface between the two voxel corners provided,
 *	and sets the internal state of the policy.
 *	@param	density0		First voxel corner value
 *	@param	density1		Second voxel corner value
 *	@param	voxelCorner0	First voxel corner position
 *	@param	voxelCorner0	Second voxel corner position
 *	@param	vertIndex		Index of the vert in the accompanying mesh
 *	@param	outputChannels	Trimesh3 channel accessors for any channel data also to be added
 *	@param	outMesh			The mesh to add the vert to.
 *	@param	meshMutex		A mutex for write locking if the function is used in a multithreaded context.
 */
void particle_zhu_bridson_is_policy::add_edge_isosurface_vertex_to_mesh(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
    frantic::geometry::trimesh3& outMesh, detail::zhu_bridson_vertex_workspace& workspace ) const {

    typedef detail::zhu_bridson_vertex_workspace::particles_t particles_t;
    particles_t& particles = workspace.particles;

    vector3f surfaceLocation = find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1, particles,
                                                              workspace.simdParticles );

    const bool gotSimdParticles = particles.size() == workspace.simdParticles.get_particle_count();

    // SIMD is slower for smaller numbers of particles
    const std::size_t simdThreshold = SIMD_PARTICLE_THRESHOLD;

    bool useSimd = ( float_v::static_size > 1 ) && gotSimdParticles;

    // add a vert to the mesh at the surface location
    outMesh.vertices_ref()[vertIndex] = surfaceLocation;

    if( outputChannels.empty() )
        return;

    // If we didnt have to fetch particles for the surface solve (ie, 0 vertex refinement)
    // we need to fetch them now.
    if( particles.empty() ) {
        m_particles.get_particles_in_range( surfaceLocation, m_kernelCompactSupport, particles );
        if( particles.size() < simdThreshold ) {
            useSimd = false;
        }
        if( useSimd ) {
            load( workspace.simdParticles, particles, m_positionAccessor, m_radiusAccessor );
        }
    }

    if( useSimd ) {
        get_sample_weights_simd( surfaceLocation, workspace.simdParticles, workspace.sampleWeights );
        filter_positive_weights( particles, workspace.sampleWeights );
        populate_vertex_channels_impl( m_inputChannels, outputChannels, vertIndex, workspace.sampleWeights, particles,
                                       workspace, m_channelMapWeightedSum );
    } else {
        populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, surfaceLocation, particles, workspace,
                                  m_channelMapWeightedSum );
    }
}

/**
 *	This function populates all vertex channels specified in the channel propagation
 *	policy using info from the surface policy.
 *	<p>
 *	All other channel info is preserved.
 *
 *	@param	outMesh		Mesh containing vertex and channel info.
 *	@param	cpp			A channel propagation policy.
 */
void particle_zhu_bridson_is_policy::populate_mesh_channels(
    frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_propagation_policy& cpp ) const {
    std::vector<frantic::tstring> previousNames;                // all channels names in the mesh
    std::vector<frantic::tstring> names;                        // all channels that are going to be populated
    std::vector<std::pair<std::size_t, data_type_t>> dataTypes; // types and arity of each channel to be populated
    std::vector<frantic::channels::channel_general_accessor> inputChannels; // channel accessors to channel_map
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>
        outputChannels; // channel accessors to vertex_channel

    // tbb::task_scheduler_init taskSchedulerInit; //please move out.

    // If channels already exist throw error.
    outMesh.get_vertex_channel_names( previousNames );
    for( unsigned i = 0; i < previousNames.size(); ++i )
        if( cpp.is_channel_included( previousNames[i] ) )
            throw runtime_error( "particle_zhu_bridson_is_policy.populate_mesh_channels: Attempt to populate a channel "
                                 "that already exists." );

    // Find the names of all channels found in both cpp and the particles.
    get_channel_names( cpp, names );

    // Prepare the input channel accessors and get the data types which
    // will be used to create new channels in the mesh.
    dataTypes.reserve( names.size() );
    inputChannels.reserve( names.size() );
    for( unsigned i = 0; i < names.size(); ++i ) {
        inputChannels.push_back( m_particles.get_channel_map().get_general_accessor( names[i] ) );
        dataTypes.push_back( std::make_pair( inputChannels.back().arity(), inputChannels.back().data_type() ) );
    }

    // Create the new vertex and get accessors to those vertex channels.
    for( size_t i = 0; i < names.size(); ++i ) {
        outMesh.add_vertex_channel_raw( names[i], dataTypes[i].first, dataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( names[i] ) );
    }

    channel_map_weighted_sum channelMapWeightedSum( m_particles.get_channel_map(), cpp );

    // Populate each vertex channel.
    tbb::parallel_for(
        tbb::blocked_range<size_t>( 0, outMesh.vertex_count() ),
        zhu_bridson_populate_vertex_channels( outputChannels, outMesh, inputChannels, this, channelMapWeightedSum ) );
}

void particle_zhu_bridson_is_policy::populate_vertex_channels(
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace ) const {
    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, vert, workspace, m_channelMapWeightedSum );
}

/**
 *	This function populates all vertex channels using data from the implicit
 *	surface policy for the input vertex.
 *
 *	@param	inputChannels	Accessors to particle channels used for population.
 *	@param	outputChannels	Accessors to each channel that should be populated.
 *	@param	vert			Location of the vertex.
 *	@param	vertIndex		Index of the vertex.
 */
void particle_zhu_bridson_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace,
    const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const {
    workspace.particles.clear();
    m_particles.get_particles_in_range( vert, m_kernelCompactSupport, workspace.particles );
    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, workspace.particles, workspace,
                              channelMapWeightedSum );
}

void particle_zhu_bridson_is_policy::populate_vertex_channels(
    std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert ) const {
    vertex_workspace_t workspace;
    workspace.particles.clear();
    m_particles.get_particles_in_range( vert, m_kernelCompactSupport, workspace.particles );
    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, workspace.particles, workspace,
                              m_channelMapWeightedSum );
}

void particle_zhu_bridson_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, const std::vector<char*>& particlesInRange, vertex_workspace_t& workspace,
    const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const {
    if( inputChannels.size() != outputChannels.size() )
        throw std::runtime_error( "particle_zhu_bridson_is_policy.populate_vertex_channels: inputChannels and "
                                  "outputChannels are not the same size" );

    // Compute the sample weights for blending between the different values
    get_sample_weights( vert, particlesInRange, workspace.sampleWeights );

    populate_vertex_channels_impl( inputChannels, outputChannels, vertIndex, workspace.sampleWeights, particlesInRange,
                                   workspace, channelMapWeightedSum );
}

namespace {

/**
 *	Data structures used internally for anistropic get_density function
 */
struct anistropic_density_user_data {
    bool getGradient;
    frantic::graphics::vector3f worldLocation;
    frantic::channels::channel_accessor<float> volumeAccessor;
    frantic::channels::channel_accessor<float> invCompactSupportVolumeAccessor;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor;
    frantic::channels::channel_general_accessor anisotropyAccessor;
    float implicitThreshold;
    float h;
    float densities[6];

    anistropic_density_user_data(
        const frantic::graphics::vector3f& worldLocation,
        const frantic::channels::channel_accessor<float>& volumeAccessor,
        const frantic::channels::channel_accessor<float>& invCompactSupportVolumeAccessor,
        const frantic::channels::channel_accessor<frantic::graphics::vector3f>& positionAccessor,
        const frantic::channels::channel_general_accessor& anisotropyAccessor, float implicitThreshold, float h = 0,
        bool getGradient = false )
        : worldLocation( worldLocation )
        , volumeAccessor( volumeAccessor )
        , invCompactSupportVolumeAccessor( invCompactSupportVolumeAccessor )
        , positionAccessor( positionAccessor )
        , anisotropyAccessor( anisotropyAccessor )
        , implicitThreshold( implicitThreshold )
        , getGradient( getGradient )
        , h( h ) {

        densities[0] = implicitThreshold;
        if( getGradient )
            for( int i = 1; i < 6; i++ )
                densities[i] = implicitThreshold;
    }
};

class anisotropic_sparse_contribution_evaluator
    : public sparse_contribution_evaluator<anisotropic_sparse_contribution_evaluator> {
  public:
    typedef detail::anisotropic_sparse_data sparse_data_t;

    typedef detail::anisotropic_sparse_data::kernel_t kernel_t;

    anisotropic_sparse_contribution_evaluator( const char* particle, sparse_data_t& data )
        : m_volume( data.volume( particle ) )
        , m_invCompactSupportVolume( data.invCompactSupportVolume( particle ) )
        , m_signedDistanceChannelNum( data.signedDistanceChannel )
        , m_signedDistanceChannel( (float*)( data.channelData[data.signedDistanceChannel] ) )
        , m_anisotropy( reinterpret_cast<const float*>( data.anisotropy.get_channel_data_pointer( particle ) ) ) {}

    BOOST_FORCEINLINE void evaluate_x_run( const char* particle, const vector3f& particlePosition, int xMin, int xMax,
                                           float particleOffsetY, float particleOffsetZ, float /*distanceSquaredYZ*/,
                                           int extentNum, int dataOffset, sparse_data_t& data ) {
        if( data.collectRunData ) {
            // Push the run into the appropriate run_tree.  This will initialize the data in if it isn't already.
            std::pair<int, int> run( xMin - data.voxelExtents.minimum().x, xMax - data.voxelExtents.minimum().x - 1 );
            if( run.second >= run.first )
                data.runTrees[extentNum].insert_run( run, data.extentData[extentNum], data.initChannelData );
        }

        for( int x = xMin; x <= xMax; x++ ) {
            float worldVoxelCenterX = data.meshingVCS.get_world_x_voxel_center( x );
            float particleOffsetX = particlePosition.x - worldVoxelCenterX;

            int i = x + dataOffset;

            const float initialSignedDistance = m_signedDistanceChannel[i];

            const frantic::graphics::vector3f dx( particleOffsetX, particleOffsetY, particleOffsetZ );
            const vector3f dxt = detail::smv33( m_anisotropy, dx );
            const float distanceSquared = dxt.get_magnitude_squared();

            const float weight =
                m_volume * kernel_t::kernel_distance_fraction2( distanceSquared, m_invCompactSupportVolume );
            m_signedDistanceChannel[i] = initialSignedDistance - weight;

            if( !data.collectRunData )
                continue;

            // add in the effect of this particle to this location in each channel
            size_t channelCount = data.channelData.size();
            for( size_t channelNum = 0; channelNum < channelCount; ++channelNum ) {

                // skip the grid particle channel
                if( channelNum == (size_t)m_signedDistanceChannelNum )
                    continue;

                // add in the weighted contribution of this particle for the channel
                frantic::channels::channel_general_accessor& ca = data.channelAccessors[channelNum];
                ca.weighted_increment( weight, ca.get_channel_data_pointer( particle ),
                                       data.channelData[channelNum] + i * ca.primitive_size() );
            }
        }
    }

  private:
    float m_volume;
    float m_invCompactSupportVolume;
    int m_signedDistanceChannelNum;
    float* m_signedDistanceChannel;
    const float* m_anisotropy;
};

void anisotropic_sparse_contribution_function( void* userData, char* treeParticle ) {
    detail::anisotropic_sparse_data* data = reinterpret_cast<detail::anisotropic_sparse_data*>( userData );

    const vector3f treeParticlePosition = data->position( treeParticle );
    const float treeParticleInteractionRadius = data->effectRadius( treeParticle );

    anisotropic_sparse_contribution_evaluator f( treeParticle, *data );
    f.evaluate( treeParticle, treeParticlePosition, treeParticleInteractionRadius, *data );

    if( f.has_particle_on_grid() ) {
        ++( data->gridParticleCount );
    }
}

class anisotropic_vert_refine_eval {
    frantic::graphics::vector3f m_voxelPosition0;

    int m_solveAxis;

    char** m_particlesBegin;
    char** m_particlesEnd;

    const particle_anisotropic_is_policy& m_anisotropicPolicy;

    anisotropic_vert_refine_eval& operator=( const anisotropic_vert_refine_eval& ); // not implemented

  public:
    anisotropic_vert_refine_eval( const particle_anisotropic_is_policy& anisotropicPolicy,
                                  const frantic::graphics::vector3f& voxelCoord0, const int solveAxis,
                                  char** particlesBegin, char** particlesEnd )
        : m_anisotropicPolicy( anisotropicPolicy )
        , m_voxelPosition0( voxelCoord0 )
        , m_solveAxis( solveAxis )
        , m_particlesBegin( particlesBegin )
        , m_particlesEnd( particlesEnd ) {}

    float operator()( float x ) const {
        vector3f vertTest( m_voxelPosition0 );
        vertTest[m_solveAxis] = x;

        // evaluate the density at this location
        return m_anisotropicPolicy.get_density( m_particlesBegin, m_particlesEnd, vertTest );
    }
};

class anisotropic_populate_vertex_channels {
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels;
    frantic::geometry::trimesh3& outMesh;
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels;
    const particle_anisotropic_is_policy* isp;
    const channel_map_weighted_sum& channelMapWeightedSum;

    anisotropic_populate_vertex_channels& operator=( const anisotropic_populate_vertex_channels& ); // not implemented

  public:
    anisotropic_populate_vertex_channels(
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
        frantic::geometry::trimesh3& outMesh,
        const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
        const particle_anisotropic_is_policy* isp, const channel_map_weighted_sum& channelMapWeightedSum )
        : outputChannels( outputChannels )
        , outMesh( outMesh )
        , inputChannels( inputChannels )
        , isp( isp )
        , channelMapWeightedSum( channelMapWeightedSum ) {}

    void operator()( const tbb::blocked_range<size_t>& r ) const {
        detail::anisotropic_vertex_workspace workspace;
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            isp->populate_vertex_channels( inputChannels, outputChannels, i, outMesh.get_vertex( i ), workspace,
                                           channelMapWeightedSum );
        }
    }
};

} // anonymous namespace

namespace detail {

void anisotropic_sparse_data::reset_for_dense_plane_evaluation(
    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
    frantic::channels::channel_accessor<float> radiusAccessor,
    frantic::channels::channel_accessor<float> volumeAccessor,
    frantic::channels::channel_accessor<float> maxDistanceAccessor,
    frantic::channels::channel_accessor<float> invCompactSupportVolumeAccessor,
    const frantic::channels::channel_general_accessor& anisotropyAccessor,
    const frantic::graphics::boundbox3& voxelExtents_, const frantic::volumetrics::voxel_coord_system& meshingVCS_,
    float* voxelCornerDensities_ ) {

    assert( voxelExtents_.zsize() == 1 );

    reset_for_dense_block_evaluation( positionAccessor, radiusAccessor, volumeAccessor, maxDistanceAccessor,
                                      invCompactSupportVolumeAccessor, anisotropyAccessor, voxelExtents_, meshingVCS_,
                                      voxelCornerDensities_ );

    // set up mutexes
    this->useMutexes = true;
    this->mutexes.reset( new tbb::spin_mutex[voxelExtents_.ysize()] );
}

void anisotropic_sparse_data::reset_for_dense_block_evaluation(
    frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAccessor,
    frantic::channels::channel_accessor<float> radiusAccessor,
    frantic::channels::channel_accessor<float> volumeAccessor,
    frantic::channels::channel_accessor<float> maxDistanceAccessor,
    frantic::channels::channel_accessor<float> invCompactSupportVolumeAccessor,
    const frantic::channels::channel_general_accessor& anisotropyAccessor,
    const frantic::graphics::boundbox3& voxelExtents_, const frantic::volumetrics::voxel_coord_system& meshingVCS_,
    float* voxelCornerDensities_ ) {

    // Set up the sparse data structure with all the info we need to determine which grid points a particle
    // will interact with
    this->position = positionAccessor;
    this->radius = radiusAccessor;
    this->volume = volumeAccessor;
    this->effectRadius = maxDistanceAccessor;
    this->invCompactSupportVolume = invCompactSupportVolumeAccessor;
    this->anisotropy = anisotropyAccessor;
    this->voxelExtents = voxelExtents_;
    this->meshingVCS = meshingVCS_;
    this->insideKernelCount = 0;
    this->outsideKernelCount = 0;

    // set up the signed distance channel
    this->signedDistanceChannel = 0;
    this->channelAccessors.clear();
    this->channelAccessors.push_back(
        frantic::channels::channel_general_accessor( 0, 1, frantic::channels::data_type_float32 ) );
    this->channelData.clear();
    this->channelData.push_back( (char*)voxelCornerDensities_ );
    this->initChannelData.clear();
    this->initChannelData.push_back( std::vector<char>( sizeof( float ) ) );
    this->gridParticleCount = 0;

    this->collectRunData = false;

    // not threaded
    this->mutexes.reset();
    this->useMutexes = false;
}

} // namespace detail

particle_anisotropic_is_policy::particle_anisotropic_is_policy(
    frantic::particles::particle_grid_tree& particles, float implicitThreshold,
    const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement )
    : particle_is_policy_base<particle_anisotropic_is_policy>( particles )
    , m_implicitThreshold( implicitThreshold )
    , m_meshingVCS( meshingVCS )
    , m_vertexRefinement( vertexRefinement ) {
    m_positionAccessor = m_particles.get_channel_map().get_accessor<vector3f>( _T("Position") );
    m_radiusAccessor = m_particles.get_channel_map().get_accessor<float>( _T("Radius") );
    m_volumeAccessor = m_particles.get_channel_map().get_accessor<float>( _T("__Volume") );
    m_invCompactSupportVolumeAccessor = m_particles.get_channel_map().get_accessor<float>( _T("__invcsv") );
    m_anisotropyAccessor = m_particles.get_channel_map().get_general_accessor( _T("__Anisotropy") );
    m_maxDistanceAccessor = m_particles.get_channel_map().get_accessor<float>( _T("__MaxDistance") );

    m_insideKernelCount = 0;
    m_outsideKernelCount = 0;

    m_maximumEffectRadius = 0;
    for( particle_grid_tree::iterator i = particles.begin(); i != particles.end(); ++i ) {
        if( m_maxDistanceAccessor( *i ) > m_maximumEffectRadius ) {
            m_maximumEffectRadius = m_maxDistanceAccessor( *i );
        }
    }

    // Compute the voxel bounds over which to do the implicit surface conversion
    boundbox3f particleWorldBounds = m_particles.compute_particle_bounds();
    particleWorldBounds.expand( m_maximumEffectRadius );
    m_particleVCSBounds = m_meshingVCS.get_voxel_bounds( particleWorldBounds );
    m_particleVCSBounds.expand( 1 );

    m_vertexRefinementEpsilon = 0.001f * m_meshingVCS.voxel_length();
}

particle_anisotropic_is_policy::~particle_anisotropic_is_policy() {
    /*
    int count = 0;
    int maxCount = 100;
    for( particle_grid_tree::iterator i = m_particles.begin(); i != m_particles.end(); ++i ) {
      FF_LOG( stats ) << boost::lexical_cast<std::string>( get_density( m_positionAccessor( * i ) ) ) << std::endl;
      ++count;
      if( count >= maxCount ) {
        break;
      }
    }
    */
    // FF_LOG( stats ) << "Inside kernel count: " << boost::lexical_cast<std::string>( m_insideKernelCount ) <<
    // std::endl; FF_LOG( stats ) << "Outside kernel count: " << boost::lexical_cast<std::string>( m_outsideKernelCount
    // ) << std::endl;
}

const voxel_coord_system& particle_anisotropic_is_policy::get_voxel_coord_system() const { return m_meshingVCS; }

frantic::graphics::boundbox3 particle_anisotropic_is_policy::get_voxel_bounds() const { return m_particleVCSBounds; }

frantic::graphics::vector3f particle_anisotropic_is_policy::corner_sample_coord_to_world(
    const frantic::graphics::vector3& voxelCornerCoord ) const {
    return m_meshingVCS.get_world_voxel_center( voxelCornerCoord );
}

frantic::graphics::vector3f particle_anisotropic_is_policy::find_edge_isosurface_location(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::vector<char*>& particles ) const {
    particles.clear();

    vector3f worldCorner0 = corner_sample_coord_to_world( voxelCorner0 );
    vector3f worldCorner1 = corner_sample_coord_to_world( voxelCorner1 );
    const int solveAxis = ( worldCorner0 - worldCorner1 ).get_largest_axis();

    vector3f surfaceLocation;

    // Refine the vertex position the requested number of times, to make it more closely approximate the isosurface.
    if( m_vertexRefinement > 0 ) {
        const float kernelCompactSupport = m_maximumEffectRadius;

        // Make a bounding box containing the extents of this voxel (this will always be an axis-aligned line)
        // This way, we only need to get the nearby particles once, though as a trade off we will get more particles
        // than we need for any individual evaluation.
        boundbox3f vertRange( worldCorner0 );
        vertRange += worldCorner1;
        m_particles.get_particles_in_range( vertRange, kernelCompactSupport, particles );

        char** particlesBegin = particles.size() > 0 ? &particles[0] : 0;
        char** particlesEnd = particlesBegin + particles.size();

        anisotropic_vert_refine_eval solver( *this, worldCorner0, solveAxis, particlesBegin, particlesEnd );

        int iterCount = 0;
        float refined = find_root( solver, worldCorner0[solveAxis], worldCorner1[solveAxis], density0, density1,
                                   m_vertexRefinementEpsilon, m_vertexRefinement, iterCount );

        surfaceLocation = worldCorner0;
        surfaceLocation[solveAxis] = refined;

        // Filter the particle array to only contain the nearby particles we need
        const float kernelCompactSupportSquared = boost::math::pow<2>( kernelCompactSupport );
        for( unsigned i = 0; i < particles.size(); ) {
            float distanceSquared = vector3f::distance_squared( surfaceLocation, m_positionAccessor( particles[i] ) );
            if( distanceSquared >= kernelCompactSupportSquared ) {
                // Remove the particle at the back
                particles[i] = particles.back();
                particles.pop_back();
            } else {
                ++i;
            }
        }
    } else {
        surfaceLocation = worldCorner0;
        float alpha = fabsf( density0 / ( density1 - density0 ) );
        surfaceLocation[solveAxis] = ( 1 - alpha ) * worldCorner0[solveAxis] + alpha * worldCorner1[solveAxis];
    }

    return surfaceLocation;
}

void particle_anisotropic_is_policy::add_edge_isosurface_vertex_to_mesh(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
    frantic::geometry::trimesh3& outMesh, vertex_workspace_t& workspace ) const {

    typedef vertex_workspace_t::particles_t particles_t;
    particles_t& particles = workspace.particles;

    particles.clear();

    vector3f surfaceLocation =
        find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1, particles );

    const float kernelCompactSupport =
        m_maximumEffectRadius; // m_maximumParticleRadius * m_effectRadiusScale * m_maxAnisotropyRadiusScale;

    {
        // add a vert to the mesh at the surface location
        // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
        outMesh.vertices_ref()[vertIndex] = surfaceLocation;
    }

    if( outputChannels.empty() )
        return;

    // If we didnt have to fetch particles for the surface solve (ie, 0 vertex refinement)
    // we need to fetch them now.
    if( particles.empty() ) {
        m_particles.get_particles_in_range( surfaceLocation, kernelCompactSupport, particles );
    }

    // In some cases with no vertex refinement, there may be no particles
    // in range of the selected vertex position.
    // Expand the search radius to pick up such distant particles.
    if( particles.empty() ) {
        m_particles.get_particles_in_range( surfaceLocation, kernelCompactSupport + m_meshingVCS.voxel_length(),
                                            particles );
    }

    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, surfaceLocation, particles, workspace,
                              m_channelMapWeightedSum );
}

/**
 *	Calculates the Density using the userData at a giving particle.
 *
 *	@param	userData				Expected anistropic_density_user_data - contains particle
 *information, point locations, and variables to holds the density information
 *	@param	p						The particle to use for density calculation
 */
void particle_anisotropic_is_policy::get_density( void* userData, char* p ) {
    anistropic_density_user_data* data = reinterpret_cast<anistropic_density_user_data*>( userData );

    const frantic::graphics::vector3f& position = data->positionAccessor( p );
    float volume = data->volumeAccessor( p );
    const float invCompactSupportVolume = data->invCompactSupportVolumeAccessor( p );
    const float* anisotropy = reinterpret_cast<const float*>( data->anisotropyAccessor.get_channel_data_pointer( p ) );
    vector3f dx = position - data->worldLocation;

    if( data->getGradient ) {
        vector3f dxh = dx - vector3f( data->h, 0, 0 );
        data->densities[0] -=
            volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dxh ).get_magnitude_squared(),
                                                          invCompactSupportVolume );

        dxh = dx + vector3f( data->h, 0, 0 );
        data->densities[1] -=
            volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dxh ).get_magnitude_squared(),
                                                          invCompactSupportVolume );

        dxh = dx - vector3f( 0, data->h, 0 );
        data->densities[2] -=
            volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dxh ).get_magnitude_squared(),
                                                          invCompactSupportVolume );

        dxh = dx + vector3f( 0, data->h, 0 );
        data->densities[3] -=
            volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dxh ).get_magnitude_squared(),
                                                          invCompactSupportVolume );

        dxh = dx - vector3f( 0, 0, data->h );
        data->densities[4] -=
            volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dxh ).get_magnitude_squared(),
                                                          invCompactSupportVolume );

        dxh = dx + vector3f( 0, 0, data->h );
        data->densities[5] -=
            volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dxh ).get_magnitude_squared(),
                                                          invCompactSupportVolume );
    } else {
        float f = volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dx ).get_magnitude_squared(),
                                                                invCompactSupportVolume );
        data->densities[0] -= f;
    }
}

float particle_anisotropic_is_policy::get_density( const char* particle,
                                                   const frantic::graphics::vector3f& testPosition ) const {
    const vector3f position = m_positionAccessor( particle );
    const vector3f dx = position - testPosition;

    const float* anisotropy =
        reinterpret_cast<const float*>( m_anisotropyAccessor.get_channel_data_pointer( particle ) );
    const float volume = m_volumeAccessor( particle );
    const float invCompactSupportVolume = m_invCompactSupportVolumeAccessor( particle );

    return volume * kernel_t::kernel_distance_fraction2( detail::smv33( anisotropy, dx ).get_magnitude_squared(),
                                                         invCompactSupportVolume );
}

float particle_anisotropic_is_policy::get_density( char** particlesBegin, char** particlesEnd,
                                                   const frantic::graphics::vector3f& testPosition ) const {
    float density = m_implicitThreshold;
    for( char** i = particlesBegin; i != particlesEnd; ++i ) {
        density -= get_density( *i, testPosition );
    }
    return density;
}

/**
 *	Calculates the density for all points in the worldBounds to the worldLocation
 *
 *	@param	worldLocation			Find the density at this point.
 */
float particle_anisotropic_is_policy::get_density( const frantic::graphics::vector3f& worldLocation ) const {
    float support = m_maximumEffectRadius;

    anistropic_density_user_data data( worldLocation, m_volumeAccessor, m_invCompactSupportVolumeAccessor,
                                       m_positionAccessor, m_anisotropyAccessor, m_implicitThreshold );
    m_particles.process_particles_in_range( &data, worldLocation, support, get_density );

    return data.densities[0];
}

/**
 *	Calculates the densities at the six points surrounding worldLocation when h is applied in all directions. Then
 *uses those densities to estimate the gradient with standard central-difference formulas.
 *
 *	@param	worldLocation			Find the density at this point.
 *	@param	h						The small step to use.
 */
frantic::graphics::vector3f
particle_anisotropic_is_policy::get_gradient( const frantic::graphics::vector3f& worldLocation, float h ) {
    float support = m_maximumEffectRadius + h;

    anistropic_density_user_data data( worldLocation, m_volumeAccessor, m_invCompactSupportVolumeAccessor,
                                       m_positionAccessor, m_anisotropyAccessor, m_implicitThreshold, h, true );
    m_particles.process_particles_in_range( &data, worldLocation, support, get_density );

    float x = data.densities[0] - data.densities[1];
    float y = data.densities[2] - data.densities[3];
    float z = data.densities[4] - data.densities[5];

    return vector3f( x, y, z ) / 2 * h;
}

/**
 *	This function populates all vertex channels specified in the channel propagation
 *	policy using info from the surface policy.
 *	<p>
 *	All other channel info is preserved.
 *
 *	@param	outMesh		Mesh containing vertex and channel info.
 *	@param	cpp			A channel propagation policy.
 */
void particle_anisotropic_is_policy::populate_mesh_channels(
    frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_propagation_policy& cpp ) const {
    std::vector<frantic::tstring> previousNames;                // all channels names in the mesh
    std::vector<frantic::tstring> names;                        // all channels that are going to be populated
    std::vector<std::pair<std::size_t, data_type_t>> dataTypes; // types and airty of each channel to be populated
    std::vector<frantic::channels::channel_general_accessor> inputChannels; // channel accessors to vertex_channel
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>
        outputChannels; // channel accessors to vertex_channel

    // If channels already exist, throw error.
    outMesh.get_vertex_channel_names( previousNames );
    for( size_t i = 0; i < previousNames.size(); ++i )
        if( cpp.is_channel_included( previousNames[i] ) )
            throw runtime_error( "particle_anisotropic_is_policy.populate_mesh_channels: Attempt to populate a channel "
                                 "that already exists." );

    // Find the names of all channels found in both cpp and the particles.
    get_channel_names( cpp, names );

    // Prepare the input channel accessors and get the data types which
    // will be used to create new channels in the mesh.
    dataTypes.reserve( names.size() );
    for( size_t i = 0; i < names.size(); ++i ) {
        inputChannels.push_back( m_particles.get_channel_map().get_general_accessor( names[i] ) );
        dataTypes.push_back( std::make_pair( inputChannels.back().arity(), inputChannels.back().data_type() ) );
    }

    // Create the new vertex and get accessors to those vertex channels.
    for( size_t i = 0; i < names.size(); ++i ) {
        outMesh.add_vertex_channel_raw( names[i], dataTypes[i].first, dataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( names[i] ) );
    }

    channel_map_weighted_sum channelMapWeightedSum( m_particles.get_channel_map(), cpp );

    // Populate each vertex channel.
    // tbb::task_scheduler_init taskSchedulerInit; //please move out.
    tbb::parallel_for(
        tbb::blocked_range<size_t>( 0, outMesh.vertex_count() ),
        anisotropic_populate_vertex_channels( outputChannels, outMesh, inputChannels, this, channelMapWeightedSum ) );
}

void particle_anisotropic_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert ) const {
    if( inputChannels.size() != outputChannels.size() )
        throw std::runtime_error( "particle_anisotropic_is_policy.populate_vertex_channels: inputChannels and "
                                  "outputChannels are not the same size" );

    // Get the particles that are near the sample location
    std::vector<char*> particlesInRange;
    m_particles.get_particles_in_range( vert, m_maximumEffectRadius, particlesInRange );

    vertex_workspace_t vertexWorkspace;
    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, particlesInRange, vertexWorkspace,
                              m_channelMapWeightedSum );
}

void particle_anisotropic_is_policy::populate_vertex_channels(
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace ) const {
    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, vert, workspace, m_channelMapWeightedSum );
}

void particle_anisotropic_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& workspace,
    const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const {
    workspace.particles.clear();
    m_particles.get_particles_in_range( vert, m_maximumEffectRadius, workspace.particles );
    populate_vertex_channels( inputChannels, outputChannels, vertIndex, vert, workspace.particles, workspace,
                              channelMapWeightedSum );
}

void particle_anisotropic_is_policy::populate_vertex_channels(
    const std::vector<frantic::channels::channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, const std::vector<char*>& particlesInRange, vertex_workspace_t& workspace,
    const frantic::channels::channel_map_weighted_sum& channelMapWeightedSum ) const {
    vector<float>& sampleWeights = workspace.sampleWeights;

    // Compute the sample weights for blending between the different values
    get_sample_weights( vert, particlesInRange, sampleWeights );

    populate_vertex_channels_impl( inputChannels, outputChannels, vertIndex, sampleWeights, particlesInRange, workspace,
                                   channelMapWeightedSum );
}

void particle_anisotropic_is_policy::get_sample_weights( const vector3f& position, const std::vector<char*>& particles,
                                                         std::vector<float>& outWeights ) const {
    const size_t particleCount = particles.size();

    outWeights.resize( particleCount );

    if( particleCount == 0 ) {
        return;
    }

    // Compute the sample weights for blending between the different values
    outWeights.resize( particles.size() );
    float weightSum = 0;
    for( size_t i = 0; i < particleCount; ++i ) {
        float weight = get_density( particles[i], position );
        outWeights[i] = weight;
        weightSum += weight;
    }

    // TODO: this is a hack to avoid infs in the output channel when there is
    // no particle in range.  What should we do in this case?
    if( weightSum == 0 ) {
        size_t nearestParticleIndex = 0;
        float nearestParticleDistance2 = std::numeric_limits<float>::max();
        for( size_t i = 0; i < particleCount; ++i ) {
            const vector3f& particlePosition = m_positionAccessor( particles[i] );
            const float distance2 = vector3f::distance_squared( position, particlePosition );
            if( distance2 < nearestParticleDistance2 ) {
                nearestParticleIndex = i;
                nearestParticleDistance2 = distance2;
            }
            outWeights[i] = 0.f;
        }
        outWeights[nearestParticleIndex] = 1.f;
        weightSum = 1.f;
    }

    // Normalize the sample weights
    for( size_t i = 0; i < particleCount; ++i ) {
        outWeights[i] /= weightSum;
    }
}

void particle_anisotropic_is_policy::get_affected_blocks(
    const frantic::graphics::size3f blockSize, std::vector<frantic::graphics::vector3>& outAffectedBlockCoordinates,
    shared_progress_logger_proxy& progressLogger ) {
    // do we need to look at every particle in this case?

    using boost::int32_t;
    using frantic::particles::particle_grid_tree;

    outAffectedBlockCoordinates.clear();

    boost::unordered_set<frantic::graphics::vector3, frantic::graphics::vector3_hasher> blockSet;

    // We need to add some length onto the max distance to ensure that we
    // get an outside voxel.  However I'm not sure if this is the right
    // value (for example, will floating point precision cause problems
    // away from the origin?).
    // It seems like it may be more appropriate to use integer math
    // and add one voxel instead.
    const float voxelDiagonalLength = 1.8f * m_meshingVCS.voxel_length();

    for( particle_grid_tree::const_node_iterator nodeIter = m_particles.const_nodes_begin(),
                                                 nodeIterEnd = m_particles.const_nodes_end();
         nodeIter != nodeIterEnd; ++nodeIter ) {
        bool first = true;
        boundbox3f bounds;

        particle_grid_tree::const_node_iterator::const_particle_iterator_pair particleIteratorPair =
            nodeIter.const_particle_iterators();
        for( const_particle_array_iterator i = particleIteratorPair.first, ie = particleIteratorPair.second; i != ie;
             ++i ) {
            const char* particle = *i;
            const vector3f position = m_positionAccessor.get( particle );

            boundbox3f particleBounds( position );
            particleBounds.expand( m_maxDistanceAccessor( particle ) + voxelDiagonalLength );

            if( first ) {
                bounds = particleBounds;
                first = false;
            } else {
                bounds += particleBounds;
            }
        }

        if( !bounds.is_empty() ) {
            const int32_t zbegin = static_cast<int32_t>( floor( bounds.minimum().z / blockSize.zsize() ) );
            const int32_t zend = static_cast<int32_t>( floor( bounds.maximum().z / blockSize.zsize() ) );
            const int32_t ybegin = static_cast<int32_t>( floor( bounds.minimum().y / blockSize.ysize() ) );
            const int32_t yend = static_cast<int32_t>( floor( bounds.maximum().y / blockSize.ysize() ) );
            const int32_t xbegin = static_cast<int32_t>( floor( bounds.minimum().x / blockSize.xsize() ) );
            const int32_t xend = static_cast<int32_t>( floor( bounds.maximum().x / blockSize.xsize() ) );
            for( int32_t z = zbegin; z <= zend; ++z ) {
                for( int32_t y = ybegin; y <= yend; ++y ) {
                    for( int32_t x = xbegin; x <= xend; ++x ) {
                        blockSet.insert( vector3( x, y, z ) );
                        if( progressLogger.is_cancelled() ) {
                            return;
                        }
                    }
                }
            }
        }
    }

    outAffectedBlockCoordinates.insert( outAffectedBlockCoordinates.begin(), blockSet.begin(), blockSet.end() );
}

void particle_anisotropic_is_policy::fill_sparse_voxel_corner_densities(
    const frantic::graphics2d::boundrect2& xyExtents, int z, float* outVoxelCornerValues,
    sparse_voxel_corner_density_workspace_t& data ) {
    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    boundbox3 voxelExtents( xyExtents.minimum().x, xyExtents.maximum().x, xyExtents.minimum().y, xyExtents.maximum().y,
                            z, z );
    boundbox3f xyzExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_plane_evaluation( m_positionAccessor, m_radiusAccessor, m_volumeAccessor,
                                           m_maxDistanceAccessor, m_invCompactSupportVolumeAccessor,
                                           m_anisotropyAccessor, voxelExtents, m_meshingVCS, outVoxelCornerValues );

    for( int i = 0, ie = xyExtents.get_area(); i != ie; ++i )
        outVoxelCornerValues[i] = m_implicitThreshold;

    // Compute the interactions
    boundbox3f worldSearchBox( xyzExtents );
    worldSearchBox.expand( m_maximumEffectRadius );

#ifndef FRANTIC_DISABLE_THREADS
    m_particles.process_particles_in_bounds_mt( &data, worldSearchBox, anisotropic_sparse_contribution_function );
#else
#pragma message( "Threads are disabled" )
    m_particles.process_particles_in_bounds( &data, worldSearchBox, anisotropic_sparse_contribution_function );
#endif

    m_insideKernelCount += data.insideKernelCount;
    m_outsideKernelCount += data.outsideKernelCount;
}

std::size_t
particle_anisotropic_is_policy::fill_voxel_corner_densities( const frantic::graphics::boundbox3& voxelExtents,
                                                             float* voxelCornerDensities,
                                                             sparse_voxel_corner_density_workspace_t& data ) const {
    // The incoming plane is in the meshing coord system, so it needs to be converted to world coords
    // so we can determine which particles will interact with it.
    const boundbox3f worldExtents = m_meshingVCS.get_world_bounds_of_voxel_centers( voxelExtents );

    data.reset_for_dense_block_evaluation( m_positionAccessor, m_radiusAccessor, m_volumeAccessor,
                                           m_maxDistanceAccessor, m_invCompactSupportVolumeAccessor,
                                           m_anisotropyAccessor, voxelExtents, m_meshingVCS, voxelCornerDensities );

    for( int i = 0, ie = voxelExtents.get_volume(); i != ie; ++i )
        voxelCornerDensities[i] = m_implicitThreshold;

    // Compute the interactions
    boundbox3f worldSearchBox( worldExtents );
    worldSearchBox.expand( m_maximumEffectRadius );

    m_particles.process_particles_in_bounds( &data, worldSearchBox, anisotropic_sparse_contribution_function );

    m_insideKernelCount += data.insideKernelCount;
    m_outsideKernelCount += data.outsideKernelCount;

    return data.gridParticleCount;
}
} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic
