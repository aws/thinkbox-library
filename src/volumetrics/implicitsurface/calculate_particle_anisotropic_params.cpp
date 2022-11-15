// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/pow.hpp>

#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_grid_tree.hpp>

#include <frantic/volumetrics/implicitsurface/calculate_particle_anisotropic_params.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>

#include <frantic/math/eigen.hpp>

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

frantic::channels::channel_map
create_channel_map_with_anisotropy_channels( const frantic::channels::channel_map& originalChannelMap,
                                             const frantic::tstring& volumeChannelName ) {
    frantic::channels::channel_map result( originalChannelMap );

    result.append_channel<float>( volumeChannelName );
    result.append_channel<float>( _T("__MaxDistance") );
    result.append_channel( _T("__Anisotropy"), 6, frantic::channels::data_type_float32 );
    // short for "inverse compact support volume"
    result.append_channel<float>( _T("__invcsv") );

    return result;
}

void calculate_anisotropy( frantic::particles::particle_array& particles, float compactSupportScale,
                           float anisotropyWindowScale, float kr, std::size_t nEps,
                           frantic::logging::progress_logger& progressLogger ) {
    frantic::volumetrics::implicitsurface::shared_progress_logger_adapter uiAdapter( progressLogger );
    frantic::volumetrics::implicitsurface::calculate_anisotropy_params params(
        particles, compactSupportScale, anisotropyWindowScale, kr, nEps, *uiAdapter.get_proxy() );
    boost::function<void( void )> f(
        boost::bind( frantic::volumetrics::implicitsurface::calculate_anisotropy, boost::ref( params ) ) );
    uiAdapter.run( f );
}

void calculate_volume_with_anisotropic_kernel( frantic::particles::particle_array& particles,
                                               const frantic::tstring& volumeChannelName,
                                               frantic::logging::progress_logger& progressLogger ) {
    frantic::volumetrics::implicitsurface::shared_progress_logger_adapter uiAdapter( progressLogger );
    frantic::volumetrics::implicitsurface::calculate_volume_with_anisotropic_kernel_params params(
        particles, volumeChannelName, *uiAdapter.get_proxy() );
    boost::function<void( void )> f( boost::bind(
        frantic::volumetrics::implicitsurface::calculate_volume_with_anisotropic_kernel, boost::ref( params ) ) );
    uiAdapter.run( f );
}

void smooth_particle_positions( frantic::particles::particle_array& particles, float effectRadiusScale, float lambda,
                                frantic::logging::progress_logger& progressLogger ) {
    frantic::volumetrics::implicitsurface::shared_progress_logger_adapter uiAdapter( progressLogger );
    frantic::volumetrics::implicitsurface::smooth_particle_positions_params params( particles, effectRadiusScale,
                                                                                    lambda, *uiAdapter.get_proxy() );
    boost::function<void( void )> f(
        boost::bind( frantic::volumetrics::implicitsurface::smooth_particle_positions, boost::ref( params ) ) );
    uiAdapter.run( f );
}

namespace detail {

// the kernel given in the reference paper has a typo -- it's given as 1 - (\|x_i - x_j\|)/r_i)^3

// 1 - ( (\x_i - x_j\|)/r_i )^3
struct diffusion_kernel {
    static float get_density2( const frantic::graphics::vector3f& xa, const frantic::graphics::vector3f& xb,
                               const float supportDistance2, const float invSupportDistance2 ) {
        const float distance2 = frantic::graphics::vector3f::distance_squared( xa, xb );
        if( distance2 < supportDistance2 ) {
            const float f = invSupportDistance2 * distance2;
            return 1.f - f * sqrt( f ); // pow( invSupportDistance2 * distance2, 1.5f );
        } else {
            return 0;
        }
    }
    static float get_uniform_cov( const float supportDistance2 ) { return 0.15f * supportDistance2; }
};

/*
// ( 1 - ( (\x_i - x_j\|)/r_i )^2 )^3
struct diffusion_kernel {
  static float get_density2( const frantic::graphics::vector3f & xa, const frantic::graphics::vector3f & xb, const float
supportDistance2, const float invSupportDistance2 ) { const float distance2 =
frantic::graphics::vector3f::distance_squared( xa, xb ); if( distance2 < supportDistance2 ) { return
boost::math::pow<3>( 1.f - ( distance2 * invSupportDistance2 ) ); } else { return 0;
    }
  }
  // diagonal entry of covariance matrix given uniform distribution
  static float get_uniform_cov( const float supportDistance2 ) {
    return supportDistance2 / 11.f;
  }
};
*/

struct particle_info_t {
    frantic::graphics::vector3f position;
    char* particle;
    std::size_t index;

    particle_info_t()
        : particle( 0 )
        , index( 0 )
        , position( 0 ) {}
    particle_info_t( std::size_t i, char* particle, const frantic::graphics::vector3f& position )
        : particle( particle )
        , index( i )
        , position( position ) {}
};

template <class ParticlePairProcessor, class voxel_info_t, class map_t>
class process_particle_pairs_in_range_impl {

    ParticlePairProcessor& m_processor;
    map_t& m_map;
    std::vector<frantic::graphics::vector3>& m_voxelCoords;
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& m_progressLogger;

    process_particle_pairs_in_range_impl& operator=( const process_particle_pairs_in_range_impl& ); // not implemented

  public:
    process_particle_pairs_in_range_impl(
        ParticlePairProcessor& processor, /*const frantic::particles::particle_array & particles,*/ map_t& map,
        std::vector<frantic::graphics::vector3>& voxelCoords,
        frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger )
        : m_processor( processor )
        , m_map( map )
        , m_voxelCoords( voxelCoords )
        , m_progressLogger( progressLogger ) {}
    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        boost::array<voxel_info_t*, 26> neighbourVoxels;
        typename ParticlePairProcessor::particle_state_t particleState;
        boost::array<vector3, 26> neighbourOffset = boost::assign::list_of( vector3( -1, -1, -1 ) )(
            vector3( 0, -1, -1 ) )( vector3( +1, -1, -1 ) )( vector3( -1, 0, -1 ) )( vector3( 0, 0, -1 ) )(
            vector3( +1, 0, -1 ) )( vector3( -1, +1, -1 ) )( vector3( 0, +1, -1 ) )( vector3( +1, +1, -1 ) )( vector3(
            -1, -1, 0 ) )( vector3( 0, -1, 0 ) )( vector3( +1, -1, 0 ) )( vector3( -1, 0, 0 ) )( vector3( +1, 0, 0 ) )(
            vector3( -1, +1, 0 ) )( vector3( 0, +1, 0 ) )( vector3( +1, +1, 0 ) )( vector3( -1, -1, +1 ) )(
            vector3( 0, -1, +1 ) )( vector3( +1, -1, +1 ) )( vector3( -1, 0, +1 ) )( vector3( 0, 0, +1 ) )(
            vector3( +1, 0, +1 ) )( vector3( -1, +1, +1 ) )( vector3( 0, +1, +1 ) )( vector3( +1, +1, +1 ) );

        for( std::size_t voxelNumber = r.begin(); voxelNumber != r.end(); ++voxelNumber ) {
            if( m_progressLogger.is_cancelled() ) {
                break;
            }

            const frantic::graphics::vector3 voxelCoord = m_voxelCoords[voxelNumber];
            typename map_t::iterator centerIter = m_map.find( voxelCoord );
            if( centerIter == m_map.end() ) {
                continue;
            }
            voxel_info_t& center = centerIter->second;

            for( std::size_t neighbourIndex = 0; neighbourIndex < neighbourOffset.size(); ++neighbourIndex ) {
                typename map_t::iterator neighbour = m_map.find( voxelCoord + neighbourOffset[neighbourIndex] );
                if( neighbour != m_map.end() ) {
                    neighbourVoxels[neighbourIndex] = &neighbour->second;
                } else {
                    neighbourVoxels[neighbourIndex] = 0;
                }
            }

            for( particle_info_t* pParticleInfo = center.begin; pParticleInfo != center.end; ++pParticleInfo ) {
                particle_info_t& particleInfo = *pParticleInfo;

                if( m_progressLogger.is_cancelled() ) {
                    break;
                }

                const std::size_t particleIndex = particleInfo.index;
                char* centerParticle = particleInfo.particle;

                m_processor.start( particleState, particleIndex, centerParticle /*m_particles[particleIndex]*/ );
                for( particle_info_t* pOtherParticleInfo = center.begin; pOtherParticleInfo != center.end;
                     ++pOtherParticleInfo ) {
                    particle_info_t& otherParticleInfo = *pOtherParticleInfo;
                    if( m_progressLogger.is_cancelled() ) {
                        break;
                    }
                    m_processor.process( particleState, otherParticleInfo.particle,
                                         otherParticleInfo.position /*m_particles[center[j].index]*/ );
                }

                BOOST_FOREACH( voxel_info_t* pNeighbourVoxel, neighbourVoxels ) {
                    if( pNeighbourVoxel ) {
                        voxel_info_t& neighbourVoxel = ( *pNeighbourVoxel );
                        for( particle_info_t *pOtherParticleInfo = neighbourVoxel.begin,
                                             *pOtherParticleInfoEnd = neighbourVoxel.end;
                             pOtherParticleInfo != pOtherParticleInfoEnd; ++pOtherParticleInfo ) {
                            particle_info_t& otherParticleInfo = *pOtherParticleInfo;
                            if( m_progressLogger.is_cancelled() ) {
                                break;
                            }
                            m_processor.process( particleState, otherParticleInfo.particle,
                                                 otherParticleInfo.position /*m_particles[neighbourVoxel[k].index]*/ );
                        }
                    }
                }
                m_processor.end( particleState );

                m_progressLogger.inc_progress();
            }
        }
    }
};

namespace {

#pragma warning( push )
// TODO : fix this instead of suppressing?
#pragma warning( disable : 4608 ) // has already been initialized by another union member in the initializer list
struct voxel_info_t {
    union {
        std::size_t particleCount;
        particle_info_t* end;
    };
    particle_info_t* begin;
    voxel_info_t()
        : begin( 0 ) {
        particleCount = 0;
        end = 0;
    }
};
#pragma warning( pop )

}; // anonymous namespace

// ParticlePairProcessor
// particle_state_t
// void start( particle_state_t & state, std::size_t particleNumber, const char * particle )
// void process( particle_state_t & state, const char * other )
// void end( particle_state_t & state )
template <class ParticlePairProcessor>
void process_particle_pairs_in_range(
    frantic::particles::particle_array& particles, const float maxDistance, ParticlePairProcessor& processor,
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger ) {
    tbb::task_scheduler_init taskScheduleInit;

    typedef std::vector<particle_info_t> voxel_particles_t;
    typedef boost::unordered_map<frantic::graphics::vector3, voxel_info_t, voxel_coord_hasher> map_t;

    map_t map;

    frantic::channels::channel_cvt_accessor<float> radiusAcc(
        particles.get_channel_map().get_cvt_accessor<float>( _T("Radius") ) );
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> positionAcc(
        particles.get_channel_map().get_cvt_accessor<frantic::graphics::vector3f>( _T("Position") ) );

    if( particles.size() == 0 ) {
        return;
    }

    frantic::volumetrics::voxel_coord_system vcs( frantic::graphics::vector3f( 0 ), maxDistance );

    for( std::size_t i = 0; i < particles.size(); ++i ) {
        frantic::graphics::vector3f pos = positionAcc( particles[i] );
        frantic::graphics::vector3f voxelCoord = vcs.get_voxel_coord( pos );
        frantic::graphics::vector3 voxelCoordInt = vector3::from_floor( voxelCoord );

        voxel_info_t& info = map[voxelCoordInt];
        ++info.particleCount;
    }

    std::vector<particle_info_t> particleInfo( particles.size() );

    std::vector<frantic::graphics::vector3> voxelCoords;
    voxelCoords.reserve( map.size() );
    std::size_t particleCumSum = 0;
    std::vector<std::size_t> voxelParticleIndex( map.size() );
    for( typename map_t::iterator i( map.begin() ); i != map.end(); ++i ) {
        voxelCoords.push_back( i->first );
        const std::size_t voxelParticleCount = i->second.particleCount;
        i->second.begin = i->second.end = &particleInfo[particleCumSum];
        particleCumSum += voxelParticleCount;
    }

    for( std::size_t i = 0; i < particles.size(); ++i ) {
        frantic::graphics::vector3f pos = positionAcc( particles[i] );
        frantic::graphics::vector3f voxelCoord = vcs.get_voxel_coord( pos );
        frantic::graphics::vector3 voxelCoordInt = vector3::from_floor( voxelCoord );

        voxel_info_t& info = map[voxelCoordInt];
        ( *info.end ) = particle_info_t( i, particles[i], positionAcc( particles[i] ) );
        ++info.end;
    }

    progressLogger.set_progress_end( particles.size() );

    process_particle_pairs_in_range_impl<ParticlePairProcessor, voxel_info_t, map_t> body( processor, map, voxelCoords,
                                                                                           progressLogger );
    tbb::blocked_range<std::size_t> range( 0, voxelCoords.size() );
#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif

    if( progressLogger.is_cancelled() ) {
        throw frantic::logging::progress_cancel_exception( "Particle-Particle Interactions" );
    }
}

template <class KernelType>
class smooth_particle_positions_impl {
    frantic::channels::channel_const_cvt_accessor<float> m_radiusAcc;
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector3f> m_positionAcc;
    std::vector<frantic::graphics::vector3f>& m_newPosition;
    float m_lambda;
    float m_effectRadiusScale;
    KernelType& m_kernel;

    smooth_particle_positions_impl& operator=( const smooth_particle_positions_impl& ); // not implemented
  public:
    struct smooth_particle_positions_particle_state {
        frantic::graphics::vector3f particlePosition;
        float effectRadius2;
        float invEffectRadius2;
        std::size_t particleNumber;
        const char* particle;
        frantic::graphics::vector3fd newPosition;
        double weightSum;

        smooth_particle_positions_particle_state() { reset( 0, 0, frantic::graphics::vector3f( 0 ), 1.f ); }

        void reset( std::size_t particleNumber, const char* particle,
                    const frantic::graphics::vector3f& particlePosition, float effectRadius ) {
            this->particlePosition = particlePosition;
            this->particleNumber = particleNumber;
            this->particle = particle;
            this->effectRadius2 = boost::math::pow<2>( effectRadius );
            this->invEffectRadius2 = 1.f / this->effectRadius2;
            newPosition.set( 0 );
            weightSum = 0;
        }
    };

    typedef smooth_particle_positions_particle_state particle_state_t;

    smooth_particle_positions_impl( const frantic::channels::channel_map& channelMap, float effectRadiusScale,
                                    float lambda, KernelType& kernel,
                                    std::vector<frantic::graphics::vector3f>& outNewPosition )
        : m_radiusAcc( channelMap.get_const_cvt_accessor<float>( _T("Radius") ) )
        , m_positionAcc( channelMap.get_const_cvt_accessor<frantic::graphics::vector3f>( _T("Position") ) )
        , m_effectRadiusScale( effectRadiusScale )
        , m_lambda( lambda )
        , m_newPosition( outNewPosition )
        , m_kernel( kernel ) {}

    void start( particle_state_t& state, std::size_t particleNumber, char* particle ) {
        state.reset( particleNumber, particle, m_positionAcc( particle ),
                     m_effectRadiusScale * m_radiusAcc( particle ) );
    }

    void process( particle_state_t& state, const char* /*other*/, const frantic::graphics::vector3f& otherPosition ) {
        if( state.effectRadius2 > 0 ) {
            const float w = m_kernel.get_density2( state.particlePosition, otherPosition, state.effectRadius2,
                                                   state.invEffectRadius2 );
            state.newPosition += frantic::graphics::vector3fd( w * otherPosition );
            state.weightSum += w;
        }
    }

    void end( particle_state_t& state ) {
        if( state.weightSum > 0 ) {
            const frantic::graphics::vector3fd temp = ( 1.0 / state.weightSum ) * state.newPosition;
            const frantic::graphics::vector3f weightedPosition( (float)temp.x, (float)temp.y, (float)temp.z );
            const frantic::graphics::vector3f newPosition =
                ( 1.f - m_lambda ) * state.particlePosition + m_lambda * weightedPosition;
            m_newPosition[state.particleNumber] = newPosition;
        } else {
            m_newPosition[state.particleNumber] = state.particlePosition;
        }
    }
};

template <class T>
class matrix33 {
    T m_elements[3][3];

  public:
    matrix33() { set_to_identity(); }
    matrix33( T a, T b, T c, T d, T e, T f, T g, T h, T i ) {
        m_elements[0][0] = a;
        m_elements[0][1] = b;
        m_elements[0][2] = c;
        m_elements[1][0] = d;
        m_elements[1][1] = e;
        m_elements[1][2] = f;
        m_elements[2][0] = g;
        m_elements[2][1] = h;
        m_elements[2][2] = i;
    }
    void set_to_identity() { set_to_diagonal( 1, 1, 1 ); }
    void set_to_diagonal( T a, T b, T c ) {
        memset( m_elements, 0, 9 * sizeof( T ) );
        m_elements[0][0] = a;
        m_elements[1][1] = b;
        m_elements[2][2] = c;
    }
    void set_from_columns( T a[3], T b[3], T c[3] ) {
        m_elements[0][0] = a[0];
        m_elements[0][1] = b[0];
        m_elements[0][2] = c[0];
        m_elements[1][0] = a[1];
        m_elements[1][1] = b[1];
        m_elements[1][2] = c[1];
        m_elements[2][0] = a[2];
        m_elements[2][1] = b[2];
        m_elements[2][2] = c[2];
    }
    T get( int row, int col ) const { return m_elements[row][col]; }
    T& get( int row, int col ) { return m_elements[row][col]; }
    void set( int row, int col, T val ) { m_elements[row][col] = val; }

    static matrix33<T> from_diagonal( T a, T b, T c ) {
        matrix33<T> m;
        m.set_to_diagonal( a, b, c );
        return m;
    }

    matrix33<T> transpose() {
        matrix33<T> m;
        for( int row = 0; row < 3; ++row ) {
            for( int col = 0; col < 3; ++col ) {
                m.set( row, col, m_elements[col][row] );
            }
        }
        return m;
    }

    matrix33<T>& operator*=( T a ) {
        for( int row = 0; row < 3; ++row ) {
            for( int col = 0; col < 3; ++col ) {
                m_elements[row][col] *= a;
            }
        }
        return *this;
    }

    matrix33<T>& operator=( const matrix33<T>& other ) {
        if( &other != this ) {
            memcpy( &m_elements[0], &other.m_elements[0], 9 * sizeof( T ) );
        }
        return *this;
    }
};

template <class T>
matrix33<T> operator*( const matrix33<T>& lhs, const matrix33<T>& rhs ) {
    matrix33<T> out;
    for( int row = 0; row < 3; ++row ) {
        for( int col = 0; col < 3; ++col ) {
            T acc = 0;
            for( int i = 0; i < 3; ++i ) {
                acc += lhs.get( row, i ) * rhs.get( i, col );
            }
            out.set( row, col, acc );
        }
    }
    return out;
}

template <class T>
matrix33<T> operator*( const T a, const matrix33<T>& m ) {
    return matrix33<T>( m ) *= a;
}

typedef matrix33<float> matrix33f;

template <class KernelType>
class anisotropy_calculator_impl {
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector3f> m_positionAcc;
    frantic::channels::channel_const_cvt_accessor<float> m_radiusAcc;
    frantic::channels::channel_accessor<float> m_maxSupportDistanceAcc;
    frantic::channels::channel_accessor<float> m_invCompactSupportVolumeAcc;
    frantic::channels::channel_general_accessor m_anisotropyAcc;
    KernelType& m_kernel;
    float m_compactSupportScale;
    float m_anisotropyRadiusScale;
    float m_kr;
    float m_kn;
    // float m_ks;
    std::size_t m_neps;
    // float m_maxmaxSigma;
    // float m_minmaxSigma;
    // float m_minvol;
    // float m_maxvol;

    anisotropy_calculator_impl& operator=( const anisotropy_calculator_impl& ); // not implemented
  public:
    struct particle_anisotropy_state {
        frantic::graphics::vector3f centrePosition;
        // float invScaledRadius2;
        float anisotropyWindowRadius2;
        float compactSupport2;
        std::size_t neighbourCount;
        float compactSupport;
        float anisotropyWindowRadius;
        std::size_t centreParticleIndex;
        std::vector<frantic::graphics::vector3f> neighbourPositions;
        std::vector<float> weights;
        char* particle;

        void reset( std::size_t index, const frantic::graphics::vector3f& position, char* particle,
                    const float scaledRadius, const float anisotropyWindowRadius ) {
            this->centreParticleIndex = index;
            this->centrePosition = position;
            this->compactSupport = scaledRadius;
            this->compactSupport2 = boost::math::pow<2>( scaledRadius );
            this->anisotropyWindowRadius = anisotropyWindowRadius;
            this->anisotropyWindowRadius2 = boost::math::pow<2>( anisotropyWindowRadius );
            this->particle = particle;
            this->neighbourPositions.clear();
            this->weights.clear();
            this->neighbourCount = 0;
        }

        float get_max_distance_squared() { return anisotropyWindowRadius2; }
        float get_neighbourhood_distance_squared() { return compactSupport2; }
    };

    typedef particle_anisotropy_state particle_state_t;

    anisotropy_calculator_impl( const frantic::channels::channel_map& channelMap, float effectRadiusScale,
                                float anisotropyRadiusScale, float kr, float kn, /*float ks,*/ std::size_t nEps,
                                KernelType& kernel, frantic::channels::channel_general_accessor& anisotropyAcc,
                                frantic::channels::channel_accessor<float>& maxSupportDistanceAcc,
                                frantic::channels::channel_accessor<float>& invCompactSupportVolumeAcc )
        : m_positionAcc( channelMap.get_const_cvt_accessor<frantic::graphics::vector3f>( _T("Position") ) )
        , m_radiusAcc( channelMap.get_const_cvt_accessor<float>( _T("Radius") ) )
        , m_compactSupportScale( effectRadiusScale )
        , m_anisotropyRadiusScale( anisotropyRadiusScale )
        , m_maxSupportDistanceAcc( maxSupportDistanceAcc )
        , m_invCompactSupportVolumeAcc( invCompactSupportVolumeAcc )
        , m_kr( kr )
        , m_kn( kn )
        ,
        // m_ks( ks ),
        m_neps( nEps )
        , m_kernel( kernel )
        , m_anisotropyAcc( anisotropyAcc ) {
        if( m_anisotropyAcc.arity() != 6 || m_anisotropyAcc.data_type() != frantic::channels::data_type_float32 ) {
            throw std::runtime_error(
                "calculate_anisotropy Error: wrong data type for anisotropy channel.  Wanted float32[6], but got " +
                frantic::strings::to_string(
                    frantic::channels::channel_data_type_str( m_anisotropyAcc.arity(), m_anisotropyAcc.data_type() ) ) +
                " instead." );
        }

        // m_maxmaxSigma = 0;
        // m_minmaxSigma = std::numeric_limits<float>::max();
        // m_minvol = std::numeric_limits<float>::max();
        // m_maxvol = 0;
    }

    void start( particle_state_t& state, std::size_t particleNumber, char* particle ) {
        const float radius = m_radiusAcc( particle );
        state.reset( particleNumber, m_positionAcc( particle ), particle, m_compactSupportScale * radius,
                     m_anisotropyRadiusScale * radius );
    }

    void process( particle_state_t& state, const char* /*other*/, const frantic::graphics::vector3f& otherPosition ) {
        const float distance2 = frantic::graphics::vector3f::distance_squared( otherPosition, state.centrePosition );
        if( distance2 < state.get_max_distance_squared() ) {
            state.neighbourPositions.push_back( otherPosition );
            if( distance2 < state.get_neighbourhood_distance_squared() ) {
                ++state.neighbourCount;
            }
        }
    }

    void end( particle_state_t& state ) {
        using boost::math::pow;

        const frantic::graphics::vector3f centrePosition = state.centrePosition;
        const float centreCompactSupport = state.compactSupport;
        const float invCentreCompactSupport = 1.f / centreCompactSupport;

        const float anisotropyWindowRadius = state.anisotropyWindowRadius;
        const float anisotropyWindowRadius2 = pow<2>( anisotropyWindowRadius );
        const float invAnisotropyWindowRadius = 1.f / anisotropyWindowRadius;
        const float invAnisotropyWindowRadius2 = pow<2>( invAnisotropyWindowRadius );

        frantic::graphics::vector3f sigma( m_kn );
        matrix33f R; // identity

        const std::size_t windowCount = state.neighbourPositions.size();

        const std::size_t nepsNeighbourhood =
            static_cast<std::size_t>( m_neps * pow<3>( centreCompactSupport ) / pow<3>( anisotropyWindowRadius ) );

        if( state.neighbourCount >= nepsNeighbourhood && windowCount >= m_neps ) {
            state.weights.resize( windowCount );
            double weightSum = 0;

            frantic::graphics::vector3f* neighbourPositions = windowCount > 0 ? &state.neighbourPositions[0] : 0;
            float* weights = windowCount > 0 ? &state.weights[0] : 0;

            for( size_t i = 0; i < windowCount; ++i ) {
                const float w = m_kernel.get_density2( centrePosition, neighbourPositions[i], anisotropyWindowRadius2,
                                                       invAnisotropyWindowRadius2 );
                weights[i] = w;
                weightSum += w;
            }

            if( weightSum > 0 ) {
                for( size_t i = 0; i < windowCount; ++i ) {
                    weights[i] = static_cast<float>( weights[i] / weightSum );
                }

                frantic::graphics::vector3f weightedMean( 0 );
                for( std::size_t i = 0; i < windowCount; ++i ) {
                    weightedMean += weights[i] * neighbourPositions[i];
                }

                // weightedMean = centrePosition;

                // C = sum_j  wij ( xj - xiw ) ( xj - xiw )^T / sum_j wij
                float c[6] = { 0, 0, 0, 0, 0, 0 };
                for( std::size_t i = 0; i < windowCount; ++i ) {
                    const float w = weights[i];
                    frantic::graphics::vector3f d = neighbourPositions[i] - weightedMean;
                    c[0] += w * pow<2>( d.x );
                    c[1] += w * d.x * d.y;
                    c[2] += w * d.x * d.z;
                    c[3] += w * pow<2>( d.y );
                    c[4] += w * d.y * d.z;
                    c[5] += w * pow<2>( d.z );
                }

                float ra[3], rb[3], rc[3]; // eigenvectors
                frantic::graphics::vector3f eigenValues;
                frantic::math::linearalgebra::eigendecompose_symmetric_3x3<float>( c, ra, rb, rc, &eigenValues.x );

                const float maxSigma = eigenValues.max_abs_component();
                const float minSigma = maxSigma / m_kr;
                // const float ks = m_ks;
                const float ks = 1.f / m_kernel.get_uniform_cov( anisotropyWindowRadius2 );

                // if( ks * maxSigma > m_kn ) {
                if( maxSigma > 0 ) {
                    // use cluster shape
                    for( int axis = 0; axis < 3; ++axis ) {
                        sigma[axis] = ks * std::max( eigenValues[axis], minSigma );
                    }

                    // TODO: symmetric so we don't need full multiplication etc
                    // G_i = 1/hi * R * inv( sigma )  * R^T;
                    R = matrix33f( ra[0], rb[0], rc[0], ra[1], rb[1], rc[1], ra[2], rb[2], rc[2] );
                }
            }
        }

        // m_maxmaxSigma = max( m_maxmaxSigma, sigma.max_abs_component() );
        // m_minmaxSigma = min( m_minmaxSigma, sigma.max_abs_component() );
        // m_minvol = min( m_minvol, sigma.x * sigma.y * sigma.z );
        // m_maxvol = max( m_maxvol, sigma.x * sigma.y * sigma.z );

        matrix33f G = invCentreCompactSupport * R * matrix33f::from_diagonal( 1 / sigma.x, 1 / sigma.y, 1 / sigma.z ) *
                      R.transpose();

        float* a = reinterpret_cast<float*>( m_anisotropyAcc.get_channel_data_pointer( state.particle ) );

        a[0] = G.get( 0, 0 );
        a[1] = G.get( 1, 0 );
        a[2] = G.get( 2, 0 );
        a[3] = G.get( 1, 1 );
        a[4] = G.get( 2, 1 );
        a[5] = G.get( 2, 2 );

        const float maxSupportDistance = centreCompactSupport * sigma.max_abs_component();
        m_maxSupportDistanceAcc( state.particle ) = maxSupportDistance;
        const float compactSupportVolume = 4.f / 3 * boost::math::constants::pi<float>() *
                                           boost::math::pow<3>( centreCompactSupport ) * sigma.x * sigma.y * sigma.z;
        const float invCompactSupportVolume = compactSupportVolume > 0 ? 1.f / compactSupportVolume : 1.f;
        m_invCompactSupportVolumeAcc( state.particle ) = invCompactSupportVolume;
#if 0
		{
			typedef detail::anisotropic_poly6_kernel KernelType;
			using boost::math::pow;
			using frantic::graphics::vector3fd;

			const vector3f position = centrePosition;
			const float testEdgeLength = 4.f * centreScaledRadius;
			const std::size_t testCount = 100000000;
			std::size_t insideCount = 0;
			double weightsum = 0;
			float maxObservedDistance = 0.f;
			for( std::size_t i = 0; i < testCount; ++i ) {
				const vector3f x = ( vector3f::from_random() - vector3f( 0.5f ) ) * testEdgeLength + position;
				const frantic::graphics::vector3f dx = position - x;
				const vector3f dxt = smv33( a, dx );
				const float distanceFraction2 = dxt.get_magnitude_squared();
				const float w = KernelType::kernel_distance_fraction2( distanceFraction2, invCompactSupportVolume );
				weightsum += w;

				if( w > 0 ) {
					++insideCount;
					const float distance = ( x - position ).get_magnitude();
					if( distance > maxObservedDistance ) {
						maxObservedDistance = distance;
					}
				}
			}
			const float observedVolume = float( insideCount ) / testCount * boost::math::pow<3>( testEdgeLength );
			const float integral = weightsum / testCount * boost::math::pow<3>( testEdgeLength );
			float dummy = 1.f;
		}
#endif
    }
};

template <class KernelType>
class volume_calculator_with_anisotropic_kernel_impl {
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAcc;
    frantic::channels::channel_accessor<float> m_radiusAcc;
    frantic::channels::channel_accessor<float> m_volumeAcc;
    frantic::channels::channel_accessor<float> m_invCompactSupportVolumeAcc;
    frantic::channels::channel_general_accessor m_anisotropyAcc;
    KernelType m_kernel;
    std::vector<float>& m_outVolume;

    volume_calculator_with_anisotropic_kernel_impl&
    operator=( const volume_calculator_with_anisotropic_kernel_impl& ); // not implemented

  public:
    struct volume_calculator_particle_state {
        std::size_t particleNumber;
        const char* particle;
        frantic::graphics::vector3f particlePosition;
        double density;
        float radius;

        volume_calculator_particle_state() { reset( 0, 0, frantic::graphics::vector3f( 0 ), 0 ); }

        void reset( std::size_t particleNumber, const char* particle,
                    const frantic::graphics::vector3f& particlePosition, float particleRadius ) {
            this->particleNumber = particleNumber;
            this->particle = particle;
            this->particlePosition = particlePosition;
            this->radius = particleRadius;
            this->density = 0;
        }

        float get_radius_squared() { return boost::math::pow<2>( radius ); }
    };

    typedef volume_calculator_particle_state particle_state_t;

    volume_calculator_with_anisotropic_kernel_impl( const frantic::channels::channel_map& channelMap,
                                                    KernelType& kernel, std::vector<float>& outVolume )
        : m_positionAcc( channelMap.get_accessor<frantic::graphics::vector3f>( _T("Position") ) )
        , m_radiusAcc( channelMap.get_accessor<float>( _T("Radius") ) )
        , m_volumeAcc( channelMap.get_accessor<float>( _T("__Volume") ) )
        , m_anisotropyAcc( channelMap.get_general_accessor( _T("__Anisotropy") ) )
        , m_invCompactSupportVolumeAcc( channelMap.get_accessor<float>( _T("__invcsv") ) )
        , m_kernel( kernel )
        , m_outVolume( outVolume ) {}

    void start( particle_state_t& state, std::size_t particleNumber, char* particle ) {
        state.reset( particleNumber, particle, m_positionAcc( particle ), m_radiusAcc( particle ) );
    }

    void process( particle_state_t& state, const char* other, const frantic::graphics::vector3f& otherPosition ) {
        const frantic::graphics::vector3f distVector = state.particlePosition - otherPosition;
        const float* otherAnisotropy =
            reinterpret_cast<const float*>( m_anisotropyAcc.get_channel_data_pointer( other ) );
        const float distanceSquared = smv33( otherAnisotropy, distVector ).get_magnitude_squared();
        const float otherRadius = m_radiusAcc( other );

        state.density += get_particle_mass( otherRadius ) *
                         m_kernel.kernel_distance_fraction2( distanceSquared, m_invCompactSupportVolumeAcc( other ) );
    }

    void end( particle_state_t& state ) {
        if( state.density > 0 ) {
            m_outVolume[state.particleNumber] = static_cast<float>( get_particle_mass( state.radius ) / state.density );
        } else {
            m_outVolume[state.particleNumber] = 0;
        }
    }
};

// density
// using the SPH kernel
// \rho_i = \sum_j m_j W( x_j - x_i, h_j )
// \rho_i = \sum_j r^3 W( x_j - x_i, h_j )

} // namespace detail

void calculate_anisotropy( frantic::particles::particle_array& particles, float compactSupportScale,
                           float anisotropyWindowScale, float kr, /*float kn, float ks,*/ std::size_t nEps,
                           frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger ) {
    // const float kn = 0.5f;
    const float kn = 1.33f / kr;
    frantic::channels::channel_cvt_accessor<float> radiusAcc(
        particles.get_channel_map().get_cvt_accessor<float>( _T("Radius") ) );
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> positionAcc(
        particles.get_channel_map().get_cvt_accessor<frantic::graphics::vector3f>( _T("Position") ) );
    frantic::channels::channel_accessor<float> maxSupportDistanceAcc =
        particles.get_channel_map().get_accessor<float>( _T("__MaxDistance") );
    frantic::channels::channel_accessor<float> invCompactSupportVolumeAcc =
        particles.get_channel_map().get_accessor<float>( _T("__invcsv") );
    frantic::channels::channel_general_accessor anisotropyAcc(
        particles.get_channel_map().get_general_accessor( _T("__Anisotropy") ) );
    if( anisotropyAcc.arity() != 6 || anisotropyAcc.data_type() != frantic::channels::data_type_float32 ) {
        throw std::runtime_error(
            "calculate_anisotropy Error: wrong data type for \'__Anisotropy\' channel.  Wanted float32[6], but got " +
            frantic::strings::to_string(
                frantic::channels::channel_data_type_str( anisotropyAcc.arity(), anisotropyAcc.data_type() ) ) +
            " instead." );
    }

    float maxRadius = 0;
    for( std::size_t i = 0; i < particles.size(); ++i ) {
        const float r = radiusAcc( particles[i] );
        if( r > maxRadius ) {
            maxRadius = r;
        }
    }

    const float maxDistance = anisotropyWindowScale * maxRadius;

    if( maxDistance > 0 ) {
        typedef detail::diffusion_kernel KernelType;
        KernelType kernel;

        detail::anisotropy_calculator_impl<KernelType> processor(
            particles.get_channel_map(), compactSupportScale, anisotropyWindowScale, kr, kn, nEps, kernel,
            anisotropyAcc, maxSupportDistanceAcc, invCompactSupportVolumeAcc );
        detail::process_particle_pairs_in_range( particles, maxDistance, processor, progressLogger );
    } else {
        float defaultAnisotropy[6] = { 1.f, 0, 0, 1.f, 0, 1.f };
        for( std::size_t i = 0; i < particles.size(); ++i ) {
            memcpy( anisotropyAcc.get_channel_data_pointer( particles[i] ), &defaultAnisotropy[0],
                    6 * sizeof( float ) );
        }
    }
}

void calculate_volume_with_anisotropic_kernel(
    frantic::particles::particle_array& particles, const frantic::tstring& volumeChannelName,
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger ) {
    typedef detail::anisotropic_poly6_kernel KernelType;
    typedef detail::volume_calculator_with_anisotropic_kernel_impl<KernelType> processor_t;

    using frantic::channels::channel_accessor;
    using frantic::channels::channel_general_accessor;

    channel_accessor<float> radiusAcc( particles.get_channel_map().get_accessor<float>( _T("Radius") ) );
    channel_accessor<frantic::graphics::vector3f> positionAcc(
        particles.get_channel_map().get_accessor<frantic::graphics::vector3f>( _T("Position") ) );
    channel_accessor<float> volumeAcc( particles.get_channel_map().get_accessor<float>( volumeChannelName ) );
    channel_general_accessor anisotropyAcc( particles.get_channel_map().get_general_accessor( _T("__Anisotropy") ) );
    if( anisotropyAcc.arity() != 6 || anisotropyAcc.data_type() != frantic::channels::data_type_float32 ) {
        throw std::runtime_error(
            "calculate_anisotropy Error: wrong data type for \'__Anisotropy\' channel.  Wanted float32[6], but got " +
            frantic::strings::to_string(
                frantic::channels::channel_data_type_str( anisotropyAcc.arity(), anisotropyAcc.data_type() ) ) +
            " instead." );
    }
    frantic::channels::channel_accessor<float> maxDistanceAcc =
        particles.get_channel_map().get_accessor<float>( _T("__MaxDistance") );

    float maxRadius = 0;
    for( std::size_t i = 0; i < particles.size(); ++i ) {
        const float r = maxDistanceAcc( particles[i] );
        if( r > maxRadius ) {
            maxRadius = r;
        }
    }
    const float maxDistance = maxRadius;

    KernelType kernel;
    std::vector<float> volume( particles.size(), 0 );

    if( maxRadius > 0 ) {
        processor_t processor( particles.get_channel_map(), kernel, volume );
        detail::process_particle_pairs_in_range( particles, maxDistance, processor, progressLogger );
    }

    for( std::size_t i = 0; i < particles.size(); ++i ) {
        volumeAcc( particles[i] ) = volume[i];
    }
}

// KernelType must have a member function .get_density( vector3f x1, vector3f x2, float radius1 )
// .get_effect_radius_scale()
template <class KernelType>
void smooth_particle_positions( frantic::particles::particle_array& particles, KernelType& kernel,
                                float effectRadiusScale, float lambda,
                                frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger ) {
    frantic::channels::channel_cvt_accessor<float> radiusAcc(
        particles.get_channel_map().get_cvt_accessor<float>( _T("Radius") ) );
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> positionAcc(
        particles.get_channel_map().get_cvt_accessor<frantic::graphics::vector3f>( _T("Position") ) );

    float maxRadius = 0;
    for( std::size_t i = 0; i < particles.size(); ++i ) {
        const float r = radiusAcc( particles[i] );
        if( r > maxRadius ) {
            maxRadius = r;
        }
    }

    const float maxDistance = effectRadiusScale * maxRadius;

    if( maxDistance <= 0 ) {
        return;
    }

    // temporary storage for new positions
    std::vector<frantic::graphics::vector3f> newPosition( particles.size() );

    detail::smooth_particle_positions_impl<KernelType> processor( particles.get_channel_map(), effectRadiusScale,
                                                                  lambda, kernel, newPosition );
    detail::process_particle_pairs_in_range( particles, maxDistance, processor, progressLogger );

    for( std::size_t i = 0; i < particles.size(); ++i ) {
        positionAcc.set( particles[i], newPosition[i] );
    }
}

void smooth_particle_positions( frantic::particles::particle_array& particles, float effectRadiusScale, float lambda,
                                frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger ) {
    typedef detail::diffusion_kernel kernel_t;
    kernel_t kernel;
    smooth_particle_positions<kernel_t>( particles, kernel, effectRadiusScale, lambda, progressLogger );
}

smooth_particle_positions_params::smooth_particle_positions_params(
    frantic::particles::particle_array& particles, float effectRadiusScale, float lambda,
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger )
    : particles( particles )
    , effectRadiusScale( effectRadiusScale )
    , lambda( lambda )
    , progressLogger( progressLogger ) {}

void smooth_particle_positions( smooth_particle_positions_params& params ) {
    smooth_particle_positions( params.particles, params.effectRadiusScale, params.lambda, params.progressLogger );
}

calculate_volume_with_anisotropic_kernel_params::calculate_volume_with_anisotropic_kernel_params(
    frantic::particles::particle_array& particles, const frantic::tstring& volumeChannelName,
    frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger )
    : particles( particles )
    , volumeChannelName( volumeChannelName )
    , progressLogger( progressLogger ) {}

void calculate_volume_with_anisotropic_kernel( calculate_volume_with_anisotropic_kernel_params& params ) {
    calculate_volume_with_anisotropic_kernel( params.particles, params.volumeChannelName, params.progressLogger );
}

calculate_anisotropy_params::calculate_anisotropy_params(
    frantic::particles::particle_array& particles, float effectRadiusScale, float anisotropyWindowRadiusScale, float kr,
    std::size_t ne, frantic::volumetrics::implicitsurface::shared_progress_logger_proxy& progressLogger )
    : particles( particles )
    , effectRadiusScale( effectRadiusScale )
    , anisotropyWindowRadiusScale( anisotropyWindowRadiusScale )
    , kr( kr )
    , ne( ne )
    , progressLogger( progressLogger ) {}

void calculate_anisotropy( calculate_anisotropy_params& params ) {
    calculate_anisotropy( params.particles, params.effectRadiusScale, params.anisotropyWindowRadiusScale, params.kr,
                          params.ne, params.progressLogger );
}

} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic
