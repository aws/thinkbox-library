// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/logging/global_progress_logger.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/math/polynomial_roots.hpp>
#include <frantic/math/utils.hpp>
#include <frantic/volumetrics/levelset/level_set_fast_marching_reinitialization.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

using namespace std;
using namespace boost;
using namespace frantic::math::polynomial_roots;
using frantic::graphics::boundbox3;
using frantic::graphics::vector3;
using frantic::graphics::vector3f;
using frantic::math::square;
using frantic::volumetrics::levelset::trilerp_float;

namespace frantic {
namespace volumetrics {
namespace levelset {

namespace detail {

/**
 * This function gets called once for each voxel inside the initial band, to reset its
 * signed distance to the nearest point on the isosurface. The return value represents
 * whether any neighbours of the voxel are on the other side of the interface.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  voxelLength  The length of one voxel.
 * @param  currentIndex  The index of the voxel in question.
 * @param  outSignedDistance  The distance from the given voxel to the closest point on the interface.
 */
bool fast_marching_signed_distance_to_closest_interface( const ris_adjacency& adj, const float* signedDistanceChannel,
                                                         const float voxelLength, const size_t currentIndex,
                                                         float& outSignedDistance ) {

    // find the neighbour indices
    int xposDataIndex = adj[currentIndex].x_pos;
    int xnegDataIndex = adj[currentIndex].x_neg;
    int yposDataIndex = adj[currentIndex].y_pos;
    int ynegDataIndex = adj[currentIndex].y_neg;
    int zposDataIndex = adj[currentIndex].z_pos;
    int znegDataIndex = adj[currentIndex].z_neg;

    // current distance value
    float currentPhi = signedDistanceChannel[currentIndex];

    // default the neighbour phis and signed_distances to infinity
    float xSignedDistanceToInterface = std::numeric_limits<float>::infinity();
    float ySignedDistanceToInterface = std::numeric_limits<float>::infinity();
    float zSignedDistanceToInterface = std::numeric_limits<float>::infinity();

    // Split this test into two cases, depending on the sign of Phi.  The code in each is almost identical,
    // just the signs of all the distance comparisons with 0 are switched.
    if( currentPhi >= 0 ) {
        // look at both neighbours in the x dimension to see if they are on the other side of the interface.
        // If one is, then we store the phi value. If both are, we store the phi value of the one that is closest
        // to the interface.

        if( xposDataIndex >= 0 && signedDistanceChannel[xposDataIndex] <= 0 ) {
            if( xnegDataIndex >= 0 && signedDistanceChannel[xnegDataIndex] <= 0 ) {
                // In this case, both directions along X cross the interface
                float maxSignedDist = max( signedDistanceChannel[xposDataIndex], signedDistanceChannel[xnegDataIndex] );
                xSignedDistanceToInterface = currentPhi / ( currentPhi - maxSignedDist );
            } else {
                // In this case, just the positive direction along X crosses the interface
                xSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[xposDataIndex] );
            }
        } else if( xnegDataIndex >= 0 && signedDistanceChannel[xnegDataIndex] <= 0 ) {
            // In this case, just the negative direction along X crosses the interface
            xSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[xnegDataIndex] );
        }

        // look at both neighbours in the y dimension to see if they are on the other side of the interface.
        // If one is, then we store the phi value. If both are, we store the phi value of the one that is closest
        // to the interface.

        if( yposDataIndex >= 0 && signedDistanceChannel[yposDataIndex] <= 0 ) {
            if( ynegDataIndex >= 0 && signedDistanceChannel[ynegDataIndex] <= 0 ) {
                // In this case, both directions along Y cross the interface
                float maxSignedDist = max( signedDistanceChannel[yposDataIndex], signedDistanceChannel[ynegDataIndex] );
                ySignedDistanceToInterface = currentPhi / ( currentPhi - maxSignedDist );
            } else {
                // In this case, just the positive direction along Y crosses the interface
                ySignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[yposDataIndex] );
            }
        } else if( ynegDataIndex >= 0 && signedDistanceChannel[ynegDataIndex] <= 0 ) {
            // In this case, just the negative direction along Y crosses the interface
            ySignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[ynegDataIndex] );
        }

        // look at both neighbours in the z dimension to see if they are on the other side of the interface.
        // If one is, then we store the phi value. If both are, we store the phi value of the one that is closest
        // to the interface.

        if( zposDataIndex >= 0 && signedDistanceChannel[zposDataIndex] <= 0 ) {
            if( znegDataIndex >= 0 && signedDistanceChannel[znegDataIndex] <= 0 ) {
                // In this case, both directions along Z cross the interface
                float maxSignedDist = max( signedDistanceChannel[zposDataIndex], signedDistanceChannel[znegDataIndex] );
                zSignedDistanceToInterface = currentPhi / ( currentPhi - maxSignedDist );
            } else {
                // In this case, just the positive direction along Z crosses the interface
                zSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[zposDataIndex] );
            }
        } else if( znegDataIndex >= 0 && signedDistanceChannel[znegDataIndex] <= 0 ) {
            // In this case, just the negative direction along Z crosses the interface
            zSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[znegDataIndex] );
        }

    }
    // Otherwise, currentPhi < 0
    else {
        // look at both neighbours in the X dimension to see if they are on the other side of the interface.
        // If one is, then we store the phi value. If both are, we store the phi value of the one that is closest
        // to the interface.

        if( xposDataIndex >= 0 && signedDistanceChannel[xposDataIndex] >= 0 ) {
            if( xnegDataIndex >= 0 && signedDistanceChannel[xnegDataIndex] >= 0 ) {
                // In this case, both directions along X cross the interface
                float minSignedDist = min( signedDistanceChannel[xposDataIndex], signedDistanceChannel[xnegDataIndex] );
                xSignedDistanceToInterface = currentPhi / ( currentPhi - minSignedDist );
            } else {
                // In this case, just the positive direction along X crosses the interface
                xSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[xposDataIndex] );
            }
        } else if( xnegDataIndex >= 0 && signedDistanceChannel[xnegDataIndex] >= 0 ) {
            // In this case, just the negative direction along X crosses the interface
            xSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[xnegDataIndex] );
        }

        // look at both neighbours in the Y dimension to see if they are on the other side of the interface.
        // If one is, then we store the phi value. If both are, we store the phi value of the one that is closest
        // to the interface.

        if( yposDataIndex >= 0 && signedDistanceChannel[yposDataIndex] >= 0 ) {
            if( ynegDataIndex >= 0 && signedDistanceChannel[ynegDataIndex] >= 0 ) {
                // In this case, both directions along Y cross the interface
                float minSignedDist = min( signedDistanceChannel[yposDataIndex], signedDistanceChannel[ynegDataIndex] );
                ySignedDistanceToInterface = currentPhi / ( currentPhi - minSignedDist );
            } else {
                // In this case, just the positive direction along Y crosses the interface
                ySignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[yposDataIndex] );
            }
        } else if( ynegDataIndex >= 0 && signedDistanceChannel[ynegDataIndex] >= 0 ) {
            // In this case, just the negative direction along Y crosses the interface
            ySignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[ynegDataIndex] );
        }

        // look at both neighbours in the Z dimension to see if they are on the other side of the interface.
        // If one is, then we store the phi value. If both are, we store the phi value of the one that is closest
        // to the interface.

        if( zposDataIndex >= 0 && signedDistanceChannel[zposDataIndex] >= 0 ) {
            if( znegDataIndex >= 0 && signedDistanceChannel[znegDataIndex] >= 0 ) {
                // In this case, both directions along Z cross the interface
                float minSignedDist = min( signedDistanceChannel[zposDataIndex], signedDistanceChannel[znegDataIndex] );
                zSignedDistanceToInterface = currentPhi / ( currentPhi - minSignedDist );
            } else {
                // In this case, just the positive direction along Z crosses the interface
                zSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[zposDataIndex] );
            }
        } else if( znegDataIndex >= 0 && signedDistanceChannel[znegDataIndex] >= 0 ) {
            // In this case, just the negative direction along Z crosses the interface
            zSignedDistanceToInterface = currentPhi / ( currentPhi - signedDistanceChannel[znegDataIndex] );
        }
    }

    // now we have 0, 1, 2 or 3 non-infinite phi values, corresponding to the number of neighbours
    // that are on the other side of the interface

    int caseBits = ( ( xSignedDistanceToInterface == std::numeric_limits<float>::infinity() ) ? 0 : 1 ) +
                   ( ( ySignedDistanceToInterface == std::numeric_limits<float>::infinity() ) ? 0 : 2 ) +
                   ( ( zSignedDistanceToInterface == std::numeric_limits<float>::infinity() ) ? 0 : 4 );

    // if we have 2 or 3 neighbours on the other side of the surface
    if( caseBits == 3 || caseBits > 4 ) {
        // if one of the phi values is less than 1e-6, set the outSignedDistance to the smallest of the phi values
        // Note that these phi values are always between 0 and 1 (ie independent of the voxel length) so comparing them
        // to a fixed epsilon.

        if( std::abs( xSignedDistanceToInterface ) < 1e-6 || std::abs( ySignedDistanceToInterface ) < 1e-6 ||
            std::abs( zSignedDistanceToInterface ) < 1e-6 ) {

            // set to the smallest of the phi values
            switch( caseBits ) {
            case 3: // X and Y
                outSignedDistance = std::abs( xSignedDistanceToInterface ) < std::abs( ySignedDistanceToInterface )
                                        ? xSignedDistanceToInterface
                                        : ySignedDistanceToInterface;
                if( currentPhi < 0 ) {
                    outSignedDistance = -outSignedDistance;
                }
                return true;
            case 5: // X and Z
                outSignedDistance = std::abs( xSignedDistanceToInterface ) < std::abs( zSignedDistanceToInterface )
                                        ? xSignedDistanceToInterface
                                        : zSignedDistanceToInterface;
                if( currentPhi < 0 ) {
                    outSignedDistance = -outSignedDistance;
                }
                return true;
            case 6: // Y and Z
                outSignedDistance = std::abs( ySignedDistanceToInterface ) < std::abs( zSignedDistanceToInterface )
                                        ? ySignedDistanceToInterface
                                        : zSignedDistanceToInterface;
                if( currentPhi < 0 ) {
                    outSignedDistance = -outSignedDistance;
                }
                return true;
            case 7: // X, Y and Z
                if( std::abs( xSignedDistanceToInterface ) < std::abs( ySignedDistanceToInterface ) ) {
                    if( std::abs( xSignedDistanceToInterface ) < std::abs( zSignedDistanceToInterface ) ) {
                        outSignedDistance = xSignedDistanceToInterface;
                        if( outSignedDistance > 0 && currentPhi < 0 ) {
                            outSignedDistance = -outSignedDistance;
                        }
                        return true;
                    } else {
                        outSignedDistance = zSignedDistanceToInterface;
                        if( outSignedDistance > 0 && currentPhi < 0 ) {
                            outSignedDistance = -outSignedDistance;
                        }
                        return true;
                    }
                } else {
                    if( std::abs( ySignedDistanceToInterface ) < std::abs( zSignedDistanceToInterface ) ) {
                        outSignedDistance = ySignedDistanceToInterface;
                        if( outSignedDistance > 0 && currentPhi < 0 ) {
                            outSignedDistance = -outSignedDistance;
                        }
                        return true;
                    } else {
                        outSignedDistance = zSignedDistanceToInterface;
                        if( outSignedDistance > 0 && currentPhi < 0 ) {
                            outSignedDistance = -outSignedDistance;
                        }
                        return true;
                    }
                }
            default:
                return false;
            }
        }
    }

    switch( caseBits ) {
    // if we have 1 neighbor on the other side of the interface, return true and just return
    // the distance we computed.
    case 1: // only X
        outSignedDistance =
            currentPhi > 0 ? voxelLength * xSignedDistanceToInterface : -voxelLength * xSignedDistanceToInterface;
        return true;
    case 2: // only Y
        outSignedDistance =
            currentPhi > 0 ? voxelLength * ySignedDistanceToInterface : -voxelLength * ySignedDistanceToInterface;
        return true;
    case 4: // only Z
        outSignedDistance =
            currentPhi > 0 ? voxelLength * zSignedDistanceToInterface : -voxelLength * zSignedDistanceToInterface;
        return true;

    // if we have 2 neighbours on the other side of the interface, return true and compute outSignedDistance
    // the derivation of this equation is in my (Dave) Labratory Notebook on pages 13 and 14.
    case 3: // X and Y
        outSignedDistance = voxelLength * xSignedDistanceToInterface * ySignedDistanceToInterface /
                            sqrtf( xSignedDistanceToInterface * xSignedDistanceToInterface +
                                   ySignedDistanceToInterface * ySignedDistanceToInterface );
        if( currentPhi < 0 ) {
            outSignedDistance = -outSignedDistance;
        }
        return true;
    case 5: // X and Z
        outSignedDistance = voxelLength * xSignedDistanceToInterface * zSignedDistanceToInterface /
                            sqrtf( xSignedDistanceToInterface * xSignedDistanceToInterface +
                                   zSignedDistanceToInterface * zSignedDistanceToInterface );
        if( currentPhi < 0 ) {
            outSignedDistance = -outSignedDistance;
        }
        return true;
    case 6: // Y and Z
        outSignedDistance = voxelLength * ySignedDistanceToInterface * zSignedDistanceToInterface /
                            sqrtf( ySignedDistanceToInterface * ySignedDistanceToInterface +
                                   zSignedDistanceToInterface * zSignedDistanceToInterface );
        if( currentPhi < 0 ) {
            outSignedDistance = -outSignedDistance;
        }
        return true;

    // otherwise we have 3 neighbours on the other side of the interface, return true and compute outSignedDistance
    // the derivation of this equation is in my (Dave) Labratory Notebook on page 15.
    case 7: // X, Y and Z
        outSignedDistance = voxelLength * xSignedDistanceToInterface * ySignedDistanceToInterface *
                            zSignedDistanceToInterface /
                            sqrtf( ySignedDistanceToInterface * ySignedDistanceToInterface *
                                       zSignedDistanceToInterface * zSignedDistanceToInterface +
                                   xSignedDistanceToInterface * xSignedDistanceToInterface *
                                       zSignedDistanceToInterface * zSignedDistanceToInterface +
                                   xSignedDistanceToInterface * xSignedDistanceToInterface *
                                       ySignedDistanceToInterface * ySignedDistanceToInterface );
        if( currentPhi < 0 ) {
            outSignedDistance = -outSignedDistance;
        }
        return true;

    default:
        return false;
    }
}

/**
 * @note this detects a zero crossing if both samples are zero, which is different
 *       from the definition in some parts of our code.
 */
bool has_zero_crossing( const float a, const float b ) { return ( a <= 0 && b >= 0 ) || ( a >= 0 && b <= 0 ); }

float get_distance_fraction_from_interface( const float centerDistance, const float adjDistance ) {
    if( centerDistance == 0 ) {
        return 0;
    } else {
        return centerDistance / ( centerDistance - adjDistance );
    }
}

float get_signed_distance_from_unsigned( float distance, float sign ) {
    if( sign <= 0 ) {
        return -distance;
    } else {
        return std::max( std::numeric_limits<float>::min(), distance );
    }
}

/**
 *  Get the signed distance from the center voxel to the interface,
 * based on the signed distance in the 6-connected neighbours.
 * Neighbours with populatedChannel == 2 are ignored.
 *
 *  This is based on frantic::volumetrics::levelset::detail::fast_marching_signed_distance_to_closest_interface()
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  populatedChannel Populated flag for each voxel.  Adjacent voxels with populated == 2 are ignored.
 * @param  voxelLength  The length of one voxel.
 * @param  currentIndex  The index of the voxel in question.
 * @param  outSignedDistance  The distance from the given voxel to the closest point on the interface.
 */
bool get_linear_signed_distance_to_interface( const ris_adjacency& adj, const float* signedDistanceChannel,
                                              const boost::uint8_t* populatedChannel, const float voxelLength,
                                              const size_t centerIndex, float& outSignedDistance ) {
    assert( signedDistanceChannel );
    assert( populatedChannel );

    const ris_adj_entry& rae = adj[centerIndex];
    const float centerPhi = signedDistanceChannel[centerIndex];

    vector3f distanceToInterface( -std::numeric_limits<float>::max() );
    for( int face = 0; face < 6; ++face ) {
        const boost::int32_t adjIndex = rae[face];
        if( adjIndex >= 0 && populatedChannel[adjIndex] != 2 ) {
            const int axis = neighbor_index_axis( face );
            const float adjPhi = signedDistanceChannel[adjIndex];
            if( has_zero_crossing( centerPhi, adjPhi ) ) {
                distanceToInterface[axis] =
                    std::max( distanceToInterface[axis], get_distance_fraction_from_interface( centerPhi, adjPhi ) );
            }
        }
    }

    int finiteDistanceCount = 0;
    float finiteDistance[3];
    float minDistanceToInterface = std::numeric_limits<float>::max();

    for( int axis = 0; axis < 3; ++axis ) {
        if( distanceToInterface[axis] >= 0 ) {
            finiteDistance[finiteDistanceCount] = distanceToInterface[axis];
            minDistanceToInterface = std::min( minDistanceToInterface, distanceToInterface[axis] );
            ++finiteDistanceCount;
        }
    }

    float distance;

    if( finiteDistanceCount > 1 && minDistanceToInterface < 1.e-6 ) {
        // In fast_marching_signed_distance_to_closest_interface,
        // this case returns the distance in voxel units (other
        // cases are scaled by voxelLength, but this case is not.)
        // This seems like it could be a problem with small voxel
        // lengths, so here I return the distance in world units
        // instead (multiplied by voxelLength).  Time will tell
        // if this is a problem...
        distance = voxelLength * minDistanceToInterface;
        outSignedDistance = get_signed_distance_from_unsigned( distance, centerPhi );
        return true;
    }

    float finiteDistanceSquared[3];

    switch( finiteDistanceCount ) {
    case 1:
        distance = voxelLength * finiteDistance[0];
        outSignedDistance = get_signed_distance_from_unsigned( distance, centerPhi );
        if( !frantic::math::is_finite( outSignedDistance ) )
            throw std::runtime_error(
                "get_linear_signed_distance_to_interface Error: outSignedDistance is not finite (1)" );
        return true;
    case 2:
        distance = voxelLength * finiteDistance[0] * finiteDistance[1] /
                   sqrtf( square( finiteDistance[0] ) + square( finiteDistance[1] ) );
        outSignedDistance = get_signed_distance_from_unsigned( distance, centerPhi );
        if( !frantic::math::is_finite( outSignedDistance ) )
            throw std::runtime_error(
                "get_linear_signed_distance_to_interface Error: outSignedDistance is not finite (2)" );
        return true;
    case 3:
        for( int i = 0; i < finiteDistanceCount; ++i ) {
            finiteDistanceSquared[i] = square( finiteDistance[i] );
        }
        distance = voxelLength * finiteDistance[0] * finiteDistance[1] * finiteDistance[2] /
                   sqrtf( finiteDistanceSquared[0] * finiteDistanceSquared[1] +
                          finiteDistanceSquared[0] * finiteDistanceSquared[2] +
                          finiteDistanceSquared[1] * finiteDistanceSquared[2] );
        outSignedDistance = get_signed_distance_from_unsigned( distance, centerPhi );
        if( !frantic::math::is_finite( outSignedDistance ) ) {
            throw std::runtime_error(
                "get_linear_signed_distance_to_interface Error: outSignedDistance is not finite (3)" );
        }
        return true;
    }

    return false;
}

/**
 *  A class to help indexing in dense 3d arrays.
 */
template <std::size_t xsize, std::size_t ysize, std::size_t zsize>
struct box_indexer {
  public:
    boost::int32_t get_index( const std::size_t x, const std::size_t y, const std::size_t z ) const {
        if( x >= xsize || y >= ysize || z >= zsize )
            return -1;

        return static_cast<boost::int32_t>( x + y * xsize + z * xsize * ysize );
    }

    boost::int32_t get_index( const vector3& v ) const { return get_index( v.x, v.y, v.z ); }

    vector3 get_coord( const std::size_t index ) const {
        if( index >= xsize * ysize * zsize )
            throw std::runtime_error( "box_indexer::get_coord Error: the index " +
                                      boost::lexical_cast<std::string>( index ) + " is greater than the box capacity " +
                                      boost::lexical_cast<std::string>( xsize * ysize * zsize - 1 ) + "." );
        return vector3( static_cast<boost::int32_t>( index % xsize ),
                        static_cast<boost::int32_t>( ( index / xsize ) % ysize ),
                        static_cast<boost::int32_t>( index / ( xsize * ysize ) ) );
    }

    boost::int32_t get_adj_index( const size_t x, const size_t y, const size_t z, const int face ) const {
        return get_index( get_adj_coord( vector3( static_cast<boost::int32_t>( x ), static_cast<boost::int32_t>( y ),
                                                  static_cast<boost::int32_t>( z ) ),
                                         face ) );
    }

    boost::int32_t get_adj_index( const vector3& coord, const int face ) const {
        return get_index( get_adj_coord( coord ), face );
    }

    boost::int32_t get_adj_index( const std::size_t index, const int face ) const {
        return get_index( get_adj_coord( get_coord( index ), face ) );
    }

    vector3 get_adj_coord( const vector3& coord, const int face ) const {
        switch( face ) {
        case rae_index_x_pos:
            return coord + vector3( 1, 0, 0 );
        case rae_index_x_neg:
            return coord + vector3( -1, 0, 0 );
        case rae_index_y_pos:
            return coord + vector3( 0, 1, 0 );
        case rae_index_y_neg:
            return coord + vector3( 0, -1, 0 );
        case rae_index_z_pos:
            return coord + vector3( 0, 0, 1 );
        case rae_index_z_neg:
            return coord + vector3( 0, 0, -1 );
        default:
            throw std::runtime_error( "box_indexer::get_adj_coord() Error: invalid face number " +
                                      boost::lexical_cast<std::string>( face ) + "." );
        }
    }

    std::size_t get_xsize( void ) const { return xsize; }
    std::size_t get_ysize( void ) const { return ysize; }
    std::size_t get_zsize( void ) const { return zsize; }

    std::size_t size( void ) const { return xsize * ysize * zsize; }
};

/**
 *  Extrapolate phi from unpopulated voxels to populated voxels,
 * while matching the gradient in grad.
 *
 * This uses a brute-force algorithm -- I assume it won't be
 * called too often.
 */
template <class box_indexer_t>
bool extrapolate_matching_gradient( float* phi, boost::uint8_t* pop, vector3f* grad, const box_indexer_t& indexer ) {
    const size_t size = indexer.size();

    size_t popCount = 0;
    for( size_t i = 0; i < size; ++i ) {
        if( pop[i] )
            popCount++;
    }

    if( popCount == 0 )
        return false;

    size_t unpopCount = size - popCount;

    boost::uint8_t popInput = 1;
    while( unpopCount > 0 ) {
        unpopCount = 0;
        for( size_t z = 0; z < indexer.get_zsize(); ++z ) {
            for( size_t y = 0; y < indexer.get_ysize(); ++y ) {
                for( size_t x = 0; x < indexer.get_xsize(); ++x ) {
                    const boost::int32_t i = indexer.get_index( x, y, z );
                    if( pop[i] == 0 ) {
                        std::size_t inputCount = 0;
                        float sum = 0;
                        for( int face = 0; face < 6; ++face ) {
                            const boost::int32_t adjIndex = indexer.get_adj_index( x, y, z, face );
                            if( adjIndex >= 0 ) {
                                if( pop[adjIndex] > 0 ) {
                                    if( pop[adjIndex] == popInput ) {
                                        ++inputCount;
                                        if( is_neighbor_index_direction_positive( face ) ) {
                                            sum +=
                                                ( 1.f - grad[adjIndex][neighbor_index_axis( face )] ) * phi[adjIndex];
                                        } else {
                                            sum += ( 1.f + grad[i][neighbor_index_axis( face )] ) * phi[adjIndex];
                                        }
                                    } else if( pop[adjIndex] < popInput ) {
                                        throw std::runtime_error(
                                            "extrapolate_matching_gradient Error: adjacent pop is not popInput" );
                                    }
                                }
                            }
                        }
                        if( inputCount > 0 ) {
                            phi[i] = sum / static_cast<float>( inputCount );
                            pop[i] = 1 + popInput;
                        } else {
                            ++unpopCount;
                        }
                    }
                }
            }
        }
        ++popInput;
    }
    return true;
}

/**
 *  Extrapolate grad from populated to unpopulated voxels.
 *
 * This uses a brute-force algorithm -- I assume it won't be
 * called too often.
 *
 */
template <class box_indexer_t>
void extrapolate_staggered( vector3f* grad, vector3* pop, const int axis, const box_indexer_t& indexer ) {
    const size_t size = indexer.size();
    if( size == 0 )
        return;

    size_t popCount = 0;
    size_t sampleCount = 0;
    for( size_t z = 0; z < indexer.get_zsize(); ++z ) {
        for( size_t y = 0; y < indexer.get_ysize(); ++y ) {
            for( size_t x = 0; x < indexer.get_xsize(); ++x ) {
                const vector3 coord( static_cast<boost::int32_t>( x ), static_cast<boost::int32_t>( y ),
                                     static_cast<boost::int32_t>( z ) );
                if( coord[axis] != 0 ) {
                    ++sampleCount;
                    const boost::int32_t i = indexer.get_index( coord );
                    if( pop[i][axis] )
                        ++popCount;
                }
            }
        }
    }

    if( popCount == 0 )
        return;

    size_t unpopCount = sampleCount - popCount;

    // Samples with this value in the pop vector will be used as
    // inputs.
    boost::int32_t popInput = 1;
    while( unpopCount > 0 ) {
        unpopCount = 0;

        for( size_t z = 0; z < indexer.get_zsize(); ++z ) {
            for( size_t y = 0; y < indexer.get_ysize(); ++y ) {
                for( size_t x = 0; x < indexer.get_xsize(); ++x ) {
                    vector3 coord( static_cast<boost::int32_t>( x ), static_cast<boost::int32_t>( y ),
                                   static_cast<boost::int32_t>( z ) );
                    if( coord[axis] == 0 )
                        continue;

                    const boost::int32_t i = indexer.get_index( coord );
                    if( pop[i][axis] == 0 ) {
                        // Look for adjacent inputs, and sum them.
                        size_t inputCount = 0;
                        float sum = 0;

                        for( int face = 0; face < 6; ++face ) {
                            vector3 adjCoord = indexer.get_adj_coord( coord, face );
                            if( adjCoord[axis] == 0 ) {
                                continue;
                            }
                            boost::int32_t adjIndex = indexer.get_index( adjCoord );
                            if( adjIndex >= 0 ) {
                                if( pop[adjIndex][axis] > 0 ) {
                                    if( pop[adjIndex][axis] == popInput ) {
                                        ++inputCount;
                                        sum += grad[adjIndex][axis];
                                    } else if( pop[adjIndex][axis] < popInput ) {
                                        throw std::runtime_error( "extrapolate_staggered() Error: voxel has input that "
                                                                  "is not current input" );
                                    }
                                }
                            }
                        }

                        if( inputCount > 0 ) {
                            pop[i][axis] = 1 + popInput;
                            grad[i][axis] = sum / static_cast<float>( inputCount );
                        } else {
                            ++unpopCount;
                        }
                    }
                }
            }
        }
        ++popInput;
    }
}

/**
 *  Perform trilinear interpolation.
 *
 * copied from frantic::graphics::trilerp_value
 * data must be an array of 8 elements of type data_type.
 * they are indexed in the order data[0] (0,0,0), data[1] (1,0,0), data[2] (0,1,0), data[3] (1,1,0), etc.
 * coord must have elements in the range 0..1.  these are the trilerp weights
 */
template <class data_type>
void trilerp( data_type* data, const vector3f& coord, data_type& out ) {
    // if( coord.x < 0 || coord.y < 0 || coord.z < 0 || coord.x > 1.f || coord.y > 1.f || coord.z > 1.f )
    // throw std::runtime_error( "trilerp Error: all components of coord must be in the range [0,1], but coord is " +
    // coord.str() + "." );

    const float dx = coord.x, dy = coord.y, dz = coord.z;
    const float dxAlt = 1.f - dx, dyAlt = 1.f - dy, dzAlt = 1.f - dz;

    out = dzAlt * ( dyAlt * ( dxAlt * data[0] + dx * data[1] ) + dy * ( dxAlt * data[2] + dx * data[3] ) ) +
          dz * ( dyAlt * ( dxAlt * data[4] + dx * data[5] ) + dy * ( dxAlt * data[6] + dx * data[7] ) );
}

/**
 *  Return true if the voxel is on the interface between inside
 * and outside voxels.
 *
 * @note tag_interface_voxels counts adjacent voxels with phi == 0
 *		as interface voxels, so I do the same here.
 *
 * @param adj adjacency information which includes the specified voxel
 *		index.
 * @param signedDistanceChannel an array of signed distance values.
 *		Signed distance <= 0 is inside, while > 0 is outside.
 * @param index the voxel to examine.
 * @return true if the index voxel has a neighboring voxel on the
 *		other side of the interface.
 */
bool is_interface_voxel( const ris_adjacency& adj, const float* signedDistanceChannel, const size_t index ) {
    const float centerPhi = signedDistanceChannel[index];

    if( centerPhi <= 0 ) {
        for( int face = 0; face < 6; ++face ) {
            const boost::int32_t adjIndex = adj[index][face];
            if( adjIndex >= 0 && signedDistanceChannel[adjIndex] >= 0 ) {
                return true;
            }
        }
    } else {
        for( int face = 0; face < 6; ++face ) {
            const boost::int32_t adjIndex = adj[index][face];
            if( adjIndex >= 0 && signedDistanceChannel[adjIndex] <= 0 ) {
                return true;
            }
        }
    }

    return false;
}

/**
 *  Return true if the voxel is on the interface between inside
 * and outside voxels.  Adjacent voxels with populated == 2 are
 * ignored.
 *
 * @note tag_interface_voxels counts adjacent voxels with phi == 0
 *		as interface voxels, so I do the same here.
 *
 * @param adj adjacency information which includes the specified voxel
 *		index.
 * @param signedDistanceChannel an array of signed distance values.
 *		Signed distance <= 0 is inside, while > 0 is outside.
 * @param populatedChannel a populated channel.  If a neighbouring
 *		voxel has populatedChannel == 2, then its distance will be
 *		ignored.
 * @param index the voxel to examine.
 * @return true if the index voxel has a neighboring voxel on the
 *		other side of the interface.
 */
bool is_interface_voxel( const ris_adjacency& adj, const float* signedDistanceChannel,
                         const boost::uint8_t* populatedChannel, const size_t index ) {
    assert( populatedChannel );

    const float centerPhi = signedDistanceChannel[index];

    if( centerPhi <= 0 ) {
        for( int face = 0; face < 6; ++face ) {
            const boost::int32_t adjIndex = adj[index][face];
            if( adjIndex >= 0 && populatedChannel[adjIndex] != 2 && signedDistanceChannel[adjIndex] >= 0 ) {
                return true;
            }
        }
    } else {
        for( int face = 0; face < 6; ++face ) {
            const boost::int32_t adjIndex = adj[index][face];
            if( adjIndex >= 0 && populatedChannel[adjIndex] != 2 && signedDistanceChannel[adjIndex] <= 0 ) {
                return true;
            }
        }
    }

    return false;
}

class trilerp_cache {
    // the cached query coordinate
    vector3f m_voxelCoordinate;

    // the bottom-left of the current box
    vector3 m_minimumCoord;

    // Can the boxes be filled for the current m_minimumCoord ?
    bool m_isPopulated;

    const rle_index_spec* m_ris;
    const float* m_signedDistance;
    const boost::uint8_t* m_populatedChannel;

    static const std::size_t m_indexCount = 27;
    float m_phiBox[m_indexCount];
    vector3f m_gradBox[m_indexCount];

    box_indexer<3, 3, 3> m_indexer;

    /**
     *  Attempt to fill m_phiBox and m_gradBox, and set m_isPopulated,
     * m_minimumCoord, and m_voxelCoord.  If the return value is true,
     * then these member variables have been set.  Otherwise these
     * member variables are undefined.
     *
     * The voxel structure is copied from
     * frantic::volumetrics::levelset::trilerp_staggered_centered_signed_distance_gradient() .
     */
    bool try_fill_box( const vector3f& voxelCoordinate ) {
        if( voxelCoordinate == m_voxelCoordinate ) {
            return m_isPopulated;
        }
        m_voxelCoordinate = voxelCoordinate;

        // We want the two closest staggered voxels, which requires the nearest 3 voxels to compute.
        // Ex. The x's are voxel centers, the dashes integer coords. We can see for both sample
        //      points(dots) that the same three voxels are closest in the x direction.
        //  0       1       2       3 <-- Staggered coordinates
        //      0       1       2     <-- Non-staggered coordinates
        //  |   x   |   x . |   x   | <-- Sample is lerped from staggered voxels 1,2
        //  |   x   | . x   |   x   | <-- Sample is lerped from staggered voxels 1,2
        //  if voxelCoordinate.x == 0.0 then xstart == -1
        //  if voxelCoordinate.x == 0.5 then xstart == -1
        //  if voxelCoordinate.x == 1.0 then xstart ==  0
        const float xStaggeredFloor = floorf( voxelCoordinate.x );
        const float yStaggeredFloor = floorf( voxelCoordinate.y );
        const float zStaggeredFloor = floorf( voxelCoordinate.z );

        // We want the 3x3x3 voxels bracketing the middle one, so subtract 1.
        const int xstart = static_cast<int>( xStaggeredFloor ) - 1;
        const int ystart = static_cast<int>( yStaggeredFloor ) - 1;
        const int zstart = static_cast<int>( zStaggeredFloor ) - 1;

        const vector3 minimum( xstart, ystart, zstart );
        // have we already filled this box ?
        if( minimum == m_minimumCoord ) {
            return m_isPopulated;
        }
        m_minimumCoord = minimum;
        m_isPopulated = false;
        vector3 maximum( xstart + 2, ystart + 2, zstart + 2 );

        boost::int32_t dataIndices[m_indexCount];
        m_ris->fill_data_index_box( boundbox3( m_minimumCoord, maximum ), dataIndices );

        // If a voxel's populatedChannel == 2, then set its dataIndex
        // to -1 so that its value will be ignored.
        if( m_populatedChannel ) {
            for( int i = 0; i < m_indexCount; ++i ) {
                const boost::int32_t dataIndex = dataIndices[i];
                if( dataIndex >= 0 ) {
                    if( m_populatedChannel[dataIndex] == 2 ) {
                        dataIndices[dataIndex] = -1;
                    }
                }
            }
        }

        int definedCount = 0;
        boost::uint8_t phiPop[m_indexCount];
        for( int i = 0; i < m_indexCount; ++i ) {
            const boost::int32_t dataIndex = dataIndices[i];
            if( dataIndex >= 0 ) {
                m_phiBox[i] = m_signedDistance[dataIndex];
                phiPop[i] = 1;
                ++definedCount;
            } else {
                phiPop[i] = 0;
            }
        }

        if( definedCount == m_indexCount ) {
            // If all of the voxels are already defined,
            // then we just need to calculate the gradient.
            // reinitializationProfilingGroup.psFillGrad.enter();
            size_t index = 0;
            for( size_t z = 0; z < 3; ++z ) {
                for( size_t y = 0; y < 3; ++y ) {
                    for( size_t x = 0; x < 3; ++x ) {
                        const size_t adjIndices[] = { index - 1, index - 3, index - 9 };
                        const size_t coord[] = { x, y, z };
                        for( int axis = 0; axis < 3; ++axis ) {
                            if( coord[axis] == 0 )
                                continue;
                            const size_t adjIndex = adjIndices[axis];

                            m_gradBox[index][axis] = m_phiBox[index] - m_phiBox[adjIndex];
                        }
                        ++index;
                    }
                }
            }
            // reinitializationProfilingGroup.psFillGrad.exit();
        } else {
            // If there are undefined voxels, then
            // we begin by calculating the gradient between
            // defined voxels.
            vector3 gradPop[m_indexCount];

            for( int i = 0; i < m_indexCount; ++i ) {
                gradPop[i].set( 0 );
            }

            { // scope for index loop
                size_t index = 0;
                for( size_t z = 0; z < 3; ++z ) {
                    for( size_t y = 0; y < 3; ++y ) {
                        for( size_t x = 0; x < 3; ++x ) {
                            const size_t adjIndices[] = { index - 1, index - 3, index - 9 };
                            const size_t coord[] = { x, y, z };
                            if( phiPop[index] ) {
                                for( int axis = 0; axis < 3; ++axis ) {
                                    if( coord[axis] == 0 )
                                        continue;
                                    const size_t adjIndex = adjIndices[axis];

                                    if( phiPop[adjIndex] ) {
                                        m_gradBox[index][axis] = m_phiBox[index] - m_phiBox[adjIndex];
                                        gradPop[index][axis] = 1;
                                    } else {
                                        m_gradBox[index][axis] = 0;
                                    }
                                }
                            }
                            ++index;
                        }
                    }
                }
            }

            // Now we extrapolate grad and phi from the
            // defined voxels into the undefined voxels.
            for( int axis = 0; axis < 3; ++axis )
                detail::extrapolate_staggered( m_gradBox, gradPop, axis, m_indexer );
            const bool success = detail::extrapolate_matching_gradient( m_phiBox, phiPop, m_gradBox, m_indexer );
            if( !success ) {
                return false;
            }
        }

        m_isPopulated = true;
        return true;
    }

    trilerp_cache();

  public:
    /**
     * @param populatedChannel a populated flag channel.  Voxels with
     *		populated flag == 2 will be ignored.  Can be NULL, in which
     *		case all voxels will be used.
     */
    trilerp_cache( const rle_index_spec* ris, const float* signedDistance, const boost::uint8_t* populatedChannel = 0 )
        : m_ris( ris )
        , m_signedDistance( signedDistance )
        , m_populatedChannel( populatedChannel )
        , m_isPopulated( false )
        , m_minimumCoord( std::numeric_limits<boost::int32_t>::max() )
        , m_voxelCoordinate( std::numeric_limits<float>::max() ) {}

    /**
     *  Estimate the signed distance at the specified voxel
     * coordinate.
     *
     * @param voxelCoordinate the signed distance is estimated at
     *		this point.
     * @param[out] outPhi if the function returns true, then this
     *		is the estimated signed distance.
     * @return true if the signed distance could be estimated, or
     *		false otherwise.
     */
    bool get_signed_distance( const vector3f& voxelCoordinate, float& outPhi ) {
        const bool success = try_fill_box( voxelCoordinate );
        if( !success )
            return false;

        const vector3f sampleCoord( voxelCoordinate - vector3f( 0.5f ) );
        const vector3f regFloor( floorf( sampleCoord.x ), floorf( sampleCoord.y ), floorf( sampleCoord.z ) );
        const vector3f regAlpha = sampleCoord - regFloor;
        const vector3f regAlphaAlt = vector3f( 1.f ) - regAlpha;
        const vector3 minCorner( int( regFloor.x ), int( regFloor.y ), int( regFloor.z ) );
        const vector3 localMinCorner( minCorner - m_minimumCoord );
        for( int axis = 0; axis < 3; ++axis ) {
            if( localMinCorner[axis] < 0 || localMinCorner[axis] > 1 )
                throw std::runtime_error( "trilerp_cache::get_signed_distance() Error: offset " + localMinCorner.str() +
                                          " is out of range for voxel coordinate " + voxelCoordinate.str() );
        }

        const std::size_t baseIndex = localMinCorner.x + 3 * localMinCorner.y + 9 * localMinCorner.z;

        outPhi =
            regAlphaAlt.z *
                ( regAlphaAlt.y * ( regAlphaAlt.x * m_phiBox[baseIndex + 0] + regAlpha.x * m_phiBox[baseIndex + 1] ) +
                  regAlpha.y * ( regAlphaAlt.x * m_phiBox[baseIndex + 3] + regAlpha.x * m_phiBox[baseIndex + 4] ) ) +
            regAlpha.z *
                ( regAlphaAlt.y * ( regAlphaAlt.x * m_phiBox[baseIndex + 9] + regAlpha.x * m_phiBox[baseIndex + 10] ) +
                  regAlpha.y * ( regAlphaAlt.x * m_phiBox[baseIndex + 12] + regAlpha.x * m_phiBox[baseIndex + 13] ) );

        return true;
    }

    /**
     *  Estimate the signed distance gradient at the specified voxel
     * coordinate.
     *
     * This is copied from
     * frantic::volumetrics::levelset::trilerp_staggered_centered_signed_distance_gradient() .
     *
     * @param voxelCoordinate the signed distance gradient is
     *		estimated at this point.
     * @param[out] outSignedDistanceGradient if the function returns
     *		true, then this	is the estimated gradient.
     * @return true if the gradient could be estimated, or false
     *		otherwise.
     */
    bool get_gradient( const vector3f& voxelCoordinate, vector3f& outSignedDistanceGradient ) {
        const bool success = try_fill_box( voxelCoordinate );
        if( !success )
            return false;

        // Compute the regular, non-staggered trilerp bounds as well.
        //     0       1       2
        //|   x   |   x . |   x   | <-- Sample is lerped from voxels 1,2
        //|   x   | . x   |   x   | <-- Sample is lerped from voxels 0,1
        float xRegFloor = floorf( voxelCoordinate.x - 0.5f );
        float yRegFloor = floorf( voxelCoordinate.y - 0.5f );
        float zRegFloor = floorf( voxelCoordinate.z - 0.5f );

        // The alpha value is (voxelCoordinate.x - 0.5f) - std::floor(voxelCoordinate.x - 0.5f)
        float xRegAlpha = voxelCoordinate.x - 0.5f - xRegFloor;
        float yRegAlpha = voxelCoordinate.y - 0.5f - yRegFloor;
        float zRegAlpha = voxelCoordinate.z - 0.5f - zRegFloor;

        // minimum is xRegFloor or xRegFloor+1, therefore xOff is 0 or 1.
        int xOff = (int)xRegFloor - m_minimumCoord.x;
        int yOff = (int)yRegFloor - m_minimumCoord.y;
        int zOff = (int)zRegFloor - m_minimumCoord.z;
        if( xOff < 0 || xOff > 1 || yOff < 0 || yOff > 1 || zOff < 0 || zOff > 1 ) {
            throw std::runtime_error(
                "trilerp_cache::get_gradient() Error: offset (" + boost::lexical_cast<std::string>( xOff ) + "," +
                boost::lexical_cast<std::string>( yOff ) + "," + boost::lexical_cast<std::string>( zOff ) +
                ") is out of range for voxel coordinate " + voxelCoordinate.str() );
        }

        float xStaggeredFloor = floorf( voxelCoordinate.x );
        float yStaggeredFloor = floorf( voxelCoordinate.y );
        float zStaggeredFloor = floorf( voxelCoordinate.z );

        float xStaggeredAlpha = voxelCoordinate.x - xStaggeredFloor;
        float yStaggeredAlpha = voxelCoordinate.y - yStaggeredFloor;
        float zStaggeredAlpha = voxelCoordinate.z - zStaggeredFloor;

        // For each component, we grab the three distance samples needed to make the
        // bracketing staggered distance differences.
        //     0       1       2     <-- Centered indices
        //|   x   |   x   |   x   | <-- Centered samples (x's)
        //|   x   d   x   d   x   | <-- Staggered samples
        //
        // Ex.
        //|   x   |   x + |   x   | <-- Sample point is '+'
        //         d-----+-d         <-- Lerp between nearest staggered differences
        //
        // We now trilerp the staggered differences. For the x-axis we compute the differences
        //  at [fl(x), fl(y-0.5)+0.5, fl(z-0.5)+0.5], and [fl(x)+1, fl(y-0.5)+0.5, fl(y-0.5)+0.5], etc.
        // Ex.
        //  The a are staggered differences in the x-axis, the ß are staggered differences in the y-axis.
        //  The + is the desired sample point. Notice that the bracketing samples are different for the
        //  x, and y computations.
        //  ------- ------- -------
        //|       |       |       |
        //|   x   a   x   a   x   |
        //|       |       |       |
        //  ------- ---ß--- ---ß---
        //|       |     + |       |
        //|   x   a   x   a   x   |
        //|       |       |       |
        //  ------- ---ß--- ---ß---
        //|       |       |       |
        //|   x   |   x   |   x   |
        //|       |       |       |
        //  ------- ------- -------
        //
        // Since we are using a 3x3x3 box of indices, we must offset the the array index
        // by one for a step in the x-direction, 3 for a step in the y-direction, and 9
        // for a step in the z-direction.
        // Ex.
        //  to get samples at [0,1,2] and [1,1,2] (relative to the index box) we would
        //  use array elements [0,1,2] = 0 + (1 * 3) + (2 * 9) = 21, and [1,1,2] = 1 + (1*3) + (2*9) = 22

        // For the x component, compute the 8 bracketing x-staggered differences.
        outSignedDistanceGradient.x = 0;
        {
            float d1, d2, w;
            int baseIndex = yOff * 3 + zOff * 9;

            d1 = m_gradBox[baseIndex + 1].x; // [x+1,  y,  z] - [x,  y,  z]
            d2 = m_gradBox[baseIndex + 2].x; // [x+2,  y,  z] - [x+1,  y,  z]
            w = ( 1.f - zRegAlpha ) * ( 1.f - yRegAlpha );
            outSignedDistanceGradient.x += w * ( d1 + xStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 4].x; // [x+1, y+1, z] - [x, y+1, z]
            d2 = m_gradBox[baseIndex + 5].x; // [x+2, y+1, z] - [x+1, y+1, z]
            w = ( 1.f - zRegAlpha ) * yRegAlpha;
            outSignedDistanceGradient.x += w * ( d1 + xStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 10].x; // [x+1, y, z+1] - [x,   y, z+1]
            d2 = m_gradBox[baseIndex + 11].x; // [x+2, y, z+1] - [x+1, y, z+1]
            w = zRegAlpha * ( 1.f - yRegAlpha );
            outSignedDistanceGradient.x += w * ( d1 + xStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 13].x; // [x+1, y+1, z+1] - [x,   y+1, z+1]
            d2 = m_gradBox[baseIndex + 14].x; // [x+2, y+1, z+1] - [x+1, y+1, z+1]
            w = zRegAlpha * yRegAlpha;
            outSignedDistanceGradient.x += w * ( d1 + xStaggeredAlpha * ( d2 - d1 ) );
        }

        // For the y component, compute the 8 bracketing y-staggered differences.
        outSignedDistanceGradient.y = 0;
        {
            float d1, d2, w;
            int baseIndex = xOff + zOff * 9;

            d1 = m_gradBox[baseIndex + 3].y; // [x,  y+1,  z] - [x,  y  ,  z]
            d2 = m_gradBox[baseIndex + 6].y; // [x,  y+2,  z] - [x,  y+1,  z]
            w = ( 1.f - zRegAlpha ) * ( 1.f - xRegAlpha );
            outSignedDistanceGradient.y += w * ( d1 + yStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 4].y; // [x+1,  y+1,  z] - [x+1,  y  ,  z]
            d2 = m_gradBox[baseIndex + 7].y; // [x+1,  y+2,  z] - [x+1,  y+1,  z]
            w = ( 1.f - zRegAlpha ) * xRegAlpha;
            outSignedDistanceGradient.y += w * ( d1 + yStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 12].y; // [x,  y+1,  z+1] - [x,  y  ,  z+1]
            d2 = m_gradBox[baseIndex + 15].y; // [x,  y+2,  z+1] - [x,  y+1,  z+1]
            w = zRegAlpha * ( 1.f - xRegAlpha );
            outSignedDistanceGradient.y += w * ( d1 + yStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 13].y; // [x+1,  y+1,  z+1] - [x+1,  y  ,  z+1]
            d2 = m_gradBox[baseIndex + 16].y; // [x+1,  y+2,  z+1] - [x+1,  y+1,  z+1]
            w = zRegAlpha * xRegAlpha;
            outSignedDistanceGradient.y += w * ( d1 + yStaggeredAlpha * ( d2 - d1 ) );
        }

        // For the z component, compute the 8 bracketing z-staggered differences.
        outSignedDistanceGradient.z = 0;
        {
            float d1, d2, w;
            int baseIndex = xOff + yOff * 3;

            d1 = m_gradBox[baseIndex + 9].z;  // [x,  y,  z+1] - [x,  y,  z]
            d2 = m_gradBox[baseIndex + 18].z; // [x,  y,  z+2] - [x,  y,  z+1]
            w = ( 1.f - yRegAlpha ) * ( 1.f - xRegAlpha );
            outSignedDistanceGradient.z += w * ( d1 + zStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 10].z; // [x+1,  y,  z+1] - [x+1,  y,  z]
            d2 = m_gradBox[baseIndex + 19].z; // [x+1,  y,  z+2] - [x+1,  y,  z+1]
            w = ( 1.f - yRegAlpha ) * xRegAlpha;
            outSignedDistanceGradient.z += w * ( d1 + zStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 12].z; // [x,  y+1,  z+1] - [x,  y+1,  z]
            d2 = m_gradBox[baseIndex + 21].z; // [x,  y+1,  z+2] - [x,  y+1,  z+1]
            w = yRegAlpha * ( 1.f - xRegAlpha );
            outSignedDistanceGradient.z += w * ( d1 + zStaggeredAlpha * ( d2 - d1 ) );

            d1 = m_gradBox[baseIndex + 13].z; // [x+1,  y+1,  z+1] - [x+1,  y+1,  z]
            d2 = m_gradBox[baseIndex + 22].z; // [x+1,  y+1,  z+2] - [x+1,  y+1,  z+1]
            w = yRegAlpha * xRegAlpha;
            outSignedDistanceGradient.z += w * ( d1 + zStaggeredAlpha * ( d2 - d1 ) );
        }

        return true;
    }
};

/**
 *  Find the distance from the specified voxelNumber to the interface.
 *
 * @param ris a filled rle_index_spec.
 * @param signedDistanceChannel an array of signed distance values.
 * @param voxelNumber the index of the voxel for which to find the
 *		distance to the interface.
 * @param voxelCoord the voxel coordinate.
 * @param[out] outSignedDistance If the function returns true, then
 *		this holds the distance from the specified voxel to the
 *		interface.
 * @param populatedChannel a populated channel.  Voxels with
 *		populated == 2 are not interface voxels and will be
 *		ignored.
 * @return true if the specified voxelNumber is on the interface and
 *		the signed distance to the interface was found.
 */
bool fast_marching_signed_distance_to_closest_interface_particle_search(
    const rle_index_spec& ris, const float* signedDistanceChannel, const boost::uint8_t* populatedChannel,
    const float voxelLength, const size_t voxelNumber, const vector3 voxelCoord, float& outSignedDistance ) {
    const ris_adjacency& adj = ris.get_cached_adjacency();
    const bool isInterfaceVoxel = is_interface_voxel( adj, signedDistanceChannel, populatedChannel, voxelNumber );

    if( isInterfaceVoxel ) {
        const float tolerance = 0.0005f * voxelLength;
        const int maxIterCount = 20;

        const float initialPhi = signedDistanceChannel[voxelNumber];
        const float invVoxelLength = 1.f / voxelLength;

        float linearSignedDistance;
        bool sanityCheck = get_linear_signed_distance_to_interface( adj, signedDistanceChannel, populatedChannel,
                                                                    voxelLength, voxelNumber, linearSignedDistance );
        if( !sanityCheck ) {
            throw std::runtime_error( "fast_marching_signed_distance_to_closest_interface_particle_search Error: is "
                                      "interface voxel, but could not measure linear distance to interface." );
        }

        if( fabsf( initialPhi ) < tolerance ) {
            outSignedDistance = linearSignedDistance;
        } else {
            vector3f voxelCenter = vector3f( voxelCoord ) + vector3f( 0.5f );
            vector3f position = voxelCenter;
            vector3f gradient;
            float positionPhi = initialPhi;

            trilerp_cache trilerpCache( &ris, signedDistanceChannel, populatedChannel );

            // Move toward the interface by moving in the gradient direction.
            // The step size is the signed distance at the current position.
            bool success = true;
            for( int iterCount = 0; iterCount < maxIterCount; ++iterCount ) {
                success = trilerpCache.get_gradient( position, gradient );
                if( !success )
                    break;

                gradient *= invVoxelLength;
                gradient.normalize();

                position -= invVoxelLength * positionPhi * gradient;
                if( !position.is_finite() ) {
                    success = false;
                    break;
                }

                success = trilerpCache.get_signed_distance( position, positionPhi );
                if( !success || fabsf( positionPhi ) < tolerance )
                    break;
            }

            const float particleVoxelDistance = ( position - voxelCenter ).get_magnitude();
            const float particleDistance = voxelLength * particleVoxelDistance;

            // Use the linear distance estimate if the particle search
            // encountered a problem.
            if( !success || particleDistance > fabs( linearSignedDistance ) || fabs( positionPhi ) > tolerance ) {
                outSignedDistance = linearSignedDistance;
            } else {
                // For positive phi, make sure we don't round down to zero,
                // which would change the voxel from outside the interface
                // to inside.
                outSignedDistance = initialPhi <= 0 ? -particleDistance
                                                    : std::max( std::numeric_limits<float>::min(), particleDistance );
            }
            if( !frantic::math::is_finite( outSignedDistance ) )
                throw std::runtime_error(
                    "fast_marching_signed_distance_to_closest_interface_particle_search Internal Error: "
                    "particle search produced infinite signed distance" );
        }

        if( ( initialPhi <= 0 && outSignedDistance > 0 ) || ( initialPhi > 0 && outSignedDistance <= 0 ) ) {
            throw std::runtime_error(
                "fast_marching_signed_distance_to_closest_interface_particle_search Internal Error: "
                "changed the sign of an interface voxel (from " +
                boost::lexical_cast<std::string>( initialPhi ) + " to " +
                boost::lexical_cast<std::string>( outSignedDistance ) + ")." );
        }
    }

    return isInterfaceVoxel;
}

/**
 * The 'Populated' channel will have 0 for unpopulated voxels and 1 for populated voxels, and on output
 * it will have 1 for each voxel that was successfully reinitialized. A value of 4 is also used for voxels
 * that are not yet reinitialized but will be when the function terminates (ie in Near). A value of 5 is
 * used to mark the initial band points so that we don't update their distances in the march.
 * There should be no 4 or 5 values in the populated channel on exit.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  voxelLength  The length of one voxel.
 * @param  updatePopulatedInterfaceDistances  Determines whether the intial band should update the distances of points
 * that are populated and on the interface.
 * @param  nearInQ  The inside portion of the initial band in the queue.
 * @param  nearOutQ  The outside portion of the initial band in the queue.
 */
void fast_marching_initial_band( const rle_index_spec& ris, unsigned char* populatedChannel,
                                 float* signedDistanceChannel, const float voxelLength,
                                 const bool updatePopulatedInterfaceDistances, queue<pair<float, int>>& nearInQ,
                                 queue<pair<float, int>>& nearOutQ ) {

    const ris_adjacency& adj = ris.get_cached_adjacency();

    // set up a vector to temporarily store new distance values
    vector<pair<float, int>> newValuesIn;
    vector<pair<float, int>> newValuesOut;

    // It seems reasonable to assume that the initial band will not exceed half of the defined voxels
    // (of course, more memory will be allocated later on if it does)
    newValuesIn.reserve( adj.data_size() / 4 );
    newValuesOut.reserve( adj.data_size() / 4 );

    if( updatePopulatedInterfaceDistances ) {
        // Note: This is split into 2 cases so that the updated distances below will only get allocated when
        // we want to update the distances

        // for each defined voxel
        for( rle_defined_iterator i( ris ), ie( ris, true ); i != ie; ++i ) {
            const size_t iterIndex = i.get_data_index();

            // check that the signedDistance value is not NAN
            if( frantic::math::is_nan( signedDistanceChannel[iterIndex] ) ) {
                throw runtime_error(
                    "fast_marching_initial_band() - The 'SignedDistance' channel contains a NAN at index " +
                    boost::lexical_cast<string>( iterIndex ) + "." );
            }

            // if we have a populated source
            if( populatedChannel[iterIndex] == 1 ) {
                const float initialSignedDistance = signedDistanceChannel[iterIndex];
                float distance = initialSignedDistance;

                // for distances that are very close to the surface, numerical error will eventually set them to 0 (when
                // they are no longer able to represent them as a floating point, ie less than 1e-45) Once it has been
                // established that the distance to the interface is 0, there is no need to update the distance (we do
                // not want to move the interface), so if the distance is exactly zero do not update the distance
                //
                // I am removing this for now.  I agree that if the distance
                // is 0, then we don't need to update it.  However, I think we
                // should still add this voxel to the priority queue so that
                // its neighbours are updated at the right time.
                // if( distance == 0.0f ) {
                //	continue;
                //}

                // if there is at least one sign change between cells, this will update the distance
                if( distance != 0 ) {
                    // const bool onInterface = get_linear_signed_distance_to_interface( adj, signedDistanceChannel,
                    // populatedChannel, voxelLength, iterIndex, distance );
                    const bool onInterface = fast_marching_signed_distance_to_closest_interface_particle_search(
                        ris, signedDistanceChannel, populatedChannel, voxelLength, iterIndex, i.get_coord(), distance );
                    if( !onInterface ) {
                        throw runtime_error(
                            "fast_marching_initial_band() - The 'Populated' channel provided includes a voxel that is "
                            "not adjacent "
                            "to the interface, an invalid condition for this function.  This is likely caused by "
                            "requesting that "
                            "the initial band distances be updated when they shouldn't be.  The voxel data index is " +
                            lexical_cast<string>( iterIndex ) + "." );
                    }
                }

                if( distance > 0 ) {
                    if( initialSignedDistance <= 0 ) {
                        throw std::runtime_error(
                            "fast_marching_initial_band() Internal Error: interface voxel signed changed from: " +
                            boost::lexical_cast<std::string>( initialSignedDistance ) +
                            " to: " + boost::lexical_cast<std::string>( distance ) );
                    }
                    // add the new distance value to new_values (and later to m_distanceData) if we are updating the
                    // distances
                    newValuesOut.push_back( make_pair( distance, (int)iterIndex ) );
                } else {
                    if( initialSignedDistance > 0 ) {
                        throw std::runtime_error(
                            "fast_marching_initial_band() Internal Error: interface voxel signed changed from: " +
                            boost::lexical_cast<std::string>( initialSignedDistance ) +
                            " to: " + boost::lexical_cast<std::string>( distance ) );
                    }
                    // add the new distance value to new_values (and later to m_distanceData) if we are updating the
                    // distances
                    newValuesIn.push_back( make_pair( -distance, (int)iterIndex ) );
                }

                // mark this point as an initial band point in the populated channel
                populatedChannel[iterIndex] = 5;
            }
            // otherwise we have one of the following cases:
            // - a populatedChannel value of 0, which we will consider a far point
            // - a populatedChannel value of 2, which is denotes that this point should be ignored in the
            // reinitialization
            // - populatedChannel values of 3 and 4 should not be passed into this function, and will be also treated as
            // points to ignore
        }

        // sort the new distance information so that they are sorted when they get added to Near
        sort( newValuesIn.begin(), newValuesIn.end() );
        sort( newValuesOut.begin(), newValuesOut.end() );

        // put the new_values into Near and update the signedDistanceChannel
        for( vector<pair<float, int>>::const_iterator i = newValuesIn.begin(), ie = newValuesIn.end(); i != ie; ++i ) {
            nearInQ.push( *i );
            signedDistanceChannel[i->second] = -i->first;
        }
        for( vector<pair<float, int>>::const_iterator i = newValuesOut.begin(), ie = newValuesOut.end(); i != ie;
             ++i ) {
            nearOutQ.push( *i );
            signedDistanceChannel[i->second] = i->first;
        }

    } else {
        // we are not updating the distances

        // for each defined voxel
        for( size_t iterIndex = 0, iterIndexEnd = adj.data_size(); iterIndex != iterIndexEnd; ++iterIndex ) {

            // check that the signedDistance value is not NAN
            if( frantic::math::is_nan( signedDistanceChannel[iterIndex] ) ) {
                throw runtime_error(
                    "fast_marching_initial_band() - The 'SignedDistance' channel contains a NAN at index " +
                    boost::lexical_cast<string>( iterIndex ) + "." );
            }

            // if we have a populated source
            if( populatedChannel[iterIndex] == 1 ) {

                float distance = signedDistanceChannel[iterIndex];

                if( distance > 0 ) {
                    // add the new distance value to new_values (and later to m_distanceData) if we are updating the
                    // distances
                    newValuesOut.push_back( make_pair( distance, (int)iterIndex ) );
                } else {
                    // add the new distance value to new_values (and later to m_distanceData) if we are updating the
                    // distances
                    newValuesIn.push_back( make_pair( -distance, (int)iterIndex ) );
                }

                // mark this point as an initial band point in the populated channel
                populatedChannel[iterIndex] = 5;
            }
            // otherwise we have one of the following cases:
            // - a populatedChannel value of 0, which we will consider a far point
            // - a populatedChannel value of 2, which is denotes that this point should be ignored in the
            // reinitialization
            // - populatedChannel values of 3 and 4 should not be passed into this function, and will be also treated as
            // points to ignore
        }

        // sort the new distance information so that they are sorted when they get added to Near
        sort( newValuesIn.begin(), newValuesIn.end() );
        sort( newValuesOut.begin(), newValuesOut.end() );

        // put the new_values into Near
        for( vector<pair<float, int>>::const_iterator i = newValuesIn.begin(), ie = newValuesIn.end(); i != ie; ++i ) {
            nearInQ.push( *i );
        }
        for( vector<pair<float, int>>::const_iterator i = newValuesOut.begin(), ie = newValuesOut.end(); i != ie;
             ++i ) {
            nearOutQ.push( *i );
        }
    }
}

/**
 * Taken from fast_marching_initial_band, and modified to handle an unsigned
 * distance field.
 *
 * The 'Populated' channel will have 0 for unpopulated voxels and 1 for populated voxels, and on output
 * it will have 1 for each voxel that was successfully reinitialized. A value of 4 is also used for voxels
 * that are not yet reinitialized but will be when the function terminates (ie in Near). A value of 5 is
 * used to mark the initial band points so that we don't update their distances in the march.
 * There should be no 4 or 5 values in the populated channel on exit.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  voxelLength  The length of one voxel.
 * @param  nearOutQ  The outside portion of the initial band in the queue.
 */
void unsigned_fast_marching_initial_band( const frantic::volumetrics::levelset::ris_adjacency& adj,
                                          unsigned char* populatedChannel, float* signedDistanceChannel,
                                          const float /*voxelLength*/, std::queue<std::pair<float, int>>& nearOutQ ) {

    // set up a vector to temporarily store new distance values
    std::vector<std::pair<float, int>> newValuesOut;

    // It seems reasonable to assume that the initial band will not exceed half of the defined voxels
    // (of course, more memory will be allocated later on if it does)
    newValuesOut.reserve( adj.data_size() / 4 );

    // for each defined voxel
    for( size_t iterIndex = 0, iterIndexEnd = adj.data_size(); iterIndex != iterIndexEnd; ++iterIndex ) {

        // check that the signedDistance value is not NAN
        if( frantic::math::is_nan( signedDistanceChannel[iterIndex] ) ) {
            throw std::runtime_error(
                "unsigned_fast_marching_initial_band() - The 'SignedDistance' channel contains a NAN at index " +
                boost::lexical_cast<std::string>( iterIndex ) + "." );
        }

        // if we have a populated source
        if( populatedChannel[iterIndex] == 1 ) {

            const float distance = signedDistanceChannel[iterIndex];

            if( distance < 0 ) {
                throw std::runtime_error( "unsigned_fast_marching_initial_band() - The 'SignedDistance' channel "
                                          "contains a negative distance " +
                                          boost::lexical_cast<std::string>( distance ) + " at index " +
                                          boost::lexical_cast<std::string>( iterIndex ) + "." );
            }

            newValuesOut.push_back( std::make_pair( distance, (int)iterIndex ) );

            // mark this point as an initial band point in the populated channel
            populatedChannel[iterIndex] = 5;
        }
        // otherwise we have one of the following cases:
        // - a populatedChannel value of 0, which we will consider a far point
        // - a populatedChannel value of 2, which is denotes that this point should be ignored in the reinitialization
        // - populatedChannel values of 3 and 4 should not be passed into this function, and will be also treated as
        // points to ignore
    }

    // sort the new distance information so that they are sorted when they get added to Near
    std::sort( newValuesOut.begin(), newValuesOut.end() );

    // put the new_values into Near
    for( std::vector<std::pair<float, int>>::const_iterator i = newValuesOut.begin(), ie = newValuesOut.end(); i != ie;
         ++i ) {
        nearOutQ.push( *i );
    }
}

/**
 * This function solves the two term quadratic equation
 *
 * @param  voxelLength  The length of one voxel.
 * @param  phi1  The distance from the given voxel to the closest point on the interface along one of axes.
 * @param  phi2  The distance from the given voxel to the closest point on the interface along another axis.
 */
float fast_marching_quad_solve2( const float voxelLength, float phi1, float phi2 ) {
    float a, b, c, result = 0;
    a = 2;
    b = -2 * ( phi1 + phi2 );
    c = ( phi1 * phi1 ) + ( phi2 * phi2 ) - square( voxelLength );

    if( get_quadratic_larger_root( a, b, c, result ) )
        return result;
    else
        return min( phi1, phi2 ) + voxelLength;
}

/**
 * This function solves the three term quadratic equation
 *
 * @param  voxelLength  The length of one voxel.
 * @param  phi1  The distance from the given voxel to the closest point on the interface along the X axis.
 * @param  phi2  The distance from the given voxel to the closest point on the interface along the Y axis.
 * @param  phi3  The distance from the given voxel to the closest point on the interface along the Z axis.
 */
float fast_marching_quad_solve3( const float voxelLength, float phi1, float phi2, float phi3 ) {
    float a, b, c, result = 0;
    a = 3;
    b = -2 * ( phi1 + phi2 + phi3 );
    c = ( phi1 * phi1 + phi2 * phi2 + phi3 * phi3 ) - square( voxelLength );

    if( get_quadratic_larger_root( a, b, c, result ) )
        return result;
    else {
        // If there was no root for the three-term equation, remove the largest phi, and use quad_solve2
        if( phi1 > phi2 ) {
            if( phi1 > phi3 ) {
                return fast_marching_quad_solve2( voxelLength, phi2, phi3 );
            } else {
                return fast_marching_quad_solve2( voxelLength, phi1, phi2 );
            }
        } else {
            if( phi2 > phi3 ) {
                return fast_marching_quad_solve2( voxelLength, phi1, phi3 );
            } else {
                return fast_marching_quad_solve2( voxelLength, phi1, phi2 );
            }
        }
    }
}

/**
 * This function chooses which quadratic function to use.
 *
 * @param  voxelLength  The length of one voxel.
 * @param  phi1  The distance from the given voxel to the closest point on the interface along the X axis.
 * @param  phi2  The distance from the given voxel to the closest point on the interface along the Y axis.
 * @param  phi3  The distance from the given voxel to the closest point on the interface along the Z axis.
 */
float fast_marching_quad_solve_chooser( const float voxelLength, float phi1, float phi2, float phi3 ) {

    int q_func = 0;
    q_func += phi1 == std::numeric_limits<float>::infinity() ? 0 : 1;
    q_func += phi2 == std::numeric_limits<float>::infinity() ? 0 : 2;
    q_func += phi3 == std::numeric_limits<float>::infinity() ? 0 : 4;

    switch( q_func ) {
    case 7:
        return fast_marching_quad_solve3( voxelLength, phi1, phi2, phi3 );
    case 6:
        return fast_marching_quad_solve2( voxelLength, phi2, phi3 );
    case 5:
        return fast_marching_quad_solve2( voxelLength, phi1, phi3 );
    case 3:
        return fast_marching_quad_solve2( voxelLength, phi1, phi2 );
    case 4:
        return phi3 + voxelLength;
    case 2:
        return phi2 + voxelLength;
    case 1:
        return phi1 + voxelLength;
    case 0:
        return std::numeric_limits<float>::infinity();
    default:
        return std::numeric_limits<float>::infinity();
    }
}

/**
 * This function computes the updated distance to the interface of the trial points' neighbours.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  pointIndex  The index of the voxel in question.
 * @param  marchInsideNotOutside  If true, we march the inside, otherwise we march the outside
 * @param  voxelLength  The length of one voxel.
 */
float fast_marching_compute_distance( const ris_adjacency& adj, unsigned char* populatedChannel,
                                      const float* signedDistanceChannel, const int pointIndex,
                                      const bool marchInsideNotOutside, const float voxelLength ) {
    vector3f phi;
    int caseBits;

    const ris_adj_entry& rae = adj[pointIndex];

    for( int axis = 0; axis < 3; ++axis ) {
        int posDataIndex = rae[get_positive_rae_neighbor_index( axis )];
        int negDataIndex = rae[get_negative_rae_neighbor_index( axis )];

        // look at the two neighbours along this axis and find the phi value from an accepted neighbour that is closest
        // to the interface
        caseBits = ( ( posDataIndex >= 0 && populatedChannel[posDataIndex] == 1 ) ? 1 : 0 ) +
                   ( ( negDataIndex >= 0 && populatedChannel[negDataIndex] == 1 ) ? 2 : 0 );

        switch( caseBits ) {
        case 0: // neither neigbour along this axis is accepted
            phi[axis] = std::numeric_limits<float>::infinity();
            break;
        case 1: // only the positive neighbour is accepted
            phi[axis] = signedDistanceChannel[posDataIndex];
            break;
        case 2: // only the negative neighbour is accepted
            phi[axis] = signedDistanceChannel[negDataIndex];
            break;
        default: // both neighbours are accepted, so pick the one that is closer
            if( marchInsideNotOutside ) {
                phi[axis] = std::max<float>( signedDistanceChannel[posDataIndex], signedDistanceChannel[negDataIndex] );
            } else {
                phi[axis] = std::min<float>( signedDistanceChannel[posDataIndex], signedDistanceChannel[negDataIndex] );
            }
        }
    }

    if( marchInsideNotOutside ) {
        return fast_marching_quad_solve_chooser( voxelLength, std::abs( phi.x ), std::abs( phi.y ), std::abs( phi.z ) );
    } else {
        return fast_marching_quad_solve_chooser( voxelLength, phi.x, phi.y, phi.z );
    }
}

/**
 * This function recomputes the distance for the given neighbour (if it exists) and adds it to near if it isn't
 * already in there.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  neighbourIndex  The index of the voxel in question.
 * @param  marchInsideNotOutside  If true, we march the inside, otherwise we march the outside
 * @param  nearPQ  The portion of the marching region corresponding to our initial band in the heap.
 * @param  nearQ  The portion of the marching region corresponding to our initial band in the queue.
 * @param  voxelLength  The length of one voxel.
 * @param  invVoxelLengthSquared  The value 1/(voxelLength^2).
 */
void fast_marching_update_neighbour(
    const ris_adjacency& adj, unsigned char* populatedChannel, const float* signedDistanceChannel,
    const int neighbourIndex, bool marchInsideNotOutside,
    priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>>& nearPQ,
    queue<pair<float, int>>& nearQ, const float voxelLength ) {

    float phi = signedDistanceChannel[neighbourIndex];
    char populatedValue = populatedChannel[neighbourIndex];

    // if the neighbour is already accepted, or is flagged to not
    // update, or is on the other side of the interface, then we
    // don't update it
    if( populatedValue == 1 || populatedValue == 2 ||
        populatedValue ==
            5 ) { // neighbour is accepted or in the initial band (and therefore, accepted) or is flagged to not update
        return;
    }

    // return immediately if the neighbour is on the other side of the interface
    if( marchInsideNotOutside ) {
        if( phi > 0 ) {
            return;
        }
    } else {
        if( phi < 0 ) {
            return;
        }
    }

    float distance = fast_marching_compute_distance( adj, populatedChannel, signedDistanceChannel, neighbourIndex,
                                                     marchInsideNotOutside, voxelLength );

    if( populatedChannel[neighbourIndex] == 0 && ( nearQ.empty() || nearQ.back().first <= distance ) ) {

        // also update the populated state
        populatedChannel[neighbourIndex] = 4;

        // add the neighbour to near, with the new distance value
        nearQ.push( std::pair<float, int>( distance, neighbourIndex ) );
    } else {
        // the neighbour is in near, but the distance needs to be updated. Or, the trial point is from the heap, so we
        // need to put it's neighbours in the heap Note: the point cannot efficiently be removed from a priority_queue,
        // so just add a point with the new distance and later we will check the points in near to make sure that they
        // have the most recent distances before we use them

        // add the point to near, with the new distance value
        nearPQ.push( std::pair<float, int>( distance, neighbourIndex ) ); // nearIn stores positive distances

        if( populatedChannel[neighbourIndex] != 4 ) {
            // also update the populated state
            populatedChannel[neighbourIndex] = 4;
        }
    }
}

/**
 * This function does most of the magic involved in reintializing a level set.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  marchInsideNotOutside  If true, we march the inside, otherwise we march the outside
 * @param  nearPQ  The portion of the marching region corresponding to our initial band in the heap.
 * @param  nearQ  The portion of the marching region corresponding to our initial band in the queue.
 * @param  voxelLength  The length of one voxel.
 * @param  stoppingDistance  The width from the interface at which to stop reinitializing the interface (-1 means no
 * stopping width)
 */
void fast_marching_march( const ris_adjacency& adj, unsigned char* populatedChannel, float* signedDistanceChannel,
                          const bool marchInsideNotOutside,
                          priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>>& nearPQ,
                          queue<pair<float, int>>& nearQ, const float voxelLength, const float stoppingDistance ) {
    std::pair<float, int> frontPair, topPair,
        nullPair = make_pair<float, int>( std::numeric_limits<float>::infinity(), -1 );
    float distance;
    int index;

    // keep going until all the points in near have been accepted
    while( !nearQ.empty() || !nearPQ.empty() ) {
        frantic::logging::check_global_abort();

        frontPair = nearQ.empty() ? nullPair : nearQ.front();
        topPair = nearPQ.empty() ? nullPair : nearPQ.top();

        // get the point from near with the smallest phi value (our trial point)
        // TODO: It might be possible to make things faster by checking if these points are populated here
        // and if so, pop until we find the next unpopulated points at the top/front (or the end)
        if( frontPair.first > topPair.first ) {
            // use the point at the top of the heap, since it is smaller than the front of the queue
            distance = topPair.first;
            index = topPair.second;
            nearPQ.pop();
        } else {
            // use the point at the front of the queue, since it is smaller than (or the same as) the top of the heap
            // Note: if both the queue and the heap are empty, the while loop above will catch it
            distance = frontPair.first;
            index = frontPair.second;
            nearQ.pop();
        }

        // Check to see that the trial point has not already been accepted.
        //	  This case will occur if a point had its distance updated while it was in near
        // because we cannot erase from a priority_queue. We know that the first occurance
        // of a point as the trial point will have the most recent distance data. This means that
        // if we find a point that is already accepted (has been a trial point before) then this
        // point should be ignored since the distance data is out of date.
        if( populatedChannel[index] != 1 ) {

            // stop marching if we are beyond the stopping width
            if( stoppingDistance >= 0 && distance > stoppingDistance ) { // distance will always be positive
                return;
            }

            // mark it as accepted
            populatedChannel[index] = 1;

            // put the accepted value into the original level set's signedDistanceChannel
            if( marchInsideNotOutside ) {
                signedDistanceChannel[index] = -distance; // update with the proper sign
            } else {
                signedDistanceChannel[index] = distance;
            }

            // find the neighbour indices
            for( int face = 0; face < 6; ++face ) {
                const int adjDataIndex = adj[index][face];
                if( adjDataIndex >= 0 ) {
                    fast_marching_update_neighbour( adj, populatedChannel, signedDistanceChannel, adjDataIndex,
                                                    marchInsideNotOutside, nearPQ, nearQ, voxelLength );
                }
            }
        }
    }
}

} // namespace detail

/**
 * This function reintializes a level set.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 * @param  signedDistanceChannel  The raw memory for the 'SignedDistance' channel.
 * @param  voxelLength  The length of one voxel.
 * @param  updatePopulatedInterfaceDistances  Determines whether the intial band should updated the distances of points
 * that are populated and on the interface.
 * @param  insideStoppingDistance  The width at which to stop reinitializing inside the interface (defaults to infinity,
 * which means no stopping width)
 * @param  outsideStoppingDistance  The width at which to stop reinitializing outside the interface (defaults to
 * infinity, which means no stopping width)
 */
void fast_marching_reinitialization( frantic::logging::progress_logger& progressLogger, const rle_index_spec& ris,
                                     unsigned char* populatedChannel, float* signedDistanceChannel,
                                     const float voxelLength, const bool updatePopulatedInterfaceDistances,
                                     const float insideStoppingDistance, const float outsideStoppingDistance ) {
    const ris_adjacency& adj = ris.get_cached_adjacency();

    // the data structures for our near points
    priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>> nearInPQ;
    priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>> nearOutPQ;
    queue<pair<float, int>> nearInQ;
    queue<pair<float, int>> nearOutQ;

    // find the initial band
    detail::fast_marching_initial_band( ris, populatedChannel, signedDistanceChannel, voxelLength,
                                        updatePopulatedInterfaceDistances, nearInQ, nearOutQ );
    progressLogger.update_progress( 20.0f );

    // march out
    detail::fast_marching_march( adj, populatedChannel, signedDistanceChannel, false, nearOutPQ, nearOutQ, voxelLength,
                                 outsideStoppingDistance );
    progressLogger.update_progress( 50.0f );

    // march in
    detail::fast_marching_march( adj, populatedChannel, signedDistanceChannel, true, nearInPQ, nearInQ, voxelLength,
                                 insideStoppingDistance );
    progressLogger.update_progress( 90.0f );

    // set anything that isn't 1 in the polulated channel to 0
    for( size_t i = 0, ie = adj.data_size(); i != ie; ++i ) {
        if( populatedChannel[i] != 1 ) {
            populatedChannel[i] = 0;
        }
    }

    // frantic::volumetrics::levelset::reinitializationProfilingGroup.print();
    progressLogger.update_progress( 100.0f );
}

/**
 * Extrapolate an unsigned distance field away from the initially populated
 * voxels.
 *
 * This function is based on fast_marching_reinitialization and
 * uses nearly the same procedure, but here the extrapolation is from the
 * initially populated voxels, and there is no check for zero crossings or
 * any changes in the initially populated band.
 *
 * The unsigned distance channel is expected to hold non-negative distances
 * in its populated voxels.
 *
 * @param  adj  A filled in ris_adjacency structure which provides the necessary adjacency information.
 * @param  populatedChannel  The raw memory for the 'Populated' channel.
 *		   This should be set to 1 for voxels with correct distance that can
 *		   be used as an extrapolation source; set to 0 for voxels with invalid
 *		   distances that should be set by extrapolation; and set to 2 for voxels that
 *		   should be ignored during extrapolation (their distance is not used
 *		   as a source and is not changed during extrapolation).
 *		   All other values are reserved.
 *		   On return, this will be 1 for voxels that were successfully extrapolated into
 *		   or were initially populated, and 0 for voxels that were not extrapolated into
 *		   or that were set to ignore.
 * @param  signedDistanceChannel  The raw memory for the unsigned distance channel.
 * @param  voxelLength  The length of one voxel.
 * @param  stoppingDistance  The distance at which to stop reinitializing the distance field (defaults to infinity,
 *which means no stopping distance)
 */
void extrapolate_unsigned_distance( const frantic::volumetrics::levelset::ris_adjacency& adj,
                                    unsigned char* populatedChannel, float* signedDistanceChannel,
                                    const float voxelLength, const float stoppingDistance ) {
    // the data structures for our near points
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>>
        nearOutPQ;
    std::queue<std::pair<float, int>> nearOutQ;

    /*
    {
      size_t populatedCount = 0, totalCount = adj.data_size();
      for( std::size_t i=0, ie = adj.data_size(); i != ie; ++i ) {
        if( populatedChannel[ i ] == 1 ) {
          ++populatedCount;
        }
      }
      if( totalCount > 0 )
        FF_LOG( stats ) << "extrapolate_unsigned_distance() Before unsigned distance extrapolation, populated fraction :
    "
    << double( populatedCount ) / totalCount << std::endl; else FF_LOG( stats ) << "extrapolate_unsigned_distance() zero
    voxels in level set" << std::endl;
    }
    */

    // find the initial band
    detail::unsigned_fast_marching_initial_band( adj, populatedChannel, signedDistanceChannel, voxelLength, nearOutQ );

    // FF_LOG( debug ) << "extrapolate_unsigned_distance() Number of voxels in the initial queue: " << nearOutQ.size()
    // << std::endl;

    // march out
    detail::fast_marching_march( adj, populatedChannel, signedDistanceChannel, false, nearOutPQ, nearOutQ, voxelLength,
                                 stoppingDistance );

    // set anything that isn't 1 in the populated channel to 0
    // size_t populatedCount = 0, totalCount = adj.data_size();
    for( std::size_t i = 0, ie = adj.data_size(); i != ie; ++i ) {
        if( populatedChannel[i] != 1 ) {
            populatedChannel[i] = 0;
        }
        // else {
        //	++populatedCount;
        // }
    }
    // if( totalCount > 0 )
    //	FF_LOG( stats ) << "extrapolate_unsigned_distance() After unsigned distance extrapolation, populated fraction :
    //"
    //<< double( populatedCount ) / totalCount << std::endl;
}

float solve_fast_marching_first_order_quad( const float voxelLength, const float phi1, const float phi2,
                                            const float phi3 ) {
    return detail::fast_marching_quad_solve_chooser( voxelLength, phi1, phi2, phi3 );
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
