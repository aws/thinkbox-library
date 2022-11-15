// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

using namespace frantic::graphics;
using namespace frantic::volumetrics::levelset;

namespace frantic {
namespace volumetrics {
namespace levelset {

/// The following are forward declarations of internal functions used in this file
namespace detail {
unsigned int get_defined_voxel_count( const boost::int32_t* indices );
void create_2x2x2_indicies_array( const vector3f& voxelCoordinate, const boost::int32_t* const inIndices,
                                  const boundbox3& indicesVoxelBounds, boost::int32_t* outIndices, float* outDeltas );
void trilerp_float_from_2x2x2_indices( const float* dataAccessor, const boost::int32_t* indices, float* deltas,
                                       float* outTrilerp, vector3f* outGradient );
void trilerp_float_from_2x2x2_indices_with_undefined( const float* dataAccessor, const boost::int32_t* indices,
                                                      float* deltas, float* outTrilerp, vector3f* outGradient );
void trilerp_vector3f_from_2x2x2_indices( const_rle_channel_accessor<vector3f> dataAccessor,
                                          const boost::int32_t* indices, float* deltas, vector3f& outTrilerp );
void trilerp_vector3f_from_2x2x2_indices_with_undefined( const_rle_channel_accessor<vector3f> dataAccessor,
                                                         const boost::int32_t* indices, float* deltas,
                                                         vector3f& outTrilerp );
void trilerp_general_channel_from_2x2x2_indices( const_rle_channel_general_accessor dataAccessor,
                                                 const boost::int32_t* indices, float* deltas, char* outTrilerp );
} // namespace detail

bool trilerp_float( const rle_index_spec& ris, const float* dataAccessor, const vector3f& voxelCoordinate,
                    float* outTrilerp, vector3f* outGradient ) {
    boost::int32_t indices[8];
    float deltas[3];
    get_trilerp_indices( ris, voxelCoordinate, deltas, indices );

    unsigned int definedVoxelCount = detail::get_defined_voxel_count( indices );
    if( definedVoxelCount == 0 ) {
        if( outTrilerp )
            *outTrilerp =
                ( indices[0] == -1 ) ? std::numeric_limits<float>::infinity() : -std::numeric_limits<float>::infinity();
        if( outGradient )
            outGradient->set( 0.0f );
        return false;
    }

    if( definedVoxelCount == 8 )
        detail::trilerp_float_from_2x2x2_indices( dataAccessor, indices, deltas, outTrilerp, outGradient );
    else
        detail::trilerp_float_from_2x2x2_indices_with_undefined( dataAccessor, indices, deltas, outTrilerp,
                                                                 outGradient );

    return true;
}

bool trilerp_float_from_indices( const float* dataAccessor, const vector3f& voxelCoordinate,
                                 const boost::int32_t* const inIndices, const boundbox3& indicesVoxelBounds,
                                 float* outTrilerp, vector3f* outGradient ) {
    boost::int32_t indices[8];
    float deltas[3];
    detail::create_2x2x2_indicies_array( voxelCoordinate, inIndices, indicesVoxelBounds, indices, deltas );

    unsigned int definedVoxelCount = detail::get_defined_voxel_count( indices );
    if( definedVoxelCount == 0 ) {
        if( outTrilerp )
            *outTrilerp =
                ( indices[0] == -1 ) ? std::numeric_limits<float>::infinity() : -std::numeric_limits<float>::infinity();
        if( outGradient )
            outGradient->set( 0.0f );
        return false;
    }

    if( definedVoxelCount == 8 )
        detail::trilerp_float_from_2x2x2_indices( dataAccessor, indices, deltas, outTrilerp, outGradient );
    else
        detail::trilerp_float_from_2x2x2_indices_with_undefined( dataAccessor, indices, deltas, outTrilerp,
                                                                 outGradient );

    return true;
}

bool trilerp_vector3f( const rle_index_spec& ris, const_rle_channel_accessor<vector3f> dataAccessor,
                       const vector3f& voxelCoordinate, vector3f& outTrilerp ) {
    boost::int32_t indices[8];
    float deltas[3];
    get_trilerp_indices( ris, voxelCoordinate, deltas, indices );

    unsigned int definedVoxelCount = detail::get_defined_voxel_count( indices );
    if( definedVoxelCount == 0 ) {
        outTrilerp.set( 0.0f );
        return false;
    }

    if( definedVoxelCount == 8 )
        detail::trilerp_vector3f_from_2x2x2_indices( dataAccessor, indices, deltas, outTrilerp );
    else
        detail::trilerp_vector3f_from_2x2x2_indices_with_undefined( dataAccessor, indices, deltas, outTrilerp );

    return true;
}

bool trilerp_vector3f_from_indices( const_rle_channel_accessor<vector3f> dataAccessor, const vector3f& voxelCoordinate,
                                    const boost::int32_t* inIndices, const boundbox3& indicesVoxelBounds,
                                    vector3f& outTrilerp ) {
    boost::int32_t indices[8];
    float deltas[3];
    detail::create_2x2x2_indicies_array( voxelCoordinate, inIndices, indicesVoxelBounds, indices, deltas );

    unsigned int definedVoxelCount = detail::get_defined_voxel_count( indices );
    if( definedVoxelCount == 0 ) {
        outTrilerp.set( 0.0f );
        return false;
    }

    if( definedVoxelCount == 8 )
        detail::trilerp_vector3f_from_2x2x2_indices( dataAccessor, indices, deltas, outTrilerp );
    else
        detail::trilerp_vector3f_from_2x2x2_indices_with_undefined( dataAccessor, indices, deltas, outTrilerp );

    return true;
}

bool trilerp_vector3f_staggered( const rle_index_spec& ris, const_rle_channel_accessor<vector3f> dataAccessor,
                                 const vector3f& voxelCoordinate, vector3f& outTrilerp ) {
    vector3 voxelCoordFloored = vector3::from_floor( voxelCoordinate );
    boundbox3 indicesBounds( voxelCoordFloored, voxelCoordFloored + vector3( 1 ) );
    if( voxelCoordinate.x - voxelCoordFloored.x < 0.5f )
        --indicesBounds.minimum().x;
    if( voxelCoordinate.y - voxelCoordFloored.y < 0.5f )
        --indicesBounds.minimum().y;
    if( voxelCoordinate.z - voxelCoordFloored.z < 0.5f )
        --indicesBounds.minimum().z;

    std::vector<boost::int32_t> dataIndices( indicesBounds.get_volume() );
    ris.fill_data_index_box( indicesBounds, &dataIndices[0] );

    return trilerp_vector3f_staggered_from_indices( dataAccessor, voxelCoordinate, &dataIndices[0], indicesBounds,
                                                    outTrilerp );
}

namespace trilerp_vector3f_staggered_detail {
void assert_valid_abc( const vector3& abc, const vector3f& voxelCoordinate, const vector3& xyz,
                       const vector3f& currentXYZMin ) {
    for( int i = 0; i < 3; ++i ) {
        if( abc[i] < 0 ) {
            logging::error << "Box Index is Invalid!\n";
            // logging::error << "Face: " << face << std::endl;
            logging::error << "Voxel Lookup: " << voxelCoordinate << std::endl;
            logging::error << "xyz: " << xyz << std::endl;
            logging::error << "currentXYZMIn: " << currentXYZMin << std::endl;
            logging::error << "abc: " << abc << std::endl;
            throw std::runtime_error( "trilerp_vector3f_staggered_from_indices: Box Index is Invalid!" );
        }
    }
}

void clamp_weights( vector3f& weights, const vector3f& voxelCoordinate, const vector3& xyz ) {
    // setting the threshold to 0.0001 * voxelsize which should be enough to capture and floating remains from the
    // offset subtraction
    const float zeroThreshold = 1e-4f;
    for( int i = 0; i < 3; ++i ) {
        float& wRef = weights[i];
        if( wRef > 1.f || wRef < -zeroThreshold ) {
            // just a little debug assertion that might be useful for a bit
            logging::error << "i: " << i << std::endl;
            logging::error << "Weight: " << wRef << std::endl;
            // logging::error << "diff: " << offsetLookup << std::endl;
            logging::error << "xyz: " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << std::endl;
            throw std::runtime_error( "trilerp_vector3f_staggered_clamp_weight: invalid weight calculated " +
                                      boost::lexical_cast<std::string>( wRef ) + " for voxel " +
                                      voxelCoordinate.str() );
        } else if( wRef < zeroThreshold ) {
            wRef = 0;
        }
    }
}
} // namespace trilerp_vector3f_staggered_detail

bool trilerp_vector3f_staggered_from_indices( const_rle_channel_accessor<vector3f> dataAccessor,
                                              const vector3f& voxelCoordinate, const boost::int32_t* indices,
                                              const boundbox3& indicesVoxelBounds, vector3f& outTrilerp ) {
    const vector3 currentXYZMin = indicesVoxelBounds.minimum();
    const size3 boxSize = indicesVoxelBounds.size();
    float weights[8];
    vector3 abc;
    vector3f w;

    bool returnValue = true;

    const vector3f voxelCoordinateFloor = vector3f::from_floor( voxelCoordinate );
    const vector3 xyzCenter( static_cast<boost::int32_t>( voxelCoordinateFloor[0] ),
                             static_cast<boost::int32_t>( voxelCoordinateFloor[1] ),
                             static_cast<boost::int32_t>( voxelCoordinateFloor[2] ) );
    const vector3 abcCenter( xyzCenter - currentXYZMin );
    trilerp_vector3f_staggered_detail::assert_valid_abc( abcCenter, voxelCoordinate, xyzCenter, currentXYZMin );

    vector3f wCenter( voxelCoordinate - voxelCoordinateFloor );
    trilerp_vector3f_staggered_detail::clamp_weights( wCenter, voxelCoordinate, xyzCenter );

    const vector3f offsetLookupFace( voxelCoordinate - vector3f( 0.5f ) );
    const vector3f offsetLookupFaceFloor = vector3f::from_floor( offsetLookupFace );

    const vector3 xyzFace( static_cast<boost::int32_t>( offsetLookupFaceFloor[0] ),
                           static_cast<boost::int32_t>( offsetLookupFaceFloor[1] ),
                           static_cast<boost::int32_t>( offsetLookupFaceFloor[2] ) );
    const vector3 abcFace( xyzFace - currentXYZMin );
    trilerp_vector3f_staggered_detail::assert_valid_abc( abcFace, voxelCoordinate, xyzFace, currentXYZMin );

    vector3f wFace( offsetLookupFace - offsetLookupFaceFloor );
    trilerp_vector3f_staggered_detail::clamp_weights( wFace, voxelCoordinate, xyzFace );

    for( int face = 0; face < 3; ++face ) {
        switch( face ) {
        case 0:
            w.set( wCenter[0], wFace[1], wFace[2] );
            abc.set( abcCenter[0], abcFace[1], abcFace[2] );
            break;
        case 1:
            w.set( wFace[0], wCenter[1], wFace[2] );
            abc.set( abcFace[0], abcCenter[1], abcFace[2] );
            break;
        case 2:
            w.set( wFace[0], wFace[1], wCenter[2] );
            abc.set( abcFace[0], abcFace[1], abcCenter[2] );
            break;
        }

        // build the 8 weights for the interpolation
        weights[0] = ( 1 - w.x ) * ( 1 - w.y ) * ( 1 - w.z );
        weights[1] = ( w.x ) * ( 1 - w.y ) * ( 1 - w.z );
        weights[2] = ( 1 - w.x ) * ( w.y ) * ( 1 - w.z );
        weights[3] = ( w.x ) * ( w.y ) * ( 1 - w.z );
        weights[4] = ( 1 - w.x ) * ( 1 - w.y ) * ( w.z );
        weights[5] = ( w.x ) * ( 1 - w.y ) * ( w.z );
        weights[6] = ( 1 - w.x ) * ( w.y ) * ( w.z );
        weights[7] = ( w.x ) * ( w.y ) * ( w.z );

        float& resultComponent = outTrilerp[face];
        float* currentW = &weights[0];

        resultComponent = 0;
        float weightSum = 0.f;
        unsigned int definedCount = 0;
        for( int k = 0; k < 2; ++k ) {
            int zIndexOffset = ( abc.z + k ) * boxSize.xsize() * boxSize.ysize();
            for( int j = 0; j < 2; ++j ) {
                int yzIndexOffset = ( abc.y + j ) * boxSize.xsize() + zIndexOffset;
                for( int i = 0; i < 2; ++i ) {
                    int boxIndex = ( abc.x + i ) + yzIndexOffset;
                    boost::int32_t dataIndex = indices[boxIndex];
                    if( dataIndex >= 0 ) {
                        const float& weight = *currentW;
                        weightSum += weight;
                        resultComponent += dataAccessor[dataIndex][face] * weight;
                        ++definedCount;
                    }
                    ++currentW;
                }
            }
        }

        // divide out the weight sum.
        // if there's no weight, it will return 0. this is a little different than the non-staggered version (which
        // returns an even weighting of all defined voxels). possibly make it consistent.
        if( weightSum < 1e-20f )
            resultComponent = 0.0f;
        else if( definedCount < 8 )
            resultComponent /= weightSum;

        // if there's no defined voxels for this component, we will alert the user by returning false.
        if( definedCount == 0 )
            returnValue = false;
    }

    return returnValue;
}

bool trilerp_general_channel( const rle_index_spec& ris, const_rle_channel_general_accessor dataAccessor,
                              const vector3f& voxelCoordinate, char* outTrilerp ) {
    boost::int32_t indices[8];
    float deltas[3];
    get_trilerp_indices( ris, voxelCoordinate, deltas, indices );

    unsigned int definedVoxelCount = detail::get_defined_voxel_count( indices );
    if( definedVoxelCount == 0 ) {
        memset( outTrilerp, 0, dataAccessor.primitive_size() );
        return false;
    }

    detail::trilerp_general_channel_from_2x2x2_indices( dataAccessor, indices, deltas, outTrilerp );
    return true;
}

bool trilerp_general_channel_from_indices( const_rle_channel_general_accessor dataAccessor,
                                           const vector3f& voxelCoordinate, const boost::int32_t* inIndices,
                                           const boundbox3& indicesVoxelBounds, char* outTrilerp ) {
    boost::int32_t indices[8];
    float deltas[3];
    detail::create_2x2x2_indicies_array( voxelCoordinate, inIndices, indicesVoxelBounds, indices, deltas );

    unsigned int definedVoxelCount = detail::get_defined_voxel_count( indices );
    if( definedVoxelCount == 0 ) {
        memset( outTrilerp, 0, dataAccessor.primitive_size() );
        return false;
    }

    detail::trilerp_general_channel_from_2x2x2_indices( dataAccessor, indices, deltas, outTrilerp );
    return true;
}

void get_trilerp_indices( const rle_index_spec& ris, const vector3f& voxelCoordinate, float* outTrilerpDeltas,
                          boost::int32_t* outDataIndices ) {
    // The data samples are at the center of each integer voxel, so we have to offset by 0.5 to compensate.
    float dx = voxelCoordinate.x - 0.5f, dy = voxelCoordinate.y - 0.5f, dz = voxelCoordinate.z - 0.5f;
    vector3f voxelMin( floorf( dx ), floorf( dy ), floorf( dz ) );
    dx -= voxelMin.x;
    dy -= voxelMin.y;
    dz -= voxelMin.z;

    // Get all 8 data index values.
    vector3 minCornerXYZ( (int)voxelMin.x, (int)voxelMin.y, (int)voxelMin.z );
    ris.fill_2x2x2_data_index_box( minCornerXYZ, outDataIndices );

    // Fill in all the multiplier weights.
    outTrilerpDeltas[0] = dx;
    outTrilerpDeltas[1] = dy;
    outTrilerpDeltas[2] = dz;
}

void get_trilerp_weights( const float* deltas, float* outWeights ) {
    // Get 1 - weights
    float dxInv = 1.0f - deltas[0];
    float dyInv = 1.0f - deltas[1];
    float dzInv = 1.0f - deltas[2];

    // Calculate 4 combinations of dx, dy, dxInv and dyInv.
    // This saves 4 float calculations every time method is called.
    float dxDy = deltas[0] * deltas[1];
    float dxInvDy = dxInv * deltas[1];
    float dxDyInv = deltas[0] * dyInv;
    float dxInvDyInv = dxInv * dyInv;

    // Fill weights array
    outWeights[0] = dxInvDyInv * dzInv;
    outWeights[1] = dxDyInv * dzInv;
    outWeights[2] = dxInvDy * dzInv;
    outWeights[3] = dxDy * dzInv;
    outWeights[4] = dxInvDyInv * deltas[2];
    outWeights[5] = dxDyInv * deltas[2];
    outWeights[6] = dxInvDy * deltas[2];
    outWeights[7] = dxDy * deltas[2];
}

void normalize_trilerp_weights_for_undefined_voxels( const boost::int32_t* dataIndices, float* outWeights ) {
    // count missing voxels
    int missingCount = 0;
    for( int i = 0; i < 8; ++i ) {
        if( dataIndices[i] < 0 ) {
            ++missingCount;
            outWeights[i] = 0.0f;
        }
    }

    // don't deal with zero missing voxels, and we can't deal when all 8 are missing.
    if( missingCount == 0 || missingCount == 8 )
        return;

    // we've zeroed out all the undefined-voxel weights. likely now the weights will not sum to 1.0f anymore, so we need
    // to fix that.
    float sum = 0.0f;
    for( int i = 0; i < 8; ++i )
        sum += outWeights[i];

    if( sum > 1e-20f ) {
        // divide by the sum, which ensures that the sum of the weights will equal 1.0f.
        for( int i = 0; i < 8; ++i )
            outWeights[i] /= sum;
    } else {
        // if the sum is infinitly small, dividing by an infinity small number will produce floating point error.
        // the following will treat each defined voxel as equal weight in the case where none of them really contribute
        // (does this make sense?).
        float weight = 1.0f / ( 8 - missingCount );
        for( int i = 0; i < 8; ++i )
            if( dataIndices[i] >= 0 )
                outWeights[i] = weight;
    }
}

void trilerp_staggered_centered_signed_distance_gradient( const rle_level_set& levelSet,
                                                          const vector3f& voxelCoordinate,
                                                          vector3f& outSignedDistanceGradient ) {
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
    float xStaggeredFloor = floorf( voxelCoordinate.x );
    float yStaggeredFloor = floorf( voxelCoordinate.y );
    float zStaggeredFloor = floorf( voxelCoordinate.z );
    float xStaggeredAlpha = voxelCoordinate.x - xStaggeredFloor;
    float yStaggeredAlpha = voxelCoordinate.y - yStaggeredFloor;
    float zStaggeredAlpha = voxelCoordinate.z - zStaggeredFloor;

    // We want the 3x3x3 voxels bracketing the middle one, so subtract 1.
    int xstart = (int)xStaggeredFloor - 1;
    int ystart = (int)yStaggeredFloor - 1;
    int zstart = (int)zStaggeredFloor - 1;

    vector3 minimum( xstart, ystart, zstart );
    vector3 maximum( xstart + 2, ystart + 2, zstart + 2 );

    boost::int32_t dataIndices[27];
    levelSet.get_rle_index_spec().fill_data_index_box( boundbox3( minimum, maximum ), dataIndices );

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

    // xstart is xRegFloor or xRegFloor+1, therefore xOff is 0 or 1.
    int xOff = (int)xRegFloor - xstart;
    int yOff = (int)yRegFloor - ystart;
    int zOff = (int)zRegFloor - zstart;

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
    //  The a are staggered differences in the x-axis, the � are staggered differences in the y-axis.
    //  The + is the desired sample point. Notice that the bracketing samples are different for the
    //  x, and y computations.
    //  ------- ------- -------
    //|       |       |       |
    //|   x   a   x   a   x   |
    //|       |       |       |
    //  ------- ---�--- ---�---
    //|       |     + |       |
    //|   x   a   x   a   x   |
    //|       |       |       |
    //  ------- ---�--- ---�---
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
        float d1, d2, d3, w;
        int baseIndex = yOff * 3 + zOff * 9;

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex] );     //  x,  y,  z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 1] ); // x+1,  y,  z
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 2] ); // x+2,  y,  z
        w = ( 1.f - zRegAlpha ) * ( 1.f - yRegAlpha );
        outSignedDistanceGradient.x += w * ( ( d2 - d1 ) + xStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 3] ); //  x, y+1, z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 4] ); // x+1, y+1, z
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 5] ); // x+2, y+1, z
        w = ( 1.f - zRegAlpha ) * yRegAlpha;
        outSignedDistanceGradient.x += w * ( ( d2 - d1 ) + xStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 9] );  // x,   y, z+1
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 10] ); // x+1, y, z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 11] ); // x+2, y, z+1
        w = zRegAlpha * ( 1.f - yRegAlpha );
        outSignedDistanceGradient.x += w * ( ( d2 - d1 ) + xStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 12] ); // x,   y+1, z+1
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 13] ); // x+1, y+1, z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 14] ); // x+2, y+1, z+1
        w = zRegAlpha * yRegAlpha;
        outSignedDistanceGradient.x += w * ( ( d2 - d1 ) + xStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );
    }

    // For the y component, compute the 8 bracketing y-staggered differences.
    outSignedDistanceGradient.y = 0;
    {
        float d1, d2, d3, w;
        int baseIndex = xOff + zOff * 9;

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex] );     // x,  y  ,  z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 3] ); // x,  y+1,  z
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 6] ); // x,  y+2,  z
        w = ( 1.f - zRegAlpha ) * ( 1.f - xRegAlpha );
        outSignedDistanceGradient.y += w * ( ( d2 - d1 ) + yStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 1] ); // x+1,  y  ,  z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 4] ); // x+1,  y+1,  z
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 7] ); // x+1,  y+2,  z
        w = ( 1.f - zRegAlpha ) * xRegAlpha;
        outSignedDistanceGradient.y += w * ( ( d2 - d1 ) + yStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 9] );  // x,  y  ,  z+1
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 12] ); // x,  y+1,  z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 15] ); // x,  y+2,  z+1
        w = zRegAlpha * ( 1.f - xRegAlpha );
        outSignedDistanceGradient.y += w * ( ( d2 - d1 ) + yStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 10] ); // x+1,  y  ,  z+1
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 13] ); // x+1,  y+1,  z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 16] ); // x+1,  y+2,  z+1
        w = zRegAlpha * xRegAlpha;
        outSignedDistanceGradient.y += w * ( ( d2 - d1 ) + yStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );
    }

    // For the z component, compute the 8 bracketing z-staggered differences.
    outSignedDistanceGradient.z = 0;
    {
        float d1, d2, d3, w;
        int baseIndex = xOff + yOff * 3;

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex] );      // x,  y,  z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 9] );  // x,  y,  z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 18] ); // x,  y,  z+2
        w = ( 1.f - yRegAlpha ) * ( 1.f - xRegAlpha );
        outSignedDistanceGradient.z += w * ( ( d2 - d1 ) + zStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 1] );  // x+1,  y,  z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 10] ); // x+1,  y,  z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 19] ); // x+1,  y,  z+2
        w = ( 1.f - yRegAlpha ) * xRegAlpha;
        outSignedDistanceGradient.z += w * ( ( d2 - d1 ) + zStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 3] );  // x,  y+1,  z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 12] ); // x,  y+1,  z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 21] ); // x,  y+1,  z+2
        w = yRegAlpha * ( 1.f - xRegAlpha );
        outSignedDistanceGradient.z += w * ( ( d2 - d1 ) + zStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );

        d1 = levelSet.get_using_data_index( dataIndices[baseIndex + 4] );  // x+1,  y+1,  z
        d2 = levelSet.get_using_data_index( dataIndices[baseIndex + 13] ); // x+1,  y+1,  z+1
        d3 = levelSet.get_using_data_index( dataIndices[baseIndex + 22] ); // x+1,  y+1,  z+2
        w = yRegAlpha * xRegAlpha;
        outSignedDistanceGradient.z += w * ( ( d2 - d1 ) + zStaggeredAlpha * ( d3 - 2 * d2 + d1 ) );
    }
}

void cached_trilerp::fill_cache( const frantic::graphics::vector3f& voxelCoordinate ) {
    if( voxelCoordinate != m_cachedVoxelCoordinate ) {
        const vector3f offsetVoxelCoordinate( voxelCoordinate.x - 0.5f, voxelCoordinate.y - 0.5f,
                                              voxelCoordinate.z - 0.5f );
        const vector3f minCornerCoordinate( floorf( offsetVoxelCoordinate.x ), floorf( offsetVoxelCoordinate.y ),
                                            floorf( offsetVoxelCoordinate.z ) );
        const vector3 minCorner( boost::int32_t( minCornerCoordinate.x ), boost::int32_t( minCornerCoordinate.y ),
                                 boost::int32_t( minCornerCoordinate.z ) );

        if( minCorner != m_cachedMinCorner ) {
            if( minCorner != m_cachedMinCorner ) {
                m_ris.fill_2x2x2_data_index_box( minCorner, m_dataIndices );
                boost::int32_t definedCount = 0;
                for( boost::int32_t* i = m_dataIndices; i != m_dataIndices + 8; ++i ) {
                    if( *i >= 0 ) {
                        m_definedDataIndices[definedCount++] = *i;
                    }
                }
                m_definedCount = definedCount;
                m_cachedMinCorner = minCorner;
            }
        }

        const vector3f deltas = offsetVoxelCoordinate - minCornerCoordinate;
        get_trilerp_weights( &deltas.x, m_weights );

        memcpy( m_definedWeights, m_weights, 8 * sizeof( float ) );

        if( m_definedCount != 8 && m_definedCount != 0 ) {
            normalize_trilerp_weights_for_undefined_voxels( m_dataIndices, m_definedWeights );

            float* iterInputWeight = m_definedWeights;
            float* iterOutputWeight = m_definedWeights;
            for( boost::int32_t* i = m_dataIndices; i != m_dataIndices + 8; ++i, ++iterInputWeight ) {
                if( *i >= 0 ) {
                    if( iterInputWeight != iterOutputWeight ) {
                        *iterOutputWeight = *iterInputWeight;
                    }
                    ++iterOutputWeight;
                }
            }
        }

        m_cachedVoxelCoordinate = voxelCoordinate;
    }
}

cached_trilerp::cached_trilerp( const rle_index_spec& ris )
    : m_ris( ris )
    , m_definedCount( 0 )
    , m_cachedVoxelCoordinate( std::numeric_limits<float>::infinity() )
    , m_cachedMinCorner( std::numeric_limits<boost::int32_t>::max() ) {}

bool cached_trilerp::get( const const_rle_channel_general_accessor& dataAccessor,
                          const frantic::graphics::vector3f& voxelCoordinate, char* outSample ) {
    assert( outSample );
    assert( dataAccessor.size() == m_ris.data_size() );

    return get( dataAccessor.data( 0 ), dataAccessor.arity(), dataAccessor.data_type(), voxelCoordinate, outSample );
}

bool cached_trilerp::get( const char* dataAccessor, const std::size_t arity, const data_type_t dataType,
                          const frantic::graphics::vector3f& voxelCoordinate, char* outSample ) {
    assert( outSample );

    fill_cache( voxelCoordinate );

    if( m_definedCount == 0 )
        return false;

    assert( dataAccessor );

    const std::size_t stride = arity * frantic::channels::sizeof_channel_data_type( dataType );
    // TODO : may be worthwhile to write a function that calculates the weighted sum from indices
    // directly, rather than require a char*
    frantic::channels::channel_weighted_sum_combine_function_t weightedSumFunction =
        channel_weighted_sum_combine_function( dataType );

    const char* data[8];
    for( size_t i = 0; i < (size_t)m_definedCount; ++i ) {
        data[i] = dataAccessor + stride * m_definedDataIndices[i];
    }

    weightedSumFunction( m_definedWeights, &data[0], m_definedCount, arity, outSample );

    return true;
}

float cached_trilerp::get( const rle_level_set& ls, const frantic::graphics::vector3f& voxelCoordinate ) {
    assert( ls.size() == m_ris.data_size() );

    fill_cache( voxelCoordinate );

    float result = 0;
    if( m_definedCount == 0 ) {
        float vote = 0;
        float* iterWeight = m_weights;
        for( boost::int32_t* iterDataIndex = m_dataIndices; iterDataIndex != m_dataIndices + 8;
             ++iterDataIndex, ++iterWeight ) {
            if( ( *iterDataIndex ) == -1 ) {
                vote += ( *iterWeight );
            } else {
                vote -= ( *iterWeight );
            }
        }
        result = ( vote <= 0 ) ? ls.get_inside_distance() : ls.get_outside_distance();
    } else {
        float* iterWeight = m_definedWeights;
        for( boost::int32_t* iterDataIndex = m_definedDataIndices;
             iterDataIndex != m_definedDataIndices + m_definedCount; ++iterDataIndex, ++iterWeight ) {
            result += ( *iterWeight ) * ls[*iterDataIndex];
        }
    }
    return result;
}

frantic::graphics::boundbox3f
cached_trilerp_vector3f_staggered::get_trilerp_vector3f_staggered_sample_range_from_indices(
    const const_rle_channel_accessor<vector3f>& dataAccessor, const vector3f& voxelCoordinate,
    const boost::int32_t* indices, const boundbox3& indicesVoxelBounds ) {
    const size3 boxSize = indicesVoxelBounds.size();
    vector3 abc;

    boundbox3f clampBounds;
    vector3f& clampMin = clampBounds.minimum();
    vector3f& clampMax = clampBounds.maximum();

    const vector3 currentXYZMin = indicesVoxelBounds.minimum();

    const vector3f voxelCoordinateFloor = vector3f::from_floor( voxelCoordinate );
    const vector3 xyzCenter( static_cast<boost::int32_t>( voxelCoordinateFloor[0] ),
                             static_cast<boost::int32_t>( voxelCoordinateFloor[1] ),
                             static_cast<boost::int32_t>( voxelCoordinateFloor[2] ) );
    const vector3 abcCenter( xyzCenter - currentXYZMin );

    const vector3f offsetLookupFace( voxelCoordinate - vector3f( 0.5f ) );
    const vector3f offsetLookupFaceFloor = vector3f::from_floor( offsetLookupFace );

    const vector3 xyzFace( static_cast<boost::int32_t>( offsetLookupFaceFloor[0] ),
                           static_cast<boost::int32_t>( offsetLookupFaceFloor[1] ),
                           static_cast<boost::int32_t>( offsetLookupFaceFloor[2] ) );
    const vector3 abcFace( xyzFace - currentXYZMin );

    for( int axis = 0; axis < 3; ++axis ) {
        switch( axis ) {
        case 0:
            abc.set( abcCenter[0], abcFace[1], abcFace[2] );
            break;
        case 1:
            abc.set( abcFace[0], abcCenter[1], abcFace[2] );
            break;
        case 2:
            abc.set( abcFace[0], abcFace[1], abcCenter[2] );
            break;
        }

        for( int k = 0; k < 2; ++k ) {
            const int zIndexOffset = ( abc.z + k ) * boxSize.xsize() * boxSize.ysize();
            for( int j = 0; j < 2; ++j ) {
                const int yzIndexOffset = ( abc.y + j ) * boxSize.xsize() + zIndexOffset;
                for( int i = 0; i < 2; ++i ) {
                    const int boxIndex = ( abc.x + i ) + yzIndexOffset;
                    const boost::int32_t dataIndex = indices[boxIndex];
                    if( dataIndex >= 0 ) {
                        const float val = dataAccessor[dataIndex][axis];
                        if( val < clampMin[axis] ) {
                            clampMin[axis] = val;
                        }
                        if( val > clampMax[axis] ) {
                            clampMax[axis] = val;
                        }
                    }
                }
            }
        }
    }

    return clampBounds;
}

bool cached_trilerp_vector3f_staggered::get( const const_rle_channel_accessor<vector3f>& dataAccessor,
                                             const vector3f& voxelCoordinate, vector3f& outSample ) {
    const vector3 voxelCoordFloored = vector3::from_floor( voxelCoordinate );
    boundbox3 indicesBounds( voxelCoordFloored, voxelCoordFloored + vector3( 1 ) );
    if( voxelCoordinate.x - voxelCoordFloored.x < 0.5f )
        --indicesBounds.minimum().x;
    if( voxelCoordinate.y - voxelCoordFloored.y < 0.5f )
        --indicesBounds.minimum().y;
    if( voxelCoordinate.z - voxelCoordFloored.z < 0.5f )
        --indicesBounds.minimum().z;

    if( indicesBounds != m_cachedBounds ) {
        m_cachedBounds = indicesBounds;
        m_ris.fill_data_index_box( indicesBounds, m_dataIndices );
    }

    return trilerp_vector3f_staggered_from_indices( dataAccessor, voxelCoordinate, m_dataIndices, indicesBounds,
                                                    outSample );
}

frantic::graphics::boundbox3f cached_trilerp_vector3f_staggered::get_sample_range(
    const const_rle_channel_accessor<frantic::graphics::vector3f>& dataAccessor,
    const frantic::graphics::vector3f& voxelCoordinate ) {
    const vector3 voxelCoordFloored = vector3::from_floor( voxelCoordinate );
    boundbox3 indicesBounds( voxelCoordFloored, voxelCoordFloored + vector3( 1 ) );
    if( voxelCoordinate.x - voxelCoordFloored.x < 0.5f )
        --indicesBounds.minimum().x;
    if( voxelCoordinate.y - voxelCoordFloored.y < 0.5f )
        --indicesBounds.minimum().y;
    if( voxelCoordinate.z - voxelCoordFloored.z < 0.5f )
        --indicesBounds.minimum().z;

    if( indicesBounds != m_cachedBounds ) {
        m_cachedBounds = indicesBounds;
        m_ris.fill_data_index_box( indicesBounds, m_dataIndices );
    }

    return get_trilerp_vector3f_staggered_sample_range_from_indices( dataAccessor, voxelCoordinate, m_dataIndices,
                                                                     indicesBounds );
}

//
//
// Internal "detail" functions
//
//

unsigned int detail::get_defined_voxel_count( const boost::int32_t* indices ) {
    unsigned int count = 0;
    for( int i = 0; i < 8; ++i )
        if( indices[i] >= 0 )
            ++count;
    return count;
}

void detail::create_2x2x2_indicies_array( const vector3f& voxelCoordinate, const boost::int32_t* const inIndices,
                                          const boundbox3& indicesVoxelBounds, boost::int32_t* outIndices,
                                          float* outDeltas ) {
    // this internal function tries to create a 2x2x2 data index array from an existing set (presumable larger set) of
    // indices. it also computes the deltas for interpolations.

    // throw an error if the data indices bound box is too small
    size3 size = indicesVoxelBounds.size();
    if( size.xsize() < 2 || size.ysize() < 2 || size.zsize() < 2 )
        throw std::runtime_error( "create_2x2x2_indicies_array: The voxel coordinate " + voxelCoordinate.str() +
                                  " can not be interpolated because the the provided dataIndicies boundbox " +
                                  indicesVoxelBounds.str() + " is too small." );

    // get the x,y,z deltas and voxelbounds minimum
    vector3f offsetVoxel = voxelCoordinate - vector3f( 0.5f );
    vector3 boxMin = vector3::from_floor( offsetVoxel );
    outDeltas[0] = offsetVoxel.x - boxMin.x;
    outDeltas[1] = offsetVoxel.y - boxMin.y;
    outDeltas[2] = offsetVoxel.z - boxMin.z;

    // check to see if voxelCoordinate is in indicesVoxelBounds. throw an error if it's not.
    // the original code also included a 1e-4 of a voxel tolerance with the justification that it should
    // be enough to capture and floating remains from the offset subtraction (although this seems large).
    for( int i = 0; i < 3; ++i ) {
        if( boxMin[i] < indicesVoxelBounds.minimum()[i] ) {
            if( indicesVoxelBounds.minimum()[i] - offsetVoxel[i] > 1e-4 ) {
                throw std::runtime_error(
                    "create_2x2x2_indicies_array: The voxel coordinate " + voxelCoordinate.str() +
                    " can not be interpolated because the interpolation samples needed are not included "
                    "within the provided dataIndicies boundbox " +
                    indicesVoxelBounds.str() + "." );
            } else {
                outDeltas[i] = 0.0f; // within tolerance
                ++boxMin[i];
            }
        }
        if( boxMin[i] + 1 > indicesVoxelBounds.maximum()[i] ) {
            if( offsetVoxel[i] - indicesVoxelBounds.maximum()[i] > 1e-4 ) {
                throw std::runtime_error(
                    "create_2x2x2_indicies_array: The voxel coordinate " + voxelCoordinate.str() +
                    " can not be interpolated because the interpolation samples needed are not included "
                    "within the provided dataIndicies boundbox " +
                    indicesVoxelBounds.str() + "." );
            } else {
                outDeltas[i] = 1.0f; // within tolerance
                --boxMin[i];
            }
        }
    }

    // build up the dataIndices 2x2x2 array
    vector3 origin = boxMin - indicesVoxelBounds.minimum();
    for( int z = origin.z, i = 0; z < origin.z + 2; ++z ) {
        for( int y = origin.y; y < origin.y + 2; ++y ) {
            for( int x = origin.x; x < origin.x + 2; ++x ) {
                outIndices[i] = inIndices[x + ( y + z * size.zsize() ) * size.ysize()];
                ++i;
            }
        }
    }
}

void detail::trilerp_float_from_2x2x2_indices( const float* dataAccessor, const boost::int32_t* indices, float* deltas,
                                               float* outTrilerp, vector3f* outGradient ) {
    // this internal function requires that the "indices" array contains no negative values. it is provided because it
    // contains less logic and will be fast. if there are negative values in the "indices" array, call the
    // trilerp_float_from_2x2x2_indices_with_undefined.

    float dx = deltas[0];
    float dy = deltas[1];
    float dz = deltas[2];

    // compute trilerp if requested
    if( outTrilerp ) {
        *outTrilerp =
            ( 1 - dz ) * ( ( 1 - dy ) * ( ( 1 - dx ) * dataAccessor[indices[0]] + dx * dataAccessor[indices[1]] ) +
                           dy * ( ( 1 - dx ) * dataAccessor[indices[2]] + dx * dataAccessor[indices[3]] ) ) +
            dz * ( ( 1 - dy ) * ( ( 1 - dx ) * dataAccessor[indices[4]] + dx * dataAccessor[indices[5]] ) +
                   dy * ( ( 1 - dx ) * dataAccessor[indices[6]] + dx * dataAccessor[indices[7]] ) );
    }

    // compute gradient if requested
    if( outGradient ) {
        outGradient->x = ( 1 - dz ) * ( ( 1 - dy ) * ( dataAccessor[indices[1]] - dataAccessor[indices[0]] ) +
                                        dy * ( dataAccessor[indices[3]] - dataAccessor[indices[2]] ) ) +
                         dz * ( ( 1 - dy ) * ( dataAccessor[indices[5]] - dataAccessor[indices[4]] ) +
                                dy * ( dataAccessor[indices[7]] - dataAccessor[indices[6]] ) );
        outGradient->y = ( 1 - dz ) * ( ( ( 1 - dx ) * dataAccessor[indices[2]] + dx * dataAccessor[indices[3]] ) -
                                        ( ( 1 - dx ) * dataAccessor[indices[0]] + dx * dataAccessor[indices[1]] ) ) +
                         dz * ( ( ( 1 - dx ) * dataAccessor[indices[6]] + dx * dataAccessor[indices[7]] ) -
                                ( ( 1 - dx ) * dataAccessor[indices[4]] + dx * dataAccessor[indices[5]] ) );
        outGradient->z = ( ( 1 - dy ) * ( ( 1 - dx ) * dataAccessor[indices[4]] + dx * dataAccessor[indices[5]] ) +
                           dy * ( ( 1 - dx ) * dataAccessor[indices[6]] + dx * dataAccessor[indices[7]] ) ) -
                         ( ( 1 - dy ) * ( ( 1 - dx ) * dataAccessor[indices[0]] + dx * dataAccessor[indices[1]] ) +
                           dy * ( ( 1 - dx ) * dataAccessor[indices[2]] + dx * dataAccessor[indices[3]] ) );
    }
}

void detail::trilerp_float_from_2x2x2_indices_with_undefined( const float* dataAccessor, const boost::int32_t* indices,
                                                              float* deltas, float* outTrilerp,
                                                              vector3f* outGradient ) {
    // this internal function is for trilerp'ing a value when one or more of the sample voxels are undefined.
    // it contains more logic than the "fully defined" version. it will ensure all the weighting of the defined voxel
    // values sum to 1.

    // compute the trilerp'd float value if requested
    if( outTrilerp != NULL ) {

        // compute the weights for each data index
        float weights[8];
        get_trilerp_weights( deltas, weights );
        normalize_trilerp_weights_for_undefined_voxels( indices, weights );

        // trilerp the output value by summing the value/weight pairs
        *outTrilerp = 0.0f;
        for( int i = 0; i < 8; ++i ) {
            boost::int32_t index = indices[i];
            if( index >= 0 )
                *outTrilerp += dataAccessor[index] * weights[i];
        }
    }

    // compute the float value's gradient if requested
    if( outGradient != NULL ) {

        float x = deltas[0], invx = 1.0f - deltas[0];
        float y = deltas[1], invy = 1.0f - deltas[1];
        float z = deltas[2], invz = 1.0f - deltas[2];

        // each gradient (in the x,y,z directions) has four samples we bi-linearly interpolate from.
        vector3f grad[4] = { vector3f( 0.0f ), vector3f( 0.0f ), vector3f( 0.0f ), vector3f( 0.0f ) };
        vector3f gradWeight[4] = { vector3f( 0.0f ), vector3f( 0.0f ), vector3f( 0.0f ), vector3f( 0.0f ) };
        vector3 gradDefinedCount( 0, 0, 0 );

        // get the four gradient samples and weights along the x-axis direction
        if( indices[1] >= 0 && indices[0] >= 0 ) {
            grad[0].x = dataAccessor[indices[1]] - dataAccessor[indices[0]];
            gradWeight[0].x = invz * invy;
            ++gradDefinedCount.x;
        }
        if( indices[3] >= 0 && indices[2] >= 0 ) {
            grad[1].x = dataAccessor[indices[3]] - dataAccessor[indices[2]];
            gradWeight[1].x = invz * y;
            ++gradDefinedCount.x;
        }
        if( indices[5] >= 0 && indices[4] >= 0 ) {
            grad[2].x = dataAccessor[indices[5]] - dataAccessor[indices[4]];
            gradWeight[2].x = z * invy;
            ++gradDefinedCount.x;
        }
        if( indices[7] >= 0 && indices[6] >= 0 ) {
            grad[3].x = dataAccessor[indices[7]] - dataAccessor[indices[6]];
            gradWeight[3].x = z * y;
            ++gradDefinedCount.x;
        }

        // get the four gradient samples and weights along the y-axis direction
        if( indices[2] >= 0 && indices[0] >= 0 ) {
            grad[0].y = dataAccessor[indices[2]] - dataAccessor[indices[0]];
            gradWeight[0].y = invz * invx;
            ++gradDefinedCount.y;
        }
        if( indices[3] >= 0 && indices[1] >= 0 ) {
            grad[1].y = dataAccessor[indices[3]] - dataAccessor[indices[1]];
            gradWeight[1].y = invz * x;
            ++gradDefinedCount.y;
        }
        if( indices[6] >= 0 && indices[4] >= 0 ) {
            grad[2].y = dataAccessor[indices[6]] - dataAccessor[indices[4]];
            gradWeight[2].y = z * invx;
            ++gradDefinedCount.y;
        }
        if( indices[7] >= 0 && indices[5] >= 0 ) {
            grad[3].y = dataAccessor[indices[7]] - dataAccessor[indices[5]];
            gradWeight[3].y = z * x;
            ++gradDefinedCount.y;
        }

        // get the four gradient samples and weights along the z-axis direction
        if( indices[4] >= 0 && indices[0] >= 0 ) {
            grad[0].z = dataAccessor[indices[4]] - dataAccessor[indices[0]];
            gradWeight[0].z = invy * invx;
            ++gradDefinedCount.z;
        }
        if( indices[5] >= 0 && indices[1] >= 0 ) {
            grad[1].z = dataAccessor[indices[5]] - dataAccessor[indices[1]];
            gradWeight[1].z = invy * x;
            ++gradDefinedCount.z;
        }
        if( indices[6] >= 0 && indices[2] >= 0 ) {
            grad[2].z = dataAccessor[indices[6]] - dataAccessor[indices[2]];
            gradWeight[2].z = y * invx;
            ++gradDefinedCount.z;
        }
        if( indices[7] >= 0 && indices[3] >= 0 ) {
            grad[3].z = dataAccessor[indices[7]] - dataAccessor[indices[3]];
            gradWeight[3].z = y * x;
            ++gradDefinedCount.z;
        }

        // set the final gradient value for each component direction x,y,z
        for( int comp = 0; comp < 3; ++comp ) {

            if( gradDefinedCount[comp] == 0 ) {
                // if there are no defined gradient pair samples for this component, we can just set the gradient to
                // zero.
                ( *outGradient )[comp] = 0.0f;
            } else {

                // if there are one or more undefined gradient pairs for this component, normalize the defined gradient
                // sample weights to sum to 1.0f.
                if( gradDefinedCount[comp] < 4 ) {
                    float sum = gradWeight[0][comp] + gradWeight[1][comp] + gradWeight[2][comp] + gradWeight[3][comp];
                    if( sum > 1e-20f ) {
                        // divide by the sum, which ensures that the sum of the weights will equal 1.0f.
                        gradWeight[0][comp] /= sum;
                        gradWeight[1][comp] /= sum;
                        gradWeight[2][comp] /= sum;
                        gradWeight[3][comp] /= sum;
                    } else {
                        // if the sum is infinitly small, dividing by an infinity small number will produce floating
                        // point error. the following will treat each defined voxel as equal weight in the case where
                        // none of them really contribute (does this make sense?). note that it doesn't matter if an
                        // undefined voxel gets a non-zero weight, because their 'grad' values are zero anyway.
                        float weight = 1.0f / gradDefinedCount[comp];
                        gradWeight[0][comp] = weight;
                        gradWeight[1][comp] = weight;
                        gradWeight[2][comp] = weight;
                        gradWeight[3][comp] = weight;
                    }
                }

                // summing the value/weight pairs to get the final component gradient value
                ( *outGradient )[comp] = grad[0][comp] * gradWeight[0][comp] + grad[1][comp] * gradWeight[1][comp] +
                                         grad[2][comp] * gradWeight[2][comp] + grad[3][comp] * gradWeight[3][comp];
            }
        }
    }
}

void detail::trilerp_vector3f_from_2x2x2_indices( const_rle_channel_accessor<vector3f> dataAccessor,
                                                  const boost::int32_t* indices, float* deltas, vector3f& outTrilerp ) {
    // this internal function requires that the "indices" array contains no negative values. it is provided because it
    // contains less logic and will be fast. if there are negative values in the "indices" array, call the
    // trilerp_float_from_2x2x2_indices_with_undefined.

    float dx = deltas[0];
    float dy = deltas[1];
    float dz = deltas[2];
    outTrilerp = ( 1 - dz ) * ( ( 1 - dy ) * ( ( 1 - dx ) * dataAccessor[indices[0]] + dx * dataAccessor[indices[1]] ) +
                                dy * ( ( 1 - dx ) * dataAccessor[indices[2]] + dx * dataAccessor[indices[3]] ) ) +
                 dz * ( ( 1 - dy ) * ( ( 1 - dx ) * dataAccessor[indices[4]] + dx * dataAccessor[indices[5]] ) +
                        dy * ( ( 1 - dx ) * dataAccessor[indices[6]] + dx * dataAccessor[indices[7]] ) );
}

void detail::trilerp_vector3f_from_2x2x2_indices_with_undefined( const_rle_channel_accessor<vector3f> dataAccessor,
                                                                 const boost::int32_t* indices, float* deltas,
                                                                 vector3f& outTrilerp ) {
    // this internal function is for trilerp'ing a value when one or more of the sample voxels are undefined.
    // it contains more logic than the "fully defined" version. it will ensure all the weighting of the defined voxel
    // values sum to 1.

    // compute the weights for each data index
    float weights[8];
    get_trilerp_weights( deltas, weights );
    normalize_trilerp_weights_for_undefined_voxels( indices, weights );

    // trilerp the output value by summing the value/weight pairs
    outTrilerp.set( 0.0f );
    for( int i = 0; i < 8; ++i ) {
        boost::int32_t index = indices[i];
        if( index >= 0 )
            outTrilerp += dataAccessor[index] * weights[i];
    }
}

void detail::trilerp_general_channel_from_2x2x2_indices( const_rle_channel_general_accessor dataAccessor,
                                                         const boost::int32_t* indices, float* deltas,
                                                         char* outTrilerp ) {
    // compute the weights for each data index
    float weights[8];
    get_trilerp_weights( deltas, weights );
    normalize_trilerp_weights_for_undefined_voxels( indices, weights );

    // the get_trilinear_interpolated functions currently does not currently re-normalize weights for undefined voxels
    // if it changes, be sure to remove the above normalization.
    dataAccessor.get_trilinear_interpolated( weights, indices, outTrilerp );
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
