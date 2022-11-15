// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <frantic/logging/global_progress_logger.hpp>
#include <frantic/volumetrics/levelset/level_set_marching_extrapolation.hpp>
#include <stdexcept>

using namespace std;
using namespace boost;
using namespace frantic::volumetrics;

namespace frantic {
namespace volumetrics {
namespace levelset {

namespace detail {

/**
 * This function implements extrapolation of level set named channels.
 *
 * Note that this function doesn't do any heap allocation, all variables are stack-allocated.  This is done
 * for performance purposes.
 *
 * This is the algorithm "marching_extrapolation" given in pseudo-code there.
 */
void marching_extrapolation_recursive( int dataIndex, const ris_adjacency& adj, float voxelLength,
                                       unsigned char* populatedChannel, const float* signedDistanceChannel,
                                       const vector<marching_extrapolation_channel>& extrapChannels,
                                       vector<int>& outDataIndexesToProcess,
                                       vector<std::pair<int, float>>& dataIndicesSorted,
                                       extrapolation_debug_info& db ) {
    // Flag this voxel as "Unpopulated, but is already on the stack to be processed"
    populatedChannel[dataIndex] = 4;
    const ris_adj_entry& rae = adj[dataIndex];
    float signedDistance = signedDistanceChannel[dataIndex];

    //	cout << debugIndent << "Started processing data index " << dataIndex << ", with signed distance " <<
    // signedDistance << endl;

    bool ambiguousExtrapolation = false;

    /////////
    // CATEGORIZE INPUTS AND OUTPUTS
    /////////

    char inputFlag, outputFlag;
    if( signedDistance >= 0 ) {
        inputFlag = 1;
        outputFlag = 2;
    } else {
        inputFlag = 2;
        outputFlag = 1;
    }

    // adjCmpFlags and adjSignedDistance specify information about the 6 neighbor voxels
    // In adjCmpFlags, 0 means undefined, 1 means less than signedDistance, 2 means greater than signedDistance
    char adjCmpFlags[6];
    float adjSignedDistance[6];

    for( int i = 0; i != 6; ++i ) {
        if( rae[i] < 0 ) {
            adjCmpFlags[i] = 0;
        } else {
            adjSignedDistance[i] = signedDistanceChannel[rae[i]];
            if( adjSignedDistance[i] < signedDistance )
                adjCmpFlags[i] = 1;
            else if( adjSignedDistance[i] > signedDistance )
                adjCmpFlags[i] = 2;
            else {
                // The weight for this as an input would be zero anyway, so just consider it as an output.
                adjCmpFlags[i] = outputFlag;
            }
        }
    }

    /////////
    // PROCESS INPUTS
    /////////

    // Go through all the inputs, and recursively populate them if necessary
    int inputCount = 0;
    for( int i = 0; i != 6; ++i ) {
        if( adjCmpFlags[i] == inputFlag ) {
            unsigned char voxelPopulated = populatedChannel[rae[i]];
            if( voxelPopulated != 1 ) {
                if( voxelPopulated == 0 ) {
                    if( ( inputFlag == 1 && signedDistanceChannel[rae[i]] < 0 ) ||
                        ( inputFlag == 2 && signedDistanceChannel[rae[i]] >= 0 ) ) {
                        // If the signed distances at the current and adjacent voxels have opposite signs,
                        // there is a dependency loop between data at these voxels.  In that case, we remove
                        // the adjacent voxel as an input, and indicate an ambiguity.  Setting voxelPopulated to 3
                        // for the if statement below.
                        voxelPopulated = 3;

                        // if( logging::is_logging_debug() ) {
                        //	std::cout <<  "InputFlag: " << (int)inputFlag;
                        //	std::cout <<  "\tphi: " << signedDistanceChannel[rae[i]] << endl;

                        //	db.reset();
                        //	db.create_debug_mesh();

                        //	throw std::runtime_error("First ambigious voxel input found");
                        //}

                    } else {
                        // Otherwise recursively fill this voxel with data
                        marching_extrapolation_recursive( rae[i], adj, voxelLength, populatedChannel,
                                                          signedDistanceChannel, extrapChannels,
                                                          outDataIndexesToProcess, dataIndicesSorted, db );

                        // db.create_debug_mesh();

                        // Get the populated channel value again
                        voxelPopulated = populatedChannel[rae[i]];
                    }
                }
                // The recursive call above may have changed voxelPopulated, so the following tests can't go in an else
                // statement.
                if( voxelPopulated == 2 ) {
                    // A 'Populated' flag of 2 indicates that it is unpopulated, but should not be used as an
                    // extrapolation target.
                    adjCmpFlags[i] = 0;
                } else if( voxelPopulated == 3 ) {
                    // A 'Populated' flag of 3 indicates that it is unpopulated, but should not be populated
                    // due to lack of data or due to ambiguity.
                    adjCmpFlags[i] = 0;
                    // ambiguousExtrapolation = true;
                } else if( voxelPopulated != 1 ) {
                    throw runtime_error( "marching_extrapolation_recursive() - Whoops, expected 1, got " +
                                         lexical_cast<string>( (int)voxelPopulated ) + " for populated value." );
                }
            }

            // If this adjacent voxel is still flagged as an input, then increment the input count;
            if( adjCmpFlags[i] == inputFlag )
                ++inputCount;
        }
    }

    /////////
    // PROCESS CURRENT VOXEL
    /////////

    if( inputCount > 1 ) {
        // Collect all the input data indexes and weights into arrays.
        float totalWeight = 0; // The sum of all the extrapolation weights
        float weights[6];      // The extrapolation weights
        int dataIndexes[6];    // The data indexes for the extrapolation data
        char* dataPointers[6]; // Pointers to the extrapolation data
        int neighborIndex[6];  // Neighbor index of the extrapolation data
        inputCount = 0;
        for( int i = 0; i != 6; i += 2 ) {
            if( adjCmpFlags[i] == inputFlag && adjCmpFlags[i + 1] == inputFlag ) {
                int index0 = rae[i], index1 = rae[i + 1];
                float weight0 = fabsf( signedDistanceChannel[index0] - signedDistance ),
                      weight1 = fabsf( signedDistanceChannel[index1] - signedDistance );
                float scaleRatio = max( weight0, weight1 ) / ( weight0 + weight1 );
                weight0 *= scaleRatio;
                weight1 *= scaleRatio;
                dataIndexes[inputCount] = index0;
                weights[inputCount] = weight0;
                neighborIndex[inputCount] = i;
                ++inputCount;
                dataIndexes[inputCount] = index1;
                weights[inputCount] = weight1;
                neighborIndex[inputCount] = i + 1;
                ++inputCount;
                totalWeight += weight0;
                totalWeight += weight1;
            } else if( adjCmpFlags[i] == inputFlag ) {
                int index = rae[i];
                float weight = fabsf( signedDistanceChannel[index] - signedDistance );
                dataIndexes[inputCount] = index;
                weights[inputCount] = weight;
                neighborIndex[inputCount] = i;
                ++inputCount;
                totalWeight += weight;
            } else if( adjCmpFlags[i + 1] == inputFlag ) {
                int index = rae[i + 1];
                float weight = fabsf( signedDistanceChannel[index] - signedDistance );
                dataIndexes[inputCount] = index;
                weights[inputCount] = weight;
                neighborIndex[inputCount] = i + 1;
                ++inputCount;
                totalWeight += weight;
            }
        }
        if( totalWeight > 0 ) {
            // Normalize the weights
            for( int i = 0; i != inputCount; ++i )
                weights[i] /= totalWeight;
        } else {
            // If no weight at all was collected, take the average of all the inputs.
            float weight = 1.f / inputCount;
            for( int i = 0; i != inputCount; ++i )
                weights[i] = weight;
        }
        // Go through all the channels and use the linear combination function to generate the extrapolated value.
        for( vector<marching_extrapolation_channel>::const_iterator ch = extrapChannels.begin(),
                                                                    che = extrapChannels.end();
             ch != che; ++ch ) {
            if( ch->gradientToMatch == 0 ) { // Equivalent to matching a zero gradient
                // First fill in the data pointers
                for( int i = 0; i != inputCount; ++i )
                    dataPointers[i] = ch->data + dataIndexes[i] * ch->primitiveSize;
                // Call the linear combination function to do the extrapolation to this voxel.
                ch->weightedSumFn( weights, dataPointers, inputCount, ch->arity,
                                   ch->data + dataIndex * ch->primitiveSize );
            } else {
                // In this case, we have to extrapolate all the input points, and then interpolate the extrapolated
                // values
                char* tempData = new char[ch->primitiveSize * inputCount];
                float extrapWeights[2] = { 1.f, voxelLength };
                const char* extrapData[2];
                // Extrapolate and fill in the data pointers
                for( int i = 0; i != inputCount; ++i ) {
                    // Pointer to the data
                    extrapData[0] = ch->data + dataIndexes[i] * ch->primitiveSize;
                    // Pointer to the element of the data gradient which corresponds to the neighbor axis direction
                    extrapData[1] =
                        ch->gradientToMatch +
                        ( 3 * dataIndexes[i] + neighbor_index_axis( neighborIndex[i] ) ) * ch->primitiveSize;
                    // The multiplication factor for the gradient is +/-voxelLength
                    extrapWeights[1] =
                        is_neighbor_index_direction_positive( neighborIndex[i] ) ? -voxelLength : voxelLength;
                    // The data pointers point inside our temporary array now
                    dataPointers[i] = tempData + i * ch->primitiveSize;
                    ch->weightedSumFn( extrapWeights, extrapData, 2, ch->arity, dataPointers[i] );
                }

                // Call the linear combination function to do the extrapolation to this voxel.
                ch->weightedSumFn( weights, dataPointers, inputCount, ch->arity,
                                   ch->data + dataIndex * ch->primitiveSize );

                delete[] tempData;
            }
        }
    } else if( inputCount == 1 ) {
        // Find the single input, and copy its value.
        for( int i = 0; i != 6; ++i ) {
            if( adjCmpFlags[i] == inputFlag ) {
                // Go through all the channels and copy this voxel value for each of them.
                for( vector<marching_extrapolation_channel>::const_iterator ch = extrapChannels.begin(),
                                                                            che = extrapChannels.end();
                     ch != che; ++ch ) {
                    memcpy( ch->data + dataIndex * ch->primitiveSize, ch->data + rae[i] * ch->primitiveSize,
                            ch->primitiveSize );
                }
                break;
            }
        }
    } else {
        // If there were no inputs, then indicate that there were was a problem with this voxel.
        ambiguousExtrapolation = true;
    }

    if( ambiguousExtrapolation ) {
        // If there was any weirdness, flag this voxel as "Was an unpopulated extrapolation target, but
        // there was insufficient or ambiguous input, so a value could not be produced"
        populatedChannel[dataIndex] = 3;
    } else {
        // Otherwise flag this voxel as populated
        populatedChannel[dataIndex] = 1;
    }

    // if(logging::is_logging_debug() && ambiguousExtrapolation ){

    //	db.create_debug_mesh();

    ////	/*if ( db.counter > 2000 )
    ////		logging::set_logging_level(3);*/
    //}
    //

    //	cout << debugIndent << "Set Populated[" << dataIndex << "] to " << (int)populatedChannel[dataIndex] << endl;

    /////////
    // PROCESS OUTPUTS
    /////////

    for( int i = 0; i != 6; ++i ) {
        if( adjCmpFlags[i] == outputFlag ) {
            int index = rae[i];
            // No need to check that index is >= 0, because if it were, adjCmpFlags[i] would have been 0 instead of
            // outputFlag.
            if( populatedChannel[index] == 0 ) {
                outDataIndexesToProcess.push_back( index );
                dataIndicesSorted.push_back( make_pair( index, signedDistanceChannel[index] ) );
            }
        }
    }
}

} // namespace detail

bool marching_predicate( const std::pair<int, float>& left, const std::pair<int, float>& right ) {
    return fabsf( left.second ) < fabsf( right.second );
}

void marching_extrapolation( const ris_adjacency& adj, float voxelLength, unsigned char* populatedChannel,
                             const float* signedDistanceChannel,
                             const std::vector<marching_extrapolation_channel>& extrapChannels,
                             extrapolation_debug_info& db ) {
    vector<int> dataIndexesToProcess;

    vector<std::pair<int, float>> dataIndicesSorted;

    // if( logging::is_logging_debug() )
    //	db.create_debug_mesh();

    for( size_t i = 0, ie = adj.data_size(); i != ie; ++i ) {
        float phi = signedDistanceChannel[i];
        switch( populatedChannel[i] ) {
        case 0: { // Unpopulated extrapolation target
            const ris_adj_entry& rae = adj[i];
            // If there's an adjacent populated voxel, then start marching from this voxel
            if( ( rae.x_pos >= 0 && populatedChannel[rae.x_pos] == 1 ) ||
                ( rae.y_pos >= 0 && populatedChannel[rae.y_pos] == 1 ) ||
                ( rae.z_pos >= 0 && populatedChannel[rae.z_pos] == 1 ) ) {
                dataIndexesToProcess.push_back( static_cast<int>( i ) );
                dataIndicesSorted.push_back( make_pair( (int)i, phi ) );
            }
            break;
        }
        case 1: { // Populated extrapolation source
            const ris_adj_entry& rae = adj[i];
            // If there's an adjacent unpopulated voxel, then start marching from that voxel
            if( rae.x_pos >= 0 && populatedChannel[rae.x_pos] == 0 ) {
                dataIndexesToProcess.push_back( rae.x_pos );

                dataIndicesSorted.push_back( pair<int, float>( rae.x_pos, signedDistanceChannel[rae.x_pos] ) );
            }
            if( rae.y_pos >= 0 && populatedChannel[rae.y_pos] == 0 ) {
                dataIndexesToProcess.push_back( rae.y_pos );
                dataIndicesSorted.push_back( pair<int, float>( rae.y_pos, signedDistanceChannel[rae.y_pos] ) );
            }

            if( rae.z_pos >= 0 && populatedChannel[rae.z_pos] == 0 ) {
                dataIndexesToProcess.push_back( rae.z_pos );
                dataIndicesSorted.push_back( pair<int, float>( rae.z_pos, signedDistanceChannel[rae.z_pos] ) );
            }
            break;
        }
        }
    }

    std::sort( dataIndicesSorted.begin(), dataIndicesSorted.end(), marching_predicate );

    while( !dataIndicesSorted.empty() ) {

        // size_t selection = rand()%dataIndexesToProcess.size();
        // std::swap( dataIndexesToProcess[selection], dataIndexesToProcess[dataIndexesToProcess.size()-1]);

        // int index = dataIndexesToProcess.back();
        //	dataIndexesToProcess.pop_back();

        int index = dataIndicesSorted[0].first;

        std::swap( dataIndicesSorted[0], dataIndicesSorted[dataIndicesSorted.size() - 1] );
        dataIndicesSorted.pop_back();

        if( populatedChannel[index] == 0 ) {
            //				cout << "Starting marching extrapolation on data index " << index << endl;
            detail::marching_extrapolation_recursive( index, adj, voxelLength, populatedChannel, signedDistanceChannel,
                                                      extrapChannels, dataIndexesToProcess, dataIndicesSorted, db );
            // if( db.channelAccessor->data(0) != (const char*)populatedChannel )
            //	throw std::runtime_error("uh oh");

            // if( logging::is_logging_debug() )
            //	db.create_debug_mesh();
        }
    }
    //}
}

/////////////////////////////////////////
// STAGGERED FIELD MARCHING EXTRAPOLATION
/////////////////////////////////////////

enum staggered_flags {
    STAGGERED_FLAG_UNTOUCHED = 0,
    STAGGERED_FLAG_VALUE_SOURCE = 1,
    STAGGERED_FLAG_VALUE_TARGET = 2,
    STAGGERED_FLAG_PROCESSING = 3,
    STAGGERED_FLAG_IGNORE = 4
};

namespace detail {

/**
 * This function sets the appropriate flag in the flagChannel, and returns it.
 */
unsigned char set_staggered_voxel_flag( unsigned char* flagChannel, int i,
                                        const boost::uint8_t* staggeredPopulatedChannel, boost::uint8_t staggeredMask )
//	const std::string& indent = "" )
{
    unsigned char flag = flagChannel[i];
    if( flag == STAGGERED_FLAG_UNTOUCHED ) {
        if( staggeredPopulatedChannel[i] & staggeredMask )
            flag = STAGGERED_FLAG_VALUE_SOURCE;
        else
            flag = STAGGERED_FLAG_VALUE_TARGET;
        flagChannel[i] = flag;
    }
    return flag;
}

static int ambiguousCount = 0;

/**
 * This function recursively extrapolates to the voxel specified by the data index.
 *
 * @param  staggeredExtrapChannel  This is the staggeredExtrapChannel from the main function, but
 *                                 with +0, +1 or +2 added depending on which of U, V, or W faces
 *                                 are being processed.
 */
void staggered_field_marching_extrapolation_recursive( int dataIndex, const ris_adjacency& adj,
                                                       unsigned char* flagChannel, const float* signedDistanceChannel,
                                                       const boost::uint8_t* staggeredPopulatedChannel,
                                                       boost::uint8_t staggeredMask, bool extrapToPositive,
                                                       boost::int32_t* dataIndexMapChannel,
                                                       float* staggeredExtrapChannel, int raeNeighborIndex,
                                                       vector<boost::int32_t>& outDataIndexesToProcess )
//	const string& indent = " " )
{
    //	cout << indent << "R rls index: " << dataIndex << "\n";

    // Flag this voxel as "Not finished, but is already on the stack to be processed"
    flagChannel[dataIndex] = STAGGERED_FLAG_PROCESSING;
    const ris_adj_entry& rae = adj[dataIndex];
    int neighborDataIndex = rae[raeNeighborIndex];

    // Get the signed distance for this voxel.  Note that we extrapolate to voxels even when there isn't an
    // adjacent neighbor to get the interpolated signed distance.
    float signedDistance;
    if( neighborDataIndex >= 0 )
        signedDistance = 0.5f * ( signedDistanceChannel[dataIndex] + signedDistanceChannel[neighborDataIndex] );
    else
        signedDistance = signedDistanceChannel[dataIndex];

    //	cout << debugIndent << "Started processing data index " << dataIndex << ", with signed distance " <<
    // signedDistance << endl;

    bool ambiguousExtrapolation = false;

    /////////
    // CATEGORIZE INPUTS AND OUTPUTS
    /////////

    char inputFlag, outputFlag;
    if( extrapToPositive ) {
        inputFlag = 1;
        outputFlag = 2;
    } else {
        inputFlag = 2;
        outputFlag = 1;
    }

    // adjCmpFlags and adjSignedDistance specify information about the 6 neighbor voxels
    // In adjCmpFlags, 0 means undefined, 1 means less than signedDistance, 2 means greater than signedDistance
    char adjCmpFlags[6];
    float adjSignedDistance[6];
    for( int i = 0; i != 6; ++i ) {
        // The first data index is the neighbor voxel in direction i
        int dataIndex0 = rae[i];
        // Exclude it if there is no distance data, or no staggered field data
        if( dataIndex0 < 0 || dataIndexMapChannel[dataIndex0] < 0 ) {
            adjCmpFlags[i] = 0;
        } else {
            // The second data index is the neighbor of the first data index in the direction of the U, V, or
            // W face being processed, as indicated by raeNeighborIndex
            int dataIndex1 = adj[dataIndex0][raeNeighborIndex];
            if( dataIndex1 < 0 ) {
                // If there is no signed distance at the neighbor, we still process it but with just this signed
                // distance
                adjSignedDistance[i] = signedDistanceChannel[dataIndex0];
            } else {
                // To find the distance value at the face, we take the average of the two closest voxel centers.
                adjSignedDistance[i] = 0.5f * ( signedDistanceChannel[dataIndex0] + signedDistanceChannel[dataIndex1] );
            }
            if( adjSignedDistance[i] < signedDistance )
                adjCmpFlags[i] = 1;
            else if( adjSignedDistance[i] > signedDistance )
                adjCmpFlags[i] = 2;
            else {
                // The weight for this as an input would be zero anyway, so just consider it as an output.
                adjCmpFlags[i] = outputFlag;
            }
        }
    }

    /////////
    // PROCESS INPUTS
    /////////

    // Go through all the inputs, and recursively populate them if necessary
    int inputCount = 0;
    for( int i = 0; i != 6; ++i ) {
        if( adjCmpFlags[i] == inputFlag ) {
            unsigned char voxelPopulated =
                set_staggered_voxel_flag( flagChannel, rae[i], staggeredPopulatedChannel, staggeredMask );
            //			cout << indent << "Input neighbor rls index: " << rae[i] << ", flag: " << (int)voxelPopulated <<
            //"\n";
            if( voxelPopulated != STAGGERED_FLAG_VALUE_SOURCE ) {
                if( voxelPopulated == STAGGERED_FLAG_VALUE_TARGET ) {
                    // Recursively fill this voxel with data
                    staggered_field_marching_extrapolation_recursive(
                        rae[i], adj, flagChannel, signedDistanceChannel, staggeredPopulatedChannel, staggeredMask,
                        extrapToPositive, dataIndexMapChannel, staggeredExtrapChannel, raeNeighborIndex,
                        outDataIndexesToProcess ); //, indent + " " );
                    // Get the populated channel value again
                    voxelPopulated = flagChannel[rae[i]];
                }
                // The recursive call above may have changed voxelPopulated, so the following tests can't go in an else
                // statement.
                if( voxelPopulated != STAGGERED_FLAG_VALUE_SOURCE ) {
                    // If the recursive call didn't fill the value, then exclude it
                    adjCmpFlags[i] = 0;
                }
            }

            // If this adjacent voxel is still flagged as an input, then increment the input count;
            if( adjCmpFlags[i] == inputFlag )
                ++inputCount;
        }
    }

    //	cout << indent << "Input count " << inputCount << "\n";

    /////////
    // PROCESS CURRENT VOXEL
    /////////

    if( inputCount > 1 ) {
        // Collect all the input data indexes and weights into arrays.
        float totalWeight = 0;
        float weights[6];
        int dataIndexes[6];
        inputCount = 0;
        // Go through i == 0, 2, 4, and process +/- directions along each axis together
        for( int i = 0; i != 6; i += 2 ) {
            // If there are inputs in both the + and - directions along this axis
            if( adjCmpFlags[i] == inputFlag && adjCmpFlags[i + 1] == inputFlag ) {
                // Get the indices of the two neighbors
                int index0 = rae[i], index1 = rae[i + 1];
                // Compute the weights, and adjust them so that this axis doesn't get over-represented because
                // of having two inputs instead of one
                float weight0 = fabsf( signedDistanceChannel[index0] - signedDistance ),
                      weight1 = fabsf( signedDistanceChannel[index1] - signedDistance );
                float scaleRatio = max( weight0, weight1 ) / ( weight0 + weight1 );
                weight0 *= scaleRatio;
                weight1 *= scaleRatio;
                // Add the weights and data indices to the inputs
                dataIndexes[inputCount] = index0;
                weights[inputCount] = weight0;
                ++inputCount;
                dataIndexes[inputCount] = index1;
                weights[inputCount] = weight1;
                ++inputCount;
                totalWeight += weight0;
                totalWeight += weight1;
            }
            // Otherwise only one or none of the directions will be an input
            else if( adjCmpFlags[i] == inputFlag ) {
                // Get the index of the neighbor
                int index = rae[i];
                // Compute the weight to use
                float weight = fabsf( signedDistanceChannel[index] - signedDistance );
                // Add the weight and data index to the inputs
                dataIndexes[inputCount] = index;
                weights[inputCount] = weight;
                ++inputCount;
                totalWeight += weight;
            } else if( adjCmpFlags[i + 1] == inputFlag ) {
                // Get the index of the neighbor
                int index = rae[i + 1];
                // Compute the weight to use
                float weight = fabsf( signedDistanceChannel[index] - signedDistance );
                // Add the weight and data index to the inputs
                dataIndexes[inputCount] = index;
                weights[inputCount] = weight;
                ++inputCount;
                totalWeight += weight;
            }
        }
        // Normalize the weights if there was any actual weight collected
        if( totalWeight > 0 ) {
            for( int i = 0; i != inputCount; ++i )
                weights[i] /= totalWeight;
        } else {
            // If no weight at all was collected, take the average of all the inputs.
            float weight = 1.f / inputCount;
            for( int i = 0; i != inputCount; ++i )
                weights[i] = weight;
        }
        // Compute the extrapolated value, and save it
        float extrapolatedValue = 0;
        for( int i = 0; i != inputCount; ++i )
            extrapolatedValue += weights[i] * staggeredExtrapChannel[3 * dataIndexMapChannel[dataIndexes[i]]];
        staggeredExtrapChannel[3 * dataIndexMapChannel[dataIndex]] = extrapolatedValue;
    } else if( inputCount == 1 ) {
        // Find the single input, and copy its value.
        for( int i = 0; i != 6; ++i ) {
            if( adjCmpFlags[i] == inputFlag ) {
                staggeredExtrapChannel[3 * dataIndexMapChannel[dataIndex]] =
                    staggeredExtrapChannel[3 * dataIndexMapChannel[rae[i]]];
                break;
            }
        }
    } else {
        // If there were no inputs, then indicate that there were was a problem with this voxel.
        ambiguousExtrapolation = true;
        // Set the value to 0
        staggeredExtrapChannel[3 * dataIndexMapChannel[dataIndex]] = 0;
    }

    if( ambiguousExtrapolation ) {
        ++ambiguousCount;
        // cout << "\t**** A rls index: " << dataIndex << " ambiguous, did not change its value! :(\n";
        //  If there was any weirdness, flag this voxel as "Was an unpopulated extrapolation target, but
        //  there was insufficient or ambiguous input, so a value could not be produced"
        flagChannel[dataIndex] = STAGGERED_FLAG_IGNORE;
    } else {
        //		cout << indent << "S rls index: " << dataIndex << " processed, set to source\n";
        // Otherwise flag this voxel as populated
        flagChannel[dataIndex] = STAGGERED_FLAG_VALUE_SOURCE;
    }

    //	cout << debugIndent << "Set Populated[" << dataIndex << "] to " << (int)populatedChannel[dataIndex] << endl;

    /////////
    // PROCESS OUTPUTS
    /////////

    /*if( logging::is_logging_debug() ) {
      std::cout << "outputing mesh: " << db.counter << std::endl;
      db.create_debug_mesh();
    }*/

    for( int i = 0; i != 6; ++i ) {
        if( adjCmpFlags[i] == outputFlag ) {
            int index = rae[i];

            // No need to check that index is >= 0, because if it were, adjCmpFlags[i] would have been 0 instead of
            // outputFlag.
            unsigned char flag =
                detail::set_staggered_voxel_flag( flagChannel, index, staggeredPopulatedChannel, staggeredMask );
            if( flag == STAGGERED_FLAG_VALUE_TARGET ) {
                //				cout << indent << "Output neighbor rls index: " << index << "\n";
                outDataIndexesToProcess.push_back( index );
            }
        }
    }
}

/**
 * This function does marching extrapolation on the U, V or W component of the staggered channel.
 *
 * @note Look at the main function definition to see what most of the other parameters mean.
 *
 * @param  staggeredExtrapChannel  This is the staggeredExtrapChannel from the main function, but
 *                                 with +0, +1 or +2 added depending on which of U, V, or W faces
 *                                 are being processed.
 * @param  raeNeighborIndex  This is the index within an rae, which will point to one of rae.x_neg,
 *                           rae.y_neg, or rae.z_neg.
 */
void staggered_field_marching_extrapolation_singleface( const ris_adjacency& adj, const float* signedDistanceChannel,
                                                        const boost::uint8_t* staggeredPopulatedChannel,
                                                        boost::uint8_t staggeredMask, bool extrapToPositive,
                                                        boost::int32_t* dataIndexMapChannel,
                                                        float* staggeredExtrapChannel, int raeNeighborIndex ) {
    // These data indexes refer to signedDistanceChannel
    vector<boost::int32_t> dataIndexesToProcess;

    // Create a flag channel (starting with all 0 values)
    vector<unsigned char> flagChannelVector( adj.data_size(), STAGGERED_FLAG_UNTOUCHED );
    unsigned char* flagChannel = &flagChannelVector[0];

    for( size_t i = 0, ie = adj.data_size(); i != ie; ++i ) {
        // Only process voxels which have a corresponding value in the staggered field
        if( dataIndexMapChannel[i] >= 0 ) {
            const ris_adj_entry& rae = adj[i];
            // Determine whether this value is a source, a target, or neither
            unsigned char flag =
                detail::set_staggered_voxel_flag( flagChannel, (int)i, staggeredPopulatedChannel, staggeredMask );
            //			cout << "M rls index: " << i << ", rvf index: " << dataIndexMapChannel[i] << ", flag: " << (int)flag
            //<<
            //"\n";
            switch( flag ) {
            case STAGGERED_FLAG_VALUE_TARGET: {
                bool crossesBoundary = false;
                if( rae.x_pos >= 0 ) {
                    flag = detail::set_staggered_voxel_flag( flagChannel, rae.x_pos, staggeredPopulatedChannel,
                                                             staggeredMask );
                    if( flag == STAGGERED_FLAG_VALUE_SOURCE )
                        crossesBoundary = true;
                }
                if( !crossesBoundary && rae.y_pos >= 0 ) {
                    flag = detail::set_staggered_voxel_flag( flagChannel, rae.y_pos, staggeredPopulatedChannel,
                                                             staggeredMask );
                    if( flag == STAGGERED_FLAG_VALUE_SOURCE )
                        crossesBoundary = true;
                }
                if( !crossesBoundary && rae.z_pos >= 0 ) {
                    flag = detail::set_staggered_voxel_flag( flagChannel, rae.z_pos, staggeredPopulatedChannel,
                                                             staggeredMask );
                    if( flag == STAGGERED_FLAG_VALUE_SOURCE )
                        crossesBoundary = true;
                }
                // If there's an adjacent source voxel, then start marching from this voxel
                if( crossesBoundary )
                    dataIndexesToProcess.push_back( static_cast<int>( i ) );
                break;
            }
            case STAGGERED_FLAG_VALUE_SOURCE: {
                // If there's an adjacent target voxel, then start marching from that voxel
                if( rae.x_pos >= 0 ) {
                    flag = detail::set_staggered_voxel_flag( flagChannel, rae.x_pos, staggeredPopulatedChannel,
                                                             staggeredMask );
                    if( flag == STAGGERED_FLAG_VALUE_TARGET && dataIndexMapChannel[rae.x_pos] >= 0 )
                        dataIndexesToProcess.push_back( rae.x_pos );
                }
                if( rae.y_pos >= 0 ) {
                    flag = detail::set_staggered_voxel_flag( flagChannel, rae.y_pos, staggeredPopulatedChannel,
                                                             staggeredMask );
                    if( flag == STAGGERED_FLAG_VALUE_TARGET && dataIndexMapChannel[rae.y_pos] >= 0 )
                        dataIndexesToProcess.push_back( rae.y_pos );
                }
                if( rae.z_pos >= 0 ) {
                    flag = detail::set_staggered_voxel_flag( flagChannel, rae.z_pos, staggeredPopulatedChannel,
                                                             staggeredMask );
                    if( flag == STAGGERED_FLAG_VALUE_TARGET && dataIndexMapChannel[rae.z_pos] >= 0 )
                        dataIndexesToProcess.push_back( rae.z_pos );
                }
                break;
            }
            }
        }
    }

    while( !dataIndexesToProcess.empty() ) {

        size_t selection = rand() % dataIndexesToProcess.size();
        std::swap( dataIndexesToProcess[selection], dataIndexesToProcess[dataIndexesToProcess.size() - 1] );

        int index = dataIndexesToProcess.back();
        dataIndexesToProcess.pop_back();

        if( flagChannel[index] == STAGGERED_FLAG_VALUE_TARGET ) {
            //				cout << "P rls index: " << index << " to be processed\n";
            //				cout << "Starting marching extrapolation on data index " << index << endl;
            staggered_field_marching_extrapolation_recursive(
                index, adj, flagChannel, signedDistanceChannel, staggeredPopulatedChannel, staggeredMask,
                extrapToPositive, dataIndexMapChannel, staggeredExtrapChannel, raeNeighborIndex, dataIndexesToProcess );

            // if( logging::is_logging_debug() ) {
            //	std::cout << "outputing mesh: " << db.counter << std::endl;
            //	db.create_debug_mesh();
            // }
        }
    }
}

} // namespace detail

void staggered_field_marching_extrapolation( const ris_adjacency& adj, const float* signedDistanceChannel,
                                             const boost::uint8_t* staggeredPopulatedChannel, int extrapDirection,
                                             boost::int32_t* dataIndexMapChannel, float* staggeredExtrapChannel )

{
    if( extrapDirection != -1 && extrapDirection != +1 )
        throw runtime_error(
            "staggered_field_marching_extrapolation() - The extrapolation direction must "
            "either be -1 for extrapolating from outside to inside, or +1 for the reverse.  The value provided was " +
            lexical_cast<string>( extrapDirection ) + "." );

    ris_adj_entry dummy;

    // detail::ambiguousCount =0;

    // logging::set_logging_level(5);

    detail::staggered_field_marching_extrapolation_singleface( adj, signedDistanceChannel, staggeredPopulatedChannel,
                                                               /*staggeredMask=*/0x01, extrapDirection == +1,
                                                               dataIndexMapChannel, staggeredExtrapChannel + 0,
                                                               int( &dummy.x_neg - &dummy[0] ) );

    // cout << "U face amibiguous count = " << detail::ambiguousCount << std::endl;
    // detail::ambiguousCount=0;
    detail::staggered_field_marching_extrapolation_singleface( adj, signedDistanceChannel, staggeredPopulatedChannel,
                                                               /*staggeredMask=*/0x02, extrapDirection == +1,
                                                               dataIndexMapChannel, staggeredExtrapChannel + 1,
                                                               int( &dummy.y_neg - &dummy[0] ) );
    // cout << "V face amibiguous count = " << detail::ambiguousCount << std::endl;
    // detail::ambiguousCount=0;
    //	cout << "neighbor index is " << int(&dummy.z_neg - &dummy[0]) << endl;

    detail::staggered_field_marching_extrapolation_singleface( adj, signedDistanceChannel, staggeredPopulatedChannel,
                                                               /*staggeredMask=*/0x04, extrapDirection == +1,
                                                               dataIndexMapChannel, staggeredExtrapChannel + 2,
                                                               int( &dummy.z_neg - &dummy[0] ) );

    // cout << "W face amibiguous count = " << detail::ambiguousCount << std::endl;
}

//
//
// Following are iterative versions of the above functions.
//
//

namespace detail {
// Store signed distance in the priority queue.
// Be careful to swap signs as necessary when you push new items into the queue !
//
// During staggered velocity extrapolation, using the abs value can switch
// the marching order of voxels close to the interface, when the queue
// contains both -ve and +ve signed distances.
// As a side effect this will also speed up the priority queue operations.
/**
 *  The predicate used to order elements in the priority
 * queue used by marching_extrapolation_iterative.
 * The smallest distance is stored at the top of the
 * priority queue.
 */
template <class index_type, class distance_type>
class marching_predicate_no_abs {
  public:
    bool operator()( const std::pair<index_type, distance_type>& left,
                     const std::pair<index_type, distance_type>& right ) {
        return left.second > right.second;
    }
};

typedef std::priority_queue<std::pair<boost::int32_t, float>, std::vector<std::pair<boost::int32_t, float>>,
                            detail::marching_predicate_no_abs<boost::int32_t, float>>
    level_set_marching_extrapolation_priority_queue;

/**
 *  Extrapolate channel values from a single adjacent
 * input voxel, matching the channel gradient if it is present.
 */
void extrapolate_channel_values_from_single_neighbor(
    const boost::int32_t destIndex, const boost::int32_t srcIndex, const int srcFace, const float voxelLength,
    const std::vector<marching_extrapolation_channel>& extrapChannels ) {
    BOOST_FOREACH( const marching_extrapolation_channel& ch, extrapChannels ) {
        if( ch.gradientToMatch ) {
            // we have a gradient to match for this channel
            // follow the gradient by taking a linear combination:
            // this value = adjacent value + gradient * voxelLength
            const float extrapWeights[2] = { 1.f, is_neighbor_index_direction_positive( srcFace ) ? -voxelLength
                                                                                                  : voxelLength };
            const char* extrapData[2];
            extrapData[0] = ch.data + srcIndex * ch.primitiveSize;
            extrapData[1] = ch.gradientToMatch + ( 3 * srcIndex + neighbor_index_axis( srcFace ) ) * ch.primitiveSize;
            ch.weightedSumFn( extrapWeights, extrapData, 2, ch.arity, ch.data + destIndex * ch.primitiveSize );
        } else {
            // no gradient matching for this voxel, so just copy
            // the adjacent value
            memcpy( ch.data + destIndex * ch.primitiveSize, ch.data + srcIndex * ch.primitiveSize, ch.primitiveSize );
        }
    }
}

/**
 *  Extrapolate channel values from multiple adjacent voxels,
 * matching the channel gradient if it is present.
 *
 * @todo consider pre-allocating the tempData buffer as a part of
 *		the marching_extrapolation_channel structure
 */
void extrapolate_channel_values_from_neighbors( const boost::int32_t destIndex, const int inputCount,
                                                const float* weights, const boost::int32_t* dataIndices,
                                                const int* neighborFace, const float voxelLength,
                                                const std::vector<marching_extrapolation_channel>& extrapChannels ) {
    if( inputCount < 1 || inputCount > 6 )
        throw std::runtime_error(
            "extrapolate_channel_values_from_neighbors Error: need 0 < inputCount <= 6, but instead inputCount is " +
            boost::lexical_cast<std::string>( inputCount ) );

    char* dataPointers[6]; // Pointers to the extrapolation data

    BOOST_FOREACH( const marching_extrapolation_channel& ch, extrapChannels ) {
        if( ch.gradientToMatch ) {
            // In this case, we have to extrapolate all the input points, and then interpolate the extrapolated values
            char* tempData = new char[ch.primitiveSize * inputCount];
            // the first extrapolation is fixed
            // the second will change sign depending on the neighbor's direction
            float extrapWeights[2] = { 1.f, voxelLength };
            const char* extrapData[2];
            // Extrapolate and fill in the data pointers
            for( int i = 0; i != inputCount; ++i ) {
                // Pointer to the data
                extrapData[0] = ch.data + dataIndices[i] * ch.primitiveSize;
                // Pointer to the element of the data gradient which corresponds to the neighbor axis direction
                extrapData[1] = ch.gradientToMatch +
                                ( 3 * dataIndices[i] + neighbor_index_axis( neighborFace[i] ) ) * ch.primitiveSize;
                // The multiplication factor for the gradient is +/-voxelLength
                extrapWeights[1] = is_neighbor_index_direction_positive( neighborFace[i] ) ? -voxelLength : voxelLength;
                // The data pointers point inside our temporary array now
                dataPointers[i] = tempData + i * ch.primitiveSize;
                ch.weightedSumFn( extrapWeights, extrapData, 2, ch.arity, dataPointers[i] );
            }

            // Call the linear combination function to do the extrapolation to this voxel.
            ch.weightedSumFn( weights, dataPointers, inputCount, ch.arity, ch.data + destIndex * ch.primitiveSize );

            delete[] tempData;
        } else { // Equivalent to matching a zero gradient
            // get the data pointer for each input
            for( int i = 0; i < inputCount; ++i )
                dataPointers[i] = ch.data + dataIndices[i] * ch.primitiveSize;
            // linear combination of inputs to extrapolate into this channel
            ch.weightedSumFn( weights, dataPointers, inputCount, ch.arity, ch.data + destIndex * ch.primitiveSize );
        }
    }
}

/**
 *  Extrapolate channel values into the specified voxel from its populated
 * neighbours.  This function assumes that information flows away from the
 * phi = 0 isosurface.
 *
 * This code is from marching_extrapolation_recursive.
 */
bool extrapolate_voxel_channels_from_populated_neighbours(
    const size_t dataIndex, const ris_adjacency& adj, const float voxelLength, unsigned char* populatedChannel,
    const float* signedDistanceChannel, const std::vector<marching_extrapolation_channel>& extrapChannels ) {
    const float signedDistance = signedDistanceChannel[dataIndex];
    const ris_adj_entry& rae = adj[dataIndex];

    // do we look for smaller or larger signed distances as inputs ?
    // const signed char inputFlag = signedDistance >= 0 ? -1 : +1;

    // inputFlag and adjSignedDistance specify information about the 6 neighbor voxels
    // inputFlag is a boolean that indicates whether the corresponding face is an input
    char inputFlag[6];
    float adjSignedDistance[6];

    int inputCount = 0;
    for( int face = 0; face < 6; ++face ) {
        const boost::int32_t adjDataIndex = rae[face];

        if( adjDataIndex >= 0 && populatedChannel[adjDataIndex] == 1 ) {
            const float faceSignedDistance = signedDistanceChannel[adjDataIndex];
            if( signedDistance >= 0 && faceSignedDistance < signedDistance ) {
                adjSignedDistance[face] = faceSignedDistance;
                inputFlag[face] = 1;
                ++inputCount;
            } else if( signedDistance < 0 && faceSignedDistance > signedDistance ) {
                adjSignedDistance[face] = faceSignedDistance;
                inputFlag[face] = 1;
                ++inputCount;
            } else {
                inputFlag[face] = 0;
            }
        } else {
            inputFlag[face] = 0;
        }
    }

    if( inputCount == 1 ) {
        // this voxel has a single input
        for( int face = 0; face < 6; ++face ) {
            if( inputFlag[face] ) {
                const boost::int32_t adjDataIndex = rae[face];
                extrapolate_channel_values_from_single_neighbor( static_cast<boost::int32_t>( dataIndex ), adjDataIndex,
                                                                 face, voxelLength, extrapChannels );
                // we're done processing our single input, so
                // we can break early
                break;
            }
        }
    } else if( inputCount > 1 ) {
        // Collect all the input data indexes and weights into arrays.
        float totalWeight = 0;         // The sum of all the extrapolation weights
        float weights[6];              // The extrapolation weights
        boost::int32_t dataIndices[6]; // The data indices for the extrapolation data
        int neighborFace[6];           // Neighbor index of the extrapolation data

        // gather the input weights and indices for all input neighbours
        int inputNumber = 0;
        for( int axis = 0; axis < 3; ++axis ) {
            const int face0 = get_rae_neighbor_index( axis, -1 );
            const int face1 = get_rae_neighbor_index( axis, +1 );

            if( inputFlag[face0] && inputFlag[face1] ) {
                const boost::int32_t dataIndex0 = rae[face0];
                const boost::int32_t dataIndex1 = rae[face1];

                float weight0 = fabs( signedDistanceChannel[dataIndex0] - signedDistance );
                float weight1 = fabs( signedDistanceChannel[dataIndex1] - signedDistance );

                // scale so that this axis is not over-represented
                // compared to axes with only on input
                const float scaleRatio = max( weight0, weight1 ) / ( weight0 + weight1 );
                weight0 *= scaleRatio;
                weight1 *= scaleRatio;

                neighborFace[inputNumber] = face0;
                dataIndices[inputNumber] = dataIndex0;
                weights[inputNumber] = weight0;
                totalWeight += weight0;
                ++inputNumber;

                neighborFace[inputNumber] = face1;
                dataIndices[inputNumber] = dataIndex1;
                weights[inputNumber] = weight1;
                totalWeight += weight1;
                ++inputNumber;
            } else if( inputFlag[face0] ) {
                const boost::int32_t dataIndex = rae[face0];

                const float weight = fabs( signedDistanceChannel[dataIndex] - signedDistance );

                neighborFace[inputNumber] = face0;
                dataIndices[inputNumber] = dataIndex;
                weights[inputNumber] = weight;
                totalWeight += weight;
                ++inputNumber;
            } else if( inputFlag[face1] ) {
                const boost::int32_t dataIndex = rae[face1];

                const float weight = fabs( signedDistanceChannel[dataIndex] - signedDistance );

                neighborFace[inputNumber] = face1;
                dataIndices[inputNumber] = dataIndex;
                weights[inputNumber] = weight;
                totalWeight += weight;
                ++inputNumber;
            }
        }

        assert( inputNumber == inputCount );

        if( totalWeight > 0 ) {
            // normalize the weights
            const float invTotalWeight = 1.f / totalWeight;
            for( int face = 0; face < 6; ++face ) {
                weights[face] *= invTotalWeight;
            }
        } else {
            // previously used 1.f/inputCount for all weights,
            // but I don't think this should come up ?
            throw std::runtime_error(
                "extrapolate_voxel_channels_from_populated_neighbours Error: voxel had inputCount > 1 "
                "but did not get a positive weight.  The total weight was " +
                boost::lexical_cast<std::string>( totalWeight ) );
        }

        extrapolate_channel_values_from_neighbors( static_cast<boost::int32_t>( dataIndex ), inputCount, weights,
                                                   dataIndices, neighborFace, voxelLength, extrapChannels );
    } else {
        return false;
    }

    return true;
}

/**
 *  Return true if information can flow from the srcIndex voxel to the
 * destIndex voxel, assuming that information flows away from the phi = 0
 * isosurface.
 */
bool has_information_flow( const size_t destIndex, const size_t srcIndex, const float* signedDistanceChannel ) {
    const float destPhi = signedDistanceChannel[destIndex];
    const float srcPhi = signedDistanceChannel[srcIndex];

    if( destPhi >= 0 ) {
        if( srcPhi < destPhi )
            return true;
    } else {
        if( srcPhi > destPhi )
            return true;
    }

    return false;
}

/**
 *  Look up the signed distance for the specified face
 * of the specified voxel.
 */
float get_face_signed_distance( const ris_adjacency& adj, const float* signedDistanceChannel,
                                const std::size_t dataIndex, const int raeNeighborIndex ) {
    const float centerPhi = signedDistanceChannel[dataIndex];
    const boost::int32_t adjDataIndex = adj[dataIndex][raeNeighborIndex];

    if( adjDataIndex >= 0 ) {
        const float adjPhi = signedDistanceChannel[adjDataIndex];
        return 0.5f * ( centerPhi + adjPhi );
    } else {
        return centerPhi;
    }
}

/**
 *  Return true if information can flow to the level set distance
 * destPhi from srcPhi.
 */
bool has_information_flow( const float destPhi, const float srcPhi, const bool extrapToPositive ) {
    if( extrapToPositive ) {
        if( destPhi > srcPhi ) {
            return true;
        }
    } else { // extrap to negative
        if( destPhi < srcPhi ) {
            return true;
        }
    }

    return false;
}

/**
 *  Assume that the specified voxel is already populated.  Look at the
 * neighbours of the specified voxel, and add them to the near queue if
 * they are unpopulated and information can propagate from the
 * specified voxel to the neighbour.  This assumes that the specified voxel
 * is already populated.
 */
void add_unpopulated_neighbours_to_near( const size_t dataIndex, const ris_adjacency& adj,
                                         unsigned char* populatedChannel, const float* signedDistanceChannel,
                                         detail::level_set_marching_extrapolation_priority_queue& nearDataIndices,
                                         const bool addPositiveNeighboursOnly = false ) {
    const ris_adj_entry& rae = adj[dataIndex];

    for( int face = 0; face < 6; ++face ) {
        if( addPositiveNeighboursOnly && !is_neighbor_index_direction_positive( face ) )
            continue;

        const boost::int32_t adjIndex = rae[face];
        if( adjIndex >= 0 && populatedChannel[adjIndex] == 0 &&
            has_information_flow( adjIndex, dataIndex, signedDistanceChannel ) ) {
            // FF_LOG( debug ) << "adding " << signedDistanceChannel[dataIndex] << " -> " <<
            // signedDistanceChannel[adjIndex]
            // << std::endl;
            nearDataIndices.push(
                std::pair<boost::int32_t, float>( adjIndex, fabsf( signedDistanceChannel[adjIndex] ) ) );
            populatedChannel[adjIndex] = 4;
        }
    }
}

} // namespace detail

void marching_extrapolation_iterative( const ris_adjacency& adj, const float voxelLength,
                                       unsigned char* populatedChannel, const float* signedDistanceChannel,
                                       const std::vector<marching_extrapolation_channel>& extrapChannels,
                                       extrapolation_debug_info& /*db*/
                                       /*std::size_t* populatedOrderDebugChannel*/ ) {
    if( voxelLength <= 0 )
        throw std::runtime_error(
            "marching_extrapolation_iterative Error: voxelLength must be positive, but instead it is " +
            boost::lexical_cast<std::string>( voxelLength ) + "." );

    if( populatedChannel == 0 )
        throw std::runtime_error( "marching_extrapolation_iterative Error: the populatedChannel is NULL." );

    if( signedDistanceChannel == 0 )
        throw std::runtime_error( "marching_extrapolation_iterative Error: the signedDistanceChannel is NULL." );

    detail::level_set_marching_extrapolation_priority_queue nearDataIndices;

    // if( logging::is_logging_debug() )
    // db.create_debug_mesh();

    // std::size_t populatedOrder = 0;

    // get the "near" voxels and place them in the priority queue
    for( size_t i = 0; i != adj.data_size(); ++i ) {
        switch( populatedChannel[i] ) {
        case 0: { // unpopulated extrapolation target
            const ris_adj_entry& rae = adj[i];
            // if there's an adjacent populated voxel, then this voxel is
            // "near", and the neighbouring voxel is "accepted"
            for( int axis = 0; axis < 3; ++axis ) {
                const int adjIndex = rae[get_positive_rae_neighbor_index( axis )];
                if( adjIndex >= 0 && populatedChannel[adjIndex] == 1 &&
                    detail::has_information_flow( i, adjIndex, signedDistanceChannel ) ) {
                    const float phi = signedDistanceChannel[i];
                    nearDataIndices.push(
                        std::pair<boost::int32_t, float>( static_cast<boost::int32_t>( i ), fabsf( phi ) ) );
                    populatedChannel[i] = 4;
                    break;
                }
            }
            break;
        }
        case 1: { // populated extrapolation source
            detail::add_unpopulated_neighbours_to_near( i, adj, populatedChannel, signedDistanceChannel,
                                                        nearDataIndices, true );
            break;
        }
        }
    }

    // go through the near voxels, in order of ascending distance
    while( !nearDataIndices.empty() ) {
        frantic::logging::check_global_abort();

        std::pair<boost::int32_t, float> trial = nearDataIndices.top();
        nearDataIndices.pop();

        const boost::int32_t trialDataIndex = trial.first;

        if( populatedChannel[trialDataIndex] == 4 ) {
            // try to set this voxel using its populated neighbours
            if( detail::extrapolate_voxel_channels_from_populated_neighbours(
                    trialDataIndex, adj, voxelLength, populatedChannel, signedDistanceChannel, extrapChannels ) ) {
                // we successfully set this voxel's channels.

                // if( populatedOrderDebugChannel )
                // populatedOrderDebugChannel[trialDataIndex] = populatedOrder++;

                // FF_LOG( debug ) << "populated voxel at distance " << signedDistanceChannel[ trialDataIndex ] <<
                // std::endl;

                // add this voxel's unpopulated neighbours to near
                detail::add_unpopulated_neighbours_to_near( trialDataIndex, adj, populatedChannel,
                                                            signedDistanceChannel, nearDataIndices );

                // flag this voxel as populated/accepted
                populatedChannel[trialDataIndex] = 1;
            } else {
                // No data was available to set this voxel's channels.
                // We indicate this situation by settings its
                // populated flag to 3.
                populatedChannel[trialDataIndex] = 3;
            }
        } else {
            throw std::runtime_error( "marching_extrapolation_iterative Error: a voxel in the near list has a "
                                      "populated flag different from 4: " +
                                      boost::lexical_cast<std::string>( populatedChannel[trialDataIndex] ) );
        }
    }

    // if( logging::is_logging_debug() )
    // db.create_debug_mesh();
}

namespace detail {
/**
 * @todo should we allow extrapolation on an isosurface ( srcPhi == destPhi ? )
 * @return true if the cell was extrapolated into, or false otherwise.
 */
bool extrapolate_staggered_channel_from_populated_neighbours_singleface(
    const boost::int32_t dataIndex, const ris_adjacency& adj, const unsigned char* flagChannel,
    const float* signedDistanceChannel, const bool extrapToPositive, const boost::int32_t* dataIndexMapChannel,
    float* staggeredExtrapChannel, const int raeNeighborIndex ) {
    const float targetSignedDistance =
        detail::get_face_signed_distance( adj, signedDistanceChannel, dataIndex, raeNeighborIndex );
    const ris_adj_entry& rae = adj[dataIndex];

    // different from centred case : specified vs. a function of signed distance
    // const signed char inputFlag = extrapToPositive ? -1 : +1;

    char inputFlag[6];
    float adjSignedDistance[6];

    int inputCount = 0;
    for( int face = 0; face < 6; ++face ) {
        const boost::int32_t adjDataIndex = rae[face];
        if( adjDataIndex >= 0 && flagChannel[adjDataIndex] == STAGGERED_FLAG_VALUE_SOURCE ) {
            const float srcSignedDistance =
                detail::get_face_signed_distance( adj, signedDistanceChannel, adjDataIndex, raeNeighborIndex );
            adjSignedDistance[face] = srcSignedDistance;
            if( extrapToPositive && srcSignedDistance < targetSignedDistance ) {
                inputFlag[face] = 1;
                ++inputCount;
            } else if( !extrapToPositive && srcSignedDistance > targetSignedDistance ) {
                inputFlag[face] = 1;
                ++inputCount;
            } else {
                inputFlag[face] = 0;
            }
        } else {
            inputFlag[face] = 0;
        }
    }

    if( inputCount == 1 ) {
        for( int face = 0; face < 6; ++face ) {
            if( inputFlag[face] ) {
                staggeredExtrapChannel[3 * dataIndexMapChannel[dataIndex]] =
                    staggeredExtrapChannel[3 * dataIndexMapChannel[rae[face]]];
            }
        }
    } else if( inputCount > 1 ) {
        // Collect all the input data indexes and weights into arrays.
        float totalWeight = 0;
        float weights[6];
        int dataIndices[6];

        int inputNumber = 0;
        // Go through i == 0, 2, 4, and process +/- directions along each axis together
        for( int axis = 0; axis < 3; ++axis ) {
            const int negFace = get_negative_rae_neighbor_index( axis );
            const int posFace = get_positive_rae_neighbor_index( axis );

            // If there are inputs in both the + and - directions along this axis
            if( inputFlag[negFace] && inputFlag[posFace] ) {
                // Get the indices of the two neighbors
                const boost::int32_t dataIndex0 = rae[negFace];
                const boost::int32_t dataIndex1 = rae[posFace];

                // Compute the weights, and adjust them so that this axis doesn't get over-represented because
                // of having two inputs instead of one
                float weight0 = fabsf( adjSignedDistance[negFace] - targetSignedDistance );
                float weight1 = fabsf( adjSignedDistance[posFace] - targetSignedDistance );

                const float scaleRatio = max( weight0, weight1 ) / ( weight0 + weight1 );
                weight0 *= scaleRatio;
                weight1 *= scaleRatio;

                // Add the weights and data indices to the inputs
                dataIndices[inputNumber] = dataIndex0;
                weights[inputNumber] = weight0;
                totalWeight += weight0;
                ++inputNumber;

                dataIndices[inputNumber] = dataIndex1;
                weights[inputNumber] = weight1;
                totalWeight += weight1;
                ++inputNumber;
            }
            // Otherwise only one or none of the directions will be an input
            else if( inputFlag[negFace] ) {
                // Get the index of the neighbor
                const boost::int32_t srcDataIndex = rae[negFace];

                // Compute the weight to use
                float weight = fabsf( adjSignedDistance[negFace] - targetSignedDistance );

                // Add the weight and data index to the inputs
                dataIndices[inputNumber] = srcDataIndex;
                weights[inputNumber] = weight;
                totalWeight += weight;
                ++inputNumber;
            } else if( inputFlag[posFace] ) {
                // Get the index of the neighbor
                const boost::int32_t srcDataIndex = rae[posFace];

                // Compute the weight to use
                float weight = fabsf( adjSignedDistance[posFace] - targetSignedDistance );

                // Add the weight and data index to the inputs
                dataIndices[inputNumber] = srcDataIndex;
                weights[inputNumber] = weight;
                totalWeight += weight;
                ++inputNumber;
            }
        }

        if( inputNumber != inputCount ) {
            throw std::runtime_error(
                "extrapolate_staggered_channel_from_populated_neighbours_singleface Error: mismatch "
                "between counted and culled extrapolation sources: " +
                boost::lexical_cast<std::string>( inputCount ) + " vs " +
                boost::lexical_cast<std::string>( inputNumber ) );
        }

        // Normalize the weights if there was any actual weight collected
        if( totalWeight > 0 ) {
            const float invTotalWeight = 1.f / totalWeight;
            for( int i = 0; i != inputCount; ++i )
                weights[i] *= invTotalWeight;
        } else {
            // If no weight at all was collected, take the average of all the inputs.
            float weight = 1.f / inputCount;
            for( int i = 0; i != inputCount; ++i )
                weights[i] = weight;
        }

        // Compute the extrapolated value, and save it
        float extrapolatedValue = 0;
        for( int i = 0; i != inputCount; ++i ) {
            extrapolatedValue += weights[i] * staggeredExtrapChannel[3 * dataIndexMapChannel[dataIndices[i]]];
        }
        staggeredExtrapChannel[3 * dataIndexMapChannel[dataIndex]] = extrapolatedValue;

    } else { // no inputs to extrapolate from
        return false;
    }

    return true;
}

/**
 *  This is called after dataIndex becomes a new source ( you have
 * already marked it as populated, or you are going to before
 * continuing extrapolation. )
 *
 *  Look for neighbours which can now be extrapolated into from
 * the new source at dataIndex.
 */
void add_unpopulated_neighbours_to_near_staggered_singleface(
    const boost::int32_t dataIndex, const ris_adjacency& adj, const float* signedDistanceChannel,
    unsigned char* flagChannel, const boost::uint8_t* staggeredPopulatedChannel, const boost::uint8_t staggeredMask,
    const bool extrapToPositive, const boost::int32_t* dataIndexMapChannel, const int raeNeighborIndex,
    detail::level_set_marching_extrapolation_priority_queue& nearDataIndices ) {
    const ris_adj_entry& rae = adj[dataIndex];
    const float srcSignedDistance =
        detail::get_face_signed_distance( adj, signedDistanceChannel, dataIndex, raeNeighborIndex );

    for( int face = 0; face < 6; ++face ) {
        const boost::int32_t adjDataIndex = rae[face];
        if( adjDataIndex >= 0 ) {
            const unsigned char flag = detail::set_staggered_voxel_flag( flagChannel, (int)adjDataIndex,
                                                                         staggeredPopulatedChannel, staggeredMask );
            if( flag == STAGGERED_FLAG_VALUE_TARGET && dataIndexMapChannel[adjDataIndex] >= 0 ) {
                const float targetSignedDistance =
                    detail::get_face_signed_distance( adj, signedDistanceChannel, adjDataIndex, raeNeighborIndex );

                if( detail::has_information_flow( targetSignedDistance, srcSignedDistance, extrapToPositive ) ) {
                    nearDataIndices.push( std::pair<boost::int32_t, float>(
                        adjDataIndex, extrapToPositive ? targetSignedDistance : -targetSignedDistance ) );
                    flagChannel[adjDataIndex] = STAGGERED_FLAG_PROCESSING;
                }
            }
        }
    }
}

/**
 * This function does marching extrapolation on the U, V or W component of the staggered channel.
 *
 * @note Look at the main function definition to see what most of the other parameters mean.
 *
 * @param  staggeredExtrapChannel  This is the staggeredExtrapChannel from the main function, but
 *                                 with +0, +1 or +2 added depending on which of U, V, or W faces
 *                                 are being processed.
 * @param  raeNeighborIndex  This is the index within an rae, which will point to one of rae.x_neg,
 *                           rae.y_neg, or rae.z_neg.
 */
void staggered_field_marching_extrapolation_singleface_iterative(
    const ris_adjacency& adj, const float* signedDistanceChannel, const boost::uint8_t* staggeredPopulatedChannel,
    boost::uint8_t staggeredMask, bool extrapToPositive, boost::int32_t* dataIndexMapChannel,
    float* staggeredExtrapChannel, int raeNeighborIndex ) {
    detail::level_set_marching_extrapolation_priority_queue nearDataIndices;

    if( adj.data_size() == 0 )
        return;

    // These data indexes refer to signedDistanceChannel
    // vector<boost::int32_t> dataIndexesToProcess;

    // Create a flag channel (starting with all 0 values)
    vector<unsigned char> flagChannelVector( adj.data_size(), STAGGERED_FLAG_UNTOUCHED );
    unsigned char* flagChannel = &flagChannelVector[0];

    // This is a little different from the original, as I
    // change a voxel's flag to STAGGERED_FLAG_PROCESSING
    // when it is enqueued.

    for( size_t i = 0, ie = adj.data_size(); i != ie; ++i ) {
        // Only process voxels which have a corresponding value in the staggered field
        if( dataIndexMapChannel[i] >= 0 ) {
            const ris_adj_entry& rae = adj[i];
            // Determine whether this value is a source, a target, or neither
            unsigned char flag =
                detail::set_staggered_voxel_flag( flagChannel, (int)i, staggeredPopulatedChannel, staggeredMask );
            //			cout << "M rls index: " << i << ", rvf index: " << dataIndexMapChannel[i] << ", flag: " << (int)flag
            //<<
            //"\n";
            switch( flag ) {
            case STAGGERED_FLAG_VALUE_TARGET: {
                bool crossesBoundary = false;
                const float targetSignedDistance =
                    detail::get_face_signed_distance( adj, signedDistanceChannel, i, raeNeighborIndex );
                for( int axis = 0; axis < 3; ++axis ) {
                    const boost::int32_t adjDataIndex = rae[get_positive_rae_neighbor_index( axis )];
                    if( adjDataIndex >= 0 ) {
                        flag = detail::set_staggered_voxel_flag( flagChannel, adjDataIndex, staggeredPopulatedChannel,
                                                                 staggeredMask );
                        if( flag == STAGGERED_FLAG_VALUE_SOURCE ) {
                            const float srcSignedDistance = detail::get_face_signed_distance(
                                adj, signedDistanceChannel, adjDataIndex, raeNeighborIndex );
                            if( detail::has_information_flow( targetSignedDistance, srcSignedDistance,
                                                              extrapToPositive ) ) {
                                crossesBoundary = true;
                                break;
                            }
                        }
                    }
                }
                // If there's an adjacent source voxel, then start marching from this voxel
                if( crossesBoundary ) {
                    nearDataIndices.push( std::pair<boost::int32_t, float>(
                        static_cast<boost::int32_t>( i ),
                        extrapToPositive ? targetSignedDistance : -targetSignedDistance ) );
                    flagChannel[i] = STAGGERED_FLAG_PROCESSING;
                }
                break;
            }
            case STAGGERED_FLAG_VALUE_SOURCE: {
                const float srcSignedDistance =
                    detail::get_face_signed_distance( adj, signedDistanceChannel, i, raeNeighborIndex );
                // If there's an adjacent target voxel, then start marching from that voxel
                for( int axis = 0; axis < 3; ++axis ) {
                    const boost::int32_t adjDataIndex = rae[get_positive_rae_neighbor_index( axis )];
                    if( adjDataIndex >= 0 ) {
                        flag = detail::set_staggered_voxel_flag( flagChannel, adjDataIndex, staggeredPopulatedChannel,
                                                                 staggeredMask );
                        if( flag == STAGGERED_FLAG_VALUE_TARGET && dataIndexMapChannel[adjDataIndex] >= 0 ) {
                            const float targetSignedDistance = detail::get_face_signed_distance(
                                adj, signedDistanceChannel, adjDataIndex, raeNeighborIndex );
                            if( detail::has_information_flow( targetSignedDistance, srcSignedDistance,
                                                              extrapToPositive ) ) {
                                nearDataIndices.push( std::pair<boost::int32_t, float>(
                                    adjDataIndex, extrapToPositive ? targetSignedDistance : -targetSignedDistance ) );
                                flagChannel[adjDataIndex] = STAGGERED_FLAG_PROCESSING;
                            }
                        }
                    }
                }
                break;
            }
            }
        }
    }

    while( !nearDataIndices.empty() ) {
        std::pair<boost::int32_t, float> trial = nearDataIndices.top();
        nearDataIndices.pop();

        const boost::int32_t trialDataIndex = trial.first;

        if( flagChannel[trialDataIndex] == STAGGERED_FLAG_PROCESSING ) {
            // try to set this voxel using its populated neighbours
            if( detail::extrapolate_staggered_channel_from_populated_neighbours_singleface(
                    trialDataIndex, adj, flagChannel, signedDistanceChannel, extrapToPositive, dataIndexMapChannel,
                    staggeredExtrapChannel, raeNeighborIndex ) ) {
                // we successfully set this voxel's channels.

                // if( populatedOrderDebugChannel )
                // populatedOrderDebugChannel[trialDataIndex] = populatedOrder++;

                // FF_LOG( debug ) << "populated voxel at distance " << signedDistanceChannel[ trialDataIndex ] <<
                // std::endl;

                // add this voxel's unpopulated neighbours to near
                detail::add_unpopulated_neighbours_to_near_staggered_singleface(
                    trialDataIndex, adj, signedDistanceChannel, flagChannel, staggeredPopulatedChannel, staggeredMask,
                    extrapToPositive, dataIndexMapChannel, raeNeighborIndex, nearDataIndices );

                // flag this voxel as populated/accepted
                flagChannel[trialDataIndex] = STAGGERED_FLAG_VALUE_SOURCE;
            } else {
                // No data was available to set this voxel's channels.
                // We indicate this situation by settings its
                // populated flag to 3.
                flagChannel[trialDataIndex] = STAGGERED_FLAG_IGNORE;
            }
        } else {
            throw std::runtime_error( "staggered_field_marching_extrapolation_singleface_iterative Error: unexpected "
                                      "flag for queued voxel: " +
                                      boost::lexical_cast<std::string>( (unsigned int)flagChannel[trialDataIndex] ) );
        }
    }
}

} // namespace detail

void staggered_field_marching_extrapolation_iterative( const ris_adjacency& adj, const float* signedDistanceChannel,
                                                       const boost::uint8_t* staggeredPopulatedChannel,
                                                       int extrapDirection, boost::int32_t* dataIndexMapChannel,
                                                       float* staggeredExtrapChannel ) {
    if( extrapDirection != -1 && extrapDirection != +1 )
        throw runtime_error(
            "staggered_field_marching_extrapolation_iterative() - The extrapolation direction must "
            "either be -1 for extrapolating from outside to inside, or +1 for the reverse.  The value provided was " +
            lexical_cast<string>( extrapDirection ) + "." );

    // detail::ambiguousCount =0;

    // logging::set_logging_level(5);

    detail::staggered_field_marching_extrapolation_singleface_iterative(
        adj, signedDistanceChannel, staggeredPopulatedChannel,
        /*staggeredMask=*/0x01, extrapDirection == +1, dataIndexMapChannel, staggeredExtrapChannel + 0,
        rae_index_x_neg );

    // cout << "U face amibiguous count = " << detail::ambiguousCount << std::endl;
    // detail::ambiguousCount=0;
    detail::staggered_field_marching_extrapolation_singleface_iterative(
        adj, signedDistanceChannel, staggeredPopulatedChannel,
        /*staggeredMask=*/0x02, extrapDirection == +1, dataIndexMapChannel, staggeredExtrapChannel + 1,
        rae_index_y_neg );
    // cout << "V face amibiguous count = " << detail::ambiguousCount << std::endl;
    // detail::ambiguousCount=0;
    //	cout << "neighbor index is " << int(&dummy.z_neg - &dummy[0]) << endl;

    detail::staggered_field_marching_extrapolation_singleface_iterative(
        adj, signedDistanceChannel, staggeredPopulatedChannel,
        /*staggeredMask=*/0x04, extrapDirection == +1, dataIndexMapChannel, staggeredExtrapChannel + 2,
        rae_index_z_neg );

    // cout << "W face amibiguous count = " << detail::ambiguousCount << std::endl;
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
