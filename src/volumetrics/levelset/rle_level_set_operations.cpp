// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/levelset/rle_pairwise_run_iterator.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::geometry;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;

namespace frantic {
namespace volumetrics {
namespace levelset {

void rle_level_set::csg_complement() {
    // Iterate over all the runs
    for( int c = 0, csize = m_rleIndex.m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_rleIndex.m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;

            int runRangeStart = m_rleIndex.m_bcToRunIndex[bcIndex],
                runRangeEnd = m_rleIndex.m_bcToRunIndex[bcIndex + 1] - 1;
            run_data *rd = &m_rleIndex.m_runData[runRangeStart], *rdEnd = &m_rleIndex.m_runData[runRangeEnd];

            // For each run go in and switch the sign of undefined regions
            for( ; rd != rdEnd; ++rd ) {
                boost::int32_t index = rd->dataIndex;
                if( index < 0 ) {
                    // if the region is undefined swap the region from inside to outside and vice versa
                    rd->dataIndex = ( index == -1 ) ? -2 : -1;
                }
            }
        }
    }

    // invert the distance function for defined voxels
    for( size_t i = 0, ie = m_distanceData.size(); i != ie; ++i )
        m_distanceData[i] = -m_distanceData[i];

    // Invert the outside area too
    m_rleIndex.m_exteriorRegionCode = ( m_rleIndex.m_exteriorRegionCode == -1 ) ? -2 : -1;

    // Switch the inside/outside values
    std::swap( m_interfaceVoxelWidthInside, m_interfaceVoxelWidthOutside );
    std::swap( m_insideDistance, m_outsideDistance );
    // invert signs
    m_insideDistance = -m_insideDistance;
    m_outsideDistance = -m_outsideDistance;
}

template <class CSGPolicy>
void rle_level_set::csg_operation( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                                   float channelBlendVoxelDistance ) {
    // Check that rleFirst and rleSecond have the same voxel coordinate systems.
    // Use a relative comparison, so that if the vcs was computed in two different ways,
    // but is trying to be the same, it will match.
    if( !rleFirst.m_voxelCoordSystem.equals( rleSecond.m_voxelCoordSystem ) )
        throw std::runtime_error(
            "rle_level_set.csg_operation: The voxel coordinate systems of the two level set operands "
            "don't match.  The first is " +
            rleFirst.m_voxelCoordSystem.str() + ", while the second is " + rleSecond.m_voxelCoordSystem.str() + "." );
    if( !m_voxelCoordSystem.equals( rleFirst.m_voxelCoordSystem ) )
        throw std::runtime_error( "rle_level_set.csg_operation: The voxel coordinate systems of the input level set "
                                  "operands doesn't match that of the output.  The input is " +
                                  rleFirst.m_voxelCoordSystem.str() + ", while the output is " +
                                  m_voxelCoordSystem.str() + "." );

    rle_index_spec& ris = m_rleIndex;

    // Get the exterior region codes
    boost::int32_t firstExteriorRegionCode = rleFirst.get_rle_index_spec().get_exterior_region_code();
    boost::int32_t secondExteriorRegionCode = rleSecond.get_rle_index_spec().get_exterior_region_code();
    // Get bools for whether the exterior codes dominate
    bool firstExteriorDominates = CSGPolicy::first_region_code_dominates( firstExteriorRegionCode );
    bool secondExteriorDominates = CSGPolicy::second_region_code_dominates( secondExteriorRegionCode );

    // get the channel handling parameters
    bool useChannelsFromFirst = CSGPolicy::use_first_op_channels();
    bool useChannelsFromSecond = CSGPolicy::use_second_op_channels();
    bool onlyChannelsFromFirst = useChannelsFromFirst && !useChannelsFromSecond;
    bool onlyChannelsFromSecond = !useChannelsFromFirst && useChannelsFromSecond;

    // Do a special case if either level set is completely empty
    if( rleFirst.m_rleIndex.is_empty() ) {
        // If rleFirst is empty, copy rleFirst or rleSecond depending on whether it's entirely inside or entirely
        // outside, respectively.
        if( firstExteriorDominates ) {
            *this = rleFirst;
            // Complement if the first operand needs it
            if( CSGPolicy::first_op_complemented() )
                csg_complement();
        } else {
            *this = rleSecond;
            // Complement if the second operand needs it
            if( CSGPolicy::second_op_complemented() )
                csg_complement();
        }
        return;
    } else if( rleSecond.m_rleIndex.is_empty() ) {
        // If rleSecond is empty, copy rleFirst or rleSecond depending on whether it's entirely outside or entirely
        // inside, respectively.
        if( secondExteriorDominates ) {
            *this = rleSecond;
            // Complement if the second operand needs it
            if( CSGPolicy::second_op_complemented() )
                csg_complement();
        } else {
            *this = rleFirst;
            // Complement if the first operand needs it
            if( CSGPolicy::first_op_complemented() )
                csg_complement();
        }
        return;
    }

    // Since the signed distances are in world space, need to convert the channel blend distance into world space
    float channelBlendWorldDistance = channelBlendVoxelDistance * m_voxelCoordSystem.voxel_length();

    ris.m_exteriorRegionCode = CSGPolicy::combine_region_codes( firstExteriorRegionCode, secondExteriorRegionCode );

    // The voxel interface width of the result will be the smaller of the voxel interface widths of the two inputs.
    m_interfaceVoxelWidthInside =
        ( std::min )( rleFirst.m_interfaceVoxelWidthInside, rleSecond.m_interfaceVoxelWidthInside );
    m_interfaceVoxelWidthOutside =
        ( std::min )( rleFirst.m_interfaceVoxelWidthOutside, rleSecond.m_interfaceVoxelWidthOutside );
    m_outsideDistance = ( std::min )( rleFirst.m_outsideDistance, rleSecond.m_outsideDistance );
    m_insideDistance = ( std::max )( rleFirst.m_insideDistance, rleSecond.m_insideDistance );

    ////////////
    // Figure out the bounding box required by the new level set.  There are 4 cases, depending on the exterior region
    // codes of the operands.
    ////////////
    if( firstExteriorDominates && secondExteriorDominates ) {
        // If both exterior region codes are dominant, then we take the intersection of the bounding boxes
        ris.m_abcCoordOrigin =
            vector3::from_max( rleFirst.m_rleIndex.m_abcCoordOrigin, rleSecond.m_rleIndex.m_abcCoordOrigin );
        ris.m_abcCoordSize = size3::from_bounds(
            ris.m_abcCoordOrigin,
            vector3::from_min( rleFirst.m_rleIndex.m_abcCoordOrigin + rleFirst.m_rleIndex.m_abcCoordSize,
                               rleSecond.m_rleIndex.m_abcCoordOrigin + rleSecond.m_rleIndex.m_abcCoordSize ) -
                vector3( 1 ) );
    } else if( firstExteriorDominates ) {
        // If just the first exterior dominates, use the first bounding box
        ris.m_abcCoordOrigin = rleFirst.m_rleIndex.m_abcCoordOrigin;
        ris.m_abcCoordSize = rleFirst.m_rleIndex.m_abcCoordSize;
    } else if( secondExteriorDominates ) {
        // If just the second exterior dominates, use the second bounding box
        ris.m_abcCoordOrigin = rleSecond.m_rleIndex.m_abcCoordOrigin;
        ris.m_abcCoordSize = rleSecond.m_rleIndex.m_abcCoordSize;
    } else {
        // If neither exterior dominates, then we take the union of the bounding boxes
        ris.m_abcCoordOrigin =
            vector3::from_min( rleFirst.m_rleIndex.m_abcCoordOrigin, rleSecond.m_rleIndex.m_abcCoordOrigin );
        ris.m_abcCoordSize = size3::from_bounds(
            ris.m_abcCoordOrigin,
            vector3::from_max( rleFirst.m_rleIndex.m_abcCoordOrigin + rleFirst.m_rleIndex.m_abcCoordSize,
                               rleSecond.m_rleIndex.m_abcCoordOrigin + rleSecond.m_rleIndex.m_abcCoordSize ) -
                vector3( 1 ) );
    }

    // clear any data already stored in ris
    ris.m_runData.clear();
    m_distanceData.clear();

    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );

    ////////////
    // Determine what named channels the CSG result will contain
    ////////////
    m_namedChannels.clear();

    // Add all the named channels from the first operand.
    if( useChannelsFromFirst ) {
        for( map<frantic::tstring, rle_channel>::const_iterator i = rleFirst.m_namedChannels.begin(),
                                                                ie = rleFirst.m_namedChannels.end();
             i != ie; ++i ) {
            add_channel( i->second.name(), i->second.arity(), i->second.data_type() );
        }
    }
    // Add all the named channels from the second operand, which weren't already added from the first
    if( useChannelsFromSecond ) {
        for( map<frantic::tstring, rle_channel>::const_iterator i = rleSecond.m_namedChannels.begin(),
                                                                ie = rleSecond.m_namedChannels.end();
             i != ie; ++i ) {
            map<frantic::tstring, rle_channel>::const_iterator j = m_namedChannels.find( i->second.name() );
            if( j != m_namedChannels.end() ) {
                // If the two input channels aren't the same, throw an error
                if( i->second.arity() != j->second.arity() || i->second.data_type() != j->second.data_type() )
                    throw runtime_error(
                        "rle_level_set.csg_operation: The two input level sets have an incompatible channel named \"" +
                        frantic::strings::to_string( i->second.name() ) + "\".  The type of the first operand is " +
                        frantic::strings::to_string(
                            channels::channel_data_type_str( i->second.arity(), i->second.data_type() ) ) +
                        ", while the type of the second operand is " +
                        frantic::strings::to_string(
                            channels::channel_data_type_str( j->second.arity(), j->second.data_type() ) ) +
                        "." );
            } else {
                add_channel( i->second.name(), i->second.arity(), i->second.data_type() );
            }
        }
    }

    // This is a zero element used for interpolation when a channel doesn't have an input
    vector<char> zero( 16 );

    // Set up all the accessors for filling in the blended channel data.
    vector<rle_channel_general_accessor> outputChannels;
    vector<const_rle_channel_general_accessor> inputChannelsFirst, inputChannelsSecond;
    // Also some flags, bit 0 indicates that the first operand has the channel, bit 1 indicates that the second operand
    // does.
    vector<char> inputChannelFlags;
    for( map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.begin(), ie = m_namedChannels.end();
         i != ie; ++i ) {
        outputChannels.push_back( get_channel_general_accessor( i->second.name() ) );
        if( onlyChannelsFromFirst ) {
            inputChannelsFirst.push_back( rleFirst.get_channel_general_accessor( i->second.name() ) );
        } else if( onlyChannelsFromSecond ) {
            inputChannelsSecond.push_back( rleSecond.get_channel_general_accessor( i->second.name() ) );
        } else {
            // Make sure that the zero vector is big enough
            if( zero.size() < i->second.primitive_size() )
                zero.resize( i->second.primitive_size() );
            char hasChannels = 0;
            if( rleFirst.has_channel( i->second.name() ) ) {
                inputChannelsFirst.push_back( rleFirst.get_channel_general_accessor( i->second.name() ) );
                hasChannels |= 1;
            } else {
                inputChannelsFirst.push_back( const_rle_channel_general_accessor() );
            }
            if( rleSecond.has_channel( i->second.name() ) ) {
                inputChannelsSecond.push_back( rleSecond.get_channel_general_accessor( i->second.name() ) );
                hasChannels |= 2;
            } else {
                inputChannelsSecond.push_back( const_rle_channel_general_accessor() );
            }
            inputChannelFlags.push_back( hasChannels );
        }
    }

    // Get a pointer to the zero data for speed
    char* zeroVectorPointer = &zero[0];

    ////////////
    // Iterate over the BCIndex values and combine the input scanlines
    ////////////

    // Allocate an array of blending alphas for combining the named channels
    vector<float> blendingAlpha;

    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0, csize = ris.m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = ris.m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int y = ris.m_abcCoordOrigin.y + b, z = ris.m_abcCoordOrigin.z + c;

            const rle_index_spec& risA = rleFirst.m_rleIndex;
            const rle_index_spec& risB = rleSecond.m_rleIndex;

            int bcIndex = b + c * bsize;
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            // The default state is an undefined region of the exterior region code
            bool creatingUndefinedRun = true;
            boost::int32_t creatingRegionCode = ris.m_exteriorRegionCode;

            // Loop through all the common sub-intervals of these scanlines
            rle_pairwise_run_iterator i( risA, risA.y_to_b( y ), risA.z_to_c( z ), risB, risB.y_to_b( y ),
                                         risB.z_to_c( z ) ),
                ie;
            for( ; i != ie; ++i ) {
                if( i.get_first_data_index() >= 0 && i.get_second_data_index() >= 0 ) {
                    // In this case, both inputs are defined, so we have to compare the data voxel-by-voxel
                    if( creatingUndefinedRun ) {
                        ris.m_runData.push_back( run_data( i.get_xmin(), (int)m_distanceData.size() ) );
                        creatingUndefinedRun = false;
                    }
                    int dataIndexFirst = i.get_first_data_index();
                    int dataIndexSecond = i.get_second_data_index();
                    if( !outputChannels.empty() ) {
                        // Resize the blending alpha array
                        blendingAlpha.resize( 2 * i.get_xsize() );
                        // Combine the distance arrays, and at the same time compute the alpha values for blending
                        // channels
                        float alpha = 0;
                        for( int d = 0, de = i.get_xsize(); d != de; ++d ) {
                            float distanceFirst = rleFirst[dataIndexFirst + d],
                                  distanceSecond = rleSecond[dataIndexSecond + d];
                            m_distanceData.push_back( CSGPolicy::combine_distances_alpha(
                                distanceFirst, distanceSecond, channelBlendWorldDistance, alpha ) );
                            // Store both alpha and 1-alpha in adjacent values, for use later with the named channel
                            // blending functions
                            blendingAlpha[2 * d] = 1.f - alpha;
                            blendingAlpha[2 * d + 1] = alpha;
                        }

                        // Blend the data for all the other channels
                        for( size_t channelIndex = 0, channelIndexEnd = outputChannels.size();
                             channelIndex != channelIndexEnd; ++channelIndex ) {

                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            size_t primitiveSize = channel.primitive_size();
                            size_t dataSize = primitiveSize * i.get_xsize();
                            char* newElements = channel.m_data->add_element( dataSize );

                            if( onlyChannelsFromFirst ) {
                                const_rle_channel_general_accessor& inputChannel = inputChannelsFirst[channelIndex];
                                memcpy( newElements, inputChannel.data( dataIndexFirst ), dataSize );
                            } else if( onlyChannelsFromSecond ) {
                                const_rle_channel_general_accessor& inputChannel = inputChannelsSecond[channelIndex];
                                memcpy( newElements, inputChannel.data( dataIndexSecond ), dataSize );
                            } else {
                                // Weighted sum style blend
                                int alphaIndex = 0;
                                channel_weighted_sum_combine_function_t ws = channel.m_weightedSumFunction;
                                size_t arity = channel.arity();
                                const char* data[2];

                                switch( inputChannelFlags[channelIndex] ) {
                                // Both the inputs are there
                                case 3: {
                                    const_rle_channel_general_accessor& inputChannelFirst =
                                        inputChannelsFirst[channelIndex];
                                    const_rle_channel_general_accessor& inputChannelSecond =
                                        inputChannelsSecond[channelIndex];
                                    data[0] = inputChannelFirst.data( dataIndexFirst );
                                    data[1] = inputChannelSecond.data( dataIndexSecond );
                                    // This writes out the weighted sum of the two inputs to the output channel
                                    for( int d = 0, de = i.get_xsize(); d != de; ++d ) {
                                        ws( &blendingAlpha[alphaIndex], data, 2, arity, newElements );
                                        alphaIndex += 2;
                                        data[0] += primitiveSize;
                                        data[1] += primitiveSize;
                                        newElements += primitiveSize;
                                    }
                                    break;
                                }
                                // Just the first input is there
                                case 1: {
                                    const_rle_channel_general_accessor& inputChannelFirst =
                                        inputChannelsFirst[channelIndex];
                                    data[0] = inputChannelFirst.data( dataIndexFirst );
                                    data[1] = zeroVectorPointer;
                                    // This writes out the weighted sum of the first input with zero to the output
                                    // channel
                                    for( int d = 0, de = i.get_xsize(); d != de; ++d ) {
                                        ws( &blendingAlpha[alphaIndex], data, 2, arity, newElements );
                                        alphaIndex += 2;
                                        data[0] += primitiveSize;
                                        newElements += primitiveSize;
                                    }
                                    break;
                                }
                                // Just the second input is there
                                case 2: {
                                    const_rle_channel_general_accessor& inputChannelSecond =
                                        inputChannelsSecond[channelIndex];
                                    data[0] = zeroVectorPointer;
                                    data[1] = inputChannelSecond.data( dataIndexSecond );
                                    // This writes out the weighted sum of the first input with zero to the output
                                    // channel
                                    for( int d = 0, de = i.get_xsize(); d != de; ++d ) {
                                        ws( &blendingAlpha[alphaIndex], data, 2, arity, newElements );
                                        alphaIndex += 2;
                                        data[1] += primitiveSize;
                                        newElements += primitiveSize;
                                    }
                                    break;
                                }
                                }
                            }
                        }
                    } else {
                        // If there are no named channels, then do a simpler loop which doesn't compute the blending
                        // alpha values
                        for( int d = 0, de = i.get_xsize(); d != de; ++d )
                            m_distanceData.push_back( CSGPolicy::combine_distances( rleFirst[dataIndexFirst + d],
                                                                                    rleSecond[dataIndexSecond + d] ) );
                    }
                } else if( i.get_first_data_index() >= 0 ) {
                    if( CSGPolicy::second_region_code_dominates( i.get_second_data_index() ) ) {
                        // In this case, the region code from the second operand dominates, so the run becomes
                        // undefined.
                        boost::int32_t regionCode = i.get_second_data_index();
                        // If necessary, complement the region code
                        if( CSGPolicy::second_op_complemented() ) {
                            regionCode = ( regionCode < -1 ) ? -1 : -2;
                        }
                        if( creatingUndefinedRun ) {
                            if( creatingRegionCode != regionCode ) {
                                ris.m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                                creatingRegionCode = regionCode;
                            }
                        } else {
                            ris.m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                            creatingUndefinedRun = true;
                            creatingRegionCode = regionCode;
                        }
                    } else {
                        // In this case, the data from the first operand is defined and dominates over the second
                        // operand region code, so we can copy the data from the first operand.
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back( run_data( i.get_xmin(), (int)m_distanceData.size() ) );
                            creatingUndefinedRun = false;
                        }
                        int dataIndex = i.get_first_data_index();
                        // Copy the data for all the other channels
                        for( size_t channelIndex = 0, channelIndexEnd = outputChannels.size();
                             channelIndex != channelIndexEnd; ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            size_t dataSize = i.get_xsize() * channel.primitive_size();
                            bool copyChannel;
                            if( onlyChannelsFromFirst || onlyChannelsFromSecond )
                                copyChannel = useChannelsFromFirst;
                            else
                                copyChannel = ( inputChannelFlags[channelIndex] & 1 ) != 0;
                            if( copyChannel ) {
                                const_rle_channel_general_accessor& inputChannel = inputChannelsFirst[channelIndex];
                                memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( dataIndex ),
                                        dataSize );
                            } else {
                                memset( channel.m_data->add_element( dataSize ), 0, dataSize );
                            }
                        }
                        if( CSGPolicy::first_op_complemented() ) {
                            for( int d = 0, de = i.get_xsize(); d != de; ++d )
                                m_distanceData.push_back( -rleFirst[dataIndex + d] );
                        } else {
                            for( int d = 0, de = i.get_xsize(); d != de; ++d )
                                m_distanceData.push_back( rleFirst[dataIndex + d] );
                        }
                    }
                } else if( i.get_second_data_index() >= 0 ) {
                    if( CSGPolicy::first_region_code_dominates( i.get_first_data_index() ) ) {
                        // In this case, the region code from the first operand dominates, so the run becomes undefined.
                        boost::int32_t regionCode = i.get_first_data_index();
                        // If necessary, complement the region code
                        if( CSGPolicy::first_op_complemented() ) {
                            regionCode = ( regionCode < -1 ) ? -1 : -2;
                        }
                        if( creatingUndefinedRun ) {
                            if( creatingRegionCode != regionCode ) {
                                ris.m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                                creatingRegionCode = regionCode;
                            }
                        } else {
                            ris.m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                            creatingUndefinedRun = true;
                            creatingRegionCode = regionCode;
                        }
                    } else {
                        // In this case, the data from the second operand is defined and dominates over the first
                        // operand region code, so we can copy the data from the second operand
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back( run_data( i.get_xmin(), (int)m_distanceData.size() ) );
                            creatingUndefinedRun = false;
                        }
                        int dataIndex = i.get_second_data_index();
                        // Copy the data for all the other channels
                        for( size_t channelIndex = 0, channelIndexEnd = outputChannels.size();
                             channelIndex != channelIndexEnd; ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            size_t dataSize = i.get_xsize() * channel.primitive_size();
                            bool copyChannel;
                            if( onlyChannelsFromFirst || onlyChannelsFromSecond )
                                copyChannel = useChannelsFromSecond;
                            else
                                copyChannel = ( inputChannelFlags[channelIndex] & 2 ) != 0;
                            if( copyChannel ) {
                                const_rle_channel_general_accessor& inputChannel = inputChannelsSecond[channelIndex];
                                memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( dataIndex ),
                                        dataSize );
                            } else {
                                memset( channel.m_data->add_element( dataSize ), 0, dataSize );
                            }
                        }
                        if( CSGPolicy::second_op_complemented() ) {
                            for( int d = 0, de = i.get_xsize(); d != de; ++d )
                                m_distanceData.push_back( -rleSecond[dataIndex + d] );
                        } else {
                            for( int d = 0, de = i.get_xsize(); d != de; ++d )
                                m_distanceData.push_back( rleSecond[dataIndex + d] );
                        }
                    }
                } else {
                    // In this case, both region codes are undefined
                    boost::int32_t regionCode =
                        CSGPolicy::combine_region_codes( i.get_first_data_index(), i.get_second_data_index() );
                    if( creatingUndefinedRun ) {
                        if( creatingRegionCode != regionCode ) {
                            ris.m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                            creatingRegionCode = regionCode;
                        }
                    } else {
                        ris.m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                        creatingUndefinedRun = true;
                        creatingRegionCode = regionCode;
                    }
                }
            }

            // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to optimize the
            // storage and access. If no runs were added, add a zero-sized run.
            if( (int)ris.m_runData.size() > ris.m_bcToRunIndex[bcIndex] ) {
                // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to optimize
                // the storage and access.
                if( !creatingUndefinedRun || creatingRegionCode != ris.m_exteriorRegionCode ) {
                    // When the rle_pairwise_run_iterator loop is done, i.get_xmin() is the value 1 past
                    // the end of the last sub-interval, so it is the correct value to close of the m_runData array.
                    ris.m_runData.push_back( run_data( i.get_xmin(), -1 ) );
                }
            } else {
                ris.m_runData.push_back( run_data( 0, -1 ) );
                ris.m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();

    // Set the data size as well
    ris.m_dataSize = m_distanceData.size();
}

namespace detail {
struct union_csg_policy {
    /**
     * For union, inside always beats outside.
     */
    static boost::int32_t combine_region_codes( boost::int32_t regionCodeFirst, boost::int32_t regionCodeSecond ) {
        if( regionCodeFirst == -1 && regionCodeSecond == -1 )
            return -1;
        else
            return -2;
    }

    /**
     * Returns whether the region codes of the operands dominate.
     */
    static bool first_region_code_dominates( boost::int32_t regionCode ) { return regionCode < -1; }
    static bool second_region_code_dominates( boost::int32_t regionCode ) { return regionCode < -1; }

    /**
     * Whether the operands are complemented.
     */
    static bool first_op_complemented() { return false; }
    static bool second_op_complemented() { return false; }

    /**
     * Union combines using the min operation.
     */
    static float combine_distances( float distanceFirst, float distanceSecond ) {
        return min( distanceFirst, distanceSecond );
    }

    /**
     * Union combines using the min operation, and this gives back an alpha for blending other channels.
     */
    static float combine_distances_alpha( float distanceFirst, float distanceSecond, float channelBlendWorldDistance,
                                          float& outAlpha ) {
        outAlpha = 0.5f + ( distanceFirst - distanceSecond ) / channelBlendWorldDistance;
        if( distanceFirst < distanceSecond ) {
            if( outAlpha < 0 )
                outAlpha = 0;
            return distanceFirst;
        } else {
            if( outAlpha > 1 )
                outAlpha = 1;
            return distanceSecond;
        }
    }

    /**
     * Union levelsets have both the first and second operand's channels
     */
    static bool use_first_op_channels() { return true; }

    /**
     * Union levelsets have both the first and second operand's channels
     */
    static bool use_second_op_channels() { return true; }
};

struct intersect_csg_policy {
    /**
     * For intersection, outside always beats inside.
     */
    static boost::int32_t combine_region_codes( boost::int32_t regionCodeFirst, boost::int32_t regionCodeSecond ) {
        if( regionCodeFirst < -1 && regionCodeSecond < -1 )
            return -2;
        else
            return -1;
    }

    /**
     * Returns whether the region codes of the operands dominate.
     */
    static bool first_region_code_dominates( boost::int32_t regionCode ) { return regionCode == -1; }
    static bool second_region_code_dominates( boost::int32_t regionCode ) { return regionCode == -1; }

    /**
     * Whether the operands are complemented.
     */
    static bool first_op_complemented() { return false; }
    static bool second_op_complemented() { return false; }

    /**
     * Intersection combines using the max operation.
     */
    static float combine_distances( float distanceFirst, float distanceSecond ) {
        return max( distanceFirst, distanceSecond );
    }

    /**
     * Intersection combines using the max operation, and this gives back an alpha for blending other channels.
     */
    static float combine_distances_alpha( float distanceFirst, float distanceSecond, float channelBlendWorldDistance,
                                          float& outAlpha ) {
        outAlpha = 0.5f + ( distanceSecond - distanceFirst ) / channelBlendWorldDistance;
        if( distanceFirst < distanceSecond ) {
            if( outAlpha > 1 )
                outAlpha = 1;
            return distanceSecond;
        } else {
            if( outAlpha < 0 )
                outAlpha = 0;
            return distanceFirst;
        }
    }

    /**
     * Intersected levelsets have both the first and second operand's channels
     */
    static bool use_first_op_channels() { return true; }

    /**
     * Intersected levelsets have both the first and second operand's channels
     */
    static bool use_second_op_channels() { return true; }
};

struct subtract_csg_policy {
    /**
     * Outside for the first operand and inside for the second operand dominate, because it is
     * A intersect ~B.
     */
    static boost::int32_t combine_region_codes( boost::int32_t regionCodeFirst, boost::int32_t regionCodeSecond ) {
        if( regionCodeFirst < -1 && regionCodeSecond == -1 )
            return -2;
        else
            return -1;
    }

    /**
     * Returns whether the region codes of the operands dominate.
     */
    static bool first_region_code_dominates( boost::int32_t regionCode ) { return regionCode == -1; }
    static bool second_region_code_dominates( boost::int32_t regionCode ) { return regionCode < -1; }

    /**
     * Whether the operands are complemented.
     */
    static bool first_op_complemented() { return false; }
    static bool second_op_complemented() { return true; }

    /**
     * Intersection combines using the max operation.
     */
    static float combine_distances( float distanceFirst, float distanceSecond ) {
        return max( distanceFirst, -distanceSecond );
    }

    /**
     * Intersection combines using the max operation, and this gives back an alpha for blending other channels.
     */
    static float combine_distances_alpha( float distanceFirst, float distanceSecond, float channelBlendWorldDistance,
                                          float& outAlpha ) {
        distanceSecond = -distanceSecond;

        outAlpha = 0.5f + ( distanceSecond - distanceFirst ) / channelBlendWorldDistance;
        if( distanceFirst < distanceSecond ) {
            if( outAlpha > 1 )
                outAlpha = 1;
            return distanceSecond;
        } else {
            if( outAlpha < 0 )
                outAlpha = 0;
            return distanceFirst;
        }
    }

    /**
     * Subtract uses only the first operand's channels
     */
    static bool use_first_op_channels() { return true; }

    /**
     * Subtract does not use the second operand's channels
     */
    static bool use_second_op_channels() { return false; }
};
} // namespace detail

void rle_level_set::csg_union( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                               float channelBlendDistance ) {
    csg_operation<detail::union_csg_policy>( rleFirst, rleSecond, channelBlendDistance );
}

void rle_level_set::csg_intersect( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                                   float channelBlendDistance ) {
    csg_operation<detail::intersect_csg_policy>( rleFirst, rleSecond, channelBlendDistance );
}

void rle_level_set::csg_subtract( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                                  float channelBlendDistance ) {
    csg_operation<detail::subtract_csg_policy>( rleFirst, rleSecond, channelBlendDistance );
}

void rle_level_set::linear_interpolate( const rle_level_set& rleFirst, const rle_level_set& rleSecond, float alpha ) {
    // check that alpha is between 0 and 1
    if( alpha > 1 || alpha < 0 ) {
        throw std::runtime_error(
            "rle_level_set.linear_interpolate: alpha should be between 0 and 1 (inclusive), but alpha is " +
            boost::lexical_cast<string>( alpha ) + ".\n" );
    }

    // check that the input level sets have the same voxel coordinate system
    if( rleFirst.get_voxel_coord_system() != rleSecond.get_voxel_coord_system() ) {
        throw std::runtime_error( "rle_level_set.linear_interpolate: the voxel coordinate systems of the two level set "
                                  "operands don't match.  The first is " +
                                  rleFirst.m_voxelCoordSystem.str() + ", while the second is " +
                                  rleSecond.m_voxelCoordSystem.str() + "." );
    }

    // Make copies of the input level sets
    rle_level_set copyFirst( rleFirst ), copySecond( rleSecond );
    // Create the rle_index_spec for blending
    rle_index_spec ris;
    ris.combine_for_blend( rleFirst.get_rle_index_spec(), rleSecond.get_rle_index_spec() );
    // Switch both level set copies to use the new rle_index_spec, creating a populated channel for them
    copyFirst.switch_rle_index_spec( ris, _T("tempPreSwitch") );
    copySecond.switch_rle_index_spec_with_swap( ris, _T("tempPreSwitch") );

    // Duplicate the populated channels so we can use them for extrapolation after the reinitialization
    copyFirst.duplicate_channel( _T("tempPreSwitchExtrap"), _T("tempPreSwitch") );
    copySecond.duplicate_channel( _T("tempPreSwitchExtrap"), _T("tempPreSwitch") );
    // Reinitialize both level sets from the previously defined voxels to the new ones
    copyFirst.reinitialize_signed_distance_from_populated( _T("tempPreSwitch") );
    copySecond.reinitialize_signed_distance_from_populated( _T("tempPreSwitch") );
    // Erase the populated channel used for reinitialization
    copyFirst.erase_channel( _T("tempPreSwitch") );
    copySecond.erase_channel( _T("tempPreSwitch") );

    // Extrapolate all the channels into the newly defined voxels
    copyFirst.extrapolate_channels( _T("tempPreSwitchExtrap") );
    copySecond.extrapolate_channels( _T("tempPreSwitchExtrap") );
    // Erase the populated channel used for extrapolation
    copyFirst.erase_channel( _T("tempPreSwitchExtrap") );
    copySecond.erase_channel( _T("tempPreSwitchExtrap") );

    // Get all the channel names from the first level set
    vector<frantic::tstring> channelNames;
    copySecond.get_channel_names( channelNames );

    vector<rle_channel_general_accessor> firstAccessors, secondAccessors;

    for( size_t i = 0, ie = channelNames.size(); i != ie; ++i ) {
        // If a channel exists in both, we interpolate, otherwise we use the data without any blending
        if( copyFirst.has_channel( channelNames[i] ) ) {
            firstAccessors.push_back( copyFirst.get_channel_general_accessor( channelNames[i] ) );
            secondAccessors.push_back( copySecond.get_channel_general_accessor( channelNames[i] ) );
            if( firstAccessors.back().arity() != secondAccessors.back().arity() )
                throw runtime_error( "rle_level_set::linear_interpolate() - Channel \"" +
                                     frantic::strings::to_string( channelNames[i] ) +
                                     "\" has different arity in the two input level sets." );
            if( firstAccessors.back().data_type() != secondAccessors.back().data_type() )
                throw runtime_error( "rle_level_set::linear_interpolate() - Channel \"" +
                                     frantic::strings::to_string( channelNames[i] ) +
                                     "\" has different data types in the two input level sets." );
        } else {
            rle_channel_general_accessor channelSecond = copySecond.get_channel_general_accessor( channelNames[i] );
            // Create the channel in the first level set
            copyFirst.add_channel( channelNames[i], channelSecond.arity(), channelSecond.data_type() );
            rle_channel_general_accessor channelFirst = copyFirst.get_channel_general_accessor( channelNames[i] );
            // Copy all the data into the newly created channel
            memcpy( channelFirst.data( 0 ), channelSecond.data( 0 ),
                    channelFirst.primitive_size() * channelFirst.size() );
        }
    }

    float weights[2] = { 1.f - alpha, alpha };
    char* data[2];
    // Go through all the channels to blend
    for( size_t j = 0, je = firstAccessors.size(); j != je; ++j ) {
        channel_weighted_sum_combine_function_t combineFunction = firstAccessors[j].get_weighted_sum_combine_function();
        data[0] = firstAccessors[j].data( 0 );
        data[1] = secondAccessors[j].data( 0 );
        //		size_t primitiveSize = firstAccessors[j].primitive_size();
        size_t arity = firstAccessors[j].arity();
        // Go through all the defined voxels of this channel
        combineFunction( weights, data, 2, arity * m_rleIndex.data_size(), data[0] );
    }

    // don't forget about the distance data either
    for( int i = 0, ie = static_cast<int>( copyFirst.m_distanceData.size() ); i != ie; ++i ) {
        copyFirst.m_distanceData[i] = ( 1 - alpha ) * copyFirst[i] + (alpha)*copySecond[i];
    }

    copyFirst.swap( *this );
}

void rle_level_set::linear_interpolate( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                                        const frantic::tstring& blendAlphaChannel ) {

    // check that the input level sets have the same voxel coordinate system
    if( rleFirst.get_voxel_coord_system() != rleSecond.get_voxel_coord_system() ) {
        throw std::runtime_error( "rle_level_set.linear_interpolate: the voxel coordinate systems of the two level set "
                                  "operands don't match.  The first is " +
                                  rleFirst.m_voxelCoordSystem.str() + ", while the second is " +
                                  rleSecond.m_voxelCoordSystem.str() + "." );
    }

    // check that alpha is between 0 and 1
    if( !rleFirst.has_channel( blendAlphaChannel ) ) {
        throw std::runtime_error(
            "rle_level_set.linear_interpolate: the first level set does not have the blend alpha channel " +
            frantic::strings::to_string( blendAlphaChannel ) + ".\n" );
    }

    // Make copies of the input level sets
    rle_level_set copyFirst( rleFirst ), copySecond( rleSecond );
    // Create the rle_index_spec for blending
    rle_index_spec ris;
    ris.combine_for_blend( rleFirst.get_rle_index_spec(), rleSecond.get_rle_index_spec() );
    // Switch both level set copies to use the new rle_index_spec, creating a populated channel for them
    copyFirst.switch_rle_index_spec( ris, _T("tempPreSwitch") );
    copySecond.switch_rle_index_spec_with_swap( ris, _T("tempPreSwitch") );

    // Duplicate the populated channels so we can use them for extrapolation after the reinitialization
    copyFirst.duplicate_channel( _T("tempPreSwitchExtrap"), _T("tempPreSwitch") );
    copySecond.duplicate_channel( _T("tempPreSwitchExtrap"), _T("tempPreSwitch") );
    // Reinitialize both level sets from the previously defined voxels to the new ones
    copyFirst.reinitialize_signed_distance_from_populated( _T("tempPreSwitch") );
    copySecond.reinitialize_signed_distance_from_populated( _T("tempPreSwitch") );
    // Erase the populated channel used for reinitialization
    copyFirst.erase_channel( _T("tempPreSwitch") );
    copySecond.erase_channel( _T("tempPreSwitch") );

    // Extrapolate all the channels into the newly defined voxels
    copyFirst.extrapolate_channels( _T("tempPreSwitchExtrap") );
    copySecond.extrapolate_channels( _T("tempPreSwitchExtrap") );
    // Erase the populated channel used for extrapolation
    copyFirst.erase_channel( _T("tempPreSwitchExtrap") );
    copySecond.erase_channel( _T("tempPreSwitchExtrap") );

    const_rle_channel_accessor<float> blendAlpha = copyFirst.get_channel_accessor<float>( blendAlphaChannel );

    // Get all the channel names from the first level set
    vector<frantic::tstring> channelNames;
    copySecond.get_channel_names( channelNames );

    vector<rle_channel_general_accessor> firstAccessors, secondAccessors;

    for( size_t i = 0, ie = channelNames.size(); i != ie; ++i ) {
        // If a channel exists in both, we interpolate, otherwise we use the data without any blending
        if( copyFirst.has_channel( channelNames[i] ) ) {
            firstAccessors.push_back( copyFirst.get_channel_general_accessor( channelNames[i] ) );
            secondAccessors.push_back( copySecond.get_channel_general_accessor( channelNames[i] ) );
            if( firstAccessors.back().arity() != secondAccessors.back().arity() )
                throw runtime_error( "rle_level_set::linear_interpolate() - Channel \"" +
                                     frantic::strings::to_string( channelNames[i] ) +
                                     "\" has different arity in the two input level sets." );
            if( firstAccessors.back().data_type() != secondAccessors.back().data_type() )
                throw runtime_error( "rle_level_set::linear_interpolate() - Channel \"" +
                                     frantic::strings::to_string( channelNames[i] ) +
                                     "\" has different data types in the two input level sets." );
        } else {
            rle_channel_general_accessor channelSecond = copySecond.get_channel_general_accessor( channelNames[i] );
            // Create the channel in the first level set
            copyFirst.add_channel( channelNames[i], channelSecond.arity(), channelSecond.data_type() );
            rle_channel_general_accessor channelFirst = copyFirst.get_channel_general_accessor( channelNames[i] );
            // Copy all the data into the newly created channel
            memcpy( channelFirst.data( 0 ), channelSecond.data( 0 ),
                    channelFirst.primitive_size() * channelFirst.size() );
        }
    }

    char* data[2];
    // Go through all the channels to blend
    for( size_t j = 0, je = firstAccessors.size(); j != je; ++j ) {
        channel_weighted_sum_combine_function_t combineFunction = firstAccessors[j].get_weighted_sum_combine_function();
        data[0] = firstAccessors[j].data( 0 );
        data[1] = secondAccessors[j].data( 0 );

        size_t arity = firstAccessors[j].arity();

        for( boost::int32_t dataIndex = 0; dataIndex < (boost::int32_t)m_rleIndex.data_size(); ++dataIndex ) {
            float alpha = blendAlpha[dataIndex];
            // we need to get the weights from the blendAlpha channel for each entry.
            float weights[2] = { 1.f - alpha, alpha };

            // Go through all the defined voxels of this channel
            combineFunction( weights, data, 2, arity, data[0] );
        }
    }

    // don't forget about the distance data either
    for( int i = 0, ie = static_cast<int>( copyFirst.m_distanceData.size() ); i != ie; ++i ) {
        float alpha = blendAlpha[i];
        copyFirst.m_distanceData[i] = ( 1 - alpha ) * copyFirst[i] + (alpha)*copySecond[i];
    }

    copyFirst.swap( *this );
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
