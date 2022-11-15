// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/level_set_fast_marching_reinitialization.hpp>
#include <frantic/volumetrics/levelset/level_set_marching_extrapolation.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::geometry;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;
using namespace frantic::fluids; // for rle_voxel_field

namespace {
#if defined( _WIN32 ) || defined( _WIN64 )
inline int strcmp( const wchar_t* str1, const wchar_t* str2 ) { return wcscmp( str1, str2 ); }
#endif
} // namespace

namespace frantic {
namespace fluids {

void rle_voxel_field::swap_with_rle_level_set( frantic::volumetrics::levelset::rle_level_set& levelSet ) {
    if( levelSet.has_channel( _T("SignedDistance") ) )
        throw std::runtime_error(
            "rle_voxel_field::swap_with_rle_level_set: The original level set cannot have a channel "
            "named SignedDistance. This channel name is reserved." );

    // save a temporary copy of the voxel field's signed distance
    std::vector<float> voxelFieldTempSignedDistChannel( size(), levelSet.m_outsideDistance );
    if( has_channel( _T("SignedDistance") ) ) {
        rle_channel_general_accessor signedDistAcc = get_channel_general_accessor( _T("SignedDistance") );
        if( signedDistAcc.data_type() != data_type_float32 || signedDistAcc.arity() != 1 )
            throw std::runtime_error( "rle_voxel_field::swap_with_rle_level_set: Could not convert vector field "
                                      "SignedDistance channel to a level set. Data type or arity is not correct" );
        if( size() > 0 )
            memcpy( &voxelFieldTempSignedDistChannel[0], signedDistAcc.data( 0 ), size() * sizeof( float ) );
        erase_channel( _T("SignedDistance") );
    }

    // swap the internal data members
    m_rleIndex.swap( levelSet.m_rleIndex );
    m_voxelCoordSystem.swap( levelSet.m_voxelCoordSystem );
    m_namedChannels.swap( levelSet.m_namedChannels );

    // copy the level set's signed distance into the voxel field
    add_channel<float>( _T("SignedDistance") );
    rle_channel_accessor<float> signedDistAcc = get_channel_accessor<float>( _T("SignedDistance") );
    if( size() > 0 )
        memcpy( &signedDistAcc[0], &levelSet[0], size() * sizeof( float ) );

    // swap the voxel field's temp signed distance into the level set
    levelSet.m_distanceData.swap( voxelFieldTempSignedDistChannel );
}

void rle_voxel_field::duplicate_channel( const frantic::tstring& destChannelName,
                                         const frantic::tstring& sourceChannelName ) {
    std::map<frantic::tstring, rle_channel>::iterator sourceChannel = m_namedChannels.find( sourceChannelName );
    if( sourceChannel == m_namedChannels.end() )
        throw std::runtime_error( "rle_level_set.duplicate_channel() - Tried to duplicate channel \"" +
                                  frantic::strings::to_string( sourceChannelName ) +
                                  "\", but no channel of that name exists." );

    if( !channels::is_valid_channel_name( destChannelName ) )
        throw std::runtime_error( "rle_level_set.duplicate_channel() - Tried to create channel \"" +
                                  frantic::strings::to_string( destChannelName ) +
                                  "\".  The specified channel name is not valid for a channel name, channel names must "
                                  "start with a letter or "
                                  "underscore, and must contain only letters, numbers, and underscores." );

    // Create the destination channel with the appropriate data type
    std::map<frantic::tstring, rle_channel>::iterator destChannel = m_namedChannels.find( destChannelName );
    if( destChannel == m_namedChannels.end() ) {
        m_namedChannels.insert(
            std::make_pair( destChannelName, rle_channel( destChannelName, sourceChannel->second.arity(),
                                                          sourceChannel->second.data_type() ) ) );
        destChannel = m_namedChannels.find( destChannelName );
        // Make sure the vertex array count matches that of the rle_index_spec's data size
        destChannel->second.m_data.resize( m_rleIndex.data_size() * destChannel->second.primitive_size() );
    } else if( destChannel->second.data_type() != sourceChannel->second.data_type() ||
               destChannel->second.arity() != sourceChannel->second.arity() ) {
        destChannel->second.set( destChannelName, sourceChannel->second.arity(), sourceChannel->second.data_type() );
        // Make sure the vertex array count matches that of the rle_index_spec's data size
        destChannel->second.m_data.resize( m_rleIndex.data_size() * destChannel->second.primitive_size() );
    }

    // Copy the memory to the duplicated channel
    memcpy( destChannel->second.m_data.begin(), sourceChannel->second.m_data.begin(),
            m_rleIndex.data_size() * destChannel->second.primitive_size() );
}

void rle_voxel_field::add_channel( const frantic::tstring& channelName, std::size_t arity, data_type_t dataType ) {
    if( !channels::is_valid_channel_name( channelName ) )
        throw std::runtime_error(
            "rle_voxel_field.add_channel() - Tried to add channel \"" + frantic::strings::to_string( channelName ) +
            "\", with data type " + frantic::strings::to_string( channel_data_type_str( arity, dataType ) ) +
            ".  "
            "The specified channel name is not valid for a channel name, channel names must start "
            "with a letter or underscore, "
            "and must contain only letters, numbers, and underscores." );

    std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( channelName );
    if( i == m_namedChannels.end() ) {
        m_namedChannels.insert( std::make_pair( channelName, rle_channel( channelName, arity, dataType ) ) );
        i = m_namedChannels.find( channelName );
        // Make sure the voxel data array count matches that of the rle_index_spec's data size
        i->second.m_data.resize( m_rleIndex.data_size() * i->second.primitive_size() );
    } else {
        // If the existing channel doesn't match arity and data type, then we tweak it so it does
        if( i->second.arity() != arity || i->second.data_type() != dataType ) {
            i->second.set( channelName, arity, dataType );
            // Make sure the voxel data array count matches that of the rle_index_spec's data size
            i->second.m_data.resize( m_rleIndex.data_size() * i->second.primitive_size() );
        }
    }
}

void rle_voxel_field::zero_channel( const frantic::tstring& channelName ) {
    std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( channelName );
    if( i == m_namedChannels.end() ) {
        throw runtime_error( "rle_voxel_field.zero_channel() - Tried to set channel \"" +
                             frantic::strings::to_string( channelName ) +
                             "\" "
                             "to all zeros, but no such channel exists in the level set." );
    } else {
        // memset the buffer for this channel to all zeros
        size_t size = i->second.m_data.size();
        if( size > 0 )
            memset( i->second.m_data.begin(), 0, size );
    }
}

float rle_voxel_field::get_channel_max_norm( const frantic::tstring& channelName ) const {
    if( !has_channel( channelName ) )
        throw runtime_error( "rle_voxel_field.get_channel_max_norm() - Tried to get the maximum L2 norm of channel \"" +
                             frantic::strings::to_string( channelName ) +
                             "\", "
                             "but no such channel exists." );

    const_rle_channel_general_accessor chanAcc = get_channel_general_accessor( channelName );
    if( chanAcc.data_type() != channels::data_type_float32 )
        throw runtime_error( "rle_voxel_field.get_channel_max_norm() - Tried to get the maximum L2 norm of channel \"" +
                             frantic::strings::to_string( channelName ) +
                             "\", "
                             "but only float32[] data types are supported currently.  Its data type is " +
                             chanAcc.type_str() + "." );

    float maximumNorm = 0;
    for( size_t i = 0, ie = chanAcc.size(); i != ie; ++i ) {
        // Compute the sum of the squares of the components (L2 norm squared)
        float norm = 0;
        const float* data = reinterpret_cast<const float*>( chanAcc.data( i ) );
        for( size_t j = 0, je = chanAcc.arity(); j != je; ++j ) {
            norm += ( *data ) * ( *data );
            ++data;
        }
        if( norm > maximumNorm )
            maximumNorm = norm;
    }

    // Take the square root to return the norm
    return sqrtf( maximumNorm );
}

void rle_voxel_field::copy_channels( const rle_voxel_field& inputRVF,
                                     const std::vector<frantic::tstring>& channelsToCopy,
                                     const frantic::tstring& populatedChannelToCreate ) {
    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;

    // The channel copying data array
    vector<rle_index_spec_channel_copying_data> channelCopyData;

    // First make sure that all the necesary channels are defined
    for( size_t i = 0, ie = channelsToCopy.size(); i != ie; ++i ) {
        if( !inputRVF.has_channel( channelsToCopy[i] ) )
            throw runtime_error(
                "rle_voxel_field::copy_channels() - The input rle_voxel_field doesn't have a channel named \"" +
                frantic::strings::to_string( channelsToCopy[i] ) + "\", which was requested as input" );
        inputAccessors.push_back( inputRVF.get_channel_general_accessor( channelsToCopy[i] ) );

        // If the output already had this data channel, then try to reuse it.
        if( has_channel( channelsToCopy[i] ) ) {
            rle_channel_general_accessor outputAccessor = get_channel_general_accessor( channelsToCopy[i] );
            if( outputAccessor.data_type() == inputAccessors.back().data_type() &&
                outputAccessor.arity() == inputAccessors.back().arity() ) {
                outputAccessors.push_back( outputAccessor );
            } else {
                erase_channel( channelsToCopy[i] );
                add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
                outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
            }
        } else {
            add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
            outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
        }

        // Initialize the output channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = inputAccessors.back().data( 0 );
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = inputAccessors.back().primitive_size();
    }

    // If the caller requested that the 'Populated' channel be made, do so
    vector<char> populatedChannelInput( populatedChannelToCreate.empty() ? 0
                                                                         : inputRVF.get_rle_index_spec().data_size() );
    if( !populatedChannelToCreate.empty() ) {
        if( !populatedChannelInput.empty() ) {
            // Initialize the input 'Populated' data to all ones
            memset( &populatedChannelInput[0], 1, populatedChannelInput.size() );
        }
        // Create the 'Populated' channel, or reuse it if it already exists
        add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        outputAccessors.push_back( get_channel_general_accessor( populatedChannelToCreate ) );
        // Initialize the output 'Populated' channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = populatedChannelInput.empty() ? 0 : &populatedChannelInput[0];
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = 1;
    }

    // Use the rle_index_spec method to copy all the defined data.
    if( !channelCopyData.empty() )
        rle_index_spec::copy_data_channels( m_rleIndex, inputRVF.get_rle_index_spec(), &channelCopyData[0],
                                            &channelCopyData[0] + channelCopyData.size() );
}

void rle_voxel_field::copy_channels( const rle_level_set& inputRLS, const std::vector<frantic::tstring>& channelsToCopy,
                                     const frantic::tstring& populatedChannelToCreate ) {
    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;

    // The channel copying data array
    vector<rle_index_spec_channel_copying_data> channelCopyData;

    // First make sure that all the necesary channels are defined
    for( size_t i = 0, ie = channelsToCopy.size(); i != ie; ++i ) {
        if( !inputRLS.has_channel( channelsToCopy[i] ) )
            throw runtime_error(
                "rle_voxel_field::copy_channels() - The input rle_voxel_field doesn't have a channel named \"" +
                frantic::strings::to_string( channelsToCopy[i] ) + "\", which was requested as input" );
        inputAccessors.push_back( inputRLS.get_channel_general_accessor( channelsToCopy[i] ) );

        // If the output already had this data channel, then try to reuse it.
        if( has_channel( channelsToCopy[i] ) ) {
            rle_channel_general_accessor outputAccessor = get_channel_general_accessor( channelsToCopy[i] );
            if( outputAccessor.data_type() == inputAccessors.back().data_type() &&
                outputAccessor.arity() == inputAccessors.back().arity() ) {
                outputAccessors.push_back( outputAccessor );
            } else {
                erase_channel( channelsToCopy[i] );
                add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
                outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
            }
        } else {
            add_channel( channelsToCopy[i], inputAccessors.back().arity(), inputAccessors.back().data_type() );
            outputAccessors.push_back( get_channel_general_accessor( channelsToCopy[i] ) );
        }

        // Initialize the output channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = inputAccessors.back().data( 0 );
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = inputAccessors.back().primitive_size();
    }

    // If the caller requested that the 'Populated' channel be made, do so
    vector<char> populatedChannelInput( populatedChannelToCreate.empty() ? 0
                                                                         : inputRLS.get_rle_index_spec().data_size() );
    if( !populatedChannelToCreate.empty() ) {
        if( !populatedChannelInput.empty() ) {
            // Initialize the input 'Populated' data to all ones
            memset( &populatedChannelInput[0], 1, populatedChannelInput.size() );
        }
        // Create the 'Populated' channel, or reuse it if it already exists
        add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        outputAccessors.push_back( get_channel_general_accessor( populatedChannelToCreate ) );
        // Initialize the output 'Populated' channel to all zeros
        memset( outputAccessors.back().data( 0 ), 0,
                outputAccessors.back().primitive_size() * outputAccessors.back().size() );

        // Add an entry to the channel copying data array
        channelCopyData.push_back( rle_index_spec_channel_copying_data() );
        rle_index_spec_channel_copying_data& cd = channelCopyData.back();
        cd.inputData = populatedChannelInput.empty() ? 0 : &populatedChannelInput[0];
        cd.outputData = outputAccessors.back().data( 0 );
        cd.primitiveSize = 1;
    }

    // Use the rle_index_spec method to copy all the defined data.
    if( !channelCopyData.empty() )
        rle_index_spec::copy_data_channels( m_rleIndex, inputRLS.get_rle_index_spec(), &channelCopyData[0],
                                            &channelCopyData[0] + channelCopyData.size() );
}

bool rle_voxel_field::check_consistency( std::ostream& out ) const {
    // check the rle index spec
    if( !m_rleIndex.check_consistency( out ) )
        return false;

    for( std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i =
             m_namedChannels.begin();
         i != m_namedChannels.end(); ++i ) {
        // check the channel
        if( !i->second.check_consistency( out ) )
            return false;

        bool pass = true;

        // check the channel against the spec
        if( i->second.size() != m_rleIndex.data_size() ) {
            out << "Channel " << frantic::strings::to_string( i->first ) << " has size " << i->second.size()
                << " which does not match index spec " << m_rleIndex.data_size() << std::endl;
            pass = false;
        }

        if( i->second.m_data.size() != i->second.m_primitiveSize * m_rleIndex.data_size() ) {
            out << "Channel " << frantic::strings::to_string( i->first ) << " has size " << i->second.m_data.size()
                << " which is not the same as " << m_rleIndex.data_size() << " elements of primitive size "
                << i->second.m_primitiveSize << std::endl;
            pass = false;
        }

        if( i->second.m_data.capacity() < i->second.m_primitiveSize * m_rleIndex.data_size() ) {
            out << "Channel " << frantic::strings::to_string( i->first ) << " has capaciy "
                << i->second.m_data.capacity() << " which cannot contain " << m_rleIndex.data_size()
                << " elements of primitive size " << i->second.m_primitiveSize << std::endl;
            pass = false;
        }

        if( !pass )
            return false;
    }

    return true;
}

// NOTE: this function was hijacked from the rle_level_set class
void rle_voxel_field::linear_interpolate( const rle_voxel_field& rleFirst, const rle_voxel_field& rleSecond,
                                          float alpha ) {
    // TODO: Refactor this function like crazy

    rle_index_spec& ris = m_rleIndex;

    // check that alpha is between 0 and 1
    if( alpha > 1 || alpha < 0 ) {
        throw std::runtime_error(
            "rle_voxel_field.linear_interpolate: alpha should be between 0 and 1 (inclusive), but alpha is " +
            boost::lexical_cast<string>( alpha ) + ".\n" );
    }

    // check that the input level sets have the same voxel coordinate system
    if( rleFirst.get_voxel_coord_system() != rleSecond.get_voxel_coord_system() ) {
        throw std::runtime_error(
            "rle_voxel_field.linear_interpolate: the voxel coordinate systems of the two level set "
            "operands don't match.  The first is " +
            rleFirst.m_voxelCoordSystem.str() + ", while the second is " + rleSecond.m_voxelCoordSystem.str() + "." );
    }

    // set the voxel coordinate system
    m_voxelCoordSystem = rleFirst.m_voxelCoordSystem;

    // set up the exterior region code
    if( rleFirst.m_rleIndex.m_exteriorRegionCode == rleSecond.m_rleIndex.m_exteriorRegionCode ) {
        ris.m_exteriorRegionCode = rleFirst.m_rleIndex.m_exteriorRegionCode;
    } else {
        // If the exterior region codes don't match, we will get weirdness, so just disallow this case.
        if( ris.m_exteriorRegionCode == -1 ) {
            throw std::runtime_error(
                "rle_voxel_field.linear_interpolate: The two exterior region codes do not match. The "
                "first operand is interior and the second is exterior." );
        } else {
            throw std::runtime_error(
                "rle_voxel_field.linear_interpolate: The two exterior region codes do not match. The "
                "first operand is exterior and the second is interior." );
        }
    }

    // set the bounding box
    ris.m_abcCoordOrigin =
        vector3::from_min( rleFirst.m_rleIndex.m_abcCoordOrigin, rleSecond.m_rleIndex.m_abcCoordOrigin );
    ris.m_abcCoordSize = size3::from_bounds(
        ris.m_abcCoordOrigin,
        vector3::from_max( rleFirst.m_rleIndex.m_abcCoordOrigin + rleFirst.m_rleIndex.m_abcCoordSize,
                           rleSecond.m_rleIndex.m_abcCoordOrigin + rleSecond.m_rleIndex.m_abcCoordSize ) -
            vector3( 1 ) );

    // clear any data already stored in ris
    ris.m_runData.clear();
    size_t totalDataSize = 0;

    // resize the bcToRunIndex vector to match our new bounding box
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );

    // oh boy, here come the channels!
    m_namedChannels.clear();

    // I think it is reasonable to expect that the two level sets we are interpolating will have the same channels
    // if they aren't the same size, throw an exception (there is a further check in the next loop)
    if( rleFirst.m_namedChannels.size() != rleSecond.m_namedChannels.size() ) {
        // find the name of a channel that only exists in one of the level sets
        vector<frantic::tstring> firstNames, secondNames;
        rleFirst.get_channel_names( firstNames );
        rleSecond.get_channel_names( secondNames );

        // look for the name of the channel that isn't in both level sets
        bool found = false;
        if( rleFirst.m_namedChannels.size() > rleSecond.m_namedChannels.size() ) {
            for( unsigned i = 0; i < rleFirst.m_namedChannels.size(); ++i ) {
                found = false;
                for( unsigned j = 0; j < rleSecond.m_namedChannels.size() && !found; ++j ) {
                    if( strcmp( firstNames[i].c_str(), secondNames[j].c_str() ) == 0 ) {
                        found = true;
                    }
                }
                if( !found ) {
                    throw std::runtime_error(
                        "rle_voxel_field.linear_interpolate: The two input level sets should contain the "
                        "same channels, but the second input level set does not contain a channel named: " +
                        frantic::strings::to_string( firstNames[i] ) + "." );
                }
            }
        } else {
            for( unsigned i = 0; i < rleSecond.m_namedChannels.size(); ++i ) {
                found = false;
                for( unsigned j = 0; j < rleFirst.m_namedChannels.size() && !found; ++j ) {
                    if( strcmp( secondNames[i].c_str(), firstNames[j].c_str() ) == 0 ) {
                        found = true;
                    }
                }
                if( !found ) {
                    throw std::runtime_error(
                        "rle_voxel_field.linear_interpolate: The two input level sets should contain the "
                        "same channels, but the first input level set does not contain a channel named: " +
                        frantic::strings::to_string( secondNames[i] ) + "." );
                }
            }
        }
    }

    // copy over all of the channels from rleFirst
    for( map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i =
             rleFirst.m_namedChannels.begin();
         i != rleFirst.m_namedChannels.end(); ++i ) {
        map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator j =
            rleSecond.m_namedChannels.find( i->second.name() );
        if( j == rleSecond.m_namedChannels.end() ) {
            throw std::runtime_error(
                "rle_voxel_field.linear_interpolate: The two input level sets should contain the same "
                "channels, but only the first one contains the channel " +
                frantic::strings::to_string( i->second.name() ) + "." );
        }
        if( i->second.arity() != j->second.arity() || i->second.data_type() != j->second.data_type() ) {
            throw runtime_error(
                "rle_voxel_field.linear_interpolate: The two input level sets have an incompatible channel named \"" +
                frantic::strings::to_string( i->second.name() ) + "\".  The type of the first operand is " +
                frantic::strings::to_string(
                    channels::channel_data_type_str( j->second.arity(), j->second.data_type() ) ) +
                ", while the type of the second operand is " +
                frantic::strings::to_string(
                    channels::channel_data_type_str( i->second.arity(), i->second.data_type() ) ) +
                "." );
        }
        add_channel( i->second.name(), i->second.arity(), i->second.data_type() );
    }

    // set up the channel accessors
    vector<rle_channel_general_accessor> outputChannels;
    vector<const_rle_channel_general_accessor> inputChannelsFirst;
    vector<const_rle_channel_general_accessor> inputChannelsSecond;
    for( map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i = m_namedChannels.begin();
         i != m_namedChannels.end(); ++i ) {
        outputChannels.push_back( get_channel_general_accessor( i->second.name() ) );
        inputChannelsFirst.push_back( rleFirst.get_channel_general_accessor( i->second.name() ) );
        inputChannelsSecond.push_back( rleSecond.get_channel_general_accessor( i->second.name() ) );
    }

    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            int y = ris.m_abcCoordOrigin.y + b, z = ris.m_abcCoordOrigin.z + c;

            int bFirst = y - rleFirst.m_rleIndex.m_abcCoordOrigin.y;
            int cFirst = z - rleFirst.m_rleIndex.m_abcCoordOrigin.z;
            int bSecond = y - rleSecond.m_rleIndex.m_abcCoordOrigin.y;
            int cSecond = z - rleSecond.m_rleIndex.m_abcCoordOrigin.z;

            int bcIndex = b + c * m_rleIndex.m_abcCoordSize.ysize();
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            int bcIndexFirst = -1;
            int bcIndexSecond = -1;
            // Get the range of runs for each
            int runRangeStartFirst = 0, runRangeEndFirst = 0;
            int runRangeStartSecond = 0, runRangeEndSecond = 0;

            if( (unsigned)bFirst < (unsigned)rleFirst.m_rleIndex.m_abcCoordSize.ysize() &&
                (unsigned)cFirst < (unsigned)rleFirst.m_rleIndex.m_abcCoordSize.zsize() ) {
                bcIndexFirst = bFirst + cFirst * rleFirst.m_rleIndex.m_abcCoordSize.ysize();
                runRangeStartFirst = rleFirst.m_rleIndex.m_bcToRunIndex[bcIndexFirst];
                runRangeEndFirst = rleFirst.m_rleIndex.m_bcToRunIndex[bcIndexFirst + 1] - 1;
                // If the scanline consists of a single zero-sized run, we set bcIndexFirst back
                // to -1, so that it will trigger one of the simpler cases below.
                if( runRangeStartFirst + 1 == runRangeEndFirst &&
                    rleFirst.m_rleIndex.m_runData[runRangeStartFirst].x ==
                        rleFirst.m_rleIndex.m_runData[runRangeStartFirst + 1].x ) {
                    bcIndexFirst = -1;
                }
            }

            if( (unsigned)bSecond < (unsigned)rleSecond.m_rleIndex.m_abcCoordSize.ysize() &&
                (unsigned)cSecond < (unsigned)rleSecond.m_rleIndex.m_abcCoordSize.zsize() ) {
                bcIndexSecond = bSecond + cSecond * rleSecond.m_rleIndex.m_abcCoordSize.ysize();
                runRangeStartSecond = rleSecond.m_rleIndex.m_bcToRunIndex[bcIndexSecond];
                runRangeEndSecond = rleSecond.m_rleIndex.m_bcToRunIndex[bcIndexSecond + 1] - 1;
                // If the scanline consists of a single zero-sized run, we set bcIndexSecond back
                // to -1, so that it will trigger one of the simpler cases below.
                if( runRangeStartSecond + 1 == runRangeEndSecond &&
                    rleSecond.m_rleIndex.m_runData[runRangeStartSecond].x ==
                        rleSecond.m_rleIndex.m_runData[runRangeStartSecond + 1].x ) {
                    bcIndexSecond = -1;
                }
            }

            if( bcIndexFirst == -1 && bcIndexSecond == -1 ) {
                // Create a zero-sized run, necessary for the coordinate queries to work properly
                ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
                ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
            } else if( bcIndexFirst == -1 ) {
                // Since the run for the first operand is zero-sized, we copy the scanline from the second operand and
                // interpolate with the outside/inside distance value from the first.

                for( int run = runRangeStartSecond; run != runRangeEndSecond; ++run ) {
                    // Copy the whole run index data value
                    ris.m_runData.push_back( rleSecond.m_rleIndex.m_runData[run] );
                    int dataIndexSecond = ris.m_runData.back().dataIndex;
                    if( dataIndexSecond >= 0 ) {
                        ris.m_runData.back().dataIndex = (int)totalDataSize;
                        int runSize = rleSecond.m_rleIndex.m_runData[run + 1].x - rleSecond.m_rleIndex.m_runData[run].x;
                        // Copy the data for all the other channels
                        for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            const_rle_channel_general_accessor& inputChannel = inputChannelsSecond[channelIndex];

                            if( inputChannel.data_type() != frantic::channels::data_type_int8 &&
                                inputChannel.data_type() != frantic::channels::data_type_int16 &&
                                inputChannel.data_type() != frantic::channels::data_type_int32 &&
                                inputChannel.data_type() != frantic::channels::data_type_int64 ) {
                                // Channels are only coming from the second input values get interpolated with 0 ( and
                                // (1-alpha)*0 is 0, so that term is not explicitly written out)

                                channel_weighted_sum_combine_function_t ws = channel.m_weightedSumFunction;
                                size_t arity = channel.arity();
                                const char* data[1];
                                float blendingAlpha[1];
                                blendingAlpha[0] = alpha;

                                // This writes out scaled values of the input data, using the input alpha values
                                for( int i = 0; i < runSize; ++i ) {
                                    data[0] = inputChannel.data( dataIndexSecond );
                                    ws( blendingAlpha, data, 1, arity, channel.add_element() );
                                }

                            } else {
                                // integer channels are just copied over from the second input

                                size_t dataSize = runSize * channel.primitive_size();
                                memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( dataIndexSecond ),
                                        dataSize );
                            }
                        }
                        totalDataSize += runSize;
                    }
                }
                // Copy the run index data for the value that finishes off the last run of the scanline
                ris.m_runData.push_back( rleSecond.m_rleIndex.m_runData[runRangeEndSecond] );
            } else if( bcIndexSecond == -1 ) {
                // Since the run for the second operand is zero-sized, we copy the scanline from the first operand and
                // interpolate with the outside/inside distance value from the second.

                if( rleSecond.m_rleIndex.m_exteriorRegionCode != -1 ) {
                    // Create a zero-sized run
                    ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
                    ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
                } else {

                    for( int run = runRangeStartFirst; run != runRangeEndFirst; ++run ) {
                        // Copy the whole run index data value
                        ris.m_runData.push_back( rleFirst.m_rleIndex.m_runData[run] );
                        int dataIndexFirst = ris.m_runData.back().dataIndex;
                        if( dataIndexFirst >= 0 ) {
                            ris.m_runData.back().dataIndex = (int)totalDataSize;
                            int runSize =
                                rleFirst.m_rleIndex.m_runData[run + 1].x - rleFirst.m_rleIndex.m_runData[run].x;
                            // Copy the data for all the other channels
                            for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                                rle_channel_general_accessor& channel = outputChannels[channelIndex];
                                const_rle_channel_general_accessor& inputChannel = inputChannelsFirst[channelIndex];
                                if( inputChannel.data_type() != frantic::channels::data_type_int8 &&
                                    inputChannel.data_type() != frantic::channels::data_type_int16 &&
                                    inputChannel.data_type() != frantic::channels::data_type_int32 &&
                                    inputChannel.data_type() != frantic::channels::data_type_int64 ) {
                                    // Channels only coming from the first input get values interpolated with 0 ( and
                                    // alpha*0 is 0, so that term is not explicitly written out)

                                    channel_weighted_sum_combine_function_t ws = channel.m_weightedSumFunction;
                                    size_t arity = channel.arity();
                                    const char* data[1];
                                    float blendingAlpha[1];
                                    blendingAlpha[0] = 1 - alpha;

                                    // This writes out scaled values of the input data, using the input alpha values
                                    for( int i = 0; i < runSize; ++i ) {
                                        data[0] = inputChannel.data( dataIndexFirst );
                                        ws( blendingAlpha, data, 1, arity, channel.add_element() );
                                    }

                                } else {
                                    // integer channels are just copied over from the first input
                                    size_t dataSize = runSize * channel.primitive_size();
                                    memcpy( channel.m_data->add_element( dataSize ),
                                            inputChannel.data( dataIndexFirst ), dataSize );
                                }
                            }
                            totalDataSize += runSize;
                        }
                    }
                    // Copy the run index data for the value that finishes off the last run of the scanline
                    ris.m_runData.push_back( rleFirst.m_rleIndex.m_runData[runRangeEndFirst] );
                }
            } else {
                // In this case, both scanlines have a set of valid runs.  The juggling about
                // we did above guarantees that when either run is zero-sized, the previous cases
                // will catch it.

                // There's an ordering of the types of regions we encounter during processing.
                // Constructing it so that defined values always win

                // The default state is an undefined region of the exterior region code
                bool creatingUndefinedRun = true;
                int creatingRegionCode = ris.m_exteriorRegionCode;
                int creatingCurrentRunX = 0;

                // These are the two run index values, which we increment side by side in a merge-sort like fashion.
                int runFirst = runRangeStartFirst;
                int runSecond = runRangeStartSecond;

                // Initially, both the inputs we're processing begin with an undefined virtual run
                // extending from negative infinity to one before the start of the first run, that has
                // the exterior region code.
                bool processingUndefinedRunFirst = true;
                int processingRegionCodeFirst = rleFirst.m_rleIndex.m_exteriorRegionCode;
                int processingNextRunXFirst = rleFirst.m_rleIndex.m_runData[runFirst].x;
                bool processingUndefinedRunSecond = true;
                int processingRegionCodeSecond = rleSecond.m_rleIndex.m_exteriorRegionCode;
                int processingNextRunXSecond = rleSecond.m_rleIndex.m_runData[runSecond].x;

                ////////////
                // Loop through all the segments created by taking the union of all the run starts in both the first and
                // second operand.
                ////////////

                for( ;; ) {
                    ////////////
                    // Set up the segment we'll be processing in the next iteration of the loop
                    ////////////

                    // To set up the next iteration of the loop, we need to advance either the run of the first scanline
                    // or the run of the second scanline depending on which is closer.
                    // If both of them are equal, we increment both.
                    bool incrementFirst = processingNextRunXFirst <= processingNextRunXSecond;
                    bool incrementSecond = processingNextRunXFirst >= processingNextRunXSecond;
                    if( incrementFirst ) {
                        creatingCurrentRunX = processingNextRunXFirst;

                        if( runFirst != runRangeEndFirst ) {
                            processingRegionCodeFirst = rleFirst.m_rleIndex.m_runData[runFirst].dataIndex;
                            processingUndefinedRunFirst = processingRegionCodeFirst < 0;
                            processingNextRunXFirst = rleFirst.m_rleIndex.m_runData[++runFirst].x;
                        } else {
                            processingRegionCodeFirst = rleFirst.m_rleIndex.m_exteriorRegionCode;
                            processingUndefinedRunFirst = true;
                            processingNextRunXFirst = ( std::numeric_limits<int>::max )();
                        }
                    }
                    if( incrementSecond ) {
                        creatingCurrentRunX = processingNextRunXSecond;

                        if( runSecond != runRangeEndSecond ) {
                            processingRegionCodeSecond = rleSecond.m_rleIndex.m_runData[runSecond].dataIndex;
                            processingUndefinedRunSecond = processingRegionCodeSecond < 0;
                            processingNextRunXSecond = rleSecond.m_rleIndex.m_runData[++runSecond].x;
                        } else {
                            processingRegionCodeSecond = rleSecond.m_rleIndex.m_exteriorRegionCode;
                            processingUndefinedRunSecond = true;
                            processingNextRunXSecond = ( std::numeric_limits<int>::max )();
                        }
                    }

                    // If we stepped past the end, stop the looping
                    if( processingNextRunXFirst == ( std::numeric_limits<int>::max )() &&
                        processingNextRunXSecond == ( std::numeric_limits<int>::max )() ) {
                        break;
                    }

                    ////////////
                    // Determine the size of the segment we're processing in this loop iteration
                    ////////////

                    int segmentSize =
                        ( std::min )( processingNextRunXFirst, processingNextRunXSecond ) - creatingCurrentRunX;

                    ////////////
                    // Process the data within this segment, either doing a jump by creating an undefined run, or by
                    // looping through the data of the segment
                    ////////////

                    if( processingUndefinedRunFirst && processingUndefinedRunSecond ) {

                        // In this case, both runs are undefined, so we can make a leap collapsing this run into an
                        // undefined run

                        // By default, use the exterior region code
                        int undefinedRegionCodeToUse = ris.m_exteriorRegionCode;

                        // If the two undefined runs have the same region code (most cases), use this shared region code
                        // instead of the default undefined outside.
                        if( processingRegionCodeFirst == processingRegionCodeSecond )
                            undefinedRegionCodeToUse = processingRegionCodeFirst;

                        if( creatingUndefinedRun ) {
                            // If the undefined region doesn't match the current region, start a new run
                            if( creatingRegionCode != undefinedRegionCodeToUse ) {
                                ris.m_runData.push_back( run_data( creatingCurrentRunX, undefinedRegionCodeToUse ) );
                                creatingRegionCode = undefinedRegionCodeToUse;
                            }
                        } else {
                            // If we were making a defined run, start making an undefined run
                            ris.m_runData.push_back( run_data( creatingCurrentRunX, undefinedRegionCodeToUse ) );
                            creatingUndefinedRun = true;
                            creatingRegionCode = undefinedRegionCodeToUse;
                        }

                    } else if( processingUndefinedRunFirst ) {
                        // In this case, the second run is defined, but the first is undefined "outside", so we can copy
                        // the data from the second run (and interpolate it with the outside distance from the first
                        // run)
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back( run_data( creatingCurrentRunX, (int)totalDataSize ) );
                            creatingUndefinedRun = false;
                        }
                        // get the proper inside or outside distance for interpolation
                        int dataIndexSecond = rleSecond.m_rleIndex.m_runData[runSecond - 1].dataIndex +
                                              ( creatingCurrentRunX - rleSecond.m_rleIndex.m_runData[runSecond - 1].x );
                        // Copy the data for all the other channels
                        for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            const_rle_channel_general_accessor& inputChannel = inputChannelsSecond[channelIndex];
                            if( inputChannel.data_type() != frantic::channels::data_type_int8 &&
                                inputChannel.data_type() != frantic::channels::data_type_int16 &&
                                inputChannel.data_type() != frantic::channels::data_type_int32 &&
                                inputChannel.data_type() != frantic::channels::data_type_int64 ) {
                                // Channels only coming from the second input get values interpolated with 0 ( and
                                // (1-alpha)*0 is 0, so that term is not explicitly written out)

                                channel_weighted_sum_combine_function_t ws = channel.m_weightedSumFunction;
                                size_t arity = channel.arity();
                                const char* data[1];
                                float blendingAlpha[1];
                                blendingAlpha[0] = alpha;

                                // This writes out scaled values of the input data, using the input alpha values
                                for( int i = 0; i < segmentSize; ++i ) {
                                    data[0] = inputChannel.data( dataIndexSecond );
                                    ws( blendingAlpha, data, 1, arity, channel.add_element() );
                                }

                                // size_t dataSize = segmentSize * channel.primitive_size();
                                // memcpy( channel.m_data->add_element(dataSize), inputChannel.data(dataIndexSecond),
                                // dataSize );
                            } else {
                                // integer channels are just copied over from the second input

                                size_t dataSize = segmentSize * channel.primitive_size();
                                memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( dataIndexSecond ),
                                        dataSize );
                            }
                        }
                        totalDataSize += segmentSize;
                    } else if( processingUndefinedRunSecond ) {
                        // In this case, the first run is defined, but the second is undefined "outside", so we can copy
                        // the data from the first run (and interpolate it with the outside distance from the second
                        // run)
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back( run_data( creatingCurrentRunX, (int)totalDataSize ) );
                            creatingUndefinedRun = false;
                        }
                        int dataIndexFirst = rleFirst.m_rleIndex.m_runData[runFirst - 1].dataIndex +
                                             ( creatingCurrentRunX - rleFirst.m_rleIndex.m_runData[runFirst - 1].x );
                        // Copy the data for all the other channels
                        for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            const_rle_channel_general_accessor& inputChannel = inputChannelsFirst[channelIndex];
                            if( inputChannel.data_type() != frantic::channels::data_type_int8 &&
                                inputChannel.data_type() != frantic::channels::data_type_int16 &&
                                inputChannel.data_type() != frantic::channels::data_type_int32 &&
                                inputChannel.data_type() != frantic::channels::data_type_int64 ) {
                                // Channels only coming from the first input get values interpolated with 0 ( and
                                // alpha*0 is 0, so that term is not explicitly written out)

                                channel_weighted_sum_combine_function_t ws = channel.m_weightedSumFunction;
                                size_t arity = channel.arity();
                                const char* data[1];
                                float blendingAlpha[1];
                                blendingAlpha[0] = 1 - alpha;

                                // This writes out scaled values of the input data, using the input alpha values
                                for( int i = 0; i < segmentSize; ++i ) {
                                    data[0] = inputChannel.data( dataIndexFirst );
                                    ws( blendingAlpha, data, 1, arity, channel.add_element() );
                                }

                                // size_t dataSize = segmentSize * channel.primitive_size();
                                // memcpy( channel.m_data->add_element(dataSize), inputChannel.data(dataIndexFirst),
                                // dataSize );
                            } else {
                                // integer channels are just copied over from the first input
                                size_t dataSize = segmentSize * channel.primitive_size();
                                memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( dataIndexFirst ),
                                        dataSize );
                            }
                        }
                        totalDataSize += segmentSize;
                    } else {
                        // In this case, both runs are defined, so we have to interpolate the data voxel-by-voxel
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back( run_data( creatingCurrentRunX, (int)totalDataSize ) );
                            creatingUndefinedRun = false;
                        }
                        int dataIndexFirstBegin =
                            rleFirst.m_rleIndex.m_runData[runFirst - 1].dataIndex +
                            ( creatingCurrentRunX - rleFirst.m_rleIndex.m_runData[runFirst - 1].x );
                        int dataIndexSecondBegin =
                            rleSecond.m_rleIndex.m_runData[runSecond - 1].dataIndex +
                            ( creatingCurrentRunX - rleSecond.m_rleIndex.m_runData[runSecond - 1].x );
                        int dataIndexFirst = dataIndexFirstBegin, dataIndexSecond = dataIndexSecondBegin;

                        totalDataSize += segmentSize;

                        // interpolate the channel values
                        for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                            dataIndexFirst = dataIndexFirstBegin;
                            dataIndexSecond = dataIndexSecondBegin;

                            if( inputChannelsFirst[channelIndex].data_type() != frantic::channels::data_type_int8 &&
                                inputChannelsFirst[channelIndex].data_type() != frantic::channels::data_type_int16 &&
                                inputChannelsFirst[channelIndex].data_type() != frantic::channels::data_type_int32 &&
                                inputChannelsFirst[channelIndex].data_type() != frantic::channels::data_type_int64 ) {

                                rle_channel_general_accessor& channel = outputChannels[channelIndex];
                                const_rle_channel_general_accessor& inputChannelFirst =
                                    inputChannelsFirst[channelIndex];
                                const_rle_channel_general_accessor& inputChannelSecond =
                                    inputChannelsSecond[channelIndex];

                                // Channels only coming from both inputs get interpolated

                                channel_weighted_sum_combine_function_t ws = channel.m_weightedSumFunction;
                                size_t arity = channel.arity();
                                const char* data[2];

                                vector<float> blendingAlpha;
                                blendingAlpha.resize( 2 );
                                blendingAlpha[0] = 1 - alpha;
                                blendingAlpha[1] = alpha;

                                // This writes out scaled values of the input data, using the input alpha values
                                for( int i = 0; i < segmentSize; ++i ) {
                                    data[0] = inputChannelFirst.data( dataIndexFirst++ );
                                    data[1] = inputChannelSecond.data( dataIndexSecond++ );
                                    ws( &blendingAlpha[0], data, 2, arity, channel.add_element() );
                                }

                                // size_t dataSize = segmentSize * channel.primitive_size();
                                // memcpy( channel.m_data->add_element(dataSize), inputChannelFirst.data(dataIndexFirst)
                                // , dataSize );

                            } else {
                                // integer channels are just copied over from the closer of the two level sets
                                rle_channel_general_accessor& channel = outputChannels[channelIndex];
                                int channelDataIndex = dataIndexFirst;
                                const_rle_channel_general_accessor& inputChannel = inputChannelsFirst[channelIndex];
                                if( alpha > 0.5 ) {
                                    inputChannel = inputChannelsSecond[channelIndex];
                                    channelDataIndex = dataIndexSecond;
                                }
                                size_t dataSize = segmentSize * channel.primitive_size();
                                memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( channelDataIndex ),
                                        dataSize );
                            }
                        }
                    }
                }
                // If any runs were flagged, finish off the last run.  Otherwise add a zero-sized run
                if( (int)ris.m_runData.size() > ris.m_bcToRunIndex[bcIndex] ) {
                    // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to
                    // optimize the storage and access.
                    if( !( creatingUndefinedRun && creatingRegionCode == ris.m_exteriorRegionCode ) )
                        ris.m_runData.push_back( run_data( creatingCurrentRunX, -1 ) );
                } else {
                    ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
                    ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
                }
            }
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();

    // Set the data size as well
    ris.m_dataSize = totalDataSize;

    // if the result is of size (0,0,0) return a single outside run of size (1,1,1)
    if( ris.m_abcCoordSize.xsize() == 0 && ris.m_abcCoordSize.ysize() == 0 && ris.m_abcCoordSize.zsize() == 0 ) {
        set_to_empty( ris.m_exteriorRegionCode );
    }
}

void rle_voxel_field::convert_staggered_velocity_to_centered(
    frantic::volumetrics::levelset::rle_level_set& destLevelSet, const frantic::tstring& staggeredChannelName,
    const frantic::tstring& centeredChannelName ) const {
    const_rle_channel_accessor<vector3f> staggeredAccessor = get_channel_accessor<vector3f>( staggeredChannelName );

    destLevelSet.add_channel<vector3f>( centeredChannelName );
    rle_channel_accessor<vector3f> centeredAccessor =
        destLevelSet.get_channel_accessor<vector3f>( centeredChannelName );

    const rle_index_spec& destRleIndexSpec = destLevelSet.get_rle_index_spec();
    //	const voxel_coord_system& vcs = destLevelSet.get_voxel_coord_system();

    rle_defined_iterator i = destRleIndexSpec.begin();
    rle_defined_iterator iend = destRleIndexSpec.end();

    for( ; i != iend; ++i ) {
        vector3 voxelCoord = i.get_coord();
        vector3f voxelCenter = voxelCoord + vector3f( 0.5f );

        if( m_rleIndex.XYZtoDataIndex( voxelCoord ) >= 0 ) {
            trilerp_vector3f_staggered( m_rleIndex, staggeredAccessor, voxelCenter,
                                        centeredAccessor[i.get_data_index()] );
        } else
            centeredAccessor[i.get_data_index()] = vector3f();
    }
}

void rle_voxel_field::convert_staggered_velocity_to_centered( const frantic::tstring& staggeredChannelName,
                                                              const frantic::tstring& centeredChannelName ) {
    const_rle_channel_accessor<vector3f> staggeredAccessor = get_channel_accessor<vector3f>( staggeredChannelName );

    add_channel<vector3f>( centeredChannelName );
    rle_channel_accessor<vector3f> centeredAccessor = get_channel_accessor<vector3f>( centeredChannelName );

    rle_defined_iterator i = m_rleIndex.begin();
    rle_defined_iterator iend = m_rleIndex.end();

    for( ; i != iend; ++i ) {
        vector3 voxelCoord = i.get_coord();
        vector3f voxelCenter = voxelCoord + vector3f( 0.5f );
        trilerp_vector3f_staggered( m_rleIndex, staggeredAccessor, voxelCenter, centeredAccessor[i.get_data_index()] );
    }
}

void rle_voxel_field::trim_to_bounds( const frantic::graphics::boundbox3& trimVoxelBounds ) {
    // Only do the trimming if the outer bounds are bigger than the trim bounds
    if( !m_rleIndex.outer_bounds().is_empty() && !trimVoxelBounds.contains( m_rleIndex.outer_bounds() ) ) {
        // To do the trimming, create a temp trimmed rle index spec, then switch the level set using
        // a swap on the rle index spec for efficiency.
        rle_index_spec risTemp;
        risTemp.build_with_trim_bounds( m_rleIndex, trimVoxelBounds );
        switch_rle_index_spec_with_swap( risTemp );
    }
}

void rle_voxel_field::switch_rle_index_spec_with_swap( rle_index_spec& ris,
                                                       const frantic::tstring& populatedChannelToCreate ) {
    // Initialize the level set we're computing from this
    rle_voxel_field result( m_voxelCoordSystem );

    // Swap in the provided rle_index_spec
    result.m_rleIndex.swap( ris );

    // Create all the data channels we need in the result
    for( std::map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.begin();
         i != m_namedChannels.end(); ++i ) {
        result.add_channel( i->second.name(), i->second.arity(), i->second.data_type() );
        // Initialize each new channel to all 0 bytes
        std::map<frantic::tstring, rle_channel>::iterator iResult = result.m_namedChannels.find( i->second.name() );
        memset( iResult->second.m_data.begin(), 0, iResult->second.m_data.size() );
    }

    // Create the populated channel.  By initializing all the values in the input to 1, and all the values
    // in the output to 0, we will automatically end up with a result that has 1 wherever there was data defined
    // and 0 whereever a new defined voxel was created.
    if( !populatedChannelToCreate.empty() ) {
        // Add the populated channel to the source level set, setting all its values to 1
        add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( populatedChannelToCreate );
        memset( i->second.m_data.begin(), 1, i->second.m_data.size() );

        // Add the populated channel to the result level set, setting all its values to 0
        result.add_channel( populatedChannelToCreate, 1, channels::data_type_uint8 );
        i = result.m_namedChannels.find( populatedChannelToCreate );
        memset( i->second.m_data.begin(), 0, i->second.m_data.size() );
    }

    // Build the accessors to all the named channels
    vector<const_rle_channel_general_accessor> inputAccessors;
    vector<rle_channel_general_accessor> outputAccessors;
    for( std::map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.begin();
         i != m_namedChannels.end(); ++i ) {
        inputAccessors.push_back( get_channel_general_accessor( i->first ) );
        outputAccessors.push_back( result.get_channel_general_accessor( i->first ) );
    }
    // Fill the channel copying data array
    vector<rle_index_spec_channel_copying_data> channelCopyData( inputAccessors.size() );
    for( size_t i = 0, ie = inputAccessors.size(); i != ie; ++i ) {
        rle_index_spec_channel_copying_data& cd = channelCopyData[i];
        cd.inputData = inputAccessors[i].data( 0 );
        cd.outputData = outputAccessors[i].data( 0 );
        cd.primitiveSize = inputAccessors[i].primitive_size();
    }
    // Use the rle_index_spec method to copy all the defined data.
    rle_index_spec::copy_data_channels( result.m_rleIndex, m_rleIndex, &channelCopyData[0],
                                        &channelCopyData[0] + channelCopyData.size() );

    // Swap the answer back into this
    result.swap( *this );
}

void rle_voxel_field::apply_axis_permutation( const vector3& axisPermutation ) {
    // Skip any processing if the permutation is the identity
    if( axisPermutation.x != 0 || axisPermutation.y != 1 || axisPermutation.z != 2 ) {
        // NOTE: this is not a vector, because vectors could copy the data and corrupt it
        deque<raw_byte_buffer> tempProcessingChannelData;
        vector<rle_index_spec_channel_copying_data> channelsToRemap;
        // Go through all the channels, and set up the channels to remap during the axis permutation
        for( map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.begin(); i != m_namedChannels.end();
             ++i ) {
            // Allocate a new buffer for this channel, and swap it so the destination is within the named channel data
            // structure
            tempProcessingChannelData.push_back( raw_byte_buffer() );
            tempProcessingChannelData.back().resize( i->second.m_data.size() );
            tempProcessingChannelData.back().swap( i->second.m_data );

            // Set up the entry for remapping this channel
            channelsToRemap.push_back( rle_index_spec_channel_copying_data() );
            channelsToRemap.back().primitiveSize = i->second.primitive_size();
            channelsToRemap.back().inputData = tempProcessingChannelData.back().begin();
            channelsToRemap.back().outputData = i->second.m_data.begin();
        }

        // Now apply the axis permutation, providing the data to remap
        if( channelsToRemap.empty() )
            m_rleIndex.apply_axis_permutation( axisPermutation );
        else
            m_rleIndex.apply_axis_permutation( axisPermutation, &channelsToRemap[0],
                                               &channelsToRemap[0] + channelsToRemap.size() );
    }

    // run a consistency check (debugging code, leave it in for now)
    stringstream sout;
    if( !m_rleIndex.check_consistency( sout ) ) {
        throw std::runtime_error(
            "rle_level_set.apply_axis_permutation: RLE Index Spec consistency check failed (permutation=" +
            axisPermutation.str() + "):\n" + sout.str() );
    }
}

void rle_voxel_field::convert_centered_velocity_to_staggered( const rle_level_set& sourceLevelSet,
                                                              const frantic::tstring& centeredChannel,
                                                              const frantic::tstring& staggeredChannel,
                                                              const frantic::tstring& dataIndexMappingChannel ) {
    // logging::debug << "Converting velocities from " << centeredChannel << " to staggered Channel " <<
    // staggeredChannel<< endl;

    if( !has_channel( staggeredChannel ) ) {
        throw std::runtime_error( "Voxel Field does not have the requested staggered velocity channel: \"" +
                                  frantic::strings::to_string( staggeredChannel ) + "\"" );
    }

    if( !m_voxelCoordSystem.equals( sourceLevelSet.get_voxel_coord_system() ) ) {
        throw std::runtime_error(
            "convert_centered_velocity_to_staggered() - The voxel coordinate systems of the voxel "
            "field and the level set must match\n" +
            m_voxelCoordSystem.str() + " vs. " + sourceLevelSet.get_voxel_coord_system().str() );
    }

    const rle_index_spec& sourceRIS = sourceLevelSet.get_rle_index_spec();
    rle_channel_accessor<vector3f> staggeredAccessor = get_channel_accessor<vector3f>( staggeredChannel );
    const_rle_channel_accessor<vector3f> centeredAccessor =
        sourceLevelSet.get_channel_accessor<vector3f>( centeredChannel );
    const_rle_channel_accessor<boost::int32_t> realDataIndexMapChannelAccessor;

    std::vector<boost::int32_t> dataIndexMapChannelArray;
    const boost::int32_t* dataIndexMapChannelAccessor;

    // if required build a temporary buffer to store the data index mapping data, otherwise use the one provided
    if( dataIndexMappingChannel.empty() ) {
        dataIndexMapChannelArray.resize( sourceLevelSet.size() );
        sourceRIS.fill_data_index_map( m_rleIndex, &dataIndexMapChannelArray[0] );
        dataIndexMapChannelAccessor = &dataIndexMapChannelArray[0];
    } else {
        realDataIndexMapChannelAccessor =
            sourceLevelSet.get_channel_accessor<boost::int32_t>( dataIndexMappingChannel );
        dataIndexMapChannelAccessor = &realDataIndexMapChannelAccessor[0];
    }

    // build the defined and adjacent iterators so that we can march over the level set
    rle_defined_and_adj_iterator i( sourceRIS );
    rle_defined_and_adj_iterator iend( sourceRIS, true );

    for( ; i != iend; ++i ) {
        boost::int32_t dataIndex = i.get_center_data_index();
        boost::int32_t staggeredIndex = dataIndexMapChannelAccessor[dataIndex];

        // we can perform an interpolation of the velocities and set the face
        if( staggeredIndex >= 0 ) {

            for( int neighbour = 0; neighbour < 6; ++neighbour ) {
                boost::int32_t adj = i.get_adjacent_data_index( neighbour );

                // -ve faces
                if( !is_neighbor_index_direction_positive( neighbour ) ) {
                    if( adj >= 0 ) {
                        // our neighbour adds a contribution to this face too
                        staggeredAccessor[staggeredIndex][neighbour / 2] =
                            0.5f *
                            ( centeredAccessor[dataIndex][neighbour / 2] + centeredAccessor[adj][neighbour / 2] );
                    } else {
                        staggeredAccessor[staggeredIndex][neighbour / 2] = centeredAccessor[dataIndex][neighbour / 2];
                    }
                }
                // +ve faces
                // these will be handled by the -ve face section when their voxel is reached
            }
        }
    }
}

} // namespace fluids
} // namespace frantic
