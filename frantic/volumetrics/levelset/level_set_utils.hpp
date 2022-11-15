// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/mpl/assert.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/none.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/levelset/rle_pairwise_run_iterator.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {
/**
 * functor for sorting level sets based on their size
 * used for sorting level sets for optimal ordering of csg boolean operations
 */
template <class LevelSet>
struct rle_size_compare {
    bool operator()( const std::pair<size_t, LevelSet>& rhs, const std::pair<size_t, LevelSet>& lhs ) {
        if( rhs.first < lhs.first ) {
            return true;
        }

        return false;
    }
};

/**
 * For a collection of levelsets, this function unions all of them together where there are any bounding box
 * intersections returns of a set of level sets that have mutually exclusive inside regions. This function allocates new
 * level sets where necessary to hold the union results.
 *
 * @param levelSets The collection of level sets to union
 * @param[out] unionResults a collection of level sets whose inside regions are exclusive (i.e have no intersections
 * with each other)
 */
void multiple_csg_union( std::vector<boost::shared_ptr<rle_level_set>>& levelSets,
                         std::vector<boost::shared_ptr<rle_level_set>>& unionResults );

/*
 * Semi lagrangian advection a level set by some frame offset using its own velocity field. Performance warning: This
 * function needs to make a copy of the level set that is later discarded.
 *
 * @param levelSet  The level set to advect
 * @param frameOffset  The amount to advect the level set
 */
void self_advect_rle_level_set( frantic::volumetrics::levelset::rle_level_set& levelSet, float frameOffset );

/*
 * Semi lagrangian advection a level set by some frame offset using its own velocity field. Performance warning: This
 * function needs to make a copy of the level set that is later discarded. Normally you would call the above version
 * unless you want to provide your own copy of the levelSet.
 *
 * @param levelSet  The level set to advect
 * @param frameOffset  The amount to advect the level set
 */
void self_advect_rle_level_set( frantic::volumetrics::levelset::rle_level_set& levelSet,
                                const frantic::volumetrics::levelset::rle_level_set& tempLevelSetCopy,
                                float frameOffset );

/*
 * Finds the maximum velocity in a level set, then dilates the level set to that distance * frameOffset.
 * This function was designed fore, and should usually be called before, "self_advect_rle_level_set".
 *
 * @param levelSet  level set to dilate
 * @param timeStep  if portion of the frame it needs to be dilated to include
 */
void expand_rle_level_set_to_max_velocity( frantic::volumetrics::levelset::rle_level_set& levelSet, float timeStep );

namespace detail {
// treat SignedDistance channel without special case in copy_channel
// this may not be worthwhile, and it could certainly be done better
bool has_channel( const frantic::fluids::rle_voxel_field& rvf, const frantic::tstring& channelName );

template <class T>
void add_channel( frantic::fluids::rle_voxel_field& rvf, const frantic::tstring& channelName ) {
    rvf.add_channel<T>( channelName );
}

void add_channel( frantic::fluids::rle_voxel_field& rvf, const frantic::tstring& channelName, std::size_t arity,
                  data_type_t dataType );

std::pair<std::size_t, data_type_t> get_channel_type( const frantic::fluids::rle_voxel_field& rvf,
                                                      const frantic::tstring& channelName );

template <class T>
T* get_channel_data_pointer( frantic::fluids::rle_voxel_field& rvf, const frantic::tstring& channelName ) {
    if( !rvf.has_channel( channelName ) )
        throw std::runtime_error(
            "get_channel_data_pointer Error: the rle_voxel_field does not have the specified channel \'" +
            frantic::strings::to_string( channelName ) + "\'." );
    else {
        rle_channel_general_accessor tempAcc = rvf.get_channel_general_accessor( channelName );
        if( frantic::channels::channel_data_type_traits<T>::arity() != tempAcc.arity() ||
            frantic::channels::channel_data_type_traits<T>::data_type() != tempAcc.data_type() ) {
            throw std::runtime_error( "get_channel_data_pointer Error: requested pointer for channel \'" +
                                      frantic::strings::to_string( channelName ) + "\' to type " +
                                      frantic::channels::channel_data_type_traits<T>::type_str() +
                                      " but the channel is of type " + tempAcc.type_str() + "." );
        }
    }

    rle_channel_accessor<T> acc = rvf.get_channel_accessor<T>( channelName );
    if( acc.size() == 0 )
        return NULL;
    return &acc[0];
}

template <class T>
const T* get_channel_data_pointer( const frantic::fluids::rle_voxel_field& rvf, const frantic::tstring& channelName ) {
    if( !rvf.has_channel( channelName ) )
        throw std::runtime_error(
            "get_channel_data_pointer Error: the rle_voxel_field does not have the specified channel \'" +
            frantic::strings::to_string( channelName ) + "\'." );
    else {
        const_rle_channel_general_accessor tempAcc = rvf.get_channel_general_accessor( channelName );
        if( frantic::channels::channel_data_type_traits<T>::arity() != tempAcc.arity() ||
            frantic::channels::channel_data_type_traits<T>::data_type() != tempAcc.data_type() ) {
            throw std::runtime_error( "get_channel_data_pointer Error: requested pointer for channel \'" +
                                      frantic::strings::to_string( channelName ) + "\' to type " +
                                      frantic::channels::channel_data_type_traits<T>::type_str() +
                                      " but the channel is of type " + tempAcc.type_str() + "." );
        }
    }

    const_rle_channel_accessor<T> acc = rvf.get_channel_accessor<T>( channelName );
    if( acc.size() == 0 )
        return NULL;
    return &acc[0];
}

bool has_channel( const frantic::volumetrics::levelset::rle_level_set& rls, const frantic::tstring& channelName );

template <class T>
void add_channel( frantic::volumetrics::levelset::rle_level_set& rls, const frantic::tstring& channelName ) {
    if( channelName == _T("SignedDistance") ) {
        if( frantic::channels::channel_data_type_traits<T>::arity() !=
                frantic::channels::channel_data_type_traits<float>::arity() ||
            frantic::channels::channel_data_type_traits<T>::data_type() !=
                frantic::channels::channel_data_type_traits<float>::data_type() ) {
            throw std::runtime_error(
                "add_channel Error: requested adding a \'SignedDistance\' channel of type " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) +
                " to a rle_level_set, but this channel must be of type " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<float>::type_str() ) + "." );
        }

    } else {
        rls.add_channel<T>( channelName );
    }
}

void add_channel( frantic::volumetrics::levelset::rle_level_set& rls, const frantic::tstring& channelName,
                  std::size_t arity, data_type_t dataType );

std::pair<std::size_t, data_type_t> get_channel_type( const frantic::volumetrics::levelset::rle_level_set& rls,
                                                      const frantic::tstring& channelName );

template <class T>
const T* get_signed_distance_data_pointer( const frantic::volumetrics::levelset::rle_level_set& /*rls*/ ) {
    BOOST_MPL_ASSERT( (boost::mpl::not_<boost::is_same<typename boost::remove_const<T>::type, float>>));
    throw std::runtime_error(
        "get_signed_distance_data_pointer Error: requested pointer for channel \'SignedDistance\' to type " +
        frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) +
        " but the channel is of type " +
        frantic::strings::to_string( frantic::channels::channel_data_type_traits<float>::type_str() ) + "." );
}

template <class T>
T* get_signed_distance_data_pointer( frantic::volumetrics::levelset::rle_level_set& /*rls*/ ) {
    BOOST_MPL_ASSERT( (boost::mpl::not_<boost::is_same<typename boost::remove_const<T>::type, float>>));
    throw std::runtime_error(
        "get_signed_distance_data_pointer Error: requested pointer for channel \'SignedDistance\' to type " +
        frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) +
        " but the channel is of type " +
        frantic::strings::to_string( frantic::channels::channel_data_type_traits<float>::type_str() ) + "." );
}

template <>
const float* get_signed_distance_data_pointer<float>( const frantic::volumetrics::levelset::rle_level_set& rls );

template <>
float* get_signed_distance_data_pointer<float>( frantic::volumetrics::levelset::rle_level_set& rls );

template <class T>
T* get_channel_data_pointer( frantic::volumetrics::levelset::rle_level_set& rls, const frantic::tstring& channelName ) {
    if( channelName == _T("SignedDistance") ) {
        return get_signed_distance_data_pointer<T>( rls );
    }
    if( !rls.has_channel( channelName ) )
        throw std::runtime_error(
            "get_channel_data_pointer Error: the rle_level_set does not have the specified channel \'" +
            frantic::strings::to_string( channelName ) + "\'." );
    else {
        rle_channel_general_accessor tempAcc = rls.get_channel_general_accessor( channelName );
        if( frantic::channels::channel_data_type_traits<T>::arity() != tempAcc.arity() ||
            frantic::channels::channel_data_type_traits<T>::data_type() != tempAcc.data_type() ) {
            throw std::runtime_error(
                "get_channel_data_pointer Error: requested pointer for channel \'" +
                frantic::strings::to_string( channelName ) + "\' to type " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) +
                " but the channel is of type " + frantic::strings::to_string( tempAcc.type_str() ) + "." );
        }
    }

    rle_channel_accessor<T> acc = rls.get_channel_accessor<T>( channelName );
    if( acc.size() == 0 )
        return NULL;
    return &acc[0];
}

template <class T>
const T* get_channel_data_pointer( const frantic::volumetrics::levelset::rle_level_set& rls,
                                   const frantic::tstring& channelName ) {
    if( channelName == _T("SignedDistance") ) {
        return get_signed_distance_data_pointer<T>( rls );
    }

    if( !rls.has_channel( channelName ) )
        throw std::runtime_error(
            "get_channel_data_pointer Error: the rle_level_set does not have the specified channel \'" +
            frantic::strings::to_string( channelName ) + "\'." );
    else {
        const_rle_channel_general_accessor tempAcc = rls.get_channel_general_accessor( channelName );
        if( frantic::channels::channel_data_type_traits<T>::arity() != tempAcc.arity() ||
            frantic::channels::channel_data_type_traits<T>::data_type() != tempAcc.data_type() ) {
            throw std::runtime_error(
                "get_channel_data_pointer Error: requested pointer for channel \'" +
                frantic::strings::to_string( channelName ) + "\' to type " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) +
                " but the channel is of type " + frantic::strings::to_string( tempAcc.type_str() ) + "." );
        }
    }

    const_rle_channel_accessor<T> acc = rls.get_channel_accessor<T>( channelName );
    if( acc.size() == 0 )
        return NULL;
    return &acc[0];
}

template <class T, class data_type>
void set_array_range_to_value( T* destChannel, boost::uint8_t* destPopulatedChannel, std::size_t startIndex,
                               std::size_t endIndex, data_type value ) {
    if( destChannel == 0 )
        return;
    if( destPopulatedChannel ) {
        for( std::size_t i = startIndex; i != endIndex; ++i ) {
            destChannel[i] = value;
            destPopulatedChannel[i] = 1;
        }
    } else {
        for( std::size_t i = startIndex; i != endIndex; ++i ) {
            destChannel[i] = value;
        }
    }
}

template <class T>
void set_array_range_to_value( T*, boost::uint8_t*, std::size_t, std::size_t, boost::none_t ) {}
} // namespace detail

/**
 *  Copy a named channel from an rle_level_set or an rle_voxel_field
 * into another rle_level_set or rle_voxel_field.
 *
 * @tparam ChannelDataType the type of the channel to copy.
 * @param[in,out] destField the rle_level_set or rle_voxel_field into which
 *		the channel will be copied.
 * @param destChannelName the name of a channel in destField that will
 *		receive the copied channel data.  If destField is an rle_level_set
 *		and this is "SignedDistance", then the source data will be copied
 *		into the level set's signed distance data.
 * @param destPopulatedChannelToCreate the name of a uint8[1] channel that
 *		will be created in destField to indicate which voxels were populated
 *		using defined data in the srcField.  This can be an empty string if
 *		you do not wish to create a populated channel.
 * @param srcField the rle_level_set or rle_voxel_field that contains the
 *		channel data to copy.
 * @param srcChannelName the name of a channel in srcField that will be
 *		copied into destField.  If the srcField is and rle_level_set and
 *		this is "SignedDistance", then the source data will be taken from
 *		the srcField's signed distance data.
 * @param defaultInsideValue undefined inside voxels in srcField will
 *		be assigned this value when they are written to destField.  Such
 *		undefined voxels will not be marked as populated.
 *		This can be boost::none if you do not want to assign a value for
 *		undefined inside voxels.
 * @param defaultOutsideValue undefined outside voxels in srcField will
 *		be assigned this value when they are written to destField.  Such
 *		undefined voxels will not be marked as populated.  This can be
 *		boost::none if you do not want to assign a value for undefined
 *		outside voxels.
 */
template <class ChannelDataType, class rle_channel_provider1, class rle_channel_provider2, class type_or_none1,
          class type_or_none2>
void copy_rle_channel( rle_channel_provider1& destField, const frantic::tstring& destChannelName,
                       const frantic::tstring& destPopulatedChannelToCreate, const rle_channel_provider2& srcField,
                       const frantic::tstring& srcChannelName, type_or_none1 defaultInsideValue,
                       type_or_none2 defaultOutsideValue ) {
    if( !detail::has_channel( srcField, srcChannelName ) ) {
        throw std::runtime_error( "copy_rle_channel Error: the source field is missing the specified channel \'" +
                                  frantic::strings::to_string( srcChannelName ) + "\'." );
    }

    std::pair<std::size_t, data_type_t> srcChannelType = detail::get_channel_type( srcField, srcChannelName );
    if( srcChannelType.first != frantic::channels::channel_data_type_traits<ChannelDataType>::arity() ||
        srcChannelType.second != frantic::channels::channel_data_type_traits<ChannelDataType>::data_type() ) {
        throw std::runtime_error(
            "copy_rle_channel Error: requested a copy of type " +
            frantic::strings::to_string( frantic::channels::channel_data_type_traits<ChannelDataType>::type_str() ) +
            ", but the source channel \'" + frantic::strings::to_string( srcChannelName ) + "\' is of type " +
            frantic::strings::to_string( channel_data_type_str( srcChannelType.first, srcChannelType.second ) ) + "." );
    }

    detail::add_channel( destField, destChannelName, srcChannelType.first, srcChannelType.second );

    boost::uint8_t* destPopulatedChannel = 0;
    if( !destPopulatedChannelToCreate.empty() ) {
        detail::add_channel<boost::uint8_t>( destField, destPopulatedChannelToCreate );
        destField.zero_channel( destPopulatedChannelToCreate );
        destPopulatedChannel =
            detail::get_channel_data_pointer<boost::uint8_t>( destField, destPopulatedChannelToCreate );
    }

    if( destField.size() == 0 )
        return;

    ChannelDataType* destChannel = detail::get_channel_data_pointer<ChannelDataType>( destField, destChannelName );

    const rle_index_spec& srcRIS = srcField.get_rle_index_spec();
    const rle_index_spec& destRIS = destField.get_rle_index_spec();

    if( srcField.size() == 0 ) {
        if( srcField.get_rle_index_spec().get_exterior_region_code() < -1 ) {
            detail::set_array_range_to_value( destChannel, destPopulatedChannel, 0, destRIS.data_size(),
                                              defaultInsideValue );
        } else if( srcField.get_rle_index_spec().get_exterior_region_code() == -1 ) {
            detail::set_array_range_to_value( destChannel, destPopulatedChannel, 0, destRIS.data_size(),
                                              defaultOutsideValue );
        }
    } else {
        const ChannelDataType* srcChannel =
            detail::get_channel_data_pointer<ChannelDataType>( srcField, srcChannelName );

        if( srcChannel == destChannel ) {
            if( destPopulatedChannel ) {
                for( std::size_t i = 0; i < destRIS.data_size(); ++i ) {
                    destPopulatedChannel[i] = 1;
                }
            }
        } else {
            const frantic::graphics::boundbox3& bounds = destRIS.outer_bounds();
            for( int z = bounds.zminimum(); z <= bounds.zmaximum(); ++z ) {
                for( int y = bounds.yminimum(); y <= bounds.ymaximum(); ++y ) {
                    rle_pairwise_run_iterator i( destRIS, destRIS.y_to_b( y ), destRIS.z_to_c( z ), srcRIS,
                                                 srcRIS.y_to_b( y ), srcRIS.z_to_c( z ) ),
                        ie;
                    // Go through all the sub-intervals shared between the destination voxel field and the level set
                    for( ; i != ie; ++i ) {
                        boost::int32_t destIndex = i.get_first_data_index();
                        boost::int32_t srcIndex = i.get_second_data_index();

                        if( destIndex >= 0 ) {
                            const boost::int32_t destIndexEnd = destIndex + i.get_xsize();

                            if( srcIndex < -1 ) {
                                detail::set_array_range_to_value( destChannel, destPopulatedChannel, destIndex,
                                                                  destIndexEnd, defaultInsideValue );
                            } else if( srcIndex == -1 ) {
                                detail::set_array_range_to_value( destChannel, destPopulatedChannel, destIndex,
                                                                  destIndexEnd, defaultOutsideValue );
                            } else if( srcIndex >= 0 ) {
                                for( ; destIndex != destIndexEnd; ++destIndex, ++srcIndex ) {
                                    destChannel[destIndex] = srcChannel[srcIndex];
                                    if( destPopulatedChannel )
                                        destPopulatedChannel[destIndex] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 *  Copy a named channel from an rle_level_set or an rle_voxel_field
 * into another rle_level_set or rle_voxel_field.
 *
 * @tparam ChannelDataType the type of the channel to copy.
 * @param[in,out] destField the rle_level_set or rle_voxel_field into which
 *		the channel will be copied.
 * @param destChannelName the name of a channel in destField that will
 *		receive the copied channel data.  If destField is an rle_level_set
 *		and this is "SignedDistance", then the source data will be copied
 *		into the level set's signed distance data.
 * @param destPopulatedChannelToCreate the name of a uint8[1] channel that
 *		will be created in destField to indicate which voxels were populated
 *		using defined data in the srcField.  This can be an empty string if
 *		you do not wish to create a populated channel.
 * @param srcField the rle_level_set or rle_voxel_field that contains the
 *		channel data to copy.
 * @param srcChannelName the name of a channel in srcField that will be
 *		copied into destField.  If the srcField is and rle_level_set and
 *		this is "SignedDistance", then the source data will be taken from
 *		the srcField's signed distance data.
 */
template <typename ChannelDataType, typename rle_channel_provider1, typename rle_channel_provider2>
void copy_rle_channel( rle_channel_provider1& destField, const frantic::tstring& destChannelName,
                       const frantic::tstring& destPopulatedChannelToCreate, const rle_channel_provider2& srcField,
                       const frantic::tstring& srcChannelName ) {
    copy_rle_channel<ChannelDataType, rle_channel_provider1, rle_channel_provider2>(
        destField, destChannelName, destPopulatedChannelToCreate, srcField, srcChannelName, boost::none, boost::none );
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
