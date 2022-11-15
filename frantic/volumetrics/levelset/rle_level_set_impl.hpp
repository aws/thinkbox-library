// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <boost/numeric/conversion/bounds.hpp>
#include <frantic/volumetrics/levelset/rle_general_block_iterator.hpp>
#include <vector>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

namespace frantic {
namespace volumetrics {
namespace levelset {

template <class ReconstructionFilter>
void rle_level_set::fill_box( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                              const frantic::graphics::boundbox3& voxelExtents,
                              std::vector<float>& outVoxelCornerValues ) const {
    frantic::graphics::size3 voxelSize = voxelExtents.size();
    outVoxelCornerValues.resize( voxelSize.volume() );

    float lsVoxelX =
        m_voxelCoordSystem.get_voxel_x_coord( destCoordSys.get_world_x_coord( voxelExtents.xminimum() + 0.5f ) ) - 0.5f;
    int xVoxelMinInit = (int)floorf( lsVoxelX - 1.f );

    float yFilterCoeff[4];
    float yzFilterCoeff[16];
    boost::int32_t yzDataIndexes[16];

    float* outVoxelCornerValuePointer = &outVoxelCornerValues[0];
    for( int z = voxelExtents.zminimum(); z <= voxelExtents.zmaximum(); ++z ) {
        for( int y = voxelExtents.yminimum(); y <= voxelExtents.ymaximum(); ++y ) {
            // Get the voxel YZ coordinate in the level set coordinate system
            float lsVoxelY = m_voxelCoordSystem.get_voxel_y_coord( destCoordSys.get_world_y_coord( y + 0.5f ) ) - 0.5f;
            float lsVoxelZ = m_voxelCoordSystem.get_voxel_z_coord( destCoordSys.get_world_z_coord( z + 0.5f ) ) - 0.5f;

            // Figure out the bounds we need in the YZ plane
            int yVoxelMin = (int)floorf( lsVoxelY - 1.f ), zVoxelMin = (int)floorf( lsVoxelZ - 1.f );
            int yVoxelMax = yVoxelMin + 3, zVoxelMax = zVoxelMin + 3;

            ////////////////////
            // Fill in the YZ reconstruction filter coefficients for this X scanline
            ////////////////////
            for( int i = 0; i < 4; ++i ) {
                yFilterCoeff[i] = (float)reconFilter( i + yVoxelMin - lsVoxelY );
            }
            for( int i = 0; i < 4; ++i ) {
                float zFilterCoeff = (float)reconFilter( i + zVoxelMin - lsVoxelZ );
                int iOffset = i << 2;
                for( int j = 0; j < 4; ++j ) {
                    yzFilterCoeff[iOffset + j] = zFilterCoeff * yFilterCoeff[j];
                }
            }

            ////////////////////
            // Create the block iterator, starting at the beginning X value
            ////////////////////
            rle_general_block_iterator iter( m_rleIndex, xVoxelMinInit, yVoxelMin, yVoxelMax, zVoxelMin, zVoxelMax );

            ////////////////////
            // Initialize the yzSliceBuffer
            ////////////////////
            float yzSliceBuffer[4];
            int yzSliceBufferNextIndex = 0;
            // This variable tracks how many slices have gone by with no change.  If >= 4 have gone by,
            // There's no need to sum up the filter along X, because all the inputs have the same value.
            int unchangedYZSliceCount = 0;

            // The first thing we need to do is fill up four values in our YZ slice buffer which we need for the
            // interpolation of the first output value.  An invariant during the loop below is that this buffer is
            // always full
            iter.get_indexes( yzDataIndexes );
            float lastReconstructedValue = 0;

            for( int i = 0; i < 16; ++i ) {
                lastReconstructedValue += yzFilterCoeff[i] * get_using_data_index( yzDataIndexes[i] );
            }

            yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
            yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;

            for( int i = 1; i < 4; ++i ) {

                // Advance the cursor to one greater X, and produce a new reconstructed value if the data indexes
                // changed at all
                if( iter.increment_x() ) {
                    iter.get_indexes( yzDataIndexes );
                    lastReconstructedValue = 0;
                    for( int i = 0; i < 16; ++i ) {
                        lastReconstructedValue += yzFilterCoeff[i] * get_using_data_index( yzDataIndexes[i] );
                    }
                    unchangedYZSliceCount = 1;
                } else {
                    ++unchangedYZSliceCount;
                }
                yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
                yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;
            }

            ////////////////////
            // Loop through to fill this output scanline
            ////////////////////

            int xVoxelMin = xVoxelMinInit;
            for( int x = voxelExtents.xminimum(); x <= voxelExtents.xmaximum(); ++x ) {

                // There are a number of cases which we have to deal with here.
                // 1) The voxel length of the output voxel coordinate system might be so large that all of the YZ slices
                //    need to be evaluated for each output value.
                // 2) It might be in between, so one or two value could be reused each time.
                // 3) It might be really small, so that the same set of YZ slices is used multiple times before
                // computing
                //    a new YZ slice.
                lsVoxelX = m_voxelCoordSystem.get_voxel_x_coord( destCoordSys.get_world_x_coord( x + 0.5f ) ) - 0.5f;
                // Update the yz slices to reflect the desired xVoxelMin
                int desiredXVoxelMin = (int)floorf( lsVoxelX - 1.f );

                if( desiredXVoxelMin > xVoxelMin ) {
                    // 4-iStart is how many values we are reusing from the last step.
                    int iStart = 4 - desiredXVoxelMin + xVoxelMin;
                    if( iStart < 0 )
                        iStart = 0;
                    for( int i = iStart; i < 4; ++i ) {
                        if( iter.jump_to_x( desiredXVoxelMin + i ) ) {
                            iter.get_indexes( yzDataIndexes );
                            lastReconstructedValue = 0;
                            for( int i = 0; i < 16; ++i ) {
                                lastReconstructedValue += yzFilterCoeff[i] * get_using_data_index( yzDataIndexes[i] );
                            }
                            unchangedYZSliceCount = 1;
                        } else {
                            ++unchangedYZSliceCount;
                        }
                        yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
                        yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;
                    }
                    xVoxelMin = desiredXVoxelMin;
                }

                // Apply the filter along X to get the final reconstructed values.
                if( unchangedYZSliceCount >= 4 ) {
                    // When there are 4 identical values, just copy one of them instead of doing the weighted sum.
                    *outVoxelCornerValuePointer++ = yzSliceBuffer[0];
                } else {
                    float result = 0;
                    for( int i = 0; i < 4; ++i ) {
                        float fortest2 = i + xVoxelMin - lsVoxelX;
                        float forest = (float)reconFilter( fortest2 );
                        result += yzSliceBuffer[( yzSliceBufferNextIndex + i ) & 0x03] * forest;
                    }
                    *outVoxelCornerValuePointer++ = result;
                }
            }
        }
    }
}

namespace detail {

template <class ReconstructionFilter>
void fill_box_mt_detail( const rle_level_set& rls, const voxel_coord_system& destCoordSys,
                         const ReconstructionFilter& reconFilter, const frantic::graphics::boundbox3& voxelExtents,
                         const frantic::graphics::boundbox3& totalVoxelExtents, float* outVoxelCornerValueStartPointer,
                         tbb::spin_mutex& /*m_mutex*/ ) {

    //{
    // tbb::spin_mutex::scoped_lock lock(m_mutex);
    // std::cout << "fill_box_mt_detail()" << std::endl;
    // std::cout << voxelExtents << std::endl << "of" << std::endl << totalVoxelExtents << std::endl;
    //}

    float lsVoxelX = rls.get_voxel_coord_system().get_voxel_x_coord(
                         destCoordSys.get_world_x_coord( voxelExtents.xminimum() + 0.5f ) ) -
                     0.5f;
    int xVoxelMinInit = (int)floorf( lsVoxelX - 1.f );

    float yFilterCoeff[4];
    float yzFilterCoeff[16];
    boost::int32_t yzDataIndexes[16];

    for( int z = voxelExtents.zminimum(); z < voxelExtents.zmaximum(); ++z ) {
        int zOffset = ( z - totalVoxelExtents.zminimum() ) * ( totalVoxelExtents.xsize() * totalVoxelExtents.ysize() );
        for( int y = voxelExtents.yminimum(); y < voxelExtents.ymaximum(); ++y ) {
            int yOffset = ( y - totalVoxelExtents.yminimum() ) * totalVoxelExtents.xsize();

            //{
            // tbb::spin_mutex::scoped_lock lock(m_mutex);
            // std::cout << "y,z: " << y << "," << z << std::endl;
            //}

            float* outVoxelCornerValuePointer = outVoxelCornerValueStartPointer + zOffset + yOffset;

            // Get the voxel YZ coordinate in the level set coordinate system
            float lsVoxelY =
                rls.get_voxel_coord_system().get_voxel_y_coord( destCoordSys.get_world_y_coord( y + 0.5f ) ) - 0.5f;
            float lsVoxelZ =
                rls.get_voxel_coord_system().get_voxel_z_coord( destCoordSys.get_world_z_coord( z + 0.5f ) ) - 0.5f;

            // Figure out the bounds we need in the YZ plane
            int yVoxelMin = (int)floorf( lsVoxelY - 1.f ), zVoxelMin = (int)floorf( lsVoxelZ - 1.f );
            int yVoxelMax = yVoxelMin + 3, zVoxelMax = zVoxelMin + 3;

            ////////////////////
            // Fill in the YZ reconstruction filter coefficients for this X scanline
            ////////////////////
            for( int i = 0; i < 4; ++i ) {
                yFilterCoeff[i] = (float)reconFilter( i + yVoxelMin - lsVoxelY );
            }
            for( int i = 0; i < 4; ++i ) {
                float zFilterCoeff = (float)reconFilter( i + zVoxelMin - lsVoxelZ );
                int iOffset = i << 2;
                for( int j = 0; j < 4; ++j ) {
                    yzFilterCoeff[iOffset + j] = zFilterCoeff * yFilterCoeff[j];
                }
            }

            ////////////////////
            // Create the block iterator, starting at the beginning X value
            ////////////////////
            rle_general_block_iterator iter( rls.get_rle_index_spec(), xVoxelMinInit, yVoxelMin, yVoxelMax, zVoxelMin,
                                             zVoxelMax );

            ////////////////////
            // Initialize the yzSliceBuffer
            ////////////////////
            float yzSliceBuffer[4];
            int yzSliceBufferNextIndex = 0;
            // This variable tracks how many slices have gone by with no change.  If >= 4 have gone by,
            // There's no need to sum up the filter along X, because all the inputs have the same value.
            int unchangedYZSliceCount = 0;

            // The first thing we need to do is fill up four values in our YZ slice buffer which we need for the
            // interpolation of the first output value.  An invariant during the loop below is that this buffer is
            // always full
            iter.get_indexes( yzDataIndexes );
            float lastReconstructedValue = 0;

            for( int i = 0; i < 16; ++i ) {
                lastReconstructedValue += yzFilterCoeff[i] * rls.get_using_data_index( yzDataIndexes[i] );
            }

            yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
            yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;

            for( int i = 1; i < 4; ++i ) {

                // Advance the cursor to one greater X, and produce a new reconstructed value if the data indexes
                // changed at all
                if( iter.increment_x() ) {
                    iter.get_indexes( yzDataIndexes );
                    lastReconstructedValue = 0;
                    for( int i = 0; i < 16; ++i ) {
                        lastReconstructedValue += yzFilterCoeff[i] * rls.get_using_data_index( yzDataIndexes[i] );
                    }
                    unchangedYZSliceCount = 1;
                } else {
                    ++unchangedYZSliceCount;
                }
                yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
                yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;
            }

            ////////////////////
            // Loop through to fill this output scanline
            ////////////////////

            int xVoxelMin = xVoxelMinInit;
            for( int x = voxelExtents.xminimum(); x <= voxelExtents.xmaximum(); ++x ) {

                // There are a number of cases which we have to deal with here.
                // 1) The voxel length of the output voxel coordinate system might be so large that all of the YZ slices
                //    need to be evaluated for each output value.
                // 2) It might be in between, so one or two value could be reused each time.
                // 3) It might be really small, so that the same set of YZ slices is used multiple times before
                // computing
                //    a new YZ slice.
                lsVoxelX =
                    rls.get_voxel_coord_system().get_voxel_x_coord( destCoordSys.get_world_x_coord( x + 0.5f ) ) - 0.5f;
                // Update the yz slices to reflect the desired xVoxelMin
                int desiredXVoxelMin = (int)floorf( lsVoxelX - 1.f );

                if( desiredXVoxelMin > xVoxelMin ) {
                    // 4-iStart is how many values we are reusing from the last step.
                    int iStart = 4 - desiredXVoxelMin + xVoxelMin;
                    if( iStart < 0 )
                        iStart = 0;
                    for( int i = iStart; i < 4; ++i ) {
                        if( iter.jump_to_x( desiredXVoxelMin + i ) ) {
                            iter.get_indexes( yzDataIndexes );
                            lastReconstructedValue = 0;
                            for( int i = 0; i < 16; ++i ) {
                                lastReconstructedValue +=
                                    yzFilterCoeff[i] * rls.get_using_data_index( yzDataIndexes[i] );
                            }
                            unchangedYZSliceCount = 1;
                        } else {
                            ++unchangedYZSliceCount;
                        }
                        yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
                        yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;
                    }
                    xVoxelMin = desiredXVoxelMin;
                }

                // Apply the filter along X to get the final reconstructed values.
                if( unchangedYZSliceCount >= 4 ) {
                    // When there are 4 identical values, just copy one of them instead of doing the weighted sum.
                    *outVoxelCornerValuePointer++ = yzSliceBuffer[0];
                } else {
                    float result = 0;
                    for( int i = 0; i < 4; ++i ) {
                        float fortest2 = i + xVoxelMin - lsVoxelX;
                        float forest = (float)reconFilter( fortest2 );
                        result += yzSliceBuffer[( yzSliceBufferNextIndex + i ) & 0x03] * forest;
                    }
                    *outVoxelCornerValuePointer++ = result;
                }
            }
        }
    }
}

// necessary for multithreading with tbb
template <class ReconstructionFilter>
class FillBoxBody {

    const rle_level_set& m_rls;
    const voxel_coord_system& m_destCoordSys;
    const ReconstructionFilter& m_reconFilter;
    const frantic::graphics::boundbox3& m_voxelExtents;
    float* m_outVoxelCornerValues;
    tbb::spin_mutex& m_mutex;

  public:
    void operator()( const tbb::blocked_range2d<int, int>& r ) const {

        //{
        //	tbb::spin_mutex::scoped_lock lock(m_mutex);
        //	std::cout << "operating on yz range: " << r.rows().begin() << "," << r.rows().end() << "  " <<
        // r.cols().begin()
        //<< "," << r.cols().end() << std::endl;
        //}

        // the 2d range is in the yz plane, so threads run on blocks of extents along the x axis
        // const frantic::graphics::boundbox3 voxelExtents(m_voxelExtents.xminimum(),
        //												m_voxelExtents.xmaximum(),
        //												r.begin(),
        //												r.end(),
        //												m_voxelExtents.zminimum(),
        //												m_voxelExtents.zminimum()+1
        //);
        const frantic::graphics::boundbox3 voxelExtents( m_voxelExtents.xminimum(), m_voxelExtents.xmaximum(),
                                                         r.rows().begin(), r.rows().end(), r.cols().begin(),
                                                         r.cols().end() );
        fill_box_mt_detail( m_rls, m_destCoordSys, m_reconFilter, voxelExtents, m_voxelExtents, m_outVoxelCornerValues,
                            m_mutex );
    }

    FillBoxBody( const rle_level_set& rls, const voxel_coord_system& destCoordSys,
                 const ReconstructionFilter& reconFilter, const frantic::graphics::boundbox3& voxelExtents,
                 float* outVoxelCornerValues, tbb::spin_mutex& mutex )
        : m_rls( rls )
        , m_destCoordSys( destCoordSys )
        , m_reconFilter( reconFilter )
        , m_voxelExtents( voxelExtents )
        , m_outVoxelCornerValues( outVoxelCornerValues )
        , m_mutex( mutex ) {}

    FillBoxBody& operator=( const FillBoxBody& rhs ) { return *( this ); } // unimplemented
};

} // namespace detail

template <class ReconstructionFilter>
void rle_level_set::fill_box_mt( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                                 const frantic::graphics::boundbox3& voxelExtents, float* outVoxelCornerValues ) const {
    // size3 voxelSize = voxelExtents.size();
    // outVoxelCornerValues.resize(voxelSize.volume());

    // std::cout << voxelExtents << std::endl;
    tbb::spin_mutex mutex;
    tbb::parallel_for( tbb::blocked_range2d<int, int>( voxelExtents.minimum().y, voxelExtents.maximum().y + 1,
                                                       voxelExtents.minimum().z, voxelExtents.maximum().z + 1 ),
                       // tbb::blocked_range2d<int,int>(voxelExtents.minimum().y,
                       //								voxelExtents.maximum().y+1,
                       //								voxelExtents.ysize(),
                       //								voxelExtents.minimum().z,
                       //								voxelExtents.maximum().z+1,
                       //								voxelExtents.zsize()),
                       detail::FillBoxBody<ReconstructionFilter>( *this, destCoordSys, reconFilter, voxelExtents,
                                                                  outVoxelCornerValues, mutex ),
                       tbb::auto_partitioner() );
}

/*********************************************************************************************************/
template <class ReconstructionFilter>
void rle_level_set::fill_plane( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                                const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                                std::vector<float>& outVoxelCornerValues ) const {
    using namespace frantic::graphics;

    fill_box( destCoordSys, reconFilter,
              boundbox3( vector3( voxelXYExtents.minimum().x, voxelXYExtents.minimum().y, voxelZ ),
                         vector3( voxelXYExtents.maximum().x, voxelXYExtents.maximum().y, voxelZ ) ),
              outVoxelCornerValues );
}

template <class ReconstructionFilter>
void rle_level_set::fill_plane_mt( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                                   const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                                   float* outVoxelCornerValues ) const {
    using namespace frantic::graphics;

    fill_box_mt( destCoordSys, reconFilter,
                 boundbox3( vector3( voxelXYExtents.minimum().x, voxelXYExtents.minimum().y, voxelZ ),
                            vector3( voxelXYExtents.maximum().x, voxelXYExtents.maximum().y, voxelZ ) ),
                 outVoxelCornerValues );
}

namespace detail {

/**
 *	A helper function to reconstruct channel data using a filter.
 *
 *	It returns true if there is defined data in this YZ place, otherwise returns false;
 *
 *	@param	channelAccessors	Accessors to the channel data to be filtered
 *	@param	yzDataIndexes		Indices to the channel data
 *	@param	yzFilterCoeff		Filter coefficients to be applied to the data at the above indices
 *	@param	filterSize			The number of data members to be filtered (length of above arrays)
 *	@param	reconstructedValues	The filter reconstructed values for each channel
 */
inline bool reconstruct_channel_data( std::vector<const_rle_channel_general_accessor>& channelAccessors,
                                      int* yzDataIndexes, float* yzFilterCoeff, int filterSize,
                                      std::vector<std::vector<char>>& reconstructedValues ) {

    // make some space for reconstruction values, will need a max of filterSize
    std::vector<char const*> data( filterSize );
    std::vector<float> weights( filterSize );

    for( size_t channelNum = 0; channelNum < channelAccessors.size(); ++channelNum ) {

        // build the data and weight arrays
        float weightSum = 0.f;
        int weightCount = 0;
        for( int i = 0; i < filterSize; ++i ) {
            // if( (unsigned int)yzDataIndexes[i] < (unsigned int)channelAccessors[channelNum].size() ) {
            if( yzDataIndexes[i] >= 0 ) {
                // std::cout << "yzDataIndexes[" << i << "]: " << yzDataIndexes[i] << std::endl;
                data[weightCount] = channelAccessors[channelNum].data( yzDataIndexes[i] );
                weights[weightCount] = yzFilterCoeff[i];
                weightSum += yzFilterCoeff[i];
                ++weightCount;
            }
        }

        if( weightCount == 0 ) {
            for( unsigned i = 0; i < channelAccessors.size(); ++i ) {
                memset( &( reconstructedValues[i][0] ), 0, channelAccessors[i].primitive_size() );
            }
            return false;
        }

        // reconstruct the channel data
        for( unsigned i = 0; i < weights.size(); ++i )
            weights[i] /= weightSum;

        channelAccessors[channelNum].get_weighted_sum_combine_function()( &weights[0], &data[0], weightCount,
                                                                          channelAccessors[channelNum].arity(),
                                                                          &( reconstructedValues[channelNum][0] ) );

        return true;
    }
    return false;
}

inline float find_last_reconstructed_value( const rle_level_set& rls, float lastReconstructedValue,
                                            float* yzFilterCoeff, std::vector<int> yzDataIndexes, int size ) {
    float sumOfCoeff = 0;
    std::vector<float> newYZFilterCoeff( 16 );
    int i;
    for( i = 0; i < size; ++i ) {
        if( yzDataIndexes[i] >= 0 )
            sumOfCoeff += yzFilterCoeff[i];
    }

    for( i = 0; i < size; ++i ) {
        if( yzDataIndexes[i] >= 0 )
            newYZFilterCoeff[i] = yzFilterCoeff[i] / sumOfCoeff;
    }

    for( i = 0; i < size; ++i ) {
        if( yzDataIndexes[i] >= 0 )
            lastReconstructedValue += newYZFilterCoeff[i] * rls.get_using_data_index( yzDataIndexes[i] );
    }
    return lastReconstructedValue;
}

} // namespace detail

/*********************************************************************************************************/
template <class ReconstructionFilter>
void rle_level_set::fill_sparse_plane_channel_data( const voxel_coord_system& destCoordSys,
                                                    const ReconstructionFilter& reconFilter,
                                                    const frantic::graphics2d::boundrect2& voxelExtents, int voxelZ,
                                                    std::vector<frantic::tstring>& channelNames,
                                                    std::vector<char*>& channelData,
                                                    frantic::volumetrics::rle_plane& outRLP ) const {
    outRLP.reset( voxelExtents );

    float lsVoxelZ = m_voxelCoordSystem.get_voxel_z_coord( destCoordSys.get_world_z_coord( voxelZ + 0.5f ) ) - 0.5f;
    float lsVoxelX =
        m_voxelCoordSystem.get_voxel_x_coord( destCoordSys.get_world_x_coord( voxelExtents.minimum().x + 0.5f ) ) -
        0.5f;
    int xVoxelMinInit = (int)floorf( lsVoxelX - 1.f );

    float yFilterCoeff[4];
    float yzFilterCoeff[16];
    std::vector<int> yzDataIndexes;

    // initialize all
    std::vector<const_rle_channel_general_accessor> channelAccessors;
    std::vector<char*> outChannelData;
    std::vector<std::vector<char>> yzSliceChannelBuffers;
    std::vector<std::vector<char>> lastReconstructedChannelValues;
    int signedDistanceChannel = -1;
    for( unsigned int i = 0; i < channelNames.size(); ++i ) {
        if( channelNames[i] == "SignedDistance" ) {
            signedDistanceChannel = i;
        } else {
            channelAccessors.push_back( get_channel_general_accessor( channelNames[i] ) );
            outChannelData.push_back( channelData[i] );
            yzSliceChannelBuffers.push_back(
                std::vector<char>( channelAccessors[channelAccessors.size() - 1].primitive_size() * 4 ) );
            lastReconstructedChannelValues.push_back(
                std::vector<char>( channelAccessors[channelAccessors.size() - 1].primitive_size() ) );
        }
    }

    for( int y = voxelExtents.minimum().y; y <= voxelExtents.maximum().y; ++y ) {
        // runs for the extent
        std::vector<std::pair<int, int>> runs;

        // Get the voxel YZ coordinate in the level set coordinate system
        float lsVoxelY = m_voxelCoordSystem.get_voxel_y_coord( destCoordSys.get_world_y_coord( y + 0.5f ) ) - 0.5f;

        // Figure out the bounds we need in the YZ plane
        int yVoxelMin = (int)floorf( lsVoxelY - 1.f ), zVoxelMin = (int)floorf( lsVoxelZ - 1.f );
        int yVoxelMax = yVoxelMin + 3, zVoxelMax = zVoxelMin + 3;

        ////////////////////
        // Fill in the YZ reconstruction filter coefficients for this X scanline
        ////////////////////
        for( int i = 0; i < 4; ++i ) {
            yFilterCoeff[i] = (float)reconFilter( i + yVoxelMin - lsVoxelY );
        }
        for( int i = 0; i < 4; ++i ) {
            float zFilterCoeff = (float)reconFilter( i + zVoxelMin - lsVoxelZ );
            int iOffset = i << 2;
            for( int j = 0; j < 4; ++j ) {
                yzFilterCoeff[iOffset + j] = zFilterCoeff * yFilterCoeff[j];
            }
        }

        ////////////////////
        // Create the block iterator, starting at the beginning X value
        ////////////////////
        rle_general_block_iterator iter( m_rleIndex, xVoxelMinInit, yVoxelMin, yVoxelMax, zVoxelMin, zVoxelMax );

        ////////////////////
        // Initialize the yzSliceBuffer
        ////////////////////
        float yzSliceBuffer[4];
        int yzSliceBufferNextIndex = 0;
        std::vector<bool> yzSliceBufferSignedDistanceDefined( 4 );
        std::vector<bool> yzSliceBufferOtherChannelDefined( 4 );
        std::vector<float> weightsForSignedDistance( 4 );
        std::vector<float> weightsForOtherChannel( 4 );

        // This variable tracks how many slices have gone by with no change.  If >= 4 have gone by,
        // There's no need to sum up the filter along X, because all the inputs have the same value.
        int unchangedYZSliceCount = 0;
        int xValueOfStartPoint = iter.current_x();

        int xOfLastDefinedDataSoFar = iter.current_x();

        // if no defined data at current X, jump to the defined data - 2
        if( !iter.is_currentX_has_defined_indexes() ) {
            if( iter.next_x() - 4 >= voxelExtents.minimum().x ) {
                iter.jump_to_x( iter.next_x() - 4 );
                if( iter.is_finished() ) {
                    outRLP.append_runs_by_extent( runs, y - voxelExtents.minimum().y, yzDataIndexes[0] );
                    continue;
                }
                xValueOfStartPoint = iter.next_x() - 2;
                xOfLastDefinedDataSoFar = iter.next_x();
            }
        }

        // The first thing we need to do is fill up four values in our YZ slice buffer which we need for the
        // interpolation of the first output value.  An invariant during the loop below is that this buffer is always
        // full
        float lastReconstructedValue = 0;
        iter.get_indexes( yzDataIndexes );

        // find out if there is defined data in this current YZ place
        yzSliceBufferSignedDistanceDefined[yzSliceBufferNextIndex] = iter.is_currentX_has_defined_indexes();

        if( yzSliceBufferSignedDistanceDefined[yzSliceBufferNextIndex] )
            lastReconstructedValue = detail::find_last_reconstructed_value( *this, lastReconstructedValue,
                                                                            yzFilterCoeff, yzDataIndexes, 16 );

        yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
        // cout << "lastReconstructedValue = " << lastReconstructedValue << endl;

        // rest of channels
        yzSliceBufferOtherChannelDefined[yzSliceBufferNextIndex] = detail::reconstruct_channel_data(
            channelAccessors, &yzDataIndexes[0], &yzFilterCoeff[0], 16, lastReconstructedChannelValues );
        for( unsigned int channelNum = 0; channelNum < channelAccessors.size(); ++channelNum )
            memcpy( &( yzSliceChannelBuffers[channelNum]
                                            [yzSliceBufferNextIndex * channelAccessors[channelNum].primitive_size()] ),
                    &( lastReconstructedChannelValues[channelNum][0] ), channelAccessors[channelNum].primitive_size() );
        // float a = *((float*)(&(lastReconstructedChannelValues[0][0])));
        // cout<< "the a = "<< a << endl;
        // std::cout << "yzSliceChannelBuffersdata: " <<
        // (*(float*)(&(yzSliceChannelBuffers[0][yzSliceBufferNextIndex*channelAccessors[0].primitive_size()]))) <<
        // std::endl;
        yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;

        for( int i = 1; i < 4; ++i ) {
            // Advance the cursor to one greater X, and produce a new reconstructed value if the data indexes changed at
            // all
            if( iter.increment_x() ) {
                if( iter.is_currentX_has_defined_indexes() ) {
                    xOfLastDefinedDataSoFar = iter.current_x();
                }
                iter.get_indexes( yzDataIndexes );

                // find out if there is defined data in this current YZ place
                yzSliceBufferSignedDistanceDefined[yzSliceBufferNextIndex] = iter.is_currentX_has_defined_indexes();
                lastReconstructedValue = 0;
                if( yzSliceBufferSignedDistanceDefined[yzSliceBufferNextIndex] )
                    lastReconstructedValue = detail::find_last_reconstructed_value( *this, lastReconstructedValue,
                                                                                    yzFilterCoeff, yzDataIndexes, 16 );

                // cout << "lastReconstructedValue = " << lastReconstructedValue << endl;

                // rest of channels
                yzSliceBufferOtherChannelDefined[yzSliceBufferNextIndex] = detail::reconstruct_channel_data(
                    channelAccessors, &yzDataIndexes[0], &yzFilterCoeff[0], 16, lastReconstructedChannelValues );

                // float a = *((float*)(&(lastReconstructedChannelValues[0][0])));
                // cout<< "the a = "<< a << endl;

                unchangedYZSliceCount = 1;
            } else {
                ++unchangedYZSliceCount;
            }
            yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
            for( unsigned int channelNum = 0; channelNum < channelAccessors.size(); ++channelNum ) {
                // std::cout << "channel " << channelNum << std::endl;
                // std::cout << "copying " << *( (float*)&(lastReconstructedChannelValues[channelNum][0])) << std::endl;
                memcpy( &( yzSliceChannelBuffers[channelNum][yzSliceBufferNextIndex *
                                                             channelAccessors[channelNum].primitive_size()] ),
                        &( lastReconstructedChannelValues[channelNum][0] ),
                        channelAccessors[channelNum].primitive_size() );
                // std::cout << "yzSliceChannelBuffersdata: " <<
                // (*(float*)(&(yzSliceChannelBuffers[0][yzSliceBufferNextIndex*channelAccessors[0].primitive_size()])))
                // << std::endl;
            }
            yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;
        }

        int offset = voxelExtents.xsize() * ( y - voxelExtents.minimum().y ) - voxelExtents.minimum().x;

        ////////////////////
        // Loop through to fill this output scanline
        ////////////////////

        bool waitingForNewRun = true;
        int xVoxelMin = xVoxelMinInit;
        // int x = voxelExtents.minimum().x;

        int x = (int)floorf(
            destCoordSys.get_voxel_x_coord( m_voxelCoordSystem.get_world_x_coord( (float)xValueOfStartPoint ) ) );

        if( x < voxelExtents.minimum().x )
            x = voxelExtents.minimum().x;

        while( x <= voxelExtents.maximum().x && x >= voxelExtents.minimum().x ) {
            if( waitingForNewRun ) {
                runs.push_back( std::pair<int, int>( offset + x, -1 ) );
                waitingForNewRun = false;
            }
            // std::cout << "running" << std::endl;
            //  There are a number of cases which we have to deal with here.
            //  1) The voxel length of the output voxel coordinate system might be so large that all of the YZ slices
            //     need to be evaluated for each output value.
            //  2) It might be in between, so one or two value could be reused each time.
            //  3) It might be really small, so that the same set of YZ slices is used multiple times before computing
            //     a new YZ slice.
            lsVoxelX = m_voxelCoordSystem.get_voxel_x_coord( destCoordSys.get_world_x_coord( x + 0.5f ) ) - 0.5f;
            // Update the yz slices to reflect the desired xVoxelMin
            int desiredXVoxelMin = (int)floorf( lsVoxelX - 1.f );

            if( desiredXVoxelMin > xVoxelMin ) {
                // 4-iStart is how many values we are reusing from the last step.
                int iStart = 4 - desiredXVoxelMin + xVoxelMin;
                if( iStart < 0 )
                    iStart = 0;
                for( int i = iStart; i < 4; ++i ) {
                    if( iter.jump_to_x( desiredXVoxelMin + i ) ) {
                        if( iter.is_currentX_has_defined_indexes() ) {
                            xOfLastDefinedDataSoFar = iter.current_x();
                        }
                        iter.get_indexes( yzDataIndexes );

                        // find out if there is defined data in this current YZ place
                        yzSliceBufferSignedDistanceDefined[yzSliceBufferNextIndex] =
                            iter.is_currentX_has_defined_indexes();
                        lastReconstructedValue = 0;
                        if( yzSliceBufferSignedDistanceDefined[yzSliceBufferNextIndex] )
                            lastReconstructedValue = detail::find_last_reconstructed_value(
                                *this, lastReconstructedValue, yzFilterCoeff, yzDataIndexes, 16 );

                        // rest of channels
                        yzSliceBufferOtherChannelDefined[yzSliceBufferNextIndex] =
                            detail::reconstruct_channel_data( channelAccessors, &yzDataIndexes[0], &yzFilterCoeff[0],
                                                              16, lastReconstructedChannelValues );

                        unchangedYZSliceCount = 1;
                    } else {
                        ++unchangedYZSliceCount;
                    }
                    yzSliceBuffer[yzSliceBufferNextIndex] = lastReconstructedValue;
                    for( unsigned int channelNum = 0; channelNum < channelAccessors.size(); ++channelNum )
                        memcpy( &( yzSliceChannelBuffers[channelNum][yzSliceBufferNextIndex *
                                                                     channelAccessors[channelNum].primitive_size()] ),
                                &( lastReconstructedChannelValues[channelNum][0] ),
                                channelAccessors[channelNum].primitive_size() );
                    yzSliceBufferNextIndex = ( yzSliceBufferNextIndex + 1 ) & 0x03;
                }
                xVoxelMin = desiredXVoxelMin;
            }

            // Apply the filter along X to get the final reconstructed values.
            if( unchangedYZSliceCount >= 4 ) {
                // When there are 4 identical values, just copy one of them instead of doing the weighted sum.
                if( signedDistanceChannel >= 0 )
                    ( (float*)( channelData[signedDistanceChannel] ) )[offset + x] = yzSliceBuffer[0];

                for( size_t channelNum = 0; channelNum < channelAccessors.size(); ++channelNum ) {
                    memcpy( outChannelData[channelNum] + ( offset + x ) * channelAccessors[channelNum].primitive_size(),
                            &( yzSliceChannelBuffers[channelNum][0] ), channelAccessors[channelNum].primitive_size() );
                }

            } else {
                float result = 0;

                // build the filter weights
                std::vector<float> weightsForSignedDistance( 4 );
                std::vector<float> weightsForOtherChannel( 4 );
                float weightsSumForSignedDistance = 0;
                float weightsSumForOtherChannel = 0;

                float w[4];
                for( int i = 0; i < 4; ++i )
                    w[i] = (float)reconFilter( i + xVoxelMin - lsVoxelX );

                // determine what is the sum of total weight
                for( int i = 0; i < 4; ++i ) {
                    if( yzSliceBufferSignedDistanceDefined[( yzSliceBufferNextIndex + i ) & 0x03] == true )
                        weightsSumForSignedDistance += w[i];

                    if( yzSliceBufferOtherChannelDefined[( yzSliceBufferNextIndex + i ) & 0x03] == true )
                        weightsSumForOtherChannel += w[i];
                }

                // renomalize weights
                if( weightsSumForSignedDistance != 0 ) {
                    for( int i = 0; i < 4; ++i ) {
                        if( yzSliceBufferSignedDistanceDefined[( yzSliceBufferNextIndex + i ) & 0x03] == true )
                            weightsForSignedDistance[i] = w[i] / weightsSumForSignedDistance;

                        if( yzSliceBufferOtherChannelDefined[( yzSliceBufferNextIndex + i ) & 0x03] == true )
                            weightsForOtherChannel[i] = w[i] / weightsSumForOtherChannel;
                    }
                }

                // reconstruct the signed distance
                for( int i = 0; i < 4; ++i )
                    result += yzSliceBuffer[( yzSliceBufferNextIndex + i ) & 0x03] * weightsForSignedDistance[i];

                if( signedDistanceChannel >= 0 )
                    ( (float*)( channelData[signedDistanceChannel] ) )[offset + x] = result;

                // reconstruct the flexible channel data
                // detail::reconstruct_channel_data(channelAccessors, &yzDataIndexes[0], &yzFilterCoeff[0], 16,
                // lastReconstructedChannelValues);
                for( size_t channelNum = 0; channelNum < channelAccessors.size(); ++channelNum ) {
                    // std::cout << "reconstructing for channel : " << channelNum << std::endl;

                    float* a = (float*)( &yzSliceChannelBuffers[channelNum][0] );

                    std::vector<const char*> data( 4 );
                    for( size_t i = 0; i < data.size(); ++i ) {
                        data[i] = (const char*)( a + ( ( yzSliceBufferNextIndex + i ) & 0x03 ) );
                        // std::cout << "data[" << i << "] " << *((float*)(data[i])) << std::endl;
                        w[i] = weightsForOtherChannel[i];
                    }

                    channelAccessors[channelNum].get_weighted_sum_combine_function()(
                        w, &data[0], 4, channelAccessors[channelNum].arity(),
                        outChannelData[channelNum] + ( offset + x ) * channelAccessors[channelNum].primitive_size() );
                    //	memcpy(outChannelData[channelNum] + (offset+x)*channelAccessors[channelNum].primitive_size(),
                    //		   &(yzSliceChannelBuffers[channelNum][0]),
                    //		   channelAccessors[channelNum].primitive_size());
                }

            } // if
            if( xOfLastDefinedDataSoFar + 3 >= iter.current_x() || xOfLastDefinedDataSoFar + 3 >= iter.next_x() - 2 )
                ++x;
            else {
                if( !waitingForNewRun ) {
                    // If we did, end the current run.
                    runs.back().second = offset + x - 1;
                    waitingForNewRun = true;
                }
                x = (int)floorf(
                    destCoordSys.get_voxel_x_coord( m_voxelCoordSystem.get_world_x_coord( iter.next_x() - 2.f ) ) );
                xOfLastDefinedDataSoFar = iter.next_x();
            }
        } // x

        // If we hit the end of the dest coord system, close off the last run
        if( !waitingForNewRun ) {
            runs.back().second = offset + x - 1;
            waitingForNewRun = false;
        }
        // Add the runs to the rle plane
        outRLP.append_runs_by_extent( runs, y - voxelExtents.minimum().y, -1 );
    } // y
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
