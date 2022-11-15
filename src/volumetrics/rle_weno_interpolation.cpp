// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/rle_weno_interpolation.hpp>

//#include <frantic/fluids/rle_voxel_field.hpp>

namespace frantic {
namespace volumetrics {

using namespace frantic::graphics;
using namespace frantic::volumetrics::levelset;

void staggered_weno3_debug_dump(
    const boost::int32_t* const dataIndices, const frantic::graphics::size3& boxSize,
    const frantic::graphics::vector3& currentXYZMin,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& /*indicatorAccessor1*/,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& /*ndicatorAccessor2*/,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& velocityAccessor,
    frantic::graphics::vector3f voxelLookup ) {

    vector3 voxelFloor = vector3::from_floor( voxelLookup );
    vector3f voxelFloorFloat( voxelFloor );

    vector3f velocity;

    vector3 xyz, abc;
    vector3f weight;

    // for each component velocity we need to do a 3d weno3 interpolation
    for( int face = 2; face < 3; ++face ) {
        // compute the basic x -xi weights for the polynomials
        for( int i = 0; i < 3; ++i ) {
            float offsetLookup = ( i == face ? voxelLookup[i] : voxelLookup[i] - 0.5f );
            // std::cout << "\t\t\tfloor: " << floor( voxelLookup[i] - offset[i] ) << std::endl;
            xyz[i] = static_cast<int>( floor( offsetLookup ) );
            // these are the local xyz coordinates of the data index block
            abc[i] = xyz[i] - currentXYZMin[i];

            // std::cout << i << std::endl;
            weight[i] = voxelLookup[i] - xyz[i];
        }

        boost::int32_t originIndex =
            dataIndices[abc.x + abc.y * boxSize.xsize() + abc.z * boxSize.xsize() * boxSize.ysize()];
        // just in case we end up looking outside of the fluid domain
        if( originIndex < 0 ) {
            velocity[face] = 0.f;
            continue;
        }

        // FF_LOG(debug)  << "indicator(1):" << indicatorAccessor1[originIndex] << "\tindicator(2):" <<
        // indicatorAccessor2[originIndex] << std::endl;
        //  compute the 1d interpolations along x first, as this is the most we need to do and should be the fastest
        for( int z = -1; z < 3; ++z ) {
            for( int y = -1; y < 3; ++y ) {
                for( int x = -1; x < 3; ++x ) {

                    boost::int32_t index = dataIndices[abc.x + x + ( abc.y + y ) * boxSize.xsize() +
                                                       ( abc.z + z ) * boxSize.xsize() * boxSize.ysize()];

                    if( index < 0 ) {
                        logging::error << "#UNDEF  ";
                    } else {
                        logging::error << "[" << index << "]:" << velocityAccessor[index] << " ";
                    }
                }
                logging::error << std::endl;
            }
            logging::error << std::endl;
        }
    }
}

frantic::graphics::vector3f staggered_weno3_lookup(
    const boost::int32_t* const dataIndices, const frantic::graphics::size3& boxSize,
    const frantic::graphics::vector3& currentXYZMin,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& /*indicatorXAccessor1*/,
    frantic::volumetrics::levelset::const_rle_channel_accessor<float>& /*indicatorXAccessor2*/,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f>& velocityAccessor,
    frantic::graphics::vector3f voxelLookup ) {

    vector3 voxelFloor = vector3::from_floor( voxelLookup );
    vector3f voxelFloorFloat( voxelFloor );

    vector3f velocity;

    vector3 xyz, abc;
    vector3f weight;

    float xlerps[16];
    float ylerps[4];

    // FF_LOG(debug) << "lookup: " << voxelLookup << std::endl;
    // FF_LOG(debug) << "CurrentMin: " << currentXYZMin << std::endl;

    bool log = false; // logging::get_logging_level() == 4;
    // for each component velocity we need to do a 3d weno3 interpolation
    for( int face = 0; face < 3; ++face ) {
        // compute the basic x -xi weights for the polynomials
        for( int i = 0; i < 3; ++i ) {
            float offsetLookup = ( i == face ? voxelLookup[i] : voxelLookup[i] - 0.5f );
            // std::cout << "\t\t\tfloor: " << floor( voxelLookup[i] - offset[i] ) << std::endl;
            xyz[i] = static_cast<int>( floor( offsetLookup ) );
            // these are the local xyz coordinates of the data index block
            abc[i] = xyz[i] - currentXYZMin[i];

            // std::cout << i << std::endl;
            weight[i] = offsetLookup - xyz[i];

            if( weight[i] > 1.f ) {
                logging::error << "i:" << i << "    face: " << face << std::endl;
                logging::error << "lookup: " << voxelLookup << std::endl;
                logging::error << "CurrentMin: " << currentXYZMin << std::endl;
                logging::error << "abc:" << abc << std::endl;
                logging::error << "xyz:" << xyz << "  offset:" << offsetLookup << std::endl;
                logging::error << "weight:" << weight << std::endl;
                logging::error << "face:" << face << "  comp:" << i << std::endl;

                throw std::runtime_error( "miscalculated weight! " );
            }
        }

        int localIndex = abc.x + abc.y * boxSize.xsize() + abc.z * boxSize.xsize() * boxSize.ysize();
        if( localIndex > boxSize.volume() ) {
            logging::error << "lookup: " << voxelLookup << std::endl;
            logging::error << "CurrentMin: " << currentXYZMin << std::endl;
            logging::error << "abc:" << abc << std::endl;
            logging::error << "xyz:" << xyz << std::endl;
            logging::error << "boxSize: " << boxSize << std::endl;
            throw std::runtime_error( "index out of bounds" );
        }
        // boost::int32_t originIndex = dataIndices[ localIndex ];
        //  just in case we end up looking outside of the fluid domain
        // if( originIndex < 0 ) {
        //	FF_LOG(debug)  << "Lookup Outside Fluid domain : " << voxelLookup << " ( abcCoord=" << abc << ",
        //dataIndex=" << originIndex << ")" << std::endl; 	velocity[face] = 0.f; 	continue;
        // }

        // FF_LOG(debug)  << "\tweight:" << weight << std::endl;
        // FF_LOG(debug)  << "indicator(1):" << indicatorAccessor1[originIndex] << "\tindicator(2):" <<
        // indicatorAccessor2[originIndex] << std::endl;
        int yUndefined[16];
        int yUndefinedCount[4] = { 0, 0, 0, 0 };

        // compute the 1d interpolations along x first, as this is the most we need to do and should be the fastest
        for( int z = -1; z < 3; ++z ) {
            for( int y = -1; y < 3; ++y ) {
                float velocities[4];
                int undefinedCount = 0;
                int undefined[4] = { 0, 0, 0, 0 };
                int resultIndex = ( y + 1 ) + ( z + 1 ) * 4;
                for( int x = -1; x < 3; ++x ) {

                    boost::int32_t index = dataIndices[abc.x + x + ( abc.y + y ) * boxSize.xsize() +
                                                       ( abc.z + z ) * boxSize.xsize() * boxSize.ysize()];

                    if( index > -1 ) {
                        // FF_LOG(debug) << "\tFACE:" << face << "\t" << vector3(x+1,y+1,z+1) <<"\tdataIndex=" << index
                        // << std::endl;

                        if( index > (int)velocityAccessor.size() ) {

                            logging::error << "lookup: " << voxelLookup << std::endl;
                            logging::error << "CurrentMin: " << currentXYZMin << std::endl;
                            logging::error << "abc:" << abc << std::endl;
                            logging::error << "xyz:" << vector3( x, y, z ) << std::endl;
                            logging::error << "boxSize: " << boxSize << std::endl;

                            logging::error << "index: " << index << " face: " << face << std::endl;
                            throw std::runtime_error( "weno3 error - bad index when collection x lookups" );
                        }

                        undefined[x + 1] = 0;
                        velocities[x + 1] = velocityAccessor[index][face];
                    } else {
                        ++undefinedCount;
                        undefined[x + 1] = 1;
                        velocities[x + 1] =
                            -987789.f; //-std::numeric_limits<float>::max();  // ? possibly handle this case better...
                    }
                }

                if( undefinedCount == 4 ) {

                    xlerps[resultIndex] = 0.f;
                    yUndefined[resultIndex] = 1;
                    ++yUndefinedCount[z + 1];
                    if( log && face == 0 && z == 2 ) {
                        logging::debug << "[y,z]: (" << y << ", " << z << ") - Undefined!" << std::endl;
                    }
                }
                // if all four samples are undefined, we should drop the whole sample
                else if( undefinedCount >= 0 ) {
                    // int sweepCount = (undefinedCount+1)/2;

                    // FF_LOG(debug) << "\tUndefined x Count: " << undefinedCount << std::endl;
                    if( logging::is_logging_debug() ) {
                        for( int i = 0; i < 4; ++i ) {
                            logging::debug << "\t" << velocities[i];
                        }
                        logging::debug << std::endl;
                    }

                    for( int iter = 0; iter < undefinedCount; ++iter ) {
                        for( int i = 3; i >= 0; --i ) {
                            // for( int  i=0; i<4; ++i ) {
                            if( log && face == 0 && z == 0 ) {
                                logging::debug << "\n**(" << y << ", " << z << ")"
                                               << "**\ni: " << i << std::endl;
                                for( int j = 0; j < 4; ++j ) {
                                    logging::debug << "\t" << velocities[j] << std::endl;
                                }

                                logging::debug << "\n undefined:" << std::endl;
                                for( int j = 0; j < 4; ++j ) {
                                    logging::debug << "\t" << undefined[j] << std::endl;
                                }

                                logging::debug << "****\n" << std::endl;
                            }

                            if( undefined[i] == 0 )
                                continue;

                            float v1 = ( i > 0 ) ? velocities[i - 1] : 0.f;
                            float v2 = ( i < 3 ) ? velocities[i + 1] : 0.f;
                            float result = 0.f;
                            int weight = 0;

                            if( i > 0 && undefined[i - 1] == 0 ) {
                                result += v1;
                                ++weight;
                            }
                            if( i < 3 && undefined[i + 1] == 0 ) {
                                result += v2;
                                ++weight;
                            }
                            if( weight != 0 ) {
                                velocities[i] = result / static_cast<float>( weight );
                                undefined[i] = 0;
                                --undefinedCount;
                            }
                        }
                    }

                    // if( log && face==0 ) {
                    //	logging::set_logging_level(5);
                    //	FF_LOG(debug) << "[y,z]: (" <<  y << ", " << z << ")" << std::endl;
                    // }

                    yUndefined[resultIndex] = 0;

                    // if( originIndex > 0 ) {
                    //	logging::error << "indicator 1: " << indicatorXAccessor1[originIndex] << " indicator 2: " <<
                    // indicatorXAccessor2[originIndex] << std::endl; 	xlerps[resultIndex] = weno3_interpolant(
                    // velocities, weight.x, indicatorXAccessor1[originIndex], indicatorXAccessor2[originIndex]);
                    // }
                    // else {
                    //	xlerps[resultIndex] = weno3_interpolant( velocities, weight.x );
                    // }

                    xlerps[resultIndex] = weno3_interpolant( velocities, weight.x );

                    // if( log && face ==0 && z ==2 ) {
                    //	logging::set_logging_level(3);
                    // }
                }

                // FF_LOG(debug) << "\txlerp[" <<  (y+1) + (z+1)*4 << "](" << y << "," << z << "):" << xlerps[(y+1) +
                // (z+1)*4]
                // << std::endl;
            }
        }

        int zUndefined[4] = { 0, 0, 0, 0 };
        int zUndefinedCount = 0;

        // compute the next set of 4 interpolations along y
        for( int z = 0; z < 4; ++z ) {
            if( yUndefinedCount[z] == 4 ) {
                ++zUndefinedCount;
                zUndefined[z] = 1;
                ylerps[z] = 0.f;
                // std::cout << "Voxel Lookup: " << voxelLookup << std::endl;

                // std::stringstream ss;
                //
                // for( int i=0;i<4;++i){
                //	ss << "yUndefined[" << i << "]: "  << yUndefined[i+z*4] << std::endl;
                // }

                // for( int i=0;i<4;++i){
                //	ss << "xlerps[" << i << "]: "  << xlerps[i+z*4] << std::endl;
                // }
                //
                // std::cout << ss.str() << std::endl;

                // throw std::runtime_error( "Unable to get a defined value during z WENO3 interpolation for z:" +
                // boost::lexical_cast<std::string>(z));
            }

            if( yUndefinedCount[z] > 0 ) {
                int zOffset = z * 4;

                for( int iter = 0; iter < yUndefinedCount[z]; ++iter ) {
                    if( log && face == 0 && z == 3 ) {

                        logging::debug << "\n****\ni: " << iter << std::endl;
                        for( int i = 0; i < 4; ++i ) {
                            int index = i + zOffset;
                            logging::debug << "\t" << xlerps[index] << std::endl;
                        }

                        logging::debug << "\n yUndefined:" << std::endl;
                        for( int i = 0; i < 4; ++i ) {
                            logging::debug << "\t" << yUndefined[i] << std::endl;
                        }

                        logging::debug << "****\n" << std::endl;
                    }
                    for( int i = 3; i >= 0; --i ) {
                        // for( int  i=0; i<4; ++i ) {

                        int index = i + zOffset;
                        if( yUndefined[index] == 0 )
                            continue;
                        float v1 = ( i > 0 ) ? xlerps[index - 1] : 0.f;
                        float v2 = ( i < 3 ) ? xlerps[index + 1] : 0.f;
                        float result = 0.f;
                        int weight = 0;

                        if( i > 0 && yUndefined[index - 1] == 0 ) {
                            result += v1;
                            ++weight;
                        }
                        if( i < 3 && yUndefined[index + 1] == 0 ) {
                            result += v2;
                            ++weight;
                        }
                        if( weight != 0 ) {
                            xlerps[index] = result / static_cast<float>( weight );
                            yUndefined[index] = 0;
                            --yUndefinedCount[z];
                        }
                    }
                }
                if( log && face == 0 && z == 3 ) {

                    logging::debug << "\n**** POST ***** \n" << std::endl;
                    for( int i = 0; i < 4; ++i ) {
                        int index = i + zOffset;
                        logging::debug << "\t" << xlerps[index] << std::endl;
                    }
                    logging::debug << "****\n" << std::endl;
                }
            }

            // if( log && face ==0 && z ==3 ) {
            //	logging::set_logging_level(5);
            //	FF_LOG(debug) << "\n[z]: (" << z << ")" << std::endl;
            //
            // }
            ylerps[z] = weno3_interpolant( xlerps + z * 4, weight.y );

            // if( log && face ==0 && z ==3 ) {
            //	logging::set_logging_level(3);
            // }

            // FF_LOG(debug) << "\t ylerp[" <<  z<< "]:" << ylerps[z] << std::endl;
        }

        if( zUndefinedCount == 4 ) {
            std::cout << "Voxel Lookup: " << voxelLookup << std::endl;

            std::stringstream ss;

            for( int i = 0; i < 4; ++i ) {
                ss << "zUndefined[" << i << "]: " << zUndefined[i] << std::endl;
            }

            for( int i = 0; i < 4; ++i ) {
                ss << "ylerps[" << i << "]: " << ylerps[i] << std::endl;
            }

            std::cout << ss.str() << std::endl;

            throw std::runtime_error(
                "Unable to get a defined value during for the final z WENO3 interpolation for z" );
        }

        if( zUndefinedCount > 0 ) {
            for( int iter = 0; iter < zUndefinedCount; ++iter ) {
                // for( int  i=0; i<4; ++i ) {
                for( int i = 3; i >= 0; --i ) {

                    int index = i;
                    if( zUndefined[index] == 0 )
                        continue;

                    float v1 = ( i > 0 ) ? ylerps[index - 1] : 0.f;
                    float v2 = ( i < 3 ) ? ylerps[index + 1] : 0.f;
                    float result = 0.f;
                    int weight = 0;

                    if( i > 0 && zUndefined[index - 1] == 0 ) {
                        result += v1;
                        ++weight;
                    }
                    if( i < 3 && zUndefined[index + 1] == 0 ) {
                        result += v2;
                        ++weight;
                    }
                    if( weight != 0 ) {
                        ylerps[index] = result / static_cast<float>( weight );
                        zUndefined[index] = 0;
                        --zUndefinedCount;
                    }
                }
            }
        }

        // if( log && face ==0 ) {
        //	logging::set_logging_level(5);
        //	FF_LOG(debug) << "final lerp:\n" ;
        // }

        // velocity[face] = (ylerps[1] + ylerps[2])*.5f;
        //  compute the final interpolation on z
        //  save the resulting value
        velocity[face] = weno3_interpolant( ylerps, weight.z );

        // FF_LOG(debug) << "\t zlerp[" <<  face << "]:" << velocity[face] << std::endl;
        // if(log &&  face ==0 ) {
        //	logging::set_logging_level(3);
        // }
    }

    return velocity;
}

float weno3_signed_distance_lookup( const boost::int32_t* const dataIndices, const frantic::graphics::size3& boxSize,
                                    const frantic::graphics::vector3& currentXYZMin,
                                    const std::vector<float>& distanceData, frantic::graphics::vector3f voxelLookup ) {
    frantic::graphics::vector3 voxelFloor = vector3::from_floor( voxelLookup );
    frantic::graphics::vector3f voxelFloorFloat( voxelFloor );
    frantic::graphics::vector3 xyz, abc;
    frantic::graphics::vector3f weight;

    // FF_LOG(debug) << "lookup: " << voxelLookup << std::endl;
    // FF_LOG(debug) << "CurrentMin: " << currentXYZMin << std::endl;

    // compute the basic x -xi weights for the polynomials
    for( int i = 0; i < 3; ++i ) {
        float offsetLookup = voxelLookup[i] - 0.5f;
        // std::cout << "\t\t\tfloor: " << floor( voxelLookup[i] - offset[i] ) << std::endl;
        xyz[i] = static_cast<int>( floor( offsetLookup ) );
        // these are the local xyz coordinates of the data index block
        abc[i] = xyz[i] - currentXYZMin[i];

        // std::cout << i << std::endl;
        weight[i] = offsetLookup - xyz[i];

        if( weight[i] > 1.f || weight[i] < 0.f ) {
            logging::error << "lookup: " << voxelLookup << std::endl;
            logging::error << "CurrentMin: " << currentXYZMin << std::endl;
            logging::error << "abc:" << abc << std::endl;
            logging::error << "xyz:" << xyz << "  offset:" << offsetLookup << std::endl;
            logging::error << "weight:" << weight << std::endl;
            throw std::runtime_error( "miscalculated weight! " );
        }
    }

    int localIndex = abc.x + abc.y * boxSize.xsize() + abc.z * boxSize.xsize() * boxSize.ysize();
    if( localIndex > boxSize.volume() ) {

        logging::error << "localIndex: " << localIndex << std::endl;
        logging::error << "lookup: " << voxelLookup << std::endl;
        logging::error << "CurrentMin: " << currentXYZMin << std::endl;
        logging::error << "abc:" << abc << std::endl;
        logging::error << "xyz:" << xyz << std::endl;
        logging::error << "boxSize: " << boxSize << std::endl;
        throw std::runtime_error( "index out of bounds" );
    }

    // boost::int32_t originIndex = dataIndices[ localIndex ];
    //  just in case we end up looking outside of the fluid domain
    // if( originIndex < 0 ) {
    //	FF_LOG(debug)  << "Lookup Outside Fluid domain : " << voxelLookup << " ( abcCoord=" << abc << ", dataIndex=" <<
    // originIndex << ")" << std::endl; 	velocity[face] = 0.f; 	continue;
    // }

    // FF_LOG(debug)  << "\tweight:" << weight << std::endl;
    // FF_LOG(debug)  << "indicator(1):" << indicatorAccessor1[originIndex] << "\tindicator(2):" <<
    // indicatorAccessor2[originIndex] << std::endl;
    float xlerps[16];
    float ylerps[4];

    // compute the 1d interpolations along x first, as this is the most we need to do and should be the fastest
    for( int z = -1; z < 3; ++z ) {
        for( int y = -1; y < 3; ++y ) {
            float values[4];

            int resultIndex = ( y + 1 ) + ( z + 1 ) * 4;

            for( int x = -1; x < 3; ++x ) {
                boost::int32_t blockIndex =
                    abc.x + x + ( abc.y + y ) * boxSize.xsize() + ( abc.z + z ) * boxSize.xsize() * boxSize.ysize();

                if( blockIndex > boxSize.volume() ) {
                    logging::error << "blockIndex: " << blockIndex << std::endl;
                    logging::error << "x: " << x << " y:" << y << " z:" << z << std::endl;
                    logging::error << "lookup: " << voxelLookup << std::endl;
                    logging::error << "CurrentMin: " << currentXYZMin << std::endl;
                    logging::error << "abc:" << abc << std::endl;
                    logging::error << "xyz:" << xyz << std::endl;
                    logging::error << "boxSize: " << boxSize << std::endl;
                    throw std::runtime_error( "Local Index out of Data Index block range" );
                }

                boost::int32_t index = dataIndices[blockIndex];

                if( index > -1 ) {
                    if( index > (int)distanceData.size() ) {

                        logging::error << "lookup: " << voxelLookup << std::endl;
                        logging::error << "CurrentMin: " << currentXYZMin << std::endl;
                        logging::error << "abc:" << abc << std::endl;
                        logging::error << "xyz:" << vector3( x, y, z ) << std::endl;
                        logging::error << "boxSize: " << boxSize << std::endl;

                        logging::error << "index: " << index << std::endl;
                        throw std::runtime_error( "weno3 error - bad index when collection x lookups" );
                    }

                    values[x + 1] = distanceData[index];

                } else {
                    throw std::runtime_error( "Undefined value in weno3_lookup() - " +
                                              boost::lexical_cast<std::string>( index ) );
                }
            }

            xlerps[resultIndex] = weno3_interpolant( values, weight.x );
        }
    }

    // compute the next set of 4 interpolations along y
    for( int z = 0; z < 4; ++z ) {
        ylerps[z] = weno3_interpolant( xlerps + z * 4, weight.y );

        // FF_LOG(debug) << "\t ylerp[" <<  z<< "]:" << ylerps[z] << std::endl;
    }

    // compute the final interpolation on z
    // save the resulting value
    return weno3_interpolant( ylerps, weight.z );
}

void create_staggered_smoothness_indicator_x_channel( frantic::fluids::rle_voxel_field& field,
                                                      const frantic::tstring& staggeredChannelName,
                                                      const frantic::tstring& indicatorChannelOneName,
                                                      const frantic::tstring& indicatorChannelTwoName ) {
    field.add_channel<float>( indicatorChannelOneName );
    field.add_channel<float>( indicatorChannelTwoName );

    const rle_index_spec& ris = field.get_rle_index_spec();

    const_rle_channel_accessor<vector3f> chanAcc = field.get_channel_accessor<vector3f>( staggeredChannelName );

    rle_channel_accessor<float> indicator1Acc = field.get_channel_accessor<float>( indicatorChannelOneName );
    rle_channel_accessor<float> indicator2Acc = field.get_channel_accessor<float>( indicatorChannelTwoName );

    const std::vector<boost::int32_t>& bcToRunIndex = ris.get_bc_to_run_index_vector();
    const std::vector<run_data>& runIndexData = ris.get_run_index_data_vector();

    // precompute the index to the center cell (2,2,2) in terms of the local 5x5x5 box
    // this is written out for clarity, but could be collapsed to save a few cyclces...
    size3 boxSize( 4, 4, 4 );
    vector3 centerOffset( 1 );

    // tbb::task_scheduler_init taskScheduleInit;

    int cellBoxIndex =
        centerOffset.x + centerOffset.y * boxSize.xsize() + centerOffset.z * boxSize.xsize() * boxSize.ysize();

    int xsize = boxSize.xsize();
    int ysize = boxSize.ysize();

    boost::int32_t baseIndex = 1 + xsize + xsize * ysize;

    boundbox3 bounds = ris.outer_bounds();
    vector3 boundsMin = bounds.minimum();
    vector3 boundsMax = bounds.maximum();

    for( int z = boundsMin.z; z <= boundsMax.z; ++z ) {
        int cIndex = ris.z_to_c( z );
        for( int y = boundsMin.y; y <= boundsMax.y; ++y ) {
            // compute the bcIndex of the scanline
            int bcIndex = ris.y_to_b( y ) + cIndex * bounds.ysize();

            // get the run range of the scanline
            int runRangeStart = bcToRunIndex[bcIndex];
            int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

            // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the runs
            // with the iterator
            if( runRangeStart != runRangeEnd || runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                // construct a block iterator for this scanline

                boundbox3 box( vector3( boundsMin.x - centerOffset.x, y - centerOffset.y, z - centerOffset.z ),
                               boxSize );

                for( levelset::rle_defined_box_iterator boxIter( ris, box );
                     !boxIter.is_xplane_finished( centerOffset.x ); ++boxIter ) {
                    const boost::int32_t* const dataIndices = boxIter.get_indices();

                    // const boundbox3& currentBox = boxIter.current_box();
                    // const vector3& currentMin = currentBox.minimum();

                    boost::int32_t cellIndex = dataIndices[cellBoxIndex];

                    if( cellIndex < 0 ) {
                        continue;
                    }

                    float f[4];

                    vector3f s1, s2;
                    for( int comp = 0; comp < 3; ++comp ) {
                        boost::int32_t offset = 0;

                        switch( comp ) {
                        case 1:
                            offset = xsize;
                            break;
                        case 2:
                            offset = xsize * ysize;
                            break;
                        }

                        for( int i = -1; i < 3; ++i ) {
                            boost::int32_t index = dataIndices[baseIndex + i * offset];

                            if( index < 0 ) {
                                f[i + 1] = 0.f;
                            } else {
                                f[i + 1] = chanAcc[index][comp];
                            }
                        }

                        s1[comp] = weno_indicator1( f[0], f[1], f[2] );
                        s2[comp] = weno_indicator2( f[1], f[2], f[3] );
                    }

                    // s1 = weno_indicator1( f[0], f[1], f[2] );
                    // s2 = weno_indicator2( f[1], f[2], f[3] );

                    indicator1Acc[cellIndex] = s1[0];
                    indicator2Acc[cellIndex] = s2[0];
                }
            }
        }
    }
}

} // namespace volumetrics
} // namespace frantic
