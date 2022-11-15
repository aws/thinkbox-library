// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/trimesh3_scan_conversion.hpp>
#include <frantic/volumetrics/levelset/geometry_to_levelset.hpp>

#include <frantic/diagnostics/profiling_section.hpp>

// Temporary include, for debugging.
#include <iomanip>

//#define _MESH_DUMP
//#define _PRINT_GTOLS_PROFILING

#ifdef _MESH_DUMP // include the mesh io code, if we are interested in saving out the mesh used in the conversion
#include <frantic/files/paths.hpp>
#include <frantic/geometry/trimesh3_file_io.hpp>
#include <frantic/volumetrics/levelset/rle_level_set_file_io.hpp>
#endif

using namespace std;
using namespace frantic::diagnostics;
using namespace frantic::channels;
using namespace frantic::geometry;
using namespace frantic::graphics;
using namespace frantic::volumetrics::levelset;

namespace frantic {
namespace volumetrics {
namespace levelset {

namespace detail {

/**
 * This is an internal function, and therefore in the detail namespace.
 *
 * Given a set of intersections (as provided by the trimesh3_scan_convert function), this converts them
 * into an rle_level_set.  A slightly confusing thing about this function is how the input and output
 * axes are a permutation of each other.  The intersections depths were created with an XY plane and contain
 * Z depths of intersections.  The output level set will have its compression axis matching the direction
 * along which the intersections were found, thus the Z of the intersectionDepths will map to the X of the
 * outResult level set.
 *
 * This axis permutation from the input intersections to the output level set is (X,Y,Z) -> (Z,X,Y).
 *
 * The outResult should already have the appropriate voxel space configured, and the voxel bounds provided
 * are in terms of the destination voxel space, not the source space where intersections are taken.  The
 * source
 *
 * @param  intersectionDepths  This is an array of intersections with the geometry, along the Z direction in the XY
 *                             plane of the input intersection coordinate system.
 * @param  geometry            This is the geometry being converted to a level set.
 * @param  interfaceVoxelWidthInside  This is the number of voxel units away from the interface inside the geometry that
 * should be defined.
 * @param  interfaceVoxelWidthOutside  This is the number of voxel units away from the interface outside the geometry
 * that should be defined.
 * @param  isHalfOpenSurface   If the input geometry was determined to be a surface which doesn't close a volume, for
 * instance a flat plane, this flag will be set to true.
 * @param  exteriorRegionCodeNeg  This should be -1 if the exterior is "outside", -2 if the exterior is "inside".  If
 * isHalfOpenSurface is false, this value is used to determine how a non-intersection is treated in the vote counting.
 * If isHalfOpenSurface is true, this is how non-intersection is treated in the negative direction of vote counting.  In
 *                             that case it should be 0 if the exterior is ambiguous.
 * @param  exteriorRegionCodePos  This should be -1 if the exterior is "outside", -2 if the exterior is "inside".  If
 *                             isHalfOpenSurface is true, this is how non-intersection is treated in the positive
 * direction of vote counting. It should be 0 if the exterior is ambiguous.  If isHalfOpenSurface is false, this must be
 *                             equal to exteriorRegionCodeNeg.
 * @param  resultVoxelBounds   This the voxel bounding box to fill for the outResult level set.
 * @param  outResult           The level set object to fill.
 *
 */
void convert_intersections_to_level_set(
    const std::vector<std::vector<geometry::scan_conversion_intersection>>& intersectionDepths,
    const geometry::trimesh3& geometry, float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside,
    bool isHalfOpenSurface, int exteriorRegionCodeNeg, int exteriorRegionCodePos, const boundbox3& resultVoxelBounds,
    rle_level_set& outResult ) {
    // Get the source 2D bounding box of the intersections, and ensure that the intersections array has the correct
    // size.
    frantic::graphics2d::vector2 sourceOrigin( resultVoxelBounds.minimum().y, resultVoxelBounds.minimum().z );
    frantic::graphics2d::size2 sourceSize( resultVoxelBounds.maximum().y - sourceOrigin.x + 1,
                                           resultVoxelBounds.maximum().z - sourceOrigin.y + 1 );
    if( (int)intersectionDepths.size() != sourceSize.get_area() )
        throw runtime_error(
            "convert_intersections_to_level_set: The input intersection depths array had the incorrect size." );

    outResult.clear();
    voxel_coord_system& vcs = outResult.m_voxelCoordSystem;
    rle_index_spec& ris = outResult.m_rleIndex;
    ris.m_abcCoordOrigin = resultVoxelBounds.minimum();
    ris.m_abcCoordSize = resultVoxelBounds.size();

    float voxelLength = vcs.voxel_length();

    // Save the untouched voxel widths into the level set
    outResult.m_interfaceVoxelWidthInside = interfaceVoxelWidthInside;
    outResult.m_interfaceVoxelWidthOutside = interfaceVoxelWidthOutside;

    // Expand the interface voxel widths by sqrt(3) so that we're guaranteed to get the voxels we need
    interfaceVoxelWidthInside *= 1.73f;
    interfaceVoxelWidthOutside *= 1.73f;

    float interfaceWorldWidthInside = interfaceVoxelWidthInside * voxelLength;
    float interfaceWorldWidthOutside = interfaceVoxelWidthOutside * voxelLength;
    float largeDistanceValue = max( interfaceWorldWidthInside, interfaceWorldWidthOutside ) + voxelLength;

    // Reset the recorded data size to 0, so when we create the named channels they start with a size of zero
    ris.m_dataSize = 0;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );
    ris.m_runData.clear();

    // Get the vertex channel names in the mesh
    vector<frantic::tstring> channelNames;
    geometry.get_vertex_channel_names( channelNames );

    // Accessors for the named channels from the mesh
    vector<const_trimesh3_vertex_channel_general_accessor> inputChannelAccessors;
    vector<rle_channel_general_accessor> outputChannelAccessors;
    int largestPrimitiveSize = 1;
    // Build the corresponding arrays of input and output channels.
    for( unsigned channelIndex = 0; channelIndex < channelNames.size(); ++channelIndex ) {
        // A small hard-coded list of channels we don't want to copy into the level set.  For example, we
        // don't want a channel of normals.
        // Also, only support propagation of floating point data types.  Unlike conversion achieved by looking up the
        // nearest point on the surface, the method of blending values doesn't work for integer types here.
        const_trimesh3_vertex_channel_general_accessor inputAccessor =
            geometry.get_vertex_channel_general_accessor( channelNames[channelIndex] );
        if( channels::is_channel_data_type_float( inputAccessor.data_type() ) &&
            channelNames[channelIndex] != _T("Normal") ) {
            inputChannelAccessors.push_back( inputAccessor );
            outResult.add_channel( channelNames[channelIndex], inputAccessor.arity(), inputAccessor.data_type() );
            outputChannelAccessors.push_back( outResult.get_channel_general_accessor( channelNames[channelIndex] ) );
            largestPrimitiveSize = max( (int)inputChannelAccessors.back().primitive_size(), largestPrimitiveSize );
        }
    }
    // Create a temp buffer to hold one intermediate named channel value.
    vector<char> tempNamedChannelValue( largestPrimitiveSize );

    // Auxiliary accessors specific to this conversion algorithm
    rle_channel_accessor<float> channelWeightAccessor, minDistanceInsideAccessor, minDistanceOutsideAccessor;
    rle_channel_accessor<int> voteCountAccessor;

    // Build the special-purpose channels and get their accessors
    outResult.add_channel<float>( _T("G2LS_TEMP_ChannelWeight") );
    outResult.add_channel<float>( _T("G2LS_TEMP_MinDistanceInside") );
    outResult.add_channel<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    outResult.add_channel<int>( _T("G2LS_TEMP_VoteCount") );
    channelWeightAccessor = outResult.get_channel_accessor<float>( _T("G2LS_TEMP_ChannelWeight") );
    minDistanceInsideAccessor = outResult.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    minDistanceOutsideAccessor = outResult.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    voteCountAccessor = outResult.get_channel_accessor<int>( _T("G2LS_TEMP_VoteCount") );

    // The exterior region code either gets no votes if the surface is half-open, or gets 2 votes to
    // the appropriate vote count. Outside votes get into the least significant 4 bits, and inside
    // votes get into the next 4 bits.
    if( isHalfOpenSurface )
        // Non-intersections outside the bounding box get no votes in this case
        ris.m_exteriorRegionCode = -0x100;
    else if( exteriorRegionCodeNeg == -1 ) {
        // Count 2 votes for "outside"
        ris.m_exteriorRegionCode = -0x102;
    } else {
        // Count 2 votes for "inside"
        ris.m_exteriorRegionCode = -0x120;
    }

    /*
    int debugB = 70; //ris.m_abcCoordSize.ysize()/2;
    int debugC = 99; //ris.m_abcCoordSize.zsize()/2;
    */

    float voxelXMin = ris.m_abcCoordOrigin.x + 0.5f,
          voxelXMax = ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() - 0.5f;

    vector3 xyz;
    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            xyz.y = b + ris.m_abcCoordOrigin.y;
            xyz.z = c + ris.m_abcCoordOrigin.z;

            int bcIndex = b + c * ris.m_abcCoordSize.ysize();
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            const vector<scan_conversion_intersection>& isect = intersectionDepths[bcIndex];

            /*
            if( b == debugB && c == debugC ) {
              cout << "DEBUG COORD" << endl;
              for( unsigned i = 0; i < isect.size(); ++i )
                cout << isect[i].z << " ";
              cout << endl;
              cout << endl;
            }
            */

            if( !isect.empty() ) {
                // Start of this variable indicating "outside" or "inside" depending on the exterior region code
                bool prevNormalFacingZPositive = exteriorRegionCodeNeg == -1;
                float prevVoxelX =
                    ris.m_abcCoordOrigin.x - 2 * ( max( interfaceVoxelWidthInside, interfaceVoxelWidthOutside ) + 1 );
                // Find the first intersection within the range we are concerned with
                float voxelX = isect[0].z / voxelLength;
                bool normalFacingZPositive = isect[0].normalFacingZPositive;

                // The default state is an undefined region of the exterior region code
                bool creatingUndefinedRun = true;
                int creatingRegionCode = ris.m_exteriorRegionCode;

                // Go through all the intersections, and build up the runs.
                for( unsigned isectIndex = 0; isectIndex <= isect.size(); ++isectIndex ) {
                    if( isectIndex < isect.size() ) {
                        normalFacingZPositive = isect[isectIndex].normalFacingZPositive;
                        voxelX = isect[isectIndex].z / voxelLength;
                    } else {
                        normalFacingZPositive = exteriorRegionCodePos != -1;
                        voxelX = ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() +
                                 2 * ( max( interfaceVoxelWidthInside, interfaceVoxelWidthOutside ) + 1 );
                    }

                    // If the interval from prevVoxelX to voxelX intersects with the voxel x extents of the runs, then
                    // we need to incorporate this
                    if( voxelXMin < voxelX && prevVoxelX < voxelXMax ) {
                        // Compute the extents of the runs
                        // First run, defined
                        int firstXStart = (int)floorf( prevVoxelX - 0.5f ) + 1;
                        int firstXEnd = (int)floorf(
                            prevVoxelX - 0.5f +
                            ( prevNormalFacingZPositive ? interfaceVoxelWidthOutside : interfaceVoxelWidthInside ) );
                        // Third run, defined
                        int thirdXStart = (int)floorf( voxelX - 0.5f -
                                                       ( normalFacingZPositive ? interfaceVoxelWidthInside
                                                                               : interfaceVoxelWidthOutside ) ) +
                                          1;
                        int thirdXEnd = (int)floorf( voxelX - 0.5f );
                        // Second run, undefined
                        int secondXStart = firstXEnd + 1, secondXEnd = thirdXStart - 1;

                        // If the first and third intervals overlap, need to adjust them so they don't
                        if( firstXEnd >= thirdXStart ) {
                            // The interval center is where the two propagating fronts would meet.
                            float centerX = 0.5f * ( prevVoxelX + voxelX );
                            firstXEnd = (int)floorf( centerX - 0.5f );
                            thirdXStart = (int)ceilf( centerX - 0.5f );
                            if( firstXEnd == thirdXStart )
                                ++thirdXStart;
                        }

                        // Make sure the intervals are contained within the voxel x extents of the runs
                        if( firstXStart < ris.m_abcCoordOrigin.x )
                            firstXStart = ris.m_abcCoordOrigin.x;
                        if( secondXStart < ris.m_abcCoordOrigin.x )
                            secondXStart = ris.m_abcCoordOrigin.x;
                        if( thirdXStart < ris.m_abcCoordOrigin.x )
                            thirdXStart = ris.m_abcCoordOrigin.x;
                        if( firstXEnd >= ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() )
                            firstXEnd = ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() - 1;
                        if( secondXEnd >= ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() )
                            secondXEnd = ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() - 1;
                        if( thirdXEnd >= ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() )
                            thirdXEnd = ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() - 1;

                        // cout << i << "/" << isect.size() << ", first: " << firstXStart << "-" << firstXEnd << ",
                        // third: " << thirdXStart << "-" << thirdXEnd << endl;

                        // Figure out the vote counts for this interval.  This deals specially with the boundary cases,
                        // which must be skipped selectively when isHalfOpenSurface is true and exteriorRegionCodeNeg
                        // and/or exteriorRegionCodePos is 0.
                        int intervalVoting = 0;
                        if( isectIndex > 0 ) {
                            if( prevNormalFacingZPositive )
                                intervalVoting += 0x001; // Cast a vote for "outside"
                            else
                                intervalVoting += 0x010; // Cast a vote for "inside"
                        } else {
                            if( !isHalfOpenSurface ) {
                                if( prevNormalFacingZPositive )
                                    intervalVoting += 0x001; // Cast a vote for "outside"
                                else
                                    intervalVoting += 0x010; // Cast a vote for "inside"
                            } else if( exteriorRegionCodeNeg != 0 ) {
                                // In the case of half-open surfaces with an exterior vote, triplicate it.  This is
                                // done because we were getting Flood:Surf surfaces with bad loops happening at extreme
                                // wave peaks.  Only doubling the vote didn't help enough.
                                if( prevNormalFacingZPositive )
                                    intervalVoting += 0x003; // Cast three votes for "outside"
                                else
                                    intervalVoting += 0x030; // Cast three votes for "inside"
                            }
                        }

                        if( isectIndex < isect.size() ) {
                            if( normalFacingZPositive )
                                intervalVoting += 0x010; // Cast a vote for "inside"
                            else
                                intervalVoting += 0x001; // Cast a vote for "outside"
                        } else {
                            if( !isHalfOpenSurface ) {
                                if( normalFacingZPositive )
                                    intervalVoting += 0x010; // Cast a vote for "inside"
                                else
                                    intervalVoting += 0x001; // Cast a vote for "outside"
                            } else if( exteriorRegionCodePos != 0 ) {
                                // In the case of half-open surfaces with an exterior vote, triplicate it.  This is
                                // done because we were getting Flood:Surf surfaces with bad loops happening at extreme
                                // wave peaks.  Only doubling the vote didn't help enough.
                                if( normalFacingZPositive )
                                    intervalVoting += 0x030; // Cast three votes for "inside"
                                else
                                    intervalVoting += 0x003; // Cast three votes for "outside"
                            }
                        }

                        /*
                        if( b == debugB && c == debugC ) {
                          cout << "first interval: " << firstXStart << ", " << firstXEnd << endl;
                          cout << "second interval: " << secondXStart << ", " << secondXEnd << endl;
                          cout << "third interval: " << thirdXStart << ", " << thirdXEnd << endl;
                        }
                        */

                        // Now add the three intervals into the rle runs
                        if( firstXStart <= firstXEnd ) {
                            const geometry::scan_conversion_intersection& isectPrevRef = isect[isectIndex - 1];
                            if( creatingUndefinedRun ) {
                                ris.m_runData.push_back(
                                    run_data( firstXStart, (int)outResult.m_distanceData.size() ) );
                                creatingUndefinedRun = false;
                            }

                            for( int x = firstXStart; x <= firstXEnd; ++x ) {
                                float voxelDistance = ( x + 0.5f - prevVoxelX );
                                float distance = voxelDistance * voxelLength;

                                // Make the channel data weight based on the voxel distance, so it's related to our
                                // sampling. Also adjust it so that its maximum value is 8.0.  this is so that
                                // interpolating half float channels avoids overflowing their numeric bounds.
                                float channelWeight = 0.08f / ( voxelDistance < 0.01f ? 0.01f : voxelDistance );
                                channelWeightAccessor.add_element( channelWeight );
                                // Add the votes for this run
                                voteCountAccessor.add_element( intervalVoting );

                                if( prevNormalFacingZPositive ) {
                                    outResult.m_distanceData.push_back( distance );
                                    minDistanceInsideAccessor.add_element( largeDistanceValue );
                                    minDistanceOutsideAccessor.add_element( distance );
                                } else {
                                    outResult.m_distanceData.push_back( -distance );
                                    minDistanceInsideAccessor.add_element( distance );
                                    minDistanceOutsideAccessor.add_element( largeDistanceValue );
                                }

                                // Add the weight times the vertex channel value which is interpolated using the
                                // barycentric coordinates
                                char* data = &tempNamedChannelValue[0];
                                for( size_t channelIndex = 0, channelIndexEnd = inputChannelAccessors.size();
                                     channelIndex != channelIndexEnd; ++channelIndex ) {
                                    const_trimesh3_vertex_channel_general_accessor& inputCA =
                                        inputChannelAccessors[channelIndex];
                                    rle_channel_general_accessor& outputCA = outputChannelAccessors[channelIndex];
                                    inputCA.get_barycentric( isectPrevRef.faceIndex, isectPrevRef.barycentricCoord,
                                                             data );
                                    outputCA.m_weightedSumFunction( &channelWeight, &data, 1, inputCA.arity(),
                                                                    outputCA.add_element() );
                                }
                            }
                        }

                        if( secondXStart <= secondXEnd ) {
                            int regionCode = -0x100 - intervalVoting;
                            if( !creatingUndefinedRun ) {
                                ris.m_runData.push_back( run_data( secondXStart, regionCode ) );
                                creatingUndefinedRun = true;
                                creatingRegionCode = regionCode;
                            } else if( regionCode != creatingRegionCode ) {
                                ris.m_runData.push_back( run_data( secondXStart, regionCode ) );
                                creatingRegionCode = regionCode;
                            }
                        }

                        if( thirdXStart <= thirdXEnd ) {
                            const geometry::scan_conversion_intersection& isectRef = isect[isectIndex];
                            if( creatingUndefinedRun ) {
                                ris.m_runData.push_back(
                                    run_data( thirdXStart, (int)outResult.m_distanceData.size() ) );
                                creatingUndefinedRun = false;
                            }

                            for( int x = thirdXStart; x <= thirdXEnd; ++x ) {
                                float voxelDistance = ( voxelX - ( x + 0.5f ) );
                                float distance = voxelDistance * voxelLength;

                                // Make the channel data weight based on the voxel distance, so it's related to our
                                // sampling. Also adjust it so that its maximum value is 8.0.  this is so that
                                // interpolating half float channels avoids overflowing their numeric bounds.
                                float channelWeight = 0.08f / ( voxelDistance < 0.01f ? 0.01f : voxelDistance );
                                channelWeightAccessor.add_element( channelWeight );
                                // Add the votes for this run
                                voteCountAccessor.add_element( intervalVoting );

                                if( normalFacingZPositive ) {
                                    outResult.m_distanceData.push_back( -distance );
                                    minDistanceInsideAccessor.add_element( distance );
                                    minDistanceOutsideAccessor.add_element( largeDistanceValue );
                                } else {
                                    outResult.m_distanceData.push_back( distance );
                                    minDistanceInsideAccessor.add_element( largeDistanceValue );
                                    minDistanceOutsideAccessor.add_element( distance );
                                }

                                // Add the weight times the vertex channel value which is interpolated using the
                                // barycentric coordinates
                                char* data = &tempNamedChannelValue[0];
                                for( size_t channelIndex = 0, channelIndexEnd = inputChannelAccessors.size();
                                     channelIndex != channelIndexEnd; ++channelIndex ) {
                                    const_trimesh3_vertex_channel_general_accessor& inputCA =
                                        inputChannelAccessors[channelIndex];
                                    rle_channel_general_accessor& outputCA = outputChannelAccessors[channelIndex];
                                    inputCA.get_barycentric( isectRef.faceIndex, isectRef.barycentricCoord, data );
                                    outputCA.m_weightedSumFunction( &channelWeight, &data, 1, inputCA.arity(),
                                                                    outputCA.add_element() );
                                }
                            }
                        }
                    }

                    prevNormalFacingZPositive = normalFacingZPositive;
                    prevVoxelX = voxelX;
                }

                // If any runs were flagged, finish off the last run.  Otherwise add a zero-sized run
                if( (int)ris.m_runData.size() > ris.m_bcToRunIndex[bcIndex] ) {
                    // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to
                    // optimize the storage and access.
                    if( !( creatingUndefinedRun && creatingRegionCode == ris.m_exteriorRegionCode ) )
                        ris.m_runData.push_back( run_data( ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize(), -1 ) );
                } else {
                    ris.m_runData.push_back( run_data( 0, -1 ) );
                    ris.m_runData.push_back( run_data( 0, -1 ) );
                }
            } else {
                // Create an empty run for this BC coordinate
                ris.m_runData.push_back( run_data( 0, -1 ) );
                ris.m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = outResult.m_distanceData.size();

    // run a consistency check
    stringstream sout;
    if( !ris.check_consistency( sout ) ) {
        ofstream fout( "c:\\temp\\badness.txt" );
        ris.dump( fout );
        throw std::runtime_error( "Intersections To Level Set: RLE Index Spec consistency check failed:\n" +
                                  sout.str() );
    }
}

/**
 * This is an internal function, and therefore in the detail namespace.
 *
 * The function removes pairs of intersections that are "inside" the surface.  If you
 * consider a ray passing through a mesh, and add +1 each time a surface faces you,
 * or -1 each time it faces away from you, you will get a count that can tell you whether
 * you are inside or outside the whole object, assuming the shape is well-structured.  By
 * removing pairs of intersections that cause the count to go above +1 or below -1, the
 * code which converts these intersections into a primitive level set doesn't have to worry
 * about such bookkeeping itself.
 *
 * @param  intersectionDepths           The arrays of intersection depths.
 *
 */
void remove_internal_intersections( vector<vector<geometry::scan_conversion_intersection>>& intersectionDepths ) {
    vector<vector<geometry::scan_conversion_intersection>>::iterator i = intersectionDepths.begin(),
                                                                     ie = intersectionDepths.end();
    vector<char> flags( 50 );
    for( ; i != ie; ++i ) {
        // Internal pairs can only be removed if there are 3 or more intersections
        if( i->size() > 2 ) {
            // Set the flags array to the size we need, and initialize it to all 0
            flags.resize( i->size() );
            memset( &flags[0], 0, flags.size() );

            // Sweep through the intersections in the forward direction
            // Skip the last index, i->size()-1, so that the boundary intersections are retained.
            int counter = 0;
            for( size_t j = 0, je = i->size() - 1; j != je; ++j ) {
                // If the counter is outside the [-1,1] interval, then mark the flag
                // Note that we have to do this both before and after adjusting the counter
                if( counter < -1 || counter > 1 )
                    flags[j] = 1;
                if( ( *i )[j].normalFacingZPositive )
                    ++counter;
                else
                    --counter;
                // If this intersection caused the counter to go outside the [-1,1] interval, then mark the flag
                if( counter < -1 || counter > 1 )
                    flags[j] = 1;
            }

            // Sweep through the intersections in the backward direction
            // Skip the last index, 0, so that the boundary intersections are retained.
            counter = 0;
            for( size_t j = i->size() - 1; j != 0; --j ) {
                // If the counter is outside the [-1,1] interval, then mark the flag
                // Note that we have to do this both before and after adjusting the counter
                if( counter < -1 || counter > 1 )
                    flags[j] = 1;
                if( ( *i )[j].normalFacingZPositive )
                    ++counter;
                else
                    --counter;
                // If this intersection caused the counter to go outside the [-1,1] interval, then mark the flag
                if( counter < -1 || counter > 1 )
                    flags[j] = 1;
            }

            // Sweep through the intersections, deleting the ones which have been flagged for deletion
            size_t destIndex = 0;
            for( size_t j = 0, je = i->size(); j != je; ++j ) {
                if( !flags[j] ) {
                    if( destIndex != j )
                        ( *i )[destIndex] = ( *i )[j];
                    ++destIndex;
                }
            }
            // Shrink the array if necessary
            if( destIndex != i->size() )
                i->resize( destIndex );
        }
    }
}

/**
 * This is an internal function, and therefore in the detail namespace.
 *
 * The function dilates the rle_index_spec region of the level set, so that voxels that have no
 * intersections in the X or Y direction but are near the surface will become defined.
 *
 * One difficulty which slightly complicates the implementation of this function is the possibility of multiple
 * undefined runs in a row, as well as undefined runs before the first defined run and after the last defined run.  The
 * dilation amount away from the end of a defined run could span several undefined runs as well as beyond the range of
 * runs into the exterior.  When doing the dilation, this function needs to convert the votes in the multiple undefined
 * runs into the vote counting channel.
 *
 * An important assumption that is made in this function is that all distances that have been computed
 * up to this point are orthogonal to the X direction.  With this assumption, Pythagoras' theorem can be
 * applied to figure out how far the dilation should go.
 *
 * @param  levelSet           The level set object to dilate.
 *
 */
void mesh_to_level_set_region_x_dilation( rle_level_set& levelSet ) {
    //	profiling_section psTotal("Total"), psInit("Initialization"), psPerNonEmptyScanline("Per Non-Empty Scanline");
    //	profiling_section psPerScanlineInit("Per-Scanline Init"), psPerDefinedRunInit("Per-Defined Run Init");
    //	profiling_section psForwardDilate("Forward Dilation"), psBackwardDilate("Backward Dilation");
    //	profiling_section psUndefinedCopy("Undefined Copy"), psDefinedCopy("Defined Copy");

    //	psTotal.enter();
    //	psInit.enter();

    float voxelLength = levelSet.get_voxel_coord_system().voxel_length();
    float interfaceWorldWidthInside = levelSet.m_interfaceVoxelWidthInside * voxelLength;
    float interfaceWorldWidthOutside = levelSet.m_interfaceVoxelWidthOutside * voxelLength;

    // Swap out to a temp variable, so we copy the result into the levelSet variable.
    rle_level_set tempLevelSet( levelSet.get_voxel_coord_system() );
    tempLevelSet.swap( levelSet );
    levelSet.m_interfaceVoxelWidthInside = tempLevelSet.m_interfaceVoxelWidthInside;
    levelSet.m_interfaceVoxelWidthOutside = tempLevelSet.m_interfaceVoxelWidthOutside;
    levelSet.m_insideDistance = tempLevelSet.m_insideDistance;
    levelSet.m_outsideDistance = tempLevelSet.m_outsideDistance;

    rle_index_spec& ris = levelSet.m_rleIndex;
    const rle_index_spec& tempRis = tempLevelSet.m_rleIndex;

    int exteriorVoteCount = ( -tempRis.m_exteriorRegionCode ) & 0xff;

    ris.m_exteriorRegionCode = tempRis.m_exteriorRegionCode;
    ris.m_abcCoordOrigin = tempRis.m_abcCoordOrigin;
    ris.m_abcCoordSize = tempRis.m_abcCoordSize;
    ris.m_dataSize = 0;

    // clear any data already stored in ris
    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );

    // Set up the named channel accessors
    levelSet.m_namedChannels.clear();

    vector<rle_channel_general_accessor> outputChannels;
    vector<const_rle_channel_general_accessor> inputChannels;
    for( map<frantic::tstring, rle_channel>::const_iterator i = tempLevelSet.m_namedChannels.begin();
         i != tempLevelSet.m_namedChannels.end(); ++i ) {
        // Add all the channels from the input to the output
        levelSet.add_channel( i->first, i->second.arity(), i->second.data_type() );
        // All channels except for three should just have their values duplicated
        if( i->first != _T("G2LS_TEMP_MinDistanceInside") && i->first != _T("G2LS_TEMP_MinDistanceOutside") &&
            i->first != _T("G2LS_TEMP_VoteCount") ) {
            outputChannels.push_back( levelSet.get_channel_general_accessor( i->second.name() ) );
            inputChannels.push_back( tempLevelSet.get_channel_general_accessor( i->second.name() ) );
        }
    }

    // Auxiliary accessors specific to this conversion algorithm
    const_rle_channel_accessor<float> inputMinDistanceInsideAccessor, inputMinDistanceOutsideAccessor;
    const_rle_channel_accessor<int> inputVoteCountAccessor;
    rle_channel_accessor<float> outputMinDistanceInsideAccessor, outputMinDistanceOutsideAccessor;
    rle_channel_accessor<int> outputVoteCountAccessor;
    inputMinDistanceInsideAccessor = tempLevelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    inputMinDistanceOutsideAccessor = tempLevelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    inputVoteCountAccessor = tempLevelSet.get_channel_accessor<int>( _T("G2LS_TEMP_VoteCount") );
    outputMinDistanceInsideAccessor = levelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    outputMinDistanceOutsideAccessor = levelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    outputVoteCountAccessor = levelSet.get_channel_accessor<int>( _T("G2LS_TEMP_VoteCount") );

    //	psInit.exit();

    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            int bcIndex = b + c * tempRis.m_abcCoordSize.ysize();
            int runRangeStart = tempRis.m_bcToRunIndex[bcIndex], runRangeEnd = tempRis.m_bcToRunIndex[bcIndex + 1] - 2;

            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            /*
            if( b >= 55 && b <= 56 && c == 5 ) {
              cout << "Processing BC " << b << ", " << c << endl;
              for( int run = runRangeStart; run <= runRangeEnd; ++run ) {
                cout << "Run " << run << ": " << tempRis.m_runData[run].x << " to " << tempRis.m_runData[run+1].x-1 <<
            ", data index " << tempRis.m_runData[run].dataIndex << endl;
              }
            }
            */

            if( runRangeStart != runRangeEnd ||
                tempRis.m_runData[runRangeStart].x != tempRis.m_runData[runRangeStart + 1].x ) {
                //				psPerNonEmptyScanline.enter();

                //				psPerScanlineInit.enter();

                int run;
                bool creatingUndefinedRun = true;

                int prevDefinedRun = runRangeStart - 1, curDefinedRun = runRangeStart;
                // Increment through the runs to find the next defined one
                while( curDefinedRun <= runRangeEnd && tempRis.m_runData[curDefinedRun].dataIndex < 0 )
                    ++curDefinedRun;

                float forwardMinDistanceInside = 0, forwardMinDistanceOutside = 0;
                float backwardMinDistanceInside = 0, backwardMinDistanceOutside = 0;
                int forwardDataIndex, backwardDataIndex;
                // forwardX is the X value at the last defined voxel in the previous defined run
                // backwardX is the X value at the first defined voxel in the current defined run
                // forwardDilateX is the largest X value contained within the forward dilation region
                // backwardDilateX is the smallest X value contained within the backward dilation region
                boost::int32_t forwardX, forwardDilateX, backwardX, backwardDilateX;

                int scanlineEndX = tempRis.m_runData[runRangeEnd + 1].x;

                //				psPerScanlineInit.exit();

                //				cout << endl;

                // Keep looping until the previous defined run passes the end.  This way we will iterate once
                // for every gap between two defined runs as well as before the first defined run and after
                // the last defined run.
                while( prevDefinedRun <= runRangeEnd ) {

                    //					psPerDefinedRunInit.enter();

                    int forwardVoxelsToDilate = 0, backwardVoxelsToDilate = 0;

                    // Get the data index from the end of the previous defined run
                    if( prevDefinedRun >= runRangeStart ) {
                        forwardX = tempRis.m_runData[prevDefinedRun + 1].x - 1;
                        forwardDataIndex = tempRis.m_runData[prevDefinedRun].dataIndex + forwardX -
                                           tempRis.m_runData[prevDefinedRun].x;

                        forwardMinDistanceInside = inputMinDistanceInsideAccessor[forwardDataIndex];
                        forwardMinDistanceOutside = inputMinDistanceOutsideAccessor[forwardDataIndex];

                        // Figure out how much to dilate forward
                        float distanceToDilate = 0;
                        if( forwardMinDistanceInside + voxelLength < interfaceWorldWidthInside )
                            distanceToDilate = sqrtf( interfaceWorldWidthInside * interfaceWorldWidthInside -
                                                      forwardMinDistanceInside * forwardMinDistanceInside );
                        if( forwardMinDistanceOutside + voxelLength < interfaceWorldWidthOutside )
                            distanceToDilate =
                                max( distanceToDilate, sqrtf( interfaceWorldWidthOutside * interfaceWorldWidthOutside -
                                                              forwardMinDistanceOutside * forwardMinDistanceOutside ) );
                        forwardVoxelsToDilate =
                            distanceToDilate <= 0 ? 0 : (int)floorf( distanceToDilate / voxelLength );

                        //						cout << endl;
                        //						cout << "interfaceWorldWidthInside: " <<
                        //interfaceWorldWidthInside << endl; 						cout << "interfaceWorldWidthOutside: " <<
                        //interfaceWorldWidthOutside << endl; 						cout << "forwardMinDistanceInside: " <<
                        //forwardMinDistanceInside << endl; 						cout << "forwardMinDistanceOutside: " <<
                        //forwardMinDistanceOutside << endl; 						cout << "voxelLength: " << voxelLength << endl; 						cout <<
                        //"forward distanceToDilate: " << distanceToDilate << " / " << forwardVoxelsToDilate << endl;

                        forwardDilateX = forwardX + forwardVoxelsToDilate;
                    } else {
                        // If the current defined run stepped passed the end, then use the bounding box to determine
                        // these values. This ensures proper clipping of the dilation later.
                        forwardX = tempRis.m_abcCoordOrigin.x - 1;
                        forwardDilateX = forwardX;
                        forwardDataIndex = -1;
                    }

                    // Get the data index from the start of the current defined run
                    if( curDefinedRun <= runRangeEnd ) {
                        backwardX = tempRis.m_runData[curDefinedRun].x;
                        backwardDataIndex = tempRis.m_runData[curDefinedRun].dataIndex;

                        backwardMinDistanceInside = inputMinDistanceInsideAccessor[backwardDataIndex];
                        backwardMinDistanceOutside = inputMinDistanceOutsideAccessor[backwardDataIndex];

                        // Figure out how much to dilate backward
                        float distanceToDilate = 0;
                        if( backwardMinDistanceInside + voxelLength < interfaceWorldWidthInside )
                            distanceToDilate = sqrtf( interfaceWorldWidthInside * interfaceWorldWidthInside -
                                                      backwardMinDistanceInside * backwardMinDistanceInside );
                        if( backwardMinDistanceOutside + voxelLength < interfaceWorldWidthOutside )
                            distanceToDilate = max( distanceToDilate,
                                                    sqrtf( interfaceWorldWidthOutside * interfaceWorldWidthOutside -
                                                           backwardMinDistanceOutside * backwardMinDistanceOutside ) );
                        backwardVoxelsToDilate =
                            distanceToDilate <= 0 ? 0 : (int)floorf( distanceToDilate / voxelLength );

                        //						cout << "backward distanceToDilate: " << distanceToDilate
                        //<< " / " << backwardVoxelsToDilate
                        //<< endl;

                        backwardDilateX = backwardX - backwardVoxelsToDilate;
                    } else {
                        // If the current defined run stepped passed the end, then use the bounding box to determine
                        // these values. This ensures proper clipping of the dilation later.
                        backwardX = tempRis.m_abcCoordOrigin.x + tempRis.m_abcCoordSize.xsize();
                        backwardDilateX = backwardX;
                        backwardDataIndex = -1;
                    }

                    // Clip the dilation extents so that they stop at the next defined run or at the bounds.
                    if( forwardDilateX >= backwardX )
                        forwardDilateX = backwardX - 1;
                    if( backwardDilateX <= forwardX )
                        backwardDilateX = forwardX + 1;

                    // If the dilation regions overlap, then we have to clip them so they don't.  Do this by
                    // taking the center position between the two and splitting them there.
                    if( forwardDilateX >= backwardDilateX ) {
                        forwardDilateX = forwardDilateX + backwardDilateX;
                        // Ensure that this gets truncates towards negative (integer divide truncates towards 0)
                        if( forwardDilateX >= 0 )
                            forwardDilateX /= 2;
                        else
                            forwardDilateX = ( forwardDilateX - 1 ) / 2;

                        backwardDilateX = forwardDilateX + 1;
                    }

                    //					psPerDefinedRunInit.exit();

                    /*
                    if( b >= 55 && b <= 56 && c == 5 ) {
                      cout << "maxX: " << tempRis.m_abcCoordOrigin.x + tempRis.m_abcCoordSize.xsize() - 1 << endl;
                      cout << "backwardDilateX: " << backwardDilateX << ", backwardX: " << backwardX << endl;
                      cout << "forwardX: " << forwardX << ", forwardDilateX: " << forwardDilateX << endl;
                    }
                    */

                    //////////
                    // First do the forward dilation
                    //////////
                    //					psForwardDilate.enter();

                    if( forwardDilateX > forwardX ) {
                        //						cout << "Forward dilation of " << forwardDilateX -
                        //forwardX << " voxels" << endl;
                        // Don't need to check if we're creating a defined run, because it's guaranteed that we are
                        // currently.
                        float lsDistance =
                            tempLevelSet
                                .m_distanceData[forwardDataIndex]; // Keeping a simple distance function in play too

                        int voteCount;
                        run = prevDefinedRun + 1;
                        if( run <= runRangeEnd )
                            voteCount = ( -tempRis.m_runData[run].dataIndex ) & 0xff;
                        else {
                            voteCount = exteriorVoteCount;
                            // In this case, the forward dilation has extended past the runs, so the terminator at the
                            // end needs to go further
                            if( forwardDilateX >= scanlineEndX )
                                scanlineEndX = forwardDilateX + 1;
                        }

                        // Go through all the undefined runs within the dilation region, and copy the votes
                        int firstForwardX = forwardX;
                        ++forwardX;
                        ++run;
                        while( forwardX <= forwardDilateX ) {
                            // runDilateX is the amount we dilate within the current undefined run. (There might be
                            // multiple undefined runs in the gap between defined runs)
                            int runDilateX;
                            if( run - 1 <= runRangeEnd )
                                runDilateX = min( tempRis.m_runData[run].x - 1, forwardDilateX );
                            else
                                runDilateX = forwardDilateX;

                            // Go through all the voxels up to the start of this run
                            while( forwardX <= runDilateX ) {
                                //								cout << "one forward dilation voxel"
                                //<< endl;
                                outputVoteCountAccessor.add_element( voteCount );

                                ++forwardX;
                            }

                            if( run <= runRangeEnd )
                                voteCount = ( -tempRis.m_runData[run].dataIndex ) & 0xff;
                            else {
                                voteCount = exteriorVoteCount;
                                // In this case, the forward dilation has extended past the runs, so the terminator at
                                // the end needs to go further
                                if( forwardDilateX >= scanlineEndX )
                                    scanlineEndX = forwardDilateX + 1;
                            }
                            ++run;
                        }
                        // Go through all the new voxels and add the minimum distance values.  These are in separate
                        // loops for better cache performance
                        float xDistance = forwardMinDistanceInside;
                        for( forwardX = firstForwardX + 1; forwardX <= forwardDilateX; ++forwardX ) {
                            xDistance += voxelLength;
                            outputMinDistanceInsideAccessor.add_element( xDistance );
                        }
                        xDistance = forwardMinDistanceOutside;
                        for( forwardX = firstForwardX + 1; forwardX <= forwardDilateX; ++forwardX ) {
                            xDistance += voxelLength;
                            outputMinDistanceOutsideAccessor.add_element( xDistance );
                        }
                        for( forwardX = firstForwardX + 1; forwardX <= forwardDilateX; ++forwardX ) {
                            /*
                            if( b >= 55 && b <= 56 && c == 5 )
                              cout << "+ forwardX " << levelSet.m_distanceData.size() << endl;
                            */
                            levelSet.m_distanceData.push_back( lsDistance );
                        }

                        // Duplicate the values for the rest of the named channels
                        for( size_t i = 0; i < outputChannels.size(); ++i ) {
                            const char* elementToCopy = inputChannels[i].data( forwardDataIndex );
                            rle_channel_general_accessor& outputChannel = outputChannels[i];
                            size_t primitiveSize = outputChannels[i].primitive_size();
                            for( forwardX = firstForwardX + 1; forwardX <= forwardDilateX; ++forwardX )
                                memcpy( outputChannel.add_element(), elementToCopy, primitiveSize );
                        }
                    }

                    //					psForwardDilate.exit();

                    //////////
                    // Second, copy the undefined runs in between
                    //////////

                    //					psUndefinedCopy.enter();

                    // Here we fill in the voxels that are in the open interval (forwardDilateX, backwardDilateX) using
                    // undefined runs
                    if( forwardDilateX + 1 < backwardDilateX ) {
                        // Find the first undefined run that intersects with this interval
                        run = prevDefinedRun + 1;
                        //						cout << "runRangeStart: " << runRangeStart << ",
                        //runRangeEnd: " << runRangeEnd << endl; 						cout << "prevDefinedRun: " << prevDefinedRun << ",
                        //curDefinedRun: " << curDefinedRun << endl; 						cout << "Gap undefined range: (" << forwardDilateX
                        //<< ", " << backwardDilateX << ")" << endl; 						cout << "Run " << run << ", x range [" <<
                        //tempRis.m_runData[run].x << ", " << tempRis.m_runData[run+1].x << ")" << endl;
                        while( run < curDefinedRun && tempRis.m_runData[run + 1].x <= forwardDilateX + 1 )
                            ++run;
                        // Add this run
                        if( run < curDefinedRun ) {
                            int x = max( forwardDilateX + 1, tempRis.m_runData[run].x );
                            if( x < backwardDilateX ) {
                                ris.m_runData.push_back( run_data( x, tempRis.m_runData[run].dataIndex ) );
                                creatingUndefinedRun = true;
                            }
                        }
                        ++run;
                        // Copy the rest of the runs until they don't intersect with this interval anymore
                        while( run < curDefinedRun && tempRis.m_runData[run].x <= backwardDilateX - 1 ) {
                            ris.m_runData.push_back( tempRis.m_runData[run] );
                            creatingUndefinedRun = true;
                            ++run;
                        }
                    }

                    // Start the defined run if necessary (i.e. an undefined run is being created and we haven't passed
                    // all the defined runs)
                    if( creatingUndefinedRun && backwardDataIndex >= 0 ) {
                        /*
                        if( b >= 55 && b <= 56 && c == 5 )
                          cout << "Starting run at backwardDilateX: " << backwardDilateX << endl;
                        */
                        ris.m_runData.push_back( run_data( backwardDilateX, (int)levelSet.m_distanceData.size() ) );
                        creatingUndefinedRun = false;
                    }

                    //					psUndefinedCopy.exit();

                    //////////
                    // Third, do the backward dilation
                    //////////

                    //					psBackwardDilate.enter();

                    if( backwardDilateX < backwardX ) {
                        //						cout << "Backward dilation of " << backwardX -
                        //backwardDilateX << " voxels" << endl;
                        float lsDistance =
                            tempLevelSet
                                .m_distanceData[backwardDataIndex]; // Keeping a simple distance function in play too

                        int voteCount;
                        run = curDefinedRun - 1;
                        while( run >= runRangeStart && tempRis.m_runData[run].x > backwardDilateX )
                            --run;
                        if( run >= runRangeStart )
                            voteCount = ( -tempRis.m_runData[run].dataIndex ) & 0xff;
                        else
                            voteCount = exteriorVoteCount;

                        // Go through all the undefined runs within the dilation region
                        int firstBackwardX = backwardDilateX;
                        ++run;
                        while( backwardDilateX < backwardX ) {
                            int runDilateX = min( tempRis.m_runData[run].x, backwardX );
                            // Go through all the voxels up to the start of this run
                            while( backwardDilateX < runDilateX ) {
                                // Create the values for all the special-cased channels
                                outputVoteCountAccessor.add_element( voteCount );

                                ++backwardDilateX;
                            }

                            voteCount = ( -tempRis.m_runData[run].dataIndex ) & 0xff;
                            ++run;
                        }
                        // Go through all the new voxels and add the minimum distance values.  These are in separate
                        // loops for better cache performance
                        float xDistance =
                            backwardMinDistanceInside + ( backwardDilateX - firstBackwardX ) * voxelLength;
                        for( backwardDilateX = firstBackwardX; backwardDilateX < backwardX; ++backwardDilateX ) {
                            outputMinDistanceInsideAccessor.add_element( xDistance );
                            xDistance -= voxelLength;
                        }
                        xDistance = backwardMinDistanceOutside + ( backwardDilateX - firstBackwardX ) * voxelLength;
                        for( backwardDilateX = firstBackwardX; backwardDilateX < backwardX; ++backwardDilateX ) {
                            outputMinDistanceOutsideAccessor.add_element( xDistance );
                            xDistance -= voxelLength;
                        }
                        for( backwardDilateX = firstBackwardX; backwardDilateX < backwardX; ++backwardDilateX ) {
                            /*
                            if( b >= 55 && b <= 56 && c == 5 )
                              cout << "+ backwardX " << levelSet.m_distanceData.size() << endl;
                            */
                            levelSet.m_distanceData.push_back( lsDistance );
                        }

                        // Duplicate the values for the rest of the named channels
                        for( size_t i = 0; i < outputChannels.size(); ++i ) {
                            const char* elementToCopy = inputChannels[i].data( backwardDataIndex );
                            rle_channel_general_accessor& outputChannel = outputChannels[i];
                            size_t primitiveSize = outputChannels[i].primitive_size();
                            for( backwardDilateX = firstBackwardX; backwardDilateX < backwardX; ++backwardDilateX )
                                memcpy( outputChannel.add_element(), elementToCopy, primitiveSize );
                        }
                    }

                    //					psBackwardDilate.exit();

                    //////////
                    // Fourth, copy the defined run
                    //////////

                    //					psDefinedCopy.enter();

                    if( backwardDataIndex >= 0 ) {
                        int runVoxelCount = tempRis.m_runData[curDefinedRun + 1].x - backwardX;
                        //						cout << "Defined run of " << runVoxelCount << " voxels" <<
                        //endl;
                        for( int i = 0; i < runVoxelCount; ++i ) {
                            /*
                            if( b >= 55 && b <= 56 && c == 5 )
                              cout << "+ defined " << levelSet.m_distanceData.size() << endl;
                            */
                            levelSet.m_distanceData.push_back( tempLevelSet.m_distanceData[backwardDataIndex + i] );
                        }
                        size_t dataSize = outputMinDistanceInsideAccessor.primitive_size() * runVoxelCount;
                        memcpy( outputMinDistanceInsideAccessor.m_data->add_element( dataSize ),
                                inputMinDistanceInsideAccessor.data( backwardDataIndex ), dataSize );
                        dataSize = outputMinDistanceOutsideAccessor.primitive_size() * runVoxelCount;
                        memcpy( outputMinDistanceOutsideAccessor.m_data->add_element( dataSize ),
                                inputMinDistanceOutsideAccessor.data( backwardDataIndex ), dataSize );
                        dataSize = outputVoteCountAccessor.primitive_size() * runVoxelCount;
                        memcpy( outputVoteCountAccessor.m_data->add_element( dataSize ),
                                inputVoteCountAccessor.data( backwardDataIndex ), dataSize );
                        for( size_t i = 0; i < outputChannels.size(); ++i ) {
                            dataSize = outputChannels[i].primitive_size() * runVoxelCount;
                            memcpy( outputChannels[i].m_data->add_element( dataSize ),
                                    inputChannels[i].data( backwardDataIndex ), dataSize );
                        }
                    }

                    //					psDefinedCopy.exit();

                    // Increment to the next gap between defined runs
                    prevDefinedRun = curDefinedRun;
                    ++curDefinedRun;
                    // Increment through the runs to find the next defined one
                    while( curDefinedRun <= runRangeEnd && tempRis.m_runData[curDefinedRun].dataIndex < 0 )
                        ++curDefinedRun;
                }

                /*
                if( b >= 55 && b <= 56 && c == 5 )
                  cout << "scanlineEndX: " << scanlineEndX << endl;
                */

                // Finish off the last run of the scanline
                ris.m_runData.push_back( run_data( scanlineEndX, -1 ) );

                //				psPerNonEmptyScanline.exit();
            } else {
                // Create an empty run for this BC coordinate
                ris.m_runData.push_back( run_data( 0, -1 ) );
                ris.m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = levelSet.m_distanceData.size();

    // run a consistency check
    stringstream sout;
    if( !ris.check_consistency( sout ) ) {
        ofstream fout( "c:\\temp\\bad.txt" );
        ris.dump( fout );
        throw std::runtime_error( "X Dilation: RLE Index Spec consistency check failed:\n" + sout.str() );
    }
    //	psTotal.exit();

    /*
    cout << endl << "X Dilation Stats" << endl;
    cout << psTotal << endl;
    cout << psInit << endl;
    cout << psPerNonEmptyScanline << endl;
    cout << psPerScanlineInit << endl;
    cout << psPerDefinedRunInit << endl;
    cout << psForwardDilate << endl;
    cout << psBackwardDilate << endl;
    cout << psUndefinedCopy << endl;
    cout << psDefinedCopy << endl;
    cout << endl;
    */
}

/**
 * This is an internal function, and therefore in the detail namespace.
 *
 * The function does a sweep along +X and -X for each run, propagating the minimum inside and minimum outside in a
 * simplified fast-sweeping like manner.
 *
 * The main reason for doing this is so that if you have multiple inconsistent intersections in a row, and voting throws
 * some of them out, the values from the other nearby intersections will be propagated and used.  This prevents getting
 * undefined voxels where you would expect a defined voxel in some cases.
 *
 * @param  levelSet           The level set object to apply the X sweeping to.
 *
 */
void mesh_to_level_set_region_x_sweeping( rle_level_set& levelSet ) {
    const rle_index_spec& ris = levelSet.get_rle_index_spec();

    float voxelLength = levelSet.get_voxel_coord_system().voxel_length();

    rle_channel_accessor<float> minDistanceInsideAccessor, minDistanceOutsideAccessor;
    minDistanceInsideAccessor = levelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    minDistanceOutsideAccessor = levelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );

    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            int bcIndex = b + c * ris.m_abcCoordSize.ysize();
            int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 2;

            // Go through all the runs to do the sweeping
            for( int run = runRangeStart; run <= runRangeEnd; ++run ) {
                int dataIndex = ris.m_runData[run].dataIndex;
                if( dataIndex >= 0 ) {
                    float* minDistanceInsideData = &minDistanceInsideAccessor[dataIndex];
                    float* minDistanceOutsideData = &minDistanceOutsideAccessor[dataIndex];

                    int xStart = ris.m_runData[run].x, xEnd = ris.m_runData[run + 1].x;
                    int x = xStart;
                    // Sweep forwards
                    float prevMinInsideValue = 1e38f, prevMinOutsideValue = 1e38f;
                    while( x < xEnd ) {
                        // Increment the min values by 1 voxel
                        prevMinInsideValue += voxelLength;
                        prevMinOutsideValue += voxelLength;

                        // Use the minimum of the current value at the voxel and the value we just computed
                        float curMinInsideValue = *minDistanceInsideData, curMinOutsideValue = *minDistanceOutsideData;
                        if( curMinInsideValue > prevMinInsideValue ) {
                            *minDistanceInsideData = prevMinInsideValue;
                        } else {
                            prevMinInsideValue = curMinInsideValue;
                        }
                        if( curMinOutsideValue > prevMinOutsideValue ) {
                            *minDistanceOutsideData = prevMinOutsideValue;
                        } else {
                            prevMinOutsideValue = curMinOutsideValue;
                        }

                        ++x;
                        ++minDistanceInsideData;
                        ++minDistanceOutsideData;
                    }
                    --x;
                    --minDistanceInsideData;
                    --minDistanceOutsideData;
                    prevMinInsideValue = 1e38f;
                    prevMinOutsideValue = 1e38f;
                    // Sweep backwards
                    while( x >= xStart ) {
                        // Increment the min values by 1 voxel
                        prevMinInsideValue += voxelLength;
                        prevMinOutsideValue += voxelLength;

                        // Use the minimum of the current value at the voxel and the value we just computed
                        float curMinInsideValue = *minDistanceInsideData, curMinOutsideValue = *minDistanceOutsideData;
                        if( curMinInsideValue > prevMinInsideValue ) {
                            *minDistanceInsideData = prevMinInsideValue;
                        } else {
                            prevMinInsideValue = curMinInsideValue;
                        }
                        if( curMinOutsideValue > prevMinOutsideValue ) {
                            *minDistanceOutsideData = prevMinOutsideValue;
                        } else {
                            prevMinOutsideValue = curMinOutsideValue;
                        }

                        --x;
                        --minDistanceInsideData;
                        --minDistanceOutsideData;
                    }
                }
            }
        }
    }
}

/**
 * This is an internal function, and therefore in the detail namespace.
 *
 * This merges two level sets together, adding together the votes being tallied for inside/outside status of a voxel,
 * as well as the named channels.  The votes counted in the rle compressed areas are also added together.
 *
 * An assumption made is that both level sets have an identical suite of named channels.  This is not a general-purpose
 * function, like the csg_union in the rle_level_set class, so it's perfectly reasonable to make that assumption.
 *
 * @param  levelSet        This is one of the input level sets, and is where the unioned output goes
 * @param  inputLevelSet   This is the other input level set.
 *
 */
void intermediate_mesh_to_level_set_union( rle_level_set& levelSet, const rle_level_set& inputLevelSet ) {
    // Swap the levelSet provided into a temporary variable, so that we've got two operands, rleFirst and rleSecond,
    // which we process into levelSet as the output.
    rle_level_set rleFirst( levelSet.get_voxel_coord_system() );
    rleFirst.swap( levelSet );
    const rle_level_set& rleSecond = inputLevelSet;

    rle_index_spec& ris = levelSet.m_rleIndex;

    // check that two operands have the same voxel coordinate systems
    if( rleFirst.m_voxelCoordSystem != rleSecond.m_voxelCoordSystem )
        throw std::runtime_error(
            "intermediate_mesh_to_level_set_union: The voxel coordinate systems of the two level set "
            "operands don't match.  The first is " +
            rleFirst.m_voxelCoordSystem.str() + ", while the second is " + rleSecond.m_voxelCoordSystem.str() + "." );

    // Add together the votes in the exterior region code.  The region code votes are ORed with 0x100, so that the lack
    // of votes is not 0, which is invalid as a region code.
    ris.m_exteriorRegionCode = -( 0x100 | ( ( ( -rleFirst.m_rleIndex.m_exteriorRegionCode ) & 0xff ) +
                                            ( ( -rleSecond.m_rleIndex.m_exteriorRegionCode ) & 0xff ) ) );

    // The voxel interface width of the result will be the smaller of the voxel interface widths of the two inputs.
    levelSet.m_interfaceVoxelWidthInside =
        ( std::min )( rleFirst.m_interfaceVoxelWidthInside, rleSecond.m_interfaceVoxelWidthInside );
    levelSet.m_interfaceVoxelWidthOutside =
        ( std::min )( rleFirst.m_interfaceVoxelWidthOutside, rleSecond.m_interfaceVoxelWidthOutside );
    levelSet.m_outsideDistance = ( std::min )( rleFirst.m_outsideDistance, rleSecond.m_outsideDistance );
    levelSet.m_insideDistance = ( std::max )( rleFirst.m_insideDistance, rleSecond.m_insideDistance );

    ////////////
    // Copy the bounding box from the first operand, and ensure that the second opeerand matches.
    ////////////
    ris.m_abcCoordOrigin = rleFirst.m_rleIndex.m_abcCoordOrigin;
    ris.m_abcCoordSize = rleFirst.m_rleIndex.m_abcCoordSize;
    if( ris.m_abcCoordOrigin != rleSecond.m_rleIndex.m_abcCoordOrigin ||
        ris.m_abcCoordSize != rleSecond.m_rleIndex.m_abcCoordSize )
        throw runtime_error(
            "intermediate_mesh_to_level_set_union: The voxel coordinate space bounding boxes of the two "
            "operands don't match.  This function requires that they match." );

    // clear any data already stored in ris
    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );

    ////////////
    // Determine what named channels the result will contain
    ////////////
    levelSet.m_namedChannels.clear();

    // Add a named channel matching each named channels in the first input
    for( map<frantic::tstring, rle_channel>::const_iterator i = rleFirst.m_namedChannels.begin();
         i != rleFirst.m_namedChannels.end(); ++i ) {
        levelSet.add_channel( i->second.name(), i->second.arity(), i->second.data_type() );
        // Check that this named channel exists in the second operand as well
        map<frantic::tstring, rle_channel>::const_iterator j = rleSecond.m_namedChannels.find( i->second.name() );
        if( j == rleSecond.m_namedChannels.end() )
            throw runtime_error(
                "intermediate_mesh_to_level_set_union: The two input level sets have mismatching named "
                "channels, including \"" +
                frantic::strings::to_string( i->second.name() ) + "\"." );
        else if( i->second.arity() != j->second.arity() || i->second.data_type() != j->second.data_type() )
            throw runtime_error(
                "intermediate_mesh_to_level_set_union: The two input level sets have an incompatible channel named \"" +
                frantic::strings::to_string( i->second.name() ) + "\".  The type of the first operand is " +
                frantic::strings::to_string(
                    channels::channel_data_type_str( j->second.arity(), j->second.data_type() ) ) +
                ", while the type of the second operand is " +
                frantic::strings::to_string(
                    channels::channel_data_type_str( i->second.arity(), i->second.data_type() ) ) +
                "." );
    }
    // The named channels in the second input should exactly match those in the first.  This finishes testing for this
    // condition.
    for( map<frantic::tstring, rle_channel>::const_iterator i = rleSecond.m_namedChannels.begin();
         i != rleSecond.m_namedChannels.end(); ++i ) {
        map<frantic::tstring, rle_channel>::const_iterator j = levelSet.m_namedChannels.find( i->second.name() );
        if( j == levelSet.m_namedChannels.end() )
            throw runtime_error(
                "intermediate_mesh_to_level_set_union: The two input level sets have mismatching named "
                "channels, including \"" +
                frantic::strings::to_string( i->second.name() ) + "\"." );
    }
    // Set up all the accessors for filling in the blended channel data.  We need three lists; the channels coming from
    // just the first or second operand, and the channels coming from both operands.
    vector<rle_channel_general_accessor> outputChannels;
    vector<const_rle_channel_general_accessor> inputChannelsFromFirst, inputChannelsFromSecond;
    for( map<frantic::tstring, rle_channel>::const_iterator i = levelSet.m_namedChannels.begin();
         i != levelSet.m_namedChannels.end(); ++i ) {
        // All channels except for three should just be added together.
        if( i->first != _T("G2LS_TEMP_MinDistanceInside") && i->first != _T("G2LS_TEMP_MinDistanceOutside") &&
            i->first != _T("G2LS_TEMP_VoteCount") ) {
            outputChannels.push_back( levelSet.get_channel_general_accessor( i->second.name() ) );
            inputChannelsFromFirst.push_back( rleFirst.get_channel_general_accessor( i->second.name() ) );
            inputChannelsFromSecond.push_back( rleSecond.get_channel_general_accessor( i->second.name() ) );
        }
    }

    // Auxiliary accessors specific to mesh to level set algorithm
    rle_channel_accessor<float> outputMinDistanceInsideAccessor, outputMinDistanceOutsideAccessor;
    rle_channel_accessor<int> outputVoteCountAccessor;
    outputMinDistanceInsideAccessor = levelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    outputMinDistanceOutsideAccessor = levelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    outputVoteCountAccessor = levelSet.get_channel_accessor<int>( _T("G2LS_TEMP_VoteCount") );
    const_rle_channel_accessor<float> inputMinDistanceInsideAccessorFirst, inputMinDistanceInsideAccessorSecond,
        inputMinDistanceOutsideAccessorFirst, inputMinDistanceOutsideAccessorSecond;
    const_rle_channel_accessor<int> inputVoteCountAccessorFirst, inputVoteCountAccessorSecond;
    inputMinDistanceInsideAccessorFirst = rleFirst.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    inputMinDistanceOutsideAccessorFirst = rleFirst.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    inputVoteCountAccessorFirst = rleFirst.get_channel_accessor<int>( _T("G2LS_TEMP_VoteCount") );
    inputMinDistanceInsideAccessorSecond = rleSecond.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    inputMinDistanceOutsideAccessorSecond = rleSecond.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    inputVoteCountAccessorSecond = rleSecond.get_channel_accessor<int>( _T("G2LS_TEMP_VoteCount") );

    ////////////
    // Iterate over the BCIndex values and combine the input scanlines
    ////////////

    // Allocate an array of blending alphas for combining the named channels
    vector<float> blendingAlpha;

    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            int bcIndex = b + c * ris.m_abcCoordSize.ysize();
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            // Get the range of runs for each.  Because the bounding boxes of the inputs line up, no fanciness is
            // required here.
            int runRangeStartFirst = rleFirst.m_rleIndex.m_bcToRunIndex[bcIndex],
                runRangeEndFirst = rleFirst.m_rleIndex.m_bcToRunIndex[bcIndex + 1] - 2;
            int runRangeStartSecond = rleSecond.m_rleIndex.m_bcToRunIndex[bcIndex],
                runRangeEndSecond = rleSecond.m_rleIndex.m_bcToRunIndex[bcIndex + 1] - 2;

            bool firstEmpty =
                runRangeStartFirst == runRangeEndFirst && rleFirst.m_rleIndex.m_runData[runRangeStartFirst].x ==
                                                              rleFirst.m_rleIndex.m_runData[runRangeStartFirst + 1].x;
            bool secondEmpty = runRangeStartSecond == runRangeEndSecond &&
                               rleSecond.m_rleIndex.m_runData[runRangeStartSecond].x ==
                                   rleSecond.m_rleIndex.m_runData[runRangeStartSecond + 1].x;

            if( firstEmpty && secondEmpty ) {
                // Create a zero-sized run, necessary for the coordinate queries to work properly
                ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
                ris.m_runData.push_back( run_data( 0, ris.m_exteriorRegionCode ) );
            } else {

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

                // However, if either of the runs is empty, we want to skip that one so that it is passed the end.
                // Note that we already dealt with the case where both are empty.
                if( firstEmpty ) {
                    runFirst++;
                    processingNextRunXFirst = ( std::numeric_limits<int>::max )();
                } else if( secondEmpty ) {
                    runSecond++;
                    processingNextRunXSecond = ( std::numeric_limits<int>::max )();
                }

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

                        if( runFirst <= runRangeEndFirst ) {
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

                        if( runSecond <= runRangeEndSecond ) {
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

                    // First deal with the case where both the inputs are undefined
                    if( processingUndefinedRunFirst && processingUndefinedRunSecond ) {
                        int regionCode = -( 0x100 | ( ( ( -processingRegionCodeFirst ) & 0xff ) +
                                                      ( ( -processingRegionCodeSecond ) & 0xff ) ) );
                        // In this case, both input runs are undefined, so we want to add the votes together
                        if( creatingUndefinedRun ) {
                            if( creatingRegionCode != regionCode ) {
                                ris.m_runData.push_back( run_data( creatingCurrentRunX, regionCode ) );
                                creatingRegionCode = regionCode;
                            }
                        } else {
                            ris.m_runData.push_back( run_data( creatingCurrentRunX, regionCode ) );
                            creatingUndefinedRun = true;
                            creatingRegionCode = regionCode;
                        }
                    }
                    // In this case, the second run is defined, but the first is undefined
                    else if( processingUndefinedRunFirst ) {
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back(
                                run_data( creatingCurrentRunX, (int)levelSet.m_distanceData.size() ) );
                            creatingUndefinedRun = false;
                        }
                        int dataIndexSecond = rleSecond.m_rleIndex.m_runData[runSecond - 1].dataIndex +
                                              ( creatingCurrentRunX - rleSecond.m_rleIndex.m_runData[runSecond - 1].x );
                        // Copy the data for all the channels in the arrays
                        for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            const_rle_channel_general_accessor& inputChannel = inputChannelsFromSecond[channelIndex];
                            // Simply copy it
                            size_t dataSize = segmentSize * channel.primitive_size();
                            memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( dataIndexSecond ),
                                    dataSize );
                        }
                        // Copy the data for the distance channel
                        int dataIndex = dataIndexSecond;
                        for( int i = 0; i < segmentSize; ++i )
                            levelSet.m_distanceData.push_back( rleSecond[dataIndex++] );
                        // Copy the minimum distance channels
                        size_t dataSize = segmentSize * inputMinDistanceInsideAccessorSecond.primitive_size();
                        memcpy( outputMinDistanceInsideAccessor.m_data->add_element( dataSize ),
                                inputMinDistanceInsideAccessorSecond.data( dataIndexSecond ), dataSize );
                        dataSize = segmentSize * inputMinDistanceOutsideAccessorSecond.primitive_size();
                        memcpy( outputMinDistanceOutsideAccessor.m_data->add_element( dataSize ),
                                inputMinDistanceOutsideAccessorSecond.data( dataIndexSecond ), dataSize );
                        // Combine the vote channel with the votes in the region code
                        int undefinedVote = ( -processingRegionCodeFirst ) & 0xff;
                        dataIndex = dataIndexSecond;
                        for( int i = 0; i < segmentSize; ++i )
                            outputVoteCountAccessor.add_element( inputVoteCountAccessorSecond[dataIndex++] +
                                                                 undefinedVote );
                    }
                    // In this case, the first run is defined, but the second is undefined
                    else if( processingUndefinedRunSecond ) {
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back(
                                run_data( creatingCurrentRunX, (int)levelSet.m_distanceData.size() ) );
                            creatingUndefinedRun = false;
                        }
                        int dataIndexFirst = rleFirst.m_rleIndex.m_runData[runFirst - 1].dataIndex +
                                             ( creatingCurrentRunX - rleFirst.m_rleIndex.m_runData[runFirst - 1].x );
                        // Copy the data for all the channels in the arrays
                        for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            const_rle_channel_general_accessor& inputChannel = inputChannelsFromFirst[channelIndex];
                            // Simply copy it
                            size_t dataSize = segmentSize * channel.primitive_size();
                            memcpy( channel.m_data->add_element( dataSize ), inputChannel.data( dataIndexFirst ),
                                    dataSize );
                        }
                        // Copy the data for the distance channel
                        int dataIndex = dataIndexFirst;
                        for( int i = 0; i < segmentSize; ++i )
                            levelSet.m_distanceData.push_back( rleFirst[dataIndex++] );
                        // Copy the minimum distance channels
                        size_t dataSize = segmentSize * inputMinDistanceInsideAccessorFirst.primitive_size();
                        memcpy( outputMinDistanceInsideAccessor.m_data->add_element( dataSize ),
                                inputMinDistanceInsideAccessorFirst.data( dataIndexFirst ), dataSize );
                        dataSize = segmentSize * inputMinDistanceOutsideAccessorFirst.primitive_size();
                        memcpy( outputMinDistanceOutsideAccessor.m_data->add_element( dataSize ),
                                inputMinDistanceOutsideAccessorFirst.data( dataIndexFirst ), dataSize );
                        // Combine the vote channel with the votes in the region code
                        int undefinedVote = ( -processingRegionCodeSecond ) & 0xff;
                        dataIndex = dataIndexFirst;
                        for( int i = 0; i < segmentSize; ++i )
                            outputVoteCountAccessor.add_element( inputVoteCountAccessorFirst[dataIndex++] +
                                                                 undefinedVote );
                    } else {
                        // In this case, both runs are defined, so we have to combine the data voxel-by-voxel
                        if( creatingUndefinedRun ) {
                            ris.m_runData.push_back(
                                run_data( creatingCurrentRunX, (int)levelSet.m_distanceData.size() ) );
                            creatingUndefinedRun = false;
                        }
                        int dataIndexFirstBegin =
                            rleFirst.m_rleIndex.m_runData[runFirst - 1].dataIndex +
                            ( creatingCurrentRunX - rleFirst.m_rleIndex.m_runData[runFirst - 1].x );
                        int dataIndexSecondBegin =
                            rleSecond.m_rleIndex.m_runData[runSecond - 1].dataIndex +
                            ( creatingCurrentRunX - rleSecond.m_rleIndex.m_runData[runSecond - 1].x );
                        float weights[2] = { 1, 1 };
                        // Add together the data for all the channels in the arrays
                        for( unsigned channelIndex = 0; channelIndex < outputChannels.size(); ++channelIndex ) {
                            rle_channel_general_accessor& channel = outputChannels[channelIndex];
                            const_rle_channel_general_accessor& inputChannelFirst =
                                inputChannelsFromFirst[channelIndex];
                            const_rle_channel_general_accessor& inputChannelSecond =
                                inputChannelsFromSecond[channelIndex];
                            channel_weighted_sum_combine_function_t ws = channel.m_weightedSumFunction;
                            // Use the weighted sum function to add the two input channels together
                            size_t dataSize = segmentSize * channel.primitive_size();
                            const char* data[2];
                            data[0] = inputChannelFirst.data( dataIndexFirstBegin );
                            data[1] = inputChannelSecond.data( dataIndexSecondBegin );
                            ws( weights, data, 2, segmentSize * channel.arity(),
                                channel.m_data->add_element( dataSize ) );
                        }
                        // Combine the data for the distance channel using the min operation
                        int dataIndexFirst = dataIndexFirstBegin, dataIndexSecond = dataIndexSecondBegin;
                        for( int i = 0; i < segmentSize; ++i )
                            levelSet.m_distanceData.push_back(
                                ( std::min )( rleFirst[dataIndexFirst++], rleSecond[dataIndexSecond++] ) );
                        // Combine the minimum distance channels using the min operation
                        dataIndexFirst = dataIndexFirstBegin;
                        dataIndexSecond = dataIndexSecondBegin;
                        for( int i = 0; i < segmentSize; ++i )
                            outputMinDistanceInsideAccessor.add_element(
                                ( std::min )( inputMinDistanceInsideAccessorFirst[dataIndexFirst++],
                                              inputMinDistanceInsideAccessorSecond[dataIndexSecond++] ) );
                        dataIndexFirst = dataIndexFirstBegin;
                        dataIndexSecond = dataIndexSecondBegin;
                        for( int i = 0; i < segmentSize; ++i )
                            outputMinDistanceOutsideAccessor.add_element(
                                ( std::min )( inputMinDistanceOutsideAccessorFirst[dataIndexFirst++],
                                              inputMinDistanceOutsideAccessorSecond[dataIndexSecond++] ) );
                        // Combine the vote channel by adding the votes together
                        dataIndexFirst = dataIndexFirstBegin;
                        dataIndexSecond = dataIndexSecondBegin;
                        for( int i = 0; i < segmentSize; ++i )
                            outputVoteCountAccessor.add_element( inputVoteCountAccessorFirst[dataIndexFirst++] +
                                                                 inputVoteCountAccessorSecond[dataIndexSecond++] );
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
    ris.m_dataSize = levelSet.m_distanceData.size();

    // run a consistency check
    stringstream sout;
    if( !ris.check_consistency( sout ) ) {
        throw std::runtime_error( "Intermediate Union: RLE Index Spec consistency check failed:\n" + sout.str() );
    }
}

/**
 * This is an internal function, and therefore in the detail namespace.
 *
 * This takes a level set which has had all its votes combined, and finalizes it into a normal distance-function based
 * level set. It expects there to be a number of auxiliary channels, as created by the function
 * convert_intersections_to_level_set.
 *
 * @param  levelSet  The level set to finalize.
 *
 */
void finalize_geometry_to_levelset_conversion( int exteriorRegionCode, rle_level_set& levelSet ) {
    // Swap the input level set into a temporary variable, so we can then write the result back into it.
    rle_level_set tempLevelSet( levelSet.get_voxel_coord_system() );
    levelSet.swap( tempLevelSet );
    levelSet.m_interfaceVoxelWidthInside = tempLevelSet.m_interfaceVoxelWidthInside;
    levelSet.m_interfaceVoxelWidthOutside = tempLevelSet.m_interfaceVoxelWidthOutside;
    levelSet.m_insideDistance = tempLevelSet.m_insideDistance;
    levelSet.m_outsideDistance = tempLevelSet.m_outsideDistance;

    // Make these thresholds a little bit liberal because this mesh to level set conversion overestimates the distance
    // values in some cases The factor is sqrt(3), because that's how much a diagonal distance could be off by.
    float interfaceWorldWidthInsideThreshold =
        1.73f * levelSet.m_interfaceVoxelWidthInside * levelSet.get_voxel_coord_system().voxel_length();
    float interfaceWorldWidthOutsideThreshold =
        1.73f * levelSet.m_interfaceVoxelWidthOutside * levelSet.get_voxel_coord_system().voxel_length();

    // The tie breaking region code is the opposite of the exterior region code.  This has gone through some debate, and
    // practical examples, and currently this seems like the best choice.
    // An exception to this rule is with a 0-0 tie, where we use the exterior region code instead.
    int tieBreakRegionCode = ( exteriorRegionCode == -1 ) ? -2 : -1;

    rle_index_spec& ris = levelSet.m_rleIndex;
    rle_index_spec& tempRis = tempLevelSet.m_rleIndex;

    // Reset the recorded data size to 0, so when we create the named channels they start with a size of zero
    ris.m_dataSize = 0;

    // Create the channels in the destination level set, and initialize accessors for them
    vector<const_rle_channel_general_accessor> inputChannels;
    vector<rle_channel_general_accessor> outputChannels;
    vector<channel_weighted_sum_combine_function_t> channelWeightedSumFunctions;
    for( map<frantic::tstring, rle_channel>::const_iterator i = tempLevelSet.m_namedChannels.begin(),
                                                            iterEnd = tempLevelSet.m_namedChannels.end();
         i != iterEnd; ++i ) {
        // Strip out all the temporary channels, which are prefixed by "G2LS_TEMP_"
        if( i->first.substr( 0, 10 ) != _T("G2LS_TEMP_") ) {
            inputChannels.push_back( tempLevelSet.get_channel_general_accessor( i->first ) );
            levelSet.add_channel( i->first, i->second.arity(), i->second.data_type() );
            outputChannels.push_back( levelSet.get_channel_general_accessor( i->first ) );
            channelWeightedSumFunctions.push_back( outputChannels.back().m_weightedSumFunction );
        }
    }

    // Auxiliary accessors specific to this conversion algorithm
    rle_channel_accessor<float> channelWeightAccessor, minDistanceInsideAccessor, minDistanceOutsideAccessor;
    rle_channel_accessor<int> voteCountAccessor;

    channelWeightAccessor = tempLevelSet.get_channel_accessor<float>( _T("G2LS_TEMP_ChannelWeight") );
    minDistanceInsideAccessor = tempLevelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceInside") );
    minDistanceOutsideAccessor = tempLevelSet.get_channel_accessor<float>( _T("G2LS_TEMP_MinDistanceOutside") );
    voteCountAccessor = tempLevelSet.get_channel_accessor<int>( _T("G2LS_TEMP_VoteCount") );

    ris.m_exteriorRegionCode = exteriorRegionCode;
    ris.m_abcCoordOrigin = tempRis.m_abcCoordOrigin;
    ris.m_abcCoordSize = tempRis.m_abcCoordSize;

    // Reserve the correct amount of space, and clear the arrays
    ris.m_bcToRunIndex.resize( ris.m_abcCoordSize.ysize() * ris.m_abcCoordSize.zsize() + 1 );

    ris.m_runData.clear();
    levelSet.m_distanceData.clear();

    vector3 xyz;
    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            xyz.y = b + ris.m_abcCoordOrigin.y;
            xyz.z = c + ris.m_abcCoordOrigin.z;

            int bcIndex = b + c * ris.m_abcCoordSize.ysize();
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            int runRangeStart = tempRis.m_bcToRunIndex[bcIndex], runRangeEnd = tempRis.m_bcToRunIndex[bcIndex + 1] - 2;
            // If it's a non-empty scanline, then process it
            if( runRangeStart != runRangeEnd ||
                tempRis.m_runData[runRangeStart].x != tempRis.m_runData[runRangeStart + 1].x ) {

                // The default state is an undefined region of the exterior region code
                bool creatingUndefinedRun = true;
                int creatingRegionCode = ris.m_exteriorRegionCode;

                for( int run = runRangeStart; run <= runRangeEnd; ++run ) {
                    const run_data& rd = tempRis.m_runData[run];
                    int dataIndex = rd.dataIndex;
                    // If it's a defined run, process it element by element
                    if( dataIndex >= 0 ) {
                        int xStart = rd.x, xEnd = tempRis.m_runData[run + 1].x;
                        for( int x = xStart; x < xEnd; ++x, ++dataIndex ) {
                            float levelSetDistance;
                            bool undefinedVoxel = false;
                            // Tally the inside/outside votes
                            int insideVotes = voteCountAccessor[dataIndex];
                            int outsideVotes = insideVotes & 0x0f;
                            insideVotes = ( insideVotes & 0xf0 ) >> 4;

                            // Using the votes, figure out the level set distance
                            if( insideVotes > outsideVotes ) {
                                levelSetDistance = -minDistanceInsideAccessor[dataIndex]; // If inside wins
                                if( levelSetDistance < -interfaceWorldWidthInsideThreshold )
                                    undefinedVoxel = true;
                            } else if( insideVotes < outsideVotes ) {
                                levelSetDistance = minDistanceOutsideAccessor[dataIndex]; // If outside wins
                                if( levelSetDistance > interfaceWorldWidthOutsideThreshold )
                                    undefinedVoxel = true;
                            } else if( tieBreakRegionCode == -1 ) {
                                // Use the tie breaker region code as the default when there's a tie, except when it's a
                                // 0-0 tie.
                                if( insideVotes > 0 ) {
                                    levelSetDistance = minDistanceOutsideAccessor[dataIndex]; // outside
                                    if( levelSetDistance > interfaceWorldWidthOutsideThreshold )
                                        undefinedVoxel = true;
                                } else {
                                    // In a 0-0 tie, make it an undefined voxel with the exterior region code
                                    levelSetDistance = exteriorRegionCode == -1 ? 1.f : -1.f;
                                    undefinedVoxel = true;
                                }
                            } else {
                                // Use the tie breaker region code as the default when there's a tie, except when it's a
                                // 0-0 tie.
                                if( insideVotes > 0 ) {
                                    levelSetDistance = -minDistanceInsideAccessor[dataIndex]; // inside
                                    if( levelSetDistance < -interfaceWorldWidthInsideThreshold )
                                        undefinedVoxel = true;
                                } else {
                                    // In a 0-0 tie, make it an undefined voxel with the exterior region code
                                    levelSetDistance = exteriorRegionCode == -1 ? 1.f : -1.f;
                                    undefinedVoxel = true;
                                }
                            }

                            if( undefinedVoxel ) {
                                int regionCode = levelSetDistance < 0 ? -2 : -1;

                                // Create the run index data entry for an undefined run
                                if( !creatingUndefinedRun ) {
                                    ris.m_runData.push_back( run_data( x, regionCode ) );
                                    creatingUndefinedRun = true;
                                    creatingRegionCode = regionCode;
                                } else if( regionCode != creatingRegionCode ) {
                                    // It's possible that adjacent runs will have different region codes
                                    ris.m_runData.push_back( run_data( x, regionCode ) );
                                    creatingRegionCode = regionCode;
                                }
                            } else {
                                // If necessary start a defined run
                                if( creatingUndefinedRun ) {
                                    ris.m_runData.push_back( run_data( x, (int)levelSet.m_distanceData.size() ) );
                                    creatingUndefinedRun = false;
                                }

                                // Set the level set distance
                                levelSet.m_distanceData.push_back( levelSetDistance );

                                // Get the normalizing factor for the named channels
                                float channelWeightInverse = 1.f / channelWeightAccessor[dataIndex];
                                const char* data;

                                // Copy all the named channel values
                                for( size_t i = 0; i < inputChannels.size(); ++i ) {
                                    data = inputChannels[i].data( dataIndex );
                                    channelWeightedSumFunctions[i]( &channelWeightInverse, &data, 1,
                                                                    inputChannels[i].arity(),
                                                                    outputChannels[i].add_element() );
                                }
                            }
                        }
                    }
                    // If it's an undefined run, process it as a single chunk
                    else {
                        // Figure out the region code based on the voting
                        int insideVotes = ( ( -dataIndex ) & 0xf0 ) >> 4, outsideVotes = ( -dataIndex ) & 0x0f;
                        int regionCode;
                        if( insideVotes > outsideVotes )
                            regionCode = -2; // If inside wins
                        else if( insideVotes < outsideVotes )
                            regionCode = -1; // If outside wins
                        else
                            regionCode =
                                tieBreakRegionCode; // Use the tie breaker region code as the default when there's a tie

                        // Create the run index data entry
                        if( !creatingUndefinedRun ) {
                            ris.m_runData.push_back( run_data( rd.x, regionCode ) );
                            creatingUndefinedRun = true;
                            creatingRegionCode = regionCode;
                        } else if( regionCode != creatingRegionCode ) {
                            // It's possible that adjacent runs will have different region codes
                            ris.m_runData.push_back( run_data( rd.x, regionCode ) );
                            creatingRegionCode = regionCode;
                        }
                    }
                }
                // If any runs were flagged, finish off the last run.  Otherwise add a zero-sized run
                if( (int)ris.m_runData.size() > ris.m_bcToRunIndex[bcIndex] ) {
                    // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to
                    // optimize the storage and access.
                    if( !( creatingUndefinedRun && creatingRegionCode == ris.m_exteriorRegionCode ) )
                        ris.m_runData.push_back( run_data( tempRis.m_runData[runRangeEnd + 1].x, -1 ) );
                } else {
                    ris.m_runData.push_back( run_data( 0, -1 ) );
                    ris.m_runData.push_back( run_data( 0, -1 ) );
                }
            }
            // Otherwise, it's an empty run
            else {
                ris.m_runData.push_back( run_data( 0, -1 ) );
                ris.m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }
    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();
    // Set the data size as well
    ris.m_dataSize = levelSet.m_distanceData.size();

    stringstream sout;
    if( !ris.check_consistency( sout ) ) {
        throw std::runtime_error( "Finalize: the RLEIndex Spec is no longer consistent:\n" + sout.str() );
    }

    // Set the inside and outside distances
    levelSet.m_outsideDistance =
        ( levelSet.m_interfaceVoxelWidthOutside + 1 ) * levelSet.get_voxel_coord_system().voxel_length();
    levelSet.m_insideDistance =
        -( levelSet.m_interfaceVoxelWidthInside + 1 ) * levelSet.get_voxel_coord_system().voxel_length();
}

/// This function prints the middle XY plane slice of the level set to cout
void debug_slice_print( ostream& out, const rle_level_set& levelSet, int z = INT_MAX ) {
    //	const_rle_channel_accessor<int> voteAccessor = levelSet.get_channel_accessor<int>("G2LS_TEMP_VoteCount");
    //	const_rle_channel_accessor<float> minDist =
    //levelSet.get_channel_accessor<float>("G2LS_TEMP_MinDistanceOutside");

    boundbox3 box = levelSet.get_rle_index_spec().outer_bounds();
    if( z == INT_MAX )
        z = (int)box.center().z;
    for( int y = box.minimum().y; y <= box.maximum().y; ++y ) {
        for( int x = box.minimum().x; x <= box.maximum().x; ++x ) {
            int dataIndex = levelSet.get_rle_index_spec().XYZtoDataIndex( vector3( x, y, z ) );
            if( dataIndex >= 0 ) {
                float distance = levelSet[dataIndex];
                out << distance << " ";
                /*
                int votes = voteAccessor[dataIndex];
                int insideVotes = ((votes&0xf0)>>4), outsideVotes = (votes&0x0f);
                char msg[3];
                if( insideVotes > outsideVotes ) {
                  msg[0] = (char)('A' + insideVotes);
                  msg[1] = (char)('A' + outsideVotes);
                } else {
                  msg[0] = (char)('a' + insideVotes);
                  msg[1] = (char)('a' + outsideVotes);
                }
                msg[2] = '\0';
                out << msg;
                */
            } else {
                out.setf( ios::hex, ios::basefield );
                out << setw( 2 ) << setfill( '0' ) << ( ( -dataIndex ) & 0xff ) << setfill( ' ' ) << setw( 0 );
                out.setf( ios::dec, ios::basefield );
            }
        }
        out << "\n";
    }
    out << endl;
}

void debug_slices_print( ostream& out, const rle_level_set& levelSet ) {
    boundbox3 box = levelSet.get_rle_index_spec().outer_bounds();
    for( int z = box.zminimum(); z <= box.zmaximum(); ++z ) {
        debug_slice_print( out, levelSet, z );
        out << "\n";
    }
}

} // namespace detail

// Converts geometry into a level set
void convert_geometry_to_levelset( const geometry::trimesh3& geometry, float interfaceVoxelDistanceInside,
                                   float interfaceVoxelDistanceOutside, const boundbox3& voxelBounds,
                                   rle_level_set& outLevelSet ) {
    if( interfaceVoxelDistanceInside >= 0 )
        throw runtime_error( "convert_geometry_to_levelset() - The interfaceVoxelDistanceInside provided must be "
                             "negative, but instead is " +
                             boost::lexical_cast<string>( interfaceVoxelDistanceInside ) + "." );
    if( interfaceVoxelDistanceOutside <= 0 )
        throw runtime_error( "convert_geometry_to_levelset() - The interfaceVoxelDistanceOutside provided must be "
                             "positive, but instead is " +
                             boost::lexical_cast<string>( interfaceVoxelDistanceOutside ) + "." );

    float interfaceVoxelWidthInside = -interfaceVoxelDistanceInside;
    float interfaceVoxelWidthOutside = interfaceVoxelDistanceOutside;

#ifdef _PRINT_GTOLS_PROFILING
    profiling_section psTotal( "Total" ), psScanConvert( "Scan Conversion" ), psScanToLS( "Scan Intersect to LS" );
    profiling_section psAxisPermute( "Axis Permutation" ), psInterUnion( "Intermediate Union" ),
        psXDilation( "X Dilation" ), psXSweeping( "X Sweeping" );
    profiling_section psFinalize( "Finalize" ), psRemoveInternal( "Remove Internal Intersections" );
#endif

#ifdef _MESH_DUMP
    static int counter = 0;
    ++counter;
    // Debugging code, save all the input state requested
    write_mesh_file( files::replace_sequence_number( "c:\\temp\\gtols\\gtoltest.xmesh", counter ), geometry );
    ofstream fout( files::replace_sequence_number( "c:\\temp\\gtols\\gtoltest.txt", counter ).c_str() );
    fout << "interfaceVoxelWidthInside: " << interfaceVoxelWidthInside << endl;
    fout << "interfaceVoxelWidthOutside: " << interfaceVoxelWidthOutside << endl;
    fout << "voxelBounds: " << voxelBounds << endl;
    fout << "outLevelSet.get_voxel_coord_system(): " << outLevelSet.get_voxel_coord_system() << endl;
    fout.close();
#endif

#ifdef _PRINT_GTOLS_PROFILING
    psTotal.enter();
#endif

    voxel_coord_system vcs = outLevelSet.get_voxel_coord_system();
    float voxelLength = vcs.voxel_length();

    boundbox3 conversionVoxelBounds;
    // If the voxel bounds we're given are empty, compute them based on the geometry
    if( voxelBounds.is_empty() ) {
        conversionVoxelBounds = vcs.get_voxel_bounds( geometry.compute_bound_box() );
        conversionVoxelBounds.expand( (int)ceil( max( interfaceVoxelWidthInside, interfaceVoxelWidthOutside ) ) );
    } else {
        conversionVoxelBounds = voxelBounds;
    }

    // Special case an empty voxel bounds
    if( conversionVoxelBounds.is_empty() ) {
        outLevelSet.set_to_empty();
        return;
    }

    //	cout << "voxel bounds " << conversionVoxelBounds << endl;
    //	cout << "volume: " << conversionVoxelBounds.get_volume() << endl;
    //	cout << "interfaceVoxelWidthInside: " << interfaceVoxelWidthInside << endl;
    //	cout << "interfaceVoxelWidthOutside: " << interfaceVoxelWidthOutside << endl;

    // Get the transform matrices for the space where the scan conversion will occur
    // NOTE: This may be confusing, because the RLE compression axis is X, but the scan conversion function scans along
    // Z.
    //       Calling the function convert_intersections_to_level_set implicitly does an axis permutation switching from
    //       Z to X.
    transform4f scanConvertTransformForZ =
        transform4f::from_translation( -conversionVoxelBounds.minimum() ) * vcs.world_to_voxel_transform();
    transform4f scanConvertTransformForX =
        transform4f::from_translation( 0, 0, conversionVoxelBounds.minimum().x * voxelLength ) *
        transform4f::from_scale( 1, 1, voxelLength ) * transform4f::from_axis_permutation( vector3( 2, 0, 1 ) ) *
        scanConvertTransformForZ;
    transform4f scanConvertTransformForY =
        transform4f::from_translation( 0, 0, conversionVoxelBounds.minimum().y * voxelLength ) *
        transform4f::from_scale( 1, 1, voxelLength ) * transform4f::from_axis_permutation( vector3( 1, 2, 0 ) ) *
        scanConvertTransformForZ;
    scanConvertTransformForZ = transform4f::from_translation( 0, 0, conversionVoxelBounds.minimum().z * voxelLength ) *
                               transform4f::from_scale( 1, 1, voxelLength ) * scanConvertTransformForZ;

    // prepare the intersection arrays
    vector<vector<scan_conversion_intersection>> intersectionDepthsX( conversionVoxelBounds.ysize() *
                                                                      conversionVoxelBounds.zsize() ),
        intersectionDepthsY( conversionVoxelBounds.xsize() * conversionVoxelBounds.zsize() ),
        intersectionDepthsZ( conversionVoxelBounds.xsize() * conversionVoxelBounds.ysize() );

#ifdef _PRINT_GTOLS_PROFILING
    psScanConvert.enter();
#endif
    // To support different kind of geometry than triangle meshes, an equivalent of trimesh3_scan_convert is the only
    // function which needs to be implemented for it.
    trimesh3_scan_convert( scanConvertTransformForX, geometry,
                           frantic::graphics2d::size2( conversionVoxelBounds.ysize(), conversionVoxelBounds.zsize() ),
                           intersectionDepthsX );
    trimesh3_scan_convert( scanConvertTransformForY, geometry,
                           frantic::graphics2d::size2( conversionVoxelBounds.zsize(), conversionVoxelBounds.xsize() ),
                           intersectionDepthsY );
    trimesh3_scan_convert( scanConvertTransformForZ, geometry,
                           frantic::graphics2d::size2( conversionVoxelBounds.xsize(), conversionVoxelBounds.ysize() ),
                           intersectionDepthsZ );
#ifdef _PRINT_GTOLS_PROFILING
    psScanConvert.exit();
#endif

#ifdef _PRINT_GTOLS_PROFILING
    psRemoveInternal.enter();
#endif
    detail::remove_internal_intersections( intersectionDepthsX );
    detail::remove_internal_intersections( intersectionDepthsY );
    detail::remove_internal_intersections( intersectionDepthsZ );
#ifdef _PRINT_GTOLS_PROFILING
    psRemoveInternal.exit();
#endif

    // Now figure out what the exterior region code should be, and whether this is a half-open surface
    int insideCountNegX = 0, outsideCountNegX = 0, insideCountPosX = 0, outsideCountPosX = 0;
    for( vector<vector<scan_conversion_intersection>>::const_iterator i = intersectionDepthsX.begin(),
                                                                      iterEnd = intersectionDepthsX.end();
         i != iterEnd; ++i ) {
        if( !i->empty() ) {
            if( i->front().normalFacingZPositive )
                ++insideCountNegX;
            else
                ++outsideCountNegX;

            if( i->back().normalFacingZPositive )
                ++outsideCountPosX;
            else
                ++insideCountPosX;
        }
    }
    int insideCountNegY = 0, outsideCountNegY = 0, insideCountPosY = 0, outsideCountPosY = 0;
    for( vector<vector<scan_conversion_intersection>>::const_iterator i = intersectionDepthsY.begin(),
                                                                      iterEnd = intersectionDepthsY.end();
         i != iterEnd; ++i ) {
        if( !i->empty() ) {
            if( i->front().normalFacingZPositive )
                ++insideCountNegY;
            else
                ++outsideCountNegY;

            if( i->back().normalFacingZPositive )
                ++outsideCountPosY;
            else
                ++insideCountPosY;
        }
    }
    int insideCountNegZ = 0, outsideCountNegZ = 0, insideCountPosZ = 0, outsideCountPosZ = 0;
    for( vector<vector<scan_conversion_intersection>>::const_iterator i = intersectionDepthsZ.begin(),
                                                                      iterEnd = intersectionDepthsZ.end();
         i != iterEnd; ++i ) {
        if( !i->empty() ) {
            if( i->front().normalFacingZPositive )
                ++insideCountNegZ;
            else
                ++outsideCountNegZ;

            if( i->back().normalFacingZPositive )
                ++outsideCountPosZ;
            else
                ++insideCountPosZ;
        }
    }

    int insideCount =
        insideCountNegX + insideCountPosX + insideCountNegY + insideCountPosY + insideCountNegZ + insideCountPosZ;
    int outsideCount =
        outsideCountNegX + outsideCountPosX + outsideCountNegY + outsideCountPosY + outsideCountNegZ + outsideCountPosZ;

    int exteriorRegionCode;
    bool isHalfOpenSurface = false;
    if( insideCount > outsideCount ) {
        exteriorRegionCode = -2;
        // If it's less than a 2:1 ratio, consider this a half-open surface
        if( insideCount / 2 < outsideCount )
            isHalfOpenSurface = true;
    } else {
        exteriorRegionCode = -1;
        // If it's less than a 2:1 ratio, consider this a half-open surface
        if( outsideCount / 2 < insideCount )
            isHalfOpenSurface = true;
    }

    // There's a special circumstance in which the exterior region code could jump between inside
    // and outside randomly.  This is if a flat plane has an ocean or noise animation on it, sometimes
    // more will be outside, sometimes more will be inside.  Most of the time, though, along one of the axes,
    // there will be an exact 50-50 split of the votes.  Here we detect that, end set the exterior region code
    // to outside in that case
    // TODO: Another option might be to just use an outside exterior code when it's a half open surface, and not give
    //       it the possibility of being interior in some circumstances.
    //	if( isHalfOpenSurface && exteriorRegionCode != -1 ) {
    //		if( insideCountX == outsideCountX || insideCountY == outsideCountY || insideCountZ == outsideCountZ )
    //			exteriorRegionCode = -1;
    //	}
    // Just going with the "outside" region code when encountering half-open surfaces for now.
    if( isHalfOpenSurface )
        exteriorRegionCode = -1;

    //	cout << "Inside votes: " << insideCount << ", outside votes: " << outsideCount << ", exterior region code: " <<
    // exteriorRegionCode << endl; 	cout << "is half open: " << isHalfOpenSurface << endl;
    rle_level_set tempLevelSet( outLevelSet.get_voxel_coord_system() );

    boundbox3 citlsBounds( conversionVoxelBounds.minimum().z, conversionVoxelBounds.maximum().z,
                           conversionVoxelBounds.minimum().x, conversionVoxelBounds.maximum().x,
                           conversionVoxelBounds.minimum().y, conversionVoxelBounds.maximum().y );
#ifdef _PRINT_GTOLS_PROFILING
    psScanToLS.enter();
#endif
    if( isHalfOpenSurface ) {
        int exteriorRegionCodeNeg = 0, exteriorRegionCodePos = 0;

        int maximumVotes = citlsBounds.ysize() * citlsBounds.zsize();
        //		cout << "Maximum Z Votes: " << maximumVotes << endl;
        //		cout << "insideCountNegZ: " << insideCountNegZ << ", outsideCountNegZ: " << outsideCountNegZ <<
        //endl; 		cout << "insideCountPosZ: " << insideCountPosZ << ", outsideCountPosZ: " << outsideCountPosZ << endl;
        // Only classify the directions if at least 50% of the votes got used
        if( insideCountNegZ + outsideCountNegZ > maximumVotes * 3 / 8 ) {
            // If we got a strong majority, use that to determine the exterior region codes for the negative direction
            // along Z
            if( insideCountNegZ > ( insideCountNegZ + outsideCountNegZ ) * 7 / 8 )
                exteriorRegionCodeNeg = -2;
            else if( outsideCountNegZ > ( insideCountNegZ + outsideCountNegZ ) * 7 / 8 )
                exteriorRegionCodeNeg = -1;

            // If we got a strong majority, use that to determine the exterior region codes for the positive direction
            // along Z
            if( insideCountPosZ > ( insideCountPosZ + outsideCountPosZ ) * 7 / 8 )
                exteriorRegionCodePos = -2;
            else if( outsideCountPosZ > ( insideCountPosZ + outsideCountPosZ ) * 7 / 8 )
                exteriorRegionCodePos = -1;

            //			cout << "ExteriorRegionCodeNeg: " << exteriorRegionCodeNeg << endl;
            //			cout << "ExteriorRegionCodePos: " << exteriorRegionCodePos << endl;
        }
        detail::convert_intersections_to_level_set(
            intersectionDepthsZ, geometry, interfaceVoxelWidthInside, interfaceVoxelWidthOutside, isHalfOpenSurface,
            exteriorRegionCodeNeg, exteriorRegionCodePos, citlsBounds, outLevelSet );
    } else {
        detail::convert_intersections_to_level_set( intersectionDepthsZ, geometry, interfaceVoxelWidthInside,
                                                    interfaceVoxelWidthOutside, isHalfOpenSurface, exteriorRegionCode,
                                                    exteriorRegionCode, citlsBounds, outLevelSet );
    }
#ifdef _PRINT_GTOLS_PROFILING
    psScanToLS.exit();
#endif

#ifdef _PRINT_GTOLS_PROFILING
    psXSweeping.enter();
#endif
    detail::mesh_to_level_set_region_x_sweeping( outLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psXSweeping.exit();
#endif

    /*
    {
      ofstream sout("c:\\temp\\slicesZ.txt");
      tempLevelSet = outLevelSet;
      tempLevelSet.apply_axis_permutation(vector3(1,2,0));
      detail::debug_slices_print(sout, tempLevelSet);
    }
    */

    //	detail::debug_slice_print(outLevelSet);

#ifdef _PRINT_GTOLS_PROFILING
    psAxisPermute.enter();
#endif
    outLevelSet.apply_axis_permutation( vector3( 2, 0, 1 ) );
#ifdef _PRINT_GTOLS_PROFILING
    psAxisPermute.exit();
#endif

//	detail::debug_slice_print(outLevelSet);
#ifdef _PRINT_GTOLS_PROFILING
    psXDilation.enter();
#endif
    detail::mesh_to_level_set_region_x_dilation( outLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psXDilation.exit();
#endif
    //	if( !outLevelSet.get_rle_index_spec().check_consistency(cout) )
    //		throw runtime_error( "inconsistency" );
    //	detail::debug_slice_print(outLevelSet);

    /*
    {
      ofstream sout("c:\\temp\\slicesZ_dilated.txt");
      tempLevelSet = outLevelSet;
      tempLevelSet.apply_axis_permutation(vector3(2,0,1));
      detail::debug_slices_print(sout, tempLevelSet);
    }
    */

    /*
    tempLevelSet = outLevelSet;
    tempLevelSet.apply_axis_permutation(vector3(2,0,1));
    detail::finalize_geometry_to_levelset_conversion( exteriorRegionCode, tempLevelSet );
    legacyflood::write_legacy_rle_level_set( "c:\\temp\\zscan.cssf", tempLevelSet );
    */

    //	detail::debug_slice_print(outLevelSet);

    citlsBounds.set( conversionVoxelBounds.minimum().y, conversionVoxelBounds.maximum().y,
                     conversionVoxelBounds.minimum().z, conversionVoxelBounds.maximum().z,
                     conversionVoxelBounds.minimum().x, conversionVoxelBounds.maximum().x );
#ifdef _PRINT_GTOLS_PROFILING
    psScanToLS.enter();
#endif
    if( isHalfOpenSurface ) {
        int exteriorRegionCodeNeg = 0, exteriorRegionCodePos = 0;
        //*
        int maximumVotes = citlsBounds.ysize() * citlsBounds.zsize();
        //		cout << "Maximum Y Votes: " << maximumVotes << endl;
        //		cout << "insideCountNegY: " << insideCountNegY << ", outsideCountNegY: " << outsideCountNegY <<
        //endl; 		cout << "insideCountPosY: " << insideCountPosY << ", outsideCountPosY: " << outsideCountPosY << endl;
        // Only classify the directions if at least 50% of the votes got used
        if( insideCountNegY + outsideCountNegY > maximumVotes * 3 / 8 ) {
            // If we got a strong majority, use that to determine the exterior region codes for the negative direction
            // along Y
            if( insideCountNegY > ( insideCountNegY + outsideCountNegY ) * 7 / 8 )
                exteriorRegionCodeNeg = -2;
            else if( outsideCountNegY > ( insideCountNegY + outsideCountNegY ) * 7 / 8 )
                exteriorRegionCodeNeg = -1;

            // If we got a strong majority, use that to determine the exterior region codes for the positive direction
            // along Y
            if( insideCountPosY > ( insideCountPosY + outsideCountPosY ) * 7 / 8 )
                exteriorRegionCodePos = -2;
            else if( outsideCountPosY > ( insideCountPosY + outsideCountPosY ) * 7 / 8 )
                exteriorRegionCodePos = -1;

            //			cout << "ExteriorRegionCodeNeg: " << exteriorRegionCodeNeg << endl;
            //			cout << "ExteriorRegionCodePos: " << exteriorRegionCodePos << endl;
        }
        //*/
        detail::convert_intersections_to_level_set(
            intersectionDepthsY, geometry, interfaceVoxelWidthInside, interfaceVoxelWidthOutside, isHalfOpenSurface,
            exteriorRegionCodeNeg, exteriorRegionCodePos, citlsBounds, tempLevelSet );
    } else {
        detail::convert_intersections_to_level_set( intersectionDepthsY, geometry, interfaceVoxelWidthInside,
                                                    interfaceVoxelWidthOutside, isHalfOpenSurface, exteriorRegionCode,
                                                    exteriorRegionCode, citlsBounds, tempLevelSet );
    }
#ifdef _PRINT_GTOLS_PROFILING
    psScanToLS.exit();
#endif

#ifdef _PRINT_GTOLS_PROFILING
    psInterUnion.enter();
#endif
    detail::intermediate_mesh_to_level_set_union( outLevelSet, tempLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psInterUnion.exit();
#endif

#ifdef _PRINT_GTOLS_PROFILING
    psXSweeping.enter();
#endif
    detail::mesh_to_level_set_region_x_sweeping( outLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psXSweeping.exit();
#endif

    /*
    {
      ofstream sout("c:\\temp\\slicesY.txt");
      tempLevelSet.apply_axis_permutation(vector3(2,0,1));
      detail::debug_slices_print(sout, tempLevelSet);
    }
    */

#ifdef _PRINT_GTOLS_PROFILING
    psAxisPermute.enter();
#endif
    outLevelSet.apply_axis_permutation( vector3( 2, 0, 1 ) );
#ifdef _PRINT_GTOLS_PROFILING
    psAxisPermute.exit();
#endif

//	detail::debug_slice_print(outLevelSet);
#ifdef _PRINT_GTOLS_PROFILING
    psXDilation.enter();
#endif
    detail::mesh_to_level_set_region_x_dilation( outLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psXDilation.exit();
#endif

    /*
    {
      ofstream sout("c:\\temp\\slicesYZ_dilated.txt");
      tempLevelSet = outLevelSet;
      detail::debug_slices_print(sout, tempLevelSet);
    }
    */

    //	outLevelSet.get_rle_index_spec().dump(cout);
    //	if( !outLevelSet.get_rle_index_spec().check_consistency(cout) )
    //		throw runtime_error( "inconsistency" );
    //	detail::debug_slice_print(outLevelSet);

    /*
    tempLevelSet = outLevelSet;
    finalize_geometry_to_levelset_conversion( exteriorRegionCode, tempLevelSet );
    legacyflood::write_legacy_rle_level_set( "c:\\temp\\zscan_unioned_yscan.cssf", tempLevelSet );
    */

#ifdef _PRINT_GTOLS_PROFILING
    psScanToLS.enter();
#endif
    if( isHalfOpenSurface ) {
        int exteriorRegionCodeNeg = 0, exteriorRegionCodePos = 0;

        int maximumVotes = conversionVoxelBounds.ysize() * conversionVoxelBounds.zsize();
        //		cout << "Maximum X Votes: " << maximumVotes << endl;
        //		cout << "insideCountNegX: " << insideCountNegX << ", outsideCountNegX: " << outsideCountNegX <<
        //endl; 		cout << "insideCountPosX: " << insideCountPosX << ", outsideCountPosX: " << outsideCountPosX << endl;
        // Only classify the directions if at least 50% of the votes got used
        if( insideCountNegX + outsideCountNegX > maximumVotes * 3 / 8 ) {
            // If we got a strong majority, use that to determine the exterior region codes for the negative direction
            // along X
            if( insideCountNegX > ( insideCountNegX + outsideCountNegX ) * 7 / 8 )
                exteriorRegionCodeNeg = -2;
            else if( outsideCountNegX > ( insideCountNegX + outsideCountNegX ) * 7 / 8 )
                exteriorRegionCodeNeg = -1;

            // If we got a strong majority, use that to determine the exterior region codes for the positive direction
            // along X
            if( insideCountPosX > ( insideCountPosX + outsideCountPosX ) * 7 / 8 )
                exteriorRegionCodePos = -2;
            else if( outsideCountPosX > ( insideCountPosX + outsideCountPosX ) * 7 / 8 )
                exteriorRegionCodePos = -1;

            //			cout << "ExteriorRegionCodeNeg: " << exteriorRegionCodeNeg << endl;
            //			cout << "ExteriorRegionCodePos: " << exteriorRegionCodePos << endl;
        }
        detail::convert_intersections_to_level_set(
            intersectionDepthsX, geometry, interfaceVoxelWidthInside, interfaceVoxelWidthOutside, isHalfOpenSurface,
            exteriorRegionCodeNeg, exteriorRegionCodePos, conversionVoxelBounds, tempLevelSet );
    } else {
        detail::convert_intersections_to_level_set( intersectionDepthsX, geometry, interfaceVoxelWidthInside,
                                                    interfaceVoxelWidthOutside, isHalfOpenSurface, exteriorRegionCode,
                                                    exteriorRegionCode, conversionVoxelBounds, tempLevelSet );
    }
#ifdef _PRINT_GTOLS_PROFILING
    psScanToLS.exit();
#endif

    /*
    {
      ofstream sout("c:\\temp\\slicesX.txt");
      detail::debug_slices_print(sout, tempLevelSet);
    }
    */

#ifdef _PRINT_GTOLS_PROFILING
    psInterUnion.enter();
#endif
    detail::intermediate_mesh_to_level_set_union( outLevelSet, tempLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psInterUnion.exit();
#endif

#ifdef _PRINT_GTOLS_PROFILING
    psXSweeping.enter();
#endif
    detail::mesh_to_level_set_region_x_sweeping( outLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psXSweeping.exit();
#endif

    /*
    {
      ofstream sout("c:\\temp\\slices.txt");
      detail::debug_slices_print(sout, outLevelSet);
    }
    */

    //	detail::debug_slice_print(outLevelSet);

    /*
    finalize_geometry_to_levelset_conversion( exteriorRegionCode, tempLevelSet );
    legacyflood::write_legacy_rle_level_set( "c:\\temp\\xscan.cssf", tempLevelSet );
    */

#ifdef _PRINT_GTOLS_PROFILING
    psFinalize.enter();
#endif
    detail::finalize_geometry_to_levelset_conversion( exteriorRegionCode, outLevelSet );
#ifdef _PRINT_GTOLS_PROFILING
    psFinalize.exit();
#endif

#ifdef _MESH_DUMP
    // Debugging code, save the resulting level set
    write_rls_rle_level_set_file( files::replace_sequence_number( "c:\\temp\\gtols\\gtoltest.rls", counter ),
                                  outLevelSet );
#endif

    /*
    {
      ofstream sout("c:\\temp\\slices_finished.txt");
      detail::debug_slices_print(sout, outLevelSet);
    }
    */

#ifdef _PRINT_GTOLS_PROFILING
    psTotal.exit();
#endif

    //	legacyflood::write_legacy_rle_level_set( "c:\\temp\\zscan_unioned_yscan_unioned_xscan.cssf", outLevelSet );

#ifdef _PRINT_GTOLS_PROFILING
    cout << endl;
    cout << "Mesh to Level Set Profiling Statistics:" << endl;
    cout << psTotal << endl;
    cout << psScanConvert << endl;
    cout << psRemoveInternal << endl;
    cout << psScanToLS << endl;
    cout << psAxisPermute << endl;
    cout << psXDilation << endl;
    cout << psXSweeping << endl;
    cout << psInterUnion << endl;
    cout << psFinalize << endl;
    cout << endl;
#endif
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
