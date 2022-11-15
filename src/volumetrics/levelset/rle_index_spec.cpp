// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_pairwise_run_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_run_iterator.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::graphics;

namespace frantic {
namespace volumetrics {
namespace levelset {

// Dumps the contents of this rle_index_spec.
void rle_index_spec::dump( std::ostream& out ) const {
    out << "----- Dumping RLE Index Spec\n";
    out << " ABC Coordinate minimum: " << m_abcCoordOrigin << "\n";
    out << " ABC Coordinate maximum: " << m_abcCoordOrigin + m_abcCoordSize - vector3( 1 ) << "\n";
    out << " ABC Coordinate size: " << m_abcCoordSize << "\n";
    out << " Outer bounding box volume: " << m_abcCoordSize.volume() << "\n";
    out << " Region code of exterior: " << m_exteriorRegionCode << " (\""
        << ( ( m_exteriorRegionCode == -1 ) ? "outside" : "inside" ) << "\")\n";
    out << " Data size: " << m_dataSize << "\n";
    out << " bcToRunIndex array size: " << m_bcToRunIndex.size() << "\n";
    out << " m_runData array size: " << m_runData.size() << "\n";
    out << " m_cachedRisAdjacency allocated: " << ( ( m_cachedRisAdjacency == 0 ) ? "no\n" : "yes\n" );
    out << "\n";

    //	for( unsigned i = 0; i < m_bcToRunIndex.size(); ++i )
    //		out << i << ": " << m_bcToRunIndex[i] << "\n";

    out << "YZ/BC coordinates and runs:\n";
    boundbox3 outerBounds = outer_bounds();
    for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
        for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
            int b = y_to_b( y ), c = z_to_c( z );
            out << " YZ (" << y << ", " << z << ") / BC (" << b << ", " << c << ") / bcIndex "
                << ( b + c * m_abcCoordSize.ysize() ) << endl;
            // Iterate through all the runs in this scanline
            for( rle_run_iterator i( *this, b, c ), ie; i != ie; ++i ) {
                out << "   Run with X range [" << i.get_xmin() << ", " << i.get_xmax() << "].\n";
                if( i.get_data_index() < 0 )
                    out << "     Undefined with region code " << i.get_data_index() << ".\n";
                else
                    out << "     Defined with data index range [" << i.get_data_index() << ", "
                        << ( i.get_data_index() + i.get_xsize() - 1 ) << "].\n";
            }
        }
    }

    out << "----- Finished Dumping RLE Index Spec" << endl;
}

// This function checks that the rle index spec is self-consistent
bool rle_index_spec::check_consistency( std::ostream& out ) const {
    // If the index spec has no runs and is zero volume, it's empty and that means consistent
    if( m_abcCoordSize.volume() == 0 && m_runData.size() == 0 ) {
        return true;
    }

    if( m_exteriorRegionCode >= 0 ) {
        out << "The m_exteriorRegionCode value cannot be 0 or positive, but its current value is "
            << m_exteriorRegionCode << ".\n";
        return false;
    }

    // The m_bcToRunIndex array should be sized exactly (m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1)
    if( (int)m_bcToRunIndex.size() != ( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 ) ) {
        out << "The m_bcToRunIndex array size is not one more than the area of the box projected onto the YZ plane + 1."
            << endl;
        out << "It is " << m_bcToRunIndex.size() << " instead of the expected "
            << ( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 ) << "." << endl;
        return false;
    }

    // The values in m_bcToRunIndex should all be within the range [0, m_runData.size()), except for the last one which
    // should be equal to m_runData.size(). These values should also be strictly increasing. The first value in the
    // m_bcToRunIndex array should be zero.
    for( unsigned i = 0; i < m_bcToRunIndex.size() - 1; ++i ) {
        if( m_bcToRunIndex[i] < 0 || m_bcToRunIndex[i] >= (int)m_runData.size() ) {
            out << "The m_bcToRunIndex array value at index " << i
                << " contains a run index value which is out of bounds." << endl;
            out << "It is " << m_bcToRunIndex[i] << " but should be within the range [0, " << ( m_runData.size() - 1 )
                << "]." << endl;
            return false;
        }
        if( m_bcToRunIndex[i] >= m_bcToRunIndex[i + 1] ) {
            out << "The m_bcToRunIndex array value at index " << i << " is not less than the array value at index "
                << i + 1 << "." << endl;
            out << "At " << i << " the value is " << m_bcToRunIndex[i] << " and at " << i + 1 << " it is "
                << m_bcToRunIndex[i + 1] << "." << endl;
            return false;
        }
    }
    if( m_bcToRunIndex.front() != 0 ) {
        out << "The first value in m_bcToRunIndex is not equal to zero." << endl;
        out << "It is " << m_bcToRunIndex.front() << ", but should be 0." << endl;
        return false;
    }
    if( m_bcToRunIndex.back() != (int)m_runData.size() ) {
        out << "The last value in m_bcToRunIndex is not equal to one past the end of m_runData.size()." << endl;
        out << "It is " << m_bcToRunIndex.back() << ", but should be " << m_runData.size() << "." << endl;
        return false;
    }

    // We use the data index tracker to make sure all the defined runs are placed end to end, which means there should
    // be no gaps.
    int dataIndexTracker = 0, dataIndexTrackerRun = -1;

    // Now iterate through all the runs, and ensure the m_runData array contains reasonable data.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = m_bcToRunIndex[bcIndex], runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;

            // Each BC coordinate must have at least one (possibly empty) run associated with it.
            if( runRangeStart + 1 > runRangeEnd ) {
                out << "BC coordinate (" << b << ", " << c << ") does not have a run associated with it." << endl;
                out << "It should have, as a minimum, a single empty run." << endl;
                return false;
            }
            // If it's not a single empty run, then do some additional validation
            if( runRangeStart + 1 != runRangeEnd || m_runData[runRangeStart].x != m_runData[runRangeStart + 1].x ) {
                int previousRunDefined = false;
                for( int run = runRangeStart; run != runRangeEnd; ++run ) {
                    // All of these runs should have non-zero size.
                    if( m_runData[run].x >= m_runData[run + 1].x ) {
                        out << "Run " << run << ", with dataIndex " << m_runData[run].dataIndex << ", in run range ["
                            << runRangeStart << ", " << runRangeEnd << "), for BC coordinate (" << b << ", " << c
                            << ") and YZ coordinate (" << b + m_abcCoordOrigin.y << ", " << c + m_abcCoordOrigin.z
                            << "), is empty or has a decreasing value." << endl;
                        out << "Its X values are the invalid half-open range [" << m_runData[run].x << ", "
                            << m_runData[run + 1].x << ")." << endl;
                        return false;
                    }
                    // Make sure that the X values fall within the bounds defined by the structure
                    if( m_runData[run].x < m_abcCoordOrigin.x ||
                        m_runData[run].x >= m_abcCoordOrigin.x + m_abcCoordSize.xsize() ) {
                        out << "Run " << run << " in run range [" << runRangeStart << ", " << runRangeEnd
                            << "), for BC coordinate (" << b << ", " << c << ") and YZ coordinate ("
                            << b + m_abcCoordOrigin.y << ", " << c + m_abcCoordOrigin.z
                            << "), has a starting X value outside the bounding box of this structure." << endl;
                        out << "The X value is " << m_runData[run].x << ", but the valid range is ["
                            << m_abcCoordOrigin.x << ", " << m_abcCoordOrigin.x + m_abcCoordSize.xsize() - 1 << "]."
                            << endl;
                        return false;
                    }
                    if( m_runData[run + 1].x - 1 < m_abcCoordOrigin.x ||
                        m_runData[run + 1].x - 1 >= m_abcCoordOrigin.x + m_abcCoordSize.xsize() ) {
                        out << "Run " << run << " in run range [" << runRangeStart << ", " << runRangeEnd
                            << "), for BC coordinate (" << b << ", " << c << ") and YZ coordinate ("
                            << b + m_abcCoordOrigin.y << ", " << c + m_abcCoordOrigin.z
                            << "), has an ending X value outside the bounding box of this structure." << endl;
                        out << "The X value is " << m_runData[run + 1].x - 1 << ", but the valid range is ["
                            << m_abcCoordOrigin.x << ", " << m_abcCoordOrigin.x + m_abcCoordSize.xsize() - 1 << "]."
                            << endl;
                        return false;
                    }
                    if( m_runData[run].dataIndex >= 0 ) {
                        // Since this run is defined, it represents an array of data values.  Make sure that it is
                        // within the range defined by m_dataSize.
                        int runEnd = m_runData[run].dataIndex + ( m_runData[run + 1].x - m_runData[run].x ) - 1;
                        if( runEnd >= (int)m_dataSize ) {
                            out << "Run " << run << " in run range [" << runRangeStart << ", " << runRangeEnd
                                << "), for BC coordinate (" << b << ", " << c << ") and YZ coordinate ("
                                << b + m_abcCoordOrigin.y << ", " << c + m_abcCoordOrigin.z
                                << "), has dataIndices extending past the end of the defined data array size." << endl;
                            out << "The last data index in this run is " << runEnd << ", but the defined data size is "
                                << m_dataSize << endl;
                            return false;
                        }
                        // Update the computedDataSize
                        if( m_runData[run].dataIndex != dataIndexTracker ) {
                            out << "Run " << run << " in run range [" << runRangeStart << ", " << runRangeEnd
                                << "), for BC coordinate (" << b << ", " << c << ") and YZ coordinate ("
                                << b + m_abcCoordOrigin.y << ", " << c + m_abcCoordOrigin.z
                                << "), has a dataIndex gap just preceding it." << endl;
                            out << "The previous defined run was " << dataIndexTrackerRun << endl;
                            return false;
                        }
                        dataIndexTracker = runEnd + 1;
                        dataIndexTrackerRun = run;

                        if( previousRunDefined ) {
                            out << "Run " << run << " in run range [" << runRangeStart << ", " << runRangeEnd
                                << "), for BC coordinate (" << b << ", " << c << ") and YZ coordinate ("
                                << b + m_abcCoordOrigin.y << ", " << c + m_abcCoordOrigin.z
                                << "), is a defined run, and has a defined run preceding it." << endl;
                            out << "Adjacent defined runs must always be coalesced into a single defined run." << endl;
                            return false;
                        }
                        previousRunDefined = true;
                    } else {
                        previousRunDefined = false;
                    }
                }
                // Also check that the terminator of the last run is in good shape
                if( m_runData[runRangeEnd].x <= m_abcCoordOrigin.x ||
                    m_runData[runRangeEnd].x > m_abcCoordOrigin.x + m_abcCoordSize.xsize() ) {
                    out << "The last run of run range [" << runRangeStart << ", " << runRangeEnd
                        << "), for BC coordinate (" << b << ", " << c << ") and YZ coordinate ("
                        << b + m_abcCoordOrigin.y << ", " << c + m_abcCoordOrigin.z
                        << "), has an invalid ending X value." << endl;
                    out << "The X value is " << m_runData[runRangeEnd].x << ", but the valid range is ["
                        << m_abcCoordOrigin.x + 1 << ", " << m_abcCoordOrigin.x + m_abcCoordSize.xsize() << "]."
                        << endl;
                    return false;
                }
            } else {
                // In the case of a single empty run, the data index should indicate undefined
                if( m_runData[runRangeStart].dataIndex >= 0 ) {
                    out << "Run " << runRangeStart
                        << ", which is a single zero-sized run for an empty scanline, had dataIndex "
                        << m_runData[runRangeStart].dataIndex
                        << ", a positive value.  It should be negative, indicating undefined.";
                    return false;
                }
            }
        }
    }

    // Ensure that the size computed by iterating over the runs matches the value in the member variable
    if( dataIndexTracker != (int)m_dataSize ) {
        out << "The data size computed by iterating over the runs doesn't match the member value data size." << endl;
        out << "The member variable says the data size is " << m_dataSize << ", but the computed size is "
            << dataIndexTracker << endl;
        return false;
    }

    return true;
}

bool rle_index_spec::check_voxels_slow( std::ostream& out ) const {
    int boundsStartX = m_abcCoordOrigin.x, boundsEndX = m_abcCoordOrigin.x + m_abcCoordSize.xsize();
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int y = b + m_abcCoordOrigin.y, z = c + m_abcCoordOrigin.z;
            int bcIndex = b + c * bsize;
            int runRangeStart = m_bcToRunIndex[bcIndex], runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;

            if( runRangeStart + 1 != runRangeEnd || m_runData[runRangeStart].x != m_runData[runRangeStart + 1].x ) {
                // First go through all the X's before the first run
                int runEndX = m_runData[runRangeStart].x;
                for( int x = boundsStartX; x < runEndX; ++x ) {
                    int dataIndex = XYZtoDataIndex( vector3( x, y, z ) );
                    if( dataIndex != m_exteriorRegionCode ) {
                        out << "XYZ coordinate (" << x << ", " << y << ", " << z << "), for BC coordinate (" << b
                            << ", " << c << ") had an incorrect lookup." << endl;
                        out << "The data index received was " << dataIndex
                            << ", but it should match the exterior region code " << m_exteriorRegionCode << endl;
                        return false;
                    }
                }
                // Go through all the runs
                for( int run = runRangeStart; run != runRangeEnd; ++run ) {
                    int runStartX = m_runData[run].x;
                    runEndX = m_runData[run + 1].x;
                    int runDataIndex = m_runData[run].dataIndex;
                    if( runDataIndex >= 0 ) {
                        for( int x = runStartX; x < runEndX; ++x ) {
                            int dataIndex = XYZtoDataIndex( vector3( x, y, z ) );
                            if( dataIndex != runDataIndex + x - runStartX ) {
                                out << "XYZ coordinate (" << x << ", " << y << ", " << z << "), for BC coordinate ("
                                    << b << ", " << c << ") had an incorrect lookup." << endl;
                                out << "The data index received was " << dataIndex
                                    << ", but it should match the data index " << ( runDataIndex + x - runStartX )
                                    << endl;
                                return false;
                            }
                        }
                    } else {
                        for( int x = runStartX; x < runEndX; ++x ) {
                            int dataIndex = XYZtoDataIndex( vector3( x, y, z ) );
                            if( dataIndex != runDataIndex ) {
                                out << "XYZ coordinate (" << x << ", " << y << ", " << z << "), for BC coordinate ("
                                    << b << ", " << c << ") had an incorrect lookup." << endl;
                                out << "The data index received was " << dataIndex
                                    << ", but it should match the region code " << runDataIndex << endl;
                                return false;
                            }
                        }
                    }
                }
                // Finally go through all the X's after the last run
                for( int x = m_runData[runRangeEnd].x; x < boundsEndX; ++x ) {
                    int dataIndex = XYZtoDataIndex( vector3( x, y, z ) );
                    if( dataIndex != m_exteriorRegionCode ) {
                        out << "XYZ coordinate (" << x << ", " << y << ", " << z << "), for BC coordinate (" << b
                            << ", " << c << ") had an incorrect lookup." << endl;
                        out << "The data index received was " << dataIndex
                            << ", but it should match the exterior region code " << m_exteriorRegionCode << endl;
                        return false;
                    }
                }
            }
            // It's a single empty run, so everything in this scanline is the exterior region code
            else {
                for( int x = boundsStartX; x < boundsEndX; ++x ) {
                    int dataIndex = XYZtoDataIndex( vector3( x, y, z ) );
                    if( dataIndex != m_exteriorRegionCode ) {
                        out << "XYZ coordinate (" << x << ", " << y << ", " << z << "), for BC coordinate (" << b
                            << ", " << c << ") had an incorrect lookup." << endl;
                        out << "The data index received was " << dataIndex
                            << ", but it should match the exterior region code " << m_exteriorRegionCode << endl;
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

const ris_adjacency& rle_index_spec::get_cached_adjacency() const {
    if( m_cachedRisAdjacency ) {
        return *m_cachedRisAdjacency;
    } else {
        m_cachedRisAdjacency = new ris_adjacency();
        m_cachedRisAdjacency->compute( *this );
        return *m_cachedRisAdjacency;
    }
}

void rle_index_spec::free_cached_adjacency() const {
    if( m_cachedRisAdjacency ) {
        delete m_cachedRisAdjacency;
        m_cachedRisAdjacency = 0;
    }
}

boundbox3 rle_index_spec::compute_defined_bounds() const {
    // Defaults to empty
    frantic::graphics::boundbox3 result;

    // Go through all the runs, and include the first and last voxel of each defined run
    boundbox3 outerBounds = outer_bounds();
    for( int z = outerBounds.zminimum(); z <= outerBounds.zmaximum(); ++z ) {
        for( int y = outerBounds.yminimum(); y <= outerBounds.ymaximum(); ++y ) {
            // Iterate through all the runs in this scanline
            for( rle_run_iterator i( *this, y_to_b( y ), z_to_c( z ) ), ie; i != ie; ++i ) {
                if( i.get_data_index() >= 0 ) {
                    result.include_x_interval( i.get_xmin(), i.get_xmax() );
                    result.include_y( y );
                    result.include_z( z );
                }
            }
        }
    }

    return result;
}

// A function object which orders vector3 values in the rle_index_spec sort order (column-major).
// NOTE: This is the same sorted order as in rle_index_spec.cpp rle_index_spec_sort_order_pair class.
class rle_index_spec_sort_order {
  public:
    bool operator()( const vector3& lhs, const vector3& rhs ) const {
        if( lhs.z != rhs.z )
            return lhs.z < rhs.z;
        else if( lhs.y != rhs.y )
            return lhs.y < rhs.y;
        else
            return lhs.x < rhs.x;
    }
};

// This function build an rle_index_spec from an array of defined voxel coordinates.  It may change the order of the
// voxelArray during processing.
void rle_index_spec::build_from_voxel_array( std::vector<vector3>& coordinateArray ) {
    // First sort the coordinate array so that it goes in the order the runs do in the data array.  After this sorting,
    // the coordinate array exactly matches what the RLE mapping from coordinates to data indices will be.
    std::sort( coordinateArray.begin(), coordinateArray.end(), rle_index_spec_sort_order() );

    //	cout << endl;
    //	for( unsigned i = 0; i < coordinateArray.size(); ++i ) {
    //		cout << i << " : " << coordinateArray[i] << endl;
    //	}
    //	cout << endl;

    // Use -1 ("outside") as the default exterior region.
    m_exteriorRegionCode = -1;

    // Get the outer bounds of the voxels
    boundbox3 coordinateBounds;
    for( std::vector<vector3>::const_iterator i = coordinateArray.begin(); i != coordinateArray.end(); ++i )
        coordinateBounds += *i;

    //	cout << "coordinateBounds: " << coordinateBounds << endl;

    m_abcCoordOrigin = coordinateBounds.minimum();
    m_abcCoordSize = coordinateBounds.size();

    //	cout << "m_abcCoordOrigin: " << m_abcCoordOrigin << endl;
    //	cout << "m_abcCoordSize: " << m_abcCoordSize << endl;

    // Initialize the vector<> classes
    // The size of the bcToRunIndex array needs to be one bigger than the 2D array, so that we can always access one
    // value greater to get the scanline range.
    m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
    m_runData.clear();

    int currentDataIndex = 0, dataIndexSize = (int)coordinateArray.size();
    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;

            int y = b + m_abcCoordOrigin.y, z = c + m_abcCoordOrigin.z;

            //			cout << "Checking yz coordinate " << y << ", " << z << endl;

            //			cout << "Its bcIndex is " << bcIndex << " (of " << m_bcToRunIndex.size() << ")" << endl;
            m_bcToRunIndex[bcIndex] = (int)m_runData.size();

            //			if( currentDataIndex < dataIndexSize ) cout << "Next voxel coordinate (index " << currentDataIndex
            //<< ") is " << coordinateArray[currentDataIndex] << endl;

            // Check that the current voxel coordinate has matching y,z coordinates
            if( currentDataIndex < dataIndexSize && coordinateArray[currentDataIndex].y == y &&
                coordinateArray[currentDataIndex].z == z ) {
                int lastXValue = coordinateArray[currentDataIndex].x;

                // Start a run of defined voxels pointing to the current data index
                m_runData.push_back( run_data( lastXValue, currentDataIndex ) );
                ++currentDataIndex;
                // Loop through all the voxels with matching y,z coordinates
                while( currentDataIndex < dataIndexSize && coordinateArray[currentDataIndex].y == y &&
                       coordinateArray[currentDataIndex].z == z ) {
                    int currentXValue = coordinateArray[currentDataIndex].x;

                    //					cout << " lastXValue " << lastXValue << ", currentXValue " << currentXValue <<
                    //endl;

                    // If the current X value is not adjacent to the last one, then we have to insert an undefined run,
                    // otherwise keep collecting defined values
                    if( currentXValue > lastXValue + 1 ) {
                        // Start an undefined run after the last seen X value
                        m_runData.push_back( run_data( lastXValue + 1, -1 ) );
                        // Start a defined run at the currently seen X value
                        m_runData.push_back( run_data( currentXValue, currentDataIndex ) );
                    } else if( currentXValue != lastXValue + 1 ) {
                        throw std::runtime_error(
                            "rle_index_spec.build_from_voxel_array: An unexpected error occurred, the X "
                            "values in the sorted coordinate array were not ascending." );
                    }
                    ++currentDataIndex;
                    lastXValue = currentXValue;
                }
                // Finish off the last run by providing its ending value
                m_runData.push_back( run_data( lastXValue + 1, -1 ) );
            } else {
                // Create a zero-sized run, necessary for the coordinate queries to work properly
                m_runData.push_back( run_data( 0, -1 ) );
                m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }
    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();

    if( currentDataIndex != dataIndexSize )
        throw std::runtime_error( "rle_index_spec.build_from_voxel_array: An unexpected error occurred, the creation "
                                  "loop failed to use all " +
                                  boost::lexical_cast<std::string>( dataIndexSize ) + " voxels, and instead used " +
                                  boost::lexical_cast<std::string>( currentDataIndex ) );

    m_dataSize = coordinateArray.size();
}

void rle_index_spec::build_from_random( const boundbox3& testBounds, float expectedRunLength,
                                        int randomRegionCodeLimit ) {
    if( randomRegionCodeLimit < 2 )
        throw runtime_error(
            "rle_index_spec::build_from_random() - The random region code limit must be 2 or greater, it was set to " +
            lexical_cast<string>( randomRegionCodeLimit ) + "." );

    float runChangeProbability = 1.f / expectedRunLength;

    // Make a randomized outer_bounds
    boundbox3 ob( testBounds.random_vector() );
    ob += testBounds.random_vector();
    while( ob.get_volume() <= 1 )
        ob += testBounds.random_vector();

    // Figure out which voxels we want to include
    vector<vector3> voxels;
    bool isDefined = false;
    for( int z = ob.zminimum(); z <= ob.zmaximum(); ++z ) {
        for( int y = ob.yminimum(); y <= ob.ymaximum(); ++y ) {
            for( int x = ob.xminimum(); x <= ob.xmaximum(); ++x ) {
                // switch the isDefined flag with a specified probability
                float chooser = (float)rand() / RAND_MAX;
                if( chooser < runChangeProbability )
                    isDefined = !isDefined;
                if( isDefined )
                    voxels.push_back( vector3( x, y, z ) );
            }
        }
    }

    // Build the rle index spec
    build_from_voxel_array( voxels );

    // Randomize the region codes of the undefined runs
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = m_bcToRunIndex[bcIndex], runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;

            run_data *rd = &m_runData[runRangeStart], *rdEnd = &m_runData[runRangeEnd];
            for( ; rd != rdEnd; ++rd ) {
                if( rd->dataIndex < 0 )
                    rd->dataIndex = -( ( rand() % randomRegionCodeLimit ) + 1 );
            }
        }
    }

    // Randomize the exterior region code
    m_exteriorRegionCode = -( ( rand() % randomRegionCodeLimit ) + 1 );
}

void rle_index_spec::apply_axis_permutation( const vector3& axisPermutation,
                                             const rle_index_spec_channel_copying_data* channelsToRemapBegin,
                                             const rle_index_spec_channel_copying_data* channelsToRemapEnd ) {
    // Only apply the permutation if it's not the identity.  If there are channels to remap, we always apply the
    // permutation, so the caller should check for the identity if it is passing in channels for remapping purposes.
    if( channelsToRemapBegin != channelsToRemapEnd || axisPermutation.x != 0 || axisPermutation.y != 1 ||
        axisPermutation.z != 2 ) {
        rle_index_spec permuteRis;
        permuteRis.build_from_axis_permutation( axisPermutation, *this, channelsToRemapBegin, channelsToRemapEnd );
        swap( permuteRis );
    }
}

void rle_index_spec::build_from_dilation( const rle_index_spec& ris, int dilationVoxels ) {
    if( dilationVoxels < 1 )
        throw std::runtime_error(
            "rle_index_spec.build_from_dilation() - Tried to dilate the voxel region with a value less than one, " +
            lexical_cast<string>( dilationVoxels ) + "." );

    // Cycle through the axes, dilating along X each time.
    rle_index_spec tempRisA, tempRisB;
    tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), ris );
    tempRisB.build_from_x_dilation( tempRisA, dilationVoxels, dilationVoxels );
    tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), tempRisB );
    tempRisB.build_from_x_dilation( tempRisA, dilationVoxels, dilationVoxels );
    tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), tempRisB );
    build_from_x_dilation( tempRisA, dilationVoxels, dilationVoxels );
}

void rle_index_spec::build_from_dilation( const rle_index_spec& ris, const boundbox3& dilationBox ) {
    if( !dilationBox.contains( vector3( 0 ) ) )
        throw std::runtime_error(
            "rle_index_spec.build_from_dilation() - Tried to dilate the voxel region with an invalid dilation box " +
            dilationBox.str() + ", that did not contain the origin." );
    if( dilationBox.maximum() - dilationBox.minimum() == vector3( 0 ) )
        throw std::runtime_error(
            "rle_index_spec.build_from_dilation() - Tried to dilate the voxel region with a dilation box " +
            dilationBox.str() + ", that specifies no dilation." );

    // If there is only dilation along X, then don't do any axis permutations
    if( dilationBox.yminimum() == dilationBox.ymaximum() && dilationBox.zminimum() == dilationBox.zmaximum() ) {
        build_from_x_dilation( ris, -dilationBox.xminimum(), dilationBox.xmaximum() );
    }
    // If there is no dilation along Y, then skip one axis permutation/dilation
    else if( dilationBox.yminimum() == dilationBox.ymaximum() ) {
        rle_index_spec tempRisA, tempRisB;
        tempRisA.build_from_axis_permutation( vector3( 2, 0, 1 ), ris );
        tempRisB.build_from_x_dilation( tempRisA, -dilationBox.zminimum(), dilationBox.zmaximum() );
        // Avoid the X dilation if possible
        if( dilationBox.xminimum() == dilationBox.xmaximum() ) {
            build_from_axis_permutation( vector3( 1, 2, 0 ), tempRisB );
        } else {
            tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), tempRisB );
            build_from_x_dilation( tempRisA, -dilationBox.xminimum(), dilationBox.xmaximum() );
        }
    }
    // If there is no dilation along Z, then skip one axis permutation/dilation
    else if( dilationBox.zminimum() == dilationBox.zmaximum() ) {
        // Cycle through the axes, dilating along X each time.
        rle_index_spec tempRisA, tempRisB;
        tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), ris );
        tempRisB.build_from_x_dilation( tempRisA, -dilationBox.yminimum(), dilationBox.ymaximum() );
        // Avoid the X dilation if possible
        if( dilationBox.xminimum() == dilationBox.xmaximum() ) {
            build_from_axis_permutation( vector3( 2, 0, 1 ), tempRisB );
        } else {
            tempRisA.build_from_axis_permutation( vector3( 2, 0, 1 ), tempRisB );
            build_from_x_dilation( tempRisA, -dilationBox.xminimum(), dilationBox.xmaximum() );
        }
    }
    // There is dilation along both Y and Z
    else {
        // Cycle through the axes, dilating along X each time.
        rle_index_spec tempRisA, tempRisB;
        tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), ris );
        tempRisB.build_from_x_dilation( tempRisA, -dilationBox.yminimum(), dilationBox.ymaximum() );
        tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), tempRisB );
        tempRisB.build_from_x_dilation( tempRisA, -dilationBox.zminimum(), dilationBox.zmaximum() );
        // Avoid the X dilation if possible
        if( dilationBox.xminimum() == dilationBox.xmaximum() ) {
            build_from_axis_permutation( vector3( 1, 2, 0 ), tempRisB );
        } else {
            tempRisA.build_from_axis_permutation( vector3( 1, 2, 0 ), tempRisB );
            build_from_x_dilation( tempRisA, -dilationBox.xminimum(), dilationBox.xmaximum() );
        }
    }
}

void rle_index_spec::build_from_x_dilation( const rle_index_spec& ris, int dilationVoxelsNeg, int dilationVoxelsPos ) {
    free_cached_adjacency();

    m_exteriorRegionCode = ris.m_exteriorRegionCode;
    m_abcCoordOrigin = ris.m_abcCoordOrigin;
    m_abcCoordSize = ris.m_abcCoordSize;
    m_dataSize = 0;

    // Expand the size by dilation voxels along X
    m_abcCoordOrigin.x -= dilationVoxelsNeg;
    m_abcCoordSize.xsize() += dilationVoxelsNeg + dilationVoxelsPos;

    // clear any data already stored in ris
    m_runData.clear();

    // Resize to match the correct size.  Note that we're building a new rle_index_spec
    // from scrach, not modifying the input one in place, therefore this operation is necessary.
    m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );

    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = ris.m_bcToRunIndex[bcIndex], runRangeEnd = ris.m_bcToRunIndex[bcIndex + 1] - 1;
            const run_data *rdStart = &ris.m_runData[runRangeStart], *rdEnd = &ris.m_runData[runRangeEnd];

            m_bcToRunIndex[bcIndex] = (int)m_runData.size();

            // Only process if this is a non-empty scanline.
            if( rdStart->x != ( rdStart + 1 )->x ) {
                bool creatingUndefinedRun = true;

                const run_data *rdPrevDefined = rdStart - 1, *rdCurDefined = rdStart;
                // Increment through the runs to find the next defined one
                while( rdCurDefined != rdEnd && rdCurDefined->dataIndex < 0 )
                    ++rdCurDefined;

                // forwardX is the X value at the last defined voxel in the previous defined run
                // backwardX is the X value at the first defined voxel in the current defined run
                // forwardDilateX is the largest X value contained within the forward dilation region
                // backwardDilateX is the smallest X value contained within the backward dilation region
                int32_t forwardX, forwardDilateX, backwardX, backwardDilateX;

                int scanlineEndX = rdEnd->x;

                // Keep looping until the previous defined run passes the end.  This way we will iterate once
                // for every gap between two defined runs as well as before the first defined run and after
                // the last defined run.
                while( rdPrevDefined != rdEnd ) {

                    // Get the X and dilated X from the end of the previous defined run
                    if( rdPrevDefined >= rdStart ) {
                        forwardX = ( rdPrevDefined + 1 )->x - 1;
                        forwardDilateX = forwardX + dilationVoxelsPos;
                    } else {
                        // This happens in the first iteration, making a virtual previous defined voxel
                        // before the outer bounds start.
                        forwardX = m_abcCoordOrigin.x - 1;
                        forwardDilateX = forwardX;
                    }

                    // Get the X and dilated X from the start of the current defined run
                    if( rdCurDefined < rdEnd ) {
                        backwardX = rdCurDefined->x;
                        backwardDilateX = backwardX - dilationVoxelsNeg;
                    } else {
                        // This happens in the last iteration, making a virtual next defined voxel
                        // afterthe outer bounds end.
                        backwardX = m_abcCoordOrigin.x + m_abcCoordSize.xsize();
                        backwardDilateX = backwardX;
                    }

                    // Clip the dilation extents so that they stop at the next defined run or at the bounds.
                    if( forwardDilateX >= backwardX )
                        forwardDilateX = backwardX - 1;
                    if( backwardDilateX <= forwardX )
                        backwardDilateX = forwardX + 1;

                    // If the dilation regions overlap, then we have to clip them so they don't.
                    if( forwardDilateX >= backwardDilateX ) {
                        backwardDilateX = forwardDilateX + 1;
                    }

                    //////////
                    // First do the forward dilation
                    //////////

                    if( forwardDilateX > forwardX ) {
                        // If the forward dilation has extended past the scanline end, increase the scanline end.
                        if( forwardDilateX >= scanlineEndX )
                            scanlineEndX = forwardDilateX + 1;

                        m_dataSize += forwardDilateX - forwardX;
                    }

                    //////////
                    // Second, copy the undefined runs in between
                    //////////

                    // Here we fill in the voxels that are in the open interval (forwardDilateX, backwardDilateX) using
                    // undefined runs
                    if( forwardDilateX + 1 < backwardDilateX ) {
                        // Find the first undefined run that intersects with this interval
                        const run_data* rd = rdPrevDefined + 1;
                        while( rd < rdCurDefined && ( rd + 1 )->x <= forwardDilateX + 1 )
                            ++rd;
                        // Add this run
                        if( rd < rdCurDefined ) {
                            int x;
                            // If this is the first run, is undefined, and started at the beginning of the
                            // outer bounds, then we extend this run backwards to the beginning of the new
                            // outer bounds.  Otherwise, the run starts either where it started before or
                            // just after the defined dilation end.
                            if( rd == rdStart && rd->x == ris.m_abcCoordOrigin.x )
                                x = m_abcCoordOrigin.x;
                            else
                                x = max( forwardDilateX + 1, rd->x );
                            if( x < backwardDilateX ) {
                                m_runData.push_back( run_data( x, rd->dataIndex ) );
                                creatingUndefinedRun = true;
                            }
                        }
                        ++rd;
                        // Copy the rest of the runs until they don't intersect with this interval anymore
                        while( rd < rdCurDefined && rd->x <= backwardDilateX - 1 ) {
                            m_runData.push_back( *rd );
                            creatingUndefinedRun = true;
                            ++rd;
                        }
                    }

                    // Start the defined run if necessary (i.e. an undefined run is being created and we haven't passed
                    // all the defined runs)
                    if( creatingUndefinedRun && rdCurDefined < rdEnd ) {
                        m_runData.push_back( run_data( backwardDilateX, (int)m_dataSize ) );
                        creatingUndefinedRun = false;
                    }

                    //////////
                    // Third, do the backward dilation
                    //////////

                    if( backwardDilateX < backwardX ) {
                        m_dataSize += backwardX - backwardDilateX;
                    }

                    //////////
                    // Fourth, copy the defined run
                    //////////

                    if( rdCurDefined < rdEnd ) {
                        m_dataSize += ( rdCurDefined + 1 )->x - backwardX;
                    }

                    // Increment to the next gap between defined runs
                    rdPrevDefined = rdCurDefined;
                    ++rdCurDefined;
                    // Increment through the runs to find the next defined one
                    while( rdCurDefined < rdEnd && rdCurDefined->dataIndex < 0 )
                        ++rdCurDefined;
                }

                // If the last run is an undefined run, and ended at the outer bounds, then extend
                // the last run so that it goes all the way to the end of the new outer bounds.
                if( ( rdEnd - 1 )->dataIndex < 0 && rdEnd->x == ris.m_abcCoordOrigin.x + ris.m_abcCoordSize.xsize() ) {
                    int newScanlineEndX = m_abcCoordOrigin.x + m_abcCoordSize.xsize();
                    // Only do this if the current scanline end doesn't stop at the outer bounds
                    if( scanlineEndX < newScanlineEndX ) {
                        // If the dilation took the defined run past this undefined run, we have to add it again.
                        if( m_runData.back().dataIndex != ( rdEnd - 1 )->dataIndex )
                            m_runData.push_back( run_data( scanlineEndX, ( rdEnd - 1 )->dataIndex ) );
                        scanlineEndX = newScanlineEndX;
                    }
                }

                // Finish off the last run of the scanline
                m_runData.push_back( run_data( scanlineEndX, -1 ) );
            } else {
                // Create an empty run for this BC coordinate
                m_runData.push_back( run_data( 0, -1 ) );
                m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();

    // run a consistency check
    stringstream sout;
    if( !ris.check_consistency( sout ) ) {
        throw std::runtime_error( "build_from_x_dilation: RLE Index Spec consistency check failed:\n" + sout.str() );
    }
}

void rle_index_spec::copy_data_channels( const rle_index_spec& outputRIS, const rle_index_spec& inputRIS,
                                         const rle_index_spec_channel_copying_data* channelsToCopyBegin,
                                         const rle_index_spec_channel_copying_data* channelsToCopyEnd ) {
    for( int c = 0, csize = outputRIS.m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = outputRIS.m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = outputRIS.m_bcToRunIndex[bcIndex],
                runRangeEnd = outputRIS.m_bcToRunIndex[bcIndex + 1] - 1;
            const run_data *rdOutput = &outputRIS.m_runData[runRangeStart],
                           *rdOutputEnd = &outputRIS.m_runData[runRangeEnd];

            // Check that the output scanline is not empty
            if( rdOutput->x != ( rdOutput + 1 )->x ) {
                // Get the bcIndex for the input rle index spec
                int bInput = b + outputRIS.m_abcCoordOrigin.y - inputRIS.m_abcCoordOrigin.y;
                int cInput = c + outputRIS.m_abcCoordOrigin.z - inputRIS.m_abcCoordOrigin.z;

                if( (unsigned)bInput < (unsigned)inputRIS.m_abcCoordSize.ysize() &&
                    (unsigned)cInput < (unsigned)inputRIS.m_abcCoordSize.zsize() ) {
                    int bcIndexInput = bInput + cInput * inputRIS.m_abcCoordSize.ysize();
                    int runRangeStartInput = inputRIS.m_bcToRunIndex[bcIndexInput],
                        runRangeEndInput = inputRIS.m_bcToRunIndex[bcIndexInput + 1] - 1;
                    const run_data *rdInput = &inputRIS.m_runData[runRangeStartInput],
                                   *rdInputEnd = &inputRIS.m_runData[runRangeEndInput];

                    // Only if this is a non-empty scanline is there anything interesting to do
                    int xBeginInputScanline = rdInput->x, xEndInputScanline = rdInputEnd->x;
                    if( xBeginInputScanline != xEndInputScanline ) {
                        // Before we can go through the runs in lockstep, we have to go through the runs which intersect
                        // with parts outside of the input scanline.
                        while( rdOutput != rdOutputEnd && ( rdOutput + 1 )->x < xBeginInputScanline )
                            ++rdOutput;

                        // Go through the runs in lockstep, copying values as necessary
                        while( rdOutput != rdOutputEnd ) {
                            int dataIndexOutput = rdOutput->dataIndex;
                            // If this is a defined run, we may input some values in the output channel
                            if( dataIndexOutput >= 0 ) {
                                // Get the half-open interval [xBegin,xEnd) over which this run is defined
                                int xBeginOutput = rdOutput->x, xEndOutput = ( rdOutput + 1 )->x;
                                // Increment the rdInput variable until it is at a run which ends at least where the
                                // interval [xBegin,xEnd) starts.
                                while( rdInput != rdInputEnd && ( rdInput + 1 )->x <= xBeginOutput )
                                    ++rdInput;
                                // Process all the intervals which intersect the interval [xBegin,xEnd).  That is, all
                                // the intervals which start before the inputination intervals end
                                while( rdInput != rdInputEnd ) {
                                    int dataIndexInput = rdInput->dataIndex;
                                    // If this is an defined run then there may be something to copy
                                    if( dataIndexInput >= 0 ) {
                                        // Get the half-open interval [xBeginInput,xEndInput) over which this input run
                                        // is defined
                                        int xBeginInput = rdInput->x, xEndInput = ( rdInput + 1 )->x;
                                        int xBeginIntersected = ( std::max )( xBeginOutput, xBeginInput );
                                        // Get the data index range
                                        int dataIndexCopy = dataIndexOutput + xBeginIntersected - xBeginOutput;
                                        int dataIndexCopySize =
                                            ( std::min )( xEndOutput, xEndInput ) - xBeginIntersected;

                                        // Offset the input data index based on the sub-interval being processed
                                        dataIndexInput += xBeginIntersected - xBeginInput;

                                        // Copy the data for each of the data channels
                                        for( const rle_index_spec_channel_copying_data* i = channelsToCopyBegin;
                                             i != channelsToCopyEnd; ++i )
                                            memcpy( i->outputData + i->primitiveSize * dataIndexCopy,
                                                    i->inputData + i->primitiveSize * dataIndexInput,
                                                    i->primitiveSize * dataIndexCopySize );
                                    }
                                    if( ( rdInput + 1 )->x < xEndOutput )
                                        ++rdInput;
                                    else
                                        break;
                                }
                            }
                            // If the inputination runs have run out, part of the current rd run we're processing may
                            // be outside the [xBeginInputScanline,xEndInputScanline) interval, so we don't increment
                            // rd in that case.
                            if( rdInput != rdInputEnd )
                                ++rdOutput;
                            else
                                break;
                        }
                    }
                }
            }
        }
    }
}

// This function builds an RLE index spec from a legacy flood rle region.
// An axis permutation will have to be applied after the region is created.  The reason a permutation is returned
// instead of being applied directly is so the caller can also read in data to remap with the axis permutation.
void rle_index_spec::build_from_legacy_flood_rle_region( int compressionAxis, const boundbox3& voxelBounds,
                                                         const std::vector<boost::int32_t>& bcOffsets,
                                                         const std::vector<boost::int16_t>& runBoundaries,
                                                         const std::vector<boost::int32_t>& runDataIndices,
                                                         vector3& outRequiredAxisPermutation ) {
    // If the legacy region was compressed on a different axis, then we need to permute the axes so it's compressed
    // along the X axis We can't do this now, because permuting the axis will reorder the data.  The caller of this
    // function needs to do the permutation once it also has the data.
    if( compressionAxis == 1 ) {
        outRequiredAxisPermutation = vector3( 1, 2, 0 );

        // Set up the bounds, applying the inverse permutation
        vector3 v = voxelBounds.minimum();
        size3 s = voxelBounds.size();
        m_abcCoordOrigin = vector3( v.z, v.x, v.y );
        m_abcCoordSize = size3( s.zsize(), s.xsize(), s.ysize() );
    } else if( compressionAxis == 2 ) {
        outRequiredAxisPermutation = vector3( 2, 1, 0 );

        // Set up the bounds, applying the inverse permutation
        vector3 v = voxelBounds.minimum();
        size3 s = voxelBounds.size();
        m_abcCoordOrigin = vector3( v.z, v.y, v.x );
        m_abcCoordSize = size3( s.zsize(), s.ysize(), s.xsize() );
    } else {
        outRequiredAxisPermutation = vector3( 0, 1, 2 );

        // Set up the bounds
        m_abcCoordOrigin = voxelBounds.minimum();
        m_abcCoordSize = voxelBounds.size();
    }

    // Use -1 ("outside") as the default exterior region.  This property didn't exist in the legacy flood structures, so
    // we set it arbitrarily.
    m_exteriorRegionCode = -1;

    if( (int)bcOffsets.size() != m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 )
        throw runtime_error(
            "rle_index_spec.build_from_legacy_flood_rle_region: The bcOffsets array is the incorrect size.  It is " +
            lexical_cast<string>( bcOffsets.size() ) + ", but should be " +
            lexical_cast<string>( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 ) + "." );

    if( runBoundaries.size() != runDataIndices.size() )
        throw runtime_error(
            "rle_index_spec.build_from_legacy_flood_rle_region: The runBoundaries and runDataIndices were "
            "different sizes, but should have been the same." );

    // The run structure is the same as in the C# class, so we could use it directly.  In legacy Flood, every voxel
    // within the region had to be contained within a run, whereas in the new structure we try to trim away the extra
    // runs at the front and the back.  Thus, we reprocess the input we get to remove the redundant runs.
    m_bcToRunIndex.resize( bcOffsets.size() );
    m_runData.clear();
    m_runData.reserve( runBoundaries.size() );

    // Go through all the scanlines.  The data index values will stay the same, so we don't need to remap them, but we
    // want to remove all the redundant runs.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = bcOffsets[bcIndex], runRangeEnd = bcOffsets[bcIndex + 1] - 2;

            m_bcToRunIndex[bcIndex] = (int)m_runData.size();

            // Special case the single run scanline
            if( runRangeStart >= runRangeEnd ) {
                int dataIndex = runDataIndices[runRangeStart];
                // Remap the region codes from the legacy values to the values we use now.
                if( dataIndex < 0 )
                    dataIndex = ( dataIndex == -2000 ) ? -2 : -1;
                if( dataIndex == m_exteriorRegionCode ) {
                    // Add an empty scanline consisting of one zero-sized run
                    m_runData.push_back( run_data( 0, -1 ) );
                    m_runData.push_back( run_data( 0, -1 ) );
                } else {
                    // Copy the single run
                    m_runData.push_back( run_data( runBoundaries[runRangeStart], dataIndex ) );
                    m_runData.push_back( run_data( runBoundaries[runRangeStart + 1], -1 ) );
                }
            }
            // Deal with multiple runs in the scanline
            else {
                // Copy the first run start if it's different from the exterior region code
                int dataIndex = runDataIndices[runRangeStart];
                // Remap the region codes from the legacy values to the values we use now.
                if( dataIndex < 0 )
                    dataIndex = ( dataIndex == -2000 ) ? -2 : -1;
                if( dataIndex != m_exteriorRegionCode )
                    m_runData.push_back( run_data( runBoundaries[runRangeStart], dataIndex ) );

                // Copy all the middle run values
                for( int run = runRangeStart + 1; run <= runRangeEnd - 1; ++run ) {
                    dataIndex = runDataIndices[run];
                    // Remap the region codes from the legacy values to the values we use now.
                    if( dataIndex < 0 )
                        dataIndex = ( dataIndex == -2000 ) ? -2 : -1;
                    m_runData.push_back( run_data( runBoundaries[run], dataIndex ) );
                }

                // Copy the last run if it's different from the exterior region code.
                dataIndex = runDataIndices[runRangeEnd];
                // Remap the region codes from the legacy values to the values we use now.
                if( dataIndex < 0 )
                    dataIndex = ( dataIndex == -2000 ) ? -2 : -1;
                // If we haven't added any runs yet, we have to deal with the potential for this turning into an empty
                // scanline, so check for that first.
                if( m_bcToRunIndex[bcIndex] == (int)m_runData.size() ) {
                    if( dataIndex != m_exteriorRegionCode ) {
                        // Copy the last run, capping off the scanline
                        m_runData.push_back( run_data( runBoundaries[runRangeEnd], dataIndex ) );
                        m_runData.push_back( run_data( runBoundaries[runRangeEnd + 1], -1 ) );
                    } else {
                        // Add an empty scanline consisting of one zero-sized run
                        m_runData.push_back( run_data( 0, -1 ) );
                        m_runData.push_back( run_data( 0, -1 ) );
                    }
                } else {
                    if( dataIndex != m_exteriorRegionCode ) {
                        // Copy the last run, capping off the scanline
                        m_runData.push_back( run_data( runBoundaries[runRangeEnd], dataIndex ) );
                        m_runData.push_back( run_data( runBoundaries[runRangeEnd + 1], -1 ) );
                    } else {
                        // Cap off the last run that was started
                        m_runData.push_back( run_data( runBoundaries[runRangeEnd], -1 ) );
                    }
                }
            }
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();

    // Loop backwards through the runs to find the last defined run, and get the data size from that
    bool foundIt = false;
    m_dataSize = 0;
    for( int c = m_abcCoordSize.zsize() - 1; c >= 0 && !foundIt; --c ) {
        for( int b = m_abcCoordSize.ysize() - 1; b >= 0 && !foundIt; --b ) {
            int bcIndex = b + c * m_abcCoordSize.ysize();
            int runRangeStart = m_bcToRunIndex[bcIndex], runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 2;

            for( int run = runRangeEnd; run >= runRangeStart && !foundIt; --run ) {
                // Check that the run value is within bounds, because this data is from an untrusted source
                if( m_runData[run].dataIndex >= 0 ) {
                    m_dataSize = m_runData[run].dataIndex + m_runData[run + 1].x - m_runData[run].x;
                    foundIt = true;
                }
            }
        }
    }

    // Check the consistency of the data structure, because you never know what might have been in this file...
    stringstream errorMessage;
    if( !check_consistency( errorMessage ) )
        throw std::runtime_error( "rle_index_spec.build_from_legacy_flood_rle_region: There was a consistency check "
                                  "failure with the data.\n" +
                                  errorMessage.str() );
}

void rle_index_spec::fill_2x2x2_data_index_box( const vector3& minCornerXYZ, int32_t* outDataIndices ) const {
    // Compute the B and C coordinates from the xyz
    int b = minCornerXYZ.y - m_abcCoordOrigin.y;
    int c = minCornerXYZ.z - m_abcCoordOrigin.z;
    int x = minCornerXYZ.x;
    // Only test the B and C coordinate values initially.  The X coordinate value will be tested after the binary search
    // is done, based on the assumption that the expected case will often be to hit defined runs.
    if( b >= 0 && b + 1 < m_abcCoordSize.ysize() && c >= 0 && c + 1 < m_abcCoordSize.zsize() ) {
        int runIndex;
        int bcIndex = b + m_abcCoordSize.ysize() * c;
        // Loop through the run index values and get the two data index values from each
        for( int cOffset = 0; cOffset < 2; ++cOffset, bcIndex += m_abcCoordSize.ysize() ) {
            for( int bOffset = 0; bOffset < 2; ++bOffset ) {
                runIndex = BCIndexXtoRunIndex( bcIndex + bOffset, x );

                int runBeginX = m_runData[runIndex].x, runEndX = m_runData[runIndex + 1].x;
                // If both points are within this run, we can fill these two indices directly from this run's data
                if( runBeginX <= x && x + 1 < runEndX ) {
                    int runDataIndex = m_runData[runIndex].dataIndex;
                    if( runDataIndex >= 0 ) {
                        runDataIndex += x - runBeginX;
                        // Fill in the two adjacent data index values
                        *outDataIndices++ = runDataIndex;
                        *outDataIndices++ = runDataIndex + 1;
                    } else {
                        // Fill in the region code twice
                        *outDataIndices++ = runDataIndex;
                        *outDataIndices++ = runDataIndex;
                    }
                }
                // If both points are outside of this run, we can fill these two indices from the exterior region code
                // (See the OPTIMIZATION NOTE for the BCIndexXtoRunIndex function.)
                else if( runBeginX > x + 1 || x >= runEndX ) {
                    *outDataIndices++ = m_exteriorRegionCode;
                    *outDataIndices++ = m_exteriorRegionCode;
                }
                // If this pair of points is straddling the end of this run.
                else if( x + 1 >= runEndX ) {
                    // First value gets retrieved from this run
                    int runDataIndex = m_runData[runIndex].dataIndex;
                    if( runDataIndex >= 0 ) {
                        *outDataIndices++ = runDataIndex + x - runBeginX;
                    } else {
                        *outDataIndices++ = runDataIndex;
                    }
                    // Check whether the next run exists, and if so, get the next value from it.
                    if( runIndex <= m_bcToRunIndex[bcIndex + bOffset + 1] - 3 ) {
                        // Second value gets retrieved from the next run
                        int runDataIndex = m_runData[runIndex + 1].dataIndex;
                        if( runDataIndex >= 0 ) {
                            *outDataIndices++ = runDataIndex + x + 1 - runEndX;
                        } else {
                            *outDataIndices++ = runDataIndex;
                        }
                    } else {
                        // Second value gets the exterior region code
                        *outDataIndices++ = m_exteriorRegionCode;
                    }
                }
                // If this pair of points is straddling the beginning of this run.
                else {
                    // First value gets the exterior region code
                    *outDataIndices++ = m_exteriorRegionCode;
                    // Second value is retrieved from the run
                    int runDataIndex = m_runData[runIndex].dataIndex;
                    if( runDataIndex >= 0 ) {
                        *outDataIndices++ = runDataIndex + x + 1 - runBeginX;
                    } else {
                        *outDataIndices++ = runDataIndex;
                    }
                }
            }
        }
    } else if( b < -1 || c < -1 || b >= m_abcCoordSize.ysize() || c >= m_abcCoordSize.zsize() ) {
        // It's outside the bounds entirely, so everything is the exterior region code
        for( int i = 0; i < 8; ++i )
            outDataIndices[i] = m_exteriorRegionCode;
    } else {
        // It's on the boundary, so do it the easier slower way
        outDataIndices[0] = XYZtoDataIndex( minCornerXYZ );
        outDataIndices[1] = XYZtoDataIndex( vector3( minCornerXYZ.x + 1, minCornerXYZ.y, minCornerXYZ.z ) );
        outDataIndices[2] = XYZtoDataIndex( vector3( minCornerXYZ.x, minCornerXYZ.y + 1, minCornerXYZ.z ) );
        outDataIndices[3] = XYZtoDataIndex( vector3( minCornerXYZ.x + 1, minCornerXYZ.y + 1, minCornerXYZ.z ) );
        outDataIndices[4] = XYZtoDataIndex( vector3( minCornerXYZ.x, minCornerXYZ.y, minCornerXYZ.z + 1 ) );
        outDataIndices[5] = XYZtoDataIndex( vector3( minCornerXYZ.x + 1, minCornerXYZ.y, minCornerXYZ.z + 1 ) );
        outDataIndices[6] = XYZtoDataIndex( vector3( minCornerXYZ.x, minCornerXYZ.y + 1, minCornerXYZ.z + 1 ) );
        outDataIndices[7] = XYZtoDataIndex( vector3( minCornerXYZ.x + 1, minCornerXYZ.y + 1, minCornerXYZ.z + 1 ) );
    }
}

void rle_index_spec::fill_data_index_box( const boundbox3& voxelExtents, boost::int32_t* outDataIndices ) const {
    ////////////////
    // INITIALIZE
    ////////////////

    int exteriorRegionCode = m_exteriorRegionCode;
    int dataSize = voxelExtents.get_volume();

    // Get the B coordinate range and the C coordinate
    int bMin = voxelExtents.minimum().y - m_abcCoordOrigin.y, bMax = voxelExtents.maximum().y - m_abcCoordOrigin.y;
    int cMin = voxelExtents.minimum().z - m_abcCoordOrigin.z, cMax = voxelExtents.maximum().z - m_abcCoordOrigin.z;

    // If the C coordinate interval or the B coordinate interval is outside the range of the bounding box, set the whole
    // array to "outside"
    if( bMax < 0 || bMin >= m_abcCoordSize.ysize() || cMax < 0 || cMin >= m_abcCoordSize.zsize() ) {
        for( int i = 0; i < dataSize; ++i )
            *outDataIndices++ = exteriorRegionCode;
        return;
    }

    // These are the number of voxels to fill in before copying values within the defined B range of the bounding box
    int bMinExteriorCount = 0, bMaxExteriorCount = 0;

    // If the array desired extends in front of the defined B range, we will have to fill the first part with exterior
    // values.
    if( bMin < 0 ) {
        bMinExteriorCount = -bMin * voxelExtents.xsize();
        bMin = 0;
    }

    // If the array desired extends beyond the defined B range, we will have to fill the last part with exterior values
    if( bMax >= m_abcCoordSize.ysize() ) {
        bMaxExteriorCount = ( bMax - m_abcCoordSize.ysize() + 1 ) * voxelExtents.xsize();
        bMax = m_abcCoordSize.ysize() - 1;
    }

    ////////////////
    // COPY LEVEL SET VALUES
    ////////////////

    // This is the value we use to step through the outVoxelCornerValues
    int outputIndex = 0;

    // If the array desired extends in front of the defined C range, fill the first part with exterior values.
    if( cMin < 0 ) {
        int blankValueCount = -cMin * voxelExtents.xsize() * voxelExtents.ysize();
        for( int i = 0; i < blankValueCount; ++i )
            outDataIndices[outputIndex++] = exteriorRegionCode;
        cMin = 0;
    }

    // Trim the C values we iterate over to within the range
    if( cMax >= m_abcCoordSize.zsize() ) {
        cMax = m_abcCoordSize.zsize() - 1;
    }

    // Now iterate through all the C coordinates and copy the planes
    for( int c = cMin; c <= cMax; ++c ) {
        // Fill in the starting B out-of-bounds scanlines with exterior values
        for( int i = 0; i < bMinExteriorCount; ++i )
            outDataIndices[outputIndex++] = exteriorRegionCode;

        // Now iterate through all the B coordinates and copy the scanlines
        for( int b = bMin; b <= bMax; ++b ) {
            // Get the run range
            int bcIndex = b + c * m_abcCoordSize.ysize();
            int nextXToFill = voxelExtents.minimum().x, lastXToFill = voxelExtents.maximum().x;
            int runRangeStart = BCIndexXtoRunIndex( bcIndex, nextXToFill ),
                runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;

            int xStart, xPastTheEnd = m_runData[runRangeStart].x;
            for( int run = runRangeStart; run != runRangeEnd; ++run ) {
                // Get the x range of this run
                xStart = xPastTheEnd;
                xPastTheEnd = m_runData[run + 1].x;

                // If there are values before xStart, fill them up with "outside" values
                while( nextXToFill < xStart && nextXToFill <= lastXToFill ) {
                    outDataIndices[outputIndex++] = exteriorRegionCode;
                    ++nextXToFill;
                }

                int dataIndex = m_runData[run].dataIndex;
                if( dataIndex >= 0 ) {
                    dataIndex += nextXToFill - m_runData[run].x;
                    // If the run is defined, copy the distance values
                    while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                        outDataIndices[outputIndex++] = dataIndex++;
                        ++nextXToFill;
                    }
                } else {
                    // If the run is undefined, copy +/-levelSetOutsideDistance
                    while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                        outDataIndices[outputIndex++] = dataIndex;
                        ++nextXToFill;
                    }
                }
            }
            // If there are any values to fill after the last run, fill them in with the exterior value
            while( nextXToFill <= lastXToFill ) {
                outDataIndices[outputIndex++] = exteriorRegionCode;
                ++nextXToFill;
            }
        }
        // Fill in the ending B out-of-bounds scanlines with exterior values
        for( int i = 0; i < bMaxExteriorCount; ++i )
            outDataIndices[outputIndex++] = exteriorRegionCode;
    }

    // If the array extended beyond in the C direction, we need to fill the rest with exterior values
    while( outputIndex < dataSize )
        outDataIndices[outputIndex++] = exteriorRegionCode;
}

void rle_index_spec::fill_data_index_box_boundary_duplicated( const boundbox3& voxelExtents,
                                                              boost::int32_t* outDataIndices ) const {
    ////////////////
    // INITIALIZE
    ////////////////

    int exteriorRegionCode = m_exteriorRegionCode;
    // int dataSize = voxelExtents.get_volume();

    // Get the B coordinate range and the C coordinate
    int bMin = voxelExtents.minimum().y - m_abcCoordOrigin.y, bMax = voxelExtents.maximum().y - m_abcCoordOrigin.y;
    int cMin = voxelExtents.minimum().z - m_abcCoordOrigin.z, cMax = voxelExtents.maximum().z - m_abcCoordOrigin.z;

    // The outer bound box size in B and C directions
    int bsize = m_abcCoordSize.ysize(), csize = m_abcCoordSize.zsize();
    // The X coordinate range of the outer bound box
    int xRangeBegin = m_abcCoordOrigin.x, xRangeEnd = m_abcCoordOrigin.x + m_abcCoordSize.xsize();

    ////////////////
    // COPY LEVEL SET VALUES
    ////////////////

    // This is the value we use to step through the outVoxelCornerValues
    int outputIndex = 0;

    // Now iterate through all the C coordinates and copy the planes
    for( int cIter = cMin; cIter <= cMax; ++cIter ) {
        // Now iterate through all the B coordinates and copy the scanlines
        for( int bIter = bMin; bIter <= bMax; ++bIter ) {
            // Clamp b and c to the outer bounds
            int b = ( bIter < 0 ) ? 0 : ( ( bIter >= bsize ) ? ( bsize - 1 ) : bIter );
            int c = ( cIter < 0 ) ? 0 : ( ( cIter >= csize ) ? ( csize - 1 ) : cIter );

            // Get the run range
            int bcIndex = b + c * m_abcCoordSize.ysize();
            int nextXToFill = voxelExtents.xminimum(), lastXToFill = voxelExtents.xmaximum();
            int runRangeStart = BCIndexXtoRunIndex( bcIndex, nextXToFill ),
                runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;
            const run_data *rd = &m_runData[runRangeStart], *rdEnd = &m_runData[runRangeEnd];

            if( rd->x != ( rd + 1 )->x ) {
                // This is a non-empty run
                int xStart, xPastTheEnd = rd->x;
                for( ; rd != rdEnd; ++rd ) {
                    // Get the x range of this run
                    xStart = xPastTheEnd;
                    xPastTheEnd = ( rd + 1 )->x;
                    int32_t dataIndex = rd->dataIndex;

                    // If there are values before xStart, fill them up with either exterior values, or
                    // if the run starts at the outer bounds, copy that value
                    if( nextXToFill < xStart && nextXToFill <= lastXToFill ) {
                        if( rd->x == xRangeBegin ) {
                            // In this case we're duplicating the first run's data index or region code
                            do {
                                outDataIndices[outputIndex++] = dataIndex;
                                ++nextXToFill;
                            } while( nextXToFill < xStart && nextXToFill <= lastXToFill );
                        } else {
                            // In this case we're duplicating the exterior region code
                            do {
                                outDataIndices[outputIndex++] = exteriorRegionCode;
                                ++nextXToFill;
                            } while( nextXToFill < xStart && nextXToFill <= lastXToFill );
                        }
                    }

                    if( dataIndex >= 0 ) {
                        dataIndex += nextXToFill - rd->x;
                        // If the run is defined, copy the distance values
                        while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                            outDataIndices[outputIndex++] = dataIndex++;
                            ++nextXToFill;
                        }
                    } else {
                        // If the run is undefined, copy +/-levelSetOutsideDistance
                        while( nextXToFill < xPastTheEnd && nextXToFill <= lastXToFill ) {
                            outDataIndices[outputIndex++] = dataIndex;
                            ++nextXToFill;
                        }
                    }
                }
                // If there are any values to fill after the last run, fill them in with the exterior value
                if( nextXToFill <= lastXToFill ) {
                    if( rdEnd->x == xRangeEnd ) {
                        // Duplicate the data index or region code of the last voxel in the run
                        int32_t dataIndex = ( rdEnd - 1 )->dataIndex;
                        if( dataIndex >= 0 )
                            dataIndex += rdEnd->x - ( rdEnd - 1 )->x - 1;
                        do {
                            // In this case we're duplicating the last run's last data index or region code
                            outDataIndices[outputIndex++] = dataIndex;
                            ++nextXToFill;
                        } while( nextXToFill <= lastXToFill );
                    } else {
                        do {
                            // In this case we're duplicating the exterior region code
                            outDataIndices[outputIndex++] = exteriorRegionCode;
                            ++nextXToFill;
                        } while( nextXToFill <= lastXToFill );
                    }
                }
            } else {
                // In the case of an empty run, everything is the exterior region code
                for( ; nextXToFill <= lastXToFill; ++nextXToFill )
                    outDataIndices[outputIndex++] = exteriorRegionCode;
            }
        }
    }
}

// If all the voxels in the box have the same region code, return it, otherwise return 0.
// This function should generally take on the order of the number of YZ coordinates in the input box to run.
boost::int32_t rle_index_spec::get_unique_region_code( const boundbox3& box ) const {
    // Get the box that's actually within the outer bounds of the rle index spec
    boundbox3 intersectedBounds = outer_bounds();
    intersectedBounds.intersect_with( box );

    if( intersectedBounds.is_empty() )
        return m_exteriorRegionCode;

    boost::int32_t uniqueRegionCode = 0;

    for( int z = intersectedBounds.zminimum(); z <= intersectedBounds.zmaximum(); ++z ) {
        for( int y = intersectedBounds.yminimum(); y <= intersectedBounds.ymaximum(); ++y ) {
            int bcIndex = ( y - m_abcCoordOrigin.y ) + ( z - m_abcCoordOrigin.z ) * m_abcCoordSize.ysize();

            int x = intersectedBounds.xminimum();
            int runIndex = BCIndexXtoRunIndex( bcIndex, x );
            const run_data* rd = &m_runData[runIndex];

            // If x isn't contained in the run, then it means the exterior region code gets involved here
            if( x < rd->x || x >= ( rd + 1 )->x ) {
                // An empty run means this scanline is made out of the exterior region code
                if( uniqueRegionCode == 0 ) {
                    uniqueRegionCode = m_exteriorRegionCode;
                } else if( uniqueRegionCode != m_exteriorRegionCode ) {
                    return 0;
                }
            }
            // If it's not an empty run, step through the runs until we get to intersectedBounds.xmaximum()
            if( rd->x != ( rd + 1 )->x ) {
                int runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;
                const run_data* rdEnd = &m_runData[runRangeEnd];
                // Loop until the end of the run range, or until the run exceeds the x bounds
                while( rd != rdEnd && rd->x <= intersectedBounds.xmaximum() ) {
                    int dataIndex = rd->dataIndex;
                    // If it's a defined run, return 0
                    if( dataIndex >= 0 )
                        return 0;
                    // Otherwise check that the region code matches throughout
                    if( uniqueRegionCode == 0 ) {
                        uniqueRegionCode = dataIndex;
                    } else if( uniqueRegionCode != dataIndex ) {
                        return 0;
                    }
                    ++rd;
                }
                // If the runs didn't extend all the way to the end of the intersectedBounds x range, then the exterior
                // region code gets involved
                if( rd != rdEnd && rdEnd->x <= intersectedBounds.xmaximum() ) {
                    if( uniqueRegionCode == 0 ) {
                        uniqueRegionCode = m_exteriorRegionCode;
                    } else if( uniqueRegionCode != m_exteriorRegionCode ) {
                        return 0;
                    }
                }
            }
        }
    }
    return uniqueRegionCode;
}

void rle_index_spec::fill_data_index_map( const rle_index_spec& risMappingTarget,
                                          boost::int32_t* outDataIndexChannel ) const {
    // Cache the target exterior region code
    int targetExteriorRegionCode = risMappingTarget.get_exterior_region_code();

    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = m_bcToRunIndex[bcIndex], runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;

            // Only if there is a target is there anything interesting to do
            if( ( runRangeStart + 1 != runRangeEnd || m_runData[runRangeStart].x != m_runData[runRangeEnd].x ) ) {
                // Get the bcIndex for the target rle index spec
                int bTarget = b + m_abcCoordOrigin.y - risMappingTarget.m_abcCoordOrigin.y;
                int cTarget = c + m_abcCoordOrigin.z - risMappingTarget.m_abcCoordOrigin.z;

                if( (unsigned)bTarget < (unsigned)risMappingTarget.m_abcCoordSize.ysize() &&
                    (unsigned)cTarget < (unsigned)risMappingTarget.m_abcCoordSize.zsize() ) {
                    int bcIndexTarget = bTarget + cTarget * risMappingTarget.m_abcCoordSize.ysize();
                    int runRangeStartTarget = risMappingTarget.m_bcToRunIndex[bcIndexTarget],
                        runRangeEndTarget = risMappingTarget.m_bcToRunIndex[bcIndexTarget + 1] - 1;

                    // Only if this is a non-empty scanline is there anything interesting to do
                    int xBeginTargetScanline = risMappingTarget.m_runData[runRangeStartTarget].x,
                        xEndTargetScanline = risMappingTarget.m_runData[runRangeEndTarget].x;
                    if( ( runRangeStartTarget + 1 != runRangeEndTarget ||
                          xBeginTargetScanline != xEndTargetScanline ) ) {
                        const run_data *rd = &m_runData[runRangeStart], *rdEnd = &m_runData[runRangeEnd];
                        // Before we can go through the runs in lockstep, we have to go through the runs which intersect
                        // with parts outside of the target scanline, and target the output channel to the target
                        // exterior region code
                        if( rd->x < xBeginTargetScanline ) {
                            while( rd != rdEnd ) {
                                int dataIndex = rd->dataIndex;
                                if( dataIndex >= 0 ) {
                                    // Go through all the data indices in the defined run (up to xBeginTargetScanline),
                                    // setting the data value to the target exterior region code
                                    int dataIndexEnd =
                                        dataIndex + ( std::min )( (int)( rd + 1 )->x, xBeginTargetScanline ) - rd->x;
                                    while( dataIndex != dataIndexEnd )
                                        outDataIndexChannel[dataIndex++] = targetExteriorRegionCode;
                                }

                                if( ( rd + 1 )->x < xBeginTargetScanline )
                                    ++rd;
                                else
                                    break;
                            }
                        }

                        // Go through the runs in lockstep, setting values as necessary
                        const run_data *rdTarget = &risMappingTarget.m_runData[runRangeStartTarget],
                                       *rdTargetEnd = &risMappingTarget.m_runData[runRangeEndTarget];
                        while( rd != rdEnd ) {
                            int dataIndex = rd->dataIndex;
                            // If this is a defined run, we may target some values in the output channel
                            if( dataIndex >= 0 ) {
                                // Get the half-open interval [xBegin,xEnd) over which this run is defined
                                int xBegin = rd->x, xEnd = ( rd + 1 )->x;
                                // Increment the rdTarget variable until it is at a run which ends at least where the
                                // interval [xBegin,xEnd) starts.
                                while( rdTarget != rdTargetEnd && ( rdTarget + 1 )->x <= xBegin )
                                    ++rdTarget;
                                // Process all the intervals which intersect the interval [xBegin,xEnd).  That is, all
                                // the intervals which start before the target intervals end
                                while( rdTarget != rdTargetEnd ) {
                                    int dataIndexTarget = rdTarget->dataIndex;
                                    // Get the half-open interval [xBeginTarget,xEndTarget) over which this input run is
                                    // defined
                                    int xBeginTarget = rdTarget->x, xEndTarget = ( rdTarget + 1 )->x;
                                    int xBeginIntersected = ( std::max )( xBegin, xBeginTarget );
                                    // Get the data index range
                                    int dataIndexCopy = dataIndex + xBeginIntersected - xBegin;
                                    int dataIndexCopyEnd = dataIndex + ( std::min )( xEnd, xEndTarget ) - xBegin;
                                    // If this is an defined run, then we need to set channel values to output data
                                    // indices
                                    if( dataIndexTarget >= 0 ) {
                                        // Offset the target data index based on the sub-interval being processed
                                        dataIndexTarget += xBeginIntersected - xBeginTarget;
                                        // Loop through and set the channel values to the target data indexes for this
                                        // run
                                        while( dataIndexCopy != dataIndexCopyEnd )
                                            outDataIndexChannel[dataIndexCopy++] = dataIndexTarget++;
                                    } else {
                                        // Loop through and set the channel values to the target region code of this
                                        // undefined run
                                        while( dataIndexCopy != dataIndexCopyEnd )
                                            outDataIndexChannel[dataIndexCopy++] = dataIndexTarget;
                                    }
                                    if( ( rdTarget + 1 )->x < xEnd )
                                        ++rdTarget;
                                    else
                                        break;
                                }
                            }
                            // If the target runs have run out, part of the current rd run we're processing may
                            // be outside the [xBeginTargetScanline,xEndTargetScanline) interval, so we don't increment
                            // rd in that case.
                            if( rdTarget != rdTargetEnd )
                                ++rd;
                            else
                                break;
                        }

                        // Go through the rest of the rd runs to set the needed values to the target exterior region
                        // code
                        while( rd != rdEnd ) {
                            int dataIndex = rd->dataIndex;
                            if( dataIndex >= 0 ) {
                                // Go through all the data indices in the defined run (up to xBeginTargetScanline),
                                // setting the data value to the target exterior region code
                                int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                                // If this rd run intersects with the end of the target scanline,
                                // we have to adjust the start data index.
                                if( xEndTargetScanline > rd->x )
                                    dataIndex += xEndTargetScanline - rd->x;
                                while( dataIndex < dataIndexEnd )
                                    outDataIndexChannel[dataIndex++] = targetExteriorRegionCode;
                            }
                            ++rd;
                        }

                    } else {
                        // If the target scanline was empty, then fill the defined voxels for this scanline
                        // with the exterior region code of risMappingTarget.
                        for( const run_data *rd = &m_runData[runRangeStart], *rdEnd = &m_runData[runRangeEnd];
                             rd != rdEnd; ++rd ) {
                            // Check for a defined run
                            int dataIndex = rd->dataIndex;
                            if( dataIndex >= 0 ) {
                                // Go through all the data indices in the defined run, setting the data value to the
                                // target exterior region code
                                int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                                while( dataIndex != dataIndexEnd )
                                    outDataIndexChannel[dataIndex++] = targetExteriorRegionCode;
                            }
                        }
                    }

                } else {
                    // If the target bcIndex was out of range, then fill the defined voxels for this scanline
                    // with the exterior region code of risMappingTarget.
                    for( const run_data *rd = &m_runData[runRangeStart], *rdEnd = &m_runData[runRangeEnd]; rd != rdEnd;
                         ++rd ) {
                        // Check for a defined run
                        int dataIndex = rd->dataIndex;
                        if( dataIndex >= 0 ) {
                            // Go through all the data indices in the defined run, setting the data value to the target
                            // exterior region code
                            int dataIndexEnd = dataIndex + ( rd + 1 )->x - rd->x;
                            while( dataIndex != dataIndexEnd )
                                outDataIndexChannel[dataIndex++] = targetExteriorRegionCode;
                        }
                    }
                }
            }
        }
    }
}

void rle_index_spec::copy_to_legacy_flood_rle_region( std::vector<boost::int32_t>& outBCOffsets,
                                                      std::vector<boost::int16_t>& outRunBoundaries,
                                                      std::vector<boost::int32_t>& outRunDataIndices ) const {
    // The bcOffsets we can copy directly
    outBCOffsets.resize( m_bcToRunIndex.size() );
    outRunBoundaries.clear();
    outRunDataIndices.clear();

    int legacyExteriorRegionCode = ( m_exteriorRegionCode == -1 ) ? -1000 : -2000;
    int dataIndex = 0;

    // We need to add in runs at the start and end to completely fill the data within the bounding box.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int bcIndex = b + c * bsize;
            int runRangeStart = m_bcToRunIndex[bcIndex], runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 1;

            outBCOffsets[bcIndex] = (int)outRunBoundaries.size();

            // If it's not a single empty run, then process it
            if( runRangeStart + 1 != runRangeEnd || m_runData[runRangeStart].x != m_runData[runRangeStart + 1].x ) {
                // If the run doesn't start at the beginning of the bounding box, add a new run at the start
                if( m_runData[runRangeStart].x > m_abcCoordOrigin.x ) {
                    outRunBoundaries.push_back( boost::int16_t( m_abcCoordOrigin.x ) );
                    outRunDataIndices.push_back( legacyExteriorRegionCode );
                }
                // Copy all the runs in the middle
                for( int run = runRangeStart; run != runRangeEnd; ++run ) {
                    outRunBoundaries.push_back( boost::int16_t( m_runData[run].x ) );
                    dataIndex = m_runData[run].dataIndex;
                    if( dataIndex < 0 )
                        dataIndex = ( dataIndex == -1 ) ? -1000 : -2000;
                    outRunDataIndices.push_back( dataIndex );
                }
                // If the run doesn't stop at the end of the bounding box, add a new run at the end.
                if( m_runData[runRangeEnd].x != m_abcCoordOrigin.x + m_abcCoordSize.xsize() ) {
                    outRunBoundaries.push_back( boost::int16_t( m_runData[runRangeEnd].x ) );
                    outRunDataIndices.push_back( legacyExteriorRegionCode );
                }
                // Cap off the scanline at the end of the bounding box
                outRunBoundaries.push_back( boost::int16_t( m_abcCoordOrigin.x + m_abcCoordSize.xsize() ) );
                outRunDataIndices.push_back( -9999 ); // Legacy Flood uses this value in the capping off run positions.
            }
            // If it is a single empty run, then we make a single run spanning the bounding box made of the exterior
            // region code
            else {
                outRunBoundaries.push_back( boost::int16_t( m_abcCoordOrigin.x ) );
                outRunDataIndices.push_back( legacyExteriorRegionCode );
                outRunBoundaries.push_back( boost::int16_t( m_abcCoordOrigin.x + m_abcCoordSize.xsize() ) );
                outRunDataIndices.push_back( -9999 ); // Legacy Flood uses this value in the capping off run positions.
            }
        }
    }

    // Finish off the outBCOffsets array, pointing past the end of the run index data to provide the range of the last
    // scanline
    outBCOffsets[outBCOffsets.size() - 1] = (int)outRunBoundaries.size();
}

// finds the first defined voxel in the rle index spec
// Note: this does not run in constant time, it must check each undefined run that comes before the first defined run
rle_defined_iterator rle_index_spec::begin() const { return rle_defined_iterator( *this ); }

// points to one past the last value of the scan line defined by bcIndex
rle_defined_iterator rle_index_spec::end() const { return rle_defined_iterator( *this, true ); }

void rle_index_spec::combine_for_blend( const rle_index_spec& risFirst, const rle_index_spec& risSecond ) {
    // Set the exterior region code
    if( risFirst.m_exteriorRegionCode == risSecond.m_exteriorRegionCode )
        m_exteriorRegionCode = risFirst.m_exteriorRegionCode;
    else
        throw runtime_error( "rle_index_spec::combine_for_blend() - The two input rle_index_spec instances have "
                             "incompatible exterior region codes.  "
                             "They must be the same, but the first operand has exterior region code " +
                             lexical_cast<string>( risFirst.m_exteriorRegionCode ) +
                             " and the second operand has exterior region code " +
                             lexical_cast<string>( risSecond.m_exteriorRegionCode ) + "." );

    // Set the bounding box
    m_abcCoordOrigin = vector3::from_min( risFirst.m_abcCoordOrigin, risSecond.m_abcCoordOrigin );
    m_abcCoordSize = size3::from_bounds( m_abcCoordOrigin,
                                         vector3::from_max( risFirst.m_abcCoordOrigin + risFirst.m_abcCoordSize,
                                                            risSecond.m_abcCoordOrigin + risSecond.m_abcCoordSize ) -
                                             vector3( 1 ) );

    // Start with empty run index data
    m_runData.clear();

    // Resize the bcToRunIndex vector to match our new bounding box
    m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );

    int32_t currentDataIndex = 0;
    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b != bsize; ++b ) {
            int y = m_abcCoordOrigin.y + b, z = m_abcCoordOrigin.z + c;

            int bcIndex = b + c * bsize;
            m_bcToRunIndex[bcIndex] = (int)m_runData.size();

            // The default state is an undefined region of the exterior region code
            bool creatingUndefinedRun = true;
            int32_t creatingRegionCode = m_exteriorRegionCode;

            // Loop through all the common sub-intervals of these scanlines
            rle_pairwise_run_iterator i( risFirst, risFirst.y_to_b( y ), risFirst.z_to_c( z ), risSecond,
                                         risSecond.y_to_b( y ), risSecond.z_to_c( z ) ),
                ie;
            for( ; i != ie; ++i ) {
                if( i.get_first_data_index() >= 0 || i.get_second_data_index() >= 0 ||
                    i.get_first_data_index() != i.get_second_data_index() ) {
                    // If either input is defined, or the two inputs have conflicting undefined region codes, the result
                    // is defined
                    if( creatingUndefinedRun ) {
                        m_runData.push_back( run_data( i.get_xmin(), currentDataIndex ) );
                        creatingUndefinedRun = false;
                    }
                    currentDataIndex += i.get_xsize();
                } else {
                    // Otherwise, we have two equal undefined region codes, so the result is undefined with that region
                    // code
                    int32_t regionCode = i.get_first_data_index();
                    if( creatingUndefinedRun ) {
                        if( creatingRegionCode != regionCode ) {
                            m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                            creatingRegionCode = regionCode;
                        }
                    } else {
                        m_runData.push_back( run_data( i.get_xmin(), regionCode ) );
                        creatingUndefinedRun = true;
                        creatingRegionCode = regionCode;
                    }
                }
            }

            // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to optimize the
            // storage and access. If no runs were added, add a zero-sized run.
            if( (int)m_runData.size() > m_bcToRunIndex[bcIndex] ) {
                // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to optimize
                // the storage and access.
                if( !creatingUndefinedRun || creatingRegionCode != m_exteriorRegionCode ) {
                    // When the rle_pairwise_run_iterator loop is done, i.get_xmin() is the value 1 past
                    // the end of the last sub-interval, so it is the correct value to close of the m_runData array.
                    m_runData.push_back( run_data( i.get_xmin(), -1 ) );
                }
            } else {
                m_runData.push_back( run_data( 0, -1 ) );
                m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }

    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();
    m_dataSize = currentDataIndex;
}

void rle_index_spec::build_by_filling( const rle_index_spec& src ) {
    // create a temp copy to swap in at the end for exception safety
    rle_index_spec ris;
    // copy the basic members first
    ris.m_abcCoordOrigin = src.m_abcCoordOrigin;
    ris.m_abcCoordSize = src.m_abcCoordSize;
    ris.m_exteriorRegionCode = src.m_exteriorRegionCode;

    ris.m_bcToRunIndex.resize( src.m_abcCoordSize.ysize() * src.m_abcCoordSize.zsize() + 1 );

    // go over the runs and convert all undefined interior runs (dataIndex <-1) into defined runs
    // combining defined runs as necessary
    for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
        for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
            int bcIndex = b + c * ris.m_abcCoordSize.ysize();

            // set the run index for the current bcToRunIndex entry
            ris.m_bcToRunIndex[bcIndex] = (int)ris.m_runData.size();

            // Get the range of runs from the source
            int runRangeStart = src.m_bcToRunIndex[bcIndex];
            int runRangeEnd = src.m_bcToRunIndex[bcIndex + 1] - 1; // one past the end

            // Special case the single run scanline
            if( runRangeStart >= runRangeEnd - 1 ) {
                run_data runIndexData = src.m_runData[runRangeStart];

                // if there is no x range for the run entry then the run is empty so we can simply enter an empty run
                // as well, if the single run matches the exterior region code, then it should be a null run as it adds
                // no information
                if( runIndexData.x == src.m_runData[runRangeStart + 1].x ||
                    runIndexData.dataIndex == ris.m_exteriorRegionCode ) {
                    ris.m_runData.push_back( run_data( 0, -1 ) );
                    ris.m_runData.push_back( run_data( 0, -1 ) );
                }
                // we have an undefined run that we need to make a defined run for
                // or we have a defined run that we need to keep (which amounts to the same thing)
                else {

                    // if the single run is inside undefined we make a new defined run
                    runIndexData.dataIndex = (int)ris.m_dataSize;

                    // grab the terminating entry to get the end x value
                    run_data runEnd = src.m_runData[runRangeStart + 1];

                    // update the data size to reflect the allocation of the new defined run
                    ris.m_dataSize += runEnd.x - runIndexData.x;

                    // push the new run and its terminating entry
                    ris.m_runData.push_back( runIndexData );
                    ris.m_runData.push_back( runEnd );
                }
            }
            // Deal with multiple runs in the scanline
            else {

                run_data runStart;

                bool buildingRun = false;
                // go over all runs in the scanline
                run_data runIndexData;
                for( int run = runRangeStart; run < runRangeEnd; ++run ) {
                    runIndexData = src.m_runData[run];

                    /*if ( run == 208 ) {
                      cout << "\n\nundefined run index = " << run << " data index=" << runIndexData.dataIndex << ", x="
                    << runIndexData.x << endl; run_data temp = ris.m_runData[206]; cout << "PREV: data index=" <<
                    temp.dataIndex <<
                    ", x=" << temp.x << endl;

                      cout << "Data Size: " << ris.m_dataSize << endl;

                    }*/

                    // the run is defined
                    if( runIndexData.dataIndex >= 0 ) {
                        // if we aren't currently building a run we need to start a new one
                        if( !buildingRun ) {
                            // cout << "starting new run" << endl;
                            //  we start a new run
                            //  create a new defined run spanning the undefined run
                            //  the data index is the current data size that we have added so far
                            runStart.dataIndex = (int)ris.m_dataSize;
                            runStart.x = runIndexData.x;

                            ris.m_runData.push_back( runStart );
                            buildingRun = true;
                        }
                    }
                    // the run is undefined
                    else {
                        // the run is undefined and is opposite the exterior region code, so we fill it
                        if( runIndexData.dataIndex != ris.m_exteriorRegionCode ) {
                            // cout << "undefined run index = " << run << " data index=" << runIndexData.dataIndex <<
                            // endl;
                            //  if we aren't currently building a run we need to start a new one
                            if( !buildingRun ) {
                                // cout << "starting new run" << endl;
                                //  we start a new run
                                //  create a new defined run spanning the undefined run
                                //  the data index is the current data size that we have added so far
                                runStart.dataIndex = (int)ris.m_dataSize;
                                runStart.x = runIndexData.x;
                                ris.m_runData.push_back( runStart );
                                buildingRun = true;
                            }
                        } else {

                            if( buildingRun ) {
                                // we close off the run
                                // cout << "closing run" << endl;
                                ris.m_runData.push_back( runIndexData );

                                runIndexData.dataIndex = (int)ris.m_dataSize;
                                int offset = runIndexData.x - runStart.x;

                                ris.m_dataSize += offset;

                                buildingRun = false;
                            } else {
                                // we just need to copy the run and its an undefined run  so the data array size isn't
                                // updated
                                ris.m_runData.push_back( runIndexData );
                            }
                        }
                    }
                }

                // cout << "closing run" << endl;
                //  close off the last run
                runIndexData = src.m_runData[runRangeEnd];
                if( buildingRun ) {
                    runIndexData.dataIndex = (int)ris.m_dataSize;
                    int offset = runIndexData.x - runStart.x;
                    ris.m_dataSize += offset; // increase the datasize of cause by the building of the run
                }

                ris.m_runData.push_back( runIndexData );
            }
        }
    }
    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    ris.m_bcToRunIndex[ris.m_bcToRunIndex.size() - 1] = (int)ris.m_runData.size();

    // swap in the result
    swap( ris );
}

void rle_index_spec::build_with_trim_bounds( const rle_index_spec& ris,
                                             const frantic::graphics::boundbox3& trimVoxelBounds ) {
    // If the bounds are already within the trimming bounding box, then just do a copy
    if( trimVoxelBounds.contains( ris.outer_bounds() ) ) {
        *this = ris;
        return;
    }

    // Because we're modifying this rle_index_spec, we need to reset the cached ris_adjacency
    free_cached_adjacency();

    m_exteriorRegionCode = ris.m_exteriorRegionCode;
    boundbox3 outerBounds = ris.outer_bounds();
    outerBounds.intersect_with( trimVoxelBounds );
    // If the resulting bounding box is empty, just clear out everything
    if( outerBounds.is_empty() ) {
        clear();
        m_exteriorRegionCode = ris.m_exteriorRegionCode;
        return;
    }

    m_abcCoordOrigin = outerBounds.minimum();
    m_abcCoordSize = outerBounds.size();

    m_dataSize = 0;

    // Reserve the correct amount of space
    m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
    m_runData.clear();
    m_runData.reserve( ris.m_runData.size() );
    int newXBoundsBegin = m_abcCoordOrigin.x, newXBoundsEnd = m_abcCoordOrigin.x + m_abcCoordSize.xsize();

    int bInput = 0, cInput = 0;

    int currentDataIndex = 0;
    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        cInput = c + m_abcCoordOrigin.z - ris.m_abcCoordOrigin.z;
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b < bsize; ++b ) {
            bInput = b + m_abcCoordOrigin.y - ris.m_abcCoordOrigin.y;
            int bcIndex = b + c * bsize;
            m_bcToRunIndex[bcIndex] = (int)m_runData.size();

            // If the input bc coordinates are within range, there will be a scanline to copy
            if( (unsigned)bInput < (unsigned)ris.m_abcCoordSize.ysize() &&
                (unsigned)cInput < (unsigned)ris.m_abcCoordSize.zsize() ) {
                int bcIndexInput = bInput + cInput * ris.m_abcCoordSize.ysize();
                int runRangeStartInput = ris.m_bcToRunIndex[bcIndexInput],
                    runRangeEndInput = ris.m_bcToRunIndex[bcIndexInput + 1] - 1;
                const run_data *rdInput = &ris.m_runData[runRangeStartInput],
                               *rdInputEnd = &ris.m_runData[runRangeEndInput];
                if( rdInput->x == ( rdInput + 1 )->x ) {
                    // Create a zero-sized run, necessary for the coordinate queries to work properly
                    m_runData.push_back( run_data( 0, -1 ) );
                    m_runData.push_back( run_data( 0, -1 ) );
                } else {
                    int runXStart = 0, runXEnd = 0;
                    // Copy the scanline, trimming to the new X bounds.
                    for( ; rdInput != rdInputEnd; ++rdInput ) {
                        runXStart = rdInput->x;
                        runXEnd = ( rdInput + 1 )->x;
                        // Clip the run to the bounds
                        if( runXStart < newXBoundsBegin )
                            runXStart = newXBoundsBegin;
                        if( runXEnd > newXBoundsEnd )
                            runXEnd = newXBoundsEnd;
                        // If this run intersects with the new X bounds, add it
                        if( runXStart < runXEnd ) {
                            int dataIndex = rdInput->dataIndex;
                            if( dataIndex < 0 ) {
                                // For undefined runs, make an exact copy of the run index data
                                m_runData.push_back( run_data( runXStart, dataIndex ) );
                            } else {
                                // For defined runs, we add on to the currentDataIndex
                                m_runData.push_back( run_data( runXStart, currentDataIndex ) );
                                currentDataIndex += runXEnd - runXStart;
                            }
                        }
                    }
                    if( m_bcToRunIndex[bcIndex] == (int)m_runData.size() ) {
                        // All runs were clipped away, so we add a zero-sized run
                        m_runData.push_back( run_data( 0, -1 ) );
                        m_runData.push_back( run_data( 0, -1 ) );
                    } else {
                        // Use the last runXEnd to cap off the last run
                        m_runData.push_back( run_data( runXEnd, -1 ) );
                    }
                }
            } else {
                // Create a zero-sized run, necessary for the coordinate queries to work properly
                m_runData.push_back( run_data( 0, -1 ) );
                m_runData.push_back( run_data( 0, -1 ) );
            }
        }
    }
    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();
    m_dataSize = currentDataIndex;
}

namespace detail {

/**
 * This is a helper function for build_from_populated, which does a flood fill algorithm on the region codes
 * for the unpopulated voxels.  Generally, already undefined voxels beat out defined voxels for determining
 * an unpopulated voxel's region code, except if that undefined voxel is outside of the outer bounds.
 *
 * All of the 0 values in the regionCodesChannel will be changed to either -1 or -2.
 *
 * @param  ris  The rle_index_spec for the channels.
 * @param  regionCodesChannel  The channel specifying the flood filled region codes.
 * @param  signedDistanceChannel  The level set signed distance values.
 */
void flood_fill_for_build_from_populated( const rle_index_spec& ris, int32_t* regionCodesChannel,
                                          const float* signedDistanceChannel ) {
    const ris_adjacency& adj = ris.get_cached_adjacency();
    boundbox3 outerBounds = ris.outer_bounds();

    // In regionCodesChannel, we have:
    //  1: populated
    //  0: unpopulated, no region code yet assigned
    // -1: undefined outside has been assigned
    // -2: undefined inside has been assigned

    // Create the voxel processing queue, and preallocate enough memory so we won't need to do another allocation during
    // processing
    vector<int32_t> voxelProcessingQueue;
    voxelProcessingQueue.reserve( ris.data_size() );

    ////////
    // First Fill Pass
    ////////

    // First do a pass through all the defined voxels, filling in all the unpopulated
    // voxels which are adjacent to a populated or undefined voxel
    for( rle_defined_iterator i = ris.begin(), ie = ris.end(); i != ie; ++i ) {
        if( regionCodesChannel[i.get_data_index()] == 0 ) {
            const ris_adj_entry& rae = adj[i.get_data_index()];

            // Collect the adjacent region codes for defined and for undefined voxels
            int32_t definedRegionCode = 0, undefinedRegionCode = 0;
            // Do -X direction
            int32_t neighborDataIndex = rae.x_neg;
            if( neighborDataIndex < 0 ) {
                // Only count the neighbor as a flood fill source if it is within the outer bounds
                if( i.get_coord().x - 1 >= outerBounds.xminimum() )
                    undefinedRegionCode = neighborDataIndex;
            } else {
                // If the neighbor is populated, then assign -1 for positive (outside) distance,
                // -2 for negative (inside) distance.
                if( regionCodesChannel[neighborDataIndex] == 1 )
                    definedRegionCode = ( signedDistanceChannel[neighborDataIndex] >= 0 ) ? -1 : -2;
            }
            // Do +X direction
            neighborDataIndex = rae.x_pos;
            if( neighborDataIndex < 0 ) {
                // Only count the neighbor as a flood fill source if it is within the outer bounds
                if( i.get_coord().x + 1 <= outerBounds.xmaximum() )
                    undefinedRegionCode = neighborDataIndex;
            } else {
                // If the neighbor is populated, then assign -1 for positive (outside) distance,
                // -2 for negative (inside) distance.
                if( regionCodesChannel[neighborDataIndex] == 1 )
                    definedRegionCode = ( signedDistanceChannel[neighborDataIndex] >= 0 ) ? -1 : -2;
            }
            // Do -Y direction
            neighborDataIndex = rae.y_neg;
            if( neighborDataIndex < 0 ) {
                // Only count the neighbor as a flood fill source if it is within the outer bounds
                if( i.get_coord().y - 1 >= outerBounds.yminimum() )
                    undefinedRegionCode = neighborDataIndex;
            } else {
                // If the neighbor is populated, then assign -1 for positive (outside) distance,
                // -2 for negative (inside) distance.
                if( regionCodesChannel[neighborDataIndex] == 1 )
                    definedRegionCode = ( signedDistanceChannel[neighborDataIndex] >= 0 ) ? -1 : -2;
            }
            // Do +Y direction
            neighborDataIndex = rae.y_pos;
            if( neighborDataIndex < 0 ) {
                // Only count the neighbor as a flood fill source if it is within the outer bounds
                if( i.get_coord().y + 1 <= outerBounds.ymaximum() )
                    undefinedRegionCode = neighborDataIndex;
            } else {
                // If the neighbor is populated, then assign -1 for positive (outside) distance,
                // -2 for negative (inside) distance.
                if( regionCodesChannel[neighborDataIndex] == 1 )
                    definedRegionCode = ( signedDistanceChannel[neighborDataIndex] >= 0 ) ? -1 : -2;
            }
            // Do -Z direction
            neighborDataIndex = rae.z_neg;
            if( neighborDataIndex < 0 ) {
                // Only count the neighbor as a flood fill source if it is within the outer bounds
                if( i.get_coord().z - 1 >= outerBounds.zminimum() )
                    undefinedRegionCode = neighborDataIndex;
            } else {
                // If the neighbor is populated, then assign -1 for positive (outside) distance,
                // -2 for negative (inside) distance.
                if( regionCodesChannel[neighborDataIndex] == 1 )
                    definedRegionCode = ( signedDistanceChannel[neighborDataIndex] >= 0 ) ? -1 : -2;
            }
            // Do +Z direction
            neighborDataIndex = rae.z_pos;
            if( neighborDataIndex < 0 ) {
                // Only count the neighbor as a flood fill source if it is within the outer bounds
                if( i.get_coord().z + 1 <= outerBounds.zmaximum() )
                    undefinedRegionCode = neighborDataIndex;
            } else {
                // If the neighbor is populated, then assign -1 for positive (outside) distance,
                // -2 for negative (inside) distance.
                if( regionCodesChannel[neighborDataIndex] == 1 )
                    definedRegionCode = ( signedDistanceChannel[neighborDataIndex] >= 0 ) ? -1 : -2;
            }

            // Now, use undefinedRegionCode in preference to definedRegionCode
            int32_t regionCode = ( undefinedRegionCode < 0 ) ? undefinedRegionCode : definedRegionCode;
            // If we got a region code to fill in, then set it and copy the indexes from the temp queue
            // to the voxelProcessingQueue.
            if( regionCode < 0 ) {
                regionCodesChannel[i.get_data_index()] = regionCode;
                voxelProcessingQueue.push_back( i.get_data_index() );
            }
        }
    }

    ////////
    // Iterative Flood Fill Of Remaining Interior
    ////////

    // Now that any voxels which are adjacent to a defined voxel or an undefined voxel with a good region code
    // are done, we can use an iterative flood fill to finish the rest.  The queue is used so that when
    // a region has different region codes at either end, they run into each other in the middle instead of
    // one of the codes filling the whole region.
    // Note that voxelProcessingQueue is growing as this loop runs, so we can't grab a copy of
    // voxelProcessingQueue.size().
    for( size_t voxelProcessingIndex = 0; voxelProcessingIndex < voxelProcessingQueue.size(); ++voxelProcessingIndex ) {
        int32_t dataIndex = voxelProcessingQueue[voxelProcessingIndex];
        int32_t regionCode = regionCodesChannel[dataIndex];

        const ris_adj_entry& rae = adj[dataIndex];
        // Go through all 6 directions and make sure they're full
        for( int i = 0; i < 6; ++i ) {
            int32_t neighborDataIndex = rae[i];
            // If a neighbor is unpopulated, then set it to this region code, and add it to the queue for processing.
            if( neighborDataIndex >= 0 && regionCodesChannel[neighborDataIndex] == 0 ) {
                regionCodesChannel[neighborDataIndex] = regionCode;
                voxelProcessingQueue.push_back( neighborDataIndex );
            }
        }
    }

    ////////
    // Non-Filled Remaining Voxels
    ////////

    // Some region codes may not have been filled.  They get set to the exterior code.
    int32_t exteriorRegionCode = ris.get_exterior_region_code();
    for( size_t i = 0, ie = ris.data_size(); i != ie; ++i ) {
        if( regionCodesChannel[i] == 0 ) {
            regionCodesChannel[i] = exteriorRegionCode;
        }
    }
}

} // namespace detail

bool rle_index_spec::build_from_populated( const rle_index_spec& ris, const unsigned char* populatedChannel,
                                           const float* signedDistanceChannel, bool shortCircuitFullyPopulated ) {
    // If there are no defined voxels in the input rle_index_spec, just copy it and return
    if( ris.data_size() == 0 ) {
        if( shortCircuitFullyPopulated ) {
            return false;
        } else {
            *this = ris;
            return true;
        }
    }

    // Copy the populated channel array into a region codes array to be used for the flood fill.
    // In the region codes array, we will use:
    //  1: populated
    //  0: unpopulated, no region code yet assigned
    // -ve: The region code has been assigned.
    vector<int32_t> regionCodesArray( ris.data_size(), 1 );
    boost::int32_t* regionCodesChannel = &regionCodesArray[0];
    bool anyUnpopulated = false;
    // All the region codes were defaulted to 1 ("populated"), so we
    // just have to set any unpopulated ones to 0 for the flood fill to use later.
    for( size_t i = 0, ie = ris.data_size(); i != ie; ++i ) {
        if( populatedChannel[i] == 0 ) {
            anyUnpopulated = true;
            regionCodesChannel[i] = 0;
        }
    }

    // If all the voxels were populated, just copy the rle_index_spec and return
    if( !anyUnpopulated ) {
        if( shortCircuitFullyPopulated ) {
            return false;
        } else {
            *this = ris;
            return true;
        }
    }

    // Do the flood fill of the region codes
    detail::flood_fill_for_build_from_populated( ris, regionCodesChannel, signedDistanceChannel );

    // Now create the new rle_index_spec using the filled in region codes
    free_cached_adjacency();

    m_exteriorRegionCode = ris.m_exteriorRegionCode;
    m_abcCoordOrigin = ris.m_abcCoordOrigin;
    m_abcCoordSize = ris.m_abcCoordSize;

    m_dataSize = 0;

    // Reserve the correct amount of space
    m_bcToRunIndex.resize( m_abcCoordSize.ysize() * m_abcCoordSize.zsize() + 1 );
    m_runData.clear();
    m_runData.reserve( ris.m_runData.size() );

    int32_t currentDataIndex = 0;
    // Iterate through the voxel coordinates and build up the arrays.
    for( int c = 0, csize = m_abcCoordSize.zsize(); c != csize; ++c ) {
        for( int b = 0, bsize = m_abcCoordSize.ysize(); b < bsize; ++b ) {
            int bcIndex = b + c * bsize;
            m_bcToRunIndex[bcIndex] = (int)m_runData.size();

            int runRangeStartInput = ris.m_bcToRunIndex[bcIndex],
                runRangeEndInput = ris.m_bcToRunIndex[bcIndex + 1] - 1;
            const run_data *rdInput = &ris.m_runData[runRangeStartInput],
                           *rdInputEnd = &ris.m_runData[runRangeEndInput];
            if( rdInput->x == ( rdInput + 1 )->x ) {
                // Create a zero-sized run, necessary for the coordinate queries to work properly
                m_runData.push_back( run_data( 0, -1 ) );
                m_runData.push_back( run_data( 0, -1 ) );
            } else {
                // These two variables let us track the run we're creating, and merge multiple runs together.
                // Initialize the default to an undefined exterior region code run, because that's how
                // things are treated outside of any runs.
                bool creatingUndefinedRun = true;
                int creatingRegionCode = m_exteriorRegionCode;

                int xStartDefinedRun = 0;

                // Loop through all the runs in this scanline.
                for( ; rdInput != rdInputEnd; ++rdInput ) {
                    // Get the data index/region code for this run
                    int dataIndexInput = rdInput->dataIndex;

                    // If it's an undefined run
                    if( dataIndexInput < 0 ) {
                        // If we aren't already creating an undefined run of this region code, start a new
                        // undefined run
                        if( !creatingUndefinedRun || creatingRegionCode != dataIndexInput ) {
                            // Can just copy the run_data from the input, because it's identical.
                            m_runData.push_back( *rdInput );

                            // If we finished off a defined run, we have to count the defined voxels
                            if( !creatingUndefinedRun )
                                currentDataIndex += rdInput->x - xStartDefinedRun;

                            creatingUndefinedRun = true;
                            creatingRegionCode = dataIndexInput;
                        }
                    }
                    // Otherwise it's a defined run
                    else {
                        // Iterate through all the defined voxels, and test them for trimming
                        int xBegin = rdInput->x, xEnd = ( rdInput + 1 )->x;
                        for( int x = xBegin; x != xEnd; ++x, ++dataIndexInput ) {
                            int32_t newRegionCode = regionCodesChannel[dataIndexInput];
                            if( newRegionCode == 1 ) {
                                // In this case, this was a populated voxel, so it should stay a defined voxel
                                if( creatingUndefinedRun ) {
                                    // Start a new defined run
                                    m_runData.push_back( run_data( x, currentDataIndex ) );
                                    // save the X coordinate for the start of this run
                                    xStartDefinedRun = x;
                                    creatingUndefinedRun = false;
                                }
                            } else {
                                if( creatingUndefinedRun ) {
                                    // If necessary switch to a run new run with a different undefined region code
                                    if( creatingRegionCode != newRegionCode ) {
                                        m_runData.push_back( run_data( x, newRegionCode ) );
                                        creatingRegionCode = newRegionCode;
                                    }
                                } else {
                                    // Here we must finish off the previous defined run, adding the number of defined
                                    // voxels into currentDataIndex.
                                    currentDataIndex += x - xStartDefinedRun;
                                    // Start the new undefined run
                                    m_runData.push_back( run_data( x, newRegionCode ) );
                                    creatingUndefinedRun = true;
                                    creatingRegionCode = newRegionCode;
                                }
                            }
                        }
                    }
                }

                // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to optimize
                // the storage and access. If no runs were added, add a zero-sized run.
                if( (int)m_runData.size() > m_bcToRunIndex[bcIndex] ) {
                    // Only finish off the last run if it isn't an undefined exterior run.  If it is, cut it off to
                    // optimize the storage and access.
                    if( !creatingUndefinedRun || creatingRegionCode != m_exteriorRegionCode ) {
                        m_runData.push_back( *rdInputEnd );
                        // If we finished off a defined run, we have to count the defined voxels
                        if( !creatingUndefinedRun )
                            currentDataIndex += rdInputEnd->x - xStartDefinedRun;
                    }
                } else {
                    m_runData.push_back( run_data( 0, -1 ) );
                    m_runData.push_back( run_data( 0, -1 ) );
                }
            }
        }
    }
    // Finish off the bcToRunIndex array, pointing past the end of the runIndexData array to provide the range of the
    // last scanline
    m_bcToRunIndex[m_bcToRunIndex.size() - 1] = (int)m_runData.size();
    m_dataSize = currentDataIndex;

    //	stringstream sout;
    //	if( !check_consistency(sout) ) {
    //		throw runtime_error( "build_from_populated: Failed consistency check:\n" + sout.str() );
    //	}

    return true;
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
