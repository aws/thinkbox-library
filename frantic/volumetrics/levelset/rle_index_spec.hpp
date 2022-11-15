// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/size3.hpp>
#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/volumetrics/levelset/rle_defined_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_adjacency.hpp>

#include <frantic/geometry/trimesh3_scan_conversion.hpp>

namespace frantic {
namespace fluids {
class rle_voxel_field;
}
} // namespace frantic

namespace frantic {
namespace volumetrics {
namespace levelset {
class rle_level_set;
}
} // namespace volumetrics
} // namespace frantic

namespace frantic {
namespace volumetrics {
namespace implicitsurface {
namespace detail {
class brls_ris_friend;
}
} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic

namespace frantic {
namespace volumetrics {
namespace levelset {

namespace detail {
// Some implementation functions that need access to the rle_index_spec and rle_level_set internals.
void convert_intersections_to_level_set(
    const std::vector<std::vector<geometry::scan_conversion_intersection>>& intersectionDepths,
    const geometry::trimesh3& geometry, float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside,
    bool isHalfOpenSurface, int exteriorRegionCodeNeg, int exteriorRegionCodePos,
    const frantic::graphics::boundbox3& resultVoxelBounds, rle_level_set& outResult );

void finalize_geometry_to_levelset_conversion( int exteriorRegionCode, rle_level_set& levelSet );

void intermediate_mesh_to_level_set_union( rle_level_set& levelSet, const rle_level_set& inputLevelSet );

void mesh_to_level_set_region_x_dilation( rle_level_set& levelSet );

void mesh_to_level_set_region_x_sweeping( rle_level_set& levelSet );

void build_test1( rle_level_set& levelSet );

void build_testSuite( rle_level_set& levelSet, const frantic::graphics::vector3 orig,
                      const frantic::graphics::size3 coordSize, const int exteriorRegionCode );

// This builds a 2d level set that we use for testing in the Test Suite
void build_test1_complement( rle_level_set& levelSet );

// This builds a 2d level set that we use for testing in the Test Suite
void build_test2( rle_level_set& levelSet );

// This builds a 2d level set that we use for testing in the Test Suite
void build_test12_union( rle_level_set& levelSet );

// This builds a 2d level set that we use for testing in the Test Suite
void build_test12_intersect( rle_level_set& levelSet );

// This builds a 3d level set that is a plane
void build_plane_test( rle_level_set& levelSet );
} // namespace detail

/**
 * This struct stores the x positions and data index/run code of each run.
 */
struct run_data {
    boost::int32_t x;
    boost::int32_t dataIndex;

    run_data()
        : x( 0 )
        , dataIndex( 0 ) {}
    run_data( boost::int32_t x_, boost::int32_t dataIndex_ )
        : x( x_ )
        , dataIndex( dataIndex_ ) {}
};

// A vector<> of these is passed into the axis permutation function.  This specifies input and output
// data that should be remapped.
struct rle_index_spec_channel_copying_data {
    std::size_t primitiveSize;
    const char* inputData;
    char* outputData;
};

// forward declaration of rle_level_set and rle_index_spec_block_iterator_x so they can be a friend of rle_index_spec
class rle_level_set;
class rle_run_iterator;
class rle_pairwise_run_iterator;
class rle_defined_and_adj_iterator;
class rle_index_spec_block_iterator_x;

class rle_index_spec {
    // These two members define the bounding box within which all the voxels lie.
    frantic::graphics::vector3 m_abcCoordOrigin;
    frantic::graphics::size3 m_abcCoordSize;
    // This value defines the region code for any voxels that aren't specified by a run
    boost::int32_t m_exteriorRegionCode;

    // Maps an BCIndex to a RunIndex.
    std::vector<boost::int32_t> m_bcToRunIndex;

    // Maps a RunIndex to an X coordinate or region code and an index in the external data array
    std::vector<run_data> m_runData;

    // Keep track of the data size separately, so we don't have to recompute it all the time
    std::size_t m_dataSize;

    /**
     * This ris_adjacency instance defaults to NULL, and is created when requested externally.
     * It gets reset to NULL by any function which modifies the rle_index_spec.  In general, you don't
     * want to create and use this if you can avoid it, because it uses storage of 24 bytes per defined voxel.
     *
     * @note This member is declared as mutable so that it can be created even in a const function.  Note
     *       that this means we have to be careful about thread safety with respect to this variable.  Why is this
     *       ok to do in this case?  Because this variable only caches information about the state of the
     * rle_index_spec, and is not used to modify the state.  Any other usage of mutable is generally a bad idea.
     */
    mutable ris_adjacency* m_cachedRisAdjacency;

    // Give the rle_level_set direct access to the members of this class, so it can do efficient construction.
    friend class rle_level_set;
    friend class frantic::fluids::rle_voxel_field;
    friend class ris_adjacency;
    friend class rle_index_spec_block_iterator_x;
    friend class rle_defined_iterator;
    friend class rle_scanline_iteration_primitive;
    friend class rle_run_iterator;
    friend class rle_pairwise_run_iterator;
    friend class rle_defined_and_adj_iterator;
    friend class rle_defined_block_iterator;
    friend class rle_general_block_iterator;

    // Some friend functions
    friend void detail::convert_intersections_to_level_set(
        const std::vector<std::vector<geometry::scan_conversion_intersection>>& intersectionDepths,
        const geometry::trimesh3& geometry, float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside,
        bool isHalfOpenSurface, int exteriorRegionCodeNeg, int exteriorRegionCodePos,
        const frantic::graphics::boundbox3& resultVoxelBounds, rle_level_set& outResult );
    friend void detail::finalize_geometry_to_levelset_conversion( int exteriorRegionCode, rle_level_set& levelSet );
    friend void detail::intermediate_mesh_to_level_set_union( rle_level_set& levelSet,
                                                              const rle_level_set& inputLevelSet );
    friend void detail::mesh_to_level_set_region_x_dilation( rle_level_set& levelSet );
    friend void detail::mesh_to_level_set_region_x_sweeping( rle_level_set& levelSet );
    friend void detail::build_test1( rle_level_set& levelSet );
    friend void detail::build_testSuite( rle_level_set& levelSet, const frantic::graphics::vector3 orig,
                                         const frantic::graphics::size3 coordSize, const int exteriorRegionCode );
    friend void detail::build_test1_complement( rle_level_set& levelSet );
    friend void detail::build_test2( rle_level_set& levelSet );
    friend void detail::build_test12_union( rle_level_set& levelSet );
    friend void detail::build_test12_intersect( rle_level_set& levelSet );
    friend void detail::build_plane_test( rle_level_set& levelSet );
    // friend void convert_voxel_state_field_debug_mesh( const frantic::fluids::rle_voxel_field& levelSet, float
    // cubeSize, frantic::geometry::trimesh3& outMesh );

    friend class ::frantic::volumetrics::implicitsurface::detail::brls_ris_friend;

    ////////////////////////////////////////////////////////
    // Internal Coordinate system conversion functions (preliminary names, TODO: make better names)
    //   Note that these functions don't do validity checking of the input.  That's for the public functions.
    ////////////////////////////////////////////////////////

    frantic::graphics::vector3 XYZtoABC( const frantic::graphics::vector3& xyz ) const {
        return xyz - m_abcCoordOrigin;
    }

    std::pair<int, int> XYZtoBCIndexX( const frantic::graphics::vector3& xyz ) const {
        return std::make_pair( ( xyz.y - m_abcCoordOrigin.y ) + m_abcCoordSize.ysize() * ( xyz.z - m_abcCoordOrigin.z ),
                               xyz.x );
    }

    std::pair<int, int> BCIndexXtoRunIndexX( const std::pair<int, int>& bcIndexX ) const {
        return std::make_pair( BCIndexXtoRunIndex( bcIndexX.first, bcIndexX.second ), bcIndexX.second );
    }

    int RunIndexXtoDataIndex( const std::pair<int, int>& runIndexX ) const {
        return m_runData[runIndexX.first].dataIndex + ( runIndexX.second - m_runData[runIndexX.first].x );
    }

    int RunIndexXtoDataIndex( int runIndex, int x ) const {
        return m_runData[runIndex].dataIndex + ( x - m_runData[runIndex].x );
    }

  public:
    rle_index_spec() {
        m_dataSize = 0;
        m_exteriorRegionCode = -1;
        m_cachedRisAdjacency = 0;
    }

    ~rle_index_spec() { free_cached_adjacency(); }

    rle_index_spec( const rle_index_spec& ris )
        : m_abcCoordOrigin( ris.m_abcCoordOrigin )
        , m_abcCoordSize( ris.m_abcCoordSize )
        , m_exteriorRegionCode( ris.m_exteriorRegionCode )
        , m_bcToRunIndex( ris.m_bcToRunIndex )
        , m_runData( ris.m_runData )
        , m_dataSize( ris.m_dataSize )
        , m_cachedRisAdjacency( 0 ) // Do not copy this mutable member!
    {}

    /**
     * This function sets the rle_index_spec by swapping the input variables into the structure.  This prevents
     * duplicate copies of the big data chunks from ever being made, so is ideal for file I/O code to use.
     *
     * If requested, it also validates the result, and produces an error message if necessary
     * using the sourceObjectName to provide some context.  The preprocessor symbol RLE_INDEX_SPEC_DEBUG
     * can be defined to force this check in all cases.
     *
     * @param  checkConsistency    Whether to do the consistency check of the resulting rle index spec.
     * @param  outerBounds         The outer bounds of the structure.
     * @param  dataSize            The data array size.
     * @param  exteriorRegionCode  The exterior region code of the structure.
     * @param  bcToRunIndex        The array mapping BC coordinates to runIndex values.
     * @param  runIndexData        The array with data for each runIndex.
     * @param  sourceObjectName    A human-readable name for use in error messages.
     */
    void set_with_swap( bool checkConsistency, const frantic::graphics::boundbox3& outerBounds, std::size_t dataSize,
                        boost::int32_t exteriorRegionCode, std::vector<boost::int32_t>& bcToRunIndex,
                        std::vector<run_data>& runIndexData, const frantic::tstring& sourceObjectName ) {
        // Reset the cached ris_adjacency.
        free_cached_adjacency();

        // Swap all the values provided into a temporary rle_index_spec
        rle_index_spec ris;
        ris.m_abcCoordOrigin = outerBounds.minimum();
        ris.m_abcCoordSize = outerBounds.size();
        ris.m_bcToRunIndex.swap( bcToRunIndex );
        ris.m_runData.swap( runIndexData );
        ris.m_dataSize = dataSize;
        ris.m_exteriorRegionCode = exteriorRegionCode;

        // If RLE_INDEX_SPEC_DEBUG is defined, always do a consistency check
#ifndef RLE_INDEX_SPEC_DEBUG
        if( checkConsistency ) {
#endif
            // Check the resulting structure for consistency, and throw an error if necessary.
            std::stringstream errorMessage;
            errorMessage << "Consistency error with RLE Index Spec from \""
                         << frantic::strings::to_string( sourceObjectName ) << "\"\n";
            if( !ris.check_consistency( errorMessage ) ) {
                throw std::runtime_error( errorMessage.str() );
            }
#ifndef RLE_INDEX_SPEC_DEBUG
        }
#endif

        // Swap the temp rle_index_spec into this
        ris.swap( *this );
    }

    void clear() {
        free_cached_adjacency();

        m_abcCoordOrigin.set( 0 );
        m_abcCoordSize.set( 0 );
        m_bcToRunIndex.clear();
        m_bcToRunIndex.push_back( 0 );
        m_runData.clear();
        m_dataSize = 0;
        m_exteriorRegionCode = -1;
    }

    bool is_empty() const {
        return m_abcCoordSize.xsize() <= 0 || m_abcCoordSize.ysize() <= 0 || m_abcCoordSize.zsize() <= 0;
    }

    rle_index_spec& operator=( const rle_index_spec& ris ) {
        free_cached_adjacency();

        m_abcCoordOrigin = ris.m_abcCoordOrigin;
        m_abcCoordSize = ris.m_abcCoordSize;
        m_bcToRunIndex = ris.m_bcToRunIndex;
        m_runData = ris.m_runData;
        m_dataSize = ris.m_dataSize;
        m_exteriorRegionCode = ris.m_exteriorRegionCode;

        return *this;
    }

    bool operator==( const rle_index_spec& ris ) const {
        if( m_abcCoordOrigin != ris.m_abcCoordOrigin || m_abcCoordSize != ris.m_abcCoordSize ||
            m_bcToRunIndex != ris.m_bcToRunIndex || m_dataSize != ris.m_dataSize ||
            m_exteriorRegionCode != ris.m_exteriorRegionCode ) {
            return false;
        }
        // We do this special checking, so that the region code in the final run index data doesn't get compared, and
        // the region codes for empty runs don't get compared.  Note that the m_bcToRunIndex arrays have already been
        // confirmed to be equal, thus the run range will be equal across the two as well.
        for( int c = 0; c < ris.m_abcCoordSize.zsize(); ++c ) {
            for( int b = 0; b < ris.m_abcCoordSize.ysize(); ++b ) {
                int bcIndex = b + c * m_abcCoordSize.ysize();
                int runRangeStart = m_bcToRunIndex[bcIndex], runRangeEnd = m_bcToRunIndex[bcIndex + 1] - 2;

                // If it's not a single zero-sized run, go through all the runs in the range
                if( runRangeStart != runRangeEnd || m_runData[runRangeStart].x != m_runData[runRangeStart + 1].x ) {
                    for( int run = runRangeStart; run <= runRangeEnd; ++run ) {
                        if( m_runData[run].x != ris.m_runData[run].x ||
                            m_runData[run].dataIndex != ris.m_runData[run].dataIndex )
                            return false;
                    }
                    // Also check the value which finishes off the last run.  Don't check the region code, though.
                    if( m_runData[runRangeEnd + 1].x != ris.m_runData[runRangeEnd + 1].x )
                        return false;
                }
                // Otherwise confirm that the run is zero-sized in the comparison rle_index_spec as well.  The region
                // code doesn't matter for zero-sized runs.
                else {
                    if( ris.m_runData[runRangeStart].x != m_runData[runRangeStart + 1].x )
                        return false;
                }
            }
        }
        return true;
    }

    bool operator!=( const rle_index_spec& ris ) const { return !operator==( ris ); }

    void swap( rle_index_spec& ris ) {
        std::swap( m_abcCoordOrigin, ris.m_abcCoordOrigin );
        std::swap( m_abcCoordSize, ris.m_abcCoordSize );
        m_bcToRunIndex.swap( ris.m_bcToRunIndex );
        m_runData.swap( ris.m_runData );
        std::swap( m_dataSize, ris.m_dataSize );
        std::swap( m_exteriorRegionCode, ris.m_exteriorRegionCode );
        std::swap( m_cachedRisAdjacency, ris.m_cachedRisAdjacency );
    }

    /**
     * This function returns a reference to an ris_adjacency structure corresponding to this
     * rle_index_spec.  The rle_index_spec keeps a cached, dynamically allocated instance of this
     * around, and only allocates it when it is requested, because it uses a LOT of memory.  When
     * a function is called that mutates an rle_index_spec, its associated cached ris_adjacency
     * gets deallocated.
     *
     * IMPORTANT: Use this only when absolutely necessary!
     */
    const ris_adjacency& get_cached_adjacency() const;

    /**
     * Returns true if the cached ris_adjacency structure is currently allocated.
     */
    bool is_cached_adjacency() const { return m_cachedRisAdjacency != 0; }

    /**
     * Frees the cached ris_adjacency structure if it has been allocated.
     * It's const because the get_cached_adjacency allocation function is const,
     * and the actual adjacency object is mutable, so we need a const deallocation.
     */
    void free_cached_adjacency() const;

    /**
     * This function computes the bounding box around all the defined voxels in the rle_index_spec.
     */
    frantic::graphics::boundbox3 compute_defined_bounds() const;

    frantic::graphics::boundbox3 outer_bounds() const {
        return frantic::graphics::boundbox3( m_abcCoordOrigin, m_abcCoordSize );
    }

    graphics2d::boundrect2 outer_xy_bounds() const {
        return graphics2d::boundrect2( graphics2d::vector2( m_abcCoordOrigin.x, m_abcCoordOrigin.y ),
                                       graphics2d::size2( m_abcCoordSize.xsize(), m_abcCoordSize.ysize() ) );
    }

    /**
     * This function converts a Y voxel coordinate into a B coordinate within
     * the outer_bounds().
     *
     * @param  y  The Y voxel coordinate to convert.
     */
    boost::int32_t y_to_b( boost::int32_t y ) const { return y - m_abcCoordOrigin.y; }

    /**
     * This function converts a Z voxel coordinate into a C coordinate within
     * the outer_bounds().
     *
     * @param  z  The Z voxel coordinate to convert.
     */
    boost::int32_t z_to_c( boost::int32_t z ) const { return z - m_abcCoordOrigin.z; }

    std::size_t data_size() const { return m_dataSize; }

    boost::int32_t get_exterior_region_code() const { return m_exteriorRegionCode; }

    void set_exterior_region_code( boost::int32_t regionCode ) {
        if( regionCode >= 0 )
            throw std::runtime_error( "rle_index_spec::set_exterior_region_code() - The region code must be negative, "
                                      "cannot set to region code " +
                                      boost::lexical_cast<std::string>( regionCode ) + "." );
        m_exteriorRegionCode = regionCode;
    }

    const std::vector<boost::int32_t>& get_bc_to_run_index_vector() const { return m_bcToRunIndex; }

    const std::vector<run_data>& get_run_index_data_vector() const { return m_runData; }

    ////////////////////////////////////////////////////////
    // Low level access, oordinate system conversion functions (preliminary names, TODO: make better names)
    ////////////////////////////////////////////////////////

    /**
     * This binary search finds the run index such that m_runData[runIndex].x <= x < m_runData[runIndex+1].x,
     * which is the largest runIndex such that m_runData[runIndex].x <= bcIndexX.second.
     *
     * NOTE: If the X value provided isn't within the range spanned by the runs, the run returned will be valid but
     * won't contain that X value.
     *
     * OPTIMIZATION NOTE: If the X value provided is lower than any of the runs, the first run is returned. This way, if
     * you are doing a query for an interval of X values, you can do one binary search, and the first X value in the
     *                    interval isn't in the run returned, you can use this property to check the rest of the X
     * values against the same run, knowing if they are all outside its range then no run in this scanline contains any
     * of them.
     *
     * @param  bcIndex  The bcIndex of the input voxel coordinate.
     * @param  x  The X component of the input voxel coordiante.
     */
    int BCIndexXtoRunIndex( int bcIndex, int x ) const {
        // Determine the possible lower and upper run index values for the given BCIndex
        int lowRunIndex = m_bcToRunIndex[bcIndex];
        int highRunIndex = m_bcToRunIndex[bcIndex + 1] - 2;

#if 1
        // Optimization for low run count
        // TODO: Measure the performance of this once we have real level set data as well as access patterns from
        // important algorithms.
        if( highRunIndex - lowRunIndex < 4 ) {
            // Perform a linear search
            while( lowRunIndex < highRunIndex && x >= m_runData[lowRunIndex + 1].x )
                ++lowRunIndex;
            return lowRunIndex;
        } else {
#else
        // No optimization for low run count
        if( lowRunIndex < highRunIndex ) {
#endif
            // Do a binary search to find C within that range
            do {
                // Get the middle index, by rounding upwards
                int middleRunIndex = ( lowRunIndex + highRunIndex + 1 ) / 2;
                if( m_runData[middleRunIndex].x <= x ) {
                    lowRunIndex = middleRunIndex;
                } else {
                    highRunIndex = middleRunIndex - 1;
                }
            } while( lowRunIndex < highRunIndex );
        }

        return lowRunIndex;
    }

    /**
     * This function, given an interval along X, returns a closed interval of run indexes that contain
     * at a minimum that range.  Any part of the input X interval which isn't covered by the runs should
     * be considered to consist of undefined voxels with the exterior region code.
     *
     * It is possible that the return value has return.second < return.first, in which case there
     * are no runs.  It is also possible that the return value has return.first == return.second, and
     * that run does not intersect with the [xMin,xMax] interval.  This check is left for the caller, to
     * avoid redundant lookups.
     *
     * @param  xMin  The start of the X interval.
     * @param  xMax  The end of the X interval (inclusive).
     * @param  y  The Y coordinate to use.
     * @param  z  The Z coordinate to use.
     */
    std::pair<int, int> x_interval_to_run_index_interval( int xMin, int xMax, int y, int z ) const {
        // Compute the B and C coordinates from the xyz
        int b = y - m_abcCoordOrigin.y;
        int c = z - m_abcCoordOrigin.z;

        if( (unsigned)b < (unsigned)m_abcCoordSize.ysize() && (unsigned)c < (unsigned)m_abcCoordSize.zsize() ) {
            int bcIndex = b + m_abcCoordSize.ysize() * c;
            return std::make_pair( BCIndexXtoRunIndex( bcIndex, xMin ), BCIndexXtoRunIndex( bcIndex, xMax ) );
        } else {
            // Return an empty interval
            return std::make_pair( 1, 0 );
        }
    }

    /**
     * This function returns the data index of the first coordinate in the run, or
     * the region code of the run as a negative number
     *
     * @param  runIndex  The run index for which to retrieve the x interval
     */
    int run_index_to_data_index( int runIndex ) const { return m_runData[runIndex].dataIndex; }

    /**
     * This function returns a half-open interval of x coordinates for this run.
     * The returned X interval is [return.first, return.second).
     *
     * @param  runIndex  The run index for which to retrieve the x interval
     */
    std::pair<int, int> run_index_to_x_interval( int runIndex ) const {
        return std::make_pair( m_runData[runIndex].x, m_runData[runIndex + 1].x );
    }

    /**
     * Converts an XYZ voxel coordinate into a DataIndex or region code.  Non-negative values mean that the voxel is
     * defined, and the return value is an index into the data array.  Negative values mean the voxel is undefined, and
     * the returned value is the region code.  The special region code -1 indicates "outside", and any other region code
     * (-2 or below) indicates "inside". To test if a voxel is defined you check that the returned value is >= 0.
     *
     * @param  xyz  The voxel XYZ coordinate to convert to a data index.
     */
    int XYZtoDataIndex( const frantic::graphics::vector3& xyz ) const {
        // Compute the B and C coordinates from the xyz
        int b = xyz.y - m_abcCoordOrigin.y;
        int c = xyz.z - m_abcCoordOrigin.z;
        // Only test the B and C coordinate values initially.  The X coordinate value will be tested after the binary
        // search is done, based on the assumption that the expected case will often be to hit defined runs.
        if( (unsigned)b < (unsigned)m_abcCoordSize.ysize() && (unsigned)c < (unsigned)m_abcCoordSize.zsize() ) {
            int runIndex = BCIndexXtoRunIndex( b + m_abcCoordSize.ysize() * c, xyz.x );
            int runBeginX = m_runData[runIndex].x;
            // Now check that the X value is within the run we found
            if( runBeginX <= xyz.x && xyz.x < m_runData[runIndex + 1].x ) {
                int runDataIndex = m_runData[runIndex].dataIndex;
                if( runDataIndex >= 0 ) {
                    return runDataIndex + ( xyz.x - runBeginX );
                } else {
                    // Return the the region code, which is exactly the value stored in the data index variable for
                    // undefined runs.
                    return runDataIndex;
                }
            }
        }
        // m_exteriorRegionCode specifies the region outside the outer_bounds or outside the range of runs within a
        // scanline
        return m_exteriorRegionCode;
    }

    ////////////////////////////////////////////////////////
    // Higher level access, which often can aggregate operations for higher performance
    ////////////////////////////////////////////////////////

    /**
     * This function fills the provided data array (which must be of size 8) with the 8
     * value block having minimum index at minCorner.
     *
     * @param  minCornerXYZ  The voxel coordinate of the minimum coordinate in the box.  The box
     *                       is from minCornerXYZ to minCornerXYZ+[1,1,1].
     * @param  outDataIndices  An array of integers into which the data indices are placed.  This array must
     *                         be at least of size 8, there is no checking done in this function.
     */
    void fill_2x2x2_data_index_box( const frantic::graphics::vector3& minCornerXYZ,
                                    boost::int32_t* outDataIndices ) const;

    /**
     * This function fills the provided data array with the data indexes in the box provided.  The
     * array provided must be at least as big as voxelExtents.get_volume().
     *
     * @param  voxelExtents  The bounding box defining the block of data to be retrieved.
     * @param  outDataIndices  An array of integers into which the data indices are placed.  This array must
     *                         be at least of size voxelExtents.get_volume(), there is no checking done in this
     * function.
     */
    void fill_data_index_box( const frantic::graphics::boundbox3& voxelExtents, boost::int32_t* outDataIndices ) const;

    /**
     * This function fills the provided data array with the data indexes in the box provided.  The
     * array provided must be at least as big as voxelExtents.get_volume().  This version duplicates the boundary
     * data indices for voxel coordinates outside the outer bounds, instead of using the exterior region code for them.
     *
     * @param  voxelExtents  The bounding box defining the block of data to be retrieved.
     * @param  outDataIndices  An array of integers into which the data indices are placed.  This array must
     *                         be at least of size voxelExtents.get_volume(), there is no checking done in this
     * function.
     */
    void fill_data_index_box_boundary_duplicated( const frantic::graphics::boundbox3& voxelExtents,
                                                  boost::int32_t* outDataIndices ) const;

    /**
     * If all the voxels in the box have the same region code, return it, otherwise return 0.
     * This function can be used to tell whether a specified box region can be skipped because it's entirely outside,
     * processed trivially because it's entirely inside, or needs to be subdivided because there's likely defined
     * data inside.
     *
     * @param  box  The voxel bounding box within which to check for a unique region code.
     */
    boost::int32_t get_unique_region_code( const frantic::graphics::boundbox3& box ) const;

    /**
     * This function fills in an int32 channel with a mapping to the destination rle_index_spec.  What this
     * means is for each voxel in this rle_index_spec, the resulting channel will either have the defined voxel
     * index in the destination rle_index_spec, or the voxel region code in the destination rle_index_spec.
     *
     * @param  risMappingTarget  The rle_index_spec which the filled in data indexes point to.
     * @param  outDataIndexChannel  An int32 array, which must contain at least this->data_size() elements, and
     *                              gets filled in with the destination data indexes and region codes.
     */
    void fill_data_index_map( const rle_index_spec& risMappingTarget, boost::int32_t* outDataIndexChannel ) const;

    ////////////////////////////////////////////////////////
    // Validation and dumping functions
    ////////////////////////////////////////////////////////

    /**
     * Dumps the contents of this rle_index_spec.
     *
     * @param  out  The output stream to which dump data is written.
     */
    void dump( std::ostream& out ) const;

    /**
     * This function checks that the rle index spec is self-consistent, and prints out
     * information about how it isn't consistent to the provided stream.
     *
     * @param  out  The output stream to which error messages are written.
     *
     * @return  True if it checks out, false otherwise.
     */
    bool check_consistency( std::ostream& out ) const;

    /**
     * This function iterates over every run within the outer bounds, going over every voxel, and confirming the value
     * with the value given by the XYZtoDataIndex function.
     *
     * @param  out  The output stream to which error messages are written.
     *
     * @return  True if it checks out, false otherwise.
     */
    bool check_voxels_slow( std::ostream& out ) const;

    ////////////////////////////////////////////////////////
    // RLE index spec manipulation functions
    ////////////////////////////////////////////////////////

    /**
     * This function applies the given axis permutation to the rle_index_spec, while at the same time as copying
     * the specified data channels to match the new configuration.
     *
     * @param  axisPermutation  The permutation to apply to the axes.
     * @param  channelsToRemapBegin  Specification of the channels to remap.  This is the begin iterator.
     *								 The iterators are pointers for performance reasons.
     * @param  channelsToRemapEnd  Specification of the channels to remap.  This is the end iterator.
     *							   The iterators are pointers for performance reasons.
     */
    void apply_axis_permutation( const frantic::graphics::vector3& axisPermutation,
                                 const rle_index_spec_channel_copying_data* channelsToRemapBegin = 0,
                                 const rle_index_spec_channel_copying_data* channelsToRemapEnd = 0 );

    /**
     * This function copies the data formatted for the inputRIS into data formatted for the outputRIS.  It doesn't
     * zero out the elements that exist in outputRIS but not in inputRIS, so you will want to memset the values to
     * zero before calling this function.
     *
     * @param  outputRIS  The rle_index_spec defining how the output data arrays are mapped.
     * @param  inputRIS  The rle_index_spec defining how the input data arrays are mapped.
     * @param  channelsToCopyBegin  Specifications of the data arrays to transfer.  This is the begin iterator.  The
     * iterators are pointers for performance reasons.
     * @param  channelsToCopyEnd  Specifications of the data arrays to transfer.  This is the end iterator.  The
     * iterators are pointers for performance reasons.
     */
    static void copy_data_channels( const rle_index_spec& outputRIS, const rle_index_spec& inputRIS,
                                    const rle_index_spec_channel_copying_data* channelsToCopyBegin,
                                    const rle_index_spec_channel_copying_data* channelsToCopyEnd );

    /**
     * These two function provide an iterator pair that goes over all the defined, providing both the
     * voxel coordinate and the data index for each defined voxel.
     *
     * @note The iterator is not guaranteed to run in constant time per increment operation,
     *       because there might be arbitrarily many undefined runs interspersed in the index spec.  In
     *       a typical rle_index_spec, however, this will be the case.
     */
    rle_defined_iterator begin() const;
    rle_defined_iterator end() const;

    ////////////////////////////////////////////////////////
    // RLE index spec creation functions
    ////////////////////////////////////////////////////////

    /**
     * This function build an rle_index_spec from an array of voxel coordinates.  It may change the order of the
     * voxelArray during processing.
     *
     * @note  This function is generally only for testing purposes
     *
     * @param  coordinateArray  An array containing all the voxel coordinates with which to make up the defined
     *                         voxels of the rle index spec.
     */
    void build_from_voxel_array( std::vector<frantic::graphics::vector3>& coordinateArray );

    /**
     * This function selects a randomized sub-box of testBounds, then generates a randomized
     * rle index spec with that as its outer_bounds.  Region codes of undefined runs
     * are randomly selected, as is the exterior region code.
     *
     * @note  This function is generally only for testing purposes
     *
     * @param  testBounds  The bounding box within which the test is being done.
     * @param  expectedRunLength  How long to try and make the runs.
     * @param  randomRegionCodeLimit  How many random region codes to generate.  Set it to 2 if you just want inside
     * (-2) and outside (-1).
     */
    void build_from_random( const frantic::graphics::boundbox3& testBounds, float expectedRunLength,
                            int randomRegionCodeLimit = 100000 );

    /**
     * This applies a permutation to the voxel coordinates of the input rle_index_spec, producing a copy.
     * It can be used to swap the X and Z coordinates, for instance, to recompress the structure on a different axis.
     *
     * @param  axisPermutation  This is a permutation vector which is used to permute the axes.
     * @param  ris  This is the input rle_index_spec.
     * @param  channelsToRemapBegin  This is the begin of the iterator range of the channel copying data.  Note that
     *                               the iterator is just a pointer, for maximum performance.
     * @param  channelsToRemapEnd  This is the end of the iterator range of the channel copying data.
     */
    void build_from_axis_permutation( const frantic::graphics::vector3& axisPermutation, const rle_index_spec& ris,
                                      const rle_index_spec_channel_copying_data* channelsToRemapBegin = 0,
                                      const rle_index_spec_channel_copying_data* channelsToRemapEnd = 0 );

    /**
     * This function will be deleted, but is still here for now to do comparison tests.
     */
    void build_from_axis_permutation_unthreaded( const frantic::graphics::vector3& axisPermutation,
                                                 const rle_index_spec& ris,
                                                 const rle_index_spec_channel_copying_data* channelsToRemapBegin = 0,
                                                 const rle_index_spec_channel_copying_data* channelsToRemapEnd = 0 );

    /**
     * This function builds this rle_index_spec from a dilation of the input rle_index_spec.  That is, it applies
     * a box filter of width 2*dilationVoxels + 1 to the defined region.
     *
     * @param  ris  This is the input rle_index_spec.
     * @param  dilationVoxels  This is the number of voxels by which to dilate.  It should be 1 or greater.
     */
    void build_from_dilation( const rle_index_spec& ris, int dilationVoxels );

    /**
     * This function builds this rle_index_spec from a dilation of the input rle_index_spec.  That is, it applies
     * the box filter specified by the dilation box.  To be equivalent to the build_from_dilation which takes
     * a single integer, you would use the dilationBox
     * boundbox3(-dilationVoxelsX,dilationVoxelsX,-dilationVoxelsY,dilationVoxelsY,-dilationVoxelsZ,dilationVoxelsZ).
     *
     * @param  ris  This is the input rle_index_spec.
     * @param  dilationBox  A bounding box, which must contain the origin, specifying the dilation to do in each
     *                      of the six dilation directions.
     */
    void build_from_dilation( const rle_index_spec& ris, const frantic::graphics::boundbox3& dilationBox );

    /**
     * This function builds this rle_index_spec from an X dilation of the input rle_index_spec.  Generally
     * not useful on its own, it is used as a building block for build_from_dilation.
     *
     * @param  ris  This is the input rle_index_spec.
     * @param  dilationVoxelsNeg  This is the number of voxels by which to dilate towards negative X.  It should be 0 or
     * greater.
     * @param  dilationVoxelsPos  This is the number of voxels by which to dilate towards positive X.  It should be 0 or
     * greater.
     */
    void build_from_x_dilation( const rle_index_spec& ris, int dilationVoxelsNeg, int dilationVoxelsPos );

    /**
     * This function combines all the defined voxels of the two input rle_index_specs.
     * Basically, if a voxel is defined in either of the inputs, it will be defined in the
     * output.  Also, if a voxel has conflicting region codes, it will become defined.  That
     * property makes this function suitable to create the region for an interpolation or other
     * type of blending between two level sets.
     *
     * If the two inputs have different exterior region codes, an error will be thrown.
     *
     * @param  risA  The first input rle_index_spec.
     * @param  risB  The second input rle_index_spec.
     */
    void combine_for_blend( const rle_index_spec& risA, const rle_index_spec& risB );

    /**
     * This function converts all undefined runs that do not match the exterior region code
     * into defined runs.
     *
     * @param  ris  The input rle_index_spec.
     */
    void build_by_filling( const rle_index_spec& ris );

    /**
     * This function builds a new rle_index_spec, as a copy of the input ris but intersecting
     * the outer bounds with trimVoxelBounds.
     *
     * @param  ris  The input rle_index_spec.
     * @param  trimVoxelBounds  The bounding box with which to intersect the outer bounds.
     */
    void build_with_trim_bounds( const rle_index_spec& ris, const frantic::graphics::boundbox3& trimVoxelBounds );

    /**
     * This function builds a new rle_index_spec, converting any defined voxels indicated
     * by a value of 0 in the populated channel to undefined voxels.  The algorithm uses
     * a flood fill to determine what region code the newly undefined voxels should have.
     *
     * @param  ris  The input rle_index_spec.
     * @param  populatedChannel  The data for a 'populated' channel, with 0 for unpopulated
     *                           and other values for populated.
     * @param  signedDistanceChannel  The data for a signed distance channel, which is used
     *                                in the flood fill part of the algorithm to determine region
     *                                codes for newly undefined voxels.
     * @param  shortCircuitFullyPopulated  If this is true, the function does nothing if the voxels were
     *                                     fully populated, and returns false in that case.  The caller
     *                                     can use this to optimize the case where there is no change
     *                                     to the rle index spec.
     *
     * @return  True if the rle_index_spec was built, false otherwise.  If shortCircuitFullyPopulated is false,
     *          it will always return true.
     */
    bool build_from_populated( const rle_index_spec& ris, const unsigned char* populatedChannel,
                               const float* signedDistanceChannel, bool shortCircuitFullyPopulated );

    /**
     * Function for compatibility with the old Flood#.  To be deleted.
     */
    void build_from_legacy_flood_rle_region( int compressionAxis, const frantic::graphics::boundbox3& voxelBounds,
                                             const std::vector<boost::int32_t>& bcOffsets,
                                             const std::vector<boost::int16_t>& runBoundaries,
                                             const std::vector<boost::int32_t>& runDataIndices,
                                             frantic::graphics::vector3& outRequiredAxisPermutation );

    ////////////////////////////////////////////////////////
    // RLE index spec output functions
    ////////////////////////////////////////////////////////

    /**
     * Function for compatibility with the old Flood#.  To be deleted.
     */
    void copy_to_legacy_flood_rle_region( std::vector<boost::int32_t>& outBCOffsets,
                                          std::vector<boost::int16_t>& outRunBoundaries,
                                          std::vector<boost::int32_t>& outRunDataIndices ) const;
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
