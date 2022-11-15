// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

// Forward declaration of the rle_index_spec classes.
class rle_index_spec;
struct run_data;

/**
 * This rle_index_spec iteration primitive is intended to provide a building block for more complex rle_index_spec
 * iterators. Its requirements, in addition to minimal storage (i.e. no redundancy when a user of this class has
 * multiple instances) and efficient execution are:
 *
 * <pre>
 *  Initialization:
 *   - initialize to a scanline, given a range of run_data pointers
 *  Operations:
 *   - increment to next defined X
 *  Available Properties:
 *   - the next X where a new run starts
 *  Values Maintained/Tracked externally:
 *   - the rle_index_spec
 *   - the data index of the current voxel
 *   - the X, Y and Z coordinates of the current voxel
 * </pre>
 *
 * This primitive is used to iterate over all the defined voxels in a scanline.
 */
class rle_scanline_defined_iteration_primitive {
    const run_data *m_rd, *m_rdEnd;
    boost::int32_t m_xEnd;

  public:
    //////
    // Constructors
    //////

    /**
     * Initializes the iteration primitive to an empty default.
     */
    rle_scanline_defined_iteration_primitive()
        : m_rd( 0 )
        , m_rdEnd( 0 )
        , m_xEnd( 0 ) {}

    /**
     * Copy constructor.
     */
    rle_scanline_defined_iteration_primitive( const rle_scanline_defined_iteration_primitive& rhs )
        : m_rd( rhs.m_rd )
        , m_rdEnd( rhs.m_rdEnd )
        , m_xEnd( rhs.m_xEnd ) {}

    /**
     * Initializes the iteration primitive to the range of run_data values provided.  Note that if the function returns
     * false, the output variables can't be assumed to have any particular values.
     *
     * @param  ris      The rle_index_spec with which to initialize the primitive.  Note that it must stay valid as long
     * as the primitive exists.
     * @param  rdBegin  A pointer to the first run_data for the scanline being processed
     * @param  rdEnd    A pointer to the one-past-the-end run_data for the scanline being processed
     * @param[out]  outXCurrent  The X coordinate of the iteration primitive starting voxel is placed here.
     * @param[out]  outDataIndex  The data index of the iteration primitive starting voxel is placed here.
     * @return  Returns true if there was a defined voxel to start at, false if not.
     */
    bool init( const rle_index_spec& ris, const run_data* rdBegin, const run_data* rdEnd, boost::int32_t& outXCurrent,
               boost::int32_t& outDataIndex );

    /**
     * Initializes the iteration primitive to the range of run_data values provided.  Note that if the function returns
     * false, the output variables can't be assumed to have any particular values.
     *
     * @param  ris      The rle_index_spec with which to initialize the primitive.  Note that it must stay valid as long
     * as the primitive exists.
     * @param  rdBegin  A pointer to the first run_data for the scanline being processed
     * @param  rdEnd    A pointer to the one-past-the-end run_data for the scanline being processed
     * @param[out]  outXCurrent  The X coordinate of the iteration primitive starting voxel is placed here.
     * @param[out]  outXNegDataIndex  This is populated with the data index of the voxel with one smaller X coordinate
     * than the one resulting from the increment.
     * @param[out]  outDataIndex  The data index of the iteration primitive starting voxel is placed here.
     * @param[out]  outXPosDataIndex  This is populated with the data index of the voxel with one greater X coordinate
     * than the one resulting from the increment.
     * @return  Returns true if there was a defined voxel to start at, false if not.
     */
    bool init( const rle_index_spec& ris, const run_data* rdBegin, const run_data* rdEnd, boost::int32_t& outXCurrent,
               boost::int32_t& outXNegDataIndex, boost::int32_t& outDataIndex, boost::int32_t& outXPosDataIndex );

    /**
     * Clears the iteration primitive.  After this call, this->finished() will return true.
     *
     */
    void clear();

    //////
    // Functions to increment X along the scanline
    //////

    /**
     * Increments the iteration primitive to the next defined X coordinate in the scanline.
     *
     * @param  ris  The rle_index_spec which was used to initialize the primitive.
     * @param[in,out]  xCurrent  The X coordinate that the iteration primitive was at previously (should be the output
     * value of the last call to an init() or increment...() function.  This value is updated to contain the new X
     * coordinate.
     * @param[in,out]  dataIndex  The data index that the iteration primitive was at previously.  This value is updated
     * to contain the new data index.
     * @return  Returns true if the increment jumped to a new defined voxel, false if it passed the end of the scanline.
     */
    bool increment_to_next_defined_x( const rle_index_spec& ris, boost::int32_t& xCurrent, boost::int32_t& dataIndex );

    /**
     * Increments the iteration primitive to the next defined X coordinate in the scanline, and provides
     * the data indices of the immediate neighbors along the negative and positive X directions.
     *
     * @param  ris  The rle_index_spec which was used to initialize the primitive.
     * @param[in,out]  xCurrent  The X coordinate that the iteration primitive was at previously (should be the output
     * value of the last call to an init() or increment...() function.  This value is updated to contain the new X
     * coordinate.
     * @param[out]  outXNegDataIndex  This is populated with the data index of the voxel with one smaller X coordinate
     * than the one resulting from the increment.
     * @param[in,out]  dataIndex  The data index that the iteration primitive was at previously.  This value is updated
     * to contain the new data index.
     * @param[out]  outXPosDataIndex  This is populated with the data index of the voxel with one greater X coordinate
     * than the one resulting from the increment.
     * @return  Returns true if the increment jumped to a new defined voxel, false if it passed the end of the scanline.
     */
    bool increment_to_next_defined_x( const rle_index_spec& ris, boost::int32_t& xCurrent,
                                      boost::int32_t& outXNegDataIndex, boost::int32_t& dataIndex,
                                      boost::int32_t& outXPosDataIndex );

    bool jump_to_next_defined_x( const rle_index_spec& /*ris*/, boost::int32_t& xTarget, boost::int32_t& xCurrent,
                                 boost::int32_t& dataIndex );

    //////
    // Properties of the iteration
    //////

    // The X coordinate and data index are stored by the client of this object, so there are no getters for them.

    /**
     * Whether the iteration primitive has gone past all of the runs it had to process.
     */
    bool finished() const { return m_rd == 0; }
};

/**
 * This rle_index_spec iteration primitive is intended to provide a building block for more complex rle_index_spec
 * iterators. Its requirements, in addition to minimal storage (i.e. no redundancy when a user of this class has
 * multiple instances) and efficient execution are:
 *
 * <pre>
 *  Initialization:
 *   - initialize to some (X,Y,Z) coordinate
 *  Operations:
 *   - skip to a specified X (greater than the current X)
 *  Available Properties:
 *   - the next X where a new run starts
 *  Values Maintained/Tracked externally:
 *   - the rle_index_spec
 *   - the data index of the current voxel
 *   - the X, Y and Z coordinates of the current voxel
 * </pre>
 */
class rle_scanline_iteration_primitive {
    const run_data *m_rd, *m_rdEnd;
    boost::int32_t m_xEnd;

  public:
    //////
    // Constructors
    //////

    /**
     * Initializes the iteration primitive to an empty default.
     */
    rle_scanline_iteration_primitive()
        : m_rd( 0 )
        , m_rdEnd( 0 )
        , m_xEnd( 0 ) {}

    /**
     * Copy constructor.
     */
    rle_scanline_iteration_primitive( const rle_scanline_iteration_primitive& rhs )
        : m_rd( rhs.m_rd )
        , m_rdEnd( rhs.m_rdEnd )
        , m_xEnd( rhs.m_xEnd ) {}

    /**
     * Initializes the iteration primitive to the specified X coordinate in the specified scanline.
     *
     * @param  ris  The rle_index_spec with which to initialize the primitive.  Note that it must stay valid as long as
     *              the primitive exists.
     * @param  b  The b coordinate of the scanline to iterate over.
     * @param  c  The c coordinate of the scanline to iterate over.
     * @param  xCurrent  The x coordinate of the voxel to start at.
     * @param[out]  outDataIndex  The data index of the iteration primitive starting voxel is placed here.
     */
    void init( const rle_index_spec& ris, boost::int32_t b, boost::int32_t c, boost::int32_t xCurrent,
               boost::int32_t& outDataIndex );

    /**
     * Initializes the iteration primitive to the specified X coordinate in the specified scanline.
     *
     * @param  ris  The rle_index_spec with which to initialize the primitive.  Note that it must stay valid as long as
     *              the primitive exists.
     * @param  b  The b coordinate of the scanline to iterate over.
     * @param  c  The c coordinate of the scanline to iterate over.
     * @param  xCurrent  The x coordinate of the voxel to start at.
     * @param[out]  outXNegDataIndex  This is populated with the data index of the voxel with one smaller X coordinate
     * than the one resulting from the increment.
     * @param[out]  outDataIndex  The data index of the iteration primitive starting voxel is placed here.
     * @param[out]  outXPosDataIndex  This is populated with the data index of the voxel with one greater X coordinate
     * than the one resulting from the increment.
     */
    void init( const rle_index_spec& ris, boost::int32_t b, boost::int32_t c, boost::int32_t xCurrent,
               boost::int32_t& outXNegDataIndex, boost::int32_t& outDataIndex, boost::int32_t& outXPosDataIndex );

    /**
     * Clears the iteration primitive.  After this call, this->finished() will return true.
     *
     */
    void clear();

    //////
    // Functions to increment X along the scanline
    //////

    /**
     * Skips ahead (forward only) to the specified X in the same scanline.
     *
     * @param  xPrevious  The x coordinate of the voxel the iteration primitive is at before this increment.
     * @param  xNew  The x coordinate of the voxel the iteration primitive is at after this increment.
     * @param  ris  The rle_index_spec which was used to initialize the primitive.
     * @param[in,out]  dataIndex  The data index that the iteration primitive was at previously.  This value is updated
     * to contain the new data index.
     */
    void increment_to_x( boost::int32_t xPrevious, boost::int32_t xNew, const rle_index_spec& ris,
                         boost::int32_t& dataIndex );

    /**
     * Skips ahead (forward only) to the specified X in the same scanline.
     *
     * @param  xPrevious  The x coordinate of the voxel the iteration primitive is at before this increment.
     * @param  xNew  The x coordinate of the voxel the iteration primitive is at after this increment.
     * @param  ris  The rle_index_spec which was used to initialize the primitive.
     * @param[out]  outXNegDataIndex  This is populated with the data index of the voxel with one smaller X coordinate
     * than the one resulting from the increment.
     * @param[in,out]  dataIndex  The data index that the iteration primitive was at previously.  This value is updated
     * to contain the new data index.
     * @param[out]  outXPosDataIndex  This is populated with the data index of the voxel with one greater X coordinate
     * than the one resulting from the increment.
     */
    void increment_to_x( boost::int32_t xPrevious, boost::int32_t xNew, const rle_index_spec& ris,
                         boost::int32_t& outXNegDataIndex, boost::int32_t& dataIndex,
                         boost::int32_t& outXPosDataIndex );

    //////
    // Properties of the iteration
    //////

    // The X coordinate and data index are stored by the client of this object, so there are no getters for them.

    /**
     * Whether the iteration primitive has gone past all of the runs it had to process.
     */
    bool finished() const { return m_rd == 0; }
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
