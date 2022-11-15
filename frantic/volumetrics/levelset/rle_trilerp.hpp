// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/size3.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_named_channels.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * Trilerp and gradient function for floating point rle map channel (or signed distance array). It takes a float*
 * for its data so that either a rle_channel_accessor or a signed distance array can be used. The two output values
 * are passed in as pointers so that the user can specify which outputs to compute (no output will be computed for
 * NULL inputs). This function does a 2x2x2 index box fetch from the provided rle_index_spec. If the user is doing
 * multiple trilerps at a single voxel coodinate, the voxel indices should be pre-computed and used with
 trilerp_float_from_indices.
 * <p>
 * Example Usage for doing a trilerp and gradient of the signed distance at a given voxel position:
 *     //given a populated rle_level_set rls, and desired vector3f voxelCoordinate
 *     float trilerpVal;
 *     vector3f gradientVal;
 *     bool success = trilerp_float( rls.get_rle_index_spec(), &rls[0], voxelCoordinate, &trilerpVal, &gradientVal );
 *
 * @param ris  The rle index spec that defines the layout of the float data array.
 * @param dataAccessor  The actual float array of data. Usually this will be the refernce to the start of the signed
 *                      distance array, or a reference to the start of an rle_channel_accessor<float> buffer.
 * @param voxelCoordinate  The coordinate in voxel space at which to do the trilerp.
 * @param outTrilerp  If the trilerp was successful, this value will be set to the resulting trilerp'ed value.
                      If the trilerp was unsuccessful and the voxels have undefined code -1, this value is set to a
 positive number (for signed distance convenience). If the trilerp was unsuccessful and the voxels have undefined code <
 -1, this value is set to a negative number (for signed distance convenience). If NULL is passed in, no trilerp will be
 computed.
 * @param outGradient  The resulting gradient value. If NULL is passed in, no gradient will be computed.
 * @return  false if there were no defined voxels in the voxelCoordinate's 2x2x2 voxel neighbourhood, and as a result,
 no meaningful
 *          output could be computed. true if one or more voxels were defined, and the function is returning meaningful
 output.
 */
bool trilerp_float( const frantic::volumetrics::levelset::rle_index_spec& ris, const float* dataAccessor,
                    const frantic::graphics::vector3f& voxelCoordinate, float* outTrilerp = NULL,
                    frantic::graphics::vector3f* outGradient = NULL );

/**
 * Trilerp and gradient function for floating point rle map channel (or signed distance array). The indices array must
 * be populated with indicesVoxelBounds.get_volume() number of values and must encapsulate the 2x2x2 data index box
 needed
 * to do the trilerp at the requested coordinate. Negative indicies will be handled correctly. The two output values are
 * passed in as pointers so that the user can specify which outputs to compute (no output will be computed for NULL
 inputs).
 * <p>
 * Example Usage for doing a trilerp of the signed distance at a given voxel position between [2.5,2.5,2.5] and
 [3.5,3.5,3.5]:
 *     //given a populated rle_level_set rls, and desired vector3f voxelCoordinate
 *     boundbox3 indicesBounds( vector3(2,2,2), vector3(3,3,3) );
 *     std::vector<boost::int32_t> indices( indicesBounds.get_volume() );
 *     rls.get_rle_index_spec().fill_data_index_box( indicesBounds, &indices[0] );
 *     float trilerp;
 *     trilerp_float_from_indices( &rls[0], voxelCoordinate, &indices[0], indicesBounds, &trilerp );
 *
 * @param dataAccessor  The float array of data. Usually this will be the refernce to the start of the signed
 *                      distance array, or a reference to the start of an rle_channel_accessor<float> buffer.
 * @param voxelCoordinate  The coordinate in voxel space at which to do the trilerp.
 *                         This coodinate must be within the interpolation bounds of the indicesVoxelBounds.
 * @param indices  An array of size indicesVoxelBounds.get_volume() containing the relevant indices into the
 dataAccessor array.
 * @param indicesVoxelBounds  The voxel bounding box that defines the size of the indices array, and were the indices
 array values correspond to in voxel space.
 * @param outTrilerp  If the trilerp was successful, this value will be set to the resulting trilerp'ed value.
                      If the trilerp was unsuccessful and the voxels have undefined code -1, this value is set to a
 positive number (for signed distance convenience). If the trilerp was unsuccessful and the voxels have undefined code <
 -1, this value is set to a negative number (for signed distance convenience). If NULL is passed in, no trilerp will be
 computed.
 * @param outGradient  The resulting gradient value. If NULL is passed in, no gradient will be computed.
 * @return  false if there were no defined voxels in the voxelCoordinate's 2x2x2 voxel neighbourhood, and as a result,
 no meaningful
 *          output could be computed. true if one or more voxels were defined, and the function is returning meaningful
 output.
 */
bool trilerp_float_from_indices( const float* dataAccessor, const frantic::graphics::vector3f& voxelCoordinate,
                                 const boost::int32_t* indices, const frantic::graphics::boundbox3& indicesVoxelBounds,
                                 float* outTrilerp = NULL, frantic::graphics::vector3f* outGradient = NULL );

/**
 * Trilerp function for a vector3f rle map channel. This function does a 2x2x2 index box fetch
 * from the provided rle_index_spec. If the user is doing multiple trilerps at a single voxel
 * coodinate, the voxel indices should be pre-computed and used with trilerp_vector3f_from_indices.
 * <p>
 * Example Usage for doing a trilerp and gradient of the "Color" map channel at a given voxel position:
 *     //given a populated rle_level_set rls with "Color" channel, and desired vector3f voxelCoordinate
 *     const_rle_channel_accessor<vector3f> colorAcc = rls.get_rle_channel_accessor<vector3f>( "Color" );
 *     vector3f trilerpColorVal;
 *     bool success = trilerp_vector3f( rls.get_rle_index_spec(), colorAcc, voxelCoordinate, trilerpColorVal );
 *
 * @param ris  The rle index spec that defines the layout of the vector3f data.
 * @param dataAccessor  The vector3f rle data accessor.
 * @param voxelCoordinate  The coordinate in voxel space at which to do the trilerp.
 * @param outTrilerp  The resulting trilerp'ed vector3f value.
 * @return  false if there were no defined voxels in the voxelCoordinate's 2x2x2 voxel neighbourhood, and as a result,
 * no meaningful output could be computed. true if one or more voxels were defined, and the function is returning
 * meaningful output.
 */
bool trilerp_vector3f(
    const frantic::volumetrics::levelset::rle_index_spec& ris,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f> dataAccessor,
    const frantic::graphics::vector3f& voxelCoordinate, frantic::graphics::vector3f& outTrilerp );

/**
 * Trilerp function for a vector3f rle map channel. The indices array must be populated with
 * indicesVoxelBounds.get_volume() number of values and must encapsulate the 2x2x2 data index box
 * needed to do the trilerp at the requested coordinate. Negative indicies will be handled correctly.
 * <p>
 * Example Usage for doing a trilerp of the "Color" channel at a given voxel position between [2.5,2.5,2.5] and
 * [3.5,3.5,3.5]:
 *     //given a populated rle_level_set rls, and desired vector3f voxelCoordinate
 *     boundbox3 indicesBounds( vector3(2,2,2), vector3(3,3,3) );
 *     std::vector<boost::int32_t> indices( indicesBounds.get_volume() );
 *     rls.get_rle_index_spec().fill_data_index_box( indicesBounds, &indices[0] );
 *     const_rle_channel_accessor<vector3f> colorAcc = rls.get_rle_channel_accessor<vector3f>( "Color" );
 *     vector3f trilerpColorVal;
 *     trilerp_vector3f_from_indices( colorAcc, voxelCoordinate, &indices[0], indicesBounds, trilerpColorVal );
 *
 * @param dataAccessor  The vector3f rle data accessor.
 * @param voxelCoordinate  The coordinate in voxel space at which to do the trilerp.
 *                         This coodinate must be within the interpolation bounds of the indicesVoxelBounds.
 * @param indices  An array of size indicesVoxelBounds.get_volume() containing the relevant indices into the
 * dataAccessor data.
 * @param indicesVoxelBounds  The voxel bounding box that defines the size of the indices array, and were the indices
 * array values correspond to in voxel space.
 * @param outTrilerp  The resulting trilerp'ed vector3f value.
 * @return  false if there were no defined voxels in the voxelCoordinate's 2x2x2 voxel neighbourhood, and as a result,
 * no meaningful output could be computed. true if one or more voxels were defined, and the function is returning
 * meaningful output.
 */
bool trilerp_vector3f_from_indices(
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f> dataAccessor,
    const frantic::graphics::vector3f& voxelCoordinate, const boost::int32_t* indices,
    const frantic::graphics::boundbox3& indicesVoxelBounds, frantic::graphics::vector3f& outTrilerp );

/**
 * Trilerp function for a vector3f rle map channel that contains "Staggered" data (data stored at voxel faces, not at
 * voxel centers). This function does a index box fetch from the provided rle_index_spec of between 2x2x2 and 3x3x3
 * voxels depending on what is needed. If the user is doing multiple trilerps for staggered data at a single voxel
 * coodinate, the voxel indices should be pre-computed and used with trilerp_vector3f_staggered_from_indices.
 * <p>
 * Example: See documentation for trilerp_vector3f.
 *
 * Parameters: See documentation for trilerp_vector3f.
 */
bool trilerp_vector3f_staggered(
    const frantic::volumetrics::levelset::rle_index_spec& ris,
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f> dataAccessor,
    const frantic::graphics::vector3f& voxelCoordinate, frantic::graphics::vector3f& outTrilerp );

/**
 * Trilerp function for a vector3f rle map channel that contains "Staggered" data (data stored at voxel faces, not at
 * voxel centers). The indices array must be populated with indicesVoxelBounds.get_volume() number of values and must
 * encapsulate the data index box required for the staggered trilerp lookup. The size of the data index box required
 * varies depending on where in a voxel the lookup is being done at, however, a 3x3x3 box centered at the
 * voxelCoordinate's voxel will work for all cases. Negative indicies will be handled correctly.
 * <p>
 * Example: See documentation for trilerp_vector3f_staggered_from_indices.
 *
 * Parameters: See documentation for trilerp_vector3f_from_indices.
 */
bool trilerp_vector3f_staggered_from_indices(
    frantic::volumetrics::levelset::const_rle_channel_accessor<frantic::graphics::vector3f> dataAccessor,
    const frantic::graphics::vector3f& voxelCoordinate, const boost::int32_t* indices,
    const frantic::graphics::boundbox3& indicesVoxelBounds, frantic::graphics::vector3f& outTrilerp );

/**
 * Trilerp function for a general channel accessor. This function does a 2x2x2 index box fetch from
 * the provided rle_index_spec. If the user is doing multiple trilerps at a single voxel coordinate,
 * the voxel indices should be pre-computed and used with trilerp_general_channel_from_indices.
 * <p>
 * Example Usage: See trilerp_vector3f comment.
 *
 * @param ris  The rle index spec that defines the layout of the vector3f data.
 * @param dataAccessor  The general channel accessor.
 * @param voxelCoordinate  The coordinate in voxel space at which to do the trilerp.
 * @param outTrilerp  The resulting trilerp'ed value.
 * @return  false if there were no defined voxels in the voxelCoordinate's 2x2x2 voxel neighbourhood, and as a result,
 * no meaningful output could be computed. true if one or more voxels were defined, and the function is returning
 * meaningful output.
 */
bool trilerp_general_channel( const frantic::volumetrics::levelset::rle_index_spec& ris,
                              frantic::volumetrics::levelset::const_rle_channel_general_accessor dataAccessor,
                              const frantic::graphics::vector3f& voxelCoordinate, char* outTrilerp );

/**
 * Trilerp function for a general channel accessor. The indices array must be populated with
 * indicesVoxelBounds.get_volume() number of values and must encapsulate the 2x2x2 data index box
 * needed to do the trilerp at the requested coordinate. Negative indicies will be handled correctly.
 * <p>
 * Example Usage: See trilerp_vector3f_from_indices comment.
 *
 * @param dataAccessor  The general channel accessor.
 * @param voxelCoordinate  The coordinate in voxel space at which to do the trilerp.
 *                         This coodinate must be within the interpolation bounds of the indicesVoxelBounds.
 * @param indices  An array of size indicesVoxelBounds.get_volume() containing the relevant indices into the
 * dataAccessor data.
 * @param indicesVoxelBounds  The voxel bounding box that defines the size of the indices array, and were the indices
 * array values correspond to in voxel space.
 * @param outTrilerp  The resulting trilerp'ed value.
 * @return  false if there were no defined voxels in the voxelCoordinate's 2x2x2 voxel neighbourhood, and as a result,
 * no meaningful output could be computed. true if one or more voxels were defined, and the function is returning
 * meaningful output.
 */
bool trilerp_general_channel_from_indices(
    frantic::volumetrics::levelset::const_rle_channel_general_accessor dataAccessor,
    const frantic::graphics::vector3f& voxelCoordinate, const boost::int32_t* indices,
    const frantic::graphics::boundbox3& indicesVoxelBounds, char* outTrilerp );

/**
 * Retrieves the weights and data indices needed to trilinearly interpolate
 * a channel for an input coordinate.  This can be used together with a
 * channel_weighted_sum_combine_function_t function to do a trilerp
 * of an arbitrary input channel.
 * <p>
 * The output outDataIndices arrays must contain at least 8
 * elements.  The outDataIndices array may contain both positive and negative
 * values.  If a value in this array is negative, then it is a region code and
 * an appropriate default value should be blended.
 *
 * @param  ris  The index spec that provides the data indices.
 * @param  voxelCoordinate  The voxel coordinate at which to do the trilerp lookup.
 * @param  outTrilerpDeltas  An output parameter, where 3 deltas will go.
 * @param  outDataIndices  An output parameter, where 8 data indices will go.
 */
void get_trilerp_indices( const frantic::volumetrics::levelset::rle_index_spec& ris,
                          const frantic::graphics::vector3f& voxelCoordinate, float* outTrilerpDeltas,
                          boost::int32_t* outDataIndices );

/**
 * This function should be used in conjunction with get_trilerp_indices.
 * <p>
 * Function will take in the three weights and return the multipliers associated with each of the 8
 * trilerp indices. Using the output weightsArray and dataIndicies (from get_weights) a trilinear
 * interpolation can be done using a weighted_sum_function.
 * <p>
 * 1) get_trilerp_indices ( voxelCoordinate, outTrilerpDeltas, outDataIndices );
 * 2) get_trilerp_weights ( outTrilerpDeltas, weightsArray );
 * 3) Create array of pointers to data;
 * 4) weighted_sum_function( weightsArray, dataPointers, 8, arity, result );
 *
 * @param  deltas  deltas given by get_trilerp_indices
 * @param  outWeights  array size 8 to hold the multipliers
 */
void get_trilerp_weights( const float* deltas, float* outWeights );

/**
 * Modifies a 2x2x2 weights array used for trilerp'ing to ensure that all undefined voxels
 * will receive a zero weight, and the remaining defined voxels weights will still sum to one.
 *
 * @param dataIndices  An array of data indices into a rle structure that are being used for a trilerp.
 * @param outWeights  An array of size 2x2x2 that represents the weights for triler'ing 2x2x2 data values.
 */
void normalize_trilerp_weights_for_undefined_voxels( const boost::int32_t* dataIndices, float* outWeights );

/**
 * This function returns the gradient of the distance field, by creating a staggered grid of the
 * signed distance derivatives around the sample point. The sample point's gradient is then interpolated
 * using a trilerp on a staggered grid for each component (as opposed to a lerp per component). It is
 * second order accurate because we are taking centered differences to create the staggered grid.
 * @note It is more expensive, but smoother that the standard forward-only differencing such
 *        as the one used in trilerp_float.
 * @note This function has been written to work only on an rle_level_set, and undefined
 *        inside/outside values are considered to be the level set's inside/outside interface distance
 *        (unlike the trilerp_float function which handles undefined voxels correctly).
 *
 * @param levelSet The rle_level_set that provides the signed distance data.
 * @param voxelCoordinate The voxel coordinate at which to do the trilerp lookup.
 * @param outSignedDistanceGradient An output value that contains the computed gradient. This gradient
 *                                   is with respect to the change in voxel coordinate (NOT WORLD COORD).
 */
void trilerp_staggered_centered_signed_distance_gradient( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                                          const frantic::graphics::vector3f& voxelCoordinate,
                                                          frantic::graphics::vector3f& outSignedDistanceGradient );

/**
 *  Accelerate consecutive trilerp lookups within the same voxel.
 */
class cached_trilerp {
    frantic::graphics::vector3f m_cachedVoxelCoordinate;
    frantic::graphics::vector3 m_cachedMinCorner;

    boost::int32_t m_definedCount;

    boost::int32_t m_definedDataIndices[8];
    float m_definedWeights[8];

    boost::int32_t m_dataIndices[8];
    float m_weights[8];

    const rle_index_spec& m_ris;

    void fill_cache( const frantic::graphics::vector3f& voxelCoordinate );

    bool get( const char* data, std::size_t arity, data_type_t dataType,
              const frantic::graphics::vector3f& voxelCoordinate, char* outSample );

    cached_trilerp& operator=( cached_trilerp& ); // not implemented

  public:
    /**
     *  Construct an object for performing trilerp lookups in the specified
     * ris.
     *
     * @param ris  The rle_index_spec that defines the layout of the defined
     *			   voxels.
     */
    cached_trilerp( const rle_index_spec& ris );

    /**
     *  Perform a trilinear interpolation of channel data.
     *
     * @todo specialize this for common data types.
     *
     * @param data  An array of channel data that will be interpolated.
     *				This array must have the same size as the rle_index_spec
     *				that was passed to the constructor.
     * @param voxelCoordinate  The voxel coordinate at which to do the trilerp
     *						   lookup.
     * @param[out] outSample  If the function returns true, then this is the
     *					 resulting trilerped value.  If the function returned
     *					 false, then this is undefined.
     * @return  false if there were no defined voxels in the voxelCoordinate's
     *			2x2x2 voxel neighbourhood, and as a result, no meaningful
     *			output could be computed. true if one or more voxels were
     *			defined, and the function is returning meaningful output.
     */
    template <class T>
    bool get( const T* data, const frantic::graphics::vector3f& voxelCoordinate, T& outSample ) {
        return get( reinterpret_cast<const char*>( data ), frantic::channels::channel_data_type_traits<T>::arity(),
                    frantic::channels::channel_data_type_traits<T>::data_type(), voxelCoordinate,
                    reinterpret_cast<char*>( &outSample ) );
    }

    /**
     *  Perform a trilinear interpolation of data in a channel accessor.
     *
     * @param dataAccessor  The channel accessor that holds the data to interpolate.
     *				This accessor must have the same size as the rle_index_spec
     *				that was passed to the constructor.
     * @param voxelCoordinate  The voxel coordinate at which to do the trilerp
     *						   lookup.
     * @param[out] outSample  If the function returns true, then this will hold
     *					 the resulting interpolated value.  If the function
     *					 returns false, then this is undefined.
     * @return  false if there were no defined voxels in the voxelCoordinate's
     *			2x2x2 voxel neighbourhood, and as a result, no meaningful
     *			output could be computed. true if one or more voxels were
     *			defined, and the function is returning meaningful output.
     */
    template <class T>
    bool get( const rle_channel_accessor<T>& dataAccessor, const frantic::graphics::vector3f& voxelCoordinate,
              T& outSample ) {
        return get( reinterpret_cast<const char*>( dataAccessor.data( 0 ) ),
                    frantic::channels::channel_data_type_traits<T>::arity(),
                    frantic::channels::channel_data_type_traits<T>::data_type(), voxelCoordinate,
                    reinterpret_cast<char*>( &outSample ) );
    }

    /**
     *  Perform a trilinear interpolation of data in a channel accessor.
     *
     * @param dataAccessor  The channel accessor that holds the data to interpolate.
     *				This accessor must have the same size as the rle_index_spec
     *				that was passed to the constructor.
     * @param voxelCoordinate  The voxel coordinate at which to do the trilerp
     *						   lookup.
     * @param[out] outSample  If the function returns true, then this will hold
     *					 the resulting interpolated value.  If the function
     *					 returns false, then this is undefined.
     * @return  false if there were no defined voxels in the voxelCoordinate's
     *			2x2x2 voxel neighbourhood, and as a result, no meaningful
     *			output could be computed. true if one or more voxels were
     *			defined, and the function is returning meaningful output.
     */
    template <class T>
    bool get( const const_rle_channel_accessor<T>& dataAccessor, const frantic::graphics::vector3f& voxelCoordinate,
              T& outSample ) {
        return get( reinterpret_cast<const char*>( dataAccessor.data( 0 ) ),
                    frantic::channels::channel_data_type_traits<T>::arity(),
                    frantic::channels::channel_data_type_traits<T>::data_type(), voxelCoordinate,
                    reinterpret_cast<char*>( &outSample ) );
    }

    /**
     *  Perform a trilinear interpolation of data in a channel accessor.
     *
     * @param dataAccessor  The channel accessor that holds the data to interpolate.
     *				This accessor must have the same size as the rle_index_spec
     *				that was passed to the constructor.
     * @param voxelCoordinate  The voxel coordinate at which to do the trilerp
     *						   lookup.
     * @param[out] outSample  If the function returns true, then this will hold
     *					 the resulting interpolated value.  If the function
     *					 returns false, then this is undefined.  This must not
     *					 be a NULL pointer, and the memory required to hold the
     *					 output must be allocated before calling.
     * @return  false if there were no defined voxels in the voxelCoordinate's
     *			2x2x2 voxel neighbourhood, and as a result, no meaningful
     *			output could be computed. true if one or more voxels were
     *			defined, and the function is returning meaningful output.
     */
    bool get( const const_rle_channel_general_accessor& dataAccessor,
              const frantic::graphics::vector3f& voxelCoordinate, char* outSample );

    /**
     *  Perform a trilinear interpolation of signed distance data in a level
     * set.
     *
     *  If at least one voxel is defined in the voxelCoordinate's 2x2x2
     * neighbourhood, then the interpolation uses input only from the
     * defined neighbours.  If there is no defined voxel in the
     * voxelCoordinate's neighbourhood, then the interpolation uses the
     * level set's default inside distance or default outside distance
     * instead.
     *
     * @param ls The level set that holds the signed distance data to
     *			 interpolate.  It should have the same rle_index_spec
     *			 as the one that was passed to the constructor.
     * @param voxelCoordinate  The voxel coordinate at which to do the
     *						   interpolation.
     * @return  The interpolated signed distance at the specified
     *			voxelCoordinate.
     */
    float get( const rle_level_set& ls, const frantic::graphics::vector3f& voxelCoordinate );
};

/**
 *  This class accelerates consecutive trilerp_vector3f_staggered lookups
 * within the same voxel.
 */
class cached_trilerp_vector3f_staggered {
    const rle_index_spec& m_ris;
    frantic::graphics::boundbox3 m_cachedBounds;
    boost::int32_t m_dataIndices[27];

    frantic::graphics::boundbox3f get_trilerp_vector3f_staggered_sample_range_from_indices(
        const const_rle_channel_accessor<frantic::graphics::vector3f>& dataAccessor,
        const frantic::graphics::vector3f& voxelCoordinate, const boost::int32_t* indices,
        const frantic::graphics::boundbox3& indicesVoxelBounds );

    cached_trilerp_vector3f_staggered& operator=( cached_trilerp_vector3f_staggered& ); // not implemented

  public:
    /**
     *  Construct an object for performing a trilerp lookup in the specified
     * rle_index_spec.
     *
     * @param ris  The rle index spec that defines the layout of the vector3f data.
     */
    cached_trilerp_vector3f_staggered( const rle_index_spec& ris )
        : m_ris( ris ) {}

    /**
     *  Sample the value of the staggered vector3f rle map channel at the
     * specified voxelCoordinate.
     *
     * @note the implementation is taken from trilerp_vector3f_staggered
     *
     * @param dataAccessor  The vector3f rle data accessor, which holds
     *      staggered data.
     * @param voxelCoordinate the voxel coordinate at which to do the trilerp
     *		lookup.
     * @param[out] outSample if the function returned true, then this is the
     *		resulting trilerped value.  If the function returned false, then
     *		this is 0.
     * @return  false if there were no defined voxels in the voxelCoordinate's 2x2x2 voxel neighbourhood, and as a
     *result, no meaningful output could be computed. true if one or more voxels were defined, and the function is
     *returning meaningful output.
     */
    bool get( const const_rle_channel_accessor<frantic::graphics::vector3f>& dataAccessor,
              const frantic::graphics::vector3f& voxelCoordinate, frantic::graphics::vector3f& outSample );

    /**
     *  Return the range of sample values used to calculate each component of
     * a trilerp result.
     *
     * @param dataAccessor the channel accessor in which to find the sample
     *		range.
     * @param voxelCoordinate the voxel coordinate at which to find the sample
     *		range.
     * @return the range of sample values used to calculate the trilerp at the
     *		specified voxelCoordinate.  If no sample is available for a
     *		component, then it is indicated by
     *		boundbox3f.maximum()[component] < boundbox3f.minimum()[component]
     */
    frantic::graphics::boundbox3f
    get_sample_range( const const_rle_channel_accessor<frantic::graphics::vector3f>& dataAccessor,
                      const frantic::graphics::vector3f& voxelCoordinate );
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
